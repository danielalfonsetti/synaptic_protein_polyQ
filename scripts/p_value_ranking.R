####################################################
# Daniel Alfonsetti, daniel.alfonsetti@gmail.com
# MIT, Littleton Lab
# 24 January 2019
# ---------------------------------
# Description: A script that takes the output of the HMM and quantifies the significance 
# of the polyQ tract found by the HMM for a given protein, and then continues
# to perform wilcoxon ranked sum tests and gene set enrichment analysis to see if
# the category of synaptic proteins is enrniched for polyQ proteins. 
####################################################
rm(list = ls())

# Load in libraries
library(dplyr)

# Biological databases and their tools
library(biomaRt)
library(biomartr)
library(STRINGdb)
library(org.Dm.eg.db)
library(org.Hs.eg.db) 
library(org.Mm.eg.db)
library(org.Ce.eg.db)
library(GO.db)

# GSEA library
library(clusterProfiler)
####################################################
# Load in the data as outputted from the HMM.
data <- read.csv("C:/UROPs/polyQ_neuronal_proteins/output/fly/Q/fly_22905prots_w_HMM_Q.csv")
data <- read.csv("C:/UROPs/polyQ_neuronal_proteins/output/human/Q/human_78328prots_w_HMM_Q.csv")

# Add a peptide length column.
data$PeptideLength <- unlist(lapply(as.character(data$peptide_seq), function(x){length(unlist(strsplit(x, "")))}))
####################################################

# GOAL: Find a way to measure sigificance of polyQ proteins outputted by our HMM while controlling for length.
# Ideally, would be able to order the proteins by significance, allowing for more resolution
# besides just a boolean of "is a polyQ"/"is not a polyQ".

# Proposed Method:
# --------------------------
# Let m be the length of a protein
# Let n be the length of the largest polyQ tract found by the HMM

# There are 20**n different possible amino acid sequences of length n.
# This implies there is a 1/(20**n) chance of seeing a specific n-mer at a specific location.
# Thus there is a 1-1/(20**n) chance of not seeing the specific n-mer at a specific location.
# There are m-(n-1) possible locations for a given n-mer to be located on a protein of length m.
# Therefore, there is a (1-1/20**n)**(m-(n-1)) chance of not seeing a specific n-mer anywhere on a protein of length m.
# Therefore, there is a 1-(1-1/20**n)**(m-(n-1)) chance of seeing the n-mer at least once on the protein.
# This is our p-value for a given protein of length m annotated by the HMM to have a maximum polyQ tract of length n.

# This p-value is testing the following for a given protein of length m:

  # Null hypothesis:  The m probability distributions for each of the m locations along the protein that define that 
  # locations' respective probabilities of being a particular amino acid are uniform and independent.
  # (i.e. each location has a 1/20 probability of being any amino acid)
  
  # Alternative hypothesis (? still debating the exact wording here?): 
  # The m probability distributions are not uniform and independent.
  # E.g. knowing that one location is a "Q" may affect the probability that another location is a "Q".

# We are able to calculate such a p-value for each protein in our dataset, 
# thus giving us a metric of determining relative significance of polyQ proteins (instead of a boolean of polyQ vs not polyQ)
# and thereby giving us a ranked order. With this ranked order, we can use certain non-parametric tests 
# such as the Wilcoxon ranked sum test as well as gene set enrichment analysis (GSEA) - as shown below.


n = ifelse(is.na(data$MaxLengthsPolyAA), 0, data$MaxLengthsPolyAA) # get max length of polyQ tract found for each protein
m = data$PeptideLength # Protein length
data$pvalues <-  1-(1-1/(20**n))**(m-(n-1)) # compute p-values

data <- data %>% arrange(pvalues) # order dataset by p-values

##################################
# WILCOXON (single sided) RANKED SUM TEST
# https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test
##################################

  # Null hypothesis: Roughly speaking, the length of 
  # the longest polyQ tract found by our HMM has no bearing on whether
  # the protein is synaptic or not. More formally, the wilcoxon test null hyp. is that the probability of a randomly selected
  # synaptic protein having a higher p-value ranking than a randomly selected non-synaptic protein is equal to the probability 
  # it is lower.

  # Alternative hypothesis: Roughly speaking, the longer the length of the maximum polyQ tract found by our HMM,
  # the more likely that protein is synaptic. More formally, the probability of a randomly selected synaptic protein having
  # a higher p-value ranking than a randomly selected non-syaptic protein is greater than the probability it is lower.

# Annotate proteins based on whether they are synaptic or not.
synapse_term1 <- "GO:0045202"
synapse_cats <- get(synapse_term1, GOCCOFFSPRING)
synapse_cats <- c(synapse_cats, synapse_term1)

data$synaptic <- unlist(
  lapply(
    lapply(as.character(data$go_ids_CC), function(x){unlist(strsplit(x, "; "))}), 
    function(x){any(x %in% synapse_cats)}
  )
)

# Do wilcoxon ranked sum test on p-values between synaptic proteins and non-synaptic proteins.
x <- data %>% filter(synaptic == TRUE)
y <- data %>% filter(synaptic == FALSE)

# Run test.
test_result <- wilcox.test(x$pvalues, y$pvalues, alternative = "less")
test_result$p.value # 7.508342e-12 (fly). Very significant.
# We reject the null hypothesis.

# Sanity check/control for wilcoxon test. Should expect a uniform distribution of p-values if
# proteins are randomly grouped into two categories (i.e. if the null hypothesis holds).
# https://stats.stackexchange.com/questions/10613/why-are-p-values-uniformly-distributed-under-the-null-hypothesis
# This is what we see. Good.

wilcox_p_values <- vector(mode='numeric', length=1000)
for (i in 1:1000){
  rows <- sample.int(nrow(data), size = nrow(x), replace = FALSE)
  x_control <- data[rows,]
  y_control <- data[-rows,]
  control_result <- wilcox.test(x_control$pvalues, y_control$pvalues, alternative = "less")
  wilcox_p_values[[i]] <- control_result$p.value
}
hist(wilcox_p_values)



##################################
# Do Kruskal Wallis test
# Test to see if there is a statistically significant difference
# between nuclear, PSD, and AZ proteins. 
##################################


############################



##################################
# Gene Set Enrichment Analysis (GSEA) using clusterProfiler library
# https://en.wikipedia.org/wiki/Gene_set_enrichment_analysis
# https://bioconductor.org/packages/devel/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#go-gene-set-enrichment-analysis


# GSEA has an advantage over standard enrichment analysis in that it can catch if many small changes (although not necessairly individually
# significant) are occuring for proteins that are in a functional group. It is normally used for differential expression,
# but I believe it can be used here since our p-values (calculated above) give us a ranking just like fold-change in DE would.
# The idea is that the 'more polyQ' proteins (as determined by a higher p-value rank) will be influencing which categories are preferential (or not) 
# to polyQ proteins better than "less polyQ" proteins, thus giving us better resolution than standard enrichment analysis would.
##################################

# Only stick to using the fly proteome for now.
species = "fly" 
if (species == "fly") {
  OrgDb = org.Dm.eg.db
  fromType = "FLYBASEPROT"
} else if (species == "human") {
  OrgDb = org.Hs.eg.db
  fromType = "ENSEMBLPROT"
} else if (species == "mouse") {
  OrgDb = org.Mm.eg.db
  fromType = "ENSEMBLPROT"
} else if (species == "worm") {
  OrgDb = org.Ce.eg.db
  fromType = "ENSEMBLPROT"
}

# Do enrichment for each subontology of GO
for (go_ont in c("BP", "MF", "CC")) {
  proteins <- data
  print(paste0("Doing enrichment analysis for ", go_ont, " ontology."))
  
  # Convert from whatever protein identifier was used in OrgDb to ENTREZID for the gseGO function.
  id_translation_df <- bitr(proteins$ensembl_peptide_id,
                      fromType = fromType,
                      toType = c("ENTREZID", "SYMBOL"),
                      OrgDb = OrgDb)
  
  # Only use proteins in our ranked list that we have a mapping for
  proteins <- proteins[proteins$ensembl_peptide_id %in% eval(parse(text=paste0("id_translation_df$", fromType))),]
  proteins <- merge(proteins, id_translation_df, by.x = "ensembl_peptide_id", by.y = fromType)
  
  # Order by p-values.
  proteins <- proteins %>% arrange(pvalues) 
  
  # Format ranked order vector as gseGO function requires.
  geneList <- 1-proteins$pvalues
  names(geneList) <- proteins$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)
  
  # Perform test
  go_gse_output_test <- gseGO(gene = geneList,
                                    OrgDb = OrgDb,
                                    ont = go_ont,
                                    nPerm = 1000,
                                    minGSSize = 5,
                                    maxGSSize = 500,
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = .3,
                                    verbose = TRUE)
  # Save results
  write.csv(go_gse_output_test, file = paste0("C:/UROPs/polyQ_neuronal_proteins/output/gsea_test_output/gsea_go_", go_ont, ".csv"))
}
# No groups enriched were enriched using GSEA. Strange, considering 
# that wilcoxon ranked sum test and standard enrichment analysis
# validate the hypothesis that synaptic proteins are enriched for polyQ tracts.
# This discrepancy may be due to the ranking metric I am using as input to the GSEA function.

##################################
##################################
# Junk/Misc code (skip)
  # kegg_gse_output_test <- gseKEGG(gene = geneList,
  #                                organism = "Fly",
  #                                ont = go_ont,
  #                                nPerm = 100,
  #                                minGSSize = 100,
  #                                maxGSSize = 500,
  #                                pAdjustMethod = "BH",
  #                                pvalueCutoff = 0.05,
  #                                verbose = TRUE)
  # 
  # mkegg_gse_output_test <- gseMKEGG(gene = geneList,
  #                               OrgDb = OrgDb,
  #                               ont = go_ont,
  #                               nPerm = 100,
  #                               minGSSize = 100,
  #                               maxGSSize = 500,
  #                               pAdjustMethod = "BH",
  #                               pvalueCutoff = 0.05,
  #                               verbose = TRUE)
  
  # go_enrich_output_test <- enrichGO(gene = proteins$ENTREZID,
  #                                   OrgDb = OrgDb,
  #                                   ont = go_ont,
  #                                   pAdjustMethod = "BH",
  #                                   pvalueCutoff = 0.05,
  #                                   qvalueCutoff = 0.05,
  #                                   readable = TRUE)
  # write.csv(go_enrich_output_test, file = paste0("C:/UROPs/polyQ_neuronal_proteins/output/gsea_test_output/overrep_go_", go_ont, ".csv"))
  # 
  #   
