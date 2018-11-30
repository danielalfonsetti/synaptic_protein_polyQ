# Daniel Alfonsetti, daniel.alfonsetti@gmail.com
# 13 Sep. 2018
# Description: polyQ analysis in synaptic proteins

##########################
library(HiddenMarkov)
library(stringr)
# library(dcGOR)
# library(dnet)
# library(supraHex)
# library(Biobase)
# library(Rgraphviz)
# library(domainsignatures)
library(grid)
library(dplyr)
library(ggplot2)
library(doParallel)

require(data.table)

library(biomaRt)
library(clusterProfiler)
library(STRINGdb)
library(HiddenMarkov)

library(org.Dm.eg.db) # Drosophila_melanogaster

##########################
# Setup environment
rm(list = ls())

setwd("C:/UROPs/polyQ_neuronal_proteins")
plot_output_dir <- "C:/UROPs/polyQ_neuronal_proteins/plots/"
data_dir <- "C:/UROPs/polyQ_neuronal_proteins/data/"
output_dir <-      "C:/UROPs/polyQ_neuronal_proteins/output/"
registerDoParallel(cores = 4)
##############################################################################


#HIDDEN MARKOV MODEL PARAMETERS
HiddenMarkov = function (proteinSet){
  
  # Debugging
    # proteinSet = b[1:10,]
  #
  
  x = NULL #string translation sequence
  Pi = matrix(c(0.95, 0.05, 0.13, 0.87), nrow = 2, ncol = 2, byrow = TRUE) #hardcode transition prob matrix
  delta = c(.99, .01) #probability of starting in either state (nonQ region, Q region)
  distn = "binom" # distribution name, average = probability
  pm = list(prob = c(0.05, 0.26), size = 1) #prob is parameter for binom dist, prob of Q in nonQ, Q states
  HMM = dthmm(x, Pi, delta, distn, pm, pn = NULL, discrete = TRUE, nonstat = TRUE) # Create the model
  
  
  # A column for the Qfraction for an entire protein.
  proteinSet$Qfraction = rep(NA, nrow(proteinSet))
  
  # A boolean indicating whether or not at least one polyQ region was found by the HMM for this protein.
  proteinSet$HMMhasPolyQ = rep(NA, nrow(proteinSet))
  
  # A column where each element is a number indicating the number of polyQ regions found on this protein by the HMM.
  proteinSet$NumberPolyQ = rep(NA, nrow(proteinSet))
  
  # A column where each element contains a list with the indicies for each region in the protein that is denoted polyQ by the HMM
  proteinSet$IndiciesPolyQ = rep(NA, nrow(proteinSet)) 
  
  # A column where each element contains a list with the AA sequence for each region in the protein that is denoted polyQ by the HMM
  proteinSet$TranslationsPolyQ = rep(NA, nrow(proteinSet)) 
  
  # A column where each contains a list of the length of the region that is denoted polyQ by the HMM
  proteinSet$LengthsPolyQ = rep(NA, nrow(proteinSet))
  proteinSet$AvgLengthsPolyQ = rep(NA, nrow(proteinSet))
  proteinSet$MaxLengthsPolyQ = rep(NA, nrow(proteinSet))
  
  # A column where each contains a list of the length of the region that is denoted polyQ by the HMM
  proteinSet$PolyQRegionQFractions = rep(NA, nrow(proteinSet))
  proteinSet$AvgPolyQRegionQFractions = rep(NA, nrow(proteinSet))
  proteinSet$MaxPolyQRegionQFractions = rep(NA, nrow(proteinSet))
  
  
  #proteinSet$binarySequence = rep(nrow(proteinSet))
  #proteinSet$Annotation = rep(nrow(proteinSet))
  for (sequenceNum in 1: nrow(proteinSet)){
    print(sequenceNum)
    # Debugging
      #sequenceNum = 3
    #
    
    proteinRow <- proteinSet[sequenceNum,]
    
    
    # Get Qfraction for the entire protein.
    nQs = str_count(proteinRow$peptide_seq, "Q")
    proteinRow$Qfraction = round(nQs/nchar(proteinRow$peptide_seq), 4)
    
    
    ### TURN TRANSLATION INTO BINARY VECTOR (coded as 1 if AA == Q, 0 otherwise)
    binarySequence = rep(1, nchar(proteinRow$peptide_seq))
    for (i in 1:length(binarySequence)){
      character = substr(proteinRow$peptide_seq, i, i)
      if (character == "Q"){
        binarySequence[i] = 1
      } else { binarySequence[i] = 0}
    }
    
    ### ADD BINARY SEQUENCE VECTOR TO HMM INPUTS
    HMM$x = binarySequence
    
    ### CREATING AN ANNOTATED (STATES) SEQUENCE FOR THE TRANSLATION
    annotatedSequence = Viterbi(HMM)

    ### DETERMINE IF THERE ARE ANY Q ISLANDS AND REPORT TRUE OR FALSE IN THE TABLE
    a = annotatedSequence == 2
    if (sum(a) >= 1 ) {
    } else {proteinRow}
    
    # Find which indices of the translation are marked as being polyQ
    polyQ_indices <- which(annotatedSequence %in% c(2))
    
    # Group the indicies into groups (so if you have multiple polyQ regions, the indices for each polyQ region will be seperated from each other)
    grouped_polyQ_indicies <- split(polyQ_indices, cumsum(c(1, diff(polyQ_indices) != 1))) 
    
    
    # Record how many polyQ regions there are
    if (length(grouped_polyQ_indicies[[1]]) == 0) {
      proteinRow$NumberPolyQ = 0
      proteinRow$HMMhasPolyQ = FALSE
    } else {
      proteinRow$NumberPolyQ = length(grouped_polyQ_indicies)
      proteinRow$HMMhasPolyQ = TRUE
      proteinRow$IndiciesPolyQ <- paste0(paste0(grouped_polyQ_indicies), collapse="; ")
      
      
      # Get the translations for each polyQ region and store them and their lengths. 
      PolyQTranslations <- character(length(proteinRow$NumberPolyQ)) # Intialize storage vectors
      PolyQLengths <- character(length(proteinRow$NumberPolyQ))
      PolyQRegionQFractions <- numeric(length(proteinRow$NumberPolyQ))
      
      for (i in 1:proteinRow$NumberPolyQ) {
        group <- grouped_polyQ_indicies[[i]] # Indexing into a list is done with double brackets
        
        PolyQTranslations[i] <- substr(proteinRow$peptide_seq, group[1], tail(group, 1))
        PolyQLengths[i] <- nchar(PolyQTranslations[i])
        
        
        nQs = str_count(PolyQTranslations[i], "Q")
        PolyQRegionQFractions[i] <- round(nQs/nchar(PolyQTranslations[i]), 4)
      }
      
      proteinRow$TranslationsPolyQ <- paste0(PolyQTranslations, collapse = "; ")
      
      PolyQLengths <- as.numeric(PolyQLengths)
      proteinRow$LengthsPolyQ <- paste0(sort(PolyQLengths, decreasing = TRUE), collapse = "; ")
      proteinRow$AvgLengthsPolyQ <- round(mean(PolyQLengths), 4)
      proteinRow$MaxLengthsPolyQ <- PolyQLengths[1]
      
      proteinRow$PolyQRegionQFractions <- paste0(sort(PolyQRegionQFractions, decreasing = TRUE), collapse = "; ")
      proteinRow$AvgPolyQRegionQFractions <- round(mean(PolyQRegionQFractions), 4)
      proteinRow$MaxPolyQRegionQFractions <- PolyQRegionQFractions[1]
    }
    
    # Add the row back
    proteinSet[sequenceNum,] <- proteinRow
  }
  
  # Return a column that gives length of largest region 
  # and the number of regions present (and the Qfraction of the largest region?)
  
  return (proteinSet)
}

########## End of helper functions ###############


# Diagnostics: Measure download time
start_time <- Sys.time() 

ensembl <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")
listFilters(ensembl)# Stuff we can filter by
listAttributes(ensembl) # stuff we can retrive


# Get protein sequences by ENTREZ ids
FLYBASEPROT <- names(as.list(org.Dm.egFLYBASEPROT2EG))
length(FLYBASEPROT)
length(unique(FLYBASEPROT))

my_df <- getBM(attributes = c('flybase_gene_id', 'ensembl_peptide_id', "external_gene_name", 'go_id', "name_1006", "namespace_1003"), 
           filters = "ensembl_peptide_id", 
           values = FLYBASEPROT, 
           mart = ensembl)
my_df <- rename(my_df, name_1006 = "go_description", namespace_1003 = "go_type")


# Paste together GO terms if they are for the same peptide and of the same GO_type.
my_df <- my_df %>%  
  distinct() %>%
  group_by(go_type, ensembl_peptide_id) %>% 
  mutate(go_list = paste(go_description, collapse = "; "), go_ids = paste(go_id, collapse = "; "))
my_df <- my_df %>% dplyr::select(-c(go_description, go_id)) 

vec <- duplicated(my_df[,c("go_type", "ensembl_peptide_id")]) # Get rows that share same GO type and peptide ID.
my_df <- my_df[!vec,]  # Remove duplicated rows.


# Do the filtering here.
seqs <- getSequence(id = FLYBASEPROT, 
                    type="ensembl_peptide_id",
                    seqType="peptide",
                    mart=ensembl) 
seqs <- rename(seqs , peptide = "peptide_seq")

markov_result <- HiddenMarkov(seqs)

my_df <- merge(my_df, markov_result, by = "ensembl_peptide_id")
my_df <- my_df %>% dplyr::select("flybase_gene_id", everything())

write.csv(my_df, file = paste0(output_dir, "fly_proteins_v2.csv"), row.names = FALSE)


############################
end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)
############################

########################################################3

# Load in protein set data from local storage if you didn't download the sets from above
  proteins <- read.csv(paste0(output_dir, "fly_proteins_v2.csv"))
#####


#########################################################
#########################################################

polyQ_proteins = proteins %>% filter(proteins$HMMhasPolyQ == TRUE)

proteins_unique <- unique(proteins$ensembl_peptide_id)
polyQ_proteins_unique <- unique(polyQ_proteins$ensembl_peptide_id)
length(proteins_unique) # 30284
length(polyQ_proteins_unique) # 9713

# Contains all genes, not just significant prpoteins
universe_df <- bitr(proteins_unique,
                    fromType = "FLYBASEPROT",
                    toType = c("ENTREZID", "UNIPROT"),
                    OrgDb = org.Dm.eg.db)


# Create df that contains only significant proteins (HMM proteins)

sig_df <- bitr(polyQ_proteins_unique,
               fromType = "FLYBASEPROT",
               toType = c("ENTREZID", "UNIPROT"),
               OrgDb = org.Dm.eg.db)



# ###########################################
# clusterProfiler Part 1 - GO
# ###########################################

for (go_ont in c("BP", "MF", "CC")) {
  print(paste0("Doing enrichment analysis for ", go_ont))
  
  
  # Control for sample size: Random subset of universe df with same size as
  # sig_df
  go_enrich_output <- enrichGO(gene = sig_df$ENTREZID,
                               universe = universe_df$ENTREZID,
                               OrgDb = org.Dm.eg.db,
                               ont = go_ont,
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.05,
                               readable = TRUE)
  
  control_df <- sample_n(universe_df, nrow(sig_df))
  go_enrich_output_control <- enrichGO(gene = control_df$ENTREZID,
                               universe = universe_df$ENTREZID,
                               OrgDb = org.Dm.eg.db,
                               ont = go_ont,
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.05,
                               readable = TRUE)
  
  go_result <- go_enrich_output@result
  go_control_result <- go_enrich_output_control@result
  if (nrow(go_enrich_output) != 0) {
    write.csv(go_enrich_output, file = paste0(output_dir, "Enriched_GO_", go_ont, ".csv"), row.names = FALSE)
    write.csv(go_enrich_output_control, file = paste0(output_dir, "Control_Enriched_GO_", go_ont, ".csv"), row.names = FALSE)
  }
  
  setwd(output_dir); getwd()
  # Control plots
  pdf(paste0("Control_Enriched_GO_", go_ont, "_plots.pdf"), width = 8.5, height = 11)
  par(mfrow = c(1, 1))
  p1 <- barplot(go_enrich_output_control, showCategory=8)
  p2 <- dotplot(go_enrich_output_control)
  p3 <- emapplot(go_enrich_output_control)
  p4 <- cnetplot(go_enrich_output_control, 
                 categorySize="pvalue", 
                 foldChange=control_df$ENTREZID,
                 node_label = FALSE,
                 colorEdge = TRUE)
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  dev.off()
  
  # Treatment
  pdf(paste0("Enriched_GO_", go_ont, "_plots.pdf"), width = 8.5, height = 11)
  par(mfrow = c(1, 1))
  p1 <- barplot(go_enrich_output, showCategory=8)
  p2 <- dotplot(go_enrich_output)
  p3 <- emapplot(go_enrich_output)
  p4 <- cnetplot(go_enrich_output, 
           categorySize="pvalue", 
           foldChange=sig_df$ENTREZID,
           node_label = FALSE,
           colorEdge = TRUE)
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  dev.off()
}



# ###########################################
# clusterProfiler Part 2 - KEGG


# kegg_enrich_output <- enrichKEGG(gene = sig_df$ENTREZID,
#                                  organism = "dme",
#                                  keyType = "ncbi-proteinid",
#                                  pvalueCutoff = 0.05)
# 
# kegg_result <- kegg_enrich_output@result
# if (nrow(go_enrich_output) != 0)
# {
#   write.csv(go_enrich_output, file = paste0(output_dir, group_name, "Enriched_KEGG.csv"), row.names = FALSE)
# }


# ###########################################

# KEGG ENRICHMENT



#########3 ClusterProfiler ##################
# Use clusterProfiler to compare proteins with polyQ to all proteins



############################################################################3
################################## Plotting ##################################
############################################################################
  
# Number of randomly chosen 50 amino acid peptide chunks [WORKS]

# Helper function
rnd_substr <- function(x, length) {
  if (nchar(x) < length) {return(x)}
  else {
    start = sample(1:(nchar(x)-length), 1, replace=T)
    end = start + (length-1)
    random_substring = substr(x, start, end)
    return(random_substring)
  }
}
###########################################
pdf(paste0(plot_output_dir, "HMM_result_plots_NP", nrow(proteins),".pdf"))
  

########################################################
#######################
# Plot type: Whole Protein Qfraction Histogram for all proteins
p <- ggplot(data = proteins, aes(x = proteins$Qfraction)) +
  geom_histogram(show.legend = FALSE, aes(fill=..count..)) +
  xlab("Whole Protein %Q") +
  ylab("# Proteins") +
  ggtitle("%Q Whole Protein") +
  theme_light()
print(p)
  
filtered_df <- proteins %>% filter(proteins$Qfraction > quantile(proteins$Qfraction, 0.95))
filtered_df <- filtered_df[-order(filtered_df$Qfraction),]
write.csv(filtered_df, file = paste0(output_dir, "Qfraction_all_wholeProts_95sig.csv"))


# Plot type: Whole Protein Qfraction Histogram split by HMMhasPolyQ
p <- ggplot(data = proteins, aes(x = Qfraction, fill = HMMhasPolyQ)) +
  geom_histogram(alpha=0.5, position="identity") +
  xlab("%Q (Whole Protein Seqs)") +
  ylab("# Proteins") +
  labs(fill = "Contains \n polyQ") +
  ggtitle("%Q Whole Proteins Histogram") +
  theme_light()
print(p)

p <- ggplot(data = proteins, aes(x = Qfraction, fill = HMMhasPolyQ)) +
  geom_density(alpha = 0.5) +
  xlab("%Q (Whole Protein Seqs)") +
  ylab("Density") +
  labs(fill = "Contains \n polyQ") +
  ggtitle("%Q Whole Proteins Density Plot") +
  theme_light()
print(p)

########################################################
# Plot type: 100000 50AA chunks, each randomly chosen from randomly chosen Protein seqs
df <- sample_n(proteins, 100000, replace = TRUE)
df$RandomSeqs <-unlist(lapply(df$Translations, rnd_substr, 50))
nQs = str_count(df$RandomSeqs, "Q")
df$QfractionRandomSeq <- nQs/nchar(df$RandomSeqs)


p <- ggplot(data = df, aes(x = df$QfractionRandomSeq)) +
  geom_histogram(show.legend = FALSE, aes(fill=..count..)) +
  xlab("%Q (Random 50AA Seqs)") +
  ylab("# of Randomly Choosen Pepitides") +
  ggtitle("%Q for 100000 50AA Random Seqs Histogram") +
  theme_light()
print(p)

filtered_df <- df %>% filter(df$QfractionRandomSeq > quantile(df$QfractionRandomSeq, 0.95))
filtered_df <- filtered_df[-order(filtered_df$QfractionRandomSeq),]
write.csv(filtered_df, file = paste0(output_dir, "Qfraction_randomAAs_95sig.csv"))

########################################################
# Plot type: Number of Proteins annotated by HMM to be "Poly Q" vs %Q

p <- ggplot(data = df, aes(x = df$AvgPolyQRegionQFractions)) +
  geom_histogram(show.legend = FALSE, aes(fill=..count..)) +
  xlab("Avg %Q of polyQ Regions") + 
  ggtitle("Avg %Q of polyQ Regions for polyQ Containing Proteins Histogram") +
  ylab("# Proteins Containing polyQ") +
  theme_light()

print(p)

df <- proteins %>% filter(HMMhasPolyQ == TRUE)
p <- ggplot(data = df, aes(x = df$MaxPolyQRegionQFractions)) +
  geom_histogram(show.legend = FALSE, aes(fill=..count..)) +
  xlab("Max %Q of polyQ Regions") +
  ggtitle("Max %Q of PolyQ Regions for polyQ Containing Proteins Histogram") +
  ylab("# Proteins Containing polyQ") +
  theme_light()

print(p)

##############
# Plot type 4: Number of Proteins annotated by HMM to be Poly Q vs Length of Region

p <- ggplot(data = df, aes(x = df$AvgLengthsPolyQ)) +
  geom_histogram(show.legend = FALSE, aes(fill=..count..)) +
  xlab("Average PolyQ Length") +
  ylab("# Proteins Containing polyQ (by HMM)") +
  ggtitle("Average PolyQ Length in PolyQ Containing Proteins Histogram") +
  theme_light()

print(p)


p <- ggplot(data = df, aes(x = df$MaxLengthsPolyQ)) +
  geom_histogram(show.legend = FALSE, aes(fill=..count..)) +
  xlab("Max PolyQ Length") +
  ylab("# Proteins Containing polyQ (by HMM)") +
  ggtitle("Max PolyQ Length in PolyQ Containing Proteins Histogram") +
  theme_light()

print(p)

dev.off()
############ End Plotting ################



#################### Misc. Code #################################


# ggplot(data = proteins, aes(x = proteins$Qfraction, y = ..count../sum(..count..))) +
#   geom_density() +
#   xlab("Q-fraction (Whole Protein Seqs)") +
#   ylab("Percent") +
#   theme_light()

# clusterProfiler needs entrenz id
#

# Add go Annotations

# head(listFilters(ensembl),20) # Stuff we can filter by
# head(listAttributes(ensembl), 20)# stuff we can retrive

# columns(org.Dm.eg.db)
# select(org.Dm.eg.db,
#        keys=as.character(protein_set$AccessionNumber), keytype="ACCNUM",
#        columns = c("GO", "SYMBOL", "ALIAS", "GENENAME", "UNIGENE", "UNIPROT"))

      
