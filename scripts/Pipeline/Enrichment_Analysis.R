# Daniel Alfonsetti
# daniel.alfonsetti@gmail.com; alfonset@mit.edu
# MIT Littleton Lab
# 7 March, 2019
# Description: Use enrichment analysis to find GO gene categories that are enriched for
# proteins that contain a polyAA track as annoated by our HMM.

#########################################
# Load libraries
#########################################
rm(list = ls())

library(doParallel)

# HMM tools
library(HiddenMarkov) 
library(seqHMM) 

# Data manipulation tools
library(TraMineR) 
library(stringr)
library(grid)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(data.table)

# Biological databases and their tools
library(biomaRt)
library(biomartr)
library(STRINGdb)
library(org.Dm.eg.db)
library(org.Hs.eg.db) 
library(org.Mm.eg.db)
library(org.Ce.eg.db)
library(GO.db)
library(protr)
library(clusterProfiler)

# Enrichment analysis tools
#########################################
# Setup environment
#########################################

registerDoParallel(cores = 4)

setwd("C:/UROPs/polyQ_neuronal_proteins")
plot_output_dir <- "C:/UROPs/polyQ_neuronal_proteins/plots/"
data_dir <- "C:/UROPs/polyQ_neuronal_proteins/data/"

output_base_dir <- "C:/UROPs/polyQ_neuronal_proteins/output/"


#########################################
# Helper function
#########################################
make_enrichment_plots_and_files <- function(go_enrich_output, type, file_dir_and_name) {
  if (!is.null(go_enrich_output)) {
    go_result <- go_enrich_output@result
    
    write.csv(go_enrich_output, file = paste0(file_dir_and_name, ".csv"), row.names = FALSE)
    pdf(paste0(file_dir_and_name, "_plots.pdf"), width = 8.5, height = 11)

    par(mfrow = c(1, 1))
    p1 <- barplot(go_enrich_output, showCategory=8)
    p2 <- dotplot(go_enrich_output)
    print(p1)
    print(p2)
    if (nrow(go_enrich_output) != 0) {
      p3 <- emapplot(go_enrich_output)
      print(p3)
      
      # p4 <- cnetplot(go_enrich_output, 
      #                categorySize="pvalue", 
      #                foldChange=sig_df$ENTREZID,
      #                node_label = FALSE,
      #                colorEdge = TRUE)
      # print(p4)
    }
    
    dev.off()
  } # if (!is.null(go_enrich_output)) {
} # make_enrichment_plots_and_files()
#########################################
# User defined variables
#########################################

species_vec = c("fly")

# candidate_AA_vec = c("Q", "D")

candidate_AA_vec = c("D", "T", "S", "E", "P", "G", "A", "C", "V", "M",
                     "I", "L", "Y", "F", "H", "K", "R", "W", "Q", "N")

models = c("adjusted", "trained")

  # debugging; species = species_vec[1]
for (model in models) {
  for (species in species_vec) {
    for (candidate_AA in candidate_AA_vec) {
      
    output_dir <- paste0(output_base_dir, model, "/", species, "/", candidate_AA, "/")
    dir.create(output_dir, recursive = TRUE)
    
    #########################################
    # ClusterProfiler module
    #########################################
    
    # Load in protein set data 
    

    proteins <- read.csv(paste0(output_dir, species, "_prots_w_HMM_", candidate_AA, "_nuclear_filt_transcriptome_filt.csv"))
    
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
    
    # Filter out nuclear proteins here using CC category
    # proteins <- 
    
    # Do enrichment for each subontology of GO
    #########################################
    # clusterProfiler Part 1 - GO
    #########################################
    
    for (go_ont in c("BP", "MF", "CC")) {
      # go_ont <- "BP"
      print(paste0("Doing enrichment analysis for ", go_ont, " ontology."))
      
      # Only allow proteins that have annotations in current ontology to be in background universe.
      proteins_has_ont_type <- proteins %>% filter(!is.na(eval(parse(text=paste0("go_ids_", go_ont)))))
      nrow(proteins); nrow(proteins_has_ont_type)
      
      
      # Get universe of proteins list
      proteins_vec <- proteins_has_ont_type$ensembl_peptide_id

      # Get significant proteins list
      polyAA_proteins_vec <- proteins_has_ont_type %>% 
                             filter(HMMhasPolyAA) %>%
                             dplyr::select(c(ensembl_peptide_id))[,1]
      
      
      # Get mapping from ensembl to entrez (entrez ids are needed for enrichment analysis)
      # Create vector of  entrez ids for significant genes.
      # Not there can be a 1 to many mapping from ensembl to entrez in the above step.
      # thus we need to use 'unique'.
      universe_df <- bitr(proteins_vec,
                          fromType = fromType,
                          toType = c("ENTREZID", "SYMBOL"),
                          OrgDb = OrgDb)
      universe_genes <- unique(universe_df$ENTREZID)
      

      sig_df <- bitr(polyAA_proteins_vec,
                     fromType = fromType,
                     toType = c("ENTREZID", "SYMBOL"),
                     OrgDb = OrgDb)
      sig_genes <- unique(sig_df$ENTREZID)
      
      # Create vector of entrez ids for a random set of genes that is the same size
      # of the set of significant genes in order to control for sample size. (FDR adjustment should do this,
      # but we are going to add this step just for kicks)
      # Create vector of random control genes with same number as significant genes.
      control_genes <- sample(universe_genes, length(sig_genes))
      # print(paste0("Length control genes: ", length(control_genes)))
      
      go_enrich_output_test <- enrichGO(gene = sig_genes,
                                        universe = universe_genes,
                                        OrgDb = OrgDb,
                                        ont = go_ont,
                                        pAdjustMethod = "BH",
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.05,
                                        readable = TRUE)
      
      # Control group
      go_enrich_output_control <- enrichGO(gene = control_genes,
                                           universe = universe_genes,
                                           OrgDb = OrgDb,
                                           ont = go_ont,
                                           pAdjustMethod = "BH",
                                           pvalueCutoff = 1,
                                           qvalueCutoff = 1,
                                           readable = TRUE)
      

      enrichment_output_dir <- paste0(paste0(output_dir, "/GO_Enrichment_", go_ont))
      dir.create(enrichment_output_dir)
      file_name = paste0(paste0(output_dir, "/GO_Enrichment_", go_ont, "/GO_", go_ont, "_results_CONTROL"))
      make_enrichment_plots_and_files(go_enrich_output_control, type = "control", file_dir_and_name = file_name)
      file_name = paste0(paste0(output_dir, "/GO_Enrichment_", go_ont, "/GO_", go_ont, "_results"))
      make_enrichment_plots_and_files(go_enrich_output_test, type = "test", file_name)
    } # for (go_ont in c("BP", "MF", "CC")) {
    
    #########################################
    # End of clusterProfiler Part 1 - GO
    #########################################
    
    #########################################
    # clusterProfiler Part 2 - KEGG
    #########################################
    
    
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
    #########################################
    # End of clusterProfiler Part 2 - KEGG
    #########################################
    
    
    #########################################
    # Plotting
    #########################################
    

    }
  }
}



###################################################################
# category_to_proteins = read.csv(text="GO_id,GO_desc,p.adjust,avg_peptide_length, number_proteins")
# # Iterate over each enriched category, seeing if the length of the proteins in that category are correlated
# # with the p-value of that category
# 
# enrich_results <- go_enrich_output@result %>% filter(p.adjust < 0.05) # Should we use a lower threshold?
# if (nrow(enrich_results) > 0) { # Make sure there is actually stuff to iterate over.
#   for (i in 1:nrow(enrich_results)) {
#     row <- enrich_results[i,]
#     cat_gene_syms <- strsplit(row$geneID, split = '/' ) # Gene symbols of genes in category.
#     
#     # Get the peptide ids of all the genes that were in this enriched category. 
#     peptide_ids <- with(universe_df[universe_df$SYMBOL %in% cat_gene_syms[[1]],], get(fromType))
#     # Now get their HMM annotations and other data.
#     df <- proteins[proteins$ensembl_peptide_id %in% peptide_ids,]
#     
#     # mean(df$peptide_length())
#     mean <- mean(unlist(lapply(as.character(df$peptide_seq), function(x){length(unlist(strsplit(x, "")))})))
#     number_proteins <- length(peptide_ids)
#     category_to_proteins[i,] <- list(row$ID, row$Description, row$p.adjust, mean, number_proteins)
#   }
#   
#   # Make Avg length of peptide in category vs adjusted p-value of categories.
#   # This is a control to make sure that categories with longer peptides aren't 
#   # more significant just based on the idea.
#   # Are synaptic categories below the line? Play with residuals to check.
#   
#   # Check if the slopes are significantly different? <- Do this
#   # Check if the sum of the residuals is significantly different?
#   # What do I do?
#   
#   # hist(category_to_proteins$number_proteins)
#   
#   # Make sure category has at least 20 proteins.
#   df <- category_to_proteins %>% filter(category_to_proteins$number_proteins > 20)
#   
#   
#   df$synaptic_cat = df$GO_id %in% synapse_cats # True or False, depending on whether or not the category is a synaptic category
#   
#   # Make two linear regression.
#   x = df$avg_peptide_length
#   y = df$p.adjust
#   fit <- lm(y ~ x)
#   
#   
#   # df$synapse <- ifelse(df$GO_id %in% synapse_cats, TRUE, FALSE)
#   p <- ggplot(df, aes(x = avg_peptide_length, y = p.adjust)) + 
#     geom_point(alpha = 0.5,  aes(size = number_proteins, color = synaptic_cat)) +
#     stat_smooth(method = "lm", col = "red") +
#     labs(title = paste("Avg Length of Peptides in Enriched Category VS p-value of Enriched Category \n",
#                        "Adj.R2 = ",signif(summary(fit)$adj.r.squared, 5),
#                        "Intercept =",signif(fit$coef[[1]],5 ),
#                        "  p.val =",signif(summary(fit)$coef[2,4], 5))) +
#     xlab("Average Peptide Length of Category") + 
#     ylab("Adjusted p-value of Category") +
#     scale_size_continuous(name = "Number of Proteins") +
#     scale_color_discrete(name = "Synaptic Category") +
#     theme_light() +
#     theme(plot.title = element_text(hjust = 0.5),
#           legend.title =) +
#     # geom_text(aes(label = GO_desc), color ='red', data = df[df$GO_id %in% synapse_cats,]) +
#     geom_text_repel(mapping = aes(label = GO_desc), 
#                     arrow = arrow(length = unit(0.03, "npc"), type = 'closed', ends ="first"),
#                     color ='#00CCCC', 
#                     data = df[df$synaptic_cat == TRUE,]) 
#   print(p)
# } # for (i in 1:nrow(enrich_results)) {
##############################################################################
