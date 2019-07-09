# Daniel Alfonsetti
# daniel.alfonsetti@gmail.com; alfonset@mit.edu
# MIT Littleton Lab
# 4 April, 2019
# Description: Use enrichment analysis to find GO gene categories that are enriched for
# proteins that contain a polyAA track as annoated by our HMM.

#########################################
# Load libraries
#########################################
rm(list = ls())
kOutputBaseDir <- "../../../output/"
source("../../ConstantsAndFunctions.R", chdir=T)
library(clusterProfiler)
#########################################
# Setup environment
#########################################

setwd("C:/UROPs/polyQ_neuronal_proteins")

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 0) {
  kModels =  unlist(strsplit(args[1],","))
  kSpecies = unlist(strsplit(args[2],","))
  kCandidateAAs = strsplit(args[3],"")[[1]]
}
#########################################
# Helper function(s)
#########################################
MakeEnrichmentPlots <- function(goEnrichOutput, type, fileName) {
  if (!is.null(goEnrichOutput)) {
    
    write.csv(goEnrichOutput, file = paste0(fileName, ".csv"), row.names = FALSE)
    pdf(paste0(fileName, "_plots.pdf"), width = 8.5, height = 11)
    
    if (nrow(goEnrichOutput) != 0) {
      par(mfrow = c(1, 1))
      p1 <- barplot(goEnrichOutput, showCategory=8)
      p2 <- dotplot(goEnrichOutput)
      print(p1)
      print(p2)
      if (nrow(goEnrichOutput) != 0) {
        p3 <- emapplot(goEnrichOutput)
        print(p3)
      }
    } else {
      plot(NA, xlim=c(0,2), ylim=c(0,2), bty='n',
           xaxt='n', yaxt='n', xlab='', ylab='')
      text(1,1.7,"No data", cex = 4)
      text(1, 1.4, "No significant categories were found", cex = 2)
      text(1, 1.2, "Cannot make enrichment plots", cex = 2)
      
    } # if (nrow(goEnrichOutput) != 0) {

    dev.off()
  } # if (!is.null(goEnrichOutput)) {
} # MakeEnrichmentPlots()

#########################################
# Main
#########################################
for (model in kModels) {
  for (species in kSpecies) {
    for (candidateAA in kCandidateAAs) {
      outputDir <- paste0(kOutputBaseDir, model, "/", species, "/", candidateAA, "/")
      dir.create(outputDir, recursive = TRUE)
      
      # Load in protein set data 
      proteins <- read.csv(paste0(outputDir, species, "_prots_w_HMM_", candidateAA, "_nuclear_filt_transcriptome_filt.csv"))
      
      if (species == "fly") {
        orgDb = org.Dm.eg.db
        fromType = "FLYBASEPROT"
      } else if (species == "human") {
        orgDb = org.Hs.eg.db
        fromType = "ENSEMBLPROT"
      } else if (species == "mouse") {
        orgDb = org.Mm.eg.db
        fromType = "ENSEMBLPROT"
      } else if (species == "worm") {
        orgDb = org.Ce.eg.db
        fromType = "ENSEMBLPROT"
      }
      
      # Do enrichment for each subontology of GO
      #########################################
      # clusterProfiler Part 1 - GO
      #########################################
      
      for (goOnt in c("BP", "MF", "CC")) {
        # goOnt <- "BP"
        print(paste0("Doing enrichment analysis for ", goOnt, " ontology."))
        
        # Only allow proteins that have annotations in current ontology to be in background universe.
        proteinsHasOntType <- proteins %>% filter(!is.na(eval(parse(text=paste0("go_ids_", goOnt)))))
        nrow(proteins); nrow(proteinsHasOntType)
        
        # Get universe of proteins list
        proteinsVec <- proteinsHasOntType$ensembl_peptide_id
  
        # Get significant proteins list
        polyAAProteinsVec <- proteinsHasOntType %>% 
                               filter(hmmHasPolyAA) %>%
                               dplyr::select(c(ensembl_peptide_id))
        
        polyAAProteinsVec <- polyAAProteinsVec$ensembl_peptide_id
        # Get mapping from ensembl to entrez (entrez ids are needed for enrichment analysis)
        # Create vector of  entrez ids for significant genes.
        # Not there can be a 1 to many mapping from ensembl to entrez in the above step.
        # thus we need to use 'unique'.
        universeDf <- bitr(proteinsVec,
                            fromType = fromType,
                            toType = c("ENTREZID", "SYMBOL"),
                            OrgDb = orgDb)
        universeGenes <- unique(universeDf$ENTREZID)
  
        sigDf <- bitr(polyAAProteinsVec,
                       fromType = fromType,
                       toType = c("ENTREZID", "SYMBOL"),
                       OrgDb = orgDb)
        sigGenes <- unique(sigDf$ENTREZID)
        
        # Create vector of entrez ids for a random set of genes that is the same size
        # of the set of significant genes in order to control for sample size. (FDR adjustment should do this,
        # but we are going to add this step just for kicks)
        # Create vector of random control genes with same number as significant genes.
        controlGenes <- sample(universeGenes, length(sigGenes))
        
        # print(paste0("Length control genes: ", length(controlGenes)))
        goEnrichOutputTest <- enrichGO(gene = sigGenes,
                                          universe = universeGenes,
                                          OrgDb = orgDb,
                                          ont = goOnt,
                                          pAdjustMethod = "BH",
                                          pvalueCutoff = 0.05,
                                          qvalueCutoff = 0.05,
                                          readable = TRUE)
        
        # Control group
        goEnrichOutputControl <- enrichGO(gene = controlGenes,
                                             universe = universeGenes,
                                             OrgDb = orgDb,
                                             ont = goOnt,
                                             pAdjustMethod = "BH",
                                             pvalueCutoff = 1,
                                             qvalueCutoff = 1,
                                             readable = TRUE)
        
        enrichmentOutputDir <- paste0(paste0(outputDir, "/GO_Enrichment_", goOnt))
        dir.create(enrichmentOutputDir)
        
        fileName = paste0(paste0(outputDir, "/GO_Enrichment_", goOnt, "/GO_", goOnt, "_results_CONTROL"))
        MakeEnrichmentPlots(goEnrichOutputControl, type = "control", fileName = fileName)
        
        fileName = paste0(paste0(outputDir, "/GO_Enrichment_", goOnt, "/GO_", goOnt, "_results"))
        MakeEnrichmentPlots(goEnrichOutputTest, type = "test", fileName = fileName)
        
      } # for (goOnt in c("BP", "MF", "CC")) {
    } # for (candidateAA in kCandidateAAs) {
  } # for (species in kSpecies) {
} # for (model in kModels) {
# END