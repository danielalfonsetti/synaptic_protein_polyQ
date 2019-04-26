####################################################
# Daniel Alfonsetti, daniel.alfonsetti@gmail.com
# MIT, Littleton Lab UROP
# 27 February 2019
# ---------------------------------
# Description: A script with functionality to take protein sequences and corresponding HMM annotations for that sequence
# and generate a plot divided by amino acid, showing the position of each amino acid of a certain type while overlayed by
# polyAA tracks denoted for that amino acid by the HMM.
####################################################
# Import libraries and set 'global' variables
####################################################
rm(list = ls())
source("C:/UROPs/polyQ_neuronal_proteins/scripts/ConstantsAndFunctions.R")

library(ggplot2)
library(reshape2)
library(dplyr)

######################################################
# Helper functions
######################################################

# Helper functions for plotting
# Helper function 1
PolyAAChart <- function(row, graphHasHMMannots = FALSE)  {

  # Get indicies of each amino acid in the protein
  peptideSeq <- as.character(row$peptideSeq)
  aaPosList <- lapply(
    kCandidateAAs,
    function(AA) 
      unlist(lapply(
        strsplit(peptideSeq, ''), 
        function(peptideSeq) 
          which(peptideSeq == AA)
      ))
  )
  
  ########## Convert list to dataframe ############3
  ## Compute maximum length
  maxLength <- max(sapply(aaPosList, length))
  
  ## Add NA values to list elements
  aaPosList <- lapply(aaPosList, function(v) { c(v, rep(NA, maxLength-length(v)))})
  
  ## Cbind
  df <- do.call(cbind, aaPosList)
  colnames(df) <- kCandidateAAs
  
  # Reorganize into a form useful for plotting
  melted_df <- melt(df)[,c(-1)]
  
  box_width = 7
  p <- ggplot() +
    geom_segment(aes(x=melted_df$value-box_width/2, xend = melted_df$value+box_width/2, 
                     y=melted_df$Var2, yend = melted_df$Var2),
                 alpha = 1/box_width,
                 size = 8) +
    xlab("Peptide Position Index") +
    ylab("Amino Acid") +
    ggtitle(paste0("Amino acid chart for ", row$external_gene_name, " (", row$ensembl_peptide_id, ")")) +
    theme_light() +
    theme(legend.position = "none")

  if (graphHasHMMannots) {
    hmmAnnotationsList <- lapply(
      kCandidateAAs,
      function(AA) 
        unlist(lapply(strsplit(as.character(eval(parse(text = paste0("row$indiciesPoly", AA)))), split = "; ")[[1]], function(i) eval(parse(text = i))))
    )
    # Remove first and last index in each list (because we plot each point as a line from its middle position to -1 and +1).
    hmmAnnotationsList <- lapply(hmmAnnotationsList, function(vec){vec[2:(length(vec)-1)]})
    
    # Find the vector with the max length, and pad the rest of the vectors so that we can convert the list to the dataframe.
    maxLength <- max(sapply(hmmAnnotationsList, length))
    hmmAnnotationsList <- lapply(hmmAnnotationsList, function(v) {if (is.na(v[1])) {rep(NA, maxLength)} else {c(v, rep(NA, maxLength-length(v)))}})
    
    df <- do.call(cbind, hmmAnnotationsList)
    colnames(df) <- kCandidateAAs
    meltedHmmAnnotsDf <- melt(df)[,c(-1)]
    if (!all(is.na(meltedHmmAnnotsDf$value))) {
      p <- p + geom_segment(aes(x=meltedHmmAnnotsDf$value-1, xend = meltedHmmAnnotsDf$value+1, 
                                  y=meltedHmmAnnotsDf$Var2, yend = meltedHmmAnnotsDf$Var2),
                              alpha = 1, color = "Green", size = 2.5) 
    }
  } #   if (graphHasHMMannots) {

  print(p)
}

# Helper function 2
PolyAAChartWrapper <- function(proteins){
  if (nrow(proteins) == 0){
    plot(NA, xlim=c(0,2), ylim=c(0,2), bty='n',
         xaxt='n', yaxt='n', xlab='', ylab='')
    text(1, 1.7,"No data", cex = 4)
    text(1, 1.4, "No proteins in this set had polyAA regions", cex = 1.5)
    text(1, 1.2, "Cannot make amino acid plots", cex = 1.5)
  } else {
    for (i in 1:nrow(proteins)){
      row <- proteins[i,]
      PolyAAChart(row, TRUE)
    }
  }
}


# Helper function 3
# Function to check if any two intervals overlap 
IsOverlap <- function(L) { 
  
  # Sort intervals in increasing order of start time 
  L = L[order(as.integer(unlist(lapply(L, function(x) x[1]))), decreasing=FALSE)]
  
  # In the sorted array, if start time of an interval 
  # is less than end of previous interval, then there 
  # is an overlap
  for (i in 2:length(L)) { 
    if (is.na(L[i-1]) | is.na(L[i])) {next}
    if (as.integer(L[[i-1]][2]) >= as.integer(L[[i]][1])){
      return(TRUE)
    }
  }
  
  #If we reach here, then no overlap
  return(FALSE)
} 
####################################################
# Preliminary Data Wrangling
####################################################

# Merge the polyAA output files for each HMM type.
for (hmmType in kModels) {
  # Important columns to take from each file are...
  # First 11 columns (general protein information, not polyAA annotation specific)
  # And 15th column (indices of polyAA regions)
  
  mergedPolyAaDf <- read.csv(paste0("C:/UROPs/polyQ_neuronal_proteins/output/",  hmmType, "/fly/", kCandidateAAs[1], "/fly_prots_w_HMM_", kCandidateAAs[1], ".csv"))
  mergedPolyAaDf <- mergedPolyAaDf[, c(1:11, 16)]
  new_name <- paste0("indiciesPoly", kCandidateAAs[1])
  colnames(mergedPolyAaDf)[colnames(mergedPolyAaDf)=="indiciesPolyAA"] <- new_name
  
  # Iterate over the rest of the AAs and append their results.
  for (candidateAA in kCandidateAAs[2:length(kCandidateAAs)]) {
    df <- read.csv(paste0("C:/UROPs/polyQ_neuronal_proteins/output/", hmmType, "/fly/", candidateAA, "/fly_prots_w_HMM_", candidateAA, ".csv"))
    if ("FBpp0070031" %in% df$ensembl_peptide_id) {
      print(candidateAA)
    }
    new_name <- paste0("indiciesPoly", candidateAA)
    colnames(df)[colnames(df)=="indiciesPolyAA"] <- new_name
    mergedPolyAaDf = merge(x = mergedPolyAaDf, y = df[, c("ensembl_peptide_id",new_name)], by = "ensembl_peptide_id")
  }
  write.csv(mergedPolyAaDf, paste0("C:/UROPs/polyQ_neuronal_proteins/output/", hmmType, "/fly/mergedPolyAaDf.csv"), row.names = FALSE)
}

####################################################
# Plotting
####################################################


####################################################
# Plotting 1 - Explore the AA charts for some hand chosen proteins
####################################################

for (hmmType in kModels) {
  proteins <- read.csv(paste0("C:/UROPs/polyQ_neuronal_proteins/output/", hmmType, "/fly/mergedPolyAaDf.csv"))
  
  dir.create(paste0("C:/UROPs/polyQ_neuronal_proteins/output/",hmmType,"/fly/AA_charts"))
  setwd(paste0("C:/UROPs/polyQ_neuronal_proteins/output/",hmmType,"/fly/AA_charts"))
  
  # Known AZ proteins
  AZproteins <- proteins %>% filter(proteins$external_gene_name %in% c("brp", "Rbp", "Rim"))
  AZproteins <- AZproteins[!duplicated(AZproteins$ensembl_gene_id),]
  
  pdf(file = paste0("AA_charts_manualAZ_proteins_Model-", hmmType, "_fly.pdf"))
  PolyAAChartWrapper(AZproteins)
  dev.off()
  
  # Known PSD proteins
  # "homer", "Grip", "Dlg", "Sh3", "shk"
  psdProteins <- proteins %>% filter(proteins$external_gene_name %in% c("homer", "Grip", "dlg1", "SH3PX1", "Prosap"))
  psdProteins <- psdProteins[!duplicated(psdProteins$ensembl_gene_id),]

  pdf(file = paste0("AA_charts_manualPSD_proteins_Model-", hmmType, "_fly.pdf"))
  PolyAAChartWrapper(psdProteins)
  dev.off()
  
  # Proteins used in HMM training set
  training_set_protein_ids <- kTrainingSetMap[["flyQ"]]
  training_set_proteins <- proteins %>% filter(proteins$ensembl_peptide_id %in% training_set_protein_ids)
  
  pdf(file = paste0("AA_charts_training_set_proteins_Model-", hmmType, "_fly.pdf"))
  PolyAAChartWrapper(training_set_proteins)
  dev.off()
}

####################################################
# Plotting 2 - Explore AA charts for proteins of certain GO categories
####################################################

kLocations = list("AZ" = kAZcats, "Synapse" =  kSynapseCats, "PSD" = kPSDcats)

for (hmmType in kModels) {
  dir.create(paste0("C:/UROPs/polyQ_neuronal_proteins/output/",hmmType,"/fly/AA_charts"))
  setwd(paste0("C:/UROPs/polyQ_neuronal_proteins/output/",hmmType,"/fly/AA_charts/"))
  
  for (location in names(kLocations)) {
    proteins <- read.csv(paste0("C:/UROPs/polyQ_neuronal_proteins/output/", hmmType, "/fly/mergedPolyAaDf.csv"))
    locationCats <- kLocations[[location]]
    proteins <- AnnotateByCategorySet(proteins, locationCats)
    proteins <- proteins %>% filter(proteins$inSet)
    
    pdf(file = paste0("AA_charts__location-", location, "_proteins__Model-", hmmType, "__fly.pdf"))
    PolyAAChartWrapper(proteins)
    dev.off()
  } #   for (location in names(kLocations)) {
} # for (hmmType in kModels) {

####################################################
# Plotting 3 -   Plot neuronally expressed proteins (based off transcriptome) that have polyAAs
####################################################
for (hmmType in kModels) {
  dir.create(paste0("C:/UROPs/polyQ_neuronal_proteins/output/", hmmType, "/fly/AA_charts"))
  setwd(paste0("C:/UROPs/polyQ_neuronal_proteins/output/", hmmType, "/fly/AA_charts/"))
  
  proteins <- read.csv(paste0("C:/UROPs/polyQ_neuronal_proteins/output/", hmmType, "/fly/mergedPolyAaDf.csv"))
  # kNeuronalTranscripts <- as.vector(read.table("C:/UROPs/polyQ_neuronal_proteins/output/fly_CNS_transcriptome_mh-l.txt", sep = "\t"))
  proteins <- proteins[proteins$ensembl_peptide_id %in% kNeuronalTranscripts$V1,]
  proteins <- proteins %>% filter(!is.na(proteins$indiciesPolyQ)) # Only plot proteins that have a polyQ region
  
  # Remove nucleus related and transcription factor proteins
  proteins <- AnnotateByCategorySet(proteins, kNuclearCats)
  proteins <- proteins %>% filter(!proteins$inSet)
  
  pdf(file = paste0("AA_charts__location-neuronal_transcriptome_proteins_(non_nuclear)__Model-", hmmType, "__fly.pdf"))
  PolyAAChartWrapper(proteins)
  dev.off()
} # for (hmmType in kModels) {


####################################################
# Plotting 4 - Plot all proteins have overlapping polyAA regions
####################################################
# Check for interval overlaps
for (hmmType in kModels) {
  dir.create(paste0("C:/UROPs/polyQ_neuronal_proteins/output/", hmmType, "/fly/AA_charts"))
  setwd(paste0("C:/UROPs/polyQ_neuronal_proteins/output/", hmmType, "/fly/AA_charts/"))
  
  proteins <- read.csv(paste0("C:/UROPs/polyQ_neuronal_proteins/output/", hmmType, "/fly/mergedPolyAaDf.csv"))

  intervalData <- proteins[,(ncol(proteins)-20+1):ncol(proteins)] # Only get the columns with polyAA intervals

  intervalsList <- apply(intervalData, 1, function(row)
    unlist(sapply(row, function(x)
      lapply(strsplit(as.character(x), "; "), function(y)
        strsplit(as.character(y), ":"))), 
      recursive = FALSE))
  
  # Find which proteins have overlapping polyAA regions
  res <- unlist(lapply(intervalsList, IsOverlap))
  proteins <- proteins[res,]
  
  # Only plot proteins that have overlapping polyAA regions.
  pdf(file = paste0("AA_charts_overlappingPolyAAs__Model-", hmmType, "__fly.pdf"))
  PolyAAChartWrapper(proteins)
  dev.off()
  
  ####################################################
  # Plot proteins with overlapping polyAAs that ALSO have neuronally expressed proteins (based off transcriptome)
  # kNeuronalTranscripts <- as.vector(read.table("C:/UROPs/polyQ_neuronal_proteins/output/fly_CNS_transcriptome_mh-l.txt", sep = "\t"))
  proteins <- proteins[proteins$ensembl_peptide_id %in% kNeuronalTranscripts$V1,]
  proteins <- proteins %>% filter(!is.na(proteins$indiciesPolyQ)) # Only plot proteins that have a polyQ region
  
  # Remove nucleus related and transcription factor proteins
  proteins <- AnnotateByCategorySet(proteins, kNuclearCats)
  proteins <- proteins %>% filter(!proteins$inSet)
  
  pdf(file = paste0("AA_charts_overlappingPolyAAs__location-neuronal_transcriptome_proteins_(non_nuclear)__Model-", hmmType, "__fly.pdf"))
  PolyAAChartWrapper(proteins)
  dev.off()
} # for (hmmType in kModels) {
# END
