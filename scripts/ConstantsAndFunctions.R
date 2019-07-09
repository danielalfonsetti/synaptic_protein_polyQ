####################################################
# Daniel Alfonsetti, daniel.alfonsetti@gmail.com
# MIT, Littleton Lab UROP
# 17 April 2019
# ---------------------------------
# Description: Some helper functions and constants that are used in multiple scripts in this project.
####################################################
print("Running ConstantsAndFunctions.R!")

# Remove everything in environment
rm(list = ls()) 

# Detach all libraries
# detach_helper <- function(x){
#   if (x %in% search()){
#     detach(x, character.only=TRUE,unload=TRUE)
#   }
# }
# invisible(lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""), detach_helper))

Packages <- c("dplyr", "GO.db", "org.Dm.eg.db", "org.Hs.eg.db", "org.Mm.eg.db", "org.Ce.eg.db")
invisible(lapply(Packages, library, character.only = TRUE))

####################################################3
# CONSTANTS
####################################################

# Default parameters
kModels = c("adjusted")
kSpecies = c("fly")
kCandidateAAs= c("D", "T", "S", "E", "P", "G", "A", "C", "V", "M",
                 "I", "L", "Y", "F", "H", "K", "R", "W", "Q", "N")


kTrainingSetMap <- list("flyQ"= c("FBpp0289769","FBpp0307700", "FBpp0086727", "FBpp0111724", "FBpp0293366",
                                  "FBpp0070830", "FBpp0309352", "FBpp0402897", "FBpp0110299","FBpp0305807"),
                        "mouseQ" = NA,
                        "wormQ" =  NA)


kNeuronalTranscripts <- as.vector(read.table("../data/fly_CNS_transcriptome_mh-l.txt", sep = "\t"))
kFlyNeuronalTranscripts <- kNeuronalTranscripts

# TODO: Download these
kWormNeuronalTranscripts <- NA 
kHumanNeuronalTranscripts <- NA
kMouseNeuronalTranscripts <- NA 


# Transcription factor categories
# "GO:0003676" is nucleic acid binding
kTranscriptionCats <- c("GO:0003676", get("GO:0003676", GOMFOFFSPRING))

# Nucleus related categories
# "GO:0009295" is nucleoid
# "GO:0005634"  is nucleus 
kNuclearCats <- c("GO:0005634", get("GO:0005634", GOCCOFFSPRING), 
                  "GO:0009295", get("GO:0009295", GOCCOFFSPRING),
                  kTranscriptionCats)

# Synaptic categories
kSynapseCats <- c("GO:0045202", get("GO:0045202", GOCCOFFSPRING))

# Active Zone categories
kAZcats <- c("GO:0048786", get("GO:0048786", GOCCOFFSPRING))

# Post synaptic density categories
kPSDcats <- c("GO:0014069", get("GO:0014069", GOCCOFFSPRING))

####################################################
# FUNCTIONS
####################################################
AnnotateByCategorySet <- function(df, cats) {
  CC <- unlist(
    lapply(
      lapply(as.character(df$go_ids_CC), function(x){unlist(strsplit(x, "; "))}), 
      function(x){any(x %in% cats)}
    )
  )
  MF <- unlist(
    lapply(
      lapply(as.character(df$go_ids_MF), function(x){unlist(strsplit(x, "; "))}), 
      function(x){any(x %in% cats)}
    )
  )
  BP <- unlist(
    lapply(
      lapply(as.character(df$go_ids_BP), function(x){unlist(strsplit(x, "; "))}), 
      function(x){any(x %in% cats)}
    )
  )
  df$inSet <- CC | MF | BP
  return(df)
}
print("Completed ConstantsAndFunctions.R")
# END


# Helper functions for plotting amino acid sequences. Used in "AminoAcidCharts.R" and "ExperimentProteins.R"
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

