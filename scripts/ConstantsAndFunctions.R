####################################################
# Daniel Alfonsetti, daniel.alfonsetti@gmail.com
# MIT, Littleton Lab UROP
# 17 April 2019
# ---------------------------------
# Description: Some helper functions and constants that are used in multiple scripts in this project.
####################################################

# Remove everything in environment
rm(list = ls()) 
# Detach all libraries
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE) 


.libPaths(c("C:/Users/danie/Documents/R/win-library/3.5",.libPaths()))
library(dplyr)
library(GO.db)

library(org.Dm.eg.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Ce.eg.db)

####################################################3
# CONSTANTS
####################################################

# Default parameters
kModels = c("adjusted")
kSpecies = c("fly")
kCandidateAAs= c("D", "T", "S", "E", "P", "G", "A", "C", "V", "M",
                 "I", "L", "Y", "F", "H", "K", "R", "W", "Q", "N")
# kCandidateAAs= c("Q")


kTrainingSetMap <- list("flyQ"= c("FBpp0289769","FBpp0307700", "FBpp0086727", "FBpp0111724", "FBpp0293366",
                                  "FBpp0070830", "FBpp0309352", "FBpp0402897", "FBpp0110299","FBpp0305807"),
                        "mouseQ" = NA,
                        "wormQ" =  NA)


kOutputBaseDir <- "C:/UROPs/polyQ_neuronal_proteins/output/"
kNeuronalTranscripts <- as.vector(read.table("C:/UROPs/polyQ_neuronal_proteins/output/fly_CNS_transcriptome_mh-l.txt", sep = "\t"))

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
# END
