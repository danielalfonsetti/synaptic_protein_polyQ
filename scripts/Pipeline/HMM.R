# Daniel Alfonsetti
# daniel.alfonsetti@gmail.com; alfonset@mit.edu
# MIT Littleton Lab
# April 4, 2019
# Description: Run HMMs on protein sequences to identify poly amino acid (polyAA) tracts. 
# Run for each of the 20 different amino acids

# Identifier Conventions: 
#   Functions: MyFunction(...)
#   Constants: kMyConstant
#   Variables: myVariables
#   https://google.github.io/styleguide/Rguide.xml

#########################################
# Load libraries
#########################################
rm(list = ls())
source("C:/UROPs/polyQ_neuronal_proteins/scripts/ConstantsAndFunctions.R")

# Efficiency tools
library(doParallel)

# HMM tools
library(HiddenMarkov)
library(seqHMM) 

# Data manipulation tools
library(dplyr)
library(stringr)

#########################################
# Set constants
#########################################
registerDoParallel(cores = 4)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 0) {
  kModels =  unlist(strsplit(args[1],","))
  kSpecies = unlist(strsplit(args[2],","))
  kCandidateAAs = strsplit(args[3],"")[[1]]
}


#########################################
# Helper Function(s)
#########################################
FormatForHmm = function(proteinDf, candidateAA) {
  ########################
  # Get data into right form for HMM
  ########################
  print("At start of FormatForHmm!")

  listSeqsSpilt <- list()
  for (i in 1:nrow(proteinDf)) {
    row <- proteinDf[i,]
    split_seq_AAs <- unlist(strsplit(row$peptideSeq, split = ""))

    split_seq <- rep(NA, length(split_seq_AAs))

    for (j in 1:length(split_seq)) {
      if (split_seq_AAs[j] == candidateAA) {
        split_seq[j] <- "TargetAA"
      } else {
        split_seq[j] <- "NonTarget"
      }
    }
    listSeqsSpilt[[length(listSeqsSpilt)+1]] <- split_seq
  }
  print("In FormatForHmm: Split the seqs")

  # Pad the list with NAs
  maxlen <- max(lengths(listSeqsSpilt))
  listSeqsSpilt_padded <- lapply(listSeqsSpilt, `length<-`, maxlen)
  print("In FormatForHmm: Padded the seqs")

  # Conver to dataframe
  split_seq_result_df <- as.data.frame(do.call(rbind, listSeqsSpilt_padded), stringsAsFactors= FALSE)
  print("In FormatForHmm: Converted list of seqs into a dataframe!")

  # Convert to a seqdef object, part of the seqHMM API.
  seqs_object <- seqdef(split_seq_result_df)
  print("In FormatForHmm: Converted seq dataframe to seqs_object (seqHMM package object)!")

  print("Finished FormatForHmm!")
  return(seqs_object)
} # FormatForHmm = function(proteinDf, candidateAA)

TrainHMM = function(trainingSet, candidateAA) {
  print("In TrainHMM!")

  trainingSeqs <- FormatForHmm(trainingSet, candidateAA)

  # Make the model
  # Initial values for emission matrices
  emissProbs <- matrix(c(0.95, 0.05,
                    0.26, 0.74),
                  nrow=2, ncol=2, byrow = TRUE)

  # Initial values for transition matrix
  transProbs <- matrix(c(0.95, 0.05,
                    0.01, 0.99),
                  nrow = 2, ncol = 2, byrow = TRUE)


  # Initial values for initial state probabilities
  initialProbs <- c(.99, .01)

  hmm <- build_hmm(observations = trainingSeqs,
                   initial_probs = initialProbs,
                   transition_probs = transProbs,
                   emission_probs = emissProbs,
                   state_names = c("NonPolyAA", "PolyAA"))

  # Fit/tune the model
  hmmFit <- fit_model(hmm,
                       em_step = TRUE,
                       global_step = TRUE,
                       local_step = FALSE,
                       control_em = list(maxeval = 5, print_level = 2),
                       control_global = list(maxtime = 1),
                       threads = 10)

  print("Finished TrainHMM!")
  return(hmmFit)
} # TrainHMM = function(trainingSet, candidateAA)

TestHMM = function(proteinDf, model, candidateAA){
  print("At start of TestHMM!")

  # Remoeve any proteins with a sequence of length 1
  proteinDf <- proteinDf[nchar(as.character(proteinDf$peptideSeq))>1,]

  # Create model using train model parameters. The 'x' variable while
  # be updated with a single sequence each time we iterate over a new row of
  # proteinDf.
  x <- NULL
  HMM <- dthmm(x,
               Pi = model$transition_probs,  # Looks right
               delta = model$initial_probs,  # Looks right
               distn = "binom",
               pm = list(prob = model$emission_probs[1:2, 2], size = 1),       #
               pn = NULL,
               discrete = TRUE,
               nonstat = TRUE)
  print("In TestHMM: Finished building the object!")

  # Columns for fractions of Single, double, and triple AAs of the particular type in the entire protein sequence.
  proteinDf$aaComp = rep(NA, nrow(proteinDf))
  # proteinDf$DipeptideComp = rep(NA, nrow(proteinDf))
  # proteinDf$TripeptideComp = rep(NA, nrow(proteinDf))

  # A column where each cell contains a bboolean indicating whether or not at least one polyQ region was found by the HMM for this protein.
  proteinDf$hmmHasPolyAA = rep(NA, nrow(proteinDf))

  # A column where each cell is a number indicating the number of polyQ regions found on this protein by the HMM.
  proteinDf$numPolyAA = rep(NA, nrow(proteinDf))

  # A column where each cell contains a list with the indicies for each region in the protein that is denoted polyQ by the HMM
  proteinDf$indiciesPolyAA = rep(NA, nrow(proteinDf))

  # A column where each cell contains a list with the AA sequence for each region in the protein that is denoted polyQ by the HMM
  proteinDf$translationsPolyAA = rep(NA, nrow(proteinDf))

  # A column where each cell contains a list of the lengths of the regions that are denoted polyQ by the HMM
  proteinDf$lengthPolyAA = rep(NA, nrow(proteinDf))

  # Average and max lengths of polyQ regions found for a given protein.
  proteinDf$avgLengthPolyAA = rep(NA, nrow(proteinDf))
  proteinDf$maxLengthPolyAA = rep(NA, nrow(proteinDf))

  # A column where each cell contains a list of the fraction of target AAs in each polyAA sequence found for a given protein.
  proteinDf$polyAARegionAAFractions = rep(NA, nrow(proteinDf))

  # Average and Max of the fractions of target AAs in each polyAA sequence found for a given protein.
  proteinDf$avgPolyAARegionAAFractions = rep(NA, nrow(proteinDf))
  proteinDf$maxPolyAARegionAAFractions = rep(NA, nrow(proteinDf))

  # A column where each cell contains a sequence of 1s and 0s of the same length of the protein
  # where strings of 1s represent polyAA regions.
  proteinDf$annotatedSequence = rep(NA, nrow(proteinDf))

  # proteinDf$pvalues = rep(0, nrow(proteinDf))

  # TODO: Add amino acid descriptor, dipeptide descriptor, and tri peptide descriptor.
  # Add autocorrelatoin descriptor

  print("In TestHMM: Extracting information about annotated paths!")
  for (sequenceNum in 1:nrow(proteinDf)) {

    if (sequenceNum %% 5000 == 0){
      print(paste0("Annotating sequence ", sequenceNum, " out of ", nrow(proteinDf)))
    }
    proteinRow <- proteinDf[sequenceNum,]

    ### TURN TRANSLATION INTO BINARY VECTOR (coded as 1 if AA == Q, 0 otherwise)
    binarySequence = rep(1, nchar(proteinRow$peptideSeq))
    for (i in 1:length(binarySequence)) {
      character = substr(proteinRow$peptideSeq, i, i)
      if (character == candidateAA) {
        binarySequence[i] = 1
      } else {
        binarySequence[i] = 0
      }
    }

    ### ADD BINARY SEQUENCE VECTOR TO HMM INPUTS
    HMM$x = binarySequence

    ### CREATING AN ANNOTATED (STATES) SEQUENCE FOR THE TRANSLATION
    annotatedSequence = Viterbi(HMM)
    proteinRow$annotatedSequence <- paste0(annotatedSequence, collapse = "")

    # Get Qfraction for the entire protein.


    # proteinRow$aaComp <- extractAAC(proteinRow$peptideSeq)[paste0(replicate(1, candidateAA), collapse = "")]
    # proteinRow$DipeptideComp <- extractDC(proteinRow$peptideSeq)[paste0(replicate(2, candidateAA), collapse = "")]
    # proteinRow$TripeptideComp <- extractTC(proteinRow$peptideSeq)[paste0(replicate(3, candidateAA), collapse = "")]
    nAAs = str_count(proteinRow$peptideSeq, candidateAA)
    proteinRow$aaComp= round(nAAs/nchar(proteinRow$peptideSeq), 4)

    # Find which indices of the translation are marked as being polyQ
    polyAA_indices <- which(annotatedSequence %in% c(2))

    # Group the indicies into groups (so if you have multiple polyQ regions, the indices for each polyQ region will be seperated from each other)
    groupedPolyAAIndicies <- split(polyAA_indices, cumsum(c(1, diff(polyAA_indices) != 1)))

    # Only keep polyAAs that are greater than 6AAs in length here
    groupedPolyAAIndicies <- groupedPolyAAIndicies[lapply(groupedPolyAAIndicies, length) > 6]

    # Record how many polyQ regions there are
    if (length(groupedPolyAAIndicies) == 0) {
      proteinRow$numPolyAA = 0
      proteinRow$hmmHasPolyAA = FALSE
    } else {
      proteinRow$numPolyAA = length(groupedPolyAAIndicies)
      proteinRow$hmmHasPolyAA = TRUE
      proteinRow$indiciesPolyAA <- paste0(paste0(groupedPolyAAIndicies), collapse="; ")


      # Get the translations for each polyQ region and store them and their lengths.
      polyAATranslations <- character(length(proteinRow$numPolyAA)) # Intialize storage vectors
      polyAALengths <- character(length(proteinRow$numPolyAA))
      polyAARegionAAFractions <- numeric(length(proteinRow$numPolyAA))

      for (i in 1:proteinRow$numPolyAA) {
        group <- groupedPolyAAIndicies[[i]] # Indexing into a list is done with double brackets

        polyAATranslations[i] <- substr(proteinRow$peptideSeq, group[1], tail(group, 1))
        polyAALengths[i] <- nchar(polyAATranslations[i])

        nAAs = str_count(polyAATranslations[i], candidateAA)
        polyAARegionAAFractions[i] <- round(nAAs/nchar(polyAATranslations[i]), 4)
      } # for (i in 1:proteinRow$numPolyAA)

      proteinRow$translationsPolyAA <- paste0(polyAATranslations, collapse = "; ")

      polyAALengths <- as.numeric(polyAALengths)
      proteinRow$lengthPolyAA <- paste0(sort(polyAALengths, decreasing = TRUE), collapse = "; ")
      proteinRow$avgLengthPolyAA <- round(mean(polyAALengths), 4)
      proteinRow$maxLengthPolyAA <- polyAALengths[1]

      proteinRow$polyAARegionAAFractions <- paste0(sort(polyAARegionAAFractions, decreasing = TRUE), collapse = "; ")
      proteinRow$avgPolyAARegionAAFractions <- round(mean(polyAARegionAAFractions), 4)
      proteinRow$maxPolyAARegionAAFractions <- polyAARegionAAFractions[1]

      # Compute p-values for protein. See NOTE 1 at bottom for details.
      # n = proteinRow$maxLengthPolyAA # get max length of polyQ tract found for each protein
      # m = proteinRow$peptideLength # Protein length
      # proteinRow$pvalues <-  1-(1-1/(20**n))**(m-(n-1)) # compute p-values

    } # if else (length(groupedPolyAAIndicies[[1]]) == 0) {

    # Add the row back into dataframe we are building up.
    proteinDf[sequenceNum,] <- proteinRow

  } # for (sequenceNum in 1:nrow(proteinDf))

  # Return a column that gives length of largest region
  # and the number of regions present (and the Qfraction of the largest region?)
  print("Finished TestHMM!")
  return (proteinDf)
} # TestHMM = function(proteinDf, model, candidateAA){


#########################################
# Main
#########################################
# Train the HMMs on hand selected polyQ
# proteins and then use them to annotate rest of protein
#########################################


# Mappings from species and amino acid repeat type to list of IDs used to train the HMM.

for (model in kModels) {
  for (species in kSpecies) {
    # Don't use the filtered list. Use all. And then filter again later if need be.
    proteinDf <- read.csv(paste0(kOutputBaseDir, species, "_prots.csv"), stringsAsFactors = FALSE)
    # proteinDf <- proteinDf[1:500,] # DEBUGGING

    outputSpeciesDir <- paste0(kOutputBaseDir, model, "/", species, "/")
    dir.create(outputSpeciesDir, recursive = TRUE)

    for (candidateAA in kCandidateAAs) {
      # debugging; candidateAA = kCandidateAAs[1]; candidateAA

      print(paste0("Species: ", species,"; Amino acid: ", candidateAA, "; Model: ", model))
      outputAADir <- paste0(outputSpeciesDir, candidateAA, "/")
      dir.create(outputAADir, recursive = TRUE)

      # Train model
      proteins <- read.csv("C:/UROPs/polyQ_neuronal_proteins/output/fly_prots.csv", stringsAsFactors = FALSE)
      print(kTrainingSetMap)
      trainingSetIDs <- kTrainingSetMap[["flyQ"]]
      print(paste0("TrainingSetIDs: ", trainingSetIDs))
      trainingSet <- proteins %>% filter(ensembl_peptide_id %in% trainingSetIDs)


      print(paste0("Look here!: ", trainingSet))
      trainedModel <- TrainHMM(trainingSet, candidateAA)

      # Annotate the proteins based on the model (either the trained version or an adjusted version of the trained version)
      if (model == "trained") {
        proteinDfHMM <- TestHMM(proteinDf, trainedModel$model, candidateAA)

        # Plot trained model
        pdf(paste0(outputAADir, species, "_", candidateAA,"_trainedHMM_graph.pdf"))
        p <- plot(trainedModel$model,
                  vertex.label = "names",
                  combine.slices = 0,
                  main = "Trained Model")
        dev.off()

      } else if (model == "adjusted") {
        # Use adjusted model
        emissProbs <- matrix(c(0.95, 0.05,
                          0.05, 0.95),
                        nrow=2, ncol=2, byrow = TRUE)
        transProbs <- matrix(c(0.995, 0.005,
                          0.0676, 0.9324),
                        nrow = 2, ncol = 2, byrow = TRUE)
        initialProbs <- c(.99, .01)
        trainingSeqs <- FormatForHmm(trainingSet, candidateAA)

        hardCodedModel <- build_hmm(observations = trainingSeqs,
                                      initial_probs = initialProbs,
                                      transition_probs = transProbs,
                                      emission_probs = emissProbs,
                                      state_names = c("NonPolyAA", "PolyAA"))

        proteinDfHMM <- TestHMM(proteinDf, hardCodedModel, candidateAA)

        # Plot adjusted model
        pdf(paste0(outputAADir, species, "_", candidateAA,"_adjustedHMM_graph.pdf"))
        p <- plot(hardCodedModel,
                  vertex.label = "names",
                  combine.slices = 0,
                  main = "Adjusted Model")
        dev.off()
      }

      # Order so most relevant proteins are at the top.
      proteinDfHMM <- proteinDfHMM %>% arrange(desc(hmmHasPolyAA), desc(maxLengthPolyAA), desc(avgPolyAARegionAAFractions))

      # Save HMM annotation results.
      write.csv(proteinDfHMM, file = paste0(outputAADir, species, "_prots_w_HMM_", candidateAA, ".csv"), row.names = FALSE)

      # Find genes where alternative splicing results in some isoforms being annoated
      # for at least one polyQ region and other isoforms not having any polyQs.
      proteinDfHMMaltSpliced <- proteinDfHMM %>%
                                      group_by(ensembl_gene_id) %>%
                                      filter(sum(hmmHasPolyAA) < n(),  1 < sum(hmmHasPolyAA)) %>%
                                      ungroup()
      # Reorder
      proteinDfHMMaltSpliced <- proteinDfHMMaltSpliced %>%
                                      arrange(ensembl_gene_id,
                                              desc(hmmHasPolyAA),
                                              desc(maxLengthPolyAA),
                                              desc(avgPolyAARegionAAFractions))

      write.csv(proteinDfHMMaltSpliced, file = paste0(outputAADir, species, "_prots_w_HMM_", candidateAA, "_alt_spliced.csv"), row.names = FALSE)

      #################################################################################33
      # Save versions without nuclear proteins.

      proteinDfHMMnucFil <- AnnotateByCategorySet(proteinDfHMM, kNuclearCats)
      proteinDfHMMnucFil <- proteinDfHMMnucFil %>% filter(!proteinDfHMMnucFil$inSet)
      write.csv(proteinDfHMMnucFil[,-ncol(proteinDfHMMnucFil)], file = paste0(outputAADir, species, "_prots_w_HMM_", candidateAA, "_nuclear_filt.csv"), row.names = FALSE)

      #######
      # Repeat filtering for alternatively spliced list.
      proteinDfHMMnucFilAltSpl <- AnnotateByCategorySet(proteinDfHMMaltSpliced, kNuclearCats)
      proteinDfHMMnucFilAltSpl <- proteinDfHMMnucFilAltSpl %>% filter(!proteinDfHMMnucFilAltSpl$inSet)
      # Save
      write.csv(proteinDfHMMaltSpliced, file = paste0(outputAADir, species, "_prots_w_HMM_", candidateAA, "_alt_spliced_nuclear_filt.csv"), row.names = FALSE)

    } # for (candidateAA in kCandidateAAs)
  } # for (species in kSpecies) {
} # for (model in kModels) {



#########################################
# Perform neuronal transcriptome filter
#########################################
for (candidateAA in kCandidateAAs) {
  for (species in kSpecies) {
    for (model in kModels) {
      outputDir <- paste0(kOutputBaseDir, model, "/", species, "/", candidateAA, "/")

      files <- list.files(outputDir, "*.csv")
      files <- files[!unlist(lapply(files, function(fileName) grepl('*transcriptome', fileName)))]
      for (file in files) {
        print(file)
        proteinDfHMM <- read.csv(paste0(outputDir, file))
        data_filtered <- proteinDfHMM %>% filter(proteinDfHMM$ensembl_peptide_id %in% kNeuronalTranscripts$V1)
        write.csv(data_filtered, paste0(outputDir, substr(file, 1, nchar(file)-4), "_transcriptome_filt.csv"), row.names = FALSE)
      }
    }
  }
}
#END
