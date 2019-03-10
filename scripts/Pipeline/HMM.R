# Daniel Alfonsetti
# daniel.alfonsetti@gmail.com; alfonset@mit.edu
# MIT Littleton Lab
# 7 March, 2019
# Description: 

#########################################
# Load libraries and set 'global' variables 
#########################################
rm(list = ls())

# Efficiency tools
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

# Enrichment analysis tools
library(clusterProfiler) 

#########################################
# Setup environment/set 'global' variables
#########################################

registerDoParallel(cores = 4)

setwd("C:/UROPs/polyQ_neuronal_proteins")
plot_output_dir <- "C:/UROPs/polyQ_neuronal_proteins/plots/"
data_dir <- "C:/UROPs/polyQ_neuronal_proteins/data/"

output_base_dir <- "C:/UROPs/polyQ_neuronal_proteins/output/"

species_vec = c("fly")

candidate_AA_vec = c("D", "T", "S", "E", "P", "G", "A", "C", "V", "M",
                     "I", "L", "Y", "F", "H", "K", "R", "W", "Q", "N")


# Categories to filter out of analysis. Anything that has to do with nucleus or DNA regulation.

# Transcription factor categories
term1 <- "GO:0003676" # nucleic acid binding
transcription_cats <- c(term1, get(term1, GOMFOFFSPRING))

term1 <- "GO:0005634" # Nucleus
term2 <- "GO:0009295" # nucleoid 
nuclear_cats <- c(term1, get(term1, GOCCOFFSPRING), 
                       term2, get(term2, GOCCOFFSPRING),
                       transcription_cats)

#########################################
# Helper Functions
#########################################

format_for_HMM = function(proteinSet, candidate_AA) {
  ########################
  # Get data into right form for HMM
  ########################
  print("At start of format_for_HMM!")
  
  list_seqs_split <- list()
  for (i in 1:nrow(proteinSet)) {
    row <- proteinSet[i,]
    split_seq_aas <- unlist(strsplit(row$peptide_seq, split = ""))
    
    split_seq <- rep(NA, length(split_seq_aas))
    
    for (j in 1:length(split_seq)) {
      if (split_seq_aas[j] == candidate_AA) {
        split_seq[j] <- "TargetAA"
      } else {
        split_seq[j] <- "NonTarget"
      }
    }
    list_seqs_split[[length(list_seqs_split)+1]] <- split_seq
  }
  print("In format_for_HMM: Split the seqs")
  
  # Pad the list with NAs
  maxlen <- max(lengths(list_seqs_split))
  list_seqs_split_padded <- lapply(list_seqs_split, `length<-`, maxlen)
  print("In format_for_HMM: Padded the seqs")
  
  # Conver to dataframe
  split_seq_result_df <- as.data.frame(do.call(rbind, list_seqs_split_padded), stringsAsFactors= FALSE)
  print("In format_for_HMM: Converted list of seqs into a dataframe!")
  
  # Convert to a seqdef object, part of the seqHMM API.
  seqs_object <- seqdef(split_seq_result_df)
  print("In format_for_HMM: Converted seq dataframe to seqs_object (seqHMM package object)!")
  
  print("Finished format_for_HMM!")
  return(seqs_object)
}

train_HMM = function(training_set, candidate_AA) {
  print("In train_HMM!")
  
  training_seqs <- format_for_HMM(training_set, candidate_AA)
  
  # Make the model
  # Initial values for emission matrices
  emiss <- matrix(c(0.95, 0.05,
                    0.26, 0.74),
                  nrow=2, ncol=2, byrow = TRUE)

  # Initial values for transition matrix
  trans <- matrix(c(0.95, 0.05, 
                    0.01, 0.99), 
                  nrow = 2, ncol = 2, byrow = TRUE)

  
  # Initial values for initial state probabilities
  initial_probs <- c(.99, .01)
  
  
  hmm <- build_hmm(observations = training_seqs,
                   initial_probs = initial_probs,
                   transition_probs = trans,
                   emission_probs = emiss,
                   state_names = c("NonPolyAA", "PolyAA"))
  
  # Fit/tune the model
  hmm_fit <- fit_model(hmm, 
                       em_step = TRUE,
                       global_step = TRUE, 
                       local_step = FALSE, 
                       control_em = list(maxeval = 5, print_level = 2),
                       control_global = list(maxtime = 1),
                       threads = 10)
  
  print("Finished train_HMM!")
  return(hmm_fit)
}

test_HMM = function(proteinSet, model, candidate_AA){
  print("At start of test_HMM!")
  
  # Remoeve any proteins with a sequence of length 1
  proteinSet <- proteinSet[nchar(as.character(proteinSet$peptide_seq))>1,]
  
  # Create model using train model parameters. The 'x' variable while
  # be updated with a single sequence each time we iterate over a new row of
  # proteinSet.
  x <- NULL
  HMM <- dthmm(x, 
               Pi = model$transition_probs,  # Looks right 
               delta = model$initial_probs,  # Looks right
               distn = "binom",             
               pm = list(prob = model$emission_probs[1:2, 2], size = 1),       # 
               pn = NULL,
               discrete = TRUE,
               nonstat = TRUE)
  print("In test_HMM: Finished building the object!")

  # Columns for fractions of Single, double, and triple AAs of the particularly type in the entire protein sequence.
  proteinSet$AAComp = rep(NA, nrow(proteinSet)) # Errors here
  # proteinSet$DipeptideComp = rep(NA, nrow(proteinSet))
  # proteinSet$TripeptideComp = rep(NA, nrow(proteinSet))
  
  # A boolean indicating whether or not at least one polyQ region was found by the HMM for this protein.
  proteinSet$HMMhasPolyAA = rep(NA, nrow(proteinSet))
  
  # A column where each element is a number indicating the number of polyQ regions found on this protein by the HMM.
  proteinSet$NumberPolyAA = rep(NA, nrow(proteinSet))
  
  # A column where each element contains a list with the indicies for each region in the protein that is denoted polyQ by the HMM
  proteinSet$IndiciesPolyAA = rep(NA, nrow(proteinSet)) 
  
  # A column where each element contains a list with the AA sequence for each region in the protein that is denoted polyQ by the HMM
  proteinSet$TranslationsPolyAA = rep(NA, nrow(proteinSet)) 
  
  # A column where each contains a list of the length of the region that is denoted polyQ by the HMM
  proteinSet$LengthsPolyAA = rep(NA, nrow(proteinSet))
  proteinSet$AvgLengthsPolyAA = rep(NA, nrow(proteinSet))
  proteinSet$MaxLengthsPolyAA = rep(NA, nrow(proteinSet))
  
  # A column where each contains a list of the length of the region that is denoted polyQ by the HMM
  proteinSet$PolyAARegionAAFractions = rep(NA, nrow(proteinSet))
  proteinSet$AvgPolyAARegionAAFractions = rep(NA, nrow(proteinSet))
  proteinSet$MaxPolyAARegionAAFractions = rep(NA, nrow(proteinSet))
  
  proteinSet$annotatedSequence = rep(NA, nrow(proteinSet))
  proteinSet$PeptideLength = rep(NA, nrow(proteinSet))
  
  
  proteinSet$pvalues = rep(0, nrow(proteinSet))
  
  # TODO: Add amino acid descriptor, dipeptide descriptor, and tri peptide descriptor.
  # Add autocorrelatoin descriptor.
  
  
  print("In test_HMM: Extracting information about annotated paths!")
  for (sequenceNum in 1:nrow(proteinSet)) {
    
    if (sequenceNum %% 5000 == 0){
      print(paste0("Annotating sequence ", sequenceNum, " out of ", nrow(proteinSet)))
    }
    proteinRow <- proteinSet[sequenceNum,]
    
    ### TURN TRANSLATION INTO BINARY VECTOR (coded as 1 if AA == Q, 0 otherwise)
    binarySequence = rep(1, nchar(proteinRow$peptide_seq))
    for (i in 1:length(binarySequence)){
      character = substr(proteinRow$peptide_seq, i, i)
      if (character == candidate_AA){
        binarySequence[i] = 1
      } else { binarySequence[i] = 0}
    }
    
    ### ADD BINARY SEQUENCE VECTOR TO HMM INPUTS
    HMM$x = binarySequence
    
    ### CREATING AN ANNOTATED (STATES) SEQUENCE FOR THE TRANSLATION
    annotatedSequence = Viterbi(HMM)
    proteinRow$annotatedSequence <- paste0(annotatedSequence, collapse = "")
    
    # Get Qfraction for the entire protein.
    

    # proteinRow$AAComp <- extractAAC(proteinRow$peptide_seq)[paste0(replicate(1, candidate_AA), collapse = "")]
    # proteinRow$DipeptideComp <- extractDC(proteinRow$peptide_seq)[paste0(replicate(2, candidate_AA), collapse = "")]
    # proteinRow$TripeptideComp <- extractTC(proteinRow$peptide_seq)[paste0(replicate(3, candidate_AA), collapse = "")]
    nAAs = str_count(proteinRow$peptide_seq, candidate_AA)
    proteinRow$AAComp= round(nAAs/nchar(proteinRow$peptide_seq), 4)
    
    # Get length of protein.
    proteinRow$PeptideLength <- nchar(proteinRow$peptide_seq)
    
    # Find which indices of the translation are marked as being polyQ
    polyAA_indices <- which(annotatedSequence %in% c(2))
    
    # Group the indicies into groups (so if you have multiple polyQ regions, the indices for each polyQ region will be seperated from each other)
    grouped_polyAA_indicies <- split(polyAA_indices, cumsum(c(1, diff(polyAA_indices) != 1))) 
    
    # Only keep polyAAs that are greater than 6AAs in length here
    grouped_polyAA_indicies <- grouped_polyAA_indicies[lapply(grouped_polyAA_indicies, length) > 6]
    
    # Record how many polyQ regions there are
    if (length(grouped_polyAA_indicies) == 0) {
      proteinRow$NumberPolyAA = 0
      proteinRow$HMMhasPolyAA = FALSE
    } else {
      proteinRow$NumberPolyAA = length(grouped_polyAA_indicies)
      proteinRow$HMMhasPolyAA = TRUE
      proteinRow$IndiciesPolyAA <- paste0(paste0(grouped_polyAA_indicies), collapse="; ")
      
      
      # Get the translations for each polyQ region and store them and their lengths. 
      PolyAATranslations <- character(length(proteinRow$NumberPolyAA)) # Intialize storage vectors
      PolyAALengths <- character(length(proteinRow$NumberPolyAA))
      PolyAARegionAAFractions <- numeric(length(proteinRow$NumberPolyAA))
      
      for (i in 1:proteinRow$NumberPolyAA) {
        group <- grouped_polyAA_indicies[[i]] # Indexing into a list is done with double brackets
        
        PolyAATranslations[i] <- substr(proteinRow$peptide_seq, group[1], tail(group, 1))
        PolyAALengths[i] <- nchar(PolyAATranslations[i])
        
        
        nAAs = str_count(PolyAATranslations[i], candidate_AA)
        PolyAARegionAAFractions[i] <- round(nAAs/nchar(PolyAATranslations[i]), 4)
      } # for (i in 1:proteinRow$NumberPolyAA)
      
      proteinRow$TranslationsPolyAA <- paste0(PolyAATranslations, collapse = "; ")
      
      PolyAALengths <- as.numeric(PolyAALengths)
      proteinRow$LengthsPolyAA <- paste0(sort(PolyAALengths, decreasing = TRUE), collapse = "; ")
      proteinRow$AvgLengthsPolyAA <- round(mean(PolyAALengths), 4)
      proteinRow$MaxLengthsPolyAA <- PolyAALengths[1]
      
      proteinRow$PolyAARegionAAFractions <- paste0(sort(PolyAARegionAAFractions, decreasing = TRUE), collapse = "; ")
      proteinRow$AvgPolyAARegionAAFractions <- round(mean(PolyAARegionAAFractions), 4)
      proteinRow$MaxPolyAARegionAAFractions <- PolyAARegionAAFractions[1]
      
      # Compute p-values for protein. See NOTE 1 at bottom for details.
      n = proteinRow$MaxLengthsPolyAA # get max length of polyQ tract found for each protein
      m = proteinRow$PeptideLength # Protein length
      proteinRow$pvalues <-  1-(1-1/(20**n))**(m-(n-1)) # compute p-values
      
    } # if else (length(grouped_polyAA_indicies[[1]]) == 0) {
    
    # Add the row back into dataframe we are building up.
    proteinSet[sequenceNum,] <- proteinRow
    
  } # for (sequenceNum in 1:nrow(proteinSet))
  
  # Return a column that gives length of largest region 
  # and the number of regions present (and the Qfraction of the largest region?)
  print("Finished test_HMM!")
  return (proteinSet)
}

#########################################
# Main - Train the HMMs on hand selected polyQ proteins and then use them to annotate rest of protein 
#########################################

# Mappings from species and amino acid repeat type to list of IDs used to train the HMM.
training_id_dict <- list("mouse_Q"= c("FBpp0289769","FBpp0307700", "FBpp0086727", "FBpp0111724", "FBpp0070830", "FBpp0309352", "FBpp0402897", "FBpp0110299","FBpp0305807"),
                         "mouse_M" = c("NA", "NA"),
                         "worm_Q" = c("Something", "Something"))

models_dict <- list("adjusted" = NA, "trained" = NA)

for (model in names(models_dict)) {
  for (species in species_vec) {
    # Don't use the filtered list. Use all. And then filter again later if need be.
    protein_df_w_seqs <- read.csv(paste0(output_base_dir, species, "_prots.csv"), stringsAsFactors = FALSE)
    
    output_species_dir <- paste0(output_base_dir, "/", model, "/", species, "/")
    dir.create(output_species_dir, recursive = TRUE)
    
    for (candidate_AA in candidate_AA_vec) {
      # debugging; candidate_AA = candidate_AA_vec[1]; candidate_AA
      
      print(paste0("Species: ", species,"; Amino acid: ", candidate_AA, "; Model: ", model))
      output_AA_dir <- paste0(output_species_dir, candidate_AA, "/")
      dir.create(output_AA_dir, recursive = TRUE)
      
      # Train model
      x <- read.csv("C:/UROPs/polyQ_neuronal_proteins/output/fly_prots.csv", stringsAsFactors = FALSE)
      training_set_ids <- training_id_dict[["mouse_Q"]]
      training_set <- x %>% filter(ensembl_peptide_id %in% training_set_ids)
      training_seqs <- format_for_HMM(training_set, candidate_AA) # We aren't actually using these. But 
      trained_model <- train_HMM(training_set, candidate_AA)

      

      # Use trained HMM model on testSet (proteinSet).
      proteinSet <- protein_df_w_seqs
      
      # Annotate the proteins based on the model (either the trained version or an adjusted version of the trained version)
      if (model == "trained") {
        protein_df_w_hmm <- test_HMM(proteinSet, trained_model$model, candidate_AA)
        
        # Plot trained model
        pdf(paste0(output_dir, species, "_", candidate_AA,"_trainedHMM_graph.pdf"))
        p <- plot(trained_model$model,
                  vertex.label = "names",
                  combine.slices = 0,
                  main = "Trained Model") # Intial hmm
        # ssplot(trained_model, plots = "both", border = NA)
        dev.off()
        
      } else if (model == "adjusted") {
        # Use adjusted model
        random_set <- sample_n(x, 1000)
        emiss <- matrix(c(0.95, 0.05,
                          0.05, 0.95),
                        nrow=2, ncol=2, byrow = TRUE)
        trans <- matrix(c(0.995, 0.005, 
                          0.0676, 0.9324), 
                        nrow = 2, ncol = 2, byrow = TRUE)
        initial_probs <- c(.99, .01)
        
        hard_coded_model <- build_hmm(observations = training_seqs,
                                      initial_probs = initial_probs,
                                      transition_probs = trans,
                                      emission_probs = emiss,
                                      state_names = c("NonPolyAA", "PolyAA"))
        
        protein_df_w_hmm <- test_HMM(proteinSet, hard_coded_model, candidate_AA)
        
        
        # Plot adjusted model
        pdf(paste0(output_dir, species, "_", candidate_AA,"_adjustedHMM_graph.pdf"))
        p <- plot(hard_coded_model,
                  vertex.label = "names",
                  combine.slices = 0,
                  main = "Adjusted Model") # Intial hmm
        # ssplot(trained_model, plots = "both", border = NA)
        dev.off()
      }
      
      # Order so most relevant proteins are at the top.
      protein_df_w_hmm <- protein_df_w_hmm %>% arrange(desc(HMMhasPolyAA), desc(MaxLengthsPolyAA), desc(AvgPolyAARegionAAFractions))
      
      # Save HMM annotation results.
      write.csv(protein_df_w_hmm, file = paste0(output_AA_dir, species, "_prots_w_HMM_", candidate_AA, ".csv"), row.names = FALSE)
      
      # Find genes where alternative splicing results in some isoforms being annoated 
      # for at least one polyQ region and other isoforms not having any polyQs.
      protein_df_w_hmm_alt_spliced <- protein_df_w_hmm %>% 
                                      group_by(ensembl_gene_id) %>% 
                                      filter(sum(HMMhasPolyAA) < n(),  1 < sum(HMMhasPolyAA))
      # Reorder
      protein_df_w_hmm_alt_spliced <- protein_df_w_hmm_alt_spliced %>%
                                      arrange(ensembl_gene_id, 
                                              desc(HMMhasPolyAA), 
                                              desc(MaxLengthsPolyAA),
                                              desc(AvgPolyAARegionAAFractions))
      write.csv(protein_df_w_hmm_alt_spliced, file = paste0(output_AA_dir, species, "_prots_w_HMM_", candidate_AA, "_alt_spliced.csv"), row.names = FALSE)
    
      
      #################################################################################33
      #################################################################################
      # Save versions without nuclear proteins.
      vec <- unlist(
        lapply(
          lapply(as.character(protein_df_w_hmm$go_ids_CC), function(x){unlist(strsplit(x, "; "))}), 
          function(x){any(x %in% nuclear_cats)}
        )
      ) # Get list of rows that have proteins that are nuclear.
      
      protein_df_w_hmm <- protein_df_w_hmm[!vec,]
      # Save
      write.csv(protein_df_w_hmm, file = paste0(output_AA_dir, species, "_prots_w_HMM_", candidate_AA, "_nuclear_filt.csv"), row.names = FALSE)
      
      #######
      # Repeat filtering for alternatively spliced list.
      vec <- unlist(
        lapply(
          lapply(as.character(protein_df_w_hmm_alt_spliced$go_ids_CC), function(x){unlist(strsplit(x, "; "))}), 
          function(x){any(x %in% nuclear_cats)}
        )
      )
      protein_df_w_hmm_alt_spliced <- protein_df_w_hmm_alt_spliced[vec,]
      # Save
      write.csv(protein_df_w_hmm_alt_spliced, file = paste0(output_AA_dir, species, "_prots_w_HMM_", candidate_AA, "_alt_spliced_nuclear_filt.csv"), row.names = FALSE)
      
    } # for (candidate_AA in candidate_AA_vec)
  } # for (species in species_vec) {
} # for (model in models) {


# Perform neuronal transcriptome filter
neuronal_transcripts <- as.vector(read.table("C:/UROPs/polyQ_neuronal_proteins/output/fly_CNS_transcriptome_mh-l.txt", sep = "\t"))
for (candidate_AA in candidate_AA_vec) {
  for (species in species_vec) {
    for (model in models) {
      output_dir <- paste0(output_base_dir, model, "/", species, "/", candidate_AA, "/")
      for (file in list.files(output_dir, "*.csv")) {
        data <- read.csv(paste0(output_dir, file))
        data_filtered <- data %>% filter(data$ensembl_peptide_id %in% neuronal_transcripts$V1)
        write.csv(data_filtered, paste0(output_dir, substr(file, 1, nchar(file)-4), "_transcriptome_filt.csv"))
      }
    }
  }
}

# Summary Plots
for (model in c("adjusted")) {
  for (species in species_vec) {
    for (candidate_AA in candidate_AA_vec) {
      output_dir <- paste0(output_base_dir, model, "/", species, "/", candidate_AA, "/")
      dir.create(output_dir, recursive = TRUE)
      
      # Plotting helper function
      # Number of randomly chosen 50 amino acid peptide chunks
      
      rnd_substr <- function(x, length) {
        if (nchar(x) < length) {return(x)}
        else {
          start = sample(1:(nchar(x)-length), 1, replace=T)
          end = start + (length-1)
          random_substring = substr(x, start, end)
          return(random_substring)
        }
      } # rnd_substr <- function(x, length)
      #########################################
      
      pdf(paste0(output_dir, "/sequence_summary_plots.pdf"))
      #########################################
      # Plot type: Whole Protein Qfraction Histogram for all proteins
      p <- ggplot(data = proteins, aes(x = proteins$AAComp)) +
        geom_histogram(show.legend = FALSE) +
        xlab(paste0("%", candidate_AA, " (Whole Protein Seqs)")) +
        ylab("# Proteins") +
        ggtitle(paste0("%", candidate_AA, " Whole Protein")) +
        theme_light()
      print(p)
      
      filtered_df <- proteins %>% filter(proteins$AAC> quantile(proteins$AAfraction, 0.95))
      filtered_df <- filtered_df %>% arrange(desc(filtered_df$AAfraction))
      write.csv(filtered_df, file = paste0(output_dir, species, "_", candidate_AA, "fraction_prots_95sig.csv"))
      
      
      # Plot type: Whole Protein AAC Histogram split by HMMhasPolyAA
      p <- ggplot(data = proteins, aes(x = AAComp, fill = HMMhasPolyAA)) +
        geom_histogram(alpha=0.5, position="identity") +
        xlab(paste0("%", candidate_AA, " (Whole Protein Seqs)")) +
        ylab("# Proteins") +
        labs(fill = paste0("Contains \n poly", candidate_AA)) +
        ggtitle(paste0("%", candidate_AA," Whole Proteins Histogram")) +
        theme_light()
      print(p)
      
      p <- ggplot(data = proteins, aes(x = AAComp, fill = HMMhasPolyAA)) +
        geom_density(alpha = 0.5) +
        xlab(paste0("%", candidate_AA, " (Whole Protein Seqs)")) +
        ylab("Density") +
        labs(fill = paste0("Contains \n poly", candidate_AA)) +
        ggtitle(paste0("%", candidate_AA, " Whole Proteins Density Plot")) +
        theme_light()
      print(p)
      
      ########################################################
      # Plot type: 100000 50AA chunks, each randomly chosen from randomly chosen Protein seqs
      df <- sample_n(proteins, 100000, replace = TRUE)
      df$RandomSeqs <-unlist(lapply(as.character(df$peptide_seq), rnd_substr, 50)) # TODO: Have length be randomly choosen from lengths of polyAA seuqneces.
      nAAs = str_count(df$RandomSeqs, candidate_AA)
      df$AAComp <- nAAs/nchar(df$RandomSeqs)
      
      
      p <- ggplot(data = df, aes(x = df$AAComp)) +
        geom_histogram(show.legend = FALSE) +
        xlab(paste0("%", candidate_AA, " (Random 50AA Seqs)")) +
        ylab("# of Randomly Choosen Pepitides") +
        ggtitle(paste0("%", candidate_AA, " for 100000 50AA Random Seqs Histogram")) +
        theme_light()
      print(p)
      
      
      
      # df <- sample_n(proteins, sum(!proteins$HMMhasPolyAA), replace = TRUE)
      # df$RandomSeqs <-unlist(lapply(as.character(df$peptide_seq), rnd_substr, 50)) # TODO: Have length be randomly choosen from lengths of polyAA seuqneces.
      # nAAs = str_count(df$RandomSeqs, candidate_AA)
      # df$AAComp <- nAAs/nchar(df$RandomSeqs)
      # 
      # p <- ggplot() +
      #   geom_histogram(data = proteins, aes(x = AAComp, fill = HMMhasPolyAA), alpha = 0.5) +
      #   geom_histogram(data = df, aes(x = df$AAComp), alpha = 0.5) +
      #   xlim(c(0, 0.5)) +
      #   ggtitle("Something")
      # print(p)
      
      filtered_df <- df %>% filter(df$AAComp > quantile(df$AAComp, 0.95))
      filtered_df <- filtered_df %>% arrange(desc(filtered_df$AAComp))
      write.csv(filtered_df, file = paste0(output_dir, species, "_", candidate_AA, "fraction_randomSeqs_95sig.csv"))
      
      # Plot type: Randomly chosen AA chunks of the same size and number as polyQ protein regions
      
      ########################################################
      # Plot type: Number of Proteins annotated by HMM to be "Poly AA" vs %AA 
      # (both avg %AA and max %AA for the polyAA regions in a given polyAA protein)
      
      p <- ggplot(data = df, aes(x = as.numeric(df$AvgPolyAARegionAAFractions))) +
        geom_histogram(show.legend = FALSE) +
        xlab(paste0("Avg %", candidate_AA, " of poly", candidate_AA, " Regions")) + 
        ggtitle(paste0("Avg %", candidate_AA, " of poly", candidate_AA," Regions for poly", candidate_AA, " Containing Proteins")) +
        ylab(paste0("# Proteins Containing poly", candidate_AA)) +
        theme_light()
      print(p)
      
      df <- proteins %>% filter(HMMhasPolyAA)
      p <- ggplot(data = df, aes(x = df$MaxPolyAARegionAAFractions, fill = "")) +
        geom_histogram(show.legend = FALSE) +
        xlab(paste0("Max %", candidate_AA, " of poly", candidate_AA, " Regions")) +
        ggtitle(paste0("Max %", candidate_AA, " of Poly", candidate_AA, " Regions for poly", candidate_AA, " Containing Proteins")) +
        ylab(paste0("# Proteins Containing poly", candidate_AA)) +
        theme_light()
      print(p)
      
      # Now compare to randomly chosen sequences of the same size
      polyAA_df <- proteins %>% filter(HMMhasPolyAA)
      random_seq_df <- data.frame(matrix(nrow = nrow(polyAA_df), ncol = 1))
      colnames(random_seq_df) <- c("AAComp")
      for (i in 1:nrow(polyAA_df)) {
        target_length <- polyAA_df[i,]$MaxLengthsPolyAA
        flag = TRUE
        while (flag){
          try = as.character(sample_n(proteins, 1)$peptide_seq)
          if (nchar(try) > target_length) {
            res <- rnd_substr(try, target_length)
            nAAs = str_count(res, candidate_AA)
            random_seq_df$AAComp[i] <- nAAs/nchar(res)
            flag = FALSE
          }
        }
      }
      
      p <- ggplot() +
        geom_histogram(data = polyAA_df, aes(x = polyAA_df$MaxPolyAARegionAAFractions, fill = "Red")) +
        geom_histogram(data = random_seq_df, aes(x= random_seq_df$AAComp, fill = "Blue")) + 
        xlab(paste0("%", candidate_AA)) +
        ggtitle(paste0("Max %", candidate_AA, " of Poly", candidate_AA, " Regions for poly", candidate_AA, " Containing Proteins \n
                       vs %Q of randomly chosen tracts")) +
        ylab(paste0("# Proteins", candidate_AA)) +
        theme_light()
      print(p)
      
      p <- ggplot() +
        geom_histogram(data = polyAA_df, aes(x = polyAA_df$AvgPolyAARegionAAFractions, fill = "Red")) +
        geom_histogram(data = random_seq_df, aes(x= random_seq_df$AAComp, fill = "Blue")) + 
        xlab(paste0("%", candidate_AA)) + 
        ggtitle(paste0("Avg %", candidate_AA, " of poly", candidate_AA," Regions for poly", candidate_AA, " Containing Proteins \n vs %Q of randomly chosen tracts")) +
        ylab(paste0("# Proteins", candidate_AA)) +
        theme_light()
      print(p)
      
      # All three on the same chart
      # Max %Q vs Avg %Q vs Random Tracts
      p <- ggplot() +
        geom_histogram(data = polyAA_df, aes(x = polyAA_df$AvgPolyAARegionAAFraction, fill= "Avg"), alpha = 0.5) +
        geom_histogram(data = polyAA_df, aes(x = polyAA_df$MaxPolyAARegionAAFractions, fill= "Max"), alpha= 0.5) +
        geom_histogram(data = random_seq_df, aes(x= random_seq_df$AAComp, fill="Random"), alpha = 0.5) + 
        xlab(paste0("%", candidate_AA)) + 
        ggtitle(paste0("Avg %", candidate_AA, " vs Max %", candidate_AA, " of poly", candidate_AA," Regions \n vs %Q of randomly chosen tracts")) +
        ylab(paste0("# Proteins", candidate_AA)) +
        scale_fill_manual(name = "", values = c("Avg" = "blue", "Max" = "green", "Random" = "red")) 
      print(p)
      
      ##############
      # Plot type: Number of Proteins annotated by HMM to be Poly Q vs Length of Region
      
      p <- ggplot(data = df, aes(x = df$AvgLengthsPolyAA)) +
        geom_histogram(show.legend = FALSE) +
        xlab(paste0("Average Poly", candidate_AA, " Length")) +
        ylab(paste0("# Proteins Containing poly", candidate_AA, " (by HMM)")) +
        ggtitle(paste0("Average Poly", candidate_AA, " Length in Poly", candidate_AA, " Containing Proteins")) +
        theme_light()
      print(p)
      
      p <- ggplot(data = df, aes(x = df$MaxLengthsPolyAA)) +
        geom_histogram(show.legend = FALSE) +
        xlab(paste0("Max poly", candidate_AA, " Length")) +
        ylab(paste0("# Proteins Containing poly", candidate_AA, " (by HMM)")) +
        ggtitle(paste0("Max poly", candidate_AA, " in poly", candidate_AA, " Containing Proteins")) +
        theme_light()
      print(p)
      
      p <- ggplot(data = df, aes(x = df$NumberPolyAA)) + 
        geom_histogram(show.legend = FALSE) +
        xlab(paste0("Number of poly", candidate_AA, " Regions")) +
        ylab(paste0("# of ", candidate_AA, " Proteins")) +
        ggtitle(paste0("Number of poly", candidate_AA, " regions in poly", candidate_AA, " Containing Proteins")) +
        theme_light()
      print(p)
      dev.off()
    } # for (model in c("adjusted", "trained")) {
  } # for (species in species_vec)
} # for (candidate_AA in candidate_AA_vec) {


# Sanity check
# data <- read.csv("C:/UROPs/polyQ_neuronal_proteins/output/adjusted/fly/Q/fly_prots_w_HMM_Q_nuclear_filt_transcriptome_filter.csv")
# training_set_protein_ids %in% data$ensembl_peptide_id
# training_set_protein_ids %in% neuronal_transcripts$V1



# NOTE 1: Find a way to measure sigificance of polyQ proteins outputted by our HMM while controlling for length.
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
