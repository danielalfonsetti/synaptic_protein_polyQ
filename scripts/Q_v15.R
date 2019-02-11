# Daniel Alfonsetti, daniel.alfonsetti@gmail.com
# 13 Sep. 2018
# Description: polyQ analysis in synaptic proteins

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

# Enrichment analysis tools
library(clusterProfiler) 
#########################################
# Setup environment
#########################################

registerDoParallel(cores = 4)

setwd("C:/UROPs/polyQ_neuronal_proteins")
plot_output_dir <- "C:/UROPs/polyQ_neuronal_proteins/plots/"
data_dir <- "C:/UROPs/polyQ_neuronal_proteins/data/"

output_base_dir <- "C:/UROPs/polyQ_neuronal_proteins/output/"

#########################################
# User defined variables
#########################################

species_vec = c("fly")

# candidate_AA_vec = c("Q", "D")

candidate_AA_vec = c("D", "T", "S", "E", "P", "G", "A", "C", "V", "M",
                     "I", "L", "Y", "F", "H", "K", "R", "W", "Q", "N")
# species_vec = c("fly", "mouse", "human", "worm")

# test_id_dict <- list(mouse_Q = c("FBpp0289769","FBpp0307700", "FBpp0086727", "FBpp0111724", "FBpp0070830", "FBpp0309352", "FBpp0402897", "FBpp0110299","FBpp0305807"),
#                  worm_Q = c("NA", "NA"))

#########################################

# Categories to filter out of analysis. ANything that has to do with nucleus or DNA regulation.
term1 <- "GO:0005634" # Nucleus
term2 <- "GO:0009295" # nucleoid 
term3 <- "GO:0003676" # nucleic acid binding

# term3 <- "GO:0019219" # regulation of nucleobase-containing compound metabolic process

nuclear_cats <- c(term1, get(term1, GOCCOFFSPRING), 
                       term2, get(term2, GOCCOFFSPRING),
                       term3, get(term3, GOMFOFFSPRING))

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
  
  ###############
  # HMM graph outputs
  print("In test_HMM: Printing HMM graphs to file")
  pdf(paste0(output_AA_dir, species, "_", candidate_AA,"_fitHMM_graph.pdf"))
  p <- plot(hmm_fit$model,
            vertex.label = "names",
            combine.slices = 0,
            main = "Trained Model") # Intial hmm
  #ssplot(hmm, plots = "both", border = NA)
  dev.off()
  print("In test_HMM: Successfully printed to output")
  #################
  
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
# End of helper functions
#########################################

#########################################
# Download proteomes
#########################################
for (species in species_vec) {
  # debugging; species = species_vec[1]
  print(paste0("Downloading proteome for ", species, "..."))
  output_species_dir <- paste0(output_base_dir, species, "/")
  dir.create(output_species_dir, recursive = TRUE)
  # Diagnostics: Measure download time
  start_time <- Sys.time() 
  
  # datasets = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org")
  # listDatasets(datasets); Relevant data sets:
  # mmusculus_gene_ensembl
  # hsapiens_gene_ensembl
  # dmelanogaster_gene_ensembl
  # celegans_gene_ensembl
  
  # Drosophila melanogaster
  # Homo sapiens
  # Mus musculus
  # Caenorhabditis elegans
  
  if (species == "fly") {
    ensembl <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")
  } else if (species == "human")  {
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  } else if (species == "mouse") {
    ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  } else if (species == "worm") {
    ensembl <- useMart("ensembl", dataset = "celegans_gene_ensembl")
  }

  listAttributes(ensembl)
  protein_df <- getBM(attributes = c('ensembl_gene_id', 'ensembl_peptide_id', "external_gene_name", 'go_id', "name_1006", "namespace_1003", "kegg_enzyme"), 
                      values = "*", 
                      mart = ensembl)
  
  
  # protein_df <- getBM(attributes = c("kegg_enzyme"),
  #                     values = "*",
  #                     mart = ensembl)
  
  # Filter here based on whether it has CC GO annot. Make two seperate files
  nuclear_prots <- protein_df[protein_df$go_id %in% nuclear_cats, ]

  
  protein_df <- dplyr::rename(protein_df, "go_description" = name_1006 , "go_type" = namespace_1003 )
  
  
  # Paste together GO terms if they are for the same peptide and of the same GO_type.
  protein_df <- protein_df %>%  
    distinct() %>%
    group_by(go_type, ensembl_peptide_id) %>% 
    mutate(go_list = paste(go_description, collapse = "; "), go_ids = paste(go_id, collapse = "; "))
  protein_df <- protein_df %>% dplyr::select(-c(go_description, go_id)) 
  
  # Put MF, BP, and CC information on same row for a given peptide.
  protein_df <- protein_df %>%
          distinct() %>%
          group_by(ensembl_peptide_id) %>%
          mutate(go_list_MF = nth(go_list, 1), go_ids_MF = nth(go_ids, 1),
                 go_list_CC = nth(go_list, 2), go_ids_CC = nth(go_ids, 2),
                 go_list_BP = nth(go_list, 3), go_ids_BP = nth(go_ids, 3))
  protein_df <- protein_df %>% dplyr::select(-c(go_type, go_list, go_ids)) 
  
  View(protein_df)
  
  vec <- duplicated(protein_df[c("ensembl_peptide_id")]) # Get rows that share the same peptide ID.
  protein_df <- protein_df[!vec,]  # Remove duplicated rows.
  
  # Remove genes that are not protein coding.
  protein_df <- protein_df[!(protein_df$ensembl_peptide_id == ""),]

  
  # Get peptide sequences
  seqs <- getSequence(id = protein_df$ensembl_peptide_id, 
                      type="ensembl_peptide_id",
                      seqType="peptide",
                      mart=ensembl) 

  # Seqs takes the form...
  
  # peptide_seq | ensembl_peptide_id
  # ------------------------------
  # DSSQSRLRSAH | FBpp0070102
  # MGKKGKGKKGK | FBpp0070068

  # Add sequences to the dataframe we are building up
  seqs <- dplyr::rename(seqs , "peptide_seq" = peptide)
  protein_df_w_seqs <- merge(protein_df, seqs, by = "ensembl_peptide_id")
  
  #  Set empty cells equal to NA
  protein_df_w_seqs[protein_df_w_seqs == ""] <- NA
  
  # Record length of each protein
  protein_df_w_seqs$peptide_length <- unlist(lapply(protein_df_w_seqs$peptide_seq, function(x){length(unlist(strsplit(x, "")))}))
  
  # Make two versions: one with all the proteins, and one with only nuclear proteins.
  write.csv(protein_df_w_seqs, file = paste0(output_species_dir, species, "_prots.csv"), row.names = FALSE)
  
  # Something is sketchy here with the nuclear filters....
  # Only include non-nuclear proteins
  protein_df_w_seqs_filt <- protein_df_w_seqs %>% filter(!(ensembl_peptide_id %in% nuclear_prots$ensembl_peptide_id))
  write.csv(protein_df_w_seqs_filt, file = paste0(output_species_dir, species, "_prots_nuclear_filt.csv"), row.names = FALSE)
  
  #####################
  print(paste0("Downloaded and saved proteome for ", species, " in ", output_species_dir, "!"))
}
##########################
# Downloads complete
##########################


# Mappings from species and amino acid repeat type to list of IDs used to train the HMM.
training_id_dict <- list("mouse_Q"= c("FBpp0289769","FBpp0307700", "FBpp0086727", "FBpp0111724", "FBpp0070830", "FBpp0309352", "FBpp0402897", "FBpp0110299","FBpp0305807"),
                         "mouse_M" = c("NA", "NA"),
                         "worm_Q" = c("Something", "Something"))

models_dict <- list("adjusted" = NA, "trained" = NA)

for (model in names(models_dict)) {
  
  output_model_dir <- paste0(output_base_dir, model, "/")
  
  # Now do the analyses...
  for (species in species_vec) {
    # debugging; species = species_vec[1]
    
    output_species_dir <- paste0(output_model_dir, species, "/")
    dir.create(output_species_dir, recursive = TRUE)
    
    # Don't use the filtered list. Use all. And then filter again later if need be.
    protein_df_w_seqs <- read.csv(paste0(output_base_dir, species, "_prots.csv"), stringsAsFactors = FALSE)
  
    for (candidate_AA in candidate_AA_vec) {
      # debugging; candidate_AA = candidate_AA_vec[1]; candidate_AA
      
      print(paste0("Species: ", species,"; Amino acid: ", candidate_AA))
      output_AA_dir <- paste0(output_species_dir, candidate_AA, "/")
      dir.create(output_AA_dir, recursive = TRUE)
      
      x <- read.csv("C:/UROPs/polyQ_neuronal_proteins/output/fly_prots.csv", stringsAsFactors = FALSE)
      
      training_set_ids <- training_id_dict[["mouse_Q"]]
      training_set <- x %>% filter(ensembl_peptide_id %in% training_set_ids)
      training_seqs <- format_for_HMM(training_set, candidate_AA) # We aren't actually using these. But 
      trained_model <- train_HMM(training_set, candidate_AA)
      
      
      # build function needs them for some reason...
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
      
      # Evaluate/inspect models (debugging)
          # Evaluate compare hard-coded model with trained model on trained model's training data.
          # hard_coded_model
          # trained_model$model
          # training_df_w_trained_hmm <- test_HMM(training_set, trained_model$model, candidate_AA)
          # View(training_df_w_trained_hmm)
          # training_df_w_hard_hmm <- test_HMM(training_set, hard_coded_model, candidate_AA)
          # View(training_df_w_hard_hmm)
          # 
          # 
          # # Evaluate compare hard-coded model with trained model on random set of protein data
          # random_df_w_trained_hmm <- test_HMM(random_set, trained_model$model, candidate_AA)
          # random_df_w_hard_hmm <- test_HMM(random_set, hard_coded_model, candidate_AA)
          # trained <- random_df_w_trained_hmm %>%
          #   filter(HMMhasPolyAA == TRUE)
          # mean(trained$MaxPolyAARegionAAFractions)
          # hard <- random_df_w_hard_hmm %>%
          #   filter(HMMhasPolyAA == TRUE)
          # mean(hard$MaxPolyAARegionAAFractions)
          # nrow(hard)
          # nrow(trained)
          # View(random_df_w_hard_hmm)
          # View(random_df_w_trained_hmm)
      
      # Use trained HMM model on testSet (proteinSet).
      proteinSet <- protein_df_w_seqs
      
      # Annotate the proteins based on the model we have trained and generate various statistics for each annotation (or lack thereof)
      if (model == "trained") {
        protein_df_w_hmm <- test_HMM(proteinSet, trained_model$model, candidate_AA)
      } else if (model == "adjusted") {
        protein_df_w_hmm <- test_HMM(proteinSet, hard_coded_model, candidate_AA)
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


###################
# Further Analyses
###################
for (candidate_AA in candidate_AA_vec) {
  # debugging; species = species_vec[1]
  output_species_dir <- paste0(output_base_dir, species, "/")
  dir.create(output_species_dir, recursive = TRUE)
  
  protein_df_w_seqs <- read.csv(paste0(output_species_dir, species, "_prots_filt.csv"), stringsAsFactors = FALSE)
  
  for (species in species_vec) {
    # debugging; candidate_AA = candidate_AA_vec[1]; candidate_AA
    
    print(paste0("Species: ", species,"; Amino acid: ", candidate_AA))
    output_AA_dir <- paste0(output_species_dir, candidate_AA, "/")
    dir.create(output_AA_dir, recursive = TRUE)
    
    #########################################
    # ClusterProfiler module
    #########################################
    
    # Load in protein set data from local drive
    # Error here
    proteins <- read.csv( paste0(output_AA_dir, species, "_", nrow(protein_df_w_hmm), "prots_w_HMM_", candidate_AA, ".csv"))

    polyAA_proteins = proteins %>% filter(proteins$HMMhasPolyAA == TRUE)
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
      

      proteins_vec <- unique(proteins_has_ont_type$ensembl_peptide_id)
      # length(proteins_vec); length(proteins_has_ont_type$ensembl_peptide_id)
      
      # Get only significant proteins
      polyAA_proteins <- proteins_has_ont_type %>% filter(HMMhasPolyAA == TRUE)
      polyAA_proteins_vec <- unique(polyAA_proteins$ensembl_peptide_id)
      length(polyAA_proteins_vec)
      
      
      
      universe_df <- bitr(proteins_vec,
                          fromType = fromType,
                          toType = c("ENTREZID", "SYMBOL"),
                          OrgDb = OrgDb)
      universe_genes <- unique(universe_df$ENTREZID)
      
      # Get mapping from ensembl to entrez (entrez ids are needed for enrichment analysis)
      sig_df <- bitr(polyAA_proteins_vec,
                     fromType = fromType,
                     toType = c("ENTREZID", "SYMBOL"),
                     OrgDb = OrgDb)
      
      # Create vector of  entrez ids for significant genes.
      # Not there can be a 1 to many mapping from ensembl to entrez in the above step.
      # thus we need to use 'unique'.
      sig_genes <- unique(sig_df$ENTREZID)
      
      
      # Create vector of entrez ids for a random set of genes that is the same size
      # of the set of significant genes in order to control for sample size. (FDR adjustment should do this,
      # but we are going to add this step just for kicks)
      # Create vector of random control genes with same number as significant genes.
      control_genes <- sample(universe_genes, length(sig_genes))
      print(paste0("Length control genes: ", length(control_genes)))
      
      # Test
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
      
      make_enrichment_plots_and_files <- function(go_enrich_output, type) {
        if (!is.null(go_enrich_output)) {
          go_result <- go_enrich_output@result
          
          if (type == "control") {
            write.csv(go_enrich_output, file = paste0(output_AA_dir,  species, "_", candidate_AA, "_Control_Enriched_GO_", go_ont, ".csv"), row.names = FALSE)
            pdf(paste0(output_AA_dir,  species, "_", candidate_AA, "_Control_Enriched_GO_", go_ont, "_plots.pdf"), width = 8.5, height = 11)
          } else {
            write.csv(go_enrich_output, file = paste0(output_AA_dir, species, "_", candidate_AA, "_Enriched_GO_", go_ont, ".csv"), row.names = FALSE)
            pdf(paste0(output_AA_dir, species, "_", candidate_AA, "_Enriched_GO_", go_ont, "_plots.pdf"), width = 8.5, height = 11)
          }
 
          par(mfrow = c(1, 1))
          p1 <- barplot(go_enrich_output, showCategory=8)
          p2 <- dotplot(go_enrich_output)
          print(p1)
          print(p2)
          if (nrow(go_enrich_output) != 0) {
            p3 <- emapplot(go_enrich_output)
            # p4 <- cnetplot(go_enrich_output, 
            #                categorySize="pvalue", 
            #                foldChange=sig_df$ENTREZID,
            #                node_label = FALSE,
            #                colorEdge = TRUE)
            print(p3)
            # print(p4)
          }
          
          
          category_to_proteins = read.csv(text="GO_id,GO_desc,p.adjust,avg_peptide_length, number_proteins")
          # Iterate over each enriched category, seeing if the length of the proteins in that category are correlated
          # with the p-value of that category
          
          enrich_results <- go_enrich_output@result %>% filter(p.adjust < 0.05) # Should we use a lower threshold?
          if (nrow(enrich_results) > 0) { # Make sure there is actually stuff to iterate over.
            for (i in 1:nrow(enrich_results)) {
              row <- enrich_results[i,]
              cat_gene_syms <- strsplit(row$geneID, split = '/' ) # Gene symbols of genes in category.
              
              # Get the peptide ids of all the genes that were in this enriched category. 
              peptide_ids <- with(universe_df[universe_df$SYMBOL %in% cat_gene_syms[[1]],], get(fromType))
              # Now get their HMM annotations and other data.
              df <- proteins[proteins$ensembl_peptide_id %in% peptide_ids,]
            
              # mean(df$peptide_length())
              mean <- mean(unlist(lapply(as.character(df$peptide_seq), function(x){length(unlist(strsplit(x, "")))})))
              number_proteins <- length(peptide_ids)
              category_to_proteins[i,] <- list(row$ID, row$Description, row$p.adjust, mean, number_proteins)
            }
            
            # Make Avg length of peptide in category vs adjusted p-value of categories.
            # This is a control to make sure that categories with longer peptides aren't 
            # more significant just based on the idea.
            # Are synaptic categories below the line? Play with residuals to check.
            
            # Check if the slopes are significantly different? <- Do this
            # Check if the sum of the residuals is significantly different?
            # What do I do?
            
            # hist(category_to_proteins$number_proteins)
            
            # Make sure category has at least 20 proteins.
            df <- category_to_proteins %>% filter(category_to_proteins$number_proteins > 20)
            
            
            df$synaptic_cat = df$GO_id %in% synapse_cats # True or False, depending on whether or not the category is a synaptic category
            
            # Make two linear regression.
            x = df$avg_peptide_length
            y = df$p.adjust
            fit <- lm(y ~ x)
            
            
            # df$synapse <- ifelse(df$GO_id %in% synapse_cats, TRUE, FALSE)
            p <- ggplot(df, aes(x = avg_peptide_length, y = p.adjust)) + 
              geom_point(alpha = 0.5,  aes(size = number_proteins, color = synaptic_cat)) +
              stat_smooth(method = "lm", col = "red") +
              labs(title = paste("Avg Length of Peptides in Enriched Category VS p-value of Enriched Category \n",
                                 "Adj.R2 = ",signif(summary(fit)$adj.r.squared, 5),
                                 "Intercept =",signif(fit$coef[[1]],5 ),
                                 "  p.val =",signif(summary(fit)$coef[2,4], 5))) +
              xlab("Average Peptide Length of Category") + 
              ylab("Adjusted p-value of Category") +
              scale_size_continuous(name = "Number of Proteins") +
              scale_color_discrete(name = "Synaptic Category") +
              theme_light() +
              theme(plot.title = element_text(hjust = 0.5),
                    legend.title =) +
              # geom_text(aes(label = GO_desc), color ='red', data = df[df$GO_id %in% synapse_cats,]) +
              geom_text_repel(mapping = aes(label = GO_desc), 
                              arrow = arrow(length = unit(0.03, "npc"), type = 'closed', ends ="first"),
                              color ='#00CCCC', 
                              data = df[df$synaptic_cat == TRUE,]) 
            print(p)
          } # for (i in 1:nrow(enrich_results)) {
          dev.off()
        } # if (!is.null(go_enrich_output)) {
      } # make_enrichment_plots_and_files()
      
      make_enrichment_plots_and_files(go_enrich_output_control, type = "control")
      make_enrichment_plots_and_files(go_enrich_output_test, type = "test")
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
    
    pdf(paste0(output_AA_dir, species, "_", candidate_AA, "_", nrow(proteins), "prots_result_plots.pdf"))
      
    #########################################
    # Plot type: Whole Protein Qfraction Histogram for all proteins
    p <- ggplot(data = proteins, aes(x = proteins$AAfraction)) +
      geom_histogram(show.legend = FALSE, aes(fill=..count..)) +
      xlab(paste0("Whole Protein %", candidate_AA)) +
      ylab("# Proteins") +
      ggtitle(paste0("%", candidate_AA, " Whole Protein")) +
      theme_light()
    print(p)
      
    filtered_df <- proteins %>% filter(proteins$AAC> quantile(proteins$AAfraction, 0.95))
    filtered_df <- filtered_df %>% arrange(desc(filtered_df$AAfraction))
    write.csv(filtered_df, file = paste0(output_AA_dir, candidate_AA, "fraction_all_wholeProts_95sig.csv"))
    
    
    # Plot type: Whole Protein AAC Histogram split by HMMhasPolyAA
    p <- ggplot(data = proteins, aes(x = AAfraction, fill = HMMhasPolyAA)) +
      geom_histogram(alpha=0.5, position="identity") +
      xlab(paste0("%", candidate_AA, " (Whole Protein Seqs)")) +
      ylab("# Proteins") +
      labs(fill = paste0("Contains \n poly", candidate_AA)) +
      ggtitle(paste0("%", candidate_AA," Whole Proteins Histogram")) +
      theme_light()
    print(p)
    
    p <- ggplot(data = proteins, aes(x = AAfraction, fill = HMMhasPolyAA)) +
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
    df$RandomSeqs <-unlist(lapply(as.character(df$peptide_seq), rnd_substr, 50))
    nAAs = str_count(df$RandomSeqs, candidate_AA)
    df$AAfractionRandomSeq <- nAAs/nchar(df$RandomSeqs)
    
    
    p <- ggplot(data = df, aes(x = df$AAfractionRandomSeq)) +
      geom_histogram(show.legend = FALSE, aes(fill=..count..)) +
      xlab(paste0("%", candidate_AA, " (Random 50AA Seqs)")) +
      ylab("# of Randomly Choosen Pepitides") +
      ggtitle(paste0("%", candidate_AA, " for 100000 50AA Random Seqs Histogram")) +
      theme_light()
    print(p)
    
    filtered_df <- df %>% filter(df$AAfractionRandomSeq > quantile(df$AAfractionRandomSeq, 0.95))
    filtered_df <- filtered_df %>% arrange(desc(filtered_df$AAfractionRandomSeq))
    write.csv(filtered_df, file = paste0(output_AA_dir, "AAfraction_randomAAs_95sig.csv"))
    
    ########################################################
    # Plot type: Number of Proteins annotated by HMM to be "Poly AA" vs %AA
    
    p <- ggplot(data = df, aes(x = as.numeric(df$AvgPolyAARegionAAFractions))) +
      geom_histogram(show.legend = FALSE, aes(fill=..count..)) +
      xlab(paste0("Avg %", candidate_AA, " of poly", candidate_AA, " Regions")) + 
      ggtitle(paste0("Avg %", candidate_AA, " of poly", candidate_AA," Regions for poly", candidate_AA, " Containing Proteins")) +
      ylab(paste0("# Proteins Containing poly", candidate_AA)) +
      theme_light()
    print(p)
    
    df <- proteins %>% filter(HMMhasPolyAA == TRUE)
    p <- ggplot(data = df, aes(x = df$MaxPolyAARegionAAFractions)) +
      geom_histogram(show.legend = FALSE, aes(fill=..count..)) +
      xlab(paste0("Max %", candidate_AA, " of poly", candidate_AA, " Regions")) +
      ggtitle(paste0("Max %", candidate_AA, " of Poly", candidate_AA, " Regions for poly", candidate_AA, " Containing Proteins")) +
      ylab(paste0("# Proteins Containing poly", candidate_AA)) +
      theme_light()
    print(p)
    
    ##############
    # Plot type 4: Number of Proteins annotated by HMM to be Poly Q vs Length of Region
    
    p <- ggplot(data = df, aes(x = df$AvgLengthsPolyAA)) +
      geom_histogram(show.legend = FALSE, aes(fill=..count..)) +
      xlab(paste0("Average Poly", candidate_AA, " Length")) +
      ylab(paste0("# Proteins Containing poly", candidate_AA, " (by HMM)")) +
      ggtitle(paste0("Average Poly", candidate_AA, " Length in Poly", candidate_AA, " Containing Proteins")) +
      theme_light()
    print(p)
    
    p <- ggplot(data = df, aes(x = df$MaxLengthsPolyAA)) +
      geom_histogram(show.legend = FALSE, aes(fill=..count..)) +
      xlab(paste0("Max poly", candidate_AA, " Length")) +
      ylab(paste0("# Proteins Containing poly", candidate_AA, " (by HMM)")) +
      ggtitle(paste0("Max poly", candidate_AA, " in poly", candidate_AA, " Containing Proteins")) +
      theme_light()
    print(p)
    
    p <- ggplot(data = df, aes(x = df$NumberPolyAA)) + 
      geom_histogram(show.legend = FALSE, aes(fill=..count..)) +
      xlab(paste0("Number of poly", candidate_AA, " Regions")) +
      ylab(paste0("# of ", candidate_AA, " Proteins")) +
      ggtitle(paste0("Number of poly", candidate_AA, " regions in poly", candidate_AA, " Containing Proteins")) +
      theme_light()
    print(p)
    dev.off()
  }
}


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
