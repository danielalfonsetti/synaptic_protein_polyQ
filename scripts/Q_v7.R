# Daniel Alfonsetti, daniel.alfonsetti@gmail.com
# 13 Sep. 2018
# Description: polyQ analysis in synaptic proteins

##########################
library(HiddenMarkov)
library(stringr)
library(grid)
library(dplyr)
library(ggplot2)


library(doParallel)
library(data.table)

library(biomaRt)
library(biomartr)

library(clusterProfiler)
library(STRINGdb)
library(HiddenMarkov)
library(org.Dm.eg.db) # Drosophila_melanogaster
library(org.Hs.eg.db) 
library(org.Mm.eg.db)
library(org.Ce.eg.db)
rm(list = ls())


##########################
# Setup environment
# 
registerDoParallel(cores = 4)

setwd("C:/UROPs/polyQ_neuronal_proteins")
plot_output_dir <- "C:/UROPs/polyQ_neuronal_proteins/plots/"
data_dir <- "C:/UROPs/polyQ_neuronal_proteins/data/"

output_base_dir <- "C:/UROPs/polyQ_neuronal_proteins/output/"
###########

# User defined variables
candidate_AA_vec = c("Q")
species_vec = c("fly", "mouse", "human", "worm") # choose from 'fly', 'human', 'mouse' and 'worm'

##########################################3
##############################################################################


#HIDDEN MARKOV MODEL PARAMETERS
HiddenMarkov = function (proteinSet, candidate_AA){
  
  # Debugging
  # proteinSet = seqs[6830:6840,]
  #
  
  x = NULL #string translation sequence
  Pi = matrix(c(0.95, 0.05, 0.13, 0.87), nrow = 2, ncol = 2, byrow = TRUE) #hardcode transition prob matrix
  delta = c(.99, .01) #probability of starting in either state (nonQ region, Q region)
  distn = "binom" # distribution name, average = probability
  pm = list(prob = c(0.05, 0.26), size = 1) #prob is parameter for binom dist, prob of Q in nonQ, Q states
  HMM = dthmm(x, Pi, delta, distn, pm, pn = NULL, discrete = TRUE, nonstat = TRUE) # Create the model
  
  
  # A column for the Qfraction for an entire protein.
  proteinSet$AAfraction = rep(NA, nrow(proteinSet))
  
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
  
  
  #proteinSet$binarySequence = rep(nrow(proteinSet))
  #proteinSet$Annotation = rep(nrow(proteinSet))
  for (sequenceNum in 1: nrow(proteinSet)){
    print(sequenceNum)
    # Debugging
    #sequenceNum = 6833
    #
    
    proteinRow <- proteinSet[sequenceNum,]
    
    # Get Qfraction for the entire protein.
    nAAs = str_count(proteinRow$peptide_seq, candidate_AA)
    proteinRow$AAfraction = round(nAAs/nchar(proteinRow$peptide_seq), 4)
    
    
    ### TURN TRANSLATION INTO BINARY VECTOR (coded as 1 if AA == Q, 0 otherwise)
    binarySequence = rep(1, nchar(proteinRow$peptide_seq))
    for (i in 1:length(binarySequence)){
      character = substr(proteinRow$peptide_seq, i, i)
      if (character == candidate_AA){
        binarySequence[i] = 1
      } else { binarySequence[i] = 0}
    }
    
    if (length(binarySequence) > 1) {
      ### ADD BINARY SEQUENCE VECTOR TO HMM INPUTS
      HMM$x = binarySequence
      ### CREATING AN ANNOTATED (STATES) SEQUENCE FOR THE TRANSLATION
      annotatedSequence = Viterbi(HMM)
    } else {
      # Viterbi doesn't work if length of the sequence is 1, so just manually do it.
        annotatedSequence = c(1)
    }

    ### DETERMINE IF THERE ARE ANY Q ISLANDS AND REPORT TRUE OR FALSE IN THE TABLE
    a = annotatedSequence == 2
    if (sum(a) >= 1 ) {
    } else {proteinRow}
    
    # Find which indices of the translation are marked as being polyQ
    polyAA_indices <- which(annotatedSequence %in% c(2))
    
    # Group the indicies into groups (so if you have multiple polyQ regions, the indices for each polyQ region will be seperated from each other)
    grouped_polyAA_indicies <- split(polyAA_indices, cumsum(c(1, diff(polyAA_indices) != 1))) 
    
    
    # Record how many polyQ regions there are
    if (length(grouped_polyAA_indicies[[1]]) == 0) {
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
      }
      
      proteinRow$TranslationsPolyAA <- paste0(PolyAATranslations, collapse = "; ")
      
      PolyAALengths <- as.numeric(PolyAALengths)
      proteinRow$LengthsPolyAA <- paste0(sort(PolyAALengths, decreasing = TRUE), collapse = "; ")
      proteinRow$AvgLengthsPolyAA <- round(mean(PolyAALengths), 4)
      proteinRow$MaxLengthsPolyAA <- PolyAALengths[1]
      
      proteinRow$PolyAARegionAAFractions <- paste0(sort(PolyAARegionAAFractions, decreasing = TRUE), collapse = "; ")
      proteinRow$AvgPolyAARegionAAFractions <- round(mean(PolyAARegionAAFractions), 4)
      proteinRow$MaxPolyAARegionAAFractions <- PolyAARegionAAFractions[1]
    }
    
    # Add the row back
    proteinSet[sequenceNum,] <- proteinRow
  }
  
  # Return a column that gives length of largest region 
  # and the number of regions present (and the Qfraction of the largest region?)
  
  return (proteinSet)
}

########## End of helper functions ###############


for (species in species_vec) {
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

  
  
  protein_df <- getBM(attributes = c('ensembl_gene_id', 'ensembl_peptide_id', "external_gene_name", 'go_id', "name_1006", "namespace_1003"), 
                      values = "*", 
                      mart = ensembl)
  protein_df <- rename(protein_df, "go_description" = name_1006 , "go_type" = namespace_1003 )
  
  # Paste together GO terms if they are for the same peptide and of the same GO_type.
  protein_df <- protein_df %>%  
    distinct() %>%
    group_by(go_type, ensembl_peptide_id) %>% 
    mutate(go_list = paste(go_description, collapse = "; "), go_ids = paste(go_id, collapse = "; "))
  protein_df <- protein_df %>% dplyr::select(-c(go_description, go_id)) 
  
  vec <- duplicated(protein_df[,c("go_type", "ensembl_peptide_id")]) # Get rows that share same GO type and peptide ID.
  protein_df <- protein_df[!vec,]  # Remove duplicated rows.
  
  
  #############
  # Do the filtering for nuclear proteins here.
  # [FILTERING CODE]
  ##############
  
  write.csv(protein_df, file = paste0(output_species_dir, species, "_", nrow(protein_df), "prots.csv"), row.names = FALSE)
  
  seqs <- getSequence(id = protein_df$ensembl_peptide_id, 
                      type="ensembl_peptide_id",
                      seqType="peptide",
                      mart=ensembl) 

  # Seqs takes the form...
  
  # peptide_seq | ensembl_peptide_id
  # ------------------------------
  # DSSQSRLRSAH | FBpp0070102
  # MGKKGKGKKGK | FBpp0070068
  
  seqs <- rename(seqs , "peptide_seq" = peptide)
  
  

  for (candidate_AA in candidate_AA_vec) {
    print(paste0("Species: ", species,"; Amino acid: ", candidate_AA))
    output_AA_dir <- paste0(output_species_dir, candidate_AA, "/")
    dir.create(output_AA_dir, recursive = TRUE)
    
    markov_result <- HiddenMarkov(seqs, candidate_AA)
    
    protein_df_w_seqs <- merge(protein_df, markov_result, by = "ensembl_peptide_id")
    protein_df_w_seqs <- protein_df_w_seqs %>% dplyr::select("ensembl_gene_id", everything())
    write.csv(protein_df_w_seqs, file = paste0(output_AA_dir, species, "_", nrow(protein_df_w_seqs), "prots_w_HMM_", candidate_AA, ".csv"), row.names = FALSE)
    
    
    ############################
    end_time <- Sys.time()
    elapsed_time <- end_time - start_time
    print(elapsed_time)
    ############################
    
    ########################################################3
    
    # Load in protein set data from local storage if you didn't download the sets from above
    proteins <- read.csv( paste0(output_AA_dir, species, "_", nrow(protein_df_w_seqs), "prots_w_HMM_", candidate_AA, ".csv"))
    #####
    
    
    #########################################################
    #########################################################
    
    polyAA_proteins = proteins %>% filter(proteins$HMMhasPolyAA == TRUE)
    
    proteins_unique <- unique(proteins$ensembl_peptide_id)
    polyAA_proteins_unique <- unique(polyAA_proteins$ensembl_peptide_id)
    length(proteins_unique) # 30284
    length(polyAA_proteins_unique) # 9713
    
  
    
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
    
    universe_df <- bitr(proteins_unique,
                        fromType = fromType,
                        toType = c("ENTREZID", "UNIPROT"),
                        OrgDb = OrgDb)
    
    # Create df that contains only significant proteins (HMM proteins)
    sig_df <- bitr(polyAA_proteins_unique,
                   fromType = fromType,
                   toType = c("ENTREZID", "UNIPROT"),
                   OrgDb = OrgDb)
    
    
    sig_df <- bitr(as.factor(ensembl_peptide_id),
                   fromType = fromType,
                   toType = c("ENTREZID", "UNIPROT"),
                   OrgDb = OrgDb)
    
    
    
    # ###########################################
    # clusterProfiler Part 1 - GO
    # ###########################################
    
    for (go_ont in c("BP", "MF", "CC")) {
      print(paste0("Doing enrichment analysis for ", go_ont))
      
      # Control for sample size: Random subset of universe df with same size as
    
      # Control  
      control_df <- sample_n(universe_df, nrow(sig_df))
      go_enrich_output_control <- enrichGO(gene = control_df$ENTREZID,
                                   universe = universe_df$ENTREZID,
                                   OrgDb = OrgDb,
                                   ont = go_ont,
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05,
                                   qvalueCutoff = 0.05,
                                   readable = TRUE)
      # Test
      go_enrich_output <- enrichGO(gene = sig_df$ENTREZID,
                                   universe = universe_df$ENTREZID,
                                   OrgDb = OrgDb,
                                   ont = go_ont,
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05,
                                   qvalueCutoff = 0.05,
                                   readable = TRUE)
      
      

      
      
      
      
      
      # Control plots
      if (!is.null(go_enrich_output_control)) {
        go_control_result <- go_enrich_output_control@result
        
        write.csv(go_enrich_output_control, file = paste0(output_AA_dir,  species, "_", candidate_AA, "_Control_Enriched_GO_", go_ont, ".csv"), row.names = FALSE)
        pdf(paste0(output_AA_dir,  species, "_", candidate_AA, "_Control_Enriched_GO_", go_ont, "_plots.pdf"), width = 8.5, height = 11)
        
        par(mfrow = c(1, 1))
        p1 <- barplot(go_enrich_output_control, showCategory=8)
        p2 <- dotplot(go_enrich_output_control)
        
        print(p1)
        print(p2)
        
        if (nrow(go_enrich_output_control) != 0) {
          p3 <- emapplot(go_enrich_output_control)
          p4 <- cnetplot(go_enrich_output_control, 
                         categorySize = "pvalue", 
                         foldChange = control_df$ENTREZID,
                         node_label = FALSE,
                         colorEdge = TRUE) 
          print(p3)
          print(p4)
        }
      

        dev.off()
      }
      
      # Treatment plots
      if (!is.null(go_enrich_output)) {
        go_result <- go_enrich_output@result
        
        write.csv(go_enrich_output, file = paste0(output_AA_dir, species, "_", candidate_AA, "_Enriched_GO_", go_ont, ".csv"), row.names = FALSE)
        pdf(paste0(output_AA_dir, species, "_", candidate_AA, "_Enriched_GO_", go_ont, "_plots.pdf"), width = 8.5, height = 11)
        
        par(mfrow = c(1, 1))
        p1 <- barplot(go_enrich_output, showCategory=8)
        p2 <- dotplot(go_enrich_output)
        print(p1)
        print(p2)
        if (nrow(go_enrich_output) != 0) {
          p3 <- emapplot(go_enrich_output)
          p4 <- cnetplot(go_enrich_output, 
                         categorySize="pvalue", 
                         foldChange=sig_df$ENTREZID,
                         node_label = FALSE,
                         colorEdge = TRUE)
          print(p3)
          print(p4)
        }
        dev.off()
      }
      
      
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
    
    
    
    ######### ClusterProfiler ##################
    # Use clusterProfiler to compare proteins with polyAA to all proteins
    
    
    
    ############################################################################3
    ################################## Plotting ##################################
    ############################################################################
      
    # Number of randomly chosen 50 amino acid peptide chunks [WORKS]
    
    # Plotting helper function
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
    pdf(paste0(output_AA_dir, species, "_", candidate_AA, "_", nrow(proteins), "prots_result_plots.pdf"))
      
    ########################################################
    #######################
    # Plot type: Whole Protein Qfraction Histogram for all proteins
    p <- ggplot(data = proteins, aes(x = proteins$AAfraction)) +
      geom_histogram(show.legend = FALSE, aes(fill=..count..)) +
      xlab(paste0("Whole Protein %", candidate_AA)) +
      ylab("# Proteins") +
      ggtitle(paste0("%", candidate_AA, " Whole Protein")) +
      theme_light()
    print(p)
      
    filtered_df <- proteins %>% filter(proteins$AAfraction > quantile(proteins$AAfraction, 0.95))
    filtered_df <- filtered_df[-order(filtered_df$AAfraction),]
    write.csv(filtered_df, file = paste0(output_AA_dir, candidate_AA, "fraction_all_wholeProts_95sig.csv"))
    
    
    # Plot type: Whole Protein AAfraction Histogram split by HMMhasPolyAA
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
    filtered_df <- filtered_df[-order(filtered_df$AAfractionRandomSeq),]
    write.csv(filtered_df, file = paste0(output_AA_dir, "AAfraction_randomAAs_95sig.csv"))
    
    ########################################################
    # Plot type: Number of Proteins annotated by HMM to be "Poly AA" vs %AA
    
    p <- ggplot(data = df, aes(x = as.numeric(df$AvgPolyAARegionAAFractions))) +
      geom_histogram(show.legend = FALSE, aes(fill=..count..)) +
      xlab(paste0("Avg %", candidate_AA, " of poly", candidate_AA, " Regions")) + 
      ggtitle(paste0("Avg %", candidate_AA, " of poly", candidate_AA," Regions for poly", candidate_AA, " Containing Proteins Histogram")) +
      ylab(paste0("# Proteins Containing poly", candidate_AA)) +
      theme_light()
    print(p)
    
    df <- proteins %>% filter(HMMhasPolyAA == TRUE)
    p <- ggplot(data = df, aes(x = df$MaxPolyAARegionAAFractions)) +
      geom_histogram(show.legend = FALSE, aes(fill=..count..)) +
      xlab(paste0("Max %", candidate_AA, " of poly", candidate_AA, " Regions")) +
      ggtitle(paste0("Max %", candidate_AA, " of Poly", candidate_AA, " Regions for poly", candidate_AA, " Containing Proteins Histogram")) +
      ylab(paste0("# Proteins Containing poly", candidate_AA)) +
      theme_light()
    print(p)
    
    ##############
    # Plot type 4: Number of Proteins annotated by HMM to be Poly Q vs Length of Region
    
    p <- ggplot(data = df, aes(x = df$AvgLengthsPolyAA)) +
      geom_histogram(show.legend = FALSE, aes(fill=..count..)) +
      xlab(paste0("Average Poly", candidate_AA, " Length")) +
      ylab(paste0("# Proteins Containing poly", candidate_AA, " (by HMM)")) +
      ggtitle(paste0("Average Poly", candidate_AA, " Length in Poly", candidate_AA, " Containing Proteins Histogram")) +
      theme_light()
    print(p)
    
    p <- ggplot(data = df, aes(x = df$MaxLengthsPolyAA)) +
      geom_histogram(show.legend = FALSE, aes(fill=..count..)) +
      xlab(paste0("Max poly", candidate_AA, " Length")) +
      ylab(paste0("# Proteins Containing poly", candidate_AA, " (by HMM)")) +
      ggtitle(paste0("Max poly", candidate_AA, " in poly", candidate_AA, " Containing Proteins Histogram")) +
      theme_light()
    print(p)
    dev.off()
  }
}


############ End Plotting ################



#################### Misc. Code #################################


# ggplot(data = proteins, aes(x = proteins$Qfraction, y = ..count../sum(..count..))) +
#   geom_density() +
#   xlab("Q-fraction (Whole Protein Seqs)") +
#   ylab("Percent") +
#   theme_light()

# clusterProfiler needs entrenz id

# Add go Annotations

# head(listFilters(ensembl),20) # Stuff we can filter by
# head(listAttributes(ensembl), 20)# stuff we can retrive

# columns(org.Dm.eg.db)
# select(org.Dm.eg.db,
#        keys=as.character(protein_set$AccessionNumber), keytype="ACCNUM",
#        columns = c("GO", "SYMBOL", "ALIAS", "GENENAME", "UNIGENE", "UNIPROT"))


# getCollection( db = "ensembl",
#                organism = "Drosophila melanogaster",
#                path = file.path("_ncbi_downloads","collection"))
# 
# getCollection( db = "refseq",
#                organism = "Drosophila melanogaster",
#                path = file.path("refseq","collection"))

# Read the protein sequence file, looks like:
# $A06852
# [1] "M" "P" "R" "L" "F" "S" "Y" "L" "L" "G" "V" "W" "L" "L" "L" "S" "Q" "L"
# ...

      
