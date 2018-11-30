# Daniel Alfonsetti, daniel.alfonsetti@gmail.com
# 13 Sep. 2018
# Description: polyQ analysis in synaptic proteins

##########################
library(rentrez)
library(HiddenMarkov)
library(stringr)
# library(dcGOR)
# library(dnet)
# library(supraHex)
# library(Biobase)
# library(Rgraphviz)
library(grid)
# library(domainsignatures)
library(dplyr)
library(ggplot2)
library(doParallel)

require(data.table)

library(biomaRt)

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

##############################################################################
# Helper functions
removeFastaHeader = function(sequence,sep = "\n"){
  b = unlist(strsplit(sequence, split = sep))
  b = b[2:length(b)]
  b = paste0(b, collapse = "")
  return(b)
}


#HIDDEN MARKOV MODEL PARAMETERS
HiddenMarkov = function (proteinSet){
  
  # Debugging
    # proteinSet = Proteins[1:10,]
  #
  
  x = NULL #string translation sequence
  Pi = matrix(c(0.95, 0.05, 0.13, 0.87), nrow = 2, ncol = 2, byrow = TRUE) #hardcode transition prob matrix
  delta = c(.99, .01) #probability of starting in either state
  distn = "binom" #distribution name, average = probability
  pm = list(prob = c(0.05, 0.26), size = 1) #prob is parameter for binom dist, prob of Q in nonQ, Q states
  HMM = dthmm(x, Pi, delta, distn, pm, pn = NULL, discrete = TRUE, nonstat = TRUE) # Create the model
  
  
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
    
    # Debugging
      #sequenceNum = 3
    #
    proteinRow <- proteinSet[sequenceNum,]
    
    
    ### TURN TRANSLATION INTO BINARY VECTOR (coded as 1 if AA == Q, 0 otherwise)
    binarySequence = rep(1, nchar(proteinRow$Translations))
    for (i in 1:length(binarySequence)){
      character = substr(proteinRow$Translations, i, i)
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
      proteinRow$IndiciesPolyQ <- paste0(paste0(grouped_polyQ_indicies), collapse=",")
      
      
      # Get the translations for each polyQ region and store them and their lengths. 
      PolyQTranslations <- character(length(proteinRow$NumberPolyQ)) # Intialize storage vectors
      PolyQLengths <- character(length(proteinRow$NumberPolyQ))
      PolyQRegionQFractions <- numeric(length(proteinRow$NumberPolyQ))
      
      for (i in 1:proteinRow$NumberPolyQ) {
        group <- grouped_polyQ_indicies[[i]] # Indexing into a list is done with double brackets
        
        PolyQTranslations[i] <- substr(proteinRow$Translations, group[1], tail(group, 1))
        PolyQLengths[i] <- nchar(PolyQTranslations[i])
        
        
        nQs = str_count(PolyQTranslations[i], "Q")
        PolyQRegionQFractions[i] <- round(nQs/nchar(PolyQTranslations[i]), 4)
      }
      
      proteinRow$TranslationsPolyQ <- paste0(PolyQTranslations, collapse = ",")
      
      PolyQLengths <- as.numeric(PolyQLengths)
      proteinRow$LengthsPolyQ <- paste0(sort(PolyQLengths, decreasing = TRUE), collapse = ",")
      proteinRow$AvgLengthsPolyQ <- round(mean(PolyQLengths), 4)
      proteinRow$MaxLengthsPolyQ <- PolyQLengths[1]
      
      proteinRow$PolyQRegionQFractions <- paste0(sort(PolyQRegionQFractions, decreasing = TRUE), collapse = ",")
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


formatter_helper_func <- function(recs) {
  
  Proteins = foreach(i = 1:length(recs), .combine = rbind)  %do% {
    x <- recs[i][[1]]
    c(x$caption, x$uid, x$title)
  }
  Proteins <- as.data.frame(Proteins)
  names(Proteins) <- c("AccessionNumber", "ENTREZID", "ProteinName") # If we can replace these three things from bioMart, we are set.
  
  
  Translations = entrez_fetch("Protein", id = Proteins$AccessionNumber, rettype = "fasta")
  
  splitTranslations <- strsplit(Translations, split = "\n>")[[1]]
  TranslationsList = unlist(lapply(splitTranslations, removeFastaHeader))
  
  Proteins$Translations = TranslationsList
  
  
  nQs = str_count(Proteins$Translations, "Q")
  Proteins$Qfraction = round(nQs/nchar(Proteins$Translations), 4)
  
  Proteins <- HiddenMarkov(Proteins)
  return(Proteins) 
}

########## End of helper functions ###############

############################################################
# Download and get the data.
############################################################

# Debugging: manual download
  # seq_start = 5
  # search <- entrez_search(db="protein", term="Drosophila melanogaster[ORGN]", retmax = 135664, web_history = TRUE)
  # recs <- entrez_summary(db="protein", id = search$ids, retmax=5, retstart=seq_start)
  # 
  # result <- formatter_helper_func(recs)
  # write.csv(result, file = "output.csv", row.names = FALSE)
#########################

  
start_time <- Sys.time() # Measure download time

# Automatic download
number_proteins = 135664  # Number of proteins to be analyzed. 135664 is all proteins in the Drosophilia melanogaster proteome.
file_size = 100

foreach(group_start = seq(1,135664, file_size)) %dopar% {
  library(rentrez)
  library(HiddenMarkov)
  library(stringr)
  library(grid)
  library(dplyr)
  library(ggplot2)
  library(doParallel)
  library(reshape2)
  library(biomaRt)
  
  search <- entrez_search(db="protein", term="Drosophila melanogaster[ORGN]", retmax = 135666, use_history = TRUE)
  
  search$ids # Entrez ids
  
  
  ensembl <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")
  seqs <- getSequence(id = search$ids[1:100], 
                      type="entrezgene",
                      seqType="peptide",
                      mart=ensembl) 
  
  
  
  ######################################################
  # http://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/Genome%20Annotation%20and%20Visualisation%20using%20R%20and%20Bioconductor.pdf
  mart <- useMart(biomart = "")

  
  ######################################################
  
  
  proteins = foreach(seq_start = seq(group_start, group_start+file_size-1, 100), .combine = rbind) %dopar% {
    recs <- entrez_summary(db="protein",  retmax=100, retstart=seq_start, web_history = search$web_history)

    # recs <- entrez_summary(db="protein", id = search$ids, retmax=100, retstart=seq_start)
    
    protein_set = formatter_helper_func(recs) # Reformats output into a dataframe, makes new columns, and runs hidden markov model.
    
    ## Testing
    
    
    protein_set
  }
  
  proteins <- as.data.frame(proteins)
  write.csv(proteins, file = paste0(data_dir, "proteins_N", group_start, "_", group_start+file_size-1, ".csv"), row.names = FALSE)
  
  
  
  # Add go Annotations
  
  head(listFilters(ensembl),20) # Stuff we can filter by
  head(listAttributes(ensembl), 20)# stuff we can retrive
  
  
  columns(org.Dm.eg.db)
  select(org.Dm.eg.db,
         keys=as.character(protein_set$AccessionNumber), keytype="ACCNUM",
         columns = c("GO", "SYMBOL", "ALIAS", "GENENAME", "UNIGENE", "UNIPROT"))
  
  
  # Using the entrenz stuff
  
  
  the_keys <- as.list(org.Dm.egACCNUM)
  
  
  ####################

  # Convert to a list
  Entrez_ids <- names(as.list(org.Dm.egFLYBASEPROT)) #Entrez ids
  length(Entrez_ids)
  
  seqs <- getSequence(id = Entrez_ids, 
              type="entrezgene",
              seqType="peptide",
              mart=ensembl) 

  as.list(org.Dm.egENSEMBL)
  
  ensemble_prot_ids <- unlist(as.list(org.Dm.egENSEMBLPROT))
  
  seqs <- getSequence(id = ensemble_prot_ids[1:10], 
                      type= "ensembl_peptide_id",
                      seqType="peptide",
                      mart=ensembl) 
  
 
  columns(org.Dm.eg.db)
  result <- select(org.Dm.eg.db, keys=as.character(seqs$entrezgene), keytype="ENTREZID", columns = c("GO"))
  result <- result[complete.cases(result),]
  
  
  #######3 Stuff below works
  #result <- select(org.Dm.eg.db, keys=as.character(proteins$AccessionNumber), keytype="ACCNUM", columns = c("GO"))
  a1 <- result %>%  
    dplyr::select(-EVIDENCE) %>%
    distinct() %>%
    group_by(ONTOLOGY, ENTREZID) %>% 
    mutate(GO_list = paste(GO, collapse = ",")) %>%
    dplyr::select(-GO) %>%
    distinct() %>%
    dcast(ENTREZID ~ ONTOLOGY)
  View(a1)
  
  final <- merge(a1, proteins, by.x = "ACCNUM", by.y = "AccessionNumber")
  
}


################################
# Stitch together files
setwd(data_dir)
file_names <- dir()


# Put all files (dataframes) into a list
df_list = list()
i <- 1
for (file_name in file_names) {
  df <- read.csv(file_name) # Read in data
  df_list[[i]] <- df # add it to the list
  i <- i + 1
}

# Merge the list
big_data = do.call(rbind, df_list)
setwd(output_dir)
write.csv(big_data, file = "proteins_v1.csv")
proteins <- big_data


proteins$Translations <- as.character(proteins$Translations)


############################
end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)
############################


##########################
# Load in protein set data from local storage if you didn't download the sets from above.

  proteins <- read.csv(paste0(output_dir, "proteins_v1.csv"))
  proteins$Translations <- as.character(proteins$Translations)
  proteins$Qfraction <- as.numeric(proteins$Qfraction)
  # Only 47569 unique translations
#####


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
      
