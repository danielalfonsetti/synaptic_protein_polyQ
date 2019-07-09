# Daniel Alfonsetti
# daniel.alfonsetti@gmail.com; alfonset@mit.edu
# MIT Littleton Lab
# 7 March, 2019
# Description: Download proteomes for human, fly, mouse, and worm for future polyQ
# analysis in "HMM.R"script.

# Identifier Conventions: 
#   Functions: MyFunction(...)
#   Constants: kMyConstant
#   Variables: myVariables
#   https://google.github.io/styleguide/Rguide.xml

#########################################
# Load libraries and set 'global' variables 
#########################################
rm(list = ls())
source("../../ConstantsAndFunctions.R", chdir=T)

library(biomaRt)
library(org.Dm.eg.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Ce.eg.db)
library(dplyr)

#########################################
# Set constants
#########################################
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 0) {
  kSpecies = unlist(strsplit(args[1],","))
}

kOutputBaseDir <- "../../../output/"
 
#########################################
# Main - Download proteomes
#########################################

for (species in kSpecies) {
  # debugging; species = kSpecies[1]
  print(paste0("Downloading proteome for ", species, "..."))
  
  # listDatasets(datasets); Relevant data sets:
  # mmusculus_gene_ensembl; Drosophila melanogaster
  # hsapiens_gene_ensembl; Homo sapiens
  # dmelanogaster_gene_ensembl; Mus musculus
  # celegans_gene_ensembl; Caenorhabditis elegans
  
  if (species == "fly") {
    ensembl <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")
  } else if (species == "human")  {
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  } else if (species == "mouse") {
    ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  } else if (species == "worm") {
    ensembl <- useMart("ensembl", dataset = "celegans_gene_ensembl")
  }
  
  # l <- listAttributes(ensembl)
  # View(l)
  proteinDf <- getBM(attributes = c('ensembl_gene_id', 'ensembl_peptide_id', "uniprotsptrembl", "external_gene_name", 'go_id', "name_1006", "namespace_1003", "kegg_enzyme"), 
                      values = "*", 
                      mart = ensembl)
  
  # Filter here based on whether a protein is a nuclear protein. Make two seperate files (will do this later)
  tmp <- proteinDf
  tmp$label <- unlist(
    lapply(
      lapply(as.character(tmp$go_id), function(x){unlist(strsplit(x, "; "))}), 
      function(x){any(x %in% kNuclearCats)}
    )
  )
  nuclearProts <- tmp %>% filter(tmp$label)
  
  # Rename columns to make intelligible.
  proteinDf <- dplyr::rename(proteinDf, "go_description" = name_1006 , "go_type" = namespace_1003)
  
  # Paste together GO terms if they are for the same peptide and of the same GO_type.
  proteinDf <- proteinDf %>%  
    distinct() %>%
    group_by(go_type, ensembl_peptide_id) %>% 
    mutate(go_list = paste(go_description, collapse = "; "), go_ids = paste(go_id, collapse = "; "))
  proteinDf <- proteinDf %>% dplyr::select(-c(go_description, go_id)) 
  
  # Put MF, BP, and CC information on same row for a given peptide.
  proteinDf <- proteinDf %>%
    distinct() %>%
    group_by(ensembl_peptide_id) %>%
    mutate(go_list_MF = nth(go_list, 1), go_ids_MF = nth(go_ids, 1),
           go_list_CC = nth(go_list, 2), go_ids_CC = nth(go_ids, 2),
           go_list_BP = nth(go_list, 3), go_ids_BP = nth(go_ids, 3))
  proteinDf <- proteinDf %>% dplyr::select(-c(go_type, go_list, go_ids)) 
  
  vec <- duplicated(proteinDf[c("ensembl_peptide_id")]) # Get rows that share the same peptide ID.
  proteinDf <- proteinDf[!vec,]  # Remove duplicated rows.
  
  # Remove genes that are not protein coding.
  proteinDf <- proteinDf[!(proteinDf$ensembl_peptide_id == ""),]
  
  # Get peptide sequences
  seqs <- getSequence(id = proteinDf$ensembl_peptide_id, 
                      type="ensembl_peptide_id",
                      seqType="peptide",
                      mart=ensembl) 
  
  # Seqs takes the form...
  # peptideSeq | ensembl_peptide_id
  # ------------------------------
  # DSSQSRLRSAH | FBpp0070102
  # MGKKGKGKKGK | FBpp0070068
  
  # Add sequences to the dataframe we are building up
  seqs <- dplyr::rename(seqs , "peptideSeq" = peptide)
  proteinDfSeqs <- merge(proteinDf, seqs, by = "ensembl_peptide_id")
  
  #  Set empty cells equal to NA
  proteinDfSeqs[proteinDfSeqs == ""] <- NA
  
  # Record length of each protein
  proteinDfSeqs$PeptideLength <- unlist(lapply(proteinDfSeqs$peptideSeq, function(x){length(unlist(strsplit(x, "")))}))
  
  # Make two versions: one with all the proteins, and one without nuclear proteins.
  
  # All proteins
  write.csv(proteinDfSeqs, file = paste0(kOutputBaseDir, species, "_prots.csv"), row.names = FALSE)
  
  # All proteins except nuclear proteins.
  proteinDfSeqsFilt <- proteinDfSeqs %>% filter(!(ensembl_peptide_id %in% nuclearProts$ensembl_peptide_id))
  write.csv(proteinDfSeqsFilt, file = paste0(kOutputBaseDir, species, "_prots_nuclear_filt.csv"), row.names = FALSE)
  
  print(paste0("Downloaded and save proteome for ", species, " in ", kOutputBaseDir, "!"))
} # for (species in kSpecies) {

# END
