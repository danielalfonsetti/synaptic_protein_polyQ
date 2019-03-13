# Daniel Alfonsetti
# daniel.alfonsetti@gmail.com; alfonset@mit.edu
# MIT Littleton Lab
# 7 March, 2019
# Description: Download proteomes for human, fly, mouse, and worm for future polyQ
# analysis in "HMM.R"script.
#########################################
# Load libraries and set 'global' variables 
#########################################
rm(list = ls())

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

output_base_dir <- "C:/UROPs/polyQ_neuronal_proteins/output/"
species_vec = c("fly")

# Transcription factor categories
term1 <- "GO:0003676" # nucleic acid binding
transcription_cats <- c(term1, get(term1, GOMFOFFSPRING))

term1 <- "GO:0005634" # Nucleus
term2 <- "GO:0009295" # nucleoid 
nuclear_cats <- c(term1, get(term1, GOCCOFFSPRING), 
                  term2, get(term2, GOCCOFFSPRING),
                  transcription_cats)

#########################################
# Helper function for filtering certain GO categories
#########################################
annotate_cats <- function(df, cats) {
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
  df$label <- CC | MF | BP
  return(df)
}


#########################################
# Main - Download proteomes
#########################################

for (species in species_vec) {
  # debugging; species = species_vec[1]
  print(paste0("Downloading proteome for ", species, "..."))
  output_species_dir <- paste0(output_base_dir, species, "/")
  dir.create(output_species_dir, recursive = TRUE)
  
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
  
  listAttributes(ensembl)
  protein_df <- getBM(attributes = c('ensembl_gene_id', 'ensembl_peptide_id', "external_gene_name", 'go_id', "name_1006", "namespace_1003", "kegg_enzyme"), 
                      values = "*", 
                      mart = ensembl)
  
  
  # protein_df <- getBM(attributes = c("kegg_enzyme"),
  #                     values = "*",
  #                     mart = ensembl)
  
  # Filter here based on whether a protein is a nuclear protein. Make two seperate files (will do this later)
  tmp <- protein_df
  tmp <- annotate_cats(tmp, nuclear_cats)
  nuclear_prots <- tmp %>% filter(tmp$label)
  
  # Rename columns.
  protein_df <- dplyr::rename(protein_df, "go_description" = name_1006 , "go_type" = namespace_1003)
  
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
  
  # Make two versions: one with all the proteins, and one without nuclear proteins.
  write.csv(protein_df_w_seqs, file = paste0(output_species_dir, species, "_prots.csv"), row.names = FALSE)
  
  # Something is sketchy here with the nuclear filters....
  # Only include non-nuclear proteins
  protein_df_w_seqs_filt <- protein_df_w_seqs %>% filter(!(ensembl_peptide_id %in% nuclear_prots$ensembl_peptide_id))
  write.csv(protein_df_w_seqs_filt, file = paste0(output_species_dir, species, "_prots_nuclear_filt.csv"), row.names = FALSE)
  
  #####################
  print(paste0("Downloaded and saved proteome for ", species, " in ", output_species_dir, "!"))
}

##########################