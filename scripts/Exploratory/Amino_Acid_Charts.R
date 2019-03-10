####################################################
# Daniel Alfonsetti, daniel.alfonsetti@gmail.com
# MIT, Littleton Lab UROP
# 27 February 2019
# ---------------------------------
# Description: 
####################################################
# Import libraries and set 'global' variables
####################################################
rm(list = ls())
library(GO.db)
library(biomaRt)
library(biomartr)
library(ggplot2)
library(reshape2)
library(dplyr)

output_base_dir <- "C:/UROPs/polyQ_neuronal_proteins/output/"

######################################################
# Helper functions
######################################################

# Helper function for plotting
plot_AA_Chart <- function(row, graph_HMM_annots = FALSE)  {
  peptide_seq <- as.character(row$peptide_seq)
  peptide_name <- row$ensembl_gene_id
  
  candidate_AA_vec = c("A", "V", "L", "I", "M", "F","W", "C", "Y", "S", "T", "G", "P", "Q","N", "E", "D", "R", "H", "K")
  
  
  # Get indicies of each amino acid in the protein
  AA_pos_list <- lapply(
    candidate_AA_vec,
    function(AA) 
      unlist(lapply(
        strsplit(peptide_seq, ''), 
        function(peptide_seq) 
          which(peptide_seq == AA)
      ))
  )
  
  ########## Convert list to dataframe ############3
  ## Compute maximum length
  max.length <- max(sapply(AA_pos_list, length))
  ## Add NA values to list elements
  AA_pos_list <- lapply(AA_pos_list, function(v) { c(v, rep(NA, max.length-length(v)))})
  ## Cbind
  df <- do.call(cbind, AA_pos_list)
  colnames(df) <- candidate_AA_vec
  
  # Reorganize into a form useful for plotting
  melted_df <- melt(df)[,c(-1)]

  
  box_width = 7
  p <- ggplot() +
    # geom_point(shape = "|", size = 4.5, alpha = 0.3, aes(x=melted_df$value, y=melted_df$Var2)) +
    geom_segment(aes(x=melted_df$value-box_width/2, xend = melted_df$value+box_width/2, 
                     y=melted_df$Var2, yend = melted_df$Var2),
                 alpha = 1/box_width,
                 size = 8) +
    xlab("Peptide Position Index") +
    ylab("Amino Acid") +
    ggtitle(paste0("Amino acid chart for ", row$external_gene_name, " (", row$ensembl_peptide_id, ")")) +
    theme_light() +
    theme(legend.position = "none")

  if (graph_HMM_annots) {
    HMM_annotations_list <- lapply(
      candidate_AA_vec,
      function(AA) 
        unlist(lapply(strsplit(as.character(eval(parse(text = paste0("row$IndiciesPoly", AA)))), split = "; ")[[1]], function(i) eval(parse(text = i))))
    )
    # Remove first and last index in each list (because we plot each point as a line from its middle position to -1 and +1).
    HMM_annotations_list <- lapply(HMM_annotations_list, function(vec){vec[2:(length(vec)-1)]})
    
    # Find the vector with the max length, and pad the rest of the vectors so that we can convert the list to the dataframe.
    max.length <- max(sapply(HMM_annotations_list, length))
    HMM_annotations_list <- lapply(HMM_annotations_list, function(v) {if (is.na(v[1])) {rep(NA, max.length)} else {c(v, rep(NA, max.length-length(v)))}})
    
    df <- do.call(cbind, HMM_annotations_list)
    colnames(df) <- candidate_AA_vec
    melted_HMM_annots_df <- melt(df)[,c(-1)]
    if (!all(is.na(melted_HMM_annots_df$value))) {
      #p <- p + geom_point(shape = "|", size = 1, color = "green", aes(x=, y))
      p <- p + geom_segment(aes(x=melted_HMM_annots_df$value-1, xend = melted_HMM_annots_df$value+1, 
                                  y=melted_HMM_annots_df$Var2, yend = melted_HMM_annots_df$Var2),
                              alpha = 1, color = "Green", size = 2.5) 
    }
    
  }
  print(p)
}


# Helper function for filtering certain GO categories
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


####################################################
# Preliminary Data Wrangling
####################################################

candidate_AA_vec = c("K", "H", "R","D", "E", "N", "Q", "P", "G", "T", "S", "Y", "C", "W", "F", "M", "I", "L", "V", "A")
# Merge the polyAA output files for each HMM type.
for (HMM_type in c("adjusted", "trained")) {
  # Important columns to take from each file are...
  # First 11 columns (general protein information, not polyAA annotation specific)
  # And 15th column (indices of polyAA regions)
  
  
  merged_polyAA_df <- read.csv(paste0("C:/UROPs/polyQ_neuronal_proteins/output/",  HMM_type, "/fly/", candidate_AA_vec[1], "/fly_prots_w_HMM_", candidate_AA_vec[1], ".csv"))
  merged_polyAA_df <- merged_polyAA_df[, c(1:11, 15)]
  new_name <- paste0("IndiciesPoly", candidate_AA_vec[1])
  colnames(merged_polyAA_df)[colnames(merged_polyAA_df)=="IndiciesPolyAA"] <- new_name
  
  # Iterate over the rest of the AAs and append their results.
  for (candidate_AA in candidate_AA_vec[2:length(candidate_AA_vec)]) {
    df <- read.csv(paste0("C:/UROPs/polyQ_neuronal_proteins/output/", HMM_type, "/fly/", candidate_AA, "/fly_prots_w_HMM_", candidate_AA, ".csv"))
    
    new_name <- paste0("IndiciesPoly", candidate_AA)
    colnames(df)[colnames(df)=="IndiciesPolyAA"] <- new_name
    merged_polyAA_df = merge(x = merged_polyAA_df, y = df[, c("ensembl_peptide_id",new_name)], by = "ensembl_peptide_id")
  }
  write.csv(merged_polyAA_df, paste0("C:/UROPs/polyQ_neuronal_proteins/output/", HMM_type, "/fly/merged_polyAA_df.csv"), row.names = FALSE)
}

####################################################
# Plotting
####################################################


####################################################
# Plotting 1 - Explore the AA charts for some hand chosen proteins
####################################################

for (HMM_type in c("adjusted")) {
  proteins <- read.csv(paste0("C:/UROPs/polyQ_neuronal_proteins/output/", HMM_type, "/fly/merged_polyAA_df.csv"))
  
  # Known AZ proteins
  AZ_proteins <- proteins %>% filter(proteins$external_gene_name %in% c("brp", "Rbp", "Rim"))
  AZ_proteins <- AZ_proteins[!duplicated(AZ_proteins$ensembl_gene_id),]
  
  dir.create(paste0("C:/UROPs/polyQ_neuronal_proteins/output/",HMM_type,"/fly/AA_charts"))
  setwd(paste0("C:/UROPs/polyQ_neuronal_proteins/output/",HMM_type,"/fly/AA_charts"))
  pdf(file = paste0("fly_AA_charts_manualAZ_protein_Model-", HMM_type, "_fly.pdf"))
  for (i in 1:nrow(AZ_proteins)){
    row <- AZ_proteins[i,]
    plot_AA_Chart(row, TRUE)
  }
  dev.off()
  
  
  # Known PSD proteins
  # "homer", "Grip", "Dlg", "Sh3", "shk"
  PSD_proteins <- proteins %>% filter(proteins$external_gene_name %in% c("homer", "Grip", "dlg1", "SH3PX1", "Prosap"))
  PSD_proteins <- PSD_proteins[!duplicated(PSD_proteins$ensembl_gene_id),]

  pdf(file = paste0("AA_charts_manualPSD_protein_Model-", HMM_type, "_fly.pdf"))
  for (i in 1:nrow(PSD_proteins)){
    row <- PSD_proteins[i,]
    plot_AA_Chart(row, TRUE)
  }
  dev.off()
  
  
  # Proteins used in HMM training set
  training_set_protein_ids <- c("FBpp0289769","FBpp0307700", "FBpp0086727", "FBpp0111724", "FBpp0293366",
                                "FBpp0070830", "FBpp0309352", "FBpp0402897", "FBpp0110299","FBpp0305807")
  training_set_proteins <- proteins %>% filter(proteins$ensembl_peptide_id %in% training_set_protein_ids)
  pdf(file = paste0("AA_charts_training_set_proteins_Model-", HMM_type, "_fly.pdf"))
  for (i in 1:nrow(training_set_proteins)){
    row <- training_set_proteins[i,]
    plot_AA_Chart(row, TRUE)
  }
  dev.off()
}


####################################################
# Plotting 2 - Explore AA charts for proteins of certain GO categories
####################################################


# Synaptic categories
synapse_term1 <- "GO:0045202"
synapse_cats <- get(synapse_term1, GOCCOFFSPRING)
synapse_cats <- c(synapse_cats, synapse_term1)

# Active Zone categories
AZ_term1 <- "GO:0048786"
AZ_cats <- get(AZ_term1, GOCCOFFSPRING)
AZ_cats <- c(AZ_cats, AZ_term1)

# Post synaptic density categories
PSD_term1 <- "GO:0014069"
PSD_cats <- get(PSD_term1, GOCCOFFSPRING)
PSD_cats <- c(PSD_cats, PSD_term1)

# Transcription factor categories
term1 <- "GO:0003676" # nucleic acid binding
transcription_cats <- c(term1, get(term1, GOMFOFFSPRING))

# Nucleus categories 
term1 <- "GO:0005634" # Nucleus
term2 <- "GO:0009295" # nucleoid 
nuclear_cats <- c(term1, get(term1, GOCCOFFSPRING), 
                  term2, get(term2, GOCCOFFSPRING),
                  transcription_cats)

location_dict = list("AZ" = AZ_cats, "Synapse" =  synapse_cats, "PSD" = PSD_cats)



for (HMM_type in c("adjusted")) {
  dir.create(paste0("C:/UROPs/polyQ_neuronal_proteins/output/",HMM_type,"/fly/AA_charts"))
  setwd(paste0("C:/UROPs/polyQ_neuronal_proteins/output/",HMM_type,"/fly/AA_charts/"))
  
  for (location in names(location_dict)) {
    proteins <- read.csv(paste0("C:/UROPs/polyQ_neuronal_proteins/output/", HMM_type, "/fly/merged_polyAA_df.csv"))
    location_cats <- location_dict[[location]]
    proteins <- annotate_cats(proteins, location_cats)
    proteins <- proteins %>% filter(proteins$label)
    
    pdf(file = paste0("AA_charts__location-", location, "_proteins__Model-", HMM_type, "__fly.pdf"))
    for (i in 1:nrow(proteins)){
      row <- proteins[i,]
      plot_AA_Chart(row, TRUE)
    }
    dev.off()
  }
}

####################################################
# Plotting 3 -   Plot neuronally expressed proteins (based off transcriptome) that have polyAAs
####################################################

for (HMM_type in c("adjusted")) {
  dir.create(paste0("C:/UROPs/polyQ_neuronal_proteins/output/",HMM_type,"/fly/AA_charts"))
  setwd(paste0("C:/UROPs/polyQ_neuronal_proteins/output/",HMM_type,"/fly/AA_charts/"))
  
  proteins <- read.csv(paste0("C:/UROPs/polyQ_neuronal_proteins/output/", HMM_type, "/fly/merged_polyAA_df.csv"))
  neuronal_transcripts <- as.vector(read.table("C:/UROPs/polyQ_neuronal_proteins/output/fly_CNS_transcriptome_mh-l.txt", sep = "\t"))
  proteins <- proteins[proteins$ensembl_peptide_id %in% neuronal_transcripts$V1,]
  proteins <- proteins %>% filter(!is.na(proteins$IndiciesPolyQ)) # Only plot proteins that have a polyQ region
  
  
  # Remove nucleus related and transcription factor proteins
  proteins <- annotate_cats(proteins, nuclear_cats)
  proteins <- proteins %>% filter(!proteins$label)
  
  pdf(file = paste0("AA_charts__location-neuronal_transcriptome_proteins__Model-", HMM_type, "__fly.pdf"))
  for (i in 1:nrow(proteins)){
    print(i)
    row <- proteins[i,]
    plot_AA_Chart(row, TRUE)
  }
  dev.off()
}
