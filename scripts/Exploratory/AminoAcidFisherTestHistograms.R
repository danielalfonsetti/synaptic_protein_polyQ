####################################################
# Daniel Alfonsetti, daniel.alfonsetti@gmail.com
# MIT, Littleton Lab UROP
# 24 January 2019
# ---------------------------------
# Description: For each amino acid, see if synaptic proteins are enriched for polyAA tracts.
# If our hypothesis is correct, synaptic proteins will only be enriched for polyQ tracts.
# We use a one sided fisher exact test for each amino acid. 
# (2x2 matrix with rows and columsn of synaptic/not-synpatic and polyAA/not-polyAA, where polyAA for a protein is determined
# by HMM output, and synaptic/not-synaptic is determined by GO categories).

# After this, as a diagnostic, we analyze the distribution of the amino acid 
# fraction of the polyAA regions of polyAA-denoted proteins
# for each of the 20 amino acids.
####################################################
rm(list = ls())
source("C:/UROPs/polyQ_neuronal_proteins/scripts/ConstantsAndFunctions.R")
# TODO: Fix isoform filter.

library(GO.db)
library(biomaRt)
library(biomartr)
library(ggplot2)
library(dplyr)
source("C:/UROPs/polyQ_neuronal_proteins/scripts/ConstantsAndFunctions.R")
####################################################

# Get lists of various neuron categories so we can label proteins.
# We will annotate proteins based on whether they are in the category or not, and see if this relation is 
# associated with whether they have a polyAA track.

####################################################
# Helper functions
####################################################
FisherTestComparison <- function(locationCats = NA, 
                                   isoform_filter = "on", nuclear_filter = "on", neuron_filter = "on", 
                                   proteinsNoHmm = NA, output_species_dir = NA) {
  

  candidate_AA_vec = c("A", "V", "L", "I", "M", "F","W", "C", "Y", "S", "T", "G", "P", "Q","N", "E", "D", "R", "H", "K")
  
  fisher_test_results <- rep(NA, length(candidate_AA_vec)) # For the p-value barplot
  aa_fraction_df <- as.data.frame(matrix(ncol = 3, nrow = 0)) # For the AA fraction histograms
  colnames(aa_fraction_df) <- c("TractFraction", "AA", "label")
  
  

  for (i in 1:length(candidate_AA_vec)) {
    candidate_AA <- candidate_AA_vec[i]
    # debugging; candidate_AA = candidate_AA_vec[1]; candidate_AA
    
    print(paste0("Species: ", species,"; Amino acid: ", candidate_AA, " Iter: ",i))
    output_AA_dir <- paste0(output_species_dir, candidate_AA, "/")
    dir.create(output_AA_dir, recursive = TRUE)
    
    proteins <- read.csv(paste0(output_AA_dir, species, "_prots_w_HMM_", candidate_AA, ".csv"))
    
    ################### 
    # Filtering Steps
    ################### 
    
    # remove isoforms unless one is polyAA and one is not, in which case, keep only one of each.
    # Why do this? Fischer test assumes independence among samples for null distribution.
    # For example, if synaptic proteins have more isoforms on average than other proteins, 
    # could screw up results.
    if (isoform_filter == "on"){
      proteins <- proteins[!duplicated(proteins[,c("ensembl_gene_id", "hmmHasPolyAA")]),]
    }
    
    if (nuclear_filter == "on") {
      proteins$nuclear <- unlist(
        lapply(
          lapply(as.character(proteins$go_ids_CC), function(x){unlist(strsplit(x, "; "))}),
          function(x){any(x %in% nuclear_cats)}
        )
      )
      proteins = proteins %>% filter(proteins$nuclear == FALSE)
    }
    
    # kNeuronalTranscripts <- as.vector(read.table("C:/UROPs/polyQ_neuronal_proteins/output/fly_CNS_transcriptome_mh-l.txt", sep = "\t"))
    if (neuron_filter == "on") {
      proteins <- proteins[proteins$ensembl_peptide_id %in% kNeuronalTranscripts$V1,]
    }
    
    ################### 
    ################### 
    # Labeling and p-value calculation steps.
    
    
    # Annotate each protein as to whether it is in the GO category list or not.
    proteins <- AnnotateByCategorySet(df = proteins, cats = locationCats)

    if (length(unique(proteins$hmmHasPolyAA)) > 1 & length(unique(proteins$inSet)) > 1) {
      p.val <- fisher.test(proteins$inSet, proteins$hmmHasPolyAA, alternative = "greater")
    } else {
      p.val = 1
    }
    
    # For the p-value histograms
    fisher_test_results[i] <- p.val
    names(fisher_test_results)[i] <- candidate_AA
    
    # For the AA fraction histograms
    polyAA_df <- proteins %>% filter(hmmHasPolyAA == TRUE)
    if (nrow(polyAA_df) == 0) {
      aa_fraction_df <- rbind(aa_fraction_df, data.frame(TractFraction = c(NA), AA = c(candidate_AA), loc_label = c(NA)))
    } else {
      vec <- rep(candidate_AA, length(polyAA_df$MaxPolyAARegionAAFractions))
      aa_fraction_df <- rbind(aa_fraction_df, data.frame(TractFraction = polyAA_df$MaxPolyAARegionAAFractions, AA = vec , loc_label = polyAA_df$inSet))
    }
  } # for (i in 1:length(candidate_AA_vec)) {
  
  # p-value histogram
  p_value_hist_df = data.frame(names = factor(names(fisher_test_results), levels = names(fisher_test_results)),
                               pvalues = unlist(unname(fisher_test_results)))
  p <- ggplot(p_value_hist_df, aes(x = names, y = -log(pvalues), fill = pvalues)) + 
    geom_col(alpha = .5, show.legend = FALSE) + 
    theme_light() +
    xlab("Amino Acid") +
    ylab("-log(Fisher Test P-value)") + 
    ggtitle(paste0("Association between polyAA tract proteins and ", location, " proteins \n",
                   "Nuclear Filter: ", nuclear_filter, ", Isoform Filter: ", isoform_filter, ", Neuronal Filter: ", neuron_filter)) +
    geom_hline(yintercept = -log(0.05), linetype="dashed", color="orange") +
    geom_hline(yintercept = -log(0.0025), linetype="dashed", color="red")
  par(mfrow = c(1,1))
  print(p)
  #############################################
  # AA fraction histograms
  
  aa_fraction_df$TractFraction <- as.numeric(aa_fraction_df$TractFraction)
  aa_fraction_df <- aa_fraction_df %>% 
    group_by(AA) %>%
    mutate(med = median(TractFraction), mean = mean(TractFraction), group_count = n(), group_count_cat = sum(loc_label))
  aa_fraction_df$group_count[is.na(aa_fraction_df$TractFraction)] = 0
  
  # Create dfs for text annotations
  tmp <- aa_fraction_df[!duplicated(aa_fraction_df[,c("AA", "group_count", "group_count_cat")]),]
  ann_text <- data.frame(label = paste0("Count: ", tmp$group_count), AA = tmp$AA)
  ann_text_2 <- data.frame(label = paste0("Cat Count: ", tmp$group_count_cat), AA = tmp$AA)
  
  
  p <- ggplot(aa_fraction_df, aes(x = TractFraction)) +
    geom_histogram() +
    facet_wrap("AA", drop = FALSE) + 
    geom_vline(aes(xintercept = mean, group = AA), color = 'red') +
    geom_text(data = ann_text, mapping = aes(x = -Inf,  y = Inf, label = label, hjust = 0, vjust = 1)) + 
    geom_text(data = ann_text_2, mapping = aes(x = -Inf,  y = Inf, label = label, hjust = 0, vjust = 2)) + 
    xlab("%AA of Largest Tract") +
    ylab("Count") + 
    scale_x_continuous(lim = c(0, 1), labels = c("0", "0.25", "0.50", "0.75", "1")) +
    theme_light()
  
  if (neuron_filter == "on") {
    p <- p + ggtitle(paste0("Histograms of %AA of polyAA tracts (", species, ") \n",
                   "# of prots in neuronal transcriptome: ", nrow(kNeuronalTranscripts), "\n",
                   "# of prots in category (", location, "): ", sum(proteinsNoHmm$inSet), "\n",
                   "Nuclear Filter: ", nuclear_filter, ", Isoform Filter: ", isoform_filter, ", Neuronal Filter: ", neuron_filter))
  } else {
    p <- p + ggtitle(paste0("Histograms of %AA of polyAA tracts (", species, ") \n",
                            "# of prots in genome: ", nrow(proteinsNoHmm), "\n",
                            "# of prots in category (", location, "): ", sum(proteinsNoHmm$inSet), "\n",
                            "Nuclear Filter: ", nuclear_filter, ", Isoform Filter: ", isoform_filter, ", Neuronal Filter: ", neuron_filter))
  }
  
  par(mfrow = c(1,1))
  print(p)
}
####################################################
# End of helper functions
####################################################

# Exploratory analysis
kLocations = list("AZ" = kAZcats, "Synapse" =  kSynapseCats, "PSD" = kPSDcats, 
                  "Transcription_factor (+ cntrl)"= kTranscriptionCats, 
                  "AZ+PSD" = c(kAZcats, kPSDcats))

for (model in kModels) {
  for (species in kSpecies) {
    for (location in names(kLocations)) {
      
      fisherOutputDir <- paste0(kOutputBaseDir, model, "/", species, "/", "AA_Fisher_Tests", "/")
      dir.create(fisherOutputDir, recursive = TRUE)
      pdf(file = paste0(fisherOutputDir, species, "_AA_Fisher_Tests__Model-", model, 
                        "__Location-", location, ".pdf"))
      
      proteinsNoHmm <- read.csv(paste0(kOutputBaseDir, species, "_prots.csv"))
      locationCats <- kLocations[[location]]
      proteinsNoHmm <- AnnotateByCategorySet(df = proteinsNoHmm, cats = locationCats)

      for (isoformFilter in c("on", "off")) {
        for (nuclearFilter in c("on", "off")) {
          for (neuronFilter in c("on", "off")) {
            print(paste0("Location: ", location))
            print(paste0("Isoform filter: ", isoform_filter))
            print(paste0("Nuclear filter: ", nuclear_filter))
            
            FisherTestComparison(locationCats, isoformFilter, nuclearFilter, neuronFilter, proteinsNoHmm, output_species_dir)
        
          } # for (neuronal_filter in c("on", "off")) {
        }  # for (nuclear_filter in c("on", "off")) {
      } # for (isoform in c("on", "off")) {
      
      dev.off()
    } #(location in names(kLocations)) {
  } #for (species in kSpecies) {
} # for (model in kModels) {
# END


