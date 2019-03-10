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
library(GO.db)
library(biomaRt)
library(biomartr)
library(ggplot2)
library(dplyr)
####################################################
output_base_dir <- "C:/UROPs/polyQ_neuronal_proteins/output/"

species_vec = c("fly")
candidate_AA_vec = c("Q", "T", "S", "E", "P", "G", "A", "C", "V", "M",                      
                     "I", "L", "Y", "F", "H", "K", "R", "W", "D", "N")


# TODO: Remove isoforms of proteins. They are too similiar. Could be screwing up results.
# TODO: Try playing with other categories, like maybe just active zone instead of synaptic active zone.
# TODO: write p-values of fisher t-tests to one csv file.

# Get a list of synaptic categories so we can label proteins.
# Annotate proteins based on whether they are synaptic or not.

synapse_term1 <- "GO:0045202"
synapse_cats <- get(synapse_term1, GOCCOFFSPRING)
synapse_cats <- c(synapse_cats, synapse_term1)


AZ_term1 <- "GO:0048786"
AZ_cats <- get(AZ_term1, GOCCOFFSPRING)
AZ_cats <- c(AZ_cats, AZ_term1)

PSD_term1 <- "GO:0014069"
PSD_cats <- get(PSD_term1, GOCCOFFSPRING)
PSD_cats <- c(PSD_cats, PSD_term1)


location_dict = list("AZ" = AZ_cats, "Synapse" =  synapse_cats)
location_dict = list("PSD" = PSD_cats)

# Categories to filter out of analysis. ANything that has to do with nucleus or DNA regulation.
term1 <- "GO:0005634" # Nucleus
term2 <- "GO:0009295" # nucleoid 
term3 <- "GO:0003676" # nucleic acid binding

nuclear_cats <- c(term1, get(term1, GOCCOFFSPRING), term2, get(term2, GOCCOFFSPRING))


# This is using all proteins. Repeat the fischer test after controlling for nuclear proteins.

models_dict <- list("trained" = NA, "adjusted" = NA)
models_dict <- list("adjusted" = NA)

for (model in names(models_dict)) {
    
  output_model_dir <- paste0(output_base_dir, model, "/")
  
  for (species in species_vec) {
    
    output_species_dir <- paste0(output_model_dir, species, "/")
    dir.create(output_species_dir, recursive = TRUE)
    
    pdf(file = paste0(output_species_dir, 
                      "fisher_test_results_log_plots.pdf"))
    
    for (location in names(location_dict)) {
      for (isoform_filter in c("on", "off")) {
        for (nuclear_filter in c("on", "off")) {
          print(paste0("Location: ", location))
          print(paste0("Isoform filter: ", isoform_filter))
          print(paste0("Nuclear filter: ", nuclear_filter))
          
          fisher_test_results <- rep(NA, length(candidate_AA_vec)) # For the p-value barplot
          aa_fraction_df <- as.data.frame(matrix(ncol = 2, nrow = 0)) # For the AA fraction histograms
          colnames(aa_fraction_df) <- c("TractFraction", "AA")
          
          for (i in 1:length(candidate_AA_vec)) {
            candidate_AA <- candidate_AA_vec[i]
            # debugging; candidate_AA = candidate_AA_vec[1]; candidate_AA
            
            print(paste0("Species: ", species,"; Amino acid: ", candidate_AA, " Iter: ",i))
            output_AA_dir <- paste0(output_species_dir, candidate_AA, "/")
            dir.create(output_AA_dir, recursive = TRUE)
            
            proteins <- read.csv(paste0(output_AA_dir, species, "_prots_w_HMM_", candidate_AA, ".csv"))
            # remove isoforms unless one is polyAA and one is not, in which case, keep only one of each.
            # Why do this? Fischer test assumes independence among samples for null distribution.
            # For example, if synaptic proteins have more isoforms on average than other proteins, 
            # could screw up results.
            if (isoform_filter == "on"){
              proteins <- proteins[!duplicated(proteins[c("ensembl_gene_id", "HMMhasPolyAA"),]),]
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
            
            
            # Annotate each protein as to whether it is synaptic or not.
            location_cats <- location_dict[[location]]
            proteins$label <- unlist(
              lapply(
                lapply(as.character(proteins$go_ids_CC), function(x){unlist(strsplit(x, "; "))}), 
                function(x){any(x %in% location_cats)}
              )
            )
            
            
            if (length(unique(proteins$HMMhasPolyAA)) > 1 & length(unique(proteins$label)) > 1) {
              p.val <- fisher.test(proteins$label, proteins$HMMhasPolyAA, alternative = "greater")
            } else {
              p.val = 1
            }
            
            # For the p-value histograms
            fisher_test_results[i] <- p.val
            names(fisher_test_results)[i] <- candidate_AA
            
            # For the AA fraction histograms
            polyAA_df <- proteins %>% filter(HMMhasPolyAA == TRUE)
            if (nrow(polyAA_df) == 0) {
              aa_fraction_df <- rbind(aa_fraction_df, data.frame(TractFraction = c(NA), AA = c(candidate_AA)))
            } else {
              vec <- rep(candidate_AA, length(polyAA_df$MaxPolyAARegionAAFractions))
              aa_fraction_df <- rbind(aa_fraction_df, data.frame(TractFraction = polyAA_df$MaxPolyAARegionAAFractions, AA = vec))
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
            ggtitle(paste0("Enrichment of polyAA tracts in ", location, " Proteins by Fisher Test \n",
                           "Nuclear Filter: ", nuclear_filter, ", Isoform Filter: ", isoform_filter)) +
            geom_hline(yintercept = -log(0.05), linetype="dashed", color="orange") +
            geom_hline(yintercept = -log(0.0025), linetype="dashed", color="red")
          par(mfrow = c(1,1))
          print(p)
          #############################################3
          
          # AA fraction histograms
          aa_fraction_df$TractFraction <- as.numeric(aa_fraction_df$TractFraction)
          aa_fraction_df <- aa_fraction_df %>% 
            group_by(AA) %>%
            mutate(med = median(TractFraction), mean = mean(TractFraction), group_count = n())
          aa_fraction_df$group_count[is.na(aa_fraction_df$TractFraction)] = 0
          
          # Create df for text annotations
          tmp <- aa_fraction_df[!duplicated(aa_fraction_df[,c("AA", "group_count")]),]
          ann_text <- data.frame(label = paste0("Count ", tmp$group_count), AA = tmp$AA)
          
          p <- ggplot(aa_fraction_df, aes(x = TractFraction)) +
            geom_histogram() +
            facet_wrap("AA", drop = FALSE) + 
            geom_vline(aes(xintercept = mean, group = AA), color = 'red') +
            geom_text(data = ann_text, mapping = aes(x = -Inf,  y = Inf, label = label, hjust = 0, vjust = 1)) + 
            xlab("%AA of Largest Tract") +
            ylab("Count") + 
            scale_x_continuous(lim = c(0, 1), labels = c("0", "0.25", "0.50", "0.75", "1")) +
            ggtitle(paste0("Histograms of %AA of polyAA tracts (", species, ") \n",
                           "Nuclear Filter: ", nuclear_filter, ", Isoform Filter: ", isoform_filter)) +
            theme_light()
          par(mfrow = c(1,1))
          print(p)
        }  # for (nuclear_filter in c("on", "off")) {
      } # for (isoform in c("remove", "keep")) {
    } #for (label in c("synapse", "AZ")) {
    # debugging; species = species_vec[1]
    dev.off()
  } #for (species in species_vec) {
}

# x <- read.csv("C:/UROPs/polyQ_neuronal_proteins/output/fly/Q/fly_prots_w_HMM_Q.csv", stringsAsFactors = FALSE)
# training_set_ids <- c("FBpp0289769","FBpp0307700", "FBpp0086727", "FBpp0111724", "FBpp0070830", "FBpp0309352", "FBpp0402897", "FBpp0110299","FBpp0305807")
# training_set <- x %>% filter(ensembl_peptide_id %in% training_set_ids)
# ggplot(training_set, aes(x = MaxPolyAARegionAAfractions)) +
#   geom_histogram()


