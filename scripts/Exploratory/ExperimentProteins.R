rm(list = ls())
source("C:/UROPs/polyQ_neuronal_proteins/scripts/ConstantsAndFunctions.R")


dir <- "C:/UROPs/polyQ_neuronal_proteins/output/exploratory"
dir.create(dir)
setwd(dir)

# Rbp6
hmmType <- 'adjusted'
proteins <- read.csv(paste0("C:/UROPs/polyQ_neuronal_proteins/output/", hmmType, "/fly/mergedPolyAaDf.csv"))

set <- c("FBpp0075045",
         "FBpp0305994",
         "FBpp0112339",
         "FBpp0112340",
         "FBpp0288505",
         "FBpp0301711",
         "FBpp0301712")

proteins <- proteins %>% filter(proteins$ensembl_peptide_id %in% set)
pdf(file = paste0("AA_charts_Rbb6_isoforms_Model-", hmmType, "_fly.pdf"))
PolyAAChartWrapper(proteins)
dev.off()
# END



hmmType <- 'adjusted'
proteins <- read.csv(paste0("C:/UROPs/polyQ_neuronal_proteins/output/", hmmType, "/fly/mergedPolyAaDf.csv"))


proteinsFilt <- read.csv("C:/UROPs/polyQ_neuronal_proteins/output/adjusted/fly/Q/fly_prots_w_HMM_Q_alt_spliced_nuclear_filt_transcriptome_filt.csv")
proteins <- proteins[proteins$ensembl_peptide_id %in% proteinsFilt$ensembl_peptide_id,]

setwd("C:/UROPs/polyQ_neuronal_proteins/output/exploratory")
pdf(file = paste0("AA_charts_for_fly_prots_w_HMM_Q_alt_spliced_nuclear_filt_transcriptome_filt.csv-", hmmType, "_fly.pdf"))
PolyAAChartWrapper(proteins)
dev.off()
