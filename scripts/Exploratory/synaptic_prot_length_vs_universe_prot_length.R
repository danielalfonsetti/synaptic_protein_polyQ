
# Libraries
library(ggplot2)
library(biomaRt)
library(biomartr)
library(STRINGdb)
library(org.Dm.eg.db)
library(org.Hs.eg.db) 
library(org.Mm.eg.db)
library(org.Ce.eg.db)
library(GO.db)


proteins <- read.csv("C:/UROPs/polyQ_neuronal_proteins/output/fly/fly_prots.csv", stringsAsFactors = FALSE)

protein_lengths <- unlist(lapply(proteins$peptide_seq, function(x){length(unlist(strsplit(x, "")))}))

# ggplot() + aes(protein_lengths)+ geom_histogram(binwidth=1, colour="black", fill="white")
# Length of universe proteins

p <- ggplot() + 
  aes(protein_lengths)+
  geom_histogram(show.legend = FALSE, aes(fill=..count..), binwidth = 50) +
  xlim(c(0, 7500)) +
  ylim(c(0, 2100)) +
  xlab("# of Amino Acids") +
  ylab("# of Proteins") +
  ggtitle(paste0("Fly Protein Length Histogram (All Proteins)")) +
  theme_light()
print(p)

# Length of synapse proteins 
synapse_term1 <- "GO:0045202"
synapse_cats <- get(synapse_term1, GOCCOFFSPRING)
synapse_cats <- c(synapse_cats, synapse_term1)

proteins_synapse_only <- proteins[unlist(
  lapply(
    lapply(proteins$go_ids_CC, function(x){unlist(strsplit(x, "; "))}), 
    function(x){any(x %in% synapse_cats)}
  )
)
,]
nrow(proteins); nrow(proteins_synapse_only)

synapse_protein_lengths <- unlist(lapply(proteins_synapse_only$peptide_seq, function(x){length(unlist(strsplit(x, "")))}))

p <- ggplot() + 
  aes(synapse_protein_lengths)+
  geom_histogram(show.legend = FALSE, aes(fill=..count..), binwidth = 50) +
  xlim(c(0, 7500)) +
  ylim(c(0, 2100)) +
  xlab("# of Amino Acids") +
  ylab("# of Proteins") +
  ggtitle(paste0("Fly Protein Length Histogram (Synaptic Proteins)")) +
  theme_light()
print(p)


# There are some outliers....
summary(synapse_protein_lengths)
summary(protein_lengths)
t_test_res <- t.test(synapse_protein_lengths, protein_lengths)
t_test_res$p.value

# Try capping to deal with outliers... 
# http://r-statistics.co/Outlier-Treatment-With-R.html
cap <- function(vec) {
  x <- vec
  qnt <- quantile(x, probs=c(.25, .75), na.rm = T)
  caps <- quantile(x, probs=c(.05, .95), na.rm = T)
  H <- 1.5 * IQR(x, na.rm = T)
  x[x < (qnt[1] - H)] <- caps[1]
  x[x > (qnt[2] + H)] <- caps[2]
  return(x)
}

synapse_protein_lengths_capped <- cap(synapse_protein_lengths)
protein_lengths_capped <- cap(protein_lengths)

summary(synapse_protein_lengths_capped)
summary(protein_lengths_capped)
t_test_res <- t.test(synapse_protein_lengths_capped, protein_lengths_capped)
t_test_res$p.value
# Still significantly different. Seems that synaptic proteins are longer on average.
# Therefore we will have to do adjusted-p value of category vs avergage peptide length of category