####################################################
# Daniel Alfonsetti, daniel.alfonsetti@gmail.com
# MIT, Littleton Lab UROP
# 17 April 2019
# ---------------------------------
# Description: Make plots summarising results from HMM
####################################################


#########################################
# Load Libraries
#########################################
rm(list = ls())
source("../../ConstantsAndFunctions.R", chdir=T)
kOutputBaseDir <- "../../../output/"

library(RColorBrewer)
library(stringr)
library(reshape2)
library(grid)

#########################################
# Command Line Args (if using command line)
#########################################
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 0) {
  kModels =  unlist(strsplit(args[1],","))
  kSpecies = unlist(strsplit(args[2],","))
  kCandidateAAs = strsplit(args[3],"")[[1]]
}
#########################################
# HELPER FUNCTIONS
#########################################

# Plotting helper function
# Number of randomly chosen 50 amino acid peptide chunks
RandomSubstr <- function(x, length) {
  if (nchar(x) < length) {return(x)}
  else {
    start = sample(1:(nchar(x)-length), 1, replace=T)
    end = start + (length-1)
    random_substring = substr(x, start, end)
    return(random_substring)
  }
} # RandomSubstr <- function(x, length)

#########################################
# MAIN (HMM summary plots)
#########################################
for (model in kModels) {
  for (species in kSpecies) {
    for (candidateAA in kCandidateAAs) {
      
      outputDir <- paste0(kOutputBaseDir, model, "/", species, "/", candidateAA, "/")
      dir.create(outputDir, recursive = TRUE)
      proteins <- read.csv(paste0(outputDir, species, "_prots_w_HMM_", candidateAA, "_nuclear_filt_transcriptome_filt.csv"))
      if (nrow(proteins) == 0) {
        pdf(paste0(outputDir, "/sequence_summary_plots.pdf"))
        dev.off()
        next
      }
      
      pdf(paste0(outputDir, "/sequence_summary_plots.pdf"))
      
      colors1 <- brewer.pal(8, "Dark2")[1]
      colors2 <- brewer.pal(8, "Dark2")[round(seq(from = 1, to = 8,length.out = 2))]
      colors3 <- brewer.pal(8, "Dark2")[round(seq(from = 1, to = 8,length.out = 3))]
      
      ####################
      # Plot: Whole Protein AAfraction Histogram for all proteins
      p <- ggplot(data = proteins, aes(x = proteins$aaComp)) +
        geom_histogram(alpha=0.5, color = colors1, 
                       fill = colors1) +
        geom_vline(aes(xintercept=mean(aaComp)),
                   color=colors1, linetype="dashed", size=1) + 
        xlab(paste0("%", candidateAA, " (Whole Protein Seqs)")) +
        ylab("# Proteins") +
        ggtitle(paste0("%", candidateAA, " Whole Protein")) +
        theme_light()
      print(p)
      
      filteredDf <- proteins %>% filter(proteins$aaComp> quantile(proteins$aaComp, 0.95))
      filteredDf <- filteredDf %>% arrange(desc(filteredDf$aaComp))
      write.csv(filteredDf, file = paste0(outputDir, species, "_", candidateAA, "fraction_prots_95sig.csv"))
      
      ####################
      # Plot: Whole Protein AAfraction Histogram split by whether the protein has at least one polyAA tract or not.
      mu <- proteins %>% group_by(hmmHasPolyAA) %>%
        summarise(hmmHasPolyAA.mean = mean(aaComp))
      p <- ggplot(data = proteins, aes(x = aaComp, fill = hmmHasPolyAA, color = hmmHasPolyAA)) +
        geom_histogram(alpha=0.5, position="identity") +
        geom_vline(data=mu, aes(xintercept= hmmHasPolyAA.mean, color = hmmHasPolyAA),
                   linetype='dashed', size=1, show.legend = FALSE) + 
        xlab(paste0("%", candidateAA, " (Whole Protein Seqs)")) +
        ylab("# Proteins") +
        labs(fill = paste0("Contains \n poly", candidateAA)) +
        ggtitle(paste0("%", candidateAA," Whole Proteins Histogram")) +
        theme_light()+
        scale_colour_manual(values = colors2)  +
        scale_fill_manual(values = colors2)  +
        guides(color = FALSE)
      print(p)
      
      # Same as previous plot, by now in density form.
      p <- ggplot(data = proteins, aes(x = aaComp, fill = hmmHasPolyAA)) +
        geom_density(alpha = 0.5) +
        xlab(paste0("%", candidateAA, " (Whole Protein Seqs)")) +
        ylab("Density") +
        labs(fill = paste0("Contains \n poly", candidateAA)) +
        ggtitle(paste0("%", candidateAA, " Whole Proteins Density Plot")) +
        scale_colour_manual(values = colors2)  +
        scale_fill_manual(values = colors2)  +
        theme_light()
      print(p)
      
      ####################
      # Plot: Histogram of 100000 50AA chunks based on %AA, each randomly chosen from randomly chosen Protein seqs
      df <- sample_n(proteins, 100000, replace = TRUE)
      df$RandomSeqs <-unlist(lapply(as.character(df$peptideSeq), RandomSubstr, 50))
      nAAs = str_count(df$RandomSeqs, candidateAA)
      df$aaComp <- nAAs/nchar(df$RandomSeqs)
      
      p <- ggplot(data = df, aes(x = df$aaComp)) +
        geom_histogram(show.legend = FALSE, alpha = 0.5,                 
                       fill = colors1, 
                       color = colors1) +
        geom_vline(aes(xintercept=mean(aaComp)), linetype="dashed", size=1, color = colors1) +
        xlab(paste0("%", candidateAA, " (Random 50AA Seqs)")) +
        ylab("# of Randomly Choosen Pepitides") +
        ggtitle(paste0("%", candidateAA, " for 100000 50AA Random Seqs Histogram")) +
        theme_light()
      print(p)
      
      #####################
      # Get proteins with the most significant amount of AA.
      filteredDf <- df %>% filter(df$aaComp > quantile(df$aaComp, 0.95))
      filteredDf <- filteredDf %>% arrange(desc(filteredDf$aaComp))
      write.csv(filteredDf, file = paste0(outputDir, species, "_", candidateAA, "fraction_randomSeqs_95sig.csv"))
      
      #####################
      # If there are no proteins with polyAA regions, the next graphs would have nothing to compare, so skip.
      if (sum(proteins$hmmHasPolyAA) == 0) {
        dev.off()
        next()
      }
      
      # Plot number of Proteins annotated by HMM to be "Poly AA" vs avg %AA of all polyAA denoted regions.
      df <- proteins %>% filter(hmmHasPolyAA)
      p <- ggplot(data = df, aes(x = as.numeric(df$avgPolyAARegionAAFractions))) +
        geom_histogram(show.legend = FALSE, alpha = 0.5,                 
                       fill = colors1, 
                       color = colors1) +
        geom_vline(aes(xintercept=mean(avgPolyAARegionAAFractions, na.rm = TRUE)), linetype="dashed", size=1, color = colors1) +        
        xlab(paste0("Avg %", candidateAA, " of poly", candidateAA, " Regions")) + 
        ggtitle(paste0("Avg %", candidateAA, " of poly", candidateAA," Regions for poly", candidateAA, " Containing Proteins")) +
        ylab(paste0("# Proteins Containing poly", candidateAA)) +
        theme_light()
      print(p)
      
      # Plot number of Proteins annotated by HMM to be "Poly AA" vs max %AA among all polyAA denoted regions.
      p <- ggplot(data = df, aes(x = df$maxPolyAARegionAAFractions)) +
        geom_histogram(show.legend = FALSE, alpha = 0.5,                 
                       fill = colors1, 
                       color = colors1) +
        geom_vline(aes(xintercept=mean(df$maxPolyAARegionAAFractions, na.rm = TRUE)), linetype="dashed", size=1, color = colors1) +              
        xlab(paste0("Max %", candidateAA, " of poly", candidateAA, " Regions")) +
        ggtitle(paste0("Max %", candidateAA, " of poly", candidateAA, " Regions for poly", candidateAA, " Containing Proteins")) +
        ylab(paste0("# Proteins Containing poly", candidateAA)) +
        theme_light()
      print(p)
      
      #####################
      # Now compare polyAA regions to a set of randomly chosen sequences whose sequence lengths
      # follow the same distribution as the distribution of lengths of the polyAA regions.
      polyAADf <- proteins %>% filter(hmmHasPolyAA)
  
      random_seq_df <- data.frame(matrix(nrow = nrow(polyAADf), ncol = 1))
      colnames(random_seq_df) <- c("aaComp")
      for (i in 1:nrow(polyAADf)) {
        target_length <- polyAADf[i,]$maxLengthPolyAA
        flag = TRUE
        while (flag){
          try = as.character(sample_n(proteins, 1)$peptideSeq)
          if (nchar(try) > target_length) {
            res <- RandomSubstr(try, target_length)

            nAAs = str_count(res, candidateAA)
            random_seq_df$aaComp[i] <- nAAs/nchar(res)
            flag = FALSE
          } # if (nchar(try) > target_length) {
        } # while (flag){
      } # for (i in 1:nrow(polyAADf)) {
      
      
      # Plot: Compare %AA of randomly chosen tracts of similiar size to max %AA among polyAA regions of a given protein
      df <- data.frame(randomAaComp = random_seq_df$aaComp, maxAaComp= polyAADf$maxPolyAARegionAAFractions)
      df <- melt(df)
      colnames(df) <- c("seqType", "aaComp")
      mu <- df %>% group_by(df$seqType) %>% summarise(aaComp.mean = mean(aaComp))
      colnames(mu) <- c("seqType", "aaComp.mean")
      
      p <- ggplot() +
        geom_histogram(data = df, aes(x = aaComp, fill = seqType, color = seqType), alpha = 0.5, position="identity") +
        geom_vline(data=mu, aes(xintercept= aaComp.mean, color = seqType),
                   linetype='dashed', size=1, show.legend = FALSE) + 
        xlab(paste0("%", candidateAA)) +
        ylab(paste0("# Proteins")) +
        ggtitle(paste0("Max %", candidateAA, " of poly", candidateAA, " Regions for poly", candidateAA, " Containing Proteins \n vs %Q of randomly chosen tracts of similiar size")) +
        scale_fill_manual(values = colors2, labels = c("Random Seqs", "Max %Q Seqs"))  +
        scale_colour_manual(values = colors2)  +
        theme_light() + 
        labs(fill = "Sequence Type") + 
        guides(color = FALSE)
      print(p)
      
      # Plot: Compare %AA of randomly chosen tracts of similiar size to avg %AA of polyAA regions of a given protein
      df <- data.frame(randomAaComp= random_seq_df$aaComp, avgAaComp= polyAADf$avgPolyAARegionAAFractions)
      df <- melt(df)
      colnames(df) <- c("seqType", "aaComp")
      mu <- df %>% group_by(df$seqType) %>% summarise(aaComp.mean = mean(aaComp))
      colnames(mu) <- c("seqType", "aaComp.mean")
      
      
      tryCatch({
        test <- t.test(random_seq_df$aaComp, polyAADf$avgPolyAARegionAAFractions, alternative = "two.sided", var.equal = FALSE)
        pval <- round(test$p.value, 3)
      }, error = function(e) {
        pval <- 1
      })
      annot_m <- grobTree(textGrob(paste0("t-test pval: ", pval), x= unit(0.5, "npc"), y = unit(0.97, "npc"),
                                   gp=gpar(col='red', fontsize=15, fontface="italic")))

      p <- ggplot() +
        geom_histogram(data = df, aes(x = aaComp, fill = seqType, color = seqType), alpha = 0.5, position="identity") +
        geom_vline(data=mu, aes(xintercept= aaComp.mean, color = seqType),
                   linetype='dashed', size=1, show.legend = FALSE) + 
        xlab(paste0("%", candidateAA)) +
        ylab(paste0("# Proteins")) +
        ggtitle(paste0("Avg %", candidateAA, " of Poly", candidateAA, " Regions for poly", candidateAA, " Containing Proteins \n vs %", candidateAA, " of Randomly Chosen Tracts")) +
        scale_fill_manual(values = colors2, labels = c("Random Seqs", "Avg %Q Seqs"))  +
        scale_colour_manual(values = colors2)  +
        theme_light() +
        labs(fill = "Sequence Type") +
        guides(color = FALSE) + 
        annotation_custom(annot_m)
      print(p)
      
      # Plot: All three on the same chart
      # Max %AA vs Avg %AA vs %AA of Random Tracts
      df <- data.frame(randomAaComp = random_seq_df$aaComp, 
                       avgAaComp = polyAADf$avgPolyAARegionAAFractions,
                       maxAaComp = polyAADf$maxPolyAARegionAAFractions)
      df <- melt(df)
      colnames(df) <- c("seqType", "aaComp")
      mu <- df %>% group_by(df$seqType) %>% summarise(aaComp.mean = mean(aaComp))
      colnames(mu) <- c("seqType", "aaComp.mean")
      p <- ggplot() +
        geom_histogram(data = df, aes(x = aaComp, fill = seqType, color = seqType), alpha = 0.5, position="identity") +
        geom_vline(data=mu, aes(xintercept= aaComp.mean, color = seqType),
                   linetype='dashed', size=1, show.legend = FALSE) + 
        xlab(paste0("%", candidateAA)) +
        ylab(paste0("# Proteins")) +
        ggtitle(paste0("Avg %", candidateAA," of poly", candidateAA, " Regions \n vs Max %", candidateAA, " of poly", candidateAA," Regions \n vs %Q of randomly chosen tracts")) +
        scale_fill_manual(values = colors3, label = c("Random Seqs", "Avg %Q Seqs", "Max %Q Seqs"))  +
        scale_colour_manual(values = colors3)   +
        theme_light() +
        labs(fill = "Sequence Type") +
        guides(color = FALSE)
      print(p)
      
      #####################
      # Plot: Number of Proteins annotated by HMM to be PolyAA vs average length of Region
      df <- proteins %>% filter(hmmHasPolyAA)
      p <- ggplot(data = df, aes(x = df$avgLengthPolyAA)) +
        geom_histogram(show.legend = FALSE, alpha = 0.5,                 
                       fill = brewer.pal(8, "Dark2")[1], 
                       color = brewer.pal(8, "Dark2")[1]) +
        geom_vline(aes(xintercept=mean(avgLengthPolyAA, na.rm = TRUE)), linetype="dashed", size=1, color = colors1) +   
        xlab(paste0("Average Poly", candidateAA, " Length")) +
        ylab(paste0("# Proteins Containing poly", candidateAA, " (by HMM)")) +
        ggtitle(paste0("Average Poly", candidateAA, " Region Length in poly", candidateAA, " Containing Proteins")) +
        theme_light()
      print(p)
      
      # Plot: Number of Proteins annotated by HMM to be PolyAA vs Max length of Region
      p <- ggplot(data = df, aes(x = df$maxLengthPolyAA)) +
        geom_histogram(show.legend = FALSE, alpha = 0.5,                 
                       fill = brewer.pal(8, "Dark2")[1], 
                       color = brewer.pal(8, "Dark2")[1]) +
        geom_vline(aes(xintercept=mean(maxLengthPolyAA, na.rm = TRUE)), linetype="dashed", size=1, color = colors1)  +
        xlab(paste0("Max poly", candidateAA, " Length")) +
        ylab(paste0("# Proteins Containing poly", candidateAA, " (by HMM)")) +
        ggtitle(paste0("Max poly", candidateAA, " Region Length in poly", candidateAA, " Containing Proteins")) +
        theme_light()
      print(p)
      
      # Plot: Number proteins by number of polyAA regions.
      p <- ggplot(data = df, aes(x = df$numPolyAA)) + 
        geom_bar(show.legend = FALSE, alpha = 0.5,                 
                 fill = brewer.pal(8, "Dark2")[1], 
                 color = brewer.pal(8, "Dark2")[1]) +
        geom_vline(aes(xintercept=mean(numPolyAA, na.rm = TRUE)), linetype="dashed", size=1, color = colors1)  +
        xlab(paste0("Number of poly", candidateAA, " Regions")) +
        ylab(paste0("# of ", candidateAA, " Proteins")) +
        ggtitle(paste0("Number of poly", candidateAA, " Regions in poly", candidateAA, " Containing Proteins")) +
        theme_light()
      print(p)
      
      dev.off()
    } # for (model in c("adjusted", "trained")) {
  } # for (species in kSpecies) {
} # for (candidateAA in kCandidateAAs) {
# END
