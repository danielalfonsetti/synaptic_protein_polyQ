# source("https://bioconductor.org/biocLite.R")
#biocLite("DECIPHER")
source("../../scripts/ConstantsAndFunctions.R", chdir=T)


library(DECIPHER)
library(bamboo)
library(reticulate)


# Helper functions for plotting
# Helper function 1
PolyAAChart <- function(row, graphHasHMMannots = FALSE, plotDisorder = FALSE)  {

  
  # Get indicies of each amino acid in the protein
  # row <- b[1,] # debugging
  peptideSeq <- as.character(row$peptideSeq)
  aaPosList <- lapply(
    kCandidateAAs,
    function(AA) 
      unlist(lapply(
        strsplit(peptideSeq, ''), 
        function(peptideSeq) 
          which(peptideSeq == AA)
      ))
  )
  
  ########## Convert list to dataframe ############3
  ## Compute maximum length
  maxLength <- max(sapply(aaPosList, length))
  
  ## Add NA values to list elements
  aaPosList <- lapply(aaPosList, function(v) { c(v, rep(NA, maxLength-length(v)))})
  
  ## Cbind
  df <- do.call(cbind, aaPosList)
  colnames(df) <- kCandidateAAs
  
  # Reorganize into a form useful for plotting
  melted_df <- melt(df)[,c(-1)]
  
  box_width = 7
  p <- ggplot() +
    geom_segment(aes(x=melted_df$value-box_width/2, xend = melted_df$value+box_width/2, 
                     y=melted_df$Var2, yend = melted_df$Var2),
                 alpha = 1/box_width,
                 size = 8) +
    xlab("Peptide Position Index") +
    ylab("Amino Acid") +
    ggtitle(paste0("Amino acid chart for ", row$external_gene_name, " (", row$ensembl_peptide_id, ")")) +
    theme_light() +
    theme(legend.position = "none")
  
  ############################
  if (graphHasHMMannots) {
    
    hmmAnnotationsList <- lapply(
      kCandidateAAs,
      function(AA) 
        unlist(lapply(strsplit(as.character(eval(parse(text = paste0("row$indiciesPoly", AA)))), split = "; ")[[1]], function(i) eval(parse(text = i))))
    )
    print("Look here")
    print(hmmAnnotationsList)
    
    # Remove first and last index in each list (because we plot each point as a line from its middle position to -1 and +1).
    hmmAnnotationsList <- lapply(hmmAnnotationsList, function(vec){vec[2:(length(vec)-1)]})
    
    # Find the vector with the max length, and pad the rest of the vectors so that we can convert the list to the dataframe.
    maxLength <- max(sapply(hmmAnnotationsList, length))
    hmmAnnotationsList <- lapply(hmmAnnotationsList, function(v) {if (is.na(v[1])) {rep(NA, maxLength)} else {c(v, rep(NA, maxLength-length(v)))}})
    
    df <- do.call(cbind, hmmAnnotationsList)
    colnames(df) <- kCandidateAAs
    meltedHmmAnnotsDf <- melt(df)[,c(-1)]
    if (!all(is.na(meltedHmmAnnotsDf$value))) {
      p <- p + geom_segment(aes(x=meltedHmmAnnotsDf$value-1, xend = meltedHmmAnnotsDf$value+1, 
                                y=meltedHmmAnnotsDf$Var2, yend = meltedHmmAnnotsDf$Var2),
                            alpha = 1, color = "Green", size = 2.5) +
        ggtitle(paste0("Amino acid chart for ", row$external_gene_name, " (", row$ensembl_peptide_id, ")", 
                       "\nwith polyAA in green."))
        
    }
  } # if (graphHasHMMannots) {
  
  ##############
  # Intrinsically disordered region predictor plotting (using DICPHER method)
  if (plotDisorder & (!all(unlist(row$disorderedPred) %in% c("", ":")))) {

    
    disorderAnnotationsList <- lapply(
      kCandidateAAs,
      function(AA) 
        unlist(lapply(strsplit(as.character(eval(parse(text = "row$disorderedPred"))), split = "; ")[[1]], function(i) eval(parse(text = i))))
    )
    print("Look here")
    print(disorderAnnotationsList)
    
    # Remove first and last index in each list (because we plot each point as a line from its middle position to -1 and +1).
    disorderAnnotationsList <- lapply(disorderAnnotationsList, function(vec){vec[2:(length(vec)-1)]})
    
    # Find the vector with the max length, and pad the rest of the vectors so that we can convert the list to the dataframe.
    maxLength <- max(sapply(disorderAnnotationsList, length))
    disorderAnnotationsList <- lapply(disorderAnnotationsList, function(v) {if (is.na(v[1])) {rep(NA, maxLength)} else {c(v, rep(NA, maxLength-length(v)))}})
    
    df <- do.call(cbind, disorderAnnotationsList)
    colnames(df) <- kCandidateAAs
    meltedDisorderAnnotsDf <- melt(df)[,c(-1)]
    if (!all(is.na(meltedDisorderAnnotsDf$value))) {
      p <- p + geom_segment(aes(x=meltedDisorderAnnotsDf$value-1, xend = meltedDisorderAnnotsDf$value+1, 
                                y=meltedDisorderAnnotsDf$Var2, yend = meltedDisorderAnnotsDf$Var2),
                            alpha = 1, color = "Red", size = 2.5) +
        ggtitle(paste0("Amino acid chart for ", row$external_gene_name, " (", row$ensembl_peptide_id, ")", 
                       "\nwith polyAA in green and disordered in red."))
    }
  } # if (plotDisorder) {
  
  print(p)
} # PolyAAChart <- function(row, graphHasHMMannots = FALSE)  {

# Helper function 2
PolyAAChartWrapper <- function(proteins, graphHasHMMannots = TRUE, plotDisorder = TRUE){
  if (nrow(proteins) == 0){
    plot(NA, xlim=c(0,2), ylim=c(0,2), bty='n',
         xaxt='n', yaxt='n', xlab='', ylab='')
    text(1, 1.7,"No data", cex = 4)
    text(1, 1.4, "No proteins in this set had polyAA regions", cex = 1.5)
    text(1, 1.2, "Cannot make amino acid plots", cex = 1.5)
  } else {
    for (i in 1:nrow(proteins)){
      print(i)
      row <- proteins[i,]
      PolyAAChart(row, graphHasHMMannots, plotDisorder)
    }
  }
}


a <- read.csv("../../output/adjusted/fly/mergedPolyAaDfDisorder.csv")

# b <- a %>% filter(!is.na(a$uniprotsptrembl))
# b <- b[1:30,]
# b$disorderedPred <-lapply(b$uniprotsptrembl, function(proteinName) paste0(lapply(getDisorder(proteinName), function(interval) paste0(interval[1], ":", interval[2])), collapse="; "))

pdf("disordered_test.pdf")
PolyAAChart(b) 
dev.off()

pdf("all_test_disorder.pdf")
PolyAAChartWrapper(b, graphHasHMMannots = TRUE, plotDisorder = TRUE)
dev.off()

       