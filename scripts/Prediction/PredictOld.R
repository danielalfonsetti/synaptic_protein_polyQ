# source("https://bioconductor.org/biocLite.R")
#biocLite("DECIPHER")
source("C:/UROPs/polyQ_neuronal_proteins/scripts/ConstantsAndFunctions.R")


library(DECIPHER)
library(bamboo)
library(reticulate)
source_python("C:/UROPs/polyQ_neuronal_proteins/scripts/Prediction/disorderPredict.py")


res <- getDisorder() 
print(res)

ssTypeDict = list('H' = "Helix", 'E' = "Strand", "T" = "Turn", "C" = "Coil", "TC" = "Disordered")

# Helper functions for plotting
# Helper function 1
PolyAAChart <- function(row, graphHasHMMannots = FALSE, plotDisorder = FALSE)  {
  # Get indicies of each amino acid in the protein
  
  # row <- a[1,] # debugging
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
  if (plotDisorder) {

    for (ssType in c("H", "E", "T", "C", "TC")) {
      regionPred <- row$regionPred
      
      if (ssType == "TC") {
        ssPosList <- lapply(
          kCandidateAAs,
          function(AA) 
            unlist(lapply(
              strsplit(regionPred, ''), 
              function(regionPred) 
                which(regionPred == "T" | regionPred == "C")
            ))
        )
      } else {
        ssPosList <- lapply(
          kCandidateAAs,
          function(AA) 
            unlist(lapply(
              strsplit(regionPred, ''), 
              function(regionPred) 
                which(regionPred == ssType)
            ))
        )
      }
      # print(ssType)
      # print(ssPosList[[1]][1:10])
      
      maxLength <- max(sapply(ssPosList, length))
      
      ssPosList <- lapply(ssPosList, function(v) { c(v, rep(NA, maxLength-length(v)))})
      
      df2 <- do.call(cbind, ssPosList)
      colnames(df2) <- kCandidateAAs
      melted_df2 <- melt(df2)[,c(-1)]
      print(melted_df2[1:5, 1:2])
      tmp <- p +  geom_segment(aes(x=melted_df2$value-box_width/2, xend = melted_df2$value+box_width/2, 
                                 y=melted_df2$Var2, yend = melted_df2$Var2, color = "red"),
                             alpha = 1/box_width,
                             size = 8) +
        ggtitle(paste0("Amino acid chart for ", row$external_gene_name, " (", row$ensembl_peptide_id, ")",
                       "\nwith ", ssType, " (", ssTypeDict[ssType], ") Secondary Structure in red, polyAA in green"))
      print(tmp)
      
    } # for (ssType in c("H", "E", "T", "C")) {
  } # if (plotDisorder) {
  
  print(p)
} # PolyAAChart <- function(row, graphHasHMMannots = FALSE)  {





# a <- read.csv("C:/UROPs/polyQ_neuronal_proteins/output/adjusted/fly/Q/fly_prots_w_HMM_Q.csv")
# a <- a[1:100,]
# seq <- AAStringSet(a$peptideSeq)
# 
# regionPred <- PredictHEC(seq,
#                          type = "states",
#                          windowSize = 7,
#                          background = c(H = -0.12, E = -0.25, C = 0.23),
#                          HEC_MI1 = NULL,
#                          HEC_MI2 = NULL)
# 
# print(regionPred[1]) # Helix (H, alpha helixes), coil (C, nothing), strand (E, beta sheets)
# 
# a['regionPred'] = regionPred
# 
# PolyAAChart(a[1:3,], TRUE)



# pondr. Check 
# a <- a %>% filter(a$ensembl_peptide_id %in% c("FBpp0289769","FBpp0307700", "FBpp0086727", "FBpp0111724", "FBpp0293366",
#                   "FBpp0070830", "FBpp0309352", "FBpp0402897", "FBpp0110299","FBpp0305807"))
# seq <- AAStringSet(a$peptideSeq)
# regionPred <- PredictHEC(seq,
#                          type = "states",
#                          windowSize = 7,
#                          background = c(H = -0.12, E = -0.25, C = 0.23),
#                          HEC_MI1 = NULL,
#                          HEC_MI2 = NULL)
# 
# View(regionPred[1])
# print(regionPred[1])
# a$peptideSeq[1]
# 
# nchar(toString(a$peptideSeq[1]))
# 
# 
# unlist(lapply(
#   strsplit(regionPred[1], ''), 
#   function(regionPred) 
#     which(regionPred == "E")
# ))




####################################################################
# Bamboo secondary structure prediction method
####################################################################
# rm(list = ls())
# data(bamboo.training,
#      bamboo.validation.casp9,
#      bamboo.validation.astral30,
#      bamboo.MSA.casp9,
#      bamboo.MSA.astral30)
# 
# 
# ## Be patient, this example can take several seconds on a fast computer.
# prior.NonInfo <- bamboo.priorNonInfo()
# 
# 
# bamboo.MSA <- c(bamboo.MSA.casp9,bamboo.MSA.astral30)
# target <- "f3rvca_0"
# 
# aa <- bamboo.validation.astral30[bamboo.validation.astral30$name==target,"primary"] # get AA sequence of target
# 
# fm.NonInfo <- bamboo.estimate(likelihood(aa),prior.NonInfo,5000,500) 
# 
# fm.MSA <- bamboo.estimate(likelihood(aa), bamboo.priorMSA(bamboo.MSA[[target]]),5000,500)
# test <- bamboo.estimate(likelihood(aa), bamboo.priorNonInfo(), 5000, 500)
# test$mpState
# 
# 
# ss <- c(
#   "Truth"=bamboo.validation.astral30[bamboo.validation.astral30$name==target,"hetc"],
#   "NonInfo-MP"=fm.NonInfo$mpState,
#   "MSA-MP"=fm.MSA$mpState
# )
# plot(fm.MSA, ss)
####################################################################


likelihood <- bamboo.likelihood(bamboo.training[,"primary"],bamboo.training[,"hetc"],force=TRUE)

a <- read.csv("C:/UROPs/polyQ_neuronal_proteins/output/adjusted/fly/mergedPolyAaDf.csv")
a <- a[1:2,]
a <- a[a$ensembl_peptide_id %in% c("FBpp0289769","FBpp0307700", "FBpp0086727", "FBpp0111724", "FBpp0293366",
                  "FBpp0070830", "FBpp0309352", "FBpp0402897", "FBpp0110299","FBpp0305807"),]

regionPred <- unlist(lapply(a$peptideSeq, function(aa)  bamboo.estimate(likelihood(as.character(aa)), bamboo.priorNonInfo(), 5000, 500)$mapState))

a['regionPred'] = regionPred 

# Helper function 2
PolyAAChartWrapper <- function(proteins){
  if (nrow(proteins) == 0){
    plot(NA, xlim=c(0,2), ylim=c(0,2), bty='n',
         xaxt='n', yaxt='n', xlab='', ylab='')
    text(1, 1.7,"No data", cex = 4)
    text(1, 1.4, "No proteins in this set had polyAA regions", cex = 1.5)
    text(1, 1.2, "Cannot make amino acid plots", cex = 1.5)
  } else {
    for (i in 1:nrow(proteins)){
      row <- proteins[i,]
      PolyAAChart(row, TRUE, TRUE)
    }
  }
}

pdf("wowIMcool.pdf")
PolyAAChart(a, graphHasHMMannots = TRUE, plotDisorder = TRUE)
dev.off()

pdf("all_test.pdf")
PolyAAChartWrapper(a)
dev.off()

       