####################################################
# Daniel Alfonsetti, daniel.alfonsetti@gmail.com
# MIT, Littleton Lab UROP
# 24 January 2019
# ---------------------
# Descripton: A script that takes the most significant 


####################################################
rm(list=ls())
library(protr)



data <- read.csv("C:/UROPs/polyQ_neuronal_proteins/output/fly/Q/fly_22905prots_w_HMM_Q.csv")

seq <- as.character(data$peptide_seq[[1]])
aac <- extractAAC(seq)
aac["Q"]
dc <- extractDC(seq)
dc["QQ"]
tc <- extractTC(seq)
tc["QQQ"]


moreau = extractMoreauBroto(seq)
head(moreau, n = 36L)


extractCTDC(seq)
