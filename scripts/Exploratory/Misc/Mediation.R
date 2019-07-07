####################################################
# Daniel Alfonsetti, daniel.alfonsetti@gmail.com
# MIT, Littleton Lab
# 24 January 2019
# -----------------
# Description: A script that performs mediation analysis to see if 
# protein length is a mediator of the correlation between 
# a protein being marked as polyQ by our HMM (a boolean var) and whether 
# the protein is in a synaptic GO category.
# If length is a mediator, then it is evidence that length 
# is something that should be controlled for.
####################################################
rm(list = ls())
library(bda) # Mediation library
####################################################

data <- read.csv("C:/UROPs/polyQ_neuronal_proteins/output/fly/Q/fly_22905prots_w_HMM_Q.csv")
#data <- read.csv("C:/UROPs/polyQ_neuronal_proteins/output/human/Q/human_78328prots_w_HMM_Q.csv")


# Annotate proteins based on whether they are a synaptic protein or not.
synapse_term1 <- "GO:0045202"
synapse_cats <- get(synapse_term1, GOCCOFFSPRING)
synapse_cats <- c(synapse_cats, synapse_term1)

data$synaptic <- unlist(
  lapply(
    lapply(as.character(data$go_ids_CC), function(x){unlist(strsplit(x, "; "))}), 
    function(x){any(x %in% synapse_cats)}
  )
)
data$peptide_length <- unlist(lapply(as.character(data$peptide_seq), 
                                     function(x){length(unlist(strsplit(x, "")))}))



# Let X, the independent variable, be whether the protein is a synaptic protein or not (bool)
# Let Y, the dependent variable, be whether the protein is polyQ or not (bool)
# Let M, the candidate mediator, be the length of the protein (num)

model.0 <- lm(data$HMMhasPolyAA ~ data$synaptic)
summary(model.0) # 7.18e-10
# X is a significant predictor of Y

model.1 <- lm(data$peptide_length ~ data$synaptic)
summary(model.1)
# X is a significant predictor of M

model.2 <- lm(data$HMMhasPolyAA ~ data$synaptic + data$peptide_length)
summary(model.2) # p-value: 7.48e-05. 
# Predicting ability of synaptic variable does drop, but is it significant?

# Let's test it.
mediation.test(mv = data$peptide_length, iv = data$synaptic, dv = data$HMMhasPolyAA)
# Mediation of peptide length is indeed significant. 
# Protein length is thus something that should be controlled for.


# What about the other way around?
mediation.test(mv = data$peptide_length, dv = data$synaptic, iv = data$HMMhasPolyAA)
# Yes, length does mediate the polyQ -> synapse effect, but not as strongly as it mediates the 
# the synapse -> polyQ effect.

####################################################
# Sanity check to make sure mediation test is actually right.
fake_mediator <- sample(c(0,1), replace=TRUE, size=nrow(data))
mediation.test(mv = fake_mediator, iv = data$synaptic, dv = data$HMMhasPolyAA)
# No significant mediation as expected. Good.


