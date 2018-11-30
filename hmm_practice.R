# ## Required library
# library(depmixS4)
# 
# ## data loading-
# physician_prescrition_data <-c(12,16,45,45,56,67,78,98,120,124,156)
# 
# ## model execution-
# HMM_model <- depmixS4::depmix(physician_prescrition_data~1, nstates = 2,ntimes=length(physician_prescrition_data))
# 
# ## model fitting
# HMM_fm <- fit(HMM_model)
# 
# ## Transition probabilties-
# HMM_fm@transition
# 
# ## posterior states-
# posterior(HMM_fm)
# plot(ts(posterior(HMM_fm)[,1]))
# 
# ## Emission probabilties-
# HMM_fm@response
# 
library(seqHMM)
library(TraMineR)

# observations = seqdef( )

data(actcal)
actcal.seq <- seqdef(actcal,13:24,
                     labels=c("> 37 hours", "19-36 hours", "1-18 hours", "no work"))

View(actcal)
View(actcal.seq)

hmm_obj <- build_hmm(observations = actcal.seq,
          n_states = 4)

hmm <- fit_model(hmm_obj)
hmm$model

plot(
  hmm$model,
  # larger vertices
  vertex.size = 45,
  # varying curvature of edges
  edge.curved = c(0, -0.7, 0.6, 0, -0.7, 0),
  # legend with two columns and less space
  ncol.legend = 2, legend.prop = 0.4,
  # new label for combined slice
  combined.slice.label = "States with probability < 0.05")

ssplot(actcal.seq)

fb <- forward_backward(hmm_obj)
apply(fb$forward_probs[, , 3], 2, which.max)


############33
data(ex1)
str(ex1)
ex1.seq <- seqdef(ex1, 1:13, right = "DEL")
ssplot(ex1.seq)
hmm <- build_hmm(ex1.seq, n_states = 2)
hmm
fit_hmm <- fit_model(hmm)
fit_hmm$model
str(fit_hmm)
plot(fit_hmm$model)

ssplot(ex1.seq)


data(hmm_biofam)
plot(hmm_biofam)
