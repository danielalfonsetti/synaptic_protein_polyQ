Folders in this directory represent outputs for different hidden markov models.

Trained - Results using HMM trained on hand selected polyQ synaptic proteins.
Adjusted - Results using manually adjusted version of the trained HMM (increased the emission frequency of the targetAA to decrease false positive rate). This version seems to do much better (see './adjusted/AA_charts')

We are currently using the mh-l version of neuronal transcriptome for transcriptome filtering.
Transcriptome data was accessed through https://flybase.org/rnaseq/profile_search.

mh-l stands for 'moderately high - low', which is in reference to the lowest peak level of expression of proteins in tissues that we want to exclude needed for us not exclude the protein (i.e. moderately high) 
and the peak level of expression of proteins in tissues that we want (CNS tissues) in order for us to not throw away (e.g. low).

So, we are only including everything that is at least 'lowly' expressed in the CNS tissues while simutaneously not being moderately expressed or more in non-CNS tissues.