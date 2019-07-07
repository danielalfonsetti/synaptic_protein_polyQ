polyQ Synaptic Protein Project
MIT Littleton Lab UROP Fall 18, Spring 19
Daniel Alfonsetti
daniel.alfonsetti@gmail.com; alfonset@mit.edu 
22 April 2019

-------------------------------
- Research Question
-------------------------------
Are proteins with polyQ tracks enriched in the pre-synaptic active zone? What about the post synaptic density?

-------------------------------
- Goal
-------------------------------
Make a computational pipeline that is general enough to answer this question and related questions such as...

For an arbitrary amino acid AA and an arbitrary cellulular location X or functional category Y...
1. Are proteins with polyAA tracks enriched at location X?
2. Are proteins with polyAA tracks enriched in the set of Y-related proteins?

Do our finds vary by species?

-------------------------------
- Procedure Overview
-------------------------------
Train hidden markov model on known polyQ proteins and then find polyAA tracks (Viterbi method).
Take the proteins annotated to have polyAA tracks and test to see if they are enriched for certain functional categories or groups (Fisher Tests, enrichment analysis).

Repeat enrichment testing procedures after filtering on only neuronally transcribed proteins and non-nulcear proteins (the latter of which are known to have polyQ regions).

-------------------------------
- Style Guide for this Project
-------------------------------
Identifier Conventions: 
  Functions: MyFunction(...)
  Constants: kMyConstant
  Variables: myVariables
  https://google.github.io/styleguide/Rguide.xml

-------------------------------
Made in...
R version 3.5.2 (2018-12-20) -- "Eggshell Igloo"


