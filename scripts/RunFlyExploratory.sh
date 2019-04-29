#!/bin/bash
set -e
echo "Running RunFlyExploratory.sh script!"

Rscript.exe -e ".libPaths()"

kModels="adjusted"
kSpecies="fly"
kCandidateAAs="DTSEPGACVMILYFHKRWQN"

echo "Running Exploratory/AminoAcidCharts.R with HMM_models"
Rscript Exploratory/AminoAcidCharts.R $kModels $kSpecies

echo "Running Exploratory/AminoAcidFisherTestHistograms.R"
Rscript Exploratory/AminoAcidFisherTestHistograms.R