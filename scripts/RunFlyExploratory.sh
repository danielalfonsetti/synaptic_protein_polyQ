#!/bin/bash
set -e
echo "Running RunFlyExploratory.sh script!"

Rscript.exe -e ".libPaths()"

echo "Running Exploratory/AminoAcidCharts.R with HMM_models"
Rscript Exploratory/AminoAcidCharts.R

echo "Running Exploratory/AminoAcidFisherTestHistograms.R"
Rscript Exploratory/AminoAcidFisherTestHistograms.R