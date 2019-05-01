#!/bin/bash
set -e
echo "Running RunFlyExploratory.sh script!"

Rscript.exe -e ".libPaths()"

kModels="adjusted"
kSpecies="fly"
kCandidateAAs="DTSEPGACVMILYFHKRWQN"


echo "Running pipeline/part2/EnrichmentAnalysis.R with HMM_models: $kModels and species: $kSpecies and AAs: $kCandidateAAs"
Rscript Pipeline/part2/EnrichmentAnalysis.R $kModels $kSpecies $kCandidateAAs

echo "Running pipeline/part2/AminoAcidCharts.R with HMM_models"
Rscript Pipeline/part2/AminoAcidCharts.R $kModels $kSpecies

echo "Running pipeline/part2/AminoAcidFisherTestHistograms.R"
Rscript Pipeline/part2/AminoAcidFisherTestHistograms.R

