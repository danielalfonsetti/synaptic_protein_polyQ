#!/bin/bash
set -e
echo "Running RunPipeline.sh script!"

Rscript.exe -e ".libPaths()"

kModels="adjusted"
kSpecies="fly"
kCandidateAAs="DTSEPGACVMILYFHKRWQN"


echo "Running pipeline/DownloadProteomes.R for species: $kSpecies"
Rscript pipeline/DownloadProteomes.R $kSpecies

echo "Running pipeline/HMM.R with HMM_models: $kModels and species: $kSpecies and AAs: $kCandidateAAs"
Rscript pipeline/HMM.R $kModels $kSpecies $kCandidateAAs

echo "Running pipeline/EnrichmentAnalysis.R with HMM_models: $kModels and species: $kSpecies and AAs: $kCandidateAAs"
Rscript pipeline/EnrichmentAnalysis.R $kModels $kSpecies $kCandidateAAs

echo "Running pipeline/HmmSummaryPlots.R with HMM_models: $kModels and species: $kSpecies and AAs: $kCandidateAAs"
Rscript pipeline/HmmSummaryPlots.R $kModels $kSpecies $kCandidateAAs

