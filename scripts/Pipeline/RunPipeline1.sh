#!/bin/bash
set -e
echo "Running RunPipeline.sh script!"

Rscript.exe -e ".libPaths()"

kModels="adjusted"
kSpecies="fly"
kCandidateAAs="DTSEPGACVMILYFHKRWQN"

echo "Running pipeline/part1/DownloadProteomes.R for species: $kSpecies"
Rscript pipeline/part1/DownloadProteomes.R $kSpecies

echo "Running pipeline/part1/HMM.R with HMM_models: $kModels and species: $kSpecies and AAs: $kCandidateAAs"
Rscript pipeline/part1/HMM.R $kModels $kSpecies $kCandidateAAs

echo "Running pipeline/part1/HmmSummaryPlots.R with HMM_models: $kModels and species: $kSpecies and AAs: $kCandidateAAs"
Rscript pipeline/part1/HmmSummaryPlots.R $kModels $kSpecies $kCandidateAAs

