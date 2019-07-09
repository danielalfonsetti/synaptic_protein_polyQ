#!/bin/bash
set -e
echo "Running RunPipeline.sh script! "

kModels="adjusted"
kSpecies="fly"
kCandidateAAs="DTSEPGACVMILYFHKRWQN"

echo "Running scripts/pipeline/part1/DownloadProteomes.R for species: $kSpecies"
(cd part1; Rscript DownloadProteomes.R $kSpecies)

echo "Running scripts/pipeline/part1/HMM.R with HMM_models: $kModels and species: $kSpecies and AAs: $kCandidateAAs"
(cd part1; Rscript HMM.R $kModels $kSpecies $kCandidateAAs)

echo "Running scripts/pipeline/part1/HmmSummaryPlots.R with HMM_models: $kModels and species: $kSpecies and AAs: $kCandidateAAs"
(cd part1; Rscript HmmSummaryPlots.R $kModels $kSpecies $kCandidateAAs)

