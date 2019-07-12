#!/bin/bash
set -e

kModels="adjusted"
kSpecies="fly"
kCandidateAAs="DTSEPGACVMILYFHKRWQN"

echo "Running pipeline/part2/EnrichmentAnalysis.R with HMM_models: $kModels and species: $kSpecies and AAs: $kCandidateAAs"
(cd part2; Rscript EnrichmentAnalysis.R $kModels $kSpecies $kCandidateAAs)

echo "Running pipeline/part2/AminoAcidCharts.R with HMM_models"
(cd part2; Rscript AminoAcidCharts.R $kModels $kSpecies)

echo "Running pipeline/part2/AminoAcidFisherTestHistograms.R"
(cd part2; Rscript AminoAcidFisherTestHistograms.R)

