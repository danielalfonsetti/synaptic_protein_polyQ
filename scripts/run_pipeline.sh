#!/bin/bash

echo "Running download_proteomes.R"
Rscript pipeline/download_proteomes.R

echo "Running Enrichment_Analysis.R"
Rscript pipeline/Enrichment_Analaysis.R

echo "Running HMM.R"
Rscript pipeline/HMM.R