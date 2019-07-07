#!/bin/bash
set -e 

echo "Running test.R script!"
Rscript test.R

echo "Running Pipeline/test2.R script!"
#Rscript Pipeline/test2.R
(cd Pipeline/ && Rscript test2.R)
