#!/usr/bin/env Rscript

# Script to rebuild package documentation
setwd("/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PDecipher")
devtools::document()
cat("Documentation rebuild complete\n")