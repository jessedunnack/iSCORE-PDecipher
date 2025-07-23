#!/bin/bash

# Wrapper script to run expression processing with correct conda environment

echo "Activating seuratv4 conda environment..."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate seuratv4

echo "R version:"
R --version | head -1

echo ""
echo "Running expression processing..."
Rscript process_all_umap_expression.R

echo ""
echo "Processing complete!"