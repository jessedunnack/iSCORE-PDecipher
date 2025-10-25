#!/bin/bash
# Quick script to re-run enrichment consolidation
# Run this after downloading new enrichment results from HPC

cd /mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PDecipher

echo "=================================================="
echo "  ENRICHMENT CONSOLIDATION SCRIPT"
echo "=================================================="
echo ""
echo "This will consolidate enrichment results from:"
echo "  - enrichment_results_FPD_p_weight"
echo "  - enrichment_results_FPD_p_weight_BH"
echo "  - enrichment_results_FPD_p_weight_bonferroni"
echo "  - enrichment_results_CRISPRi_p_weight"
echo "  - enrichment_results_CRISPRi_p_weight_BH"
echo "  - enrichment_results_CRISPRi_p_weight_bonferroni"
echo ""
echo "Output: all_enrichment_padj005_complete_with_direction.rds in each directory"
echo ""
echo "Estimated time: 1-2 hours for all 6 directories"
echo ""
read -p "Press ENTER to continue or Ctrl+C to cancel..."

echo ""
echo "Starting consolidation with seuratv4 environment..."
echo "Log: consolidate_enrichment_log.txt"
echo ""

conda run -n seuratv4 Rscript consolidate_pooled_enrichment.R 2>&1 | tee consolidate_enrichment_log.txt

echo ""
echo "=================================================="
echo "  CONSOLIDATION COMPLETE!"
echo "=================================================="
echo ""
echo "Created files:"
ls -lh ../final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/enrichment_results_*/all_enrichment_padj005_complete_with_direction.rds 2>/dev/null
echo ""
