#!/bin/bash
# Inventory all enrichment files to check completeness before consolidation

BASE="/mnt/e/ASAP/scRNASeq/PerturbSeq/final/final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper"

echo "========================================================================"
echo "  ENRICHMENT FILES INVENTORY"
echo "========================================================================"
echo ""

# Expected methods
EXPECTED_METHODS=("GO_ALL" "GO_BP" "GO_CC" "GO_MF" "KEGG" "Reactome" "WikiPathways" "STRING" "GSEA")
EXPECTED_DIRECTIONS=("ALL" "UP" "DOWN")

for ENRICH_DIR in "$BASE"/enrichment_results_*; do
    DIR_NAME=$(basename "$ENRICH_DIR")

    echo "========================================================================"
    echo "DATASET: $DIR_NAME"
    echo "========================================================================"

    # Count clusters
    CLUSTER_DIRS=$(find "$ENRICH_DIR" -maxdepth 1 -type d -name "all_*_clust_*" | sort)
    CLUSTER_COUNT=$(echo "$CLUSTER_DIRS" | grep -c "clust_")

    echo "Clusters found: $CLUSTER_COUNT"
    echo ""

    TOTAL_PERTS=0
    TOTAL_METHODS=0
    TOTAL_RDS=0
    INCOMPLETE_PERTS=0
    EMPTY_PERTS=0

    for CLUSTER_DIR in $CLUSTER_DIRS; do
        CLUSTER_NAME=$(basename "$CLUSTER_DIR" | sed 's/.*_clust_/cluster_/')

        # Count perturbations in this cluster
        PERT_DIRS=$(find "$CLUSTER_DIR" -maxdepth 1 -type d | grep -v "^$CLUSTER_DIR$")
        PERT_COUNT=$(echo "$PERT_DIRS" | grep -v "^$" | wc -l)
        TOTAL_PERTS=$((TOTAL_PERTS + PERT_COUNT))

        echo "  $CLUSTER_NAME: $PERT_COUNT perturbations"

        # Check each perturbation
        CLUSTER_INCOMPLETE=0
        CLUSTER_EMPTY=0

        for PERT_DIR in $PERT_DIRS; do
            if [ -z "$PERT_DIR" ] || [ ! -d "$PERT_DIR" ]; then
                continue
            fi

            # Count method directories
            METHOD_COUNT=$(find "$PERT_DIR" -maxdepth 1 -type d | grep -v diagnostics | grep -v "^$PERT_DIR$" | wc -l)
            TOTAL_METHODS=$((TOTAL_METHODS + METHOD_COUNT))

            # Count RDS files
            RDS_COUNT=$(find "$PERT_DIR" -name "*.rds" | wc -l)
            TOTAL_RDS=$((TOTAL_RDS + RDS_COUNT))

            # Check completeness
            if [ "$METHOD_COUNT" -eq 0 ]; then
                CLUSTER_EMPTY=$((CLUSTER_EMPTY + 1))
                EMPTY_PERTS=$((EMPTY_PERTS + 1))
            elif [ "$METHOD_COUNT" -lt 9 ]; then
                CLUSTER_INCOMPLETE=$((CLUSTER_INCOMPLETE + 1))
                INCOMPLETE_PERTS=$((INCOMPLETE_PERTS + 1))
            fi
        done

        if [ "$CLUSTER_INCOMPLETE" -gt 0 ] || [ "$CLUSTER_EMPTY" -gt 0 ]; then
            echo "    ⚠️  Incomplete: $CLUSTER_INCOMPLETE | Empty: $CLUSTER_EMPTY"
        fi
    done

    echo ""
    echo "SUMMARY:"
    echo "  Total perturbations: $TOTAL_PERTS"
    echo "  Total method directories: $TOTAL_METHODS"
    echo "  Total RDS files: $TOTAL_RDS"
    echo "  Incomplete perturbations: $INCOMPLETE_PERTS"
    echo "  Empty perturbations: $EMPTY_PERTS"

    # Calculate expected vs actual
    EXPECTED_RDS=$((TOTAL_PERTS * 9 * 3))  # 9 methods × 3 directions
    COMPLETENESS=$((TOTAL_RDS * 100 / EXPECTED_RDS))

    echo "  Expected RDS files: $EXPECTED_RDS"
    echo "  Completeness: ${COMPLETENESS}%"

    if [ "$COMPLETENESS" -ge 90 ]; then
        echo "  ✅ READY for consolidation"
    elif [ "$COMPLETENESS" -ge 50 ]; then
        echo "  ⚠️  PARTIAL - may be missing some HPC results"
    else
        echo "  ❌ INCOMPLETE - many HPC results missing"
    fi

    echo ""
done

echo "========================================================================"
echo "INVENTORY COMPLETE"
echo "========================================================================"
