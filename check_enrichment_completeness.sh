#!/bin/bash
# Check enrichment completeness across all directories

BASE="/mnt/e/ASAP/scRNASeq/PerturbSeq/final/final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper"

echo "==============================================="
echo "CHECKING ENRICHMENT DIRECTORY COMPLETENESS"
echo "==============================================="
echo ""

# Expected methods
EXPECTED_METHODS=("GO_ALL" "GO_BP" "GO_CC" "GO_MF" "KEGG" "Reactome" "WikiPathways" "STRING" "GSEA")
EXPECTED_COUNT=9

for ENRICH_DIR in "$BASE"/enrichment_results_*; do
    DIR_NAME=$(basename "$ENRICH_DIR")
    echo "=== $DIR_NAME ==="

    # Count clusters
    CLUSTER_COUNT=$(find "$ENRICH_DIR" -maxdepth 1 -type d -name "all_*_clust_*" | wc -l)
    echo "Clusters: $CLUSTER_COUNT"

    # Check each cluster
    INCOMPLETE_COUNT=0
    TOTAL_PERTS=0

    for CLUSTER_DIR in "$ENRICH_DIR"/all_*_clust_*; do
        if [ ! -d "$CLUSTER_DIR" ]; then continue; fi

        CLUSTER_NAME=$(basename "$CLUSTER_DIR" | sed 's/.*_clust_/cluster_/')

        # Count perturbations
        PERT_COUNT=$(find "$CLUSTER_DIR" -maxdepth 1 -type d | grep -v "^$CLUSTER_DIR$" | wc -l)
        TOTAL_PERTS=$((TOTAL_PERTS + PERT_COUNT))

        # Check each perturbation for complete methods
        for PERT_DIR in "$CLUSTER_DIR"/*; do
            if [ ! -d "$PERT_DIR" ]; then continue; fi

            PERT_NAME=$(basename "$PERT_DIR")

            # Count method directories (exclude diagnostics)
            METHOD_COUNT=$(find "$PERT_DIR" -maxdepth 1 -type d | grep -v diagnostics | grep -v "^$PERT_DIR$" | wc -l)

            if [ "$METHOD_COUNT" -ne "$EXPECTED_COUNT" ]; then
                echo "  ⚠️  $CLUSTER_NAME/$PERT_NAME: $METHOD_COUNT methods (expected $EXPECTED_COUNT)"
                INCOMPLETE_COUNT=$((INCOMPLETE_COUNT + 1))

                # Show which methods are missing
                for METHOD in "${EXPECTED_METHODS[@]}"; do
                    if [ ! -d "$PERT_DIR/$METHOD" ]; then
                        echo "      Missing: $METHOD"
                    fi
                done
            fi
        done
    done

    echo "Total perturbations: $TOTAL_PERTS"
    if [ "$INCOMPLETE_COUNT" -eq 0 ]; then
        echo "✓ All perturbations have complete enrichment methods"
    else
        echo "⚠️  $INCOMPLETE_COUNT perturbations have incomplete methods"
    fi
    echo ""
done

echo "==============================================="
echo "CHECK COMPLETE"
echo "==============================================="
