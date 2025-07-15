#!/usr/bin/env Rscript

#' Diagnostic test to check why we're getting 0 overlaps
#' Tests with increasingly lenient thresholds

# Load DE data
de_data_path <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/full_DE_results.rds"
full_de <- readRDS(de_data_path)

# Test LRRK2 in cluster_0
mast_data <- full_de$iSCORE_PD_MAST[["LRRK2"]][["cluster_0"]]$results
crispri_data <- full_de$CRISPRi_Mixscale[["LRRK2"]][["cluster_0"]]$results

cat("=== DIAGNOSTIC CHECK FOR LRRK2 CLUSTER_0 ===\n\n")

# Check MAST data structure
cat("MAST data structure:\n")
cat("  Dimensions:", nrow(mast_data), "x", ncol(mast_data), "\n")
cat("  Column names:", paste(head(names(mast_data), 10), collapse=", "), "\n")
cat("  Sample rows:\n")
print(head(mast_data[, c("avg_log2FC", "p_val_adj")], 5))

# Check CRISPRi data structure
cat("\n\nCRISPRi data structure:\n")
cat("  Dimensions:", nrow(crispri_data), "x", ncol(crispri_data), "\n")
cat("  Column names:", paste(head(names(crispri_data), 10), collapse=", "), "\n")

# Find log2FC columns
log2fc_cols <- grep("^log2FC_", names(crispri_data), value = TRUE)
cat("  Log2FC columns found:", paste(log2fc_cols, collapse=", "), "\n")

# Check for different thresholds
thresholds <- c(0.05, 0.1, 0.25, 0.5, 1.0)

cat("\n\nSignificant genes at different thresholds:\n")
cat("=========================================\n")

for (p_thresh in c(0.01, 0.05, 0.1)) {
  for (lfc_thresh in thresholds) {
    # MAST significant genes
    mast_sig_up <- sum(mast_data$avg_log2FC > lfc_thresh & mast_data$p_val_adj < p_thresh, na.rm = TRUE)
    mast_sig_down <- sum(mast_data$avg_log2FC < -lfc_thresh & mast_data$p_val_adj < p_thresh, na.rm = TRUE)
    
    # CRISPRi significant genes (use first experiment)
    if (length(log2fc_cols) > 0) {
      lfc_col <- log2fc_cols[1]
      crispri_sig_up <- sum(crispri_data[[lfc_col]] > lfc_thresh & crispri_data$p_val_adj < p_thresh, na.rm = TRUE)
      crispri_sig_down <- sum(crispri_data[[lfc_col]] < -lfc_thresh & crispri_data$p_val_adj < p_thresh, na.rm = TRUE)
    } else {
      crispri_sig_up <- crispri_sig_down <- 0
    }
    
    cat(sprintf("p < %.2f, |log2FC| > %.2f: MAST (up:%d, down:%d), CRISPRi (up:%d, down:%d)\n",
                p_thresh, lfc_thresh, mast_sig_up, mast_sig_down, crispri_sig_up, crispri_sig_down))
  }
  cat("\n")
}

# Check top DE genes
cat("\nTop 10 DE genes in MAST (by p-value):\n")
mast_top <- mast_data[order(mast_data$p_val_adj), ]
print(head(mast_top[, c("avg_log2FC", "p_val_adj")], 10))

cat("\n\nTop 10 DE genes in CRISPRi (by p-value):\n")
if ("p_val_adj" %in% names(crispri_data)) {
  crispri_top <- crispri_data[order(crispri_data$p_val_adj), ]
  if (length(log2fc_cols) > 0) {
    print(head(crispri_top[, c(log2fc_cols[1], "p_val_adj")], 10))
  }
} else {
  cat("No p_val_adj column found in CRISPRi data\n")
}

# Check gene name overlap
cat("\n\nGene name overlap check:\n")
mast_genes <- rownames(mast_data)
crispri_genes <- if (!is.null(rownames(crispri_data))) rownames(crispri_data) else crispri_data$gene
overlap <- intersect(mast_genes, crispri_genes)
cat("  MAST genes:", length(mast_genes), "\n")
cat("  CRISPRi genes:", length(crispri_genes), "\n")
cat("  Overlapping genes:", length(overlap), "\n")
cat("  Example overlapping genes:", paste(head(overlap, 10), collapse=", "), "\n")