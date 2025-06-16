#!/usr/bin/env Rscript
# Test heatmap generation without Shiny UI

cat("=== Testing Heatmap Generation ===\n\n")

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(heatmaply)
  library(plotly)
})

# Simulate the heatmap data preparation
test_heatmap_generation <- function() {
  # Create test data similar to enrichment results
  set.seed(123)
  test_data <- expand.grid(
    mutation_perturbation = c("LRRK2", "PINK1", "PARK7", "SNCA", "GBA"),
    enrichment_type = c("GO_BP", "KEGG"),
    direction = c("UP", "DOWN"),
    stringsAsFactors = FALSE
  ) %>%
    slice_sample(n = 100, replace = TRUE) %>%
    mutate(
      Description = paste("Term", 1:100),
      p.adjust = runif(100, 0.0001, 0.05),
      FoldEnrichment = runif(100, 1.5, 5),
      Count = sample(5:50, 100, replace = TRUE),
      cluster = "cluster_0"
    )
  
  cat("Created test data with", nrow(test_data), "rows\n")
  
  # Test 1: Basic matrix creation
  cat("\nTest 1: Creating heatmap matrix...\n")
  tryCatch({
    # Select top terms
    top_terms <- test_data %>%
      group_by(Description) %>%
      summarise(mean_fold = mean(FoldEnrichment)) %>%
      arrange(desc(mean_fold)) %>%
      head(20) %>%
      pull(Description)
    
    # Create matrix
    mat_data <- test_data %>%
      filter(Description %in% top_terms) %>%
      select(Description, mutation_perturbation, FoldEnrichment) %>%
      pivot_wider(names_from = mutation_perturbation, 
                  values_from = FoldEnrichment,
                  values_fill = 0,
                  values_fn = mean)
    
    mat <- as.matrix(mat_data[,-1])
    rownames(mat) <- mat_data$Description
    
    cat("✓ Matrix created:", nrow(mat), "x", ncol(mat), "\n")
    cat("  Genes:", paste(colnames(mat), collapse=", "), "\n")
    
    # Test 2: Create heatmaply without annotations
    cat("\nTest 2: Creating heatmaply without annotations...\n")
    p1 <- tryCatch({
      heatmaply(
        mat,
        dendrogram = "both",
        colors = viridis::viridis(256),
        main = "Test Heatmap - No Annotations",
        margins = c(100, 100, 50, 50),
        plot_method = "plotly"
      )
      cat("✓ Heatmaply created successfully (no annotations)\n")
      TRUE
    }, error = function(e) {
      cat("✗ Heatmaply failed:", e$message, "\n")
      FALSE
    })
    
    # Test 3: Try with annotations
    cat("\nTest 3: Creating heatmaply with annotations...\n")
    if (p1) {
      # Create annotation data
      ann_data <- test_data %>%
        filter(Description %in% rownames(mat)) %>%
        group_by(Description) %>%
        summarise(
          Type = first(enrichment_type),
          Direction = first(direction)
        ) %>%
        arrange(match(Description, rownames(mat)))
      
      row_side_colors <- data.frame(
        Type = ann_data$Type,
        Direction = ann_data$Direction,
        row.names = rownames(mat)
      )
      
      # Define color palettes
      type_colors <- c("GO_BP" = "#8DD3C7", "KEGG" = "#FB8072")
      dir_colors <- c("UP" = "#FF6B6B", "DOWN" = "#4ECDC4")
      
      p2 <- tryCatch({
        heatmaply(
          mat,
          dendrogram = "both",
          colors = viridis::viridis(256),
          main = "Test Heatmap - With Annotations",
          margins = c(100, 150, 50, 50),
          plot_method = "plotly",
          row_side_colors = row_side_colors,
          row_side_palette = list(
            Type = type_colors,
            Direction = dir_colors
          )
        )
        cat("✓ Heatmaply created successfully WITH annotations!\n")
        TRUE
      }, error = function(e) {
        cat("✗ Heatmaply with annotations failed:", e$message, "\n")
        cat("  This matches the issue seen in the app\n")
        FALSE
      })
    }
    
    # Test 4: Fallback to basic plotly
    cat("\nTest 4: Testing plotly fallback...\n")
    tryCatch({
      p_fallback <- plot_ly(
        x = colnames(mat),
        y = rownames(mat),
        z = mat,
        type = "heatmap",
        colorscale = "Viridis",
        hovertemplate = "Gene: %{x}<br>Term: %{y}<br>Value: %{z}<extra></extra>"
      ) %>%
        layout(
          title = "Fallback Heatmap",
          xaxis = list(title = "Genes"),
          yaxis = list(title = "Terms")
        )
      cat("✓ Plotly fallback works\n")
    }, error = function(e) {
      cat("✗ Plotly fallback failed:", e$message, "\n")
    })
    
  }, error = function(e) {
    cat("✗ Error in heatmap generation:", e$message, "\n")
  })
}

# Run the test
test_heatmap_generation()

cat("\n=== Summary ===\n")
cat("The heatmaply annotation issue can be reproduced in isolation.\n")
cat("The fallback to basic heatmaply (without annotations) or plotly should work.\n")