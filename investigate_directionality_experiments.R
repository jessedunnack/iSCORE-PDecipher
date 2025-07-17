# Investigate Directionality vs Multiple CRISPRi Experiments
# This script checks if directionality "inflation" is due to multiple experiments

library(dplyr)

# Find consolidated data
find_consolidated_data <- function() {
  data_paths <- c(
    "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi",
    "../../iSCORE-PD_plus_CRISPRi",
    "../iSCORE-PD_plus_CRISPRi"
  )
  
  for (path in data_paths) {
    test_file <- file.path(path, "all_enrichment_padj005_complete_with_direction.rds")
    if (file.exists(test_file)) {
      return(test_file)
    }
  }
  return(NULL)
}

cat("=== INVESTIGATING DIRECTIONALITY vs EXPERIMENTS ===\n")

# Load data
consolidated_file <- find_consolidated_data()
if (is.null(consolidated_file)) {
  stop("Could not find consolidated data file!")
}

cat("Loading data from:", consolidated_file, "\n")
data <- readRDS(consolidated_file)

# Check for experiment column
if (!"experiment" %in% colnames(data)) {
  cat("❌ No 'experiment' column found in data\n")
  cat("Available columns:", paste(colnames(data), collapse = ", "), "\n")
  stop("Cannot investigate experiment-related directionality")
}

cat("✅ Found experiment column\n")

# Check unique experiments
cat("\nUnique experiments in data:\n")
experiments <- unique(data$experiment)
print(experiments)

# Focus on CRISPRi data for this analysis
crispri_data <- data %>%
  filter(grepl("MixScale", method) | grepl("CRISPRi", method))

cat("\nCRISPRi data dimensions:", nrow(crispri_data), "rows\n")

# 1. Check if same pathway has different directions across experiments
cat("\n1. ANALYZING DIRECTIONALITY ACROSS EXPERIMENTS\n")
cat("=============================================\n")

# Look for cases where same gene/cluster/pathway has different directions across experiments
direction_analysis <- crispri_data %>%
  filter(!is.na(mutation_perturbation), !is.na(cluster), !is.na(Description), !is.na(experiment)) %>%
  group_by(mutation_perturbation, cluster, enrichment_type, Description) %>%
  summarise(
    n_experiments = n_distinct(experiment),
    experiments = paste(unique(experiment), collapse = ","),
    n_directions = n_distinct(direction),
    directions = paste(unique(direction), collapse = ","),
    experiment_direction_pairs = paste(paste(experiment, direction, sep = ":"), collapse = ","),
    .groups = "drop"
  ) %>%
  filter(n_directions > 1)  # Multiple directions for same pathway

cat("Pathways with multiple directions:", nrow(direction_analysis), "\n")

if (nrow(direction_analysis) > 0) {
  # Sample some cases to examine
  cat("\nSample cases of directional differences:\n")
  cat("======================================\n")
  
  sample_cases <- head(direction_analysis, 5)
  for (i in 1:nrow(sample_cases)) {
    case <- sample_cases[i, ]
    cat(sprintf("%d. Gene: %s | Cluster: %s\n", i, case$mutation_perturbation, case$cluster))
    cat(sprintf("   Pathway: %s\n", case$Description))
    cat(sprintf("   Experiments: %s\n", case$experiments))
    cat(sprintf("   Directions: %s\n", case$directions))
    cat(sprintf("   Experiment:Direction pairs: %s\n", case$experiment_direction_pairs))
    cat("\n")
  }
  
  # Check if different directions come from different experiments
  direction_experiment_check <- crispri_data %>%
    filter(!is.na(mutation_perturbation), !is.na(cluster), !is.na(Description), !is.na(experiment)) %>%
    group_by(mutation_perturbation, cluster, enrichment_type, Description, direction) %>%
    summarise(
      experiments_for_this_direction = paste(unique(experiment), collapse = ","),
      n_experiments_for_this_direction = n_distinct(experiment),
      .groups = "drop"
    )
  
  # Check for pathways where same direction appears in multiple experiments
  same_direction_multiple_experiments <- direction_experiment_check %>%
    filter(n_experiments_for_this_direction > 1)
  
  cat("Cases where SAME direction appears in MULTIPLE experiments:", nrow(same_direction_multiple_experiments), "\n")
  
  if (nrow(same_direction_multiple_experiments) > 0) {
    cat("This suggests some inflation may still be due to data processing issues.\n")
    cat("Sample cases:\n")
    print(head(same_direction_multiple_experiments, 3))
  }
  
} else {
  cat("✅ No pathways found with multiple directions!\n")
  cat("This suggests the directionality 'inflation' is NOT due to experiment differences.\n")
}

# 2. Check experiment-specific direction patterns
cat("\n\n2. EXPERIMENT-SPECIFIC DIRECTION PATTERNS\n")
cat("=========================================\n")

experiment_direction_summary <- crispri_data %>%
  group_by(experiment, direction) %>%
  summarise(
    n_terms = n(),
    n_genes = n_distinct(mutation_perturbation),
    n_clusters = n_distinct(cluster),
    .groups = "drop"
  ) %>%
  arrange(experiment, direction)

cat("Direction distribution by experiment:\n")
print(experiment_direction_summary)

# 3. Look for specific examples with LRRK2 (user's test case)
cat("\n\n3. LRRK2 EXAMPLE ANALYSIS\n")
cat("========================\n")

lrrk2_crispri <- crispri_data %>%
  filter(mutation_perturbation == "LRRK2", cluster == "cluster_0") %>%
  group_by(enrichment_type, Description, direction) %>%
  summarise(
    experiments = paste(unique(experiment), collapse = ","),
    n_experiments = n_distinct(experiment),
    .groups = "drop"
  ) %>%
  arrange(Description, direction)

cat("LRRK2 cluster_0 CRISPRi examples:\n")
print(head(lrrk2_crispri, 10))

# Check for LRRK2 pathways appearing in multiple directions
lrrk2_multi_direction <- lrrk2_crispri %>%
  group_by(enrichment_type, Description) %>%
  summarise(
    n_directions = n_distinct(direction),
    directions = paste(unique(direction), collapse = ","),
    .groups = "drop"
  ) %>%
  filter(n_directions > 1)

cat("\nLRRK2 pathways with multiple directions:", nrow(lrrk2_multi_direction), "\n")

if (nrow(lrrk2_multi_direction) > 0) {
  cat("Examples:\n")
  print(head(lrrk2_multi_direction, 5))
}

# 4. CONCLUSION AND RECOMMENDATIONS
cat("\n\n4. CONCLUSION AND RECOMMENDATIONS\n")
cat("=================================\n")

if (nrow(direction_analysis) > 0) {
  cat("🔍 FINDINGS:\n")
  cat("• Multiple directions found for same pathways\n")
  cat("• This could be due to:\n")
  cat("  1. Different CRISPRi experiments having different effects\n")
  cat("  2. Data processing combining results inappropriately\n")
  cat("  3. ALL direction including both UP and DOWN terms\n")
  
  # Check if the issue is experiment-specific
  exp_specific_directions <- direction_analysis %>%
    mutate(
      has_experiment_specific_directions = !grepl(",", experiment_direction_pairs)
    )
  
  n_exp_specific <- sum(exp_specific_directions$has_experiment_specific_directions)
  
  if (n_exp_specific > 0) {
    cat(sprintf("\n✅ GOOD NEWS: %d/%d cases have experiment-specific directions\n", 
                n_exp_specific, nrow(direction_analysis)))
    cat("This suggests biological variability across experiments is real.\n")
  }
  
  cat("\n💡 RECOMMENDATIONS:\n")
  cat("1. For Fisher's tests, filter by BOTH direction AND experiment\n")
  cat("2. Consider experiment as a grouping factor in analyses\n")
  cat("3. The 'inflation' may be legitimate biological signal\n")
  
} else {
  cat("🎯 CONCLUSION: No directional inflation found in CRISPRi data\n")
  cat("The issue may be in MAST data or data consolidation process\n")
}

cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("Analysis complete. Please review findings above.\n")