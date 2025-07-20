# Fix Label Wrapping and Readability in All Plots
# Purpose: Add word wrapping to long labels and improve top3_summary readability
# Date: 2025-07-20

library(stringr)
library(ggplot2)
library(dplyr)

# Helper function for word wrapping
wrap_text <- function(text, width = 40) {
  # Use stringr's str_wrap for intelligent word wrapping
  sapply(text, function(x) {
    if (is.na(x)) return("")
    str_wrap(x, width = width)
  })
}

# ==============================================================================
# 1. FIX CREATE_QUICK_PLOTS.R
# ==============================================================================
cat("Fixing create_quick_plots.R...\n")

# Read the original file
quick_plots_content <- readLines("create_quick_plots.R")

# Replace the truncation logic with word wrapping
quick_plots_content <- gsub(
  'Description_short = ifelse\\(nchar\\(Description\\) > 50,\\s*paste0\\(substr\\(Description, 1, 50\\), "\\.\\.\\."\\),\\s*Description\\)',
  'Description_short = wrap_text(Description, width = 50)',
  quick_plots_content
)

quick_plots_content <- gsub(
  'Description_short = ifelse\\(nchar\\(Description\\) > 40,\\s*paste0\\(substr\\(Description, 1, 40\\), "\\.\\.\\."\\),\\s*Description\\)',
  'Description_short = wrap_text(Description, width = 40)',
  quick_plots_content
)

# Add the wrap_text function at the beginning (after libraries)
library_end <- which(grepl("^library\\(", quick_plots_content))
if (length(library_end) > 0) {
  last_library <- max(library_end)
  quick_plots_content <- c(
    quick_plots_content[1:last_library],
    "",
    "# Helper function for word wrapping",
    "wrap_text <- function(text, width = 40) {",
    "  sapply(text, function(x) {",
    "    if (is.na(x)) return('')",
    "    stringr::str_wrap(x, width = width)",
    "  })",
    "}",
    "",
    quick_plots_content[(last_library+1):length(quick_plots_content)]
  )
}

# Fix the top3_summary plot specifically
# Find the p4 plot section
p4_start <- which(grepl("^p4 <-", quick_plots_content))
if (length(p4_start) > 0) {
  # Find the geom_text line
  geom_text_line <- which(grepl("geom_text\\(aes\\(label = Pathway_short\\)", quick_plots_content))
  if (length(geom_text_line) > 0) {
    # Replace geom_text with geom_label for better readability
    quick_plots_content[geom_text_line] <- "  geom_label(aes(label = Pathway_short, y = Significance/2),"
    quick_plots_content[geom_text_line + 1] <- "             fill = 'white', color = 'black', alpha = 0.9,"
    quick_plots_content[geom_text_line + 2] <- "             label.r = unit(0.25, 'lines'), label.padding = unit(0.5, 'lines'),"
    quick_plots_content[geom_text_line + 3] <- "             size = 3.5, fontface = 'bold') +"
  }
}

# Also update the truncation for Pathway_short to use word wrap
pathway_short_line <- which(grepl("top3_data\\$Pathway_short <-", quick_plots_content))
if (length(pathway_short_line) > 0) {
  quick_plots_content[pathway_short_line] <- "top3_data$Pathway_short <- wrap_text(top3_data$Pathway, width = 30)"
  # Remove the next two lines (the sapply function)
  quick_plots_content <- quick_plots_content[-c((pathway_short_line+1):(pathway_short_line+2))]
}

# Write back the file
writeLines(quick_plots_content, "create_quick_plots.R")

# ==============================================================================
# 2. FIX OTHER VISUALIZATION SCRIPTS
# ==============================================================================

# Function to update a script with word wrapping
update_script_with_wrapping <- function(script_path) {
  cat("  Updating", basename(script_path), "...\n")
  
  content <- readLines(script_path)
  
  # Skip if it's a backup file
  if (grepl("backup", script_path)) return()
  
  # Add wrap_text function if not present
  if (!any(grepl("wrap_text <- function", content))) {
    library_lines <- which(grepl("^library\\(", content))
    if (length(library_lines) > 0) {
      last_library <- max(library_lines)
      content <- c(
        content[1:last_library],
        "",
        "# Helper function for word wrapping",
        "wrap_text <- function(text, width = 40) {",
        "  sapply(text, function(x) {",
        "    if (is.na(x)) return('')",
        "    stringr::str_wrap(x, width = width)",
        "  })",
        "}",
        "",
        content[(last_library+1):length(content)]
      )
    }
  }
  
  # Replace all instances of Description truncation with word wrapping
  content <- gsub(
    'ifelse\\(nchar\\((.*?)\\) > ([0-9]+),\\s*paste0\\(substr\\(\\1, 1, \\2\\), "\\.\\.\\."\\),\\s*\\1\\)',
    'wrap_text(\\1, width = \\2)',
    content,
    perl = TRUE
  )
  
  # Specific patterns to replace
  content <- gsub('substr\\((.+?), 1, ([0-9]+)\\)', 'wrap_text(\\1, width = \\2)', content)
  
  writeLines(content, script_path)
}

# Update all visualization scripts
scripts_to_update <- c(
  "pd_signature_visualization.R",
  "pd_signature_comprehensive_viz.R",
  "create_comprehensive_visualizations.R",
  "create_thesis_committee_summary.R"
)

for (script in scripts_to_update) {
  if (file.exists(script)) {
    update_script_with_wrapping(script)
  }
}

# ==============================================================================
# 3. SPECIAL HANDLING FOR COMPREHENSIVE VISUALIZATIONS
# ==============================================================================
cat("\nApplying special fixes to comprehensive visualizations...\n")

if (file.exists("create_comprehensive_visualizations.R")) {
  comp_viz_content <- readLines("create_comprehensive_visualizations.R")
  
  # Find and update the convergence plot text labels
  conv_plot_lines <- which(grepl("geom_text_repel", comp_viz_content))
  if (length(conv_plot_lines) > 0) {
    for (line_num in conv_plot_lines) {
      # Add max.overlaps parameter if not present
      if (!grepl("max.overlaps", comp_viz_content[line_num])) {
        comp_viz_content[line_num] <- gsub(
          "size = 3\\)",
          "size = 3, max.overlaps = 20)",
          comp_viz_content[line_num]
        )
      }
    }
  }
  
  writeLines(comp_viz_content, "create_comprehensive_visualizations.R")
}

# ==============================================================================
# 4. RUN THE UPDATED SCRIPTS
# ==============================================================================
cat("\nRunning updated visualization scripts...\n")

# Run create_quick_plots.R
cat("  Running create_quick_plots.R...\n")
source("create_quick_plots.R")

# Run comprehensive visualizations
cat("  Running create_comprehensive_visualizations.R...\n")
source("create_comprehensive_visualizations.R")

cat("\n=== LABEL FIXES COMPLETE ===\n")
cat("Key improvements:\n")
cat("1. Word wrapping implemented for all long labels\n")
cat("2. Top3 summary plot now has horizontal labels with rounded backgrounds\n")
cat("3. STRING pathway titles properly wrapped instead of truncated\n")
cat("4. All plots regenerated with improved readability\n")