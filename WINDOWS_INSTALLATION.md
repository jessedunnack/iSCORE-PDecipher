# iSCORE-PDecipher Windows Installation Guide

This guide helps you install and use the iSCORE-PDecipher R package on Windows with RStudio.

## Prerequisites

### 1. Required Software
- **R** version 4.0.0 or higher ([Download R](https://cran.r-project.org/bin/windows/base/))
- **RStudio** ([Download RStudio](https://posit.co/download/rstudio-desktop/))
- **Rtools** (for building packages from source) ([Download Rtools](https://cran.r-project.org/bin/windows/Rtools/))

### 2. Data Location
Ensure your iSCORE datasets are accessible on Windows. Common locations:
- `E:/ASAP/scRNASeq/PerturbSeq/final/`
- `C:/ASAP/scRNASeq/PerturbSeq/final/`
- `D:/ASAP/scRNASeq/PerturbSeq/final/`

Your dataset directory should contain subdirectories like:
- `iSCORE-PD/`
- `iSCORE-PD_plus_CRISPRi/`

## Installation

### Method 1: Install from GitHub (Recommended)

1. **Open RStudio** on your Windows machine

2. **Install devtools** (if not already installed):
```r
install.packages("devtools")
```

3. **Install required Bioconductor packages**:
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "clusterProfiler", "ReactomePA", "DOSE", "org.Hs.eg.db", 
    "pathview", "enrichplot", "ComplexHeatmap", "SingleCellExperiment"
))
```

4. **Install iSCORE-PDecipher from GitHub**:
```r
devtools::install_github("jessedunnack/iSCORE-PDecipher")
```

### Method 2: Install from Local Files

If you have the package files locally:

1. **Navigate to the package directory** in R:
```r
# Adjust path to your local copy
setwd("E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PDecipher")
```

2. **Install dependencies**:
```r
devtools::install_deps()
```

3. **Install the package**:
```r
devtools::install()
```

## First-Time Setup

### 1. Load the Package
```r
library(iSCORE.PDecipher)
```

### 2. Configure Data Location
On first launch, the package will prompt you to configure your data directory:

```r
# This will open a setup dialog
launch_app()
```

**Or configure manually**:
```r
# Set your parent data directory (one-time setup)
set_parent_data_dir("E:/ASAP/scRNASeq/PerturbSeq/final")

# Verify configuration
get_parent_data_dir()
```

## Launching the Application

### Quick Launch (Recommended)
```r
# Simple launch with interactive dataset selection
launch_app()
```

### Launch with Specific Dataset
```r
# For mutations only
launch_app("E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD")

# For mutations + CRISPRi
launch_app("E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi")
```

### Advanced Launch Options
```r
# Custom port and browser settings
launch_app(port = 8080, launch.browser = TRUE)

# Full function name (equivalent to launch_app)
launch_iscore_app()
```

## Dataset Validation

The package will automatically check for required files:
- `full_DE_results.rds` - Differential expression results
- `all_enrichment_padj005_complete_with_direction.rds` - Enrichment data
- `enrichment_results/` directory - Individual enrichment files

If files are missing, the app will offer to generate them (requires source data).

## Troubleshooting

### Common Issues

#### 1. Package Installation Errors
```r
# If you encounter compilation errors, try installing binaries
install.packages("devtools", type = "binary")
```

#### 2. Path Issues
```r
# Check your current configuration
get_parent_data_dir()

# Reset configuration if needed
set_parent_data_dir("E:/your/actual/path")
```

#### 3. Missing Dependencies
```r
# Check for missing packages
check_required_packages()

# Install specific missing packages
install.packages("package_name")
```

#### 4. Data Files Not Found
```r
# Check what datasets are available
get_dataset_options()

# Validate a specific directory
validate_dataset_directory("E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD")
```

### Memory Issues
For large datasets (200K+ cells), ensure you have:
- At least 16GB RAM
- Close other memory-intensive applications

### Performance Tips
- Use static UMAP visualizations for large datasets
- Enable caching in the app settings
- Consider using smaller sample percentages for initial exploration

## Complete Example Session

```r
# 1. Load package
library(iSCORE.PDecipher)

# 2. Set up data location (one-time)
set_parent_data_dir("E:/ASAP/scRNASeq/PerturbSeq/final")

# 3. Check available datasets
get_dataset_options()

# 4. Launch with interactive selection
launch_app()

# The app should now open in your default web browser
```

## Support

If you encounter issues:

1. **Check the console** for error messages in RStudio
2. **Verify file paths** using Windows Explorer
3. **Ensure adequate memory** (16GB+ recommended)
4. **Report issues** at: https://github.com/jessedunnack/iSCORE-PDecipher/issues

## Version Information

- Package version: 0.4.0
- R version required: ≥ 4.0.0
- Platform: Windows 10/11 (tested)
- RStudio: Any recent version

---

*Last updated: August 2025*