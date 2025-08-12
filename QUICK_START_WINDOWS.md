# 🚀 Quick Start for Windows RStudio Users

## Step 1: Install from GitHub in RStudio

Open RStudio and run:

```r
# Install the package from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("jessedunnack/iSCORE-PDecipher")
```

## Step 2: Load the Package

```r
library(iSCORE.PDecipher)
```

## Step 3: Configure Data Location (First Time Only)

The package will automatically look for your data in:
- `E:/ASAP/scRNASeq/PerturbSeq/final/`
- `C:/ASAP/scRNASeq/PerturbSeq/final/`
- `D:/ASAP/scRNASeq/PerturbSeq/final/`

If your data is elsewhere, set it:

```r
# Only needed if data is in a different location
set_parent_data_dir("E:/ASAP/scRNASeq/PerturbSeq/final")
```

## Step 4: Launch the App

```r
# Simple launch - will prompt for dataset selection
launch_app()
```

### What You'll See:

1. **Console Output**: Shows available datasets
   - ✅ iSCORE-PD (Mutations Only) - 210K cells
   - ✅ iSCORE-PD + CRISPRi - 232K cells

2. **Dataset Selection**: Choose 1 or 2 in the console

3. **App Opens**: The Shiny app launches in your browser with:
   - Purple banner showing current dataset
   - "Change Dataset" button to switch anytime
   - All analysis tabs ready to use

## Alternative: Direct Launch with Dataset Pre-Selected

```r
# Launch with iSCORE-PD only
launch_app(dataset = "iSCORE-PD")

# Launch with iSCORE-PD + CRISPRi
launch_app(dataset = "iSCORE-PD_plus_CRISPRi")
```

## Troubleshooting

### If the app doesn't launch:

1. **Check your working directory**:
```r
getwd()
# Should show your project directory
```

2. **Verify data files exist**:
```r
check_data_availability()
```

3. **Manual path setting**:
```r
# If auto-detection fails
options(iscore.data_dir = "E:/ASAP/scRNASeq/PerturbSeq/final")
launch_app()
```

### If you get package errors:

```r
# Reinstall with dependencies
devtools::install_github("jessedunnack/iSCORE-PDecipher", 
                        dependencies = TRUE,
                        force = TRUE)
```

## Features Available

Once the app is running, you can:

- **Switch datasets** without restarting (Change Dataset button)
- **View DE results** with cached UMAP plots (40-200x faster)
- **Explore enrichment** across all pathways
- **Generate heatmaps** with interactive controls
- **Use preview mode** for faster loading with large datasets

## Best Workflow for RStudio

1. Open RStudio
2. Set working directory: `setwd("E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PDecipher")`
3. Run: `library(iSCORE.PDecipher); launch_app()`
4. Select dataset when prompted
5. App opens in browser
6. Use "Change Dataset" button to switch between datasets

That's it! The app handles all the path conversions and data loading automatically.