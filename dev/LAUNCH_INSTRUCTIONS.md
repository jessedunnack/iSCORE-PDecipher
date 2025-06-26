# How to Launch iSCORE-PDecipher App

## Issue: Function Not Found
You're getting `Error: could not find function "launch_iscore_app"` because:
1. The R version is too old (3.6.3, but package requires R >= 4.0.0)
2. Dependencies are missing or incompatible

## Quick Solutions

### Option 1: Launch Directly (RECOMMENDED)
```r
# In R, navigate to the package directory first:
setwd("/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher")

# Then run the quick launch script:
source("quick_launch.R")
```

### Option 2: Manual Function Loading
```r
# Load functions directly
source("R/launch_app.R")
source("R/config_manager.R")

# Then launch
launch_iscore_app()
```

### Option 3: GitHub Installation (if you have R >= 4.0.0)
```r
# Update the package from GitHub
remotes::install_github("jessedunnack/iSCORE-PDecipher", force = TRUE)
library(iSCORE.PDecipher)
launch_iscore_app()
```

## What to Expect
When you run `launch_iscore_app()`, you should see:
1. **First-time setup**: Prompt to select parent data directory
2. **Dataset selection**: Choose from available datasets
3. **App launch**: Shiny app opens in browser

## Testing the Signature Nomination Fixes
Once the app is running:
1. Navigate to **Signature Nomination** tab
2. Click **"Discover Signatures"** button  
3. Check these tabs:
   - ✅ **Pan-Cluster Signatures**: Should show table with gene pairs
   - ✅ **PD Biology Focus**: Should show biological analysis (not empty)
   - ⚠️ **Cluster-Specific Signatures**: Check if dropdown populates
   - ⚠️ **Signature Heatmap**: Check if it loads or spins infinitely

## Next Steps After Testing
Please report back:
- [ ] Does the app launch successfully?
- [ ] Are the signature tabs now showing data?
- [ ] Which specific features still need work?

The core PD analysis bug is fixed, but some UI elements may need additional work depending on what you see.