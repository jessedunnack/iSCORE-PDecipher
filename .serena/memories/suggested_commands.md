# Essential Commands for iSCORE-PDecipher Development

## Environment Setup
```bash
# Activate R environment
conda activate seuratv4  # R 4.3.3 environment
```

## Package Development Workflow
```r
# Load development tools
library(devtools)
library(testthat)

# Load package in development
devtools::load_all()           # Ctrl+Shift+L in RStudio
devtools::document()           # Update documentation
devtools::check()              # Comprehensive package check

# Install package locally
devtools::install()
```

## Testing Commands (CRITICAL - Currently Missing!)
```r
# Set up testing infrastructure (NEEDED)
usethis::use_testthat()        # Create tests/ directory structure

# Run tests (once implemented)
devtools::test()               # Ctrl+Shift+T in RStudio
testthat::test_dir("tests/testthat/")

# Check test coverage
covr::report()
```

## Git Commands
```bash
git status                     # Check repository status
git add .                      # Stage changes
git commit -m "message"        # Commit changes
git push origin main           # Push to remote
```

## Analysis Commands
```r
# Launch Shiny app
launch_iscore_app()

# Run core analysis (if exported properly)
# run_mast_analysis(mutation, seurat_path, output_dir)
# run_mixscale_analysis(experiment_path, output_dir, modality)
```

## System Utilities (Linux WSL2)
```bash
ls -la                         # List files with details
find . -name "*.R" -type f     # Find R files
grep -r "pattern" R/           # Search in R files
cd /mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PDecipher/
```

## Critical Priority Actions
1. **IMMEDIATE**: Set up formal testing with `usethis::use_testthat()`  
2. **HIGH**: Fix NAMESPACE synchronization with @export tags
3. **HIGH**: Create comprehensive tests for core analysis functions
4. **MEDIUM**: Implement golden standard test data
5. **MEDIUM**: Set up automated testing workflow