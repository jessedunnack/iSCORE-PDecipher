# iSCORE-PDecipher Folder Structure
**Last Updated:** January 20, 2025 | **Version:** 0.3.5

This document describes the organized folder structure of the iSCORE-PDecipher package.

## Root Directory
The root directory contains only essential package files:
- `DESCRIPTION` - Package metadata, dependencies, and version information
- `NAMESPACE` - Package exports and imports
- `NEWS.md` - Detailed changelog of all package versions
- `README.md` - Main package documentation and quick start guide

## Directory Structure

### `/R`
Core R package functions, organized by functionality:
- `launch_app.R` - Main app launcher function
- `*_analysis.R` - Analysis functions (MAST, MixScale, enrichment, etc.)
- `config_manager.R` - Cross-platform configuration management
- `signature_*.R` - Signature discovery and visualization functions
- `data_*.R` - Data import and processing functions
- Other utility and helper functions

### `/inst`
Files to be installed with the package:

#### `/inst/shiny`
Complete Shiny application:
- `app.R` - Main application file
- `global.R` - Global settings and libraries
- `/modules` - Shiny modules for different features
- `/R` - Shiny-specific helper functions
- `/www` - Web assets (CSS, JavaScript)

#### `/inst/analysis_scripts`
Standalone analysis scripts for manuscript preparation:
- PD signature discovery scripts
- Visualization generation scripts
- Comprehensive report templates

#### `/inst/extdata`
External data files:
- Gene association databases
- UMAP data files for visualization
- Reference datasets

#### `/inst/results`
Analysis results and outputs:
- Fisher test results
- Correlation analysis results
- PD biological relevance data

#### `/inst/scripts`
Utility scripts for data processing:
- UMAP data extraction
- Cluster marker calculation
- Correlation analysis tools

### `/man`
R documentation files (auto-generated)

### `/dev`
Development and testing files:
- `/debug` - Debugging scripts
- `/tests` - Test scripts and validation
- `/temp` - Temporary working files
- Quick launch scripts and utilities

### `/docs`
Comprehensive documentation:
- Analysis guides and methodologies
- Bug fix documentation
- Cross-platform setup guides
- Statistical method improvements
- Archive of older documentation

### `/scripts`
Main analysis scripts:
- PD signature analysis scripts
- Visualization creation scripts
- Data processing utilities

### `/results`
Analysis outputs (not tracked in git):
- `/pd_signatures` - PD signature discovery results
- Visualization outputs
- Summary reports

### `/vignettes`
Package vignettes and tutorials (if any)

## File Organization Principles

1. **Root Cleanliness**: Only essential package files in root
2. **Logical Grouping**: Files grouped by function and purpose
3. **Clear Naming**: Descriptive file names indicating content
4. **Documentation**: Each major directory has its own README
5. **Version Control**: Large data files excluded via .gitignore

## Quick Navigation

- **To run the app**: See `/inst/shiny/app.R`
- **For analysis scripts**: Check `/inst/analysis_scripts/`
- **For documentation**: Browse `/docs/`
- **For development**: Use `/dev/`
- **For results**: Look in `/results/` (local only)

## Recent Changes (v0.3.5)

- Moved all documentation files from root to `/docs/`
- Moved test scripts from root to `/dev/tests/`
- Moved analysis scripts from root to `/scripts/`
- Moved data files from root to `/inst/results/`
- Created this folder structure documentation