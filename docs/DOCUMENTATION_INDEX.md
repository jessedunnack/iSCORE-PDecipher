# iSCORE-PDecipher Documentation Index

## 📚 Complete Documentation Directory

### 🚀 User Guides
- **[Main README](../README.md)** - Primary package documentation with installation and quick start
- **[LABMATE_QUICKSTART](LABMATE_QUICKSTART.md)** - Step-by-step guide for new users
- **[DATASET_SHARING_PLAN](DATASET_SHARING_PLAN.md)** - Instructions for sharing datasets with collaborators

### 🖥️ Cross-Platform Setup  
- **[MAC_COMPATIBILITY_GUIDE](MAC_COMPATIBILITY_GUIDE.md)** - Complete Mac setup instructions
- **[MAC_SETUP_CHECKLIST](MAC_SETUP_CHECKLIST.md)** - Quick Mac setup reference
- **[CROSS_PLATFORM_DEVELOPMENT_GUIDELINES](CROSS_PLATFORM_DEVELOPMENT_GUIDELINES.md)** - Developer guidelines for cross-platform compatibility

### 📊 Analysis Guides
- **[PD_SIGNATURE_ANALYSIS_GUIDE](PD_SIGNATURE_ANALYSIS_GUIDE.md)** - Guide to signature nomination and PD biology analysis
- **[GSEA_VISUALIZATION_SUMMARY](GSEA_VISUALIZATION_SUMMARY.md)** - GSEA analysis and visualization guide

### 🔬 Technical Documentation
- **[STATISTICAL_OVERLAP_ANALYSIS_REPORT](STATISTICAL_OVERLAP_ANALYSIS_REPORT.md)** - Statistical methods and validation
- **[SHINY_APP_TESTING_CHECKLIST](SHINY_APP_TESTING_CHECKLIST.md)** - Quality assurance testing procedures

### 📈 Development & Status
- **[SIGNATURE_NOMINATION_STATUS](SIGNATURE_NOMINATION_STATUS.md)** - Current status of signature analysis features
- **[REORGANIZATION_STATUS_LOG](REORGANIZATION_STATUS_LOG.md)** - App reorganization progress
- **[APP_REORGANIZATION_SUMMARY](APP_REORGANIZATION_SUMMARY.md)** - Summary of app structure changes
- **[FUTURE_ENHANCEMENTS](FUTURE_ENHANCEMENTS.md)** - Planned features and improvements

### 🛠️ Development Plans
- **[BIDIRECTIONAL_SYNC_IMPLEMENTATION_PLAN](BIDIRECTIONAL_SYNC_IMPLEMENTATION_PLAN.md)** - Data synchronization planning
- **[DE_GENE_EXPANSION_PLAN](DE_GENE_EXPANSION_PLAN.md)** - Differential expression analysis expansion
- **[CHAT_CONTEXT_SUMMARY](CHAT_CONTEXT_SUMMARY.md)** - Development conversation summaries

### 🎯 Specialized Features
- **[CHECKPOINT_v0.1.4](CHECKPOINT_v0.1.4.md)** - Version 0.1.4 feature checkpoint
- **[GOOGLE_DRIVE_SHARING_CHECKLIST](GOOGLE_DRIVE_SHARING_CHECKLIST.md)** - Google Drive integration checklist

## 🔧 Development Files

### Development Tools (dev/)
- **Launch utilities**: `INSTALL_LATEST_VERSION.R`, `quick_launch.R`
- **Setup guides**: `LAUNCH_INSTRUCTIONS.md`

### Debug & Testing (dev/debug/)
- **Debug scripts**: Various `debug_*.R` files for troubleshooting
- **Test scripts**: Various `test_*.R` files for validation
- **Analysis scripts**: `analyze_current_results.R`, `explore_*.R`

### Development Testing (dev/)
- **Function tests**: `test_app_functions.R`, `test_heatmap_generation.R`
- **Module tests**: `test_reactive_flow.R`, `test_reliability.R`
- **Integration tests**: `test_launch.R`, `test_de_results.R`

## 📂 Package Structure Overview

```
iSCORE-PDecipher/
├── 📄 Core Package Files
│   ├── DESCRIPTION          # Package metadata and dependencies
│   ├── NAMESPACE           # Exported functions
│   └── README.md           # Main documentation
│
├── 🔧 R Functions (R/)
│   ├── Core analysis functions
│   ├── Signature discovery functions  
│   ├── Data import/export functions
│   └── Configuration management
│
├── 📱 Shiny Application (inst/shiny/)
│   ├── app.R               # Main application
│   ├── modules/           # Modular UI components
│   ├── R/                 # App-specific functions
│   └── www/              # Web assets (CSS, JS)
│
├── 📊 Data & Scripts (inst/)
│   ├── extdata/           # Package data files
│   ├── scripts/           # Utility scripts
│   └── examples/          # Usage examples
│
├── 📚 Documentation (docs/)
│   ├── User guides and tutorials
│   ├── Technical documentation
│   └── Development status files
│
├── 🛠️ Development (dev/)
│   ├── debug/             # Debug scripts
│   ├── Installation utilities
│   └── Testing tools
│
└── 📝 Manual Pages (man/)
    └── Function documentation
```

## 🎯 Quick Navigation

**For Users**: Start with [Main README](../README.md) → [LABMATE_QUICKSTART](LABMATE_QUICKSTART.md)

**For Mac Users**: [MAC_COMPATIBILITY_GUIDE](MAC_COMPATIBILITY_GUIDE.md) → [MAC_SETUP_CHECKLIST](MAC_SETUP_CHECKLIST.md)

**For Developers**: [CROSS_PLATFORM_DEVELOPMENT_GUIDELINES](CROSS_PLATFORM_DEVELOPMENT_GUIDELINES.md) → [Development files in dev/](../dev/)

**For Analysis**: [PD_SIGNATURE_ANALYSIS_GUIDE](PD_SIGNATURE_ANALYSIS_GUIDE.md) → [Signature Nomination features](../README.md#signature-nomination-module--new-v020)

---

*Last Updated: January 2025 (v0.2.0)*