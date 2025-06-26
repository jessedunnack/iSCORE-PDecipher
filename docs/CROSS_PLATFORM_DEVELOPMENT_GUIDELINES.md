# Cross-Platform Development Guidelines for iSCORE-PDecipher

## 🎯 Overview

This document provides essential guidelines for maintaining cross-platform compatibility in the iSCORE-PDecipher R package. **All future development must follow these principles** to ensure the package works seamlessly across Windows, Mac, and Linux systems.

## 🚨 Critical Rules - NEVER Break These

### 1. **Path Management Rules**

#### ✅ **DO - Use Cross-Platform Path Functions**
```r
# Correct - cross-platform path construction
file.path("data", "results", "analysis.rds")
file.path(base_dir, "enrichment_results", paste0(gene, "_", cluster, ".rds"))

# Correct - normalize paths when storing
config$parent_dir <- normalizePath(user_path, mustWork = FALSE)

# Correct - expand home directory 
expanded_path <- path.expand("~/Documents/data")
```

#### ❌ **DON'T - Hard-code Path Separators or Drive Letters**
```r
# WRONG - will break on Mac/Linux
file_path <- "E:/ASAP/data/results/analysis.rds"
relative_path <- "data\\results\\analysis.rds"  # Windows-only separators

# WRONG - platform-specific assumptions
if (file.exists("C:/Program Files/R/")) { ... }
```

### 2. **Configuration Storage Rules**

#### ✅ **DO - Use Platform-Appropriate Config Directories**
```r
# Correct - uses rappdirs for cross-platform config
config_dir <- rappdirs::user_config_dir("iSCORE.PDecipher", "jessedunnack")
config_file <- file.path(config_dir, "config.json")

# Correct - create directories if needed
if (!dir.exists(config_dir)) {
  dir.create(config_dir, recursive = TRUE)
}
```

#### ❌ **DON'T - Assume Writable Locations**
```r
# WRONG - package directory may not be writable
config_file <- file.path(system.file(package = "iSCORE.PDecipher"), "config.json")

# WRONG - working directory assumptions
config_file <- "config.json"  # Current directory may not be writable
```

### 3. **Package Function Export Rules**

#### ✅ **DO - Proper NAMESPACE Exports**
```r
# Correct - export all functions used by other package functions
#' @export
my_function <- function() { ... }

# Correct - update NAMESPACE file
export(my_function)
```

#### ❌ **DON'T - Source Files in Installed Packages**
```r
# WRONG - R files don't exist in installed packages
source(file.path(pkg_root, "R", "helper_functions.R"))

# WRONG - assuming development directory structure
if (file.exists("R/config.R")) { source("R/config.R") }
```

### 4. **User Interface Rules**

#### ✅ **DO - Provide Cross-Platform File Selection**
```r
# Correct - cross-platform file selection with fallback
select_directory <- function() {
  if (interactive() && requireNamespace("tcltk", quietly = TRUE)) {
    tryCatch({
      tcltk::tk_choose.dir(caption = "Select directory")
    }, error = function(e) {
      # Fallback to manual entry
      readline("Enter directory path: ")
    })
  } else {
    readline("Enter directory path: ")
  }
}
```

#### ❌ **DON'T - Assume Specific File Dialogs**
```r
# WRONG - Windows-specific file dialog
path <- choose.dir()  # May not work on all platforms

# WRONG - assuming GUI availability
file.choose()  # May fail in non-interactive environments
```

## 📋 Development Checklist

Before adding any new code, verify:

### **Path Handling Checklist**
- [ ] All paths constructed with `file.path()`
- [ ] No hard-coded drive letters (E:/, C:/, etc.)
- [ ] No hard-coded path separators (`\`, `/`)
- [ ] Paths normalized with `normalizePath()` when stored
- [ ] Home directory expansion supported with `path.expand()`

### **Configuration Checklist**
- [ ] Config stored in platform-appropriate location (`rappdirs`)
- [ ] Config format is cross-platform (JSON, not R-specific)
- [ ] Graceful handling of missing/corrupted config files
- [ ] No assumptions about writable working directory

### **Function Export Checklist**
- [ ] All public functions have `@export` documentation
- [ ] NAMESPACE file updated with new exports
- [ ] No source() calls for package-internal files
- [ ] Package installation tested with `remotes::install_github()`

### **User Experience Checklist**
- [ ] Clear error messages for platform-specific issues
- [ ] Fallback options for platform-specific features
- [ ] First-time setup guidance provided
- [ ] Documentation covers platform differences

## 🧪 Testing Protocol

### **Required Testing Before Any Commit**

1. **Test Package Installation**
```r
# Test fresh installation
remove.packages("iSCORE.PDecipher")
remotes::install_github("jessedunnack/iSCORE-PDecipher", force = TRUE)
library(iSCORE.PDecipher)
```

2. **Test Cross-Platform Functions**
```r
# Test config system
is_first_launch()
get_config_path()

# Test path handling
generate_transfer_report("/test/path")
```

3. **Test App Launch**
```r
# Test app launch without errors
launch_iscore_app()
```

### **Platform-Specific Testing**

#### **Windows Testing**
- Test with both forward and backslash paths
- Verify drive letter handling
- Test with spaces in paths

#### **Mac Testing**
- Test with `~/` home directory paths
- Verify case-sensitive filesystem handling
- Test with special characters in paths

#### **Linux Testing**
- Test with various mount points
- Verify symbolic link handling
- Test with permission-restricted directories

## 🔧 Common Platform Issues & Solutions

### **Issue 1: Path Construction**
```r
# Problem: Platform-specific separators
bad_path <- paste("data", "files", "result.rds", sep = "/")  # Fails on Windows

# Solution: Use file.path()
good_path <- file.path("data", "files", "result.rds")  # Works everywhere
```

### **Issue 2: Config Storage**
```r
# Problem: Non-portable config location  
bad_config <- "~/.myapp/config.json"  # May not work on Windows

# Solution: Use rappdirs
good_config <- file.path(
  rappdirs::user_config_dir("myapp", "author"), 
  "config.json"
)
```

### **Issue 3: File Existence Checks**
```r
# Problem: Case sensitivity assumptions
if (file.exists("Data/Results.RDS")) { ... }  # May fail on case-sensitive systems

# Solution: Store exact filenames or use case-insensitive logic
stored_files <- list.files(dirname(path), full.names = TRUE)
target_file <- stored_files[tolower(basename(stored_files)) == tolower("Results.RDS")]
```

### **Issue 4: Package Dependencies**
```r
# Problem: Source file dependencies
source("R/helper.R")  # Fails in installed packages

# Solution: Proper exports
#' @export
helper_function <- function() { ... }
# And add to NAMESPACE: export(helper_function)
```

## 📚 Reference Implementation

### **Config Manager Pattern**
See `R/config_manager.R` for complete cross-platform config implementation:
- Platform-appropriate directory detection
- JSON serialization for portability
- Graceful error handling
- First-launch detection

### **Path Validation Pattern**
See `R/prepare_mac_transfer.R` for path validation examples:
- Cross-platform directory checking
- Smart fallback directory detection
- User-friendly error messages

### **Function Export Pattern**
See `NAMESPACE` for proper export declarations:
- All user-facing functions exported
- Internal helper functions not exported
- No source() dependencies

## 🚨 Emergency Troubleshooting

### **If Package Breaks on Mac/Linux**

1. **Check for Hard-coded Paths**
```bash
grep -r "E:/" R/
grep -r "C:/" R/
grep -r "\\\\" R/
```

2. **Check for Missing Exports**
```r
# Test function availability after installation
library(iSCORE.PDecipher)
exists("problematic_function")
```

3. **Check for Source Dependencies**
```bash
grep -r "source(" R/
```

### **If Config System Fails**

1. **Check Config Directory Creation**
```r
config_dir <- rappdirs::user_config_dir("iSCORE.PDecipher", "jessedunnack")
dir.exists(config_dir)
file.access(config_dir, mode = 2)  # Check write permissions
```

2. **Reset User Config**
```r
# For testing - remove user config
unlink(get_config_path())
is_first_launch()  # Should return TRUE
```

## 🔄 Maintenance Tasks

### **Regular Reviews**
- **Monthly**: Review new code for platform assumptions
- **Before Releases**: Full cross-platform testing
- **After Dependency Updates**: Verify compatibility maintained

### **Documentation Updates**
- Keep this guide updated with new patterns
- Update user guides when adding platform-specific features
- Document any known platform limitations

## 🏆 Success Metrics

The package is considered cross-platform compatible when:

- ✅ **Installation**: Works via `remotes::install_github()` on all platforms
- ✅ **First Launch**: Succeeds without errors on fresh systems
- ✅ **Config Persistence**: Settings survive R restarts and package updates
- ✅ **Data Transfer**: Files can be moved between platforms and work
- ✅ **Error Handling**: Clear messages for platform-specific issues

## 📞 Support & Contact

**For Cross-Platform Issues:**
- Check this guide first
- Test with minimal reproducible example
- Document specific platform and R version
- Include error messages and config file contents

**Remember**: Cross-platform compatibility is not optional - it's a core requirement for the iSCORE-PDecipher package. Every line of code should work on Windows, Mac, and Linux systems.