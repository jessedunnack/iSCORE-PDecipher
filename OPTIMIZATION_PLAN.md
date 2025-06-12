# iSCORE-PDecipher Performance Optimization Plan

## Issues Identified

### Critical Performance Problems
1. **Multiple Data Loading** - 227K records loaded 3-5 times during startup
2. **Reactive Cascade Storms** - UI updates triggering cycles of more updates 
3. **File Redundancy** - Multiple versions of same modules/configs
4. **Library Loading Bloat** - Same packages loaded multiple times
5. **Synchronous Startup** - Everything loads at once instead of progressively

### Flickering Root Cause
App ready → data loads → reactive cascade → UI updates → repeat 5-10x until stable

## Optimization Strategy

### Phase 1: Emergency Performance Fixes (IMMEDIATE)
- [ ] Implement centralized data loading with caching
- [ ] Fix reactive cascades with proper isolate() and req()
- [ ] Remove immediate execution from startup_manager.R
- [ ] Implement progressive loading (basic UI first)

### Phase 2: Code Consolidation (WEEK 1)
- [ ] Merge global.R and global_minimal.R
- [ ] Remove unused module versions
- [ ] Create shared utility functions
- [ ] Implement centralized data management

### Phase 3: Library Optimization (WEEK 2)  
- [ ] Implement lazy library loading
- [ ] Remove redundant library() calls
- [ ] Conditional loading for heavy packages
- [ ] Create dependency management system

## Files Flagged for Removal/Consolidation

### Immediate Removal
- `app_reference.R` (backup file)
- `mod_landing_page.R` (superseded)
- `mod_landing_page_with_umap.R` (superseded)
- `mod_precomputed.R` (superseded by reactive version)
- `mod_heatmap.R` (superseded by unified version)

### Consolidation Required
- `global.R` + `global_minimal.R` → single `global.R`
- Duplicate library loading across files
- Repeated data loading functions
- Redundant configuration settings

### Simplification Opportunities
- Startup sequence (too complex)
- Module sourcing chain (circular dependencies)
- Reactive expression patterns (cascading updates)
- Environment variable handling (duplicated)

## Expected Results
- **80%+ reduction** in startup time
- **Elimination** of flickering behavior  
- **50%+ reduction** in codebase size
- **Much more maintainable** architecture

## Implementation Priority
1. Fix data loading (biggest performance impact)
2. Remove file redundancy (immediate cleanup)
3. Fix reactive patterns (eliminate flickering)
4. Optimize library loading (final polish)