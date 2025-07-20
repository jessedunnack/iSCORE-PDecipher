# iSCORE.PDecipher 0.3.3

## Enhanced Visualizations and Documentation (2025-07-20)

### 🎨 **Visualization Improvements**
- **Fixed**: Label visibility in convergence strength plots
- **Corrected**: Pathway count display to show actual totals instead of filtered top 30
- **Updated**: Summary documentation with corrected pathway counts and key insights

---

# iSCORE.PDecipher 0.3.2

## Improved Plot Readability (2025-07-20)

### 📊 **Visualization Enhancements**
- **Implemented**: Word wrapping for long pathway names in all PD signature plots
- **Improved**: Label readability across bar charts and visualizations
- **Enhanced**: Overall plot aesthetics for publication quality

---

# iSCORE.PDecipher 0.3.1

## Comprehensive Fixes Applied In-Place (2025-07-20)

This patch release applies critical fixes to all visualization outputs without creating separate directories.

### 🔧 **Fixes Applied**

#### Gene Harmonization
- **Fixed**: Gene naming mismatches between MAST and CRISPRi datasets
- **Mappings**: PRKN→PARK2, VPS13C variants→VPS13C, SNCA variants→SNCA
- **Impact**: PARK2 now shows 106 convergent pathways (was 0), proper cross-method analysis

#### Heatmap Clustering
- **Fixed**: All heatmaps now use hierarchical clustering
- **Applied to**: Gene pathway heatmaps, cluster heatmaps
- **Benefit**: Patterns and relationships more clearly visible

#### Natural Cluster Sorting
- **Fixed**: Clusters now sort numerically (cluster_10 after cluster_9, not cluster_1)
- **Applied to**: All cluster visualizations and data files
- **Order**: cluster_0 through cluster_14 in proper sequence

### 📁 **Important Note**
All fixes applied directly to original file locations:
- `/results/pd_signatures/by_gene/`
- `/results/pd_signatures/by_cluster/`
- `/results/pd_signatures/visualizations/comprehensive/`

No "_fixed" directories created - all updates in place.

---

# iSCORE.PDecipher 0.3.0

## Comprehensive PD Signature Analysis Suite (2025-07-20)

This major release introduces a comprehensive analysis pipeline for identifying compelling Parkinson's Disease signatures across mutation and CRISPRi perturbation methods, preparing for manuscript publication.

### 🎯 **New Analysis Features**

#### PD Signature Discovery Pipeline
- **Fast discovery script**: Efficiently filters 663K+ enrichments to identify PD-relevant pathways
- **Method-specific analysis**: Identifies mutation-only, CRISPRi-only, and convergent signatures
- **Statistical rigor**: All pathways filtered with p.adjust < 0.05

#### Gene-by-Gene Analysis
- **Individual gene reports**: Analyzes all 16 PD genes (including SNCA variants)
- **Pattern classification**: Identifies mutation-dominant, CRISPRi-dominant, and convergent genes
- **Visual summaries**: Creates gene-specific signature plots and heatmaps

#### Cluster Analysis Framework
- **Comprehensive coverage**: Analyzes all 15 clusters (cluster_0 through cluster_14)
- **Specificity analysis**: Identifies cluster-specific vs ubiquitous pathways
- **Method preferences**: Determines which clusters prefer mutation vs CRISPRi signatures

#### Advanced Visualizations
- **Multi-panel figures**: Comprehensive summary visualizations for manuscripts
- **Gene × pathway matrices**: Shows pathway distribution across all genes
- **Cluster heatmaps**: Visualizes enrichment patterns across clusters
- **Updated labeling**: "Mutation - iSCORE-PD" and "CRISPRi Perturbation" throughout

### 📊 **Key Discoveries**

- **30 mutation-only pathways**: Including vesicle-mediated transport in synapse
- **30 CRISPRi-only pathways**: Including mitochondrial small ribosomal subunit
- **30 convergent pathways**: Led by synapse pathway (p < 10^-17)
- **Gene patterns**: SNCA variants show mutation-dominance, ATP13A2 shows convergence
- **Cluster insights**: Cluster_4 has highest enrichment (1,152 pathways)

### 🔧 **Technical Improvements**

#### Enhanced Heatmap Module (from concurrent development)
- **CRISPRi experiment separation**: Individual columns for C12_FPD-24, C12_FPD-23, C18_FPD-23
- **Column annotations**: Visual method distinction (MAST blue, CRISPRi orange)
- **Advanced UI**: Collapsible sections, gene filters, term search, p-value threshold
- **Biological categorization**: Transparent keyword-based pathway classification

#### Analysis Scripts Organization
- **Location**: `inst/analysis_scripts/` for package distribution
- **Documentation**: Comprehensive README with usage instructions
- **Reproducibility**: All paths and parameters documented

### 📝 **Reports and Documentation**

- **Comprehensive R Markdown report**: Full analysis with interactive elements
- **Executive summary**: Key findings for thesis committee presentation
- **Visualization gallery**: Publication-ready figures in multiple formats

### 🐛 **Bug Fixes** (from 0.2.9 branch)
- Fixed gene details modal to show actual genes
- Fixed CRISPRi dotplots with modality column fallback
- Restricted DE table search to Gene column only
- Maintained CRISPRa exclusion from analyses

---

# iSCORE.PDecipher 0.2.9

## Minimal Bug Fixes from Clean Baseline (2025-07-17)

This release applies minimal, targeted fixes to address 9 bugs while avoiding the introduction of new issues. Started from baseline commit d6c08c0 and applied only necessary changes.

### 🐛 **Bug Fixes Applied**

#### Bug #1: Gene Details Modal
- Added Details button to gene pair analysis table
- Implemented modal showing actual DE genes and enrichment terms
- Fixed early returns that prevented gene list extraction

#### Bug #3: Cross-Cluster DE Gene Analysis 
- Added new Cross-Cluster Analysis tab
- Implemented proper column handling to prevent [object Object] errors
- Shows genes DE across multiple clusters with statistics

#### Bug #5: CRISPRi Dotplot Hover Tooltips
- Fixed gene lookup for CRISPRi data - now uses actual experiment ID
- Added conservative fallback for missing data with transparency
- Hover tooltips now show actual DE genes instead of "No genes found"

#### Bug #6: Interactive DE Gene Overlap Heatmap
- Added DE Gene Overlap Heatmap tab with interactive visualization
- Implemented using plotly/heatmaply for full interactivity
- Shows Fisher's exact test p-values for gene overlap significance

#### Bug #8: CRISPRa Exclusion
- Ensured CRISPRa data is always excluded from analyses
- Modified detect_available_methods to always return FALSE for CRISPRa

#### Bug #9: Gene Dropdown Cleanup
- Clean gene names without "(X variants)" text (avoided by using baseline)

### 🔧 **Technical Notes**
- Based on clean baseline (d6c08c0) to avoid tangential issues
- Each fix applied independently and minimally
- Conservative approach to prevent introduction of new bugs
- See BUG_FIX_COMPREHENSIVE_ANALYSIS.md for detailed documentation

---

# iSCORE.PDecipher 0.2.7

## Major Correlation Plot Visualization Improvements (2025-01-15)

### 🎯 **Enhanced Interactive Correlation Plots**

#### Visual Reference Lines and Statistics Display
- **Added bold reference lines**: x=0 and y=0 lines (black, 0.8 linewidth, 0.7 alpha) for easy identification of positive/negative effects
- **Correlation statistics in titles**: Both cluster-specific and pan-cluster plots now display:
  - Pearson correlation coefficient (r = 0.73)
  - Statistical significance (p = 0.002 or p < 0.001)
  - Sample size (n = 53 genes)
- **Format**: "Cluster 3\nr = 0.81, p < 0.001 (n = 45)" for easy interpretation

#### Revolutionary Layout Change: Grid to Vertical Stacking
- **Problem Solved**: Compressed grid layout made text unreadable
- **New Approach**: Vertical stacking of cluster-specific plots
- **Benefits**:
  - Much larger individual plots with readable axis labels
  - Clear correlation statistics and gene names
  - Natural scrolling through clusters (0, 1, 2, ...)
  - Independent axis scaling per cluster
  - Better focus on individual cluster correlations

#### Technical Implementation Details
- **Dynamic height**: 400px × number of clusters for optimal viewing
- **Independent scaling**: `shareX = FALSE, shareY = FALSE` for proper axis ranges
- **Enhanced spacing**: 8% margin between plots for clear separation
- **User guidance**: "Scroll to view all clusters" in main title
- **UI updates**: Changed from "Cluster grid view" to "Cluster-specific view (vertical)"

### 🔧 **Critical Bug Fixes**
- **Fixed "condition has length > 1" error**: Replaced vectorized `if()` with `ifelse()` in hover text generation
- **Resolved plot rendering failures**: Proper vectorized operations for gene highlighting

### 📊 **Gene Selection Methodology Confirmed**
**Question**: Are genes significant in both approaches?
**Answer**: Yes, but by effect size rather than p-value significance:
- **MAST filtering**: Top N genes by absolute log2FC (highest effect size changes)
- **CRISPRi filtering**: Top N genes by absolute log2FC (highest effect size changes)
- **Final correlation**: Only genes present in BOTH filtered datasets (via merge by gene name)
- **Default**: Top 200 most changed genes per method
- **Options**: Top 100, 200, 500, 1000, or all genes

### 🎨 **User Experience Improvements**
- **Increased plot height**: From 450px to 700px for better spacing
- **Reference lines**: Easy identification of origin (0,0) and effect quadrants
- **Scrollable design**: Natural interaction pattern for cluster exploration
- **Statistical transparency**: All correlation statistics visible at a glance

### 📝 **Code Structure Updates**
- **File**: `inst/shiny/modules/mod_signature_nomination.R`
- **Functions enhanced**: `renderPlotly` correlation plot generation
- **Layout logic**: Replaced `plotly::subplot()` grid with vertical stacking
- **Statistical calculation**: Added real-time correlation computation per cluster

---

# iSCORE.PDecipher 0.2.6

## Enhanced Cross-Method Correlation Analysis (2025-01-15)

### New Features

#### Gene Filtering for Correlation Analysis
- **Improvement**: Gene filtering increases correlation strength by 11x (mean |r| = 0.593 vs 0.053)
- **Methods**: Top 100/200/500 genes vs all genes comparison
- **UI**: Interactive dropdown in Gene Pair Analysis with explanatory text
- **Validation**: Tested across 180 combinations (12 gene pairs × 3 experiments × 5 clusters)

#### Interactive Correlation Plots
- **Visualization**: Interactive plotly plots with gene hover information
- **Features**: Trend lines, experiment comparison, gene selection details
- **Default**: C12_FPD-24 experiment set as default
- **Location**: Signature Nomination module, Gene Pair Analysis tab

#### Analysis Scripts
- **Scripts**: Comprehensive correlation analysis pipeline
- **Location**: `inst/scripts/correlation_analysis/`
- **Results**: `inst/results/comprehensive_correlation_results.csv`

#### Statistical Framework Enhancements
- **Direction-aware analysis**: Same vs opposite direction overlap detection
- **Hierarchical FDR correction**: Multi-level multiple comparison correction
- **Experiment weighting**: Cell count-based weighting system

### Strong Correlations Identified
- **DNAJC6**: r = 0.99 (strongest correlation)
- **VPS13C variants**: Consistent strong correlations
- **SNCA variants**: Strong directional concordance
- **Total**: 61 gene pairs with |r| ≥ 0.5

### Documentation
- **Guide**: `docs/correlation_analysis_guide.md`
- **Methodology**: Complete description of gene filtering approach
- **Results**: Comprehensive validation results

---

# iSCORE.PDecipher 0.2.3

## Critical Volcano Plot Fix (2025-01-13)

This hotfix release resolves a critical discrepancy in volcano plot data display for gene variants.

### 🔥 **CRITICAL BUG FIX: VPS13C and SNCA Volcano Plot Gene Harmonization**
- **Issue**: VPS13C_W395C volcano plots showed "no results available" while summary statistics showed overlap data
- **Root Cause**: Missing gene harmonization in volcano plot rendering vs summary statistics
- **Fix Applied**: Added gene harmonization to both static and interactive MixScale volcano plots
- **Impact**: All gene variants now show consistent data between volcano plots and statistics

### 🧬 **Gene Variants Fixed**
- **VPS13C_W395C** → Now correctly maps to VPS13C MixScale data
- **VPS13C_A444P** → Now correctly maps to VPS13C MixScale data  
- **SNCA_A30P** → Now correctly maps to SNCA MixScale data
- **SNCA_A53T** → Now correctly maps to SNCA MixScale data
- **PRKN** → Confirmed working mapping to PARK2 MixScale data

### ✅ **Testing & Validation**
- Added comprehensive test suite (`test_vps13c_volcano_fix.R`)
- Verified all 11 clusters have data for each variant
- Confirmed gene harmonization logic working correctly
- All target genes available in MixScale dataset

### 🔧 **Technical Implementation**
- Modified `inst/shiny/modules/mod_de_results.R` volcano plot renderers
- Applied same harmonization logic used in Fisher's exact test calculations
- Ensures consistency across all DE Results page visualizations

**Resolution**: Volcano plots now display correctly for all gene variants, matching the behavior of summary statistics.

---

# iSCORE.PDecipher 0.2.2

## Major Statistical Improvements (2025-01-13)

This release introduces critical statistical improvements to cross-method gene overlap analysis and comprehensive code organization.

### 🔬 Fisher's Exact Test Enhancements
- **Fixed**: Fundamental statistical flaw in Fisher's exact test implementation
- **Issue**: Background gene universe used union of significant genes instead of all tested genes
- **Solution**: Dual approach implementation with proper background extraction
  - **Intersection approach** (conservative): Genes tested in BOTH methods (~15-18K genes)
  - **Union approach** (liberal): Genes tested in EITHER method (~20-22K genes)
- **Impact**: Statistically valid cross-method comparisons between MAST and CRISPRi

### 📊 New Analysis Functions
- **Added**: `run_comprehensive_fishers_analysis()` function for systematic overlap testing
- **Features**: 
  - Tests all gene-cluster combinations automatically
  - Provides both intersection and union statistics
  - Exports detailed results and gene-level summaries
  - Includes gene harmonization (PRKN→PARK2, SNCA variants→SNCA)

### 🧹 Code Organization
- **Created**: `.gitignore` file to prevent tracking of large data files
- **Organized**: Moved test scripts to `dev/tests/`
- **Organized**: Moved temporary working scripts to `dev/temp/`
- **Cleaned**: Removed temporary gene association files (*_temp_*.rds)
- **Refactored**: Converted analysis scripts to proper package functions

### 🐛 Additional Bug Fixes
- **Fixed**: Stale data bug in DE Results "Overlapping DE Genes" section
- **Fixed**: Infinite reactive loop when changing global mutation settings
- **Fixed**: Gene harmonization now properly applied in DE Results module
- **Fixed**: P-value formatting for better readability (decimals for p ≥ 0.001)

## Technical Details
- Background genes now properly extracted from current selection (not all comparisons)
- Fisher's exact test results displayed with both approaches in UI
- Comprehensive documentation in `FISHER_EXACT_TEST_IMPROVEMENTS.md`
- Example analysis of 130 gene-cluster combinations completed successfully

---

# iSCORE.PDecipher 0.2.1

## Critical Fixes (2025-01-13)

This release addresses three critical issues that prevented core functionality:

### 🧬 Gene Association Loading Fix (HIGH-PRIORITY)
- **Fixed**: Gene display feature completely non-functional due to locked binding error
- **Issue**: `cannot change value of locked binding for '.gene_associations'`
- **Solution**: Implemented environment-based storage replacing global variables
- **Impact**: Gene lists now visible in enrichment tables and hover tooltips (24,000 associations)

### 🔧 DE Results Module Namespace Fix  
- **Fixed**: DE Results page crashes with "ns function not found" error
- **Issue**: `DT::dataTableOutput(ns("cluster_markers_table"))` namespace error at line 771
- **Solution**: Added explicit namespace capture in renderUI functions
- **Impact**: DE Results page now loads correctly without crashes

### ⚡ DE Heatmap Performance Optimization
- **Fixed**: DE Heatmap generation stalling/hanging indefinitely  
- **Issue**: O(13 × 14 × 200K) = 36M operations with inefficient nested loops
- **Solution**: Performance optimization with progress indicators
  - Pre-filtering by significance before processing
  - Vectorized operations replace nested loops  
  - Real-time progress feedback for users
  - Configurable parameters for different dataset sizes
- **Impact**: Fast heatmap generation with user feedback, no more perceived "hanging"

## Technical Improvements

### Performance
- 52% faster processing on small datasets
- Progress indicators eliminate user confusion during long operations
- Memory usage optimized with pre-allocated storage
- Vectorized filtering operations

### User Experience  
- Gene lists visible in Plot Details tables
- Interactive hover tooltips show associated genes
- All pages load without crashes
- Clear progress feedback during operations

### Development
- Environment-based storage eliminates package namespace conflicts
- Backwards compatibility maintained for all existing functions
- Comprehensive test suite added for critical functionality
- Enhanced error handling and fallback mechanisms

## Testing & Validation

- ✅ Live Shiny app testing: All modules load without errors
- ✅ User acceptance testing: 3/3 scenarios passed  
- ✅ Performance benchmarking: Significant improvements verified
- ✅ Integration testing: All critical fixes working together

## Deployment
- All fixes tested and validated in live environment
- Rollback procedures documented with checkpoint tags
- Ready for immediate production deployment

---

# iSCORE.PDecipher 0.2.0

## Previous Release Notes
[Previous changelog content would go here if it existed]