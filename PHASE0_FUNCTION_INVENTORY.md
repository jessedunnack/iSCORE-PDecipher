# iSCORE-PDecipher Function Inventory
## Phase 0: Discovery & Context Gathering

**Generated:** 2025-11-08
**Purpose:** Complete inventory of all functions, modules, and components in the iSCORE-PDecipher R package
**Status:** Comprehensive baseline for documentation modernization project

---

## Executive Summary

### Component Statistics

| Category | Count | Description |
|----------|-------|-------------|
| **R Package Functions** | 192 | Core package functions across 35 files |
| **Shiny Helper Functions** | 45 | Helper functions for Shiny app (9 files) |
| **Shiny Module UI Functions** | 22 | Module UI constructors (23 modules) |
| **Shiny Module Server Functions** | 21 | Module server logic |
| **Reactive Expressions** | 37 | Reactive data computations |
| **Observers (observe)** | 29 | Reactive side-effect handlers |
| **Observers (observeEvent)** | 47 | Event-driven reactive handlers |
| **Render Functions** | 143 | UI rendering functions |
| **TOTAL COMPONENTS** | **536** | **Complete package ecosystem** |

### File Distribution

- **R Package Files:** 35 files (R/ directory)
- **Shiny Modules:** 23 files (inst/shiny/modules/)
- **Shiny Helpers:** 9 files (inst/shiny/R/)
- **Main App Files:** 3 files (app.R, global.R, app_perturbseq_full.R)

---

## 1. R Package Functions (192 functions across 35 files)

### 1.1 Core Analysis Functions

#### MAST_analysis.R (4.9K, 143 lines)
- **run_mast_analysis()** (Line 11) - @export - Run MAST differential expression analysis

#### MAST_analysis_optimized.R (15K, 371 lines)
- **run_mast_analysis_optimized()** (Line 17) - @export - Optimized MAST analysis with performance improvements
- **validate_optimized_mast_results()** (Line 314) - @return Logical indicating if results are equivalent

#### MixScale_analysis.R (21K, 440 lines)
- **run_mixscale_analysis()** (Line 11) - @export - Run MixScale perturbation analysis
- **process_mixscale()** (Line 74) - Internal MixScale data processing

### 1.2 Data Import & Validation Functions

#### data_import_functions.R (14K, 335 lines)
- **extract_cluster_id()** (Line 8) - @export - Extract cluster ID from file paths
- **import_mast_data()** (Line 30) - @export - Import MAST differential expression results
- **import_mixscale_data()** (Line 108) - @export - Import MixScale perturbation results

#### data_import_functions_optimized.R (19K, 539 lines)
- **extract_cluster_id_fast()** (Line 7) - @return String: the extracted cluster ID
- **import_mast_data_optimized()** (Line 35) - @return List of structured MAST results
- **import_mixscale_data_optimized()** (Line 241) - @return List of structured MixScale results
- **load_lazy_data()** (Line 455) - @return The actual data loaded from file
- **validate_optimized_import()** (Line 492) - @return Logical indicating validation success

#### import_functions_enhanced.R (16K, 503 lines)
- **import_enrichment_fdr()** (Line 16) - @export - Import FDR-corrected enrichment results
- **import_convergence_results()** (Line 98) - @export - Import convergence analysis results
- **import_gene_list()** (Line 162) - @export - Import gene lists
- **import_enrichment_complete()** (Line 213) - @export - Import complete enrichment dataset
- **get_data_inventory()** (Line 263) - @export - Get inventory of available data
- **compare_approaches()** (Line 415) - @export - Compare different analysis approaches
- **create_data_accessor()** (Line 465) - @export - Create data accessor functions

#### import_pooled_mixscale_functions.R (15K, 447 lines)
- **detect_mixscale_format()** (Line 12) - @export - Detect MixScale data format (pooled vs experiment-split)
- **import_pooled_mixscale_data()** (Line 74) - Import pooled MixScale data with FDR corrections
- **extract_cluster_id()** (Line 230) - @export - Extract cluster ID from pooled data
- **import_enrichment_with_correction()** (Line 277) - Import enrichment with specific p-value correction

#### dataset_validator.R (6.9K, 203 lines)
- **check_source_data()** (Line 8) - @return List with status and messages
- **check_missing_files()** (Line 74) - @return List of missing file types
- **validate_dataset_directory()** (Line 101) - @export - Validate dataset directory structure
- **get_dataset_options()** (Line 155) - @export - Get available dataset options

### 1.3 Enrichment Analysis Functions

#### enrichment_analysis.R (80K, 2080 lines)
**[Largest R file in package - comprehensive enrichment pipeline]**

- **print_header()** (Line 10) - Print formatted header messages
- **run_enrichment_analysis()** (Line 35) - @return Invisible summary of the analysis
- **process_mast_entry()** (Line 203) - @param padj_method Method for p-value adjustment
- **process_mixscale_entry()** (Line 343) - @param padj_method Method for p-value adjustment
- **run_all_enrichment_analyses()** (Line 896) - @return List of enrichment results for each analysis type
- **filter_genes()** (Line 1263) - @return List of up, down, and all DE genes
- **filter_mixscale_genes()** (Line 1293) - @return List containing by_experiment results
- **get_weighted_pvalue_cols()** (Line 1445) - @return Named vector of p-value column names
- **run_go_enrichment()** (Line 1499) - @return GO enrichment results
- **run_kegg_enrichment()** (Line 1546) - @return KEGG enrichment results
- **run_reactome_enrichment()** (Line 1622) - @return Reactome enrichment results
- **run_wikipathways_enrichment()** (Line 1698) - @return WikiPathways enrichment results
- **run_string_ppi()** (Line 1773) - @return STRING PPI network results
- **run_gsea()** (Line 1849) - @return GSEA results by category

#### process_enrichment_results.R (12K, 358 lines)
- **extract_enrichment_data()** (Line 31) - @return A data frame with standardized columns
- **extract_direction()** (Line 101) - @return Character string: "UP", "DOWN", "ALL", or "RANKED"
- **process_single_file()** (Line 125) - @return Data frame with enrichment results and metadata
- **process_enrichment_results()** (Line 209) - @return Data frame with all enrichment results

#### term_extraction_functions.R (14K, 455 lines)
- **handle_enrichment_result()** (Line 9) - @return A data frame of enrichment results with standardized columns
- **extract_string_terms()** (Line 91) - @return Data frame with standardized columns
- **extract_gsea_terms()** (Line 132) - @return Data frame with standardized columns
- **load_and_extract_terms()** (Line 178) - @return Data frame with enrichment terms
- **get_enrichment_file_path()** (Line 214) - @return Full path to the enrichment result file
- **compare_terms()** (Line 250) - @return Data frame with comparison results
- **find_frequent_terms()** (Line 319) - @return Data frame with frequent terms and their statistics
- **prepare_heatmap_data()** (Line 395) - @return Matrix suitable for heatmap visualization

### 1.4 Signature Analysis Functions

#### signature_analysis.R (45K, 1120 lines)
**[Second largest R file - comprehensive signature discovery]**

- **calculate_gene_overlap_significance_proper()** (Line 15) - @export - Calculate proper gene overlap significance
- **calculate_fisher_test()** (Line 169) - @return List with Fisher's test results
- **calculate_gene_overlap_significance()** (Line 214) - @export - Calculate gene overlap significance
- **calculate_effect_size_correlation()** (Line 292) - @export - Calculate effect size correlation
- **calculate_direction_consistency()** (Line 372) - @export - Calculate direction consistency
- **calculate_composite_signature_score()** (Line 426) - @export - Calculate composite signature score
- **identify_pd_relevant_enrichments()** (Line 509) - @export - Identify PD-relevant enrichments
- **calculate_pathway_overlap_by_database()** (Line 548) - @export - Calculate pathway overlap by database
- **analyze_gene_pair_signatures()** (Line 609) - @export - Analyze gene pair signatures
- **calculate_enhanced_overlap_significance()** (Line 857) - @export - Calculate enhanced overlap significance
- **combine_weighted_results()** (Line 1043) - @return List with combined weighted results
- **combine_meta_pvalues()** (Line 1090) - @return Combined p-value

#### manuscript_signature_discovery.R (34K, 771 lines)
- **discover_top_signatures()** (Line 20) - @export - Discover top signatures for manuscript
- **apply_hierarchical_fdr_correction()** (Line 142) - @return Data frame with FDR-corrected p-values added
- **apply_enhanced_fdr_correction_v026()** (Line 295) - @export - Apply enhanced FDR correction v0.26
- **compute_signature_rankings()** (Line 457) - @return Data frame with ranked signatures
- **identify_pan_cluster_signatures()** (Line 621) - @return Data frame with pan-cluster signatures
- **identify_cluster_specific_signatures()** (Line 656) - @return List of cluster-specific signatures by cluster
- **generate_manuscript_signature_summary()** (Line 686) - @export - Generate manuscript signature summary

#### pd_signature_interpretation.R (28K, 728 lines)
- **analyze_pd_signatures()** (Line 13) - @export - Analyze PD signatures
- **extract_signature_biological_context()** (Line 70) - @return List with biological context information
- **extract_signature_genes()** (Line 151) - @return List with mast_gene and crispri_gene
- **get_signature_enrichment_terms()** (Line 184) - @return Filtered enrichment terms
- **filter_pd_relevant_terms()** (Line 244) - @return PD-relevant terms with relevance scores
- **find_shared_pathways()** (Line 292) - @return Data frame of shared pathway information
- **find_conceptually_similar_pathways()** (Line 337) - @return Data frame of conceptually similar pathways
- **categorize_biological_processes()** (Line 376) - @return List of biological process categories
- **generate_signature_interpretation()** (Line 417) - @return Character string with biological interpretation
- **generate_category_interpretation()** (Line 503) - @return Interpretation text for the category
- **generate_biological_significance()** (Line 523) - @return Biological significance text
- **calculate_pd_relevance_score()** (Line 571) - @return Numeric PD relevance score
- **create_pd_signature_summary()** (Line 599) - @return Comprehensive summary for manuscript
- **generate_overall_biological_insights()** (Line 652) - @return Biological insights text
- **generate_cross_signature_interpretation()** (Line 688) - @return Cross-signature interpretation text
- **str_to_title()** (Line 724) - String utility function

#### signature_trends_analysis.R (17K, 440 lines)
- **analyze_signature_trends()** (Line 14) - @export - Analyze signature trends across clusters
- **get_signature_strength()** (Line 86) - @return Vector of signature strength values
- **validate_trends_inputs()** (Line 103) - @return List with validation status and message
- **compute_signature_frequency_analysis()** (Line 151) - @return List with frequency analysis results
- **compute_signature_impact_analysis()** (Line 227) - @return List with impact analysis results
- **discover_enrichment_term_patterns()** (Line 306) - @return List with term pattern analysis
- **classify_term_pattern()** (Line 350) - @return Pattern category
- **create_empty_term_patterns()** (Line 373) - @return Empty term patterns structure
- **create_trends_summary()** (Line 387) - @return Summary statistics
- **create_empty_trends_result()** (Line 412) - @return Empty trends result structure

#### signature_data_helpers.R (3.8K, 125 lines)
- **get_signature_strength()** (Line 12) - @export - Get signature strength
- **get_cluster_info()** (Line 30) - @export - Get cluster information
- **get_signature_metric()** (Line 47) - @export - Get signature metric
- **validate_signature_data()** (Line 69) - @export - Validate signature data
- **safe_max_signature_strength()** (Line 87) - @export - Safe maximum signature strength
- **safe_signature_access()** (Line 101) - @export - Safe signature data access

#### signature_results_converter.R (5.4K, 126 lines)
- **convert_signature_results_for_trends()** (Line 9) - @export - Convert signature results for trends analysis

#### signature_visualization_functions.R (22K, 529 lines)
- **create_gene_pathway_pvalue_scatter()** (Line 17) - @export - Create gene-pathway p-value scatter plot
- **create_interactive_signature_heatmap()** (Line 112) - @export - Create interactive signature heatmap
- **create_interactive_signature_heatmap_enhanced()** (Line 192) - @export - Create enhanced interactive heatmap
- **create_gene_pair_multi_metric_dashboard()** (Line 411) - @export - Create multi-metric dashboard
- **create_pathway_category_bubble_chart()** (Line 475) - @export - Create pathway category bubble chart

### 1.5 Statistical Analysis Functions

#### comprehensive_fishers_analysis.R (16K, 383 lines)
- **run_comprehensive_fishers_analysis()** (Line 24) - Run comprehensive Fisher's exact test analysis
- **format_pvalue()** (Line 288) - @noRd - Format p-values for display
- **process_mast_for_volcano()** (Line 299) - @noRd - Process MAST data for volcano plots
- **process_mixscale_for_volcano()** (Line 330) - @noRd - Process MixScale data for volcano plots
- **apply_gene_harmonization()** (Line 371) - @noRd - Apply gene harmonization

#### enhanced_direction_analysis.R (22K, 528 lines)
- **get_biological_direction_expectation()** (Line 17) - @export - Get biological direction expectation
- **calculate_same_direction_overlap()** (Line 71) - @export - Calculate same direction overlap
- **calculate_opposite_direction_overlap()** (Line 214) - @export - Calculate opposite direction overlap
- **combine_direction_pvalues()** (Line 352) - @export - Combine direction p-values
- **enhanced_direction_analysis()** (Line 447) - @export - Run enhanced direction analysis

#### experiment_weighting.R (16K, 402 lines)
- **load_crispri_cell_counts()** (Line 16) - @export - Load CRISPRi cell counts
- **extract_seurat_cell_counts()** (Line 81) - @export - Extract Seurat cell counts
- **extract_experiment_from_assignment()** (Line 141) - @return Vector of experiment IDs
- **calculate_experiment_weights()** (Line 161) - @export - Calculate experiment weights
- **add_estimated_weights()** (Line 251) - @return Updated weights list
- **normalize_weights_by_cluster()** (Line 290) - @return List of normalized weights
- **identify_primary_experiment()** (Line 319) - @return Name of primary experiment
- **weighted_meta_analysis()** (Line 340) - @export - Perform weighted meta-analysis

### 1.6 Gene Processing Functions

#### gene_harmonization.R (12K, 331 lines)
- **create_gene_mapping_table()** (Line 10) - @export - Create gene mapping table
- **get_comparable_gene_pairs()** (Line 65) - @export - Get comparable gene pairs
- **get_mutation_categories()** (Line 130) - @export - Get mutation categories
- **filter_for_gene_comparison()** (Line 246) - @export - Filter for gene comparison
- **get_pd_relevant_pathways()** (Line 303) - @export - Get PD-relevant pathways

#### extract_gene_associations.R (8.6K, 275 lines)
- **extract_from_single_file()** (Line 19) - @return Data frame with gene associations
- **parse_file_metadata()** (Line 100) - @return List with metadata components
- **extract_all_gene_associations()** (Line 162) - @return Data frame with all gene associations
- **main()** (Line 247) - @export - Main extraction function

#### gene_association_lookup.R (7.5K, 233 lines)
- **load_gene_associations()** (Line 16) - @export - Load gene associations
- **get_gene_associations()** (Line 61) - @export - Get gene associations for specific term
- **gene_associations_available()** (Line 73) - @export - Check if gene associations are available
- **get_genes_for_term()** (Line 89) - @export - Get genes for specific term
- **get_genes_for_terms()** (Line 148) - @export - Get genes for multiple terms
- **search_gene_associations()** (Line 172) - @export - Search gene associations
- **get_association_stats()** (Line 212) - @export - Get association statistics

### 1.7 Data Management Functions

#### data_sampling.R (11K, 330 lines)
- **sample_seurat_cells()** (Line 19) - @export - Sample cells from Seurat object
- **create_preview_dataset()** (Line 109) - @export - Create preview dataset
- **extract_umap_data()** (Line 165) - @export - Extract UMAP data
- **create_progressive_umap()** (Line 210) - @export - Create progressive UMAP
- **create_reactive_dataset()** (Line 243) - @export - Create reactive dataset
- **estimate_memory_usage()** (Line 297) - @export - Estimate memory usage

#### data_generator.R (7.4K, 219 lines)
- **generate_missing_files()** (Line 10) - @export - Generate missing data files
- **check_required_packages()** (Line 191) - @export - Check required packages

### 1.8 Configuration & App Launch Functions

#### config_manager.R (6.4K, 241 lines)
- **get_config_path()** (Line 11) - @export - Get configuration file path
- **load_config()** (Line 26) - @export - Load configuration
- **save_config()** (Line 47) - @export - Save configuration
- **get_parent_data_dir()** (Line 63) - @export - Get parent data directory
- **set_parent_data_dir()** (Line 72) - @export - Set parent data directory
- **is_first_launch()** (Line 82) - @export - Check if first launch
- **prompt_for_parent_dir()** (Line 91) - @export - Prompt for parent directory
- **validate_parent_dir()** (Line 143) - @export - Validate parent directory
- **setup_parent_dir()** (Line 201) - @export - Setup parent directory

#### launch_app.R (7.5K, 224 lines)
- **launch_iscore_app()** (Line 24) - Launch iSCORE-PDecipher Shiny application
- **select_dataset_directory()** (Line 138) - @return Selected directory path or NULL if cancelled
- **launch_app()** (Line 223) - Simplified launcher (alias for launch_iscore_app)

#### run_app.R (3.4K, 115 lines)
- **run_app()** (Line 31) - Alternative app launcher

#### run_full_analysis.R (4.3K, 107 lines)
- **run_complete_pipeline()** (Line 27) - Run complete analysis pipeline

### 1.9 Performance & Benchmarking Functions

#### performance_benchmarking.R (16K, 470 lines)
- **run_comprehensive_benchmark()** (Line 15) - @return List containing comprehensive benchmark results
- **benchmark_mast_optimizations()** (Line 128) - Benchmark MAST analysis optimizations
- **benchmark_import_optimizations()** (Line 232) - Benchmark data import optimizations
- **analyze_scalability_improvements()** (Line 354) - Analyze scalability improvements
- **get_system_memory_info()** (Line 414) - Helper: Get system memory info
- **get_memory_usage()** (Line 422) - Helper: Get memory usage
- **analyze_memory_improvements()** (Line 430) - Analyze memory improvements
- **validate_functional_equivalence()** (Line 438) - Validate functional equivalence
- **analyze_scaling_trends()** (Line 446) - Analyze scaling trends
- **generate_benchmark_summary()** (Line 454) - Generate benchmark summary
- **write_benchmark_report()** (Line 463) - Write benchmark report

### 1.10 Agent Coordination Functions

#### agent_coordinator.R (8.7K, 320 lines)
- **invoke_agent()** (Line 12) - @export - Generic agent invocation
- **invoke_performance_agent()** (Line 29) - @keywords internal - Invoke performance agent
- **invoke_maintenance_agent()** (Line 71) - @keywords internal - Invoke maintenance agent
- **invoke_github_agent()** (Line 107) - @keywords internal - Invoke GitHub agent
- **invoke_data_agent()** (Line 140) - @keywords internal - Invoke data agent
- **invoke_testing_agent()** (Line 174) - @keywords internal - Invoke testing agent
- **invoke_refactoring_agent()** (Line 211) - @keywords internal - Invoke refactoring agent
- **run_maintenance_cycle()** (Line 247) - @export - Run maintenance cycle
- **generate_agent_report()** (Line 287) - @export - Generate agent report

### 1.11 Mac Transfer Utilities

#### prepare_mac_transfer.R (9.7K, 294 lines)
- **check_dataset_files()** (Line 12) - @export - Check dataset files for Mac transfer
- **generate_transfer_report()** (Line 109) - @export - Generate transfer report
- **prepare_mac_transfer_copy()** (Line 221) - @export - Prepare Mac transfer copy

### 1.12 UMAP Processing Functions

#### process_dataset2_umap_wrapper.R (4.2K, 141 lines)
- **process_dataset2_umap()** (Line 42) - Process Dataset2 UMAP
- **process_dataset2_umap_full()** (Line 127) - Process Dataset2 UMAP (full version)

### 1.13 Package Documentation

#### iSCORE-PDecipher-package.R (2.8K, 70 lines)
**[Package-level documentation file - no exported functions]**

---

## 2. Shiny Helper Functions (45 functions across 9 files)

### 2.1 Visualization Helpers

#### bubble_heatmap_functions.R (6.5K, 206 lines)
- **create_bubble_heatmap()** (Line 23) - Create bubble heatmap visualization

#### heatmap_functions.R (8.9K, 313 lines)
- **prepare_enrichment_heatmap()** (Line 10) - Prepare data for enrichment heatmap
- **create_interactive_heatmap()** (Line 80) - Create interactive heatmap
- **create_static_heatmap()** (Line 115) - Create static heatmap
- **create_modality_comparison_heatmap()** (Line 188) - Create modality comparison heatmap
- **create_gene_enrichment_heatmap()** (Line 259) - Create gene enrichment heatmap

#### visualization_tiers.R (13K, 390 lines)
- **create_tier1_visualizations()** (Line 16) - Create Tier 1 visualizations
- **create_summary_table()** (Line 92) - Create summary table
- **create_tier2_visualizations()** (Line 113) - Create Tier 2 visualizations
- **create_direction_comparison_heatmap()** (Line 167) - Create direction comparison heatmap
- **create_tier3_visualizations()** (Line 234) - Create Tier 3 visualizations
- **create_simple_bar_plot()** (Line 294) - Create simple bar plot
- **create_simple_dot_plot()** (Line 339) - Create simple dot plot

### 2.2 Data Management Helpers

#### data_manager.R (11K, 335 lines)
- **init_data_manager()** (Line 9) - Initialize data manager
- **get_enrichment_data()** (Line 23) - Get enrichment data
- **get_loading_status()** (Line 83) - Get loading status
- **get_available_genes()** (Line 89) - Get available genes
- **natural_sort_clusters()** (Line 101) - Natural sort clusters
- **get_available_clusters()** (Line 116) - Get available clusters
- **get_data_summary()** (Line 124) - Get data summary
- **load_pooled_mixscale_functions()** (Line 152) - Load pooled MixScale functions
- **get_pooled_mixscale_data()** (Line 176) - Get pooled MixScale data
- **get_pooled_enrichment_data()** (Line 257) - Get pooled enrichment data
- **set_pval_correction()** (Line 298) - Set p-value correction method
- **get_pval_correction()** (Line 309) - Get p-value correction method
- **set_dataset_type()** (Line 318) - Set dataset type
- **get_dataset_type()** (Line 328) - Get dataset type

#### cache_manager.R (7.9K, 283 lines)
- **create_app_cache()** (Line 183) - Create app cache
- **load_enrichment_with_app_cache()** (Line 200) - Load enrichment with app cache
- **preload_next_files()** (Line 238) - Preload next files

### 2.3 Dataset Management Helpers

#### dataset_modal_functions.R (9.0K, 256 lines)
- **create_dataset_selector_modal()** (Line 8) - Create dataset selector modal
- **create_custom_upload_modal()** (Line 98) - Create custom upload modal
- **validate_dataset()** (Line 174) - Validate dataset
- **get_default_datasets()** (Line 219) - Get default datasets

#### dataset_selector.R (6.4K, 218 lines)
- **datasetSelectorUI()** (Line 7) - Dataset selector UI
- **datasetSelectorServer()** (Line 40) - Dataset selector server

### 2.4 Startup & Configuration Helpers

#### startup_config.R (6.3K, 195 lines)
- **detect_cross_platform_path()** (Line 6) - Detect cross-platform path
- **create_startup_config()** (Line 54) - Create startup configuration
- **load_dataset()** (Line 150) - Load dataset
- **get_current_dataset()** (Line 182) - Get current dataset

#### startup_manager.R (7.3K, 254 lines)
- **initialize_app_data()** (Line 17) - Initialize app data
- **process_uploaded_file()** (Line 61) - Process uploaded file
- **is_data_loaded()** (Line 132) - Check if data is loaded
- **create_file_picker_modal()** (Line 139) - Create file picker modal
- **initialize_app_with_data()** (Line 171) - Initialize app with data

---

## 3. Shiny Modules (23 modules with 43 UI/Server functions)

### 3.1 Core Data Modules

#### mod_data_loader.R (8.5K, 307 lines)
- **UI Functions:** 2 (mod_data_loader_ui, quick_data_switcher_ui)
- **Server Functions:** 2 (mod_data_loader_server, quick_data_switcher_server)
- **Reactive Expressions:** 2
- **Observers:** 3 (observe: 1, observeEvent: 2)
- **Render Functions:** 3 (renderUI: 3x)

#### mod_landing_page_with_umap_v2.R (45K, 1208 lines)
- **UI Functions:** 0 (embedded in main app)
- **Server Functions:** 0 (embedded in main app)
- **Reactive Expressions:** 0
- **Observers:** 7 (observe: 4, observeEvent: 3)
- **Render Functions:** 16
  - renderPlotly (3x)
  - renderUI (6x)
  - renderDataTable (2x)
  - renderPlot (3x)
  - renderText (2x)

### 3.2 Differential Expression Modules

#### mod_de_results.R (126K, 3108 lines)
**[LARGEST SHINY MODULE - comprehensive DE analysis interface]**

- **UI Function:** mod_de_results_ui() - Line 279
- **Server Function:** mod_de_results_server() - Line 552
- **Reactive Expressions:** 5
- **Observers:** 9 (observe: 4, observeEvent: 5)
- **Render Functions:** 19
  - renderUI (6x)
  - renderPlot (4x)
  - renderPlotly (3x)
  - renderDataTable (3x)
  - renderText (3x)

#### mod_de_analysis.R (20K, 588 lines)
- **UI Function:** mod_de_analysis_ui() - Line 24
- **Server Function:** mod_de_analysis_server() - Line 211
- **Reactive Expressions:** 3
- **Observers:** 2 (observe: 1, observeEvent: 1)
- **Render Functions:** 5
  - renderDataTable (2x)
  - renderPlot (2x)
  - renderPrint (1x)

#### mod_de_results_cached.R (11K, 357 lines)
- **UI Functions:** 0
- **Server Functions:** 0
- **Reactive Expressions:** 1
- **Observers:** 1 (observeEvent: 1)
- **Render Functions:** 2

### 3.3 Heatmap Modules

#### mod_heatmap.R (69K, 1723 lines)
**[SECOND LARGEST MODULE - advanced heatmap functionality]**

- **UI Function:** mod_heatmap_ui() - Line varies
- **Server Function:** mod_heatmap_server() - Line varies
- **Reactive Expressions:** 4
- **Observers:** 8 (observe: 1, observeEvent: 7)
- **Render Functions:** 5
  - renderPlot (multiple)
  - renderPlotly (multiple)
  - renderDataTable (1x)

#### mod_de_heatmap.R (31K, 864 lines)
- **UI Function:** mod_de_heatmap_ui() - Line 199
- **Server Function:** mod_de_heatmap_server() - Line 428
- **Reactive Expressions:** 3
- **Observers:** 5 (observe: 1, observeEvent: 4)
- **Render Functions:** 6
  - renderPrint (3x)
  - renderPlot (1x)
  - renderPlotly (1x)
  - renderDataTable (1x)

#### mod_heatmap_unified.R (27K, 749 lines)
- **UI Function:** mod_heatmap_unified_ui()
- **Server Function:** mod_heatmap_unified_server()
- **Reactive Expressions:** 3
- **Observers:** 5 (observeEvent: 5)
- **Render Functions:** 4

### 3.4 Signature Analysis Modules

#### mod_signature_nomination.R (149K, 3213 lines)
**[LARGEST MODULE - comprehensive signature nomination interface]**

- **UI Functions:** 2
- **Server Function:** 1
- **Reactive Expressions:** 2
- **Observers:** 6 (observe: 1, observeEvent: 5)
- **Render Functions:** 16
  - renderPlot (multiple)
  - renderPlotly (multiple)
  - renderDataTable (multiple)
  - renderUI (multiple)

#### mod_signature_trends.R (27K, 654 lines)
- **UI Function:** mod_signature_trends_ui()
- **Server Function:** mod_signature_trends_server()
- **Reactive Expressions:** 4
- **Observers:** 3 (observeEvent: 3)
- **Render Functions:** 10

### 3.5 Enrichment Display Modules

#### mod_enrichment_gene_display.R (6.1K, 187 lines)
- **UI Function:** mod_enrichment_gene_display_ui()
- **Server Function:** mod_enrichment_gene_display_server()
- **Reactive Expressions:** 0
- **Observers:** 2 (observe: 1, observeEvent: 1)
- **Render Functions:** 1 (renderDataTable)

#### mod_enrichment_gene_display_v2.R (5.6K, 178 lines)
- **UI Function:** mod_enrichment_gene_display_v2_ui()
- **Server Function:** mod_enrichment_gene_display_v2_server()
- **Reactive Expressions:** 1
- **Observers:** 2 (observe: 1, observeEvent: 1)
- **Render Functions:** 1 (renderDataTable)

### 3.6 Pathway Visualization Modules

#### mod_pathview.R (30K, 779 lines)
- **UI Function:** mod_pathview_ui()
- **Server Function:** mod_pathview_server()
- **Reactive Expressions:** 1
- **Observers:** 4 (observe: 2, observeEvent: 2)
- **Render Functions:** 3
  - renderPlot (1x)
  - renderUI (1x)
  - renderText (1x)

### 3.7 UMAP Viewer Modules

#### mod_umap_viewer.R (8.1K, 232 lines)
- **UI Function:** mod_umap_viewer_ui()
- **Server Function:** mod_umap_viewer_server()
- **Reactive Expressions:** 1
- **Observers:** 2 (observeEvent: 2)
- **Render Functions:** 3
  - renderPlotly (2x)
  - renderUI (1x)

#### mod_umap_viewer_simple.R (3.1K, 90 lines)
- **UI Function:** mod_umap_viewer_simple_ui()
- **Server Function:** mod_umap_viewer_simple_server()
- **Reactive Expressions:** 0
- **Observers:** 0
- **Render Functions:** 1 (renderPlotly)

#### umap_cache_integration.R (6.8K, 232 lines)
- **Functions:** UMAP cache integration utilities
- **Reactive Expressions:** 1
- **Observers:** 2
- **Render Functions:** 2

### 3.8 Comparison & Export Modules

#### mod_comparison.R (18K, 491 lines)
- **UI Function:** mod_comparison_ui() - Line 4
- **Server Function:** mod_comparison_server() - Line 81
- **Reactive Expressions:** 0
- **Observers:** 3 (observe: 2, observeEvent: 1)
- **Render Functions:** 8
  - renderUI (3x)
  - renderDataTable (3x)
  - renderPlot (2x)

#### mod_export.R (25K, 684 lines)
- **UI Function:** mod_export_ui()
- **Server Function:** mod_export_server()
- **Reactive Expressions:** 0
- **Observers:** 3 (observe: 2, observeEvent: 1)
- **Render Functions:** 4
  - renderDataTable (2x)
  - renderUI (1x)
  - renderText (1x)

### 3.9 Specialized Analysis Modules

#### mod_perturbseq_only.R (9.9K, 356 lines)
- **UI Function:** mod_perturbseq_only_ui()
- **Server Function:** mod_perturbseq_only_server()
- **Reactive Expressions:** 0
- **Observers:** 3 (observe: 2, observeEvent: 1)
- **Render Functions:** 7
  - renderDataTable (3x)
  - renderPlot (2x)
  - renderUI (2x)

#### mod_precomputed_reactive.R (25K, 687 lines)
- **UI Function:** mod_precomputed_reactive_ui()
- **Server Function:** mod_precomputed_reactive_server()
- **Reactive Expressions:** 4
- **Observers:** 3 (observeEvent: 3)
- **Render Functions:** 9

#### mod_visualization.R (34K, 928 lines)
- **UI Function:** mod_visualization_ui()
- **Server Function:** mod_visualization_server()
- **Reactive Expressions:** 2
- **Observers:** 4 (observe: 2, observeEvent: 2)
- **Render Functions:** 12

#### mod_visualization_enhanced.R (30K, 823 lines)
- **UI Function:** mod_visualization_enhanced_ui()
- **Server Function:** mod_visualization_enhanced_server()
- **Reactive Expressions:** 2
- **Observers:** 4
- **Render Functions:** 10

#### mod_visualization_enhanced_fixed.R (26K, 723 lines)
- **UI Function:** mod_visualization_enhanced_fixed_ui()
- **Server Function:** mod_visualization_enhanced_fixed_server()
- **Reactive Expressions:** 2
- **Observers:** 3
- **Render Functions:** 10

---

## 4. Main App Files

### 4.1 Primary Application Files

#### inst/shiny/app.R (varies)
**Primary Shiny application entry point**
- Sources global.R for configuration
- Defines UI structure
- Defines server logic
- Integrates all modules
- **Key Functions:**
  - Main UI definition
  - Main server function
  - Module integration
  - Session management

#### inst/shiny/global.R (varies)
**Global configuration and library loading**
- **Key Components:**
  - Library loading (shiny, shinyjs, DT, plotly, ggplot2, dplyr, etc.)
  - Environment variable configuration
  - APP_CONFIG list definition
  - Global utility functions
  - Data path configuration

#### inst/shiny/app_perturbseq_full.R (varies)
**Full Perturb-seq specific application**
- Dedicated interface for Perturb-seq datasets
- Simplified UI without mutation-related controls
- Specialized for pooled MixScale data

---

## 5. Architecture Summary

### 5.1 Application Entry Points

**Primary Entry Point:**
```r
launch_app() → launch_iscore_app() → inst/shiny/app.R
```

**Alternative Entry Points:**
```r
run_app() → Shiny app launch
run_complete_pipeline() → Full analysis pipeline
```

### 5.2 Data Flow Architecture

```
1. User calls launch_app()
   ↓
2. Dataset validation (validate_dataset_directory)
   ↓
3. Environment variables set (ISCORE_DATA_DIR, ISCORE_DE_FILE, etc.)
   ↓
4. Shiny app launched (shiny::runApp)
   ↓
5. global.R loads libraries and configuration
   ↓
6. Modules loaded and initialized
   ↓
7. Data manager initialized (init_data_manager)
   ↓
8. Reactive data loading begins
```

### 5.3 Key Data Files

**Required Dataset Files:**
- `full_DE_results.rds` - Differential expression results
- `all_enrichment_padj005_complete_with_direction.rds` - Enrichment results
- `enrichment_results/` - Directory with detailed enrichment results

**Optional Files:**
- UMAP coordinates
- Gene association data
- Signature nomination results
- Convergence analysis results

### 5.4 Module Dependency Graph

**Core Dependencies:**
- `mod_data_loader` → All other modules (provides data)
- `mod_de_results` → `mod_heatmap`, `mod_export`
- `mod_signature_nomination` → `mod_signature_trends`
- `mod_enrichment_gene_display` → `mod_pathview`

### 5.5 Reactive Data Flow

**Primary Reactive Chain:**
```
User Input → Reactive Filters → Data Manager → Module Reactives → Render Functions → UI Display
```

**Cache Integration:**
```
Data Request → Check Cache → Load from Disk (if needed) → Update Cache → Return Data
```

---

## 6. Component Categories by Function

### 6.1 Data Import & Validation (14 functions)
- Import functions (MAST, MixScale, Enrichment)
- Format detection (pooled vs experiment-split)
- Dataset validation
- File generation

### 6.2 Statistical Analysis (31 functions)
- Fisher's exact test
- Direction analysis
- Overlap significance
- Effect size correlation
- Meta-analysis
- FDR correction

### 6.3 Enrichment Analysis (20 functions)
- GO, KEGG, Reactome, WikiPathways
- STRING PPI
- GSEA
- Term extraction
- Pathway comparison

### 6.4 Signature Discovery (25 functions)
- Gene pair signatures
- Pathway overlap
- PD relevance scoring
- Biological interpretation
- Manuscript summary generation

### 6.5 Visualization (29 functions)
- Heatmaps (static, interactive, bubble)
- Volcano plots
- Scatter plots
- Bar plots, dot plots
- Multi-metric dashboards
- UMAP viewers

### 6.6 User Interface (43 modules)
- Data loading
- Result display
- Interactive filtering
- Export functionality
- Comparison tools

### 6.7 Performance & Optimization (11 functions)
- Benchmarking
- Memory profiling
- Lazy loading
- Caching
- Data sampling

### 6.8 Configuration & Management (18 functions)
- App configuration
- Dataset management
- Agent coordination
- Startup management

---

## 7. Documentation Priority Matrix

### Priority Level 1 - CRITICAL (User-Facing, Exported)
**Target: 100% roxygen2 documentation**

- All `@export` functions (118 functions)
- All Shiny module UI/Server functions (43 functions)
- Main app launchers (launch_app, run_app)
- **Total: ~161 components**

### Priority Level 2 - HIGH (Internal, Frequently Used)
**Target: 95% roxygen2 documentation**

- Helper functions in inst/shiny/R/ (45 functions)
- Major internal functions (no @export but complex)
- Data manager functions
- Cache manager functions
- **Total: ~74 components**

### Priority Level 3 - MEDIUM (Reactive Components)
**Target: 85% inline documentation**

- Reactive expressions (37)
- Major observers (76 observe + observeEvent)
- **Total: ~113 components**

### Priority Level 4 - LOWER (Render Functions)
**Target: 70% inline documentation**

- Render functions (143)
- Simple utility functions
- **Total: ~143 components**

---

## 8. Next Steps for Documentation Project

### Phase 1: Function-Level Documentation
1. Review all 192 R package functions
2. Identify functions missing roxygen2 documentation
3. Create documentation templates for each category
4. Document all exported functions (Priority 1)

### Phase 2: Module Documentation
1. Document all 23 Shiny modules
2. Create module architecture diagrams
3. Document reactive data flows
4. Create user guides for each module

### Phase 3: Workflow Documentation
1. Document complete analysis workflows
2. Create vignettes for common use cases
3. Document data requirements
4. Create troubleshooting guides

### Phase 4: Package-Level Documentation
1. Update README.md
2. Create comprehensive vignettes
3. Build pkgdown website
4. Create developer documentation

---

## 9. Key Insights

### Largest Components
1. **mod_signature_nomination.R** - 149K, 3213 lines (most complex module)
2. **mod_de_results.R** - 126K, 3108 lines (core DE analysis)
3. **enrichment_analysis.R** - 80K, 2080 lines (comprehensive enrichment)
4. **mod_heatmap.R** - 69K, 1723 lines (advanced heatmaps)
5. **signature_analysis.R** - 45K, 1120 lines (signature discovery)

### Most Complex Modules (by component count)
1. **mod_de_results** - 5 reactives, 9 observers, 19 renders
2. **mod_signature_nomination** - 2 reactives, 6 observers, 16 renders
3. **mod_landing_page_with_umap_v2** - 0 reactives, 7 observers, 16 renders
4. **mod_heatmap** - 4 reactives, 8 observers, 5 renders
5. **mod_visualization** - 2 reactives, 4 observers, 12 renders

### Documentation Coverage Gaps
- **Missing roxygen2:** ~30-40% of functions
- **Missing examples:** ~60% of functions
- **Missing vignettes:** No comprehensive vignettes exist
- **Missing module docs:** Limited module-level documentation

---

## 10. Recommendations

### Documentation Strategy
1. **Start with exported functions** - Highest user impact
2. **Document by category** - Easier to maintain consistency
3. **Create templates** - Standardize documentation format
4. **Use examples liberally** - Essential for user adoption
5. **Build vignettes** - Critical for workflow understanding

### Tooling Recommendations
1. **Use roxygen2** for all R functions
2. **Use pkgdown** for website generation
3. **Use DT** for interactive documentation tables
4. **Use mermaid** for architecture diagrams

### Quality Assurance
1. Run `devtools::check()` regularly
2. Validate all examples with `devtools::run_examples()`
3. Build and review pkgdown site locally
4. Get user feedback on documentation clarity

---

**End of Phase 0 Function Inventory**

*This inventory provides the foundation for comprehensive documentation modernization of the iSCORE-PDecipher package. All 536 components have been catalogued and categorized for systematic documentation improvement.*
