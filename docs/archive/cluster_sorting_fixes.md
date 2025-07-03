# Cluster Sorting Fixes

## Problem
Clusters were being sorted alphabetically as strings, causing cluster_10, cluster_11, etc. to appear immediately after cluster_1 instead of after cluster_9.

## Solution
Added a `natural_sort_clusters()` function that extracts numeric parts and sorts clusters numerically.

## Files Modified

1. **inst/shiny/R/data_manager.R**
   - Added `natural_sort_clusters()` function
   - Updated `get_available_clusters()` to use natural sorting

2. **inst/shiny/app.R**
   - Line 912: Changed `sort(unique(filtered_data$cluster))` to `natural_sort_clusters(unique(filtered_data$cluster))`

3. **inst/shiny/modules/mod_de_results.R**
   - Line 340: Updated cluster sorting in UMAP data loading
   - Line 490: Updated cluster sorting in cluster selector update

4. **inst/shiny/modules/mod_signature_nomination.R**
   - Line 377: Changed `sort(unique_clusters[!is.na(unique_clusters)])` to `natural_sort_clusters(unique_clusters[!is.na(unique_clusters)])`

## Hardcoded Issues Found
- **mod_de_analysis.R** line 51: Has hardcoded `choices = paste0("cluster_", 0:13)` - This may need to be made dynamic if cluster counts vary between datasets.

## Testing
After these changes, cluster selection dropdowns throughout the app should show:
- cluster_0, cluster_1, ..., cluster_9, cluster_10, cluster_11, ..., cluster_14
Instead of:
- cluster_0, cluster_1, cluster_10, cluster_11, ..., cluster_14, cluster_2, ...