# Shiny Integration Plan for PD Signatures
## Date: 2025-07-19

## Overview
Plan to add PD signature exploration capabilities to the existing Shiny app without conflicting with ongoing heatmap module enhancements.

## Key Findings to Highlight

### 1. **MAST-only Signature**: Vesicle-mediated transport in synapse
- Found in 8 genes (LRRK2, PARK7, PRKN, etc.)
- Mean -log10(p) = 6.34
- **Biological Story**: Genetic mutations converge on synaptic vesicle dysfunction

### 2. **MixScale-only Signature**: Mitochondrial small ribosomal subunit
- Found in 9 genes (ATP13A2, DNAJC6, FBXO7, etc.)
- Mean -log10(p) = 1.88
- **Biological Story**: Acute knockdowns reveal compensatory mitochondrial protein synthesis response

### 3. **Convergent Signature**: Synapse (STRING database)
- MAST: 13 genes, MixScale: 10 genes
- Mean -log10(p) = 4.75
- **Special Note**: "(2018) Dopamine perturbation" signature shows extreme significance (p = 17.05!)
- **Biological Story**: Core synaptic dysfunction is a universal feature across perturbation methods

## Proposed Implementation Strategy

### Option 1: Add to Existing Heatmap Module (Minimal Conflict)
Since another agent is enhancing the heatmap module, we could:
1. Add a "PD Signatures" preset button that automatically:
   - Filters for PD-relevant terms
   - Sets appropriate parameters
   - Highlights convergent signatures
2. Add these keywords to the term search functionality already being implemented

### Option 2: Create New Tab in Signature Nomination Section (Recommended)
Create `mod_pd_signature_explorer.R` with:
1. **Pre-computed Results Display**
   - Show top 3 signature categories
   - Interactive tables with gene lists
   - Biological interpretation text
2. **Visualization Panel**
   - Bar plots of top signatures
   - Method comparison plots
   - Download options for presentations
3. **Deep Dive Analysis**
   - Click on signature → see all genes/clusters
   - Export gene lists for further analysis

### Option 3: Quick Implementation - Add to Landing Page
Add a "Key PD Signatures" summary box showing:
- Top signature from each category
- Quick stats (genes, p-values)
- Link to full analysis

## Files to Create/Modify

### New Files:
1. `inst/shiny/modules/mod_pd_signature_explorer.R` - Main module
2. `inst/shiny/www/pd_signatures_precomputed.rds` - Pre-computed results

### Files to Modify (AFTER git sync):
1. `inst/shiny/app.R` - Add new tab/module
2. `inst/shiny/modules/mod_heatmap.R` - Add PD preset (if Option 1)
3. `inst/shiny/global.R` - Add any required functions

## Git Coordination Strategy
1. Check for recent commits: `git pull`
2. Review changes to `mod_heatmap.R`
3. Create feature branch: `git checkout -b feature/pd-signatures`
4. Implement changes
5. Test with existing heatmap enhancements
6. Create PR with clear description

## Quick Win for Tomorrow's Meeting
For immediate use without code changes:
1. Use existing heatmap module with manual filtering:
   - Term search: "synap|mitochond|autophagy|lysosom"
   - Method filter: "intersection" for convergent
   - Direction: "ALL" to see full picture
2. Generate static plots using visualization script
3. Create PowerPoint slides with key findings

## Key Talking Points for Committee
1. **Method Validation**: Convergent signatures validate our approach
2. **Novel Finding**: Dopamine perturbation network (p = 1e-17!)
3. **Therapeutic Targets**: Synaptic vesicle pathway genes
4. **Method-Specific Insights**: Why mutations vs knockdowns show different signatures