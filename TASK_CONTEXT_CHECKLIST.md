# Task Context Checklist
**Always consult before starting significant work**

## 🔍 Pre-Task Context Check

### Step 1: Read Main Knowledge Base
- [ ] **Review [CLAUDE.md](../CLAUDE.md)** - Main knowledge base and current status
- [ ] **Check Executive Summary** - Latest achievements and manuscript candidates
- [ ] **Review User Preferences** - Environment setup and workflow requirements

### Step 2: Task-Specific Context

#### For Gene Analysis Tasks:
- [ ] **Check [COMPREHENSIVE_PD_GENE_ANALYSIS_2025.md](../COMPREHENSIVE_PD_GENE_ANALYSIS_2025.md)**
  - Current gene rankings (ATP13A2, PARK7, FBXO7, LRRK2 top tier)
  - Convergence patterns (4 types identified)
  - Biological mechanisms and mutation types
- [ ] **Review manuscript candidates** - Don't duplicate existing analysis

#### For Statistical/Correlation Tasks:
- [ ] **Check [ENHANCED_METHODOLOGY_AND_MIGRATION_GUIDE.md](../ENHANCED_METHODOLOGY_AND_MIGRATION_GUIDE.md)**
  - Use enhanced Fisher's exact tests (intersection/union backgrounds)
  - Apply gene filtering methodology (Top 200 genes = 11x improvement)
  - Use hierarchical FDR correction
  - Direction-aware analysis expectations

#### For Data Processing Tasks:
- [ ] **Check data column names in CLAUDE.md**:
  - MAST: `mutation_tidy`, `eWT` control, `batch`
  - MixScale: `scMAGeCK_gene_assignment`, `Non-Targeting` control
  - Experiments: `log2FC_<experiment>`, `p_cell_type<experiment>:weight`
- [ ] **Check existing data files**:
  - `full_DE_results.rds` (190 MB)
  - `all_enrichment_padj005_complete_with_direction.rds` (266 MB)

#### For App/Visualization Tasks:
- [ ] **Check Recent Updates & Fixes** in CLAUDE.md
  - Correlation plot improvements (vertical layout)
  - Heatmap module status (interactive heatmaply)
  - UMAP layout optimization (950×700px)
- [ ] **Review Common Tasks & Solutions** - Use existing patterns

#### For Cross-Platform Tasks:
- [ ] **Check config_manager.R** functionality
- [ ] **Use `file.path()`, never hardcode paths**
- [ ] **Check platform-specific guidelines**

### Step 3: Check for Existing Solutions
- [ ] **Review Technical Architecture** - Don't recreate existing modules
- [ ] **Check Common Tasks section** - Use established patterns
- [ ] **Review Development Guidelines** - Follow best practices

### Step 4: Validate Approach
- [ ] **Does this align with user preferences?** (batch processing, comprehensive reports)
- [ ] **Are we using the latest statistical methods?** (v0.2.6+ enhancements)
- [ ] **Will this work cross-platform?** (Mac/Windows/Linux)
- [ ] **Are we building on recent fixes?** (avoid regression)

## ⚠️ Red Flags - STOP and Re-check

- **Hardcoding paths** → Use config_manager.R
- **Using old correlation methods** → Apply gene filtering (11x improvement)
- **Wrong Fisher's exact tests** → Use intersection/union backgrounds
- **Ignoring direction expectations** → Use direction-aware analysis
- **Recreating existing functionality** → Check existing modules first
- **Not following data structure** → Check column names in CLAUDE.md

## 🎯 Quick Reference Links

### Primary Documents:
- **[CLAUDE.md](../CLAUDE.md)** - Main knowledge base
- **[COMPREHENSIVE_PD_GENE_ANALYSIS_2025.md](../COMPREHENSIVE_PD_GENE_ANALYSIS_2025.md)** - Gene analysis details
- **[ENHANCED_METHODOLOGY_AND_MIGRATION_GUIDE.md](../ENHANCED_METHODOLOGY_AND_MIGRATION_GUIDE.md)** - Statistical methods

### Technical References:
- **Common Tasks & Solutions** (CLAUDE.md lines 119-141)
- **User Workflow Preferences** (CLAUDE.md lines 64-89)
- **Recent Updates & Fixes** (CLAUDE.md lines 127-145)
- **Development Guidelines** (CLAUDE.md lines 163-187)

## 📝 Context Notes Template

When starting a task, fill this out:

```
Task: [Brief description]
Context Checked: [✓] CLAUDE.md [✓] Gene Analysis [✓] Statistical Methods
Relevant Info Found:
- 
- 
- 
Existing Solutions to Build On:
- 
- 
Potential Conflicts/Issues:
- 
- 
```

---

**Remember**: 5 minutes of context checking can save hours of redundant work!