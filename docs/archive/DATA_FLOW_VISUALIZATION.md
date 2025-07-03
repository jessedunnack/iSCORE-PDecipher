# iSCORE-PDecipher Data Flow Visualization

## Core Data Processing Pipeline

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                              RAW DATA SOURCES                                  │
├─────────────────────────────────────────────────────────────────────────────────┤
│                                                                                 │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐                     │
│  │ MAST Results │    │CRISPRi Results│   │CRISPRa Results│                     │
│  │              │    │               │   │               │                     │
│  │ 13 mutations │    │ 10 genes      │   │ 10 genes      │                     │
│  │ 14 clusters  │    │ 10 clusters   │   │ 10 clusters   │                     │
│  │ .xlsx files  │    │ 3 experiments │   │ 3 experiments │                     │
│  └──────┬───────┘    └───────┬───────┘   └───────┬───────┘                     │
│         │                    │                   │                             │
└─────────┼────────────────────┼───────────────────┼─────────────────────────────┘
          │                    │                   │
          ▼                    ▼                   ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                           DATA IMPORT LAYER                                    │
├─────────────────────────────────────────────────────────────────────────────────┤
│                                                                                 │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐                     │
│  │import_mast_  │    │import_mixscale│   │import_mixscale│                     │
│  │data()        │    │_data()        │   │_data()        │                     │
│  │              │    │ (CRISPRi)     │   │ (CRISPRa)     │                     │
│  │• Gene symbol │    │• Exp detection│   │• Modality     │                     │
│  │  validation  │    │• Background   │   │  detection    │                     │
│  │• Cluster ID  │    │  gene tracking│   │• Column       │                     │
│  │  extraction  │    │• Log2FC cols  │   │  mapping      │                     │
│  │• Statistical │    │• P-value cols │   │• Statistical  │                     │
│  │  formatting  │    │• Metadata     │   │  formatting   │                     │
│  └──────┬───────┘    └───────┬───────┘   └───────┬───────┘                     │
│         │                    │                   │                             │
│         └────────────────────┼───────────────────┘                             │
│                              │                                                 │
└──────────────────────────────┼─────────────────────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                      UNIFIED DATA STRUCTURE                                    │
├─────────────────────────────────────────────────────────────────────────────────┤
│                    full_DE_results.rds (190MB)                                 │
│                                                                                 │
│  ┌─────────────────────────────────────────────────────────────────────────┐   │
│  │ list(                                                                   │   │
│  │   iSCORE_PD_MAST = list(                                               │   │
│  │     "LRRK2" = list(                                                    │   │
│  │       "cluster_0" = list(results = data.frame, background_genes = ..) │   │
│  │       "cluster_1" = ...                                               │   │
│  │     ),                                                                 │   │
│  │     "DNAJC6" = ...,                                                   │   │
│  │     ...                                                               │   │
│  │   ),                                                                  │   │
│  │   CRISPRi_Mixscale = list(...),                                      │   │
│  │   CRISPRa_Mixscale = list(...)                                       │   │
│  │ )                                                                     │   │
│  └─────────────────────────────────────────────────────────────────────────┘   │
└──────────────────────────────┼─────────────────────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                       ENRICHMENT ANALYSIS LAYER                                │
├─────────────────────────────────────────────────────────────────────────────────┤
│                                                                                 │
│  ┌──────────────┐ ┌──────────────┐ ┌──────────────┐ ┌──────────────┐          │
│  │   GO Terms   │ │ KEGG Paths   │ │  Reactome    │ │   STRING     │          │
│  │              │ │              │ │              │ │              │          │
│  │• Biological  │ │• Metabolic   │ │• Pathways    │ │• Protein     │          │
│  │  Process     │ │  pathways    │ │• Reactions   │ │  interactions│          │
│  │• Cellular    │ │• Disease     │ │• Complexes   │ │• Functional  │          │
│  │  Component   │ │  pathways    │ │              │ │  networks    │          │
│  │• Molecular   │ │              │ │              │ │              │          │
│  │  Function    │ │              │ │              │ │              │          │
│  └──────┬───────┘ └──────┬───────┘ └──────┬───────┘ └──────┬───────┘          │
│         │                │                │                │                   │
│         └────────────────┼────────────────┼────────────────┘                   │
│                          │                │                                    │
│  ┌──────────────┐ ┌──────┼──────┐ ┌───────┼──────┐                            │
│  │ WikiPathways │ │   GSEA      │ │ Fisher Exact │                            │
│  │              │ │             │ │   Testing    │                            │
│  │• Curated     │ │• Hallmark   │ │              │                            │
│  │  biological  │ │  gene sets  │ │• Hypergeom   │                            │
│  │  pathways    │ │• C2/C5 sets │ │  testing     │                            │
│  │              │ │• Running    │ │• P-value     │                            │
│  │              │ │  enrichment │ │  adjustment  │                            │
│  └──────┬───────┘ └──────┬──────┘ └──────┬───────┘                            │
│         │                │               │                                    │
│         └────────────────┼───────────────┘                                    │
│                          │                                                    │
└──────────────────────────┼────────────────────────────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                     RESULTS CONSOLIDATION                                      │
├─────────────────────────────────────────────────────────────────────────────────┤
│                                                                                 │
│  Processing 14,052 individual enrichment result files                          │
│  Standardizing column names and data formats                                   │
│  Adding direction metadata (UP/DOWN/ALL/RANKED)                                │
│  Including method information (MAST/MixScale/CRISPRa)                          │
│  Filtering to significant results (p.adjust < 0.05)                            │
│                                                                                 │
│  ┌─────────────────────────────────────────────────────────────────────────┐   │
│  │          all_enrichment_padj005_complete_with_direction.rds             │   │
│  │                            (266MB)                                     │   │
│  │                                                                         │   │
│  │ 767,337 significant enrichment terms:                                  │   │
│  │ • STRING: 661,608 terms (86.2%)                                        │   │
│  │ • GO: 78,524 terms (10.2%)                                             │   │
│  │ • Reactome: 19,234 terms (2.5%)                                        │   │
│  │ • GSEA: 8,691 terms (1.1%)                                             │   │
│  │ • KEGG: 3,088 terms (0.4%)                                             │   │
│  │ • WikiPathways: 2,329 terms (0.3%)                                     │   │
│  │                                                                         │   │
│  │ Direction distribution:                                                 │   │
│  │ • UP-regulated: 298,916 terms (39.0%)                                  │   │
│  │ • DOWN-regulated: 226,779 terms (29.5%)                                │   │
│  │ • ALL genes: 232,951 terms (30.4%)                                     │   │
│  │ • RANKED (GSEA): 8,691 terms (1.1%)                                    │   │
│  └─────────────────────────────────────────────────────────────────────────┘   │
│                                                                                 │
└──────────────────────────────┼─────────────────────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                        SHINY VISUALIZATION LAYER                               │
├─────────────────────────────────────────────────────────────────────────────────┤
│                                                                                 │
│  ┌───────────────┐    ┌───────────────┐    ┌───────────────┐                  │
│  │ Global Filter │────▶│ Data Reactive │────▶│ Module Render │                  │
│  │               │    │               │    │               │                  │
│  │• Gene/Mutation│    │• Real-time    │    │• Plot         │                  │
│  │• Cluster      │    │  filtering    │    │  generation   │                  │
│  │• Enrichment   │    │• Validation   │    │• Interactive  │                  │
│  │  Type         │    │• Caching      │    │  features     │                  │
│  │• Direction    │    │• Performance  │    │• Export       │                  │
│  │• Analysis     │    │  optimization │    │  capabilities │                  │
│  │  Method       │    │               │    │               │                  │
│  └───────────────┘    └───────────────┘    └───────────────┘                  │
│                                                                                 │
└─────────────────────────────────────────────────────────────────────────────────┘

```

## Module Interaction Diagram

```
┌────────────────────────────────────────────────────────────────────────────────┐
│                          SHINY APP ARCHITECTURE                               │
├────────────────────────────────────────────────────────────────────────────────┤
│                                                                                │
│  ┌─────────────┐         ┌─────────────┐         ┌─────────────┐              │
│  │   app.R     │◄────────┤  global.R   │────────▶│ startup_    │              │
│  │             │         │             │         │ manager.R   │              │
│  │• UI layout  │         │• Libraries  │         │• Data       │              │
│  │• Server     │         │• Config     │         │  loading    │              │
│  │  logic      │         │• Utilities  │         │• Validation │              │
│  │• Module     │         │• Data       │         │• Path       │              │
│  │  loading    │         │  structures │         │  detection  │              │
│  └──────┬──────┘         └─────────────┘         └─────────────┘              │
│         │                                                                     │
│         ▼                                                                     │
│  ┌─────────────────────────────────────────────────────────────────────────┐ │
│  │                            MODULE LAYER                                │ │
│  ├─────────────────────────────────────────────────────────────────────────┤ │
│  │                                                                         │ │
│  │ ┌─────────────┐  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐     │ │
│  │ │ Landing     │  │Functional   │  │ Heatmaps    │  │ Signature   │     │ │
│  │ │ Page        │  │Enrichment   │  │             │  │ Nomination  │     │ │
│  │ │             │  │             │  │             │  │             │     │ │
│  │ │• UMAP       │  │• Dot plots  │  │• Interactive│  │• Cross-     │     │ │
│  │ │• Summary    │  │• Bar plots  │  │  heatmaply  │  │  method     │     │ │
│  │ │• Statistics │  │• Lollipops  │  │• Clustering │  │  comparison │     │ │
│  │ │• Markers    │  │• GSEA plots │  │• Annotations│  │• Statistical│     │ │
│  │ └─────────────┘  └─────────────┘  └─────────────┘  └─────────────┘     │ │
│  │                                                                         │ │
│  │ ┌─────────────┐  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐     │ │
│  │ │ DE Results  │  │ Export      │  │ Comparison  │  │ Help/Docs   │     │ │
│  │ │             │  │             │  │             │  │             │     │ │
│  │ │• UMAP       │  │• Multi-     │  │• Method     │  │• Tutorials  │     │ │
│  │ │• Volcano    │  │  format     │  │  comparison │  │• Context    │     │ │
│  │ │  plots      │  │  download   │  │• Statistical│  │  help       │     │ │
│  │ │• Static/    │  │• Metadata   │  │  testing    │  │• Examples   │     │ │
│  │ │  Interactive│  │  inclusion  │  │• Overlaps   │  │             │     │ │
│  │ └─────────────┘  └─────────────┘  └─────────────┘  └─────────────┘     │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                                                                │
│  ┌─────────────────────────────────────────────────────────────────────────┐ │
│  │                         REACTIVE DATA FLOW                             │ │
│  ├─────────────────────────────────────────────────────────────────────────┤ │
│  │                                                                         │ │
│  │  Global Settings ────▶ Data Filter ────▶ Module Processing ────▶ Render │ │
│  │        │                    │                      │                │  │ │
│  │        ▼                    ▼                      ▼                ▼  │ │
│  │  • Gene select        • p-value filter      • Plot data prep   • ggplot │ │
│  │  • Cluster select     • Direction filter    • Title generation • plotly │ │
│  │  • Enrichment type    • Method filter       • Color mapping    • DT     │ │
│  │  • Analysis method    • Experiment filter   • Statistical calc • Export │ │
│  │                                                                         │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
└────────────────────────────────────────────────────────────────────────────────┘
```

## Performance Characteristics

### Data Processing Timeline
```
Data Import:        5-10 minutes   (One-time setup)
Enrichment Analysis: 1-2 hours     (One-time computation)
Results Processing: 10-30 minutes  (One-time consolidation)
App Startup:        10-15 seconds  (Data loading)
Interactive Response: <1 second    (Filtering/plotting)
Export Generation:  2-10 seconds   (Format dependent)
```

### Memory Usage Profile
```
Peak Memory Usage: ~1.2 GB (during enrichment analysis)
Steady State:      ~400 MB (Shiny app running)
Consolidated Data: 266 MB  (Primary results file)
UMAP Data:         14 MB   (Visualization coordinates)
Cache Storage:     50-100 MB (Filtered results cache)
```

### Bottleneck Analysis
1. **Initial Enrichment**: Most time-consuming step (hours)
2. **Data Loading**: App startup delay (10-15s)
3. **Large Result Sets**: Heatmap generation with >1000 terms
4. **Interactive Response**: Complex filtering operations
5. **Export Operations**: Large file generation

This visualization shows how data flows from raw differential expression results through comprehensive enrichment analysis to interactive visualization, highlighting the major processing steps and performance characteristics of each component.