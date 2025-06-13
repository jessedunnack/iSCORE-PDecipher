# iSCORE-PDecipher Dataset Sharing Plan

## 📊 Available Dataset Combinations

You have **three pre-computed dataset options** for sharing with labmates:

### Option 1: iSCORE-PD Only (76 MB + 4.1 MB = 80.1 MB)
- **Enrichment Data**: `/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD/all_enrichment_padj005_complete_with_direction.rds` (76 MB)
- **UMAP Data**: `/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/inst/extdata/umap_data/iSCORE_PD_umap_data.rds` (4.1 MB)
- **Contains**: 13 genetic mutations only (LRRK2, GBA, PINK1, PRKN, etc.)
- **Best for**: Labmates focused only on genetic mutations

### Option 2: iSCORE-PD + CRISPRi (221 MB + 4.6 MB = 225.6 MB)
- **Enrichment Data**: `/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds` (221 MB)
- **UMAP Data**: `/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/inst/extdata/umap_data/iSCORE_PD_CRISPRi_umap_data.rds` (4.6 MB)
- **Contains**: Mutations + CRISPRi knockdowns
- **Best for**: Comparing mutations with gene knockdowns

### Option 3: Complete Dataset (266 MB + 14 MB = 280 MB)
- **Enrichment Data**: `/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/data_files/all_enrichment_padj005_complete_with_direction.rds` (266 MB)
- **UMAP Data**: `/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/inst/extdata/umap_data/all_umap_data_combined.rds` (14 MB)
- **Contains**: All mutations + CRISPRi + CRISPRa (767,337 terms)
- **Best for**: Complete analysis with all experimental approaches

## 📁 Recommended Google Drive Structure

Create folders for each option:

```
iSCORE-PDecipher-Datasets/
├── Option1_iSCORE-PD-Only/
│   ├── all_enrichment_padj005_complete_with_direction.rds (76 MB)
│   └── iSCORE_PD_umap_data.rds (4.1 MB)
│   
├── Option2_iSCORE-PD-plus-CRISPRi/
│   ├── all_enrichment_padj005_complete_with_direction.rds (221 MB)
│   └── iSCORE_PD_CRISPRi_umap_data.rds (4.6 MB)
│   
├── Option3_Complete-Dataset/
│   ├── all_enrichment_padj005_complete_with_direction.rds (266 MB)
│   └── all_umap_data_combined.rds (14 MB)
│   
└── Quickstart_Guide.pdf
```

## 🎯 Which Option to Choose?

**For most labmates**: **Option 3 (Complete Dataset)** - only 280 MB total and gives access to everything.

**For specific use cases**:
- **Genetics focus only**: Option 1 (80 MB)
- **Mutations vs knockdowns**: Option 2 (226 MB)
- **Complete analysis**: Option 3 (280 MB)

## 📧 Email Template for Labmates

```
Subject: iSCORE-PDecipher Analysis Tool - Ready to Use!

Hi [Name],

I've prepared the iSCORE-PDecipher Shiny app for analyzing our Parkinson's disease data. You can now explore 767K+ enrichment terms across mutations, CRISPRi, and CRISPRa experiments.

**What you get:**
- Interactive visualizations (UMAP, heatmaps, dotplots)
- 13 genetic mutations + 10 CRISPR perturbations
- 6 enrichment databases (GO, KEGG, Reactome, etc.)
- Publication-quality figure export

**Setup (5 minutes):**
1. Download files from Google Drive: [LINK]
2. Follow the Quickstart Guide (attached PDF)
3. Install from GitHub and launch

**Google Drive folder:**
[INSERT YOUR GOOGLE DRIVE LINK]

**Need help?** Slack me or check the troubleshooting section in the guide.

Happy analyzing!
Jesse
```

## 🚀 Next Steps

1. **Choose your sharing strategy** (recommend Option 3 for most users)
2. **Upload files to Google Drive** using the folder structure above
3. **Convert the quickstart guide to PDF** (I'll help with this next)
4. **Test the download/install process** with one labmate first
5. **Send the email template** to your lab

## 📋 File Copy Commands

To prepare files for upload:

```bash
# Option 1: iSCORE-PD Only
mkdir -p ~/Desktop/iSCORE_Sharing/Option1_iSCORE-PD-Only/
cp "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD/all_enrichment_padj005_complete_with_direction.rds" ~/Desktop/iSCORE_Sharing/Option1_iSCORE-PD-Only/
cp "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/inst/extdata/umap_data/iSCORE_PD_umap_data.rds" ~/Desktop/iSCORE_Sharing/Option1_iSCORE-PD-Only/

# Option 2: iSCORE-PD + CRISPRi  
mkdir -p ~/Desktop/iSCORE_Sharing/Option2_iSCORE-PD-plus-CRISPRi/
cp "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds" ~/Desktop/iSCORE_Sharing/Option2_iSCORE-PD-plus-CRISPRi/
cp "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/inst/extdata/umap_data/iSCORE_PD_CRISPRi_umap_data.rds" ~/Desktop/iSCORE_Sharing/Option2_iSCORE-PD-plus-CRISPRi/

# Option 3: Complete Dataset
mkdir -p ~/Desktop/iSCORE_Sharing/Option3_Complete-Dataset/
cp "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/data_files/all_enrichment_padj005_complete_with_direction.rds" ~/Desktop/iSCORE_Sharing/Option3_Complete-Dataset/
cp "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/inst/extdata/umap_data/all_umap_data_combined.rds" ~/Desktop/iSCORE_Sharing/Option3_Complete-Dataset/
```

Now your labmates can choose exactly what they need!