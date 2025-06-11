#' Run MixScale Differential Expression Analysis
#'
#' This function performs MixScale analysis for PerturbSeq experiments.
#'
#' @param experiment_path Path to the experiment directory
#' @param output_dir Directory to save results
#' @param modality Experiment modality ("CRISPRi" or "CRISPRa")
#'
#' @return List containing MixScale analysis results
#' @export
run_mixscale_analysis <- function(experiment_path, output_dir = ".", modality = "CRISPRi") {
  
  # Check required packages
  required_packages <- c("Seurat", "dplyr", "ggplot2", "Matrix", "magrittr", "mixtools")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }
  
  # Check additional packages
  additional_packages <- c("reshape2", "glmGamPoi", "gridExtra", "scMAGeCK")
  for (pkg in additional_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning(paste("Optional package", pkg, "is not installed. Some functionality may be limited."))
    }
  }
  
  # Set options
  options(future.globals.maxSize = 64 * 10^9)


## This is used to exclude genes that are not found in both MUT v eWT and PerturbSeq experiments
exclude <- c("SSR1", "SLC6A3", "CTSB", "PITX3", "APOE", "LRP1B", "GPNMB", "RAB7L1", "NEUROD1", "VPS35", "POLG", 
             "FOXA2", "NEUROG1", "NEUROG2", "PAX6", "EIF4G1", "RHEB", "GIGYF2", 
             "TMEM230", "HTRA2", "NEUROD2", "NR4A2", "DNAJC13", "RPTOR", "UCHL1", 
             "LRP10", "B2M", "SOX2", "ASCL1", 
             "PLA2G6", "LMX1A")



## These are genes that are unique to GWAS library G5O13. This omits any entries seen also in FPD Library G5O11.
gwas_genes <- c("Non-Targeting", "C8orf58", "SEC23IP", "NUPR1", "AREL1", "CASR", "ARIH2OS", "ASXL3", "PSD", "TMEM175", 
                "MAP4K4", "ZNF668", "CCDC58", "PARP9", "TOB2P1", "RIMS1", "ZSCAN16-AS1", "SAR1B", "HIST1H1B", "SPNS1", 
                "PAM", "SLC50A1", "INTS2", "RNF141", "FBRSL1", "ZNF2", "SATB1", "NCKIPSD", "HSD17B1", "BST1", 
                "MBNL2", "DDX46", "LAT", "CUEDC2", "USP4", "SCAF11", "GBA", "POLR2A", "ZNRD1ASP", "HLA-DRA", 
                "GALC", "DALRD3", "GUSB", "PGS1", "LCORL", "IP6K2", "BTNL2", "CTF1", "NFATC2IP", "HISLA", 
                "PPIP5K2", "MIR4697HG", "PCBD2", "FAM171A2", "MED12L", "ZNF629", "FAM47E", "KLHL7", "TUBG1", "MEX3C", 
                "COASY", "ZBTB7B", "PDLIM2", "CRHR1", "FCGR2A", "SLC45A3", "TUFM", "MRPS5", "CYLD", "RETREG3", 
                "IGSF9B", "CASC16", "UBQLN4", "RPS12", "MTRNR2L8", "SPPL2B", "DGKQ", "RNF39", "ATP6V0A1", "ZNF165", 
                "SYT4", "NDUFAF2", "PPM1L", "DNAH17", "LINC00693", "HLA-DQB1", "PMVK", "HIST1H3J", "HLA-DQA2", "MLX", 
                "CAB39L", "RIT2", "TNFSF12", "DCAF16", "LAMB2P1", "KCNS3", "CD19", "SCARB2", "EEF1AKNMT", "SH3GL2", 
                "RABGEF1", "KCNN3", "FOXA1", "C5orf30", "BRIP1", "HCG23", "HSD3B7", "GRN", "AMPD3", "C5orf24", 
                "BAG3", "BCL7C", "EIF3C", "GS1-124K5.11", "SLC18B1", "TMEM248", "NOL4", "GCH1", "PRSS3", "SH3RF1", 
                "CCAR2", "EGR3", "WNT3", "P4HTM", "TPST1", "SEC24A", "NCAPG", "UBE2R2", "CCDC71", "KCNIP3", 
                "MAPT", "ART3", "NDUFAF3", "CCDC62", "FAM49B", "TUBG2", "FBXL19", "GAK", "UBAP2", "LINC01012", 
                "SLC26A1", "TMEM163", "MIR146B", "SLC25A21-AS1", "FDFT1", "ITGA2B", "MED13", "QRICH1", "CRCP", "LAMB2", 
                "CAVIN1", "TYW1", "WDR5B", "PLEKHH3", "FAM47E-STBD1", "DTX3L", "KPNA1", "SETD1A", "ANKRD20A8P", "ZNF514", 
                "PROM2", "NUPL2", "DLST", "PRKAR2A", "WDR6", "TRIM40", "VKORC1L1", "P2RY12", "ARID2", "UBTF", 
                "TXNDC15", "FAM184B", "HIST1H2BN", "IGF2BP3", "GPATCH8", "STK39", "HCG17", "FCF1", "USP19", "DCUN1D1", 
                "SYT17", "NUCKS1", "PPP1R11", "NOLC1", "SBDS", "ITGA8", "MCCC1", "CLCN3", "BIN3", "GBF1", 
                "HIST1H2BL", "RABEP2", "HIST1H2AK", "HIP1R", "VAMP4", "GIN1", "FBXL19-AS1", "CATSPER3", "ZBTB4", "PBXIP1", 
                "KRTCAP2", "ARIH2", "NOD2", "CAMK2D", "TRIM26", "HIST1H2BO", "MAL", "NEK1", "PRSS53", "ZSCAN16", 
                "PGF", "ERCC8", "CPLX1", "CCT6P1", "OR2B6", "NFKB2", "CHD9", "MALSU1", "TOX3", "SULT1A1", 
                "HLA-DRB5", "ATP2A1", "FGF20", "CHRNB1", "KLHDC8B", "SH2B1", "CAMLG", "SNCA-AS1", "LOC442028", "TRIM15", 
                "ELOVL7", "ITPKB", "SLC25A20", "DYRK1A", "B3GALNT1", "ASH1L", "ASL", "PART1", "DLG2", "QARS", 
                "HIST1H3H", "SGF29", "DCAF12", "WDHD1", "RAB29", "ZNF646", "LRRN4", "DPM3", "DEPDC1B", "GTF2IRD1P1", 
                "NAGLU", "KLHL7-DT", "TRIM31", "FAM162A", "HIST1H2BM", "LINC02067", "MMRN1", "JADE2", "NMD3", "ZKSCAN8", 
                "INTS4P2", "STX4", "MIPOL1", "ELOVL3", "GBAP1", "HLA-DRB6", "CSTA", "PROX2", "CCDC36", "SMIM15", 
                "HIST1H2AM", "SIPA1L2", "C3orf84", "HLA-DQA1", "YLPM1", "INPP5F", "CRLS1", "RPS6KL1", "ZNF192P1", "LOC339862", 
                "VKORC1", "PTENP1", "SNX20", "LINC00174", "STX1B", "LOC100131289", "ATXN2L", "PMS2P4", "ZSCAN12P1", "FYN", 
                "HIST1H4J", "GPR65", "OR2B2", "HIST1H4L", "TSBP1", "TIAL1", "SLC4A1", "SPTSSB", "441242", "TRIM10")

process_mixscale <- function(obj,
                              phrase,
                              split_by,
                              neighbors,
                              omit,
                              omit_genes = NULL,
                              min_nt_cells = 5,
                              pca_dims = 50,
                              var_features = 3000,
                              logfc_threshold = 0.25,
                              max_de_genes = 500,
                              modality = "both") {  # New parameter: "CRISPRi", "CRISPRa", or "both"

    # Input validation
    if (!inherits(obj, "Seurat")) {
      stop("Input must be a Seurat object")
    }

    required_cols <- c("scMAGeCK_multiplet_id",
                       "cellranger_multiplet_id",
                       "scMAGeCK_gene_assignment",
                       "experiment",
                       "crispr_modality")  # Added requirement for modality column
    if (!all(required_cols %in% colnames(obj@meta.data))) {
      stop("Required metadata columns missing from object")
    }

    # Validate modality parameter
    if (!modality %in% c("CRISPRi", "CRISPRa", "both")) {
      stop("modality must be 'CRISPRi', 'CRISPRa', or 'both'")
    }

    # Check which modalities are present in the data
    available_modalities <- unique(obj@meta.data$crispr_modality)
    available_modalities <- available_modalities[!is.na(available_modalities)]
    
    message(sprintf("Available modalities in data: %s", paste(available_modalities, collapse = ", ")))

    # Determine which modalities to process
    if (modality == "both") {
      modalities_to_process <- intersect(c("CRISPRi", "CRISPRa"), available_modalities)
      if (length(modalities_to_process) == 0) {
        stop("No valid CRISPR modalities found in data")
      }
    } else {
      if (!modality %in% available_modalities) {
        stop(sprintf("Requested modality '%s' not found in data. Available: %s", 
                     modality, paste(available_modalities, collapse = ", ")))
      }
      modalities_to_process <- modality
    }

    # Create base output directory
    base_output_dir <- paste0("PerturbSeq_MixScale_analysis_", phrase)
    dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)

    # Process each modality
    for (current_modality in modalities_to_process) {
      message(sprintf("\n========== Processing %s modality ==========", current_modality))
      
      # Create modality-specific directory
      modality_dir <- file.path(base_output_dir, current_modality)
      dir.create(modality_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Log file setup for this modality
      log_file <- file.path(modality_dir, sprintf("%s_analysis_log.txt", current_modality))
      cat(sprintf("%s Analysis started: %s\n", current_modality, Sys.time()), file = log_file)
      
      # Filter object by modality
      obj_modality <- subset(obj, crispr_modality == current_modality)
      
      message(sprintf("Processing %d cells with %s", ncol(obj_modality), current_modality))
      
      # Store original cluster identities
      original_clusters <- obj_modality[[split_by]]
      obj_modality$seurat_clusters <- original_clusters

      # Process valid clusters
      clusters <- unique(obj_modality$seurat_clusters)
      clusters <- clusters[!clusters %in% omit]
      clusters <- as.character(sort(as.numeric(as.character(clusters))))

      # Initialize analysis tracking for this modality
      analysis_log <- data.frame()

      # Process each cluster
      for (clust in clusters) {
        message(sprintf("\nProcessing %s - Cluster %s", current_modality, clust))

        # Create cluster-specific directory within modality directory
        cluster_dir <- file.path(modality_dir,
                                 sprintf("all_%s_%s_Cluster%s", current_modality, phrase, clust))
        dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)

        tryCatch({
          # Subset cluster and filter multiplets
          cluster_obj <- subset(obj_modality,
                                seurat_clusters == clust &
                                  scMAGeCK_multiplet_id != "multiplet: DIFFERENT targets" &
                                  cellranger_multiplet_id != "multiplet: DIFFERENT targets")

          # Filter Non-Targeting cells
          message("Checking Non-Targeting cells by experiment:")
          nts_tab <- table(cluster_obj$scMAGeCK_gene_assignment,
                           cluster_obj$experiment)["Non-Targeting", ]
          print(nts_tab)

          # Remove experiments with insufficient NT cells
          toss <- names(nts_tab)[nts_tab <= neighbors]
          if (length(toss) > 0) {
            message(sprintf("Removing experiments with too few NT cells: %s",
                            paste(toss, collapse = ", ")))
            cluster_obj <- subset(cluster_obj, experiment %in% toss, invert = TRUE)
          }

          # Update NT table after filtering
          nts_tab <- table(cluster_obj$scMAGeCK_gene_assignment,
                           cluster_obj$experiment)["Non-Targeting",]

          if (length(nts_tab) == 0) {
            message(sprintf("No valid experiments remain for %s cluster %s", current_modality, clust))
            next
          }

          # Process data
          message("Processing expression data...")
          cluster_obj[["RNA"]] <- split(cluster_obj[["RNA"]],
                                        f = cluster_obj$experiment)

          cluster_obj <- NormalizeData(cluster_obj, assay = "RNA")
          cluster_obj <- ScaleData(cluster_obj,
                                   features = rownames(cluster_obj@assays$RNA),
                                   assay = "RNA")
          cluster_obj <- FindVariableFeatures(cluster_obj,
                                              assay = "RNA",
                                              nfeatures = var_features)
          cluster_obj <- RunPCA(cluster_obj,
                                assay = "RNA",
                                features = VariableFeatures(cluster_obj,
                                                            assay = "RNA"))
          cluster_obj <- FindNeighbors(cluster_obj,
                                       assay = "RNA",
                                       dims = 1:pca_dims,
                                       features = VariableFeatures(cluster_obj,
                                                                   assay = "RNA"))
          cluster_obj <- JoinLayers(cluster_obj, assay = "RNA")

          # Calculate perturbation signatures
          message("Calculating perturbation signatures...")
          nts <- table(cluster_obj$scMAGeCK_gene_assignment)["Non-Targeting"]
          message(sprintf("%d total Non-Targeting guides in %s cluster %s",
                          nts, current_modality, clust))

          cluster_obj <- CalcPerturbSig(
            object = cluster_obj,
            assay = "RNA",
            features = rownames(cluster_obj@assays$RNA),
            slot = "data",
            gd.class = "scMAGeCK_gene_assignment",
            nt.cell.class = "Non-Targeting",
            reduction = "pca",
            ndims = pca_dims,
            num.neighbors = neighbors,
            new.assay.name = "PRTB",
            split.by = "experiment")

          # Run MixScale
          message("Running MixScale analysis...")
          cluster_obj <- RunMixscale(
            object = cluster_obj,
            assay = "PRTB",
            slot = "scale.data",
            labels = "scMAGeCK_gene_assignment",
            nt.class.name = "Non-Targeting",
            min.de.genes = 1,
            logfc.threshold = logfc_threshold,
            de.assay = "RNA",
            max.de.genes = max_de_genes,
            new.class.name = "mixscale_score",
            fine.mode = TRUE,
            fine.mode.labels = "scMAGeCK_gRNA_assignment",
            harmonize = TRUE,
            verbose = TRUE,
            split.by = "experiment")

          # Save processed object
          saveRDS(cluster_obj,
                  file.path(cluster_dir,
                            sprintf("%s_%s_Cluster%s_mixscale_object.rds",
                                    current_modality, phrase, clust)))

          # Run differential expression
          message("Running differential expression analysis...")
          genes <- unique(cluster_obj$scMAGeCK_gene_assignment)
          genes <- genes[!genes %in% c("Non-Targeting", omit_genes)]

          de_res <- tryCatch(
            Run_wmvRegDE(
              object = cluster_obj,
              assay = "RNA",
              slot = "counts",
              labels = "scMAGeCK_gene_assignment",
              nt.class.name = "Non-Targeting",
              logfc.threshold = 0.1,
              PRTB_list = genes,
              full.results = TRUE,
              split.by = "experiment"
            ),
            error = function(e) {
              message(sprintf("DE analysis failed: %s", e$message))
              return(NULL)
            }
          )

          # Restructure results with metadata
          if (!is.null(de_res)) {
            # Extract all experiment information (dataset level)
            all_experiments <- unique(cluster_obj$experiment)
            message(sprintf("Total experiments in this %s cluster: %s",
                          current_modality, paste(all_experiments, collapse=", ")))

            # Get gene to experiment mapping from the experiment metadata
            gene_experiment_map <- table(cluster_obj$scMAGeCK_gene_assignment,
                                         cluster_obj$experiment)

            # Create a map between genes and their experiments
            gene_to_exps <- list()
            for (gene in rownames(gene_experiment_map)) {
              if (gene != "Non-Targeting" && !gene %in% omit_genes) {
                # Only include experiments where this gene has cells
                exps_with_gene <- colnames(gene_experiment_map)[gene_experiment_map[gene, ] > 0]
                gene_to_exps[[gene]] <- exps_with_gene
              }
            }

            # Initialize restructured results
            modified_de_res <- list()

            # Process each gene entry in the original results
            for (gene in names(de_res)) {
              if (gene != "metadata" && !is.null(de_res[[gene]])) {
                # Determine if this result has weighted p-values
                has_weighted <- FALSE
                if (is.data.frame(de_res[[gene]])) {
                  has_weighted <- any(grepl(":weight", colnames(de_res[[gene]]))) ||
                                  "p_weight" %in% colnames(de_res[[gene]])
                }

                # Get the experiments where this gene was tested
                gene_experiments <- gene_to_exps[[gene]]
                if (is.null(gene_experiments)) {
                  # Fallback: try to detect experiments from column names in results
                  if (is.data.frame(de_res[[gene]])) {
                    # Look for log2FC columns to detect experiments
                    fc_cols <- grep("^log2FC_", colnames(de_res[[gene]]), value = TRUE)
                    if (length(fc_cols) > 0) {
                      gene_experiments <- gsub("^log2FC_", "", fc_cols)
                    } else {
                      # If no log2FC columns, use all experiments as default
                      gene_experiments <- all_experiments
                    }
                  } else {
                    gene_experiments <- all_experiments
                  }
                }

                # Create metadata for this gene's results
                gene_metadata <- list(
                  gene = gene,
                  cluster = clust,
                  modality = current_modality,  # Add modality to metadata
                  experiments = gene_experiments,
                  is_weighted = has_weighted,
                  analysis_date = Sys.time(),
                  parameters = list(
                    neighbors = neighbors,
                    pca_dims = pca_dims,
                    logfc_threshold = logfc_threshold
                  )
                )

                # Create new structure with both results and metadata
                modified_de_res[[gene]] <- list(
                  results = de_res[[gene]],  # Original results
                  metadata = gene_metadata   # Gene-specific metadata
                )
              }
            }

            # Add global metadata to the entire result object
            modified_de_res$metadata <- list(
              cluster = clust,
              modality = current_modality,  # Add modality to global metadata
              n_perturbations = length(genes),
              experiments = all_experiments,
              analysis_date = Sys.time(),
              parameters = list(
                neighbors = neighbors,
                pca_dims = pca_dims,
                var_features = var_features,
                logfc_threshold = logfc_threshold,
                max_de_genes = max_de_genes
              ),
              stats = list(
                n_cells = ncol(cluster_obj),
                n_features = nrow(cluster_obj),
                n_nt_cells = as.numeric(nts),
                n_experiments = length(all_experiments)
              )
            )

            # Save the restructured results
            saveRDS(modified_de_res,
                    file.path(cluster_dir,
                              sprintf("%s_%s_clust_%s_mixscale_DEGs.rds",
                                      current_modality, phrase, clust)))
          }

          # Update analysis log
          analysis_log <- rbind(
            analysis_log,
            data.frame(
              modality = current_modality,
              cluster = clust,
              n_cells = ncol(cluster_obj),
              n_features = nrow(cluster_obj),
              n_nt_cells = nts,
              n_experiments = length(unique(cluster_obj$experiment)),
              timestamp = Sys.time()
            )
          )

        }, error = function(e) {
          message(sprintf("Error processing %s cluster %s: %s", current_modality, clust, e$message))
          cat(sprintf("Error processing %s cluster %s: %s\n",
                      current_modality, clust, e$message),
              file = log_file, append = TRUE)
        })
      }

      # Save analysis log for this modality
      write.csv(analysis_log,
                file.path(modality_dir, sprintf("%s_analysis_log.csv", current_modality)),
                row.names = FALSE)

      # Log completion
      cat(sprintf("\n%s Analysis completed: %s\n", current_modality, Sys.time()),
          file = log_file, append = TRUE)
    }

    # Create a summary log for all modalities
    summary_log <- file.path(base_output_dir, "analysis_summary.txt")
    cat(sprintf("MixScale analysis completed for modalities: %s\n", 
                paste(modalities_to_process, collapse = ", ")), 
        file = summary_log)
    cat(sprintf("Completed at: %s\n", Sys.time()), file = summary_log, append = TRUE)
  }



  # Return results summary
  return(list(
    output_directory = output_dir,
    modality = modality,
    completed_at = Sys.time()
  ))
}
