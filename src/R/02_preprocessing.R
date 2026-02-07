# =============================================================================
# 02_preprocessing.R - Standardized Preprocessing Pipeline
# scRNA-seq Clustering Comparative Study
# =============================================================================

source("src/R/00_setup.R")

# -----------------------------------------------------------------------------
# Main Preprocessing Function
# -----------------------------------------------------------------------------

preprocess_seurat <- function(
sobj,
min_genes = 200,
max_genes = 5000,
max_pct_mt = 5,
min_cells = 3,
n_features = 2000,
n_pcs = 30,
scale_factor = 10000,
seed = 42
) {
set.seed(seed)
message("Starting preprocessing pipeline...")

# Record initial cell count
n_initial <- ncol(sobj)
message(sprintf("Initial cells: %d", n_initial))

# 1. Calculate mitochondrial percentage if not present
if (!"percent.mt" %in% colnames(sobj@meta.data)) {
  message("Calculating mitochondrial percentage...")
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-|^mt-")
}

# 2. Quality Control Filtering
message("Applying QC filters...")
sobj <- subset(
  sobj,
  subset = nFeature_RNA > min_genes &
           nFeature_RNA < max_genes &
           percent.mt < max_pct_mt
)
n_after_qc <- ncol(sobj)
message(sprintf("Cells after QC: %d (removed %d)", n_after_qc, n_initial - n_after_qc))

# 3. Filter genes by minimum cells
message("Filtering genes...")
gene_counts <- rowSums(GetAssayData(sobj, layer = "counts") > 0)
genes_keep <- names(gene_counts[gene_counts >= min_cells])
sobj <- subset(sobj, features = genes_keep)
message(sprintf("Genes retained: %d", nrow(sobj)))

# 4. Normalization
message("Normalizing data (LogNormalize)...")
sobj <- NormalizeData(
  sobj,
  normalization.method = "LogNormalize",
  scale.factor = scale_factor,
  verbose = FALSE
)

# 5. Feature Selection (HVGs)
message(sprintf("Selecting top %d variable features (VST)...", n_features))
sobj <- FindVariableFeatures(
  sobj,
  selection.method = "vst",
  nfeatures = n_features,
  verbose = FALSE
)

# 6. Scaling
message("Scaling data...")
sobj <- ScaleData(sobj, verbose = FALSE)

# 7. PCA
message(sprintf("Running PCA (n_pcs = %d)...", n_pcs))
sobj <- RunPCA(
  sobj,
  npcs = n_pcs,
  verbose = FALSE,
  seed.use = seed
)

# 8. UMAP for visualization
message("Computing UMAP...")
sobj <- RunUMAP(
  sobj,
  dims = 1:n_pcs,
  verbose = FALSE,
  seed.use = seed
)

# 9. Build neighbor graph (for clustering methods)
message("Building neighbor graph...")
sobj <- FindNeighbors(
  sobj,
  dims = 1:n_pcs,
  verbose = FALSE
)

message("Preprocessing complete.")
return(sobj)
}

# -----------------------------------------------------------------------------
# Preprocessing from Config
# -----------------------------------------------------------------------------

preprocess_from_config <- function(sobj, config_path = "config/preprocessing_params.yaml") {
config <- yaml::read_yaml(config_path)

preprocess_seurat(
  sobj,
  min_genes = config$qc$min_genes_per_cell,
  max_genes = config$qc$max_genes_per_cell,
  max_pct_mt = config$qc$max_pct_mt,
  min_cells = config$qc$min_cells_per_gene,
  n_features = config$feature_selection$n_top_genes,
  n_pcs = config$dimensionality_reduction$n_pcs,
  scale_factor = config$normalization$scale_factor,
  seed = config$random_seed
)
}

# -----------------------------------------------------------------------------
# Preprocess All Datasets
# -----------------------------------------------------------------------------

preprocess_all_datasets <- function(
input_dir = "data/raw",
output_dir = "data/processed",
config_path = "config/preprocessing_params.yaml"
) {
source("src/R/01_download_data.R")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

datasets <- c("pbmc3k", "zeisel")  # scmixology handled separately

results <- list()

for (dataset in datasets) {
  message("\n", strrep("=", 60))
  message("Processing: ", dataset)
  message(strrep("=", 60))

  tryCatch({
    # Load raw data
    sobj <- load_dataset(dataset, input_dir)

    # Convert SingleCellExperiment to Seurat if needed
    if (inherits(sobj, "SingleCellExperiment")) {
      suppressPackageStartupMessages(library(SingleCellExperiment))
      count_matrix <- SummarizedExperiment::assay(sobj, "counts")
      sobj <- CreateSeuratObject(
        counts = count_matrix,
        meta.data = as.data.frame(colData(sobj))
      )
    }

    # Ensure Seurat v5 format compatibility
    if (inherits(sobj, "Seurat")) {
      # Update object if needed
      tryCatch({
        sobj <- UpdateSeuratObject(sobj)
      }, error = function(e) NULL)
    }

    # Set ground truth column based on dataset
    if (dataset == "pbmc3k") {
      # PBMC 3k may have seurat_annotations or cell type labels
      possible_cols <- c("seurat_annotations", "celltype", "cell_type",
                         "CellType", "cluster", "louvain")
      for (col in possible_cols) {
        if (col %in% colnames(sobj@meta.data)) {
          sobj$ground_truth <- sobj@meta.data[[col]]
          message("Ground truth column set from: ", col)
          break
        }
      }
    } else if (dataset == "zeisel") {
      # Zeisel brain data has level1class or similar
      possible_cols <- c("level1class", "level2class", "cell_type1",
                         "tissue", "CellType", "cell_type")
      for (col in possible_cols) {
        if (col %in% colnames(sobj@meta.data)) {
          sobj$ground_truth <- sobj@meta.data[[col]]
          message("Ground truth column set from: ", col)
          break
        }
      }
    }

    # List available metadata columns
    message("Available metadata columns: ", paste(colnames(sobj@meta.data), collapse = ", "))

    # Preprocess
    sobj <- preprocess_from_config(sobj, config_path)

    # Save
    output_file <- file.path(output_dir, paste0(dataset, "_processed.rds"))
    saveRDS(sobj, output_file)
    message("Saved to: ", output_file)

    results[[dataset]] <- output_file
  }, error = function(e) {
    warning("Failed to process ", dataset, ": ", e$message)
    results[[dataset]] <- NULL
  })
}

# Handle scMixology separately (multiple platforms)
message("\n", strrep("=", 60))
message("Processing: scMixology (10x)")
message(strrep("=", 60))

tryCatch({
  # Load SingleCellExperiment for counts() function
  suppressPackageStartupMessages(library(SingleCellExperiment))

  env <- new.env()
  load(file.path(input_dir, "scMixology.RData"), envir = env)

  # Find available objects
  available_objects <- ls(envir = env)
  message("Available scMixology objects: ", paste(available_objects, collapse = ", "))

  # Use 10x data - find the right object
  sce <- NULL
  if ("sce10x_qc" %in% available_objects) {
    sce <- env$sce10x_qc
  } else if ("sce_sc_10x_qc" %in% available_objects) {
    sce <- env$sce_sc_10x_qc
  } else if (any(grepl("10x", available_objects, ignore.case = TRUE))) {
    obj_name <- available_objects[grepl("10x", available_objects, ignore.case = TRUE)][1]
    sce <- env[[obj_name]]
  } else if (length(available_objects) > 0) {
    sce <- env[[available_objects[1]]]
  }

  if (is.null(sce)) {
    stop("Could not find scMixology data object")
  }

  # Convert SCE to Seurat - use SummarizedExperiment::assay for robustness
  count_matrix <- SummarizedExperiment::assay(sce, "counts")
  sobj <- CreateSeuratObject(
    counts = count_matrix,
    meta.data = as.data.frame(colData(sce))
  )

  # List available metadata
  message("Available metadata columns: ", paste(colnames(sobj@meta.data), collapse = ", "))

  # Store ground truth - try multiple possible column names
  gt_cols <- c("cell_line", "cell_line_demuxlet", "cellline", "CellLine",
               "cell.line", "sample", "group")
  gt_found <- FALSE
  for (col in gt_cols) {
    if (col %in% colnames(sobj@meta.data)) {
      sobj$ground_truth <- sobj@meta.data[[col]]
      message("Ground truth column set from: ", col)
      gt_found <- TRUE
      break
    }
  }
  if (!gt_found) {
    warning("No ground truth column found for scMixology. Available: ",
            paste(colnames(sobj@meta.data), collapse = ", "))
  }

  # scMixology has higher gene counts, use relaxed thresholds
  sobj <- preprocess_seurat(
    sobj,
    min_genes = 500,      # More lenient for this dataset
    max_genes = 15000,    # scMixology cells have 6000-11000 genes
    max_pct_mt = 20,      # No MT genes in this dataset anyway
    min_cells = 3,
    n_features = 2000,
    n_pcs = 30
  )

  output_file <- file.path(output_dir, "scmixology_processed.rds")
  saveRDS(sobj, output_file)
  message("Saved to: ", output_file)

  results[["scmixology"]] <- output_file
}, error = function(e) {
  warning("Failed to process scMixology: ", e$message)
  results[["scmixology"]] <- NULL
})

message("\n=== Preprocessing Summary ===")
for (name in names(results)) {
  status <- if (!is.null(results[[name]])) "SUCCESS" else "FAILED"
  message(sprintf("%-15s: %s", name, status))
}

results
}

# -----------------------------------------------------------------------------
# QC Visualization
# -----------------------------------------------------------------------------

plot_qc <- function(sobj, title = "") {
p1 <- VlnPlot(sobj, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
p2 <- VlnPlot(sobj, features = "nCount_RNA", pt.size = 0.1) + NoLegend()
p3 <- VlnPlot(sobj, features = "percent.mt", pt.size = 0.1) + NoLegend()

(p1 | p2 | p3) + plot_annotation(title = title)
}

# -----------------------------------------------------------------------------
# Main execution
# -----------------------------------------------------------------------------

if (sys.nframe() == 0) {
setup_environment()
preprocess_all_datasets()
}
