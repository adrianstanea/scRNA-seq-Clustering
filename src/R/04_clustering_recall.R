# =============================================================================
# 04_clustering_recall.R - recall (formerly callback) Clustering
# Knockoff-based Calibration for FDR Control
# =============================================================================

source("src/R/00_setup.R")

# -----------------------------------------------------------------------------
# recall Clustering Function
# -----------------------------------------------------------------------------

run_recall_clustering <- function(
sobj,
dims = 30,
algorithm = "louvain",
resolution_start = 0.8,
num_clusters_start = 20,
null_method = "ZIP",
n_cores = 4,
verbose = TRUE
) {
# Load recall
if (!requireNamespace("recall", quietly = TRUE)) {
  stop("recall not installed. Run: devtools::install_github('lcrawlab/recall')")
}
library(recall)

# Load presto for speed if available
if (requireNamespace("presto", quietly = TRUE)) {
  library(presto)
  message("presto loaded for accelerated DE testing")
}

message("Running recall clustering...")
message(sprintf("  - Dimensions: 1:%d", dims))
message(sprintf("  - Algorithm: %s", algorithm))
message(sprintf("  - Resolution start: %.2f", resolution_start))
message(sprintf("  - Null method: %s", null_method))

start_time <- Sys.time()

# Run recall using the correct API
# FindClustersRecall(seurat_obj, resolution_start, reduction_percentage,
#                    num_clusters_start, dims, algorithm, null_method,
#                    assay, cores, shared_memory_max, verbose)
sobj <- recall::FindClustersRecall(
  seurat_obj = sobj,
  resolution_start = resolution_start,
  num_clusters_start = num_clusters_start,
  dims = 1:dims,
  algorithm = algorithm,
  null_method = null_method,
  assay = "RNA",
  cores = n_cores,
  verbose = verbose
)

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

# Store runtime in metadata
sobj@misc$recall_runtime <- runtime
message(sprintf("recall completed in %.1f seconds", runtime))

# Find the recall cluster column and standardize it
recall_cols <- grep("recall|Recall", colnames(sobj@meta.data), value = TRUE)
if (length(recall_cols) > 0) {
  sobj$recall_clusters <- sobj@meta.data[[recall_cols[1]]]
  n_clusters <- length(unique(sobj$recall_clusters))
  message(sprintf("Found %d FDR-controlled clusters", n_clusters))
}

return(sobj)
}

# -----------------------------------------------------------------------------
# Run recall from Config
# -----------------------------------------------------------------------------

run_recall_from_config <- function(sobj, config_path = "config/clustering_params.yaml") {
config <- yaml::read_yaml(config_path)

run_recall_clustering(
  sobj,
  dims = if (!is.null(config$recall$dims)) config$recall$dims else 30,
  algorithm = if (!is.null(config$recall$algorithm)) config$recall$algorithm else "louvain",
  resolution_start = if (!is.null(config$recall$resolution_start)) config$recall$resolution_start else 0.8,
  num_clusters_start = if (!is.null(config$recall$num_clusters_start)) config$recall$num_clusters_start else 20,
  null_method = if (!is.null(config$recall$null_method)) config$recall$null_method else "ZIP",
  n_cores = if (!is.null(config$recall$n_cores)) config$recall$n_cores else 4
)
}

# -----------------------------------------------------------------------------
# Process All Datasets with recall
# -----------------------------------------------------------------------------

run_recall_all_datasets <- function(
input_dir = "data/processed",
output_dir = "data/results",
config_path = "config/clustering_params.yaml"
) {
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

datasets <- c("scmixology", "pbmc3k", "zeisel")
results <- list()

for (dataset in datasets) {
  message("\n", strrep("=", 60))
  message("recall clustering: ", dataset)
  message(strrep("=", 60))

  input_file <- file.path(input_dir, paste0(dataset, "_processed.rds"))

  if (!file.exists(input_file)) {
    warning("Processed file not found: ", input_file)
    next
  }

  tryCatch({
    sobj <- readRDS(input_file)
    sobj <- run_recall_from_config(sobj, config_path)

    # Save results
    output_file <- file.path(output_dir, paste0(dataset, "_recall.rds"))
    saveRDS(sobj, output_file)
    message("Saved to: ", output_file)

    # Store summary
    results[[dataset]] <- list(
      n_clusters = length(unique(sobj$recall_clusters)),
      runtime = sobj@misc$recall_runtime,
      file = output_file
    )
  }, error = function(e) {
    warning("recall failed for ", dataset, ": ", e$message)
    results[[dataset]] <- NULL
  })
}

# Summary
message("\n=== recall Clustering Summary ===")
for (name in names(results)) {
  if (!is.null(results[[name]])) {
    message(sprintf("%-15s: %d clusters (%.1f sec)",
                    name,
                    results[[name]]$n_clusters,
                    results[[name]]$runtime))
  } else {
    message(sprintf("%-15s: FAILED", name))
  }
}

results
}

# -----------------------------------------------------------------------------
# Main execution
# -----------------------------------------------------------------------------

if (sys.nframe() == 0) {
setup_environment()
run_recall_all_datasets()
}
