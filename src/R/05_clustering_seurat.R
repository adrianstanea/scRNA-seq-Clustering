# =============================================================================
# 05_clustering_seurat.R - Baseline Seurat/Leiden Clustering
# Standard graph-based clustering at multiple resolutions
# =============================================================================

source("src/R/00_setup.R")

# -----------------------------------------------------------------------------
# Seurat Leiden Clustering Function
# -----------------------------------------------------------------------------

run_seurat_clustering <- function(
sobj,
resolutions = c(0.2, 0.5, 0.8, 1.0, 1.5),
algorithm = 4,  # 4 = Leiden
seed = 42,
verbose = TRUE
) {
message("Running Seurat Leiden clustering...")
message(sprintf("  - Resolutions: %s", paste(resolutions, collapse = ", ")))
message(sprintf("  - Random seed: %d", seed))

start_time <- Sys.time()

# Run clustering at multiple resolutions
for (res in resolutions) {
  if (verbose) message(sprintf("  Clustering at resolution %.1f...", res))

  sobj <- FindClusters(
    sobj,
    resolution = res,
    algorithm = algorithm,
    random.seed = seed,
    verbose = FALSE
  )

  # Rename cluster column for clarity
  col_name <- paste0("leiden_res", gsub("\\.", "_", as.character(res)))
  sobj@meta.data[[col_name]] <- sobj$seurat_clusters
}

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

# Store runtime
sobj@misc$seurat_runtime <- runtime
message(sprintf("Seurat clustering completed in %.1f seconds", runtime))

# Report cluster counts
message("Cluster counts by resolution:")
for (res in resolutions) {
  col_name <- paste0("leiden_res", gsub("\\.", "_", as.character(res)))
  n_clusters <- length(unique(sobj@meta.data[[col_name]]))
  message(sprintf("  Resolution %.1f: %d clusters", res, n_clusters))
}

return(sobj)
}

# -----------------------------------------------------------------------------
# Run Seurat from Config
# -----------------------------------------------------------------------------

run_seurat_from_config <- function(sobj, config_path = "config/clustering_params.yaml") {
config <- yaml::read_yaml(config_path)

run_seurat_clustering(
  sobj,
  resolutions = config$leiden$resolutions,
  seed = config$random_seed
)
}

# -----------------------------------------------------------------------------
# Process All Datasets with Seurat
# -----------------------------------------------------------------------------

run_seurat_all_datasets <- function(
input_dir = "data/processed",
output_dir = "data/results",
config_path = "config/clustering_params.yaml"
) {
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

config <- yaml::read_yaml(config_path)
resolutions <- config$leiden$resolutions

datasets <- c("scmixology", "pbmc3k", "zeisel")
results <- list()

for (dataset in datasets) {
  message("\n", strrep("=", 60))
  message("Seurat Leiden clustering: ", dataset)
  message(strrep("=", 60))

  input_file <- file.path(input_dir, paste0(dataset, "_processed.rds"))

  if (!file.exists(input_file)) {
    warning("Processed file not found: ", input_file)
    next
  }

  tryCatch({
    sobj <- readRDS(input_file)
    sobj <- run_seurat_from_config(sobj, config_path)

    # Save results
    output_file <- file.path(output_dir, paste0(dataset, "_seurat.rds"))
    saveRDS(sobj, output_file)
    message("Saved to: ", output_file)

    # Store summary - cluster counts at each resolution
    cluster_counts <- sapply(resolutions, function(res) {
      col_name <- paste0("leiden_res", gsub("\\.", "_", as.character(res)))
      length(unique(sobj@meta.data[[col_name]]))
    })
    names(cluster_counts) <- paste0("res_", resolutions)

    results[[dataset]] <- list(
      cluster_counts = cluster_counts,
      runtime = sobj@misc$seurat_runtime,
      file = output_file
    )
  }, error = function(e) {
    warning("Seurat failed for ", dataset, ": ", e$message)
    results[[dataset]] <- NULL
  })
}

# Summary
message("\n=== Seurat Clustering Summary ===")
for (name in names(results)) {
  if (!is.null(results[[name]])) {
    message(sprintf("%-15s: %s (%.1f sec)",
                    name,
                    paste(results[[name]]$cluster_counts, collapse = "/"),
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
run_seurat_all_datasets()
}
