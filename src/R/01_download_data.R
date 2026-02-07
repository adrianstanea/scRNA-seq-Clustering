# =============================================================================
# 01_download_data.R - Dataset Acquisition
# scRNA-seq Clustering Comparative Study
# =============================================================================

source("src/R/00_setup.R")

# -----------------------------------------------------------------------------
# Download scMixology Dataset
# Ground truth: 3 cell lines (HCC827, H1975, H2228)
# -----------------------------------------------------------------------------

download_scmixology <- function(output_dir = "data/raw") {
message("Downloading scMixology dataset...")

output_file <- file.path(output_dir, "scMixology.RData")

if (file.exists(output_file)) {
  message("scMixology already exists at: ", output_file)
  return(output_file)
}

options(timeout = 600)
url <- "https://github.com/LuyiTian/sc_mixology/raw/master/data/sincell_with_class.RData"

tryCatch({
  download.file(url, destfile = output_file, mode = "wb")
  message("Downloaded to: ", output_file)
  output_file
}, error = function(e) {
  warning("Failed to download scMixology: ", e$message)
  NULL
})
}

# -----------------------------------------------------------------------------
# Download PBMC 3k Dataset
# Community standard: ~9 immune cell types
# -----------------------------------------------------------------------------

download_pbmc3k <- function(output_dir = "data/raw") {
message("Downloading PBMC 3k dataset...")

output_file <- file.path(output_dir, "pbmc3k.rds")

if (file.exists(output_file)) {
  message("PBMC 3k already exists at: ", output_file)
  return(output_file)
}

# Load SeuratData
if (!requireNamespace("SeuratData", quietly = TRUE)) {
  remotes::install_github("satijalab/seurat-data")
}
library(SeuratData)

# Install and load the dataset
tryCatch({
  InstallData("pbmc3k")
  data("pbmc3k")

  # Save to file
  saveRDS(pbmc3k, output_file)
  message("Saved to: ", output_file)
  output_file
}, error = function(e) {
  warning("Failed to download PBMC 3k: ", e$message)
  NULL
})
}

# -----------------------------------------------------------------------------
# Download Zeisel Brain Dataset
# Complex hierarchical structure: Mouse cortex/hippocampus
# -----------------------------------------------------------------------------

download_zeisel <- function(output_dir = "data/raw") {
message("Downloading Zeisel Brain dataset...")

output_file <- file.path(output_dir, "zeisel_brain.rds")

if (file.exists(output_file)) {
  message("Zeisel Brain already exists at: ", output_file)
  return(output_file)
}

# Load scRNAseq package
if (!requireNamespace("scRNAseq", quietly = TRUE)) {
  BiocManager::install("scRNAseq")
}
library(scRNAseq)

tryCatch({
  sce <- ZeiselBrainData()

  # Convert to Seurat for consistency
  library(Seurat)
  sobj <- CreateSeuratObject(
    counts = counts(sce),
    meta.data = as.data.frame(colData(sce))
  )

  saveRDS(sobj, output_file)
  message("Saved to: ", output_file)
  output_file
}, error = function(e) {
  warning("Failed to download Zeisel Brain: ", e$message)
  NULL
})
}

# -----------------------------------------------------------------------------
# Download All Datasets
# -----------------------------------------------------------------------------
download_all_datasets <- function(output_dir = "data/raw") {
# Create output directory if needed
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

results <- list(
  scmixology = download_scmixology(output_dir),
  pbmc3k = download_pbmc3k(output_dir),
  zeisel = download_zeisel(output_dir)
)

# Summary
message("\n=== Download Summary ===")
for (name in names(results)) {
  status <- if (!is.null(results[[name]])) "SUCCESS" else "FAILED"
  message(sprintf("%-15s: %s", name, status))
}

results
}

# -----------------------------------------------------------------------------
# Load Downloaded Datasets
# -----------------------------------------------------------------------------

load_dataset <- function(dataset_name, data_dir = "data/raw") {
file_map <- list(
  scmixology = "scMixology.RData",
  pbmc3k = "pbmc3k.rds",
  zeisel = "zeisel_brain.rds"
)

if (!dataset_name %in% names(file_map)) {
  stop("Unknown dataset: ", dataset_name,
       ". Available: ", paste(names(file_map), collapse = ", "))
}

file_path <- file.path(data_dir, file_map[[dataset_name]])

if (!file.exists(file_path)) {
  stop("Dataset not found at: ", file_path, ". Run download_all_datasets() first.")
}

if (dataset_name == "scmixology") {
  # Load RData file
  env <- new.env()
  load(file_path, envir = env)
  # Return the 10x data by default
  # Check available objects and return the appropriate one
  available_objects <- ls(envir = env)
  message("Available objects in scMixology: ", paste(available_objects, collapse = ", "))

  if ("sce10x_qc" %in% available_objects) {
    return(env$sce10x_qc)
  } else if ("sce_sc_10x_qc" %in% available_objects) {
    return(env$sce_sc_10x_qc)
  } else if (any(grepl("10x", available_objects, ignore.case = TRUE))) {
    # Try to find any 10x object
    obj_name <- available_objects[grepl("10x", available_objects, ignore.case = TRUE)][1]
    return(env[[obj_name]])
  } else if (length(available_objects) > 0) {
    # Return the first object
    message("Using first available object: ", available_objects[1])
    return(env[[available_objects[1]]])
  } else {
    stop("No data objects found in scMixology RData file")
  }
} else {
  return(readRDS(file_path))
}
}

# -----------------------------------------------------------------------------
# Main execution
# -----------------------------------------------------------------------------

if (sys.nframe() == 0) {
# Script is being run directly
message("Starting dataset download...")
download_all_datasets()
}
