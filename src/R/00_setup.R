# =============================================================================
# 00_setup.R - Environment Setup and Package Loading
# scRNA-seq Clustering Comparative Study
# =============================================================================

# Disable renv in Docker container - use system library instead
if (file.exists("/.dockerenv") || Sys.getenv("DOCKER_CONTAINER") == "true") {
    # Deactivate renv if active
    if ("renv" %in% loadedNamespaces()) {
        try(renv::deactivate(), silent = TRUE)
    }
    # Use system library
    .libPaths("/usr/local/lib/R/site-library")
}

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# -----------------------------------------------------------------------------
# Required Packages
# -----------------------------------------------------------------------------

required_packages <- c(
  "Seurat",
  "ggplot2",
  "dplyr",
  "tidyr",
  "patchwork",
  "yaml",
  "future",
  "RColorBrewer",
  "viridis",
  "pheatmap"
)

bioc_packages <- c(
  "SingleCellExperiment",
  "scRNAseq",
  "scater",
  "scran"
)

github_packages <- list(
  presto = "immunogenomics/presto",
  SeuratData = "satijalab/seurat-data"
)

# -----------------------------------------------------------------------------
# Installation Functions
# -----------------------------------------------------------------------------

install_if_missing <- function(packages) {
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) {
  message("Installing: ", paste(new_packages, collapse = ", "))
  install.packages(new_packages)
}
}

install_bioc_if_missing <- function(packages) {
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) {
  message("Installing Bioconductor packages: ", paste(new_packages, collapse = ", "))
  BiocManager::install(new_packages, ask = FALSE)
}
}

install_github_if_missing <- function(packages) {
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
for (pkg_name in names(packages)) {
  if (!(pkg_name %in% installed.packages()[, "Package"])) {
    message("Installing from GitHub: ", packages[[pkg_name]])
    remotes::install_github(packages[[pkg_name]],
                            repos = BiocManager::repositories(),
                            upgrade = "never")
  }
}
}

# -----------------------------------------------------------------------------
# Run Installation
# -----------------------------------------------------------------------------

setup_environment <- function(install_missing = TRUE) {
if (install_missing) {
  message("Checking and installing required packages...")
  install_if_missing(required_packages)
  install_bioc_if_missing(bioc_packages)
  install_github_if_missing(github_packages)
}

# Load core packages
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(yaml)
})

message("Environment setup complete.")
}

# -----------------------------------------------------------------------------
# Configuration Loading
# -----------------------------------------------------------------------------

load_config <- function(config_dir = "config") {
preprocessing <- yaml::read_yaml(file.path(config_dir, "preprocessing_params.yaml"))
clustering <- yaml::read_yaml(file.path(config_dir, "clustering_params.yaml"))

list(
  preprocessing = preprocessing,
  clustering = clustering
)
}

# -----------------------------------------------------------------------------
# Reproducibility Settings
# -----------------------------------------------------------------------------

set_seed <- function(seed = 42) {
set.seed(seed)
message("Random seed set to: ", seed)
}

# -----------------------------------------------------------------------------
# Parallel Processing Setup
# -----------------------------------------------------------------------------

setup_parallel <- function(n_cores = NULL) {
if (is.null(n_cores)) {
  n_cores <- max(1, parallel::detectCores() - 2)
}
future::plan("multicore", workers = n_cores)
message("Parallel processing enabled with ", n_cores, " cores")
}

# -----------------------------------------------------------------------------
# Project Paths
# -----------------------------------------------------------------------------

get_paths <- function(project_root = ".") {
list(
  root = project_root,
  data_raw = file.path(project_root, "data", "raw"),
  data_processed = file.path(project_root, "data", "processed"),
  data_results = file.path(project_root, "data", "results"),
  src_r = file.path(project_root, "src", "R"),
  config = file.path(project_root, "config"),
  figures = file.path(project_root, "docs", "latex", "figures")
)
}

# -----------------------------------------------------------------------------
# Auto-run if sourced directly
# -----------------------------------------------------------------------------

if (interactive()) {
message("Run setup_environment() to install packages and load libraries.")
message("Run load_config() to load configuration parameters.")
}
