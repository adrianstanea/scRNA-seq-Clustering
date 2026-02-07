# =============================================================================
# 03_clustering_sc3.R - SC3 Consensus Clustering
# Single-Cell Consensus Clustering (replaces CHOIR)
# =============================================================================

source("src/R/00_setup.R")

# -----------------------------------------------------------------------------
# SC3 Clustering Function
# -----------------------------------------------------------------------------

run_sc3_clustering <- function(
    sobj,
    k_range = NULL,
    n_cores = 8,
    biology = TRUE,
    verbose = TRUE
) {
    # Load required packages
    if (!requireNamespace("SC3", quietly = TRUE)) {
        stop("SC3 not installed. Run: BiocManager::install('SC3')")
    }
    library(SC3)
    library(SingleCellExperiment)

    message("Running SC3 consensus clustering...")

    start_time <- Sys.time()

    # Convert Seurat to SingleCellExperiment
    # Get counts matrix directly from Seurat
    counts_mat <- GetAssayData(sobj, layer = "counts")

    # Get normalized data
    data_mat <- tryCatch({
        GetAssayData(sobj, layer = "data")
    }, error = function(e) {
        log1p(counts_mat)
    })

    # SC3 requires dense matrices for its internal operations
    # Convert sparse to dense (this may use more memory but is required)
    message("Converting to dense matrices for SC3...")
    counts_mat <- as.matrix(counts_mat)
    data_mat <- as.matrix(data_mat)

    # Create SCE with explicit assays
    sce <- SingleCellExperiment(
        assays = list(
            counts = counts_mat,
            logcounts = data_mat
        ),
        colData = sobj@meta.data
    )

    # Set rowData for gene filtering
    rowData(sce)$feature_symbol <- rownames(sce)

    # Determine k range if not provided
    if (is.null(k_range)) {
        # Estimate reasonable k range based on data size
        n_cells <- ncol(sce)
        # Use a more conservative range that always works
        k_min <- max(2, min(5, round(sqrt(n_cells) / 4)))
        k_max <- min(15, max(k_min + 3, round(sqrt(n_cells) / 2)))
        # Ensure k_min <= k_max
        if (k_min > k_max) {
            k_min <- 2
            k_max <- max(k_min + 1, min(10, round(sqrt(n_cells))))
        }
        k_range <- seq(k_min, k_max, by = 1)
    }

    message(sprintf("  - k range: %d to %d", min(k_range), max(k_range)))
    message(sprintf("  - Cores: %d", n_cores))
    message(sprintf("  - Biology: %s", biology))

    # Prepare SC3 - compute distance matrices and transformations
    message("Setting SC3 parameters...")
    # svm_num_cells must be strictly less than number of cells
    n_cells <- ncol(sce)
    svm_cells <- min(5000, floor(n_cells * 0.8))  # Use 80% of cells max
    if (svm_cells >= n_cells) {
        svm_cells <- max(100, n_cells - 10)  # At least 10 less than total
    }
    sce <- sc3_prepare(sce, gene_filter = TRUE, n_cores = n_cores,
                       svm_num_cells = svm_cells)

    # Estimate optimal k
    message("Estimating optimal k...")
    sce <- sc3_estimate_k(sce)
    k_estimated <- metadata(sce)$sc3$k_estimation
    if (is.null(k_estimated) || is.na(k_estimated)) {
        # Fallback to middle of k_range if estimation fails
        k_estimated <- round(median(k_range))
        message(sprintf("  - k estimation failed, using median of range: %d", k_estimated))
    } else {
        message(sprintf("  - Estimated optimal k: %d", k_estimated))
    }

    # Run SC3 clustering for the estimated k
    # Also run for a range around the estimated k for comparison
    k_to_run <- unique(c(k_estimated, k_range[k_range >= 2]))
    k_to_run <- k_to_run[k_to_run <= 20]  # Cap at 20 clusters
    k_to_run <- sort(k_to_run)

    message(sprintf("Running SC3 for k values: %s", paste(k_to_run, collapse=", ")))
    sce <- sc3(sce, ks = k_to_run, biology = biology, n_cores = n_cores)

    end_time <- Sys.time()
    runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # Extract cluster assignments at estimated k
    cluster_col <- paste0("sc3_", k_estimated, "_clusters")
    if (cluster_col %in% colnames(colData(sce))) {
        sc3_clusters <- colData(sce)[[cluster_col]]
    } else {
        # Fall back to first available
        sc3_cols <- grep("^sc3_\\d+_clusters$", colnames(colData(sce)), value = TRUE)
        if (length(sc3_cols) > 0) {
            sc3_clusters <- colData(sce)[[sc3_cols[1]]]
            cluster_col <- sc3_cols[1]
        } else {
            stop("No SC3 cluster columns found")
        }
    }

    # Add SC3 clusters back to Seurat object
    sobj$SC3_clusters <- sc3_clusters
    sobj$SC3_k_estimated <- k_estimated

    # Store additional metadata
    sobj@misc$sc3_runtime <- runtime
    sobj@misc$sc3_k_estimated <- k_estimated
    sobj@misc$sc3_sce <- sce  # Store SCE for later analysis

    n_clusters <- length(unique(sc3_clusters))
    message(sprintf("SC3 completed in %.1f seconds", runtime))
    message(sprintf("Found %d clusters (k=%d)", n_clusters, k_estimated))

    return(sobj)
}

# -----------------------------------------------------------------------------
# Run SC3 from Config
# -----------------------------------------------------------------------------

run_sc3_from_config <- function(sobj, config_path = "config/clustering_params.yaml") {
    config <- yaml::read_yaml(config_path)

    # SC3 config (use defaults if not specified)
    k_range <- if (!is.null(config$sc3$k_range)) config$sc3$k_range else NULL
    n_cores <- if (!is.null(config$sc3$n_cores)) config$sc3$n_cores else 8
    biology <- if (!is.null(config$sc3$biology)) config$sc3$biology else TRUE

    run_sc3_clustering(
        sobj,
        k_range = k_range,
        n_cores = n_cores,
        biology = biology
    )
}

# -----------------------------------------------------------------------------
# Process All Datasets with SC3
# -----------------------------------------------------------------------------

run_sc3_all_datasets <- function(
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
        message("SC3 clustering: ", dataset)
        message(strrep("=", 60))

        input_file <- file.path(input_dir, paste0(dataset, "_processed.rds"))

        if (!file.exists(input_file)) {
            warning("Processed file not found: ", input_file)
            next
        }

        tryCatch({
            sobj <- readRDS(input_file)
            sobj <- run_sc3_from_config(sobj, config_path)

            # Save results
            output_file <- file.path(output_dir, paste0(dataset, "_sc3.rds"))
            saveRDS(sobj, output_file)
            message("Saved to: ", output_file)

            # Store summary
            results[[dataset]] <- list(
                n_clusters = length(unique(sobj$SC3_clusters)),
                k_estimated = sobj@misc$sc3_k_estimated,
                runtime = sobj@misc$sc3_runtime,
                file = output_file
            )
        }, error = function(e) {
            warning("SC3 failed for ", dataset, ": ", e$message)
            results[[dataset]] <- NULL
        })
    }

    # Summary
    message("\n=== SC3 Clustering Summary ===")
    for (name in names(results)) {
        if (!is.null(results[[name]])) {
            message(sprintf("%-15s: %d clusters (k=%d, %.1f sec)",
                            name,
                            results[[name]]$n_clusters,
                            results[[name]]$k_estimated,
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
    run_sc3_all_datasets()
}
