# =============================================================================
# run_pipeline.R - Master Pipeline Script
# scRNA-seq Clustering Comparative Study
#
# This script orchestrates the full analysis pipeline from data download
# through clustering, evaluation, and visualization.
#
# Usage:
#   Rscript src/R/run_pipeline.R [--skip-download] [--skip-preprocessing]
#                                [--methods=seurat,choir,recall] [--datasets=all]
# =============================================================================

# -----------------------------------------------------------------------------
# Parse Command Line Arguments
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

# Default options
SKIP_DOWNLOAD <- "--skip-download" %in% args
SKIP_PREPROCESSING <- "--skip-preprocessing" %in% args
SKIP_CLUSTERING <- "--skip-clustering" %in% args
SKIP_EVALUATION <- "--skip-evaluation" %in% args
SKIP_VISUALIZATION <- "--skip-visualization" %in% args
DRY_RUN <- "--dry-run" %in% args
VERBOSE <- !("--quiet" %in% args)

# Parse methods argument
methods_arg <- args[grepl("--methods=", args)]
if (length(methods_arg) > 0) {
  METHODS <- unlist(strsplit(gsub("--methods=", "", methods_arg), ","))
} else {
  METHODS <- c("seurat", "choir", "recall")
}

# Parse datasets argument
datasets_arg <- args[grepl("--datasets=", args)]
if (length(datasets_arg) > 0) {
  datasets_str <- gsub("--datasets=", "", datasets_arg)
  if (datasets_str == "all") {
    DATASETS <- c("scmixology", "pbmc3k", "zeisel")
  } else {
    DATASETS <- unlist(strsplit(datasets_str, ","))
  }
} else {
  DATASETS <- c("scmixology", "pbmc3k", "zeisel")
}

# -----------------------------------------------------------------------------
# Utility Functions
# -----------------------------------------------------------------------------

log_step <- function(msg, level = "INFO") {
  if (VERBOSE) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] [%s] %s\n", timestamp, level, msg))
  }
}

log_section <- function(title) {
  if (VERBOSE) {
    cat("\n")
    cat(strrep("=", 70), "\n")
    cat("  ", title, "\n")
    cat(strrep("=", 70), "\n\n")
  }
}

run_with_timing <- function(expr, name) {
  start_time <- Sys.time()
  tryCatch({
    result <- expr
    end_time <- Sys.time()
    runtime <- as.numeric(difftime(end_time, start_time, units = "mins"))
    log_step(sprintf("%s completed in %.2f minutes", name, runtime))
    return(list(success = TRUE, result = result, runtime = runtime))
  }, error = function(e) {
    log_step(sprintf("%s FAILED: %s", name, e$message), "ERROR")
    return(list(success = FALSE, error = e$message, runtime = NA))
  })
}

check_file_exists <- function(file_path, description = NULL) {
  if (file.exists(file_path)) {
    log_step(sprintf("Found: %s", file_path))
    return(TRUE)
  } else {
    desc <- if (!is.null(description)) paste0(" (", description, ")") else ""
    log_step(sprintf("Missing: %s%s", file_path, desc), "WARN")
    return(FALSE)
  }
}

# -----------------------------------------------------------------------------
# Main Pipeline
# -----------------------------------------------------------------------------

run_pipeline <- function() {
  pipeline_start <- Sys.time()

  log_section("scRNA-seq Clustering Comparative Study - Pipeline")
  log_step(sprintf("Working directory: %s", getwd()))
  log_step(sprintf("Datasets: %s", paste(DATASETS, collapse = ", ")))
  log_step(sprintf("Methods: %s", paste(METHODS, collapse = ", ")))

  if (DRY_RUN) {
    log_step("DRY RUN MODE - No actual execution", "WARN")
  }

  # Track results
  results <- list(
    download = list(),
    preprocessing = list(),
    clustering = list(),
    evaluation = NULL
  )

  # -------------------------------------------------------------------------
  # Step 0: Environment Setup
  # -------------------------------------------------------------------------
  log_section("Step 0: Environment Setup")

  log_step("Loading setup script...")
  source("src/R/00_setup.R")

  if (!DRY_RUN) {
    setup_result <- run_with_timing({
      setup_environment(install_missing = TRUE)
      set_seed(42)
    }, "Environment setup")

    if (!setup_result$success) {
      stop("Environment setup failed. Cannot continue.")
    }
  }

  # -------------------------------------------------------------------------
  # Step 1: Data Download
  # -------------------------------------------------------------------------
  if (!SKIP_DOWNLOAD) {
    log_section("Step 1: Data Download")

    log_step("Loading download script...")
    source("src/R/01_download_data.R")

    if (!DRY_RUN) {
      download_result <- run_with_timing({
        download_all_datasets("data/raw")
      }, "Data download")

      results$download <- download_result
    } else {
      log_step("Would download datasets to data/raw/")
    }
  } else {
    log_step("Skipping data download (--skip-download)")

    # Verify data exists
    for (dataset in DATASETS) {
      if (dataset == "scmixology") {
        check_file_exists("data/raw/scMixology.RData")
      } else if (dataset == "pbmc3k") {
        check_file_exists("data/raw/pbmc3k.rds")
      } else if (dataset == "zeisel") {
        check_file_exists("data/raw/zeisel_brain.rds")
      }
    }
  }

  # -------------------------------------------------------------------------
  # Step 2: Preprocessing
  # -------------------------------------------------------------------------
  if (!SKIP_PREPROCESSING) {
    log_section("Step 2: Preprocessing")

    log_step("Loading preprocessing script...")
    source("src/R/02_preprocessing.R")

    if (!DRY_RUN) {
      preprocess_result <- run_with_timing({
        preprocess_all_datasets(
          input_dir = "data/raw",
          output_dir = "data/processed",
          config_path = "config/preprocessing_params.yaml"
        )
      }, "Preprocessing")

      results$preprocessing <- preprocess_result
    } else {
      log_step("Would preprocess datasets to data/processed/")
    }
  } else {
    log_step("Skipping preprocessing (--skip-preprocessing)")

    # Verify preprocessed data exists
    for (dataset in DATASETS) {
      check_file_exists(
        file.path("data/processed", paste0(dataset, "_processed.rds"))
      )
    }
  }

  # -------------------------------------------------------------------------
  # Step 3: Clustering
  # -------------------------------------------------------------------------
  if (!SKIP_CLUSTERING) {
    log_section("Step 3: Clustering")

    # Seurat Leiden clustering
    if ("seurat" %in% METHODS) {
      log_step("Running Seurat Leiden clustering...")
      source("src/R/05_clustering_seurat.R")

      if (!DRY_RUN) {
        seurat_result <- run_with_timing({
          run_seurat_all_datasets(
            input_dir = "data/processed",
            output_dir = "data/results",
            config_path = "config/clustering_params.yaml"
          )
        }, "Seurat clustering")

        results$clustering$seurat <- seurat_result
      }
    }

    # CHOIR clustering
    if ("choir" %in% METHODS) {
      log_step("Running CHOIR clustering...")
      source("src/R/03_clustering_choir.R")

      if (!DRY_RUN) {
        choir_result <- run_with_timing({
          run_choir_all_datasets(
            input_dir = "data/processed",
            output_dir = "data/results",
            config_path = "config/clustering_params.yaml"
          )
        }, "CHOIR clustering")

        results$clustering$choir <- choir_result
      }
    }

    # recall clustering
    if ("recall" %in% METHODS) {
      log_step("Running recall clustering...")
      source("src/R/04_clustering_recall.R")

      if (!DRY_RUN) {
        recall_result <- run_with_timing({
          run_recall_all_datasets(
            input_dir = "data/processed",
            output_dir = "data/results",
            config_path = "config/clustering_params.yaml"
          )
        }, "recall clustering")

        results$clustering$recall <- recall_result
      }
    }
  } else {
    log_step("Skipping clustering (--skip-clustering)")
  }

  # -------------------------------------------------------------------------
  # Step 4: Evaluation
  # -------------------------------------------------------------------------
  if (!SKIP_EVALUATION) {
    log_section("Step 4: Evaluation")

    log_step("Loading evaluation script...")
    source("src/R/06_evaluation.R")

    if (!DRY_RUN) {
      eval_result <- run_with_timing({
        evaluate_all(
          results_dir = "data/results",
          output_file = "data/results/metrics_full.csv"
        )
      }, "Evaluation")

      results$evaluation <- eval_result

      # Print summary if successful
      if (eval_result$success && nrow(eval_result$result) > 0) {
        log_step("Evaluation Summary:")
        print(generate_summary_table(eval_result$result))
      }
    }
  } else {
    log_step("Skipping evaluation (--skip-evaluation)")
  }

  # -------------------------------------------------------------------------
  # Step 5: Visualization
  # -------------------------------------------------------------------------
  if (!SKIP_VISUALIZATION) {
    log_section("Step 5: Visualization")

    log_step("Loading visualization script...")
    source("src/R/07_visualization.R")

    if (!DRY_RUN) {
      viz_result <- run_with_timing({
        generate_all_figures(
          results_dir = "data/results",
          output_dir = "docs/latex/figures",
          metrics_file = "data/results/metrics_full.csv"
        )
        generate_umap_grids(
          results_dir = "data/results",
          output_dir = "docs/latex/figures"
        )
      }, "Visualization")

      results$visualization <- viz_result
    }
  } else {
    log_step("Skipping visualization (--skip-visualization)")
  }

  # -------------------------------------------------------------------------
  # Summary
  # -------------------------------------------------------------------------
  log_section("Pipeline Complete")

  pipeline_end <- Sys.time()
  total_runtime <- as.numeric(difftime(pipeline_end, pipeline_start, units = "mins"))
  log_step(sprintf("Total runtime: %.2f minutes", total_runtime))

  # Summary of results
  log_step("Results summary:")

  if (!is.null(results$evaluation) && results$evaluation$success) {
    metrics_file <- "data/results/metrics_full.csv"
    if (file.exists(metrics_file)) {
      log_step(sprintf("Metrics saved to: %s", metrics_file))
    }
  }

  figures_dir <- "docs/latex/figures"
  if (dir.exists(figures_dir)) {
    n_figures <- length(list.files(figures_dir, pattern = "\\.(pdf|png)$"))
    log_step(sprintf("Generated %d figures in: %s", n_figures, figures_dir))
  }

  # Return results for programmatic use
  invisible(results)
}

# -----------------------------------------------------------------------------
# Help Message
# -----------------------------------------------------------------------------

print_help <- function() {
  cat("
scRNA-seq Clustering Pipeline

Usage:
  Rscript src/R/run_pipeline.R [options]

Options:
  --skip-download       Skip data download step
  --skip-preprocessing  Skip preprocessing step
  --skip-clustering     Skip clustering step
  --skip-evaluation     Skip evaluation step
  --skip-visualization  Skip visualization step
  --methods=METHOD1,... Methods to run (default: seurat,choir,recall)
  --datasets=DATASET1,... Datasets to process (default: all)
  --dry-run             Show what would be done without executing
  --quiet               Suppress progress messages
  --help                Show this help message

Examples:
  # Run full pipeline
  Rscript src/R/run_pipeline.R

  # Run only Seurat clustering on PBMC dataset
  Rscript src/R/run_pipeline.R --methods=seurat --datasets=pbmc3k

  # Skip download and preprocessing (data already prepared)
  Rscript src/R/run_pipeline.R --skip-download --skip-preprocessing

  # Dry run to see what would happen
  Rscript src/R/run_pipeline.R --dry-run
")
}

# -----------------------------------------------------------------------------
# Main Execution
# -----------------------------------------------------------------------------

if ("--help" %in% args) {
  print_help()
} else if (sys.nframe() == 0) {
  # Script is being run directly
  run_pipeline()
} else {
  # Script is being sourced - make functions available
  message("Pipeline functions loaded. Run run_pipeline() to execute.")
}
