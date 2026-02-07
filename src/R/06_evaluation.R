# =============================================================================
# 06_evaluation.R - Metrics Computation
# Evaluation of clustering performance against ground truth
# =============================================================================

source("src/R/00_setup.R")

# -----------------------------------------------------------------------------
# Load Required Packages
# -----------------------------------------------------------------------------

load_eval_packages <- function() {
if (!requireNamespace("mclust", quietly = TRUE)) {
  install.packages("mclust")
}
if (!requireNamespace("cluster", quietly = TRUE)) {
  install.packages("cluster")
}
if (!requireNamespace("aricode", quietly = TRUE)) {
  install.packages("aricode")
}
library(mclust)
library(cluster)
library(aricode)
}

# -----------------------------------------------------------------------------
# Metric Functions
# -----------------------------------------------------------------------------

# Adjusted Rand Index
compute_ari <- function(true_labels, pred_labels) {
mclust::adjustedRandIndex(true_labels, pred_labels)
}

# Normalized Mutual Information
compute_nmi <- function(true_labels, pred_labels) {
aricode::NMI(true_labels, pred_labels)
}

# Homogeneity (all members of a cluster belong to the same class)
compute_homogeneity <- function(true_labels, pred_labels) {
# Entropy of classes given clusters
n <- length(true_labels)
clusters <- unique(pred_labels)

h_c_given_k <- 0
for (k in clusters) {
  mask <- pred_labels == k
  n_k <- sum(mask)
  if (n_k == 0) next

  class_counts <- table(true_labels[mask])
  probs <- class_counts / n_k
  probs <- probs[probs > 0]
  h_c_given_k <- h_c_given_k - (n_k / n) * sum(probs * log(probs))
}

# Entropy of classes
class_counts <- table(true_labels)
probs <- class_counts / n
probs <- probs[probs > 0]
h_c <- -sum(probs * log(probs))

if (h_c == 0) return(1.0)
1 - h_c_given_k / h_c
}

# Completeness (all members of a class are assigned to the same cluster)
compute_completeness <- function(true_labels, pred_labels) {
compute_homogeneity(pred_labels, true_labels)
}

# Silhouette Score
compute_silhouette <- function(sobj, cluster_col, reduction = "pca", dims = 30) {
embeddings <- Embeddings(sobj, reduction)[, 1:dims]
clusters <- as.numeric(factor(sobj@meta.data[[cluster_col]]))

if (length(unique(clusters)) < 2) {
  return(NA)
}

sil <- cluster::silhouette(clusters, dist(embeddings))
mean(sil[, "sil_width"])
}

# Over-clustering Index
compute_overclustering_index <- function(k_pred, k_true) {
if (k_true == 0) {
  return(NA_real_)
}
(k_pred - k_true) / k_true
}

# -----------------------------------------------------------------------------
# Evaluate Single Method
# -----------------------------------------------------------------------------

evaluate_clustering <- function(
sobj,
cluster_col,
ground_truth_col,
reduction = "pca",
dims = 30
) {
load_eval_packages()

pred_labels <- sobj@meta.data[[cluster_col]]
true_labels <- sobj@meta.data[[ground_truth_col]]

# Remove NAs
valid_idx <- !is.na(pred_labels) & !is.na(true_labels)
pred_labels <- pred_labels[valid_idx]
true_labels <- true_labels[valid_idx]

k_pred <- length(unique(pred_labels))
k_true <- length(unique(true_labels))

metrics <- list(
  method = cluster_col,
  k_predicted = k_pred,
  k_true = k_true,
  ARI = compute_ari(true_labels, pred_labels),
  NMI = compute_nmi(true_labels, pred_labels),
  Homogeneity = compute_homogeneity(true_labels, pred_labels),
  Completeness = compute_completeness(true_labels, pred_labels),
  Silhouette = compute_silhouette(sobj, cluster_col, reduction, dims),
  Overclustering_Index = compute_overclustering_index(k_pred, k_true)
)

as.data.frame(metrics)
}

# -----------------------------------------------------------------------------
# Evaluate All Methods for a Dataset
# -----------------------------------------------------------------------------

evaluate_dataset <- function(
dataset_name,
results_dir = "data/results",
ground_truth_col = "ground_truth"
) {
message("Evaluating: ", dataset_name)

results <- data.frame()

# Load and evaluate each method
methods <- list(
  seurat = list(
    file = paste0(dataset_name, "_seurat.rds"),
    cols = c("leiden_res0_2", "leiden_res0_5", "leiden_res0_8", "leiden_res1", "leiden_res1_5")
  ),
  sc3 = list(
    file = paste0(dataset_name, "_sc3.rds"),
    cols = c("SC3_clusters")
  ),
  recall = list(
    file = paste0(dataset_name, "_recall.rds"),
    cols = c("recall_clusters")
  )
)

for (method_name in names(methods)) {
  method_info <- methods[[method_name]]
  file_path <- file.path(results_dir, method_info$file)

  if (!file.exists(file_path)) {
    message("  Skipping ", method_name, " (file not found)")
    next
  }

  sobj <- readRDS(file_path)

  # Check ground truth exists
  if (!ground_truth_col %in% colnames(sobj@meta.data)) {
    message("  Ground truth column '", ground_truth_col, "' not found")
    next
  }

  for (col in method_info$cols) {
    if (!col %in% colnames(sobj@meta.data)) {
      message("  Skipping ", col, " (column not found)")
      next
    }

    tryCatch({
      metrics <- evaluate_clustering(sobj, col, ground_truth_col)
      metrics$dataset <- dataset_name
      results <- rbind(results, metrics)
    }, error = function(e) {
      message("  Error evaluating ", col, ": ", e$message)
    })
  }
}

results
}

# -----------------------------------------------------------------------------
# Evaluate All Datasets
# -----------------------------------------------------------------------------

evaluate_all <- function(
results_dir = "data/results",
output_file = "data/results/metrics_full.csv"
) {
# Ground truth column names for each dataset
# These should match what was set in preprocessing
gt_cols <- list(
  scmixology = "ground_truth",
  pbmc3k = "ground_truth",
  zeisel = "ground_truth"
)

# Alternative ground truth column names to try if ground_truth is missing
alt_gt_cols <- list(
  scmixology = c("cell_line", "cell_line_demuxlet"),
  pbmc3k = c("seurat_annotations", "celltype", "louvain"),
  zeisel = c("level1class", "level2class", "cell_type1")
)

all_results <- data.frame()

for (dataset in names(gt_cols)) {
  tryCatch({
    # Try primary ground truth column first
    gt_col <- gt_cols[[dataset]]

    # Check if we need to use alternative column
    # Load one result file to check available columns
    test_file <- file.path(results_dir, paste0(dataset, "_seurat.rds"))
    if (file.exists(test_file)) {
      sobj <- readRDS(test_file)
      if (!gt_col %in% colnames(sobj@meta.data)) {
        # Try alternatives
        for (alt_col in alt_gt_cols[[dataset]]) {
          if (alt_col %in% colnames(sobj@meta.data)) {
            gt_col <- alt_col
            message("Using alternative ground truth column: ", alt_col, " for ", dataset)
            break
          }
        }
      }
    }

    results <- evaluate_dataset(dataset, results_dir, gt_col)
    all_results <- rbind(all_results, results)
  }, error = function(e) {
    message("Failed to evaluate ", dataset, ": ", e$message)
  })
}

if (nrow(all_results) > 0) {
  write.csv(all_results, output_file, row.names = FALSE)
  message("Results saved to: ", output_file)
}

all_results
}

# -----------------------------------------------------------------------------
# Summary Table Generation
# -----------------------------------------------------------------------------

generate_summary_table <- function(metrics_df) {
summary_df <- metrics_df %>%
  group_by(dataset, method) %>%
  summarize(
    K = first(k_predicted),
    ARI = round(mean(ARI, na.rm = TRUE), 3),
    NMI = round(mean(NMI, na.rm = TRUE), 3),
    Silhouette = round(mean(Silhouette, na.rm = TRUE), 3),
    .groups = "drop"
  )

summary_df
}

# -----------------------------------------------------------------------------
# Export to LaTeX
# -----------------------------------------------------------------------------

export_to_latex <- function(df, output_file = "docs/latex/tables/metrics_full.tex") {
if (!requireNamespace("xtable", quietly = TRUE)) {
  install.packages("xtable")
}
library(xtable)

latex_table <- xtable(
  df,
  caption = "Clustering performance metrics across datasets and methods",
  label = "tab:metrics"
)

print(
  latex_table,
  file = output_file,
  include.rownames = FALSE,
  booktabs = TRUE
)

message("LaTeX table saved to: ", output_file)
}

# -----------------------------------------------------------------------------
# Main execution
# -----------------------------------------------------------------------------

if (sys.nframe() == 0) {
setup_environment()
results <- evaluate_all()
if (nrow(results) > 0) {
  print(generate_summary_table(results))
}
}
