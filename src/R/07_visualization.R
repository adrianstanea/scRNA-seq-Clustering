# =============================================================================
# 07_visualization.R - Figures and Plots
# Publication-quality visualizations for comparative analysis
# =============================================================================

source("src/R/00_setup.R")

# -----------------------------------------------------------------------------
# Theme and Color Settings
# -----------------------------------------------------------------------------

# Publication-ready theme
theme_publication <- function(base_size = 12) {
theme_minimal(base_size = base_size) +
  theme(
    plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = base_size, hjust = 0.5),
    axis.title = element_text(size = base_size),
    axis.text = element_text(size = base_size - 1),
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = base_size - 1),
    strip.text = element_text(size = base_size, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA)
  )
}

# Colorblind-friendly palette
get_cluster_colors <- function(n) {
if (n <= 8) {
  RColorBrewer::brewer.pal(max(3, n), "Set2")[1:n]
} else if (n <= 12) {
  RColorBrewer::brewer.pal(n, "Set3")
} else {
  viridis::viridis(n)
}
}

# -----------------------------------------------------------------------------
# UMAP Plotting
# -----------------------------------------------------------------------------

plot_umap <- function(
sobj,
group_by,
title = NULL,
colors = NULL,
pt_size = 0.5,
label = TRUE
) {
p <- DimPlot(
  sobj,
  reduction = "umap",
  group.by = group_by,
  pt.size = pt_size,
  label = label,
  label.size = 3
) +
  theme_publication() +
  theme(legend.position = "right")

if (!is.null(title)) {
  p <- p + ggtitle(title)
}

if (!is.null(colors)) {
  p <- p + scale_color_manual(values = colors)
}

p
}

# -----------------------------------------------------------------------------
# UMAP Comparison Grid
# -----------------------------------------------------------------------------

plot_umap_comparison <- function(
sobj,
cluster_cols,
titles = NULL,
ncol = 3
) {
if (is.null(titles)) {
  titles <- cluster_cols
}

plots <- lapply(seq_along(cluster_cols), function(i) {
  plot_umap(sobj, cluster_cols[i], titles[i])
})

wrap_plots(plots, ncol = ncol) +
  plot_layout(guides = "collect")
}

# -----------------------------------------------------------------------------
# Metrics Bar Plot
# -----------------------------------------------------------------------------

plot_metrics_barplot <- function(
metrics_df,
metric = "ARI",
title = NULL
) {
if (is.null(title)) {
  title <- paste(metric, "by Method and Dataset")
}

ggplot(metrics_df, aes(x = method, y = .data[[metric]], fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(label = round(.data[[metric]], 2)),
    position = position_dodge(width = 0.8),
    vjust = -0.5,
    size = 3
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = title,
    x = "Method",
    y = metric,
    fill = "Dataset"
  ) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# -----------------------------------------------------------------------------
# Cluster Count Comparison
# -----------------------------------------------------------------------------

plot_cluster_counts <- function(
metrics_df,
true_k_col = "k_true",
pred_k_col = "k_predicted"
) {
# Reshape data
df_long <- metrics_df %>%
  select(dataset, method, k_true = all_of(true_k_col), k_predicted = all_of(pred_k_col)) %>%
  pivot_longer(
    cols = c(k_true, k_predicted),
    names_to = "type",
    values_to = "k"
  ) %>%
  mutate(type = ifelse(type == "k_true", "Ground Truth", "Predicted"))

ggplot(df_long, aes(x = method, y = k, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~dataset, scales = "free_y") +
  scale_fill_manual(values = c("Ground Truth" = "#4DAF4A", "Predicted" = "#377EB8")) +
  labs(
    title = "Cluster Count: Predicted vs Ground Truth",
    x = "Method",
    y = "Number of Clusters",
    fill = ""
  ) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# -----------------------------------------------------------------------------
# Overclustering Analysis
# -----------------------------------------------------------------------------

plot_overclustering <- function(metrics_df) {
ggplot(metrics_df, aes(x = method, y = Overclustering_Index, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Overclustering Index by Method",
    subtitle = "Positive = overclustering, Negative = underclustering",
    x = "Method",
    y = "Overclustering Index",
    fill = "Dataset"
  ) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# -----------------------------------------------------------------------------
# Runtime Comparison
# -----------------------------------------------------------------------------

plot_runtime <- function(runtime_df) {
ggplot(runtime_df, aes(x = method, y = runtime, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Runtime Comparison",
    x = "Method",
    y = "Runtime (seconds)",
    fill = "Dataset"
  ) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# -----------------------------------------------------------------------------
# Generate All Figures
# -----------------------------------------------------------------------------

generate_all_figures <- function(
results_dir = "data/results",
output_dir = "docs/latex/figures",
metrics_file = "data/results/metrics_full.csv"
) {
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load metrics
if (!file.exists(metrics_file)) {
  warning("Metrics file not found: ", metrics_file)
  return(NULL)
}

metrics_df <- read.csv(metrics_file)

# 1. ARI Bar Plot
message("Generating ARI barplot...")
p_ari <- plot_metrics_barplot(metrics_df, "ARI", "Adjusted Rand Index by Method")
ggsave(
  file.path(output_dir, "ari_comparison.pdf"),
  p_ari, width = 8, height = 5
)

# 2. NMI Bar Plot
message("Generating NMI barplot...")
p_nmi <- plot_metrics_barplot(metrics_df, "NMI", "Normalized Mutual Information by Method")
ggsave(
  file.path(output_dir, "nmi_comparison.pdf"),
  p_nmi, width = 8, height = 5
)

# 3. Cluster Count Comparison
message("Generating cluster count comparison...")
p_k <- plot_cluster_counts(metrics_df)
ggsave(
  file.path(output_dir, "cluster_counts.pdf"),
  p_k, width = 10, height = 6
)

# 4. Overclustering Analysis
message("Generating overclustering analysis...")
p_oc <- plot_overclustering(metrics_df)
ggsave(
  file.path(output_dir, "overclustering_analysis.pdf"),
  p_oc, width = 8, height = 5
)

message("All figures saved to: ", output_dir)
}

# -----------------------------------------------------------------------------
# UMAP Grid for All Datasets
# -----------------------------------------------------------------------------

generate_umap_grids <- function(
results_dir = "data/results",
output_dir = "docs/latex/figures"
) {
datasets <- c("scmixology", "pbmc3k", "zeisel")

for (dataset in datasets) {
  message("Generating UMAP grid for: ", dataset)

  # Try to load results
  files <- list(
    seurat = file.path(results_dir, paste0(dataset, "_seurat.rds")),
    choir = file.path(results_dir, paste0(dataset, "_choir.rds")),
    recall = file.path(results_dir, paste0(dataset, "_recall.rds"))
  )

  plots <- list()

  # Seurat baseline (res 0.8)
  if (file.exists(files$seurat)) {
    sobj <- readRDS(files$seurat)
    if ("leiden_res0_8" %in% colnames(sobj@meta.data)) {
      plots$seurat <- plot_umap(sobj, "leiden_res0_8", "Leiden (res=0.8)")
    }
  }

  # CHOIR
  if (file.exists(files$choir)) {
    sobj <- readRDS(files$choir)
    if ("CHOIR_clusters" %in% colnames(sobj@meta.data)) {
      plots$choir <- plot_umap(sobj, "CHOIR_clusters", "CHOIR")
    }
  }

  # recall
  if (file.exists(files$recall)) {
    sobj <- readRDS(files$recall)
    if ("recall_clusters" %in% colnames(sobj@meta.data)) {
      plots$recall <- plot_umap(sobj, "recall_clusters", "recall")
    }
  }

  if (length(plots) > 0) {
    combined <- wrap_plots(plots, ncol = length(plots)) +
      plot_annotation(title = paste("Clustering Comparison:", dataset))

    ggsave(
      file.path(output_dir, paste0("umap_", dataset, ".pdf")),
      combined, width = 4 * length(plots), height = 4
    )
  }
}
}

# -----------------------------------------------------------------------------
# Main execution
# -----------------------------------------------------------------------------

if (sys.nframe() == 0) {
setup_environment()
generate_all_figures()
generate_umap_grids()
}
