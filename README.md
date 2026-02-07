# scRNA-seq Clustering Comparative Study

A comparative study of clustering methods for single-cell RNA sequencing data, evaluating statistical calibration approaches against standard community detection algorithms.

## Overview

This project benchmarks clustering methods for scRNA-seq data to assess their ability to recover known cell types. We compare:

- **Seurat/Leiden**: Standard graph-based community detection at multiple resolutions (baseline)
- **SC3**: Single-Cell Consensus Clustering using ensemble k-means
- **recall**: Knockoff-based FDR-controlled clustering

### Key Findings

| Dataset | Best Method | Clusters | ARI | NMI |
|---------|-------------|----------|-----|-----|
| scMixology (3 cell lines) | Leiden (res=0.2) | 4 | **0.856** | 0.812 |
| PBMC 3k (9 cell types) | Leiden (res=1.0) | 8 | **0.781** | 0.837 |
| Zeisel Brain (7 classes) | Leiden (res=0.2) | 4 | **0.797** | 0.692 |

**Observations:**
- Leiden clustering at low resolution (0.2-0.5) generally achieves best ARI scores
- Higher resolutions lead to overclustering, reducing agreement with ground truth
- recall provides FDR-controlled clustering with competitive performance (ARI 0.48-0.68)
- SC3 achieves moderate performance (ARI=0.74 on scMixology) but has high memory requirements

## Quick Start (Docker - Recommended)

Docker provides a fully reproducible environment with all R packages pre-installed.

### Prerequisites

- Docker and Docker Compose installed

### Run the Complete Pipeline

```bash
# Clone the repository
git clone https://github.com/adrianstanea/scRNA-seq-Clustering.git
cd scRNA-seq-Clustering

# Build the Docker image
docker-compose build

# Run the full pipeline
docker-compose run scrna-clustering ./scripts/run_pipeline.sh
```

### Run Individual Steps

```bash
# Start an interactive container
docker-compose run scrna-clustering bash

# Inside the container:
Rscript src/R/01_download_data.R      # Download datasets
Rscript src/R/02_preprocessing.R      # Preprocess data
Rscript src/R/05_clustering_seurat.R  # Seurat/Leiden clustering
Rscript src/R/03_clustering_sc3.R     # SC3 clustering
Rscript src/R/04_clustering_recall.R  # recall clustering
Rscript src/R/06_evaluation.R         # Compute metrics
Rscript src/R/07_visualization.R      # Generate figures
```

## Project Structure

```
scRNA-seq-Clustering/
├── config/
│   ├── preprocessing_params.yaml    # QC and normalization parameters
│   └── clustering_params.yaml       # Clustering method parameters
├── data/
│   ├── raw/                         # Downloaded raw datasets
│   ├── processed/                   # Preprocessed Seurat objects
│   └── results/                     # Clustering results and metrics
├── docs/
│   └── latex/
│       └── figures/                 # Publication-ready figures
├── scripts/
│   └── run_pipeline.sh              # Master pipeline script
├── src/
│   └── R/
│       ├── 00_setup.R               # Environment setup
│       ├── 01_download_data.R       # Dataset acquisition
│       ├── 02_preprocessing.R       # Standardized preprocessing
│       ├── 03_clustering_sc3.R      # SC3 consensus clustering
│       ├── 04_clustering_recall.R   # recall FDR-controlled clustering
│       ├── 05_clustering_seurat.R   # Seurat/Leiden baseline
│       ├── 06_evaluation.R          # Metrics computation
│       └── 07_visualization.R       # Figure generation
├── Dockerfile                       # R 4.4 + Bioconductor 3.19
├── docker-compose.yml               # Container configuration
└── README.md
```

## Datasets

| Dataset | Cells | Features | Ground Truth | Source |
|---------|-------|----------|--------------|--------|
| scMixology | 902 | 16,468 | 3 cell lines | [Tian et al. 2019](https://github.com/LuyiTian/sc_mixology) |
| PBMC 3k | 2,643 | 13,697 | 9 cell types | [10x Genomics](https://support.10xgenomics.com/) |
| Zeisel Brain | 680 | 15,500 | 7 level1 classes | [Zeisel et al. 2015](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60361) |

Note: Cell/feature counts are after QC filtering.

## Methods

### Seurat/Leiden (Baseline)

Standard graph-based clustering using the Leiden algorithm at multiple resolutions.

```yaml
# config/clustering_params.yaml
leiden:
  resolutions: [0.2, 0.5, 0.8, 1.0, 1.5]
  algorithm: "leiden"
```

### SC3 (Consensus Clustering)

Single-Cell Consensus Clustering combines multiple k-means runs with different distance metrics and transformations.

```r
# Install
BiocManager::install("SC3")
```

**Parameters:**
- Automatic k estimation via Tracy-Widom theory
- Gene filtering enabled
- Biology computation for marker genes

### recall (FDR-Controlled)

Knockoff-based calibration method that controls false discovery rate during cluster identification.

```r
# Install
devtools::install_github("lcrawlab/recall")
```

**Parameters:**
- Algorithm: Louvain
- Null method: ZIP (Zero-Inflated Poisson)
- Iterative resolution adjustment

## Evaluation Metrics

- **ARI** (Adjusted Rand Index): Measures agreement with ground truth, adjusted for chance
- **NMI** (Normalized Mutual Information): Information-theoretic measure of clustering quality
- **Silhouette**: Internal validation based on cluster separation in PCA space

Results are saved to `data/results/metrics_full.csv`.

## Output Files

### Results
- `data/results/{dataset}_seurat.rds` - Seurat object with Leiden clusters
- `data/results/{dataset}_sc3.rds` - Seurat object with SC3 clusters
- `data/results/{dataset}_recall.rds` - Seurat object with recall clusters
- `data/results/metrics_full.csv` - All evaluation metrics

### Figures
- `docs/latex/figures/ari_comparison.pdf` - ARI barplot by method
- `docs/latex/figures/nmi_comparison.pdf` - NMI barplot by method
- `docs/latex/figures/cluster_counts.pdf` - Cluster count comparison
- `docs/latex/figures/umap_{dataset}.pdf` - UMAP visualizations

## Configuration

### Preprocessing Parameters (`config/preprocessing_params.yaml`)

```yaml
qc:
  min_genes_per_cell: 200
  max_genes_per_cell: 5000
  max_pct_mt: 5
  min_cells_per_gene: 3

normalization:
  method: "LogNormalize"
  scale_factor: 10000

feature_selection:
  method: "vst"
  n_top_genes: 2000

dimensionality_reduction:
  n_pcs: 30
```

### Clustering Parameters (`config/clustering_params.yaml`)

```yaml
leiden:
  resolutions: [0.2, 0.5, 0.8, 1.0, 1.5]

sc3:
  k_range: null  # Auto-estimate
  n_cores: 8
  biology: true

recall:
  dims: 30
  algorithm: "louvain"
  null_method: "ZIP"
```

## Known Limitations

### SC3 Memory Requirements
SC3 requires dense matrices for its internal computations:

- Successfully completes on smaller datasets (scMixology: 902 cells)
- Runs out of memory on larger datasets (PBMC3k: 2,643 cells; Zeisel: 680 cells with 15k genes)
