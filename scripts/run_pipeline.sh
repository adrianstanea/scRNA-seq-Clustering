#!/bin/bash
# =============================================================================
# run_pipeline.sh - Execute full scRNA-seq clustering comparative study
# Run inside Docker container: docker-compose run scrna-clustering ./scripts/run_pipeline.sh
# =============================================================================

set -e  # Exit on error

echo "=============================================="
echo "scRNA-seq Clustering Comparative Study"
echo "=============================================="
echo ""

# Create directories
mkdir -p data/raw data/processed data/results docs/latex/figures

# Step 1: Download datasets
echo "[1/7] Downloading datasets..."
Rscript src/R/01_download_data.R
echo "Done."
echo ""

# Step 2: Preprocess datasets
echo "[2/7] Preprocessing datasets..."
Rscript src/R/02_preprocessing.R
echo "Done."
echo ""

# Step 3: Run Seurat/Leiden baseline clustering
echo "[3/7] Running Seurat/Leiden clustering..."
Rscript src/R/05_clustering_seurat.R
echo "Done."
echo ""

# Step 4: Run SC3 consensus clustering
echo "[4/7] Running SC3 clustering..."
Rscript src/R/03_clustering_sc3.R
echo "Done."
echo ""

# Step 5: Run recall clustering
echo "[5/7] Running recall clustering..."
Rscript src/R/04_clustering_recall.R
echo "Done."
echo ""

# Step 6: Compute evaluation metrics
echo "[6/7] Computing evaluation metrics..."
Rscript src/R/06_evaluation.R
echo "Done."
echo ""

# Step 7: Generate visualizations
echo "[7/7] Generating visualizations..."
Rscript src/R/07_visualization.R
echo "Done."
echo ""

echo "=============================================="
echo "Pipeline completed successfully!"
echo "Results saved to: data/results/"
echo "Figures saved to: docs/latex/figures/"
echo "=============================================="
