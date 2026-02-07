# =============================================================================
# Dockerfile for scRNA-seq Clustering Comparative Study
# R 4.4 + Python 3.12 + Seurat 5 + SC3 + recall
# =============================================================================

FROM rocker/r-ver:4.4.0

LABEL maintainer="adrianstanea1@gmail.com"
LABEL description="scRNA-seq Clustering Comparative Study Environment"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Handle CA certificates (optional - build succeeds without certs/)
COPY cert[s]/ /tmp/certs/
RUN if [ -f /tmp/certs/nscacert_combined.crt ]; then \
    cp /tmp/certs/nscacert_combined.crt /usr/local/share/ca-certificates/nscacert_combined.crt && \
    update-ca-certificates; \
    fi

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Build tools
    build-essential \
    gfortran \
    cmake \
    # Python
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    # R package dependencies
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libgit2-dev \
    libgmp-dev \
    libmpfr-dev \
    libglpk-dev \
    libhdf5-dev \
    libmagick++-dev \
    libgsl-dev \
    libcairo2-dev \
    libxt-dev \
    # Other utilities
    git \
    wget \
    curl \
    vim \
    ca-certificates \
    && update-ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install uv for Python package management
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:$PATH"

# Set up R library path and configure repos
ENV R_LIBS_USER=/usr/local/lib/R/site-library
RUN mkdir -p $R_LIBS_USER

# Configure R to use Posit Package Manager (fast binary packages for Ubuntu)
RUN echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))' >> /usr/local/lib/R/etc/Rprofile.site

# Install BiocManager first
RUN Rscript -e "install.packages('BiocManager')"

# Set Bioconductor version
RUN Rscript -e "BiocManager::install(version = '3.19', ask = FALSE, update = FALSE)"

# Install core R packages
RUN Rscript -e "install.packages(c('remotes', 'devtools', 'renv'))"

# Install Seurat and dependencies (using Posit binaries for speed)
RUN Rscript -e "install.packages(c( \
    'Seurat', \
    'SeuratObject', \
    'sctransform', \
    'leiden', \
    'igraph', \
    'harmony', \
    'clustree', \
    'ggplot2', \
    'dplyr', \
    'tidyr', \
    'patchwork', \
    'cowplot', \
    'ggrepel', \
    'pheatmap', \
    'viridis', \
    'RColorBrewer' \
    ))"

# Install Bioconductor packages
RUN Rscript -e "BiocManager::install(c( \
    'SingleCellExperiment', \
    'scran', \
    'scater', \
    'scuttle', \
    'bluster', \
    'BiocNeighbors', \
    'BiocSingular', \
    'scRNAseq', \
    'AnnotationHub', \
    'ensembldb' \
    ), ask = FALSE, update = FALSE)"

# Install SC3 (consensus clustering - replaces CHOIR)
RUN Rscript -e "BiocManager::install('SC3', ask = FALSE, update = FALSE)"

# Install recall dependencies
RUN Rscript -e "install.packages(c('knockoff', 'glmnet', 'Matrix'))"

# Install leidenbase for Seurat Leiden clustering
RUN Rscript -e "install.packages('leidenbase')"

# Install evaluation packages
RUN Rscript -e "install.packages(c('mclust', 'cluster', 'aricode', 'xtable', 'yaml'))"

# Install GitHub packages (presto for fast DE, SeuratData for datasets, recall for FDR clustering)
RUN Rscript -e "remotes::install_github('immunogenomics/presto', repos = BiocManager::repositories(), upgrade = 'never')"
RUN Rscript -e "remotes::install_github('satijalab/seurat-data', repos = BiocManager::repositories(), upgrade = 'never')"
RUN Rscript -e "remotes::install_github('lcrawlab/recall', repos = BiocManager::repositories(), upgrade = 'never')"

# Set up working directory
WORKDIR /workspace

# Copy project files
COPY . /workspace/

# Create Python virtual environment with uv
RUN cd /workspace && uv sync || true

# Set default command
CMD ["R"]
