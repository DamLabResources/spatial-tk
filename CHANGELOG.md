# Changelog

All notable changes to the Scanpy Single-Cell RNA-seq Clustering and Annotation Tool will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2025-10-22

### Added
- **Differential Expression Analysis**: Automated identification of marker genes for each cluster
  - Uses `sc.tl.rank_genes_groups()` with Wilcoxon rank-sum test
  - Runs automatically for all clustering resolutions
  - Supports resume functionality to skip already computed analyses
- **Differential Expression Output Files**:
  - `deg_all_clusters_res{resolution}.csv`: Complete DE results for all genes across all clusters
  - `deg_top100_per_cluster_res{resolution}.csv`: Top 100 marker genes per cluster with statistics
  - Files organized in `differential_expression/` subdirectory
- **Differential Expression Visualizations**:
  - Dotplot showing top 5 DE genes per cluster (`deg_dotplot_res{resolution}.png`)
  - Heatmap showing top 10 DE genes with expression patterns (`deg_heatmap_res{resolution}.png`)
  - Both visualizations generated automatically when `--save-plots` flag is used
- **Enhanced Documentation**:
  - Added detailed documentation for differential expression features
  - Updated script docstring to reflect new capabilities

### Changed
- Updated workflow to integrate DE analysis after clustering and before cell type annotation
- Enhanced `save_results()` function to automatically save DE results for each resolution
- Expanded `save_plots()` function to generate DE visualizations

## [0.1.0] - 2025-10-21

### Added
- **Initial Release**: Complete scRNA-seq analysis pipeline
- **Data Loading and Preprocessing**:
  - Support for h5ad format input files
  - Automatic handling of variable and observation name uniqueness
  - Optional downsampling for quick testing
- **Quality Control**:
  - Calculation of QC metrics (mitochondrial, ribosomal, hemoglobin gene percentages)
  - Configurable cell and gene filtering thresholds
  - QC visualization plots (violin plots, scatter plots)
- **Normalization and Feature Selection**:
  - Median total count normalization
  - Log transformation
  - Highly variable gene selection (default: 2000 genes)
  - Batch-aware feature selection when sample information available
- **Dimensionality Reduction**:
  - PCA analysis with variance ratio plots
  - Neighborhood graph computation
  - UMAP embedding for visualization
- **Clustering**:
  - Leiden clustering with configurable resolution(s)
  - Support for multiple clustering resolutions in a single run
  - Results stored with unique keys for each resolution
- **Cell Type Annotation**:
  - Marker-based annotation using decoupler's MLM (multivariate linear model)
  - CSV format for marker gene input (cell_type, gene columns)
  - Automatic enrichment score calculation per cell type
  - Cluster-level annotation based on top enrichment scores
- **Visualizations** (with `--save-plots` flag):
  - QC plots: violin plots, scatter plots
  - Highly variable genes plot
  - PCA variance ratio plot
  - UMAP colored by sample (if available)
  - UMAP colored by clusters for each resolution
  - UMAP colored by cell type annotations for each resolution
  - Marker gene dotplots for each resolution
  - Enrichment score heatmaps for each resolution
- **Output Files**:
  - Processed AnnData object (`processed_data.h5ad`)
  - Cell type annotation CSVs for each resolution
  - Organized plot directory structure
- **Command-Line Interface**:
  - Required arguments: `--input`, `--output-dir`
  - Optional arguments: `--markers`, `--save-plots`, `--min-genes`, `--min-cells`, 
    `--n-top-genes`, `--leiden-resolution`, `--downsample`, `--resume`
  - Comprehensive help documentation and usage examples
- **Resume Functionality**:
  - Ability to resume from existing analysis
  - Skips already computed steps (QC, normalization, PCA, UMAP, clustering, annotation)
  - Useful for adding new markers or resolutions without recomputing everything
- **Logging**:
  - Detailed logging with timestamps
  - Progress tracking for all major steps
  - Error handling with informative messages
- **Performance Features**:
  - Non-interactive backend for plot generation
  - Efficient processing of large datasets
  - Graceful handling of missing marker genes

### Technical Details
- Based on Scanpy and the Scverse ecosystem
- Follows best practices from Scverse tutorials
- Uses igraph-based Leiden algorithm for clustering
- Implements decoupler for marker-based enrichment analysis
- Compatible with standard h5ad file formats

