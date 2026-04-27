# Changelog

All notable changes to the Xenium Spatial Clustering and Annotation Tool will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed (breaking)

- Renamed the Python package from `xenium_process` to `spatial_tk`, the distribution/PyPI name from `xenium-process` to `spatial-tk`, and the CLI entry point from `xenium_process` to `spatial-tk`. Run `pip install spatial-tk` (or install from source) and invoke `spatial-tk ...` or `python -m spatial_tk.cli ...`.

## [1.3.0] - 2025-11-09

### Added
- **Memory-Efficient Table I/O**:
  - `load_table_only()`: Direct AnnData loading from `zarr_path/tables/table` without loading SpatialData
  - `save_table_only()`: Direct AnnData writing to zarr without loading SpatialData
  - Automatic table name detection (handles custom table names)
  - Fallback to SpatialData read if direct path fails (for robustness)
- **Selective Image Loading**:
  - `load_existing_spatial_data()`: Added `load_images=False` parameter to skip image loading
  - `load_spatial_datasets()`: Added `load_images=True` parameter (default True for concat)
  - Images automatically removed from memory when not needed
  - On-demand image loading in visualization script (`load_images_on_demand()`)
- **Raw Xenium Dataset Support**:
  - Automatic detection of `.zarr` vs raw Xenium dataset directories
  - `load_xenium_dataset()`: Loads raw Xenium datasets using `spatialdata_io.xenium()`
  - Automatic image loading from `morphology_focus` directory
  - `setup_squidpy_structure()`: Configures `adata.uns['spatial']` for squidpy compatibility

### Changed
- **Memory Optimization for Processing Commands**:
  - `normalize` and `cluster` commands now use `load_table_only()` for all operations
  - Inplace operations (`--inplace`) use `save_table_only()` - no SpatialData loading at all
  - Non-inplace operations load SpatialData without images for efficiency
  - Significant memory savings (90%+ for image-heavy datasets) during clustering/normalization
- **Visualization Script**:
  - Images only loaded when `--overlay-image` flag is used
  - On-demand image loading when overlay is requested
  - Automatic image path detection and loading from zarr
- **Data Loading Strategy**:
  - Concat command: Loads images by default (may be needed for downstream visualization)
  - Normalize/Cluster commands: Skip images entirely (not needed for processing)
  - Visualization: Load images only when explicitly requested

### Performance Improvements
- **Memory Usage**:
  - Table-only operations avoid loading images, shapes, and SpatialData metadata
  - Direct AnnData I/O eliminates SpatialData object overhead
  - Inplace operations are now most efficient (no SpatialData loading)
- **Loading Speed**:
  - Direct table I/O faster than loading full SpatialData and extracting table
  - Image loading deferred until actually needed
  - Reduced I/O overhead for processing workflows

## [1.2.0] - 2025-11-05

### Changed
- **Enrichment Scoring Method**: Switched from Univariate Linear Model (ULM) to Multivariate Linear Model (MLM)
  - All enrichment scoring now uses `dc.mt.mlm()` instead of `dc.mt.ulm()`
  - Function renamed: `calculate_ulm_scores()` → `calculate_mlm_scores()`
  - Obsm keys updated: `score_ulm_{resource}` → `score_mlm_{resource}`
  - MLM provides more accurate activity estimates by accounting for correlated regulators and their interactions
  - CLI flag `--calculate-ulm` retained for backward compatibility (still calculates MLM scores)
- **Plotting**: Updated enrichment heatmap function to correctly detect MLM score keys in obsm

### Documentation
- Added `notes/ULM_vs_MLM.md` with detailed explanation of differences and toy example
- Updated all documentation references from ULM to MLM
- Updated configuration file examples and help text

## [1.1.0] - 2025-11-04

### Added
- **TOML Configuration File Support**:
  - All commands now accept an optional `--config` argument for TOML configuration files
  - Each command reads its own section from the config file (`[concat]`, `[normalize]`, `[cluster]`, `[annotate]`, `[differential]`)
  - Config values can be overridden by CLI arguments for maximum flexibility
  - Automatic type conversion for config values (strings to ints/floats/bools)
  - Support for underscore/hyphen key mapping (e.g., `min_genes` in config maps to `--min-genes` CLI arg)
- **Makefile Run Target**:
  - New `make run ROOT=/path/to/directory` target for running full pipeline
  - Executes all five pipeline steps sequentially using a single config file
  - Validates config file existence and stops on any step failure
  - Changes working directory to ROOT for relative path resolution
- **Configuration Documentation**:
  - `example_config.toml` - Complete example with all command sections documented
  - `notes/TOML_CONFIG_GUIDE.md` - Comprehensive guide to TOML config files
  - `notes/MAKEFILE_RUN_TARGET.md` - Documentation for the make run target
- **Config Utility Module**:
  - `spatial_tk/utils/config.py` with `load_config()` and `merge_config_with_args()` functions
  - Uses Python 3.11+ built-in `tomllib` module (no additional dependencies)
  - Handles config loading, merging with CLI args, and type conversion

### Changed
- **Command Argument Handling**:
  - Required arguments (`--input`, `--output`, etc.) are now optional when `--config` is provided
  - Arguments are validated after config merge instead of during argparse parsing
  - Clear error messages when required values are missing from both CLI and config
- **CLI Behavior**:
  - Commands can now be run with only `--config` argument when all required parameters are in config file
  - CLI arguments always override config values (standard precedence)
  - Config files support relative paths resolved from config file location

### Testing
- **Unit Tests**:
  - `tests/unit/test_config.py` - Tests for config loading, merging, and type conversion
- **Functional Tests**:
  - `tests/functional/test_config_integration.py` - Integration tests for config files with all commands
  - Tests verify config values override defaults, CLI overrides config, and error handling

### Documentation
- Updated README.md with configuration file section
- Added examples of using config files in workflows
- Documented config key naming conventions (underscores vs hyphens)

## [1.0.0] - 2025-10-28

### Major Refactor - Xenium Spatial Data Support

This is a major rewrite transforming the tool from single-cell RNA-seq to Xenium spatial transcriptomics analysis.

### Added
- **Xenium Spatial Data Support**:
  - Load and concatenate multiple Xenium .zarr datasets
  - CSV input format with sample metadata (sample, path, optional columns)
  - SpatialData integration for spatial coordinate preservation
  - Output as .zarr format preserving spatial information
- **Modular Architecture**:
  - Refactored ~1000 line monolithic script into organized modules:
    - `data_io.py`: Data loading, concatenation, and saving
    - `preprocessing.py`: QC, filtering, normalization, HVG selection
    - `clustering.py`: PCA, neighbors, UMAP, Leiden clustering
    - `annotation.py`: Marker loading, cell type annotation, ULM scores, DE analysis
    - `plotting.py`: All visualization functions
    - `main.py`: CLI and workflow orchestration
- **ULM Enrichment Scoring**:
  - Pre-calculate ULM scores for pathway/TF resources via `--calculate-ulm` flag
  - Resources included: hallmark, collectri, dorothea, progeny, PanglaoDB
  - PanglaoDB markers filtered by canonical status and sensitivity (default: >0.5)
  - Scores stored in `adata.obsm['score_ulm_{resource}']`
  - Configurable PanglaoDB sensitivity via `--panglao-min-sensitivity`
- **Enhanced Dependencies**:
  - Added `spatialdata` for spatial data handling
  - Added `squidpy` for future spatial analysis features
- **Multi-sample Support**:
  - Automatic sample concatenation with metadata preservation
  - Sample-aware batch correction in HVG selection
  - UMAP visualization colored by sample

### Changed
- **Input Format**: Now accepts CSV file with sample paths instead of single h5ad file
- **Output Format**: Saves processed data as .zarr SpatialData object instead of h5ad
- **Main Script**: Renamed from `scanpy_cluster.py` to `main.py` (legacy script preserved)
- **Cell Type Annotation**: Updated to use latest decoupler API:
  - `dc.run_mlm()` instead of `dc.mt.mlm()`
  - `dc.get_acts()` instead of `dc.pp.get_obsm()`
  - `dc.rank_sources_groups()` instead of `dc.tl.rankby_group()`
- **Documentation**: Completely rewritten README for Xenium workflow

### Maintained
- All existing features from v0.2.0:
  - Quality control and filtering
  - Normalization and feature selection
  - PCA, UMAP, Leiden clustering
  - Differential expression analysis
  - Cell type annotation with markers
  - Resume functionality
  - Downsampling
  - Multiple clustering resolutions
  - Comprehensive plotting

### Notes
- Legacy `scanpy_cluster.py` remains available for h5ad single-cell workflows
- Spatial-specific analyses (spatial autocorrelation, niche detection, etc.) planned for future releases

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

