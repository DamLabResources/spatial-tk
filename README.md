# Scanpy Clustering Tool

A comprehensive single-file Python tool for single-cell RNA-seq analysis. This script takes a single-cell sequencing dataset in h5 format and runs through a complete analysis workflow:

- **Quality Control** - Calculate QC metrics, filter cells and genes
- **Preprocessing** - Normalization, log transformation, feature selection
- **Dimensionality Reduction** - PCA and UMAP embedding
- **Clustering** - Leiden clustering at multiple resolutions
- **Cell Type Annotation** - Marker gene-based cell type assignment

The tool outputs a processed AnnData object (h5ad) and CSV files with cell type annotations for downstream analysis.

This tool is an operational version of the Scverse basic tutorial:
https://scverse-tutorials.readthedocs.io/en/latest/notebooks/basic-scrna-tutorial.html

## Installation

Install dependencies using pip:

```bash
pip install -r requirements.txt
```

## Usage

### Basic Usage

Process a dataset with default parameters:

```bash
python scanpy_cluster.py --input data.h5 --output-dir results/
```

### With Marker-Based Annotation

Annotate cell types using a marker gene CSV file:

```bash
python scanpy_cluster.py \
  --input data.h5 \
  --output-dir results/ \
  --markers example_markers.csv \
  --save-plots
```

### Multiple Clustering Resolutions

Run clustering at multiple resolutions (comma-separated):

```bash
python scanpy_cluster.py \
  --input data.h5 \
  --output-dir results/ \
  --leiden-resolution 0.2,0.5,1.0,2.0 \
  --markers example_markers.csv
```

### Quick Testing with Downsampling

Test the pipeline quickly on a subset of cells:

```bash
python scanpy_cluster.py \
  --input data.h5 \
  --output-dir results/ \
  --downsample 0.1 \
  --save-plots
```

### Resume from Existing Analysis

Resume from a previously saved analysis, skipping already-computed steps:

```bash
python scanpy_cluster.py \
  --input data.h5 \
  --output-dir results/ \
  --markers new_markers.csv \
  --resume \
  --save-plots
```

### All Options

```bash
python scanpy_cluster.py \
  --input data.h5 \
  --output-dir results/ \
  --markers example_markers.csv \
  --save-plots \
  --min-genes 100 \
  --min-cells 3 \
  --n-top-genes 2000 \
  --leiden-resolution 0.5 \
  --downsample 1.0
```

## Command-Line Arguments

### Required Arguments
- `--input`: Path to input h5 or h5ad file
- `--output-dir`: Directory where all output files will be saved

### Optional Arguments
- `--markers`: Path to CSV file with marker genes (format: `cell_type,gene`)
- `--save-plots`: Flag to generate and save QC/analysis plots
- `--min-genes`: Minimum genes per cell for filtering (default: 100)
- `--min-cells`: Minimum cells per gene for filtering (default: 3)
- `--n-top-genes`: Number of highly variable genes (default: 2000)
- `--leiden-resolution`: Clustering resolution(s), comma-separated (default: 0.5)
- `--downsample`: Fraction of cells to keep for analysis, 0-1 (default: 1.0 = no downsampling)
- `--resume`: Resume from existing analysis, skip already-computed steps

## Marker Gene CSV Format

The marker gene file should be a CSV with two columns: `cell_type` and `gene`. Each row specifies one marker gene for a cell type. See `example_markers.csv` for a template:

```csv
cell_type,gene
B cells,MS4A1
B cells,CD19
T cells,CD3D
T cells,CD3E
NK cells,GNLY
NK cells,NKG7
```

## Output Files

The tool creates the following outputs in the specified output directory:

### Data Files
- `processed_data.h5ad`: Processed AnnData object with all analysis results
- `cell_annotations_res{resolution}.csv`: Cell annotations for each resolution
  - Columns: `cell_barcode`, `cluster`, `cell_type`
  - One file per clustering resolution

### Plot Files (if `--save-plots` is used)
- `plots/qc_violin.png`: QC metrics violin plots
- `plots/qc_scatter.png`: Total counts vs genes scatter plot
- `plots/highly_variable_genes.png`: Highly variable genes dispersion
- `plots/pca_variance_ratio.png`: PCA explained variance
- `plots/umap_by_sample.png`: UMAP colored by sample (if available)
- `plots/umap_leiden_res{resolution}.png`: UMAP colored by clusters
- `plots/umap_celltype_res{resolution}.png`: UMAP colored by cell type
- `plots/marker_dotplot_res{resolution}.png`: Marker gene expression dotplot
- `plots/enrichment_scores_res{resolution}.png`: MLM enrichment score heatmap per cluster

## Example Workflow

1. Prepare your h5 data file
2. (Optional) Create a marker gene CSV file
3. Run the tool:
   ```bash
   python scanpy_cluster.py \
     --input compiled.h5 \
     --output-dir results/ \
     --markers example_markers.csv \
     --leiden-resolution 0.2,0.5,1.0 \
     --save-plots
   ```
4. Review the output:
   - Check QC plots in `results/plots/`
   - Load `results/processed_data.h5ad` in Python for further analysis
   - Use `results/cell_annotations_res*.csv` for downstream applications

## Features

### Quality Control
- Optional downsampling for quick testing on subset of cells
- Calculates mitochondrial, ribosomal, and hemoglobin gene percentages
- Filters low-quality cells and lowly-expressed genes

### Preprocessing
- Median normalization (log1p transformation)
- Highly variable gene selection (batch-aware if multiple samples)

### Analysis
- PCA dimensionality reduction
- Nearest neighbor graph construction
- UMAP embedding for visualization
- Leiden clustering at customizable resolutions

### Cell Type Annotation
- Marker gene-based annotation using decoupler's MLM (multivariate linear model)
- Enrichment analysis approach tests if marker gene sets are enriched in cells
- Uses decoupler's rankby_group to assign cell type with highest enrichment per cluster
- Supports multiple clustering resolutions

### Logging
- Detailed progress logging with timestamps
- Summary statistics at each step
- Error handling with informative messages

### Resume Capability
- Can resume from existing `processed_data.h5ad` file
- Skips already-computed steps (QC, normalization, PCA, UMAP, clustering, annotation)
- Useful for adding new clustering resolutions or annotations
- Saves time when experimenting with parameters

## Requirements

See `requirements.txt` for full dependencies. Key packages:
- scanpy >= 1.9.0
- anndata >= 0.9.0
- decoupler >= 1.4.0
- pandas >= 1.5.0
- numpy >= 1.23.0
- matplotlib >= 3.6.0
- igraph >= 0.10.0
- leidenalg >= 0.9.0

