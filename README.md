# spatial-tk: Spatial Transcriptomics Analysis Toolkit

A comprehensive modular Python toolkit for Xenium spatial transcriptomics analysis. This package provides a command-line interface with separate subcommands for each stage of the analysis pipeline, enabling flexible and efficient processing of spatial transcriptomics data.

## Features

- **Modular Pipeline**: Each analysis step is a separate command for maximum flexibility
- **Inplace Processing**: Optionally modify datasets without duplication to save disk space
- **Multiple Resolutions**: Support for multi-resolution clustering analysis
- **Rich Annotations**: Marker-based and MLM enrichment-based cell type annotation
- **Flexible Differential Analysis**: Compare groups or find cluster markers
- **Configuration Files**: TOML config files for reproducible pipelines
- **Comprehensive Testing**: Unit and functional tests ensure reliability

## Installation

### From Source

```bash
# Clone or navigate to the repository
cd clustering-tools

# Install in development mode (recommended)
pip install -e .

# Or install normally
pip install .
```

### Dependencies

The package requires Python ≥3.9 and includes dependencies for:
- Spatial data handling (spatialdata, squidpy)
- Single-cell analysis (scanpy, anndata)
- Cell type annotation (decoupler)
- Visualization (matplotlib, seaborn)
- Clustering (igraph, leidenalg)

See `pyproject.toml` for complete dependency list.

## Configuration Files

All commands support TOML configuration files for reproducible pipelines. Each command has its own section in the config file, and CLI arguments override config values when both are provided.

### Basic Usage

```bash
# Use a config file
spatial-tk concat --config config.toml --input samples.csv --output merged.zarr
spatial-tk normalize --config config.toml --input merged.zarr --inplace
```

### Config File Format

Create a `config.toml` file with sections for each command:

```toml
[concat]
input = "samples.csv"
output = "merged.zarr"
downsample = 1.0

[normalize]
input = "merged.zarr"
inplace = true
min_genes = 100
min_cells = 3
n_top_genes = 2000
save_plots = false

[cluster]
input = "merged.zarr"
inplace = true
leiden_resolution = "0.2,0.5,1.0"
save_plots = true

[annotate]
input = "merged.zarr"
inplace = true
markers = "markers.csv"
calculate_ulm = true
panglao_min_sensitivity = 0.5
tmin = 2
save_plots = true

[differential]
input = "merged.zarr"
output_dir = "results/"
groupby = "leiden_res0p5"
method = "wilcoxon"
n_genes = 100
save_plots = false
```

### Config Key Naming

Config keys use underscores (e.g., `min_genes`, `n_top_genes`), which correspond to CLI arguments with hyphens (`--min-genes`, `--n-top-genes`). The config system automatically handles this conversion.

### CLI Arguments Override Config

When both a config file and CLI arguments are provided, CLI arguments take precedence:

```bash
# Config specifies downsample = 0.5, but CLI overrides it to 0.8
spatial-tk concat --config config.toml --input samples.csv --output merged.zarr --downsample 0.8
```

### Example Config File

See `example_config.toml` in the repository root for a complete example with all available options documented.

```bash
# 1. Concatenate multiple samples
spatial-tk concat --input samples.csv --output merged.zarr

# 2. Normalize (inplace to save space)
spatial-tk normalize --input merged.zarr --inplace

# 3. Cluster with multiple resolutions
spatial-tk cluster --input merged.zarr --inplace --leiden-resolution 0.2,0.5,1.0

# 4. Annotate cell types
spatial-tk annotate --input merged.zarr --inplace --markers markers.csv

# 5. Differential expression analysis
spatial-tk differential --input merged.zarr --output-dir results/ --groupby leiden_res0p5
```

## Commands

### `spatial-tk concat`

Concatenate multiple Xenium .zarr files into a single dataset.

```bash
spatial-tk concat --input samples.csv --output merged.zarr

# With downsampling for testing
spatial-tk concat --input samples.csv --output merged.zarr --downsample 0.1
```

**Arguments:**
- `--input`: Path to CSV file with columns: `sample`, `path`, [optional metadata]
- `--output`: Path to output .zarr file
- `--downsample`: Fraction of cells to keep (0-1, default: 1.0)
- `--config`: Path to TOML configuration file (optional)

**CSV Format:**
```csv
sample,path,status,location
sample1,/path/to/sample1.zarr,HIV,Drexel
sample2,/path/to/sample2.zarr,NEG,OSU
```

### `spatial-tk normalize`

Perform QC, filtering, normalization, and feature selection.

```bash
# Save to new file
spatial-tk normalize --input data.zarr --output normalized.zarr

# Modify in place
spatial-tk normalize --input data.zarr --inplace

# With custom parameters and plots
spatial-tk normalize --input data.zarr --inplace \
  --min-genes 200 \
  --min-cells 5 \
  --n-top-genes 3000 \
  --save-plots
```

**Arguments:**
- `--input`: Input .zarr file
- `--output`: Output .zarr file (mutually exclusive with --inplace)
- `--inplace`: Modify input file in place
- `--min-genes`: Minimum genes per cell (default: 100)
- `--min-cells`: Minimum cells per gene (default: 3)
- `--n-top-genes`: Number of highly variable genes (default: 2000)
- `--save-plots`: Generate QC plots
- `--config`: Path to TOML configuration file (optional)

### `spatial-tk cluster`

Perform PCA, neighbor graph computation, UMAP, and Leiden clustering.

```bash
# Single resolution
spatial-tk cluster --input data.zarr --inplace --leiden-resolution 0.5

# Multiple resolutions with plots
spatial-tk cluster --input data.zarr --inplace \
  --leiden-resolution 0.2,0.5,1.0,2.0 \
  --save-plots
```

**Arguments:**
- `--input`: Input normalized .zarr file
- `--output`: Output .zarr file (mutually exclusive with --inplace)
- `--inplace`: Modify input file in place
- `--leiden-resolution`: Clustering resolution(s), comma-separated (default: 0.5)
- `--save-plots`: Generate UMAP plots
- `--config`: Path to TOML configuration file (optional)

### `spatial-tk spatial_neighbors`

Build a spatial graph on coordinates using Squidpy.

```bash
# kNN graph on obsm['spatial']
spatial-tk spatial_neighbors --input data.zarr --inplace \
  --spatial-key spatial --n-neighs 8

# Radius-based graph with cosine transform, writing to new file
spatial-tk spatial_neighbors --input data.zarr --output neighbors.zarr \
  --coord-type generic --radius 50,200 --transform cosine
```

**Arguments:**
- `--input`: Input .zarr file
- `--output`: Output .zarr file (mutually exclusive with --inplace)
- `--inplace`: Modify input file in place
- `--table-key`: Optional table key in `SpatialData.tables`
- `--spatial-key`: Coordinate key in `adata.obsm` (default: `spatial`)
- `--library-key`: Optional obs column containing library ids
- `--library-id`: Optional single-library convenience value
- `--coord-type`: `grid` or `generic` (default: inferred by Squidpy)
- `--n-neighs`: Number of neighbors (default: 6)
- `--radius`: Scalar radius or `min,max` interval
- `--transform`: `spectral`, `cosine`, or `none`
- `--key-added`: Output prefix in `adata.obsp`/`adata.uns` (default: `spatial`)
- `--config`: Path to TOML configuration file (optional)

### `spatial-tk annotate`

Annotate cell types using marker genes and/or MLM scoring.

```bash
# Basic annotation with markers
spatial-tk annotate --input data.zarr --inplace --markers markers.csv

# With MLM enrichment scores
spatial-tk annotate --input data.zarr --inplace \
  --markers markers.csv \
  --calculate-ulm \
  --save-plots

# Annotate specific clustering
spatial-tk annotate --input data.zarr --inplace \
  --markers markers.csv \
  --cluster-key leiden_res1p0
```

**Arguments:**
- `--input`: Input clustered .zarr file
- `--output`: Output .zarr file (mutually exclusive with --inplace)
- `--inplace`: Modify input file in place
- `--markers`: Path to marker genes CSV (columns: `cell_type`, `gene`)
- `--cluster-key`: Specific cluster column to annotate (default: all leiden_res*)
- `--calculate-ulm`: Calculate MLM enrichment scores for pathways/TFs
- `--panglao-min-sensitivity`: Min sensitivity for PanglaoDB markers (default: 0.5)
- `--tmin`: Minimum marker genes per cell type (default: 2)
- `--save-plots`: Generate annotation plots
- `--config`: Path to TOML configuration file (optional)

**MLM Resources:**
- **hallmark**: MSigDB Hallmark gene sets
- **collectri**: CollecTRI TF regulons
- **dorothea**: DoRothEA TF activities
- **progeny**: PROGENy pathway activities
- **PanglaoDB**: Filtered cell type markers

### `spatial-tk differential`

Differential expression analysis with two modes:

**Mode A**: Compare two specific groups (e.g., HIV vs NEG)
**Mode B**: Find marker genes for all groups/clusters

```bash
# Mode B: Find markers for all clusters
spatial-tk differential \
  --input data.zarr \
  --output-dir results/ \
  --groupby leiden_res0p5

# Mode A: Compare two groups
spatial-tk differential \
  --input data.zarr \
  --output-dir results/ \
  --groupby status \
  --compare-groups HIV,NEG

# With obsm enrichment scores
spatial-tk differential \
  --input data.zarr \
  --output-dir results/ \
  --groupby status \
  --compare-groups HIV,NEG \
  --obsm-layer score_mlm_PanglaoDB \
  --save-plots

# Compare cell types
spatial-tk differential \
  --input data.zarr \
  --output-dir results/ \
  --groupby cell_type_res0p5 \
  --n-genes 50
```

**Arguments:**
- `--input`: Input .zarr file with annotations
- `--output-dir`: Directory for results
- `--groupby`: Column in obs to group by (e.g., "leiden_res0p5", "status", "cell_type")
- `--compare-groups`: Two groups to compare (Mode A), comma-separated
- `--obsm-layer`: Optional obsm layer for enrichment analysis (e.g., "score_mlm_PanglaoDB")
- `--method`: Statistical test method (default: wilcoxon)
- `--layer`: Layer to use for expression (default: None uses .X)
- `--n-genes`: Number of top genes to save (default: 100)
- `--save-plots`: Generate differential analysis plots
- `--config`: Path to TOML configuration file (optional)

## Example Workflows

### Full Pipeline with Config File

```bash
# Create config.toml with your settings
# Then run pipeline with config
spatial-tk concat --config config.toml --input samples.csv --output data.zarr
spatial-tk normalize --config config.toml --input data.zarr --inplace
spatial-tk cluster --config config.toml --input data.zarr --inplace
spatial-tk annotate --config config.toml --input data.zarr --inplace
spatial-tk differential --config config.toml --input data.zarr --output-dir results/
```

### Full Pipeline (In-place to Save Space)

```bash
# Step 1: Concatenate samples
spatial-tk concat --input samples.csv --output data.zarr

# Step 2-5: Process in place
spatial-tk normalize --input data.zarr --inplace --save-plots
spatial-tk cluster --input data.zarr --inplace --leiden-resolution 0.5,1.0 --save-plots
spatial-tk annotate --input data.zarr --inplace --markers markers.csv --calculate-ulm --save-plots
spatial-tk differential --input data.zarr --output-dir results/ --groupby leiden_res0p5 --save-plots
```

### Separate Files for Each Step

```bash
spatial-tk concat --input samples.csv --output step1_concat.zarr
spatial-tk normalize --input step1_concat.zarr --output step2_normalized.zarr
spatial-tk cluster --input step2_normalized.zarr --output step3_clustered.zarr
spatial-tk annotate --input step3_clustered.zarr --output step4_annotated.zarr
spatial-tk differential --input step4_annotated.zarr --output-dir results/
```

### Compare Disease Status

```bash
# Process and normalize
spatial-tk concat --input samples.csv --output data.zarr
spatial-tk normalize --input data.zarr --inplace

# Compare HIV vs NEG
spatial-tk differential \
  --input data.zarr \
  --output-dir hiv_vs_neg/ \
  --groupby status \
  --compare-groups HIV,NEG \
  --save-plots
```

### Multi-Resolution Analysis

```bash
spatial-tk concat --input samples.csv --output data.zarr
spatial-tk normalize --input data.zarr --inplace
spatial-tk cluster --input data.zarr --inplace --leiden-resolution 0.2,0.5,1.0,2.0

# Annotate all resolutions
spatial-tk annotate --input data.zarr --inplace --markers markers.csv --save-plots

# Differential analysis for each resolution
for res in 0p2 0p5 1p0 2p0; do
  spatial-tk differential \
    --input data.zarr \
    --output-dir results_res${res}/ \
    --groupby leiden_res${res}
done
```

## Output Files

### Concat
- `{output}.zarr`: Concatenated spatial dataset

### Normalize
- `{output}.zarr`: Normalized dataset with QC metrics
- `plots/qc_*.png`: QC plots (if --save-plots)

### Cluster
- `{output}.zarr`: Dataset with clustering results
- `plots/umap_leiden_res*.png`: UMAP plots (if --save-plots)

### Annotate
- `{output}.zarr`: Dataset with cell type annotations
- `plots/umap_celltype_res*.png`: Annotated UMAP plots (if --save-plots)
- `plots/marker_dotplot_res*.png`: Marker expression dotplots
- `plots/deg_*.png`: Differential expression plots

### Differential
- `de_genes_*.csv`: Differential expression results
- `de_{obsm_layer}_*.csv`: obsm enrichment results (if --obsm-layer used)
- `plots/`: Visualization plots (if --save-plots)

## Development

### Running Tests

```bash
# Install with dev dependencies
pip install -e ".[dev]"

# Run all tests
pytest

# Run only unit tests (fast)
pytest tests/unit/

# Run only functional tests (slower, requires test data)
pytest tests/functional/

# Run with coverage
pytest --cov=spatial_tk --cov-report=html
```

### Creating Test Data

```bash
python scripts/create_test_data.py \
  --input-csv example.csv \
  --output-dir tests/test_data \
  --n-cells 500
```

### Building Package

```bash
# Build distribution
python -m build

# Install locally
pip install dist/spatial_tk-*.whl
```

## Marker Gene CSV Format

```csv
cell_type,gene
T cells,CD3D
T cells,CD3E
B cells,MS4A1
B cells,CD19
Macrophages,CD68
Macrophages,CD14
```

## Advanced Usage

### Python API

The package can also be used programmatically:

```python
from spatial_tk.core import data_io, preprocessing, clustering, annotation
from spatial_tk.utils.helpers import get_table, set_table

# Load data
sdata = data_io.load_existing_spatial_data("data.zarr")
adata = get_table(sdata)

# Process
adata = preprocessing.normalize_and_log(adata)
adata = clustering.run_pca(adata)
adata = clustering.compute_neighbors_and_umap(adata)
adata = clustering.cluster_leiden(adata, resolution=0.5)

# Save
set_table(sdata, adata)
data_io.save_spatial_data(sdata, "processed.zarr")
```

## Citation

This tool is based on the Scverse ecosystem and follows best practices from:
- [Scverse Basic Tutorial](https://scverse-tutorials.readthedocs.io/en/latest/notebooks/basic-scrna-tutorial.html)
- [Decoupler documentation](https://decoupler.readthedocs.io/)

## License

MIT License

## Support

For issues, questions, or contributions, please contact the Hope Lab or open an issue on GitHub.
