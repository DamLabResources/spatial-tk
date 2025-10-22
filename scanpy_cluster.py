#!/usr/bin/env python3
"""
Scanpy Single-Cell RNA-seq Clustering and Annotation Tool

This script implements a complete scRNA-seq analysis workflow including:
- Quality control and filtering
- Normalization and feature selection
- Dimensionality reduction (PCA, UMAP)
- Leiden clustering at multiple resolutions
- Differential expression analysis for each cluster
- Cell type annotation using marker genes

The workflow is based on the Scverse basic tutorial:
https://scverse-tutorials.readthedocs.io/en/latest/notebooks/basic-scrna-tutorial.html

Example Usage:
    # Basic usage with default parameters
    python scanpy_cluster.py --input data.h5 --output-dir results/

    # With marker-based annotation and plots
    python scanpy_cluster.py --input data.h5 --output-dir results/ \\
        --markers cell_markers.csv --save-plots

    # Multiple clustering resolutions
    python scanpy_cluster.py --input data.h5 --output-dir results/ \\
        --leiden-resolution 0.2,0.5,1.0 --markers markers.csv

    # Downsample to 10% of cells for quick testing
    python scanpy_cluster.py --input data.h5 --output-dir results/ \\
        --downsample 0.1 --save-plots

    # Resume from existing analysis and add new annotations
    python scanpy_cluster.py --input data.h5 --output-dir results/ \\
        --markers updated_markers.csv --resume --save-plots

Marker CSV Format:
    cell_type,gene
    B cells,MS4A1
    B cells,CD19
    T cells,CD3D
    T cells,CD3E
"""
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

import argparse
import logging
import os
import sys
import time

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for saving plots
import matplotlib.pyplot as plt

# Suppress FutureWarnings


# Configure scanpy settings
sc.settings.verbosity = 1  # Reduce scanpy's verbosity
sc.settings.set_figure_params(dpi=100, facecolor="white", frameon=False)


def setup_logging():
    """Configure logging to show INFO level messages with timestamps."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def load_data(input_path: str) -> ad.AnnData:
    """
    Load single-cell data from h5 file.
    
    Args:
        input_path: Path to h5 or h5ad file
        
    Returns:
        AnnData object containing the count matrix
    """
    logging.info(f"Loading data from {input_path}")
    
    adata = ad.read_h5ad(input_path)
    
    # Make variable names unique
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    
    logging.info(f"Loaded data: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


def downsample_cells(adata: ad.AnnData, fraction: float) -> ad.AnnData:
    """
    Randomly downsample cells to specified fraction.
    
    Args:
        adata: AnnData object
        fraction: Fraction of cells to keep (0-1)
        
    Returns:
        Downsampled AnnData object
    """
    if fraction <= 0 or fraction > 1:
        raise ValueError("Downsample fraction must be between 0 and 1")
    
    if fraction == 1.0:
        logging.info("No downsampling (fraction=1.0)")
        return adata
    
    n_cells_original = adata.n_obs
    n_cells_keep = int(n_cells_original * fraction)
    
    logging.info(f"Downsampling from {n_cells_original} to {n_cells_keep} cells (fraction={fraction})")
    
    # Randomly sample cells
    sc.pp.subsample(adata, fraction=fraction)
    
    logging.info(f"Downsampled: {adata.n_obs} cells remaining")
    
    return adata


def calculate_qc_metrics(adata: ad.AnnData, resume: bool = False) -> ad.AnnData:
    """
    Calculate quality control metrics including mitochondrial, ribosomal,
    and hemoglobin gene percentages.
    
    Args:
        adata: AnnData object
        resume: If True, skip if QC metrics already exist
        
    Returns:
        AnnData object with QC metrics added
    """
    if resume and "pct_counts_mt" in adata.obs.columns:
        logging.info("QC metrics already calculated (resuming)")
        return adata
    
    logging.info("Calculating QC metrics")
    
    # Identify mitochondrial genes (MT- for human, Mt- for mouse)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    
    # Identify ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    
    # Identify hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        inplace=True,
        log1p=True
    )
    
    logging.info(f"QC metrics calculated - Median genes/cell: {np.median(adata.obs['n_genes_by_counts']):.0f}")
    logging.info(f"Median UMI/cell: {np.median(adata.obs['total_counts']):.0f}")
    logging.info(f"Median MT%: {np.median(adata.obs['pct_counts_mt']):.2f}%")
    
    return adata


def filter_cells_and_genes(
    adata: ad.AnnData,
    min_genes: int = 100,
    min_cells: int = 3
) -> ad.AnnData:
    """
    Filter cells and genes based on minimum thresholds.
    
    Args:
        adata: AnnData object
        min_genes: Minimum number of genes expressed per cell
        min_cells: Minimum number of cells expressing a gene
        
    Returns:
        Filtered AnnData object
    """
    logging.info(f"Filtering cells (min_genes={min_genes}) and genes (min_cells={min_cells})")
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars
    
    # Filter cells with too few genes
    sc.pp.filter_cells(adata, min_genes=min_genes)
    
    # Filter genes expressed in too few cells
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    logging.info(f"Filtered: {n_cells_before - adata.n_obs} cells, {n_genes_before - adata.n_vars} genes")
    logging.info(f"Remaining: {adata.n_obs} cells × {adata.n_vars} genes")
    
    return adata




def normalize_and_log(adata: ad.AnnData, resume: bool = False) -> ad.AnnData:
    """
    Normalize to median total counts and apply log transformation.
    
    Args:
        adata: AnnData object
        resume: If True, skip if normalization already done
        
    Returns:
        Normalized AnnData object
    """
    if resume and "counts" in adata.layers:
        logging.info("Normalization already done (resuming)")
        return adata
    
    logging.info("Normalizing and log-transforming data")
    
    # Save raw counts in layers
    adata.layers["counts"] = adata.X.copy()
    
    # Normalize to median total counts
    sc.pp.normalize_total(adata)
    
    # Log transform
    sc.pp.log1p(adata)
    
    logging.info("Normalization complete")
    return adata


def select_variable_genes(adata: ad.AnnData, n_top_genes: int = 2000, resume: bool = False) -> ad.AnnData:
    """
    Select highly variable genes for downstream analysis.
    
    Args:
        adata: AnnData object
        n_top_genes: Number of highly variable genes to select
        resume: If True, skip if highly variable genes already computed
        
    Returns:
        AnnData object with highly variable genes annotated
    """
    if resume and "highly_variable" in adata.var.columns:
        logging.info(f"Highly variable genes already selected (resuming)")
        return adata
    
    logging.info(f"Selecting {n_top_genes} highly variable genes")
    
    # Use batch_key if available for batch-aware feature selection
    batch_key = "sample" if "sample" in adata.obs.columns else None
    
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        batch_key=batch_key
    )
    
    n_hvg = adata.var["highly_variable"].sum()
    logging.info(f"Selected {n_hvg} highly variable genes")
    
    return adata


def run_pca(adata: ad.AnnData, resume: bool = False) -> ad.AnnData:
    """
    Perform PCA dimensionality reduction.
    
    Args:
        adata: AnnData object
        resume: If True, skip if PCA already computed
        
    Returns:
        AnnData object with PCA computed
    """
    if resume and "X_pca" in adata.obsm:
        logging.info("PCA already computed (resuming)")
        return adata
    
    logging.info("Running PCA")
    sc.tl.pca(adata)
    logging.info("PCA complete")
    return adata


def compute_neighbors_and_umap(adata: ad.AnnData, resume: bool = False) -> ad.AnnData:
    """
    Compute neighborhood graph and UMAP embedding.
    
    Args:
        adata: AnnData object
        resume: If True, skip if neighbors and UMAP already computed
        
    Returns:
        AnnData object with neighbors and UMAP computed
    """
    if resume and "X_umap" in adata.obsm:
        logging.info("Neighbors and UMAP already computed (resuming)")
        return adata
    
    logging.info("Computing neighborhood graph")
    sc.pp.neighbors(adata)
    
    logging.info("Computing UMAP embedding")
    sc.tl.umap(adata)
    
    logging.info("Neighbors and UMAP complete")
    return adata


def cluster_leiden(
    adata: ad.AnnData,
    resolution: float,
    key_added: str = "leiden",
    resume: bool = False
) -> ad.AnnData:
    """
    Perform Leiden clustering at specified resolution.
    
    Args:
        adata: AnnData object
        resolution: Clustering resolution parameter
        key_added: Key name for storing clustering results in adata.obs
        resume: If True, skip if clustering already exists
        
    Returns:
        AnnData object with clustering results added
    """
    if resume and key_added in adata.obs.columns:
        n_clusters = adata.obs[key_added].nunique()
        logging.info(f"Leiden clustering (resolution={resolution}) already exists with {n_clusters} clusters (resuming)")
        return adata
    
    logging.info(f"Running Leiden clustering (resolution={resolution})")
    
    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added=key_added,
        flavor="igraph",
        n_iterations=2,
        directed=False
    )
    
    n_clusters = adata.obs[key_added].nunique()
    logging.info(f"Found {n_clusters} clusters at resolution {resolution}")
    
    return adata


def load_marker_genes(marker_path: str) -> Dict[str, List[str]]:
    """
    Load marker genes from CSV file.
    
    Args:
        marker_path: Path to CSV file with columns: cell_type, gene
        
    Returns:
        Dictionary mapping cell type to list of marker genes
    """
    logging.info(f"Loading marker genes from {marker_path}")
    
    df = pd.read_csv(marker_path)
    
    if not all(col in df.columns for col in ["cell_type", "gene"]):
        raise ValueError("Marker CSV must have 'cell_type' and 'gene' columns")
    
    # Group by cell type
    markers = df.groupby("cell_type")["gene"].apply(list).to_dict()
    
    total_markers = sum(len(genes) for genes in markers.values())
    logging.info(f"Loaded {len(markers)} cell types with {total_markers} total marker genes")
    
    return markers


def annotate_with_markers(
    adata: ad.AnnData,
    markers: Dict[str, List[str]],
    cluster_key: str = "leiden",
    annotation_key: str = "cell_type",
    resume: bool = False
) -> ad.AnnData:
    """
    Annotate clusters with cell types based on marker gene expression using
    decoupler's multivariate linear model (MLM) approach.
    
    This method uses enrichment analysis to test if marker gene collections
    are enriched in cells, similar to the approach in the Scverse tutorial.
    
    Args:
        adata: AnnData object
        markers: Dictionary mapping cell type to list of marker genes
        cluster_key: Key in adata.obs containing cluster assignments
        annotation_key: Key name for storing cell type annotations
        resume: If True, skip if annotation already exists
        
    Returns:
        AnnData object with cell type annotations added
    """
    if resume and annotation_key in adata.obs.columns:
        logging.info(f"Cell type annotation already exists (resuming)")
        return adata
    
    logging.info(f"Annotating cell types using decoupler MLM (cluster_key={cluster_key})")
    
    # Check which marker genes are present
    all_marker_genes = set()
    for genes in markers.values():
        all_marker_genes.update(genes)
    
    missing_genes = all_marker_genes - set(adata.var_names)
    if missing_genes:
        logging.info(f"Note: {len(missing_genes)} marker genes not found in dataset")
    
    # Convert markers dictionary to DataFrame format expected by decoupler
    # Format: columns 'source' (cell_type) and 'target' (gene)
    marker_rows = []
    for cell_type, genes in markers.items():
        for gene in genes:
            marker_rows.append({"source": cell_type, "target": gene})
    
    marker_df = pd.DataFrame(marker_rows)
    
    # Add weight column (all weights = 1)
    marker_df["weight"] = 1
    
    logging.info(f"Running MLM with {len(marker_df)} marker gene entries across {len(markers)} cell types")
    
    # Run multivariate linear model
    # This calculates enrichment scores for each cell type in each cell
    dc.mt.mlm(adata, net=marker_df, verbose=False)
    
    # Extract the MLM scores from adata.obsm
    # This creates a new AnnData-like object with cells x cell_types
    acts = dc.pp.get_obsm(adata, "score_mlm")
    
    # For each cluster, find the cell type with highest enrichment score
    # Use decoupler's rankby_group to get top scoring cell type per cluster
    enr = dc.tl.rankby_group(acts, groupby=cluster_key)
    
    # Get the top cell type (highest stat) for each cluster
    # Filter to positive stats only (stat > 0)
    annotation_dict = (
        enr[enr["stat"] > 0]
        .groupby("group", observed=True)
        .head(1)
        .set_index("group")["name"]
        .to_dict()
    )
    
    # Handle clusters that may not have positive enrichment scores
    all_clusters = adata.obs[cluster_key].unique()
    for cluster in all_clusters:
        if cluster not in annotation_dict:
            annotation_dict[cluster] = "Unknown"
    
    # Map cluster annotations to cells
    adata.obs[annotation_key] = adata.obs[cluster_key].map(annotation_dict)
    adata.obs[annotation_key] = adata.obs[annotation_key].astype("category")
    
    # Log annotation summary
    annotation_counts = adata.obs[annotation_key].value_counts()
    logging.info("Cell type annotation summary:")
    for cell_type, count in annotation_counts.items():
        logging.info(f"  {cell_type}: {count} cells")
    
    return adata


def create_enrichment_dotplot(
    adata: ad.AnnData,
    cluster_key: str,
    output_path: Path
) -> None:
    """
    Create a dotplot showing MLM enrichment scores per cluster for each cell type.
    
    Args:
        adata: AnnData object with MLM scores in obsm
        cluster_key: Key in adata.obs containing cluster assignments
        output_path: Path to save the plot
    """
    if "score_mlm" not in adata.obsm:
        logging.warning("MLM scores not found in adata.obsm, skipping enrichment dotplot")
        return
    
    try:
        # Extract MLM scores directly from obsm
        mlm_scores = adata.obsm["score_mlm"]
        
        # Get cell types (column names from MLM scores)
        if hasattr(mlm_scores, 'dtype') and hasattr(mlm_scores.dtype, 'names'):
            # Structured array
            cell_types = list(mlm_scores.dtype.names)
            score_matrix = np.column_stack([mlm_scores[ct] for ct in cell_types])
        else:
            # Regular array - need to get cell types from decoupler result
            acts = dc.pp.get_obsm(adata, "score_mlm")
            cell_types = acts.var_names.tolist()
            score_matrix = mlm_scores
        
        # Create DataFrame with scores
        scores_df = pd.DataFrame(
            score_matrix,
            index=adata.obs_names,
            columns=cell_types
        )
        
        # Add cluster information
        scores_df[cluster_key] = adata.obs[cluster_key].values
        
        # Calculate mean score per cluster
        cluster_means = scores_df.groupby(cluster_key).mean()
        
        # Ensure all values are numeric
        cluster_means = cluster_means.astype(float)
        
        # Create heatmap
        import seaborn as sns
        fig, ax = plt.subplots(figsize=(max(12, len(cell_types) * 0.6), max(6, len(cluster_means) * 0.4)))
        
        sns.heatmap(
            cluster_means, 
            cmap="RdBu_r", 
            center=0,
            cbar_kws={'label': 'Mean MLM Enrichment Score'},
            linewidths=0.5, 
            linecolor='lightgray',
            ax=ax,
            robust=True
        )
        
        ax.set_xlabel("Cell Type", fontsize=12)
        ax.set_ylabel("Cluster", fontsize=12)
        ax.set_title("Cell Type Enrichment Scores by Cluster", fontsize=14)
        
        # Rotate x labels for better readability
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        
        plt.tight_layout()
        plt.savefig(output_path, bbox_inches="tight", dpi=150)
        plt.close()
        logging.info(f"  Saved enrichment heatmap to {output_path}")
        
    except Exception as e:
        logging.warning(f"  Could not create enrichment dotplot: {e}")
        import traceback
        logging.debug(traceback.format_exc())


def save_plots(
    adata: ad.AnnData,
    output_dir: Path,
    resolutions: List[float],
    markers: Optional[Dict[str, List[str]]] = None
) -> None:
    """
    Generate and save QC and analysis plots.
    
    Args:
        adata: AnnData object
        output_dir: Directory to save plots
        resolutions: List of clustering resolutions used
        markers: Optional marker gene dictionary for dotplots
    """
    logging.info("Generating plots")
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)
    
    # QC violin plots
    try:
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
            show=False
        )
        plt.savefig(plots_dir / "qc_violin.png", bbox_inches="tight", dpi=150)
        plt.close()
        logging.info("  Saved QC violin plot")
    except Exception as e:
        logging.warning(f"  Could not generate QC violin plot: {e}")
    
    # QC scatter plot
    try:
        sc.pl.scatter(
            adata,
            "total_counts",
            "n_genes_by_counts",
            color="pct_counts_mt",
            show=False
        )
        plt.savefig(plots_dir / "qc_scatter.png", bbox_inches="tight", dpi=150)
        plt.close()
        logging.info("  Saved QC scatter plot")
    except Exception as e:
        logging.warning(f"  Could not generate QC scatter plot: {e}")
    
    # Highly variable genes
    try:
        sc.pl.highly_variable_genes(adata, show=False)
        plt.savefig(plots_dir / "highly_variable_genes.png", bbox_inches="tight", dpi=150)
        plt.close()
        logging.info("  Saved highly variable genes plot")
    except Exception as e:
        logging.warning(f"  Could not generate highly variable genes plot: {e}")
    
    # PCA variance ratio
    try:
        sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
        plt.savefig(plots_dir / "pca_variance_ratio.png", bbox_inches="tight", dpi=150)
        plt.close()
        logging.info("  Saved PCA variance ratio plot")
    except Exception as e:
        logging.warning(f"  Could not generate PCA variance plot: {e}")
    
    # UMAP by sample if available
    if "sample" in adata.obs.columns:
        try:
            sc.pl.umap(adata, color="sample", show=False)
            plt.savefig(plots_dir / "umap_by_sample.png", bbox_inches="tight", dpi=150)
            plt.close()
            logging.info("  Saved UMAP by sample plot")
        except Exception as e:
            logging.warning(f"  Could not generate UMAP by sample plot: {e}")
    
    # UMAP by leiden clustering and annotations for each resolution
    for resolution in resolutions:
        res_str = str(resolution).replace(".", "p")
        cluster_key = f"leiden_res{res_str}"
        annotation_key = f"cell_type_res{res_str}"
        
        # UMAP by cluster
        if cluster_key in adata.obs.columns:
            try:
                sc.pl.umap(adata, color=cluster_key, legend_loc="on data", show=False)
                plt.savefig(plots_dir / f"umap_leiden_res{res_str}.png", bbox_inches="tight", dpi=150)
                plt.close()
                logging.info(f"  Saved UMAP leiden resolution {resolution} plot")
            except Exception as e:
                logging.warning(f"  Could not generate UMAP leiden plot: {e}")
        
        # UMAP by cell type annotation
        if annotation_key in adata.obs.columns:
            try:
                sc.pl.umap(adata, color=annotation_key, show=False)
                plt.savefig(plots_dir / f"umap_celltype_res{res_str}.png", bbox_inches="tight", dpi=150)
                plt.close()
                logging.info(f"  Saved UMAP cell type annotation resolution {resolution} plot")
            except Exception as e:
                logging.warning(f"  Could not generate UMAP cell type plot: {e}")
        
        # Dotplot of markers if available
        if markers and cluster_key in adata.obs.columns:
            try:
                # Filter markers to only include genes present in the dataset
                filtered_markers = {}
                for cell_type, genes in markers.items():
                    present_genes = [g for g in genes if g in adata.var_names]
                    if present_genes:
                        filtered_markers[cell_type] = present_genes
                
                if filtered_markers:
                    sc.pl.dotplot(
                        adata,
                        filtered_markers,
                        groupby=cluster_key,
                        show=False
                    )
                    plt.savefig(plots_dir / f"marker_dotplot_res{res_str}.png", bbox_inches="tight", dpi=150)
                    plt.close()
                    logging.info(f"  Saved marker dotplot resolution {resolution}")
                else:
                    logging.warning(f"  No marker genes found in dataset for dotplot")
            except Exception as e:
                logging.warning(f"  Could not generate marker dotplot: {e}")
        
        # Enrichment score heatmap if annotation was done
        if annotation_key in adata.obs.columns and "score_mlm" in adata.obsm:
            enrichment_path = plots_dir / f"enrichment_scores_res{res_str}.png"
            create_enrichment_dotplot(adata, cluster_key, enrichment_path)
        
        # Differential expression dotplot
        rank_key = f"rank_genes_{cluster_key}"
        if rank_key in adata.uns:
            try:
                # Create dotplot of top differentially expressed genes
                sc.pl.rank_genes_groups_dotplot(
                    adata,
                    n_genes=5,  # Top 5 genes per cluster
                    key=rank_key,
                    groupby=cluster_key,
                    show=False,
                    dendrogram=False
                )
                plt.savefig(plots_dir / f"deg_dotplot_res{res_str}.png", bbox_inches="tight", dpi=150)
                plt.close()
                logging.info(f"  Saved differential expression dotplot resolution {resolution}")
            except Exception as e:
                logging.warning(f"  Could not generate differential expression dotplot: {e}")
            
            # Also create a heatmap of top genes
            try:
                sc.pl.rank_genes_groups_heatmap(
                    adata,
                    n_genes=10,  # Top 10 genes per cluster
                    key=rank_key,
                    groupby=cluster_key,
                    show=False,
                    show_gene_labels=True
                )
                plt.savefig(plots_dir / f"deg_heatmap_res{res_str}.png", bbox_inches="tight", dpi=150)
                plt.close()
                logging.info(f"  Saved differential expression heatmap resolution {resolution}")
            except Exception as e:
                logging.warning(f"  Could not generate differential expression heatmap: {e}")
    
    logging.info(f"All plots saved to {plots_dir}")


def run_differential_expression(
    adata: ad.AnnData,
    cluster_key: str,
    method: str = "wilcoxon",
    resume: bool = False
) -> ad.AnnData:
    """
    Run differential expression analysis to find marker genes for each cluster.
    
    Args:
        adata: AnnData object
        cluster_key: Key in adata.obs containing cluster assignments
        method: Statistical test to use (default: wilcoxon)
        resume: If True, skip if differential expression already computed
        
    Returns:
        AnnData object with differential expression results added
    """
    rank_key = f"rank_genes_{cluster_key}"
    
    if resume and "rank_genes_groups" in adata.uns and adata.uns.get("rank_genes_groups_key") == rank_key:
        logging.info(f"Differential expression already computed for {cluster_key} (resuming)")
        return adata
    
    logging.info(f"Running differential expression analysis for {cluster_key}")
    
    # Run rank_genes_groups
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method=method,
        use_raw=False,  # Use normalized data in .X
        key_added=rank_key,
        layer=None
    )
    
    # Store which key was used
    adata.uns["rank_genes_groups_key"] = rank_key
    
    n_clusters = adata.obs[cluster_key].nunique()
    logging.info(f"  Differential expression completed for {n_clusters} clusters")
    
    return adata


def save_differential_expression_results(
    adata: ad.AnnData,
    cluster_key: str,
    output_dir: Path,
    n_genes: int = 100
) -> None:
    """
    Save differential expression results to CSV files.
    
    Args:
        adata: AnnData object with differential expression results
        cluster_key: Key in adata.obs containing cluster assignments
        output_dir: Directory to save output files
        n_genes: Number of top genes to save per cluster
    """
    rank_key = f"rank_genes_{cluster_key}"
    
    if rank_key not in adata.uns:
        logging.warning(f"  No differential expression results found for {cluster_key}")
        return
    
    logging.info(f"  Saving differential expression results for {cluster_key}")
    
    # Get the differential expression results as a DataFrame
    result = sc.get.rank_genes_groups_df(adata, group=None, key=rank_key)
    
    # Save all results
    de_dir = output_dir / "differential_expression"
    de_dir.mkdir(exist_ok=True)
    
    res_str = cluster_key.replace("leiden_res", "")
    
    # Save complete results
    all_results_path = de_dir / f"deg_all_clusters_res{res_str}.csv"
    result.to_csv(all_results_path, index=False)
    logging.info(f"    Saved all DE genes to {all_results_path}")
    
    # Save top N genes per cluster
    top_results_path = de_dir / f"deg_top{n_genes}_per_cluster_res{res_str}.csv"
    top_result = result.groupby('group').head(n_genes)
    top_result.to_csv(top_results_path, index=False)
    logging.info(f"    Saved top {n_genes} DE genes per cluster to {top_results_path}")


def save_results(
    adata: ad.AnnData,
    output_dir: Path,
    resolutions: List[float]
) -> None:
    """
    Save processed AnnData object and cell type annotation CSVs.
    
    Args:
        adata: Processed AnnData object
        output_dir: Directory to save output files
        resolutions: List of clustering resolutions used
    """
    logging.info(f"Saving results to {output_dir}")
    
    # Save h5ad file
    h5ad_path = output_dir / "processed_data.h5ad"
    adata.write_h5ad(h5ad_path)
    logging.info(f"  Saved AnnData object to {h5ad_path}")
    
    # Save annotation CSVs for each resolution
    for resolution in resolutions:
        res_str = str(resolution).replace(".", "p")
        cluster_key = f"leiden_res{res_str}"
        annotation_key = f"cell_type_res{res_str}"
        
        if annotation_key in adata.obs.columns:
            # Create DataFrame with cell barcode, cluster, and cell type
            df = pd.DataFrame({
                "cell_barcode": adata.obs_names,
                "cluster": adata.obs[cluster_key],
                "cell_type": adata.obs[annotation_key]
            })
            
            csv_path = output_dir / f"cell_annotations_res{res_str}.csv"
            df.to_csv(csv_path, index=False)
            logging.info(f"  Saved annotations (resolution {resolution}) to {csv_path}")
        
        # Save differential expression results
        save_differential_expression_results(adata, cluster_key, output_dir)
    
    logging.info("Results saved successfully")


def parse_resolutions(resolution_str: str) -> List[float]:
    """
    Parse comma-separated resolution string into list of floats.
    
    Args:
        resolution_str: Comma-separated string of resolutions (e.g., "0.2,0.5,1.0")
        
    Returns:
        List of resolution values
    """
    resolutions = []
    for res in resolution_str.split(","):
        try:
            resolutions.append(float(res.strip()))
        except ValueError:
            raise ValueError(f"Invalid resolution value: {res}")
    return resolutions


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="scRNA-seq clustering and cell type annotation tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  %(prog)s --input data.h5 --output-dir results/
  %(prog)s --input data.h5 --output-dir results/ --markers markers.csv --save-plots
  %(prog)s --input data.h5 --output-dir results/ --leiden-resolution 0.2,0.5,1.0
  %(prog)s --input data.h5 --output-dir results/ --downsample 0.1 --save-plots
  %(prog)s --input data.h5 --output-dir results/ --markers new_markers.csv --resume
        """
    )
    
    # Required arguments
    parser.add_argument(
        "--input",
        required=True,
        help="Path to input h5 or h5ad file"
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory path where all output files will be placed"
    )
    
    # Optional arguments
    parser.add_argument(
        "--markers",
        help="Path to CSV file with marker genes (columns: cell_type, gene)"
    )
    parser.add_argument(
        "--save-plots",
        action="store_true",
        help="Generate and save QC and analysis plots"
    )
    parser.add_argument(
        "--min-genes",
        type=int,
        default=100,
        help="Minimum number of genes expressed per cell (default: 100)"
    )
    parser.add_argument(
        "--min-cells",
        type=int,
        default=3,
        help="Minimum number of cells expressing a gene (default: 3)"
    )
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=2000,
        help="Number of highly variable genes to select (default: 2000)"
    )
    parser.add_argument(
        "--leiden-resolution",
        type=str,
        default="0.5",
        help="Leiden clustering resolution(s), comma-separated for multiple (default: 0.5)"
    )
    parser.add_argument(
        "--downsample",
        type=float,
        default=1.0,
        help="Fraction of cells to keep (0-1, default: 1.0 = no downsampling)"
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from existing analysis, skip already computed steps"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Parse resolutions
    try:
        resolutions = parse_resolutions(args.leiden_resolution)
    except ValueError as e:
        logging.error(f"Error parsing resolutions: {e}")
        sys.exit(1)
    
    logging.info("="*60)
    logging.info("scRNA-seq Clustering and Annotation Pipeline")
    logging.info("="*60)
    
    start_time = time.time()
    
    try:
        # Load data - check for existing processed file if resuming
        processed_path = output_dir / "processed_data.h5ad"
        if args.resume and processed_path.exists():
            logging.info(f"Resuming from existing file: {processed_path}")
            adata = ad.read_h5ad(processed_path)
            logging.info(f"Loaded data: {adata.n_obs} cells × {adata.n_vars} genes")
        else:
            adata = load_data(args.input)
        
        # Downsample if requested
        if args.downsample < 1.0:
            adata = downsample_cells(adata, args.downsample)
        
        # QC and filtering
        adata = calculate_qc_metrics(adata, resume=args.resume)
        adata = filter_cells_and_genes(adata, args.min_genes, args.min_cells)
        
        # Normalization and feature selection
        adata = normalize_and_log(adata, resume=args.resume)
        adata = select_variable_genes(adata, args.n_top_genes, resume=args.resume)
        
        # Dimensionality reduction
        adata = run_pca(adata, resume=args.resume)
        adata = compute_neighbors_and_umap(adata, resume=args.resume)
        
        # Clustering at multiple resolutions
        for resolution in resolutions:
            res_str = str(resolution).replace(".", "p")
            cluster_key = f"leiden_res{res_str}"
            adata = cluster_leiden(adata, resolution, key_added=cluster_key, resume=args.resume)
        
        # Run differential expression analysis for each resolution
        for resolution in resolutions:
            res_str = str(resolution).replace(".", "p")
            cluster_key = f"leiden_res{res_str}"
            adata = run_differential_expression(adata, cluster_key, resume=args.resume)
        
        # Cell type annotation if markers provided
        markers = None
        if args.markers:
            markers = load_marker_genes(args.markers)
            
            for resolution in resolutions:
                res_str = str(resolution).replace(".", "p")
                cluster_key = f"leiden_res{res_str}"
                annotation_key = f"cell_type_res{res_str}"
                adata = annotate_with_markers(
                    adata,
                    markers,
                    cluster_key=cluster_key,
                    annotation_key=annotation_key,
                    resume=args.resume
                )
        
        # Save results
        save_results(adata, output_dir, resolutions)
        
        # Generate plots if requested
        if args.save_plots:
            save_plots(adata, output_dir, resolutions, markers)
        
        # Summary
        elapsed_time = time.time() - start_time
        logging.info("="*60)
        logging.info(f"Pipeline completed successfully in {elapsed_time:.1f} seconds")
        logging.info(f"Final dataset: {adata.n_obs} cells × {adata.n_vars} genes")
        logging.info(f"Results saved to: {output_dir}")
        logging.info("="*60)
        
    except Exception as e:
        logging.error(f"Pipeline failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()

