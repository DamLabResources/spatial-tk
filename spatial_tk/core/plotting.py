#!/usr/bin/env python3
"""
Visualization functions for spatial transcriptomics analysis.

This module handles generation of QC plots, UMAP visualizations,
marker dotplots, and enrichment heatmaps.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional

import anndata as ad
import decoupler as dc
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for saving plots
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


def create_plots_directory(output_dir: Path) -> Path:
    """
    Create plots directory if it doesn't exist.
    
    Args:
        output_dir: Base output directory
        
    Returns:
        Path to plots directory
    """
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)
    return plots_dir


def save_qc_plots(adata: ad.AnnData, plots_dir: Path) -> None:
    """
    Generate and save QC plots (violin and scatter plots).
    
    Args:
        adata: AnnData object with QC metrics
        plots_dir: Directory to save plots
    """
    logging.info("Generating QC plots")
    
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


def save_umap_plots(
    adata: ad.AnnData,
    plots_dir: Path,
    cluster_key: str,
    annotation_key: Optional[str] = None,
    resolution: float = None
) -> None:
    """
    Generate and save UMAP plots colored by sample, cluster, and cell type.
    
    Args:
        adata: AnnData object with UMAP coordinates
        plots_dir: Directory to save plots
        cluster_key: Key for cluster assignments in adata.obs
        annotation_key: Optional key for cell type annotations in adata.obs
        resolution: Clustering resolution (for filename)
    """
    res_str = str(resolution).replace(".", "p") if resolution is not None else ""
    
    # UMAP by sample if available
    if "sample" in adata.obs.columns:
        try:
            sc.pl.umap(adata, color="sample", show=False)
            plt.savefig(plots_dir / "umap_by_sample.png", bbox_inches="tight", dpi=150)
            plt.close()
            logging.info("  Saved UMAP by sample plot")
        except Exception as e:
            logging.warning(f"  Could not generate UMAP by sample plot: {e}")
    
    # UMAP by cluster
    if cluster_key in adata.obs.columns:
        try:
            sc.pl.umap(adata, color=cluster_key, legend_loc="on data", show=False)
            filename = f"umap_leiden_res{res_str}.png" if res_str else "umap_leiden.png"
            plt.savefig(plots_dir / filename, bbox_inches="tight", dpi=150)
            plt.close()
            logging.info(f"  Saved UMAP leiden resolution {resolution} plot")
        except Exception as e:
            logging.warning(f"  Could not generate UMAP leiden plot: {e}")
    
    # UMAP by cell type annotation
    if annotation_key and annotation_key in adata.obs.columns:
        try:
            sc.pl.umap(adata, color=annotation_key, show=False)
            filename = f"umap_celltype_res{res_str}.png" if res_str else "umap_celltype.png"
            plt.savefig(plots_dir / filename, bbox_inches="tight", dpi=150)
            plt.close()
            logging.info(f"  Saved UMAP cell type annotation resolution {resolution} plot")
        except Exception as e:
            logging.warning(f"  Could not generate UMAP cell type plot: {e}")


def save_marker_dotplot(
    adata: ad.AnnData,
    plots_dir: Path,
    markers: Dict[str, List[str]],
    cluster_key: str,
    resolution: float = None
) -> None:
    """
    Generate and save dotplot of marker genes by cluster.
    
    Args:
        adata: AnnData object
        plots_dir: Directory to save plots
        markers: Dictionary mapping cell type to marker genes
        cluster_key: Key for cluster assignments in adata.obs
        resolution: Clustering resolution (for filename)
    """
    res_str = str(resolution).replace(".", "p") if resolution is not None else ""
    
    if cluster_key not in adata.obs.columns:
        return
    
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
            filename = f"marker_dotplot_res{res_str}.png" if res_str else "marker_dotplot.png"
            plt.savefig(plots_dir / filename, bbox_inches="tight", dpi=150)
            plt.close()
            logging.info(f"  Saved marker dotplot resolution {resolution}")
        else:
            logging.warning(f"  No marker genes found in dataset for dotplot")
    except Exception as e:
        logging.warning(f"  Could not generate marker dotplot: {e}")


def save_de_plots(
    adata: ad.AnnData,
    plots_dir: Path,
    cluster_key: str,
    resolution: float = None
) -> None:
    """
    Generate and save differential expression plots (dotplot and heatmap).
    
    Args:
        adata: AnnData object with DE results
        plots_dir: Directory to save plots
        cluster_key: Key for cluster assignments in adata.obs
        resolution: Clustering resolution (for filename)
    """
    res_str = str(resolution).replace(".", "p") if resolution is not None else ""
    rank_key = f"rank_genes_{cluster_key}"
    
    if rank_key not in adata.uns:
        return
    
    # Differential expression dotplot
    try:
        sc.pl.rank_genes_groups_dotplot(
            adata,
            n_genes=5,  # Top 5 genes per cluster
            key=rank_key,
            groupby=cluster_key,
            show=False,
            dendrogram=False
        )
        filename = f"deg_dotplot_res{res_str}.png" if res_str else "deg_dotplot.png"
        plt.savefig(plots_dir / filename, bbox_inches="tight", dpi=150)
        plt.close()
        logging.info(f"  Saved differential expression dotplot resolution {resolution}")
    except Exception as e:
        logging.warning(f"  Could not generate differential expression dotplot: {e}")
    
    # Differential expression heatmap
    try:
        sc.pl.rank_genes_groups_heatmap(
            adata,
            n_genes=10,  # Top 10 genes per cluster
            key=rank_key,
            groupby=cluster_key,
            show=False,
            show_gene_labels=True
        )
        filename = f"deg_heatmap_res{res_str}.png" if res_str else "deg_heatmap.png"
        plt.savefig(plots_dir / filename, bbox_inches="tight", dpi=150)
        plt.close()
        logging.info(f"  Saved differential expression heatmap resolution {resolution}")
    except Exception as e:
        logging.warning(f"  Could not generate differential expression heatmap: {e}")


def create_enrichment_heatmap(
    adata: ad.AnnData,
    plots_dir: Path,
    cluster_key: str,
    resolution: float = None
) -> None:
    """
    Create a heatmap showing MLM enrichment scores per cluster for each cell type.
    
    Args:
        adata: AnnData object with MLM scores in obsm
        plots_dir: Directory to save plot
        cluster_key: Key in adata.obs containing cluster assignments
        resolution: Clustering resolution (for filename)
    """
    # Check for MLM scores from marker-based annotation (score_mlm)
    # or from resource-based scoring (score_mlm_*)
    mlm_key = None
    if "score_mlm" in adata.obsm:
        mlm_key = "score_mlm"
    else:
        # Look for any score_mlm_* key
        mlm_keys = [k for k in adata.obsm.keys() if k.startswith("score_mlm_")]
        if mlm_keys:
            mlm_key = mlm_keys[0]  # Use first available
    
    if mlm_key is None:
        logging.warning("MLM scores not found in adata.obsm, skipping enrichment heatmap")
        return
    
    res_str = str(resolution).replace(".", "p") if resolution is not None else ""
    
    try:
        # Extract MLM scores
        acts = dc.pp.get_obsm(adata, mlm_key)
        
        # Get cell types (column names from MLM scores)
        cell_types = acts.var_names.tolist()
        score_matrix = acts.X
        
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
        filename = f"enrichment_scores_res{res_str}.png" if res_str else "enrichment_scores.png"
        plt.savefig(plots_dir / filename, bbox_inches="tight", dpi=150)
        plt.close()
        logging.info(f"  Saved enrichment heatmap")
        
    except Exception as e:
        logging.warning(f"  Could not create enrichment heatmap: {e}")
        import traceback
        logging.debug(traceback.format_exc())


def save_all_plots(
    adata: ad.AnnData,
    output_dir: Path,
    resolutions: List[float],
    markers: Optional[Dict[str, List[str]]] = None
) -> None:
    """
    Generate and save all plots for the analysis.
    
    Args:
        adata: AnnData object
        output_dir: Base output directory
        resolutions: List of clustering resolutions used
        markers: Optional marker gene dictionary for dotplots
    """
    logging.info("Generating plots")
    plots_dir = create_plots_directory(output_dir)
    
    # Save QC plots
    save_qc_plots(adata, plots_dir)
    
    # Save plots for each resolution
    for resolution in resolutions:
        res_str = str(resolution).replace(".", "p")
        cluster_key = f"leiden_res{res_str}"
        annotation_key = f"cell_type_res{res_str}"
        
        # UMAP plots
        save_umap_plots(
            adata,
            plots_dir,
            cluster_key,
            annotation_key if annotation_key in adata.obs.columns else None,
            resolution
        )
        
        # Marker dotplot if markers provided
        if markers:
            save_marker_dotplot(adata, plots_dir, markers, cluster_key, resolution)
        
        # Differential expression plots
        save_de_plots(adata, plots_dir, cluster_key, resolution)
        
        # Enrichment heatmap if annotation was done
        if annotation_key in adata.obs.columns:
            create_enrichment_heatmap(adata, plots_dir, cluster_key, resolution)
    
    logging.info(f"All plots saved to {plots_dir}")

