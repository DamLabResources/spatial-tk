#!/usr/bin/env python3
"""
Clustering and dimensionality reduction functions.

This module handles PCA, neighbor graph computation, UMAP,
and Leiden clustering for spatial transcriptomics data.
"""

import logging
from typing import Optional

import anndata as ad
import scanpy as sc


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

