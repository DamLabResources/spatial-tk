"""
Unit tests for clustering module.
"""

import pytest
from spatial_tk.core import clustering, preprocessing


def test_run_pca(mock_adata):
    """Test PCA computation."""
    # Normalize first
    adata = preprocessing.normalize_and_log(mock_adata)
    adata = preprocessing.select_variable_genes(adata, n_top_genes=20)
    
    adata = clustering.run_pca(adata)
    
    # Check that PCA was computed
    assert 'X_pca' in adata.obsm
    assert adata.obsm['X_pca'].shape[0] == adata.n_obs


def test_compute_neighbors_and_umap(mock_adata):
    """Test neighbor graph and UMAP computation."""
    # Need to run PCA first
    adata = preprocessing.normalize_and_log(mock_adata)
    adata = preprocessing.select_variable_genes(adata, n_top_genes=20)
    adata = clustering.run_pca(adata)
    
    adata = clustering.compute_neighbors_and_umap(adata)
    
    # Check that neighbors and UMAP were computed
    assert 'X_umap' in adata.obsm
    assert adata.obsm['X_umap'].shape == (adata.n_obs, 2)


def test_cluster_leiden(mock_adata):
    """Test Leiden clustering."""
    # Prepare data
    adata = preprocessing.normalize_and_log(mock_adata)
    adata = preprocessing.select_variable_genes(adata, n_top_genes=20)
    adata = clustering.run_pca(adata)
    adata = clustering.compute_neighbors_and_umap(adata)
    
    # Run clustering
    resolution = 0.5
    cluster_key = "leiden_test"
    adata = clustering.cluster_leiden(adata, resolution, key_added=cluster_key)
    
    # Check that clustering was added
    assert cluster_key in adata.obs.columns
    
    # Check that we have at least 1 cluster (sometimes random data produces only 1 cluster)
    n_clusters = adata.obs[cluster_key].nunique()
    assert n_clusters >= 1
    assert n_clusters <= adata.n_obs  # Can't have more clusters than cells


def test_cluster_leiden_resume(mock_adata_with_clusters):
    """Test that resume skips already computed clustering."""
    adata = mock_adata_with_clusters
    
    # Run clustering with resume
    adata = clustering.cluster_leiden(
        adata, 
        resolution=0.5, 
        key_added='leiden_res0p5',
        resume=True
    )
    
    # Should still have the column
    assert 'leiden_res0p5' in adata.obs.columns

