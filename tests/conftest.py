"""
Pytest fixtures for spatial_tk tests.
"""

import pytest
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import spatialdata as sd
import shutil
import gc


@pytest.fixture(scope="function", autouse=True)
def cleanup_after_test(request):
    """Automatically clean up after each test to free disk space."""
    yield
    # Force garbage collection after each test
    gc.collect()


@pytest.fixture
def test_data_dir():
    """Return path to test data directory."""
    return Path(__file__).parent / "test_data"


@pytest.fixture
def test_samples_csv(test_data_dir):
    """Return path to test samples CSV."""
    return test_data_dir / "test_samples.csv"


@pytest.fixture
def test_markers_csv(test_data_dir):
    """Return path to test markers CSV."""
    return test_data_dir / "test_markers.csv"


@pytest.fixture
def mock_adata():
    """Create a mock AnnData object for unit tests."""
    n_obs = 100
    n_vars = 300  # Increased to avoid scanpy QC issues with percent_top (needs >200)
    
    # Create random expression matrix
    np.random.seed(42)
    X = np.random.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(np.float32)
    
    # Create obs DataFrame
    obs = pd.DataFrame({
        'sample': np.random.choice(['sample1', 'sample2'], n_obs),
        'status': np.random.choice(['HIV', 'NEG'], n_obs),
        'location': np.random.choice(['Drexel', 'OSU'], n_obs)
    })
    obs.index = [f'cell_{i}' for i in range(n_obs)]
    
    # Create var DataFrame with some real gene names for testing annotation
    # Include marker genes that will be used in tests
    real_genes = ['CD3D', 'CD3E', 'MS4A1', 'CD19', 'CD68', 'CD14']
    generic_genes = [f'gene_{i}' for i in range(n_vars - len(real_genes))]
    all_genes = real_genes + generic_genes
    
    var = pd.DataFrame({
        'gene_name': all_genes
    })
    var.index = all_genes
    
    # Create AnnData object
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    return adata


@pytest.fixture
def mock_adata_with_clusters(mock_adata):
    """Create a mock AnnData object with clustering results."""
    # Add PCA
    np.random.seed(42)
    mock_adata.obsm['X_pca'] = np.random.randn(mock_adata.n_obs, 50)
    
    # Add UMAP
    mock_adata.obsm['X_umap'] = np.random.randn(mock_adata.n_obs, 2)
    
    # Add clustering - ensure we have multiple clusters by design
    # Divide cells into 5 groups deterministically
    n_obs = mock_adata.n_obs
    cluster_labels = np.array([str(i % 5) for i in range(n_obs)])
    mock_adata.obs['leiden_res0p5'] = pd.Categorical(cluster_labels)
    
    return mock_adata


@pytest.fixture
def mock_markers():
    """Return a mock markers dictionary."""
    return {
        'T cells': ['CD3D', 'CD3E'],
        'B cells': ['MS4A1', 'CD19'],
        'Macrophages': ['CD68', 'CD14']
    }


@pytest.fixture
def subsampled_zarr_path(test_data_dir):
    """Return path to first subsampled zarr file."""
    zarr_files = list(test_data_dir.glob("subsampled_*.zarr"))
    if zarr_files:
        return zarr_files[0]
    return None


@pytest.fixture
def tmp_zarr_cleanup(tmp_path):
    """
    Provide a temp directory that aggressively cleans up .zarr files.
    Use this instead of tmp_path for tests that create large files.
    """
    yield tmp_path
    # Clean up all .zarr directories immediately after test
    for zarr_dir in tmp_path.glob("*.zarr"):
        if zarr_dir.is_dir():
            shutil.rmtree(zarr_dir, ignore_errors=True)
    gc.collect()

