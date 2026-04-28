"""
Functional tests for spatial_cluster command.
"""

import subprocess
import sys

import pandas as pd
import pytest


@pytest.mark.slow
def test_spatial_cluster_single_sample_persists_results(test_samples_csv, test_markers_csv, tmp_zarr_cleanup):
    """Run spatial_cluster on a one-sample pipeline for faster functional coverage."""
    if not test_samples_csv.exists():
        pytest.skip("Test samples CSV not found")
    if not test_markers_csv.exists():
        pytest.skip("Test markers CSV not found")

    # Use one sample for faster runtime.
    one_sample_csv = tmp_zarr_cleanup / "one_sample.csv"
    sample_df = pd.read_csv(test_samples_csv).iloc[[0]]
    sample_df.to_csv(one_sample_csv, index=False)

    work_zarr = tmp_zarr_cleanup / "single_sample_work.zarr"

    # concat (single sample)
    result = subprocess.run([
        sys.executable, "-m", "spatial_tk.cli",
        "concat",
        "--input", str(one_sample_csv),
        "--output", str(work_zarr),
    ], capture_output=True, text=True)
    assert result.returncode == 0, f"Concat failed: {result.stderr}"

    # normalize
    result = subprocess.run([
        sys.executable, "-m", "spatial_tk.cli",
        "normalize",
        "--input", str(work_zarr),
        "--inplace",
        "--n-top-genes", "500",
    ], capture_output=True, text=True)
    assert result.returncode == 0, f"Normalize failed: {result.stderr}"

    # expression clustering
    result = subprocess.run([
        sys.executable, "-m", "spatial_tk.cli",
        "cluster",
        "--input", str(work_zarr),
        "--inplace",
        "--leiden-resolution", "0.5",
    ], capture_output=True, text=True)
    assert result.returncode == 0, f"Cluster failed: {result.stderr}"

    # quantitate and assign so cell_type_res0p5 exists
    result = subprocess.run([
        sys.executable, "-m", "spatial_tk.cli",
        "quantitate",
        "--input", str(work_zarr),
        "--inplace",
        "--markers", str(test_markers_csv),
        "--tmin", "1",
    ], capture_output=True, text=True)
    assert result.returncode == 0, f"Quantitate failed: {result.stderr}"

    result = subprocess.run([
        sys.executable, "-m", "spatial_tk.cli",
        "assign",
        "--input", str(work_zarr),
        "--inplace",
        "--score-key", "score_mlm_custom",
        "--cluster-key", "leiden_res0p5",
    ], capture_output=True, text=True)
    assert result.returncode == 0, f"Assign failed: {result.stderr}"

    # spatial neighbors
    result = subprocess.run([
        sys.executable, "-m", "spatial_tk.cli",
        "spatial_neighbors",
        "--input", str(work_zarr),
        "--inplace",
        "--spatial-key", "spatial",
        "--n-neighs", "6",
    ], capture_output=True, text=True)
    assert result.returncode == 0, f"spatial_neighbors failed: {result.stderr}"

    # command under test
    result = subprocess.run([
        sys.executable, "-m", "spatial_tk.cli",
        "spatial_cluster",
        "--input", str(work_zarr),
        "--inplace",
        "--cell-type-key", "cell_type_res0p5",
        "--max-clusters", "10",
        "--output-key", "spatial_cluster_res",
        "--results-key", "spatial_cluster_results",
    ], capture_output=True, text=True)
    assert result.returncode == 0, f"spatial_cluster failed: {result.stderr}"

    # verify persisted outputs
    import spatialdata as sd
    from spatial_tk.utils.helpers import get_table

    sdata = sd.read_zarr(work_zarr)
    table = get_table(sdata)
    assert table is not None
    assert "spatial_cluster_res" in table.obs.columns
    assert "spatial_cluster_results" in table.uns
    result_uns = table.uns["spatial_cluster_results"]
    assert "silhouette_scores" in result_uns
    assert "labels_by_n_clusters" in result_uns
    assert "best_n_clusters" in result_uns


@pytest.mark.slow
def test_spatial_cluster_hdbscan_mode_single_sample(test_samples_csv, test_markers_csv, tmp_zarr_cleanup):
    """Run HDBSCAN mode and verify persisted HDBSCAN-specific uns fields."""
    if not test_samples_csv.exists():
        pytest.skip("Test samples CSV not found")
    if not test_markers_csv.exists():
        pytest.skip("Test markers CSV not found")

    one_sample_csv = tmp_zarr_cleanup / "one_sample_hdbscan.csv"
    sample_df = pd.read_csv(test_samples_csv).iloc[[0]]
    sample_df.to_csv(one_sample_csv, index=False)
    work_zarr = tmp_zarr_cleanup / "single_sample_hdbscan_work.zarr"

    cmds = [
        [
            sys.executable, "-m", "spatial_tk.cli", "concat",
            "--input", str(one_sample_csv), "--output", str(work_zarr),
        ],
        [
            sys.executable, "-m", "spatial_tk.cli", "normalize",
            "--input", str(work_zarr), "--inplace", "--n-top-genes", "500",
        ],
        [
            sys.executable, "-m", "spatial_tk.cli", "cluster",
            "--input", str(work_zarr), "--inplace", "--leiden-resolution", "0.5",
        ],
        [
            sys.executable, "-m", "spatial_tk.cli", "quantitate",
            "--input", str(work_zarr), "--inplace", "--markers", str(test_markers_csv), "--tmin", "1",
        ],
        [
            sys.executable, "-m", "spatial_tk.cli", "assign",
            "--input", str(work_zarr), "--inplace", "--score-key", "score_mlm_custom", "--cluster-key", "leiden_res0p5",
        ],
        [
            sys.executable, "-m", "spatial_tk.cli", "spatial_neighbors",
            "--input", str(work_zarr), "--inplace", "--spatial-key", "spatial", "--n-neighs", "6",
        ],
    ]

    for cmd in cmds:
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Step failed: {' '.join(cmd)}\n{result.stderr}"

    result = subprocess.run([
        sys.executable, "-m", "spatial_tk.cli",
        "spatial_cluster",
        "--input", str(work_zarr),
        "--inplace",
        "--cell-type-key", "cell_type_res0p5",
        "--mode", "hdbscan",
        "--hdbscan-min-cluster-size", "5",
        "--output-key", "spatial_cluster_hdbscan",
        "--results-key", "spatial_cluster_hdbscan_results",
    ], capture_output=True, text=True)
    assert result.returncode == 0, f"spatial_cluster hdbscan failed: {result.stderr}"

    import spatialdata as sd
    from spatial_tk.utils.helpers import get_table

    sdata = sd.read_zarr(work_zarr)
    table = get_table(sdata)
    assert table is not None
    assert "spatial_cluster_hdbscan" in table.obs.columns
    assert "spatial_cluster_hdbscan_results" in table.uns
    result_uns = table.uns["spatial_cluster_hdbscan_results"]
    assert result_uns["mode"] == "hdbscan"
    assert "n_clusters_found" in result_uns
    assert "n_noise" in result_uns
    assert "noise_fraction" in result_uns
