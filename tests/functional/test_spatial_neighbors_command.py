"""
Functional tests for spatial_neighbors command.
"""

import subprocess
import sys

import pytest


@pytest.mark.slow
def test_spatial_neighbors_command_persists_obsp_fields(test_samples_csv, tmp_zarr_cleanup):
    """After running spatial_neighbors, obsp graph matrices should persist after reload."""
    if not test_samples_csv.exists():
        pytest.skip("Test samples CSV not found")

    concat_output = tmp_zarr_cleanup / "concat_neighbors.zarr"

    # Prepare input by running the minimum upstream pipeline steps.
    result = subprocess.run([
        sys.executable, "-m", "spatial_tk.cli",
        "concat",
        "--input", str(test_samples_csv),
        "--output", str(concat_output),
    ], capture_output=True, text=True)
    assert result.returncode == 0, f"Concat failed: {result.stderr}"

    result = subprocess.run([
        sys.executable, "-m", "spatial_tk.cli",
        "normalize",
        "--input", str(concat_output),
        "--inplace",
        "--n-top-genes", "500",
    ], capture_output=True, text=True)
    assert result.returncode == 0, f"Normalize failed: {result.stderr}"

    result = subprocess.run([
        sys.executable, "-m", "spatial_tk.cli",
        "cluster",
        "--input", str(concat_output),
        "--inplace",
        "--leiden-resolution", "0.5",
    ], capture_output=True, text=True)
    assert result.returncode == 0, f"Cluster failed: {result.stderr}"

    # Command under test.
    result = subprocess.run([
        sys.executable, "-m", "spatial_tk.cli",
        "spatial_neighbors",
        "--input", str(concat_output),
        "--inplace",
        "--spatial-key", "spatial",
        "--n-neighs", "6",
        "--key-added", "spatial",
    ], capture_output=True, text=True)
    assert result.returncode == 0, f"spatial_neighbors failed: {result.stderr}"

    # Reload and verify persisted graph outputs.
    import spatialdata as sd
    from spatial_tk.utils.helpers import get_table

    sdata = sd.read_zarr(concat_output)
    table = get_table(sdata)
    assert table is not None
    assert "spatial_connectivities" in table.obsp
    assert "spatial_distances" in table.obsp

    connectivities = table.obsp["spatial_connectivities"]
    distances = table.obsp["spatial_distances"]
    assert connectivities.shape == (table.n_obs, table.n_obs)
    assert distances.shape == (table.n_obs, table.n_obs)
