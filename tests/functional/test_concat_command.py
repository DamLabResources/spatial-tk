"""
Functional tests for concat command.
"""

import pytest
import subprocess
import sys
from pathlib import Path


def test_concat_command(test_samples_csv, tmp_zarr_cleanup):
    """Test concat command with test data."""
    if not test_samples_csv.exists():
        pytest.skip("Test samples CSV not found")
    
    output_path = tmp_zarr_cleanup / "concat_output.zarr"
    
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'concat',
        '--input', str(test_samples_csv),
        '--output', str(output_path)
    ], capture_output=True, text=True)
    
    # Check that command succeeded
    assert result.returncode == 0, f"Command failed: {result.stderr}"
    assert "Expected an iterable of integers" not in result.stderr
    
    # Check that output file was created
    assert output_path.exists()
    
    # Verify the output is valid
    import spatialdata as sd
    sdata = sd.read_zarr(output_path)
    assert sdata is not None
    
    from spatial_tk.utils.helpers import get_table
    table = get_table(sdata)
    assert table is not None
    assert table.n_obs > 0
    assert "sample" in table.obs
    assert table.obs["sample"].nunique() >= 1

    # Labels are optional in concat output, but if present they should be readable.
    if getattr(sdata, "labels", None):
        first_label = next(iter(sdata.labels.values()))
        first_scale = next(iter(first_label.keys()))
        image_da = first_label[first_scale].ds["image"]
        assert image_da.shape[0] > 0
        assert image_da.shape[1] > 0


def test_concat_with_downsampling(test_samples_csv, tmp_zarr_cleanup):
    """Test concat command with downsampling."""
    if not test_samples_csv.exists():
        pytest.skip("Test samples CSV not found")
    
    output_path = tmp_zarr_cleanup / "concat_downsampled.zarr"
    
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'concat',
        '--input', str(test_samples_csv),
        '--output', str(output_path),
        '--downsample', '0.5'
    ], capture_output=True, text=True)
    
    assert result.returncode == 0
    assert output_path.exists()

