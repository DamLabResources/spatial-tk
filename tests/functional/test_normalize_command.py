"""
Functional tests for normalize command.
"""

import pytest
import subprocess
import sys
from pathlib import Path
import shutil


def test_normalize_command(test_samples_csv, tmp_zarr_cleanup):
    """Test normalize command."""
    if not test_samples_csv.exists():
        pytest.skip("Test samples CSV not found")
    
    # First concat the data
    concat_output = tmp_zarr_cleanup / "concat.zarr"
    
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'concat',
        '--input', str(test_samples_csv),
        '--output', str(concat_output)
    ], capture_output=True, text=True)
    
    assert result.returncode == 0
    
    # Now normalize
    normalize_output = tmp_zarr_cleanup / "normalized.zarr"
    
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'normalize',
        '--input', str(concat_output),
        '--output', str(normalize_output),
        '--n-top-genes', '500'
    ], capture_output=True, text=True)
    
    assert result.returncode == 0, f"Normalize failed: {result.stderr}"
    assert normalize_output.exists()
    
    # Verify normalization was applied
    import spatialdata as sd
    from spatial_tk.utils.helpers import get_table
    
    sdata = sd.read_zarr(normalize_output)
    table = get_table(sdata)
    
    assert 'highly_variable' in table.var.columns
    assert 'counts' in table.layers


def test_normalize_inplace(test_samples_csv, tmp_zarr_cleanup):
    """Test normalize command with --inplace."""
    if not test_samples_csv.exists():
        pytest.skip("Test samples CSV not found")
    
    # First concat the data
    concat_output = tmp_zarr_cleanup / "concat_inplace.zarr"
    
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'concat',
        '--input', str(test_samples_csv),
        '--output', str(concat_output)
    ], capture_output=True, text=True)
    
    assert result.returncode == 0
    
    # Normalize in place
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'normalize',
        '--input', str(concat_output),
        '--inplace',
        '--n-top-genes', '500'
    ], capture_output=True, text=True)
    
    assert result.returncode == 0
    
    # Verify that the same file was modified
    import spatialdata as sd
    from spatial_tk.utils.helpers import get_table
    
    sdata = sd.read_zarr(concat_output)
    table = get_table(sdata)
    
    assert 'highly_variable' in table.var.columns

