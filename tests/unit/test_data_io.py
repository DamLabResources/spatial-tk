"""
Unit tests for data_io module.
"""

import pytest
import pandas as pd
import tempfile
from pathlib import Path
from spatial_tk.core import data_io


def test_load_sample_metadata(test_samples_csv):
    """Test loading sample metadata from CSV."""
    df = data_io.load_sample_metadata(str(test_samples_csv))
    
    # Check that required columns exist
    assert 'sample' in df.columns
    assert 'path' in df.columns
    
    # Check that data was loaded
    assert len(df) > 0


def test_load_sample_metadata_missing_columns():
    """Test that missing columns raise an error."""
    # Create a temporary CSV without required columns
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("wrong_column\n")
        f.write("value\n")
        temp_path = f.name
    
    try:
        with pytest.raises(ValueError, match="missing required columns"):
            data_io.load_sample_metadata(temp_path)
    finally:
        Path(temp_path).unlink()


def test_load_existing_spatial_data(subsampled_zarr_path):
    """Test loading existing spatial data."""
    if subsampled_zarr_path is None:
        pytest.skip("No subsampled zarr file available")
    
    sdata = data_io.load_existing_spatial_data(subsampled_zarr_path)
    
    # Check that spatial data was loaded
    assert sdata is not None
    
    # Check that table exists
    from spatial_tk.utils.helpers import get_table
    table = get_table(sdata)
    assert table is not None
    assert table.n_obs > 0
    assert table.n_vars > 0

