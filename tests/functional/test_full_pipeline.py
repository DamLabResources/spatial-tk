"""
Functional test for full pipeline end-to-end.
"""

import pytest
import subprocess
import sys
from pathlib import Path
import shutil


@pytest.mark.slow
def test_full_pipeline_end_to_end(test_samples_csv, test_markers_csv, tmp_zarr_cleanup):
    """
    Test the full pipeline from concat through differential analysis.
    
    This is an integration test that runs all commands sequentially.
    """
    # Skip if test data doesn't exist
    if not test_samples_csv.exists():
        pytest.skip("Test samples CSV not found")
    
    if not test_markers_csv.exists():
        pytest.skip("Test markers CSV not found")
    
    # Paths for intermediate outputs
    concat_output = tmp_zarr_cleanup / "step1_concat.zarr"
    
    # Step 1: Concat
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'concat',
        '--input', str(test_samples_csv),
        '--output', str(concat_output)
    ], capture_output=True, text=True)
    
    assert result.returncode == 0, f"Concat failed: {result.stderr}"
    assert concat_output.exists()
    
    # Step 2: Normalize (inplace)
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'normalize',
        '--input', str(concat_output),
        '--inplace',
        '--n-top-genes', '500'
    ], capture_output=True, text=True)
    
    assert result.returncode == 0, f"Normalize failed: {result.stderr}"
    
    # Step 3: Cluster (inplace)
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'cluster',
        '--input', str(concat_output),
        '--inplace',
        '--leiden-resolution', '0.5'
    ], capture_output=True, text=True)
    
    assert result.returncode == 0, f"Cluster failed: {result.stderr}"
    
    # Step 4: Quantitate – score cells against the marker gene list (inplace)
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'quantitate',
        '--input', str(concat_output),
        '--inplace',
        '--markers', str(test_markers_csv),
        '--tmin', '1',
    ], capture_output=True, text=True)

    assert result.returncode == 0, f"Quantitate failed: {result.stderr}"

    # Step 5: Assign – label clusters from the scored obsm matrix (inplace)
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'assign',
        '--input', str(concat_output),
        '--inplace',
        '--score-key', 'score_mlm_custom',
        '--cluster-key', 'leiden_res0p5',
    ], capture_output=True, text=True)

    assert result.returncode == 0, f"Assign failed: {result.stderr}"

    # Step 6: Differential analysis
    diff_output_dir = tmp_zarr_cleanup / "differential"
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'differential',
        '--input', str(concat_output),
        '--output-dir', str(diff_output_dir),
        '--groupby', 'leiden_res0p5'
    ], capture_output=True, text=True)
    
    assert result.returncode == 0, f"Differential failed: {result.stderr}"
    assert diff_output_dir.exists()
    
    # Check that output files were created
    csv_files = list(diff_output_dir.glob("*.csv"))
    assert len(csv_files) > 0, "No differential expression results found"


@pytest.mark.slow
def test_pipeline_with_group_comparison(test_samples_csv, tmp_zarr_cleanup):
    """
    Test pipeline with group comparison (Mode A differential).
    """
    if not test_samples_csv.exists():
        pytest.skip("Test samples CSV not found")
    
    # Run concat and normalize
    concat_output = tmp_zarr_cleanup / "concat.zarr"
    
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'concat',
        '--input', str(test_samples_csv),
        '--output', str(concat_output)
    ], capture_output=True, text=True)
    
    assert result.returncode == 0
    
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'normalize',
        '--input', str(concat_output),
        '--inplace',
        '--n-top-genes', '500'
    ], capture_output=True, text=True)
    
    assert result.returncode == 0
    
    # Run differential analysis comparing status groups
    diff_output_dir = tmp_zarr_cleanup / "differential_status"
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'differential',
        '--input', str(concat_output),
        '--output-dir', str(diff_output_dir),
        '--groupby', 'status',
        '--compare-groups', 'HIV,NEG'
    ], capture_output=True, text=True)
    
    assert result.returncode == 0
    assert diff_output_dir.exists()
    
    # Check for comparison output files
    comparison_files = list(diff_output_dir.glob("*_vs_*.csv"))
    assert len(comparison_files) > 0

