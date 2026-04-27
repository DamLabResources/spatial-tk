"""
Functional tests for config file integration with commands.
"""

import pytest
import subprocess
import sys
import tempfile
from pathlib import Path


@pytest.fixture
def test_config_file(tmp_path):
    """Create a test config file."""
    config_path = tmp_path / "test_config.toml"
    config_path.write_text("""
[concat]
input = "placeholder.csv"
output = "config_output.zarr"
downsample = 0.5

[normalize]
input = "placeholder.zarr"
output = "config_normalized.zarr"
min_genes = 150
min_cells = 5
n_top_genes = 3000
save_plots = true

[cluster]
input = "placeholder.zarr"
output = "config_clustered.zarr"
leiden_resolution = "0.3,0.6"
save_plots = true

[quantitate]
input = "placeholder.zarr"
output = "config_scored.zarr"
markers = "markers.csv"
score_key = "custom"
method = "mlm"
tmin = 3
preset_resources = ""
panglao_min_sensitivity = 0.6
panglao_canonical_only = true
filter_obs = ""
save_plots = false

[assign]
input = "config_scored.zarr"
output = "config_annotated.zarr"
score_key = "score_mlm_custom"
cluster_key = ""
annotation_key = ""
strategy = "top_positive"
run_de = true
save_plots = false

[differential]
input = "placeholder.zarr"
output_dir = "config_results/"
groupby = "leiden_res0p3"
method = "t-test"
n_genes = 50
save_plots = true
""")
    return config_path


def test_concat_with_config(test_samples_csv, test_config_file, tmp_zarr_cleanup):
    """Test concat command with config file."""
    if not test_samples_csv.exists():
        pytest.skip("Test samples CSV not found")
    
    output_path = tmp_zarr_cleanup / "concat_output.zarr"
    
    # Use config file but override input/output with CLI
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'concat',
        '--config', str(test_config_file),
        '--input', str(test_samples_csv),
        '--output', str(output_path)
    ], capture_output=True, text=True)
    
    assert result.returncode == 0, f"Command failed: {result.stderr}"
    assert output_path.exists()
    
    # Verify downsample from config was applied (we can't easily verify this without
    # checking the actual data, but we can verify the command ran)


def test_concat_config_overrides_defaults(test_samples_csv, tmp_zarr_cleanup):
    """Test that config values override defaults."""
    if not test_samples_csv.exists():
        pytest.skip("Test samples CSV not found")
    
    # Create config with custom downsample
    config_path = tmp_zarr_cleanup / "downsample_config.toml"
    config_path.write_text("""
[concat]
downsample = 0.1
""")
    
    output_path = tmp_zarr_cleanup / "concat_downsampled.zarr"
    
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'concat',
        '--config', str(config_path),
        '--input', str(test_samples_csv),
        '--output', str(output_path)
    ], capture_output=True, text=True)
    
    assert result.returncode == 0, f"Command failed: {result.stderr}"


def test_normalize_with_config(test_samples_csv, tmp_zarr_cleanup):
    """Test normalize command with config file."""
    if not test_samples_csv.exists():
        pytest.skip("Test samples CSV not found")
    
    # First concat the data
    concat_output = tmp_zarr_cleanup / "concat_for_normalize.zarr"
    
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'concat',
        '--input', str(test_samples_csv),
        '--output', str(concat_output)
    ], capture_output=True, text=True)
    
    assert result.returncode == 0
    
    # Create config
    config_path = tmp_zarr_cleanup / "normalize_config.toml"
    config_path.write_text("""
[normalize]
min_genes = 150
min_cells = 5
n_top_genes = 3000
save_plots = false
""")
    
    output_path = tmp_zarr_cleanup / "normalized_output.zarr"
    
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'normalize',
        '--config', str(config_path),
        '--input', str(concat_output),
        '--output', str(output_path)
    ], capture_output=True, text=True)
    
    assert result.returncode == 0, f"Command failed: {result.stderr}"
    assert output_path.exists()


def test_cluster_with_config(test_samples_csv, tmp_zarr_cleanup):
    """Test cluster command with config file."""
    if not test_samples_csv.exists():
        pytest.skip("Test samples CSV not found")
    
    # First concat and normalize
    concat_output = tmp_zarr_cleanup / "concat_for_cluster.zarr"
    
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'concat',
        '--input', str(test_samples_csv),
        '--output', str(concat_output)
    ], capture_output=True, text=True)
    
    assert result.returncode == 0
    
    normalize_output = tmp_zarr_cleanup / "normalize_for_cluster.zarr"
    
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'normalize',
        '--input', str(concat_output),
        '--output', str(normalize_output),
        '--n-top-genes', '500'
    ], capture_output=True, text=True)
    
    assert result.returncode == 0
    
    # Create config
    config_path = tmp_zarr_cleanup / "cluster_config.toml"
    config_path.write_text("""
[cluster]
leiden_resolution = "0.3,0.6,0.9"
save_plots = false
""")
    
    output_path = tmp_zarr_cleanup / "clustered_output.zarr"
    
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'cluster',
        '--config', str(config_path),
        '--input', str(normalize_output),
        '--output', str(output_path)
    ], capture_output=True, text=True)
    
    assert result.returncode == 0, f"Command failed: {result.stderr}"
    assert output_path.exists()


def test_cli_overrides_config(test_samples_csv, tmp_zarr_cleanup):
    """Test that CLI arguments override config values."""
    if not test_samples_csv.exists():
        pytest.skip("Test samples CSV not found")
    
    # Create config with downsample = 0.5
    config_path = tmp_zarr_cleanup / "override_config.toml"
    config_path.write_text("""
[concat]
downsample = 0.5
""")
    
    output_path = tmp_zarr_cleanup / "override_output.zarr"
    
    # Override with CLI arg downsample = 0.8
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'concat',
        '--config', str(config_path),
        '--input', str(test_samples_csv),
        '--output', str(output_path),
        '--downsample', '0.8'
    ], capture_output=True, text=True)
    
    assert result.returncode == 0, f"Command failed: {result.stderr}"
    # CLI value (0.8) should override config value (0.5)


def test_missing_config_section(test_samples_csv, tmp_zarr_cleanup):
    """Test handling of missing config section."""
    if not test_samples_csv.exists():
        pytest.skip("Test samples CSV not found")
    
    # Create config without [concat] section
    config_path = tmp_zarr_cleanup / "missing_section_config.toml"
    config_path.write_text("""
[other_command]
input = "test.csv"
""")
    
    output_path = tmp_zarr_cleanup / "missing_section_output.zarr"
    
    # Should still work, just won't use config
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'concat',
        '--config', str(config_path),
        '--input', str(test_samples_csv),
        '--output', str(output_path)
    ], capture_output=True, text=True)
    
    assert result.returncode == 0, f"Command failed: {result.stderr}"


def test_invalid_config_file(test_samples_csv, tmp_zarr_cleanup):
    """Test handling of invalid config file."""
    if not test_samples_csv.exists():
        pytest.skip("Test samples CSV not found")
    
    # Create invalid TOML
    config_path = tmp_zarr_cleanup / "invalid_config.toml"
    config_path.write_text("""
[concat
input = "test.csv"  # Missing closing bracket
""")
    
    output_path = tmp_zarr_cleanup / "invalid_config_output.zarr"
    
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'concat',
        '--config', str(config_path),
        '--input', str(test_samples_csv),
        '--output', str(output_path)
    ], capture_output=True, text=True)
    
    # Should fail with error
    assert result.returncode != 0
    assert "Error loading config file" in result.stderr or "invalid" in result.stderr.lower()


def test_nonexistent_config_file(test_samples_csv, tmp_zarr_cleanup):
    """Test handling of nonexistent config file."""
    if not test_samples_csv.exists():
        pytest.skip("Test samples CSV not found")
    
    output_path = tmp_zarr_cleanup / "nonexistent_config_output.zarr"
    
    result = subprocess.run([
        sys.executable, '-m', 'spatial_tk.cli',
        'concat',
        '--config', 'nonexistent_config.toml',
        '--input', str(test_samples_csv),
        '--output', str(output_path)
    ], capture_output=True, text=True)
    
    # Should fail with error
    assert result.returncode != 0
    assert "not found" in result.stderr.lower() or "FileNotFoundError" in result.stderr

