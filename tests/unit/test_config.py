"""
Unit tests for config utilities.
"""

import pytest
import argparse
import tempfile
from pathlib import Path

from xenium_process.utils.config import load_config, merge_config_with_args, convert_value


class TestLoadConfig:
    """Test load_config function."""
    
    def test_load_valid_config(self):
        """Test loading a valid TOML config file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.toml', delete=False) as f:
            f.write("""
[concat]
input = "test.csv"
output = "test.zarr"
downsample = 0.5
""")
            f.flush()
            config = load_config(f.name)
            assert 'concat' in config
            assert config['concat']['input'] == "test.csv"
            assert config['concat']['output'] == "test.zarr"
            assert config['concat']['downsample'] == 0.5
            Path(f.name).unlink()
    
    def test_load_nonexistent_config(self):
        """Test loading a nonexistent config file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            load_config("nonexistent.toml")
    
    def test_load_invalid_toml(self):
        """Test loading an invalid TOML file raises ValueError."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.toml', delete=False) as f:
            f.write("""
[concat
input = "test.csv"  # Missing closing bracket
""")
            f.flush()
            with pytest.raises(ValueError):
                load_config(f.name)
            Path(f.name).unlink()


class TestConvertValue:
    """Test convert_value function."""
    
    def test_convert_bool(self):
        """Test bool conversion."""
        assert convert_value(True, bool) == True
        assert convert_value(False, bool) == False
        assert convert_value("true", bool) == True
        assert convert_value("True", bool) == True
        assert convert_value("1", bool) == True
        assert convert_value("yes", bool) == True
        assert convert_value("false", bool) == False
        assert convert_value("0", bool) == False
    
    def test_convert_int(self):
        """Test int conversion."""
        assert convert_value(1, int) == 1
        assert convert_value(1.0, int) == 1
        assert convert_value("1", int) == 1
        assert convert_value("1.0", int) == 1
    
    def test_convert_float(self):
        """Test float conversion."""
        assert convert_value(1, float) == 1.0
        assert convert_value(1.5, float) == 1.5
        assert convert_value("1.5", float) == 1.5
        assert convert_value("1", float) == 1.0
    
    def test_convert_none(self):
        """Test None handling."""
        assert convert_value(None, int) is None
        assert convert_value(None, str) is None


class TestMergeConfigWithArgs:
    """Test merge_config_with_args function."""
    
    def test_merge_config_values(self):
        """Test merging config values into args."""
        parser = argparse.ArgumentParser()
        parser.add_argument('--input', required=True)
        parser.add_argument('--output', required=True)
        parser.add_argument('--downsample', type=float, default=1.0)
        
        config_dict = {
            'concat': {
                'input': 'config_input.csv',
                'output': 'config_output.zarr',
                'downsample': 0.5
            }
        }
        
        # Parse empty args (will use defaults)
        args = parser.parse_args(['--input', 'cli_input.csv', '--output', 'cli_output.zarr'])
        
        # Merge config - CLI args should override config
        merged = merge_config_with_args('concat', config_dict, args, parser)
        
        # CLI args take precedence
        assert merged.input == 'cli_input.csv'
        assert merged.output == 'cli_output.zarr'
        # downsample wasn't provided, so config value should be applied
        assert merged.downsample == 0.5
    
    def test_config_overrides_defaults(self):
        """Test that config values override argparse defaults."""
        parser = argparse.ArgumentParser()
        parser.add_argument('--min-genes', type=int, default=100)
        parser.add_argument('--n-top-genes', type=int, default=2000)
        
        config_dict = {
            'normalize': {
                'min_genes': 150,
                'n_top_genes': 3000
            }
        }
        
        # Parse with no args (uses defaults)
        args = parser.parse_args([])
        
        # Merge config
        merged = merge_config_with_args('normalize', config_dict, args, parser)
        
        # Config values should override defaults
        assert merged.min_genes == 150
        assert merged.n_top_genes == 3000
    
    def test_cli_overrides_config(self):
        """Test that CLI args override config values."""
        parser = argparse.ArgumentParser()
        parser.add_argument('--min-genes', type=int, default=100)
        
        config_dict = {
            'normalize': {
                'min_genes': 150
            }
        }
        
        # Parse with CLI arg
        args = parser.parse_args(['--min-genes', '200'])
        
        # Merge config
        merged = merge_config_with_args('normalize', config_dict, args, parser)
        
        # CLI arg should override config
        assert merged.min_genes == 200
    
    def test_missing_section(self):
        """Test handling of missing config section."""
        parser = argparse.ArgumentParser()
        parser.add_argument('--input', default='default.zarr')
        
        config_dict = {
            'other_command': {
                'input': 'other.zarr'
            }
        }
        
        args = parser.parse_args([])
        merged = merge_config_with_args('normalize', config_dict, args, parser)
        
        # Should return args unchanged
        assert merged.input == 'default.zarr'
    
    def test_type_conversion(self):
        """Test type conversion from config."""
        parser = argparse.ArgumentParser()
        parser.add_argument('--min-genes', type=int, default=100)
        parser.add_argument('--save-plots', action='store_true')
        
        config_dict = {
            'normalize': {
                'min_genes': '150',  # String that should convert to int
                'save_plots': 'true'  # String that should convert to bool
            }
        }
        
        args = parser.parse_args([])
        merged = merge_config_with_args('normalize', config_dict, args, parser)
        
        assert merged.min_genes == 150
        assert merged.save_plots == True
    
    def test_underscore_hyphen_mapping(self):
        """Test that underscore config keys map to hyphen CLI args."""
        parser = argparse.ArgumentParser()
        parser.add_argument('--min-genes', type=int, default=100)
        parser.add_argument('--n-top-genes', type=int, default=2000)
        
        config_dict = {
            'normalize': {
                'min_genes': 150,      # Underscore in config
                'n_top_genes': 3000     # Underscore in config
            }
        }
        
        args = parser.parse_args([])
        merged = merge_config_with_args('normalize', config_dict, args, parser)
        
        # Should map correctly (argparse stores as min_genes with underscore)
        assert merged.min_genes == 150
        assert merged.n_top_genes == 3000

