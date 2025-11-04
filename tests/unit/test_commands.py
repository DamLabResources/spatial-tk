"""
Unit tests for CLI command modules.

These tests verify that command-line arguments are properly passed to core functions.
"""

import pytest
from unittest.mock import patch, MagicMock
from argparse import Namespace
from pathlib import Path


def test_annotate_command_passes_tmin_default():
    """Test that annotate command uses default tmin value of 2."""
    from xenium_process.commands import annotate
    
    # Create mock args with tmin default
    args = Namespace(
        input="test.zarr",
        output="output.zarr",
        inplace=False,
        markers="markers.csv",
        cluster_key=None,
        calculate_ulm=False,
        panglao_min_sensitivity=0.5,
        tmin=2,  # Default value
        save_plots=False,
        config=None  # No config file
    )
    
    # Mock all the functions that would be called
    with patch('xenium_process.commands.annotate.load_existing_spatial_data') as mock_load, \
         patch('xenium_process.commands.annotate.get_table') as mock_get_table, \
         patch('xenium_process.commands.annotate.annotation.load_marker_genes') as mock_load_markers, \
         patch('xenium_process.commands.annotate.annotation.annotate_with_markers') as mock_annotate, \
         patch('xenium_process.commands.annotate.annotation.run_differential_expression'), \
         patch('xenium_process.commands.annotate.save_spatial_data'), \
         patch('xenium_process.commands.annotate.set_table'), \
         patch('xenium_process.commands.annotate.prepare_spatial_data_for_save'), \
         patch('xenium_process.commands.annotate.Path'):
        
        # Setup mocks
        mock_sdata = MagicMock()
        mock_adata = MagicMock()
        mock_adata.obs.columns = ['leiden_res0p5']
        mock_load.return_value = mock_sdata
        mock_get_table.return_value = mock_adata
        mock_load_markers.return_value = {'T cells': ['CD3D']}
        
        # Make Path(markers).exists() return True
        with patch('xenium_process.commands.annotate.Path') as mock_path_class:
            mock_path_obj = MagicMock()
            mock_path_obj.exists.return_value = True
            mock_path_class.return_value = mock_path_obj
            
            # Run the command
            try:
                annotate.main(args)
            except SystemExit:
                pass  # Expected when paths don't exist
            
            # Verify annotate_with_markers was called with tmin=2
            assert mock_annotate.called, "annotate_with_markers should have been called"
            call_kwargs = mock_annotate.call_args[1]
            assert 'tmin' in call_kwargs, "tmin parameter should be passed"
            assert call_kwargs['tmin'] == 2, "tmin should be 2 (default)"


def test_annotate_command_passes_custom_tmin():
    """Test that annotate command respects custom tmin value."""
    from xenium_process.commands import annotate
    
    # Create mock args with custom tmin
    args = Namespace(
        input="test.zarr",
        output="output.zarr",
        inplace=False,
        markers="markers.csv",
        cluster_key=None,
        calculate_ulm=False,
        panglao_min_sensitivity=0.5,
        tmin=1,  # Custom value for small marker sets
        save_plots=False,
        config=None  # No config file
    )
    
    # Mock all the functions
    with patch('xenium_process.commands.annotate.load_existing_spatial_data') as mock_load, \
         patch('xenium_process.commands.annotate.get_table') as mock_get_table, \
         patch('xenium_process.commands.annotate.annotation.load_marker_genes') as mock_load_markers, \
         patch('xenium_process.commands.annotate.annotation.annotate_with_markers') as mock_annotate, \
         patch('xenium_process.commands.annotate.annotation.run_differential_expression'), \
         patch('xenium_process.commands.annotate.save_spatial_data'), \
         patch('xenium_process.commands.annotate.set_table'), \
         patch('xenium_process.commands.annotate.prepare_spatial_data_for_save'):
        
        # Setup mocks
        mock_sdata = MagicMock()
        mock_adata = MagicMock()
        mock_adata.obs.columns = ['leiden_res0p5']
        mock_load.return_value = mock_sdata
        mock_get_table.return_value = mock_adata
        mock_load_markers.return_value = {'T cells': ['CD3D']}
        
        # Make Path(markers).exists() return True
        with patch('xenium_process.commands.annotate.Path') as mock_path_class:
            mock_path_obj = MagicMock()
            mock_path_obj.exists.return_value = True
            mock_path_class.return_value = mock_path_obj
            
            # Run the command
            try:
                annotate.main(args)
            except SystemExit:
                pass
            
            # Verify annotate_with_markers was called with tmin=1
            assert mock_annotate.called
            call_kwargs = mock_annotate.call_args[1]
            assert call_kwargs['tmin'] == 1, "tmin should be 1 (custom value)"


def test_annotate_command_without_markers_no_tmin_error():
    """Test that annotate command works without markers (no tmin needed)."""
    from xenium_process.commands import annotate
    
    # Create mock args without markers
    args = Namespace(
        input="test.zarr",
        output="output.zarr",
        inplace=False,
        markers=None,  # No markers
        cluster_key="leiden_res0p5",
        calculate_ulm=False,
        panglao_min_sensitivity=0.5,
        save_plots=False,
        config=None  # No config file
    )
    
    # Mock all the functions
    with patch('xenium_process.commands.annotate.load_existing_spatial_data') as mock_load, \
         patch('xenium_process.commands.annotate.get_table') as mock_get_table, \
         patch('xenium_process.commands.annotate.annotation.run_differential_expression'), \
         patch('xenium_process.commands.annotate.save_spatial_data'), \
         patch('xenium_process.commands.annotate.set_table'), \
         patch('xenium_process.commands.annotate.prepare_spatial_data_for_save'), \
         patch('xenium_process.commands.annotate.Path'):
        
        mock_sdata = MagicMock()
        mock_adata = MagicMock()
        mock_adata.obs.columns = ['leiden_res0p5']
        mock_load.return_value = mock_sdata
        mock_get_table.return_value = mock_adata
        
        # Should not raise an error
        try:
            annotate.main(args)
        except SystemExit:
            pass  # Expected when paths don't exist

