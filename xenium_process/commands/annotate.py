#!/usr/bin/env python3
"""
Annotate command: Perform cell type annotation using marker genes and/or ULM scoring.
"""

import argparse
import logging
import sys
from pathlib import Path

from xenium_process.core.data_io import load_existing_spatial_data, save_spatial_data
from xenium_process.core import annotation
from xenium_process.core import plotting
from xenium_process.utils.helpers import (
    get_table, set_table, get_output_path, 
    prepare_spatial_data_for_save, parse_resolutions
)
from xenium_process.utils.config import load_config, merge_config_with_args


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """
    Add arguments for the annotate command.
    
    Args:
        parser: ArgumentParser to add arguments to
    """
    parser.add_argument(
        '--input',
        required=True,
        help='Path to input clustered .zarr file'
    )
    parser.add_argument(
        '--output',
        help='Path to output .zarr file (required unless --inplace is used)'
    )
    parser.add_argument(
        '--inplace',
        action='store_true',
        help='Modify the input file in place instead of creating a new file'
    )
    parser.add_argument(
        '--markers',
        help='Path to CSV file with marker genes (columns: cell_type, gene)'
    )
    parser.add_argument(
        '--cluster-key',
        default=None,
        help='Cluster column key to use for annotation (e.g., leiden_res0p5). If not specified, will use all leiden_res* columns found'
    )
    parser.add_argument(
        '--calculate-ulm',
        action='store_true',
        help='Pre-calculate ULM enrichment scores for pathway/TF resources'
    )
    parser.add_argument(
        '--panglao-min-sensitivity',
        type=float,
        default=0.5,
        help='Minimum sensitivity for PanglaoDB markers in ULM (default: 0.5)'
    )
    parser.add_argument(
        '--tmin',
        type=int,
        default=2,
        help='Minimum number of marker genes per cell type for ULM annotation (default: 2)'
    )
    parser.add_argument(
        '--save-plots',
        action='store_true',
        help='Generate and save annotation plots'
    )
    parser.add_argument(
        '--config',
        help='Path to TOML configuration file (optional)'
    )


def main(args: argparse.Namespace) -> None:
    """
    Execute the annotate command.
    
    Args:
        args: Parsed command-line arguments
    """
    # Load and merge config if provided
    if args.config:
        try:
            config_dict = load_config(args.config)
            # Create a temporary parser to get defaults
            temp_parser = argparse.ArgumentParser()
            add_arguments(temp_parser)
            args = merge_config_with_args('annotate', config_dict, args, temp_parser)
        except Exception as e:
            logging.error(f"Error loading config file: {e}")
            sys.exit(1)
    
    logging.info("="*60)
    logging.info("Xenium Process: Cell Type Annotation")
    logging.info("="*60)
    
    # Validate inputs
    input_path = Path(args.input)
    if not input_path.exists():
        logging.error(f"Input file not found: {input_path}")
        sys.exit(1)
    
    try:
        output_path = get_output_path(args.input, args.output, args.inplace)
    except ValueError as e:
        logging.error(str(e))
        sys.exit(1)
    
    try:
        # Load spatial data
        sdata = load_existing_spatial_data(input_path)
        adata = get_table(sdata)
        
        if adata is None:
            raise ValueError("No expression table found in spatial data")
        
        logging.info(f"Starting annotation: {adata.n_obs} cells × {adata.n_vars} genes")
        
        # Calculate ULM enrichment scores if requested (before annotation)
        if args.calculate_ulm:
            adata = annotation.calculate_ulm_scores(
                adata,
                use_panglao=True,
                panglao_min_sensitivity=args.panglao_min_sensitivity
            )
        
        # Determine which cluster keys to annotate
        if args.cluster_key:
            cluster_keys = [args.cluster_key]
        else:
            # Find all leiden_res* columns
            cluster_keys = [col for col in adata.obs.columns if col.startswith('leiden_res')]
            if not cluster_keys:
                logging.warning("No leiden clustering columns found. Run cluster command first.")
        
        # Cell type annotation if markers provided
        markers = None
        if args.markers:
            markers_path = Path(args.markers)
            if not markers_path.exists():
                logging.error(f"Markers file not found: {markers_path}")
                sys.exit(1)
            
            markers = annotation.load_marker_genes(str(markers_path))
            
            # Annotate each clustering resolution
            for cluster_key in cluster_keys:
                # Extract resolution from cluster key for annotation key naming
                res_str = cluster_key.replace("leiden_res", "")
                annotation_key = f"cell_type_res{res_str}"
                
                adata = annotation.annotate_with_markers(
                    adata,
                    markers,
                    cluster_key=cluster_key,
                    annotation_key=annotation_key,
                    tmin=args.tmin
                )
        
        # Run differential expression for each cluster key if not already done
        for cluster_key in cluster_keys:
            adata = annotation.run_differential_expression(adata, cluster_key)
        
        # Update the SpatialData table with processed AnnData
        prepare_spatial_data_for_save(adata)
        set_table(sdata, adata)
        
        # Save results
        if args.inplace:
            logging.info(f"Saving results in place: {output_path}")
        else:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            logging.info(f"Saving results to: {output_path}")
        
        save_spatial_data(sdata, output_path, overwrite=args.inplace)
        
        # Generate plots if requested
        if args.save_plots:
            plots_dir = output_path.parent / "plots"
            plots_dir.mkdir(exist_ok=True)
            
            for cluster_key in cluster_keys:
                # Extract resolution
                res_str = cluster_key.replace("leiden_res", "")
                try:
                    resolution = float(res_str.replace("p", "."))
                except ValueError:
                    resolution = None
                
                annotation_key = f"cell_type_res{res_str}"
                
                # UMAP with annotation
                if annotation_key in adata.obs.columns:
                    plotting.save_umap_plots(
                        adata, plots_dir, cluster_key, annotation_key, resolution
                    )
                
                # Marker dotplots
                if markers:
                    plotting.save_marker_dotplot(adata, plots_dir, markers, cluster_key, resolution)
                
                # Differential expression plots
                plotting.save_de_plots(adata, plots_dir, cluster_key, resolution)
                
                # Enrichment heatmap
                if annotation_key in adata.obs.columns:
                    plotting.create_enrichment_heatmap(adata, plots_dir, cluster_key, resolution)
        
        logging.info("="*60)
        logging.info(f"Annotation complete: {output_path}")
        logging.info("="*60)
        
    except Exception as e:
        logging.error(f"Annotation failed: {e}", exc_info=True)
        sys.exit(1)

