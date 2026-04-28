#!/usr/bin/env python3
"""
Main CLI entry point for spatial-tk.

This module provides the command-line interface for the Xenium spatial
transcriptomics processing pipeline.
"""

import argparse
import sys
import warnings

# Suppress warnings
warnings.filterwarnings('ignore', category=FutureWarning)

from spatial_tk.commands import (
    concat,
    normalize,
    cluster,
    quantitate,
    spatial_neighbors,
    assign,
    differential,
)
from spatial_tk.utils.helpers import setup_logging


def create_parser() -> argparse.ArgumentParser:
    """
    Create the main argument parser with all subcommands.
    
    Returns:
        Configured ArgumentParser
    """
    parser = argparse.ArgumentParser(
        prog='spatial-tk',
        description='Xenium Spatial Transcriptomics Processing Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Concatenate multiple samples
  spatial-tk concat --input samples.csv --output merged.zarr
  
  # Normalize (in place to save space)
  spatial-tk normalize --input merged.zarr --inplace --save-plots
  
  # Cluster with multiple resolutions
  spatial-tk cluster --input merged.zarr --inplace --leiden-resolution 0.2,0.5,1.0
  
  # Score enrichment with a custom marker list (all cells)
  spatial-tk quantitate --input clustered.zarr --inplace --markers markers.csv

  # Score only fibroblasts against a custom list, plus built-in PanglaoDB
  spatial-tk quantitate --input clustered.zarr --inplace \\
      --markers markers.csv --filter-obs "cell_type==Fibroblast" \\
      --preset-resources panglao

  # Assign cell type labels to clusters from the computed scores
  spatial-tk assign --input clustered.zarr --inplace \\
      --score-key score_mlm_custom

  # Build a Squidpy spatial neighbors graph
  spatial-tk spatial_neighbors --input clustered.zarr --inplace \\
      --spatial-key spatial --n-neighs 8 --transform cosine

  # Differential analysis between groups
  spatial-tk differential --input annotated.zarr --output-dir results/ \\
      --groupby status --compare-groups HIV,NEG
  
  # Full pipeline (separate files)
  spatial-tk concat --input samples.csv --output step1_concat.zarr
  spatial-tk normalize --input step1_concat.zarr --output step2_normalized.zarr
  spatial-tk cluster --input step2_normalized.zarr --output step3_clustered.zarr
  spatial-tk quantitate --input step3_clustered.zarr --output step4_scored.zarr \\
      --markers markers.csv
  spatial-tk assign --input step4_scored.zarr --output step5_annotated.zarr \\
      --score-key score_mlm_custom
  spatial-tk differential --input step5_annotated.zarr --output-dir results/
        """
    )
    
    # Create subparsers
    subparsers = parser.add_subparsers(
        dest='command',
        help='Available commands',
        required=True
    )
    
    # Add concat subcommand
    concat_parser = subparsers.add_parser(
        'concat',
        help='Concatenate multiple Xenium .zarr files',
        description='Join multiple Xenium spatial datasets into a single .zarr file'
    )
    concat.add_arguments(concat_parser)
    concat_parser.set_defaults(func=concat.main)
    
    # Add normalize subcommand
    normalize_parser = subparsers.add_parser(
        'normalize',
        help='Normalize and preprocess data',
        description='Perform QC, filtering, normalization, and feature selection'
    )
    normalize.add_arguments(normalize_parser)
    normalize_parser.set_defaults(func=normalize.main)
    
    # Add cluster subcommand
    cluster_parser = subparsers.add_parser(
        'cluster',
        help='Perform clustering analysis',
        description='Run PCA, compute neighbors, UMAP, and Leiden clustering'
    )
    cluster.add_arguments(cluster_parser)
    cluster_parser.set_defaults(func=cluster.main)
    
    # Add quantitate subcommand
    quantitate_parser = subparsers.add_parser(
        'quantitate',
        help='Run enrichment scoring (MLM/ULM) on a gene list or built-in resources',
        description=(
            'Run MLM or ULM enrichment scoring using a custom marker gene list, '
            'decoupler built-in resources (panglao, hallmark, collectri, dorothea, progeny), '
            'or both. Supports optional cell filtering via --filter-obs.'
        ),
    )
    quantitate.add_arguments(quantitate_parser)
    quantitate_parser.set_defaults(func=quantitate.main)

    # Add spatial_neighbors subcommand
    spatial_neighbors_parser = subparsers.add_parser(
        'spatial_neighbors',
        help='Compute spatial neighbor graph with Squidpy',
        description=(
            'Build spatial connectivities/distances with squidpy.gr.spatial_neighbors '
            'using configurable spatial key, neighbor definition, and transform.'
        ),
    )
    spatial_neighbors.add_arguments(spatial_neighbors_parser)
    spatial_neighbors_parser.set_defaults(func=spatial_neighbors.main)

    # Add assign subcommand
    assign_parser = subparsers.add_parser(
        'assign',
        help='Assign cell type labels to clusters from enrichment scores',
        description=(
            'Read an enrichment score matrix from obsm (produced by quantitate) '
            'and assign a cell type label to each cluster using a configurable strategy. '
            'Optionally runs per-cluster differential expression.'
        ),
    )
    assign.add_arguments(assign_parser)
    assign_parser.set_defaults(func=assign.main)

    # Add differential subcommand
    differential_parser = subparsers.add_parser(
        'differential',
        help='Differential expression analysis',
        description='Perform differential analysis between groups or find cluster markers'
    )
    differential.add_arguments(differential_parser)
    differential_parser.set_defaults(func=differential.main)
    
    return parser


def main():
    """Main entry point for the CLI."""
    # Setup logging
    setup_logging()
    
    # Parse arguments
    parser = create_parser()
    args = parser.parse_args()
    
    # Execute the appropriate command
    try:
        args.func(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

