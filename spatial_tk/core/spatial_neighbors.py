#!/usr/bin/env python3
"""
Spatial neighbor graph utilities built on Squidpy.
"""

import logging
from typing import Optional, Tuple, Union

import anndata as ad
import squidpy as sq

RadiusType = Optional[Union[float, Tuple[float, float]]]
TransformType = Optional[str]


def parse_radius(radius: Optional[str]) -> RadiusType:
    """
    Parse radius CLI input into Squidpy-compatible value.

    Accepted formats:
    - "100" -> 100.0
    - "50,200" -> (50.0, 200.0)
    """
    if radius is None:
        return None

    value = radius.strip()
    if not value:
        return None

    if "," in value:
        parts = [p.strip() for p in value.split(",") if p.strip()]
        if len(parts) != 2:
            raise ValueError(
                "Invalid --radius format. Use a float (e.g. '100') or 'min,max' (e.g. '50,200')."
            )
        min_radius = float(parts[0])
        max_radius = float(parts[1])
        if min_radius > max_radius:
            raise ValueError("Invalid --radius interval: min must be <= max.")
        return (min_radius, max_radius)

    return float(value)


def normalize_transform(transform: Optional[str]) -> TransformType:
    """Normalize transform argument for Squidpy."""
    if transform is None:
        return None
    if transform == "none":
        return None
    return transform


def compute_spatial_neighbors(
    adata: ad.AnnData,
    spatial_key: str = "spatial",
    library_key: Optional[str] = None,
    coord_type: Optional[str] = None,
    n_neighs: int = 6,
    radius: RadiusType = None,
    transform: TransformType = None,
    key_added: str = "spatial",
) -> ad.AnnData:
    """Compute spatial neighbors in-place using squidpy.gr.spatial_neighbors."""
    logging.info("Computing spatial neighbors with squidpy")
    sq.gr.spatial_neighbors(
        adata,
        spatial_key=spatial_key,
        library_key=library_key,
        coord_type=coord_type,
        n_neighs=n_neighs,
        radius=radius,
        transform=transform,
        key_added=key_added,
        copy=False,
    )
    logging.info(
        "Spatial neighbors complete. Stored keys: obsp['%s_connectivities'], obsp['%s_distances']",
        key_added,
        key_added,
    )
    return adata
