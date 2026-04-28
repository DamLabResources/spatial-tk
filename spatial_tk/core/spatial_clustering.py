#!/usr/bin/env python3
"""
Spatial neighborhood composition and clustering utilities.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.cluster import KMeans, HDBSCAN
from sklearn.metrics import silhouette_score


def _as_binary_connectivity(matrix: Any, include_self: bool) -> sparse.csr_matrix:
    """Convert connectivity matrix to binary CSR neighborhood membership matrix."""
    if sparse.issparse(matrix):
        conn = matrix.tocsr().copy()
    else:
        conn = sparse.csr_matrix(matrix)

    conn.data = np.ones_like(conn.data, dtype=np.float64)
    if include_self:
        conn = conn + sparse.eye(conn.shape[0], format="csr", dtype=np.float64)
        conn.data = np.ones_like(conn.data, dtype=np.float64)
    return conn


def build_neighborhood_composition(
    adata: ad.AnnData,
    connectivities_key: str,
    cell_type_key: str,
    include_self: bool = True,
    normalize: bool = True,
) -> Dict[str, Any]:
    """
    Build per-cell neighborhood composition vectors from spatial connectivities.
    """
    if connectivities_key not in adata.obsp:
        raise KeyError(f"Connectivity key not found in adata.obsp: {connectivities_key}")
    if cell_type_key not in adata.obs.columns:
        raise KeyError(f"Cell-type key not found in adata.obs: {cell_type_key}")

    cell_types = adata.obs[cell_type_key].astype("category")
    categories = list(cell_types.cat.categories)
    codes = cell_types.cat.codes.to_numpy()
    n_cells = adata.n_obs
    n_types = len(categories)

    one_hot = np.zeros((n_cells, n_types), dtype=np.float64)
    valid_mask = codes >= 0
    one_hot[np.arange(n_cells)[valid_mask], codes[valid_mask]] = 1.0

    connectivity = _as_binary_connectivity(adata.obsp[connectivities_key], include_self=include_self)
    composition = connectivity @ one_hot

    row_sums = np.asarray(composition.sum(axis=1)).reshape(-1)
    if normalize:
        nonzero = row_sums > 0
        composition[nonzero] = composition[nonzero] / row_sums[nonzero, None]

    return {
        "composition": composition,
        "cell_type_categories": categories,
        "neighbor_counts": row_sums.astype(np.float64),
    }


def run_spatial_kmeans(
    composition: np.ndarray,
    min_clusters: int = 2,
    max_clusters: int = 20,
    random_state: int = 0,
    force_n_clusters: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Run k-means for cluster counts in range and score with silhouette.
    """
    if min_clusters < 2:
        raise ValueError("--min-clusters must be >= 2")
    if max_clusters < min_clusters:
        raise ValueError("--max-clusters must be >= --min-clusters")
    if force_n_clusters is not None and not (min_clusters <= force_n_clusters <= max_clusters):
        raise ValueError("--force-n-clusters must be within [min-clusters, max-clusters]")

    n_samples = composition.shape[0]
    tested_clusters: list[int] = []
    silhouette_scores: list[float] = []
    inertia_values: list[float] = []
    labels_by_n_clusters: Dict[str, list[int]] = {}

    for n_clusters in range(min_clusters, max_clusters + 1):
        if n_clusters >= n_samples:
            break

        model = KMeans(n_clusters=n_clusters, random_state=random_state, n_init=10)
        labels = model.fit_predict(composition)
        score = silhouette_score(composition, labels)

        tested_clusters.append(n_clusters)
        silhouette_scores.append(float(score))
        inertia_values.append(float(model.inertia_))
        labels_by_n_clusters[str(n_clusters)] = labels.astype(int).tolist()

    if not tested_clusters:
        raise ValueError("No valid cluster counts were tested. Check sample size and cluster bounds.")

    best_idx = int(np.argmax(silhouette_scores))
    silhouette_best_n_clusters = tested_clusters[best_idx]
    best_n_clusters = force_n_clusters if force_n_clusters is not None else silhouette_best_n_clusters
    best_labels = labels_by_n_clusters[str(best_n_clusters)]

    return {
        "mode": "kmeans",
        "n_clusters": tested_clusters,
        "silhouette_scores": silhouette_scores,
        "inertia": inertia_values,
        "labels_by_n_clusters": labels_by_n_clusters,
        "silhouette_best_n_clusters": int(silhouette_best_n_clusters),
        "best_n_clusters": int(best_n_clusters),
        "best_silhouette_score": float(silhouette_scores[best_idx]),
        "selection_method": "forced" if force_n_clusters is not None else "silhouette_max",
        "force_n_clusters": int(force_n_clusters) if force_n_clusters is not None else None,
        "best_labels": best_labels,
    }


def run_spatial_hdbscan(
    composition: np.ndarray,
    min_cluster_size: int = 5,
    min_samples: Optional[int] = None,
    cluster_selection_epsilon: float = 0.0,
    metric: str = "euclidean",
    allow_single_cluster: bool = False,
) -> Dict[str, Any]:
    """
    Run sklearn HDBSCAN on neighborhood composition vectors.
    """
    model = HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_epsilon=cluster_selection_epsilon,
        metric=metric,
        allow_single_cluster=allow_single_cluster,
    )
    labels = model.fit_predict(composition).astype(int)

    non_noise = labels != -1
    non_noise_labels = labels[non_noise]
    unique_non_noise = np.unique(non_noise_labels) if np.any(non_noise) else np.array([])
    n_clusters_found = int(len(unique_non_noise))
    n_noise = int(np.sum(labels == -1))
    noise_fraction = float(n_noise / len(labels)) if len(labels) else 0.0

    silhouette = None
    if n_clusters_found >= 2 and np.sum(non_noise) >= 2:
        silhouette = float(silhouette_score(composition[non_noise], non_noise_labels))

    return {
        "mode": "hdbscan",
        "best_labels": labels.tolist(),
        "labels": labels.tolist(),
        "n_clusters_found": n_clusters_found,
        "n_noise": n_noise,
        "noise_fraction": noise_fraction,
        "silhouette_score": silhouette,
        "hdbscan_params": {
            "min_cluster_size": int(min_cluster_size),
            "min_samples": None if min_samples is None else int(min_samples),
            "cluster_selection_epsilon": float(cluster_selection_epsilon),
            "metric": metric,
            "allow_single_cluster": bool(allow_single_cluster),
        },
    }


def cluster_cell_type_composition(
    composition: np.ndarray,
    best_labels: np.ndarray,
    categories: list[str],
) -> Dict[str, list[float]]:
    """Compute mean neighborhood composition per selected cluster."""
    df = pd.DataFrame(composition, columns=categories)
    clusters = pd.Series(best_labels.astype(int), name="cluster")
    grouped = df.groupby(clusters).mean()
    return {str(idx): grouped.loc[idx].astype(float).tolist() for idx in grouped.index}


def store_spatial_cluster_results(
    adata: ad.AnnData,
    output_key: str,
    results_key: str,
    params: Dict[str, Any],
    composition: np.ndarray,
    categories: list[str],
    cluster_results: Dict[str, Any],
    store_composition_in_obsm: bool = True,
) -> ad.AnnData:
    """Store selected labels in obs and detailed run outputs in uns."""
    mode = cluster_results.get("mode", "kmeans")
    best_labels = np.asarray(cluster_results["best_labels"], dtype=int)
    adata.obs[output_key] = pd.Categorical(best_labels.astype(str))

    composition_key = f"{results_key}_composition"
    if store_composition_in_obsm:
        adata.obsm[composition_key] = composition
    else:
        composition_key = None

    common_payload: Dict[str, Any] = {
        "mode": mode,
        "params": params,
        "cell_type_categories": categories,
        "composition_key": composition_key,
        "cluster_cell_type_composition": cluster_cell_type_composition(
            composition=composition,
            best_labels=best_labels,
            categories=categories,
        ),
    }
    if mode == "kmeans":
        common_payload.update(
            {
                "n_clusters": cluster_results["n_clusters"],
                "silhouette_scores": cluster_results["silhouette_scores"],
                "inertia": cluster_results["inertia"],
                "best_n_clusters": cluster_results["best_n_clusters"],
                "best_silhouette_score": cluster_results["best_silhouette_score"],
                "selection_method": cluster_results["selection_method"],
                "force_n_clusters": cluster_results["force_n_clusters"],
                "silhouette_best_n_clusters": cluster_results["silhouette_best_n_clusters"],
                "labels_by_n_clusters": cluster_results["labels_by_n_clusters"],
            }
        )
    else:
        common_payload.update(
            {
                "labels": cluster_results["labels"],
                "n_clusters_found": cluster_results["n_clusters_found"],
                "n_noise": cluster_results["n_noise"],
                "noise_fraction": cluster_results["noise_fraction"],
                "silhouette_score": cluster_results["silhouette_score"],
            }
        )

    adata.uns[results_key] = common_payload
    return adata
