"""
Unit tests for spatial clustering core utilities.
"""

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse


def _make_adata():
    X = np.ones((6, 4), dtype=float)
    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(["A", "A", "B", "B", "C", "C"]),
        },
        index=[f"cell_{i}" for i in range(6)],
    )
    var = pd.DataFrame(index=[f"g{i}" for i in range(4)])
    adata = ad.AnnData(X=X, obs=obs, var=var)

    rows = np.array([0, 0, 1, 1, 2, 3, 4, 5])
    cols = np.array([1, 2, 0, 2, 0, 4, 3, 4])
    data = np.ones(len(rows), dtype=float)
    adata.obsp["spatial_connectivities"] = sparse.csr_matrix((data, (rows, cols)), shape=(6, 6))
    return adata


def test_build_neighborhood_composition_shapes_and_categories():
    from spatial_tk.core import spatial_clustering

    adata = _make_adata()
    out = spatial_clustering.build_neighborhood_composition(
        adata,
        connectivities_key="spatial_connectivities",
        cell_type_key="cell_type",
        include_self=True,
        normalize=False,
    )

    composition = out["composition"]
    assert composition.shape == (adata.n_obs, 3)
    assert out["cell_type_categories"] == ["A", "B", "C"]
    assert out["neighbor_counts"].shape == (adata.n_obs,)


def test_include_self_changes_neighbor_counts():
    from spatial_tk.core import spatial_clustering

    adata = _make_adata()
    out_include = spatial_clustering.build_neighborhood_composition(
        adata, "spatial_connectivities", "cell_type", include_self=True, normalize=False
    )
    out_exclude = spatial_clustering.build_neighborhood_composition(
        adata, "spatial_connectivities", "cell_type", include_self=False, normalize=False
    )
    assert np.all(out_include["neighbor_counts"] >= out_exclude["neighbor_counts"])
    assert np.any(out_include["neighbor_counts"] > out_exclude["neighbor_counts"])


def test_normalize_composition_row_sums_one_for_nonzero_rows():
    from spatial_tk.core import spatial_clustering

    adata = _make_adata()
    out = spatial_clustering.build_neighborhood_composition(
        adata, "spatial_connectivities", "cell_type", include_self=True, normalize=True
    )
    row_sums = out["composition"].sum(axis=1)
    assert np.allclose(row_sums, np.ones(adata.n_obs))


def test_run_spatial_kmeans_returns_full_sweep():
    from spatial_tk.core import spatial_clustering

    rng = np.random.default_rng(7)
    composition = rng.normal(size=(40, 5))
    out = spatial_clustering.run_spatial_kmeans(
        composition=composition,
        min_clusters=2,
        max_clusters=6,
        random_state=0,
    )
    assert out["n_clusters"] == [2, 3, 4, 5, 6]
    assert len(out["silhouette_scores"]) == 5
    assert len(out["labels_by_n_clusters"]["2"]) == 40
    assert out["selection_method"] == "silhouette_max"


def test_run_spatial_kmeans_force_n_clusters_overrides_selection():
    from spatial_tk.core import spatial_clustering

    rng = np.random.default_rng(11)
    composition = rng.normal(size=(40, 6))
    out = spatial_clustering.run_spatial_kmeans(
        composition=composition,
        min_clusters=2,
        max_clusters=8,
        random_state=0,
        force_n_clusters=5,
    )
    assert out["best_n_clusters"] == 5
    assert out["force_n_clusters"] == 5
    assert out["selection_method"] == "forced"
    assert out["silhouette_best_n_clusters"] in [2, 3, 4, 5, 6, 7, 8]


def test_run_spatial_hdbscan_returns_labels_and_noise_counts():
    from spatial_tk.core import spatial_clustering

    rng = np.random.default_rng(21)
    cluster1 = rng.normal(loc=0.0, scale=0.2, size=(30, 5))
    cluster2 = rng.normal(loc=3.0, scale=0.2, size=(30, 5))
    composition = np.vstack([cluster1, cluster2])

    out = spatial_clustering.run_spatial_hdbscan(
        composition=composition,
        min_cluster_size=5,
        min_samples=2,
        cluster_selection_epsilon=0.0,
        metric="euclidean",
        allow_single_cluster=False,
    )
    assert out["mode"] == "hdbscan"
    assert len(out["labels"]) == composition.shape[0]
    assert out["n_clusters_found"] >= 1
    assert out["n_noise"] >= 0
    assert 0.0 <= out["noise_fraction"] <= 1.0


def test_run_spatial_hdbscan_silhouette_can_be_none():
    from spatial_tk.core import spatial_clustering

    # Constant vectors tend to yield one cluster or all noise.
    composition = np.zeros((20, 4), dtype=float)
    out = spatial_clustering.run_spatial_hdbscan(
        composition=composition,
        min_cluster_size=5,
        min_samples=None,
        cluster_selection_epsilon=0.0,
        metric="euclidean",
        allow_single_cluster=True,
    )
    assert "silhouette_score" in out
    assert out["silhouette_score"] is None or isinstance(out["silhouette_score"], float)

