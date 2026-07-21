#!/usr/bin/env python3
"""Create a coarse per-sample Leiden view targeted to 20-30 clusters."""
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

PROJECT_HOME = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
ROOT = Path(os.environ.get("SP_OFFICIAL_RESULTS", PROJECT_HOME / "results" / "segmented_official_v1"))
OUT_ROOT = ROOT / "analysis_coarse"
SAMPLES = ["SC000895-R1", "SC000895-R2", "SC000895-R4", "SC000895-R5", "SC000895-R6", "SC000895-R7", "SC000895-R8", "SC000895-R9"]
RESOLUTIONS = [0.40, 0.25, 0.15, 0.10, 0.075, 0.05, 0.035, 0.025, 0.015, 0.01, 0.005]
TARGET_MIN, TARGET_MAX, TARGET = 20, 30, 25


def choose_resolution(graph) -> tuple[float, dict[str, int]]:
    counts = {}
    labels = {}
    for resolution in RESOLUTIONS:
        key = f"coarse_{str(resolution).replace('.', '_')}"
        sc.tl.leiden(graph, resolution=resolution, flavor="igraph", directed=False, key_added=key, random_state=0)
        labels[resolution] = graph.obs[key].astype(str).copy()
        counts[str(resolution)] = int(graph.obs[key].nunique())
    candidates = [(abs(n - TARGET), -resolution, resolution) for resolution, n in ((float(k), v) for k, v in counts.items()) if TARGET_MIN <= n <= TARGET_MAX]
    if not candidates:
        candidates = [(abs(n - TARGET), -resolution, resolution) for resolution, n in ((float(k), v) for k, v in counts.items())]
    chosen = min(candidates)[2]
    graph.obs["coarse_leiden"] = labels[chosen].values
    return chosen, counts


def main() -> None:
    index = int(os.environ.get("LSB_JOBINDEX", "1"))
    sample = SAMPLES[index - 1]
    input_path = ROOT / "analysis" / "samples" / sample / f"{sample}_official_cell_level_analysis.h5ad"
    out = OUT_ROOT / "samples" / sample
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(input_path)
    if "counts" not in adata.layers:
        raise ValueError(f"{sample}: counts layer missing")
    nonzero = adata.obs["analysis_nonzero_counts"].to_numpy(dtype=bool)
    graph = adata[nonzero].copy()
    sc.pp.highly_variable_genes(graph, n_top_genes=2000, flavor="cell_ranger", subset=True)
    n_pcs = max(2, min(30, graph.n_vars - 1, graph.n_obs - 1))
    sc.pp.pca(graph, n_comps=n_pcs, zero_center=False, random_state=0)
    sc.pp.neighbors(graph, n_neighbors=15, n_pcs=n_pcs, random_state=0)
    chosen, candidate_counts = choose_resolution(graph)
    sc.tl.umap(graph, random_state=0)
    adata.obs["cell_level_leiden_coarse"] = "unassigned_zero_counts"
    adata.obs.loc[graph.obs_names, "cell_level_leiden_coarse"] = graph.obs["coarse_leiden"].astype(str)
    adata.obs["cell_level_leiden_coarse"] = adata.obs["cell_level_leiden_coarse"].astype("category")
    adata.obsm["X_umap_coarse"] = np.full((adata.n_obs, 2), np.nan, dtype=np.float32)
    adata.obsm["X_umap_coarse"][nonzero] = graph.obsm["X_umap"]
    adata.uns["cell_level_analysis_coarse"] = {
        "target_min_clusters": TARGET_MIN, "target_max_clusters": TARGET_MAX, "target_clusters": TARGET,
        "candidate_resolutions": RESOLUTIONS, "candidate_cluster_counts": candidate_counts,
        "chosen_resolution": chosen, "chosen_cluster_count": int(graph.obs["coarse_leiden"].nunique()),
        "graph_backend": "scanpy_igraph_neighbors", "leiden_backend": "scanpy_igraph",
        "v1_preserved": True, "zero_counts_retained": True,
    }
    sc.pl.umap(graph, color="coarse_leiden", show=False, frameon=False)
    plt.savefig(figures / f"{sample}_coarse_umap_leiden.png", dpi=180, bbox_inches="tight")
    plt.close("all")
    if "spatial" in adata.obsm:
        positions = {cell: i for i, cell in enumerate(graph.obs_names)}
        keep = np.arange(graph.n_obs)
        if len(keep) > 100000:
            keep = np.random.default_rng(0).choice(keep, 100000, replace=False)
        labels = graph.obs["coarse_leiden"].astype("category").cat.codes.to_numpy()
        plt.figure(figsize=(8, 8))
        plt.scatter(graph.obsm["spatial"][keep, 0] if "spatial" in graph.obsm else adata.obsm["spatial"][np.asarray([positions[c] for c in graph.obs_names[keep]]), 0], graph.obsm["spatial"][keep, 1] if "spatial" in graph.obsm else adata.obsm["spatial"][np.asarray([positions[c] for c in graph.obs_names[keep]]), 1], c=labels[keep], s=0.15, linewidths=0, rasterized=True, cmap="tab20")
        plt.gca().invert_yaxis(); plt.axis("equal"); plt.axis("off")
        plt.title(f"{sample}: coarse Leiden spatial overview")
        plt.savefig(figures / f"{sample}_coarse_spatial_leiden.png", dpi=180, bbox_inches="tight")
        plt.close("all")
    output = out / f"{sample}_official_cell_level_analysis_coarse.h5ad"
    adata.write_h5ad(output, compression="lzf")
    summary = {
        "sample_id": sample, "status": "pass", "cells_total": int(adata.n_obs),
        "cells_nonzero_graph": int(nonzero.sum()), "zero_counts_retained": int((~nonzero).sum()),
        "chosen_resolution": chosen, "chosen_cluster_count": int(graph.obs["coarse_leiden"].nunique()),
        "target_range": [TARGET_MIN, TARGET_MAX], "candidate_cluster_counts": candidate_counts,
        "output_h5ad": str(output), "v1_preserved": True,
    }
    (out / "coarse_analysis_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
