#!/usr/bin/env python3
"""Repair coarse Leiden fragmentation for samples above 30 clusters."""
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
SAMPLES = ["SC000895-R2", "SC000895-R4", "SC000895-R6", "SC000895-R9"]
NEIGHBORS = [30, 50, 80]
RESOLUTIONS = [0.4, 0.25, 0.15, 0.1, 0.075, 0.05, 0.035, 0.025, 0.015, 0.01, 0.005, 0.002, 0.001]
TARGET_MIN, TARGET_MAX, TARGET = 20, 30, 25


def main() -> None:
    sample = SAMPLES[int(os.environ.get("LSB_JOBINDEX", "1")) - 1]
    input_path = ROOT / "analysis" / "samples" / sample / f"{sample}_official_cell_level_analysis.h5ad"
    out = ROOT / "analysis_coarse_repaired" / "samples" / sample
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(input_path)
    nonzero = adata.obs["analysis_nonzero_counts"].to_numpy(dtype=bool)
    graph_base = adata[nonzero].copy()
    sc.pp.highly_variable_genes(graph_base, n_top_genes=2000, flavor="cell_ranger", subset=True)
    n_pcs = max(2, min(30, graph_base.n_vars - 1, graph_base.n_obs - 1))
    sc.pp.pca(graph_base, n_comps=n_pcs, zero_center=False, random_state=0)
    candidates = []
    saved_labels = {}
    for neighbors in NEIGHBORS:
        graph = graph_base.copy()
        sc.pp.neighbors(graph, n_neighbors=neighbors, n_pcs=n_pcs, random_state=0)
        for resolution in RESOLUTIONS:
            key = f"repair_{neighbors}_{str(resolution).replace('.', '_')}"
            sc.tl.leiden(graph, resolution=resolution, flavor="igraph", directed=False, key_added=key, random_state=0)
            count = int(graph.obs[key].nunique())
            candidates.append({"n_neighbors": neighbors, "resolution": resolution, "clusters": count})
            saved_labels[(neighbors, resolution)] = graph.obs[key].astype(str).copy()
    eligible = [x for x in candidates if TARGET_MIN <= x["clusters"] <= TARGET_MAX]
    pool = eligible or candidates
    chosen = min(pool, key=lambda x: (abs(x["clusters"] - TARGET), x["n_neighbors"], -x["resolution"]))
    graph = graph_base.copy()
    sc.pp.neighbors(graph, n_neighbors=chosen["n_neighbors"], n_pcs=n_pcs, random_state=0)
    graph.obs["coarse_leiden_repaired"] = saved_labels[(chosen["n_neighbors"], chosen["resolution"])].values
    sc.tl.umap(graph, random_state=0)
    adata.obs["cell_level_leiden_coarse_repaired"] = "unassigned_zero_counts"
    adata.obs.loc[graph.obs_names, "cell_level_leiden_coarse_repaired"] = graph.obs["coarse_leiden_repaired"].astype(str)
    adata.obs["cell_level_leiden_coarse_repaired"] = adata.obs["cell_level_leiden_coarse_repaired"].astype("category")
    adata.obsm["X_umap_coarse_repaired"] = np.full((adata.n_obs, 2), np.nan, dtype=np.float32)
    adata.obsm["X_umap_coarse_repaired"][nonzero] = graph.obsm["X_umap"]
    adata.uns["cell_level_analysis_coarse_repaired"] = {"target_range": [TARGET_MIN, TARGET_MAX], "target_clusters": TARGET, "neighbor_candidates": NEIGHBORS, "resolution_candidates": RESOLUTIONS, "chosen_n_neighbors": chosen["n_neighbors"], "chosen_resolution": chosen["resolution"], "chosen_cluster_count": chosen["clusters"], "candidate_results_json": json.dumps(candidates, sort_keys=True), "leiden_backend": "scanpy_igraph", "v1_and_coarse_v2_preserved": True, "zero_counts_retained": True}
    sc.pl.umap(graph, color="coarse_leiden_repaired", show=False, frameon=False)
    plt.savefig(figures / f"{sample}_repaired_umap_leiden.png", dpi=180, bbox_inches="tight")
    plt.close("all")
    if "spatial" in graph.obsm:
        keep = np.arange(graph.n_obs)
        if len(keep) > 100000:
            keep = np.random.default_rng(0).choice(keep, 100000, replace=False)
        labels = graph.obs["coarse_leiden_repaired"].astype("category").cat.codes.to_numpy()
        plt.figure(figsize=(8, 8))
        plt.scatter(graph.obsm["spatial"][keep, 0], graph.obsm["spatial"][keep, 1], c=labels[keep], s=0.15, linewidths=0, rasterized=True, cmap="tab20")
        plt.gca().invert_yaxis(); plt.axis("equal"); plt.axis("off")
        plt.title(f"{sample}: repaired coarse Leiden spatial overview")
        plt.savefig(figures / f"{sample}_repaired_spatial_leiden.png", dpi=180, bbox_inches="tight")
        plt.close("all")
    output = out / f"{sample}_official_cell_level_analysis_coarse_repaired.h5ad"
    adata.write_h5ad(output, compression="lzf")
    summary = {"sample_id": sample, "status": "pass", "cells_total": int(adata.n_obs), "chosen_n_neighbors": chosen["n_neighbors"], "chosen_resolution": chosen["resolution"], "chosen_cluster_count": chosen["clusters"], "candidate_results": candidates, "output_h5ad": str(output)}
    (out / "coarse_repaired_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
