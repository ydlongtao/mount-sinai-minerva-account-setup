#!/usr/bin/env python3
"""Per-sample first-pass analysis of official Space Ranger segmented cells.

保留全部官方细胞；零 counts 对象只在构建图时排除，不从输出 H5AD 删除。
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

PROJECT_HOME = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
MANIFEST = Path(os.environ.get("SP_PROJECT_MANIFEST", PROJECT_HOME / "config" / "sample_manifest_v2_corrected.csv"))
ROOT = Path(os.environ.get("SP_OFFICIAL_RESULTS", PROJECT_HOME / "results" / "segmented_official_v1"))
MARKERS = {
    "luminal": ["KLK3", "ACPP", "AR", "NKX3-1", "KRT8", "KRT18"],
    "basal": ["KRT5", "KRT14", "TP63", "KRT15"],
    "tumor": ["AMACR", "ERG", "MYC", "MKI67", "TOP2A"],
    "immune": ["PTPRC", "CD3D", "CD3E", "MS4A1", "CD68", "LYZ"],
    "stromal": ["COL1A1", "COL1A2", "DCN", "LUM", "ACTA2"],
    "endothelial": ["PECAM1", "VWF", "KDR", "EMCN"],
}


def manifest_row(index: int) -> dict[str, str]:
    rows = list(csv.DictReader(MANIFEST.open(newline="")))
    return rows[index - 1]


def marker_scores(adata: sc.AnnData) -> dict[str, list[str]]:
    upper = {str(g).upper(): str(g) for g in adata.var_names}
    totals = np.asarray(adata.X.sum(axis=1)).ravel().astype(float)
    scale = np.divide(1e4, totals, out=np.zeros_like(totals), where=totals > 0)
    present = {}
    for group, genes in MARKERS.items():
        found = [upper[g] for g in genes if g in upper]
        present[group] = found
        if found:
            values = np.asarray(adata[:, found].X.sum(axis=1)).ravel()
            adata.obs[f"score_{group}"] = np.log1p(values * scale).astype(np.float32)
    return present


def main() -> None:
    index = int(os.environ.get("LSB_JOBINDEX", "1"))
    row = manifest_row(index)
    sample = row["sample_id"]
    input_path = ROOT / "samples" / sample / f"{sample}_official_segmented_cells.h5ad"
    out = ROOT / "analysis" / "samples" / sample
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(input_path)
    adata.var_names_make_unique()
    if "counts" not in adata.layers:
        raise ValueError(f"{sample}: counts layer is required")
    counts = adata.layers["counts"]
    if hasattr(counts, "data") and not np.allclose(counts.data, np.rint(counts.data)):
        raise ValueError(f"{sample}: counts layer is not integer-valued")
    adata.obs["sample_id"] = sample
    adata.obs["sample_group"] = row.get("group", "unknown")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    marker_present = marker_scores(adata)
    nonzero = np.asarray(counts.sum(axis=1)).ravel() > 0
    adata.obs["analysis_nonzero_counts"] = nonzero
    graph = adata[nonzero].copy()
    if graph.n_obs < 100:
        raise ValueError(f"{sample}: too few nonzero cells for graph analysis")
    sc.pp.highly_variable_genes(graph, n_top_genes=2000, flavor="cell_ranger", subset=True)
    n_pcs = min(30, graph.n_vars - 1, graph.n_obs - 1)
    sc.pp.pca(graph, n_comps=max(2, n_pcs), zero_center=False, random_state=0)
    sc.pp.neighbors(graph, n_neighbors=15, n_pcs=max(2, n_pcs), random_state=0)
    sc.tl.leiden(graph, resolution=0.4, flavor="igraph", directed=False, key_added="leiden", random_state=0)
    sc.tl.umap(graph, random_state=0)
    adata.obs["cell_level_leiden"] = "unassigned_zero_counts"
    adata.obs.loc[graph.obs_names, "cell_level_leiden"] = graph.obs["leiden"].astype(str)
    adata.obs["cell_level_leiden"] = adata.obs["cell_level_leiden"].astype("category")
    adata.obsm["X_pca"] = np.full((adata.n_obs, graph.obsm["X_pca"].shape[1]), np.nan, dtype=np.float32)
    adata.obsm["X_umap"] = np.full((adata.n_obs, 2), np.nan, dtype=np.float32)
    adata.obsm["X_pca"][nonzero] = graph.obsm["X_pca"]
    adata.obsm["X_umap"][nonzero] = graph.obsm["X_umap"]
    adata.uns["cell_level_analysis"] = {
        "normalization": "normalize_total_target_sum_1e4_then_log1p",
        "hvg": "cell_ranger_2000", "pca_components": int(max(2, n_pcs)),
        "neighbors": 15, "leiden_resolution": 0.4,
        "leiden_backend": "scanpy_igraph", "zero_counts_retained": True,
        "zero_counts_graph_status": "excluded_from_graph_only",
    }
    sc.pl.umap(graph, color="leiden", show=False, frameon=False)
    plt.savefig(figures / f"{sample}_cell_level_umap_leiden.png", dpi=180, bbox_inches="tight")
    plt.close("all")
    if "spatial" in adata.obsm:
        keep = np.flatnonzero(nonzero)
        if len(keep) > 100000:
            keep = np.random.default_rng(0).choice(keep, 100000, replace=False)
        labels = pd.Categorical(adata.obs.iloc[keep]["cell_level_leiden"]).codes
        plt.figure(figsize=(8, 8))
        plt.scatter(adata.obsm["spatial"][keep, 0], adata.obsm["spatial"][keep, 1], c=labels, s=0.15, linewidths=0, rasterized=True, cmap="tab20")
        plt.gca().invert_yaxis(); plt.axis("equal"); plt.axis("off")
        plt.title(f"{sample}: cell-level Leiden spatial overview")
        plt.savefig(figures / f"{sample}_cell_level_spatial_leiden.png", dpi=180, bbox_inches="tight")
        plt.close("all")
    output = out / f"{sample}_official_cell_level_analysis.h5ad"
    adata.write_h5ad(output, compression="lzf")
    summary = {
        "sample_id": sample, "status": "pass", "cells_total": int(adata.n_obs),
        "cells_nonzero_graph": int(nonzero.sum()), "zero_counts_retained": int((~nonzero).sum()),
        "genes_hvg": int(graph.n_vars), "clusters": {str(k): int(v) for k, v in adata.obs.loc[graph.obs_names, "cell_level_leiden"].value_counts().sort_index().items()},
        "markers_present": marker_present, "output_h5ad": str(output), "leiden_backend": "scanpy_igraph",
    }
    (out / "cell_level_analysis_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
