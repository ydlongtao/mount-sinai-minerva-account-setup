#!/usr/bin/env python3
"""中文说明：生成 cell-level QC 和前列腺癌 marker scoring 报告。

当前版本只生成统计和图，不自动删除低质量细胞，避免在人工复核前
改变 Cellpose/bin2cell 的原始结果。
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from batch_utils import marker_groups, sample_by_index, sample_by_id, sample_root
from sp_v2_utils import write_h5ad_atomic


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate cell-level QC and marker report for one sample.")
    parser.add_argument("--sample-index", type=int, default=int(os.environ.get("LSB_JOBINDEX", "1")))
    parser.add_argument("--sample-id")
    return parser.parse_args()


def sparse_sum(x, axis: int) -> np.ndarray:
    return np.asarray(x.sum(axis=axis)).ravel()


def score_raw_markers(adata) -> pd.DataFrame:
    x = adata.X
    totals = sparse_sum(x, 1)
    genes = np.asarray((x > 0).sum(axis=1)).ravel()
    out = pd.DataFrame(index=adata.obs_names)
    out["total_counts"] = totals
    out["n_genes"] = genes
    upper = {g.upper(): g for g in adata.var_names}
    for name, geneset in marker_groups().items():
        present = [upper[g.upper()] for g in geneset if g.upper() in upper]
        col = f"{name}_score"
        if present:
            sub = adata[:, present].X
            out[col] = np.asarray(sub.mean(axis=1)).ravel()
            out[f"{name}_genes"] = ",".join(present)
        else:
            out[col] = 0.0
            out[f"{name}_genes"] = ""
    return out


def scatter_spatial(df: pd.DataFrame, coords: np.ndarray, color: str, path: Path, title: str) -> None:
    values = df[color].to_numpy()
    plt.figure(figsize=(7, 7))
    if np.issubdtype(values.dtype, np.number):
        vmax = np.nanpercentile(values, 99) if values.size else None
        plt.scatter(coords[:, 0], coords[:, 1], c=values, s=6, cmap="viridis", vmax=vmax)
        plt.colorbar(label=color)
    else:
        plt.scatter(coords[:, 0], coords[:, 1], s=6)
    plt.gca().invert_yaxis()
    plt.title(title)
    plt.xlabel("spatial x")
    plt.ylabel("spatial y")
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def main() -> None:
    args = parse_args()
    row = sample_by_index(args.sample_index) if args.sample_id is None else sample_by_id(args.sample_id)
    sample_id = row["sample_id"]
    cellpose_dir = sample_root(sample_id) / "cellpose"
    out_dir = sample_root(sample_id) / "cell_qc_marker"
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    raw_cell_h5ad = cellpose_dir / f"{sample_id}_cell_level_cellpose_raw_counts.h5ad"
    if not raw_cell_h5ad.exists():
        raise FileNotFoundError(f"Run 23_bin2cell_raw_counts_each_sample.py first: {raw_cell_h5ad}")
    adata = sc.read_h5ad(raw_cell_h5ad)
    metrics = score_raw_markers(adata)
    for col in ["object_id", "bin_count", "labels_joint_source", "cellid"]:
        if col in adata.obs:
            metrics[col] = adata.obs[col].astype(str).to_numpy()
    metrics.to_csv(out_dir / f"{sample_id}_cell_qc_marker_table.csv")

    adata.obs["total_counts"] = metrics["total_counts"].to_numpy()
    adata.obs["n_genes"] = metrics["n_genes"].to_numpy()
    for col in metrics.columns:
        if col.endswith("_score"):
            adata.obs[col] = metrics[col].to_numpy()
    write_h5ad_atomic(adata, out_dir / f"{sample_id}_cell_level_qc_marker.h5ad")

    plt.figure(figsize=(8, 5))
    sns.histplot(metrics["total_counts"], bins=50)
    plt.title(f"{sample_id} cell-level total raw counts")
    plt.tight_layout()
    plt.savefig(fig_dir / "qc_total_counts_hist.png", dpi=180)
    plt.close()

    plt.figure(figsize=(8, 5))
    sns.histplot(metrics["n_genes"], bins=50)
    plt.title(f"{sample_id} cell-level detected genes")
    plt.tight_layout()
    plt.savefig(fig_dir / "qc_n_genes_hist.png", dpi=180)
    plt.close()

    score_cols = [c for c in metrics.columns if c.endswith("_score")]
    if score_cols:
        plt.figure(figsize=(10, max(4, 0.35 * len(score_cols))))
        sns.heatmap(metrics[score_cols].corr(), cmap="vlag", center=0)
        plt.title(f"{sample_id} marker score correlation")
        plt.tight_layout()
        plt.savefig(fig_dir / "marker_score_correlation.png", dpi=180)
        plt.close()
    coords = np.asarray(adata.obsm["spatial"])
    scatter_spatial(metrics, coords, "total_counts", fig_dir / "spatial_total_counts.png", f"{sample_id} spatial total counts")
    for col in score_cols:
        scatter_spatial(metrics, coords, col, fig_dir / f"spatial_{col}.png", f"{sample_id} {col}")

    summary = {
        "sample_id": sample_id,
        "status": "pass",
        "cell_h5ad": str(raw_cell_h5ad),
        "qc_marker_h5ad": str(out_dir / f"{sample_id}_cell_level_qc_marker.h5ad"),
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "total_counts_sum": float(metrics["total_counts"].sum()),
        "median_total_counts": float(metrics["total_counts"].median()),
        "median_n_genes": float(metrics["n_genes"].median()),
        "marker_scores": {col: float(metrics[col].mean()) for col in score_cols},
    }
    (out_dir / "cell_qc_marker_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    imgs = "\n".join(f'<img src="figures/{p.name}">' for p in sorted(fig_dir.glob("*.png")))
    html = f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><title>{sample_id} Cell-Level QC Marker Report</title>
<style>body{{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;margin:28px;line-height:1.55;color:#1f2933}}img{{max-width:48%;margin:6px;border:1px solid #ddd}}code{{background:#f3f4f6;padding:2px 4px}}</style></head>
<body><h1>{sample_id} Cell-Level QC / Marker Report</h1>
<p>Cells: <strong>{adata.n_obs}</strong>; genes: <strong>{adata.n_vars}</strong>; median counts: <strong>{summary['median_total_counts']:.1f}</strong>; median genes: <strong>{summary['median_n_genes']:.1f}</strong>.</p>
<p>QC table: <code>{sample_id}_cell_qc_marker_table.csv</code></p>
{imgs}
</body></html>"""
    (out_dir / "cell_qc_marker_report.html").write_text(html)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
