#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


MARKERS = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "Prostate_luminal": ["KLK3", "ACPP", "AR", "NKX3-1"],
    "Basal": ["KRT5", "KRT14", "TP63"],
    "Tumor_stress": ["AMACR", "MKI67", "TOP2A", "EGR1"],
    "Immune": ["PTPRC", "CD3D", "CD3E", "MS4A1", "LYZ"],
    "Stromal": ["COL1A1", "COL1A2", "DCN", "LUM", "ACTA2"],
    "Endothelial": ["PECAM1", "VWF", "KDR"],
}


def score_markers(adata: sc.AnnData) -> pd.DataFrame:
    x = adata.X
    if hasattr(x, "toarray"):
        totals = np.asarray(x.sum(axis=1)).ravel()
        genes = np.asarray((x > 0).sum(axis=1)).ravel()
    else:
        totals = np.asarray(x.sum(axis=1)).ravel()
        genes = np.asarray((x > 0).sum(axis=1)).ravel()
    out = pd.DataFrame(index=adata.obs_names)
    out["total_counts"] = totals
    out["n_genes"] = genes
    upper = {g.upper(): g for g in adata.var_names}
    for name, geneset in MARKERS.items():
        present = [upper[g] for g in geneset if g in upper]
        if present:
            sub = adata[:, present].X
            out[f"{name}_score"] = np.asarray(sub.mean(axis=1)).ravel()
            out[f"{name}_genes"] = ",".join(present)
        else:
            out[f"{name}_score"] = 0.0
            out[f"{name}_genes"] = ""
    return out


def scatter_spatial(df: pd.DataFrame, x: np.ndarray, y: np.ndarray, color: str, path: Path, title: str) -> None:
    plt.figure(figsize=(7, 7))
    values = df[color].to_numpy()
    if np.issubdtype(values.dtype, np.number):
        vmax = np.nanpercentile(values, 99)
        plt.scatter(x, y, c=values, s=6, cmap="viridis", vmax=vmax)
        plt.colorbar(label=color)
    else:
        plt.scatter(x, y, s=6)
    plt.gca().invert_yaxis()
    plt.title(title)
    plt.xlabel("spatial x")
    plt.ylabel("spatial y")
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate QC, marker, and spatial overlay report for Cellpose/bin2cell output.")
    parser.add_argument("--cell-h5ad", type=Path, required=True)
    parser.add_argument("--raw-cell-h5ad", type=Path)
    parser.add_argument("--gex-sweep-json", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = args.output_dir / "figures"
    fig_dir.mkdir(exist_ok=True)

    adata = sc.read_h5ad(args.raw_cell_h5ad if args.raw_cell_h5ad and args.raw_cell_h5ad.exists() else args.cell_h5ad)
    metrics = score_markers(adata)
    for col in adata.obs.columns:
        if col in ["object_id", "bin_count", "labels_joint_source", "cellid"]:
            metrics[col] = adata.obs[col].astype(str).to_numpy()
    metrics.to_csv(args.output_dir / "cellpose_cell_qc_marker_table.csv")

    plt.figure(figsize=(8, 5))
    sns.histplot(metrics["total_counts"], bins=50)
    plt.title("Cell-level total raw counts")
    plt.tight_layout()
    plt.savefig(fig_dir / "qc_total_counts_hist.png", dpi=180)
    plt.close()

    plt.figure(figsize=(8, 5))
    sns.histplot(metrics["n_genes"], bins=50)
    plt.title("Cell-level detected genes")
    plt.tight_layout()
    plt.savefig(fig_dir / "qc_n_genes_hist.png", dpi=180)
    plt.close()

    score_cols = [c for c in metrics.columns if c.endswith("_score")]
    plt.figure(figsize=(10, max(4, 0.35 * len(score_cols))))
    sns.heatmap(metrics[score_cols].corr(), cmap="vlag", center=0)
    plt.title("Marker score correlation")
    plt.tight_layout()
    plt.savefig(fig_dir / "marker_score_correlation.png", dpi=180)
    plt.close()

    coords = np.asarray(adata.obsm["spatial"])
    scatter_spatial(metrics, coords[:, 0], coords[:, 1], "total_counts", fig_dir / "spatial_total_counts.png", "Spatial total counts")
    for col in score_cols[:]:
        scatter_spatial(metrics, coords[:, 0], coords[:, 1], col, fig_dir / f"spatial_{col}.png", f"Spatial {col}")

    gex_sweep = None
    if args.gex_sweep_json and args.gex_sweep_json.exists():
        gex_sweep = json.loads(args.gex_sweep_json.read_text())

    summary = {
        "cell_h5ad_used": str(args.raw_cell_h5ad if args.raw_cell_h5ad and args.raw_cell_h5ad.exists() else args.cell_h5ad),
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "total_counts_sum": float(metrics["total_counts"].sum()),
        "median_total_counts": float(metrics["total_counts"].median()),
        "median_n_genes": float(metrics["n_genes"].median()),
        "marker_scores": {col: float(metrics[col].mean()) for col in score_cols},
        "gex_sweep": gex_sweep,
    }
    (args.output_dir / "cellpose_qc_marker_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    sweep_rows = ""
    if gex_sweep:
        for row in gex_sweep.get("results", []):
            sweep_rows += (
                f"<tr><td>{row['name']}</td><td>{row['status']}</td><td>{row['sigma']}</td>"
                f"<td>{row['prob_thresh']}</td><td>{row['nms_thresh']}</td>"
                f"<td>{row['labels_gex_nonzero_bins']}</td><td>{row['labels_gex_unique_objects']}</td></tr>"
            )
    html = f"""<!doctype html>
<html lang=\"zh-CN\"><head><meta charset=\"utf-8\"><title>SC000895-R4 Cellpose QC Marker Report</title>
<style>body{{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;margin:28px;line-height:1.55;color:#1f2933}}img{{max-width:48%;margin:6px;border:1px solid #ddd}}table{{border-collapse:collapse}}td,th{{border:1px solid #ddd;padding:6px 8px}}th{{background:#eef3f8}}code{{background:#f3f4f6;padding:2px 4px}}</style></head>
<body><h1>SC000895-R4 Cellpose QC / Marker / Spatial Report</h1>
<p>Cells: <strong>{adata.n_obs}</strong>; genes: <strong>{adata.n_vars}</strong>; median counts: <strong>{summary['median_total_counts']:.1f}</strong>; median genes: <strong>{summary['median_n_genes']:.1f}</strong>.</p>
<h2>QC</h2><img src=\"figures/qc_total_counts_hist.png\"><img src=\"figures/qc_n_genes_hist.png\">
<h2>Marker Spatial Overlays</h2>
{''.join(f'<img src=\"figures/spatial_{c}.png\">' for c in score_cols)}
<h2>Marker Score Correlation</h2><img src=\"figures/marker_score_correlation.png\">
<h2>GEX Cellpose Sweep</h2><table><tr><th>name</th><th>status</th><th>sigma</th><th>prob</th><th>nms</th><th>nonzero bins</th><th>objects</th></tr>{sweep_rows}</table>
<h2>Files</h2><p>QC table: <code>cellpose_cell_qc_marker_table.csv</code></p><p>Summary: <code>cellpose_qc_marker_summary.json</code></p>
</body></html>"""
    (args.output_dir / "cellpose_qc_marker_report.html").write_text(html)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
