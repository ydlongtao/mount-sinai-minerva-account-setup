#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


QUANTILES = [0.0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1.0]


def load_marker_genes(path: Path) -> list[str]:
    with path.open(newline="") as handle:
        return sorted({row["gene"] for row in csv.DictReader(handle)})


def array_summary(values) -> dict[str, object]:
    arr = np.asarray(values)
    finite = np.isfinite(arr)
    return {
        "shape": list(arr.shape),
        "finite_fraction": float(finite.mean()),
        "min": np.nanmin(arr, axis=0).tolist(),
        "max": np.nanmax(arr, axis=0).tolist(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--pilot-report", type=Path, required=True)
    parser.add_argument("--markers", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pilot_report = json.loads(args.pilot_report.read_text())
    marker_genes = load_marker_genes(args.markers)

    adata = ad.read_h5ad(args.h5ad, backed="r")
    cluster_sizes = adata.obs["leiden"].astype(str).value_counts().rename_axis("leiden").reset_index(name="n_bins")
    cluster_sizes["fraction"] = cluster_sizes["n_bins"] / adata.n_obs
    cluster_sizes.to_csv(args.output_dir / "cluster_sizes.csv", index=False)

    qc_columns = [column for column in ["total_counts", "n_genes_by_counts", "pct_counts_mt"] if column in adata.obs]
    qc_quantiles = adata.obs[qc_columns].quantile(QUANTILES)
    qc_quantiles.index.name = "quantile"
    qc_quantiles.to_csv(args.output_dir / "qc_quantiles.csv")

    present_markers = [gene for gene in marker_genes if gene in adata.var_names]
    metrics = {
        "sample_id": pilot_report.get("sample_id"),
        "status": pilot_report.get("status"),
        "shape": [adata.n_obs, adata.n_vars],
        "raw_counts_preserved": "counts" in adata.layers or adata.raw is not None,
        "layers": list(adata.layers.keys()),
        "obsm": list(adata.obsm.keys()),
        "obsp": list(adata.obsp.keys()),
        "spatial_coordinate_source": adata.uns.get("spatial_coordinate_source"),
        "leiden_backend": pilot_report.get("datasets", {}).get("008um", {}).get("leiden_backend"),
        "n_clusters": int(len(cluster_sizes)),
        "cluster_size": {
            "min": int(cluster_sizes["n_bins"].min()),
            "median": float(cluster_sizes["n_bins"].median()),
            "max": int(cluster_sizes["n_bins"].max()),
            "clusters_lt_10_bins": int((cluster_sizes["n_bins"] < 10).sum()),
            "clusters_lt_20_bins": int((cluster_sizes["n_bins"] < 20).sum()),
            "clusters_lt_50_bins": int((cluster_sizes["n_bins"] < 50).sum()),
            "clusters_lt_100_bins": int((cluster_sizes["n_bins"] < 100).sum()),
            "fraction_in_largest_cluster": float(cluster_sizes.iloc[0]["fraction"]),
            "fraction_in_top_2_clusters": float(cluster_sizes.iloc[:2]["fraction"].sum()),
            "fraction_in_top_10_clusters": float(cluster_sizes.iloc[:10]["fraction"].sum()),
        },
        "qc": {
            "zero_total_count_bins": int((adata.obs["total_counts"] <= 0).sum()),
            "pct_mt_ge_20_bins": int((adata.obs["pct_counts_mt"] >= 20).sum()),
            "pct_mt_ge_50_bins": int((adata.obs["pct_counts_mt"] >= 50).sum()),
            "pct_mt_eq_100_bins": int((adata.obs["pct_counts_mt"] >= 99.999).sum()),
        },
        "marker_panel": {
            "configured": len(marker_genes),
            "present_in_hvg_object": len(present_markers),
            "present": present_markers,
            "missing": [gene for gene in marker_genes if gene not in adata.var_names],
        },
        "stages": pilot_report.get("stages", {}),
    }

    if "X_umap" in adata.obsm:
        metrics["umap"] = array_summary(adata.obsm["X_umap"])
    if "spatial" in adata.obsm:
        metrics["spatial"] = array_summary(adata.obsm["spatial"])

    (args.output_dir / "review_metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    adata.file.close()


if __name__ == "__main__":
    main()
