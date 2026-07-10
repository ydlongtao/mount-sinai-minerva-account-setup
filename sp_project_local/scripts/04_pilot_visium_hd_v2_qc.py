#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

from sp_v2_utils import (
    PROJECT_HOME,
    RESULTS_V2,
    WORK_V2,
    apply_qc_flags,
    attach_barcode_spatial,
    candidate_qc_thresholds,
    deterministic_indices,
    file_fingerprint,
    load_config,
    read_manifest_sample,
    write_h5ad_atomic,
    write_json_atomic,
)


def make_qc_plots(adata, outdir: Path, max_points: int, seed: int) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    obs = adata.obs
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    for ax, column, color in zip(
        axes,
        ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
        ["#2563eb", "#0f766e", "#dc2626"],
    ):
        upper = obs[column].quantile(0.995)
        ax.hist(obs.loc[obs[column] <= upper, column], bins=80, color=color, alpha=0.85)
        ax.set_title(column)
        ax.set_ylabel("Bins")
    fig.tight_layout()
    fig.savefig(outdir / "qc_distributions.png", dpi=180)
    plt.close(fig)

    idx = deterministic_indices(adata.n_obs, max_points, seed)
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    hb = ax.hexbin(
        obs["total_counts"].to_numpy()[idx],
        obs["n_genes_by_counts"].to_numpy()[idx],
        gridsize=90,
        bins="log",
        mincnt=1,
        cmap="viridis",
    )
    fig.colorbar(hb, ax=ax, label="log10 bins")
    ax.set_xlabel("Total counts")
    ax.set_ylabel("Detected genes")
    ax.set_title("16 µm bin complexity")
    fig.tight_layout()
    fig.savefig(outdir / "counts_vs_genes_hexbin.png", dpi=180)
    plt.close(fig)

    coords = adata.obsm["spatial"]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for ax, column, cmap in zip(
        axes,
        ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
        ["magma", "viridis", "Reds"],
    ):
        values = obs[column].to_numpy()[idx]
        vmax = np.nanquantile(values, 0.99)
        points = ax.scatter(
            coords[idx, 0],
            coords[idx, 1],
            c=values,
            s=1,
            rasterized=True,
            cmap=cmap,
            vmin=0,
            vmax=vmax,
        )
        ax.invert_yaxis()
        ax.set_aspect("equal")
        ax.set_title(column)
        ax.axis("off")
        fig.colorbar(points, ax=ax, fraction=0.04)
    fig.tight_layout()
    fig.savefig(outdir / "spatial_qc.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=PROJECT_HOME / "config" / "pilot_v2.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    sample_id = config["sample_id"]
    outdir = RESULTS_V2 / sample_id / "qc"
    workdir = WORK_V2 / sample_id
    state_path = workdir / "run_state.json"
    outdir.mkdir(parents=True, exist_ok=True)
    workdir.mkdir(parents=True, exist_ok=True)

    sample = read_manifest_sample(PROJECT_HOME / "config" / "sample_manifest.csv", sample_id)
    matrix_path = Path(sample["matrix_016um"])
    started = time.time()
    state = {
        "sample_id": sample_id,
        "stage": "qc",
        "status": "running",
        "config": config,
        "input": file_fingerprint(matrix_path),
    }
    write_json_atomic(state_path, state)

    try:
        adata = sc.read_10x_h5(matrix_path, gex_only=True)
        adata.var_names_make_unique()
        adata.obs["sample_id"] = sample_id
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)
        attach_barcode_spatial(adata)
        thresholds = candidate_qc_thresholds(adata.obs)
        apply_qc_flags(adata, thresholds)
        adata.layers["counts"] = adata.X.copy()
        adata.uns["pilot_v2_config"] = json.dumps(config, sort_keys=True)
        adata.uns["candidate_qc_thresholds"] = thresholds

        qc_columns = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
        quantiles = adata.obs[qc_columns].quantile([0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.995, 1])
        quantiles.index.name = "quantile"
        quantiles.to_csv(outdir / "qc_quantiles.csv")
        make_qc_plots(adata, outdir, config["plot_max_points"], config["random_state"])

        checkpoint = workdir / "016um_raw_qc.h5ad"
        write_h5ad_atomic(adata, checkpoint)
        summary = {
            "sample_id": sample_id,
            "shape": [adata.n_obs, adata.n_vars],
            "qc_mode": config["qc_mode"],
            "candidate_thresholds": thresholds,
            "candidate_pass_bins": int(adata.obs["passes_qc_candidate"].sum()),
            "candidate_pass_fraction": float(adata.obs["passes_qc_candidate"].mean()),
            "checkpoint": str(checkpoint),
        }
        write_json_atomic(outdir / "qc_summary.json", summary)
        state.update({"status": "complete", "elapsed_seconds": round(time.time() - started, 2), "qc": summary})
        write_json_atomic(state_path, state)
        print(json.dumps(summary, indent=2))
    except Exception as exc:
        state.update({"status": "failed", "error": repr(exc), "elapsed_seconds": round(time.time() - started, 2)})
        write_json_atomic(state_path, state)
        raise


if __name__ == "__main__":
    main()
