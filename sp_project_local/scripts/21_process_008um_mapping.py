#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import scanpy as sc

from batch_utils import (
    BATCH_RESULTS,
    load_batch_config,
    make_basic_spatial_plot,
    marker_groups,
    sample_by_index,
    sample_root,
)
from sp_v2_utils import attach_barcode_spatial, write_h5ad_atomic


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Map 16 um labels and marker scores onto 8 um Visium HD bins.")
    parser.add_argument("--sample-index", type=int, default=int(__import__("os").environ.get("LSB_JOBINDEX", "1")))
    parser.add_argument("--sample-id")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    row = sample_by_index(args.sample_index) if args.sample_id is None else {"sample_id": args.sample_id}
    sample_id = row["sample_id"]
    config = load_batch_config()
    model_cfg = config.get("model", {})
    random_state = int(model_cfg.get("random_state", 0))

    if "outs_dir" not in row:
        from batch_utils import sample_by_id
        row = sample_by_id(sample_id)

    out_dir = sample_root(sample_id) / "008um"
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    h5_008 = Path(row["matrix_008um"])
    mapping_path = Path(row["outs_dir"]) / "barcode_mappings.parquet"
    h5_016_processed = sample_root(sample_id) / "016um" / f"{sample_id}_016um_processed.h5ad"
    if not h5_016_processed.exists():
        raise FileNotFoundError(f"Run 20_process_016um_each_sample.py first: {h5_016_processed}")

    adata8 = sc.read_10x_h5(h5_008)
    adata8.var_names_make_unique()
    attach_barcode_spatial(adata8)
    adata8.obs["sample_id"] = sample_id
    adata8.obs["bin_size_um"] = 8

    adata16 = sc.read_h5ad(h5_016_processed)
    if not mapping_path.exists():
        sc.pp.calculate_qc_metrics(adata8, inplace=True, percent_top=None)
        output_h5ad = out_dir / f"{sample_id}_008um_mapped.h5ad"
        adata8.uns["mapping_source"] = {
            "mapping_path": str(mapping_path),
            "processed_016um_h5ad": str(h5_016_processed),
            "mapped_008um_bins": 0,
            "total_008um_bins": int(adata8.n_obs),
            "status": "missing_barcode_mappings",
        }
        write_h5ad_atomic(adata8, output_h5ad)
        make_basic_spatial_plot(adata8, "total_counts", fig_dir / f"{sample_id}_008um_total_counts.png", title=f"{sample_id} 8um total counts")
        make_basic_spatial_plot(adata8, "n_genes_by_counts", fig_dir / f"{sample_id}_008um_n_genes_by_counts.png", title=f"{sample_id} 8um detected genes")
        summary = {
            "sample_id": sample_id,
            "status": "warning_missing_barcode_mappings",
            "h5ad": str(output_h5ad),
            "n_008um_bins": int(adata8.n_obs),
            "n_genes": int(adata8.n_vars),
            "mapped_008um_bins": 0,
            "mapping_fraction": 0.0,
            "missing_mapping_path": str(mapping_path),
        }
        (out_dir / "008um_mapping_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        print(json.dumps(summary, indent=2, sort_keys=True))
        return

    map_df = pd.read_parquet(mapping_path)
    if not {"square_008um", "square_016um"}.issubset(map_df.columns):
        raise ValueError(f"Unexpected mapping columns in {mapping_path}: {list(map_df.columns)}")
    map_df = map_df.loc[map_df["square_008um"].isin(adata8.obs_names)]
    map_df = map_df.loc[map_df["square_016um"].isin(adata16.obs_names)]
    map_df = map_df.drop_duplicates("square_008um").set_index("square_008um")

    columns = ["leiden", "total_counts", "n_genes_by_counts"] + [
        f"score_{name}" for name in marker_groups().keys() if f"score_{name}" in adata16.obs
    ]
    mapped = adata16.obs.loc[map_df["square_016um"], [c for c in columns if c in adata16.obs]].copy()
    mapped.index = map_df.index
    for col in mapped.columns:
        adata8.obs[f"mapped_016um_{col}"] = mapped[col].reindex(adata8.obs_names)

    adata8.uns["mapping_source"] = {
        "mapping_path": str(mapping_path),
        "processed_016um_h5ad": str(h5_016_processed),
        "mapped_008um_bins": int(map_df.shape[0]),
        "total_008um_bins": int(adata8.n_obs),
    }

    output_h5ad = out_dir / f"{sample_id}_008um_mapped.h5ad"
    write_h5ad_atomic(adata8, output_h5ad)

    plot_cols = ["mapped_016um_leiden"] + [
        c for c in adata8.obs.columns if c.startswith("mapped_016um_score_")
    ][:6]
    for col in plot_cols:
        if col in adata8.obs and adata8.obs[col].notna().any():
            make_basic_spatial_plot(adata8, col, fig_dir / f"{sample_id}_008um_{col}.png", title=f"{sample_id} 8um {col}")

    summary = {
        "sample_id": sample_id,
        "status": "pass",
        "h5ad": str(output_h5ad),
        "n_008um_bins": int(adata8.n_obs),
        "n_genes": int(adata8.n_vars),
        "mapped_008um_bins": int(map_df.shape[0]),
        "mapping_fraction": float(map_df.shape[0] / max(1, adata8.n_obs)),
    }
    (out_dir / "008um_mapping_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
