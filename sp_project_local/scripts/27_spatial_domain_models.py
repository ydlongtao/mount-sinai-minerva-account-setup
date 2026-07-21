#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os

import numpy as np
import scanpy as sc
from scipy import sparse
from sklearn.preprocessing import StandardScaler

from batch_utils import load_batch_config, make_basic_spatial_plot, sample_by_index, sample_by_id, sample_root, score_markers
from sp_v2_utils import write_h5ad_atomic


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a stable expression-plus-spatial Leiden domain model for one sample.")
    parser.add_argument("--sample-index", type=int, default=int(os.environ.get("LSB_JOBINDEX", "1")))
    parser.add_argument("--sample-id")
    parser.add_argument("--spatial-weight", type=float, default=3.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    row = sample_by_index(args.sample_index) if args.sample_id is None else sample_by_id(args.sample_id)
    sample_id = row["sample_id"]
    config = load_batch_config()
    model_cfg = config.get("model", {})
    out_dir = sample_root(sample_id) / "spatial_domains"
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    cell_qc_h5ad = sample_root(sample_id) / "cell_qc_marker" / f"{sample_id}_cell_level_qc_marker.h5ad"
    if cell_qc_h5ad.exists():
        adata = sc.read_h5ad(cell_qc_h5ad)
        input_h5ad = cell_qc_h5ad
    else:
        input_h5ad = sample_root(sample_id) / "016um" / f"{sample_id}_016um_processed.h5ad"
        if not input_h5ad.exists():
            raise FileNotFoundError(f"No cell-level or 16um input found for {sample_id}")
        adata = sc.read_h5ad(input_h5ad)

    if sparse.issparse(adata.X):
        counts = np.asarray(adata.X.sum(axis=1)).ravel()
    else:
        counts = np.asarray(adata.X.sum(axis=1)).ravel()
    if "counts" not in adata.layers and np.nanmax(counts) > 50:
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=min(int(model_cfg.get("n_top_genes", 3000)), max(1, adata.n_vars - 1)))
    model = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(model, max_value=10)
    n_comps = min(int(model_cfg.get("n_pcs", 30)), max(2, model.n_vars - 1), max(2, model.n_obs - 1))
    sc.tl.pca(model, n_comps=n_comps, svd_solver="arpack")
    coords = np.asarray(adata.obsm["spatial"], dtype=np.float32)
    coords_scaled = StandardScaler().fit_transform(coords)
    model.obsm["X_expr_spatial"] = np.column_stack([model.obsm["X_pca"], coords_scaled * float(args.spatial_weight)])
    sc.pp.neighbors(model, n_neighbors=int(model_cfg.get("n_neighbors", 15)), use_rep="X_expr_spatial")
    sc.tl.umap(model, random_state=int(model_cfg.get("random_state", 0)))
    sc.tl.leiden(model, resolution=float(model_cfg.get("leiden_resolution", 0.4)), flavor="igraph", directed=False, key_added="spatial_domain")

    adata.obsm["X_pca_spatial_domain"] = model.obsm["X_pca"]
    adata.obsm["X_umap_spatial_domain"] = model.obsm["X_umap"]
    adata.obs["spatial_domain"] = model.obs["spatial_domain"].astype("category")
    score_markers(adata, random_state=int(model_cfg.get("random_state", 0)))

    output_h5ad = out_dir / f"{sample_id}_spatial_domains.h5ad"
    write_h5ad_atomic(adata, output_h5ad)
    make_basic_spatial_plot(adata, "spatial_domain", fig_dir / f"{sample_id}_spatial_domain.png", title=f"{sample_id} spatial domains")
    for col in [c for c in adata.obs.columns if c.startswith("score_")][:6]:
        make_basic_spatial_plot(adata, col, fig_dir / f"{sample_id}_{col}.png", title=f"{sample_id} {col}")

    summary = {
        "sample_id": sample_id,
        "status": "pass",
        "input_h5ad": str(input_h5ad),
        "h5ad": str(output_h5ad),
        "shape": list(adata.shape),
        "spatial_weight": float(args.spatial_weight),
        "domains": adata.obs["spatial_domain"].value_counts().sort_index().astype(int).to_dict(),
        "method_note": "Leiden on PCA embedding augmented with scaled spatial coordinates; stable first-pass substitute for heavier spatial graph models.",
    }
    (out_dir / "spatial_domain_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
