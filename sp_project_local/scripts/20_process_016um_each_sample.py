#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

import scanpy as sc

from batch_utils import BATCH_RESULTS, load_batch_config, make_basic_spatial_plot, marker_groups, sample_by_index, sample_root, score_markers, write_json
from sp_v2_utils import attach_barcode_spatial, write_h5ad_atomic


def main() -> None:
    parser = argparse.ArgumentParser(description="Process one sample at 16um resolution for batch v1.")
    parser.add_argument("--sample-index", type=int, required=True)
    args = parser.parse_args()
    sample = sample_by_index(args.sample_index)
    sample_id = sample["sample_id"]
    config = load_batch_config()
    model_cfg = config["model"]
    outdir = sample_root(sample_id) / "016um"
    figdir = outdir / "figures"
    outdir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_10x_h5(sample["matrix_016um"], gex_only=True)
    adata.var_names_make_unique()
    adata.obs["sample_id"] = sample_id
    adata.obs["sample_index"] = args.sample_index
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)
    attach_barcode_spatial(adata)
    adata.layers["counts"] = adata.X.copy()
    adata = adata[adata.obs["total_counts"] > 0].copy()
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=model_cfg["n_top_genes"], flavor="seurat")
    model = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(model, max_value=10, zero_center=False)
    n_pcs = min(model_cfg["n_pcs"], model.n_obs - 1, model.n_vars - 1)
    sc.tl.pca(model, n_comps=n_pcs)
    sc.pp.neighbors(model, n_neighbors=model_cfg["n_neighbors"], n_pcs=n_pcs)
    sc.tl.umap(model, init_pos="random", random_state=model_cfg["random_state"])
    sc.tl.leiden(model, resolution=model_cfg["leiden_resolution"], key_added="leiden", flavor="igraph", directed=False)
    adata.obsm["X_pca"] = model.obsm["X_pca"]
    adata.obsm["X_umap"] = model.obsm["X_umap"]
    adata.obsp["connectivities"] = model.obsp["connectivities"]
    adata.obsp["distances"] = model.obsp["distances"]
    adata.obs["leiden"] = model.obs["leiden"].astype("category")
    present = score_markers(adata, marker_groups(), model_cfg["random_state"])
    adata.uns["batch_v1"] = {"sample": sample, "config": config, "marker_genes_present": present}

    output = outdir / f"{sample_id}_016um_processed.h5ad"
    write_h5ad_atomic(adata, output)
    for color in ["leiden", "total_counts", "n_genes_by_counts"]:
        make_basic_spatial_plot(adata, color, figdir / f"{sample_id}_spatial_{color}.png", f"{sample_id} {color}")
    for color in [c for c in adata.obs.columns if c.startswith("score_")]:
        make_basic_spatial_plot(adata, color, figdir / f"{sample_id}_spatial_{color}.png", f"{sample_id} {color}")
    summary = {
        "sample_id": sample_id,
        "shape": [adata.n_obs, adata.n_vars],
        "clusters": int(adata.obs["leiden"].nunique()),
        "output": str(output),
        "figures": str(figdir),
    }
    write_json(outdir / "016um_summary.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
