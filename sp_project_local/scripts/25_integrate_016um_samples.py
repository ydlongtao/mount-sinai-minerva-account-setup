#!/usr/bin/env python3
from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

from batch_utils import BATCH_RESULTS, load_batch_config, make_basic_spatial_plot, read_manifest, sample_root, score_markers
from sp_v2_utils import write_h5ad_atomic


def main() -> None:
    config = load_batch_config()
    sample_ids = config.get("sample_ids") or [row["sample_id"] for row in read_manifest()]
    model_cfg = config.get("model", {})
    out_dir = BATCH_RESULTS / "integration" / "016um"
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    adatas = []
    used_samples = []
    missing = []
    for sample_id in sample_ids:
        path = sample_root(sample_id) / "016um" / f"{sample_id}_016um_processed.h5ad"
        if not path.exists():
            missing.append(str(path))
            continue
        adata = sc.read_h5ad(path)
        adata.obs_names = [f"{sample_id}:{x}" for x in adata.obs_names]
        adata.obs["sample_id"] = sample_id
        adatas.append(adata)
        used_samples.append(sample_id)
    if not adatas:
        raise RuntimeError("No 16um processed sample H5AD files were found.")

    integrated = sc.concat(adatas, join="inner", label="sample_id_from_concat", keys=used_samples, index_unique=None)
    integrated.var_names_make_unique()
    if "counts" in integrated.layers:
        integrated.X = integrated.layers["counts"].copy()
    sc.pp.calculate_qc_metrics(integrated, inplace=True, percent_top=None)
    sc.pp.normalize_total(integrated, target_sum=1e4)
    sc.pp.log1p(integrated)
    sc.pp.highly_variable_genes(integrated, n_top_genes=int(model_cfg.get("n_top_genes", 3000)), batch_key="sample_id")
    model = integrated[:, integrated.var["highly_variable"]].copy()
    sc.pp.scale(model, max_value=10)
    sc.tl.pca(model, n_comps=int(model_cfg.get("n_pcs", 30)), svd_solver="arpack")
    sc.pp.neighbors(model, n_neighbors=int(model_cfg.get("n_neighbors", 15)), n_pcs=int(model_cfg.get("n_pcs", 30)))
    sc.tl.umap(model, random_state=int(model_cfg.get("random_state", 0)))
    sc.tl.leiden(model, resolution=float(model_cfg.get("leiden_resolution", 0.4)), flavor="igraph", directed=False)
    integrated.obsm["X_pca"] = model.obsm["X_pca"]
    integrated.obsm["X_umap"] = model.obsm["X_umap"]
    integrated.obs["integrated_leiden"] = model.obs["leiden"].astype("category")
    score_markers(integrated, random_state=int(model_cfg.get("random_state", 0)))

    output_h5ad = out_dir / "integrated_016um.h5ad"
    write_h5ad_atomic(integrated, output_h5ad)

    sc.pl.umap(integrated, color=["sample_id", "integrated_leiden"], wspace=0.35, show=False)
    plt.savefig(fig_dir / "integrated_016um_umap.png", dpi=180, bbox_inches="tight")
    plt.close()
    for sample_id in used_samples[:]:
        subset = integrated[integrated.obs["sample_id"] == sample_id].copy()
        make_basic_spatial_plot(subset, "integrated_leiden", fig_dir / f"{sample_id}_016um_integrated_leiden.png", title=f"{sample_id} 16um integrated Leiden")

    summary = {
        "status": "pass",
        "h5ad": str(output_h5ad),
        "samples": used_samples,
        "missing_inputs": missing,
        "shape": list(integrated.shape),
        "clusters": integrated.obs["integrated_leiden"].value_counts().sort_index().astype(int).to_dict(),
    }
    (out_dir / "integrated_016um_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
