#!/usr/bin/env python3
"""Run an infercnvpy pilot for R1/R7 using marker-selected main references.

中文说明：
1. 只把非恶性候选的内皮、基质/成纤维和平滑肌/周细胞作为主参考。
2. 免疫细胞不混入主参考，后续作为敏感性分析。
3. 保留 counts layer；infercnvpy 使用 counts 重新归一化并 log1p 后的 X。
4. 默认每个样本最多抽取 40,000 个 query cells，保留全部 selected reference cells。
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad
import infercnvpy as cnv
import numpy as np
import pandas as pd
import scanpy as sc


PROJECT = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
ROOT = PROJECT / "results" / "segmented_official_v1"
REF_DIR = ROOT / "infercnv_reference_r1_r7"
HG38 = Path(os.environ.get("SP_HG38_REF", "/sc/arion/work/huangl21/references/hg38/gencode38"))
GTF = HG38 / "gencode.v38.primary_assembly.annotation.gtf.gz"
INPUT_DIR = ROOT / "analysis_coarse" / "samples"
OUT_ROOT = ROOT / "infercnv_pilot_r1_r7"
MAIN_TYPES = {"endothelial", "fibroblast_stromal", "smooth_muscle_pericyte"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True, choices=["SC000895-R1", "SC000895-R7"])
    p.add_argument("--max-query", type=int, default=40000)
    p.add_argument("--seed", type=int, default=20260722)
    p.add_argument("--output-root", default=str(PROJECT / "results" / "segmented_official_v1" / "infercnv_pilot_r1_r7"))
    return p.parse_args()


def input_path(sample: str) -> Path:
    return INPUT_DIR / sample / f"{sample}_official_cell_level_analysis_coarse.h5ad"


def main() -> None:
    args = parse_args()
    sample = args.sample
    out = Path(args.output_root) / sample
    out.mkdir(parents=True, exist_ok=True)

    selected_path = REF_DIR / "reference_cells_for_infercnvpy.csv.gz"
    selected = pd.read_csv(selected_path)
    selected = selected[(selected["sample_id"] == sample) &
                        (selected["reference_group"] == "main_reference") &
                        (selected["provisional_type"].isin(MAIN_TYPES)) &
                        selected["selected_for_pilot"].astype(bool)]
    selected_ids = set(selected["cell_id"].astype(str))
    if not selected_ids:
        raise RuntimeError(f"{sample}: no selected main-reference cells")

    adata = sc.read_h5ad(input_path(sample))
    adata.var_names = adata.var_names.astype(str)
    adata.var_names_make_unique()
    adata.obs_names = adata.obs_names.astype(str)
    if "counts" not in adata.layers:
        raise RuntimeError(f"{sample}: counts layer is required for CNV pilot")

    present_refs = np.array([cell in selected_ids for cell in adata.obs_names], dtype=bool)
    if present_refs.sum() < 1000:
        raise RuntimeError(f"{sample}: only {present_refs.sum()} reference cells found in H5AD")

    query_idx = np.flatnonzero(~present_refs)
    if args.max_query > 0 and len(query_idx) > args.max_query:
        rng = np.random.default_rng(args.seed)
        query_idx = np.sort(rng.choice(query_idx, args.max_query, replace=False))
    keep = np.sort(np.concatenate([np.flatnonzero(present_refs), query_idx]))
    adata = adata[keep].copy()
    adata.obs["cnv_reference_label"] = np.where(
        adata.obs_names.isin(selected_ids), "main_reference", "query"
    ).astype(str)
    adata.obs["cnv_reference_label"] = pd.Categorical(
        adata.obs["cnv_reference_label"], categories=["main_reference", "query"]
    )

    # Use the preserved raw counts layer as input; infercnvpy operates on X.
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
    sc.pp.log1p(adata)
    cnv.io.genomic_position_from_gtf(
        GTF, adata=adata, gtf_gene_id="gene_name", adata_gene_id=None, inplace=True
    )
    if "chromosome" not in adata.var:
        raise RuntimeError(f"{sample}: GTF did not add chromosome positions")
    positioned = adata.var["chromosome"].notna().to_numpy()
    if positioned.sum() < 10000:
        raise RuntimeError(f"{sample}: only {positioned.sum()} genes have genomic positions")
    adata = adata[:, positioned].copy()

    cnv.tl.infercnv(
        adata,
        reference_key="cnv_reference_label",
        reference_cat="main_reference",
        window_size=100,
        step=10,
        dynamic_threshold=1.5,
        exclude_chromosomes=("chrX", "chrY"),
        chunksize=2000,
        n_jobs=max(1, int(os.environ.get("LSB_DJOB_NUMPROC", "4"))),
        key_added="cnv",
        calculate_gene_values=False,
    )

    # Group-level heatmap is compact and avoids plotting tens of thousands of cells.
    try:
        cnv.pl.chromosome_heatmap(
            adata,
            groupby="cnv_reference_label",
            use_rep="cnv",
            show=False,
            figsize=(18, 4),
        )
        plt.savefig(out / f"{sample}_cnv_reference_heatmap.png", dpi=180, bbox_inches="tight")
        plt.close("all")
    except Exception as exc:
        (out / "heatmap_error.txt").write_text(repr(exc) + "\n")

    result_h5ad = out / f"{sample}_infercnv_pilot.h5ad"
    adata.uns["infercnv_pilot_metadata"] = {
        "sample": sample,
        "reference_types": sorted(MAIN_TYPES),
        "reference_selection": "balanced marker-based non-ambiguous endothelial/stromal/smooth-muscle cells",
        "input_counts_layer": "counts",
        "normalization": "normalize_total_target_sum_1e4_then_log1p",
        "gtf": str(GTF),
        "max_query": args.max_query,
        "seed": args.seed,
        "exclude_chromosomes": ["chrX", "chrY"],
    }
    adata.write_h5ad(result_h5ad, compression="lzf")
    cnv_matrix = adata.obsm.get("X_cnv")
    summary = {
        "status": "pass",
        "sample": sample,
        "cells_total": int(adata.n_obs),
        "reference_cells": int((adata.obs["cnv_reference_label"] == "main_reference").sum()),
        "query_cells": int((adata.obs["cnv_reference_label"] == "query").sum()),
        "genes_with_positions": int(adata.n_vars),
        "cnv_matrix_shape": list(cnv_matrix.shape) if cnv_matrix is not None else None,
        "output_h5ad": str(result_h5ad),
        "heatmap": str(out / f"{sample}_cnv_reference_heatmap.png"),
        "reference_types": sorted(MAIN_TYPES),
    }
    (out / "pilot_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
