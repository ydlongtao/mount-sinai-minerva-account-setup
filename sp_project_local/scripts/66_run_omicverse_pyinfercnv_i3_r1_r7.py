#!/usr/bin/env python3
"""Run OmicVerse py-inferCNV with the i3 HMM mode for R1/R7."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import omicverse as ov
import pandas as pd
from pyinfercnv.io.genome import load_gene_positions
import scanpy as sc

PROJECT = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
ROOT = PROJECT / "results" / "segmented_official_v1"
INPUT_DIR = ROOT / "analysis_coarse" / "samples"
REF_TABLE = ROOT / "infercnv_reference_r1_r7" / "reference_cells_for_infercnvpy.csv.gz"
MAIN_TYPES = {"endothelial", "fibroblast_stromal", "smooth_muscle_pericyte"}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True, choices=["SC000895-R1", "SC000895-R7"])
    p.add_argument("--output-root", default=str(ROOT / "omicverse_pyinfercnv_i3_r1_r7"))
    args = p.parse_args()
    sample = args.sample
    out = Path(args.output_root) / sample
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(INPUT_DIR / sample / f"{sample}_official_cell_level_analysis_coarse.h5ad")
    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)
    if "counts" not in adata.layers:
        raise RuntimeError(f"{sample}: counts layer is required")

    refs = pd.read_csv(REF_TABLE)
    refs = refs[(refs.sample_id == sample) &
                (refs.reference_group == "main_reference") &
                refs.provisional_type.isin(MAIN_TYPES) &
                refs.selected_for_pilot.astype(bool)]
    reference_ids = set(refs.cell_id.astype(str))

    marker_path = ROOT / "preannotation_prostate" / "samples" / sample / "cluster_prostate_marker_spatial_preannotation.csv"
    marker = pd.read_csv(marker_path, dtype={"cluster": str}).set_index("cluster")
    clusters = adata.obs["cell_level_leiden_coarse"].astype(str)
    malignant_clusters = set(marker.index[marker["predicted_type"].eq("malignant_epithelial_candidate")])
    is_reference = adata.obs_names.isin(reference_ids)
    is_malignant = clusters.isin(malignant_clusters).to_numpy()
    adata.obs["cnv_group"] = np.where(is_reference, "main_reference", np.where(is_malignant, "malignant_epithelial_candidate", "other"))
    adata.obs["cnv_group"] = pd.Categorical(adata.obs["cnv_group"], categories=["main_reference", "malignant_epithelial_candidate", "other"])

    # Use py-inferCNV's bundled GRCh38 coordinates. The counts layer remains raw.
    positions = load_gene_positions("hg38").drop_duplicates("gene_symbol").set_index("gene_symbol")
    gene_symbols = adata.var_names.astype(str)
    matched = positions.reindex(gene_symbols)
    adata.var["chromosome"] = matched["chromosome"].to_numpy()
    adata.var["start"] = matched["start"].to_numpy()
    adata.var["end"] = matched["end"].to_numpy()
    adata.var_names_make_unique()
    required = adata.var[["chromosome", "start", "end"]].notna().all(axis=1).to_numpy()
    required &= ~adata.var["chromosome"].isin(["chrX", "chrY", "chrM"]).to_numpy()
    adata = adata[:, required].copy()

    cnv = ov.single.CNV(adata, method="infercnv", layer="counts")
    cnv.run(
        reference_key="cnv_group",
        reference_cat="main_reference",
        platform="10x",
        exclude_chromosomes=("chrX", "chrY", "chrM"),
        HMM=True,
        HMM_type="i3",
        HMM_report_by="cell",
        analysis_mode="subclusters",
        denoise=True,
        cluster_by_groups=True,
        num_threads=max(1, int(os.environ.get("LSB_DJOB_NUMPROC", "8"))),
        random_state=20260723,
        verbose=True,
    )

    # Keep an explicit run record inside AnnData for downstream audits.
    adata.uns.setdefault("omicverse_pyinfercnv_i3", {})
    adata.uns["omicverse_pyinfercnv_i3"].update({
        "sample": sample, "backend": "omicverse.ov.single.CNV -> pyinfercnv",
        "pyinfercnv_hmm": "i3", "reference_group": "main_reference",
        "tumor_candidate_group": "malignant_epithelial_candidate",
        "reference_types": sorted(MAIN_TYPES), "gene_position_source": "pyinfercnv bundled hg38",
        "counts_layer": "counts", "exclude_chromosomes": ["chrX", "chrY", "chrM"],
    })

    # OmicVerse plotting helpers are optional; preserve a compact group heatmap when available.
    try:
        ov.pl.cnv_heatmap(adata, groupby="cnv_group", figsize=(12, 5), title=f"{sample} OmicVerse py-inferCNV i3")
        plt.savefig(figures / f"{sample}_pyinfercnv_i3_heatmap.png", dpi=180, bbox_inches="tight")
        plt.close("all")
    except Exception as exc:
        (out / "heatmap_error.txt").write_text(repr(exc) + "\n")

    try:
        ov.pl.cnv_summary(adata, groupby="cnv_group", subset="malignant_epithelial_candidate", figsize=(12, 3), title=f"{sample} malignant epithelial candidate CNV")
        plt.savefig(figures / f"{sample}_malignant_candidate_cnv_summary.png", dpi=180, bbox_inches="tight")
        plt.close("all")
    except Exception as exc:
        (out / "summary_plot_error.txt").write_text(repr(exc) + "\n")

    output = out / f"{sample}_omicverse_pyinfercnv_i3.h5ad"
    adata.write_h5ad(output, compression="lzf")
    hmm = adata.obsm.get("X_cnv_hmm_states_i3", adata.obsm.get("X_cnv_hmm_states"))
    summary = {
        "status": "pass", "sample": sample, "backend": "omicverse_pyinfercnv",
        "pyinfercnv_version": "0.2.0", "hmm_type": "i3",
        "cells_total": int(adata.n_obs), "reference_cells": int(is_reference.sum()),
        "malignant_epithelial_candidate_cells": int(is_malignant.sum()),
        "other_cells": int((~is_reference & ~is_malignant).sum()),
        "genes_with_positions": int(adata.n_vars),
        "cnv_matrix_shape": list(adata.obsm["X_cnv"].shape),
        "hmm_state_matrix_shape": list(hmm.shape) if hmm is not None else None,
        "cnv_regions_present": bool(adata.uns.get("cnv", {}).get("cnv_regions") is not None),
        "output_h5ad": str(output),
    }
    (out / "i3_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
