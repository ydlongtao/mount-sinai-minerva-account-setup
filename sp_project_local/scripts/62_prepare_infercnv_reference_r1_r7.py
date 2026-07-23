#!/usr/bin/env python3
"""Prepare balanced normal-reference candidates for infercnvpy."""
from __future__ import annotations

import json
import os
from pathlib import Path

import anndata as ad
import infercnvpy as cnv
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


PROJECT = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
ROOT = PROJECT / "results" / "segmented_official_v1"
REF = Path(os.environ.get("SP_HG38_REF", "/sc/arion/work/huangl21/references/hg38/gencode38"))
GTF = REF / "gencode.v38.primary_assembly.annotation.gtf.gz"
OUT = ROOT / "infercnv_reference_r1_r7"
SAMPLES = ["SC000895-R1", "SC000895-R7"]
MAIN_TYPES = {"fibroblast_stromal", "smooth_muscle_pericyte", "endothelial"}
AUX_TYPES = {"T_NK", "B_plasma", "myeloid"}
MAX_MAIN_PER_TYPE_SAMPLE = 3000
MAX_AUX_PER_TYPE_SAMPLE = 1000


def h5ad_path(sample: str) -> Path:
    return ROOT / "analysis_coarse" / "samples" / sample / f"{sample}_official_cell_level_analysis_coarse.h5ad"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    all_rows = []
    selected = []
    gtf_reports = []
    for sample in SAMPLES:
        path = h5ad_path(sample)
        a = sc.read_h5ad(path, backed="r")
        marker_path = ROOT / "preannotation_prostate" / "samples" / sample / "cluster_prostate_marker_spatial_preannotation.csv"
        marker = pd.read_csv(marker_path, dtype={"cluster": str})
        marker["cluster"] = marker["cluster"].astype(str)
        marker = marker.set_index("cluster")
        clusters = a.obs["cell_level_leiden_coarse"].astype(str).to_numpy()
        provisional = np.array([marker["predicted_type"].get(c, "ambiguous_or_low_signal") for c in clusters], dtype=object)
        confidence = np.array([marker["confidence"].get(c, "ambiguous") for c in clusters], dtype=object)
        niche_table = pd.read_csv(ROOT / "spatial_niche_r1_r7" / "samples" / sample / "niche_cells.csv.gz", usecols=["niche"])
        if len(niche_table) != a.n_obs:
            raise ValueError(f"{sample}: niche/H5AD row mismatch")
        qc = a.obs
        quality = ((qc["total_counts"].to_numpy(float) >= 20) &
                   (qc["n_genes_by_counts"].to_numpy(float) >= 10) &
                   (qc["pct_counts_mt"].to_numpy(float) <= 30))
        base = pd.DataFrame({"sample_id": sample, "cell_id": a.obs_names.astype(str),
                             "cluster": clusters, "niche": niche_table["niche"].to_numpy(int),
                             "provisional_type": provisional, "confidence": confidence,
                             "total_counts": qc["total_counts"].to_numpy(float),
                             "n_genes_by_counts": qc["n_genes_by_counts"].to_numpy(float),
                             "pct_counts_mt": qc["pct_counts_mt"].to_numpy(float),
                             "quality_pass": quality})
        base["reference_group"] = "exclude"
        base.loc[base["provisional_type"].isin(MAIN_TYPES) & base["quality_pass"] & base["confidence"].ne("ambiguous"), "reference_group"] = "main_reference"
        base.loc[base["provisional_type"].isin(AUX_TYPES) & base["quality_pass"] & base["confidence"].ne("ambiguous"), "reference_group"] = "auxiliary_reference"
        base["selected_for_pilot"] = False
        rng = np.random.default_rng(20260722)
        for group, cap in [("main_reference", MAX_MAIN_PER_TYPE_SAMPLE), ("auxiliary_reference", MAX_AUX_PER_TYPE_SAMPLE)]:
            for typ in sorted(MAIN_TYPES | AUX_TYPES):
                idx = base.index[(base["reference_group"] == group) & (base["provisional_type"] == typ)].to_numpy()
                if len(idx):
                    take = idx if len(idx) <= cap else np.sort(rng.choice(idx, cap, replace=False))
                    base.loc[take, "selected_for_pilot"] = True
                    selected.append(base.loc[take])
        all_rows.append(base)
        # Validate that the selected GTF can annotate this H5AD's gene symbols.
        var = pd.DataFrame(index=pd.Index(a.var_names.astype(str), name="gene"))
        mini = ad.AnnData(X=sparse.csr_matrix((1, len(var))), var=var)
        cnv.io.genomic_position_from_gtf(GTF, adata=mini, gtf_gene_id="gene_name", inplace=True)
        annotated = int(mini.var["chromosome"].notna().sum()) if "chromosome" in mini.var else 0
        gtf_reports.append({"sample_id": sample, "genes": int(a.n_vars), "genes_with_hg38_position": annotated,
                            "coverage_fraction": float(annotated / a.n_vars), "gtf": str(GTF)})
    candidates = pd.concat(all_rows, ignore_index=True)
    selected_df = pd.concat(selected, ignore_index=True) if selected else candidates.iloc[0:0].copy()
    candidates.to_csv(OUT / "all_reference_candidates.csv.gz", index=False, compression="gzip")
    selected_df.to_csv(OUT / "reference_cells_for_infercnvpy.csv.gz", index=False, compression="gzip")
    summary = candidates.groupby(["sample_id", "reference_group", "provisional_type"], dropna=False).size().reset_index(name="candidate_cells")
    selected_summary = selected_df.groupby(["sample_id", "reference_group", "provisional_type"], dropna=False).size().reset_index(name="selected_cells") if len(selected_df) else pd.DataFrame()
    summary.to_csv(OUT / "reference_candidate_summary.csv", index=False)
    if len(selected_summary):
        selected_summary.to_csv(OUT / "reference_selected_summary.csv", index=False)
    report = {"status": "pass", "gtf_validation": gtf_reports, "main_reference_types": sorted(MAIN_TYPES),
              "auxiliary_reference_types": sorted(AUX_TYPES), "max_main_per_type_sample": MAX_MAIN_PER_TYPE_SAMPLE,
              "max_aux_per_type_sample": MAX_AUX_PER_TYPE_SAMPLE, "candidate_cells": int(len(candidates)),
              "selected_cells": int(len(selected_df)), "excluded_cells": int((candidates.reference_group == "exclude").sum()),
              "selection_rule": "non-ambiguous stromal/contractile/endothelial as main; non-ambiguous immune as auxiliary; QC total_counts>=20, genes>=10, pct_mt<=30; deterministic balanced cap"}
    (OUT / "infercnv_reference_preparation.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
