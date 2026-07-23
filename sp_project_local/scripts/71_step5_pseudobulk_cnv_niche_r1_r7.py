#!/usr/bin/env python3
"""Step 4 continuation + Step 5: niche CNV burden, pseudobulk DE, pathways."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


PROJECT = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
ROOT = PROJECT / "results" / "segmented_official_v1"
FINAL = ROOT / "final_epithelial_groups_r1_r7"
I3 = ROOT / "omicverse_pyinfercnv_i3_epithelial_v2_r1_r7"
STEP4 = ROOT / "step4_final_groups_spatial_niche_r1_r7"
OUT = ROOT / "step5_niche_cnv_pseudobulk_r1_r7"
SAMPLES = ["SC000895-R1", "SC000895-R7"]
GROUPS = ["malignant_epithelial", "normal_epithelial"]
PROGRAMS = {
    "epithelial_identity": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "prostate_luminal": ["AR", "KLK3", "KLK2", "ACPP", "NKX3-1", "MSMB"],
    "malignant_program": ["AMACR", "ERG", "PCA3", "MYC", "EZH2", "TPD52", "GOLM1", "FASN", "GSTP1"],
    "cell_cycle": ["MKI67", "TOP2A", "PCNA", "TYMS", "UBE2C", "CENPF", "BIRC5", "MCM2", "MCM6", "CDK1", "CCNB1"],
    "DNA_repair": ["BRCA1", "BRCA2", "RAD51", "ATM", "ATR", "CHEK1", "CHEK2", "BARD1", "MRE11"],
    "EMT_invasion": ["VIM", "COL1A1", "COL3A1", "MMP2", "MMP9", "ITGA6", "ITGB1", "CDH2", "SNAI1", "SNAI2"],
    "hypoxia": ["HIF1A", "EPAS1", "CA9", "VEGFA", "LDHA", "SLC2A1", "ADM", "NDRG1"],
    "inflammatory": ["IL6", "IL1B", "TNF", "NFKB1", "RELA", "CXCL8", "CCL2", "STAT3"],
    "neuroendocrine": ["CHGA", "SYP", "NCAM1", "INSM1", "ASCL1", "DLL3"],
}


def read_inputs(sample: str):
    final = pd.read_csv(FINAL / sample / f"{sample}_final_epithelial_groups.csv.gz")
    niche = pd.read_csv(STEP4 / sample / f"{sample}_step4_niche_cells.csv.gz")
    if len(final) != len(niche):
        raise ValueError(f"{sample}: final/niche row mismatch")
    groups = final["final_epithelial_group"].astype(str).to_numpy()
    niches = niche["niche"].to_numpy(int)
    h5ad = I3 / sample / f"{sample}_omicverse_pyinfercnv_i3.h5ad"
    adata = sc.read_h5ad(h5ad)
    if len(adata) != len(final):
        raise ValueError(f"{sample}: H5AD/final row mismatch")
    return adata, final, groups, niches


def pseudobulk_de(adata, groups: np.ndarray, sample: str) -> pd.DataFrame:
    X = adata.layers["counts"]
    names = np.asarray(adata.var_names.astype(str))
    masks = {g: groups == g for g in GROUPS}
    sums = {g: np.asarray(X[masks[g]].sum(axis=0)).ravel().astype(float) for g in GROUPS}
    totals = {g: max(v.sum(), 1.0) for g, v in sums.items()}
    cpm = {g: sums[g] / totals[g] * 1e6 for g in GROUPS}
    log2fc = np.log2((cpm["malignant_epithelial"] + 0.5) / (cpm["normal_epithelial"] + 0.5))
    df = pd.DataFrame({"sample": sample, "gene": names, "malignant_cells": int(masks["malignant_epithelial"].sum()), "normal_cells": int(masks["normal_epithelial"].sum()), "malignant_cpm": cpm["malignant_epithelial"], "normal_cpm": cpm["normal_epithelial"], "log2fc_malignant_vs_normal": log2fc})
    df = df[(df.malignant_cpm + df.normal_cpm) > 1].copy()
    df["abs_log2fc"] = df.log2fc_malignant_vs_normal.abs()
    return df.sort_values("abs_log2fc", ascending=False)


def program_scores(adata, groups: np.ndarray, sample: str) -> pd.DataFrame:
    names = {str(g).upper(): str(g) for g in adata.var_names}
    needed = sorted({names[g] for panel in PROGRAMS.values() for g in panel if g in names})
    x = adata[:, needed].layers["counts"]
    if sparse.issparse(x): x = x.toarray()
    x = np.asarray(x, dtype=float); totals = x.sum(axis=1); logx = np.log1p(np.divide(x, totals[:, None], out=np.zeros_like(x), where=totals[:, None] > 0) * 1e4)
    pos = {g: i for i, g in enumerate(needed)}; rows = []
    for group in GROUPS:
        mask = groups == group
        row = {"sample": sample, "group": group, "cells": int(mask.sum())}
        for pname, panel in PROGRAMS.items():
            found = [names[g] for g in panel if g in names]
            row[pname] = float(logx[np.ix_(mask, [pos[g] for g in found])].mean()) if found and mask.any() else 0.0
        rows.append(row)
    return pd.DataFrame(rows)


def niche_tables(adata, groups: np.ndarray, niches: np.ndarray, sample: str):
    rows = []
    cnv = np.asarray(adata.obsm["X_cnv"], dtype=float)
    burden = np.nanmean(np.abs(cnv), axis=1)
    high_fraction = np.mean(np.abs(cnv) >= 0.25, axis=1)
    for niche in sorted(np.unique(niches[niches >= 0])):
        in_niche = niches == niche; total = max(int(in_niche.sum()), 1)
        epi = in_niche & np.isin(groups, GROUPS); epi_total = max(int(epi.sum()), 1)
        for group in GROUPS:
            mask = in_niche & (groups == group)
            rows.append({"sample": sample, "niche": int(niche), "group": group, "cells": int(mask.sum()), "fraction_of_niche": float(mask.sum() / total), "fraction_of_epithelial_in_niche": float(mask.sum() / epi_total), "cnv_burden_mean_abs": float(burden[mask].mean()) if mask.any() else np.nan, "cnv_high_bin_fraction": float(high_fraction[mask].mean()) if mask.any() else np.nan})
    return pd.DataFrame(rows), burden


def plots(niche_df: pd.DataFrame, program_df: pd.DataFrame):
    OUT.mkdir(parents=True, exist_ok=True)
    pivot = niche_df.pivot_table(index=["sample", "niche"], columns="group", values="fraction_of_niche", fill_value=0)
    ax = pivot.plot(kind="bar", figsize=(12, 5), color=["#d73027", "#1a9850"]); ax.set_ylabel("fraction of niche"); ax.set_title("Malignant/normal epithelial fraction by niche"); ax.figure.tight_layout(); ax.figure.savefig(OUT / "niche_epithelial_composition.png", dpi=180); plt.close(ax.figure)
    b = niche_df.pivot_table(index=["sample", "niche"], columns="group", values="cnv_burden_mean_abs")
    ax = b.plot(kind="bar", figsize=(12, 5), color=["#d73027", "#1a9850"]); ax.set_ylabel("mean abs CNV"); ax.set_title("CNV burden by niche and final epithelial group"); ax.figure.tight_layout(); ax.figure.savefig(OUT / "niche_cnv_burden.png", dpi=180); plt.close(ax.figure)
    p = program_df.set_index(["sample", "group"])[list(PROGRAMS)]
    fig, ax = plt.subplots(figsize=(13, 4)); im = ax.imshow(p.to_numpy(), aspect="auto", cmap="RdBu_r"); ax.set_xticks(range(len(PROGRAMS)), list(PROGRAMS), rotation=45, ha="right"); ax.set_yticks(range(len(p)), [f"{a} {b}" for a,b in p.index]); fig.colorbar(im, ax=ax, label="mean log1p normalized score"); fig.tight_layout(); fig.savefig(OUT / "malignant_normal_program_scores.png", dpi=180); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument("--sample", choices=SAMPLES); args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); de_all=[]; programs=[]; niche_all=[]; summaries=[]
    for sample in ([args.sample] if args.sample else SAMPLES):
        adata, final, groups, niches = read_inputs(sample)
        de = pseudobulk_de(adata, groups, sample); de.to_csv(OUT / f"{sample}_malignant_vs_normal_pseudobulk.csv", index=False); de.head(100).to_csv(OUT / f"{sample}_malignant_vs_normal_top100.csv", index=False); de_all.append(de)
        prog = program_scores(adata, groups, sample); programs.append(prog)
        nt, burden = niche_tables(adata, groups, niches, sample); niche_all.append(nt)
        summaries.append({"sample": sample, "cells": int(len(groups)), "malignant_cells": int((groups == GROUPS[0]).sum()), "normal_cells": int((groups == GROUPS[1]).sum()), "niche_count": int((niches >= 0).sum()), "de_genes": int(len(de)), "top_up_genes": de.head(20).query("log2fc_malignant_vs_normal > 0").gene.tolist(), "top_down_genes": de.sort_values("log2fc_malignant_vs_normal").head(20).gene.tolist(), "mean_cnv_burden_malignant": float(burden[groups == GROUPS[0]].mean()), "mean_cnv_burden_normal": float(burden[groups == GROUPS[1]].mean())})
    niche_df = pd.concat(niche_all, ignore_index=True); program_df = pd.concat(programs, ignore_index=True); niche_df.to_csv(OUT / "niche_epithelial_cnv_burden.csv", index=False); program_df.to_csv(OUT / "malignant_normal_program_scores.csv", index=False); plots(niche_df, program_df)
    report = {"status": "pass", "samples": summaries, "note": "Pseudobulk is one aggregate per final group per sample; with one patient per condition, inferential patient-level p-values are not estimated.", "programs": PROGRAMS, "cnv_metric": "mean absolute X_cnv and fraction of bins with abs(X_cnv)>=0.25", "niche_metric": "12-neighbor local composition; pooled 5-cluster niche model from Step 4"}
    (OUT / "step5_summary.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__": main()
