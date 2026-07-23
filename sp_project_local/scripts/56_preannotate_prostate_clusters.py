#!/usr/bin/env python3
"""Prostate-cancer-aware preliminary annotation from cluster marker evidence."""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

PROJECT_HOME = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
ROOT = Path(os.environ.get("SP_OFFICIAL_RESULTS", PROJECT_HOME / "results" / "segmented_official_v1"))
SAMPLES = ["SC000895-R1", "SC000895-R2", "SC000895-R4", "SC000895-R5", "SC000895-R6", "SC000895-R7", "SC000895-R8", "SC000895-R9"]
PANELS = {
    "luminal_epithelial": ["KLK3", "ACPP", "AR", "NKX3-1", "KRT8", "KRT18", "KRT19", "EPCAM"],
    "basal_epithelial": ["KRT5", "KRT14", "TP63", "KRT15", "KRT17", "KRT23"],
    "malignant_epithelial_candidate": ["AMACR", "ERG", "PCA3", "MYC", "MKI67", "TOP2A", "EZH2", "TPD52"],
    "T_NK": ["CD3D", "CD3E", "TRBC1", "TRBC2", "CD4", "CD8A", "CD8B", "IL7R", "NKG7", "GNLY"],
    "B_plasma": ["CD79A", "MS4A1", "CD74", "HLA-DRA", "CD37", "MZB1", "JCHAIN", "CD38"],
    "myeloid": ["LYZ", "TYROBP", "FCER1G", "CTSS", "CD68", "LGALS3", "CSF1R"],
    "fibroblast_stromal": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "CFD", "C7"],
    "smooth_muscle_pericyte": ["ACTA2", "TAGLN", "MYH11", "RGS5", "CSPG4", "MCAM"],
    "endothelial": ["PECAM1", "VWF", "KDR", "EMCN", "ESAM", "RAMP2", "PLVAP"],
}


def main() -> None:
    sample = SAMPLES[int(os.environ.get("LSB_JOBINDEX", "1")) - 1]
    repaired = ROOT / "analysis_coarse_repaired" / "samples" / sample / f"{sample}_official_cell_level_analysis_coarse_repaired.h5ad"
    if repaired.exists():
        path, label_key, version = repaired, "cell_level_leiden_coarse_repaired", "coarse_repaired_v2"
    else:
        path, label_key, version = ROOT / "analysis_coarse" / "samples" / sample / f"{sample}_official_cell_level_analysis_coarse.h5ad", "cell_level_leiden_coarse", "coarse_v2"
    adata = sc.read_h5ad(path)
    labels = adata.obs[label_key].astype(str)
    upper = {str(g).upper(): str(g) for g in adata.var_names}
    genes = sorted({upper[g] for panel in PANELS.values() for g in panel if g in upper})
    if len(genes) < 15:
        raise ValueError(f"{sample}: only {len(genes)} prostate marker genes found")
    values = adata[:, genes].X
    if hasattr(values, "toarray"):
        values = values.toarray()
    expr = pd.DataFrame(np.asarray(values), columns=genes)
    expr["cluster"] = labels.to_numpy()
    means = expr.groupby("cluster", observed=True)[genes].mean()
    detection = expr.assign(**{g: expr[g] > 0 for g in genes}).groupby("cluster", observed=True)[genes].mean()
    z = (means - means.mean(axis=0)) / means.std(axis=0).replace(0, 1)
    rows = []
    for cluster in means.index:
        scores = {}
        present = {}
        for name, panel in PANELS.items():
            found = [upper[g] for g in panel if g in upper]
            scores[name] = float(z.loc[cluster, found].mean()) if found else 0.0
            present[name] = int(len(found))
        ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
        predicted, top = ranked[0]
        second = ranked[1][1]
        margin = top - second
        epithelial = max(scores["luminal_epithelial"], scores["basal_epithelial"], scores["malignant_epithelial_candidate"])
        if cluster == "unassigned_zero_counts" or top < 0.35 or margin < 0.15:
            predicted, confidence = "ambiguous_or_low_signal", "ambiguous"
        elif predicted == "malignant_epithelial_candidate" and epithelial < 0.25:
            predicted, confidence = "ambiguous_malignant_candidate", "provisional"
        elif margin < 0.30:
            confidence = "provisional"
        else:
            confidence = "provisional_strong"
        item = {"sample_id": sample, "result_version": version, "cluster": cluster, "n_cells": int((labels == cluster).sum()), "predicted_type": predicted, "confidence": confidence, "top_score": top, "second_score": second, "score_margin": margin, "epithelial_context_score": epithelial}
        item.update({f"score_{name}": value for name, value in scores.items()})
        item.update({f"marker_genes_{name}": present[name] for name in PANELS})
        for name, panel in PANELS.items():
            found = [upper[g] for g in panel if g in upper]
            item[f"detect_{name}"] = float(detection.loc[cluster, found].mean()) if found else 0.0
        if "spatial" in adata.obsm:
            idx = labels == cluster
            coords = np.asarray(adata.obsm["spatial"])[idx.to_numpy()]
            item.update({"spatial_x_mean": float(coords[:, 0].mean()), "spatial_y_mean": float(coords[:, 1].mean()), "spatial_x_sd": float(coords[:, 0].std()), "spatial_y_sd": float(coords[:, 1].std())})
        rows.append(item)
    frame = pd.DataFrame(rows).sort_values(["sample_id", "n_cells"], ascending=[True, False])
    out = ROOT / "preannotation_prostate" / "samples" / sample
    out.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out / "cluster_prostate_marker_spatial_preannotation.csv", index=False)
    summary = {"sample_id": sample, "result_version": version, "input_h5ad": str(path), "clusters": int(len(frame)), "candidate_type_counts": frame["predicted_type"].value_counts().to_dict(), "ambiguous_clusters": int((frame["confidence"] == "ambiguous").sum()), "panel_genes_present": genes, "panels": PANELS}
    (out / "prostate_preannotation_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
