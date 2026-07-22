#!/usr/bin/env python3
"""Extract marker-score and spatial evidence for preliminary cell-type labels."""
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
GROUPS = ["luminal", "basal", "tumor", "immune", "stromal", "endothelial"]
MARKERS = {"luminal": ["KLK3", "ACPP", "AR", "NKX3-1", "KRT8", "KRT18"], "basal": ["KRT5", "KRT14", "TP63", "KRT15"], "tumor": ["AMACR", "ERG", "MYC", "MKI67", "TOP2A"], "immune": ["PTPRC", "CD3D", "CD3E", "MS4A1", "CD68", "LYZ"], "stromal": ["COL1A1", "COL1A2", "DCN", "LUM", "ACTA2"], "endothelial": ["PECAM1", "VWF", "KDR", "EMCN"]}


def main() -> None:
    sample = SAMPLES[int(os.environ.get("LSB_JOBINDEX", "1")) - 1]
    repaired = ROOT / "analysis_coarse_repaired" / "samples" / sample / f"{sample}_official_cell_level_analysis_coarse_repaired.h5ad"
    if repaired.exists():
        path, label_key, version = repaired, "cell_level_leiden_coarse_repaired", "coarse_repaired_v2"
    else:
        path, label_key, version = ROOT / "analysis_coarse" / "samples" / sample / f"{sample}_official_cell_level_analysis_coarse.h5ad", "cell_level_leiden_coarse", "coarse_v2"
    adata = sc.read_h5ad(path, backed="r")
    if label_key not in adata.obs:
        raise ValueError(f"{sample}: missing {label_key}")
    obs = adata.obs
    labels = obs[label_key].astype(str)
    rows = []
    for cluster, idx in labels.groupby(labels).groups.items():
        sub = obs.loc[idx]
        means = {g: float(sub[f"score_{g}"].mean()) if f"score_{g}" in sub else 0.0 for g in GROUPS}
        ranked = sorted(means.items(), key=lambda x: x[1], reverse=True)
        predicted, top_score = ranked[0]
        second_score = ranked[1][1]
        margin = top_score - second_score
        if cluster == "unassigned_zero_counts" or top_score <= 0 or margin < 0.05:
            confidence, predicted = "ambiguous", "ambiguous_or_low_signal"
        elif margin < 0.15:
            confidence = "provisional"
        else:
            confidence = "provisional_strong"
        item = {"sample_id": sample, "result_version": version, "cluster": cluster, "n_cells": int(len(sub)), "predicted_type": predicted, "confidence": confidence, "top_score": top_score, "second_score": second_score, "score_margin": margin}
        item.update({f"score_{g}": means[g] for g in GROUPS})
        if "spatial" in adata.obsm:
            coords = pd.DataFrame(np.asarray(adata.obsm["spatial"]), index=adata.obs_names, columns=["x", "y"]).loc[idx]
            item.update({"spatial_x_mean": float(coords.x.mean()), "spatial_y_mean": float(coords.y.mean()), "spatial_x_sd": float(coords.x.std()), "spatial_y_sd": float(coords.y.std())})
        rows.append(item)
    frame = pd.DataFrame(rows).sort_values(["sample_id", "n_cells"], ascending=[True, False])
    out = ROOT / "preannotation" / "samples" / sample
    out.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out / "cluster_marker_spatial_preannotation.csv", index=False)
    summary = {"sample_id": sample, "result_version": version, "input_h5ad": str(path), "clusters": int(len(frame)), "candidate_type_counts": frame["predicted_type"].value_counts().to_dict(), "ambiguous_clusters": int((frame["confidence"] == "ambiguous").sum()), "marker_groups": MARKERS}
    (out / "preannotation_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
