#!/usr/bin/env python3
"""Reclassify R1/R7 epithelial cells using core markers plus OmicVerse i3 candidates."""
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
from PIL import Image
import scanpy as sc


PROJECT = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
ROOT = PROJECT / "results" / "segmented_official_v1"
I3_ROOT = ROOT / "omicverse_pyinfercnv_i3_epithelial_v2_r1_r7"
MARKER_ROOT = ROOT / "preannotation_prostate" / "samples"
RAW = Path(os.environ.get("SP_RAW", "/sc/arion/projects/DiseaseGeneCell/Huang_lab_data/SpatialTranscriptome"))
OUT = ROOT / "epithelial_cnv_corepanel_reclassification_r1_r7"
SAMPLES = ["SC000895-R1", "SC000895-R7"]
PANELS = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "luminal": ["AR", "KLK3", "KLK2", "ACPP", "NKX3-1", "MSMB"],
    "malignant": ["AMACR", "ERG", "PCA3", "MYC", "EZH2", "TPD52", "GOLM1", "FASN", "GSTP1"],
    "proliferation": ["MKI67", "TOP2A", "PCNA", "TYMS", "UBE2C", "CENPF", "BIRC5", "MCM2", "MCM6"],
    "basal": ["KRT5", "KRT14", "TP63", "KRT15", "KRT17"],
    "neuroendocrine": ["CHGA", "SYP", "NCAM1", "INSM1", "ASCL1", "DLL3"],
}
SCORE_THRESHOLD = 0.25
EPITHELIAL_DETECTION_THRESHOLD = 0.25


def hires_image(sample: str) -> Path:
    return RAW / "wangy33.u.hpc.mssm.edu" / "10X_Single_Cell_RNA" / "SC000895_Kuan_lin_Huang_T643" / sample / "outs" / "segmented_outputs" / "spatial" / "tissue_hires_image.png"


def panel_scores(adata, panels: dict[str, list[str]]) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    upper = {str(g).upper(): str(g) for g in adata.var_names}
    found = {name: [upper[g] for g in genes if g in upper] for name, genes in panels.items()}
    all_genes = sorted(set(g for genes in found.values() for g in genes))
    if "counts" not in adata.layers:
        raise RuntimeError("counts layer missing")
    x = adata[:, all_genes].layers["counts"][:]
    if hasattr(x, "toarray"):
        x = x.toarray()
    x = np.asarray(x, dtype=np.float32)
    totals = x.sum(axis=1)
    norm = np.divide(x, totals[:, None], out=np.zeros_like(x), where=totals[:, None] > 0) * 1e4
    logx = np.log1p(norm)
    pos = {g: i for i, g in enumerate(all_genes)}
    result = {}
    for name, genes in found.items():
        idx = [pos[g] for g in genes]
        result[f"score_{name}"] = logx[:, idx].mean(axis=1) if idx else np.zeros(adata.n_obs, dtype=np.float32)
        result[f"detect_{name}"] = (x[:, idx] > 0).mean(axis=1) if idx else np.zeros(adata.n_obs, dtype=np.float32)
    return pd.DataFrame(result, index=adata.obs_names.astype(str)), found


def classify(adata, scores: pd.DataFrame, malignant_clusters: set[str]) -> pd.DataFrame:
    clusters = adata.obs["cell_level_leiden_coarse"].astype(str).to_numpy()
    cnv = adata.obs["cnv_group"].astype(str).to_numpy()
    cnv_expanded = cnv == "epithelial_candidate"
    cnv_high = cnv_expanded & np.isin(clusters, list(malignant_clusters))
    epithelial = (scores.score_epithelial >= SCORE_THRESHOLD) & (scores.detect_epithelial >= EPITHELIAL_DETECTION_THRESHOLD)
    marker_malignant = epithelial & ((scores.score_malignant >= SCORE_THRESHOLD) | (scores.score_proliferation >= SCORE_THRESHOLD))
    marker_ne = epithelial & (scores.score_neuroendocrine >= SCORE_THRESHOLD)
    marker_luminal = epithelial & (scores.score_luminal >= SCORE_THRESHOLD)
    marker_basal = epithelial & (scores.score_basal >= SCORE_THRESHOLD)
    labels = np.full(adata.n_obs, "non_epithelial_or_low_signal", dtype=object)
    labels[epithelial] = "epithelial_ambiguous"
    labels[epithelial & marker_luminal] = "luminal_epithelial"
    labels[epithelial & marker_basal & ~marker_luminal] = "basal_epithelial"
    labels[marker_ne] = "neuroendocrine_like"
    labels[marker_malignant] = "malignant_marker_candidate"
    labels[cnv_expanded & marker_malignant] = "malignant_marker_CNV_candidate"
    labels[cnv_expanded & epithelial & ~marker_malignant] = "epithelial_CNV_candidate"
    labels[cnv_high & epithelial] = "high_conf_malignant_cluster_CNV_candidate"
    labels[cnv_high & marker_malignant] = "malignant_marker_CNV_candidate"
    out = scores.copy()
    out["cluster"] = clusters
    out["cnv_expanded"] = cnv_expanded
    out["cnv_high_confidence_cluster"] = cnv_high
    out["core_epithelial_positive"] = epithelial
    out["core_malignant_marker_positive"] = marker_malignant
    out["classification"] = labels
    return out


def overlay(image: np.ndarray, coords: np.ndarray, labels: np.ndarray, output: Path, sample: str) -> None:
    colors = {
        "malignant_marker_CNV_candidate": "#e31a1c",
        "high_conf_malignant_cluster_CNV_candidate": "#ff7f00",
        "epithelial_CNV_candidate": "#00d7d7",
        "malignant_marker_candidate": "#8e24aa",
        "luminal_epithelial": "#1565c0",
        "basal_epithelial": "#2e7d32",
        "neuroendocrine_like": "#6a1b9a",
        "epithelial_ambiguous": "#90a4ae",
    }
    fig, ax = plt.subplots(figsize=(13, 13))
    ax.imshow(image)
    order = ["epithelial_ambiguous", "luminal_epithelial", "basal_epithelial", "neuroendocrine_like", "epithelial_CNV_candidate", "high_conf_malignant_cluster_CNV_candidate", "malignant_marker_CNV_candidate", "malignant_marker_candidate"]
    for label in order:
        mask = labels == label
        if mask.any():
            size = 1.0 if label in {"epithelial_ambiguous", "epithelial_CNV_candidate", "luminal_epithelial"} else 3.2
            alpha = 0.34 if label == "epithelial_ambiguous" else 0.72
            ax.scatter(coords[mask, 0], coords[mask, 1], s=size, c=colors[label], alpha=alpha, linewidths=0, rasterized=True, label=f"{label} (n={int(mask.sum()):,})")
    ax.set_xlim(0, image.shape[1]); ax.set_ylim(image.shape[0], 0); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"{sample}: core-panel + i3 CNV epithelial reclassification")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), markerscale=4, fontsize=8)
    fig.tight_layout(); fig.savefig(output, dpi=220, bbox_inches="tight"); plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", choices=SAMPLES)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    summaries = []
    for sample in ([args.sample] if args.sample else SAMPLES):
        path = I3_ROOT / sample / f"{sample}_omicverse_pyinfercnv_i3.h5ad"
        adata = sc.read_h5ad(path)
        marker = pd.read_csv(MARKER_ROOT / sample / "cluster_prostate_marker_spatial_preannotation.csv", dtype={"cluster": str})
        malignant_clusters = set(marker.loc[marker.predicted_type.eq("malignant_epithelial_candidate"), "cluster"].astype(str))
        scores, found = panel_scores(adata, PANELS)
        table = classify(adata, scores, malignant_clusters)
        table["x_hires"] = np.asarray(adata.obsm["spatial_hires"])[:, 0]
        table["y_hires"] = np.asarray(adata.obsm["spatial_hires"])[:, 1]
        table.insert(0, "sample", sample)
        sample_out = OUT / sample; sample_out.mkdir(parents=True, exist_ok=True)
        table.to_csv(sample_out / f"{sample}_epithelial_cnv_corepanel_cells.csv.gz", index=True, compression="gzip")
        image = np.asarray(Image.open(hires_image(sample)).convert("RGB"))
        coords = table[["x_hires", "y_hires"]].to_numpy(float)
        inside = np.isfinite(coords).all(axis=1) & (coords[:, 0] >= 0) & (coords[:, 0] < image.shape[1]) & (coords[:, 1] >= 0) & (coords[:, 1] < image.shape[0])
        overlay(image, coords[inside], table.classification.to_numpy()[inside], sample_out / f"{sample}_corepanel_cnv_reclassification_he_overlay.png", sample)
        counts = table.classification.value_counts().to_dict()
        summaries.append({"sample": sample, "cells": int(len(table)), "cnv_expanded": int(table.cnv_expanded.sum()), "cnv_high_confidence_cluster": int(table.cnv_high_confidence_cluster.sum()), "core_epithelial": int(table.core_epithelial_positive.sum()), "core_malignant_marker": int(table.core_malignant_marker_positive.sum()), "overlap_expanded_core_malignant": int((table.cnv_expanded & table.core_malignant_marker_positive).sum()), "overlap_high_cnv_core_malignant": int((table.cnv_high_confidence_cluster & table.core_malignant_marker_positive).sum()), "class_counts": counts, "panel_genes_present": found, "thresholds": {"score_threshold": SCORE_THRESHOLD, "epithelial_detection_threshold": EPITHELIAL_DETECTION_THRESHOLD}})
    summary_path = OUT / "reclassification_summary.json"
    prior = json.loads(summary_path.read_text()) if summary_path.exists() else {"samples": []}
    by_sample = {item["sample"]: item for item in prior.get("samples", [])}
    by_sample.update({item["sample"]: item for item in summaries})
    summary_path.write_text(json.dumps({"status": "pass", "samples": [by_sample[s] for s in sorted(by_sample)], "logic": "high-confidence malignant cluster CNV and expanded epithelial CNV were combined with cell-level core-panel scores; marker-positive and CNV-supported classes are reported separately", "panels": PANELS}, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summaries, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
