#!/usr/bin/env python3
"""Create H&E overlays for high-confidence malignant and expanded epithelial candidates."""
from __future__ import annotations

import json
import os
import argparse
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
OUT = ROOT / "omicverse_pyinfercnv_i3_epithelial_v2_r1_r7_he_overlay"
SAMPLES = ["SC000895-R1", "SC000895-R7"]


def image_path(sample: str) -> Path:
    return RAW / "wangy33.u.hpc.mssm.edu" / "10X_Single_Cell_RNA" / "SC000895_Kuan_lin_Huang_T643" / sample / "outs" / "segmented_outputs" / "spatial" / "tissue_hires_image.png"


def save_overlay(image: np.ndarray, coords: np.ndarray, mask: np.ndarray, output: Path, title: str, color: str, label: str, size: float, alpha: float) -> None:
    fig, ax = plt.subplots(figsize=(13, 13))
    ax.imshow(image)
    if mask.any():
        ax.scatter(coords[mask, 0], coords[mask, 1], s=size, c=color, alpha=alpha, linewidths=0, rasterized=True, label=f"{label} (n={int(mask.sum()):,})")
    ax.set_xlim(0, image.shape[1]); ax.set_ylim(image.shape[0], 0)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), markerscale=5, frameon=True)
    fig.tight_layout(); fig.savefig(output, dpi=220, bbox_inches="tight"); plt.close(fig)


def save_combined(image: np.ndarray, coords: np.ndarray, malignant: np.ndarray, candidate: np.ndarray, output: Path, sample: str) -> None:
    fig, ax = plt.subplots(figsize=(13, 13))
    ax.imshow(image)
    expanded_only = candidate & ~malignant
    if expanded_only.any():
        ax.scatter(coords[expanded_only, 0], coords[expanded_only, 1], s=1.1, c="#00d7d7", alpha=0.42, linewidths=0, rasterized=True, label=f"expanded epithelial only (n={int(expanded_only.sum()):,})")
    if malignant.any():
        ax.scatter(coords[malignant, 0], coords[malignant, 1], s=4.0, c="#e31a1c", alpha=0.90, linewidths=0, rasterized=True, label=f"high-confidence malignant (n={int(malignant.sum()):,})")
    ax.set_xlim(0, image.shape[1]); ax.set_ylim(image.shape[0], 0)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(f"{sample}: H&E + i3 epithelial candidate overlays")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), markerscale=4, frameon=True)
    fig.tight_layout(); fig.savefig(output, dpi=220, bbox_inches="tight"); plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", choices=SAMPLES)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for sample in ([args.sample] if args.sample else SAMPLES):
        sample_root = I3_ROOT / sample
        h5ad = sample_root / f"{sample}_omicverse_pyinfercnv_i3.h5ad"
        adata = sc.read_h5ad(h5ad, backed="r")
        coords = np.asarray(adata.obsm["spatial_hires"], dtype=np.float32)
        candidate = adata.obs["cnv_group"].astype(str).to_numpy() == "epithelial_candidate"
        marker = pd.read_csv(MARKER_ROOT / sample / "cluster_prostate_marker_spatial_preannotation.csv", dtype={"cluster": str})
        malignant_clusters = set(marker.loc[marker["predicted_type"].eq("malignant_epithelial_candidate"), "cluster"].astype(str))
        clusters = adata.obs["cell_level_leiden_coarse"].astype(str).to_numpy()
        malignant = np.isin(clusters, list(malignant_clusters)) & candidate
        image = np.asarray(Image.open(image_path(sample)).convert("RGB"))
        inside = np.isfinite(coords).all(axis=1)
        inside &= (coords[:, 0] >= 0) & (coords[:, 0] < image.shape[1]) & (coords[:, 1] >= 0) & (coords[:, 1] < image.shape[0])
        coords = coords[inside]; candidate = candidate[inside]; malignant = malignant[inside]
        sample_out = OUT / sample; sample_out.mkdir(parents=True, exist_ok=True)
        save_overlay(image, coords, malignant, sample_out / f"{sample}_he_high_confidence_malignant_overlay.png", f"{sample}: H&E + high-confidence malignant candidates", "#e31a1c", "high-confidence malignant", 3.5, 0.88)
        save_overlay(image, coords, candidate, sample_out / f"{sample}_he_expanded_epithelial_overlay.png", f"{sample}: H&E + expanded epithelial candidates", "#00d7d7", "expanded epithelial", 1.1, 0.46)
        save_combined(image, coords, malignant, candidate, sample_out / f"{sample}_he_i3_epithelial_combined_overlay.png", sample)
        rows.append({"sample": sample, "image_width": int(image.shape[1]), "image_height": int(image.shape[0]), "cells_in_image": int(inside.sum()), "malignant_cells": int(malignant.sum()), "expanded_epithelial_cells": int(candidate.sum()), "expanded_only_cells": int((candidate & ~malignant).sum()), "malignant_clusters": sorted(malignant_clusters)})
    summary_path = OUT / "overlay_summary.json"
    prior = json.loads(summary_path.read_text()) if summary_path.exists() else {"samples": []}
    by_sample = {row["sample"]: row for row in prior.get("samples", [])}
    by_sample.update({row["sample"]: row for row in rows})
    summary_path.write_text(json.dumps({"status": "pass", "samples": [by_sample[s] for s in sorted(by_sample)], "coordinate_key": "spatial_hires", "image": "official segmented_outputs/spatial/tissue_hires_image.png", "malignant_definition": "marker predicted_type == malignant_epithelial_candidate, restricted to cnv_group epithelial_candidate", "expanded_definition": "cnv_group == epithelial_candidate from OmicVerse i3 epithelial_v2"}, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
