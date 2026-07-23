#!/usr/bin/env python3
"""Apply the user-approved final epithelial grouping to R1/R7."""
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


PROJECT = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
ROOT = PROJECT / "results" / "segmented_official_v1"
SOURCE = ROOT / "epithelial_cnv_corepanel_reclassification_r1_r7"
OUT = ROOT / "final_epithelial_groups_r1_r7"
RAW = Path(os.environ.get("SP_RAW", "/sc/arion/projects/DiseaseGeneCell/Huang_lab_data/SpatialTranscriptome"))
SAMPLES = ["SC000895-R1", "SC000895-R7"]
MALIGNANT = {"malignant_marker_candidate", "malignant_marker_CNV_candidate"}
NORMAL = {"epithelial_CNV_candidate", "epithelial_ambiguous", "basal_epithelial"}


def image_path(sample: str) -> Path:
    return RAW / "wangy33.u.hpc.mssm.edu" / "10X_Single_Cell_RNA" / "SC000895_Kuan_lin_Huang_T643" / sample / "outs" / "segmented_outputs" / "spatial" / "tissue_hires_image.png"


def final_label(value: str) -> str:
    if value in MALIGNANT:
        return "malignant_epithelial"
    if value in NORMAL:
        return "normal_epithelial"
    return value


def overlay(image: np.ndarray, table: pd.DataFrame, output: Path, sample: str) -> None:
    colors = {"malignant_epithelial": "#e31a1c", "normal_epithelial": "#2ca25f", "luminal_epithelial": "#3182bd", "neuroendocrine_like": "#88419d", "epithelial_ambiguous": "#9e9e9e"}
    fig, ax = plt.subplots(figsize=(13, 13)); ax.imshow(image)
    for label, color in colors.items():
        mask = table.final_epithelial_group.to_numpy() == label
        if mask.any():
            size = 2.3 if label in {"malignant_epithelial", "normal_epithelial"} else 1.2
            alpha = 0.70 if label in {"malignant_epithelial", "normal_epithelial"} else 0.48
            ax.scatter(table.loc[mask, "x_hires"], table.loc[mask, "y_hires"], s=size, c=color, alpha=alpha, linewidths=0, rasterized=True, label=f"{label} (n={int(mask.sum()):,})")
    ax.set_xlim(0, image.shape[1]); ax.set_ylim(image.shape[0], 0); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"{sample}: final epithelial grouping")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), markerscale=4, fontsize=9)
    fig.tight_layout(); fig.savefig(output, dpi=220, bbox_inches="tight"); plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--sample", choices=SAMPLES); args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); summaries = []
    for sample in ([args.sample] if args.sample else SAMPLES):
        source = SOURCE / sample / f"{sample}_epithelial_cnv_corepanel_cells.csv.gz"
        table = pd.read_csv(source, index_col=0)
        table["final_epithelial_group"] = table["classification"].map(final_label)
        out_dir = OUT / sample; out_dir.mkdir(parents=True, exist_ok=True)
        table.to_csv(out_dir / f"{sample}_final_epithelial_groups.csv.gz", compression="gzip")
        image = np.asarray(Image.open(image_path(sample)).convert("RGB"))
        coords = table[["x_hires", "y_hires"]].to_numpy(float)
        inside = np.isfinite(coords).all(axis=1) & (coords[:, 0] >= 0) & (coords[:, 0] < image.shape[1]) & (coords[:, 1] >= 0) & (coords[:, 1] < image.shape[0])
        overlay(image, table.loc[inside].copy(), out_dir / f"{sample}_final_epithelial_groups_he_overlay.png", sample)
        summaries.append({"sample": sample, "cells": int(len(table)), "final_group_counts": table["final_epithelial_group"].value_counts().to_dict(), "merged_malignant_source": sorted(MALIGNANT), "merged_normal_source": sorted(NORMAL)})
    summary_path = OUT / "final_epithelial_groups_summary.json"
    prior = json.loads(summary_path.read_text()) if summary_path.exists() else {"samples": []}
    by_sample = {x["sample"]: x for x in prior.get("samples", [])}; by_sample.update({x["sample"]: x for x in summaries})
    summary_path.write_text(json.dumps({"status": "pass", "samples": [by_sample[x] for x in sorted(by_sample)], "final_group_rules": {"malignant_epithelial": sorted(MALIGNANT), "normal_epithelial": sorted(NORMAL), "preserved": ["luminal_epithelial", "neuroendocrine_like", "non_epithelial_or_low_signal"]}}, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summaries, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
