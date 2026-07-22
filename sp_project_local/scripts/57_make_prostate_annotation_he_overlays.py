#!/usr/bin/env python3
from __future__ import annotations

"""将前列腺癌细胞类型预注释叠加到官方 H&E hires 图像。

这些图用于人工核查，不是最终细胞类型标签。灰色表示低信号或暂不能判断。
"""
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
RESULTS = PROJECT / "results" / "segmented_official_v1"
RAW = Path(os.environ.get("SP_RAW", "/sc/arion/projects/DiseaseGeneCell/Huang_lab_data/SpatialTranscriptome"))
SAMPLES = ["SC000895-R1", "SC000895-R7"]
COLORS = {
    "malignant_epithelial_candidate": "#e31a1c", "luminal_epithelial": "#1f78b4",
    "basal_epithelial": "#33a02c", "T_NK": "#6a3d9a", "B_plasma": "#ff7f00",
    "myeloid": "#b15928", "fibroblast_stromal": "#a6cee3",
    "smooth_muscle_pericyte": "#fb9a99", "endothelial": "#cab2d6",
    "ambiguous_or_low_signal": "#bdbdbd",
}
DISPLAY = {
    "malignant_epithelial_candidate": "候选恶性上皮", "luminal_epithelial": "腔面上皮",
    "basal_epithelial": "基底上皮", "T_NK": "T/NK", "B_plasma": "B/浆细胞",
    "myeloid": "髓系", "fibroblast_stromal": "成纤维/间质",
    "smooth_muscle_pericyte": "平滑肌/周细胞", "endothelial": "内皮",
    "ambiguous_or_low_signal": "低信号/待核查",
}


def input_paths(sample: str) -> tuple[Path, str]:
    repaired = RESULTS / "analysis_coarse_repaired" / "samples" / sample / f"{sample}_official_cell_level_analysis_coarse_repaired.h5ad"
    if repaired.is_file():
        return repaired, "cell_level_leiden_coarse_repaired"
    return (RESULTS / "analysis_coarse" / "samples" / sample / f"{sample}_official_cell_level_analysis_coarse.h5ad",
            "cell_level_leiden_coarse")


def get_hires_coords(adata, image: np.ndarray) -> np.ndarray:
    if "spatial_hires" in adata.obsm:
        coords = np.asarray(adata.obsm["spatial_hires"], dtype=float)
    elif "spatial" in adata.obsm:
        coords = np.asarray(adata.obsm["spatial"], dtype=float)
        scale = None
        for payload in (adata.uns.get("spatial", {}) or {}).values():
            scale = (payload.get("scalefactors", {}) or {}).get("tissue_hires_scalef")
            if scale is not None:
                break
        if scale is not None and (coords[:, 0].max() > image.shape[1] or coords[:, 1].max() > image.shape[0]):
            coords = coords * float(scale)
    else:
        raise ValueError("H5AD lacks spatial coordinates")
    if coords.shape != (adata.n_obs, 2) or not np.isfinite(coords).all():
        raise ValueError("Invalid spatial coordinate array")
    return coords


def save_overlay(image, coords, labels, sample, output, title, focus=None):
    fig, ax = plt.subplots(figsize=(12, 12))
    ax.imshow(image)
    if focus is None:
        keep, size, alpha = np.ones(len(labels), dtype=bool), 1.0, 0.42
    else:
        keep, size, alpha = np.isin(labels, list(focus)), 2.2, 0.78
        ax.scatter(coords[:, 0], coords[:, 1], s=0.20, c="#777777", alpha=0.10,
                   linewidths=0, rasterized=True)
    for label, color in COLORS.items():
        mask = keep & (labels == label)
        if mask.any():
            ax.scatter(coords[mask, 0], coords[mask, 1], s=size, c=color, alpha=alpha,
                       linewidths=0, rasterized=True, label=DISPLAY[label])
    ax.set_xlim(0, image.shape[1]); ax.set_ylim(image.shape[0], 0)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"{sample}: {title}", fontsize=14)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=9, markerscale=4)
    fig.tight_layout(); fig.savefig(output, dpi=220, bbox_inches="tight"); plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", choices=SAMPLES, required=True)
    sample = parser.parse_args().sample
    h5ad, label_key = input_paths(sample)
    marker_csv = RESULTS / "preannotation_prostate" / "samples" / sample / "cluster_prostate_marker_spatial_preannotation.csv"
    image_path = (RAW / "wangy33.u.hpc.mssm.edu" / "10X_Single_Cell_RNA" /
                  "SC000895_Kuan_lin_Huang_T643" / sample / "outs" / "segmented_outputs" / "spatial" /
                  "tissue_hires_image.png")
    out = RESULTS / "preannotation_prostate" / "he_overlays" / sample
    out.mkdir(parents=True, exist_ok=True)
    for path in (h5ad, marker_csv, image_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    adata = sc.read_h5ad(h5ad, backed="r")
    if label_key not in adata.obs:
        raise KeyError(f"Missing cluster label: {label_key}")
    marker = pd.read_csv(marker_csv, dtype={"cluster": str})
    marker["cluster"] = marker["cluster"].astype(str)
    mapping = marker.set_index("cluster")["predicted_type"].to_dict()
    clusters = adata.obs[label_key].astype(str).to_numpy()
    labels = np.array([mapping.get(x, "ambiguous_or_low_signal") for x in clusters], dtype=object)
    image = np.asarray(Image.open(image_path).convert("RGB"))
    coords = get_hires_coords(adata, image)
    in_bounds = ((coords[:, 0] >= 0) & (coords[:, 0] < image.shape[1]) &
                 (coords[:, 1] >= 0) & (coords[:, 1] < image.shape[0]))
    out_of_bounds = int((~in_bounds).sum())
    coords = coords[in_bounds]
    labels = labels[in_bounds]
    save_overlay(image, coords, labels, sample, out / f"{sample}_he_overlay_all_labels.png", "H&E + prostate cancer preliminary labels")
    save_overlay(image, coords, labels, sample, out / f"{sample}_he_overlay_epithelial_focus.png", "H&E + epithelial focus",
                 {"malignant_epithelial_candidate", "luminal_epithelial", "basal_epithelial"})
    save_overlay(image, coords, labels, sample, out / f"{sample}_he_overlay_immune_focus.png", "H&E + immune focus",
                 {"T_NK", "B_plasma", "myeloid"})
    marker.to_csv(out / f"{sample}_cluster_label_review_table.csv", index=False)
    summary = {
        "sample_id": sample, "h5ad": str(h5ad), "h5ad_cluster_key": label_key,
        "he_image": str(image_path), "image_shape_yx": list(image.shape[:2]),
        "n_cells_total": int(adata.n_obs), "n_cells_overlay": int(len(labels)),
        "n_cells_outside_he_bounds": out_of_bounds,
        "coordinate_x_range": [float(coords[:, 0].min()), float(coords[:, 0].max())],
        "coordinate_y_range": [float(coords[:, 1].min()), float(coords[:, 1].max())],
        "label_counts_cells": {str(k): int(v) for k, v in pd.Series(labels).value_counts().items()},
        "outputs": [str(p) for p in sorted(out.glob("*.png"))],
        "interpretation_status": "preliminary_marker_spatial_only; H&E_requires_manual_review",
    }
    (out / f"{sample}_he_overlay_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
