#!/usr/bin/env python3
"""Step 4: spatial niche analysis using the final epithelial grouping."""
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
from scipy import sparse
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors


PROJECT = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
ROOT = PROJECT / "results" / "segmented_official_v1"
FINAL = ROOT / "final_epithelial_groups_r1_r7"
OLD_NICHE = ROOT / "spatial_niche_r1_r7"
RAW = Path(os.environ.get("SP_RAW", "/sc/arion/projects/DiseaseGeneCell/Huang_lab_data/SpatialTranscriptome"))
OUT = ROOT / "step4_final_groups_spatial_niche_r1_r7"
SAMPLES = ["SC000895-R1", "SC000895-R7"]
TYPES = ["malignant_epithelial", "normal_epithelial", "luminal_epithelial", "neuroendocrine_like", "T_NK", "B_plasma", "myeloid", "fibroblast_stromal", "smooth_muscle_pericyte", "endothelial"]
COLORS = ["#d73027", "#1a9850", "#4575b4", "#762a83", "#e66101", "#5e3c99", "#b2182b", "#66a61e", "#a6761d", "#1b9e77"]


def image_path(sample: str) -> Path:
    return RAW / "wangy33.u.hpc.mssm.edu" / "10X_Single_Cell_RNA" / "SC000895_Kuan_lin_Huang_T643" / sample / "outs" / "segmented_outputs" / "spatial" / "tissue_hires_image.png"


def effective_labels(final: pd.DataFrame, old: pd.DataFrame) -> np.ndarray:
    labels = old["cell_type_provisional"].astype(str).to_numpy()
    final_labels = final["final_epithelial_group"].astype(str).to_numpy()
    keep = {"malignant_epithelial", "normal_epithelial", "luminal_epithelial", "neuroendocrine_like"}
    labels[np.isin(final_labels, list(keep))] = final_labels[np.isin(final_labels, list(keep))]
    return labels


def composition(coords: np.ndarray, labels: np.ndarray, k: int = 12) -> tuple[np.ndarray, np.ndarray]:
    nn = NearestNeighbors(n_neighbors=k + 1, n_jobs=-1).fit(coords)
    idx = nn.kneighbors(return_distance=False)[:, 1:]
    codes = {x: i for i, x in enumerate(TYPES)}
    code = np.array([codes.get(x, -1) for x in labels])
    row = np.repeat(np.arange(len(labels)), k); col = idx.ravel()
    adj = sparse.csr_matrix((np.ones(len(row)), (row, col)), shape=(len(labels), len(labels)))
    valid = code >= 0
    rr = np.flatnonzero(valid); one = sparse.csr_matrix((np.ones(len(rr)), (rr, code[rr])), shape=(len(labels), len(TYPES)))
    counts = (adj @ one).toarray(); totals = counts.sum(axis=1)
    valid = totals >= 3
    prop = np.divide(counts, totals[:, None], out=np.zeros_like(counts), where=totals[:, None] > 0)
    return prop, valid


def overlay(image: np.ndarray, coords: np.ndarray, niche: np.ndarray, labels: np.ndarray, output: Path, sample: str) -> None:
    fig, ax = plt.subplots(figsize=(13, 13)); ax.imshow(image)
    niche_colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"]
    for n in sorted(set(niche)):
        mask = niche == n
        color = niche_colors[int(n)] if 0 <= int(n) < len(niche_colors) else "#bdbdbd"
        ax.scatter(coords[mask, 0], coords[mask, 1], s=0.8, c=color, alpha=0.32, linewidths=0, rasterized=True)
    malignant = labels == "malignant_epithelial"
    normal = labels == "normal_epithelial"
    ax.scatter(coords[normal, 0], coords[normal, 1], s=1.8, facecolors="none", edgecolors="#00ff66", linewidths=0.18, rasterized=True, label=f"normal_epithelial (n={int(normal.sum()):,})")
    ax.scatter(coords[malignant, 0], coords[malignant, 1], s=3.0, c="#ff1f1f", alpha=0.88, linewidths=0, rasterized=True, label=f"malignant_epithelial (n={int(malignant.sum()):,})")
    ax.set_xlim(0, image.shape[1]); ax.set_ylim(image.shape[0], 0); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(f"{sample}: final groups + spatial niche")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), markerscale=4); fig.tight_layout(); fig.savefig(output, dpi=220, bbox_inches="tight"); plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--sample", choices=SAMPLES); args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); summaries = []
    # Fit one pooled niche model so R1/R7 niche IDs are comparable.
    frames = []
    for sample in ([args.sample] if args.sample else SAMPLES):
        final = pd.read_csv(FINAL / sample / f"{sample}_final_epithelial_groups.csv.gz")
        old = pd.read_csv(OLD_NICHE / "samples" / sample / "niche_cells.csv.gz")
        if len(final) != len(old): raise ValueError(f"{sample}: final and niche rows differ")
        labels = effective_labels(final, old); coords = final[["x_hires", "y_hires"]].to_numpy(float)
        prop, valid = composition(coords, labels)
        frames.append((sample, final, old, labels, coords, prop, valid))
    model = KMeans(n_clusters=5, n_init=20, random_state=20260723).fit(np.vstack([x[5][x[6]] for x in frames]))
    for sample, final, old, labels, coords, prop, valid in frames:
        niche = np.full(len(labels), -1, dtype=int); niche[valid] = model.predict(prop[valid])
        out = OUT / sample; out.mkdir(parents=True, exist_ok=True)
        cell = pd.DataFrame({"sample": sample, "x_hires": coords[:, 0], "y_hires": coords[:, 1], "final_epithelial_group": final.final_epithelial_group, "effective_cell_type": labels, "niche": niche})
        cell.to_csv(out / f"{sample}_step4_niche_cells.csv.gz", index=False, compression="gzip")
        comp = pd.DataFrame(prop, columns=TYPES); comp["niche"] = niche
        means = comp[comp.niche >= 0].groupby("niche")[TYPES].mean(); means.to_csv(out / f"{sample}_step4_niche_composition.csv")
        image = np.asarray(Image.open(image_path(sample)).convert("RGB")); inside = (coords[:, 0] >= 0) & (coords[:, 0] < image.shape[1]) & (coords[:, 1] >= 0) & (coords[:, 1] < image.shape[0])
        overlay(image, coords[inside], niche[inside], labels[inside], out / f"{sample}_step4_niche_he_overlay.png", sample)
        summaries.append({"sample": sample, "cells": len(labels), "valid_niche_cells": int(valid.sum()), "niches": int(niche[niche >= 0].max() + 1 if (niche >= 0).any() else 0), "final_group_counts": final.final_epithelial_group.value_counts().to_dict(), "effective_type_counts": pd.Series(labels).value_counts().to_dict()})
    summary_path = OUT / "step4_summary.json"; prior = json.loads(summary_path.read_text()) if summary_path.exists() else {"samples": []}; by = {x["sample"]: x for x in prior.get("samples", [])}; by.update({x["sample"]: x for x in summaries})
    summary_path.write_text(json.dumps({"status": "pass", "samples": [by[x] for x in sorted(by)], "types": TYPES, "k_neighbors": 12, "niche_clusters": 5}, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summaries, indent=2, ensure_ascii=False))


if __name__ == "__main__": main()
