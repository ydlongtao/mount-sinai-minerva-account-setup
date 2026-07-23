#!/usr/bin/env python3
"""R1/R7 prostate cancer spatial niche analysis.

Uses spatial k-nearest-neighbor composition and H&E overlays. Niche labels are
exploratory ecological states, not final cell-type annotations.
"""
from __future__ import annotations

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
import scanpy as sc


PROJECT = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
RAW = Path(os.environ.get("SP_RAW", "/sc/arion/projects/DiseaseGeneCell/Huang_lab_data/SpatialTranscriptome"))
ROOT = PROJECT / "results" / "segmented_official_v1"
OUT = ROOT / "spatial_niche_r1_r7"
SAMPLES = ["SC000895-R1", "SC000895-R7"]
TYPES = ["luminal_epithelial", "basal_epithelial", "malignant_epithelial_candidate",
         "T_NK", "B_plasma", "myeloid", "fibroblast_stromal",
         "smooth_muscle_pericyte", "endothelial"]
TYPE_DISPLAY = {
    "luminal_epithelial": "luminal epithelium", "basal_epithelial": "basal epithelium",
    "malignant_epithelial_candidate": "candidate malignant epithelium", "T_NK": "T/NK",
    "B_plasma": "B/plasma", "myeloid": "myeloid", "fibroblast_stromal": "fibroblast/stroma",
    "smooth_muscle_pericyte": "smooth muscle/pericyte", "endothelial": "endothelial",
}
NICHE_COLORS = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#bdbdbd"]


def input_paths(sample: str) -> tuple[Path, str]:
    repaired = ROOT / "analysis_coarse_repaired" / "samples" / sample / f"{sample}_official_cell_level_analysis_coarse_repaired.h5ad"
    if repaired.is_file():
        return repaired, "cell_level_leiden_coarse_repaired"
    return ROOT / "analysis_coarse" / "samples" / sample / f"{sample}_official_cell_level_analysis_coarse.h5ad", "cell_level_leiden_coarse"


def load_sample(sample: str) -> dict:
    h5ad, label_key = input_paths(sample)
    marker_path = ROOT / "preannotation_prostate" / "samples" / sample / "cluster_prostate_marker_spatial_preannotation.csv"
    image_path = RAW / "wangy33.u.hpc.mssm.edu" / "10X_Single_Cell_RNA" / "SC000895_Kuan_lin_Huang_T643" / sample / "outs" / "segmented_outputs" / "spatial" / "tissue_hires_image.png"
    adata = sc.read_h5ad(h5ad, backed="r")
    marker = pd.read_csv(marker_path, dtype={"cluster": str})
    mapping = marker.assign(cluster=marker["cluster"].astype(str)).set_index("cluster")["predicted_type"].to_dict()
    clusters = adata.obs[label_key].astype(str).to_numpy()
    labels = np.array([mapping.get(x, "ambiguous_or_low_signal") for x in clusters], dtype=object)
    coords = np.asarray(adata.obsm["spatial_hires"], dtype=np.float32)
    image = np.asarray(Image.open(image_path).convert("RGB"))
    inside = ((coords[:, 0] >= 0) & (coords[:, 0] < image.shape[1]) &
              (coords[:, 1] >= 0) & (coords[:, 1] < image.shape[0]))
    if inside.mean() < 0.99:
        raise ValueError(f"{sample}: only {inside.mean():.3%} cells are inside H&E bounds")
    return {"sample": sample, "coords": coords, "labels": labels, "image": image,
            "inside": inside, "h5ad": str(h5ad), "label_key": label_key}


def spatial_composition(coords: np.ndarray, labels: np.ndarray, k: int = 12) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = len(labels)
    nn = NearestNeighbors(n_neighbors=k + 1, n_jobs=-1).fit(coords)
    indices = nn.kneighbors(return_distance=False)[:, 1:]
    rows = np.repeat(np.arange(n), k)
    cols = indices.ravel()
    adj = sparse.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)), shape=(n, n))
    codes = {name: i for i, name in enumerate(TYPES)}
    informative = np.array([codes.get(x, -1) for x in labels], dtype=np.int16)
    r = np.flatnonzero(informative >= 0)
    one_hot = sparse.csr_matrix((np.ones(len(r), dtype=np.float32), (r, informative[r])), shape=(n, len(TYPES)))
    counts = (adj @ one_hot).toarray()
    valid = counts.sum(axis=1) >= 3
    composition = np.divide(counts, counts.sum(axis=1, keepdims=True), out=np.zeros_like(counts), where=counts.sum(axis=1, keepdims=True) > 0)
    return composition, valid, informative


def edge_enrichment(labels: np.ndarray, coords: np.ndarray, k: int = 12) -> pd.DataFrame:
    nn = NearestNeighbors(n_neighbors=k + 1, n_jobs=-1).fit(coords)
    idx = nn.kneighbors(return_distance=False)[:, 1:]
    codes = {name: i for i, name in enumerate(TYPES)}
    code = np.array([codes.get(x, -1) for x in labels])
    valid = code >= 0
    obs = np.zeros((len(TYPES), len(TYPES)), dtype=np.int64)
    for i in np.flatnonzero(valid):
        js = idx[i][code[idx[i]] >= 0]
        for j in js:
            obs[code[i], code[j]] += 1
    freq = np.bincount(code[valid], minlength=len(TYPES)).astype(float) / max(valid.sum(), 1)
    expected = max(obs.sum(), 1) * np.outer(freq, freq)
    log2 = np.log2((obs + 0.5) / (expected + 0.5))
    return pd.DataFrame(log2, index=TYPES, columns=TYPES)


def overlay(sample: str, image: np.ndarray, coords: np.ndarray, niche: np.ndarray, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 12))
    ax.imshow(image)
    for n in sorted(set(niche)):
        mask = niche == n
        color = NICHE_COLORS[int(n)] if 0 <= int(n) < len(NICHE_COLORS) else "#bdbdbd"
        ax.scatter(coords[mask, 0], coords[mask, 1], s=1.0, c=color, alpha=0.50,
                   linewidths=0, rasterized=True,
                   label=(f"Niche {int(n) + 1}" if int(n) >= 0 else "low-signal / unassigned"))
    ax.set_xlim(0, image.shape[1]); ax.set_ylim(image.shape[0], 0)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"{sample}: H&E + spatial niches")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), markerscale=4)
    fig.tight_layout(); fig.savefig(out, dpi=220, bbox_inches="tight"); plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    data = [load_sample(s) for s in SAMPLES]
    features = []
    for d in data:
        comp, valid, informative = spatial_composition(d["coords"], d["labels"])
        d.update({"composition": comp, "valid": valid, "informative": informative})
        features.append(comp[valid])
    pooled = np.vstack(features)
    n_clusters = 5
    model = KMeans(n_clusters=n_clusters, random_state=42, n_init=20).fit(pooled)
    global_rows = []
    pair_tables = []
    for d in data:
        sample = d["sample"]
        niche = np.full(len(d["labels"]), -1, dtype=np.int16)
        niche[d["valid"]] = model.predict(d["composition"][d["valid"]])
        d["niche"] = niche
        sample_out = OUT / "samples" / sample
        sample_out.mkdir(parents=True, exist_ok=True)
        inside = d["inside"]
        table = pd.DataFrame({"sample_id": sample, "x_hires": d["coords"][:, 0], "y_hires": d["coords"][:, 1],
                              "cell_type_provisional": d["labels"], "niche": niche, "informative_neighbors": d["composition"].sum(axis=1) * 12,
                              "inside_he": inside})
        table.to_csv(sample_out / "niche_cells.csv.gz", index=False, compression="gzip")
        comp = pd.DataFrame(d["composition"], columns=TYPES)
        comp["niche"] = niche
        means = comp[comp.niche >= 0].groupby("niche")[TYPES].mean()
        means.to_csv(sample_out / "niche_composition.csv")
        pairs = edge_enrichment(d["labels"], d["coords"])
        pairs.to_csv(sample_out / "cell_type_neighbor_log2_enrichment.csv")
        pair_tables.append(pairs)
        overlay(sample, d["image"], d["coords"][inside], niche[inside], sample_out / f"{sample}_he_niche_overlay.png")
        for n, row in means.iterrows():
            dominant = row.sort_values(ascending=False).head(3)
            global_rows.append({"sample_id": sample, "niche": int(n + 1), "n_cells": int((niche == n).sum()),
                                "dominant_types": "; ".join(f"{TYPE_DISPLAY[k]} ({v:.2f})" for k, v in dominant.items())})
        (sample_out / "niche_summary.json").write_text(json.dumps({"sample_id": sample, "n_cells": len(niche),
            "niche_clusters": n_clusters, "valid_niche_cells": int((niche >= 0).sum()),
            "low_signal_niche_cells": int((niche < 0).sum()), "inside_he_cells": int(inside.sum()),
            "outside_he_cells": int((~inside).sum()), "k_neighbors": 12}, indent=2) + "\n")
    summary = pd.DataFrame(global_rows).sort_values(["sample_id", "niche"])
    summary.to_csv(OUT / "r1_r7_niche_summary.csv", index=False)
    # HTML intentionally links to files so large PNGs remain inspectable locally.
    html = ["<!doctype html><html><head><meta charset='utf-8'><title>R1/R7 Spatial Niche Review</title>",
            "<style>body{font-family:Arial,sans-serif;max-width:1200px;margin:2rem auto} img{max-width:100%;border:1px solid #ccc} table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:5px}</style></head><body>",
            "<h1>R1/R7 前列腺癌空间 Niche 探索性分析</h1>",
            "<p>R1 为未复发原发灶，R7 为复发原发灶。每组目前只有一个患者，以下结果用于空间模式审阅，不作患者层面统计推断。</p>",
            "<h2>Niche 组成</h2>", summary.to_html(index=False), "<h2>H&amp;E Overlay</h2>"]
    for sample in SAMPLES:
        rel = f"samples/{sample}/{sample}_he_niche_overlay.png"
        html += [f"<h3>{sample}</h3><p><img src='{rel}' alt='{sample} niche overlay'></p>"]
    html += ["<h2>解释限制</h2><ul><li>Niche 是基于局部空间细胞类型组成的探索性状态，不等同于病理诊断。</li><li>低信号或未确认细胞类型未被强行解释。</li><li>候选恶性上皮仍需结合 CNV、AMACR/ERG/PCA3 和 H&amp;E 复核。</li></ul></body></html>"]
    (OUT / "r1_r7_spatial_niche_report.html").write_text("\n".join(html), encoding="utf-8")
    print(json.dumps({"status": "pass", "output": str(OUT), "report": str(OUT / "r1_r7_spatial_niche_report.html"), "samples": SAMPLES}, indent=2))


if __name__ == "__main__":
    main()
