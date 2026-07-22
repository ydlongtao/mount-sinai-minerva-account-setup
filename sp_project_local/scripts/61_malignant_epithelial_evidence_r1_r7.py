#!/usr/bin/env python3
"""Step 3: evidence-based screening of candidate malignant epithelium."""
from __future__ import annotations

import importlib.util
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
import scanpy as sc


PROJECT = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
RAW = Path(os.environ.get("SP_RAW", "/sc/arion/projects/DiseaseGeneCell/Huang_lab_data/SpatialTranscriptome"))
ROOT = PROJECT / "results" / "segmented_official_v1"
NICHE_ROOT = ROOT / "spatial_niche_r1_r7"
OUT = ROOT / "malignant_epithelial_evidence_r1_r7"
SAMPLES = ["SC000895-R1", "SC000895-R7"]
PANELS = {
    "epithelial_context": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT14", "TP63"],
    "luminal_context": ["AR", "KLK3", "ACPP", "NKX3-1", "KRT8", "KRT18", "KRT19"],
    "malignant_candidate": ["AMACR", "ERG", "PCA3", "MYC", "EZH2", "TPD52"],
    "proliferation": ["MKI67", "TOP2A", "PCNA", "TYMS", "UBE2C", "MCM2", "MCM6"],
    "immune_context": ["PTPRC", "CD3D", "CD3E", "MS4A1", "CD79A", "LYZ", "TYROBP", "FCER1G", "CD68"],
    "stromal_context": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "ACTA2", "TAGLN"],
}


def robust_z(values: np.ndarray) -> np.ndarray:
    med = float(np.nanmedian(values))
    mad = float(np.nanmedian(np.abs(values - med)))
    scale = 1.4826 * mad
    if scale < 1e-8:
        scale = float(np.nanstd(values)) or 1.0
    return (values - med) / scale


def input_path(sample: str) -> Path:
    return ROOT / "analysis_coarse" / "samples" / sample / f"{sample}_official_cell_level_analysis_coarse.h5ad"


def score_panel(x, genes, index):
    present = [index[g.upper()][0] for g in genes if g.upper() in index]
    if not present:
        return np.zeros(x.shape[0], dtype=np.float32), []
    sub = x[:, present]
    score = np.asarray(sub.mean(axis=1)).ravel() if sparse.issparse(sub) else np.asarray(sub).mean(axis=1)
    return score.astype(np.float32), [index[g.upper()][1] for g in genes if g.upper() in index]


def overlay(image, coords, evidence, output, sample):
    colors = {"high_candidate": "#e31a1c", "moderate_candidate": "#ff7f00",
              "epithelial_reference": "#1f78b4", "low_or_non_epithelial": "#bdbdbd"}
    names = {"high_candidate": "high candidate", "moderate_candidate": "moderate candidate",
             "epithelial_reference": "epithelial reference", "low_or_non_epithelial": "low/non-epithelial"}
    fig, ax = plt.subplots(figsize=(12, 12)); ax.imshow(image)
    for group in colors:
        mask = evidence == group
        if mask.any():
            ax.scatter(coords[mask, 0], coords[mask, 1], s=1.3 if group.startswith("high") else 0.9,
                       c=colors[group], alpha=0.65, linewidths=0, rasterized=True, label=names[group])
    ax.set_xlim(0, image.shape[1]); ax.set_ylim(image.shape[0], 0); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"{sample}: H&E + malignant epithelial evidence tiers")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), markerscale=4)
    fig.tight_layout(); fig.savefig(output, dpi=220, bbox_inches="tight"); plt.close(fig)


def cnv_status() -> dict:
    result = {"infercnvpy_importable": False, "note": "CNV not run in this step; requires genomic positions and a reference population."}
    if importlib.util.find_spec("infercnvpy") is not None:
        result["infercnvpy_importable"] = True
        result["note"] = "infercnvpy is importable, but CNV was not run without a validated gene-position reference and reference-cell definition."
    return result


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cluster_rows, cell_rows, sample_sections = [], [], []
    for sample in SAMPLES:
        adata = sc.read_h5ad(input_path(sample))
        niche_table = pd.read_csv(NICHE_ROOT / "samples" / sample / "niche_cells.csv.gz")
        marker_table = pd.read_csv(ROOT / "preannotation_prostate" / "samples" / sample / "cluster_prostate_marker_spatial_preannotation.csv", dtype={"cluster": str})
        if len(niche_table) != adata.n_obs:
            raise ValueError(f"{sample}: niche/H5AD row mismatch")
        label_key = "cell_level_leiden_coarse"
        clusters = adata.obs[label_key].astype(str).to_numpy()
        x = adata.X.tocsr() if sparse.issparse(adata.X) else np.asarray(adata.X)
        genes = [str(g) for g in adata.var_names]
        index = {g.upper(): (i, g) for i, g in enumerate(genes)}
        scores = {}
        present = {}
        for name, panel in PANELS.items():
            scores[name], present[name] = score_panel(x, panel, index)
        z_ep = robust_z(scores["epithelial_context"])
        z_mal = robust_z(scores["malignant_candidate"])
        z_pro = robust_z(scores["proliferation"])
        z_lum = robust_z(scores["luminal_context"])
        # Use cluster-level evidence to avoid amplifying noisy cell-level deviations.
        marker_map = marker_table.assign(cluster=marker_table["cluster"].astype(str)).set_index("cluster")["predicted_type"].to_dict()
        score_table = pd.DataFrame({"cluster": clusters, "epithelial": scores["epithelial_context"],
                                    "malignant": scores["malignant_candidate"]}).groupby("cluster").mean()
        malignant_threshold = float(score_table["malignant"].quantile(0.75))
        cluster_evidence = {}
        for cluster, row in score_table.iterrows():
            marker_label = marker_map.get(str(cluster), "ambiguous_or_low_signal")
            if marker_label == "malignant_epithelial_candidate" and row["epithelial"] > 0.03:
                cluster_evidence[str(cluster)] = "high_candidate"
            elif marker_label in {"malignant_epithelial_candidate", "luminal_epithelial"} and row["malignant"] >= malignant_threshold and row["epithelial"] > 0.03:
                cluster_evidence[str(cluster)] = "moderate_candidate"
            elif marker_label in {"luminal_epithelial", "basal_epithelial"}:
                cluster_evidence[str(cluster)] = "epithelial_reference"
            else:
                cluster_evidence[str(cluster)] = "low_or_non_epithelial"
        evidence = np.array([cluster_evidence.get(str(c), "low_or_non_epithelial") for c in clusters], dtype=object)
        out = OUT / "samples" / sample; out.mkdir(parents=True, exist_ok=True)
        coords = niche_table[["x_hires", "y_hires"]].to_numpy(float)
        image_path = RAW / "wangy33.u.hpc.mssm.edu" / "10X_Single_Cell_RNA" / "SC000895_Kuan_lin_Huang_T643" / sample / "outs" / "segmented_outputs" / "spatial" / "tissue_hires_image.png"
        image = np.asarray(Image.open(image_path).convert("RGB"))
        inside = ((coords[:, 0] >= 0) & (coords[:, 0] < image.shape[1]) & (coords[:, 1] >= 0) & (coords[:, 1] < image.shape[0]))
        overlay(image, coords[inside], evidence[inside], out / f"{sample}_he_malignant_evidence_overlay.png", sample)
        cells = pd.DataFrame({"sample_id": sample, "cluster": clusters, "niche": niche_table["niche"].to_numpy(),
                              "x_hires": coords[:, 0], "y_hires": coords[:, 1], **{f"score_{k}": v for k, v in scores.items()},
                              "z_epithelial": z_ep, "z_malignant_candidate": z_mal, "z_proliferation": z_pro,
                              "z_luminal": z_lum, "evidence_tier": evidence})
        cells.to_csv(out / "malignant_evidence_cells.csv.gz", index=False, compression="gzip")
        for cluster, group in cells.groupby("cluster", sort=False):
            m = marker_table[marker_table["cluster"].astype(str) == str(cluster)]
            tier_counts = group["evidence_tier"].value_counts().to_dict()
            cluster_rows.append({"sample_id": sample, "cluster": str(cluster), "n_cells": len(group),
                                 "provisional_marker_label": m["predicted_type"].iloc[0] if len(m) else "unknown",
                                 "mean_epithelial": group.score_epithelial_context.mean(),
                                 "mean_luminal": group.score_luminal_context.mean(),
                                 "mean_malignant_candidate": group.score_malignant_candidate.mean(),
                                 "mean_proliferation": group.score_proliferation.mean(),
                                 "evidence_tier": cluster_evidence.get(str(cluster), "low_or_non_epithelial"),
                                 "high_candidate_cells": int(tier_counts.get("high_candidate", 0)),
                                 "moderate_candidate_cells": int(tier_counts.get("moderate_candidate", 0)),
                                 "evidence_counts": json.dumps(tier_counts)})
        sample_sections.append((sample, out, present))
    cluster_df = pd.DataFrame(cluster_rows)
    cluster_df.to_csv(OUT / "malignant_evidence_cluster_summary.csv", index=False)
    (OUT / "cnv_status.json").write_text(json.dumps(cnv_status(), indent=2) + "\n")
    html = ["<!doctype html><html><head><meta charset='utf-8'><title>R1/R7 Malignant Evidence</title>",
            "<style>body{font-family:Arial,sans-serif;max-width:1250px;margin:2rem auto}img{max-width:100%;border:1px solid #ccc}table{border-collapse:collapse;font-size:12px}td,th{border:1px solid #ccc;padding:4px;vertical-align:top}</style></head><body>",
            "<h1>R1/R7 候选恶性上皮证据报告</h1>",
            "<p>证据分层整合上皮背景、前列腺腔面、AMACR/ERG/PCA3 等候选恶性 marker、增殖信号、空间 Niche 和 H&amp;E。分层不是癌症诊断。</p>",
            "<h2>Cluster summary</h2>", cluster_df.to_html(index=False)]
    for sample, out, present in sample_sections:
        rel = f"samples/{sample}/{sample}_he_malignant_evidence_overlay.png"
        html += [f"<h2>{sample}</h2><p>Panel genes present: {json.dumps(present, ensure_ascii=False)}</p><img src='{rel}' alt='{sample} malignant evidence overlay'>"]
    html += ["<h2>Next confirmation</h2><ol><li>优先复核 high_candidate cluster 的 H&amp;E 腺体形态。</li><li>检查 EPCAM/KRT 与 AMACR/ERG/PCA3 是否共同出现。</li><li>补充可靠的基因组位置和正常参考细胞后，再运行 infercnvpy/CNV 分析。</li><li>只有多证据一致时，才将候选标签升级为 tumor epithelial。</li></ol></body></html>"]
    (OUT / "r1_r7_malignant_epithelial_evidence_report.html").write_text("\n".join(html), encoding="utf-8")
    (OUT / "malignant_evidence_summary.json").write_text(json.dumps({"status": "pass", "samples": SAMPLES, "panels": PANELS, "cnv": cnv_status()}, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"status": "pass", "report": str(OUT / "r1_r7_malignant_epithelial_evidence_report.html")}, indent=2))


if __name__ == "__main__":
    main()
