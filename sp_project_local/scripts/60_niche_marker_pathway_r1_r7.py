#!/usr/bin/env python3
"""Step 2: molecular characterization of fixed k=12 R1/R7 niches."""
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
import scanpy as sc


PROJECT = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
RAW = Path(os.environ.get("SP_RAW", "/sc/arion/projects/DiseaseGeneCell/Huang_lab_data/SpatialTranscriptome"))
ROOT = PROJECT / "results" / "segmented_official_v1"
NICHE_ROOT = ROOT / "spatial_niche_r1_r7"
OUT = ROOT / "niche_marker_pathway_r1_r7"
SAMPLES = ["SC000895-R1", "SC000895-R7"]
NICHES = {
    0: "fibroblast-luminal epithelial interface",
    1: "contractile smooth muscle-pericyte stroma",
    2: "basal epithelial-contractile stroma interface",
    3: "luminal epithelial-enriched niche",
    4: "fibroblast-rich stroma",
}
GENESETS = {
    "AR_luminal": ["AR", "KLK3", "ACPP", "NKX3-1", "KRT8", "KRT18", "KRT19", "EPCAM"],
    "basal_epithelial": ["KRT5", "KRT14", "TP63", "KRT15", "KRT17", "KRT23"],
    "malignant_candidate": ["AMACR", "ERG", "PCA3", "MYC", "MKI67", "TOP2A", "EZH2", "TPD52"],
    "fibroblast_stroma": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "CFD", "C7"],
    "contractile_perivascular": ["ACTA2", "TAGLN", "MYH11", "RGS5", "CSPG4", "MCAM"],
    "immune": ["PTPRC", "CD3D", "CD3E", "TRBC1", "MS4A1", "CD79A", "LYZ", "TYROBP", "FCER1G", "CD68"],
    "endothelial": ["PECAM1", "VWF", "KDR", "EMCN", "ESAM", "RAMP2", "PLVAP"],
    "proliferation": ["MKI67", "TOP2A", "PCNA", "TYMS", "UBE2C", "MCM2", "MCM6"],
    "inflammatory": ["HLA-DRA", "CD74", "IL7R", "NKG7", "CCL2", "CCL5", "CXCL9", "CXCL10"],
}


def input_path(sample: str) -> Path:
    return ROOT / "analysis_coarse" / "samples" / sample / f"{sample}_official_cell_level_analysis_coarse.h5ad"


def sparse_mean_and_detect(x):
    if sparse.issparse(x):
        mean = np.asarray(x.mean(axis=0)).ravel()
        detect = np.asarray((x > 0).mean(axis=0)).ravel()
    else:
        mean = np.asarray(x).mean(axis=0)
        detect = (np.asarray(x) > 0).mean(axis=0)
    return mean.astype(float), detect.astype(float)


def score_gene_set(x, genes, gene_index):
    present = [gene_index[g.upper()] for g in genes if g.upper() in gene_index]
    if not present:
        return np.zeros(x.shape[0], dtype=np.float32), []
    sub = x[:, present]
    if sparse.issparse(sub):
        score = np.asarray(sub.mean(axis=1)).ravel()
    else:
        score = np.asarray(sub).mean(axis=1)
    return score.astype(np.float32), [str(gene_index[g.upper()][1]) if isinstance(gene_index[g.upper()], tuple) else str(g) for g in genes if g.upper() in gene_index]


def make_spatial_overlay(image, coords, values, title, output, cmap="viridis"):
    fig, ax = plt.subplots(figsize=(11, 11))
    ax.imshow(image)
    finite = np.isfinite(values)
    v = values[finite]
    lo, hi = np.percentile(v, [1, 99]) if len(v) else (0, 1)
    if hi <= lo:
        hi = lo + 1e-6
    sca = ax.scatter(coords[finite, 0], coords[finite, 1], c=values[finite], s=0.8,
                     alpha=0.65, linewidths=0, cmap=cmap, vmin=lo, vmax=hi, rasterized=True)
    ax.set_xlim(0, image.shape[1]); ax.set_ylim(image.shape[0], 0); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title)
    fig.colorbar(sca, ax=ax, fraction=0.035, pad=0.01, label="log1p expression score")
    fig.tight_layout(); fig.savefig(output, dpi=200, bbox_inches="tight"); plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    marker_rows = []
    pathway_rows = []
    report_sections = []
    for sample in SAMPLES:
        niche_dir = NICHE_ROOT / "samples" / sample
        cell_table = pd.read_csv(niche_dir / "niche_cells.csv.gz")
        niche = cell_table["niche"].to_numpy(dtype=int)
        adata = sc.read_h5ad(input_path(sample))
        if len(niche) != adata.n_obs:
            raise ValueError(f"{sample}: niche cells ({len(niche)}) != H5AD cells ({adata.n_obs})")
        x = adata.X.tocsr() if sparse.issparse(adata.X) else np.asarray(adata.X)
        genes = [str(g) for g in adata.var_names]
        gene_index = {g.upper(): (i, g) for i, g in enumerate(genes)}
        sample_out = OUT / "samples" / sample
        sample_out.mkdir(parents=True, exist_ok=True)
        coords = cell_table[["x_hires", "y_hires"]].to_numpy(dtype=float)
        image_path = RAW / "wangy33.u.hpc.mssm.edu" / "10X_Single_Cell_RNA" / "SC000895_Kuan_lin_Huang_T643" / sample / "outs" / "segmented_outputs" / "spatial" / "tissue_hires_image.png"
        image = np.asarray(Image.open(image_path).convert("RGB"))
        inside = ((coords[:, 0] >= 0) & (coords[:, 0] < image.shape[1]) & (coords[:, 1] >= 0) & (coords[:, 1] < image.shape[0])) & (niche >= 0)
        # Descriptive marker statistics: mean and detection within niche vs all other valid niches.
        for n, nname in NICHES.items():
            mask = niche == n
            other = (niche >= 0) & ~mask
            if not mask.any():
                continue
            mean_n, det_n = sparse_mean_and_detect(x[mask])
            mean_o, det_o = sparse_mean_and_detect(x[other]) if other.any() else (np.zeros_like(mean_n), np.zeros_like(det_n))
            logfc = mean_n - mean_o
            top = np.argsort(logfc)[::-1]
            kept = 0
            for idx in top:
                if kept >= 40 or logfc[idx] <= 0:
                    break
                marker_rows.append({"sample_id": sample, "niche": n + 1, "niche_name": nname,
                                    "gene": genes[idx], "mean_in_niche": mean_n[idx], "mean_outside": mean_o[idx],
                                    "log1p_difference": logfc[idx], "detection_in_niche": det_n[idx],
                                    "detection_outside": det_o[idx]})
                kept += 1
        # Pathway scores and per-niche summaries.
        score_frame = pd.DataFrame(index=np.arange(adata.n_obs))
        present_sets = {}
        for name, gene_set in GENESETS.items():
            present = [gene_index[g.upper()][0] for g in gene_set if g.upper() in gene_index]
            present_sets[name] = [genes[i] for i in present]
            if present:
                sub = x[:, present]
                score_frame[name] = np.asarray(sub.mean(axis=1)).ravel() if sparse.issparse(sub) else np.asarray(sub).mean(axis=1)
            else:
                score_frame[name] = 0.0
            for n, nname in NICHES.items():
                mask = niche == n
                if mask.any():
                    pathway_rows.append({"sample_id": sample, "niche": n + 1, "niche_name": nname,
                                         "pathway": name, "genes_present": ";".join(present_sets[name]),
                                         "mean_score": float(score_frame.loc[mask, name].mean()),
                                         "median_score": float(score_frame.loc[mask, name].median())})
        score_frame.insert(0, "sample_id", sample); score_frame.insert(1, "niche", niche)
        score_frame.to_csv(sample_out / "niche_cell_pathway_scores.csv.gz", index=False, compression="gzip")
        for name in GENESETS:
            if inside.any():
                make_spatial_overlay(image, coords[inside], score_frame.loc[inside, name].to_numpy(dtype=float),
                                     f"{sample}: {name} score", sample_out / f"{sample}_spatial_{name}.png")
        means = pd.DataFrame(pathway_rows).query("sample_id == @sample").pivot_table(index="niche_name", columns="pathway", values="mean_score")
        means.to_csv(sample_out / "niche_pathway_mean_matrix.csv")
        report_sections.append((sample, sample_out, present_sets))
    marker_df = pd.DataFrame(marker_rows)
    marker_df.to_csv(OUT / "niche_top_markers.csv", index=False)
    pathway_df = pd.DataFrame(pathway_rows)
    pathway_df.to_csv(OUT / "niche_pathway_scores.csv", index=False)
    html = ["<!doctype html><html><head><meta charset='utf-8'><title>R1/R7 Niche Marker Pathway</title>",
            "<style>body{font-family:Arial,sans-serif;max-width:1250px;margin:2rem auto}img{max-width:100%;border:1px solid #ccc}table{border-collapse:collapse;font-size:12px}td,th{border:1px solid #ccc;padding:4px;vertical-align:top}</style></head><body>",
            "<h1>R1/R7 Niche Marker 与通路特征</h1>",
            "<p>使用固定 k=12 Niche。marker 使用 Niche 内与其他 Niche 的 log1p 平均表达差异排序；通路为基因集平均表达分数。结果为描述性分析，未使用细胞级 p 值替代患者重复。</p>",
            "<h2>Pathway scores</h2>", pathway_df.to_html(index=False), "<h2>Top markers</h2>", marker_df.groupby(["sample_id", "niche", "niche_name"]).head(10).to_html(index=False)]
    for sample, sample_out, present_sets in report_sections:
        html.append(f"<h2>{sample}</h2><p>Genes present: {json.dumps(present_sets, ensure_ascii=False)}</p>")
        for name in GENESETS:
            rel = f"samples/{sample}/{sample}_spatial_{name}.png"
            html.append(f"<h3>{name}</h3><img src='{rel}' alt='{sample} {name}'>")
    html += ["<h2>Interpretation limits</h2><ul><li>候选恶性上皮仍不是最终癌细胞标签，需要 CNV、EPCAM/KRT 和 H&amp;E 共同验证。</li><li>R1/R7 各只有一个患者，不能进行正式复发组推断。</li><li>pathway score 是表达描述，不等于通路活性因果证据。</li></ul></body></html>"]
    (OUT / "r1_r7_niche_marker_pathway_report.html").write_text("\n".join(html), encoding="utf-8")
    (OUT / "niche_marker_pathway_summary.json").write_text(json.dumps({"status": "pass", "samples": SAMPLES, "niches": NICHES, "genesets": GENESETS}, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"status": "pass", "report": str(OUT / "r1_r7_niche_marker_pathway_report.html")}, indent=2))


if __name__ == "__main__":
    main()
