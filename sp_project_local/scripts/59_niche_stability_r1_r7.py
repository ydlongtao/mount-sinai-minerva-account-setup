#!/usr/bin/env python3
"""Step 1: sensitivity analysis for R1/R7 spatial niche definitions."""
from __future__ import annotations

import json
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

_base_path = Path(__file__).with_name("58_spatial_niche_r1_r7.py")
_spec = importlib.util.spec_from_file_location("spatial_niche_base", _base_path)
_base = importlib.util.module_from_spec(_spec)
assert _spec and _spec.loader
_spec.loader.exec_module(_base)
BASE_OUT = _base.OUT
SAMPLES = _base.SAMPLES
load_sample = _base.load_sample
overlay = _base.overlay
spatial_composition = _base.spatial_composition
TYPES = _base.TYPES
TYPE_DISPLAY = _base.TYPE_DISPLAY


OUT = BASE_OUT.parent / "niche_stability_r1_r7"
K_VALUES = [8, 12, 20]
N_NICHES = 5


def align_mapping(reference: np.ndarray, current: np.ndarray, n: int) -> dict[int, int]:
    mask = (reference >= 0) & (current >= 0)
    table = np.zeros((n, n), dtype=np.int64)
    np.add.at(table, (reference[mask], current[mask]), 1)
    rows, cols = linear_sum_assignment(-table)
    return {int(c): int(r) for r, c in zip(rows, cols)}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    data = [load_sample(s) for s in SAMPLES]
    results: dict[int, list[dict]] = {}
    pooled_features: dict[int, np.ndarray] = {}
    pooled_valid: dict[int, np.ndarray] = {}
    for k in K_VALUES:
        chunks = []
        for d in data:
            comp, valid, _ = spatial_composition(d["coords"], d["labels"], k=k)
            d[f"composition_{k}"] = comp
            d[f"valid_{k}"] = valid
            chunks.append(comp[valid])
        pooled_features[k] = np.vstack(chunks)
        pooled_valid[k] = np.concatenate([d[f"valid_{k}"] for d in data])
        model = KMeans(n_clusters=N_NICHES, random_state=42, n_init=20).fit(pooled_features[k])
        offset = 0
        for d in data:
            labels = np.full(len(d["labels"]), -1, dtype=np.int16)
            labels[d[f"valid_{k}"]] = model.predict(d[f"composition_{k}"][d[f"valid_{k}"]])
            d[f"niche_{k}"] = labels
            offset += len(labels)

    reference = {d["sample"]: d["niche_12"] for d in data}
    rows = []
    for k in K_VALUES:
        k_out = OUT / f"k{k}"
        for d in data:
            sample = d["sample"]
            niche = d[f"niche_{k}"]
            valid = niche >= 0
            ari = float(adjusted_rand_score(reference[sample][valid & (reference[sample] >= 0)], niche[valid & (reference[sample] >= 0)])) if k != 12 else 1.0
            comp = pd.DataFrame(d[f"composition_{k}"], columns=TYPES)
            comp["niche"] = niche
            means = comp[comp.niche >= 0].groupby("niche")[TYPES].mean()
            sample_out = k_out / "samples" / sample
            sample_out.mkdir(parents=True, exist_ok=True)
            means.to_csv(sample_out / "niche_composition.csv")
            overlay(sample, d["image"], d["coords"][d["inside"]], niche[d["inside"]], sample_out / f"{sample}_he_niche_overlay_k{k}.png")
            dominant = []
            for n, row in means.iterrows():
                top = row.sort_values(ascending=False).head(3)
                dominant.append({"niche": int(n + 1), "n_cells": int((niche == n).sum()),
                                 "dominant_types": "; ".join(f"{TYPE_DISPLAY[x]} ({v:.2f})" for x, v in top.items())})
            rows.append({"sample_id": sample, "k_neighbors": k, "n_niches": N_NICHES,
                         "valid_niche_cells": int(valid.sum()), "valid_fraction": float(valid.mean()),
                         "ari_vs_k12": ari, "niche_sizes": json.dumps(dominant)})
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "niche_stability_summary.csv", index=False)
    html = ["<!doctype html><html><head><meta charset='utf-8'><title>R1/R7 Niche Stability</title>",
            "<style>body{font-family:Arial,sans-serif;max-width:1250px;margin:2rem auto}img{max-width:100%;border:1px solid #ccc}table{border-collapse:collapse;font-size:13px}td,th{border:1px solid #ccc;padding:5px;vertical-align:top}</style></head><body>",
            "<h1>R1/R7 Niche 稳定性报告</h1>",
            "<p>比较空间邻域 k=8、12、20。Niche 数固定为 5；ARI 以 k=12 作为参考。结果用于判断参数敏感性，不是患者层面统计推断。</p>",
            "<h2>参数稳定性</h2>", summary.drop(columns=["niche_sizes"]).to_html(index=False),
            "<h2>各参数 H&amp;E Overlay</h2>"]
    for sample in SAMPLES:
        html.append(f"<h3>{sample}</h3>")
        for k in K_VALUES:
            rel = f"k{k}/samples/{sample}/{sample}_he_niche_overlay_k{k}.png"
            html.append(f"<h4>k={k}</h4><img src='{rel}' alt='{sample} k={k} niche overlay'>")
        sample_rows = summary[summary.sample_id == sample][["k_neighbors", "niche_sizes"]]
        html.append(sample_rows.to_html(index=False))
    html += ["<h2>判定建议</h2><ul><li>优先选择在 k=8/12/20 下空间边界和主要组成均稳定的 Niche。</li><li>若 ARI 较低但 H&amp;E 结构一致，应优先依据病理边界和组成解释，而不是单独依据 ARI。</li><li>R1/R7 各只有一个患者，不能据此进行复发组统计推断。</li></ul></body></html>"]
    (OUT / "r1_r7_niche_stability_report.html").write_text("\n".join(html), encoding="utf-8")
    (OUT / "niche_stability_summary.json").write_text(json.dumps({"status": "pass", "k_values": K_VALUES, "n_niches": N_NICHES}, indent=2) + "\n")
    print(json.dumps({"status": "pass", "output": str(OUT), "report": str(OUT / "r1_r7_niche_stability_report.html")}, indent=2))


if __name__ == "__main__":
    main()
