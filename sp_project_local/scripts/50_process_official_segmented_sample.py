#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd
from PIL import Image
import scanpy as sc


PROJECT_HOME = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
MANIFEST = Path(os.environ.get("SP_PROJECT_MANIFEST", PROJECT_HOME / "config" / "sample_manifest_v2_corrected.csv"))
RESULTS = Path(os.environ.get("SP_OFFICIAL_RESULTS", PROJECT_HOME / "results" / "segmented_official_v1"))

MARKERS = {
    "luminal": ["KLK3", "ACPP", "AR", "NKX3-1", "KRT8", "KRT18"],
    "basal": ["KRT5", "KRT14", "TP63", "KRT15"],
    "tumor": ["AMACR", "ERG", "MYC", "MKI67", "TOP2A"],
    "immune": ["PTPRC", "CD3D", "CD3E", "MS4A1", "CD68", "LYZ"],
    "stromal": ["COL1A1", "COL1A2", "DCN", "LUM", "ACTA2"],
    "endothelial": ["PECAM1", "VWF", "KDR", "EMCN"],
}


def manifest_row(index: int) -> dict[str, str]:
    rows = list(csv.DictReader(MANIFEST.open(newline="")))
    if not 1 <= index <= len(rows):
        raise IndexError(index)
    return rows[index - 1]


def polygon_stats(ring: np.ndarray) -> tuple[float, float, float]:
    if ring.shape[0] < 3:
        return float(ring[:, 0].mean()), float(ring[:, 1].mean()), 0.0
    if np.allclose(ring[0], ring[-1]):
        ring = ring[:-1]
    x, y = ring[:, 0], ring[:, 1]
    cross = x * np.roll(y, -1) - np.roll(x, -1) * y
    signed = cross.sum() / 2.0
    if abs(signed) < 1e-8:
        return float(x.mean()), float(y.mean()), 0.0
    cx = ((x + np.roll(x, -1)) * cross).sum() / (6.0 * signed)
    cy = ((y + np.roll(y, -1)) * cross).sum() / (6.0 * signed)
    return float(cx), float(cy), float(abs(signed))


def read_cell_geojson(path: Path) -> tuple[pd.DataFrame, list[np.ndarray]]:
    # Space Ranger GeoJSON is a single FeatureCollection. Loading one sample at
    # a time is bounded and keeps the protected source directory read-only.
    data = json.loads(path.read_text())
    records = []
    rings = []
    for feature in data.get("features", []):
        geometry = feature.get("geometry") or {}
        if geometry.get("type") != "Polygon" or not geometry.get("coordinates"):
            continue
        ring = np.asarray(geometry["coordinates"][0], dtype=np.float32)
        if ring.ndim != 2 or ring.shape[1] != 2 or ring.shape[0] < 3:
            continue
        numeric_id = int(feature["properties"]["cell_id"])
        cx, cy, area = polygon_stats(ring)
        records.append((f"cellid_{numeric_id:09d}-1", cx, cy, area))
        rings.append(ring)
    frame = pd.DataFrame(records, columns=["cell_id", "x_fullres", "y_fullres", "area_fullres_px2"]).set_index("cell_id")
    if frame.index.has_duplicates:
        raise ValueError(f"Duplicate cell IDs in {path}")
    return frame, rings


def marker_scores(adata) -> dict[str, list[str]]:
    upper = {str(g).upper(): g for g in adata.var_names}
    totals = np.asarray(adata.X.sum(axis=1)).ravel().astype(np.float64)
    scale = np.divide(1e4, totals, out=np.zeros_like(totals), where=totals > 0)
    present = {}
    for group, genes in MARKERS.items():
        found = [upper[g] for g in genes if g in upper]
        present[group] = [str(g) for g in found]
        if found:
            values = np.asarray(adata[:, found].X.sum(axis=1)).ravel()
            adata.obs[f"score_{group}"] = np.log1p(values * scale).astype(np.float32)
    return present


def overview_overlay(image_path: Path, rings: list[np.ndarray], scalef: float, output: Path) -> None:
    image = np.asarray(Image.open(image_path).convert("RGB"))
    rng = np.random.default_rng(0)
    keep = np.arange(len(rings))
    if len(keep) > 30000:
        keep = np.sort(rng.choice(keep, 30000, replace=False))
    lines = [rings[i] * scalef for i in keep]
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(image)
    ax.add_collection(LineCollection(lines, colors="#00e5ff", linewidths=0.18, alpha=0.55, rasterized=True))
    ax.set_xlim(0, image.shape[1]); ax.set_ylim(image.shape[0], 0)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"Official Space Ranger cell boundaries (sampled {len(lines):,}/{len(rings):,})")
    fig.tight_layout(); fig.savefig(output, dpi=220, bbox_inches="tight"); plt.close(fig)


def qc_figure(adata, output: Path, sample: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, col, title in [
        (axes[0, 0], "total_counts", "Total counts"),
        (axes[0, 1], "n_genes_by_counts", "Detected genes"),
        (axes[1, 0], "pct_counts_mt", "Mitochondrial fraction (%)"),
        (axes[1, 1], "cell_area_um2", "Cell area (um²)"),
    ]:
        values = adata.obs[col].to_numpy(dtype=float)
        cap = np.nanpercentile(values, 99.5)
        ax.hist(values[np.isfinite(values) & (values <= cap)], bins=80, color="#2563eb", alpha=0.85)
        ax.set_title(title); ax.set_xlabel(col); ax.set_ylabel("Cells")
    fig.suptitle(f"{sample}: official segmented-cell QC")
    fig.tight_layout(); fig.savefig(output, dpi=180, bbox_inches="tight"); plt.close(fig)


def main() -> None:
    index = int(os.environ.get("LSB_JOBINDEX", os.environ.get("SP_SAMPLE_INDEX", "1")))
    row = manifest_row(index)
    sample = row["sample_id"]
    outs = Path(row["outs_dir"])
    segmented = outs / "segmented_outputs"
    out = RESULTS / "samples" / sample
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    adata = sc.read_10x_h5(row["segmented_h5"])
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    if adata.X.nnz and not np.allclose(adata.X.data, np.rint(adata.X.data)):
        raise ValueError("Official feature-cell matrix does not contain integer counts")
    adata.layers["counts"] = adata.X.copy()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    polygons, rings = read_cell_geojson(segmented / "cell_segmentations.geojson")
    matched = polygons.reindex(adata.obs_names)
    if matched[["x_fullres", "y_fullres"]].isna().any(axis=None):
        missing = int(matched["x_fullres"].isna().sum())
        raise ValueError(f"{sample}: {missing} matrix cells have no official polygon")
    scalefactors = json.loads((segmented / "spatial" / "scalefactors_json.json").read_text())
    mpp = float(scalefactors["microns_per_pixel"])
    hires_scalef = float(scalefactors["tissue_hires_scalef"])
    adata.obs["x_fullres"] = matched["x_fullres"].to_numpy(dtype=np.float32)
    adata.obs["y_fullres"] = matched["y_fullres"].to_numpy(dtype=np.float32)
    adata.obs["cell_area_um2"] = (matched["area_fullres_px2"].to_numpy() * mpp * mpp).astype(np.float32)
    adata.obsm["spatial"] = matched[["x_fullres", "y_fullres"]].to_numpy(dtype=np.float32)
    adata.obsm["spatial_hires"] = adata.obsm["spatial"] * hires_scalef
    adata.uns["spatial_coordinate_source"] = "Space Ranger cell_segmentations.geojson polygon centroid"
    adata.uns["spatial_coordinate_unit"] = "full-resolution pixels"
    adata.uns["official_segmentation"] = {
        "cell_geojson": str(segmented / "cell_segmentations.geojson"),
        "nucleus_geojson": str(segmented / "nucleus_segmentations.geojson"),
        "barcode_mappings": str(outs / "barcode_mappings.parquet"),
        "microns_per_pixel": mpp,
        "tissue_hires_scalef": hires_scalef,
    }
    present_markers = marker_scores(adata)

    # 保留全部官方 filtered cells，仅添加候选 QC 标记，等待人工复核后再筛选。
    adata.obs["passes_qc_candidate"] = (
        (adata.obs["total_counts"] >= 20)
        & (adata.obs["n_genes_by_counts"] >= 10)
        & (adata.obs["pct_counts_mt"] <= 30)
        & (adata.obs["cell_area_um2"] >= 10)
        & (adata.obs["cell_area_um2"] <= 1000)
    )

    overview_overlay(
        segmented / "spatial" / "tissue_hires_image.png",
        rings,
        hires_scalef,
        figures / "official_cell_boundaries_he_overview.png",
    )
    qc_figure(adata, figures / "official_cell_qc_distributions.png", sample)

    h5ad = out / f"{sample}_official_segmented_cells.h5ad"
    adata.write_h5ad(h5ad, compression="lzf")
    summary = {
        "sample_id": sample,
        "status": "pass",
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "polygon_features": int(len(polygons)),
        "matched_polygon_percent": 100.0,
        "counts_layer": "counts",
        "integer_counts": True,
        "median_total_counts": float(adata.obs["total_counts"].median()),
        "median_genes": float(adata.obs["n_genes_by_counts"].median()),
        "median_pct_mt": float(adata.obs["pct_counts_mt"].median()),
        "median_cell_area_um2": float(adata.obs["cell_area_um2"].median()),
        "passes_qc_candidate": int(adata.obs["passes_qc_candidate"].sum()),
        "passes_qc_candidate_percent": float(100 * adata.obs["passes_qc_candidate"].mean()),
        "markers_present": present_markers,
        "h5ad": str(h5ad),
        "overlay": str(figures / "official_cell_boundaries_he_overview.png"),
    }
    (out / "official_segmented_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    adata.obs.to_csv(out / "official_segmented_cell_qc.csv.gz", compression="gzip")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
