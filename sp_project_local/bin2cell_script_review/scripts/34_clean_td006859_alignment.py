#!/usr/bin/env python3
"""Make uncluttered alignment figures for the two TD006859 samples."""
from __future__ import annotations

import html
import json
import re
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from PIL import Image
from scipy.spatial import cKDTree

from batch_utils import BATCH_RESULTS, MANIFEST, read_manifest, write_json

TARGETS = {"TD006859-B408", "TD006859-B573"}
BARCODE = re.compile(r"^s_(\d+)um_(\d+)_(\d+)-\d+$")
CROP_X = (2000.0, 3000.0)
CROP_Y = (2000.0, 3000.0)


def decode(values):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in values]


def raw_2um(path: Path):
    with h5py.File(path, "r") as h5:
        matrix = h5["matrix"]
        barcodes = decode(matrix["barcodes"][:])
        indptr = np.asarray(matrix["indptr"][:], dtype=np.int64)
        data = matrix["data"]
        counts = np.zeros(len(barcodes), dtype=np.float64)
        for start in range(0, len(barcodes), 50000):
            stop = min(start + 50000, len(barcodes))
            ds, de = int(indptr[start]), int(indptr[stop])
            vals = np.asarray(data[ds:de], dtype=np.float64)
            offsets = indptr[start:stop + 1] - ds
            cumulative = np.concatenate(([0.0], np.cumsum(vals)))
            counts[start:stop] = cumulative[offsets[1:]] - cumulative[offsets[:-1]]
    matches = [BARCODE.match(x) for x in barcodes]
    size = int(matches[0].group(1))
    x_um = np.asarray([int(m.group(3)) * size + size / 2 for m in matches], dtype=float)
    y_um = np.asarray([int(m.group(2)) * size + size / 2 for m in matches], dtype=float)
    return np.asarray(barcodes), x_um, y_um, counts


def get_crop(image_path: Path, x_um, y_um):
    image = np.asarray(Image.open(image_path).convert("RGB"))
    px_per_um = min((image.shape[1] - 1) / max(x_um), (image.shape[0] - 1) / max(y_um))
    x0, y0 = int(np.floor(CROP_X[0] * px_per_um)), int(np.floor(CROP_Y[0] * px_per_um))
    x1, y1 = min(image.shape[1], int(np.ceil(CROP_X[1] * px_per_um))), min(image.shape[0], int(np.ceil(CROP_Y[1] * px_per_um)))
    return image[y0:y1, x0:x1], px_per_um, x0, y0


def plot(image, output, title, points=None, colors=None, sizes=3, alpha=0.8):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(image)
    if points is not None and len(points):
        if colors is None:
            colors = "#00e5ff"
        ax.scatter(points[:, 0], points[:, 1], s=sizes, c=colors, alpha=alpha, linewidths=0, rasterized=True)
    ax.set_title(title)
    ax.set_xlim(0, image.shape[1]); ax.set_ylim(image.shape[0], 0); ax.axis("off")
    fig.tight_layout(); fig.savefig(output, dpi=240, bbox_inches="tight"); plt.close(fig)


def main():
    out = BATCH_RESULTS / "debug_td006859" / "alignment_clean"
    out.mkdir(parents=True, exist_ok=True)
    summaries = []
    for row in read_manifest(MANIFEST):
        sample = row["sample_id"]
        if sample not in TARGETS:
            continue
        sample_out = out / sample; sample_out.mkdir(parents=True, exist_ok=True)
        barcodes, x_um, y_um, counts = raw_2um(Path(row["matrix_002um"]))
        in_crop = (x_um >= CROP_X[0]) & (x_um < CROP_X[1]) & (y_um >= CROP_Y[0]) & (y_um < CROP_Y[1])
        image, px, x0, y0 = get_crop(Path(row["outs_dir"]) / "spatial" / "tissue_hires_image.png", x_um, y_um)
        crop_coords = np.column_stack((x_um[in_crop] * px - x0, y_um[in_crop] * px - y0))
        crop_counts = counts[in_crop]
        threshold = float(np.percentile(crop_counts[crop_counts > 0], 99.5))
        high = crop_counts >= threshold
        high_points = crop_coords[high]
        if len(high_points) > 8000:
            high_points = high_points[np.linspace(0, len(high_points) - 1, 8000).astype(int)]

        h5ad = BATCH_RESULTS / "samples" / sample / "cellpose" / f"{sample}_002um_cellpose_he.h5ad"
        adata = sc.read_h5ad(h5ad, backed="r")
        coords = np.asarray(adata.obsm["spatial"])
        labels = np.asarray(adata.obs["labels_he"]) > 0
        label_points = coords[labels]
        names = np.asarray(adata.obs_names).astype(str)
        high_names = set(barcodes[in_crop][crop_counts >= threshold])
        label_names = set(names[labels])
        overlap = len(high_names & label_names)
        transform_scores = {}
        width, height = image.shape[1], image.shape[0]
        transforms = {
            "identity": label_points,
            "flip_x": np.column_stack((width - label_points[:, 0], label_points[:, 1])),
            "flip_y": np.column_stack((label_points[:, 0], height - label_points[:, 1])),
            "flip_xy": np.column_stack((width - label_points[:, 0], height - label_points[:, 1])),
            "swap_xy": label_points[:, ::-1],
        }
        if len(label_points) and len(high_points):
            for name, transformed in transforms.items():
                nearest = cKDTree(transformed).query(high_points, k=1)[0]
                transform_scores[name] = {"median_px": float(np.median(nearest)), "p90_px": float(np.percentile(nearest, 90)), "within_10px": float(np.mean(nearest <= 10)), "within_25px": float(np.mean(nearest <= 25))}
        plot(image, sample_out / "00_he_only.png", f"{sample}: H&E only")
        plot(image, sample_out / "01_he_highcount_bins.png", f"{sample}: top 0.5% count bins", high_points, "#00e5ff", 2.5, 0.75)
        plot(image, sample_out / "02_he_labels.png", f"{sample}: Cellpose labels", label_points, "#ff3b30", 4, 0.9)
        if len(high_points) and len(label_points):
            fig, ax = plt.subplots(figsize=(8, 8)); ax.imshow(image)
            ax.scatter(high_points[:, 0], high_points[:, 1], s=2.5, c="#00e5ff", alpha=0.7, label="top 0.5% counts", linewidths=0, rasterized=True)
            ax.scatter(label_points[:, 0], label_points[:, 1], s=5, facecolors="none", edgecolors="#ff3b30", alpha=0.95, label="Cellpose labels", linewidths=0.8, rasterized=True)
            ax.legend(loc="upper right"); ax.set_title(f"{sample}: high-count bins vs Cellpose"); ax.set_xlim(0, image.shape[1]); ax.set_ylim(image.shape[0], 0); ax.axis("off"); fig.tight_layout(); fig.savefig(sample_out / "03_highcount_vs_labels.png", dpi=240, bbox_inches="tight"); plt.close(fig)
        summaries.append({"sample_id": sample, "image_shape_px": [int(image.shape[1]), int(image.shape[0])], "crop_shape_px": [int(image.shape[1]), int(image.shape[0])], "highcount_threshold": threshold, "highcount_bins": int(high.sum()), "cellpose_labels": int(labels.sum()), "highcount_label_barcode_overlap": int(overlap), "highcount_label_overlap_fraction": float(overlap / max(1, len(high_names))), "nearest_transform_scores": transform_scores, "figures": [str(p) for p in sorted(sample_out.glob("*.png"))]})
        adata.file.close()
    write_json(out / "alignment_clean_summary.json", summaries)
    rows = "".join(f"<tr><td>{html.escape(r['sample_id'])}</td><td>{r['highcount_threshold']:.1f}</td><td>{r['highcount_bins']}</td><td>{r['cellpose_labels']}</td><td>{r['highcount_label_barcode_overlap']}</td><td>{r['highcount_label_overlap_fraction']:.3f}</td></tr>" for r in summaries)
    sections = "".join(f"<h2>{html.escape(r['sample_id'])}</h2>" + "".join(f'<img src="{Path(p).relative_to(out / r["sample_id"])}" style="max-width:48%;margin:6px;border:1px solid #ddd">' for p in r["figures"]) for r in summaries)
    report = f"""<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><title>TD006859 clean alignment</title><style>body{{font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;margin:24px;color:#1f2933;background:#f7f9fc}}table{{border-collapse:collapse;background:#fff}}th,td{{border:1px solid #d9e2ec;padding:8px}}th{{background:#eef3f8}}</style></head><body><h1>TD006859 clean alignment</h1><p>使用 H&amp;E 原图、top 0.5% counts bins 和 Cellpose labels，避免全部 positive bins 遮盖背景。</p><table><tr><th>sample</th><th>high-count threshold</th><th>high-count bins</th><th>Cellpose labels</th><th>barcode overlap</th><th>overlap fraction</th></tr>{rows}</table>{sections}</body></html>"""
    (out / "alignment_clean.html").write_text(report)
    print(json.dumps({"status": "pass", "report": str(out / "alignment_clean.html"), "summary": str(out / "alignment_clean_summary.json")}, indent=2))


if __name__ == "__main__":
    main()
