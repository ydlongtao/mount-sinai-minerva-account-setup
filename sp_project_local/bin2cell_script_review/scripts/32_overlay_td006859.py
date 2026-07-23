#!/usr/bin/env python3
"""Create read-only H&E/Visium/Cellpose overlay diagnostics for TD006859."""
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

from batch_utils import BATCH_RESULTS, MANIFEST, read_manifest, write_json


TARGETS = {"TD006859-B408", "TD006859-B573"}
BARCODE = re.compile(r"^s_(\d+)um_(\d+)_(\d+)-\d+$")
CROP_X = (2000.0, 3000.0)
CROP_Y = (2000.0, 3000.0)


def decode(values) -> list[str]:
    return [x.decode() if isinstance(x, bytes) else str(x) for x in values]


def barcode_coordinates(path: Path) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5:
        matrix = h5["matrix"]
        barcodes = decode(matrix["barcodes"][:])
        indptr = np.asarray(matrix["indptr"][:], dtype=np.int64)
        data = matrix["data"]
        n_bins = len(barcodes)
        counts = np.zeros(n_bins, dtype=np.float64)
        for start in range(0, n_bins, 50000):
            stop = min(start + 50000, n_bins)
            data_start, data_stop = int(indptr[start]), int(indptr[stop])
            values = np.asarray(data[data_start:data_stop], dtype=np.float64)
            offsets = indptr[start:stop + 1] - data_start
            cumulative = np.concatenate(([0.0], np.cumsum(values)))
            counts[start:stop] = cumulative[offsets[1:]] - cumulative[offsets[:-1]]
    matches = [BARCODE.match(x) for x in barcodes]
    if any(x is None for x in matches):
        raise ValueError(f"Unparseable barcode in {path}")
    size = int(matches[0].group(1))
    x_um = np.asarray([int(m.group(3)) * size + size / 2 for m in matches], dtype=np.float64)
    y_um = np.asarray([int(m.group(2)) * size + size / 2 for m in matches], dtype=np.float64)
    return barcodes, x_um, y_um, counts


def crop_image(source: Path, x_um: np.ndarray, y_um: np.ndarray) -> tuple[np.ndarray, float, tuple[int, int, int, int]]:
    image = np.asarray(Image.open(source).convert("RGB"))
    pixel_per_um = min((image.shape[1] - 1) / max(x_um), (image.shape[0] - 1) / max(y_um))
    x0 = int(np.floor(CROP_X[0] * pixel_per_um))
    y0 = int(np.floor(CROP_Y[0] * pixel_per_um))
    x1 = min(image.shape[1], int(np.ceil(CROP_X[1] * pixel_per_um)))
    y1 = min(image.shape[0], int(np.ceil(CROP_Y[1] * pixel_per_um)))
    return image[y0:y1, x0:x1], float(pixel_per_um), (x0, y0, x1, y1)


def save_overlay(image: np.ndarray, x_px: np.ndarray, y_px: np.ndarray, values: np.ndarray | None, output: Path, title: str, cmap: str = "turbo", size: float = 1.2) -> None:
    fig, ax = plt.subplots(figsize=(9, 8))
    ax.imshow(image)
    if values is None:
        ax.scatter(x_px, y_px, s=size, c="#00d084", alpha=0.35, linewidths=0, rasterized=True)
    else:
        scatter = ax.scatter(x_px, y_px, s=size, c=values, cmap=cmap, alpha=0.55, linewidths=0, rasterized=True)
        fig.colorbar(scatter, ax=ax, fraction=0.035, pad=0.02)
    ax.set_title(title)
    ax.set_xlim(0, image.shape[1])
    ax.set_ylim(image.shape[0], 0)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def sample_overlay(sample: dict[str, str], out_dir: Path) -> dict[str, object]:
    sample_id = sample["sample_id"]
    sample_dir = out_dir / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)
    raw_h5 = Path(sample["matrix_002um"])
    source_image = Path(sample["outs_dir"]) / "spatial" / "tissue_hires_image.png"
    barcodes, x_um, y_um, counts = barcode_coordinates(raw_h5)
    image, pixel_per_um, crop_box = crop_image(source_image, x_um, y_um)
    in_crop = (x_um >= CROP_X[0]) & (x_um < CROP_X[1]) & (y_um >= CROP_Y[0]) & (y_um < CROP_Y[1])
    positive = in_crop & (counts > 0)
    crop_x_px = x_um[in_crop] * pixel_per_um - crop_box[0]
    crop_y_px = y_um[in_crop] * pixel_per_um - crop_box[1]
    positive_idx = np.flatnonzero(positive)
    if len(positive_idx) > 120000:
        rng = np.random.default_rng(0)
        positive_idx = np.sort(rng.choice(positive_idx, 120000, replace=False))
    positive_crop_x = x_um[positive_idx] * pixel_per_um - crop_box[0]
    positive_crop_y = y_um[positive_idx] * pixel_per_um - crop_box[1]
    save_overlay(image, positive_crop_x, positive_crop_y, None, sample_dir / "01_he_positive_2um_bins.png", f"{sample_id}: H&E + positive 2um bins", size=1.0)

    he_h5ad = BATCH_RESULTS / "samples" / sample_id / "cellpose" / f"{sample_id}_002um_cellpose_he.h5ad"
    labels = sc.read_h5ad(he_h5ad, backed="r")
    label_x = np.asarray(labels.obsm["spatial"][:, 0])
    label_y = np.asarray(labels.obsm["spatial"][:, 1])
    he_labels = np.asarray(labels.obs["labels_he"].to_numpy()) if "labels_he" in labels.obs else np.zeros(labels.n_obs, dtype=int)
    expanded = np.asarray(labels.obs["labels_he_expanded"].to_numpy()) if "labels_he_expanded" in labels.obs else np.zeros(labels.n_obs, dtype=int)
    label_positive = he_labels > 0
    expanded_positive = expanded > 0
    save_overlay(image, label_x[label_positive], label_y[label_positive], he_labels[label_positive], sample_dir / "02_he_cellpose_labels.png", f"{sample_id}: H&E Cellpose labels", size=1.4)
    save_overlay(image, label_x[expanded_positive], label_y[expanded_positive], expanded[expanded_positive], sample_dir / "03_he_expanded_labels.png", f"{sample_id}: expanded labels", size=1.2)

    cell_h5ad = BATCH_RESULTS / "samples" / sample_id / "cellpose" / f"{sample_id}_cell_level_cellpose_raw_counts.h5ad"
    cells = sc.read_h5ad(cell_h5ad, backed="r")
    cell_coords = np.asarray(cells.obsm["spatial"])
    save_overlay(image, cell_coords[:, 0], cell_coords[:, 1], None, sample_dir / "04_bin2cell_centroids.png", f"{sample_id}: bin2cell cell centroids", size=10.0)
    labels.file.close()
    cells.file.close()

    # Match by barcode string, never by row position: read_10x_h5 and the
    # Cellpose-derived H5AD are not required to preserve the same ordering.
    label_barcodes = np.asarray(labels.obs_names).astype(str)
    crop_positive_barcodes = set(np.asarray(barcodes)[positive])
    label_nonzero_barcodes = set(label_barcodes[label_positive])
    overlap = len(crop_positive_barcodes & label_nonzero_barcodes)
    unique_he = np.unique(he_labels[label_positive])
    unique_expanded = np.unique(expanded[expanded_positive])
    area_values = np.bincount(he_labels[label_positive].astype(int)) if label_positive.any() else np.array([], dtype=int)
    area_values = area_values[area_values > 0]
    return {
        "sample_id": sample_id,
        "image_shape_px": [int(image.shape[1]), int(image.shape[0])],
        "pixel_per_um": pixel_per_um,
        "crop_box_px": list(crop_box),
        "crop_bins": int(in_crop.sum()),
        "positive_crop_bins": int(positive.sum()),
        "he_h5ad": str(he_h5ad),
        "cell_h5ad": str(cell_h5ad),
        "labels_he_nonzero_bins": int(label_positive.sum()),
        "labels_he_expanded_nonzero_bins": int(expanded_positive.sum()),
        "labels_he_unique_objects": int(len(unique_he)),
        "labels_he_expanded_unique_objects": int(len(unique_expanded)),
        "label_positive_fraction": float(label_positive.mean()),
        "expanded_positive_fraction": float(expanded_positive.mean()),
        "positive_label_overlap_bins": int(overlap),
        "positive_label_overlap_fraction": float(overlap / max(1, positive.sum())),
        "cell_centroids": int(cells.n_obs),
        "label_area_quantiles_bins": {k: float(v) for k, v in zip(["q50", "q95", "q995"], np.percentile(area_values, [50, 95, 99.5]))} if area_values.size else {},
        "figures": [str(p) for p in sorted(sample_dir.glob("*.png"))],
    }


def main() -> None:
    out_dir = BATCH_RESULTS / "debug_td006859" / "overlay"
    out_dir.mkdir(parents=True, exist_ok=True)
    results = [sample_overlay(row, out_dir) for row in read_manifest(MANIFEST) if row["sample_id"] in TARGETS]
    write_json(out_dir / "td006859_overlay_summary.json", results)
    rows = "".join(f"<tr><td>{html.escape(r['sample_id'])}</td><td>{r['positive_crop_bins']}</td><td>{r['labels_he_nonzero_bins']}</td><td>{r['labels_he_expanded_nonzero_bins']}</td><td>{r['labels_he_unique_objects']}</td><td>{r['cell_centroids']}</td><td>{r['positive_label_overlap_fraction']:.3f}</td></tr>" for r in results)
    sections = []
    for r in results:
        imgs = "".join(f'<figure><img src="{Path(p).name}" style="max-width:95%"><figcaption>{html.escape(Path(p).stem)}</figcaption></figure>' for p in r["figures"])
        sections.append(f"<h2>{html.escape(r['sample_id'])}</h2><p><code>{html.escape(json.dumps(r, indent=2, ensure_ascii=False))}</code></p><div class=\"grid\">{imgs}</div>")
    report = f"""<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>TD006859 overlay</title><style>body{{font:14px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;color:#1f2933;background:#f7f9fc;margin:0}}main{{max-width:1400px;margin:auto;padding:28px 18px 60px}}table{{width:100%;border-collapse:collapse;background:#fff}}th,td{{border:1px solid #d9e2ec;padding:8px;text-align:left}}th{{background:#eef3f8}}.grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}}figure{{background:#fff;margin:0;padding:8px;border:1px solid #d9e2ec}}img{{display:block;width:100%}}figcaption{{color:#52606d;padding:5px}}code{{white-space:pre-wrap}}</style></head><body><main><h1>TD006859-B408 / B573 H&amp;E overlay 诊断</h1><p>本报告使用现有 Cellpose 结果和原始 2um positive bins，只读生成，不重跑分割。</p><table><tr><th>sample</th><th>positive crop bins</th><th>HE label bins</th><th>expanded bins</th><th>HE objects</th><th>cell centroids</th><th>positive-label overlap</th></tr>{rows}</table>{''.join(sections)}</main></body></html>"""
    (out_dir / "td006859_overlay.html").write_text(report)
    print(json.dumps({"status": "pass", "report": str(out_dir / "td006859_overlay.html"), "summary": str(out_dir / "td006859_overlay_summary.json")}, indent=2))


if __name__ == "__main__":
    main()
