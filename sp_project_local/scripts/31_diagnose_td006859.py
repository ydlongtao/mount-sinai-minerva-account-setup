#!/usr/bin/env python3
"""Read-only diagnosis for TD006859-B408/B573 Visium HD inputs."""
from __future__ import annotations

import html
import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from PIL import Image

from batch_utils import BATCH_RESULTS, MANIFEST, read_manifest, write_json


TARGETS = {"TD006859-B408", "TD006859-B573"}
BARCODE = re.compile(r"^s_(\d+)um_(\d+)_(\d+)-\d+$")
CROP_X = (2000.0, 3000.0)
CROP_Y = (2000.0, 3000.0)


def decode(values) -> list[str]:
    return [x.decode() if isinstance(x, bytes) else str(x) for x in values]


def summarize_h5(path: Path) -> dict[str, object]:
    result: dict[str, object] = {"path": str(path), "exists": path.exists()}
    if not path.exists():
        result["status"] = "missing"
        return result
    with h5py.File(path, "r") as h5:
        matrix = h5["matrix"]
        shape = [int(x) for x in matrix["shape"][:]]
        barcodes = decode(matrix["barcodes"][:])
        n_genes, n_bins = shape
        matches = [BARCODE.match(x) for x in barcodes]
        if any(x is None for x in matches):
            raise ValueError(f"unparseable barcode in {path}")
        size = int(matches[0].group(1))
        x = np.asarray([int(m.group(3)) * size + size / 2 for m in matches], dtype=np.float64)
        y = np.asarray([int(m.group(2)) * size + size / 2 for m in matches], dtype=np.float64)
        indptr = np.asarray(matrix["indptr"][:], dtype=np.int64)
        data = matrix["data"]
        total_counts = np.zeros(n_bins, dtype=np.float64)
        for start in range(0, n_bins, 50000):
            stop = min(start + 50000, n_bins)
            data_start, data_stop = int(indptr[start]), int(indptr[stop])
            values = np.asarray(data[data_start:data_stop], dtype=np.float64)
            offsets = indptr[start:stop + 1] - data_start
            cumulative = np.concatenate(([0.0], np.cumsum(values)))
            total_counts[start:stop] = cumulative[offsets[1:]] - cumulative[offsets[:-1]]
        positive = total_counts > 0
        result.update({
            "status": "pass",
            "shape": [n_bins, n_genes],
            "bin_size_um": size,
            "barcode_count": len(barcodes),
            "coordinate_range_um": {"x": [float(x.min()), float(x.max())], "y": [float(y.min()), float(y.max())]},
            "positive_bins": int(positive.sum()),
            "positive_fraction": float(positive.mean()),
            "positive_coordinate_range_um": {
                "x": [float(x[positive].min()), float(x[positive].max())] if positive.any() else None,
                "y": [float(y[positive].min()), float(y[positive].max())] if positive.any() else None,
            },
            "positive_median_counts": float(np.median(total_counts[positive])) if positive.any() else 0.0,
            "current_crop_positive_bins": int((positive & (x >= CROP_X[0]) & (x < CROP_X[1]) & (y >= CROP_Y[0]) & (y < CROP_Y[1])).sum()),
            "current_crop_all_bins": int(((x >= CROP_X[0]) & (x < CROP_X[1]) & (y >= CROP_Y[0]) & (y < CROP_Y[1])).sum()),
            "barcodes": barcodes,
        })
    return result


def image_summary(path: Path) -> dict[str, object]:
    result = {"path": str(path), "exists": path.exists()}
    if path.exists():
        with Image.open(path) as image:
            result.update({"size_px": [int(image.width), int(image.height)], "mode": image.mode})
    return result


def mapping_summary(path: Path, barcodes8: list[str], barcodes16: list[str]) -> dict[str, object]:
    result: dict[str, object] = {"path": str(path), "exists": path.exists()}
    if not path.exists():
        result.update({"status": "missing", "review_flags": ["barcode_mappings_missing"]})
        return result
    df = pd.read_parquet(path)
    result["status"] = "pass"
    result["shape"] = [int(df.shape[0]), int(df.shape[1])]
    result["columns"] = list(df.columns)
    required = {"square_008um", "square_016um"}
    result["required_columns_present"] = required.issubset(df.columns)
    flags = []
    if not required.issubset(df.columns):
        flags.append("mapping_columns_missing")
    else:
        set8, set16 = set(barcodes8), set(barcodes16)
        result["square_008um_unique"] = int(df["square_008um"].nunique())
        result["square_016um_unique"] = int(df["square_016um"].nunique())
        result["8um_intersection"] = int(df["square_008um"].isin(set8).sum())
        result["16um_intersection"] = int(df["square_016um"].isin(set16).sum())
        if result["8um_intersection"] == 0 or result["16um_intersection"] == 0:
            flags.append("barcode_format_or_sample_mismatch")
    result["review_flags"] = flags
    return result


def main() -> None:
    out_dir = BATCH_RESULTS / "debug_td006859" / "diagnostic"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for row in read_manifest(MANIFEST):
        sample = row["sample_id"]
        if sample not in TARGETS:
            continue
        h5 = {}
        for stage in ["002um", "008um", "016um"]:
            h5[stage] = summarize_h5(Path(row[f"matrix_{stage}"]))
        image = image_summary(Path(row["outs_dir"]) / "spatial" / "tissue_hires_image.png")
        mapping = mapping_summary(Path(row["outs_dir"]) / "barcode_mappings.parquet", h5["008um"]["barcodes"], h5["016um"]["barcodes"])
        rows.append({"sample_id": sample, "group": row.get("group", ""), "h5": h5, "image": image, "mapping": mapping})
        for stage in h5:
            h5[stage].pop("barcodes", None)
    # A compact decision-oriented table for the HTML report.
    table_rows = []
    for row in rows:
        for stage, info in row["h5"].items():
            pos = info.get("positive_coordinate_range_um") or {}
            table_rows.append({
                "sample_id": row["sample_id"], "stage": stage, "shape": info.get("shape"),
                "positive_bins": info.get("positive_bins"), "positive_fraction": info.get("positive_fraction"),
                "positive_x": pos.get("x"), "positive_y": pos.get("y"),
                "crop_positive_bins": info.get("current_crop_positive_bins"), "crop_all_bins": info.get("current_crop_all_bins"),
            })
    payload = {"crop_um": {"x": CROP_X, "y": CROP_Y}, "samples": rows, "table": table_rows}
    write_json(out_dir / "td006859_diagnostic.json", payload)
    columns = ["sample_id", "stage", "shape", "positive_bins", "positive_fraction", "positive_x", "positive_y", "crop_positive_bins", "crop_all_bins"]
    pd.DataFrame(table_rows, columns=columns).to_csv(out_dir / "td006859_coordinate_summary.csv", index=False)

    html_rows = "".join("<tr>" + "".join(f"<td>{html.escape(str(item.get(c, '')))}</td>" for c in columns) + "</tr>" for item in table_rows)
    map_items = []
    for row in rows:
        m = row["mapping"]
        map_items.append(f"<h3>{html.escape(row['sample_id'])}</h3><pre>{html.escape(json.dumps(m, indent=2, ensure_ascii=False))}</pre>")
    image_items = "".join(f"<li>{html.escape(row['sample_id'])}: {html.escape(json.dumps(row['image'], ensure_ascii=False))}</li>" for row in rows)
    head = "".join(f"<th>{c}</th>" for c in columns)
    report = f"""<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>TD006859 diagnostic</title>
<style>body{{font:14px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;color:#1f2933;background:#f7f9fc;margin:0}}main{{max-width:1400px;margin:auto;padding:28px 18px 60px}}table{{width:100%;border-collapse:collapse;background:#fff;font-size:12px}}th,td{{border:1px solid #d9e2ec;padding:7px;text-align:left;vertical-align:top}}th{{background:#eef3f8}}.box{{background:#fff;border-left:5px solid #0b5cad;padding:12px 16px;margin:12px 0}}pre{{background:#fff;border:1px solid #d9e2ec;padding:12px;overflow:auto}}code{{background:#eef1f4;padding:1px 4px}}</style></head><body><main>
<h1>TD006859-B408 / B573 原始数据只读诊断</h1><div class="box"><strong>当前 crop：</strong>x={CROP_X[0]}–{CROP_X[1]} um；y={CROP_Y[0]}–{CROP_Y[1]} um。若 crop 内 positive bins 很少而全图 positive bins 较多，优先判定为裁剪位置问题。</div>
<h2>坐标与 crop 覆盖</h2><table><tr>{head}</tr>{html_rows}</table><h2>H&amp;E 图像</h2><ul>{image_items}</ul><h2>Barcode mapping</h2>{''.join(map_items)}<h2>下一步判定</h2><ol><li>先检查 positive coordinate range 与当前 crop 的关系。</li><li>若组织在 crop 外，进入重新裁剪步骤。</li><li>若组织在 crop 内，再检查 H&amp;E 与空间坐标 overlay。</li><li>只有 mapping 文件存在且 barcode 可匹配时，才重建 8um mapped 结果。</li></ol>
</main></body></html>"""
    (out_dir / "td006859_diagnostic.html").write_text(report)
    print(json.dumps({"status": "pass", "samples": [x["sample_id"] for x in rows], "report": str(out_dir / "td006859_diagnostic.html"), "output_dir": str(out_dir)}, indent=2))


if __name__ == "__main__":
    main()
