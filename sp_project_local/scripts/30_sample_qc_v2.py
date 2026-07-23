#!/usr/bin/env python3
"""Compute sample-level QC distributions for the formal Visium HD workflow."""
from __future__ import annotations

import csv
import html
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from batch_utils import BATCH_RESULTS, read_manifest, sample_root, write_json


STAGES = {
    "016um": ("016um", "{sample}_016um_processed.h5ad"),
    "008um": ("008um", "{sample}_008um_mapped.h5ad"),
    "cell_level": ("cell_qc_marker", "{sample}_cell_level_qc_marker.h5ad"),
}


def row_qc(adata, chunk_size: int = 10000) -> tuple[np.ndarray, np.ndarray]:
    counts = np.empty(adata.n_obs, dtype=np.float64)
    genes = np.empty(adata.n_obs, dtype=np.int32)
    for start in range(0, adata.n_obs, chunk_size):
        stop = min(start + chunk_size, adata.n_obs)
        block = adata.X[start:stop]
        if sparse.issparse(block):
            counts[start:stop] = np.asarray(block.sum(axis=1)).ravel()
            genes[start:stop] = np.diff(block.indptr)
        else:
            dense = np.asarray(block)
            counts[start:stop] = dense.sum(axis=1)
            genes[start:stop] = (dense > 0).sum(axis=1)
    return counts, genes


def quantiles(values: np.ndarray) -> dict[str, float]:
    positive = values[np.isfinite(values) & (values > 0)]
    if positive.size == 0:
        return {"q01": 0.0, "q05": 0.0, "median": 0.0, "q95": 0.0, "q995": 0.0}
    q = np.percentile(positive, [1, 5, 50, 95, 99.5])
    return {key: float(value) for key, value in zip(["q01", "q05", "median", "q95", "q995"], q)}


def stage_qc(sample_id: str, stage: str) -> dict[str, object]:
    folder, pattern = STAGES[stage]
    path = sample_root(sample_id) / folder / pattern.format(sample=sample_id)
    result: dict[str, object] = {"sample_id": sample_id, "stage": stage, "path": str(path), "exists": path.exists()}
    if not path.exists():
        result.update(status="missing", review_flags=["missing_result"])
        return result
    try:
        adata = sc.read_h5ad(path, backed="r")
        counts, genes = row_qc(adata)
        positive = counts > 0
        c_q = quantiles(counts)
        g_q = quantiles(genes.astype(float))
        result.update({
            "status": "pass",
            "shape": [int(adata.n_obs), int(adata.n_vars)],
            "positive_obs": int(positive.sum()),
            "positive_fraction": float(positive.mean()),
            "zero_count_obs": int((~positive).sum()),
            "counts_quantiles": c_q,
            "genes_quantiles": g_q,
            "candidate_min_counts": float(max(1.0, c_q["q01"])),
            "candidate_max_counts": float(c_q["q995"]),
            "candidate_min_genes": float(max(1.0, g_q["q01"])),
            "candidate_max_genes": float(g_q["q995"]),
            "counts_layer_present": "counts" in adata.layers,
            "var_names_unique": bool(adata.var_names.is_unique),
            "spatial_complete": bool("spatial" in adata.obsm and np.isfinite(np.asarray(adata.obsm["spatial"])).all()),
            "review_flags": [],
        })
        if not result["counts_layer_present"]:
            result["review_flags"].append("counts_layer_absent_X_is_source_check_needed")
        if not result["var_names_unique"]:
            result["review_flags"].append("duplicate_gene_names")
        if not result["spatial_complete"]:
            result["review_flags"].append("spatial_coordinates_incomplete")
        if stage == "cell_level":
            if adata.n_obs < 100:
                result["review_flags"].append("too_few_cells_for_cell_level_statistics")
            if g_q["median"] < 50:
                result["review_flags"].append("low_median_detected_genes")
        if stage == "008um":
            mapping_summary = path.parent / "008um_mapping_summary.json"
            if mapping_summary.exists():
                mapping = json.loads(mapping_summary.read_text())
                result["mapping_fraction"] = mapping.get("mapping_fraction")
                result["mapping_status"] = mapping.get("status")
                if float(mapping.get("mapping_fraction", 0.0)) < 0.99:
                    result["review_flags"].append("mapping_fraction_below_99_percent")
        adata.file.close()
        return result
    except Exception as exc:  # noqa: BLE001
        result.update(status="fail", error=f"{type(exc).__name__}: {exc}", review_flags=["read_or_qc_failed"])
        return result


def html_report(rows: list[dict[str, object]], out_path: Path) -> None:
    headers = ["sample_id", "stage", "status", "shape", "positive_fraction", "counts_median", "genes_median", "counts_layer_present", "mapping_fraction", "review_flags"]
    body = []
    for row in rows:
        body.append("<tr>" + "".join(f"<td>{html.escape(str(row.get(key, '')))}</td>" for key in headers) + "</tr>")
    flag_counts: dict[str, int] = {}
    for row in rows:
        for flag in row.get("review_flags", []):
            flag_counts[flag] = flag_counts.get(flag, 0) + 1
    flags_html = "".join(f"<li><code>{html.escape(k)}</code>: {v}</li>" for k, v in sorted(flag_counts.items())) or "<li>无</li>"
    table_head = "".join(f"<th>{html.escape(x)}</th>" for x in headers)
    report = f"""<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Visium HD Sample QC v2</title>
<style>body{{font:14px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;color:#1f2933;background:#f7f9fc;margin:0}}main{{max-width:1500px;margin:auto;padding:28px 18px 60px}}table{{border-collapse:collapse;width:100%;background:#fff;font-size:12px}}th,td{{border:1px solid #d9e2ec;padding:7px;text-align:left;vertical-align:top}}th{{background:#eef3f8;position:sticky;top:0}}.box{{background:#fff;border-left:5px solid #0b5cad;padding:12px 16px;margin:12px 0}}code{{background:#eef1f4;padding:1px 4px}}</style></head><body><main>
<h1>Visium HD 样本级 QC v2</h1><div class="box"><strong>用途：</strong>候选 QC 阈值和异常样本复核，不自动删除任何 bin 或 cell。阈值按每个样本的阳性观测分布计算。</div>
<h2>复核标记</h2><ul>{flags_html}</ul><h2>逐样本逐分辨率结果</h2><table><tr>{table_head}</tr>{''.join(body)}</table>
<h2>解释</h2><ul><li><code>positive_fraction</code>：X 中 total counts 大于 0 的观测比例。</li><li><code>counts_median</code> 和 <code>genes_median</code>：阳性观测的中位数。</li><li><code>candidate_min/max_*</code>：1% 和 99.5% 分位数形成的候选范围，需结合 H&amp;E 和空间覆盖确认。</li><li>TD006859-B408/B573 的细胞数不足时，只做描述性展示，不进入正式组间统计。</li></ul>
</main></body></html>"""
    out_path.write_text(report)


def main() -> None:
    out_dir = BATCH_RESULTS / "qc_v2"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = [stage_qc(row["sample_id"], stage) for row in read_manifest() for stage in STAGES]
    write_json(out_dir / "sample_qc_v2.json", rows)
    fields = ["sample_id", "stage", "status", "shape", "positive_fraction", "zero_count_obs", "counts_quantiles", "genes_quantiles", "candidate_min_counts", "candidate_max_counts", "candidate_min_genes", "candidate_max_genes", "counts_layer_present", "var_names_unique", "spatial_complete", "mapping_fraction", "mapping_status", "review_flags", "error"]
    with (out_dir / "sample_qc_v2.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(row.get(key), ensure_ascii=False) if isinstance(row.get(key), (list, dict)) else row.get(key, "") for key in fields})
    report_path = BATCH_RESULTS / "reports" / "sample_qc_v2.html"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    html_report(rows, report_path)
    failed = sum(row.get("status") == "fail" for row in rows)
    print(json.dumps({"status": "fail" if failed else "pass", "rows": len(rows), "failed": failed, "report": str(report_path), "output_dir": str(out_dir)}, indent=2))


if __name__ == "__main__":
    main()
