#!/usr/bin/env python3
"""Validate batch v1 result files and write a corrected summary report.

This script is read-only with respect to H5AD inputs. It records structural
issues and regenerates completion counts from the actual result files.
"""
from __future__ import annotations

import csv
import html
import json
import os
from pathlib import Path

import numpy as np
import scanpy as sc
from scipy import sparse

from batch_utils import BATCH_RESULTS, read_manifest, sample_root, write_json


STAGES = {
    "016um": ("016um", "{sample}_016um_processed.h5ad", True),
    "008um": ("008um", "{sample}_008um_mapped.h5ad", False),
    "cellpose_he": ("cellpose", "{sample}_002um_cellpose_he.h5ad", True),
    "cell_raw_counts": ("cellpose", "{sample}_cell_level_cellpose_raw_counts.h5ad", True),
    "cell_qc": ("cell_qc_marker", "{sample}_cell_level_qc_marker.h5ad", True),
    "spatial_domains": ("spatial_domains", "{sample}_spatial_domains.h5ad", True),
}


def check_h5ad(path: Path) -> dict[str, object]:
    result: dict[str, object] = {"path": str(path), "exists": path.exists()}
    if not path.exists():
        result.update(status="missing", errors=["file does not exist"])
        return result
    try:
        adata = sc.read_h5ad(path, backed="r")
        errors: list[str] = []
        warnings: list[str] = []
        result["shape"] = [int(adata.n_obs), int(adata.n_vars)]
        result["obs_names_unique"] = bool(adata.obs_names.is_unique)
        result["var_names_unique"] = bool(adata.var_names.is_unique)
        if not adata.obs_names.is_unique:
            errors.append("obs_names are not unique")
        if not adata.var_names.is_unique:
            warnings.append("var_names are not unique; make unique before downstream modeling")

        spatial = adata.obsm.get("spatial")
        if spatial is None:
            errors.append("obsm['spatial'] is missing")
        else:
            coords = np.asarray(spatial)
            result["spatial_shape"] = list(coords.shape)
            result["spatial_finite"] = bool(np.isfinite(coords).all())
            if coords.ndim != 2 or coords.shape != (adata.n_obs, 2):
                errors.append(f"spatial shape is {list(coords.shape)}, expected [{adata.n_obs}, 2]")
            elif not np.isfinite(coords).all():
                errors.append("spatial contains non-finite values")

        layers = list(adata.layers.keys())
        result["layers"] = layers
        result["counts_layer_present"] = "counts" in layers
        if "counts" not in layers:
            warnings.append("counts layer is absent; X must be treated as the raw-count source only if confirmed")
        else:
            counts = adata.layers["counts"]
            result["counts_shape"] = list(counts.shape)
            if tuple(counts.shape) != (adata.n_obs, adata.n_vars):
                errors.append("counts layer shape does not match X")

        # Backed sparse H5AD matrices expose a CSRDataset object. Inspect a
        # bounded row sample instead of applying NumPy reductions to the
        # dataset proxy itself.
        x_sample = adata.X[: min(256, adata.n_obs)]
        if sparse.issparse(x_sample):
            x_values = np.asarray(x_sample.data)
        else:
            x_values = np.asarray(x_sample)
        result["x_dtype"] = str(getattr(adata.X, "dtype", x_values.dtype))
        result["x_nonnegative"] = bool(x_values.size == 0 or np.nanmin(x_values) >= 0)
        if not result["x_nonnegative"]:
            errors.append("X contains negative values")
        result["status"] = "fail" if errors else ("warn" if warnings else "pass")
        result["errors"] = errors
        result["warnings"] = warnings
        adata.file.close()
        return result
    except Exception as exc:  # noqa: BLE001 - report all unreadable files
        result.update(status="fail", errors=[f"read failed: {type(exc).__name__}: {exc}"], warnings=[])
        return result


def percent(value: int, total: int) -> str:
    return f"{100.0 * value / total:.1f}%" if total else "0.0%"


def main() -> None:
    out_dir = BATCH_RESULTS / "validation"
    report_dir = BATCH_RESULTS / "reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    sample_rows = []
    all_checks = []
    manifest = read_manifest()
    for row in manifest:
        sample_id = row["sample_id"]
        stage_status: dict[str, str] = {}
        stage_checks: dict[str, dict[str, object]] = {}
        for stage, (folder, pattern, required) in STAGES.items():
            path = sample_root(sample_id) / folder / pattern.format(sample=sample_id)
            check = check_h5ad(path)
            check["required"] = required
            stage_checks[stage] = check
            all_checks.append({"sample_id": sample_id, "stage": stage, **check})
            stage_status[stage] = str(check["status"])
        sample_rows.append({"sample_id": sample_id, "group": row.get("group", ""), "stages": stage_status, "checks": stage_checks})

    summary = {}
    for stage in STAGES:
        checks = [s["checks"][stage] for s in sample_rows]
        present = sum(bool(c.get("exists")) for c in checks)
        readable = sum(c.get("status") in {"pass", "warn"} for c in checks)
        summary[stage] = {"present": present, "readable": readable, "total": len(checks), "fraction": percent(present, len(checks))}

    integration_paths = {
        "integrated_016um": BATCH_RESULTS / "integration/016um/integrated_016um.h5ad",
        "integrated_cell_level": BATCH_RESULTS / "integration/cell_level/integrated_cell_level.h5ad",
    }
    integration_checks = {name: check_h5ad(path) for name, path in integration_paths.items()}
    validation = {"status": "pass" if not any(c["status"] == "fail" for c in all_checks) else "fail", "samples": sample_rows, "stage_summary": summary, "integration": integration_checks}
    write_json(out_dir / "batch_validation.json", validation)

    csv_path = out_dir / "batch_validation_table.csv"
    fields = ["sample_id", "stage", "status", "exists", "shape", "spatial_shape", "spatial_finite", "var_names_unique", "counts_layer_present", "layers", "errors", "warnings"]
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for item in all_checks:
            writer.writerow({key: json.dumps(item.get(key), ensure_ascii=False) if isinstance(item.get(key), (list, dict)) else item.get(key, "") for key in fields})

    rows_html = []
    for sample in sample_rows:
        stage_bits = "; ".join(f"{stage}: {sample['stages'][stage]}" for stage in STAGES)
        rows_html.append(f"<tr><td>{html.escape(sample['sample_id'])}</td><td>{html.escape(sample['group'])}</td><td>{html.escape(stage_bits)}</td></tr>")
    stage_html = "".join(f"<tr><td>{stage}</td><td>{v['present']}/{v['total']}</td><td>{v['readable']}/{v['total']}</td><td>{v['fraction']}</td></tr>" for stage, v in summary.items())
    integration_html = "".join(f"<li><code>{html.escape(name)}</code>: {html.escape(str(check['status']))}; shape={html.escape(str(check.get('shape', 'unknown')))}</li>" for name, check in integration_checks.items())
    report = f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Visium HD Batch v1 Step 1 Validation</title>
<style>body{{font:15px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;color:#1f2933;background:#f7f9fc;margin:0}}main{{max-width:1100px;margin:auto;padding:30px 20px 60px}}h1{{font-size:28px}}table{{width:100%;border-collapse:collapse;background:#fff;margin:12px 0 22px}}th,td{{border:1px solid #d9e2ec;padding:8px;text-align:left;vertical-align:top}}th{{background:#eef3f8}}.box{{background:#fff;border-left:5px solid #0b5cad;padding:12px 16px;margin:14px 0}}code{{background:#eef1f4;padding:1px 4px}}</style></head>
<body><main><h1>Visium HD 批量 v1：Step 1 结果验证报告</h1>
<p>验证时间：{html.escape(str(os.environ.get('LSB_JOBID', 'manual')))}；状态：<strong>{validation['status']}</strong></p>
<div class="box"><strong>统计修正：</strong>完成数由实际 H5AD 文件重新计算。16um 单样本结果为 <strong>{summary['016um']['present']}/{summary['016um']['total']}</strong>，不是此前本地报告中的 0/10。</div>
<h2>阶段文件统计</h2><table><tr><th>阶段</th><th>文件存在</th><th>可读取</th><th>存在率</th></tr>{stage_html}</table>
<h2>整合文件</h2><ul>{integration_html}</ul>
<h2>样本级状态</h2><table><tr><th>样本</th><th>组别</th><th>阶段检查</th></tr>{''.join(rows_html)}</table>
<h2>检查定义</h2><ul><li>H5AD：使用 backed 模式读取，避免把大型矩阵全部载入内存。</li><li>空间坐标：检查 <code>obsm['spatial']</code> 是否存在、形状是否为 n_obs × 2、是否全为有限值。</li><li>gene names：检查 <code>var_names</code> 是否唯一；重复名称记录为 warning，不静默修改原始 H5AD。</li><li>counts：记录 <code>layers['counts']</code> 是否存在并检查形状；没有 counts layer 的文件要求后续确认 X 是否为 raw counts。</li><li>表达矩阵：检查 X 是否包含负值。</li></ul>
<p>详细机器可读结果：<code>validation/batch_validation.json</code>；逐文件表：<code>validation/batch_validation_table.csv</code>。</p>
</main></body></html>"""
    (report_dir / "spatial_batch_summary_corrected.html").write_text(report)
    print(json.dumps({"status": validation["status"], "stage_summary": summary, "report": str(report_dir / "spatial_batch_summary_corrected.html")}, indent=2))


if __name__ == "__main__":
    main()
