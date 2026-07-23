#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from batch_utils import BATCH_RESULTS, load_batch_config, read_json, read_manifest, sample_root


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(BATCH_RESULTS))
    except ValueError:
        return str(path)


def status_cell(path: Path) -> str:
    if path.exists():
        return f'<span class="ok">done</span><br><code>{rel(path)}</code>'
    return '<span class="missing">missing</span>'


def main() -> None:
    config = load_batch_config()
    sample_ids = config.get("sample_ids") or [row["sample_id"] for row in read_manifest()]
    out_dir = BATCH_RESULTS / "reports"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for sample_id in sample_ids:
        root = sample_root(sample_id)
        rows.append({
            "sample_id": sample_id,
            "016um": root / "016um" / f"{sample_id}_016um_processed.h5ad",
            "008um": root / "008um" / f"{sample_id}_008um_mapped.h5ad",
            "cellpose": root / "cellpose" / f"{sample_id}_cell_level_cellpose_raw_counts.h5ad",
            "cell_qc": root / "cell_qc_marker" / "cell_qc_marker_summary.json",
            "domains": root / "spatial_domains" / "spatial_domain_summary.json",
        })

    qc_path = BATCH_RESULTS / "qc" / "sample_qc_summary.csv"
    qc_html = ""
    if qc_path.exists():
        qc_df = pd.read_csv(qc_path)
        qc_html = qc_df.to_html(index=False, classes="dataframe")

    integration_016 = read_json(BATCH_RESULTS / "integration" / "016um" / "integrated_016um_summary.json", default={})
    integration_cell = read_json(BATCH_RESULTS / "integration" / "cell_level" / "integrated_cell_level_summary.json", default={})

    table_rows = "\n".join(
        "<tr>"
        f"<td>{row['sample_id']}</td>"
        f"<td>{status_cell(row['016um'])}</td>"
        f"<td>{status_cell(row['008um'])}</td>"
        f"<td>{status_cell(row['cellpose'])}</td>"
        f"<td>{status_cell(row['cell_qc'])}</td>"
        f"<td>{status_cell(row['domains'])}</td>"
        "</tr>"
        for row in rows
    )
    html = f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><title>Visium HD Batch v1 Report</title>
<style>
body{{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;margin:30px;line-height:1.55;color:#1f2933}}
table{{border-collapse:collapse;margin:16px 0;width:100%}}td,th{{border:1px solid #d7dde5;padding:7px 9px;vertical-align:top}}th{{background:#eef3f8}}
code{{background:#f3f4f6;padding:2px 4px}}.ok{{color:#0f766e;font-weight:700}}.missing{{color:#b91c1c;font-weight:700}}
img{{max-width:48%;border:1px solid #ddd;margin:6px}}
</style></head><body>
<h1>Visium HD Prostate Cancer Batch v1 Report</h1>
<p>本报告汇总正式批量流程的文件状态。v1 重点包括样本 QC、16um 全局建模、8um 映射复核、H&E Cellpose 细胞分割、raw-count cell-level 矩阵、marker/QC 和空间域初步识别。空间通讯和 transition tensor 暂不纳入 v1。</p>
<h2>Run Configuration</h2>
<pre>{json.dumps(config, indent=2, sort_keys=True)}</pre>
<h2>Sample Status</h2>
<table><tr><th>sample</th><th>16um model</th><th>8um mapped</th><th>cellpose raw cells</th><th>cell QC</th><th>spatial domains</th></tr>{table_rows}</table>
<h2>Sample QC Summary</h2>
{qc_html or '<p class="missing">QC summary not available yet.</p>'}
<h2>Integration</h2>
<h3>16um</h3><pre>{json.dumps(integration_016, indent=2, sort_keys=True)}</pre>
<h3>Cell Level</h3><pre>{json.dumps(integration_cell, indent=2, sort_keys=True)}</pre>
<h2>GPU/CPU Policy</h2>
<p>GPU 节点仅用于 H&E Cellpose 分割。16um/8um 矩阵处理、raw-count 聚合、marker/QC、整合和空间域模型默认使用 CPU 队列，减少 GPU 队列等待并避免非 GPU 步骤占用 GPU。</p>
</body></html>"""
    report_html = out_dir / "spatial_batch_summary.html"
    report_html.write_text(html)

    md = ["# Visium HD Batch v1 Report", "", "## Sample Status", ""]
    for row in rows:
        done = [name for name in ["016um", "008um", "cellpose", "cell_qc", "domains"] if row[name].exists()]
        md.append(f"- {row['sample_id']}: {', '.join(done) if done else 'no outputs yet'}")
    md.append("")
    md.append(f"HTML report: `{report_html}`")
    (out_dir / "spatial_batch_summary.md").write_text("\n".join(md) + "\n")
    print(json.dumps({"status": "pass", "html": str(report_html), "samples": sample_ids}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
