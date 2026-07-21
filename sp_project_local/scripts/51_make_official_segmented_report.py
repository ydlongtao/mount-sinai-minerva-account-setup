#!/usr/bin/env python3
from __future__ import annotations

import csv
import html
import json
import os
from pathlib import Path


PROJECT_HOME = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
MANIFEST = Path(os.environ.get("SP_PROJECT_MANIFEST", PROJECT_HOME / "config" / "sample_manifest_v2_corrected.csv"))
RESULTS = Path(os.environ.get("SP_OFFICIAL_RESULTS", PROJECT_HOME / "results" / "segmented_official_v1"))


def main() -> None:
    samples = [row["sample_id"] for row in csv.DictReader(MANIFEST.open(newline=""))]
    summaries = []
    for sample in samples:
        path = RESULTS / "samples" / sample / "official_segmented_summary.json"
        if not path.is_file():
            raise FileNotFoundError(path)
        summaries.append(json.loads(path.read_text()))

    reports = RESULTS / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    (reports / "official_segmented_qc_summary.json").write_text(
        json.dumps(summaries, indent=2, sort_keys=True) + "\n"
    )
    with (reports / "official_segmented_qc_summary.csv").open("w", newline="") as handle:
        fields = [
            "sample_id", "cells", "genes", "polygon_features", "matched_polygon_percent",
            "median_total_counts", "median_genes", "median_pct_mt", "median_cell_area_um2",
            "passes_qc_candidate", "passes_qc_candidate_percent",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for item in summaries:
            writer.writerow({key: item[key] for key in fields})

    rows = []
    panels = []
    for item in summaries:
        sample = item["sample_id"]
        rel = f"../samples/{sample}/figures/official_cell_boundaries_he_overview.png"
        rows.append(
            "<tr>"
            f"<td>{html.escape(sample)}</td><td>{item['cells']:,}</td><td>{item['median_total_counts']:.1f}</td>"
            f"<td>{item['median_genes']:.1f}</td><td>{item['median_pct_mt']:.2f}</td>"
            f"<td>{item['median_cell_area_um2']:.1f}</td><td>{item['passes_qc_candidate_percent']:.1f}%</td>"
            "</tr>"
        )
        panels.append(
            f'<section><h2>{html.escape(sample)}</h2><img src="{rel}" alt="{html.escape(sample)} official cell boundaries"></section>'
        )

    page = f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>8 样本 Space Ranger 官方细胞分割复核</title>
<style>
body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;max-width:1500px;margin:28px auto;padding:0 18px;color:#202124;line-height:1.5}}
table{{border-collapse:collapse;width:100%;margin:18px 0 28px}}th,td{{border:1px solid #d9dce1;padding:8px;text-align:right}}th:first-child,td:first-child{{text-align:left}}th{{background:#f3f5f7}}
.note{{border-left:4px solid #2563eb;padding:10px 14px;background:#f6f8fa}}.grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px}}section{{min-width:0}}img{{width:100%;height:auto;border:1px solid #d9dce1}}h2{{font-size:18px}}code{{background:#f1f3f4;padding:2px 5px}}@media(max-width:900px){{.grid{{grid-template-columns:1fr}}}}
</style></head><body>
<h1>8 样本 Space Ranger 官方细胞分割复核</h1>
<p class="note">本阶段保留全部官方 filtered cells，只添加候选 QC 标记；尚未据此删除细胞。蓝绿色边界为随机抽样最多 30,000 个官方 cell polygons，用于检查 H&E 配准、组织覆盖和明显碎片化。</p>
<table><thead><tr><th>样本</th><th>细胞</th><th>中位 counts</th><th>中位基因</th><th>中位 mt%</th><th>中位面积 µm²</th><th>候选 QC 通过</th></tr></thead><tbody>{''.join(rows)}</tbody></table>
<div class="grid">{''.join(panels)}</div>
</body></html>"""
    output = reports / "official_segmented_qc_review.html"
    output.write_text(page)
    print(json.dumps({"status": "pass", "samples": len(summaries), "report": str(output)}, indent=2))


if __name__ == "__main__":
    main()
