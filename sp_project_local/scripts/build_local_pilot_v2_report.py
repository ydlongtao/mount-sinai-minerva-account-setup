#!/usr/bin/env python3
from __future__ import annotations

import csv
import html
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "local_results" / "SC000895-R4_v2_20260710"


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def table(rows: list[dict[str, str]]) -> str:
    if not rows:
        return "<p>无记录。</p>"
    columns = list(rows[0])
    head = "".join(f"<th>{html.escape(column)}</th>" for column in columns)
    body = []
    for row in rows:
        body.append("<tr>" + "".join(f"<td>{html.escape(str(row.get(column, '')))}</td>" for column in columns) + "</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def findings(qc: dict, model: dict, graph: list[dict], clusters: list[dict], markers: list[dict], h5ad: dict) -> list[str]:
    qc_fraction = qc["candidate_pass_fraction"] * 100
    mapped = model["mapping_008um"]["mapped_fraction"] * 100
    present = sum(row["present"].lower() == "true" for row in markers)
    total = len(markers)
    selected = model["selected_neighbors"]
    largest = graph[0].get("largest_component_fraction", "NA") if graph else "NA"
    return [
        f"16 µm 主分析保留 {model['shape_016um'][0]:,} 个 bins、{model['shape_016um'][1]:,} 个基因；QC 候选通过比例为 {qc_fraction:.2f}%。",
        f"CAGRA 在 n_neighbors={selected} 时通过连通性门槛，最大连通分量比例为 {largest}；该结果适合继续进行空间聚类解释。",
        f"16 µm 聚类标签已映射到 8 µm，映射覆盖率为 {mapped:.2f}%；8 µm 对象用于高分辨率空间定位，不作为本轮全局聚类依据。",
        f"marker panel 中 {present}/{total} 个 marker 在完整表达矩阵中可用；marker 结果仍属于候选区室提示，不应视为最终细胞类型注释。",
        f"最终 16 µm H5AD 的结构检查结果：{h5ad['files'][0]['shape'][0]:,} 个 observations，包含空间坐标、PCA、UMAP 和邻居图字段。",
        "8 µm H5AD 保留完整表达矩阵、空间坐标和映射标签，但没有独立 PCA、UMAP 或邻居图；后续若需要 8 µm 局部聚类，应单独设计局部建模流程。",
        "GPU 建模在 H100 上完成；GPU Leiden 因 Dask 依赖不匹配回退到 Scanpy igraph，因此当前结果是 GPU 邻居图/降维加 CPU Leiden 的混合流程。",
    ]


def main() -> None:
    qc = read_json(RESULTS / "qc" / "qc_summary.json")
    model = read_json(RESULTS / "model" / "model_summary.json")
    graph = read_json(RESULTS / "model" / "graph_diagnostics.json")
    clusters = read_csv(RESULTS / "model" / "clustering_resolution_summary.csv")
    markers = read_csv(RESULTS / "model" / "marker_availability.csv")
    h5ad = read_json(RESULTS / "h5ad_structure.json")
    paragraphs = "".join(f"<li>{html.escape(item)}</li>" for item in findings(qc, model, graph, clusters, markers, h5ad))
    report = f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><title>SC000895-R4 Visium HD v2 分析报告</title>
<style>
body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;max-width:1400px;margin:32px auto;padding:0 24px;color:#17202a;line-height:1.55}}
h1,h2{{color:#174a5b}} .note{{background:#fff7e6;border-left:4px solid #d97706;padding:12px 16px}}
table{{border-collapse:collapse;font-size:13px;margin:12px 0 28px;max-width:100%}}
th,td{{border:1px solid #c9d4d8;padding:5px 8px;text-align:left}} th{{background:#eaf1f3}}
img{{max-width:100%;height:auto;border:1px solid #d6dee2}} code,pre{{background:#f5f7f8;padding:2px 5px}}
</style></head><body>
<h1>SC000895-R4 Visium HD 第二版分析报告</h1>
<p>样本：人前列腺癌；主分辨率：16 µm；8 µm 用于标签映射和高分辨率定位。</p>
<h2>主要结论</h2><ul>{paragraphs}</ul>
<div class="note"><strong>解释边界：</strong>本报告是单样本 pilot 的 QC、图结构和候选 marker 审阅，不是最终病理区域或细胞类型注释。</div>
<h2>QC</h2>{table([qc])}
<h2>图连通性</h2>{table(graph)}
<h2>Leiden resolution 比较</h2>{table(clusters)}
<h2>Marker 可用性</h2>{table(markers)}
<h2>H5AD 结构检查</h2>{table(h5ad["files"])}
<h2>QC 图</h2>
<p><img src="qc/qc_distributions.png" alt="QC distributions"></p>
<p><img src="qc/counts_vs_genes_hexbin.png" alt="counts versus genes"></p>
<p><img src="qc/spatial_qc.png" alt="spatial QC"></p>
<h2>结果文件</h2>
<ul>
<li><code>model/SC000895-R4_016um_v2_full.h5ad</code></li>
<li><code>model/SC000895-R4_008um_v2_mapped.h5ad</code></li>
<li><code>pilot_v2_review.html</code>：服务器生成的原始 v2 表格报告</li>
</ul>
</body></html>
"""
    output = RESULTS / "SC000895-R4_v2_analysis_report.html"
    output.write_text(report)
    print(output)


if __name__ == "__main__":
    main()
