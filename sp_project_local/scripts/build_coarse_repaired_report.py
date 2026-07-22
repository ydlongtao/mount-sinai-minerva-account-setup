#!/usr/bin/env python3
"""Build a detailed local report for the final 20-30 cluster coarse results."""
from __future__ import annotations

import csv
import html
import json
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE / "local_results" / "segmented_official_v1"
COARSE = ROOT / "analysis_coarse" / "samples"
REPAIRED = ROOT / "analysis_coarse_repaired" / "samples"
OUT = ROOT / "coarse_repaired_report"
SRC = OUT / "sources"
SAMPLES = ["SC000895-R1", "SC000895-R2", "SC000895-R4", "SC000895-R5", "SC000895-R6", "SC000895-R7", "SC000895-R8", "SC000895-R9"]


def load() -> list[dict]:
    rows = []
    for sample in SAMPLES:
        repaired = REPAIRED / sample / "coarse_repaired_summary.json"
        if repaired.exists():
            d = json.loads(repaired.read_text())
            base = json.loads((COARSE / sample / "coarse_analysis_summary.json").read_text())
            for key in ("cells_total", "cells_nonzero_graph", "zero_counts_retained"):
                d.setdefault(key, base[key])
            d["result_version"] = "coarse_repaired_v2"
            d["umap"] = f"../analysis_coarse_repaired/samples/{sample}/figures/{sample}_repaired_umap_leiden.png"
            d["spatial"] = f"../analysis_coarse_repaired/samples/{sample}/figures/{sample}_repaired_spatial_leiden.png"
        else:
            d = json.loads((COARSE / sample / "coarse_analysis_summary.json").read_text())
            d["chosen_n_neighbors"] = 15
            d["result_version"] = "coarse_v2"
            d["umap"] = f"../analysis_coarse/samples/{sample}/figures/{sample}_coarse_umap_leiden.png"
            d["spatial"] = f"../analysis_coarse/samples/{sample}/figures/{sample}_coarse_spatial_leiden.png"
        rows.append(d)
    if len(rows) != 8 or any(not (20 <= r["chosen_cluster_count"] <= 30) for r in rows):
        raise RuntimeError("Final coarse results are incomplete or outside the 20-30 cluster target")
    return rows


def write_sources(rows: list[dict]) -> None:
    SRC.mkdir(parents=True, exist_ok=True)
    fields = ["sample_id", "result_version", "cells_total", "cells_nonzero_graph", "zero_counts_retained", "chosen_n_neighbors", "chosen_resolution", "chosen_cluster_count", "output_h5ad"]
    with (SRC / "final_coarse_cluster_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def write_images(rows: list[dict]) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    sections = []
    for row in rows:
        sections.append(f"<section><h2>{html.escape(row['sample_id'])}</h2><p>{row['chosen_cluster_count']} clusters; n_neighbors={row['chosen_n_neighbors']}; resolution={row['chosen_resolution']}; {html.escape(row['result_version'])}.</p><div class='pair'><figure><figcaption>UMAP / Leiden</figcaption><img src='{row['umap']}' alt='{row['sample_id']} UMAP Leiden'></figure><figure><figcaption>Spatial / Leiden</figcaption><img src='{row['spatial']}' alt='{row['sample_id']} spatial Leiden'></figure></div></section>")
    path = OUT / "final_coarse_image_review.html"
    path.write_text("<!doctype html><html lang='zh-CN'><meta charset='utf-8'><title>Final coarse Leiden image review</title><style>body{font:15px system-ui,sans-serif;max-width:1200px;margin:30px auto;padding:0 20px;color:#172033}section{border-top:1px solid #ccd3df;padding:20px 0}.pair{display:grid;grid-template-columns:1fr 1fr;gap:20px}figure{margin:0}img{max-width:100%;border:1px solid #ccd3df}@media(max-width:800px){.pair{grid-template-columns:1fr}}</style><h1>8 个样本最终 coarse Leiden 图像审阅</h1><p>所有样本均控制在 20–30 个 cluster；cluster 仍是表达结构标签，不等同于最终细胞类型。</p>" + "".join(sections) + "</html>\n")
    return path


def build(rows: list[dict]) -> Path:
    total = sum(r["cells_total"] for r in rows)
    nonzero = sum(r["cells_nonzero_graph"] for r in rows)
    source_id = "final_coarse_summary"
    sample_rows = [{"sample_id": r["sample_id"], "result_version": r["result_version"], "cells_total": r["cells_total"], "cells_nonzero_graph": r["cells_nonzero_graph"], "zero_counts_retained": r["zero_counts_retained"], "chosen_n_neighbors": r["chosen_n_neighbors"], "chosen_resolution": r["chosen_resolution"], "chosen_cluster_count": r["chosen_cluster_count"]} for r in rows]
    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1,
            "title": "Visium HD 8样本最终粗粒度聚类结果与分析计划",
            "description": "8 个官方 segmented-cell 样本的 20–30 cluster coarse Leiden 结果、参数修复和后续空间分析计划。",
            "generatedAt": "2026-07-22T00:00:00Z",
            "sources": [{"id": source_id, "label": "Final coarse Leiden summaries", "path": "sources/final_coarse_cluster_summary.csv", "query": {"sql": "SELECT * FROM final_coarse_cluster_summary ORDER BY sample_id", "description": "Final coarse per-sample Leiden summaries", "engine": "CSV snapshot", "language": "sql", "tables_used": ["final_coarse_cluster_summary"]}}],
            "cards": [
                {"id": "samples", "dataset": "headline", "metrics": [{"label": "完成样本", "field": "sample_count", "format": "number"}], "sourceId": source_id},
                {"id": "cells", "dataset": "headline", "metrics": [{"label": "官方细胞", "field": "total_cells", "format": "compact"}], "sourceId": source_id},
                {"id": "nonzero", "dataset": "headline", "metrics": [{"label": "非零 counts 比例", "field": "nonzero_share", "format": "percent"}], "sourceId": source_id},
                {"id": "target", "dataset": "headline", "metrics": [{"label": "目标范围内样本", "field": "target_samples", "format": "number"}], "sourceId": source_id},
            ],
            "charts": [
                {"id": "cluster_count", "type": "bar", "title": "最终 Leiden cluster 数量", "subtitle": "8 个样本均位于 20–30 的预设人工审阅范围内。", "dataset": "sample_summary", "encodings": {"x": {"field": "sample_id", "type": "nominal", "label": "样本"}, "y": {"field": "chosen_cluster_count", "type": "quantitative", "label": "clusters"}}, "sourceId": source_id},
                {"id": "resolution", "type": "bar", "title": "每样本选择的 Leiden resolution", "subtitle": "不同样本使用不同 resolution，以适应图结构差异；不能跨样本比较 cluster 编号。", "dataset": "sample_summary", "encodings": {"x": {"field": "sample_id", "type": "nominal", "label": "样本"}, "y": {"field": "chosen_resolution", "type": "quantitative", "label": "resolution"}}, "sourceId": source_id},
                {"id": "neighbors", "type": "bar", "title": "每样本使用的邻居数", "subtitle": "R2、R4、R6、R9 使用更大的邻居数以减少图碎片化。", "dataset": "sample_summary", "encodings": {"x": {"field": "sample_id", "type": "nominal", "label": "样本"}, "y": {"field": "chosen_n_neighbors", "type": "quantitative", "label": "n_neighbors"}}, "sourceId": source_id},
            ],
            "tables": [{"id": "summary", "title": "最终 coarse Leiden 参数与结果", "description": "每行是一个样本；v1/v2 结果均保留，表中为用于下一步人工审阅的最终版本。", "dataset": "sample_summary", "density": "spacious", "defaultSort": {"field": "sample_id", "direction": "asc"}, "columns": [{"field": "sample_id", "label": "样本", "type": "text"}, {"field": "result_version", "label": "结果版本", "type": "text"}, {"field": "cells_total", "label": "官方细胞", "format": "number"}, {"field": "chosen_n_neighbors", "label": "n_neighbors", "format": "number"}, {"field": "chosen_resolution", "label": "resolution", "format": "number"}, {"field": "chosen_cluster_count", "label": "clusters", "format": "number"} ], "sourceId": source_id}],
            "blocks": [
                {"id": "title", "type": "markdown", "body": "# Visium HD 8样本最终粗粒度聚类结果与分析计划"},
                {"id": "summary", "type": "markdown", "sourceId": source_id, "body": f"## 技术摘要\n\n**8 个样本的 coarse Leiden 结果均已进入 20–30 cluster 目标范围。** 共覆盖 {total:,} 个官方细胞，非零 counts 对象占 {nonzero / total:.2%}。\n\n**R2、R4、R6、R9 通过增大 n_neighbors 修复了图碎片化。** 最终参数分别为 R2: 80/0.005、R4: 50/0.002、R6: 80/0.002、R9: 50/0.015；R1、R5、R7、R8 使用原 coarse v2 结果。\n\n**该版本适合人工细胞类型和空间结构审阅。** cluster 数量被控制为可读范围，但 cluster 仍不是最终细胞类型，也不能把不同样本的 cluster 编号直接对应。"},
                {"id": "cluster_result", "type": "markdown", "sourceId": source_id, "body": "## 目标已实现，但参数不可跨样本直接比较\n\n所有样本最终 cluster 数为 22–30。R4 达到上限 30，R6 为 22，其余样本为 23–27。这个范围适合在 UMAP、空间图和 marker table 中做人工复核。\n\n需要注意，R2/R4/R6/R9 使用更大的 n_neighbors，说明这些样本原先的 15-neighbor graph 存在较强碎片化。它们的 cluster 粒度因此是参数适配后的结果，不应被解释为与其他样本具有相同的无监督结构。"},
                {"id": "cluster_chart", "type": "chart", "chartId": "cluster_count"},
                {"id": "parameter_result", "type": "markdown", "sourceId": source_id, "body": "## 参数修复保留了原始数据和旧版本\n\n本轮只改变邻居图和 Leiden 参数，未删除细胞、未修改 counts layer、未改变空间坐标。v1 高分辨率结果和 coarse v2 结果仍保留，可用于敏感性对照。\n\n选择不同 resolution 是为了达到预设的审阅粒度，不是为了制造预先指定的生物学类别。正式细胞类型仍需 marker 和空间证据支持。"},
                {"id": "resolution_chart", "type": "chart", "chartId": "resolution"},
                {"id": "neighbor_chart", "type": "chart", "chartId": "neighbors"},
                {"id": "table", "type": "table", "tableId": "summary"},
                {"id": "methods", "type": "markdown", "body": "## 方法和结果版本定义\n\n- **输入：** 官方 segmented-cell H5AD；整数原始矩阵保留在 `layers[\"counts\"]`。\n- **图构建：** 2,000 HVG、最多 30 个 PCA 成分；固定 metric 和 random seed；邻居数按样本在 15、30、50、80 中搜索。\n- **Leiden：** `scanpy_igraph`，resolution 在 0.4 到 0.001 的候选集合中搜索，优先选择 20–30 个 cluster 且最接近 25 的组合。\n- **零 counts：** 保留在 H5AD，仅从 PCA/邻居图/Leiden 中排除。\n- **版本：** v1 为高分辨率探索；coarse v2 为首次粗粒度结果；repaired v2 仅覆盖 R2、R4、R6、R9。"},
                {"id": "limitations", "type": "markdown", "body": "## 局限性与稳健性检查\n\n- 20–30 是人工审阅的显示粒度，不是统计学或病理学的真实细胞类型数。\n- 不同样本使用不同 n_neighbors/resolution，因此 cluster 编号、cluster 数和 cluster 大小不能直接做跨样本差异比较。\n- 当前没有进行 batch correction、整合、差异表达、细胞类型最终注释或空间通讯。\n- 需要检查 R4/R6 在增大邻居数后是否过度合并局部组织结构；这一步必须通过 UMAP、空间图和 marker 共同判断。\n- 细胞不是独立生物学重复；正式统计应以患者/样本为重复单位。"},
                {"id": "next", "type": "markdown", "body": "## 下一步分析计划\n\n1. **图像与 marker 审阅：** 查看 8 个样本的最终 coarse UMAP、空间 Leiden 图和 marker scoring；重点检查 R4/R6 是否出现过度合并。\n2. **候选细胞类型注释：** 对 luminal、basal、tumor、immune、stromal、endothelial 等候选类型，结合 marker、差异基因、空间位置和 H&E 进行确认。\n3. **跨样本整合：** 保留 raw counts，使用样本/患者/区域 metadata 进行 batch-aware PCA/整合；比较整合前后结构，避免过度校正。\n4. **整合 cluster 与细胞类型审阅：** 以 20–30 个 coarse cluster 作为人工审阅入口，但建立跨样本稳定的细胞类型标签。\n5. **空间域与邻域：** 在官方细胞坐标和 8 um bin 结果之间交叉验证空间域，叠加 H&E 检查腺体、间质和炎症区域。\n6. **差异分析：** 以患者/样本为生物学重复，针对临床和区域分组开展 pseudobulk 或 mixed-effects 分析；R2/R4/R6 不因低复杂度单独剔除。\n7. **空间通讯：** 固定 ligand-receptor 数据库、邻域半径和最小支持数，以样本级通讯分数进行统计，禁止把 cell/bin 当作独立重复。\n8. **最终报告：** 输出细胞类型组成、空间域、差异表达、通路和空间通讯结果，并记录所有参数敏感性分析。"},
                {"id": "images", "type": "markdown", "body": "## 图像审阅\n\n逐样本 UMAP 和空间图位于配套页面：[打开最终 coarse 图像审阅](final_coarse_image_review.html)。图像页用于确认参数调整后的空间结构是否仍符合 H&E 和组织学预期。"},
                {"id": "questions", "type": "markdown", "body": "## 进入整合前需要确认的问题\n\n- R4/R6 增大 n_neighbors 后是否把相邻但生物学不同的区域过度合并？\n- 哪些 coarse cluster 可以被稳定映射为细胞类型，而不是纯粹的深度/批次效应？\n- 患者、区域和临床分组在整合模型中使用哪些固定协变量？\n- 空间通讯采用何种 ligand-receptor 数据库、邻域半径和样本级支持阈值？"},
            ],
        },
        "snapshot": {"version": 1, "status": "ready", "generatedAt": "2026-07-22T00:00:00Z", "datasets": {"headline": [{"sample_count": 8, "total_cells": total, "nonzero_share": nonzero / total, "target_samples": 8}], "sample_summary": sample_rows}},
        "package_info": {"report_notes": "Figures are available in the companion local image-review HTML."},
    }
    path = OUT / "artifact.json"
    path.write_text(json.dumps(artifact, ensure_ascii=False, indent=2) + "\n")
    (OUT / "report_notes.json").write_text(json.dumps({"sample_count": 8, "source_csv": "sources/final_coarse_cluster_summary.csv", "all_samples_in_target": True}, ensure_ascii=False, indent=2) + "\n")
    return path


def main() -> None:
    rows = load()
    OUT.mkdir(parents=True, exist_ok=True)
    write_sources(rows)
    write_images(rows)
    print(build(rows))


if __name__ == "__main__":
    main()
