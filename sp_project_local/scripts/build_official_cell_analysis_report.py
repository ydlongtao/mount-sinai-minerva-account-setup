#!/usr/bin/env python3
"""Build the local HTML review report for the official cell-level analysis."""
from __future__ import annotations

import csv
import html
import json
from pathlib import Path


BASE = Path(__file__).resolve().parents[1]
RESULTS = BASE / "local_results" / "segmented_official_v1" / "analysis"
REPORT_DIR = RESULTS / "detailed_report"
SOURCE_DIR = REPORT_DIR / "sources"
SAMPLES = RESULTS / "samples"


def load_rows() -> list[dict]:
    rows = []
    for path in sorted(SAMPLES.glob("*/cell_level_analysis_summary.json")):
        data = json.loads(path.read_text())
        data["n_clusters"] = len(data["clusters"])
        data["zero_rate"] = data["zero_counts_retained"] / data["cells_total"]
        largest = max(data["clusters"].values())
        data["largest_cluster"] = largest
        data["largest_cluster_share"] = largest / data["cells_nonzero_graph"]
        data["marker_groups_present"] = sum(bool(v) for v in data.get("markers_present", {}).values())
        rows.append(data)
    if len(rows) != 8:
        raise RuntimeError(f"Expected 8 sample summaries, found {len(rows)}")
    return rows


def write_sources(rows: list[dict]) -> Path:
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    source = SOURCE_DIR / "official_cell_level_summary.csv"
    fields = ["sample_id", "cells_total", "cells_nonzero_graph", "zero_counts_retained", "zero_rate", "genes_hvg", "n_clusters", "largest_cluster", "largest_cluster_share", "marker_groups_present"]
    with source.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({key: row[key] for key in fields} for row in rows)
    return source


def write_image_review(rows: list[dict]) -> Path:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    cards = []
    for row in rows:
        sample = row["sample_id"]
        fig = f"../samples/{sample}/figures/"
        cards.append(f"""
        <section><h2>{html.escape(sample)}</h2>
        <p>{row['cells_total']:,} cells; {row['n_clusters']} Leiden clusters; zero-count objects retained: {row['zero_counts_retained']:,}.</p>
        <div class="pair"><figure><figcaption>UMAP / Leiden</figcaption><img src="{fig}{sample}_cell_level_umap_leiden.png" alt="{sample} UMAP Leiden"></figure>
        <figure><figcaption>Spatial / Leiden</figcaption><img src="{fig}{sample}_cell_level_spatial_leiden.png" alt="{sample} spatial Leiden"></figure></div></section>""")
    output = REPORT_DIR / "cell_level_image_review.html"
    output.write_text("""<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>Official cell-level image review</title>
    <style>body{font:15px system-ui,sans-serif;max-width:1200px;margin:30px auto;padding:0 20px;color:#172033}section{border-top:1px solid #ccd3df;padding:20px 0}.pair{display:grid;grid-template-columns:1fr 1fr;gap:20px}figure{margin:0}figcaption{font-weight:600;margin-bottom:8px}img{max-width:100%;border:1px solid #ccd3df}@media(max-width:800px){.pair{grid-template-columns:1fr}}</style>
    <h1>8 个样本官方细胞级结果图像审阅</h1><p>UMAP 显示每样本 Leiden 结构，空间图显示相同 cluster 在组织坐标中的分布。cluster 不是最终细胞类型标签。</p>""" + "".join(cards) + "</html>\n")
    return output


def build_artifact(rows: list[dict], source: Path, image_review: Path) -> Path:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    total_cells = sum(r["cells_total"] for r in rows)
    total_nonzero = sum(r["cells_nonzero_graph"] for r in rows)
    total_zero = sum(r["zero_counts_retained"] for r in rows)
    median_clusters = sorted(r["n_clusters"] for r in rows)[3:5]
    median_clusters = sum(median_clusters) / 2
    source_id = "official_cell_level_summary"
    sample_qc = [{k: r[k] for k in ["sample_id", "cells_total", "cells_nonzero_graph", "zero_counts_retained", "zero_rate", "genes_hvg", "n_clusters", "largest_cluster", "largest_cluster_share", "marker_groups_present"]} for r in rows]
    depth = [{"sample_id": r["sample_id"], "metric": "total cells", "value": r["cells_total"]} for r in rows]
    cluster = [{"sample_id": r["sample_id"], "metric": "Leiden clusters", "value": r["n_clusters"]} for r in rows]
    zero = [{"sample_id": r["sample_id"], "metric": "zero-count share", "value": r["zero_rate"]} for r in rows]
    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1,
            "title": "8 个样本 Visium HD 细胞级分析结果与下一步计划",
            "description": "官方 Space Ranger 分割细胞的逐样本标准化、聚类、marker scoring 和空间分布审阅。",
            "generatedAt": "2026-07-21T00:00:00Z",
            "blocks": [
                {"id": "title", "type": "markdown", "body": "# 8 个样本 Visium HD 细胞级分析结果与下一步计划"},
                {"id": "summary", "type": "markdown", "sourceId": source_id, "body": f"## 技术摘要\n\n**8 个样本的逐样本 cell-level 首轮分析已全部完成。** 共保留 {total_cells:,} 个官方 Space Ranger 细胞；其中 {total_nonzero:,} 个非零 counts 对象进入 PCA、邻居图和 Leiden 聚类，{total_zero:,} 个零 counts 对象仍保存在输出 H5AD 中，仅从图计算中排除。\n\n**当前结果适合进入跨样本整合前的人工审阅。** 每个样本均已完成 normalize_total、log1p、2,000 HVG、PCA、15-neighbor graph、Leiden resolution 0.4、UMAP 和前列腺 marker scoring。样本级 cluster 数量的中位数为 {median_clusters:.0f}，但 cluster 不是最终细胞类型标签。\n\n**下一步不应直接把 cluster 当作细胞类型或空间域。** 需要先比较 marker score、cluster 的空间位置和样本组成，再进行跨样本整合与正式注释。"},
                {"id": "finding_cells", "type": "markdown", "sourceId": source_id, "body": "## 细胞数量和非零 counts 结果支持继续分析\n\n各样本均成功生成 cell-level H5AD 和空间聚类图。R2 的 218,052 个细胞中只有 11 个为零 counts；其余样本也保持了极高的非零比例。因此，零 counts 不是当前批次的主要质量瓶颈。\n\n该结果支持继续使用官方分割细胞作为主分析对象，同时保留零 counts 标记以便审计。低复杂度样本 R2、R4、R6 不因复杂度单独剔除，后续在样本/患者层面建模时纳入其测序深度和组织背景。"},
                {"id": "cells_chart", "type": "chart", "chartId": "cells_by_sample"},
                {"id": "cluster_finding", "type": "markdown", "sourceId": source_id, "body": "## Cluster 数量存在样本间差异，需先做 marker 与空间复核\n\nLeiden cluster 数量从约 50 到 340 不等，说明同一 resolution 在不同样本的局部密度和表达复杂度下产生了不同粒度。这个差异本身不能解释为生物学差异，也不能直接用于比较细胞类型丰度。\n\n下一轮审阅应重点看：cluster 是否由已知 marker 支持、是否集中在合理的腺体/间质/炎症区域、是否由极少数细胞组成，以及同一候选细胞类型是否跨样本重复出现。"},
                {"id": "cluster_chart", "type": "chart", "chartId": "clusters_by_sample"},
                {"id": "zero_finding", "type": "markdown", "sourceId": source_id, "body": "## 零 counts 对象比例很低，处理方式保持为图分析排除而非数据删除\n\n零 counts 对象在全部官方细胞中的比例很低；它们保留在 H5AD 中，便于回溯官方分割与矩阵匹配。它们不进入 PCA、邻居图、Leiden 和 UMAP，是计算稳定性的处理，不是额外 QC 阈值。"},
                {"id": "zero_chart", "type": "chart", "chartId": "zero_rate_by_sample"},
                {"id": "sample_table", "type": "table", "tableId": "sample_summary"},
                {"id": "methods", "type": "markdown", "body": "## 方法、数据范围和指标定义\n\n- **分析对象：** 8 个保留样本的官方 Space Ranger segmented-cell H5AD；原始整数矩阵保存在 `layers[\"counts\"]`。\n- **标准化：** `normalize_total(target_sum=1e4)` 后 `log1p`；未对 R2/R4/R6 设置特殊最小 counts/genes 阈值。\n- **特征和图：** 每样本选 2,000 个 HVG，PCA 最多 30 个成分，15 个近邻，Leiden resolution 0.4，后端为 `scanpy_igraph`。\n- **零 counts：** 保留在输出对象；只在构建图时排除。\n- **空间结果：** 使用官方 polygon 质心对应的 `obsm[\"spatial\"]` 坐标绘制 cluster 分布。\n- **marker scoring：** 使用 luminal、basal、tumor、immune、stromal、endothelial 候选 marker panel；评分用于候选注释，不等同于最终标签。"},
                {"id": "limitations", "type": "markdown", "body": "## 局限性和稳健性检查\n\n- 当前是逐样本分析，不是跨样本整合；cluster 编号和数量不能跨样本直接对应。\n- Leiden resolution 会影响 cluster 粒度，需要在整合前做 resolution sensitivity review。\n- marker score 不是差异表达检验，也没有替代病理专家或参考单细胞数据的人工确认。\n- 当前报告没有使用患者、区域和临床变量做正式差异分析；这些变量应在整合后以患者/样本为重复单位建模。\n- 空间通讯不能把每个 cell 当作独立生物学重复，必须进行样本级汇总或混合效应/置换框架。"},
                {"id": "next_steps", "type": "markdown", "body": "## 下一步计划\n\n1. **逐样本结果审阅：** 查看 8 个样本的 UMAP、空间 Leiden 图、marker score 和小 cluster；确认候选细胞类型与组织位置相符。\n2. **跨样本整合：** 以样本、患者、区域和临床分组为 metadata，使用可审计的 batch-aware HVG/PCA 整合；同时保留未整合表达用于差异分析。\n3. **细胞类型注释：** 采用 marker score、cluster marker、空间位置和外部参考联合确认，不把 Leiden domain 直接命名为细胞类型。\n4. **空间域分析：** 在确认细胞类型后，基于 8 um/官方细胞坐标构建空间邻域和组织域，并叠加 H&E 进行复核。\n5. **差异分析：** 按患者/样本作为生物学重复，比较区域和临床分组；低复杂度样本保留并纳入协变量或敏感性分析。\n6. **空间通讯：** 先确定配体-受体数据库、邻域半径和最小支持规则，再以样本级通讯分数进行统计；不把单个 cell/bin 当作独立重复。\n7. **最终报告：** 汇总细胞类型比例、空间域、marker、差异结果和空间通讯，并附完整参数与失败记录。"},
                {"id": "questions", "type": "markdown", "body": "## 当前需要确认的问题\n\n- 哪些 Leiden cluster 可以稳定映射到 luminal、basal、tumor、immune、stromal 和 endothelial 候选类型？\n- 同一候选类型在 8 个样本之间是否具有一致 marker 和空间位置？\n- 跨样本整合后是否仍保持 R2/R4/R6 的组织结构，而不是被 batch correction 过度抹平？\n- 正式空间通讯中，邻域半径、配体-受体数据库和最小样本支持数采用什么固定规则？"},
            ],
            "sources": [{"id": source_id, "label": "Official cell-level analysis summaries", "path": "sources/official_cell_level_summary.csv", "query": {"sql": "SELECT * FROM official_cell_level_summary ORDER BY sample_id", "description": "Reviewed sample-level cell analysis summaries", "engine": "CSV snapshot", "language": "sql", "tables_used": ["official_cell_level_summary"]}}],
            "cards": [
                {"id": "n_samples", "dataset": "headline", "metrics": [{"label": "完成样本", "field": "sample_count", "format": "number"}], "sourceId": source_id},
                {"id": "n_cells", "dataset": "headline", "metrics": [{"label": "官方细胞数", "field": "total_cells", "format": "compact"}], "sourceId": source_id},
                {"id": "nonzero_share", "dataset": "headline", "metrics": [{"label": "非零 counts 比例", "field": "nonzero_share", "format": "percent"}], "sourceId": source_id},
                {"id": "median_clusters", "dataset": "headline", "metrics": [{"label": "样本 cluster 数中位数", "field": "median_clusters", "format": "number"}], "sourceId": source_id},
            ],
            "charts": [
                {"id": "cells_by_sample", "type": "bar", "title": "官方细胞数（按样本）", "subtitle": "每个样本官方 Space Ranger segmented-cell 数量；柱形从零开始。", "dataset": "sample_qc", "encodings": {"x": {"field": "sample_id", "type": "nominal", "label": "样本"}, "y": {"field": "cells_total", "type": "quantitative", "label": "官方细胞"}}, "sourceId": source_id},
                {"id": "clusters_by_sample", "type": "bar", "title": "Leiden cluster 数量（按样本）", "subtitle": "固定 resolution 0.4；数量差异不能直接解释为细胞类型差异。", "dataset": "cluster_summary", "encodings": {"x": {"field": "sample_id", "type": "nominal", "label": "样本"}, "y": {"field": "value", "type": "quantitative", "label": "Leiden clusters"}}, "sourceId": source_id},
                {"id": "zero_rate_by_sample", "type": "bar", "title": "零 counts 比例（按样本）", "subtitle": "零 counts 对象保留在 H5AD，仅从图计算中排除。", "dataset": "zero_summary", "valueFormat": "percent", "unit": "%", "encodings": {"x": {"field": "sample_id", "type": "nominal", "label": "样本"}, "y": {"field": "value", "type": "quantitative", "label": "零 counts 比例", "unit": "%"}}, "sourceId": source_id},
            ],
            "tables": [{"id": "sample_summary", "title": "8 个样本 cell-level 结果汇总", "description": "官方细胞、图分析对象、HVG 数量、Leiden cluster 数量和零 counts 比例。", "dataset": "sample_qc", "columns": [{"field": "sample_id", "label": "样本", "type": "string"}, {"field": "cells_total", "label": "官方细胞", "type": "number"}, {"field": "cells_nonzero_graph", "label": "图分析细胞", "type": "number"}, {"field": "zero_counts_retained", "label": "零 counts", "type": "number"}, {"field": "zero_rate", "label": "零 counts 比例", "type": "percent"}, {"field": "n_clusters", "label": "Leiden clusters", "type": "number"}, {"field": "largest_cluster_share", "label": "最大 cluster 占比", "type": "percent"}], "defaultSort": {"field": "cells_total", "direction": "desc"}, "density": "spacious", "sourceId": source_id}],
        },
        "snapshot": {"version": 1, "status": "ready", "generatedAt": "2026-07-21T00:00:00Z", "datasets": {"headline": [{"sample_count": 8, "total_cells": total_cells, "nonzero_share": total_nonzero / total_cells, "median_clusters": median_clusters}], "sample_qc": sample_qc, "cluster_summary": cluster, "zero_summary": zero, "depth_metrics": depth}},
        "package_info": {"report_notes": "Image figures are provided in the companion local image review page."},
    }
    path = REPORT_DIR / "artifact.json"
    path.write_text(json.dumps(artifact, ensure_ascii=False, indent=2) + "\n")
    (REPORT_DIR / "report_notes.json").write_text(json.dumps({"source": str(source), "image_review": str(image_review), "samples": [r["sample_id"] for r in rows]}, ensure_ascii=False, indent=2) + "\n")
    return path


def main() -> None:
    rows = load_rows()
    source = write_sources(rows)
    image_review = write_image_review(rows)
    print(build_artifact(rows, source, image_review))


if __name__ == "__main__":
    main()
