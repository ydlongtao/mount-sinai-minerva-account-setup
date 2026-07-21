#!/usr/bin/env python3
from __future__ import annotations

import base64
import csv
import io
import json
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from PIL import Image, ImageDraw, ImageOps


PROJECT = Path(__file__).resolve().parents[1]
INPUT = PROJECT / "local_results" / "segmented_official_v1"
OUTPUT = INPUT / "detailed_report"
TITLE = "8 样本 Visium HD 官方细胞分割：结果与技术审阅"


def image_data_uri(path: Path, max_width: int, quality: int = 76) -> str:
    # 使用 macOS 自带 sips 压缩审阅图，避免要求本地 Python 额外安装 Pillow。
    with tempfile.NamedTemporaryFile(suffix=".jpg") as handle:
        subprocess.run(
            ["sips", "-s", "format", "jpeg", "-s", "formatOptions", str(quality),
             "-Z", str(max_width), str(path), "--out", handle.name],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )
        data = Path(handle.name).read_bytes()
    return "data:image/jpeg;base64," + base64.b64encode(data).decode("ascii")


def contact_sheet_data_uri(items: list[tuple[str, Path]], cell_width: int, cell_height: int, quality: int = 65) -> str:
    columns = 2
    rows = (len(items) + columns - 1) // columns
    label_height = 34
    sheet = Image.new("RGB", (columns * cell_width, rows * (cell_height + label_height)), "white")
    draw = ImageDraw.Draw(sheet)
    for index, (label, path) in enumerate(items):
        with Image.open(path) as source:
            source = ImageOps.contain(source.convert("RGB"), (cell_width - 12, cell_height - 12), Image.Resampling.LANCZOS)
        col, row = index % columns, index // columns
        x = col * cell_width + (cell_width - source.width) // 2
        y = row * (cell_height + label_height) + label_height + (cell_height - source.height) // 2
        sheet.paste(source, (x, y))
        draw.text((col * cell_width + 12, row * (cell_height + label_height) + 9), label, fill="#202124")
    buffer = io.BytesIO()
    sheet.save(buffer, format="JPEG", quality=quality, optimize=True, progressive=True)
    return "data:image/jpeg;base64," + base64.b64encode(buffer.getvalue()).decode("ascii")


def fmt1(value: float) -> float:
    return round(float(value), 1)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    summary_csv = INPUT / "reports" / "official_segmented_qc_summary.csv"
    rows = list(csv.DictReader(summary_csv.open(newline="")))
    if len(rows) != 8:
        raise ValueError(f"Expected 8 sample rows, found {len(rows)}")

    sample_rows = []
    total_cells = 0
    total_pass = 0
    for row in rows:
        item = {
            "sample_id": row["sample_id"],
            "cells": int(row["cells"]),
            "genes": int(row["genes"]),
            "polygon_features": int(row["polygon_features"]),
            "matched_polygon_percent": float(row["matched_polygon_percent"]),
            "median_total_counts": float(row["median_total_counts"]),
            "median_genes": float(row["median_genes"]),
            "median_pct_mt": float(row["median_pct_mt"]),
            "median_cell_area_um2": float(row["median_cell_area_um2"]),
            "passes_qc_candidate": int(row["passes_qc_candidate"]),
            "passes_qc_candidate_percent": float(row["passes_qc_candidate_percent"]),
        }
        total_cells += item["cells"]
        total_pass += item["passes_qc_candidate"]
        sample_rows.append(item)

    weighted_qc = total_pass / total_cells
    depth_rows = []
    for item in sample_rows:
        depth_rows.extend([
            {"sample_id": item["sample_id"], "metric": "Median counts", "value": item["median_total_counts"], "cells": item["cells"]},
            {"sample_id": item["sample_id"], "metric": "Median genes", "value": item["median_genes"], "cells": item["cells"]},
        ])
    generated = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    source_id = "official_segmented_qc"
    headline = [{
        "sample_count": 8,
        "total_cells": total_cells,
        "gene_count": 18085,
        "polygon_match_fraction": 1.0,
        "candidate_qc_fraction": weighted_qc,
    }]

    cards = [
        {"id": "samples", "dataset": "headline", "sourceId": source_id, "description": "保留并进入正式分析的 Visium HD 样本数。", "metrics": [{"label": "正式样本", "field": "sample_count", "format": "number"}]},
        {"id": "cells", "dataset": "headline", "sourceId": source_id, "description": "8 个官方 filtered feature-cell matrices 中的细胞总数。", "metrics": [{"label": "官方细胞", "field": "total_cells", "format": "compact"}]},
        {"id": "match", "dataset": "headline", "sourceId": source_id, "description": "矩阵 cell IDs 能在官方 cell segmentation polygons 中找到的比例。", "metrics": [{"label": "细胞边界匹配", "field": "polygon_match_fraction", "format": "percent"}]},
        {"id": "qc", "dataset": "headline", "sourceId": source_id, "description": "按细胞数加权的候选 QC 通过率；该标记尚未用于删除细胞。", "metrics": [{"label": "候选 QC 通过", "field": "candidate_qc_fraction", "format": "percent"}]},
        {"id": "genes", "dataset": "headline", "sourceId": source_id, "description": "每个样本官方矩阵中的基因特征数。", "metrics": [{"label": "基因特征", "field": "gene_count", "format": "compact"}]},
    ]

    charts = [
        {
            "id": "cells_by_sample", "title": "官方分割细胞数", "subtitle": "8 个保留样本；每行对应一个 filtered feature-cell matrix",
            "type": "bar", "dataset": "sample_qc", "sourceId": source_id, "valueFormat": "compact",
            "encodings": {"x": {"field": "sample_id", "type": "nominal", "label": "Sample"}, "y": {"field": "cells", "type": "quantitative", "label": "Cells"}, "tooltip": [{"field": "polygon_features", "type": "quantitative", "label": "GeoJSON polygons"}]},
        },
        {
            "id": "depth_by_sample", "title": "单细胞表达复杂度", "subtitle": "各样本单细胞总 counts 与检出基因数的中位数",
            "type": "bar", "dataset": "depth_metrics", "sourceId": source_id, "valueFormat": "number",
            "encodings": {"x": {"field": "sample_id", "type": "nominal", "label": "Sample"}, "y": {"field": "value", "type": "quantitative", "label": "Median per cell"}, "color": {"field": "metric", "type": "nominal", "label": "Metric"}, "tooltip": [{"field": "cells", "type": "quantitative", "label": "Cells"}]},
        },
        {
            "id": "qc_by_sample", "title": "候选 QC 通过率", "subtitle": "counts、genes、mt% 与细胞面积四项条件的交集；未用于实际过滤",
            "type": "bar", "dataset": "sample_qc", "sourceId": source_id, "valueFormat": "number", "unit": "%",
            "encodings": {"x": {"field": "sample_id", "type": "nominal", "label": "Sample"}, "y": {"field": "passes_qc_candidate_percent", "type": "quantitative", "label": "Candidate QC pass", "unit": "%"}, "tooltip": [{"field": "passes_qc_candidate", "type": "quantitative", "label": "Passing cells"}, {"field": "cells", "type": "quantitative", "label": "All cells"}]},
        },
    ]

    tables = [{
        "id": "expression_qc_table", "title": "逐样本表达与候选 QC", "subtitle": "全部 8 个保留样本；精确值用于阈值审阅",
        "dataset": "sample_qc", "sourceId": source_id, "density": "spacious",
        "defaultSort": {"field": "sample_id", "direction": "asc"},
        "columns": [
            {"field": "sample_id", "label": "样本", "type": "text"},
            {"field": "cells", "label": "细胞", "format": "number"},
            {"field": "median_total_counts", "label": "中位 counts", "format": "number"},
            {"field": "median_genes", "label": "中位 genes", "format": "number"},
            {"field": "median_pct_mt", "label": "中位 mt%", "format": "number"},
            {"field": "passes_qc_candidate_percent", "label": "候选 QC %", "format": "number"},
        ],
    }, {
        "id": "segmentation_qc_table", "title": "逐样本分割与形态指标", "subtitle": "官方矩阵、GeoJSON 边界匹配与细胞面积",
        "dataset": "sample_qc", "sourceId": source_id, "density": "spacious",
        "defaultSort": {"field": "sample_id", "direction": "asc"},
        "columns": [
            {"field": "sample_id", "label": "样本", "type": "text"},
            {"field": "cells", "label": "矩阵细胞", "format": "number"},
            {"field": "polygon_features", "label": "Polygons", "format": "number"},
            {"field": "matched_polygon_percent", "label": "匹配率 %", "format": "number"},
            {"field": "median_cell_area_um2", "label": "中位面积 µm²", "format": "number"},
        ],
    }]

    overlay_items = []
    qc_items = []
    sample_notes = []
    for item in sample_rows:
        sample = item["sample_id"]
        fig_dir = INPUT / "samples" / sample / "figures"
        status = "重点复核" if item["passes_qc_candidate_percent"] < 90 else "常规复核"
        overlay_items.append((sample, fig_dir / "official_cell_boundaries_he_overview.png"))
        qc_items.append((sample, fig_dir / "official_cell_qc_distributions.png"))
        sample_notes.append(f"<li><strong>{sample} · {status}</strong>：细胞 {item['cells']:,}；中位 counts/genes {item['median_total_counts']:.0f}/{item['median_genes']:.0f}；候选 QC {item['passes_qc_candidate_percent']:.1f}%。</li>")
    # 阅读器正文宽度约 760 px；联系表保持在 720 px 内，避免横向溢出。
    visual_review_markdown = "## 逐样本 H&E 边界与 QC 分布\n\n完整分辨率图像集中在伴随审阅页：[打开 8 样本图像复核](image_review.html)。主报告保留定量比较与方法结论，图像页保留每个样本的官方 polygon overlay 和四项 QC 分布。\n\n" + "\n".join(
        f"- **{item['sample_id']}**：{('重点复核' if item['passes_qc_candidate_percent'] < 90 else '常规复核')}；候选 QC {item['passes_qc_candidate_percent']:.1f}%"
        for item in sample_rows
    )

    blocks = [
        {"id": "title", "type": "markdown", "body": f"# {TITLE}"},
        {"id": "technical_summary", "type": "markdown", "sourceId": source_id, "body": "## 技术摘要\n\n**官方 Space Ranger 细胞分割可作为后续细胞级主分析输入。** 8 个样本共包含 1,402,720 个细胞，表达矩阵中的 cell IDs 与官方 cell polygons 匹配率均为 100%，并已在 H5AD 中保存原始整数 `counts` layer 和官方 polygon 质心。\n\n**样本间测序复杂度差异明显。** R4 的中位 counts/genes 为 39/38，候选 QC 通过率 79.6%，是当前最需要重点复核的样本；R5 的对应值为 143/129 和 98.5%。这属于描述性 QC 差异，尚不能解释为生物学差异。\n\n**当前没有删除细胞。** 候选 QC 标记仅用于阈值审阅；在确认 H&E 边界覆盖和低深度样本后，再固定正式筛选规则。"},
        {"id": "headline_metrics", "type": "metric-strip", "cardIds": ["samples", "cells", "match", "qc", "genes"]},
        {"id": "finding_cells", "type": "markdown", "sourceId": source_id, "body": "## 官方细胞产量覆盖全部 8 个样本\n\nR6 的细胞数最低（95,354），R5 最高（253,596），约相差 2.7 倍。细胞产量受组织面积、细胞密度、切片完整性和分割结果共同影响，不能单独作为样本质量排序。"},
        {"id": "chart_cells", "type": "chart", "chartId": "cells_by_sample"},
        {"id": "finding_depth", "type": "markdown", "sourceId": source_id, "body": "## R4、R2 与 R6 的表达复杂度偏低\n\nR4 的中位 counts 和检出基因最低，其次是 R2；R6 的深度高于这两例，但候选 QC 通过率仍偏低，提示细胞面积或低计数尾部也在影响结果。正式聚类前应采用逐样本 QC，而不是一个未经检查的全局硬阈值。"},
        {"id": "chart_depth", "type": "chart", "chartId": "depth_by_sample"},
        {"id": "finding_qc", "type": "markdown", "sourceId": source_id, "body": f"## 候选 QC 保留总体较高，但低深度样本需要单独审阅\n\n按细胞数加权，候选 QC 通过率为 {100 * weighted_qc:.1f}%。R4 为 79.6%，R6 为 84.2%，R2 为 88.1%；其余样本均超过 90%。该差异由四项规则共同产生，报告不将未通过细胞直接等同于低质量细胞。"},
        {"id": "chart_qc", "type": "chart", "chartId": "qc_by_sample"},
        {"id": "exact_table_intro", "type": "markdown", "body": "## 精确值支持阈值复核\n\n下表保留每个样本的精确分母和中位指标。图表用于比较形状，表格用于审核具体阈值。"},
        {"id": "expression_table", "type": "table", "tableId": "expression_qc_table"},
        {"id": "segmentation_table", "type": "table", "tableId": "segmentation_qc_table"},
        {"id": "image_review", "type": "markdown", "sourceId": source_id, "body": visual_review_markdown},
        {"id": "scope", "type": "markdown", "sourceId": source_id, "body": "## 分析范围与指标定义\n\n- **分析队列：** R1、R2、R4、R5、R6、R7、R8、R9；旧版 Visium 样本 B408/B573 已排除。\n- **细胞：** `segmented_outputs/filtered_feature_cell_matrix.h5` 中的官方 filtered cell barcode。\n- **边界匹配率：** 矩阵 cell ID 能在 `cell_segmentations.geojson` 中找到 polygon 的比例。\n- **候选 QC：** `total_counts ≥ 20`、`n_genes ≥ 10`、`mt% ≤ 30`、`10 ≤ cell_area ≤ 1000 µm²` 四项同时满足。\n- **比较基础：** 当前批次 8 个样本的描述性横向比较，无条件组、时间点或患者级推断。"},
        {"id": "method", "type": "markdown", "body": "## 数据处理与空间映射方法\n\n官方 H5 被读取为稀疏原始计数矩阵；重复 gene symbols 经唯一化后，原始整数矩阵复制到 `layers['counts']`。细胞质心与面积由官方 GeoJSON polygon 计算，质心保存为 full-resolution pixel 坐标；H&E 叠加使用每个样本的官方 `tissue_hires_scalef`，并保持图像坐标的 y 轴方向。Marker score 已基于每细胞 library-size 标准化后的指定前列腺癌 marker 集计算，但本报告不据此赋予最终细胞类型。"},
        {"id": "limitations", "type": "markdown", "body": "## 限制、不确定性与稳健性检查\n\n- 100% ID 匹配证明矩阵和边界可连接，但不证明每个细胞边界在病理学上都正确。\n- 全景 overlay 为确定性随机抽样，不显示全部 polygon；局部边界精度仍需在高倍区域检查。\n- 候选 QC 阈值是操作性起点，未按样本深度或细胞类型调整。低 RNA 含量的免疫细胞可能更容易被误删。\n- 当前结果是描述性 QC，不包含聚类、批次校正、细胞类型推断、差异表达或细胞通讯结论。\n- 细胞不是独立生物学重复；后续统计必须以样本或患者为重复单位。"},
        {"id": "next_steps", "type": "markdown", "body": "## 推荐的下一阶段\n\n1. 人工确认 8 张 H&E overlay，重点检查 R4、R2 和 R6。\n2. 基于每样本分布固定正式 QC 规则，并保留 `qc_reason`，避免静默删除。\n3. 在 QC 通过细胞上进行单样本 PCA/邻接图/Leiden 与 marker 初注释。\n4. 将细胞类型和主要 marker 映射回 H&E，并与现有 8 µm 空间域交叉验证。\n5. 通过样本级伪批量和留一法检查整合稳定性，再进入空间邻域与配体–受体通讯。"},
        {"id": "questions", "type": "markdown", "body": "## 需要在进入聚类前确认的问题\n\n- R4、R2、R6 的低复杂度是否与病理区域、组织质量或测序深度一致？\n- 官方边界是否在腺体、间质和炎症密集区都保持合理覆盖？\n- 是否需要对低 RNA 免疫细胞采用不同的最小 counts/genes 阈值？\n- 后续差异分析可用的患者、区域和临床分组信息是否完整？"},
    ]

    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": TITLE,
            "description": "8 个 Visium HD 前列腺癌样本的 Space Ranger 官方细胞分割、表达复杂度、候选 QC 与 H&E 边界技术审阅。",
            "generatedAt": generated,
            "cards": cards,
            "charts": charts,
            "tables": tables,
            "sources": [{"id": source_id, "label": "Official segmented-cell QC summary", "path": "sources/official_segmented_qc_summary.csv"}],
            "blocks": blocks,
        },
        "snapshot": {"version": 1, "generatedAt": generated, "status": "ready", "datasets": {"headline": headline, "sample_qc": sample_rows, "depth_metrics": depth_rows}},
        "sources": [{
            "id": source_id,
            "label": "Official segmented-cell QC summary",
            "path": "sources/official_segmented_qc_summary.csv",
            "query": {
                "engine": "sqlite",
                "language": "sql",
                "sql": "SELECT sample_id, cells, genes, polygon_features, matched_polygon_percent, median_total_counts, median_genes, median_pct_mt, median_cell_area_um2, passes_qc_candidate, passes_qc_candidate_percent FROM official_segmented_qc_summary ORDER BY sample_id",
                "description": "Aggregated from eight Space Ranger filtered feature-cell matrices and matched cell segmentation GeoJSON files.",
                "tables_used": ["official_segmented_qc_summary"],
                "filters": ["Retained samples: R1, R2, R4, R5, R6, R7, R8, R9", "B408 and B573 excluded", "No cells deleted"],
                "metric_definitions": ["Candidate QC requires counts>=20, genes>=10, mt%<=30, and cell area 10-1000 um2", "Weighted QC pass rate is total passing cells divided by total official filtered cells"],
            },
        }],
    }
    artifact_path = OUTPUT / "artifact.json"
    artifact_path.write_text(json.dumps(artifact, ensure_ascii=False, indent=2) + "\n")
    notes = {
        "audience": "technical",
        "delivery_mode": "html",
        "chart_map": [
            {"section": "cell yield", "question": "How many official cells are available per sample?", "family": "comparison", "type": "horizontal bar", "dataset": "sample_qc", "fields": ["sample_id", "cells"]},
            {"section": "expression complexity", "question": "Which samples have lower median expression complexity?", "family": "comparison", "type": "grouped bar", "dataset": "sample_qc", "fields": ["sample_id", "median_total_counts", "median_genes"]},
            {"section": "candidate QC", "question": "Which samples require threshold review?", "family": "comparison", "type": "horizontal bar", "dataset": "sample_qc", "fields": ["sample_id", "passes_qc_candidate_percent"]},
        ],
        "omitted_visuals": [{"type": "scatter", "reason": "Only eight sample-level observations; a relationship chart would be underpowered and easy to over-interpret."}],
        "source": str(summary_csv),
    }
    (OUTPUT / "report_notes.json").write_text(json.dumps(notes, ensure_ascii=False, indent=2) + "\n")
    print(artifact_path)


if __name__ == "__main__":
    main()
