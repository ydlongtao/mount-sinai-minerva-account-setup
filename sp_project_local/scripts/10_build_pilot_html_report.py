#!/usr/bin/env python3
from __future__ import annotations

import argparse
import base64
import csv
import html
import json
import re
from pathlib import Path


def source_value(value: str, tooltip_id: str, source_file: str) -> str:
    return (
        f'<span class="source-tooltip" tabindex="0" aria-describedby="{tooltip_id}">{html.escape(value)}'
        f'<span class="source-tooltip-content" id="{tooltip_id}" role="tooltip">'
        f'Source: SC000895-R4 GPU pilot review<br>File: {html.escape(source_file)}</span></span>'
    )


def image_uri(path: Path) -> str:
    mime = "image/png" if path.suffix.lower() == ".png" else "image/jpeg"
    return f"data:{mime};base64,{base64.b64encode(path.read_bytes()).decode('ascii')}"


def horizontal_bar_svg(rows: list[dict[str, object]], label_key: str, value_key: str, color: str) -> str:
    width, height, left, right, top = 960, 430, 150, 70, 30
    plot_width = width - left - right
    max_value = max(float(row[value_key]) for row in rows)
    bar_height, gap = 24, 12
    parts = [f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="Horizontal bar chart">']
    for index, row in enumerate(rows):
        y = top + index * (bar_height + gap)
        value = float(row[value_key])
        bar_width = plot_width * value / max_value
        parts.append(f'<text x="{left - 12}" y="{y + 17}" text-anchor="end" fill="currentColor" font-size="13">{html.escape(str(row[label_key]))}</text>')
        parts.append(f'<rect x="{left}" y="{y}" width="{bar_width:.1f}" height="{bar_height}" rx="3" fill="{color}"/>')
        parts.append(f'<text x="{min(left + bar_width + 8, width - right)}" y="{y + 17}" fill="currentColor" font-size="12">{value:,.0f}</text>')
    parts.append('</svg>')
    return "".join(parts)


def qc_line_svg(rows: list[dict[str, object]]) -> str:
    width, height, left, right, top, bottom = 960, 380, 80, 40, 28, 58
    plot_width, plot_height = width - left - right, height - top - bottom
    max_value = max(float(row["value"]) for row in rows) * 1.08
    colors = {"Total counts": "#0169cc", "Detected genes": "#e25507"}
    grouped: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault(str(row["metric"]), []).append(row)
    parts = [f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="QC percentile curves">']
    for tick in [0, 25, 50, 75, 100]:
        x = left + plot_width * tick / 100
        parts.append(f'<line x1="{x}" y1="{top}" x2="{x}" y2="{top + plot_height}" stroke="var(--grid)"/>')
        parts.append(f'<text x="{x}" y="{height - 24}" text-anchor="middle" fill="currentColor" font-size="12">P{tick}</text>')
    for tick in [0, 20, 40, 60, 80]:
        y = top + plot_height * (1 - tick / max_value)
        parts.append(f'<line x1="{left}" y1="{y}" x2="{left + plot_width}" y2="{y}" stroke="var(--grid)"/>')
        parts.append(f'<text x="{left - 12}" y="{y + 4}" text-anchor="end" fill="currentColor" font-size="12">{tick}</text>')
    for metric, metric_rows in grouped.items():
        points = []
        for row in metric_rows:
            x = left + plot_width * float(row["percentile"]) / 100
            y = top + plot_height * (1 - float(row["value"]) / max_value)
            points.append((x, y))
        point_text = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
        parts.append(f'<polyline points="{point_text}" fill="none" stroke="{colors[metric]}" stroke-width="3"/>')
        for x, y in points:
            parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{colors[metric]}"/>')
    parts.append('<text x="80" y="16" fill="#0169cc" font-size="12">Total counts</text>')
    parts.append('<text x="190" y="16" fill="#e25507" font-size="12">Detected genes</text>')
    parts.append('</svg>')
    return "".join(parts)


def chart_figure(chart_id: str, title: str, subtitle: str, fallback: str, source_file: str, note: str) -> str:
    tooltip_id = f"{chart_id}-source"
    return f'''<div class="wide"><figure class="card source-figure">
      <div class="card-head"><h3>{html.escape(title)}</h3><p>{html.escape(subtitle)}</p></div>
      <div class="chart-wrap"><div data-recharts-chart="{chart_id}">
        <div class="chart-fallback" data-recharts-fallback>{fallback}</div>
        <div data-recharts-live aria-hidden="true"></div>
      </div></div>
      <figcaption class="chart-note">{html.escape(note)}</figcaption>
      <button type="button" class="source-tooltip" aria-describedby="{tooltip_id}">Source<span class="source-tooltip-content" id="{tooltip_id}" role="tooltip">Source: SC000895-R4 GPU pilot review<br>File: {html.escape(source_file)}</span></button>
    </figure></div>'''


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--shell", type=Path, required=True)
    parser.add_argument("--output-shell", type=Path, required=True)
    parser.add_argument("--payload", type=Path, required=True)
    args = parser.parse_args()

    metrics = json.loads((args.input_dir / "review" / "review_metrics.json").read_text())
    with (args.input_dir / "review" / "cluster_sizes.csv").open(newline="") as handle:
        clusters = list(csv.DictReader(handle))
    with (args.input_dir / "review" / "qc_quantiles.csv").open(newline="") as handle:
        qc_rows = list(csv.DictReader(handle))

    top_clusters = [{"cluster": f"Cluster {row['leiden']}", "n_bins": int(row["n_bins"])} for row in clusters[:10]]
    stage_labels = {
        "008um:neighbors": "CAGRA neighbors",
        "008um:umap": "GPU UMAP",
        "008um:leiden": "CPU Leiden",
        "008um:plots_final_write": "Plots + H5AD write",
    }
    stage_rows = [
        {"stage": label, "minutes": round(metrics["stages"][key]["seconds"] / 60, 2)}
        for key, label in stage_labels.items()
    ]
    selected_quantiles = {"0.01", "0.05", "0.25", "0.5", "0.75", "0.95", "0.99"}
    qc_chart_rows = []
    for row in qc_rows:
        if row["quantile"] not in selected_quantiles:
            continue
        percentile = float(row["quantile"]) * 100
        qc_chart_rows.extend([
            {"percentile": percentile, "metric": "Total counts", "value": float(row["total_counts"])},
            {"percentile": percentile, "metric": "Detected genes", "value": float(row["n_genes_by_counts"])},
        ])

    top_chart = chart_figure(
        "top-clusters",
        "Top 10 Leiden clusters",
        "Exact bin counts; the long tail contains 1,442 additional clusters",
        horizontal_bar_svg(top_clusters, "cluster", "n_bins", "#0169cc"),
        "review/cluster_sizes.csv",
        "The two largest clusters contain 35.2% of bins, while hundreds of small clusters create an unusable long tail.",
    )
    qc_chart = chart_figure(
        "qc-percentiles",
        "008um bin complexity by percentile",
        "Counts and detected genes per retained nonzero bin",
        qc_line_svg(qc_chart_rows),
        "review/qc_quantiles.csv",
        "Median complexity is only 12 counts and 11 detected genes, so bin-level clustering is dominated by sparse observations.",
    )
    runtime_chart = chart_figure(
        "stage-runtime",
        "Runtime by major stage",
        "Minutes; checkpoint loading is omitted because it took under three seconds",
        horizontal_bar_svg(stage_rows, "stage", "minutes", "#8046d9"),
        "results/pilot_gpu_report.json",
        "CPU Leiden consumed 84.14 minutes and dominates runtime after GPU Leiden failed on the incompatible Dask stack.",
    )

    src_metrics = "review/review_metrics.json"
    src_qc = "review/qc_quantiles.csv"
    s = lambda value, key, source=src_metrics: source_value(str(value), key, source)

    qc_table_rows = []
    quantile_labels = {"0.25": "P25", "0.5": "P50", "0.75": "P75", "0.95": "P95", "0.99": "P99"}
    for row in qc_rows:
        if row["quantile"] not in quantile_labels:
            continue
        key = row["quantile"].replace(".", "-")
        qc_table_rows.append(
            f'<tr><td>{quantile_labels[row["quantile"]]}</td><td>{s(f"{float(row["total_counts"]):.0f}", f"qc-counts-{key}", src_qc)}</td>'
            f'<td>{s(f"{float(row["n_genes_by_counts"]):.0f}", f"qc-genes-{key}", src_qc)}</td>'
            f'<td>{s(f"{float(row["pct_counts_mt"]):.2f}%", f"qc-mt-{key}", src_qc)}</td></tr>'
        )

    marker_present = ", ".join(metrics["marker_panel"]["present"])
    marker_missing = ", ".join(metrics["marker_panel"]["missing"])
    image_specs = [
        ("QC violin", "SC000895-R4_008um_gpu_qc_violin.png", "65.8万个散点导致严重遮挡，不适合阈值判断。"),
        ("UMAP + Leiden", "SC000895-R4_008um_gpu_umap_leiden.png", "1,452项图例将主图压缩到几乎不可读。"),
        ("空间 Leiden", "SC000895-R4_008um_gpu_spatial_leiden.png", "图例不可读，且未叠加病理图像。"),
        ("Marker UMAP", "SC000895-R4_008um_gpu_umap_marker_genes.png", "34个预设 marker 中仅保留3个，不能支持前列腺癌注释。"),
    ]
    image_cards = []
    for index, (title, filename, caption) in enumerate(image_specs, start=1):
        uri = image_uri(args.input_dir / "results" / filename)
        image_cards.append(f'''<figure class="image-panel source-figure">
          <div class="card-head"><h3>{html.escape(title)}</h3></div>
          <div class="image-stage"><img src="{uri}" alt="{html.escape(title)}" loading="lazy"></div>
          <figcaption>{html.escape(caption)}</figcaption>
          <button type="button" class="source-tooltip" aria-describedby="image-source-{index}">Source<span class="source-tooltip-content" id="image-source-{index}" role="tooltip">Source: SC000895-R4 GPU pilot output<br>File: {html.escape(filename)}</span></button>
        </figure>''')

    main_html = f'''<main data-report-audience="technical">
      <article class="reading">
        <div class="kicker">10x Visium HD · Human prostate cancer · 008um pilot</div>
        <header data-contract-section="title"><h1>SC000895-R4 GPU Pilot 分析审阅报告</h1></header>
        <section class="technical-summary" data-contract-section="technical-summary">
          <h2>技术跑通，但当前结果不能扩展到全部样本</h2>
          <p><strong>GPU 预处理、CAGRA 邻居图、UMAP、空间坐标和 H5AD 写出均已完成。</strong>但分析对象丢失 raw counts 和大多数 marker，CPU Leiden 生成了 {s('1,452', 'summary-clusters')} 个 clusters，图形也无法人工审阅。</p>
          <ul>
            <li>建议决策：<span class="status stop">暂停 10 样本 array</span>，先完成第二版 R4 pilot。</li>
            <li>主要质量问题：中位 008um bin 仅 {s('12 counts / 11 genes', 'summary-qc', src_qc)}，且 {s('31/34', 'summary-markers')} 个预设 marker 不在最终对象中。</li>
            <li>性能瓶颈：GPU 步骤已明显加速，但 GPU Leiden 环境不兼容导致 CPU 回退，单独耗时 {s('84.14 min', 'summary-leiden', 'results/pilot_gpu_report.json')}。</li>
          </ul>
        </section>
      </article>

      <article class="reading" data-contract-section="key-findings">
        <section class="metrics">
          <div class="metric"><div class="metric-label">QC 后 bins</div><div class="metric-value">{s('658,342', 'metric-bins')}</div><div class="metric-note">008um nonzero bins</div></div>
          <div class="metric"><div class="metric-label">Leiden clusters</div><div class="metric-value danger">{s('1,452', 'metric-clusters')}</div><div class="metric-note">805 clusters &lt;100 bins</div></div>
          <div class="metric"><div class="metric-label">Marker 保留</div><div class="metric-value danger">{s('3 / 34', 'metric-markers')}</div><div class="metric-note">HVG 裁剪后</div></div>
          <div class="metric"><div class="metric-label">Median complexity</div><div class="metric-value">{s('12 / 11', 'metric-complexity', src_qc)}</div><div class="metric-note">counts / genes</div></div>
        </section>
        <section class="narrative"><h2>1,452 个 clusters 表明图结构过度碎片化</h2><p>最大两类容纳 {s('35.2%', 'cluster-top2')} bins，但剩余观测被分散到大量小类；cluster 中位规模仅 {s('86 bins', 'cluster-median')}，最小为 {s('5 bins', 'cluster-min')}。这不支持稳定的 marker 检验或细胞类型注释。</p></section>
      </article>
      {top_chart}

      <article class="reading">
        <section class="narrative"><h2>008um bins 过于稀疏，当前 QC 只删除零计数远远不够</h2><p>{s('P50', 'qc-p50-label', src_qc)} 仅为 {s('12 counts', 'qc-p50-counts', src_qc)} 和 {s('11 genes', 'qc-p50-genes', src_qc)}；同时有 {s('10,881', 'qc-mt20')} 个 bins 的线粒体比例≥20%。应先识别组织内 bins，再根据分布定义 `passes_QC`。</p></section>
      </article>
      {qc_chart}
      <article class="reading"><section class="card table-card"><div class="card-head"><h3>QC 分位数</h3><p>658,342 个非零 008um bins</p></div><div class="table-scroll"><table><thead><tr><th>分位数</th><th>Total counts</th><th>Detected genes</th><th>MT %</th></tr></thead><tbody>{''.join(qc_table_rows)}</tbody></table></div></section></article>

      <article class="reading"><section class="narrative"><h2>对象结构阻断了 marker-based 前列腺癌解读</h2><p>最终 `.h5ad` 没有 `.raw` 或 counts layer，并从 QC 后 {s('17,870 genes', 'full-genes')} 裁剪为 {s('3,000 HVGs', 'hvg-genes')}。预设 marker 仅保留 <strong>{html.escape(marker_present)}</strong>；缺失的 31 个 marker 包括 <span class="muted-list">{html.escape(marker_missing)}</span>。</p></section></article>

      <article class="reading"><section class="narrative"><h2>GPU 已加速邻居图和 UMAP，CPU Leiden 是当前主要瓶颈</h2><p>CAGRA 邻居图仅需 {s('30.58 s', 'runtime-neighbors', 'results/pilot_gpu_report.json')}，GPU UMAP 需 {s('10.28 min', 'runtime-umap', 'results/pilot_gpu_report.json')}。由于 `dask.dataframe.dask_expr` 缺失，GPU Leiden 回退至 CPU，占了绝大多数时间。</p></section></article>
      {runtime_chart}

      <article class="reading"><section class="narrative"><h2>原始图像证实当前输出不适合人工审阅</h2><p>QC 散点过密，UMAP 和空间图图例过长，marker 面板严重不完整。下一版需使用 histogram/hexbin、主要 cluster 或合并标签，并叠加病理图像。</p></section></article>
      <div class="figure-grid">{''.join(image_cards)}</div>

      <article class="reading" data-contract-section="scope-data-and-metric-definitions">
        <section class="narrative"><h2>范围、数据和指标定义</h2><p>本报告仅审阅人前列腺癌样本 SC000895-R4 的 square_008um GPU pilot。“bin complexity”指每个保留 bin 的 total counts 和 detected genes；“cluster size”指 Leiden label 中的 bin 数。本报告不将当前 cluster 解读为细胞类型或病理区域。</p></section>
      </article>

      <article class="reading" data-contract-section="methodology">
        <section class="narrative"><h2>方法和完整性检查</h2><p>流程使用 rapids-singlecell 0.15.2 完成 GPU 归一化、HVG、PCA、CAGRA neighbors 和 UMAP；Leiden 回退至 Scanpy/leidenalg。空间坐标从 Visium HD barcode 解析，全部为有限值。本地 H5AD 与远程 SHA-256 一致：<code>0d083e84ef5b80cf1e200e22d7c6306c2c6994ade468e77009d2f004015bf855</code>。</p></section>
      </article>

      <article class="reading" data-contract-section="limitations-uncertainty-and-robustness-checks">
        <section class="narrative"><h2>局限性会直接改变当前结论</h2><ul><li>未保留 raw counts 和全基因，无法做可追溯的 marker 验证。</li><li>未输出邻接图 connected components，暂不能区分是低深度、CAGRA 参数还是 Leiden 实现导致碎片化。</li><li>空间坐标虽可用，但未与病理图像做可视化配准。</li><li>全项目 Space Ranger QC 字段仍为 NaN，样本间质量尚无法比较。</li></ul></section>
      </article>

      <article class="reading" data-contract-section="recommended-next-steps">
        <section class="narrative"><h2>第二版 pilot 应先修正对象结构和聚类诊断</h2><ol><li>保留全部 QC 后 genes 和 raw counts layer，HVG 仅作 PCA mask。</li><li>根据组织内 bins 的 counts/genes/MT 分布建立 `passes_QC`。</li><li>在 Leiden 前计算 connected components，并测试更低 resolution 和最小 cluster 规模。</li><li>将 RAPIDS GPU 核心环境与 Squidpy/SpatialData 绘图环境分开，修复 GPU Leiden。</li><li>用 016um 建立全样本概览，008um 仅用于局部细化；8 个有 segmented 输出的样本并行做 cell-level 对照。</li><li>重做无超长图例的 QC、UMAP、空间及 marker 图。</li></ol></section>
      </article>

      <article class="reading" data-contract-section="further-questions">
        <section class="narrative"><h2>下一轮需要回答的问题</h2><ul><li>R4 在 016um 下能否形成可解读的主要组织区域？</li><li>当前 1,452 个 clusters 中有多少来自图的独立 connected components？</li><li>组织内 bin 的 QC 谷值和合理下限在哪里？</li><li>segmented 输出与 016um/008um 在上皮、免疫、基质和内皮 marker 上是否一致？</li></ul></section>
      </article>
    </main>'''

    extra_css = '''
    .shell, .card, .metric, .image-panel { border-radius: 8px; }
    .mark { border-radius: 4px; background: var(--blue); }
    h1, h2, .metric-value { letter-spacing: 0; }
    .technical-summary { margin: 28px 0 34px; padding: 22px 0; border-top: 1px solid var(--border-strong); border-bottom: 1px solid var(--border-strong); }
    .technical-summary h2 { font-size: 22px; }
    .technical-summary p, .technical-summary li, .narrative li { color: var(--secondary); }
    .status { display: inline-block; padding: 2px 8px; border-radius: 4px; font-weight: 650; }
    .status.stop { color: var(--warning); background: var(--warning-bg); }
    .metric-value.danger { color: var(--warning); }
    .figure-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; width: min(1080px, calc(100% - 32px)); margin: 22px auto 52px; }
    .image-panel { position: relative; overflow: hidden; border: 1px solid var(--border); background: var(--surface); }
    .image-stage { display: flex; height: 360px; align-items: center; justify-content: center; padding: 12px; overflow: hidden; background: #fff; }
    .image-stage img { display: block; max-width: 100%; max-height: 336px; object-fit: contain; }
    .image-panel figcaption { min-height: 64px; padding: 12px 16px; color: var(--secondary); font-size: 12px; line-height: 19px; }
    .muted-list { color: var(--secondary); }
    code { overflow-wrap: anywhere; }
    @media (max-width: 760px) { .figure-grid { grid-template-columns: 1fr; width: calc(100% - 24px); } .image-stage { height: 300px; } }
    @media print { .shell { box-shadow: none; } .source-tooltip-content { display: none !important; } }
    '''

    shell = args.shell.read_text()
    shell = shell.replace("{{TITLE}}", "SC000895-R4 GPU Pilot 分析审阅报告")
    shell = shell.replace("{{SOURCE_AND_DATE}}", "SC000895-R4 · 10 Jul 2026")
    shell = shell.replace("Data Analytics", "Visium HD Technical Review")
    shell = shell.replace("</style>", extra_css + "</style>")
    shell = re.sub(r'<main data-report-audience=.*?</main>', main_html, shell, count=1, flags=re.DOTALL)
    args.output_shell.write_text(shell)

    payload = {
        "charts": [
            {
                "id": "top-clusters", "height": 420, "type": "bar",
                "dataset": {"id": "top-clusters", "title": "Top 10 Leiden clusters", "data": top_clusters,
                    "chart_spec": {"id": "top-clusters", "dataset": "top-clusters", "title": "Top 10 Leiden clusters", "type": "bar",
                        "encodings": {"x": {"field": "cluster", "type": "nominal"}, "y": {"field": "n_bins", "type": "quantitative", "label": "Bins"}},
                        "settings": {"orientation": "horizontal", "groupMode": "grouped"}, "xAxisTitle": "", "yAxisTitle": "Bins", "valueFormat": "number"}},
            },
            {
                "id": "qc-percentiles", "height": 360, "type": "line",
                "dataset": {"id": "qc-percentiles", "title": "008um bin complexity by percentile", "data": qc_chart_rows,
                    "chart_spec": {"id": "qc-percentiles", "dataset": "qc-percentiles", "title": "008um bin complexity by percentile", "type": "line",
                        "encodings": {"x": {"field": "percentile", "type": "quantitative", "label": "Percentile"}, "y": {"field": "value", "type": "quantitative", "label": "Value"}, "color": {"field": "metric", "type": "nominal"}},
                        "xAxisTitle": "Percentile", "yAxisTitle": "Value", "valueFormat": "number"}},
            },
            {
                "id": "stage-runtime", "height": 320, "type": "bar",
                "dataset": {"id": "stage-runtime", "title": "Runtime by major stage", "data": stage_rows,
                    "chart_spec": {"id": "stage-runtime", "dataset": "stage-runtime", "title": "Runtime by major stage", "type": "bar",
                        "encodings": {"x": {"field": "stage", "type": "nominal"}, "y": {"field": "minutes", "type": "quantitative", "label": "Minutes"}},
                        "settings": {"orientation": "horizontal", "groupMode": "grouped"}, "xAxisTitle": "", "yAxisTitle": "Minutes", "valueFormat": "decimal"}},
            },
        ]
    }
    args.payload.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
