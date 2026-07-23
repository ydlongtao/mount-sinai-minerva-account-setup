#!/usr/bin/env python3
from __future__ import annotations

import html
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "local_results" / "batch_v2_no_legacy"
OUT = ROOT / "docs" / "step8_review_v2_no_legacy.html"


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def esc(value: object) -> str:
    return html.escape(str(value))


def main() -> None:
    scope = read_json(DATA / "config" / "scope_manifest.json")
    int16 = read_json(DATA / "integration" / "integrated_016um_summary.json")
    intcell = read_json(DATA / "integration" / "integrated_cell_level_summary.json")
    summaries = []
    for path in sorted((DATA / "sample_summaries").glob("*/spatial_domain_summary.json")):
        summaries.append(read_json(path))

    rows = []
    for d in summaries:
        domains = d.get("domains", {})
        rows.append(
            "<tr>"
            f"<td>{esc(d.get('sample_id'))}</td>"
            f"<td>{esc(d.get('status'))}</td>"
            f"<td>{esc(d.get('shape', ['?', '?'])[0])}</td>"
            f"<td>{esc(d.get('shape', ['?', '?'])[1])}</td>"
            f"<td>{esc(len(domains))}</td>"
            f"<td>{esc(d.get('spatial_weight'))}</td>"
            "</tr>"
        )

    cards = []
    for d in summaries:
        sample = d["sample_id"]
        image = DATA / "figures" / "overlay" / f"{sample}_spatial_domain_he_overlay.png"
        image_html = f'<img src="../local_results/batch_v2_no_legacy/figures/overlay/{sample}_spatial_domain_he_overlay.png" alt="{sample} H&E spatial domain overlay">' if image.exists() else "<p>图像未下载</p>"
        domain_counts = ", ".join(f"{k}: {v}" for k, v in d.get("domains", {}).items())
        cards.append(f"""
        <article class="sample-card">
          <h3>{esc(sample)}</h3>
          {image_html}
          <p><b>对象数：</b>{esc(d.get('shape', ['?', '?'])[0])}；<b>基因数：</b>{esc(d.get('shape', ['?', '?'])[1])}；<b>空间域数：</b>{esc(len(d.get('domains', {})))}</p>
          <p><b>域组成：</b>{esc(domain_counts)}</p>
          <p class="muted">输入：{esc(Path(d.get('input_h5ad', '')).name)}</p>
        </article>
        """)

    html_text = f"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<title>Step 8 前审阅报告：Visium HD 前列腺癌 v2</title>
<style>
body {{ font-family: -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; margin: 0; color: #1f2933; background: #f5f7fa; line-height: 1.55; }}
main {{ max-width: 1180px; margin: auto; padding: 32px 24px 60px; }}
h1 {{ margin-bottom: 4px; }} h2 {{ margin-top: 34px; border-bottom: 2px solid #d9e2ec; padding-bottom: 6px; }}
.subtitle,.muted {{ color: #52606d; }}
.status {{ display: inline-block; padding: 3px 10px; border-radius: 12px; background: #d9f2e6; color: #176b45; font-weight: 600; }}
.warning {{ background: #fff4cc; border-left: 5px solid #d99e00; padding: 12px 16px; }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fit,minmax(270px,1fr)); gap: 16px; }}
.sample-card {{ background: white; border: 1px solid #d9e2ec; border-radius: 6px; padding: 14px; }}
.sample-card img {{ width: 100%; height: 230px; object-fit: contain; background: #f8fafc; }}
table {{ width: 100%; border-collapse: collapse; background: white; }} th,td {{ border: 1px solid #d9e2ec; padding: 8px; text-align: left; }} th {{ background: #eaf0f6; }}
li {{ margin: 6px 0; }} label {{ display: block; margin: 8px 0; }}
.metric {{ background: white; border: 1px solid #d9e2ec; padding: 14px; border-radius: 6px; }}
code {{ background: #eef2f7; padding: 2px 4px; border-radius: 3px; }}
</style>
</head>
<body><main>
<h1>Step 8 前审阅报告</h1>
<p class="subtitle">Visium HD 前列腺癌正式分析 v2；排除 TD006859-B408 与 TD006859-B573</p>
<p><span class="status">Step 7 已完成</span> 本报告用于确认是否进入空间邻域、共定位和空间通讯分析，不替代正式统计报告。</p>

<div class="warning"><b>解释边界：</b>当前空间域模型是基于表达 PCA 与空间坐标的第一版 Leiden 模型；marker 结果是初步评分，不是人工确认的细胞类型标签。空间通讯结果必须以样本为重复单位进行汇总。</div>

<h2>1. 队列与完成状态</h2>
<div class="grid">
<div class="metric"><b>纳入样本</b><br>{esc(len(scope['included_sample_ids']))} 个 SC000895 Visium HD 样本</div>
<div class="metric"><b>排除样本</b><br>{esc(', '.join(scope['excluded_sample_ids']))}</div>
<div class="metric"><b>016um 整合</b><br>{esc(int16.get('shape'))}；状态 {esc(int16.get('status'))}</div>
<div class="metric"><b>细胞级整合</b><br>{esc(intcell.get('shape'))}；状态 {esc(intcell.get('status'))}</div>
</div>
<p>纳入样本：{esc(', '.join(scope['included_sample_ids']))}</p>

<h2>2. Step 7 样本级空间域结果</h2>
<table><thead><tr><th>样本</th><th>状态</th><th>对象数</th><th>基因数</th><th>空间域数</th><th>空间权重</th></tr></thead><tbody>{''.join(rows)}</tbody></table>
<div class="grid" style="margin-top:16px">{''.join(cards)}</div>

<h2>3. 进入 Step 8 前的确认清单</h2>
<label><input type="checkbox"> 8 个样本的空间域图均能与组织覆盖和表达分布对应，没有明显空域、孤立域或坐标错位。</label>
<label><input type="checkbox"> 细胞类型标签将采用 marker + QC + 空间位置综合确认，不直接把 Leiden domain 当作细胞类型。</label>
<label><input type="checkbox"> 通讯分析的基本统计单位确定为“样本内的细胞类型/空间域对”，而不是单个 cell/bin。</label>
<label><input type="checkbox"> 低细胞数或低 counts 样本只做敏感性分析，不作为独立强证据。</label>
<label><input type="checkbox"> 已确认配体-受体数据库版本和物种为 human，并记录数据库版本。</label>
<label><input type="checkbox"> 已确认空间邻域距离/邻居数阈值，并计划做至少一组敏感性分析。</label>

<h2>4. 建议的 Step 8 执行顺序</h2>
<ol>
<li>先在每个样本内建立细胞级空间邻接图，计算细胞类型/空间域的邻域富集。</li>
<li>再计算肿瘤-CAF、肿瘤-免疫、肿瘤-血管等候选关系的样本级统计。</li>
<li>使用 human 配体-受体数据库开展通讯候选筛选，并报告每个样本的效应方向和一致性。</li>
<li>最后进行跨样本汇总和敏感性分析；不把所有细胞拼接后直接当作独立重复。</li>
</ol>

<h2>5. 当前尚未执行的内容</h2>
<ul><li>空间邻域富集与共定位</li><li>配体-受体空间通讯</li><li>样本级统计检验与敏感性分析</li><li>Spatial transition tensor</li></ul>
<p class="muted">原始结果位置：<code>/sc/arion/work/huangl21/sp_project/results/batch_v2_no_legacy/</code><br>本地摘要目录：<code>local_results/batch_v2_no_legacy/</code></p>
</main></body></html>
"""
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(html_text)
    print(OUT)


if __name__ == "__main__":
    main()
