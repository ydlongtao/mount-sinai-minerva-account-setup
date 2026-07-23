#!/usr/bin/env python3
"""Build a local HTML review report for R1/R7 H&E candidate overlays."""
from __future__ import annotations

import html
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "local_results" / "segmented_official_v1" / "omicverse_pyinfercnv_i3_epithelial_v2_r1_r7_he_overlay"
OUT = DATA / "r1_r7_i3_epithelial_he_overlay_report.html"


def e(x: object) -> str:
    return html.escape(str(x))


def main() -> None:
    summary = json.loads((DATA / "overlay_summary.json").read_text())
    sections = []
    for row in summary["samples"]:
        sample = row["sample"]
        imgs = {
            "高置信恶性候选": f"{sample}_he_high_confidence_malignant_overlay.png",
            "扩大上皮候选": f"{sample}_he_expanded_epithelial_overlay.png",
            "联合图（红色恶性，青色扩大候选）": f"{sample}_he_i3_epithelial_combined_overlay.png",
        }
        figures = "".join(
            f'<figure><img src="{e(sample + "/" + filename)}" alt="{e(sample + " " + title)}"><figcaption>{e(title)}</figcaption></figure>'
            for title, filename in imgs.items()
        )
        sections.append(f"""
        <section class="sample">
          <h2>{e(sample)}</h2>
          <p><b>高置信恶性候选：</b>{row['malignant_cells']:,}；<b>扩大上皮候选：</b>{row['expanded_epithelial_cells']:,}；<b>仅扩大候选：</b>{row['expanded_only_cells']:,}；图像尺寸：{row['image_width']} × {row['image_height']}。</p>
          <p>高置信恶性 cluster：<code>{e(', '.join(row['malignant_clusters']) or '无')}</code></p>
          <div class="grid">{figures}</div>
        </section>""")
    report = f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>R1/R7 i3 H&amp;E overlay</title>
<style>
body{{font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;color:#20252b;max-width:1500px;margin:28px auto;padding:0 18px}}
h1{{margin-bottom:4px}} .muted{{color:#68737d}} .note{{background:#f3f6f8;border-left:4px solid #3b6f8f;padding:12px 16px;margin:18px 0}}
.sample{{border-top:1px solid #ccd5dc;margin-top:28px;padding-top:22px}} .grid{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:16px}}
figure{{margin:0;border:1px solid #d7dee4;padding:8px;background:#fff}} img{{display:block;width:100%;height:auto}} figcaption{{font-weight:600;padding:7px 2px 2px}}
code{{background:#f1f4f6;padding:2px 5px}} li{{margin:6px 0}}
@media(max-width:1000px){{.grid{{grid-template-columns:1fr 1fr}}}} @media(max-width:650px){{.grid{{grid-template-columns:1fr}};body{{margin:14px}}}}
</style></head><body>
<h1>R1/R7：i3 候选标签 H&amp;E overlay</h1>
<p class="muted">基于 OmicVerse py-inferCNV i3 v2 输出和官方 segmented H&amp;E hires 图。</p>
<div class="note"><b>判读规则：</b>红色是 marker 预注释为 <code>malignant_epithelial_candidate</code> 的高置信候选；青色是 i3 v2 的 <code>epithelial_candidate</code> 集合，包括 malignant、luminal、basal 和符合阈值的 epithelial-like ambiguous clusters。青色不等于恶性。R7 没有红色点表示当前 marker 规则没有高置信 malignant cluster，不表示 R7 没有肿瘤细胞。</div>
<h2>建议核查</h2><ol><li>红色点是否集中在腺体上皮、腺腔周围或形态异常的上皮区域。</li><li>青色区域是否覆盖主要上皮组织，同时避免大量落在间质、炎症或空白区域。</li><li>R1 与 R7 的候选区域是否存在空间结构差异，并与 CNV heatmap 中的连续异常相互支持。</li><li>注意点是细胞中心坐标而不是边界，因此密集区域会显示为点云。</li></ol>
{''.join(sections)}
</body></html>"""
    OUT.write_text(report, encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
