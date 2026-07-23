#!/usr/bin/env python3
"""Build a local HTML report for core-panel/CNV epithelial reclassification."""
from __future__ import annotations

import html
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "local_results" / "segmented_official_v1" / "epithelial_cnv_corepanel_reclassification_r1_r7"
OUT = DATA / "r1_r7_epithelial_cnv_corepanel_reclassification_report.html"


def e(x: object) -> str:
    return html.escape(str(x))


def main() -> None:
    summary = json.loads((DATA / "reclassification_summary.json").read_text())
    cards = []
    for row in summary["samples"]:
        sample = row["sample"]
        counts = row["class_counts"]
        label_rows = "".join(f"<tr><td>{e(k)}</td><td>{v:,}</td><td>{100*v/row['cells']:.2f}%</td></tr>" for k, v in counts.items())
        image = f"{sample}/{sample}_corepanel_cnv_reclassification_he_overlay.png"
        cards.append(f"""
        <section class="sample"><h2>{e(sample)}</h2>
          <div class="metrics">
            <div><b>{row['cells']:,}</b><span>cells</span></div>
            <div><b>{row['core_epithelial']:,}</b><span>core epithelial</span></div>
            <div><b>{row['core_malignant_marker']:,}</b><span>core malignant-marker positive</span></div>
            <div><b>{row['cnv_expanded']:,}</b><span>expanded CNV candidate</span></div>
            <div><b>{row['overlap_expanded_core_malignant']:,}</b><span>expanded CNV + malignant marker overlap</span></div>
          </div>
          <table><tr><th>分类</th><th>细胞数</th><th>比例</th></tr>{label_rows}</table>
          <p><b>high-confidence cluster CNV ∩ core malignant marker：</b>{row['overlap_high_cnv_core_malignant']:,}</p>
          <figure><img src="{e(image)}" alt="{e(sample)} core panel CNV reclassification H&E overlay"><figcaption>H&amp;E merge：红色为 marker + CNV candidate；青色为 epithelial CNV candidate；蓝/绿为 luminal/basal。</figcaption></figure>
        </section>""")
    report = f"""<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>R1/R7 epithelial CNV core panel</title>
<style>
body{{font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;color:#20252b;max-width:1500px;margin:28px auto;padding:0 18px}} h1{{margin-bottom:4px}} .muted{{color:#68737d}}
.note{{background:#f3f6f8;border-left:4px solid #3b6f8f;padding:12px 16px;margin:18px 0}} .sample{{border-top:1px solid #ccd5dc;margin-top:30px;padding-top:22px}}
.metrics{{display:grid;grid-template-columns:repeat(5,1fr);gap:10px;margin:15px 0}} .metrics div{{background:#f7f8f9;border:1px solid #e1e6ea;padding:11px}} .metrics b,.metrics span{{display:block}} .metrics b{{font-size:19px}} .metrics span{{font-size:12px;color:#5f6b76}}
table{{border-collapse:collapse;width:100%;margin:12px 0 18px}}th,td{{border:1px solid #dfe4e8;padding:6px 8px;text-align:left}}th{{background:#f7f8f9}} figure{{margin:0}} img{{width:100%;height:auto;border:1px solid #dfe4e8}} figcaption{{color:#5f6b76;font-size:13px;padding-top:5px}} code{{background:#f1f4f6;padding:2px 5px}}
@media(max-width:900px){{.metrics{{grid-template-columns:repeat(2,1fr)}}}} @media(max-width:600px){{.metrics{{grid-template-columns:1fr}};body{{margin:14px}}}}
</style></head><body>
<h1>R1/R7 上皮细胞：核心 marker + i3 CNV 重新分类</h1><p class="muted">细胞级分类报告；R1 为未复发原发灶，R7 为复发原发灶。</p>
<div class="note"><b>判定逻辑：</b>先用 EPCAM/KRT8/KRT18/KRT19 识别核心上皮细胞；再用前列腺腔面、恶性相关、增殖、基底和神经内分泌 panel 打分；最后与 OmicVerse i3 的 <code>epithelial_candidate</code> 和高置信 malignant cluster 合并。高置信 malignant cluster 是 expanded CNV candidate 的子集。阈值：panel score ≥ 0.25；核心上皮 marker detection fraction ≥ 0.25。该结果是候选分层，不是病理诊断。</div>
<h2>类别说明</h2><ul><li><b>malignant_marker_CNV_candidate：</b>核心上皮细胞同时具有恶性/增殖 marker 信号和 expanded CNV candidate 标签。</li><li><b>epithelial_CNV_candidate：</b>核心上皮细胞有 CNV candidate 标签，但未达到恶性/增殖 marker 阈值。</li><li><b>malignant_marker_candidate：</b>具有恶性/增殖 marker 信号，但没有被当前 expanded CNV 标签覆盖。</li><li><b>luminal/basal/neuroendocrine：</b>相应谱系程序阳性的核心上皮细胞。</li></ul>
{''.join(cards)}
<h2>解释重点</h2><ol><li>R1 的 expanded CNV candidate 与 core malignant marker overlap 为 {summary['samples'][0]['overlap_expanded_core_malignant']:,} 个细胞；其中高置信 malignant cluster CNV 与核心 malignant marker overlap 为 {summary['samples'][0]['overlap_high_cnv_core_malignant']:,}。</li><li>R7 的对应 overlap 为 {summary['samples'][1]['overlap_expanded_core_malignant']:,} 个细胞；当前 R7 没有高置信 malignant cluster 标签，因此后一项为 0。</li><li>下一步应优先查看 <code>malignant_marker_CNV_candidate</code> 的连续 CNV 区段、H&amp;E 腺体形态和空间 niche，而不是把所有 expanded epithelial cells 直接称作 malignant。</li></ol>
</body></html>"""
    OUT.write_text(report, encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
