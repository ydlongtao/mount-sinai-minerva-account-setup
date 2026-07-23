#!/usr/bin/env python3
"""Build the local report for the final user-approved epithelial grouping."""
from __future__ import annotations
import html, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "local_results" / "segmented_official_v1" / "final_epithelial_groups_r1_r7"
OUT = DATA / "r1_r7_final_epithelial_groups_report.html"

def e(x): return html.escape(str(x))

def main():
    summary = json.loads((DATA / "final_epithelial_groups_summary.json").read_text())
    sections = []
    for row in summary["samples"]:
        sample = row["sample"]; counts = row["final_group_counts"]
        rows = "".join(f"<tr><td>{e(k)}</td><td>{v:,}</td><td>{100*v/row['cells']:.2f}%</td></tr>" for k,v in counts.items())
        image = f"{sample}/{sample}_final_epithelial_groups_he_overlay.png"
        sections.append(f"<section class='sample'><h2>{e(sample)}</h2><table><tr><th>最终分群</th><th>细胞数</th><th>比例</th></tr>{rows}</table><figure><img src='{e(image)}' alt='{e(sample)} final epithelial grouping H&E overlay'><figcaption>红色：malignant_epithelial；绿色：normal_epithelial；蓝色/紫色：保留的 luminal/neuroendocrine 类别。</figcaption></figure></section>")
    report = f"""<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>R1/R7 final epithelial groups</title><style>body{{font:15px/1.55 -apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;color:#20252b;max-width:1500px;margin:28px auto;padding:0 18px}}h1{{margin-bottom:4px}}.muted{{color:#68737d}}.note{{background:#f3f6f8;border-left:4px solid #3b6f8f;padding:12px 16px;margin:18px 0}}.sample{{border-top:1px solid #ccd5dc;margin-top:28px;padding-top:22px}}table{{border-collapse:collapse;width:100%;margin:12px 0 20px}}th,td{{border:1px solid #dfe4e8;padding:7px 9px;text-align:left}}th{{background:#f7f8f9}}img{{width:100%;height:auto;border:1px solid #dfe4e8}}figcaption{{color:#5f6b76;font-size:13px;padding-top:5px}}code{{background:#f1f4f6;padding:2px 5px}}</style></head><body><h1>R1/R7 最终上皮细胞分群</h1><p class='muted'>基于核心 marker 与 OmicVerse i3 CNV 候选的用户确认版。</p><div class='note'><b>最终合并规则：</b><code>malignant_marker_candidate</code> + <code>malignant_marker_CNV_candidate</code> → <b>malignant_epithelial</b>；<code>epithelial_CNV_candidate</code> + <code>epithelial_ambiguous</code> + <code>basal_epithelial</code> → <b>normal_epithelial</b>。<code>luminal_epithelial</code>、<code>neuroendocrine_like</code> 和 <code>non_epithelial_or_low_signal</code> 保留原标签，没有被静默合并。</div>{''.join(sections)}<h2>使用建议</h2><ol><li>下游差异分析可将 malignant_epithelial 与 normal_epithelial 作为主要对照分组。</li><li>luminal_epithelial 和 neuroendocrine_like 建议先作为单独的保留类别，后续结合病理和 CNV 决定是否纳入 normal 或 malignant。</li><li>最终标签文件位于每个样本目录的 <code>*_final_epithelial_groups.csv.gz</code>。</li></ol></body></html>"""
    OUT.write_text(report, encoding='utf-8'); print(OUT)

if __name__ == '__main__': main()
