#!/usr/bin/env python3
from __future__ import annotations
import html, json
from pathlib import Path

ROOT = Path("/sc/arion/work/huangl21/sp_project/results/batch_v2_corrected")
OUT = ROOT / "reports" / "r4_r5_gex_segmentation_report.html"
SAMPLES = ["SC000895-R4", "SC000895-R5"]
TAGS = ["flow02_probm2", "flow04_prob0", "flow06_prob1"]

def e(x): return html.escape(str(x))

def main():
    rows, cards = [], []
    for sample in SAMPLES:
        for tag in TAGS:
            base = ROOT / "samples" / sample / "gex_compare" / tag
            s = json.loads((base / "gex_segmentation_summary.json").read_text())
            rel = Path("../samples") / sample / "gex_compare" / tag / f"{sample}_{tag}_gex_salvage_overlay.png"
            rows.append(f"<tr><td>{e(sample)}</td><td>{e(tag)}</td><td>{s['he_labeled_bins']}</td><td>{s['gex_labeled_bins']}</td><td>{s['salvaged_bins']}</td><td>{s['gex_cells']}</td><td>{s['joint_cells']}</td><td>{s['joint_coverage_percent']:.1f}%</td></tr>")
            cards.append(f"<article><h3>{e(sample)} / {e(tag)}</h3><img src='{e(rel.as_posix())}'><p>H&E 标签 bins：<b>{s['he_labeled_bins']}</b>；GEX 新增 bins：<b>{s['salvaged_bins']}</b>；联合细胞数：<b>{s['joint_cells']}</b>；联合覆盖率：<b>{s['joint_coverage_percent']:.1f}%</b>。</p></article>")
    report = f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><title>R4/R5 GEX 细胞分割报告</title><style>body{{font:14px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;color:#243b53;background:#f5f7fa;margin:0}}main{{max-width:1400px;margin:auto;padding:28px}}table{{border-collapse:collapse;width:100%;background:#fff}}th,td{{border:1px solid #d9e2ec;padding:7px}}th{{background:#eaf0f6}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(460px,1fr));gap:16px}}article{{background:#fff;border:1px solid #d9e2ec;border-radius:6px;padding:12px}}article img{{width:100%;height:auto}}.note{{background:#fff4cc;border-left:5px solid #d99e00;padding:12px}}</style></head><body><main><h1>R4/R5 GEX 细胞分割与 H&E salvage 比较</h1><div class="note">紫色表示 GEX 标签区域，绿色表示只由 GEX 补回的 H&E 未覆盖区域。GEX 结果用于补充密集区域，不应仅根据细胞数量判断优劣。</div><table><tr><th>样本</th><th>参数组</th><th>H&E bins</th><th>GEX bins</th><th>salvaged bins</th><th>GEX cells</th><th>联合 cells</th><th>联合覆盖率</th></tr>{''.join(rows)}</table><h2>Overlay</h2><div class="grid">{''.join(cards)}</div></main></body></html>'''
    OUT.parent.mkdir(parents=True, exist_ok=True); OUT.write_text(report); print(OUT)

if __name__ == '__main__': main()
