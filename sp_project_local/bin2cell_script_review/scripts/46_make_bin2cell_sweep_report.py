#!/usr/bin/env python3
from __future__ import annotations
import html, json
from pathlib import Path

ROOT = Path("/sc/arion/work/huangl21/sp_project/results/batch_v2_corrected")
OUT = ROOT / "reports" / "r4_r5_cellpose_bin2cell_sweep_report.html"
SAMPLES = ["SC000895-R4", "SC000895-R5"]
TAGS = ["flow02_probm2", "flow04_prob0", "flow06_prob1"]

def esc(x): return html.escape(str(x))

def main():
    rows, cards = [], []
    for sample in SAMPLES:
        for tag in TAGS:
            base = ROOT / "samples" / sample / "cellpose_sweep" / tag
            cp = json.loads((base / "cellpose_he_summary.json").read_text())
            bc = json.loads((base / "bin2cell" / "bin2cell_sweep_summary.json").read_text())
            rel = Path("../samples") / sample / "cellpose_sweep" / tag / "bin2cell" / "figures" / f"{sample}_{tag}_bin2cell_overlay.png"
            rows.append(f"<tr><td>{esc(sample)}</td><td>{esc(tag)}</td><td>{cp.get('cellpose_flow_threshold')}</td><td>{cp.get('cellpose_prob_thresh')}</td><td>{cp.get('source_mpp')}</td><td>{bc['cells']}</td><td>{bc['labeled_bins']}</td><td>{bc['coverage_percent']:.1f}%</td></tr>")
            cards.append(f"<article><h3>{esc(sample)} / {esc(tag)}</h3><img src='{esc(rel.as_posix())}'><p>细胞数：<b>{bc['cells']}</b>；扩展标签覆盖：<b>{bc['coverage_percent']:.1f}%</b>；source mpp：<b>{cp.get('source_mpp')}</b>。</p></article>")
    report = f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><title>R4/R5 Cellpose bin2cell 参数比较</title><style>body{{font:14px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;color:#243b53;background:#f5f7fa;margin:0}}main{{max-width:1400px;margin:auto;padding:28px}}table{{border-collapse:collapse;width:100%;background:#fff}}th,td{{border:1px solid #d9e2ec;padding:7px}}th{{background:#eaf0f6}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(460px,1fr));gap:16px}}article{{background:#fff;border:1px solid #d9e2ec;border-radius:6px;padding:12px}}article img{{width:100%;height:auto}}.note{{background:#fff4cc;border-left:5px solid #d99e00;padding:12px}}</style></head><body><main><h1>R4/R5 Cellpose + bin2cell 参数比较</h1><p>本报告比较 6 组参数；黄色为扩展后有标签的 2um bins，蓝色为 bin2cell 细胞中心。</p><div class="note">细胞数和覆盖率只能作为筛选指标，最终选择必须结合 H&E 形态、细胞边界是否合理以及转录本是否落入细胞。</div><h2>汇总</h2><table><tr><th>样本</th><th>参数组</th><th>flow</th><th>prob</th><th>mpp</th><th>cells</th><th>labeled bins</th><th>coverage</th></tr>{''.join(rows)}</table><h2>H&E / bin2cell 叠加图</h2><div class="grid">{''.join(cards)}</div></main></body></html>'''
    OUT.parent.mkdir(parents=True, exist_ok=True); OUT.write_text(report); print(OUT)

if __name__ == '__main__': main()
