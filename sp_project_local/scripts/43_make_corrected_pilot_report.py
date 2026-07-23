#!/usr/bin/env python3
from __future__ import annotations

import html
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "local_results" / "batch_v2_corrected"
OUT = ROOT / "docs" / "corrected_cellpose_pilot_report.html"
SAMPLES = ["SC000895-R4", "SC000895-R5", "SC000895-R9"]


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def e(x: object) -> str:
    return html.escape(str(x))


def main() -> None:
    cards = []
    rows = []
    for sample in SAMPLES:
        summary_dir = DATA / "summaries" / sample
        he_path = summary_dir / "cellpose_he_summary.json"
        raw_path = summary_dir / "raw_counts_cellseg_summary.json"
        qc_path = summary_dir / "cell_qc_marker_summary.json"
        overlay = DATA / "figures" / f"{sample}_corrected_he_coverage.png"
        if not he_path.exists():
            status = "CPU 重试失败：缺少 cv2，尚未生成校正 Cellpose H5AD"
            rows.append(f"<tr><td>{e(sample)}</td><td class=\"fail\">{e(status)}</td><td colspan=5>待安装 OpenCV 后重跑</td></tr>")
            cards.append(f"<article class=\"card\"><h3>{e(sample)}</h3><p class=\"fail\">{e(status)}</p><p>本样本没有被纳入校正通过结论。</p></article>")
            continue
        he = load(he_path); raw = load(raw_path); qc = load(qc_path)
        expanded = he["labels_he_expanded_nonzero_bins"]
        bins = he["n_bins"]
        cells = qc["cells"]
        coverage = 100 * expanded / bins
        rows.append(f"<tr><td>{e(sample)}</td><td class=\"pass\">pass</td><td>{e(he.get('source_mpp'))}</td><td>{e(he.get('cellpose_mpp_used'))}</td><td>{e(he.get('labels_he_nonzero_bins'))}</td><td>{e(expanded)}</td><td>{coverage:.1f}%</td><td>{e(cells)}</td></tr>")
        img = f'<img src="../local_results/batch_v2_corrected/figures/{sample}_corrected_he_coverage.png" alt="{sample} corrected H&E coverage">' if overlay.exists() else "<p>叠加图未下载</p>"
        cards.append(f"<article class=\"card\"><h3>{e(sample)}</h3>{img}<p><b>Cellpose source mpp：</b>{e(he.get('source_mpp'))}；<b>扩展标签覆盖：</b>{coverage:.1f}%；<b>bin2cell cells：</b>{e(cells)}</p><p><b>median counts：</b>{e(qc.get('median_total_counts'))}；<b>median genes：</b>{e(qc.get('median_n_genes'))}</p></article>")

    report = f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><title>校正 Cellpose 试点报告</title>
<style>body{{font:14px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;color:#1f2933;background:#f5f7fa;margin:0}}main{{max-width:1200px;margin:auto;padding:30px 22px 60px}}h1{{margin-bottom:4px}}h2{{margin-top:32px;border-bottom:2px solid #d9e2ec;padding-bottom:6px}}table{{border-collapse:collapse;width:100%;background:#fff}}th,td{{border:1px solid #d9e2ec;padding:8px;text-align:left}}th{{background:#eaf0f6}}.pass{{color:#176b45;font-weight:600}}.fail{{color:#b42318;font-weight:600}}.warning{{background:#fff4cc;border-left:5px solid #d99e00;padding:12px 16px}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(430px,1fr));gap:16px}}.card{{background:#fff;border:1px solid #d9e2ec;border-radius:6px;padding:14px}}.card img{{width:100%;height:auto}}code{{background:#eef2f7;padding:2px 4px;border-radius:3px}}</style></head><body><main>
<h1>校正 Cellpose / bin2cell 试点报告</h1><p>Visium HD 前列腺癌；试点样本：R4、R5、R9。结果目录：<code>batch_v2_corrected</code>。</p>
<div class="warning"><b>本报告的判断：</b>R4/R5 已使用实际 source mpp 完成 Cellpose、bin2cell 和 cell-level QC；R9 的 CPU 重试因 omicverse CPU 环境缺少 <code>cv2</code> 未完成。当前不启动全量 Cellpose，也不进入 Step 8。</div>
<h2>1. 校正结果对照</h2>
<table><thead><tr><th>样本</th><th>状态</th><th>source mpp</th><th>使用 mpp</th><th>HE 标签 bins</th><th>扩展标签 bins</th><th>扩展覆盖率</th><th>bin2cell cells</th></tr></thead><tbody>{''.join(rows)}</tbody></table>
<h2>2. H&E 覆盖与细胞中心</h2><div class="grid">{''.join(cards)}</div>
<h2>3. 当前建议</h2><ol><li>先安装/确认 CPU 环境的 OpenCV，完成 R9 source-mpp CPU 重试。</li><li>审阅 R4/R5 叠加图中橙色扩展 bins 是否覆盖 H&E 组织，青色 bin2cell 中心是否落在可识别细胞区域。</li><li>若 R4/R5/R9 均通过，再用 source mpp 批量重跑其余 5 个样本。</li><li>8 µm bins 的空间域分析仍作为主空间域路线；Cellpose 只作为细胞级辅助路线。</li></ol>
<p class="muted">旧 v1/v2 结果未删除；本试点独立输出到 <code>/sc/arion/work/huangl21/sp_project/results/batch_v2_corrected/</code>。</p>
</main></body></html>"""
    OUT.write_text(report)
    print(OUT)


if __name__ == "__main__": main()
