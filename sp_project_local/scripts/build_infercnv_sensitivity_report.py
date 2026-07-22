#!/usr/bin/env python3
"""Build a local comparison report for infercnvpy reference sensitivity runs."""
from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

MODES = ["main_reference", "endothelial", "fibroblast_stromal", "smooth_muscle_pericyte", "auxiliary_reference"]
SAMPLES = ["SC000895-R1", "SC000895-R7"]
LABELS = {
    "main_reference": "综合主参考",
    "endothelial": "仅内皮参考",
    "fibroblast_stromal": "仅成纤维/基质参考",
    "smooth_muscle_pericyte": "仅平滑肌/周细胞参考",
    "auxiliary_reference": "免疫辅助参考",
}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    base = args.input_dir.resolve()
    cards, table_rows = [], []
    for mode in MODES:
        for sample in SAMPLES:
            folder = base / mode / sample
            summary = json.loads((folder / "pilot_summary.json").read_text())
            image = f"{mode}/{sample}/{sample}_cnv_reference_heatmap.png"
            cards.append(f"""
<section class="card"><h3>{html.escape(LABELS[mode])} · {html.escape(sample)}</h3>
<img src="{image}" alt="{html.escape(mode)} {html.escape(sample)} CNV heatmap">
<p>参考细胞：<b>{summary['reference_cells']:,}</b>；query：{summary['query_cells']:,}；基因：{summary['genes_with_positions']:,}</p>
</section>""")
            table_rows.append(f"<tr><td>{html.escape(LABELS[mode])}</td><td>{html.escape(sample)}</td><td>{summary['reference_cells']:,}</td><td>{summary['query_cells']:,}</td><td>{'pass' if summary['status'] == 'pass' else summary['status']}</td></tr>")

    page = f"""<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">
<title>infercnvpy R1/R7 Reference Sensitivity</title>
<style>
body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;background:#f4f6f8;color:#17202a;line-height:1.5;margin:0}}
main{{max-width:1500px;margin:auto;padding:28px}}h1{{margin:0 0 6px}}h2{{color:#1f4e79;margin-top:28px}}
.subtitle{{color:#5d6873;margin-bottom:22px}}.panel,.card{{background:#fff;border:1px solid #d8dee5;border-radius:8px;padding:16px;box-shadow:0 2px 8px #17202a12}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(560px,1fr));gap:18px}}.card img{{display:block;width:100%;border:1px solid #d8dee5;margin:8px 0}}
table{{border-collapse:collapse;width:100%;background:#fff}}th,td{{border:1px solid #d8dee5;padding:8px;text-align:left}}th{{background:#eef2f5}}
code{{background:#eef1f4;padding:2px 5px;border-radius:4px}}.pass{{color:#177245;font-weight:600}}
@media(max-width:700px){{main{{padding:14px}}.grid{{grid-template-columns:1fr}}}}
</style></head><body><main>
<h1>infercnvpy R1/R7 参考细胞敏感性分析</h1>
<div class="subtitle">比较不同候选二倍体参考对前列腺癌空间转录组 CNV 推断的影响</div>
<div class="panel"><h2>分析设计</h2>
<p>每个参考模式均使用 R1/R7 的全部 query cells，并分别以单一参考类型或综合主参考为 infercnvpy 的 reference。综合主参考包含 endothelial、fibroblast/stromal、smooth muscle/pericyte；免疫参考仅用于敏感性分析。</p>
<p>判定稳定 CNV 时，应优先寻找在多个参考模式下方向一致、覆盖相邻多个基因和染色体区段的信号。单一参考模式特异的信号应标记为低置信度，不能直接解释为恶性 CNV。</p></div>
<h2>运行汇总</h2><table><thead><tr><th>参考模式</th><th>样本</th><th>参考细胞</th><th>Query 细胞</th><th>状态</th></tr></thead><tbody>{''.join(table_rows)}</tbody></table>
<h2>热图比较</h2><div class="grid">{''.join(cards)}</div>
<div class="panel"><h2>审阅要点</h2><ul>
<li>先确认不同参考模式下参考组本身是否接近中性。</li>
<li>再比较 R1 与 R7 是否出现相同的连续染色体异常，区分共同事件与复发相关事件。</li>
<li>若异常只在免疫参考下出现，优先考虑参考细胞组成、低 RNA 或免疫表达程序造成的偏移。</li>
<li>最终应将稳定 CNV 区段与恶性上皮候选标签、空间位置和 H&amp;E overlay 叠加验证。</li>
</ul></div>
<p class="subtitle">本报告只包含热图和汇总 JSON；大型 H5AD 文件保留在服务器。</p>
</main></body></html>"""
    args.output.write_text(page, encoding="utf-8")
    print(args.output.resolve())


if __name__ == "__main__":
    main()
