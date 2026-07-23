#!/usr/bin/env python3
"""Build a local HTML review report for the R1/R7 infercnvpy pilot."""
from __future__ import annotations

import html
import json
import argparse
from pathlib import Path

SAMPLES = ["SC000895-R1", "SC000895-R7"]


def read_summary(base: Path, sample: str) -> dict:
    return json.loads((base / sample / "pilot_summary.json").read_text())


def row(label: str, value: object) -> str:
    return f"<tr><th>{html.escape(label)}</th><td>{html.escape(str(value))}</td></tr>"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    base = args.input_dir.resolve()
    out = args.output.resolve()
    cards = []
    for sample in SAMPLES:
        s = read_summary(base, sample)
        image = f"{sample}/{sample}_cnv_reference_heatmap.png"
        cards.append(f"""
<section class="sample">
<h2>{html.escape(sample)}</h2>
<img src="{image}" alt="{html.escape(sample)} CNV reference heatmap">
<table>
{row('状态', s['status'])}
{row('总细胞数', f"{s['cells_total']:,}")}
{row('主参考细胞', f"{s['reference_cells']:,}")}
{row('Query 细胞', f"{s['query_cells']:,}")}
{row('有 Hg38 坐标的基因', f"{s['genes_with_positions']:,}")}
{row('CNV 矩阵', ' x '.join(map(str, s['cnv_matrix_shape'])))}
</table>
</section>""")

    page = f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8">
<title>infercnvpy R1/R7 CNV Pilot Report</title>
<style>
body {{ font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; margin:0; background:#f4f6f8; color:#17202a; line-height:1.55; }}
main {{ max-width:1400px; margin:0 auto; padding:32px; }}
h1 {{ margin:0 0 8px; font-size:30px; }} h2 {{ margin-top:0; color:#1f4e79; }}
.subtitle {{ color:#5d6873; margin-bottom:24px; }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(520px,1fr)); gap:24px; }}
.sample,.panel {{ background:#fff; border:1px solid #d8dee5; border-radius:8px; padding:20px; box-shadow:0 2px 8px #17202a12; }}
.sample img {{ width:100%; height:auto; border:1px solid #d8dee5; background:#fff; margin:6px 0 16px; }}
table {{ width:100%; border-collapse:collapse; font-size:14px; }} th,td {{ text-align:left; border-top:1px solid #e5e9ed; padding:8px; }} th {{ width:42%; color:#53606d; }}
.panel {{ margin-top:24px; }} code {{ background:#eef1f4; padding:2px 5px; border-radius:4px; }}
@media(max-width:700px) {{ main {{ padding:16px; }} .grid {{ grid-template-columns:1fr; }} .sample {{ padding:12px; }} }}
</style></head><body><main>
<h1>infercnvpy R1/R7 CNV Pilot</h1>
<div class="subtitle">人前列腺癌 Visium HD 单细胞级空间转录组初步 CNV 评估</div>
<div class="panel"><h2>本轮设置</h2>
<p>本 pilot 使用 GENCODE v38 的 GRCh38 primary assembly 基因注释。主参考细胞为非恶性候选的 endothelial、fibroblast/stromal 和 smooth muscle/pericyte。免疫细胞未混入主参考，后续作为敏感性分析。</p>
<p>每个样本保留全部选中的主参考细胞，并随机保留最多 40,000 个 query cells。CNV 输入由 counts layer 重新进行总量归一化和 log1p，排除 chrX/chrY。</p>
</div><div class="grid">{''.join(cards)}</div>
<div class="panel"><h2>如何阅读热图</h2>
<p>蓝白红信号表示相对于主参考的平滑表达偏移，不应直接等同于经过 DNA 测序验证的拷贝数。优先关注连续覆盖多个相邻基因的大片段信号，而不是单个基因的孤立变化。</p>
<p>首先检查 <code>main_reference</code> 行是否整体接近中性且内部一致；随后观察 query 中是否存在跨染色体连续的增益/缺失模式，并与恶性上皮候选、空间位置和 H&amp;E overlay 对照。</p>
</div><div class="panel"><h2>当前限制与下一步</h2>
<ul><li>GENCODE v38 对本批 H5AD 的基因坐标覆盖率为 98.55%；232 个基因未匹配，30 个重复标识符被跳过。</li>
<li>参考细胞来自 marker-based 预注释，不是配对正常组织；因此本结果属于探索性 CNV 证据。</li>
<li>下一步应运行全量 R1/R7，并分别使用 endothelial、stromal、smooth-muscle 及 immune auxiliary reference 做敏感性比较。</li>
<li>最终 CNV 阳性候选需同时满足：连续染色体区段、多个参考组合方向一致、空间上富集于恶性上皮候选区域。</li></ul>
</div><p class="subtitle">报告生成自本地下载的 LSF pilot 输出；原始大 H5AD 未复制到本地。</p>
</main></body></html>"""
    out.write_text(page, encoding="utf-8")
    print(out)


if __name__ == "__main__":
    main()
