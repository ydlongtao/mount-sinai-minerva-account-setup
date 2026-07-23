#!/usr/bin/env python3
"""Build a local HTML review report for the OmicVerse py-inferCNV i3 v2 run."""

from __future__ import annotations

import html
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "local_results" / "segmented_official_v1" / "omicverse_pyinfercnv_i3_epithelial_v2_r1_r7"
OUT = RESULTS / "r1_r7_omicverse_pyinfercnv_i3_epithelial_v2_report.html"


def esc(value: object) -> str:
    return html.escape(str(value))


def pct(value: int, total: int) -> str:
    return f"{100.0 * value / total:.2f}%" if total else "NA"


def img(path: Path, alt: str) -> str:
    return f'<img src="{esc(path.relative_to(OUT.parent))}" alt="{esc(alt)}">'


def card(summary: dict, sample_dir: Path) -> str:
    total = int(summary["cells_total"])
    malignant = int(summary["malignant_epithelial_candidate_cells"])
    candidate = int(summary["epithelial_candidate_cells"])
    ambiguous = int(summary["epithelial_like_ambiguous_cells"])
    ref = int(summary["reference_cells"])
    cnv_rows, cnv_cols = summary["cnv_matrix_shape"]
    heatmap = sample_dir / "figures" / f"{summary['sample']}_pyinfercnv_i3_heatmap.png"
    cnv_summary = sample_dir / "figures" / f"{summary['sample']}_epithelial_candidate_cnv_summary.png"
    return f"""
    <section class="sample">
      <div class="sample-head"><h2>{esc(summary['sample'])}</h2><span class="status">{esc(summary['status'])}</span></div>
      <div class="metrics">
        <div><b>{total:,}</b><span>total cells</span></div>
        <div><b>{candidate:,}</b><span>epithelial candidates ({pct(candidate, total)})</span></div>
        <div><b>{malignant:,}</b><span>high-confidence malignant ({pct(malignant, total)})</span></div>
        <div><b>{ambiguous:,}</b><span>epithelial-like ambiguous ({pct(ambiguous, total)})</span></div>
        <div><b>{ref:,}</b><span>reference cells ({pct(ref, total)})</span></div>
      </div>
      <table>
        <tr><th>Backend</th><td>{esc(summary['backend'])} {esc(summary['pyinfercnv_version'])}</td><th>HMM</th><td>{esc(summary['hmm_type'])}</td></tr>
        <tr><th>Genes with hg38 positions</th><td>{int(summary['genes_with_positions']):,}</td><th>CNV matrix</th><td>{int(cnv_rows):,} cells x {int(cnv_cols):,} genomic bins</td></tr>
        <tr><th>CNV regions</th><td>{esc(summary['cnv_regions_present'])}</td><th>Definition</th><td>{esc(summary['candidate_definition'])}</td></tr>
      </table>
      <div class="figures">
        <figure>{img(heatmap, summary['sample'] + ' i3 CNV heatmap')}<figcaption>OmicVerse py-inferCNV i3 heatmap</figcaption></figure>
        <figure>{img(cnv_summary, summary['sample'] + ' epithelial candidate CNV summary')}<figcaption>Expanded epithelial-candidate CNV summary</figcaption></figure>
      </div>
    </section>
    """


def main() -> None:
    summaries = []
    for sample in ("SC000895-R1", "SC000895-R7"):
        sample_dir = RESULTS / sample
        with (sample_dir / "i3_summary.json").open() as handle:
            summaries.append((json.load(handle), sample_dir))

    generated = "2026-07-23"
    body = "\n".join(card(summary, sample_dir) for summary, sample_dir in summaries)
    html_text = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>R1/R7 OmicVerse py-inferCNV i3 v2</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 32px auto; max-width: 1500px; color: #20252b; line-height: 1.45; }}
h1 {{ margin-bottom: 4px; }} h2 {{ margin: 0; }}
.muted {{ color: #66717d; }} .note {{ background: #f3f6f8; border-left: 4px solid #3b6f8f; padding: 12px 16px; margin: 20px 0; }}
.sample {{ border-top: 1px solid #cfd7dd; padding-top: 24px; margin-top: 30px; }}
.sample-head {{ display: flex; align-items: center; justify-content: space-between; }}
.status {{ color: #176b42; border: 1px solid #8bc6a7; padding: 3px 10px; border-radius: 4px; font-weight: 600; }}
.metrics {{ display: grid; grid-template-columns: repeat(5, 1fr); gap: 10px; margin: 16px 0; }}
.metrics div {{ background: #f7f8f9; border: 1px solid #e1e6ea; padding: 12px; min-height: 58px; }}
.metrics b, .metrics span {{ display: block; }} .metrics b {{ font-size: 20px; }} .metrics span {{ color: #5f6b76; font-size: 13px; }}
table {{ border-collapse: collapse; width: 100%; margin: 12px 0 20px; }} th, td {{ border: 1px solid #dfe4e8; padding: 7px 9px; text-align: left; }} th {{ background: #f7f8f9; width: 20%; }}
.figures {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }} figure {{ margin: 0; }} figure img {{ width: 100%; height: auto; border: 1px solid #dfe4e8; }} figcaption {{ color: #5f6b76; font-size: 13px; margin-top: 5px; }}
li {{ margin: 5px 0; }} code {{ background: #f3f6f8; padding: 1px 4px; }}
@media (max-width: 900px) {{ .metrics, .figures {{ grid-template-columns: 1fr 1fr; }} }}
@media (max-width: 600px) {{ .metrics, .figures {{ grid-template-columns: 1fr; }} body {{ margin: 16px; }} }}
</style></head><body>
<h1>R1/R7 OmicVerse py-inferCNV i3 v2</h1>
<p class="muted">本地复核报告，生成日期 {generated}。R1：未复发原发灶；R7：复发原发灶。</p>
<div class="note"><b>本版目的：</b>在上一版 <code>malignant_only</code> 过于保守的基础上，使用 <code>epithelial_v2</code> 重新运行 i3。候选组包含 malignant、luminal、basal 三类明确上皮标签，以及同时具备上皮上下文和上皮 marker 信号的 ambiguous cluster。<b>epithelial candidate 不等于恶性细胞。</b>真正的恶性判定仍需结合 CNV 模式、H&E、marker 和空间位置。</div>
<h2>如何阅读</h2>
<ul>
 <li><b>High-confidence malignant：</b>上一版 marker 规则直接标出的 <code>malignant_epithelial_candidate</code>。</li>
 <li><b>Expanded epithelial candidates：</b>本版用于 CNV 筛查的扩大上皮集合；R7 即使没有高置信 malignant 标签，也可以在该集合中检查克隆性 CNV。</li>
 <li><b>Reference cells：</b>主要使用 endothelial、fibroblast/stromal 和 smooth-muscle/pericyte 细胞；参考细胞不代表绝对正常组织，只作为相对表达基线。</li>
 <li><b>i3：</b>OmicVerse 的 py-inferCNV HMM 三状态模式。热图中的大尺度、连续、跨多个相邻 genomic bins 的改变比单个 bin 的噪声更值得关注。</li>
</ul>
{body}
<h2>本版结论与下一步</h2>
<ol>
 <li>R1 和 R7 的 i3 计算均成功，hg38 gene-position 映射和 CNV/HMM 矩阵均已生成。</li>
 <li>R7 现在有 63,130 个扩大上皮候选细胞可用于 CNV 检查；这解决了“高置信 malignant 标签为 0 时无法运行候选组汇总”的技术限制，但不等同于证明 R7 没有恶性细胞。</li>
 <li>下一步应按 cluster/空间区域汇总 CNV，再与 H&E、AR/KRT8/KRT18/EPCAM、KRT5/KRT14、PTEN/TP53 等 marker 和 R1/R7 对照关系联合判断，而不是按单细胞逐个下结论。</li>
 <li>建议保留本版作为敏感性分析，并把 <code>malignant_only</code> 与 <code>epithelial_v2</code> 两版结果并列报告。</li>
</ol>
<p class="muted">源文件：本地下载的 JSON summary 与 PNG；远端完整 H5AD 位于服务器结果目录，未下载到本地报告包。</p>
</body></html>"""
    OUT.write_text(html_text, encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
