#!/usr/bin/env python3
"""Build a complete local R1/R7 analysis report through Step 5."""
from __future__ import annotations

import csv
import html
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "local_results" / "segmented_official_v1"
STEP5 = BASE / "step5_niche_cnv_pseudobulk_r1_r7"
OUT = BASE / "r1_r7_complete_analysis_report.html"


def e(x: object) -> str:
    return html.escape(str(x))


def table(path: Path, columns: list[str] | None = None, limit: int | None = None) -> str:
    if not path.exists():
        return f"<p>Missing: <code>{e(path.name)}</code></p>"
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if columns:
        rows = [{k: row.get(k, "") for k in columns} for row in rows]
    if limit:
        rows = rows[:limit]
    if not rows:
        return "<p>Empty table</p>"
    cols = list(rows[0])
    head = "".join(f"<th>{e(x)}</th>" for x in cols)
    body = "".join("<tr>" + "".join(f"<td>{e(row.get(x, ''))}</td>" for x in cols) + "</tr>" for row in rows)
    return f"<table><tr>{head}</tr>{body}</table>"


def img(path: Path, alt: str) -> str:
    if not path.exists():
        return f"<p>Missing image: <code>{e(path.name)}</code></p>"
    return f'<img src="{e(path.relative_to(OUT.parent))}" alt="{e(alt)}">'


def main() -> None:
    step5 = json.loads((STEP5 / "step5_summary.json").read_text())
    final = json.loads((BASE / "final_epithelial_groups_r1_r7" / "final_epithelial_groups_summary.json").read_text())
    step4 = json.loads((BASE / "step4_final_groups_spatial_niche_r1_r7" / "step4_summary.json").read_text())
    sample_cards = []
    for row in step5["samples"]:
        sample = row["sample"]
        top_path = STEP5 / f"{sample}_malignant_vs_normal_top100.csv"
        top_html = table(top_path, ["gene", "malignant_cpm", "normal_cpm", "log2fc_malignant_vs_normal"], 20)
        final_row = next(x for x in final["samples"] if x["sample"] == sample)
        niche_row = next(x for x in step4["samples"] if x["sample"] == sample)
        sample_cards.append(f"""
        <section class="sample"><h2>{e(sample)}</h2>
          <div class="metrics">
            <div><b>{row['malignant_cells']:,}</b><span>malignant epithelial</span></div>
            <div><b>{row['normal_cells']:,}</b><span>normal epithelial</span></div>
            <div><b>{niche_row['valid_niche_cells']:,}</b><span>valid niche cells</span></div>
            <div><b>{row['mean_cnv_burden_malignant']:.4f}</b><span>malignant mean abs CNV</span></div>
            <div><b>{row['mean_cnv_burden_normal']:.4f}</b><span>normal mean abs CNV</span></div>
          </div>
          <p><b>Top malignant-up genes:</b> <code>{e(', '.join(row['top_up_genes'][:12]))}</code></p>
          <p><b>Top malignant-down genes:</b> <code>{e(', '.join(row['top_down_genes'][:12]))}</code></p>
          <p>最终分组：<code>{e(json.dumps(final_row['final_group_counts'], ensure_ascii=False))}</code></p>
          <h3>Top pseudobulk effect-size genes</h3>{top_html}
          <div class="grid2">
            <figure>{img(BASE / "final_epithelial_groups_r1_r7" / sample / f"{sample}_final_epithelial_groups_he_overlay.png", sample + " final epithelial grouping") }<figcaption>最终上皮分群 H&amp;E overlay</figcaption></figure>
            <figure>{img(BASE / "step4_final_groups_spatial_niche_r1_r7" / sample / f"{sample}_step4_niche_he_overlay.png", sample + " Step 4 niche overlay") }<figcaption>最终分群与空间 niche overlay</figcaption></figure>
          </div>
        </section>""")
    report = f"""<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>R1/R7 Complete Analysis</title>
<style>
body{{font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;color:#20252b;max-width:1600px;margin:28px auto;padding:0 18px}}h1{{margin-bottom:4px}}h2{{margin-top:28px}}.muted{{color:#68737d}}.note{{background:#f3f6f8;border-left:4px solid #3b6f8f;padding:12px 16px;margin:18px 0}}.warn{{background:#fff4cc;border-left:4px solid #c58b00;padding:12px 16px;margin:18px 0}}.sample{{border-top:1px solid #ccd5dc;margin-top:30px;padding-top:22px}}.metrics{{display:grid;grid-template-columns:repeat(5,1fr);gap:10px;margin:15px 0}}.metrics div{{background:#f7f8f9;border:1px solid #e1e6ea;padding:11px}}.metrics b,.metrics span{{display:block}}.metrics b{{font-size:19px}}.metrics span{{font-size:12px;color:#5f6b76}}table{{border-collapse:collapse;width:100%;font-size:12px;margin:12px 0 20px}}th,td{{border:1px solid #dfe4e8;padding:5px 7px;text-align:left}}th{{background:#f7f8f9}}.grid2,.grid3{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}}.grid3{{grid-template-columns:repeat(3,minmax(0,1fr))}}figure{{margin:0}}img{{width:100%;height:auto;border:1px solid #dfe4e8}}figcaption{{font-size:13px;color:#5f6b76;padding-top:5px}}code{{background:#f1f4f6;padding:2px 5px}}li{{margin:6px 0}}a{{color:#155e9c}}@media(max-width:900px){{.metrics,.grid2,.grid3{{grid-template-columns:1fr 1fr}}}}@media(max-width:600px){{.metrics,.grid2,.grid3{{grid-template-columns:1fr}};body{{margin:14px}}}}
</style></head><body>
<h1>R1/R7 前列腺癌空间转录组完整分析报告</h1>
<p class="muted">R1：未复发患者原发灶；R7：复发患者原发灶。报告覆盖最终上皮分群、i3 CNV、空间 niche、pseudobulk 表达和恶性程序比较。</p>
<div class="note"><b>数据与样本范围：</b>本报告只分析 R1 和 R7。原始数据保持只读；所有计算在 Minerva LSF compute node 完成；本地仅下载报告、图像、汇总表和压缩细胞级标签。</div>
<h2>一、完整处理流程</h2><ol>
<li><b>输入与结构验证：</b>读取 official segmented cell H5AD，确认空间坐标、counts layer、gene names 和 H&amp;E hires 图。</li>
<li><b>预注释：</b>基于前列腺腔面、基底、恶性候选、免疫、间质和内皮 marker 进行 cluster-level preliminary annotation。</li>
<li><b>OmicVerse py-inferCNV i3：</b>使用 main reference（endothelial、fibroblast/stromal、smooth muscle/pericyte）和 hg38 gene positions 进行 CNV 推断。</li>
<li><b>核心 panel + CNV 细胞级重分类：</b>结合 epithelial、luminal、malignant、proliferation、basal、neuroendocrine panel。</li>
<li><b>最终分群：</b><code>malignant_marker_candidate</code> 与 <code>malignant_marker_CNV_candidate</code> 合并为 <code>malignant_epithelial</code>；<code>epithelial_CNV_candidate</code>、<code>epithelial_ambiguous</code>、<code>basal_epithelial</code> 合并为 <code>normal_epithelial</code>。</li>
<li><b>Step 4 spatial niche：</b>以每个细胞最近 12 个邻居的细胞组成建立 pooled 5-niche 模型。</li>
<li><b>Step 5：</b>按 niche 汇总 malignant/normal 比例、CNV burden；进行样本内 pseudobulk malignant vs normal 表达效应量比较和程序评分。</li></ol>
<div class="warn"><b>统计限制：</b>R1/R7 各只有一个患者样本。pseudobulk 目前是每个样本、每个最终分组一个聚合矩阵，不能据此计算可靠的患者层面 FDR 或显著性；本报告使用 log2FC、CNV burden 和程序分数进行探索性比较。</div>
<h2>二、最终分群与样本结果</h2>{''.join(sample_cards)}
<h2>三、Niche 组成与 CNV burden</h2>
<p>按 Step 4 的 pooled 5-niche ID 汇总最终上皮分群。表中 fraction of niche 表示该组占 niche 全部细胞的比例；CNV burden 使用 <code>mean(abs(X_cnv))</code>。</p>
{table(STEP5 / "niche_epithelial_cnv_burden.csv")}
<div class="grid3"><figure>{img(STEP5 / "niche_epithelial_composition.png", "niche epithelial composition")}<figcaption>每个 niche 的 malignant/normal epithelial 比例</figcaption></figure><figure>{img(STEP5 / "niche_cnv_burden.png", "niche CNV burden")}<figcaption>每个 niche 的 CNV burden</figcaption></figure><figure>{img(STEP5 / "malignant_normal_program_scores.png", "program scores")}<figcaption>恶性与正常上皮程序分数</figcaption></figure></div>
<h2>四、Pseudobulk 与恶性相关程序</h2>
<p>每个样本分别将 malignant_epithelial 和 normal_epithelial 的 raw counts 聚合，计算 CPM 后的 log2 fold change。该结果用于候选基因排序，不等价于有生物学重复的差异检验。</p>
{table(STEP5 / "malignant_normal_program_scores.csv")}
<h2>五、CNV、marker 和空间结果的综合结论</h2><ol>
<li>R1 的 malignant_epithelial 为 6,406 个细胞，R7 为 13,578 个细胞；R7 的 marker-defined malignant epithelial 数量更多，但这本身不代表肿瘤负荷更高。</li>
<li>R1 的 malignant/normal mean abs CNV 分别为 {step5['samples'][0]['mean_cnv_burden_malignant']:.4f}/{step5['samples'][0]['mean_cnv_burden_normal']:.4f}；R7 分别为 {step5['samples'][1]['mean_cnv_burden_malignant']:.4f}/{step5['samples'][1]['mean_cnv_burden_normal']:.4f}。两个样本中 normal 组反而略高，说明当前 CNV burden metric 和最终 marker 分群之间不完全一致，必须进一步检查 CNV reference、CNV 区段连续性和细胞组成。</li>
<li>R1 与 R7 的 malignant-up pseudobulk genes 均包含 AMACR、GOLM1、FASN、TPD52、MYC 等肿瘤相关程序；R1 还出现 MCM2、MKI67、TOP2A、UBE2C 等增殖基因。</li>
<li>下一步不应直接把当前标签作为最终恶性诊断。应在 niche 内重新计算 CNV 的连续染色体区段、比较 malignant/normal 的 CNV state 分布，并将结果与 H&amp;E 腺体形态联合复核。</li></ol>
<h2>六、文件索引</h2><ul>
<li><a href="final_epithelial_groups_r1_r7/r1_r7_final_epithelial_groups_report.html">最终上皮分群报告</a></li>
<li><a href="step4_final_groups_spatial_niche_r1_r7/r1_r7_step4_final_groups_spatial_niche_report.html">Step 4 niche 报告</a></li>
<li><a href="omicverse_pyinfercnv_i3_epithelial_v2_r1_r7/r1_r7_omicverse_pyinfercnv_i3_epithelial_v2_report.html">OmicVerse i3 CNV v2 报告</a></li>
<li><a href="epithelial_cnv_corepanel_reclassification_r1_r7/r1_r7_epithelial_cnv_corepanel_reclassification_report.html">核心 panel + CNV 重分类报告</a></li></ul>
</body></html>"""
    OUT.write_text(report, encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
