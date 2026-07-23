#!/usr/bin/env python3
"""Build a local HTML report for Step 4 spatial niche analysis."""
from __future__ import annotations
import csv, html, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "local_results" / "segmented_official_v1" / "step4_final_groups_spatial_niche_r1_r7"
OUT = DATA / "r1_r7_step4_final_groups_spatial_niche_report.html"

def e(x): return html.escape(str(x))

def csv_table(path: Path) -> str:
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle))
    if not rows:
        return "<p>Composition table empty</p>"
    head = "".join(f"<th>{e(x)}</th>" for x in rows[0])
    body = "".join("<tr>" + "".join(f"<td>{e(x)}</td>" for x in row) + "</tr>" for row in rows[1:])
    return f"<table><tr>{head}</tr>{body}</table>"

def main():
    summary = json.loads((DATA / "step4_summary.json").read_text())
    sections = []
    for row in summary["samples"]:
        sample = row["sample"]
        comp = DATA / sample / f"{sample}_step4_niche_composition.csv"
        image = f"{sample}/{sample}_step4_niche_he_overlay.png"
        table = "<p>Composition table unavailable</p>"
        if comp.exists():
            table = csv_table(comp)
        sections.append(f"<section class='sample'><h2>{e(sample)}</h2><p><b>有效 niche 细胞：</b>{row['valid_niche_cells']:,}；<b>niche 数量：</b>{row['niches']}。最终分群计数：{e(json.dumps(row['final_group_counts'], ensure_ascii=False))}</p><h3>Niche 平均组成</h3>{table}<figure><img src='{e(image)}' alt='{e(sample)} Step 4 niche H&E overlay'><figcaption>背景颜色：niche；红色：malignant_epithelial；绿色空心：normal_epithelial。</figcaption></figure></section>")
    report = f"""<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>Step 4 R1/R7 spatial niche</title><style>body{{font:15px/1.55 -apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;color:#20252b;max-width:1500px;margin:28px auto;padding:0 18px}}.note{{background:#f3f6f8;border-left:4px solid #3b6f8f;padding:12px 16px;margin:18px 0}}.sample{{border-top:1px solid #ccd5dc;margin-top:30px;padding-top:22px}}table{{border-collapse:collapse;width:100%;font-size:13px}}th,td{{border:1px solid #dfe4e8;padding:6px 8px;text-align:left}}th{{background:#f7f8f9}}img{{width:100%;height:auto;border:1px solid #dfe4e8}}figcaption{{color:#5f6b76;font-size:13px;padding-top:5px}}code{{background:#f1f4f6;padding:2px 5px}}</style></head><body><h1>Step 4：最终上皮分群 × 空间 niche</h1><p>R1 为未复发原发灶，R7 为复发原发灶。niche 使用两样本 pooled 的局部 12-neighbor 组成并聚类为 5 个可比较空间状态。</p><div class='note'><b>解释限制：</b>每个样本目前只有一个患者，niche 结果用于空间模式和后续假设生成，不作患者层面统计推断。niche 是局部组成状态，不是病理诊断；红色 malignant_epithelial 仍需与连续 CNV 区段、H&E 和 marker 联合解释。</div>{''.join(sections)}<h2>下一步</h2><ol><li>按 niche 汇总 malignant/normal epithelial 的比例和 CNV burden。</li><li>进入 Step 5：在 malignant_epithelial 与 normal_epithelial 之间进行样本内 pseudobulk 差异表达。</li><li>随后比较 R1/R7 的 niche composition 和恶性上皮相关通路。</li></ol></body></html>"""
    OUT.write_text(report, encoding='utf-8'); print(OUT)

if __name__ == '__main__': main()
