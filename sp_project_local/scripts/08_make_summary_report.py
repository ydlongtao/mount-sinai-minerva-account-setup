#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from html import escape

import pandas as pd

from sp_utils import PROJECT_HOME, read_manifest


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT_HOME))
    except Exception:
        return str(path)


def main() -> None:
    results = PROJECT_HOME / "results"
    rows = read_manifest(PROJECT_HOME / "config" / "sample_manifest.csv")
    qc = results / "sample_qc_summary.csv"
    marker = results / "prostate_marker_scores.csv"
    lines = [
        "# Visium HD Prostate Cancer Summary Report",
        "",
        f"- Project: `{PROJECT_HOME}`",
        f"- Samples: {len(rows)}",
        f"- QC summary: `{rel(qc)}`" if qc.exists() else "- QC summary: missing",
        f"- Marker scores: `{rel(marker)}`" if marker.exists() else "- Marker scores: missing",
        "",
        "## Samples",
        "",
    ]
    for row in rows:
        sample_id = row["sample_id"]
        lines.append(f"### {sample_id}")
        sample_dir = results / "samples" / sample_id
        for image in sorted(sample_dir.glob("*.png")):
            lines.append(f"- `{rel(image)}`")
        lines.append("")
    if qc.exists():
        df = pd.read_csv(qc)
        lines += ["## QC Table Preview", "", "```csv", df.head(20).to_csv(index=False).strip(), "```", ""]
    md = results / "summary_report.md"
    html = results / "summary_report.html"
    md.write_text("\n".join(lines) + "\n")
    body = "\n".join(f"<pre>{escape(line)}</pre>" for line in lines)
    html.write_text(f"<html><body>{body}</body></html>\n")
    print(f"Wrote {md}")
    print(f"Wrote {html}")


if __name__ == "__main__":
    main()
