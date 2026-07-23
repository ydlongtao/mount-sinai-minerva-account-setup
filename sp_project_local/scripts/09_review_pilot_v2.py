#!/usr/bin/env python3
from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

import pandas as pd

from sp_v2_utils import PROJECT_HOME, RESULTS_V2, load_config


def table(path: Path) -> str:
    if not path.exists():
        return f"<p class='missing'>Missing: {html.escape(str(path))}</p>"
    return pd.read_csv(path).to_html(index=False, classes="data", border=0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=PROJECT_HOME / "config" / "pilot_v2.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    sample_id = config["sample_id"]
    root = RESULTS_V2 / sample_id
    qc = root / "qc"
    model = root / "model"
    summary = {}
    for name, path in [("qc", qc / "qc_summary.json"), ("model", model / "model_summary.json")]:
        if path.exists():
            summary[name] = json.loads(path.read_text())
    title = f"Visium HD pilot v2 review: {sample_id}"
    sections = [
        f"<h1>{html.escape(title)}</h1>",
        "<p>This report is a QC and workflow review, not a final biological annotation.</p>",
        "<h2>Stage summaries</h2>",
        f"<pre>{html.escape(json.dumps(summary, indent=2, sort_keys=True))}</pre>",
        "<h2>QC quantiles</h2>", table(qc / "qc_quantiles.csv"),
        "<h2>Graph diagnostics</h2>", table(model / "graph_diagnostics.csv"),
        "<h2>Cluster resolution comparison</h2>", table(model / "clustering_resolution_summary.csv"),
        "<h2>Marker availability</h2>", table(model / "marker_availability.csv"),
    ]
    output = root / "pilot_v2_review.html"
    output.parent.mkdir(parents=True, exist_ok=True)
    content = """<!doctype html><html><head><meta charset="utf-8"><title>"""
    content += html.escape(title)
    content += """</title><style>
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;max-width:1400px;margin:32px auto;padding:0 24px;color:#17202a}
h1,h2{color:#174a5b}.data{border-collapse:collapse;font-size:13px;margin:12px 0 28px}
.data th,.data td{border:1px solid #c9d4d8;padding:5px 8px;text-align:left}.data th{background:#eaf1f3}
pre{background:#f5f7f8;padding:14px;overflow:auto}.missing{color:#a33}
</style></head><body>"""
    content += "\n".join(sections) + "</body></html>\n"
    output.write_text(content)
    print(output)


if __name__ == "__main__":
    main()
