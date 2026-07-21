#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import pandas as pd

from batch_utils import BATCH_RESULTS, MANIFEST, read_manifest, write_json


def h5_shape(path: str) -> list[int | None]:
    p = Path(path)
    if not p.exists():
        return [None, None]
    with h5py.File(p, "r") as h5:
        shape = h5["matrix/shape"][:]
        return [int(shape[1]), int(shape[0])]


def main() -> None:
    parser = argparse.ArgumentParser(description="Collect batch sample QC from Space Ranger outputs.")
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--output-dir", type=Path, default=BATCH_RESULTS / "qc")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for sample in read_manifest(args.manifest):
        row = {
            "sample_index": int(sample["sample_index"]),
            "sample_id": sample["sample_id"],
            "group": sample.get("group", ""),
            "metrics_summary": sample.get("metrics_summary", ""),
            "web_summary": sample.get("web_summary", ""),
        }
        for key, col in [("002um", "matrix_002um"), ("008um", "matrix_008um"), ("016um", "matrix_016um"), ("segmented", "segmented_h5")]:
            shape = h5_shape(sample[col])
            row[f"{key}_bins_or_cells"] = shape[0]
            row[f"{key}_genes"] = shape[1]
            row[f"{key}_exists"] = Path(sample[col]).exists()
        metrics_path = Path(sample["metrics_summary"])
        if metrics_path.exists():
            try:
                metrics = pd.read_csv(metrics_path)
                if metrics.shape[1] >= 2:
                    metric_map = dict(zip(metrics.iloc[:, 0].astype(str), metrics.iloc[:, 1].astype(str)))
                    for name in ["Number of bins", "Mean reads per bin", "Median genes per bin"]:
                        if name in metric_map:
                            row[name] = metric_map[name]
            except Exception as exc:
                row["metrics_read_error"] = repr(exc)
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(args.output_dir / "sample_qc_summary.csv", index=False)
    write_json(args.output_dir / "sample_qc_summary.json", rows)
    html = "<html><body><h1>Batch QC Summary</h1>" + df.to_html(index=False) + "</body></html>"
    (args.output_dir / "sample_qc_report.html").write_text(html)
    print(json.dumps({"samples": len(rows), "output": str(args.output_dir)}, indent=2))


if __name__ == "__main__":
    main()
