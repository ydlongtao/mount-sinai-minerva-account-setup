#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

from sp_utils import PROJECT_HOME, h5_matrix_shape, read_manifest


def read_metrics(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    try:
        df = pd.read_csv(path)
        if df.shape[1] >= 2:
            return {str(k): str(v) for k, v in zip(df.iloc[:, 0], df.iloc[:, 1])}
    except Exception:
        pass
    return {}


def main() -> None:
    manifest = PROJECT_HOME / "config" / "sample_manifest.csv"
    rows = []
    for sample in read_manifest(manifest):
        metrics = read_metrics(Path(sample["metrics_summary"]))
        row = {
            "sample_index": sample["sample_index"],
            "sample_id": sample["sample_id"],
            "group": sample["group"],
            "metrics_summary_exists": Path(sample["metrics_summary"]).exists(),
            "web_summary_exists": Path(sample["web_summary"]).exists(),
        }
        for key in ["matrix_008um", "matrix_016um", "segmented_h5"]:
            path = Path(sample[key])
            row[f"{key}_exists"] = path.exists()
            if path.exists():
                genes, barcodes_or_nnz = h5_matrix_shape(path)
                row[f"{key}_shape_0"] = genes
                row[f"{key}_shape_1"] = barcodes_or_nnz
                row[f"{key}_size_bytes"] = path.stat().st_size
        for metric_key in [
            "Estimated Number of Cells",
            "Mean Reads per Cell",
            "Median Genes per Cell",
            "Number of Reads",
            "Valid Barcodes",
            "Sequencing Saturation",
        ]:
            row[metric_key] = metrics.get(metric_key, "")
        rows.append(row)

    out = PROJECT_HOME / "results" / "sample_qc_summary.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted({k for row in rows for k in row}))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote QC summary to {out}")


if __name__ == "__main__":
    main()
