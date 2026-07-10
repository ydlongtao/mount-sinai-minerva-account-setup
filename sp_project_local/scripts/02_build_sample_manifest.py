#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

from sp_utils import PROJECT_HOME, RAW_ROOT, sample_group


def main() -> None:
    rows = []
    for outs in sorted(RAW_ROOT.glob("**/outs")):
        sample_dir = outs.parent
        sample_id = sample_dir.name
        if not (sample_id.startswith("SC000895") or sample_id.startswith("TD006859")):
            continue
        rows.append(
            {
                "sample_index": len(rows) + 1,
                "sample_id": sample_id,
                "group": sample_group(sample_id),
                "sample_dir": str(sample_dir),
                "outs_dir": str(outs),
                "matrix_002um": str(outs / "binned_outputs" / "square_002um" / "filtered_feature_bc_matrix.h5"),
                "matrix_008um": str(outs / "binned_outputs" / "square_008um" / "filtered_feature_bc_matrix.h5"),
                "matrix_016um": str(outs / "binned_outputs" / "square_016um" / "filtered_feature_bc_matrix.h5"),
                "segmented_h5": str(outs / "segmented_outputs" / "filtered_feature_cell_matrix.h5"),
                "metrics_summary": str(outs / "metrics_summary.csv"),
                "web_summary": str(outs / "web_summary.html"),
                "spatial_dir": str(outs / "binned_outputs" / "square_008um" / "spatial"),
                "segmented_spatial_dir": str(outs / "segmented_outputs" / "spatial"),
            }
        )
    if not rows:
        raise SystemExit(f"No expected samples found under {RAW_ROOT}")
    out = PROJECT_HOME / "config" / "sample_manifest.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} samples to {out}")


if __name__ == "__main__":
    main()
