#!/usr/bin/env python3
"""Prepare an isolated corrected workspace for the Cellpose mpp pilot."""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path


PROJECT_HOME = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
SOURCE_MANIFEST = PROJECT_HOME / "config" / "sample_manifest_v2_no_legacy.csv"
MANIFEST = PROJECT_HOME / "config" / "sample_manifest_v2_corrected.csv"
RESULTS = PROJECT_HOME / "results" / "batch_v2_corrected"
WORK = PROJECT_HOME / "work" / "batch_v2_corrected"
SOURCE_RESULTS = PROJECT_HOME / "results" / "batch"
PILOT = {"SC000895-R4", "SC000895-R5", "SC000895-R9"}


def main() -> None:
    rows = list(csv.DictReader(SOURCE_MANIFEST.open(newline="")))
    with MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)
    (PROJECT_HOME / "config" / "batch_analysis_v2_corrected.yaml").write_text(json.dumps({
        "analysis_version": "batch_v2_corrected_mpp_pilot",
        "sample_ids": [row["sample_id"] for row in rows],
        "n_samples": len(rows),
        "mpp_mode": "source",
        "crop_defaults": {"x_range_um": [2000, 3000], "y_range_um": [2000, 3000], "mpp": 0.3239732935, "mpp_mode": "source", "buffer": 150, "expand_max_bin_distance": 4},
        "crop_by_sample": {},
        "model": {"n_top_genes": 3000, "n_pcs": 30, "n_neighbors": 15, "leiden_resolution": 0.4, "random_state": 0}
    }, indent=2) + "\n")
    for row in rows:
        sample = row["sample_id"]
        target = RESULTS / "samples" / sample
        target.mkdir(parents=True, exist_ok=True)
        source = SOURCE_RESULTS / "samples" / sample
        for name in ("002um", "008um", "016um"):
            link = target / name
            if not link.exists(): link.symlink_to(source / name, target_is_directory=True)
    scope = {"analysis_version": "batch_v2_corrected_mpp_pilot", "pilot_samples": sorted(PILOT), "all_samples": [r["sample_id"] for r in rows], "cellpose_mpp": "source_mpp derived from morphology image and barcode coordinate extent", "source_results": str(SOURCE_RESULTS), "results": str(RESULTS)}
    out = RESULTS / "config" / "corrected_scope.json"; out.parent.mkdir(parents=True, exist_ok=True); out.write_text(json.dumps(scope, indent=2) + "\n")
    print(json.dumps(scope, indent=2))


if __name__ == "__main__": main()
