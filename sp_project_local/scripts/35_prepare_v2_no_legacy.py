#!/usr/bin/env python3
"""Create an auditable 8-sample manifest excluding legacy Visium samples."""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path


PROJECT_HOME = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
RAW_MANIFEST = PROJECT_HOME / "config" / "sample_manifest.csv"
V2_MANIFEST = PROJECT_HOME / "config" / "sample_manifest_v2_no_legacy.csv"
V2_RESULTS = PROJECT_HOME / "results" / "batch_v2_no_legacy"
EXCLUDED = {"TD006859-B408", "TD006859-B573"}


def main() -> None:
    with RAW_MANIFEST.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    excluded = [row for row in rows if row["sample_id"] in EXCLUDED]
    kept = [row for row in rows if row["sample_id"] not in EXCLUDED]
    if len(kept) != 8 or {row["sample_id"] for row in excluded} != EXCLUDED:
        raise RuntimeError("Unexpected manifest scope; refusing to write v2 manifest")

    for index, row in enumerate(kept, 1):
        row["sample_index"] = str(index)
        if not Path(row["matrix_016um"]).is_file():
            raise FileNotFoundError(row["matrix_016um"])

    V2_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with V2_MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=kept[0].keys())
        writer.writeheader()
        writer.writerows(kept)

    scope = {
        "analysis_version": "batch_v2_no_legacy",
        "included_sample_ids": [row["sample_id"] for row in kept],
        "excluded_sample_ids": sorted(EXCLUDED),
        "excluded_reason": "Legacy Visium chemistry/format differs from the SC000895 Visium HD cohort.",
        "source_manifest": str(RAW_MANIFEST),
        "manifest": str(V2_MANIFEST),
        "raw_data_read_only": True,
    }
    out = V2_RESULTS / "config" / "scope_manifest.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(scope, indent=2) + "\n")
    print(json.dumps(scope, indent=2))


if __name__ == "__main__":
    main()
