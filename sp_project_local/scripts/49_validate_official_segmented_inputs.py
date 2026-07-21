#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import h5py


PROJECT_HOME = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
MANIFEST = Path(os.environ.get(
    "SP_PROJECT_MANIFEST",
    PROJECT_HOME / "config" / "sample_manifest_v2_corrected.csv",
))
RESULTS = Path(os.environ.get(
    "SP_OFFICIAL_RESULTS",
    PROJECT_HOME / "results" / "segmented_official_v1",
))


def main() -> None:
    rows = list(csv.DictReader(MANIFEST.open(newline="")))
    if len(rows) != 8:
        raise ValueError(f"Expected 8 retained samples, found {len(rows)} in {MANIFEST}")

    report = []
    errors = []
    for row in rows:
        sample = row["sample_id"]
        outs = Path(row["outs_dir"])
        segmented = outs / "segmented_outputs"
        files = {
            "filtered_cell_h5": Path(row["segmented_h5"]),
            "raw_cell_h5": segmented / "raw_feature_cell_matrix.h5",
            "cell_geojson": segmented / "cell_segmentations.geojson",
            "nucleus_geojson": segmented / "nucleus_segmentations.geojson",
            "hires_image": segmented / "spatial" / "tissue_hires_image.png",
            "scalefactors": segmented / "spatial" / "scalefactors_json.json",
            "barcode_mappings": outs / "barcode_mappings.parquet",
            "matrix_008um": Path(row["matrix_008um"]),
        }
        missing = [name for name, path in files.items() if not path.is_file()]
        item = {
            "sample_id": sample,
            "status": "fail" if missing else "pass",
            "missing": missing,
            "files": {name: str(path) for name, path in files.items()},
        }
        if not missing:
            with h5py.File(files["filtered_cell_h5"], "r") as handle:
                shape = tuple(int(x) for x in handle["matrix/shape"][:])
                item["genes"] = shape[0]
                item["cells"] = shape[1]
                item["matrix_nnz"] = int(handle["matrix/data"].shape[0])
            sf = json.loads(files["scalefactors"].read_text())
            item["microns_per_pixel"] = float(sf["microns_per_pixel"])
            item["tissue_hires_scalef"] = float(sf["tissue_hires_scalef"])
        else:
            errors.append(f"{sample}: missing {', '.join(missing)}")
        report.append(item)

    out = RESULTS / "validation"
    out.mkdir(parents=True, exist_ok=True)
    (out / "official_segmented_input_validation.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n"
    )
    with (out / "official_segmented_input_validation.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "sample_id", "status", "cells", "genes", "matrix_nnz",
            "microns_per_pixel", "tissue_hires_scalef", "missing",
        ])
        writer.writeheader()
        for item in report:
            writer.writerow({
                "sample_id": item["sample_id"],
                "status": item["status"],
                "cells": item.get("cells", ""),
                "genes": item.get("genes", ""),
                "matrix_nnz": item.get("matrix_nnz", ""),
                "microns_per_pixel": item.get("microns_per_pixel", ""),
                "tissue_hires_scalef": item.get("tissue_hires_scalef", ""),
                "missing": ";".join(item["missing"]),
            })
    if errors:
        raise RuntimeError("; ".join(errors))
    print(json.dumps({"status": "pass", "samples": len(report), "results": str(out)}, indent=2))


if __name__ == "__main__":
    main()
