#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from sp_utils import (
    PILOT_SAMPLE,
    PROJECT_HOME,
    PROJECT_WORK,
    attach_spatial_from_10x,
    load_markers,
    preprocess_basic,
    qc_basic,
    read_10x_h5,
    read_manifest,
    save_basic_plots,
    write_json,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", default=PILOT_SAMPLE)
    args = parser.parse_args()

    sample = next(r for r in read_manifest(PROJECT_HOME / "config" / "sample_manifest.csv") if r["sample_id"] == args.sample_id)
    outdir = PROJECT_HOME / "results" / "pilot" / args.sample_id
    workdir = PROJECT_WORK / args.sample_id
    outdir.mkdir(parents=True, exist_ok=True)
    workdir.mkdir(parents=True, exist_ok=True)

    markers = load_markers(PROJECT_HOME / "config" / "prostate_markers.csv")
    marker_flat = [g for genes in markers.values() for g in genes]

    report = {"sample_id": args.sample_id, "warnings": []}
    adata = read_10x_h5(Path(sample["matrix_008um"]), args.sample_id)
    report["warnings"].extend(attach_spatial_from_10x(adata, Path(sample["matrix_008um"]).parent))
    qc_basic(adata)
    adata = preprocess_basic(adata)
    report["warnings"].extend(save_basic_plots(adata, outdir, f"{args.sample_id}_008um", marker_flat))
    adata.write_h5ad(workdir / f"{args.sample_id}_008um_pilot.h5ad")

    seg = read_10x_h5(Path(sample["segmented_h5"]), args.sample_id)
    qc_basic(seg)
    seg = preprocess_basic(seg)
    report["warnings"].extend(save_basic_plots(seg, outdir, f"{args.sample_id}_segmented", marker_flat))
    seg.write_h5ad(workdir / f"{args.sample_id}_segmented_pilot.h5ad")

    write_json(outdir / "pilot_report.json", report)
    print(f"Pilot complete for {args.sample_id}")


if __name__ == "__main__":
    main()
