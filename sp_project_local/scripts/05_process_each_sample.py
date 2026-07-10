#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd
import scanpy as sc

from sp_utils import (
    PROJECT_HOME,
    PROJECT_WORK,
    attach_spatial_from_10x,
    load_markers,
    preprocess_basic,
    qc_basic,
    read_10x_h5,
    sample_from_index,
    save_basic_plots,
    write_json,
)


def marker_table(adata, out_csv: Path) -> None:
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
    df = sc.get.rank_genes_groups_df(adata, group=None)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)


def process_one(sample: dict[str, str]) -> None:
    sample_id = sample["sample_id"]
    outdir = PROJECT_HOME / "results" / "samples" / sample_id
    workdir = PROJECT_WORK / sample_id
    outdir.mkdir(parents=True, exist_ok=True)
    workdir.mkdir(parents=True, exist_ok=True)
    markers = load_markers(PROJECT_HOME / "config" / "prostate_markers.csv")
    marker_flat = [g for genes in markers.values() for g in genes]
    report = {"sample_id": sample_id, "warnings": []}

    for label, key in [("008um", "matrix_008um"), ("segmented", "segmented_h5")]:
        adata = read_10x_h5(Path(sample[key]), sample_id)
        adata.obs["assay_resolution"] = label
        if label == "008um":
            report["warnings"].extend(attach_spatial_from_10x(adata, Path(sample[key]).parent))
        qc_basic(adata)
        adata = preprocess_basic(adata)
        save_basic_plots(adata, outdir, f"{sample_id}_{label}", marker_flat)
        marker_table(adata, outdir / f"{sample_id}_{label}_rank_genes_groups.csv")
        adata.write_h5ad(workdir / f"{sample_id}_{label}.h5ad")

    write_json(outdir / "sample_report.json", report)
    print(f"Processed {sample_id}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-index", type=int, default=int(os.environ.get("LSB_JOBINDEX", "0") or "0"))
    args = parser.parse_args()
    if args.sample_index < 1:
        raise SystemExit("Provide --sample-index or run as an LSF array job with LSB_JOBINDEX")
    sample = sample_from_index(args.sample_index, PROJECT_HOME / "config" / "sample_manifest.csv")
    process_one(sample)


if __name__ == "__main__":
    main()
