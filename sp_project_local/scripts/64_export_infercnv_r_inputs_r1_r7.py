#!/usr/bin/env python3
"""Export counts, annotations and GRCh38 gene order for R inferCNV i3."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import shutil
from pathlib import Path

import infercnvpy as cnv
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmwrite

PROJECT = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
ROOT = PROJECT / "results" / "segmented_official_v1"
GTF = Path(os.environ.get("SP_HG38_REF", "/sc/arion/work/huangl21/references/hg38/gencode38")) / "gencode.v38.primary_assembly.annotation.gtf.gz"
INPUT_DIR = ROOT / "analysis_coarse" / "samples"
REF_TABLE = ROOT / "infercnv_reference_r1_r7" / "reference_cells_for_infercnvpy.csv.gz"
OUT_ROOT = ROOT / "infercnv_i3_inputs_r1_r7"
MALIGNANT_TYPES = {"malignant_epithelial_candidate"}


def marker_labels(sample: str, obs_names: pd.Index) -> np.ndarray:
    marker_path = ROOT / "preannotation_prostate" / "samples" / sample / "cluster_prostate_marker_spatial_preannotation.csv"
    marker = pd.read_csv(marker_path, dtype={"cluster": str}).set_index("cluster")
    cluster = pd.Series(obs_names.astype(str), index=obs_names)
    # The coarse analysis stores the cluster label used by the marker table.
    return np.array(["malignant_epithelial_candidate" if marker["predicted_type"].get(str(c), "") in MALIGNANT_TYPES else "other" for c in []], dtype=object)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True, choices=["SC000895-R1", "SC000895-R7"])
    p.add_argument("--output-root", default=str(OUT_ROOT))
    args = p.parse_args()
    sample = args.sample
    out = Path(args.output_root) / sample
    out.mkdir(parents=True, exist_ok=True)
    input_path = INPUT_DIR / sample / f"{sample}_official_cell_level_analysis_coarse.h5ad"
    adata = sc.read_h5ad(input_path)
    adata.obs_names = adata.obs_names.astype(str)
    if "counts" not in adata.layers:
        raise RuntimeError(f"{sample}: counts layer missing")

    selected = pd.read_csv(REF_TABLE)
    ref_ids = set(selected.loc[(selected.sample_id == sample) &
                               (selected.reference_group == "main_reference") &
                               selected.selected_for_pilot.astype(bool), "cell_id"].astype(str))
    marker_path = ROOT / "preannotation_prostate" / "samples" / sample / "cluster_prostate_marker_spatial_preannotation.csv"
    marker = pd.read_csv(marker_path, dtype={"cluster": str}).set_index("cluster")
    clusters = adata.obs["cell_level_leiden_coarse"].astype(str)
    malignant_clusters = set(marker.index[marker["predicted_type"].eq("malignant_epithelial_candidate")])
    malignant = clusters.isin(malignant_clusters).to_numpy()
    reference = adata.obs_names.isin(ref_ids)
    labels = np.where(reference, "main_reference", np.where(malignant, "malignant_epithelial_candidate", "other"))

    # Add gene positions using the same GENCODE v38 annotation validated for infercnvpy.
    cnv.io.genomic_position_from_gtf(GTF, adata=adata, gtf_gene_id="gene_name", inplace=True)
    valid = adata.var["chromosome"].notna().to_numpy()
    valid &= ~adata.var["chromosome"].isin(["chrX", "chrY"]).to_numpy()
    adata = adata[:, valid].copy()
    counts = adata.layers["counts"]
    if not hasattr(counts, "tocsr"):
        from scipy import sparse
        counts = sparse.csr_matrix(counts)
    else:
        counts = counts.tocsr()
    mm_path = out / "counts.mtx"
    mmwrite(mm_path, counts)
    with open(mm_path, "rb") as src, gzip.open(str(mm_path) + ".gz", "wb") as dst:
        shutil.copyfileobj(src, dst)
    mm_path.unlink()

    with gzip.open(out / "genes.tsv.gz", "wt") as fh:
        for gene in adata.var_names.astype(str):
            fh.write(gene + "\n")
    with gzip.open(out / "cells.tsv.gz", "wt") as fh:
        for cell in adata.obs_names.astype(str):
            fh.write(cell + "\n")
    with gzip.open(out / "annotations.tsv.gz", "wt") as fh:
        for cell, label in zip(adata.obs_names.astype(str), labels):
            fh.write(f"{cell}\t{label}\n")
    with gzip.open(out / "gene_order.tsv.gz", "wt") as fh:
        for gene, row in adata.var.iterrows():
            fh.write(f"{gene}\t{row['chromosome']}\t{int(row['start'])}\t{int(row['end'])}\n")
    report = {
        "sample": sample, "cells": int(adata.n_obs), "genes": int(adata.n_vars),
        "main_reference_cells": int((labels == "main_reference").sum()),
        "malignant_epithelial_candidate_cells": int((labels == "malignant_epithelial_candidate").sum()),
        "other_cells": int((labels == "other").sum()), "excluded_chrX_chrY": True,
        "gtf": str(GTF), "counts_input": str(out / "counts.mtx.gz"),
    }
    (out / "export_summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
