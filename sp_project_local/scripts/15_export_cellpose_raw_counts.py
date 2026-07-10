#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import omicverse as ov
import scanpy as sc
from scipy import sparse
from shapely import wkt
from shapely.geometry import mapping


def write_complete_cell_geojson(cdata: ad.AnnData, export_dir: Path) -> Path:
    geojson_path = export_dir / "graphclust_annotated_cell_segmentations.geojson"
    backup_path = export_dir / "graphclust_annotated_cell_segmentations.omicsverse_simple.geojson"
    if geojson_path.exists() and not backup_path.exists():
        geojson_path.replace(backup_path)

    features = []
    skipped = []
    for i, (idx, row) in enumerate(cdata.obs.iterrows(), start=1):
        try:
            geom = wkt.loads(str(row.get("geometry", "")))
        except Exception as exc:
            skipped.append({"cell": str(idx), "reason": f"parse_error:{exc!r}"})
            continue
        if geom.is_empty or not geom.is_valid:
            skipped.append({"cell": str(idx), "reason": "empty_or_invalid"})
            continue
        features.append({
            "type": "Feature",
            "properties": {
                "cell_id": i,
                "cellid": str(row.get("cellid", idx)),
                "object_id": str(row.get("object_id", idx)),
            },
            "geometry": mapping(geom),
        })

    geojson_path.write_text(json.dumps({"type": "FeatureCollection", "features": features}, separators=(",", ":")) + "\n")
    report = {
        "geojson": str(geojson_path),
        "features": len(features),
        "skipped": len(skipped),
        "skipped_examples": skipped[:10],
        "backup": str(backup_path) if backup_path.exists() else None,
    }
    (export_dir / "complete_geojson_report.json").write_text(json.dumps(report, indent=2) + "\n")
    return geojson_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Aggregate raw Visium HD 2um bin counts to Cellpose/bin2cell cells.")
    parser.add_argument("--labels-h5ad", type=Path, required=True)
    parser.add_argument("--cell-h5ad", type=Path, required=True)
    parser.add_argument("--raw-2um-h5", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    labels = sc.read_h5ad(args.labels_h5ad)
    cells = sc.read_h5ad(args.cell_h5ad)
    if "labels_joint" not in labels.obs:
        raise ValueError("labels_h5ad must contain obs['labels_joint']")
    if "geometry" not in cells.obs:
        raise ValueError("cell_h5ad must contain obs['geometry']")

    raw = sc.read_10x_h5(args.raw_2um_h5)
    raw.var_names_make_unique()

    missing_bins = labels.obs_names.difference(raw.obs_names)
    if len(missing_bins):
        raise ValueError(f"{len(missing_bins)} label bins are missing from raw 2um matrix")
    common_genes = labels.var_names.intersection(raw.var_names)
    if len(common_genes) == 0:
        raise ValueError("No overlapping genes between labels and raw 2um matrix")

    labels = labels[:, common_genes].copy()
    raw_subset = raw[labels.obs_names, common_genes].copy()
    raw_subset.X = raw_subset.X.astype(np.int32)

    label_ids = labels.obs["labels_joint"].astype(int).to_numpy()
    positive = label_ids > 0
    object_ids = cells.obs["object_id"].astype(int).to_numpy() if "object_id" in cells.obs else np.arange(1, cells.n_obs + 1)
    id_to_row = {int(obj): i for i, obj in enumerate(object_ids)}
    row_ids = np.array([id_to_row.get(int(x), -1) for x in label_ids[positive]], dtype=np.int64)
    keep = row_ids >= 0
    if not np.all(keep):
        positive_idx = np.flatnonzero(positive)
        positive[positive_idx[~keep]] = False
        row_ids = row_ids[keep]

    col_ids = np.flatnonzero(positive)
    weights = sparse.csr_matrix(
        (np.ones(len(row_ids), dtype=np.int32), (row_ids, col_ids)),
        shape=(cells.n_obs, raw_subset.n_obs),
    )
    raw_counts = (weights @ raw_subset.X).tocsr().astype(np.int32)

    raw_cells = ad.AnnData(X=raw_counts, obs=cells.obs.copy(), var=labels.var.copy(), uns=cells.uns.copy(), obsm={k: v.copy() for k, v in cells.obsm.items()})
    raw_cells.uns.setdefault("cellpose_raw_counts", {})
    raw_cells.uns["cellpose_raw_counts"].update({
        "raw_2um_h5": str(args.raw_2um_h5),
        "labels_h5ad": str(args.labels_h5ad),
        "cell_h5ad": str(args.cell_h5ad),
        "positive_labeled_bins": int(positive.sum()),
        "total_counts": int(raw_counts.sum()),
    })

    raw_cell_h5ad = args.output_dir / "SC000895-R4_cell_level_cellpose_raw_counts.h5ad"
    raw_cells.write_h5ad(raw_cell_h5ad, compression="lzf")

    export_dir = args.output_dir / "cellpose_spaceranger_raw_counts_output"
    ov.space.write_visium_hd_cellseg(raw_cells, str(export_dir))
    full_geojson = write_complete_cell_geojson(raw_cells, export_dir)

    summary = {
        "raw_cell_h5ad": str(raw_cell_h5ad),
        "export_dir": str(export_dir),
        "complete_geojson": str(full_geojson),
        "cell_shape": list(raw_cells.shape),
        "positive_labeled_bins": int(positive.sum()),
        "total_raw_counts": int(raw_counts.sum()),
        "max_raw_count": int(raw_counts.max()),
        "genes": int(raw_cells.n_vars),
        "cells": int(raw_cells.n_obs),
    }
    (args.output_dir / "raw_counts_cellseg_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
