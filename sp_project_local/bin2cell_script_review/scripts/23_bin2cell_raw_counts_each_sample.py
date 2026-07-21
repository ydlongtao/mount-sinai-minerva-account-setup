#!/usr/bin/env python3
"""中文说明：将 Cellpose 标签映射到 2um bins 并聚合 raw counts。

本脚本不重新分割细胞；它只使用 labels_joint 将原始 2um counts
确定性地加总到细胞对象，并导出 H5AD/GeoJSON。
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import omicverse as ov
import scanpy as sc
from scipy import sparse
from shapely import wkt
from shapely.geometry import mapping

from batch_utils import crop_config_for_sample, load_batch_config, sample_by_index, sample_by_id, sample_root
from sp_v2_utils import write_h5ad_atomic


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run bin2cell and aggregate raw 2 um counts to Cellpose cells.")
    parser.add_argument("--sample-index", type=int, default=int(os.environ.get("LSB_JOBINDEX", "1")))
    parser.add_argument("--sample-id")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def ensure_spatial_scalefactors(adata, sample_id: str, bin_size_um: float = 2.0) -> None:
    spatial = adata.uns.get("spatial", {})
    if not isinstance(spatial, dict) or not spatial:
        adata.uns["spatial"] = {sample_id: {"scalefactors": {}}}
        spatial = adata.uns["spatial"]
    for library, payload in list(spatial.items()):
        if not isinstance(payload, dict):
            spatial[library] = {"scalefactors": {}}
            payload = spatial[library]
        scalefactors = payload.setdefault("scalefactors", {})
        if "spot_diameter_fullres" not in scalefactors:
            source_mpp = adata.uns.get("visium_hd_source_mpp", scalefactors.get("microns_per_pixel"))
            try:
                source_mpp = float(source_mpp)
            except (TypeError, ValueError):
                source_mpp = None
            scalefactors["spot_diameter_fullres"] = float(bin_size_um / source_mpp) if source_mpp and source_mpp > 0 else float(bin_size_um)
        scalefactors.setdefault("tissue_hires_scalef", 1.0)
        scalefactors.setdefault("tissue_lowres_scalef", 0.1)


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
    (export_dir / "complete_geojson_report.json").write_text(json.dumps({
        "geojson": str(geojson_path),
        "features": len(features),
        "skipped": len(skipped),
        "skipped_examples": skipped[:10],
        "backup": str(backup_path) if backup_path.exists() else None,
    }, indent=2, sort_keys=True) + "\n")
    return geojson_path


def aggregate_raw_counts(labels: ad.AnnData, cells: ad.AnnData, raw_2um_h5: Path) -> ad.AnnData:
    raw = sc.read_10x_h5(raw_2um_h5)
    raw.var_names_make_unique()
    missing_bins = labels.obs_names.difference(raw.obs_names)
    if len(missing_bins):
        raise ValueError(f"{len(missing_bins)} label bins are missing from raw 2um matrix")
    common_genes = labels.var_names.intersection(raw.var_names)
    if len(common_genes) == 0:
        raise ValueError("No overlapping genes between labeled bins and raw 2um matrix")

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
    raw_cells = ad.AnnData(
        X=raw_counts,
        obs=cells.obs.copy(),
        var=labels.var.copy(),
        uns=cells.uns.copy(),
        obsm={k: v.copy() for k, v in cells.obsm.items()},
    )
    raw_cells.uns["cellpose_raw_counts"] = {
        "raw_2um_h5": str(raw_2um_h5),
        "positive_labeled_bins": int(positive.sum()),
        "total_counts": int(raw_counts.sum()),
    }
    return raw_cells


def main() -> None:
    args = parse_args()
    row = sample_by_index(args.sample_index) if args.sample_id is None else sample_by_id(args.sample_id)
    sample_id = row["sample_id"]
    config = load_batch_config()
    crop_cfg = crop_config_for_sample(sample_id, config)
    out_dir = sample_root(sample_id) / "cellpose"
    out_dir.mkdir(parents=True, exist_ok=True)

    labels_input = out_dir / f"{sample_id}_002um_cellpose_he.h5ad"
    if not labels_input.exists():
        raise FileNotFoundError(f"Run 22_cellpose_he_each_sample.py first: {labels_input}")
    labels_output = out_dir / f"{sample_id}_002um_cellpose_labels.h5ad"
    cell_output = out_dir / f"{sample_id}_cell_level_cellpose.h5ad"
    raw_cell_output = out_dir / f"{sample_id}_cell_level_cellpose_raw_counts.h5ad"
    if raw_cell_output.exists() and not args.force:
        print(f"{raw_cell_output} exists; use --force to overwrite.")
        return

    labels = sc.read_h5ad(labels_input)
    ensure_spatial_scalefactors(labels, sample_id)
    if "labels_he_expanded" not in labels.obs:
        raise ValueError("Missing obs['labels_he_expanded']; rerun H&E Cellpose")
    labels.obs["labels_joint"] = labels.obs["labels_he_expanded"].astype(int)
    labels.obs["labels_joint_source"] = np.where(labels.obs["labels_joint"].astype(int) > 0, "he_expanded", "background")
    write_h5ad_atomic(labels, labels_output)

    cells = ov.space.bin2cell(
        labels,
        labels_key="labels_joint",
        spatial_keys=["spatial", "spatial_cropped_150_buffer"],
    )
    write_h5ad_atomic(cells, cell_output)

    raw_cells = aggregate_raw_counts(labels, cells, Path(row["matrix_002um"]))
    write_h5ad_atomic(raw_cells, raw_cell_output)

    export_dir = out_dir / "cellpose_spaceranger_raw_counts_output"
    export_data = ov.space.bin2cell(
        labels,
        labels_key="labels_joint",
        spatial_keys=["spatial"],
        add_geometry=True,
        geometry_spatial_key="spatial",
    )
    export_data = aggregate_raw_counts(labels, export_data, Path(row["matrix_002um"]))
    ov.space.write_visium_hd_cellseg(export_data, str(export_dir))
    geojson = write_complete_cell_geojson(export_data, export_dir)

    summary = {
        "sample_id": sample_id,
        "status": "pass",
        "labels_h5ad": str(labels_output),
        "cell_h5ad": str(cell_output),
        "raw_cell_h5ad": str(raw_cell_output),
        "spaceranger_export": str(export_dir),
        "complete_geojson": str(geojson),
        "cell_shape": list(raw_cells.shape),
        "positive_labeled_bins": int(raw_cells.uns["cellpose_raw_counts"]["positive_labeled_bins"]),
        "total_raw_counts": int(raw_cells.uns["cellpose_raw_counts"]["total_counts"]),
        "crop_config": crop_cfg,
    }
    (out_dir / "raw_counts_cellseg_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
