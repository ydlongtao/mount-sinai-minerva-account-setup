#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from shapely import wkt
from shapely.geometry import mapping
import omicverse as ov
import scanpy as sc



def ensure_spatial_scalefactors(adata, bin_size_um: float = 2.0) -> None:
    spatial = adata.uns.get("spatial", {})
    if not isinstance(spatial, dict) or not spatial:
        adata.uns["spatial"] = {"SC000895-R4": {"scalefactors": {}}}
        spatial = adata.uns["spatial"]
    for library, payload in spatial.items():
        if not isinstance(payload, dict):
            spatial[library] = {"scalefactors": {}}
            payload = spatial[library]
        scalefactors = payload.setdefault("scalefactors", {})
        if "spot_diameter_fullres" not in scalefactors:
            source_mpp = adata.uns.get("visium_hd_source_mpp", scalefactors.get("microns_per_pixel", None))
            try:
                source_mpp = float(source_mpp)
            except (TypeError, ValueError):
                source_mpp = None
            if source_mpp and source_mpp > 0:
                scalefactors["spot_diameter_fullres"] = float(bin_size_um / source_mpp)
            else:
                scalefactors["spot_diameter_fullres"] = float(bin_size_um)
        scalefactors.setdefault("tissue_hires_scalef", 1.0)
        scalefactors.setdefault("tissue_lowres_scalef", 0.1)



def write_complete_cell_geojson(cdata, export_dir: Path) -> Path:
    geojson_path = export_dir / "graphclust_annotated_cell_segmentations.geojson"
    backup_path = export_dir / "graphclust_annotated_cell_segmentations.omicsverse_simple.geojson"
    if geojson_path.exists() and not backup_path.exists():
        geojson_path.replace(backup_path)

    features = []
    skipped = []
    for i, (idx, row) in enumerate(cdata.obs.iterrows(), start=1):
        geom_wkt = row.get("geometry", "")
        try:
            geom = wkt.loads(str(geom_wkt))
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

    geojson = {
        "type": "FeatureCollection",
        "features": features,
    }
    geojson_path.write_text(json.dumps(geojson, separators=(",", ":")) + "\n")
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
    parser = argparse.ArgumentParser(description="Add GEX Cellpose labels and aggregate Visium HD bins to cells.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mpp", type=float, default=0.3239732935)
    parser.add_argument("--buffer", type=int, default=150)
    parser.add_argument("--sigma", type=float, default=5)
    parser.add_argument("--gpu", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--export-spaceranger", action="store_true")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    gex_image = args.output_dir / "gex_cellpose.tiff"
    bin_output = args.output_dir / "SC000895-R4_002um_cellpose_labels.h5ad"
    if bin_output.exists():
        adata = sc.read_h5ad(bin_output)
        print(f"Resuming from existing bin-level labels: {bin_output}")
    else:
        adata = sc.read_h5ad(args.input)
        ensure_spatial_scalefactors(adata)
        if "labels_he_expanded" not in adata.obs:
            raise ValueError("Missing required primary label: labels_he_expanded")
        if "total_counts" not in adata.obs:
            adata.obs["total_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
        ov.space.visium_10x_hd_cellpose_gex(
            adata, obs_key="total_counts", log1p=False, mpp=args.mpp,
            sigma=args.sigma, gex_save_path=str(gex_image), prob_thresh=0.01,
            nms_thresh=0.1, gpu=args.gpu, buffer=args.buffer,
        )
        ov.space.salvage_secondary_labels(
            adata, primary_label="labels_he_expanded",
            secondary_label="labels_gex", labels_key="labels_joint",
        )
        adata.write_h5ad(bin_output, compression="lzf")
    ensure_spatial_scalefactors(adata)
    adata.write_h5ad(bin_output, compression="lzf")
    cdata = ov.space.bin2cell(
        adata,
        labels_key="labels_joint",
        spatial_keys=["spatial", "spatial_cropped_150_buffer"],
    )
    cell_output = args.output_dir / "SC000895-R4_cell_level_cellpose.h5ad"
    cdata.write_h5ad(cell_output, compression="lzf")

    export_dir = None
    if args.export_spaceranger:
        export_data = ov.space.bin2cell(
            adata,
            labels_key="labels_joint",
            spatial_keys=["spatial"],
            add_geometry=True,
            geometry_spatial_key="spatial",
        )
        export_dir = args.output_dir / "cellpose_spaceranger_output"
        ov.space.write_visium_hd_cellseg(export_data, str(export_dir))
        full_geojson = write_complete_cell_geojson(export_data, export_dir)
        print(f"[VisiumHD][OK] Complete GeoJSON written: {full_geojson}")

    summary = {
        "input": str(args.input),
        "bin_label_shape": list(adata.shape),
        "cell_shape": list(cdata.shape),
        "labels": [key for key in ["labels_he", "labels_he_expanded", "labels_gex", "labels_joint"] if key in adata.obs],
        "bin_output": str(bin_output),
        "cell_output": str(cell_output),
        "spaceranger_export": str(export_dir) if export_dir else None,
        "complete_geojson_report": str(export_dir / "complete_geojson_report.json") if export_dir else None,
        "gpu": args.gpu,
    }
    (args.output_dir / "cellpose_bin2cell_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
