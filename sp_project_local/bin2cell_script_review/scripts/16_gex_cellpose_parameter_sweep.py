#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import omicverse as ov
import scanpy as sc


def ensure_inputs(adata) -> None:
    if "labels_he_expanded" not in adata.obs:
        raise ValueError("Missing labels_he_expanded")
    if "total_counts" not in adata.obs:
        adata.obs["total_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    spatial = adata.uns.get("spatial", {})
    for payload in spatial.values():
        scalefactors = payload.setdefault("scalefactors", {})
        if "spot_diameter_fullres" not in scalefactors:
            mpp = float(adata.uns.get("visium_hd_source_mpp", scalefactors.get("microns_per_pixel", 1.0)))
            scalefactors["spot_diameter_fullres"] = 2.0 / mpp


def main() -> None:
    parser = argparse.ArgumentParser(description="Small GEX Cellpose parameter sweep for Visium HD pilot crop.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mpp", type=float, default=0.3239732935)
    parser.add_argument("--buffer", type=int, default=150)
    parser.add_argument("--gpu", action=argparse.BooleanOptionalAction, default=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    grid = [
        {"name": "sigma2_prob-2e-2_nms10", "sigma": 2.0, "prob_thresh": 0.02, "nms_thresh": 0.10},
        {"name": "sigma2_prob5e-3_nms10", "sigma": 2.0, "prob_thresh": 0.005, "nms_thresh": 0.10},
        {"name": "sigma5_prob5e-3_nms10", "sigma": 5.0, "prob_thresh": 0.005, "nms_thresh": 0.10},
        {"name": "sigma8_prob5e-3_nms05", "sigma": 8.0, "prob_thresh": 0.005, "nms_thresh": 0.05},
    ]
    rows = []
    for params in grid:
        adata = sc.read_h5ad(args.input)
        ensure_inputs(adata)
        out_tiff = args.output_dir / f"{params['name']}.tiff"
        try:
            ov.space.visium_10x_hd_cellpose_gex(
                adata,
                obs_key="total_counts",
                log1p=False,
                mpp=args.mpp,
                sigma=params["sigma"],
                gex_save_path=str(out_tiff),
                prob_thresh=params["prob_thresh"],
                nms_thresh=params["nms_thresh"],
                gpu=args.gpu,
                buffer=args.buffer,
            )
            labels = adata.obs.get("labels_gex")
            label_values = labels.astype(int).to_numpy() if labels is not None else np.zeros(adata.n_obs, dtype=int)
            nonzero = int((label_values > 0).sum())
            unique = int(len(set(label_values[label_values > 0])))
            status = "ok"
            error = None
        except Exception as exc:
            nonzero = 0
            unique = 0
            status = "fail"
            error = repr(exc)
        rows.append({
            **params,
            "status": status,
            "labels_gex_nonzero_bins": nonzero,
            "labels_gex_unique_objects": unique,
            "output_tiff": str(out_tiff),
            "error": error,
        })
        print(json.dumps(rows[-1], ensure_ascii=False))

    report = {
        "input": str(args.input),
        "mpp": args.mpp,
        "buffer": args.buffer,
        "gpu": args.gpu,
        "results": rows,
    }
    (args.output_dir / "gex_cellpose_parameter_sweep.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
