#!/usr/bin/env python3
"""中文说明：对单个 Visium HD 样本执行 H&E Cellpose 和标签扩展。

输入是 2um bins 与 H&E 图像；输出是 labels_he/labels_he_expanded，
后续由 23_bin2cell_raw_counts_each_sample.py 聚合为细胞矩阵。
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import omicverse as ov

from batch_utils import BATCH_WORK, crop_config_for_sample, load_batch_config, sample_by_index, sample_by_id, sample_root
from sp_v2_utils import crop_visium_hd_um, read_visium_hd_compat, write_h5ad_atomic


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run H&E Cellpose on one cropped 2 um Visium HD sample.")
    parser.add_argument("--sample-index", type=int, default=int(os.environ.get("LSB_JOBINDEX", "1")))
    parser.add_argument("--sample-id")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def total_counts(adata) -> float:
    return float(np.asarray(adata.X.sum()).ravel()[0])


def crop_or_auto(adata, crop_cfg: dict, work_dir: Path):
    x_range = crop_cfg.get("x_range_um", [2000, 3000])
    y_range = crop_cfg.get("y_range_um", [2000, 3000])
    cropped = crop_visium_hd_um(adata, x_range, y_range, work_dir / "crop")
    if cropped.n_obs > 0 and total_counts(cropped) > 0:
        cropped.uns["cellpose_crop_mode"] = "configured"
        return cropped, {"x_range_um": x_range, "y_range_um": y_range, "mode": "configured"}

    coords = np.asarray(adata.obsm["spatial_um"])
    counts = np.asarray(adata.X.sum(axis=1)).ravel()
    positive = counts > 0
    if not positive.any():
        raise ValueError("No positive-count 2um bins are available for automatic Cellpose crop selection")
    width = float(x_range[1] - x_range[0])
    height = float(y_range[1] - y_range[0])
    center_x = float(np.median(coords[positive, 0]))
    center_y = float(np.median(coords[positive, 1]))
    auto_x = [max(0.0, center_x - width / 2), center_x + width / 2]
    auto_y = [max(0.0, center_y - height / 2), center_y + height / 2]
    cropped = crop_visium_hd_um(adata, auto_x, auto_y, work_dir / "crop_auto")
    if cropped.n_obs == 0 or total_counts(cropped) == 0:
        raise ValueError(f"Automatic crop still has no positive bins: x={auto_x}, y={auto_y}")
    cropped.uns["cellpose_crop_mode"] = "auto_positive_median"
    return cropped, {"x_range_um": auto_x, "y_range_um": auto_y, "mode": "auto_positive_median"}


def main() -> None:
    args = parse_args()
    row = sample_by_index(args.sample_index) if args.sample_id is None else sample_by_id(args.sample_id)
    sample_id = row["sample_id"]
    config = load_batch_config()
    crop_cfg = crop_config_for_sample(sample_id, config)

    output_subdir = os.environ.get("SP_CELLPOSE_OUTPUT_SUBDIR", "cellpose")
    out_dir = sample_root(sample_id) / output_subdir
    work_dir = BATCH_WORK / sample_id / output_subdir
    out_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    output_h5ad = out_dir / f"{sample_id}_002um_cellpose_he.h5ad"
    if output_h5ad.exists() and not args.force:
        print(f"{output_h5ad} exists; use --force to overwrite.")
        return

    ov.settings.cpu_gpu_mixed_init()
    bin_dir = Path(row["matrix_002um"]).parent
    source_image = Path(row["outs_dir"]) / "spatial" / "tissue_hires_image.png"
    adata = read_visium_hd_compat(bin_dir, source_image, work_dir)
    ov.pp.filter_genes(adata, min_cells=3)
    ov.pp.filter_cells(adata, min_counts=1)
    cropped, effective_crop = crop_or_auto(adata, crop_cfg, work_dir)
    cropped.var_names_make_unique()
    he_image = out_dir / "he_cellpose.tiff"
    configured_mpp = float(crop_cfg.get("mpp", 0.3239732935))
    source_mpp = float(cropped.uns.get("visium_hd_source_mpp", configured_mpp))
    used_mpp = source_mpp if crop_cfg.get("mpp_mode") == "source" else configured_mpp
    use_gpu = os.environ.get("SP_CELLPOSE_USE_GPU", "1").lower() not in {"0", "false", "no"}
    prob_thresh = float(os.environ.get("SP_CELLPOSE_PROB_THRESH", "0"))
    flow_threshold = float(os.environ.get("SP_CELLPOSE_FLOW_THRESHOLD", "0.4"))
    ov.space.visium_10x_hd_cellpose_he(
        cropped,
        mpp=used_mpp,
        he_save_path=str(he_image),
        prob_thresh=prob_thresh,
        flow_threshold=flow_threshold,
        gpu=use_gpu,
        buffer=int(crop_cfg.get("buffer", 150)),
        backend="pil",
    )
    ov.space.visium_10x_hd_cellpose_expand(
        cropped,
        labels_key="labels_he",
        expanded_labels_key="labels_he_expanded",
        max_bin_distance=int(crop_cfg.get("expand_max_bin_distance", 4)),
    )
    write_h5ad_atomic(cropped, output_h5ad)
    summary = {
        "sample_id": sample_id,
        "status": "pass",
        "h5ad": str(output_h5ad),
        "n_bins": int(cropped.n_obs),
        "n_genes": int(cropped.n_vars),
        "labels_he_nonzero_bins": int((cropped.obs.get("labels_he", 0).astype(int) > 0).sum()),
        "labels_he_expanded_nonzero_bins": int((cropped.obs.get("labels_he_expanded", 0).astype(int) > 0).sum()),
        "crop_config": crop_cfg,
        "effective_crop": effective_crop,
        "source_mpp": source_mpp,
        "cellpose_mpp_used": used_mpp,
        "cellpose_gpu_used": use_gpu,
        "cellpose_prob_thresh": prob_thresh,
        "cellpose_flow_threshold": flow_threshold,
        "output_subdir": output_subdir,
        "svg_step": "skipped_for_he_cellpose_batch_v1",
    }
    (out_dir / "cellpose_he_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
