#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import omicverse as ov
import scanpy as sc
import tifffile
from PIL import Image

from batch_utils import sample_by_id, sample_root
from sp_v2_utils import write_h5ad_atomic

PROJECT_HOME = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))


def helpers():
    path = PROJECT_HOME / "scripts" / "23_bin2cell_raw_counts_each_sample.py"
    spec = importlib.util.spec_from_file_location("bin2cell_helpers", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def overlay(adata, image_path, output, sample, tag):
    image = np.asarray(Image.open(image_path).convert("RGB"))
    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    fig, ax = plt.subplots(figsize=(9, 9)); ax.imshow(image)
    gex = adata.obs["labels_gex"].astype(int).to_numpy() > 0
    joint = adata.obs["labels_joint_source"].astype(str).to_numpy()
    ax.scatter(coords[gex, 0], coords[gex, 1], s=1.6, c="#d946ef", alpha=.18,
               linewidths=0, rasterized=True, label="GEX labeled 2um bins")
    sec = joint == "secondary"
    ax.scatter(coords[sec, 0], coords[sec, 1], s=2.0, c="#22c55e", alpha=.35,
               linewidths=0, rasterized=True, label="GEX-salvaged bins")
    ax.set_xlim(0, image.shape[1]); ax.set_ylim(image.shape[0], 0); ax.axis("off")
    ax.set_title(f"{sample} {tag}: GEX labels and salvage")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=8)
    fig.tight_layout(); fig.savefig(output, dpi=220, bbox_inches="tight"); plt.close(fig)


def main():
    sample, tag = os.environ["SP_GEX_SAMPLE"], os.environ["SP_GEX_TAG"]
    variant = os.environ.get("SP_GEX_VARIANT", "baseline")
    # 允许在 GPU 上做小规模参数扫描，避免每次修改脚本后重新上传。
    mpp = float(os.environ.get("SP_GEX_MPP", "0.3"))
    sigma = float(os.environ.get("SP_GEX_SIGMA", "5"))
    norm_lo = float(os.environ.get("SP_GEX_NORM_LO", "1"))
    norm_hi = float(os.environ.get("SP_GEX_NORM_HI", "99.5"))
    mask_threshold = float(os.environ.get("SP_GEX_MASK_THRESHOLD", "0.01"))
    flow_threshold = float(os.environ.get("SP_GEX_FLOW_THRESHOLD", "0.4"))
    diameter_raw = os.environ.get("SP_GEX_DIAMETER", "")
    diameter = float(diameter_raw) if diameter_raw else None
    root = sample_root(sample)
    source = root / "cellpose_sweep" / tag
    out = root / "gex_compare" / f"{tag}__{variant}"
    out.mkdir(parents=True, exist_ok=True)
    labels = sc.read_h5ad(source / f"{sample}_002um_cellpose_he.h5ad")
    h = helpers(); h.ensure_spatial_scalefactors(labels, sample)
    labels.obs["n_counts_adjusted"] = np.asarray(labels.X.sum(axis=1)).ravel().astype(float)
    # 对 GEX 图像做 log1p 和百分位强度归一化；原始 uint8 counts 图像
    # 在本批样本上过于稀疏，Cellpose 会返回空标签。
    from omicverse.external.bin2cell import grid_image, cellseg, insert_labels
    raw_gex_image = out / "gex_counts_log1p_raw.tiff"
    normalized_gex_image = out / "gex_counts_log1p_norm.tiff"
    grid_image(labels, "n_counts_adjusted", log1p=True, mpp=mpp, sigma=sigma,
               save_path=str(raw_gex_image))
    image = tifffile.imread(raw_gex_image).astype(np.float32)
    positive = image[image > 0]
    if positive.size == 0:
        raise ValueError("log1p GEX image has no positive pixels")
    lo, hi = np.percentile(positive, [norm_lo, norm_hi])
    image = np.clip((image - lo) / max(hi - lo, 1e-6), 0, 1) * 255
    tifffile.imwrite(normalized_gex_image, image.astype(np.uint8))
    labels_npz = out / "gex_counts_log1p_norm.npz"
    use_gpu = os.environ.get("SP_GEX_USE_GPU", "0").lower() not in {"0", "false", "no"}
    cellseg(str(normalized_gex_image), str(labels_npz),
            stardist_model="2D_versatile_he", backend="cellpose",
            gpu=use_gpu, block_size=2048, min_overlap=256, context=64,
            diameter=diameter, mask_threshold=mask_threshold,
            flow_threshold=flow_threshold, show_progress=True)
    # GEX 网格图像的坐标来自 array_row/array_col，不能用 H&E 的 spatial
    # 坐标回写；basis=array 会按 bin2cell 的 array 坐标规则匹配 GEX 像素。
    insert_labels(labels, labels_npz_path=str(labels_npz), basis="array", mpp=mpp,
                  labels_key="labels_gex")
    ov.space.salvage_secondary_labels(labels, primary_label="labels_he_expanded",
                                      secondary_label="labels_gex", labels_key="labels_joint")
    labels.obs["labels_joint_source"] = labels.obs["labels_joint_source"].astype(str)
    write_h5ad_atomic(labels, out / f"{sample}_{tag}_gex_labels.h5ad")
    gex_cells = ov.space.bin2cell(labels, labels_key="labels_gex",
                                  spatial_keys=["spatial", "spatial_cropped_150_buffer"])
    joint_cells = ov.space.bin2cell(labels, labels_key="labels_joint",
                                    spatial_keys=["spatial", "spatial_cropped_150_buffer"])
    raw_matrix = Path(sample_by_id(sample)["matrix_002um"])
    gex_raw = h.aggregate_raw_counts(labels, gex_cells, raw_matrix)
    joint_raw = h.aggregate_raw_counts(labels, joint_cells, raw_matrix)
    write_h5ad_atomic(gex_raw, out / f"{sample}_{tag}_gex_cell_raw_counts.h5ad")
    write_h5ad_atomic(joint_raw, out / f"{sample}_{tag}_joint_cell_raw_counts.h5ad")
    image = PROJECT_HOME / "work" / "batch_v2_corrected" / sample / "cellpose_sweep" / tag / "crop" / "cropped_source_image.png"
    fig = out / f"{sample}_{tag}_gex_salvage_overlay.png"
    if image.exists(): overlay(labels, image, fig, sample, tag)
    primary = labels.obs["labels_he_expanded"].astype(int) > 0
    secondary = labels.obs["labels_joint_source"] == "secondary"
    summary = {
        "sample_id": sample, "tag": tag, "variant": variant, "status": "pass", "mpp": mpp, "sigma": sigma,
        "he_labeled_bins": int(primary.sum()), "gex_labeled_bins": int((labels.obs["labels_gex"].astype(int) > 0).sum()),
        "salvaged_bins": int(secondary.sum()), "total_bins": int(labels.n_obs),
        "gex_cells": int(gex_raw.n_obs), "joint_cells": int(joint_raw.n_obs),
        "gex_coverage_percent": float(100 * (labels.obs["labels_gex"].astype(int) > 0).sum() / labels.n_obs),
        "joint_coverage_percent": float(100 * (labels.obs["labels_joint"].astype(int) > 0).sum() / labels.n_obs),
        "overlay": str(fig) if fig.exists() else None,
    }
    (out / "gex_segmentation_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__": main()
