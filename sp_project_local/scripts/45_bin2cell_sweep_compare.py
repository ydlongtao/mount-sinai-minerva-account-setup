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
import scanpy as sc
from PIL import Image

from batch_utils import sample_by_id, sample_root
from sp_v2_utils import write_h5ad_atomic

PROJECT_HOME = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))


def load_bin2cell_helpers():
    path = PROJECT_HOME / "scripts" / "23_bin2cell_raw_counts_each_sample.py"
    spec = importlib.util.spec_from_file_location("bin2cell_helpers", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def make_overlay(labels, cells, image_path: Path, output: Path, title: str) -> None:
    image = np.asarray(Image.open(image_path).convert("RGB"))
    coords = np.asarray(labels.obsm["spatial"], dtype=float)
    expanded = labels.obs["labels_he_expanded"].astype(int).to_numpy() > 0
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.imshow(image)
    ax.scatter(coords[expanded, 0], coords[expanded, 1], s=1.2, c="#ffb000", alpha=0.18,
               linewidths=0, rasterized=True, label="expanded labeled 2um bins")
    if cells is not None and cells.n_obs:
        cc = np.asarray(cells.obsm["spatial"], dtype=float)
        ax.scatter(cc[:, 0], cc[:, 1], s=9, c="#00d4ff", alpha=0.85,
                   linewidths=0, rasterized=True, label="bin2cell centroids")
    ax.set_xlim(0, image.shape[1]); ax.set_ylim(image.shape[0], 0); ax.axis("off")
    ax.set_title(title)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=8)
    fig.tight_layout(); fig.savefig(output, dpi=220, bbox_inches="tight"); plt.close(fig)


def main() -> None:
    sample = os.environ["SP_SWEEP_SAMPLE"]
    tag = os.environ["SP_SWEEP_TAG"]
    root = sample_root(sample)
    sweep_dir = root / "cellpose_sweep" / tag
    labels_path = sweep_dir / f"{sample}_002um_cellpose_he.h5ad"
    raw_matrix = Path(sample_by_id(sample)["matrix_002um"])
    out_dir = sweep_dir / "bin2cell"
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True); fig_dir.mkdir(parents=True, exist_ok=True)
    labels_output = out_dir / f"{sample}_002um_cellpose_labels.h5ad"
    cell_output = out_dir / f"{sample}_cell_level_cellpose.h5ad"
    raw_output = out_dir / f"{sample}_cell_level_cellpose_raw_counts.h5ad"
    if not labels_path.exists():
        raise FileNotFoundError(labels_path)

    helpers = load_bin2cell_helpers()
    labels = sc.read_h5ad(labels_path)
    helpers.ensure_spatial_scalefactors(labels, sample)
    labels.obs["labels_joint"] = labels.obs["labels_he_expanded"].astype(int)
    labels.obs["labels_joint_source"] = np.where(labels.obs["labels_joint"] > 0, "he_expanded", "background")
    write_h5ad_atomic(labels, labels_output)
    cells = __import__("omicverse").space.bin2cell(labels, labels_key="labels_joint",
        spatial_keys=["spatial", "spatial_cropped_150_buffer"])
    write_h5ad_atomic(cells, cell_output)
    raw_cells = helpers.aggregate_raw_counts(labels, cells, raw_matrix)
    write_h5ad_atomic(raw_cells, raw_output)

    image_path = PROJECT_HOME / "work" / "batch_v2_corrected" / sample / "cellpose_sweep" / tag / "crop" / "cropped_source_image.png"
    overlay = fig_dir / f"{sample}_{tag}_bin2cell_overlay.png"
    if image_path.exists():
        make_overlay(labels, cells, image_path, overlay, f"{sample} {tag}: H&E + bin2cell")
    summary = {
        "sample_id": sample, "tag": tag, "status": "pass",
        "labels_h5ad": str(labels_path), "cell_h5ad": str(cell_output), "raw_cell_h5ad": str(raw_output),
        "n_bins": int(labels.n_obs), "n_genes": int(labels.n_vars), "cells": int(raw_cells.n_obs),
        "labeled_bins": int((labels.obs["labels_joint"] > 0).sum()),
        "coverage_percent": float(100 * (labels.obs["labels_joint"] > 0).sum() / labels.n_obs),
        "median_total_counts": float(np.median(np.asarray(raw_cells.X.sum(axis=1)).ravel())) if raw_cells.n_obs else 0.0,
        "overlay": str(overlay) if overlay.exists() else None,
    }
    (out_dir / "bin2cell_sweep_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
