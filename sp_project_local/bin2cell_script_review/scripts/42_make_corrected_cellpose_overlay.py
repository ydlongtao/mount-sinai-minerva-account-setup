#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from PIL import Image

from batch_utils import BATCH_RESULTS, MANIFEST, read_manifest


PROJECT_HOME = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
SOURCE_WORK = Path(os.environ.get("SP_PROJECT_INPUT_BATCH_WORK", str(PROJECT_HOME / "work" / "batch")))


def main() -> None:
    for row in read_manifest(MANIFEST):
        sample = row["sample_id"]
        cellpose = BATCH_RESULTS / "samples" / sample / "cellpose"
        labels_path = cellpose / f"{sample}_002um_cellpose_he.h5ad"
        cells_path = cellpose / f"{sample}_cell_level_cellpose_raw_counts.h5ad"
        image_path = SOURCE_WORK / sample / "cellpose_he" / "crop" / "cropped_source_image.png"
        if not labels_path.is_file() or not image_path.is_file():
            continue
        labels = sc.read_h5ad(labels_path, backed="r")
        image = np.asarray(Image.open(image_path).convert("RGB"))
        coords = np.asarray(labels.obsm["spatial"], dtype=float)
        expanded = np.asarray(labels.obs["labels_he_expanded"], dtype=int) > 0
        out = BATCH_RESULTS / "samples" / sample / "cellpose" / f"{sample}_corrected_he_coverage.png"
        fig, ax = plt.subplots(figsize=(9, 9))
        ax.imshow(image)
        ax.scatter(coords[expanded, 0], coords[expanded, 1], s=1.2, c="#ffb000", alpha=0.18, linewidths=0, rasterized=True, label="expanded labeled 2um bins")
        if cells_path.is_file():
            cells = sc.read_h5ad(cells_path, backed="r")
            cc = np.asarray(cells.obsm["spatial"], dtype=float)
            ax.scatter(cc[:, 0], cc[:, 1], s=9, c="#00d4ff", alpha=0.85, linewidths=0, rasterized=True, label="bin2cell centroids")
            cells.file.close()
        ax.set_xlim(0, image.shape[1]); ax.set_ylim(image.shape[0], 0); ax.axis("off")
        ax.set_title(f"{sample}: corrected source-mpp H&E coverage")
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=8)
        fig.tight_layout(); fig.savefig(out, dpi=220, bbox_inches="tight"); plt.close(fig)
        labels.file.close()
        print(out)


if __name__ == "__main__": main()
