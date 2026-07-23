#!/usr/bin/env python3
"""Overlay v2 spatial domains on the Cellpose crop in the matching pixel frame."""
from __future__ import annotations

import json
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
    summary = []
    for row in read_manifest(MANIFEST):
        sample = row["sample_id"]
        h5ad_path = BATCH_RESULTS / "samples" / sample / "spatial_domains" / f"{sample}_spatial_domains.h5ad"
        image_path = SOURCE_WORK / sample / "cellpose_he" / "crop" / "cropped_source_image.png"
        if not h5ad_path.is_file() or not image_path.is_file():
            raise FileNotFoundError(f"{sample}: missing {h5ad_path} or {image_path}")
        adata = sc.read_h5ad(h5ad_path, backed="r")
        coords = np.asarray(adata.obsm["spatial"], dtype=float)
        domains = np.asarray(adata.obs["spatial_domain"].astype(str))
        image = np.asarray(Image.open(image_path).convert("RGB"))
        if coords[:, 0].max() > image.shape[1] or coords[:, 1].max() > image.shape[0]:
            raise ValueError(f"{sample}: spatial coordinates exceed cropped H&E image bounds")
        out_dir = BATCH_RESULTS / "samples" / sample / "spatial_domains" / "figures"
        out_dir.mkdir(parents=True, exist_ok=True)
        categories = sorted(set(domains), key=lambda x: int(x) if x.isdigit() else x)
        colors = plt.get_cmap("tab20")(np.linspace(0, 1, max(1, len(categories))))
        color_map = dict(zip(categories, colors))
        fig, ax = plt.subplots(figsize=(9, 9))
        ax.imshow(image)
        for category in categories:
            mask = domains == category
            ax.scatter(coords[mask, 0], coords[mask, 1], s=7, c=[color_map[category]], alpha=0.72, label=f"Domain {category}", linewidths=0, rasterized=True)
        ax.set_xlim(0, image.shape[1]); ax.set_ylim(image.shape[0], 0)
        ax.set_title(f"{sample}: H&E + spatial domains")
        ax.axis("off")
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=True, fontsize=8, title="Spatial domain")
        fig.tight_layout()
        output = out_dir / f"{sample}_spatial_domain_he_overlay.png"
        fig.savefig(output, dpi=220, bbox_inches="tight")
        plt.close(fig)
        summary.append({"sample_id": sample, "image": str(image_path), "h5ad": str(h5ad_path), "image_shape": list(image.shape[:2]), "n_objects": int(adata.n_obs), "domains": {x: int((domains == x).sum()) for x in categories}, "output": str(output), "coordinate_frame": "Cellpose cropped_source_image pixel frame"})
        adata.file.close()
    out = BATCH_RESULTS / "reports" / "step7_he_domain_overlay_summary.json"
    out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"status": "pass", "n_samples": len(summary), "summary": str(out)}, indent=2))


if __name__ == "__main__":
    main()
