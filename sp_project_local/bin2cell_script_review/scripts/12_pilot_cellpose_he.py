#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import omicverse as ov
from sp_v2_utils import crop_visium_hd_um, read_visium_hd_compat


def main() -> None:
    parser = argparse.ArgumentParser(description="Pilot Cellpose H&E segmentation on a cropped 2 um Visium HD region.")
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--source-image", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mpp", type=float, default=0.3239732935)
    parser.add_argument("--buffer", type=int, default=150)
    parser.add_argument("--x-range", type=float, nargs=2, required=True)
    parser.add_argument("--y-range", type=float, nargs=2, required=True)
    parser.add_argument("--gpu", action=argparse.BooleanOptionalAction, default=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    he_image = args.output_dir / "he_cellpose.tiff"
    adata = read_visium_hd_compat(args.bin_dir, args.source_image, args.output_dir / "compat")
    ov.pp.filter_genes(adata, min_cells=3)
    ov.pp.filter_cells(adata, min_counts=1)
    adata = crop_visium_hd_um(adata, args.x_range, args.y_range, args.output_dir)

    ov.space.visium_10x_hd_cellpose_he(
        adata,
        mpp=args.mpp,
        he_save_path=str(he_image),
        prob_thresh=0,
        flow_threshold=0.4,
        gpu=args.gpu,
        buffer=args.buffer,
        backend="pil",
    )
    ov.space.visium_10x_hd_cellpose_expand(
        adata,
        labels_key="labels_he",
        expanded_labels_key="labels_he_expanded",
        max_bin_distance=4,
    )
    output = args.output_dir / "SC000895-R4_002um_cellpose_he.h5ad"
    adata.write_h5ad(output, compression="lzf")
    summary = {
        "bin_dir": str(args.bin_dir),
        "source_image": str(args.source_image),
        "input_shape_after_crop": list(adata.shape),
        "mpp": args.mpp,
        "buffer": args.buffer,
        "gpu": args.gpu,
        "label_columns": [key for key in ["labels_he", "labels_he_expanded"] if key in adata.obs],
        "output": str(output),
    }
    (args.output_dir / "cellpose_he_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
