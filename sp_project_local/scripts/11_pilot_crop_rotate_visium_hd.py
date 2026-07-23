#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import omicverse as ov
from sp_v2_utils import crop_visium_hd_um, read_visium_hd_compat


def main() -> None:
    parser = argparse.ArgumentParser(description="Pilot OmicVerse image-aware crop and rotate on Visium HD.")
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--source-image", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--x-range", type=float, nargs=2, required=True)
    parser.add_argument("--y-range", type=float, nargs=2, required=True)
    parser.add_argument("--angle", type=float, default=0)
    parser.add_argument("--map-after-rotate", action="store_true")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    adata = read_visium_hd_compat(args.bin_dir, args.source_image, args.output_dir / "compat")
    cropped = crop_visium_hd_um(adata, args.x_range, args.y_range, args.output_dir)
    library_id = list(cropped.uns["spatial"].keys())[0]
    cropped_path = args.output_dir / "SC000895-R4_008um_crop.h5ad"
    cropped.write_h5ad(cropped_path, compression="lzf")

    rotated = cropped
    if args.angle:
        rotated = ov.space.rotate_space_visium(
            cropped, angle=args.angle, library_id=library_id,
            interpolation_order=1,
        )
    map_skipped_reason = None
    if args.map_after_rotate and "Anno_manual" in rotated.obs:
        ov.space.map_spatial_auto(rotated, method="phase")
    elif args.map_after_rotate:
        map_skipped_reason = "Anno_manual is absent; map_spatial_auto was skipped"
    if "spatial" not in rotated.obsm:
        raise RuntimeError("OmicVerse rotation did not return obsm['spatial']")
    rotated_path = args.output_dir / "SC000895-R4_008um_crop_rotated.h5ad"
    rotated.write_h5ad(rotated_path, compression="lzf")

    summary = {
        "bin_dir": str(args.bin_dir),
        "source_image": str(args.source_image),
        "x_range_um": list(args.x_range),
        "y_range_um": list(args.y_range),
        "angle": args.angle,
        "map_after_rotate": args.map_after_rotate,
        "map_skipped_reason": map_skipped_reason,
        "input_shape": list(adata.shape),
        "cropped_shape": list(cropped.shape),
        "rotated_shape": list(rotated.shape),
        "cropped_output": str(cropped_path),
        "rotated_output": str(rotated_path),
    }
    (args.output_dir / "crop_rotate_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
