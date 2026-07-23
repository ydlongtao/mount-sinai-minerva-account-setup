#!/usr/bin/env python3
from pathlib import Path

import openslide


SVS_PATH = Path("/Volumes/Extreme SSD/关键work/Data/XeniumofGC/BS06-9313-8_Tumor__2.svs")
OUT_DIR = Path("/Volumes/Extreme SSD/关键work/Data/XeniumofGC/svs_previews")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main() -> None:
    slide = openslide.OpenSlide(str(SVS_PATH))

    print(f"File: {SVS_PATH}")
    print(f"Vendor: {slide.properties.get(openslide.PROPERTY_NAME_VENDOR, 'unknown')}")
    print(f"Level count: {slide.level_count}")
    print("Levels:")
    for idx, (dims, downsample) in enumerate(zip(slide.level_dimensions, slide.level_downsamples)):
        print(f"  level {idx}: {dims[0]} x {dims[1]}, downsample={downsample:.2f}")

    # Choose a small pyramid level for a quick whole-slide overview.
    target_max_side = 2500
    level = min(
        range(slide.level_count),
        key=lambda i: abs(max(slide.level_dimensions[i]) - target_max_side),
    )
    lowres_size = slide.level_dimensions[level]
    overview = slide.read_region((0, 0), level, lowres_size).convert("RGB")
    overview_path = OUT_DIR / f"{SVS_PATH.stem}_overview_level{level}_{lowres_size[0]}x{lowres_size[1]}.jpg"
    overview.save(overview_path, quality=90)
    print(f"Saved overview: {overview_path}")

    # Also export a low-resolution center tile from level 0 coordinates.
    w0, h0 = slide.level_dimensions[0]
    tile_size = 2048
    x0 = max((w0 - tile_size) // 2, 0)
    y0 = max((h0 - tile_size) // 2, 0)
    tile = slide.read_region((x0, y0), 0, (tile_size, tile_size)).convert("RGB")
    tile = tile.resize((1024, 1024))
    tile_path = OUT_DIR / f"{SVS_PATH.stem}_center_tile_level0_1024.jpg"
    tile.save(tile_path, quality=90)
    print(f"Saved center tile: {tile_path}")

    slide.close()


if __name__ == "__main__":
    main()
