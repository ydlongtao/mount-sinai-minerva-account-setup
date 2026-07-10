#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
from pathlib import Path
from typing import Any

import numpy as np


PROJECT_HOME = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
WORK_V2 = Path(os.environ.get("SP_PROJECT_WORK_V2", str(PROJECT_HOME / "work_v2")))
RESULTS_V2 = PROJECT_HOME / "results" / "pilot_v2"
BARCODE_PATTERN = re.compile(r"^s_(\d+)um_(\d+)_(\d+)-\d+$")


def load_config(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write_json_atomic(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def write_h5ad_atomic(adata, path: Path, compression: str = "lzf") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.stem}.{os.getpid()}.tmp.h5ad")
    adata.write_h5ad(tmp, compression=compression)
    os.replace(tmp, path)


def file_fingerprint(path: Path) -> dict[str, Any]:
    stat = path.stat()
    payload = f"{path.resolve()}:{stat.st_size}:{stat.st_mtime_ns}"
    return {
        "path": str(path),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "fingerprint": hashlib.sha256(payload.encode()).hexdigest(),
    }


def read_manifest_sample(path: Path, sample_id: str) -> dict[str, str]:
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["sample_id"] == sample_id:
                return row
    raise KeyError(f"Sample {sample_id!r} was not found in {path}")


def attach_barcode_spatial(adata) -> None:
    parsed = [BARCODE_PATTERN.match(str(barcode)) for barcode in adata.obs_names]
    if not parsed or any(match is None for match in parsed):
        raise ValueError("Visium HD barcodes could not be converted to spatial coordinates")
    sizes = {int(match.group(1)) for match in parsed if match is not None}
    if len(sizes) != 1:
        raise ValueError(f"Mixed Visium HD bin sizes: {sorted(sizes)}")
    size = sizes.pop()
    coords = np.empty((adata.n_obs, 2), dtype=np.float32)
    coords[:, 0] = [int(match.group(3)) * size + size / 2 for match in parsed if match is not None]
    coords[:, 1] = [int(match.group(2)) * size + size / 2 for match in parsed if match is not None]
    adata.obsm["spatial"] = coords
    adata.uns["spatial_coordinate_source"] = "visium_hd_bin_barcode"
    adata.uns["spatial_coordinate_unit"] = "micrometer"


def read_visium_hd_compat(bin_dir: Path, source_image: Path, work_dir: Path):
    """Load Visium HD output through OmicVerse's bin2cell-compatible reader.

    Some Cell Ranger HD outputs contain barcode-encoded coordinates but omit
    the legacy ``tissue_positions_list.csv`` and scale-factor files.  This
    adapter materializes those small metadata files in scratch, leaving the
    raw directory untouched.
    """
    import h5py
    import pandas as pd
    from PIL import Image
    import omicverse as ov

    bin_dir = Path(bin_dir)
    source_image = Path(source_image)
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    h5_path = bin_dir / "filtered_feature_bc_matrix.h5"
    if not h5_path.exists():
        raise FileNotFoundError(h5_path)
    if not source_image.exists():
        raise FileNotFoundError(source_image)

    with h5py.File(h5_path, "r") as handle:
        barcodes = [x.decode() if isinstance(x, bytes) else str(x)
                    for x in handle["matrix/barcodes"][:]]
    parsed = [BARCODE_PATTERN.match(barcode) for barcode in barcodes]
    if any(match is None for match in parsed):
        raise ValueError("Visium HD barcode format is not recognized")
    sizes = {int(match.group(1)) for match in parsed if match is not None}
    if len(sizes) != 1:
        raise ValueError(f"Mixed Visium HD bin sizes: {sorted(sizes)}")
    bin_size = sizes.pop()
    x_um = np.array([int(match.group(3)) * bin_size + bin_size / 2
                     for match in parsed], dtype=np.float64)
    y_um = np.array([int(match.group(2)) * bin_size + bin_size / 2
                     for match in parsed], dtype=np.float64)

    with Image.open(source_image) as image:
        width, height = image.size
        hires = np.asarray(image.convert("RGB"))
        lowres = np.asarray(image.resize((600, 600)))
    # Fit the barcode coordinate extent inside the supplied morphology image.
    # The resulting scale is recorded explicitly and used by bin2cell.
    pixel_per_um = min((width - 1) / max(x_um), (height - 1) / max(y_um))
    mpp_source = 1.0 / pixel_per_um
    px_x = x_um * pixel_per_um
    px_y = y_um * pixel_per_um

    compat_dir = work_dir / f"visium_hd_compat_{bin_size:03d}um"
    spatial_dir = compat_dir / "spatial"
    spatial_dir.mkdir(parents=True, exist_ok=True)
    link_h5 = compat_dir / "filtered_feature_bc_matrix.h5"
    if not link_h5.exists():
        link_h5.symlink_to(h5_path)
    link_hires = spatial_dir / "tissue_hires_image.png"
    if not link_hires.exists():
        link_hires.symlink_to(source_image)
    lowres_path = spatial_dir / "tissue_lowres_image.png"
    if not lowres_path.exists():
        Image.fromarray(lowres).save(lowres_path)

    positions = pd.DataFrame({
        "barcode": barcodes,
        "in_tissue": 1,
        "array_row": [int(match.group(2)) for match in parsed],
        "array_col": [int(match.group(3)) for match in parsed],
        "pxl_col_in_fullres": px_x,
        "pxl_row_in_fullres": px_y,
    }).set_index("barcode")
    positions.to_csv(spatial_dir / "tissue_positions.csv")
    (spatial_dir / "scalefactors_json.json").write_text(json.dumps({
        "microns_per_pixel": mpp_source,
        "tissue_hires_scalef": 1.0,
        "tissue_lowres_scalef": 600.0 / max(width, height),
        "spot_diameter_fullres": float(bin_size * pixel_per_um),
    }) + "\n")

    adata = ov.space.read_visium_10x(
        str(compat_dir),
        source_image_path=str(source_image),
    )
    adata.var_names_make_unique()
    adata.obsm["spatial_um"] = np.column_stack([x_um, y_um]).astype(np.float32)
    adata.uns["spatial_coordinate_source"] = "visium_hd_barcode_scaled_to_morphology"
    adata.uns["spatial_coordinate_unit"] = "fullres_pixels"
    adata.uns["visium_hd_source_mpp"] = float(mpp_source)
    return adata


def crop_visium_hd_um(adata, x_range_um, y_range_um, output_dir: Path):
    """Crop a loaded HD object and its morphology image using micrometer bounds."""
    from PIL import Image

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if "spatial_um" not in adata.obsm:
        raise KeyError("spatial_um is required for micrometer-based cropping")
    coords_um = np.asarray(adata.obsm["spatial_um"])
    x0, x1 = map(float, x_range_um)
    y0, y1 = map(float, y_range_um)
    mask = ((coords_um[:, 0] >= x0) & (coords_um[:, 0] < x1) &
            (coords_um[:, 1] >= y0) & (coords_um[:, 1] < y1))
    cropped = adata[mask].copy()
    library = list(cropped.uns["spatial"].keys())[0]
    spatial = np.asarray(cropped.obsm["spatial"], dtype=np.float32)
    source_mpp = float(cropped.uns["spatial"][library]["scalefactors"]["microns_per_pixel"])
    px_per_um = 1.0 / source_mpp
    px_x0, px_y0 = int(np.floor(x0 * px_per_um)), int(np.floor(y0 * px_per_um))
    px_x1, px_y1 = int(np.ceil(x1 * px_per_um)), int(np.ceil(y1 * px_per_um))
    spatial[:, 1] -= px_x0
    spatial[:, 0] -= px_y0
    cropped.obsm["spatial"] = spatial
    cropped.obsm["spatial_um"] = coords_um[mask] - np.array([x0, y0])

    metadata = cropped.uns["spatial"][library]
    source_path = Path(metadata["metadata"]["source_image_path"])
    with Image.open(source_path) as image:
        image = image.convert("RGB")
        px_x1 = min(px_x1, image.width)
        px_y1 = min(px_y1, image.height)
        crop_image = image.crop((px_x0, px_y0, px_x1, px_y1))
    crop_path = output_dir / "cropped_source_image.png"
    crop_image.save(crop_path)
    cropped.uns["spatial"][library]["images"]["hires"] = np.asarray(crop_image)
    cropped.uns["spatial"][library]["images"]["lowres"] = np.asarray(
        crop_image.resize((max(1, crop_image.width // 10), max(1, crop_image.height // 10)))
    )
    cropped.uns["spatial"][library]["scalefactors"]["tissue_hires_scalef"] = 1.0
    cropped.uns["spatial"][library]["scalefactors"]["tissue_lowres_scalef"] = 0.1
    cropped.uns["spatial"][library]["metadata"]["source_image_path"] = str(crop_path)
    cropped.uns["crop_um"] = {"x_range": [x0, x1], "y_range": [y0, y1]}
    return cropped


def candidate_qc_thresholds(obs) -> dict[str, float]:
    positive = obs.loc[obs["total_counts"] > 0]
    if positive.empty:
        raise ValueError("No positive-count bins were found")
    return {
        "min_counts": float(max(1, np.floor(positive["total_counts"].quantile(0.01)))),
        "max_counts": float(np.ceil(positive["total_counts"].quantile(0.995))),
        "min_genes": float(max(1, np.floor(positive["n_genes_by_counts"].quantile(0.01)))),
        "max_genes": float(np.ceil(positive["n_genes_by_counts"].quantile(0.995))),
        "max_pct_mt": float(min(100, max(20, positive["pct_counts_mt"].quantile(0.995)))),
    }


def apply_qc_flags(adata, thresholds: dict[str, float]) -> None:
    obs = adata.obs
    obs["passes_qc_candidate_counts"] = (
        (obs["total_counts"] >= thresholds["min_counts"])
        & (obs["total_counts"] <= thresholds["max_counts"])
    )
    obs["passes_qc_candidate_genes"] = (
        (obs["n_genes_by_counts"] >= thresholds["min_genes"])
        & (obs["n_genes_by_counts"] <= thresholds["max_genes"])
    )
    obs["passes_qc_candidate_mt"] = obs["pct_counts_mt"] <= thresholds["max_pct_mt"]
    obs["passes_qc_candidate"] = (
        obs["passes_qc_candidate_counts"]
        & obs["passes_qc_candidate_genes"]
        & obs["passes_qc_candidate_mt"]
    )


def graph_component_diagnostics(connectivities) -> tuple[dict[str, Any], np.ndarray]:
    from scipy.sparse.csgraph import connected_components

    graph = connectivities.get() if hasattr(connectivities, "get") else connectivities
    n_components, labels = connected_components(graph, directed=False)
    sizes = np.bincount(labels, minlength=n_components)
    diagnostics = {
        "n_components": int(n_components),
        "largest_component_size": int(sizes.max()),
        "largest_component_fraction": float(sizes.max() / graph.shape[0]),
        "components_lt_10": int((sizes < 10).sum()),
        "components_lt_50": int((sizes < 50).sum()),
        "components_lt_100": int((sizes < 100).sum()),
    }
    return diagnostics, sizes


def load_marker_groups(path: Path) -> dict[str, list[str]]:
    groups: dict[str, list[str]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            groups.setdefault(row["category"], []).append(row["gene"])
    return groups


def deterministic_indices(n_obs: int, max_points: int, seed: int = 0) -> np.ndarray:
    if n_obs <= max_points:
        return np.arange(n_obs)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n_obs, max_points, replace=False))
