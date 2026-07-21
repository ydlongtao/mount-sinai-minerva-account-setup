from __future__ import annotations

import csv
import json
import os
from pathlib import Path
from typing import Any

import numpy as np


PROJECT_HOME = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
BATCH_RESULTS = Path(os.environ.get("SP_PROJECT_BATCH_RESULTS", str(PROJECT_HOME / "results" / "batch")))
BATCH_WORK = Path(os.environ.get("SP_PROJECT_BATCH_WORK", str(PROJECT_HOME / "work" / "batch")))
SCRATCH = Path(os.environ.get("SP_PROJECT_SCRATCH", "/sc/arion/scratch/huangl21/sp_project"))
MANIFEST = Path(os.environ.get("SP_PROJECT_MANIFEST", str(PROJECT_HOME / "config" / "sample_manifest.csv")))
CONFIG = Path(os.environ.get("SP_PROJECT_CONFIG", str(PROJECT_HOME / "config" / "batch_analysis.yaml")))


DEFAULT_MARKERS = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "Prostate_luminal": ["KLK3", "ACPP", "AR", "NKX3-1"],
    "Basal": ["KRT5", "KRT14", "TP63"],
    "Tumor_stress": ["AMACR", "MKI67", "TOP2A", "EGR1"],
    "Immune": ["PTPRC", "CD3D", "CD3E", "MS4A1", "LYZ"],
    "Stromal": ["COL1A1", "COL1A2", "DCN", "LUM", "ACTA2"],
    "Endothelial": ["PECAM1", "VWF", "KDR"],
}


def read_json(path: Path, default: Any | None = None) -> Any:
    if not path.exists():
        return default
    return json.loads(path.read_text())


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def read_manifest(path: Path = MANIFEST) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def sample_by_index(index: int, manifest: Path = MANIFEST) -> dict[str, str]:
    rows = read_manifest(manifest)
    if index < 1 or index > len(rows):
        raise IndexError(f"sample index {index} is outside 1..{len(rows)}")
    return rows[index - 1]


def sample_by_id(sample_id: str, manifest: Path = MANIFEST) -> dict[str, str]:
    for row in read_manifest(manifest):
        if row["sample_id"] == sample_id:
            return row
    raise KeyError(sample_id)


def sample_root(sample_id: str) -> Path:
    return BATCH_RESULTS / "samples" / sample_id


def load_batch_config(path: Path = CONFIG) -> dict[str, Any]:
    data = read_json(path, default={})
    if not data:
        data = {
            "analysis_version": "batch_v1",
            "sample_ids": [row["sample_id"] for row in read_manifest()],
            "crop_defaults": {
                "x_range_um": [2000, 3000],
                "y_range_um": [2000, 3000],
                "mpp": 0.3239732935,
                "buffer": 150,
                "expand_max_bin_distance": 4,
            },
            "crop_by_sample": {},
            "model": {
                "n_top_genes": 3000,
                "n_pcs": 30,
                "n_neighbors": 15,
                "leiden_resolution": 0.4,
                "random_state": 0,
            },
        }
    return data


def crop_config_for_sample(sample_id: str, config: dict[str, Any]) -> dict[str, Any]:
    defaults = dict(config.get("crop_defaults", {}))
    specific = config.get("crop_by_sample", {}).get(sample_id, {})
    defaults.update(specific)
    return defaults


def marker_groups(marker_csv: Path | None = None) -> dict[str, list[str]]:
    if marker_csv is None:
        marker_csv = PROJECT_HOME / "config" / "prostate_markers_v2.csv"
    if not marker_csv.exists():
        return DEFAULT_MARKERS
    groups: dict[str, list[str]] = {}
    with marker_csv.open(newline="") as handle:
        for row in csv.DictReader(handle):
            category = row.get("category") or row.get("group") or row.get("cell_type")
            gene = row.get("gene") or row.get("symbol")
            if category and gene:
                groups.setdefault(category, []).append(gene)
    return groups or DEFAULT_MARKERS


def barcode_spatial(adata) -> None:
    import re

    pattern = re.compile(r"^s_(\d+)um_(\d+)_(\d+)-\d+$")
    parsed = [pattern.match(str(x)) for x in adata.obs_names]
    if any(x is None for x in parsed):
        raise ValueError("Visium HD barcodes are not parseable")
    sizes = {int(x.group(1)) for x in parsed if x is not None}
    if len(sizes) != 1:
        raise ValueError(f"mixed bin sizes: {sorted(sizes)}")
    size = sizes.pop()
    coords = np.empty((adata.n_obs, 2), dtype=np.float32)
    coords[:, 0] = [int(x.group(3)) * size + size / 2 for x in parsed if x is not None]
    coords[:, 1] = [int(x.group(2)) * size + size / 2 for x in parsed if x is not None]
    adata.obsm["spatial"] = coords
    adata.uns["spatial_coordinate_source"] = "visium_hd_barcode"
    adata.uns["spatial_coordinate_unit"] = "micrometer"


def score_markers(adata, groups: dict[str, list[str]] | None = None, random_state: int = 0) -> dict[str, list[str]]:
    import scanpy as sc

    groups = groups or DEFAULT_MARKERS
    upper = {g.upper(): g for g in adata.var_names}
    present: dict[str, list[str]] = {}
    for name, genes in groups.items():
        found = [upper[g.upper()] for g in genes if g.upper() in upper]
        present[name] = found
        if found:
            sc.tl.score_genes(adata, found, score_name=f"score_{name}", random_state=random_state)
    return present


def make_basic_spatial_plot(adata, color: str, output: Path, title: str | None = None, max_points: int = 150000) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    output.parent.mkdir(parents=True, exist_ok=True)
    coords = np.asarray(adata.obsm["spatial"])
    idx = np.arange(adata.n_obs)
    if adata.n_obs > max_points:
        rng = np.random.default_rng(0)
        idx = np.sort(rng.choice(idx, max_points, replace=False))
    values = adata.obs[color].iloc[idx] if color in adata.obs else adata[:, color].X
    fig, ax = plt.subplots(figsize=(7, 7))
    if color in adata.obs and pd.api.types.is_categorical_dtype(adata.obs[color]):
        codes = adata.obs[color].cat.codes.to_numpy()[idx]
        sca = ax.scatter(coords[idx, 0], coords[idx, 1], c=codes, s=1, cmap="tab20", rasterized=True)
    else:
        vals = np.asarray(values).ravel()
        vmax = np.nanpercentile(vals, 99) if vals.size else None
        sca = ax.scatter(coords[idx, 0], coords[idx, 1], c=vals, s=1, cmap="viridis", vmax=vmax, rasterized=True)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title or color)
    fig.colorbar(sca, ax=ax, fraction=0.035)
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
