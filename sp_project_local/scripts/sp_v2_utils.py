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
