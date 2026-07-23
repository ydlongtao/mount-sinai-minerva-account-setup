#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


PROJECT_HOME = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
PROJECT_WORK = Path(os.environ.get("SP_PROJECT_WORK", str(PROJECT_HOME / "work")))
PROJECT_SCRATCH = Path(os.environ.get("SP_PROJECT_SCRATCH", "/sc/arion/scratch/huangl21/sp_project"))
RAW_ROOT = Path(os.environ.get("SP_RAW", "/sc/arion/projects/DiseaseGeneCell/Huang_lab_data/SpatialTranscriptome"))
DEFAULT_RESOLUTION = os.environ.get("SP_DEFAULT_RESOLUTION", "square_008um")
PILOT_SAMPLE = os.environ.get("SP_PILOT_SAMPLE", "SC000895-R4")


def project_path(*parts: str) -> Path:
    return PROJECT_HOME.joinpath(*parts)


def ensure_project_dirs() -> None:
    for path in [
        PROJECT_HOME,
        PROJECT_HOME / "scripts",
        PROJECT_HOME / "lsf",
        PROJECT_HOME / "config",
        PROJECT_HOME / "notebooks",
        PROJECT_HOME / "docs",
        PROJECT_HOME / "results",
        PROJECT_HOME / "logs",
        PROJECT_WORK,
        PROJECT_SCRATCH,
    ]:
        path.mkdir(parents=True, exist_ok=True)


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def write_markdown(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def run_cmd(cmd: list[str]) -> dict[str, Any]:
    proc = subprocess.run(cmd, text=True, capture_output=True)
    return {
        "cmd": cmd,
        "returncode": proc.returncode,
        "stdout": proc.stdout.strip(),
        "stderr": proc.stderr.strip(),
    }


def import_status(modules: list[str]) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for name in modules:
        try:
            mod = __import__(name)
            out[name] = {
                "ok": True,
                "version": getattr(mod, "__version__", "unknown"),
                "error": "",
            }
        except Exception as exc:
            out[name] = {"ok": False, "version": "", "error": repr(exc)}
    return out


def h5_matrix_shape(path: Path) -> tuple[int | None, int | None]:
    import h5py

    with h5py.File(path, "r") as h5:
        if "matrix" in h5 and "shape" in h5["matrix"]:
            shape = h5["matrix"]["shape"][:]
            if len(shape) >= 2:
                return int(shape[0]), int(shape[1])
        if "matrix" in h5 and "data" in h5["matrix"]:
            return None, int(len(h5["matrix"]["data"]))
    return None, None


def read_10x_h5(path: Path, sample_id: str | None = None):
    import scanpy as sc

    adata = sc.read_10x_h5(str(path), gex_only=True)
    adata.var_names_make_unique()
    if sample_id:
        adata.obs["sample_id"] = sample_id
    return adata


def _candidate_spatial_files(base: Path) -> list[Path]:
    return [
        base / "spatial" / "tissue_positions.parquet",
        base / "spatial" / "tissue_positions.csv",
        base / "spatial" / "tissue_positions_list.csv",
        base / "tissue_positions.parquet",
        base / "tissue_positions.csv",
    ]


def attach_spatial_from_10x(adata, base: Path) -> list[str]:
    import numpy as np
    import pandas as pd

    warnings: list[str] = []
    pos = next((p for p in _candidate_spatial_files(base) if p.exists()), None)
    if pos is None:
        # Visium HD bin barcodes encode bin size, row, and column directly.
        pattern = re.compile(r"^s_(\d+)um_(\d+)_(\d+)-\d+$")
        parsed = [pattern.match(str(barcode)) for barcode in adata.obs_names]
        if parsed and all(match is not None for match in parsed):
            bin_sizes = {int(match.group(1)) for match in parsed if match is not None}
            if len(bin_sizes) == 1:
                bin_size_um = bin_sizes.pop()
                coords = np.empty((adata.n_obs, 2), dtype=float)
                coords[:, 0] = [int(match.group(3)) * bin_size_um + bin_size_um / 2 for match in parsed if match is not None]
                coords[:, 1] = [int(match.group(2)) * bin_size_um + bin_size_um / 2 for match in parsed if match is not None]
                adata.obsm["spatial"] = coords
                adata.obs["has_spatial"] = True
                adata.uns["spatial_coordinate_source"] = "visium_hd_bin_barcode"
                adata.uns["spatial_coordinate_unit"] = "micrometer"
                return warnings
        warnings.append(f"no tissue_positions file and barcodes were not parseable under {base}")
        return warnings

    if pos.suffix == ".parquet":
        df = pd.read_parquet(pos)
    else:
        df = pd.read_csv(pos, header=None if pos.name == "tissue_positions_list.csv" else "infer")

    cols = list(df.columns)
    barcode_col = "barcode" if "barcode" in cols else cols[0]
    x_candidates = ["pxl_col_in_fullres", "array_col", "x", "imagecol"]
    y_candidates = ["pxl_row_in_fullres", "array_row", "y", "imagerow"]
    x_col = next((c for c in x_candidates if c in cols), None)
    y_col = next((c for c in y_candidates if c in cols), None)
    if x_col is None or y_col is None:
        if len(cols) >= 6:
            y_col, x_col = cols[-2], cols[-1]
        elif len(cols) >= 3:
            x_col, y_col = cols[-2], cols[-1]
        else:
            warnings.append(f"cannot infer spatial columns from {pos}")
            return warnings

    df = df.set_index(barcode_col)
    shared = adata.obs_names.intersection(df.index)
    if len(shared) == 0:
        warnings.append(f"no shared barcodes between AnnData and {pos}")
        return warnings
    coords = np.full((adata.n_obs, 2), np.nan)
    loc = adata.obs_names.get_indexer(shared)
    coords[loc, 0] = df.loc[shared, x_col].to_numpy()
    coords[loc, 1] = df.loc[shared, y_col].to_numpy()
    adata.obsm["spatial"] = coords
    adata.obs["has_spatial"] = ~np.isnan(coords[:, 0])
    return warnings


def qc_basic(adata) -> None:
    import scanpy as sc

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)


def preprocess_basic(adata, n_top_genes: int = 3000, n_pcs: int = 30, resolution: float = 0.6):
    import numpy as np
    import scanpy as sc

    if "total_counts" in adata.obs:
        nonzero = np.asarray(adata.obs["total_counts"] > 0)
        if not nonzero.all():
            adata = adata[nonzero].copy()
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="seurat")
        adata = adata[:, adata.var["highly_variable"]].copy()
    except Exception:
        pass
    # Keep sparse Visium HD matrices from being accidentally densified.
    sc.pp.scale(adata, max_value=10, zero_center=False)
    sc.tl.pca(adata, n_comps=min(n_pcs, max(2, min(adata.n_obs, adata.n_vars) - 1)))
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=min(n_pcs, adata.obsm["X_pca"].shape[1]))
    # Random initialization avoids the costly spectral embedding seen on large HD bin graphs.
    sc.tl.umap(adata, init_pos="random", random_state=0)
    run_leiden(adata, resolution=resolution, key_added="leiden")
    return adata


def run_leiden(adata, resolution: float = 0.6, key_added: str = "leiden") -> str:
    import scanpy as sc

    import leidenalg  # noqa: F401

    sc.tl.leiden(adata, resolution=resolution, key_added=key_added, flavor="leidenalg")
    adata.uns[f"{key_added}_backend"] = "leidenalg"
    return "leidenalg"


def save_basic_plots(adata, outdir: Path, sample_id: str, marker_genes: list[str] | None = None) -> list[str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc

    outdir.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []

    sc.pl.umap(adata, color=["leiden"], show=False)
    plt.savefig(outdir / f"{sample_id}_umap_leiden.png", dpi=180, bbox_inches="tight")
    plt.close("all")

    qc_cols = [c for c in ["total_counts", "n_genes_by_counts", "pct_counts_mt"] if c in adata.obs.columns]
    if qc_cols:
        sc.pl.violin(adata, qc_cols, jitter=0.2, multi_panel=True, show=False)
        plt.savefig(outdir / f"{sample_id}_qc_violin.png", dpi=180, bbox_inches="tight")
        plt.close("all")

    if "spatial" in adata.obsm:
        coords = adata.obsm["spatial"]
        fig, ax = plt.subplots(figsize=(7, 6))
        groups = adata.obs["leiden"].astype(str)
        for group in sorted(groups.unique()):
            mask = groups == group
            ax.scatter(coords[mask, 0], coords[mask, 1], s=1, label=group, alpha=0.8)
        ax.invert_yaxis()
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"{sample_id} spatial Leiden")
        ax.legend(markerscale=5, fontsize=6, bbox_to_anchor=(1.02, 1), loc="upper left")
        fig.savefig(outdir / f"{sample_id}_spatial_leiden.png", dpi=180, bbox_inches="tight")
        plt.close(fig)
    else:
        warnings.append(f"{sample_id}: no spatial coordinates attached")

    if marker_genes:
        present = [g for g in marker_genes if g in adata.var_names]
        if present:
            sc.pl.umap(adata, color=present[:12], show=False)
            plt.savefig(outdir / f"{sample_id}_umap_marker_genes.png", dpi=180, bbox_inches="tight")
            plt.close("all")
    return warnings


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def sample_from_index(index: int, manifest_path: Path) -> dict[str, str]:
    rows = read_manifest(manifest_path)
    if index < 1 or index > len(rows):
        raise SystemExit(f"LSB_JOBINDEX/sample index {index} out of range 1..{len(rows)}")
    return rows[index - 1]


def load_markers(path: Path) -> dict[str, list[str]]:
    markers: dict[str, list[str]] = {}
    if not path.exists():
        return markers
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            markers.setdefault(row["category"], []).append(row["gene"])
    return markers


def sample_group(sample_id: str) -> str:
    if sample_id.startswith("SC000895"):
        return "SC000895_Kuan_lin_Huang_T643"
    if sample_id.startswith("TD006859"):
        return "TD006859_KHuang_T598"
    return "unknown"


def sanitize(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
