#!/usr/bin/env python3
from __future__ import annotations

import os
import platform
import sys
from pathlib import Path

from sp_utils import (
    PILOT_SAMPLE,
    PROJECT_HOME,
    PROJECT_SCRATCH,
    PROJECT_WORK,
    RAW_ROOT,
    attach_spatial_from_10x,
    h5_matrix_shape,
    import_status,
    read_10x_h5,
    run_leiden,
    write_json,
    write_markdown,
)


def find_sample(sample_id: str) -> Path:
    hits = sorted(RAW_ROOT.glob(f"**/{sample_id}/outs"))
    if not hits:
        raise FileNotFoundError(f"Cannot find {sample_id}/outs under {RAW_ROOT}")
    return hits[0]


def main() -> None:
    outdir = PROJECT_HOME / "results" / "preflight"
    outdir.mkdir(parents=True, exist_ok=True)
    report: dict = {
        "status": "pass",
        "warnings": [],
        "errors": [],
        "python": sys.executable,
        "python_version": sys.version,
        "platform": platform.platform(),
        "hostname": platform.node(),
        "conda_default_env": os.environ.get("CONDA_DEFAULT_ENV", ""),
        "lsf_jobid": os.environ.get("LSB_JOBID", ""),
        "paths": {
            "project_home": str(PROJECT_HOME),
            "project_work": str(PROJECT_WORK),
            "project_scratch": str(PROJECT_SCRATCH),
            "raw_root": str(RAW_ROOT),
        },
    }

    core = [
        "omicverse",
        "scanpy",
        "anndata",
        "numpy",
        "pandas",
        "scipy",
        "sklearn",
        "matplotlib",
        "seaborn",
        "h5py",
        "igraph",
        "leidenalg",
        "squidpy",
        "spatialdata",
    ]
    optional = ["tifffile", "PIL"]
    report["imports"] = import_status(core + optional)
    missing_core = [m for m in core if not report["imports"][m]["ok"]]
    if missing_core:
        report["status"] = "fail"
        report["errors"].append(f"Missing core modules: {', '.join(missing_core)}")
    missing_optional = [m for m in optional if not report["imports"][m]["ok"]]
    if missing_optional:
        report["warnings"].append(f"Missing optional modules: {', '.join(missing_optional)}")

    for path_name, path in [("work", PROJECT_WORK), ("scratch", PROJECT_SCRATCH), ("results", PROJECT_HOME / "results")]:
        try:
            path.mkdir(parents=True, exist_ok=True)
            probe = path / f"preflight_write_{os.getpid()}.txt"
            probe.write_text("ok\n")
            probe.unlink()
        except Exception as exc:
            report["status"] = "fail"
            report["errors"].append(f"{path_name} not writable: {path}: {exc!r}")

    try:
        outs = find_sample(PILOT_SAMPLE)
        h5_008 = outs / "binned_outputs" / "square_008um" / "filtered_feature_bc_matrix.h5"
        seg = outs / "segmented_outputs" / "filtered_feature_cell_matrix.h5"
        report["pilot_outs"] = str(outs)
        report["h5_shapes"] = {
            "square_008um": h5_matrix_shape(h5_008),
            "segmented": h5_matrix_shape(seg),
        }
        adata = read_10x_h5(h5_008, PILOT_SAMPLE)
        report["read_008um"] = {"n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars)}
        report["spatial_warnings"] = attach_spatial_from_10x(adata, h5_008.parent)
        smoke = adata[: min(200, adata.n_obs), : min(500, adata.n_vars)].copy()
        import scanpy as sc

        sc.pp.filter_genes(smoke, min_cells=1)
        sc.pp.normalize_total(smoke, target_sum=1e4)
        sc.pp.log1p(smoke)
        sc.tl.pca(smoke, n_comps=min(10, max(2, min(smoke.n_obs, smoke.n_vars) - 1)))
        sc.pp.neighbors(smoke, n_neighbors=min(10, max(2, smoke.n_obs - 1)), n_pcs=smoke.obsm["X_pca"].shape[1])
        report["leiden_backend"] = run_leiden(smoke, resolution=0.2, key_added="leiden_smoke")
        test_h5ad = PROJECT_WORK / "preflight" / "preflight_test.h5ad"
        test_h5ad.parent.mkdir(parents=True, exist_ok=True)
        adata[: min(100, adata.n_obs), : min(100, adata.n_vars)].write_h5ad(test_h5ad)

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1])
        fig.savefig(outdir / "matplotlib_smoke_test.png")
        plt.close(fig)
    except Exception as exc:
        report["status"] = "fail"
        report["errors"].append(f"Data read or plotting preflight failed: {exc!r}")

    write_json(outdir / "environment_report.json", report)
    lines = [
        "# Environment Preflight Report",
        "",
        f"- Status: **{report['status']}**",
        f"- Host: `{report['hostname']}`",
        f"- Python: `{report['python']}`",
        f"- Conda env: `{report['conda_default_env']}`",
        "",
        "## Errors",
        "",
        *(f"- {e}" for e in report["errors"]),
        "",
        "## Warnings",
        "",
        *(f"- {w}" for w in report["warnings"]),
        "",
        "## Module Imports",
        "",
    ]
    for name, info in report["imports"].items():
        lines.append(f"- `{name}`: {'OK' if info['ok'] else 'MISSING'} {info.get('version', '')}")
    write_markdown(outdir / "environment_report.md", lines)
    print(f"Preflight status: {report['status']}")
    if report["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
