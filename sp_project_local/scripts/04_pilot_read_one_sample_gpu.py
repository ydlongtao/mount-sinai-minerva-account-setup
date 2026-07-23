#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gc
import os
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator

from sp_utils import (
    PILOT_SAMPLE,
    PROJECT_HOME,
    attach_spatial_from_10x,
    load_markers,
    qc_basic,
    read_10x_h5,
    read_manifest,
    save_basic_plots,
    write_json,
)


GPU_WORK = Path(os.environ.get("SP_PROJECT_WORK_GPU", str(PROJECT_HOME / "work_gpu")))


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def gpu_memory() -> dict[str, float]:
    import cupy as cp

    free, total = cp.cuda.runtime.memGetInfo()
    return {
        "free_gib": round(free / 1024**3, 2),
        "used_gib": round((total - free) / 1024**3, 2),
        "total_gib": round(total / 1024**3, 2),
    }


def init_gpu_backend() -> dict[str, Any]:
    import cupy as cp
    import rapids_singlecell as rsc

    device_count = cp.cuda.runtime.getDeviceCount()
    if device_count < 1:
        raise RuntimeError("No CUDA device is visible to CuPy")

    status: dict[str, Any] = {
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES", ""),
        "cupy_devices": device_count,
        "cupy": getattr(cp, "__version__", "unknown"),
        "rapids_singlecell": getattr(rsc, "__version__", "unknown"),
        "gpu_memory_at_start": gpu_memory(),
    }
    try:
        import omicverse as ov

        ov.settings.cpu_gpu_mixed_init()
        status["omicverse_gpu_init"] = "ok"
    except Exception as exc:
        status["omicverse_gpu_init"] = f"failed: {exc!r}"

    log(f"GPU backend ready: {status}")
    return status


@contextmanager
def timed_stage(report: dict[str, Any], report_path: Path, name: str) -> Iterator[None]:
    started = time.perf_counter()
    log(f"START {name}")
    try:
        yield
    except Exception as exc:
        elapsed = round(time.perf_counter() - started, 2)
        report.setdefault("stages", {})[name] = {
            "status": "failed",
            "seconds": elapsed,
            "error": repr(exc),
        }
        write_json(report_path, report)
        log(f"FAILED {name} after {elapsed}s: {exc!r}")
        raise
    else:
        elapsed = round(time.perf_counter() - started, 2)
        report.setdefault("stages", {})[name] = {"status": "complete", "seconds": elapsed}
        write_json(report_path, report)
        log(f"DONE {name} in {elapsed}s")


def write_gpu_checkpoint(adata, path: Path) -> None:
    import rapids_singlecell as rsc

    log(f"Writing checkpoint {path}")
    cpu_copy = rsc.get.anndata_to_CPU(adata, copy=True)
    cpu_copy.write_h5ad(path, compression="lzf")
    del cpu_copy
    gc.collect()


def preprocess_gpu(
    adata,
    report: dict[str, Any],
    report_path: Path,
    checkpoint_dir: Path,
    prefix: str,
    n_top_genes: int,
    n_pcs: int,
    neighbor_algorithm: str,
    resolution: float,
):
    import numpy as np
    import rapids_singlecell as rsc
    import scanpy as sc
    from scipy import sparse

    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    with timed_stage(report, report_path, f"{prefix}:cpu_qc_filter"):
        adata.var_names_make_unique()
        qc_basic(adata)
        nonzero = np.asarray(adata.obs["total_counts"] > 0)
        report.setdefault("datasets", {}).setdefault(prefix, {})["zero_count_bins_removed"] = int((~nonzero).sum())
        if not nonzero.all():
            adata = adata[nonzero].copy()
        sc.pp.filter_genes(adata, min_cells=3)
        adata.X = adata.X.astype(np.float32)
        if sparse.issparse(adata.X):
            adata.X = adata.X.tocsr()
            adata.X.sum_duplicates()
            adata.X.sort_indices()
        report["datasets"][prefix].update({"n_obs_after_qc": adata.n_obs, "n_vars_after_qc": adata.n_vars})
        adata.write_h5ad(checkpoint_dir / f"{prefix}_qc_checkpoint.h5ad", compression="lzf")

    with timed_stage(report, report_path, f"{prefix}:to_gpu"):
        rsc.get.anndata_to_GPU(adata)
        report["datasets"][prefix]["gpu_memory_after_transfer"] = gpu_memory()

    with timed_stage(report, report_path, f"{prefix}:normalize_log1p"):
        rsc.pp.normalize_total(adata, target_sum=1e4)
        rsc.pp.log1p(adata)

    with timed_stage(report, report_path, f"{prefix}:hvg"):
        rsc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="seurat")
        hvg_mask = np.asarray(adata.var["highly_variable"].to_numpy(), dtype=bool)
        if hvg_mask.sum() < 2:
            raise RuntimeError(f"Only {hvg_mask.sum()} highly variable genes were selected")
        adata = adata[:, hvg_mask].copy()
        report["datasets"][prefix]["n_hvg"] = int(hvg_mask.sum())

    with timed_stage(report, report_path, f"{prefix}:scale_pca"):
        rsc.pp.scale(adata, max_value=10, zero_center=False)
        actual_pcs = min(n_pcs, max(2, min(adata.n_obs, adata.n_vars) - 1))
        rsc.pp.pca(adata, n_comps=actual_pcs, zero_center=False)
        report["datasets"][prefix]["n_pcs"] = actual_pcs
        report["datasets"][prefix]["gpu_memory_after_pca"] = gpu_memory()
        write_gpu_checkpoint(adata, checkpoint_dir / f"{prefix}_pca_checkpoint.h5ad")

    return finish_gpu_embedding(
        adata,
        report,
        report_path,
        prefix,
        actual_pcs,
        neighbor_algorithm,
        resolution,
    )


def finish_gpu_embedding(
    adata,
    report: dict[str, Any],
    report_path: Path,
    prefix: str,
    actual_pcs: int,
    neighbor_algorithm: str,
    resolution: float,
):
    import rapids_singlecell as rsc
    import scanpy as sc

    with timed_stage(report, report_path, f"{prefix}:neighbors"):
        rsc.pp.neighbors(
            adata,
            n_neighbors=15,
            n_pcs=actual_pcs,
            algorithm=neighbor_algorithm,
        )
        report["datasets"][prefix]["neighbor_algorithm"] = neighbor_algorithm

    with timed_stage(report, report_path, f"{prefix}:umap"):
        rsc.tl.umap(adata, random_state=0)

    with timed_stage(report, report_path, f"{prefix}:leiden"):
        on_gpu = True
        try:
            rsc.tl.leiden(adata, resolution=resolution, key_added="leiden", random_state=0)
            report["datasets"][prefix]["leiden_backend"] = "rapids_singlecell"
        except Exception as exc:
            report.setdefault("warnings", []).append(f"GPU Leiden failed; used CPU leidenalg: {exc!r}")
            rsc.get.anndata_to_CPU(adata)
            on_gpu = False
            sc.tl.leiden(adata, resolution=resolution, key_added="leiden", flavor="leidenalg", random_state=0)
            report["datasets"][prefix]["leiden_backend"] = "scanpy_leidenalg"

    with timed_stage(report, report_path, f"{prefix}:to_cpu"):
        if on_gpu:
            rsc.get.anndata_to_CPU(adata)

    return adata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", default=PILOT_SAMPLE)
    parser.add_argument("--include-segmented", action="store_true")
    parser.add_argument("--resume-pca-checkpoint", action="store_true")
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--neighbor-algorithm", choices=["brute", "ivfflat", "ivfpq", "cagra", "nn_descent"], default="ivfflat")
    parser.add_argument("--resolution", type=float, default=0.6)
    args = parser.parse_args()

    outdir = PROJECT_HOME / "results" / "pilot_gpu" / args.sample_id
    workdir = GPU_WORK / args.sample_id
    outdir.mkdir(parents=True, exist_ok=True)
    workdir.mkdir(parents=True, exist_ok=True)
    report_path = outdir / "pilot_gpu_report.json"

    report: dict[str, Any] = {
        "sample_id": args.sample_id,
        "status": "running",
        "include_segmented": args.include_segmented,
        "parameters": vars(args),
        "warnings": [],
        "stages": {},
        "datasets": {},
    }
    write_json(report_path, report)

    try:
        report["gpu_status"] = init_gpu_backend()
        sample = next(
            row
            for row in read_manifest(PROJECT_HOME / "config" / "sample_manifest.csv")
            if row["sample_id"] == args.sample_id
        )
        markers = load_markers(PROJECT_HOME / "config" / "prostate_markers.csv")
        marker_flat = [gene for genes in markers.values() for gene in genes]

        pca_checkpoint = workdir / "008um_pca_checkpoint.h5ad"
        if args.resume_pca_checkpoint:
            if not pca_checkpoint.exists():
                raise FileNotFoundError(f"PCA checkpoint not found: {pca_checkpoint}")
            with timed_stage(report, report_path, "008um:resume_pca_checkpoint"):
                import anndata as ad

                adata = ad.read_h5ad(pca_checkpoint)
                report["warnings"].extend(attach_spatial_from_10x(adata, Path(sample["matrix_008um"]).parent))
                if "X_pca" not in adata.obsm:
                    raise KeyError(f"X_pca missing from {pca_checkpoint}")
                actual_pcs = min(args.n_pcs, adata.obsm["X_pca"].shape[1])
                report["datasets"]["008um"] = {
                    "resume_checkpoint": str(pca_checkpoint),
                    "checkpoint_shape": list(adata.shape),
                    "n_pcs": actual_pcs,
                    "spatial_coordinate_source": adata.uns.get("spatial_coordinate_source", "missing"),
                }
            with timed_stage(report, report_path, "008um:resume_to_gpu"):
                import rapids_singlecell as rsc

                rsc.get.anndata_to_GPU(adata)
                report["datasets"]["008um"]["gpu_memory_after_resume"] = gpu_memory()
            adata = finish_gpu_embedding(
                adata,
                report,
                report_path,
                "008um",
                actual_pcs,
                args.neighbor_algorithm,
                args.resolution,
            )
        else:
            with timed_stage(report, report_path, "008um:read_attach_spatial"):
                adata = read_10x_h5(Path(sample["matrix_008um"]), args.sample_id)
                report["warnings"].extend(attach_spatial_from_10x(adata, Path(sample["matrix_008um"]).parent))
                report["datasets"]["008um"] = {"input_shape": list(adata.shape)}

            adata = preprocess_gpu(
                adata,
                report,
                report_path,
                workdir,
                "008um",
                args.n_top_genes,
                args.n_pcs,
                args.neighbor_algorithm,
                args.resolution,
            )

        with timed_stage(report, report_path, "008um:plots_final_write"):
            report["warnings"].extend(save_basic_plots(adata, outdir, f"{args.sample_id}_008um_gpu", marker_flat))
            final_path = outdir / f"{args.sample_id}_008um_gpu_pilot.h5ad"
            adata.write_h5ad(final_path, compression="lzf")
            report["datasets"]["008um"]["final_h5ad"] = str(final_path)
        del adata
        gc.collect()

        if args.include_segmented:
            with timed_stage(report, report_path, "segmented:read"):
                seg = read_10x_h5(Path(sample["segmented_h5"]), args.sample_id)
                report["datasets"]["segmented"] = {"input_shape": list(seg.shape)}
            seg = preprocess_gpu(
                seg,
                report,
                report_path,
                workdir,
                "segmented",
                args.n_top_genes,
                args.n_pcs,
                args.neighbor_algorithm,
                args.resolution,
            )
            with timed_stage(report, report_path, "segmented:plots_final_write"):
                report["warnings"].extend(save_basic_plots(seg, outdir, f"{args.sample_id}_segmented_gpu", marker_flat))
                final_path = outdir / f"{args.sample_id}_segmented_gpu_pilot.h5ad"
                seg.write_h5ad(final_path, compression="lzf")
                report["datasets"]["segmented"]["final_h5ad"] = str(final_path)

        report["status"] = "complete"
        write_json(report_path, report)
        log(f"GPU pilot complete for {args.sample_id}")
    except Exception:
        report["status"] = "failed"
        write_json(report_path, report)
        raise


if __name__ == "__main__":
    main()
