#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gc
import json
import os
import time
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from sp_v2_utils import (
    PROJECT_HOME,
    RESULTS_V2,
    WORK_V2,
    attach_barcode_spatial,
    file_fingerprint,
    graph_component_diagnostics,
    load_config,
    load_marker_groups,
    read_manifest_sample,
    write_h5ad_atomic,
    write_json_atomic,
)


def package_versions(names: list[str]) -> dict[str, str]:
    result = {}
    for name in names:
        try:
            result[name] = version(name)
        except PackageNotFoundError:
            result[name] = "missing"
    return result


def gpu_memory() -> dict[str, float]:
    import cupy as cp

    free, total = cp.cuda.runtime.memGetInfo()
    return {
        "free_gib": round(free / 1024**3, 2),
        "used_gib": round((total - free) / 1024**3, 2),
        "total_gib": round(total / 1024**3, 2),
    }


def normalize_full(adata):
    nonzero = np.asarray(adata.obs["total_counts"] > 0)
    if not nonzero.all():
        adata = adata[nonzero].copy()
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat")
    return adata


def map_labels_to_008um(sample: dict[str, str], full16, label_columns: list[str], output: Path):
    mapping_path = Path(sample["outs_dir"]) / "barcode_mappings.parquet"
    mapping = pd.read_parquet(mapping_path, columns=["square_008um", "square_016um"]).drop_duplicates()
    conflict = mapping.groupby("square_008um", observed=True)["square_016um"].nunique()
    if int((conflict > 1).sum()) > 0:
        raise RuntimeError("At least one 8 µm barcode maps to multiple 16 µm barcodes")
    mapping = mapping.drop_duplicates("square_008um").set_index("square_008um")
    labels = full16.obs[label_columns].copy()
    labels.index.name = "square_016um"
    mapped = mapping.join(labels, on="square_016um")

    adata8 = sc.read_10x_h5(sample["matrix_008um"], gex_only=True)
    adata8.var_names_make_unique()
    adata8.obs["sample_id"] = full16.obs["sample_id"].iloc[0]
    attach_barcode_spatial(adata8)
    for column in label_columns:
        adata8.obs[column] = mapped[column].reindex(adata8.obs_names).astype("category")
    adata8.obs["mapped_square_016um"] = mapped["square_016um"].reindex(adata8.obs_names)
    adata8.uns["mapping_source"] = str(mapping_path)
    adata8.uns["mapping_strategy"] = "deduplicated square_008um to square_016um"
    write_h5ad_atomic(adata8, output)
    mapped.reset_index().to_csv(output.with_suffix(".labels.csv.gz"), index=False, compression="gzip")
    return {
        "n_008um_bins": adata8.n_obs,
        "mapped_fraction": float(adata8.obs[label_columns[0]].notna().mean()),
        "output": str(output),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=PROJECT_HOME / "config" / "pilot_v2.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    sample_id = config["sample_id"]
    outdir = RESULTS_V2 / sample_id / "model"
    workdir = WORK_V2 / sample_id
    state_path = workdir / "run_state.json"
    checkpoint = workdir / "016um_raw_qc.h5ad"
    outdir.mkdir(parents=True, exist_ok=True)
    if not checkpoint.exists():
        raise FileNotFoundError(f"QC checkpoint is missing: {checkpoint}")

    started = time.time()
    state = json.loads(state_path.read_text()) if state_path.exists() else {}
    state.update({
        "stage": "model",
        "status": "running",
        "checkpoint": file_fingerprint(checkpoint),
        "versions": package_versions([
            "scanpy", "anndata", "cupy", "rapids-singlecell", "cudf", "cuml",
            "dask", "distributed", "python-igraph", "leidenalg",
        ]),
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES", ""),
    })
    write_json_atomic(state_path, state)

    try:
        import cupy as cp
        import rapids_singlecell as rsc

        if cp.cuda.runtime.getDeviceCount() < 1:
            raise RuntimeError("No CUDA GPU is visible")
        state["gpu_memory_start"] = gpu_memory()

        full = sc.read_h5ad(checkpoint)
        full = normalize_full(full)
        full.uns["pilot_v2_config"] = json.dumps(config, sort_keys=True)
        normalized_checkpoint = workdir / "016um_normalized_full.h5ad"
        write_h5ad_atomic(full, normalized_checkpoint)

        hvg_mask = full.var["highly_variable"].to_numpy(dtype=bool)
        model = full[:, hvg_mask].copy()
        model.layers.clear()
        rsc.get.anndata_to_GPU(model)
        rsc.pp.scale(model, max_value=10, zero_center=False)
        n_pcs = min(config["n_pcs"], model.n_vars - 1, model.n_obs - 1)
        rsc.pp.pca(model, n_comps=n_pcs, zero_center=False)

        graph_rows = []
        component_rows = []
        selected_neighbors = None
        for n_neighbors in config["neighbor_candidates"]:
            for key in ["connectivities", "distances"]:
                if key in model.obsp:
                    del model.obsp[key]
            rsc.pp.neighbors(
                model,
                n_neighbors=n_neighbors,
                n_pcs=n_pcs,
                algorithm=config["neighbor_algorithm"],
            )
            diagnostics, sizes = graph_component_diagnostics(model.obsp["connectivities"])
            diagnostics["n_neighbors"] = n_neighbors
            graph_rows.append(diagnostics)
            component_rows.extend(
                {"n_neighbors": n_neighbors, "component": i, "size": int(size)}
                for i, size in enumerate(np.sort(sizes)[::-1])
            )
            if diagnostics["largest_component_fraction"] >= config["largest_component_min_fraction"]:
                selected_neighbors = n_neighbors
                break
        pd.DataFrame(graph_rows).to_csv(outdir / "graph_diagnostics.csv", index=False)
        pd.DataFrame(component_rows).to_csv(outdir / "component_sizes.csv", index=False)
        write_json_atomic(outdir / "graph_diagnostics.json", graph_rows)
        if selected_neighbors is None:
            raise RuntimeError("No neighbor graph passed the 99% largest-component quality gate")

        rsc.tl.umap(model, init_pos="random", random_state=config["random_state"])
        leiden_backend = "rapids_singlecell"
        try:
            for resolution in config["leiden_resolutions"]:
                key = f"leiden_r{str(resolution).replace('.', '_')}"
                rsc.tl.leiden(model, resolution=resolution, key_added=key, random_state=config["random_state"])
        except Exception as exc:
            state["gpu_leiden_error"] = repr(exc)
            rsc.get.anndata_to_CPU(model)
            leiden_backend = "scanpy_igraph"
            for resolution in config["leiden_resolutions"]:
                key = f"leiden_r{str(resolution).replace('.', '_')}"
                sc.tl.leiden(
                    model,
                    resolution=resolution,
                    key_added=key,
                    flavor="igraph",
                    directed=False,
                    n_iterations=2,
                    random_state=config["random_state"],
                )
        else:
            rsc.get.anndata_to_CPU(model)

        label_columns = [f"leiden_r{str(value).replace('.', '_')}" for value in config["leiden_resolutions"]]
        full.obsm["X_pca"] = np.asarray(model.obsm["X_pca"])
        full.obsm["X_umap"] = np.asarray(model.obsm["X_umap"])
        full.obsp["connectivities"] = model.obsp["connectivities"]
        full.obsp["distances"] = model.obsp["distances"]
        for column in label_columns:
            full.obs[column] = model.obs[column].astype("category")
        full.uns["pilot_v2_selected_neighbors"] = selected_neighbors
        full.uns["pilot_v2_leiden_backend"] = leiden_backend

        markers = load_marker_groups(PROJECT_HOME / "config" / "prostate_markers_v2.csv")
        marker_rows = []
        for category, genes in markers.items():
            present = [gene for gene in genes if gene in full.var_names]
            marker_rows.extend({"category": category, "gene": gene, "present": gene in present} for gene in genes)
            if present:
                sc.tl.score_genes(full, present, score_name=f"score_{category}", random_state=config["random_state"])
        pd.DataFrame(marker_rows).to_csv(outdir / "marker_availability.csv", index=False)

        final16 = outdir / f"{sample_id}_016um_v2_full.h5ad"
        write_h5ad_atomic(full, final16)
        sample = read_manifest_sample(PROJECT_HOME / "config" / "sample_manifest.csv", sample_id)
        mapping_summary = map_labels_to_008um(
            sample,
            full,
            label_columns,
            outdir / f"{sample_id}_008um_v2_mapped.h5ad",
        )

        cluster_rows = []
        for resolution, column in zip(config["leiden_resolutions"], label_columns):
            sizes = full.obs[column].value_counts()
            cluster_rows.append({
                "resolution": resolution,
                "column": column,
                "n_clusters": int(len(sizes)),
                "min_cluster_size": int(sizes.min()),
                "median_cluster_size": float(sizes.median()),
                "clusters_lt_100": int((sizes < 100).sum()),
                "fraction_bins_in_clusters_lt_100": float(sizes[sizes < 100].sum() / full.n_obs),
            })
        pd.DataFrame(cluster_rows).to_csv(outdir / "clustering_resolution_summary.csv", index=False)
        summary = {
            "sample_id": sample_id,
            "status": "complete",
            "shape_016um": [full.n_obs, full.n_vars],
            "n_hvg": int(hvg_mask.sum()),
            "selected_neighbors": selected_neighbors,
            "leiden_backend": leiden_backend,
            "final_016um": str(final16),
            "mapping_008um": mapping_summary,
            "gpu_memory_end": gpu_memory(),
        }
        write_json_atomic(outdir / "model_summary.json", summary)
        state.update({"status": "complete", "elapsed_seconds": round(time.time() - started, 2), "model": summary})
        write_json_atomic(state_path, state)
        print(json.dumps(summary, indent=2))
        del model, full
        gc.collect()
    except Exception as exc:
        state.update({"status": "failed", "error": repr(exc), "elapsed_seconds": round(time.time() - started, 2)})
        write_json_atomic(state_path, state)
        raise


if __name__ == "__main__":
    main()
