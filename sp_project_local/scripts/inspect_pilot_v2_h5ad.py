#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad


def inspect(path: Path) -> dict:
    obj = ad.read_h5ad(path, backed="r")
    return {
        "path": str(path),
        "shape": [int(obj.n_obs), int(obj.n_vars)],
        "obs_columns": sorted(map(str, obj.obs_keys())),
        "layers": sorted(map(str, obj.layers.keys())),
        "obsm": sorted(map(str, obj.obsm_keys())),
        "obsp": sorted(map(str, obj.obsp.keys())),
        "has_raw": obj.raw is not None,
        "spatial_shape": list(obj.obsm["spatial"].shape) if "spatial" in obj.obsm else None,
        "umap_shape": list(obj.obsm["X_umap"].shape) if "X_umap" in obj.obsm else None,
        "pca_shape": list(obj.obsm["X_pca"].shape) if "X_pca" in obj.obsm else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = {"files": [inspect(path) for path in args.input]}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
