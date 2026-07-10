#!/usr/bin/env python3
from __future__ import annotations

import scanpy as sc

from sp_utils import PROJECT_HOME, PROJECT_WORK, preprocess_basic, read_manifest


def main() -> None:
    adatas = []
    for sample in read_manifest(PROJECT_HOME / "config" / "sample_manifest.csv"):
        sample_id = sample["sample_id"]
        path = PROJECT_WORK / sample_id / f"{sample_id}_008um.h5ad"
        if path.exists():
            adata = sc.read_h5ad(path)
            adata.obs["sample_id"] = sample_id
            adata.obs["sample_group"] = sample["group"]
            adatas.append(adata)
    if not adatas:
        raise SystemExit("No per-sample 008um h5ad files found")
    merged = sc.concat(adatas, join="outer", label="sample_id_from_concat", fill_value=0)
    merged.var_names_make_unique()
    merged = preprocess_basic(merged, n_top_genes=4000, n_pcs=40, resolution=0.8)
    out = PROJECT_HOME / "results" / "integrated_visium_hd.h5ad"
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5ad(out)
    import matplotlib

    matplotlib.use("Agg")
    sc.pl.umap(merged, color=["sample_id", "sample_group", "leiden"], show=False)
    import matplotlib.pyplot as plt

    plt.savefig(PROJECT_HOME / "results" / "integrated_umap.png", dpi=180, bbox_inches="tight")
    print(f"Wrote integrated object to {out}")


if __name__ == "__main__":
    main()
