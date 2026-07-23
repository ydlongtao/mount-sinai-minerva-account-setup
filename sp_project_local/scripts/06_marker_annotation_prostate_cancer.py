#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import pandas as pd
import scanpy as sc

from sp_utils import PROJECT_HOME, PROJECT_WORK, load_markers, read_manifest


def score_markers(adata, markers: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    for category, genes in markers.items():
        present = [g for g in genes if g in adata.var_names]
        if not present:
            rows.append({"category": category, "present_genes": 0, "top_cluster": "", "mean_score": ""})
            continue
        key = f"score_{category}"
        sc.tl.score_genes(adata, gene_list=present, score_name=key)
        means = adata.obs.groupby("leiden", observed=True)[key].mean().sort_values(ascending=False)
        rows.append(
            {
                "category": category,
                "present_genes": len(present),
                "genes": ";".join(present),
                "top_cluster": str(means.index[0]),
                "mean_score": float(means.iloc[0]),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    markers = load_markers(PROJECT_HOME / "config" / "prostate_markers.csv")
    all_rows = []
    for sample in read_manifest(PROJECT_HOME / "config" / "sample_manifest.csv"):
        sample_id = sample["sample_id"]
        for label in ["008um", "segmented"]:
            path = PROJECT_WORK / sample_id / f"{sample_id}_{label}.h5ad"
            if not path.exists():
                continue
            adata = sc.read_h5ad(path)
            df = score_markers(adata, markers)
            df.insert(0, "assay_resolution", label)
            df.insert(0, "sample_id", sample_id)
            all_rows.append(df)
    out = PROJECT_HOME / "results" / "prostate_marker_scores.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    if all_rows:
        pd.concat(all_rows, ignore_index=True).to_csv(out, index=False)
    else:
        pd.DataFrame(
            columns=["sample_id", "assay_resolution", "category", "present_genes", "genes", "top_cluster", "mean_score"]
        ).to_csv(out, index=False)
    print(f"Wrote marker scores to {out}")


if __name__ == "__main__":
    main()
