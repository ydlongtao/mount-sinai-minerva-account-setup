#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from batch_utils import BATCH_RESULTS, CONFIG, PROJECT_HOME, load_batch_config, read_manifest, write_json


def main() -> None:
    parser = argparse.ArgumentParser(description="Create/validate a Visium HD batch configuration.")
    parser.add_argument("--output", type=Path, default=CONFIG)
    args = parser.parse_args()

    rows = read_manifest()
    config = load_batch_config(args.output)
    config["sample_ids"] = [row["sample_id"] for row in rows]
    config["n_samples"] = len(rows)
    config.setdefault("analysis_version", "batch_v1")
    config.setdefault("crop_defaults", {
        "x_range_um": [2000, 3000],
        "y_range_um": [2000, 3000],
        "mpp": 0.3239732935,
        "buffer": 150,
        "expand_max_bin_distance": 4,
    })
    config.setdefault("crop_by_sample", {})
    config.setdefault("model", {
        "n_top_genes": 3000,
        "n_pcs": 30,
        "n_neighbors": 15,
        "leiden_resolution": 0.4,
        "random_state": 0,
    })
    write_json(args.output, config)
    write_json(BATCH_RESULTS / "config" / "config_snapshot.json", config)
    write_json(BATCH_RESULTS / "config" / "manifest_snapshot.json", rows)
    print(f"Wrote {args.output}")
    print(f"Samples: {len(rows)}")


if __name__ == "__main__":
    main()
