#!/usr/bin/env python3
"""Stage validated v1 sample outputs read-only into the 8-sample v2 namespace."""
from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

from batch_utils import BATCH_RESULTS, read_manifest


PROJECT_HOME = Path(os.environ.get("SP_PROJECT_HOME", "/sc/arion/work/huangl21/sp_project"))
SOURCE_RESULTS = Path(os.environ.get("SP_PROJECT_INPUT_BATCH_RESULTS", str(PROJECT_HOME / "results" / "batch")))
EXCLUDED = {"TD006859-B408", "TD006859-B573"}


def main() -> None:
    rows = read_manifest()
    if any(row["sample_id"] in EXCLUDED for row in rows):
        raise RuntimeError("The active v2 manifest still contains an excluded legacy sample")
    staged = []
    for row in rows:
        sample_id = row["sample_id"]
        source = SOURCE_RESULTS / "samples" / sample_id
        target = BATCH_RESULTS / "samples" / sample_id
        if not source.is_dir():
            raise FileNotFoundError(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.is_symlink() and target.resolve() == source.resolve():
            target.unlink()
        elif target.exists():
            raise FileExistsError(f"Refusing to replace existing v2 path: {target}")
        target.mkdir()
        for name in ("002um", "008um", "016um", "cellpose", "cell_qc_marker"):
            input_dir = source / name
            if input_dir.is_dir():
                (target / name).symlink_to(input_dir, target_is_directory=True)
        source_domains = source / "spatial_domains"
        if source_domains.is_dir() and not (target / "spatial_domains").exists():
            shutil.copytree(source_domains, target / "spatial_domains")
        staged.append({"sample_id": sample_id, "source": str(source), "target": str(target), "read_only_stage": True})

    out = BATCH_RESULTS / "config" / "staged_samples.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"source_results": str(SOURCE_RESULTS), "samples": staged}, indent=2) + "\n")
    print(json.dumps({"n_samples": len(staged), "excluded": sorted(EXCLUDED), "staged": staged}, indent=2))


if __name__ == "__main__":
    main()
