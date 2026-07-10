#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--source-image", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    modules = ["cellpose", "tifffile", "zarr", "skimage", "shapely"]
    status = {}
    for name in modules:
        try:
            module = importlib.import_module(name)
            status[name] = {"ok": True, "version": getattr(module, "__version__", "unknown")}
        except Exception as exc:
            status[name] = {"ok": False, "error": repr(exc)}

    paths = {
        "bin_dir": {"path": str(args.bin_dir), "exists": args.bin_dir.is_dir()},
        "source_image": {"path": str(args.source_image), "exists": args.source_image.is_file()},
        "output_dir": {"path": str(args.output_dir), "writable_parent": args.output_dir.parent.exists()},
    }
    result = {"modules": status, "paths": paths}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "cellpose_preflight.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    if not all(item["ok"] for item in status.values()) or not all(
        item.get("exists", item.get("writable_parent", False)) for item in paths.values()
    ):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
