cat > /tmp/make_spatial_structure_report.py <<'PY'
from pathlib import Path
from datetime import datetime
import os
import socket


def human(n):
    n = float(n)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if n < 1024 or unit == "TB":
            return f"{int(n)} B" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024


def find_root():
    candidates = [
        Path.cwd(),
        Path.home() / "DiseaseGeneCell" / "Huang_lab_project" / "SpatialTranscriptome",
        Path("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/SpatialTranscriptome"),
        Path.home() / "DiseaseGeneCell" / "Huang_lab_data" / "SpatialTranscriptome",
        Path("/sc/arion/projects/DiseaseGeneCell/Huang_lab_data/SpatialTranscriptome"),
    ]
    for path in candidates:
        if path.exists() and path.is_dir() and (
            path.name == "SpatialTranscriptome" or (path / "wangy33.u.hpc.mssm.edu").exists()
        ):
            return path.resolve()

    base = Path.home() / "DiseaseGeneCell"
    if base.exists():
        hits = sorted(base.glob("**/SpatialTranscriptome"))
        if hits:
            return hits[0].resolve()
    return Path.cwd().resolve()


root = find_root()
all_dirs = []
all_files = []

for dirpath, dirnames, filenames in os.walk(root):
    p = Path(dirpath)
    dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
    all_dirs.append(p)
    for fn in sorted(filenames):
        if fn.startswith("."):
            continue
        all_files.append(p / fn)

lines = []
lines.append("# SpatialTranscriptome Directory Structure")
lines.append("")
lines.append(f"- Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
lines.append(f"- Host: {socket.gethostname()}")
lines.append(f"- User: {os.environ.get('USER', 'unknown')}")
lines.append(f"- Root: `{root}`")
lines.append("")
lines.append("## Summary")
lines.append("")
lines.append(f"- Directories: {len(all_dirs)}")
lines.append(f"- Files: {len(all_files)}")
try:
    total_size = sum(f.stat().st_size for f in all_files if f.exists())
    lines.append(f"- Total visible file size under root: {human(total_size)}")
except Exception as e:
    lines.append(f"- Total visible file size under root: unavailable ({e})")
lines.append("")

lines.append("## Directory Tree (directories, max depth 5)")
lines.append("")
for p in sorted(all_dirs, key=lambda x: str(x.relative_to(root))):
    rel = p.relative_to(root)
    depth = 0 if str(rel) == "." else len(rel.parts)
    if depth > 5:
        continue
    indent = "  " * depth
    name = root.name if depth == 0 else p.name
    lines.append(f"{indent}- `{name}/`")
lines.append("")

key_ext = {
    ".h5", ".html", ".csv", ".cloupe", ".bam", ".bai", ".parquet",
    ".tgz", ".gz", ".json", ".mtx", ".tsv", ".txt", ".h5ad", ".rds"
}
key_names = {
    "metrics_summary.csv", "web_summary.html", "feature_slice.h5",
    "molecule_info.h5", "barcode_mappings.parquet", "probe_set.csv",
    "possorted_genome_bam.bam", "possorted_genome_bam.bam.bai",
}

lines.append("## Key Files")
lines.append("")
for f in sorted(all_files, key=lambda x: str(x.relative_to(root))):
    suffix_combo = "".join(f.suffixes[-2:]).lower()
    if f.name in key_names or f.suffix.lower() in key_ext or suffix_combo in {".tsv.gz", ".mtx.gz"}:
        try:
            size = human(f.stat().st_size)
        except Exception:
            size = "unavailable"
        lines.append(f"- `{f.relative_to(root)}` ({size})")
lines.append("")

lines.append("## 10X Outs Directories")
lines.append("")
for outs in sorted(root.glob("**/outs")):
    lines.append(f"### `{outs.relative_to(root)}`")
    for child in sorted(outs.iterdir(), key=lambda x: x.name):
        try:
            if child.is_dir():
                n_dirs = sum(1 for item in child.rglob("*") if item.is_dir())
                n_files = sum(1 for item in child.rglob("*") if item.is_file())
                lines.append(f"- `{child.name}/` ({n_dirs} subdirs, {n_files} files)")
            elif child.is_file():
                lines.append(f"- `{child.name}` ({human(child.stat().st_size)})")
        except Exception as e:
            lines.append(f"- `{child.name}` (unavailable: {e})")
    lines.append("")

out = Path.home() / f"SpatialTranscriptome_structure_{datetime.now().strftime('%Y%m%d_%H%M%S')}.md"
out.write_text("\n".join(lines) + "\n")
print(f"REPORT_PATH={out}")
print("BEGIN_SPATIAL_REPORT")
print(out.read_text())
print("END_SPATIAL_REPORT")
PY
python /tmp/make_spatial_structure_report.py
