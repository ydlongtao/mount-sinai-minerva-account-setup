# AGENTS.md

This repository documents a working Mount Sinai Minerva account setup for
single-cell transcriptomics analysis.

## Scope

Maintain this repository as an operational setup guide, not as a project data
repository. Do not add:

- passwords
- MFA screenshots or prompts
- API tokens
- SSH private keys
- protected patient or project data
- private server-side analysis paths

## Target Platform

- Mount Sinai Minerva HPC
- user-space installation under `~/miniconda3`
- bash shell on the remote server
- conda/mamba environments:
  - `base`
  - `omicverse`
  - `scop_r`
  - `scop_env`

## Editing Guidelines

- Keep commands copy-paste runnable.
- Prefer user-space commands over root-level changes.
- Separate Mac-side commands from server-side commands when both are involved.
- Record observed package versions when they are verified.
- Include troubleshooting notes for errors that were actually encountered.
- Keep public documentation free of secrets and account-specific credentials.

## Validation Checklist

Before publishing updates, verify the relevant environment:

```bash
source ~/.bashrc
conda env list
```

For OmicVerse:

```bash
conda activate omicverse
python - <<'PY'
import torch
import torch_geometric
import omicverse as ov
print(torch.__version__)
print(torch_geometric.__version__)
print(ov.__version__)
PY
```

For scop R:

```bash
conda activate scop_r
Rscript - <<'RSCRIPT'
library(scop)
cat("scop=", as.character(packageVersion("scop")), "\n")
cat("Seurat=", as.character(packageVersion("Seurat")), "\n")
cat("Signac=", as.character(packageVersion("Signac")), "\n")
RSCRIPT
```

For the Python environment created by scop:

```bash
conda activate scop_env
python - <<'PY'
import scanpy
import scvelo
import anndata
print(scanpy.__version__)
print(scvelo.__version__)
print(anndata.__version__)
PY
```

## Known Decisions

- Use `conda-forge` plus `nodefaults` to avoid non-interactive Anaconda defaults
  channel Terms of Service failures.
- Use `CONDA_SOLVER=classic` for large R environment creation if `mamba` fails
  with transaction resource errors.
- Install `scop` from a local GitHub clone with `R CMD INSTALL` when `pak` hides
  build details.
- Keep OmicVerse and scop in separate environments to reduce dependency
  conflicts.
