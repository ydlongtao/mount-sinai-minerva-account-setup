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

## Visium HD Archive Rules

This repository also contains the reproducible Visium HD prostate cancer pilot
under sp_project_local/. Maintain code and operational documentation, not data.

- Keep v1 and v2 scripts versioned; do not overwrite a validated workflow
  silently.
- Use 16 um for global pilot modeling and 8 um for mapped high-resolution
  review unless the analysis plan is explicitly revised.
- Keep #BSUB -P acc_DiseaseGeneCell in every production LSF script.
- Keep raw data read-only and put temporary files under
  /sc/arion/scratch/huangl21/sp_project/.
- Treat H5AD structure checks, graph connectivity, and report generation as
  required validation steps.
- Do not commit local_results/, results/, logs/, work/, H5AD, BAM/CRAM, SVS,
  or credentials.

The successful v2 GPU pilot ran on an H100 node with CUDA 13 RAPIDS. GPU
Leiden may fall back to Scanpy igraph when the GPU Dask stack is incompatible;
this behavior must be recorded in run_state.json.
