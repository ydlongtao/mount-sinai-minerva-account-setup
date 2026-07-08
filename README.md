# Mount Sinai Minerva Account Setup

Practical setup notes for a Mount Sinai Minerva HPC user account focused on
single-cell transcriptomics analysis.

This guide records a working user-space setup for:

- Miniconda3 and `mamba`
- a Python OmicVerse environment
- an R-based `scop` environment
- the Python environment created by `scop::PrepareEnv()`

It is intentionally account-safe: no passwords, MFA details, tokens, SSH keys, or
private project paths are included.

## Tested Context

- HPC: Mount Sinai Minerva
- Login node observed during setup: `li04e04`
- User-space installation prefix: `~/miniconda3`
- Shell: bash on the remote server
- Package strategy: conda-forge first, with `nodefaults` to avoid Anaconda
  defaults channel Terms of Service prompts in non-interactive installs

## Environment Summary

| Environment | Purpose | Path |
| --- | --- | --- |
| `base` | Miniconda and mamba bootstrap | `~/miniconda3` |
| `omicverse` | Python OmicVerse analysis | `~/miniconda3/envs/omicverse` |
| `scop_r` | R/Seurat/Signac/scop analysis | `~/miniconda3/envs/scop_r` |
| `scop_env` | Python dependencies managed by `scop::PrepareEnv()` | `~/miniconda3/envs/scop_env` |

## 1. Install Miniconda3 and mamba

Install Miniconda3 under the user home directory:

```bash
cd "$HOME"
curl -L --fail --retry 3 \
  -o miniconda3-installer.sh \
  https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash miniconda3-installer.sh -b -p "$HOME/miniconda3"
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate base
```

Configure conda to prefer conda-forge and avoid the defaults channel:

```bash
conda config --set auto_activate_base true
conda config --set channel_priority strict
conda config --add channels nodefaults
conda config --add channels conda-forge
```

If the installer-level config still includes `defaults`, remove it there too:

```bash
conda config --file "$HOME/miniconda3/.condarc" --remove channels defaults || true
conda config --file "$HOME/miniconda3/.condarc" --add channels nodefaults
conda config --file "$HOME/miniconda3/.condarc" --add channels conda-forge
conda config --file "$HOME/miniconda3/.condarc" --set channel_priority strict
```

Install `mamba`:

```bash
conda install -n base --override-channels -c conda-forge -y mamba
conda init bash
source ~/.bashrc
```

Verify:

```bash
conda --version
mamba --version
conda config --show channels channel_priority
```

Expected working state:

```text
channels:
  - conda-forge
  - nodefaults
channel_priority: strict
```

## 2. Install OmicVerse

The conda-only installation path can be heavy on a login node. A lighter and
working approach is to create a Python environment with conda, then install
PyTorch and OmicVerse with `uv pip`.

```bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate base

mamba create -n omicverse --override-channels -c conda-forge -y \
  python=3.10 pip uv ipykernel jupyterlab

conda activate omicverse

uv pip install --index-url https://download.pytorch.org/whl/cpu \
  torch torchvision torchaudio

uv pip install torch_geometric
uv pip install omicverse

python -m ipykernel install --user \
  --name omicverse \
  --display-name "Python (omicverse)"
```

Verify:

```bash
python - <<'PY'
import sys
import torch
import torch_geometric
import omicverse as ov

print("python", sys.version)
print("torch", torch.__version__, "cuda_available", torch.cuda.is_available())
print("torch_geometric", torch_geometric.__version__)
print("omicverse", ov.__version__)
PY
```

Observed working versions:

```text
Python 3.10.20
torch 2.12.1+cpu
torch_geometric 2.8.0
omicverse 2.2.3
```

## 3. Install R, Seurat, Signac, and scop

`scop` is an R package for spatial and single-cell omics analysis. The package
requires R >= 4.1 and is installed from GitHub.

Create the R environment. On Minerva login nodes, the classic conda solver was
more stable than libmamba for this large R transaction.

```bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate base

CONDA_SOLVER=classic conda create -n scop_r \
  --override-channels -c conda-forge -c bioconda -y \
  r-base=4.4 r-pak r-remotes r-devtools r-irkernel r-reticulate \
  r-seurat r-seuratobject r-signac r-uwot \
  r-rcpp r-rcppparallel r-rcppannoy r-rcppeigen r-rcpparmadillo \
  r-dplyr r-ggplot2 r-ggrepel r-ggforce r-ggnewscale \
  r-igraph r-matrix r-rlang r-gtable r-cli \
  bioconductor-complexheatmap \
  bioconductor-genomicranges \
  bioconductor-s4vectors \
  bioconductor-summarizedexperiment
```

Install GitHub dependencies and `scop`:

```bash
conda activate scop_r

Rscript - <<'RSCRIPT'
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 1200)
Sys.setenv(MAKEFLAGS = "-j2")

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(
  "mengxu98/thisutils",
  upgrade = "never",
  dependencies = c("Depends", "Imports", "LinkingTo")
)

remotes::install_github(
  "mengxu98/thisplot",
  upgrade = "never",
  dependencies = c("Depends", "Imports", "LinkingTo")
)
RSCRIPT

mkdir -p "$HOME/src"
git clone https://github.com/mengxu98/scop.git "$HOME/src/scop"
R CMD INSTALL --no-multiarch --with-keep.source "$HOME/src/scop"
```

Register an R Jupyter kernel:

```bash
Rscript -e 'IRkernel::installspec(name = "scop_r", displayname = "R (scop_r)", user = TRUE)'
```

Verify:

```bash
Rscript - <<'RSCRIPT'
library(scop)
cat("scop=", as.character(packageVersion("scop")), "\n")
cat("Seurat=", as.character(packageVersion("Seurat")), "\n")
cat("Signac=", as.character(packageVersion("Signac")), "\n")
cat("reticulate=", as.character(packageVersion("reticulate")), "\n")
RSCRIPT
```

Observed working versions:

```text
scop 0.8.9
Seurat 5.5.1
Signac 1.16.0
reticulate 1.46.0
```

## 4. Prepare scop Python Environment

Some `scop` functions, such as PAGA and scVelo workflows, require a separate
Python environment. `scop::PrepareEnv()` creates this default environment as
`scop_env`.

```bash
conda activate scop_r

Rscript - <<'RSCRIPT'
options(scop_env_init = TRUE)
options(reticulate.conda_binary = file.path(Sys.getenv("HOME"), "miniconda3", "bin", "mamba"))
options(timeout = 1200)

library(scop)
scop::PrepareEnv(conda = file.path(Sys.getenv("HOME"), "miniconda3", "bin", "mamba"))
RSCRIPT
```

Verify the Python side:

```bash
conda activate scop_env

python - <<'PY'
import scanpy
import scvelo
import anndata
import pandas
import numpy

print("scanpy", scanpy.__version__)
print("scvelo", scvelo.__version__)
print("anndata", anndata.__version__)
print("pandas", pandas.__version__)
print("numpy", numpy.__version__)
PY
```

Observed working versions:

```text
scanpy 1.11.3
scvelo 0.3.3
anndata 0.11.4
pandas 2.2.0
numpy 1.26.4
```

## 5. Daily Use

Use `scop` in R:

```bash
source ~/.bashrc
conda activate scop_r
R
```

```r
library(scop)
```

Use the Python environment created for `scop`:

```bash
conda activate scop_env
python
```

Use OmicVerse:

```bash
conda activate omicverse
python
```

List available Jupyter kernels:

```bash
jupyter kernelspec list
```

## 6. Troubleshooting Notes

### Anaconda defaults Terms of Service

If conda reports that `https://repo.anaconda.com/pkgs/main` or
`https://repo.anaconda.com/pkgs/r` Terms of Service have not been accepted, avoid
defaults by using:

```bash
conda config --show-sources
conda config --remove channels defaults || true
conda config --add channels nodefaults
conda config --add channels conda-forge
```

Also check installer-level config:

```bash
cat "$HOME/miniconda3/.condarc"
```

### libmamba transaction resource errors

Large conda transactions on login nodes may fail with errors such as:

```text
Resource temporarily unavailable
std::bad_alloc
```

Mitigations:

- reduce `mamba` concurrency:

```bash
MAMBA_EXTRACT_THREADS=1 MAMBA_DOWNLOAD_THREADS=1 mamba install ...
```

- for large R environments, use conda classic:

```bash
CONDA_SOLVER=classic conda create ...
```

- split installation into smaller stages.

### scop PrepareEnv bedtools channel issue

During `scop::PrepareEnv()`, the conda substep may report that `bedtools` is not
available in the selected channel. In the observed setup, `PrepareEnv()` then
continued with `uv pip` and successfully installed the Python packages.

## References

- Mount Sinai Minerva HPC documentation: https://labs.icahn.mssm.edu/minervalab/
- OmicVerse installation guide: https://omicverse.readthedocs.io/en/latest/Installation_guild.html
- scop GitHub repository: https://github.com/mengxu98/scop
- scop documentation: https://mengxu98.github.io/scop/
