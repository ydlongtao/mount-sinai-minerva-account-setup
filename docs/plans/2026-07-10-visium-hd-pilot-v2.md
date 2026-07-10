# Visium HD Pilot V2 Implementation Plan

> Review draft only. Do not upload or submit jobs until the user approves it.

## Goal

Prepare a second pilot workflow for the human prostate cancer Visium HD sample
\`SC000895-R4\`. The workflow must preserve full-gene expression, avoid unstable
global modeling of sparse 8 µm bins, validate graph connectivity before
clustering, and produce an auditable HTML report.

## Version Boundary

- Keep all v1 scripts and results unchanged.
- New scripts use a \`v2\` suffix or a dedicated \`v2\` directory.
- Remote result root: \`results/pilot_v2/SC000895-R4/\`.
- Remote checkpoint root: \`work_v2/SC000895-R4/\`.
- No raw-data writes.

## Proposed Files

- \`config/pilot_v2.yaml\`
- \`scripts/sp_v2_utils.py\`
- \`scripts/04_pilot_visium_hd_v2_qc.py\`
- \`scripts/04_pilot_visium_hd_v2_model.py\`
- \`scripts/09_review_pilot_v2.py\`
- \`lsf/04a_pilot_v2_qc.lsf\`
- \`lsf/04b_pilot_v2_gpu_model.lsf\`
- \`lsf/04c_pilot_v2_report.lsf\`
- focused tests under \`tests/\`

## Main Changes

### 1. Resolution Strategy

- Use the 16 µm matrix for global QC, PCA, neighbors, UMAP, and Leiden.
- Keep the 8 µm matrix for high-resolution expression display and local
  refinement, not global clustering.
- Map 16 µm labels to 8 µm bins with \`barcode_mappings.parquet\`.
- Defer segmented-output modeling until this 16/8 µm workflow passes.

Reason: the v1 8 µm pilot had median 12 counts and 11 genes per bin, producing
1,452 clusters, including 805 clusters with fewer than 100 bins.

### 2. Preserve Full-Gene Data

- Store original counts in \`layers["counts"]\`.
- Keep normalized, log-transformed full-gene data in the final AnnData object.
- Mark HVGs in \`var["highly_variable"]\`; do not replace the full object with
  an HVG-only slice.
- Create a temporary HVG-only model object for GPU calculations.
- Merge PCA, UMAP, graph metadata, and cluster labels back into the full object.

### 3. Separate QC and Modeling

- Stage 04a generates diagnostics and candidate QC thresholds only.
- Default QC mode is \`flag_only\`; candidate failures remain in the output.
- Stage 04b runs only after threshold review.
- QC report includes count, gene, mitochondrial fraction, and spatial maps.
- Do not run single-cell doublet detection on square Visium HD bins.
- Treat ambient RNA signatures as review signals, not automatic exclusions.

### 4. Tissue Registration Check

- Generate a barcode-coordinate and histology overlay before applying an image
  mask.
- Record the affine transform and image dimensions in the report.
- Apply a tissue-image mask only after visual confirmation of registration.

### 5. Graph Connectivity Guardrail

- Use RAPIDS CAGRA; do not use the previously unstable IVFFlat path.
- Test \`n_neighbors\` sequentially at 15, 30, and 50.
- For each graph, report connected-component count, largest-component fraction,
  and counts of components smaller than 10, 50, and 100 observations.
- Require the largest component to contain at least 99% of observations before
  UMAP and Leiden run.
- Stop with a non-zero exit status if no candidate graph passes.

### 6. Clustering Grid

- Run Leiden resolutions 0.02, 0.05, 0.1, 0.2, and 0.4.
- Preserve every label column; do not automatically declare one final.
- Report cluster count, median cluster size, minimum cluster size, and fractions
  in small clusters for every resolution.
- Use GPU Leiden when the preflight succeeds.
- CPU fallback uses Scanpy's igraph implementation rather than the slow direct
  leidenalg path used by v1.

### 7. GPU Environment Validation

- Keep RAPIDS core dependencies internally compatible, especially Dask and
  Distributed.
- GPU preflight must execute a synthetic CAGRA graph and GPU Leiden, not only
  import packages.
- Keep CPU plotting/report dependencies separate if Squidpy or SpatialData pins
  conflict with the RAPIDS stack.
- Record package versions and GPU model in \`run_state.json\`.

### 8. Marker Panel

- Retain all existing prostate cancer markers.
- Add candidates for review: \`MSMB\`, \`AZGP1\`, \`KRT15\`, \`LST1\`,
  \`TYROBP\`, \`COL3A1\`, \`PDGFRA\`, \`MYL9\`, \`DES\`, \`EMCN\`, \`CHGA\`,
  \`CHGB\`, \`SYP\`, and \`INSM1\`.
- Plot markers from the full normalized object.
- V2 produces coarse compartment scores and candidate labels only; it does not
  claim final cell-type annotation.

### 9. Plot and HTML Report Changes

- Replace full-data jitter violins with histograms, density plots, and hexbin
  plots.
- Deterministically subsample at most 100,000 observations for UMAP rendering;
  retain all coordinates in H5AD.
- Limit plot legends to the 20 largest clusters plus \`Other\`.
- Rasterize dense spatial layers and include histology overlays.
- Create a standalone HTML report with explicit pass/fail gates, graph
  diagnostics, QC thresholds, resolution comparison, marker panels, and output
  checksums.

### 10. Metrics Parser Repair

- Update \`03_collect_qc_metrics.py\` to parse Space Ranger's two-row wide
  \`metrics_summary.csv\` format.
- Extract both 8 µm and 16 µm bin counts, genes/UMIs per bin, sequencing
  saturation, reads, and segmented-cell metrics when available.
- Report missing fields instead of silently returning incorrect values.

### 11. Checkpoint and Resume

- Checkpoints: \`raw_qc\`, \`normalized_full\`, \`hvg_pca\`, \`graph\`, and
  \`embedding_labels\`.
- Write checkpoints atomically through a temporary file and rename.
- Store input hashes, parameter hashes, package versions, and stage status in
  \`run_state.json\`.
- Resume only when input and parameter hashes match.

## Proposed LSF Jobs

- \`04a_pilot_v2_qc.lsf\`: express, CPU, 2-4 cores, 8 GB/core, 1 hour.
- \`04b_pilot_v2_gpu_model.lsf\`: gpuexpress, 1 GPU, 4 cores, 16 GB/core,
  4 hours.
- \`04c_pilot_v2_report.lsf\`: express, CPU, 4 cores, 12 GB/core, 2 hours.
- Every job uses \`#BSUB -P acc_DiseaseGeneCell\` and
  \`#BSUB -R "span[hosts=1]"\`.

## Acceptance Tests

- Final 16 µm and mapped 8 µm H5AD files reopen successfully.
- Full gene set and \`layers["counts"]\` are retained.
- All configured marker genes present in the source matrix remain available.
- Largest graph component is at least 99%.
- No plotting legend contains more than 21 displayed categories.
- Metrics parser matches known R4 values from the Space Ranger CSV.
- HTML report links to every generated table and figure.
- Any failed quality gate causes the relevant LSF job to fail visibly.

## Decisions Required Before Implementation

1. Approve 16 µm global modeling with label mapping to 8 µm.
2. Approve separate QC review and GPU modeling jobs.
3. Approve deferring segmented-output modeling.
4. Approve the expanded marker list.
5. Approve a separate RAPIDS core environment when dependency pins conflict.
6. Confirm whether a pathologist-provided region annotation or tissue mask is
   available for this sample.
