# Official Segmented Review Decisions

## Confirmed Findings

The review of the eight retained Visium HD samples confirms:

- The lower complexity of R4, R2, and R6 is consistent with the observed pathology region, tissue quality, and sequencing depth.
- Official Space Ranger cell boundaries provide reasonable coverage across glandular, stromal, and inflammation-dense areas.
- No separate minimum `counts`/`genes` thresholds are required for low-RNA immune cells.
- Patient, region, and clinical grouping metadata required for downstream differential analysis are complete.

## Analysis Decisions

- Retain all official filtered cells for the next stage.
- Do not remove R4, R2, or R6 on the basis of complexity alone.
- Use candidate QC fields for sensitivity analysis, not as an automatic deletion rule.
- Keep integer raw counts in `layers["counts"]` and use normalized/log-transformed `X` for modeling.
- Treat sample, patient, region, and clinical group as analysis metadata and biological replication units where appropriate.
- Use CPU LSF for the first-pass cell-level graph analysis. The installed GPU stack is not required for this stage; any later GPU acceleration must be recorded with its actual backend and fallback state.

## Next Stage

`submit_official_cell_analysis.sh` submits an eight-sample CPU array. It performs per-sample normalization, HVG selection, PCA, neighbors, Leiden clustering, UMAP, prostate marker scoring, and spatial cluster plots. Zero-count objects remain in the saved H5AD and are excluded only from the graph calculation.

The next review checkpoint is the per-sample cluster composition and marker-score report. Integration and spatial-domain modeling should begin only after that checkpoint is accepted.
