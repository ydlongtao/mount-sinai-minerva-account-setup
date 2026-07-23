#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: 65_run_infercnv_i3.R <input_dir> <output_dir>")
input_dir <- normalizePath(args[[1]])
output_dir <- normalizePath(args[[2]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(infercnv))

counts <- file.path(input_dir, "counts.mtx.gz")
annotations <- file.path(input_dir, "annotations.tsv.gz")
gene_order <- file.path(input_dir, "gene_order.tsv.gz")

obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file = annotations,
  delim = "\t",
  gene_order_file = gene_order,
  ref_group_names = c("main_reference")
)

res <- infercnv::run(
  obj,
  cutoff = 0.1,
  out_dir = output_dir,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = TRUE,
  HMM_type = "i3",
  analysis_mode = "samples",
  num_threads = as.integer(Sys.getenv("LSB_DJOB_NUMPROC", "8")),
  output_format = "pdf",
  plot_steps = FALSE,
  resume_mode = FALSE
)

summary <- list(
  status = "pass",
  hmm_type = "i3",
  reference_group = "main_reference",
  tumor_group = "malignant_epithelial_candidate",
  input_dir = input_dir,
  output_dir = output_dir,
  infercnv_version = as.character(packageVersion("infercnv"))
)
jsonlite::write_json(summary, file.path(output_dir, "i3_run_summary.json"), pretty = TRUE, auto_unbox = TRUE)
cat(jsonlite::toJSON(summary, pretty = TRUE, auto_unbox = TRUE), "\n")
