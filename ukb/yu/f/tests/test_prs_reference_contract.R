#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else normalizePath(".")
project_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

sources <- fread(file.path(project_dir, "config", "prs_gwas_sources.tsv"))
policy <- fread(file.path(project_dir, "config", "prs_weight_source_policy.tsv"))
thresholds <- fread(file.path(project_dir, "config", "prs_thresholds.tsv"))

stopifnot(
  nrow(sources) == 13L,
  uniqueN(sources$outcome_id) == 13L,
  nrow(policy) == 13L,
  uniqueN(policy$outcome_id) == 13L,
  setequal(sources$outcome_id, policy$outcome_id),
  nrow(thresholds) == 5L,
  uniqueN(thresholds$threshold) == 5L,
  uniqueN(thresholds$score_column) == 5L,
  all(policy$author_weight_file_public %in% c(FALSE, "FALSE"))
)

hf <- policy[outcome_id == "heart_failure"]
non_hf <- policy[outcome_id != "heart_failure"]
stopifnot(
  nrow(hf) == 1L,
  grepl("not reported|provisional", paste(hf$weight_basis, hf$interpretation), ignore.case = TRUE),
  nrow(non_hf) == 12L,
  all(grepl("S26", non_hf$weight_basis, fixed = TRUE))
)

cat("PASS: 13 outcomes, 5 thresholds, 12 S26 GWAS reference-weight sources, and provisional HF source.\n")
