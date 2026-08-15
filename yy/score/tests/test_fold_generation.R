#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
project <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
source(file.path(project, "R", "00_common.R"))

eid <- sprintf("%07d", seq_len(30L))
baseline <- as.Date("2010-01-01")
event_date <- rep(as.Date(NA), 30L)
event_date[1:10] <- as.Date("2008-01-01")
event_date[11:20] <- as.Date("2012-01-01")
all_data <- data.table(
  eid = eid,
  date_attend = baseline,
  fod_icd10_cvd_cad = event_date,
  date_lost = as.Date(NA),
  date_death = as.Date(NA),
  ethnic.c = "White"
)
protein_data <- data.table(eid = eid, protein_1 = seq_along(eid))
config <- list(
  follow_end = "2024-01-01",
  expected = list(yin_n = 20L, yin_events = 10L, yang_n = 10L, outer_folds = 5L),
  restricted_reference = list(yin_md5 = "", yang_md5 = ""),
  fold_generation = list(
    ethnicity_column = "ethnic.c", ethnicity_pattern = "White",
    include_missing_ethnicity = TRUE, yin_seed = 2026L, yang_seed = 2027L,
    yang_duration_breaks = list(0, 1, 3, 5, 10, "Inf")
  )
)

first <- score_generate_fold_tables(copy(all_data), copy(protein_data), config)
second <- score_generate_fold_tables(copy(all_data), copy(protein_data), config)
score_validate_fold_tables(first, config)
stopifnot(
  identical(first$yin, second$yin),
  identical(first$yang, second$yang),
  identical(sort(unique(first$yin$fold)), 1:5),
  identical(sort(unique(first$yang$fold)), 1:5)
)

test_root <- tempfile("yy_fold_generation_test_")
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)
all_path <- file.path(test_root, "all.rds")
prot_path <- file.path(test_root, "prot.rds")
saveRDS(all_data, all_path)
saveRDS(protein_data, prot_path)
hash_root <- file.path(test_root, "hash")
hash_files <- score_write_fold_candidates(first, hash_root)
config$restricted_reference$yin_canonical_md5 <- score_canonical_fold_md5(hash_files$yin, "yin")
config$restricted_reference$yang_canonical_md5 <- score_canonical_fold_md5(hash_files$yang, "yang")
fold_root <- file.path(test_root, "folds")
paths <- list(
  all_rds = all_path, prot_rds = prot_path, fold_root = fold_root,
  fold_yin = file.path(fold_root, "fold_assignment_yin.csv"),
  fold_yang = file.path(fold_root, "fold_assignment_yang.csv")
)

score_ensure_fold_manifests(paths, config, copy(all_data), copy(protein_data), install = FALSE)
stopifnot(!file.exists(paths$fold_yin), !file.exists(paths$fold_yang))
score_ensure_fold_manifests(paths, config, copy(all_data), copy(protein_data), install = TRUE)
stopifnot(
  file.exists(paths$fold_yin), file.exists(paths$fold_yang),
  file.exists(file.path(fold_root, "fold_generation_manifest.csv"))
)
score_ensure_fold_manifests(paths, config, copy(all_data), copy(protein_data), install = TRUE)

partial_root <- file.path(test_root, "partial")
dir.create(partial_root)
partial_paths <- paths
partial_paths$fold_root <- partial_root
partial_paths$fold_yin <- file.path(partial_root, "fold_assignment_yin.csv")
partial_paths$fold_yang <- file.path(partial_root, "fold_assignment_yang.csv")
invisible(file.copy(paths$fold_yin, partial_paths$fold_yin))
partial_error <- tryCatch({
  score_ensure_fold_manifests(
    partial_paths, config, copy(all_data), copy(protein_data), install = TRUE
  )
  ""
}, error = function(error) conditionMessage(error))
stopifnot(grepl("Partial fold reference", partial_error, fixed = TRUE))
cat("test_fold_generation.R: PASS\n")
