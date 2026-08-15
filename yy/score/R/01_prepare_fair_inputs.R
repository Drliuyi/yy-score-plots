#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
project_dir <- dirname(dirname(gsub("\\\\", "/", script_file)))
source(file.path(project_dir, "R", "00_common.R"))
parsed <- score_args()
paths <- score_resolve(project_dir, parsed)
cfg <- jsonlite::fromJSON(paths$config, simplifyVector = TRUE)
score_print_roots(paths)

if (parsed$stage == "help") {
  cat("Usage: Rscript 01_prepare_fair_inputs.R --stage=preflight|prepare|status [--dir0=/mnt/d ...]\n")
  quit(status = 0L)
}

required <- c(paths$all_rds, paths$prot_rds, paths$config)
if (parsed$stage == "status") {
  marker <- file.path(paths$common_root, "COMPLETE")
  cat(if (file.exists(marker)) "COMPLETE\n" else "NOT_COMPLETE\n")
  if (file.exists(marker)) cat("  ", paths$common_root, "\n", sep = "")
  quit(status = if (file.exists(marker)) 0L else 1L)
}
if (!parsed$stage %in% c("preflight", "prepare")) stop("--stage must be preflight, prepare, or status.", call. = FALSE)
score_require(required)

all_data <- data.table::as.data.table(readRDS(paths$all_rds))
protein_data <- data.table::as.data.table(readRDS(paths$prot_rds))
if (anyDuplicated(as.character(all_data$eid)) || anyDuplicated(as.character(protein_data$eid))) {
  stop("all.rds/prot.rds eid must be unique.", call. = FALSE)
}
proteins <- setdiff(names(protein_data), "eid")
if (length(proteins) != cfg$expected$protein_n || anyDuplicated(proteins)) {
  stop("Raw prot.rds protein contract failed: n=", length(proteins), call. = FALSE)
}
folds <- score_ensure_fold_manifests(
  paths, cfg, all_data, protein_data, install = identical(parsed$stage, "prepare")
)
fold_yin <- folds$yin
fold_yang <- folds$yang
endpoint <- score_build_endpoint(all_data, cfg$follow_end)

all_index_yin <- match(fold_yin$eid, as.character(all_data$eid))
all_index_yang <- match(fold_yang$eid, as.character(all_data$eid))
prot_index_yin <- match(fold_yin$eid, as.character(protein_data$eid))
prot_index_yang <- match(fold_yang$eid, as.character(protein_data$eid))
endpoint_index_yin <- match(fold_yin$eid, endpoint$eid)
endpoint_index_yang <- match(fold_yang$eid, endpoint$eid)
if (anyNA(c(all_index_yin, all_index_yang, prot_index_yin, prot_index_yang,
            endpoint_index_yin, endpoint_index_yang))) {
  stop("Raw all.rds/prot.rds do not cover every locked Yin/Yang EID.", call. = FALSE)
}
yin_endpoint <- endpoint[endpoint_index_yin]
yang_endpoint <- endpoint[endpoint_index_yang]
if (!identical(as.integer(yin_endpoint$event), as.integer(fold_yin$event)) ||
    any(yin_endpoint$prevalent != 0L) || any(!is.finite(yin_endpoint$time) | yin_endpoint$time <= 0) ||
    any(yang_endpoint$prevalent != 1L) || any(!is.finite(yang_endpoint$b2e) | yang_endpoint$b2e >= 0)) {
  stop("Recomputed endpoint disagrees with the locked fair-cohort contract.", call. = FALSE)
}
if (nrow(fold_yin) != cfg$expected$yin_n || sum(fold_yin$event) != cfg$expected$yin_events ||
    nrow(fold_yang) != cfg$expected$yang_n) {
  stop("Locked fair-cohort counts disagree with config.", call. = FALSE)
}
required_clinical <- c("age", "sex")
if (length(setdiff(required_clinical, names(all_data)))) {
  stop("all.rds lacks age/sex for fair score trajectories.", call. = FALSE)
}
yin <- data.table(
  eid = fold_yin$eid, time = yin_endpoint$time, event = yin_endpoint$event,
  outer_fold = as.integer(fold_yin$fold),
  age = as.numeric(all_data$age[all_index_yin]), sex = as.integer(all_data$sex[all_index_yin])
)
yang <- data.table(
  eid = fold_yang$eid, yang_fold = as.integer(fold_yang$fold),
  disease_duration_years = -yang_endpoint$b2e,
  age = as.numeric(all_data$age[all_index_yang]), sex = as.integer(all_data$sex[all_index_yang])
)
if (any(!is.finite(yin$age)) || any(!yin$sex %in% c(0L, 1L)) ||
    any(!is.finite(yang$age)) || any(!yang$sex %in% c(0L, 1L))) {
  stop("Age/sex contract failed in the locked fair cohort.", call. = FALSE)
}

cat("PREPARE_QC_PASS yin=", nrow(yin), " events=", sum(yin$event),
    " yang=", nrow(yang), " proteins=", length(proteins), "\n", sep = "")
if (parsed$stage == "preflight") quit(status = 0L)

marker <- file.path(paths$common_root, "COMPLETE")
if (file.exists(marker)) {
  expected_paths <- file.path(paths$common_root, c(
    "participants_yin.csv.gz", "participants_yang.csv.gz", "protein_features.csv",
    "protein_yin.f32", "protein_yang.f32", "locked_yin_target.rds",
    "locked_yang_target.rds", "input_manifest.csv"
  ))
  score_require(expected_paths, "completed fair input")
  cat("Existing complete fair input preserved: ", paths$common_root, "\n", sep = "")
  quit(status = 0L)
}
if (dir.exists(paths$common_root)) {
  existing <- list.files(paths$common_root, all.files = TRUE, no.. = TRUE)
  if (length(existing)) stop("Incomplete fair-input directory exists; inspect it: ", paths$common_root, call. = FALSE)
}
dir.create(paths$common_root, recursive = TRUE, showWarnings = FALSE)

protein_matrix_yin <- as.matrix(protein_data[prot_index_yin, ..proteins])
protein_matrix_yang <- as.matrix(protein_data[prot_index_yang, ..proteins])
storage.mode(protein_matrix_yin) <- "double"
storage.mode(protein_matrix_yang) <- "double"
score_atomic_csv(yin, file.path(paths$common_root, "participants_yin.csv.gz"))
score_atomic_csv(yang, file.path(paths$common_root, "participants_yang.csv.gz"))
score_atomic_rds(
  yin[, .(eid, time_years = time, event)],
  file.path(paths$common_root, "locked_yin_target.rds")
)
score_atomic_rds(
  yang[, .(eid, disease_duration_years)],
  file.path(paths$common_root, "locked_yang_target.rds")
)
score_atomic_csv(data.table(feature_index = seq_along(proteins) - 1L, feature = proteins),
                 file.path(paths$common_root, "protein_features.csv"))
score_write_f32(protein_matrix_yin, file.path(paths$common_root, "protein_yin.f32"))
score_write_f32(protein_matrix_yang, file.path(paths$common_root, "protein_yang.f32"))

manifest_inputs <- c(required, paths$fold_yin, paths$fold_yang)
manifest <- data.table(
  role = c("all_rds", "prot_rds", "config", "fold_yin", "fold_yang"),
  path = normalizePath(manifest_inputs, winslash = "/", mustWork = TRUE),
  bytes = as.numeric(file.info(manifest_inputs)$size),
  md5 = unname(tools::md5sum(manifest_inputs))
)
score_atomic_csv(manifest, file.path(paths$common_root, "input_manifest.csv"))
score_atomic_text(c(
  "status=COMPLETE", paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  paste0("yin_n=", nrow(yin)), paste0("yin_events=", sum(yin$event)),
  paste0("yang_n=", nrow(yang)), paste0("protein_n=", length(proteins)),
  paste0("follow_end=", cfg$follow_end)
), marker)
cat("FAIR_INPUT_COMPLETE\nOUTPUT_ROOT=", paths$common_root, "\n", sep = "")
