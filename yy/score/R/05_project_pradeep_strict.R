#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
project_dir <- dirname(dirname(gsub("\\\\", "/", script_file)))
source(file.path(project_dir, "R", "00_common.R"))
parsed <- score_args()
paths <- score_resolve(project_dir, parsed, "pradeep-strict")
source_root <- score_norm(parsed$source_root %||% score_env(
  "PRADEEP_STRICT_MODEL_ROOT", file.path(paths$analysis_root, "ukb", "pradeep_strict")
))
score_print_roots(paths)
cat("  STRICT_SOURCE_ROOT=", source_root, "\n", sep = "")

find_participant <- function(side) {
  candidates <- file.path(paths$common_root, paste0("participants_", side, c(".csv.gz", ".csv")))
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop("Common participant file missing: ", paste(candidates, collapse = "; "), call. = FALSE)
  hit[[1L]]
}

base_file <- file.path(source_root, "outputs", "ukbppp_cardiac_analysis_base.rds")
coefficient_file <- file.path(source_root, "outputs", "lasso", "cad_coefficients.csv")
prediction_file <- file.path(source_root, "outputs", "lasso", "cad_predictions.csv")
common_files <- c(
  yin = find_participant("yin"), yang = find_participant("yang"),
  features = file.path(paths$common_root, "protein_features.csv"),
  protein_yin = file.path(paths$common_root, "protein_yin.f32"),
  protein_yang = file.path(paths$common_root, "protein_yang.f32")
)
required <- c(base_file, coefficient_file, prediction_file, common_files)

if (parsed$stage == "status") {
  cat(if (file.exists(file.path(paths$output_root, "COMPLETE"))) "COMPLETE\n" else "NOT_COMPLETE\n")
  quit(status = if (file.exists(file.path(paths$output_root, "COMPLETE"))) 0L else 1L)
}
if (!parsed$stage %in% c("preflight", "project")) {
  stop("Use --stage=preflight|project|status.", call. = FALSE)
}
score_require(required, "Pradeep strict projection input")

yin <- fread(common_files[["yin"]], colClasses = list(character = "eid"))
yang <- fread(common_files[["yang"]], colClasses = list(character = "eid"))
feature_table <- fread(common_files[["features"]])
features <- as.character(feature_table$feature)
if (nrow(yin) != 37127L || nrow(yang) != 1766L || length(features) != 2910L ||
    anyDuplicated(yin$eid) || anyDuplicated(yang$eid) || anyDuplicated(features)) {
  stop("Common Yin/Yang projection contract failed.", call. = FALSE)
}

coefficients_all <- fread(coefficient_file)
coefficients <- coefficients_all[model == "Proteins" & outc == "cad",
  .(feature = as.character(feature), coefficient = as.numeric(coefficient))]
intercept <- coefficients[feature == "(Intercept)", coefficient]
selected <- coefficients[feature != "(Intercept)" & coefficient != 0]
if (length(intercept) != 1L || !nrow(selected) || anyDuplicated(selected$feature) ||
    any(!is.finite(selected$coefficient))) {
  stop("Pradeep strict coefficient contract failed.", call. = FALSE)
}

key <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))
map_columns <- function(required_names, available_names, label) {
  available_key <- key(available_names)
  duplicated_key <- unique(available_key[duplicated(available_key)])
  if (length(intersect(key(required_names), duplicated_key))) {
    stop(label, " has ambiguous normalized feature names.", call. = FALSE)
  }
  index <- match(key(required_names), available_key)
  if (anyNA(index)) stop(label, " is missing: ", paste(required_names[is.na(index)], collapse = ", "), call. = FALSE)
  data.table(index = index, name = available_names[index])
}

base <- readRDS(base_file)
if (!all(c("dat", "protein_cols") %in% names(base))) stop("Pradeep strict base RDS contract failed.", call. = FALSE)
base_dat <- as.data.frame(base$dat)
if (!"eid" %in% names(base_dat) || anyDuplicated(as.character(base_dat$eid))) {
  stop("Pradeep strict base EID contract failed.", call. = FALSE)
}
base_map <- map_columns(selected$feature, intersect(as.character(base$protein_cols), names(base_dat)), "Pradeep base")
common_map <- map_columns(selected$feature, features, "Common 2,910-protein matrix")

## Reproduce the frozen Pradeep shuffle/split and derive projection imputation
## only from its native training partition.
set.seed(1234)
base_dat <- base_dat[sample(seq_len(nrow(base_dat))), , drop = FALSE]
base_dat$.strict_split_id <- seq_len(nrow(base_dat))
set.seed(1)
train_ids <- sample(base_dat$.strict_split_id, floor(0.80 * nrow(base_dat)))
eligible <- !is.na(base_dat$cad_inc) & is.finite(base_dat$cad_fu) & base_dat$cad_fu > 0
strict_dat <- base_dat[eligible, , drop = FALSE]
train_index <- which(strict_dat$.strict_split_id %in% train_ids)
training_median <- vapply(base_map$name, function(column) {
  value <- suppressWarnings(as.numeric(strict_dat[[column]][train_index]))
  result <- stats::median(value[is.finite(value)])
  if (!is.finite(result)) stop("No finite strict-training value for ", column, call. = FALSE)
  result
}, numeric(1L))

cat("PRADEEP_STRICT_PROJECTION_PREFLIGHT_PASS selected=", nrow(selected),
    " yin=", nrow(yin), " yang=", nrow(yang), "\n", sep = "")
if (parsed$stage == "preflight") quit(status = 0L)

output_files <- file.path(paths$output_root, c(
  "scores_yin.csv.gz", "scores_yang.csv.gz", "roc_native.csv.gz", "metrics.csv",
  "parameter_audit.csv", "input_manifest.csv", "COMPLETE"
))
if (all(file.exists(output_files))) {
  cat("PRADEEP_STRICT_PROJECTION_PRESERVED\n")
  quit(status = 0L)
}
if (any(file.exists(output_files))) {
  stop("Partial Pradeep strict projection exists; refusing overwrite: ", paths$output_root, call. = FALSE)
}

core_file <- file.path(paths$script_root, "yy", "R", "core.R")
score_require(core_file, "YY float32 helper")
source(core_file, local = FALSE)
x_yin <- read_f32_fortran_columns(
  common_files[["protein_yin"]], nrow(yin), length(features), common_map$index
)
x_yang <- read_f32_fortran_columns(
  common_files[["protein_yang"]], nrow(yang), length(features), common_map$index
)
common_eid <- c(yin$eid, yang$eid)
common_selected <- rbind(x_yin, x_yang)
strict_train_eid <- as.character(strict_dat$eid[train_index])
common_train_index <- match(strict_train_eid, common_eid)
scale_rows <- lapply(seq_len(nrow(selected)), function(column) {
  native_value <- suppressWarnings(as.numeric(strict_dat[[base_map$name[[column]]]][train_index]))
  projected_value <- common_selected[common_train_index, column]
  ok <- !is.na(common_train_index) & is.finite(native_value) & is.finite(projected_value)
  if (sum(ok) < 1000L) stop("Too few scale-audit rows for ", selected$feature[[column]], call. = FALSE)
  correlation <- suppressWarnings(stats::cor(native_value[ok], projected_value[ok]))
  fit <- stats::lm.fit(cbind(1, projected_value[ok]), native_value[ok])
  coefficient <- as.numeric(fit$coefficients)
  if (!is.finite(correlation) || correlation < 0.98 || length(coefficient) != 2L ||
      any(!is.finite(coefficient)) || coefficient[[2L]] <= 0) {
    stop("Pradeep scale bridge failed for ", selected$feature[[column]],
         ": correlation=", format(correlation, digits = 5), call. = FALSE)
  }
  data.table(
    feature = selected$feature[[column]], scale_overlap_n = sum(ok),
    scale_correlation = correlation, scale_intercept = coefficient[[1L]],
    scale_slope = coefficient[[2L]]
  )
})
scale_audit <- rbindlist(scale_rows)
project_one <- function(matrix_value) {
  matrix_value <- sweep(matrix_value, 2L, scale_audit$scale_slope, `*`)
  matrix_value <- sweep(matrix_value, 2L, scale_audit$scale_intercept, `+`)
  for (column in seq_len(ncol(matrix_value))) {
    bad <- !is.finite(matrix_value[, column])
    if (any(bad)) matrix_value[bad, column] <- training_median[[column]]
  }
  as.numeric(intercept[[1L]] + matrix_value %*% selected$coefficient)
}
yin_raw <- project_one(x_yin)
yang_raw <- project_one(x_yang)
center <- mean(yin_raw)
scale <- sd(yin_raw)
if (!is.finite(center) || !is.finite(scale) || scale <= 0 ||
    any(!is.finite(yin_raw)) || any(!is.finite(yang_raw))) {
  stop("Pradeep strict projected score scale failed.", call. = FALSE)
}
yin_output <- data.table(
  eid = yin$eid, method = "pradeep-strict", cohort_side = "Yin",
  score_raw = yin_raw, score_z = (yin_raw - center) / scale,
  score_source = "fixed native Pradeep protein model projected to common Yin"
)
yang_output <- data.table(
  eid = yang$eid, method = "pradeep-strict", cohort_side = "Yang",
  score_raw = yang_raw, score_z = (yang_raw - center) / scale,
  score_source = "fixed native Pradeep protein model projected to common Yang"
)

native <- fread(prediction_file, colClasses = list(character = "eid"))[
  model == "Proteins" & outc == "cad",
  .(eid, y = as.integer(y), score = as.numeric(score_link))
]
if (nrow(native) < 1000L || length(unique(native$y)) != 2L || anyDuplicated(native$eid)) {
  stop("Pradeep strict native prediction contract failed.", call. = FALSE)
}
replay <- merge(
  native,
  data.table(eid = common_eid, replayed_score = c(yin_raw, yang_raw)),
  by = "eid", all = FALSE, sort = FALSE
)
replay_correlation <- suppressWarnings(cor(replay$score, replay$replayed_score))
replay_rank_correlation <- suppressWarnings(cor(replay$score, replay$replayed_score, method = "spearman"))
replay_rmse <- sqrt(mean((replay$score - replay$replayed_score)^2))
if (nrow(replay) < 1000L || !is.finite(replay_correlation) || replay_correlation < 0.95) {
  stop("Pradeep frozen-score replay validation failed: n=", nrow(replay),
       " correlation=", format(replay_correlation, digits = 5), call. = FALSE)
}
roc_object <- pROC::roc(native$y, native$score, levels = c(0, 1), direction = "<", quiet = TRUE)
roc_curve <- data.table(
  method = "pradeep-strict", series = "Pradeep strict",
  false_positive_rate = 1 - as.numeric(roc_object$specificities),
  true_positive_rate = as.numeric(roc_object$sensitivities)
)
interval <- as.numeric(pROC::ci.auc(roc_object, method = "delong"))
metrics <- data.table(
  method = "pradeep-strict", level = "native_heldout", n = nrow(native),
  events = sum(native$y), auc = as.numeric(pROC::auc(roc_object)),
  auc_ci_low = interval[[1L]], auc_ci_high = interval[[3L]],
  estimand = "native held-out eventual incident CAD binary AUC",
  projection_validation_n = nrow(replay),
  projection_validation_pearson = replay_correlation,
  projection_validation_spearman = replay_rank_correlation,
  projection_validation_rmse = replay_rmse
)
parameter_audit <- cbind(data.table(
  feature = selected$feature, coefficient = selected$coefficient,
  strict_base_feature = base_map$name, common_feature = common_map$name,
  strict_training_median = training_median
), scale_audit[, -"feature"])
manifest <- data.table(
  role = c("strict_base", "strict_coefficients", "strict_predictions", names(common_files), "core"),
  path = c(base_file, coefficient_file, prediction_file, unname(common_files), core_file)
)
manifest[, `:=`(bytes = as.numeric(file.info(path)$size), md5 = unname(tools::md5sum(path)))]

dir.create(paths$output_root, recursive = TRUE, showWarnings = FALSE)
score_atomic_csv(yin_output, output_files[[1L]])
score_atomic_csv(yang_output, output_files[[2L]])
score_atomic_csv(roc_curve, output_files[[3L]])
score_atomic_csv(metrics, output_files[[4L]])
score_atomic_csv(parameter_audit, output_files[[5L]])
score_atomic_csv(manifest, output_files[[6L]])
score_atomic_text(c(
  "status=COMPLETE", "method=pradeep-strict",
  paste0("source_root=", source_root), paste0("selected_proteins=", nrow(selected)),
  paste0("yin_n=", nrow(yin)), paste0("yang_n=", nrow(yang)),
  paste0("yin_center=", format(center, digits = 16)),
  paste0("yin_scale=", format(scale, digits = 16)),
  paste0("native_auc=", format(metrics$auc[[1L]], digits = 16)),
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
), output_files[[7L]])
cat("PRADEEP_STRICT_PROJECTION_COMPLETE\nOUTPUT_ROOT=", paths$output_root,
    "\nNATIVE_AUC=", sprintf("%.6f", metrics$auc[[1L]]), "\n", sep = "")
