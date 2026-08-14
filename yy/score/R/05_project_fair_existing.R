#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = "") {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (length(hit)) sub("^[^=]+=?", "", hit[[length(hit)]]) else default
}

method <- tolower(arg("method"))
source_root <- gsub("\\\\", "/", path.expand(arg("source-root")))
output_root <- arg("output-root")
if (!method %in% c("pradeep-fair", "yu-fair")) {
  stop("--method must be pradeep-fair or yu-fair.", call. = FALSE)
}
if (!nzchar(output_root)) stop("--output-root is required.", call. = FALSE)
if (!dir.exists(source_root)) stop("--source-root does not exist: ", source_root, call. = FALSE)

source_label <- if (method == "pradeep-fair") {
  "Pradeep-style LASSO-logistic"
} else {
  "Yu-style LightGBM"
}
source_data <- file.path(source_root, "03_source_data")
yin_file <- file.path(source_data, "unified_yin_oof_scores.csv.gz")
yang_file <- file.path(source_data, "unified_yang_ensemble_scores.csv.gz")
metric_file <- file.path(source_data, "panel_b_unified_model_metrics.csv")
fold_file <- file.path(source_data, "panel_b_unified_fold_auc5.csv")
roc_file <- file.path(source_data, "panel_b_unified_ipcw_roc_curves.csv")
required <- c(yin_file, yang_file, metric_file, fold_file, roc_file)
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop("Completed fair-comparison project is incomplete: ",
       paste(missing, collapse = "; "), call. = FALSE)
}

outputs <- file.path(output_root, c(
  "scores_yin_oof.csv.gz", "scores_yang.csv.gz", "metrics.csv",
  "fold_metrics.csv", "roc_curves.csv.gz", "input_manifest.csv", "COMPLETE"
))
if (all(file.exists(outputs))) {
  cat("FAIR_PROJECT_PRESERVED method=", method, "\n", sep = "")
  quit(status = 0L)
}
if (any(file.exists(outputs))) {
  stop("Partial projected fair output exists; refusing overwrite: ",
       output_root, call. = FALSE)
}

yin <- fread(yin_file, colClasses = list(character = "eid"))[series == source_label]
yang <- fread(yang_file, colClasses = list(character = "eid"))[series == source_label]
metrics <- fread(metric_file)[series == source_label]
fold_metrics <- fread(fold_file)[series == source_label]
roc <- fread(roc_file)[series == source_label]

if (nrow(yin) != 37127L || nrow(yang) != 1766L ||
    sum(yin$event) != 3442L || anyDuplicated(yin$eid) ||
    anyDuplicated(yang$eid) || nrow(metrics) != 1L ||
    nrow(fold_metrics) != 5L || !nrow(roc) ||
    any(!is.finite(yin$score_raw)) || any(!is.finite(yin$score_z)) ||
    any(!is.finite(yang$score_raw)) || any(!is.finite(yang$score_z))) {
  stop("Existing fair-comparison score contract failed for ", method,
       call. = FALSE)
}

yin_output <- yin[, .(
  eid, outer_fold = as.integer(outer_fold), time = as.numeric(time),
  event = as.integer(event), method = method, cohort_side = "Yin",
  score_raw = as.numeric(score_raw), score_z = as.numeric(score_z),
  score_source = paste0("existing locked common-fold project: ", source_label)
)]
yang_output <- yang[, .(
  eid, method = method, cohort_side = "Yang",
  score_raw = as.numeric(score_raw), score_z = as.numeric(score_z),
  score_source = paste0("existing locked common-fold ensemble: ", source_label)
)]
metrics_output <- metrics[, .(
  method = method, level = "pooled_oof", n = as.integer(n),
  events = as.integer(incident_events), AUC_5y = as.numeric(AUC_5y),
  estimand = "locked common five-fold 5-year IPCW AUC"
)]
fold_output <- fold_metrics[, .(
  method = method, outer_fold = as.integer(outer_fold),
  AUC_5y = as.numeric(AUC_5y)
)]
roc_output <- roc[, .(
  method = method, false_positive_rate = as.numeric(false_positive_rate),
  true_positive_rate = as.numeric(true_positive_rate)
)]
manifest <- data.table(
  role = c("yin_oof_scores", "yang_ensemble_scores", "metrics", "fold_metrics", "roc"),
  source_path = required,
  source_project = source_root,
  projection_action = "column selection and method relabelling only; no fitting or tuning"
)

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
atomic_write <- function(object, path) {
  temporary <- if (grepl("\\.gz$", path)) {
    paste0(sub("\\.gz$", "", path), ".tmp.", Sys.getpid(), ".gz")
  } else {
    paste0(path, ".tmp.", Sys.getpid())
  }
  fwrite(object, temporary)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Atomic output rename failed: ", path, call. = FALSE)
  }
}
atomic_write(yin_output, outputs[[1L]])
atomic_write(yang_output, outputs[[2L]])
atomic_write(metrics_output, outputs[[3L]])
atomic_write(fold_output, outputs[[4L]])
atomic_write(roc_output, outputs[[5L]])
atomic_write(manifest, outputs[[6L]])
writeLines(c(
  paste0("method=", method),
  paste0("source_project=", source_root),
  "action=project_existing_scores_without_refit"
), outputs[[7L]])
cat("FAIR_PROJECT_COMPLETE method=", method,
    " AUC5=", format(metrics_output$AUC_5y, digits = 8), "\n", sep = "")
