#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(survival)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
project_dir <- dirname(dirname(gsub("\\\\", "/", script_file)))
source(file.path(project_dir, "R", "00_common.R"))
parsed <- score_args()
method_arg <- commandArgs(trailingOnly = TRUE)
method_hit <- method_arg[startsWith(method_arg, "--method=")]
if (length(method_hit) != 1L) stop("Required: --method=pradeep-fair|yu-fair", call. = FALSE)
method <- sub("^--method=", "", method_hit[[1L]])
if (!method %in% c("pradeep-fair", "yu-fair")) stop("Unsupported fair method: ", method, call. = FALSE)
paths <- score_resolve(project_dir, parsed, method)
cfg <- jsonlite::fromJSON(paths$config, simplifyVector = TRUE)
score_print_roots(paths)

fold_yin <- file.path(paths$output_root, "predictions", sprintf("fold%02d_yin.csv.gz", 1:5))
fold_yang <- file.path(paths$output_root, "predictions", sprintf("fold%02d_yang.csv.gz", 1:5))
markers <- file.path(paths$output_root, "models", sprintf("fold%02d_COMPLETE", 1:5))
score_require(c(fold_yin, fold_yang, markers), paste(method, "fold output"))
final_paths <- c(
  file.path(paths$output_root, "scores_yin_oof.csv.gz"),
  file.path(paths$output_root, "scores_yang.csv.gz"),
  file.path(paths$output_root, "metrics.csv"),
  file.path(paths$output_root, "COMPLETE")
)
if (all(file.exists(final_paths))) {
  cat("FAIR_METHOD_PRESERVED method=", method, "\n", sep = "")
  quit(status = 0L)
}
if (any(file.exists(final_paths))) stop("Partial finalized method output exists: ", paths$output_root, call. = FALSE)

yin <- rbindlist(lapply(fold_yin, fread, colClasses = list(character = "eid")))
yang_members <- rbindlist(lapply(fold_yang, fread, colClasses = list(character = "eid")))
if (nrow(yin) != cfg$expected$yin_n || sum(yin$event) != cfg$expected$yin_events ||
    anyDuplicated(yin$eid) || nrow(yang_members) != 5L * cfg$expected$yang_n ||
    anyDuplicated(yang_members[, .(eid, outer_fold)]) ||
    any(yin$method != method) || any(yang_members$method != method)) {
  stop("Final fair-score coverage contract failed: ", method, call. = FALSE)
}
center <- mean(yin$score_raw)
scale <- stats::sd(yin$score_raw)
if (!is.finite(center) || !is.finite(scale) || scale <= 0) stop("Fair score scale is invalid.", call. = FALSE)
yin[, score_z := (score_raw - center) / scale]
yang <- yang_members[, .(score_raw = mean(score_raw)), by = .(eid, method)]
yang[, score_z := (score_raw - center) / scale]
if (nrow(yang) != cfg$expected$yang_n || anyDuplicated(yang$eid)) stop("Yang mean-score contract failed.", call. = FALSE)

ipcw_auc <- function(time, event, score, horizon) {
  case <- which(event == 1L & time <= horizon)
  control <- which(time > horizon)
  if (!length(case) || !length(control)) return(NA_real_)
  censor_fit <- survfit(Surv(time, 1L - event) ~ 1)
  censor_survival <- function(value, left = FALSE) {
    query <- if (left) pmax(0, value - 1e-10) else value
    pmax(summary(censor_fit, times = query, extend = TRUE)$surv, 1e-6)
  }
  case_weight <- 1 / censor_survival(time[case], TRUE)
  control_weight <- rep(1 / censor_survival(horizon), length(control))
  control_score <- score[control]
  order_index <- order(control_score)
  control_score <- control_score[order_index]
  control_weight <- control_weight[order_index]
  cumulative <- cumsum(control_weight)
  concordant <- vapply(score[case], function(value) {
    less <- findInterval(value, control_score, left.open = TRUE)
    less_equal <- findInterval(value, control_score)
    weight_less <- if (less > 0L) cumulative[[less]] else 0
    weight_equal <- if (less_equal > less) sum(control_weight[(less + 1L):less_equal]) else 0
    weight_less + 0.5 * weight_equal
  }, numeric(1L))
  sum(case_weight * concordant) / (sum(case_weight) * sum(control_weight))
}
fold_metrics <- yin[, .(
  n = .N, events = sum(event),
  AUC_5y = ipcw_auc(time, event, score_z, cfg$horizon_years)
), by = outer_fold]
metrics <- rbindlist(list(
  fold_metrics[, .(method = method, level = "outer_fold", outer_fold, n, events, AUC_5y)],
  fold_metrics[, .(
    method = method, level = "mean_outer_fold", outer_fold = NA_integer_,
    n = sum(n), events = sum(events), AUC_5y = mean(AUC_5y)
  )]
), use.names = TRUE)
score_atomic_csv(yin, final_paths[[1L]])
score_atomic_csv(yang, final_paths[[2L]])
score_atomic_csv(metrics, final_paths[[3L]])
score_atomic_text(c(
  "status=COMPLETE", paste0("method=", method),
  paste0("yin_center=", format(center, digits = 16)),
  paste0("yin_scale=", format(scale, digits = 16)),
  paste0("mean_outer_fold_AUC5=", format(metrics[level == "mean_outer_fold", AUC_5y], digits = 16)),
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
), final_paths[[4L]])
cat("FAIR_METHOD_COMPLETE method=", method,
    " AUC5=", sprintf("%.6f", metrics[level == "mean_outer_fold", AUC_5y]), "\n", sep = "")
