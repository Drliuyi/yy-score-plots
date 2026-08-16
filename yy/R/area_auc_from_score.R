#!/usr/bin/env Rscript

## Cross-fitted association between CAD Yin/Yang trajectory area and AUC5.
##
## This script reuses the locked common-cohort products created by
## `yy score --main`. It never refits the four Pradeep/Yu score models.

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit)) return(sub(prefix, "", hit[[length(hit)]], fixed = TRUE))
  index <- match(paste0("--", name), args)
  if (!is.na(index) && index < length(args)) return(args[[index + 1L]])
  default
}
has_flag <- function(name) paste0("--", name) %in% args
env_or <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

if (has_flag("h") || has_flag("help")) {
  cat(paste(
    "Usage:",
    "  Rscript area_auc_from_score.R --stage=preflight",
    "  Rscript area_auc_from_score.R --stage=compute --workers=10 --resume",
    "  Rscript area_auc_from_score.R --stage=status",
    "",
    "Optional smoke test:",
    "  --max-proteins=6 --output-root=/mnt/d/analysis/yy/area-auc-smoke",
    "",
    "Default input:",
    "  D:/analysis/yy/score/common-fair-inputs",
    "",
    "Default output:",
    "  D:/analysis/yy/area-auc",
    "",
    "The compute stage also writes one four-panel figure as PNG, PDF and SVG",
    "under 05_figures, plus the exact figure source tables under 02_summary.",
    sep = "\n"
  ))
  quit(status = 0L)
}

stage <- tolower(arg_value("stage", "preflight"))
if (!stage %in% c("preflight", "compute", "status")) {
  stop("--stage must be preflight, compute or status.", call. = FALSE)
}
workers <- as.integer(arg_value("workers", "10"))
if (!is.finite(workers) || workers < 1L) stop("--workers must be positive.", call. = FALSE)
resume <- has_flag("resume")
max_proteins <- as.integer(arg_value("max-proteins", "2910"))
if (!is.finite(max_proteins) || max_proteins < 1L) {
  stop("--max-proteins must be positive.", call. = FALSE)
}

dir0 <- arg_value("dir0", env_or("DIR0", "/mnt/d"))
script_root <- arg_value("script-root", env_or("SCRIPT_ROOT", file.path(dir0, "scripts")))
analysis_root <- arg_value("analysis-root", env_or("ANALYSIS_ROOT", file.path(dir0, "analysis")))
yy_outdir <- arg_value("yy-outdir", env_or("YY_OUTDIR", file.path(analysis_root, "yy")))
score_root <- arg_value("score-root", file.path(yy_outdir, "score"))
common_root <- arg_value("common-root", file.path(score_root, "common-fair-inputs"))
output_root <- arg_value("output-root", file.path(yy_outdir, "area-auc"))
core_file <- file.path(script_root, "yy", "R", "core.R")

paths <- list(
  yin = file.path(common_root, "participants_yin.csv.gz"),
  yang = file.path(common_root, "participants_yang.csv.gz"),
  yin_target = file.path(common_root, "locked_yin_target.rds"),
  yang_target = file.path(common_root, "locked_yang_target.rds"),
  features = file.path(common_root, "protein_features.csv"),
  protein_yin = file.path(common_root, "protein_yin.f32"),
  protein_yang = file.path(common_root, "protein_yang.f32"),
  score_complete = file.path(score_root, "MAIN_CAD_COMPLETE.tsv"),
  core = core_file
)

cat("Resolved Huang-lab roots:\n")
cat("  DIR0=", dir0,
    "\n  SCRIPT_ROOT=", script_root,
    "\n  ANALYSIS_ROOT=", analysis_root,
    "\n  YY_OUTDIR=", yy_outdir,
    "\n  SCORE_ROOT=", score_root,
    "\n  COMMON_ROOT=", common_root,
    "\n  OUTPUT_ROOT=", output_root, "\n", sep = "")

complete_file <- file.path(output_root, "COMPLETE")
fold_dir <- file.path(output_root, "01_fold_metrics")
summary_dir <- file.path(output_root, "02_summary")
qc_dir <- file.path(output_root, "03_qc")
report_dir <- file.path(output_root, "04_report")
figure_dir <- file.path(output_root, "05_figures")
figure_stem <- file.path(figure_dir, "Figure_CAD_YinYang_area_AUC_four_panel_crossfitted")
expected_figure_files <- paste0(figure_stem, c(".png", ".pdf", ".svg"))

if (stage == "status") {
  fold_files <- file.path(fold_dir, sprintf("fold%02d.csv.gz", 1:5))
  status_ok <- file.exists(complete_file) && all(file.exists(expected_figure_files))
  cat("STATUS=", if (status_ok) "COMPLETE" else "NOT_COMPLETE", "\n", sep = "")
  cat(paste0("fold", 1:5, "=", ifelse(file.exists(fold_files), "COMPLETE", "MISSING")), sep = "\n")
  cat("\n")
  cat(paste0("figure_", c("png", "pdf", "svg"), "=",
             ifelse(file.exists(expected_figure_files), "COMPLETE", "MISSING")), sep = "\n")
  cat("\n")
  quit(status = if (status_ok) 0L else 1L)
}

missing <- names(paths)[!file.exists(unlist(paths, use.names = FALSE))]
if (length(missing)) {
  stop("Missing yy score input(s): ",
       paste(sprintf("%s=%s", missing, unlist(paths)[missing]), collapse = "; "),
       call. = FALSE)
}
source(core_file, local = FALSE)

yin <- fread(paths$yin, colClasses = list(character = "eid"))
yang <- fread(paths$yang, colClasses = list(character = "eid"))
yin_target <- as.data.table(readRDS(paths$yin_target)); yin_target[, eid := as.character(eid)]
yang_target <- as.data.table(readRDS(paths$yang_target)); yang_target[, eid := as.character(eid)]
feature_table <- fread(paths$features)
features <- as.character(feature_table$feature)

required_yin <- c("eid", "time", "event", "outer_fold")
required_yang <- c("eid", "disease_duration_years")
if (!all(required_yin %in% names(yin)) || !all(required_yang %in% names(yang)) ||
    nrow(yin) != 37127L || sum(yin$event) != 3442L || nrow(yang) != 1766L ||
    nrow(yin_target) != nrow(yin) || nrow(yang_target) != nrow(yang) ||
    length(features) != 2910L || anyDuplicated(yin$eid) || anyDuplicated(yang$eid) ||
    anyDuplicated(yin_target$eid) || anyDuplicated(yang_target$eid) ||
    !setequal(sort(unique(yin$outer_fold)), 1:5)) {
  stop("Locked yy score cohort contract failed.", call. = FALSE)
}
if (!identical(yin$eid, yin_target$eid)) {
  index <- match(yin$eid, yin_target$eid)
  if (anyNA(index)) stop("Yin target EID alignment failed.", call. = FALSE)
  yin_target <- yin_target[index]
}
if (!identical(yang$eid, yang_target$eid)) {
  index <- match(yang$eid, yang_target$eid)
  if (anyNA(index)) stop("Yang target EID alignment failed.", call. = FALSE)
  yang_target <- yang_target[index]
}
if (any(abs(yin$time - yin_target$time_years) > 1e-8) ||
    any(abs(yang$disease_duration_years - yang_target$disease_duration_years) > 1e-8)) {
  stop("Locked target time alignment failed.", call. = FALSE)
}

expected_yin_bytes <- as.double(nrow(yin)) * length(features) * 4
expected_yang_bytes <- as.double(nrow(yang)) * length(features) * 4
if (file.info(paths$protein_yin)$size != expected_yin_bytes ||
    file.info(paths$protein_yang)$size != expected_yang_bytes) {
  stop("Float32 protein matrix contract failed.", call. = FALSE)
}

max_proteins <- min(max_proteins, length(features))
selected <- seq_len(max_proteins)
cat(sprintf(
  "AREA_AUC_PREFLIGHT_PASS yin=%d events=%d yang=%d proteins=%d workers=%d\n",
  nrow(yin), sum(yin$event), nrow(yang), length(selected), workers
))
if (stage == "preflight") quit(status = 0L)

dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
fold_files <- file.path(fold_dir, sprintf("fold%02d.csv.gz", 1:5))
reuse_completed_folds <- FALSE
if (file.exists(complete_file)) {
  complete_lines <- readLines(complete_file, warn = FALSE)
  expected_line <- paste0("proteins=", length(selected))
  if (expected_line %in% complete_lines) {
    if (all(file.exists(expected_figure_files))) {
      cat("AREA_AUC_COMPLETE preserved\n")
      quit(status = 0L)
    }
    if (all(file.exists(fold_files))) {
      reuse_completed_folds <- TRUE
      resume <- TRUE
      cat("Completed fold metrics found; rebuilding summaries and figures only.\n")
    } else {
      stop("COMPLETE exists but fold outputs needed for figure rebuilding are missing: ",
           output_root, call. = FALSE)
    }
  } else {
    stop("Completed output was created with a different protein count: ",
         output_root, call. = FALSE)
  }
}
if (!reuse_completed_folds && !resume && any(file.exists(fold_files))) {
  stop("Partial fold outputs exist; rerun with --resume: ", output_root,
       call. = FALSE)
}

## Read the locked matrices once. Under WSL/Linux, mclapply workers share these
## pages copy-on-write, avoiding ten independent 400-MB input copies.
if (!reuse_completed_folds) {
  cat("Loading locked score matrices...\n")
  if (length(selected) == length(features)) {
    protein_yin <- read_f32_fortran(paths$protein_yin, nrow(yin), length(features))
    protein_yang <- read_f32_fortran(paths$protein_yang, nrow(yang), length(features))
  } else {
    protein_yin <- read_f32_fortran_columns(
      paths$protein_yin, nrow(yin), length(features), selected
    )
    protein_yang <- read_f32_fortran_columns(
      paths$protein_yang, nrow(yang), length(features), selected
    )
  }
}

## The plotted trajectory uses two-year bins centred at 1, 3 and 5 years.
## For the exact 0-5-year window their overlap widths are 2, 2 and 1 years.
area_0_5 <- function(time, value) {
  keep <- is.finite(time) & time >= 0 & time <= 5 & is.finite(value)
  if (!any(keep)) return(c(signed = NA_real_, absolute = NA_real_, mean1 = NA_real_, mean3 = NA_real_, mean5 = NA_real_))
  bin <- baseline_bin(time[keep])
  d <- data.table(bin = bin, value = value[keep])[bin %in% c(1L, 3L, 5L),
                                                        .(mean = mean(value)), by = bin]
  means <- setNames(rep(NA_real_, 3L), c("1", "3", "5"))
  means[as.character(d$bin)] <- d$mean
  if (any(!is.finite(means))) {
    return(c(signed = NA_real_, absolute = NA_real_, mean1 = means[[1L]], mean3 = means[[2L]], mean5 = means[[3L]]))
  }
  weight <- c(2, 2, 1)
  c(signed = sum(weight * means), absolute = sum(weight * abs(means)),
    mean1 = means[[1L]], mean3 = means[[2L]], mean5 = means[[3L]])
}

one_protein <- function(j, fold) {
  train <- yin$outer_fold != fold
  test <- !train
  raw_train <- protein_yin[train, j]
  center <- mean(raw_train, na.rm = TRUE)
  if (!is.finite(center)) center <- 0
  raw_train[!is.finite(raw_train)] <- center
  scale <- sd(raw_train)
  if (!is.finite(scale) || scale <= 1e-8) return(NULL)
  z_train <- (raw_train - center) / scale

  fit <- tryCatch(
    suppressWarnings(coxph(Surv(yin$time[train], yin$event[train]) ~ z_train,
                             ties = "efron", singular.ok = TRUE)),
    error = function(e) NULL
  )
  beta <- if (is.null(fit)) NA_real_ else suppressWarnings(as.numeric(coef(fit)[[1L]]))
  if (!is.finite(beta)) return(NULL)
  direction <- if (beta < 0) -1 else 1

  raw_test <- protein_yin[test, j]
  raw_test[!is.finite(raw_test)] <- center
  z_test <- (raw_test - center) / scale
  auc <- tryCatch(
    ipcw_roc(yin$time[test], yin$event[test], beta * z_test, horizon = 5)$auc,
    error = function(e) NA_real_
  )

  incident_train <- yin$event[train] == 1L
  yin_area <- area_0_5(yin$time[train][incident_train], z_train[incident_train])
  raw_yang <- protein_yang[, j]
  raw_yang[!is.finite(raw_yang)] <- center
  z_yang <- (raw_yang - center) / scale
  yang_area <- area_0_5(yang$disease_duration_years, z_yang)

  data.table(
    protein = features[selected[[j]]], outer_fold = fold,
    beta_train = beta, risk_direction = direction, AUC_5y = auc,
    AUC_5y_margin = auc - 0.5,
    yin_signed_area_0_5 = yin_area[["signed"]],
    yin_absolute_area_0_5 = yin_area[["absolute"]],
    yin_risk_aligned_area_0_5 = direction * yin_area[["signed"]],
    yang_signed_area_0_5 = yang_area[["signed"]],
    yang_absolute_area_0_5 = yang_area[["absolute"]],
    yang_risk_aligned_area_0_5 = direction * yang_area[["signed"]],
    yin_mean_bin1 = yin_area[["mean1"]], yin_mean_bin3 = yin_area[["mean3"]],
    yin_mean_bin5 = yin_area[["mean5"]], yang_mean_bin1 = yang_area[["mean1"]],
    yang_mean_bin3 = yang_area[["mean3"]], yang_mean_bin5 = yang_area[["mean5"]]
  )
}

for (fold in 1:5) {
  fold_file <- fold_files[[fold]]
  if (resume && file.exists(fold_file)) {
    existing <- tryCatch(fread(fold_file), error = function(e) NULL)
    if (!is.null(existing) && nrow(existing) == length(selected) &&
        uniqueN(existing$protein) == length(selected)) {
      cat(sprintf("fold%d=COMPLETE preserved\n", fold))
      next
    }
    stop("Existing fold output is incomplete; refusing overwrite: ", fold_file, call. = FALSE)
  }
  cat(sprintf("Computing fold %d/5 across %d proteins with %d workers...\n",
              fold, length(selected), workers))
  if (.Platform$OS.type == "unix" && workers > 1L) {
    rows <- parallel::mclapply(seq_along(selected), one_protein, fold = fold,
                               mc.cores = workers, mc.preschedule = TRUE)
  } else {
    rows <- lapply(seq_along(selected), one_protein, fold = fold)
  }
  if (any(vapply(rows, is.null, logical(1L)))) {
    failed <- which(vapply(rows, is.null, logical(1L)))
    stop("Protein calculation failed in fold ", fold, ": ",
         paste(features[selected[failed]][seq_len(min(10L, length(failed)))], collapse = ", "),
         call. = FALSE)
  }
  fold_result <- rbindlist(rows)
  if (nrow(fold_result) != length(selected) ||
      any(!is.finite(as.matrix(fold_result[, .(
        beta_train, AUC_5y, yin_signed_area_0_5, yin_absolute_area_0_5,
        yin_risk_aligned_area_0_5, yang_signed_area_0_5,
        yang_absolute_area_0_5, yang_risk_aligned_area_0_5
      )])))) {
    stop("Fold result contract failed: ", fold, call. = FALSE)
  }
  atomic_csv(fold_result, fold_file)
  cat(sprintf("fold%d=COMPLETE\n", fold))
}

fold_metrics <- rbindlist(lapply(fold_files, fread))
if (nrow(fold_metrics) != 5L * length(selected) ||
    anyDuplicated(fold_metrics[, .(protein, outer_fold)])) {
  stop("Merged fold metric contract failed.", call. = FALSE)
}

rank_spline_residual <- function(value, control) {
  result <- rep(NA_real_, length(value))
  ok <- complete.cases(value, control)
  if (sum(ok) < 6L) return(result)
  value_rank <- rank(value[ok], ties.method = "average")
  control_rank <- rank(control[ok], ties.method = "average")
  result[ok] <- residuals(lm(value_rank ~ splines::ns(control_rank, df = 3)))
  result
}

partial_spearman <- function(x, y, control) {
  ok <- complete.cases(x, y, control)
  rx <- rank_spline_residual(x[ok], control[ok])
  ry <- rank_spline_residual(y[ok], control[ok])
  cor(rx, ry, method = "spearman")
}

spearman_result <- function(x, y) {
  ok <- complete.cases(x, y)
  test <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(n = sum(ok), estimate = unname(test$estimate), p_value = test$p.value)
}

fold_correlations <- fold_metrics[, .(
  n_proteins = .N,
  yin_risk_area_vs_auc5 = cor(yin_risk_aligned_area_0_5, AUC_5y_margin, method = "spearman"),
  yin_absolute_area_vs_auc5 = cor(yin_absolute_area_0_5, AUC_5y_margin, method = "spearman"),
  yin_area_vs_yang_area = cor(yin_risk_aligned_area_0_5, yang_risk_aligned_area_0_5,
                              method = "spearman"),
  yang_area_vs_auc5 = cor(yang_risk_aligned_area_0_5, AUC_5y_margin, method = "spearman"),
  yang_partial_vs_auc5_given_yin = partial_spearman(
    yang_risk_aligned_area_0_5, AUC_5y_margin, yin_risk_aligned_area_0_5
  )
), by = outer_fold]

protein_summary <- fold_metrics[, .(
  beta_train_mean = mean(beta_train), direction_consistency = abs(mean(risk_direction)),
  AUC_5y = mean(AUC_5y), AUC_5y_sd = sd(AUC_5y),
  yin_signed_area_0_5 = mean(yin_signed_area_0_5),
  yin_absolute_area_0_5 = mean(yin_absolute_area_0_5),
  yin_risk_aligned_area_0_5 = mean(yin_risk_aligned_area_0_5),
  yang_signed_area_0_5 = mean(yang_signed_area_0_5),
  yang_absolute_area_0_5 = mean(yang_absolute_area_0_5),
  yang_risk_aligned_area_0_5 = mean(yang_risk_aligned_area_0_5)
), by = protein]

protein_summary[, continuum_class := fcase(
  yin_risk_aligned_area_0_5 * yang_risk_aligned_area_0_5 < 0, "Reversal",
  abs(yin_risk_aligned_area_0_5) > 1.5 * abs(yang_risk_aligned_area_0_5), "Yin-dominant",
  abs(yang_risk_aligned_area_0_5) > 1.5 * abs(yin_risk_aligned_area_0_5), "Yang-dominant",
  default = "Continuum"
)]
protein_summary[, `:=`(
  yin_rank = rank(yin_risk_aligned_area_0_5, ties.method = "average"),
  yang_rank = rank(yang_risk_aligned_area_0_5, ties.method = "average"),
  auc5_rank = rank(AUC_5y, ties.method = "average")
)]
protein_summary[, `:=`(
  yang_residual_given_yin = rank_spline_residual(
    yang_risk_aligned_area_0_5, yin_risk_aligned_area_0_5
  ),
  auc5_residual_given_yin = rank_spline_residual(AUC_5y, yin_risk_aligned_area_0_5),
  log2_yin_yang_absolute_ratio = log2(
    (abs(yin_risk_aligned_area_0_5) + 1e-8) /
      (abs(yang_risk_aligned_area_0_5) + 1e-8)
  )
)]

trajectory_fold <- rbindlist(list(
  fold_metrics[, .(protein, outer_fold, side = "Yang", years_from_baseline = -5,
                   mean_risk_aligned = risk_direction * yang_mean_bin5)],
  fold_metrics[, .(protein, outer_fold, side = "Yang", years_from_baseline = -3,
                   mean_risk_aligned = risk_direction * yang_mean_bin3)],
  fold_metrics[, .(protein, outer_fold, side = "Yang", years_from_baseline = -1,
                   mean_risk_aligned = risk_direction * yang_mean_bin1)],
  fold_metrics[, .(protein, outer_fold, side = "Yin", years_from_baseline = 1,
                   mean_risk_aligned = risk_direction * yin_mean_bin1)],
  fold_metrics[, .(protein, outer_fold, side = "Yin", years_from_baseline = 3,
                   mean_risk_aligned = risk_direction * yin_mean_bin3)],
  fold_metrics[, .(protein, outer_fold, side = "Yin", years_from_baseline = 5,
                   mean_risk_aligned = risk_direction * yin_mean_bin5)]
))
trajectory_summary <- trajectory_fold[, .(
  mean_risk_aligned = mean(mean_risk_aligned),
  fold_sd = sd(mean_risk_aligned),
  contributing_folds = sum(is.finite(mean_risk_aligned))
), by = .(protein, side, years_from_baseline)]

preferred_proteins <- c("NTPROBNP", "MMP12", "GDF15", "HAVCR1", "EDA2R", "PCSK9")
protein_keys <- protein_key(protein_summary$protein)
preferred_index <- match(preferred_proteins, protein_keys)
representative_proteins <- protein_summary$protein[preferred_index[!is.na(preferred_index)]]
if (length(representative_proteins) < min(6L, nrow(protein_summary))) {
  fill <- protein_summary[order(-AUC_5y), protein]
  representative_proteins <- unique(c(representative_proteins, fill))
}
representative_proteins <- head(representative_proteins, min(6L, nrow(protein_summary)))

display_name <- function(protein) {
  labels <- as.character(protein)
  labels[protein_key(labels) == "NTPROBNP"] <- "NT-proBNP"
  labels
}
trajectory_summary[, display_label := display_name(protein)]
protein_summary[, display_label := display_name(protein)]

pooled_estimates <- list(
  spearman_result(protein_summary$yin_risk_aligned_area_0_5, protein_summary$AUC_5y - 0.5),
  spearman_result(protein_summary$yin_absolute_area_0_5, protein_summary$AUC_5y - 0.5),
  spearman_result(protein_summary$yin_risk_aligned_area_0_5,
                  protein_summary$yang_risk_aligned_area_0_5),
  spearman_result(protein_summary$yang_risk_aligned_area_0_5, protein_summary$AUC_5y - 0.5),
  spearman_result(protein_summary$yang_residual_given_yin,
                  protein_summary$auc5_residual_given_yin)
)

pooled <- data.table(
  comparison = c(
    "Yin risk-aligned area vs AUC5",
    "Yin absolute area vs AUC5",
    "Yin risk-aligned area vs Yang risk-aligned area",
    "Yang risk-aligned area vs AUC5",
    "Yang area vs AUC5 controlling Yin area"
  ),
  method = c(rep("Spearman on five-fold protein means", 4L),
             "Partial Spearman after spline rank residualization"),
  n = vapply(pooled_estimates, `[[`, numeric(1L), "n"),
  estimate = vapply(pooled_estimates, `[[`, numeric(1L), "estimate"),
  p_value = vapply(pooled_estimates, `[[`, numeric(1L), "p_value")
)

atomic_csv(fold_metrics, file.path(summary_dir, "protein_fold_area_auc.csv.gz"))
atomic_csv(protein_summary, file.path(summary_dir, "protein_area_auc_summary.csv.gz"))
atomic_csv(trajectory_summary, file.path(summary_dir, "protein_trajectory_bin_summary.csv.gz"))
atomic_csv(fold_correlations, file.path(summary_dir, "fold_correlations.csv"))
atomic_csv(pooled, file.path(summary_dir, "correlation_summary.csv"))

class_levels <- c("Continuum", "Yin-dominant", "Yang-dominant", "Reversal")
class_counts <- merge(
  data.table(continuum_class = class_levels),
  protein_summary[, .N, by = continuum_class],
  by = "continuum_class", all.x = TRUE, sort = FALSE
)
class_counts[is.na(N), N := 0L]
atomic_csv(class_counts, file.path(summary_dir, "continuum_class_counts.csv"))

format_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 2.2e-16) return("<2.2e-16")
  formatC(p, digits = 3L, format = "g")
}
stat_label <- function(row) {
  sprintf("Spearman rho = %.3f; P %s; n = %d",
          row$estimate, if (row$p_value < 2.2e-16) format_p(row$p_value)
          else paste0("= ", format_p(row$p_value)), row$n)
}

representative_trajectory <- trajectory_summary[protein %in% representative_proteins]
representative_points <- protein_summary[protein %in% representative_proteins]
atomic_csv(representative_trajectory,
           file.path(summary_dir, "figure_panel_a_representative_trajectories.csv"))
atomic_csv(protein_summary,
           file.path(summary_dir, "figure_panels_bcd_protein_source.csv.gz"))
atomic_csv(pooled[c(1L, 3L, 5L)],
           file.path(summary_dir, "figure_displayed_statistics.csv"))

protein_palette <- c(
  "NT-proBNP" = "#C95B47", MMP12 = "#B38B00", GDF15 = "#159A62",
  HAVCR1 = "#1C93A5", EDA2R = "#337AC4", PCSK9 = "#BF3DBB"
)
missing_colours <- setdiff(unique(representative_trajectory$display_label), names(protein_palette))
if (length(missing_colours)) {
  protein_palette <- c(protein_palette,
                       setNames(scales::hue_pal()(length(missing_colours)), missing_colours))
}
class_palette <- c(
  Continuum = "#29967F", `Yin-dominant` = "#2D7FC1",
  `Yang-dominant` = "#F08A34", Reversal = "#8061B4"
)
class_labels <- setNames(
  sprintf("%s (n=%s)", class_counts$continuum_class,
          format(class_counts$N, big.mark = ",")),
  class_counts$continuum_class
)
common_theme <- theme_classic(base_size = 11, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 12.5),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(linewidth = 0.55),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8.5)
  )

panel_a <- ggplot(representative_trajectory,
                  aes(years_from_baseline, mean_risk_aligned,
                      colour = display_label, group = display_label)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "#F2EFEC", alpha = 0.70) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#EAF3FB", alpha = 0.70) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50", linewidth = 0.45) +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey30", linewidth = 0.55) +
  geom_line(linewidth = 1.05) +
  geom_point(shape = 21, fill = "white", stroke = 1.05, size = 3.1) +
  annotate("text", x = -4.3, y = Inf, label = "Yang", vjust = 1.7,
           colour = "#645D57", fontface = "bold", size = 4.2) +
  annotate("text", x = 4.3, y = Inf, label = "Yin", vjust = 1.7,
           colour = "#2878B8", fontface = "bold", size = 4.2) +
  scale_colour_manual(values = protein_palette, breaks = unique(representative_trajectory$display_label)) +
  scale_x_continuous(breaks = c(-5, -3, -1, 0, 1, 3, 5), limits = c(-5.4, 5.4)) +
  guides(colour = guide_legend(ncol = 3, byrow = TRUE)) +
  labs(title = "Baseline-centred Yin-Yang trajectories",
       x = "Years from baseline", y = "Mean risk-aligned protein z-score") +
  common_theme

highlight_b <- representative_points
panel_b <- ggplot(protein_summary,
                  aes(yin_risk_aligned_area_0_5, AUC_5y)) +
  geom_point(shape = 21, fill = "#4C86A8", colour = "#4C86A8",
             alpha = 0.28, size = 1.55) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "#BA4438", linewidth = 1.15) +
  geom_point(data = highlight_b, shape = 21, fill = "white", colour = "black",
             stroke = 0.9, size = 3) +
  geom_text(data = highlight_b, aes(label = display_label),
            hjust = -0.08, vjust = -0.45, size = 3.2, check_overlap = TRUE) +
  annotate("text", x = -Inf, y = -Inf, hjust = -0.02, vjust = -0.55,
           label = stat_label(pooled[1L]), size = 3.35) +
  labs(title = "Yin area and single-protein AUC5",
       x = "Yin risk-aligned area, 0-5 years", y = "Five-year IPCW AUC") +
  coord_cartesian(clip = "off") + common_theme +
  theme(plot.margin = margin(7, 19, 20, 7))

protein_summary[, continuum_class_plot := factor(
  continuum_class, levels = class_levels
)]
panel_c <- ggplot(protein_summary,
                  aes(yin_risk_aligned_area_0_5, yang_risk_aligned_area_0_5,
                      colour = continuum_class_plot)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.4) +
  geom_abline(slope = 1.5, intercept = 0, linetype = 3, colour = "grey65") +
  geom_abline(slope = 1 / 1.5, intercept = 0, linetype = 3, colour = "grey65") +
  geom_point(alpha = 0.42, size = 1.35) +
  geom_point(data = representative_points, inherit.aes = FALSE,
             aes(yin_risk_aligned_area_0_5, yang_risk_aligned_area_0_5),
             shape = 21, fill = "white", colour = "black", stroke = 0.9, size = 3) +
  geom_text(data = representative_points, inherit.aes = FALSE,
            aes(yin_risk_aligned_area_0_5, yang_risk_aligned_area_0_5,
                label = display_label),
            hjust = -0.08, vjust = -0.45, size = 3.2, check_overlap = TRUE) +
  annotate("text", x = -Inf, y = -Inf, hjust = -0.02, vjust = -0.55,
           label = stat_label(pooled[3L]), size = 3.35) +
  scale_colour_manual(values = class_palette, labels = class_labels, drop = FALSE) +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE)) +
  labs(title = "Four Yin-Yang area phenotypes",
       x = "Yin risk-aligned area, 0-5 years",
       y = "Yang risk-aligned area, 0-5 years") +
  coord_cartesian(clip = "off") + common_theme +
  theme(plot.margin = margin(7, 19, 20, 7))

panel_d <- ggplot(protein_summary,
                  aes(yang_residual_given_yin, auc5_residual_given_yin)) +
  geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.45) +
  geom_vline(xintercept = 0, colour = "grey65", linewidth = 0.45) +
  geom_point(shape = 21, fill = "grey55", colour = "grey45",
             alpha = 0.30, size = 1.45) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "#BA4438", linewidth = 1.15) +
  annotate("text", x = -Inf, y = -Inf, hjust = -0.02, vjust = -0.55,
           label = stat_label(pooled[5L]), size = 3.35) +
  labs(title = "Residual Yang area and residual AUC5",
       x = "Yang-area rank residual", y = "AUC5-rank residual") +
  coord_cartesian(clip = "off") + common_theme +
  theme(legend.position = "none", plot.margin = margin(7, 7, 20, 7))

four_panel <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 16))
figure_stem <- file.path(figure_dir, "Figure_CAD_YinYang_area_AUC_four_panel_crossfitted")
ggsave(paste0(figure_stem, ".png"), four_panel, width = 15.8, height = 11.8,
       dpi = 300, device = ragg::agg_png, bg = "white")
ggsave(paste0(figure_stem, ".pdf"), four_panel, width = 15.8, height = 11.8,
       device = cairo_pdf, bg = "white")
ggsave(paste0(figure_stem, ".svg"), four_panel, width = 15.8, height = 11.8,
       device = svglite::svglite, bg = "white")
figure_files <- paste0(figure_stem, c(".png", ".pdf", ".svg"))
if (any(!file.exists(figure_files)) || any(file.info(figure_files)$size <= 1000)) {
  stop("Four-panel figure output contract failed.", call. = FALSE)
}

manifest <- data.table(
  role = names(paths), path = normalizePath(unlist(paths), winslash = "/", mustWork = TRUE),
  bytes = as.double(file.info(unlist(paths))$size),
  modified = format(file.info(unlist(paths))$mtime, "%Y-%m-%dT%H:%M:%S%z")
)
atomic_csv(manifest, file.path(qc_dir, "input_manifest.csv"))
qc <- data.table(
  check = c("yin_n", "yin_events", "yang_n", "protein_n", "fold_rows",
            "finite_primary_metrics", "five_folds", "trajectory_rows",
            "figure_png", "figure_pdf", "figure_svg"),
  observed = c(nrow(yin), sum(yin$event), nrow(yang), nrow(protein_summary),
               nrow(fold_metrics),
               as.integer(all(is.finite(pooled$estimate))), uniqueN(fold_metrics$outer_fold),
               nrow(trajectory_summary), as.integer(file.exists(figure_files[[1L]])),
               as.integer(file.exists(figure_files[[2L]])),
               as.integer(file.exists(figure_files[[3L]]))),
  expected = c(37127, 3442, 1766, length(selected), 5L * length(selected), 1, 5,
               6L * length(selected), 1, 1, 1)
)
qc[, status := fifelse(observed == expected, "PASS", "FAIL")]
atomic_csv(qc, file.path(qc_dir, "qc.csv"))
if (any(qc$status != "PASS")) stop("Area-AUC QC failed.", call. = FALSE)

report <- c(
  "# Cross-fitted Yin-Yang trajectory area and AUC5",
  "",
  sprintf("- Yin: %d participants; %d incident CAD events.", nrow(yin), sum(yin$event)),
  sprintf("- Yang: %d prevalent CAD participants.", nrow(yang)),
  sprintf("- Proteins: %d; outer folds: 5.", length(selected)),
  "- Areas are learned in each outer-training fold from baseline-centred two-year bins.",
  "- AUC5 is evaluated only in the corresponding outer-test fold using IPCW.",
  "- Risk alignment uses only the sign of the training-fold Cox coefficient.",
  "- Yang partial correlation is evaluated after flexible rank control for Yin area.",
  "- The four-panel figure is generated from these same cross-fitted summaries.",
  "",
  "## Correlations",
  "",
  paste(sprintf("- %s: rho = %.6f", pooled$comparison, pooled$estimate), collapse = "\n"),
  "",
  "## Interpretation boundary",
  "",
  "These are cross-sectional case-timing pseudo-trajectories from baseline measurements, not within-person longitudinal trajectories. Correlation does not show that area weighting improves a multivariable prediction model.",
  "",
  "## Figure outputs",
  "",
  "- 05_figures/Figure_CAD_YinYang_area_AUC_four_panel_crossfitted.png",
  "- 05_figures/Figure_CAD_YinYang_area_AUC_four_panel_crossfitted.pdf",
  "- 05_figures/Figure_CAD_YinYang_area_AUC_four_panel_crossfitted.svg",
  "- Exact panel source tables are stored under 02_summary."
)
writeLines(report, file.path(report_dir, "RESULTS_SUMMARY.md"), useBytes = TRUE)
writeLines(c(
  "status=COMPLETE",
  paste0("proteins=", length(selected)),
  paste0("workers=", workers),
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
), complete_file)
cat("AREA_AUC_COMPLETE\n")
print(pooled)
