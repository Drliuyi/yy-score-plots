#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = "") {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (length(hit)) sub("^[^=]+=", "", hit[[length(hit)]]) else default
}
as_flag <- function(name) tolower(get_arg(name, "0")) %in% c("1", "true", "yes")
split_values <- function(x) {
  if (!nzchar(x)) return(character())
  value <- trimws(unlist(strsplit(x, "[,|]", perl = TRUE), use.names = FALSE))
  unique(value[nzchar(value)])
}
env_or <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}
first_existing <- function(paths, label) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) stop(label, " not found. Checked: ", paste(paths, collapse = "; "), call. = FALSE)
  normalizePath(hit[[1L]], winslash = "/", mustWork = TRUE)
}
slug <- function(x) gsub("(^-+|-+$)", "", tolower(gsub("[^A-Za-z0-9]+", "-", x)))
hash_text <- function(x) {
  result <- 0
  for (value in utf8ToInt(enc2utf8(x))) result <- (result * 131 + value) %% 2147483647
  sprintf("%08x", as.integer(result))
}

dir0 <- env_or("DIR0", "/mnt/d")
script_root <- env_or("SCRIPT_ROOT", file.path(dir0, "scripts"))
analysis_root <- env_or("ANALYSIS_ROOT", file.path(dir0, "analysis"))
yy_outdir <- env_or("YY_OUTDIR", file.path(analysis_root, "yy"))
output_root <- get_arg("output-root", file.path(yy_outdir, "plot"))

dynamic_helper <- file.path(script_root, "yy", "R", "dynamic_proteins.R")
score_helper <- file.path(script_root, "yy", "R", "dynamic_scores.R")
if (!file.exists(dynamic_helper) || !file.exists(score_helper)) {
  stop("YY plot helper missing under SCRIPT_ROOT: ", script_root, call. = FALSE)
}
source(dynamic_helper, local = FALSE)

method_table <- data.table(
  method = c("pradeep-strict", "yu-strict", "pradeep-fair", "yu-fair"),
  label = c("Pradeep strict", "Yu strict", "Pradeep fair", "Yu fair"),
  family = c("LASSO-logistic", "LightGBM", "LASSO-logistic", "LightGBM"),
  target = c("native eventual incident CAD", "native eventual incident CAD",
             "common 5-year incident CAD", "common 5-year incident CAD"),
  roc_estimand = c("native held-out binary AUC", "native held-out binary AUC",
                   "common five-fold IPCW AUC5", "common five-fold IPCW AUC5")
)
method_alias <- c(
  pradeep = "pradeep-fair", yu = "yu-fair",
  pradeepstrict = "pradeep-strict", yustrict = "yu-strict",
  pradeepfair = "pradeep-fair", yufair = "yu-fair"
)
normalize_method <- function(value) {
  key <- tolower(gsub("[^A-Za-z0-9]", "", value))
  result <- unname(method_alias[[key]])
  if (is.null(result)) {
    direct <- method_table$method[tolower(method_table$method) == tolower(value)]
    if (length(direct)) result <- direct[[1L]]
  }
  if (is.null(result) || !length(result)) {
    stop("Unknown score method: ", value, ". Use: ", paste(method_table$method, collapse = ", "), call. = FALSE)
  }
  result
}

if (any(args == "--status")) {
  cat("REGISTERED SCORE METHODS\n")
  print(method_table)
  cat("ADJUSTMENT SYNTAX\n")
  cat("  raw\n")
  cat("  age+sex\n")
  cat("  age+sex+pc+center+tdi\n")
  cat("  age+sex+le8+medications\n")
  cat("  tokens are independent; omitted tokens are not added silently\n")
  cat("SCORE OUTPUT STATUS\n")
  for (method_value in method_table$method) {
    marker <- file.path(yy_outdir, "score", method_value, "COMPLETE")
    cat("  ", method_value, "=", if (file.exists(marker)) "COMPLETE" else "MISSING", "\n", sep = "")
  }
  quit(status = 0L)
}

side_mode <- tolower(get_arg("side", "yy"))
anchor <- tolower(get_arg("anchor", "baseline"))
adjustment_spec <- yy_adjustment_spec(get_arg("adj", "raw"))
adjustment <- adjustment_spec$id
show_traj <- as_flag("traj")
show_roc <- as_flag("roc")
show_bar <- as_flag("bar")
show_mean <- as_flag("mean")
show_sd <- as_flag("sd")
show_sd_line <- as_flag("sd-line")
recompute <- as_flag("recompute")
require_score <- as_flag("require-score")
requested_proteins <- split_values(get_arg("proteins"))
requested_methods <- vapply(split_values(get_arg("scores")), normalize_method, character(1L))
requested_methods <- unique(requested_methods)

if (!side_mode %in% c("yin", "yang", "yy")) stop("--side must be yin, yang, or yy.", call. = FALSE)
if (!anchor %in% c("baseline", "diagnosis")) stop("--anchor must be baseline or diagnosis.", call. = FALSE)
if (!show_traj && !show_roc) stop("Select --traj and/or --roc.", call. = FALSE)
if (show_roc && side_mode == "yang") stop("Incident ROC is not valid with Yang-only display.", call. = FALSE)
if ((show_bar || show_mean || show_sd || show_sd_line) && !show_traj) stop("Trajectory annotations require --traj.", call. = FALSE)
if ((show_sd || show_sd_line) && adjustment != "raw") {
  stop("Observed SD is defined only for --adj=raw; adjusted marginal-mean uncertainty is SE, not SD.", call. = FALSE)
}
if (!length(requested_proteins) && !length(requested_methods)) stop("Select --proteins and/or --score.", call. = FALSE)

dynamic_paths <- yy_dynamic_paths(dir0, analysis_root, yy_outdir, script_root)
source(dynamic_paths$core, local = FALSE)
yin_target <- as.data.table(readRDS(dynamic_paths$yin_target)); yin_target[, eid := as.character(eid)]
yang_target <- as.data.table(readRDS(dynamic_paths$yang_target)); yang_target[, eid := as.character(eid)]
if (nrow(yin_target) != 37127L || sum(yin_target$event) != 3442L || nrow(yang_target) != 1766L) {
  stop("Locked Yin/Yang target contract failed.", call. = FALSE)
}
bins <- rbindlist(list(
  yin_target[event == 1L, .(eid, relative_bin = baseline_bin(as.numeric(time_years)))],
  yang_target[, .(eid, relative_bin = baseline_bin(-as.numeric(disease_duration_years)))]
))
if (anchor == "diagnosis") bins[, relative_bin := -relative_bin]
bins[, side := if (anchor == "diagnosis") fifelse(relative_bin < 0, "Yin", "Yang") else fifelse(relative_bin < 0, "Yang", "Yin")]
counts <- bins[, .(n = .N), by = .(side, relative_bin)]
setorder(counts, relative_bin)

load_cached_scores <- function(methods, require_complete = FALSE) {
  rows <- list()
  input_paths <- character()
  origins <- list()
  strict_methods <- intersect(methods, c("pradeep-strict", "yu-strict"))
  old_input_root <- NULL
  if ("pradeep-strict" %in% strict_methods) {
    new_root <- file.path(yy_outdir, "score", "pradeep-strict")
    if (file.exists(file.path(new_root, "COMPLETE"))) {
      yin_file <- file.path(new_root, "scores_yin.csv.gz")
      yang_file <- file.path(new_root, "scores_yang.csv.gz")
      strict_yin <- fread(yin_file, colClasses = list(character = "eid"))[, .(eid, method, cohort_side, score_raw, score_z)]
      strict_yang <- fread(yang_file, colClasses = list(character = "eid"))[, .(eid, method, cohort_side, score_raw, score_z)]
      origins[["pradeep-strict"]] <- data.table(method = "pradeep-strict", origin = "completed frozen-project projection")
    } else {
      stop("Completed public score projection missing for pradeep-strict. Run: yy score pradeep-strict --project", call. = FALSE)
    }
    rows[["pradeep-strict-yin"]] <- strict_yin
    rows[["pradeep-strict-yang"]] <- strict_yang
    input_paths <- c(input_paths, yin_file, yang_file)
  }
  if ("yu-strict" %in% strict_methods) {
    new_root <- file.path(yy_outdir, "score", "yu-strict")
    if (file.exists(file.path(new_root, "COMPLETE"))) {
      yin_file <- file.path(new_root, "scores_yin.csv.gz")
      yang_file <- file.path(new_root, "scores_yang.csv.gz")
      strict_yin <- fread(yin_file, colClasses = list(character = "eid"))[, .(eid, method, cohort_side, score_raw, score_z)]
      strict_yang <- fread(yang_file, colClasses = list(character = "eid"))[, .(eid, method, cohort_side, score_raw, score_z)]
      origins[["yu-strict"]] <- data.table(method = "yu-strict", origin = "completed frozen-project projection")
    } else {
      stop("Completed public score projection missing for yu-strict. Run: yy score yu-strict --project", call. = FALSE)
    }
    rows[["yu-strict-yin"]] <- strict_yin
    rows[["yu-strict-yang"]] <- strict_yang
    input_paths <- c(input_paths, yin_file, yang_file)
  }
  fair_methods <- intersect(methods, c("pradeep-fair", "yu-fair"))
  for (method in fair_methods) {
    new_root <- file.path(yy_outdir, "score", method)
    if (file.exists(file.path(new_root, "COMPLETE"))) {
      yin_file <- file.path(new_root, "scores_yin_oof.csv.gz")
      yang_file <- file.path(new_root, "scores_yang.csv.gz")
      fair_yin <- fread(yin_file, colClasses = list(character = "eid"))[, .(eid, method, cohort_side = "Yin", score_raw, score_z)]
      fair_yang <- fread(yang_file, colClasses = list(character = "eid"))[, .(eid, method, cohort_side = "Yang", score_raw, score_z)]
    } else {
      stop("Completed public score projection missing for ", method,
           ". Run: yy score ", method, " --project", call. = FALSE)
    }
    if (file.exists(file.path(new_root, "COMPLETE"))) origins[[method]] <- data.table(method = method, origin = "completed common-fold model")
    rows[[paste0(method, "-yin")]] <- fair_yin
    rows[[paste0(method, "-yang")]] <- fair_yang
    input_paths <- c(input_paths, yin_file, yang_file)
  }
  scores <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (nrow(scores) != (37127L + 1766L) * length(methods) ||
      anyDuplicated(scores[, .(eid, method)]) || any(!is.finite(scores$score_z))) {
    stop("Cached individual score coverage failed.", call. = FALSE)
  }
  list(scores = scores, inputs = unique(input_paths), origins = rbindlist(origins, use.names = TRUE, fill = TRUE))
}

score_cache <- if (length(requested_methods)) load_cached_scores(requested_methods, require_score) else list(scores = data.table(), inputs = character(), origins = data.table())
if (recompute && length(requested_methods)) {
  source(score_helper, local = FALSE)
  message("--recompute validates frozen model parameters; cached plotting values remain the display contract.")
  strict <- intersect(requested_methods, c("pradeep-strict", "yu-strict"))
  fair <- intersect(requested_methods, c("pradeep-fair", "yu-fair"))
  if (length(strict)) {
    strict_labels <- c(`pradeep-strict` = "Pradeep local LASSO", `yu-strict` = "Yu local LightGBM")[strict]
    replay <- yy_legacy_score_bundle(strict_labels, "raw", anchor, show_roc, dir0, analysis_root, yy_outdir, script_root)
    score_cache$inputs <- unique(c(score_cache$inputs, replay$input_paths))
  }
  if (length(fair)) {
    fair_labels <- c(`pradeep-fair` = "Pradeep-style LASSO-logistic", `yu-fair` = "Yu-style LightGBM")[fair]
    replay <- yy_unified_score_bundle(fair_labels, "raw", anchor, show_roc, dir0, analysis_root, yy_outdir, script_root)
    score_cache$inputs <- unique(c(score_cache$inputs, replay$input_paths))
  }
}

trajectory_parts <- list()
roc_parts <- list()
metric_parts <- list()
input_paths <- c(dynamic_helper, dynamic_paths$core, dynamic_paths$yin_target, dynamic_paths$yang_target, score_cache$inputs)

if (length(requested_methods)) {
  score_cases <- score_cache$scores[, .(eid, method, value = score_z)]
  score_participants <- merge(bins, score_cases, by = "eid", allow.cartesian = TRUE, sort = FALSE)
  if (adjustment == "raw") {
    score_trajectory <- score_participants[, {
      n_value <- .N; mean_value <- mean(value); se_value <- sd(value) / sqrt(n_value)
      list(n = n_value, mean = mean_value, se = se_value,
           low = mean_value - 1.96 * se_value, high = mean_value + 1.96 * se_value)
    }, by = .(side, relative_bin, method)]
  } else {
    covariates <- yy_load_stepwise_covariates(dynamic_paths, adjustment_spec)
    input_paths <- c(input_paths, attr(covariates, "input_path"))
    score_participants <- merge(score_participants, covariates, by = "eid", all = FALSE, sort = FALSE)
    score_trajectory <- rbindlist(lapply(requested_methods, function(method_value) {
      marginal_adjusted_bins(
        score_participants[method == method_value], method_value, "Model score",
        adjustment = adjustment, rhs_terms = adjustment_spec$rhs
      )
    }))
  }
  if ("method" %in% names(score_trajectory)) {
    score_trajectory[, series := method_table$label[match(as.character(method), method_table$method)]]
    score_trajectory[, method := NULL]
  } else {
    score_trajectory[, series := method_table$label[match(as.character(series), method_table$method)]]
  }
  score_trajectory[, `:=`(series_type = "Model score", bin_n = n, displayed = !relative_bin %in% range(bins$relative_bin))]
  trajectory_parts[["scores"]] <- score_trajectory

  if (show_roc) {
    fair <- intersect(requested_methods, c("pradeep-fair", "yu-fair"))
    if (length(fair)) {
      yin_participant <- fread(dynamic_paths$participants_yin, colClasses = list(character = "eid"))
      for (method_value in fair) {
        one <- merge(
          yin_participant[, .(eid, outer_fold, time, event)],
          score_cache$scores[method == method_value & cohort_side == "Yin", .(eid, score = score_z)],
          by = "eid", all = FALSE, sort = FALSE
        )
        if (adjustment != "raw") {
          covariates <- yy_load_stepwise_covariates(dynamic_paths, adjustment_spec)
          one <- merge(one, covariates, by = "eid", all = FALSE, sort = FALSE)
          one[, score := cross_fitted_covariate_residual(.SD, adjustment_spec$rhs)]
        }
        result <- mean_fold_roc(one[, .(eid, outer_fold, time, event, score)], 5)
        label <- method_table[method == method_value, label]
        result$curve[, series := label]
        roc_parts[[method_value]] <- result$curve
        metric_parts[[method_value]] <- data.table(
          series = label, AUC_5y = result$mean_fold_auc,
          cohort = "Locked common Yin five-fold OOF cohort",
          model_definition = paste0("common five-year IPCW AUC; adjustment=", adjustment)
        )
      }
    }
    strict <- intersect(requested_methods, c("pradeep-strict", "yu-strict"))
    if (length(strict)) {
      strict_missing <- character()
      for (method_value in strict) {
        method_root <- file.path(yy_outdir, "score", method_value)
        if (file.exists(file.path(method_root, "COMPLETE"))) {
          curve_file <- file.path(method_root, "roc_native.csv.gz")
          metric_file <- file.path(method_root, "metrics.csv")
          curves <- fread(curve_file)
          metrics <- fread(metric_file)
          label <- method_table[method == method_value, label]
          curves[, series := label]
          metrics[, `:=`(
            series = label,
            cohort = if (method_value == "pradeep-strict") "Pradeep native held-out" else "Yu native held-out",
            model_definition = paste0("native eventual-event binary AUC; trajectory adjustment=", adjustment,
                                      " is not applied to the strict native ROC")
          )]
          roc_parts[[paste0("strict-", method_value)]] <- curves
          metric_parts[[paste0("strict-", method_value)]] <- metrics
          input_paths <- c(input_paths, curve_file, metric_file)
        } else {
          strict_missing <- c(strict_missing, method_value)
        }
      }
      if (length(strict_missing)) {
        stop("Completed public strict ROC projection missing for: ",
             paste(strict_missing, collapse = ", "),
             ". Run the corresponding yy score METHOD --project command.",
             call. = FALSE)
      }
    }
  }
}

selected_protein_labels <- character()
if (length(requested_proteins)) {
  protein_bundle <- yy_dynamic_protein_bundle(
    requested = requested_proteins, adjustment = adjustment, version = "unified",
    anchor = anchor, show_roc = show_roc,
    dir0 = dir0, analysis_root = analysis_root, yy_outdir = yy_outdir, script_root = script_root
  )
  trajectory_parts[["proteins"]] <- protein_bundle$trajectories
  selected_protein_labels <- protein_bundle$mapping$series
  input_paths <- c(input_paths, protein_bundle$input_paths)
  if (show_roc) {
    protein_bundle$metrics[, `:=`(
      cohort = "Locked common Yin five-fold OOF cohort",
      model_definition = paste0("single-protein Cox; five-year IPCW AUC; adjustment=", adjustment)
    )]
    roc_parts[["proteins"]] <- protein_bundle$roc
    metric_parts[["proteins"]] <- protein_bundle$metrics
  }
}

trajectories <- rbindlist(trajectory_parts, use.names = TRUE, fill = TRUE)
if (side_mode != "yy") {
  wanted <- if (side_mode == "yin") "Yin" else "Yang"
  trajectories <- trajectories[side == wanted]
  counts <- counts[side == wanted]
}
if ("displayed" %in% names(trajectories)) trajectories <- trajectories[is.na(displayed) | displayed]
if (!nrow(trajectories)) stop("No trajectory rows remain.", call. = FALSE)
method_labels <- method_table[method %in% requested_methods, label]
selected_series <- c(method_labels, selected_protein_labels)
trajectories[, series := factor(as.character(series), levels = selected_series)]
if (show_sd || show_sd_line) trajectories[, `:=`(sd = se * sqrt(n), sd_low = mean - se * sqrt(n), sd_high = mean + se * sqrt(n))]
if (show_mean || show_sd) {
  trajectories[, point_label := if (show_mean && show_sd) sprintf("mean %.2f; SD %.2f", mean, sd) else if (show_mean) sprintf("%.2f", mean) else sprintf("SD %.2f", sd)]
}

palette <- c(
  "Pradeep strict" = "#0B3C5D", "Yu strict" = "#D55E00",
  "Pradeep fair" = "#4C78A8", "Yu fair" = "#F2A541"
)
if (length(selected_protein_labels)) {
  protein_colours <- grDevices::hcl(seq(15, 375, length.out = length(selected_protein_labels) + 1L)[seq_along(selected_protein_labels)], c = 90, l = 48)
  names(protein_colours) <- selected_protein_labels
  palette <- c(palette, protein_colours)
}
palette <- palette[selected_series]

make_trajectory <- function() {
  z_values <- c(0, trajectories$mean)
  if (show_sd_line) z_values <- c(z_values, trajectories$sd_low, trajectories$sd_high)
  z_range <- range(z_values, finite = TRUE); z_pad <- max(0.12, diff(z_range) * 0.08)
  z_limits <- z_range + c(-z_pad, z_pad)
  if (show_bar) {
    count_max <- ceiling(max(counts$n) * 1.15 / 50) * 50
    z_scale <- count_max / diff(z_limits)
    transform <- function(x) (x - z_limits[[1L]]) * z_scale
    trajectories[, plot_y := transform(mean)]
    if (show_sd_line) trajectories[, `:=`(plot_low = transform(sd_low), plot_high = transform(sd_high))]
    plot <- ggplot() +
      geom_col(data = counts, aes(relative_bin, n, fill = side), width = 1.8, alpha = 0.58, colour = "white") +
      geom_text(data = counts, aes(relative_bin, n, label = n), vjust = 1.35,
        family = "Helvetica", size = 2.05, colour = "#59636C") +
      scale_fill_manual(values = c(Yin = "#B8CBE2", Yang = "#E5B3AA"), guide = "none") +
      scale_y_continuous(name = "Patient count", limits = c(0, count_max),
        sec.axis = sec_axis(~ . / z_scale + z_limits[[1L]], name = "Standardized level (z; Yin reference)"),
        expand = expansion(mult = c(0, 0))) +
      geom_hline(yintercept = transform(0), linewidth = 0.35, linetype = 2, colour = "#5A6168")
  } else {
    trajectories[, plot_y := mean]
    if (show_sd_line) trajectories[, `:=`(plot_low = sd_low, plot_high = sd_high)]
    plot <- ggplot() + scale_y_continuous(name = "Standardized level (z; Yin reference)") +
      geom_hline(yintercept = 0, linewidth = 0.35, linetype = 2, colour = "#5A6168")
  }
  if (show_sd_line) plot <- plot + geom_errorbar(data = trajectories,
    aes(relative_bin, ymin = plot_low, ymax = plot_high, colour = series), width = 0.25, linewidth = 0.55)
  protein_data <- trajectories[!as.character(series) %in% method_labels]
  score_data <- trajectories[as.character(series) %in% method_labels]
  if (nrow(protein_data)) plot <- plot +
    geom_line(data = protein_data, aes(relative_bin, plot_y, colour = series, group = series), linewidth = 1.05) +
    geom_point(data = protein_data, aes(relative_bin, plot_y, colour = series), size = 3)
  if (nrow(score_data)) plot <- plot +
    geom_line(data = score_data, aes(relative_bin, plot_y, colour = series, group = series), linewidth = 1.9) +
    geom_point(data = score_data, aes(relative_bin, plot_y, colour = series), size = 4, shape = 21, fill = "white", stroke = 1.1)
  if (show_mean || show_sd) plot <- plot + geom_text(data = trajectories,
    aes(relative_bin, plot_y, label = point_label, colour = series), vjust = -0.85, size = 2.45, check_overlap = TRUE)
  if (min(trajectories$relative_bin) <= 0 && max(trajectories$relative_bin) >= 0) {
    plot <- plot + geom_vline(xintercept = 0, linewidth = 0.65, linetype = 3, colour = "#4C4C4C")
  }
  plot + scale_colour_manual(values = palette, breaks = selected_series, drop = FALSE, name = NULL) +
    labs(x = if (anchor == "baseline") "CAD diagnosis timing relative to baseline (years)" else "Years relative to CAD diagnosis") +
    theme_classic(base_size = 11.5, base_family = "Helvetica") +
    theme(axis.title = element_text(face = "bold"), legend.position = "bottom", legend.text = element_text(size = 8.3, face = "bold"))
}

make_roc <- function() {
  curves <- rbindlist(roc_parts, use.names = TRUE, fill = TRUE)
  metrics <- rbindlist(metric_parts, use.names = TRUE, fill = TRUE)
  if (!nrow(curves) || !nrow(metrics)) stop("No ROC results were resolved.", call. = FALSE)
  metrics[, metric_value := NA_real_]
  if ("AUC_5y" %in% names(metrics)) metrics[is.finite(AUC_5y), metric_value := AUC_5y]
  if ("auc" %in% names(metrics)) metrics[!is.finite(metric_value) & is.finite(auc), metric_value := auc]
  metrics[, metric_label := fifelse(grepl("common|five-year", paste(cohort, model_definition), ignore.case = TRUE), "AUC5", "AUC")]
  labels <- metrics[, .(series, legend = sprintf("%s  %s %.3f", series, metric_label, metric_value))]
  labels <- unique(labels, by = "series")
  curves <- merge(curves, labels, by = "series", all.x = TRUE)
  order <- labels[match(selected_series, series), legend]; order <- order[!is.na(order)]
  colours <- setNames(unname(palette[labels$series]), labels$legend)
  curves[, legend := factor(legend, levels = order)]
  mixed <- any(requested_methods %in% c("pradeep-strict", "yu-strict")) &&
    (length(requested_proteins) || any(requested_methods %in% c("pradeep-fair", "yu-fair")))
  plot <- ggplot(curves, aes(false_positive_rate, true_positive_rate, colour = legend)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") + geom_line(linewidth = 1.1) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    scale_colour_manual(values = colours, breaks = order, drop = FALSE, name = NULL) +
    guides(colour = guide_legend(ncol = 2, byrow = TRUE)) +
    labs(x = "1 - Specificity", y = "Sensitivity",
      title = if (mixed) "ROC audit across native and common estimands" else "Incident CAD ROC",
      caption = if (mixed) "Strict native binary AUC and common five-year IPCW AUC are different estimands; visual juxtaposition is not a paired comparison." else NULL) +
    theme_classic(base_size = 11.5, base_family = "Helvetica") +
    theme(axis.title = element_text(face = "bold"), legend.position = "bottom", legend.box = "vertical",
      legend.text = element_text(size = 7.2, face = "bold"), plot.caption = element_text(size = 6.8, colour = "#7A2E2E"), aspect.ratio = 1)
  list(plot = plot, curves = curves, metrics = metrics)
}

trajectory_plot <- if (show_traj) make_trajectory() else NULL
roc_result <- if (show_roc) make_roc() else NULL
figure_mode <- if (show_traj && show_roc) "trajectory_roc_combined" else if (show_traj) "trajectory" else "roc"
final_plot <- if (show_traj && show_roc) {
  trajectory_plot + roc_result$plot + plot_layout(widths = c(1.35, 1)) + plot_annotation(tag_levels = "a")
} else if (show_traj) trajectory_plot else roc_result$plot

parameter_text <- paste(side_mode, anchor, adjustment, figure_mode,
  paste(requested_methods, collapse = "|"), paste(selected_protein_labels, collapse = "|"), sep = ";")
run_name <- paste0("yyplot_", slug(side_mode), "_", slug(anchor), "_", slug(adjustment), "_", figure_mode, "_", hash_text(parameter_text))
run_dir <- file.path(output_root, "custom", run_name)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
stem <- file.path(run_dir, run_name)
width <- if (show_traj && show_roc) 13.2 else if (show_roc) 8.5 else 11.5
height <- if (show_traj && show_roc) 6.4 else 8.5
ggsave(paste0(stem, ".png"), final_plot, width = width, height = height, dpi = 300, device = ragg::agg_png, bg = "white")
ggsave(paste0(stem, ".pdf"), final_plot, width = width, height = height, device = cairo_pdf)
ggsave(paste0(stem, ".svg"), final_plot, width = width, height = height, device = svglite::svglite)
if (show_traj) {
  fwrite(trajectories, paste0(stem, ".trajectory_source.csv"))
  fwrite(counts, paste0(stem, ".count_source.csv"))
}
if (show_roc) {
  fwrite(roc_result$curves, paste0(stem, ".roc_source.csv"))
  fwrite(roc_result$metrics, paste0(stem, ".roc_metrics.csv"))
}
fwrite(method_table[method %in% requested_methods], paste0(stem, ".score_methods.csv"))
parameters <- data.table(
  parameter = c("side", "anchor", "adjustment", "adjustment_formula", "proteins", "score_methods", "score_origins", "require_score", "recompute", "roc_warning"),
  value = c(side_mode, anchor, adjustment, adjustment_spec$rhs,
    paste(selected_protein_labels, collapse = "|"), paste(requested_methods, collapse = "|"),
    paste(score_cache$origins[, paste(method, origin, sep = "=")], collapse = "|"), require_score, recompute,
    if (show_roc && any(requested_methods %in% c("pradeep-strict", "yu-strict")))
      "strict native binary AUC is not the same estimand as common five-year IPCW AUC5" else "none")
)
fwrite(parameters, paste0(stem, ".parameters.csv"))
input_paths <- unique(input_paths[file.exists(input_paths)])
fwrite(data.table(path = normalizePath(input_paths, winslash = "/", mustWork = TRUE),
  bytes = file.info(input_paths)$size, md5 = unname(tools::md5sum(input_paths))), paste0(stem, ".input_manifest.csv"))
writeLines(c("status\tadjustment\tcompleted_at",
  paste("COMPLETE", adjustment, format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), sep = "\t")), file.path(run_dir, "run_status.tsv"))
file.create(file.path(run_dir, "COMPLETE"))
cat("YY_PLOT_COMPLETE\nOUTPUT_DIR=", run_dir, "\nFIGURE_PNG=", paste0(stem, ".png"),
    "\nADJUSTMENT=", adjustment, "\nSERIES=", paste(selected_series, collapse = " | "), "\n", sep = "")
