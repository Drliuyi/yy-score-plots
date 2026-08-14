# Generate figures from the local UKB-PPP cardiac reproduction outputs.

locate_repro_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")))
  frames <- sys.frames()
  ofiles <- vapply(frames, function(frame) if (!is.null(frame$ofile)) frame$ofile else NA_character_, character(1))
  ofiles <- ofiles[!is.na(ofiles)]
  if (length(ofiles) > 0) return(dirname(normalizePath(ofiles[length(ofiles)], winslash = "/")))
  project_dir <- Sys.getenv("PRADEEP_PROJECT_DIR", unset = "")
  candidates <- if (nzchar(project_dir)) file.path(project_dir, "f") else getwd()
  hit <- candidates[file.exists(file.path(candidates, "00_config.R"))]
  if (length(hit) > 0) return(normalizePath(hit[1], winslash = "/"))
  stop("Cannot locate pradeep/f. Run this step through pradeep.sh.", call. = FALSE)
}

this_dir <- locate_repro_script_dir()
source(file.path(this_dir, "00_config.R"))
source(file.path(script_dir, "R", "repro_utils.R"))
repro_load_packages(c("data.table", "ggplot2", "ggrepel", "patchwork", "survival"))

message("Step 05: making figures")

fig_dir <- file.path(output_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
make_legacy_figures <- repro_bool_env("UKBPPP_MAKE_LEGACY_FIGURES", default = FALSE)
legacy_figure_files <- c(
  "primary_association_volcano.png",
  "primary_association_significant_counts.png",
  "primary_association_summary.csv",
  "primary_top25_forest_af.png",
  "primary_top25_forest_as.png",
  "primary_top25_forest_cad.png",
  "primary_top25_forest_hf.png",
  "sex_interaction_volcano.png",
  "sex_interaction_top20_forest.png",
  "sex_interaction_summary.csv",
  "lasso_auc_by_outcome.png",
  "lasso_auc_summary_for_plot.csv",
  "lasso_roc_curves.png"
)
if (!make_legacy_figures) {
  unlink(file.path(fig_dir, legacy_figure_files), force = TRUE)
}

primary_path <- file.path(output_dir, "primary_association_cox.csv")
sex_path <- file.path(output_dir, "sex_interaction_cox.csv")
auc_path <- file.path(output_dir, "lasso_risk_score_auc.csv")
if (!file.exists(primary_path)) stop("Missing primary association results: ", primary_path, call. = FALSE)
if (!file.exists(sex_path)) stop("Missing sex interaction results: ", sex_path, call. = FALSE)
if (!file.exists(auc_path)) stop("Missing LASSO AUC results: ", auc_path, call. = FALSE)

safe_neglog10 <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  p[!is.finite(p) | p <= 0] <- .Machine$double.xmin
  -log10(p)
}

outcome_levels <- c("CAD", "Afib", "HF", "AS")
outcome_labels <- c(CAD = "CAD", Afib = "AF", HF = "HF", AS = "AS")
outcome_full_labels <- c(
  CAD = "Coronary artery disease",
  Afib = "Atrial fibrillation",
  HF = "Heart failure",
  AS = "Aortic stenosis"
)
lasso_outcome_full_labels <- c(
  cad = "Coronary artery disease",
  afib = "Atrial fibrillation",
  hfail = "Heart failure",
  ao_sten = "Aortic stenosis"
)
model_levels <- c("Clinical", "Proteins", "Combination")
model_cols <- c(Clinical = "#4E79A7", Proteins = "#59A14F", Combination = "#E15759")
sig_cols <- c("Not significant" = "grey78", "FDR < 0.05" = "#F28E2B", "Bonferroni" = "#D62728")

save_plot <- function(plot, filename, width, height) {
  if (!make_legacy_figures && !startsWith(filename, "original_style_")) {
    return(invisible(NA_character_))
  }
  path <- file.path(fig_dir, filename)
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = 320, bg = "white")
  invisible(path)
}

primary <- data.table::fread(primary_path)
repro_assert_primary_complete(primary, audit_dir, outcome_map, context = "Primary association results used for figures")
primary <- primary[is.finite(P_Value) & is.finite(beta)]
primary[, Outcome := factor(Outcome, levels = outcome_levels)]
primary[, outcome_label := factor(outcome_labels[as.character(Outcome)], levels = outcome_labels[outcome_levels])]
primary[, neglog10p := safe_neglog10(P_Value)]
primary[, significance := data.table::fifelse(
  Bonferroni %in% TRUE, "Bonferroni",
  data.table::fifelse(FDR < 0.05, "FDR < 0.05", "Not significant")
)]
primary[, significance := factor(significance, levels = names(sig_cols))]
primary[, logHR := beta]
n_proteins <- length(unique(primary$Protein))
primary_bonferroni_threshold <- if ("Bonferroni_Threshold" %in% names(primary)) {
  unique(primary$Bonferroni_Threshold[is.finite(primary$Bonferroni_Threshold)])[1]
} else {
  0.05 / (n_proteins * length(outcome_levels))
}
bonf_line <- -log10(primary_bonferroni_threshold)

primary_labels <- primary[order(P_Value), head(.SD, 8), by = Outcome]
primary_counts <- primary[, .(
  n_proteins = .N,
  n_fdr = sum(FDR < 0.05, na.rm = TRUE),
  n_bonferroni = sum(Bonferroni %in% TRUE, na.rm = TRUE),
  top_protein = Protein[which.min(P_Value)],
  top_p = min(P_Value, na.rm = TRUE)
), by = .(Outcome, outcome_label)]
if (make_legacy_figures) {
  data.table::fwrite(primary_counts, file.path(fig_dir, "primary_association_summary.csv"))
}

p_primary_volcano <- ggplot2::ggplot(primary, ggplot2::aes(logHR, neglog10p, color = significance)) +
  ggplot2::geom_hline(yintercept = bonf_line, linewidth = 0.3, linetype = "dashed", color = "grey35") +
  ggplot2::geom_point(size = 1.15, alpha = 0.78) +
  ggrepel::geom_text_repel(
    data = primary_labels,
    ggplot2::aes(label = Protein),
    size = 2.4,
    max.overlaps = Inf,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~ outcome_label, scales = "free_y", nrow = 2) +
  ggplot2::scale_color_manual(values = sig_cols, drop = FALSE) +
  ggplot2::labs(
    title = "Primary protein associations with incident cardiac diseases",
    subtitle = paste0("Adjusted time-varying Cox models; dashed line = Bonferroni threshold across ", n_proteins, " proteins"),
    x = "Log hazard ratio per 1 SD protein",
    y = expression(-log[10](P)),
    color = NULL
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
save_plot(p_primary_volcano, "primary_association_volcano.png", 11, 8)

primary_count_long <- data.table::melt(
  primary_counts,
  id.vars = c("Outcome", "outcome_label"),
  measure.vars = c("n_fdr", "n_bonferroni"),
  variable.name = "threshold",
  value.name = "n"
)
primary_count_long[, threshold := factor(threshold, levels = c("n_fdr", "n_bonferroni"), labels = c("FDR < 0.05", "Bonferroni"))]
p_primary_counts <- ggplot2::ggplot(primary_count_long, ggplot2::aes(outcome_label, n, fill = threshold)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.65) +
  ggplot2::geom_text(ggplot2::aes(label = n), position = ggplot2::position_dodge(width = 0.75), vjust = -0.25, size = 3) +
  ggplot2::scale_fill_manual(values = c("FDR < 0.05" = "#F28E2B", Bonferroni = "#D62728")) +
  ggplot2::labs(title = "Number of significant protein associations", x = NULL, y = "Proteins", fill = NULL) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
save_plot(p_primary_counts, "primary_association_significant_counts.png", 7, 5)

for (outc in outcome_levels) {
  d <- primary[Outcome == outc][order(P_Value)][seq_len(min(.N, 25))]
  d[, Protein_plot := factor(Protein, levels = rev(Protein))]
  p <- ggplot2::ggplot(d, ggplot2::aes(HR, Protein_plot)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.35) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_Lower, xmax = CI_Upper, color = significance), height = 0.18, linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(color = significance), size = 2) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(values = sig_cols, drop = FALSE) +
    ggplot2::labs(
      title = paste0("Top protein associations: ", outcome_labels[[outc]]),
      x = "Hazard ratio per 1 SD protein",
      y = NULL,
      color = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
  save_plot(p, paste0("primary_top25_forest_", tolower(as.character(outcome_labels[[outc]])), ".png"), 7.2, 7.5)
}

sex <- data.table::fread(sex_path)
sex <- sex[is.finite(P_Val_Int) & is.finite(HR_Int)]
sex[, Outcome := factor(Outcome, levels = outcome_levels)]
sex[, outcome_label := factor(outcome_labels[as.character(Outcome)], levels = outcome_labels[outcome_levels])]
sex[, logHR_Int := log(HR_Int)]
sex[, neglog10p := safe_neglog10(P_Val_Int)]

# Paper rule: first define protein-outcome candidates that are Bonferroni
# significant in either sex, then Bonferroni-correct interaction P values only
# across those candidates. Recompute here as a guard for older result files.
n_sex_tests <- sum(is.finite(sex$P_Val_M) | is.finite(sex$P_Val_F))
sex_association_threshold <- if (n_sex_tests > 0L) 0.05 / n_sex_tests else NA_real_
sex[, Paper_Candidate := (
  is.finite(P_Val_M) & P_Val_M < sex_association_threshold
) | (
  is.finite(P_Val_F) & P_Val_F < sex_association_threshold
)]
n_paper_candidates <- sex[Paper_Candidate %in% TRUE, .N]
paper_interaction_threshold <- if (n_paper_candidates > 0L) 0.05 / n_paper_candidates else NA_real_
sex[, Paper_Interaction_Threshold := paper_interaction_threshold]
sex[, Paper_Interaction_P_Adjusted := NA_real_]
if (n_paper_candidates > 0L) {
  sex[Paper_Candidate %in% TRUE, Paper_Interaction_P_Adjusted := p.adjust(
    P_Val_Int,
    method = "bonferroni",
    n = n_paper_candidates
  )]
}
sex[, Paper_Interaction_Bonferroni := Paper_Candidate %in% TRUE &
  is.finite(P_Val_Int) & P_Val_Int < paper_interaction_threshold]
sex[, significance := data.table::fifelse(
  Paper_Interaction_Bonferroni %in% TRUE,
  "Bonferroni interaction",
  data.table::fifelse(Paper_Candidate %in% TRUE, "Paper candidate", "Not a paper candidate")
)]
sex[, significance := factor(
  significance,
  levels = c("Not a paper candidate", "Paper candidate", "Bonferroni interaction")
)]
sex_labels <- sex[order(-as.integer(Paper_Interaction_Bonferroni), P_Val_Int), head(.SD, 8), by = Outcome]
sex_counts <- sex[, .(
  n_candidates = sum(Paper_Candidate, na.rm = TRUE),
  n_bonferroni_interactions = sum(Paper_Interaction_Bonferroni, na.rm = TRUE),
  top_protein = Protein[which.min(P_Val_Int)],
  top_p = min(P_Val_Int, na.rm = TRUE),
  sex_association_threshold = sex_association_threshold,
  paper_interaction_threshold = paper_interaction_threshold
), by = .(Outcome, outcome_label)]
data.table::fwrite(sex_counts, file.path(fig_dir, "sex_interaction_summary.csv"))
data.table::fwrite(
  sex[Paper_Interaction_Bonferroni %in% TRUE][order(P_Val_Int)],
  file.path(fig_dir, "sex_interaction_paper_significant_source_data.csv")
)

p_sex_volcano <- ggplot2::ggplot(sex, ggplot2::aes(logHR_Int, neglog10p, color = significance)) +
  ggplot2::geom_point(size = 1.15, alpha = 0.78) +
  ggplot2::geom_hline(
    yintercept = -log10(paper_interaction_threshold),
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.35
  ) +
  ggrepel::geom_text_repel(
    data = sex_labels,
    ggplot2::aes(label = Protein),
    size = 2.4,
    max.overlaps = Inf,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~ outcome_label, scales = "free_y", nrow = 2) +
  ggplot2::scale_color_manual(values = c(
    "Not a paper candidate" = "grey82",
    "Paper candidate" = "#80B1D3",
    "Bonferroni interaction" = "#7B3294"
  ), drop = FALSE) +
  ggplot2::labs(
    title = "Protein-by-sex interaction associations",
    x = "Log interaction HR",
    y = expression(-log[10](P[interaction])),
    color = NULL
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
save_plot(p_sex_volcano, "sex_interaction_volcano.png", 11, 8)

sex_top <- sex[order(-as.integer(Paper_Interaction_Bonferroni), -as.integer(Paper_Candidate), P_Val_Int), head(.SD, 20), by = Outcome]
sex_top[, Protein_plot := factor(paste(outcome_label, Protein, sep = " - "), levels = rev(paste(outcome_label, Protein, sep = " - ")))]
p_sex_forest <- ggplot2::ggplot(sex_top, ggplot2::aes(HR_Int, Protein_plot, color = significance)) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_Low_Int, xmax = CI_High_Int), height = 0.18, linewidth = 0.5) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_color_manual(values = c(
    "Not a paper candidate" = "grey65",
    "Paper candidate" = "#80B1D3",
    "Bonferroni interaction" = "#7B3294"
  ), drop = FALSE) +
  ggplot2::labs(
    title = "Top protein-by-sex interaction estimates",
    x = "Interaction HR",
    y = NULL,
    color = NULL
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
save_plot(p_sex_forest, "sex_interaction_top20_forest.png", 8.2, 12)

auc <- data.table::fread(auc_path)
auc[, outc_label := factor(c(cad = "CAD", afib = "AF", hfail = "HF", ao_sten = "AS")[outc], levels = c("CAD", "AF", "HF", "AS"))]
auc[, model := factor(model, levels = model_levels)]
if (make_legacy_figures) {
  data.table::fwrite(auc, file.path(fig_dir, "lasso_auc_summary_for_plot.csv"))
}
p_auc <- ggplot2::ggplot(auc, ggplot2::aes(outc_label, auc, color = model, group = model)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = auc_lci, ymax = auc_uci), width = 0.14, position = ggplot2::position_dodge(width = 0.55), linewidth = 0.55) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.55), size = 2.4) +
  ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.55), linewidth = 0.45) +
  ggplot2::scale_color_manual(values = model_cols, drop = FALSE) +
  ggplot2::coord_cartesian(ylim = c(max(0.5, min(auc$auc_lci, na.rm = TRUE) - 0.03), min(0.95, max(auc$auc_uci, na.rm = TRUE) + 0.03))) +
  ggplot2::labs(
    title = "LASSO risk-score discrimination",
    x = NULL,
    y = "AUC with 95% CI",
    color = NULL
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
save_plot(p_auc, "lasso_auc_by_outcome.png", 8, 5.5)

roc_files <- list.files(file.path(output_dir, "lasso"), pattern = "_roc_rawdata\\.csv$", full.names = TRUE)
if (length(roc_files) > 0) {
  roc <- data.table::rbindlist(lapply(roc_files, data.table::fread), fill = TRUE)
  roc[, outc_label := factor(c(cad = "CAD", afib = "AF", hfail = "HF", ao_sten = "AS")[outc], levels = c("CAD", "AF", "HF", "AS"))]
  roc[, model := factor(model, levels = model_levels)]
  auc_lab <- auc[, .(outc_label, model, auc_label = sprintf("%s %.3f", model, auc))]
  roc <- merge(roc, auc_lab, by = c("outc_label", "model"), all.x = TRUE)
  p_roc <- ggplot2::ggplot(roc, ggplot2::aes(FPR, TPR, color = model)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.35) +
    ggplot2::geom_line(linewidth = 0.65) +
    ggplot2::facet_wrap(~ outc_label, nrow = 2) +
    ggplot2::scale_color_manual(values = model_cols, drop = FALSE) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::labs(title = "LASSO ROC curves in held-out test sets", x = "False positive rate", y = "True positive rate", color = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
  save_plot(p_roc, "lasso_roc_curves.png", 8, 7)
}

metrics_path <- file.path(output_dir, "lasso_prediction_metrics.csv")
if (file.exists(metrics_path)) {
  metrics <- data.table::fread(metrics_path)
  metrics[, outc_label := factor(c(cad = "CAD", afib = "AF", hfail = "HF", ao_sten = "AS")[outc], levels = c("CAD", "AF", "HF", "AS"))]
  metrics[, model := factor(model, levels = model_levels)]
  metric_cols <- intersect(
    c("auc", "sensitivity", "specificity", "ppv", "npv", "accuracy", "f1", "brier", "calibration_slope"),
    names(metrics)
  )
  metric_long <- data.table::melt(
    metrics,
    id.vars = c("outc_label", "model"),
    measure.vars = metric_cols,
    variable.name = "metric",
    value.name = "value"
  )
  metric_long[, metric_label := factor(
    metric,
    levels = metric_cols,
    labels = c(
      auc = "AUC",
      sensitivity = "Sensitivity",
      specificity = "Specificity",
      ppv = "PPV",
      npv = "NPV",
      accuracy = "Accuracy",
      f1 = "F1",
      brier = "Brier",
      calibration_slope = "Calibration slope"
    )[metric_cols]
  )]
  metric_long[, value_label := ifelse(is.finite(value), sprintf("%.2f", value), "NA")]
  p_metrics <- ggplot2::ggplot(metric_long, ggplot2::aes(model, outc_label, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = value_label), size = 2.7) +
    ggplot2::facet_wrap(~ metric_label, ncol = 3, scales = "free") +
    ggplot2::scale_fill_gradient(low = "white", high = "#3C5488", na.value = "grey95", name = "Value") +
    ggplot2::labs(title = "Held-out test-set prediction metrics", x = NULL, y = NULL) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA)
    )
  save_plot(p_metrics, "original_style_lasso_prediction_metrics.png", 10.5, 8.0)
}

calibration_path <- file.path(output_dir, "lasso_calibration.csv")
if (file.exists(calibration_path)) {
  calibration <- data.table::fread(calibration_path)
  calibration[, outc_label := factor(c(cad = "CAD", afib = "AF", hfail = "HF", ao_sten = "AS")[outc], levels = c("CAD", "AF", "HF", "AS"))]
  calibration[, model := factor(model, levels = model_levels)]
  cal_max <- max(calibration$observed_rate, calibration$mean_predicted_risk, na.rm = TRUE)
  cal_max <- max(0.05, min(1, cal_max * 1.15))
  p_calibration <- ggplot2::ggplot(calibration, ggplot2::aes(mean_predicted_risk, observed_rate, color = model, group = model)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.35) +
    ggplot2::geom_line(linewidth = 0.55) +
    ggplot2::geom_point(ggplot2::aes(size = n), alpha = 0.9) +
    ggplot2::facet_wrap(~ outc_label, nrow = 2) +
    ggplot2::scale_color_manual(values = model_cols, drop = FALSE, name = "") +
    ggplot2::scale_size_continuous(range = c(1.4, 3.2), name = "Bin n") +
    ggplot2::guides(size = "none") +
    ggplot2::coord_equal(xlim = c(0, cal_max), ylim = c(0, cal_max), expand = FALSE) +
    ggplot2::labs(title = "Held-out test-set calibration by risk decile", x = "Mean predicted risk", y = "Observed event rate") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")
  save_plot(p_calibration, "original_style_lasso_calibration.png", 8.5, 7.0)
}

# Paper-style figures adapted from the original 8_figures_tables.R script.
# These keep the original visual grammar while using the local result files.
# The paper caption defines the Miami/Manhattan colored points as
# Bonferroni-corrected P < 0.05 across all protein-outcome tests.
paper_primary_sig <- primary_bonferroni_threshold
paper_primary_line <- -log10(paper_primary_sig)
paper_chr_levels <- c(as.character(1:22), "X")
paper_chr_cols <- c(
  stats::setNames(rep(c("grey90", "grey80"), length.out = length(paper_chr_levels)), paper_chr_levels),
  "Associated with more than one outcome" = "#3C5488",
  "Associated with one outcome" = "#6BB56B"
)

require_paper_figures <- repro_bool_env("UKBPPP_REQUIRE_PAPER_FIGURES", default = TRUE)
paper_density_spec <- data.table::data.table(
  outcome = c("CAD", "HF", "Afib", "AS"),
  protein = c("GDF15", "WFDC2", "NTproBNP", "GDF15"),
  protein_label = c("GDF15", "WFDC2", "NT-proBNP", "GDF15"),
  event_col = c("cad_inc", "hfail_inc", "afib_inc", "ao_sten_inc"),
  inset_left = c(0.51, 0.47, 0.19, 0.68),
  inset_bottom = c(0.64, 0.64, 0.64, 0.64),
  inset_right = c(0.74, 0.73, 0.47, 0.95),
  inset_top = c(0.98, 0.98, 0.98, 0.98)
)
paper_density_cols <- c(Controls = "#9EAAC4", Cases = "#3C5488")

if (!file.exists(analysis_base_file)) {
  stop(
    "PAPER_FIGURE_GATE_FAILED: missing analysis base for primary density insets: ",
    analysis_base_file,
    call. = FALSE
  )
}
base_for_density_insets <- readRDS(analysis_base_file)
if (!is.list(base_for_density_insets) || is.null(base_for_density_insets$dat)) {
  stop("PAPER_FIGURE_GATE_FAILED: analysis base does not contain $dat.", call. = FALSE)
}
density_base <- data.table::as.data.table(base_for_density_insets$dat)
ntprobnp_candidates <- c("NTproBNP", "NPPB", "NTPROBNP")
ntprobnp_column <- ntprobnp_candidates[ntprobnp_candidates %in% names(density_base)][1]
if (!is.na(ntprobnp_column)) {
  paper_density_spec[outcome == "Afib", protein := ntprobnp_column]
}

density_qc_rows <- vector("list", nrow(paper_density_spec))
density_source_rows <- list()
for (i in seq_len(nrow(paper_density_spec))) {
  spec <- paper_density_spec[i]
  required_cols <- c(spec$protein, spec$event_col)
  columns_present <- all(required_cols %in% names(density_base))
  d <- data.table::data.table(value = numeric(), event = integer())
  if (columns_present) {
    d <- density_base[, .(
      value = suppressWarnings(as.numeric(get(spec$protein))),
      event = suppressWarnings(as.integer(get(spec$event_col)))
    )]
    d <- d[is.finite(value) & event %in% c(0L, 1L)]
  }
  n_controls <- sum(d$event == 0L)
  n_cases <- sum(d$event == 1L)
  controls_variable <- n_controls >= 2L && data.table::uniqueN(d[event == 0L, value]) >= 2L
  cases_variable <- n_cases >= 2L && data.table::uniqueN(d[event == 1L, value]) >= 2L
  status <- if (columns_present && controls_variable && cases_variable) "PASS" else "FAIL"
  density_qc_rows[[i]] <- data.table::data.table(
    outcome = spec$outcome,
    protein = spec$protein,
    event_col = spec$event_col,
    columns_present = columns_present,
    n_controls = n_controls,
    n_cases = n_cases,
    controls_variable = controls_variable,
    cases_variable = cases_variable,
    status = status
  )

  if (identical(status, "PASS")) {
    common_range <- range(d$value, finite = TRUE)
    for (event_value in c(0L, 1L)) {
      group_label <- if (event_value == 0L) "Controls" else "Cases"
      density_fit <- stats::density(
        d[event == event_value, value],
        n = 512L,
        from = common_range[1],
        to = common_range[2],
        na.rm = TRUE
      )
      density_source_rows[[length(density_source_rows) + 1L]] <- data.table::data.table(
        outcome = spec$outcome,
        protein = spec$protein,
        protein_label = spec$protein_label,
        event_col = spec$event_col,
        group = group_label,
        n = if (event_value == 0L) n_controls else n_cases,
        x = density_fit$x,
        density = density_fit$y
      )
    }
  }
}
paper_density_qc <- data.table::rbindlist(density_qc_rows, fill = TRUE)
paper_density_source <- data.table::rbindlist(density_source_rows, fill = TRUE)
if (ncol(paper_density_source) == 0L) {
  paper_density_source <- data.table::data.table(
    outcome = character(), protein = character(), protein_label = character(),
    event_col = character(), group = character(), n = integer(),
    x = numeric(), density = numeric()
  )
}
data.table::fwrite(
  paper_density_qc,
  file.path(fig_dir, "original_style_primary_density_inset_qc.csv")
)
data.table::fwrite(
  paper_density_source,
  file.path(fig_dir, "original_style_primary_density_inset_source_data.csv")
)
if (require_paper_figures && any(paper_density_qc$status != "PASS")) {
  failed_density <- paper_density_qc[status != "PASS", paste(outcome, protein, event_col, sep = "/")]
  stop(
    "PAPER_FIGURE_GATE_FAILED: primary density inset QC failed for ",
    paste(failed_density, collapse = ", "),
    call. = FALSE
  )
}

make_primary_density_plot <- function(outc, compact = FALSE) {
  spec <- paper_density_spec[outcome == outc]
  d <- copy(paper_density_source[outcome == outc])
  d[, group := factor(group, levels = c("Controls", "Cases"))]
  p <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = density, fill = group, color = group)) +
    ggplot2::geom_area(position = "identity", alpha = 0.4, color = NA) +
    ggplot2::geom_line(linewidth = if (compact) 0.28 else 0.45) +
    ggplot2::scale_fill_manual(values = paper_density_cols, drop = FALSE, name = NULL) +
    ggplot2::scale_color_manual(values = paper_density_cols, drop = FALSE, guide = "none") +
    ggplot2::labs(
      title = if (compact) NULL else outcome_full_labels[[outc]],
      x = paste0("Circulating ", spec$protein_label, " (SD)"),
      y = "Density"
    ) +
    ggplot2::theme_classic(base_size = if (compact) 5.2 else 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
      legend.position = if (compact) c(0.76, 0.79) else "bottom",
      legend.justification = c(0, 0),
      legend.direction = "vertical",
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key.height = grid::unit(if (compact) 3.0 else 8.0, "pt"),
      legend.key.width = grid::unit(if (compact) 4.0 else 11.0, "pt"),
      legend.text = ggplot2::element_text(size = if (compact) 3.8 else 8),
      axis.title = ggplot2::element_text(size = if (compact) 4.2 else 8.5),
      axis.text = ggplot2::element_text(size = if (compact) 3.6 else 7.5),
      axis.ticks = ggplot2::element_line(linewidth = 0.2),
      panel.border = if (compact) {
        ggplot2::element_rect(color = "grey60", fill = NA, linetype = "dashed", linewidth = 0.3)
      } else {
        ggplot2::element_blank()
      },
      plot.margin = ggplot2::margin(if (compact) 1.5 else 5, if (compact) 1.5 else 5,
                                    if (compact) 1.5 else 5, if (compact) 1.5 else 5)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(
      override.aes = list(alpha = 1, color = NA),
      ncol = 1
    ))
  p
}

paper_density_order <- c("CAD", "HF", "Afib", "AS")
paper_density_plots <- stats::setNames(
  lapply(paper_density_order, make_primary_density_plot, compact = FALSE),
  paper_density_order
)
paper_density_insets <- stats::setNames(
  lapply(paper_density_order, make_primary_density_plot, compact = TRUE),
  paper_density_order
)
save_plot(
  paper_density_plots[["CAD"]] | paper_density_plots[["HF"]] |
    paper_density_plots[["Afib"]] | paper_density_plots[["AS"]],
  "original_style_primary_density_insets.png",
  12.0,
  3.2
)

bed_path <- ppp_bed_file
if (is.na(bed_path) && require_paper_figures) {
  stop(
    "PAPER_FIGURE_GATE_FAILED: ppp.b38.bed was not found. Set UKBPPP_PPP_BED_FILE.",
    call. = FALSE
  )
}
if (file.exists(bed_path)) {
  bed <- data.table::fread(bed_path, header = FALSE)
  if (ncol(bed) >= 4) {
    bed <- bed[, 1:4, with = FALSE]
    data.table::setnames(bed, c("chr", "gene_start", "gene_end", "Protein_pos"))
    bed[, chr := gsub("^chr", "", as.character(chr), ignore.case = TRUE)]
    bed[, gene_start := suppressWarnings(as.numeric(gene_start))]
    bed[, gene_end := suppressWarnings(as.numeric(gene_end))]
    bed <- bed[chr %in% paper_chr_levels & is.finite(gene_start)]
    bed[, chr := factor(chr, levels = paper_chr_levels)]
    bed <- bed[order(chr, gene_start)]
    bed_duplicate_positions <- bed[duplicated(Protein_pos) | duplicated(Protein_pos, fromLast = TRUE)]
    if (nrow(bed_duplicate_positions) > 0) {
      data.table::fwrite(
        bed_duplicate_positions,
        file.path(fig_dir, "original_style_primary_manhattan_duplicate_positions.csv")
      )
    }
    bed <- bed[!duplicated(Protein_pos)]
    chr_offsets <- bed[, .(chr_len = max(gene_start, gene_end, na.rm = TRUE)), by = chr][order(chr)]
    chr_offsets[, bp_add := data.table::shift(cumsum(chr_len), fill = 0)]
    bed <- merge(bed, chr_offsets[, .(chr, bp_add)], by = "chr", all.x = TRUE)
    bed[, bp_cum := gene_start + bp_add]
    axis_set <- bed[, .(center = (min(bp_cum, na.rm = TRUE) + max(bp_cum, na.rm = TRUE)) / 2), by = chr][order(chr)]

    manhattan <- copy(primary)
    manhattan[, Protein_pos := data.table::fifelse(Protein == "NTPROBNP", "NPPB", Protein)]
    manhattan <- merge(
      manhattan,
      bed[, .(Protein_pos, chr, gene_start, bp_cum)],
      by = "Protein_pos",
      all.x = FALSE,
      allow.cartesian = TRUE
    )
    manhattan <- manhattan[is.finite(bp_cum)]
    manhattan[, p_min := sign(beta) * safe_neglog10(P_Value)]
    manhattan[, paper_sig := P_Value < paper_primary_sig]
    sig_outcome_counts <- manhattan[paper_sig %in% TRUE, .(n_sig_outcomes = data.table::uniqueN(Outcome)), by = Protein]
    multi_outcome_sig_proteins <- sig_outcome_counts[n_sig_outcomes > 1, Protein]
    paper_manhattan_y_limits <- suppressWarnings(as.numeric(strsplit(
      Sys.getenv("UKBPPP_MIAMI_Y_LIMITS", unset = "-15,65"),
      "[,;[:space:]]+"
    )[[1]]))
    if (length(paper_manhattan_y_limits) != 2 ||
        any(!is.finite(paper_manhattan_y_limits)) ||
        paper_manhattan_y_limits[1] >= paper_manhattan_y_limits[2]) {
      paper_manhattan_y_limits <- c(-15, 65)
    }
    manhattan_clipped <- manhattan[
      p_min < paper_manhattan_y_limits[1] | p_min > paper_manhattan_y_limits[2],
      .(Protein, Outcome, HR, P_Value, p_min, chr, gene_start)
    ]
    if (nrow(manhattan_clipped) > 0) {
      data.table::fwrite(
        manhattan_clipped[order(Outcome, -abs(p_min))],
        file.path(fig_dir, "original_style_primary_manhattan_clipped_extreme_points.csv")
      )
    }
    manhattan_y_limits <- c(
      min(paper_manhattan_y_limits),
      max(paper_manhattan_y_limits)
    )

    make_paper_manhattan <- function(outc, show_y = TRUE) {
      d <- copy(manhattan[as.character(Outcome) == outc])
      d[, current_outcome_sig := paper_sig %in% TRUE]
      d[, multi_outcome_sig := current_outcome_sig & (Protein %in% multi_outcome_sig_proteins)]
      d[, single_outcome_sig := current_outcome_sig & !multi_outcome_sig]
      d[, col := data.table::fifelse(
        multi_outcome_sig,
        "Associated with more than one outcome",
        data.table::fifelse(single_outcome_sig, "Associated with one outcome", as.character(chr))
      )]
      d[, col := factor(col, levels = names(paper_chr_cols))]
      label_data <- d[current_outcome_sig == TRUE][order(P_Value)][seq_len(min(.N, 25))]
      label_data <- label_data[p_min >= manhattan_y_limits[1] & p_min <= manhattan_y_limits[2]]
      if (nrow(label_data) == 0) label_data <- d[order(P_Value)][seq_len(min(.N, 5))]

      p <- ggplot2::ggplot(d, ggplot2::aes(x = bp_cum, y = p_min, color = col)) +
        ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +
        ggplot2::geom_hline(yintercept = paper_primary_line, color = "#DC0000B2", linetype = "solid", linewidth = 0.35) +
        ggplot2::geom_hline(yintercept = -paper_primary_line, color = "#DC0000B2", linetype = "solid", linewidth = 0.35) +
        ggplot2::geom_point(size = 0.7, alpha = 0.9) +
        ggrepel::geom_text_repel(
          data = label_data,
          ggplot2::aes(label = Protein),
          color = "black",
          size = 2.5,
          max.overlaps = 300,
          min.segment.length = 0,
          show.legend = FALSE
        ) +
        ggplot2::scale_color_manual(values = paper_chr_cols, drop = FALSE) +
        ggplot2::scale_x_continuous(
          labels = as.character(axis_set$chr),
          breaks = axis_set$center,
          expand = ggplot2::expansion(mult = c(0.005, 0.005))
        ) +
        ggplot2::scale_y_continuous(limits = manhattan_y_limits, expand = ggplot2::expansion(mult = c(0.01, 0.04))) +
        ggplot2::labs(
          title = outcome_full_labels[[outc]],
          x = "",
          y = if (show_y) expression(sgn(beta) %*% -log[10](P)) else ""
        ) +
        ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(
          legend.position = "none",
          axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 7),
          axis.title.y = ggplot2::element_text(size = 9),
          plot.title = ggplot2::element_text(hjust = 0.5, size = 11)
        )
      inset_spec <- paper_density_spec[outcome == outc]
      if (nrow(inset_spec) == 1L && outc %in% names(paper_density_insets)) {
        p <- p + patchwork::inset_element(
          paper_density_insets[[outc]],
          left = inset_spec$inset_left,
          bottom = inset_spec$inset_bottom,
          right = inset_spec$inset_right,
          top = inset_spec$inset_top,
          align_to = "panel",
          on_top = TRUE,
          clip = TRUE
        )
      }
      p
    }

    paper_manhattan_order <- c("CAD", "HF", "Afib", "AS")
    paper_manhattan_source <- data.table::rbindlist(lapply(paper_manhattan_order, function(outc) {
      d <- copy(manhattan[as.character(Outcome) == outc])
      d[, current_outcome_sig := paper_sig %in% TRUE]
      d[, multi_outcome_sig := current_outcome_sig & (Protein %in% multi_outcome_sig_proteins)]
      d[, single_outcome_sig := current_outcome_sig & !multi_outcome_sig]
      d[, plot_group := data.table::fifelse(
        multi_outcome_sig,
        "Associated with more than one outcome",
        data.table::fifelse(single_outcome_sig, "Associated with one outcome", "Not significant")
      )]
      d[, .(
        Protein, Outcome, chr, gene_start, bp_cum, beta, HR, P_Value, FDR,
        Bonferroni, current_outcome_sig, multi_outcome_sig, single_outcome_sig, plot_group, p_min
      )]
    }), fill = TRUE)
    data.table::fwrite(
      paper_manhattan_source,
      file.path(fig_dir, "original_style_primary_manhattan_source_data.csv")
    )
    data.table::fwrite(
      paper_manhattan_source[, .(
        n_points = .N,
        n_proteins = data.table::uniqueN(Protein),
        n_current_outcome_sig = sum(current_outcome_sig, na.rm = TRUE),
        n_multi_outcome_sig = sum(multi_outcome_sig, na.rm = TRUE),
        n_single_outcome_sig = sum(single_outcome_sig, na.rm = TRUE),
        n_clipped_outside_main_y_axis = sum(p_min < manhattan_y_limits[1] | p_min > manhattan_y_limits[2], na.rm = TRUE),
        main_y_axis_lower = manhattan_y_limits[1],
        main_y_axis_upper = manhattan_y_limits[2],
        n_strict_bonferroni = sum(P_Value < primary_bonferroni_threshold, na.rm = TRUE),
        strict_bonferroni_threshold = primary_bonferroni_threshold,
        miami_plot_threshold = paper_primary_sig,
        n_fdr = sum(FDR < 0.05, na.rm = TRUE),
        n_p_lt_0_05 = sum(P_Value < 0.05, na.rm = TRUE)
      ), by = Outcome][order(Outcome)],
      file.path(fig_dir, "original_style_primary_manhattan_counts.csv")
    )
    paper_manhattan_plots <- lapply(seq_along(paper_manhattan_order), function(i) {
      make_paper_manhattan(paper_manhattan_order[[i]], show_y = i %in% c(1, 3))
    })
    names(paper_manhattan_plots) <- paper_manhattan_order
    for (outc in paper_manhattan_order) {
      save_plot(
        paper_manhattan_plots[[outc]],
        paste0("original_style_primary_manhattan_", tolower(outcome_labels[[outc]]), ".png"),
        7.0,
        4.8
      )
    }
    save_plot(
      (paper_manhattan_plots[["CAD"]] | paper_manhattan_plots[["HF"]]) /
        (paper_manhattan_plots[["Afib"]] | paper_manhattan_plots[["AS"]]),
      "original_style_primary_manhattan_all.png",
      13.5,
      8.8
    )

    missing_positions <- unique(primary[!(Protein %in% unique(manhattan$Protein)), Protein])
    if (length(missing_positions) > 0) {
      data.table::fwrite(
        data.table::data.table(Protein = missing_positions),
        file.path(fig_dir, "original_style_primary_manhattan_missing_positions.csv")
      )
    }
  } else {
    warning("Could not make paper-style Manhattan plots: ppp.b38.bed has fewer than 4 columns.")
  }
} else {
  warning("Could not make paper-style Manhattan plots: missing ", bed_path)
}

paper_sex_sig <- 0.05 / (n_proteins * length(outcome_levels))
sex_paper <- copy(sex)
sex_paper[, diff_M_F := log(HR_M) - log(HR_F)]
sex_paper <- sex_paper[is.finite(diff_M_F)]
# The paper's lollipop colors use suggestive interaction evidence (P_int <
# 0.05) among protein-outcome pairs significant in at least one sex. The
# stricter Bonferroni interaction flag remains the formal inference result,
# but must not be used to suppress the descriptive Figure 4 colors.
sex_paper[, Figure_Suggestive_Interaction :=
  Paper_Candidate %in% TRUE & is.finite(P_Val_Int) & P_Val_Int < 0.05]
sex_paper[, sex_specific := (
  ((P_Val_F > 0.05 | log(HR_F) * log(HR_M) < 0) & P_Val_M < paper_sex_sig) |
    ((P_Val_M > 0.05 | log(HR_F) * log(HR_M) < 0) & P_Val_F < paper_sex_sig)
)]
sex_paper[, col := data.table::fifelse(
  (P_Val_F > 0.05 | log(HR_F) * log(HR_M) < 0) & P_Val_M < paper_sex_sig & Figure_Suggestive_Interaction,
  "Male-specific biomarker",
  data.table::fifelse(
    (P_Val_M > 0.05 | log(HR_F) * log(HR_M) < 0) & P_Val_F < paper_sex_sig & Figure_Suggestive_Interaction,
    "Female-specific biomarker",
    "Other"
  )
)]
sex_paper[, col_3 := data.table::fifelse(
  sex_specific & Figure_Suggestive_Interaction & diff_M_F > 0,
  "Sex-specific biomarker (higher in males)",
  data.table::fifelse(
    sex_specific & Figure_Suggestive_Interaction & diff_M_F < 0,
    "Sex-specific biomarker (higher in females)",
    data.table::fifelse(
      diff_M_F > 0 & Figure_Suggestive_Interaction,
      "Higher in males",
      data.table::fifelse(
        diff_M_F < 0 & Figure_Suggestive_Interaction,
        "Higher in females",
        "Other"
      )
    )
  )
)]
sex_paper <- sex_paper[order(-diff_M_F)]
sex_paper[, num_peroutc := seq_len(.N), by = Outcome]
sex_paper[, col_3 := factor(
  col_3,
  levels = c(
    "Sex-specific biomarker (higher in females)",
    "Higher in females",
    "Other",
    "Higher in males",
    "Sex-specific biomarker (higher in males)"
  )
)]

# Keep the descriptive lollipop readable: label no more than seven proteins
# per outcome. Formal interaction findings take priority, followed by the
# smallest interaction P value and then the largest absolute sex difference.
sex_label_max <- 7L
sex_paper[, `:=`(Label_Selected = FALSE, Label_Rank = NA_integer_)]
sex_paper[col != "Other", Label_Rank := {
  label_order <- order(
    -as.integer(Paper_Interaction_Bonferroni %in% TRUE),
    P_Val_Int,
    -abs(diff_M_F),
    Protein,
    na.last = TRUE
  )
  label_rank <- integer(.N)
  label_rank[label_order] <- seq_along(label_order)
  label_rank
}, by = Outcome]
sex_paper[!is.na(Label_Rank) & Label_Rank <= sex_label_max, Label_Selected := TRUE]
sex_diff_cols <- c(
  "Sex-specific biomarker (higher in females)" = "#fa6f8f",
  "Higher in females" = "#FDB7C7",
  "Other" = "grey90",
  "Higher in males" = "#BAB1D1",
  "Sex-specific biomarker (higher in males)" = "#7563a3"
)
sex_diff_labels <- c(
  "Sex-specific biomarker (higher in females)" = "Sex-specific; log(HR) higher in female participants",
  "Higher in females" = "log(HR) higher in female participants",
  "Other" = "Other",
  "Higher in males" = "log(HR) higher in male participants",
  "Sex-specific biomarker (higher in males)" = "Sex-specific; log(HR) higher in male participants"
)

data.table::fwrite(
  sex_paper[, .(
    Protein, Outcome, HR_M, HR_F, P_Val_M, P_Val_F, P_Val_Int,
    Paper_Candidate, Paper_Interaction_P_Adjusted,
    Paper_Interaction_Bonferroni, Figure_Suggestive_Interaction,
    diff_M_F, num_peroutc, col_3, Label_Selected, Label_Rank
  )],
  file.path(fig_dir, "original_style_sex_difference_source_data.csv")
)

make_paper_sex_diff <- function(outc, show_y = TRUE) {
  d <- copy(sex_paper[as.character(Outcome) == outc])
  y_limit <- max(0.45, ceiling(max(abs(d$diff_M_F), na.rm = TRUE) * 10) / 10)
  y_limit <- min(y_limit, 1.0)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = num_peroutc, y = diff_M_F)) +
    ggplot2::geom_segment(ggplot2::aes(xend = num_peroutc, y = 0, yend = diff_M_F), color = "grey90", linewidth = 0.25) +
    ggplot2::geom_point(data = d[col_3 == "Other"], color = "grey90", size = 1.6)
  if (nrow(d[col_3 != "Other"]) > 0) {
    p <- p +
      ggplot2::geom_point(data = d[col_3 != "Other"], ggplot2::aes(color = col_3), size = 1.8) +
      ggplot2::scale_color_manual(values = sex_diff_cols, labels = sex_diff_labels, drop = FALSE, name = "")
  }
  if (nrow(d[Label_Selected %in% TRUE]) > 0) {
    p <- p +
      ggplot2::geom_point(data = d[Label_Selected %in% TRUE], ggplot2::aes(color = col_3), size = 2.2) +
      ggrepel::geom_label_repel(
        data = d[Label_Selected %in% TRUE],
        ggplot2::aes(label = Protein),
        color = "black",
        fill = "white",
        label.size = 0.15,
        size = 2.4,
        max.overlaps = Inf,
        force = 5,
        show.legend = FALSE
      )
  }
  p +
    ggplot2::scale_y_continuous(limits = c(-y_limit, y_limit), breaks = pretty(c(-y_limit, y_limit), n = 5)) +
    ggplot2::labs(
      title = outcome_full_labels[[outc]],
      x = "",
      y = if (show_y) "log(HR)males - log(HR)females" else "",
      color = ""
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 11)
    )
}

paper_sex_order <- c("CAD", "HF", "Afib", "AS")
paper_sex_diff_plots <- lapply(seq_along(paper_sex_order), function(i) {
  make_paper_sex_diff(paper_sex_order[[i]], show_y = i %in% c(1, 3))
})
names(paper_sex_diff_plots) <- paper_sex_order
save_plot(
  ((paper_sex_diff_plots[["CAD"]] | paper_sex_diff_plots[["HF"]]) /
    (paper_sex_diff_plots[["Afib"]] | paper_sex_diff_plots[["AS"]])) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom"),
  "original_style_sex_difference_alloutcomes.png",
  10.5,
  7.4
)

sex_top5 <- sex[Paper_Candidate %in% TRUE][order(P_Val_Int), head(.SD, 5), by = Outcome]
sex_top5[, rank := seq_len(.N), by = Outcome]
sex_hr_long <- data.table::rbindlist(list(
  sex_top5[, .(Protein, sex = "Female", hr = HR_F, lci = CI_Low_F, uci = CI_High_F, rank, Outcome)],
  sex_top5[, .(Protein, sex = "Male", hr = HR_M, lci = CI_Low_M, uci = CI_High_M, rank, Outcome)]
))
sex_hr_long[, sex := factor(sex, levels = c("Female", "Male"))]
data.table::fwrite(sex_hr_long, file.path(fig_dir, "original_style_sex_hr_top5_source_data.csv"))

make_paper_sex_hr <- function(outc, show_y = TRUE) {
  d <- copy(sex_hr_long[as.character(Outcome) == outc])
  source_order <- sex_top5[as.character(Outcome) == outc][order(rank), Protein]
  rank_levels <- rev(as.character(seq_along(source_order)))
  label_values <- stats::setNames(rev(source_order), rank_levels)
  ggplot2::ggplot(d, ggplot2::aes(x = hr, y = as.factor(rank), color = sex, xmin = lci, xmax = uci)) +
    ggplot2::geom_pointrange(size = 0.3, shape = 20, position = ggplot2::position_dodge2(width = 0.65, padding = 0.5)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3) +
    ggplot2::scale_color_manual(values = c("Female" = "#fa6f8f", "Male" = "#7563a3"), name = "") +
    ggplot2::scale_y_discrete(limits = rank_levels, labels = label_values) +
    ggplot2::scale_x_continuous(trans = "log", breaks = c(0.8, 1.2, 1.6, 2.0)) +
    ggplot2::coord_cartesian(xlim = c(0.7, 2.1)) +
    ggplot2::labs(
      title = outcome_full_labels[[outc]],
      x = "HR (95% CI)",
      y = if (show_y) "" else NULL
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 11)
    )
}

paper_sex_hr_plots <- lapply(seq_along(paper_sex_order), function(i) {
  make_paper_sex_hr(paper_sex_order[[i]], show_y = i %in% c(1, 3))
})
names(paper_sex_hr_plots) <- paper_sex_order
save_plot(
  (paper_sex_hr_plots[["CAD"]] | paper_sex_hr_plots[["HF"]]) /
    (paper_sex_hr_plots[["Afib"]] | paper_sex_hr_plots[["AS"]]),
  "original_style_sex_hr_top5.png",
  8.0,
  6.0
)

make_paper_sex_panel_with_inset <- function(outc, show_y = TRUE) {
  base_panel <- make_paper_sex_diff(outc, show_y = show_y) +
    ggplot2::theme(legend.position = "none")
  inset_panel <- make_paper_sex_hr(outc, show_y = FALSE) +
    ggplot2::labs(title = NULL, y = NULL) +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(1, 1, 1, 1),
      axis.text.y = ggplot2::element_text(size = 6),
      axis.text.x = ggplot2::element_text(size = 6),
      axis.title.x = ggplot2::element_text(size = 7)
    )
  base_panel + patchwork::inset_element(inset_panel, left = 0.50, bottom = 0.55, right = 0.98, top = 0.98)
}

paper_sex_inset_plots <- lapply(seq_along(paper_sex_order), function(i) {
  make_paper_sex_panel_with_inset(paper_sex_order[[i]], show_y = i %in% c(1, 3))
})
names(paper_sex_inset_plots) <- paper_sex_order
save_plot(
  (paper_sex_inset_plots[["CAD"]] | paper_sex_inset_plots[["HF"]]) /
    (paper_sex_inset_plots[["Afib"]] | paper_sex_inset_plots[["AS"]]),
  "original_style_sex_difference_with_hr_insets.png",
  11.0,
  7.6
)

if (length(roc_files) > 0) {
  roc_paper <- data.table::rbindlist(lapply(roc_files, data.table::fread), fill = TRUE)
  roc_paper[, outc_full := factor(lasso_outcome_full_labels[outc], levels = lasso_outcome_full_labels[c("cad", "hfail", "afib", "ao_sten")])]
  roc_paper[, model_paper := factor(
    c(Clinical = "Clinical parameters", Proteins = "Proteins", Combination = "Combined")[as.character(model)],
    levels = c("Clinical parameters", "Proteins", "Combined")
  )]
  paper_roc_cols <- c("Clinical parameters" = "black", "Combined" = "#1EC000", "Proteins" = "#C2EEB9")
  make_paper_roc <- function(outc_code, show_y = FALSE) {
    d <- roc_paper[roc_paper$outc == outc_code]
    auc_labels <- copy(auc[outc == outc_code])
    auc_labels[, model_paper := factor(
      c(Clinical = "Clinical parameters", Proteins = "Proteins", Combination = "Combined")[as.character(model)],
      levels = c("Clinical parameters", "Proteins", "Combined")
    )]
    auc_labels <- auc_labels[match(c("Clinical parameters", "Proteins", "Combined"), as.character(model_paper))]
    auc_text <- paste(
      "AUC (95% CI)",
      paste0(
        as.character(auc_labels$model_paper), ": ",
        sprintf("%.3f (%.3f-%.3f)", auc_labels$auc, auc_labels$auc_lci, auc_labels$auc_uci)
      ),
      collapse = "\n"
    )
    ggplot2::ggplot(d, ggplot2::aes(x = FPR, y = TPR, color = model_paper)) +
      ggplot2::geom_line(linewidth = 0.55) +
      ggplot2::geom_abline(linetype = "dotted", color = "grey90", linewidth = 0.35) +
      ggplot2::annotate(
        "text",
        x = 0.98,
        y = 0.10,
        label = auc_text,
        hjust = 1,
        vjust = 0,
        size = 2.35,
        color = "grey38",
        lineheight = 0.95
      ) +
      ggplot2::scale_color_manual(values = paper_roc_cols, name = "", drop = FALSE) +
      ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      ggplot2::labs(
        title = lasso_outcome_full_labels[[outc_code]],
        x = "False positive rate",
        y = if (show_y) "True positive rate" else ""
      ) +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 11)
      )
  }
  paper_roc_order <- c("cad", "hfail", "afib", "ao_sten")
  paper_roc_plots <- lapply(seq_along(paper_roc_order), function(i) {
    make_paper_roc(paper_roc_order[[i]], show_y = i == 1)
  })
  names(paper_roc_plots) <- paper_roc_order
  paper_roc_combined <- (paper_roc_plots[["cad"]] | paper_roc_plots[["hfail"]] |
    paper_roc_plots[["afib"]] | paper_roc_plots[["ao_sten"]]) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  save_plot(
    paper_roc_combined,
    "original_style_lasso_rocs.png",
    14.5,
    4.4
  )

  prediction_path <- file.path(output_dir, "lasso_predictions.csv")
  if (file.exists(prediction_path)) {
    pred_paper <- data.table::fread(prediction_path)
    pred_paper[, outc := as.character(outc)]
    pred_paper[, model := as.character(model)]
    pred_paper[, original_order := .I]
    pred_paper[, pred_row := seq_len(.N), by = .(outc, model)]

    has_eid <- "eid" %in% names(pred_paper)
    has_followup <- "followup_years" %in% names(pred_paper)
    if ((!has_eid || !has_followup) && file.exists(analysis_base_file)) {
      base_for_pred <- readRDS(analysis_base_file)
      pred_base <- data.table::as.data.table(base_for_pred$dat)
      lookup_rows <- lapply(seq_len(nrow(base_for_pred$outcome_map)), function(i) {
        key <- base_for_pred$outcome_map$outcome_key[i]
        event_col <- paste0(key, "_inc")
        fu_col <- paste0(key, "_fu")
        d0 <- pred_base[
          !is.na(get(event_col)) & !is.na(get(fu_col)) &
            is.finite(get(fu_col)) & get(fu_col) > 0
        ]
        y0 <- as.integer(d0[[event_col]] == 1)
        set.seed(1234)
        test_idx <- unlist(tapply(seq_along(y0), y0, function(idx) {
          sample(idx, max(1, floor(length(idx) * 0.2)))
        }))
        data.table::data.table(
          outc = key,
          pred_row = seq_along(test_idx),
          eid_lookup = if ("eid" %in% names(d0)) d0$eid[test_idx] else if ("id" %in% names(d0)) d0$id[test_idx] else test_idx,
          followup_years_lookup = as.numeric(d0[[fu_col]][test_idx]),
          y_lookup = y0[test_idx]
        )
      })
      pred_lookup <- data.table::rbindlist(lookup_rows, fill = TRUE)
      pred_paper <- merge(pred_paper, pred_lookup, by = c("outc", "pred_row"), all.x = TRUE, sort = FALSE)
      if (!has_eid) {
        pred_paper[, eid := eid_lookup]
      } else {
        pred_paper[is.na(eid), eid := eid_lookup]
      }
      if (!has_followup) {
        pred_paper[, followup_years := followup_years_lookup]
      } else {
        pred_paper[!is.finite(followup_years), followup_years := followup_years_lookup]
      }
      if (!"y" %in% names(pred_paper)) pred_paper[, y := y_lookup]
      pred_paper[, c("eid_lookup", "followup_years_lookup", "y_lookup") := NULL]
      pred_paper <- pred_paper[order(original_order)]
    }

    if (!"score_sd" %in% names(pred_paper)) {
      pred_paper[, score_sd := as.numeric(scale(stats::qlogis(pmin(pmax(predicted_risk, 1e-6), 1 - 1e-6)))),
                 by = .(outc, model)]
    }
    pred_paper[, outc_full := factor(lasso_outcome_full_labels[outc], levels = lasso_outcome_full_labels[c("cad", "hfail", "afib", "ao_sten")])]
    pred_paper[, model_paper := factor(
      c(Clinical = "Clinical parameters", Proteins = "Proteins", Combination = "Combined")[as.character(model)],
      levels = c("Clinical parameters", "Proteins", "Combined")
    )]

    protein_score <- pred_paper[model == "Proteins" & is.finite(score_sd)]
    if (nrow(protein_score) > 0) {
      protein_score[, event_status := factor(ifelse(y == 1, "Cases", "Controls"), levels = c("Controls", "Cases"))]
      protein_score[, score_percentile := 100 * (data.table::frank(score_sd, ties.method = "average") - 0.5) / .N, by = outc]
      protein_score[, score_group := cut(
        score_percentile,
        breaks = c(0, 20, 40, 60, 80, 100),
        include.lowest = TRUE,
        labels = c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%")
      )]
      protein_score[, score_decile := cut(
        score_percentile,
        breaks = seq(0, 100, by = 10),
        include.lowest = TRUE,
        labels = paste0(seq(0, 90, by = 10), "-", seq(10, 100, by = 10), "%")
      )]

      score_group_cols <- c(
        "0-20%" = "#FDE0C5",
        "20-40%" = "#FDBF7E",
        "40-60%" = "#F28E2B",
        "60-80%" = "#D07A13",
        "80-100%" = "#9C5A05"
      )
      cuminc_ymax <- protein_score[!is.na(score_group), .(event_rate = mean(y == 1, na.rm = TRUE)), by = .(outc, score_group)]
      cuminc_ymax <- min(0.35, max(0.05, max(cuminc_ymax$event_rate, na.rm = TRUE) * 1.35))
      density_fills <- c(Controls = "#F6E8D5", Cases = "#C9A37A")
      blank_prediction_panel <- function(title, message, tag = NULL) {
        p <- ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5, label = message, size = 3.2, color = "grey40") +
          ggplot2::xlim(0, 1) +
          ggplot2::ylim(0, 1) +
          ggplot2::labs(title = title) +
          ggplot2::theme_void(base_size = 10) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 11))
        if (!is.null(tag)) {
          p <- p + ggplot2::labs(tag = tag) +
            ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 13))
        }
        p
      }

      make_density_panel <- function(outc_code, show_y = FALSE, tag = NULL) {
        d <- protein_score[outc == outc_code]
        if (nrow(d) < 20 || length(unique(d$event_status)) < 2) {
          return(blank_prediction_panel(lasso_outcome_full_labels[[outc_code]], "Insufficient test-set data", tag))
        }
        threshold <- stats::quantile(d[event_status == "Controls", score_sd], probs = 0.95, na.rm = TRUE, type = 8)
        dr <- mean(d[event_status == "Cases", score_sd] >= threshold, na.rm = TRUE)
        fpr <- mean(d[event_status == "Controls", score_sd] >= threshold, na.rm = TRUE)
        p <- ggplot2::ggplot(d, ggplot2::aes(score_sd, fill = event_status, color = event_status)) +
          ggplot2::geom_density(alpha = 0.65, linewidth = 0.45, adjust = 1.05) +
          ggplot2::geom_vline(xintercept = threshold, color = "grey35", linewidth = 0.35) +
          ggplot2::annotate(
            "text",
            x = Inf,
            y = Inf,
            label = sprintf("DR = %.1f%%\nFPR = %.1f%%", 100 * dr, 100 * fpr),
            hjust = 1.04,
            vjust = 1.12,
            size = 2.6,
            color = "grey35"
          ) +
          ggplot2::scale_fill_manual(values = density_fills, name = "") +
          ggplot2::scale_color_manual(values = c(Controls = "grey25", Cases = "grey25"), guide = "none") +
          ggplot2::labs(
            title = lasso_outcome_full_labels[[outc_code]],
            x = "Protein score (s.d.)",
            y = if (show_y) "Density" else ""
          ) +
          ggplot2::theme_classic(base_size = 10) +
          ggplot2::theme(
            legend.position = if (outc_code == "cad") c(0.78, 0.80) else "none",
            legend.background = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(hjust = 0.5, size = 11)
          )
        if (!is.null(tag)) {
          p <- p + ggplot2::labs(tag = tag) +
            ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 13))
        }
        p
      }

      make_cuminc_panel <- function(outc_code, show_y = FALSE, tag = NULL) {
        d <- protein_score[outc == outc_code & is.finite(followup_years) & followup_years > 0 & !is.na(score_group)]
        if (nrow(d) < 20 || sum(d$y == 1, na.rm = TRUE) == 0) {
          return(blank_prediction_panel(lasso_outcome_full_labels[[outc_code]], "Follow-up data unavailable", tag))
        }
        fit <- tryCatch(
          survival::survfit(survival::Surv(followup_years, y == 1) ~ score_group, data = d),
          error = function(e) NULL
        )
        if (is.null(fit)) return(blank_prediction_panel(lasso_outcome_full_labels[[outc_code]], "Survival fit failed", tag))
        s <- summary(fit)
        cum_dt <- data.table::data.table(
          time = s$time,
          cumulative_incidence = 1 - s$surv,
          score_group = sub("^score_group=", "", s$strata)
        )
        cum_dt[, score_group := factor(score_group, levels = names(score_group_cols))]
        p <- ggplot2::ggplot(cum_dt, ggplot2::aes(time, cumulative_incidence, color = score_group)) +
          ggplot2::geom_step(linewidth = 0.6) +
          ggplot2::scale_color_manual(values = score_group_cols, breaks = rev(names(score_group_cols)), name = "Protein score percentile:") +
          ggplot2::coord_cartesian(xlim = c(0, min(15, max(cum_dt$time, na.rm = TRUE))), ylim = c(0, cuminc_ymax)) +
          ggplot2::labs(
            title = lasso_outcome_full_labels[[outc_code]],
            x = "Follow-up time (years)",
            y = if (show_y) "Cumulative incidence" else ""
          ) +
          ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
          ggplot2::theme_classic(base_size = 10) +
          ggplot2::theme(
            legend.position = if (outc_code == "cad") c(0.56, 0.88) else "none",
            legend.background = ggplot2::element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0.75), color = NA),
            legend.title = ggplot2::element_text(size = 6.5),
            legend.text = ggplot2::element_text(size = 6.5),
            legend.key.width = grid::unit(0.9, "lines"),
            legend.key.height = grid::unit(0.55, "lines"),
            plot.title = ggplot2::element_text(hjust = 0.5, size = 11)
          )
        if (!is.null(tag)) {
          p <- p + ggplot2::labs(tag = tag) +
            ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 13))
        }
        p
      }

      make_rate_panel <- function(outc_code, show_y = FALSE, tag = NULL) {
        d <- protein_score[outc == outc_code & is.finite(followup_years) & followup_years > 0 & !is.na(score_decile)]
        if (nrow(d) < 20) return(blank_prediction_panel(lasso_outcome_full_labels[[outc_code]], "Follow-up data unavailable", tag))
        rate_dt <- d[, .(
          x_mid = mean(score_percentile, na.rm = TRUE),
          events = sum(y == 1, na.rm = TRUE),
          person_years = sum(followup_years, na.rm = TRUE)
        ), by = score_decile]
        rate_dt[, incidence_rate := 1000 * events / person_years]
        rate_dt[incidence_rate <= 0 | !is.finite(incidence_rate), incidence_rate := NA_real_]
        p <- ggplot2::ggplot(rate_dt, ggplot2::aes(x_mid, incidence_rate, color = x_mid)) +
          ggplot2::geom_point(size = 1.7, na.rm = TRUE) +
          ggplot2::scale_color_gradient(low = "#FDE0C5", high = "#9C5A05", guide = "none") +
          ggplot2::scale_y_log10(breaks = c(0.25, 2, 16), labels = c("0.25", "2.00", "16.00")) +
          ggplot2::coord_cartesian(xlim = c(0, 100)) +
          ggplot2::labs(
            title = lasso_outcome_full_labels[[outc_code]],
            x = "Protein score percentile",
            y = if (show_y) "Incidence rate per 1,000 years" else ""
          ) +
          ggplot2::theme_classic(base_size = 10) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 11))
        if (!is.null(tag)) {
          p <- p + ggplot2::labs(tag = tag) +
            ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 13))
        }
        p
      }

      make_prediction_roc_panel <- function(outc_code, show_y = FALSE, tag = NULL) {
        p <- make_paper_roc(outc_code, show_y = show_y) +
          ggplot2::labs(x = "FPR", y = if (show_y) "True positive rate" else "") +
          ggplot2::theme(legend.position = "none")
        if (!is.null(tag)) {
          p <- p + ggplot2::labs(tag = tag) +
            ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 13))
        }
        p
      }

      prediction_order <- c("cad", "hfail", "afib", "ao_sten")
      density_row <- patchwork::wrap_plots(lapply(seq_along(prediction_order), function(i) {
        make_density_panel(prediction_order[[i]], show_y = i == 1, tag = if (i == 1) "a" else NULL)
      }), nrow = 1)
      cuminc_row <- patchwork::wrap_plots(lapply(seq_along(prediction_order), function(i) {
        make_cuminc_panel(prediction_order[[i]], show_y = i == 1, tag = if (i == 1) "b" else NULL)
      }), nrow = 1)
      rate_row <- patchwork::wrap_plots(lapply(seq_along(prediction_order), function(i) {
        make_rate_panel(prediction_order[[i]], show_y = i == 1, tag = if (i == 1) "c" else NULL)
      }), nrow = 1)
      roc_row <- patchwork::wrap_plots(lapply(seq_along(prediction_order), function(i) {
        make_prediction_roc_panel(prediction_order[[i]], show_y = i == 1, tag = if (i == 1) "d" else NULL)
      }), nrow = 1)
      prediction_figure <- density_row / cuminc_row / rate_row / roc_row +
        patchwork::plot_layout(heights = c(1.0, 1.1, 1.1, 1.25))
      save_plot(
        prediction_figure,
        "original_style_lasso_prediction_figure.png",
        14.5,
        14.2
      )
    }
  }
}

coef_files <- list.files(file.path(output_dir, "lasso"), pattern = "_coefficients\\.csv$", full.names = TRUE)
if (length(coef_files) > 0) {
  coef_dt <- data.table::rbindlist(lapply(coef_files, data.table::fread), fill = TRUE)
  coef_dt <- coef_dt[model == "Proteins" & feature != "(Intercept)" & is.finite(coefficient) & coefficient != 0]
  if (nrow(coef_dt) > 0) {
    coef_dt[, outc_full := factor(lasso_outcome_full_labels[outc], levels = lasso_outcome_full_labels[c("cad", "hfail", "afib", "ao_sten")])]
    coef_dt[, coefficient_abs := abs(coefficient)]
    coef_dt[, feature_label := data.table::fifelse(feature == "NTPROBNP", "NT-proBNP", feature)]
    coef_top <- coef_dt[order(-coefficient_abs), head(.SD, 40), by = outc]
    data.table::fwrite(
      coef_top[order(outc, -coefficient_abs), .(outc, feature = feature_label, coefficient, coefficient_abs)],
      file.path(fig_dir, "original_style_lasso_protein_features_source_data.csv")
    )

    make_paper_feature_bar <- function(outc_code) {
      d <- copy(coef_top[coef_top$outc == outc_code][order(coefficient_abs)])
      d[, feature_plot := factor(feature_label, levels = feature_label)]
      ggplot2::ggplot(d, ggplot2::aes(x = feature_plot, y = coefficient_abs)) +
        ggplot2::geom_col(color = "black", fill = "black", linewidth = 0.2) +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = lasso_outcome_full_labels[[outc_code]],
          x = "",
          y = "Regression coefficient in prediction model"
        ) +
        ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(
          axis.text.y = ggplot2::element_text(size = 5),
          plot.title = ggplot2::element_text(hjust = 0.5, size = 11)
        )
    }
    paper_feature_order <- c("cad", "hfail", "afib", "ao_sten")
    paper_feature_plots <- lapply(paper_feature_order, make_paper_feature_bar)
    names(paper_feature_plots) <- paper_feature_order
    save_plot(
      (paper_feature_plots[["cad"]] | paper_feature_plots[["hfail"]]) /
        (paper_feature_plots[["afib"]] | paper_feature_plots[["ao_sten"]]),
      "original_style_lasso_protein_features.png",
      10.5,
      12.0
    )
  }
}

required_paper_figure_files <- c(
  "original_style_primary_manhattan_all.png",
  "original_style_primary_density_insets.png",
  "original_style_primary_manhattan_cad.png",
  "original_style_primary_manhattan_hf.png",
  "original_style_primary_manhattan_af.png",
  "original_style_primary_manhattan_as.png",
  "original_style_sex_difference_alloutcomes.png",
  "original_style_sex_hr_top5.png",
  "original_style_sex_difference_with_hr_insets.png",
  "original_style_lasso_rocs.png",
  "original_style_lasso_prediction_figure.png",
  "original_style_lasso_prediction_metrics.png",
  "original_style_lasso_calibration.png",
  "original_style_lasso_protein_features.png"
)
missing_paper_figures <- required_paper_figure_files[
  !file.exists(file.path(fig_dir, required_paper_figure_files))
]
figure_qc <- data.table::data.table(
  file = required_paper_figure_files,
  exists = file.exists(file.path(fig_dir, required_paper_figure_files)),
  path = file.path(fig_dir, required_paper_figure_files)
)
data.table::fwrite(figure_qc, file.path(fig_dir, "paper_figure_qc.csv"))
if (require_paper_figures && length(missing_paper_figures) > 0L) {
  stop(
    "PAPER_FIGURE_GATE_FAILED: missing ",
    paste(missing_paper_figures, collapse = ", "),
    call. = FALSE
  )
}

figure_manifest <- data.table::data.table(
  file = list.files(fig_dir, pattern = "\\.(png|csv)$", full.names = FALSE),
  path = list.files(fig_dir, pattern = "\\.(png|csv)$", full.names = TRUE)
)
data.table::fwrite(figure_manifest, file.path(fig_dir, "figure_manifest.csv"))
message("Figures saved to: ", fig_dir)
print(figure_manifest)
