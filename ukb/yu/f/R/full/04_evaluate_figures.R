if (identical(Sys.info()[["sysname"]], "Darwin") && !"Arial" %in% names(grDevices::pdfFonts())) {
  grDevices::pdfFonts(Arial = grDevices::pdfFonts("Helvetica")[[1]])
}

yur_auc <- function(y, p) {
  ok <- is.finite(y) & is.finite(p)
  y <- as.integer(y[ok]); p <- p[ok]
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  if (!n1 || !n0) return(NA_real_)
  ranks <- rank(p, ties.method = "average")
  (sum(ranks[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

yur_calibration_metrics <- function(y, p) {
  p <- pmin(1 - 1e-9, pmax(1e-9, p))
  linear_predictor <- qlogis(p)
  intercept <- tryCatch(
    unname(coef(suppressWarnings(glm(y ~ offset(linear_predictor), family = binomial())))[[1]]),
    error = function(e) NA_real_
  )
  slope <- tryCatch(
    unname(coef(suppressWarnings(glm(y ~ linear_predictor, family = binomial())))[[2]]),
    error = function(e) NA_real_
  )
  c(calibration_intercept = intercept, calibration_slope = slope)
}

yur_binary_metrics <- function(y, p, threshold) {
  p <- pmin(1 - 1e-9, pmax(1e-9, p))
  pred <- as.integer(p >= threshold)
  tp <- sum(pred == 1L & y == 1L); tn <- sum(pred == 0L & y == 0L)
  fp <- sum(pred == 1L & y == 0L); fn <- sum(pred == 0L & y == 1L)
  c(
    auc = yur_auc(y, p), accuracy = (tp + tn) / length(y),
    sensitivity = tp / pmax(tp + fn, 1), specificity = tn / pmax(tn + fp, 1),
    f1 = 2 * tp / pmax(2 * tp + fp + fn, 1), brier = mean((y - p)^2),
    yur_calibration_metrics(y, p)
  )
}

yur_nri_idi <- function(y, p_new, p_old) {
  up <- p_new > p_old; down <- p_new < p_old
  cases <- y == 1L; controls <- y == 0L
  c(
    nri = mean(up[cases]) - mean(down[cases]) + mean(down[controls]) - mean(up[controls]),
    idi = (mean(p_new[cases]) - mean(p_new[controls])) -
      (mean(p_old[cases]) - mean(p_old[controls]))
  )
}

yur_evaluate_prediction <- function(cfg) {
  prediction_file <- file.path(cfg$paths$models, "test_predictions.csv.gz")
  if (!file.exists(prediction_file)) stop("Run train first: ", prediction_file)
  prediction <- fread(prediction_file)
  required_models <- c("SCORE2", "Protein", "Protein_SCORE2")
  if (!all(required_models %in% prediction$model_id)) stop("Final prediction file lacks the three published model families.")
  set.seed(cfg$bootstrap_seed)
  metric_rows <- list(); comparison_rows <- list(); bootstrap_rows <- list()
  endpoints <- unique(prediction$outcome_id)
  for (endpoint in endpoints) {
    endpoint_data <- prediction[outcome_id == endpoint]
    wide <- dcast(endpoint_data, eid + y ~ model_id, value.var = "prediction")
    thresholds <- unique(endpoint_data[, .(model_id, threshold)])
    n <- nrow(wide)
    indices <- replicate(cfg$bootstrap_n, sample.int(n, n, replace = TRUE), simplify = FALSE)

    endpoint_models <- intersect(required_models, names(wide))
    for (model in endpoint_models) {
      threshold <- thresholds[model_id == model, threshold][[1]]
      observed <- yur_binary_metrics(wide$y, wide[[model]], threshold)
      boot <- rbindlist(lapply(seq_along(indices), function(b) {
        idx <- indices[[b]]
        values <- tryCatch(
          yur_binary_metrics(wide$y[idx], wide[[model]][idx], threshold),
          error = function(e) setNames(rep(NA_real_, 8), c("auc", "accuracy", "sensitivity", "specificity", "f1", "brier", "calibration_intercept", "calibration_slope"))
        )
        data.table(bootstrap = b, metric = names(values), value = as.numeric(values))
      }))
      quantiles <- boot[, .(
        q025 = quantile(value, .025, na.rm = TRUE), q25 = quantile(value, .25, na.rm = TRUE),
        q50 = quantile(value, .50, na.rm = TRUE), q75 = quantile(value, .75, na.rm = TRUE),
        q975 = quantile(value, .975, na.rm = TRUE)
      ), by = metric]
      row <- data.table(outcome_id = endpoint, model_id = model, n = n, events = sum(wide$y), threshold = threshold)
      for (metric_name in names(observed)) {
        row[[metric_name]] <- observed[[metric_name]]
        q <- quantiles[metric == metric_name]
        if (nrow(q) != 1L) stop("Bootstrap quantiles are not unique for metric=", metric_name)
        row[[paste0(metric_name, "_q025")]] <- q$q025
        row[[paste0(metric_name, "_q25")]] <- q$q25
        row[[paste0(metric_name, "_q50")]] <- q$q50
        row[[paste0(metric_name, "_q75")]] <- q$q75
        row[[paste0(metric_name, "_q975")]] <- q$q975
      }
      metric_rows[[length(metric_rows) + 1L]] <- row
      boot[, `:=`(outcome_id = endpoint, model_id = model)]
      bootstrap_rows[[length(bootstrap_rows) + 1L]] <- boot
    }

    pairs <- data.table(
      comparison_id = c("Protein_vs_SCORE2", "Protein_SCORE2_vs_SCORE2", "Protein_SCORE2_vs_Protein"),
      new_model = c("Protein", "Protein_SCORE2", "Protein_SCORE2"),
      old_model = c("SCORE2", "SCORE2", "Protein")
    )
    for (i in seq_len(nrow(pairs))) {
      new_model <- pairs$new_model[[i]]; old_model <- pairs$old_model[[i]]
      observed_new <- yur_auc(wide$y, wide[[new_model]])
      observed_old <- yur_auc(wide$y, wide[[old_model]])
      delong <- pROC::roc.test(
        pROC::roc(wide$y, wide[[new_model]], quiet = TRUE, direction = "<"),
        pROC::roc(wide$y, wide[[old_model]], quiet = TRUE, direction = "<"),
        paired = TRUE, method = "delong"
      )
      boot_pair <- rbindlist(lapply(seq_along(indices), function(b) {
        idx <- indices[[b]]
        delta <- tryCatch(yur_auc(wide$y[idx], wide[[new_model]][idx]) - yur_auc(wide$y[idx], wide[[old_model]][idx]), error = function(e) NA_real_)
        nri_idi <- tryCatch(yur_nri_idi(wide$y[idx], wide[[new_model]][idx], wide[[old_model]][idx]), error = function(e) c(nri = NA, idi = NA))
        data.table(bootstrap = b, delta_auc = delta, nri = nri_idi[["nri"]], idi = nri_idi[["idi"]])
      }))
      observed_ni <- yur_nri_idi(wide$y, wide[[new_model]], wide[[old_model]])
      comparison_rows[[length(comparison_rows) + 1L]] <- data.table(
        outcome_id = endpoint, comparison_id = pairs$comparison_id[[i]],
        new_model = new_model, old_model = old_model, n = n, events = sum(wide$y),
        auc_new = observed_new, auc_old = observed_old, delta_auc = observed_new - observed_old,
        delta_auc_ci_low = quantile(boot_pair$delta_auc, .025, na.rm = TRUE),
        delta_auc_ci_high = quantile(boot_pair$delta_auc, .975, na.rm = TRUE),
        delong_p = delong$p.value, nri = observed_ni[["nri"]],
        nri_ci_low = quantile(boot_pair$nri, .025, na.rm = TRUE),
        nri_ci_high = quantile(boot_pair$nri, .975, na.rm = TRUE),
        idi = observed_ni[["idi"]], idi_ci_low = quantile(boot_pair$idi, .025, na.rm = TRUE),
        idi_ci_high = quantile(boot_pair$idi, .975, na.rm = TRUE)
      )
    }
  }
  metrics <- rbindlist(metric_rows, fill = TRUE)
  comparisons <- rbindlist(comparison_rows, fill = TRUE)
  bootstrap <- rbindlist(bootstrap_rows, fill = TRUE)
  yur_write_csv(metrics, file.path(cfg$paths$evaluation, "table_s10_prediction_metrics.csv"))
  yur_write_csv(comparisons, file.path(cfg$paths$evaluation, "table_s11_model_comparisons_nri_idi.csv"))
  cad_reference_file <- file.path(cfg$project_dir, "f", "config", "yu_cad_published_metrics.csv")
  if (file.exists(cad_reference_file) && "cad" %in% metrics$outcome_id) {
    published <- fread(cad_reference_file)
    published[, auc_published_point := suppressWarnings(as.numeric(sub(" .*", "", auc_published)))]
    published[, model_local := model]
    published[model == "SCORE.V2", model_local := "SCORE2"]
    published[model == "SCORE.V2+Protein", model_local := "Protein_SCORE2"]
    local_cad <- metrics[outcome_id == "cad", .(
      model_local = model_id, auc_local = auc, auc_local_ci_low = auc_q025, auc_local_ci_high = auc_q975
    )]
    replication <- merge(published, local_cad, by = "model_local", all.x = TRUE, sort = FALSE)
    setcolorder(replication, c("model", "model_local", setdiff(names(replication), c("model", "model_local"))))
    replication[, auc_difference_local_minus_published := auc_local - auc_published_point]
    replication[, status := fifelse(
      is.na(auc_local), "FAIL_MISSING_LOCAL_MODEL",
      fifelse(abs(auc_difference_local_minus_published) <= .05, "PASS_WITHIN_0.05", "WARN_DIFFERENCE_GT_0.05")
    )]
    yur_write_csv(replication, file.path(cfg$paths$evaluation, "cad_prediction_replication_qc.csv"))
  }
  fwrite(bootstrap, file.path(cfg$paths$evaluation, "prediction_bootstrap_replicates.csv.gz"), na = "")
  importance <- fread(file.path(cfg$paths$models, "final_model_importance.csv.gz"))
  yur_write_csv(importance, file.path(cfg$paths$evaluation, "table_s12_model_importance.csv.gz"))
  yur_write_json(list(
    status = "PASS", endpoints = endpoints, endpoint_n = length(endpoints),
    bootstrap_n = cfg$bootstrap_n, interval_exports = c("IQR", "95CI"),
    prediction_panel_mode = cfg$prediction_panel_mode,
    score2_qc = "Cohort construction hard-fails on implausible lipid units or saturated SCORE2.",
    note = "Figure 4 caption reports IQR while methods report 95% CI; both are retained."
  ), file.path(cfg$paths$evaluation, "evaluation_summary.json"))
}

yur_theme <- function() {
  ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      strip.background = ggplot2::element_rect(fill = "#F2F4F3", color = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}

yur_save_plot <- function(plot, stem, cfg, width, height, directory = cfg$paths$figures) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(directory, paste0(stem, ".pdf"))
  png_path <- file.path(directory, paste0(stem, ".png"))
  tiff_path <- file.path(directory, paste0(stem, ".tiff"))
  cairo_ok <- tryCatch({
    suppressWarnings(ggplot2::ggsave(
      pdf_path, plot, width = width, height = height, units = "in",
      device = grDevices::cairo_pdf
    ))
    file.exists(pdf_path)
  }, error = function(e) FALSE)
  if (!cairo_ok) {
    if (file.exists(pdf_path)) unlink(pdf_path)
    ggplot2::ggsave(pdf_path, plot, width = width, height = height, units = "in", device = grDevices::pdf)
  }
  png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else grDevices::png
  tiff_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_tiff else grDevices::tiff
  ggplot2::ggsave(
    png_path, plot, width = width, height = height, units = "in", dpi = 400,
    device = png_device
  )
  ggplot2::ggsave(
    tiff_path, plot, width = width, height = height, units = "in", dpi = 600,
    device = tiff_device, compression = "lzw"
  )
}

yur_reference_sheet <- function(cfg, sheet) {
  path <- file.path(cfg$paths$source_tables, paste0(tolower(sheet), "_official.csv.gz"))
  if (file.exists(path)) fread(path) else NULL
}

yur_build_figures <- function(cfg) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  colors <- c(SCORE2 = "#9A9A9A", Protein = "#2C7FB8", Protein_SCORE2 = "#D95F0E")

  flow_file <- file.path(cfg$paths$cohort, "cohort_flow.csv")
  event_file <- file.path(cfg$paths$cohort, "endpoint_event_summary.csv")
  if (file.exists(flow_file)) {
    flow <- fread(flow_file)
    flow[, x := seq_len(.N)]
    flow[, label := paste0(gsub("_", " ", step), "\n", format(n, big.mark = ","))]
    flow[, fill_color := rep(c("#F4F7F6", "#EDF5F2", "#E2F0EB", "#EAF0F7", "#F7ECE8", "#EFEAF7"), length.out = .N)]
    label_layer_args <- list(
      mapping = ggplot2::aes(label = label, fill = fill_color),
      size = 3.1,
      lineheight = .95,
      label.r = grid::unit(1.2, "mm"),
      label.size = .3,
      show.legend = FALSE
    )
    flow_label_layer <- suppressWarnings(do.call(ggplot2::geom_label, label_layer_args))
    p1 <- ggplot2::ggplot(flow, ggplot2::aes(x, 1)) +
      ggplot2::geom_segment(
        data = flow[-nrow(flow)], ggplot2::aes(x = x + .34, xend = x + .66, y = 1, yend = 1),
        linewidth = .55, color = "#7A8C87", arrow = grid::arrow(length = grid::unit(2.2, "mm"))
      ) +
      flow_label_layer +
      ggplot2::scale_fill_identity() +
      ggplot2::annotate(
        "text", x = mean(range(flow$x)), y = .72,
        label = "Baseline Olink panel -> incident association -> derivation-only selection -> locked hold-out prediction",
        size = 3.1, color = "#40524D"
      ) +
      ggplot2::coord_cartesian(xlim = c(.55, max(flow$x) + .45), ylim = c(.62, 1.16), clip = "off") +
      ggplot2::labs(title = "Yu/Chen 2025 local reproduction workflow") +
      ggplot2::theme_void(base_size = 10) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = .5, margin = ggplot2::margin(b = 14)))
    yur_save_plot(p1, "figure1_workflow", cfg, 10.5, 3.1)
    yur_write_csv(flow, file.path(cfg$paths$figures, "figure1_source_data.csv"))
  }

  if (file.exists(event_file)) {
    events <- fread(event_file)
    events_long <- melt(
      events,
      id.vars = "outcome_id",
      measure.vars = c("derivation_events", "test_events"),
      variable.name = "split", value.name = "events"
    )
    events_long[, split := factor(split, levels = c("derivation_events", "test_events"), labels = c("Derivation", "Hold-out"))]
    ps1 <- ggplot2::ggplot(events_long, ggplot2::aes(events, reorder(outcome_id, events), color = split)) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = events, yend = reorder(outcome_id, events)), linewidth = .35, color = "#C5CECB") +
      ggplot2::geom_point(size = 2) +
      ggplot2::facet_wrap(~ split, scales = "free_x") +
      ggplot2::scale_color_manual(values = c(Derivation = "#287A67", `Hold-out` = "#C9523D"), guide = "none") +
      ggplot2::labs(x = "Incident events", y = NULL, title = "Incident event distribution") + yur_theme()
    yur_save_plot(ps1, "figure_s1_event_distribution", cfg, 8.5, 5.5, cfg$paths$supplement_figures)
    yur_write_csv(events_long, file.path(cfg$paths$supplement_figures, "figure_s1_source_data.csv"))
  }

  cox_file <- file.path(cfg$paths$cox, "table_s2_incident_associations.csv.gz")
  if (file.exists(cox_file)) {
    cox <- fread(cox_file)
    figure2_outcome_order <- c(
      "abdominal_aneurysm", "aortic_valve_stenosis", "atrial_fibrillation",
      "cardiomyopathy", "cad", "deep_vein_thrombosis", "heart_failure",
      "intracerebral_hemorrhage", "ischemic_stroke", "peripheral_arterial_disease",
      "pulmonary_embolism", "subarachnoid_hemorrhage", "thoracic_aneurysm",
      "transient_ischemic_attack"
    )
    figure2_outcome_labels <- c(
      abdominal_aneurysm = "Abdominal aneurysm",
      aortic_valve_stenosis = "Aortic valve stenosis",
      atrial_fibrillation = "Atrial fibrillation",
      cardiomyopathy = "Cardiomyopathy",
      cad = "Coronary artery disease",
      deep_vein_thrombosis = "Deep vein thrombosis",
      heart_failure = "Heart failure",
      intracerebral_hemorrhage = "Intracerebral hemorrhage",
      ischemic_stroke = "Ischemic stroke",
      peripheral_arterial_disease = "Peripheral arterial disease",
      pulmonary_embolism = "Pulmonary embolism",
      subarachnoid_hemorrhage = "Subarachnoid hemorrhage",
      thoracic_aneurysm = "Thoracic aneurysm",
      transient_ischemic_attack = "Transient ischemic attack"
    )
    panel_letters <- setNames(LETTERS[2:15], figure2_outcome_order)
    positive_color <- "#4F78A8"
    negative_color <- "#E5B86B"
    nonsignificant_color <- "#CFCFCF"
    risk_color <- "#C94C45"

    counts <- cox[, .(n_positive = sum(bonferroni_significant & beta > 0, na.rm = TRUE),
                      n_negative = sum(bonferroni_significant & beta < 0, na.rm = TRUE)),
                  by = .(outcome_id, outcome_label)]
    counts[, outcome_id := factor(outcome_id, levels = figure2_outcome_order)]
    setorder(counts, outcome_id)
    counts[, total := n_positive + n_negative]
    counts_long <- melt(
      counts,
      id.vars = c("outcome_id", "outcome_label"),
      measure.vars = c("n_positive", "n_negative"),
      variable.name = "direction", value.name = "count"
    )
    counts_long[, outcome_id := factor(outcome_id, levels = figure2_outcome_order)]
    # ggplot stacks the first factor level on top. Keep positive below and
    # negative above, matching the article's visual hierarchy.
    counts_long[, direction := factor(direction, levels = c("n_negative", "n_positive"))]
    p2a <- ggplot2::ggplot(counts_long, ggplot2::aes(outcome_id, count, fill = direction)) +
      ggplot2::geom_col(width = .72, color = "white", linewidth = .12) +
      ggplot2::geom_text(
        data = counts, ggplot2::aes(outcome_id, total, label = total),
        inherit.aes = FALSE, vjust = -.28, size = 2.65,
        color = "#202020", fontface = "bold"
      ) +
      ggplot2::scale_fill_manual(
        values = c(n_positive = positive_color, n_negative = negative_color),
        breaks = c("n_negative", "n_positive"),
        labels = c("Negative associated", "Positive associated")
      ) +
      ggplot2::scale_x_discrete(labels = unname(figure2_outcome_labels[figure2_outcome_order])) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .12))) +
      ggplot2::labs(x = NULL, y = "No. significant\nassociations", title = "A", fill = NULL) +
      ggplot2::theme_classic(base_size = 8.4, base_family = "Arial") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 10.5, hjust = -.08, margin = ggplot2::margin(b = 1)),
        axis.text.x = ggplot2::element_text(angle = 52, hjust = 1, vjust = 1, size = 6.8, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 7, face = "bold"),
        axis.title.y = ggplot2::element_text(size = 7.8, face = "bold", margin = ggplot2::margin(r = 2)),
        axis.line = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_line(linewidth = .25),
        panel.border = ggplot2::element_rect(color = "#202020", fill = NA, linewidth = .45),
        legend.position = "inside",
        legend.position.inside = c(.985, .965),
        legend.justification = c(1, 1),
        legend.direction = "horizontal",
        legend.background = ggplot2::element_rect(fill = scales::alpha("white", .88), color = NA),
        legend.key.size = grid::unit(2.5, "mm"),
        legend.text = ggplot2::element_text(size = 6.8, face = "bold"),
        plot.margin = ggplot2::margin(.3, .5, 0, .3)
      )
    yur_save_plot(p2a, "figure2a_significant_counts", cfg, 7.2, 3.0)

    cox[, neglog10p := -log10(pmax(p, .Machine$double.xmin))]
    cox[, direction := fifelse(bonferroni_significant & beta > 0, "Risk", fifelse(bonferroni_significant & beta < 0, "Protective", "Not significant"))]
    cox[, display_protein := fifelse(
      tolower(feature_id) == "ntprobnp", "NT-proBNP",
      fifelse(tolower(feature_id) == "nppb", "NPPB",
              fifelse(is.na(protein) | !nzchar(protein), toupper(feature_id), toupper(protein)))
    )]
    cox[, label_rank := frank(p, ties.method = "first"), by = outcome_id]
    cox[, significant_n := sum(bonferroni_significant, na.rm = TRUE), by = outcome_id]
    cox[, label := fifelse(
      bonferroni_significant & label_rank <= 10L,
      display_protein,
      NA_character_
    )]

    volcano_plots <- lapply(figure2_outcome_order, function(endpoint) {
      dat <- cox[outcome_id == endpoint]
      dat <- dat[order(direction != "Not significant")]
      letter <- panel_letters[[endpoint]]
      outcome_label <- figure2_outcome_labels[[endpoint]]
      threshold <- unique(dat$bonferroni_threshold)
      threshold <- threshold[is.finite(threshold) & threshold > 0][1]
      finite_hr <- dat$hr[is.finite(dat$hr) & dat$hr > 0]
      finite_y <- dat$neglog10p[is.finite(dat$neglog10p) & dat$neglog10p >= 0]
      hr_range <- diff(range(finite_hr))
      if (!is.finite(hr_range) || hr_range <= 0) hr_range <- max(abs(finite_hr - 1), .1)
      x_lower <- max(0, min(finite_hr) - .04 * hr_range)
      x_data_upper <- max(finite_hr) + .07 * hr_range
      # Reserve more horizontal room for the risk side while retaining every
      # observed HR. This places HR=1 near the left third whenever feasible.
      x_target_upper <- 1 + 2.15 * max(1 - x_lower, .08)
      x_upper <- max(x_data_upper, x_target_upper)
      y_upper <- max(finite_y) * 1.08
      if (!is.finite(y_upper) || y_upper <= 0) y_upper <- 1
      p <- ggplot2::ggplot(dat, ggplot2::aes(hr, neglog10p)) +
        ggplot2::geom_point(
          data = dat[direction == "Not significant"],
          color = nonsignificant_color, size = 1.08, alpha = .58, stroke = 0
        ) +
        ggplot2::geom_point(
          data = dat[direction != "Not significant"],
          ggplot2::aes(color = direction),
          size = 1.28, alpha = .92, stroke = 0
        ) +
        ggplot2::geom_vline(xintercept = 1, linewidth = .38, linetype = 2, color = "#404040") +
        ggplot2::geom_hline(yintercept = -log10(threshold), linewidth = .36, color = "#D8B77A") +
        ggplot2::scale_color_manual(
          values = c(Risk = risk_color, Protective = positive_color, `Not significant` = nonsignificant_color),
          breaks = c("Not significant", "Protective", "Risk")
        ) +
        ggplot2::scale_x_continuous(
          limits = c(x_lower, x_upper),
          expand = ggplot2::expansion(mult = c(0, 0))
        ) +
        ggplot2::scale_y_continuous(
          limits = c(0, y_upper),
          expand = ggplot2::expansion(mult = c(0, 0))
        ) +
        ggplot2::labs(
          x = outcome_label,
          y = expression(-Log[10](italic(P) * "-value")),
          title = letter,
          color = NULL
        ) +
        ggplot2::theme_classic(base_size = 8.4, base_family = "Arial") +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 10.5, hjust = -.12, margin = ggplot2::margin(b = 1)),
          axis.title.x = ggplot2::element_text(size = 7.2, face = "bold", margin = ggplot2::margin(t = 1.5)),
          axis.title.y = ggplot2::element_text(size = 7.2, face = "bold", margin = ggplot2::margin(r = 1.5)),
          axis.text = ggplot2::element_text(size = 6.8, face = "bold", color = "#202020"),
          axis.line = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_line(linewidth = .3),
          panel.border = ggplot2::element_rect(color = "#202020", fill = NA, linewidth = .45),
          aspect.ratio = 1,
          legend.position = "none",
          plot.margin = ggplot2::margin(.15, .3, .15, .3)
        )
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
          data = dat[!is.na(label)], ggplot2::aes(label = label),
          size = 2.25, color = "black", fontface = "bold",
          box.padding = .16, point.padding = .07,
          min.segment.length = 0, segment.size = .22, segment.color = "#555555",
          seed = 20250715, max.overlaps = Inf, force = .9, show.legend = FALSE
        )
      } else {
        p <- p + ggplot2::geom_text(
          data = dat[!is.na(label)], ggplot2::aes(label = label),
          size = 2.1, color = "black", fontface = "bold",
          check_overlap = TRUE, vjust = -.35, show.legend = FALSE
        )
      }
      p
    })
    names(volcano_plots) <- LETTERS[2:15]
    if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required for Figure 2.")
    p2b <- patchwork::wrap_plots(volcano_plots, ncol = 4)
    yur_save_plot(p2b, "figure2b_o_cox_volcano", cfg, 11.2, 11.2)

    figure2_plots <- c(list(A = p2a), volcano_plots)
    figure2_design <- paste(c("AABC", "DEFG", "HIJK", "LMNO"), collapse = "\n")
    p2_complete <- patchwork::wrap_plots(
      figure2_plots,
      design = figure2_design,
      widths = rep(1, 4),
      heights = rep(1, 4)
    )
    yur_save_plot(p2_complete, "figure2_incident_protein_associations", cfg, 12, 10.9)
    yur_write_csv(counts, file.path(cfg$paths$figures, "figure2_source_data_counts.csv"))
    fwrite(cox, file.path(cfg$paths$figures, "figure2_source_data_volcano.csv.gz"), na = "")

    breadth <- cox[bonferroni_significant == TRUE, .(
      associated_outcomes = uniqueN(outcome_id),
      strongest_log10p = max(neglog10p, na.rm = TRUE)
    ), by = .(feature_id, protein)]
    breadth[, protein_label := fifelse(is.na(protein) | !nzchar(protein), feature_id, protein)]
    breadth <- breadth[order(-associated_outcomes, -strongest_log10p)]
    ps3 <- ggplot2::ggplot(breadth, ggplot2::aes(associated_outcomes, reorder(protein_label, associated_outcomes))) +
      ggplot2::geom_point(ggplot2::aes(size = strongest_log10p, color = associated_outcomes), alpha = .82) +
      ggplot2::scale_color_viridis_c(option = "C", direction = -1) +
      ggplot2::labs(x = "Number of associated CVD outcomes", y = NULL, size = expression(max(-log[10](P))), color = "Outcomes", title = "Cross-CVD association breadth") +
      yur_theme() + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
    yur_save_plot(ps3, "figure_s3_multimorbidity_breadth", cfg, 7.5, 6.2, cfg$paths$supplement_figures)
    yur_write_csv(breadth, file.path(cfg$paths$supplement_figures, "figure_s3_source_data.csv"))
  }

  metrics_file <- file.path(cfg$paths$evaluation, "table_s10_prediction_metrics.csv")
  importance_file <- file.path(cfg$paths$evaluation, "table_s12_model_importance.csv.gz")
  model_design_file <- file.path(cfg$paths$models, "model_design_contract.csv")
  article_outcome_order <- c(
    "abdominal_aneurysm", "aortic_valve_stenosis", "atrial_fibrillation",
    "cardiomyopathy", "cad", "deep_vein_thrombosis", "heart_failure",
    "intracerebral_hemorrhage", "ischemic_stroke", "peripheral_arterial_disease",
    "pulmonary_embolism", "subarachnoid_hemorrhage", "thoracic_aneurysm",
    "transient_ischemic_attack"
  )
  article_outcome_labels <- c(
    abdominal_aneurysm = "Abdominal aneurysm",
    aortic_valve_stenosis = "Aortic valve stenosis",
    atrial_fibrillation = "Atrial fibrillation",
    cardiomyopathy = "Cardiomyopathy",
    cad = "Coronary artery disease",
    deep_vein_thrombosis = "Deep vein thrombosis",
    heart_failure = "Heart failure",
    intracerebral_hemorrhage = "Intracerebral hemorrhage",
    ischemic_stroke = "Ischemic stroke",
    peripheral_arterial_disease = "Peripheral arterial disease",
    pulmonary_embolism = "Pulmonary embolism",
    subarachnoid_hemorrhage = "Subarachnoid hemorrhage",
    thoracic_aneurysm = "Thoracic aneurysm",
    transient_ischemic_attack = "Transient ischemic attack"
  )
  article_outcome_colors <- c(
    abdominal_aneurysm = "#264B6A", aortic_valve_stenosis = "#3B7198",
    atrial_fibrillation = "#70A4C2", cardiomyopathy = "#B9D6D2",
    cad = "#2E6E62", deep_vein_thrombosis = "#75B7A5",
    heart_failure = "#789A5A", intracerebral_hemorrhage = "#BE4B43",
    ischemic_stroke = "#D67848", peripheral_arterial_disease = "#E9A34F",
    pulmonary_embolism = "#DF7084", subarachnoid_hemorrhage = "#E9A4B6",
    thoracic_aneurysm = "#65537D", transient_ischemic_attack = "#A99CC0"
  )
  article_model_colors <- c(
    "Proteins + SCORE2" = "#D96652", Proteins = "#4B78A6", SCORE2 = "#E5C36B"
  )
  parse_figure4_extras <- function(x) {
    x <- trimws(as.character(x %||% ""))
    if (!nzchar(x)) return(character())
    trimws(strsplit(x, ";", fixed = TRUE)[[1]])
  }
  figure4_extra_projects <- parse_figure4_extras(cfg$figure4_extra_projects)
  figure4_extra_outcomes <- parse_figure4_extras(cfg$figure4_extra_outcomes)
  figure4_extra_labels <- parse_figure4_extras(cfg$figure4_extra_labels)
  for (values in list(figure4_extra_outcomes, figure4_extra_labels)) {
    if (length(values) && length(values) != length(figure4_extra_projects)) {
      stop("Figure 4 extra outcome/label counts must match the extra project count.")
    }
  }
  figure4_parts <- list()
  if (file.exists(metrics_file)) {
    if (!file.exists(model_design_file)) {
      stop("Figure 4 requires the model design contract: ", model_design_file)
    }
    model_design <- fread(model_design_file)
    required_design_columns <- c("outcome_id", "model_id", "feature_n")
    missing_design_columns <- setdiff(required_design_columns, names(model_design))
    if (length(missing_design_columns)) {
      stop(
        "Figure 4 model design contract lacks required columns: ",
        paste(missing_design_columns, collapse = ", ")
      )
    }
    model_design[, feature_n := suppressWarnings(as.integer(feature_n))]
    protein_design <- model_design[model_id == "Protein"]
    if (!nrow(protein_design) || anyNA(protein_design$feature_n) || any(protein_design$feature_n < 1L)) {
      stop("Figure 4 could not resolve valid protein feature counts from the Protein model rows.")
    }
    inconsistent_counts <- protein_design[, .(n_counts = uniqueN(feature_n)), by = outcome_id][n_counts != 1L]
    if (nrow(inconsistent_counts)) {
      stop(
        "Figure 4 found more than one Protein feature count for outcome(s): ",
        paste(inconsistent_counts$outcome_id, collapse = ", ")
      )
    }
    protein_counts <- unique(protein_design[, .(
      outcome_id,
      protein_n = feature_n
    )])

    metrics <- fread(metrics_file)
    metrics_public <- metrics[model_id %in% c("Protein_SCORE2", "Protein", "SCORE2")]
    main_outcomes <- article_outcome_order[article_outcome_order %chin% unique(metrics_public$outcome_id)]
    if (!length(main_outcomes)) stop("Figure 4 found no configured outcomes in the prediction metrics.")
    protein_counts <- protein_counts[outcome_id %chin% main_outcomes]
    comparison_tables <- list()
    comparison_path <- file.path(cfg$paths$evaluation, "table_s11_model_comparisons_nri_idi.csv")
    if (file.exists(comparison_path)) comparison_tables[[1L]] <- fread(comparison_path)
    extra_manifest <- list()
    extra_colors <- c("#8A6A9B", "#B56576", "#5D8A82", "#A67845", "#536B90", "#7A8F45")
    for (extra_index in seq_along(figure4_extra_projects)) {
      project_spec <- figure4_extra_projects[[extra_index]]
      extra_dir <- if (grepl("^([A-Za-z]:[/\\\\]|/)", project_spec)) {
        normalizePath(project_spec, winslash = "/", mustWork = FALSE)
      } else {
        file.path(cfg$analysis_root, project_spec)
      }
      extra_metrics_path <- file.path(extra_dir, "10_evaluation", "table_s10_prediction_metrics.csv")
      extra_design_path <- file.path(extra_dir, "09_models", "model_design_contract.csv")
      if (!file.exists(extra_metrics_path) || !file.exists(extra_design_path)) {
        stop(
          "Figure 4 extra project is incomplete: ", project_spec,
          ". Required: ", extra_metrics_path, " and ", extra_design_path
        )
      }
      extra_metrics <- fread(extra_metrics_path)[
        model_id %in% c("Protein_SCORE2", "Protein", "SCORE2")
      ]
      available_extra_outcomes <- unique(extra_metrics$outcome_id)
      requested_outcome <- if (length(figure4_extra_outcomes)) {
        figure4_extra_outcomes[[extra_index]]
      } else if (length(available_extra_outcomes) == 1L) {
        available_extra_outcomes[[1L]]
      } else {
        stop(
          "Figure 4 extra project ", project_spec,
          " contains multiple outcomes; specify -Figure4ExtraOutcome."
        )
      }
      if (!requested_outcome %chin% available_extra_outcomes) {
        stop(
          "Figure 4 extra outcome ", requested_outcome, " is absent from ",
          project_spec, ". Available: ", paste(available_extra_outcomes, collapse = ", ")
        )
      }
      display_id <- paste0(
        "extra_", yur_norm_name(basename(extra_dir)), "_", yur_norm_name(requested_outcome)
      )
      if (display_id %chin% c(main_outcomes, vapply(extra_manifest, `[[`, character(1), "display_id"))) {
        stop("Figure 4 extra display ID collision: ", display_id)
      }
      display_label <- if (length(figure4_extra_labels)) {
        figure4_extra_labels[[extra_index]]
      } else {
        paste0(article_outcome_labels[[requested_outcome]] %||% requested_outcome, " (custom panel)")
      }
      extra_metrics <- extra_metrics[outcome_id == requested_outcome]
      extra_metrics[, outcome_id := display_id]
      metrics_public <- rbind(metrics_public, extra_metrics, fill = TRUE)

      extra_design <- fread(extra_design_path)
      extra_count <- unique(extra_design[
        outcome_id == requested_outcome & model_id == "Protein",
        suppressWarnings(as.integer(feature_n))
      ])
      extra_count <- extra_count[is.finite(extra_count) & extra_count > 0]
      if (length(extra_count) != 1L) {
        stop("Figure 4 could not resolve one protein count for ", project_spec, "/", requested_outcome)
      }
      protein_counts <- rbind(
        protein_counts,
        data.table(outcome_id = display_id, protein_n = extra_count[[1L]])
      )
      extra_comparison_path <- file.path(
        extra_dir, "10_evaluation", "table_s11_model_comparisons_nri_idi.csv"
      )
      if (file.exists(extra_comparison_path)) {
        extra_comparisons <- fread(extra_comparison_path)[outcome_id == requested_outcome]
        extra_comparisons[, outcome_id := display_id]
        comparison_tables[[length(comparison_tables) + 1L]] <- extra_comparisons
      }
      extra_manifest[[length(extra_manifest) + 1L]] <- list(
        display_id = display_id, display_label = display_label,
        source_project = project_spec, source_analysis_dir = extra_dir,
        source_outcome_id = requested_outcome,
        protein_n = extra_count[[1L]],
        metrics_sha256 = yur_sha_file(extra_metrics_path),
        design_sha256 = yur_sha_file(extra_design_path)
      )
    }
    extra_ids <- if (length(extra_manifest)) {
      vapply(extra_manifest, `[[`, character(1), "display_id")
    } else character()
    figure4a_outcome_order <- c(main_outcomes, extra_ids)
    figure4a_outcome_labels <- article_outcome_labels[main_outcomes]
    figure4a_outcome_colors <- article_outcome_colors[main_outcomes]
    if (length(extra_manifest)) {
      extra_color_values <- rep(extra_colors, length.out = length(extra_ids))
      names(extra_color_values) <- extra_ids
      figure4a_outcome_labels <- c(
        figure4a_outcome_labels,
        setNames(vapply(extra_manifest, `[[`, character(1), "display_label"), extra_ids)
      )
      figure4a_outcome_colors <- c(
        figure4a_outcome_colors,
        extra_color_values
      )
      yur_write_csv(
        rbindlist(lapply(extra_manifest, as.data.table), fill = TRUE),
        file.path(cfg$paths$figures, "figure4a_extra_project_manifest.csv")
      )
    }
    metrics_public[, outcome_index := match(outcome_id, figure4a_outcome_order)]
    metrics_public[, model_label := factor(
      model_id,
      levels = c("Protein_SCORE2", "Protein", "SCORE2"),
      labels = c("Proteins + SCORE2", "Proteins", "SCORE2")
    )]
    metrics_public <- metrics_public[!is.na(outcome_index)]
    metrics_public[, label_base := max(auc_q75, na.rm = TRUE) + .034, by = outcome_id]
    metrics_public[, label_y := pmin(
      label_base + fifelse(model_id == "Protein", .038, 0),
      1.055
    )]
    strip_data <- data.table(
      outcome_id = figure4a_outcome_order,
      outcome_index = seq_along(figure4a_outcome_order)
    )
    strip_data <- merge(strip_data, protein_counts, by = "outcome_id", all.x = TRUE, sort = FALSE)
    setorder(strip_data, outcome_index)
    if (anyNA(strip_data$protein_n)) {
      stop(
        "Figure 4 lacks Protein feature counts for outcome(s): ",
        paste(strip_data[is.na(protein_n), outcome_id], collapse = ", ")
      )
    }
    strip_data[, `:=`(
      outcome_id = factor(outcome_id, levels = figure4a_outcome_order),
      protein_n_label = sprintf("n=%d", protein_n),
      strip_fill = unname(figure4a_outcome_colors[as.character(outcome_id)])
    )]
    strip_rgb <- t(grDevices::col2rgb(strip_data$strip_fill)) / 255
    strip_rgb_linear <- ifelse(
      strip_rgb <= .03928,
      strip_rgb / 12.92,
      ((strip_rgb + .055) / 1.055)^2.4
    )
    strip_luminance <-
      .2126 * strip_rgb_linear[, 1] +
      .7152 * strip_rgb_linear[, 2] +
      .0722 * strip_rgb_linear[, 3]
    strip_data[, strip_text_color := fifelse(strip_luminance > .42, "dark", "light")]

    brackets <- data.table()
    if (length(comparison_tables)) {
      comparisons <- rbindlist(comparison_tables, fill = TRUE)[delong_p < .05]
      bracket_map <- data.table(
        comparison_id = c("Protein_SCORE2_vs_Protein", "Protein_vs_SCORE2", "Protein_SCORE2_vs_SCORE2"),
        x_offset_start = c(-.24, 0, -.24), x_offset_end = c(0, .24, .24), bracket_rank = c(1L, 2L, 3L)
      )
      brackets <- merge(comparisons, bracket_map, by = "comparison_id")
      brackets[, outcome_index := match(outcome_id, figure4a_outcome_order)]
      brackets <- brackets[!is.na(outcome_index)]
      upper <- metrics_public[, .(upper = max(auc_q75, na.rm = TRUE)), by = outcome_id]
      brackets <- merge(brackets, upper, by = "outcome_id")
      brackets[, `:=`(
        x_start = outcome_index + x_offset_start,
        x_end = outcome_index + x_offset_end,
        y = pmin(1.085, upper + .102 + (bracket_rank - 1L) * .035),
        star = fifelse(delong_p < .001, "***", fifelse(delong_p < .01, "**", "*"))
      )]
    }

    dodge <- ggplot2::position_dodge(width = .72)
    p4a <- ggplot2::ggplot(metrics_public, ggplot2::aes(outcome_index, auc_q50, fill = model_label, group = model_label)) +
      ggplot2::geom_col(position = dodge, width = .66, color = "white", linewidth = .18) +
      ggplot2::geom_rect(
        data = strip_data,
        ggplot2::aes(
          xmin = outcome_index - .49, xmax = outcome_index + .49,
          ymin = -.050, ymax = -.010
        ),
        inherit.aes = FALSE,
        fill = strip_data$strip_fill,
        color = "#333333", linewidth = .18
      ) +
      ggplot2::geom_text(
        data = strip_data[strip_text_color == "light"],
        ggplot2::aes(x = outcome_index, y = -.030, label = protein_n_label),
        inherit.aes = FALSE, color = "#FFFFFF", family = "Arial",
        fontface = "bold", size = 2.15
      ) +
      ggplot2::geom_text(
        data = strip_data[strip_text_color == "dark"],
        ggplot2::aes(x = outcome_index, y = -.030, label = protein_n_label),
        inherit.aes = FALSE, color = "#171717", family = "Arial",
        fontface = "bold", size = 2.15
      ) +
      ggplot2::geom_hline(yintercept = 0, linewidth = .34, color = "#222222") +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = auc_q25, ymax = auc_q75), position = dodge,
        width = .11, linewidth = .32, color = "#333333"
      ) +
      ggplot2::geom_text(
        ggplot2::aes(y = label_y, label = sprintf("%.2f", auc_q50), color = model_label),
        position = dodge, vjust = 0, size = 2.65, fontface = "bold", show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(values = article_model_colors, drop = FALSE) +
      ggplot2::scale_color_manual(values = article_model_colors, drop = FALSE) +
      ggplot2::scale_x_continuous(breaks = seq_along(figure4a_outcome_order), labels = NULL, expand = c(.015, .015)) +
      ggplot2::scale_y_continuous(
        breaks = seq(0, 1, .2), limits = c(-.055, 1.13), expand = c(0, 0)
      ) +
      ggplot2::labs(x = NULL, y = "AUC", fill = NULL, tag = "A") +
      ggplot2::theme_classic(base_size = 9.2, base_family = "Arial") +
      ggplot2::theme(
        axis.ticks.x = ggplot2::element_blank(), legend.position = "inside",
        legend.position.inside = c(.86, .96),
        legend.direction = "horizontal", legend.justification = c(1, 1),
        axis.title.y = ggplot2::element_text(face = "bold", size = 9.4),
        axis.text.y = ggplot2::element_text(face = "bold", size = 8.2),
        panel.border = ggplot2::element_rect(fill = NA, color = "#222222", linewidth = .42),
        legend.key.width = grid::unit(4.8, "mm"),
        legend.text = ggplot2::element_text(face = "bold", size = 7.8),
        plot.tag = ggplot2::element_text(face = "bold", size = 11),
        plot.tag.position = c(.005, .99), plot.margin = ggplot2::margin(3, 4, 0, 4)
      )
    if (nrow(brackets)) {
      p4a <- p4a +
        ggplot2::geom_segment(
          data = brackets, ggplot2::aes(x = x_start, xend = x_end, y = y, yend = y),
          inherit.aes = FALSE, linewidth = .28, color = "#343434"
        ) +
        ggplot2::geom_segment(
          data = brackets, ggplot2::aes(x = x_start, xend = x_start, y = y, yend = y - .012),
          inherit.aes = FALSE, linewidth = .28, color = "#343434"
        ) +
        ggplot2::geom_segment(
          data = brackets, ggplot2::aes(x = x_end, xend = x_end, y = y, yend = y - .012),
          inherit.aes = FALSE, linewidth = .28, color = "#343434"
        ) +
        ggplot2::geom_text(
          data = brackets, ggplot2::aes(x = (x_start + x_end) / 2, y = y + .008, label = star),
          inherit.aes = FALSE, size = 2.7, fontface = "bold", vjust = 0, color = "#222222"
        )
    }

    legend_data <- data.table(
      outcome_id = factor(figure4a_outcome_order, levels = figure4a_outcome_order),
      label = unname(figure4a_outcome_labels[figure4a_outcome_order])
    )
    legend_n_rows <- max(1L, ceiling(nrow(legend_data) / 5L))
    legend_data[, `:=`(
      legend_row = legend_n_rows - ((seq_len(.N) - 1L) %/% 5L),
      legend_col = 1L + ((seq_len(.N) - 1L) %% 5L)
    )]
    legend_data[, item_width := .42 + nchar(label) * .057]
    legend_col_widths <- legend_data[, .(col_width = max(item_width)), by = legend_col][order(legend_col)]
    legend_col_widths[, x := cumsum(data.table::shift(col_width, fill = 0))]
    legend_width <- legend_col_widths[, max(x + col_width)]
    legend_data <- merge(legend_data, legend_col_widths, by = "legend_col", sort = FALSE)
    legend_data[, y := legend_row - .35]

    p4legend <- ggplot2::ggplot(legend_data, ggplot2::aes(x, y, color = outcome_id)) +
      ggplot2::annotate(
        "text", x = legend_width / 2, y = legend_n_rows + .78, label = "Incident CVD outcomes",
        family = "Arial", fontface = "bold", size = 3.0
      ) +
      ggplot2::geom_point(shape = 15, size = 3.6) +
      ggplot2::geom_text(
        ggplot2::aes(x = x + .11, label = label),
        hjust = 0, color = "#222222", family = "Arial", fontface = "bold", size = 2.65
      ) +
      ggplot2::scale_color_manual(values = figure4a_outcome_colors, guide = "none") +
      ggplot2::coord_cartesian(
        xlim = c(-.05, legend_width + .05), ylim = c(.25, legend_n_rows + 1.55), clip = "off"
      ) +
      ggplot2::theme_void(base_family = "Arial") +
      ggplot2::theme(plot.margin = ggplot2::margin(0, 3, 0, 3))

    yur_save_plot(p4a, "figure4a_prediction_auc", cfg, 10.8, 4.0)
    figure4a_source <- merge(metrics_public, protein_counts, by = "outcome_id", all.x = TRUE, sort = FALSE)
    yur_write_csv(figure4a_source, file.path(cfg$paths$figures, "figure4a_source_data.csv"))
    yur_write_csv(
      strip_data[, .(outcome_id, outcome_index, protein_n, protein_n_label)],
      file.path(cfg$paths$figures, "figure4a_model_protein_counts.csv")
    )
    figure4_parts$legend <- p4legend
    figure4_parts$a <- p4a
  }
  if (file.exists(importance_file)) {
    importance <- fread(importance_file)[model_id == "Protein"]
    importance[, rank := frank(-standardized_gain, ties.method = "first"), by = outcome_id]
    top15 <- importance[rank <= 15]
    recurrent <- top15[, .(times_top15 = uniqueN(outcome_id)), by = feature][times_top15 >= 2]
    repeated <- merge(top15, recurrent, by = "feature")
    panel_map <- fread(file.path(cfg$paths$cox, "retained_panel.csv"))[, .(feature = feature_id, protein)]
    repeated <- merge(repeated, panel_map, by = "feature", all.x = TRUE)
    repeated[, protein_label := fifelse(
      tolower(feature) == "ntprobnp", "NT-proBNP",
      fifelse(tolower(feature) == "nppb", "NPPB",
              toupper(fifelse(is.na(protein) | !nzchar(protein), feature, protein)))
    )]
    duplicate_labels <- unique(repeated[, .(feature, protein_label)])[, .N, by = protein_label][N > 1L, protein_label]
    repeated[protein_label %in% duplicate_labels, protein_label := paste0(protein_label, " (", toupper(feature), ")")]
    repeated[, overall_importance := sum(standardized_gain, na.rm = TRUE), by = feature]
    protein_order <- unique(repeated[, .(overall_importance = unique(overall_importance)), by = protein_label][order(-overall_importance), protein_label])
    repeated[, protein_label := factor(protein_label, levels = protein_order)]
    repeated[, outcome_id := factor(outcome_id, levels = article_outcome_order)]
    totals <- repeated[, .(
      overall_importance = unique(overall_importance), times_top15 = unique(times_top15)
    ), by = protein_label][order(-overall_importance)]
    repeated[, importance_panel := factor(
      fifelse(as.character(protein_label) == "GDF15", "GDF15", "Other proteins"),
      levels = c("GDF15", "Other proteins")
    )]
    totals[, importance_panel := factor(
      fifelse(as.character(protein_label) == "GDF15", "GDF15", "Other proteins"),
      levels = c("GDF15", "Other proteins")
    )]
    importance_y_max <- max(totals$overall_importance) * 1.12
    importance_right_axis <- data.table(
      importance_panel = factor("Other proteins", levels = c("GDF15", "Other proteins")),
      axis_x = .50
    )
    importance_right_ticks <- data.table(
      importance_panel = factor("Other proteins", levels = c("GDF15", "Other proteins")),
      tick_y = pretty(c(0, importance_y_max), n = 5)
    )[tick_y >= 0 & tick_y <= importance_y_max]
    importance_right_ticks[, tick_label := sprintf("%.1f", tick_y)]
    b_theme <- ggplot2::theme_classic(base_size = 8.4, base_family = "Arial") + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 52, hjust = 1, vjust = 1, size = 7.1, face = "bold"
      ),
      axis.text.y = ggplot2::element_text(size = 7.5, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 8.5, face = "bold"),
      axis.title.x = ggplot2::element_blank(), legend.position = "none",
      strip.text = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      panel.spacing.x = grid::unit(6.0, "mm"),
      plot.margin = ggplot2::margin(2, 0, 3, 0)
    )
    p4b <- ggplot2::ggplot(repeated, ggplot2::aes(protein_label, standardized_gain, fill = outcome_id)) +
      ggplot2::geom_col(width = .82) +
      ggplot2::geom_vline(
        data = importance_right_axis, ggplot2::aes(xintercept = axis_x),
        linewidth = .42, color = "#222222"
      ) +
      ggplot2::geom_segment(
        data = importance_right_ticks,
        ggplot2::aes(x = .50, xend = .58, y = tick_y, yend = tick_y),
        inherit.aes = FALSE, linewidth = .34, color = "#222222"
      ) +
      ggplot2::geom_text(
        data = importance_right_ticks,
        ggplot2::aes(x = .45, y = tick_y, label = tick_label),
        inherit.aes = FALSE, hjust = 1, size = 2.15,
        color = "#222222", fontface = "bold"
      ) +
      ggplot2::geom_text(
        data = totals, ggplot2::aes(x = protein_label, y = overall_importance, label = times_top15),
        inherit.aes = FALSE, vjust = -.35, size = 2.45, fontface = "bold"
      ) +
      ggplot2::facet_grid(
        cols = ggplot2::vars(importance_panel), scales = "free_x", space = "free_x"
      ) +
      ggplot2::scale_fill_manual(values = article_outcome_colors, drop = FALSE) +
      ggplot2::scale_x_discrete(drop = TRUE, expand = ggplot2::expansion(add = c(.5, .5))) +
      ggplot2::scale_y_continuous(limits = c(0, importance_y_max), expand = c(0, 0)) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(y = "Stacked protein importance", tag = "B") + b_theme +
      ggplot2::theme(
        plot.tag = ggplot2::element_text(face = "bold", size = 11),
        plot.tag.position = c(-.022, 1.018),
        plot.margin = ggplot2::margin(4, 1, 3, 5)
      )
    yur_save_plot(p4b, "figure4b_recurrent_importance", cfg, 10.8, 4.3)
    yur_write_csv(repeated, file.path(cfg$paths$figures, "figure4b_source_data.csv"))
    figure4_parts$b <- p4b
  }
  if (all(c("legend", "a", "b") %in% names(figure4_parts))) {
    if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required for article-style Figure 4.")
    p4_complete <- patchwork::wrap_plots(
      list(figure4_parts$legend, figure4_parts$a, figure4_parts$b),
      ncol = 1, heights = c(.16, 1.02, .88)
    )
    yur_save_plot(p4_complete, "figure4_prediction_and_importance", cfg, 11, 9.1)
  }

  local_cmr_file <- file.path(cfg$paths$cmr, "cmr_associations.csv.gz")
  cmr_source_mode <- if (file.exists(local_cmr_file)) "local_complete_matrix" else "official_s9_reference"
  s9 <- if (cmr_source_mode == "local_complete_matrix") {
    fread(local_cmr_file)
  } else {
    yur_reference_sheet(cfg, "S9")
  }
  if (!is.null(s9)) {
    if (!"p" %in% names(s9)) {
      p_col <- intersect(c("P value", "P.value"), names(s9))
      if (!length(p_col)) stop("Figure 3 association source lacks a P-value column.")
      setnames(s9, old = p_col[[1]], new = "p")
    }
    if (!"beta" %in% names(s9)) {
      beta_col <- intersect(c("β", "Beta"), names(s9))
      if (!length(beta_col)) stop("Figure 3 association source lacks a beta column.")
      setnames(s9, old = beta_col[[1]], new = "beta")
    }
    if (!"Protein" %in% names(s9) && "protein" %in% names(s9)) setnames(s9, "protein", "Protein")
    s9[, `:=`(p = as.numeric(p), beta = as.numeric(beta))]
    s9 <- s9[is.finite(p) & is.finite(beta)]

    cmr_outcome_order <- c(
      "LVEDV (mL)", "LVESV (mL)", "LVSV (mL)", "LVEF (%)", "LVCO (L/min)",
      "LVM (g)", "Global wall thickness for the LV mean myocardial (mm)",
      "RVEDV (mL)", "RVESV (mL)", "RVSV (mL)", "RVEF (%)",
      "LAV max (mL)", "LAV min (mL)", "LASV (mL)", "LAEF (%)",
      "RAV max (mL)", "RAV min (mL)", "RASV (mL)", "RAEF (%)"
    )
    cmr_outcome_labels <- c(
      "LVEDV", "LVESV", "LVSV", "LVEF", "LVCO", "LVM", "Global wall thickness",
      "RVEDV", "RVESV", "RVSV", "RVEF", "LAVmax", "LAVmin", "LASV", "LAEF",
      "RAVmax", "RAVmin", "RASV", "RAEF"
    )
    names(cmr_outcome_labels) <- cmr_outcome_order
    cmr_anatomy <- c(rep("LV", 7), rep("RV", 4), rep("LA", 4), rep("RA", 4))
    names(cmr_anatomy) <- cmr_outcome_order
    anatomy_colors <- c(LV = "#32345D", RV = "#D45B54", LA = "#4FA6B7", RA = "#7B68A7")
    direction_colors <- c(
      `Non-significant` = "#C9C9C9",
      Negative = "#4F78A8",
      Positive = "#D45B54"
    )
    cmr_protein_n <- if ("feature_id" %in% names(s9)) uniqueN(s9$feature_id) else uniqueN(s9$Protein)
    threshold_values <- if ("bonferroni_threshold" %in% names(s9)) {
      unique(as.numeric(s9$bonferroni_threshold[is.finite(s9$bonferroni_threshold)]))
    } else numeric()
    cmr_threshold <- if (length(threshold_values) == 1L) threshold_values[[1]] else .05 / cmr_protein_n

    s9 <- s9[Outcome %in% cmr_outcome_order]
    s9[, outcome_idx := match(Outcome, cmr_outcome_order)]
    s9[, protein_idx := frank(Protein, ties.method = "dense"), by = Outcome]
    s9[, proteins_in_outcome := uniqueN(Protein), by = Outcome]
    s9[, manhattan_x := outcome_idx +
         (protein_idx - (proteins_in_outcome + 1) / 2) / proteins_in_outcome * .78]
    s9[, neglog10p := -log10(pmax(p, .Machine$double.xmin))]
    s9[, significant := p < cmr_threshold]
    s9[, effect_direction := fcase(
      !significant, "Non-significant",
      beta < 0, "Negative",
      default = "Positive"
    )]
    s9[, effect_direction := factor(
      effect_direction,
      levels = c("Non-significant", "Negative", "Positive")
    )]
    s9[, effect_size_bin := cut(
      abs(beta),
      breaks = c(-Inf, .015, .04, .08, Inf),
      labels = c("|beta| <= 0.015", "0.015 < |beta| <= 0.04", "0.04 < |beta| <= 0.08", "|beta| > 0.08"),
      right = TRUE
    )]
    setorder(s9, significant)

    max_cmr_y <- max(s9$neglog10p, na.rm = TRUE) * 1.06
    cmr_band_height <- max_cmr_y * .026
    cmr_y_breaks <- pretty(c(0, max_cmr_y), n = 5)
    cmr_y_breaks <- cmr_y_breaks[cmr_y_breaks >= 0 & cmr_y_breaks <= max_cmr_y]
    anatomy_legend <- data.table(
      anatomy = factor(names(anatomy_colors), levels = names(anatomy_colors)),
      x = -Inf,
      y = -Inf
    )
    p3a_top <- ggplot2::ggplot(
      s9,
      ggplot2::aes(manhattan_x, neglog10p, fill = effect_direction, size = effect_size_bin)
    ) +
      ggplot2::annotate("rect", xmin = .5, xmax = 7.5, ymin = -cmr_band_height, ymax = 0, fill = anatomy_colors[["LV"]]) +
      ggplot2::annotate("rect", xmin = 7.5, xmax = 11.5, ymin = -cmr_band_height, ymax = 0, fill = anatomy_colors[["RV"]]) +
      ggplot2::annotate("rect", xmin = 11.5, xmax = 15.5, ymin = -cmr_band_height, ymax = 0, fill = anatomy_colors[["LA"]]) +
      ggplot2::annotate("rect", xmin = 15.5, xmax = 19.5, ymin = -cmr_band_height, ymax = 0, fill = anatomy_colors[["RA"]]) +
      ggplot2::geom_hline(yintercept = 0, linewidth = .28, color = "#303030") +
      ggplot2::geom_hline(
        yintercept = -log10(cmr_threshold), linetype = 2,
        linewidth = .42, color = "#C7433E"
      ) +
      ggplot2::geom_vline(
        xintercept = c(7.5, 11.5, 15.5), linetype = 2,
        linewidth = .34, color = "#777777"
      ) +
      ggplot2::geom_point(shape = 21, color = "white", stroke = .34, alpha = .94) +
      ggplot2::geom_point(
        data = anatomy_legend,
        ggplot2::aes(x = x, y = y, color = anatomy),
        inherit.aes = FALSE, shape = 15, size = 0, alpha = 0, na.rm = TRUE
      ) +
      ggplot2::scale_fill_manual(values = direction_colors, drop = FALSE) +
      ggplot2::scale_color_manual(values = anatomy_colors, drop = FALSE) +
      ggplot2::scale_size_manual(
        values = c(1.7, 2.2, 2.85, 3.6),
        drop = FALSE
      ) +
      ggplot2::scale_x_continuous(
        limits = c(.5, 19.5), breaks = seq_along(cmr_outcome_order),
        labels = NULL, expand = c(0, 0)
      ) +
      ggplot2::scale_y_continuous(
        limits = c(-cmr_band_height, max_cmr_y), breaks = cmr_y_breaks,
        expand = c(0, 0)
      ) +
      ggplot2::labs(
        x = NULL,
        y = expression(-Log[10](italic(P) * "-value")),
        title = "A",
        size = "Effect size",
        fill = "Effect direction",
        color = "Anatomy"
      ) +
      ggplot2::guides(
        size = ggplot2::guide_legend(order = 1, nrow = 1),
        fill = ggplot2::guide_legend(
          order = 2, nrow = 1,
          override.aes = list(shape = 21, color = "white", size = 3, alpha = 1),
          theme = ggplot2::theme(legend.margin = ggplot2::margin(l = 5.5))
        ),
        color = ggplot2::guide_legend(
          order = 3, nrow = 1,
          override.aes = list(shape = 15, size = 3.3, alpha = 1),
          theme = ggplot2::theme(legend.margin = ggplot2::margin(l = 5.5))
        )
      ) +
      ggplot2::theme_classic(base_size = 9.2, base_family = "Arial") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = -.035),
        axis.title.y = ggplot2::element_text(
          face = "bold", size = 9.4, margin = ggplot2::margin(r = 4)
        ),
        axis.text.y = ggplot2::element_text(
          face = "bold", size = 8.3, margin = ggplot2::margin(r = .6)
        ),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(fill = NA, color = "#222222", linewidth = .45),
        legend.position = "top",
        legend.box = "horizontal",
        legend.title = ggplot2::element_text(face = "bold", size = 8.8),
        legend.text = ggplot2::element_text(size = 7.7, face = "bold"),
        legend.key.height = grid::unit(3.2, "mm"),
        legend.key.width = grid::unit(3.2, "mm"),
        legend.spacing.x = grid::unit(1.2, "mm"),
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.box.margin = ggplot2::margin(0, 0, 0, 0),
        plot.margin = ggplot2::margin(.4, .8, 0, .2)
      )

    cmr_counts <- s9[significant == TRUE, .(
      positive = sum(beta > 0),
      negative = sum(beta < 0)
    ), by = .(Outcome, outcome_idx)]
    cmr_counts <- merge(
      data.table(Outcome = cmr_outcome_order, outcome_idx = seq_along(cmr_outcome_order)),
      cmr_counts,
      by = c("Outcome", "outcome_idx"), all.x = TRUE
    )
    cmr_counts[is.na(positive), positive := 0L]
    cmr_counts[is.na(negative), negative := 0L]
    cmr_counts[, total := positive + negative]
    cmr_counts[, `:=`(
      positive_ymin = -positive,
      positive_ymax = 0,
      negative_ymin = -(positive + negative),
      negative_ymax = -positive
    )]
    count_limit <- max(cmr_counts$total) * 1.2
    p3a_counts <- ggplot2::ggplot(cmr_counts) +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = outcome_idx - .28, xmax = outcome_idx + .28,
          ymin = positive_ymin, ymax = positive_ymax
        ),
        fill = direction_colors[["Positive"]], color = "white", linewidth = .1
      ) +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = outcome_idx - .28, xmax = outcome_idx + .28,
          ymin = negative_ymin, ymax = negative_ymax
        ),
        fill = direction_colors[["Negative"]], color = "white", linewidth = .1
      ) +
      ggplot2::geom_text(
        ggplot2::aes(outcome_idx, -total, label = total),
        vjust = 1.25, size = 2.55, fontface = "bold", color = "#242424"
      ) +
      ggplot2::geom_vline(
        xintercept = c(7.5, 11.5, 15.5), linetype = 2,
        linewidth = .34, color = "#777777"
      ) +
      ggplot2::scale_x_continuous(
        limits = c(.5, 19.5), breaks = seq_along(cmr_outcome_order),
        labels = cmr_outcome_labels[cmr_outcome_order], expand = c(0, 0)
      ) +
      ggplot2::scale_y_continuous(
        limits = c(-count_limit, 0),
        breaks = -c(0, 150, 300), labels = c("0", "150", "300"),
        expand = c(0, 0)
      ) +
      ggplot2::labs(x = NULL, y = "No. significant\nassociations") +
      ggplot2::theme_classic(base_size = 8.3, base_family = "Arial") +
      ggplot2::theme(
        axis.title.y = ggplot2::element_text(
          face = "bold", size = 8.1, margin = ggplot2::margin(r = 4)
        ),
        axis.text.y = ggplot2::element_text(face = "bold", size = 7.1),
        axis.text.x = ggplot2::element_text(angle = 58, hjust = 1, vjust = 1, face = "bold", size = 6.8),
        axis.line = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_line(linewidth = .3),
        panel.border = ggplot2::element_rect(fill = NA, color = "#222222", linewidth = .45),
        plot.margin = ggplot2::margin(0, .8, .5, .2)
      )

    if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required for Figure 3.")
    p3a <- patchwork::wrap_plots(
      list(
        patchwork::free(p3a_top, type = "label", side = "l"),
        patchwork::free(p3a_counts, type = "label", side = "l")
      ),
      ncol = 1, heights = c(3.25, 1)
    )
    figure3_prefix <- if (cmr_source_mode == "local_complete_matrix") "figure3_local_cmr" else "reference_figure3_cmr_from_official_s9"
    yur_save_plot(
      p3a,
      if (cmr_source_mode == "local_complete_matrix") "figure3a_local_cmr" else "reference_figure3a_cmr_from_official_s9",
      cfg, 13.2, 5.9
    )

    top15 <- s9[order(outcome_idx, p), head(.SD, 15), by = Outcome]
    top_protein_order <- unique(top15[order(outcome_idx, p), Protein])
    heat <- s9[Protein %in% top_protein_order]
    heat[, protein_factor := factor(Protein, levels = top_protein_order)]
    heat_layout <- data.table(
      Outcome = cmr_outcome_order,
      anatomy = unname(cmr_anatomy[cmr_outcome_order])
    )
    heat_layout[, group_start := anatomy != shift(anatomy, fill = first(anatomy))]
    heat_layout[, display_idx := seq_len(.N) + cumsum(group_start) * .72]
    heat_layout[, heat_y := max(display_idx) + 1 - display_idx]
    heat <- merge(heat, heat_layout[, .(Outcome, anatomy, heat_y)], by = "Outcome", all.x = TRUE)
    heat[, significance_label := fcase(
      p < .001 / cmr_protein_n, "***",
      p < .01 / cmr_protein_n, "**",
      p < .05 / cmr_protein_n, "*",
      default = ""
    )]
    heat[, label_color := fifelse(abs(beta) >= .30, "white", "#333333")]
    n_heat_proteins <- length(top_protein_order)
    group_boxes <- heat_layout[, .(
      xmin = .5,
      xmax = n_heat_proteins + .5,
      ymin = min(heat_y) - .48,
      ymax = max(heat_y) + .48
    ), by = anatomy]
    group_boxes[, anatomy := factor(anatomy, levels = names(anatomy_colors))]
    heat_anatomy_legend <- data.table(
      anatomy = factor(names(anatomy_colors), levels = names(anatomy_colors)),
      protein_factor = factor(top_protein_order[[1]], levels = top_protein_order),
      heat_y = -Inf
    )
    heat_y_breaks <- heat_layout$heat_y
    heat_y_labels <- cmr_outcome_labels[heat_layout$Outcome]
    p3b <- ggplot2::ggplot(heat, ggplot2::aes(protein_factor, heat_y, fill = beta)) +
      ggplot2::geom_tile(height = .92, color = "white", linewidth = .08) +
      ggplot2::geom_text(
        ggplot2::aes(label = significance_label, color = label_color),
        size = 1.15, fontface = "bold", show.legend = FALSE
      ) +
      ggplot2::geom_rect(
        data = group_boxes,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = anatomy),
        inherit.aes = FALSE, fill = NA, linewidth = .65, show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = heat_anatomy_legend,
        ggplot2::aes(protein_factor, heat_y, color = anatomy),
        inherit.aes = FALSE, shape = 15, size = 0, alpha = 0, na.rm = TRUE
      ) +
      ggplot2::scale_fill_gradient2(
        low = "#4F7FA5", mid = "#FAF7F1", high = "#D87B4D",
        midpoint = 0, limits = c(-.50, .50),
        breaks = c(-.50, -.30, -.10, .10, .30, .50), oob = scales::squish,
        name = "Effect size"
      ) +
      ggplot2::scale_color_manual(
        values = c(anatomy_colors, white = "white", `#333333` = "#333333"),
        breaks = names(anatomy_colors), name = "Anatomy"
      ) +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(
        breaks = heat_y_breaks, labels = heat_y_labels,
        limits = range(heat_layout$heat_y) + c(-.52, .52), expand = c(0, 0)
      ) +
      ggplot2::labs(x = NULL, y = NULL, title = "B") +
      ggplot2::guides(
        fill = ggplot2::guide_colorbar(
          order = 1, title.position = "top", title.hjust = .5,
          barwidth = grid::unit(33, "mm"), barheight = grid::unit(4.2, "mm"),
          ticks = FALSE, frame.colour = NA
        ),
        color = ggplot2::guide_legend(
          order = 2, nrow = 1, title.position = "top", title.hjust = .5,
          override.aes = list(shape = 15, size = 4.8, alpha = 1, linewidth = 0),
          theme = ggplot2::theme(legend.margin = ggplot2::margin(l = 4.5))
        )
      ) +
      ggplot2::theme_minimal(base_size = 8.2, base_family = "Arial") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 11, hjust = -.035),
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(fill = NA, color = "#222222", linewidth = .45),
        axis.text.x = ggplot2::element_text(angle = 58, hjust = 1, vjust = 1, face = "bold", size = 5.15),
        axis.text.y = ggplot2::element_text(face = "bold", size = 7.2),
        legend.position = "top",
        legend.box = "horizontal",
        legend.title = ggplot2::element_text(face = "bold", size = 8.5),
        legend.text = ggplot2::element_text(size = 7.6, face = "bold"),
        legend.spacing.x = grid::unit(2.4, "mm"),
        plot.margin = ggplot2::margin(.6, .8, .5, .2)
      )
    yur_save_plot(
      p3b,
      if (cmr_source_mode == "local_complete_matrix") "figure3b_local_cmr_top15_heatmap" else "reference_figure3b_cmr_top15_heatmap_from_official_s9",
      cfg, 13.2, 6.2
    )

    p3_complete <- patchwork::wrap_plots(
      list(
        patchwork::free(p3a_top, type = "label", side = "l"),
        patchwork::free(p3a_counts, type = "label", side = "l"),
        p3b
      ),
      ncol = 1,
      heights = c(.98 * 3.25 / 4.25, .98 / 4.25, 1.12)
    )
    yur_save_plot(p3_complete, figure3_prefix, cfg, 13.2, 12.0)
    yur_write_csv(cmr_counts, file.path(cfg$paths$figures, paste0(figure3_prefix, "_source_counts.csv")))
    fwrite(s9, file.path(cfg$paths$figures, paste0(figure3_prefix, "_source_associations.csv.gz")), na = "")
    yur_write_csv(
      heat[, .(Outcome, Protein, beta, p, significance_label)],
      file.path(cfg$paths$figures, paste0(figure3_prefix, "_source_heatmap.csv"))
    )
    yur_write_json(list(
      source_mode = cmr_source_mode, participants = if ("n" %in% names(s9)) max(s9$n, na.rm = TRUE) else NA_integer_,
      metrics = uniqueN(s9$Outcome), proteins = cmr_protein_n, rows = nrow(s9),
      bonferroni_threshold = cmr_threshold,
      source_file = if (cmr_source_mode == "local_complete_matrix") local_cmr_file else "official supplementary Table S9"
    ), file.path(cfg$paths$figures, paste0(figure3_prefix, "_manifest.json")))
    if (cmr_source_mode == "local_complete_matrix") {
      stale_reference <- list.files(
        cfg$paths$figures, pattern = "^reference_figure3", full.names = TRUE
      )
      if (length(stale_reference)) unlink(stale_reference, force = TRUE)
    }
  }

  local_mr_file <- file.path(cfg$paths$mr, "mr_results.csv")
  cmest_mediation_file <- file.path(cfg$paths$mediation, "mediation_cmest_results.csv")
  cmest_qc_file <- file.path(cfg$paths$mediation, "mediation_cmest_article_qc.json")
  delta_mediation_file <- file.path(cfg$paths$mediation, "mediation_results.csv")
  local_mediation_file <- if (file.exists(cmest_mediation_file)) {
    cmest_mediation_file
  } else {
    delta_mediation_file
  }
  mediation_source_mode <- if (identical(local_mediation_file, cmest_mediation_file)) {
    "CMAVERSE_BOOTSTRAP_MAIN"
  } else {
    "DELTA_METHOD_AUDIT"
  }
  s13 <- NULL
  s16 <- NULL
  if (file.exists(local_mr_file)) {
    mr_local <- fread(local_mr_file)
    mr_local <- mr_local[primary_method %in% TRUE]
    mr_name_map <- c(
      abdominal_aneurysm = "Abdominal_aneurysm", aortic_valve_stenosis = "Aortic_valve_stenosis",
      atrial_fibrillation = "AF", cardiomyopathy = "Cardiomyopathy", cad = "CAD",
      deep_vein_thrombosis = "Deep_vein_thrombosis", heart_failure = "HF",
      intracerebral_hemorrhage = "Intracerebral_hemorrhage", ischemic_stroke = "Ischemic_stroke",
      peripheral_arterial_disease = "Peripheral_arterial_disease", pulmonary_embolism = "Pulmonary_embolism",
      thoracic_aneurysm = "Thoracic_aneurysm", transient_ischemic_attack = "Transient_ischemic_attack"
    )
    s13 <- mr_local[, .(
      Outcome = unname(mr_name_map[outcome_id]), Exposure = protein,
      Method = method, `P value` = p_value, FDR = fdr, OR = or,
      `95% CI low` = ci_low, `95% CI high` = ci_high, Nsnp = nsnp
    )]
  }
  if (file.exists(local_mediation_file)) {
    med_local <- fread(local_mediation_file)
    if (mediation_source_mode == "CMAVERSE_BOOTSTRAP_MAIN") {
      med_local <- med_local[
        status == "PASS" & significant_bonferroni %in% TRUE &
          is.finite(proportion_mediated_pct) & proportion_mediated_pct >= 0 &
          proportion_mediated_pct <= 100
      ]
      med_local[, `:=`(p_value = indirect_p, cox_warning = "")]
    } else {
      if (!"cox_warning" %in% names(med_local)) med_local[, cox_warning := ""]
      med_local <- med_local[
        !nzchar(cox_warning) & p_value < .05 & is.finite(proportion_mediated_pct) &
          proportion_mediated_pct >= 0 & proportion_mediated_pct <= 100
      ]
    }
    med_name_map <- c(
      abdominal_aneurysm = "Abdominal aneurysm", aortic_valve_stenosis = "Aortic valve stenosis",
      atrial_fibrillation = "AF", cad = "Coronary artery disease",
      deep_vein_thrombosis = "Deep vein thrombosis", heart_failure = "HF",
      ischemic_stroke = "Ischemic stroke", peripheral_arterial_disease = "Peripheral arterial disease",
      pulmonary_embolism = "Pulmonary embolism"
    )
    s16 <- med_local[outcome_id %chin% names(med_name_map), .(
      Exposure = exposure_label, Mediator = protein,
      Outcome = unname(med_name_map[outcome_id]),
      `Proportion mediated   (%)` = proportion_mediated_pct,
      `P value` = p_value, FDR = fdr
    )]
  }
  figure5_parts <- list()
  mr_outcome_map <- c(
    Abdominal_aneurysm = "abdominal_aneurysm",
    Aortic_valve_stenosis = "aortic_valve_stenosis",
    AF = "atrial_fibrillation",
    Cardiomyopathy = "cardiomyopathy",
    CAD = "cad",
    Deep_vein_thrombosis = "deep_vein_thrombosis",
    HF = "heart_failure",
    Intracerebral_hemorrhage = "intracerebral_hemorrhage",
    Ischemic_stroke = "ischemic_stroke",
    Peripheral_arterial_disease = "peripheral_arterial_disease",
    Pulmonary_embolism = "pulmonary_embolism",
    Thoracic_aneurysm = "thoracic_aneurysm",
    Transient_ischemic_attack = "transient_ischemic_attack"
  )
  mr_outcome_order <- article_outcome_order
  if (!is.null(s13)) {
    method_col <- intersect(c("Method", "method"), names(s13))[[1]]
    p_col <- intersect(c("P value", "P.value"), names(s13))[[1]]
    mr <- copy(s13[grepl("inverse variance|ivw", get(method_col), ignore.case = TRUE)])
    mr[, `:=`(
      p = as.numeric(get(p_col)),
      outcome_id = unname(mr_outcome_map[Outcome])
    )]
    mr <- mr[is.finite(p) & p < .05 & !is.na(outcome_id)]
    mr[, outcome_id := factor(outcome_id, levels = mr_outcome_order)]
    mr_present_order <- mr_outcome_order[mr_outcome_order %chin% unique(as.character(mr$outcome_id))]
    mr[, outcome_index := match(as.character(outcome_id), mr_present_order)]
    mr[, within_index := frank(Exposure, ties.method = "first"), by = outcome_id]
    mr[, x := outcome_index + if (.N == 1L) 0 else -.43 + .86 * (within_index - 1) / (.N - 1), by = outcome_id]
    mr[, neglog10p := -log10(pmax(p, .Machine$double.xmin))]
    mr_labels <- mr[order(p), head(.SD, 16)]

    p5a <- ggplot2::ggplot(mr, ggplot2::aes(x, neglog10p, color = outcome_id)) +
      ggplot2::geom_hline(
        yintercept = -log10(.05), linetype = 2,
        linewidth = .42, color = "#C7433E"
      ) +
      ggplot2::geom_point(size = 2.55, alpha = .82, stroke = 0) +
      ggplot2::scale_color_manual(
        values = article_outcome_colors,
        breaks = mr_outcome_order,
        labels = unname(article_outcome_labels[mr_outcome_order]),
        name = "Incident CVD outcomes",
        drop = FALSE
      ) +
      ggplot2::scale_x_continuous(
        breaks = seq_along(mr_present_order), labels = NULL,
        expand = ggplot2::expansion(mult = c(.015, .015))
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, max(mr$neglog10p) * 1.10),
        breaks = seq(0, ceiling(max(mr$neglog10p) / 5) * 5, by = 5),
        expand = ggplot2::expansion(mult = c(0, .01))
      ) +
      ggplot2::labs(x = NULL, y = expression(-log[10](P)), tag = "A") +
      ggplot2::theme_classic(base_family = "Arial", base_size = 9.5) +
      ggplot2::theme(
        axis.title.y = ggplot2::element_text(face = "bold", size = 10),
        axis.text.y = ggplot2::element_text(face = "bold", color = "#222222"),
        axis.line = ggplot2::element_line(linewidth = .35, color = "#222222"),
        panel.border = ggplot2::element_rect(fill = NA, linewidth = .35, color = "#222222"),
        aspect.ratio = 1,
        legend.position = "none",
        plot.tag = ggplot2::element_text(family = "Arial", face = "bold", size = 12),
        plot.tag.position = c(.01, .99),
        plot.margin = ggplot2::margin(2, 3, 1, 2)
      )
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p5a <- p5a + ggrepel::geom_text_repel(
        data = mr_labels, ggplot2::aes(label = Exposure),
        color = "#222222", family = "Arial", fontface = "bold", size = 2.45,
        seed = 20260716, box.padding = .22, point.padding = .12,
        min.segment.length = 0, max.overlaps = Inf, segment.size = .22,
        show.legend = FALSE
      )
    }
    yur_save_plot(p5a, "figure5a_local_mr", cfg, 7.0, 4.3)
    yur_write_csv(mr, file.path(cfg$paths$figures, "figure5a_local_mr_source_data.csv"))
    legend5_data <- data.table(
      outcome_id = factor(mr_outcome_order, levels = mr_outcome_order),
      label = unname(article_outcome_labels[mr_outcome_order]),
      legend_row = c(rep(4, 4), rep(3, 4), rep(2, 4), rep(1, 2)),
      legend_col = c(rep(1:4, 3), 1:2)
    )
    legend5_data[, item_width := .42 + nchar(label) * .057]
    legend5_col_widths <- legend5_data[, .(col_width = max(item_width)), by = legend_col][order(legend_col)]
    legend5_col_widths[, x := cumsum(data.table::shift(col_width, fill = 0))]
    legend5_width <- legend5_col_widths[, max(x + col_width)]
    legend5_data <- merge(legend5_data, legend5_col_widths, by = "legend_col", sort = FALSE)
    # Keep the heading fixed while moving the legend keys down to create a
    # clearer title-to-legend gap.
    legend5_data[, y := legend_row - .90]
    p5legend <- ggplot2::ggplot(legend5_data, ggplot2::aes(x, y, color = outcome_id)) +
      ggplot2::annotate(
        "text", x = legend5_width / 2, y = 4.72, label = "Incident CVD outcomes",
        family = "Arial", fontface = "bold", size = 3.0
      ) +
      ggplot2::geom_point(shape = 15, size = 3.5) +
      ggplot2::geom_text(
        ggplot2::aes(x = x + .11, label = label),
        hjust = 0, color = "#222222", family = "Arial", fontface = "bold", size = 2.5
      ) +
      ggplot2::scale_color_manual(values = article_outcome_colors, guide = "none") +
      ggplot2::coord_cartesian(
        xlim = c(-.05, legend5_width + .05), ylim = c(-.05, 5.38), clip = "off"
      ) +
      ggplot2::theme_void(base_family = "Arial") +
      ggplot2::theme(plot.margin = ggplot2::margin(0, 2, 0, 2))
    figure5_parts$legend <- p5legend
    figure5_parts$a <- p5a
  }
  if (!is.null(s16)) {
    prop_col <- intersect(c("Proportion mediated   (%)", "Proportion mediated (%)"), names(s16))
    if (length(prop_col)) {
      med <- copy(s16)
      med[, source_row := .I]
      med[, proportion := as.numeric(sub(" .*", "", get(prop_col[[1]])))]
      med <- med[is.finite(proportion) & proportion >= 0]
      med_summary <- med[, .(
        n_outcomes = uniqueN(Outcome),
        n_exposures = uniqueN(Exposure),
        median_proportion = median(proportion, na.rm = TRUE)
      ), by = Mediator]
      # Label the two strongest mediators within every outcome-breadth stratum.
      # This keeps annotation balanced across x=1..9 instead of concentrating
      # labels among proteins associated with the largest number of outcomes.
      med_summary[, label_rank_within_breadth := data.table::frank(
        -median_proportion, ties.method = "first"
      ), by = n_outcomes]
      med_summary[, label := fifelse(
        label_rank_within_breadth <= 2L, Mediator, NA_character_
      )]
      med_summary[, label_rule := "top 2 median mediation proportions within each outcome-breadth stratum"]
      exposure_colors <- c(`1` = "#6D7080", `2` = "#8B7B87", `3` = "#A87B86", `4` = "#C37780", `5` = "#DA706F")
      p5b <- ggplot2::ggplot(
        med_summary,
        ggplot2::aes(n_outcomes, median_proportion, size = factor(n_exposures), color = factor(n_exposures))
      ) +
        ggplot2::geom_point(alpha = .72) +
        ggplot2::scale_size_manual(values = c(`1` = 1.8, `2` = 2.3, `3` = 2.8, `4` = 3.3, `5` = 3.8), name = "No. of exposures") +
        ggplot2::scale_color_manual(values = exposure_colors, name = "No. of exposures") +
        ggplot2::scale_x_continuous(breaks = 1:9, limits = c(.75, 9.25), expand = c(0, 0)) +
        ggplot2::scale_y_continuous(
          limits = c(-3, max(52, max(med_summary$median_proportion) * 1.08)),
          breaks = seq(0, 50, by = 10), labels = function(x) paste0(x, "%"), expand = c(0, 0)
        ) +
        ggplot2::labs(x = "No. of associated outcomes", y = "Median proportion mediated by proteins (%)", tag = "B") +
        ggplot2::theme_classic(base_family = "Arial", base_size = 9.5) +
        ggplot2::theme(
          axis.title = ggplot2::element_text(face = "bold", size = 9.2),
          axis.text = ggplot2::element_text(face = "bold", color = "#222222"),
          panel.border = ggplot2::element_rect(fill = NA, linewidth = .35, color = "#222222"),
          aspect.ratio = 1,
          legend.position = "top",
          legend.direction = "horizontal",
          legend.title = ggplot2::element_text(face = "bold", size = 7.4),
          legend.text = ggplot2::element_text(face = "bold", size = 7),
          legend.key.width = grid::unit(3.2, "mm"),
          plot.tag = ggplot2::element_text(family = "Arial", face = "bold", size = 12),
          plot.tag.position = c(.01, .96),
          plot.margin = ggplot2::margin(1, 3, 2, 2)
        ) +
        ggplot2::guides(
          color = ggplot2::guide_legend(order = 1, nrow = 1),
          size = ggplot2::guide_legend(order = 1, nrow = 1)
        )
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p5b <- p5b + ggrepel::geom_text_repel(
          data = med_summary[!is.na(label)], ggplot2::aes(label = label),
          color = "#222222", family = "Arial", fontface = "bold", size = 2.45,
          seed = 20260716, box.padding = .18, point.padding = .1,
          min.segment.length = 0, max.overlaps = Inf, segment.size = .2,
          show.legend = FALSE
        )
      }
      yur_save_plot(p5b, "figure5b_local_mediation_breadth", cfg, 7.0, 4.1)

      med_outcome_map <- c(
        `Abdominal aneurysm` = "abdominal_aneurysm",
        `Aortic valve stenosis` = "aortic_valve_stenosis",
        AF = "atrial_fibrillation",
        `Coronary artery disease` = "cad",
        `Deep vein thrombosis` = "deep_vein_thrombosis",
        HF = "heart_failure",
        `Ischemic stroke` = "ischemic_stroke",
        `Peripheral arterial disease` = "peripheral_arterial_disease",
        `Pulmonary embolism` = "pulmonary_embolism"
      )
      med[, outcome_id := unname(med_outcome_map[Outcome])]
      med <- med[!is.na(outcome_id)]
      med[, outcome_id := factor(outcome_id, levels = mr_outcome_order)]
      # Select the five largest local mediation proportions per outcome, then
      # retain source-table order within each outcome for a stable circular plot.
      top <- med[order(-proportion), head(.SD, 5), by = outcome_id]
      top[, outcome_index := as.integer(outcome_id)]
      setorder(top, outcome_index, source_row)
      top[, group_index := match(as.character(outcome_id), unique(as.character(outcome_id)))]
      gap_width <- 1.2
      outer_gap_width <- 2.0
      top[, bar_index := seq_len(.N) + (group_index - 1L) * gap_width + outer_gap_width / 2]
      top[, exposure_short := fifelse(
        Exposure == "Body mass index", "BMI",
        fifelse(Exposure == "Smoking status", "Smoking",
                fifelse(Exposure == "Systolic blood pressure", "SBP",
                        fifelse(Exposure == "Triglycerides", "TG", Exposure)))
      )]
      top[, phenotype_label := Exposure]
      total_bars <- max(top$bar_index) + outer_gap_width / 2
      top[, angle := 90 - 360 * (bar_index - .5) / total_bars]
      top[, hjust := fifelse(angle < -90, 1, 0)]
      top[angle < -90, angle := angle + 180]
      group_labels <- top[, .(
        bar_index = mean(range(bar_index)),
        x_start = min(bar_index) - .48,
        x_end = max(bar_index) + .48,
        outcome_label = c(
          abdominal_aneurysm = "AAA", aortic_valve_stenosis = "AS",
          atrial_fibrillation = "AF", cad = "CAD", deep_vein_thrombosis = "DVT",
          heart_failure = "HF", ischemic_stroke = "IS",
          peripheral_arterial_disease = "PAD", pulmonary_embolism = "PE"
        )[as.character(outcome_id)][1]
      ), by = outcome_id]
      group_labels[, angle := 90 - 360 * (bar_index - .5) / total_bars]
      group_labels[, hjust := .5]
      group_labels[angle < -90, angle := angle + 180]
      radial_axis_max <- max(80, ceiling(max(top$proportion, na.rm = TRUE) / 20) * 20)
      radial_breaks <- seq(0, radial_axis_max, by = 20)
      radial_scale <- 42 / radial_axis_max
      radial_base <- 29
      protein_radius <- 21.5
      category_arc_inner <- 12.4
      category_arc_outer <- 13.2
      radial_max <- radial_base + radial_axis_max * radial_scale
      separator_arcs <- group_labels[seq_len(.N - 1L), .(gap_x = x_end + gap_width / 2)]
      separator_arcs <- separator_arcs[, .(
        x = gap_x - .22, xend = gap_x + .22,
        y = radial_base + radial_breaks * radial_scale,
        yend = radial_base + radial_breaks * radial_scale
      ), by = gap_x]
      radial_ticks <- data.table(
        x = 0,
        y = radial_base + radial_breaks * radial_scale,
        label = as.character(radial_breaks)
      )

      p5c <- ggplot2::ggplot() +
        ggplot2::geom_rect(
          data = group_labels,
          ggplot2::aes(xmin = x_start, xmax = x_end, ymin = category_arc_inner, ymax = category_arc_outer, fill = outcome_id),
          color = "white", linewidth = .25, show.legend = FALSE
        ) +
        ggplot2::geom_rect(
          data = top,
          ggplot2::aes(
            xmin = bar_index - .43, xmax = bar_index + .43,
            ymin = radial_base, ymax = radial_base + proportion * radial_scale,
            fill = outcome_id
          ),
          color = "white", linewidth = .24, show.legend = FALSE
        ) +
        ggplot2::geom_segment(
          data = separator_arcs,
          ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
          inherit.aes = FALSE, color = "#89979A", linewidth = .72, lineend = "round"
        ) +
        ggplot2::geom_text(
          data = radial_ticks,
          ggplot2::aes(x = x, y = y, label = label),
          inherit.aes = FALSE, size = 1.62, color = "#333333",
          hjust = .5, angle = 0
        ) +
        ggplot2::geom_text(
          data = group_labels,
          ggplot2::aes(x = bar_index, y = 15.5, label = outcome_label, angle = angle),
          inherit.aes = FALSE, size = 1.85, color = "#333333", fontface = "bold"
        ) +
        ggplot2::geom_text(
          data = top,
          ggplot2::aes(x = bar_index, y = protein_radius, label = Mediator, angle = angle),
          inherit.aes = FALSE, size = 1.62, color = "#333333"
        ) +
        ggplot2::geom_text(
          data = top,
          ggplot2::aes(
            x = bar_index,
            y = radial_base + proportion * radial_scale * .54,
            label = sprintf("%.1f", proportion), angle = angle
          ),
          inherit.aes = FALSE, size = 1.72, color = "#111111", fontface = "bold"
        ) +
        ggplot2::geom_text(
          data = top,
          ggplot2::aes(
            x = bar_index, y = radial_base + proportion * radial_scale + 2.2,
            label = phenotype_label, angle = angle, hjust = hjust, color = outcome_id
          ),
          inherit.aes = FALSE, size = 1.72, show.legend = FALSE
        ) +
        ggplot2::annotate(
          "rect", xmin = 0, xmax = total_bars, ymin = 0, ymax = 8.5,
          fill = "white", color = NA
        ) +
        ggplot2::annotate(
          "text", x = total_bars / 2, y = 0, label = "Mediation\nproportion (%)",
          size = 2.35, color = "#222222", hjust = .5, vjust = .5, lineheight = .92
        ) +
        ggplot2::scale_fill_manual(values = article_outcome_colors, drop = FALSE) +
        ggplot2::scale_color_manual(values = article_outcome_colors, drop = FALSE) +
        ggplot2::scale_x_continuous(limits = c(0, total_bars), expand = c(0, 0)) +
        ggplot2::scale_y_continuous(limits = c(0, radial_max + 5.5), expand = c(0, 0)) +
        ggplot2::coord_polar(start = 0, clip = "off") +
        ggplot2::labs(tag = "C") +
        ggplot2::theme_void(base_family = "Arial") +
        ggplot2::theme(
          plot.tag = ggplot2::element_text(family = "Arial", face = "bold", size = 12),
          plot.tag.position = c(.105, .90),
          plot.margin = ggplot2::margin(1.5, 2.5, 1.5, 2.5)
        )
      yur_save_plot(p5c, "figure5c_local_mediation", cfg, 8.2, 8.2)
      yur_write_csv(med_summary, file.path(cfg$paths$figures, "figure5b_local_mediation_source_data.csv"))
      yur_write_csv(top, file.path(cfg$paths$figures, "figure5c_local_mediation_source_data.csv"))
      figure5_parts$b <- p5b
      figure5_parts$c <- p5c
    }
  }
  if (all(c("legend", "a", "b", "c") %in% names(figure5_parts))) {
    if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required for article-style Figure 5.")
    p5_left <- patchwork::wrap_plots(
      list(patchwork::plot_spacer(), figure5_parts$a, figure5_parts$b),
      ncol = 1, heights = c(.06, .47, .47)
    )
    p5_right <- patchwork::wrap_plots(
      list(patchwork::plot_spacer(), figure5_parts$legend, figure5_parts$c),
      ncol = 1, heights = c(.06, .09, .85)
    )
    p5_complete <- patchwork::wrap_plots(
      list(p5_left, p5_right), ncol = 2, widths = c(.68, 1.32)
    )
    yur_save_plot(p5_complete, "figure5_local_mr_and_mediation", cfg, 12.2, 8.8)
    yur_write_json(list(
      source_mode = "LOCAL_UKB_RECOMPUTATION",
      mr_file = local_mr_file, mr_sha256 = yur_sha_file(local_mr_file),
      mediation_file = local_mediation_file, mediation_sha256 = yur_sha_file(local_mediation_file),
      mediation_source_mode = mediation_source_mode,
      mediation_figure_inclusion_rule = if (mediation_source_mode == "CMAVERSE_BOOTSTRAP_MAIN") {
        "CMAverse indirect-effect bootstrap P<0.05/local candidate count"
      } else {
        "delta-method nominal mediation P<0.05; audit figure only"
      },
      mediation_figure_proportion_range = "0 to 100%; suppression paths outside this range remain in mediation_results.csv",
      mediation_inference = if (mediation_source_mode == "CMAVERSE_BOOTSTRAP_MAIN") {
        "CMAverse cmest regression-based paramfunc; linear mediator; Cox outcome; percentile bootstrap"
      } else {
        "product-of-coefficients delta method"
      },
      exact_article_cma_method = FALSE,
      exact_article_cma_status = if (mediation_source_mode == "CMAVERSE_BOOTSTRAP_MAIN") {
        "MATCHED_REPORTED_COMPONENTS_WITH_RECORDED_UNREPORTED_ASSUMPTIONS"
      } else {
        "CMAverse main track not yet complete"
      },
      mediation_numeric_replication_qc = if (file.exists(cmest_qc_file)) cmest_qc_file else NA_character_,
      mediation_numeric_replication_status = if (file.exists(cmest_qc_file)) {
        as.character(read_json(cmest_qc_file, simplifyVector = TRUE)$status)
      } else {
        "NOT_AVAILABLE"
      },
      official_s13_s16_used = FALSE
    ), file.path(cfg$paths$figures, "figure5_local_manifest.json"))
    stale_reference <- list.files(cfg$paths$figures, pattern = "^reference_figure5", full.names = TRUE)
    if (length(stale_reference)) unlink(stale_reference, force = TRUE)
  }

  s17 <- yur_reference_sheet(cfg, "S17")
  s18 <- yur_reference_sheet(cfg, "S18")
  local_prs_file <- file.path(cfg$paths$prs, "figure6a_source_data.csv")
  prs <- NULL
  prs_source_mode <- NULL
  if (file.exists(local_prs_file)) {
    prs <- fread(local_prs_file)
    required_prs <- c(
      "outcome_id", "outcome_label", "protein", "threshold_label", "beta", "p_value",
      "bonferroni_significant"
    )
    missing_prs <- setdiff(required_prs, names(prs))
    if (length(missing_prs)) {
      stop("Local Figure 6A source is missing: ", paste(missing_prs, collapse = ", "), call. = FALSE)
    }
    prs[, `:=`(
      source_row = .I,
      disease = outcome_label,
      Protein = protein,
      threshold = as.character(threshold_label),
      beta = as.numeric(beta),
      p_value = as.numeric(p_value),
      cell_significant = as.logical(bonferroni_significant)
    )]
    prs_source_mode <- "local_complete_matrix"
  } else if (!is.null(s17)) {
    prs <- copy(s17)
    prs[, source_row := .I]
    prs[, disease := trimws(sub(" \\(Pt.*", "", Outcome))]
    prs[, threshold := sub(".*\\(Pt ([^)]+)\\).*", "\\1", Outcome)]
    prs[, beta := as.numeric(get(intersect(c("beta", "Beta", "β"), names(prs))[[1]]))]
    prs[, p_value := as.numeric(get(intersect(c("P value", "P.value"), names(prs))[[1]]))]
    prs_disease_map <- c(
      `Abdominal aneurysm` = "abdominal_aneurysm",
      `Aortic valve stenosis` = "aortic_valve_stenosis",
      `Atrial fibrillation` = "atrial_fibrillation",
      `Coronary artery disease` = "cad",
      `Deep vein thrombosis` = "deep_vein_thrombosis",
      `Heart failure` = "heart_failure",
      `Ischemic stroke` = "ischemic_stroke",
      `Peripheral arterial disease` = "peripheral_arterial_disease",
      `Pulmonary embolism` = "pulmonary_embolism"
    )
    prs[, outcome_id := unname(prs_disease_map[disease])]
    prs[, cell_significant := TRUE]
    prs_source_mode <- "official_significant_rows_only"
  }
  if (!is.null(prs)) {
    prs <- prs[!is.na(outcome_id) & is.finite(beta)]
    threshold_contract <- fread(file.path(cfg$project_dir, "f", "config", "prs_thresholds.tsv"))
    threshold_levels <- as.character(threshold_contract$label)
    threshold_labels <- paste("Pt", threshold_levels)
    prs[, threshold := factor(threshold, levels = threshold_levels)]
    prs[, pair_id := paste(disease, Protein, sep = "||")]
    coef_cap <- max(1e-4, as.numeric(stats::quantile(abs(prs$beta), .98, na.rm = TRUE)))
    disease_alpha <- sort(unique(prs$disease))

    make_prs_block <- function(block_pair_ids, y_top, x_label, x_disease, x_cells) {
      pairs <- prs[pair_id %chin% block_pair_ids, .(
        source_row = min(source_row),
        outcome_id = outcome_id[[1]],
        Protein = Protein[[1]]
      ), by = .(pair_id, disease)]
      pairs[, disease_rank := match(disease, disease_alpha)]
      setorder(pairs, disease_rank, source_row)
      pairs[, y := y_top - seq_len(.N) + 1]
      pairs[, `:=`(x_disease = x_disease, x_label = x_label)]
      grid <- merge(
        data.table::CJ(pair_id = pairs$pair_id, threshold = threshold_levels, unique = TRUE),
        pairs[, .(pair_id, disease, outcome_id, Protein, source_row, y)],
        by = "pair_id", all.x = TRUE, sort = FALSE
      )
      grid <- merge(
        grid,
        prs[, .(pair_id, threshold = as.character(threshold), beta, p_value, cell_significant)],
        by = c("pair_id", "threshold"), all.x = TRUE, sort = FALSE
      )
      grid[, threshold_index := match(threshold, threshold_levels)]
      grid[, x := x_cells + threshold_index - 1]
      list(pairs = pairs, grid = grid)
    }

    disease_counts <- prs[, .(pair_n = uniqueN(pair_id)), by = disease]
    disease_counts[, disease_rank := match(disease, disease_alpha)]
    setorder(disease_counts, disease_rank)
    split_candidates <- data.table(split_at = seq_len(nrow(disease_counts) - 1L))
    split_candidates[, left_n := vapply(split_at, function(k) sum(disease_counts$pair_n[seq_len(k)]), numeric(1))]
    split_candidates[, right_n := sum(disease_counts$pair_n) - left_n]
    split_candidates[, legend_gap := right_n - left_n]
    legend_rows_required <- length(disease_alpha) * 2L + 8L
    valid_splits <- split_candidates[right_n > left_n & legend_gap >= legend_rows_required]
    if (nrow(valid_splits)) {
      split_at <- valid_splits[which.min(legend_gap), split_at]
    } else {
      split_at <- split_candidates[which.min(abs(legend_gap)), split_at]
    }
    left_diseases <- disease_counts$disease[seq_len(split_at)]
    right_diseases <- disease_counts$disease[seq.int(split_at + 1L, nrow(disease_counts))]
    left_pair_ids <- unique(prs[disease %chin% left_diseases, pair_id])
    right_pair_ids <- unique(prs[disease %chin% right_diseases, pair_id])
    left_n <- length(left_pair_ids)
    right_n <- length(right_pair_ids)
    plot_top <- max(right_n, left_n + legend_rows_required, 30L)
    heatmap_cell_width <- .78
    heatmap_cell_gap <- 1 - heatmap_cell_width
    disease_marker_width <- .16
    disease_to_cell_offset <- heatmap_cell_width / 2 + heatmap_cell_gap + disease_marker_width / 2
    left_cell_start <- 3.35
    right_cell_start <- 12.50
    left_disease_x <- left_cell_start - disease_to_cell_offset
    right_disease_x <- right_cell_start - disease_to_cell_offset
    protein_to_disease_gap <- .24
    left_block <- make_prs_block(
      left_pair_ids, left_n, left_disease_x - protein_to_disease_gap, left_disease_x, left_cell_start
    )
    right_block <- make_prs_block(
      right_pair_ids, right_n, right_disease_x - protein_to_disease_gap, right_disease_x, right_cell_start
    )
    left_block$pairs[, label_hjust := 1]
    right_block$pairs[, `:=`(
      x_label = right_cell_start + 4 + heatmap_cell_width / 2 + protein_to_disease_gap,
      label_hjust = 0
    )]
    pairs_plot <- data.table::rbindlist(list(left_block$pairs, right_block$pairs), fill = TRUE)
    grid_plot <- data.table::rbindlist(list(left_block$grid, right_block$grid), fill = TRUE)

    coefficient_legend <- data.table(
      x = .35,
      y = seq(plot_top - 9.4, plot_top - 2.90, length.out = 101),
      beta = seq(-coef_cap, coef_cap, length.out = 101)
    )
    coefficient_ticks <- data.table(
      x = .76,
      y = c(plot_top - 9.4, plot_top - 6.15, plot_top - 2.90),
      label = formatC(c(-coef_cap, 0, coef_cap), format = "e", digits = 1)
    )
    legend_diseases <- disease_alpha
    legend_ids <- vapply(legend_diseases, function(label) prs[disease == label, outcome_id][[1]], character(1))
    disease_legend <- data.table(
      disease = legend_diseases,
      outcome_id = legend_ids,
      x = 0,
      y = plot_top - 13.5 - (seq_along(legend_ids) - 1L) * 1.40
    )
    disease_legend[, label := unname(article_outcome_labels[outcome_id])]

    threshold_axis <- data.table::rbindlist(list(
      data.table(x = c(left_disease_x, left_cell_start + 0:4), label = c("CVDs", threshold_labels)),
      data.table(x = c(right_disease_x, right_cell_start + 0:4), label = c("CVDs", threshold_labels))
    ))

    p6a <- ggplot2::ggplot() +
      ggplot2::geom_tile(
        data = coefficient_legend,
        ggplot2::aes(x = x, y = y, fill = beta),
        width = .46, height = .095, linewidth = 0
      ) +
      ggplot2::geom_text(
        data = coefficient_ticks,
        ggplot2::aes(x = x, y = y, label = label),
        hjust = 0, family = "sans", fontface = "bold", size = 2.30, color = "#252525"
      ) +
      ggplot2::annotate(
        "text", x = 1.05, y = plot_top + .30, label = "Standardized\ncoefficient",
        family = "sans", fontface = "bold", size = 2.95, lineheight = .9,
        hjust = 0, vjust = 1
      ) +
      ggplot2::annotate(
        "text", x = 0, y = plot_top - 11.8, label = "Incident CVD outcomes",
        family = "sans", fontface = "bold", hjust = 0, size = 2.65
      ) +
      ggplot2::geom_segment(
        data = disease_legend,
        ggplot2::aes(x = x, xend = x, y = y - .53, yend = y + .53, color = outcome_id),
        linewidth = 2.2, lineend = "butt", show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = disease_legend,
        ggplot2::aes(x = x + .42, y = y, label = label),
        family = "sans", fontface = "bold", hjust = 0, size = 2.35, color = "#252525"
      ) +
      ggplot2::geom_segment(
        data = pairs_plot,
        ggplot2::aes(x = x_disease, xend = x_disease, y = y - .39, yend = y + .39, color = outcome_id),
        linewidth = 1.8, lineend = "butt", show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = left_block$pairs,
        ggplot2::aes(x = x_label, y = y, label = Protein),
        hjust = 1, family = "sans", fontface = "bold", size = 2.30, color = "#252525"
      ) +
      ggplot2::geom_text(
        data = right_block$pairs,
        ggplot2::aes(y = y, label = Protein),
        x = right_cell_start + 4 + heatmap_cell_width / 2 + protein_to_disease_gap,
        hjust = 0, family = "sans", fontface = "bold", size = 2.30, color = "#252525"
      ) +
      ggplot2::geom_tile(
        data = grid_plot,
        ggplot2::aes(x = x, y = y, fill = beta),
        width = heatmap_cell_width, height = .78, color = "white", linewidth = .10
      ) +
      ggplot2::geom_text(
        data = grid_plot[cell_significant %in% TRUE],
        ggplot2::aes(x = x, y = y, label = "*"),
        family = "sans", fontface = "bold", size = 1.48, color = "#202020", vjust = .72
      ) +
      ggplot2::geom_text(
        data = threshold_axis,
        ggplot2::aes(x = x, y = -1.0, label = label),
        family = "sans", fontface = "bold", size = 2.30,
        angle = 90, hjust = 1, vjust = .5, color = "#252525"
      ) +
      ggplot2::scale_fill_gradient2(
        low = "#3D74A6", mid = "#F7F7F7", high = "#C94E45", midpoint = 0,
        limits = c(-coef_cap, coef_cap), oob = scales::squish,
        na.value = "white", guide = "none"
      ) +
      ggplot2::scale_color_manual(values = article_outcome_colors, guide = "none") +
      ggplot2::coord_fixed(
        ratio = 1, xlim = c(-.65, 19.0), ylim = c(-5.5, plot_top + 1.5), clip = "off"
      ) +
      ggplot2::labs(tag = "A") +
      ggplot2::theme_void(base_family = "sans") +
      ggplot2::theme(
        plot.tag = ggplot2::element_text(family = "sans", face = "bold", size = 12),
        plot.tag.position = c(.005, .995),
        plot.margin = ggplot2::margin(1, 1, 1, 1)
      )

    figure6a_name <- if (prs_source_mode == "local_complete_matrix") {
      "figure6a_local_prs_protein_heatmap"
    } else {
      "reference_figure6a_prs_from_official_s17"
    }
    yur_save_plot(p6a, figure6a_name, cfg, 4.8, max(8.8, plot_top * .17))
    yur_write_csv(prs, file.path(cfg$paths$figures, paste0(figure6a_name, "_source_data.csv")))
    yur_write_json(list(
      source_mode = prs_source_mode,
      complete_non_significant_cells_available = identical(prs_source_mode, "local_complete_matrix"),
      outcomes = uniqueN(prs$outcome_id), pairs = uniqueN(prs$pair_id), rows = nrow(prs)
    ), file.path(cfg$paths$figures, paste0(figure6a_name, "_manifest.json")))
  }
  if (!is.null(s18)) {
    s18[, qlog := as.numeric(`-Log10(q)`)]
    s18[, count := as.numeric(`Count(bubble size)`)]
    p6b <- ggplot2::ggplot(s18, ggplot2::aes(qlog, reorder(Description, qlog), size = count, color = `Category(bubble color)`)) +
      ggplot2::geom_point(alpha = .85) + ggplot2::labs(x = expression(-log[10](q)), y = NULL, size = "Count", color = "Source", title = "Reference enrichment from Table S18") + yur_theme()
    yur_save_plot(p6b, "reference_figure6b_enrichment_from_official_s18", cfg, 8.5, 6.5)
  }

  gates <- data.table(
    figure = c("Figure 1 local workflow", "Figure 2 local Cox", "Figure 3 local CMR", "Figure 4 local prediction", "Figure 5 local MR/mediation", "Figure 6 local PRS/systems"),
    status = c("CODE_READY", ifelse(file.exists(cox_file), "GENERATED", "WAITING_COX"),
               ifelse(file.exists(file.path(cfg$paths$cmr, "cmr_associations.csv.gz")), "READY", "WAITING_LOCAL_CMR"),
               ifelse(file.exists(metrics_file), "GENERATED", "WAITING_PREDICTION"),
               ifelse(file.exists(file.path(cfg$paths$mr, "mr_results.csv")) &&
                        file.exists(file.path(cfg$paths$mediation, "mediation_results.csv")), "READY", "WAITING_LOCAL_MR_MEDIATION"),
               ifelse(any(file.exists(c(
                 file.path(cfg$paths$prs, "prs_protein_associations.csv"),
                 file.path(cfg$paths$prs, "prs_protein_associations.csv.gz"),
                 file.path(cfg$paths$prs, "figure6a_source_data.csv")
               ))), "READY", "WAITING_EXTERNAL_PRS_NETWORK"))
  )
  yur_write_csv(gates, file.path(cfg$paths$figures, "figure_generation_gates.csv"))
  artifacts <- c(
    list.files(cfg$paths$figures, full.names = TRUE),
    list.files(cfg$paths$supplement_figures, full.names = TRUE)
  )
  artifacts <- artifacts[file.exists(artifacts) & !dir.exists(artifacts)]
  if (length(artifacts)) {
    provenance <- data.table(
      file = basename(artifacts),
      path = normalizePath(artifacts, winslash = "/", mustWork = TRUE),
      bytes = file.info(artifacts)$size,
      sha256 = vapply(artifacts, yur_sha_file, character(1))
    )
    provenance[, evidence_scope := fifelse(
      startsWith(file, "reference_"), "OFFICIAL_SUPPLEMENT_TABLE_RECONSTRUCTION",
      "LOCAL_SOURCE_LOCKED_RECOMPUTATION"
    )]
    yur_write_csv(provenance, file.path(cfg$paths$figures, "figure_provenance_manifest.csv"))
  }
}

yur_build_report <- function(cfg) {
  files <- list.files(cfg$analysis_dir, recursive = TRUE, full.names = TRUE)
  files <- files[file.exists(files) & !dir.exists(files)]
  manifest <- data.table(
    path = normalizePath(files, winslash = "/", mustWork = TRUE),
    bytes = file.info(files)$size,
    sha256 = vapply(files, yur_sha_file, character(1))
  )
  yur_write_csv(manifest, file.path(cfg$paths$report, "result_file_manifest.csv"))
  cox_summary <- file.path(cfg$paths$cox, "cox_summary.json")
  selection_summary <- file.path(cfg$paths$selection, "selection_summary.json")
  evaluation_summary <- file.path(cfg$paths$evaluation, "evaluation_summary.json")
  selection <- if (file.exists(selection_summary)) read_json(selection_summary, simplifyVector = TRUE) else NULL
  score2_qc_file <- file.path(cfg$paths$cohort, "score2_input_output_qc.csv")
  score2_qc_status <- if (file.exists(score2_qc_file) && all(fread(score2_qc_file)$status == "PASS")) "PASS" else "MISSING_OR_FAIL"
  status <- data.table(
    module = c("sources", "cohort", "cox", "selection", "prediction", "figures", "cmr", "mr", "mediation", "prs", "enrichment"),
    status = c(
      ifelse(file.exists(file.path(cfg$paths$sources, "code_availability_audit.json")), "PASS", "MISSING"),
      ifelse(file.exists(file.path(cfg$paths$cohort, "cohort_summary.json")), "PASS", "MISSING"),
      ifelse(file.exists(cox_summary), "PASS", "MISSING"),
      ifelse(file.exists(selection_summary), "PASS", "MISSING"),
      ifelse(file.exists(evaluation_summary), "PASS", "MISSING"),
      ifelse(file.exists(file.path(cfg$paths$figures, "figure_provenance_manifest.csv")), "PASS", "NOT_RUN"),
      ifelse(file.exists(file.path(cfg$paths$cmr, "cmr_associations.csv.gz")), "PASS", "NOT_RUN"),
      ifelse(file.exists(file.path(cfg$paths$mr, "mr_results.csv")), "PASS", "EXTERNAL_INPUT_REQUIRED"),
      ifelse(file.exists(file.path(cfg$paths$mediation, "mediation_results.csv")), "PASS", "NOT_RUN"),
      ifelse(any(file.exists(c(
        file.path(cfg$paths$prs, "prs_protein_associations.csv"),
        file.path(cfg$paths$prs, "prs_protein_associations.csv.gz"),
        file.path(cfg$paths$prs, "figure6a_source_data.csv")
      ))), "PASS", "EXTERNAL_INPUT_REQUIRED"),
      ifelse(file.exists(file.path(cfg$paths$enrichment, "enrichment_results.csv")), "PASS", "EXTERNAL_WEB_EXPORT_REQUIRED")
    )
  )
  yur_write_csv(status, file.path(cfg$paths$report, "module_status.csv"))
  lines <- c(
    "# Yu/Chen 2025 reproduction report", "",
    paste0("Prediction panel mode: `", cfg$prediction_panel_mode, "`."), "",
    "This report distinguishes locally recomputed results from official-table reference reconstructions.", "",
    "## Code availability", "",
    "No author source-code repository was identified. The article lists software but does not release scripts, split EIDs, or tuning code.", "",
    "## Module status", "",
    "|Module|Status|", "|---|---|",
    vapply(seq_len(nrow(status)), function(i) paste0("|", status$module[[i]], "|", status$status[[i]], "|"), character(1)), "",
    "## Interpretation guardrails", "",
    paste0("- SCORE2 source and saturation QC: `", score2_qc_status, "`."),
    if (!is.null(selection)) paste0(
      "- Local derivation reselection yielded ", selection$local_reselected_union_n,
      " proteins and is the primary prediction panel; published S12 contains ", selection$published_union_n,
      " reference proteins. Overlap: ", selection$local_published_overlap_n, "."
    ) else "- Selection summary is not available.",
    "- Only an all-14-endpoint selection run can be compared with the published 671-candidate and 257-protein anchors.",
    "- The local derivation-only panel is the primary prediction panel; the published 257 panel is reference-only.",
    "- Figures whose names start with `reference_` are reconstructions from official supplementary tables, not local participant-level results.",
    "- The reported prediction AUC is an incident-event classification AUC over observed follow-up, not a fixed 1/3/5/10-year AUC."
  )
  writeLines(lines, file.path(cfg$paths$report, "RESULTS_AND_QC.md"))
}
