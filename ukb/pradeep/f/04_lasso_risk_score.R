# Protein, clinical and combined LASSO risk-score construction/testing.

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

message("Step 04: LASSO risk-score construction/testing")
if (!file.exists(analysis_base_file)) {
  stop("Missing analysis base. Run 01_build_analysis_base.R first: ", analysis_base_file, call. = FALSE)
}

base <- readRDS(analysis_base_file)
dt <- as.data.frame(base$dat)
proteins <- base$protein_cols
if (max_proteins > 0L) proteins <- head(proteins, max_proteins)
outcomes <- base$outcome_map[match(analysis_outcome_keys, base$outcome_map$outcome_key)]

lasso_dir <- file.path(output_dir, "lasso")
dir.create(lasso_dir, recursive = TRUE, showWarnings = FALSE)
force_lasso <- repro_bool_env("UKBPPP_FORCE_LASSO", default = FALSE)
lasso_nfolds <- max(3L, repro_int_env("UKBPPP_LASSO_NFOLDS", 10L))
lasso_type_measure <- Sys.getenv("UKBPPP_LASSO_TYPE_MEASURE", unset = "auc")
glmnet_maxit <- max(100000L, repro_int_env("UKBPPP_GLMNET_MAXIT", 1000000L))

impute_for_model <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    x <- suppressWarnings(as.numeric(x))
    med <- stats::median(x, na.rm = TRUE)
    if (!is.finite(med)) med <- 0
    x[is.na(x)] <- med
    return(x)
  }
  x <- as.character(x)
  mode_val <- repro_mode_value(x)
  x[is.na(x)] <- mode_val
  factor(x)
}

make_x <- function(dat, vars) {
  dat <- dat[, vars, drop = FALSE]
  dat[] <- lapply(dat, impute_for_model)
  stats::model.matrix(~ . - 1, data = dat)
}

safe_divide <- function(num, den) {
  ifelse(is.finite(den) & den > 0, num / den, NA_real_)
}

calc_prediction_metrics <- function(pred, y, roc_obj) {
  pred <- pmin(pmax(as.numeric(pred), 1e-6), 1 - 1e-6)
  y <- as.integer(y)
  coords <- tryCatch(
    pROC::coords(roc_obj, "best", best.method = "youden",
                 ret = c("threshold", "sensitivity", "specificity"),
                 transpose = FALSE),
    error = function(e) NULL
  )
  threshold <- if (!is.null(coords) && "threshold" %in% names(coords)) as.numeric(coords$threshold[1]) else NA_real_
  if (!is.finite(threshold)) threshold <- mean(y == 1, na.rm = TRUE)
  class_pred <- as.integer(pred >= threshold)
  tp <- sum(class_pred == 1 & y == 1, na.rm = TRUE)
  tn <- sum(class_pred == 0 & y == 0, na.rm = TRUE)
  fp <- sum(class_pred == 1 & y == 0, na.rm = TRUE)
  fn <- sum(class_pred == 0 & y == 1, na.rm = TRUE)

  logit_pred <- stats::qlogis(pred)
  cal <- tryCatch(stats::glm(y ~ logit_pred, family = stats::binomial()), error = function(e) NULL)
  cal_intercept <- if (!is.null(cal) && length(stats::coef(cal)) >= 1) unname(stats::coef(cal)[1]) else NA_real_
  cal_slope <- if (!is.null(cal) && length(stats::coef(cal)) >= 2) unname(stats::coef(cal)[2]) else NA_real_

  data.frame(
    threshold = threshold,
    brier = mean((pred - y)^2, na.rm = TRUE),
    observed_event_rate = mean(y == 1, na.rm = TRUE),
    mean_predicted_risk = mean(pred, na.rm = TRUE),
    accuracy = safe_divide(tp + tn, tp + tn + fp + fn),
    sensitivity = safe_divide(tp, tp + fn),
    specificity = safe_divide(tn, tn + fp),
    ppv = safe_divide(tp, tp + fp),
    npv = safe_divide(tn, tn + fn),
    f1 = safe_divide(2 * tp, 2 * tp + fp + fn),
    calibration_intercept = cal_intercept,
    calibration_slope = cal_slope,
    tp = tp,
    fp = fp,
    tn = tn,
    fn = fn,
    stringsAsFactors = FALSE
  )
}

make_calibration_data <- function(pred, y, n_bins = 10L) {
  d <- data.table::data.table(pred = as.numeric(pred), y = as.integer(y))
  d <- d[is.finite(pred) & !is.na(y)]
  if (nrow(d) == 0) return(NULL)
  breaks <- unique(stats::quantile(d$pred, probs = seq(0, 1, length.out = n_bins + 1L), na.rm = TRUE, type = 8))
  if (length(breaks) <= 2L) {
    d[, risk_bin := "1"]
  } else {
    d[, risk_bin := as.character(cut(pred, breaks = breaks, include.lowest = TRUE, labels = FALSE))]
  }
  d[, .(
    n = .N,
    observed_rate = mean(y == 1, na.rm = TRUE),
    mean_predicted_risk = mean(pred, na.rm = TRUE),
    min_predicted_risk = min(pred, na.rm = TRUE),
    max_predicted_risk = max(pred, na.rm = TRUE)
  ), by = risk_bin][order(as.integer(risk_bin))]
}

fit_model_auc <- function(x_train, y_train, x_test, y_test, outcome, model_name, lambda_rule,
                          id_test = NULL, fu_test = NULL) {
  if (length(unique(y_train)) < 2 || length(unique(y_test)) < 2 || ncol(x_train) == 0) {
    return(list(summary = data.frame(
      model = model_name, outc = outcome, auc = NA_real_, auc_lci = NA_real_, auc_uci = NA_real_,
      n_train = length(y_train), n_test = length(y_test), n_features = ncol(x_train),
      error = "single outcome class or no predictors", stringsAsFactors = FALSE
    ), roc = NULL, coef = NULL, lambda = NULL, pred = NULL, metrics = NULL, calibration = NULL))
  }
  set.seed(1234)
  cvfit <- tryCatch(glmnet::cv.glmnet(
    x_train, y_train, alpha = 1, family = "binomial", nfolds = lasso_nfolds,
    type.measure = lasso_type_measure, keep = TRUE, maxit = glmnet_maxit
  ), error = function(e) e)
  if (inherits(cvfit, "error")) {
    return(list(summary = data.frame(
      model = model_name, outc = outcome, auc = NA_real_, auc_lci = NA_real_, auc_uci = NA_real_,
      n_train = length(y_train), n_test = length(y_test), n_features = ncol(x_train),
      error = conditionMessage(cvfit), stringsAsFactors = FALSE
    ), roc = NULL, coef = NULL, lambda = NULL, pred = NULL, metrics = NULL, calibration = NULL))
  }
  lambda <- if (lambda_rule == "lambda.min") cvfit$lambda.min else cvfit$lambda.1se
  pred <- as.numeric(stats::predict(cvfit, newx = x_test, s = lambda, type = "response"))
  score_link <- as.numeric(stats::predict(cvfit, newx = x_test, s = lambda, type = "link"))
  score_sd <- as.numeric(scale(score_link))
  if (all(!is.finite(score_sd))) score_sd <- score_link
  roc_obj <- pROC::roc(response = as.numeric(y_test), predictor = pred, levels = c(0, 1), ci = TRUE, quiet = TRUE)
  ci <- as.numeric(roc_obj$ci)
  roc_df <- data.frame(
    FPR = rev(1 - roc_obj$specificities),
    TPR = rev(roc_obj$sensitivities),
    model = model_name,
    outc = outcome,
    stringsAsFactors = FALSE
  )
  coef_mat <- as.matrix(stats::coef(cvfit, s = lambda))
  coef_df <- data.frame(
    feature = rownames(coef_mat),
    coefficient = as.numeric(coef_mat[, 1]),
    model = model_name,
    outc = outcome,
    stringsAsFactors = FALSE
  )
  lambda_df <- data.frame(
    lambda = cvfit$lambda,
    cvm = cvfit$cvm,
    cvsd = cvfit$cvsd,
    cvup = cvfit$cvup,
    cvlo = cvfit$cvlo,
    nzero = cvfit$nzero,
    model = model_name,
    outc = outcome,
    stringsAsFactors = FALSE
  )
  summary_df <- data.frame(
    model = model_name,
    outc = outcome,
    auc = ci[2],
    auc_lci = ci[1],
    auc_uci = ci[3],
    n_train = length(y_train),
    n_test = length(y_test),
    n_features = ncol(x_train),
    error = NA_character_,
    stringsAsFactors = FALSE
  )
  if (is.null(id_test)) id_test <- seq_along(y_test)
  if (is.null(fu_test)) fu_test <- rep(NA_real_, length(y_test))
  pred_df <- data.frame(
    model = model_name,
    outc = outcome,
    eid = id_test,
    followup_years = as.numeric(fu_test),
    y = as.integer(y_test),
    predicted_risk = pred,
    score_link = score_link,
    score_sd = score_sd,
    stringsAsFactors = FALSE
  )
  metrics_df <- calc_prediction_metrics(pred, y_test, roc_obj)
  metrics_df$model <- model_name
  metrics_df$outc <- outcome
  metrics_df$auc <- ci[2]
  metrics_df$auc_lci <- ci[1]
  metrics_df$auc_uci <- ci[3]
  metrics_df$n_test <- length(y_test)
  metrics_df$n_events_test <- sum(y_test == 1)
  metrics_df <- metrics_df[, c(
    "model", "outc", "auc", "auc_lci", "auc_uci", "n_test", "n_events_test",
    setdiff(names(metrics_df), c("model", "outc", "auc", "auc_lci", "auc_uci", "n_test", "n_events_test"))
  )]
  calibration_df <- make_calibration_data(pred, y_test)
  if (!is.null(calibration_df)) {
    calibration_df$model <- model_name
    calibration_df$outc <- outcome
    data.table::setcolorder(calibration_df, c("model", "outc", setdiff(names(calibration_df), c("model", "outc"))))
  }
  list(summary = summary_df, roc = roc_df, coef = coef_df, lambda = lambda_df,
       pred = pred_df, metrics = metrics_df, calibration = calibration_df)
}

all_auc <- list()
all_qc <- list()
all_metrics <- list()
all_calibration <- list()
all_predictions <- list()

set.seed(1234)
dt <- dt[sample(seq_len(nrow(dt))), , drop = FALSE]
dt$.lasso_split_id <- seq_len(nrow(dt))
set.seed(1)
train_split_ids <- sample(dt$.lasso_split_id, size = floor(0.80 * nrow(dt)))

for (i in seq_len(nrow(outcomes))) {
  key <- outcomes$outcome_key[i]
  label <- outcomes$label[i]
  message("Outcome: ", label, " (", key, ")")
  existing_files <- file.path(
    lasso_dir,
    paste0(key, c("_roc_rawdata.csv", "_coefficients.csv", "_lambda_path.csv",
                  "_predictions.csv", "_prediction_metrics.csv", "_calibration.csv"))
  )
  event_col <- paste0(key, "_inc")
  fu_col <- paste0(key, "_fu")
  pred_clinical <- repro_usable_covariates(dt, base$clinical_predictors_for_lasso)
  pred_protein <- intersect(proteins, names(dt))
  pred_combined <- unique(c(pred_clinical, pred_protein))
  d <- dt[!is.na(dt[[event_col]]) & !is.na(dt[[fu_col]]) & is.finite(dt[[fu_col]]) & dt[[fu_col]] > 0, , drop = FALSE]
  y <- as.integer(d[[event_col]] == 1)
  if (length(unique(y)) < 2) {
    warning("Skipping ", key, ": only one outcome class.")
    next
  }
  train_idx <- which(d$.lasso_split_id %in% train_split_ids)
  test_idx <- setdiff(seq_len(nrow(d)), train_idx)
  y_train <- y[train_idx]
  y_test <- y[test_idx]
  id_test <- if ("eid" %in% names(d)) d$eid[test_idx] else if ("id" %in% names(d)) d$id[test_idx] else test_idx
  fu_test <- d[[fu_col]][test_idx]

  if (!force_lasso && all(file.exists(existing_files))) {
    message("  existing LASSO outputs found; reusing ", key)
    existing_metrics <- data.table::fread(file.path(lasso_dir, paste0(key, "_prediction_metrics.csv")))
    existing_calibration <- data.table::fread(file.path(lasso_dir, paste0(key, "_calibration.csv")))
    existing_predictions <- data.table::fread(file.path(lasso_dir, paste0(key, "_predictions.csv")))
    all_metrics[[key]] <- existing_metrics
    all_calibration[[key]] <- existing_calibration
    all_predictions[[key]] <- existing_predictions
    all_auc[[key]] <- existing_metrics[, .(
      model, outc, auc, auc_lci, auc_uci,
      n_train = length(y_train),
      n_test,
      n_features = NA_integer_,
      error = NA_character_
    )]
    all_qc[[key]] <- data.frame(
      outcome = key,
      n_total = nrow(d),
      n_events = sum(y == 1),
      n_train = length(train_idx),
      n_test = length(test_idx),
      n_clinical_features = NA_integer_,
      n_protein_features = length(pred_protein),
      n_combined_features = NA_integer_,
      lasso_nfolds = lasso_nfolds,
      lasso_type_measure = lasso_type_measure,
      glmnet_maxit = glmnet_maxit,
      reused_existing_outputs = TRUE,
      stringsAsFactors = FALSE
    )
    next
  }

  x_cl_train <- make_x(d[train_idx, , drop = FALSE], pred_clinical)
  x_cl_test <- make_x(d[test_idx, , drop = FALSE], pred_clinical)
  x_pr_train <- make_x(d[train_idx, , drop = FALSE], pred_protein)
  x_pr_test <- make_x(d[test_idx, , drop = FALSE], pred_protein)
  x_cb_train <- make_x(d[train_idx, , drop = FALSE], pred_combined)
  x_cb_test <- make_x(d[test_idx, , drop = FALSE], pred_combined)

  fits <- list(
    clinical = fit_model_auc(x_cl_train, y_train, x_cl_test, y_test, key, "Clinical", "lambda.min", id_test = id_test, fu_test = fu_test),
    proteins = fit_model_auc(x_pr_train, y_train, x_pr_test, y_test, key, "Proteins", "lambda.1se", id_test = id_test, fu_test = fu_test),
    combined = fit_model_auc(x_cb_train, y_train, x_cb_test, y_test, key, "Combination", "lambda.1se", id_test = id_test, fu_test = fu_test)
  )
  auc <- do.call(rbind, lapply(fits, `[[`, "summary"))
  all_auc[[key]] <- auc

  roc <- do.call(rbind, Filter(Negate(is.null), lapply(fits, `[[`, "roc")))
  coef <- do.call(rbind, Filter(Negate(is.null), lapply(fits, `[[`, "coef")))
  lambda <- do.call(rbind, Filter(Negate(is.null), lapply(fits, `[[`, "lambda")))
  pred <- do.call(rbind, Filter(Negate(is.null), lapply(fits, `[[`, "pred")))
  metrics <- do.call(rbind, Filter(Negate(is.null), lapply(fits, `[[`, "metrics")))
  calibration <- data.table::rbindlist(Filter(Negate(is.null), lapply(fits, `[[`, "calibration")), fill = TRUE)
  if (!is.null(roc)) repro_write_csv(roc, file.path(lasso_dir, paste0(key, "_roc_rawdata.csv")))
  if (!is.null(coef)) repro_write_csv(coef, file.path(lasso_dir, paste0(key, "_coefficients.csv")))
  if (!is.null(lambda)) repro_write_csv(lambda, file.path(lasso_dir, paste0(key, "_lambda_path.csv")))
  if (!is.null(pred)) repro_write_csv(pred, file.path(lasso_dir, paste0(key, "_predictions.csv")))
  if (!is.null(metrics)) repro_write_csv(metrics, file.path(lasso_dir, paste0(key, "_prediction_metrics.csv")))
  if (nrow(calibration) > 0) repro_write_csv(calibration, file.path(lasso_dir, paste0(key, "_calibration.csv")))
  all_metrics[[key]] <- metrics
  all_calibration[[key]] <- calibration
  all_predictions[[key]] <- pred

  all_qc[[key]] <- data.frame(
    outcome = key,
    n_total = nrow(d),
    n_events = sum(y == 1),
    n_train = length(train_idx),
    n_test = length(test_idx),
    n_clinical_features = ncol(x_cl_train),
    n_protein_features = ncol(x_pr_train),
    n_combined_features = ncol(x_cb_train),
    lasso_nfolds = lasso_nfolds,
    lasso_type_measure = lasso_type_measure,
    glmnet_maxit = glmnet_maxit,
    reused_existing_outputs = FALSE,
    stringsAsFactors = FALSE
  )
  gc()
}

auc_full <- do.call(rbind, all_auc)
qc <- do.call(rbind, all_qc)
metrics_full <- do.call(rbind, all_metrics)
calibration_full <- data.table::rbindlist(all_calibration, fill = TRUE)
predictions_full <- do.call(rbind, all_predictions)

delong_pairs <- data.table::data.table(
  model = c("Combination", "Combination", "Proteins"),
  reference = c("Clinical", "Proteins", "Clinical")
)
delong_rows <- list()
if (!is.null(predictions_full) && nrow(predictions_full) > 0L) {
  pred_dt <- data.table::as.data.table(predictions_full)
  if (pred_dt[, anyDuplicated(paste(outc, model, eid, sep = "|"))] > 0L) {
    stop("Paired DeLong gate failed: duplicate outcome/model/EID predictions.", call. = FALSE)
  }
  for (outc_i in unique(pred_dt$outc)) {
    for (pair_i in seq_len(nrow(delong_pairs))) {
      model_i <- delong_pairs$model[[pair_i]]
      reference_i <- delong_pairs$reference[[pair_i]]
      a <- pred_dt[outc == outc_i & model == model_i, .(
        eid = as.character(eid), y, pred_model = predicted_risk
      )]
      b <- pred_dt[outc == outc_i & model == reference_i, .(
        eid = as.character(eid), y_ref = y, pred_reference = predicted_risk
      )]
      paired <- merge(a, b, by = "eid", all = FALSE)
      if (nrow(paired) == 0L || any(paired$y != paired$y_ref)) {
        stop("Paired DeLong gate failed: unmatched participants or outcomes for ", outc_i, ".", call. = FALSE)
      }
      roc_model <- pROC::roc(paired$y, paired$pred_model, levels = c(0, 1), quiet = TRUE)
      roc_reference <- pROC::roc(paired$y, paired$pred_reference, levels = c(0, 1), quiet = TRUE)
      test <- pROC::roc.test(roc_model, roc_reference, paired = TRUE, method = "delong")
      delong_rows[[paste(outc_i, model_i, reference_i, sep = "|")]] <- data.table::data.table(
        outc = outc_i,
        model = model_i,
        reference = reference_i,
        n_test = nrow(paired),
        n_events_test = sum(paired$y == 1L),
        auc_model = as.numeric(pROC::auc(roc_model)),
        auc_reference = as.numeric(pROC::auc(roc_reference)),
        delta_auc = as.numeric(pROC::auc(roc_model) - pROC::auc(roc_reference)),
        p_delong = as.numeric(test$p.value),
        method = "paired DeLong"
      )
    }
  }
}
delong_full <- data.table::rbindlist(delong_rows, fill = TRUE)

repro_write_csv(auc_full, lasso_auc_file)
if (!is.null(metrics_full)) repro_write_csv(metrics_full, file.path(output_dir, "lasso_prediction_metrics.csv"))
if (nrow(calibration_full) > 0) repro_write_csv(calibration_full, file.path(output_dir, "lasso_calibration.csv"))
if (!is.null(predictions_full)) repro_write_csv(predictions_full, file.path(output_dir, "lasso_predictions.csv"))
if (nrow(delong_full) > 0L) repro_write_csv(delong_full, file.path(output_dir, "lasso_delong_comparisons.csv"))
repro_write_csv(qc, file.path(audit_dir, "lasso_qc.csv"))

message("LASSO AUC summary saved: ", lasso_auc_file)
print(auc_full)
