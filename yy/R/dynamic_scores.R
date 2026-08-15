## Recalculate the frozen legacy Pradeep and Yu scores from learned model
## parameters plus participant-level protein data.  Pre-binned score curves
## and precomputed participant scores are never used as prediction inputs.

yy_score_first_existing <- function(paths, label) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) {
    stop(label, " not found. Checked: ", paste(paths, collapse = "; "), call. = FALSE)
  }
  normalizePath(hit[[1L]], winslash = "/", mustWork = TRUE)
}

yy_score_key <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

yy_score_resolve_yu_python <- function(dir0 = Sys.getenv("DIR0", unset = "/mnt/d")) {
  explicit <- Sys.getenv("YU_PYTHON", unset = "")
  candidates <- unique(c(
    explicit,
    file.path(dir0, "software", "conda", "envs", "yu_proteomic_repo_py39", "bin", "python"),
    file.path(dir0, "software", "conda", "envs", "yu_proteomic_repo_py39", "python.exe"),
    file.path(dir0, "software", "python", "yu_proteomic_repo_py39", "python.exe"),
    Sys.glob("/mnt/c/Users/*/anaconda3/envs/yu_proteomic_repo_py39/python.exe"),
    Sys.glob("/mnt/c/Users/*/miniconda3/envs/yu_proteomic_repo_py39/python.exe"),
    Sys.which("python3.9")
  ))
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  for (candidate in candidates) {
    version <- suppressWarnings(system2(
      candidate,
      c("-c", shQuote("import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")),
      stdout = TRUE, stderr = FALSE
    ))
    if (length(version) && identical(trimws(version[[1L]]), "3.9")) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  stop(
    "Yu Python 3.9 runtime not found. Set YU_PYTHON or install the environment under <DIR0>/software.",
    call. = FALSE
  )
}

yy_score_paths <- function(dir0, analysis_root, yy_outdir, script_root) {
  source_root <- file.path(yy_outdir, "score", "source-projects")
  pradeep_root <- yy_score_first_existing(
    c(Sys.getenv("PRADEEP_STRICT_MODEL_ROOT", unset = ""),
      file.path(source_root, "pradeep-strict"),
      file.path(analysis_root, "ukb", "pradeep_strict")),
    "Pradeep completed derived project"
  )
  yu_root <- yy_score_first_existing(
    c(Sys.getenv("YU_STRICT_MODEL_ROOT", unset = ""),
      file.path(source_root, "yu-strict"),
      file.path(analysis_root, "yu_strict")),
    "Yu completed derived project"
  )
  old_input_root <- yy_score_first_existing(
    file.path(source_root, "strict-projection-inputs"),
    "Public strict projection input project"
  )
  analysis_bundle <- yy_score_first_existing(
    file.path(source_root, "analysis-bundle", "analysis_bundle.rds"),
    "Locked Yin-Yang analysis bundle"
  )
  python_script <- yy_score_first_existing(
    file.path(script_root, "yy", "python", "replay_strict_yu.py"),
    "Yu frozen-score Python helper"
  )
  yu_raw_protein <- yy_score_first_existing(
    file.path(dir0, "data", "ukb", "phe", "raw", "prot_full_unimputed.tsv"),
    "Yu participant-level raw protein source"
  )
  list(
    analysis_bundle = analysis_bundle,
    pradeep_base = file.path(pradeep_root, "outputs", "ukbppp_cardiac_analysis_base.rds"),
    pradeep_coefficients = file.path(pradeep_root, "outputs", "lasso", "cad_coefficients.csv"),
    pradeep_predictions = file.path(pradeep_root, "outputs", "lasso", "cad_predictions.csv"),
    yu_model = file.path(yu_root, "09_models", "cad__Protein.txt"),
    yu_test_predictions = file.path(yu_root, "09_models", "test_predictions.csv.gz"),
    yu_yin_matrix = file.path(old_input_root, "yu_target_matrix.csv.gz"),
    yu_yang_matrix = file.path(old_input_root, "yu_yang_matrix.csv.gz"),
    yu_raw_protein = yu_raw_protein,
    validation_pradeep_yin = file.path(old_input_root, "locked_yin_target.rds"),
    validation_yu_yin = file.path(old_input_root, "yu_target_scores.csv.gz"),
    python_script = python_script,
    all_rds_candidates = file.path(dir0, "data", "ukb", "phe", "Rdata", "all.rds")
  )
}

yy_score_require_paths <- function(paths, names_required) {
  required <- unlist(paths[names_required], use.names = FALSE)
  missing <- required[!file.exists(required)]
  if (length(missing)) {
    stop("Locked score recalculation input missing: ", paste(missing, collapse = "; "), call. = FALSE)
  }
}

yy_score_map_names <- function(required, available, label) {
  required_key <- yy_score_key(required)
  available_key <- yy_score_key(available)
  if (anyDuplicated(available_key)) {
    duplicated_key <- unique(available_key[duplicated(available_key)])
    ambiguous <- intersect(required_key, duplicated_key)
    if (length(ambiguous)) {
      stop(label, " contains ambiguous names: ", paste(required[required_key %in% ambiguous], collapse = ", "), call. = FALSE)
    }
  }
  index <- match(required_key, available_key)
  if (anyNA(index)) {
    stop(label, " is missing: ", paste(required[is.na(index)], collapse = ", "), call. = FALSE)
  }
  available[index]
}

yy_score_linear <- function(dat, columns, coefficients, medians, intercept) {
  if (length(columns) != length(coefficients) || length(columns) != length(medians)) {
    stop("Linear score vectors have unequal lengths.", call. = FALSE)
  }
  score <- rep(as.numeric(intercept), nrow(dat))
  for (i in seq_along(columns)) {
    value <- suppressWarnings(as.numeric(dat[[columns[[i]]]]))
    value[!is.finite(value)] <- medians[[i]]
    score <- score + coefficients[[i]] * value
  }
  score
}

yy_score_z_reference <- function(x, reference) {
  reference <- as.numeric(reference)
  reference <- reference[is.finite(reference)]
  mu <- mean(reference)
  sigma <- stats::sd(reference)
  if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) {
    stop("Score reference distribution is invalid.", call. = FALSE)
  }
  (as.numeric(x) - mu) / sigma
}

yy_score_to_windows <- function(path) {
  ## Use mixed/forward-slash Windows paths.  Backslashes passed through
  ## system2() can be consumed as shell escapes before reaching python.exe.
  result <- system2("wslpath", c("-m", shQuote(normalizePath(path, winslash = "/", mustWork = FALSE))), stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) stop("wslpath failed for: ", path, call. = FALSE)
  trimws(result[[1L]])
}

yy_score_run_yu <- function(paths, temporary_root) {
  win_python <- yy_score_resolve_yu_python()
  dir.create(temporary_root, recursive = TRUE, showWarnings = FALSE)
  feature_out <- file.path(temporary_root, "yu_model_features.csv")
  yin_out <- file.path(temporary_root, "yu_yin_scores.csv.gz")
  yang_out <- file.path(temporary_root, "yu_yang_scores.csv.gz")
  test_out <- file.path(temporary_root, "yu_test_scores.csv.gz")
  yin_qc <- file.path(temporary_root, "yu_yin_validation.json")
  yang_qc <- file.path(temporary_root, "yu_yang_validation.json")
  test_qc <- file.path(temporary_root, "yu_test_validation.json")

  run_python <- function(arguments, label) {
    command_arguments <- as.character(arguments)
    names(command_arguments) <- names(arguments)
    path_arguments <- intersect(
      names(command_arguments),
      c(
        "model", "features-out", "matrix", "raw-source", "eid-filter",
        "scores-out", "test-predictions", "trusted-validation-qc", "qc-out"
      )
    )
    command_arguments[path_arguments] <- vapply(
      command_arguments[path_arguments], yy_score_to_windows, character(1L)
    )
    flattened <- unlist(Map(function(name, value) c(paste0("--", name), value), names(command_arguments), command_arguments), use.names = FALSE)
    output <- system2(win_python, args = c(yy_score_to_windows(paths$python_script), shQuote(flattened)), stdout = TRUE, stderr = TRUE)
    status <- attr(output, "status")
    if (!is.null(status) && status != 0L) {
      stop(label, " failed (exit ", status, "):\n", paste(output, collapse = "\n"), call. = FALSE)
    }
    output
  }

  ## Feature order is read from the fitted booster itself on every run.
  run_python(c(stage = "features", model = paths$yu_model, `features-out` = feature_out), "Yu feature extraction")
  run_python(c(
    stage = "score", model = paths$yu_model, matrix = paths$yu_yin_matrix,
    `scores-out` = yin_out, `test-predictions` = paths$yu_test_predictions,
    `qc-out` = yin_qc
  ), "Yu Yin score recalculation")
  run_python(c(
    stage = "score", model = paths$yu_model, matrix = paths$yu_yang_matrix,
    `scores-out` = yang_out, `trusted-validation-qc` = yin_qc,
    `qc-out` = yang_qc
  ), "Yu Yang score recalculation")
  ## Recalculate the complete Yu held-out cohort from the original participant-
  ## level raw protein file; the locked Yin matrix is a narrower YY cohort.
  run_python(c(
    stage = "score", model = paths$yu_model, `raw-source` = paths$yu_raw_protein,
    `eid-filter` = paths$yu_test_predictions,
    `scores-out` = test_out, `test-predictions` = paths$yu_test_predictions,
    `qc-out` = test_qc
  ), "Yu full held-out score recalculation")

  list(
    yin = data.table::fread(yin_out, colClasses = list(character = "eid")),
    yang = data.table::fread(yang_out, colClasses = list(character = "eid")),
    test = data.table::fread(test_out, colClasses = list(character = "eid")),
    features = data.table::fread(feature_out),
    yin_qc = jsonlite::fromJSON(yin_qc),
    yang_qc = jsonlite::fromJSON(yang_qc),
    test_qc = jsonlite::fromJSON(test_qc),
    generated_paths = c(feature_out, yin_out, yang_out, test_out, yin_qc, yang_qc, test_qc)
  )
}

yy_score_binary_roc <- function(y, score, series, cohort, model_definition) {
  y <- as.integer(y)
  score <- as.numeric(score)
  keep <- y %in% c(0L, 1L) & is.finite(score)
  y <- y[keep]
  score <- score[keep]
  if (length(unique(y)) != 2L) stop("Binary ROC has fewer than two classes: ", series, call. = FALSE)
  roc_object <- pROC::roc(y, score, levels = c(0, 1), direction = "<", quiet = TRUE)
  interval <- as.numeric(pROC::ci.auc(roc_object, method = "delong"))
  list(
    curve = data.table::data.table(
      series = series, cohort = cohort,
      false_positive_rate = 1 - as.numeric(roc_object$specificities),
      true_positive_rate = as.numeric(roc_object$sensitivities)
    ),
    metric = data.table::data.table(
      series = series, cohort = cohort, auc = as.numeric(pROC::auc(roc_object)),
      auc_ci_low = interval[[1L]], auc_ci_high = interval[[3L]],
      n = length(y), events = sum(y == 1L), controls = sum(y == 0L),
      model_definition = model_definition
    )
  )
}

yy_score_summarize_trajectory <- function(participants, labels, adjustment, paths) {
  spec <- yy_adjustment_spec(adjustment)
  extra_inputs <- character()
  if (spec$id == "raw") {
    trajectories <- participants[, {
      n_value <- .N
      mean_value <- mean(value)
      se_value <- stats::sd(value) / sqrt(n_value)
      list(
        n = n_value, mean = mean_value, se = se_value,
        low = mean_value - 1.96 * se_value,
        high = mean_value + 1.96 * se_value
      )
    }, by = .(side, relative_bin, series)]
  } else if (spec$id == "age+sex") {
    trajectories <- data.table::rbindlist(lapply(labels, function(label) {
      age_sex_adjusted_bins(participants[series == label], label, "Model score")
    }))
  } else {
    covariates <- yy_load_stepwise_covariates(
      list(all_rds_candidates = paths$all_rds_candidates), spec
    )
    extra_inputs <- attr(covariates, "input_path")
    participants <- merge(participants, covariates, by = "eid", all = FALSE, sort = FALSE)
    trajectories <- data.table::rbindlist(lapply(labels, function(label) {
      marginal_adjusted_bins(
        participants[series == label], label, "Model score",
        adjustment = spec$id, rhs_terms = spec$rhs
      )
    }))
  }
  trajectories[, `:=`(series_type = "Model score", bin_n = n)]
  list(trajectories = trajectories, extra_inputs = extra_inputs)
}

yy_legacy_score_bundle <- function(selected_score_names, adjustment, anchor, show_roc,
                                   dir0, analysis_root, yy_outdir, script_root) {
  allowed <- c("Pradeep local LASSO", "Yu local LightGBM")
  selected_score_names <- intersect(allowed, unique(as.character(selected_score_names)))
  if (!length(selected_score_names)) {
    return(list(
      trajectories = data.table::data.table(), roc = data.table::data.table(),
      metrics = data.table::data.table(), individual_scores = data.table::data.table(),
      parameter_audit = data.table::data.table(), qc = data.table::data.table(),
      input_paths = character()
    ))
  }
  if (!anchor %in% c("baseline", "diagnosis")) stop("Score anchor must be baseline or diagnosis.", call. = FALSE)
  paths <- yy_score_paths(dir0, analysis_root, yy_outdir, script_root)
  if (!exists("baseline_bin", mode = "function") || !exists("marginal_adjusted_bins", mode = "function")) {
    dynamic_paths <- yy_dynamic_paths(dir0, analysis_root, yy_outdir, script_root)
    source(dynamic_paths$core, local = FALSE)
  }
  needed <- c("analysis_bundle")
  if ("Pradeep local LASSO" %in% selected_score_names) {
    needed <- c(needed, "pradeep_base", "pradeep_coefficients", "pradeep_predictions")
  }
  if ("Yu local LightGBM" %in% selected_score_names) {
    needed <- c(
      needed, "yu_model", "yu_test_predictions", "yu_yin_matrix",
      "yu_yang_matrix", "yu_raw_protein", "python_script"
    )
  }
  yy_score_require_paths(paths, unique(needed))

  bundle <- readRDS(paths$analysis_bundle)
  yin <- data.table::as.data.table(bundle[["dat_yin"]])
  yang <- data.table::as.data.table(bundle[["dat_yang"]])
  protein_names <- as.character(bundle[["proteins"]])
  time_col <- as.character(bundle[["time_col"]])
  event_col <- as.character(bundle[["event_col"]])
  if (nrow(yin) != 37127L || sum(as.integer(yin[[event_col]]), na.rm = TRUE) != 3442L ||
      nrow(yang) != 1766L || anyDuplicated(as.character(yin$eid)) || anyDuplicated(as.character(yang$eid))) {
    stop("Locked Yin-Yang score cohort contract failed.", call. = FALSE)
  }
  yin[, eid := as.character(eid)]
  yang[, eid := as.character(eid)]
  if (!all(c("age", "sex") %in% names(yin)) || !all(c("age", "sex", "cvd_cad.b2e") %in% names(yang))) {
    stop("Locked score cohort lacks age, sex, or Yang timing.", call. = FALSE)
  }

  score_rows <- list()
  parameter_rows <- list()
  qc_rows <- list()
  roc_rows <- list()
  metric_rows <- list()
  input_paths <- c(paths$analysis_bundle)

  if ("Pradeep local LASSO" %in% selected_score_names) {
    coefficients_all <- data.table::fread(paths$pradeep_coefficients)
    coefficients <- coefficients_all[model == "Proteins" & outc == "cad",
      .(feature = as.character(feature), coefficient = as.numeric(coefficient))]
    intercept_rows <- coefficients[feature == "(Intercept)"]
    selected <- coefficients[feature != "(Intercept)" & coefficient != 0]
    if (nrow(intercept_rows) != 1L || nrow(selected) != 49L || anyDuplicated(selected$feature)) {
      stop("Pradeep learned-parameter contract failed.", call. = FALSE)
    }
    intercept <- intercept_rows$coefficient[[1L]]
    base <- readRDS(paths$pradeep_base)
    base_dat <- as.data.frame(base$dat)
    selected[, base_column := yy_score_map_names(feature, intersect(as.character(base$protein_cols), names(base_dat)), "Pradeep base")]
    selected[, yin_column := yy_score_map_names(feature, intersect(protein_names, names(yin)), "Pradeep Yin")]
    selected[, yang_column := yy_score_map_names(feature, intersect(protein_names, names(yang)), "Pradeep Yang")]

    set.seed(1234)
    base_dat <- base_dat[sample(seq_len(nrow(base_dat))), , drop = FALSE]
    base_dat$.lasso_split_id <- seq_len(nrow(base_dat))
    set.seed(1)
    train_split_ids <- sample(base_dat$.lasso_split_id, size = floor(0.80 * nrow(base_dat)))
    eligible <- !is.na(base_dat$cad_inc) & !is.na(base_dat$cad_fu) & is.finite(base_dat$cad_fu) & base_dat$cad_fu > 0
    d <- base_dat[eligible, , drop = FALSE]
    train_index <- which(d$.lasso_split_id %in% train_split_ids)
    test_index <- setdiff(seq_len(nrow(d)), train_index)
    training_medians <- vapply(selected$base_column, function(column) {
      value <- stats::median(suppressWarnings(as.numeric(d[[column]][train_index])), na.rm = TRUE)
      if (is.finite(value)) value else 0
    }, numeric(1L))
    test_medians <- vapply(selected$base_column, function(column) {
      value <- stats::median(suppressWarnings(as.numeric(d[[column]][test_index])), na.rm = TRUE)
      if (is.finite(value)) value else 0
    }, numeric(1L))
    test_score <- yy_score_linear(d[test_index, , drop = FALSE], selected$base_column, selected$coefficient, test_medians, intercept)
    frozen_test <- data.table::fread(paths$pradeep_predictions, select = c("model", "outc", "eid", "score_link"))
    frozen_test <- frozen_test[model == "Proteins" & outc == "cad", .(eid = as.character(eid), frozen_score = as.numeric(score_link))]
    validation <- merge(
      data.table::data.table(eid = as.character(d$eid[test_index]), recalculated_score = test_score),
      frozen_test, by = "eid", all = FALSE
    )
    validation[, delta := recalculated_score - frozen_score]
    if (nrow(validation) != length(test_index) || stats::cor(validation$recalculated_score, validation$frozen_score) < 0.999999 ||
        max(abs(validation$delta)) > 1e-6) {
      stop("Pradeep parameter-to-prediction validation failed.", call. = FALSE)
    }
    yin_raw <- yy_score_linear(yin, selected$yin_column, selected$coefficient, training_medians, intercept)
    yang_raw <- yy_score_linear(yang, selected$yang_column, selected$coefficient, training_medians, intercept)
    yin_z <- yy_score_z_reference(yin_raw, yin_raw)
    yang_z <- yy_score_z_reference(yang_raw, yin_raw)
    score_rows[["pradeep_yin"]] <- data.table::data.table(
      eid = yin$eid, cohort_side = "Yin", series = "Pradeep local LASSO",
      raw_value = yin_raw, model_scale_value = yin_raw, standardized_value = yin_z
    )
    score_rows[["pradeep_yang"]] <- data.table::data.table(
      eid = yang$eid, cohort_side = "Yang", series = "Pradeep local LASSO",
      raw_value = yang_raw, model_scale_value = yang_raw, standardized_value = yang_z
    )
    parameter_rows[["pradeep"]] <- data.table::rbindlist(list(
      data.table::data.table(
        model = "Pradeep local LASSO", parameter_type = "intercept", feature_order = 0L,
        feature = "(Intercept)", learned_value = intercept, imputation_value = NA_real_,
        definition = "Frozen protein-only binomial LASSO link score"
      ),
      selected[, .(
        model = "Pradeep local LASSO", parameter_type = "nonzero_coefficient",
        feature_order = seq_len(.N), feature, learned_value = coefficient,
        imputation_value = training_medians,
        definition = "Coefficient times protein value; target missing values use frozen-training median"
      )]
    ), use.names = TRUE)
    qc_rows[["pradeep_test"]] <- data.table::data.table(
      model = "Pradeep local LASSO", check = c("heldout_n", "heldout_pearson", "heldout_mae", "heldout_max_abs"),
      value = c(nrow(validation), stats::cor(validation$recalculated_score, validation$frozen_score),
                mean(abs(validation$delta)), max(abs(validation$delta))), status = "PASS"
    )
    if (file.exists(paths$validation_pradeep_yin)) {
      old <- data.table::as.data.table(readRDS(paths$validation_pradeep_yin))[, .(eid = as.character(eid), old_score = as.numeric(pradeep_score_raw))]
      check <- merge(data.table::data.table(eid = yin$eid, new_score = yin_raw), old, by = "eid", all = FALSE)
      check[, delta := new_score - old_score]
      if (nrow(check) != nrow(yin) || stats::cor(check$new_score, check$old_score) < 0.999999 || max(abs(check$delta)) > 1e-6) {
        stop("Pradeep locked-Yin implementation validation failed.", call. = FALSE)
      }
      qc_rows[["pradeep_yin"]] <- data.table::data.table(
        model = "Pradeep local LASSO", check = c("locked_yin_n", "locked_yin_pearson", "locked_yin_max_abs"),
        value = c(nrow(check), stats::cor(check$new_score, check$old_score), max(abs(check$delta))), status = "PASS"
      )
      input_paths <- c(input_paths, paths$validation_pradeep_yin)
    }
    if (isTRUE(show_roc)) {
      one_roc <- yy_score_binary_roc(
        d$cad_inc[test_index], test_score, "Pradeep local LASSO", "Pradeep legacy held-out test",
        "recalculated from 49 nonzero LASSO coefficients plus held-out individual proteins; conventional eventual-event ROC"
      )
      roc_rows[["pradeep"]] <- one_roc$curve
      metric_rows[["pradeep"]] <- one_roc$metric
    }
    input_paths <- c(input_paths, paths$pradeep_base, paths$pradeep_coefficients, paths$pradeep_predictions)
    rm(base, base_dat, d, frozen_test, validation); invisible(gc(FALSE))
  }

  if ("Yu local LightGBM" %in% selected_score_names) {
    temporary_root <- tempfile("yy-yu-recalculate-")
    yu <- yy_score_run_yu(paths, temporary_root)
    if (nrow(yu$features) != 359L || !identical(as.integer(yu$features$feature_order), seq_len(359L))) {
      stop("Yu fitted-booster feature contract failed.", call. = FALSE)
    }
    yin_index <- match(yin$eid, yu$yin$eid)
    yang_index <- match(yang$eid, yu$yang$eid)
    if (anyNA(yin_index) || anyNA(yang_index)) stop("Yu recalculated score EIDs do not align.", call. = FALSE)
    yin_probability <- as.numeric(yu$yin$yu_probability[yin_index])
    yang_probability <- as.numeric(yu$yang$yu_probability[yang_index])
    epsilon <- 1e-6
    yin_link <- stats::qlogis(pmin(pmax(yin_probability, epsilon), 1 - epsilon))
    yang_link <- stats::qlogis(pmin(pmax(yang_probability, epsilon), 1 - epsilon))
    yin_z <- yy_score_z_reference(yin_link, yin_link)
    yang_z <- yy_score_z_reference(yang_link, yin_link)
    score_rows[["yu_yin"]] <- data.table::data.table(
      eid = yin$eid, cohort_side = "Yin", series = "Yu local LightGBM",
      raw_value = yin_probability, model_scale_value = yin_link, standardized_value = yin_z
    )
    score_rows[["yu_yang"]] <- data.table::data.table(
      eid = yang$eid, cohort_side = "Yang", series = "Yu local LightGBM",
      raw_value = yang_probability, model_scale_value = yang_link, standardized_value = yang_z
    )
    parameter_rows[["yu"]] <- yu$features[, .(
      model = "Yu local LightGBM", parameter_type = "booster_feature_order",
      feature_order = as.integer(feature_order), feature = as.character(feature),
      learned_value = NA_real_, imputation_value = NA_real_,
      definition = "Feature order read directly from frozen fitted LightGBM booster; native missing-value routing"
    )]
    qc_rows[["yu"]] <- data.table::data.table(
      model = "Yu local LightGBM",
      check = c("feature_n", "validation_overlap_n", "validation_pearson", "validation_mae", "validation_max_abs"),
      value = c(
        yu$test_qc$feature_n, yu$test_qc$validation_overlap_n,
        yu$test_qc$validation_pearson, yu$test_qc$validation_mae, yu$test_qc$validation_max_abs
      ), status = "PASS"
    )
    if (file.exists(paths$validation_yu_yin)) {
      old <- data.table::fread(paths$validation_yu_yin, colClasses = list(character = "eid"))
      check <- merge(
        data.table::data.table(eid = yin$eid, new_score = yin_probability),
        old[, .(eid, old_score = as.numeric(yu_probability))], by = "eid", all = FALSE
      )
      check[, delta := new_score - old_score]
      if (nrow(check) != nrow(yin) || stats::cor(check$new_score, check$old_score) < 0.999999 || max(abs(check$delta)) > 1e-4) {
        stop("Yu locked-Yin implementation validation failed.", call. = FALSE)
      }
      qc_rows[["yu_yin"]] <- data.table::data.table(
        model = "Yu local LightGBM", check = c("locked_yin_n", "locked_yin_pearson", "locked_yin_max_abs"),
        value = c(nrow(check), stats::cor(check$new_score, check$old_score), max(abs(check$delta))), status = "PASS"
      )
      input_paths <- c(input_paths, paths$validation_yu_yin)
    }
    if (isTRUE(show_roc)) {
      labels <- data.table::fread(paths$yu_test_predictions, select = c("eid", "outcome_id", "y", "model_id"), colClasses = list(character = "eid"))
      labels <- labels[outcome_id == "cad" & model_id == "Protein", .(eid, y = as.integer(y))]
      evaluation <- merge(yu$test[, .(eid, score = as.numeric(yu_probability))], labels, by = "eid", all = FALSE)
      if (nrow(evaluation) != nrow(labels) || nrow(evaluation) < 15000L || anyDuplicated(evaluation$eid)) {
        stop("Yu fresh-score held-out overlap contract failed: n=", nrow(evaluation), call. = FALSE)
      }
      one_roc <- yy_score_binary_roc(
        evaluation$y, evaluation$score, "Yu local LightGBM", "Yu legacy held-out test",
        "recalculated from the frozen 359-feature LightGBM booster plus held-out individual proteins; conventional eventual-event ROC"
      )
      roc_rows[["yu"]] <- one_roc$curve
      metric_rows[["yu"]] <- one_roc$metric
    }
    input_paths <- c(
      input_paths, paths$yu_model, paths$yu_test_predictions,
      paths$yu_yin_matrix, paths$yu_yang_matrix, paths$yu_raw_protein, paths$python_script
    )
    unlink(temporary_root, recursive = TRUE, force = TRUE)
  }

  individual_scores <- data.table::rbindlist(score_rows, use.names = TRUE, fill = TRUE)
  labels <- selected_score_names
  duration <- -as.numeric(yang[["cvd_cad.b2e"]])
  if (any(!is.finite(duration)) || any(duration <= 0)) stop("Yang duration is invalid.", call. = FALSE)
  bins <- data.table::rbindlist(list(
    yin[as.integer(get(event_col)) == 1L, .(
      eid, side = "Yin", relative_bin = baseline_bin(as.numeric(get(time_col))),
      age = as.numeric(age), sex = as.integer(sex)
    )],
    yang[, .(
      eid, side = "Yang", relative_bin = baseline_bin(-duration),
      age = as.numeric(age), sex = as.integer(sex)
    )]
  ))
  if (anchor == "diagnosis") bins[, relative_bin := -relative_bin]
  bins[, side := if (anchor == "diagnosis") {
    data.table::fifelse(relative_bin < 0, "Yin", "Yang")
  } else {
    data.table::fifelse(relative_bin < 0, "Yang", "Yin")
  }]
  participants <- merge(
    bins,
    individual_scores[, .(eid, series, value = standardized_value)],
    by = "eid", all.x = TRUE, allow.cartesian = TRUE, sort = FALSE
  )
  participants <- participants[series %in% labels]
  if (nrow(participants) != nrow(bins) * length(labels) || any(!is.finite(participants$value))) {
    stop("Recalculated score trajectory participant contract failed.", call. = FALSE)
  }
  summary_result <- yy_score_summarize_trajectory(participants, labels, adjustment, paths)
  trajectories <- summary_result$trajectories
  trajectories[, displayed := !relative_bin %in% range(bins$relative_bin)]
  input_paths <- unique(c(input_paths, summary_result$extra_inputs))

  list(
    trajectories = trajectories,
    roc = data.table::rbindlist(roc_rows, use.names = TRUE, fill = TRUE),
    metrics = data.table::rbindlist(metric_rows, use.names = TRUE, fill = TRUE),
    individual_scores = individual_scores,
    parameter_audit = data.table::rbindlist(parameter_rows, use.names = TRUE, fill = TRUE),
    qc = data.table::rbindlist(qc_rows, use.names = TRUE, fill = TRUE),
    input_paths = unique(input_paths[file.exists(input_paths)])
  )
}

## Replay the frozen five-fold unified/fair-comparison models.  This is a
## parameter-to-person prediction audit only: no model is fitted, no lambda is
## selected, and no outcome is used to alter a saved model.
yy_unified_score_paths <- function(dir0, analysis_root, yy_outdir, script_root) {
  source_root <- file.path(yy_outdir, "score", "source-projects")
  unified_root <- yy_score_first_existing(
    file.path(source_root, "fair-unified"),
    "Public unified Yin-Yang trajectory/ROC project"
  )
  benchmark_root <- yy_score_first_existing(
    file.path(source_root, "fair-yu-benchmark"),
    "Public unified Yu benchmark project"
  )
  cache_root <- yy_score_first_existing(
    file.path(yy_outdir, "score", "common-fair-inputs"),
    "Public common Yin-Yang protein cache"
  )
  python_script <- yy_score_first_existing(
    file.path(script_root, "yy", "python", "replay_unified_yu.py"),
    "Unified Yu five-model replay helper"
  )
  dynamic <- yy_dynamic_paths(dir0, analysis_root, yy_outdir, script_root)
  list(
    unified_root = unified_root,
    benchmark_root = benchmark_root,
    cache_root = cache_root,
    participants_yin = file.path(cache_root, "participants_yin.csv"),
    participants_yang = file.path(cache_root, "participants_yang.csv"),
    features = file.path(cache_root, "protein_features.csv"),
    protein_yin = file.path(cache_root, "protein_yin.f32"),
    protein_yang = file.path(cache_root, "protein_yang.f32"),
    pradeep_models = file.path(
      unified_root, "02_models", sprintf("fold%02d_pradeep_single_cox.rds", 1:5)
    ),
    yu_models = file.path(
      benchmark_root, "03_models", sprintf("fold%02d_yu_models", 1:5),
      "Yu5_ProteinAll.txt"
    ),
    validation_yin = file.path(unified_root, "03_source_data", "unified_yin_oof_scores.csv.gz"),
    validation_yang = file.path(unified_root, "03_source_data", "unified_yang_ensemble_scores.csv.gz"),
    yin_target = dynamic$yin_target,
    yang_target = dynamic$yang_target,
    core = dynamic$core,
    all_rds_candidates = dynamic$all_rds_candidates,
    python_script = python_script
  )
}

yy_unified_replay_pradeep <- function(paths, yin, yang, features) {
  model_objects <- lapply(paths$pradeep_models, readRDS)
  for (fold in 1:5) {
    object <- model_objects[[fold]]
    if (!identical(as.integer(object$outer_fold), as.integer(fold)) ||
        is.null(object$pradeep$coefficient) || is.null(object$pradeep$center) ||
        is.null(object$pradeep$scale)) {
      stop("Unified Pradeep fold-model contract failed in fold ", fold, call. = FALSE)
    }
  }
  selected_by_fold <- lapply(model_objects, function(object) {
    coefficient <- object$pradeep$coefficient
    names(coefficient)[names(coefficient) != "(Intercept)" & coefficient != 0]
  })
  selected_features <- unique(unlist(selected_by_fold, use.names = FALSE))
  feature_index <- match(selected_features, features)
  if (!length(selected_features) || anyNA(feature_index) || anyDuplicated(feature_index)) {
    stop("Unified Pradeep selected-feature mapping failed.", call. = FALSE)
  }
  x_yin <- read_f32_fortran_columns(
    paths$protein_yin, nrow(yin), length(features), feature_index
  )
  x_yang <- read_f32_fortran_columns(
    paths$protein_yang, nrow(yang), length(features), feature_index
  )
  colnames(x_yin) <- colnames(x_yang) <- selected_features

  predict_one <- function(x, object) {
    coefficient <- object$pradeep$coefficient
    selected <- names(coefficient)[names(coefficient) != "(Intercept)" & coefficient != 0]
    value <- x[, selected, drop = FALSE]
    center <- as.numeric(object$pradeep$center[selected])
    scale <- as.numeric(object$pradeep$scale[selected])
    beta <- as.numeric(coefficient[selected])
    if (any(!is.finite(center)) || any(!is.finite(scale)) || any(scale <= 0) ||
        any(!is.finite(beta))) {
      stop("Unified Pradeep learned parameter is invalid.", call. = FALSE)
    }
    for (j in seq_along(selected)) {
      missing <- !is.finite(value[, j])
      if (any(missing)) value[missing, j] <- center[[j]]
    }
    value <- sweep(sweep(value, 2L, center, "-"), 2L, scale, "/")
    as.numeric(coefficient[["(Intercept)"]]) + drop(value %*% beta)
  }

  yin_rows <- vector("list", 5L)
  yang_rows <- vector("list", 5L)
  audit_rows <- vector("list", 5L)
  for (fold in 1:5) {
    object <- model_objects[[fold]]
    test_index <- which(as.integer(yin$outer_fold) == fold)
    yin_score <- predict_one(x_yin[test_index, , drop = FALSE], object)
    yang_score <- predict_one(x_yang, object)
    if (any(!is.finite(yin_score)) || any(!is.finite(yang_score))) {
      stop("Unified Pradeep replay produced non-finite scores in fold ", fold, call. = FALSE)
    }
    yin_rows[[fold]] <- data.table::data.table(
      eid = yin$eid[test_index], outer_fold = fold,
      time = as.numeric(yin$time[test_index]), event = as.integer(yin$event[test_index]),
      series = "Pradeep-style LASSO-logistic", score_raw = yin_score
    )
    yang_rows[[fold]] <- data.table::data.table(
      eid = yang$eid, outer_fold = fold,
      series = "Pradeep-style LASSO-logistic", score_raw = yang_score
    )
    nonzero <- selected_by_fold[[fold]]
    audit_rows[[fold]] <- data.table::data.table(
      model = "Pradeep-style LASSO-logistic", outer_fold = fold,
      model_path = paths$pradeep_models[[fold]],
      lambda = as.numeric(object$pradeep$lambda),
      training_n = as.integer(object$pradeep$training_n),
      training_cases5 = as.integer(object$pradeep$training_cases5),
      nonzero_feature_n = length(nonzero),
      intercept = as.numeric(object$pradeep$coefficient[["(Intercept)"]])
    )
  }
  yin_result <- data.table::rbindlist(yin_rows)
  yang_result <- data.table::rbindlist(yang_rows)
  if (nrow(yin_result) != nrow(yin) || anyDuplicated(yin_result$eid) ||
      nrow(yang_result) != 5L * nrow(yang) ||
      anyDuplicated(yang_result[, .(eid, outer_fold)])) {
    stop("Unified Pradeep replay coverage contract failed.", call. = FALSE)
  }
  rm(x_yin, x_yang, model_objects); invisible(gc(FALSE))
  list(
    yin = yin_result, yang = yang_result,
    audit = data.table::rbindlist(audit_rows),
    selected_features = selected_features
  )
}

yy_unified_replay_yu <- function(paths, temporary_root) {
  win_python <- yy_score_resolve_yu_python()
  dir.create(temporary_root, recursive = TRUE, showWarnings = FALSE)
  output_paths <- list(
    yin = file.path(temporary_root, "unified_yu_yin_oof.csv.gz"),
    yang = file.path(temporary_root, "unified_yu_yang_members.csv.gz"),
    metadata = file.path(temporary_root, "unified_yu_model_metadata.csv"),
    qc = file.path(temporary_root, "unified_yu_replay_qc.json")
  )
  arguments <- c(
    "--cache-root", yy_score_to_windows(paths$cache_root),
    "--benchmark-root", yy_score_to_windows(paths$benchmark_root),
    "--yin-output", yy_score_to_windows(output_paths$yin),
    "--yang-output", yy_score_to_windows(output_paths$yang),
    "--metadata-output", yy_score_to_windows(output_paths$metadata),
    "--qc-output", yy_score_to_windows(output_paths$qc)
  )
  output <- system2(
    win_python,
    args = c(shQuote(yy_score_to_windows(paths$python_script)), shQuote(arguments)),
    stdout = TRUE, stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("Unified Yu five-model replay failed (exit ", status, "):\n",
         paste(output, collapse = "\n"), call. = FALSE)
  }
  yin <- data.table::fread(output_paths$yin, colClasses = list(character = "eid"))
  yang <- data.table::fread(output_paths$yang, colClasses = list(character = "eid"))
  data.table::setnames(yin, "score", "score_raw")
  data.table::setnames(yang, "score", "score_raw")
  metadata <- data.table::fread(output_paths$metadata)
  metadata[, `:=`(
    model = "Yu-style LightGBM", model_path = as.character(model_path),
    training_n = NA_integer_, training_cases5 = NA_integer_,
    nonzero_feature_n = as.integer(feature_n), lambda = NA_real_, intercept = NA_real_
  )]
  audit <- metadata[, .(
    model, outer_fold, model_path, lambda, training_n, training_cases5,
    nonzero_feature_n, intercept, num_trees, best_iteration, lightgbm_version
  )]
  list(
    yin = yin, yang = yang, audit = audit,
    python_qc = jsonlite::fromJSON(output_paths$qc),
    generated_paths = unlist(output_paths, use.names = FALSE),
    console = output
  )
}

yy_unified_validate_replay <- function(yin_scores, yang_scores, paths, labels) {
  cached_yin <- data.table::fread(
    paths$validation_yin, colClasses = list(character = "eid")
  )[series %in% labels]
  cached_yang <- data.table::fread(
    paths$validation_yang, colClasses = list(character = "eid")
  )[series %in% labels]
  required <- c("eid", "series", "score_raw", "score_z")
  if (!all(required %in% names(cached_yin)) || !all(required %in% names(cached_yang))) {
    stop("Unified locked score cache lacks raw/z validation columns.", call. = FALSE)
  }
  validate_side <- function(replayed, cached, expected_n, side) {
    check <- merge(
      replayed[, .(eid, series, replay_raw = score_raw, replay_z = score_z)],
      cached[, .(eid, series, cache_raw = as.numeric(score_raw), cache_z = as.numeric(score_z))],
      by = c("eid", "series"), all = FALSE
    )
    if (nrow(check) != expected_n * length(labels) || anyDuplicated(check[, .(eid, series)])) {
      stop("Unified ", side, " replay-validation coverage failed.", call. = FALSE)
    }
    data.table::rbindlist(lapply(labels, function(label) {
      one <- check[series == label]
      raw_delta <- one$replay_raw - one$cache_raw
      z_delta <- one$replay_z - one$cache_z
      raw_cor <- stats::cor(one$replay_raw, one$cache_raw)
      z_cor <- stats::cor(one$replay_z, one$cache_z)
      tolerance <- if (label == "Yu-style LightGBM") 1e-5 else 1e-7
      if (!is.finite(raw_cor) || !is.finite(z_cor) || raw_cor < 0.999999 ||
          z_cor < 0.999999 || max(abs(raw_delta)) > tolerance ||
          max(abs(z_delta)) > tolerance) {
        stop("Unified ", label, " ", side,
             " parameter-to-cache validation failed.", call. = FALSE)
      }
      data.table::data.table(
        model = label, side = side,
        check = c("n", "raw_pearson", "raw_max_abs", "z_pearson", "z_max_abs"),
        value = c(nrow(one), raw_cor, max(abs(raw_delta)), z_cor, max(abs(z_delta))),
        status = "PASS"
      )
    }))
  }
  data.table::rbindlist(list(
    validate_side(yin_scores, cached_yin, 37127L, "Yin OOF"),
    validate_side(yang_scores, cached_yang, 1766L, "Yang five-model mean")
  ))
}

yy_unified_score_bundle <- function(selected_score_names, adjustment, anchor, show_roc,
                                    dir0, analysis_root, yy_outdir, script_root) {
  allowed <- c("Pradeep-style LASSO-logistic", "Yu-style LightGBM")
  labels <- intersect(allowed, unique(as.character(selected_score_names)))
  if (!length(labels)) {
    return(list(
      trajectories = data.table::data.table(), roc = data.table::data.table(),
      metrics = data.table::data.table(), individual_scores = data.table::data.table(),
      parameter_audit = data.table::data.table(), qc = data.table::data.table(),
      input_paths = character()
    ))
  }
  if (!anchor %in% c("baseline", "diagnosis")) {
    stop("Unified score anchor must be baseline or diagnosis.", call. = FALSE)
  }
  paths <- yy_unified_score_paths(dir0, analysis_root, yy_outdir, script_root)
  required <- c(
    paths$participants_yin, paths$participants_yang, paths$features,
    paths$protein_yin, paths$protein_yang, paths$validation_yin,
    paths$validation_yang, paths$yin_target, paths$yang_target, paths$core
  )
  if ("Pradeep-style LASSO-logistic" %in% labels) required <- c(required, paths$pradeep_models)
  if ("Yu-style LightGBM" %in% labels) required <- c(required, paths$yu_models, paths$python_script)
  missing <- required[!file.exists(required)]
  if (length(missing)) {
    stop("Unified frozen score replay input missing: ", paste(missing, collapse = "; "), call. = FALSE)
  }
  source(paths$core, local = FALSE)
  yin <- data.table::fread(paths$participants_yin, colClasses = list(character = "eid"))
  yang <- data.table::fread(paths$participants_yang, colClasses = list(character = "eid"))
  features <- as.character(data.table::fread(paths$features)$feature)
  if (nrow(yin) != 37127L || sum(as.integer(yin$event)) != 3442L ||
      nrow(yang) != 1766L || length(features) != 2910L ||
      anyDuplicated(yin$eid) || anyDuplicated(yang$eid) ||
      !identical(sort(unique(as.integer(yin$outer_fold))), 1:5)) {
    stop("Unified replay participant/feature contract failed.", call. = FALSE)
  }

  yin_parts <- yang_parts <- audit_parts <- list()
  input_paths <- c(
    paths$participants_yin, paths$participants_yang, paths$features,
    paths$protein_yin, paths$protein_yang, paths$validation_yin,
    paths$validation_yang, paths$yin_target, paths$yang_target, paths$core
  )
  temporary_root <- tempfile("yy-unified-score-replay-")
  on.exit(unlink(temporary_root, recursive = TRUE, force = TRUE), add = TRUE)
  if ("Pradeep-style LASSO-logistic" %in% labels) {
    message("Replaying five frozen Pradeep outer-fold models (no fitting).")
    replay <- yy_unified_replay_pradeep(paths, yin, yang, features)
    yin_parts[["pradeep"]] <- replay$yin
    yang_parts[["pradeep"]] <- replay$yang
    audit_parts[["pradeep"]] <- replay$audit
    input_paths <- c(input_paths, paths$pradeep_models)
  }
  if ("Yu-style LightGBM" %in% labels) {
    message("Replaying five frozen Yu outer-fold boosters (no fitting).")
    replay <- yy_unified_replay_yu(paths, temporary_root)
    yin_parts[["yu"]] <- replay$yin
    yang_parts[["yu"]] <- replay$yang
    audit_parts[["yu"]] <- replay$audit
    input_paths <- c(input_paths, paths$yu_models, paths$python_script)
  }
  yin_scores <- data.table::rbindlist(yin_parts, use.names = TRUE, fill = TRUE)
  yang_members <- data.table::rbindlist(yang_parts, use.names = TRUE, fill = TRUE)
  if (nrow(yin_scores) != nrow(yin) * length(labels) ||
      anyDuplicated(yin_scores[, .(eid, series)]) ||
      nrow(yang_members) != 5L * nrow(yang) * length(labels) ||
      anyDuplicated(yang_members[, .(eid, outer_fold, series)])) {
    stop("Unified replay combined coverage contract failed.", call. = FALSE)
  }
  reference <- yin_scores[, .(center = mean(score_raw), scale = stats::sd(score_raw)), by = series]
  if (any(!is.finite(reference$center)) || any(!is.finite(reference$scale)) || any(reference$scale <= 0)) {
    stop("Unified replay Yin reference standardization failed.", call. = FALSE)
  }
  yin_scores <- reference[yin_scores, on = "series"]
  yin_scores[, score_z := (score_raw - center) / scale]
  yang_scores <- yang_members[, .(score_raw = mean(score_raw)), by = .(eid, series)]
  yang_scores <- reference[yang_scores, on = "series"]
  yang_scores[, score_z := (score_raw - center) / scale]
  qc <- yy_unified_validate_replay(yin_scores, yang_scores, paths, labels)

  yin_target <- data.table::as.data.table(readRDS(paths$yin_target)); yin_target[, eid := as.character(eid)]
  yang_target <- data.table::as.data.table(readRDS(paths$yang_target)); yang_target[, eid := as.character(eid)]
  if (nrow(yin_target) != 37127L || sum(as.integer(yin_target$event)) != 3442L ||
      nrow(yang_target) != 1766L || anyDuplicated(yin_target$eid) || anyDuplicated(yang_target$eid)) {
    stop("Unified replay target contract failed.", call. = FALSE)
  }
  base_bins <- data.table::rbindlist(list(
    yin_target[event == 1L, .(
      eid, side = "Yin", relative_bin = baseline_bin(as.numeric(time_years))
    )],
    yang_target[, .(
      eid, side = "Yang", relative_bin = baseline_bin(-as.numeric(disease_duration_years))
    )]
  ))
  if (anchor == "diagnosis") base_bins[, relative_bin := -relative_bin]
  base_bins[, side := if (anchor == "diagnosis") {
    data.table::fifelse(relative_bin < 0, "Yin", "Yang")
  } else {
    data.table::fifelse(relative_bin < 0, "Yang", "Yin")
  }]
  score_cases <- data.table::rbindlist(list(
    yin_scores[, .(eid, series, value = score_z)],
    yang_scores[, .(eid, series, value = score_z)]
  ))
  covariates <- data.table::rbindlist(list(
    yin[, .(eid, age = as.numeric(age), sex = as.integer(sex))],
    yang[, .(eid, age = as.numeric(age), sex = as.integer(sex))]
  ))
  participants <- merge(base_bins, score_cases, by = "eid", allow.cartesian = TRUE, sort = FALSE)
  participants <- merge(participants, covariates, by = "eid", all.x = TRUE, sort = FALSE)
  if (nrow(participants) != nrow(base_bins) * length(labels) ||
      any(!is.finite(participants$value)) || anyNA(participants$age) || anyNA(participants$sex)) {
    stop("Unified replay trajectory participant contract failed.", call. = FALSE)
  }
  trajectory_result <- yy_score_summarize_trajectory(participants, labels, adjustment, paths)
  trajectories <- trajectory_result$trajectories
  trajectories[, displayed := !relative_bin %in% range(base_bins$relative_bin)]
  input_paths <- c(input_paths, trajectory_result$extra_inputs)

  roc <- metrics <- data.table::data.table()
  if (isTRUE(show_roc)) {
    adjustment_spec <- yy_adjustment_spec(adjustment)
    roc_covariates <- NULL
    if (adjustment_spec$id == "age+sex") {
      roc_covariates <- yin[, .(eid, age = as.numeric(age), sex = as.integer(sex))]
      roc_covariates[, sex_factor := factor(
        sex, levels = c(0L, 1L), labels = c("Female", "Male")
      )]
    } else if (adjustment_spec$id != "raw") {
      roc_covariates <- yy_load_stepwise_covariates(paths, adjustment_spec)
      input_paths <- c(input_paths, attr(roc_covariates, "input_path"))
    }
    roc_objects <- lapply(labels, function(label) {
      one <- yin_scores[series == label, .(
        eid, outer_fold, time, event, score = score_z
      )]
      if (!is.null(roc_covariates)) {
        one <- merge(one, roc_covariates, by = "eid", all = FALSE, sort = FALSE)
        one[, score := cross_fitted_covariate_residual(.SD, adjustment_spec$rhs)]
      }
      result <- mean_fold_roc(one[, .(eid, outer_fold, time, event, score)], 5)
      result$curve[, series := label]
      list(
        curve = result$curve,
        metric = data.table::data.table(
          series = label, AUC_5y = result$mean_fold_auc,
          n = nrow(one), incident_events = sum(as.integer(one$event)),
          model_definition = paste0(
            "five frozen outer-fold models replayed; matched five-year IPCW ROC; adjustment=",
            adjustment_spec$id, "; no refitting"
          )
        )
      )
    })
    roc <- data.table::rbindlist(lapply(roc_objects, `[[`, "curve"))
    metrics <- data.table::rbindlist(lapply(roc_objects, `[[`, "metric"))
  }

  individual_scores <- data.table::rbindlist(list(
    yin_scores[, .(
      eid, cohort_side = "Yin", outer_fold, time, event, series,
      raw_value = score_raw, standardized_value = score_z
    )],
    yang_scores[, .(
      eid, cohort_side = "Yang", outer_fold = NA_integer_,
      time = NA_real_, event = NA_integer_, series,
      raw_value = score_raw, standardized_value = score_z
    )]
  ), use.names = TRUE, fill = TRUE)
  list(
    trajectories = trajectories, roc = roc, metrics = metrics,
    individual_scores = individual_scores,
    parameter_audit = data.table::rbindlist(audit_parts, use.names = TRUE, fill = TRUE),
    qc = qc,
    input_paths = unique(input_paths[file.exists(input_paths)])
  )
}
