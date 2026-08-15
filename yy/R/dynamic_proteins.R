## On-demand protein trajectories and single-protein Cox ROC for `yy plot`.
## This file intentionally contains no command-line side effects.

yy_first_existing <- function(paths, label) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) {
    stop(label, " not found. Checked: ", paste(paths, collapse = "; "), call. = FALSE)
  }
  normalizePath(hit[[1L]], winslash = "/", mustWork = TRUE)
}

yy_protein_key <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

yy_legacy_pradeep_paths <- function(dir0, analysis_root) {
  project_root_default <- file.path(
    analysis_root, "yy", "score", "source-projects", "pradeep-strict"
  )
  project_root_override <- Sys.getenv("PRADEEP_STRICT_MODEL_ROOT", unset = "")
  project_root <- yy_first_existing(
    c(project_root_override[nzchar(project_root_override)], project_root_default,
      file.path(analysis_root, "ukb", "pradeep_strict")),
    "Pradeep completed derived project"
  )
  list(
    project_root = project_root,
    base = file.path(project_root, "outputs", "ukbppp_cardiac_analysis_base.rds"),
    predictions = file.path(project_root, "outputs", "lasso", "cad_predictions.csv")
  )
}

yy_resolve_legacy_pradeep_proteins <- function(requested, protein_cols) {
  requested_key <- yy_protein_key(requested)
  protein_key <- yy_protein_key(protein_cols)
  duplicate_key <- unique(protein_key[duplicated(protein_key)])
  ambiguous <- intersect(requested_key, duplicate_key)
  if (length(ambiguous)) {
    stop(
      "Protein alias is ambiguous in the 1,459-protein Pradeep legacy panel: ",
      paste(requested[requested_key %in% ambiguous], collapse = ", "),
      call. = FALSE
    )
  }
  index <- match(requested_key, protein_key)
  if (anyNA(index)) {
    stop(
      "Legacy Cox ROC cannot be computed because the original Pradeep panel lacks: ",
      paste(requested[is.na(index)], collapse = ", "),
      ". The strict Pradeep model contains 1,459 mapped proteins. ",
      "Use the common single-protein Cox layer for proteins outside that strict panel.",
      call. = FALSE
    )
  }
  data.table::data.table(
    requested = requested,
    feature = protein_cols[index],
    series = requested
  )
}

yy_legacy_pradeep_cox_roc <- function(requested, adjustment, dir0, analysis_root) {
  if (!length(requested)) {
    return(list(
      roc = data.table::data.table(), metrics = data.table::data.table(),
      mapping = data.table::data.table(), input_paths = character()
    ))
  }
  if (!adjustment %in% c("raw", "age+sex", "age_sex")) {
    stop(
      "Legacy on-demand Cox ROC supports --adj raw or --adj=age+sex only. ",
      "Higher adjustment layers are available for trajectories but were not part of the original Pradeep split.",
      call. = FALSE
    )
  }
  paths <- yy_legacy_pradeep_paths(dir0, analysis_root)
  required <- unlist(paths[c("base", "predictions")], use.names = FALSE)
  missing <- required[!file.exists(required)]
  if (length(missing)) {
    stop("Pradeep legacy Cox input missing: ", paste(missing, collapse = "; "), call. = FALSE)
  }

  bundle <- readRDS(paths$base)
  if (!all(c("dat", "protein_cols") %in% names(bundle))) {
    stop("Pradeep legacy base RDS contract failed.", call. = FALSE)
  }
  mapping <- yy_resolve_legacy_pradeep_proteins(requested, bundle$protein_cols)
  required_columns <- unique(c(
    "eid", "cad_fu", "cad_inc", "age", "Sex_numeric", mapping$feature
  ))
  if (!all(required_columns %in% names(bundle$dat))) {
    stop(
      "Pradeep legacy base lacks Cox columns: ",
      paste(setdiff(required_columns, names(bundle$dat)), collapse = ", "),
      call. = FALSE
    )
  }
  d <- data.table::as.data.table(bundle$dat)[, ..required_columns]
  rm(bundle); invisible(gc(FALSE))
  d[, eid := as.character(eid)]

  ## Reproduce 04_lasso_risk_score.R exactly: global shuffle, global 80/20 split,
  ## then outcome-specific validity filtering. The stored prediction eids below
  ## are used as a hard contract, not to train the protein Cox models.
  set.seed(1234)
  d <- d[sample(seq_len(nrow(d)))]
  d[, legacy_split_id := .I]
  set.seed(1)
  train_split_ids <- sample(d$legacy_split_id, floor(0.80 * nrow(d)))
  d <- d[is.finite(cad_fu) & cad_fu > 0 & !is.na(cad_inc)]
  train <- d[legacy_split_id %in% train_split_ids]
  test <- d[!legacy_split_id %in% train_split_ids]
  if (anyDuplicated(train$eid) || anyDuplicated(test$eid) || length(intersect(train$eid, test$eid))) {
    stop("Pradeep legacy train/test partition is not disjoint and unique.", call. = FALSE)
  }

  stored <- data.table::fread(paths$predictions, colClasses = list(character = "eid"))
  stored <- stored[model == "Proteins" & outc == "cad"]
  if (
    nrow(test) != 8983L || sum(test$cad_inc) != 590L ||
    nrow(stored) != nrow(test) || sum(stored$y) != sum(test$cad_inc) ||
    length(setdiff(test$eid, stored$eid)) || length(setdiff(stored$eid, test$eid))
  ) {
    stop("Pradeep legacy held-out eid/event contract failed.", call. = FALSE)
  }

  fit_one <- function(feature, series) {
    raw_train <- as.numeric(train[[feature]])
    median_train <- stats::median(raw_train[is.finite(raw_train)])
    if (!is.finite(median_train)) stop("All training values are missing: ", series, call. = FALSE)
    raw_train[!is.finite(raw_train)] <- median_train
    center_train <- mean(raw_train)
    scale_train <- stats::sd(raw_train)
    if (!is.finite(scale_train) || scale_train <= 0) {
      stop("Invalid legacy training SD: ", series, call. = FALSE)
    }
    z_train <- (raw_train - center_train) / scale_train
    fit <- survival::coxph(
      survival::Surv(train$cad_fu, train$cad_inc) ~ z_train,
      ties = "efron", x = FALSE, y = FALSE, model = FALSE
    )
    beta <- as.numeric(stats::coef(fit)[[1L]])
    if (!is.finite(beta)) stop("Legacy Cox coefficient is not finite: ", series, call. = FALSE)

    raw_test <- as.numeric(test[[feature]])
    raw_test[!is.finite(raw_test)] <- median_train
    score_train <- beta * z_train
    score_test <- beta * ((raw_test - center_train) / scale_train)
    if (adjustment %in% c("age+sex", "age_sex")) {
      residual_train <- data.frame(
        score = score_train,
        age = train$age,
        sex_factor = factor(train$Sex_numeric, levels = c(0, 1))
      )
      residual_test <- data.frame(
        age = test$age,
        sex_factor = factor(test$Sex_numeric, levels = c(0, 1))
      )
      residual_model <- stats::lm(
        score ~ splines::ns(age, df = 3) + sex_factor,
        data = residual_train,
        na.action = stats::na.fail
      )
      expected_test <- stats::predict(residual_model, newdata = residual_test)
      score_test <- score_test - as.numeric(expected_test)
    }
    if (any(!is.finite(score_test))) {
      stop("Legacy held-out Cox score is not finite: ", series, call. = FALSE)
    }
    roc_object <- pROC::roc(
      response = test$cad_inc, predictor = score_test,
      levels = c(0, 1), direction = "<", quiet = TRUE
    )
    auc_ci <- as.numeric(pROC::ci.auc(roc_object, method = "delong"))
    list(
      curve = data.table::data.table(
        series = series,
        cohort = "Pradeep legacy held-out test",
        false_positive_rate = 1 - as.numeric(roc_object$specificities),
        true_positive_rate = as.numeric(roc_object$sensitivities)
      ),
      metric = data.table::data.table(
        series = series,
        cohort = "Pradeep legacy held-out test",
        auc = as.numeric(pROC::auc(roc_object)),
        auc_ci_low = auc_ci[[1L]], auc_ci_high = auc_ci[[3L]],
        n = nrow(test), events = sum(test$cad_inc), controls = sum(test$cad_inc == 0L),
        cox_beta_per_training_sd = beta,
        cox_hr_per_training_sd = exp(beta),
        model_definition = paste0(
          "single-protein Cox trained on legacy 80%; ",
          if (adjustment %in% c("age+sex", "age_sex")) "held-out score residualized for age spline and sex" else "unadjusted",
          "; conventional eventual-event ROC on held-out 20%"
        )
      )
    )
  }

  results <- Map(fit_one, mapping$feature, mapping$series)
  list(
    roc = data.table::rbindlist(lapply(results, `[[`, "curve")),
    metrics = data.table::rbindlist(lapply(results, `[[`, "metric")),
    mapping = mapping,
    input_paths = c(paths$base, paths$predictions)
  )
}

yy_dynamic_paths <- function(dir0, analysis_root, yy_outdir, script_root) {
  public_common <- yy_first_existing(
    c(
      file.path(yy_outdir, "score", "common-fair-inputs", "cad"),
      file.path(yy_outdir, "score", "common-fair-inputs")
    ),
    "Public common Yin-Yang protein cache"
  )
  use_public_common <- file.exists(file.path(public_common, "COMPLETE"))
  if (use_public_common) {
    cache_root <- public_common
    figure_input_root <- yy_first_existing(
      c(public_common,
        file.path(yy_outdir, "score", "source-projects", "strict-projection-inputs")),
      "Public locked Yin-Yang target root"
    )
    if (!all(file.exists(file.path(
      figure_input_root, c("locked_yin_target.rds", "locked_yang_target.rds")
    )))) {
      figure_input_root <- yy_first_existing(
        file.path(yy_outdir, "score", "source-projects", "strict-projection-inputs"),
        "Public locked Yin-Yang target root"
      )
    }
    if (all(file.exists(file.path(cache_root, c(
      "participants_yin.csv.gz", "participants_yang.csv.gz"
    ))))) {
      participant_suffix <- ".csv.gz"
    } else if (all(file.exists(file.path(cache_root, c(
      "participants_yin.csv", "participants_yang.csv"
    ))))) {
      participant_suffix <- ".csv"
    } else {
      stop("Public common input has no matched Yin/Yang participant pair: ",
           cache_root, call. = FALSE)
    }
  } else {
    cache_root <- yy_first_existing(
      file.path(yy_outdir, "score", "source-projects", "common-yinyang-cache"),
      "Public Yin-Yang derived protein cache"
    )
    figure_input_root <- yy_first_existing(
      file.path(yy_outdir, "score", "source-projects", "strict-projection-inputs"),
      "Public locked Yin-Yang target root"
    )
    participant_suffix <- ".csv"
  }
  core_path <- yy_first_existing(
    c(
      file.path(script_root, "yy", "R", "core.R"),
      file.path(script_root, "yy", "cad_yinyang_trajectory_roc_unified_v1", "R", "core.R")
    ),
    "Yin-Yang trajectory core"
  )
  all_rds_candidates <- file.path(dir0, "data", "ukb", "phe", "Rdata", "all.rds")
  list(
    cache_root = cache_root,
    cache_source = if (use_public_common) "public_common_contract" else "public_derived_source_project",
    participants_yin = file.path(cache_root, paste0("participants_yin", participant_suffix)),
    participants_yang = file.path(cache_root, paste0("participants_yang", participant_suffix)),
    features = file.path(cache_root, "protein_features.csv"),
    protein_yin = file.path(cache_root, "protein_yin.f32"),
    protein_yang = file.path(cache_root, "protein_yang.f32"),
    yin_target = file.path(figure_input_root, "locked_yin_target.rds"),
    yang_target = file.path(figure_input_root, "locked_yang_target.rds"),
    all_rds_candidates = all_rds_candidates,
    core = core_path
  )
}

yy_dynamic_feature_catalog <- function(paths) {
  required <- unlist(paths[c(
    "participants_yin", "participants_yang", "features", "protein_yin",
    "protein_yang", "yin_target", "yang_target", "core"
  )], use.names = FALSE)
  missing <- required[!file.exists(required)]
  if (length(missing)) {
    stop("On-demand protein input missing: ", paste(missing, collapse = "; "), call. = FALSE)
  }
  data.table::fread(paths$features)$feature
}

yy_dynamic_display_name <- function(feature) {
  result <- as.character(feature)
  result[protein_key(result) == "NTPROBNP"] <- "NT-proBNP"
  result
}

yy_resolve_dynamic_proteins <- function(requested, features) {
  feature_key <- protein_key(features)
  requested_key <- protein_key(requested)
  if (any(!nzchar(requested_key))) stop("Protein names cannot be empty.", call. = FALSE)
  if (anyDuplicated(feature_key)) {
    duplicate_key <- unique(feature_key[duplicated(feature_key)])
    ambiguous <- intersect(requested_key, duplicate_key)
    if (length(ambiguous)) {
      stop("Protein alias is ambiguous in the 2,910-feature map: ",
           paste(requested[requested_key %in% ambiguous], collapse = ", "), call. = FALSE)
    }
  }
  index <- match(requested_key, feature_key)
  if (anyNA(index)) {
    unavailable <- requested[is.na(index)]
    stop(
      "Protein not found in the locked 2,910-feature matrix: ",
      paste(unavailable, collapse = ", "),
      ". Run `yy plot --status` or check the Olink feature name.",
      call. = FALSE
    )
  }
  data.table::data.table(
    requested = requested,
    feature = features[index],
    feature_index = as.integer(index),
    series = yy_dynamic_display_name(features[index])
  )
}

yy_adjustment_spec <- function(adjustment) {
  value <- tolower(trimws(as.character(adjustment[[1L]])))
  value <- gsub("[ -]", "_", value)
  aliases <- c(
    none = "raw", unadjusted = "raw", age_sex = "age+sex", agesex = "age+sex",
    age_sex_pc = "age+sex+pc",
    center_tdi = "age+sex+pc+center+tdi"
  )
  if (value %in% names(aliases)) value <- unname(aliases[[value]])
  if (value == "raw") return(list(id = "raw", rhs = "", tokens = character(), raw_columns = character()))
  tokens <- unique(trimws(unlist(strsplit(value, "\\+", perl = TRUE), use.names = FALSE)))
  tokens <- tokens[nzchar(tokens)]
  allowed <- c("age", "sex", "pc", "center", "tdi", "le8", "medications")
  invalid <- setdiff(tokens, allowed)
  if (!length(tokens) || length(invalid)) {
    stop(
      "Unsupported --adj=", adjustment,
      ". Use raw or a '+' combination of: ", paste(allowed, collapse = ", "),
      call. = FALSE
    )
  }
  canonical_order <- allowed[allowed %in% tokens]
  terms <- character()
  raw_columns <- character()
  if ("age" %in% tokens) {
    terms <- c(terms, "splines::ns(age, df=3)")
    raw_columns <- c(raw_columns, "age")
  }
  if ("sex" %in% tokens) {
    terms <- c(terms, "sex_factor")
    raw_columns <- c(raw_columns, "sex")
  }
  if ("pc" %in% tokens) {
    terms <- c(terms, "PC1", "PC2")
    raw_columns <- c(raw_columns, "PC1", "PC2")
  }
  if ("center" %in% tokens) {
    terms <- c(terms, "center_factor")
    raw_columns <- c(raw_columns, "center")
  }
  if ("tdi" %in% tokens) {
    terms <- c(terms, "tdi")
    raw_columns <- c(raw_columns, "tdi")
  }
  if ("le8" %in% tokens) {
    le8 <- c("diet.pts", "pa.pts", "smoke.pts", "bmi.pts", "nonhdl.pts",
             "hba1c.pts", "bp.pts", "sleep.pts")
    terms <- c(terms, le8)
    raw_columns <- c(raw_columns, le8)
  }
  if ("medications" %in% tokens) {
    medications <- c("drug.lipid", "drug.dm", "drug.htn")
    terms <- c(terms, medications)
    raw_columns <- c(raw_columns, medications)
  }
  list(
    id = paste(canonical_order, collapse = "+"),
    rhs = paste(terms, collapse = " + "),
    tokens = canonical_order,
    raw_columns = unique(raw_columns)
  )
}

yy_load_stepwise_covariates <- function(paths, adjustment_spec = NULL) {
  all_rds <- yy_first_existing(paths$all_rds_candidates, "Phenotype all.rds for stepwise adjustment")
  if (is.null(adjustment_spec)) {
    raw_covariates <- c(
      "age", "sex", "PC1", "PC2", "center", "tdi",
      "diet.pts", "pa.pts", "smoke.pts", "bmi.pts", "nonhdl.pts",
      "hba1c.pts", "bp.pts", "sleep.pts", "drug.lipid", "drug.dm", "drug.htn"
    )
  } else {
    raw_covariates <- as.character(adjustment_spec$raw_columns)
  }
  all_dat <- data.table::as.data.table(readRDS(all_rds))
  required <- c("eid", raw_covariates)
  if (!all(required %in% names(all_dat))) {
    stop("Phenotype all.rds lacks dynamic-adjustment columns: ",
         paste(setdiff(required, names(all_dat)), collapse = ", "), call. = FALSE)
  }
  covariates <- all_dat[, ..required]
  rm(all_dat); invisible(gc(FALSE))
  covariates[, eid := as.character(eid)]
  covariates <- covariates[stats::complete.cases(covariates)]
  if ("sex" %in% names(covariates)) {
    covariates[, sex_factor := factor(sex, levels = c(0L, 1L), labels = c("Female", "Male"))]
  }
  if ("center" %in% names(covariates)) covariates[, center_factor := factor(center)]
  if (anyDuplicated(covariates$eid)) stop("Complete covariate eid is not unique.", call. = FALSE)
  attr(covariates, "input_path") <- all_rds
  covariates
}

yy_fit_single_protein_oof <- function(yin, protein_yin, labels) {
  result <- vector("list", length(labels))
  for (j in seq_along(labels)) {
    fold_rows <- vector("list", 5L)
    for (fold in 1:5) {
      train_index <- which(yin$outer_fold != fold)
      test_index <- which(yin$outer_fold == fold)
      raw_train <- protein_yin[train_index, j]
      center <- mean(raw_train, na.rm = TRUE)
      if (!is.finite(center)) center <- 0
      raw_train[!is.finite(raw_train)] <- center
      scale <- stats::sd(raw_train)
      if (!is.finite(scale) || scale <= 0) stop("Invalid training SD: ", labels[[j]], call. = FALSE)
      fit <- survival::coxph(
        survival::Surv(yin$time[train_index], yin$event[train_index]) ~
          I((raw_train - center) / scale),
        ties = "efron"
      )
      beta <- as.numeric(stats::coef(fit)[[1L]])
      raw_test <- protein_yin[test_index, j]
      raw_test[!is.finite(raw_test)] <- center
      fold_rows[[fold]] <- data.table::data.table(
        eid = yin$eid[test_index], outer_fold = fold,
        time = yin$time[test_index], event = yin$event[test_index],
        score = beta * ((raw_test - center) / scale)
      )
    }
    result[[j]] <- data.table::rbindlist(fold_rows)[, series := labels[[j]]]
  }
  data.table::rbindlist(result)
}

yy_dynamic_protein_bundle <- function(requested, adjustment, version, anchor, show_roc,
                                      dir0, analysis_root, yy_outdir, script_root) {
  paths <- yy_dynamic_paths(dir0, analysis_root, yy_outdir, script_root)
  source(paths$core, local = FALSE)
  features <- yy_dynamic_feature_catalog(paths)
  mapping <- yy_resolve_dynamic_proteins(requested, features)

  yin <- data.table::fread(paths$participants_yin, colClasses = list(character = "eid"))
  yang <- data.table::fread(paths$participants_yang, colClasses = list(character = "eid"))
  yin_target <- data.table::as.data.table(readRDS(paths$yin_target)); yin_target[, eid := as.character(eid)]
  yang_target <- data.table::as.data.table(readRDS(paths$yang_target)); yang_target[, eid := as.character(eid)]
  if (nrow(yin) != 37127L || sum(yin$event) != 3442L || nrow(yang) != 1766L ||
      length(features) != 2910L || anyDuplicated(yin$eid) || anyDuplicated(yang$eid) ||
      nrow(yin_target) != 37127L || nrow(yang_target) != 1766L) {
    stop("Locked dynamic-protein cohort contract failed.", call. = FALSE)
  }

  protein_yin <- read_f32_fortran_columns(
    paths$protein_yin, nrow(yin), length(features), mapping$feature_index
  )
  protein_yang <- read_f32_fortran_columns(
    paths$protein_yang, nrow(yang), length(features), mapping$feature_index
  )
  yin_case_index <- match(yin_target[event == 1L, eid], yin$eid)
  yang_case_index <- match(yang_target$eid, yang$eid)
  if (anyNA(yin_case_index) || anyNA(yang_case_index)) stop("Dynamic protein eid alignment failed.", call. = FALSE)

  base_bins <- data.table::rbindlist(list(
    yin_target[event == 1L, .(
      eid, side = "Yin", relative_bin = baseline_bin(as.numeric(time_years))
    )],
    yang_target[, .(
      eid, side = "Yang", relative_bin = baseline_bin(-as.numeric(disease_duration_years))
    )]
  ))
  if (!anchor %in% c("baseline", "diagnosis")) {
    stop("Dynamic protein anchor must be baseline or diagnosis.", call. = FALSE)
  }
  if (anchor == "diagnosis") base_bins[, relative_bin := -relative_bin]
  base_bins[, side := if (anchor == "diagnosis") {
    data.table::fifelse(relative_bin < 0, "Yin", "Yang")
  } else {
    data.table::fifelse(relative_bin < 0, "Yang", "Yin")
  }]

  participant_rows <- vector("list", nrow(mapping))
  for (j in seq_len(nrow(mapping))) {
    yin_all <- protein_yin[, j]
    yin_median <- stats::median(yin_all[is.finite(yin_all)])
    if (!is.finite(yin_median)) stop("All Yin values are missing: ", mapping$series[[j]], call. = FALSE)
    yin_all[!is.finite(yin_all)] <- yin_median
    yang_all <- protein_yang[, j]
    yang_all[!is.finite(yang_all)] <- yin_median
    center <- mean(yin_all)
    scale <- stats::sd(yin_all)
    if (!is.finite(scale) || scale <= 0) stop("Invalid Yin scale: ", mapping$series[[j]], call. = FALSE)
    values <- data.table::rbindlist(list(
      data.table::data.table(eid = yin_target[event == 1L, eid], value = (yin_all[yin_case_index] - center) / scale),
      data.table::data.table(eid = yang_target$eid, value = (yang_all[yang_case_index] - center) / scale)
    ))
    participant_rows[[j]] <- merge(base_bins, values, by = "eid", all.x = TRUE, sort = FALSE)[,
      series := mapping$series[[j]]]
  }
  participants <- data.table::rbindlist(participant_rows)
  adjustment_spec <- yy_adjustment_spec(adjustment)
  extra_inputs <- character()

  if (adjustment_spec$id == "raw") {
    trajectories <- participants[, {
      n_value <- .N
      mean_value <- mean(value)
      se_value <- stats::sd(value) / sqrt(n_value)
      list(n = n_value, mean = mean_value, se = se_value,
           low = mean_value - 1.96 * se_value, high = mean_value + 1.96 * se_value)
    }, by = .(side, relative_bin, series)]
    trajectories[, `:=`(series_type = "Protein", bin_n = n)]
  } else if (adjustment_spec$id == "age+sex") {
    yin_covariates <- yin[, .(eid, age, sex)]
    yang_covariates <- yang[, .(eid, age, sex)]
    covariates <- data.table::rbindlist(list(yin_covariates, yang_covariates))
    participants <- merge(participants, covariates, by = "eid", all.x = TRUE, sort = FALSE)
    trajectories <- data.table::rbindlist(lapply(mapping$series, function(label) {
      age_sex_adjusted_bins(participants[series == label], label, "Protein")
    }))
  } else {
    covariates <- yy_load_stepwise_covariates(paths, adjustment_spec)
    extra_inputs <- attr(covariates, "input_path")
    participants <- merge(participants, covariates, by = "eid", all = FALSE, sort = FALSE)
    trajectories <- data.table::rbindlist(lapply(mapping$series, function(label) {
      marginal_adjusted_bins(
        participants[series == label], label, "Protein",
        adjustment = adjustment_spec$id, rhs_terms = adjustment_spec$rhs
      )
    }))
  }
  trajectories[, displayed := !relative_bin %in% range(base_bins$relative_bin)]

  roc <- metrics <- NULL
  if (isTRUE(show_roc)) {
    if (version != "unified") stop("Internal common-cohort ROC contract failed.", call. = FALSE)
    predictions <- yy_fit_single_protein_oof(yin, protein_yin, mapping$series)
    if (adjustment_spec$id != "raw") {
      if (adjustment_spec$id == "age+sex") {
        covariates <- yin[, .(eid, age, sex)]
        covariates[, `:=`(
          sex_factor = factor(sex, levels = c(0L, 1L), labels = c("Female", "Male"))
        )]
      } else {
        covariates <- yy_load_stepwise_covariates(paths, adjustment_spec)
        extra_inputs <- unique(c(extra_inputs, attr(covariates, "input_path")))
      }
      predictions <- merge(predictions, covariates, by = "eid", all = FALSE, sort = FALSE)
      predictions[, score := cross_fitted_covariate_residual(.SD, adjustment_spec$rhs), by = series]
    }
    roc_objects <- lapply(mapping$series, function(label) {
      one <- predictions[series == label, .(eid, outer_fold, time, event, score)]
      result <- mean_fold_roc(one, 5)
      result$curve[, series := label]
      list(
        curve = result$curve,
        metric = data.table::data.table(
          series = label, AUC_5y = result$mean_fold_auc,
          n = nrow(one), incident_events = sum(one$event)
        )
      )
    })
    roc <- data.table::rbindlist(lapply(roc_objects, `[[`, "curve"))
    metrics <- data.table::rbindlist(lapply(roc_objects, `[[`, "metric"))
  }

  input_paths <- unique(c(
    paths$participants_yin, paths$participants_yang, paths$features,
    paths$protein_yin, paths$protein_yang, paths$yin_target, paths$yang_target,
    paths$core, extra_inputs
  ))
  list(
    trajectories = trajectories,
    roc = roc,
    metrics = metrics,
    mapping = mapping,
    input_paths = input_paths
  )
}
