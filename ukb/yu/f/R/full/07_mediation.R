yur_mediation_dataset_file <- function(cfg) file.path(cfg$paths$mediation, "mediation_analysis_dataset.rds")

yur_mediation_mode_impute <- function(x) {
  if (is.character(x) || is.factor(x)) {
    x <- trimws(as.character(x))
    x[!nzchar(x) | toupper(x) == "NA"] <- NA_character_
  }
  replacement <- yur_mode(x)
  x[is.na(x)] <- replacement
  x
}

yur_mediation_factor <- function(x) factor(yur_mediation_mode_impute(x))

yur_mediation_multiselect <- function(x, codes, names_out) {
  value <- trimws(as.character(x))
  value[!nzchar(value) | value %in% c("-3", "NA")] <- NA_character_
  value <- yur_mediation_mode_impute(value)
  out <- as.data.table(lapply(codes, function(code) {
    as.integer(grepl(paste0("(^|\\|)", code, "(\\||$)"), value))
  }))
  setnames(out, names_out)
  out
}

yur_mediation_full_panel_n <- function(cfg, cox = NULL) {
  summary_file <- file.path(cfg$paths$cox, "cox_summary.json")
  summary_n <- NA_integer_
  if (file.exists(summary_file)) {
    summary <- read_json(summary_file, simplifyVector = TRUE)
    summary_n <- as.integer(summary$retained_protein_n)
  }
  threshold_n <- NA_integer_
  if (!is.null(cox) && "bonferroni_threshold" %in% names(cox)) {
    threshold <- unique(as.numeric(cox$bonferroni_threshold))
    threshold <- threshold[is.finite(threshold) & threshold > 0]
    if (length(threshold)) threshold_n <- as.integer(round(.05 / threshold[[1]]))
  }
  candidates <- unique(c(summary_n, threshold_n))
  candidates <- candidates[is.finite(candidates) & candidates > 0]
  if (!length(candidates) || length(candidates) != 1L) {
    stop("Full-panel Bonferroni denominator is missing or inconsistent.", call. = FALSE)
  }
  candidates[[1]]
}

yur_write_mediation_method_comparison <- function(cfg, participants, proteins, full_panel_n) {
  comparison <- data.table(
    component = c(
      "analysis cohort", "protein panel", "candidate path screening",
      "covariate set", "missing-data handling", "education/employment coding",
      "mediation estimator", "final Figure 5B/5C inclusion"
    ),
    article_method = c(
      "53,026 participants with complete proteomics and no baseline CVD",
      "2,920 proteins",
      "protein-outcome P<0.05/2920; factor-outcome P<0.05; factor-protein P<0.05/2920; direction concordant",
      "age, sex, ethnicity, education, employment, TDI, smoking, alcohol, BMI, SBP",
      "continuous median; categorical mode; TDI site median; ethnicity White; BMI/SBP sex-specific median",
      "UK Biobank education and employment categories; exact collapse not reported in the article",
      "CMAverse cmest causal mediation; linear mediator model plus Cox outcome model; bootstrap B=1000",
      "article text uses indirect-effect P<0.05/candidate count; official Table S16 is a reported subset that does not fully reconcile with the article count"
    ),
    local_implementation = c(
      sprintf("%s participants in the locally derived no-baseline-CVD cohort", format(participants, big.mark = ",")),
      sprintf(
        "%s locally Bonferroni-associated mediators selected from the %s-protein QC panel",
        format(proteins, big.mark = ","), format(full_panel_n, big.mark = ",")
      ),
      "same three association thresholds and direction-concordance screen",
      "same named covariates; smoking and alcohol retained as never/former/current factors",
      "implemented as stated, before standardization and model fitting",
      "raw p6138_i0 and p6142_i0 mode-imputed then expanded to option-level multi-hot indicators",
      "CMAverse cmest regression-based paramfunc; linear mediator model plus Cox outcome model; percentile bootstrap B=1000",
      "indirect-effect bootstrap P<0.05/local frozen candidate count; valid 0-100% mediation proportions enter Figure 5"
    ),
    status = c(
      "DIFFERENT_LOCAL_COHORT", "LOCAL_PANEL_DERIVED", "MATCHED",
      "MATCHED_NAMED_COVARIATES", "MATCHED", "PARTIAL_EXACT_COLLAPSE_UNREPORTED",
      "MATCHED_REPORTED_COMPONENTS_WITH_UNREPORTED_ASSUMPTIONS",
      "MATCHED_ARTICLE_TEXT_NUMERIC_REPLICATION_DIFFERENT"
    )
  )
  yur_write_csv(comparison, file.path(cfg$paths$mediation, "mediation_article_method_comparison.csv"))
}

yur_write_cmest_article_qc <- function(cfg, result, nboot) {
  pass <- result[status == "PASS"]
  positive_p <- sort(unique(pass$indirect_p[is.finite(pass$indirect_p) & pass$indirect_p > 0]))
  valid <- pass[
    is.finite(proportion_mediated_pct) &
      proportion_mediated_pct >= 0 & proportion_mediated_pct <= 100
  ]
  official_rows <- official_mediators <- official_outcomes <- official_exposures <- NA_integer_
  workbook <- cfg$supplement_workbook_file
  if (file.exists(workbook) && requireNamespace("readxl", quietly = TRUE)) {
    official <- as.data.table(readxl::read_excel(workbook, sheet = "S16", skip = 1))
    setnames(official, trimws(names(official)))
    official_rows <- nrow(official)
    if ("Mediator" %in% names(official)) official_mediators <- uniqueN(official$Mediator)
    if ("Outcome" %in% names(official)) official_outcomes <- uniqueN(official$Outcome)
    if ("Exposure" %in% names(official)) official_exposures <- uniqueN(official$Exposure)
  }

  article_candidates <- 6665L
  article_significant <- 4757L
  local_candidates <- nrow(result)
  local_significant <- pass[significant_bonferroni %in% TRUE, .N]
  qc <- data.table(
    item = c(
      "candidate_paths", "successful_paths", "bonferroni_significant_paths",
      "valid_0_100pct_significant_paths", "minimum_positive_bootstrap_p",
      "theoretical_bootstrap_p_resolution", "official_table_s16_rows",
      "official_table_s16_unique_mediators", "official_table_s16_outcomes",
      "official_table_s16_exposures"
    ),
    article_or_expected = c(
      article_candidates, article_candidates, article_significant, NA, NA,
      2 / nboot, official_rows, official_mediators, official_outcomes, official_exposures
    ),
    local = c(
      local_candidates, nrow(pass), local_significant,
      valid[significant_bonferroni %in% TRUE, .N],
      if (length(positive_p)) positive_p[[1]] else NA_real_,
      2 / nboot, official_rows, official_mediators, official_outcomes, official_exposures
    )
  )
  qc[, difference := local - article_or_expected]
  qc[, status := fifelse(
    item %chin% c("candidate_paths", "bonferroni_significant_paths") & difference != 0,
    "NUMERICALLY_DIFFERENT",
    fifelse(is.na(difference) | difference == 0, "MATCHED_OR_DESCRIPTIVE", "DIFFERENT")
  )]
  yur_write_csv(qc, file.path(cfg$paths$mediation, "mediation_cmest_article_qc.csv"))
  yur_write_json(list(
    status = "METHOD_MATCHED_NUMERIC_REPLICATION_DIFFERENT",
    article_candidate_paths = article_candidates,
    local_candidate_paths = local_candidates,
    article_bonferroni_significant_paths = article_significant,
    local_bonferroni_significant_paths = local_significant,
    local_valid_0_100pct_significant_paths = valid[significant_bonferroni %in% TRUE, .N],
    candidate_count_difference = local_candidates - article_candidates,
    significant_count_difference = local_significant - article_significant,
    bootstrap_replicates = nboot,
    minimum_positive_bootstrap_p = if (length(positive_p)) positive_p[[1]] else NA_real_,
    theoretical_two_sided_sign_count_resolution = 2 / nboot,
    interpretation = paste(
      "The reported CMAverse components and multiplicity rule were implemented.",
      "The local candidate count nearly matches the article, but the significant-path count does not.",
      "No threshold or model setting was changed after inspecting this difference."
    ),
    likely_sources_of_difference = c(
      "local cohort size and endpoint reconstruction",
      "local covariate coding and missing-data implementation",
      "article-unreported cmest estimation and bootstrap details",
      "official Table S16 does not numerically reconcile with the article-text count"
    ),
    article_source = "https://academic.oup.com/proteincell/article/17/3/231/8250438"
  ), file.path(cfg$paths$mediation, "mediation_cmest_article_qc.json"))

  comparison_file <- file.path(cfg$paths$mediation, "mediation_article_method_comparison.csv")
  if (file.exists(comparison_file)) {
    comparison <- fread(comparison_file)
    comparison[
      component == "mediation estimator",
      `:=`(
        local_implementation = paste(
          "CMAverse cmest regression-based paramfunc; linear mediator model plus",
          "Cox outcome model; percentile bootstrap B=1000"
        ),
        status = "MATCHED_REPORTED_COMPONENTS_WITH_UNREPORTED_ASSUMPTIONS"
      )
    ]
    comparison[
      component == "final Figure 5B/5C inclusion",
      `:=`(
        local_implementation = paste(
          "indirect-effect bootstrap P<0.05/local frozen candidate count;",
          "valid 0-100% mediation proportions enter Figure 5"
        ),
        status = "MATCHED_ARTICLE_TEXT_NUMERIC_REPLICATION_DIFFERENT"
      )
    ]
    yur_write_csv(comparison, comparison_file)
  }
}

yur_write_mediation_numeric_comparison <- function(cfg, res) {
  workbook <- cfg$supplement_workbook_file
  if (!file.exists(workbook) || !requireNamespace("readxl", quietly = TRUE)) return(invisible(NULL))

  official <- as.data.table(readxl::read_excel(workbook, sheet = "S16", skip = 1))
  setnames(official, trimws(names(official)))
  prop_col <- grep("Proportion mediated", names(official), value = TRUE, fixed = TRUE)[[1]]
  official[, proportion := suppressWarnings(as.numeric(sub(" .*", "", get(prop_col))))]

  article_outcomes <- c(
    "abdominal_aneurysm", "aortic_valve_stenosis", "atrial_fibrillation", "cad",
    "deep_vein_thrombosis", "heart_failure", "ischemic_stroke",
    "peripheral_arterial_disease", "pulmonary_embolism"
  )
  local <- res[
    !nzchar(cox_warning) & p_value < .05 & is.finite(proportion_mediated_pct) &
      proportion_mediated_pct >= 0 & proportion_mediated_pct <= 100 &
      outcome_id %chin% article_outcomes
  ]

  comparison <- data.table(
    metric = c(
      "reported_rows", "unique_mediators", "outcomes", "risk_factor_exposures",
      "minimum_mediation_proportion_pct", "maximum_mediation_proportion_pct"
    ),
    article_S16 = c(
      nrow(official), uniqueN(official$Mediator), uniqueN(official$Outcome),
      uniqueN(official$Exposure), min(official$proportion, na.rm = TRUE),
      max(official$proportion, na.rm = TRUE)
    ),
    local_reconstruction = c(
      nrow(local), uniqueN(local$protein), uniqueN(local$outcome_id),
      uniqueN(local$exposure), min(local$proportion_mediated_pct, na.rm = TRUE),
      max(local$proportion_mediated_pct, na.rm = TRUE)
    )
  )
  comparison[, difference := local_reconstruction - article_S16]
  comparison[, status := fifelse(
    metric %chin% c("outcomes", "risk_factor_exposures") & difference == 0,
    "MATCHED", "DIFFERENT"
  )]
  yur_write_csv(
    comparison,
    file.path(cfg$paths$mediation, "mediation_original_numeric_comparison.csv")
  )

  official_breadth <- official[, .(outcomes_supported = uniqueN(Outcome)), by = Mediator][
    , .(article_mediators = .N), by = outcomes_supported
  ]
  local_breadth <- local[, .(outcomes_supported = uniqueN(outcome_id)), by = protein][
    , .(local_mediators = .N), by = outcomes_supported
  ]
  breadth <- merge(
    data.table(outcomes_supported = seq_len(9L)), official_breadth,
    by = "outcomes_supported", all.x = TRUE
  )
  breadth <- merge(breadth, local_breadth, by = "outcomes_supported", all.x = TRUE)
  breadth[is.na(article_mediators), article_mediators := 0L]
  breadth[is.na(local_mediators), local_mediators := 0L]
  breadth[, difference := local_mediators - article_mediators]
  yur_write_csv(
    breadth,
    file.path(cfg$paths$mediation, "mediation_breadth_original_vs_local.csv")
  )
}

yur_read_education_employment <- function(cfg, eids) {
  cache <- file.path(cfg$paths$cache, "education_employment.csv.gz")
  if (file.exists(cache)) {
    x <- fread(cache, showProgress = FALSE)
    x[, eid := yur_norm_eid(eid)]
    return(x[eid %chin% eids])
  }
  raw <- cfg$raw_phenotype_file
  if (!file.exists(raw)) stop("Raw phenotype file is required for education/employment: ", raw, call. = FALSE)
  header <- names(fread(raw, nrows = 0, showProgress = FALSE))
  eid_col <- intersect(c("eid", "f.eid", "participant_id"), header)
  education_col <- intersect(c("p6138_i0", "p6138_i0_a0"), header)
  employment_col <- intersect(c("p6142_i0", "p6142_i0_a0"), header)
  if (!length(eid_col) || !length(education_col) || !length(employment_col)) {
    stop("Required raw fields p6138_i0 (education) and p6142_i0 (employment) are unresolved.", call. = FALSE)
  }
  x <- fread(raw, select = c(eid_col[[1]], education_col[[1]], employment_col[[1]]), showProgress = TRUE)
  setnames(x, c(eid_col[[1]], education_col[[1]], employment_col[[1]]), c("eid", "education", "employment"))
  x[, eid := yur_norm_eid(eid)]
  x <- x[eid %chin% eids]
  fwrite(x, cache, na = "")
  x
}

yur_prepare_local_mediation <- function(cfg) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("survival package is required.", call. = FALSE)
  cox_file <- file.path(cfg$paths$cox, "table_s2_incident_associations.csv.gz")
  cohort_files <- file.path(cfg$paths$cohort, c("derivation_cohort.csv.gz", "test_cohort.csv.gz"))
  if (!file.exists(cox_file) || any(!file.exists(cohort_files))) {
    stop("Local Cox and cohort outputs are required before mediation.", call. = FALSE)
  }
  cox <- fread(cox_file, showProgress = FALSE)
  full_panel_n <- yur_mediation_full_panel_n(cfg, cox)
  significant <- cox[bonferroni_significant %in% TRUE]
  features <- sort(unique(significant$feature_id))
  if (!length(features)) stop("No Bonferroni-significant local proteins for mediation.", call. = FALSE)

  cohort <- rbindlist(lapply(cohort_files, fread, showProgress = FALSE), use.names = TRUE)
  cohort[, eid := yur_norm_eid(eid)]
  all <- as.data.table(readRDS(cfg$phenotype_rds))
  if (!"eid" %in% names(all)) stop("all.rds has no eid column.", call. = FALSE)
  all[, eid := yur_norm_eid(eid)]
  available <- intersect(c("eid", "bb_TG", "bb_HBA1C", "hba1c_ngsp", "occupation"), names(all))
  cohort <- merge(cohort, all[, ..available], by = "eid", all.x = TRUE)
  ee <- yur_read_education_employment(cfg, cohort$eid)
  cohort <- merge(cohort, ee, by = "eid", all.x = TRUE)

  protein_header <- names(fread(cfg$raw_protein_file, nrows = 0, showProgress = FALSE))
  features <- intersect(features, protein_header)
  if (!length(features)) stop("None of the mediation proteins exist in raw protein data.", call. = FALSE)
  proteins <- fread(cfg$raw_protein_file, select = c("eid", features), showProgress = TRUE)
  proteins[, eid := yur_norm_eid(eid)]
  cohort <- merge(cohort, proteins, by = "eid", all.x = FALSE)
  rm(all, proteins); invisible(gc())

  cohort[, tg := as.numeric(bb_TG)]
  cohort[, hba1c := as.numeric(fifelse(is.finite(as.numeric(hba1c_ngsp)), hba1c_ngsp, bb_HBA1C))]
  cohort[, smoking_binary := as.numeric(smoking_current)]

  imputation_variables <- c(
    "ethnicity_white", "tdi", "smoking", "alcohol", "education", "employment",
    "bmi", "sbp", "tg", "hba1c", "smoking_binary"
  )
  missing_before <- data.table(
    variable = imputation_variables,
    missing_before_n = vapply(imputation_variables, function(v) {
      value <- cohort[[v]]
      if (is.character(value) || is.factor(value)) {
        sum(is.na(value) | !nzchar(trimws(as.character(value))) | toupper(trimws(as.character(value))) == "NA")
      } else {
        sum(is.na(value))
      }
    }, integer(1))
  )

  cohort[is.na(ethnicity_white), ethnicity_white := 1L]
  if (anyNA(cohort$tdi)) {
    site_medians <- cohort[, .(site_median = median(as.numeric(tdi), na.rm = TRUE)), by = assessment_center]
    site_medians[!is.finite(site_median), site_median := NA_real_]
    cohort <- site_medians[cohort, on = "assessment_center"]
    cohort[is.na(tdi), tdi := site_median]
    cohort[, site_median := NULL]
    tdi_median <- median(as.numeric(cohort$tdi), na.rm = TRUE)
    cohort[is.na(tdi), tdi := tdi_median]
  }
  for (v in c("bmi", "sbp")) {
    sex_medians <- cohort[, .(sex_median = median(as.numeric(get(v)), na.rm = TRUE)), by = sex]
    sex_medians[!is.finite(sex_median), sex_median := NA_real_]
    for (sx in sex_medians$sex) {
      value <- sex_medians[sex == sx, sex_median]
      cohort[sex == sx & is.na(get(v)), (v) := value]
    }
    overall_median <- median(as.numeric(cohort[[v]]), na.rm = TRUE)
    cohort[is.na(get(v)), (v) := overall_median]
  }
  for (v in c("tg", "hba1c")) {
    value <- median(as.numeric(cohort[[v]]), na.rm = TRUE)
    cohort[is.na(get(v)), (v) := value]
  }
  cohort[, smoking_binary := as.numeric(yur_mediation_mode_impute(smoking_binary))]
  cohort[, smoking_factor := yur_mediation_factor(smoking)]
  cohort[, alcohol_factor := yur_mediation_factor(alcohol)]
  education_columns <- c(paste0("education_", 1:6), "education_none")
  employment_columns <- c(paste0("employment_", 1:7), "employment_none")
  cohort <- cbind(
    cohort,
    yur_mediation_multiselect(cohort$education, c(as.character(1:6), "-7"), education_columns),
    yur_mediation_multiselect(cohort$employment, c(as.character(1:7), "-7"), employment_columns)
  )
  cohort[, sex_factor := factor(sex)]
  cohort[, ethnicity_factor := factor(ethnicity_white)]
  for (v in c("bmi", "sbp", "tg", "hba1c")) {
    mu <- mean(cohort[[v]], na.rm = TRUE); sig <- sd(cohort[[v]], na.rm = TRUE)
    cohort[, (v) := (get(v) - mu) / sig]
  }
  for (v in features) {
    mu <- mean(cohort[[v]], na.rm = TRUE); sig <- sd(cohort[[v]], na.rm = TRUE)
    if (is.finite(sig) && sig > 0) set(cohort, j = v, value = (cohort[[v]] - mu) / sig)
  }
  saveRDS(cohort, yur_mediation_dataset_file(cfg), compress = FALSE)
  yur_write_csv(data.table(feature_id = features), file.path(cfg$paths$mediation, "mediation_protein_panel.csv"))
  missingness <- missing_before[, .(
    variable, missing_before_n,
    missing_after_n = vapply(variable, function(v) {
      target <- switch(v,
        smoking = "smoking_factor", alcohol = "alcohol_factor",
        education = education_columns[[1]], employment = employment_columns[[1]], v
      )
      sum(is.na(cohort[[target]]))
    }, integer(1)),
    n = nrow(cohort)
  )]
  missingness[, `:=`(
    missing_before_rate = missing_before_n / n,
    missing_after_rate = missing_after_n / n
  )]
  yur_write_csv(missingness, file.path(cfg$paths$mediation, "mediation_input_missingness.csv"))
  yur_write_mediation_method_comparison(cfg, nrow(cohort), length(features), full_panel_n)
  yur_write_json(list(
    status = "PASS", participants = nrow(cohort), proteins = length(features), full_panel_n = full_panel_n,
    risk_factors = c("Body mass index", "Smoking status", "Systolic blood pressure", "Triglycerides", "HbA1c"),
    method = "local observational mediation; product of coefficients with Cox outcome model",
    covariate_imputation = "article-aligned median/mode; TDI site median; ethnicity White; BMI/SBP sex-specific median",
    smoking_alcohol_coding = "never/former/current categorical covariates",
    exact_article_cma_method = FALSE,
    exact_article_cma_status = if (requireNamespace("CMAverse", quietly = TRUE)) {
      "READY_CMAVERSE_SEPARATE_MAIN_TRACK"
    } else {
      "BLOCKED_CMAVERSE_NOT_INSTALLED"
    },
    dataset_sha256 = yur_sha_file(yur_mediation_dataset_file(cfg))
  ), file.path(cfg$paths$mediation, "mediation_prepare_summary.json"))
}

yur_fit_term <- function(formula, data, term, family = c("lm", "cox")) {
  family <- match.arg(family)
  warning_messages <- character()
  fit <- withCallingHandlers(
    tryCatch(
      if (family == "lm") stats::lm(formula, data = data, na.action = na.omit)
      else survival::coxph(formula, data = data, na.action = na.omit, ties = "efron"),
      error = function(e) NULL
    ),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  warning_text <- paste(unique(warning_messages), collapse = " | ")
  if (is.null(fit)) {
    return(data.table(beta = NA_real_, se = NA_real_, p = NA_real_, n = 0L, events = NA_integer_, warning = warning_text))
  }
  sm <- summary(fit)
  co <- if (family == "lm") sm$coefficients else sm$coefficients
  fit_n <- if (family == "cox") as.integer(fit$n) else as.integer(stats::nobs(fit))
  fit_events <- if (family == "cox") as.integer(fit$nevent) else NA_integer_
  if (!term %in% rownames(co)) {
    return(data.table(beta = NA_real_, se = NA_real_, p = NA_real_, n = fit_n, events = fit_events, warning = warning_text))
  }
  beta_col <- if (family == "lm") "Estimate" else "coef"
  se_col <- if (family == "lm") "Std. Error" else "se(coef)"
  pcol <- if (family == "lm") "Pr(>|t|)" else "Pr(>|z|)"
  if (!all(c(beta_col, se_col, pcol) %in% colnames(co))) {
    return(data.table(beta = NA_real_, se = NA_real_, p = NA_real_, n = fit_n, events = fit_events, warning = warning_text))
  }
  data.table(
    beta = unname(co[term, beta_col]), se = unname(co[term, se_col]),
    p = unname(co[term, pcol]), n = fit_n, events = fit_events, warning = warning_text
  )
}

yur_fit_cox_terms <- function(formula, data, terms) {
  fit <- tryCatch(
    suppressWarnings(survival::coxph(formula, data = data, na.action = na.omit, ties = "efron")),
    error = function(e) NULL
  )
  empty <- data.table(term = terms, beta = NA_real_, se = NA_real_, p = NA_real_, n = 0L, events = NA_integer_)
  if (is.null(fit)) return(empty)
  co <- summary(fit)$coefficients
  empty[, `:=`(n = as.integer(fit$n), events = as.integer(fit$nevent))]
  for (term in terms) {
    if (term %in% rownames(co) && all(c("coef", "se(coef)", "Pr(>|z|)") %in% colnames(co))) {
      empty[term == ..term, `:=`(
        beta = unname(co[term, "coef"]), se = unname(co[term, "se(coef)"]),
        p = unname(co[term, "Pr(>|z|)"])
      )]
    }
  }
  empty
}

yur_fast_lm_panel <- function(data, proteins, exposure, covars) {
  design <- stats::model.matrix(
    stats::reformulate(c(exposure, covars)), data = data,
    contrasts.arg = NULL
  )
  term_index <- match(exposure, colnames(design))
  if (is.na(term_index)) stop("Exposure is absent from linear-model design: ", exposure, call. = FALSE)
  design_ok <- rowSums(!is.finite(design)) == 0L
  out <- vector("list", length(proteins))
  for (i in seq_along(proteins)) {
    protein <- proteins[[i]]
    y <- as.numeric(data[[protein]])
    ok <- design_ok & is.finite(y)
    n <- sum(ok)
    row <- data.table(feature_id = protein, beta = NA_real_, se = NA_real_, p = NA_real_, n = n, events = NA_integer_)
    if (n > ncol(design) + 2L && stats::sd(y[ok]) > 0) {
      fit <- tryCatch(stats::lm.fit(design[ok, , drop = FALSE], y[ok]), error = function(e) NULL)
      if (!is.null(fit) && fit$rank > 0L) {
        pivot <- fit$qr$pivot[seq_len(fit$rank)]
        position <- match(term_index, pivot)
        if (!is.na(position) && is.finite(fit$coefficients[[term_index]])) {
          residual_df <- n - fit$rank
          r_matrix <- qr.R(fit$qr)[seq_len(fit$rank), seq_len(fit$rank), drop = FALSE]
          variance <- sum(fit$residuals^2) / residual_df
          unscaled <- tryCatch(chol2inv(r_matrix), error = function(e) NULL)
          if (!is.null(unscaled) && is.finite(unscaled[position, position])) {
            beta_value <- unname(fit$coefficients[[term_index]])
            se_value <- sqrt(variance * unscaled[position, position])
            p_value <- 2 * stats::pt(abs(beta_value / se_value), df = residual_df, lower.tail = FALSE)
            row[, `:=`(beta = beta_value, se = se_value, p = p_value)]
          }
        }
      }
    }
    out[[i]] <- row
  }
  rbindlist(out)
}

yur_fast_cox_exposure_mediator <- function(time, event, exposure, mediator, covariates) {
  design <- cbind(exposure = exposure, mediator = mediator, covariates)
  ok <- is.finite(time) & is.finite(event) & rowSums(!is.finite(design)) == 0L
  n <- sum(ok)
  events <- sum(event[ok] == 1L)
  empty <- data.table(
    term = c("exposure", "mediator"), beta = NA_real_, se = NA_real_, p = NA_real_,
    n = n, events = events, warning = ""
  )
  if (n < 100L || events < 5L) return(empty)
  warning_messages <- character()
  fit <- withCallingHandlers(
    tryCatch(
      survival::coxph.fit(
        x = design[ok, , drop = FALSE], y = survival::Surv(time[ok], event[ok]),
        strata = NULL, offset = rep(0, n), init = NULL, weights = rep(1, n),
        control = survival::coxph.control(iter.max = 20), method = "efron",
        rownames = as.character(seq_len(n))
      ),
      error = function(e) NULL
    ),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  empty[, warning := paste(unique(warning_messages), collapse = " | ")]
  if (is.null(fit) || length(fit$coefficients) < 2L) return(empty)
  for (i in seq_len(2L)) {
    beta_value <- unname(fit$coefficients[[i]])
    se_value <- if (nrow(fit$var) >= i) sqrt(fit$var[i, i]) else NA_real_
    if (is.finite(beta_value) && is.finite(se_value) && se_value > 0) {
      p_value <- 2 * stats::pnorm(abs(beta_value / se_value), lower.tail = FALSE)
      empty[i, `:=`(beta = beta_value, se = se_value, p = p_value)]
    }
  }
  empty
}

yur_run_local_mediation <- function(cfg) {
  dataset_file <- yur_mediation_dataset_file(cfg)
  cox_file <- file.path(cfg$paths$cox, "table_s2_incident_associations.csv.gz")
  panel_file <- file.path(cfg$paths$mediation, "mediation_protein_panel.csv")
  if (!file.exists(dataset_file) || !file.exists(panel_file)) stop("Run mediation_prepare first.", call. = FALSE)
  d <- as.data.table(readRDS(dataset_file))
  panel_features <- fread(panel_file)$feature_id
  cox <- fread(cox_file, showProgress = FALSE)
  full_panel_n <- yur_mediation_full_panel_n(cfg, cox)
  cox <- cox[feature_id %chin% panel_features & bonferroni_significant %in% TRUE]
  endpoints <- unique(cox$outcome_id)
  exposure_map <- c(
    bmi = "Body mass index", smoking_binary = "Smoking status",
    sbp = "Systolic blood pressure", tg = "Triglycerides", hba1c = "HbA1c"
  )
  education_columns <- c(paste0("education_", 1:6), "education_none")
  employment_columns <- c(paste0("employment_", 1:7), "employment_none")
  required_design_columns <- c(education_columns, employment_columns, "smoking_factor", "alcohol_factor")
  if (any(!required_design_columns %in% names(d))) {
    stop("Mediation dataset predates article-aligned categorical coding; rerun mediation_prepare.", call. = FALSE)
  }
  base_covars <- c(
    "age", "sex_factor", "ethnicity_factor", education_columns, employment_columns,
    "tdi", "smoking_factor", "alcohol_factor", "bmi", "sbp"
  )

  rf_outcome <- list(); rf_protein <- list(); outcome_idx <- 0L; protein_idx <- 0L
  for (exposure in names(exposure_map)) {
    covars <- setdiff(base_covars, exposure)
    if (exposure == "smoking_binary") covars <- setdiff(covars, "smoking_factor")
    rhs <- paste(c(exposure, covars), collapse = " + ")
    for (outcome in endpoints) {
      event <- paste0("event_", outcome); time <- paste0("time_", outcome)
      if (!all(c(event, time) %in% names(d))) next
      fit <- yur_fit_term(as.formula(paste0("survival::Surv(", time, ", ", event, ") ~ ", rhs)), d, exposure, "cox")
      outcome_idx <- outcome_idx + 1L
      rf_outcome[[outcome_idx]] <- cbind(data.table(exposure, exposure_label = exposure_map[[exposure]], outcome_id = outcome), fit)
    }
    protein_idx <- protein_idx + 1L
    fast_panel <- yur_fast_lm_panel(d, panel_features, exposure, covars)
    fast_panel[, `:=`(exposure = exposure, exposure_label = exposure_map[[exposure]])]
    setcolorder(fast_panel, c("exposure", "exposure_label", "feature_id", "beta", "se", "p", "n", "events"))
    rf_protein[[protein_idx]] <- fast_panel
    yur_log(cfg, "Risk-factor/protein screening complete: ", exposure, " proteins=", nrow(fast_panel))
  }
  rfo <- rbindlist(rf_outcome, fill = TRUE)
  rfp <- rbindlist(rf_protein, fill = TRUE)
  yur_write_csv(rfo, file.path(cfg$paths$mediation, "risk_factor_outcome_associations.csv"))
  yur_write_csv(rfp, file.path(cfg$paths$mediation, "risk_factor_protein_associations.csv"))

  candidates <- merge(
    cox[, .(outcome_id, outcome_label, feature_id, protein, protein_outcome_beta = beta,
            protein_outcome_p = p)],
    rfp[p < .05 / full_panel_n, .(exposure, exposure_label, feature_id,
                                   exposure_protein_beta = beta, exposure_protein_se = se,
                                   exposure_protein_p = p, exposure_protein_n = n)],
    by = "feature_id", allow.cartesian = TRUE
  )
  candidates <- merge(candidates, rfo[p < .05, .(exposure, outcome_id,
                                                  exposure_outcome_beta = beta,
                                                  exposure_outcome_p = p)],
                      by = c("exposure", "outcome_id"))
  candidates[, direction_consistent := sign(exposure_protein_beta * protein_outcome_beta) == sign(exposure_outcome_beta)]
  candidates <- candidates[direction_consistent %in% TRUE]
  yur_write_csv(candidates, file.path(cfg$paths$mediation, "mediation_candidate_triangles.csv"))
  if (!nrow(candidates)) stop("No directionally consistent mediation triangles passed screening.", call. = FALSE)

  cox_covariates <- lapply(names(exposure_map), function(exposure) {
    covars <- setdiff(base_covars, exposure)
    if (exposure == "smoking_binary") covars <- setdiff(covars, "smoking_factor")
    design <- stats::model.matrix(stats::reformulate(covars), data = d)
    design[, colnames(design) != "(Intercept)", drop = FALSE]
  })
  names(cox_covariates) <- names(exposure_map)

  result <- vector("list", nrow(candidates))
  for (i in seq_len(nrow(candidates))) {
    z <- candidates[i]; exposure <- z$exposure; mediator <- z$feature_id; outcome <- z$outcome_id
    a <- data.table(
      beta = z$exposure_protein_beta, se = z$exposure_protein_se,
      n = z$exposure_protein_n
    )
    out_terms <- yur_fast_cox_exposure_mediator(
      time = as.numeric(d[[paste0("time_", outcome)]]),
      event = as.numeric(d[[paste0("event_", outcome)]]),
      exposure = as.numeric(d[[exposure]]), mediator = as.numeric(d[[mediator]]),
      covariates = cox_covariates[[exposure]]
    )
    direct <- out_terms[term == "exposure"]
    b <- out_terms[term == "mediator"]
    indirect <- a$beta * b$beta
    indirect_se <- sqrt((b$beta^2 * a$se^2) + (a$beta^2 * b$se^2))
    p_indirect <- 2 * stats::pnorm(-abs(indirect / indirect_se))
    total <- direct$beta + indirect
    proportion <- 100 * indirect / total
    result[[i]] <- data.table(
      exposure, exposure_label = z$exposure_label, mediator, protein = z$protein,
      outcome_id = outcome, outcome_label = z$outcome_label,
      n = min(a$n, b$n), events = b$events,
      path_a = a$beta, path_a_se = a$se, path_b = b$beta, path_b_se = b$se,
      direct_effect = direct$beta, indirect_effect = indirect, indirect_se = indirect_se,
      total_effect = total, proportion_mediated_pct = proportion,
      p_value = p_indirect, cox_warning = b$warning,
      method = "product_of_coefficients_delta"
    )
    if (i %% 100L == 0L || i == nrow(candidates)) {
      yur_log(cfg, "Mediation progress ", i, "/", nrow(candidates))
      fwrite(rbindlist(result[seq_len(i)], fill = TRUE),
             file.path(cfg$paths$mediation, "mediation_results.partial.csv"), na = "")
    }
  }
  res <- rbindlist(result, fill = TRUE)
  res[, candidate_threshold := .05 / nrow(candidates)]
  res[, significant_nominal := p_value < .05]
  res[, significant_bonferroni := p_value < candidate_threshold]
  res[, fdr := p.adjust(p_value, method = "BH")]
  res[, significant_fdr := fdr < .05]
  res <- res[is.finite(proportion_mediated_pct)]
  yur_write_csv(res, file.path(cfg$paths$mediation, "mediation_results.csv"))
  warning_paths <- res[nzchar(cox_warning)]
  yur_write_csv(
    warning_paths,
    file.path(cfg$paths$mediation, "mediation_cox_warning_paths.csv")
  )
  yur_write_mediation_numeric_comparison(cfg, res)
  partial <- file.path(cfg$paths$mediation, "mediation_results.partial.csv")
  if (file.exists(partial)) unlink(partial)
  valid <- res[!nzchar(cox_warning)]
  yur_write_json(list(
    status = "PASS", candidate_triangles = nrow(candidates), completed_triangles = nrow(res),
    full_panel_bonferroni_n = full_panel_n,
    risk_factor_protein_threshold = format(.05 / full_panel_n, scientific = TRUE, digits = 8),
    valid_model_paths = nrow(valid),
    nominal_significant = res[significant_nominal %in% TRUE, .N],
    nominal_significant_valid = valid[significant_nominal %in% TRUE, .N],
    fdr_significant = res[significant_fdr %in% TRUE, .N],
    fdr_significant_valid = valid[significant_fdr %in% TRUE, .N],
    bonferroni_significant = res[significant_bonferroni %in% TRUE, .N],
    bonferroni_significant_valid = valid[significant_bonferroni %in% TRUE, .N],
    cox_warning_paths = res[nzchar(cox_warning), .N],
    proteins = uniqueN(res$mediator), outcomes = uniqueN(res$outcome_id),
    method = "product-of-coefficients with linear mediator and Cox outcome; delta-method inference",
    exact_article_cma_method = FALSE,
    exact_article_cma_status = if (requireNamespace("CMAverse", quietly = TRUE)) {
      "READY_CMAVERSE_SEPARATE_MAIN_TRACK"
    } else {
      "BLOCKED_CMAVERSE_NOT_INSTALLED"
    },
    figure5_inclusion_rule = "nominal mediation P<0.05, matching official Table S16",
    limitation = "Fast audit estimator only; the CMAverse bootstrap main track is stored separately."
  ), file.path(cfg$paths$mediation, "mediation_summary.json"))
}

yur_cmest_shard_dir <- function(cfg) {
  path <- file.path(cfg$paths$mediation, "cmest_shards")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

yur_cmest_shard_file <- function(cfg, index = cfg$cmest_shard_index, count = cfg$cmest_shard_count) {
  file.path(yur_cmest_shard_dir(cfg), sprintf("cmest_shard_%03d_of_%03d.csv", index, count))
}

yur_cmest_contract_file <- function(cfg, index = cfg$cmest_shard_index, count = cfg$cmest_shard_count) {
  file.path(yur_cmest_shard_dir(cfg), sprintf("cmest_shard_%03d_of_%03d.contract.json", index, count))
}

yur_cmest_base_covars <- function(exposure) {
  covars <- c(
    "age", "sex_factor", "ethnicity_factor", paste0("education_", 1:6), "education_none",
    paste0("employment_", 1:7), "employment_none", "tdi", "smoking_factor",
    "alcohol_factor", "bmi", "sbp"
  )
  covars <- setdiff(covars, exposure)
  if (identical(exposure, "smoking_binary")) covars <- setdiff(covars, "smoking_factor")
  covars
}

yur_cmest_requirements <- function(cfg) {
  if (!requireNamespace("CMAverse", quietly = TRUE)) {
    stop("CMAverse is required. Run tools/install_cmaverse.R first.", call. = FALSE)
  }
  dataset_file <- yur_mediation_dataset_file(cfg)
  candidate_file <- file.path(cfg$paths$mediation, "mediation_candidate_triangles.csv")
  if (!file.exists(dataset_file) || !file.exists(candidate_file)) {
    stop("Run mediation_prepare and mediation_run before CMAverse.", call. = FALSE)
  }
  candidates <- fread(candidate_file, showProgress = FALSE)
  candidates[, row_id := seq_len(.N)]
  list(
    dataset_file = dataset_file,
    candidate_file = candidate_file,
    candidates = candidates,
    dataset_sha256 = yur_sha_file(dataset_file),
    candidate_sha256 = yur_sha_file(candidate_file),
    cmaverse_version = as.character(utils::packageVersion("CMAverse"))
  )
}

yur_cmest_effect <- function(effect, fit) {
  value <- function(x) {
    if (is.null(x) || !effect %in% names(x)) return(NA_real_)
    as.numeric(x[[effect]])
  }
  data.table(
    effect = effect,
    estimate = value(fit$effect.pe),
    std_error = value(fit$effect.se),
    ci_low = value(fit$effect.ci.low),
    ci_high = value(fit$effect.ci.high),
    p_value = value(fit$effect.pval)
  )
}

yur_cmest_fit_path <- function(cfg, d, candidate, nboot) {
  exposure <- as.character(candidate$exposure)
  mediator <- as.character(candidate$feature_id)
  outcome <- as.character(candidate$outcome_id)
  time_col <- paste0("time_", outcome)
  event_col <- paste0("event_", outcome)
  covars <- yur_cmest_base_covars(exposure)
  required <- unique(c(time_col, event_col, exposure, mediator, covars))
  missing_columns <- setdiff(required, names(d))
  if (length(missing_columns)) {
    stop("Missing CMAverse columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  path_data <- as.data.frame(d[, ..required])
  path_data <- path_data[stats::complete.cases(path_data), , drop = FALSE]
  n <- nrow(path_data)
  events <- sum(path_data[[event_col]] == 1L)
  if (n < 100L || events < 5L) stop("Insufficient complete observations/events.", call. = FALSE)
  if (stats::sd(path_data[[exposure]]) <= 0 || stats::sd(path_data[[mediator]]) <= 0) {
    stop("Exposure or mediator has zero variance.", call. = FALSE)
  }

  seed <- as.integer((cfg$bootstrap_seed + candidate$row_id * 104729L) %% .Machine$integer.max)
  set.seed(seed)
  warnings <- character()
  started <- Sys.time()
  invisible(utils::capture.output(
    fit <- withCallingHandlers(
      CMAverse::cmest(
        data = path_data,
        model = "rb",
        full = TRUE,
        estimation = "paramfunc",
        inference = "bootstrap",
        outcome = time_col,
        event = event_col,
        exposure = exposure,
        mediator = mediator,
        EMint = FALSE,
        basec = covars,
        yreg = "coxph",
        mreg = list("linear"),
        astar = 0,
        a = 1,
        mval = list(0),
        nboot = as.integer(nboot),
        boot.ci.type = "per"
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    type = "output"
  ))
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  effects <- rbindlist(lapply(c("Rpnde", "Rtnie", "Rte", "pm"), yur_cmest_effect, fit = fit))
  effect_value <- function(effect_name, column) {
    value <- effects[effects[["effect"]] == effect_name, get(column)]
    if (length(value)) as.numeric(value[[1]]) else NA_real_
  }
  data.table(
    row_id = as.integer(candidate$row_id),
    exposure = exposure,
    exposure_label = as.character(candidate$exposure_label),
    mediator = mediator,
    protein = as.character(candidate$protein),
    outcome_id = outcome,
    outcome_label = as.character(candidate$outcome_label),
    n = n,
    events = events,
    direct_hr = effect_value("Rpnde", "estimate"),
    direct_ci_low = effect_value("Rpnde", "ci_low"),
    direct_ci_high = effect_value("Rpnde", "ci_high"),
    direct_p = effect_value("Rpnde", "p_value"),
    indirect_hr = effect_value("Rtnie", "estimate"),
    indirect_ci_low = effect_value("Rtnie", "ci_low"),
    indirect_ci_high = effect_value("Rtnie", "ci_high"),
    indirect_p = effect_value("Rtnie", "p_value"),
    total_hr = effect_value("Rte", "estimate"),
    total_ci_low = effect_value("Rte", "ci_low"),
    total_ci_high = effect_value("Rte", "ci_high"),
    total_p = effect_value("Rte", "p_value"),
    proportion_mediated_pct = 100 * effect_value("pm", "estimate"),
    proportion_ci_low_pct = 100 * effect_value("pm", "ci_low"),
    proportion_ci_high_pct = 100 * effect_value("pm", "ci_high"),
    proportion_p = effect_value("pm", "p_value"),
    nboot = as.integer(nboot),
    seed = seed,
    elapsed_seconds = elapsed,
    warning = paste(unique(warnings), collapse = " | "),
    error = "",
    status = "PASS",
    method = "CMAverse_cmest_rb_paramfunc_coxph_linear_bootstrap_percentile"
  )
}

yur_cmest_error_row <- function(candidate, nboot, message, elapsed = NA_real_) {
  data.table(
    row_id = as.integer(candidate$row_id), exposure = as.character(candidate$exposure),
    exposure_label = as.character(candidate$exposure_label),
    mediator = as.character(candidate$feature_id), protein = as.character(candidate$protein),
    outcome_id = as.character(candidate$outcome_id), outcome_label = as.character(candidate$outcome_label),
    n = NA_integer_, events = NA_integer_, direct_hr = NA_real_, direct_ci_low = NA_real_,
    direct_ci_high = NA_real_, direct_p = NA_real_, indirect_hr = NA_real_,
    indirect_ci_low = NA_real_, indirect_ci_high = NA_real_, indirect_p = NA_real_,
    total_hr = NA_real_, total_ci_low = NA_real_, total_ci_high = NA_real_, total_p = NA_real_,
    proportion_mediated_pct = NA_real_, proportion_ci_low_pct = NA_real_,
    proportion_ci_high_pct = NA_real_, proportion_p = NA_real_, nboot = as.integer(nboot),
    seed = NA_integer_, elapsed_seconds = elapsed, warning = "", error = as.character(message),
    status = "ERROR", method = "CMAverse_cmest_rb_paramfunc_coxph_linear_bootstrap_percentile"
  )
}

yur_run_cmest_pilot <- function(cfg) {
  req <- yur_cmest_requirements(cfg)
  candidate <- req$candidates[
    exposure == "bmi" & toupper(protein) == "GDF15" & outcome_id == "cad"
  ][1]
  if (!nrow(candidate)) candidate <- req$candidates[1]
  d <- as.data.table(readRDS(req$dataset_file))
  started <- Sys.time()
  result <- tryCatch(
    yur_cmest_fit_path(cfg, d, candidate, cfg$cmest_pilot_nboot),
    error = function(e) yur_cmest_error_row(
      candidate, cfg$cmest_pilot_nboot, conditionMessage(e),
      as.numeric(difftime(Sys.time(), started, units = "secs"))
    )
  )
  output <- file.path(cfg$paths$mediation, "mediation_cmest_pilot.csv")
  yur_write_csv(result, output)
  seconds_per_boot <- result$elapsed_seconds / max(1L, cfg$cmest_pilot_nboot)
  projected_serial_hours <- seconds_per_boot * cfg$bootstrap_n * nrow(req$candidates) / 3600
  yur_write_json(list(
    status = result$status[[1]], row_id = result$row_id[[1]], path = paste(
      result$exposure[[1]], result$protein[[1]], result$outcome_id[[1]], sep = " -> "
    ),
    pilot_nboot = cfg$cmest_pilot_nboot,
    elapsed_seconds = result$elapsed_seconds[[1]],
    seconds_per_bootstrap = seconds_per_boot,
    full_candidate_n = nrow(req$candidates),
    full_nboot = cfg$bootstrap_n,
    projected_serial_hours = projected_serial_hours,
    projected_hours_at_8_jobs = projected_serial_hours / 8,
    dataset_sha256 = req$dataset_sha256,
    candidate_sha256 = req$candidate_sha256,
    cmaverse_version = req$cmaverse_version,
    estimation = "paramfunc",
    inference = "bootstrap",
    boot_ci_type = "percentile"
  ), file.path(cfg$paths$mediation, "mediation_cmest_pilot_summary.json"))
  if (result$status[[1]] != "PASS") stop("CMAverse pilot failed: ", result$error[[1]], call. = FALSE)
}

yur_cmest_shard_complete <- function(cfg) {
  output <- yur_cmest_shard_file(cfg)
  contract <- yur_cmest_contract_file(cfg)
  if (!file.exists(output) || !file.exists(contract)) return(FALSE)
  req <- yur_cmest_requirements(cfg)
  expected_ids <- req$candidates[((row_id - 1L) %% cfg$cmest_shard_count) + 1L == cfg$cmest_shard_index, row_id]
  observed <- tryCatch(fread(output, showProgress = FALSE)$row_id, error = function(e) integer())
  meta <- tryCatch(read_json(contract, simplifyVector = TRUE), error = function(e) NULL)
  identical(sort(as.integer(observed)), sort(as.integer(expected_ids))) &&
    !is.null(meta) && identical(as.character(meta$dataset_sha256), req$dataset_sha256) &&
    identical(as.character(meta$candidate_sha256), req$candidate_sha256) &&
    identical(as.integer(meta$nboot), as.integer(cfg$bootstrap_n))
}

yur_run_cmest_shard <- function(cfg) {
  if (cfg$cmest_shard_count < 1L || cfg$cmest_shard_index < 1L ||
      cfg$cmest_shard_index > cfg$cmest_shard_count) {
    stop("Invalid CMAverse shard index/count.", call. = FALSE)
  }
  req <- yur_cmest_requirements(cfg)
  candidates <- req$candidates[
    ((row_id - 1L) %% cfg$cmest_shard_count) + 1L == cfg$cmest_shard_index
  ]
  output <- yur_cmest_shard_file(cfg)
  completed <- integer()
  result <- list()
  if (cfg$resume && file.exists(output)) {
    previous <- fread(output, showProgress = FALSE)
    result <- split(previous, seq_len(nrow(previous)))
    completed <- as.integer(previous$row_id)
  } else if (file.exists(output)) {
    unlink(output)
  }
  d <- as.data.table(readRDS(req$dataset_file))
  pending <- candidates[!row_id %in% completed]
  yur_log(
    cfg, "CMAverse shard ", cfg$cmest_shard_index, "/", cfg$cmest_shard_count,
    " paths=", nrow(candidates), " pending=", nrow(pending), " nboot=", cfg$bootstrap_n
  )
  for (i in seq_len(nrow(pending))) {
    candidate <- pending[i]
    started <- Sys.time()
    row <- tryCatch(
      yur_cmest_fit_path(cfg, d, candidate, cfg$bootstrap_n),
      error = function(e) yur_cmest_error_row(
        candidate, cfg$bootstrap_n, conditionMessage(e),
        as.numeric(difftime(Sys.time(), started, units = "secs"))
      )
    )
    result[[length(result) + 1L]] <- row
    combined <- rbindlist(result, fill = TRUE)
    setorder(combined, row_id)
    fwrite(combined, output, na = "")
    yur_log(
      cfg, "CMAverse shard ", cfg$cmest_shard_index, " progress ", i, "/", nrow(pending),
      " row_id=", candidate$row_id, " status=", row$status,
      " seconds=", round(row$elapsed_seconds, 1)
    )
  }
  final <- fread(output, showProgress = FALSE)
  expected_ids <- sort(candidates$row_id)
  if (!identical(sort(as.integer(final$row_id)), as.integer(expected_ids))) {
    stop("CMAverse shard row-ID contract failed.", call. = FALSE)
  }
  yur_write_json(list(
    status = "PASS", shard_index = cfg$cmest_shard_index, shard_count = cfg$cmest_shard_count,
    candidate_rows = nrow(candidates), pass_rows = final[status == "PASS", .N],
    error_rows = final[status == "ERROR", .N], nboot = cfg$bootstrap_n,
    dataset_sha256 = req$dataset_sha256, candidate_sha256 = req$candidate_sha256,
    cmaverse_version = req$cmaverse_version,
    method = "CMAverse cmest regression-based paramfunc; linear mediator; Cox outcome; percentile bootstrap"
  ), yur_cmest_contract_file(cfg))
}

yur_merge_cmest_shards <- function(cfg) {
  req <- yur_cmest_requirements(cfg)
  parts <- vector("list", cfg$cmest_shard_count)
  contracts <- vector("list", cfg$cmest_shard_count)
  for (index in seq_len(cfg$cmest_shard_count)) {
    output <- yur_cmest_shard_file(cfg, index, cfg$cmest_shard_count)
    contract <- yur_cmest_contract_file(cfg, index, cfg$cmest_shard_count)
    if (!file.exists(output) || !file.exists(contract)) stop("Missing CMAverse shard: ", index, call. = FALSE)
    meta <- read_json(contract, simplifyVector = TRUE)
    if (!identical(as.character(meta$dataset_sha256), req$dataset_sha256) ||
        !identical(as.character(meta$candidate_sha256), req$candidate_sha256) ||
        !identical(as.integer(meta$nboot), as.integer(cfg$bootstrap_n))) {
      stop("Stale or incompatible CMAverse shard contract: ", index, call. = FALSE)
    }
    parts[[index]] <- fread(output, showProgress = FALSE)
    contracts[[index]] <- meta
  }
  result <- rbindlist(parts, fill = TRUE)
  setorder(result, row_id)
  if (anyDuplicated(result$row_id) ||
      !identical(as.integer(result$row_id), as.integer(req$candidates$row_id))) {
    stop("Merged CMAverse row-ID coverage does not match frozen candidates.", call. = FALSE)
  }
  result[, candidate_threshold := .05 / .N]
  result[, significant_nominal := status == "PASS" & indirect_p < .05]
  result[, significant_bonferroni := status == "PASS" & indirect_p < candidate_threshold]
  result[, fdr := NA_real_]
  result[status == "PASS" & is.finite(indirect_p), fdr := p.adjust(indirect_p, method = "BH")]
  result[, significant_fdr := status == "PASS" & fdr < .05]
  output <- file.path(cfg$paths$mediation, "mediation_cmest_results.csv")
  yur_write_csv(result, output)
  yur_write_cmest_article_qc(cfg, result, cfg$bootstrap_n)
  pass <- result[status == "PASS"]
  yur_write_json(list(
    status = if (nrow(pass) == nrow(result)) "PASS" else "PASS_WITH_PATH_ERRORS",
    candidate_triangles = nrow(req$candidates), completed_triangles = nrow(result),
    successful_triangles = nrow(pass), error_triangles = result[status == "ERROR", .N],
    nominal_significant = pass[indirect_p < .05, .N],
    bonferroni_threshold = .05 / nrow(result),
    bonferroni_significant = pass[indirect_p < .05 / nrow(result), .N],
    fdr_significant = pass[fdr < .05, .N],
    nboot = cfg$bootstrap_n, shard_count = cfg$cmest_shard_count,
    dataset_sha256 = req$dataset_sha256, candidate_sha256 = req$candidate_sha256,
    output_sha256 = yur_sha_file(output), cmaverse_version = req$cmaverse_version,
    article_candidate_count = 6665,
    local_candidate_count = nrow(req$candidates),
    article_method_match = "CMAverse v0.1.0 cmest; regression-based; linear mediator; Cox outcome; bootstrap",
    explicit_assumptions = c(
      "paramfunc selected because the paper did not report paramfunc versus imputation",
      "percentile bootstrap CI selected because it is CMAverse default and the paper did not report CI type"
    )
  ), file.path(cfg$paths$mediation, "mediation_cmest_summary.json"))
  yur_write_json(list(
    package = "CMAverse", version = req$cmaverse_version,
    model = "rb", estimation = "paramfunc", inference = "bootstrap",
    outcome_model = "coxph", mediator_model = "linear", exposure_mediator_interaction = FALSE,
    astar = 0, a = 1, mediator_control_value = 0,
    nboot = cfg$bootstrap_n, bootstrap_ci = "percentile",
    primary_mediation_test = "Rtnie bootstrap P value",
    proportion_mediated = "100 * pm",
    multiplicity = "Bonferroni 0.05 divided by local frozen candidate count",
    source_article_reported = "CMAverse 0.1.0; cmest; linear mediator; Cox outcome; bootstrap 1000",
    source_article_unreported = c("paramfunc versus imputation", "percentile versus BCa CI")
  ), file.path(cfg$paths$mediation, "mediation_cmest_method_manifest.json"))
}
