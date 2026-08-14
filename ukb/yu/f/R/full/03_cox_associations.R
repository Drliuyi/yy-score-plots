yur_impute_association_covariates <- function(meta, lag_variable = "protein_sampling_lag_days") {
  out <- copy(meta)
  if (!lag_variable %in% names(out)) lag_variable <- "protein_sampling_lag_days"
  if (!lag_variable %in% names(out)) stop("No protein sampling-lag covariate is available.", call. = FALSE)
  if (lag_variable != "protein_sampling_lag_days") out[, protein_sampling_lag_days := get(lag_variable)]
  numeric_covariates <- intersect(
    c("age", "protein_sampling_lag_days", "fasting_time", "bmi", "sbp", "tdi"),
    names(out)
  )
  for (v in numeric_covariates) set(out, j = v, value = as.numeric(out[[v]]))
  continuous <- c("age", "protein_sampling_lag_days", "fasting_time")
  for (v in continuous) {
    value <- median(out[[v]], na.rm = TRUE)
    out[is.na(get(v)), (v) := value]
  }
  for (v in c("bmi", "sbp")) {
    values <- out[, .(value = as.numeric(median(get(v), na.rm = TRUE))), by = sex]
    for (sx in values$sex) out[sex == sx & is.na(get(v)), (v) := values[sex == sx, value]]
  }
  if (anyNA(out$tdi)) {
    site <- out[, .(value = as.numeric(median(tdi, na.rm = TRUE))), by = assessment_center]
    for (center in site$assessment_center) {
      out[assessment_center == center & is.na(tdi), tdi := site[assessment_center == center, value]]
    }
    out[is.na(tdi), tdi := median(out$tdi, na.rm = TRUE)]
  }
  out[is.na(ethnicity_white), ethnicity_white := 1L]
  for (v in c("smoking", "alcohol", "blood_collection_season")) {
    out[, (v) := as.character(get(v))]
    replacement <- yur_mode(out[[v]])
    out[is.na(get(v)) | !nzchar(as.character(get(v))), (v) := replacement]
    out[, (v) := as.factor(get(v))]
  }
  out
}

yur_covariate_matrix <- function(meta, protein_panel = NULL) {
  panel_key <- yur_panel_key(protein_panel %||% "")
  requested_lag <- if (nzchar(panel_key)) paste0("protein_sampling_lag_days__", panel_key) else "protein_sampling_lag_days"
  lag_variable <- if (requested_lag %in% names(meta)) requested_lag else "protein_sampling_lag_days"
  x <- yur_impute_association_covariates(meta, lag_variable = lag_variable)
  mm <- model.matrix(
    ~ age + sex + ethnicity_white + tdi + blood_collection_season +
      protein_sampling_lag_days + fasting_time + sbp + bmi + smoking + alcohol,
    data = x
  )
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  attr(mm, "lag_variable") <- lag_variable
  mm
}

yur_cox_one <- function(protein, time, event, covariates) {
  cov_ok <- rowSums(!is.finite(covariates)) == 0L
  ok <- is.finite(protein) & is.finite(time) & is.finite(event) & cov_ok
  n <- sum(ok)
  events <- sum(event[ok] == 1L)
  if (n < 100L || events < 5L || stats::sd(protein[ok]) <= 0) {
    return(c(n = n, events = events, beta = NA, se = NA, z = NA, p = NA, hr = NA, ci_low = NA, ci_high = NA))
  }
  fit <- tryCatch({
    survival::coxph.fit(
      x = cbind(protein = protein[ok], covariates[ok, , drop = FALSE]),
      y = survival::Surv(time[ok], event[ok]),
      strata = NULL, offset = rep(0, n), init = NULL,
      weights = rep(1, n),
      control = survival::coxph.control(iter.max = 20), method = "efron",
      rownames = as.character(seq_len(n))
    )
  }, error = function(e) NULL)
  if (is.null(fit) || !length(fit$coefficients) || !is.finite(fit$coefficients[[1]])) {
    return(c(n = n, events = events, beta = NA, se = NA, z = NA, p = NA, hr = NA, ci_low = NA, ci_high = NA))
  }
  beta <- unname(fit$coefficients[[1]])
  se <- sqrt(fit$var[1, 1])
  z <- beta / se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  c(n = n, events = events, beta = beta, se = se, z = z, p = p,
    hr = exp(beta), ci_low = exp(beta - 1.96 * se), ci_high = exp(beta + 1.96 * se))
}

yur_cox_result_row <- function(feature_id, protein_mean, protein_sd, lag_variable,
                               lag_panel_fallback, stat) {
  stat_table <- as.data.table(as.list(stat))
  if (nrow(stat_table) != 1L) stop("Cox statistic assembly must produce exactly one row.", call. = FALSE)
  cbind(
    data.table(
      feature_id = feature_id, exposure_scale = "per_scope_SD",
      protein_mean = protein_mean, protein_sd = protein_sd,
      lag_covariate = lag_variable, lag_panel_fallback = lag_panel_fallback
    ),
    stat_table
  )
}

yur_panel_mapping <- function(features, cfg) {
  mapping <- fread(cfg$panel_mapping_file)
  candidates <- intersect(c("Assay", "assay", "HGNC.symbol", "gene_symbol", "symbol", "OlinkID"), names(mapping))
  if (!length(candidates)) {
    return(data.table(feature_id = features, protein = toupper(features), panel = NA_character_, olink_id = NA_character_, mapping_status = "UNMAPPED"))
  }
  long <- rbindlist(lapply(candidates, function(v) {
    data.table(norm = yur_norm_name(mapping[[v]]), row_id = seq_len(nrow(mapping)))
  }))
  raw <- data.table(feature_id = features, norm = yur_norm_name(features))
  hit <- merge(raw, long, by = "norm", all.x = TRUE, allow.cartesian = TRUE)
  hit <- hit[order(feature_id, row_id)][, .SD[1], by = feature_id]
  symbol_col <- intersect(c("HGNC.symbol", "gene_symbol", "symbol", "Assay", "assay"), names(mapping))
  panel_col <- intersect(c("Panel", "panel"), names(mapping))
  olink_col <- intersect(c("OlinkID", "olink_id"), names(mapping))
  hit[, protein := if (length(symbol_col)) toupper(as.character(mapping[[symbol_col[[1]]]][row_id])) else toupper(feature_id)]
  hit[is.na(protein) | !nzchar(protein), protein := toupper(feature_id)]
  hit[, panel := if (length(panel_col)) as.character(mapping[[panel_col[[1]]]][row_id]) else NA_character_]
  hit[, olink_id := if (length(olink_col)) as.character(mapping[[olink_col[[1]]]][row_id]) else NA_character_]
  hit[, .(feature_id, protein, panel, olink_id, mapping_status = fifelse(is.na(row_id), "UNMAPPED", "MAPPED"))]
}

yur_cox_contract_values <- function(cfg, outcome_row, panel_hash = NULL) {
  cohort_hash_file <- file.path(cfg$paths$cohort, "cohort_contract_hash.txt")
  definition_hash_file <- file.path(cfg$paths$cohort, "outcome_definition_hash.txt")
  if (!file.exists(cohort_hash_file) || !file.exists(definition_hash_file)) {
    stop("Cohort contract files are missing; rerun cohort with the current endpoint dictionary.")
  }
  if (is.null(panel_hash)) {
    panel_hash_file <- file.path(cfg$paths$cox, "retained_panel_hash.txt")
    if (!file.exists(panel_hash_file)) stop("Retained panel hash is missing; rerun cox_prepare.")
    panel_hash <- readLines(panel_hash_file, warn = FALSE)[[1]]
  }
  list(
    outcome_id = outcome_row$outcome_id[[1]],
    outcome_definition_hash = outcome_row$definition_hash[[1]],
    all_outcome_definitions_hash = readLines(definition_hash_file, warn = FALSE)[[1]],
    cohort_contract_hash = readLines(cohort_hash_file, warn = FALSE)[[1]],
    retained_panel_hash = panel_hash
  )
}

yur_cox_contract_matches <- function(path, expected) {
  if (!file.exists(path)) return(FALSE)
  observed <- tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(observed)) return(FALSE)
  keys <- names(expected)
  all(vapply(keys, function(k) identical(as.character(observed[[k]]), as.character(expected[[k]])), logical(1)))
}

yur_resolved_outcomes_for_cox <- function(cfg) {
  path <- file.path(cfg$paths$preflight, "outcome_field_resolution.csv")
  if (!file.exists(path)) stop("Resolved outcome dictionary is missing; rerun preflight.")
  x <- fread(path)
  requested <- trimws(strsplit(as.character(cfg$endpoint_subset), ",", fixed = TRUE)[[1]])
  if (!(length(requested) == 1L && tolower(requested) == "all")) x <- x[outcome_id %in% requested]
  if (!nrow(x)) stop("No resolved outcomes selected for Cox analysis.")
  x
}

yur_cox_shard_complete <- function(cfg) {
  fail <- function(reason) {
    yur_log(cfg, "Cox shard resume validation: ", reason, level = "WARN")
    FALSE
  }
  retained_file <- file.path(cfg$paths$cox, "retained_panel.csv")
  if (!file.exists(retained_file)) return(fail("retained_panel.csv is missing"))
  retained <- tryCatch(fread(retained_file, select = "feature_id"), error = function(e) NULL)
  if (is.null(retained) || !nrow(retained) || anyDuplicated(retained$feature_id)) {
    return(fail("retained panel is unreadable, empty, or duplicated"))
  }
  outcomes <- tryCatch(yur_resolved_outcomes_for_cox(cfg), error = function(e) NULL)
  if (is.null(outcomes) || !nrow(outcomes)) return(fail("resolved outcome contract is unavailable"))
  for (i in seq_len(nrow(outcomes))) {
    outcome_id <- outcomes$outcome_id[[i]]
    expected_contract <- tryCatch(yur_cox_contract_values(cfg, outcomes[i]), error = function(e) NULL)
    if (is.null(expected_contract)) return(fail(paste0("contract inputs are missing for ", outcome_id)))
    for (scope in c("full_incident", "derivation")) {
      out_file <- file.path(cfg$paths$cox, paste0(scope, "_", outcome_id, "_cox.csv.gz"))
      contract_file <- file.path(cfg$paths$cox, paste0(scope, "_", outcome_id, "_cox.contract.json"))
      if (!file.exists(out_file)) return(fail(paste0("missing ", basename(out_file))))
      if (!yur_cox_contract_matches(contract_file, expected_contract)) {
        return(fail(paste0("missing/stale contract for ", basename(out_file))))
      }
      observed <- tryCatch(fread(out_file, select = "feature_id"), error = function(e) NULL)
      if (is.null(observed)) return(fail(paste0("unreadable ", basename(out_file))))
      if (nrow(observed) != nrow(retained)) {
        return(fail(paste0("row-count mismatch in ", basename(out_file))))
      }
      if (!identical(sort(observed$feature_id), sort(retained$feature_id))) {
        return(fail(paste0("feature-set mismatch in ", basename(out_file))))
      }
    }
  }
  TRUE
}

yur_run_cox_scope <- function(cfg, meta, proteins, mapping, outcomes, scope) {
  mapping <- copy(mapping)
  mapping[, panel_key := yur_panel_key(panel)]
  mapping[is.na(panel_key), panel_key := ""]
  panel_keys <- unique(mapping$panel_key)
  covariates_by_panel <- setNames(lapply(panel_keys, function(key) {
    yur_covariate_matrix(meta, protein_panel = key)
  }), panel_keys)
  rows <- vector("list", nrow(outcomes))
  for (i in seq_len(nrow(outcomes))) {
    outcome_id <- outcomes$outcome_id[[i]]
    out_file <- file.path(cfg$paths$cox, paste0(scope, "_", outcome_id, "_cox.csv.gz"))
    contract_file <- file.path(cfg$paths$cox, paste0(scope, "_", outcome_id, "_cox.contract.json"))
    expected_contract <- yur_cox_contract_values(cfg, outcomes[i])
    if (cfg$resume && file.exists(out_file) && !cfg$force &&
        yur_cox_contract_matches(contract_file, expected_contract)) {
      yur_log(cfg, "Resume: reuse ", basename(out_file))
      rows[[i]] <- fread(out_file)
      next
    }
    if (cfg$resume && file.exists(out_file) && !cfg$force) {
      yur_log(cfg, "Resume: stale contract; recompute ", basename(out_file), level = "WARN")
    }
    time <- meta[[paste0("time_", outcome_id)]]
    event <- meta[[paste0("event_", outcome_id)]]
    yur_log(cfg, "Cox scope=", scope, " outcome=", outcome_id,
            " n=", nrow(meta), " events=", sum(event), " proteins=", ncol(proteins))
    result <- vector("list", ncol(proteins))
    for (j in seq_len(ncol(proteins))) {
      panel_hit <- mapping[feature_id == names(proteins)[[j]], panel_key]
      panel_key <- if (length(panel_hit)) panel_hit[[1]] %||% "" else ""
      covariates <- covariates_by_panel[[panel_key]]
      lag_variable <- attr(covariates, "lag_variable")
      raw_value <- proteins[[j]]
      protein_mean <- mean(raw_value, na.rm = TRUE)
      protein_sd <- stats::sd(raw_value, na.rm = TRUE)
      scaled_value <- if (is.finite(protein_sd) && protein_sd > 0) {
        (raw_value - protein_mean) / protein_sd
      } else {
        rep(NA_real_, length(raw_value))
      }
      stat <- yur_cox_one(scaled_value, time, event, covariates)
      result[[j]] <- yur_cox_result_row(
        feature_id = names(proteins)[[j]], protein_mean = protein_mean,
        protein_sd = protein_sd, lag_variable = lag_variable,
        lag_panel_fallback = lag_variable == "protein_sampling_lag_days" && nzchar(panel_key),
        stat = stat
      )
      if (j %% 250L == 0L) yur_log(cfg, "Cox progress ", scope, "/", outcome_id, " ", j, "/", ncol(proteins))
    }
    result <- rbindlist(result)
    result <- merge(result, mapping, by = "feature_id", all.x = TRUE)
    result[, `:=`(outcome_id = outcome_id, outcome_label = outcomes$outcome_label[[i]], scope = scope)]
    result[, bonferroni_threshold := 0.05 / ncol(proteins)]
    result[, bonferroni_significant := is.finite(p) & p < bonferroni_threshold]
    setcolorder(result, c("scope", "outcome_id", "outcome_label", "feature_id", "protein", "panel", "olink_id",
                          "exposure_scale", "protein_mean", "protein_sd", "lag_covariate", "lag_panel_fallback",
                          "n", "events", "beta", "se", "z", "p", "hr", "ci_low", "ci_high",
                          "bonferroni_threshold", "bonferroni_significant", "mapping_status"))
    fwrite(result, out_file, na = "")
    expected_contract$scope <- scope
    expected_contract$rows <- nrow(result)
    expected_contract$generated <- yur_now()
    yur_write_json(expected_contract, contract_file)
    rows[[i]] <- result
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

yur_load_cox_cohorts <- function(cfg) {
  cohort_file <- file.path(cfg$paths$cohort, "derivation_cohort.csv.gz")
  test_file <- file.path(cfg$paths$cohort, "test_cohort.csv.gz")
  if (!file.exists(cohort_file) || !file.exists(test_file)) stop("Run cohort first.")
  derivation <- fread(cohort_file)
  test <- fread(test_file)
  derivation <- yur_normalize_eid_column(derivation)
  test <- yur_normalize_eid_column(test)
  list(
    derivation = derivation,
    test = test,
    all_meta = rbindlist(list(derivation, test), use.names = TRUE, fill = TRUE)
  )
}

yur_prepare_cox_panel <- function(cfg) {
  cohorts <- yur_load_cox_cohorts(cfg)
  all_meta <- cohorts$all_meta
  header <- names(fread(cfg$raw_protein_file, nrows = 0, showProgress = FALSE))
  eid_col <- intersect(c("eid", "id", "f.eid", "participant_id"), header)[[1]]
  features <- setdiff(header, eid_col)
  raw <- fread(cfg$raw_protein_file, select = c(eid_col, features), showProgress = TRUE, nThread = cfg$workers)
  setnames(raw, eid_col, "eid")
  raw[, eid := yur_norm_eid(eid)]
  setkey(raw, eid)
  missing_eids <- setdiff(all_meta$eid, raw$eid)
  if (length(missing_eids)) {
    stop("Raw protein input is missing ", length(missing_eids), " cohort EIDs; first=",
         paste(head(missing_eids, 10L), collapse = ","))
  }
  raw_all <- raw[J(all_meta$eid)]
  if (!identical(raw_all$eid, all_meta$eid)) {
    stop("Raw protein rows are not aligned to the requested cohort EID order.")
  }
  missing_rate <- vapply(raw_all[, ..features], function(x) mean(is.na(x)), numeric(1))
  panel_qc <- data.table(feature_id = features, missing_rate = missing_rate,
                         keep = missing_rate <= cfg$protein_missingness_max)
  mapping_all <- yur_panel_mapping(features, cfg)
  panel_qc <- merge(panel_qc, mapping_all, by = "feature_id", all.x = TRUE)
  yur_write_csv(panel_qc, file.path(cfg$paths$cox, "full_panel_missingness_qc.csv"))
  keep <- panel_qc[keep == TRUE, feature_id]
  if (length(keep) < 1000L) stop("Only ", length(keep), " proteins survived 30% missingness QC.")
  yur_write_csv(panel_qc[keep == TRUE], file.path(cfg$paths$cox, "retained_panel.csv"))
  writeLines(yur_sha_text(keep), file.path(cfg$paths$cox, "retained_panel_hash.txt"))
  yur_log(cfg, "Panel QC raw=", length(features), " retained=", length(keep))
  yur_write_json(list(
    status = "PASS", raw_protein_n = length(features), retained_protein_n = length(keep),
    panel_hash = yur_sha_text(keep), missingness_scope = "full incident cohort before derivation/hold-out modeling"
  ), file.path(cfg$paths$cox, "panel_prepare_summary.json"))
}

yur_run_cox_shard <- function(cfg) {
  retained_file <- file.path(cfg$paths$cox, "retained_panel.csv")
  if (!file.exists(retained_file)) stop("Run cox_prepare first: ", retained_file)
  retained <- fread(retained_file)
  keep <- retained$feature_id
  expected_hash <- readLines(file.path(cfg$paths$cox, "retained_panel_hash.txt"), warn = FALSE)[[1]]
  if (!identical(yur_sha_text(keep), expected_hash)) stop("Retained panel hash mismatch.")

  cohorts <- yur_load_cox_cohorts(cfg)
  derivation <- cohorts$derivation
  all_meta <- cohorts$all_meta
  outcomes <- yur_resolved_outcomes_for_cox(cfg)
  header <- names(fread(cfg$raw_protein_file, nrows = 0, showProgress = FALSE))
  eid_col <- intersect(c("eid", "id", "f.eid", "participant_id"), header)[[1]]
  missing_features <- setdiff(keep, header)
  if (length(missing_features)) stop("Raw table no longer matches retained panel; missing ", length(missing_features), " features.")
  raw <- fread(cfg$raw_protein_file, select = c(eid_col, keep), showProgress = TRUE, nThread = cfg$workers)
  setnames(raw, eid_col, "eid")
  raw[, eid := yur_norm_eid(eid)]
  setkey(raw, eid)
  missing_eids <- setdiff(all_meta$eid, raw$eid)
  if (length(missing_eids)) stop("Raw protein input is missing ", length(missing_eids), " cohort EIDs.")
  protein_all <- raw[J(all_meta$eid), ..keep]
  protein_derivation <- protein_all[match(derivation$eid, all_meta$eid)]
  mapping <- retained[, .(feature_id, protein, panel, olink_id, mapping_status)]
  full_result <- yur_run_cox_scope(cfg, all_meta, protein_all, mapping, outcomes, "full_incident")
  derivation_result <- yur_run_cox_scope(cfg, derivation, protein_derivation, mapping, outcomes, "derivation")
  yur_log(cfg, "Cox shard complete endpoints=", paste(outcomes$outcome_id, collapse = ","),
          " retained_panel_hash=", expected_hash)
  invisible(list(full = full_result, derivation = derivation_result))
}

yur_merge_cox <- function(cfg) {
  outcomes <- yur_resolved_outcomes_for_cox(cfg)
  retained_file <- file.path(cfg$paths$cox, "retained_panel.csv")
  if (!file.exists(retained_file)) stop("Run cox_prepare first.")
  retained <- fread(retained_file)
  read_scope <- function(scope) {
    rbindlist(lapply(outcomes$outcome_id, function(outcome_id) {
      id_value <- outcome_id
      path <- file.path(cfg$paths$cox, paste0(scope, "_", id_value, "_cox.csv.gz"))
      contract_path <- file.path(cfg$paths$cox, paste0(scope, "_", id_value, "_cox.contract.json"))
      if (!file.exists(path)) stop("Missing Cox shard: ", path)
      outcome_row <- outcomes[outcome_id == id_value]
      expected_contract <- yur_cox_contract_values(cfg, outcome_row)
      if (!yur_cox_contract_matches(contract_path, expected_contract)) {
        stop("Cox shard contract mismatch for ", scope, "/", id_value,
             "; stale results cannot be merged. Recompute this shard.")
      }
      x <- fread(path)
      if (nrow(x) != nrow(retained)) {
        stop("Cox shard row count mismatch for ", scope, "/", id_value,
             ": expected ", nrow(retained), ", found ", nrow(x))
      }
      if (!identical(sort(x$feature_id), sort(retained$feature_id))) {
        stop("Cox shard feature set mismatch for ", scope, "/", id_value)
      }
      x
    }), use.names = TRUE, fill = TRUE)
  }
  full_result <- read_scope("full_incident")
  derivation_result <- read_scope("derivation")
  fwrite(full_result, file.path(cfg$paths$cox, "table_s2_incident_associations.csv.gz"), na = "")
  fwrite(derivation_result, file.path(cfg$paths$cox, "derivation_associations.csv.gz"), na = "")

  summary <- rbindlist(lapply(c("full_incident", "derivation"), function(scope) {
    z <- if (scope == "full_incident") full_result else derivation_result
    z[, .(tested = .N, significant = sum(bonferroni_significant, na.rm = TRUE),
          unique_significant_proteins = uniqueN(feature_id[bonferroni_significant])),
      by = .(scope, outcome_id, outcome_label)]
  }))
  yur_write_csv(summary, file.path(cfg$paths$cox, "association_summary.csv"))
  candidates <- unique(derivation_result[bonferroni_significant == TRUE, feature_id])
  yur_write_csv(data.table(feature_id = candidates), file.path(cfg$paths$selection, "derivation_bonferroni_candidate_union.csv"))
  yur_write_json(list(
    status = "PASS", raw_protein_n = nrow(fread(file.path(cfg$paths$cox, "full_panel_missingness_qc.csv"))),
    retained_protein_n = nrow(retained),
    derivation_candidate_union_n = length(candidates), selected_outcomes = outcomes$outcome_id,
    panel_hash = readLines(file.path(cfg$paths$cox, "retained_panel_hash.txt")),
    published_anchors = list(retained_protein_n = 2920, derivation_candidate_union_n = 671)
  ), file.path(cfg$paths$cox, "cox_summary.json"))
}

yur_run_cox <- function(cfg) {
  yur_prepare_cox_panel(cfg)
  yur_run_cox_shard(cfg)
  yur_merge_cox(cfg)
}
