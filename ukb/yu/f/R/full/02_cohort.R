yur_binary_sex <- function(x) {
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(z))
  out[z %in% c("female", "f", "woman", "0")] <- 0L
  out[z %in% c("male", "m", "man", "1")] <- 1L
  out
}

yur_binary_current <- function(x) {
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(z))
  out[grepl("current|daily|almost daily", z)] <- 1L
  out[grepl("never|previous|former", z)] <- 0L
  suppressWarnings({
    num <- as.numeric(z)
    idx <- is.na(out) & num %in% c(0, 1)
    out[idx] <- as.integer(num[idx])
  })
  out
}

yur_binary_diabetes <- function(x) {
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(z))
  out[z %in% c("yes", "true", "1", "type 2", "type2")] <- 1L
  out[z %in% c("no", "false", "0", "none")] <- 0L
  suppressWarnings({
    num <- as.numeric(z)
    idx <- is.na(out) & num %in% c(0, 1)
    out[idx] <- as.integer(num[idx])
  })
  out
}

yur_icd10_tokens <- function(x) {
  x <- toupper(gsub("\\.", "", trimws(as.character(x))))
  lapply(strsplit(x, "[^A-Z0-9]+", perl = TRUE), function(z) {
    z[!is.na(z) & nzchar(z) & z != "NA"]
  })
}

yur_death_date_for_codes <- function(dt, prefixes, primary = "death_icd10",
                                     secondary = "death_icd10sec", death_date = "death_date",
                                     primary_tokens = NULL, secondary_tokens = NULL) {
  prefixes <- toupper(gsub("[^A-Z0-9]", "", trimws(as.character(prefixes))))
  prefixes <- unique(prefixes[nzchar(prefixes) & prefixes != "NA"])
  if (!length(prefixes)) return(as.Date(rep(NA_real_, nrow(dt)), origin = "1970-01-01"))
  required <- c(primary, secondary, death_date)
  if (any(!required %in% names(dt))) {
    stop("Death-registry endpoint derivation is missing fields: ",
         paste(setdiff(required, names(dt)), collapse = ", "))
  }
  if (is.null(primary_tokens)) primary_tokens <- yur_icd10_tokens(dt[[primary]])
  if (is.null(secondary_tokens)) secondary_tokens <- yur_icd10_tokens(dt[[secondary]])
  matched <- vapply(seq_len(nrow(dt)), function(i) {
    tokens <- c(primary_tokens[[i]], secondary_tokens[[i]])
    any(vapply(prefixes, function(code) any(startsWith(tokens, code)), logical(1)))
  }, logical(1))
  out <- yur_as_date(dt[[death_date]])
  out[!matched] <- as.Date(NA)
  out
}

yur_outcome_pairwise_qc <- function(dt, outcome_ids, prefix = "date_", as_date = TRUE) {
  pairs <- utils::combn(outcome_ids, 2L, simplify = FALSE)
  rbindlist(lapply(pairs, function(pair) {
    a <- dt[[paste0(prefix, pair[[1]])]]
    b <- dt[[paste0(prefix, pair[[2]])]]
    if (as_date) {
      a <- yur_as_date(a)
      b <- yur_as_date(b)
    }
    same_missingness <- identical(is.na(a), is.na(b))
    same_values <- same_missingness && identical(as.numeric(a), as.numeric(b))
    data.table(
      outcome_a = pair[[1]], outcome_b = pair[[2]],
      nonmissing_a = sum(!is.na(a)), nonmissing_b = sum(!is.na(b)),
      exact_duplicate = same_values
    )
  }))
}

yur_validate_endpoint_events <- function(event_summary, outcomes) {
  qc <- merge(
    event_summary,
    outcomes[, .(outcome_id, min_derivation_events, min_test_events, min_total_events)],
    by = "outcome_id", all.x = TRUE, sort = FALSE
  )
  qc[, total_events := derivation_events + test_events]
  qc[, status := fifelse(
    derivation_events >= min_derivation_events & test_events >= min_test_events &
      total_events >= min_total_events,
    "PASS", "FAIL_LOW_EVENTS"
  )]
  qc
}

yur_outcome_definition_hash <- function(outcomes) {
  ordered <- copy(outcomes)[order(outcome_id)]
  yur_sha_text(apply(ordered[, .(
    outcome_id, resolved_date_fields, death_icd10_prefixes, source_policy,
    min_derivation_events, min_test_events, min_total_events, definition_hash
  )], 1L, paste, collapse = "|"))
}

yur_mode <- function(x) {
  z <- x[!is.na(x)]
  if (!length(z)) return(NA)
  names(sort(table(z), decreasing = TRUE))[1]
}

yur_score2_vector <- function(age, sex, smoker, sbp, diabetes, tc, hdl) {
  out <- rep(NA_real_, length(age))
  for (i in seq_along(out)) {
    values <- c(age[[i]], sex[[i]], smoker[[i]], sbp[[i]], diabetes[[i]], tc[[i]], hdl[[i]])
    if (any(!is.finite(values))) next
    male <- sex[[i]] == 1
    a <- age[[i]]; sm <- smoker[[i]]; s <- sbp[[i]]; dm <- diabetes[[i]]
    t <- tc[[i]]; h <- hdl[[i]]
    if (a < 70) {
      if (male) {
        lp <- .3742*(a-60)/5 + .6012*sm + .2777*(s-120)/20 + .6457*dm +
          .1458*(t-6) - .2698*(h-1.3)/.5 - .0755*(a-60)/5*sm -
          .0255*(a-60)/5*(s-120)/20 - .0281*(a-60)/5*(t-6) +
          .0426*(a-60)/5*(h-1.3)/.5 - .0983*(a-60)/5*dm
        base <- 1 - .9605^exp(lp); sc1 <- -.5699; sc2 <- .7476
      } else {
        lp <- .4648*(a-60)/5 + .7744*sm + .3131*(s-120)/20 + .8096*dm +
          .1002*(t-6) - .2606*(h-1.3)/.5 - .1088*(a-60)/5*sm -
          .0277*(a-60)/5*(s-120)/20 - .0226*(a-60)/5*(t-6) +
          .0613*(a-60)/5*(h-1.3)/.5 - .1272*(a-60)/5*dm
        base <- 1 - .9776^exp(lp); sc1 <- -.7380; sc2 <- .7019
      }
    } else {
      if (male) {
        lp <- .0634*(a-73) + .4245*dm + .3524*sm + .0094*(s-150) + .0850*(t-6) -
          .3564*(h-1.4) - .0174*(a-73)*dm - .0247*(a-73)*sm -
          .0005*(a-73)*(s-150) + .0073*(a-73)*(t-6) + .0091*(a-73)*(h-1.4)
        base <- 1 - .7576^exp(lp-.0929); sc1 <- -.34; sc2 <- 1.19
      } else {
        lp <- .0789*(a-73) + .6010*dm + .4921*sm + .0102*(s-150) + .0605*(t-6) -
          .3040*(h-1.4) - .0107*(a-73)*dm - .0255*(a-73)*sm -
          .0004*(a-73)*(s-150) + .0067*(a-73)*(t-6) + .0094*(a-73)*(h-1.4)
        base <- 1 - .8082^exp(lp-.229); sc1 <- -.52; sc2 <- 1.01
      }
    }
    out[[i]] <- 1 - exp(-exp(sc1 + sc2 * log(-log(1-base))))
  }
  pmin(1, pmax(0, out))
}

yur_impute_score2 <- function(train, test) {
  params <- list()
  for (v in c("total_cholesterol", "hdl")) {
    params[[v]] <- median(train[[v]], na.rm = TRUE)
    train[is.na(get(v)), (v) := params[[v]]]
    test[is.na(get(v)), (v) := params[[v]]]
  }
  for (v in c("smoking_current", "diabetes")) {
    params[[v]] <- as.numeric(yur_mode(train[[v]]))
    train[is.na(get(v)), (v) := params[[v]]]
    test[is.na(get(v)), (v) := params[[v]]]
  }
  sex_medians <- train[, .(sbp_median = median(sbp, na.rm = TRUE)), by = sex]
  params$sbp_by_sex <- split(sex_medians$sbp_median, sex_medians$sex)
  for (sx in sex_medians$sex) {
    value <- sex_medians[sex == sx, sbp_median]
    train[sex == sx & is.na(sbp), sbp := value]
    test[sex == sx & is.na(sbp), sbp := value]
  }
  train[, score2_raw := yur_score2_vector(age, sex, smoking_current, sbp, diabetes, total_cholesterol, hdl)]
  test[, score2_raw := yur_score2_vector(age, sex, smoking_current, sbp, diabetes, total_cholesterol, hdl)]
  list(train = train, test = test, parameters = params)
}

yur_score2_distribution_qc <- function(train, test, cfg) {
  qc <- cfg$score2_qc
  rbindlist(lapply(names(list(derivation = train, test = test)), function(scope) {
    x <- list(derivation = train, test = test)[[scope]]
    input_rows <- rbindlist(lapply(c("total_cholesterol", "hdl", "sbp"), function(variable) {
      values <- suppressWarnings(as.numeric(x[[variable]]))
      finite <- values[is.finite(values)]
      med <- median(finite, na.rm = TRUE)
      pass <- switch(
        variable,
        total_cholesterol = med >= qc$total_cholesterol_median_min &&
          med <= qc$total_cholesterol_median_max &&
          mean(finite > 20) <= qc$total_cholesterol_gt20_fraction_max,
        hdl = med >= qc$hdl_median_min && med <= qc$hdl_median_max,
        sbp = med >= qc$sbp_median_min && med <= qc$sbp_median_max
      )
      data.table(
        scope = scope, variable = variable, nonmissing_n = length(finite),
        q01 = unname(quantile(finite, .01)), median = med, q99 = unname(quantile(finite, .99)),
        boundary_fraction = if (variable == "total_cholesterol") mean(finite > 20) else NA_real_,
        distinct_n = uniqueN(finite), status = if (pass) "PASS" else "FAIL_IMPLAUSIBLE_RANGE"
      )
    }))
    score <- suppressWarnings(as.numeric(x$score2_raw))
    finite_score <- score[is.finite(score)]
    boundary_fraction <- mean(finite_score <= 1e-12 | finite_score >= 1 - 1e-12)
    score_pass <- length(finite_score) == nrow(x) &&
      boundary_fraction <= qc$score_boundary_fraction_max &&
      uniqueN(round(finite_score, 12)) >= qc$score_distinct_min
    score_row <- data.table(
      scope = scope, variable = "score2_raw", nonmissing_n = length(finite_score),
      q01 = unname(quantile(finite_score, .01)), median = median(finite_score),
      q99 = unname(quantile(finite_score, .99)), boundary_fraction = boundary_fraction,
      distinct_n = uniqueN(round(finite_score, 12)),
      status = if (score_pass) "PASS" else "FAIL_SCORE2_SATURATION"
    )
    rbind(input_rows, score_row, fill = TRUE)
  }))
}

yur_cohort_complete <- function(cfg) {
  summary_file <- file.path(cfg$paths$cohort, "cohort_summary.json")
  source_file <- file.path(cfg$paths$cohort, "score2_source_manifest.csv")
  qc_file <- file.path(cfg$paths$cohort, "score2_input_output_qc.csv")
  contract_file <- file.path(cfg$paths$cohort, "cohort_contract_hash.txt")
  if (!all(file.exists(c(summary_file, source_file, qc_file, contract_file)))) return(FALSE)
  summary <- read_json(summary_file, simplifyVector = TRUE)
  source <- fread(source_file)
  qc <- fread(qc_file)
  identical(summary$score2_qc_status, "PASS") &&
    identical(source[canonical == "total_cholesterol", source_field], "bb_TC") &&
    identical(source[canonical == "hdl", source_field], "bb_HDL") &&
    all(qc$status == "PASS") && length(readLines(contract_file, warn = FALSE)) == 1L
}

yur_build_cohort <- function(cfg) {
  preflight <- file.path(cfg$paths$preflight, "preflight_summary.json")
  if (!file.exists(preflight)) stop("Run preflight first.")

  ph <- as.data.table(readRDS(cfg$phenotype_rds))
  fields <- fread(file.path(cfg$paths$preflight, "local_field_resolution.csv"))
  fields <- fields[!is.na(resolved) & source_type == "phenotype_rds"]
  rename_map <- setNames(fields$resolved, fields$canonical)
  outcomes <- fread(file.path(cfg$paths$preflight, "outcome_field_resolution.csv"))

  raw_header <- names(fread(cfg$raw_protein_file, nrows = 0, showProgress = FALSE))
  raw_eid <- intersect(c("eid", "id", "f.eid", "participant_id"), raw_header)[[1]]
  protein_ids <- fread(cfg$raw_protein_file, select = raw_eid, showProgress = FALSE)
  setnames(protein_ids, raw_eid, "eid")
  protein_ids[, eid := yur_norm_eid(eid)]

  death_source_fields <- c("death_icd10", "death_icd10sec")
  needed <- unique(c(
    unname(rename_map),
    unlist(strsplit(outcomes$resolved_date_fields, ";", fixed = TRUE)),
    death_source_fields
  ))
  needed <- needed[needed %in% names(ph)]
  x <- copy(ph[, c(unique(rename_map[["eid"]]), setdiff(needed, rename_map[["eid"]])), with = FALSE])
  for (canonical in names(rename_map)) {
    source <- rename_map[[canonical]]
    if (source %in% names(x) && canonical != source) setnames(x, source, canonical)
  }
  x[, eid := yur_norm_eid(eid)]
  x <- x[eid %in% protein_ids$eid]
  x[, baseline_date := yur_as_date(baseline_date)]
  if ("death_date" %in% names(x)) x[, death_date := yur_as_date(death_date)]

  supplemental_file <- file.path(cfg$paths$preflight, "supplemental_covariates.csv.gz")
  if (!file.exists(supplemental_file)) stop("Run preflight first; supplemental fasting cache is missing.")
  supplemental <- fread(supplemental_file, showProgress = FALSE)
  supplemental <- yur_normalize_eid_column(supplemental)
  x <- merge(x, supplemental[, .(eid, fasting_time)], by = "eid", all.x = TRUE, sort = FALSE)
  x[, blood_collection_season := yur_blood_collection_season(baseline_date)]
  technical <- yur_derive_technical_covariates(
    x[, .(eid, baseline_date, protein_plate)], cfg$olink_processing_start_date_file
  )
  technical$values <- yur_normalize_eid_column(technical$values)
  x <- merge(x, technical$values, by = "eid", all.x = TRUE, sort = FALSE)
  technical_qc <- copy(technical$coverage)
  technical_qc[, `:=`(
    processing_resource = cfg$olink_processing_start_date_file,
    processing_resource_sha256 = yur_sha_file(cfg$olink_processing_start_date_file),
    lag_definition = "panel processing start date minus baseline blood collection date"
  )]
  yur_write_csv(technical_qc, file.path(cfg$paths$cohort, "technical_covariate_qc.csv"))

  x[, sex := yur_binary_sex(sex)]
  x[, ethnicity_white := fifelse(is.na(ethnicity), NA_integer_, as.integer(tolower(as.character(ethnicity)) %in% c("white", "white european", "1")))]
  x[, smoking_current := yur_binary_current(smoking)]
  x[, alcohol_current := yur_binary_current(alcohol)]
  x[, diabetes := yur_binary_diabetes(diabetes)]

  cutoff <- as.Date(cfg$followup_cutoff)
  outcome_ids <- outcomes$outcome_id
  death_primary_tokens <- yur_icd10_tokens(x$death_icd10)
  death_secondary_tokens <- yur_icd10_tokens(x$death_icd10sec)
  source_coverage <- vector("list", nrow(outcomes))
  for (i in seq_len(nrow(outcomes))) {
    id <- outcomes$outcome_id[[i]]
    date_fields <- strsplit(outcomes$resolved_date_fields[[i]], ";", fixed = TRUE)[[1]]
    hospital_or_first_date <- yur_min_date(x, date_fields)
    death_codes <- strsplit(outcomes$death_icd10_prefixes[[i]], ";", fixed = TRUE)[[1]]
    death_event_date <- yur_death_date_for_codes(
      x, death_codes,
      primary_tokens = death_primary_tokens,
      secondary_tokens = death_secondary_tokens
    )
    event_date <- yur_min_date(
      data.table(hospital_or_first_date = hospital_or_first_date, death_event_date = death_event_date),
      c("hospital_or_first_date", "death_event_date")
    )
    x[[paste0("date_", id)]] <- event_date
    x[[paste0("prev_", id)]] <- !is.na(event_date) & event_date <= x$baseline_date
    source_coverage[[i]] <- data.table(
      outcome_id = id, outcome_label = outcomes$outcome_label[[i]],
      source_policy = outcomes$source_policy[[i]],
      resolved_date_fields = outcomes$resolved_date_fields[[i]],
      death_icd10_prefixes = outcomes$death_icd10_prefixes[[i]],
      hospital_or_first_nonmissing = sum(!is.na(hospital_or_first_date)),
      death_registry_nonmissing = sum(!is.na(death_event_date)),
      composite_nonmissing = sum(!is.na(event_date)),
      definition_hash = outcomes$definition_hash[[i]]
    )
  }
  source_coverage <- rbindlist(source_coverage)
  yur_write_csv(source_coverage, file.path(cfg$paths$cohort, "outcome_source_coverage.csv"))
  pairwise_date_qc <- yur_outcome_pairwise_qc(x, outcome_ids)
  yur_write_csv(pairwise_date_qc, file.path(cfg$paths$cohort, "outcome_pairwise_duplicate_qc.csv"))
  if (any(pairwise_date_qc$exact_duplicate)) {
    bad <- pairwise_date_qc[exact_duplicate == TRUE]
    stop(
      "Endpoint definition QC failed: exact duplicate date vectors: ",
      paste(paste(bad$outcome_a, bad$outcome_b, sep = "="), collapse = ", ")
    )
  }
  prev_cols <- paste0("prev_", outcome_ids)
  x[, any_baseline_cvd := rowSums(as.matrix(.SD), na.rm = TRUE) > 0, .SDcols = prev_cols]
  other_prev_cols <- setdiff(prev_cols, "prev_cad")
  x[, other_baseline_cvd := rowSums(as.matrix(.SD), na.rm = TRUE) > 0, .SDcols = other_prev_cols]
  x <- x[!is.na(baseline_date) & !is.na(sex)]

  censor <- rep(cutoff, nrow(x))
  if ("death_date" %in% names(x)) {
    death_num <- as.numeric(x$death_date)
    censor_num <- pmin(as.numeric(censor), death_num, na.rm = TRUE)
    censor_num[!is.finite(censor_num)] <- as.numeric(cutoff)
    censor <- as.Date(censor_num, origin = "1970-01-01")
  }
  x[, censor_date := censor]

  for (id in outcome_ids) {
    event_date <- x[[paste0("date_", id)]]
    event <- !x[[paste0("prev_", id)]] & !is.na(event_date) & event_date > x$baseline_date & event_date <= x$censor_date
    end_date <- as.Date(ifelse(event, as.numeric(event_date), as.numeric(x$censor_date)), origin = "1970-01-01")
    x[[paste0("event_", id)]] <- as.integer(event)
    x[[paste0("time_", id)]] <- pmax(1, as.numeric(end_date - x$baseline_date))
  }

  incident <- x[any_baseline_cvd == FALSE]
  yang <- x[prev_cad == TRUE & other_baseline_cvd == FALSE]
  yang[, years_since_cad := as.numeric(baseline_date - date_cad) / 365.25]
  yang <- yang[is.finite(years_since_cad) & years_since_cad >= 0]
  set.seed(cfg$split_seed)
  shuffled <- sample(incident$eid)
  n_train <- floor(length(shuffled) * cfg$train_fraction)
  train_ids <- shuffled[seq_len(n_train)]
  test_ids <- shuffled[-seq_len(n_train)]
  train <- incident[match(train_ids, eid)]
  test <- incident[match(test_ids, eid)]
  score2 <- yur_impute_score2(train, test)
  train <- score2$train; test <- score2$test
  score2_qc <- yur_score2_distribution_qc(train, test, cfg)
  yur_write_csv(score2_qc, file.path(cfg$paths$cohort, "score2_input_output_qc.csv"))
  if (any(score2_qc$status != "PASS")) {
    bad <- score2_qc[status != "PASS"]
    stop(
      "SCORE2 cohort QC failed: ",
      paste(sprintf("%s/%s median=%.5f status=%s", bad$scope, bad$variable, bad$median, bad$status), collapse = "; "),
      call. = FALSE
    )
  }
  score2_sources <- fields[canonical %in% c("age", "sex", "smoking", "sbp", "total_cholesterol", "hdl", "diabetes"),
                           .(canonical, source_field = resolved, source_type)]
  yur_write_csv(score2_sources, file.path(cfg$paths$cohort, "score2_source_manifest.csv"))
  score2$parameters$source_fields <- setNames(score2_sources$source_field, score2_sources$canonical)
  score2$parameters$source_manifest_hash <- yur_sha_text(apply(score2_sources, 1L, paste, collapse = "|"))

  set.seed(cfg$inner_fold_seed)
  foldid <- sample(rep(seq_len(cfg$inner_folds), length.out = nrow(train)))
  fold_table <- data.table(eid = train$eid, foldid = foldid)
  split_table <- rbind(
    data.table(eid = train$eid, split = "derivation"),
    data.table(eid = test$eid, split = "test")
  )

  panel_lag_columns <- grep("^protein_sampling_lag_days__", names(x), value = TRUE)
  all_columns <- unique(c(
    "eid", "baseline_date", "censor_date", "age", "sex", "ethnicity_white", "tdi", "assessment_center",
    "smoking", "smoking_current", "alcohol", "alcohol_current", "bmi", "sbp",
    "blood_collection_season", "protein_plate", "protein_sampling_lag_days", panel_lag_columns, "fasting_time",
    "total_cholesterol", "hdl", "diabetes", "score2_raw",
    paste0("event_", outcome_ids), paste0("time_", outcome_ids)
  ))
  train <- train[, intersect(all_columns, names(train)), with = FALSE]
  test <- test[, intersect(all_columns, names(test)), with = FALSE]
  fwrite(train, file.path(cfg$paths$cohort, "derivation_cohort.csv.gz"), na = "")
  fwrite(test, file.path(cfg$paths$cohort, "test_cohort.csv.gz"), na = "")
  yang_columns <- intersect(
    c("eid", "baseline_date", "date_cad", "years_since_cad", "age", "sex", "ethnicity_white",
      "smoking_current", "bmi", "sbp"),
    names(yang)
  )
  fwrite(yang[, ..yang_columns], file.path(cfg$paths$cohort, "yang_cad_auxiliary.csv.gz"), na = "")
  yur_write_csv(split_table, file.path(cfg$paths$cohort, "split_eid.csv"))
  yur_write_csv(fold_table, file.path(cfg$paths$cohort, "foldid.csv"))
  writeLines(yur_sha_text(paste(split_table$eid, split_table$split, sep = ":")), file.path(cfg$paths$cohort, "split_hash.txt"))
  writeLines(yur_sha_text(paste(fold_table$eid, fold_table$foldid, sep = ":")), file.path(cfg$paths$cohort, "foldid_hash.txt"))
  saveRDS(score2$parameters, file.path(cfg$paths$cache, "score2_imputation_parameters.rds"))

  event_summary <- rbindlist(lapply(outcome_ids, function(id) {
    data.table(
      outcome_id = id,
      derivation_n = nrow(train), derivation_events = sum(train[[paste0("event_", id)]], na.rm = TRUE),
      test_n = nrow(test), test_events = sum(test[[paste0("event_", id)]], na.rm = TRUE)
    )
  }))
  event_qc <- yur_validate_endpoint_events(event_summary, outcomes)
  yur_write_csv(event_qc, file.path(cfg$paths$cohort, "endpoint_event_summary.csv"))
  if (any(event_qc$status != "PASS")) {
    bad <- event_qc[status != "PASS"]
    stop(
      "Endpoint event-count QC failed: ",
      paste(sprintf(
        "%s derivation=%d test=%d total=%d",
        bad$outcome_id, bad$derivation_events, bad$test_events, bad$total_events
      ), collapse = "; ")
    )
  }

  event_pairwise_qc <- yur_outcome_pairwise_qc(
    rbindlist(list(train, test), use.names = TRUE, fill = TRUE), outcome_ids,
    prefix = "event_", as_date = FALSE
  )
  setnames(event_pairwise_qc, "exact_duplicate", "exact_duplicate_event_vector")
  yur_write_csv(event_pairwise_qc, file.path(cfg$paths$cohort, "outcome_event_pairwise_qc.csv"))
  if (any(event_pairwise_qc$exact_duplicate_event_vector)) {
    bad <- event_pairwise_qc[exact_duplicate_event_vector == TRUE]
    stop(
      "Endpoint event-vector QC failed: exact duplicates: ",
      paste(paste(bad$outcome_a, bad$outcome_b, sep = "="), collapse = ", ")
    )
  }

  definition_hash <- yur_outcome_definition_hash(outcomes)
  cohort_contract_hash <- yur_sha_text(c(
    definition_hash,
    score2$parameters$source_manifest_hash,
    apply(score2_qc, 1L, paste, collapse = "|"),
    paste(split_table$eid, split_table$split, sep = ":"),
    unlist(lapply(outcome_ids, function(id) {
      c(
        paste0(id, ":event:", c(train[[paste0("event_", id)]], test[[paste0("event_", id)]])),
        paste0(id, ":time:", c(train[[paste0("time_", id)]], test[[paste0("time_", id)]]))
      )
    }))
  ))
  writeLines(definition_hash, file.path(cfg$paths$cohort, "outcome_definition_hash.txt"))
  writeLines(cohort_contract_hash, file.path(cfg$paths$cohort, "cohort_contract_hash.txt"))
  flow <- data.table(
    step = c("phenotype_rows", "with_protein_and_baseline", "free_of_all_14_baseline_CVD", "derivation", "holdout", "cad_yang_auxiliary"),
    n = c(nrow(ph), nrow(x), nrow(incident), nrow(train), nrow(test), nrow(yang))
  )
  yur_write_csv(flow, file.path(cfg$paths$cohort, "cohort_flow.csv"))

  characteristics <- data.table(
    characteristic = c("Participants", "Age mean", "Women", "White", "Current smoking", "BMI mean", "SBP mean"),
    value = c(
      as.character(nrow(incident)), sprintf("%.2f", mean(incident$age, na.rm = TRUE)),
      sprintf("%d (%.1f%%)", sum(incident$sex == 0, na.rm = TRUE), 100 * mean(incident$sex == 0, na.rm = TRUE)),
      sprintf("%d (%.1f%%)", sum(incident$ethnicity_white == 1, na.rm = TRUE), 100 * mean(incident$ethnicity_white == 1, na.rm = TRUE)),
      sprintf("%d (%.1f%%)", sum(incident$smoking_current == 1, na.rm = TRUE), 100 * mean(incident$smoking_current == 1, na.rm = TRUE)),
      sprintf("%.2f", mean(incident$bmi, na.rm = TRUE)), sprintf("%.2f", mean(incident$sbp, na.rm = TRUE))
    )
  )
  yur_write_csv(characteristics, file.path(cfg$paths$cohort, "table_s1_characteristics.csv"))
  yur_write_json(list(
    status = "PASS", outcome_n = length(outcome_ids), incident_cohort_n = nrow(incident),
    derivation_n = nrow(train), test_n = nrow(test), cad_yang_auxiliary_n = nrow(yang), split_seed = cfg$split_seed,
    split_hash = readLines(file.path(cfg$paths$cohort, "split_hash.txt")),
    score2_source_manifest_hash = score2$parameters$source_manifest_hash,
    score2_qc_status = "PASS",
    technical_covariates = list(
      blood_collection_season = "derived from baseline assessment date",
      fasting_time = "UKB field 74 from raw pheno.tsv.gz",
      protein_sampling_lag_days = "panel-specific UKB Resource 1019 processing date minus baseline date",
      panel_lag_columns = panel_lag_columns
    ),
    note = "Split is unstratified, matching the reported two-thirds/one-third design; seed is a frozen local adaptation."
  ), file.path(cfg$paths$cohort, "cohort_summary.json"))
}
