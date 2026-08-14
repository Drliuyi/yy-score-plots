# Shared helpers for the local UKB-PPP cardiac disease reproduction workflow.

repro_file_status <- function(path, role = NA_character_) {
  info <- file.info(path)
  data.frame(
    role = role,
    file = path,
    exists = file.exists(path),
    size_mb = ifelse(file.exists(path), round(info$size / 1024^2, 3), NA_real_),
    modified = ifelse(file.exists(path), as.character(info$mtime), NA_character_),
    stringsAsFactors = FALSE
  )
}

repro_first_existing <- function(dt, candidates, required = TRUE) {
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) > 0) return(hit[1])
  if (required) stop("Missing required columns: ", paste(candidates, collapse = ", "), call. = FALSE)
  NA_character_
}

repro_rank_norm <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ok <- !is.na(x)
  out <- rep(NA_real_, length(x))
  if (sum(ok) > 1) out[ok] <- qnorm((rank(x[ok], ties.method = "average") - 0.5) / sum(ok))
  out
}

repro_mode_value <- function(x) {
  tab <- sort(table(x, useNA = "no"), decreasing = TRUE)
  if (length(tab) == 0) return(NA)
  names(tab)[1]
}

repro_has_variation <- function(x) {
  if (is.factor(x) || is.character(x)) return(length(unique(na.omit(as.character(x)))) > 1)
  x <- suppressWarnings(as.numeric(x))
  stats::sd(x, na.rm = TRUE) > 0
}

repro_usable_covariates <- function(dt, covariates) {
  covariates <- covariates[covariates %in% names(dt)]
  covariates[vapply(covariates, function(v) repro_has_variation(dt[[v]]), logical(1))]
}

repro_clean_model_data <- function(dt, cols) {
  d <- as.data.frame(dt[, cols, drop = FALSE])
  d <- d[stats::complete.cases(d), , drop = FALSE]
  d
}

repro_resolve_outcome_date_col <- function(dt, outcome_key, date_source, outcome_map) {
  suffix <- outcome_map[["local_suffix"]][match(outcome_key, outcome_map[["outcome_key"]])]
  if (is.na(suffix)) stop("Unknown outcome key: ", outcome_key, call. = FALSE)
  candidates <- unique(c(
    paste0(date_source, "_", suffix),
    paste0("fod_ref_", suffix),
    paste0("fod_icd10_", suffix),
    paste0("fod_gp_", suffix),
    paste0("fod_srd_", suffix)
  ))
  repro_first_existing(dt, candidates, required = TRUE)
}

repro_make_birth_date <- function(dt) {
  if ("birth_date" %in% names(dt)) {
    dt$birth_date <- as.Date(dt$birth_date)
    return(dt)
  }
  if (!all(c("birth_year", "birth_month") %in% names(dt))) {
    stop("Need birth_date or birth_year + birth_month in phenotype data.", call. = FALSE)
  }
  dt$birth_date <- as.Date(sprintf("%04d-%02d-15", as.integer(dt$birth_year), as.integer(dt$birth_month)))
  dt
}

# Build one outcome axis with explicit administrative/loss/death censoring.
# This is intentionally local to the reproduction project and does not rely on
# the legacy phe.f.R t2e() helper, whose event flag ignores administrative end.
repro_construct_outcome_axis <- function(dt, outcome_date_col, baseline_col,
                                         birth_date_col, lost_date_col,
                                         death_date_col, admin_end,
                                         domain = NA_character_,
                                         same_day_policy = "exclude") {
  stopifnot(outcome_date_col %in% names(dt), baseline_col %in% names(dt))
  if (!same_day_policy %in% c("exclude", "prevalent")) {
    stop("same_day_policy must be 'exclude' or 'prevalent'.", call. = FALSE)
  }

  n <- nrow(dt)
  baseline <- as.Date(dt[[baseline_col]])
  event_date <- as.Date(dt[[outcome_date_col]])
  birth_date <- if (birth_date_col %in% names(dt)) as.Date(dt[[birth_date_col]]) else rep(as.Date(NA), n)
  lost_date <- if (lost_date_col %in% names(dt)) as.Date(dt[[lost_date_col]]) else rep(as.Date(NA), n)
  death_date <- if (death_date_col %in% names(dt)) as.Date(dt[[death_date_col]]) else rep(as.Date(NA), n)
  admin_end <- as.Date(admin_end)
  if (length(admin_end) != 1L || is.na(admin_end)) stop("Invalid administrative end date.", call. = FALSE)

  valid_lost <- !is.na(lost_date) & !is.na(baseline) & lost_date > baseline & lost_date <= admin_end
  valid_death <- !is.na(death_date) & !is.na(baseline) & death_date > baseline & death_date <= admin_end
  censor_num <- rep(as.numeric(admin_end), n)
  censor_num[valid_lost] <- pmin(censor_num[valid_lost], as.numeric(lost_date[valid_lost]))
  censor_num[valid_death] <- pmin(censor_num[valid_death], as.numeric(death_date[valid_death]))
  censor_date <- as.Date(censor_num, origin = "1970-01-01")

  same_day <- !is.na(event_date) & !is.na(baseline) & event_date == baseline
  prevalent <- !is.na(event_date) & !is.na(baseline) & event_date < baseline
  if (same_day_policy == "prevalent") prevalent <- prevalent | same_day
  incident <- !is.na(event_date) & !is.na(baseline) & event_date > baseline & event_date <= censor_date
  post_censor_event <- !is.na(event_date) & !is.na(baseline) & event_date > censor_date

  domain_happened <- rep(FALSE, n)
  if (!is.na(domain) && nzchar(domain)) {
    domain_col <- if (domain %in% names(dt)) domain else {
      count_col <- paste0("icd10Ct_", domain)
      if (count_col %in% names(dt)) count_col else NA_character_
    }
    if (!is.na(domain_col)) {
      if (inherits(dt[[domain_col]], "Date")) {
        domain_date <- as.Date(dt[[domain_col]])
        domain_happened <- !is.na(domain_date) & !is.na(baseline) & domain_date > baseline
      } else {
        domain_happened <- suppressWarnings(as.numeric(dt[[domain_col]])) > 0
        domain_happened[is.na(domain_happened)] <- FALSE
      }
    }
  }
  ambiguous_missing_date <- is.na(event_date) & domain_happened
  ambiguous_same_day <- same_day & same_day_policy == "exclude"
  ambiguous <- ambiguous_missing_date | ambiguous_same_day | is.na(baseline)

  incident_flag <- ifelse(ambiguous, NA_integer_, ifelse(incident, 1L, 0L))
  prevalent_flag <- ifelse(ambiguous, NA_integer_, ifelse(prevalent, 1L, 0L))
  final_date <- ifelse(incident, as.numeric(event_date), as.numeric(censor_date))
  final_date <- as.Date(final_date, origin = "1970-01-01")
  followup_years <- ifelse(
    !ambiguous & !prevalent & !is.na(final_date) & !is.na(baseline),
    as.numeric(final_date - baseline) / 365.25,
    NA_real_
  )
  r2e_years <- ifelse(
    prevalent,
    as.numeric(baseline - event_date) / 365.25,
    ifelse(!is.na(birth_date) & !is.na(baseline), as.numeric(baseline - birth_date) / 365.25, NA_real_)
  )
  b2e_years <- ifelse(prevalent, -r2e_years, ifelse(incident, followup_years, NA_real_))

  invalid_incident <- !is.na(incident_flag) & incident_flag == 1L &
    (is.na(event_date) | event_date > admin_end | event_date > censor_date)
  invalid_followup <- !is.na(followup_years) & !is.na(baseline) &
    followup_years > (as.numeric(admin_end - baseline) / 365.25 + 1e-8)
  if (any(invalid_incident, na.rm = TRUE) || any(invalid_followup, na.rm = TRUE)) {
    stop("OUTCOME_CENSORING_GATE_FAILED inside outcome construction.", call. = FALSE)
  }

  audit <- data.table::data.table(
    outcome_date_column = outcome_date_col,
    same_day_policy = same_day_policy,
    administrative_end = format(admin_end),
    n_nonmissing_event_date = sum(!is.na(event_date)),
    n_event_date_after_admin_end = sum(!is.na(event_date) & event_date > admin_end),
    n_event_after_person_censor = sum(post_censor_event),
    n_event_after_loss_or_death_before_admin = sum(
      post_censor_event & censor_date < admin_end,
      na.rm = TRUE
    ),
    n_same_day_excluded = sum(ambiguous_same_day),
    n_ambiguous_missing_specific_date = sum(ambiguous_missing_date),
    n_prevalent = sum(prevalent_flag == 1L, na.rm = TRUE),
    n_incident = sum(incident_flag == 1L, na.rm = TRUE),
    n_censored = sum(incident_flag == 0L & prevalent_flag == 0L, na.rm = TRUE),
    n_incident_after_admin_end = sum(incident_flag == 1L & event_date > admin_end, na.rm = TRUE),
    n_incident_after_person_censor = sum(incident_flag == 1L & event_date > censor_date, na.rm = TRUE),
    n_followup_after_admin_end = sum(invalid_followup, na.rm = TRUE),
    max_followup_years = suppressWarnings(max(followup_years, na.rm = TRUE))
  )
  if (!is.finite(audit$max_followup_years)) audit$max_followup_years <- NA_real_

  list(
    incident = incident_flag,
    prevalent = prevalent_flag,
    followup = followup_years,
    r2e = r2e_years,
    b2e = b2e_years,
    event_date = event_date,
    censor_date = censor_date,
    audit = audit
  )
}

repro_read_kinship_pairs_file <- function(path, cutoff = 0.0884, keep_ids = NULL) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(data.table::data.table())
  rel <- data.table::fread(path)
  if (ncol(rel) < 2) return(data.table::data.table())
  nms <- names(rel)
  id1_col <- grep("^(ID1|id1|eid1|f\\.eid1|IID1)$", nms, ignore.case = TRUE, value = TRUE)[1]
  id2_col <- grep("^(ID2|id2|eid2|f\\.eid2|IID2)$", nms, ignore.case = TRUE, value = TRUE)[1]
  if (is.na(id1_col)) id1_col <- nms[1]
  if (is.na(id2_col)) id2_col <- nms[2]
  kin_col <- grep("kinship|kin", nms, ignore.case = TRUE, value = TRUE)[1]
  if (is.na(kin_col) && ncol(rel) >= 5) kin_col <- nms[5]
  out <- data.table::data.table(
    id1 = as.character(rel[[id1_col]]),
    id2 = as.character(rel[[id2_col]])
  )
  if (!is.na(kin_col)) {
    out[, kinship := suppressWarnings(as.numeric(rel[[kin_col]]))]
    out <- out[is.na(kinship) | kinship > cutoff]
  }
  out <- out[!is.na(id1) & !is.na(id2) & id1 != id2]
  if (!is.null(keep_ids)) {
    keep_ids <- as.character(keep_ids)
    out <- out[id1 %in% keep_ids & id2 %in% keep_ids]
  }
  unique(out[, .(id1, id2)])
}

repro_pairs_from_related_pair_fields <- function(dt, keep_ids = NULL) {
  dt <- data.table::as.data.table(dt)
  id_col <- repro_first_existing(dt, c("id", "eid"), required = FALSE)
  if (is.na(id_col)) return(data.table::data.table())
  pair_cols <- grep("^related_pair_a[0-9]+$", names(dt), value = TRUE)
  if (length(pair_cols) == 0) pair_cols <- grep("^p22011(_a[0-9]+)?$", names(dt), value = TRUE)
  if (length(pair_cols) == 0) return(data.table::data.table())
  keep_ids <- if (is.null(keep_ids)) as.character(dt[[id_col]]) else as.character(keep_ids)
  long <- data.table::rbindlist(lapply(pair_cols, function(col) {
    x <- data.table::data.table(
      id = as.character(dt[[id_col]]),
      pair_id = as.character(dt[[col]])
    )
    x <- x[!is.na(pair_id) & pair_id != "" & pair_id != "." & id %in% keep_ids]
    x[, pair_id := paste(col, pair_id, sep = ":")]
    x
  }), use.names = TRUE, fill = TRUE)
  if (nrow(long) == 0) return(data.table::data.table())
  pairs <- long[, {
    ids <- unique(id)
    if (length(ids) < 2L) {
      NULL
    } else {
      cmb <- utils::combn(sort(ids), 2L)
      data.table::data.table(id1 = cmb[1, ], id2 = cmb[2, ])
    }
  }, by = pair_id]
  unique(pairs[, .(id1, id2)])
}

repro_samples_to_remove_from_pairs <- function(pairs, seed = 1L) {
  pairs <- data.table::as.data.table(pairs)
  if (nrow(pairs) == 0) return(character())
  pairs <- unique(pairs[!is.na(id1) & !is.na(id2) & id1 != id2, .(
    id1 = pmin(as.character(id1), as.character(id2)),
    id2 = pmax(as.character(id1), as.character(id2))
  )])
  ids <- sort(unique(c(pairs$id1, pairs$id2)))
  edge_i <- match(pairs$id1, ids)
  edge_j <- match(pairs$id2, ids)
  active <- rep(TRUE, length(edge_i))
  removed <- rep(FALSE, length(ids))
  set.seed(seed)
  while (any(active)) {
    deg <- tabulate(c(edge_i[active], edge_j[active]), nbins = length(ids))
    candidates <- which(deg == max(deg))
    pick <- if (length(candidates) == 1L) candidates else sample(candidates[order(ids[candidates])], 1L)
    removed[pick] <- TRUE
    active <- active & edge_i != pick & edge_j != pick
  }
  ids[removed]
}

repro_apply_relatedness_exclusion <- function(dt, kinship_file = NA_character_,
                                             cutoff = 0.0884, seed = 1L,
                                             skip = FALSE) {
  dt <- data.table::as.data.table(dt)
  id_col <- repro_first_existing(dt, c("id", "eid"), required = FALSE)
  if (skip || is.na(id_col)) {
    return(list(
      dt = dt,
      removed_ids = character(),
      pairs = data.table::data.table(),
      audit = data.table::data.table(
        method = ifelse(skip, "skipped_by_env", "missing_id_column"),
        source = NA_character_,
        cutoff = cutoff,
        n_pairs = 0L,
        n_related_participants = 0L,
        n_removed = 0L
      )
    ))
  }
  keep_ids <- as.character(dt[[id_col]])
  pairs <- data.table::data.table()
  method <- "none"
  source <- NA_character_
  if (!is.na(kinship_file) && nzchar(kinship_file) && file.exists(kinship_file)) {
    pairs <- repro_read_kinship_pairs_file(kinship_file, cutoff = cutoff, keep_ids = keep_ids)
    method <- "kinship_file_greedy_vertex_cover"
    source <- kinship_file
  }
  if (nrow(pairs) == 0) {
    pairs <- repro_pairs_from_related_pair_fields(dt, keep_ids = keep_ids)
    method <- if (nrow(pairs) > 0) "ukb_22011_related_pair_fields_greedy_vertex_cover" else "no_relatedness_source_found"
    source <- if (nrow(pairs) > 0) paste(grep("^related_pair_a[0-9]+$", names(dt), value = TRUE), collapse = ";") else NA_character_
  }
  removed <- repro_samples_to_remove_from_pairs(pairs, seed = seed)
  out <- dt[!as.character(get(id_col)) %in% removed]
  audit <- data.table::data.table(
    method = method,
    source = source,
    cutoff = cutoff,
    n_pairs = nrow(pairs),
    n_related_participants = length(unique(c(pairs$id1, pairs$id2))),
    n_removed = length(removed)
  )
  list(dt = out, removed_ids = removed, pairs = pairs, audit = audit)
}

repro_add_plate_well <- function(dt, plate_well_file = NA_character_) {
  dt <- data.table::as.data.table(dt)
  audit <- data.table::data.table(
    method = "not_found",
    source = NA_character_,
    plate_col = NA_character_,
    well_col = NA_character_,
    n_plate_nonmissing = 0L,
    n_well_nonmissing = 0L
  )
  plate_candidates <- c("plate", "prot.plate", "met.plate", "gen_plate", "f.30901.0.0", "p30901_i0", "p30901")
  well_candidates <- c("well", "prot.well", "met.well", "gen_well", "f.30902.0.0", "p30902_i0", "p30902")

  if (!is.na(plate_well_file) && nzchar(plate_well_file) && file.exists(plate_well_file)) {
    hdr <- names(data.table::fread(plate_well_file, nrows = 0))
    id_col <- grep("^(f\\.eid|eid|id)$", hdr, ignore.case = TRUE, value = TRUE)[1]
    plate_col <- plate_candidates[plate_candidates %in% hdr][1]
    well_col <- well_candidates[well_candidates %in% hdr][1]
    if (!is.na(id_col) && !is.na(plate_col) && !is.na(well_col)) {
      plate_source_col <- plate_col
      well_source_col <- well_col
      pw <- data.table::fread(plate_well_file, select = c(id_col, plate_col, well_col))
      data.table::setnames(pw, c(id_col, plate_col, well_col), c("eid", "plate", "well"))
      pw[, eid := as.character(eid)]
      pw[, plate := repro_as_character_preserve_integer64(plate)]
      pw[, well := repro_as_character_preserve_integer64(well)]
      pw[plate %in% c("", "0", "NA"), plate := NA_character_]
      pw[well %in% c("", "0", "NA"), well := NA_character_]
      if ("eid" %in% names(dt)) {
        dt <- merge(dt, pw, by = "eid", all.x = TRUE)
      } else if ("id" %in% names(dt)) {
        data.table::setnames(pw, "eid", "id")
        dt <- merge(dt, pw, by = "id", all.x = TRUE)
      }
      audit[, `:=`(
        method = "external_file",
        source = plate_well_file,
        plate_col = plate_source_col,
        well_col = well_source_col,
        n_plate_nonmissing = sum(!is.na(dt$plate)),
        n_well_nonmissing = sum(!is.na(dt$well))
      )]
      return(list(dt = dt, audit = audit))
    }
  }

  plate_col <- plate_candidates[plate_candidates %in% names(dt)][1]
  well_col <- well_candidates[well_candidates %in% names(dt)][1]
  if (!is.na(plate_col) && !is.na(well_col)) {
    plate_source_col <- plate_col
    well_source_col <- well_col
    dt[, plate := repro_as_character_preserve_integer64(get(plate_col))]
    dt[, well := repro_as_character_preserve_integer64(get(well_col))]
    dt[plate %in% c("", "0", "NA"), plate := NA_character_]
    dt[well %in% c("", "0", "NA"), well := NA_character_]
    audit[, `:=`(
      method = "existing_columns",
      source = "analysis_base_input",
      plate_col = plate_source_col,
      well_col = well_source_col,
      n_plate_nonmissing = sum(!is.na(dt$plate)),
      n_well_nonmissing = sum(!is.na(dt$well))
    )]
  }
  list(dt = dt, audit = audit)
}

repro_safe_factor <- function(x) {
  factor(ifelse(is.na(x), NA, as.character(x)))
}

repro_as_character_preserve_integer64 <- function(x) {
  if (inherits(x, "integer64") && requireNamespace("bit64", quietly = TRUE)) {
    return(bit64::as.character.integer64(x))
  }
  as.character(x)
}

repro_impute_continuous <- function(dt, vars, predictors) {
  predictors <- repro_usable_covariates(dt, predictors)
  for (v in vars) {
    if (!v %in% names(dt)) next
    x <- suppressWarnings(as.numeric(dt[[v]]))
    miss_name <- paste0(v, "_missing")
    final_name <- paste0(v, "_final")
    dt[[miss_name]] <- as.integer(is.na(x))
    dt[[final_name]] <- x
    miss_idx <- which(is.na(x))
    if (length(miss_idx) == 0) next
    complete_idx <- which(!is.na(x) & stats::complete.cases(dt[, predictors, drop = FALSE]))
    if (length(complete_idx) < 50 || length(predictors) == 0) {
      med <- stats::median(x, na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      dt[[final_name]][miss_idx] <- med
      next
    }
    form <- stats::reformulate(predictors, response = v)
    fit <- tryCatch(stats::lm(form, data = dt[complete_idx, , drop = FALSE]), error = function(e) e)
    if (inherits(fit, "error")) {
      med <- stats::median(x, na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      dt[[final_name]][miss_idx] <- med
    } else {
      pred <- suppressWarnings(stats::predict(fit, newdata = dt[miss_idx, predictors, drop = FALSE]))
      med <- stats::median(x, na.rm = TRUE)
      pred[!is.finite(pred)] <- med
      dt[[final_name]][miss_idx] <- pred
    }
  }
  dt
}

repro_median_impute_matrix <- function(mat) {
  for (j in seq_len(ncol(mat))) {
    miss <- is.na(mat[, j])
    if (any(miss)) {
      med <- stats::median(mat[, j], na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      mat[miss, j] <- med
    }
  }
  mat
}

repro_scale_matrix <- function(mat) {
  mat <- scale(mat)
  mat[, colSums(is.na(mat)) < nrow(mat), drop = FALSE]
}

repro_make_timevarying_data <- function(dt, outcome_key, outcome_map) {
  others <- setdiff(outcome_map$outcome_key, outcome_key)
  fu_col <- paste0(outcome_key, "_fu")
  inc_col <- paste0(outcome_key, "_inc")
  need <- unique(c("id", fu_col, inc_col, paste0(others, "_fu"), paste0(others, "_inc"), paste0(others, "_prev")))
  missing <- setdiff(need, names(dt))
  if (length(missing) > 0) stop("Missing time-to-event columns: ", paste(missing, collapse = ", "), call. = FALSE)

  tmp <- as.data.frame(dt)
  tmp$tstop_main <- suppressWarnings(as.numeric(tmp[[fu_col]]))
  tmp$event_main <- as.integer(tmp[[inc_col]] == 1)
  tmp <- tmp[!is.na(tmp$tstop_main) & is.finite(tmp$tstop_main) & tmp$tstop_main > 0 & !is.na(tmp$event_main), , drop = FALSE]
  tv <- do.call(
    survival::tmerge,
    list(data1 = tmp, data2 = tmp, id = quote(id), tstop = quote(tstop_main), event = quote(event(tstop_main, event_main)))
  )
  for (other in others) {
    time_var <- paste0(other, "_tdc_time")
    tmp[[time_var]] <- ifelse(
      tmp[[paste0(other, "_prev")]] == 1,
      0,
      ifelse(tmp[[paste0(other, "_inc")]] == 1, suppressWarnings(as.numeric(tmp[[paste0(other, "_fu")]])), NA_real_)
    )
    if (any(is.finite(tmp[[time_var]]), na.rm = TRUE)) {
      args <- list(data1 = tv, data2 = tmp, id = quote(id))
      args[[other]] <- as.call(list(as.name("tdc"), as.name(time_var)))
      tv <- do.call(survival::tmerge, args)
    } else {
      tv[[other]] <- 0
    }
  }
  tv
}

repro_empty_cox_result <- function(protein, outcome, error, n = NA_integer_, events = NA_integer_) {
  data.frame(
    Protein = protein, Outcome = outcome, beta = NA_real_, HR = NA_real_,
    CI_Lower = NA_real_, CI_Upper = NA_real_, SE = NA_real_, z = NA_real_,
    P_Value = NA_real_, n = n, events = events, error = error,
    stringsAsFactors = FALSE
  )
}

repro_fit_one_cox <- function(dt, protein, outcome, covariates, min_events = 10L) {
  cols <- unique(c("tstart", "tstop", "event", protein, covariates))
  d <- repro_clean_model_data(dt, cols)
  n_events <- sum(d$event == 1, na.rm = TRUE)
  if (nrow(d) < 50 || n_events < min_events) {
    return(repro_empty_cox_result(protein, outcome, "too few rows/events", nrow(d), n_events))
  }
  p <- suppressWarnings(as.numeric(d[[protein]]))
  if (stats::sd(p, na.rm = TRUE) <= 0 || all(is.na(p))) {
    return(repro_empty_cox_result(protein, outcome, "zero protein variance", nrow(d), n_events))
  }
  d$protein_z <- as.numeric(scale(p))
  covariates <- repro_usable_covariates(d, covariates)
  form <- stats::reformulate(c("protein_z", covariates), response = "survival::Surv(tstart, tstop, event)")
  fit <- tryCatch(survival::coxph(form, data = d, ties = "efron"), error = function(e) e)
  if (inherits(fit, "error")) {
    return(repro_empty_cox_result(protein, outcome, conditionMessage(fit), nrow(d), n_events))
  }
  sm <- summary(fit)
  if (!"protein_z" %in% rownames(sm$coefficients)) {
    return(repro_empty_cox_result(protein, outcome, "protein term dropped", nrow(d), n_events))
  }
  beta <- unname(sm$coefficients["protein_z", "coef"])
  se <- unname(sm$coefficients["protein_z", "se(coef)"])
  z <- unname(sm$coefficients["protein_z", "z"])
  pval <- unname(sm$coefficients["protein_z", "Pr(>|z|)"])
  data.frame(
    Protein = protein,
    Outcome = outcome,
    beta = beta,
    HR = exp(beta),
    CI_Lower = exp(beta - 1.96 * se),
    CI_Upper = exp(beta + 1.96 * se),
    SE = se,
    z = z,
    P_Value = pval,
    n = nrow(d),
    events = n_events,
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

repro_write_csv <- function(x, path) {
  data.table::fwrite(as.data.frame(x), path)
  invisible(path)
}

repro_file_signature <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  info <- file.info(path)
  paste(normalizePath(path, winslash = "/", mustWork = FALSE), info$size, as.numeric(info$mtime), sep = "|")
}

repro_prepare_chunk_dir <- function(chunk_dir, manifest, allow_reuse = TRUE) {
  dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)
  manifest_file <- file.path(chunk_dir, "_chunk_manifest.csv")
  expected <- data.table::data.table(
    setting = names(manifest),
    value = vapply(manifest, function(x) paste(as.character(x), collapse = "|"), character(1))
  )
  reset <- TRUE
  if (allow_reuse && file.exists(manifest_file)) {
    current <- tryCatch(data.table::fread(manifest_file), error = function(e) data.table::data.table())
    reset <- !(nrow(current) == nrow(expected) &&
      identical(as.character(current$setting), as.character(expected$setting)) &&
      identical(as.character(current$value), as.character(expected$value)))
  }
  if (reset) {
    stale <- list.files(chunk_dir, pattern = "\\.csv$", full.names = TRUE)
    stale <- stale[basename(stale) != "_chunk_manifest.csv"]
    if (length(stale) > 0) unlink(stale, force = TRUE)
    data.table::fwrite(expected, manifest_file)
    return("reset")
  }
  "reuse"
}

repro_chunk_file_matches <- function(path, protein_chunk, outcome_label) {
  if (!file.exists(path) || file.info(path)$size <= 0) return(FALSE)
  hdr <- tryCatch(names(data.table::fread(path, nrows = 0)), error = function(e) character())
  if (!all(c("Protein", "Outcome") %in% hdr)) return(FALSE)
  d <- tryCatch(data.table::fread(path, select = c("Protein", "Outcome")), error = function(e) NULL)
  if (is.null(d) || nrow(d) != length(protein_chunk)) return(FALSE)
  identical(as.character(d$Protein), as.character(protein_chunk)) &&
    all(as.character(d$Outcome) == as.character(outcome_label))
}

repro_assert_primary_complete <- function(primary, audit_dir, outcome_map, context = "primary association results") {
  protein_audit <- file.path(audit_dir, "protein_columns_used.csv")
  if (!file.exists(protein_audit)) return(invisible(TRUE))
  n_proteins <- nrow(data.table::fread(protein_audit, select = "protein"))
  expected_total <- n_proteins * nrow(outcome_map)
  if (nrow(primary) != expected_total) {
    stop(
      context, " are incomplete: found ", nrow(primary), " rows, expected ",
      expected_total, " rows (", n_proteins, " proteins x ", nrow(outcome_map),
      " outcomes). Remove stale chunks and rerun 02_primary_association_cox.R.",
      call. = FALSE
    )
  }
  expected_labels <- outcome_map$label
  counts <- primary[, .N, by = Outcome]
  bad <- counts[!Outcome %in% expected_labels | N != n_proteins]
  missing <- setdiff(expected_labels, counts$Outcome)
  if (nrow(bad) > 0 || length(missing) > 0) {
    stop(context, " have incomplete per-outcome rows. Rerun 02_primary_association_cox.R.", call. = FALSE)
  }
  invisible(TRUE)
}
