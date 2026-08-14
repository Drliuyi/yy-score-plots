yur_sources <- function(cfg) {
  required <- c("readxl", "data.table", "digest", "jsonlite")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing source-audit packages: ", paste(missing, collapse = ", "))

  workbook <- cfg$supplement_workbook_file
  methods_pdf <- cfg$supplement_methods_file
  if (!file.exists(workbook) || !file.exists(methods_pdf)) {
    stop(
      "Official supplementary workbook/PDF is missing. Resolved files: ",
      workbook, "; ", methods_pdf, ". Run ./yu.sh setup or supply explicit paths."
    )
  }

  expected_manifest <- file.path(cfg$project_dir, "references", "SOURCE_MANIFEST.csv")
  if (!file.exists(expected_manifest)) stop("Reference source manifest is missing: ", expected_manifest)
  expected <- fread(expected_manifest)
  expected <- expected[source_id %in% c("supplement_tables", "supplement_methods")]
  if (nrow(expected) != 2L || any(!nzchar(expected$sha256))) {
    stop("Reference source manifest does not contain both frozen supplement hashes.")
  }
  actual_hash <- c(
    supplement_tables = yur_sha_file(workbook),
    supplement_methods = yur_sha_file(methods_pdf)
  )
  expected_hash <- setNames(tolower(expected$sha256), expected$source_id)
  mismatched <- names(actual_hash)[tolower(actual_hash) != expected_hash[names(actual_hash)]]
  if (length(mismatched)) {
    stop(
      "Official supplement SHA256 validation failed: ", paste(mismatched, collapse = ", "),
      ". Restore the tracked publication files before analysis."
    )
  }

  sheets <- readxl::excel_sheets(workbook)
  inventory <- rbindlist(lapply(sheets, function(sheet) {
    x <- readxl::read_excel(workbook, sheet = sheet, n_max = 5, col_names = FALSE)
    data.table(sheet = sheet, preview = paste(na.omit(as.character(unlist(x[1, ]))), collapse = " | "))
  }))
  yur_write_csv(inventory, file.path(cfg$paths$sources, "supplement_sheet_inventory.csv"))

  data_sheets <- grep("^S[0-9]+$", sheets, value = TRUE)
  for (sheet in data_sheets) {
    x <- as.data.table(readxl::read_excel(workbook, sheet = sheet, skip = 1))
    fwrite(x, file.path(cfg$paths$source_tables, paste0(tolower(sheet), "_official.csv.gz")), na = "")
  }

  for (sheet in c("S23", "S24", "S25", "S26")) {
    x <- as.data.table(readxl::read_excel(workbook, sheet = sheet, skip = 1))
    yur_write_csv(x, file.path(cfg$paths$source_tables, paste0(tolower(sheet), "_official.csv")))
  }

  s10 <- as.data.table(readxl::read_excel(workbook, sheet = "S10", skip = 1))
  s11 <- as.data.table(readxl::read_excel(workbook, sheet = "S11", skip = 1))
  s12 <- as.data.table(readxl::read_excel(workbook, sheet = "S12", skip = 1))
  s12[, rank_numeric__ := suppressWarnings(as.numeric(Rank))]
  s12 <- s12[
    !is.na(Proteins) & !is.na(Disease) & is.finite(rank_numeric__) &
      rank_numeric__ >= 1 & rank_numeric__ == floor(rank_numeric__)
  ]
  s12[, rank_numeric__ := NULL]
  yur_write_csv(s10, file.path(cfg$paths$source_tables, "s10_prediction_reference.csv"))
  yur_write_csv(s11, file.path(cfg$paths$source_tables, "s11_nri_idi_reference.csv"))
  yur_write_csv(s12, file.path(cfg$paths$source_tables, "s12_importance_reference.csv.gz"))

  cad <- s12[Disease == "CAD" & !is.na(Proteins)]
  published <- list(
    cohort_n = 53026L,
    incident_free_n = 46818L,
    any_incident_cvd_n = 9096L,
    raw_unique_protein_n = 2923L,
    retained_protein_n = 2920L,
    significant_association_n = 3089L,
    associated_unique_protein_n = 892L,
    derivation_candidate_union_n = 671L,
    prediction_union_n = 257L,
    cad_importance_rows = nrow(cad),
    split = "two-thirds derivation / one-third hold-out",
    bootstrap_n = 1000L
  )
  yur_write_json(published, file.path(cfg$paths$sources, "published_numeric_anchors.json"))

  manifest <- data.table(
    source_id = c("article", "supplement_workbook", "supplement_methods"),
    title = c(
      "Systematic analyses uncover plasma proteins linked to incident cardiovascular diseases",
      "Official Tables S1-S26", "Official methods and Figures S1-S5"
    ),
    doi = c("10.1093/procel/pwaf072", "10.1093/procel/pwaf072", "10.1093/procel/pwaf072"),
    path = c("PMC12987571", workbook, methods_pdf),
    sha256 = c(NA_character_, yur_sha_file(workbook), yur_sha_file(methods_pdf)),
    code_repository = c("NOT_RELEASED", "NOT_APPLICABLE", "NOT_APPLICABLE")
  )
  yur_write_csv(manifest, file.path(cfg$paths$sources, "official_source_manifest.csv"))
  yur_write_json(list(
    status = "PASS",
    author_code_repository_found = FALSE,
    evidence = "Article Code availability lists software only; exact-title, DOI, article-code and author repository searches returned no repository.",
    workbook_sheets = length(sheets),
    important_limitation = "Independent source-locked implementation; not line-for-line author-code replication."
  ), file.path(cfg$paths$sources, "code_availability_audit.json"))
}

yur_supplemental_source_identity <- function(path) {
  info <- file.info(path)
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  prefix <- readBin(con, what = "raw", n = 65536L)
  yur_sha_text(c(normalizePath(path, winslash = "/", mustWork = TRUE), info$size,
                 as.numeric(info$mtime), digest(prefix, algo = "sha256", serialize = FALSE)))
}

yur_prepare_supplemental_covariates <- function(cfg) {
  cache_file <- file.path(cfg$paths$preflight, "supplemental_covariates.csv.gz")
  manifest_file <- file.path(cfg$paths$preflight, "supplemental_covariates_manifest.json")
  source_identity <- yur_supplemental_source_identity(cfg$raw_phenotype_file)
  if (file.exists(cache_file) && file.exists(manifest_file) && !cfg$force) {
    manifest <- read_json(manifest_file, simplifyVector = TRUE)
    if (identical(as.character(manifest$source_identity), source_identity)) {
      cached <- fread(cache_file, showProgress = FALSE)
      return(yur_normalize_eid_column(cached))
    }
  }
  yur_log(cfg, "Extracting one-time fasting-time cache from ", cfg$raw_phenotype_file)
  x <- tryCatch(
    fread(
      cfg$raw_phenotype_file, select = c("eid", "p74_i0"),
      showProgress = TRUE, nThread = cfg$workers
    ),
    error = function(e) stop(
      "Cannot extract UKB baseline fasting time; raw phenotype must contain eid and p74_i0. ",
      conditionMessage(e), call. = FALSE
    )
  )
  setnames(x, c("eid", "p74_i0"), c("eid", "fasting_time"))
  x <- yur_normalize_eid_column(x)
  x[, fasting_time := suppressWarnings(as.numeric(fasting_time))]
  x[fasting_time < 0, fasting_time := NA_real_]
  x <- unique(x[!is.na(eid)], by = "eid")
  fwrite(x, cache_file, na = "")
  yur_write_json(list(
    status = "PASS", source = cfg$raw_phenotype_file, source_field = "p74_i0",
    source_identity = source_identity, identity_method = "path_size_mtime_and_prefix_sha256",
    participants = nrow(x), nonmissing = sum(!is.na(x$fasting_time)), generated = yur_now()
  ), manifest_file)
  x
}

yur_numeric_source_summary <- function(values, canonical, source, selected, cfg) {
  values <- suppressWarnings(as.numeric(values))
  finite <- values[is.finite(values)]
  quantile_or_na <- function(p) if (length(finite)) unname(quantile(finite, p, na.rm = TRUE)) else NA_real_
  median_value <- if (length(finite)) median(finite) else NA_real_
  gt20_fraction <- if (length(finite)) mean(finite > 20) else NA_real_
  qc <- cfg$score2_qc
  pass <- switch(
    canonical,
    total_cholesterol = is.finite(median_value) &&
      median_value >= qc$total_cholesterol_median_min &&
      median_value <= qc$total_cholesterol_median_max &&
      gt20_fraction <= qc$total_cholesterol_gt20_fraction_max,
    hdl = is.finite(median_value) && median_value >= qc$hdl_median_min && median_value <= qc$hdl_median_max,
    sbp = is.finite(median_value) && median_value >= qc$sbp_median_min && median_value <= qc$sbp_median_max,
    age = is.finite(median_value) && median_value >= 35 && median_value <= 85,
    FALSE
  )
  data.table(
    canonical = canonical, source = source, selected = selected,
    total_n = length(values), nonmissing_n = length(finite), missing_rate = 1 - length(finite) / length(values),
    q01 = quantile_or_na(.01), median = median_value, q99 = quantile_or_na(.99),
    greater_than_20_fraction = gt20_fraction,
    expected_unit = fifelse(canonical %in% c("total_cholesterol", "hdl"), "mmol/L",
                            fifelse(canonical == "sbp", "mmHg", "years")),
    status = if (pass) "PASS" else "FAIL_IMPLAUSIBLE_RANGE"
  )
}

yur_score2_source_qc <- function(ph, ph_index, field_resolution, cfg) {
  field_map <- fread(cfg$field_map_file)
  variables <- c("age", "sbp", "total_cholesterol", "hdl")
  rows <- rbindlist(lapply(variables, function(var_name) {
    candidates <- strsplit(field_map[canonical == var_name, candidates][[1]], ";", fixed = TRUE)[[1]]
    candidates <- candidates[candidates %in% names(ph)]
    selected_source <- field_resolution[canonical == var_name, resolved][[1]]
    rbindlist(lapply(candidates, function(source) {
      yur_numeric_source_summary(ph[[source]][ph_index], var_name, source, identical(source, selected_source), cfg)
    }))
  }))
  rows
}

yur_preflight_complete <- function(cfg) {
  summary_file <- file.path(cfg$paths$preflight, "preflight_summary.json")
  source_qc_file <- file.path(cfg$paths$preflight, "score2_candidate_source_qc.csv")
  manifest_file <- file.path(cfg$paths$preflight, "input_manifest.csv")
  if (!all(file.exists(c(summary_file, source_qc_file, manifest_file)))) return(FALSE)
  summary <- read_json(summary_file, simplifyVector = TRUE)
  source_qc <- fread(source_qc_file)
  manifest <- fread(manifest_file)
  field_hash <- manifest[role == "field_dictionary", sha256]
  identical(summary$score2_source_qc_status, "PASS") &&
    identical(summary$score2_total_cholesterol_source, "bb_TC") &&
    identical(summary$score2_hdl_source, "bb_HDL") &&
    nrow(source_qc[selected == TRUE]) == 4L &&
    all(source_qc[selected == TRUE, status] == "PASS") &&
    length(field_hash) == 1L && identical(field_hash, yur_sha_file(cfg$field_map_file))
}

yur_preflight <- function(cfg) {
  required <- c("data.table", "jsonlite", "digest", "survival", "pROC", "ggplot2", "bit64")
  pkg <- data.table(package = required, available = vapply(required, requireNamespace, logical(1), quietly = TRUE))
  yur_write_csv(pkg, file.path(cfg$paths$preflight, "r_package_check.csv"))
  if (any(!pkg$available)) {
    stop("Missing R packages: ", paste(pkg[available == FALSE, package], collapse = ", "))
  }

  inputs <- data.table(
    role = c("phenotype", "raw_phenotype", "raw_protein", "panel_mapping",
             "supplement_workbook", "supplement_methods", "olink_processing_dates",
             "outcome_dictionary", "field_dictionary", "method_provenance"),
    path = c(cfg$phenotype_rds, cfg$raw_phenotype_file, cfg$raw_protein_file, cfg$panel_mapping_file,
             cfg$supplement_workbook_file, cfg$supplement_methods_file,
             cfg$olink_processing_start_date_file, cfg$outcomes_file, cfg$field_map_file, cfg$method_provenance_file)
  )
  inputs[, exists := file.exists(path)]
  inputs[, bytes := fifelse(exists, file.info(path)$size, NA_real_)]
  inputs[, identity_method := fifelse(role == "raw_phenotype", "path_size_mtime_and_prefix_sha256", "full_sha256")]
  inputs[, sha256 := NA_character_]
  inputs[exists & role != "raw_phenotype", sha256 := vapply(path, yur_sha_file, character(1))]
  inputs[exists & role == "raw_phenotype", sha256 := yur_supplemental_source_identity(path[[1]])]
  yur_write_csv(inputs, file.path(cfg$paths$preflight, "input_manifest.csv"))
  if (any(!inputs$exists)) {
    stop("Missing required inputs: ", paste(inputs[exists == FALSE, path], collapse = ", "))
  }

  raw_header <- names(fread(cfg$raw_protein_file, nrows = 0, showProgress = FALSE))
  eid_hits <- intersect(c("eid", "id", "f.eid", "participant_id"), raw_header)
  if (!length(eid_hits)) stop("No EID column found in raw protein table.")
  raw_features <- setdiff(raw_header, eid_hits[[1]])
  if (length(raw_features) < 1000L) stop("Raw protein table has only ", length(raw_features), " feature columns; full Explore panel required.")
  protein_ids <- fread(
    cfg$raw_protein_file, select = eid_hits[[1]], showProgress = FALSE, nThread = cfg$workers
  )[[eid_hits[[1]]]]
  protein_ids <- unique(yur_norm_eid(protein_ids))

  ph <- readRDS(cfg$phenotype_rds)
  supplemental <- yur_prepare_supplemental_covariates(cfg)
  field_resolution <- yur_resolve_fields(names(ph), cfg$field_map_file)
  field_resolution[, `:=`(source_type = "phenotype_rds", derivation = "direct")]
  baseline_source <- field_resolution[canonical == "baseline_date", resolved]
  plate_source <- field_resolution[canonical == "protein_plate", resolved]
  if (length(baseline_source) && !is.na(baseline_source)) {
    field_resolution[canonical == "blood_collection_season", `:=`(
      resolved = baseline_source, status = "PASS", source_type = "derived",
      derivation = "meteorological season from baseline assessment date"
    )]
  }
  field_resolution[canonical == "fasting_time", `:=`(
    resolved = "p74_i0", status = "PASS", source_type = "supplemental_raw_phenotype",
    derivation = "UKB field 74 baseline participant-reported fasting time"
  )]
  processing <- yur_read_processing_dates(cfg$olink_processing_start_date_file)
  if (length(plate_source) && !is.na(plate_source) && nrow(processing)) {
    field_resolution[canonical == "protein_sampling_lag_days", `:=`(
      resolved = paste(plate_source, "Panel", "Processing_StartDate", sep = "+"),
      status = "PASS", source_type = "derived_panel_specific",
      derivation = "panel-specific Olink processing start date minus baseline blood collection date"
    )]
  }
  yur_write_csv(field_resolution, file.path(cfg$paths$preflight, "local_field_resolution.csv"))
  if (any(field_resolution$status == "FAIL")) {
    stop("Required phenotype fields unresolved: ", paste(field_resolution[status == "FAIL", canonical], collapse = ", "))
  }

  ph_eid_source <- field_resolution[canonical == "eid", resolved]
  ph_eid <- yur_norm_eid(ph[[ph_eid_source]])
  ph_index <- match(protein_ids, ph_eid)
  score2_source_qc <- yur_score2_source_qc(ph, ph_index, field_resolution, cfg)
  yur_write_csv(score2_source_qc, file.path(cfg$paths$preflight, "score2_candidate_source_qc.csv"))
  selected_score2_qc <- score2_source_qc[selected == TRUE]
  if (nrow(selected_score2_qc) != 4L || any(selected_score2_qc$status != "PASS")) {
    bad <- selected_score2_qc[status != "PASS"]
    stop(
      "SCORE2 source QC failed; expected plausible age, SBP, total cholesterol and HDL. ",
      paste(sprintf("%s=%s median=%s", bad$canonical, bad$source, format(bad$median, digits = 5)), collapse = "; "),
      call. = FALSE
    )
  }
  fasting_match <- supplemental[match(protein_ids, eid), fasting_time]
  baseline_values <- ph[[baseline_source]][ph_index]
  plate_values <- if (length(plate_source) && !is.na(plate_source)) {
    yur_plate_character(ph[[plate_source]][ph_index])
  } else {
    rep(NA_character_, length(protein_ids))
  }
  matched_plates <- unique(processing$plate_id)
  supplemental_qc <- data.table(
    scope = "raw_protein_participants",
    variable = c("blood_collection_season", "fasting_time", "protein_plate", "protein_sampling_lag_days"),
    source = c(baseline_source, "pheno.tsv.gz:p74_i0", plate_source, "UKB resource 1019 panel processing dates"),
    nonmissing_n = c(
      sum(!is.na(yur_blood_collection_season(baseline_values))),
      sum(!is.na(fasting_match)), sum(!is.na(plate_values)), sum(!is.na(plate_values) & plate_values %in% matched_plates)
    ),
    total_n = length(protein_ids)
  )
  supplemental_qc[, missing_rate := 1 - nonmissing_n / total_n]
  yur_write_csv(supplemental_qc, file.path(cfg$paths$preflight, "supplemental_covariate_qc.csv"))
  provenance <- field_resolution[canonical %in% c(
    "blood_collection_season", "fasting_time", "protein_plate", "protein_sampling_lag_days"
  ), .(variable = canonical, source_type, source_field = resolved, derivation, status)]
  provenance[, official_reference := fifelse(
    variable == "protein_sampling_lag_days",
    "UK Biobank Resource 1019: olink_processing_start_date.dat",
    fifelse(variable == "fasting_time", "UK Biobank Field 74", "local baseline phenotype")
  )]
  yur_write_csv(provenance, file.path(cfg$paths$preflight, "covariate_provenance.csv"))

  outcomes <- fread(cfg$outcomes_file)
  required_outcome_columns <- c(
    "outcome_id", "outcome_label", "local_date_candidates", "death_icd10_prefixes",
    "source_policy", "min_derivation_events", "min_test_events", "min_total_events"
  )
  missing_outcome_columns <- setdiff(required_outcome_columns, names(outcomes))
  if (length(missing_outcome_columns)) {
    stop("Outcome dictionary is missing columns: ", paste(missing_outcome_columns, collapse = ", "))
  }
  if (anyDuplicated(outcomes$outcome_id)) stop("Outcome dictionary contains duplicate outcome_id values.")
  death_required <- any(nzchar(trimws(outcomes$death_icd10_prefixes)))
  death_fields <- c("death_icd10", "death_icd10sec", "date_death")
  if (death_required && any(!death_fields %in% names(ph))) {
    stop(
      "S24 source-locked death linkage requires phenotype fields: ",
      paste(setdiff(death_fields, names(ph)), collapse = ", ")
    )
  }
  outcome_resolution <- outcomes[, {
    candidates <- strsplit(local_date_candidates, ";", fixed = TRUE)[[1]]
    hit <- candidates[candidates %in% names(ph)]
    .(resolved_date_fields = paste(hit, collapse = ";"), resolved_n = length(hit),
      status = if (length(hit)) "PASS" else "FAIL")
  }, by = .(
    outcome_id, outcome_label, local_date_candidates, death_icd10_prefixes,
    source_policy, min_derivation_events, min_test_events, min_total_events
  )]
  outcome_resolution[, definition_hash := vapply(seq_len(.N), function(i) {
    yur_sha_text(c(
      outcome_id[[i]], resolved_date_fields[[i]], death_icd10_prefixes[[i]], source_policy[[i]],
      min_derivation_events[[i]], min_test_events[[i]], min_total_events[[i]]
    ))
  }, character(1))]
  yur_write_csv(outcome_resolution, file.path(cfg$paths$preflight, "outcome_field_resolution.csv"))
  if (any(outcome_resolution$status == "FAIL")) {
    stop("Outcome fields unresolved: ", paste(outcome_resolution[status == "FAIL", outcome_id], collapse = ", "))
  }

  mapping <- fread(cfg$panel_mapping_file)
  yur_write_json(list(
    status = "PASS",
    raw_participant_id_column = eid_hits[[1]],
    raw_participant_n = length(protein_ids),
    raw_feature_columns = length(raw_features),
    raw_duplicate_column_names = sum(duplicated(raw_header)),
    local_mapping_rows = nrow(mapping),
    supplemental_fasting_nonmissing_in_raw_protein_participants = sum(!is.na(fasting_match)),
    olink_processing_plate_n = uniqueN(processing$plate_id),
    olink_processing_panel_n = uniqueN(processing$panel_key),
    phenotype_plate_match_n_in_raw_protein_participants = sum(!is.na(plate_values) & plate_values %in% matched_plates),
    outcome_count = nrow(outcomes),
    selected_outcome_count = nrow(yur_select_outcomes(cfg)),
    score2_total_cholesterol_source = selected_score2_qc[canonical == "total_cholesterol", source],
    score2_hdl_source = selected_score2_qc[canonical == "hdl", source],
    score2_source_qc_status = "PASS",
    followup_cutoff = cfg$followup_cutoff,
    split_seed = cfg$split_seed,
    warning = "The publication reported 2,923 raw and 2,920 retained proteins; local counts are audited and never forced."
  ), file.path(cfg$paths$preflight, "preflight_summary.json"))
}
