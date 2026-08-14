yur_cmr_metrics <- function(cfg) {
  metrics <- fread(cfg$cmr_metrics_file)
  required <- c("metric_id", "Outcome", "anatomy", "source_column", "source_field")
  missing <- setdiff(required, names(metrics))
  if (length(missing)) stop("CMR metric dictionary lacks: ", paste(missing, collapse = ", "))
  if (nrow(metrics) != 19L || anyDuplicated(metrics$metric_id) || anyDuplicated(metrics$Outcome)) {
    stop("CMR metric dictionary must contain 19 unique metric IDs and labels.")
  }
  metrics
}

yur_cmr_requested_metrics <- function(cfg) {
  metrics <- yur_cmr_metrics(cfg)
  requested <- trimws(strsplit(as.character(cfg$cmr_metric_subset), ",", fixed = TRUE)[[1]])
  if (length(requested) == 1L && tolower(requested) == "all") return(metrics)
  unknown <- setdiff(requested, metrics$metric_id)
  if (length(unknown)) stop("Unknown CMR metric IDs: ", paste(unknown, collapse = ", "))
  metrics[metric_id %in% requested]
}

yur_cmr_prepare_complete <- function(cfg) {
  input_file <- file.path(cfg$paths$cmr, "cmr_analysis_input.rds")
  manifest_file <- file.path(cfg$paths$cmr, "cmr_prepare_manifest.json")
  if (!file.exists(input_file) || !file.exists(manifest_file)) return(FALSE)
  manifest <- tryCatch(read_json(manifest_file, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(manifest) || !identical(manifest$status, "PASS")) return(FALSE)
  identical(as.character(manifest$retained_panel_hash),
            readLines(file.path(cfg$paths$cox, "retained_panel_hash.txt"), warn = FALSE)[[1]]) &&
    identical(as.character(manifest$cmr_source_sha256), yur_sha_file(cfg$cmr_feature_file))
}

yur_prepare_cmr <- function(cfg) {
  retained_file <- file.path(cfg$paths$cox, "retained_panel.csv")
  panel_hash_file <- file.path(cfg$paths$cox, "retained_panel_hash.txt")
  if (!file.exists(retained_file) || !file.exists(panel_hash_file)) {
    stop("Run cox_prepare first so CMR and incident analyses share the retained protein panel.")
  }
  if (!file.exists(cfg$cmr_feature_file)) stop("CMR feature file is missing: ", cfg$cmr_feature_file)
  metrics <- yur_cmr_metrics(cfg)
  cmr_header <- names(fread(cfg$cmr_feature_file, nrows = 0, showProgress = FALSE))
  eid_col <- intersect(c("participant_id", "eid", "id"), cmr_header)
  if (!length(eid_col)) stop("No participant ID column found in CMR feature file.")
  missing_metrics <- setdiff(metrics$source_column, cmr_header)
  if (length(missing_metrics)) stop("CMR source lacks metrics: ", paste(missing_metrics, collapse = ", "))

  cmr <- fread(
    cfg$cmr_feature_file, select = c(eid_col[[1]], metrics$source_column),
    showProgress = TRUE, nThread = cfg$workers
  )
  setnames(cmr, eid_col[[1]], "eid")
  cmr <- yur_normalize_eid_column(cmr)
  cmr <- unique(cmr[!is.na(eid)], by = "eid")
  for (v in metrics$source_column) set(cmr, j = v, value = suppressWarnings(as.numeric(cmr[[v]])))

  cohorts <- yur_load_cox_cohorts(cfg)
  meta <- copy(cohorts$all_meta)
  meta <- yur_normalize_eid_column(meta)
  cmr_meta <- merge(meta, cmr, by = "eid", all = FALSE, sort = FALSE)
  metric_columns <- metrics$source_column
  complete_metric_n <- rowSums(is.finite(as.matrix(cmr_meta[, ..metric_columns])))
  cmr_meta <- cmr_meta[complete_metric_n > 0L]
  if (nrow(cmr_meta) < 3000L) stop("Only ", nrow(cmr_meta), " CMR participants overlap the baseline-CVD-free proteomics cohort.")

  retained <- fread(retained_file)
  keep <- retained$feature_id
  expected_hash <- readLines(panel_hash_file, warn = FALSE)[[1]]
  if (!identical(yur_sha_text(keep), expected_hash)) stop("Retained protein panel hash mismatch.")
  protein_header <- names(fread(cfg$raw_protein_file, nrows = 0, showProgress = FALSE))
  protein_eid <- intersect(c("eid", "id", "f.eid", "participant_id"), protein_header)
  if (!length(protein_eid)) stop("No EID column found in raw protein table.")
  absent <- setdiff(keep, protein_header)
  if (length(absent)) stop("Raw protein table lacks ", length(absent), " retained proteins.")
  raw <- fread(
    cfg$raw_protein_file, select = c(protein_eid[[1]], keep),
    showProgress = TRUE, nThread = cfg$workers
  )
  setnames(raw, protein_eid[[1]], "eid")
  raw <- yur_normalize_eid_column(raw)
  setkey(raw, eid)
  missing_eids <- setdiff(cmr_meta$eid, raw$eid)
  if (length(missing_eids)) stop("Raw protein table lacks ", length(missing_eids), " CMR cohort EIDs.")
  proteins <- raw[J(cmr_meta$eid), ..keep]
  if (nrow(proteins) != nrow(cmr_meta)) stop("CMR protein alignment failed.")

  metric_qc <- rbindlist(lapply(seq_len(nrow(metrics)), function(i) {
    values <- cmr_meta[[metrics$source_column[[i]]]]
    finite <- values[is.finite(values)]
    data.table(
      metric_id = metrics$metric_id[[i]], Outcome = metrics$Outcome[[i]],
      anatomy = metrics$anatomy[[i]], source_column = metrics$source_column[[i]],
      source_field = metrics$source_field[[i]], n = length(finite),
      missing_rate = 1 - length(finite) / nrow(cmr_meta), mean = mean(finite),
      sd = sd(finite), q01 = unname(quantile(finite, .01)),
      median = median(finite), q99 = unname(quantile(finite, .99)),
      status = if (length(finite) >= 3000L && is.finite(sd(finite)) && sd(finite) > 0) "PASS" else "FAIL"
    )
  }))
  yur_write_csv(metric_qc, file.path(cfg$paths$cmr, "cmr_metric_qc.csv"))
  if (any(metric_qc$status != "PASS")) {
    stop("CMR metric QC failed: ", paste(metric_qc[status != "PASS", metric_id], collapse = ", "))
  }

  input <- list(
    meta = cmr_meta,
    proteins = as.data.frame(proteins),
    mapping = retained[, .(feature_id, protein, panel, olink_id, mapping_status)],
    metrics = metrics,
    retained_panel_hash = expected_hash
  )
  saveRDS(input, file.path(cfg$paths$cmr, "cmr_analysis_input.rds"), compress = FALSE)
  fwrite(cmr_meta[, c("eid", metrics$source_column), with = FALSE],
         file.path(cfg$paths$cmr, "cmr_phenotypes.csv.gz"), na = "")
  yur_write_csv(data.table(eid = cmr_meta$eid), file.path(cfg$paths$cmr, "cmr_cohort_eid.csv"))
  yur_write_json(list(
    status = "PASS", participants = nrow(cmr_meta), metrics = nrow(metrics),
    retained_proteins = length(keep), retained_panel_hash = expected_hash,
    cmr_source = cfg$cmr_feature_file, cmr_source_sha256 = yur_sha_file(cfg$cmr_feature_file),
    cohort_definition = "UKB-PPP participants free of all 14 CVDs at baseline with at least one CMR metric",
    generated = yur_now()
  ), file.path(cfg$paths$cmr, "cmr_prepare_manifest.json"))
  yur_log(cfg, "CMR prepared participants=", nrow(cmr_meta), " metrics=", nrow(metrics),
          " proteins=", length(keep))
}

yur_linear_one <- function(protein, outcome, covariates) {
  cov_ok <- rowSums(!is.finite(covariates)) == 0L
  ok <- is.finite(protein) & is.finite(outcome) & cov_ok
  n <- sum(ok)
  if (n < 100L) return(c(n = n, beta = NA, se = NA, t = NA, p = NA))
  protein_mean <- mean(protein[ok])
  protein_sd <- sd(protein[ok])
  outcome_mean <- mean(outcome[ok])
  outcome_sd <- sd(outcome[ok])
  if (!is.finite(protein_sd) || protein_sd <= 0 || !is.finite(outcome_sd) || outcome_sd <= 0) {
    return(c(n = n, beta = NA, se = NA, t = NA, p = NA))
  }
  z_protein <- (protein[ok] - protein_mean) / protein_sd
  z_outcome <- (outcome[ok] - outcome_mean) / outcome_sd
  design <- cbind(`(Intercept)` = 1, protein = z_protein, covariates[ok, , drop = FALSE])
  fit <- tryCatch(lm.fit(design, z_outcome), error = function(e) NULL)
  if (is.null(fit) || fit$rank < 2L || !is.finite(fit$coefficients[[2]])) {
    return(c(n = n, beta = NA, se = NA, t = NA, p = NA))
  }
  df <- n - fit$rank
  if (df <= 0L) return(c(n = n, beta = NA, se = NA, t = NA, p = NA))
  r <- qr.R(fit$qr)[seq_len(fit$rank), seq_len(fit$rank), drop = FALSE]
  unscaled <- tryCatch(chol2inv(r), error = function(e) NULL)
  if (is.null(unscaled) || !is.finite(unscaled[2, 2])) {
    return(c(n = n, beta = NA, se = NA, t = NA, p = NA))
  }
  beta <- unname(fit$coefficients[[2]])
  se <- sqrt(sum(fit$residuals^2) / df * unscaled[2, 2])
  statistic <- beta / se
  c(n = n, beta = beta, se = se, t = statistic,
    p = 2 * pt(abs(statistic), df = df, lower.tail = FALSE))
}

yur_cmr_contract <- function(cfg, metric_id) {
  prepare <- read_json(file.path(cfg$paths$cmr, "cmr_prepare_manifest.json"), simplifyVector = TRUE)
  list(
    metric_id = metric_id,
    retained_panel_hash = as.character(prepare$retained_panel_hash),
    cmr_source_sha256 = as.character(prepare$cmr_source_sha256),
    participants = as.integer(prepare$participants)
  )
}

yur_cmr_shard_complete <- function(cfg) {
  metrics <- yur_cmr_requested_metrics(cfg)
  all(vapply(metrics$metric_id, function(metric_id) {
    out <- file.path(cfg$paths$cmr, paste0("cmr_", metric_id, ".csv.gz"))
    contract <- file.path(cfg$paths$cmr, paste0("cmr_", metric_id, ".contract.json"))
    if (!file.exists(out) || !file.exists(contract)) return(FALSE)
    observed <- tryCatch(fread(out, select = "feature_id"), error = function(e) NULL)
    expected_n <- nrow(fread(file.path(cfg$paths$cox, "retained_panel.csv")))
    !is.null(observed) && nrow(observed) == expected_n &&
      yur_cox_contract_matches(contract, yur_cmr_contract(cfg, metric_id))
  }, logical(1)))
}

yur_run_cmr_shard <- function(cfg) {
  input_file <- file.path(cfg$paths$cmr, "cmr_analysis_input.rds")
  if (!file.exists(input_file)) stop("Run cmr_prepare first.")
  input <- readRDS(input_file)
  metrics <- yur_cmr_requested_metrics(cfg)
  mapping <- copy(input$mapping)
  mapping[, panel_key := yur_panel_key(panel)]
  mapping[is.na(panel_key), panel_key := ""]
  panel_keys <- unique(mapping$panel_key)
  covariates_by_panel <- setNames(lapply(panel_keys, function(key) {
    yur_covariate_matrix(input$meta, protein_panel = key)
  }), panel_keys)
  protein_n <- ncol(input$proteins)
  threshold <- .05 / protein_n

  for (i in seq_len(nrow(metrics))) {
    metric <- metrics[i]
    out_file <- file.path(cfg$paths$cmr, paste0("cmr_", metric$metric_id, ".csv.gz"))
    contract_file <- file.path(cfg$paths$cmr, paste0("cmr_", metric$metric_id, ".contract.json"))
    expected <- yur_cmr_contract(cfg, metric$metric_id)
    if (cfg$resume && file.exists(out_file) && file.exists(contract_file) && !cfg$force &&
        yur_cox_contract_matches(contract_file, expected)) {
      yur_log(cfg, "Resume: reuse ", basename(out_file))
      next
    }
    outcome <- input$meta[[metric$source_column]]
    outcome_all <- outcome[is.finite(outcome)]
    result <- vector("list", protein_n)
    yur_log(cfg, "CMR metric=", metric$metric_id, " n=", length(outcome_all), " proteins=", protein_n)
    for (j in seq_len(protein_n)) {
      current_feature <- names(input$proteins)[[j]]
      panel_hit <- mapping[feature_id == current_feature, panel_key]
      panel_key <- if (length(panel_hit)) panel_hit[[1]] %||% "" else ""
      covariates <- covariates_by_panel[[panel_key]]
      stat <- yur_linear_one(input$proteins[[j]], outcome, covariates)
      result[[j]] <- cbind(data.table(feature_id = current_feature), as.data.table(as.list(stat)))
      if (j %% 250L == 0L) yur_log(cfg, "CMR progress ", metric$metric_id, " ", j, "/", protein_n)
    }
    result <- rbindlist(result)
    result <- merge(result, mapping[, .(feature_id, protein, panel, olink_id, mapping_status)],
                    by = "feature_id", all.x = TRUE, sort = FALSE)
    result[, `:=`(
      metric_id = metric$metric_id, Outcome = metric$Outcome, anatomy = metric$anatomy,
      source_column = metric$source_column, source_field = metric$source_field,
      effect_scale = "SD_CMR_per_SD_protein", bonferroni_threshold = threshold,
      bonferroni_significant = is.finite(p) & p < threshold
    )]
    setnames(result, "protein", "Protein")
    setcolorder(result, c(
      "metric_id", "Outcome", "anatomy", "source_column", "source_field",
      "feature_id", "Protein", "panel", "olink_id", "effect_scale",
      "n", "beta", "se", "t", "p", "bonferroni_threshold",
      "bonferroni_significant", "mapping_status"
    ))
    fwrite(result, out_file, na = "")
    expected$rows <- nrow(result)
    expected$generated <- yur_now()
    yur_write_json(expected, contract_file)
  }
}

yur_merge_cmr <- function(cfg) {
  metrics <- yur_cmr_metrics(cfg)
  retained <- fread(file.path(cfg$paths$cox, "retained_panel.csv"))
  rows <- lapply(metrics$metric_id, function(metric_id) {
    path <- file.path(cfg$paths$cmr, paste0("cmr_", metric_id, ".csv.gz"))
    contract <- file.path(cfg$paths$cmr, paste0("cmr_", metric_id, ".contract.json"))
    if (!file.exists(path) || !file.exists(contract)) stop("Missing CMR shard: ", metric_id)
    expected <- yur_cmr_contract(cfg, metric_id)
    if (!yur_cox_contract_matches(contract, expected)) stop("Stale CMR shard contract: ", metric_id)
    x <- fread(path)
    if (nrow(x) != nrow(retained) || !identical(sort(x$feature_id), sort(retained$feature_id))) {
      stop("CMR shard panel mismatch: ", metric_id)
    }
    x
  })
  result <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  fwrite(result, file.path(cfg$paths$cmr, "cmr_associations.csv.gz"), na = "")
  summary <- result[, .(
    tested = .N, significant = sum(bonferroni_significant, na.rm = TRUE),
    positive = sum(bonferroni_significant & beta > 0, na.rm = TRUE),
    negative = sum(bonferroni_significant & beta < 0, na.rm = TRUE),
    unique_significant_proteins = uniqueN(feature_id[bonferroni_significant])
  ), by = .(metric_id, Outcome, anatomy)]
  yur_write_csv(summary, file.path(cfg$paths$cmr, "cmr_association_summary.csv"))
  prepare <- read_json(file.path(cfg$paths$cmr, "cmr_prepare_manifest.json"), simplifyVector = TRUE)
  yur_write_json(list(
    status = "PASS", participants = prepare$participants, metrics = nrow(metrics),
    retained_proteins = nrow(retained), rows = nrow(result),
    significant_associations = sum(result$bonferroni_significant, na.rm = TRUE),
    significant_metrics = uniqueN(result[bonferroni_significant == TRUE, metric_id]),
    unique_significant_proteins = uniqueN(result[bonferroni_significant == TRUE, feature_id]),
    retained_panel_hash = prepare$retained_panel_hash, generated = yur_now()
  ), file.path(cfg$paths$cmr, "cmr_summary.json"))
}

yur_run_cmr <- function(cfg) {
  yur_prepare_cmr(cfg)
  yur_run_cmr_shard(cfg)
  yur_merge_cmr(cfg)
}
