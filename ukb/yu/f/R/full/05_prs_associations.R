yur_prs_thresholds <- function(cfg) {
  x <- fread(file.path(cfg$project_dir, "f", "config", "prs_thresholds.tsv"))
  x[, threshold := as.numeric(threshold)]
  if (nrow(x) != 5L || anyDuplicated(x$threshold) || anyDuplicated(x$score_column)) {
    stop("PRS threshold contract must contain five unique thresholds and score columns.", call. = FALSE)
  }
  setorder(x, threshold)
  x
}

yur_prs_sources <- function(cfg) {
  x <- fread(file.path(cfg$project_dir, "f", "config", "prs_gwas_sources.tsv"))
  required <- c("outcome_id", "outcome_label", "source_kind", "source_id", "build", "url", "provenance", "status")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("PRS GWAS manifest is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (nrow(x) != 13L || anyDuplicated(x$outcome_id)) {
    stop("PRS GWAS manifest must contain 13 unique outcome rows.", call. = FALSE)
  }
  x
}

yur_prs_weight_policy <- function(cfg) {
  path <- file.path(cfg$project_dir, "f", "config", "prs_weight_source_policy.tsv")
  if (!file.exists(path)) stop("PRS weight-source policy is missing: ", path, call. = FALSE)
  x <- fread(path)
  required <- c("outcome_id", "weight_basis", "author_weight_file_public", "interpretation")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("PRS weight-source policy is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (nrow(x) != 13L || anyDuplicated(x$outcome_id)) {
    stop("PRS weight-source policy must contain 13 unique outcome rows.", call. = FALSE)
  }
  x
}

yur_prs_candidate_table <- function(cfg) {
  path <- file.path(cfg$paths$cox, "table_s2_incident_associations.csv.gz")
  if (!file.exists(path)) stop("Local full-cohort Cox table is missing: ", path, call. = FALSE)
  sources <- yur_prs_sources(cfg)
  x <- fread(path)
  required <- c("scope", "outcome_id", "feature_id", "protein", "panel", "bonferroni_significant")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Local Cox table is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  x <- x[scope == "full_incident" & outcome_id %chin% sources$outcome_id & bonferroni_significant == TRUE]
  x <- unique(x[, .(outcome_id, feature_id, protein, panel, olink_id)])
  if (!nrow(x)) stop("No local outcome-specific Bonferroni proteins are available for PRS analysis.", call. = FALSE)
  x[, candidate_n := .N, by = outcome_id]
  x
}

yur_prs_target_cohort <- function(cfg) {
  files <- file.path(cfg$paths$cohort, c("derivation_cohort.csv.gz", "test_cohort.csv.gz"))
  if (any(!file.exists(files))) stop("Yu cohort files are missing; run the cohort stage first.", call. = FALSE)
  x <- rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE)
  x <- yur_normalize_eid_column(x)
  if (anyDuplicated(x$eid)) stop("Yu target cohort has duplicated EIDs.", call. = FALSE)
  x
}

yur_prs_preflight <- function(cfg) {
  sources <- yur_prs_sources(cfg)
  weight_policy <- yur_prs_weight_policy(cfg)
  thresholds <- yur_prs_thresholds(cfg)
  cohort <- yur_prs_target_cohort(cfg)
  candidates <- yur_prs_candidate_table(cfg)
  required_files <- data.table(
    role = c("raw_protein", "phenotype_rds", "cox_results", "retained_panel", "split_eid"),
    path = c(
      cfg$raw_protein_file, cfg$phenotype_rds,
      file.path(cfg$paths$cox, "table_s2_incident_associations.csv.gz"),
      file.path(cfg$paths$cox, "retained_panel.csv"),
      file.path(cfg$paths$cohort, "split_eid.csv")
    )
  )
  required_files[, exists := file.exists(path)]
  if (any(!required_files$exists)) {
    stop("PRS preflight missing files: ", paste(required_files[exists == FALSE, path], collapse = "; "), call. = FALSE)
  }

  raw_header <- names(fread(cfg$raw_protein_file, nrows = 0, showProgress = FALSE))
  eid_col <- intersect(c("eid", "id", "f.eid", "participant_id"), raw_header)
  if (length(eid_col) != 1L) stop("Raw protein table must have exactly one resolvable EID column.", call. = FALSE)
  missing_features <- setdiff(unique(candidates$feature_id), raw_header)
  if (length(missing_features)) {
    stop("Raw protein table is missing ", length(missing_features), " PRS candidate proteins; first=",
         paste(head(missing_features, 10L), collapse = ", "), call. = FALSE)
  }

  keep_dir <- file.path(cfg$paths$prs, "inputs")
  dir.create(keep_dir, recursive = TRUE, showWarnings = FALSE)
  keep_file <- file.path(keep_dir, "target.keep.tsv")
  fwrite(cohort[, .(`#FID` = eid, IID = eid)], keep_file, sep = "\t", quote = FALSE)
  writeLines(yur_sha_text(sort(cohort$eid)), file.path(keep_dir, "target_eid_hash.txt"))

  expected <- candidates[, .(candidate_proteins = uniqueN(feature_id)), by = outcome_id]
  expected <- merge(sources[, .(outcome_id, outcome_label, source_id, provenance, status)], expected,
                    by = "outcome_id", all.x = TRUE)
  expected <- merge(expected, weight_policy, by = "outcome_id", all.x = TRUE, sort = FALSE)
  expected[is.na(candidate_proteins), candidate_proteins := 0L]
  expected[, expected_regression_rows := candidate_proteins * nrow(thresholds)]
  yur_write_csv(expected, file.path(cfg$paths$prs, "expected_association_rows.csv"))
  yur_write_csv(required_files, file.path(cfg$paths$prs, "input_file_qc.csv"))
  yur_write_json(list(
    status = "PASS", target_participants = nrow(cohort), outcomes = nrow(sources),
    thresholds = as.character(thresholds$label), candidate_pairs = nrow(candidates),
    expected_complete_rows = sum(expected$expected_regression_rows),
    source_manifest_hash = yur_sha_file(file.path(cfg$project_dir, "f", "config", "prs_gwas_sources.tsv")),
    weight_policy_hash = yur_sha_file(file.path(cfg$project_dir, "f", "config", "prs_weight_source_policy.tsv")),
    threshold_manifest_hash = yur_sha_file(file.path(cfg$project_dir, "f", "config", "prs_thresholds.tsv")),
    target_eid_hash = yur_sha_text(sort(cohort$eid)),
    reconstruction_note = paste(
      "GWAS effect alleles and beta values are used as reference weights.",
      "The article did not publish final per-SNP PRS files or its complete clumping command.",
      "This is therefore an article-source reconstruction with frozen PRSice-1.25-compatible clumping."
    )
  ), file.path(cfg$paths$prs, "prs_preflight_summary.json"))
  invisible(expected)
}

yur_prs_score_column <- function(names_available, score_column) {
  normalized <- yur_norm_name(names_available)
  target <- yur_norm_name(paste0(score_column, "_SUM"))
  hit <- names_available[normalized == target]
  if (!length(hit)) {
    target <- yur_norm_name(score_column)
    hit <- names_available[normalized == target]
  }
  if (length(hit) != 1L) {
    stop("Could not resolve exactly one score column for ", score_column,
         "; available=", paste(names_available, collapse = ", "), call. = FALSE)
  }
  hit[[1]]
}

yur_prs_merge_scores <- function(cfg) {
  sources <- yur_prs_sources(cfg)
  thresholds <- yur_prs_thresholds(cfg)
  cohort <- yur_prs_target_cohort(cfg)
  reference_outcome <- sources$outcome_id[[1L]]
  reference_path <- file.path(
    cfg$paths$prs, "scores", reference_outcome, "chr1",
    paste0(reference_outcome, "_chr1.sscore")
  )
  if (!file.exists(reference_path)) stop("Missing reference PRS score: ", reference_path, call. = FALSE)
  reference <- fread(reference_path, showProgress = FALSE)
  reference_iid <- intersect(c("#IID", "IID"), names(reference))
  if (length(reference_iid) != 1L) stop("Cannot resolve reference PRS IID: ", reference_path, call. = FALSE)
  reference[, eid := yur_norm_eid(get(reference_iid))]
  if (anyDuplicated(reference$eid)) stop("Reference PRS score has duplicated EIDs.", call. = FALSE)

  target <- cohort[eid %chin% reference$eid, .(eid)]
  excluded <- cohort[!eid %chin% reference$eid, .(eid)]
  excluded[, reason := "missing_genotype_score"]
  if (!nrow(target)) stop("No Yu target participants have genotype scores.", call. = FALSE)
  if (nrow(reference) != nrow(target) || !setequal(reference$eid, target$eid)) {
    stop("Reference PRS sample set differs from the genotype-available Yu target cohort.", call. = FALSE)
  }
  yur_write_csv(target, file.path(cfg$paths$prs, "prs_genotype_available_eid.csv"))
  yur_write_csv(excluded, file.path(cfg$paths$prs, "prs_excluded_missing_genotype_eid.csv"))
  yur_write_csv(data.table(
    requested_target_participants = nrow(cohort),
    genotype_available_participants = nrow(target),
    excluded_missing_genotype = nrow(excluded),
    availability_rate = nrow(target) / nrow(cohort),
    requested_target_eid_hash = yur_sha_text(sort(cohort$eid)),
    genotype_available_eid_hash = yur_sha_text(sort(target$eid))
  ), file.path(cfg$paths$prs, "prs_sample_qc.csv"))
  variant_path <- file.path(cfg$paths$prs, "score_variant_counts.tsv")
  variants <- if (file.exists(variant_path)) fread(variant_path) else NULL
  if (!is.null(variants)) variants <- variants[, .SD[.N], by = .(outcome_id, chr, threshold, score_column)]

  rows <- vector("list", nrow(sources) * nrow(thresholds))
  row_index <- 0L
  merge_qc <- list()
  for (source_index in seq_len(nrow(sources))) {
    outcome_id <- sources$outcome_id[[source_index]]
    source_outcome <- outcome_id
    sums <- setNames(lapply(thresholds$score_column, function(x) numeric(nrow(target))), thresholds$score_column)
    for (chr in 1:22) {
      path <- file.path(cfg$paths$prs, "scores", outcome_id, paste0("chr", chr),
                        paste0(outcome_id, "_chr", chr, ".sscore"))
      if (!file.exists(path)) stop("Missing PRS chromosome score: ", path, call. = FALSE)
      score <- fread(path, showProgress = FALSE)
      iid_col <- intersect(c("#IID", "IID"), names(score))
      if (length(iid_col) != 1L) stop("Cannot resolve IID in ", path, call. = FALSE)
      score[, eid := yur_norm_eid(get(iid_col))]
      if (anyDuplicated(score$eid)) stop("Duplicated score EIDs in ", path, call. = FALSE)
      if (nrow(score) != nrow(target) || !setequal(score$eid, target$eid)) {
        stop("Inconsistent genotype-available sample set in ", path, call. = FALSE)
      }
      index <- match(target$eid, score$eid)
      if (anyNA(index)) stop("PRS score is missing target EIDs in ", path, call. = FALSE)
      for (score_name in thresholds$score_column) {
        observed <- yur_prs_score_column(names(score), score_name)
        value <- as.numeric(score[[observed]][index])
        if (any(!is.finite(value))) stop("Non-finite PRS chromosome score in ", path, call. = FALSE)
        sums[[score_name]] <- sums[[score_name]] + value
      }
    }
    for (threshold_index in seq_len(nrow(thresholds))) {
      row_index <- row_index + 1L
      score_name <- thresholds$score_column[[threshold_index]]
      raw <- sums[[score_name]]
      score_sd <- sd(raw)
      if (!is.finite(score_sd) || score_sd <= 0) stop("Zero-variance PRS: ", outcome_id, "/", score_name, call. = FALSE)
      variant_count <- NA_integer_
      if (!is.null(variants)) {
        hit <- variants[outcome_id == source_outcome &
                          abs(as.numeric(threshold) - thresholds$threshold[[threshold_index]]) < 1e-12,
                        sum(as.integer(variant_count), na.rm = TRUE)]
        if (length(hit)) variant_count <- hit[[1]]
      }
      rows[[row_index]] <- data.table(
        eid = target$eid, outcome_id = outcome_id,
        outcome_label = sources$outcome_label[[source_index]],
        threshold = thresholds$threshold[[threshold_index]],
        threshold_label = thresholds$label[[threshold_index]],
        score_column = score_name, prs_raw = raw,
        prs_z = (raw - mean(raw)) / score_sd,
        variant_count = variant_count
      )
      merge_qc[[row_index]] <- data.table(
        outcome_id = outcome_id, threshold = thresholds$threshold[[threshold_index]],
        n = length(raw), mean = mean(raw), sd = score_sd,
        min = min(raw), max = max(raw), variant_count = variant_count
      )
    }
  }
  scores <- rbindlist(rows)
  expected_rows <- nrow(target) * nrow(sources) * nrow(thresholds)
  if (nrow(scores) != expected_rows || anyDuplicated(scores[, .(eid, outcome_id, threshold)])) {
    stop("Merged PRS score contract failed.", call. = FALSE)
  }
  fwrite(scores, file.path(cfg$paths$prs, "prs_scores_long.csv.gz"), na = "")
  yur_write_csv(rbindlist(merge_qc), file.path(cfg$paths$prs, "prs_score_qc.csv"))
  yur_write_json(list(
    status = "PASS", rows = nrow(scores), participants = nrow(target),
    requested_target_participants = nrow(cohort),
    excluded_missing_genotype = nrow(excluded),
    outcomes = nrow(sources), thresholds = nrow(thresholds),
    requested_target_eid_hash = yur_sha_text(sort(cohort$eid)),
    target_eid_hash = yur_sha_text(sort(target$eid))
  ), file.path(cfg$paths$prs, "prs_merge_summary.json"))
  invisible(scores)
}

yur_prs_covariate_matrix <- function(meta, protein_panel = NULL) {
  panel_key <- yur_panel_key(protein_panel %||% "")
  requested_lag <- if (nzchar(panel_key)) paste0("protein_sampling_lag_days__", panel_key) else "protein_sampling_lag_days"
  lag_variable <- if (requested_lag %in% names(meta)) requested_lag else "protein_sampling_lag_days"
  x <- yur_impute_association_covariates(meta, lag_variable = lag_variable)
  for (v in paste0("PC", 1:10)) {
    x[, (v) := as.numeric(get(v))]
    x[!is.finite(get(v)), (v) := median(get(v), na.rm = TRUE)]
  }
  formula <- stats::as.formula(paste(
    "~ age + sex + ethnicity_white + tdi + blood_collection_season +",
    "protein_sampling_lag_days + fasting_time + sbp + bmi + smoking + alcohol +",
    paste(paste0("PC", 1:10), collapse = " + ")
  ))
  mm <- model.matrix(formula, data = x)
  attr(mm, "lag_variable") <- lag_variable
  mm
}

yur_prs_partial_regression <- function(protein, prs_matrix, covariates) {
  protein <- as.numeric(protein)
  ok <- is.finite(protein) & rowSums(!is.finite(prs_matrix)) == 0L & rowSums(!is.finite(covariates)) == 0L
  n <- sum(ok)
  if (n < 100L || sd(protein[ok]) <= 0) {
    return(data.table(threshold_index = seq_len(ncol(prs_matrix)), n = n, beta = NA_real_, se = NA_real_, p = NA_real_))
  }
  y <- as.numeric(scale(protein[ok]))
  cmat <- covariates[ok, , drop = FALSE]
  xprs <- prs_matrix[ok, , drop = FALSE]
  q <- qr(cmat, LAPACK = FALSE)
  residualized <- qr.resid(q, cbind(y, xprs))
  ry <- residualized[, 1]
  rank_c <- q$rank
  rows <- lapply(seq_len(ncol(xprs)), function(j) {
    rx <- residualized[, j + 1L]
    denominator <- sum(rx^2)
    if (!is.finite(denominator) || denominator <= 0) {
      return(data.table(threshold_index = j, n = n, beta = NA_real_, se = NA_real_, p = NA_real_))
    }
    beta <- sum(rx * ry) / denominator
    df <- n - rank_c - 1L
    sigma2 <- sum((ry - beta * rx)^2) / df
    se <- sqrt(sigma2 / denominator)
    p <- 2 * pnorm(abs(beta / se), lower.tail = FALSE)
    data.table(threshold_index = j, n = n, beta = beta, se = se, p = p)
  })
  rbindlist(rows)
}

yur_prs_load_meta_with_pcs <- function(cfg, target_eids = NULL) {
  meta <- yur_prs_target_cohort(cfg)
  if (!is.null(target_eids)) {
    target_eids <- yur_norm_eid(target_eids)
    if (anyDuplicated(target_eids)) stop("Requested PRS target EIDs are duplicated.", call. = FALSE)
    meta <- meta[eid %chin% target_eids]
    if (nrow(meta) != length(target_eids) || !setequal(meta$eid, target_eids)) {
      stop("PRS score EIDs do not match the Yu target cohort.", call. = FALSE)
    }
  }
  meta[, prs_order__ := .I]
  all <- as.data.table(readRDS(cfg$phenotype_rds))
  if (!"eid" %in% names(all)) stop("Phenotype RDS lacks eid.", call. = FALSE)
  pc_names <- paste0("PC", 1:10)
  missing <- setdiff(pc_names, names(all))
  if (length(missing)) stop("Phenotype RDS is missing PCs: ", paste(missing, collapse = ", "), call. = FALSE)
  pc <- all[, c("eid", pc_names), with = FALSE]
  pc <- yur_normalize_eid_column(pc)
  meta <- merge(meta, pc, by = "eid", all.x = TRUE, sort = FALSE)
  setorder(meta, prs_order__)
  meta[, prs_order__ := NULL]
  if (anyNA(meta[, ..pc_names])) stop("Some PRS target participants lack PC1-PC10.", call. = FALSE)
  meta
}

yur_prs_associate_shard <- function(cfg) {
  requested <- trimws(strsplit(as.character(cfg$endpoint_subset), ",", fixed = TRUE)[[1]])
  if (length(requested) != 1L || tolower(requested) == "all") {
    stop("PRS association shard requires exactly one --endpoint_subset.", call. = FALSE)
  }
  outcome_id <- requested[[1]]
  requested_outcome <- outcome_id
  sources <- yur_prs_sources(cfg)
  if (!outcome_id %chin% sources$outcome_id) stop("Unknown PRS outcome: ", outcome_id, call. = FALSE)
  thresholds <- yur_prs_thresholds(cfg)
  candidates <- yur_prs_candidate_table(cfg)[outcome_id == requested_outcome]
  out_file <- file.path(cfg$paths$prs, paste0("prs_protein_associations_", outcome_id, ".csv.gz"))
  if (!nrow(candidates)) {
    fwrite(data.table(), out_file)
    return(invisible(NULL))
  }
  score_file <- file.path(cfg$paths$prs, "prs_scores_long.csv.gz")
  if (!file.exists(score_file)) stop("Run PRS merge first.", call. = FALSE)
  scores <- fread(score_file)[outcome_id == requested_outcome]
  score_wide <- dcast(scores, eid ~ score_column, value.var = "prs_z")
  score_wide[, eid := yur_norm_eid(eid)]
  expected_score_columns <- thresholds$score_column
  if (!all(expected_score_columns %in% names(score_wide))) stop("PRS score matrix is incomplete for ", outcome_id, call. = FALSE)

  meta <- yur_prs_load_meta_with_pcs(cfg, score_wide$eid)
  meta <- merge(meta, score_wide, by = "eid", all = FALSE, sort = FALSE)
  raw_header <- names(fread(cfg$raw_protein_file, nrows = 0, showProgress = FALSE))
  eid_col <- intersect(c("eid", "id", "f.eid", "participant_id"), raw_header)[[1]]
  raw <- fread(cfg$raw_protein_file, select = c(eid_col, candidates$feature_id),
               showProgress = TRUE, nThread = cfg$workers)
  setnames(raw, eid_col, "eid")
  raw <- yur_normalize_eid_column(raw)
  raw <- raw[match(meta$eid, raw$eid)]
  if (!identical(raw$eid, meta$eid)) stop("Raw protein and PRS metadata EIDs are misaligned.", call. = FALSE)

  panel_keys <- unique(yur_panel_key(candidates$panel))
  covariates <- setNames(lapply(panel_keys, function(panel) yur_prs_covariate_matrix(meta, panel)), panel_keys)
  prs_matrix <- as.matrix(meta[, ..expected_score_columns])
  results <- vector("list", nrow(candidates))
  for (i in seq_len(nrow(candidates))) {
    panel_key <- yur_panel_key(candidates$panel[[i]])
    stat <- yur_prs_partial_regression(raw[[candidates$feature_id[[i]]]], prs_matrix, covariates[[panel_key]])
    stat[, `:=`(
      outcome_id = outcome_id,
      outcome_label = sources[outcome_id == requested_outcome, outcome_label][[1]],
      feature_id = candidates$feature_id[[i]], protein = candidates$protein[[i]],
      panel = candidates$panel[[i]], olink_id = candidates$olink_id[[i]],
      threshold = thresholds$threshold[threshold_index],
      threshold_label = thresholds$label[threshold_index],
      score_column = thresholds$score_column[threshold_index],
      lag_covariate = attr(covariates[[panel_key]], "lag_variable"),
      protein_scale = "per_analysis_SD", prs_scale = "per_analysis_SD",
      candidate_n = nrow(candidates), bonferroni_threshold = 0.05 / nrow(candidates)
    )]
    stat[, bonferroni_significant := is.finite(p) & p < bonferroni_threshold]
    results[[i]] <- stat
    if (i %% 100L == 0L) yur_log(cfg, "PRS association ", outcome_id, " ", i, "/", nrow(candidates))
  }
  result <- rbindlist(results, use.names = TRUE, fill = TRUE)
  setcolorder(result, c(
    "outcome_id", "outcome_label", "feature_id", "protein", "panel", "olink_id",
    "threshold", "threshold_label", "score_column", "n", "beta", "se", "p",
    "bonferroni_threshold", "bonferroni_significant", "candidate_n",
    "protein_scale", "prs_scale", "lag_covariate"
  ))
  expected_rows <- nrow(candidates) * nrow(thresholds)
  if (nrow(result) != expected_rows || anyDuplicated(result[, .(feature_id, threshold)])) {
    stop("PRS association shard is incomplete for ", outcome_id, call. = FALSE)
  }
  fwrite(result, out_file, na = "")
  yur_write_json(list(
    status = "PASS", outcome_id = outcome_id, candidates = nrow(candidates),
    thresholds = nrow(thresholds), rows = nrow(result), significant_rows = sum(result$bonferroni_significant),
    complete_non_significant_rows_retained = TRUE
  ), file.path(cfg$paths$prs, paste0("prs_association_", outcome_id, ".done.json")))
  invisible(result)
}

yur_prs_merge_associations <- function(cfg) {
  sources <- yur_prs_sources(cfg)
  thresholds <- yur_prs_thresholds(cfg)
  candidates <- yur_prs_candidate_table(cfg)
  rows <- lapply(sources$outcome_id, function(id_value) {
    path <- file.path(cfg$paths$prs, paste0("prs_protein_associations_", id_value, ".csv.gz"))
    if (!file.exists(path)) stop("Missing PRS association shard: ", path, call. = FALSE)
    x <- fread(path)
    expected <- candidates[outcome_id == id_value, uniqueN(feature_id)] * nrow(thresholds)
    if (nrow(x) != expected) stop("PRS association row-count mismatch for ", id_value, call. = FALSE)
    x
  })
  all <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  all[, pair_significant_any := any(bonferroni_significant), by = .(outcome_id, feature_id)]
  expected_rows <- nrow(candidates) * nrow(thresholds)
  if (nrow(all) != expected_rows || anyDuplicated(all[, .(outcome_id, feature_id, threshold)])) {
    stop("Complete PRS-protein association matrix failed its row contract.", call. = FALSE)
  }
  if (any(all[, .N, by = .(outcome_id, feature_id)]$N != nrow(thresholds))) {
    stop("At least one PRS-protein pair does not contain all five thresholds.", call. = FALSE)
  }
  fwrite(all, file.path(cfg$paths$prs, "prs_protein_associations.csv.gz"), na = "")
  yur_write_csv(all[bonferroni_significant == TRUE],
                file.path(cfg$paths$prs, "prs_protein_associations_significant.csv"))
  figure <- all[pair_significant_any == TRUE]
  figure[, `:=`(disease = outcome_label, Protein = protein, p_value = p)]
  yur_write_csv(figure, file.path(cfg$paths$prs, "figure6a_source_data.csv"))
  summary <- all[, .(
    candidate_proteins = uniqueN(feature_id), complete_rows = .N,
    significant_threshold_rows = sum(bonferroni_significant),
    proteins_significant_any_threshold = uniqueN(feature_id[pair_significant_any])
  ), by = .(outcome_id, outcome_label)]
  yur_write_csv(summary, file.path(cfg$paths$prs, "prs_association_summary.csv"))
  yur_write_json(list(
    status = "PASS", rows = nrow(all), expected_rows = expected_rows,
    all_non_significant_rows_retained = TRUE,
    significant_rows = sum(all$bonferroni_significant),
    significant_proteins = uniqueN(all[bonferroni_significant == TRUE, feature_id]),
    figure_rows = nrow(figure)
  ), file.path(cfg$paths$prs, "prs_association_merge_summary.json"))
  invisible(all)
}

yur_prs_report <- function(cfg) {
  required <- c("prs_preflight_summary.json", "prs_merge_summary.json", "prs_association_merge_summary.json")
  paths <- file.path(cfg$paths$prs, required)
  if (any(!file.exists(paths))) stop("PRS report prerequisites are incomplete.", call. = FALSE)
  preflight <- read_json(paths[[1]], simplifyVector = TRUE)
  score <- read_json(paths[[2]], simplifyVector = TRUE)
  association <- read_json(paths[[3]], simplifyVector = TRUE)
  summary <- fread(file.path(cfg$paths$prs, "prs_association_summary.csv"))
  lines <- c(
    "# Yu/Chen PRS-protein local reconstruction",
    "",
    "## Scope",
    "",
    "GWAS effect alleles and beta values from the frozen article-linked sources are used as reference weights.",
    "This remains a local reconstruction because final author PRS files and complete clumping parameters were not published.",
    "All five threshold-specific regression coefficients are retained, including non-significant results.",
    "",
    sprintf("- Target participants: %s", format(preflight$target_participants, big.mark = ",")),
    sprintf("- PRS outcomes: %s", score$outcomes),
    sprintf("- PRS thresholds: %s", score$thresholds),
    sprintf("- Complete PRS-protein rows: %s", format(association$rows, big.mark = ",")),
    sprintf("- Significant threshold rows: %s", format(association$significant_rows, big.mark = ",")),
    "",
    "## Outcome summary",
    "",
    paste(c("|Outcome|Candidates|Complete rows|Significant rows|Proteins significant at any threshold|",
            "|---|---:|---:|---:|---:|",
            vapply(seq_len(nrow(summary)), function(i) sprintf(
              "|%s|%d|%d|%d|%d|", summary$outcome_label[[i]], summary$candidate_proteins[[i]],
              summary$complete_rows[[i]], summary$significant_threshold_rows[[i]],
              summary$proteins_significant_any_threshold[[i]]
            ), character(1))), collapse = "\n")
  )
  writeLines(lines, file.path(cfg$paths$prs, "PRS_LOCAL_RECONSTRUCTION_REPORT.md"))
}
