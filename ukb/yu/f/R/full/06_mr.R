yur_mr_candidate_file <- function(cfg) file.path(cfg$paths$mr, "mr_candidates.csv")
yur_mr_instrument_file <- function(cfg) file.path(cfg$paths$mr, "mr_exposure_instruments.csv.gz")

yur_prepare_local_mr <- function(cfg) {
  cox_file <- file.path(cfg$paths$cox, "table_s2_incident_associations.csv.gz")
  if (!file.exists(cox_file)) stop("Local Cox results are required before MR: ", cox_file, call. = FALSE)
  pqtl_root <- cfg$pqtl_root %||% file.path(cfg$dir0, "ppp", "clean")
  if (!dir.exists(pqtl_root)) stop("Local pQTL root is missing: ", pqtl_root, call. = FALSE)

  cox <- fread(cox_file, showProgress = FALSE)
  needed <- c("outcome_id", "outcome_label", "feature_id", "protein", "p", "bonferroni_significant")
  miss <- setdiff(needed, names(cox))
  if (length(miss)) stop("Cox results miss MR fields: ", paste(miss, collapse = ", "), call. = FALSE)
  candidates <- unique(cox[bonferroni_significant %in% TRUE, .(
    outcome_id, outcome_label, feature_id, protein,
    observational_beta = beta, observational_p = p
  )])
  setorder(candidates, protein, outcome_id, observational_p, feature_id)
  candidates <- candidates[, .SD[1L], by = .(protein, outcome_id)]
  if (!nrow(candidates)) stop("No local Bonferroni-significant protein-outcome candidates.", call. = FALSE)

  proteins <- unique(candidates$protein)
  rows <- vector("list", length(proteins))
  audit <- vector("list", length(proteins))
  for (i in seq_along(proteins)) {
    protein <- proteins[[i]]
    f <- file.path(pqtl_root, protein, paste0(protein, ".clump.assoc"))
    if (!file.exists(f) || is.na(file.info(f)$size) || file.info(f)$size == 0) {
      audit[[i]] <- data.table(protein, status = "MISSING_OR_EMPTY_CLUMP", file = f, n_instruments = 0L)
      next
    }
    x <- tryCatch(fread(f, showProgress = FALSE), error = function(e) NULL)
    req <- c("CHR", "POS", "NEA", "EA", "EAF", "N", "BETA", "SE", "P", "SNP")
    if (is.null(x) || length(setdiff(req, names(x)))) {
      audit[[i]] <- data.table(protein, status = "INVALID_CLUMP_SCHEMA", file = f, n_instruments = 0L)
      next
    }
    x <- x[P < 5e-8 & is.finite(BETA) & is.finite(SE) & SE > 0]
    x[, `:=`(
      protein = protein, snp = sub(",.*$", "", as.character(SNP)),
      chr = as.character(CHR), pos = as.integer(POS),
      effect_allele = toupper(as.character(EA)), other_allele = toupper(as.character(NEA)),
      eaf = as.numeric(EAF), beta_exposure = as.numeric(BETA),
      se_exposure = as.numeric(SE), p_exposure = as.numeric(P), n_exposure = as.numeric(N),
      f_statistic = (as.numeric(BETA) / as.numeric(SE))^2,
      instrument_rule = "p<5e-8; pre-clumped local UKB-PPP resource"
    )]
    x <- x[nzchar(snp) & !duplicated(snp)]
    rows[[i]] <- x[, .(
      protein, snp, chr, pos, effect_allele, other_allele, eaf,
      beta_exposure, se_exposure, p_exposure, n_exposure, f_statistic, instrument_rule
    )]
    audit[[i]] <- data.table(
      protein, status = ifelse(nrow(x) > 2L, "PASS", "RELAXED_P_CLUMP_REQUIRED"),
      file = f, n_instruments = nrow(x)
    )
  }
  instruments <- rbindlist(rows, fill = TRUE)
  instrument_audit <- rbindlist(audit, fill = TRUE)
  if (!nrow(instruments)) stop("No usable local pQTL instruments were found.", call. = FALSE)

  candidates <- merge(candidates, instrument_audit[, .(protein, instrument_status = status, n_instruments)],
                      by = "protein", all.x = TRUE)
  fwrite(instruments, yur_mr_instrument_file(cfg), na = "")
  yur_write_csv(candidates, yur_mr_candidate_file(cfg))
  yur_write_csv(instrument_audit, file.path(cfg$paths$mr, "mr_instrument_qc.csv"))
  union <- unique(instruments[, .(snp, chr, pos)])
  yur_write_csv(union, file.path(cfg$paths$mr, "mr_union_variants.csv"))
  yur_write_json(list(
    status = "PASS", candidate_pairs = nrow(candidates), candidate_proteins = uniqueN(candidates$protein),
    proteins_with_instruments = uniqueN(instruments$protein), instruments = nrow(instruments),
    relaxed_clumping_required = instrument_audit[status == "RELAXED_P_CLUMP_REQUIRED", .N],
    pqtl_root = normalizePath(pqtl_root, winslash = "/", mustWork = TRUE),
    cox_sha256 = yur_sha_file(cox_file), instrument_sha256 = yur_sha_file(yur_mr_instrument_file(cfg))
  ), file.path(cfg$paths$mr, "mr_prepare_summary.json"))
}

yur_harmonise_mr_pair <- function(exp, out) {
  if (!nrow(exp) || !nrow(out)) return(data.table())
  out <- copy(out)
  out[, `:=`(snp = as.character(snp), chr = as.character(chr), pos = as.integer(pos))]
  z <- merge(exp, out, by = "snp", allow.cartesian = TRUE)
  if (!nrow(z)) {
    exp[, key := paste(chr, pos, sep = ":")]
    out[, key := paste(chr, pos, sep = ":")]
    z <- merge(exp, out, by = "key", allow.cartesian = TRUE, suffixes = c("", "_out"))
  }
  if (!nrow(z)) return(data.table())
  z[, `:=`(
    effect_allele = toupper(effect_allele), other_allele = toupper(other_allele),
    effect_allele_outcome = toupper(effect_allele_outcome),
    other_allele_outcome = toupper(other_allele_outcome)
  )]
  z[, alignment := fifelse(
    effect_allele == effect_allele_outcome & other_allele == other_allele_outcome, "direct",
    fifelse(effect_allele == other_allele_outcome & other_allele == effect_allele_outcome, "flip", "mismatch")
  )]
  z <- z[alignment != "mismatch"]
  if (!nrow(z)) return(data.table())
  z[alignment == "flip", `:=`(beta_outcome = -beta_outcome, eaf_outcome = 1 - eaf_outcome)]
  z[, palindromic := paste0(effect_allele, other_allele) %chin% c("AT", "TA", "CG", "GC")]
  z[, ambiguous_palindrome := palindromic & is.finite(eaf) & is.finite(eaf_outcome) &
       abs(eaf - .5) < .08 & abs(eaf_outcome - .5) < .08]
  z <- z[!ambiguous_palindrome & is.finite(beta_exposure) & is.finite(se_exposure) & se_exposure > 0 &
           is.finite(beta_outcome) & is.finite(se_outcome) & se_outcome > 0]
  setorder(z, p_exposure, p_outcome)
  z[!duplicated(snp)]
}

yur_run_mr_methods <- function(h, protein, outcome_id) {
  if (!nrow(h)) return(data.table())
  d <- data.frame(
    SNP = h$snp, beta.exposure = h$beta_exposure, beta.outcome = h$beta_outcome,
    se.exposure = h$se_exposure, se.outcome = h$se_outcome,
    effect_allele.exposure = h$effect_allele,
    other_allele.exposure = h$other_allele,
    effect_allele.outcome = h$effect_allele_outcome,
    other_allele.outcome = h$other_allele_outcome,
    eaf.exposure = h$eaf, eaf.outcome = h$eaf_outcome,
    pval.exposure = h$p_exposure, pval.outcome = h$p_outcome,
    exposure = protein, outcome = outcome_id,
    id.exposure = protein, id.outcome = outcome_id, mr_keep = TRUE,
    stringsAsFactors = FALSE
  )
  methods <- if (nrow(d) == 1L) "mr_wald_ratio" else c("mr_ivw", "mr_weighted_median")
  if (nrow(d) >= 3L) methods <- c(methods, "mr_egger_regression")
  fit <- tryCatch(suppressMessages(TwoSampleMR::mr(d, method_list = methods)), error = function(e) NULL)
  if (is.null(fit) || !nrow(fit)) return(data.table())
  as.data.table(fit)[, .(
    protein, outcome_id, method, nsnp = as.integer(nsnp), beta = as.numeric(b),
    se = as.numeric(se), p_value = as.numeric(pval), or = exp(as.numeric(b)),
    ci_low = exp(as.numeric(b) - 1.96 * as.numeric(se)),
    ci_high = exp(as.numeric(b) + 1.96 * as.numeric(se))
  )]
}

yur_run_local_mr <- function(cfg) {
  if (!requireNamespace("TwoSampleMR", quietly = TRUE)) stop("TwoSampleMR is required.", call. = FALSE)
  candidates_file <- yur_mr_candidate_file(cfg)
  instruments_file <- yur_mr_instrument_file(cfg)
  lookup_dir <- cfg$mr_outcome_lookup_dir %||% file.path(cfg$paths$mr, "outcome_lookup")
  if (!file.exists(candidates_file) || !file.exists(instruments_file)) {
    stop("Run mode=mr_prepare before mr_run.", call. = FALSE)
  }
  lookup_files <- list.files(lookup_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)
  if (!length(lookup_files)) stop("No extracted outcome GWAS lookups in: ", lookup_dir, call. = FALSE)
  candidates <- fread(candidates_file)
  instruments <- fread(instruments_file)
  lookup <- rbindlist(lapply(lookup_files, fread, showProgress = FALSE), fill = TRUE)
  required <- c("outcome_id", "snp", "chr", "pos", "effect_allele_outcome", "other_allele_outcome",
                "eaf_outcome", "beta_outcome", "se_outcome", "p_outcome")
  miss <- setdiff(required, names(lookup))
  if (length(miss)) stop("Outcome lookup misses: ", paste(miss, collapse = ", "), call. = FALSE)

  results <- list(); harmonised <- list(); qc <- list(); k <- 0L
  for (i in seq_len(nrow(candidates))) {
    cand <- candidates[i]
    e <- instruments[protein == cand$protein]
    o <- lookup[outcome_id == cand$outcome_id]
    h <- yur_harmonise_mr_pair(e, o)
    k <- k + 1L
    qc[[k]] <- data.table(
      protein = cand$protein, outcome_id = cand$outcome_id,
      n_exposure_instruments = nrow(e), n_outcome_matches = nrow(o), n_harmonised = nrow(h),
      status = ifelse(nrow(h), "PASS", "NO_HARMONISED_INSTRUMENT")
    )
    if (!nrow(h)) next
    h[, `:=`(protein = cand$protein, outcome_id = cand$outcome_id)]
    harmonised[[k]] <- h
    results[[k]] <- yur_run_mr_methods(h, cand$protein, cand$outcome_id)
  }
  res <- rbindlist(results, fill = TRUE)
  harm <- rbindlist(harmonised, fill = TRUE)
  qc_dt <- rbindlist(qc, fill = TRUE)
  if (!nrow(res)) stop("No MR model produced a result; inspect mr_pair_qc.csv.", call. = FALSE)
  res[, primary_method := grepl("Inverse variance weighted|Wald ratio", method, ignore.case = TRUE)]
  res[, fdr := NA_real_]
  res[primary_method %in% TRUE, fdr := p.adjust(p_value, method = "BH")]
  res[, fdr_within_method := p.adjust(p_value, method = "BH"), by = method]
  res[, significant_nominal := p_value < .05]
  res[, significant_fdr := fdr < .05]
  res <- merge(res, unique(candidates[, .(protein, outcome_id, outcome_label, feature_id,
                                          observational_beta, observational_p)]),
               by = c("protein", "outcome_id"), all.x = TRUE)
  setcolorder(res, c("outcome_id", "outcome_label", "protein", "feature_id", "method", "nsnp",
                     "beta", "se", "or", "ci_low", "ci_high", "p_value", "fdr", "fdr_within_method",
                     "significant_nominal", "significant_fdr", "primary_method",
                     "observational_beta", "observational_p"))
  yur_write_csv(res, file.path(cfg$paths$mr, "mr_results.csv"))
  fwrite(harm, file.path(cfg$paths$mr, "mr_instruments_harmonised.csv.gz"), na = "")
  yur_write_csv(qc_dt, file.path(cfg$paths$mr, "mr_pair_qc.csv"))
  yur_write_json(list(
    status = "PASS", candidate_pairs = nrow(candidates), tested_pairs = uniqueN(res[, paste(protein, outcome_id)]),
    primary_nominal_pairs = res[primary_method & significant_nominal, .N],
    primary_fdr_pairs = res[primary_method & significant_fdr, .N],
    harmonised_instruments = nrow(harm), lookup_files = length(lookup_files),
    method = "TwoSampleMR: IVW or Wald ratio; weighted median and MR-Egger sensitivity"
  ), file.path(cfg$paths$mr, "mr_summary.json"))
}
