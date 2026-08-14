# Sex interaction and sex-stratified Cox analyses.

locate_repro_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")))
  frames <- sys.frames()
  ofiles <- vapply(frames, function(frame) if (!is.null(frame$ofile)) frame$ofile else NA_character_, character(1))
  ofiles <- ofiles[!is.na(ofiles)]
  if (length(ofiles) > 0) return(dirname(normalizePath(ofiles[length(ofiles)], winslash = "/")))
  project_dir <- Sys.getenv("PRADEEP_PROJECT_DIR", unset = "")
  candidates <- if (nzchar(project_dir)) file.path(project_dir, "f") else getwd()
  hit <- candidates[file.exists(file.path(candidates, "00_config.R"))]
  if (length(hit) > 0) return(normalizePath(hit[1], winslash = "/"))
  stop("Cannot locate pradeep/f. Run this step through pradeep.sh.", call. = FALSE)
}
this_dir <- locate_repro_script_dir()
source(file.path(this_dir, "00_config.R"))
source(file.path(script_dir, "R", "repro_utils.R"))

message("Step 03: sex interaction and sex-stratified Cox")
if (!file.exists(analysis_base_file)) {
  stop("Missing analysis base. Run 01_build_analysis_base.R first: ", analysis_base_file, call. = FALSE)
}

base <- readRDS(analysis_base_file)
dt <- base$dat
proteins <- base$protein_cols
if (max_proteins > 0L) proteins <- head(proteins, max_proteins)
outcome_family <- base$outcome_map
outcomes <- outcome_family[match(analysis_outcome_keys, outcome_family$outcome_key)]
sex_chunk_size <- max(10L, repro_int_env("UKBPPP_SEX_CHUNK_SIZE", 50L))
sex_chunk_dir <- file.path(output_dir, "sex_interaction_chunks")
chunk_reuse <- repro_bool_env("UKBPPP_REUSE_CHUNKS", default = TRUE)
chunk_state <- repro_prepare_chunk_dir(
  sex_chunk_dir,
  list(
    script = "03_sex_stratified_interaction.R",
    analysis_base = repro_file_signature(analysis_base_file),
    n_proteins = length(proteins),
    proteins = paste(proteins, collapse = "|"),
    outcome_keys = paste(outcomes$outcome_key, collapse = "|"),
    outcome_labels = paste(outcomes$label, collapse = "|"),
    chunk_size = sex_chunk_size,
    max_proteins = max_proteins,
    min_events_per_model = min_events_per_model
  ),
  allow_reuse = chunk_reuse
)
message("  chunk cache: ", chunk_state, " (", sex_chunk_dir, ")")

empty_interaction <- function(protein, outcome, error, n = NA_integer_, events = NA_integer_) {
  data.frame(
    Protein = protein, Outcome = outcome,
    HR_Int = NA_real_, CI_Low_Int = NA_real_, CI_High_Int = NA_real_,
    SE_Int = NA_real_, P_Val_Int = NA_real_,
    HR_M = NA_real_, CI_Low_M = NA_real_, CI_High_M = NA_real_, P_Val_M = NA_real_,
    HR_F = NA_real_, CI_Low_F = NA_real_, CI_High_F = NA_real_, P_Val_F = NA_real_,
    n = n, events = events, error = error,
    stringsAsFactors = FALSE
  )
}

fit_interaction_one <- function(tv, protein, outcome_label, covars) {
  cols <- unique(c("tstart", "tstop", "event", protein, "Sex_numeric", covars))
  d <- repro_clean_model_data(tv, cols)
  n_events <- sum(d$event == 1, na.rm = TRUE)
  if (nrow(d) < 50 || n_events < min_events_per_model) {
    return(empty_interaction(protein, outcome_label, "too few rows/events", nrow(d), n_events))
  }
  if (length(unique(d$Sex_numeric)) < 2) {
    return(empty_interaction(protein, outcome_label, "only one sex level", nrow(d), n_events))
  }
  p <- suppressWarnings(as.numeric(d[[protein]]))
  if (stats::sd(p, na.rm = TRUE) <= 0 || all(is.na(p))) {
    return(empty_interaction(protein, outcome_label, "zero protein variance", nrow(d), n_events))
  }
  d$protein_z <- as.numeric(scale(p))
  covars <- setdiff(repro_usable_covariates(d, covars), "Sex_numeric")
  form <- stats::reformulate(c("protein_z * Sex_numeric", covars), response = "survival::Surv(tstart, tstop, event)")
  fit <- tryCatch(survival::coxph(form, data = d, ties = "efron"), error = function(e) e)
  if (inherits(fit, "error")) {
    return(empty_interaction(protein, outcome_label, conditionMessage(fit), nrow(d), n_events))
  }
  sm <- summary(fit)
  term <- grep("^protein_z:Sex_numeric|^Sex_numeric:protein_z", rownames(sm$coefficients), value = TRUE)
  if (length(term) == 0) {
    return(empty_interaction(protein, outcome_label, "interaction term dropped", nrow(d), n_events))
  }
  beta <- unname(sm$coefficients[term[1], "coef"])
  se <- unname(sm$coefficients[term[1], "se(coef)"])
  pval <- unname(sm$coefficients[term[1], "Pr(>|z|)"])

  male <- repro_fit_one_cox(d[d$Sex_numeric == 1, , drop = FALSE], protein, outcome_label, setdiff(covars, "Sex_numeric"), min_events_per_model)
  female <- repro_fit_one_cox(d[d$Sex_numeric == 0, , drop = FALSE], protein, outcome_label, setdiff(covars, "Sex_numeric"), min_events_per_model)

  data.frame(
    Protein = protein,
    Outcome = outcome_label,
    HR_Int = exp(beta),
    CI_Low_Int = exp(beta - 1.96 * se),
    CI_High_Int = exp(beta + 1.96 * se),
    SE_Int = se,
    P_Val_Int = pval,
    HR_M = male$HR[1],
    CI_Low_M = male$CI_Lower[1],
    CI_High_M = male$CI_Upper[1],
    P_Val_M = male$P_Value[1],
    HR_F = female$HR[1],
    CI_Low_F = female$CI_Lower[1],
    CI_High_F = female$CI_Upper[1],
    P_Val_F = female$P_Value[1],
    n = nrow(d),
    events = n_events,
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

all_results <- list()
all_qc <- list()
sex_worker_cap <- max(1L, repro_int_env("UKBPPP_SEX_WORKER_CAP", 4L))
sex_workers <- min(n_workers, sex_worker_cap)
if (sex_workers < n_workers) {
  message("  sex-interaction worker cap: requested=", n_workers, "; effective=", sex_workers,
          " to avoid duplicating the full time-varying table across too many workers")
}

for (i in seq_len(nrow(outcomes))) {
  key <- outcomes$outcome_key[i]
  label <- outcomes$label[i]
  message("Outcome: ", label, " (", key, ")")
  tv <- repro_make_timevarying_data(dt, key, outcome_family)
  tv_covars <- setdiff(outcome_family$outcome_key, key)
  covars <- unique(c(base$clinical_covariates_base, tv_covars))
  covars <- repro_usable_covariates(tv, covars)

  chunks <- split(proteins, ceiling(seq_along(proteins) / sex_chunk_size))
  cl <- NULL
  if (sex_workers > 1L) {
    message("  using ", sex_workers, " parallel workers for sex-interaction chunks")
    cl <- parallel::makeCluster(sex_workers)
    parallel::clusterEvalQ(cl, {
      library(survival)
      NULL
    })
    parallel::clusterExport(
      cl,
      varlist = c(
        "tv", "label", "covars", "min_events_per_model",
        "fit_interaction_one", "empty_interaction",
        "repro_fit_one_cox", "repro_clean_model_data", "repro_empty_cox_result",
        "repro_usable_covariates", "repro_has_variation"
      ),
      envir = environment()
    )
  }
  chunk_files <- character(length(chunks))
  for (chunk_idx in seq_along(chunks)) {
    chunk_file <- file.path(sex_chunk_dir, sprintf("%s_chunk_%03d.csv", key, chunk_idx))
    chunk_files[[chunk_idx]] <- chunk_file
    if (repro_chunk_file_matches(chunk_file, chunks[[chunk_idx]], label)) {
      message("  chunk ", chunk_idx, "/", length(chunks), " already exists; reusing")
      next
    }
    message("  chunk ", chunk_idx, "/", length(chunks), " (", length(chunks[[chunk_idx]]), " proteins)")
    if (!is.null(cl)) {
      chunk_res <- tryCatch(
        do.call(rbind, parallel::parLapply(cl, chunks[[chunk_idx]], function(p) {
          fit_interaction_one(tv, p, label, covars)
        })),
        error = function(e) {
          message("  parallel worker failure: ", conditionMessage(e), "; retrying this chunk serially")
          NULL
        }
      )
      if (is.null(chunk_res)) {
        try(parallel::stopCluster(cl), silent = TRUE)
        cl <- NULL
        chunk_res <- do.call(rbind, lapply(chunks[[chunk_idx]], function(p) {
          fit_interaction_one(tv, p, label, covars)
        }))
      }
    } else {
      chunk_res <- do.call(rbind, lapply(chunks[[chunk_idx]], function(p) {
        fit_interaction_one(tv, p, label, covars)
      }))
    }
    repro_write_csv(chunk_res, chunk_file)
    rm(chunk_res)
    gc()
  }
  if (!is.null(cl)) parallel::stopCluster(cl)
  res <- do.call(rbind, lapply(chunk_files, data.table::fread))
  res$FDR_Int <- p.adjust(res$P_Val_Int, method = "BH")
  all_results[[key]] <- res
  all_qc[[key]] <- data.frame(
    outcome = key,
    label = label,
    n_rows_timevarying = nrow(tv),
    n_participants = length(unique(tv$id)),
    n_events = sum(tv$event == 1, na.rm = TRUE),
    n_proteins = length(proteins),
    stringsAsFactors = FALSE
  )
  rm(tv)
  gc()
}

sex_res <- do.call(rbind, all_results)
sex_res <- sex_res[order(sex_res$Outcome, sex_res$P_Val_Int), ]
qc <- do.call(rbind, all_qc)

# Match the paper's two-stage sex-interaction multiplicity rule. Candidate
# protein-outcome pairs first need a sex-stratified association passing the
# Bonferroni threshold across all tested protein-outcome pairs. Interaction
# P values are then Bonferroni-corrected within that candidate set.
n_sex_tests <- sum(is.finite(sex_res$P_Val_M) | is.finite(sex_res$P_Val_F))
sex_association_threshold <- if (n_sex_tests > 0L) 0.05 / n_sex_tests else NA_real_
sex_res$Paper_Candidate <- is.finite(sex_res$P_Val_M) & sex_res$P_Val_M < sex_association_threshold |
  is.finite(sex_res$P_Val_F) & sex_res$P_Val_F < sex_association_threshold
n_paper_candidates <- sum(sex_res$Paper_Candidate, na.rm = TRUE)
paper_interaction_threshold <- if (n_paper_candidates > 0L) 0.05 / n_paper_candidates else NA_real_
sex_res$Paper_Interaction_Threshold <- paper_interaction_threshold
sex_res$Paper_Interaction_P_Adjusted <- NA_real_
if (n_paper_candidates > 0L) {
  candidate_idx <- which(sex_res$Paper_Candidate & is.finite(sex_res$P_Val_Int))
  sex_res$Paper_Interaction_P_Adjusted[candidate_idx] <- p.adjust(
    sex_res$P_Val_Int[candidate_idx],
    method = "bonferroni",
    n = n_paper_candidates
  )
}
sex_res$Paper_Interaction_Bonferroni <- sex_res$Paper_Candidate &
  is.finite(sex_res$P_Val_Int) & sex_res$P_Val_Int < paper_interaction_threshold

paper_summary <- do.call(rbind, lapply(c(unique(as.character(sex_res$Outcome)), "ALL"), function(outcome) {
  d <- if (identical(outcome, "ALL")) sex_res else sex_res[as.character(sex_res$Outcome) == outcome, , drop = FALSE]
  data.frame(
    outcome = outcome,
    n_tests = nrow(d),
    n_candidates = sum(d$Paper_Candidate, na.rm = TRUE),
    n_significant_interactions = sum(d$Paper_Interaction_Bonferroni, na.rm = TRUE),
    sex_association_threshold = sex_association_threshold,
    paper_interaction_threshold = paper_interaction_threshold,
    stringsAsFactors = FALSE
  )
}))
paper_significant <- sex_res[sex_res$Paper_Interaction_Bonferroni %in% TRUE, , drop = FALSE]
paper_significant <- paper_significant[order(paper_significant$P_Val_Int), , drop = FALSE]

repro_write_csv(sex_res, sex_interaction_file)
repro_write_csv(qc, file.path(audit_dir, "sex_interaction_qc.csv"))
repro_write_csv(paper_summary, file.path(audit_dir, "sex_interaction_paper_rule_summary.csv"))
repro_write_csv(paper_significant, file.path(output_dir, "sex_interaction_paper_significant.csv"))

message("Sex interaction results saved: ", sex_interaction_file)
message(
  "Paper-rule candidates/significant interactions: ",
  n_paper_candidates, "/", nrow(paper_significant),
  " (interaction threshold ", format(paper_interaction_threshold, scientific = TRUE), ")"
)
print(qc)
print(paper_summary)
