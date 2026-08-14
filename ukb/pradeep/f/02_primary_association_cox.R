# Primary proteomic association analysis.
# Reproduces the paper's adjusted Cox workflow for CAD, AF, HF and AS.

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

message("Step 02: primary adjusted Cox association")
if (!file.exists(analysis_base_file)) {
  stop("Missing analysis base. Run 01_build_analysis_base.R first: ", analysis_base_file, call. = FALSE)
}

base <- readRDS(analysis_base_file)
dt <- base$dat
proteins <- base$protein_cols
if (max_proteins > 0L) proteins <- head(proteins, max_proteins)
outcome_family <- base$outcome_map
outcomes <- outcome_family[match(analysis_outcome_keys, outcome_family$outcome_key)]
bonferroni_threshold <- 0.05 / (length(proteins) * nrow(outcome_family))
message("  Bonferroni threshold: ", signif(bonferroni_threshold, 4),
        " (0.05 / ", length(proteins), " proteins / ", nrow(outcome_family),
        " paper-family outcomes; fitted=", paste(outcomes$outcome_key, collapse = ","), ")")
cox_chunk_size <- max(10L, repro_int_env("UKBPPP_COX_CHUNK_SIZE", 100L))
cox_chunk_dir <- file.path(output_dir, "primary_cox_chunks")
chunk_reuse <- repro_bool_env("UKBPPP_REUSE_CHUNKS", default = TRUE)
chunk_state <- repro_prepare_chunk_dir(
  cox_chunk_dir,
  list(
    script = "02_primary_association_cox.R",
    analysis_base = repro_file_signature(analysis_base_file),
    n_proteins = length(proteins),
    proteins = paste(proteins, collapse = "|"),
    outcome_keys = paste(outcomes$outcome_key, collapse = "|"),
    outcome_labels = paste(outcomes$label, collapse = "|"),
    chunk_size = cox_chunk_size,
    max_proteins = max_proteins,
    min_events_per_model = min_events_per_model
  ),
  allow_reuse = chunk_reuse
)
message("  chunk cache: ", chunk_state, " (", cox_chunk_dir, ")")
cox_worker_cap <- max(1L, repro_int_env("UKBPPP_COX_WORKER_CAP", 4L))
cox_workers <- min(n_workers, cox_worker_cap)
if (cox_workers < n_workers) {
  message("  Cox worker cap: requested=", n_workers, "; effective=", cox_workers,
          " to avoid duplicating the full time-varying table across too many workers")
}

fit_protein_chunk <- function(protein_chunk, tv, label, covars) {
  do.call(rbind, lapply(protein_chunk, function(p) {
    repro_fit_one_cox(tv, p, label, covars, min_events = min_events_per_model)
  }))
}

all_results <- list()
all_qc <- list()

for (i in seq_len(nrow(outcomes))) {
  key <- outcomes$outcome_key[i]
  label <- outcomes$label[i]
  message("Outcome: ", label, " (", key, ")")
  tv <- repro_make_timevarying_data(dt, key, outcome_family)
  tv_covars <- setdiff(outcome_family$outcome_key, key)
  covars <- unique(c(base$clinical_covariates_base, tv_covars))
  covars <- repro_usable_covariates(tv, covars)
  message("  covariates: ", paste(covars, collapse = ", "))

  chunks <- split(proteins, ceiling(seq_along(proteins) / cox_chunk_size))
  cl <- NULL
  if (cox_workers > 1L) {
    message("  using ", cox_workers, " parallel workers for Cox chunks")
    cl <- parallel::makeCluster(cox_workers)
    parallel::clusterEvalQ(cl, {
      library(survival)
      NULL
    })
    parallel::clusterExport(
      cl,
      varlist = c(
        "tv", "label", "covars", "min_events_per_model",
        "repro_fit_one_cox", "repro_clean_model_data", "repro_empty_cox_result",
        "repro_usable_covariates", "repro_has_variation"
      ),
      envir = environment()
    )
  }
  chunk_files <- character(length(chunks))
  for (chunk_idx in seq_along(chunks)) {
    chunk_file <- file.path(cox_chunk_dir, sprintf("%s_chunk_%03d.csv", key, chunk_idx))
    chunk_files[[chunk_idx]] <- chunk_file
    if (repro_chunk_file_matches(chunk_file, chunks[[chunk_idx]], label)) {
      message("  chunk ", chunk_idx, "/", length(chunks), " already exists; reusing")
      next
    }
    message("  chunk ", chunk_idx, "/", length(chunks), " (", length(chunks[[chunk_idx]]), " proteins)")
    if (!is.null(cl)) {
      chunk_res <- tryCatch(
        do.call(rbind, parallel::parLapply(cl, chunks[[chunk_idx]], function(p) {
          repro_fit_one_cox(tv, p, label, covars, min_events = min_events_per_model)
        })),
        error = function(e) {
          message("  parallel worker failure: ", conditionMessage(e), "; retrying this chunk serially")
          NULL
        }
      )
      if (is.null(chunk_res)) {
        try(parallel::stopCluster(cl), silent = TRUE)
        cl <- NULL
        chunk_res <- fit_protein_chunk(chunks[[chunk_idx]], tv, label, covars)
      }
    } else {
      chunk_res <- fit_protein_chunk(chunks[[chunk_idx]], tv, label, covars)
    }
    repro_write_csv(chunk_res, chunk_file)
    rm(chunk_res)
    gc()
  }
  if (!is.null(cl)) parallel::stopCluster(cl)
  res <- do.call(rbind, lapply(chunk_files, data.table::fread))
  res$FDR <- p.adjust(res$P_Value, method = "BH")
  res$Bonferroni_Threshold <- bonferroni_threshold
  res$Bonferroni <- res$P_Value < bonferroni_threshold
  all_results[[key]] <- res
  all_qc[[key]] <- data.frame(
    outcome = key,
    label = label,
    n_rows_timevarying = nrow(tv),
    n_participants = length(unique(tv$id)),
    n_events = sum(tv$event == 1, na.rm = TRUE),
    n_proteins = length(proteins),
    bonferroni_threshold = bonferroni_threshold,
    n_covariates = length(covars),
    covariates = paste(covars, collapse = ";"),
    stringsAsFactors = FALSE
  )
  rm(tv)
  gc()
}

assoc <- do.call(rbind, all_results)
assoc <- assoc[order(assoc$Outcome, assoc$P_Value), ]
qc <- do.call(rbind, all_qc)

repro_write_csv(assoc, primary_assoc_file)
repro_write_csv(qc, file.path(audit_dir, "primary_association_qc.csv"))

message("Primary association results saved: ", primary_assoc_file)
print(qc)
