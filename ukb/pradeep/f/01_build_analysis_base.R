# Build the local analysis base used by all reproduction scripts.

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

message("Step 01: building analysis base")

followup_cutoff_audit <- data.table::data.table(
  analysis_project = analysis_project,
  resolved_follow_end = format(follow_end),
  follow_end_source = follow_end_source,
  required_follow_end = ifelse(is.na(required_follow_end), NA_character_, format(required_follow_end)),
  gate_pass = is.na(required_follow_end) || identical(follow_end, required_follow_end)
)
repro_write_csv(followup_cutoff_audit, file.path(audit_dir, "followup_cutoff_audit.csv"))
if (!isTRUE(followup_cutoff_audit$gate_pass[[1]])) {
  stop("FOLLOW_END_GATE_FAILED before outcome construction.", call. = FALSE)
}

input_audit <- do.call(rbind, Map(repro_file_status, input_files, names(input_files)))
repro_write_csv(input_audit, file.path(audit_dir, "input_audit.csv"))
print(input_audit)

if (!file.exists(input_files$all_rds)) stop("Missing all.rds: ", input_files$all_rds, call. = FALSE)
if (!file.exists(input_files$prot_rds)) stop("Missing prot.rds: ", input_files$prot_rds, call. = FALSE)

message("Reading phenotype: ", input_files$all_rds)
phe <- readRDS(input_files$all_rds)
message("Reading proteins: ", input_files$prot_rds)
prot <- readRDS(input_files$prot_rds)

phe <- data.table::as.data.table(phe)
prot <- data.table::as.data.table(prot)
phe[, eid := as.character(eid)]
prot[, eid := as.character(eid)]

if (!is.na(analysis_population)) {
  pop_col <- repro_first_existing(phe, c("ethnic.c", "gen_ethnicity", "ethnicity"), required = FALSE)
  if (!is.na(pop_col)) {
    phe <- phe[as.character(get(pop_col)) == analysis_population]
  } else {
    warning("Population filter requested but no ethnicity column found; using all participants.")
  }
}

message("Merging phenotype and protein data by eid")
dt <- merge(phe, prot, by = "eid")
rm(phe, prot)
gc()

dt <- repro_make_birth_date(dt)
if (!"date_lost" %in% names(dt)) dt[, date_lost := as.Date(NA)]
if (!"date_death" %in% names(dt)) dt[, date_death := as.Date(NA)]

qc <- data.table::data.table(stage = "after_merge_population_filter", n = nrow(dt), n_proteins = ncol(dt) - 1L)

message("Constructing outcome time axes")
outcome_source_audit <- list()
for (i in seq_len(nrow(outcome_map))) {
  key <- outcome_map$outcome_key[i]
  source_for_key <- outcome_date_sources[[key]]
  if (is.null(source_for_key) || is.na(source_for_key) || !nzchar(source_for_key)) source_for_key <- date_source
  date_col <- repro_resolve_outcome_date_col(dt, key, source_for_key, outcome_map)
  message("  ", key, " <- ", date_col, " (", source_for_key, ")")
  axis <- repro_construct_outcome_axis(
    dt = dt,
    outcome_date_col = date_col,
    baseline_col = "date_attend",
    birth_date_col = "birth_date",
    lost_date_col = "date_lost",
    death_date_col = "date_death",
    admin_end = follow_end,
    domain = "cvd",
    same_day_policy = "exclude"
  )
  dt[[paste0(key, "_inc")]] <- axis$incident
  dt[[paste0(key, "_prev")]] <- axis$prevalent
  dt[[paste0(key, "_fu")]] <- axis$followup
  dt[[paste0(key, "_r2e")]] <- axis$r2e
  dt[[paste0(key, "_b2e")]] <- axis$b2e
  dt[[paste0(key, "_event_date")]] <- axis$event_date
  dt[[paste0(key, "_censor_date")]] <- axis$censor_date
  outcome_source_audit[[key]] <- cbind(data.table::data.table(
    outcome = key,
    label = outcome_map$label[i],
    requested_date_source = source_for_key,
    date_column = date_col,
    outcome_logic_version = "administrative_loss_death_censor_v2"
  ), axis$audit)
}
outcome_source_audit <- data.table::rbindlist(outcome_source_audit, fill = TRUE)
if (any(outcome_source_audit$n_incident_after_admin_end > 0L) ||
    any(outcome_source_audit$n_incident_after_person_censor > 0L) ||
    any(outcome_source_audit$n_followup_after_admin_end > 0L)) {
  stop("OUTCOME_CENSORING_GATE_FAILED before cohort construction.", call. = FALSE)
}

message("Mapping local covariates to paper-style names")
for (target in names(covariate_source_map)) {
  source_col <- covariate_source_map[[target]]
  if (source_col %in% names(dt)) {
    dt[[target]] <- dt[[source_col]]
  } else {
    dt[[target]] <- NA
  }
}
for (pc in paste0("PC", 1:10)) {
  if (!pc %in% names(dt)) dt[[pc]] <- NA_real_
}

dt[, id := eid]
dt[, Sex_numeric := suppressWarnings(as.integer(Sex_numeric))]
dt[, age := suppressWarnings(as.numeric(age))]
dt[, age2 := age^2]
dt[, mergedrace := repro_safe_factor(mergedrace)]
dt[, ever_smoked := repro_safe_factor(ifelse(suppressWarnings(as.numeric(ever_smoked)) > 0, 1, 0))]
dt[, antihtnbase := repro_safe_factor(ifelse(suppressWarnings(as.numeric(antihtnbase)) > 0, 1, 0))]
dt[, cholmed := repro_safe_factor(ifelse(suppressWarnings(as.numeric(cholmed)) > 0, 1, 0))]
dt[, dm2_prev := as.integer(suppressWarnings(as.numeric(dm2_prev)) > 0)]

dt[, tdi_log := as.numeric(scale(repro_rank_norm(tdi)))]

protein_cols_all <- setdiff(names(readRDS(input_files$prot_rds)), "eid")
protein_panel_source <- "all_local_proteins"
if (isTRUE(use_original_15k_panel) && file.exists(protein_map_15k_file)) {
  protein_map_15k <- data.table::fread(protein_map_15k_file, select = "Assay")
  original_panel <- unique(protein_map_15k$Assay)
  original_panel <- setdiff(original_panel, c("PCOLCE", "TACSTD2", "CTSS", "NPM1"))
  original_panel <- original_panel[original_panel %in% protein_cols_all]
  if (length(original_panel) > 0) {
    protein_cols <- original_panel
    protein_panel_source <- protein_map_15k_file
  } else {
    warning("Original 1.5k panel requested but no Assay symbols matched prot.rds; using all local proteins.")
    protein_cols <- protein_cols_all
  }
} else {
  protein_cols <- protein_cols_all
}
protein_cols <- intersect(protein_cols, names(dt))
if (max_proteins > 0L) protein_cols <- head(protein_cols, max_proteins)
if (length(protein_cols) == 0) stop("No protein columns found after merge.", call. = FALSE)
message("Protein panel source: ", protein_panel_source)
message("Proteins entering missingness filter: ", length(protein_cols))

message("Protein missingness filter")
for (p in protein_cols) dt[[p]] <- suppressWarnings(as.numeric(dt[[p]]))
protein_missing <- colMeans(is.na(as.data.frame(dt[, ..protein_cols])))
protein_keep <- names(protein_missing)[protein_missing <= max_missing_protein]
protein_drop <- setdiff(protein_cols, protein_keep)
protein_cols <- protein_keep
qc <- rbind(qc, data.table::data.table(stage = "after_protein_column_missingness_filter", n = nrow(dt), n_proteins = length(protein_cols)))

row_missing <- rowMeans(is.na(as.matrix(dt[, ..protein_cols])))
dt <- dt[row_missing <= max_missing_individual]
qc <- rbind(qc, data.table::data.table(stage = "after_participant_protein_missingness_filter", n = nrow(dt), n_proteins = length(protein_cols)))

protein_var <- vapply(protein_cols, function(p) stats::var(dt[[p]], na.rm = TRUE), numeric(1))
protein_cols <- names(protein_var)[is.finite(protein_var) & protein_var > 1e-8]
qc <- rbind(qc, data.table::data.table(stage = "after_zero_variance_protein_filter", n = nrow(dt), n_proteins = length(protein_cols)))

essential_covars <- c("Sex_numeric", "age", paste0("PC", 1:10))
essential_covars <- essential_covars[essential_covars %in% names(dt)]
dt <- dt[stats::complete.cases(dt[, ..essential_covars])]
qc <- rbind(qc, data.table::data.table(stage = "after_essential_covariate_filter", n = nrow(dt), n_proteins = length(protein_cols)))

message("Excluding related participants")
relatedness <- repro_apply_relatedness_exclusion(
  dt,
  kinship_file = kinship_file,
  cutoff = relatedness_cutoff,
  seed = 1L,
  skip = skip_relatedness
)
dt <- relatedness$dt
qc <- rbind(qc, data.table::data.table(stage = "after_relatedness_exclusion", n = nrow(dt), n_proteins = length(protein_cols)))

prev_cols <- paste0(outcome_map$outcome_key, "_prev")
inc_cols <- paste0(outcome_map$outcome_key, "_inc")
dt <- dt[stats::complete.cases(dt[, ..prev_cols])]
dt <- dt[rowSums(as.data.frame(dt[, ..prev_cols]) == 1, na.rm = TRUE) == 0]
dt <- dt[stats::complete.cases(dt[, ..inc_cols])]
qc <- rbind(qc, data.table::data.table(stage = "after_prevalent_cardiac_disease_exclusion", n = nrow(dt), n_proteins = length(protein_cols)))

message("Imputing continuous clinical covariates")
predictor_covars <- c("Sex_numeric", "age", "mergedrace", paste0("PC", 1:10))
dt <- repro_impute_continuous(as.data.frame(dt), continuous_impute_vars, predictor_covars)
dt <- data.table::as.data.table(dt)

for (v in c("ever_smoked", "antihtnbase", "cholmed")) {
  if (v %in% names(dt) && any(is.na(dt[[v]]))) {
    mode_val <- repro_mode_value(dt[[v]])
    dt[[v]][is.na(dt[[v]])] <- mode_val
    dt[[v]] <- repro_safe_factor(dt[[v]])
  }
}
if ("dm2_prev" %in% names(dt)) dt[is.na(dm2_prev), dm2_prev := 0L]

message("Adding plate/well variables")
plate_well <- repro_add_plate_well(dt, plate_well_file = plate_well_file)
dt <- plate_well$dt

message("Imputing and scaling proteins")
mat <- as.matrix(dt[, ..protein_cols])
if (!skip_knn && length(protein_cols) > 1L) {
  imp <- impute::impute.knn(mat, rowmax = max_missing_individual, colmax = max_missing_protein, k = 10, rng.seed = 1)
  mat <- imp$data
} else {
  message("KNN imputation skipped; using column medians for missing protein values.")
  mat <- repro_median_impute_matrix(mat)
}
mat <- repro_scale_matrix(mat)
protein_cols <- colnames(mat)
dt[, (protein_cols) := as.data.table(mat)]
rm(mat)
gc()

endpoint_cols <- as.vector(outer(
  outcome_map$outcome_key,
  c("_inc", "_prev", "_fu", "_r2e", "_b2e", "_event_date", "_censor_date"),
  paste0
))
clinical_cols <- unique(c(
  "eid", "id", "date_attend", endpoint_cols, "Sex_numeric", "age", "age2", "center", "mergedrace", "plate", "well",
  "ever_smoked", "antihtnbase", "cholmed", "dm2_prev",
  paste0("PC", 1:10),
  continuous_impute_vars,
  paste0(continuous_impute_vars, "_missing"),
  paste0(continuous_impute_vars, "_final")
))
clinical_cols <- intersect(clinical_cols, names(dt))
dt_out <- dt[, c(clinical_cols, protein_cols), with = FALSE]

outcome_qc <- data.table::rbindlist(lapply(outcome_map$outcome_key, function(key) {
  event_date <- as.Date(dt_out[[paste0(key, "_event_date")]])
  censor_date <- as.Date(dt_out[[paste0(key, "_censor_date")]])
  inc <- dt_out[[paste0(key, "_inc")]]
  fu <- dt_out[[paste0(key, "_fu")]]
  data.table::data.table(
    outcome = key,
    outcome_logic_version = "administrative_loss_death_censor_v2",
    incident = sum(inc == 1, na.rm = TRUE),
    censored = sum(inc == 0, na.rm = TRUE),
    median_followup_year = stats::median(fu, na.rm = TRUE),
    max_followup_year = max(fu, na.rm = TRUE),
    incident_after_admin_end = sum(inc == 1 & event_date > follow_end, na.rm = TRUE),
    incident_after_person_censor = sum(inc == 1 & event_date > censor_date, na.rm = TRUE),
    followup_after_admin_end = sum(
      !is.na(fu) & fu > as.numeric(follow_end - as.Date(dt_out$date_attend)) / 365.25 + 1e-8,
      na.rm = TRUE
    )
  )
}))
if (any(outcome_qc$incident_after_admin_end > 0L) ||
    any(outcome_qc$incident_after_person_censor > 0L) ||
    any(outcome_qc$followup_after_admin_end > 0L)) {
  stop("OUTCOME_CENSORING_GATE_FAILED after cohort construction.", call. = FALSE)
}

missingness <- data.table::data.table(
  protein = names(protein_missing),
  missing_rate_before_filter = as.numeric(protein_missing),
  kept = names(protein_missing) %in% protein_cols
)

analysis_base <- list(
  dat = dt_out,
  protein_cols = protein_cols,
  outcome_map = outcome_map,
  clinical_covariates_base = clinical_covariates_base,
  clinical_predictors_for_lasso = clinical_predictors_for_lasso,
  qc = qc,
  outcome_qc = outcome_qc,
  protein_missingness = missingness,
  settings = list(
    analysis_project = analysis_project,
    analysis_population = analysis_population,
    date_source = date_source,
    outcome_date_sources = outcome_date_sources,
    follow_end = follow_end,
    follow_end_source = follow_end_source,
    required_follow_end = required_follow_end,
    outcome_logic_version = "administrative_loss_death_censor_v2",
    max_missing_protein = max_missing_protein,
    max_missing_individual = max_missing_individual,
    skip_knn = skip_knn,
    skip_relatedness = skip_relatedness,
    relatedness_cutoff = relatedness_cutoff,
    kinship_file = kinship_file,
    plate_well_file = plate_well_file,
    max_proteins = max_proteins,
    use_original_15k_panel = use_original_15k_panel,
    protein_panel_source = protein_panel_source
  )
)

saveRDS(analysis_base, analysis_base_file, compress = "gzip")
repro_write_csv(qc, file.path(audit_dir, "build_base_qc.csv"))
repro_write_csv(relatedness$audit, file.path(audit_dir, "relatedness_exclusion_audit.csv"))
repro_write_csv(data.frame(eid = relatedness$removed_ids), file.path(audit_dir, "relatedness_removed_ids.csv"))
if (nrow(relatedness$pairs) > 0) repro_write_csv(relatedness$pairs, file.path(audit_dir, "relatedness_pairs_used.csv"))
repro_write_csv(plate_well$audit, file.path(audit_dir, "plate_well_audit.csv"))
repro_write_csv(outcome_qc, file.path(audit_dir, "outcome_qc.csv"))
repro_write_csv(outcome_source_audit, file.path(audit_dir, "outcome_date_source_audit.csv"))
repro_write_csv(missingness, file.path(audit_dir, "protein_missingness.csv"))
repro_write_csv(data.frame(protein = protein_cols), file.path(audit_dir, "protein_columns_used.csv"))

message("Analysis base saved: ", analysis_base_file)
print(qc)
print(outcome_qc)
