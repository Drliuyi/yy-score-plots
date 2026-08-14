# Configuration for the command-driven Pradeep/Schuermans reproduction.

options(stringsAsFactors = FALSE)

repro_get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")))
  }
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(dirname(normalizePath(ofile, winslash = "/")))
  project_dir <- Sys.getenv("PRADEEP_PROJECT_DIR", unset = "")
  if (nzchar(project_dir)) return(file.path(project_dir, "f"))
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

repro_env_path <- function(name, default = "", required = FALSE) {
  value <- trimws(Sys.getenv(name, unset = default))
  value <- gsub("\\\\", "/", value)
  if (required && !nzchar(value)) {
    stop("Missing path setting ", name, ". Run pradeep.sh --h for examples.", call. = FALSE)
  }
  value
}

repro_bool_env <- function(name, default = FALSE) {
  val <- tolower(Sys.getenv(name, unset = if (default) "1" else "0"))
  val %in% c("1", "true", "t", "yes", "y")
}

repro_int_env <- function(name, default = 0L) {
  val <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  ifelse(is.na(val), default, val)
}

repro_load_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         ". Install them before running this workflow.", call. = FALSE)
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

repro_load_packages(c(
  "data.table", "dplyr", "survival", "impute", "glmnet", "pROC", "stats"
))

script_dir <- repro_get_script_dir()
project_dir <- repro_env_path("PRADEEP_PROJECT_DIR", dirname(script_dir), required = TRUE)
root_dir <- repro_env_path("PRADEEP_DIR0", if (.Platform$OS.type == "windows") "D:/" else "/mnt/d", required = TRUE)
phe_dir <- repro_env_path("PRADEEP_PHEDIR", file.path(root_dir, "data", "ukb", "phe"), required = TRUE)
script_root <- repro_env_path("PRADEEP_SCRIPT_ROOT", file.path(root_dir, "scripts"), required = TRUE)
helper_dir <- repro_env_path("PRADEEP_HELPER_DIR", file.path(script_root, "0f"), required = TRUE)
output_root <- repro_env_path("PRADEEP_OUTPUT_ROOT", file.path(root_dir, "analysis", "ukb"), required = TRUE)
analysis_project <- trimws(Sys.getenv(
  "UKBPPP_ANALYSIS_PROJECT",
  unset = "pradeep"
))
if (!nzchar(analysis_project) || grepl("[/\\\\]", analysis_project)) {
  stop("UKBPPP_ANALYSIS_PROJECT must be a single directory name.", call. = FALSE)
}
analysis_dir_override <- repro_env_path("PRADEEP_ANALYSIS_DIR", "")
analysis_dir <- if (nzchar(analysis_dir_override)) analysis_dir_override else file.path(output_root, analysis_project)
output_dir <- file.path(analysis_dir, "outputs")
audit_dir <- file.path(analysis_dir, "audit")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- list(
  all_rds = repro_env_path("PRADEEP_ALL_RDS", file.path(phe_dir, "Rdata", "all.rds")),
  prot_rds = repro_env_path("PRADEEP_PROT_RDS", file.path(phe_dir, "Rdata", "prot.rds")),
  pheno_tsv_gz = repro_env_path("PRADEEP_PHENO_TSV", file.path(phe_dir, "pheno.tsv.gz")),
  prot_tab_gz = repro_env_path("PRADEEP_RAW_PROTEIN", file.path(root_dir, "data.BIG", "gwas", "ppp", "prot.tab.gz"))
)

helper_file <- file.path(helper_dir, "phe.f.R")
if (file.exists(helper_file)) {
  source(helper_file)
}

follow_end_env <- trimws(Sys.getenv("UKBPPP_FOLLOW_END", unset = ""))
if (nzchar(follow_end_env)) {
  follow_end <- as.Date(follow_end_env)
  follow_end_source <- "UKBPPP_FOLLOW_END"
} else if (exists("date_follow_end", inherits = TRUE)) {
  follow_end <- as.Date(get("date_follow_end", inherits = TRUE))
  follow_end_source <- "upstream_date_follow_end"
} else {
  follow_end <- as.Date("2020-03-31")
  follow_end_source <- "Pradeep_operational_default"
}
if (length(follow_end) != 1L || is.na(follow_end)) {
  stop("Invalid follow-up cutoff. Set UKBPPP_FOLLOW_END as YYYY-MM-DD.", call. = FALSE)
}

required_follow_end_env <- trimws(Sys.getenv("UKBPPP_REQUIRE_FOLLOW_END", unset = ""))
if (nzchar(required_follow_end_env)) {
  required_follow_end <- as.Date(required_follow_end_env)
  if (is.na(required_follow_end) || !identical(follow_end, required_follow_end)) {
    stop(
      "FOLLOW_END_GATE_FAILED: resolved=", format(follow_end),
      "; required=", required_follow_end_env,
      "; source=", follow_end_source,
      call. = FALSE
    )
  }
} else {
  required_follow_end <- as.Date(NA)
}

analysis_population <- Sys.getenv("UKBPPP_POPULATION", unset = "White")
if (toupper(analysis_population) %in% c("ALL", "NONE", "")) analysis_population <- NA_character_

date_source <- Sys.getenv("UKBPPP_DATE_SOURCE", unset = "fod_ref")
outcome_date_sources <- c(
  cad = Sys.getenv("UKBPPP_DATE_SOURCE_CAD", unset = date_source),
  afib = Sys.getenv("UKBPPP_DATE_SOURCE_AF", unset = date_source),
  hfail = Sys.getenv("UKBPPP_DATE_SOURCE_HF", unset = date_source),
  ao_sten = Sys.getenv("UKBPPP_DATE_SOURCE_AS", unset = "fod_icd10")
)
max_proteins <- repro_int_env("UKBPPP_MAX_PROTEINS", 0L)
max_missing_protein <- as.numeric(Sys.getenv("UKBPPP_MAX_MISSING_PROTEIN", unset = "0.10"))
max_missing_individual <- as.numeric(Sys.getenv("UKBPPP_MAX_MISSING_INDIVIDUAL", unset = "0.10"))
skip_knn <- repro_bool_env("UKBPPP_SKIP_KNN", default = FALSE)
n_workers <- max(1L, repro_int_env("UKBPPP_WORKERS", 1L))
min_events_per_model <- max(10L, repro_int_env("UKBPPP_MIN_EVENTS", 10L))
use_original_15k_panel <- repro_bool_env("UKBPPP_ORIGINAL_15K_PANEL", default = TRUE)
protein_map_15k_file <- repro_env_path(
  "PRADEEP_PROTEIN_MAP",
  file.path(root_dir, "data.BIG", "gwas", "ppp", "map.raw", "olink_protein_map_1.5k_v1.tsv")
)
skip_relatedness <- repro_bool_env("UKBPPP_SKIP_RELATEDNESS", default = FALSE)
relatedness_cutoff <- as.numeric(Sys.getenv("UKBPPP_RELATEDNESS_CUTOFF", unset = "0.0884"))
if (!is.finite(relatedness_cutoff) || relatedness_cutoff <= 0) relatedness_cutoff <- 0.0884
kinship_file_env <- Sys.getenv("UKBPPP_KINSHIP_FILE", unset = "")
kinship_file_candidates <- unique(c(
  kinship_file_env,
  file.path(phe_dir, "common", "ukb7089_rel_s488363.dat"),
  file.path(phe_dir, "ukb7089_rel_s488363.dat"),
  file.path(root_dir, "data", "ukb", "ukb7089_rel_s488363.dat"),
  file.path(root_dir, "ukb7089_rel_s488363.dat")
))
kinship_file_candidates <- kinship_file_candidates[nzchar(kinship_file_candidates)]
kinship_file <- kinship_file_candidates[file.exists(kinship_file_candidates)][1]
if (is.na(kinship_file) && nzchar(kinship_file_env)) kinship_file <- kinship_file_env
if (length(kinship_file) == 0 || is.na(kinship_file)) kinship_file <- NA_character_
plate_well_file <- Sys.getenv("UKBPPP_PLATE_WELL_FILE", unset = "")
if (!nzchar(plate_well_file)) plate_well_file <- NA_character_

ppp_bed_env <- trimws(Sys.getenv("UKBPPP_PPP_BED_FILE", unset = ""))
ppp_bed_candidates <- unique(c(
  ppp_bed_env,
  file.path(root_dir, "data.BIG", "gwas", "ppp", "ppp_3k_b38.bed"),
  file.path(root_dir, "data.BIG", "gwas", "ppp", "ppp.b38.bed")
))
ppp_bed_candidates <- ppp_bed_candidates[nzchar(ppp_bed_candidates)]
ppp_bed_file <- ppp_bed_candidates[file.exists(ppp_bed_candidates)][1]
if (length(ppp_bed_file) == 0L || is.na(ppp_bed_file)) ppp_bed_file <- NA_character_

outcome_map <- data.table::data.table(
  outcome_key = c("cad", "afib", "hfail", "ao_sten"),
  label = c("CAD", "Afib", "HF", "AS"),
  local_suffix = c("cvd_cad", "cvd_afib", "cvd_hfail", "cvd_aosten"),
  paper_name = c("Coronary artery disease", "Atrial fibrillation", "Heart failure", "Aortic stenosis")
)

# The full four-outcome map remains the cohort-control and multiplicity family.
# UKBPPP_OUTCOME_SUBSET only controls which endpoints are fitted and emitted by
# Stages 2-4, so a CAD-only run preserves the paper cohort definition and its
# original four-outcome Bonferroni family.
analysis_outcome_keys <- trimws(unlist(strsplit(
  tolower(Sys.getenv("UKBPPP_OUTCOME_SUBSET", unset = "cad,afib,hfail,ao_sten")),
  "[,;+]"
)))
analysis_outcome_keys <- unique(analysis_outcome_keys[nzchar(analysis_outcome_keys)])
unknown_analysis_outcomes <- setdiff(analysis_outcome_keys, outcome_map$outcome_key)
if (length(analysis_outcome_keys) == 0L || length(unknown_analysis_outcomes) > 0L) {
  stop(
    "UKBPPP_OUTCOME_SUBSET contains unsupported outcomes: ",
    paste(unknown_analysis_outcomes, collapse = ","),
    call. = FALSE
  )
}
analysis_outcome_map <- outcome_map[match(analysis_outcome_keys, outcome_map$outcome_key)]

covariate_source_map <- c(
  Sex_numeric = "sex",
  age = "age",
  center = "center",
  mergedrace = "ethnic.c",
  ever_smoked = "smoke_ever",
  BMI = "bmi",
  SBP = "sbp",
  DBP = "dbp",
  tchol = "bb_TC",
  ldl = "bb_LDL",
  hdl = "bb_HDL",
  tg = "bb_TG",
  creat = "bb_CRE",
  tdi = "tdi",
  alc_int_num = "alcohol_freq",
  antihtnbase = "drug.htn",
  cholmed = "drug.lipid",
  dm2_prev = "dm.yes"
)

continuous_impute_vars <- c(
  "BMI", "tdi_log", "alc_int_num", "SBP", "DBP", "tchol", "ldl", "hdl", "tg", "creat"
)

clinical_covariates_base <- c(
  "age", "age2", "Sex_numeric", "mergedrace",
  paste0("PC", 1:10),
  "ever_smoked", "BMI_final", "SBP_final", "antihtnbase",
  "tchol_final", "hdl_final", "cholmed",
  "dm2_prev", "tdi_log_final", "creat_final"
)

clinical_predictors_for_lasso <- c(
  "age", "Sex_numeric", "mergedrace", "ever_smoked",
  "BMI_final", "SBP_final", "antihtnbase", "tchol_final",
  "hdl_final", "cholmed", "dm2_prev", "creat_final"
)

analysis_base_file <- file.path(output_dir, "ukbppp_cardiac_analysis_base.rds")
primary_assoc_file <- file.path(output_dir, "primary_association_cox.csv")
sex_interaction_file <- file.path(output_dir, "sex_interaction_cox.csv")
lasso_auc_file <- file.path(output_dir, "lasso_risk_score_auc.csv")

message("UKB-PPP cardiac reproduction config loaded")
message("  DIR0: ", root_dir)
message("  PHEDIR: ", phe_dir)
message("  SCRIPT_ROOT: ", script_root)
message("  HELPER_DIR: ", helper_dir)
message("  OUTPUT_ROOT: ", output_root)
message("  PROJECT_DIR: ", project_dir)
message("  analysis_project: ", analysis_project)
message("  analysis_dir: ", analysis_dir)
message("  population: ", ifelse(is.na(analysis_population), "ALL", analysis_population))
message("  follow_end: ", format(follow_end), " (", follow_end_source, ")")
message("  required_follow_end: ", ifelse(is.na(required_follow_end), "not set", format(required_follow_end)))
message("  date_source_default: ", date_source)
message("  date_source_by_outcome: ", paste(names(outcome_date_sources), outcome_date_sources, sep = "=", collapse = "; "))
message("  original_15k_panel: ", use_original_15k_panel)
message("  relatedness_exclusion: ", !skip_relatedness)
message("  ppp_bed_file: ", ifelse(is.na(ppp_bed_file), "not found", ppp_bed_file))
