args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))))

source(file.path(root, "R", "full", "00_core.R"))
source(file.path(root, "R", "full", "01_sources_preflight.R"))
source(file.path(root, "R", "full", "02_cohort.R"))
source(file.path(root, "R", "full", "03_cox_associations.R"))
source(file.path(root, "R", "full", "03b_cmr_associations.R"))
source(file.path(root, "R", "full", "04_evaluate_figures.R"))

defaults <- jsonlite::read_json(file.path(root, "config", "full_reproduction_defaults.json"), simplifyVector = TRUE)
stopifnot(identical(defaults$prediction_panel_mode, "local_reselected"))
stopifnot(
  identical(defaults$figure4_extra_projects, ""),
  identical(defaults$figure4_extra_outcomes, ""),
  identical(defaults$figure4_extra_labels, "")
)

tc_good <- yur_numeric_source_summary(c(3.5, 5.5, 7.5), "total_cholesterol", "bb_TC", TRUE, defaults)
tc_bad <- yur_numeric_source_summary(c(90, 107, 130), "total_cholesterol", "total_cholesterol", FALSE, defaults)
hdl_good <- yur_numeric_source_summary(c(.8, 1.4, 2.2), "hdl", "bb_HDL", TRUE, defaults)
stopifnot(tc_good$status == "PASS", tc_bad$status == "FAIL_IMPLAUSIBLE_RANGE", hdl_good$status == "PASS")

score2_synthetic <- data.table::data.table(
  total_cholesterol = seq(3.2, 8.2, length.out = 220),
  hdl = seq(.7, 2.4, length.out = 220), sbp = seq(100, 180, length.out = 220),
  score2_raw = seq(.01, .25, length.out = 220)
)
score2_qc_good <- yur_score2_distribution_qc(score2_synthetic, score2_synthetic, defaults)
stopifnot(all(score2_qc_good$status == "PASS"))
score2_synthetic[, score2_raw := 1]
score2_qc_bad <- yur_score2_distribution_qc(score2_synthetic, score2_synthetic, defaults)
stopifnot(all(score2_qc_bad[variable == "score2_raw", status] == "FAIL_SCORE2_SATURATION"))

eid_cache_file <- tempfile(fileext = ".csv.gz")
data.table::fwrite(data.table::data.table(eid = c("1001", "1002"), fasting_time = c(2, 4)), eid_cache_file)
eid_cache <- data.table::fread(eid_cache_file)
stopifnot(!is.character(eid_cache$eid))
eid_cache <- yur_normalize_eid_column(eid_cache)
eid_master <- data.table::data.table(eid = c("1001", "1002"), value = 1:2)
eid_joined <- merge(eid_master, eid_cache, by = "eid", all.x = TRUE)
stopifnot(is.character(eid_cache$eid), identical(eid_joined$fasting_time, c(2L, 4L)))
unlink(eid_cache_file)

cohort_test_dir <- tempfile(pattern = "cohort_eid_test_")
dir.create(cohort_test_dir)
data.table::fwrite(
  data.table::data.table(eid = c("2001", "2002"), event_cad = c(0L, 1L)),
  file.path(cohort_test_dir, "derivation_cohort.csv.gz")
)
data.table::fwrite(
  data.table::data.table(eid = c("2003", "2004"), event_cad = c(1L, 0L)),
  file.path(cohort_test_dir, "test_cohort.csv.gz")
)
cohort_loaded <- yur_load_cox_cohorts(list(paths = list(cohort = cohort_test_dir)))
stopifnot(
  is.character(cohort_loaded$derivation$eid),
  is.character(cohort_loaded$test$eid),
  is.character(cohort_loaded$all_meta$eid),
  identical(cohort_loaded$all_meta$eid, c("2001", "2002", "2003", "2004"))
)
unlink(cohort_test_dir, recursive = TRUE)

official_outcomes <- c(
  "abdominal_aneurysm", "atrial_fibrillation", "aortic_valve_stenosis", "cad",
  "cardiomyopathy", "deep_vein_thrombosis", "heart_failure",
  "intracerebral_hemorrhage", "ischemic_stroke", "peripheral_arterial_disease",
  "pulmonary_embolism", "subarachnoid_hemorrhage", "thoracic_aneurysm",
  "transient_ischemic_attack"
)
outcomes <- data.table::fread(file.path(root, "config", "outcomes.csv"))
stopifnot(identical(sort(outcomes$outcome_id), sort(official_outcomes)))
stopifnot(
  identical(
    outcomes[outcome_id == "abdominal_aneurysm", local_date_candidates],
    "fod_icd10_cvd_abaneu;fod_icd9_cvd_abaneu"
  ),
  identical(
    outcomes[outcome_id == "thoracic_aneurysm", local_date_candidates],
    "fod_icd10_cvd_thaneu;fod_icd9_cvd_thaneu"
  ),
  identical(
    outcomes[outcome_id == "deep_vein_thrombosis", local_date_candidates],
    "fod_icd10_cvd_dvt;fod_icd9_cvd_dvt"
  ),
  !any(grepl("fod_ref_cvd_(abaneu|thaneu|dvt)", outcomes$local_date_candidates))
)

death_test <- data.table::data.table(
  death_icd10 = c("I714", "I219", NA_character_, "I802"),
  death_icd10sec = c(NA_character_, "I713 I10", "I712", NA_character_),
  death_date = as.Date(c("2020-01-01", "2020-02-01", "2020-03-01", "2020-04-01"))
)
aaa_death <- yur_death_date_for_codes(death_test, c("I713", "I714"))
taa_death <- yur_death_date_for_codes(death_test, c("I711", "I712"))
dvt_death <- yur_death_date_for_codes(death_test, "I802")
stopifnot(
  identical(which(!is.na(aaa_death)), c(1L, 2L)),
  identical(which(!is.na(taa_death)), 3L),
  identical(which(!is.na(dvt_death)), 4L)
)

duplicate_test <- data.table::data.table(
  date_a = as.Date(c("2020-01-01", NA, "2021-01-01")),
  date_b = as.Date(c("2020-01-01", NA, "2021-01-01")),
  date_c = as.Date(c("2020-01-02", NA, "2021-01-01"))
)
duplicate_qc <- yur_outcome_pairwise_qc(duplicate_test, c("a", "b", "c"))
stopifnot(
  duplicate_qc[outcome_a == "a" & outcome_b == "b", exact_duplicate],
  !duplicate_qc[outcome_a == "a" & outcome_b == "c", exact_duplicate]
)

event_test <- data.table::data.table(
  outcome_id = c("a", "b"), derivation_n = 100L, derivation_events = c(25L, 19L),
  test_n = 50L, test_events = c(12L, 9L)
)
event_rules <- data.table::data.table(
  outcome_id = c("a", "b"), min_derivation_events = 20L,
  min_test_events = 10L, min_total_events = 30L
)
event_qc <- yur_validate_endpoint_events(event_test, event_rules)
stopifnot(event_qc[outcome_id == "a", status] == "PASS")
stopifnot(event_qc[outcome_id == "b", status] == "FAIL_LOW_EVENTS")

contract_dir <- tempfile(pattern = "cox_contract_test_")
dir.create(file.path(contract_dir, "cohort"), recursive = TRUE)
dir.create(file.path(contract_dir, "cox"), recursive = TRUE)
dir.create(file.path(contract_dir, "preflight"), recursive = TRUE)
dir.create(file.path(contract_dir, "logs"), recursive = TRUE)
writeLines("cohort-hash", file.path(contract_dir, "cohort", "cohort_contract_hash.txt"))
writeLines("definitions-hash", file.path(contract_dir, "cohort", "outcome_definition_hash.txt"))
writeLines("panel-hash", file.path(contract_dir, "cox", "retained_panel_hash.txt"))
contract_cfg <- list(
  endpoint_subset = "cad",
  paths = list(
    cohort = file.path(contract_dir, "cohort"), cox = file.path(contract_dir, "cox"),
    preflight = file.path(contract_dir, "preflight"), logs = file.path(contract_dir, "logs"),
    run_log = file.path(contract_dir, "logs", "run.log")
  )
)
contract_outcome <- data.table::data.table(outcome_id = "cad", definition_hash = "cad-definition-hash")
contract_expected <- yur_cox_contract_values(contract_cfg, contract_outcome)
contract_file <- file.path(contract_dir, "contract.json")
yur_write_json(contract_expected, contract_file)
stopifnot(yur_cox_contract_matches(contract_file, contract_expected))
contract_stale <- contract_expected
contract_stale$cohort_contract_hash <- "different-cohort"
stopifnot(!yur_cox_contract_matches(contract_file, contract_stale))

data.table::fwrite(contract_outcome, file.path(contract_dir, "preflight", "outcome_field_resolution.csv"))
data.table::fwrite(
  data.table::data.table(feature_id = c("P1", "P2")),
  file.path(contract_dir, "cox", "retained_panel.csv")
)
for (scope in c("full_incident", "derivation")) {
  shard_file <- file.path(contract_dir, "cox", paste0(scope, "_cad_cox.csv.gz"))
  shard_contract <- file.path(contract_dir, "cox", paste0(scope, "_cad_cox.contract.json"))
  data.table::fwrite(data.table::data.table(feature_id = c("P1", "P2")), shard_file)
  yur_write_json(contract_expected, shard_contract)
}
stopifnot(yur_cox_shard_complete(contract_cfg))
unlink(file.path(contract_dir, "cox", "full_incident_cad_cox.csv.gz"))
stopifnot(!yur_cox_shard_complete(contract_cfg))
unlink(contract_dir, recursive = TRUE)

mixed_integer_covariates <- data.table::data.table(
  age = c(60L, 61L, 62L, 63L, 64L),
  protein_sampling_lag_days = c(10L, 11L, 12L, 13L, 14L),
  fasting_time = c(4L, 5L, 6L, 7L, 8L),
  bmi = c(25L, NA, 27L, 28L, NA), sbp = c(120L, 121L, NA, 123L, 124L),
  tdi = c(1L, NA, 2L, 3L, NA), sex = c("F", "F", "F", "M", "M"),
  assessment_center = c(1L, 1L, 2L, 2L, 2L), ethnicity_white = c(1L, 1L, 1L, 0L, 0L),
  smoking = c("No", "No", "Yes", "No", "Yes"),
  alcohol = c("No", "Yes", "Yes", "No", "Yes"),
  blood_collection_season = c("Winter", "Spring", "Summer", "Autumn", "Winter")
)
imputed_covariates <- yur_impute_association_covariates(mixed_integer_covariates)
stopifnot(
  all(vapply(imputed_covariates[, .(age, protein_sampling_lag_days, fasting_time, bmi, sbp, tdi)], is.double, logical(1))),
  !anyNA(imputed_covariates$bmi), !anyNA(imputed_covariates$sbp), !anyNA(imputed_covariates$tdi)
)

stage_dir <- tempfile(pattern = "resume_validator_test_")
dir.create(stage_dir, recursive = TRUE)
stage_cfg <- list(
  resume = TRUE, force = FALSE,
  paths = list(logs = stage_dir, run_log = file.path(stage_dir, "run.log"))
)
yur_write_json(list(stage = "synthetic", status = "PASS"), file.path(stage_dir, "synthetic.done.json"))
stage_ran <- FALSE
yur_run_stage(stage_cfg, "synthetic", function() stage_ran <<- TRUE, validate_done = function() FALSE)
stopifnot(stage_ran)
stage_ran <- FALSE
yur_run_stage(stage_cfg, "synthetic", function() stage_ran <<- TRUE, validate_done = function() TRUE)
stopifnot(!stage_ran)
unlink(stage_dir, recursive = TRUE)

season <- as.character(yur_blood_collection_season(as.Date(c(
  "2026-01-15", "2026-04-15", "2026-07-15", "2026-10-15"
))))
stopifnot(identical(season, c("Winter", "Spring", "Summer", "Autumn")))

processing_file <- tempfile(fileext = ".dat")
data.table::fwrite(data.table::data.table(
  PlateID = c("890000000001", "890000000001"),
  Panel = c("Cardiometabolic", "Cardiometabolic II"),
  Processing_StartDate = c("2021-04-28", "2021-05-03")
), processing_file, sep = "\t")
technical <- yur_derive_technical_covariates(data.table::data.table(
  eid = "1", baseline_date = as.Date("2021-04-01"), protein_plate = "890000000001"
), processing_file)
stopifnot(
  technical$values$protein_sampling_lag_days__cardiometabolic == 27,
  technical$values$protein_sampling_lag_days__cardiometabolic_ii == 32,
  technical$values$protein_sampling_lag_days == 29.5
)
unlink(processing_file)

set.seed(11)
n <- 450L
protein <- scale(rnorm(n))[, 1]
event_time <- rexp(n, rate = exp(.35 * protein) / 1200)
censor_time <- rexp(n, rate = 1 / 1700)
event <- as.integer(event_time <= censor_time)
time <- pmax(1, pmin(event_time, censor_time))
covariates <- cbind(age = rnorm(n), sex = rbinom(n, 1, .5))
fit <- yur_cox_one(protein, time, event, covariates)
stopifnot(is.finite(fit[["beta"]]), is.finite(fit[["se"]]), fit[["events"]] >= 5)
fit_row <- yur_cox_result_row(
  feature_id = "TEST_PROTEIN", protein_mean = 0, protein_sd = 1,
  lag_variable = "protein_sampling_lag_days", lag_panel_fallback = FALSE,
  stat = fit
)
stopifnot(
  nrow(fit_row) == 1L,
  all(c("n", "events", "beta", "se", "z", "p", "hr", "ci_low", "ci_high") %in% names(fit_row)),
  is.finite(fit_row$p)
)

cmr_metrics <- data.table::fread(file.path(root, "config", "cmr_metrics.csv"))
stopifnot(
  nrow(cmr_metrics) == 19L,
  !anyDuplicated(cmr_metrics$metric_id),
  !anyDuplicated(cmr_metrics$Outcome),
  identical(as.integer(table(cmr_metrics$anatomy)), c(4L, 7L, 4L, 4L)),
  identical(names(table(cmr_metrics$anatomy)), c("LA", "LV", "RA", "RV"))
)
cmr_column_selection_test <- data.table::as.data.table(
  setNames(replicate(19L, c(1, NA_real_), simplify = FALSE), cmr_metrics$source_column)
)
cmr_metric_columns <- cmr_metrics$source_column
stopifnot(identical(
  rowSums(is.finite(as.matrix(cmr_column_selection_test[, ..cmr_metric_columns]))),
  c(19, 0)
))
set.seed(20260718)
cmr_n <- 600L
cmr_covariates <- cbind(
  age = scale(rnorm(cmr_n, 62, 7))[, 1],
  male = rbinom(cmr_n, 1, .48),
  bmi = scale(rnorm(cmr_n, 27, 4))[, 1]
)
cmr_protein <- rnorm(cmr_n)
cmr_outcome <- .28 * cmr_protein + .12 * cmr_covariates[, "age"] +
  .18 * cmr_covariates[, "male"] + rnorm(cmr_n)
cmr_protein[c(5, 17, 44)] <- NA_real_
cmr_outcome[c(8, 29)] <- NA_real_
cmr_stat <- yur_linear_one(cmr_protein, cmr_outcome, cmr_covariates)
cmr_ok <- is.finite(cmr_protein) & is.finite(cmr_outcome) &
  rowSums(!is.finite(cmr_covariates)) == 0L
cmr_reference <- summary(lm(
  scale(cmr_outcome[cmr_ok])[, 1] ~ scale(cmr_protein[cmr_ok])[, 1] +
    cmr_covariates[cmr_ok, "age"] + cmr_covariates[cmr_ok, "male"] +
    cmr_covariates[cmr_ok, "bmi"]
))$coefficients[2, ]
stopifnot(
  cmr_stat[["n"]] == sum(cmr_ok),
  isTRUE(all.equal(unname(cmr_stat[["beta"]]), unname(cmr_reference[["Estimate"]]), tolerance = 1e-10)),
  isTRUE(all.equal(unname(cmr_stat[["se"]]), unname(cmr_reference[["Std. Error"]]), tolerance = 1e-10)),
  isTRUE(all.equal(unname(cmr_stat[["p"]]), unname(cmr_reference[["Pr(>|t|)"]]), tolerance = 1e-10))
)

y <- c(rep(0L, 100), rep(1L, 100))
p_good <- c(seq(.01, .45, length.out = 100), seq(.55, .99, length.out = 100))
p_bad <- rev(p_good)
stopifnot(yur_auc(y, p_good) > .99, yur_auc(y, p_bad) < .01)
metrics <- yur_binary_metrics(y, p_good, .5)
stopifnot(all(c("auc", "accuracy", "sensitivity", "specificity", "f1", "brier") %in% names(metrics)))

cat("ALL FULL-REPRODUCTION R TESTS PASSED\n")
