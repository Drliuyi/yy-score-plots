args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))))
source(file.path(root, "R", "full", "00_core.R"))
source(file.path(root, "R", "full", "04_evaluate_figures.R"))

tmp <- Sys.getenv("YU_SMOKE_DIR")
if (!nzchar(tmp)) tmp <- tempfile("yu_figure_smoke_")
folders <- c(cohort = "04_cohort", cox = "05_cox", models = "09_models", evaluation = "10_evaluation",
             figures = "15_figures", supplement_figures = "16_supplementary_figures",
             source_tables = "03_source_tables", cmr = "07_cmr", mr = "11_mr",
             mediation = "12_mediation", prs = "13_prs", enrichment = "14_enrichment")
paths <- lapply(folders, function(x) file.path(tmp, x))
for (path in paths) dir.create(path, recursive = TRUE, showWarnings = FALSE)
cfg <- list(paths = paths)

outcomes <- data.table::fread(file.path(root, "config", "outcomes.csv"))
data.table::fwrite(
  data.table::data.table(step = c("phenotype", "incident", "derivation", "holdout"), n = c(53026, 46818, 31212, 15606)),
  file.path(paths$cohort, "cohort_flow.csv")
)
data.table::fwrite(
  outcomes[, .(outcome_id, derivation_events = seq(200, 1500, length.out = .N), test_events = seq(100, 750, length.out = .N))],
  file.path(paths$cohort, "endpoint_event_summary.csv")
)

set.seed(20260715)
cox <- outcomes[, {
  n <- 80L
  beta <- rnorm(n, 0, .18)
  p <- pmin(1, 10^(-runif(n, 0, 8)))
  data.table::data.table(
    feature_id = paste0("P", seq_len(n)), protein = paste0("GENE", seq_len(n)),
    beta = beta, p = p, hr = exp(beta), bonferroni_threshold = .05 / 80,
    bonferroni_significant = p < .05 / 80
  )
}, by = .(outcome_id, outcome_label)]
data.table::fwrite(cox, file.path(paths$cox, "table_s2_incident_associations.csv.gz"))
data.table::fwrite(
  data.table::data.table(feature_id = paste0("P", 1:80), protein = paste0("GENE", 1:80)),
  file.path(paths$cox, "retained_panel.csv")
)

models <- c("SCORE2", "Protein", "Protein_SCORE2")
metrics <- data.table::CJ(outcome_id = outcomes$outcome_id, model_id = models)
metrics[, auc_q50 := runif(.N, .62, .83)]
metrics[, `:=`(auc_q25 = auc_q50 - .015, auc_q75 = auc_q50 + .015)]
data.table::fwrite(metrics, file.path(paths$evaluation, "table_s10_prediction_metrics.csv"))
model_design <- data.table::CJ(outcome_id = outcomes$outcome_id, model_id = models)
model_design[, feature_n := data.table::fcase(
  model_id == "Protein", 35L,
  model_id == "Protein_SCORE2", 36L,
  default = 1L
)]
data.table::fwrite(model_design, file.path(paths$models, "model_design_contract.csv"))

importance <- data.table::CJ(outcome_id = outcomes$outcome_id, feature = paste0("P", 1:35))
importance[, model_id := "Protein"]
importance[, gain := rexp(.N)]
importance[, standardized_gain := gain / sum(gain), by = outcome_id]
data.table::fwrite(importance, file.path(paths$evaluation, "table_s12_model_importance.csv.gz"))

yur_build_figures(cfg)
required <- c(
  "figure2a_significant_counts.pdf", "figure2b_o_cox_volcano.tiff",
  "figure4a_prediction_auc.pdf", "figure4b_recurrent_importance.tiff"
)
missing <- required[!file.exists(file.path(paths$figures, required))]
stopifnot(!length(missing))
stopifnot(all(file.info(file.path(paths$figures, required))$size > 0))
if (!nzchar(Sys.getenv("YU_KEEP_SMOKE"))) {
  unlink(tmp, recursive = TRUE)
} else {
  cat("FIGURE_SMOKE_DIR=", tmp, "\n", sep = "")
}
cat("ALL FIGURE SMOKE TESTS PASSED\n")
