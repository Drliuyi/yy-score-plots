#!/usr/bin/env Rscript

# The macOS PDF device does not expose Arial by default. Register a metric-
# compatible alias so publication figures render identically across Mac and
# Windows instead of failing before the PNG/TIFF outputs are written.
if (identical(Sys.info()[["sysname"]], "Darwin") && !"Arial" %in% names(grDevices::pdfFonts())) {
  grDevices::pdfFonts(Arial = grDevices::pdfFonts("Helvetica")[[1]])
}

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
entry_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  getwd()
}
project_dir <- normalizePath(file.path(entry_dir, "..", ".."), winslash = "/", mustWork = TRUE)

source(file.path(project_dir, "f", "R", "full", "00_core.R"))
source(file.path(project_dir, "f", "R", "full", "01_sources_preflight.R"))
source(file.path(project_dir, "f", "R", "full", "02_cohort.R"))
source(file.path(project_dir, "f", "R", "full", "03_cox_associations.R"))
source(file.path(project_dir, "f", "R", "full", "03b_cmr_associations.R"))
source(file.path(project_dir, "f", "R", "full", "04_evaluate_figures.R"))
source(file.path(project_dir, "f", "R", "full", "06_mr.R"))
source(file.path(project_dir, "f", "R", "full", "07_mediation.R"))
source(file.path(project_dir, "f", "R", "full", "08_systems_biology.R"))

cli <- yur_parse_cli(commandArgs(trailingOnly = TRUE))
cfg <- yur_init_config(project_dir, cli)

if (cfg$mode == "help") {
  yur_print_help()
  quit(save = "no", status = 0)
}

if (cfg$mode %in% c("all", "all_fast")) {
  stop(
    "Mode=all must be dispatched by yu.sh ",
    "because the workflow alternates between R and Python stages.",
    call. = FALSE
  )
}

yur_session_snapshot(cfg)
yur_log(cfg, "Mode=", cfg$mode, " endpoint_subset=", cfg$endpoint_subset)

run <- function(stage, fun, validate_done = NULL) {
  yur_run_stage(cfg, stage, fun, validate_done = validate_done)
}
if (cfg$mode == "sources") run("sources", function() yur_sources(cfg))
if (cfg$mode == "preflight") run(
  "preflight", function() yur_preflight(cfg), validate_done = function() yur_preflight_complete(cfg)
)
if (cfg$mode == "cohort") run(
  "cohort", function() yur_build_cohort(cfg), validate_done = function() yur_cohort_complete(cfg)
)
if (cfg$mode == "cox") run("cox", function() yur_run_cox(cfg))
if (cfg$mode == "cox_prepare") run("cox_prepare", function() yur_prepare_cox_panel(cfg))
if (cfg$mode == "cox_shard") run(
  paste0("cox_shard_", gsub("[^A-Za-z0-9]+", "_", cfg$endpoint_subset)),
  function() yur_run_cox_shard(cfg),
  validate_done = function() yur_cox_shard_complete(cfg)
)
if (cfg$mode == "cox_merge") run("cox_merge", function() yur_merge_cox(cfg))
if (cfg$mode == "cmr_prepare") run(
  "cmr_prepare", function() yur_prepare_cmr(cfg), validate_done = function() yur_cmr_prepare_complete(cfg)
)
if (cfg$mode == "cmr_shard") run(
  paste0("cmr_shard_", gsub("[^A-Za-z0-9]+", "_", cfg$cmr_metric_subset)),
  function() yur_run_cmr_shard(cfg), validate_done = function() yur_cmr_shard_complete(cfg)
)
if (cfg$mode == "cmr_merge") run("cmr_merge", function() yur_merge_cmr(cfg))
if (cfg$mode == "cmr") run("cmr", function() yur_run_cmr(cfg))
if (cfg$mode == "mr_prepare") run("mr_prepare", function() yur_prepare_local_mr(cfg))
if (cfg$mode == "mr_run") run("mr_run", function() yur_run_local_mr(cfg))
if (cfg$mode == "mediation_prepare") run("mediation_prepare", function() yur_prepare_local_mediation(cfg))
if (cfg$mode == "mediation_run") run("mediation_run", function() yur_run_local_mediation(cfg))
if (cfg$mode == "mediation_cmest_pilot") {
  run("mediation_cmest_pilot", function() yur_run_cmest_pilot(cfg))
}
if (cfg$mode == "mediation_cmest_shard") {
  shard_name <- sprintf(
    "mediation_cmest_shard_%03d_of_%03d", cfg$cmest_shard_index, cfg$cmest_shard_count
  )
  run(
    shard_name, function() yur_run_cmest_shard(cfg),
    validate_done = function() yur_cmest_shard_complete(cfg)
  )
}
if (cfg$mode == "mediation_cmest_merge") {
  run("mediation_cmest_merge", function() yur_merge_cmest_shards(cfg))
}
if (cfg$mode == "systems_prepare") run("systems_prepare", function() yur_prepare_figure6_systems(cfg))
if (cfg$mode == "systems_enrichment") run("systems_enrichment", function() yur_run_figure6_enrichment(cfg))
if (cfg$mode == "systems_tf") run("systems_tf", function() yur_run_figure6_tf(cfg))
if (cfg$mode == "systems_ppi") run("systems_ppi", function() yur_run_figure6_ppi(cfg))
if (cfg$mode == "systems_figures") run("systems_figures", function() yur_plot_figure6_systems(cfg))
if (cfg$mode == "figure6_systems") run("figure6_systems", function() yur_run_figure6_systems_all(cfg))
if (cfg$mode == "evaluate") run("evaluate", function() yur_evaluate_prediction(cfg))
if (cfg$mode == "figures") run("figures", function() yur_build_figures(cfg))
if (cfg$mode == "report") run("report", function() yur_build_report(cfg))

python_modes <- c("select", "train")
if (cfg$mode %in% python_modes) {
  stop("Mode=", cfg$mode, " is dispatched by yu.sh", call. = FALSE)
}

valid <- c("sources", "preflight", "cohort", "cox_prepare", "cox_shard", "cox_merge", "cox",
           "cmr_prepare", "cmr_shard", "cmr_merge", "cmr",
           "mr_prepare", "mr_run", "mediation_prepare", "mediation_run",
           "mediation_cmest_pilot", "mediation_cmest_shard", "mediation_cmest_merge",
           "systems_prepare", "systems_enrichment", "systems_tf", "systems_ppi",
           "systems_figures", "figure6_systems",
           "select", "train", "evaluate", "figures", "report", "all", "all_fast")
if (!cfg$mode %in% valid) stop("Unsupported mode: ", cfg$mode, call. = FALSE)
