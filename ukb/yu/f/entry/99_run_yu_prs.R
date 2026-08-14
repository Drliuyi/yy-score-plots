#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
entry_dir <- normalizePath(dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])), winslash = "/")
project_dir <- normalizePath(file.path(entry_dir, "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(project_dir, "f", "R", "full", "00_core.R"))
source(file.path(project_dir, "f", "R", "full", "02_cohort.R"))
source(file.path(project_dir, "f", "R", "full", "03_cox_associations.R"))
source(file.path(project_dir, "f", "R", "full", "05_prs_associations.R"))

cli <- yur_parse_cli(args)
cfg <- yur_init_config(project_dir, cli)
cfg$endpoint_subset <- as.character(cli$endpoint_subset %||% "all")
cfg$workers <- as.integer(cli$workers %||% cfg$workers)

mode <- tolower(as.character(cli$mode %||% "help"))
if (mode == "help") {
  cat(paste0(
    "Yu/Chen local 13-outcome PRS-protein reconstruction\n\n",
    "Modes:\n",
    "  preflight          Freeze target EIDs, candidate proteins, manifests and expected row counts.\n",
    "  merge_scores       Merge 22 chromosome .sscore files into 13 x 5 participant PRSs.\n",
    "  associate_shard    Regress all candidate proteins for one --endpoint_subset.\n",
    "  merge_associations Validate all 13 shards; retain significant and non-significant rows.\n",
    "  report             Write a source-locked PRS reconstruction report.\n\n",
    "GWAS preparation and chromosome scoring are dispatched by yu.sh.\n"
  ))
  quit(status = 0)
}

yur_session_snapshot(cfg)
switch(mode,
  preflight = yur_run_stage(cfg, "prs_preflight", function() yur_prs_preflight(cfg)),
  merge_scores = yur_run_stage(cfg, "prs_merge_scores", function() yur_prs_merge_scores(cfg)),
  associate_shard = yur_run_stage(
    cfg, paste0("prs_associate_", cfg$endpoint_subset), function() yur_prs_associate_shard(cfg)
  ),
  merge_associations = yur_run_stage(cfg, "prs_merge_associations", function() yur_prs_merge_associations(cfg)),
  report = yur_run_stage(cfg, "prs_report", function() yur_prs_report(cfg)),
  stop("Unknown PRS mode: ", mode, call. = FALSE)
)
