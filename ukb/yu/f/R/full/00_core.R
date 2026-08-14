suppressWarnings(suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(digest)
}))

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

yur_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
yur_bool <- function(x) isTRUE(x) || tolower(as.character(x %||% "false")) %in% c("1", "true", "yes", "y")
yur_norm_eid <- function(x) sub("\\.0$", "", trimws(as.character(x)))
yur_normalize_eid_column <- function(x, column = "eid") {
  if (!data.table::is.data.table(x)) x <- data.table::as.data.table(x)
  if (!column %in% names(x)) stop("Missing EID column: ", column, call. = FALSE)
  x[, (column) := yur_norm_eid(get(column))]
  x
}
yur_norm_name <- function(x) gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", tolower(trimws(as.character(x)))))

yur_parse_cli <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) stop("Unexpected argument: ", token, call. = FALSE)
    token <- sub("^--", "", token)
    if (grepl("=", token, fixed = TRUE)) {
      bits <- strsplit(token, "=", fixed = TRUE)[[1]]
      key <- bits[[1]]
      value <- paste(bits[-1], collapse = "=")
    } else {
      key <- token
      if (i < length(args) && !startsWith(args[[i + 1L]], "--")) {
        i <- i + 1L
        value <- args[[i]]
      } else {
        value <- TRUE
      }
    }
    out[[gsub("-", "_", key, fixed = TRUE)]] <- value
    i <- i + 1L
  }
  out
}

yur_abs_path <- function(path, dir0, project_dir = NULL) {
  path <- as.character(path %||% "")
  if (!nzchar(path)) return("")
  if (grepl("^([A-Za-z]:[/\\\\]|/)", path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }
  root <- if (!is.null(project_dir) && startsWith(path, "config/")) {
    file.path(project_dir, "f")
  } else if (!is.null(project_dir) && startsWith(path, "references/")) {
    project_dir
  } else {
    dir0
  }
  normalizePath(file.path(root, path), winslash = "/", mustWork = FALSE)
}

yur_init_config <- function(project_dir, cli) {
  cfg_file <- file.path(project_dir, "f", "config", "full_reproduction_defaults.json")
  cfg <- read_json(cfg_file, simplifyVector = TRUE)
  for (nm in names(cli)) cfg[[nm]] <- cli[[nm]]
  cfg$project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  cfg$dir0 <- normalizePath(cfg$dir0 %||% getwd(), winslash = "/", mustWork = FALSE)
  cfg$mode <- tolower(as.character(cfg$mode %||% "help"))
  cfg$resume <- yur_bool(cfg$resume %||% FALSE)
  cfg$force <- yur_bool(cfg$force %||% FALSE)
  for (v in c(
    "workers", "cox_parallel_jobs", "cmr_parallel_jobs", "bootstrap_n", "bootstrap_seed",
    "split_seed", "inner_fold_seed", "inner_folds",
    "cmest_shard_index", "cmest_shard_count", "cmest_pilot_nboot",
    "string_required_score", "systems_top_n_per_outcome", "systems_max_tf"
  )) {
    cfg[[v]] <- as.integer(cfg[[v]])
  }
  cfg$protein_missingness_max <- as.numeric(cfg$protein_missingness_max)
  cfg$importance_cumulative_fraction <- as.numeric(cfg$importance_cumulative_fraction)
  cfg$systems_enrichment_fdr <- as.numeric(cfg$systems_enrichment_fdr %||% 0.05)
  cfg$figure4_extra_projects <- as.character(cfg$figure4_extra_projects %||% "")
  cfg$figure4_extra_outcomes <- as.character(cfg$figure4_extra_outcomes %||% "")
  cfg$figure4_extra_labels <- as.character(cfg$figure4_extra_labels %||% "")
  cfg$phenotype_rds <- yur_abs_path(cfg$phenotype_rds, cfg$dir0)
  cfg$raw_phenotype_file <- yur_abs_path(cfg$raw_phenotype_file, cfg$dir0, cfg$project_dir)
  cfg$raw_protein_file <- yur_abs_path(cfg$raw_protein_file, cfg$dir0)
  cfg$cmr_feature_file <- yur_abs_path(cfg$cmr_feature_file, cfg$dir0)
  cfg$panel_mapping_file <- yur_abs_path(cfg$panel_mapping_file, cfg$dir0)
  cfg$pqtl_root <- yur_abs_path(cfg$pqtl_root %||% "ppp/clean", cfg$dir0)
  cfg$supplement_workbook_file <- yur_abs_path(
    cfg$supplement_workbook_file %||% "references/raw/pwaf072_supplementary_table_1.xlsx",
    cfg$dir0, cfg$project_dir
  )
  cfg$supplement_methods_file <- yur_abs_path(
    cfg$supplement_methods_file %||% "references/raw/pwaf072_supplementary_figure_1.pdf",
    cfg$dir0, cfg$project_dir
  )
  cfg$olink_processing_start_date_file <- yur_abs_path(
    cfg$olink_processing_start_date_file, cfg$dir0, cfg$project_dir
  )
  cfg$outcomes_file <- file.path(cfg$project_dir, "f", "config", "outcomes.csv")
  cfg$field_map_file <- file.path(cfg$project_dir, "f", "config", "local_field_map.csv")
  cfg$method_provenance_file <- file.path(cfg$project_dir, "f", "config", "method_provenance.csv")
  cfg$cmr_metrics_file <- file.path(cfg$project_dir, "f", "config", "cmr_metrics.csv")
  cfg$analysis_root <- yur_abs_path(
    cfg$analysis_root %||% file.path(cfg$dir0, "analysis"), cfg$dir0
  )
  cfg$analysis_dir <- normalizePath(
    file.path(cfg$analysis_root, cfg$analysis_project), winslash = "/", mustWork = FALSE
  )
  folders <- c(
    logs = "00_logs", sources = "01_sources", preflight = "02_preflight",
    source_tables = "03_source_tables", cohort = "04_cohort",
    cox = "05_cox", sensitivity = "06_sensitivity",
    cmr = "07_cmr", selection = "08_selection", models = "09_models",
    evaluation = "10_evaluation", mr = "11_mr", mediation = "12_mediation",
    prs = "13_prs", enrichment = "14_enrichment", figures = "15_figures",
    supplement_figures = "16_supplementary_figures", report = "17_report", cache = "90_cache"
  )
  cfg$paths <- lapply(folders, function(x) file.path(cfg$analysis_dir, x))
  for (p in cfg$paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  cfg$paths$run_log <- file.path(cfg$paths$logs, "run.log")
  cfg$mr_outcome_lookup_dir <- yur_abs_path(
    cfg$mr_outcome_lookup_dir %||% file.path("analysis", cfg$analysis_project, "11_mr", "outcome_lookup"),
    cfg$dir0
  )
  cfg
}

yur_log <- function(cfg, ..., level = "INFO") {
  line <- sprintf("%s | %s | %s", yur_now(), level, paste(..., collapse = ""))
  cat(line, "\n")
  cat(line, "\n", file = cfg$paths$run_log, append = TRUE)
  invisible(line)
}

yur_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(as.data.table(x), path, na = "")
  invisible(path)
}

yur_write_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_json(x, path, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null")
  invisible(path)
}

yur_sha_file <- function(path) {
  if (file.exists(path)) digest(path, algo = "sha256", file = TRUE, serialize = FALSE) else NA_character_
}

yur_sha_text <- function(x) digest(paste(as.character(x), collapse = "\n"), algo = "sha256", serialize = FALSE)

yur_run_stage <- function(cfg, stage, fun, validate_done = NULL) {
  done_file <- file.path(cfg$paths$logs, paste0(stage, ".done.json"))
  error_file <- file.path(cfg$paths$logs, paste0(stage, ".error.json"))
  if (cfg$resume && file.exists(done_file) && !cfg$force) {
    validated <- is.null(validate_done) || isTRUE(tryCatch(
      validate_done(),
      error = function(e) {
        yur_log(cfg, "Resume validation failed for ", stage, ": ", conditionMessage(e), level = "WARN")
        FALSE
      }
    ))
    if (validated) {
      yur_log(cfg, "Resume: skip ", stage)
      return(invisible(NULL))
    }
    yur_log(cfg, "Resume marker is stale/incomplete; rerun ", stage, level = "WARN")
  }
  started <- Sys.time()
  yur_log(cfg, "START ", stage)
  tryCatch({
    fun()
    if (file.exists(error_file)) unlink(error_file)
    yur_write_json(list(stage = stage, status = "PASS", started = format(started),
                        ended = yur_now(), elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))),
                   done_file)
    yur_log(cfg, "DONE ", stage)
  }, error = function(e) {
    yur_write_json(list(stage = stage, status = "ERROR", message = conditionMessage(e), ended = yur_now()),
                   error_file)
    stop(e)
  })
}

yur_as_date <- function(x) {
  if (inherits(x, "IDate")) return(as.Date(as.character(x)))
  if (inherits(x, "Date")) return(as.Date(x))
  suppressWarnings(as.Date(x))
}

yur_blood_collection_season <- function(x) {
  d <- yur_as_date(x)
  month <- as.integer(format(d, "%m"))
  out <- rep(NA_character_, length(month))
  out[month %in% c(12L, 1L, 2L)] <- "Winter"
  out[month %in% 3:5] <- "Spring"
  out[month %in% 6:8] <- "Summer"
  out[month %in% 9:11] <- "Autumn"
  factor(out, levels = c("Winter", "Spring", "Summer", "Autumn"))
}

yur_plate_character <- function(x) {
  if (inherits(x, "integer64")) {
    if (!requireNamespace("bit64", quietly = TRUE)) {
      stop("Package bit64 is required to decode Olink plate identifiers.", call. = FALSE)
    }
    x <- as.character(x)
  } else {
    x <- as.character(x)
  }
  x <- trimws(x)
  x <- sub("\\.0+$", "", x)
  x[x %in% c("", "NA", "NaN", "<NA>")] <- NA_character_
  x
}

yur_panel_key <- function(x) {
  out <- tolower(trimws(as.character(x)))
  out <- gsub("[[:space:]-]+", "_", out)
  out <- gsub("_+", "_", out)
  out
}

yur_read_processing_dates <- function(path) {
  x <- fread(path, showProgress = FALSE)
  required <- c("PlateID", "Panel", "Processing_StartDate")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("Olink processing-date resource is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  x[, plate_id := yur_plate_character(PlateID)]
  x[, panel_key := yur_panel_key(Panel)]
  x[, processing_date := yur_as_date(Processing_StartDate)]
  x <- x[!is.na(plate_id) & !is.na(panel_key) & !is.na(processing_date)]
  x[, .(processing_date = min(processing_date)), by = .(plate_id, panel_key)]
}

yur_derive_technical_covariates <- function(x, processing_file) {
  required <- c("eid", "baseline_date", "protein_plate")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Cannot derive technical covariates; missing: ", paste(missing, collapse = ", "), call. = FALSE)
  participant <- x[, .(
    eid = yur_norm_eid(eid),
    baseline_date = yur_as_date(baseline_date),
    plate_id = yur_plate_character(protein_plate)
  )]
  processing <- yur_read_processing_dates(processing_file)
  joined <- merge(participant, processing, by = "plate_id", all.x = TRUE, allow.cartesian = TRUE)
  joined[, lag_days := as.numeric(processing_date - baseline_date)]
  wide <- dcast(joined, eid ~ panel_key, value.var = "lag_days", fun.aggregate = function(z) {
    z <- z[is.finite(z)]
    if (length(z)) median(z) else NA_real_
  })
  panel_columns <- setdiff(names(wide), "eid")
  if (length(panel_columns)) {
    setnames(wide, panel_columns, paste0("protein_sampling_lag_days__", panel_columns))
    lag_columns <- setdiff(names(wide), "eid")
    wide[, protein_sampling_lag_days := apply(.SD, 1L, function(z) {
      z <- as.numeric(z)
      z <- z[is.finite(z)]
      if (length(z)) median(z) else NA_real_
    }), .SDcols = lag_columns]
  } else {
    wide[, protein_sampling_lag_days := NA_real_]
  }
  setkey(wide, eid)
  setkey(participant, eid)
  coverage <- data.table(
    participants = nrow(participant),
    plate_nonmissing = sum(!is.na(participant$plate_id)),
    plate_matched_any_panel = uniqueN(joined[!is.na(processing_date), eid]),
    panel_lag_columns = sum(startsWith(names(wide), "protein_sampling_lag_days__")),
    median_lag_days = median(wide$protein_sampling_lag_days, na.rm = TRUE),
    min_lag_days = suppressWarnings(min(wide$protein_sampling_lag_days, na.rm = TRUE)),
    max_lag_days = suppressWarnings(max(wide$protein_sampling_lag_days, na.rm = TRUE))
  )
  list(values = wide, coverage = coverage, processing = processing)
}

yur_min_date <- function(dt, columns) {
  columns <- intersect(columns, names(dt))
  if (!length(columns)) return(as.Date(rep(NA_real_, nrow(dt)), origin = "1970-01-01"))
  values <- lapply(columns, function(v) as.numeric(yur_as_date(dt[[v]])))
  out <- do.call(pmin, c(values, list(na.rm = TRUE)))
  out[!is.finite(out)] <- NA_real_
  as.Date(out, origin = "1970-01-01")
}

yur_resolve_fields <- function(names_available, map_file) {
  m <- fread(map_file)
  out <- m[, {
    choices <- strsplit(candidates, ";", fixed = TRUE)[[1]]
    hit <- choices[choices %in% names_available]
    .(resolved = if (length(hit)) hit[[1]] else NA_character_)
  }, by = .(canonical, required, candidates)]
  out[, status := fifelse(!is.na(resolved), "PASS", fifelse(required, "FAIL", "OPTIONAL_MISSING"))]
  out
}

yur_select_outcomes <- function(cfg) {
  x <- fread(cfg$outcomes_file)
  requested <- trimws(strsplit(as.character(cfg$endpoint_subset), ",", fixed = TRUE)[[1]])
  if (length(requested) == 1L && tolower(requested) == "all") return(x)
  missing <- setdiff(requested, x$outcome_id)
  if (length(missing)) stop("Unknown endpoint_subset values: ", paste(missing, collapse = ", "), call. = FALSE)
  x[outcome_id %in% requested]
}

yur_print_help <- function() {
  cat(paste0(
    "Yu/Chen 2025 full-article reproduction\n\n",
    "Modes:\n",
    "  sources    Audit article, supplement workbook/PDF and export Tables S23-S26.\n",
    "  preflight  Validate local phenotype/protein inputs, fields, panel and packages.\n",
    "  cohort     Build the 14-CVD-free incident cohort and fixed 2/3-1/3 split.\n",
    "  cox_prepare  Freeze the retained full panel before parallel Cox shards.\n",
    "  cox_shard    Run one or more endpoint Cox shards (use --endpoint_subset).\n",
    "  cox_merge    Validate and merge all requested Cox endpoint shards.\n",
    "  cox          Sequential compatibility wrapper for prepare/shard/merge.\n",
    "  cmr_prepare  Freeze the baseline-CVD-free CMR cohort and shared protein matrix.\n",
    "  cmr_shard    Run one or more CMR metric association shards.\n",
    "  cmr_merge    Validate and merge all 19 local CMR association shards.\n",
    "  cmr          Sequential compatibility wrapper for local CMR associations.\n",
    "  mr_prepare   Select local Cox hits and freeze local pQTL instruments.\n",
    "  mr_run       Harmonise extracted 13-outcome GWAS rows and run local MR.\n",
    "  mediation_prepare  Build the local five-risk-factor mediation dataset.\n",
    "  mediation_run      Screen paths and run the fast delta-method audit track.\n",
    "  mediation_cmest_pilot  Run one BMI-GDF15-CAD CMAverse timing/QC pilot.\n",
    "  mediation_cmest_shard  Run one resumable CMAverse path shard.\n",
    "  mediation_cmest_merge  Validate and merge all CMAverse shards.\n",
    "  systems_prepare     Build Figure 6B-D inputs from local significant Cox proteins.\n",
    "  systems_enrichment  Run STRING pathway enrichment against the measured panel.\n",
    "  systems_tf          Download TRRUST and build CVD-protein-TF edges.\n",
    "  systems_ppi         Retrieve STRING PPI, cluster it, and calculate MNC hubs.\n",
    "  systems_figures     Draw local Figure 6B, 6C, 6D and the combined panel.\n",
    "  figure6_systems     Run all Figure 6B-D systems stages in sequence.\n",
    "  select     Derivation-only preliminary LightGBM gain selection to cumulative 30%.\n",
    "  train      Fit final SCORE2, Protein and Protein+SCORE2 models.\n",
    "  evaluate   Hold-out metrics, DeLong, NRI/IDI and paired bootstrap.\n",
    "  figures    Generate source-data tables and manuscript figures.\n",
    "  report     Build a source-locked result/QC report.\n",
    "  all_fast   Parallel Cox plus source-locked benchmark.\n",
    "  all        Full workflow; use yu.sh for cross-language dispatch.\n\n",
    "Key options:\n",
    "  --endpoint_subset=all or comma-separated IDs such as cad,heart_failure\n",
    "  --raw_protein_file=<unimputed baseline NPX table>\n",
    "  --raw_phenotype_file=<UKB raw phenotype table containing p74_i0>\n",
    "  --cmr_feature_file=<local UKB cardiac MRI feature table>\n",
    "  --cmr_metric_subset=all or comma-separated metric IDs such as LVEDV,LVEF\n",
    "  --pqtl_root=<local UKB-PPP protein pQTL directories>\n",
    "  --mr_outcome_lookup_dir=<13 extracted full-GWAS lookup files>\n",
    "  --supplement_workbook_file=<official Tables S1-S26 XLSX>\n",
    "  --supplement_methods_file=<official supplementary methods PDF>\n",
    "  --olink_processing_start_date_file=<UKB Resource 1019 table>\n",
    "  --workers=16 --bootstrap_n=1000 --resume=true\n",
    "  --cmest_shard_index=1 --cmest_shard_count=8 --cmest_pilot_nboot=20\n",
    "  --string_required_score=700 --systems_top_n_per_outcome=15 --systems_max_tf=46 --systems_enrichment_fdr=0.05\n",
    "The paper did not release author scripts, participant split EIDs, or tuning code.\n",
    "Required technical covariates are source locked: season from baseline date,\n",
    "fasting time from UKB field 74, and panel-specific sampling-to-processing lag\n",
    "from Olink plate plus UKB Resource 1019.\n\n",
    "This implementation records every published versus inferred parameter.\n"
  ))
}

yur_session_snapshot <- function(cfg) {
  capture.output(sessionInfo(), file = file.path(cfg$paths$logs, "sessionInfo.txt"))
  capture.output(RNGkind(), file = file.path(cfg$paths$logs, "RNGkind.txt"))
}
