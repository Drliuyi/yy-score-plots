`%||%` <- function(x, y) if (is.null(x) || !length(x) || !nzchar(as.character(x[[1L]]))) y else x

score_env <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

score_norm <- function(path, must_work = FALSE) {
  value <- gsub("\\\\", "/", path.expand(as.character(path[[1L]])))
  if (!grepl("^(/|[A-Za-z]:/)", value)) value <- file.path(getwd(), value)
  value <- if (grepl("^[A-Za-z]:/$", value)) value else sub("/+$", "", value)
  if (isTRUE(must_work) && !file.exists(value) && !dir.exists(value)) {
    stop("Required path does not exist: ", value, call. = FALSE)
  }
  value
}

score_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  result <- list(stage = "help", fold = NA_integer_, workers = 1L)
  for (arg in args) {
    if (arg %in% c("--h", "-h", "--help")) {
      result$stage <- "help"
    } else if (grepl("^--[^=]+=", arg)) {
      key <- gsub("-", "_", sub("^--([^=]+)=.*$", "\\1", arg))
      value <- sub("^--[^=]+=", "", arg)
      if (!key %in% c(
        "stage", "fold", "workers", "dir0", "phe_dir", "script_root",
        "analysis_root", "yy_outdir", "output_root", "common_root", "config",
        "method", "source_root", "fold_root"
      )) stop("Unknown argument: ", arg, call. = FALSE)
      result[[key]] <- value
    } else {
      stop("Arguments must use --name=value: ", arg, call. = FALSE)
    }
  }
  result$fold <- suppressWarnings(as.integer(result$fold))
  result$workers <- suppressWarnings(as.integer(result$workers))
  result
}

score_resolve <- function(project_dir, parsed, method = "") {
  default_dir0 <- if (.Platform$OS.type == "windows") "D:/" else "/mnt/d"
  dir0 <- parsed$dir0 %||% score_env("DIR0", default_dir0)
  phe_dir <- parsed$phe_dir %||% score_env("PHEDIR", file.path(dir0, "data", "ukb", "phe"))
  script_root <- parsed$script_root %||% score_env("SCRIPT_ROOT", file.path(dir0, "scripts"))
  analysis_root <- parsed$analysis_root %||% score_env("ANALYSIS_ROOT", file.path(dir0, "analysis"))
  yy_outdir <- parsed$yy_outdir %||% score_env("YY_OUTDIR", file.path(analysis_root, "yy"))
  config <- parsed$config %||% file.path(project_dir, "config", "fair.json")
  common_root <- parsed$common_root %||% score_env(
    "YY_SCORE_COMMON_ROOT", file.path(yy_outdir, "score", "common-fair-inputs", "cad")
  )
  fold_root <- parsed$fold_root %||% score_env(
    "YY_SCORE_FOLD_ROOT", file.path(yy_outdir, "reference", "cad_fivefold_v1")
  )
  output_root <- parsed$output_root %||% if (nzchar(method)) {
    score_env("YY_SCORE_OUTPUT_ROOT", file.path(yy_outdir, "score", method))
  } else common_root
  list(
    project_dir = score_norm(project_dir, TRUE),
    dir0 = score_norm(dir0), phe_dir = score_norm(phe_dir),
    script_root = score_norm(script_root), analysis_root = score_norm(analysis_root),
    yy_outdir = score_norm(yy_outdir), config = score_norm(config, TRUE),
    common_root = score_norm(common_root), output_root = score_norm(output_root),
    fold_root = score_norm(fold_root),
    all_rds = score_norm(file.path(phe_dir, "Rdata", "all.rds")),
    prot_rds = score_norm(file.path(phe_dir, "Rdata", "prot.rds")),
    fold_yin = score_norm(file.path(fold_root, "fold_assignment_yin.csv")),
    fold_yang = score_norm(file.path(fold_root, "fold_assignment_yang.csv"))
  )
}

score_print_roots <- function(paths) {
  cat(
    "Resolved Huang-lab roots:\n",
    "  DIR0=", paths$dir0, "\n",
    "  PHEDIR=", paths$phe_dir, "\n",
    "  SCRIPT_ROOT=", paths$script_root, "\n",
    "  ANALYSIS_ROOT=", paths$analysis_root, "\n",
    "  YY_OUTDIR=", paths$yy_outdir, "\n",
    "Resolved inputs/outputs:\n",
    "  ALL_RDS=", paths$all_rds, "\n",
    "  PROT_RDS=", paths$prot_rds, "\n",
    "  FOLD_ROOT=", paths$fold_root, "\n",
    "  COMMON_ROOT=", paths$common_root, "\n",
    "  OUTPUT_ROOT=", paths$output_root, "\n",
    sep = ""
  )
}

score_require <- function(paths, label = "required input") {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop(label, " missing: ", paste(missing, collapse = "; "), call. = FALSE)
}

score_as_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1970-01-01"))
  suppressWarnings(as.Date(as.character(x)))
}

score_row_min_date <- function(...) {
  value <- do.call(cbind, lapply(list(...), function(x) as.numeric(score_as_date(x))))
  result <- apply(value, 1L, function(row) if (all(!is.finite(row))) NA_real_ else min(row, na.rm = TRUE))
  as.Date(result, origin = "1970-01-01")
}

score_build_endpoint <- function(data, follow_end) {
  required <- c("eid", "date_attend", "fod_icd10_cvd_cad", "date_lost", "date_death")
  score_require_columns <- setdiff(required, names(data))
  if (length(score_require_columns)) {
    stop("all.rds is missing endpoint columns: ", paste(score_require_columns, collapse = ", "), call. = FALSE)
  }
  baseline <- score_as_date(data$date_attend)
  event_date <- score_as_date(data$fod_icd10_cvd_cad)
  censor <- score_row_min_date(data$date_lost, data$date_death, rep(as.Date(follow_end), nrow(data)))
  censor[is.na(censor)] <- as.Date(follow_end)
  prevalent <- !is.na(event_date) & !is.na(baseline) & event_date < baseline
  incident <- !is.na(event_date) & !is.na(baseline) & event_date >= baseline & event_date <= censor
  event_time <- as.numeric(event_date - baseline) / 365.25
  follow_time <- as.numeric(censor - baseline) / 365.25
  time <- ifelse(incident, event_time, follow_time)
  time[!is.finite(time) | time <= 0] <- NA_real_
  data.table::data.table(
    eid = as.character(data$eid), time = time, event = as.integer(incident),
    prevalent = as.integer(prevalent), b2e = event_time
  )
}

score_stratified_folds <- function(stratum, folds, seed) {
  stratum <- as.integer(stratum)
  folds <- as.integer(folds)
  seed <- as.integer(seed)
  if (!length(stratum) || anyNA(stratum) || !is.finite(folds) || folds < 2L ||
      !is.finite(seed)) {
    stop("Invalid stratified-fold inputs.", call. = FALSE)
  }
  set.seed(seed)
  result <- integer(length(stratum))
  for (value in unique(stratum)) {
    index <- which(stratum == value)
    result[index] <- sample(rep(seq_len(folds), length.out = length(index)))
  }
  result
}

score_fold_generation_settings <- function(config) {
  generation <- config$fold_generation
  if (is.null(generation)) stop("fair.json lacks fold_generation settings.", call. = FALSE)
  breaks <- unlist(generation$yang_duration_breaks, use.names = FALSE)
  breaks <- vapply(breaks, function(value) {
    if (identical(as.character(value), "Inf")) Inf else as.numeric(value)
  }, numeric(1L))
  if (length(breaks) < 2L || anyNA(breaks) || is.unsorted(breaks, strictly = TRUE)) {
    stop("Invalid Yang duration breaks in fair.json.", call. = FALSE)
  }
  list(
    ethnicity_column = as.character(generation$ethnicity_column),
    ethnicity_pattern = as.character(generation$ethnicity_pattern),
    include_missing_ethnicity = isTRUE(generation$include_missing_ethnicity),
    yin_seed = as.integer(generation$yin_seed),
    yang_seed = as.integer(generation$yang_seed),
    yang_duration_breaks = breaks,
    outer_folds = as.integer(config$expected$outer_folds)
  )
}

score_generate_fold_tables <- function(all_data, protein_data, config) {
  settings <- score_fold_generation_settings(config)
  all_data <- data.table::as.data.table(all_data)
  protein_data <- data.table::as.data.table(protein_data)
  if (!"eid" %in% names(all_data) || !"eid" %in% names(protein_data)) {
    stop("all.rds and prot.rds must contain eid before fold generation.", call. = FALSE)
  }
  if (!settings$ethnicity_column %in% names(all_data)) {
    stop("all.rds lacks locked ethnicity column: ", settings$ethnicity_column, call. = FALSE)
  }
  all_data[, eid := as.character(eid)]
  protein_ids <- unique(as.character(protein_data$eid))
  eligible <- all_data[eid %chin% protein_ids]
  ethnicity <- as.character(eligible[[settings$ethnicity_column]])
  keep_ethnicity <- grepl(settings$ethnicity_pattern, ethnicity, ignore.case = TRUE)
  if (settings$include_missing_ethnicity) keep_ethnicity <- is.na(ethnicity) | keep_ethnicity
  eligible <- eligible[keep_ethnicity]
  data.table::setorder(eligible, eid)

  endpoint <- score_build_endpoint(eligible, config$follow_end)
  yin <- endpoint[prevalent == 0L & is.finite(time) & time > 0,
                  .(eid, event = as.integer(event))]
  yang <- endpoint[prevalent == 1L & is.finite(b2e) & b2e < 0,
                   .(eid, disease_duration_years = abs(as.numeric(b2e)))]
  data.table::setorder(yin, eid)
  data.table::setorder(yang, eid)
  yin[, fold := score_stratified_folds(event, settings$outer_folds, settings$yin_seed)]
  duration_group <- cut(
    yang$disease_duration_years,
    breaks = settings$yang_duration_breaks,
    include.lowest = TRUE
  )
  if (anyNA(duration_group)) stop("Yang duration stratification produced missing groups.", call. = FALSE)
  yang[, fold := score_stratified_folds(
    as.integer(duration_group), settings$outer_folds, settings$yang_seed
  )]
  list(
    yin = yin[, .(eid, event, fold)],
    yang = yang[, .(eid, fold)],
    settings = settings
  )
}

score_validate_fold_tables <- function(folds, config) {
  yin <- data.table::as.data.table(folds$yin)
  yang <- data.table::as.data.table(folds$yang)
  required_yin <- c("eid", "event", "fold")
  required_yang <- c("eid", "fold")
  if (length(setdiff(required_yin, names(yin))) ||
      length(setdiff(required_yang, names(yang)))) {
    stop("Fold manifest columns are invalid.", call. = FALSE)
  }
  expected <- config$expected
  checks <- c(
    nrow(yin) == as.integer(expected$yin_n),
    sum(as.integer(yin$event)) == as.integer(expected$yin_events),
    nrow(yang) == as.integer(expected$yang_n),
    !anyDuplicated(as.character(yin$eid)),
    !anyDuplicated(as.character(yang$eid)),
    !length(intersect(as.character(yin$eid), as.character(yang$eid))),
    identical(sort(unique(as.integer(yin$fold))), seq_len(as.integer(expected$outer_folds))),
    identical(sort(unique(as.integer(yang$fold))), seq_len(as.integer(expected$outer_folds)))
  )
  if (!all(checks)) {
    stop(
      "Generated fold contract failed: yin=", nrow(yin),
      " events=", sum(as.integer(yin$event)), " yang=", nrow(yang),
      ". Check the all.rds/prot.rds versions and CAD endpoint fields.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

score_expected_fold_md5 <- function(config) {
  c(
    yin = as.character(config$restricted_reference$yin_canonical_md5),
    yang = as.character(config$restricted_reference$yang_canonical_md5)
  )
}

score_canonical_fold_md5 <- function(path, side) {
  columns <- if (identical(side, "yin")) c("eid", "event", "fold") else c("eid", "fold")
  value <- data.table::fread(path, colClasses = list(character = "eid"))
  if (length(setdiff(columns, names(value)))) {
    stop("Fold file lacks canonical columns: ", path, call. = FALSE)
  }
  value <- value[, ..columns]
  data.table::setorder(value, eid)
  temporary <- tempfile("canonical_fold_", fileext = ".csv")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  lines <- c(
    paste(columns, collapse = ","),
    do.call(paste, c(value, sep = ","))
  )
  writeLines(lines, temporary, useBytes = TRUE)
  unname(tools::md5sum(temporary))
}

score_write_fold_candidates <- function(folds, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  yin <- tempfile("fold_assignment_yin.", tmpdir = directory, fileext = ".csv")
  yang <- tempfile("fold_assignment_yang.", tmpdir = directory, fileext = ".csv")
  data.table::fwrite(folds$yin, yin)
  data.table::fwrite(folds$yang, yang)
  list(yin = yin, yang = yang)
}

score_validate_fold_file_identity <- function(files, config) {
  observed <- c(
    score_canonical_fold_md5(files$yin, "yin"),
    score_canonical_fold_md5(files$yang, "yang")
  )
  expected <- unname(score_expected_fold_md5(config))
  if (!identical(observed, expected)) {
    stop(
      "Canonical fold content differs from the frozen CAD reference. Expected ",
      paste(expected, collapse = ", "), "; generated ", paste(observed, collapse = ", "),
      ". No fold files were installed. Check source versions and cohort rules.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

score_ensure_fold_manifests <- function(paths, config, all_data, protein_data, install = FALSE) {
  existing <- file.exists(c(paths$fold_yin, paths$fold_yang))
  if (xor(existing[[1L]], existing[[2L]])) {
    stop(
      "Partial fold reference detected; refusing regeneration. Present: ",
      paste(c(paths$fold_yin, paths$fold_yang)[existing], collapse = "; "),
      call. = FALSE
    )
  }
  if (all(existing)) {
    files <- list(yin = paths$fold_yin, yang = paths$fold_yang)
    score_validate_fold_file_identity(files, config)
    folds <- list(
      yin = data.table::fread(paths$fold_yin, colClasses = list(character = "eid")),
      yang = data.table::fread(paths$fold_yang, colClasses = list(character = "eid"))
    )
    score_validate_fold_tables(folds, config)
    cat("FOLD_MANIFESTS_VALIDATED existing\n")
    return(invisible(folds))
  }

  folds <- score_generate_fold_tables(all_data, protein_data, config)
  score_validate_fold_tables(folds, config)
  temporary_directory <- tempfile("cad_fivefold_generation_")
  dir.create(temporary_directory, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(temporary_directory, recursive = TRUE, force = TRUE), add = TRUE)
  temporary <- score_write_fold_candidates(folds, temporary_directory)
  score_validate_fold_file_identity(temporary, config)
  if (!isTRUE(install)) {
    cat("FOLD_MANIFESTS_GENERATABLE exact_frozen_identity\n")
    return(invisible(folds))
  }

  dir.create(paths$fold_root, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(paths$fold_yin) || file.exists(paths$fold_yang)) {
    stop("Fold files appeared during generation; refusing to overwrite.", call. = FALSE)
  }
  if (!file.rename(temporary$yin, paths$fold_yin)) {
    stop("Cannot install generated Yin fold manifest: ", paths$fold_yin, call. = FALSE)
  }
  if (!file.rename(temporary$yang, paths$fold_yang)) {
    unlink(paths$fold_yin, force = TRUE)
    stop("Cannot install generated Yang fold manifest: ", paths$fold_yang, call. = FALSE)
  }
  manifest <- data.table::data.table(
    item = c("status", "source_all_rds", "source_prot_rds", "yin_seed", "yang_seed",
            "yin_file_md5", "yang_file_md5", "yin_canonical_md5",
            "yang_canonical_md5", "generated_at"),
    value = c(
      "GENERATED_FROM_HUANG_RAW_INPUTS",
      normalizePath(paths$all_rds, winslash = "/", mustWork = TRUE),
      normalizePath(paths$prot_rds, winslash = "/", mustWork = TRUE),
      as.character(folds$settings$yin_seed), as.character(folds$settings$yang_seed),
      unname(tools::md5sum(paths$fold_yin)), unname(tools::md5sum(paths$fold_yang)),
      score_canonical_fold_md5(paths$fold_yin, "yin"),
      score_canonical_fold_md5(paths$fold_yang, "yang"),
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
    )
  )
  score_atomic_csv(manifest, file.path(paths$fold_root, "fold_generation_manifest.csv"))
  cat("FOLD_MANIFESTS_GENERATED ", paths$fold_root, "\n", sep = "")
  invisible(folds)
}

score_known_horizon <- function(time, event, horizon = 5) {
  known <- (event == 1L & time <= horizon) | time > horizon
  list(known = known, label = as.integer(event == 1L & time <= horizon))
}

score_write_f32 <- function(matrix_value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  connection <- file(temporary, open = "wb")
  on.exit({
    try(close(connection), silent = TRUE)
    if (file.exists(temporary)) unlink(temporary)
  }, add = TRUE)
  writeBin(as.numeric(matrix_value), connection, size = 4L, endian = "little")
  close(connection)
  expected <- as.double(nrow(matrix_value)) * as.double(ncol(matrix_value)) * 4
  if (file.info(temporary)$size != expected) stop("Float32 write size failed: ", path, call. = FALSE)
  if (!file.rename(temporary, path)) stop("Atomic float32 rename failed: ", path, call. = FALSE)
  invisible(path)
}

score_read_f32 <- function(path, rows, columns) {
  expected <- as.double(rows) * as.double(columns) * 4
  if (!file.exists(path) || file.info(path)$size != expected) {
    stop("Float32 matrix contract failed: ", path, call. = FALSE)
  }
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  value <- readBin(connection, what = "numeric", n = rows * columns, size = 4L, endian = "little")
  if (length(value) != rows * columns) stop("Float32 read incomplete: ", path, call. = FALSE)
  matrix(value, nrow = rows, ncol = columns)
}

score_atomic_csv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- if (grepl("[.]gz$", path)) {
    paste0(sub("[.]gz$", "", path), ".tmp.", Sys.getpid(), ".gz")
  } else paste0(path, ".tmp.", Sys.getpid())
  data.table::fwrite(data, temporary)
  if (!file.rename(temporary, path)) stop("Atomic CSV rename failed: ", path, call. = FALSE)
  invisible(path)
}

score_atomic_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  saveRDS(object, temporary, compress = "xz")
  if (!file.rename(temporary, path)) stop("Atomic RDS rename failed: ", path, call. = FALSE)
  invisible(path)
}

score_atomic_text <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  writeLines(as.character(lines), temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) stop("Atomic text rename failed: ", path, call. = FALSE)
  invisible(path)
}

score_inner_foldid <- function(label, folds, seed) {
  set.seed(as.integer(seed))
  result <- integer(length(label))
  for (value in sort(unique(label))) {
    index <- sample(which(label == value))
    result[index] <- rep(seq_len(folds), length.out = length(index))
  }
  result
}

score_fit_preprocessor <- function(x) {
  center <- apply(x, 2L, function(value) {
    finite <- value[is.finite(value)]
    if (length(finite)) stats::median(finite) else 0
  })
  for (j in seq_len(ncol(x))) x[!is.finite(x[, j]), j] <- center[[j]]
  mean_value <- colMeans(x)
  x <- sweep(x, 2L, mean_value, "-")
  scale_value <- sqrt(colSums(x * x) / pmax(1, nrow(x) - 1L))
  scale_value[!is.finite(scale_value) | scale_value < 1e-8] <- 1
  list(x = sweep(x, 2L, scale_value, "/"), median = center, center = mean_value, scale = scale_value)
}

score_apply_preprocessor <- function(x, spec) {
  for (j in seq_len(ncol(x))) x[!is.finite(x[, j]), j] <- spec$median[[j]]
  sweep(sweep(x, 2L, spec$center, "-"), 2L, spec$scale, "/")
}
