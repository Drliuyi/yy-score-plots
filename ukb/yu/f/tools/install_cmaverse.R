#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

args <- commandArgs(trailingOnly = TRUE)
manifest <- if (length(args)) args[[1]] else file.path(getwd(), "cmaverse_install_manifest.txt")
source_tarball <- if (length(args) >= 2L) args[[2]] else ""
commit <- "b2ce0598ed362ddfa0db8fff5486ee43d3e73f54"

dir.create(dirname(manifest), recursive = TRUE, showWarnings = FALSE)

imports <- c(
  "ggdag", "stringr", "simex", "dplyr", "mice", "nnet", "MASS", "survey",
  "survival", "SuppDists", "boot", "msm", "Matrix", "EValue", "ggplot2",
  "doSNOW", "foreach", "medflex", "mstate", "predint", "survminer"
)

missing <- imports[!vapply(imports, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing CMAverse imports: ", paste(missing, collapse = ", "))
  install.packages(
    missing,
    dependencies = c("Depends", "Imports", "LinkingTo"),
    Ncpus = max(1L, min(6L, parallel::detectCores(logical = FALSE)))
  )
}

still_missing <- imports[!vapply(imports, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing)) {
  stop("CMAverse dependencies failed to install: ", paste(still_missing, collapse = ", "), call. = FALSE)
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", dependencies = c("Depends", "Imports", "LinkingTo"))
}

installed_commit <- tryCatch(utils::packageDescription("CMAverse")$RemoteSha, error = function(e) NA_character_)
if (!requireNamespace("CMAverse", quietly = TRUE) || !identical(installed_commit, commit)) {
  if (nzchar(source_tarball)) {
    if (!file.exists(source_tarball)) stop("CMAverse source tarball does not exist: ", source_tarball, call. = FALSE)
    install.packages(source_tarball, repos = NULL, type = "source")
  } else {
    remotes::install_github(
      paste0("BS1125/CMAverse@", commit),
      dependencies = FALSE,
      upgrade = "never",
      build_vignettes = FALSE,
      force = TRUE
    )
  }
}

stopifnot(
  requireNamespace("CMAverse", quietly = TRUE),
  identical(as.character(utils::packageVersion("CMAverse")), "0.1.0"),
  exists("cmest", envir = asNamespace("CMAverse"), inherits = FALSE)
)

description <- utils::packageDescription("CMAverse")
dependency_versions <- vapply(
  imports,
  function(package) as.character(utils::packageVersion(package)),
  character(1)
)

lines <- c(
  paste0("installed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  paste0("R_version=", R.version.string),
  paste0("library_paths=", paste(.libPaths(), collapse = ";")),
  paste0("CMAverse_version=", as.character(utils::packageVersion("CMAverse"))),
  paste0(
    "CMAverse_remote_sha=",
    if (is.null(description$RemoteSha) || is.na(description$RemoteSha)) commit else description$RemoteSha
  ),
  paste0("CMAverse_source=https://github.com/BS1125/CMAverse/commit/", commit),
  paste0("CMAverse_source_tarball=", if (nzchar(source_tarball)) normalizePath(source_tarball, winslash = "/") else "GitHub API"),
  paste0("CMAverse_source_md5=", if (nzchar(source_tarball)) unname(tools::md5sum(source_tarball)) else "not_recorded"),
  paste0(names(dependency_versions), "=", unname(dependency_versions))
)
writeLines(lines, manifest, useBytes = TRUE)
cat("CMAVERSE_INSTALL_PASS\n", paste(lines, collapse = "\n"), "\n", sep = "")
