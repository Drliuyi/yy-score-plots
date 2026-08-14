#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else getwd()
manifest <- if (length(args) >= 2L) args[[2]] else file.path(project_dir, "dependency_install_manifest.txt")
dir.create(dirname(manifest), recursive = TRUE, showWarnings = FALSE)

cran_packages <- c(
  "data.table", "R.utils", "jsonlite", "digest", "readxl", "survival", "pROC", "ggplot2",
  "bit64", "patchwork", "ggrepel", "ragg", "httr2", "curl", "igraph",
  "ggalluvial", "scales", "msigdbr", "png", "remotes", "BiocManager"
)

install_missing <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Installing CRAN packages: ", paste(missing, collapse = ", "))
    install.packages(
      missing,
      dependencies = c("Depends", "Imports", "LinkingTo"),
      Ncpus = max(1L, min(6L, parallel::detectCores(logical = FALSE)))
    )
  }
  packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
}

still_missing <- install_missing(cran_packages)
if (length(still_missing)) {
  stop("CRAN packages failed to install: ", paste(still_missing, collapse = ", "), call. = FALSE)
}

bioc_packages <- c("AnnotationDbi", "org.Hs.eg.db", "GO.db")
missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}
still_missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing_bioc)) {
  stop("Bioconductor packages failed to install: ", paste(still_missing_bioc, collapse = ", "), call. = FALSE)
}

if (!requireNamespace("TwoSampleMR", quietly = TRUE) ||
    as.character(utils::packageVersion("TwoSampleMR")) != "0.5.8") {
  remotes::install_github(
    "MRCIEU/TwoSampleMR@v0.5.8",
    dependencies = c("Depends", "Imports", "LinkingTo"),
    upgrade = "never",
    build_vignettes = FALSE,
    force = TRUE
  )
}
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  stop("TwoSampleMR 0.5.8 installation failed.", call. = FALSE)
}

cmaverse_installer <- file.path(project_dir, "f", "tools", "install_cmaverse.R")
cmaverse_tarball <- file.path(project_dir, "references", "raw", "CMAverse_b2ce0598.tar.gz")
cmaverse_manifest <- file.path(dirname(manifest), "cmaverse_install_manifest.txt")
cmaverse_args <- c(cmaverse_installer, cmaverse_manifest)
if (file.exists(cmaverse_tarball)) cmaverse_args <- c(cmaverse_args, cmaverse_tarball)
rscript_bin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
status <- system2(rscript_bin, c("--vanilla", shQuote(cmaverse_args)), stdout = "", stderr = "")
if (!identical(status, 0L) || !requireNamespace("CMAverse", quietly = TRUE)) {
  stop("Pinned CMAverse installation failed.", call. = FALSE)
}

required <- c(cran_packages, bioc_packages, "TwoSampleMR", "CMAverse")
versions <- vapply(required, function(x) as.character(utils::packageVersion(x)), character(1))
lines <- c(
  paste0("installed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  paste0("R_version=", R.version.string),
  paste0("library_paths=", paste(.libPaths(), collapse = ";")),
  paste0(names(versions), "=", unname(versions))
)
writeLines(lines, manifest, useBytes = TRUE)
cat("YU R DEPENDENCIES PASS\n", paste(lines, collapse = "\n"), "\n", sep = "")
