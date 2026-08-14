options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "tests/test_cli.R"
project_dir <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
entry <- file.path(project_dir, "f", "cli.R")
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
old_project_dir <- Sys.getenv("PRADEEP_PROJECT_DIR", unset = NA_character_)
Sys.setenv(PRADEEP_PROJECT_DIR = project_dir)
on.exit({
  if (is.na(old_project_dir)) Sys.unsetenv("PRADEEP_PROJECT_DIR") else Sys.setenv(PRADEEP_PROJECT_DIR = old_project_dir)
}, add = TRUE)

help_out <- system2(
  rscript,
  c(shQuote(entry), "--h"),
  stdout = TRUE,
  stderr = TRUE
)
stopifnot(is.null(attr(help_out, "status")))
stopifnot(any(grepl("Pradeep/Schuermans", help_out, fixed = TRUE)))
stopifnot(any(grepl("--step core", help_out, fixed = TRUE)))
stopifnot(any(grepl("--panel 1.5k|3k", help_out, fixed = TRUE)))
stopifnot(any(grepl("--outcome NAME", help_out, fixed = TRUE)))
stopifnot(any(grepl("<OUTPUT_ROOT>/pradeep", help_out, fixed = TRUE)))
stopifnot(!any(grepl("[\u4e00-\u9fff]", help_out, perl = TRUE)))

unknown_out <- suppressWarnings(system2(
  rscript,
  c(shQuote(entry), "--does-not-exist"),
  stdout = TRUE,
  stderr = TRUE
))
stopifnot(identical(as.integer(attr(unknown_out, "status")), 1L))
stopifnot(any(grepl("Unknown option", unknown_out, fixed = TRUE)))

bad_panel_out <- suppressWarnings(system2(
  rscript,
  c(shQuote(entry), "--step", "preflight", "--panel", "invalid"),
  stdout = TRUE,
  stderr = TRUE
))
stopifnot(identical(as.integer(attr(bad_panel_out, "status")), 1L))
stopifnot(any(grepl("--panel must be 1.5k or 3k", bad_panel_out, fixed = TRUE)))

bad_outcome_out <- suppressWarnings(system2(
  rscript,
  c(shQuote(entry), "--step", "preflight", "--outcome", "stroke"),
  stdout = TRUE,
  stderr = TRUE
))
stopifnot(identical(as.integer(attr(bad_outcome_out, "status")), 1L))
stopifnot(any(grepl("--outcome must be cad", bad_outcome_out, fixed = TRUE)))

english_surface_files <- c(
  file.path(project_dir, "README.md"),
  file.path(project_dir, "pradeep.sh"),
  list.files(file.path(project_dir, "f"), pattern = "\\.R$", full.names = TRUE)
)
surface_lines <- unlist(lapply(english_surface_files, readLines, warn = FALSE), use.names = FALSE)
stopifnot(!any(grepl("[\u4e00-\u9fff]", surface_lines, perl = TRUE)))

cat("CLI_TESTS_PASSED\n")
