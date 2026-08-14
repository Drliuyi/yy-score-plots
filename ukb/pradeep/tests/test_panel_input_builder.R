options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else "tests/test_panel_input_builder.R"
project_dir <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
builder <- file.path(project_dir, "f", "00_prepare_panel_inputs.R")
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

tmp <- tempfile("pradeep_panel_builder_")
dir.create(tmp, recursive = TRUE)
on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)

removed <- c("CTSS", "NPM1", "PCOLCE", "TACSTD2")
wanted <- sprintf("P%04d", seq_len(1459L))
map <- data.table(Assay = c(wanted, removed))
map_file <- file.path(tmp, "map.tsv")
fwrite(map, map_file, sep = "\t")

prot <- data.table(eid = as.character(1:5))
for (protein in wanted) prot[, (protein) := seq_len(.N)]
set(prot, i = 1L, j = wanted[[1L]], value = NA_real_)
raw_file <- file.path(tmp, "prot.tsv")
fwrite(prot, raw_file, sep = "\t")

all <- data.frame(eid = as.character(1:6), ethnic.c = "White")
for (pc in paste0("PC", 1:10)) all[[pc]] <- seq_len(nrow(all))
all$PC1[[6L]] <- NA_real_
source_all <- file.path(tmp, "all.rds")
saveRDS(all, source_all)

out_dir <- file.path(tmp, "generated")
old_env <- Sys.getenv(c(
  "PRADEEP_PROJECT_DIR", "PRADEEP_SOURCE_ALL_RDS", "PRADEEP_RAW_PROTEIN",
  "PRADEEP_PROTEIN_MAP", "PRADEEP_PANEL_INPUT_DIR"
), unset = NA_character_)
on.exit({
  for (key in names(old_env)) {
    if (is.na(old_env[[key]])) Sys.unsetenv(key) else do.call(Sys.setenv, setNames(list(old_env[[key]]), key))
  }
}, add = TRUE)

Sys.setenv(
  PRADEEP_PROJECT_DIR = project_dir,
  PRADEEP_SOURCE_ALL_RDS = source_all,
  PRADEEP_RAW_PROTEIN = raw_file,
  PRADEEP_PROTEIN_MAP = map_file,
  PRADEEP_PANEL_INPUT_DIR = out_dir
)

status <- system2(rscript, c("--vanilla", shQuote(builder)))
stopifnot(identical(status, 0L))
generated_prot <- readRDS(file.path(out_dir, "prot.rds"))
generated_all <- readRDS(file.path(out_dir, "all.rds"))
stopifnot(nrow(generated_prot) == 5L, ncol(generated_prot) == 1460L)
stopifnot(sum(is.na(generated_prot)) == 1L)
stopifnot(nrow(generated_all) == 5L)
stopifnot(nrow(fread(file.path(out_dir, "strict_1459_assay_mapping.csv"))) == 1459L)
reuse_status <- system2(rscript, c("--vanilla", shQuote(builder)))
stopifnot(identical(reuse_status, 0L))

cat("PANEL_INPUT_BUILDER_TESTS_PASSED\n")
