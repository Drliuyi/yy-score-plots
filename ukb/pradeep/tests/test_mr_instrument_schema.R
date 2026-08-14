options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "tests/test_mr_instrument_schema.R"
project_dir <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
mr_file <- file.path(project_dir, "f", "08_local_pqtl_mr_coloc.R")

exprs <- parse(mr_file)
target <- NULL
for (expr in exprs) {
  if (
    is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("run_mr_for_outcome"))
  ) {
    target <- expr[[3]]
    break
  }
}

stopifnot(!is.null(target))
body_text <- paste(deparse(target, width.cutoff = 500L), collapse = "\n")
stopifnot(grepl("exp = protein", body_text, fixed = TRUE))
stopifnot(grepl("outc = info$outc_lower", body_text, fixed = TRUE))
stopifnot(grepl("pvalthreshold = mr_p_threshold", body_text, fixed = TRUE))
stopifnot(grepl("rsqthreshold = mr_r2", body_text, fixed = TRUE))

cat("MR_INSTRUMENT_SCHEMA_TEST_PASSED\n")
