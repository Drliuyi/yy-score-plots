#!/usr/bin/env Rscript

project <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
config <- jsonlite::fromJSON(file.path(project, "config", "fair.json"), simplifyVector = TRUE)
stopifnot(
  config$expected$yin_n == 37127L,
  config$expected$yin_events == 3442L,
  config$expected$yang_n == 1766L,
  config$expected$protein_n == 2910L,
  config$expected$outer_folds == 5L,
  config$pradeep_fair$lambda_rule == "lambda.1se",
  config$yu_fair$device_type == "gpu",
  config$restricted_reference$yin_md5 == "befca4f9c0c1c2b90afd00e1c30a812e",
  config$restricted_reference$yang_md5 == "14f882efef58ec6754b586c0f4307ce1"
)

# Participant-level fold manifests are restricted UKB-derived inputs and are
# intentionally absent from the public source tree. If a reviewed protected
# root is exposed, validate it without printing any EID.
fold_root <- Sys.getenv("YY_SCORE_FOLD_ROOT", unset = "")
if (nzchar(fold_root)) {
  yin_file <- file.path(fold_root, config$restricted_reference$yin_file)
  yang_file <- file.path(fold_root, config$restricted_reference$yang_file)
  stopifnot(file.exists(yin_file), file.exists(yang_file))
  stopifnot(
    unname(tools::md5sum(yin_file)) == config$restricted_reference$yin_md5,
    unname(tools::md5sum(yang_file)) == config$restricted_reference$yang_md5
  )
  yin <- data.table::fread(yin_file, colClasses = list(character = "eid"))
  yang <- data.table::fread(yang_file, colClasses = list(character = "eid"))
  stopifnot(
    nrow(yin) == 37127L,
    sum(yin$event) == 3442L,
    nrow(yang) == 1766L,
    !anyDuplicated(yin$eid),
    !anyDuplicated(yang$eid),
    !length(intersect(yin$eid, yang$eid)),
    identical(sort(unique(yin$fold)), 1:5),
    identical(sort(unique(yang$fold)), 1:5)
  )
}
cat("test_contracts.R: PASS\n")
