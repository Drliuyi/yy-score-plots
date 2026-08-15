#!/usr/bin/env Rscript

project <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
config <- jsonlite::fromJSON(file.path(project, "config", "fair.json"), simplifyVector = TRUE)
source(file.path(project, "R", "00_common.R"))
stopifnot(
  config$expected$yin_n == 37127L,
  config$expected$yin_events == 3442L,
  config$expected$yang_n == 1766L,
  config$expected$protein_n == 2910L,
  config$expected$outer_folds == 5L,
  config$fold_generation$yin_seed == 2026L,
  config$fold_generation$yang_seed == 2027L,
  config$pradeep_fair$lambda_rule == "lambda.1se",
  config$yu_fair$device_type == "gpu",
  config$restricted_reference$yin_md5 == "befca4f9c0c1c2b90afd00e1c30a812e",
  config$restricted_reference$yang_md5 == "14f882efef58ec6754b586c0f4307ce1",
  config$restricted_reference$yin_canonical_md5 == "d2d6aae96c662bfdf0756ac0611a40a2",
  config$restricted_reference$yang_canonical_md5 == "e56061bb00588ff936a8525309d1fc34"
)

# Participant-level fold manifests are generated runtime outputs and remain
# absent from the public source tree. If a generated root is exposed, validate
# it without printing any EID.
fold_root <- Sys.getenv("YY_SCORE_FOLD_ROOT", unset = "")
if (nzchar(fold_root)) {
  yin_file <- file.path(fold_root, config$restricted_reference$yin_file)
  yang_file <- file.path(fold_root, config$restricted_reference$yang_file)
  present <- file.exists(c(yin_file, yang_file))
  stopifnot(!xor(present[[1L]], present[[2L]]))
  if (all(present)) {
    stopifnot(
      score_canonical_fold_md5(yin_file, "yin") == config$restricted_reference$yin_canonical_md5,
      score_canonical_fold_md5(yang_file, "yang") == config$restricted_reference$yang_canonical_md5
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
}
cat("test_contracts.R: PASS\n")
