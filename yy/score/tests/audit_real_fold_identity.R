#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
project <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
source(file.path(project, "R", "00_common.R"))
phe_dir <- Sys.getenv("PHEDIR", unset = "/mnt/d/data/ukb/phe")
locked_root <- Sys.getenv(
  "YY_SCORE_FOLD_ROOT", unset = "/mnt/d/analysis/yy/reference/cad_fivefold_v1"
)
config <- fromJSON(file.path(project, "config", "fair.json"), simplifyVector = TRUE)
all_data <- as.data.table(readRDS(file.path(phe_dir, "Rdata", "all.rds")))
protein_data <- as.data.table(readRDS(file.path(phe_dir, "Rdata", "prot.rds")))
generated <- score_generate_fold_tables(all_data, protein_data, config)
locked <- list(
  yin = fread(file.path(locked_root, "fold_assignment_yin.csv"), colClasses = list(character = "eid")),
  yang = fread(file.path(locked_root, "fold_assignment_yang.csv"), colClasses = list(character = "eid"))
)
yin <- merge(
  locked$yin, generated$yin, by = c("eid", "event"), suffixes = c("_locked", "_generated")
)
yang <- merge(
  locked$yang, generated$yang, by = "eid", suffixes = c("_locked", "_generated")
)
cat(
  "REAL_FOLD_AUDIT\n",
  "yin_locked=", nrow(locked$yin), " yin_generated=", nrow(generated$yin),
  " yin_eid_event_overlap=", nrow(yin),
  " yin_fold_match=", sprintf("%.8f", mean(yin$fold_locked == yin$fold_generated)), "\n",
  "yang_locked=", nrow(locked$yang), " yang_generated=", nrow(generated$yang),
  " yang_eid_overlap=", nrow(yang),
  " yang_fold_match=", sprintf("%.8f", mean(yang$fold_locked == yang$fold_generated)), "\n",
  sep = ""
)
