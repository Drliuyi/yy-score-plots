#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

pkg <- data.table(package = c("present", "missing"), available = c(TRUE, FALSE))
inputs <- data.table(path = c("present", "missing"), exists = c(TRUE, FALSE))

stopifnot(identical(pkg[available == FALSE, package], "missing"))
stopifnot(identical(inputs[exists == FALSE, path], "missing"))

cat("test_preflight_data_table_scope.R: PASS\n")
