#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(glmnet)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
project_dir <- dirname(dirname(gsub("\\\\", "/", script_file)))
source(file.path(project_dir, "R", "00_common.R"))
parsed <- score_args()
paths <- score_resolve(project_dir, parsed, "pradeep-fair")
cfg <- jsonlite::fromJSON(paths$config, simplifyVector = TRUE)
score_print_roots(paths)

if (parsed$stage == "help") {
  cat("Usage: Rscript 02_pradeep_fair.R --stage=fold|status --fold=1..5 [root overrides]\n")
  quit(status = 0L)
}
if (parsed$stage == "status") {
  markers <- file.path(paths$output_root, "models", sprintf("fold%02d_COMPLETE", 1:5))
  cat(sprintf("fold%02d\t%s\n", 1:5, ifelse(file.exists(markers), "COMPLETE", "MISSING")))
  quit(status = if (all(file.exists(markers))) 0L else 1L)
}
if (parsed$stage != "fold" || !parsed$fold %in% 1:5) {
  stop("Use --stage=fold --fold=1..5.", call. = FALSE)
}
fold <- parsed$fold
common <- paths$common_root
required <- file.path(common, c(
  "COMPLETE", "participants_yin.csv.gz", "participants_yang.csv.gz",
  "protein_features.csv", "protein_yin.f32", "protein_yang.f32"
))
score_require(required, "fair common input")

prediction_path <- file.path(paths$output_root, "predictions", sprintf("fold%02d_yin.csv.gz", fold))
yang_path <- file.path(paths$output_root, "predictions", sprintf("fold%02d_yang.csv.gz", fold))
model_path <- file.path(paths$output_root, "models", sprintf("fold%02d_model.rds", fold))
marker <- file.path(paths$output_root, "models", sprintf("fold%02d_COMPLETE", fold))
outputs <- c(prediction_path, yang_path, model_path, marker)
if (all(file.exists(outputs))) {
  cat("PRADEEP_FAIR_FOLD_PRESERVED fold=", fold, "\n", sep = "")
  quit(status = 0L)
}
if (any(file.exists(outputs))) {
  stop("Partial Pradeep fair fold exists; refusing overwrite: fold ", fold, call. = FALSE)
}

yin <- fread(file.path(common, "participants_yin.csv.gz"), colClasses = list(character = "eid"))
yang <- fread(file.path(common, "participants_yang.csv.gz"), colClasses = list(character = "eid"))
features <- fread(file.path(common, "protein_features.csv"))$feature
if (nrow(yin) != cfg$expected$yin_n || sum(yin$event) != cfg$expected$yin_events ||
    nrow(yang) != cfg$expected$yang_n || length(features) != cfg$expected$protein_n) {
  stop("Pradeep fair input count contract failed.", call. = FALSE)
}
x_all <- score_read_f32(file.path(common, "protein_yin.f32"), nrow(yin), length(features))
x_yang_raw <- score_read_f32(file.path(common, "protein_yang.f32"), nrow(yang), length(features))
colnames(x_all) <- colnames(x_yang_raw) <- features
train_index <- which(yin$outer_fold != fold)
test_index <- which(yin$outer_fold == fold)
known <- score_known_horizon(yin$time[train_index], yin$event[train_index], cfg$horizon_years)
known_train_index <- train_index[known$known]
y_train <- known$label[known$known]
if (length(unique(y_train)) != 2L) stop("Pradeep fair training label has fewer than two classes.", call. = FALSE)

preprocess <- score_fit_preprocessor(x_all[known_train_index, , drop = FALSE])
x_train <- preprocess$x
x_test <- score_apply_preprocessor(x_all[test_index, , drop = FALSE], preprocess)
x_yang <- score_apply_preprocessor(x_yang_raw, preprocess)
inner_foldid <- score_inner_foldid(
  y_train, as.integer(cfg$pradeep_fair$inner_folds),
  as.integer(cfg$seed) + fold * 1000L
)
cat("PRADEEP_FAIR_FIT fold=", fold, " train_known_n=", length(y_train),
    " cases=", sum(y_train), " test_n=", length(test_index), "\n", sep = "")
fit <- cv.glmnet(
  x = x_train, y = y_train, family = "binomial",
  alpha = as.numeric(cfg$pradeep_fair$alpha),
  foldid = inner_foldid, type.measure = cfg$pradeep_fair$type_measure,
  standardize = FALSE, intercept = TRUE, parallel = FALSE,
  maxit = 1000000L
)
lambda <- unname(fit[[cfg$pradeep_fair$lambda_rule]])
yin_score <- as.numeric(predict(fit, newx = x_test, s = lambda, type = "link"))
yang_score <- as.numeric(predict(fit, newx = x_yang, s = lambda, type = "link"))
coefficient <- as.matrix(coef(fit, s = lambda))[, 1L]
if (any(!is.finite(yin_score)) || any(!is.finite(yang_score)) || any(!is.finite(coefficient))) {
  stop("Pradeep fair produced non-finite output.", call. = FALSE)
}
yin_output <- yin[test_index, .(eid, outer_fold, time, event)]
yin_output[, `:=`(method = "pradeep-fair", score_raw = yin_score)]
yang_output <- data.table(
  eid = yang$eid, outer_fold = fold, method = "pradeep-fair", score_raw = yang_score
)
model <- list(
  method = "pradeep-fair", outer_fold = fold,
  target = "five-year incident CAD", family = "binomial LASSO",
  lambda = lambda, lambda_min = unname(fit$lambda.min), lambda_1se = unname(fit$lambda.1se),
  coefficient = coefficient, feature_order = features,
  preprocessing = preprocess[c("median", "center", "scale")],
  training_n = length(y_train), training_cases = sum(y_train),
  seed = as.integer(cfg$seed) + fold * 1000L
)
score_atomic_csv(yin_output, prediction_path)
score_atomic_csv(yang_output, yang_path)
score_atomic_rds(model, model_path)
score_atomic_text(c(
  "status=COMPLETE", paste0("fold=", fold), paste0("lambda=", format(lambda, digits = 16)),
  paste0("selected_proteins=", sum(coefficient[names(coefficient) != "(Intercept)"] != 0)),
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
), marker)
cat("PRADEEP_FAIR_FOLD_COMPLETE fold=", fold, " selected=",
    sum(coefficient[names(coefficient) != "(Intercept)"] != 0), "\n", sep = "")
