env_or <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

protein_key <- function(x) gsub("[^A-Z0-9]", "", toupper(as.character(x)))

atomic_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- if (grepl("\\.gz$", path)) {
    paste0(sub("\\.gz$", "", path), ".tmp.", Sys.getpid(), ".gz")
  } else {
    paste0(path, ".tmp.", Sys.getpid())
  }
  data.table::fwrite(x, temporary)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Atomic CSV rename failed: ", path, call. = FALSE)
  }
  invisible(path)
}

atomic_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  saveRDS(x, temporary, compress = "xz")
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Atomic RDS rename failed: ", path, call. = FALSE)
  }
  invisible(path)
}

read_f32_fortran <- function(path, rows, columns) {
  expected <- as.double(rows) * as.double(columns) * 4
  if (!file.exists(path) || file.info(path)$size != expected) {
    stop("Float32 matrix contract failed: ", path, call. = FALSE)
  }
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  value <- readBin(connection, what = "numeric", n = rows * columns,
                   size = 4L, endian = "little")
  if (length(value) != rows * columns) stop("Incomplete float32 read.", call. = FALSE)
  matrix(value, nrow = rows, ncol = columns)
}

read_f32_fortran_columns <- function(path, rows, columns, column_index) {
  expected <- as.double(rows) * as.double(columns) * 4
  if (!file.exists(path) || file.info(path)$size != expected) {
    stop("Float32 matrix contract failed: ", path, call. = FALSE)
  }
  column_index <- as.integer(column_index)
  if (!length(column_index) || anyNA(column_index) ||
      any(column_index < 1L | column_index > columns) || anyDuplicated(column_index)) {
    stop("Selected float32 column contract failed.", call. = FALSE)
  }
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  value <- matrix(NA_real_, nrow = rows, ncol = length(column_index))
  for (j in seq_along(column_index)) {
    seek(connection, where = as.double(column_index[[j]] - 1L) * rows * 4,
         origin = "start", rw = "read")
    one <- readBin(connection, what = "numeric", n = rows, size = 4L,
                   endian = "little")
    if (length(one) != rows) stop("Incomplete selected float32 read.", call. = FALSE)
    value[, j] <- one
  }
  value
}

age_sex_covariates <- function(age, sex) {
  age <- suppressWarnings(as.numeric(age))
  sex <- suppressWarnings(as.integer(sex))
  if (any(!is.finite(age))) stop("Age must be complete and finite.", call. = FALSE)
  if (any(!sex %in% c(0L, 1L))) stop("Sex must be complete and coded 0/1.", call. = FALSE)
  data.table::data.table(
    age = age,
    sex_factor = factor(sex, levels = c(0L, 1L), labels = c("Female", "Male"))
  )
}

age_sex_adjusted_bins <- function(dat, label, series_type) {
  required <- c("value", "relative_bin", "side", "age", "sex")
  if (!all(required %in% names(dat))) {
    stop("Age/sex-adjusted bin input is missing: ",
         paste(setdiff(required, names(dat)), collapse = ", "), call. = FALSE)
  }
  covariates <- age_sex_covariates(dat$age, dat$sex)
  d <- data.table::data.table(
    value = suppressWarnings(as.numeric(dat$value)),
    relative_bin = suppressWarnings(as.numeric(dat$relative_bin)),
    side = as.character(dat$side), age = covariates$age,
    sex_factor = covariates$sex_factor
  )
  if (any(!is.finite(d$value)) || any(!is.finite(d$relative_bin))) {
    stop("Age/sex-adjusted bin input contains non-finite values.", call. = FALSE)
  }
  bin_levels <- sort(unique(d$relative_bin))
  if (length(bin_levels) < 2L) stop("At least two bins are required.", call. = FALSE)
  d[, bin_factor := factor(relative_bin, levels = bin_levels)]
  fit <- stats::lm(
    value ~ bin_factor + splines::ns(age, df = 3) + sex_factor,
    data = d, na.action = stats::na.fail
  )
  coefficient <- stats::coef(fit)
  covariance <- stats::vcov(fit)
  if (any(!is.finite(coefficient)) || any(!is.finite(covariance))) {
    stop("Age/sex-adjusted bin model is rank deficient for ", label, call. = FALSE)
  }
  model_terms <- stats::delete.response(stats::terms(fit))
  reference <- d[, .(age, sex_factor)]
  counts <- d[, .(n = .N, side = side[[1L]]), by = relative_bin]
  if (any(d[, .(side_n = data.table::uniqueN(side)), by = relative_bin]$side_n != 1L)) {
    stop("Each adjusted bin must map to one Yin/Yang side.", call. = FALSE)
  }
  adjusted <- data.table::rbindlist(lapply(bin_levels, function(bin_value) {
    newdata <- data.frame(
      bin_factor = factor(rep(bin_value, nrow(reference)), levels = bin_levels),
      age = reference$age,
      sex_factor = factor(reference$sex_factor, levels = levels(d$sex_factor))
    )
    design <- stats::model.matrix(model_terms, newdata, contrasts.arg = fit$contrasts)
    if (!identical(colnames(design), names(coefficient))) {
      stop("Adjusted marginal-mean design does not match the fitted model.", call. = FALSE)
    }
    design_mean <- colMeans(design)
    estimate <- sum(design_mean * coefficient)
    standard_error <- sqrt(drop(t(design_mean) %*% covariance %*% design_mean))
    data.table::data.table(
      relative_bin = bin_value, mean = estimate, se = standard_error,
      low = estimate - 1.96 * standard_error,
      high = estimate + 1.96 * standard_error
    )
  }))
  adjusted <- merge(adjusted, counts, by = "relative_bin", all.x = TRUE, sort = FALSE)
  adjusted[, `:=`(
    series = label, series_type = series_type,
    adjustment = "age natural cubic spline (df=3) + sex",
    standardization_population_n = nrow(d), bin_n = n
  )]
  adjusted[]
}

cross_fitted_age_sex_residual <- function(prediction) {
  required <- c("outer_fold", "score", "age", "sex")
  if (!all(required %in% names(prediction))) {
    stop("Age/sex residualization input is incomplete.", call. = FALSE)
  }
  age_sex_covariates(prediction$age, prediction$sex)
  result <- rep(NA_real_, nrow(prediction))
  for (fold in sort(unique(prediction$outer_fold))) {
    train_index <- which(prediction$outer_fold != fold)
    test_index <- which(prediction$outer_fold == fold)
    train <- data.frame(
      score = as.numeric(prediction$score[train_index]),
      age = as.numeric(prediction$age[train_index]),
      sex_factor = factor(prediction$sex[train_index], levels = c(0L, 1L),
                          labels = c("Female", "Male"))
    )
    test <- data.frame(
      age = as.numeric(prediction$age[test_index]),
      sex_factor = factor(prediction$sex[test_index], levels = c(0L, 1L),
                          labels = c("Female", "Male"))
    )
    fit <- stats::lm(score ~ splines::ns(age, df = 3) + sex_factor,
                     data = train, na.action = stats::na.fail)
    expected <- stats::predict(fit, newdata = test)
    result[test_index] <- as.numeric(prediction$score[test_index]) - expected
  }
  if (any(!is.finite(result))) stop("Cross-fitted age/sex residualization failed.", call. = FALSE)
  result
}

marginal_adjusted_bins <- function(dat, label, series_type, adjustment,
                                   rhs_terms) {
  required <- c("value", "relative_bin", "side")
  if (!all(required %in% names(dat))) {
    stop("Marginal-bin input is missing: ",
         paste(setdiff(required, names(dat)), collapse = ", "), call. = FALSE)
  }
  d <- data.table::copy(data.table::as.data.table(dat))
  d[, `:=`(
    value = suppressWarnings(as.numeric(value)),
    relative_bin = suppressWarnings(as.numeric(relative_bin)),
    side = as.character(side)
  )]
  if (any(!is.finite(d$value)) || any(!is.finite(d$relative_bin)) ||
      anyNA(d$side) || any(!stats::complete.cases(d))) {
    stop("Marginal-bin input must be complete and finite.", call. = FALSE)
  }
  bin_levels <- sort(unique(d$relative_bin))
  if (length(bin_levels) < 2L) stop("At least two bins are required.", call. = FALSE)
  d[, bin_factor := factor(relative_bin, levels = bin_levels)]
  rhs_terms <- trimws(as.character(rhs_terms[[1L]]))
  formula_text <- paste("value ~ bin_factor", if (nzchar(rhs_terms)) paste("+", rhs_terms) else "")
  fit <- stats::lm(stats::as.formula(formula_text), data = d,
                   na.action = stats::na.fail)
  coefficient <- stats::coef(fit)
  covariance <- stats::vcov(fit)
  if (any(!is.finite(coefficient)) || any(!is.finite(covariance))) {
    stop("Marginal-bin model is rank deficient for ", label,
         " at ", adjustment, call. = FALSE)
  }
  model_terms <- stats::delete.response(stats::terms(fit))
  counts <- d[, .(n = .N, side = side[[1L]]), by = relative_bin]
  if (any(d[, .(side_n = data.table::uniqueN(side)), by = relative_bin]$side_n != 1L)) {
    stop("Each adjusted bin must map to one Yin/Yang side.", call. = FALSE)
  }
  adjusted <- data.table::rbindlist(lapply(bin_levels, function(bin_value) {
    newdata <- as.data.frame(d)
    newdata$bin_factor <- factor(rep(bin_value, nrow(d)), levels = bin_levels)
    design <- stats::model.matrix(model_terms, newdata,
                                  contrasts.arg = fit$contrasts)
    if (!identical(colnames(design), names(coefficient))) {
      stop("Marginal-bin design does not match the fitted model.", call. = FALSE)
    }
    design_mean <- colMeans(design)
    estimate <- sum(design_mean * coefficient)
    standard_error <- sqrt(drop(t(design_mean) %*% covariance %*% design_mean))
    data.table::data.table(
      relative_bin = bin_value, mean = estimate, se = standard_error,
      low = estimate - 1.96 * standard_error,
      high = estimate + 1.96 * standard_error
    )
  }))
  adjusted <- merge(adjusted, counts, by = "relative_bin", all.x = TRUE,
                    sort = FALSE)
  adjusted[, `:=`(
    series = label, series_type = series_type, adjustment = adjustment,
    adjustment_formula = formula_text,
    standardization_population_n = nrow(d), bin_n = n
  )]
  adjusted[]
}

cross_fitted_covariate_residual <- function(prediction, rhs_terms) {
  required <- c("outer_fold", "score")
  if (!all(required %in% names(prediction))) {
    stop("Cross-fitted residualization input is incomplete.", call. = FALSE)
  }
  d <- data.table::copy(data.table::as.data.table(prediction))
  if (any(!stats::complete.cases(d))) {
    stop("Cross-fitted residualization requires a complete-case cohort.", call. = FALSE)
  }
  rhs_terms <- trimws(as.character(rhs_terms[[1L]]))
  if (!nzchar(rhs_terms)) return(as.numeric(d$score))
  formula_text <- paste("score ~", rhs_terms)
  result <- rep(NA_real_, nrow(d))
  for (fold in sort(unique(d$outer_fold))) {
    train_index <- which(d$outer_fold != fold)
    test_index <- which(d$outer_fold == fold)
    fit <- stats::lm(stats::as.formula(formula_text),
                     data = as.data.frame(d[train_index]),
                     na.action = stats::na.fail)
    coefficient <- stats::coef(fit)
    if (any(!is.finite(coefficient))) {
      stop("Cross-fitted covariate model is rank deficient in fold ", fold,
           call. = FALSE)
    }
    expected <- stats::predict(fit, newdata = as.data.frame(d[test_index]))
    result[test_index] <- as.numeric(d$score[test_index]) - expected
  }
  if (any(!is.finite(result))) {
    stop("Cross-fitted covariate residualization failed.", call. = FALSE)
  }
  result
}

fit_preprocessor <- function(x) {
  center <- colMeans(x, na.rm = TRUE)
  center[!is.finite(center)] <- 0
  for (j in seq_len(ncol(x))) {
    bad <- !is.finite(x[, j])
    if (any(bad)) x[bad, j] <- center[[j]]
  }
  x <- sweep(x, 2L, center, "-")
  scale <- sqrt(colSums(x * x) / pmax(1, nrow(x) - 1L))
  scale[!is.finite(scale) | scale < 1e-8] <- 1
  list(x = sweep(x, 2L, scale, "/"), center = center, scale = scale)
}

apply_preprocessor <- function(x, spec) {
  for (j in seq_len(ncol(x))) {
    bad <- !is.finite(x[, j])
    if (any(bad)) x[bad, j] <- spec$center[[j]]
  }
  sweep(sweep(x, 2L, spec$center, "-"), 2L, spec$scale, "/")
}

stratified_foldid <- function(y, k, seed) {
  y <- as.integer(y)
  if (!all(y %in% c(0L, 1L)) || min(table(y)) < k) {
    stop("Stratified-fold outcome contract failed.", call. = FALSE)
  }
  set.seed(as.integer(seed))
  foldid <- integer(length(y))
  for (value in c(0L, 1L)) {
    index <- sample(which(y == value))
    foldid[index] <- rep(seq_len(k), length.out = length(index))
  }
  foldid
}

known_horizon_outcome <- function(time, event, horizon = 5) {
  known <- (event == 1L & time <= horizon) | time > horizon
  label <- as.integer(event == 1L & time <= horizon)
  list(known = known, label = label)
}

km_censor <- function(time, event) {
  survival::survfit(survival::Surv(as.numeric(time), 1L - as.integer(event)) ~ 1)
}

km_value <- function(km, times, left = FALSE) {
  requested <- if (left) pmax(0, as.numeric(times) - 1e-10) else as.numeric(times)
  pmax(summary(km, times = requested, extend = TRUE)$surv, 1e-6)
}

ipcw_roc <- function(time, event, score, horizon = 5) {
  ok <- is.finite(time) & is.finite(score) & event %in% c(0L, 1L)
  time <- as.numeric(time[ok]); event <- as.integer(event[ok]); score <- as.numeric(score[ok])
  case <- event == 1L & time <= horizon
  control <- time > horizon
  if (!any(case) || !any(control)) stop("IPCW ROC lacks cases or controls.", call. = FALSE)
  km <- km_censor(time, event)
  weight <- numeric(length(time))
  weight[case] <- 1 / km_value(km, time[case], left = TRUE)
  weight[control] <- 1 / km_value(km, horizon)
  keep <- case | control
  d <- data.table::data.table(score = score[keep], case = case[keep], weight = weight[keep])
  grouped <- d[, .(case_weight = sum(weight * case), control_weight = sum(weight * !case)), by = score]
  data.table::setorder(grouped, -score)
  grouped[, `:=`(
    true_positive_rate = cumsum(case_weight) / sum(case_weight),
    false_positive_rate = cumsum(control_weight) / sum(control_weight)
  )]
  curve <- data.table::rbindlist(list(
    data.table::data.table(false_positive_rate = 0, true_positive_rate = 0),
    grouped[, .(false_positive_rate, true_positive_rate)],
    data.table::data.table(false_positive_rate = 1, true_positive_rate = 1)
  ))
  curve <- unique(curve, by = c("false_positive_rate", "true_positive_rate"))
  data.table::setorder(curve, false_positive_rate, true_positive_rate)
  auc <- sum(diff(curve$false_positive_rate) *
               (head(curve$true_positive_rate, -1L) + tail(curve$true_positive_rate, -1L)) / 2)
  list(curve = curve, auc = auc,
       cases5 = sum(case), controls5 = sum(control), censored_before5 = sum(!case & !control))
}

mean_fold_roc <- function(prediction, horizon = 5, grid = seq(0, 1, length.out = 2001L)) {
  folds <- sort(unique(prediction$outer_fold))
  fold_objects <- lapply(folds, function(fold) {
    one <- prediction[outer_fold == fold]
    roc <- ipcw_roc(one$time, one$event, one$score, horizon)
    interpolated <- stats::approx(
      x = roc$curve$false_positive_rate,
      y = roc$curve$true_positive_rate,
      xout = grid, method = "linear", ties = max, rule = 2
    )$y
    list(fold = fold, auc = roc$auc, tpr = interpolated)
  })
  mean_tpr <- rowMeans(do.call(cbind, lapply(fold_objects, `[[`, "tpr")))
  curve <- data.table::data.table(false_positive_rate = grid, true_positive_rate = mean_tpr)
  list(curve = curve, mean_fold_auc = mean(vapply(fold_objects, `[[`, numeric(1), "auc")),
       fold_auc = data.table::data.table(
         outer_fold = folds,
         AUC_5y = vapply(fold_objects, `[[`, numeric(1), "auc")
       ))
}

standardize_from_reference <- function(value, reference) {
  center <- mean(reference)
  scale <- stats::sd(reference)
  if (!is.finite(center) || !is.finite(scale) || scale <= 0) {
    stop("Invalid score reference moments.", call. = FALSE)
  }
  list(value = (as.numeric(value) - center) / scale, center = center, scale = scale)
}

event_year <- function(years) pmax(1L, as.integer(floor(as.numeric(years)) + 1L))

baseline_bin <- function(time_from_baseline) {
  direction <- ifelse(time_from_baseline < 0, -1L, 1L)
  absolute_year <- pmin(event_year(abs(time_from_baseline)), 16L)
  as.integer(direction * (2L * ((absolute_year - 1L) %/% 2L) + 1L))
}
