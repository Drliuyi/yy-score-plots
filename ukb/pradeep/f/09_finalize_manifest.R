# Build a single result and figure inventory after completed analysis steps.

analysis_dir <- Sys.getenv("PRADEEP_ANALYSIS_DIR", unset = "")
if (!nzchar(analysis_dir)) stop("PRADEEP_ANALYSIS_DIR is not set.", call. = FALSE)
output_dir <- file.path(analysis_dir, "outputs")
audit_dir <- file.path(analysis_dir, "audit")
fig_dir <- file.path(output_dir, "figures")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

all_files <- if (dir.exists(output_dir)) {
  list.files(output_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
} else {
  character()
}
all_files <- all_files[file.info(all_files)$isdir %in% FALSE]
relative_path <- function(path, root) {
  path <- gsub("\\\\", "/", path)
  prefix <- paste0(sub("/+$", "", gsub("\\\\", "/", root)), "/")
  ifelse(startsWith(path, prefix), substring(path, nchar(prefix) + 1L), path)
}

if (length(all_files) > 0L) {
  info <- file.info(all_files)
  output_manifest <- data.frame(
    file = basename(all_files),
    relative_path = relative_path(all_files, output_dir),
    size_bytes = info$size,
    modified = format(info$mtime, "%Y-%m-%d %H:%M:%S %z"),
    md5 = unname(tools::md5sum(all_files)),
    stringsAsFactors = FALSE
  )
  output_manifest <- output_manifest[order(output_manifest$relative_path), , drop = FALSE]
} else {
  output_manifest <- data.frame(
    file = character(), relative_path = character(), size_bytes = numeric(),
    modified = character(), md5 = character(), stringsAsFactors = FALSE
  )
}
write.csv(output_manifest, file.path(audit_dir, "output_manifest.csv"), row.names = FALSE, na = "")

figure_files <- all_files[grepl("\\.(png|pdf|svg|tif|tiff)$", all_files, ignore.case = TRUE)]
if (length(figure_files) > 0L) {
  finfo <- file.info(figure_files)
  figure_manifest <- data.frame(
    file = basename(figure_files),
    relative_path = relative_path(figure_files, output_dir),
    size_bytes = finfo$size,
    modified = format(finfo$mtime, "%Y-%m-%d %H:%M:%S %z"),
    stringsAsFactors = FALSE
  )
  figure_manifest <- figure_manifest[order(figure_manifest$relative_path), , drop = FALSE]
} else {
  figure_manifest <- data.frame(
    file = character(), relative_path = character(), size_bytes = numeric(),
    modified = character(), stringsAsFactors = FALSE
  )
}
write.csv(figure_manifest, file.path(audit_dir, "figure_manifest_complete.csv"), row.names = FALSE, na = "")

main_required_all <- c(
  "original_style_primary_manhattan_all.png",
  "original_style_primary_density_insets.png",
  "original_style_primary_manhattan_cad.png",
  "original_style_primary_manhattan_hf.png",
  "original_style_primary_manhattan_af.png",
  "original_style_primary_manhattan_as.png",
  "original_style_sex_difference_alloutcomes.png",
  "original_style_sex_hr_top5.png",
  "original_style_sex_difference_with_hr_insets.png",
  "original_style_lasso_rocs.png",
  "original_style_lasso_prediction_figure.png",
  "original_style_lasso_prediction_metrics.png",
  "original_style_lasso_calibration.png",
  "original_style_lasso_protein_features.png"
)
main_required <- if (file.exists(file.path(fig_dir, "paper_figure_qc.csv"))) main_required_all else character()
go_required <- if (file.exists(file.path(output_dir, "enrichment", "go_enrichment_all_terms.csv"))) {
  "original_style_go_enrichment_all.png"
} else {
  character()
}
mr_required <- if (file.exists(file.path(output_dir, "mr", "local_pqtl_mr_all_outcomes.csv"))) {
  "original_style_local_pqtl_mr_all.png"
} else {
  character()
}
coloc_file <- file.path(output_dir, "coloc", "local_pqtl_coloc_all_outcomes.csv")
coloc_required <- character()
if (file.exists(coloc_file)) {
  coloc_data <- tryCatch(utils::read.csv(coloc_file, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(coloc_data) && "status" %in% names(coloc_data) && any(coloc_data$status == "ok", na.rm = TRUE)) {
    coloc_required <- "original_style_local_pqtl_coloc_pph4_all.png"
  }
}
required <- unique(c(main_required, go_required, mr_required, coloc_required))
figure_qc <- data.frame(
  figure = required,
  exists = file.exists(file.path(fig_dir, required)),
  path = file.path(fig_dir, required),
  stringsAsFactors = FALSE
)
write.csv(figure_qc, file.path(audit_dir, "figure_completeness.csv"), row.names = FALSE, na = "")

cat("\nResult inventory\n")
cat("  analysis_dir: ", analysis_dir, "\n", sep = "")
cat("  output_files: ", nrow(output_manifest), "\n", sep = "")
cat("  figure_files: ", nrow(figure_manifest), "\n", sep = "")
cat("  required_figures_present: ", sum(figure_qc$exists), "/", nrow(figure_qc), "\n", sep = "")
if (any(!figure_qc$exists)) {
  cat("  missing_figures: ", paste(figure_qc$figure[!figure_qc$exists], collapse = ", "), "\n", sep = "")
  stop("FIGURE_COMPLETENESS_FAILED: required result-backed figures are missing.", call. = FALSE)
}
