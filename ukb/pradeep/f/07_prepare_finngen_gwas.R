# Prepare local FinnGen R12 GWAS files for MR/coloc.

locate_repro_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")))
  frames <- sys.frames()
  ofiles <- vapply(frames, function(frame) if (!is.null(frame$ofile)) frame$ofile else NA_character_, character(1))
  ofiles <- ofiles[!is.na(ofiles)]
  if (length(ofiles) > 0) return(dirname(normalizePath(ofiles[length(ofiles)], winslash = "/")))
  project_dir <- Sys.getenv("PRADEEP_PROJECT_DIR", unset = "")
  if (nzchar(project_dir)) return(file.path(project_dir, "f"))
  stop("Cannot locate pradeep/f. Run this step through pradeep.sh.", call. = FALSE)
}

this_dir <- locate_repro_script_dir()
source(file.path(this_dir, "00_config.R"))
source(file.path(script_dir, "R", "repro_utils.R"))
repro_load_packages(c("data.table"))

message("Step 7: preparing FinnGen R12 GWAS files")

download_dir <- normalizePath(
  Sys.getenv(
    "PRADEEP_FINNGEN_DOWNLOAD_DIR",
    unset = Sys.getenv("UKBPPP_FINNGEN_DOWNLOAD_DIR", unset = file.path(root_dir, "data", "gwas", "downloads"))
  ),
  winslash = "/",
  mustWork = FALSE
)
gwas_root <- repro_env_path("PRADEEP_GWAS_DIR", file.path(root_dir, "data", "gwas"))
force_standardize <- repro_bool_env("UKBPPP_FORCE_GWAS_STANDARDIZE", default = FALSE)
finngen_n_total <- as.numeric(Sys.getenv("UKBPPP_FINNGEN_N_TOTAL", unset = "500348"))
if (!is.finite(finngen_n_total) || finngen_n_total <= 0) finngen_n_total <- 500348

gwas_manifest <- data.table::data.table(
  outcome = c("CAD", "HF", "AF", "AS"),
  local_outcome_key = c("cad", "hfail", "afib", "ao_sten"),
  finngen_r12_endpoint = c("I9_CHD", "I9_HEARTFAIL", "I9_AF", "I9_CAVS_OPERATED"),
  phenotype = c(
    "Major coronary heart disease event",
    "Heart failure",
    "Atrial fibrillation and flutter",
    "Calcific aortic valve stenosis, operated"
  ),
  local_dir = file.path(gwas_root, c("cad_eur", "hfail_eur", "afib_eur", "ao_sten_eur"))
)
gwas_manifest[, file_name := paste0("finngen_R12_", finngen_r12_endpoint, ".gz")]
gwas_manifest[, raw_local_path := file.path(local_dir, file_name)]
gwas_manifest[, raw_download_path := file.path(download_dir, file_name)]
gwas_manifest[, standard_path := file.path(local_dir, sub("\\.gz$", ".standard.tsv.gz", file_name))]
gwas_manifest[, source_tsv := file.path(local_dir, sub("\\.gz$", ".source.tsv", file_name))]

standardize_one_finngen <- function(info) {
  dir.create(info$local_dir, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(info$standard_path) && !force_standardize) {
    message("  existing standard file; skipping ", info$outcome, ": ", info$standard_path)
    return(data.table::data.table(
      outcome = info$outcome,
      raw_path = NA_character_,
      standard_path = info$standard_path,
      n_rows = NA_integer_,
      n_case = NA_integer_,
      n_control = NA_integer_,
      status = "existing_standard_file"
    ))
  }

  raw_candidates <- unique(c(info$raw_local_path, info$raw_download_path))
  raw_path <- raw_candidates[file.exists(raw_candidates)][1]
  if (is.na(raw_path) || length(raw_path) == 0) {
    stop(
      "Missing FinnGen raw GWAS for ", info$outcome, ". Expected one of:\n  ",
      paste(raw_candidates, collapse = "\n  "),
      call. = FALSE
    )
  }

  message("  reading ", info$outcome, ": ", raw_path)
  raw_cols <- c("#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls")
  x <- data.table::fread(raw_path, select = raw_cols, showProgress = TRUE, nThread = max(1L, data.table::getDTthreads()))
  missing_cols <- setdiff(raw_cols, names(x))
  if (length(missing_cols) > 0) {
    stop("FinnGen file is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  denom <- x$af_alt_cases - x$af_alt_controls
  case_fraction <- (x$af_alt - x$af_alt_controls) / denom
  case_fraction <- case_fraction[is.finite(case_fraction) & case_fraction > 0 & case_fraction < 1 & abs(denom) > 1e-5]
  case_fraction <- stats::median(case_fraction, na.rm = TRUE)
  if (!is.finite(case_fraction)) {
    warning("Could not estimate case fraction for ", info$outcome, "; using 0.10.")
    case_fraction <- 0.10
  }
  n_case <- as.integer(round(finngen_n_total * case_fraction))
  n_control <- as.integer(finngen_n_total - n_case)

  x[, SNP := data.table::tstrsplit(as.character(rsids), ",", fixed = TRUE, keep = 1L)[[1L]]]
  x[is.na(SNP) | SNP == "" | SNP == ".", SNP := paste(`#chrom`, pos, ref, alt, sep = ":")]
  out <- x[, .(
    SNP = as.character(SNP),
    CHR = as.character(`#chrom`),
    POS = as.integer(pos),
    EA = as.character(alt),
    NEA = as.character(ref),
    EAF = as.numeric(af_alt),
    BETA = as.numeric(beta),
    SE = as.numeric(sebeta),
    P = as.numeric(pval)
  )]
  out[CHR == "23", CHR := "X"]
  out <- out[
    !is.na(SNP) & SNP != "" &
      is.finite(POS) & is.finite(BETA) & is.finite(SE) & SE > 0 &
      is.finite(P) & P > 0 & P <= 1 &
      is.finite(EAF) & EAF > 0 & EAF < 1
  ]
  out <- out[!duplicated(SNP)]
  out[, `:=`(
    N = as.integer(finngen_n_total),
    N_CASE = n_case,
    N_CONTROL = n_control
  )]

  message("  writing ", info$outcome, ": ", info$standard_path)
  data.table::fwrite(out, info$standard_path)
  source_dt <- data.table::data.table(
    field = c("source", "phenocode", "phenotype", "raw_path", "standard_path", "n_total", "n_case", "n_control", "case_fraction_estimation"),
    value = c(
      "FinnGen R12",
      info$finngen_r12_endpoint,
      info$phenotype,
      raw_path,
      info$standard_path,
      as.character(finngen_n_total),
      as.character(n_case),
      as.character(n_control),
      "Estimated from af_alt, af_alt_cases and af_alt_controls; N total default 500348."
    )
  )
  data.table::fwrite(source_dt, info$source_tsv, sep = "\t")
  data.table::data.table(
    outcome = info$outcome,
    raw_path = raw_path,
    standard_path = info$standard_path,
    n_rows = nrow(out),
    n_case = n_case,
    n_control = n_control,
    status = "standardized"
  )
}

status <- data.table::rbindlist(lapply(seq_len(nrow(gwas_manifest)), function(i) {
  standardize_one_finngen(gwas_manifest[i])
}), fill = TRUE)

status_path <- file.path(output_dir, "finngen_r12_gwas_prepare_status.csv")
data.table::fwrite(status, status_path)
message("FinnGen prepare status saved: ", status_path)
print(status)
