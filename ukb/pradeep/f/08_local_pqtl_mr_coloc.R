# Paper-style local MR and coloc for all four cardiac outcomes.
# This follows the upstream scripts as closely as local inputs allow:
#   - observational Bonferroni hits are used as protein-disease candidates
#   - cis pQTL instruments are selected at P < 5e-6
#   - PLINK2 clumping uses 1000G EUR with kb=10000 and r2=0.1
#   - MR uses Wald ratio for one SNP and correlated IVW for >=2 SNPs
#   - Egger and Egger intercept are also reported when >=3 SNPs are available
#   - coloc is run for MR-significant protein-disease pairs with an LD matrix

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
repro_load_packages(c(
  "data.table", "dplyr", "tidyr", "ggplot2", "ggrepel", "patchwork",
  "TwoSampleMR", "MendelianRandomization", "coloc"
))

message("Step 8: paper-style local pQTL MR/coloc for all outcomes")

mr_dir <- file.path(output_dir, "mr")
coloc_dir <- file.path(output_dir, "coloc")
fig_dir <- file.path(output_dir, "figures")
ld_tmp_dir <- file.path(output_dir, "tmp_plink2")
dir.create(mr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(coloc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ld_tmp_dir, recursive = TRUE, showWarnings = FALSE)
gwas_root <- repro_env_path("PRADEEP_GWAS_DIR", file.path(root_dir, "data", "gwas"))
pqtl_root <- repro_env_path("PRADEEP_PQTL_DIR", file.path(root_dir, "data.BIG", "gwas", "ppp", "clean"))

outcome_info <- data.table::data.table(
  Outcome = c("CAD", "HF", "Afib", "AS"),
  outc = c("CAD", "HF", "AF", "AS"),
  outc_lower = c("cad", "hf", "af", "as"),
  title = c("Coronary artery disease", "Heart failure", "Atrial fibrillation", "Aortic stenosis"),
  gwas_path = file.path(gwas_root, c("cad_eur", "hfail_eur", "afib_eur", "ao_sten_eur"), c(
    "finngen_R12_I9_CHD.standard.tsv.gz",
    "finngen_R12_I9_HEARTFAIL.standard.tsv.gz",
    "finngen_R12_I9_AF.standard.tsv.gz",
    "finngen_R12_I9_CAVS_OPERATED.standard.tsv.gz"
  )),
  paper_N = c(377277, 377277, 237690, 377277),
  paper_s = c(0.115347609, 0.072371229, 0.192544911, 0.024260689)
)

missing_gwas <- outcome_info[!file.exists(gwas_path)]
if (nrow(missing_gwas) > 0) {
  stop(
    "Missing standardized GWAS files. Run 10_prepare_finngen_gwas.R first.\n  ",
    paste(missing_gwas$gwas_path, collapse = "\n  "),
    call. = FALSE
  )
}

primary <- data.table::fread(primary_assoc_file)
repro_assert_primary_complete(primary, audit_dir, outcome_map, context = "Primary association results used for MR/coloc")
primary <- primary[Outcome %in% outcome_info$Outcome]
primary[, outc := outcome_info$outc[match(as.character(Outcome), outcome_info$Outcome)]]

mr_p_threshold <- as.numeric(Sys.getenv("UKBPPP_LOCAL_MR_P_THRESHOLD", unset = "5e-6"))
if (!is.finite(mr_p_threshold) || mr_p_threshold <= 0) mr_p_threshold <- 5e-6
mr_r2 <- as.numeric(Sys.getenv("UKBPPP_LOCAL_MR_R2", unset = "0.1"))
if (!is.finite(mr_r2) || mr_r2 <= 0 || mr_r2 > 1) mr_r2 <- 0.1
mr_kb <- repro_int_env("UKBPPP_LOCAL_MR_KB", 10000L)
mr_max_proteins <- repro_int_env("UKBPPP_LOCAL_MR_MAX_PROTEINS_PER_OUTCOME", 0L)
coloc_max_proteins <- repro_int_env("UKBPPP_COLOC_MAX_PROTEINS", 0L)
run_mr_sensitivity <- repro_bool_env("UKBPPP_RUN_MR_SENSITIVITY", default = TRUE)
parse_numeric_vector_env <- function(name, default) {
  raw <- Sys.getenv(name, unset = paste(default, collapse = ","))
  vals <- suppressWarnings(as.numeric(strsplit(raw, "[,;[:space:]]+")[[1]]))
  vals <- vals[is.finite(vals) & vals > 0]
  if (length(vals) == 0) default else unique(vals)
}
mr_sensitivity_p_thresholds <- parse_numeric_vector_env(
  "UKBPPP_MR_SENSITIVITY_P_THRESHOLDS",
  c(5e-4, 5e-6, 5e-8)
)
mr_sensitivity_r2_thresholds <- parse_numeric_vector_env(
  "UKBPPP_MR_SENSITIVITY_R2_THRESHOLDS",
  c(0.001, 0.01, 0.1, 0.2, 0.4)
)
plink2_wsl <- Sys.getenv("UKBPPP_PLINK2_WSL", unset = unname(Sys.which("plink2")))
ld_pfile_template <- Sys.getenv(
  "UKBPPP_LD_PFILE_TEMPLATE",
  unset = ""
)
ld_bfile <- trimws(Sys.getenv("UKBPPP_LD_BFILE", unset = ""))

windows_to_wsl_path <- function(path) {
  if (.Platform$OS.type != "windows") return(path)
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  out <- suppressWarnings(system2("wsl.exe", c("wslpath", "-a", path), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop("wslpath failed for path: ", path, "\n", paste(out, collapse = "\n"), call. = FALSE)
  }
  trimws(out[1])
}

run_wsl_plink2 <- function(args) {
  if (.Platform$OS.type != "windows") {
    out <- suppressWarnings(system2(plink2_wsl, args, stdout = TRUE, stderr = TRUE))
  } else {
    out <- suppressWarnings(system2("wsl.exe", c(plink2_wsl, args), stdout = TRUE, stderr = TRUE))
  }
  status <- attr(out, "status")
  list(status = ifelse(is.null(status), 0L, as.integer(status)), output = out)
}

check_plink2 <- function() {
  res <- run_wsl_plink2("--version")
  if (res$status != 0) {
    stop(
      "Cannot run PLINK2. Set UKBPPP_PLINK2_WSL to the WSL plink2 path.\n",
      paste(res$output, collapse = "\n"),
      call. = FALSE
    )
  }
  message("Using PLINK2: ", paste(res$output[1], collapse = ""))
}
check_plink2()

clean_variant_id <- function(x) {
  x <- as.character(x)
  x <- sub(",.*$", "", x)
  trimws(x)
}

pfile_prefix_for_chr <- function(chr) {
  chr <- as.character(chr)
  chr <- ifelse(chr %in% c("23", "X", "x"), "X", chr)
  prefix <- gsub("{CHR}", chr, ld_pfile_template, fixed = TRUE)
  if (.Platform$OS.type == "windows") {
    miss <- !file.exists(paste0(prefix, ".pgen")) ||
      !file.exists(paste0(prefix, ".pvar")) ||
      !file.exists(paste0(prefix, ".psam"))
    if (miss) stop("Missing LD pfile prefix for chr ", chr, ": ", prefix, call. = FALSE)
  }
  windows_to_wsl_path(prefix)
}

ld_reference_args <- function(chr) {
  if (nzchar(ld_bfile)) {
    if (.Platform$OS.type == "windows") {
      missing <- !file.exists(paste0(ld_bfile, ".bed")) ||
        !file.exists(paste0(ld_bfile, ".bim")) ||
        !file.exists(paste0(ld_bfile, ".fam"))
      if (missing) stop("Missing LD bfile prefix: ", ld_bfile, call. = FALSE)
    }
    return(c("--bfile", windows_to_wsl_path(ld_bfile)))
  }
  c("--pfile", pfile_prefix_for_chr(chr))
}

read_pqtl_cis <- function(protein) {
  pdir <- file.path(pqtl_root, protein)
  cma_file <- file.path(pdir, paste0(protein, ".cma.cojo"))
  gcta_file <- file.path(pdir, paste0(protein, ".4gcta"))
  if (!file.exists(cma_file) || !file.exists(gcta_file)) return(NULL)

  cma <- data.table::fread(cma_file)
  gcta <- data.table::fread(gcta_file)
  if (!all(c("Chr", "SNP", "bp", "freq", "b", "se", "p", "n") %in% names(cma))) return(NULL)
  if (!all(c("SNP", "EA", "NEA", "EAF", "BETA", "SE", "P", "N") %in% names(gcta))) return(NULL)

  cma <- cma[, .(
    SNP = clean_variant_id(SNP),
    CHR = as.character(Chr),
    POS = as.integer(bp),
    EAF_cma = as.numeric(freq),
    BETA_cma = as.numeric(b),
    SE_cma = as.numeric(se),
    P_cma = as.numeric(p),
    N_cma = as.numeric(n)
  )]
  gcta <- gcta[, .(
    SNP = clean_variant_id(SNP),
    EA = toupper(as.character(EA)),
    NEA = toupper(as.character(NEA)),
    EAF = as.numeric(EAF),
    BETA = as.numeric(BETA),
    SE = as.numeric(SE),
    P = as.numeric(P),
    N = as.numeric(N)
  )]
  dat <- merge(cma, gcta, by = "SNP", all.x = TRUE)
  dat[is.na(EAF), EAF := EAF_cma]
  dat[is.na(BETA), BETA := BETA_cma]
  dat[is.na(SE), SE := SE_cma]
  dat[is.na(P), P := P_cma]
  dat[is.na(N), N := N_cma]
  dat <- dat[
    !is.na(SNP) & SNP != "" & SNP != "." &
      !is.na(CHR) & is.finite(POS) &
      !is.na(EA) & !is.na(NEA) & EA != "" & NEA != "" &
      is.finite(BETA) & is.finite(SE) & SE > 0 &
      is.finite(EAF) & EAF > 0 & EAF < 1 &
      is.finite(P) & P >= 0 & is.finite(N)
  ]
  dat[, `:=`(
    Protein = protein,
    P_fmt = pmax(P, .Machine$double.xmin)
  )]
  dat <- dat[order(P_fmt)]
  dat <- dat[!duplicated(SNP)]
  if (nrow(dat) == 0) return(NULL)
  dat
}

candidate_hits <- function(outcome_label, max_n = 0L) {
  hits <- primary[Outcome == outcome_label & Bonferroni %in% TRUE][order(P_Value)]
  hits <- hits[file.exists(file.path(pqtl_root, Protein, paste0(Protein, ".cma.cojo"))) &
    file.exists(file.path(pqtl_root, Protein, paste0(Protein, ".4gcta")))]
  if (max_n > 0L) hits <- hits[seq_len(min(.N, max_n))]
  hits
}

read_gwas_subset <- function(gwas_path, needed_keys = character(), needed_snps = character()) {
  cols <- c("SNP", "CHR", "POS", "EA", "NEA", "EAF", "BETA", "SE", "P", "N", "N_CASE", "N_CONTROL")
  gwas <- data.table::fread(gwas_path, select = cols, showProgress = TRUE, nThread = max(1L, data.table::getDTthreads()))
  missing_cols <- setdiff(cols, names(gwas))
  if (length(missing_cols) > 0) stop("GWAS missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  gwas[, `:=`(
    SNP = clean_variant_id(SNP),
    CHR = as.character(CHR),
    POS = as.integer(POS),
    EA = toupper(as.character(EA)),
    NEA = toupper(as.character(NEA)),
    key = paste(CHR, POS, sep = ":"),
    P_fmt = pmax(as.numeric(P), .Machine$double.xmin)
  )]
  gwas <- gwas[key %chin% needed_keys | SNP %chin% needed_snps]
  gwas <- gwas[!duplicated(paste(CHR, POS, EA, NEA, sep = ":"))]
  gwas
}

make_harmonised <- function(pdat, gwas, outcome_lower) {
  if (is.null(pdat) || nrow(pdat) == 0 || nrow(gwas) == 0) return(data.table::data.table())
  g <- gwas[CHR %in% unique(pdat$CHR) & (POS %in% pdat$POS | SNP %chin% pdat$SNP)]
  if (nrow(g) == 0) return(data.table::data.table())

  outcome_overlap <- copy(g)
  outcome_overlap[, `:=`(
    phen = outcome_lower,
    id = paste(CHR, POS, EA, NEA, sep = ":")
  )]
  out_fmt <- TwoSampleMR::format_data(
    as.data.frame(outcome_overlap),
    type = "outcome",
    phenotype_col = "phen",
    snp_col = "id",
    beta_col = "BETA",
    se_col = "SE",
    eaf_col = "EAF",
    effect_allele_col = "EA",
    other_allele_col = "NEA",
    pval_col = "P_fmt",
    chr_col = "CHR",
    pos_col = "POS",
    samplesize_col = "N"
  )

  chrom_overlap <- pdat[POS %in% outcome_overlap$POS | SNP %chin% outcome_overlap$SNP]
  if (nrow(chrom_overlap) == 0) return(data.table::data.table())
  chrom_overlap_2 <- copy(chrom_overlap)
  old_ea <- chrom_overlap_2$EA
  old_nea <- chrom_overlap_2$NEA
  chrom_overlap_2[, `:=`(
    BETA = -BETA,
    EAF = 1 - EAF,
    EA = old_nea,
    NEA = old_ea
  )]
  chrom_overlap <- data.table::rbindlist(list(chrom_overlap, chrom_overlap_2), use.names = TRUE, fill = TRUE)
  chrom_overlap[, `:=`(
    id = paste(CHR, POS, EA, NEA, sep = ":"),
    phen = unique(Protein)[1]
  )]
  exp_fmt <- TwoSampleMR::format_data(
    as.data.frame(chrom_overlap),
    type = "exposure",
    phenotype_col = "phen",
    snp_col = "id",
    beta_col = "BETA",
    se_col = "SE",
    eaf_col = "EAF",
    effect_allele_col = "EA",
    other_allele_col = "NEA",
    pval_col = "P_fmt",
    chr_col = "CHR",
    pos_col = "POS",
    samplesize_col = "N"
  )

  dat <- suppressMessages(TwoSampleMR::harmonise_data(exposure_dat = exp_fmt, outcome_dat = out_fmt))
  dat <- data.table::as.data.table(dat)
  if (nrow(dat) == 0) return(dat)
  rsid_map <- unique(outcome_overlap[, .(
    pos.exposure = POS,
    rsids = SNP,
    outcome_key = key,
    N_CASE,
    N_CONTROL
  )])
  dat <- merge(dat, rsid_map, by = "pos.exposure", all.x = TRUE)
  data.table::setnames(dat, "SNP", "pos_id")
  data.table::setnames(dat, "rsids", "SNP")
  dat <- dat[!is.na(SNP) & SNP != "" & SNP != "."]
  dat <- dat[order(pval.exposure)]
  dat <- dat[!duplicated(SNP)]
  dat <- dat[!duplicated(pos_id)]
  dat
}

plink2_clump <- function(dat_u, p_threshold = mr_p_threshold, r2 = mr_r2, kb = mr_kb) {
  dat_u <- data.table::as.data.table(dat_u)
  dat_u <- dat_u[
    !is.na(SNP) & SNP != "" & SNP != "." &
      is.finite(pval.exposure) &
      !is.na(effect_allele.exposure) & effect_allele.exposure != ""
  ]
  if (nrow(dat_u) == 0) return(character())
  clumped <- character()
  for (chr in unique(dat_u$chr.exposure)) {
    dchr <- dat_u[chr.exposure == chr]
    if (nrow(dchr) == 0) next
    tmp_prefix <- file.path(ld_tmp_dir, paste0("clump_chr", chr, "_", as.integer(stats::runif(1, 1, 1e9))))
    clump_file <- paste0(tmp_prefix, ".txt")
    clump_input <- dchr[, .(
      SNP,
      P = pval.exposure,
      A1 = toupper(effect_allele.exposure)
    )]
    data.table::fwrite(clump_input, clump_file, sep = "\t", quote = FALSE)
    args <- c(
      ld_reference_args(chr),
      "--rm-dup", "force-first",
      "--clump", windows_to_wsl_path(clump_file),
      "--clump-p1", format(p_threshold, scientific = TRUE),
      "--clump-r2", as.character(r2),
      "--clump-kb", as.character(kb),
      "--out", windows_to_wsl_path(tmp_prefix)
    )
    res <- run_wsl_plink2(args)
    if (res$status != 0) {
      stop("PLINK2 clump failed for chr ", chr, "\n", paste(res$output, collapse = "\n"), call. = FALSE)
    }
    out_file <- paste0(tmp_prefix, ".clumps")
    if (file.exists(out_file)) {
      cl <- data.table::fread(out_file)
      if ("ID" %in% names(cl)) clumped <- c(clumped, as.character(cl$ID))
    }
  }
  unique(clumped)
}

plink2_ld_raw <- function(snps, chr) {
  snps <- unique(snps[!is.na(snps) & snps != "" & snps != "."])
  if (length(snps) == 0) return(list(ld = matrix(numeric(0), 0, 0), pvar = data.table::data.table()))
  tmp_prefix <- file.path(ld_tmp_dir, paste0("ld_chr", chr, "_", as.integer(stats::runif(1, 1, 1e9))))
  extract_file <- paste0(tmp_prefix, ".extract")
  data.table::fwrite(data.table::data.table(SNP = snps), extract_file, col.names = FALSE)
  args <- c(
    ld_reference_args(chr),
    "--rm-dup", "force-first",
    "--extract", windows_to_wsl_path(extract_file),
    "--r-unphased", "square",
    "--make-just-pvar",
    "--out", windows_to_wsl_path(tmp_prefix)
  )
  res <- run_wsl_plink2(args)
  if (res$status != 0) {
    stop("PLINK2 LD matrix failed for chr ", chr, "\n", paste(res$output, collapse = "\n"), call. = FALSE)
  }
  vars_file <- paste0(tmp_prefix, ".unphased.vcor1.vars")
  mat_file <- paste0(tmp_prefix, ".unphased.vcor1")
  pvar_file <- paste0(tmp_prefix, ".pvar")
  if (!file.exists(vars_file) || !file.exists(mat_file)) {
    return(list(ld = matrix(numeric(0), 0, 0), pvar = data.table::data.table()))
  }
  vars <- readLines(vars_file, warn = FALSE)
  vars <- vars[nzchar(vars)]
  if (length(vars) == 0) return(list(ld = matrix(numeric(0), 0, 0), pvar = data.table::data.table()))
  mat <- as.matrix(data.table::fread(mat_file, header = FALSE))
  if (length(vars) == 1L) mat <- matrix(as.numeric(mat[1, 1]), nrow = 1L)
  storage.mode(mat) <- "numeric"
  rownames(mat) <- vars
  colnames(mat) <- vars
  pvar <- data.table::fread(pvar_file, skip = "#CHROM")
  data.table::setnames(pvar, old = names(pvar)[1:5], new = c("CHR", "POS37", "SNP", "REF", "ALT"))
  pvar[, `:=`(
    SNP = as.character(SNP),
    REF = toupper(as.character(REF)),
    ALT = toupper(as.character(ALT))
  )]
  list(ld = mat, pvar = pvar)
}

build_ld_matrix <- function(dat) {
  dat <- data.table::as.data.table(dat)
  dat <- dat[!is.na(SNP) & !is.na(chr.exposure)]
  if (nrow(dat) == 0) return(list(dat = dat, ld = matrix(numeric(0), 0, 0)))
  blocks <- list()
  kept <- list()
  for (chr in unique(dat$chr.exposure)) {
    dchr <- dat[chr.exposure == chr]
    raw <- plink2_ld_raw(dchr$SNP, chr)
    ld <- raw$ld
    if (nrow(ld) == 0) next
    dchr <- dchr[SNP %in% rownames(ld)]
    dchr <- dchr[match(rownames(ld), SNP)]
    pvar <- raw$pvar[match(dchr$SNP, SNP)]
    sign <- rep(1, nrow(dchr))
    sign[toupper(dchr$effect_allele.exposure) == pvar$REF] <- -1
    sign[toupper(dchr$effect_allele.exposure) == pvar$ALT] <- 1
    ld <- sweep(sweep(ld, 1, sign, `*`), 2, sign, `*`)
    blocks[[as.character(chr)]] <- ld
    kept[[as.character(chr)]] <- dchr
  }
  if (length(blocks) == 0) return(list(dat = dat[0], ld = matrix(numeric(0), 0, 0)))
  dat_kept <- data.table::rbindlist(kept, fill = TRUE)
  if (length(blocks) == 1L) return(list(dat = dat_kept, ld = blocks[[1L]]))
  n <- sum(vapply(blocks, nrow, integer(1)))
  out <- matrix(0, n, n)
  rn <- character()
  start <- 1L
  for (b in blocks) {
    idx <- start:(start + nrow(b) - 1L)
    out[idx, idx] <- b
    rn <- c(rn, rownames(b))
    start <- start + nrow(b)
  }
  rownames(out) <- rn
  colnames(out) <- rn
  list(dat = dat_kept, ld = out)
}

run_mr_estimates <- function(dat, protein, outc_lower, p_threshold = mr_p_threshold, r2 = mr_r2) {
  dat <- data.table::as.data.table(dat)
  dat <- dat[mr_keep %in% TRUE]
  if (nrow(dat) == 0) return(list(results = data.table::data.table(), instruments = data.table::data.table()))

  if (nrow(dat) == 1L) {
    results_mr <- TwoSampleMR::mr(as.data.frame(dat), method_list = c("mr_wald_ratio"))
    res <- data.table::data.table(
      exp = protein,
      outc = outc_lower,
      pvalthreshold = p_threshold,
      rsqthreshold = r2,
      nsnp = results_mr$nsnp,
      method = results_mr$method,
      b = results_mr$b,
      se = results_mr$se,
      pval = results_mr$pval
    )
    dat[, marker_ld := paste(SNP, effect_allele.exposure, other_allele.exposure, sep = "_")]
    return(list(results = res, instruments = dat))
  }

  ld_obj <- build_ld_matrix(dat)
  dat <- ld_obj$dat
  ld <- ld_obj$ld
  if (nrow(dat) == 0 || nrow(ld) == 0) return(list(results = data.table::data.table(), instruments = data.table::data.table()))
  if (nrow(dat) == 1L) return(run_mr_estimates(dat, protein, outc_lower, p_threshold = p_threshold, r2 = r2))

  dat <- dat[match(rownames(ld), SNP)]
  dat[, marker_ld := paste(SNP, effect_allele.exposure, other_allele.exposure, sep = "_")]
  dat2 <- MendelianRandomization::mr_input(
    bx = dat$beta.exposure,
    bxse = dat$se.exposure,
    by = dat$beta.outcome,
    byse = dat$se.outcome,
    correlation = ld
  )
  output_mr_ivw_corr <- tryCatch(
    MendelianRandomization::mr_ivw(dat2, correl = TRUE),
    error = function(e) e
  )
  if (inherits(output_mr_ivw_corr, "error")) {
    return(list(results = data.table::data.table(), instruments = dat))
  }
  results <- data.table::data.table(
    exp = protein,
    outc = outc_lower,
    pvalthreshold = p_threshold,
    rsqthreshold = r2,
    nsnp = output_mr_ivw_corr@SNPs,
    method = "Inverse variance weighted (correlation inc)",
    b = output_mr_ivw_corr@Estimate,
    se = output_mr_ivw_corr@StdError,
    pval = output_mr_ivw_corr@Pvalue
  )
  if (nrow(dat) >= 3L) {
    output_mr_egger_corr <- tryCatch(
      MendelianRandomization::mr_egger(dat2, correl = TRUE),
      error = function(e) e
    )
    if (!inherits(output_mr_egger_corr, "error")) {
      results <- data.table::rbindlist(list(
        results,
        data.table::data.table(
          exp = protein,
          outc = outc_lower,
          pvalthreshold = p_threshold,
          rsqthreshold = r2,
          nsnp = output_mr_egger_corr@SNPs,
          method = "Egger (correlation inc)",
          b = output_mr_egger_corr@Estimate,
          se = output_mr_egger_corr@StdError.Est,
          pval = output_mr_egger_corr@Pvalue.Est
        ),
        data.table::data.table(
          exp = protein,
          outc = outc_lower,
          pvalthreshold = p_threshold,
          rsqthreshold = r2,
          nsnp = output_mr_egger_corr@SNPs,
          method = "Egger intercept (correlation inc)",
          b = output_mr_egger_corr@Intercept,
          se = output_mr_egger_corr@StdError.Int,
          pval = output_mr_egger_corr@Pvalue.Int
        )
      ), fill = TRUE)
    }
  }
  list(results = results, instruments = dat)
}

run_mr_for_outcome <- function(info) {
  hits <- candidate_hits(info$Outcome, mr_max_proteins)
  if (nrow(hits) == 0) return(list(results = data.table::data.table(), instruments = data.table::data.table()))

  pqtls <- list()
  needed_keys <- character()
  needed_snps <- character()
  for (protein in hits$Protein) {
    pdat <- read_pqtl_cis(protein)
    if (!is.null(pdat)) {
      pdat_sig <- pdat[P < mr_p_threshold]
      if (nrow(pdat_sig) > 0) {
        pqtls[[protein]] <- pdat_sig
        needed_keys <- union(needed_keys, pdat_sig[, paste(CHR, POS, sep = ":")])
        needed_snps <- union(needed_snps, pdat_sig$SNP)
      }
    }
  }
  if (length(pqtls) == 0) return(list(results = data.table::data.table(), instruments = data.table::data.table()))
  message("  MR ", info$outc, ": reading GWAS overlap for ", length(needed_keys), " cis pQTL positions")
  gwas <- read_gwas_subset(info$gwas_path, needed_keys, needed_snps)

  result_rows <- list()
  instrument_rows <- list()
  for (protein in names(pqtls)) {
    dat_u <- make_harmonised(pqtls[[protein]], gwas, info$outc_lower)
    if (nrow(dat_u) == 0) next
    clumped <- plink2_clump(dat_u, p_threshold = mr_p_threshold, r2 = mr_r2, kb = mr_kb)
    dat <- dat_u[SNP %in% clumped]
    if (nrow(dat) == 0) next
    mr <- run_mr_estimates(dat, protein, info$outc_lower, p_threshold = mr_p_threshold, r2 = mr_r2)
    if (nrow(mr$results) > 0) result_rows[[protein]] <- mr$results
    if (nrow(mr$instruments) > 0) {
      instr <- copy(mr$instruments)
      instr[, `:=`(
        exp = protein,
        outc = info$outc_lower,
        pvalthreshold = mr_p_threshold,
        rsqthreshold = mr_r2
      )]
      instrument_rows[[protein]] <- instr
    }
  }
  list(
    results = if (length(result_rows) > 0) data.table::rbindlist(result_rows, fill = TRUE) else data.table::data.table(),
    instruments = if (length(instrument_rows) > 0) data.table::rbindlist(instrument_rows, fill = TRUE) else data.table::data.table()
  )
}

run_mr_sensitivity_grid <- function(primary_pairs) {
  primary_pairs <- data.table::as.data.table(primary_pairs)
  if (nrow(primary_pairs) == 0) {
    return(list(results = data.table::data.table(), instruments = data.table::data.table()))
  }
  primary_pairs <- unique(primary_pairs[, .(exp, outc)])
  p_thresholds <- sort(mr_sensitivity_p_thresholds, decreasing = TRUE)
  r2_thresholds <- sort(mr_sensitivity_r2_thresholds)
  max_p <- max(p_thresholds)
  result_rows <- list()
  instrument_rows <- list()

  for (i in seq_len(nrow(outcome_info))) {
    info <- outcome_info[i]
    pairs <- primary_pairs[outc == info$outc]
    if (nrow(pairs) == 0) next

    pqtls <- list()
    needed_keys <- character()
    needed_snps <- character()
    for (protein in pairs$exp) {
      pdat <- read_pqtl_cis(protein)
      if (!is.null(pdat)) {
        pdat <- pdat[P < max_p]
        if (nrow(pdat) > 0) {
          pqtls[[protein]] <- pdat
          needed_keys <- union(needed_keys, pdat[, paste(CHR, POS, sep = ":")])
          needed_snps <- union(needed_snps, pdat$SNP)
        }
      }
    }
    if (length(pqtls) == 0) next
    message("  sensitivity ", info$outc, ": reading GWAS overlap for ", length(needed_keys), " cis pQTL positions")
    gwas <- read_gwas_subset(info$gwas_path, needed_keys, needed_snps)

    for (protein in names(pqtls)) {
      for (p_threshold in p_thresholds) {
        pdat_sig <- pqtls[[protein]][P < p_threshold]
        if (nrow(pdat_sig) == 0) next
        dat_u <- make_harmonised(pdat_sig, gwas, info$outc_lower)
        if (nrow(dat_u) == 0) next
        for (r2_threshold in r2_thresholds) {
          clumped <- plink2_clump(dat_u, p_threshold = p_threshold, r2 = r2_threshold, kb = mr_kb)
          dat <- dat_u[SNP %in% clumped]
          if (nrow(dat) == 0) next
          mr <- run_mr_estimates(dat, protein, info$outc_lower, p_threshold = p_threshold, r2 = r2_threshold)
          key <- paste(info$outc, protein, p_threshold, r2_threshold, sep = "_")
          if (nrow(mr$results) > 0) result_rows[[key]] <- mr$results
          if (nrow(mr$instruments) > 0) {
            instr <- copy(mr$instruments)
            instr[, `:=`(
              exp = protein,
              outc = info$outc_lower,
              pvalthreshold = p_threshold,
              rsqthreshold = r2_threshold
            )]
            instrument_rows[[key]] <- instr
          }
        }
      }
    }
  }

  list(
    results = if (length(result_rows) > 0) data.table::rbindlist(result_rows, fill = TRUE) else data.table::data.table(),
    instruments = if (length(instrument_rows) > 0) data.table::rbindlist(instrument_rows, fill = TRUE) else data.table::data.table()
  )
}

run_coloc_for_pair <- function(protein, info, gwas = NULL) {
  pdat <- read_pqtl_cis(protein)
  if (is.null(pdat) || nrow(pdat) == 0) return(NULL)
  if (is.null(gwas)) {
    gwas <- read_gwas_subset(info$gwas_path, pdat[, paste(CHR, POS, sep = ":")], pdat$SNP)
  }
  dat <- make_harmonised(pdat, gwas, info$outc_lower)
  dat <- dat[mr_keep %in% TRUE]
  if (nrow(dat) < 10L) {
    return(data.table::data.table(
      exp = protein, outc = info$outc_lower, nsnps = nrow(dat),
      pp_h0 = NA_real_, pp_h1 = NA_real_, pp_h2 = NA_real_, pp_h3 = NA_real_, pp_h4 = NA_real_,
      status = "too_few_snps"
    ))
  }
  ld_obj <- build_ld_matrix(dat)
  dat <- ld_obj$dat
  ld <- ld_obj$ld
  if (nrow(dat) < 10L || nrow(ld) < 10L) {
    return(data.table::data.table(
      exp = protein, outc = info$outc_lower, nsnps = nrow(dat),
      pp_h0 = NA_real_, pp_h1 = NA_real_, pp_h2 = NA_real_, pp_h3 = NA_real_, pp_h4 = NA_real_,
      status = "too_few_ld_snps"
    ))
  }
  dat <- dat[match(rownames(ld), SNP)]
  dat_exp <- list(
    beta = dat$beta.exposure,
    varbeta = dat$se.exposure^2,
    snp = dat$SNP,
    position = dat$pos.exposure,
    type = "quant",
    sdY = 1,
    LD = ld,
    N = unique(dat$samplesize.exposure)[1]
  )
  dat_outc <- list(
    beta = dat$beta.outcome,
    varbeta = dat$se.outcome^2,
    snp = dat$SNP,
    position = dat$pos.outcome,
    type = "cc",
    LD = ld,
    N = info$paper_N,
    s = info$paper_s
  )
  coloc_outc <- tryCatch(coloc::coloc.abf(dataset1 = dat_exp, dataset2 = dat_outc), error = function(e) e)
  if (inherits(coloc_outc, "error")) {
    return(data.table::data.table(
      exp = protein, outc = info$outc_lower, nsnps = nrow(dat),
      pp_h0 = NA_real_, pp_h1 = NA_real_, pp_h2 = NA_real_, pp_h3 = NA_real_, pp_h4 = NA_real_,
      status = paste("error", coloc_outc$message, sep = ": ")
    ))
  }
  data.table::data.table(
    exp = protein,
    outc = info$outc_lower,
    nsnps = unname(coloc_outc$summary["nsnps"]),
    pp_h0 = unname(coloc_outc$summary["PP.H0.abf"]),
    pp_h1 = unname(coloc_outc$summary["PP.H1.abf"]),
    pp_h2 = unname(coloc_outc$summary["PP.H2.abf"]),
    pp_h3 = unname(coloc_outc$summary["PP.H3.abf"]),
    pp_h4 = unname(coloc_outc$summary["PP.H4.abf"]),
    status = "ok"
  )
}

message("Running paper-style primary MR")
mr_runs <- lapply(seq_len(nrow(outcome_info)), function(i) run_mr_for_outcome(outcome_info[i]))
mr_res <- data.table::rbindlist(lapply(mr_runs, `[[`, "results"), fill = TRUE)
mr_instr <- data.table::rbindlist(lapply(mr_runs, `[[`, "instruments"), fill = TRUE)

if (nrow(mr_res) > 0) {
  mr_res[, outc := toupper(outc)]
  obs <- primary[, .(exp = Protein, outc, beta_prim_obs = beta, HR_obs = HR, P_obs = P_Value, FDR_obs = FDR, Bonferroni_obs = Bonferroni)]
  mr_res <- merge(mr_res, obs, by = c("exp", "outc"), all.x = TRUE)
  mr_res[, pval_fdr_local := p.adjust(pval, method = "BH")]
  mr_res[, directional_obs := data.table::fifelse(b * beta_prim_obs > 0, "Consistent", "Not consistent")]
  mr_res[, col := data.table::fifelse(
    pval < 0.05 & b < 0, "Protective",
    data.table::fifelse(pval < 0.05 & b > 0, "Detrimental", "Other")
  )]
  mr_res[, group := data.table::fcase(
    pval < 0.05 & directional_obs == "Consistent", "Nominal MR / consistent with\nobservational findings",
    pval < 0.05 & directional_obs == "Not consistent", "Nominal MR / not consistent\nwith observational findings",
    default = "Nonsignificant"
  )]
  mr_res[, source := paste0(
    "paper-style local cis MR; P<", mr_p_threshold,
    "; PLINK2 1000G EUR clump kb=", mr_kb,
    ", r2=", mr_r2,
    "; correlated IVW via MendelianRandomization"
  )]
}
if (nrow(mr_instr) > 0) {
  required_instrument_columns <- c("exp", "outc", "pvalthreshold", "rsqthreshold")
  missing_instrument_columns <- setdiff(required_instrument_columns, names(mr_instr))
  if (length(missing_instrument_columns) > 0) {
    stop(
      "Internal MR instrument schema error; missing columns: ",
      paste(missing_instrument_columns, collapse = ", "),
      call. = FALSE
    )
  }
  mr_instr[, outc := toupper(outc)]
}
if (nrow(mr_res) == 0) {
  mr_res <- data.table::data.table(
    exp = character(), outc = character(), pvalthreshold = numeric(), rsqthreshold = numeric(),
    nsnp = integer(), method = character(), b = numeric(), se = numeric(), pval = numeric()
  )
}
if (nrow(mr_instr) == 0) {
  mr_instr <- data.table::data.table()
}

data.table::fwrite(mr_res, file.path(mr_dir, "primaryanalyses_loop_corr.csv"))
data.table::fwrite(mr_instr, file.path(mr_dir, "primaryanalyses_loop_corr_instr.csv"))
data.table::fwrite(mr_res, file.path(mr_dir, "local_pqtl_mr_all_outcomes.csv"))
data.table::fwrite(mr_instr, file.path(mr_dir, "local_pqtl_mr_instruments.csv"))
if (!"pval_fdr_local" %in% names(mr_res)) mr_res[, pval_fdr_local := numeric()]
mr_summary <- mr_res[, .(
  n_rows = .N,
  n_primary_methods = sum(method %in% c("Wald ratio", "Inverse variance weighted (correlation inc)"), na.rm = TRUE),
  n_nominal_primary = sum(method %in% c("Wald ratio", "Inverse variance weighted (correlation inc)") & pval < 0.05, na.rm = TRUE),
  n_fdr_primary = sum(method %in% c("Wald ratio", "Inverse variance weighted (correlation inc)") & pval_fdr_local < 0.05, na.rm = TRUE)
), by = outc]
data.table::fwrite(mr_summary, file.path(mr_dir, "local_pqtl_mr_summary.csv"))

message("Preparing MR-significant pairs")
primary_mr <- mr_res[
  method %in% c("Wald ratio", "Inverse variance weighted (correlation inc)") &
    pval < 0.05
]
primary_mr <- primary_mr[order(pval)]
primary_mr <- primary_mr[!duplicated(primary_mr[, .(exp, outc)])]
if (coloc_max_proteins > 0L) {
  primary_mr <- primary_mr[, head(.SD, coloc_max_proteins), by = outc]
}

message("Running paper-style MR sensitivity grid")
if (run_mr_sensitivity && nrow(primary_mr) > 0) {
  sens_runs <- run_mr_sensitivity_grid(primary_mr)
  sens_res <- sens_runs$results
  sens_instr <- sens_runs$instruments
} else {
  sens_res <- data.table::data.table()
  sens_instr <- data.table::data.table()
}
if (nrow(sens_res) > 0) {
  sens_res[, outc := toupper(outc)]
  obs_sens <- primary[, .(exp = Protein, outc, beta_prim_obs = beta, HR_obs = HR, P_obs = P_Value, FDR_obs = FDR, Bonferroni_obs = Bonferroni)]
  sens_res <- merge(sens_res, obs_sens, by = c("exp", "outc"), all.x = TRUE)
  sens_res[, pval_fdr_sensitivity := p.adjust(pval, method = "BH")]
  sens_res[, directional_obs := data.table::fifelse(b * beta_prim_obs > 0, "Consistent", "Not consistent")]
  sens_res[, source := paste0(
    "paper-style local cis MR sensitivity grid; pQTL P in ",
    paste(mr_sensitivity_p_thresholds, collapse = "/"),
    "; PLINK2 1000G EUR clump kb=", mr_kb,
    ", r2 in ", paste(mr_sensitivity_r2_thresholds, collapse = "/"),
    "; correlated IVW via MendelianRandomization"
  )]
}
if (nrow(sens_instr) > 0) {
  sens_instr[, outc := toupper(outc)]
}
if (nrow(sens_res) == 0) {
  sens_res <- data.table::data.table(
    exp = character(), outc = character(), pvalthreshold = numeric(), rsqthreshold = numeric(),
    nsnp = integer(), method = character(), b = numeric(), se = numeric(), pval = numeric()
  )
}
if (nrow(sens_instr) == 0) sens_instr <- data.table::data.table()
data.table::fwrite(sens_res, file.path(mr_dir, "sensitivityanalyses_loop_corr.csv"))
data.table::fwrite(sens_instr, file.path(mr_dir, "sensitivityanalyses_loop_corr_instr.csv"))
data.table::fwrite(sens_res, file.path(mr_dir, "local_pqtl_mr_sensitivity_grid.csv"))
if (nrow(sens_res) > 0) {
  sens_summary <- sens_res[, .(
    n_rows = .N,
    n_primary_methods = sum(method %in% c("Wald ratio", "Inverse variance weighted (correlation inc)"), na.rm = TRUE),
    n_nominal_primary = sum(method %in% c("Wald ratio", "Inverse variance weighted (correlation inc)") & pval < 0.05, na.rm = TRUE)
  ), by = .(outc, pvalthreshold, rsqthreshold)]
} else {
  sens_summary <- data.table::data.table()
}
data.table::fwrite(sens_summary, file.path(mr_dir, "local_pqtl_mr_sensitivity_summary.csv"))

message("Running paper-style coloc on MR-significant pairs")
coloc_rows <- list()
if (nrow(primary_mr) > 0) {
  for (i in seq_len(nrow(outcome_info))) {
    info <- outcome_info[i]
    pairs <- primary_mr[outc == info$outc]
    if (nrow(pairs) == 0) next
    needed_keys <- character()
    needed_snps <- character()
    pqtls_full <- list()
    for (protein in pairs$exp) {
      pdat <- read_pqtl_cis(protein)
      if (!is.null(pdat)) {
        pqtls_full[[protein]] <- pdat
        needed_keys <- union(needed_keys, pdat[, paste(CHR, POS, sep = ":")])
        needed_snps <- union(needed_snps, pdat$SNP)
      }
    }
    if (length(pqtls_full) == 0) next
    message("  coloc ", info$outc, ": reading GWAS overlap for ", length(needed_keys), " cis positions")
    gwas <- read_gwas_subset(info$gwas_path, needed_keys, needed_snps)
    for (protein in names(pqtls_full)) {
      coloc_rows[[paste(info$outc, protein, sep = "_")]] <- run_coloc_for_pair(protein, info, gwas)
    }
  }
}
coloc_res <- data.table::rbindlist(coloc_rows, fill = TRUE)
if (nrow(coloc_res) == 0) {
  coloc_res <- data.table::data.table(
    exp = character(), outc = character(), nsnps = integer(),
    pp_h0 = numeric(), pp_h1 = numeric(), pp_h2 = numeric(), pp_h3 = numeric(), pp_h4 = numeric(),
    status = character()
  )
}
if (nrow(coloc_res) > 0) {
  coloc_res[, outc := toupper(outc)]
  obs <- primary[, .(exp = Protein, outc, obs_beta = beta, obs_p = P_Value, obs_hr = HR)]
  coloc_res <- merge(coloc_res, obs, by = c("exp", "outc"), all.x = TRUE)
  coloc_res[, source := "paper-style local coloc; MR-significant pairs; coloc.abf with PLINK2 1000G EUR LD matrix"]
}
data.table::fwrite(coloc_res, file.path(coloc_dir, "coloc_loop.csv"))
data.table::fwrite(coloc_res, file.path(coloc_dir, "local_pqtl_coloc_all_outcomes.csv"))

mr_cols <- c(
  "Nonsignificant" = "grey90",
  "Nominal MR / consistent with\nobservational findings" = "#E8B900",
  "Nominal MR / not consistent\nwith observational findings" = "red3"
)
make_mr_panel <- function(outc_code) {
  d <- copy(mr_res[outc == outc_code & method %in% c("Wald ratio", "Inverse variance weighted (correlation inc)")])
  title <- outcome_info$title[match(outc_code, outcome_info$outc)]
  if (nrow(d) == 0) {
    return(ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No MR result", color = "grey40") +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void() + ggplot2::ggtitle(title))
  }
  d[, exp_label := data.table::fifelse(exp == "NTPROBNP", "NT-proBNP", exp)]
  label_dt <- d[pval < 0.05][order(pval)][seq_len(min(.N, 12))]
  ggplot2::ggplot(d, ggplot2::aes(b, -log10(pval), color = group)) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_color_manual(values = mr_cols, name = "", drop = FALSE) +
    ggrepel::geom_text_repel(
      data = label_dt,
      ggplot2::aes(label = exp_label),
      size = 2.7,
      force = 20,
      max.overlaps = 80,
      show.legend = FALSE
    ) +
    ggplot2::labs(title = title, x = "log(OR)", y = expression(-log[10](P))) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(hjust = 0.5, size = 11))
}

mr_order <- c("CAD", "HF", "AF", "AS")
mr_panels <- lapply(mr_order, make_mr_panel)
names(mr_panels) <- mr_order
mr_plot <- (mr_panels[["CAD"]] | mr_panels[["HF"]]) / (mr_panels[["AF"]] | mr_panels[["AS"]]) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(fig_dir, "original_style_local_pqtl_mr_all.png"), mr_plot, width = 11, height = 8.4, dpi = 320, bg = "white")
ggplot2::ggsave(file.path(mr_dir, "local_pqtl_mr_all.png"), mr_plot, width = 11, height = 8.4, dpi = 320, bg = "white")

if ("status" %in% names(coloc_res) && nrow(coloc_res[status == "ok"]) > 0) {
  pdat <- coloc_res[status == "ok"][order(outc, -pp_h4)]
  pdat[, exp_label := data.table::fifelse(exp == "NTPROBNP", "NT-proBNP", exp)]
  pdat[, exp_plot := factor(exp_label, levels = rev(unique(exp_label)))]
  pdat[, outc_title := factor(outcome_info$title[match(outc, outcome_info$outc)], levels = outcome_info$title[match(mr_order, outcome_info$outc)])]
  p_coloc <- ggplot2::ggplot(pdat, ggplot2::aes(exp_plot, pp_h4)) +
    ggplot2::geom_col(color = "black", fill = "#3C5488", linewidth = 0.2) +
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "red3", linewidth = 0.35) +
    ggplot2::annotate("text", x = Inf, y = 0.805, label = "PP.H4 = 0.8", hjust = 1.05, vjust = -0.2, size = 2.7, color = "red3") +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ outc_title, scales = "free_y", ncol = 2) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0.01, 0.06))) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::labs(title = "Local coloc by outcome", x = "", y = "PP.H4")
  ggplot2::ggsave(file.path(fig_dir, "original_style_local_pqtl_coloc_pph4_all.png"), p_coloc, width = 11, height = 9.0, dpi = 320, bg = "white")
  ggplot2::ggsave(file.path(coloc_dir, "local_pqtl_coloc_pph4_all.png"), p_coloc, width = 11, height = 9.0, dpi = 320, bg = "white")
}

figure_manifest <- data.table::data.table(
  file = list.files(fig_dir, pattern = "\\.(png|csv)$", full.names = FALSE),
  path = list.files(fig_dir, pattern = "\\.(png|csv)$", full.names = TRUE)
)
data.table::fwrite(figure_manifest, file.path(fig_dir, "figure_manifest.csv"))

message("MR results: ", file.path(mr_dir, "primaryanalyses_loop_corr.csv"))
message("MR instruments: ", file.path(mr_dir, "primaryanalyses_loop_corr_instr.csv"))
message("MR sensitivity results: ", file.path(mr_dir, "sensitivityanalyses_loop_corr.csv"))
message("Coloc results: ", file.path(coloc_dir, "coloc_loop.csv"))
print(mr_summary)
