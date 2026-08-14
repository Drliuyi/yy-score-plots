#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0L) {
    return(normalizePath(sub("^--file=", "", hit[[1]]), winslash = "/", mustWork = FALSE))
  }
  normalizePath(file.path("f", "cli.R"), winslash = "/", mustWork = FALSE)
}

project_dir <- Sys.getenv("PRADEEP_PROJECT_DIR", unset = "")
if (!nzchar(project_dir)) project_dir <- dirname(dirname(script_path()))
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)
f_dir <- file.path(project_dir, "f")

default_dir0 <- function() {
  explicit <- Sys.getenv("DIR0", unset = Sys.getenv("PRADEEP_DIR0", unset = ""))
  if (nzchar(explicit)) return(explicit)
  if (.Platform$OS.type == "windows") "D:/" else "/mnt/d"
}

clean_path <- function(x) {
  x <- gsub("\\\\", "/", trimws(as.character(x)))
  if (x %in% c("/", "")) return(x)
  if (grepl("^[A-Za-z]:/$", x)) return(x)
  sub("/+$", "", x)
}

as_flag <- function(x) {
  if (is.null(x) || length(x) == 0L) return(FALSE)
  if (is.logical(x)) return(isTRUE(x))
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

show_help <- function() {
  cat(paste0(
    "Pradeep/Schuermans UKB-PPP reproduction\n",
    "\n",
    "Usage\n",
    "  bash pradeep.sh --h\n",
    "  bash pradeep.sh --step <stage> [paths and options]\n",
    "\n",
    "Stages\n",
    "  preflight   Validate paths, inputs, R packages, and runtime only\n",
    "  1           Build the cohort, four outcome controls, protein QC, and analysis base\n",
    "  2           Run protein Cox associations for the selected outcome(s)\n",
    "  3           Run sex-stratified analyses and interaction tests\n",
    "  4           Fit and evaluate Clinical, Protein, and Combined LASSO models\n",
    "  5           Generate paper-style main, ROC, calibration, and feature figures\n",
    "  6           Run GO enrichment and generate paper-style enrichment figures\n",
    "  7           Standardize downloaded FinnGen files for MR/coloc input\n",
    "  8           Run local cis-pQTL MR/coloc and generate result figures\n",
    "  core        Run Stages 1-6 in order\n",
    "  figures     Run Stages 5-6 and rebuild the complete figure inventory\n",
    "  mr          Run Stage 8; run Stage 7 first only if GWAS files are not standardized\n",
    "  all         Run Stages 1-8 in order\n",
    "  status      Read-only result and figure completeness check\n",
    "  combinations  --step 2,3,4 or --step core,mr\n",
    "\n",
    "Validated environment: WSL Ubuntu 22.04.5, Bash 5.1.16, and R 4.3.2.\n",
    "Stage 8 uses PLINK2 2.0.0-a.6.9LM. See README.md for all R package versions.\n",
    "\n",
    "Huang-lab path defaults\n",
    "  DIR0          D:/ (or /mnt/d in WSL)\n",
    "  PHEDIR        <DIR0>/data/ukb/phe\n",
    "  SCRIPT_ROOT   <DIR0>/scripts\n",
    "  HELPER_DIR    <DIR0>/scripts/0f\n",
    "  OUTPUT_ROOT   <DIR0>/analysis/ukb\n",
    "  RESULT_DIR    <OUTPUT_ROOT>/pradeep\n",
    "  PROTEIN_BED   <DIR0>/data.BIG/gwas/ppp/ppp_3k_b38.bed\n",
    "\n",
    "Common options\n",
    "  --profile huang|winpc      Path preset; UKB_data paths are detected as winpc\n",
    "  --panel 1.5k|3k            Paper 1.5k panel or the existing full 3k panel\n",
    "  --dir0 PATH                 Logical D-drive root\n",
    "  --phe-dir PATH              Phenotype directory\n",
    "  --out-dir PATH              Final result directory; highest priority\n",
    "  --output-root PATH          Parent output directory\n",
    "  --analysis-project NAME     Result subdirectory name; default: pradeep\n",
    "  --outcome NAME              Outcome to fit: cad, afib, hfail, ao_sten, or all\n",
    "  --all-rds FILE              all.rds\n",
    "  --prot-rds FILE             prot.rds\n",
    "  --raw-protein FILE          Unimputed prot.tab.gz; required for the 1.5k panel\n",
    "  --protein-map FILE          1.5k assay mapping\n",
    "  --kinship-file FILE         Pairwise kinship file\n",
    "  --bed-file FILE             PPP protein genomic-coordinate BED\n",
    "  --follow-end YYYY-MM-DD     Follow-up cutoff; default: 2020-03-31\n",
    "  --workers N                 Requested parallel workers; default: 6\n",
    "  --resume                    Skip stages with accepted existing outputs\n",
    "  --allow-kinship-fallback    Allow fallback to UKB field 22011 in all.rds\n",
    "  --dry-run                   Validate and print the plan without computation\n",
    "\n",
    "MR/coloc options\n",
    "  --gwas-dir PATH             Root of standardized GWAS files\n",
    "  --pqtl-dir PATH             Protein cis-pQTL directory\n",
    "  --ld-bfile PREFIX           PLINK1 LD-reference prefix; may point to Zspace\n",
    "  --ld-pfile-template PATH    PLINK2 template containing {CHR}\n",
    "  --plink2 PATH               PLINK2 executable in WSL/Linux\n",
    "  --finngen-download-dir PATH FinnGen download directory\n",
    "\n",
    "Example: Huang-lab layout, CAD-only reproduction\n",
    "  bash /mnt/d/scripts/ukb/pradeep/pradeep.sh --step 1,2,3,4 --outcome cad --panel 1.5k --workers 12\n",
    "\n",
    "Example: Huang-lab layout, existing 3k all.rds/prot.rds\n",
    "  bash /mnt/d/scripts/ukb/pradeep/pradeep.sh --step core --panel 3k --workers 12\n",
    "\n",
    "Example: from an open WinPC WSL terminal (use backslash for Bash continuation)\n",
    "  cd /mnt/d/scripts/ukb/pradeep\n",
    "  ./pradeep.sh \\\n",
    "    --step core \\\n",
    "    --profile winpc \\\n",
    "    --panel 1.5k \\\n",
    "    --allow-kinship-fallback \\\n",
    "    --workers 12\n",
    "\n",
    "Example: call WSL from WinPC PowerShell on one line\n",
    "  wsl bash /mnt/d/scripts/ukb/pradeep/pradeep.sh --step core --profile winpc --panel 1.5k --allow-kinship-fallback --workers 12\n",
    "\n",
    "Use --kinship-file FILE when a genuine pairwise kinship file is available.\n",
    "Use --allow-kinship-fallback only when UKB field 22011 fallback is accepted.\n",
    "Path mismatches are reported before computation with the missing item,\n",
    "default lookup location, and required override; data are never substituted silently.\n"
  ))
}

option_aliases <- c(
  h = "help", help = "help", step = "step", steps = "step",
  profile = "profile", panel = "panel", `panel-mode` = "panel",
  dir0 = "dir0", `phe-dir` = "phe_dir", phedir = "phe_dir",
  `script-root` = "script_root", `helper-dir` = "helper_dir",
  `output-root` = "output_root", `out-dir` = "out_dir",
  `analysis-project` = "analysis_project", `all-rds` = "all_rds",
  outcome = "outcome", outcomes = "outcome",
  `prot-rds` = "prot_rds", `raw-protein` = "raw_protein",
  `pheno-tsv` = "pheno_tsv", `protein-map` = "protein_map",
  `kinship-file` = "kinship_file", `bed-file` = "bed_file",
  `follow-end` = "follow_end", workers = "workers", resume = "resume",
  `dry-run` = "dry_run", `allow-kinship-fallback` = "allow_kinship_fallback",
  `gwas-dir` = "gwas_dir", `pqtl-dir` = "pqtl_dir",
  `ld-bfile` = "ld_bfile", `ld-pfile-template` = "ld_pfile_template",
  plink2 = "plink2", `finngen-download-dir` = "finngen_download_dir"
)
flag_options <- c("help", "resume", "dry_run", "allow_kinship_fallback")

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "-")) stop("Unexpected positional argument: ", token, call. = FALSE)
    token <- sub("^-+", "", token)
    if (grepl("=", token, fixed = TRUE)) {
      key <- sub("=.*$", "", token)
      value <- sub("^[^=]*=", "", token)
    } else {
      key <- token
      canonical <- unname(option_aliases[key])
      if (length(canonical) == 0L || is.na(canonical)) stop("Unknown option --", key, ". Run --h.", call. = FALSE)
      if (canonical %in% flag_options) {
        value <- TRUE
      } else {
        if (i == length(args) || startsWith(args[[i + 1L]], "-")) {
          stop("Option --", key, " requires a value.", call. = FALSE)
        }
        i <- i + 1L
        value <- args[[i]]
      }
    }
    canonical <- unname(option_aliases[key])
    if (length(canonical) == 0L || is.na(canonical)) stop("Unknown option --", key, ". Run --h.", call. = FALSE)
    out[[canonical]] <- value
    i <- i + 1L
  }
  out
}

cli <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(cli$help) || length(commandArgs(trailingOnly = TRUE)) == 0L) {
  show_help()
  quit(status = 0L)
}

auto_profile <- if (grepl("^(D:|/mnt/d)/UKB_data(/|$)", project_dir, ignore.case = TRUE)) "winpc" else "huang"
profile <- tolower(if (!is.null(cli$profile)) as.character(cli$profile) else auto_profile)
if (!profile %in% c("huang", "winpc")) stop("--profile must be huang or winpc.", call. = FALSE)
panel_input <- tolower(if (!is.null(cli$panel)) as.character(cli$panel) else "1.5k")
panel_compact <- gsub("[^a-z0-9]", "", panel_input)
panel_mode <- if (panel_compact %in% c("15k", "1536")) {
  "1.5k"
} else if (panel_compact %in% c("3k", "3072", "full", "fullpanel")) {
  "3k"
} else {
  stop("--panel must be 1.5k or 3k.", call. = FALSE)
}
profile_dir0 <- if (profile == "winpc") {
  if (.Platform$OS.type == "windows") "D:/UKB_data" else "/mnt/d/UKB_data"
} else {
  default_dir0()
}
dir0 <- clean_path(if (!is.null(cli$dir0)) cli$dir0 else profile_dir0)
default_phe_dir <- if (profile == "winpc") file.path(dir0, "phe") else file.path(dir0, "data", "ukb", "phe")
phe_dir <- clean_path(if (!is.null(cli$phe_dir)) cli$phe_dir else Sys.getenv("PHEDIR", unset = default_phe_dir))
script_root <- clean_path(if (!is.null(cli$script_root)) cli$script_root else Sys.getenv("SCRIPT_ROOT", unset = file.path(dir0, "scripts")))
helper_dir <- clean_path(if (!is.null(cli$helper_dir)) cli$helper_dir else Sys.getenv("HELPER_DIR", unset = file.path(script_root, "0f")))
default_output_root <- if (profile == "winpc") file.path(dir0, "analysis") else file.path(dir0, "analysis", "ukb")
output_root <- clean_path(if (!is.null(cli$output_root)) cli$output_root else Sys.getenv("UKB_OUT", unset = default_output_root))
analysis_project <- if (!is.null(cli$analysis_project)) cli$analysis_project else "pradeep"
if (!grepl("^[A-Za-z0-9_.-]+$", analysis_project)) {
  stop("--analysis-project must be a directory name, not a path.", call. = FALSE)
}
outcome_input <- tolower(trimws(if (!is.null(cli$outcome)) as.character(cli$outcome) else "all"))
outcome_subset <- if (outcome_input %in% c("all", "*")) {
  c("cad", "afib", "hfail", "ao_sten")
} else {
  unique(trimws(unlist(strsplit(outcome_input, "[,;+]"))))
}
valid_outcomes <- c("cad", "afib", "hfail", "ao_sten")
unknown_outcomes <- setdiff(outcome_subset, valid_outcomes)
if (length(outcome_subset) == 0L || length(unknown_outcomes) > 0L) {
  stop("--outcome must be cad, afib, hfail, ao_sten, or all.", call. = FALSE)
}
analysis_dir <- clean_path(if (!is.null(cli$out_dir)) cli$out_dir else file.path(output_root, analysis_project))
source_all_rds <- clean_path(if (!is.null(cli$all_rds)) cli$all_rds else file.path(phe_dir, "Rdata", "all.rds"))
source_prot_rds <- clean_path(if (!is.null(cli$prot_rds)) cli$prot_rds else file.path(phe_dir, "Rdata", "prot.rds"))
default_ppp_dir <- if (profile == "winpc") file.path(dir0, "ppp") else file.path(dir0, "data.BIG", "gwas", "ppp")
raw_protein <- clean_path(if (!is.null(cli$raw_protein)) cli$raw_protein else if (profile == "winpc") file.path(dir0, "prot.tab.gz") else file.path(phe_dir, "rap", "raw", "prot.tab.gz"))
pheno_tsv <- clean_path(if (!is.null(cli$pheno_tsv)) cli$pheno_tsv else if (profile == "winpc") file.path(dir0, "pheno.tsv.gz") else file.path(phe_dir, "pheno.tsv.gz"))
protein_map <- clean_path(if (!is.null(cli$protein_map)) cli$protein_map else file.path(default_ppp_dir, "map.raw", "olink_protein_map_1.5k_v1.tsv"))
panel_input_dir <- clean_path(file.path(analysis_dir, "inputs", if (panel_mode == "1.5k") "pradeep_15k" else "existing_3k"))
all_rds <- if (panel_mode == "1.5k") file.path(panel_input_dir, "all.rds") else source_all_rds
prot_rds <- if (panel_mode == "1.5k") file.path(panel_input_dir, "prot.rds") else source_prot_rds

detect_existing_panel <- function(project_dir) {
  manifest_file <- file.path(project_dir, "audit", "command_manifest.csv")
  if (file.exists(manifest_file)) {
    manifest <- tryCatch(utils::read.csv(manifest_file, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(manifest) && all(c("key", "value") %in% names(manifest))) {
      explicit <- manifest$value[manifest$key == "PRADEEP_PANEL_MODE"]
      if (length(explicit) > 0L && explicit[[1L]] %in% c("1.5k", "3k")) return(explicit[[1L]])
      legacy <- manifest$value[manifest$key == "UKBPPP_ORIGINAL_15K_PANEL"]
      if (length(legacy) > 0L) return(if (legacy[[1L]] %in% c("1", "TRUE", "true")) "1.5k" else "3k")
    }
  }
  proteins_file <- file.path(project_dir, "audit", "protein_columns_used.csv")
  if (file.exists(proteins_file)) {
    protein_n <- tryCatch(nrow(utils::read.csv(proteins_file, stringsAsFactors = FALSE)), error = function(e) NA_integer_)
    if (is.finite(protein_n) && protein_n > 0L) return(if (protein_n <= 1600L) "1.5k" else "3k")
  }
  NA_character_
}
existing_panel <- detect_existing_panel(analysis_dir)
bed_file <- clean_path(if (!is.null(cli$bed_file)) cli$bed_file else file.path(default_ppp_dir, if (profile == "winpc") "ppp.b38.bed" else "ppp_3k_b38.bed"))
kinship_candidates <- c(
  if (!is.null(cli$kinship_file)) clean_path(cli$kinship_file) else character(),
  file.path(phe_dir, "common", "ukb7089_rel_s488363.dat"),
  file.path(phe_dir, "ukb7089_rel_s488363.dat")
)
kinship_hit <- kinship_candidates[file.exists(kinship_candidates)]
kinship_file <- if (length(kinship_hit) > 0L) clean_path(kinship_hit[[1]]) else if (!is.null(cli$kinship_file)) clean_path(cli$kinship_file) else ""
follow_end <- if (!is.null(cli$follow_end)) as.character(cli$follow_end) else "2020-03-31"
if (is.na(as.Date(follow_end))) stop("--follow-end must be YYYY-MM-DD.", call. = FALSE)
workers <- suppressWarnings(as.integer(if (!is.null(cli$workers)) cli$workers else 6L))
if (is.na(workers) || workers < 1L) stop("--workers must be a positive integer.", call. = FALSE)
resume <- as_flag(cli$resume)
dry_run <- as_flag(cli$dry_run)
allow_kinship_fallback <- as_flag(cli$allow_kinship_fallback)
gwas_dir <- clean_path(if (!is.null(cli$gwas_dir)) cli$gwas_dir else if (profile == "winpc") file.path(dir0, "gwas") else file.path(dir0, "data", "gwas"))
pqtl_dir <- clean_path(if (!is.null(cli$pqtl_dir)) cli$pqtl_dir else file.path(default_ppp_dir, "clean"))
ld_bfile <- clean_path(if (!is.null(cli$ld_bfile)) cli$ld_bfile else "")
ld_pfile_template <- clean_path(if (!is.null(cli$ld_pfile_template)) cli$ld_pfile_template else Sys.getenv("PRADEEP_LD_PFILE_TEMPLATE", unset = ""))
plink2 <- if (!is.null(cli$plink2)) {
  as.character(cli$plink2)
} else {
  Sys.getenv("UKBPPP_PLINK2_WSL", unset = unname(Sys.which("plink2")))
}
finngen_download_dir <- clean_path(if (!is.null(cli$finngen_download_dir)) cli$finngen_download_dir else file.path(gwas_dir, "downloads"))

step_aliases <- list(
  `1` = "1", `2` = "2", `3` = "3", `4` = "4", `5` = "5", `6` = "6", `7` = "7", `8` = "8",
  step01 = "1", step02 = "2", step03 = "3", step04 = "4", step05 = "5", step06 = "6", step07 = "7", step08 = "8",
  associations = c("2", "3"), models = "4", prediction = "4", figures = c("5", "6"), enrichment = "6",
  gwas = "7", mr = "8", core = as.character(1:6), all = as.character(1:8), preflight = character(), status = character()
)
requested_step <- tolower(if (!is.null(cli$step)) cli$step else "preflight")
requested_tokens <- trimws(unlist(strsplit(requested_step, "[,;+]")))
unknown_steps <- requested_tokens[!requested_tokens %in% names(step_aliases)]
if (length(unknown_steps) > 0L) stop("Unknown step: ", paste(unknown_steps, collapse = ", "), ". Run --h.", call. = FALSE)
selected_steps <- unique(unlist(step_aliases[requested_tokens], use.names = FALSE))
is_preflight <- identical(requested_tokens, "preflight")
is_status <- identical(requested_tokens, "status")

step_defs <- list(
  `1` = list(label = "Cohort and analysis base", script = "01_build_analysis_base.R", output = "ukbppp_cardiac_analysis_base.rds"),
  `2` = list(label = "Protein Cox associations", script = "02_primary_association_cox.R", output = "primary_association_cox.csv"),
  `3` = list(label = "Sex-stratified and interaction analyses", script = "03_sex_stratified_interaction.R", output = "sex_interaction_cox.csv"),
  `4` = list(label = "LASSO prediction", script = "04_lasso_risk_score.R", output = "lasso_delong_comparisons.csv"),
  `5` = list(label = "Main and prediction figures", script = "05_make_figures.R", output = "figures/paper_figure_qc.csv"),
  `6` = list(label = "GO enrichment and figures", script = "06_go_enrichment.R", output = "enrichment/go_enrichment_all_terms.csv"),
  `7` = list(label = "FinnGen GWAS standardization", script = "07_prepare_finngen_gwas.R", output = "finngen_r12_gwas_prepare_status.csv"),
  `8` = list(label = "cis-pQTL MR and coloc", script = "08_local_pqtl_mr_coloc.R", output = "mr/local_pqtl_mr_all_outcomes.csv")
)

output_dir <- file.path(analysis_dir, "outputs")
audit_dir <- file.path(analysis_dir, "audit")
base_file <- file.path(output_dir, "ukbppp_cardiac_analysis_base.rds")

path_rows <- data.frame(
  item = c(
    "PROJECT_DIR", "DIR0", "PHEDIR", "SCRIPT_ROOT", "HELPER_DIR", "ANALYSIS_DIR",
    "PANEL_MODE", "EXISTING_PANEL_LOCK", "SOURCE_ALL_RDS", "SOURCE_PROT_RDS", "RAW_PROTEIN",
    "EFFECTIVE_ALL_RDS", "EFFECTIVE_PROT_RDS", "PROTEIN_MAP", "KINSHIP",
    "PPP_BED", "GWAS_DIR", "PQTL_DIR", "LD_BFILE", "LD_PFILE_TEMPLATE"
  ),
  path = c(
    project_dir, dir0, phe_dir, script_root, helper_dir, analysis_dir,
    panel_mode, if (is.na(existing_panel)) "<not set>" else existing_panel,
    source_all_rds, source_prot_rds, raw_protein, all_rds, prot_rds,
    protein_map, if (nzchar(kinship_file)) kinship_file else "<not found>",
    bed_file, gwas_dir, pqtl_dir, if (nzchar(ld_bfile)) ld_bfile else "<not set>",
    if (nzchar(ld_pfile_template)) ld_pfile_template else "<not set>"
  ),
  stringsAsFactors = FALSE
)
cat("\nResolved paths\n")
print(path_rows, row.names = FALSE, right = FALSE)
cat("\nProfile=", profile, "\n", sep = "")
cat("Panel=", panel_mode, "\n", sep = "")
cat("Fitted outcomes=", paste(outcome_subset, collapse = ","), "\n", sep = "")
cat("Requested=", requested_step, "\n", sep = "")
cat("Selected steps=", if (length(selected_steps)) paste(selected_steps, collapse = ",") else "none", "\n", sep = "")
cat("Follow end=", follow_end, "; workers=", workers, "; resume=", resume, "\n", sep = "")

missing <- character()
if (!is.na(existing_panel) && !identical(existing_panel, panel_mode)) {
  missing <- c(missing, paste0(
    "- panel selection conflicts with this existing project\n",
    "  project: ", analysis_dir, "\n",
    "  existing panel: ", existing_panel, "\n",
    "  requested panel: ", panel_mode, "\n",
    "  use the panel already selected for this project, or explicitly choose another --analysis-project"
  ))
}
add_missing <- function(label, path, option) {
  missing <<- c(missing, paste0("- ", label, "\n  expected: ", path, "\n  set with: ", option))
}
need_step <- function(x) x %in% selected_steps
if (is_preflight || need_step("1")) {
  if (!file.exists(source_all_rds)) add_missing("source all.rds", source_all_rds, "--all-rds FILE or --phe-dir DIR")
  if (panel_mode == "1.5k") {
    if (!file.exists(raw_protein)) add_missing("un-imputed raw protein table", raw_protein, "--raw-protein FILE")
    if (!file.exists(protein_map)) add_missing("1.5k protein mapping", protein_map, "--protein-map FILE")
  } else if (!file.exists(source_prot_rds)) {
    add_missing("existing 3k prot.rds", source_prot_rds, "--prot-rds FILE or --phe-dir DIR")
  }
  if (!nzchar(kinship_file) && !allow_kinship_fallback) {
    add_missing("pairwise kinship", file.path(phe_dir, "common", "ukb7089_rel_s488363.dat"), "--kinship-file FILE; or explicitly --allow-kinship-fallback")
  }
}
if (any(c("2", "3", "4", "5") %in% selected_steps) && !need_step("1") && !file.exists(base_file)) {
  add_missing("Step 1 analysis base", base_file, "run --step 1 first")
}
if ((is_preflight || need_step("5")) && !file.exists(bed_file)) add_missing("PPP genomic BED", bed_file, "--bed-file FILE")
if (need_step("5") && !need_step("2") && !file.exists(file.path(output_dir, "primary_association_cox.csv"))) add_missing("Step 2 result", file.path(output_dir, "primary_association_cox.csv"), "run --step 2 first")
if (need_step("5") && !need_step("3") && !file.exists(file.path(output_dir, "sex_interaction_cox.csv"))) add_missing("Step 3 result", file.path(output_dir, "sex_interaction_cox.csv"), "run --step 3 first")
if (need_step("5") && !need_step("4") && !file.exists(file.path(output_dir, "lasso_risk_score_auc.csv"))) add_missing("Step 4 result", file.path(output_dir, "lasso_risk_score_auc.csv"), "run --step 4 first")
if (need_step("6") && !need_step("2") && !file.exists(file.path(output_dir, "primary_association_cox.csv"))) add_missing("Step 2 result", file.path(output_dir, "primary_association_cox.csv"), "run --step 2 first")
if (need_step("8")) {
  if (!need_step("2") && !file.exists(file.path(output_dir, "primary_association_cox.csv"))) add_missing("Step 2 result", file.path(output_dir, "primary_association_cox.csv"), "run --step 2 first")
  if (!dir.exists(pqtl_dir)) add_missing("cis-pQTL directory", pqtl_dir, "--pqtl-dir DIR")
  gwas_expected <- file.path(gwas_dir, c("cad_eur", "hfail_eur", "afib_eur", "ao_sten_eur"), c("finngen_R12_I9_CHD.standard.tsv.gz", "finngen_R12_I9_HEARTFAIL.standard.tsv.gz", "finngen_R12_I9_AF.standard.tsv.gz", "finngen_R12_I9_CAVS_OPERATED.standard.tsv.gz"))
  if (!need_step("7")) {
    for (g in gwas_expected[!file.exists(gwas_expected)]) add_missing("standardized FinnGen GWAS", g, "run --step 7 or set --gwas-dir DIR")
  }
  if (!nzchar(ld_bfile) && !nzchar(ld_pfile_template)) add_missing("LD reference", "not set", "--ld-bfile PREFIX or --ld-pfile-template PATH_WITH_{CHR}; genotype may remain on Zspace")
  if (!nzchar(plink2)) add_missing("PLINK2 executable", "not found", "--plink2 PATH or UKBPPP_PLINK2_WSL")
  if (nzchar(ld_bfile)) {
    for (ext in c(".bed", ".bim", ".fam")) if (!file.exists(paste0(ld_bfile, ext))) add_missing("LD bfile component", paste0(ld_bfile, ext), "--ld-bfile PREFIX")
  }
}

package_map <- list(
  `1` = c("data.table", "dplyr", "survival", "impute", "glmnet", "pROC"),
  `2` = c("data.table", "survival"), `3` = c("data.table", "survival"),
  `4` = c("data.table", "glmnet", "pROC"),
  `5` = c("data.table", "ggplot2", "ggrepel", "patchwork", "survival"),
  `6` = c("data.table", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "ggplot2", "ggtext", "ggpubr", "gtools", "openxlsx"),
  `7` = "data.table",
  `8` = c("data.table", "dplyr", "tidyr", "ggplot2", "ggrepel", "patchwork", "TwoSampleMR", "MendelianRandomization", "coloc")
)
common_packages <- c("data.table", "dplyr", "survival", "impute", "glmnet", "pROC")
packages_to_check <- unique(c(
  if (length(selected_steps) > 0L || is_preflight) common_packages else character(),
  unlist(package_map[if (is_preflight) as.character(1:6) else selected_steps], use.names = FALSE)
))
missing_packages <- packages_to_check[!vapply(packages_to_check, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  missing <- c(missing, paste0("- missing R packages: ", paste(missing_packages, collapse = ", "), "\n  install these packages before running the selected step"))
}

if (length(missing) > 0L) {
  cat("\nPATH_OR_DEPENDENCY_GATE_FAILED\n", paste(missing, collapse = "\n"), "\n", sep = "")
  cat("\nWinPC legacy example: --dir0 D:/UKB_data --phe-dir D:/UKB_data/phe --output-root D:/UKB_data/analysis\n")
  quit(status = 2L)
}

if (is_preflight) {
  cat("\nPREFLIGHT_PASSED\n")
  quit(status = 0L)
}

if (is_status) {
  output_dir_status <- file.path(analysis_dir, "outputs")
  output_files <- if (dir.exists(output_dir_status)) {
    listed <- list.files(output_dir_status, recursive = TRUE, full.names = TRUE)
    listed[file.info(listed)$isdir %in% FALSE]
  } else character()
  figure_files <- output_files[grepl("[.](png|pdf|svg|tif|tiff)$", output_files, ignore.case = TRUE)]
  run_status_file <- file.path(analysis_dir, "audit", "run_status.csv")
  figure_qc_file <- file.path(analysis_dir, "audit", "figure_completeness.csv")
  cat("\nRead-only result status\n")
  cat("  analysis_exists: ", dir.exists(analysis_dir), "\n", sep = "")
  cat("  output_files: ", length(output_files), "\n", sep = "")
  cat("  figure_files: ", length(figure_files), "\n", sep = "")
  cat("  run_status: ", if (file.exists(run_status_file)) run_status_file else "<not found>", "\n", sep = "")
  cat("  figure_completeness: ", if (file.exists(figure_qc_file)) figure_qc_file else "<not found>", "\n", sep = "")
  cat("\nSTATUS_COMPLETE\nResults: ", analysis_dir, "\n", sep = "")
  quit(status = 0L)
}

if (dry_run) {
  cat("\nDRY_RUN_PASSED: no analysis was started.\n")
  quit(status = 0L)
}

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(analysis_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(analysis_dir, "audit"), recursive = TRUE, showWarnings = FALSE)

env_values <- c(
  PRADEEP_PROJECT_DIR = project_dir,
  PRADEEP_DIR0 = dir0,
  PRADEEP_PHEDIR = phe_dir,
  PRADEEP_SCRIPT_ROOT = script_root,
  PRADEEP_HELPER_DIR = helper_dir,
  PRADEEP_OUTPUT_ROOT = output_root,
  PRADEEP_ANALYSIS_DIR = analysis_dir,
  PRADEEP_PANEL_MODE = panel_mode,
  PRADEEP_PANEL_INPUT_DIR = panel_input_dir,
  PRADEEP_SOURCE_ALL_RDS = source_all_rds,
  PRADEEP_SOURCE_PROT_RDS = source_prot_rds,
  PRADEEP_ALL_RDS = all_rds,
  PRADEEP_PROT_RDS = prot_rds,
  PRADEEP_RAW_PROTEIN = raw_protein,
  PRADEEP_PHENO_TSV = pheno_tsv,
  PRADEEP_PROTEIN_MAP = protein_map,
  PRADEEP_GWAS_DIR = gwas_dir,
  PRADEEP_PQTL_DIR = pqtl_dir,
  PRADEEP_FINNGEN_DOWNLOAD_DIR = finngen_download_dir,
  UKBPPP_ANALYSIS_PROJECT = analysis_project,
  UKBPPP_OUTCOME_SUBSET = paste(outcome_subset, collapse = ","),
  UKBPPP_FOLLOW_END = follow_end,
  UKBPPP_REQUIRE_FOLLOW_END = follow_end,
  UKBPPP_PPP_BED_FILE = bed_file,
  UKBPPP_POPULATION = "ALL",
  UKBPPP_ORIGINAL_15K_PANEL = if (panel_mode == "1.5k") "1" else "0",
  UKBPPP_MAX_PROTEINS = "0",
  UKBPPP_MAX_MISSING_PROTEIN = "0.10",
  UKBPPP_MAX_MISSING_INDIVIDUAL = "0.10",
  UKBPPP_SKIP_KNN = "0",
  # A missing pairwise file falls back to UKB field 22011 related-pair groups.
  # --allow-kinship-fallback only acknowledges that fallback; it must not skip exclusion.
  UKBPPP_SKIP_RELATEDNESS = "0",
  UKBPPP_RELATEDNESS_CUTOFF = "0.0884",
  UKBPPP_DATE_SOURCE = "fod_ref",
  UKBPPP_DATE_SOURCE_CAD = "fod_ref",
  UKBPPP_DATE_SOURCE_AF = "fod_ref",
  UKBPPP_DATE_SOURCE_HF = "fod_ref",
  UKBPPP_DATE_SOURCE_AS = "fod_icd10",
  UKBPPP_LASSO_NFOLDS = "10",
  UKBPPP_GLMNET_MAXIT = "1000000",
  UKBPPP_FORCE_LASSO = "1",
  UKBPPP_REQUIRE_PAPER_FIGURES = "1",
  UKBPPP_WORKERS = as.character(workers),
  UKBPPP_PLINK2_WSL = plink2,
  UKBPPP_LD_BFILE = ld_bfile,
  UKBPPP_LD_PFILE_TEMPLATE = ld_pfile_template
)
if (nzchar(kinship_file)) env_values <- c(env_values, UKBPPP_KINSHIP_FILE = kinship_file)
do.call(Sys.setenv, as.list(env_values))

manifest <- data.frame(
  key = c(names(env_values), "PROFILE", "REQUESTED_STEP", "SELECTED_STEPS", "STARTED_AT"),
  value = c(unname(env_values), profile, requested_step, paste(selected_steps, collapse = ","), format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(analysis_dir, "audit", "command_manifest.csv"), row.names = FALSE, na = "")

if ("1" %in% selected_steps && file.exists(base_file) && !resume) {
  stop("Existing Step 1 output found: ", base_file, "\nUse --resume or choose a new --analysis-project.", call. = FALSE)
}

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
if ("1" %in% selected_steps && panel_mode == "1.5k" && !(resume && file.exists(base_file))) {
  cat("\n===== PREPARE 1.5K INPUTS =====\n")
  prepare_script <- file.path(f_dir, "00_prepare_panel_inputs.R")
  prepare_exit <- system2(rscript, c("--vanilla", shQuote(prepare_script)))
  if (!identical(prepare_exit, 0L)) {
    stop("1.5k input preparation failed with exit code ", prepare_exit, ".", call. = FALSE)
  }
  if (!all(file.exists(c(all_rds, prot_rds)))) {
    stop("1.5k input preparation finished without generated all.rds/prot.rds.", call. = FALSE)
  }
}
run_status <- list()
for (step_id in selected_steps) {
  def <- step_defs[[step_id]]
  script <- file.path(f_dir, def$script)
  expected <- file.path(output_dir, def$output)
  started <- Sys.time()
  if (resume && file.exists(expected)) {
    cat("\nSKIP step ", step_id, " (output exists): ", expected, "\n", sep = "")
    exit_code <- 0L
    state <- "SKIPPED_RESUME"
  } else {
    cat("\n===== STEP ", step_id, ": ", def$label, " =====\n", sep = "")
    exit_code <- system2(rscript, c("--vanilla", shQuote(script)))
    state <- if (identical(exit_code, 0L) && file.exists(expected)) "COMPLETED" else "FAILED"
  }
  ended <- Sys.time()
  run_status[[length(run_status) + 1L]] <- data.frame(
    step = step_id, label = def$label, status = state, exit_code = exit_code,
    expected_output = expected, started_at = format(started, "%Y-%m-%d %H:%M:%S %z"),
    ended_at = format(ended, "%Y-%m-%d %H:%M:%S %z"), stringsAsFactors = FALSE
  )
  write.csv(do.call(rbind, run_status), file.path(analysis_dir, "audit", "run_status.csv"), row.names = FALSE, na = "")
  if (!identical(exit_code, 0L)) stop("Step ", step_id, " failed with exit code ", exit_code, ".", call. = FALSE)
  if (!file.exists(expected)) stop("Step ", step_id, " finished without expected output: ", expected, call. = FALSE)
}

if (any(c("5", "6", "8") %in% selected_steps)) {
  source(file.path(f_dir, "09_finalize_manifest.R"), chdir = TRUE)
}

cat("\nPRADEEP_SELECTED_STEPS_COMPLETE\n")
cat("Results: ", analysis_dir, "\n", sep = "")
cat("Figures: ", file.path(output_dir, "figures"), "\n", sep = "")
