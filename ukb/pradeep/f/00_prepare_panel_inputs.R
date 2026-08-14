locate_repro_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")))
  }
  project_dir <- Sys.getenv("PRADEEP_PROJECT_DIR", unset = "")
  if (nzchar(project_dir)) return(file.path(project_dir, "f"))
  stop("Cannot locate pradeep/f. Run this step through pradeep.sh.", call. = FALSE)
}

this_dir <- locate_repro_script_dir()
suppressPackageStartupMessages(library(data.table))

source_all <- Sys.getenv("PRADEEP_SOURCE_ALL_RDS", unset = "")
raw_protein <- Sys.getenv("PRADEEP_RAW_PROTEIN", unset = "")
protein_map <- Sys.getenv("PRADEEP_PROTEIN_MAP", unset = "")
panel_input_dir <- Sys.getenv("PRADEEP_PANEL_INPUT_DIR", unset = "")

required <- c(source_all = source_all, raw_protein = raw_protein,
              protein_map = protein_map, panel_input_dir = panel_input_dir)
missing_values <- names(required)[!nzchar(required)]
if (length(missing_values) > 0L) {
  stop("Missing panel-builder settings: ", paste(missing_values, collapse = ", "), call. = FALSE)
}
missing_files <- names(required)[names(required) != "panel_input_dir" & !file.exists(required)]
if (length(missing_files) > 0L) {
  stop("Missing 1.5k source files: ", paste(required[missing_files], collapse = ", "), call. = FALSE)
}

dir.create(panel_input_dir, recursive = TRUE, showWarnings = FALSE)
out_all <- file.path(panel_input_dir, "all.rds")
out_prot <- file.path(panel_input_dir, "prot.rds")
out_map <- file.path(panel_input_dir, "strict_1459_assay_mapping.csv")
out_manifest <- file.path(panel_input_dir, "panel_input_manifest.csv")
out_qc <- file.path(panel_input_dir, "panel_input_qc.csv")

file_signature <- function(path) {
  info <- file.info(path)
  paste(normalizePath(path, winslash = "/", mustWork = TRUE), info$size,
        format(info$mtime, "%Y-%m-%d %H:%M:%OS6 %z"), sep = "|")
}

source_signatures <- c(
  source_all = file_signature(source_all),
  raw_protein = file_signature(raw_protein),
  protein_map = file_signature(protein_map)
)

if (all(file.exists(c(out_all, out_prot, out_map, out_manifest, out_qc)))) {
  old_manifest <- fread(out_manifest)
  old_signatures <- setNames(old_manifest$value[old_manifest$key %in% names(source_signatures)],
                             old_manifest$key[old_manifest$key %in% names(source_signatures)])
  if (identical(unname(old_signatures[names(source_signatures)]), unname(source_signatures))) {
    message("1.5k generated inputs match the current sources; reusing: ", panel_input_dir)
    quit(status = 0L)
  }
}

norm_id <- function(x) gsub("[^A-Z0-9]", "", toupper(x))
removed_by_paper <- c("CTSS", "NPM1", "PCOLCE", "TACSTD2")

message("Preparing the Pradeep 1.5k assay panel from the un-imputed protein table.")
map <- fread(protein_map)
if (!"Assay" %in% names(map)) stop("The 1.5k protein map lacks an Assay column.", call. = FALSE)
wanted <- setdiff(unique(as.character(map$Assay)), removed_by_paper)
if (length(wanted) != 1459L) {
  stop("Expected 1459 paper proteins after exclusions; obtained ", length(wanted), call. = FALSE)
}

header <- names(fread(raw_protein, nrows = 0L))
eid_candidates <- c("eid", "EID", "f.eid")
eid_col <- eid_candidates[eid_candidates %in% header][1L]
if (is.na(eid_col)) stop("The raw protein table has no eid column.", call. = FALSE)
raw_proteins <- setdiff(header, eid_col)
idx <- match(norm_id(wanted), norm_id(raw_proteins))
missing_assays <- wanted[is.na(idx)]
if (length(missing_assays) > 0L) {
  fwrite(data.table(missing_assay = missing_assays),
         file.path(panel_input_dir, "missing_1.5k_assays.csv"))
  stop("Missing Pradeep assays: ", paste(missing_assays, collapse = ", "), call. = FALSE)
}
matched_raw <- raw_proteins[idx]
if (anyDuplicated(norm_id(matched_raw))) stop("Ambiguous normalized assay mapping detected.", call. = FALSE)

mapping <- data.table(
  paper_assay = wanted,
  raw_variable = matched_raw,
  normalized_id = norm_id(wanted)
)
fwrite(mapping, out_map)

prot <- fread(raw_protein, select = c(eid_col, matched_raw))
setnames(prot, c(eid_col, matched_raw), c("eid", wanted))
prot[, eid := as.character(eid)]
generated_protein_cols <- setdiff(names(prot), "eid")
protein_missing_n <- sum(vapply(
  generated_protein_cols,
  function(protein) sum(is.na(prot[[protein]])),
  numeric(1L)
))
if (protein_missing_n == 0L) {
  stop("The raw protein table contains no missing values and appears pre-imputed.", call. = FALSE)
}
saveRDS(prot, out_prot, compress = FALSE)

message("Preparing the 1.5k phenotype input from the source all.rds.")
all <- as.data.table(readRDS(source_all))
if (!"eid" %in% names(all)) stop("The source all.rds has no eid column.", call. = FALSE)
all[, eid := as.character(eid)]
race_candidates <- c("ethnic.c", "gen_ethnicity", "ethnicity")
race_col <- race_candidates[race_candidates %in% names(all)][1L]
if (is.na(race_col)) stop("Race variable not found in source all.rds.", call. = FALSE)
pcs <- paste0("PC", 1:10)
if (!base::all(pcs %in% names(all))) stop("PC1-PC10 are incomplete in source all.rds.", call. = FALSE)
source_all_n <- nrow(all)
all <- all[
  !is.na(get(race_col)) & trimws(as.character(get(race_col))) != "" &
    complete.cases(all[, ..pcs])
]
saveRDS(all, out_all, compress = FALSE)

qc <- data.table(
  item = c(
    "panel_mode", "source_all_n", "generated_all_n", "raw_protein_n",
    "generated_protein_n", "generated_protein_count", "raw_missing_values",
    "generated_missing_values", "matched_assays"
  ),
  value = c(
    "1.5k", source_all_n, nrow(all), nrow(prot), nrow(prot), ncol(prot) - 1L,
    protein_missing_n, protein_missing_n, nrow(mapping)
  )
)
fwrite(qc, out_qc)

manifest <- as.data.table(data.frame(
  key = c(names(source_signatures), "panel_mode", "effective_all_rds", "effective_prot_rds",
          "mapping_file", "created_at"),
  value = c(unname(source_signatures), "1.5k", normalizePath(out_all, winslash = "/"),
            normalizePath(out_prot, winslash = "/"), normalizePath(out_map, winslash = "/"),
            format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  stringsAsFactors = FALSE
))
fwrite(manifest, out_manifest)

message("PRADEEP_15K_INPUT_PREPARATION_COMPLETE")
message("  participants in generated protein input: ", nrow(prot))
message("  participants in generated phenotype input: ", nrow(all))
message("  proteins: ", ncol(prot) - 1L)
message("  remaining raw protein missing values: ", protein_missing_n)
message("  output: ", panel_input_dir)
