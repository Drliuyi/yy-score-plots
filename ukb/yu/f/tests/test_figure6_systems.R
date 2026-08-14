suppressWarnings(suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
}))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))))
source(file.path(root, "R", "full", "00_core.R"))
source(file.path(root, "R", "full", "08_systems_biology.R"))

specification <- list(collection = "C2", subcollection = "CP:REACTOME")
stopifnot(
  identical(names(yur_msigdb_arguments(specification, c("species", "collection", "subcollection"))),
            c("species", "collection", "subcollection")),
  identical(names(yur_msigdb_arguments(specification, c("species", "category", "subcategory"))),
            c("species", "category", "subcategory"))
)
version_test <- yur_msigdb_version(data.table(gs_name = "TEST"), package_version = "test")
stopifnot("db_version" %in% names(version_test), nzchar(version_test$db_version[[1]]))
old_specs <- yur_msigdb_specifications(c("species", "category", "subcategory"))
new_specs <- yur_msigdb_specifications(c("species", "collection", "subcollection"))
stopifnot(
  identical(vapply(old_specs, `[[`, character(1), "subcollection")[2], "CP:KEGG"),
  all(c("CP:KEGG_LEGACY", "CP:KEGG_MEDICUS") %in%
        vapply(new_specs, `[[`, character(1), "subcollection"))
)

tmp <- tempfile("figure6_systems_test_")
dir.create(file.path(tmp, "05_cox"), recursive = TRUE)
dir.create(file.path(tmp, "14_enrichment"), recursive = TRUE)
dir.create(file.path(tmp, "15_figures"), recursive = TRUE)

beta_values <- c(.2, .1, -.2, .01, .3, -.1, .01, .02)
cox <- data.table(
  scope = "full_incident",
  outcome_id = rep(c("cad", "heart_failure"), each = 4),
  outcome_label = rep(c("Coronary artery disease", "Heart failure"), each = 4),
  feature_id = rep(c("A1", "A2", "B1", "C1"), 2),
  protein = rep(c("IL6", "IL6", "TNF", "CRP"), 2),
  panel = "Cardiometabolic", olink_id = NA_character_, mapping_status = "PASS",
  p = c(1e-9, 2e-8, 3e-7, .2, 1e-10, 2e-6, .3, .4),
  beta = beta_values, hr = exp(beta_values), bonferroni_threshold = .05 / 4
)
cox[, bonferroni_significant := p < bonferroni_threshold]
fwrite(cox, file.path(tmp, "05_cox", "table_s2_incident_associations.csv.gz"))

cfg <- list(
  paths = list(
    cox = file.path(tmp, "05_cox"), enrichment = file.path(tmp, "14_enrichment"),
    figures = file.path(tmp, "15_figures")
  ),
  systems_top_n_per_outcome = 15L, force = TRUE
)
yur_prepare_figure6_systems(cfg)

manifest <- jsonlite::read_json(
  file.path(tmp, "14_enrichment", "local_cox_systems", "figure6_systems_input_manifest.json"),
  simplifyVector = TRUE
)
stopifnot(
  manifest$tested_assays == 4L,
  manifest$tested_gene_symbols == 3L,
  manifest$significant_assay_level_associations == 5L,
  manifest$significant_gene_level_associations == 3L,
  manifest$significant_unique_proteins == 2L,
  identical(manifest$figure6a_input_is_used, FALSE)
)

graph <- graph_from_edgelist(matrix(c(
  "A", "B", "A", "C", "B", "C", "C", "D"
), ncol = 2, byrow = TRUE), directed = FALSE)
mnc <- yur_mnc_score(graph)
stopifnot(mnc[gene_symbol == "C", mnc] >= mnc[gene_symbol == "D", mnc])

source_nodes <- yur_sankey_node_layout(
  data.table(label = c("Outcome A", "Outcome B"), weight = c(2, 4)),
  c("Outcome A", "Outcome B"), x = 0
)
target_nodes <- yur_sankey_node_layout(
  data.table(label = c("Protein 1", "Protein 2"), weight = c(3, 3)),
  c("Protein 1", "Protein 2"), x = 1
)
stopifnot(
  abs(source_nodes[label == "Outcome B", node_ymax - node_ymin] /
        source_nodes[label == "Outcome A", node_ymax - node_ymin] - 2) < 1e-8
)
weighted_links <- data.table(
  outcome = c("Outcome A", "Outcome A", "Outcome B"),
  protein = c("Protein 1", "Protein 2", "Protein 1"),
  weight = c(1, 1, 2), edge_color = c("Outcome A", "Outcome A", "Outcome B")
)
allocated <- yur_sankey_allocate_links(weighted_links, source_nodes, target_nodes, "outcome", "protein")
ribbons <- yur_sankey_ribbons(allocated, .03, .97, "test")
stopifnot(nrow(ribbons) > 0L, uniqueN(ribbons$ribbon_group) == nrow(weighted_links))

cat("FIGURE 6 SYSTEMS TESTS PASSED\n")
