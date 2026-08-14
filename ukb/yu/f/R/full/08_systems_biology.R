suppressWarnings(suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(digest)
}))

yur_systems_dir <- function(cfg) {
  path <- file.path(cfg$paths$enrichment, "local_cox_systems")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

yur_systems_required_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Figure 6 systems stage requires R packages: ", paste(missing, collapse = ", "),
      ". Install them before rerunning.", call. = FALSE
    )
  }
}

yur_systems_font_family <- function() {
  if (.Platform$OS.type == "windows") "Arial" else "Helvetica"
}

yur_systems_cox_file <- function(cfg) {
  candidates <- c(
    file.path(cfg$paths$cox, "table_s2_incident_associations.csv.gz"),
    file.path(cfg$paths$cox, "table_s2_incident_associations.csv")
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop("Merged local Cox table is missing; run cox_merge first.", call. = FALSE)
  hit[[1]]
}

yur_clean_gene_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "<NA>")] <- NA_character_
  x
}

yur_prepare_figure6_systems <- function(cfg) {
  out_dir <- yur_systems_dir(cfg)
  cox_file <- yur_systems_cox_file(cfg)
  cox <- fread(cox_file, showProgress = FALSE)
  required <- c(
    "scope", "outcome_id", "outcome_label", "feature_id", "protein", "p", "beta",
    "hr", "bonferroni_threshold", "bonferroni_significant"
  )
  missing <- setdiff(required, names(cox))
  if (length(missing)) stop("Cox table is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  cox <- cox[scope == "full_incident"]
  cox[, gene_symbol := yur_clean_gene_symbol(protein)]
  cox <- cox[!is.na(gene_symbol) & is.finite(p)]
  if (!nrow(cox)) stop("No valid full-incident Cox rows were found.", call. = FALSE)

  background <- unique(cox[, .(gene_symbol, feature_id, protein, panel, olink_id, mapping_status)])
  setorder(background, gene_symbol, feature_id)
  background[, assay_n_for_symbol := .N, by = gene_symbol]

  significant_assay <- cox[bonferroni_significant %in% TRUE]
  if (!nrow(significant_assay)) stop("No Bonferroni-significant local Cox associations were found.", call. = FALSE)
  setorder(significant_assay, outcome_id, gene_symbol, p)
  significant <- significant_assay[, .SD[1L], by = .(outcome_id, gene_symbol)]
  foreground <- significant[, .(
    association_n = .N,
    outcome_n = uniqueN(outcome_id),
    min_p = min(p),
    max_abs_beta = max(abs(beta)),
    direction_summary = if (all(beta > 0)) "positive" else if (all(beta < 0)) "negative" else "mixed"
  ), by = gene_symbol]
  setorder(foreground, min_p)

  top_n <- as.integer(cfg$systems_top_n_per_outcome %||% 15L)
  top <- significant[order(outcome_id, p), head(.SD, top_n), by = .(outcome_id, outcome_label)]
  top[, rank_within_outcome := seq_len(.N), by = outcome_id]

  outcome_counts <- significant[, .(
    significant_associations = .N,
    unique_proteins = uniqueN(gene_symbol),
    positive = sum(beta > 0),
    negative = sum(beta < 0),
    top_n_exported = min(.N, top_n)
  ), by = .(outcome_id, outcome_label)]
  setorder(outcome_counts, -significant_associations)

  yur_write_csv(background, file.path(out_dir, "figure6_systems_background_all_tested_assays.csv"))
  yur_write_csv(unique(background[, .(gene_symbol)]), file.path(out_dir, "figure6_systems_background_unique_genes.csv"))
  yur_write_csv(significant_assay, file.path(out_dir, "figure6_systems_significant_associations_assay_level.csv.gz"))
  yur_write_csv(significant, file.path(out_dir, "figure6_systems_significant_associations.csv.gz"))
  yur_write_csv(foreground, file.path(out_dir, "figure6_systems_foreground_unique_genes.csv"))
  yur_write_csv(top, file.path(out_dir, "figure6c_top15_by_outcome.csv"))
  yur_write_csv(outcome_counts, file.path(out_dir, "figure6_systems_outcome_counts.csv"))
  yur_write_json(list(
    status = "PASS",
    source = normalizePath(cox_file, winslash = "/", mustWork = TRUE),
    source_sha256 = yur_sha_file(cox_file),
    tested_assays = uniqueN(cox$feature_id),
    tested_gene_symbols = uniqueN(cox$gene_symbol),
    significant_assay_level_associations = nrow(significant_assay),
    significant_gene_level_associations = nrow(significant),
    significant_unique_proteins = nrow(foreground),
    outcomes_with_significant_proteins = uniqueN(significant$outcome_id),
    top_n_per_outcome = top_n,
    inference = "Local full-incident Cox results; Bonferroni P < 0.05 / retained assay count.",
    figure6a_input_is_used = FALSE
  ), file.path(out_dir, "figure6_systems_input_manifest.json"))
}

yur_http_post_tsv <- function(url, fields, raw_path, timeout_seconds = 300L, retries = 3L) {
  yur_systems_required_packages(c("httr2"))
  dir.create(dirname(raw_path), recursive = TRUE, showWarnings = FALSE)
  last_error <- NULL
  for (attempt in seq_len(retries)) {
    response <- tryCatch({
      request <- httr2::request(url)
      request <- do.call(httr2::req_body_form, c(list(request), fields))
      request <- httr2::req_timeout(request, seconds = timeout_seconds)
      request <- httr2::req_retry(request, max_tries = 2L)
      httr2::req_perform(request)
    }, error = function(e) e)
    if (!inherits(response, "error")) {
      body <- httr2::resp_body_string(response)
      if (httr2::resp_status(response) == 200L && nzchar(body)) {
        writeLines(body, raw_path, useBytes = TRUE)
        return(fread(text = body, sep = "\t", header = TRUE, showProgress = FALSE))
      }
      last_error <- paste0("HTTP ", httr2::resp_status(response))
    } else {
      last_error <- conditionMessage(response)
    }
    Sys.sleep(attempt * 2)
  }
  stop("External API request failed after retries: ", url, "; ", last_error, call. = FALSE)
}

yur_string_map <- function(genes, cfg, stem) {
  genes <- sort(unique(yur_clean_gene_symbol(genes)))
  genes <- genes[!is.na(genes)]
  out_dir <- yur_systems_dir(cfg)
  url <- paste0(cfg$string_api_base, "/api/tsv/get_string_ids")
  mapped <- yur_http_post_tsv(
    url,
    list(
      identifiers = paste(genes, collapse = "\r"), species = "9606", limit = "1",
      echo_query = "1", caller_identity = "yu_protein_analysis"
    ),
    file.path(out_dir, paste0(stem, "_string_mapping_raw.tsv")),
    timeout_seconds = 300L
  )
  required <- c("queryItem", "stringId", "preferredName")
  missing <- setdiff(required, names(mapped))
  if (length(missing)) stop("STRING mapping response is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  mapped[, query_gene := yur_clean_gene_symbol(queryItem)]
  mapped[, preferred_gene := yur_clean_gene_symbol(preferredName)]
  setorder(mapped, query_gene, stringId)
  mapped <- mapped[, .SD[1L], by = query_gene]
  submitted <- data.table(query_gene = genes)
  mapped <- merge(submitted, mapped, by = "query_gene", all.x = TRUE, sort = FALSE)
  mapped[, mapping_status := fifelse(!is.na(stringId) & nzchar(stringId), "MAPPED", "UNMAPPED")]
  yur_write_csv(mapped, file.path(out_dir, paste0(stem, "_string_mapping.csv")))
  mapped
}

yur_clean_pathway_label <- function(x) {
  x <- gsub("^(WP|KEGG|REACTOME|BIOCARTA|PID|GOBP)_", "", x)
  x <- gsub("_", " ", x, fixed = TRUE)
  tools::toTitleCase(tolower(x))
}

yur_msigdb_arguments <- function(specification, formal_names = names(formals(msigdbr::msigdbr))) {
  args <- list(species = "Homo sapiens")
  if ("collection" %in% formal_names) {
    args$collection <- specification$collection
    args$subcollection <- specification$subcollection
  } else {
    args$category <- specification$collection
    args$subcategory <- specification$subcollection
  }
  args
}

yur_msigdb_version <- function(value, package_version = NULL) {
  if (!"db_version" %in% names(value)) {
    package_version <- package_version %||% as.character(utils::packageVersion("msigdbr"))
    value[, db_version := paste0("msigdbr-", package_version)]
  }
  value
}

yur_msigdb_specifications <- function(formal_names = names(formals(msigdbr::msigdbr))) {
  kegg <- if ("collection" %in% formal_names) c("CP:KEGG_LEGACY", "CP:KEGG_MEDICUS") else "CP:KEGG"
  c(
    list(list(collection = "C2", subcollection = "CP:WIKIPATHWAYS", source = "Wiki")),
    lapply(kegg, function(x) list(collection = "C2", subcollection = x, source = "KEGG")),
    list(
      list(collection = "C2", subcollection = "CP:BIOCARTA", source = "Canonical"),
      list(collection = "C2", subcollection = "CP:PID", source = "Canonical"),
      list(collection = "C2", subcollection = "CP:REACTOME", source = "Reactome")
    )
  )
}

yur_msigdb_pathways <- function(background) {
  yur_systems_required_packages(c("msigdbr", "AnnotationDbi", "org.Hs.eg.db", "GO.db"))
  specifications <- yur_msigdb_specifications()
  pathway_parts <- lapply(specifications, function(specification) {
    value <- yur_msigdb_version(as.data.table(do.call(msigdbr::msigdbr, yur_msigdb_arguments(specification))))
    value <- value[yur_clean_gene_symbol(gene_symbol) %chin% background]
    unique(value[, .(
      source = specification$source, term = gs_name,
      description = yur_clean_pathway_label(gs_name),
      gene_symbol = yur_clean_gene_symbol(gene_symbol),
      msigdb_version = db_version
    )])
  })
  go <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db, keys = background, keytype = "SYMBOL",
    columns = c("GOALL", "ONTOLOGYALL")
  ))
  go <- as.data.table(go)[ONTOLOGYALL == "BP" & !is.na(GOALL)]
  go_terms <- suppressMessages(AnnotationDbi::select(
    GO.db::GO.db, keys = unique(go$GOALL), keytype = "GOID", columns = "TERM"
  ))
  go <- merge(go, as.data.table(go_terms), by.x = "GOALL", by.y = "GOID", all.x = TRUE)
  pathway_parts[[length(pathway_parts) + 1L]] <- unique(go[, .(
    source = "GO:BP", term = GOALL, description = TERM,
    gene_symbol = yur_clean_gene_symbol(SYMBOL),
    msigdb_version = paste0("org.Hs.eg.db-", as.character(utils::packageVersion("org.Hs.eg.db")))
  )])
  pathways <- rbindlist(pathway_parts, fill = TRUE)
  unique(pathways[!is.na(gene_symbol)])
}

yur_run_figure6_enrichment <- function(cfg) {
  out_dir <- yur_systems_dir(cfg)
  fg_file <- file.path(out_dir, "figure6_systems_foreground_unique_genes.csv")
  bg_file <- file.path(out_dir, "figure6_systems_background_unique_genes.csv")
  if (!all(file.exists(c(fg_file, bg_file)))) stop("Run systems_prepare first.", call. = FALSE)
  foreground <- unique(yur_clean_gene_symbol(fread(fg_file)$gene_symbol))
  background <- unique(yur_clean_gene_symbol(fread(bg_file)$gene_symbol))
  foreground <- intersect(foreground, background)
  pathways <- yur_msigdb_pathways(background)
  set_sizes <- pathways[, .(gene_set_size = uniqueN(gene_symbol)), by = .(source, term, description, msigdb_version)]
  overlaps <- pathways[gene_symbol %chin% foreground, .(
    number_of_genes = uniqueN(gene_symbol),
    input_genes = paste(sort(unique(gene_symbol)), collapse = ";")
  ), by = .(source, term, description, msigdb_version)]
  enrichment <- merge(set_sizes, overlaps, by = c("source", "term", "description", "msigdb_version"), all.x = TRUE)
  enrichment[is.na(number_of_genes), `:=`(number_of_genes = 0L, input_genes = "")]
  enrichment <- enrichment[gene_set_size >= 3L & number_of_genes >= 2L]
  enrichment[, p_value := phyper(
    number_of_genes - 1L, gene_set_size, length(background) - gene_set_size,
    length(foreground), lower.tail = FALSE
  )]
  enrichment[, fdr := p.adjust(p_value, method = "BH")]
  enrichment[, significant := is.finite(fdr) & fdr < as.numeric(cfg$systems_enrichment_fdr %||% 0.05)]
  enrichment[, minus_log10_fdr := -log10(pmax(fdr, .Machine$double.xmin))]
  setorder(enrichment, fdr, -number_of_genes)
  yur_write_csv(pathways, file.path(out_dir, "msigdb_pathway_gene_sets.csv.gz"))
  yur_write_csv(enrichment, file.path(out_dir, "enrichment_results_all_msigdb_categories.csv.gz"))
  main <- enrichment[significant == TRUE]
  missing_sources <- setdiff(c("Wiki", "KEGG", "Canonical", "GO:BP", "Reactome"), unique(main$source))
  if (length(missing_sources)) {
    stop("No FDR-significant local enrichment result for: ", paste(missing_sources, collapse = ", "), call. = FALSE)
  }
  yur_write_csv(main, file.path(cfg$paths$enrichment, "enrichment_results.csv"))
  yur_write_json(list(
    status = "PASS", method = "MSigDB over-representation analysis with local measured-protein background",
    msigdb_version = paste(sort(unique(pathways$msigdb_version)), collapse = ";"),
    foreground_proteins = length(foreground), background_proteins = length(background),
    fdr_threshold = as.numeric(cfg$systems_enrichment_fdr %||% 0.05),
    significant_terms = nrow(main),
    source_categories = sort(unique(main$source)),
    gene_set_sha256 = yur_sha_file(file.path(out_dir, "msigdb_pathway_gene_sets.csv.gz"))
  ), file.path(out_dir, "enrichment_manifest.json"))
}

yur_download_trrust <- function(cfg) {
  yur_systems_required_packages("curl")
  out_dir <- yur_systems_dir(cfg)
  path <- file.path(out_dir, "trrust_rawdata_human.tsv")
  if (!file.exists(path) || cfg$force) {
    curl::curl_download(cfg$trrust_human_url, path, quiet = FALSE, mode = "wb")
  }
  trrust <- fread(path, header = FALSE, sep = "\t", quote = "", showProgress = FALSE)
  if (ncol(trrust) < 4L || nrow(trrust) < 1000L) stop("Downloaded TRRUST file failed structural QC.", call. = FALSE)
  setnames(trrust, names(trrust)[1:4], c("tf", "target", "regulation", "pmid"))
  trrust[, `:=`(tf = yur_clean_gene_symbol(tf), target = yur_clean_gene_symbol(target))]
  trrust[!is.na(tf) & !is.na(target)]
}

yur_run_figure6_tf <- function(cfg) {
  out_dir <- yur_systems_dir(cfg)
  top_file <- file.path(out_dir, "figure6c_top15_by_outcome.csv")
  if (!file.exists(top_file)) stop("Run systems_prepare first.", call. = FALSE)
  top <- fread(top_file)
  trrust <- yur_download_trrust(cfg)
  edges <- merge(top, trrust, by.x = "gene_symbol", by.y = "target", allow.cartesian = TRUE)
  if (!nrow(edges)) stop("No local top proteins mapped to TRRUST.", call. = FALSE)
  setnames(edges, "gene_symbol", "protein")
  setorder(edges, outcome_id, rank_within_outcome, tf)
  edges <- unique(edges, by = c("outcome_id", "protein", "tf", "regulation", "pmid"))
  compact <- edges[, .(
    supporting_records = .N,
    regulation = paste(sort(unique(na.omit(regulation))), collapse = ";"),
    pmid = paste(sort(unique(na.omit(as.character(pmid)))), collapse = ";")
  ), by = .(outcome_id, outcome_label, protein, tf, rank_within_outcome, p, beta, hr)]
  tf_summary <- compact[, .(
    outcome_n = uniqueN(outcome_id), protein_n = uniqueN(protein), edge_n = .N,
    supporting_records = sum(supporting_records)
  ), by = tf][order(-outcome_n, -protein_n, -supporting_records, tf)]
  max_tf <- as.integer(cfg$systems_max_tf %||% 46L)
  display_tfs <- head(tf_summary$tf, max_tf)
  compact_display <- compact[tf %chin% display_tfs]
  yur_write_csv(compact, file.path(out_dir, "figure6c_cvd_protein_tf_edges_all.csv"))
  yur_write_csv(compact_display, file.path(out_dir, "figure6c_cvd_protein_tf_edges.csv"))
  yur_write_csv(tf_summary, file.path(out_dir, "figure6c_tf_summary.csv"))
  yur_write_json(list(
    status = "PASS", database = "TRRUST v2 human downloadable relationship table",
    source_url = cfg$trrust_human_url,
    source_sha256 = yur_sha_file(file.path(out_dir, "trrust_rawdata_human.tsv")),
    outcomes = uniqueN(compact$outcome_id), proteins = uniqueN(compact$protein),
    transcription_factors_total = uniqueN(compact$tf), edges_total = nrow(compact),
    display_rule = "Top TFs by outcome count, protein count, supporting records, then symbol",
    transcription_factors_displayed = uniqueN(compact_display$tf),
    edges_displayed = nrow(compact_display)
  ), file.path(out_dir, "figure6c_tf_manifest.json"))
}

yur_mnc_score <- function(graph) {
  vertices <- igraph::V(graph)$name
  score <- vapply(vertices, function(vertex) {
    neighbors <- igraph::neighbors(graph, vertex, mode = "all")$name
    if (!length(neighbors)) return(0)
    subgraph <- igraph::induced_subgraph(graph, vids = neighbors)
    max(igraph::components(subgraph)$csize)
  }, numeric(1))
  data.table(gene_symbol = vertices, mnc = score)
}

yur_mcode_node_weights <- function(graph) {
  vertices <- igraph::V(graph)$name
  weights <- vapply(vertices, function(vertex) {
    neighborhood <- unique(c(vertex, igraph::neighbors(graph, vertex, mode = "all")$name))
    subgraph <- igraph::induced_subgraph(graph, vids = neighborhood)
    if (igraph::vcount(subgraph) < 3L) return(0)
    core <- igraph::coreness(subgraph)
    highest <- max(core)
    if (!is.finite(highest) || highest < 2) return(0)
    core_graph <- igraph::induced_subgraph(subgraph, vids = names(core)[core >= highest])
    highest * igraph::edge_density(core_graph, loops = FALSE)
  }, numeric(1))
  setNames(weights, vertices)
}

yur_haircut_cluster <- function(graph, vertices) {
  vertices <- unique(vertices)
  repeat {
    subgraph <- igraph::induced_subgraph(graph, vids = vertices)
    keep <- names(igraph::degree(subgraph))[igraph::degree(subgraph) > 1L]
    if (length(keep) == length(vertices) || length(keep) < 3L) break
    vertices <- keep
  }
  vertices
}

yur_mcode_clusters <- function(graph, node_score_cutoff = 0.2, max_depth = 100L) {
  weights <- yur_mcode_node_weights(graph)
  seeds <- names(sort(weights, decreasing = TRUE))
  claimed <- character()
  clusters <- list()
  for (seed in seeds) {
    if (seed %chin% claimed || weights[[seed]] <= 0) next
    cutoff <- weights[[seed]] * (1 - node_score_cutoff)
    members <- seed
    frontier <- seed
    depth <- 0L
    while (length(frontier) && depth < max_depth) {
      candidates <- unique(unlist(lapply(frontier, function(vertex) {
        igraph::neighbors(graph, vertex, mode = "all")$name
      }), use.names = FALSE))
      candidates <- setdiff(candidates, members)
      candidates <- candidates[weights[candidates] >= cutoff]
      if (!length(candidates)) break
      members <- unique(c(members, candidates))
      frontier <- candidates
      depth <- depth + 1L
    }
    members <- yur_haircut_cluster(graph, members)
    if (length(members) < 3L) next
    subgraph <- igraph::induced_subgraph(graph, vids = members)
    if (max(igraph::coreness(subgraph)) < 2L) next
    clusters[[length(clusters) + 1L]] <- list(
      seed = seed, members = members,
      score = igraph::vcount(subgraph) * igraph::edge_density(subgraph, loops = FALSE),
      nodes = igraph::vcount(subgraph), edges = igraph::ecount(subgraph)
    )
    claimed <- unique(c(claimed, members))
  }
  if (!length(clusters)) stop("No MCODE-compatible PPI cluster was identified.", call. = FALSE)
  clusters[order(vapply(clusters, `[[`, numeric(1), "score"), decreasing = TRUE)]
}

yur_run_figure6_ppi <- function(cfg) {
  yur_systems_required_packages("igraph")
  out_dir <- yur_systems_dir(cfg)
  fg_file <- file.path(out_dir, "figure6_systems_foreground_unique_genes.csv")
  if (!file.exists(fg_file)) stop("Run systems_prepare first.", call. = FALSE)
  foreground <- fread(fg_file)$gene_symbol
  map_file <- file.path(out_dir, "background_string_mapping.csv")
  mapping <- if (file.exists(map_file)) fread(map_file) else yur_string_map(foreground, cfg, "ppi_foreground")
  mapping <- mapping[query_gene %chin% foreground & mapping_status == "MAPPED"]
  if (nrow(mapping) < 10L) stop("Too few Cox-associated proteins mapped to STRING.", call. = FALSE)

  edges <- yur_http_post_tsv(
    paste0(cfg$string_api_base, "/api/tsv/network"),
    list(
      identifiers = paste(unique(mapping$stringId), collapse = "\r"), species = "9606",
      required_score = as.character(cfg$string_required_score %||% 150L),
      network_type = "functional", add_nodes = "0", caller_identity = "yu_protein_analysis"
    ),
    file.path(out_dir, "string_network_raw.tsv"), timeout_seconds = 600L
  )
  required <- c("preferredName_A", "preferredName_B", "score")
  missing <- setdiff(required, names(edges))
  if (length(missing)) stop("STRING network response is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  edges[, `:=`(
    from = yur_clean_gene_symbol(preferredName_A),
    to = yur_clean_gene_symbol(preferredName_B),
    score = as.numeric(score)
  )]
  edges <- edges[!is.na(from) & !is.na(to) & from != to & is.finite(score)]
  edges[, pair := fifelse(from < to, paste(from, to, sep = "||"), paste(to, from, sep = "||"))]
  setorder(edges, pair, -score)
  edges <- edges[, .SD[1L], by = pair]
  graph <- igraph::graph_from_data_frame(edges[, .(from, to, score)], directed = FALSE)
  graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE,
                            edge.attr.comb = list(score = "max", "ignore"))
  clusters <- yur_mcode_clusters(graph)
  membership <- rbindlist(lapply(seq_along(clusters), function(index) {
    cluster <- clusters[[index]]
    data.table(
      cluster = index, seed = cluster$seed, cluster_score = cluster$score,
      cluster_nodes = cluster$nodes, cluster_edges = cluster$edges,
      gene_symbol = cluster$members
    )
  }))
  yur_write_csv(membership, file.path(out_dir, "figure6d_mcode_compatible_clusters.csv"))

  main_members <- clusters[[1]]$members
  main_graph <- igraph::induced_subgraph(graph, vids = main_members)
  mnc <- yur_mnc_score(main_graph)
  degree <- data.table(gene_symbol = igraph::V(main_graph)$name, degree = as.numeric(igraph::degree(main_graph)))
  nodes <- merge(mnc, degree, by = "gene_symbol", all = TRUE)
  nodes[, mnc_rank := frank(-mnc, ties.method = "min")]
  nodes[, cluster := 1L]
  main_edges <- igraph::as_data_frame(main_graph, what = "edges") |> as.data.table()
  setnames(main_edges, intersect(c("from", "to"), names(main_edges)), c("from", "to"))
  yur_write_csv(nodes[order(-mnc, -degree)], file.path(out_dir, "figure6d_main_cluster_nodes.csv"))
  yur_write_csv(main_edges, file.path(out_dir, "figure6d_main_cluster_edges.csv.gz"))
  yur_write_csv(edges, file.path(out_dir, "figure6d_all_string_edges.csv.gz"))
  yur_write_csv(nodes, file.path(cfg$paths$enrichment, "network_results.csv"))
  yur_write_json(list(
    status = "PASS", database = "STRING v12.0", api_base = cfg$string_api_base,
    species = 9606, network_type = "functional",
    required_score = as.integer(cfg$string_required_score %||% 150L),
    mapped_input_proteins = nrow(mapping), all_network_nodes = igraph::vcount(graph),
    all_network_edges = igraph::ecount(graph), clusters = length(clusters),
    main_cluster_nodes = igraph::vcount(main_graph), main_cluster_edges = igraph::ecount(main_graph),
    cluster_method = "MCODE-compatible defaults: degree cutoff 2, node score cutoff 0.2, k-core 2, max depth 100, haircut on, fluff off",
    hub_method = "Maximum Neighborhood Component implemented from the cytoHubba definition",
    raw_response_sha256 = yur_sha_file(file.path(out_dir, "string_network_raw.tsv"))
  ), file.path(out_dir, "figure6d_ppi_manifest.json"))
}

yur_figure6_outcome_colors <- function(ids) {
  palette <- c(
    "#4E79A7", "#76B7B2", "#59A14F", "#EDC948", "#F28E2B", "#E15759", "#B07AA1",
    "#9C755F", "#BAB0AC", "#2F6B5F", "#7B8DAA", "#C77C8A", "#6F4E7C", "#86A873"
  )
  setNames(rep(palette, length.out = length(ids)), ids)
}

yur_figure6_outcome_short <- function(ids) {
  labels <- c(
    abdominal_aneurysm = "AAA", aortic_valve_stenosis = "AS",
    atrial_fibrillation = "AF", cad = "CAD", cardiomyopathy = "CM",
    deep_vein_thrombosis = "DVT", heart_failure = "HF",
    intracerebral_hemorrhage = "ICH", ischemic_stroke = "IS",
    peripheral_arterial_disease = "PAD", pulmonary_embolism = "PE",
    subarachnoid_hemorrhage = "SAH", thoracic_aneurysm = "TAA",
    transient_ischemic_attack = "TIA"
  )
  unname(labels[ids])
}

yur_bezier_paths <- function(edges, stage) {
  if (!nrow(edges)) return(data.table())
  steps <- 35L
  expanded <- edges[rep(seq_len(.N), each = steps)]
  expanded[, curve_t := rep(seq(0, 1, length.out = steps), times = nrow(edges))]
  expanded[, curve_group := paste(stage, edge_id, sep = "_")]
  expanded[, x := {
    t <- curve_t
    (1 - t)^3 * x0 + 3 * (1 - t)^2 * t * (x0 + 0.38 * (x1 - x0)) +
      3 * (1 - t) * t^2 * (x0 + 0.62 * (x1 - x0)) + t^3 * x1
  }]
  expanded[, y := {
    t <- curve_t
    (1 - t)^3 * y0 + 3 * (1 - t)^2 * t * y0 +
      3 * (1 - t) * t^2 * y1 + t^3 * y1
  }]
  expanded
}

yur_sankey_node_layout <- function(counts, order_labels, x, ymin = 0.02, ymax = 0.98, gap = 0.006) {
  counts <- copy(counts)
  setnames(counts, names(counts), c("label", "weight"))
  nodes <- merge(data.table(label = order_labels), counts, by = "label", all.x = TRUE, sort = FALSE)
  nodes[, order_id := match(label, order_labels)]
  setorder(nodes, order_id)
  nodes[is.na(weight), weight := 0]
  available <- ymax - ymin - gap * max(0, nrow(nodes) - 1L)
  unit_height <- available / sum(nodes$weight)
  node_heights <- nodes$weight * unit_height
  top_offsets <- if (length(node_heights) > 1L) c(0, cumsum(head(node_heights + gap, -1L))) else 0
  nodes[, `:=`(
    x = x,
    node_ymax = ymax - top_offsets,
    node_ymin = ymax - top_offsets - node_heights,
    unit_height = unit_height
  )]
  nodes[, y := (node_ymin + node_ymax) / 2]
  nodes
}

yur_sankey_allocate_links <- function(links, source_nodes, target_nodes, source_col, target_col) {
  links <- copy(links)
  links[, source_label := get(source_col)]
  links[, target_label := get(target_col)]
  source_map <- source_nodes[, .(
    source_label = label, source_ymin = node_ymin,
    source_unit = unit_height, source_center = y
  )]
  target_map <- target_nodes[, .(
    target_label = label, target_ymin = node_ymin,
    target_unit = unit_height, target_center = y
  )]
  links <- merge(links, source_map, by = "source_label", all.x = TRUE, sort = FALSE)
  links <- merge(links, target_map, by = "target_label", all.x = TRUE, sort = FALSE)

  # Offsets are accumulated from node_ymin (bottom) upward. Sorting centers in
  # ascending order keeps bottom links at the bottom and top links at the top,
  # minimizing avoidable crossings between adjacent Sankey stages.
  setorder(links, source_label, target_center, target_label)
  links[, source_offset := cumsum(weight) - weight, by = source_label]
  links[, `:=`(
    source_lower = source_ymin + source_offset * source_unit,
    source_upper = source_ymin + (source_offset + weight) * source_unit
  )]

  setorder(links, target_label, source_center, source_label)
  links[, target_offset := cumsum(weight) - weight, by = target_label]
  links[, `:=`(
    target_lower = target_ymin + target_offset * target_unit,
    target_upper = target_ymin + (target_offset + weight) * target_unit
  )]
  links
}

yur_sankey_ribbons <- function(links, x0, x1, stage) {
  if (!nrow(links)) return(data.table())
  steps <- 28L
  rbindlist(lapply(seq_len(nrow(links)), function(i) {
    t <- seq(0, 1, length.out = steps)
    ease <- t * t * (3 - 2 * t)
    x <- x0 + (x1 - x0) * t
    lower <- links$source_lower[[i]] + (links$target_lower[[i]] - links$source_lower[[i]]) * ease
    upper <- links$source_upper[[i]] + (links$target_upper[[i]] - links$source_upper[[i]]) * ease
    data.table(
      x = c(x, rev(x)), y = c(lower, rev(upper)),
      ribbon_group = paste(stage, i, sep = "_"),
      edge_color = links$edge_color[[i]]
    )
  }))
}

yur_trimmed_png_grob <- function(path) {
  image <- png::readPNG(path)
  if (length(dim(image)) == 3L && dim(image)[[3]] == 4L) {
    # Figure 6A is supplied as a transparent PNG with its own panel tag.
    # Remove that tag before trimming so the heatmap, rather than the old tag,
    # determines the visible top and left bounds in the composite figure.
    clear_rows <- seq_len(max(1L, round(dim(image)[[1]] * .02)))
    clear_cols <- seq_len(max(1L, round(dim(image)[[2]] * .25)))
    image[clear_rows, clear_cols, 4] <- 0
    alpha <- image[, , 4]
    rows <- which(rowSums(alpha > 0.005) > 0)
    cols <- which(colSums(alpha > 0.005) > 0)
    if (length(rows) && length(cols)) {
      image <- image[min(rows):max(rows), min(cols):max(cols), , drop = FALSE]
    }
  }
  list(
    grob = grid::rasterGrob(image, interpolate = TRUE),
    aspect = dim(image)[[1]] / dim(image)[[2]]
  )
}

yur_circle_positions <- function(labels, radius, start = pi / 2) {
  labels <- unique(as.character(labels))
  angles <- start - seq(0, 2 * pi, length.out = length(labels) + 1L)[seq_along(labels)]
  data.table(
    label = labels, angle = angles, x = radius * cos(angles), y = radius * sin(angles)
  )
}

yur_plot_figure6_systems <- function(cfg) {
  yur_systems_required_packages(c("ggplot2", "ggalluvial", "igraph", "ggrepel", "patchwork", "scales", "png"))
  font_family <- yur_systems_font_family()
  out_dir <- yur_systems_dir(cfg)
  enrichment_file <- file.path(cfg$paths$enrichment, "enrichment_results.csv")
  tf_file <- file.path(out_dir, "figure6c_cvd_protein_tf_edges.csv")
  node_file <- file.path(out_dir, "figure6d_main_cluster_nodes.csv")
  edge_file <- file.path(out_dir, "figure6d_main_cluster_edges.csv.gz")
  figure6a_file <- file.path(cfg$paths$figures, "figure6a_local_prs_protein_heatmap.png")
  required_files <- c(enrichment_file, tf_file, node_file, edge_file, figure6a_file)
  if (!all(file.exists(required_files))) {
    stop("Systems result files are incomplete; run systems_enrichment, systems_tf and systems_ppi first.", call. = FALSE)
  }

  enrichment <- fread(enrichment_file)
  category_order <- c("Wiki", "KEGG", "Canonical", "GO:BP", "Reactome")
  top_terms <- enrichment[order(fdr), head(.SD, 2L), by = source]
  top_terms <- top_terms[order(minus_log10_fdr)]
  top_terms[, description_wrapped := vapply(
    description, function(value) paste(strwrap(value, width = 34), collapse = "\n"), character(1)
  )]
  top_terms[, term_label := factor(description_wrapped, levels = unique(description_wrapped))]
  top_terms[, source := factor(source, levels = category_order)]
  source_colors <- c(
    "Wiki" = "#31466E", "KEGG" = "#4F7F99", "Canonical" = "#80A5A5",
    "GO:BP" = "#E97E4D", "Reactome" = "#E9B64F"
  )
  p6b <- ggplot2::ggplot(top_terms, ggplot2::aes(minus_log10_fdr, term_label)) +
    ggplot2::geom_point(ggplot2::aes(size = number_of_genes, color = source), alpha = .95) +
    ggplot2::scale_color_manual(values = source_colors, drop = FALSE) +
    ggplot2::scale_size_continuous(range = c(3.2, 8.6)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(.02, .05))) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 1, byrow = TRUE, order = 1),
      size = "none"
    ) +
    ggplot2::labs(x = expression(-log[10](q)), y = NULL, color = NULL, size = "Proteins", tag = "B") +
    ggplot2::theme_classic(base_size = 11, base_family = font_family) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(face = "bold", color = "#222222", size = 11.2),
      axis.text.y = ggplot2::element_text(face = "bold", color = "#222222", size = 9.2, lineheight = .84),
      axis.title = ggplot2::element_text(face = "bold", size = 12),
      legend.position = "top", legend.justification = "center", legend.box = "horizontal",
      legend.box.just = "center", legend.spacing.y = grid::unit(1.5, "pt"),
      legend.text = ggplot2::element_text(size = 8.8, face = "bold"),
      legend.key.width = grid::unit(6, "pt"),
      legend.key.height = grid::unit(7, "pt"),
      legend.spacing.x = grid::unit(1, "pt"),
      plot.tag = ggplot2::element_text(face = "bold", size = 16),
      plot.margin = ggplot2::margin(4, 8, 4, 4)
    )

  tf <- fread(tf_file)
  tf[, outcome_label_full := outcome_label]
  tf <- unique(tf, by = c("outcome_id", "outcome_label_full", "protein", "tf"))
  setnames(tf, "outcome_label_full", "outcome")

  outcome_order <- tf[, .N, by = .(outcome_id, outcome)][order(-N, outcome)]$outcome
  outcome_rank <- data.table(outcome = outcome_order, outcome_rank = seq_along(outcome_order))
  tf <- merge(tf, outcome_rank, by = "outcome", all.x = TRUE, sort = FALSE)
  protein_dominant <- tf[, .N, by = .(protein, outcome)][order(protein, -N, outcome)][, .SD[1L], by = protein]
  protein_total <- tf[, .(total = .N), by = protein]
  protein_dominant <- merge(protein_dominant, protein_total, by = "protein")
  protein_dominant[, outcome_order := match(outcome, outcome_order)]
  protein_order <- tf[, .(
    barycenter = mean(outcome_rank), total = .N
  ), by = protein][order(barycenter, -total, protein)]$protein
  protein_rank <- data.table(protein = protein_order, protein_rank = seq_along(protein_order))
  tf <- merge(tf, protein_rank, by = "protein", all.x = TRUE, sort = FALSE)
  tf_dominant <- tf[, .N, by = .(tf, outcome)][order(tf, -N, outcome)][, .SD[1L], by = tf]
  tf_total <- tf[, .(total = .N), by = tf]
  tf_dominant <- merge(tf_dominant, tf_total, by = "tf")
  tf_dominant[, outcome_order := match(outcome, outcome_order)]
  tf_order <- tf[, .(
    barycenter = mean(protein_rank), total = .N
  ), by = tf][order(barycenter, -total, tf)]$tf

  outcome_ids <- tf[, outcome_id[[1]], by = outcome]$V1
  outcome_colors <- yur_figure6_outcome_colors(outcome_ids)
  names(outcome_colors) <- tf[, outcome[[1]], by = outcome_id]$V1

  outcome_counts <- tf[, .(weight = .N), by = outcome][, .(label = outcome, weight)]
  protein_counts <- tf[, .(weight = .N), by = protein][, .(label = protein, weight)]
  tf_counts <- tf[, .(weight = .N), by = tf][, .(label = tf, weight)]
  outcome_nodes <- yur_sankey_node_layout(
    outcome_counts, outcome_order, x = 0, ymin = .08, ymax = .92, gap = .008
  )
  protein_nodes <- yur_sankey_node_layout(
    protein_counts, protein_order, x = 1.04, ymin = .035, ymax = .965, gap = .0035
  )
  tf_nodes <- yur_sankey_node_layout(
    tf_counts, tf_order, x = 2.08, ymin = .005, ymax = .995, gap = .0025
  )
  outcome_nodes[, fill_key := label]
  protein_nodes <- merge(protein_nodes, protein_dominant[, .(label = protein, fill_key = outcome)], by = "label", all.x = TRUE, sort = FALSE)
  tf_nodes <- merge(tf_nodes, tf_dominant[, .(label = tf, fill_key = outcome)], by = "label", all.x = TRUE, sort = FALSE)
  protein_nodes[, order_id := match(label, protein_order)]
  tf_nodes[, order_id := match(label, tf_order)]
  setorder(protein_nodes, order_id)
  setorder(tf_nodes, order_id)

  first_links <- tf[, .(weight = .N, edge_color = outcome[[1]]), by = .(outcome, protein)]
  first_links <- yur_sankey_allocate_links(first_links, outcome_nodes, protein_nodes, "outcome", "protein")
  second_links <- tf[, .(weight = 1, edge_color = outcome), by = .(outcome, protein, tf)]
  second_links <- yur_sankey_allocate_links(second_links, protein_nodes, tf_nodes, "protein", "tf")
  ribbon_data <- rbindlist(list(
    yur_sankey_ribbons(first_links, x0 = .016, x1 = 1.024, stage = "cvd_protein"),
    yur_sankey_ribbons(second_links, x0 = 1.056, x1 = 2.064, stage = "protein_tf")
  ))

  outcome_nodes[, label_display := vapply(label, function(value) paste(strwrap(value, width = 24), collapse = "\n"), character(1))]
  protein_nodes[, label_side := rep(c("left", "right"), length.out = .N)]
  tf_nodes[, label_side := rep(c("left", "right"), length.out = .N)]
  protein_nodes[, label_x := x + fifelse(label_side == "left", -.034, .034)]
  tf_nodes[, label_x := x + fifelse(label_side == "left", -.034, .034)]
  p6c <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = ribbon_data, ggplot2::aes(x, y, group = ribbon_group, fill = edge_color),
      color = NA, alpha = .38
    ) +
    ggplot2::geom_rect(
      data = outcome_nodes,
      ggplot2::aes(xmin = x - .016, xmax = x + .016, ymin = node_ymin, ymax = node_ymax, fill = fill_key),
      color = "white", linewidth = .22
    ) +
    ggplot2::geom_rect(
      data = protein_nodes,
      ggplot2::aes(xmin = x - .016, xmax = x + .016, ymin = node_ymin, ymax = node_ymax, fill = fill_key),
      color = "white", linewidth = .22
    ) +
    ggplot2::geom_rect(
      data = tf_nodes,
      ggplot2::aes(xmin = x - .016, xmax = x + .016, ymin = node_ymin, ymax = node_ymax, fill = fill_key),
      color = "white", linewidth = .22
    ) +
    ggplot2::geom_text(
      data = outcome_nodes,
      ggplot2::aes(x = x - .035, y = y, label = label_display),
      hjust = 1, family = font_family, fontface = "bold",
      size = 2.55, lineheight = .90, color = "#222222"
    ) +
    ggplot2::geom_text(
      data = protein_nodes[label_side == "left"], ggplot2::aes(label_x, y, label = label),
      hjust = 1, family = font_family, fontface = "bold", size = 2.15
    ) +
    ggplot2::geom_text(
      data = protein_nodes[label_side == "right"], ggplot2::aes(label_x, y, label = label),
      hjust = 0, family = font_family, fontface = "bold", size = 2.15
    ) +
    ggplot2::geom_text(
      data = tf_nodes[label_side == "left"], ggplot2::aes(label_x, y, label = label),
      hjust = 1, family = font_family, fontface = "bold", size = 1.90
    ) +
    ggplot2::geom_text(
      data = tf_nodes[label_side == "right"], ggplot2::aes(label_x, y, label = label),
      hjust = 0, family = font_family, fontface = "bold", size = 1.90
    ) +
    ggplot2::annotate("text", x = c(0, 1.04, 2.08), y = -0.035, label = c("CVDs", "Proteins", "TFs"),
                      family = font_family, fontface = "bold", size = 4.3) +
    ggplot2::scale_fill_manual(values = outcome_colors, guide = "none") +
    ggplot2::coord_cartesian(xlim = c(-.15, 2.46), ylim = c(-.06, 1.01), clip = "off") +
    ggplot2::labs(tag = "C") +
    ggplot2::theme_void(base_family = font_family) +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 16),
      plot.margin = ggplot2::margin(4, 20, 12, 35),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )

  nodes <- fread(node_file)
  edges <- fread(edge_file)
  graph <- igraph::graph_from_data_frame(edges[, .(from, to, score)], directed = FALSE, vertices = nodes)
  node_plot <- copy(nodes[order(-mnc, -degree)])
  inner_n <- min(3L, nrow(node_plot))
  middle1_n <- min(12L, nrow(node_plot))
  middle2_n <- min(24L, nrow(node_plot))
  node_plot[, ring := fifelse(
    seq_len(.N) <= inner_n, "Inner",
    fifelse(
      seq_len(.N) <= middle1_n, "Middle1",
      fifelse(seq_len(.N) <= middle2_n, "Middle2", "Outer")
    )
  )]
  inner <- yur_circle_positions(node_plot[ring == "Inner"]$gene_symbol, radius = .24, start = pi / 2)
  middle1 <- yur_circle_positions(node_plot[ring == "Middle1"]$gene_symbol, radius = .53, start = pi / 2)
  middle2 <- yur_circle_positions(node_plot[ring == "Middle2"]$gene_symbol, radius = .86, start = pi / 2)
  outer <- yur_circle_positions(node_plot[ring == "Outer"]$gene_symbol, radius = 1.22, start = pi / 2)
  node_plot <- merge(
    node_plot,
    rbindlist(list(inner, middle1, middle2, outer))[, .(gene_symbol = label, angle, x, y)],
    by = "gene_symbol"
  )
  node_plot[, label_color := fifelse(mnc >= quantile(mnc, .72, na.rm = TRUE), "white", "black")]
  edge_ends <- igraph::ends(graph, igraph::E(graph), names = TRUE)
  edge_plot <- data.table(from = edge_ends[, 1], to = edge_ends[, 2], score = igraph::E(graph)$score)
  edge_plot <- merge(edge_plot, node_plot[, .(from = gene_symbol, x, y)], by = "from")
  setnames(edge_plot, c("x", "y"), c("x_from", "y_from"))
  edge_plot <- merge(edge_plot, node_plot[, .(to = gene_symbol, x, y)], by = "to")
  setnames(edge_plot, c("x", "y"), c("x_to", "y_to"))
  p6d <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_plot,
      ggplot2::aes(x = x_from, y = y_from, xend = x_to, yend = y_to, alpha = score, linewidth = score),
      color = "#667078", lineend = "round"
    ) +
    ggplot2::geom_point(
      data = node_plot,
      ggplot2::aes(x, y, size = degree, fill = mnc), shape = 21, color = "white", stroke = .55, alpha = .99
    ) +
    ggplot2::geom_text(
      data = node_plot[label_color == "black"], ggplot2::aes(x, y, label = gene_symbol),
      family = font_family, fontface = "bold", size = 2.05, color = "black", check_overlap = FALSE
    ) +
    ggplot2::geom_text(
      data = node_plot[label_color == "white"], ggplot2::aes(x, y, label = gene_symbol),
      family = font_family, fontface = "bold", size = 2.05, color = "white", check_overlap = FALSE
    ) +
    ggplot2::scale_fill_gradientn(colors = c("#F6D891", "#EEA246", "#E45B36", "#B91F2D")) +
    ggplot2::scale_size_continuous(range = c(5.8, 16.2)) +
    ggplot2::scale_alpha_continuous(range = c(.04, .30), guide = "none") +
    ggplot2::scale_linewidth_continuous(range = c(.15, .65), guide = "none") +
    ggplot2::coord_equal(xlim = c(-1.28, 1.28), ylim = c(-1.28, 1.28), clip = "off") +
    ggplot2::labs(tag = "D", size = "Degree", fill = "MNC") +
    ggplot2::theme_void(base_family = font_family) +
    ggplot2::theme(
      legend.position = "none", plot.tag = ggplot2::element_text(face = "bold", size = 16),
      plot.margin = ggplot2::margin(2, 2, 2, 2),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )

  regulatory_edges <- unique(tf[, .(protein, tf)])
  regulatory_proteins <- tf[, .N, by = protein][order(-N)]$protein
  regulatory_tfs <- tf[, .N, by = tf][order(-N)]$tf
  regulatory_inner <- yur_circle_positions(regulatory_proteins, radius = .78, start = pi / 2)
  regulatory_outer <- yur_circle_positions(regulatory_tfs, radius = 1.62, start = pi / 2)
  regulatory_edges <- merge(regulatory_edges, regulatory_inner[, .(protein = label, x0 = x, y0 = y)], by = "protein")
  regulatory_edges <- merge(regulatory_edges, regulatory_outer[, .(tf = label, x1 = x, y1 = y)], by = "tf")
  regulatory_outer[, `:=`(
    label_x = 1.84 * cos(angle), label_y = 1.84 * sin(angle),
    label_angle = angle * 180 / pi, hjust = fifelse(cos(angle) < 0, 1, 0)
  )]
  regulatory_outer[label_angle < -90, label_angle := label_angle + 180]
  p6c_circle <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = regulatory_edges, ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1),
      color = "#8E9698", alpha = .18, linewidth = .28
    ) +
    ggplot2::geom_point(data = regulatory_inner, ggplot2::aes(x, y), shape = 21, size = 6,
                        fill = "#D95F45", color = "white", stroke = .45) +
    ggplot2::geom_text(data = regulatory_inner, ggplot2::aes(x, y, label = label),
                       family = font_family, fontface = "bold", size = 2.25) +
    ggplot2::geom_point(data = regulatory_outer, ggplot2::aes(x, y), shape = 21, size = 3.8,
                        fill = "#EAB654", color = "white", stroke = .4) +
    ggplot2::geom_text(
      data = regulatory_outer,
      ggplot2::aes(label_x, label_y, label = label, angle = label_angle, hjust = hjust),
      family = font_family, fontface = "bold", size = 2.4, color = "#333333"
    ) +
    ggplot2::coord_equal(xlim = c(-2.25, 2.25), ylim = c(-2.25, 2.25), clip = "off") +
    ggplot2::theme_void(base_family = font_family) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA))

  figure6a_raster <- yur_trimmed_png_grob(figure6a_file)
  p6a <- ggplot2::ggplot() +
    ggplot2::annotation_custom(figure6a_raster$grob, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    ggplot2::labs(tag = "A") +
    ggplot2::theme_void(base_family = font_family) +
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 16),
                   plot.margin = ggplot2::margin(2, 2, 2, 2),
                   plot.background = ggplot2::element_rect(fill = "white", color = NA))

  yur_save_plot(p6b, "figure6b_local_cvd_protein_enrichment", cfg, 7.2, 7.6)
  yur_save_plot(p6c, "figure6c_local_cvd_protein_tf_sankey", cfg, 7.8, 12.2)
  yur_save_plot(p6d, "figure6d_local_string_ppi_mnc", cfg, 8.2, 7.4)
  yur_save_plot(p6c_circle, "figure6c_supp_local_protein_tf_concentric", cfg, 9.0, 9.0)
  b_row <- (p6b | patchwork::plot_spacer()) + patchwork::plot_layout(widths = c(.527, .473))
  c_row <- patchwork::wrap_elements(full = p6c)
  d_row <- (patchwork::wrap_elements(full = p6d) | patchwork::plot_spacer()) +
    patchwork::plot_layout(widths = c(.878, .122))
  right_panel <- b_row / c_row / d_row + patchwork::plot_layout(heights = c(.31, .37, .32))
  combined <- (p6a | right_panel) + patchwork::plot_layout(widths = c(3.96, 6.84)) &
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA))
  yur_save_plot(combined, "figure6abcd_local_systems_biology", cfg, 10.8, 18.0)
  yur_save_plot(combined, "figure6bcd_local_systems_biology", cfg, 10.8, 18.0)

  yur_write_csv(top_terms, file.path(cfg$paths$figures, "figure6b_local_source_data.csv"))
  yur_write_csv(tf, file.path(cfg$paths$figures, "figure6c_local_source_data.csv.gz"))
  yur_write_csv(node_plot, file.path(cfg$paths$figures, "figure6d_local_node_source_data.csv"))
  yur_write_csv(edge_plot, file.path(cfg$paths$figures, "figure6d_local_edge_source_data.csv.gz"))
  yur_write_csv(regulatory_inner, file.path(cfg$paths$figures, "figure6c_supp_protein_nodes.csv"))
  yur_write_csv(regulatory_outer, file.path(cfg$paths$figures, "figure6c_supp_tf_nodes.csv"))
  yur_write_csv(regulatory_edges, file.path(cfg$paths$figures, "figure6c_supp_protein_tf_edges.csv.gz"))
  yur_write_json(list(
    status = "PASS", evidence = "LOCAL_COX_PLUS_VERSION_LOCKED_EXTERNAL_ANNOTATION",
    figure6a = "Local PRS-protein coefficient heat map included as the full-height left column",
    figure6b = "MSigDB 2026.1 over-representation analysis with local measured-protein background",
    figure6c = "Local per-outcome top-15 Cox proteins mapped to downloaded TRRUST human TF-target edges",
    figure6d = "STRING v12 functional PPI, MCODE-compatible cluster 1, MNC hub score",
    figure6c_supplement = "Concentric regulatory network with proteins in the inner ring and TFs in the outer ring",
    caution = "Numerically reproducible for the frozen API/database versions; not a reuse of the article's unavailable Metascape/Cytoscape sessions."
  ), file.path(out_dir, "figure6bcd_manifest.json"))
}

yur_run_figure6_systems_all <- function(cfg) {
  yur_prepare_figure6_systems(cfg)
  yur_run_figure6_enrichment(cfg)
  yur_run_figure6_tf(cfg)
  yur_run_figure6_ppi(cfg)
  yur_plot_figure6_systems(cfg)
}
