# GO enrichment using the local Cox-significant proteins.
# This follows the paper's enrichment logic at a runnable local level:
# foreground = significant disease-associated proteins; universe = tested proteins.

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
repro_load_packages(c("data.table", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "ggplot2", "ggtext", "ggpubr", "gtools", "openxlsx"))

message("Step 06: GO enrichment")

enrich_dir <- file.path(output_dir, "enrichment")
fig_dir <- file.path(output_dir, "figures")
dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

primary <- data.table::fread(primary_assoc_file)
repro_assert_primary_complete(primary, audit_dir, outcome_map, context = "Primary association results used for GO enrichment")
primary <- primary[is.finite(P_Value)]
outcome_levels <- c("CAD", "HF", "Afib", "AS")
outcome_full <- c(CAD = "Coronary artery disease", HF = "Heart failure", Afib = "Atrial fibrillation", AS = "Aortic stenosis")
domain_labels <- c(BP = "Biological process", MF = "Molecular function", CC = "Cellular component")
capitalize <- function(x) gsub("\\b([[:lower:]])", "\\U\\1", x, perl = TRUE)

clean_symbol <- function(x) {
  x <- as.character(x)
  x <- ifelse(x == "NTPROBNP", "NPPB", x)
  x <- ifelse(grepl("^HLA_", x), gsub("_", "-", x), x)
  x
}

map_symbols <- function(symbols) {
  data.table::data.table(Protein = unique(symbols), SYMBOL_QUERY = clean_symbol(unique(symbols)))[
    , {
      hit_symbol <- suppressMessages(AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = unique(SYMBOL_QUERY),
        keytype = "SYMBOL",
        columns = c("SYMBOL", "ENTREZID")
      ))
      hit_alias <- suppressMessages(AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = setdiff(unique(SYMBOL_QUERY), hit_symbol$SYMBOL[!is.na(hit_symbol$ENTREZID)]),
        keytype = "ALIAS",
        columns = c("SYMBOL", "ENTREZID")
      ))
      data.table::rbindlist(list(
        data.table::as.data.table(hit_symbol)[, .(SYMBOL_QUERY = SYMBOL, SYMBOL, ENTREZID)],
        data.table::as.data.table(hit_alias)[, .(SYMBOL_QUERY = ALIAS, SYMBOL, ENTREZID)]
      ), fill = TRUE)
    }
  ][!is.na(ENTREZID)][!duplicated(SYMBOL_QUERY)]
}

symbol_map <- merge(
  data.table::data.table(Protein = unique(primary$Protein), SYMBOL_QUERY = clean_symbol(unique(primary$Protein))),
  map_symbols(unique(primary$Protein)),
  by = "SYMBOL_QUERY",
  all.x = TRUE
)
data.table::fwrite(symbol_map, file.path(enrich_dir, "protein_entrez_mapping.csv"))

universe_entrez <- unique(symbol_map[!is.na(ENTREZID), ENTREZID])
enrich_rows <- list()

for (outc in outcome_levels) {
  foreground <- unique(primary[Outcome == outc & Bonferroni %in% TRUE, Protein])
  foreground_rule <- "Bonferroni"
  if (length(foreground) < 5) {
    foreground_fdr <- unique(primary[Outcome == outc & FDR < 0.05, Protein])
    if (length(foreground_fdr) >= 5) {
      foreground <- foreground_fdr
      foreground_rule <- "FDR < 0.05 fallback because Bonferroni foreground had fewer than 5 mapped proteins"
    }
  }
  foreground_map <- symbol_map[Protein %in% foreground & !is.na(ENTREZID)]
  if (nrow(foreground_map) < 5) {
    next
  }
  for (ont in c("BP", "MF", "CC")) {
    ego <- tryCatch(
      clusterProfiler::enrichGO(
        gene = unique(foreground_map$ENTREZID),
        universe = universe_entrez,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        readable = TRUE
      ),
      error = function(e) NULL
    )
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      d <- data.table::as.data.table(as.data.frame(ego))
      d[, `:=`(
        Outcome = outc,
        Disease = outcome_full[[outc]],
        Domain = domain_labels[[ont]],
        Ontology = ont,
        Foreground_Rule = foreground_rule,
        Foreground_N = length(unique(foreground_map$Protein)),
        Number = Count,
        Ont = Description,
        P = pvalue,
        FDR_P = p.adjust
      )]
      enrich_rows[[paste(outc, ont, sep = "_")]] <- d
    }
  }
}

go_all <- data.table::rbindlist(enrich_rows, fill = TRUE)
if (nrow(go_all) == 0) {
  stop("No GO enrichment rows were produced. Check significant protein mapping.", call. = FALSE)
}
data.table::fwrite(go_all, file.path(enrich_dir, "go_enrichment_all_terms.csv"))
openxlsx::write.xlsx(go_all, file.path(enrich_dir, "2_GO.local_enrichment.xlsx"), overwrite = TRUE)

top_terms_per_domain <- max(2L, repro_int_env("UKBPPP_GO_TOP_PER_DOMAIN", 4L))
plot_go <- go_all[order(FDR_P, P), head(.SD, top_terms_per_domain), by = .(Disease, Domain)]
plot_go[, Domain := factor(Domain, levels = c("Biological process", "Molecular function", "Cellular component"))]
plot_go[, Ont_clean := vapply(strwrap(capitalize(tolower(sub("\\s*\\(GO:[^)]+\\)", "", Ont))), 38, simplify = FALSE), paste, character(1), collapse = "\n")]
plot_go[, alpha := ifelse(FDR_P < 0.05, 0.9, 0.8)]
data.table::fwrite(plot_go[, .(Disease, Domain, Foreground_Rule, Foreground_N, Ont, Number, P, FDR_P)], file.path(enrich_dir, "go_enrichment_plot_terms.csv"))

make_go_plot <- function(disease) {
  d <- plot_go[Disease == disease]
  d <- d[order(Domain, FDR_P, P)]
  d[, term_plot := factor(Ont_clean, levels = rev(unique(Ont_clean)))]
  ggplot2::ggplot(d, ggplot2::aes(x = term_plot, y = -log10(P), fill = Domain)) +
    ggplot2::geom_col(ggplot2::aes(alpha = alpha)) +
    ggplot2::scale_alpha_continuous(guide = "none", range = c(0.5, 0.9)) +
    ggplot2::scale_fill_manual(
      values = c("Biological process" = "#F39B7F", "Molecular function" = "#B97E64", "Cellular component" = "#7E6148"),
      name = "GO domain",
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(limits = c(0, max(6, ceiling(max(-log10(plot_go$P), na.rm = TRUE)))), expand = ggplot2::expansion(mult = c(0, 0.04))) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = disease, y = "-log<sub>10</sub>(P)", x = NULL) +
    ggplot2::coord_flip() +
    ggplot2::theme(
      axis.title.y = ggtext::element_markdown(),
      axis.title.x = ggtext::element_markdown(),
      legend.position = "bottom",
      axis.text.y = ggplot2::element_text(size = 7),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 11)
    )
}

plots <- lapply(outcome_full[outcome_levels], make_go_plot)
names(plots) <- outcome_levels
for (outc in names(plots)) {
  ggplot2::ggsave(
    file.path(fig_dir, paste0("original_style_go_enrichment_", tolower(outc), ".png")),
    plots[[outc]], width = 7.5, height = 6.5, dpi = 320, bg = "white"
  )
}
combined <- ggpubr::ggarrange(
  plots[["CAD"]], plots[["HF"]], plots[["Afib"]], plots[["AS"]],
  common.legend = TRUE, ncol = 2, nrow = 2, legend = "bottom"
)
ggplot2::ggsave(file.path(fig_dir, "original_style_go_enrichment_all.png"), combined, width = 15, height = 12, dpi = 320, bg = "white")
ggplot2::ggsave(file.path(enrich_dir, "2_GO.local_enrichment.png"), combined, width = 15, height = 12, dpi = 320, bg = "white")

required_go_figures <- c(
  "original_style_go_enrichment_all.png",
  "original_style_go_enrichment_cad.png",
  "original_style_go_enrichment_hf.png",
  "original_style_go_enrichment_afib.png",
  "original_style_go_enrichment_as.png"
)
go_figure_qc <- data.table::data.table(
  file = required_go_figures,
  exists = file.exists(file.path(fig_dir, required_go_figures)),
  path = file.path(fig_dir, required_go_figures)
)
data.table::fwrite(go_figure_qc, file.path(fig_dir, "go_figure_qc.csv"))
if (any(!go_figure_qc$exists)) {
  stop(
    "GO_FIGURE_GATE_FAILED: missing ",
    paste(go_figure_qc[exists == FALSE, file], collapse = ", "),
    call. = FALSE
  )
}

# Step 06 creates the enrichment figures after Step 05, so refresh the final
# manifest here to keep the website/report inventory complete.
figure_manifest <- data.table::data.table(
  file = list.files(fig_dir, pattern = "\\.(png|csv)$", full.names = FALSE),
  path = list.files(fig_dir, pattern = "\\.(png|csv)$", full.names = TRUE)
)
data.table::fwrite(figure_manifest, file.path(fig_dir, "figure_manifest.csv"))

message("GO enrichment outputs saved to: ", enrich_dir)
