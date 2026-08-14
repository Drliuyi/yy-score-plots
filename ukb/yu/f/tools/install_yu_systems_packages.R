required <- c(
  "data.table", "jsonlite", "digest", "httr2", "curl", "igraph",
  "ggplot2", "ggalluvial", "ggrepel", "patchwork", "scales", "msigdbr"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing missing Figure 6 systems packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org", dependencies = TRUE)
}
still_missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing)) {
  stop("Required Figure 6 systems packages remain unavailable: ", paste(still_missing, collapse = ", "))
}
bioconductor <- c("AnnotationDbi", "org.Hs.eg.db", "GO.db")
missing_bioc <- bioconductor[!vapply(bioconductor, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}
still_missing_bioc <- bioconductor[!vapply(bioconductor, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing_bioc)) {
  stop("Required Bioconductor packages remain unavailable: ", paste(still_missing_bioc, collapse = ", "))
}
cat("FIGURE 6 SYSTEMS PACKAGE CHECK PASSED\n")
