# Export JASPAR2020 PFMs (default CORE vertebrates, non-redundant latest
# versions) in JASPAR format for TOBIAS, so footprinting and chromVAR use the
# same motif set.

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")
suppressPackageStartupMessages({
  library(JASPAR2020)
  library(TFBSTools)
})

pfms <- getMatrixSet(JASPAR2020, opts = list(
  collection = snakemake@params[["collection"]],
  tax_group = snakemake@params[["tax_group"]],
  all_versions = FALSE))
cat("Exporting", length(pfms), "matrices\n")

con <- file(snakemake@output[[1]], "w")
for (i in seq_along(pfms)) {
  m <- Matrix(pfms[[i]])
  writeLines(sprintf(">%s\t%s", ID(pfms[[i]]), name(pfms[[i]])), con)
  for (b in c("A", "C", "G", "T")) {
    writeLines(sprintf("%s  [ %s ]", b, paste(sprintf("%6d", as.integer(m[b, ])), collapse = " ")), con)
  }
}
close(con)
