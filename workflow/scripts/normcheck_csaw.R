# Re-run the diff contrasts under csaw background-bin normalization
# (port of 08b NB01). Reuses the counted DBA; only normalize + analyze re-run.

source(snakemake@params[["helpers"]])
start_log()
suppressPackageStartupMessages({
  library(DiffBind)
  library(csaw)
})

fdr <- as.numeric(snakemake@params[["fdr"]])
lfc <- as.numeric(snakemake@params[["lfc"]])
tab_dir <- snakemake@params[["tabdir"]]
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

dba_obj <- readRDS(snakemake@input[["dba"]])
dba_obj$config$cores <- snakemake@threads
ct <- read_contrasts(snakemake@input[["contrasts"]])

dba_csaw <- dba.normalize(dba_obj, method = DBA_DESEQ2,
                          normalize = DBA_NORM_NATIVE, background = TRUE)
sf <- dba.normalize(dba_csaw, bRetrieve = TRUE)
sf_tbl <- data.frame(
  SampleID = dba_csaw$samples$SampleID,
  Condition = dba_csaw$samples$Condition,
  Replicate = dba_csaw$samples$Replicate,
  lib_size = sf$lib.sizes,
  csaw_norm_factor = sf$norm.factors
)
write_tsv(sf_tbl, snakemake@output[["sizefactors"]])
print(sf_tbl)

dba_csaw <- add_contrasts(dba_csaw, ct)
dba_csaw <- dba.analyze(dba_csaw)

summaries <- list()
for (i in seq_len(nrow(ct))) {
  df <- report_df(dba_csaw, i)
  write_tsv(df, file.path(tab_dir, paste0(ct$label[i], "_all.tsv")))
  summaries[[i]] <- summarise_contrast(df, ct$label[i], fdr, lfc)
}
summary_df <- do.call(rbind, summaries)
summary_df$method <- "csaw"
write_tsv(summary_df, snakemake@output[["summary"]])
print(summary_df)

saveRDS(dba_csaw, snakemake@output[["rds"]])
cat("\n"); print(sessionInfo())
