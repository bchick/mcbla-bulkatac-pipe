# Re-run the diff contrasts under quantile normalization + limma
# (port of 08b NB02).
#   raw reads per consensus peak -> edgeR::filterByExpr -> log2(count + 0.5)
#   -> preprocessCore::normalize.quantiles -> lmFit(~0 + condition)
#   -> contrasts.fit -> eBayes -> topTable
# Deviation from the notebook: counts are retrieved explicitly as raw reads
# (DBA_SCORE_READS). Quantile normalization removes per-sample scale, so this
# changes little, but it makes the input unambiguous.

source(snakemake@params[["helpers"]])
start_log()
suppressPackageStartupMessages({
  library(DiffBind)
  library(limma)
  library(edgeR)
  library(preprocessCore)
})

fdr <- as.numeric(snakemake@params[["fdr"]])
lfc <- as.numeric(snakemake@params[["lfc"]])
tab_dir <- snakemake@params[["tabdir"]]
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

dba_obj <- readRDS(snakemake@input[["dba"]])
ct <- read_contrasts(snakemake@input[["contrasts"]])

dba_obj <- dba.count(dba_obj, peaks = NULL, score = DBA_SCORE_READS)
gr <- dba.peakset(dba_obj, bRetrieve = TRUE, DataType = DBA_DATA_GRANGES)
ids <- dba_obj$samples$SampleID
count_mat <- as.matrix(GenomicRanges::mcols(gr)[, ids])
rownames(count_mat) <- sprintf("%s:%d-%d", GenomicRanges::seqnames(gr),
                               GenomicRanges::start(gr), GenomicRanges::end(gr))
cond <- dba_obj$samples$Condition

keep <- filterByExpr(DGEList(counts = count_mat, group = cond), group = cond)
cat("Keep", sum(keep), "/", length(keep), "peaks after filterByExpr\n")
log_counts <- log2(count_mat[keep, , drop = FALSE] + 0.5)
# preprocessCore (as in 08b NB02). Some conda builds of preprocessCore fail
# with "return code from pthread_create() is 22"; fall back to
# limma::normalizeQuantiles, the same algorithm with ties averaged.
qn <- tryCatch(
  normalize.quantiles(log_counts, keep.names = TRUE),
  error = function(e) {
    message("preprocessCore::normalize.quantiles failed (", conditionMessage(e),
            "); using limma::normalizeQuantiles")
    limma::normalizeQuantiles(log_counts)
  })

write_tsv(data.frame(sample = ids,
                     mean_pre = colMeans(log_counts), mean_post = colMeans(qn),
                     sd_pre = apply(log_counts, 2, sd), sd_post = apply(qn, 2, sd)),
          snakemake@output[["diagnostics"]])

safe <- make.names(unique(cond))
names(safe) <- unique(cond)
f <- factor(safe[cond], levels = safe)
design <- model.matrix(~ 0 + f)
colnames(design) <- levels(f)
fit <- lmFit(qn, design)
cm <- makeContrasts(contrasts = paste0(safe[ct$group1], "-", safe[ct$group2]),
                    levels = design)
colnames(cm) <- ct$label
fit2 <- eBayes(contrasts.fit(fit, cm))

summaries <- list()
for (i in seq_len(nrow(ct))) {
  tt <- topTable(fit2, coef = ct$label[i], number = Inf, sort.by = "none")
  parts <- do.call(rbind, regmatches(rownames(tt),
                   regexec("^(.*):([0-9]+)-([0-9]+)$", rownames(tt))))
  df <- data.frame(seqnames = parts[, 2], start = as.integer(parts[, 3]),
                   end = as.integer(parts[, 4]), Fold = tt$logFC,
                   AveExpr = tt$AveExpr, t = tt$t, p.value = tt$P.Value,
                   FDR = tt$adj.P.Val, B = tt$B, peak_id = rownames(tt))
  write_tsv(df, file.path(tab_dir, paste0(ct$label[i], "_all.tsv")))
  summaries[[i]] <- summarise_contrast(df, ct$label[i], fdr, lfc)
}
summary_df <- do.call(rbind, summaries)
summary_df$method <- "quantile"
write_tsv(summary_df, snakemake@output[["summary"]])
print(summary_df)
cat("\n"); print(sessionInfo())
