# Time-course clustering for one treatment series
# (port of 05_atac_temporal_clustering 00_preprocessing + 01_degpatterns).

source(snakemake@params[["helpers"]])
start_log()
suppressPackageStartupMessages({
  library(DESeq2)
  library(DEGreport)
  library(ggplot2)
})
p <- snakemake@params
series <- p$series

counts <- read.delim(snakemake@input[["counts"]], row.names = 1, check.names = FALSE)
meta <- data.frame(sample = unlist(p$libs), treatment = unlist(p$treatments),
                   time = as.numeric(unlist(p$times)),
                   condition = unlist(p$conditions), stringsAsFactors = FALSE)
rownames(meta) <- meta$sample
counts <- as.matrix(counts[, meta$sample])

# 00_preprocessing: DESeq2 object on all included samples, mean-count filter,
# blind VST.
dds_all <- DESeqDataSetFromMatrix(counts, colData = meta, design = ~ condition)
keep <- rowMeans(counts(dds_all)) >= p$min_mean_counts
dds_all <- dds_all[keep, ]
cat(sprintf("Peaks after rowMeans >= %s: %d / %d\n", p$min_mean_counts, sum(keep), length(keep)))
# vst() subsamples 1000 rows to fit the dispersion trend; fall back to the
# full transformation for small peak sets (e.g. the test data).
vst_all <- if (nrow(dds_all) >= 1000) {
  assay(vst(dds_all, blind = TRUE))
} else {
  assay(varianceStabilizingTransformation(dds_all, blind = TRUE))
}

# Series = this treatment + shared baseline samples
idx <- meta$treatment == series | meta$treatment %in% unlist(p$baseline)
meta_s <- meta[idx, ]
times <- sort(unique(meta_s$time))
cat("Series", series, ":", nrow(meta_s), "samples; times =", paste(times, collapse = ", "), "\n")
if (length(times) < 2) stop("Series ", series, " has fewer than two time points")
meta_s$time_f <- factor(meta_s$time, levels = times)

dds <- DESeqDataSetFromMatrix(counts(dds_all)[, meta_s$sample],
                              colData = meta_s, design = ~ time_f)
dds <- DESeq(dds, test = "LRT", reduced = ~ 1, parallel = FALSE)
res <- results(dds)

vst_s <- vst_all[, meta_s$sample, drop = FALSE]
tmeans <- sapply(split(seq_len(ncol(vst_s)), meta_s$time_f),
                 function(i) rowMeans(vst_s[, i, drop = FALSE]))
vrange <- apply(tmeans, 1, function(x) max(x) - min(x))
sig <- !is.na(res$padj) & res$padj < p$lrt_padj
dynamic <- sig & vrange >= p$vst_range
cat(sprintf("LRT padj < %s: %d; with VST range >= %s: %d dynamic peaks\n",
            p$lrt_padj, sum(sig), p$vst_range, sum(dynamic)))

lrt <- data.frame(peak_id = rownames(res), baseMean = res$baseMean, stat = res$stat,
                  pvalue = res$pvalue, padj = res$padj, vst_range = vrange,
                  dynamic = dynamic)
write_tsv(lrt, snakemake@output[["lrt"]])

vst_dyn <- vst_s[dynamic, , drop = FALSE]
vst_dyn <- vst_dyn[apply(vst_dyn, 1, var, na.rm = TRUE) > 0, , drop = FALSE]
md <- data.frame(time = meta_s$time_f, row.names = meta_s$sample)

set.seed(p$seed)
degp <- NULL
if (nrow(vst_dyn) >= p$minc) {
  degp <- tryCatch(
    degPatterns(vst_dyn, metadata = md, time = "time", col = NULL,
                minc = p$minc, cutoff = p$cutoff, plot = FALSE),
    error = function(e) { message("degPatterns failed: ", conditionMessage(e)); NULL })
} else {
  message(sprintf("Only %d dynamic peaks (< minc = %s); skipping degPatterns",
                  nrow(vst_dyn), p$minc))
}

if (!is.null(degp)) {
  cl <- unique(degp$df[, c("genes", "cluster")])
  colnames(cl) <- c("peak_id", "cluster")
  norm <- degp$normalized
  prof <- aggregate(value ~ cluster + time, data = norm, FUN = mean)
  g <- ggplot(norm, aes(time, value)) +
    geom_boxplot(outlier.shape = NA, fill = "grey90") +
    stat_summary(aes(group = 1), fun = mean, geom = "line", colour = "#1B7837") +
    facet_wrap(~ cluster, scales = "free_y") +
    labs(x = "Time", y = "Z-score (VST)", title = paste("degPatterns clusters:", series)) +
    theme_bw()
} else {
  cl <- data.frame(peak_id = character(), cluster = integer())
  prof <- data.frame(cluster = integer(), time = character(), value = numeric())
  g <- ggplot() + annotate("text", x = 0, y = 0, label = "No clusters") + theme_void()
}
write_tsv(cl, snakemake@output[["clusters"]])
write_tsv(prof, snakemake@output[["profiles"]])
ggsave(snakemake@output[["plot"]], g, width = 12, height = 8)
saveRDS(list(lrt = lrt, degpatterns = degp, meta = meta_s), snakemake@output[["rds"]])
cat("\n"); print(sessionInfo())
