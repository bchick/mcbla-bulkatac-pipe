# DiffBind differential accessibility under DiffBind's default normalization
# (port of 08 NB01: dba.normalize -> contrasts -> dba.analyze -> dba.report).

source(snakemake@params[["helpers"]])
start_log()
suppressPackageStartupMessages({
  library(DiffBind)
  library(ggplot2)
})

fdr <- as.numeric(snakemake@params[["fdr"]])
lfc <- as.numeric(snakemake@params[["lfc"]])
tab_dir <- snakemake@params[["tabdir"]]
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

dba_obj <- readRDS(snakemake@input[["dba"]])
dba_obj$config$cores <- snakemake@threads
ct <- read_contrasts(snakemake@input[["contrasts"]])

dba_obj <- dba.normalize(dba_obj)
print(dba.normalize(dba_obj, bRetrieve = TRUE))

pdf(snakemake@output[["heatmap"]], width = 10, height = 9)
dba.plotHeatmap(dba_obj, main = "ATAC-seq sample correlation")
dev.off()
pdf(snakemake@output[["pca"]], width = 9, height = 6)
dba.plotPCA(dba_obj, label = DBA_CONDITION, attributes = DBA_CONDITION)
dev.off()

dba_obj <- add_contrasts(dba_obj, ct)
dba_obj <- dba.analyze(dba_obj)
print(dba.show(dba_obj, bContrasts = TRUE))

summaries <- list()
for (i in seq_len(nrow(ct))) {
  label <- ct$label[i]
  df <- report_df(dba_obj, i)
  write_tsv(df, file.path(tab_dir, paste0(label, "_all.tsv")))
  write_tsv(df[!is.na(df$FDR) & df$FDR < fdr, ], file.path(tab_dir, paste0(label, "_sig.tsv")))
  summaries[[label]] <- summarise_contrast(df, label, fdr, lfc)
}
summary_df <- do.call(rbind, summaries)
summary_df$method <- "depth"
write_tsv(summary_df, snakemake@output[["summary"]])
print(summary_df)

pdf(snakemake@output[["ma"]], width = 7, height = 6)
for (i in seq_len(nrow(ct))) dba.plotMA(dba_obj, contrast = i, sub = ct$label[i])
dev.off()

long <- rbind(
  data.frame(Contrast = summary_df$Contrast, Direction = "Gained", N = summary_df$Gained),
  data.frame(Contrast = summary_df$Contrast, Direction = "Lost", N = -summary_df$Lost)
)
long$Contrast <- factor(long$Contrast, levels = rev(ct$label))
p <- ggplot(long, aes(Contrast, N, fill = Direction)) +
  geom_col() + geom_hline(yintercept = 0) + coord_flip() +
  scale_fill_manual(values = c(Gained = "#1B7837", Lost = "#762A83")) +
  labs(x = NULL, y = sprintf("Differential peaks (FDR < %g)", fdr),
       title = "Differential accessibility (DiffBind default normalization)") +
  theme_bw()
ggsave(snakemake@output[["barplot"]], p, width = 9, height = 1.5 + 0.35 * nrow(ct))

saveRDS(dba_obj, snakemake@output[["rds"]])
cat("\n"); print(sessionInfo())
