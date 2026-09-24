# chromVAR TF deviation analysis (port of 05/07_chromvar.Rmd and the motif
# matching in 05/05_chisq_motif_enrichment.Rmd).

source(snakemake@params[["helpers"]])
start_log()
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(Rsamtools)
  library(chromVAR)
  library(motifmatchr)
  library(JASPAR2020)
  library(TFBSTools)
  library(BiocParallel)
})
p <- snakemake@params
register(MulticoreParam(snakemake@threads))

counts <- read.delim(snakemake@input[["counts"]], row.names = 1, check.names = FALSE)
libs <- unlist(p$libs)
counts <- as.matrix(counts[, libs])
m <- regmatches(rownames(counts), regexec("^(.*):([0-9]+)-([0-9]+)$", rownames(counts)))
m <- do.call(rbind, m)
gr <- GRanges(m[, 2], IRanges(as.integer(m[, 3]), as.integer(m[, 4])))
gr$peak_id <- rownames(counts)

if (nzchar(p$bsgenome)) {
  suppressPackageStartupMessages(library(p$bsgenome, character.only = TRUE))
  genome <- get(p$bsgenome)
  cat("Genome: BSgenome", p$bsgenome, "\n")
} else {
  genome <- FaFile(snakemake@input[["fasta"]])
  cat("Genome: FASTA", snakemake@input[["fasta"]], "\n")
}
seqs <- if (is(genome, "FaFile")) as.character(seqnames(scanFaIndex(genome))) else seqnames(genome)
ok <- as.character(seqnames(gr)) %in% seqs
gr <- gr[ok]; counts <- counts[ok, , drop = FALSE]

se <- SummarizedExperiment(assays = list(counts = counts), rowRanges = gr,
                           colData = DataFrame(sample = libs,
                                               condition = unlist(p$conditions)))
se <- se[rowSums(assay(se, "counts")) >= p$min_fragments, ]
cat(sprintf("Peaks with >= %s fragments: %d\n", p$min_fragments, nrow(se)))

set.seed(p$seed)
se <- addGCBias(se, genome = genome)

pwms <- getMatrixSet(JASPAR2020, opts = list(collection = p$collection,
                                             tax_group = p$tax_group,
                                             all_versions = FALSE))
cat("JASPAR2020 motifs:", length(pwms), "\n")
windows <- resize(rowRanges(se), width = p$peak_width, fix = "center")
motif_ix <- matchMotifs(pwms, windows, genome = genome, out = "matches")
# chromVAR needs annotations aligned to the counts object
motif_se <- SummarizedExperiment(assays = list(motifMatches = motifMatches(motif_ix)),
                                 rowRanges = rowRanges(se), colData = colData(motif_ix))

set.seed(p$seed)
bg <- getBackgroundPeaks(se, niterations = p$iterations, w = 0.1)
dev <- computeDeviations(object = se, annotations = motif_se, background_peaks = bg)

name_map <- setNames(TFBSTools::name(pwms), TFBSTools::ID(pwms))
label <- paste0(rownames(dev), "_", name_map[rownames(dev)])
z <- assay(dev, "z"); d <- assay(dev, "deviations")
write_tsv(data.frame(motif = label, z, check.names = FALSE), snakemake@output[["z"]])
write_tsv(data.frame(motif = label, d, check.names = FALSE), snakemake@output[["dev"]])
v <- computeVariability(dev)
v$motif <- paste0(rownames(v), "_", name_map[rownames(v)])
v <- v[order(-v$variability), ]
write_tsv(v, snakemake@output[["var"]])

top <- head(rownames(v), p$n_top)
zt <- z[top, , drop = FALSE]
rownames(zt) <- v$motif[match(top, rownames(v))]
ord <- order(unlist(p$conditions))
pdf(snakemake@output[["heatmap"]], width = 9, height = max(6, 0.18 * nrow(zt) + 3))
heatmap(zt[, ord, drop = FALSE], Colv = NA, scale = "none",
        col = colorRampPalette(c("#762A83", "white", "#1B7837"))(101),
        margins = c(10, 14), main = "chromVAR deviation Z (top variable motifs)")
dev.off()

saveRDS(dev, snakemake@output[["rds"]])
cat("\n"); print(sessionInfo())
