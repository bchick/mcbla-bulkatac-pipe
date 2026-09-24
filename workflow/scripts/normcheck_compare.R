# Compare DiffBind default (depth) results with each alternative
# normalization and flag normalization-sensitive contrasts
# (08b NB03 comparison + the >20% rule from 22_normalization_decision).
#
# A contrast is normalization-sensitive under a method when its gained OR
# lost count changes by more than `threshold` (fraction, default 0.20) of
# the depth count AND by at least `min_abs` peaks (so 2 -> 3 peaks is not
# flagged). When the depth count is 0, only the min_abs rule applies.

source(snakemake@params[["helpers"]])
start_log()
suppressPackageStartupMessages(library(ggplot2))

methods <- snakemake@params[["methods"]]
labels <- snakemake@params[["labels"]]
fdr <- as.numeric(snakemake@params[["fdr"]])
thr <- as.numeric(snakemake@params[["threshold"]])
min_abs <- as.numeric(snakemake@params[["min_abs"]])

read_tab <- function(p) read.delim(p, stringsAsFactors = FALSE)
depth <- read_tab(snakemake@input[["depth_summary"]])
depth_tabs <- setNames(snakemake@input[["depth_tables"]], labels)

sig_ids <- function(df) df$peak_id[!is.na(df$FDR) & df$FDR < fdr]
moved <- function(new, old) {
  d <- abs(new - old)
  d >= min_abs & (old == 0 | d > thr * old)
}

rows <- list()
for (k in seq_along(methods)) {
  m <- methods[k]
  summ <- read_tab(snakemake@input[["summaries"]][k])
  for (lab in labels) {
    d <- depth[depth$Contrast == lab, ]
    s <- summ[summ$Contrast == lab, ]
    dt <- read_tab(depth_tabs[[lab]])
    mt <- read_tab(file.path(dirname(snakemake@input[["summaries"]][k]), "tables",
                             paste0(lab, "_all.tsv")))
    a <- sig_ids(dt); b <- sig_ids(mt)
    jac <- if (length(union(a, b)) == 0) NA else length(intersect(a, b)) / length(union(a, b))
    shared <- intersect(dt$peak_id, mt$peak_id)
    r <- if (length(shared) > 2) {
      cor(dt$Fold[match(shared, dt$peak_id)], mt$Fold[match(shared, mt$peak_id)])
    } else NA
    rows[[length(rows) + 1]] <- data.frame(
      Contrast = lab, method = m,
      Gained_depth = d$Gained, Lost_depth = d$Lost,
      Gained = s$Gained, Lost = s$Lost,
      dGained = s$Gained - d$Gained, dLost = s$Lost - d$Lost,
      pct_Gained = ifelse(d$Gained > 0, round(100 * (s$Gained - d$Gained) / d$Gained, 1), NA),
      pct_Lost = ifelse(d$Lost > 0, round(100 * (s$Lost - d$Lost) / d$Lost, 1), NA),
      jaccard_sig = round(jac, 3), lfc_pearson = round(r, 4),
      sensitive = moved(s$Gained, d$Gained) | moved(s$Lost, d$Lost),
      stringsAsFactors = FALSE)
  }
}
tab <- do.call(rbind, rows)
write_tsv(tab, snakemake@output[["table"]])
print(tab)

verdict <- do.call(rbind, lapply(labels, function(lab) {
  x <- tab[tab$Contrast == lab, ]
  sens <- x$method[x$sensitive]
  rec <- if ("csaw" %in% sens) {
    "normalization-sensitive under background bins: global shift suspected; report csaw alongside depth"
  } else if (length(sens) > 0) {
    "sensitive only to quantile normalization (strong equal-distribution assumption); depth default retained"
  } else {
    "robust: depth and alternatives agree"
  }
  data.frame(Contrast = lab,
             normalization_sensitive = length(sens) > 0,
             sensitive_methods = if (length(sens)) paste(sens, collapse = ",") else "",
             recommendation = rec, stringsAsFactors = FALSE)
}))
write_tsv(verdict, snakemake@output[["verdict"]])
print(verdict)

long <- rbind(
  data.frame(Contrast = depth$Contrast, method = "depth", Direction = "Gained", N = depth$Gained),
  data.frame(Contrast = depth$Contrast, method = "depth", Direction = "Lost", N = -depth$Lost),
  data.frame(Contrast = tab$Contrast, method = tab$method, Direction = "Gained", N = tab$Gained),
  data.frame(Contrast = tab$Contrast, method = tab$method, Direction = "Lost", N = -tab$Lost)
)
long$Contrast <- factor(long$Contrast, levels = rev(labels))
flag <- verdict$Contrast[verdict$normalization_sensitive]
p <- ggplot(long, aes(Contrast, N, fill = method, alpha = Direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_hline(yintercept = 0) + coord_flip() +
  scale_alpha_manual(values = c(Gained = 1, Lost = 0.5)) +
  labs(x = NULL, y = sprintf("Differential peaks (FDR < %g); lost plotted negative", fdr),
       title = "Normalization sensitivity",
       subtitle = if (length(flag)) paste("Sensitive:", paste(flag, collapse = ", ")) else
         "No contrast is normalization-sensitive") +
  theme_bw()
ggsave(snakemake@output[["barplot"]], p, width = 10, height = 2 + 0.45 * length(labels))
