# Shared helpers for the R steps (sourced via snakemake@params$helpers).

start_log <- function() {
  log <- file(snakemake@log[[1]], open = "wt")
  sink(log)
  sink(log, type = "message")
  invisible(log)
}

read_contrasts <- function(path) {
  ct <- read.delim(path, stringsAsFactors = FALSE, comment.char = "#")
  stopifnot(all(c("group1", "group2", "label") %in% colnames(ct)))
  ct
}

# Add explicit two-group contrasts (group1 vs group2 by Condition), in the
# order of contrasts.tsv, exactly as 08 NB01 did with dba$masks.
add_contrasts <- function(dba_obj, ct) {
  dba_obj$contrasts <- NULL
  for (i in seq_len(nrow(ct))) {
    g1 <- DiffBind::dba.mask(dba_obj, DiffBind::DBA_CONDITION, ct$group1[i])
    g2 <- DiffBind::dba.mask(dba_obj, DiffBind::DBA_CONDITION, ct$group2[i])
    if (sum(g1) == 0 || sum(g2) == 0) {
      stop(sprintf("Contrast %s: no samples for %s or %s",
                   ct$label[i], ct$group1[i], ct$group2[i]))
    }
    dba_obj <- DiffBind::dba.contrast(dba_obj, group1 = g1, group2 = g2,
                                      name1 = ct$group1[i], name2 = ct$group2[i])
  }
  dba_obj
}

# Full (th = 1) report for contrast i as a data.frame with stable columns.
report_df <- function(dba_obj, i) {
  res <- DiffBind::dba.report(dba_obj, contrast = i, th = 1, bCounts = FALSE)
  df <- as.data.frame(res)
  conc_cols <- grep("^Conc_", colnames(df))
  if (length(conc_cols) == 2) colnames(df)[conc_cols] <- c("Conc_group1", "Conc_group2")
  df$peak_id <- sprintf("%s:%d-%d", df$seqnames, df$start, df$end)
  df
}

summarise_contrast <- function(df, label, fdr, lfc) {
  sig <- !is.na(df$FDR) & df$FDR < fdr
  data.frame(
    Contrast   = label,
    Total_peaks = nrow(df),
    Gained     = sum(sig & df$Fold > 0),
    Lost       = sum(sig & df$Fold < 0),
    Sig_total  = sum(sig),
    Gained_lfc = sum(sig & df$Fold >= lfc),
    Lost_lfc   = sum(sig & df$Fold <= -lfc),
    stringsAsFactors = FALSE
  )
}

write_tsv <- function(df, path) {
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}
