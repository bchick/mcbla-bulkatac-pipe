# Outputs

All paths are relative to the Snakemake working directory. `<lib>` is
`<sample>_REP<replicate>`, `<cond>` is a `condition` value, `<label>` a
contrast label and `<series>` a time-course treatment. Logs for every step are
under `logs/`, mirroring this layout.

## Reference (`results/reference/`)

| File | Description |
|---|---|
| `genome.fa`, `genome.fa.fai`, `chrom.sizes` | link to (or decompressed copy of) `reference.fasta`, indexed |
| `bowtie2/genome.*.bt2` | only when `reference.bowtie2_index` is empty |
| `jaspar_motifs.jaspar` | only when `footprint.motifs` is empty (JASPAR2020 export) |

## Alignment (`results/bam/`, `results/qc/`)

| File | Description |
|---|---|
| `bam/<lib>.final.bam` (+ `.bai`) | MAPQ ≥ 30, proper pairs, no chrM, blacklist-filtered, duplicates removed. In nfcore mode, a link to `<aligner>/merged_library/<lib>.mLb.clN.sorted.bam` |
| `bam/merged/<cond>.merged.bam` (+ `.bai`) | replicates merged per condition (MACS2 merged calls, TOBIAS) |
| `qc/alignment_qc_report.tsv` | Sample, Raw_Reads, Trimmed_Reads, Aligned_Reads, Aligned_Pct, ChrM_Reads, ChrM_Pct, Blacklist_Removed, Final_Reads, Dup_Pct, Mean_FragSize (fastq mode) |
| `qc/filter_stats/<lib>.filter_stats.tsv` | read counts before/after the mito and blacklist filters |
| `qc/markdup/<lib>.markdup.txt` | `samtools markdup -f` statistics |
| `qc/flagstat/<lib>.flagstat.txt` | `samtools flagstat` of the final BAM |
| `qc/fragment_sizes/<lib>_fragment_sizes.tsv` | fragment size (bp) and count, < 1 kb |

## QC (`results/bigwig/`, `results/qc/deeptools/`, `results/qc/multiqc/`)

| File | Description |
|---|---|
| `bigwig/<lib>.bw` | RPGC-normalized coverage, 10 bp bins |
| `qc/deeptools/tss_enrichment_profile.png` / `.tab` | TSS ±2 kb profile (from `tss_matrix.gz`) |
| `qc/deeptools/fragment_size_distribution.png`, `fragment_size_table.tsv`, `fragment_size_raw.tsv` | bamPEFragmentSize |
| `qc/deeptools/multiBamSummary.npz` | 500 bp bin counts |
| `qc/deeptools/correlation_spearman.png` / `.tsv` | Spearman correlation heatmap and matrix |
| `qc/deeptools/pca_plot.png`, `pca_data.tsv` | PCA of the bin counts |
| `qc/deeptools/fingerprint.png`, `fingerprint_metrics.tsv`, `fingerprint_counts.tsv` | plotFingerprint |
| `qc/multiqc/multiqc_report.html` | MultiQC: cutadapt, Bowtie 2, samtools, deepTools, MACS2, plus pipeline tables (alignment QC, peak summaries, IDR) |

## Peaks (`results/peaks/`)

| File | Description |
|---|---|
| `individual/<lib>_peaks.narrowPeak` (+ `.xls`, `_summits.bed`) | per-replicate MACS2 `-p 0.01`, blacklist-filtered |
| `merged_relaxed/<cond>_peaks.narrowPeak` | merged BAM, `-p 0.01` (IDR oracle) |
| `merged_stringent/<cond>_peaks.narrowPeak` | merged BAM, `-q 0.05` (production peaks) |
| `<set>/peak_summary.tsv` | Peaks, FRiP (reads in peaks / reads, `-F 2304`), Median_Width per file |
| `idr/pairs/<cond>/<a>__<b>.idr_all.narrowPeak` | raw IDR output for one replicate pair (plus `.png` plot) |
| `idr/pairs/<cond>/<a>__<b>.idr.narrowPeak` | pair peaks passing scaled IDR ≥ 540 |
| `idr/<cond>_idr.narrowPeak` | reproducible peaks per condition (the best pair when > 2 reps) |
| `idr/idr_summary.tsv` | Condition, N_Reps, Rep1/Rep2/Oracle/IDR peaks, Reprod_Rate, Selected_Pair, All_Pairs, Status (`OK`, `PARTIAL`, `FAILED`, `single_rep:*`) |
| `consensus/consensus_idr.bed` (+ `.summary.tsv`) | merge of per-condition reproducible peaks |
| `consensus/union_stringent.bed` | merge of per-condition stringent peaks |

## Counts (`results/counts/`, analyses only)

Built only for the analyses that are switched on.

| File | Description |
|---|---|
| `<peak set>/peaks.saf`, `<peak set>/peak_counts.tsv` | peaks (width-filtered) and a fragments × libraries matrix (featureCounts) for each peak set named in `timecourse.peaks` / `chromvar.peaks`; a BED path gets the directory `custom_<file stem>` |
| `diffbind_samplesheet.csv`, `dba_counted.rds` | DiffBind sample sheet and the counted DBA shared by diff and normcheck |

## Differential accessibility (`results/diff/depth/`)

| File | Description |
|---|---|
| `tables/<label>_all.tsv` | every peak: seqnames, start, end, Conc, Conc_group1, Conc_group2, Fold (log2 group1/group2), p.value, FDR, peak_id |
| `tables/<label>_sig.tsv` | FDR < `diff.fdr` |
| `summary.tsv` | Contrast, Total_peaks, Gained, Lost, Sig_total, Gained_lfc, Lost_lfc |
| `plots/ma_plots.pdf`, `pca.pdf`, `correlation_heatmap.pdf`, `summary_barplot.pdf` | DiffBind plots |
| `dba_analyzed.rds` | analysed DBA object |

## Normalization check (`results/normcheck/`)

| File | Description |
|---|---|
| `csaw/size_factors.tsv` | library sizes and csaw background-bin normalization factors |
| `csaw/tables/<label>_all.tsv`, `csaw/summary.tsv`, `csaw/dba_csaw.rds` | contrasts re-run with background-bin normalization |
| `quantile/tables/<label>_all.tsv`, `quantile/summary.tsv`, `quantile/quantile_diagnostics.tsv` | quantile normalization + limma |
| `norm_comparison.tsv` | per contrast × method: gained/lost, deltas and % change vs depth, Jaccard of significant sets, log2FC Pearson r, `sensitive` |
| `norm_verdict.tsv` | per contrast: `normalization_sensitive`, `sensitive_methods`, recommendation |
| `norm_comparison_barplot.pdf` | gained/lost under each method |

## Time course (`results/timecourse/<series>/`)

| File | Description |
|---|---|
| `lrt_results.tsv` | per peak: baseMean, LRT stat, pvalue, padj, vst_range, dynamic |
| `degpatterns_clusters.tsv` | peak_id → degPatterns cluster |
| `degpatterns_profiles.tsv` | mean Z-scored VST per cluster and time |
| `degpatterns_clusters.pdf`, `timecourse.rds` | cluster plot; LRT + degPatterns objects |

## chromVAR (`results/chromvar/`)

| File | Description |
|---|---|
| `deviation_zscores.tsv`, `deviations.tsv` | motif × library deviation Z-scores and bias-corrected deviations |
| `variability.tsv` | per-motif variability (sorted) |
| `top_variable_heatmap.pdf`, `chromvar_deviations.rds` | heatmap of the top variable motifs; chromVARDeviations object |

## Footprinting (`results/footprint/`)

| File | Description |
|---|---|
| `atacorrect/<cond>_{corrected,bias,expected,uncorrected}.bw` | TOBIAS ATACorrect tracks |
| `footprintscores/<cond>_footprints.bw` | TOBIAS ScoreBigwig |
| `bindetect/bindetect_results.txt` (+ `.xlsx`, per-motif folders, figures) | TOBIAS BINDetect across `footprint.conditions` |
