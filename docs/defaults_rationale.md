# Defaults and the evidence behind them

Every default in `config/config.yaml` comes from the lab's AS28 ATAC-seq
analysis: MCF7 cells, unstimulated or stimulated with EGF (transient ERK) or
HRG (sustained ERK) for 15/30/60/120/240 min, two biological replicates per
condition, 22 libraries. Where a default was chosen by benchmarking, the
benchmark and its result are given here. Where it is an inherited convention,
it says so.

The analyses (diff, normcheck, timecourse, chromvar, footprint) are off by
default and have no default peak set. The peak sets listed below are the
ones AS28 used, and `config/salk_example.yaml` sets them.

## Trimming and alignment

| Default | Value | Rationale |
|---|---|---|
| Adapter | Nextera `CTGTCTCTTATACACATCT`, both reads | Tn5 library adapters. |
| cutadapt quality / length | `-q 20,20 -m 30 --pair-filter=any` | Trims low-quality 3' and 5' ends; drops the pair if either mate falls below 30 bp. Short NFR fragments still map uniquely at 30 bp. |
| bowtie2 | `--maxins 2000 --no-mixed --no-discordant` | Keeps mono-, di- and tri-nucleosome fragments (the default `--maxins 500` truncates the nucleosomal ladder) and reports only concordant pairs, which is what `-f BAMPE` peak calling needs. |
| MAPQ | `-q 30` | Raised from 10 in an earlier lab pipeline version; removes multi-mapping reads that inflate signal in repeats. |
| Flags | `-f 2 -F 2828` | Proper pairs only; drops unmapped, mate-unmapped, secondary, QC-fail and supplementary records. |
| Mitochondrial reads | removed (`mito_chrom: chrM`) | Often 30 to 50% of ATAC reads. They carry no nuclear accessibility information and distort depth normalization. |
| Blacklist | ENCODE v2 (Amemiya et al. 2019); empty disables | Removes artefact high-signal regions before deduplication and again from every peak file. Set it to empty for genomes without a blacklist (for example rn6, as in the lab's PC12 AS08 run). |
| Duplicates | `samtools markdup -r` (removed) | ATAC duplicates are PCR duplicates. This differs from the lab's CUT&RUN convention of marking without removing, which exists because CUT&RUN cleavage produces real identical fragments. |

## QC

| Default | Value | Rationale |
|---|---|---|
| bigWig | deepTools RPGC, bin 10 bp | Base-level resolution for browsing and for comparing samples across depths. |
| Effective genome size | 2,913,022,398 (GRCh38) | Used for every human dataset in the lab. For rn6 use 2,375,372,135. |
| TSS profile | ±2 kb | Standard ATAC enrichment check. |
| Correlation / PCA | multiBamSummary 500 bp bins, Spearman | Genome-wide and independent of peak calls, so batch effects show up before any peak-based analysis. |

The QC table is what found the AS28 15 min technical batch. The four 15 min
libraries aligned at 80.1 to 82.3% (the other 18: 89.1 to 92.3%), had a mean
fragment of 129 to 138 bp (148 to 167 bp) and a FRiP of 0.285 to 0.303 (0.099
to 0.177). They were PCA outliers and were excluded from the differential and
time-course analyses. The `exclude_conditions` key exists for this situation:
the libraries still get peaks, QC and footprints, but they are kept out of the
statistics. `config/salk_example.yaml` sets it for AS28.

## Peak calling: MACS2, four phases

| Default | Value | Rationale |
|---|---|---|
| Caller | MACS2 2.2.9.1 `-f BAMPE --keep-dup all` | See the benchmark below. BAMPE uses the real fragment span. `--keep-dup all` because duplicates are already removed. |
| Per-replicate peaks | `-p 0.01` (relaxed) | IDR needs a relaxed, noisy peak list to model the reproducible and irreproducible components. |
| Merged relaxed peaks | `-p 0.01` on the merged BAM | The IDR `--peak-list` oracle (ENCODE practice). |
| Merged stringent peaks | `-q 0.05` on the merged BAM | Production per-condition peaks (FRiP, union for footprinting). |

### Peak-caller benchmark (AS28)

AS28 has high background, so we tested whether ATAC-specific callers or
nucleosome-free filtering recover peaks that MACS2 misses. The rule for
switching was a caller calling at least 30% more peaks with a better FRiP.

Sanity condition HRG 60 min (both replicates):

| Method | Peaks | FRiP (NF) |
|---|---|---|
| MACS2, replicate 1 | 90,327 | 0.144 |
| MACS2, replicate union | 61,928 | 0.159 |
| MACS2 on nucleosome-free (<120 bp) pool | 99,142 | 0.119 |
| MACS2 NF pool `--nolambda` | 371,993 | 0.253 |
| Genrich, NF pool | 102,366 | 0.200 |
| HMMRATAC | 31,448 | 0.166 |

What we concluded:

* **Nucleosome-free filtering does not help.** FRiP dropped from 0.159 to
  0.119. The background problem is not a fragment-size problem.
* **`--nolambda` over-calls** (4 to 6× more peaks). The higher FRiP reflects
  more peaks, not recovered signal.
* **Genrich** looked better on this one condition, but only because pooled
  Genrich was being compared with single-replicate MACS2. Compared per
  replicate across the whole dataset, it called 30 to 70% fewer peaks (median
  about 50% fewer) with no consistent FRiP gain. Rejected.
* **HMMRATAC** was very conservative (31k vs 90k peaks), gave only a modest
  FRiP and was the slowest. Rejected.
* **ROCCO** gives one clean consensus across all samples (AS28: 325k peaks,
  FRiP 0.249), which is what it is designed for. Per condition it was erratic
  (collapsing to 11 to 15k peaks at the highest-signal time points). It is
  suitable only for an all-sample union and is not in this pipeline yet
  (listed under "Known gaps" in the README).

MACS2 at these settings therefore stays the default.

## IDR

| Default | Value | Rationale |
|---|---|---|
| Threshold | 0.05, keep column 5 ≥ 540 | ENCODE standard. The scaled IDR is `min(int(-125·log2(IDR)), 1000)`, and −125·log2(0.05) = 540.2. |
| Ranking | `--rank p.value`, inputs `sort -k8,8rn` | MACS2 p-value is the ENCODE ranking measure for relaxed calls. |
| 2 replicates | rep1 vs rep2 with the merged relaxed oracle | Same as the lab script. |
| >2 replicates | `encode_max_pair`: IDR on every pair, keep the pair with the most reproducible peaks | IDR is pairwise. ENCODE's ATAC pipeline takes the maximum over true-replicate pairs as its "conservative" set. We do not generate pseudo-replicates, so there is no "optimal" set or rescue ratio. `first_two` reproduces the old rep1-vs-rep2 behaviour and logs a warning. |
| 1 replicate | IDR skipped; merged stringent peaks enter the consensus (`single_rep_fallback: stringent`) | A single replicate cannot be assessed for reproducibility. The fallback keeps the condition represented and is labelled as a fallback in `idr_summary.tsv`. Use `skip` to leave it out. |
| Failed IDR pair | recorded as `FAILED`, run continues (`allow_failure: true`) | Matches the lab script, where one bad condition does not stop the batch. Set it to false to fail hard. |
| Consensus | `cat` IDR peaks → `sort` → `bedtools merge` | As in the lab script. Unlike nf-core's merge of all peak calls, only reproducible peaks enter. |

## Differential accessibility (DiffBind)

| Default | Value | Rationale |
|---|---|---|
| Peak set | per-replicate relaxed peaks, `minOverlap = 2` | The AS28 analysis. DiffBind builds its own consensus of peaks present in at least 2 libraries and re-centres them on summits (DiffBind default). Set with `diff.peaks: individual`; `idr_consensus`, `stringent_union` or a BED path count over a fixed set instead. |
| Normalization | DiffBind default (library size, DESeq2) | Supported by the normalization benchmark below. |
| Significance | FDR < 0.05 | Tables also carry `Gained_lfc` / `Lost_lfc` counts at \|log2FC\| ≥ 1. |
| Contrasts | explicit group1 vs group2 by condition | Contrasts are defined within one experiment, so they are not confounded by batch. |

## Normalization check

This check is the main selling point. For DiffBind-style analyses the
normalization method matters more than the test, and its assumption (that
most peaks do not change) is exactly what a strong stimulus can violate.

**Evidence that the choice can matter.** In the lab's CUT&RUN data (the SD51
AP-1 degrader experiment), depth normalization understated HRG-induced AP-1
binding by about 2× compared with greenlist (invariant-region) size factors,
while the degrader contrast itself was robust. The lab adopted the rule that a
contrast whose gained-peak count moves by more than 20% under an alternative
normalization is "load-bearing" on normalization and must be reported with
both. That rule is `normcheck.sensitivity_threshold: 0.20`.

**Evidence for the ATAC default.** For AS28 ATAC, the same 10 contrasts were
re-run under csaw background bins (the ATAC equivalent of greenlist: 15 kb
bins, DESeq2 native normalization of the bin counts) and under quantile
normalization + limma:

| Contrast | Depth gained / lost | csaw gained / lost | Quantile gained / lost |
|---|---|---|---|
| EGF 30m vs US | 17 / 0 | 17 / 0 | 0 / 1 |
| HRG 30m vs US | 42 / 80 | 42 / 80 | 14 / 4 |
| EGF 60m vs US | 34 / 20 | 34 / 20 | 0 / 0 |
| HRG 60m vs US | 1,643 / 386 | 1,643 / 386 | 3,306 / 1,205 |
| EGF 120m vs US | 26 / 4 | 26 / 4 | 111 / 2 |
| HRG 120m vs US | 7,137 / 3,185 | 7,137 / 3,185 | 10,222 / 4,601 |
| EGF 240m vs US | 0 / 0 | 0 / 0 | 0 / 0 |
| HRG 240m vs US | 7,552 / 701 | 7,552 / 701 | 7,819 / 1,868 |
| HRG 60m vs EGF 60m | 233 / 5 | 233 / 5 | 0 / 0 |
| HRG 240m vs EGF 240m | 3,681 / 189 | 3,681 / 189 | 4,281 / 975 |

csaw size factors were within 2 to 3% of the depth factors, per-peak log2FC
shifted by less than 0.03, and the significant sets were identical across all
10 contrasts. So library-size normalization is appropriate for AS28, and the
global-shift failure seen in CUT&RUN did not occur in bulk ATAC here.
Quantile normalization moved many contrasts in both directions. Its
assumption that every sample has the same distribution is stronger than
depth's, so a quantile-only change is reported but is not treated as a
reason to switch.

Because the answer depends on the biology, the pipeline does not assume it.
Every run repeats the check and writes `norm_verdict.tsv`:

* sensitive under csaw → "global shift suspected; report csaw alongside depth"
* sensitive only under quantile → "depth default retained"
* otherwise → "robust"

`min_abs_change: 10` stops tiny contrasts (for example 2 → 3 peaks) from being
flagged by the percentage rule alone. This value is new in the pipeline and
has not been benchmarked.

## Time course

| Default | Value | Rationale |
|---|---|---|
| Counted peaks | `timecourse.peaks: idr_consensus` in AS28, widths 100 to 5,000 bp, featureCounts fragments | As in the AS28 preprocessing (width filter on the consensus union). |
| Low-count filter | rowMeans ≥ 10 | AS28 preprocessing. |
| VST | `vst(blind = TRUE)` on all included samples | AS28 preprocessing. |
| Dynamic peaks | DESeq2 LRT `~time` vs `~1` per series, padj < 0.01 **and** range of per-time VST means ≥ 0.5 | The LRT alone selects many peaks with tiny effects at 22-library depth. The effect-size floor keeps clustering on peaks that actually move. |
| Clustering | DEGreport `degPatterns`, `minc = 50`, `cutoff = 0.5`, `set.seed(42)` | Produced the interpretable AS28 HRG solution of five superclusters (decreasing, transient, transient-increasing, sustained-increasing, late-increasing). |
| Shared baseline | `baseline_treatments` (e.g. `unstim`) included in every series at its time | The unstimulated sample is time 0 for both ligands. |

## chromVAR

| Default | Value | Rationale |
|---|---|---|
| Motifs | JASPAR2020 CORE vertebrates, latest versions (746 matrices) | Same set as the TOBIAS run, so TF activity and footprints are comparable. |
| Counted peaks | `chromvar.peaks: idr_consensus` in AS28 | Same counts as the time course. |
| Motif windows | 200 bp centred on each peak | Motif matches in wide merged peaks dilute the per-peak signal. |
| Peak filter | ≥ 10 fragments summed over samples | chromVAR recommendation for low-coverage peaks. |
| GC bias | `addGCBias` with the genome FASTA (or a BSgenome) | FASTA via `Rsamtools::FaFile` avoids requiring a BSgenome package. Set `chromvar.bsgenome` to use one. |
| Background | `getBackgroundPeaks(niterations = 200, w = 0.1)`, `set.seed(2025)` | AS28 settings. |

## Footprinting (TOBIAS)

| Default | Value | Rationale |
|---|---|---|
| Signal | per-condition merged BAMs | Footprints need depth. Replicates are pooled per condition, as in the AS28 run. |
| Regions | `footprint.peaks: stringent_union` in AS28 | One region set across all conditions keeps BINDetect comparisons like for like. |
| Steps | ATACorrect (Tn5 bias) → ScoreBigwig → one BINDetect across all conditions | TOBIAS 0.16.1, as used for AS28 (746 motifs × 11 conditions). |
| Motifs | `footprint.motifs`; if empty, JASPAR2020 CORE vertebrates exported with TFBSTools | Same set as chromVAR. |

## Tool versions

Pinned per stage in `workflow/envs/` to the versions of the lab toolchain:
bowtie2 2.5.2, samtools 1.13, bedtools 2.30.0, cutadapt 4.6, macs2 2.2.9.1,
deeptools 3.5.4, multiqc 1.17, idr 2.0.4.2, subread 2.0.3, TOBIAS 0.16.1, and
Bioconductor 3.17 with DiffBind 3.10, DESeq2 1.40 and DEGreport 1.36. They
are separate environments because the old Python tools do not co-solve.
