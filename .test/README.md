# Test dataset

A synthetic paired-end ATAC-seq dataset small enough to run the pipeline in
minutes. Nothing here is committed except the scripts and config; the data
are rebuilt with

```bash
pixi run build-test                         # downloads chr22 + chrM from UCSC (~12 MB)
bash .test/scripts/build_test_data.sh --fasta /path/to/hg38.fa   # or use a local FASTA
```

`scripts/simulate_atac.py` (standard-library Python) takes a 4 Mb window of
chr22 plus chrM and simulates 50 bp reads: 600 peaks (20% induced over time,
about 5% repressed) with a nucleosome-free/nucleosomal fragment mixture,
background, 8% chrM, two blacklist artefact pileups, 5% PCR duplicates and
Nextera read-through on short fragments. A synthetic GTF puts TSSs on 40% of
the peaks.

Design (`config/samples.tsv`):

| condition | replicates | exercises |
|---|---|---|
| ctrl_0m | 2 | shared time-course baseline |
| stim_30m | 2 | standard two-replicate IDR |
| stim_60m | 3 (replicate 3 in two runs) | all-pairs IDR selection, FASTQ merging |
| stim_120m | 1 | single-replicate IDR fallback |

Replicate 2 libraries are in batch `b2`, the rest in `b1`, and the config sets
`batch: true`, so diff, normcheck and the time course run with the batch
covariate.

This directory is the Snakemake working directory:

```bash
pixi run dry        # snakemake -n -s workflow/Snakefile --directory .test
pixi run test       # alignment -> peaks -> IDR with --sdm conda
pixi run test-all   # everything (also builds the R and TOBIAS environments)
```

Coordinates are relative to the extracted window, not to hg38.

## Head-to-head with nf-core/atacseq

`nfcore/` runs nf-core/atacseq 2.1.2 on the same simulated reads and scores
both pipelines against the simulated truth. It needs Nextflow and Docker
(`NFCORE_PROFILE=singularity` to switch) and the `test-all` outputs:

```bash
pixi run test-all
pixi run test-nfcore      # nf-core run (~30 min at 8 cores, then -resume) + compare
pixi run compare-nfcore   # comparison only, on existing outputs
```

`nfcore/compare_nfcore.py` writes `results/nfcore_compare/{metrics.tsv,report.md}`.
Only these checks can fail it (tolerances in the script docstring):

- precision and F1 of the IDR consensus (`consensus_idr.bed`) vs
  `truth_peaks.bed`, compared with nf-core's merged-library consensus;
- per-library duplicate rate and chrM fraction vs the simulated 5% and 8%,
  compared with nf-core's Picard and idxstats numbers.

Recall, per-condition peak sets, FRiP and the peak-set Jaccard are reported
only. The IDR consensus is stricter than nf-core's merge of every peak call,
so it trades some recall for precision on purpose.
