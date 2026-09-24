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

This directory is the Snakemake working directory:

```bash
pixi run dry        # snakemake -n -s workflow/Snakefile --directory .test
pixi run test       # alignment -> peaks -> IDR with --sdm conda
pixi run test-all   # everything (also builds the R and TOBIAS environments)
```

Coordinates are relative to the extracted window, not to hg38.
