# Changelog

All notable changes to this workflow are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

Results produced with 0.1.0 change after this release: differential
accessibility p-values and FDRs (all contrasts), the normalization check
verdicts, and, slightly, every final BAM and the peaks and counts built from
it. Reprocess and re-run the analyses rather than mixing outputs from the two
versions.

### Fixed

- **Differential accessibility now uses the normalization it reports.**
  Contrasts were given to DiffBind as group masks, which puts DiffBind 3.x in
  its legacy per-contrast mode: the DESeq2 test set its own size factors and
  ignored `dba.normalize()`, while the reported `Conc` and `Fold` used the
  stored factors. As a result, the depth p-values were not computed under
  library-size normalization, and the csaw arm of the normalization check gave
  p-values identical to depth, so it could never flag a contrast. Contrasts now
  use a `~Condition` design, and every diff and csaw run stops if the size
  factors DESeq2 used differ from the stored ones. The AS28 csaw benchmark in
  `docs/defaults_rationale.md` is marked invalid until it is re-run.
- **Orphaned mates are removed from the final BAMs.** The per-read MAPQ filter
  and the blacklist intersect removed reads one at a time, leaving the partner
  of a low-MAPQ or blacklisted mate behind. A new `fix_pairs` step
  (`fixmate -r -m`, then proper pairs only) drops them before deduplication, as
  the ENCODE ATAC pipeline does. The count is reported as `Orphans_Removed` in
  `alignment_qc_report.tsv`.
- DiffBind's own blacklist step is no longer applied inside `dba.analyze`:
  reads and peaks are already filtered against `reference.blacklist`, or
  deliberately not when it is empty.

### Added

- **Per-library QC summary with ENCODE thresholds**
  (`results/qc/qc_summary.tsv`, also in MultiQC). Fragments, alignment rate,
  NRF, PBC1, PBC2, FRiP and TSS enrichment are each flagged PASS / WARN / FAIL
  against `qc.thresholds` (ENCODE ATAC standards by default). Flags are
  reported, never enforced. It comes with:
  - a TSS enrichment score from Tn5 cut sites, using the ENCODE definition;
  - library complexity (NRF, PBC1, PBC2) before deduplication, also added to
    `alignment_qc_report.tsv` (FASTQ mode);
  - FRiP per library over the condition's merged stringent peaks;
  - FastQC on the raw reads (FASTQ mode).
- **Opt-in batch covariate** (`batch: true` plus a `batch` column in
  `samples.tsv`). diff and normcheck fit `~batch + condition`; the time course
  fits an LRT of `~batch + time` against `~batch` and removes batch from the
  VST before clustering. The run stops before any job if a library has no
  batch, or if batch is single-level or confounded with condition or time.
- **Summit-centred fixed-width peak set** (`peaks: fixed_width`,
  `results/peaks/consensus/fixed_width.bed`): 501 bp windows around each
  per-replicate summit, built with the iterative-overlap method of Corces
  et al. 2018 (score per million ≥ 5, supported by ≥ 2 libraries). It keeps
  neighbouring sites separate instead of chaining them into multi-kb regions,
  so chromVAR's motif window sits on a summit. It is now the README's
  suggested set for chromVAR.
- `.test/scripts/check_results.py`, run at the end of `pixi run test-all` (or
  on its own with `pixi run check-test`). It fails if the csaw p-values match
  depth while the folds differ, or if any final BAM holds records that are not
  properly paired.

### Changed

- `alignment_qc_report.tsv` has new columns `Orphans_Removed`, `NRF`, `PBC1`
  and `PBC2`. Scripts that read it by column position need updating.
- The DiffBind sample sheet's `Factor` column holds the batch when
  `batch: true` (otherwise the constant `ATAC`, as before).
- The test dataset now has two batches with `batch: true`, and runs chromVAR
  on `fixed_width`.

### Known issues

- The depth and csaw columns of the AS28 normalization benchmark in
  `docs/defaults_rationale.md` were produced by the bug fixed above and have
  not been re-run yet.
- Snakemake can report "Software environment definition has changed" and plan
  a full rerun even when no file in `workflow/envs/` changed. Until this is
  understood, `--rerun-triggers mtime params input code` avoids it.
- The `multiqc` rule reruns on every invocation, because its config param
  resolves to a temporary cache path.
- The nf-core import is still checked only against the documented 2.x
  layout, not a real nf-core run.

## [0.1.0] - 2026-09-24

First version: a Snakemake port of the McBla lab's AS28 bulk ATAC-seq
pipeline, generalised to any design.

- FASTQ mode (cutadapt, Bowtie 2, MAPQ/flag/chrM/blacklist filtering,
  duplicate removal) or import of filtered BAMs from an nf-core/atacseq 2.x
  outdir.
- MACS2 in four phases, IDR per condition for any number of replicates, and
  two consensus peak sets (`idr_consensus`, `stringent_union`).
- deepTools QC and MultiQC.
- Opt-in analyses, each with an explicit peak set: DiffBind contrasts, a
  normalization check (csaw background bins, quantile + limma), time-course
  clustering (DESeq2 LRT + DEGreport), chromVAR and TOBIAS footprinting.
- Pinned per-stage conda environments, a pixi launcher, Slurm and local
  profiles, a synthetic test dataset, and CI lint and dry runs.

[Unreleased]: https://github.com/bchick/mcbla-bulkatac-pipe/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/bchick/mcbla-bulkatac-pipe/releases/tag/v0.1.0
