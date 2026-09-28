# AGENTS.md — running mcbla-bulkatac-pipe for a lab member

Instructions for AI agents (Claude Code, Codex, …) asked to run this bulk
ATAC-seq pipeline on the McBla lab server. People can follow them too. For
pipeline internals, read `README.md` and `docs/`.

## 1. Before running anything, ask the user

Ask these questions and wait for the answers. **Do not guess or pick defaults
without asking**, even if the FASTQ names seem to hint at an answer.

1. **Which genome?** Run `pixi run init --list` and show the options
   (`hg38`, `mm10`, `mm39`, `rn6`). mm10 and mm39 are different builds, so the
   user should confirm which one their other data uses.
2. **Blacklist?** The genome's default (ENCODE v2 for hg38/mm10, excluderanges
   for mm39), `none`, or a path to their own BED. rn6 has no blacklist. The
   blacklist is applied to reads, peaks and TOBIAS footprints.

Also ask where the FASTQs are, where the project directory should go (usually
under their own `/data/<user>/`, never inside this repo or `/data/resource`),
what the conditions and replicates are, and which comparisons they want
(contrasts, time course).

This is ATAC-seq, so there are no IgG controls or spike-ins. If the user asks
for spike-in normalization, explain that it doesn't apply here. The
normalization check (`normcheck`) compares depth, csaw background and quantile
normalization instead.

## 2. Write the project config with `init`

Pass the user's answers as flags:

```bash
cd /path/to/mcbla-bulkatac-pipe
pixi run init --dir /data/<user>/<project> --genome hg38 --blacklist encode_v2 \
    --fastq-dir /path/to/fastqs
```

`init` takes the reference paths from the lab manifest,
`/data/resource/manifest.yaml`. It checks that each file exists, then writes
`<project>/project.yaml` and copies example `samples.tsv` / `contrasts.tsv`
(the AS28 design) into the project. It exits with status 2 if an answer is
missing; that means go back and ask the user.

- **Never download or build a genome, index or blacklist yourself**, and never
  point the config at files in someone's personal directory. If something the
  user needs isn't in the manifest, stop and tell them. New resources go into
  `/data/resource` through its maintainer (see
  `/data/resource/_admin/docs/README.md`).
- Entries marked `unverified` in the manifest are hidden. Only use
  `--allow-unverified` if the user explicitly asks for it after you explain
  why the entry is unverified.

## 3. Fill in the samplesheet

Replace the example rows in `<project>/samples.tsv` (columns are described in
the README, "Describe your samples"). Show the finished sheet to the user
before running.

- `condition` groups replicates.
- `treatment` + `time` are only needed for the time-course module.
- Add a `batch` column only if the user wants a batch covariate
  (`batch: true`).
- `contrasts.tsv` lists `group1 group2 label` (log2FC = group1 / group2).

## 4. Dry run, then run

```bash
pixi run snakemake -s workflow/Snakefile --directory /data/<user>/<project> \
    --configfile /data/<user>/<project>/project.yaml -n            # dry run: fix any errors first
pixi run snakemake -s workflow/Snakefile --directory /data/<user>/<project> \
    --configfile /data/<user>/<project>/project.yaml --profile profiles/local --cores 32
```

- This server has **no Slurm**, so use `profiles/local`. It is shared, so keep
  `--cores` at 32 or fewer unless the user says otherwise.
- The run takes hours. Start it in the background or under `tmux`/`nohup`, and
  tell the user where the log is (`<project>/.snakemake/log/`).

## 5. Things that need the user's decision

- **Downstream analyses** (`diff`, `normcheck`, `timecourse`, `chromvar`,
  `footprint`) are off by default. Each needs a peak set (`peaks:`) and turns
  on only when the user asks for it. `footprint` also needs a motif file
  (`footprint.motifs`, JASPAR or MEME format); the lab's motif databases are
  in `/data/resource/meme_motif_files/`.
- **Excluding conditions** (`exclude_conditions`), for example a failed
  batch, is the user's call.
