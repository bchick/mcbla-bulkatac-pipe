#!/usr/bin/env bash
# run_nfcore.sh -- run nf-core/atacseq on the synthetic .test dataset.
#
# Writes an nf-core samplesheet from .test/config/samples.tsv (absolute FASTQ
# paths; rows sharing sample + replicate are runs, merged as in our pipeline)
# and runs the pinned nf-core revision with params.yaml + nfcore.config.
# Everything lands in .test/nfcore/run/ (gitignored); reruns use -resume.
# compare_nfcore.py then scores both pipelines against the simulated truth.
#
# Usage (repo root):  pixi run test-nfcore-run
#   NFCORE_PROFILE=singularity pixi run test-nfcore-run   (default: docker)
#   NFCORE_REVISION=2.1.2                                 (pinned default)
#   NXF_VER=25.10.0                                       (pinned default; downloaded on first use)
set -euo pipefail
REPO=$(cd "$(dirname "$0")/../.." && pwd)
HERE=$REPO/.test/nfcore
RUN=$HERE/run
PROFILE=${NFCORE_PROFILE:-docker}
REVISION=${NFCORE_REVISION:-2.1.2}
DATA=$REPO/.test/data
export NXF_VER=${NXF_VER:-25.10.0}

command -v nextflow >/dev/null || { echo "nextflow not on PATH (see .test/README.md)" >&2; exit 1; }
[[ -s $DATA/ref/genome.fa ]] || { echo "no test data; run: pixi run build-test" >&2; exit 1; }

mkdir -p "$RUN"
awk -F'\t' -v d="$DATA/fastq" 'BEGIN{OFS=","; print "sample","fastq_1","fastq_2","replicate"}
    NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    {print $c["sample"], d "/" $c["fastq_1"], d "/" $c["fastq_2"], $c["replicate"]}' \
    "$REPO/.test/config/samples.tsv" > "$RUN/samplesheet.csv"

cd "$RUN"
nextflow run nf-core/atacseq -r "$REVISION" -profile "$PROFILE" \
    -params-file "$HERE/params.yaml" -c "$HERE/nfcore.config" \
    --input "$RUN/samplesheet.csv" --outdir "$RUN/results" \
    --fasta "$DATA/ref/genome.fa" --gtf "$DATA/ref/genes.gtf" \
    --blacklist "$DATA/ref/blacklist.bed" \
    -resume
