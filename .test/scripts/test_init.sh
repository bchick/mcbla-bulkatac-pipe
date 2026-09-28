#!/usr/bin/env bash
# test_init.sh -- check workflow/scripts/init_project.py against a fixture manifest.
#
# Builds a throwaway manifest whose "resources" are empty placeholder files,
# runs init non-interactively for each genome x blacklist choice, and dry-runs
# the Snakemake DAG of every generated project. Also checks that init refuses
# to guess when answers are missing.
#
# Usage (repo root):  pixi run test-init
set -euo pipefail
REPO=$(cd "$(dirname "$0")/../.." && pwd)
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT

# --- placeholder resources ---------------------------------------------------
R=$T/resource
mkdir -p "$R"
for g in hs mm rn; do
    touch "$R/$g.fa" "$R/$g.fa.fai" "$R/$g.gtf" "$R/$g.blacklist.bed"
    for s in 1 2 3 4 rev.1 rev.2; do touch "$R/$g.$s.bt2"; done
done
cat > "$T/manifest.yaml" <<YAML
manifest_version: 1
updated: test
genomes:
  hg38:
    description: "fixture human"
    aliases: [GRCh38]
    fasta: $R/hs.fa
    bowtie2_index: $R/hs
    gtf: $R/hs.gtf
    effective_genome_size: 2913022398
    macs2_gsize: hs
    mito_chrom: chrM
    blacklists:
      encode_v2: {path: $R/hs.blacklist.bed, default: true}
      old: {path: $R/hs.blacklist.bed, status: unverified}
  mm39:
    description: "fixture mouse"
    fasta: $R/mm.fa
    bowtie2_index: $R/mm
    gtf: $R/mm.gtf
    effective_genome_size: 2654621783
    macs2_gsize: mm
    blacklists:
      excluderanges: {path: $R/mm.blacklist.bed, default: true}
  rn6:
    description: "fixture rat, no blacklist, no GTF"
    fasta: $R/rn.fa
    bowtie2_index: $R/rn
    effective_genome_size: 2651731494
    macs2_gsize: "2651731494"
    blacklists: {}
YAML

init() { python "$REPO/workflow/scripts/init_project.py" --manifest "$T/manifest.yaml" "$@" < /dev/null; }
fail() { echo "FAIL: $*" >&2; exit 1; }

# --- refusals ------------------------------------------------------------------
set +e
init --dir "$T/x" --genome hg38 > /dev/null 2>&1; [[ $? == 2 ]] || fail "missing answers should exit 2"
init --dir "$T/x" --genome hg38 --blacklist old > /dev/null 2>&1
[[ $? == 1 ]] || fail "unverified blacklist should be refused"
init --dir "$T/x" --genome hg19 --blacklist none > /dev/null 2>&1
[[ $? == 1 ]] || fail "unknown genome should be refused"
set -e
init --list > /dev/null
init --list --json | python -c "import json,sys; d=json.load(sys.stdin); assert 'old' not in d['hg38']['blacklists']"

# --- every combination dry-runs -------------------------------------------------
n=0
for combo in hg38:encode_v2 hg38:none mm39:excluderanges mm39:none rn6:none; do
    genome=${combo%%:*}; bl=${combo#*:}
    d=$T/proj_${genome}_${bl}
    init --dir "$d" --genome "$genome" --blacklist "$bl" > /dev/null
    mkdir -p "$d/fastq"
    for f in $(tail -n +2 "$d/samples.tsv" | cut -f2,3); do touch "$d/fastq/$f"; done
    snakemake -n -s "$REPO/workflow/Snakefile" --directory "$d" \
        --configfile "$d/project.yaml" > "$d.log" 2>&1 \
        || { cat "$d.log"; fail "dry run $genome blacklist=$bl"; }
    [[ $genome == rn6 ]] && grep -q "tss_profile" "$d.log" && fail "TSS profile without a GTF"
    n=$((n + 1))
done
echo "test_init: refusals OK, $n generated projects dry-run OK"
