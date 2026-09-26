"""TSS enrichment score for one library (ENCODE ATAC-seq definition).

Tn5 cut sites (read 5' ends shifted +4 on the forward strand, -5 on the
reverse strand) are piled up, strand-aware, in a +/- flank window around every
unique TSS. The aggregate profile is divided by the mean of its outer `edge` bp
at each end (background), smoothed with a `smooth` bp moving average, and the
score is the maximum of the smoothed profile.

Outputs:
  tsv     : Sample, TSS_Enrichment, TSS_Used, Cut_Sites, Background
  profile : position (bp from TSS), cut sites, enrichment (smoothed)
"""

import sys

import numpy as np
import pysam

sm = snakemake  # noqa: F821  (injected by Snakemake)
sys.stderr = open(sm.log[0], "w")

flank = int(sm.params.flank)
edge = int(sm.params.edge)
smooth = int(sm.params.smooth)
mito = sm.params.mito

bam = pysam.AlignmentFile(sm.input.bam)
refs = set(bam.references)
profile = np.zeros(2 * flank + 1)
n_tss = 0

with open(sm.input.tss) as fh:
    for line in fh:
        chrom, start, _end, _name, _score, strand = line.rstrip("\n").split("\t")[:6]
        if chrom not in refs or chrom == mito:
            continue
        tss = int(start)
        n_tss += 1
        for read in bam.fetch(chrom, max(0, tss - flank - 5), tss + flank + 5):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_reverse:
                cut = read.reference_end - 1 - 5
            else:
                cut = read.reference_start + 4
            rel = cut - tss if strand != "-" else tss - cut
            if -flank <= rel <= flank:
                profile[rel + flank] += 1

background = np.concatenate([profile[:edge], profile[-edge:]]).mean()
if background > 0:
    enrichment = profile / background
    kernel = np.ones(smooth) / smooth
    smoothed = np.convolve(enrichment, kernel, mode="same")
    score = f"{smoothed.max():.2f}"
else:
    smoothed = np.full_like(profile, np.nan)
    score = "NA"
    print(f"WARNING: no background cut sites for {sm.wildcards.lib}", file=sys.stderr)

with open(sm.output.tsv, "w") as out:
    out.write("Sample\tTSS_Enrichment\tTSS_Used\tCut_Sites\tBackground\n")
    out.write(f"{sm.wildcards.lib}\t{score}\t{n_tss}\t{int(profile.sum())}\t{background:.3f}\n")

with open(sm.output.profile, "w") as out:
    out.write("position\tcut_sites\tenrichment\n")
    for i in range(len(profile)):
        out.write(f"{i - flank}\t{int(profile[i])}\t{smoothed[i]:.4f}\n")

print(f"{sm.wildcards.lib}: TSS enrichment {score} over {n_tss} TSSs", file=sys.stderr)
