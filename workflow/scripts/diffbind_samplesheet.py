"""Write the DiffBind sample sheet (one row per library)."""

import csv
import sys

sm = snakemake  # noqa: F821
sys.stderr = open(sm.log[0], "w")

with open(sm.output[0], "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["SampleID", "Tissue", "Factor", "Condition", "Treatment",
                "Replicate", "bamReads", "Peaks", "PeakCaller"])
    for lib, cond, trt, rep, bam, peaks in zip(
            sm.params.libs, sm.params.conditions, sm.params.treatments,
            sm.params.replicates, sm.input.bams, sm.input.peaks):
        w.writerow([lib, sm.params.tissue, "ATAC", cond, trt, rep, bam, peaks,
                    "narrow"])
