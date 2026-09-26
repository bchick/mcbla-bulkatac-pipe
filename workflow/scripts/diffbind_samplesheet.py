"""Write the DiffBind sample sheet (one row per library).

Factor is the batch when config `batch` is true (design ~Factor + Condition),
otherwise the constant "ATAC".
"""

import csv
import sys

sm = snakemake  # noqa: F821
sys.stderr = open(sm.log[0], "w")

with open(sm.output[0], "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["SampleID", "Tissue", "Factor", "Condition", "Treatment",
                "Replicate", "bamReads", "Peaks", "PeakCaller"])
    for lib, cond, trt, rep, factor, bam, peaks in zip(
            sm.params.libs, sm.params.conditions, sm.params.treatments,
            sm.params.replicates, sm.params.factors, sm.input.bams,
            sm.input.peaks):
        w.writerow([lib, sm.params.tissue, factor, cond, trt, rep, bam, peaks,
                    "narrow"])
