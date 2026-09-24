"""Assemble the per-library alignment QC table (port of the QC TSV row written
by 1.1_atac_align_cc.sh).

Columns: Sample, Raw_Reads, Trimmed_Reads, Aligned_Reads, Aligned_Pct,
ChrM_Reads, ChrM_Pct, Blacklist_Removed, Final_Reads, Dup_Pct, Mean_FragSize

Differences from the shell original (deliberate):
  * Raw_Reads comes from cutadapt's "Total read pairs processed" instead of a
    separate zcat pass over R1 (same number, one fewer full read of the FASTQ).
  * Dup_Pct is DUPLICATE TOTAL / EXAMINED from `samtools markdup -f`; the
    original awk matched several "DUPLICATE ..." lines and printed garbage.
  * Final_Reads is the flagstat total (= `samtools view -c` on the final BAM).
"""

import re
import sys

sm = snakemake  # noqa: F821  (injected by Snakemake)
sys.stderr = open(sm.log[0], "w")


def grab(path, pattern, cast=int, default="NA"):
    with open(path) as fh:
        text = fh.read()
    m = re.search(pattern, text, flags=re.MULTILINE)
    if not m:
        return default
    return cast(m.group(1).replace(",", ""))


def kv(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                out[parts[0]] = parts[1]
    return out


def pct(num, den):
    try:
        return f"{100.0 * float(num) / float(den):.1f}%"
    except (ValueError, ZeroDivisionError, TypeError):
        return "NA"


header = [
    "Sample", "Raw_Reads", "Trimmed_Reads", "Aligned_Reads", "Aligned_Pct",
    "ChrM_Reads", "ChrM_Pct", "Blacklist_Removed", "Final_Reads", "Dup_Pct",
    "Mean_FragSize",
]

rows = []
for i, lib in enumerate(sm.params.libs):
    cut = sm.input.cutadapt[i]
    bt2 = sm.input.bowtie2[i]
    raw = grab(cut, r"Total read pairs processed:\s+([\d,]+)")
    trimmed = grab(cut, r"Pairs written \(passing filters\):\s+([\d,]+)")
    aligned_pct = grab(bt2, r"([\d.]+)% overall alignment rate", cast=str)
    filt = kv(sm.input.filt[i])
    total = filt.get("total_aligned", "NA")
    mito = filt.get("mito_reads", "NA")
    bl = filt.get("blacklist_removed", "NA")
    final = grab(sm.input.flagstat[i], r"^(\d+) \+ \d+ in total")
    examined = grab(sm.input.markdup[i], r"^EXAMINED:\s+(\d+)")
    dups = grab(sm.input.markdup[i], r"^DUPLICATE TOTAL:\s+(\d+)")
    s = n = 0
    with open(sm.input.frag[i]) as fh:
        for line in fh:
            size, count = line.split()
            s += int(size) * int(count)
            n += int(count)
    mean_frag = f"{s / n:.0f}" if n else "NA"
    rows.append([
        lib, raw, trimmed, total,
        f"{aligned_pct}%" if aligned_pct != "NA" else "NA",
        mito, pct(mito, total), bl, final, pct(dups, examined), mean_frag,
    ])

with open(sm.output[0], "w") as out:
    out.write("\t".join(header) + "\n")
    for r in rows:
        out.write("\t".join(str(x) for x in r) + "\n")
