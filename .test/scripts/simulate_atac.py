#!/usr/bin/env python3
"""Simulate a tiny paired-end ATAC-seq dataset for the pipeline test.

Standard library only. Produces, under --outdir:
  ref/genome.fa        a 4 Mb window of chr22 (renamed "chr22") + chrM
  ref/genes.gtf        synthetic transcripts whose TSSs sit on simulated peaks
  ref/blacklist.bed    two synthetic blacklist regions (with artefact pileups)
  ref/truth_peaks.bed  simulated peak centres (for eyeballing)
  fastq/*.fastq.gz     paired-end 50 bp reads, one pair of files per run

Model: fragments come from peaks (nucleosome-free + mono/di-nucleosome length
mixture), from the genome background, from chrM and from blacklist artefacts.
A subset of peaks is induced over time in the "stim" series so that the
differential/time-course modules have signal. Fragments shorter than the read
length carry Nextera adapter read-through, so cutadapt has work to do.
Coordinates are relative to the extracted window, not to hg38.
"""

import argparse
import gzip
import math
import os
import random

ADAPTER = "CTGTCTCTTATACACATCT"
READ_LEN = 50
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# (sample, replicate label, run label, time)
RUNS = [
    ("ctrl_0m", "r1", "r1", 0),
    ("ctrl_0m", "r2", "r2", 0),
    ("stim_30m", "r1", "r1", 30),
    ("stim_30m", "r2", "r2", 30),
    ("stim_60m", "r1", "r1", 60),
    ("stim_60m", "r2", "r2", 60),
    ("stim_60m", "r3", "r3a", 60),   # replicate 3 sequenced in two runs
    ("stim_60m", "r3", "r3b", 60),
    ("stim_120m", "r1", "r1", 120),
]


def read_fasta_region(path, contig, start=None, end=None):
    """Return sequence of contig[start:end]; uses .fai when present."""
    fai = path + ".fai"
    if os.path.exists(fai) and not path.endswith(".gz"):
        with open(fai) as fh:
            for line in fh:
                name, length, offset, lb, lw = line.split("\t")[:5]
                if name == contig:
                    length, offset, lb, lw = int(length), int(offset), int(lb), int(lw)
                    break
            else:
                raise SystemExit(f"{contig} not in {fai}")
        s = 0 if start is None else start
        e = length if end is None else min(end, length)
        first = offset + (s // lb) * lw + s % lb
        last = offset + (e // lb) * lw + e % lb
        with open(path, "rb") as fh:
            fh.seek(first)
            raw = fh.read(last - first).decode()
        return raw.replace("\n", "").replace("\r", "").upper()
    opener = gzip.open if path.endswith(".gz") else open
    seq, keep = [], False
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if keep:
                    break
                keep = line[1:].split()[0] == contig
                continue
            if keep:
                seq.append(line.strip())
    if not seq:
        raise SystemExit(f"{contig} not found in {path}")
    s = "".join(seq).upper()
    return s[start:end]


def revcomp(s):
    return s.translate(COMP)[::-1]


def mutate(s, rng, rate=0.001):
    if rate <= 0:
        return s
    out = list(s)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice("ACGT")
    return "".join(out)


def frag_len(rng, kind):
    if kind == "peak":
        u = rng.random()
        if u < 0.6:
            mu, sd = 80, 20
        elif u < 0.9:
            mu, sd = 200, 25
        else:
            mu, sd = 380, 30
    else:
        mu, sd = 200, 60
    return int(min(800, max(40, rng.gauss(mu, sd))))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, help="FASTA with chr22 and chrM (plain+.fai, or .gz)")
    ap.add_argument("--chrm-fasta", default=None, help="separate FASTA holding chrM (optional)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--start", type=int, default=20_000_000)
    ap.add_argument("--length", type=int, default=4_000_000)
    ap.add_argument("--pairs", type=int, default=80_000, help="read pairs per run")
    ap.add_argument("--peaks", type=int, default=600)
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args()

    rng = random.Random(a.seed)
    ref = os.path.join(a.outdir, "ref")
    fq = os.path.join(a.outdir, "fastq")
    os.makedirs(ref, exist_ok=True)
    os.makedirs(fq, exist_ok=True)

    chrom = read_fasta_region(a.fasta, "chr22", a.start, a.start + a.length)
    chrm = read_fasta_region(a.chrm_fasta or a.fasta, "chrM")
    L = len(chrom)
    with open(os.path.join(ref, "genome.fa"), "w") as out:
        for name, s in (("chr22", chrom), ("chrM", chrm)):
            out.write(f">{name}\n")
            for i in range(0, len(s), 60):
                out.write(s[i:i + 60] + "\n")

    def ok(s):
        return s.count("N") < 0.05 * len(s)

    # peak centres on non-N sequence, >= 2 kb apart
    centres = []
    while len(centres) < a.peaks:
        c = rng.randrange(5_000, L - 5_000)
        if ok(chrom[c - 500:c + 500]) and all(abs(c - x) > 2_000 for x in centres):
            centres.append(c)
    centres.sort()
    strength = [math.exp(rng.gauss(0, 0.7)) for _ in centres]
    # 20% of peaks induced by stimulation (rising with time), 5% repressed
    kind = ["induced" if rng.random() < 0.2 else "repressed" if rng.random() < 0.06 else "static"
            for _ in centres]

    with open(os.path.join(ref, "truth_peaks.bed"), "w") as out:
        for c, s, k in zip(centres, strength, kind):
            out.write(f"chr22\t{c - 150}\t{c + 150}\t{k}\t{s:.3f}\n")

    # blacklist: two 5 kb regions, each with an artefact pileup
    bl = []
    while len(bl) < 2:
        c = rng.randrange(10_000, L - 10_000)
        clear = all(abs(c - x) > 4_000 for x in centres) and all(abs(c - x) > 50_000 for x in bl)
        if clear and ok(chrom[c - 2_500:c + 2_500]):
            bl.append(c)
    with open(os.path.join(ref, "blacklist.bed"), "w") as out:
        for c in sorted(bl):
            out.write(f"chr22\t{c - 2_500}\t{c + 2_500}\n")

    # GTF: transcripts starting at 40% of peaks, random strand
    with open(os.path.join(ref, "genes.gtf"), "w") as out:
        for i, c in enumerate(centres):
            if rng.random() > 0.4:
                continue
            strand = rng.choice("+-")
            tss = c + 1
            s, e = (tss, min(L, tss + 5_000)) if strand == "+" else (max(1, tss - 5_000), tss)
            gid, tid = f"GENE{i:04d}", f"TX{i:04d}"
            attr = f'gene_id "{gid}"; transcript_id "{tid}"; gene_name "{gid}";'
            for feat in ("gene", "transcript", "exon"):
                a_ = f'gene_id "{gid}"; gene_name "{gid}";' if feat == "gene" else attr
                out.write(f"chr22\tsim\t{feat}\t{s}\t{e}\t.\t{strand}\t.\t{a_}\n")

    rep_noise = {}
    for sample, rep, run, time in RUNS:
        key = (sample, rep)
        if key not in rep_noise:
            rep_noise[key] = [math.exp(rng.gauss(0, 0.2)) for _ in centres]
        noise = rep_noise[key]
        stim = 0.0 if sample.startswith("ctrl") else min(1.0, time / 60)
        w = []
        for s, k, n in zip(strength, kind, noise):
            f = 1.0
            if k == "induced":
                f = 0.15 + 3.0 * stim
            elif k == "repressed":
                f = 1.0 - 0.7 * stim
            w.append(s * f * n)
        cum, tot = [], 0.0
        for x in w:
            tot += x
            cum.append(tot)

        def pick_peak():
            r = rng.random() * tot
            lo, hi = 0, len(cum) - 1
            while lo < hi:
                mid = (lo + hi) // 2
                if cum[mid] < r:
                    lo = mid + 1
                else:
                    hi = mid
            return centres[lo]

        r1p = os.path.join(fq, f"{sample}_{run}_R1.fastq.gz")
        r2p = os.path.join(fq, f"{sample}_{run}_R2.fastq.gz")
        with gzip.open(r1p, "wt", compresslevel=3) as o1, gzip.open(r2p, "wt", compresslevel=3) as o2:
            n = 0
            prev = None
            while n < a.pairs:
                u = rng.random()
                if prev is not None and rng.random() < 0.05:
                    src, fs, fl = prev            # PCR duplicate
                elif u < 0.40:
                    src, fl = chrom, frag_len(rng, "peak")
                    fs = int(pick_peak() + rng.gauss(0, 60) - fl / 2)
                elif u < 0.48:
                    src, fl = chrm, frag_len(rng, "peak")
                    fs = rng.randrange(0, len(chrm) - fl)
                elif u < 0.50:
                    src, fl = chrom, frag_len(rng, "peak")
                    fs = int(rng.choice(bl) + rng.gauss(0, 300) - fl / 2)
                else:
                    src, fl = chrom, frag_len(rng, "bg")
                    fs = rng.randrange(0, L - fl)
                if fs < 0 or fs + fl > len(src):
                    continue
                frag = src[fs:fs + fl]
                if "N" in frag:
                    continue
                prev = (src, fs, fl)
                if rng.random() < 0.5:
                    frag = revcomp(frag)
                r1 = (frag + ADAPTER + "A" * READ_LEN)[:READ_LEN]
                r2 = (revcomp(frag) + ADAPTER + "A" * READ_LEN)[:READ_LEN]
                r1, r2 = mutate(r1, rng), mutate(r2, rng)
                q = "I" * READ_LEN
                name = f"@{sample}_{run}_{n}"
                o1.write(f"{name}/1\n{r1}\n+\n{q}\n")
                o2.write(f"{name}/2\n{r2}\n+\n{q}\n")
                n += 1
        print(f"  {sample} {run}: {a.pairs} pairs")


if __name__ == "__main__":
    main()
