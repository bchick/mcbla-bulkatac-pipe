#!/usr/bin/env python3
"""Head-to-head of this pipeline and nf-core/atacseq on the synthetic .test data.

Run from the repo root after `pixi run test-all` and `pixi run test-nfcore-run`:

    python .test/nfcore/compare_nfcore.py

Both pipelines are scored against the simulated truth (simulate_atac.py):
600 peaks in truth_peaks.bed, 5% PCR duplicates, 8% chrM fragments. Two
kinds of check:

  gate    fails the test (exit 1) when this pipeline is meaningfully worse
          than nf-core against the truth:
            * final consensus peak set: F1 or precision more than 0.05 below
              nf-core's;
            * per-library duplicate rate and chrM fraction: further from the
              simulated value than nf-core's, by more than 0.02.
  report  written to the report only; differences are expected because the
          defaults differ on purpose (IDR-filtered consensus vs a merge of
          every peak call; cutadapt Nextera vs Trim Galore).

Writes .test/results/nfcore_compare/{metrics.tsv,report.md}.
"""

import bisect
import csv
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OURS = os.path.join(HERE, "results")
NF = os.path.join(HERE, "nfcore", "run", "results")
REF = os.path.join(HERE, "data", "ref")
OUT = os.path.join(OURS, "nfcore_compare")

SIM_DUP = 0.05
SIM_MITO = 0.08
TOL_PEAKS = 0.05
TOL_QC = 0.02

metrics = []  # (section, metric, ours, nfcore, kind, status)
failures = []


def fmt(x):
    return "NA" if x is None else (f"{x:.3f}" if isinstance(x, float) else str(x))


def record(section, metric, ours, nfcore, gate=None):
    """gate: None (report) or a bool (True = pass)."""
    kind = "report" if gate is None else "gate"
    status = "-" if gate is None else ("PASS" if gate else "FAIL")
    metrics.append((section, metric, ours, nfcore, kind, status))
    print(f"{status:4} {section}: {metric}  ours={fmt(ours)}  nf-core={fmt(nfcore)}")
    if gate is False:
        failures.append(f"{section}: {metric}")


def need(pattern, what):
    hits = sorted(glob.glob(pattern))
    if not hits:
        sys.exit(f"missing {what}: {pattern}\n(run `pixi run test-all` and `pixi run test-nfcore-run` first)")
    return hits


# --- intervals ----------------------------------------------------------------
def bed(path):
    with open(path) as fh:
        return [(f[0], int(f[1]), int(f[2])) for f in (l.split("\t") for l in fh)
                if len(f) >= 3 and not f[0].startswith(("#", "track"))]


def merged(ivs):
    """Disjoint sorted intervals per chromosome."""
    by = {}
    for c, s, e in sorted(ivs):
        m = by.setdefault(c, [])
        if m and s <= m[-1][1]:
            m[-1][1] = max(m[-1][1], e)
        else:
            m.append([s, e])
    return {c: ([s for s, _ in m], [e for _, e in m]) for c, m in by.items()}


def hit_fraction(a, b):
    """Fraction of intervals in a overlapping any interval in b."""
    if not a:
        return 0.0
    mb = merged(b)
    hit = 0
    for c, s, e in a:
        if c in mb:
            starts, ends = mb[c]
            i = bisect.bisect_left(starts, e) - 1
            hit += i >= 0 and ends[i] > s
    return hit / len(a)


def jaccard(a, b):
    """Base-pair Jaccard of two interval sets."""
    def total(m):
        return sum(e - s for st, en in m.values() for s, e in zip(st, en))
    ma, mb = merged(a), merged(b)
    inter = 0
    for c in ma.keys() & mb.keys():
        (sa, ea), (sb, eb) = ma[c], mb[c]
        i = j = 0
        while i < len(sa) and j < len(sb):
            inter += max(0, min(ea[i], eb[j]) - max(sa[i], sb[j]))
            if ea[i] < eb[j]:
                i += 1
            else:
                j += 1
    union = total(ma) + total(mb) - inter
    return inter / union if union else 0.0


def prf(peaks, truth):
    p, r = hit_fraction(peaks, truth), hit_fraction(truth, peaks)
    return p, r, (2 * p * r / (p + r) if p + r else 0.0)


# --- peaks --------------------------------------------------------------------
truth = bed(os.path.join(REF, "truth_peaks.bed"))
ours_cons = bed(need(os.path.join(OURS, "peaks/consensus/consensus_idr.bed"), "our consensus")[0])
nf_cons = bed(need(os.path.join(NF, "*/merged_library/macs2/narrow_peak/consensus/consensus_peaks.mLb.clN.bed"),
                   "nf-core consensus peaks")[0])

po, ro, fo = prf(ours_cons, truth)
pn, rn, fn = prf(nf_cons, truth)
record("consensus", "peaks", len(ours_cons), len(nf_cons))
record("consensus", "precision vs truth", po, pn, po >= pn - TOL_PEAKS)
record("consensus", "recall vs truth", ro, rn)
record("consensus", "F1 vs truth", fo, fn, fo >= fn - TOL_PEAKS)
record("consensus", "bp Jaccard ours vs nf-core", jaccard(ours_cons, nf_cons), None)
mrp = glob.glob(os.path.join(NF, "*/merged_replicate/macs2/narrow_peak/consensus/consensus_peaks.mRp.clN.bed"))
if mrp:
    nf_mrp = bed(mrp[0])
    record("consensus", "nf-core merged-replicate consensus: peaks / F1 vs truth",
           None, f"{len(nf_mrp)} / {prf(nf_mrp, truth)[2]:.3f}")

# Per condition: our IDR set vs nf-core's merged-replicate peaks (conditions
# with >= 2 replicates; the single-replicate condition has no IDR set).
for path in sorted(glob.glob(os.path.join(OURS, "peaks/idr/*_idr.narrowPeak"))):
    cond = os.path.basename(path)[: -len("_idr.narrowPeak")]
    nf_hits = glob.glob(os.path.join(NF, f"*/merged_replicate/macs2/narrow_peak/{cond}.mRp.clN_peaks.narrowPeak"))
    if not nf_hits:
        continue
    _, _, f_o = prf(bed(path), truth)
    _, _, f_n = prf(bed(nf_hits[0]), truth)
    record(f"condition {cond}", "F1 vs truth (IDR vs merged-replicate peaks)", f_o, f_n)


# --- per-library QC -----------------------------------------------------------
def kv(path):
    with open(path) as fh:
        return dict(l.rstrip("\n").split("\t")[:2] for l in fh if "\t" in l)


def ours_dup(lib):
    with open(os.path.join(OURS, f"qc/markdup/{lib}.markdup.txt")) as fh:
        t = fh.read()
    n = lambda k: int(re.search(rf"^{k}: (\d+)", t, re.M).group(1))
    return (n("DUPLICATE PAIR") + n("DUPLICATE SINGLE")) / n("EXAMINED")


def ours_mito(lib):
    s = kv(os.path.join(OURS, f"qc/filter_stats/{lib}.filter_stats.tsv"))
    return int(s["mito_reads"]) / int(s["total_aligned"])


def nf_dup(lib):
    path = need(os.path.join(NF, f"*/merged_library/picard_metrics/{lib}.mLb.mkD.sorted.MarkDuplicates.metrics.txt"),
                f"nf-core MarkDuplicates metrics for {lib}")[0]
    with open(path) as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh]
    i = next(k for k, l in enumerate(lines) if l[0] == "LIBRARY")
    return float(dict(zip(lines[i], lines[i + 1]))["PERCENT_DUPLICATION"])


def nf_mito(lib):
    path = need(os.path.join(NF, f"*/merged_library/samtools_stats/{lib}.mLb.mkD.sorted.bam.idxstats"),
                f"nf-core idxstats for {lib}")[0]
    with open(path) as fh:
        rows = [l.split("\t") for l in fh]
    mapped = {r[0]: int(r[2]) for r in rows}
    return mapped.get("chrM", 0) / sum(mapped.values())


def frip(lib):
    with open(os.path.join(OURS, f"qc/frip/{lib}.frip.tsv")) as fh:
        return float(next(csv.DictReader(fh, delimiter="\t"))["FRiP"])


def nf_frip(lib):
    hits = glob.glob(os.path.join(NF, f"*/merged_library/macs2/narrow_peak/qc/{lib}.mLb.clN_peaks.FRiP_mqc.tsv"))
    if not hits:
        return None
    with open(hits[0]) as fh:
        return float(fh.read().split()[-1])


libs = sorted(os.path.basename(p)[: -len(".markdup.txt")]
              for p in need(os.path.join(OURS, "qc/markdup/*.markdup.txt"), "our markdup stats"))
for lib in libs:
    do, dn = ours_dup(lib), nf_dup(lib)
    record(lib, f"duplicate rate (simulated {SIM_DUP})", do, dn,
           abs(do - SIM_DUP) <= abs(dn - SIM_DUP) + TOL_QC)
    mo, mn = ours_mito(lib), nf_mito(lib)
    record(lib, f"chrM fraction (simulated {SIM_MITO})", mo, mn,
           abs(mo - SIM_MITO) <= abs(mn - SIM_MITO) + TOL_QC)
    record(lib, "FRiP (own peaks)", frip(lib), nf_frip(lib))


# --- write --------------------------------------------------------------------
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "metrics.tsv"), "w") as fh:
    fh.write("section\tmetric\tours\tnfcore\tkind\tstatus\n")
    for m in metrics:
        fh.write("\t".join(fmt(x) for x in m) + "\n")
with open(os.path.join(OUT, "report.md"), "w") as fh:
    fh.write("# mcbla-bulkatac-pipe vs nf-core/atacseq on the .test dataset\n\n")
    fh.write("Truth: `.test/data/ref/truth_peaks.bed` (600 peaks), 5% PCR duplicates, 8% chrM.\n")
    fh.write("Gates fail the test; report rows are for reading only (the defaults differ on purpose: "
             "our consensus is IDR-filtered, nf-core merges every peak call).\n\n")
    fh.write(f"**{len(failures)} gate failure(s).**\n\n")
    fh.write("| section | metric | ours | nf-core | kind | status |\n|---|---|---|---|---|---|\n")
    for m in metrics:
        fh.write("| " + " | ".join(fmt(x) for x in m) + " |\n")
print(f"\nwrote {OUT}/metrics.tsv and report.md")
print(f"{len(failures)} gate failure(s)")
sys.exit(1 if failures else 0)
