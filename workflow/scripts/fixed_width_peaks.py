"""Summit-centred fixed-width consensus peak set ("iterative overlap").

Follows Corces et al. 2018 (Science 362:eaav1898) and Grandi et al. 2022
(Nat Protoc 17:1518):

  1. per library: every MACS2 summit (narrowPeak start + column 10) becomes a
     `width` bp window; windows running off a chromosome end are dropped.
  2. per library, iterative overlap removal: windows are ranked by MACS2
     -log10(p) (column 8), the most significant is kept and every window
     overlapping it is discarded, and so on down the list.
  3. scores become score-per-million (SPM): -log10(p) / (sum over the
     library's kept windows / 1e6), so libraries of different depth and peak
     number are comparable.
  4. windows with SPM >= min_spm from all libraries are pooled and the same
     iterative overlap removal is run on SPM.
  5. a kept window is reproducible if windows from >= min_samples distinct
     libraries (SPM >= min_spm) overlap it.

All windows have the same width, so two overlap exactly when their summits are
less than `width` bp apart; positions are bucketed by `width`, so each check
looks at three buckets only.

Output (BED, sorted): chrom, start, end, name, SPM, number of supporting
libraries. Blacklist filtering is done by the rule afterwards.
"""

import sys
from collections import defaultdict

sm = snakemake  # noqa: F821  (injected by Snakemake)
sys.stderr = open(sm.log[0], "w")

width = int(sm.params.width)
half = width // 2
min_spm = float(sm.params.min_spm)
min_samples = int(sm.params.min_samples)

sizes = {}
with open(sm.input.sizes) as fh:
    for line in fh:
        chrom, size = line.split("\t")[:2]
        sizes[chrom] = int(size)


def overlaps(buckets, chrom, summit):
    """Kept summits (at most one per bucket) closer than `width` to summit."""
    b = summit // width
    for k in (b - 1, b, b + 1):
        other = buckets.get((chrom, k))
        if other is not None and abs(other - summit) < width:
            return True
    return False


def iterative_overlap(peaks):
    """peaks: [(score, chrom, summit, ...)]; returns the kept subset."""
    kept, buckets = [], {}
    for peak in sorted(peaks, key=lambda p: -p[0]):
        _score, chrom, summit = peak[:3]
        if not overlaps(buckets, chrom, summit):
            buckets[(chrom, summit // width)] = summit
            kept.append(peak)
    return kept


pooled = []
for lib, path in zip(sm.params.libs, sm.input.peaks):
    windows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom, start, offset, pval = f[0], int(f[1]), int(f[9]), float(f[7])
            summit = start + offset
            if chrom not in sizes or summit - half < 0 or summit + half + 1 > sizes[chrom]:
                continue
            windows.append((pval, chrom, summit))
    kept = iterative_overlap(windows)
    total = sum(p[0] for p in kept)
    n_pass = 0
    for pval, chrom, summit in kept:
        spm = pval / (total / 1e6) if total > 0 else 0.0
        if spm >= min_spm:
            pooled.append((spm, chrom, summit, lib))
            n_pass += 1
    print(f"{lib}: {len(windows)} windows, {len(kept)} after iterative overlap, "
          f"{n_pass} with SPM >= {min_spm}", file=sys.stderr)

# supporting libraries per kept window: all pooled windows within `width`
index = defaultdict(list)
for spm, chrom, summit, lib in pooled:
    index[(chrom, summit // width)].append((summit, lib))

consensus = iterative_overlap(pooled)
rows = []
for spm, chrom, summit, _lib in consensus:
    b = summit // width
    libs = {
        lib
        for k in (b - 1, b, b + 1)
        for other, lib in index.get((chrom, k), ())
        if abs(other - summit) < width
    }
    if len(libs) >= min_samples:
        rows.append((chrom, summit - half, summit + half + 1, spm, len(libs)))

rows.sort(key=lambda r: (r[0], r[1]))
with open(sm.output[0], "w") as out:
    for i, (chrom, start, end, spm, n) in enumerate(rows, 1):
        out.write(f"{chrom}\t{start}\t{end}\tfixed_{i}\t{spm:.2f}\t{n}\n")

print(f"pooled: {len(pooled)} windows, {len(consensus)} after iterative overlap, "
      f"{len(rows)} supported by >= {min_samples} libraries", file=sys.stderr)
