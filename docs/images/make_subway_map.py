"""Generate the nf-core style subway map SVGs (light + dark) used in the README.

Usage: python docs/images/make_subway_map.py docs/images
"""
import sys
from html import escape
from pathlib import Path

W, H = 1520, 660
Y = 310  # trunk
Y_NFC, Y_TC, Y_NORM, Y_CV, Y_FP, Y_QC = 200, 100, 190, 250, 410, 510

LINES = {
    "core":  ("#24B064", "Processing: FASTQ to consensus peak sets"),
    "nfc":   ("#8C8C8C", "nf-core/atacseq BAM input"),
    "qc":    ("#F2B138", "QC"),
    "diff":  ("#1F6FEB", "Differential accessibility"),
    "norm":  ("#8E44AD", "Normalization check"),
    "tc":    ("#E8702A", "Time course"),
    "cv":    ("#16A2B8", "TF activity (chromVAR)"),
    "fp":    ("#D83B3B", "TF footprinting (TOBIAS)"),
}

X0, X_BAM, X_MERGED, X_PEAKS, X_COUNTS = 90, 590, 790, 1000, 1090
# fan-out from the counts hub climbs one 45-degree diagonal
FAN = [(X_COUNTS, Y), (X_COUNTS + 60, Y_CV), (X_COUNTS + 120, Y_NORM), (X_COUNTS + 210, Y_TC)]
PATHS = {
    "nfc":  [(X0, Y_NFC), (X_BAM - (Y - Y_NFC), Y_NFC), (X_BAM, Y)],
    "qc":   [(X_BAM, Y), (X_BAM + 60, Y + 60), (X_BAM + 60, Y_QC - 60), (X_BAM + 120, Y_QC), (1025, Y_QC)],
    "fp":   [(X_PEAKS, Y), (X_PEAKS + (Y_FP - Y), Y_FP), (1370, Y_FP)],
    "tc":   FAN + [(1460, Y_TC)],
    "norm": FAN[:3] + [(1460, Y_NORM)],
    "cv":   FAN[:2] + [(1270, Y_CV)],
    "diff": [(X_COUNTS, Y), (1450, Y)],
    "core": [(X0, Y), (X_COUNTS, Y)],
}

# stations: (x, y, label, label side, kind); kind: stop | hub | start | end
S = [
    (X0, Y, "FASTQ", "below", "start"),
    (190, Y, "cutadapt\n(Nextera trim)", "below", "stop"),
    (295, Y, "Bowtie 2", "below", "stop"),
    (400, Y, "Filter\nMAPQ30, pairs,\nchrM, blacklist", "below", "stop"),
    (500, Y, "Picard\nMarkDuplicates", "below", "stop"),
    (X_BAM, Y, "Final BAMs", "above", "hub"),
    (X0, Y_NFC, "nf-core/atacseq\n*.mLb.clN.bam", "above", "start"),
    (680, Y, "MACS2\nper replicate", "above", "stop"),
    (X_MERGED, Y, "MACS2 merged\nrelaxed + stringent", "above", "stop"),
    (900, Y, "IDR\nENCODE 0.05", "above", "stop"),
    (X_PEAKS, Y, "Consensus\npeak sets", "above", "stop"),
    (X_COUNTS, Y, "Counts over\nchosen peaks", "below", "hub"),
    # QC
    (710, Y_QC, "deepTools\nbigWig", "below", "stop"),
    (790, Y_QC, "TSS\nenrichment", "below", "stop"),
    (870, Y_QC, "Fingerprint", "below", "stop"),
    (950, Y_QC, "Correlation\n+ PCA", "above", "stop"),
    (1025, Y_QC, "MultiQC", "above", "end"),
    # footprinting (merged BAMs + chosen peak set)
    (1150, Y_FP, "ATACorrect", "below", "stop"),
    (1260, Y_FP, "ScoreBigwig", "below", "stop"),
    (1370, Y_FP, "BINDetect", "below", "end"),
    # differential
    (1270, Y, "DiffBind +\nDESeq2 contrasts", "below", "stop"),
    (1450, Y, "Tables\n+ MA plots", "below", "end"),
    # chromVAR
    (1270, Y_CV, "chromVAR\n(JASPAR2020)", "right", "end"),
    # normalization check
    (1345, Y_NORM, "csaw bins,\nquantile + limma", "above", "stop"),
    (1460, Y_NORM, "Sensitivity\nverdict", "above", "end"),
    # time course
    (1345, Y_TC, "DESeq2 LRT\nover time", "above", "stop"),
    (1460, Y_TC, "degPatterns\nclusters", "above", "end"),
]

STAGES = [(40, X_BAM + 30, "1", "Preprocessing"), (X_BAM + 30, X_PEAKS + 52, "2", "Peaks & QC"),
          (X_PEAKS + 52, W - 20, "3", "Analyses (opt-in)")]


def svg(theme):
    fg = "#1F2328" if theme == "light" else "#E6EDF3"
    muted = "#59636E" if theme == "light" else "#9198A1"
    band = "#F6F8FA" if theme == "light" else "#161B22"
    stfill = "#FFFFFF"
    stroke = "#1F2328"
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
         f'font-family="Helvetica Neue, Helvetica, Arial, sans-serif">',
         '<title>mcbla-bulkatac-pipe subway map</title>']
    # stage bands
    for x0, x1, n, name in STAGES:
        o.append(f'<rect x="{x0}" y="20" width="{x1 - x0 - 8}" height="{H - 100}" rx="14" fill="{band}"/>')
        o.append(f'<circle cx="{x0 + 22}" cy="44" r="13" fill="{fg}"/>'
                 f'<text x="{x0 + 22}" y="49" font-size="14" font-weight="700" fill="{band}" text-anchor="middle">{n}</text>'
                 f'<text x="{x0 + 42}" y="49" font-size="15" font-weight="700" fill="{fg}">{escape(name)}</text>')
    # lines
    for key in ["qc", "fp", "tc", "norm", "cv", "diff", "nfc", "core"]:
        pts = " ".join(f"{x},{y}" for x, y in PATHS[key])
        cap = ' stroke-dasharray="16 8" stroke-linecap="butt"' if key == "nfc" else ' stroke-linecap="round"'
        o.append(f'<polyline points="{pts}" fill="none" stroke="{LINES[key][0]}" stroke-width="9" '
                 f'stroke-linejoin="round"{cap}/>')
    # stations
    for x, y, label, side, kind in S:
        if kind == "hub":
            o.append(f'<rect x="{x - 13}" y="{y - 13}" width="26" height="26" rx="13" fill="{stfill}" stroke="{stroke}" stroke-width="3.5"/>')
        elif kind in ("start", "end"):
            o.append(f'<rect x="{x - 10}" y="{y - 12}" width="20" height="24" rx="3" fill="{stfill}" stroke="{stroke}" stroke-width="3"/>'
                     f'<line x1="{x - 5}" y1="{y - 4}" x2="{x + 5}" y2="{y - 4}" stroke="{stroke}" stroke-width="2"/>'
                     f'<line x1="{x - 5}" y1="{y + 2}" x2="{x + 5}" y2="{y + 2}" stroke="{stroke}" stroke-width="2"/>')
        else:
            o.append(f'<circle cx="{x}" cy="{y}" r="9" fill="{stfill}" stroke="{stroke}" stroke-width="3"/>')
        lines = label.split("\n")
        lh = 15
        if side == "below":
            ys, anchor, lx = [y + 32 + i * lh for i in range(len(lines))], "middle", x
        elif side == "above":
            ys, anchor, lx = [y - 22 - (len(lines) - 1 - i) * lh for i in range(len(lines))], "middle", x
        else:
            ys, anchor, lx = [y + 5 + (i - (len(lines) - 1) / 2) * lh for i in range(len(lines))], "start", x + 18
        for i, (t, ty) in enumerate(zip(lines, ys)):
            weight = "700" if i == 0 else "400"
            col = fg if i == 0 else muted
            size = 13 if i == 0 else 11.5
            o.append(f'<text x="{lx}" y="{ty:.0f}" font-size="{size}" font-weight="{weight}" fill="{col}" text-anchor="{anchor}">{escape(t)}</text>')
    # legend
    ly = H - 58
    order = ["core", "nfc", "qc", "diff", "norm", "tc", "cv", "fp"]
    x = 40
    for i, key in enumerate(order):
        col, name = LINES[key]
        cx = 40 + (i % 4) * 350
        cy = ly + (i // 4) * 28
        cap = ' stroke-dasharray="10 5" stroke-linecap="butt"' if key == "nfc" else ' stroke-linecap="round"'
        o.append(f'<line x1="{cx}" y1="{cy}" x2="{cx + 36}" y2="{cy}" stroke="{col}" stroke-width="8"{cap}/>'
                 f'<text x="{cx + 50}" y="{cy + 5}" font-size="13" fill="{fg}">{escape(name)}</text>')
    o.append("</svg>")
    return "\n".join(o)


out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
for t in ("light", "dark"):
    (out / f"subway_map_{t}.svg").write_text(svg(t) + "\n")
