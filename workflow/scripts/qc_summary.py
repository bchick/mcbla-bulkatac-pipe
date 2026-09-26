"""Per-library QC summary with PASS / WARN / FAIL flags.

Collects fragments (flagstat), alignment rate, NRF / PBC1 / PBC2 (fastq mode,
from alignment_qc_report.tsv), FRiP over the condition's merged stringent peaks
and the TSS enrichment score, and flags each metric against qc.thresholds
([pass, warn]: >= pass PASS, >= warn WARN, else FAIL). Missing metrics are NA
and not flagged. Overall is the worst flag of the library.
"""

import csv
import re
import sys

sm = snakemake  # noqa: F821  (injected by Snakemake)
sys.stderr = open(sm.log[0], "w")

METRICS = ["fragments", "aligned_pct", "nrf", "pbc1", "pbc2", "frip", "tss_enrichment"]
COLUMNS = {
    "fragments": "Fragments",
    "aligned_pct": "Aligned_Pct",
    "nrf": "NRF",
    "pbc1": "PBC1",
    "pbc2": "PBC2",
    "frip": "FRiP",
    "tss_enrichment": "TSS_Enrichment",
}
RANK = {"PASS": 0, "WARN": 1, "FAIL": 2}
thresholds = sm.params.thresholds
libs = list(sm.params.libs)


def read_rows(path):
    with open(path) as fh:
        return {row["Sample"]: row for row in csv.DictReader(fh, delimiter="\t")}


def number(x):
    try:
        return float(str(x).rstrip("%"))
    except (TypeError, ValueError):
        return None


values = {lib: dict.fromkeys(METRICS) for lib in libs}

for lib, path in zip(libs, sm.input.flagstat):
    with open(path) as fh:
        m = re.search(r"^(\d+) \+ \d+ paired in sequencing", fh.read(), flags=re.MULTILINE)
    values[lib]["fragments"] = int(m.group(1)) // 2 if m else None

if sm.input.get("report"):
    report = read_rows(sm.input.report)
    for lib in libs:
        row = report.get(lib, {})
        for key in ("aligned_pct", "nrf", "pbc1", "pbc2"):
            values[lib][key] = number(row.get(COLUMNS[key]))

for path in sm.input.frip:
    for lib, row in read_rows(path).items():
        values[lib]["frip"] = number(row["FRiP"])

for path in sm.input.get("tss", []):
    for lib, row in read_rows(path).items():
        values[lib]["tss_enrichment"] = number(row["TSS_Enrichment"])


def flag(metric, value):
    if value is None or metric not in thresholds:
        return "NA"
    pass_at, warn_at = thresholds[metric]
    if value >= pass_at:
        return "PASS"
    return "WARN" if value >= warn_at else "FAIL"


def fmt(metric, value):
    if value is None:
        return "NA"
    if metric == "fragments":
        return str(int(value))
    if metric == "aligned_pct":
        return f"{value:.2f}"
    return f"{value:.4g}"


header = ["Sample", "Condition"]
for m in METRICS:
    header += [COLUMNS[m], f"{COLUMNS[m]}_flag"]
header += ["Overall", "Not_Passing"]

with open(sm.output[0], "w") as out:
    out.write("\t".join(header) + "\n")
    for lib, cond in zip(libs, sm.params.conditions):
        row = [lib, cond]
        flags = {}
        for m in METRICS:
            flags[m] = flag(m, values[lib][m])
            row += [fmt(m, values[lib][m]), flags[m]]
        scored = [f for f in flags.values() if f != "NA"]
        overall = max(scored, key=RANK.get) if scored else "NA"
        not_passing = [f"{COLUMNS[m]}:{f}" for m, f in flags.items() if f in ("WARN", "FAIL")]
        row += [overall, ";".join(not_passing) or "-"]
        out.write("\t".join(row) + "\n")
        if overall == "FAIL":
            print(f"QC FAIL {lib}: {';'.join(not_passing)}", file=sys.stderr)
