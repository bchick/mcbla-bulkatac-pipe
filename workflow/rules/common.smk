# Shared helpers: config and samplesheet loading, validation, input-mode
# resolution (fastq vs nf-core/atacseq outdir) and target lists.
#
# All Python functions live here, so the rule files contain only rules.

import os
import re
import sys
import glob
import itertools
from pathlib import Path

import pandas as pd
from snakemake.utils import validate

# ---------------------------------------------------------------------------
# Config + samplesheet
# ---------------------------------------------------------------------------
validate(config, schema="../schemas/config.schema.yaml")

INPUT_MODE = config["input_mode"]
MODULES = config["modules"]


def _warn(msg):
    print(f"[mcbla-bulkatac-pipe] WARNING: {msg}", file=sys.stderr)


def _read_table(path):
    return pd.read_csv(path, sep=None, engine="python", dtype=str, comment="#").fillna(
        ""
    )


samples_raw = _read_table(config["samples"])
# `condition` is optional in the nf-core samplesheet; default it to `sample`.
if "condition" not in samples_raw.columns:
    samples_raw["condition"] = samples_raw["sample"]
samples_raw["condition"] = [
    c if c else s for c, s in zip(samples_raw["condition"], samples_raw["sample"])
]
for col in ("fastq_1", "fastq_2", "treatment", "time"):
    if col not in samples_raw.columns:
        samples_raw[col] = ""
validate(samples_raw, schema="../schemas/samples.schema.yaml")
if INPUT_MODE == "fastq":
    _nofq = samples_raw[(samples_raw["fastq_1"] == "") | (samples_raw["fastq_2"] == "")]
    if len(_nofq):
        raise ValueError(
            "input_mode=fastq needs paired-end fastq_1 and fastq_2 for every row; "
            f"missing for sample(s): {', '.join(sorted(set(_nofq['sample'])))}"
        )

# nf-core/atacseq names each library <sample>_REP<replicate>; rows that share
# sample + replicate are sequencing runs of one library (merged before trimming).
samples_raw["lib"] = [
    f"{s}_REP{r}" for s, r in zip(samples_raw["sample"], samples_raw["replicate"])
]

_per_lib = samples_raw.groupby("lib", sort=False).agg(
    {
        "sample": "first",
        "replicate": "first",
        "condition": lambda x: ";".join(sorted(set(x))),
        "treatment": lambda x: ";".join(sorted(set(x))),
        "time": lambda x: ";".join(sorted(set(x))),
    }
)
for col in ("condition", "treatment", "time"):
    bad = _per_lib[_per_lib[col].str.contains(";")]
    if len(bad):
        raise ValueError(
            f"Samplesheet: runs of the same library disagree on '{col}': "
            f"{', '.join(bad.index)}"
        )
LIBRARIES = _per_lib
LIBS = list(LIBRARIES.index)
CONDITIONS = list(dict.fromkeys(LIBRARIES["condition"]))
LIBS_BY_COND = {
    c: list(LIBRARIES.index[LIBRARIES["condition"] == c]) for c in CONDITIONS
}

contrasts = _read_table(config["contrasts"]) if config.get("contrasts") else None
if contrasts is not None and len(contrasts):
    validate(contrasts, schema="../schemas/contrasts.schema.yaml")
    _unknown = set(contrasts["group1"]) | set(contrasts["group2"])
    _unknown -= set(CONDITIONS)
    if _unknown:
        raise ValueError(
            f"contrasts.tsv references unknown conditions: {sorted(_unknown)}"
        )
    CONTRASTS = list(contrasts["label"])
    if len(set(CONTRASTS)) != len(CONTRASTS):
        raise ValueError("contrasts.tsv: labels must be unique")
else:
    CONTRASTS = []

# Conditions excluded from count-based statistics (diff, normcheck,
# timecourse, chromvar). Peaks, QC and footprinting still use every condition.
EXCLUDED = set(config.get("exclude_conditions") or [])
STAT_LIBS = [l for l in LIBS if LIBRARIES.loc[l, "condition"] not in EXCLUDED]
if CONTRASTS:
    _hit = contrasts[
        contrasts["group1"].isin(EXCLUDED) | contrasts["group2"].isin(EXCLUDED)
    ]
    if len(_hit):
        raise ValueError(
            "contrasts.tsv uses excluded conditions: " + ", ".join(_hit["label"])
        )


# ---------------------------------------------------------------------------
# Wildcard constraints
# ---------------------------------------------------------------------------
def _alt(values):
    return "|".join(re.escape(v) for v in values) if values else "__none__"


wildcard_constraints:
    lib=_alt(LIBS),
    condition=_alt(CONDITIONS),
    peakset="individual|merged_relaxed|merged_stringent",
    name=_alt(LIBS + CONDITIONS),


# ---------------------------------------------------------------------------
# Reference helpers
# ---------------------------------------------------------------------------
REF = config["reference"]
BLACKLIST = REF.get("blacklist") or ""
USE_BLACKLIST = bool(BLACKLIST)
MITO = REF.get("mito_chrom", "chrM")


def bowtie2_index_prefix():
    """User-supplied index prefix, or one built by the pipeline from the FASTA."""
    if REF.get("bowtie2_index"):
        return REF["bowtie2_index"]
    return "results/reference/bowtie2/genome"


def bowtie2_index_files():
    prefix = bowtie2_index_prefix()
    if REF.get("bowtie2_index"):
        # small (.bt2) or large (.bt2l) index
        found = glob.glob(f"{prefix}.1.bt2*")
        ext = ".bt2l" if found and found[0].endswith(".bt2l") else ".bt2"
    else:
        ext = ".bt2"
    return [f"{prefix}.{s}{ext}" for s in ("1", "2", "3", "4", "rev.1", "rev.2")]


def blacklist_input():
    return [BLACKLIST] if USE_BLACKLIST else []


# ---------------------------------------------------------------------------
# FASTQ-mode inputs
# ---------------------------------------------------------------------------
def _fastq_path(p):
    base = config.get("fastq_dir") or ""
    if not p:
        return p
    if base and not os.path.isabs(p):
        return os.path.join(base, p)
    return p


def lib_runs(lib):
    rows = samples_raw[samples_raw["lib"] == lib]
    return [(_fastq_path(r.fastq_1), _fastq_path(r.fastq_2)) for r in rows.itertuples()]


def trim_inputs(wildcards):
    runs = lib_runs(wildcards.lib)
    if len(runs) == 1:
        return {"r1": runs[0][0], "r2": runs[0][1]}
    return {
        "r1": f"results/fastq/{wildcards.lib}_R1.merged.fastq.gz",
        "r2": f"results/fastq/{wildcards.lib}_R2.merged.fastq.gz",
    }


def merge_fastq_inputs(wildcards):
    idx = 0 if wildcards.read == "1" else 1
    return [r[idx] for r in lib_runs(wildcards.lib)]


# ---------------------------------------------------------------------------
# nf-core/atacseq-mode inputs
# ---------------------------------------------------------------------------
NFCORE = config.get("nfcore", {})
NFCORE_ALIGNERS = ["bowtie2", "bwa", "chromap", "star"]


def nfcore_aligner():
    outdir = NFCORE.get("outdir") or ""
    want = NFCORE.get("aligner", "auto")
    if want != "auto":
        return want
    for a in NFCORE_ALIGNERS:
        if os.path.isdir(os.path.join(outdir, a, "merged_library")):
            return a
    raise ValueError(
        f"input_mode=nfcore: no <aligner>/merged_library/ directory found under "
        f"nfcore.outdir='{outdir}' (looked for {', '.join(NFCORE_ALIGNERS)})."
    )


def nfcore_lib_bam(wildcards):
    d = os.path.join(NFCORE["outdir"], nfcore_aligner(), "merged_library")
    return os.path.join(d, f"{wildcards.lib}.mLb.clN.sorted.bam")


def nfcore_mrp_bam(condition):
    """nf-core merged-replicate BAM for a condition, if usable, else None.

    Only used when every library of the condition comes from a single nf-core
    `sample` group (nf-core merges replicates per `sample`).
    """
    if INPUT_MODE != "nfcore" or not NFCORE.get("use_merged_replicate", False):
        return None
    groups = set(LIBRARIES.loc[LIBS_BY_COND[condition], "sample"])
    if len(groups) != 1:
        return None
    d = os.path.join(NFCORE["outdir"], nfcore_aligner(), "merged_replicate")
    return os.path.join(d, f"{groups.pop()}.mRp.clN.sorted.bam")


# ---------------------------------------------------------------------------
# BAM accessors (identical downstream of both input modes)
# ---------------------------------------------------------------------------
def lib_bam(lib):
    return f"results/bam/{lib}.final.bam"


def all_lib_bams(libs=None):
    return [lib_bam(l) for l in (libs or LIBS)]


def cond_libs(wildcards):
    return [lib_bam(l) for l in LIBS_BY_COND[wildcards.condition]]


def cond_bam_inputs(wildcards):
    """nf-core merged-replicate BAM if configured and usable, else rep BAMs."""
    mrp = nfcore_mrp_bam(wildcards.condition)
    return [mrp] if mrp else cond_libs(wildcards)


# Shared by the per-replicate and merged MACS2 rules (2.1_atac_peaks.sh):
# MACS2 -f BAMPE --keep-dup all, then blacklist-filter the narrowPeak in place.
MACS2_SHELL = """
(
set -euo pipefail
macs2 callpeak -t {input.bam} -f BAMPE -g {params.gsize} --keep-dup all \
    {params.threshold} --outdir {params.outdir} --name {params.name} {params.extra}
if [ -n "{input.blacklist}" ]; then
    raw=$(wc -l < {output.peaks})
    mv {output.peaks} {output.peaks}.unfiltered
    bedtools intersect -v -a {output.peaks}.unfiltered -b {input.blacklist} > {output.peaks}
    rm -f {output.peaks}.unfiltered
    echo "Blacklist: $raw -> $(wc -l < {output.peaks})"
fi
) > {log} 2>&1
"""


# ---------------------------------------------------------------------------
# IDR pairing
# ---------------------------------------------------------------------------
IDR = config["idr"]


def idr_pairs(condition):
    libs = LIBS_BY_COND[condition]
    if len(libs) < 2:
        return []
    if len(libs) > 2 and IDR.get("multi_rep_strategy") == "first_two":
        _warn(
            f"IDR: condition '{condition}' has {len(libs)} replicates; "
            "multi_rep_strategy=first_two uses only the first two."
        )
        return [(libs[0], libs[1])]
    return list(itertools.combinations(libs, 2))


IDR_CONDITIONS = [c for c in CONDITIONS if len(LIBS_BY_COND[c]) >= 2]
SINGLE_REP_CONDITIONS = [c for c in CONDITIONS if len(LIBS_BY_COND[c]) < 2]
if SINGLE_REP_CONDITIONS and MODULES.get("idr", True):
    _warn(
        "IDR skipped for single-replicate condition(s): "
        + ", ".join(SINGLE_REP_CONDITIONS)
        + f" (fallback: {IDR.get('single_rep_fallback', 'stringent')})."
    )


def idr_pair_files(wildcards):
    return [
        f"results/peaks/idr/pairs/{wildcards.condition}/{a}__{b}.idr.narrowPeak"
        for a, b in idr_pairs(wildcards.condition)
    ]


def idr_scaled_threshold():
    # IDR column 5 = min(int(-125 * log2(IDR)), 1000); 0.05 -> 540 (as in source)
    import math

    return min(int(-125 * math.log2(float(IDR["threshold"]))), 1000)


def condition_reproducible_peaks(condition):
    if len(LIBS_BY_COND[condition]) >= 2:
        return f"results/peaks/idr/{condition}_idr.narrowPeak"
    if IDR.get("single_rep_fallback", "stringent") == "stringent":
        return f"results/peaks/merged_stringent/{condition}_peaks.narrowPeak"
    return None


def consensus_inputs(wildcards):
    files = [condition_reproducible_peaks(c) for c in CONDITIONS]
    return [f for f in files if f]


# ---------------------------------------------------------------------------
# Downstream peak-set selection
# ---------------------------------------------------------------------------
def consensus_bed():
    """Peak set used for counting (timecourse, chromVAR, optional DiffBind)."""
    if MODULES.get("idr", True) and config["consensus"]["source"] == "idr":
        return "results/peaks/consensus/consensus_idr.bed"
    return "results/peaks/consensus/union_stringent.bed"


def footprint_peaks():
    if config["footprint"]["peaks"] == "consensus":
        return consensus_bed()
    return "results/peaks/consensus/union_stringent.bed"


def peakset_members(wildcards):
    if wildcards.peakset == "individual":
        return [f"results/peaks/individual/{l}.stats.tsv" for l in LIBS]
    return [f"results/peaks/{wildcards.peakset}/{c}.stats.tsv" for c in CONDITIONS]


def peak_stats_bam(wildcards):
    if wildcards.peakset == "individual":
        return lib_bam(wildcards.name)
    return f"results/bam/merged/{wildcards.name}.merged.bam"


# ---------------------------------------------------------------------------
# Time-course series
# ---------------------------------------------------------------------------
TC = config["timecourse"]
HAS_TIME = bool((LIBRARIES["time"] != "").all()) and bool(
    (LIBRARIES["treatment"] != "").all()
)
RUN_TIMECOURSE = bool(MODULES.get("timecourse", False)) and HAS_TIME
if MODULES.get("timecourse", False) and not HAS_TIME:
    _warn(
        "modules.timecourse is on, but the samplesheet lacks complete "
        "'treatment' and 'time' columns; the time-course module is skipped."
    )


def timecourse_series():
    if not RUN_TIMECOURSE:
        return []
    base = set(TC.get("baseline_treatments") or [])
    keep = LIBRARIES.loc[STAT_LIBS]
    return [t for t in dict.fromkeys(keep["treatment"]) if t not in base]


def multiqc_inputs(wildcards):
    files = expand("results/qc/flagstat/{lib}.flagstat.txt", lib=LIBS)
    if INPUT_MODE == "fastq":
        files += expand("logs/align/{lib}.cutadapt.log", lib=LIBS)
        files += expand("logs/align/{lib}.bowtie2.log", lib=LIBS)
        files += expand("results/qc/markdup/{lib}.markdup.txt", lib=LIBS)
        files.append("results/qc/alignment_qc_report.tsv")
    if MODULES.get("qc", True):
        files += [
            "results/qc/deeptools/tss_enrichment_profile.tab",
            "results/qc/deeptools/fragment_size_table.tsv",
            "results/qc/deeptools/fragment_size_raw.tsv",
            "results/qc/deeptools/correlation_spearman.tsv",
            "results/qc/deeptools/pca_data.tsv",
            "results/qc/deeptools/fingerprint_metrics.tsv",
            "results/qc/deeptools/fingerprint_counts.tsv",
        ]
    files += expand("results/peaks/individual/{lib}_peaks.xls", lib=LIBS)
    files += [
        "results/peaks/individual/peak_summary.tsv",
        "results/peaks/merged_stringent/peak_summary.tsv",
    ]
    if MODULES.get("idr", True):
        files.append("results/peaks/idr/idr_summary.tsv")
    return files


# ---------------------------------------------------------------------------
# Motifs
# ---------------------------------------------------------------------------
def tobias_motifs():
    return config["footprint"].get("motifs") or "results/reference/jaspar_motifs.jaspar"


# ---------------------------------------------------------------------------
# Resources
# ---------------------------------------------------------------------------
def threads(key, default=4):
    return int(config.get("threads", {}).get(key, default))


# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------
def core_targets():
    t = all_lib_bams()
    t += [f"results/bam/{l}.final.bam.bai" for l in LIBS]
    if INPUT_MODE == "fastq":
        t.append("results/qc/alignment_qc_report.tsv")
        t += [f"results/qc/fragment_sizes/{l}_fragment_sizes.tsv" for l in LIBS]
    t += expand("results/bam/merged/{c}.merged.bam", c=CONDITIONS)
    t += [
        "results/peaks/individual/peak_summary.tsv",
        "results/peaks/merged_relaxed/peak_summary.tsv",
        "results/peaks/merged_stringent/peak_summary.tsv",
        "results/peaks/consensus/union_stringent.bed",
    ]
    if MODULES.get("idr", True):
        t += [
            "results/peaks/idr/idr_summary.tsv",
            "results/peaks/consensus/consensus_idr.bed",
        ]
    return t


def all_targets():
    t = core_targets()
    if MODULES.get("qc", True):
        t += expand("results/bigwig/{l}.bw", l=LIBS)
        t += [
            "results/qc/deeptools/tss_enrichment_profile.png",
            "results/qc/deeptools/fragment_size_distribution.png",
            "results/qc/deeptools/correlation_spearman.png",
            "results/qc/deeptools/pca_plot.png",
            "results/qc/deeptools/fingerprint.png",
        ]
    if MODULES.get("multiqc", True):
        t.append("results/qc/multiqc/multiqc_report.html")
    if MODULES.get("diff", True) and CONTRASTS:
        t += [
            "results/diff/depth/dba_analyzed.rds",
            "results/diff/depth/summary.tsv",
        ]
    if MODULES.get("normcheck", True) and CONTRASTS:
        t += [
            "results/normcheck/norm_comparison.tsv",
            "results/normcheck/norm_comparison_barplot.pdf",
        ]
    for s in timecourse_series():
        t += [
            f"results/timecourse/{s}/lrt_results.tsv",
            f"results/timecourse/{s}/degpatterns_clusters.tsv",
        ]
    if MODULES.get("chromvar", True):
        t += [
            "results/chromvar/deviation_zscores.tsv",
            "results/chromvar/variability.tsv",
        ]
    if MODULES.get("footprint", True):
        t.append("results/footprint/bindetect/bindetect_results.txt")
    return t


wildcard_constraints:
    series=_alt(timecourse_series()),
