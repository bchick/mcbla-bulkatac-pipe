#!/usr/bin/env python3
"""Set up an ATAC-seq project config from the lab's reference manifest.

Asks which genome and which blacklist, then writes <project>/project.yaml with the matching paths from the
resource manifest (/data/resource/manifest.yaml on the Salk server), ready for

    pixi run snakemake -s workflow/Snakefile --directory <project> \
        --configfile <project>/project.yaml --profile profiles/local -n

Run it interactively (`pixi run init`), or give every answer as a flag, which is
how agents should call it after asking the user:

    pixi run init --dir /path/to/project --genome hg38 --blacklist encode_v2

`--list` prints the choices the manifest offers. The manifest is found at
--manifest, else $MCBLA_RESOURCE_MANIFEST, else /data/resource/manifest.yaml.
"""

import argparse
import datetime
import json
import os
import shutil
import sys

import yaml

DEFAULT_MANIFEST = "/data/resource/manifest.yaml"
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
# pixi runs tasks from the repo root; INIT_CWD is where the user typed the command
CWD = os.environ.get("INIT_CWD") or os.getcwd()


def die(msg, code=1):
    print(f"init: {msg}", file=sys.stderr)
    sys.exit(code)


def find_manifest(arg):
    for p in (arg, os.environ.get("MCBLA_RESOURCE_MANIFEST"), DEFAULT_MANIFEST):
        if p and os.path.exists(p):
            return p
    die(
        "no resource manifest found (looked at --manifest, $MCBLA_RESOURCE_MANIFEST, "
        f"{DEFAULT_MANIFEST}). Outside the Salk server, copy config/salk_example.yaml "
        "and fill in your reference paths by hand (see README)."
    )


def usable(entry, allow_unverified):
    return entry.get("status", "ok") == "ok" or allow_unverified


def choices(genome, allow_unverified):
    """Blacklist keys offered for one genome."""
    return {
        k: v
        for k, v in (genome.get("blacklists") or {}).items()
        if usable(v, allow_unverified)
    }


def default_blacklist(bl):
    for k, v in bl.items():
        if v.get("default"):
            return k
    return next(iter(bl), "none")


def resolve_genome(genomes, answer):
    a = answer.lower()
    for key, g in genomes.items():
        if a == key.lower() or a in (x.lower() for x in g.get("aliases", [])):
            return key
    die(f"unknown genome '{answer}'; choose one of: {', '.join(genomes)}")


def ask(question, options, default=None):
    """Numbered-choice prompt. options: list of (value, description)."""
    print(f"\n{question}")
    for i, (val, desc) in enumerate(options, 1):
        mark = "  (default)" if val == default else ""
        print(f"  {i}. {val:<16} {desc}{mark}")
    while True:
        raw = input(
            f"Choice [1-{len(options)}{', Enter = default' if default else ''}]: "
        ).strip()
        if not raw and default:
            return default
        if raw.isdigit() and 1 <= int(raw) <= len(options):
            return options[int(raw) - 1][0]
        for val, _ in options:
            if raw.lower() == val.lower():
                return val
        if os.path.sep in raw or raw.endswith(".bed"):  # custom BED path
            return raw
        print("  Please enter one of the numbers above.")


def ask_text(question, default):
    raw = input(f"\n{question} [{default}]: ").strip()
    return raw or default


def list_choices(man, genomes, allow_unverified, as_json):
    out = {}
    for key, g in genomes.items():
        bl = choices(g, allow_unverified)
        out[key] = {
            "description": g.get("description", ""),
            "blacklists": {k: v.get("description", "") for k, v in bl.items()},
            "default_blacklist": default_blacklist(bl),
        }
    if as_json:
        print(json.dumps(out, indent=2))
        return
    print(
        f"Manifest version {man['manifest_version']} (updated {man.get('updated')})\n"
    )
    for key, o in out.items():
        print(f"{key}: {o['description']}")
        bls = [
            f"{k} (default)" if k == o["default_blacklist"] else k
            for k in o["blacklists"]
        ]
        print("  blacklists: " + ", ".join(bls + ["none"]) + "\n")


def check_exists(paths):
    missing = [p for p in paths if p and not os.path.exists(p)]
    if missing:
        die(
            "these resources are missing on disk (tell the manifest maintainer):\n  "
            + "\n  ".join(missing)
        )


def bt2(prefix):
    for ext in (".bt2", ".bt2l"):
        if os.path.exists(prefix + ".1" + ext):
            return prefix + ".1" + ext
    return prefix + ".1.bt2"


def q(v):
    """YAML scalar."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return str(v)
    return json.dumps(str(v))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--dir", help="project directory (created if needed)")
    ap.add_argument(
        "--genome", help="genome key or alias from the manifest (hg38, mm10, mm39, ...)"
    )
    ap.add_argument(
        "--blacklist", help="blacklist key from the manifest, none, or a BED path"
    )
    ap.add_argument(
        "--fastq-dir", help="directory holding the FASTQs (default: <project>/fastq)"
    )
    ap.add_argument(
        "--manifest", help=f"resource manifest (default {DEFAULT_MANIFEST})"
    )
    ap.add_argument(
        "--allow-unverified",
        action="store_true",
        help="also offer entries not marked status: ok",
    )
    ap.add_argument(
        "--force", action="store_true", help="overwrite an existing project.yaml"
    )
    ap.add_argument(
        "--list", action="store_true", help="print the manifest's choices and exit"
    )
    ap.add_argument(
        "--json", action="store_true", help="with --list: machine-readable output"
    )
    a = ap.parse_args()

    mpath = find_manifest(a.manifest)
    with open(mpath) as fh:
        man = yaml.safe_load(fh)
    genomes = {
        k: g
        for k, g in man["genomes"].items()
        if g.get("bowtie2_index") and usable(g, a.allow_unverified)
    }
    if a.list:
        list_choices(man, genomes, a.allow_unverified, a.json)
        return

    interactive = sys.stdin.isatty()
    needed = [
        f
        for f, v in (
            ("--dir", a.dir),
            ("--genome", a.genome),
            ("--blacklist", a.blacklist),
        )
        if v is None
    ]
    if needed and not interactive:
        die(
            "not a terminal, so every answer must be a flag. Missing: "
            + ", ".join(needed)
            + ".\nAsk the user (genome; blacklist) -- do not guess. "
            "`pixi run init --list` shows the options.",
            2,
        )

    if a.dir is None:
        a.dir = ask_text("Project directory (created if needed)", ".")
    proj = os.path.abspath(os.path.join(CWD, a.dir))
    out = os.path.join(proj, "project.yaml")
    if os.path.exists(out) and not a.force:
        die(f"{out} exists; use --force to overwrite")

    if a.genome is None:
        a.genome = ask(
            "Which genome are the samples from?",
            [(k, g.get("description", "")) for k, g in genomes.items()],
            "hg38" if "hg38" in genomes else None,
        )
    gkey = resolve_genome(genomes, a.genome)
    g = genomes[gkey]
    bl = choices(g, a.allow_unverified)

    if a.blacklist is None:
        a.blacklist = ask(
            "Blacklist to filter reads, peaks and footprints? (or type a BED path)",
            [(k, v.get("description", "")) for k, v in bl.items()]
            + [("none", "no blacklist filtering")],
            default_blacklist(bl),
        )
    if a.blacklist in bl:
        blacklist = bl[a.blacklist]["path"]
    elif a.blacklist == "none":
        blacklist = ""
    elif os.path.exists(os.path.join(CWD, a.blacklist)):
        blacklist = os.path.abspath(os.path.join(CWD, a.blacklist))
    else:
        die(
            f"blacklist '{a.blacklist}' is neither a {gkey} manifest key "
            f"({', '.join(bl) or 'none'}), 'none', nor an existing file"
        )

    if a.fastq_dir is None:
        a.fastq_dir = (
            ask_text(
                "Directory with the FASTQs (relative to the project directory, "
                "or absolute)",
                "fastq",
            )
            if interactive
            else "fastq"
        )

    gtf = g.get("gtf") or ""
    check_exists(
        [g["fasta"], g["fasta"] + ".fai", bt2(g["bowtie2_index"]), gtf, blacklist]
    )

    run = (
        f"pixi run snakemake -s workflow/Snakefile --directory {proj} "
        f"--configfile {out} --profile profiles/local -n"
    )
    lines = [
        "# =============================================================================",
        f"# ATAC-seq project config -- written by `pixi run init` on "
        f"{datetime.date.today()}",
        f"# from {mpath} (manifest version {man['manifest_version']}, "
        f"updated {man.get('updated')})",
        f"# Choices: genome={gkey}  blacklist={a.blacklist}",
        "#",
        "# Holds only the keys that differ from the pipeline's config/config.yaml.",
        "# Relative paths are relative to this directory. Run from the pipeline repo:",
        f"#   {run}",
        "# =============================================================================",
        "",
        "# Fill these in (see the pipeline README, 'Samplesheet').",
        "samples: samples.tsv",
        "contrasts: contrasts.tsv",
        f"fastq_dir: {q(a.fastq_dir)}",
        "",
        f"# {g.get('description', gkey)}",
        "reference:",
        f"  fasta: {q(g['fasta'])}",
        f"  bowtie2_index: {q(g['bowtie2_index'])}",
        f"  gtf: {q(gtf)}",
        f"  blacklist: {q(blacklist)}"
        + (f"    # {a.blacklist}" if a.blacklist in bl else ""),
        f"  effective_genome_size: {q(g['effective_genome_size'])}",
        f"  macs2_gsize: {q(g['macs2_gsize'])}",
        f"  mito_chrom: {q(g.get('mito_chrom', 'chrM'))}",
        "",
    ]
    if not gtf:
        lines += [
            "# No GTF for this genome: the qc module (TSS profile, tracks) needs one.",
            "modules:",
            "  qc: false",
            "",
        ]

    os.makedirs(proj, exist_ok=True)
    with open(out, "w") as fh:
        fh.write("\n".join(lines))
    copied = []
    for name in ("samples.tsv", "contrasts.tsv"):
        dst = os.path.join(proj, name)
        if not os.path.exists(dst):
            shutil.copy(os.path.join(REPO, "config", name), dst)
            copied.append(dst)

    print(f"\nWrote {out}")
    for c in copied:
        what = "libraries" if c.endswith("samples.tsv") else "contrasts"
        print(
            f"Copied the example {os.path.basename(c)} to {c} -- replace its rows "
            f"with your {what}."
        )
    print(f"\nNext, a dry run from {REPO}:\n  {run}")


if __name__ == "__main__":
    main()
