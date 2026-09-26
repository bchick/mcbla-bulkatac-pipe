"""Regression checks on the outputs of `pixi run test-all`.

Run from the repo root after the full test workflow:

    python .test/scripts/check_results.py [.test/results]

Checks:
  * normcheck/csaw: the csaw contrasts are really re-tested under csaw size
    factors. DiffBind's legacy per-contrast mode used to ignore
    dba.normalize() in the DESeq2 test, so csaw and depth gave identical
    p-values while only Fold moved. When the reported folds differ, the
    p-values must differ too.
"""

import csv
import os
import sys

RESULTS = sys.argv[1] if len(sys.argv) > 1 else ".test/results"
failures = []


def read_tsv(path):
    with open(path) as fh:
        return {row["peak_id"]: row for row in csv.DictReader(fh, delimiter="\t")}


def check_csaw():
    depth_dir = os.path.join(RESULTS, "diff/depth/tables")
    csaw_dir = os.path.join(RESULTS, "normcheck/csaw/tables")
    if not os.path.isdir(csaw_dir):
        failures.append(f"missing {csaw_dir} (run `pixi run test-all` first)")
        return
    tables = sorted(f for f in os.listdir(csaw_dir) if f.endswith("_all.tsv"))
    if not tables:
        failures.append(f"no *_all.tsv tables in {csaw_dir}")
    for name in tables:
        depth = read_tsv(os.path.join(depth_dir, name))
        csaw = read_tsv(os.path.join(csaw_dir, name))
        shared = depth.keys() & csaw.keys()
        if not shared:
            failures.append(f"{name}: no peaks shared between depth and csaw")
            continue
        fold_moved = any(
            abs(float(depth[p]["Fold"]) - float(csaw[p]["Fold"])) > 1e-6 for p in shared
        )
        p_moved = any(depth[p]["p.value"] != csaw[p]["p.value"] for p in shared)
        if fold_moved and not p_moved:
            failures.append(
                f"{name}: csaw Fold differs from depth but every p-value is "
                "identical; the csaw size factors are not used by the test"
            )
        else:
            print(f"ok  csaw {name}: fold moved={fold_moved}, p-values moved={p_moved}")


check_csaw()

if failures:
    for f in failures:
        print(f"FAIL {f}", file=sys.stderr)
    sys.exit(1)
print("all checks passed")
