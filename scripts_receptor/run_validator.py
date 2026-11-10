#!/usr/bin/env python3
"""Run check_receptor.py over every structure under the option sets the manuscript reports.

  python3 run_validator.py <pdb_dir> [--out validator_runs.csv]

check_receptor.py checks one file. The counts in the manuscript (which rule fires on how many
structures, and for how many structures the exit status is non-zero) come from running it over the
whole sample, and the option set decides them:

  default            W5 (incomplete side chain) is an ERROR, both altloc A and B are examined
  --altloc A         as the default, with only altloc A examined
  --allow-bad-res    W5 drops to a warning; this is the setting for a Meeko run that passes -a

Under --allow-bad-res the errors that remain are E1 and E2, the two rules for runs that stop even
with -a. Under the default the exit status is also non-zero for every structure with an incomplete
side chain, because without -a Meeko stops on those.

One row is written per structure, with the exit status and the rule identifiers raised under each
option set. E2 is written with the altloc it refers to, as E2[A] or E2[B], when the validator's
message names it.
"""

import argparse
import collections
import csv
import glob
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SETS = (("default", []), ("altlocA", ["--altloc", "A"]), ("allow_bad_res", ["--allow-bad-res"]))


def run(path, extra):
    p = subprocess.run(
        [sys.executable, os.path.join(HERE, "check_receptor.py"), path, "--json", *extra],
        capture_output=True,
        text=True,
    )
    ids = []
    try:
        for f in json.loads(p.stdout)["findings"]:
            fid = f["id"]
            if fid == "E2":
                m = re.search(r"altloc\s+([A-Z])\b", f["message"])
                if m:
                    fid = f"E2[{m.group(1)}]"
            ids.append((fid, f["level"]))
    except (ValueError, KeyError):
        ids.append(("unparsed", "ERROR"))
    return p.returncode, ids


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdb_dir")
    ap.add_argument("--out", help="CSV to write one row per structure to")
    a = ap.parse_args()

    files = sorted(glob.glob(os.path.join(a.pdb_dir, "*.pdb")))
    if not files:
        sys.exit(f"no *.pdb in {a.pdb_dir}")
    fields = ["pdb_id"] + [f"{k}_{s}" for s, _ in SETS for k in ("exit", "errors", "warnings")]
    rows = []
    fired = {s: collections.Counter() for s, _ in SETS}
    nonzero = {s: [] for s, _ in SETS}
    for f in files:
        pid = os.path.basename(f)[:4].upper()
        row = {"pdb_id": pid}
        for s, extra in SETS:
            rc, ids = run(f, extra)
            row[f"exit_{s}"] = rc
            row[f"errors_{s}"] = "+".join(sorted({i for i, lv in ids if lv == "ERROR"}))
            row[f"warnings_{s}"] = "+".join(sorted({i for i, lv in ids if lv == "WARN"}))
            fired[s].update({re.sub(r"\[.\]", "", i) for i, _lv in ids})
            if rc:
                nonzero[s].append(pid)
        rows.append(row)
    if a.out:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)

    print(f"{len(files)} structures")
    for s, extra in SETS:
        print(f"\noption set: {s} ({' '.join(extra) or 'no options'})")
        print(f"  non-zero exit status: {len(nonzero[s])} structures")
        if len(nonzero[s]) <= 10:
            print(f"    {', '.join(nonzero[s])}")
        print("  structures per rule: " + ", ".join(f"{k} {v}" for k, v in sorted(fired[s].items())))
    if a.out:
        print(f"\n-> {a.out}")


if __name__ == "__main__":
    main()
