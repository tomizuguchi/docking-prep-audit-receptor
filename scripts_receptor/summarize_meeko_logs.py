#!/usr/bin/env python3
"""Tabulate what the logs of the Meeko runs say, one row per structure and run.

  python3 summarize_meeko_logs.py <meeko_dir> [<meeko_dir> ...] [--prof prof100.csv] [--out meeko_runs.csv]

run_meeko.py writes <ID>.log (standard output and standard error of mk_prepare_receptor.py together)
and, on success, <ID>.pdbqt. The numbers in the manuscript that come from those logs were first read
off by hand. This script reads them, so that they can be regenerated:

  - whether the run completed (a pdbqt with at least one ATOM record)
  - the residues listed as "Ignored due to allow_bad_res", the line Meeko prints when the flag deletes
  - the residues listed as "Template matching failed for" when the run stops instead
  - the residues reported with "matched with excess inter-residue bond(s)"
  - the input residues reported as "not in residue_templates", and whether a download of a residue
    definition from rcsb.org was attempted, and for which residue it failed
  - whether the run stopped on alternate locations, or on a residue that could not be built
  - whether the log mentions alternate locations at all

Meeko's own recommendation for batch processing is part of the error message of a failed run. The
number of logs that carry it is counted, because the manuscript cites it as the source of the
option pair "--default_altloc A -a".

With --prof, the number of excess-bond residues is compared with incomplete_sc in prof100.csv.

The logs contain the paths of the machine they were written on and are not deposited; the table this
script writes is.
"""

import argparse
import ast
import collections
import csv
import glob
import os
import re

FIELDS = [
    "run",
    "pdb_id",
    "completed",
    "n_ignored",
    "n_template_failed",
    "n_excess_bond",
    "n_no_template",
    "no_template_codes",
    "download_attempted",
    "download_failed_for",
    "stopped_on_altloc",
    "build_failed_for",
    "mentions_altloc",
    "batch_recommendation",
]


def has_atoms(path):
    return os.path.exists(path) and any(l.startswith("ATOM") for l in open(path, errors="replace"))


def keys_in(text):
    return re.findall(r"'([^']+)'", text)


def read_log(path):
    t = open(path, errors="replace").read()
    rec = dict.fromkeys(FIELDS, 0)
    m = re.search(r"Template matching failed for: (\[.*?\])\s*Ignored due to allow_bad_res", t, re.S)
    if m:
        rec["n_ignored"] = len(keys_in(m.group(1)))
    else:
        m = re.search(r"Template matching failed for: (\[.*?\])", t, re.S)
        if m:
            rec["n_template_failed"] = len(keys_in(m.group(1)))
    rec["n_excess_bond"] = len(set(re.findall(r"matched with excess inter-residue bond\(s\): (\S+)", t)))
    m = re.search(r"Input residues (\{.*?\}) not in residue_templates", t, re.S)
    codes = sorted(set(ast.literal_eval(m.group(1)).values())) if m else []
    rec["n_no_template"] = len(ast.literal_eval(m.group(1))) if m else 0
    rec["no_template_codes"] = "+".join(codes)
    rec["download_attempted"] = int(bool(re.search(r"rcsb\.org|urlopen|download CIF", t)))
    m = re.search(r"Could not download CIF file for (\w+)", t)
    rec["download_failed_for"] = m.group(1) if m else ""
    rec["stopped_on_altloc"] = int(bool(re.search(r"Requested altlocs not found|Either specify an altloc", t)))
    m = re.search(r"unable to build rdkit mol for residue (\w+)", t)
    rec["build_failed_for"] = m.group(1) if m else ""
    rec["mentions_altloc"] = int(bool(re.search(r"altloc|alternate", t, re.I)))
    rec["batch_recommendation"] = int("(for batch processing)" in t)
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("meeko_dirs", nargs="+")
    ap.add_argument("--prof", help="prof100.csv, to compare excess-bond counts with incomplete_sc")
    ap.add_argument("--out", help="CSV to write one row per structure and run to")
    a = ap.parse_args()

    prof = {}
    if a.prof:
        prof = {r["id"].upper(): int(r["incomplete_sc"] or 0) for r in csv.DictReader(open(a.prof))}

    rows = []
    for d in a.meeko_dirs:
        run = os.path.basename(os.path.normpath(d))
        logs = sorted(glob.glob(os.path.join(d, "*.log")))
        mine = []
        for f in logs:
            pid = os.path.basename(f)[:4].upper()
            rec = read_log(f)
            rec["run"], rec["pdb_id"] = run, pid
            rec["completed"] = int(has_atoms(os.path.join(d, pid + ".pdbqt")))
            mine.append(rec)
        rows += mine

        n = len(mine)
        done = [r for r in mine if r["completed"]]
        fail = [r for r in mine if not r["completed"]]
        print(f"{run}: {n} logs, {len(done)} completed, {len(fail)} without output")
        ign = [r for r in done if r["n_ignored"]]
        if ign:
            print(f"  completed with an 'Ignored due to allow_bad_res' line: {len(ign)} structures,"
                  f" {sum(r['n_ignored'] for r in ign)} residues listed")
        print(f"  completed runs whose log mentions alternate locations: {sum(r['mentions_altloc'] for r in done)}")
        if fail:
            kinds = collections.Counter()
            for r in fail:
                if r["download_attempted"]:
                    kinds["download of a residue definition attempted"] += 1
                elif r["stopped_on_altloc"]:
                    kinds["alternate locations"] += 1
                elif r["build_failed_for"]:
                    kinds[f"residue could not be built ({r['build_failed_for']})"] += 1
                elif r["n_excess_bond"]:
                    kinds["template matching failed, excess inter-residue bond reported"] += 1
                elif r["n_template_failed"]:
                    kinds["template matching failed, no excess bond reported"] += 1
                else:
                    kinds["other"] += 1
            for k, v in kinds.most_common():
                print(f"  without output, {k}: {v}")
            print(f"  failed runs whose message recommends -a and --default_altloc for batch processing:"
                  f" {sum(r['batch_recommendation'] for r in fail)}")
        dl = [r for r in mine if r["download_attempted"]]
        print(f"  runs that attempted a download: {len(dl)}")
        if dl:
            codes = collections.Counter(c for r in dl for c in r["no_template_codes"].split("+") if c)
            per = collections.Counter(r["n_no_template"] for r in dl)
            print(f"    distinct residue codes without a template: {len(codes)};"
                  f" most frequent (structures): " + ", ".join(f"{k} {v}" for k, v in codes.most_common(6)))
            print("    residues without a template per structure: "
                  + ", ".join(f"{k} in {per[k]}" for k in sorted(per)))
            first = collections.Counter(r["download_failed_for"] for r in dl)
            print(f"    the download that failed, by residue code: {len(first)} distinct codes"
                  " (which one is tried first is not fixed between runs)")
        if prof:
            eb = [r for r in mine if r["n_excess_bond"] and not r["completed"]]
            if eb:
                d_ = [r["n_excess_bond"] - prof.get(r["pdb_id"], 0) for r in eb]
                print(f"  {len(eb)} failed runs report excess inter-residue bonds; against incomplete_sc:")
                print(f"    equal {sum(1 for x in d_ if x == 0)}, fewer reported {sum(1 for x in d_ if x < 0)}"
                      f" ({dict(sorted(collections.Counter(x for x in d_ if x < 0).items()))}),"
                      f" more reported {sum(1 for x in d_ if x > 0)}"
                      f" ({dict(sorted(collections.Counter(x for x in d_ if x > 0).items()))})")
                more = [(r["pdb_id"], x) for r, x in zip(eb, d_, strict=True) if x > 0]
                print("    more reported in: " + ", ".join(f"{p} +{x}" for p, x in more))
            inc = {p for p, v in prof.items() if v}
            if inc and fail:
                fset = {r["pdb_id"] for r in fail}
                print(f"  structures with an incomplete side chain: {len(inc)}, of which failed {len(inc & fset)};"
                      f" with none: {n - len(inc)}, of which completed {n - len(inc) - len(fset - inc)}")
        print()

    if a.out:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS)
            w.writeheader()
            w.writerows(rows)
        print(f"-> {a.out}")


if __name__ == "__main__":
    main()
