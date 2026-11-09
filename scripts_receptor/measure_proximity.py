#!/usr/bin/env python3
"""Measure how close the residues whose titration state moved are to the metal or the ligand.

  python3 measure_proximity.py <diff.csv> <structure_dir> [--col ligand_removal_changed]

Takes `arms_diff_3way.csv` (residues whose state moved with the order of ligand removal) or `tool_diff.csv`
(residues whose formal charge disagrees between the tools) and, for each residue, gives

  - the distance to the nearest metal
  - the distance to the nearest HETATM other than water

The distances are measured from all of the titratable atoms of the residue.

Why this was written
  This measurement was done earlier with one-off code, and no script was left.
  "48% within 2.6 A of a metal" and "all 64 residues within 9.4 A of the ligand" in the record
  were in a state where they could not be regenerated.

Why it is measured again
  W1 of check_receptor.py looked at only the first of the titratable atoms and stopped there.
  For HIS it used only ND1 and did not look at NE2, and for ARG it used only NE. The decision on
  the coordination sphere changes, so this table has to be measured again with the same definition.
  With --legacy, the earlier "first atom only" behaviour is reported alongside.

Note
  The diff CSV has no chain column. When the same residue number is in several chains, all of the
  chains are measured and the minimum is taken, and how many cases were ambiguous is reported separately.
"""

import argparse
import collections
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_receptor import TITR_ATOM, dist, read_pdb


def min_dist(e, rn, targets, first_only=False):
    """Shortest distance from the titratable atoms of the residue to targets, and the name of that atom.

    first_only=True is the behaviour before the fix (stop at the first atom that exists)."""
    best, who = None, None
    for a in TITR_ATOM.get(rn, ()):
        if a not in e["atoms"]:
            continue
        d = min(dist(e["atoms"][a], t) for t in targets)
        if best is None or d < best:
            best, who = d, a
        if first_only:
            break
    return best, who


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("diff_csv")
    ap.add_argument("structure_dir")
    ap.add_argument("--col", default=None, help="use only the rows where this column is '1' (default is all rows)")
    ap.add_argument("--legacy", action="store_true", help="also report the 'first atom only' values before the fix")
    ap.add_argument("--out", help="CSV to write the per-row results to")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.diff_csv)))
    if a.col:
        rows = [r for r in rows if r.get(a.col) == "1"]
    print(f"{a.diff_csv}: {len(rows)} rows selected{f' ({a.col}=1 only)' if a.col else ''}")

    cache, out, ambiguous, missing = {}, [], 0, []
    for r in rows:
        pid, num, base = r["pdb_id"], r["resnum"], r["base"]
        if pid not in cache:
            p = os.path.join(a.structure_dir, f"{pid}.pdb")
            cache[pid] = read_pdb(p) if os.path.isfile(p) else None
        got = cache[pid]
        if got is None:
            missing.append(pid)
            continue
        res, het, metals, _nH, _g, _inc = got
        # The chain is not in the CSV, so pick up the same residue number and residue name from all chains
        hits = [e for (ch, n, rn), e in res.items() if n == str(num) and rn == base]
        if not hits:
            missing.append(f"{pid}:{base}{num}")
            continue
        if len(hits) > 1:
            ambiguous += 1
        mx = [xyz for _, _, xyz in metals]
        hx = [xyz for _, _, xyz in het]
        rec = dict(
            pdb_id=pid, resnum=num, base=base, n_chains=len(hits), has_metal=int(bool(mx)), has_het=int(bool(hx))
        )
        for tag, tg in (("metal", mx), ("het", hx)):
            if not tg:
                rec[f"d_{tag}"] = ""
                rec[f"atom_{tag}"] = ""
                if a.legacy:
                    rec[f"d_{tag}_legacy"] = ""
                continue
            cands = [min_dist(e, base, tg) for e in hits]
            cands = [c for c in cands if c[0] is not None]
            if not cands:
                rec[f"d_{tag}"] = ""
                rec[f"atom_{tag}"] = ""
                if a.legacy:
                    rec[f"d_{tag}_legacy"] = ""
                continue
            d, who = min(cands, key=lambda c: c[0])
            rec[f"d_{tag}"] = round(d, 2)
            rec[f"atom_{tag}"] = who
            if a.legacy:
                lc = [min_dist(e, base, tg, first_only=True) for e in hits]
                lc = [c for c in lc if c[0] is not None]
                rec[f"d_{tag}_legacy"] = round(min(c[0] for c in lc), 2) if lc else ""
        out.append(rec)

    if missing:
        print(f"  !! {len(missing)} not found: {missing[:6]}")
    if ambiguous:
        print(f"  !! {ambiguous} rows whose residue number is in several chains (the minimum over all chains taken)")

    def summarize(key, label, cut):
        vals = [r[key] for r in out if r[key] != ""]
        if not vals:
            print(f"\n{label}: none")
            return
        vals = sorted(vals)
        n = len(vals)
        within = sum(1 for v in vals if v < cut)
        print(f"\n{label}  {n} residues")
        print(f"  median {vals[n // 2]:.2f} A   min {vals[0]:.2f}   max {vals[-1]:.2f}")
        print(f"  under {cut} A {within} cases ({100 * within / n:.0f}%)")
        for t in (4.0, 6.0, 10.0):
            print(
                f"  under {t:>4.0f} A {sum(1 for v in vals if v < t):>3} cases"
                f" ({100 * sum(1 for v in vals if v < t) / n:>3.0f}%)"
            )

    summarize("d_metal", "distance to the metal", 2.6)
    summarize("d_het", "distance to the HETATM other than water", 2.6)

    if a.legacy:
        for tag, cut in (("metal", 2.6),):
            new = [r for r in out if r[f"d_{tag}"] != "" and r[f"d_{tag}"] < cut]
            old = [r for r in out if r.get(f"d_{tag}_legacy") != "" and r[f"d_{tag}_legacy"] < cut]
            print(f"\nnumber of cases under {cut} A for {tag}")
            print(f"  before the fix (first atom only) {len(old)}")
            print(f"  after the fix (all atoms)        {len(new)}")
            gained = [r for r in new if r not in old]
            if gained:
                c = collections.Counter((r["base"], r[f"atom_{tag}"]) for r in gained)
                print(
                    "  residues newly inside the range: "
                    + ", ".join(f"{at} of {b} {n} cases" for (b, at), n in c.most_common())
                )

    print("\nBreakdown of the atoms (the titratable atom nearest the metal):")
    c = collections.Counter(r["atom_metal"] for r in out if r["atom_metal"])
    for (at), n in c.most_common():
        print(f"  {at:<5} {n}")

    if a.out and out:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0]))
            w.writeheader()
            w.writerows(out)
        print(f"\n-> {a.out}")


if __name__ == "__main__":
    main()
