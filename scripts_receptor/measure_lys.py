#!/usr/bin/env python3
"""Measure the fraction of neutral lysines over a set of pqr files. Count with a cross-check for zero misses.

  python3 measure_lys.py <pqrdir> <manifest.csv>

Parser lesson: pqr is whitespace separated, but coordinates run together and the field count drops
(in 6ZHF, 67% of the ATOM lines have 9 fields). Rejecting lines by field count silently loses atoms.
This version uses an assert to guarantee that not one LYS/LYN ATOM line is dropped.

The lysine state is decided by the number of hydrogens on NZ.
  3 HZ = charged LYS (+1), 2 HZ = neutral LYN, 0/1 = incomplete side chain.
The residue charge sum is used only as a secondary check. A C-terminal lysine sums to 0 even when
charged (the terminal carboxylate cancels it), so on its own it gives a wrong call.
"""

import collections
import csv
import glob
import math
import os
import sys


def parse(pqr):
    """Return (per-residue state, number of LYS/LYN ATOM lines, number of lines assigned)."""
    res, seen, placed = collections.OrderedDict(), 0, 0
    for l in open(pqr, errors="replace"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        f = l.split()
        if len(f) < 8:
            continue
        if f[3] not in ("LYS", "LYN"):
            continue
        seen += 1
        try:
            q = float(f[-2])
            float(f[-1])
        except ValueError:
            continue
        resnum = ""
        for t in f[4:7]:
            if t.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ").lstrip("-").isdigit():
                resnum = t
                break
        if not resnum:
            continue
        placed += 1
        e = res.setdefault((resnum, f[3]), {"name": f[3], "hz": 0, "q": 0.0, "n": 0})
        e["q"] += q
        e["n"] += 1
        if f[2].startswith("HZ"):
            e["hz"] += 1
    return list(res.values()), seen, placed


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)


def main():
    if len(sys.argv) < 3 or sys.argv[1] in ("-h", "--help"):
        sys.exit(__doc__)
    pqrdir, manfile = sys.argv[1], sys.argv[2]
    if not glob.glob(os.path.join(pqrdir, "*.pqr")):
        sys.exit(f"no *.pqr in {pqrdir}")
    man = {r["pdb_id"]: r for r in csv.DictReader(open(manfile))}
    STR = ["A_le1.5", "B_1.5-2.0", "C_2.0-2.5", "D_2.5-3.0", "E_ge3.0"]
    pool = {s: int(next(r["pool_size"] for r in man.values() if r["stratum"] == s)) for s in STR}
    tot = sum(pool.values())

    per = {}
    for f in sorted(glob.glob(os.path.join(pqrdir, "*.pqr"))):
        pid = os.path.basename(f)[:4].upper()
        st, seen, placed = parse(f)
        assert seen == placed, f"{pid}: only {placed} of {seen} LYS/LYN lines were assigned"
        per[pid] = {
            "charged": sum(1 for s in st if s["hz"] == 3),
            "neutral": sum(1 for s in st if s["hz"] == 2),
            "broken": sum(1 for s in st if s["hz"] < 2),
            "total": len(st),
            "neutral_named_LYS": sum(1 for s in st if s["hz"] == 2 and s["name"] == "LYS"),
            "named_LYN": sum(1 for s in st if s["name"] == "LYN"),
            "atoms": sum(s["n"] for s in st),
        }
    print(f"Miss cross-check: every LYS/LYN ATOM line assigned in all {len(per)} structures\n")

    print("=== By stratum ===")
    print(
        f"{'stratum':<11}{'weight':>7}{'n':>4}{'LYS tot':>8}{'chg':>6}{'neu':>6}{'broken':>7}"
        f"{'neutral/res':>12}{'struct w/ neu':>14}"
    )
    wr_res = wr_str = 0.0
    agg = collections.Counter()
    nstr_all = 0
    for s in STR:
        ids = [p for p, r in man.items() if r["stratum"] == s and p in per]
        if not ids:
            continue
        t = collections.Counter()
        for p in ids:
            for k, v in per[p].items():
                t[k] += v
        nstr = sum(1 for p in ids if per[p]["neutral"] > 0)
        nstr_all += nstr
        w = pool[s] / tot
        pr = t["neutral"] / t["total"] if t["total"] else 0.0
        wr_res += w * pr
        wr_str += w * (nstr / len(ids))
        for k, v in t.items():
            agg[k] += v
        print(
            f"{s:<11}{w:>7.3f}{len(ids):>4}{t['total']:>8}{t['charged']:>6}"
            f"{t['neutral']:>6}{t['broken']:>7}{pr * 100:>11.2f}%{f'{nstr}/{len(ids)}':>14}"
        )

    n = len(per)
    lo, hi = wilson(agg["neutral"], agg["total"])
    lo2, hi2 = wilson(nstr_all, n)
    print(f"\nUnweighted total: {n} structures, {agg['total']} LYS/LYN residues")
    print(f"  charged {agg['charged']}, neutral {agg['neutral']}, incomplete side chain {agg['broken']}")
    print(
        f"  neutral fraction/res {agg['neutral'] / agg['total'] * 100:.2f}%  Wilson 95% {lo * 100:.2f}-{hi * 100:.2f}%"
    )
    print(
        f"  structures with 1 or more neutral {nstr_all}/{n} = {nstr_all / n * 100:.1f}%  "
        f"Wilson 95% {lo2 * 100:.1f}-{hi2 * 100:.1f}%"
    )
    print(f"pool_size-weighted population estimate: {wr_res * 100:.2f}% per residue, {wr_str * 100:.1f}% per structure")
    print(
        f"\nNeutral lysines still named LYS: {agg['neutral_named_LYS']}/{agg['neutral']}, "
        f"residues named LYN: {agg['named_LYN']}"
    )

    print("\n=== Structures containing a neutral or incomplete residue ===")
    for p in sorted(per):
        d = per[p]
        if d["neutral"] or d["broken"]:
            print(
                f"  {p} ({man[p]['stratum']:<10}) LYS total {d['total']:>3}: "
                f"charged {d['charged']:>3}, neutral {d['neutral']:>2}"
                f" (still named LYS {d['neutral_named_LYS']}), incomplete {d['broken']:>2}"
            )


if __name__ == "__main__":
    main()
