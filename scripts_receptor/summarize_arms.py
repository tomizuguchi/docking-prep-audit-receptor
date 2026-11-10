#!/usr/bin/env python3
"""Net charge per structure across the three derived inputs, and the split by metal and ligand.

  python3 summarize_arms.py <pqr_hetall> <pqr_nowat> <pqr_recept> <structure_dir> <arms_diff_3way.csv> [out.csv]

compare_arms.py reports which residues change state when the HETATM records are removed. This script
reports what follows from that per structure:

  - the sum of the atomic charges of each pqr file, and how far it moves between the input that keeps
    the non-water HETATM records (nowat) and the input that keeps none (recept)
  - which structures contain a non-water HETATM record (a "ligand") and which contain a metal
  - the number of structures and residues affected in each of those groups
  - for the changed residues in structures with a metal, the distance from the titratable atoms to the
    nearest metal, the element of that metal, and the direction of the change

The differences in charge are rounded to the nearest integer before they are counted. The sums carry
floating-point noise (0.9999999999999991 for a difference of one), and a comparison such as ">= 1.0"
made without rounding puts such a structure in the wrong class.

<structure_dir> is the directory of PDB files the distances are measured in. It must contain the
HETATM records, so pass full100 or in_nowat, not recept. The distance definition is the one in
check_receptor.py and measure_proximity.py: the minimum over all titratable atoms of the residue.
"""

import collections
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_receptor import TITR_ATOM, dist, read_pdb

SPHERE = 2.6  # A. A titratable atom at or under this distance from a metal is counted as coordinating.

FIELDS = [
    "pdb_id",
    "has_ligand",
    "has_metal",
    "metal_elements",
    "sumq_hetall",
    "sumq_nowat",
    "sumq_recept",
    "d_sumq_ligand",
    "n_residues_changed",
]


def sumq(pqr):
    """Sum of the atomic charges of a pqr file. Charge and radius are the last two fields of a line.

    The coordinates of a pqr line can run together, so the fields are counted from the right."""
    total = 0.0
    for l in open(pqr, errors="replace"):
        if l.startswith(("ATOM", "HETATM")):
            total += float(l.split()[-2])
    return total


def ids(d):
    return {os.path.basename(f)[:4].upper() for f in glob.glob(os.path.join(d, "*.pqr"))}


def main():
    if len(sys.argv) < 6 or sys.argv[1] in ("-h", "--help"):
        sys.exit(__doc__)
    d_het, d_now, d_rec, sdir, diff_csv = sys.argv[1:6]
    out = sys.argv[6] if len(sys.argv) > 6 else None

    pids = sorted(ids(d_het) & ids(d_now) & ids(d_rec))
    changed = collections.defaultdict(list)
    for r in csv.DictReader(open(diff_csv)):
        if r["ligand_removal_changed"] == "1":
            changed[r["pdb_id"]].append(r)

    rows, struct = [], {}
    for pid in pids:
        res, het, metals, _nh, _g, _inc = read_pdb(os.path.join(sdir, pid + ".pdb"))
        struct[pid] = (res, metals)
        s = [sumq(os.path.join(d, pid + ".pqr")) for d in (d_het, d_now, d_rec)]
        rows.append(
            {
                "pdb_id": pid,
                "has_ligand": int(bool(het)),
                "has_metal": int(bool(metals)),
                "metal_elements": "+".join(sorted({el for _rn, el, _xyz in metals})),
                "sumq_hetall": round(s[0], 3),
                "sumq_nowat": round(s[1], 3),
                "sumq_recept": round(s[2], 3),
                "d_sumq_ligand": round(s[2] - s[1]),
                "n_residues_changed": len(changed.get(pid, ())),
            }
        )
    if out:
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS)
            w.writeheader()
            w.writerows(rows)

    n = len(rows)
    print(f"{n} structures with a pqr file under all three inputs")
    water = sum(1 for r in rows if round(r["sumq_hetall"] - r["sumq_nowat"]) != 0)
    print(f"  net charge differs between hetall and nowat (water removed): {water}")
    moved = [r for r in rows if r["d_sumq_ligand"] != 0]
    print(f"  net charge differs between nowat and recept (ligand removed): {len(moved)}/{n}")
    hist = collections.Counter(abs(r["d_sumq_ligand"]) for r in moved)
    print("  size of the change: " + ", ".join(f"{k} unit(s) in {hist[k]}" for k in sorted(hist)))
    print(f"  one or more: {len(moved)}   two or more: {sum(v for k, v in hist.items() if k >= 2)}")

    groups = (
        ("with a metal", lambda r: r["has_metal"]),
        ("with a ligand but no metal", lambda r: r["has_ligand"] and not r["has_metal"]),
        ("with no ligand", lambda r: not r["has_ligand"]),
        ("without a metal (the two groups above)", lambda r: not r["has_metal"]),
    )
    print(
        f"\n  {'group':<42}{'structures':>10}{'affected':>10}{'residues':>10}{'per struct':>12}"
        f"{'charge moved':>14}{'sum|d|':>8}{'max|d|':>8}"
    )
    for label, test in groups:
        g = [r for r in rows if test(r)]
        aff = [r for r in g if r["n_residues_changed"]]
        nres = sum(r["n_residues_changed"] for r in g)
        mv = [abs(r["d_sumq_ligand"]) for r in g if r["d_sumq_ligand"]]
        print(
            f"  {label:<42}{len(g):>10}{len(aff):>10}{nres:>10}{nres / len(g) if g else 0:>12.2f}"
            f"{len(mv):>14}{sum(mv):>8}{max(mv, default=0):>8}"
        )

    # Metal content of the whole structure directory, not only the structures compared above.
    allfiles = sorted(glob.glob(os.path.join(sdir, "*.pdb")))
    el_count, n_lig, n_met = collections.Counter(), 0, 0
    for f in allfiles:
        _res, het, metals, _nh, _g, _inc = read_pdb(f)
        n_lig += bool(het)
        n_met += bool(metals)
        el_count.update({el for _rn, el, _xyz in metals})
    print(f"\n{len(allfiles)} structures in {sdir}: {n_lig} with a ligand, {n_met} with a metal")
    print("  structures per metal element (a structure counts once for each element it contains):")
    print("  " + ", ".join(f"{el} {c}" for el, c in sorted(el_count.items(), key=lambda kv: (-kv[1], kv[0]))))

    # Changed residues in structures with a metal: distance to the nearest metal, and what changed.
    bins = collections.Counter()
    by_el = collections.defaultdict(collections.Counter)
    n_in_metal = 0
    for pid, rs in sorted(changed.items()):
        if pid not in struct:
            continue
        res, metals = struct[pid]
        if not metals:
            continue
        for r in rs:
            n_in_metal += 1
            hits = [e for (_ch, num, rn), e in res.items() if num == str(r["resnum"]) and rn == r["base"]]
            best = None
            for e in hits:
                for a in TITR_ATOM.get(r["base"], ()):
                    if a not in e["atoms"]:
                        continue
                    for _rn, el, xyz in metals:
                        d = dist(e["atoms"][a], xyz)
                        if best is None or d < best[0]:
                            best = (d, el)
            if best is None:
                bins["no titratable atom found"] += 1
                continue
            d, el = best
            if d <= SPHERE:
                bins[f"within {SPHERE} A"] += 1
                by_el[el][f"{r['base']} {r['state_nowat']} -> {r['state_none']}"] += 1
            elif d <= 4.0:
                bins[f"{SPHERE} to 4 A"] += 1
            elif d <= 6.0:
                bins["4 to 6 A"] += 1
            else:
                bins["beyond 6 A"] += 1
    print(f"\n{n_in_metal} changed residues in structures with a metal, by distance to the nearest metal:")
    for k in (f"within {SPHERE} A", f"{SPHERE} to 4 A", "4 to 6 A", "beyond 6 A", "no titratable atom found"):
        if bins[k]:
            print(f"  {k:<26}{bins[k]:>4}")
    print(f"  within {SPHERE} A, by the element of the nearest metal (state with the ligand -> state without):")
    for el, c in sorted(by_el.items(), key=lambda kv: (-sum(kv[1].values()), kv[0])):
        print(f"    {el:<3}{sum(c.values()):>3}   " + "; ".join(f"{k} x{v}" for k, v in c.most_common()))
    if out:
        print(f"\n-> {out}")


if __name__ == "__main__":
    main()
