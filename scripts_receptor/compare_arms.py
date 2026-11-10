#!/usr/bin/env python3
"""Compare residue by residue how the titration state changes with and without HETATM.

  python3 compare_arms.py <pqr_hetall> <pqr_nowat> <pqr_none> [out.csv]

The 3 conditions have exactly the same protein ATOM records and differ only in the presence of HETATM.
Hydrogens are dropped in all 3 conditions, so deposited hydrogens have no influence.

  hetall : heavy atoms + all HETATM (ligand + water)
  nowat  : heavy atoms + HETATM other than water (water removed only)
  none   : heavy atoms only (both water and ligand removed)

The hetall -> nowat difference is the change from removing water, nowat -> none from removing ligands.
States are classified by classify_titr.py (by the presence of hydrogens, not by residue name).
"""

import collections
import csv
import glob
import importlib.util
import os
import sys

spec = importlib.util.spec_from_file_location(
    "ct", os.path.join(os.path.dirname(os.path.abspath(__file__)), "classify_titr.py")
)
ct = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ct)


def states(pqr):
    res, seen, placed = ct.parse(pqr)
    assert seen == placed, f"{pqr}: only {placed} of {seen} titratable residue lines were assigned"
    # A residue with the same atom name twice means same-numbered residues of different chains were
    # merged, or altLoc. Passing it silently creates a fictitious HIP(+1) from one HD1 and another HE2.
    dup = sum(1 for v in res.values() if v.get("dup"))
    return {k: v for k, v in ct.classify(res).items()}, dup


def main():
    if len(sys.argv) < 4 or sys.argv[1] in ("-h", "--help"):
        sys.exit(__doc__)
    dirs = sys.argv[1:4]
    for d in dirs:
        if not glob.glob(os.path.join(d, "*.pqr")):
            sys.exit(f"no *.pqr in {d}")
    out = sys.argv[4] if len(sys.argv) > 4 else None
    labels = ["hetall", "nowat", "none"]
    pids = sorted(
        set.intersection(
            *[{os.path.basename(f)[:4].upper() for f in glob.glob(os.path.join(d, "*.pqr"))} for d in dirs]
        )
    )
    rows, keymismatch, dup_struct = [], 0, []
    for pid in pids:
        got3 = [states(os.path.join(d, pid + ".pqr")) for d in dirs]
        S = [g[0] for g in got3]
        if any(g[1] for g in got3):
            dup_struct.append((pid, max(g[1] for g in got3)))
        keys = set(S[0]) | set(S[1]) | set(S[2])
        if not (set(S[0]) == set(S[1]) == set(S[2])):
            keymismatch += 1
        for k in sorted(keys):
            got = [s.get(k) for s in S]
            if any(g is None for g in got):
                continue
            base = got[0][0]
            st = [g[2] for g in got]
            if st[0] == st[1] == st[2]:
                continue
            # The key is (chain, residue number, residue name). PDB2PQR drops the chain ID by default,
            # so chain is often empty. Leave it empty so later stages can tell whether a chain is present.
            rows.append(
                {
                    "pdb_id": pid,
                    "chain": k[0],
                    "resnum": k[1],
                    "base": base,
                    "state_hetall": st[0],
                    "state_nowat": st[1],
                    "state_none": st[2],
                    "water_removal_changed": int(st[0] != st[1]),
                    "ligand_removal_changed": int(st[1] != st[2]),
                }
            )
    print(
        f"Compared {len(pids)} structures"
        + (f" ({keymismatch} structures with mismatched residue keys)" if keymismatch
           else " (residue keys match in all structures)")
    )
    nochain = sum(1 for r in rows if not r["chain"])
    if nochain:
        print(f"  Note: no chain ID in the pqr for {nochain}/{len(rows)} rows. "
              "PDB2PQR drops chains by default (--keep-chain keeps them)")
    if dup_struct:
        print(f"  !! {len(dup_struct)} structures have a residue with the same atom name twice: {dup_struct[:5]}")
        print("     Same-numbered residues of different chains may have been merged. Rebuild with --keep-chain")
    else:
        print("  No duplicate atom names. No residue merging occurred")
    tot = sum(len(states(os.path.join(dirs[2], p + ".pqr"))[0]) for p in pids)
    print(f"Of {tot} titratable residues, {len(rows)} do not share the same state across the 3 conditions "
          f"({len(rows) / tot * 100:.2f}%)")
    if out and rows:
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    print("\n=== Breakdown of changes (by residue type) ===")
    print(f"{'base':<9}{'water removal':>16}{'ligand removal':>20}{'total':>7}")
    byb = collections.defaultdict(lambda: [0, 0, 0])
    for r in rows:
        b = byb[r["base"]]
        b[0] += r["water_removal_changed"]
        b[1] += r["ligand_removal_changed"]
        b[2] += 1
    for b, (w, l, t) in sorted(byb.items(), key=lambda kv: -kv[1][2]):
        print(f"{b:<9}{w:>16}{l:>20}{t:>7}")
    w_all = sum(r["water_removal_changed"] for r in rows)
    l_all = sum(r["ligand_removal_changed"] for r in rows)
    print(f"{'total':<9}{w_all:>16}{l_all:>20}{len(rows):>7}")
    nstr = len({r["pdb_id"] for r in rows})
    print(f"\nStructures with 1 or more residues that change state {nstr}/{len(pids)} = {nstr / len(pids) * 100:.0f}%")

    print("\n=== From which state to which state (top) ===")
    tr = collections.Counter()
    for r in rows:
        if r["state_hetall"] != r["state_none"]:
            tr[(r["base"], r["state_hetall"], r["state_none"])] += 1
    for (b, a, z), n in tr.most_common(12):
        print(f"  {b}: {a} -> {z}   {n} residues")
    if out:
        print(f"\n-> {out}")


if __name__ == "__main__":
    main()
