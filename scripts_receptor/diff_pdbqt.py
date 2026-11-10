#!/usr/bin/env python3
"""Compare two receptor pdbqt files made from the same structure with only the altLoc setting changed.

  python3 diff_pdbqt.py <dirA> <dirB> [out.csv]

This measures not whether a run fails, but how much the contents change in structures that both
settings got through. Compared: atom count, total charge, residues gained or lost, the charge
difference of matched atoms, atom types, and the coordinates of matched atoms.

Charges and coordinates are reported separately because on this sample they do not move together.
The altLoc setting moves coordinates while leaving every charge, atom name, residue name and
residue number untouched, so a comparison that reads only charges reports the two files as the
same when they are not. The identical_charges column says charges and atom identity agree; it
does not say the files agree.

The pdbqt columns are PDB conforming + charge (70-76) + AutoDock atom type (77-79).
Read by fixed columns. Whitespace splitting breaks when coordinates run together (as hit with pqr).
"""

import csv
import glob
import math
import os
import sys


def read(fn):
    """Return a dict (residue name, residue number, atom name) -> (charge, atom type, xyz)."""
    d = {}
    for l in open(fn, errors="replace"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        key = (l[17:20].strip(), l[22:27].strip(), l[12:16].strip())
        try:
            q = float(l[70:76])
        except ValueError:
            q = float("nan")
        try:
            xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        except ValueError:
            xyz = None
        d[key] = (q, l[77:79].strip(), xyz)
    return d


def main():
    if len(sys.argv) < 3 or sys.argv[1] in ("-h", "--help"):
        sys.exit(__doc__)
    dA, dB = sys.argv[1], sys.argv[2]
    out = sys.argv[3] if len(sys.argv) > 3 else None
    rows = []
    for fa in sorted(glob.glob(os.path.join(dA, "*.pdbqt"))):
        pid = os.path.basename(fa)[:4].upper()
        fb = os.path.join(dB, pid + ".pdbqt")
        if not os.path.exists(fb):
            continue
        A, B = read(fa), read(fb)
        onlyA, onlyB = set(A) - set(B), set(B) - set(A)
        common = set(A) & set(B)
        dq = {k: A[k][0] - B[k][0] for k in common}
        moved = {k: v for k, v in dq.items() if abs(v) > 1e-9}
        tdiff = [k for k in common if A[k][1] != B[k][1]]
        # Coordinates of the matched atoms. A displacement is counted when the two files do not
        # print the same string in columns 31-54, which is the resolution the format itself has.
        shifted = {}
        for k in common:
            pa, pb = A[k][2], B[k][2]
            if pa is None or pb is None or pa == pb:
                continue
            shifted[k] = math.dist(pa, pb)
        resA = {(k[0], k[1]) for k in A}
        resB = {(k[0], k[1]) for k in B}
        rows.append(
            {
                "pdb_id": pid,
                "n_atom_A": len(A),
                "n_atom_B": len(B),
                "d_n_atom": len(A) - len(B),
                "sumq_A": round(sum(v[0] for v in A.values()), 4),
                "sumq_B": round(sum(v[0] for v in B.values()), 4),
                "d_sumq": round(sum(v[0] for v in A.values()) - sum(v[0] for v in B.values()), 4),
                "atoms_only_A": len(onlyA),
                "atoms_only_B": len(onlyB),
                "res_only_A": len(resA - resB),
                "res_only_B": len(resB - resA),
                "n_charge_moved": len(moved),
                "max_abs_dq": round(max((abs(v) for v in moved.values()), default=0.0), 4),
                "n_type_diff": len(tdiff),
                "n_coord_moved": len(shifted),
                "max_d_xyz": round(max(shifted.values(), default=0.0), 3),
                "identical_charges": int(not onlyA and not onlyB and not moved and not tdiff),
                "identical_all": int(
                    not onlyA and not onlyB and not moved and not tdiff and not shifted
                ),
            }
        )
    if not rows:
        print("No structures available for comparison")
        return
    if out:
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    n = len(rows)
    same_q = sum(r["identical_charges"] for r in rows)
    same_all = sum(r["identical_all"] for r in rows)
    print(f"{n} structures produced a pdbqt under both settings")
    print(f"  identical in every respect {same_all}, different somewhere {n - same_all}")
    print(f"  identical in charges and atom identity {same_q}, of which coordinates still move"
          f" {same_q - same_all}")
    print(f"  atom count differs {sum(1 for r in rows if r['d_n_atom'])}")
    print(f"  total charge differs {sum(1 for r in rows if abs(r['d_sumq']) > 1e-9)}")
    print(f"  contains an atom whose atom type differs {sum(1 for r in rows if r['n_type_diff'])}")
    print(f"  residues gained or lost {sum(1 for r in rows if r['res_only_A'] or r['res_only_B'])}")
    print(f"  contains an atom whose coordinates move {sum(1 for r in rows if r['n_coord_moved'])}")
    print(f"  atoms moved in total {sum(r['n_coord_moved'] for r in rows)}, "
          f"largest displacement of one atom {max((r['max_d_xyz'] for r in rows), default=0):.2f} A")
    diff = [r for r in rows if not r["identical_all"]]
    if diff:
        print(
            f"\n{'struct':<7}{'atoms A':>7}{'atoms B':>7}{'d atoms':>7}{'sumq A':>9}{'sumq B':>9}"
            f"{'d sumq':>9}{'q moved':>9}{'max|dq|':>9}{'types':>7}{'xyz moved':>11}{'max dxyz':>10}"
        )
        for r in sorted(diff, key=lambda r: -r["n_coord_moved"]):
            print(
                f"{r['pdb_id']:<7}{r['n_atom_A']:>7}{r['n_atom_B']:>7}{r['d_n_atom']:>+7}"
                f"{r['sumq_A']:>9.3f}{r['sumq_B']:>9.3f}{r['d_sumq']:>+9.3f}"
                f"{r['n_charge_moved']:>9}{r['max_abs_dq']:>9.3f}{r['n_type_diff']:>7}"
                f"{r['n_coord_moved']:>11}{r['max_d_xyz']:>10.2f}"
            )
    if out:
        print(f"\n-> {out}")


if __name__ == "__main__":
    main()
