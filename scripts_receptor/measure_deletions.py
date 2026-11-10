#!/usr/bin/env python3
"""List the residues that Meeko's -a (--allow_bad_res / --delete_bad_res) flag deleted, and summarize them.

  python3 measure_deletions.py <recept_dir> <meeko_dir> <prof100.csv> [--noa <meeko_noa_dir>] [--out deleted.csv]

A residue counts as deleted when it has ATOM records in the input PDB file and no atom in the pdbqt
file Meeko wrote from it. Residues are matched on chain, residue number with insertion code, and
residue name.

Reported:
  - the number of residues deleted, the number of structures affected, and the largest count
  - the residue types deleted
  - the share of charged residues (ARG, LYS, ASP, GLU) among the deleted residues, against their share
    among all residues of the inputs on which Meeko completed
  - the formal charge removed from each structure. The sign convention matters: the column
    q_deleted is the sum of the formal charges of the DELETED residues (ARG and LYS +1, ASP and GLU -1),
    so the charge of the receptor that remains moves by -q_deleted. Deleting more lysines than
    aspartates makes q_deleted positive and the receptor more negative.
  - the number deleted against the number of incomplete side chains in prof100.csv (incomplete_sc)
  - with --noa, whether the structures that fail without the flag and complete with it are the same
    set as the structures in which residues are deleted

The pdbqt columns are read by fixed position, as in diff_pdbqt.py.
"""

import argparse
import collections
import csv
import glob
import os

CHARGE = {"ARG": 1, "LYS": 1, "ASP": -1, "GLU": -1}
FIELDS = ["pdb_id", "chain", "resnum", "resname", "formal_charge"]


def residues(path, het=False):
    """Ordered dict (chain, resnum+icode, resname) -> True for the residues with a coordinate record."""
    d = collections.OrderedDict()
    for l in open(path, errors="replace"):
        if l.startswith("ATOM") or (het and l.startswith("HETATM")):
            d[(l[21], l[22:27].strip(), l[17:20].strip())] = True
    return d


def has_atoms(path):
    return os.path.exists(path) and any(l.startswith("ATOM") for l in open(path, errors="replace"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("recept_dir", help="the input PDB files given to Meeko")
    ap.add_argument("meeko_dir", help="the pdbqt files written with the flag")
    ap.add_argument("prof_csv", help="prof100.csv, for the incomplete_sc column")
    ap.add_argument("--noa", help="the pdbqt directory of the run without the flag")
    ap.add_argument("--out", help="CSV to write one row per deleted residue to")
    a = ap.parse_args()

    prof = {r["id"].upper(): int(r["incomplete_sc"] or 0) for r in csv.DictReader(open(a.prof_csv))}
    inputs = sorted(os.path.basename(f)[:4].upper() for f in glob.glob(os.path.join(a.recept_dir, "*.pdb")))
    done = [p for p in inputs if has_atoms(os.path.join(a.meeko_dir, p + ".pdbqt"))]

    rows, per, n_res, n_chg = [], {}, 0, 0
    for pid in done:
        inp = residues(os.path.join(a.recept_dir, pid + ".pdb"))
        kept = residues(os.path.join(a.meeko_dir, pid + ".pdbqt"), het=True)
        n_res += len(inp)
        n_chg += sum(1 for k in inp if k[2] in CHARGE)
        gone = [k for k in inp if k not in kept]
        per[pid] = gone
        for ch, num, rn in gone:
            rows.append(
                {"pdb_id": pid, "chain": ch.strip(), "resnum": num, "resname": rn, "formal_charge": CHARGE.get(rn, 0)}
            )
    if a.out:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS)
            w.writeheader()
            w.writerows(rows)

    affected = {p: g for p, g in per.items() if g}
    n_del = len(rows)
    print(f"{len(inputs)} inputs, Meeko completed on {len(done)}")
    print(f"residues deleted: {n_del} in {len(affected)}/{len(done)} structures")
    if not affected:
        return
    big = max(affected, key=lambda p: len(affected[p]))
    tot_big = len(residues(os.path.join(a.recept_dir, big + ".pdb")))
    share = 100 * len(affected[big]) / tot_big
    print(f"  largest: {len(affected[big])} in {big} ({share:.1f}% of its {tot_big} residues)")
    types = collections.Counter(r["resname"] for r in rows)
    print("  by residue type: " + ", ".join(f"{k} {v}" for k, v in types.most_common()))

    del_chg = sum(1 for r in rows if r["resname"] in CHARGE)
    pos = types["ARG"] + types["LYS"]
    neg = types["ASP"] + types["GLU"]
    share_all = n_chg / n_res
    share_del = del_chg / n_del
    print(f"\ncharged residues (ARG, LYS, ASP, GLU) among all residues of the {len(done)} inputs: "
          f"{n_chg}/{n_res} = {100 * share_all:.1f}%")
    print(f"charged residues among the deleted: {del_chg}/{n_del} = {100 * share_del:.1f}%"
          f"   enrichment {share_del / share_all:.2f}")
    print(f"  LYS + ARG deleted {pos}, ASP + GLU deleted {neg}")

    q = {p: sum(CHARGE.get(k[2], 0) for k in g) for p, g in affected.items()}
    moved = {p: v for p, v in q.items() if v}
    print(f"\nformal charge of the receptor moves in {len(moved)}/{len(done)} structures")
    print("  q_deleted = charge of the deleted residues; the receptor moves by the opposite amount")
    for p, v in sorted(moved.items(), key=lambda kv: -abs(kv[1]))[:5]:
        print(f"    {p}  q_deleted {v:+d}   receptor {-v:+d}")
    toward_neg = sum(1 for v in moved.values() if v > 0)
    print(f"  receptor moves towards negative in {toward_neg}, towards positive in {len(moved) - toward_neg}")

    same = [p for p in affected if len(affected[p]) == prof.get(p, 0)]
    more = {p: len(affected[p]) - prof.get(p, 0) for p in affected if len(affected[p]) != prof.get(p, 0)}
    n_inc = sum(1 for p in inputs if prof.get(p, 0))
    print(f"\nincomplete side chains in the inputs (prof100.csv): {sum(prof.get(p, 0) for p in inputs)}"
          f" in {n_inc} structures")
    print(f"  structures in which the number deleted equals incomplete_sc: {len(same)}/{len(affected)}")
    if more:
        print("  differs in: " + ", ".join(f"{p} {v:+d}" for p, v in sorted(more.items())))
        print(f"  sum of the differences: {sum(more.values()):+d}")

    if a.noa:
        ok_noa = {p for p in inputs if has_atoms(os.path.join(a.noa, p + ".pdbqt"))}
        rescued = set(done) - ok_noa
        print(f"\nwithout the flag: {len(ok_noa)} complete, {len(inputs) - len(ok_noa)} fail")
        print(f"  fail without the flag and complete with it: {len(rescued)}")
        print(f"  same set as the structures with deletions: {rescued == set(affected)}")
        still = set(inputs) - set(done)
        print(f"  fail with the flag as well: {sorted(still)}")
    if a.out:
        print(f"\n-> {a.out}")


if __name__ == "__main__":
    main()
