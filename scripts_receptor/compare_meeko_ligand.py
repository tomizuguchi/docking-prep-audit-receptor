#!/usr/bin/env python3
"""Compare Meeko's output on the input that keeps the ligand with its output on the input without it.

  python3 compare_meeko_ligand.py <meeko_recept_dir> <meeko_nowat_dir> <recept_dir> <in_nowat_dir>

For every structure on which both runs wrote a pdbqt file, the protein part of the two files is
compared: the atoms that the file from the ligand-free input contains, matched on chain, residue
number, residue name and atom name. Reported per structure and in total:

  - whether the two input PDB files are the same file (a structure with no non-water HETATM record
    has identical recept and in_nowat inputs, so for it the comparison is not a test)
  - whether the two pdbqt files are byte-identical
  - otherwise, the number of protein atoms missing from either side, the number whose charge or
    AutoDock atom type differs, and the number of atoms present only in the ligand-kept output

This is the comparison behind the statement that, where it could be made, keeping the ligand did not
change what Meeko assigned to the protein.
"""

import glob
import hashlib
import os
import sys


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def atoms(path):
    d = {}
    for l in open(path, errors="replace"):
        if l.startswith(("ATOM", "HETATM")):
            d[(l[21], l[22:27].strip(), l[17:20].strip(), l[12:16].strip())] = (l[70:76].strip(), l[77:79].strip())
    return d


def main():
    if len(sys.argv) != 5 or sys.argv[1] in ("-h", "--help"):
        sys.exit(__doc__)
    m_rec, m_now, d_rec, d_now = sys.argv[1:5]
    both = sorted(
        {os.path.basename(f)[:4] for f in glob.glob(os.path.join(m_rec, "*.pdbqt")) if os.path.getsize(f)}
        & {os.path.basename(f)[:4] for f in glob.glob(os.path.join(m_now, "*.pdbqt")) if os.path.getsize(f)}
    )
    same_input = identical = protein_same = 0
    with_ligand = []
    print(f"{'struct':<7}{'same input':>11}{'same pdbqt':>11}{'missing':>9}{'q/type differ':>15}{'extra atoms':>13}")
    for pid in both:
        si = sha(os.path.join(d_rec, pid + ".pdb")) == sha(os.path.join(d_now, pid + ".pdb"))
        fa, fb = os.path.join(m_rec, pid + ".pdbqt"), os.path.join(m_now, pid + ".pdbqt")
        so = sha(fa) == sha(fb)
        A, B = atoms(fa), atoms(fb)
        missing = len(set(A) - set(B))
        differ = sum(1 for k in set(A) & set(B) if A[k] != B[k])
        extra = len(set(B) - set(A))
        same_input += si
        identical += so
        protein_same += missing == 0 and differ == 0
        if not si:
            with_ligand.append(pid)
        print(f"{pid:<7}{si!s:>11}{so!s:>11}{missing:>9}{differ:>15}{extra:>13}")
    print(f"\n{len(both)} structures completed under both inputs")
    print(f"  identical input files (no ligand in the structure): {same_input}")
    print(f"  input files differ (a ligand is present): {len(with_ligand)}  {', '.join(with_ligand)}")
    print(f"  byte-identical pdbqt files: {identical}")
    print(f"  protein part identical in atoms, charges and atom types: {protein_same}/{len(both)}")


if __name__ == "__main__":
    main()
