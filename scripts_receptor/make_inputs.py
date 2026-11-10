#!/usr/bin/env python3
"""Build the inputs for the 3 preparation conditions from the PDB files as downloaded.

  python3 make_inputs.py <srcdir> <outbase>

Put `*.pdb` in `<srcdir>` (the 100 structures of `full100_manifest.csv` as downloaded
from RCSB). Three directories are created under `<outbase>`.

  in_hetall/   heavy atoms + all HETATM (ligand + water)
  in_nowat/    heavy atoms + HETATM other than water
  recept/      heavy atoms only (both water and ligand removed)

Rules
  Keep only the ATOM and HETATM lines. Drop those whose element column (columns 77-78) is H or D.
  Keep TER and END as they are. Drop every other line (HEADER, REMARK, ANISOU, CONECT and so on).
  Keep the altloc indicator as it is. Do not rewrite the atom numbering.

  The only difference between the conditions is the treatment of HETATM.
    in_hetall  keep all HETATM
    in_nowat   drop HETATM whose residue name is HOH
    recept     drop all HETATM

For a line whose element column is empty, the element is decided from the leading alphabetic character
of the atom name (columns 13-16). A PDB atom name starts at column 13 only when the element symbol has
2 characters, so if column 13 is blank the alphabetic character in column 14 is taken as the element.
"""

import argparse
import glob
import os
import sys

KEEP_OTHER = ("TER", "END")
DROP_ELEM = ("H", "D")

ARMS = (
    ("in_hetall", "all"),
    ("in_nowat", "nowat"),
    ("recept", "none"),
)


def element_of(line):
    """Element symbol. If columns 77-78 are empty, take it from the atom name."""
    e = line[76:78].strip().upper()
    if e:
        return e
    name = line[12:16]
    if name[:1] != " ":
        return name[:2].strip().upper()
    return name[1:2].strip().upper()


def is_hydrogen(line):
    return element_of(line) in DROP_ELEM


def filter_lines(lines, het):
    """het is 'all' | 'nowat' | 'none'."""
    out = []
    for l in lines:
        rec = l[:6]
        if rec.startswith(("ATOM", "HETATM")):
            if is_hydrogen(l):
                continue
            if rec.startswith("HETATM"):
                if het == "none":
                    continue
                if het == "nowat" and l[17:20].strip() == "HOH":
                    continue
            out.append(l)
        elif l.startswith(KEEP_OTHER):
            out.append(l)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("srcdir", help="directory holding the PDB files as downloaded")
    ap.add_argument("outbase", help="where the 3 directories are created")
    a = ap.parse_args()

    files = sorted(glob.glob(os.path.join(a.srcdir, "*.pdb")))
    # glob does not pick up the macOS companion files (._*), but exclude them just in case
    files = [f for f in files if not os.path.basename(f).startswith("._")]
    if not files:
        sys.exit(f"no *.pdb in {a.srcdir}")

    dirs = {}
    for name, _ in ARMS:
        d = os.path.join(a.outbase, name)
        os.makedirs(d, exist_ok=True)
        dirs[name] = d

    n = 0
    for f in files:
        with open(f, errors="replace") as fh:
            lines = fh.readlines()
        for name, het in ARMS:
            out = filter_lines(lines, het)
            with open(os.path.join(dirs[name], os.path.basename(f)), "w") as fh:
                fh.writelines(out)
        n += 1

    print(f"Expanded {n} structures into the 3 conditions")
    for name, _ in ARMS:
        print(f"  {dirs[name]}")


if __name__ == "__main__":
    main()
