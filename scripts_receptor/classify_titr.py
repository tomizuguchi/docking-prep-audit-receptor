#!/usr/bin/env python3
"""Decide the state of the titratable residues in a pqr from the presence of hydrogens, not from the residue name.

  python3 classify_titr.py <pqrdir> [out.csv]

Why the residue name is not trusted: PDB2PQR leaves the name as LYS even when it assigns the LYN
parameters to a neutral lysine (18/18 in the raw files of the 100 structures).
Read by name, not a single neutral lysine is visible. Hydrogens are read instead.

Decision (AMBER hydrogen naming)
  ASP  HD2 present       -> neutral (ASH)    absent -> charged (-1)
  GLU  HE2 present       -> neutral (GLH)    absent -> charged (-1)
  HIS  both HD1 and HE2  -> HIP(+1)
       HD1 only          -> HID (neutral)    HE2 only -> HIE (neutral)
  CYS  HG present        -> neutral thiol    absent -> SS bond or thiolate (-1)
       SS bond when it has a partner with SG-SG under 2.5 A, thiolate when it does not
  TYR  HH present        -> neutral          absent -> tyrosinate (-1)
  LYS  3 HZ              -> charged (+1)     2 -> neutral (LYN)
  ARG  5 H on NE/NH1/NH2 -> charged (+1)     4 -> neutral (AR0)
Residues that could not be classified are not thrown away; they are counted as unknown and always reported.

Parser note: pqr is whitespace-separated, but the field count drops when coordinates run together
(in 6ZHF, 67% of the ATOM lines have 9 fields). Lines must not be rejected by field count.
Coordinates are read from fixed columns (30-54). An assert guarantees that not one ATOM line of a
titratable residue is dropped.
"""

import collections
import csv
import glob
import os
import sys

import numpy as np

TITR = {
    "ASP",
    "ASH",
    "GLU",
    "GLH",
    "HIS",
    "HID",
    "HIE",
    "HIP",
    "CYS",
    "CYM",
    "CYX",
    "LYS",
    "LYN",
    "TYR",
    "TYM",
    "ARG",
    "AR0",
}
# Name to the base residue type
BASE = {
    "ASH": "ASP",
    "GLH": "GLU",
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "CYM": "CYS",
    "CYX": "CYS",
    "LYN": "LYS",
    "TYM": "TYR",
    "AR0": "ARG",
}


def parse(pqr):
    """Return the information for each titratable residue and the line counts for the drop cross-check."""
    res, seen, placed = collections.OrderedDict(), 0, 0
    for l in open(pqr, errors="replace"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        f = l.split()
        if len(f) < 8:
            continue
        rn = f[3]
        if rn not in TITR:
            continue
        seen += 1
        try:
            q = float(f[-2])
            float(f[-1])
        except ValueError:
            continue
        # The chain ID goes between the residue name and the residue number. PDB2PQR drops it by
        # default, but it stays when --keep-chain is given. Put it in the key when present.
        chain, num = "", ""
        for j in range(4, min(7, len(f))):
            if f[j].rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ").lstrip("-").isdigit():
                num = f[j]
                if j > 4:
                    chain = f[j - 1]
                break
        if not num:
            continue
        placed += 1
        name = f[2]
        e = res.setdefault(
            (chain, num, rn), {"base": BASE.get(rn, rn), "name": rn, "at": set(), "q": 0.0, "sg": None, "dup": 0}
        )
        # The same atom name coming twice means same-numbered residues from different chains were
        # merged, or altLoc. Putting it in the set and dropping it silently makes HD1 from one and
        # HE2 from the other look like one residue, and a spurious HIP(+1) is created. Count them and
        # report to the caller.
        if name in e["at"]:
            e["dup"] += 1
        e["at"].add(name)
        e["q"] += q
        if name == "SG":
            try:
                e["sg"] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            except ValueError:
                pass
    return res, seen, placed


def classify(res):
    """Decide the state of each residue. CYS needs the SG-SG distance, so the whole set is looked at."""
    sgs = [(k, v["sg"]) for k, v in res.items() if v["base"] == "CYS" and v["sg"] is not None]
    ss = set()
    for i in range(len(sgs)):
        for j in range(i + 1, len(sgs)):
            if np.linalg.norm(sgs[i][1] - sgs[j][1]) < 2.5:
                ss.add(sgs[i][0])
                ss.add(sgs[j][0])
    out = {}
    for k, v in res.items():
        a, b = v["at"], v["base"]
        if b == "ASP":
            s = "neutral" if "HD2" in a else "charged-"
        elif b == "GLU":
            s = "neutral" if "HE2" in a else "charged-"
        elif b == "HIS":
            d, e = "HD1" in a, "HE2" in a
            s = "HIP+" if (d and e) else "HID" if d else "HIE" if e else "unknown"
        elif b == "CYS":
            s = "neutral-SH" if "HG" in a else ("disulfide" if k in ss else "thiolate-")
        elif b == "TYR":
            s = "neutral" if "HH" in a else "tyrosinate-"
        elif b == "LYS":
            n = sum(1 for x in a if x.startswith("HZ"))
            s = "charged+" if n == 3 else "neutral" if n == 2 else f"unknown(HZ={n})"
        elif b == "ARG":
            n = sum(1 for x in a if x in ("HE", "HH11", "HH12", "HH21", "HH22"))
            s = "charged+" if n == 5 else "neutral" if n == 4 else f"unknown(H={n})"
        else:
            s = "unknown"
        out[k] = (b, v["name"], s)
    return out


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        sys.exit(__doc__)
    pqrdir = sys.argv[1]
    if not glob.glob(os.path.join(pqrdir, "*.pqr")):
        sys.exit(f"no *.pqr in {pqrdir}")
    out = sys.argv[2] if len(sys.argv) > 2 else None
    rows = []
    for f in sorted(glob.glob(os.path.join(pqrdir, "*.pqr"))):
        pid = os.path.basename(f)[:4].upper()
        res, seen, placed = parse(f)
        assert seen == placed, f"{pid}: only {placed} of {seen} titratable-residue lines were assigned"
        # The key carries the chain ID as well, which is empty unless PDB2PQR was given --keep-chain.
        # It is not written to the CSV; the columns below are the ones the deposited table has.
        for (_chain, num, rn), (base, name, state) in classify(res).items():
            rows.append({"pdb_id": pid, "resnum": num, "base": base, "name_in_pqr": name, "state": state})
    if out:
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    print(f"{len({r['pdb_id'] for r in rows})} structures, {len(rows)} titratable residues (zero drops cross-checked)")
    unk = [r for r in rows if r["state"].startswith("unknown")]
    print(
        f"{len(unk)} residues could not be classified"
        + (f": {collections.Counter((r['base'], r['state']) for r in unk).most_common(6)}" if unk else "")
    )
    print("\n=== Whether PDB2PQR reflects the state in the name ===")
    print(f"{'base':<9}{'state':<14}{'count':>7}  name in pqr")
    tab = collections.defaultdict(collections.Counter)
    for r in rows:
        tab[(r["base"], r["state"])][r["name_in_pqr"]] += 1
    for (base, state), names in sorted(tab.items()):
        shown = ", ".join(f"{n} {c}" for n, c in names.most_common())
        hidden = (
            "  <- not distinguishable by name"
            if len(names) == 1 and state not in ("charged+", "charged-", "disulfide")
            else ""
        )
        print(f"{base:<9}{state:<14}{sum(names.values()):>7}  {shown}{hidden}")
    if out:
        print(f"\n-> {out}")


if __name__ == "__main__":
    main()
