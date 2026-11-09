#!/usr/bin/env python3
"""Run PDB2PQR over the sample and count how often a neutral lysine comes out named LYS.

  python3 run_pqr.py <pdbdir> <outdir>

The protocol is the same as in the survey proper:
  pdb2pqr --ff=AMBER --with-ph=7.4 --titration-state-method=propka

A neutral lysine is identified by two tests.
  (1) The number of hydrogens on NZ. A charged LYS has 3, HZ1/HZ2/HZ3, a neutral LYN has 2.
  (2) The residue charge sum. A charged LYS is +1, a neutral LYN is 0.
The AMBER N/CA/C charges could also tell them apart, but the values change at a terminal residue, so they are not used.
Residues where the two tests disagree are counted and reported separately (this validates the tests).
"""

import collections
import csv
import glob
import os
import subprocess
import sys
import time

CMD = ["pdb2pqr", "--ff=AMBER", "--with-ph=7.4", "--titration-state-method=propka"]


def lys_states(pqr):
    """Return (name, HZ count, charge sum) for each LYS/LYN residue in the .pqr."""
    res = collections.OrderedDict()
    for l in open(pqr, errors="replace"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        f = l.split()
        # pqr is whitespace separated rather than column fixed, and the field count drops when coordinates
        # run together. Example: "20.691  54.380-100.128" gives 2 fields, not 3.
        # In 6ZHF, 9137 of the 13614 ATOM lines actually had 9 fields.
        # So the field count must not be used to decide. Take the name and residue from the front, the charge
        # and radius from the back, and take the residue number as the first field from f[4] on that converts
        # to an integer (the chain ID column is not always present, and being alphabetic it fails that test).
        if len(f) < 8:
            continue
        try:
            q = float(f[-2])
            float(f[-1])
        except ValueError:
            continue
        name, rn = f[2], f[3]
        if rn not in ("LYS", "LYN"):
            continue
        resnum = ""
        for t in f[4:7]:
            t2 = t.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")  # insertion code
            if t2.lstrip("-").isdigit():
                resnum = t
                break
        if not resnum:
            continue
        e = res.setdefault((resnum, rn), {"name": rn, "hz": 0, "q": 0.0})
        e["q"] += q
        if name.startswith("HZ"):
            e["hz"] += 1
    return list(res.values())


def main():
    pdbdir, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    rows, fails = [], []
    files = sorted(glob.glob(os.path.join(pdbdir, "*.pdb")))
    for i, f in enumerate(files, 1):
        pid = os.path.basename(f)[:4].upper()
        pqr = os.path.join(outdir, f"{pid}.pqr")
        log = os.path.join(outdir, f"{pid}.pqr.log")
        t0 = time.time()
        if not os.path.exists(pqr):
            p = subprocess.run(CMD + [f, pqr], capture_output=True, text=True, timeout=1800)
            open(log, "w").write(p.stdout + "\n--- stderr ---\n" + p.stderr)
            if p.returncode != 0 or not os.path.exists(pqr):
                tail = [x for x in (p.stderr or p.stdout).strip().split("\n") if x.strip()]
                fails.append(
                    {"pdb_id": pid, "returncode": p.returncode, "error": tail[-1][:300] if tail else "(no output)"}
                )
                print(f"[{i:2}/{len(files)}] {pid} failed rc={p.returncode}", flush=True)
                continue
        st = lys_states(pqr)
        # Test 1: 2 HZ.  Test 2: charge sum near 0
        hz_neutral = [s for s in st if s["hz"] == 2]
        q_neutral = [s for s in st if abs(s["q"]) < 0.5]
        named_lyn = [s for s in st if s["name"] == "LYN"]
        disagree = [s for s in st if (s["hz"] == 2) != (abs(s["q"]) < 0.5)]
        rows.append(
            {
                "pdb_id": pid,
                "n_lys_total": len(st),
                "n_neutral_hz": len(hz_neutral),
                "n_neutral_q": len(q_neutral),
                "n_named_LYN": len(named_lyn),
                "n_neutral_but_named_LYS": len([s for s in hz_neutral if s["name"] == "LYS"]),
                "n_disagree": len(disagree),
                "secs": round(time.time() - t0, 1),
            }
        )
        print(
            f"[{i:2}/{len(files)}] {pid} LYS {len(st):3} neutral {len(hz_neutral):2} "
            f"of which named LYS {len([s for s in hz_neutral if s['name'] == 'LYS']):2} "
            f"({time.time() - t0:.0f}s)",
            flush=True,
        )
    if rows:
        with open(os.path.join(outdir, "lys.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    if fails:
        with open(os.path.join(outdir, "pqr_failed.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(fails[0]))
            w.writeheader()
            w.writerows(fails)
    print(f"\nsucceeded {len(rows)} / failed {len(fails)} -> {outdir}")


if __name__ == "__main__":
    main()
