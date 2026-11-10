#!/usr/bin/env python3
"""Check that PDB2PQR run with --ffout=AMBER names each titratable residue for its protonation state.

  python3 check_ffout.py <pdb_dir> <out_dir> <ID> [<ID> ...]

For each structure, PDB2PQR is run as in run_pqr.py with --ffout=AMBER added:

  pdb2pqr --ff=AMBER --ffout=AMBER --with-ph=7.4 --titration-state-method=propka <in.pdb> <out.pqr>

The state of every titratable residue is then read from the hydrogens present (classify_titr.py),
independently of the name, and compared with the name the file gives the residue. The names expected
are the AMBER ones: ASH and GLH for the neutral acids, LYN for neutral lysine, HID, HIE and HIP for
histidine, CYX for a disulfide cysteine and CYM for a thiolate.

The manuscript reports this check for 1T5S, 3DX2 and 1EA7, and for the thiolate name on 4NV6.
"""

import collections
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import classify_titr

EXPECT = {
    ("ASP", "neutral"): "ASH",
    ("ASP", "charged-"): "ASP",
    ("GLU", "neutral"): "GLH",
    ("GLU", "charged-"): "GLU",
    ("LYS", "neutral"): "LYN",
    ("LYS", "charged+"): "LYS",
    ("HIS", "HIP+"): "HIP",
    ("HIS", "HID"): "HID",
    ("HIS", "HIE"): "HIE",
    ("CYS", "disulfide"): "CYX",
    ("CYS", "thiolate-"): "CYM",
    ("CYS", "neutral-SH"): "CYS",
}


def main():
    if len(sys.argv) < 4 or sys.argv[1] in ("-h", "--help"):
        sys.exit(__doc__)
    pdb_dir, out_dir, pids = sys.argv[1], sys.argv[2], sys.argv[3:]
    os.makedirs(out_dir, exist_ok=True)
    bad_total = 0
    for pid in pids:
        src = os.path.join(pdb_dir, pid + ".pdb")
        dst = os.path.join(out_dir, pid + ".pqr")
        p = subprocess.run(
            ["pdb2pqr", "--ff=AMBER", "--ffout=AMBER", "--with-ph=7.4", "--titration-state-method=propka", src, dst],
            capture_output=True,
            text=True,
        )
        if not os.path.exists(dst):
            print(f"{pid}: pdb2pqr wrote no file (rc={p.returncode})")
            bad_total += 1
            continue
        res = classify_titr.parse(dst)[0]
        seen, bad = collections.Counter(), []
        for (_chain, num, _rn), (base, name, state) in classify_titr.classify(res).items():
            want = EXPECT.get((base, state))
            if want is None:
                continue
            seen[(base, state, name)] += 1
            if name != want:
                bad.append(f"{base}{num} state {state} named {name}, expected {want}")
        bad_total += len(bad)
        print(f"{pid}: {sum(seen.values())} residues checked, {len(bad)} named otherwise than their state")
        for (base, state, name), c in sorted(seen.items()):
            print(f"    {base} {state:<11} -> {name}  x{c}")
        for b in bad[:10]:
            print(f"    !! {b}")
    sys.exit(1 if bad_total else 0)


if __name__ == "__main__":
    main()
