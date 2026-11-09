#!/usr/bin/env python3
"""Inspect the input before docking receptor preparation and cross-check the output afterwards.

  python3 check_receptor.py <input.pdb> [--altloc A] [--allow-bad-res]
                              [--prepared <out.pqr|out.pdbqt>] [--json]

  Specify --allow-bad-res when -a is going to be given to Meeko.
  The report of incomplete side chains drops from ERROR to WARN (nothing crashes, but residues disappear).

Every rule checked here is backed by measurements on a stratified sample of 100 structures.
No dependencies (standard library only).

Rules where preparation crashes (zero false positives and zero misses in the measurements)
  E1  A residue lacks both C and O -> PDB2PQR cannot cap it and crashes (3/3)
  E2  A heavy atom has neither X nor blank for --default_altloc X
      -> Meeko crashes (1/1 with X=A, 2/2 with X=B)

Rules where preparation silently changes the receptor
  W1  A metal is present. Dropping HETATM and then titrating flips the charge of the coordinating residues
      (sum q moves in 23 of the 41 structures with a metal, by up to 4 units; without a metal at most 1 unit.
      41 is the subset of the 44 metal-containing structures that PDB2PQR produced output for in both arms,
      which is where sum q can be compared)
      All titratable atoms are examined. When the scan stopped at the first one, only 73 residues inside the
      coordination sphere were picked up in the 44 structures containing a metal (the correct count is 124.
      The misses break down as 23 NE2 of HIS, 18 OD2 of ASP, 9 OE2 of GLU and 1 NH2 of ARG. Only ND1 of HIS
      was examined)
  W2  HETATM other than water is present. The removal order changes the titration state (64 residues / 34 structures)
  W3  Deposited hydrogens are present. PDB2PQR can misreport a chain break and crash (1EA7)
  W4  A modified residue inside the chain comes on a HETATM line. Removing it breaks the main chain (6/100)
      The decision uses the C-N distance to the neighbouring residue. When a list of residue names was used,
      only MSE was picked up and OSE and CSS were missed, and conversely there was a risk of picking up
      SAM or free PRO, which merely carry main chain atom names
  W5  A residue has an incomplete side chain. Without -a Meeko crashes (37/37, 0 false alarms).
      With -a it passes but this residue is deleted, with a 2.9-fold bias toward charged residues
      (306 residues / 38 structures, the formal charge moves by up to 9 units).
      Of the 39 failures without -a, 37 are this, 1 is altLoc (E2) and 1 (4LYZ) is unexplained

Cross-checks on the output
  C1  The residue names in the pqr do not express the state (--ffout was forgotten)
  C2  Residues are missing from the pdbqt (deletion by Meeko's -a). The formal charge of the loss is also given
  C3  Residues whose titration state differs between two pqr files (--prepared A.pqr --against B.pqr)
      Run PDB2PQR on both the input with HETATM and the input without it and compare.
      The pqr records no provenance (no header and no REMARK), so one file alone cannot tell
      which condition the titration was run under. The only trace is a single hydrogen on the changed
      residue, and it is indistinguishable from the correct state.
      In the sample of 100 structures the state moved for 64 residues / 34 structures (matches arms_diff_3way.csv)
"""

import collections
import json
import math
import os
import sys

METALS = set(
    [
        "LI",
        "NA",
        "K",
        "RB",
        "CS",
        "MG",
        "CA",
        "SR",
        "BA",
        "SC",
        "Y",
        "LA",
        "CE",
        "TI",
        "ZR",
        "HF",
        "V",
        "NB",
        "TA",
        "CR",
        "MO",
        "W",
        "MN",
        "TC",
        "RE",
        "FE",
        "RU",
        "OS",
        "CO",
        "RH",
        "IR",
        "NI",
        "PD",
        "PT",
        "CU",
        "AG",
        "AU",
        "ZN",
        "CD",
        "HG",
        "AL",
        "GA",
        "IN",
        "TL",
        "SN",
        "PB",
        "BI",
        "U",
        "TH",
        "EU",
        "GD",
        "TB",
        "DY",
        "HO",
        "ER",
        "TM",
        "YB",
        "LU",
        "SM",
        "ND",
        "PR",
    ]
)
STD = set(
    [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    ]
)
# Modified residues that can appear inside a chain (they often come on HETATM lines)
INCHAIN_MOD = {"MSE", "SEP", "TPO", "PTR", "CSO", "CME", "KCX", "MLY", "ALY", "HYP", "PCA", "LLP"}
BB = {"N", "CA", "C", "O"}
# Heavy atom count of a complete side chain (includes main chain N CA C O)
_SIDE = {
    "ALA": "CB",
    "ARG": "CB CG CD NE CZ NH1 NH2",
    "ASN": "CB CG OD1 ND2",
    "ASP": "CB CG OD1 OD2",
    "CYS": "CB SG",
    "GLN": "CB CG CD OE1 NE2",
    "GLU": "CB CG CD OE1 OE2",
    "GLY": "",
    "HIS": "CB CG ND1 CD2 CE1 NE2",
    "ILE": "CB CG1 CG2 CD1",
    "LEU": "CB CG CD1 CD2",
    "LYS": "CB CG CD CE NZ",
    "MET": "CB CG SD CE",
    "PHE": "CB CG CD1 CD2 CE1 CE2 CZ",
    "PRO": "CB CG CD",
    "SER": "CB OG",
    "THR": "CB OG1 CG2",
    "TRP": "CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
    "TYR": "CB CG CD1 CD2 CE1 CE2 CZ OH",
    "VAL": "CB CG1 CG2",
}
# Set of heavy atom names of a complete residue (includes main chain N CA C O).
# Compare by set, not by count. With a count, a LYS that lacks NZ but carries OXT
# is judged "complete" at 9 atoms and W5 gives a false alarm.
ATOMSET = {rn: BB | set(s.split()) for rn, s in _SIDE.items()}
NHEAVY = {rn: len(s) for rn, s in ATOMSET.items()}
TITR_ATOM = {
    "ASP": ("OD1", "OD2"),
    "GLU": ("OE1", "OE2"),
    "HIS": ("ND1", "NE2"),
    "LYS": ("NZ",),
    "CYS": ("SG",),
    "TYR": ("OH",),
    "ARG": ("NE", "NH1", "NH2"),
}
FORMAL = {"LYS": 1, "ARG": 1, "ASP": -1, "GLU": -1}


def dist(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def guess_element(l):
    """(element symbol, whether it could be read from the element columns).

    In a PDB where columns 77-78 are empty there is no choice but to infer from the atom name, but taking
    the first character turns ZN into 'Z' and loses the metal. When the 2 characters right-justified in
    columns 13-14 match a metal name, those are taken instead (" CA " is C-alpha, "CA  " is calcium).
    The number of atoms that fell back to inference is counted and warned about by the caller."""
    el = l[76:78].strip().upper()
    if el:
        return el, True
    name = l[12:16] if len(l) >= 16 else l[12:].ljust(4)
    two = name[0:2].strip().upper()
    if len(two) == 2 and two in METALS:
        return two, False
    raw = name.strip().upper()
    if raw and raw[0].isdigit():
        raw = raw[1:]  # PDB v2 hydrogen name such as 1HB
    return (raw[0] if raw else ""), False


def find_inchain(hetres, chainres):
    """Return those residues coming on HETATM lines that are linked to the chain by a peptide bond.

    A fixed list of residue names is not enough. In the sample of 100 structures the list picked up only
    MSE and missed OSE / 0QF / DSN / CSS. Conversely SAM and free PRO carry main chain atom names
    (N, CA, C) but are not linked to the chain, and the name or the atom names alone cannot tell them
    apart. Deciding by the C-N distance separates both cases correctly."""
    PEP = 1.5  # the peptide bond C-N is about 1.33 A
    out = {}
    for k, e in hetres.items():
        n, c = e["atoms"].get("N"), e["atoms"].get("C")
        linked = False
        for k2, e2 in chainres.items():
            if k2 == k:
                continue
            c2, n2 = e2["atoms"].get("C"), e2["atoms"].get("N")
            if n and c2 and dist(n, c2) < PEP:
                linked = True
                break
            if c and n2 and dist(c, n2) < PEP:
                linked = True
                break
        if linked:
            out[k] = e
    return out


def read_pdb(path):
    """Read ATOM/HETATM. altLoc is held as a set per atom name."""
    res = collections.OrderedDict()  # (chain, resseq+icode, resname) -> info
    hetres = collections.OrderedDict()
    het, metals, nH, guessed = [], [], 0, 0
    for l in open(path, errors="replace"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        el, from_col = guess_element(l)
        if not from_col:
            guessed += 1
        name = l[12:16].strip()
        alt = l[16] if len(l) > 16 else " "
        rn = l[17:20].strip()
        key = (l[21] if len(l) > 21 else " ", l[22:27].strip(), rn)
        try:
            xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        except ValueError:
            continue
        if el == "H" or el == "D":
            nH += 1
            continue  # count heavy atoms only (the receptor file drops H)
        if l.startswith("HETATM"):
            if rn != "HOH":
                het.append((rn, name, xyz))
                if el in METALS:
                    metals.append((rn, el, xyz))
            # A modified residue inside the chain comes on a HETATM line but is part of the main chain.
            # They are all collected separately here first and sorted out later by bond distance.
            if rn != "HOH":
                e = hetres.setdefault(key, {"atoms": {}, "alt": collections.defaultdict(set)})
                e["atoms"][name] = xyz
                e["alt"][name].add(alt)
            continue
        e = res.setdefault(key, {"atoms": {}, "alt": collections.defaultdict(set)})
        e["atoms"][name] = xyz
        e["alt"][name].add(alt)
    # Merge the HETATM residues linked to the chain into res. Without putting them in res,
    # the checks for E1 (lacking both C and O) and W5 (incomplete side chain) do not reach them.
    inchain = find_inchain(hetres, res)
    for k, e in inchain.items():
        res[k] = e
    return res, het, metals, nH, guessed, inchain


def _min_titr_dist(e, rn, targets):
    """For all titratable atoms of a residue, the shortest distance to the targets and the name of that atom.

    Previously only the first atom found was examined and the loop broke. For HIS only ND1 was examined
    and not NE2, so a histidine coordinating a metal through NE2 was classified as outside the
    coordination sphere. For ARG only NE was examined and not NH1 / NH2."""
    best, who = None, None
    for a in TITR_ATOM.get(rn, ()):
        if a not in e["atoms"]:
            continue
        d = min(dist(e["atoms"][a], xyz) for xyz in targets)
        if best is None or d < best:
            best, who = d, a
    return best, who


def check_input(path, altlocs=("A", "B"), allow_bad_res=False):
    res, het, metals, nH, guessed, inchain = read_pdb(path)
    out = []
    if guessed:
        out.append(
            (
                "W0",
                "WARN",
                f"{guessed} atoms have an empty element column (digits 77-78). The element was inferred from the name",
                ["inference can miss a metal"],
                "Use a file with the element column filled in. The metal detection below is less reliable",
            )
        )

    # --- E1: residues lacking both C and O ---
    cap = []
    for (ch, num, rn), e in res.items():
        if rn not in STD and rn not in INCHAIN_MOD:
            continue
        miss = BB - set(e["atoms"])
        if {"C", "O"} <= miss:
            cap.append(f"{rn}{num}({ch}) missing {'+'.join(sorted(miss))}")
    if cap:
        out.append(
            (
                "E1",
                "ERROR",
                "PDB2PQR cannot cap a residue that lacks both C and O and crashes",
                cap,
                "Delete that residue, or complete the main chain with an external tool",
            )
        )

    # --- E2: heavy atoms with no altLoc variant ---
    for X in altlocs:
        bad = [
            f"{rn}{num}({ch}):{a}"
            for (ch, num, rn), e in res.items()
            for a, s in e["alt"].items()
            if not (s & {" ", X})
        ]
        if bad:
            out.append(
                (
                    f"E2[{X}]",
                    "ERROR",
                    f"Some atoms are not found with --default_altloc {X}. Meeko crashes",
                    bad,
                    "Specify them individually with --wanted_altloc, or make another altLoc the default",
                )
            )
    ok = [X for X in altlocs if not any(not (s & {" ", X}) for e in res.values() for s in e["alt"].values())]
    if len(altlocs) > 1:
        out.append(
            (
                "I0",
                "INFO",
                "safe --default_altloc",
                [", ".join(ok) if ok else "this structure has no safe default altLoc"],
                "",
            )
        )

    # --- W1 / W2: HETATM and metals ---
    if metals:
        kinds = collections.Counter(e for _, e, _ in metals)
        coord, nearby = [], []
        mxyz = [xyz for _, _, xyz in metals]
        for (ch, num, rn), e in res.items():
            d, a = _min_titr_dist(e, rn, mxyz)
            if d is None:
                continue
            if d < 2.6:
                coord.append(f"{rn}{num}({ch}) {a} {d:.2f} A to the metal")
            elif d < 6.0:
                nearby.append(f"{rn}{num}({ch}) {a} {d:.1f} A")
        out.append(
            (
                "W1",
                "WARN",
                f"Metals {dict(kinds)} are present. Dropping HETATM then titrating flips coordinating residue charges",
                [f"{len(coord)} titratable residues inside the coordination sphere (under 2.6 A)"]
                + coord
                + [f"{len(nearby)} titratable residues at 2.6-6 A"],
                "Run PDB2PQR with the metal left in the input. "
                "No removal step is needed. PDB2PQR does not write HETATM other than water to the output, "
                "so leaving it in the input does not make it appear in the pqr",
            )
        )
    if het:
        ligs = sorted({rn for rn, _, _ in het})
        near = []
        hxyz = [xyz for _, _, xyz in het]
        for (ch, num, rn), e in res.items():
            d, a = _min_titr_dist(e, rn, hxyz)
            if d is not None and d < 6.0:
                near.append(f"{rn}{num}({ch}) {a} {d:.1f} A")
        out.append(
            (
                "W2",
                "WARN",
                f"HETATM other than water is present ({', '.join(ligs[:8])}"
                f"{' and others' if len(ligs) > 8 else ''}). The removal order changes the titration state",
                [f"{len(near)} titratable residues within 6 A of a ligand"],
                "Run PDB2PQR on both the input with HETATM and the input without it, "
                "and check the difference with --prepared A.pqr --against B.pqr. "
                "The pqr records no provenance, so one file alone cannot decide",
            )
        )

    # --- W3: deposited hydrogens ---
    if nH:
        out.append(
            (
                "W3",
                "WARN",
                f"{nH} deposited hydrogens are present. PDB2PQR can misreport a chain break and crash",
                [],
                "Strip the hydrogens before passing the file. Each tool adds them back on its own",
            )
        )

    # --- W4: modified residues inside the chain on HETATM lines ---
    if inchain:
        kinds = sorted({k[2] for k in inchain})
        out.append(
            (
                "W4",
                "WARN",
                f"Modified residues inside the chain come on HETATM lines ({', '.join(kinds)})",
                [f"{k[2]}{k[1]}({k[0]})" for k in list(inchain)[:12]]
                + ([f"... and {len(inchain) - 12} more"] if len(inchain) > 12 else []),
                "Deleting all HETATM at once creates a break in the main chain. Keep these residues",
            )
        )

    # --- W5: residues with an incomplete side chain ---
    inc, q_at_risk = [], 0
    for (ch, num, rn), e in res.items():
        want = ATOMSET.get(rn)
        if not want:
            continue
        # Compare by set, not by count. This keeps an extra atom such as OXT from compensating
        # for a missing one and causing a miss.
        miss = want - set(e["atoms"])
        if miss:
            inc.append(f"{rn}{num}({ch}) missing {'+'.join(sorted(miss))}")
            q_at_risk += FORMAL.get(rn, 0)
    if inc:
        # Without -a Meeko crashes (37/37 over 100 structures, 0 false alarms). With -a it passes, but
        # this residue is deleted (306 residues over 38 structures, 2.9-fold bias toward charged residues).
        lvl = "WARN" if allow_bad_res else "ERROR"
        msg = (
            f"{len(inc)} residues have an incomplete side chain. Meeko's -a deletes them"
            if allow_bad_res
            else f"{len(inc)} residues have an incomplete side chain. Without -a Meeko crashes"
            f" (reported as excess inter-residue bond)"
        )
        out.append(
            (
                "W5",
                lvl,
                msg,
                inc[:12]
                + ([f"... and {len(inc) - 12} more"] if len(inc) > 12 else [])
                + [f"formal charge change if deleted by -a: {q_at_risk:+d} units"],
                "Complete the side chains before passing the file (running PDB2PQR fills them in). "
                "If deletion with -a is allowed, cross-check in the output which residues disappeared",
            )
        )
    return out, res


PQR_BASE = {
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


def pqr_states(path):
    """The state of each titratable residue in the pqr. Decided by the presence of hydrogens, not by name.

    PDB2PQR leaves the name as LYS even when it applies the LYN parameters to a neutral lysine
    (when --ffout is not given). Reading by name shows no neutral lysine at all."""
    res = parse_pqr(path)
    sg = {k: v["sg"] for k, v in res.items() if PQR_BASE.get(v["rn"], v["rn"]) == "CYS" and v["sg"]}
    ss = set()
    ks = list(sg)
    for i in range(len(ks)):
        for j in range(i + 1, len(ks)):
            if dist(sg[ks[i]], sg[ks[j]]) < 2.5:
                ss.add(ks[i])
                ss.add(ks[j])
    out = {}
    for k, v in res.items():
        a, b = v["at"], PQR_BASE.get(v["rn"], v["rn"])
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
            continue
        out[k] = (b, s)
    return out


def check_against(a_path, b_path):
    """C3: residues whose titration state differs between two pqr files.

    The pqr records no provenance. Whether the titration was run with the ligand in or with it out cannot
    be decided from a single file (the only trace is a single hydrogen on the changed residue, and it is
    indistinguishable from the "correct state"). The only way is to match the two against each other."""
    A, B = pqr_states(a_path), pqr_states(b_path)
    common = set(A) & set(B)
    diff = sorted(k for k in common if A[k][1] != B[k][1])
    only = set(A) ^ set(B)
    det = [f"{A[k][0]}{k[1]}{('(' + k[0] + ')') if k[0] else ''} {A[k][1]} -> {B[k][1]}" for k in diff[:20]]
    if len(diff) > 20:
        det.append(f"... and {len(diff) - 20} more")
    if only:
        det.append(f"{len(only)} titratable residues present in only one of them")
    if not diff and not only:
        return [
            ("C3", "INFO", "the titration states agree completely in the two pqr",
             [f"{len(common)} titratable residues compared"], "")
        ]
    return [
        (
            "C3",
            "WARN",
            f"{len(diff)} residues differ in titration state between the two pqr",
            [f"{len(common)} titratable residues compared"] + det,
            "Which one is used changes the receptor charge. Check whether a residue at the binding site is included",
        )
    ]


def parse_pqr(path):
    """Read the pqr. The key is chain, residue number and residue name; atom names and SG coordinates are kept."""
    res = collections.OrderedDict()
    for l in open(path, errors="replace"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        f = l.split()
        if len(f) < 8:
            continue
        rn = f[3]
        chain, num = "", ""
        for j in range(4, min(7, len(f))):
            if f[j].rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ").lstrip("-").isdigit():
                num = f[j]
                if j > 4:
                    chain = f[j - 1]  # the chain ID sits between residue name and number
                break
        if not num:
            continue
        # Include the chain in the key. Without it, residues with the same number on different chains merge,
        # and HD1 from one and HE2 from the other look like one residue, giving a fictitious "+1 HIS".
        e = res.setdefault((chain, num, rn), {"rn": rn, "at": set(), "sg": None})
        e["at"].add(f[2])
        if f[2] == "SG":
            try:
                e["sg"] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            except ValueError:
                pass
    return res


def check_pqr(path):
    """C1: whether the residue names express the state. The state is read from the number of hydrogens."""
    res = parse_pqr(path)
    hidden = collections.Counter()
    for (_ch, _num, rn), e in res.items():
        at = e["at"]
        if rn in ("ASP", "ASH") and "HD2" in at and rn == "ASP":
            hidden["neutral ASP still named ASP"] += 1
        if rn in ("GLU", "GLH") and "HE2" in at and rn == "GLU":
            hidden["neutral GLU still named GLU"] += 1
        if rn == "LYS" and sum(1 for x in at if x.startswith("HZ")) == 2:
            hidden["neutral LYS still named LYS"] += 1
        if rn == "HIS" and "HD1" in at and "HE2" in at:
            hidden["+1 HIS still named HIS"] += 1
        if rn == "HIS" and (("HD1" in at) != ("HE2" in at)):
            hidden["tautomer still named HIS"] += 1
        if rn == "CYS" and "HG" not in at:
            hidden["SS bond or thiolate still named CYS"] += 1
    if hidden:
        return [
            (
                "C1",
                "WARN",
                "the residue names in the pqr do not express the state. --ffout=AMBER was forgotten",
                [f"{k}: {v} residues" for k, v in hidden.most_common()],
                "Add --ffout=AMBER again. A downstream tool that reads the state from the name gets it wrong",
            )
        ]
    return []


def check_pdbqt(path, src_res):
    """C2: residues present in the input but absent from the pdbqt, and their formal charge."""
    # Do not include the residue name in the key. The preparation tools rename HIS to HIE / HID / HIP
    # and ASP to ASH, so matching with the name included would report every renamed residue
    # as "deleted".
    got = {(l[21], l[22:27].strip()) for l in open(path, errors="replace") if l.startswith(("ATOM", "HETATM"))}
    gone = [k for k in src_res if (k[0], k[1]) not in got]
    if not gone:
        return []
    q = sum(FORMAL.get(k[2], 0) for k in gone)
    ch = sum(1 for k in gone if k[2] in FORMAL)
    return [
        (
            "C2",
            "WARN",
            f"{len(gone)} residues are present in the input but absent from the pdbqt",
            [f"{k[2]}{k[1]}({k[0]})" for k in gone[:12]]
            + ([f"... and {len(gone) - 12} more"] if len(gone) > 12 else [])
            + [f"{ch} of them charged residues, formal charge change {q:+d} units"],
            "Check whether a deleted residue overlaps the binding site",
        )
    ]


def main():
    # The hand-written argv scan is dropped and argparse does the work. Previously
    #   args = [a for a in sys.argv[1:] if not a.startswith("--")]
    # removed only the flags and left "their values", so writing
    #   check_receptor.py --prepared out.pdbqt in.pdb
    # made out.pdbqt the file under inspection, and an input that should have been an ERROR
    # passed with 0 ERROR and exit 0. Putting --altloc last also crashed with IndexError.
    import argparse

    ap = argparse.ArgumentParser(
        description="Inspect the input before docking receptor preparation and cross-check the output afterwards",
        epilog="exit code: 0=no blocker  1=ERROR present  2=usage or input/output error",
    )
    ap.add_argument("input", help="PDB file of the receptor")
    ap.add_argument("--altloc", metavar="X", help="the altLoc that is going to be passed to --default_altloc")
    ap.add_argument(
        "--allow-bad-res",
        action="store_true",
        help="when -a is going to be given to Meeko. Drops incomplete side chains from ERROR to WARN",
    )
    ap.add_argument("--prepared", metavar="FILE", help="output to cross-check (.pqr or .pdbqt)")
    ap.add_argument(
        "--against",
        metavar="FILE",
        help="a second pqr to match against the pqr of --prepared. "
        "Run PDB2PQR on both the input with HETATM and the input without it, "
        "and report the residues whose titration state moved",
    )
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()

    src, as_json = a.input, a.json
    for p in (src, a.prepared, a.against):
        if p and not os.path.isfile(p):
            print(f"no such file: {p}", file=sys.stderr)
            return 2
    altlocs = (a.altloc,) if a.altloc else ("A", "B")
    try:
        findings, res = check_input(src, altlocs, a.allow_bad_res)
    except Exception as e:
        print(f"failed to read {src}: {type(e).__name__}: {e}", file=sys.stderr)
        return 2
    prep = a.prepared
    if prep:
        if prep.endswith(".pqr"):
            findings += check_pqr(prep)
        elif prep.endswith(".pdbqt"):
            findings += check_pdbqt(prep, set(res))
        else:
            print(f"--prepared only accepts .pqr or .pdbqt: {prep}", file=sys.stderr)
            return 2
    if a.against:
        if not prep or not prep.endswith(".pqr"):
            print("--against can only be used when a .pqr is given to --prepared", file=sys.stderr)
            return 2
        if not a.against.endswith(".pqr"):
            print(f"--against only accepts .pqr: {a.against}", file=sys.stderr)
            return 2
        findings += check_against(prep, a.against)
    nerr = sum(1 for f in findings if f[1] == "ERROR")
    nwarn = sum(1 for f in findings if f[1] == "WARN")
    if as_json:
        print(
            json.dumps(
                {
                    "input": src,
                    "prepared": prep,
                    "against": a.against,
                    "n_error": nerr,
                    "n_warn": nwarn,
                    "findings": [
                        {"id": fid, "level": lvl, "message": msg, "detail": det, "action": act}
                        for fid, lvl, msg, det, act in findings
                    ],
                },
                ensure_ascii=False,
                indent=1,
            )
        )
        return 1 if nerr else 0
    print(f"=== {os.path.basename(src)} ===")
    print(f"{nerr} ERROR, {nwarn} WARN\n")
    for fid, lvl, msg, det, act in findings:
        print(f"[{lvl:<5} {fid}] {msg}")
        for d in det:
            print(f"        {d}")
        if act:
            print(f"        -> {act}")
        print()
    if not findings:
        print("None of the rules checked here apply.")
    return 1 if nerr else 0


if __name__ == "__main__":
    sys.exit(main())
