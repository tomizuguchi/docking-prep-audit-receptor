#!/usr/bin/env python3
"""
Profiler for PDB structures. Enumerates the code paths of the preprocessing pipeline.

  python3 profile_pdb.py prof  <pdbdir> [out.csv]   tabulate the flags of each structure
  python3 profile_pdb.py pick  <out.csv> [N]        choose the smallest set that covers all flags

Residue numbers are handled as the pair (resSeq, iCode). They are not collapsed with int().
"""

import collections
import csv
import glob
import os
import sys

AA = {
    "ALA": 5,
    "ARG": 11,
    "ASN": 8,
    "ASP": 8,
    "CYS": 6,
    "GLN": 9,
    "GLU": 9,
    "GLY": 4,
    "HIS": 10,
    "ILE": 8,
    "LEU": 8,
    "LYS": 9,
    "MET": 8,
    "PHE": 11,
    "PRO": 7,
    "SER": 6,
    "THR": 7,
    "TRP": 14,
    "TYR": 12,
    "VAL": 7,
}
NONSTD = {"MSE", "SEP", "TPO", "PTR", "CSO", "CME", "KCX", "MLY", "ALY", "HYP", "PCA", "LLP", "SAH", "SAM"}
TITRATABLE = {"LYS", "ARG", "ASP", "GLU", "HIS", "CYS", "TYR"}


def parse(path):
    r = {"file": os.path.basename(path), "id": os.path.basename(path).split(".")[0][-4:].upper()}
    method = resol = rwork = None
    res = collections.OrderedDict()  # (chain, resSeq, iCode) -> [resName, {atomName: altLoc}]
    het = collections.Counter()
    altlocs = collections.Counter()
    # Keep the set of altLocs for each atom name. The {atomName: altLoc} in res holds only one value,
    # so it cannot express "no A and no blank, only B" (the CA of 4MLK MET59).
    altsets = collections.defaultdict(set)
    ssbond = 0
    nH = 0
    nmodel = 0
    for l in open(path, errors="replace"):
        rec = l[:6]
        if rec == "EXPDTA":
            method = l[10:].strip()[:40]
        elif l.startswith("REMARK   2 RESOLUTION"):
            # Take the 1.30 out of "REMARK   2 RESOLUTION.    1.30 ANGSTROMS.".
            # split() on the whole line picks up the 2 of the REMARK number first,
            # so only what follows "RESOLUTION." is examined. A NOT APPLICABLE line has no value.
            for tok in l.split("RESOLUTION.")[-1].split():
                try:
                    v = float(tok)
                    if 0.3 < v < 30:
                        resol = v
                        break
                except ValueError:
                    pass
        elif l.startswith("REMARK   3   R VALUE            (WORKING SET)"):
            try:
                rwork = float(l.split(":")[1])
            except Exception:
                pass
        elif rec == "SSBOND":
            ssbond += 1
        elif rec == "MODEL ":
            nmodel += 1
        elif rec == "ATOM  " or rec == "HETATM":
            el = l[76:78].strip() or l[12:16].strip()[0]
            if el == "H":
                nH += 1
            alt = l[16]
            if alt != " ":
                altlocs[alt] += 1
            if rec == "HETATM":
                het[l[17:20].strip()] += 1
            else:
                k = (l[21], int(l[22:26]), l[26])
                res.setdefault(k, [l[17:20].strip(), {}])[1][l[12:16].strip()] = alt
                # Put only heavy atoms into the altLoc set. The receptor file passed to Meeko has its
                # hydrogens removed, so counting the altLocs of deposited hydrogens gives a false positive
                # (1EA7 has 2036 hydrogens with an altLoc, but Meeko passes with both A and B).
                if el != "H":
                    altsets[(k, l[12:16].strip())].add(alt)
    r["method"] = method or ""
    r["xray"] = 1 if method and "X-RAY" in method.upper() else 0
    r["resol"] = resol if resol else ""
    r["rwork"] = rwork if rwork else ""
    chains = sorted({k[0] for k in res})
    r["n_chain"] = len(chains)
    r["chains"] = "".join(chains)[:12]
    r["n_res"] = len(res)
    r["n_H"] = nH
    r["n_model"] = nmodel
    # --- per-chain, using the first chain with >=50 residues
    tgt = next((c for c in chains if sum(1 for k in res if k[0] == c) >= 50), chains[0] if chains else None)
    r["target_chain"] = tgt or ""
    sub = {k: v for k, v in res.items() if k[0] == tgt}
    nums = sorted({k[1] for k in sub})
    r["chain_res"] = len(sub)
    r["gaps"] = sum(1 for a, b in zip(nums, nums[1:]) if b != a + 1)
    r["icode"] = sum(1 for k in sub if k[2] != " ")
    r["neg_resnum"] = sum(1 for k in sub if k[1] <= 0)
    r["altloc_atoms"] = sum(altlocs.values())
    r["altloc_only_B"] = sum(1 for k, v in sub.items() if v[1] and all(a not in (" ", "A") for a in v[1].values()))
    inc = [(k, v[0], len(v[1]), AA[v[0]]) for k, v in sub.items() if v[0] in AA and len(v[1]) < AA[v[0]]]
    r["incomplete_sc"] = len(inc)
    r["incomplete_atoms"] = sum(e - n for _, _, n, e in inc)
    r["incomplete_titr"] = sum(1 for _, rn, _, _ in inc if rn in TITRATABLE)
    # A modified residue can arrive on an ATOM line or on a HETATM line. MSE and the like are
    # normally HETATM, so looking only at ATOM lines drops them (4/40 missed in the pilot).
    # Missing main-chain atoms. What makes PDB2PQR fail is a missing main chain, not a missing side chain
    # (pilot: the 2 that failed both had a residue missing C/CA/O at once, and
    #  with only 1 atom missing, repair_heavy rebuilt it and it passed).
    BB = {"N", "CA", "C", "O"}
    bbmiss = [(k, v[0], BB - set(v[1])) for k, v in sub.items() if BB - set(v[1])]
    r["bb_missing_res"] = len(bbmiss)
    r["bb_max_missing"] = max((len(m) for _, _, m in bbmiss), default=0)
    # Residues missing C and O at once. The PDB2PQR "Too few atoms present to
    # reconstruct or cap residue" appears under this condition (3/3 in 100 entries, 0 false positives).
    r["bb_cap"] = sum(1 for _, _, m in bbmiss if {"C", "O"} <= m)
    r["bb_missing_detail"] = ";".join(f"{rn}{k[1]}:{'+'.join(sorted(m))}" for k, rn, m in bbmiss[:4])

    # Per-atom altLoc: the case where an atom exists only under an altLoc that is neither A nor blank.
    # The per-residue F_altloc_onlyB cannot catch it (in 4MLK MET59, N/C/O are blank and
    # only CA has altLoc B). Meeko's --default_altloc A cannot find this atom and fails.
    # An atom with neither A nor blank present, existing only under another altLoc. --default_altloc A does not find it.
    # Atoms that have neither X nor blank when the default altLoc is X. Meeko
    # cannot find this atom with --default_altloc X and fails (4MLK with A, 3I6H/7USK with B).
    def no_variant(X):
        return [(k, a) for (k, a), st in altsets.items() if k[0] == tgt and not (st & {" ", X})]

    for X in ("A", "B"):
        nv = no_variant(X)
        r[f"atoms_no_{X}"] = len(nv)
        r[f"atoms_no_{X}_bb"] = sum(1 for _, a in nv if a in ("N", "CA", "C", "O"))
        r[f"atoms_no_{X}_detail"] = ";".join(f"{k[1]}:{a}" for k, a in nv[:4])
    r["nonstd"] = sum(1 for k, v in sub.items() if v[0] in NONSTD)
    r["nonstd_het"] = sum(1 for n in het if n in NONSTD)
    r["nonstd_names"] = ",".join(
        sorted(set([v[0] for v in sub.values() if v[0] in NONSTD] + [n for n in het if n in NONSTD]))
    )
    r["ssbond"] = ssbond
    r["n_titr"] = sum(1 for v in sub.values() if v[0] in TITRATABLE)
    r["n_lys"] = sum(1 for v in sub.values() if v[0] == "LYS")
    r["het_kinds"] = len(het)
    r["het_top"] = ",".join(f"{n}" for n, _ in het.most_common(4))
    r["has_nuc"] = 1 if any(v[0] in ("DA", "DT", "DG", "DC", "A", "U", "G", "C") for v in sub.values()) else 0
    return r


# Flags that stay at zero however much the sample is increased, because the search conditions of
# fetch_strata.py remove them from the population. Distinguish these from a detection failure.
FRAME_EXCLUDED = {
    "F_multichain": "deposited_polymer_entity_instance_count == 1",
    "F_nucleic": 'selected_polymer_entity_types == "Protein (only)"',
    "F_no_resol": "range specified on resolution_combined",
}

FLAGS = {
    "F_clean_sc": lambda r: r["incomplete_sc"] == 0,
    "F_many_inc_sc": lambda r: r["incomplete_sc"] >= 20,
    "F_inc_titr": lambda r: r["incomplete_titr"] >= 5,
    "F_gaps": lambda r: r["gaps"] >= 1,
    "F_altloc": lambda r: r["altloc_atoms"] >= 1,
    "F_altloc_onlyB": lambda r: r["altloc_only_B"] >= 1,
    "F_icode": lambda r: r["icode"] >= 1,
    "F_neg_resnum": lambda r: r["neg_resnum"] >= 1,
    "F_nonstd": lambda r: r["nonstd"] >= 1 or r["nonstd_het"] >= 1,
    "F_bb_missing": lambda r: r["bb_missing_res"] >= 1,
    "F_bb_stub": lambda r: r["bb_max_missing"] >= 2,
    "F_atom_noA": lambda r: r["atoms_no_A"] >= 1,
    "F_atom_noA_bb": lambda r: r["atoms_no_A_bb"] >= 1,
    "F_atom_noB": lambda r: r["atoms_no_B"] >= 1,
    "F_bb_cap": lambda r: r["bb_cap"] >= 1,
    "F_ssbond": lambda r: r["ssbond"] >= 1,
    "F_het": lambda r: r["het_kinds"] >= 2,
    "F_has_H": lambda r: r["n_H"] >= 1,
    "F_multimodel": lambda r: r["n_model"] >= 2,
    "F_multichain": lambda r: r["n_chain"] >= 3,
    "F_nucleic": lambda r: r["has_nuc"] == 1,
    "F_hires": lambda r: r["resol"] != "" and float(r["resol"]) <= 1.5,
    "F_lores": lambda r: r["resol"] != "" and float(r["resol"]) >= 2.8,
    "F_big": lambda r: r["chain_res"] >= 500,
    "F_small": lambda r: r["chain_res"] <= 150,
    "F_no_resol": lambda r: r["resol"] == "",
}


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    mode = sys.argv[1]
    if mode == "prof":
        files = sorted(glob.glob(os.path.join(sys.argv[2], "*.pdb")) + glob.glob(os.path.join(sys.argv[2], "*.ent")))
        out = sys.argv[3] if len(sys.argv) > 3 else "profile.csv"
        rows = []
        for f in files:
            try:
                rows.append(parse(f))
            except Exception as e:
                print(f"  !! {os.path.basename(f)}: {type(e).__name__}: {e}")
        for r in rows:
            for name, fn in FLAGS.items():
                try:
                    r[name] = int(bool(fn(r)))
                except Exception:
                    r[name] = ""
        cols = [c for c in rows[0] if not c.startswith("F_")] + list(FLAGS)
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        print(f"{len(rows)} structures -> {out}\n")
        print("number of matching structures per flag")
        for name in FLAGS:
            n = sum(r.get(name) == 1 for r in rows)
            if n:
                mark = ""
            elif name in FRAME_EXCLUDED:
                mark = f"  <- no match. Cannot occur under the search conditions ({FRAME_EXCLUDED[name]})"
            else:
                mark = "  <- no match. Fetch additional structures with this condition"
            print(f"  {name:<16} {n:>3}{mark}")
    elif mode == "pick":
        rows = list(csv.DictReader(open(sys.argv[2])))
        N = int(sys.argv[3]) if len(sys.argv) > 3 else 15
        need = {f for f in FLAGS if any(r.get(f) == "1" for r in rows)}
        chosen, covered = [], set()
        while need - covered and len(chosen) < N:
            best = max(
                (r for r in rows if r not in chosen),
                key=lambda r: (len({f for f in need if r.get(f) == "1"} - covered), -float(r["n_res"] or 0)),
            )
            gain = {f for f in need if best.get(f) == "1"} - covered
            if not gain:
                break
            chosen.append(best)
            covered |= gain
            print(f"  {best['id']}  newly covered flags: {' '.join(sorted(gain))}")
        print(f"\npilot {len(chosen)} entries: " + " ".join(r["id"] for r in chosen))
        miss = need - covered
        print(
            f"flags covered {len(covered)}/{len(need)}"
            + (f", not reached: {' '.join(sorted(miss))}" if miss else " (all flags achieved)")
        )
    else:
        print(__doc__)


if __name__ == "__main__":
    main()
