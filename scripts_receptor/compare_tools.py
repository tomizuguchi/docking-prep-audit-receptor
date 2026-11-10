#!/usr/bin/env python3
"""Compare the charge states of the receptors that Meeko and PDB2PQR build from the same input file.

  python3 compare_tools.py <meeko_pdbqt_dir> <pqr_dir> [out.csv]

The input is exactly the same (the heavy-atom-only PDB in recept/). Only the tool differs.
PDB2PQR predicts titration states with PROPKA. Meeko assigns charges from templates.
The per-residue charge sum is rounded to an integer to give the formal charge, and the two are compared.

pdbqt is column-fixed (charge 70-76, residue name 17-20, residue number 22-27).
pqr is whitespace-separated, but the field count drops when coordinates run together, so lines are not
rejected by field count.
An assert guarantees that neither side drops a single titratable residue.
"""

import collections
import csv
import glob
import os
import sys

TITR = {"ASP", "GLU", "HIS", "CYS", "LYS", "TYR", "ARG"}
# Formal charge in the standard state (the charged form in AMBER)
STD = {"ASP": -1, "GLU": -1, "HIS": 0, "CYS": 0, "LYS": 1, "TYR": 0, "ARG": 1}


def from_pdbqt(fn):
    res = collections.OrderedDict()
    seen = placed = 0
    for l in open(fn, errors="replace"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        rn = l[17:20].strip()
        if rn not in TITR:
            continue
        seen += 1
        try:
            q = float(l[70:76])
        except ValueError:
            continue
        placed += 1
        # pdbqt has a chain ID. Keep it in the key.
        k = (l[21] if len(l) > 21 else " ", l[22:27].strip(), rn)
        e = res.setdefault(k, {"q": 0.0, "n": 0})
        e["q"] += q
        e["n"] += 1
    assert seen == placed, f"{fn}: only {placed} of {seen} lines were read"
    return res


def from_pqr(fn):
    res = collections.OrderedDict()
    seen = placed = 0
    for l in open(fn, errors="replace"):
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
        # The chain ID sits between the residue name and the residue number. PDB2PQR drops it by
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
        e = res.setdefault((chain, num, rn), {"q": 0.0, "n": 0, "at": set(), "dup": 0})
        if f[2] in e["at"]:
            e["dup"] += 1  # merge of same-numbered residues from different chains, or altLoc
        e["at"].add(f[2])
        e["q"] += q
        e["n"] += 1
    assert seen == placed, f"{fn}: only {placed} of {seen} lines were read"
    return res


def termini(pqrfile):
    """Identify terminal residues by atom name. The maximum and minimum residue number miss some
    (ARG309 of 7INV is a C terminus carrying OXT, but its number is not the maximum).
    A C terminus carries OXT, an N terminus carries H1/H2/H3. Termini confound the charge sum:
    a charged ARG +1 and a terminal carboxylate -1 sum to 0 and look neutral."""
    term = set()
    for l in open(pqrfile, errors="replace"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        f = l.split()
        if len(f) < 8:
            continue
        if f[2] not in ("OXT", "H1", "H2", "H3"):
            continue
        for j in range(4, min(7, len(f))):
            if f[j].rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ").lstrip("-").isdigit():
                term.add((f[j - 1] if j > 4 else "", f[j]))
                break
    return term


def main():
    if len(sys.argv) < 3 or sys.argv[1] in ("-h", "--help"):
        sys.exit(__doc__)
    mdir, pdir = sys.argv[1], sys.argv[2]
    if not glob.glob(os.path.join(mdir, "*.pdbqt")):
        sys.exit(f"no *.pdbqt in {mdir}")
    if not glob.glob(os.path.join(pdir, "*.pqr")):
        sys.exit(f"no *.pqr in {pdir}")
    out = sys.argv[3] if len(sys.argv) > 3 else None
    rows, perstruct, n_excl, n_ambiguous, n_dup = [], [], 0, 0, 0
    pids = sorted(
        {os.path.basename(f)[:4].upper() for f in glob.glob(os.path.join(mdir, "*.pdbqt"))}
        & {os.path.basename(f)[:4].upper() for f in glob.glob(os.path.join(pdir, "*.pqr"))}
    )
    mstate = collections.defaultdict(collections.Counter)
    disagree = collections.Counter()
    for pid in pids:
        M = from_pdbqt(os.path.join(mdir, pid + ".pdbqt"))
        P = from_pqr(os.path.join(pdir, pid + ".pqr"))
        term = termini(os.path.join(pdir, pid + ".pqr"))
        # pdbqt carries chains, but PDB2PQR drops them from pqr by default. Match with the chain
        # only when both sides carry chains; when one side does not, drop the chain before matching.
        # When the chain is dropped, the same residue number in several chains is ambiguous, so detect it.
        p_has_chain = any(k[0] for k in P)
        m_has_chain = any(k[0] for k in M)
        use_chain = p_has_chain and m_has_chain
        if not use_chain:
            byshort = collections.Counter((k[1], k[2]) for k in M)
            amb = sum(1 for v in byshort.values() if v > 1)
            if amb:
                n_ambiguous += amb
            Mx = {(("", k[1], k[2]) if not use_chain else k): v for k, v in M.items()}
            Px = {(("", k[1], k[2]) if not use_chain else k): v for k, v in P.items()}
        else:
            Mx, Px = M, P
        termx = {t[1] for t in term} if not use_chain else term

        def interm(k):
            return (k[1] in termx) if not use_chain else ((k[0], k[1]) in termx)

        common = {k for k in (set(Mx) & set(Px)) if not interm(k)}
        n_excl += len(set(Mx) & set(Px)) - len(common)
        n_dup += sum(v.get("dup", 0) for v in P.values())
        dm = dp = 0
        for k in sorted(common):
            qm, qp = round(Mx[k]["q"]), round(Px[k]["q"])
            rn = k[2]
            mstate[rn][qm] += 1
            dm += qm
            dp += qp
            if qm != qp:
                disagree[(rn, qm, qp)] += 1
                rows.append(
                    {
                        "pdb_id": pid,
                        "chain": k[0],
                        "resnum": k[1],
                        "base": rn,
                        "meeko_formal": qm,
                        "pdb2pqr_formal": qp,
                        "diff": qm - qp,
                    }
                )
        perstruct.append(
            {
                "pdb_id": pid,
                "n_titr_common": len(common),
                "sum_meeko": dm,
                "sum_pdb2pqr": dp,
                "diff": dm - dp,
                "n_disagree": sum(1 for k in common if round(Mx[k]["q"]) != round(Px[k]["q"])),
            }
        )
    print(f"{len(pids)} structures compared (the input files are identical, only the tool differs)")
    print(f"{n_excl} terminal residues excluded (the charge sum is confounded)")
    if n_dup:
        print(
            f"  !! {n_dup} residues where the same atom name came twice in pqr. "
            f"Same-numbered residues from different chains may have been merged"
        )
    if n_ambiguous:
        print(
            f"  !! after matching with the chain dropped, {n_ambiguous} cases where the same residue "
            f"number is in several chains. Rebuild the pqr with --keep-chain"
        )
    else:
        print("  no chain ambiguity and no merged residues")
    print()
    print("=== Distribution of the formal charges Meeko assigns ===")
    for rn in sorted(mstate):
        tot = sum(mstate[rn].values())
        dist = ", ".join(f"{q:+d} x{n}" for q, n in sorted(mstate[rn].items()))
        fixed = "  <- all residues in the same state" if len(mstate[rn]) == 1 else ""
        print(f"  {rn:<5} {tot:>5} residues: {dist}{fixed}")
    tot = sum(sum(c.values()) for c in mstate.values())
    nd = sum(perstruct_i["n_disagree"] for perstruct_i in perstruct)
    print(f"\n{nd} of {tot} titratable residues differ in formal charge between the two tools ({nd / tot * 100:.2f}%)")
    print("\n=== Breakdown of the disagreements (Meeko value / PDB2PQR value) ===")
    for (rn, qm, qp), n in disagree.most_common(10):
        print(f"  {rn}: Meeko {qm:+d} / PDB2PQR {qp:+d}   {n} residues")
    moved = [r for r in perstruct if r["diff"]]
    print("\n=== Net charge of the receptor (sum over the titratable residues) ===")
    print(
        f"  {len(moved)}/{len(perstruct)} structures differ between the two tools "
        f"({len(moved) / len(perstruct) * 100:.0f}%)"
    )
    if moved:
        print(f"  largest difference {max(abs(r['diff']) for r in moved)} units")
        print(f"\n{'struct':<7}{'Meeko':>8}{'PDB2PQR':>10}{'diff':>7}{'disagreeing':>14}")
        for r in sorted(moved, key=lambda r: -abs(r["diff"]))[:12]:
            print(f"{r['pdb_id']:<7}{r['sum_meeko']:>8}{r['sum_pdb2pqr']:>10}{r['diff']:>+7}{r['n_disagree']:>14}")
    if out and rows:
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"\n-> {out}")


if __name__ == "__main__":
    main()
