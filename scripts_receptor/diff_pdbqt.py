#!/usr/bin/env python3
"""同じ構造から altLoc 設定だけを変えて作った 2 つの receptor pdbqt を比べる。

  python3 diff_pdbqt.py <dirA> <dirB> [out.csv]

落ちるか落ちないかではなく、両方通った構造で中身がどれだけ変わるかを測る。
比べるのは原子数、電荷の総和、残基の出入り、対応づけた原子の電荷差、原子型。

pdbqt の列は PDB 準拠 + 電荷(70-76) + AutoDock 原子型(77-79)。
列固定で読む。空白区切りは座標が繋がると壊れる（pqr で踏んだ）。
"""

import csv
import glob
import os
import sys


def read(fn):
    """(残基名, 残基番号, 原子名) -> (電荷, 原子型) の辞書を返す。"""
    d = {}
    for l in open(fn, errors="replace"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        key = (l[17:20].strip(), l[22:27].strip(), l[12:16].strip())
        try:
            q = float(l[70:76])
        except ValueError:
            q = float("nan")
        d[key] = (q, l[77:79].strip())
    return d


def main():
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
        resA = {(k[0], k[1]) for k in A}
        resB = {(k[0], k[1]) for k in B}
        rows.append(
            {
                "pdb_id": pid,
                "n_atom_A": len(A),
                "n_atom_B": len(B),
                "d_n_atom": len(A) - len(B),
                "sumq_A": round(sum(q for q, _ in A.values()), 4),
                "sumq_B": round(sum(q for q, _ in B.values()), 4),
                "d_sumq": round(sum(q for q, _ in A.values()) - sum(q for q, _ in B.values()), 4),
                "atoms_only_A": len(onlyA),
                "atoms_only_B": len(onlyB),
                "res_only_A": len(resA - resB),
                "res_only_B": len(resB - resA),
                "n_charge_moved": len(moved),
                "max_abs_dq": round(max((abs(v) for v in moved.values()), default=0.0), 4),
                "n_type_diff": len(tdiff),
                "identical": int(not onlyA and not onlyB and not moved and not tdiff),
            }
        )
    if not rows:
        print("比較できる構造がありません")
        return
    if out:
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    n = len(rows)
    same = sum(r["identical"] for r in rows)
    print(f"両設定で pdbqt が出た構造 {n} 本")
    print(f"  完全に同一 {same} 本、どこかが違う {n - same} 本")
    print(f"  原子数が違う {sum(1 for r in rows if r['d_n_atom'])} 本")
    print(f"  電荷の総和が違う {sum(1 for r in rows if abs(r['d_sumq']) > 1e-9)} 本")
    print(f"  原子型が違う原子を含む {sum(1 for r in rows if r['n_type_diff'])} 本")
    print(f"  残基の出入りがある {sum(1 for r in rows if r['res_only_A'] or r['res_only_B'])} 本")
    diff = [r for r in rows if not r["identical"]]
    if diff:
        print(
            f"\n{'構造':<7}{'原子 A':>7}{'原子 B':>7}{'Δ原子':>7}{'Σq A':>9}{'Σq B':>9}"
            f"{'ΔΣq':>9}{'電荷が動いた原子':>17}{'max|Δq|':>9}{'型違い':>7}"
        )
        for r in sorted(diff, key=lambda r: -abs(r["d_sumq"])):
            print(
                f"{r['pdb_id']:<7}{r['n_atom_A']:>7}{r['n_atom_B']:>7}{r['d_n_atom']:>+7}"
                f"{r['sumq_A']:>9.3f}{r['sumq_B']:>9.3f}{r['d_sumq']:>+9.3f}"
                f"{r['n_charge_moved']:>17}{r['max_abs_dq']:>9.3f}{r['n_type_diff']:>7}"
            )
    if out:
        print(f"\n→ {out}")


if __name__ == "__main__":
    main()
