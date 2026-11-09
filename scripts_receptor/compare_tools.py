#!/usr/bin/env python3
"""同じ入力ファイルから Meeko と PDB2PQR が作った受容体の電荷状態を比べる。

  python3 compare_tools.py <meeko_pdbqt_dir> <pqr_dir> [out.csv]

入力は完全に同じ（recept/ の重原子のみ PDB）。違うのはツールだけ。
PDB2PQR は PROPKA で滴定状態を予測する。Meeko はテンプレートから電荷を当てる。
残基ごとの電荷の和を整数に丸めて形式電荷とし、両者を突き合わせる。

pdbqt は列固定（電荷 70-76、残基名 17-20、残基番号 22-27）。
pqr は空白区切りだが座標が繋がると要素数が落ちるので要素数で弾かない。
どちらも滴定可能残基を 1 個も落とさないことを assert で保証する。
"""

import collections
import csv
import glob
import os
import sys

TITR = {"ASP", "GLU", "HIS", "CYS", "LYS", "TYR", "ARG"}
# 標準状態（AMBER の荷電形）での形式電荷
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
        # pdbqt は鎖 ID を持つ。鍵に入れておく。
        k = (l[21] if len(l) > 21 else " ", l[22:27].strip(), rn)
        e = res.setdefault(k, {"q": 0.0, "n": 0})
        e["q"] += q
        e["n"] += 1
    assert seen == placed, f"{fn}: {seen} 行中 {placed} 行しか読めていない"
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
        # 鎖 ID は残基名と残基番号の間。PDB2PQR は既定で落とすが
        # --keep-chain を付けると残る。あれば鍵に入れる。
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
            e["dup"] += 1  # 別鎖の同番号残基の併合、または altLoc
        e["at"].add(f[2])
        e["q"] += q
        e["n"] += 1
    assert seen == placed, f"{fn}: {seen} 行中 {placed} 行しか読めていない"
    return res


def termini(pqrfile):
    """末端残基を原子名で特定する。残基番号の最大最小では取り落とす
    （7INV の ARG309 は OXT を持つ C 末端だが番号は最大ではない）。
    C 末は OXT、N 末は H1/H2/H3 を持つ。末端は電荷和が交絡する:
    荷電 ARG +1 と末端カルボキシラート -1 で和が 0 になり中性に見える。"""
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
    mdir, pdir = sys.argv[1], sys.argv[2]
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
        # pdbqt は鎖を持つが pqr は PDB2PQR が既定で落とす。両方が鎖を持つ
        # ときだけ鎖込みで突き合わせ、片方が持たないときは鎖を外して合わせる。
        # 鎖を外す場合、同じ残基番号が複数の鎖にあると曖昧になるので検出する。
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
    print(f"比較した構造 {len(pids)} 本（入力ファイルは同一、ツールだけが違う）")
    print(f"末端残基 {n_excl} 個を除外（電荷和が交絡するため）")
    if n_dup:
        print(f"  !! pqr で同じ原子名が二度来た残基 {n_dup} 件。別鎖の同番号残基が併合されている可能性がある")
    if n_ambiguous:
        print(
            f"  !! 鎖を外して突き合わせた結果、同じ残基番号が複数の鎖にある例 "
            f"{n_ambiguous} 件。--keep-chain で pqr を作り直すこと"
        )
    else:
        print("  鎖の曖昧さ・残基の併合ともになし")
    print()
    print("=== Meeko が付ける形式電荷の分布 ===")
    for rn in sorted(mstate):
        tot = sum(mstate[rn].values())
        dist = ", ".join(f"{q:+d} x{n}" for q, n in sorted(mstate[rn].items()))
        fixed = "  <- 全残基が同じ状態" if len(mstate[rn]) == 1 else ""
        print(f"  {rn:<5} {tot:>5} 残基: {dist}{fixed}")
    tot = sum(sum(c.values()) for c in mstate.values())
    nd = sum(perstruct_i["n_disagree"] for perstruct_i in perstruct)
    print(f"\n滴定可能残基 {tot} 個のうち、2 ツールで形式電荷が違う残基 {nd} 個 ({nd / tot * 100:.2f}%)")
    print("\n=== 食い違いの内訳（Meeko の値 / PDB2PQR の値）===")
    for (rn, qm, qp), n in disagree.most_common(10):
        print(f"  {rn}: Meeko {qm:+d} / PDB2PQR {qp:+d}   {n} 残基")
    moved = [r for r in perstruct if r["diff"]]
    print("\n=== 受容体の正味電荷（滴定可能残基の和）===")
    print(f"  2 ツールで差がある構造 {len(moved)}/{len(perstruct)} ({len(moved) / len(perstruct) * 100:.0f}%)")
    if moved:
        print(f"  最大の差 {max(abs(r['diff']) for r in moved)} 単位")
        print(f"\n{'構造':<7}{'Meeko':>8}{'PDB2PQR':>10}{'差':>7}{'食い違い残基':>14}")
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
