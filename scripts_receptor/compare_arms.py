#!/usr/bin/env python3
"""HETATM の有無で滴定状態がどう変わるかを残基ごとに比べる。

  python3 compare_arms.py <pqr_hetall> <pqr_nowat> <pqr_none> [out.csv]

3 条件はタンパクの ATOM 記録が完全に同じで、HETATM の有無だけが違う。
水素は 3 条件すべてで落としてあるので、寄託水素の影響は入らない。

  hetall : 重原子 + HETATM 全部（配位子 + 水）
  nowat  : 重原子 + 水以外の HETATM（水だけ除去）
  none   : 重原子のみ（水も配位子も除去）

hetall -> nowat の差が水の除去、nowat -> none の差が配位子の除去による変化。
状態の判定は classify_titr.py（残基名ではなく水素の有無）。
"""

import collections
import csv
import glob
import importlib.util
import os
import sys

spec = importlib.util.spec_from_file_location(
    "ct", os.path.join(os.path.dirname(os.path.abspath(__file__)), "classify_titr.py")
)
ct = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ct)


def states(pqr):
    res, seen, placed = ct.parse(pqr)
    assert seen == placed, f"{pqr}: 滴定可能残基の行 {seen} 中 {placed} しか割り当てていない"
    # 同じ原子名が二度来た残基は、別鎖の同番号残基が併合されたか altLoc。
    # 黙って通すと片方の HD1 と他方の HE2 で架空の HIP(+1) が生まれる。
    dup = sum(1 for v in res.values() if v.get("dup"))
    return {k: v for k, v in ct.classify(res).items()}, dup


def main():
    dirs = sys.argv[1:4]
    out = sys.argv[4] if len(sys.argv) > 4 else None
    labels = ["hetall", "nowat", "none"]
    pids = sorted(
        set.intersection(
            *[{os.path.basename(f)[:4].upper() for f in glob.glob(os.path.join(d, "*.pqr"))} for d in dirs]
        )
    )
    rows, keymismatch, dup_struct = [], 0, []
    for pid in pids:
        got3 = [states(os.path.join(d, pid + ".pqr")) for d in dirs]
        S = [g[0] for g in got3]
        if any(g[1] for g in got3):
            dup_struct.append((pid, max(g[1] for g in got3)))
        keys = set(S[0]) | set(S[1]) | set(S[2])
        if not (set(S[0]) == set(S[1]) == set(S[2])):
            keymismatch += 1
        for k in sorted(keys):
            got = [s.get(k) for s in S]
            if any(g is None for g in got):
                continue
            base = got[0][0]
            st = [g[2] for g in got]
            if st[0] == st[1] == st[2]:
                continue
            # 鍵は (鎖, 残基番号, 残基名)。PDB2PQR は既定で鎖 ID を落とすので
            # chain は空になることが多い。空のまま残し、後段が鎖の有無を知れるようにする。
            rows.append(
                {
                    "pdb_id": pid,
                    "chain": k[0],
                    "resnum": k[1],
                    "base": base,
                    "state_hetall": st[0],
                    "state_nowat": st[1],
                    "state_none": st[2],
                    "water_removal_changed": int(st[0] != st[1]),
                    "ligand_removal_changed": int(st[1] != st[2]),
                }
            )
    print(
        f"比較した構造 {len(pids)} 本"
        + (f"（残基キーが一致しない構造 {keymismatch} 本）" if keymismatch else "（残基キーは全構造で一致）")
    )
    nochain = sum(1 for r in rows if not r["chain"])
    if nochain:
        print(f"  注意: pqr に鎖 ID が無い行 {nochain}/{len(rows)}。PDB2PQR は既定で鎖を落とす（--keep-chain で残る）")
    if dup_struct:
        print(f"  !! 同じ原子名が二度来た残基を持つ構造 {len(dup_struct)} 本: {dup_struct[:5]}")
        print("     別鎖の同番号残基が併合されている可能性がある。--keep-chain を付けて作り直すこと")
    else:
        print("  原子名の重複なし。残基の併合は起きていない")
    tot = sum(len(states(os.path.join(dirs[2], p + ".pqr"))[0]) for p in pids)
    print(f"滴定可能残基 {tot} 個のうち、3 条件で状態が揃わない残基 {len(rows)} 個 ({len(rows) / tot * 100:.2f}%)")
    if out and rows:
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    print("\n=== 変化の内訳（残基種ごと）===")
    print(f"{'素の残基':<9}{'水の除去で変化':>16}{'配位子の除去で変化':>20}{'合計':>7}")
    byb = collections.defaultdict(lambda: [0, 0, 0])
    for r in rows:
        b = byb[r["base"]]
        b[0] += r["water_removal_changed"]
        b[1] += r["ligand_removal_changed"]
        b[2] += 1
    for b, (w, l, t) in sorted(byb.items(), key=lambda kv: -kv[1][2]):
        print(f"{b:<9}{w:>16}{l:>20}{t:>7}")
    w_all = sum(r["water_removal_changed"] for r in rows)
    l_all = sum(r["ligand_removal_changed"] for r in rows)
    print(f"{'計':<9}{w_all:>16}{l_all:>20}{len(rows):>7}")
    nstr = len({r["pdb_id"] for r in rows})
    print(f"\n状態が変わる残基を 1 個以上持つ構造 {nstr}/{len(pids)} = {nstr / len(pids) * 100:.0f}%")

    print("\n=== どの状態からどの状態へ（上位）===")
    tr = collections.Counter()
    for r in rows:
        if r["state_hetall"] != r["state_none"]:
            tr[(r["base"], r["state_hetall"], r["state_none"])] += 1
    for (b, a, z), n in tr.most_common(12):
        print(f"  {b}: {a} -> {z}   {n} 残基")
    if out:
        print(f"\n-> {out}")


if __name__ == "__main__":
    main()
