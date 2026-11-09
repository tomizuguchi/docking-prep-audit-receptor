#!/usr/bin/env python3
"""pqr 群から中性リジン率を測る。取りこぼしゼロを検算しながら数える。

  python3 measure_lys.py <pqrdir> <manifest.csv>

パーサの教訓: pqr は空白区切りだが座標が繋がって要素数が落ちる
（6ZHF では ATOM 行の 67% が 9 要素）。要素数で行を弾くと黙って原子が消える。
この版は LYS/LYN の ATOM 行を 1 本も落とさないことを assert で保証する。

リジンの状態は NZ に付く水素の数で決める。
  HZ 3 本 = 荷電 LYS (+1)、HZ 2 本 = 中性 LYN、0/1 本 = 側鎖が不完全。
残基の電荷の和は補助にしか使わない。C 末リジンは荷電でも和が 0 になる
（末端カルボキシラートが打ち消す）ので、単独では誤判定する。
"""

import collections
import csv
import glob
import math
import os
import sys


def parse(pqr):
    """(残基ごとの状態, LYS/LYN の ATOM 行数, 割り当てられた行数) を返す。"""
    res, seen, placed = collections.OrderedDict(), 0, 0
    for l in open(pqr, errors="replace"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        f = l.split()
        if len(f) < 8:
            continue
        if f[3] not in ("LYS", "LYN"):
            continue
        seen += 1
        try:
            q = float(f[-2])
            float(f[-1])
        except ValueError:
            continue
        resnum = ""
        for t in f[4:7]:
            if t.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ").lstrip("-").isdigit():
                resnum = t
                break
        if not resnum:
            continue
        placed += 1
        e = res.setdefault((resnum, f[3]), {"name": f[3], "hz": 0, "q": 0.0, "n": 0})
        e["q"] += q
        e["n"] += 1
        if f[2].startswith("HZ"):
            e["hz"] += 1
    return list(res.values()), seen, placed


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)


def main():
    pqrdir, manfile = sys.argv[1], sys.argv[2]
    man = {r["pdb_id"]: r for r in csv.DictReader(open(manfile))}
    STR = ["A_le1.5", "B_1.5-2.0", "C_2.0-2.5", "D_2.5-3.0", "E_ge3.0"]
    pool = {s: int(next(r["pool_size"] for r in man.values() if r["stratum"] == s)) for s in STR}
    tot = sum(pool.values())

    per = {}
    for f in sorted(glob.glob(os.path.join(pqrdir, "*.pqr"))):
        pid = os.path.basename(f)[:4].upper()
        st, seen, placed = parse(f)
        assert seen == placed, f"{pid}: LYS/LYN の行 {seen} 本中 {placed} 本しか割り当てられていない"
        per[pid] = {
            "charged": sum(1 for s in st if s["hz"] == 3),
            "neutral": sum(1 for s in st if s["hz"] == 2),
            "broken": sum(1 for s in st if s["hz"] < 2),
            "total": len(st),
            "neutral_named_LYS": sum(1 for s in st if s["hz"] == 2 and s["name"] == "LYS"),
            "named_LYN": sum(1 for s in st if s["name"] == "LYN"),
            "atoms": sum(s["n"] for s in st),
        }
    print(f"取りこぼし検算: {len(per)} 本すべてで LYS/LYN の ATOM 行を全数割り当て済み\n")

    print("=== 層ごと ===")
    print(
        f"{'層':<11}{'重み':>7}{'本':>4}{'LYS 計':>8}{'荷電':>6}{'中性':>6}{'不完全':>7}"
        f"{'中性率/残基':>12}{'中性を含む構造':>14}"
    )
    wr_res = wr_str = 0.0
    agg = collections.Counter()
    nstr_all = 0
    for s in STR:
        ids = [p for p, r in man.items() if r["stratum"] == s and p in per]
        if not ids:
            continue
        t = collections.Counter()
        for p in ids:
            for k, v in per[p].items():
                t[k] += v
        nstr = sum(1 for p in ids if per[p]["neutral"] > 0)
        nstr_all += nstr
        w = pool[s] / tot
        pr = t["neutral"] / t["total"] if t["total"] else 0.0
        wr_res += w * pr
        wr_str += w * (nstr / len(ids))
        for k, v in t.items():
            agg[k] += v
        print(
            f"{s:<11}{w:>7.3f}{len(ids):>4}{t['total']:>8}{t['charged']:>6}"
            f"{t['neutral']:>6}{t['broken']:>7}{pr * 100:>11.2f}%{f'{nstr}/{len(ids)}':>14}"
        )

    n = len(per)
    lo, hi = wilson(agg["neutral"], agg["total"])
    lo2, hi2 = wilson(nstr_all, n)
    print(f"\n単純合計: {n} 本、LYS/LYN {agg['total']} 残基")
    print(f"  荷電 {agg['charged']}、中性 {agg['neutral']}、側鎖が不完全 {agg['broken']}")
    print(
        f"  中性率（残基あたり）{agg['neutral'] / agg['total'] * 100:.2f}%  Wilson 95% {lo * 100:.2f}–{hi * 100:.2f}%"
    )
    print(
        f"  中性を 1 本以上含む構造 {nstr_all}/{n} = {nstr_all / n * 100:.1f}%  "
        f"Wilson 95% {lo2 * 100:.1f}–{hi2 * 100:.1f}%"
    )
    print(f"pool_size で重み付けした母集団推定: 残基あたり {wr_res * 100:.2f}%、構造あたり {wr_str * 100:.1f}%")
    print(
        f"\n中性リジンのうち LYS の名前のまま: {agg['neutral_named_LYS']}/{agg['neutral']}、"
        f"LYN と命名された残基: {agg['named_LYN']}"
    )

    print("\n=== 中性または不完全を含む構造 ===")
    for p in sorted(per):
        d = per[p]
        if d["neutral"] or d["broken"]:
            print(
                f"  {p} ({man[p]['stratum']:<10}) LYS 計 {d['total']:>3}: "
                f"荷電 {d['charged']:>3}、中性 {d['neutral']:>2}"
                f"（LYS の名前のまま {d['neutral_named_LYS']}）、不完全 {d['broken']:>2}"
            )


if __name__ == "__main__":
    main()
