#!/usr/bin/env python3
"""滴定状態が動いた残基が、金属や配位子からどれだけ近いかを測る。

  python3 measure_proximity.py <diff.csv> <structure_dir> [--col ligand_removal_changed]

`arms_diff.csv`（配位子の除去順で状態が動いた残基）や `tool_diff.csv`
（ツール間で形式電荷が食い違う残基）を受け取り、各残基について

  - 最も近い金属までの距離
  - 最も近い水以外の HETATM までの距離

を出す。距離は残基の滴定可能原子すべてから測る。

なぜ書いたか
  この測定は以前その場限りのコードで行われ、スクリプトが残っていなかった。
  記録にある「48% が金属から 2.6 Å 以内」「64 残基すべてが配位子から 9.4 Å 以内」
  は再生成できない状態だった。

なぜ測り直すか
  check_receptor.py の W1 が、滴定可能原子の最初の 1 つしか見ずに打ち切って
  いた。HIS は ND1 だけで NE2 を見ず、ARG は NE だけだった。配位圏の判定が
  変わるので、この表も同じ定義で測り直す必要がある。
  --legacy を付けると以前の「最初の 1 原子だけ」の挙動も併記する。

注意
  diff の CSV に鎖の列が無い。同じ残基番号が複数の鎖にある場合は全部の鎖で
  測って最小値を採り、何件が曖昧だったかを別に報告する。
"""

import argparse
import collections
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_receptor import TITR_ATOM, dist, read_pdb


def min_dist(e, rn, targets, first_only=False):
    """残基の滴定可能原子から targets までの最短距離と、その原子名。

    first_only=True は修正前の挙動（存在する最初の 1 原子で打ち切る）。"""
    best, who = None, None
    for a in TITR_ATOM.get(rn, ()):
        if a not in e["atoms"]:
            continue
        d = min(dist(e["atoms"][a], t) for t in targets)
        if best is None or d < best:
            best, who = d, a
        if first_only:
            break
    return best, who


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("diff_csv")
    ap.add_argument("structure_dir")
    ap.add_argument("--col", default=None, help="この列が '1' の行だけを使う（既定は全行）")
    ap.add_argument("--legacy", action="store_true", help="修正前の『最初の1原子だけ』の値も併記する")
    ap.add_argument("--out", help="行ごとの結果を書き出す CSV")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.diff_csv)))
    if a.col:
        rows = [r for r in rows if r.get(a.col) == "1"]
    print(f"{a.diff_csv}: 対象 {len(rows)} 行{f'（{a.col}=1 のみ）' if a.col else ''}")

    cache, out, ambiguous, missing = {}, [], 0, []
    for r in rows:
        pid, num, base = r["pdb_id"], r["resnum"], r["base"]
        if pid not in cache:
            p = os.path.join(a.structure_dir, f"{pid}.pdb")
            cache[pid] = read_pdb(p) if os.path.isfile(p) else None
        got = cache[pid]
        if got is None:
            missing.append(pid)
            continue
        res, het, metals, _nH, _g, _inc = got
        # 鎖が CSV に無いので、同じ残基番号・同じ残基名を全鎖から拾う
        hits = [e for (ch, n, rn), e in res.items() if n == str(num) and rn == base]
        if not hits:
            missing.append(f"{pid}:{base}{num}")
            continue
        if len(hits) > 1:
            ambiguous += 1
        mx = [xyz for _, _, xyz in metals]
        hx = [xyz for _, _, xyz in het]
        rec = dict(
            pdb_id=pid, resnum=num, base=base, n_chains=len(hits), has_metal=int(bool(mx)), has_het=int(bool(hx))
        )
        for tag, tg in (("metal", mx), ("het", hx)):
            if not tg:
                rec[f"d_{tag}"] = ""
                rec[f"atom_{tag}"] = ""
                if a.legacy:
                    rec[f"d_{tag}_legacy"] = ""
                continue
            cands = [min_dist(e, base, tg) for e in hits]
            cands = [c for c in cands if c[0] is not None]
            if not cands:
                rec[f"d_{tag}"] = ""
                rec[f"atom_{tag}"] = ""
                if a.legacy:
                    rec[f"d_{tag}_legacy"] = ""
                continue
            d, who = min(cands, key=lambda c: c[0])
            rec[f"d_{tag}"] = round(d, 2)
            rec[f"atom_{tag}"] = who
            if a.legacy:
                lc = [min_dist(e, base, tg, first_only=True) for e in hits]
                lc = [c for c in lc if c[0] is not None]
                rec[f"d_{tag}_legacy"] = round(min(c[0] for c in lc), 2) if lc else ""
        out.append(rec)

    if missing:
        print(f"  !! 見つからなかった {len(missing)} 件: {missing[:6]}")
    if ambiguous:
        print(f"  !! 残基番号が複数の鎖にある行 {ambiguous} 件（全鎖の最小値を採った）")

    def summarize(key, label, cut):
        vals = [r[key] for r in out if r[key] != ""]
        if not vals:
            print(f"\n{label}: 該当なし")
            return
        vals = sorted(vals)
        n = len(vals)
        within = sum(1 for v in vals if v < cut)
        print(f"\n{label}  対象 {n} 残基")
        print(f"  中央値 {vals[n // 2]:.2f} Å   最小 {vals[0]:.2f}   最大 {vals[-1]:.2f}")
        print(f"  {cut} Å 未満 {within} 件 ({100 * within / n:.0f}%)")
        for t in (4.0, 6.0, 10.0):
            print(
                f"  {t:>4.0f} Å 未満 {sum(1 for v in vals if v < t):>3} 件"
                f" ({100 * sum(1 for v in vals if v < t) / n:>3.0f}%)"
            )

    summarize("d_metal", "金属までの距離", 2.6)
    summarize("d_het", "水以外の HETATM までの距離", 2.6)

    if a.legacy:
        for tag, cut in (("metal", 2.6),):
            new = [r for r in out if r[f"d_{tag}"] != "" and r[f"d_{tag}"] < cut]
            old = [r for r in out if r.get(f"d_{tag}_legacy") != "" and r[f"d_{tag}_legacy"] < cut]
            print(f"\n{tag} {cut} Å 未満の件数")
            print(f"  修正前（最初の1原子だけ）{len(old)}")
            print(f"  修正後（全原子）        {len(new)}")
            gained = [r for r in new if r not in old]
            if gained:
                c = collections.Counter((r["base"], r[f"atom_{tag}"]) for r in gained)
                print("  新たに圏内になった残基: " + ", ".join(f"{b} の {at} {n} 件" for (b, at), n in c.most_common()))

    print("\n原子の内訳（金属に最も近い滴定可能原子）:")
    c = collections.Counter(r["atom_metal"] for r in out if r["atom_metal"])
    for (at), n in c.most_common():
        print(f"  {at:<5} {n}")

    if a.out and out:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0]))
            w.writeheader()
            w.writerows(out)
        print(f"\n-> {a.out}")


if __name__ == "__main__":
    main()
