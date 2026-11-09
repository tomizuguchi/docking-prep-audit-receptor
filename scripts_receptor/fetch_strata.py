#!/usr/bin/env python3
"""
分解能で層化した PDB 構造の取得。乱数シードとクエリを manifest に残す。

  python3 fetch_strata.py list                        各層の件数だけ確認（落とさない）
  python3 fetch_strata.py list  --group               配列90%で重複排除した版の件数
  python3 fetch_strata.py fetch <N> <outdir> [seed] [--group]

同じ seed なら N を増やしたときに前の標本が先頭にそのまま残る
（パイロット 8 本/層 ⊂ 本番 20 本/層）。manifest の draw_rank が抽選順。

既定は return_type=entry、重複排除なし（＝エントリ母集団の標本）。
--group を付けると配列同一性90%の代表のみ（＝タンパク母集団に近い標本）。
どちらを使ったかは manifest.csv と query.json に記録される。
"""

import csv
import json
import os
import random
import sys
import time
import urllib.error
import urllib.request

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
FILES = "https://files.rcsb.org/download/{}.pdb"
STRATA = [
    ("A_le1.5", 0.0, 1.5),
    ("B_1.5-2.0", 1.5, 2.0),
    ("C_2.0-2.5", 2.0, 2.5),
    ("D_2.5-3.0", 2.5, 3.0),
    ("E_ge3.0", 3.0, 9.9),
]


def T(attr, op, val):
    return {"type": "terminal", "service": "text", "parameters": {"attribute": attr, "operator": op, "value": val}}


def query(lo, hi, group=False):
    nodes = [
        T("exptl.method", "exact_match", "X-RAY DIFFRACTION"),
        T("rcsb_entry_info.selected_polymer_entity_types", "exact_match", "Protein (only)"),
        T("rcsb_entry_info.deposited_polymer_entity_instance_count", "equals", 1),
        T(
            "rcsb_entry_info.resolution_combined",
            "range",
            {"from": lo, "include_lower": True, "to": hi, "include_upper": False},
        ),
    ]
    q = {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "polymer_entity" if group else "entry",
        "request_options": {"paginate": {"start": 0, "rows": 10000}, "results_content_type": ["experimental"]},
    }
    if group:
        q["request_options"]["group_by"] = {"aggregation_method": "sequence_identity", "similarity_cutoff": 90}
        q["request_options"]["group_by_return_type"] = "representatives"
    return q


def post(q, label=""):
    req = urllib.request.Request(SEARCH, data=json.dumps(q).encode(), headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=180) as r:
            return json.loads(r.read())
    except urllib.error.HTTPError as e:
        print(f"  !! {label} HTTP {e.code}")
        print("  --- サーバの返答 ---")
        print("  " + e.read().decode()[:1500].replace("\n", "\n  "))
        print("  --- 送ったクエリ ---")
        print("  " + json.dumps(q)[:900])
        return None
    except Exception as e:
        print(f"  !! {label} {type(e).__name__}: {e}")
        return None


def ids_for(lo, hi, group, verbose=True):
    """paginate の上限(10000)を超える層は start を進めて全件集める。"""
    PAGE, got, total, start = 10000, set(), None, 0
    while True:
        q = query(lo, hi, group)
        q["request_options"]["paginate"] = {"start": start, "rows": PAGE}
        d = post(q, f"{lo}-{hi} start={start}")
        if d is None:
            return None, 0
        if total is None:
            total = d.get("total_count", 0)
        rs = d.get("result_set", [])
        if not rs:
            break
        got |= {r["identifier"].split("_")[0] for r in rs}
        start += PAGE
        if verbose and total > PAGE:
            print(f"      ... {len(got)}/{total}")
        if start >= total:
            break
        time.sleep(0.5)
    return sorted(got), total


def pools_for(group):
    """全層の ID を集め、2 層以上に現れるエントリを全層から除いて返す。

    resolution_combined は配列属性で、range 検索は「配列のどれかが範囲内」で
    当たる。分解能値を複数持つエントリは複数の層に同時に入り、層化が崩れる。
    どちらの値が「その構造の分解能」かをここでは決められないので、
    またいだエントリは母集団から外す。外した ID は excluded.csv に残す。
    """
    pools, totals = {}, {}
    for name, lo, hi in STRATA:
        ids, total = ids_for(lo, hi, group)
        if ids is None:
            return None, None, None
        pools[name], totals[name] = set(ids), total
        time.sleep(1)
    names = [n for n, _, _ in STRATA]
    seen, multi, where = set(), set(), {}
    for n in names:
        multi |= seen & pools[n]
        seen |= pools[n]
    for n in names:
        for pid in pools[n] & multi:
            where.setdefault(pid, []).append(n)
    if multi:
        for a in range(len(names)):
            for b in range(a + 1, len(names)):
                ov = pools[names[a]] & pools[names[b]]
                if ov:
                    print(f"  層の重複 {names[a]} ∩ {names[b]}: {len(ov)} 件")
        print(
            f"  分解能値を複数持ち層をまたぐエントリ {len(multi)} 件を全層から除外"
            f"（母集団 {len(seen)} 本の {len(multi) / len(seen) * 100:.2f}%）"
        )
        for n in names:
            pools[n] -= multi
    else:
        print("  層間の重複なし")
    print(f"  有効母集団 {sum(len(p) for p in pools.values())} 本")
    return pools, totals, where


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    group = "--group" in sys.argv
    if not args:
        print(__doc__)
        return
    mode = args[0]
    if mode == "list":
        print(f"重複排除: {'配列同一性90%の代表のみ' if group else 'なし（エントリ母集団）'}")
        pools, totals, where = pools_for(group)
        if pools is None:
            print("\n上のサーバ返答を貼っていただければ直します。")
            return
        for name, lo, hi in STRATA:
            print(f"  {name:<10} {lo}–{hi} Å : 検索 {totals[name]}、重複除外後 {len(pools[name])}")
        print("\nfetch に進めます。")
        return
    if mode != "fetch":
        print(__doc__)
        return
    N, outdir = int(args[1]), args[2]
    seed = int(args[3]) if len(args) > 3 else 20261108
    os.makedirs(outdir, exist_ok=True)
    pools, totals, where = pools_for(group)
    if pools is None:
        print("  取得失敗。中止")
        return
    with open(os.path.join(outdir, "excluded.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["pdb_id", "reason", "strata"])
        for pid, ns in sorted(where.items()):
            w.writerow([pid, "resolution_combined が複数層にまたがる", "|".join(ns)])
    man, failed = [], []
    for name, lo, hi in STRATA:
        ids = sorted(pools[name])
        # 全順列を作って先頭から取る。同じ seed なら N を増やしても
        # 前の標本がそのまま先頭に残る（パイロット ⊂ 本番）。
        order = random.Random(f"{seed}:{name}").sample(ids, len(ids))
        pick, rank, i, bad = [], {}, 0, 0
        while len(pick) < min(N, len(ids)) and i < len(order):
            pid = order[i]
            i += 1
            dst = os.path.join(outdir, f"{pid}.pdb")
            if not os.path.exists(dst):
                try:
                    urllib.request.urlretrieve(FILES.format(pid), dst)
                    time.sleep(0.2)
                except Exception as e:
                    print(f"    !! {pid}: {type(e).__name__} {e}")
                    bad += 1
                    failed.append({"stratum": name, "pdb_id": pid, "draw_rank": i, "reason": f"{type(e).__name__} {e}"})
                    continue
            # 404 や mmCIF のみのエントリで HTML/空ファイルが残ることがある
            try:
                head = open(dst, errors="replace").read(200)
                nl = sum(1 for _ in open(dst, errors="replace"))
            except Exception:
                head, nl = "", 0
            if not head.startswith("HEADER") or nl < 100:
                print(f"    !! {pid}: PDB 形式でない（{nl} 行）。除外")
                os.remove(dst)
                bad += 1
                failed.append({"stratum": name, "pdb_id": pid, "draw_rank": i, "reason": f"PDB 形式でない（{nl} 行）"})
                continue
            rank[pid] = i
            pick.append(pid)
        print(
            f"  {name:<10} 母集団 {len(ids)} から {len(pick)} 本" + (f"（{bad} 本は取得できず飛ばした）" if bad else "")
        )
        for pid in pick:
            man.append(
                {
                    "stratum": name,
                    "res_lo": lo,
                    "res_hi": hi,
                    "pdb_id": pid,
                    "draw_rank": rank[pid],
                    "pool_size": len(ids),
                    "seed": seed,
                    "dedup90": int(group),
                }
            )
    if not man:
        print("  1 本も取得できませんでした")
        return
    if failed:
        with open(os.path.join(outdir, "failed.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(failed[0]))
            w.writeheader()
            w.writerows(failed)
        print(f"  取得できなかった {len(failed)} 本は failed.csv に記録（次順位で補充済み）")
    with open(os.path.join(outdir, "manifest.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(man[0]))
        w.writeheader()
        w.writerows(man)
    with open(os.path.join(outdir, "query.json"), "w") as fh:
        json.dump(
            {
                "seed": seed,
                "dedup90": group,
                "endpoint": SEARCH,
                "exclusion_rule": (
                    "resolution_combined は配列属性のため複数層に当たる"
                    "エントリがある。2 層以上に現れた ID は全層から除外した"
                    "（excluded.csv）。層ごとの母集団は除外後の数で、"
                    "manifest.csv の pool_size がそれ。"
                ),
                "n_excluded": len(where),
                "strata": [
                    {
                        "name": n,
                        "lo": lo,
                        "hi": hi,
                        "query": query(lo, hi, group),
                        "search_total": totals[n],
                        "pool_after_exclusion": len(pools[n]),
                    }
                    for n, lo, hi in STRATA
                ],
            },
            fh,
            indent=1,
        )
    print(f"\n{len(man)} 本 → {outdir}（manifest.csv と query.json 同梱）")


if __name__ == "__main__":
    main()
