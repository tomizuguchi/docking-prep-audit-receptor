#!/usr/bin/env python3
"""PDB2PQR を標本に通し、中性リジンが LYS の名前で出る頻度を数える。

  python3 run_pqr.py <pdbdir> <outdir>

プロトコルはサーベイ本体と同じ:
  pdb2pqr --ff=AMBER --with-ph=7.4 --titration-state-method=propka

中性リジンの判定は 2 本立てにする。
  (1) NZ に付く水素の数。荷電 LYS は HZ1/HZ2/HZ3 の 3 本、中性 LYN は 2 本。
  (2) 残基の電荷の和。荷電 LYS は +1、中性 LYN は 0。
AMBER の N/CA/C の電荷で見分ける手もあるが、末端残基で値が変わるので使わない。
2 つの判定が食い違った残基は別に数えて報告する（判定自体の検証になる）。
"""

import collections
import csv
import glob
import os
import subprocess
import sys
import time

CMD = ["pdb2pqr", "--ff=AMBER", "--with-ph=7.4", "--titration-state-method=propka"]


def lys_states(pqr):
    """.pqr の LYS/LYN 残基ごとに (名前, HZ 本数, 電荷の和) を返す。"""
    res = collections.OrderedDict()
    for l in open(pqr, errors="replace"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        f = l.split()
        # pqr は列固定でなく空白区切りだが、座標が繋がると要素数が落ちる。
        # 例: "20.691  54.380-100.128" は 3 要素でなく 2 要素になる。
        # 実際に 6ZHF では ATOM 行 13614 本のうち 9137 本が 9 要素だった。
        # よって要素数で判定してはいけない。前から名前と残基、後ろから
        # 電荷と半径を取り、残基番号は f[4] 以降で最初に整数になる要素とする
        # （鎖 ID 列はあるときとないときがあり、英字なので整数化で弾ける）。
        if len(f) < 8:
            continue
        try:
            q = float(f[-2])
            float(f[-1])
        except ValueError:
            continue
        name, rn = f[2], f[3]
        if rn not in ("LYS", "LYN"):
            continue
        resnum = ""
        for t in f[4:7]:
            t2 = t.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")  # 挿入コード
            if t2.lstrip("-").isdigit():
                resnum = t
                break
        if not resnum:
            continue
        e = res.setdefault((resnum, rn), {"name": rn, "hz": 0, "q": 0.0})
        e["q"] += q
        if name.startswith("HZ"):
            e["hz"] += 1
    return list(res.values())


def main():
    pdbdir, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    rows, fails = [], []
    files = sorted(glob.glob(os.path.join(pdbdir, "*.pdb")))
    for i, f in enumerate(files, 1):
        pid = os.path.basename(f)[:4].upper()
        pqr = os.path.join(outdir, f"{pid}.pqr")
        log = os.path.join(outdir, f"{pid}.pqr.log")
        t0 = time.time()
        if not os.path.exists(pqr):
            p = subprocess.run(CMD + [f, pqr], capture_output=True, text=True, timeout=1800)
            open(log, "w").write(p.stdout + "\n--- stderr ---\n" + p.stderr)
            if p.returncode != 0 or not os.path.exists(pqr):
                tail = [x for x in (p.stderr or p.stdout).strip().split("\n") if x.strip()]
                fails.append(
                    {"pdb_id": pid, "returncode": p.returncode, "error": tail[-1][:300] if tail else "(出力なし)"}
                )
                print(f"[{i:2}/{len(files)}] {pid} 失敗 rc={p.returncode}", flush=True)
                continue
        st = lys_states(pqr)
        # 判定 1: HZ が 2 本、 判定 2: 電荷の和が 0 近傍
        hz_neutral = [s for s in st if s["hz"] == 2]
        q_neutral = [s for s in st if abs(s["q"]) < 0.5]
        named_lyn = [s for s in st if s["name"] == "LYN"]
        disagree = [s for s in st if (s["hz"] == 2) != (abs(s["q"]) < 0.5)]
        rows.append(
            {
                "pdb_id": pid,
                "n_lys_total": len(st),
                "n_neutral_hz": len(hz_neutral),
                "n_neutral_q": len(q_neutral),
                "n_named_LYN": len(named_lyn),
                "n_neutral_but_named_LYS": len([s for s in hz_neutral if s["name"] == "LYS"]),
                "n_disagree": len(disagree),
                "secs": round(time.time() - t0, 1),
            }
        )
        print(
            f"[{i:2}/{len(files)}] {pid} LYS {len(st):3} 中性 {len(hz_neutral):2} "
            f"うち LYS の名前 {len([s for s in hz_neutral if s['name'] == 'LYS']):2} "
            f"({time.time() - t0:.0f}s)",
            flush=True,
        )
    if rows:
        with open(os.path.join(outdir, "lys.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    if fails:
        with open(os.path.join(outdir, "pqr_failed.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(fails[0]))
            w.writeheader()
            w.writerows(fails)
    print(f"\n成功 {len(rows)} / 失敗 {len(fails)} → {outdir}")


if __name__ == "__main__":
    main()
