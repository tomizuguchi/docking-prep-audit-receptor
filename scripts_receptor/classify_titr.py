#!/usr/bin/env python3
"""pqr の滴定可能残基の状態を、残基名ではなく水素の有無から決める。

  python3 classify_titr.py <pqrdir> [out.csv]

残基名を信用しない理由: PDB2PQR は中性リジンに LYN のパラメータを当てても
名前を LYS のままにする（100 構造で 18/18、HSA 3 構造で 6/6）。
名前で読むと中性リジンは 1 本も見えない。だから水素で読む。

判定（AMBER の水素命名）
  ASP  HD2 あり        -> 中性(ASH)        なし -> 荷電(-1)
  GLU  HE2 あり        -> 中性(GLH)        なし -> 荷電(-1)
  HIS  HD1 と HE2 両方 -> HIP(+1)
       HD1 のみ        -> HID(中性)        HE2 のみ -> HIE(中性)
  CYS  HG あり         -> 中性チオール      なし -> SS 結合 か チオラート(-1)
       SG-SG が 2.5 A 未満の相手を持てば SS 結合、なければチオラート
  TYR  HH あり         -> 中性             なし -> チロシナート(-1)
  LYS  HZ 3 本         -> 荷電(+1)         2 本 -> 中性(LYN)
  ARG  NE/NH1/NH2 の H 5 本 -> 荷電(+1)    4 本 -> 中性(AR0)
分類できなかった残基は捨てずに unknown として数え、必ず報告する。

パーサの注意: pqr は空白区切りだが座標が繋がると要素数が落ちる
（6ZHF では ATOM 行の 67% が 9 要素）。要素数で行を弾いてはいけない。
座標は列固定(30-54)で読む。滴定可能残基の ATOM 行を 1 本も落とさないことを
assert で保証する。
"""

import collections
import csv
import glob
import os
import sys

import numpy as np

TITR = {
    "ASP",
    "ASH",
    "GLU",
    "GLH",
    "HIS",
    "HID",
    "HIE",
    "HIP",
    "CYS",
    "CYM",
    "CYX",
    "LYS",
    "LYN",
    "TYR",
    "TYM",
    "ARG",
    "AR0",
}
# 名前から素の残基種へ
BASE = {
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


def parse(pqr):
    """滴定可能残基ごとの情報と、取りこぼし検算用の行数を返す。"""
    res, seen, placed = collections.OrderedDict(), 0, 0
    for l in open(pqr, errors="replace"):
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
        # 鎖 ID は残基名と残基番号の間に入る。PDB2PQR は既定で落とすが
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
        name = f[2]
        e = res.setdefault(
            (chain, num, rn), {"base": BASE.get(rn, rn), "name": rn, "at": set(), "q": 0.0, "sg": None, "dup": 0}
        )
        # 同じ原子名が二度来たら、別の鎖の同番号残基が併合されたか altLoc である。
        # set に入れて黙って消すと、片方の HD1 と他方の HE2 が 1 残基に見え、
        # 架空の HIP(+1) が生まれる。数えて呼び出し側に報告する。
        if name in e["at"]:
            e["dup"] += 1
        e["at"].add(name)
        e["q"] += q
        if name == "SG":
            try:
                e["sg"] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            except ValueError:
                pass
    return res, seen, placed


def classify(res):
    """残基ごとに状態を決める。CYS は SG-SG 距離が要るので全体を見る。"""
    sgs = [(k, v["sg"]) for k, v in res.items() if v["base"] == "CYS" and v["sg"] is not None]
    ss = set()
    for i in range(len(sgs)):
        for j in range(i + 1, len(sgs)):
            if np.linalg.norm(sgs[i][1] - sgs[j][1]) < 2.5:
                ss.add(sgs[i][0])
                ss.add(sgs[j][0])
    out = {}
    for k, v in res.items():
        a, b = v["at"], v["base"]
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
            s = "unknown"
        out[k] = (b, v["name"], s)
    return out


def main():
    pqrdir = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else None
    rows = []
    for f in sorted(glob.glob(os.path.join(pqrdir, "*.pqr"))):
        pid = os.path.basename(f)[:4].upper()
        res, seen, placed = parse(f)
        assert seen == placed, f"{pid}: 滴定可能残基の行 {seen} 本中 {placed} 本しか割り当てていない"
        for (num, rn), (base, name, state) in classify(res).items():
            rows.append({"pdb_id": pid, "resnum": num, "base": base, "name_in_pqr": name, "state": state})
    if out:
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    print(f"{len({r['pdb_id'] for r in rows})} 本、滴定可能残基 {len(rows)} 個 （取りこぼしゼロを検算済み）")
    unk = [r for r in rows if r["state"].startswith("unknown")]
    print(
        f"分類できなかった残基 {len(unk)} 個"
        + (f": {collections.Counter((r['base'], r['state']) for r in unk).most_common(6)}" if unk else "")
    )
    print("\n=== PDB2PQR が名前に状態を反映するか ===")
    print(f"{'素の残基':<9}{'状態':<14}{'個数':>7}  pqr 中の名前")
    tab = collections.defaultdict(collections.Counter)
    for r in rows:
        tab[(r["base"], r["state"])][r["name_in_pqr"]] += 1
    for (base, state), names in sorted(tab.items()):
        shown = ", ".join(f"{n} {c}" for n, c in names.most_common())
        hidden = (
            "  <- 名前で区別できない" if len(names) == 1 and state not in ("charged+", "charged-", "disulfide") else ""
        )
        print(f"{base:<9}{state:<14}{sum(names.values()):>7}  {shown}{hidden}")
    if out:
        print(f"\n-> {out}")


if __name__ == "__main__":
    main()
