#!/usr/bin/env python3
"""
PDB構造のプロファイラ。前処理パイプラインのコードパスを洗い出す。

  python3 profile_pdb.py prof  <pdbdir> [out.csv]   各構造のフラグを表にする
  python3 profile_pdb.py pick  <out.csv> [N]        全フラグを覆う最小集合を選ぶ

残基番号は (resSeq, iCode) の組で扱う。int() で潰さない。
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
    # 原子名ごとに altLoc の集合を持つ。res の {atomName: altLoc} は 1 個しか
    # 保持しないので「A も空白も無く B のみ」を表現できない（4MLK MET59 の CA）。
    altsets = collections.defaultdict(set)
    ssbond = 0
    nH = 0
    nmodel = 0
    for l in open(path, errors="replace"):
        rec = l[:6]
        if rec == "EXPDTA":
            method = l[10:].strip()[:40]
        elif l.startswith("REMARK   2 RESOLUTION"):
            # "REMARK   2 RESOLUTION.    1.30 ANGSTROMS." の 1.30 を取る。
            # 行全体を split() すると REMARK 番号の 2 を先に拾うので、
            # "RESOLUTION." より後ろだけを見る。NOT APPLICABLE の行は値なし。
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
                # altLoc の集合には重原子だけ入れる。Meeko に渡す受容体ファイルは
                # 水素を落としてあるので、寄託水素の altLoc を数えると空振りする
                # （1EA7 は altLoc 付き水素 2036 個を持つが Meeko は A/B 両方通る）。
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
    # 修飾残基は ATOM 行で来ることも HETATM 行で来ることもある。MSE などは
    # 普通 HETATM なので、ATOM 行だけ見ていると取り落とす（パイロットで 4/40 見落とし）。
    # 主鎖原子の欠損。PDB2PQR が落ちる原因は側鎖ではなく主鎖の欠損である
    # （パイロット: 落ちた 2 本はどちらも C/CA/O を同時に欠く残基を持ち、
    #  1 原子だけの欠損なら repair_heavy が再構築して通った）。
    BB = {"N", "CA", "C", "O"}
    bbmiss = [(k, v[0], BB - set(v[1])) for k, v in sub.items() if BB - set(v[1])]
    r["bb_missing_res"] = len(bbmiss)
    r["bb_max_missing"] = max((len(m) for _, _, m in bbmiss), default=0)
    # C と O を同時に欠く残基。PDB2PQR の『Too few atoms present to
    # reconstruct or cap residue』はこの条件で出る（100 本で 3/3、空振り 0）。
    r["bb_cap"] = sum(1 for _, _, m in bbmiss if {"C", "O"} <= m)
    r["bb_missing_detail"] = ";".join(f"{rn}{k[1]}:{'+'.join(sorted(m))}" for k, rn, m in bbmiss[:4])

    # 原子単位の altLoc: ある原子が A でも空白でもない altLoc にしか存在しない場合。
    # 残基単位の F_altloc_onlyB では捉えられない（4MLK MET59 は N/C/O が空白で
    # CA だけ altLoc B）。Meeko の --default_altloc A はこの原子を見つけられず落ちる。
    # A も空白も存在せず、別の altLoc にしか無い原子。--default_altloc A では見つからない。
    # 既定 altLoc を X にしたとき、X も空白も持たない原子。Meeko は
    # --default_altloc X でこの原子を見つけられず落ちる（A で 4MLK、B で 3I6H/7USK）。
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


# fetch_strata.py の検索条件で母集団から除かれているため、
# 標本をいくら増やしてもゼロのままになるフラグ。検出漏れと区別する。
FRAME_EXCLUDED = {
    "F_multichain": "deposited_polymer_entity_instance_count == 1",
    "F_nucleic": 'selected_polymer_entity_types == "Protein (only)"',
    "F_no_resol": "resolution_combined に範囲指定",
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
        print(f"{len(rows)} 構造 → {out}\n")
        print("フラグごとの該当構造数")
        for name in FLAGS:
            n = sum(r.get(name) == 1 for r in rows)
            if n:
                mark = ""
            elif name in FRAME_EXCLUDED:
                mark = f"  ← 該当なし。検索条件（{FRAME_EXCLUDED[name]}）で出得ない"
            else:
                mark = "  ← 該当なし。この条件の構造を追加で落とすこと"
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
            print(f"  {best['id']}  新たに覆うフラグ: {' '.join(sorted(gain))}")
        print(f"\nパイロット {len(chosen)} 本: " + " ".join(r["id"] for r in chosen))
        miss = need - covered
        print(
            f"覆えたフラグ {len(covered)}/{len(need)}"
            + (f"、未達: {' '.join(sorted(miss))}" if miss else "（全フラグ達成）")
        )
    else:
        print(__doc__)


if __name__ == "__main__":
    main()
