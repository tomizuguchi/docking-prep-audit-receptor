#!/usr/bin/env python3
"""ドッキング受容体の前処理の前に入力を検査し、後で出力を検算する。

  python3 check_receptor.py <input.pdb> [--altloc A] [--allow-bad-res]
                              [--prepared <out.pqr|out.pdbqt>] [--json]

  --allow-bad-res は Meeko に -a を付けるつもりのときに指定する。
  不完全側鎖の指摘が ERROR から WARN に下がる（落ちはしないが残基が消える）。

検査する規則はすべて 100 構造の層化標本での測定に裏付けがある。
依存はなし（標準ライブラリのみ）。

前処理が落ちる規則（測定では誤検出ゼロ・見逃しゼロ）
  E1  C と O を同時に欠く残基がある -> PDB2PQR が cap できず落ちる（3/3）
  E2  --default_altloc X に対し、X も空白も持たない重原子がある
      -> Meeko が落ちる（X=A で 1/1、X=B で 2/2）

前処理が黙って受容体を変える規則
  W1  金属がある。HETATM を落としてから滴定させると配位残基の電荷が反転する
      （金属あり 41 本中 23 本で Σq が動き最大 4 単位、金属なしは最大 1 単位）
      滴定可能原子は全部見る。最初の 1 つで打ち切っていたときは、金属を含む
      44 構造で配位圏内の残基を 73 個しか拾えていなかった（正しくは 124 個。
      取りこぼしの内訳は HIS の NE2 が 23、ASP の OD2 が 18、GLU の OE2 が 9、
      ARG の NH2 が 1。HIS は ND1 しか見ていなかった）
  W2  水以外の HETATM がある。落とす順で滴定状態が変わる（64 残基 / 34 本）
  W3  寄託水素がある。PDB2PQR が鎖切断と誤報告して落ちることがある（1EA7）
  W4  鎖の中の修飾残基が HETATM 行で来ている。落とすと主鎖に切断ができる（6/100）
      隣接残基との C–N 距離で判定する。残基名の一覧で判定していたときは
      MSE しか拾えず OSE と CSS を取りこぼし、逆に主鎖の原子名を持つだけの
      SAM や遊離 PRO を拾う危険があった
  W5  側鎖が不完全な残基がある。-a なしでは Meeko が落ちる（37/37、空振り 0）。
      -a を付けると通るがこの残基が削除され、荷電残基に 2.9 倍偏る
      （306 残基 / 38 本、形式電荷が最大 9 単位動く）。
      -a なしの失敗 39 本のうち 37 本がこれ、1 本は altLoc（E2）、1 本（4LYZ）は未解明

出力の検算
  C1  pqr の残基名が状態を表していない（--ffout を付け忘れた）
  C2  pdbqt から残基が消えている（Meeko の -a による削除）。消えた分の形式電荷も出す
  C3  2 本の pqr で滴定状態が違う残基（--prepared A.pqr --against B.pqr）
      HETATM を入れた入力と抜いた入力の両方で PDB2PQR を走らせて比べる。
      pqr には出自が記録されない（ヘッダも REMARK も無い）ので、片方だけでは
      どちらの条件で滴定したか判定できない。痕跡は変化した残基の水素 1 個だけで、
      それは正しい状態と見分けがつかない。
      100 構造の標本では 64 残基 / 34 構造で状態が動いた（arms_diff.csv と一致）
"""

import collections
import json
import math
import os
import sys

METALS = set(
    [
        "LI",
        "NA",
        "K",
        "RB",
        "CS",
        "MG",
        "CA",
        "SR",
        "BA",
        "SC",
        "Y",
        "LA",
        "CE",
        "TI",
        "ZR",
        "HF",
        "V",
        "NB",
        "TA",
        "CR",
        "MO",
        "W",
        "MN",
        "TC",
        "RE",
        "FE",
        "RU",
        "OS",
        "CO",
        "RH",
        "IR",
        "NI",
        "PD",
        "PT",
        "CU",
        "AG",
        "AU",
        "ZN",
        "CD",
        "HG",
        "AL",
        "GA",
        "IN",
        "TL",
        "SN",
        "PB",
        "BI",
        "U",
        "TH",
        "EU",
        "GD",
        "TB",
        "DY",
        "HO",
        "ER",
        "TM",
        "YB",
        "LU",
        "SM",
        "ND",
        "PR",
    ]
)
STD = set(
    [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    ]
)
# 鎖の中に入りうる修飾残基（HETATM 行で来ることが多い）
INCHAIN_MOD = {"MSE", "SEP", "TPO", "PTR", "CSO", "CME", "KCX", "MLY", "ALY", "HYP", "PCA", "LLP"}
BB = {"N", "CA", "C", "O"}
# 完全な側鎖の重原子数（主鎖 N CA C O を含む）
_SIDE = {
    "ALA": "CB",
    "ARG": "CB CG CD NE CZ NH1 NH2",
    "ASN": "CB CG OD1 ND2",
    "ASP": "CB CG OD1 OD2",
    "CYS": "CB SG",
    "GLN": "CB CG CD OE1 NE2",
    "GLU": "CB CG CD OE1 OE2",
    "GLY": "",
    "HIS": "CB CG ND1 CD2 CE1 NE2",
    "ILE": "CB CG1 CG2 CD1",
    "LEU": "CB CG CD1 CD2",
    "LYS": "CB CG CD CE NZ",
    "MET": "CB CG SD CE",
    "PHE": "CB CG CD1 CD2 CE1 CE2 CZ",
    "PRO": "CB CG CD",
    "SER": "CB OG",
    "THR": "CB OG1 CG2",
    "TRP": "CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
    "TYR": "CB CG CD1 CD2 CE1 CE2 CZ OH",
    "VAL": "CB CG1 CG2",
}
# 完全な残基の重原子名の集合（主鎖 N CA C O を含む）。
# 個数ではなく集合で比べる。個数だと、NZ を欠く LYS に OXT が付いていると
# 9 原子で「完全」と判定され、W5 が空振りする。
ATOMSET = {rn: BB | set(s.split()) for rn, s in _SIDE.items()}
NHEAVY = {rn: len(s) for rn, s in ATOMSET.items()}
TITR_ATOM = {
    "ASP": ("OD1", "OD2"),
    "GLU": ("OE1", "OE2"),
    "HIS": ("ND1", "NE2"),
    "LYS": ("NZ",),
    "CYS": ("SG",),
    "TYR": ("OH",),
    "ARG": ("NE", "NH1", "NH2"),
}
FORMAL = {"LYS": 1, "ARG": 1, "ASP": -1, "GLU": -1}


def dist(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def guess_element(l):
    """(元素記号, 元素列から読めたか)。

    77-78 列が空の PDB では原子名から推定するほかないが、先頭 1 文字を取ると
    ZN が 'Z' になり金属を取りこぼす。13-14 列に右詰めされた 2 文字が金属名に
    一致する場合はそちらを採る（" CA " は Cα、"CA  " はカルシウム）。
    推定に落ちた件数は呼び出し側で数えて警告する。"""
    el = l[76:78].strip().upper()
    if el:
        return el, True
    name = l[12:16] if len(l) >= 16 else l[12:].ljust(4)
    two = name[0:2].strip().upper()
    if len(two) == 2 and two in METALS:
        return two, False
    raw = name.strip().upper()
    if raw and raw[0].isdigit():
        raw = raw[1:]  # 1HB のような PDB v2 の水素名
    return (raw[0] if raw else ""), False


def find_inchain(hetres, chainres):
    """HETATM 行で来ている残基のうち、ペプチド結合で鎖に繋がっているものを返す。

    固定の残基名一覧では足りない。100 構造の標本では、一覧が MSE しか拾わず
    OSE / 0QF / DSN / CSS を取りこぼしていた。逆に SAM や遊離 PRO は主鎖の
    原子名（N, CA, C）を持つが鎖には繋がっておらず、名前や原子名だけでは
    区別できない。C–N 距離で判定すればどちらも正しく分かれる。"""
    PEP = 1.5  # ペプチド結合 C–N は約 1.33 Å
    out = {}
    for k, e in hetres.items():
        n, c = e["atoms"].get("N"), e["atoms"].get("C")
        linked = False
        for k2, e2 in chainres.items():
            if k2 == k:
                continue
            c2, n2 = e2["atoms"].get("C"), e2["atoms"].get("N")
            if n and c2 and dist(n, c2) < PEP:
                linked = True
                break
            if c and n2 and dist(c, n2) < PEP:
                linked = True
                break
        if linked:
            out[k] = e
    return out


def read_pdb(path):
    """ATOM/HETATM を読む。altLoc は原子名ごとに集合で持つ。"""
    res = collections.OrderedDict()  # (chain, resseq+icode, resname) -> info
    hetres = collections.OrderedDict()
    het, metals, nH, guessed = [], [], 0, 0
    for l in open(path, errors="replace"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        el, from_col = guess_element(l)
        if not from_col:
            guessed += 1
        name = l[12:16].strip()
        alt = l[16] if len(l) > 16 else " "
        rn = l[17:20].strip()
        key = (l[21] if len(l) > 21 else " ", l[22:27].strip(), rn)
        try:
            xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        except ValueError:
            continue
        if el == "H" or el == "D":
            nH += 1
            continue  # 重原子だけを数える（受容体ファイルは H を落とす）
        if l.startswith("HETATM"):
            if rn != "HOH":
                het.append((rn, name, xyz))
                if el in METALS:
                    metals.append((rn, el, xyz))
            # 鎖の中の修飾残基は HETATM 行で来るが、主鎖の一部である。
            # ここでは一旦すべて別に溜め、あとで結合距離で選り分ける。
            if rn != "HOH":
                e = hetres.setdefault(key, {"atoms": {}, "alt": collections.defaultdict(set)})
                e["atoms"][name] = xyz
                e["alt"][name].add(alt)
            continue
        e = res.setdefault(key, {"atoms": {}, "alt": collections.defaultdict(set)})
        e["atoms"][name] = xyz
        e["alt"][name].add(alt)
    # 鎖に繋がっている HETATM 残基を res に合流させる。res に入れないと
    # E1（C と O を同時に欠く）や W5（不完全側鎖）の検査が届かない。
    inchain = find_inchain(hetres, res)
    for k, e in inchain.items():
        res[k] = e
    return res, het, metals, nH, guessed, inchain


def _min_titr_dist(e, rn, targets):
    """残基の滴定可能原子すべてについて、target 群までの最短距離と、その原子名。

    以前は最初に見つかった 1 原子だけを見て break していた。HIS は ND1 しか
    見ず NE2 を見ないので、NE2 で金属に配位するヒスチジンを配位圏外に分類して
    いた。ARG も NE だけで NH1 / NH2 を見ていなかった。"""
    best, who = None, None
    for a in TITR_ATOM.get(rn, ()):
        if a not in e["atoms"]:
            continue
        d = min(dist(e["atoms"][a], xyz) for xyz in targets)
        if best is None or d < best:
            best, who = d, a
    return best, who


def check_input(path, altlocs=("A", "B"), allow_bad_res=False):
    res, het, metals, nH, guessed, inchain = read_pdb(path)
    out = []
    if guessed:
        out.append(
            (
                "W0",
                "WARN",
                f"元素列（77-78 桁）が空の原子が {guessed} 個ある。原子名から推定した",
                ["推定は金属の取りこぼしを起こしうる"],
                "元素列を埋めたファイルを使う。以降の金属判定の信頼度が下がる",
            )
        )

    # --- E1: C と O を同時に欠く残基 ---
    cap = []
    for (ch, num, rn), e in res.items():
        if rn not in STD and rn not in INCHAIN_MOD:
            continue
        miss = BB - set(e["atoms"])
        if {"C", "O"} <= miss:
            cap.append(f"{rn}{num}({ch}) 欠損 {'+'.join(sorted(miss))}")
    if cap:
        out.append(
            (
                "E1",
                "ERROR",
                "PDB2PQR は C と O を同時に欠く残基を cap できず落ちる",
                cap,
                "その残基を削除するか、主鎖を外部ツールで補完する",
            )
        )

    # --- E2: altLoc 変種が無い重原子 ---
    for X in altlocs:
        bad = [
            f"{rn}{num}({ch}):{a}"
            for (ch, num, rn), e in res.items()
            for a, s in e["alt"].items()
            if not (s & {" ", X})
        ]
        if bad:
            out.append(
                (
                    f"E2[{X}]",
                    "ERROR",
                    f"--default_altloc {X} では見つからない原子がある。Meeko が落ちる",
                    bad,
                    "--wanted_altloc で個別に指定するか、別の altLoc を既定にする",
                )
            )
    ok = [X for X in altlocs if not any(not (s & {" ", X}) for e in res.values() for s in e["alt"].values())]
    if len(altlocs) > 1:
        out.append(
            (
                "I0",
                "INFO",
                "安全な --default_altloc",
                [", ".join(ok) if ok else "この構造には安全な既定 altLoc がない"],
                "",
            )
        )

    # --- W1 / W2: HETATM と金属 ---
    if metals:
        kinds = collections.Counter(e for _, e, _ in metals)
        coord, nearby = [], []
        mxyz = [xyz for _, _, xyz in metals]
        for (ch, num, rn), e in res.items():
            d, a = _min_titr_dist(e, rn, mxyz)
            if d is None:
                continue
            if d < 2.6:
                coord.append(f"{rn}{num}({ch}) {a} 金属まで {d:.2f} Å")
            elif d < 6.0:
                nearby.append(f"{rn}{num}({ch}) {a} {d:.1f} Å")
        out.append(
            (
                "W1",
                "WARN",
                f"金属 {dict(kinds)} がある。HETATM を落としてから滴定させると配位残基の電荷が反転する",
                [f"配位圏内（2.6 Å 未満）の滴定可能残基 {len(coord)} 個"]
                + coord
                + [f"2.6–6 Å の滴定可能残基 {len(nearby)} 個"],
                "金属を入力に残したまま PDB2PQR にかける。"
                "落とす操作は要らない。PDB2PQR は水以外の HETATM を出力に書かないので、"
                "入力に残しても pqr には現れない",
            )
        )
    if het:
        ligs = sorted({rn for rn, _, _ in het})
        near = []
        hxyz = [xyz for _, _, xyz in het]
        for (ch, num, rn), e in res.items():
            d, a = _min_titr_dist(e, rn, hxyz)
            if d is not None and d < 6.0:
                near.append(f"{rn}{num}({ch}) {a} {d:.1f} Å")
        out.append(
            (
                "W2",
                "WARN",
                f"水以外の HETATM がある（{', '.join(ligs[:8])}"
                f"{' ほか' if len(ligs) > 8 else ''}）。落とす順で滴定状態が変わる",
                [f"配位子から 6 Å 以内の滴定可能残基 {len(near)} 個"],
                "HETATM を入れた入力と抜いた入力の両方で PDB2PQR を走らせ、"
                "--prepared A.pqr --against B.pqr で差を確認する。"
                "pqr には出自が記録されないので、片方だけでは判定できない",
            )
        )

    # --- W3: 寄託水素 ---
    if nH:
        out.append(
            (
                "W3",
                "WARN",
                f"寄託水素が {nH} 個ある。PDB2PQR が鎖切断と誤報告して落ちることがある",
                [],
                "水素を落としてから渡す。各ツールが自前で付け直す",
            )
        )

    # --- W4: 鎖の中の修飾残基が HETATM 行 ---
    if inchain:
        kinds = sorted({k[2] for k in inchain})
        out.append(
            (
                "W4",
                "WARN",
                f"鎖の中の修飾残基が HETATM 行で来ている（{', '.join(kinds)}）",
                [f"{k[2]}{k[1]}({k[0]})" for k in list(inchain)[:12]]
                + ([f"... 他 {len(inchain) - 12} 個"] if len(inchain) > 12 else []),
                "HETATM を一括削除すると主鎖に切断ができる。この残基は残す",
            )
        )

    # --- W5: 側鎖が不完全な残基 ---
    inc, q_at_risk = [], 0
    for (ch, num, rn), e in res.items():
        want = ATOMSET.get(rn)
        if not want:
            continue
        # 個数ではなく集合で比べる。OXT のような余分な原子が欠損を埋め合わせて
        # 見逃しになるのを防ぐ。
        miss = want - set(e["atoms"])
        if miss:
            inc.append(f"{rn}{num}({ch}) 欠損 {'+'.join(sorted(miss))}")
            q_at_risk += FORMAL.get(rn, 0)
    if inc:
        # -a を付けなければ Meeko は落ちる（100 本で 37/37、空振り 0）。
        # 付ければ通るが、この残基が削除される（38 本で 306 残基、荷電残基に 2.9 倍偏る）。
        lvl = "WARN" if allow_bad_res else "ERROR"
        msg = (
            f"側鎖が不完全な残基が {len(inc)} 個。Meeko の -a はこれを削除する"
            if allow_bad_res
            else f"側鎖が不完全な残基が {len(inc)} 個。-a なしでは Meeko が落ちる"
            f"（excess inter-residue bond として報告される）"
        )
        out.append(
            (
                "W5",
                lvl,
                msg,
                inc[:12]
                + ([f"... 他 {len(inc) - 12} 個"] if len(inc) > 12 else [])
                + [f"-a で削除された場合の形式電荷の変化 {q_at_risk:+d} 単位"],
                "側鎖を補完してから渡す（PDB2PQR を通すと埋まる）。"
                "-a で削除を許すなら、どの残基が消えたかを出力で検算する",
            )
        )
    return out, res


PQR_BASE = {
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


def pqr_states(path):
    """pqr の滴定可能残基ごとの状態。残基名ではなく水素の有無で決める。

    PDB2PQR は中性リジンに LYN のパラメータを当てても名前を LYS のままにする
    （--ffout を付けない場合）。名前で読むと中性リジンが 1 本も見えない。"""
    res = parse_pqr(path)
    sg = {k: v["sg"] for k, v in res.items() if PQR_BASE.get(v["rn"], v["rn"]) == "CYS" and v["sg"]}
    ss = set()
    ks = list(sg)
    for i in range(len(ks)):
        for j in range(i + 1, len(ks)):
            if dist(sg[ks[i]], sg[ks[j]]) < 2.5:
                ss.add(ks[i])
                ss.add(ks[j])
    out = {}
    for k, v in res.items():
        a, b = v["at"], PQR_BASE.get(v["rn"], v["rn"])
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
            continue
        out[k] = (b, s)
    return out


def check_against(a_path, b_path):
    """C3: 2 本の pqr で滴定状態が違う残基。

    pqr には出自が記録されない。配位子を入れて滴定したのか抜いて滴定したのか、
    ファイル 1 本からは判定できない（痕跡は変化した残基の水素 1 個だけで、
    それは「正しい状態」と見分けがつかない）。両方を突き合わせるしかない。"""
    A, B = pqr_states(a_path), pqr_states(b_path)
    common = set(A) & set(B)
    diff = sorted(k for k in common if A[k][1] != B[k][1])
    only = set(A) ^ set(B)
    det = [f"{A[k][0]}{k[1]}{('(' + k[0] + ')') if k[0] else ''} {A[k][1]} -> {B[k][1]}" for k in diff[:20]]
    if len(diff) > 20:
        det.append(f"... 他 {len(diff) - 20} 個")
    if only:
        det.append(f"片方にしかない滴定可能残基 {len(only)} 個")
    if not diff and not only:
        return [("C3", "INFO", "2 本の pqr で滴定状態は完全に一致した", [f"比較した滴定可能残基 {len(common)} 個"], "")]
    return [
        (
            "C3",
            "WARN",
            f"2 本の pqr で滴定状態が違う残基が {len(diff)} 個",
            [f"比較した滴定可能残基 {len(common)} 個"] + det,
            "どちらを使うかで受容体の電荷が変わる。結合部位にかかる残基が含まれていないか確認する",
        )
    ]


def parse_pqr(path):
    """pqr を読む。鎖・残基番号・残基名で鍵を作り、原子名と SG 座標を持つ。"""
    res = collections.OrderedDict()
    for l in open(path, errors="replace"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        f = l.split()
        if len(f) < 8:
            continue
        rn = f[3]
        chain, num = "", ""
        for j in range(4, min(7, len(f))):
            if f[j].rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ").lstrip("-").isdigit():
                num = f[j]
                if j > 4:
                    chain = f[j - 1]  # 残基名と残基番号の間にあるのが鎖 ID
                break
        if not num:
            continue
        # 鎖をキーに含める。含めないと、別の鎖の同じ番号の残基が併合され、
        # 片方の HD1 と他方の HE2 が 1 残基に見えて架空の「+1 の HIS」になる。
        e = res.setdefault((chain, num, rn), {"rn": rn, "at": set(), "sg": None})
        e["at"].add(f[2])
        if f[2] == "SG":
            try:
                e["sg"] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            except ValueError:
                pass
    return res


def check_pqr(path):
    """C1: 残基名が状態を表しているか。状態は水素の数で読む。"""
    res = parse_pqr(path)
    hidden = collections.Counter()
    for (_ch, _num, rn), e in res.items():
        at = e["at"]
        if rn in ("ASP", "ASH") and "HD2" in at and rn == "ASP":
            hidden["中性 ASP が ASP のまま"] += 1
        if rn in ("GLU", "GLH") and "HE2" in at and rn == "GLU":
            hidden["中性 GLU が GLU のまま"] += 1
        if rn == "LYS" and sum(1 for x in at if x.startswith("HZ")) == 2:
            hidden["中性 LYS が LYS のまま"] += 1
        if rn == "HIS" and "HD1" in at and "HE2" in at:
            hidden["+1 の HIS が HIS のまま"] += 1
        if rn == "HIS" and (("HD1" in at) != ("HE2" in at)):
            hidden["互変異性体が HIS のまま"] += 1
        if rn == "CYS" and "HG" not in at:
            hidden["SS 結合かチオラートが CYS のまま"] += 1
    if hidden:
        return [
            (
                "C1",
                "WARN",
                "pqr の残基名が状態を表していない。--ffout=AMBER を付け忘れている",
                [f"{k}: {v} 残基" for k, v in hidden.most_common()],
                "--ffout=AMBER を付け直す。名前から状態を読む下流ツールが取り違える",
            )
        ]
    return []


def check_pdbqt(path, src_res):
    """C2: 入力にあって pdbqt に無い残基と、その形式電荷。"""
    # 残基名をキーに含めない。前処理ツールは HIS を HIE / HID / HIP に、
    # ASP を ASH に改名するので、名前込みで突き合わせると改名された残基を
    # すべて「削除された」と報告してしまう。
    got = {(l[21], l[22:27].strip()) for l in open(path, errors="replace") if l.startswith(("ATOM", "HETATM"))}
    gone = [k for k in src_res if (k[0], k[1]) not in got]
    if not gone:
        return []
    q = sum(FORMAL.get(k[2], 0) for k in gone)
    ch = sum(1 for k in gone if k[2] in FORMAL)
    return [
        (
            "C2",
            "WARN",
            f"入力にあって pdbqt に無い残基が {len(gone)} 個",
            [f"{k[2]}{k[1]}({k[0]})" for k in gone[:12]]
            + ([f"... 他 {len(gone) - 12} 個"] if len(gone) > 12 else [])
            + [f"うち荷電残基 {ch} 個、形式電荷の変化 {q:+d} 単位"],
            "削除された残基が結合部位にかかっていないか確認する",
        )
    ]


def main():
    # 自前の argv 走査をやめて argparse に任せる。以前は
    #   args = [a for a in sys.argv[1:] if not a.startswith("--")]
    # としてフラグだけ除き「その値」を残していたため、
    #   check_receptor.py --prepared out.pdbqt in.pdb
    # と書くと out.pdbqt が検査対象になり、本来 ERROR になる入力が
    # ERROR 0 件・exit 0 で通った。--altloc を末尾に置くと IndexError でも落ちた。
    import argparse

    ap = argparse.ArgumentParser(
        description="ドッキング受容体の前処理の前に入力を検査し、後で出力を検算する",
        epilog="終了コード: 0=阻害要因なし  1=ERROR あり  2=使い方または入出力の誤り",
    )
    ap.add_argument("input", help="受容体の PDB ファイル")
    ap.add_argument("--altloc", metavar="X", help="--default_altloc に渡すつもりの altLoc")
    ap.add_argument(
        "--allow-bad-res",
        action="store_true",
        help="Meeko に -a を付けるつもりのとき。不完全側鎖を ERROR から WARN に下げる",
    )
    ap.add_argument("--prepared", metavar="FILE", help="検算する出力（.pqr または .pdbqt）")
    ap.add_argument(
        "--against",
        metavar="FILE",
        help="--prepared の pqr と突き合わせるもう 1 本の pqr。"
        "HETATM を入れた入力と抜いた入力の両方で PDB2PQR を走らせ、"
        "滴定状態が動いた残基を出す",
    )
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()

    src, as_json = a.input, a.json
    for p in (src, a.prepared, a.against):
        if p and not os.path.isfile(p):
            print(f"ファイルが無い: {p}", file=sys.stderr)
            return 2
    altlocs = (a.altloc,) if a.altloc else ("A", "B")
    try:
        findings, res = check_input(src, altlocs, a.allow_bad_res)
    except Exception as e:
        print(f"{src} の読み取りに失敗: {type(e).__name__}: {e}", file=sys.stderr)
        return 2
    prep = a.prepared
    if prep:
        if prep.endswith(".pqr"):
            findings += check_pqr(prep)
        elif prep.endswith(".pdbqt"):
            findings += check_pdbqt(prep, set(res))
        else:
            print(f"--prepared は .pqr か .pdbqt のみ扱える: {prep}", file=sys.stderr)
            return 2
    if a.against:
        if not prep or not prep.endswith(".pqr"):
            print("--against は --prepared に .pqr を指定したときだけ使える", file=sys.stderr)
            return 2
        if not a.against.endswith(".pqr"):
            print(f"--against は .pqr のみ扱える: {a.against}", file=sys.stderr)
            return 2
        findings += check_against(prep, a.against)
    nerr = sum(1 for f in findings if f[1] == "ERROR")
    nwarn = sum(1 for f in findings if f[1] == "WARN")
    if as_json:
        print(
            json.dumps(
                {
                    "input": src,
                    "prepared": prep,
                    "against": a.against,
                    "n_error": nerr,
                    "n_warn": nwarn,
                    "findings": [
                        {"id": fid, "level": lvl, "message": msg, "detail": det, "action": act}
                        for fid, lvl, msg, det, act in findings
                    ],
                },
                ensure_ascii=False,
                indent=1,
            )
        )
        return 1 if nerr else 0
    print(f"=== {os.path.basename(src)} ===")
    print(f"ERROR {nerr} 件、WARN {nwarn} 件\n")
    for fid, lvl, msg, det, act in findings:
        print(f"[{lvl:<5} {fid}] {msg}")
        for d in det:
            print(f"        {d}")
        if act:
            print(f"        → {act}")
        print()
    if not findings:
        print("検査した規則には該当なし。")
    return 1 if nerr else 0


if __name__ == "__main__":
    sys.exit(main())
