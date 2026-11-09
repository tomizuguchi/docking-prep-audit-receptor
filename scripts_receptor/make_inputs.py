#!/usr/bin/env python3
"""取得したままの PDB から、準備の 3 条件の入力を作る。

  python3 make_inputs.py <srcdir> <outbase>

`<srcdir>` に `*.pdb` を置く（`full100_manifest.csv` の 100 構造を RCSB から
取得したもの）。`<outbase>` の下に 3 つのディレクトリができる。

  in_hetall/   重原子 + HETATM 全部（配位子 + 水）
  in_nowat/    重原子 + 水以外の HETATM
  recept/      重原子のみ（水も配位子も除去）

規則
  ATOM と HETATM の行だけを残す。元素欄（77-78 桁）が H または D のものは落とす。
  TER と END はそのまま残す。それ以外の行（HEADER, REMARK, ANISOU, CONECT など）は
  落とす。altloc 指示子はそのまま残す。原子番号は書き換えない。

  条件ごとの違いは HETATM の扱いだけ。
    in_hetall  HETATM をすべて残す
    in_nowat   残基名が HOH の HETATM を落とす
    recept     HETATM をすべて落とす

元素欄が空の行は、原子名（13-16 桁）の先頭の英字で判定する。PDB の原子名は
元素記号が 2 文字のときだけ 13 桁目から始まるので、13 桁目が空白なら 14 桁目の
英字を元素とみなす。
"""

import argparse
import glob
import os
import sys

KEEP_OTHER = ("TER", "END")
DROP_ELEM = ("H", "D")

ARMS = (
    ("in_hetall", "all"),
    ("in_nowat", "nowat"),
    ("recept", "none"),
)


def element_of(line):
    """元素記号。77-78 桁が空なら原子名から拾う。"""
    e = line[76:78].strip().upper()
    if e:
        return e
    name = line[12:16]
    if name[:1] != " ":
        return name[:2].strip().upper()
    return name[1:2].strip().upper()


def is_hydrogen(line):
    return element_of(line) in DROP_ELEM


def filter_lines(lines, het):
    """het は 'all' | 'nowat' | 'none'。"""
    out = []
    for l in lines:
        rec = l[:6]
        if rec.startswith(("ATOM", "HETATM")):
            if is_hydrogen(l):
                continue
            if rec.startswith("HETATM"):
                if het == "none":
                    continue
                if het == "nowat" and l[17:20].strip() == "HOH":
                    continue
            out.append(l)
        elif l.startswith(KEEP_OTHER):
            out.append(l)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("srcdir", help="取得したままの PDB が入ったディレクトリ")
    ap.add_argument("outbase", help="3 つのディレクトリを作る場所")
    a = ap.parse_args()

    files = sorted(glob.glob(os.path.join(a.srcdir, "*.pdb")))
    # macOS の付随ファイル（._*）は glob が拾わないが、念のため除く
    files = [f for f in files if not os.path.basename(f).startswith("._")]
    if not files:
        sys.exit(f"{a.srcdir} に *.pdb が無い")

    dirs = {}
    for name, _ in ARMS:
        d = os.path.join(a.outbase, name)
        os.makedirs(d, exist_ok=True)
        dirs[name] = d

    n = 0
    for f in files:
        with open(f, errors="replace") as fh:
            lines = fh.readlines()
        for name, het in ARMS:
            out = filter_lines(lines, het)
            with open(os.path.join(dirs[name], os.path.basename(f)), "w") as fh:
                fh.writelines(out)
        n += 1

    print(f"{n} 構造を 3 条件に展開した")
    for name, _ in ARMS:
        print(f"  {dirs[name]}")


if __name__ == "__main__":
    main()
