#!/usr/bin/env python3
"""準備ツールの違いが、受容体表面の静電ポテンシャルをどれだけ変えるかを測る。

  python3 measure_field.py <入力pdbdir> <電荷dirA> <電荷dirB> [--out field_diff.csv]
  python3 measure_field.py --selftest          （解析解との照合）

二つの準備ツールの出力は原子集合が違う（Meeko の pdbqt は非極性水素を併合し、
PDB2PQR の pqr は水素を全部書く）。そのため「同じ原子の電荷を比べる」ことはできない。
かわりに、どちらのツールもリガンドに見せることになる場そのものを比べる。

  - 点集合は入力 PDB の重原子から作る。どちらのツールの出力にも依存しないので、
    残基を削除したツールでも同じ点で評価される（削除はその点での場の差として出る）。
  - 場は各ファイルの電荷を全部使って計算する。そのツールが実際に書き出した
    電荷セット（追加した水素を含む）が、ドッキングに渡るものだから。

方法
  接触面   各重原子を中心に半径 (vdW 半径 + 1.7 Å) の球を置き、Fibonacci 球で
           140 点/原子を散らす。他の原子の同じ半径の 0.92 倍より内側に入る点は捨てる。
           接触してくる重原子の核が来る位置の近似。
  vdW 半径 C 1.70 / N 1.55 / O 1.52 / S 1.80 / B 1.92 / H 1.20、その他 2.00
  電位     V = 332.0637 × Σ(qᵢ/rᵢ)。単位は kcal/mol 毎単位電荷。r < 0.5 Å は 0.5 で打ち切る。
  出力     二つの電荷セットの V の差の RMS と最大。各セットの Σq も併記する。

読める電荷ファイル
  .pdbqt   電荷は 71-76 桁、原子型は 78-79 桁
  .pqr     座標は PDB と同じ列位置（空白区切りにすると座標が繋がって要素数が落ちる
           ファイルがある。実測で 9138 行）。電荷と半径は 55 桁以降を分割して取る。

Fibonacci 球の作り方は実装によって違うので、ここで使ったものを書いておく。
  y = 1 - 2(i+0.5)/N、r = sqrt(1-y²)、φ = i π (3-√5)、x = r cos φ、z = r sin φ
"""

import argparse
import csv
import glob
import math
import os
import sys

import numpy as np

COUL = 332.0637  # kcal/mol·Å/e²
PROBE = 1.7  # 接触してくる重原子の核までの上乗せ
OCCL = 0.92  # 他原子の半径のこの倍より内側の点は捨てる
NPTS = 140  # 原子あたりの点数
RMIN = 0.5  # これより近い距離はこの値で打ち切る
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "B": 1.92, "H": 1.20}
VDW_DEFAULT = 2.00
MAXEL = 30_000_000  # 一度に確保する距離行列の要素数の上限（約 240 MB）


def chunk_for(n_atoms):
    """原子数に応じて一度に処理する点の数を決める。

    点×原子の 2 次元行列しか作らないので、要素数を MAXEL で抑える。
    3 次元の中間配列（点×原子×3）を作る実装では、重原子 7671・pqr 15440 原子の
    構造で 7.4 GB を要求して落ちた。
    """
    return max(256, int(MAXEL // max(1, n_atoms)))


def _dist2(p, a, pn, an):
    """|p-a|² を 2 次元のまま作る。r² = |p|² + |a|² - 2p·a。"""
    d2 = pn[:, None] + an[None, :]
    d2 -= 2.0 * (p @ a.T)
    np.maximum(d2, 0.0, out=d2)  # 丸めで負になることがある
    return d2


def radius(el):
    return VDW.get((el or "").strip().upper()[:1], VDW_DEFAULT)


# ---------- 読み込み ----------
def read_pdb_heavy(path):
    """入力 PDB の重原子。座標と元素を返す。"""
    xyz, el = [], []
    for l in open(path, errors="replace"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        e = l[76:78].strip() or l[12:16].strip()[:1]
        if e.upper().startswith(("H", "D")):
            continue
        xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
        el.append(e)
    return np.asarray(xyz, float), el


def read_pdbqt(path):
    xyz, q = [], []
    for l in open(path, errors="replace"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        try:
            c = float(l[70:76])
        except ValueError:
            tail = l[54:].split()
            if not tail:
                continue
            c = float(tail[-2]) if len(tail) >= 2 else float(tail[-1])
        xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
        q.append(c)
    return np.asarray(xyz, float), np.asarray(q, float)


def read_pqr(path):
    """座標は列位置で取る。電荷と半径は 55 桁以降を分割して取る。"""
    xyz, q = [], []
    for l in open(path, errors="replace"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        tail = l[54:].split()
        if len(tail) < 2:
            continue
        xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
        q.append(float(tail[0]))
    return np.asarray(xyz, float), np.asarray(q, float)


def read_charges(path):
    if path.endswith(".pdbqt"):
        return read_pdbqt(path)
    if path.endswith(".pqr"):
        return read_pqr(path)
    raise SystemExit(f"読めない拡張子: {path}")


# ---------- 接触面 ----------
def fibonacci(n):
    i = np.arange(n, dtype=float)
    y = 1.0 - 2.0 * (i + 0.5) / n
    r = np.sqrt(np.clip(1.0 - y * y, 0.0, None))
    phi = i * math.pi * (3.0 - math.sqrt(5.0))
    return np.column_stack((r * np.cos(phi), y, r * np.sin(phi)))


def contact_points(xyz, els, npts=NPTS, probe=PROBE, occl=OCCL):
    """各原子の球面に点を散らし、他原子に埋もれる点を捨てる。"""
    rad = np.asarray([radius(e) + probe for e in els], float)
    unit = fibonacci(npts)
    pts = (xyz[:, None, :] + rad[:, None, None] * unit[None, :, :]).reshape(-1, 3)
    cut = (occl * rad) ** 2  # 原子ごとの打ち切り（二乗）
    keep = np.ones(len(pts), bool)
    an = (xyz**2).sum(1)
    ch = chunk_for(len(xyz))
    for s in range(0, len(pts), ch):
        p = pts[s : s + ch]
        d2 = _dist2(p, xyz, (p**2).sum(1), an)
        keep[s : s + ch] = ~(d2 < cut[None, :]).any(1)
    return pts[keep]


# ---------- 電位 ----------
def potential(pts, xyz, q):
    """各点での V = COUL Σ qᵢ/rᵢ。"""
    nz = q != 0.0
    xyz, q = xyz[nz], q[nz]
    out = np.empty(len(pts), float)
    an = (xyz**2).sum(1)
    ch = chunk_for(len(xyz))
    for s in range(0, len(pts), ch):
        p = pts[s : s + ch]
        r = np.sqrt(_dist2(p, xyz, (p**2).sum(1), an))
        np.maximum(r, RMIN, out=r)
        out[s : s + ch] = COUL * (q[None, :] / r).sum(1)
    return out


# ---------- 検算 ----------
def selftest():
    ok = True

    def chk(name, got, want, tol=1e-9):
        nonlocal ok
        good = abs(got - want) <= tol
        ok &= good
        print(f"  {'OK ' if good else 'NG '} {name}: {got:.6f}（期待 {want:.6f}）")

    # 1 点電荷。1 点だけを置いて解析解と比べる
    pts = np.array([[10.0, 0.0, 0.0]])
    chk("点電荷 +1 を 10 Å 先から", potential(pts, np.zeros((1, 3)), np.array([1.0]))[0], COUL / 10.0)
    # 打ち切り
    chk(
        "r=0.1 は 0.5 に打ち切る",
        potential(np.array([[0.1, 0.0, 0.0]]), np.zeros((1, 3)), np.array([1.0]))[0],
        COUL / 0.5,
    )
    # 同じ電荷セットなら差は 0
    xyz = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
    q = np.array([0.3, -0.3])
    v = potential(pts, xyz, q)
    chk("同じ電荷セットの差", float(np.abs(v - potential(pts, xyz, q)).max()), 0.0)
    # 単極子を足すと遠方で COUL·Δq/r 増える
    v2 = potential(np.array([[1000.0, 0.0, 0.0]]), xyz, q + np.array([1.0, 0.0]))
    v1 = potential(np.array([[1000.0, 0.0, 0.0]]), xyz, q)
    chk("+1 を足すと 1000 Å 先で", float(v2[0] - v1[0]), COUL / 1000.0, tol=1e-4)
    # 点集合: 孤立した 1 原子なら 140 点すべて残り、半径は vdW+probe
    p1 = contact_points(np.zeros((1, 3)), ["C"])
    chk("孤立原子の点数", float(len(p1)), float(NPTS))
    chk("その点の半径", float(np.linalg.norm(p1, axis=1).max()), VDW["C"] + PROBE, tol=1e-9)
    # 遮蔽: 同じ位置に 2 原子を置くと、互いの 0.92 倍より内側…には入らない（半径が同じ）。
    # 近接させると減る
    p2 = contact_points(np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]), ["C", "C"])
    print(f"  -- 1 Å 離した 2 原子の残存点 {len(p2)} / {2 * NPTS}（減ることの確認）")
    ok &= len(p2) < 2 * NPTS
    print("\n検算: " + ("すべて通った" if ok else "失敗あり"))
    return 0 if ok else 1


# ---------- 本体 ----------
def find(d, pid):
    for ext in (".pdbqt", ".pqr"):
        p = os.path.join(d, pid + ext)
        if os.path.exists(p) and os.path.getsize(p) > 0:
            return p
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdbdir", nargs="?", help="入力 PDB（重原子のみ）のディレクトリ")
    ap.add_argument("dirA", nargs="?", help="電荷セット A のディレクトリ")
    ap.add_argument("dirB", nargs="?", help="電荷セット B のディレクトリ")
    ap.add_argument("--out", default="field_diff.csv")
    ap.add_argument("--resume", action="store_true", help="既存の出力にある構造を飛ばして追記する")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()

    if a.selftest:
        sys.exit(selftest())
    if not (a.pdbdir and a.dirA and a.dirB):
        ap.error("pdbdir dirA dirB が必要（または --selftest）")

    files = sorted(glob.glob(os.path.join(a.pdbdir, "*.pdb")))
    files = [f for f in files if not os.path.basename(f).startswith("._")]
    if not files:
        sys.exit(f"{a.pdbdir} に *.pdb が無い")

    cols = [
        "pdb_id",
        "n_heavy",
        "n_points",
        "n_atoms_a",
        "n_atoms_b",
        "sum_q_a",
        "sum_q_b",
        "d_sum_q",
        "rms_dv",
        "max_abs_dv",
        "mean_dv",
        "rms_v_a",
        "rms_v_b",
    ]
    rows, done = [], set()
    if a.resume and os.path.exists(a.out):
        with open(a.out, newline="") as fh:
            for r in csv.DictReader(fh):
                rows.append({k: r[k] for k in cols})
                done.add(r["pdb_id"].upper())
        print(f"既存の出力に {len(done)} 構造。これを飛ばして追記する")

    with open(a.out, "a" if done else "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        if not done:
            w.writeheader()
        for i, f in enumerate(files, 1):
            pid = os.path.basename(f)[:4].upper()
            if pid in done:
                continue
            pa, pb = find(a.dirA, pid), find(a.dirB, pid)
            if not (pa and pb):
                print(f"[{i:3}/{len(files)}] {pid} 片方が無いので飛ばす", flush=True)
                continue
            xyz, els = read_pdb_heavy(f)
            pts = contact_points(xyz, els)
            xa, qa = read_charges(pa)
            xb, qb = read_charges(pb)
            va, vb = potential(pts, xa, qa), potential(pts, xb, qb)
            d = va - vb
            rec = dict(
                pdb_id=pid,
                n_heavy=len(xyz),
                n_points=len(pts),
                n_atoms_a=len(xa),
                n_atoms_b=len(xb),
                sum_q_a=round(float(qa.sum()), 3),
                sum_q_b=round(float(qb.sum()), 3),
                d_sum_q=round(float(qa.sum() - qb.sum()), 3),
                rms_dv=round(float(np.sqrt((d**2).mean())), 3),
                max_abs_dv=round(float(np.abs(d).max()), 3),
                mean_dv=round(float(d.mean()), 3),
                rms_v_a=round(float(np.sqrt((va**2).mean())), 3),
                rms_v_b=round(float(np.sqrt((vb**2).mean())), 3),
            )
            rows.append(rec)
            w.writerow(rec)
            fh.flush()
            print(
                f"[{i:3}/{len(files)}] {pid} 点 {len(pts):6d}  "
                f"ΔV rms {rec['rms_dv']:8.3f}  max {rec['max_abs_dv']:9.3f}  "
                f"Σq {rec['sum_q_a']:+.2f} / {rec['sum_q_b']:+.2f}",
                flush=True,
            )

    if rows:
        r = np.array([float(x["rms_dv"]) for x in rows])
        m = np.array([float(x["max_abs_dv"]) for x in rows])
        print(f"\n{len(rows)} 構造")
        print(
            f"  ΔV の RMS: 中央値 {np.median(r):.2f}、四分位 {np.percentile(r, 25):.2f}"
            f"–{np.percentile(r, 75):.2f}、範囲 {r.min():.2f}–{r.max():.2f} kcal/mol/e"
        )
        print(f"  ΔV の最大: 中央値 {np.median(m):.2f}、範囲 {m.min():.2f}–{m.max():.2f}")
    print(f"\n-> {a.out}")


if __name__ == "__main__":
    main()
