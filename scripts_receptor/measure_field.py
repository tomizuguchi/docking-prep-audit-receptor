#!/usr/bin/env python3
"""Measure how much the difference between preparation tools changes the receptor surface electrostatic potential.

  python3 measure_field.py <input pdbdir> <charge dirA> <charge dirB> [--out field_diff.csv]
  python3 measure_field.py --selftest          (cross-check against the analytic solution)

The outputs of the two preparation tools have different atom sets (the Meeko pdbqt merges non-polar
hydrogens, the PDB2PQR pqr writes every hydrogen). "Compare the charge of the same atom" is therefore
impossible. Instead, compare the field itself, which is what each of the two tools ends up showing to
the ligand.

  - The point set is built from the heavy atoms of the input PDB. It does not depend on the output of
    either tool, so a tool that deleted a residue is evaluated at the same points (the deletion shows
    up as a difference in the field at those points).
  - The field is computed using all of the charges in each file, because the charge set the tool
    actually wrote out (including the hydrogens it added) is what is passed on to docking.

Method
  contact surface  Place a sphere of radius (vdW radius + 1.7 A) around each heavy atom and scatter
                   140 points/atom over it with a Fibonacci sphere. Points that fall inside 0.92 times
                   the same radius of another atom are discarded. This approximates the positions
                   reached by the nucleus of an approaching heavy atom.
  vdW radius       C 1.70 / N 1.55 / O 1.52 / S 1.80 / B 1.92 / H 1.20, 2.00 otherwise
  potential        V = 332.0637 x sum(q_i/r_i). Units are kcal/mol per unit charge. r < 0.5 A is clamped to 0.5.
  output           The RMS and the maximum of the difference in V between the two charge sets. The
                   sum of q of each set is reported alongside.

Readable charge files
  .pdbqt   charge in columns 71-76, atom type in columns 78-79
  .pqr     coordinates are at the same column positions as in the PDB (there are files where the
           coordinates run together, so splitting on whitespace drops the field count; 9138 lines
           by measurement). Charge and radius are taken by splitting from column 55 onwards.

Implementations differ in how they build a Fibonacci sphere, so the one used here is written out.
  y = 1 - 2(i+0.5)/N, r = sqrt(1-y^2), phi = i pi (3-sqrt(5)), x = r cos phi, z = r sin phi
"""

import argparse
import csv
import glob
import math
import os
import sys

import numpy as np

COUL = 332.0637  # kcal/mol*A/e^2
PROBE = 1.7  # added on to reach the nucleus of an approaching heavy atom
OCCL = 0.92  # points inside this multiple of another atom's radius are discarded
NPTS = 140  # number of points per atom
RMIN = 0.5  # distances closer than this are clamped to this value
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "B": 1.92, "H": 1.20}
VDW_DEFAULT = 2.00
MAXEL = 30_000_000  # upper bound on the element count of the distance matrix allocated at once (about 240 MB)


def chunk_for(n_atoms):
    """Decide how many points to process at once according to the number of atoms.

    Only the 2-dimensional points x atoms matrix is built, so the element count is held down by MAXEL.
    An implementation that builds a 3-dimensional intermediate array (points x atoms x 3) asked for
    7.4 GB and died on a structure with 7671 heavy atoms and 15440 pqr atoms.
    """
    return max(256, int(MAXEL // max(1, n_atoms)))


def _dist2(p, a, pn, an):
    """Build |p-a|^2 while staying 2-dimensional. r^2 = |p|^2 + |a|^2 - 2p.a."""
    d2 = pn[:, None] + an[None, :]
    d2 -= 2.0 * (p @ a.T)
    np.maximum(d2, 0.0, out=d2)  # rounding can make this negative
    return d2


def radius(el):
    return VDW.get((el or "").strip().upper()[:1], VDW_DEFAULT)


# ---------- reading ----------
def read_pdb_heavy(path):
    """Heavy atoms of the input PDB. Returns coordinates and elements."""
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
    """Coordinates are taken by column position. Charge and radius are taken by splitting from column 55 onwards."""
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
    raise SystemExit(f"unreadable extension: {path}")


# ---------- contact surface ----------
def fibonacci(n):
    i = np.arange(n, dtype=float)
    y = 1.0 - 2.0 * (i + 0.5) / n
    r = np.sqrt(np.clip(1.0 - y * y, 0.0, None))
    phi = i * math.pi * (3.0 - math.sqrt(5.0))
    return np.column_stack((r * np.cos(phi), y, r * np.sin(phi)))


def contact_points(xyz, els, npts=NPTS, probe=PROBE, occl=OCCL):
    """Scatter points over the sphere of each atom and discard the points buried in other atoms."""
    rad = np.asarray([radius(e) + probe for e in els], float)
    unit = fibonacci(npts)
    pts = (xyz[:, None, :] + rad[:, None, None] * unit[None, :, :]).reshape(-1, 3)
    cut = (occl * rad) ** 2  # per-atom cutoff (squared)
    keep = np.ones(len(pts), bool)
    an = (xyz**2).sum(1)
    ch = chunk_for(len(xyz))
    for s in range(0, len(pts), ch):
        p = pts[s : s + ch]
        d2 = _dist2(p, xyz, (p**2).sum(1), an)
        keep[s : s + ch] = ~(d2 < cut[None, :]).any(1)
    return pts[keep]


# ---------- potential ----------
def potential(pts, xyz, q):
    """V = COUL sum(q_i/r_i) at each point."""
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


# ---------- cross-check ----------
def selftest():
    ok = True

    def chk(name, got, want, tol=1e-9):
        nonlocal ok
        good = abs(got - want) <= tol
        ok &= good
        print(f"  {'OK ' if good else 'NG '} {name}: {got:.6f} (expected {want:.6f})")

    # 1 point charge. Place a single point and compare with the analytic solution
    pts = np.array([[10.0, 0.0, 0.0]])
    chk("point charge +1 from 10 A away", potential(pts, np.zeros((1, 3)), np.array([1.0]))[0], COUL / 10.0)
    # clamping
    chk(
        "r=0.1 is clamped to 0.5",
        potential(np.array([[0.1, 0.0, 0.0]]), np.zeros((1, 3)), np.array([1.0]))[0],
        COUL / 0.5,
    )
    # with the same charge set the difference is 0
    xyz = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
    q = np.array([0.3, -0.3])
    v = potential(pts, xyz, q)
    chk("difference for the same charge set", float(np.abs(v - potential(pts, xyz, q)).max()), 0.0)
    # adding a monopole increases the far field by COUL*dq/r
    v2 = potential(np.array([[1000.0, 0.0, 0.0]]), xyz, q + np.array([1.0, 0.0]))
    v1 = potential(np.array([[1000.0, 0.0, 0.0]]), xyz, q)
    chk("adding +1, at 1000 A away", float(v2[0] - v1[0]), COUL / 1000.0, tol=1e-4)
    # point set: for a single isolated atom all 140 points survive and the radius is vdW+probe
    p1 = contact_points(np.zeros((1, 3)), ["C"])
    chk("point count for an isolated atom", float(len(p1)), float(NPTS))
    chk("radius of those points", float(np.linalg.norm(p1, axis=1).max()), VDW["C"] + PROBE, tol=1e-9)
    # occlusion: with 2 atoms at the same position, neither falls inside 0.92 times the other's
    # radius (the radii are equal). Bringing them close together reduces the count
    p2 = contact_points(np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]), ["C", "C"])
    print(f"  -- surviving points for 2 atoms 1 A apart {len(p2)} / {2 * NPTS} (check that the count drops)")
    ok &= len(p2) < 2 * NPTS
    print("\ncross-check: " + ("all passed" if ok else "failures present"))
    return 0 if ok else 1


# ---------- main ----------
def find(d, pid):
    for ext in (".pdbqt", ".pqr"):
        p = os.path.join(d, pid + ext)
        if os.path.exists(p) and os.path.getsize(p) > 0:
            return p
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdbdir", nargs="?", help="directory of input PDBs (heavy atoms only)")
    ap.add_argument("dirA", nargs="?", help="directory of charge set A")
    ap.add_argument("dirB", nargs="?", help="directory of charge set B")
    ap.add_argument("--out", default="field_diff.csv")
    ap.add_argument("--resume", action="store_true", help="skip structures already in the output and append")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()

    if a.selftest:
        sys.exit(selftest())
    if not (a.pdbdir and a.dirA and a.dirB):
        ap.error("pdbdir dirA dirB are required (or --selftest)")

    files = sorted(glob.glob(os.path.join(a.pdbdir, "*.pdb")))
    files = [f for f in files if not os.path.basename(f).startswith("._")]
    if not files:
        sys.exit(f"no *.pdb in {a.pdbdir}")

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
        print(f"{len(done)} structures in the existing output. Skip these and append")

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
                print(f"[{i:3}/{len(files)}] {pid} one of the two is missing, so skip", flush=True)
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
                f"[{i:3}/{len(files)}] {pid} points {len(pts):6d}  "
                f"dV rms {rec['rms_dv']:8.3f}  max {rec['max_abs_dv']:9.3f}  "
                f"sum_q {rec['sum_q_a']:+.2f} / {rec['sum_q_b']:+.2f}",
                flush=True,
            )

    if rows:
        r = np.array([float(x["rms_dv"]) for x in rows])
        m = np.array([float(x["max_abs_dv"]) for x in rows])
        print(f"\n{len(rows)} structures")
        print(
            f"  RMS of dV: median {np.median(r):.2f}, quartiles {np.percentile(r, 25):.2f}"
            f"-{np.percentile(r, 75):.2f}, range {r.min():.2f}-{r.max():.2f} kcal/mol/e"
        )
        print(f"  maximum of dV: median {np.median(m):.2f}, range {m.min():.2f}-{m.max():.2f}")
    print(f"\n-> {a.out}")


if __name__ == "__main__":
    main()
