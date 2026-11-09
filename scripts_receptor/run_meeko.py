#!/usr/bin/env python3
"""Meeko の mk_prepare_receptor を標本に通し、落ち方を分類する。

  python3 run_meeko.py <pdbdir> <outdir> [--rerun | --skip-existing]
                       [-- 追加オプション...]

既定は Meeko 自身が「バッチ処理向け」として案内する組み合わせ:
  --default_altloc A -a
これで通らない構造は、人が中を見て --wanted_altloc を決めない限り
パイプラインに乗らない。何本がそうなるかを数えるのが目的。

成否の判定
  mk_prepare_receptor は失敗しても rc=0 を返すことがある。そのため returncode は
  列に残すが判定には使わない。

  ファイルが存在するかどうかでも足りない。テンプレート照合に失敗した構造でも
  0 バイトの pdbqt を置くことを実測した（1 原子の PDB で確認。log には
  "Template matching failed for" が出て、rc=0、pdbqt は 0 バイト）。
  そのため判定は ATOM 記録が 1 行以上あることとする。ファイルは出たが
  ATOM が無い場合は failure_kind=empty_output として区別する。

出力がすでにある場合（重要）
  判定がファイルの存在に依っているので、前回の出力が残っていると、今回の
  実行が失敗しても成功と記録される。初版はこれを塞いでいなかった。
  出力がある状態では既定で走らせず、どちらかを選ばせる。
    --rerun          その構造の pdbqt を消してから走らせる（測り直す）
    --skip-existing  走らせず、既存ファイルをそのまま測る（途中から再開する）
  --skip-existing で飛ばした行は measured_now=0 として記録する。どの行が
  今回の実行で得られたものかを、出力した CSV から追えるようにするため。

途中で落ちても結果を失わない
  構造ごとに CSV へ追記して flush する。初版は最後に一括で書いていたため、
  タイムアウトや例外でラン全体（数百行）が失われた。タイムアウト自体も
  所見なので、例外で止めず failure_kind に記録して次へ進む。
"""

import argparse
import collections
import csv
import glob
import os
import re
import subprocess
import sys
import time

BASE = ["mk_prepare_receptor.py", "-p"]
DEFAULT_EXTRA = ["--default_altloc", "A", "-a"]
TIMEOUT = 1800

PATTERNS = [
    ("altloc_not_found", r"Requested altlocs not found"),
    ("altloc_unset", r"Either specify an altloc for each"),
    ("template_mismatch", r"did not match any template|does not match templates"),
    ("excess_bonds", r"excess inter-residue bond"),
    ("missing_padding", r"Expected \d+ paddings"),
    ("unknown_residue", r"[Uu]nknown residue|not in the templates"),
    ("gap_or_blunt", r"blunt|chain break|gap"),
]

FIELDS = ["pdb_id", "ok", "measured_now", "preexisting", "returncode", "n_atom", "failure_kind", "secs"]


def classify(txt):
    hits = [name for name, pat in PATTERNS if re.search(pat, txt)]
    return "|".join(hits) if hits else "other"


def count_atoms(path):
    """pdbqt の ATOM 記録の行数。ファイルが無ければ 0。"""
    if not os.path.exists(path):
        return 0
    with open(path, errors="replace") as fh:
        return sum(1 for l in fh if l.startswith("ATOM"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdbdir")
    ap.add_argument("outdir")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--rerun", action="store_true", help="既存の pdbqt を消してから走らせる")
    g.add_argument("--skip-existing", action="store_true", help="既存の pdbqt があれば走らせず、それを測る")
    ap.add_argument(
        "extra", nargs="*", help=f"mk_prepare_receptor への追加オプション（既定 {' '.join(DEFAULT_EXTRA)}）"
    )
    a = ap.parse_args()
    extra = a.extra or DEFAULT_EXTRA

    os.makedirs(a.outdir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(a.pdbdir, "*.pdb")))
    if not files:
        sys.exit(f"{a.pdbdir} に *.pdb が無い")

    # 先に出力の有無を数える。判定がファイルの存在に依るので、黙って走らせない。
    existing = [
        os.path.basename(f)[:4].upper()
        for f in files
        if os.path.exists(os.path.join(a.outdir, os.path.basename(f)[:4].upper() + ".pdbqt"))
    ]
    if existing and not (a.rerun or a.skip_existing):
        sys.exit(
            f"{a.outdir} に pdbqt が {len(existing)} 本ある"
            f"（{', '.join(existing[:5])}{' ...' if len(existing) > 5 else ''}）。\n"
            f"成否は pdbqt の ATOM 記録で判定するので、残っていると今回の実行が"
            f"失敗しても前回のファイルで成功と記録される。\n"
            f"  --rerun          消してから測り直す\n"
            f"  --skip-existing  走らせず既存のものを測る（measured_now=0 で記録）"
        )
    if existing:
        print(
            f"{a.outdir} に pdbqt が {len(existing)} 本ある -> "
            f"{'消して測り直す' if a.rerun else '走らせず既存のものを測る'}"
        )

    csv_path = os.path.join(a.outdir, "meeko.csv")
    rows = []
    with open(csv_path, "w", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=FIELDS)
        writer.writeheader()
        for i, f in enumerate(files, 1):
            pid = os.path.basename(f)[:4].upper()
            base = os.path.join(a.outdir, pid)
            pdbqt = base + ".pdbqt"
            pre = os.path.exists(pdbqt)
            t0 = time.time()

            if pre and a.skip_existing:
                na = count_atoms(pdbqt)
                rec = dict(
                    pdb_id=pid,
                    ok=int(na > 0),
                    measured_now=0,
                    preexisting=1,
                    returncode="",
                    n_atom=na,
                    failure_kind="" if na else "empty_output",
                    secs=0,
                )
                print(f"[{i:2}/{len(files)}] {pid} 既存（走らせていない）{'' if na else ' ATOM が無い'}", flush=True)
            else:
                if pre:
                    # 消してから走らせる。これをしないと、今回失敗しても
                    # 前回のファイルが残って成功と記録される。
                    os.remove(pdbqt)
                try:
                    p = subprocess.run(
                        BASE + ["--read_pdb", f, "-o", base] + extra, capture_output=True, text=True, timeout=TIMEOUT
                    )
                    txt = (p.stdout or "") + "\n" + (p.stderr or "")
                    rc = p.returncode
                except subprocess.TimeoutExpired as e:
                    txt = (
                        (e.stdout or "") if isinstance(e.stdout, str) else (e.stdout or b"").decode("utf-8", "replace")
                    )
                    txt += f"\n--- timeout after {TIMEOUT}s ---"
                    rc = "timeout"
                except FileNotFoundError:
                    sys.exit(f"{BASE[0]} が見つからない")
                with open(base + ".log", "w") as fh:
                    fh.write(txt)
                # 存在だけでは足りない。0 バイトの pdbqt を置いて rc=0 で
                # 終わることがあるので、ATOM が 1 行以上あることを見る。
                na = count_atoms(pdbqt)
                got = na > 0
                if got:
                    kind = ""
                elif rc == "timeout":
                    kind = "timeout"
                elif os.path.exists(pdbqt):
                    kind = "empty_output"
                else:
                    kind = classify(txt)
                rec = dict(
                    pdb_id=pid,
                    ok=int(got),
                    measured_now=1,
                    preexisting=int(pre),
                    returncode=rc,
                    n_atom=na,
                    failure_kind=kind,
                    secs=round(time.time() - t0, 1),
                )
                print(
                    f"[{i:2}/{len(files)}] {pid} {'成功' if got else '失敗 ' + kind} ({time.time() - t0:.0f}s)",
                    flush=True,
                )

            rows.append(rec)
            writer.writerow(rec)
            fout.flush()

    ok = sum(r["ok"] for r in rows)
    now = sum(r["measured_now"] for r in rows)
    print(f"\n成功 {ok} / 失敗 {len(rows) - ok}  （オプション: {' '.join(extra)}）")
    print(f"このランで実際に走らせた構造 {now}/{len(rows)}")
    for k, n in collections.Counter(r["failure_kind"] for r in rows if not r["ok"]).most_common():
        print(f"  {k}: {n} 本")
    # returncode は判定に使っていないが、rc=0 のまま失敗した件数は所見として出す
    liars = sum(1 for r in rows if not r["ok"] and r["returncode"] == 0)
    if liars:
        print(f"\n失敗したが returncode が 0 だった構造 {liars} 本（rc で判定できないことの実測）")
    print(f"\n-> {csv_path}")


if __name__ == "__main__":
    main()
