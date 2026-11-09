#!/usr/bin/env python3
"""マニフェストに載っている構造を、その ID で RCSB から取得する。

  python3 fetch_manifest.py <manifest.csv> <outdir>

`fetch_strata.py` を再実行すると、検索の母集団がその日の RCSB の内容になるので、
同じシードでも標本が変わりうる（`pool_size` 列が当時の母集団の大きさ）。
測定に使ったのと同じ 100 構造を取るには、マニフェストの `pdb_id` 列から
直接取る。こちらを使うこと。

取得後、`data_receptor/full100_sha256.txt` と照合できる。

  cd <outdir の親>
  shasum -a 256 -c data_receptor/full100_sha256.txt

RCSB が同じ ID のファイルを改訂した場合、照合は通らない。そのときは
通らなかった ID を控えること（測定時のファイルと中身が違うという事実になる）。
"""

import argparse
import csv
import os
import sys
import time
import urllib.error
import urllib.request

URL = "https://files.rcsb.org/download/{}.pdb"
RETRY = 3


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("manifest")
    ap.add_argument("outdir")
    ap.add_argument("--col", default="pdb_id", help="ID の列名（既定 pdb_id）")
    a = ap.parse_args()

    with open(a.manifest, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows or a.col not in rows[0]:
        sys.exit(f"{a.manifest} に {a.col} 列が無い")
    ids = []
    for r in rows:
        p = (r[a.col] or "").strip().upper()
        if p and p not in ids:
            ids.append(p)
    os.makedirs(a.outdir, exist_ok=True)

    got, skipped, failed = 0, 0, []
    for i, p in enumerate(ids, 1):
        out = os.path.join(a.outdir, p + ".pdb")
        if os.path.exists(out) and os.path.getsize(out) > 0:
            skipped += 1
            continue
        for k in range(RETRY):
            try:
                with urllib.request.urlopen(URL.format(p), timeout=120) as r:
                    body = r.read()
                with open(out, "wb") as fh:
                    fh.write(body)
                got += 1
                print(f"[{i}/{len(ids)}] {p} {len(body)} bytes", flush=True)
                break
            except Exception as e:
                if k == RETRY - 1:
                    failed.append((p, f"{type(e).__name__}: {e}"))
                    print(f"[{i}/{len(ids)}] {p} 失敗 {type(e).__name__}", flush=True)
                else:
                    time.sleep(2 * (k + 1))

    print(f"\n取得 {got} / 既存を飛ばした {skipped} / 失敗 {len(failed)}（ID {len(ids)} 件）")
    for p, why in failed:
        print(f"  {p}: {why}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
