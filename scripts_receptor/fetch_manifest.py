#!/usr/bin/env python3
"""Download the structures listed in the manifest from RCSB by their IDs.

  python3 fetch_manifest.py <manifest.csv> <outdir>

Re-running `fetch_strata.py` makes the search population the contents of RCSB on that day, so the
sample can change even with the same seed (the `pool_size` column is the population size at the time).
To get the same 100 structures that were used for the measurement, take them directly from the
`pdb_id` column of the manifest. Use this script.

After downloading, the files can be checked against `data_receptor/full100_sha256.txt`.

  cd <parent of outdir>
  shasum -a 256 -c data_receptor/full100_sha256.txt

If RCSB has revised the file for a given ID, the check does not pass. In that case, note down the IDs
that did not pass (it is the fact that the contents differ from the files used for the measurement).
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
    ap.add_argument("--col", default="pdb_id", help="name of the ID column (default pdb_id)")
    a = ap.parse_args()

    with open(a.manifest, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows or a.col not in rows[0]:
        sys.exit(f"{a.manifest} has no {a.col} column")
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
                    print(f"[{i}/{len(ids)}] {p} failed {type(e).__name__}", flush=True)
                else:
                    time.sleep(2 * (k + 1))

    print(f"\nfetched {got} / skipped existing {skipped} / failed {len(failed)} ({len(ids)} IDs)")
    for p, why in failed:
        print(f"  {p}: {why}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
