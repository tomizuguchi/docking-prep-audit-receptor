#!/usr/bin/env python3
"""
Fetch PDB structures stratified by resolution. The random seed and the query are left in the manifest.

  python3 fetch_strata.py list                        check only the count of each stratum (no download)
  python3 fetch_strata.py list  --group               counts for the version deduplicated at 90% sequence identity
  python3 fetch_strata.py fetch <N> <outdir> [seed] [--group]

With the same seed, increasing N leaves the previous sample at the head unchanged
(pilot 8 per stratum is a subset of production 20 per stratum). draw_rank in the manifest is the draw order.

The default is return_type=entry with no deduplication (= a sample of the entry population).
With --group, only the representatives at 90% sequence identity are used (= a sample close to the
protein population). Which of the two was used is recorded in manifest.csv and query.json.
"""

import csv
import json
import os
import random
import sys
import time
import urllib.error
import urllib.request

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
FILES = "https://files.rcsb.org/download/{}.pdb"
STRATA = [
    ("A_le1.5", 0.0, 1.5),
    ("B_1.5-2.0", 1.5, 2.0),
    ("C_2.0-2.5", 2.0, 2.5),
    ("D_2.5-3.0", 2.5, 3.0),
    ("E_ge3.0", 3.0, 9.9),
]


def T(attr, op, val):
    return {"type": "terminal", "service": "text", "parameters": {"attribute": attr, "operator": op, "value": val}}


def query(lo, hi, group=False):
    nodes = [
        T("exptl.method", "exact_match", "X-RAY DIFFRACTION"),
        T("rcsb_entry_info.selected_polymer_entity_types", "exact_match", "Protein (only)"),
        T("rcsb_entry_info.deposited_polymer_entity_instance_count", "equals", 1),
        T(
            "rcsb_entry_info.resolution_combined",
            "range",
            {"from": lo, "include_lower": True, "to": hi, "include_upper": False},
        ),
    ]
    q = {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "polymer_entity" if group else "entry",
        "request_options": {"paginate": {"start": 0, "rows": 10000}, "results_content_type": ["experimental"]},
    }
    if group:
        q["request_options"]["group_by"] = {"aggregation_method": "sequence_identity", "similarity_cutoff": 90}
        q["request_options"]["group_by_return_type"] = "representatives"
    return q


def post(q, label=""):
    req = urllib.request.Request(SEARCH, data=json.dumps(q).encode(), headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=180) as r:
            return json.loads(r.read())
    except urllib.error.HTTPError as e:
        print(f"  !! {label} HTTP {e.code}")
        print("  --- server response ---")
        print("  " + e.read().decode()[:1500].replace("\n", "\n  "))
        print("  --- query sent ---")
        print("  " + json.dumps(q)[:900])
        return None
    except Exception as e:
        print(f"  !! {label} {type(e).__name__}: {e}")
        return None


def ids_for(lo, hi, group, verbose=True):
    """For a stratum above the paginate limit (10000), advance start and collect every entry."""
    PAGE, got, total, start = 10000, set(), None, 0
    while True:
        q = query(lo, hi, group)
        q["request_options"]["paginate"] = {"start": start, "rows": PAGE}
        d = post(q, f"{lo}-{hi} start={start}")
        if d is None:
            return None, 0
        if total is None:
            total = d.get("total_count", 0)
        rs = d.get("result_set", [])
        if not rs:
            break
        got |= {r["identifier"].split("_")[0] for r in rs}
        start += PAGE
        if verbose and total > PAGE:
            print(f"      ... {len(got)}/{total}")
        if start >= total:
            break
        time.sleep(0.5)
    return sorted(got), total


def pools_for(group):
    """Collect the IDs of all strata and return them with entries in 2 or more strata removed from every stratum.

    resolution_combined is an array attribute, and a range search matches when "any element of the
    array is in the range". An entry with several resolution values enters several strata at once, and
    the stratification breaks down. Which of the values is "the resolution of that structure" cannot be
    decided here, so entries that straddle strata are removed from the population. The removed IDs are
    left in excluded.csv.
    """
    pools, totals = {}, {}
    for name, lo, hi in STRATA:
        ids, total = ids_for(lo, hi, group)
        if ids is None:
            return None, None, None
        pools[name], totals[name] = set(ids), total
        time.sleep(1)
    names = [n for n, _, _ in STRATA]
    seen, multi, where = set(), set(), {}
    for n in names:
        multi |= seen & pools[n]
        seen |= pools[n]
    for n in names:
        for pid in pools[n] & multi:
            where.setdefault(pid, []).append(n)
    if multi:
        for a in range(len(names)):
            for b in range(a + 1, len(names)):
                ov = pools[names[a]] & pools[names[b]]
                if ov:
                    print(f"  stratum overlap {names[a]} & {names[b]}: {len(ov)} entries")
        print(
            f"  excluding {len(multi)} entries with several resolution values that straddle strata from all strata"
            f" ({len(multi) / len(seen) * 100:.2f}% of the population of {len(seen)})"
        )
        for n in names:
            pools[n] -= multi
    else:
        print("  no overlap between strata")
    print(f"  effective population {sum(len(p) for p in pools.values())} entries")
    return pools, totals, where


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    group = "--group" in sys.argv
    if not args:
        print(__doc__)
        return
    mode = args[0]
    if mode == "list":
        print(f"dedup: {'representatives at 90% sequence identity only' if group else 'none (entry population)'}")
        pools, totals, where = pools_for(group)
        if pools is None:
            print("\nPaste the server response above and this can be fixed.")
            return
        for name, lo, hi in STRATA:
            print(f"  {name:<10} {lo}-{hi} A : search {totals[name]}, after overlap exclusion {len(pools[name])}")
        print("\nReady to proceed to fetch.")
        return
    if mode != "fetch":
        print(__doc__)
        return
    N, outdir = int(args[1]), args[2]
    seed = int(args[3]) if len(args) > 3 else 20261108
    os.makedirs(outdir, exist_ok=True)
    pools, totals, where = pools_for(group)
    if pools is None:
        print("  fetch failed. Abort")
        return
    with open(os.path.join(outdir, "excluded.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["pdb_id", "reason", "strata"])
        for pid, ns in sorted(where.items()):
            w.writerow([pid, "resolution_combined straddles several strata", "|".join(ns)])
    man, failed = [], []
    for name, lo, hi in STRATA:
        ids = sorted(pools[name])
        # Build the full permutation and take from the head. With the same seed, increasing N
        # leaves the previous sample at the head unchanged (pilot is a subset of production).
        order = random.Random(f"{seed}:{name}").sample(ids, len(ids))
        pick, rank, i, bad = [], {}, 0, 0
        while len(pick) < min(N, len(ids)) and i < len(order):
            pid = order[i]
            i += 1
            dst = os.path.join(outdir, f"{pid}.pdb")
            if not os.path.exists(dst):
                try:
                    urllib.request.urlretrieve(FILES.format(pid), dst)
                    time.sleep(0.2)
                except Exception as e:
                    print(f"    !! {pid}: {type(e).__name__} {e}")
                    bad += 1
                    failed.append({"stratum": name, "pdb_id": pid, "draw_rank": i, "reason": f"{type(e).__name__} {e}"})
                    continue
            # A 404, or an entry that exists only as mmCIF, can leave an HTML or empty file behind
            try:
                head = open(dst, errors="replace").read(200)
                nl = sum(1 for _ in open(dst, errors="replace"))
            except Exception:
                head, nl = "", 0
            if not head.startswith("HEADER") or nl < 100:
                print(f"    !! {pid}: not PDB format ({nl} lines). Excluded")
                os.remove(dst)
                bad += 1
                failed.append(
                    {"stratum": name, "pdb_id": pid, "draw_rank": i, "reason": f"not PDB format ({nl} lines)"}
                )
                continue
            rank[pid] = i
            pick.append(pid)
        print(
            f"  {name:<10} {len(pick)} entries from a population of {len(ids)}"
            + (f" ({bad} could not be fetched and were skipped)" if bad else "")
        )
        for pid in pick:
            man.append(
                {
                    "stratum": name,
                    "res_lo": lo,
                    "res_hi": hi,
                    "pdb_id": pid,
                    "draw_rank": rank[pid],
                    "pool_size": len(ids),
                    "seed": seed,
                    "dedup90": int(group),
                }
            )
    if not man:
        print("  not a single entry could be fetched")
        return
    if failed:
        with open(os.path.join(outdir, "failed.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(failed[0]))
            w.writeheader()
            w.writerows(failed)
        print(
            f"  the {len(failed)} entries that could not be fetched are recorded in failed.csv"
            f" (already refilled from the next rank)"
        )
    with open(os.path.join(outdir, "manifest.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(man[0]))
        w.writeheader()
        w.writerows(man)
    with open(os.path.join(outdir, "query.json"), "w") as fh:
        json.dump(
            {
                "seed": seed,
                "dedup90": group,
                "endpoint": SEARCH,
                "exclusion_rule": (
                    "resolution_combined is an array attribute, so some entries match "
                    "several strata. IDs that appeared in 2 or more strata were excluded from all strata "
                    "(excluded.csv). The population of each stratum is the count after exclusion, "
                    "which is pool_size in manifest.csv."
                ),
                "n_excluded": len(where),
                "strata": [
                    {
                        "name": n,
                        "lo": lo,
                        "hi": hi,
                        "query": query(lo, hi, group),
                        "search_total": totals[n],
                        "pool_after_exclusion": len(pools[n]),
                    }
                    for n, lo, hi in STRATA
                ],
            },
            fh,
            indent=1,
        )
    print(f"\n{len(man)} entries -> {outdir} (manifest.csv and query.json included)")


if __name__ == "__main__":
    main()
