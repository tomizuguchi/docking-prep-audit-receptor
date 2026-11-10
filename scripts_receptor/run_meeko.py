#!/usr/bin/env python3
"""Run Meeko's mk_prepare_receptor over the sample and classify how it fails.

  python3 run_meeko.py <pdbdir> <outdir> [--rerun | --skip-existing]
                       [-- extra options...]

The default here is
  --default_altloc A -a
These two were made the default so that 100 structures go through without a person intervening.
The pair comes from Meeko itself: when mk_prepare_receptor.py stops on a residue it cannot match or
on alternate locations, its error message ends with
  "(for batch processing) Use -a/--allow_bad_res to automatically remove residues
   that do not match templates, and --default_altloc to set a default altloc variant.
   Use these at your own risk."
(Meeko 0.8.0; the text is in the log of every run that fails that way.) In 0.8.0 --help lists the
flag as -x/--delete_bad_res; -a and --allow_bad_res are still accepted.
A structure that does not pass with these does not get onto the pipeline unless a person looks
inside it and decides --wanted_altloc. The purpose is to count how many structures are like that.

Deciding success or failure
  mk_prepare_receptor can return rc=0 even when it fails. returncode is therefore kept as a
  column but is not used for the decision.

  Whether the file exists is not enough either. It was measured that a structure whose template
  matching failed still gets a 0-byte pdbqt (confirmed with a 1-atom PDB; the log shows
  "Template matching failed for", rc=0, and the pdbqt is 0 bytes).
  The decision is therefore that there is at least one ATOM record. When the file appeared but
  has no ATOM, it is distinguished as failure_kind=empty_output.

When output already exists (important)
  Because the decision depends on the existence of the file, output left from the previous run makes
  this run be recorded as a success even when it fails. The first version did not block this.
  With output present, nothing runs by default; one of the two has to be chosen.
    --rerun          delete that structure's pdbqt and then run (measure again)
    --skip-existing  do not run and measure the existing file as it is (resume part-way)
  Rows skipped by --skip-existing are recorded as measured_now=0. This is so that which rows came
  from this run can be followed from the CSV that is written.

Results are not lost when it dies part-way
  The CSV is appended and flushed per structure. The first version wrote everything at the end, so
  a timeout or an exception lost the whole run (several hundred rows). A timeout is itself a finding,
  so it is recorded in failure_kind and the run goes on to the next instead of stopping on an exception.
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

# All patterns that match are joined with "|" in failure_kind. A failed download of a residue
# definition from rcsb.org also matches unknown_residue, so those rows read
# "template_download|unknown_residue"; template_download is what tells them from a residue that
# is unknown with the network available.
PATTERNS = [
    ("template_download", r"Could not download CIF|urlopen error"),
    ("residue_build", r"unable to build rdkit mol"),
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
    """Number of ATOM record lines in the pdbqt. 0 when the file does not exist."""
    if not os.path.exists(path):
        return 0
    with open(path, errors="replace") as fh:
        return sum(1 for l in fh if l.startswith("ATOM"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdbdir")
    ap.add_argument("outdir")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--rerun", action="store_true", help="delete the existing pdbqt and then run")
    g.add_argument("--skip-existing", action="store_true", help="if a pdbqt exists, do not run and measure it")
    ap.add_argument(
        "extra", nargs="*", help=f"extra options for mk_prepare_receptor (default {' '.join(DEFAULT_EXTRA)})"
    )
    a = ap.parse_args()
    extra = a.extra or DEFAULT_EXTRA

    os.makedirs(a.outdir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(a.pdbdir, "*.pdb")))
    if not files:
        sys.exit(f"no *.pdb in {a.pdbdir}")

    # Count the existing output first. The decision depends on the existence of the file, so do not
    # run silently.
    existing = [
        os.path.basename(f)[:4].upper()
        for f in files
        if os.path.exists(os.path.join(a.outdir, os.path.basename(f)[:4].upper() + ".pdbqt"))
    ]
    if existing and not (a.rerun or a.skip_existing):
        sys.exit(
            f"{a.outdir} holds {len(existing)} pdbqt files"
            f" ({', '.join(existing[:5])}{' ...' if len(existing) > 5 else ''}).\n"
            f"Success is decided from the ATOM records of the pdbqt, so if they are left in place this run"
            f" is recorded as a success from the previous files even when it fails.\n"
            f"  --rerun          delete them and measure again\n"
            f"  --skip-existing  do not run and measure the existing ones (recorded as measured_now=0)"
        )
    if existing:
        print(
            f"{a.outdir} holds {len(existing)} pdbqt files -> "
            f"{'delete them and measure again' if a.rerun else 'do not run and measure the existing ones'}"
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
                print(f"[{i:2}/{len(files)}] {pid} existing (not run){'' if na else ' no ATOM'}", flush=True)
            else:
                if pre:
                    # Delete, then run. Without this, even when this run fails, the previous
                    # file is left in place and recorded as a success.
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
                    sys.exit(f"{BASE[0]} not found")
                with open(base + ".log", "w") as fh:
                    fh.write(txt)
                # Existence alone is not enough. It can leave a 0-byte pdbqt and finish with rc=0,
                # so look for at least one ATOM line.
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
                    f"[{i:2}/{len(files)}] {pid} {'success' if got else 'failure ' + kind} ({time.time() - t0:.0f}s)",
                    flush=True,
                )

            rows.append(rec)
            writer.writerow(rec)
            fout.flush()

    ok = sum(r["ok"] for r in rows)
    now = sum(r["measured_now"] for r in rows)
    print(f"\nsuccess {ok} / failure {len(rows) - ok}  (options: {' '.join(extra)})")
    print(f"structures actually run in this run {now}/{len(rows)}")
    for k, n in collections.Counter(r["failure_kind"] for r in rows if not r["ok"]).most_common():
        print(f"  {k}: {n} structures")
    # returncode is not used for the decision, but the number that failed with rc=0 is reported as a finding
    liars = sum(1 for r in rows if not r["ok"] and r["returncode"] == 0)
    if liars:
        print(f"\n{liars} structures failed but had returncode 0 (measured evidence that rc cannot decide this)")
    print(f"\n-> {csv_path}")


if __name__ == "__main__":
    main()
