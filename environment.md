# Run environment

## Versions

| Tool | Version | Where this comes from |
|---|---|---|
| PDB2PQR | 3.8.0 | `pdb2pqr --version` |
| propka | 3.5.1 | called internally by PDB2PQR |
| Meeko | 0.8.0 | `python3 -c "import meeko; print(meeko.__version__)"`; see below |
| RDKit | 2026.09.1 | used by Meeko |
| gemmi | 0.7.5 | used by Meeko |
| Python | 3.13.18 | the interpreter path in the logs is `python3.13` |
| numpy | 2.5.3 | analysis only |

The Meeko version in this deposit is 0.8.0. It differs from the 0.7.x series in two
ways.

- The contents of the built-in templates differ.
- `-a` / `--allow_bad_res` was renamed to `-x` / `--delete_bad_res` in 0.8.0. The
  old names are still accepted but do not appear in `--help`. `run_meeko.py`
  passes the old `-a`, which 0.8.0 accepts.

Meeko does not write its version into its log. That every prepared file came from
0.8.0 was established by running the whole preparation again under 0.8.0 and
comparing checksums; see "What was confirmed by re-running" below.

Run under the 0.7.x series, the same scripts may not produce identical output.
Match the version before checking against `data_receptor/prepared_sha256.txt`.

## Two machines were involved

| Stage | Where it ran | Reaches rcsb.org |
|---|---|---|
| Fetching the 100 structures | the machine the structures were collected on | yes |
| Every Meeko and PDB2PQR condition | a machine with restricted network access | no |

The PDB files in `full100/` were fetched on the machine with network access, which
is why that directory carries macOS sidecar files (`._*`). The preparation and the
measurements that followed ran in an environment that cannot reach rcsb.org. HTTP
connections are refused by a proxy with 403.

```
$ python3 -c "import urllib.request; urllib.request.urlopen('https://files.rcsb.org/download/1BGP.pdb')"
urllib.error.URLError: <urlopen error Tunnel connection failed: 403 Forbidden>
```

## The restriction bites only in the condition that keeps the ligands

On meeting a residue absent from its built-in templates, `mk_prepare_receptor`
tries to fetch that residue's CIF from rcsb.org at run time. This happens only for
inputs in which the HETATM records were kept.

| Condition | logs | fetch attempts | pdbqt | why there is no output |
|---|---|---|---|---|
| `meeko_nowat` (ligands present) | 100 | 77 | 21 | network 77 / other than network 2 (4MLK on alternate locations, 6ZHF on a residue, BEF, that could not be built) |
| `meekoA` | 100 | 0 | 99 | other than network 1 |
| `meekoA_noa` | 100 | 0 | 61 | other than network 39 |
| `meekoB` | 100 | 0 | 98 | other than network 2 |
| PDB2PQR, 4 conditions (`pqr_raw`, `pqr_hetall`, `pqr_nowat`, `pqr100r`) | 724 | 0 | — | — |

A fetch attempt was counted when the body of a log contains `rcsb.org`, `urlopen`
or `download CIF`. The string left in a log takes this form.

```
Could not download CIF file for HEM. If this is a certificate error and you trust
rcsb.org is not spoofed, you could try the --ignore_https_cert option.
Error: <urlopen error Tunnel connection failed: 403 Forbidden>
```

Each of the 77 logs carries, before the failed download, a warning that lists every
input residue without a template:

```
RuntimeWarning: Input residues {'A:502': 'NA', 'A:400': 'HEM'} not in residue_templates
```

The run then stops at the first download that fails, so the error itself names one
residue. Across the 77 logs the warnings name 97 distinct residue codes; the most
frequent are SO4 in 15 structures, NA in 6, GOL in 5 and MSE in 4. Which residue
the download is attempted for first is not fixed between runs: in the original run
the failed downloads were for 62 distinct codes, and in a later run of the same
inputs for 64. The pdbqt files written were the same in both.
`scripts_receptor/summarize_meeko_logs.py` reads all of this from the logs, and its
table is `data_receptor/meeko_runs.csv`. The logs themselves contain the paths of
the machine they were written on and are not deposited.

## Notes for the reader

The 21 pdbqt files in `meeko_nowat` are not a pass rate for a machine with network
access. Running `run_meeko.py` on a machine that reaches rcsb.org will give a
different count. If a reproduction attempt produces 403, that is the restricted
network access of the environment, not a failure of Meeko.

`mk_prepare_receptor.py --help` offers `--cache_templates` and `--add_templates`.
These are the route to holding fetched templates locally, and neither was used in
these measurements. `--cache_templates` exists from 0.7.0 onward and is off by
default.

## What was confirmed by re-running

| Stage | How it was checked | Result |
|---|---|---|
| The three conditions from `make_inputs.py` | all 300 files rebuilt and compared byte for byte | identical |
| PDB2PQR, 3 conditions | 6 files rebuilt and compared by md5 | identical |
| Meeko, 4 conditions | 6 files rebuilt and compared by md5 | identical |
| `measure_field.py` | `--selftest`, 6 analytic checks | identical |
| The whole of Stage 1 and Stage 2, 2026-11-10 | run again from a copy of the deposit, in a separate working directory, by a process given only the manuscript and the deposit | see below |

In the run of 2026-11-10 the 100 PDB files could not be fetched (the same 403 as
above) and were taken from a local copy, which matched `full100_sha256.txt` 100 of
100. The preparation then reproduced every prepared file: the 300 inputs, 291 pqr
files from the three derived inputs and 279 pdbqt files matched
`prepared_sha256.txt`, and the 96 pqr files from the raw inputs were byte-identical
to the original ones (their checksums were added to the list afterwards). The six
measurement tables that existed at that date were regenerated byte for byte. The
software was the same as in the table at the top of this file.

That run also showed what the deposit lacked, and the following were added after
it: the Stage 1 command for the raw inputs, the scripts `summarize_arms.py`,
`measure_deletions.py`, `summarize_meeko_logs.py`, `compare_meeko_ligand.py`,
`summarize_field.py`, `run_validator.py`, `check_ffout.py` and `stage2.sh`, the
tables they write, and the recorded standard output in `data_receptor/stdout/`.

Which directory was built from what is recorded in `data_receptor/provenance.csv`.
