# An audit of receptor preparation — scripts and measurements

What receptor preparation does to a protein structure before docking, measured on
100 protein crystal structures sampled stratified by resolution. The scripts, the
measurements and the steps needed to reproduce them are here.

Companion paper: [to be filled in on acceptance: citation and DOI]

## What was measured

| Finding | What is compared |
|---|---|
| Ligand present at titration | PDB2PQR given the HETATM records, and PDB2PQR not given them |
| Choice of tool | PDB2PQR and Meeko on the same input |
| Meeko's `-a` flag | Meeko with it and without it, and the residues it deletes |
| Default altloc | `--default_altloc A` and `--default_altloc B` |

In the first three the receptor's charges change; in the fourth its coordinates
change and no charge does. In none of them does the run stop or return a non-zero
exit status, and in none does the prepared file record what was changed. Only the
`-a` flag produces a message: Meeko prints the numbers of the residues it drops to
standard error.

## Licences

| Covered | Licence |
|---|---|
| Software (everything in `scripts_receptor/`) | MIT (`LICENSE`) |
| Data and prose (`data_receptor/`, `README.md`, `environment.md`, and the metadata files) | CC BY 4.0 (`LICENSE-DATA`) |

`LICENSE` is the MIT text itself; the table above is what defines its scope.
`LICENSE-DATA` states the scope of CC BY 4.0 and the provenance of third-party
material.

The form to cite is in `CITATION.cff`.

## Installing

The analysis needs only Python 3 and numpy.

```
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Reproducing the preparation itself needs these as well. All three install with pip.

```
pip install "pdb2pqr==3.8.0" "propka==3.5.1" "meeko==0.8.0"
```

Match the versions exactly. Meeko 0.7.x and 0.8.0 differ in their built-in
templates and in the names of their flags. See `environment.md`.

## Reproducing, in two stages

The 100 PDB files (38 MB) and the prepared intermediate files (the pqr and pdbqt
files, 229 MB together) are not included. Their sha256 checksums are:
`data_receptor/full100_sha256.txt` for the 100 PDB files, and
`data_receptor/prepared_sha256.txt` for 966 prepared files (300 input PDB files
under the three derived conditions, 387 pqr files and 279 pdbqt files). Use those
to check whether the files you regenerate are the ones measured here.

The checks below print one line for every file that does not match and nothing
when all match. Where `shasum` is not installed, `sha256sum -c` takes the same
files.

### Stage 1 — fetch the structures and prepare them

```
python3 scripts_receptor/fetch_manifest.py data_receptor/full100_manifest.csv full100
shasum -a 256 -c data_receptor/full100_sha256.txt | grep -v ': OK$'

python3 scripts_receptor/make_inputs.py full100 .

python3 scripts_receptor/run_pqr.py full100   pqr_raw
python3 scripts_receptor/run_pqr.py in_hetall pqr_hetall
python3 scripts_receptor/run_pqr.py in_nowat  pqr_nowat
python3 scripts_receptor/run_pqr.py recept    pqr100r

python3 scripts_receptor/run_meeko.py recept   meekoA      -- --default_altloc A -a
python3 scripts_receptor/run_meeko.py recept   meekoB      -- --default_altloc B -a
python3 scripts_receptor/run_meeko.py recept   meekoA_noa  -- --default_altloc A
python3 scripts_receptor/run_meeko.py in_nowat meeko_nowat -- --default_altloc A -a

shasum -a 256 -c data_receptor/prepared_sha256.txt | grep -v ': OK$'
```

If the final check prints nothing, you have the same files that were measured here.
On two cores the eight preparation runs take about an hour in all.

`fetch_manifest.py` needs to reach files.rcsb.org. If an entry has been revised at
RCSB since the files were fetched its checksum will not match; the file measured here is then
not recoverable from this deposit, and everything downstream of it will differ.

Only `meeko_nowat` needs network access during preparation.
`mk_prepare_receptor.py` fetches the CIF of any residue absent from its built-in
templates from rcsb.org at run time, and the measurements here were made in an
environment that cannot reach rcsb.org. Only 21 of the 100 therefore went through.
On a machine with network access the count will be different, and the 21 pdbqt
checksums of that directory are the only ones that can be expected to match. See
`environment.md`.

### Stage 2 — rebuild the measurements from the prepared files

```
bash scripts_receptor/stage2.sh
```

This runs every analysis script, writes the tables to `regenerated/` and the
standard output of each command to `regenerated/stdout/`, and compares each file
with the deposited one in `data_receptor/`. Nothing in `data_receptor/` is
overwritten. The commands it runs are in the script, one per line. Several numbers
in the paper are printed by a script and not written to a table, which is why the
standard output is kept and compared as well.

`measure_field.py` takes about half an hour on two cores; `SKIP_FIELD=1 bash
scripts_receptor/stage2.sh` leaves it out. `check_ffout.py` runs PDB2PQR on four
structures.

`measure_field.py --selftest` checks the computation against analytic solutions and
needs no files.

## The scripts

```
fetch_strata.py   draws 100 structures stratified by resolution and fetches them from RCSB (redraws the sample)
fetch_manifest.py fetches the same 100 structures by the IDs in the manifest (for reproduction; use this one)
make_inputs.py    builds the three derived input conditions
profile_pdb.py    counts features of the inputs
run_pqr.py        runs PDB2PQR
run_meeko.py      runs Meeko
classify_titr.py  decides a pqr file's titration states from the hydrogens present, not from residue names
compare_arms.py   compares, per residue, how the titration states change with and without the HETATM records
summarize_arms.py net charge per structure across the inputs, and the split by metal and ligand
compare_tools.py  compares Meeko and PDB2PQR results from the same input
diff_pdbqt.py     compares two pdbqt files that differ only in the altloc setting
measure_proximity.py  distance from the metal or the ligand for residues whose state moved
measure_deletions.py  the residues Meeko's -a flag deleted, their types and charge
summarize_meeko_logs.py  what the Meeko logs say: deletions, failures, download attempts
compare_meeko_ligand.py  Meeko's output with the ligand kept against its output without it
measure_lys.py        fraction of neutral lysines
measure_field.py      difference in electrostatic potential that two charge sets produce on the receptor surface
summarize_field.py    the field table by difference in total charge, and the rank correlation
check_receptor.py     checks one input and cross-checks its outputs (the validator)
run_validator.py      runs check_receptor.py over the sample under each option set
check_ffout.py        checks that --ffout=AMBER names residues for their state
stage2.sh             runs all of the analysis and compares with the deposited tables
```

Use `fetch_strata.py` to redraw the sample, but note that the search population is
whatever RCSB holds that day, so the same seed can yield a different sample. To
obtain the same 100 structures that were measured, use `fetch_manifest.py`.
The original draw was `python3 scripts_receptor/fetch_strata.py fetch 20 full100`,
with the default seed (20261108). `fetch_strata.py` writes `manifest.csv`,
`query.json`, `excluded.csv` and `failed.csv` into the output directory; they are
deposited under the names `full100_manifest.csv` and so on, and the text inside
`full100_query.json` still refers to them by the undecorated names.

## The data

```
full100_manifest.csv   stratum, resolution, draw order and seed for the 100 structures
full100_query.json     the per-stratum search queries, the seed, and the exclusion rule
full100_excluded.csv   the 74 identifiers excluded for falling in more than one stratum
full100_failed.csv     the 4 identifiers that returned HTTP 404 and were replaced from the next ranks
full100_sha256.txt     sha256 of the 100 PDB files as fetched
prepared_sha256.txt    sha256 of the 966 prepared files (300 inputs, 387 pqr, 279 pdbqt)
provenance.csv         which directory was built from what
prof100.csv            features of the inputs
titr_none.csv          titration states read from the pqr files of the recept inputs
arms_diff_3way.csv     residues whose state moved, with the three derived inputs compared
netcharge_arms.csv     net charge of each structure under the three inputs; ligand and metal content
proximity_arms.csv     distance to the nearest ligand atom and metal for the residues in arms_diff_3way.csv
tool_diff.csv          residues on which Meeko and PDB2PQR assign different formal charges
proximity_tools.csv    the same distances for the residues in tool_diff.csv
deleted_residues.csv   the 306 residues Meeko's -a flag deleted
meeko_runs.csv         one row per structure and Meeko run: completion, deletions, failures, downloads
altloc_diff.csv        the difference the altloc setting makes; charges and coordinates are
                       reported separately, and identical_charges does not mean the files agree
field_diff.csv         difference in electrostatic potential on the contact surface, 96 structures
validator_runs.csv     exit status and rules raised by check_receptor.py under three option sets
stdout/                the standard output of each analysis command, as stage2.sh writes it
```

## Notes for the reader

Read `provenance.csv` first. It states which directory was built from which input,
and how that was checked.

The 21 files in `meeko_nowat` are not a pass rate for a machine with network
access. They are what remained after 77 of the 100 failed on the restricted
network access of the machine the measurements were run on. See the run
environment in `environment.md`.

The Meeko version here is 0.8.0. Its flag names and built-in templates differ from
the 0.7.x series.

The sign of a charge change caused by `-a` needs care. `deleted_residues.csv` and
`measure_deletions.py` give the charge of the residues that were deleted; the
receptor that remains moves by the opposite amount. More lysines and arginines are
deleted than aspartates and glutamates, so the receptor becomes more negative.

`prof100.csv` counts alternate locations over ATOM and HETATM records together
(46 structures). In ATOM records alone the count is 45; 1JZE has them only on its
ligands.

`profile_pdb.py` and the W5 rule of `check_receptor.py` count incomplete side
chains differently, by number of atoms and by atom names respectively, and give
288 and 299 over the sample. The figures in the paper are the former.

## How MANIFEST.sha256 was built

```
cd <deposit tree>
find . -type f ! -name MANIFEST.sha256 \
     ! -path './.git/*' ! -path '*/__pycache__/*' ! -path '*/.ruff_cache/*' \
     ! -name '*.pyc' ! -name '._*' ! -name '.DS_Store' \
  | sort | xargs shasum -a 256 > MANIFEST.sha256
```

Check it with `shasum -a 256 -c MANIFEST.sha256`.

## The linter

`ruff.toml` holds this deposit's configuration. `ruff check .` reports nothing
under it. Every ignored rule carries a written reason; each one is a finding that
was looked at and left as it stands.
