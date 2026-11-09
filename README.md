# An audit of receptor preparation — scripts and measurements

What receptor preparation does to a protein structure before docking, measured on
100 protein crystal structures sampled stratified by resolution. The scripts, the
measurements and the steps needed to reproduce them are here.

Companion paper: [to be filled in on acceptance: citation and DOI]

## What was measured

| Condition | What is compared |
|---|---|
| Order of ligand removal | PDB2PQR given the HETATM records, and PDB2PQR not given them |
| Choice of tool | PDB2PQR and Meeko on the same input |
| Meeko's `-a` flag | Meeko with it and without it |
| Default altloc | `--default_altloc A` and `--default_altloc B` |

In each case the receptor's charges change while the output carries no error and no
warning.

## Licences

| Covered | Licence |
|---|---|
| Software (the 14 scripts in `scripts_receptor/`) | MIT (`LICENSE`) |
| Data and prose (`data_receptor/`, `README.md`, `environment.md`) | CC BY 4.0 (`LICENSE-DATA`) |

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

The prepared intermediate files (the pqr and pdbqt files, 205 MB together) are not
included. Their sha256 checksums are, in
`data_receptor/prepared_sha256.txt`, 870 of them. Use those to check whether the
files you regenerate are the ones measured here.

### Stage 1 — fetch the structures and prepare them

```
python3 scripts_receptor/fetch_manifest.py data_receptor/full100_manifest.csv full100
shasum -a 256 -c data_receptor/full100_sha256.txt

python3 scripts_receptor/make_inputs.py full100 .

python3 scripts_receptor/run_pqr.py in_hetall pqr_hetall
python3 scripts_receptor/run_pqr.py in_nowat  pqr_nowat
python3 scripts_receptor/run_pqr.py recept    pqr100r

python3 scripts_receptor/run_meeko.py recept   meekoA      -- --default_altloc A -a
python3 scripts_receptor/run_meeko.py recept   meekoB      -- --default_altloc B -a
python3 scripts_receptor/run_meeko.py recept   meekoA_noa  -- --default_altloc A
python3 scripts_receptor/run_meeko.py in_nowat meeko_nowat -- --default_altloc A -a

shasum -a 256 -c data_receptor/prepared_sha256.txt
```

If the final check reports nothing, you have the same files that were measured here.

Only `meeko_nowat` needs network access. `mk_prepare_receptor` fetches the CIF of
any residue absent from its built-in templates from rcsb.org at run time, and the
measurements here were made in an environment that cannot reach rcsb.org. Only 21
of the 100 therefore went through. On a machine with network access the count will
be different. See `environment.md`.

### Stage 2 — rebuild the measurements from the prepared files

```
python3 scripts_receptor/profile_pdb.py prof full100 data_receptor/prof100.csv
python3 scripts_receptor/classify_titr.py pqr100r data_receptor/titr_none.csv
python3 scripts_receptor/compare_arms.py pqr_hetall pqr_nowat pqr100r data_receptor/arms_diff_3way.csv
python3 scripts_receptor/compare_tools.py meekoA pqr100r data_receptor/tool_diff.csv
python3 scripts_receptor/diff_pdbqt.py meekoA meekoB data_receptor/altloc_diff.csv
python3 scripts_receptor/measure_field.py recept meekoA pqr100r --out data_receptor/field_diff.csv
```

`measure_field.py --selftest` checks the computation against analytic solutions and
needs no files.

## The scripts

```
fetch_strata.py   draws 100 structures stratified by resolution and fetches them from RCSB (redraws the sample)
fetch_manifest.py fetches the same 100 structures by the IDs in the manifest (for reproduction; use this one)
make_inputs.py    builds the three input conditions
profile_pdb.py    counts features of the inputs
run_pqr.py        runs PDB2PQR
run_meeko.py      runs Meeko
classify_titr.py  decides a pqr file's titration states from the hydrogens present, not from residue names
compare_arms.py   compares, per residue, how the titration states change with and without the HETATM records
compare_tools.py  compares Meeko and PDB2PQR results from the same input
diff_pdbqt.py     compares two pdbqt files that differ only in the altloc setting
measure_proximity.py  distance from the metal or the ligand for residues whose state moved
measure_lys.py        fraction of neutral lysines
measure_field.py      difference in electrostatic potential that two charge sets produce on the receptor surface
check_receptor.py     checks the inputs and cross-checks the outputs
```

Use `fetch_strata.py` to redraw the sample, but note that the search population is
whatever RCSB holds that day, so the same seed can yield a different sample. To
obtain the same 100 structures that were measured, use `fetch_manifest.py`.

## The data

```
full100_manifest.csv   stratum, resolution, draw order and seed for the 100 structures
full100_query.json     the per-stratum search queries, the seed, and the exclusion rule
full100_excluded.csv   the 74 identifiers excluded for falling in more than one stratum
full100_failed.csv     the 4 identifiers that returned HTTP 404 and were replaced from the next ranks
full100_sha256.txt     sha256 of the 100 PDB files as fetched
prepared_sha256.txt    sha256 of the 870 prepared intermediate files
prof100.csv            features of the inputs
titr_none.csv          titration states
arms_diff_3way.csv     residues whose state moved, with the three heteroatom conditions compared
tool_diff.csv          residues on which Meeko and PDB2PQR assign different formal charges
altloc_diff.csv        the difference the altloc setting makes; charges and coordinates are
                       reported separately, and identical_charges does not mean the files agree
field_diff.csv         difference in electrostatic potential on the contact surface, 96 structures
provenance.csv         which directory was built from what
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
