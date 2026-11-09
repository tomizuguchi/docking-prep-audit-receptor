# Run environment

## Versions

| Tool | Version | Where this comes from |
|---|---|---|
| PDB2PQR | 3.8.0 | `pdb2pqr --version` |
| propka | 3.5.1 | called internally by PDB2PQR |
| Meeko | 0.8.0 | `mk_prepare_receptor.py --help` and the text of the logs |
| Python | 3.11 | |
| numpy | 2.x | analysis only |

The Meeko version in this deposit is 0.8.0. It differs from the 0.7.x series in two
ways.

- The contents of the built-in templates differ.
- `-a` / `--allow_bad_res` was renamed to `-x` / `--delete_bad_res` in 0.8.0. The
  old names are still accepted but do not appear in `--help`. `run_meeko.py`
  passes the old `-a`, which 0.8.0 accepts.

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
| `meeko_nowat` (ligands present) | 100 | 77 | 21 | network 77 / other than network 2 |
| `meekoA` | 100 | 0 | 99 | other than network 1 |
| `meekoA_noa` | 100 | 0 | 61 | other than network 39 |
| `meekoB` | 100 | 0 | 98 | other than network 2 |
| PDB2PQR, 4 conditions | 724 | 0 | — | — |

A fetch attempt was counted when the body of a log contains `rcsb.org`, `urlopen`
or `download CIF`. The string left in a log takes this form.

```
Could not download CIF file for HEM. If this is a certificate error and you trust
rcsb.org is not spoofed, you could try the --ignore_https_cert option.
Error: <urlopen error Tunnel connection failed: 403 Forbidden>
```

Each log names exactly one residue (77 of 77). The run stops at the first residue it
cannot match and the rest are never reported. The 77 logs name 62 distinct
residues; the most frequent are SO4 in 8 structures, NA in 3, and HEM, GTP, MRD,
GOL, ANP and ACT in 2 each. For 22 structures the run stopped on a
crystallization additive alone.

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

Which directory was built from what is recorded in `data_receptor/provenance.csv`.
