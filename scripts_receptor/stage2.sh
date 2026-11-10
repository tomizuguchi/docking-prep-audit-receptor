#!/usr/bin/env bash
# Stage 2: rebuild every measurement table from the prepared files, and compare with the deposit.
#
#   bash scripts_receptor/stage2.sh [outdir]
#
# Run from the root of the deposit, after Stage 1 has produced the working directories
# (full100, in_hetall, in_nowat, recept, pqr_raw, pqr_hetall, pqr_nowat, pqr100r,
#  meekoA, meekoB, meekoA_noa, meeko_nowat). See README.md.
#
# Tables are written to <outdir> (default: regenerated/), never to data_receptor/, so that the
# deposited tables stay in place to be compared against. The standard output of every command is
# kept in <outdir>/stdout/, because several of the numbers reported in the paper are printed and
# not tabulated. At the end each regenerated file is compared with the deposited one. The recorded
# output names the output directory, so the comparison of stdout/ holds for the default <outdir> only.
#
# measure_field.py takes about half an hour on two cores. Set SKIP_FIELD=1 to leave it out; the
# summary of the field table is then made from the deposited field_diff.csv.
# check_ffout.py runs PDB2PQR on four structures and needs it installed.

set -u
S=scripts_receptor
D=data_receptor
R=${1:-regenerated}
mkdir -p "$R/stdout"

run() {
    name=$1
    shift
    if ! "$@" > "$R/stdout/$name.txt" 2>&1; then
        echo "!! $name exited with a non-zero status; see $R/stdout/$name.txt"
    fi
}

run profile_pdb             python3 $S/profile_pdb.py prof full100 "$R/prof100.csv"
run classify_titr           python3 $S/classify_titr.py pqr100r "$R/titr_none.csv"
run compare_arms            python3 $S/compare_arms.py pqr_hetall pqr_nowat pqr100r "$R/arms_diff_3way.csv"
run compare_tools           python3 $S/compare_tools.py meekoA pqr100r "$R/tool_diff.csv"
run diff_pdbqt              python3 $S/diff_pdbqt.py meekoA meekoB "$R/altloc_diff.csv"
run summarize_arms          python3 $S/summarize_arms.py pqr_hetall pqr_nowat pqr100r full100 "$R/arms_diff_3way.csv" "$R/netcharge_arms.csv"
run measure_proximity_arms  python3 $S/measure_proximity.py "$R/arms_diff_3way.csv" full100 --col ligand_removal_changed --out "$R/proximity_arms.csv"
run measure_proximity_tools python3 $S/measure_proximity.py "$R/tool_diff.csv" full100 --out "$R/proximity_tools.csv"
run measure_deletions       python3 $S/measure_deletions.py recept meekoA "$R/prof100.csv" --noa meekoA_noa --out "$R/deleted_residues.csv"
run summarize_meeko_logs    python3 $S/summarize_meeko_logs.py meekoA meekoB meekoA_noa meeko_nowat --prof "$R/prof100.csv" --out "$R/meeko_runs.csv"
run compare_meeko_ligand    python3 $S/compare_meeko_ligand.py meekoA meeko_nowat recept in_nowat
run measure_lys_raw         python3 $S/measure_lys.py pqr_raw $D/full100_manifest.csv
run measure_lys_recept      python3 $S/measure_lys.py pqr100r $D/full100_manifest.csv
run run_validator           python3 $S/run_validator.py full100 --out "$R/validator_runs.csv"
run check_ffout             python3 $S/check_ffout.py recept "$R/ffout" 1T5S 3DX2 1EA7 4NV6

if [ "${SKIP_FIELD:-0}" = 1 ]; then
    run summarize_field     python3 $S/summarize_field.py $D/field_diff.csv
else
    run measure_field       python3 $S/measure_field.py recept meekoA pqr100r --out "$R/field_diff.csv"
    run summarize_field     python3 $S/summarize_field.py "$R/field_diff.csv"
fi

echo
echo "Comparison with the deposited files"
same=0
differ=0
absent=0
for f in "$R"/*.csv "$R"/stdout/*.txt; do
    rel=${f#"$R"/}
    if [ ! -f "$D/$rel" ]; then
        echo "  not in the deposit: $rel"
        absent=$((absent + 1))
    elif cmp -s "$f" "$D/$rel"; then
        same=$((same + 1))
    else
        echo "  DIFFERS: $rel"
        differ=$((differ + 1))
    fi
done
echo "  identical $same, different $differ, not in the deposit $absent"
[ "$differ" -eq 0 ]
