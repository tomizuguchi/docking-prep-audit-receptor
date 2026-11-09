# 受容体準備の監査 — スクリプトとデータ

分解能で層化した 100 構造のタンパク質結晶構造に対し、ドッキングの前の受容体準備
が何をするかを測った。スクリプト、測定結果、再現の手順をここに置く。

対応する論文: 【要記入：受理後に書誌情報と DOI】

## 何を測ったか

| 条件 | 比べたもの |
|---|---|
| 配位子の除去順序 | PDB2PQR に HETATM を入れた場合と入れない場合 |
| ツールの選択 | 同一入力に対する PDB2PQR と Meeko |
| `-a` フラグ | Meeko でこれを付けた場合と付けない場合 |
| altloc の設定 | `--default_altloc A` と `B` |

いずれも出力にエラーも警告も出ないまま、受容体の電荷が変わる。

## ライセンス

| 対象 | ライセンス |
|---|---|
| ソフトウェア（`scripts_receptor/` の 14 本） | MIT（`LICENSE`） |
| データと文章（`data_receptor/`, `README.md`, `environment.md`） | CC BY 4.0（`LICENSE-DATA`） |

`LICENSE` は MIT の条文そのもので、適用範囲はこの表が定める。`LICENSE-DATA` には
CC BY 4.0 の適用範囲と、第三者由来のものの出自を書いてある。

引用の形は `CITATION.cff` にある。

## 導入

解析に要るのは Python 3 と numpy だけ。

```
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

準備そのものを再現する場合は、これに加えて次が要る。いずれも pip で入る。

```
pip install "pdb2pqr==3.8.0" "propka==3.5.1" "meeko==0.8.0"
```

バージョンは厳密に合わせること。Meeko は 0.7 系と 0.8.0 で既定のテンプレートと
フラグ名が違う。詳細は `environment.md`。

## 再現の 2 段階

寄託物には準備済みの中間ファイル（pqr と pdbqt、合わせて 205 MB）を入れていない。
かわりに、その sha256 を `data_receptor/prepared_sha256.txt` に置いた。870 件。
自分で作り直したものがこちらのものと同じかを、これで照合できる。

### 第 1 段階 — 構造を取得して準備する

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

最後の照合で何も出なければ、こちらが測ったのと同じファイルが得られている。

`meeko_nowat` だけは通信を要する。`mk_prepare_receptor` が内蔵テンプレートに無い
残基の CIF を rcsb.org から取りに行くためで、こちらの測定は rcsb.org に到達でき
ない環境で行った。そのため 100 本中 21 本しか通っていない。通信のある環境では
別の本数になる。`environment.md` を見ること。

### 第 2 段階 — 準備済みのファイルから測定を作り直す

```
python3 scripts_receptor/profile_pdb.py prof full100 data_receptor/prof100.csv
python3 scripts_receptor/classify_titr.py pqr100r data_receptor/titr_none.csv
python3 scripts_receptor/compare_arms.py pqr_hetall pqr_nowat pqr100r data_receptor/arms_diff.csv
python3 scripts_receptor/compare_tools.py meekoA pqr100r data_receptor/tool_diff.csv
python3 scripts_receptor/diff_pdbqt.py meekoA meekoB data_receptor/altloc_diff.csv
python3 scripts_receptor/measure_field.py recept meekoA pqr100r --out data_receptor/field_diff.csv
```

`measure_field.py` は `--selftest` で解析解と照合できる。ファイルを要らない。

## スクリプト

```
fetch_strata.py   分解能で層化して 100 構造を選び RCSB から取得する（標本を選び直す）
fetch_manifest.py マニフェストの ID で同じ 100 構造を取得する（再現用。こちらを使う）
make_inputs.py    3 条件の入力を作る
profile_pdb.py    入力側の素性を数える
run_pqr.py        PDB2PQR を走らせる
run_meeko.py      Meeko を走らせる
classify_titr.py  pqr の滴定状態を残基名ではなく水素の有無から決める
compare_arms.py   HETATM の有無で滴定状態がどう変わるかを残基ごとに比べる
compare_tools.py  同じ入力から Meeko と PDB2PQR の結果を比べる
diff_pdbqt.py     altloc 設定だけを変えた 2 つの pdbqt を比べる
measure_proximity.py  状態が動いた残基の金属・配位子からの距離
measure_lys.py        中性リジン率
measure_field.py      2 つの電荷セットが受容体表面に作る静電ポテンシャルの差
check_receptor.py     入力の検査と出力の検算
```

標本を選び直す場合は `fetch_strata.py` を使うが、検索の母集団がその日の RCSB の
内容になるので、同じシードでも標本が変わりうる。測定に使ったのと同じ 100 構造を
取るには `fetch_manifest.py` を使うこと。

## データ

```
full100_manifest.csv   100 構造の層・分解能・抽出順・シード
full100_sha256.txt     取得した PDB 100 本の sha256
prepared_sha256.txt    準備済み中間ファイル 870 件の sha256
prof100.csv            入力側の素性
titr_none.csv          滴定状態
arms_diff.csv          配位子除去順で動いた残基
arms_diff_3way.csv     3 条件の比較
tool_diff.csv          Meeko と PDB2PQR で形式電荷が食い違う残基
altloc_diff.csv        altloc 設定の違い
field_diff.csv         接触面の静電ポテンシャルの差（96 構造）
provenance.csv         どのディレクトリが何から作られたか
```

## 読むときの注意

`provenance.csv` を先に読むこと。どのディレクトリがどの入力から作られたかと、
それがどう確かめられたかが書いてある。

`meeko_nowat` の 21 本は、通信のある環境での通過率ではない。走行環境の通信制限
で 77 本が落ちた結果である。`environment.md` の「走行環境」を見ること。

Meeko の版は 0.8.0。0.7 系とはフラグ名と既定テンプレートが違う。

## MANIFEST.sha256 の作り方

```
cd <寄託ツリー>
find . -type f ! -name MANIFEST.sha256 \
     ! -path './.git/*' ! -path '*/__pycache__/*' ! -path '*/.ruff_cache/*' \
     ! -name '*.pyc' ! -name '._*' ! -name '.DS_Store' \
  | sort | xargs shasum -a 256 > MANIFEST.sha256
```

照合は `shasum -a 256 -c MANIFEST.sha256`。

## linter

`ruff.toml` にこの寄託物の設定がある。`ruff check .` で何も出ない。除外している
規則には 1 件ずつ理由を書いた。見たうえで直さないと決めたものだけである。
