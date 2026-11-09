# 走行環境

## バージョン

| ツール | バージョン | 出所 |
|---|---|---|
| PDB2PQR | 3.8.0 | `pdb2pqr --version` |
| propka | 3.5.1 | PDB2PQR が内部で呼ぶ |
| Meeko | 0.8.0 | `mk_prepare_receptor.py --help` と log 本文 |
| Python | 3.11 | |
| numpy | 2.x | 解析のみ |

Meeko の版はこの寄託物では 0.8.0 である。0.7 系との違いは次の 2 点。

- 既定のテンプレートの内容が違う。
- `-a` / `--allow_bad_res` が 0.8.0 では `-x` / `--delete_bad_res` に改名された。
  旧名も受け付けるが `--help` には出ない。`run_meeko.py` は旧名の `-a` を
  渡しており、0.8.0 でも通る。

0.7 系で走らせた場合、同じスクリプトでも出力が一致しない可能性がある。
`data_receptor/prepared_sha256.txt` との照合はバージョンを合わせてから行うこと。

## 2 つの計算機にまたがる

| 段 | 走らせた場所 | rcsb.org への到達 |
|---|---|---|
| 100 構造の取得 | 構造を集めた計算機 | あり |
| Meeko・PDB2PQR の全条件 | 通信制限のある計算機 | なし |

`full100/` の PDB は通信のある計算機で取得した（ディレクトリに macOS の付随
ファイル `._*` が混ざっているのはそのため）。その後の準備と測定は rcsb.org に
到達できない環境で走らせた。HTTP の接続はプロキシが 403 を返す。

```
$ python3 -c "import urllib.request; urllib.request.urlopen('https://files.rcsb.org/download/1BGP.pdb')"
urllib.error.URLError: <urlopen error Tunnel connection failed: 403 Forbidden>
```

## この制限がかかるのは配位子を残した条件だけ

`mk_prepare_receptor` は内蔵テンプレートに無い残基に出会うと、その残基の CIF を
実行時に rcsb.org から取得しようとする。HETATM を残した入力でだけこれが起きる。

| 条件 | log | 取得試行 | pdbqt | 出力なしの理由 |
|---|---|---|---|---|
| `meeko_nowat`（配位子あり） | 100 | 77 | 21 | 通信 77 / 通信以外 2 |
| `meekoA` | 100 | 0 | 99 | 通信以外 1 |
| `meekoA_noa` | 100 | 0 | 61 | 通信以外 39 |
| `meekoB` | 100 | 0 | 98 | 通信以外 2 |
| PDB2PQR 4 条件 | 724 | 0 | — | — |

取得試行の判定は、log 本文に `rcsb.org`、`urlopen`、`download CIF` のいずれかが
あるかで行った。log に残る文字列は次の形。

```
Could not download CIF file for HEM. If this is a certificate error and you trust
rcsb.org is not spoofed, you could try the --ignore_https_cert option.
Error: <urlopen error Tunnel connection failed: 403 Forbidden>
```

1 本の log に出てくる残基名は 1 つだけ（77/77）。照合できない最初の残基で止まり、
残りは報告されない。77 本で挙がった残基名は 62 種類で、構造数の多いものは
SO4 8 本、NA 3 本、HEM・GTP・MRD・GOL・ANP・ACT が各 2 本。結晶化添加物だけで
止まった構造は 22 本。

## 読むときの注意

`meeko_nowat` の pdbqt 21 本は、通信のある環境での通過率ではない。rcsb.org に
到達できる計算機で `run_meeko.py` を走らせれば別の本数になる。再現を試みて 403 が
出た場合、それは走行環境の通信制限であって Meeko の失敗ではない。

`mk_prepare_receptor.py --help` に `--cache_templates` と `--add_templates` がある。
取得したテンプレートをローカルに持たせる経路だが、この測定では使っていない。
`--cache_templates` は 0.7.0 以降にあり、既定では無効である。

## 再実行で確かめた範囲

| 段 | 確かめ方 | 結果 |
|---|---|---|
| `make_inputs.py` の 3 条件 | 全 300 ファイルを作り直して byte 比較 | 一致 |
| PDB2PQR 3 条件 | 6 本を作り直して md5 比較 | 一致 |
| Meeko 4 条件 | 6 本を作り直して md5 比較 | 一致 |
| `measure_field.py` | `--selftest`（解析解 6 件） | 一致 |

どのディレクトリが何から作られたかは `data_receptor/provenance.csv` にある。
