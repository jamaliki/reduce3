# Reference data for validating Reduce3

This directory holds the scripts that produce Reduce2 reference data. The data itself is not
in git (the cctbx environment and the outputs come to several hundred MB). To recreate it:

```bash
cd reference
micromamba create -p env -c conda-forge python=3.12 cctbx-base=2026.9 chem_data=2026.9
mkdir -p pdbs out dumps
for id in 1a28 1crn 1d3z 1dfu 1ehz 1ubq 1xso 2oob 3gfh 3vyk 4fen 6oge 7c31; do
  curl -sSfo pdbs/$id.pdb https://files.rcsb.org/download/$id.pdb
done
./run_ref.sh 1a28 1crn 1d3z 1dfu 1ehz 1ubq 1xso 2oob 3gfh 3vyk 4fen 6oge 7c31
harness/run_all.sh
```

What each step produces:

* `run_ref.sh` runs `mmtbx.reduce2` on each structure, with and without flips. It writes
  `out/<id>H.pdb`, `out/<id>FH.pdb`, the description files, and wall-clock times in
  `out/timings.txt`.
* `harness/run_all.sh` runs `harness/dump_ref.py`. That script replays Reduce2 step by step
  and writes its intermediate state (atoms after placement, riding parameters, atom info,
  Movers, scores, deletions) to `dumps/<id>[_flips].json`. `harness/README.md` describes the
  format.

With the data in place:

```bash
tools/compare_outputs.sh --compat                   # end-to-end, byte for byte
cargo build --release --features refcheck
./target/release/reduce3 refcheck reference/dumps/*.json              # optimizer
./target/release/reduce3 hcheck reference/pdbs/1crn.pdb reference/dumps/1crn.json reference/env/lib/python3.12/site-packages/chem_data
./target/release/reduce3 wcheck reference/pdbs/1crn.pdb reference/dumps/1crn.json reference/env/lib/python3.12/site-packages/chem_data
```

`hcheck` compares hydrogen placement, and `wcheck` compares the optimizer inputs Reduce3
builds itself (bonded lists and atom info). `../docs/INTERPRETATION_SPEC.md` documents the
cctbx interpretation behavior that Reduce3 reproduces.
