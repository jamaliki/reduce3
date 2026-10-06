#!/bin/bash
# Run reference reduce2 on test files, with and without flips; record timings
cd "$(dirname "$0")/out"
for id in "$@"; do
  for mode in H FH; do
    extra=""; [ "$mode" = FH ] && extra="add_flip_movers=True"
    start=$(date +%s.%N)
    ../env/bin/mmtbx.reduce2 ../pdbs/$id.pdb overwrite=True $extra output.filename=${id}${mode}.pdb > ${id}${mode}.log 2>&1
    rc=$?
    end=$(date +%s.%N)
    echo "$id $mode rc=$rc wall=$(echo "$end - $start" | bc)" >> timings.txt
  done
done
