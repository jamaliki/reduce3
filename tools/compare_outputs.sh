#!/bin/bash
# Run reduce3 on the reference PDB files and diff against reduce2's outputs.
# Usage: tools/compare_outputs.sh [--compat|--fixed] [names...]
# Needs the reference data described in reference/README.md (or REDUCE3_REF
# pointing at a directory with the same layout).
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
ref="${REDUCE3_REF:-$here/reference}"
mode="--compat"
if [ "${1:-}" = "--compat" ] || [ "${1:-}" = "--fixed" ]; then mode="$1"; shift; fi
names=("$@")
if [ ${#names[@]} -eq 0 ]; then
  names=(); for f in "$ref"/pdbs/*.pdb; do b=$(basename "$f" .pdb); names+=("$b"); done
fi
work="${TMPDIR:-/tmp}/reduce3_compare"; mkdir -p "$work"
bin="$here/target/release/reduce3"
chem="$ref/env/lib/python3.12/site-packages/chem_data"
# timing lines and the header differ by nature; deletions are listed in
# Python set order, so compare them sorted
filter() {
  local t; t=$(grep -v -e "Time to" -e "^reduce[23] v" -e "^ /" -e "^ --" "$1")
  printf '%s\n' "$t" | grep -v "^  Deleting "
  printf '%s\n' "$t" | grep "^  Deleting " | sort
}
for b in "${names[@]}"; do
  for kind in H FH; do
    flips=False; [ $kind = FH ] && flips=True
    out="$work/$b$kind.pdb"
    t0=$(python3 -c 'import time; print(time.time())')
    "$bin" $mode -q --chem-data "$chem" "$ref/pdbs/$b.pdb" add_flip_movers=$flips output.filename="$out" 2>"$work/$b$kind.err"
    rc=$?
    t1=$(python3 -c 'import time; print(time.time())')
    secs=$(python3 -c "print('%.3f' % ($t1-$t0))")
    if [ $rc -ne 0 ]; then echo "$b $kind: FAILED ($(head -1 "$work/$b$kind.err"))"; continue; fi
    pdbdiff=$(diff "$ref/out/$b$kind.pdb" "$out" | grep -c '^[<>]')
    txtdiff=$(diff <(filter "$ref/out/$b$kind.txt") <(filter "$work/$b$kind.txt") | grep -c '^[<>]')
    echo "$b $kind: pdb diff lines $pdbdiff, description diff lines $txtdiff, ${secs}s"
  done
done
