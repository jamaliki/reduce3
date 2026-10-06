#!/bin/bash
# Run dump_ref.py on all reference structures, with and without flips, then compare the
# final atoms with the outputs of the real mmtbx.reduce2 in reference/out.
# Usage: reference/harness/run_all.sh [id ...]      (default: all reference/pdbs/*.pdb)
#        JOBS=4 reference/harness/run_all.sh        (parallel; timings in the dumps get noisier)
HERE="$(cd "$(dirname "$0")" && pwd)"
REF="$(dirname "$HERE")"
PY="$REF/env/bin/python"
LOGS="$REF/dumps/logs"
mkdir -p "$LOGS"
if [ $# -gt 0 ]; then ids="$*"; else ids=$(cd "$REF/pdbs" && ls *.pdb | sed 's/\.pdb$//'); fi
JOBS=${JOBS:-1}  # sequential by default so the recorded timings are not inflated

one() {  # $1 = id, $2 = optional --flips
  local id=$1 flag=$2 tag=$1
  [ -n "$flag" ] && tag=${id}_flips
  if "$PY" "$HERE/dump_ref.py" "$REF/pdbs/$id.pdb" $flag > "$LOGS/$tag.log" 2>&1; then
    echo "ok   $tag: $(tail -1 "$LOGS/$tag.log" | sed 's/.* in /in /')"
  else
    echo "FAIL $tag (see $LOGS/$tag.log): $(tail -1 "$LOGS/$tag.log")"
  fi
}

for id in $ids; do
  for flag in "" "--flips"; do
    while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 0.2; done
    one "$id" "$flag" &
  done
done
wait
"$PY" "$HERE/compare_out.py" "$@"
