#!/bin/sh
# usage: ladder.sh <qasm> <label> chi1 chi2 ...   (restartable: skips finished chi, resumes from checkpoint)
q=$1; lab=$2; shift 2
for chi in "$@"; do
  out=/res/${lab}_chi${chi}
  if [ -f "$out.done" ]; then continue; fi
  echo "=== start $(date -u +%FT%TZ) chi=$chi" >> "$out.log"
  python /code/run_mps.py --qasm "$q" --chi "$chi" --out "$out.json" --ckpt "/ckpt/${lab}_chi${chi}.pkl" >> "$out.log" 2>&1 && touch "$out.done"
done
echo "ladder finished $(date -u +%FT%TZ)" >> "/res/${lab}_ladder.finished"
exec sleep infinity
