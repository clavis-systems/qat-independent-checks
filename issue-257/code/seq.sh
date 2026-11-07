#!/bin/sh
# usage: seq.sh <chi> <qasm1> <label1> <qasm2> <label2> ...  (one circuit after the other; restartable)
chi=$1; shift
while [ $# -ge 2 ]; do
  q=$1; lab=$2; shift 2
  out=/res/${lab}_chi${chi}
  if [ -f "$out.done" ]; then continue; fi
  echo "=== start $(date -u +%FT%TZ) chi=$chi" >> "$out.log"
  python /code/run_mps.py --qasm "$q" --chi "$chi" --out "$out.json" --ckpt "/ckpt/${lab}_chi${chi}.pkl" >> "$out.log" 2>&1 && touch "$out.done"
done
echo "seq finished $(date -u +%FT%TZ)" >> "/res/seq_chi${chi}.finished"
exec sleep infinity
