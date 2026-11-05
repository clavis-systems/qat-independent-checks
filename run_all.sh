#!/usr/bin/env bash
# Re-run all checks with one command.
#   bash run_all.sh PATH_TO_TRACKER_CLONE           fast checks: #240 and #248 (about 5 minutes)
#   bash run_all.sh PATH_TO_TRACKER_CLONE --dmrg    also the #239 DMRG (about 8 more minutes on 4 cores, about 8 GB)
# Results go to ./out/. Compare them with ./logs/.
set -euo pipefail
PY="${PYTHON:-python3}"
T="$(cd "$1" && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"
C="$HERE/code"
OUT="$HERE/out"
mkdir -p "$OUT"
FCI="$T/data/variational-problems/hamiltonians/2fe_2s/2fe_2s_30e_20o.fcidump"
H49="$T/data/variational-problems/hamiltonians/guided_sparse_ground_state_problem/49Q_v1/49Q_v1_hamiltonian.npy"

echo "== input fingerprints (must match README)"
sha256sum "$FCI" "$H49"
echo "bd6ef31773f7b082171f73faa6c98924948e2c250b6fbec74504d64361c8c6c7  expected for the FCIDUMP"
echo "c15ad8d74b5a41508a7cf99c389ad8a45ac50206a0eb7a992a72a3e11bf07566  expected for 49Q_v1_hamiltonian.npy"

cd "$OUT"   # run from a folder that holds no data files
echo "== #240 structure";        "$PY" -I "$C/struttura_240.py" "$H49"            | tee struttura_240.log
echo "== #240 symmetries";       "$PY" -I "$C/symmetries_240.py" "$H49"           | tee symmetries_240.log
echo "== #240 product state";    "$PY" -I "$C/prodotto_240.py" "$H49" prodotto_240.json | tee prodotto_240.log
echo "== #240 lower bound";      "$PY" -I "$C/lower_bound_240.py" "$H49" lower_bound_240.json 1.0 | tee lower_bound_240.log
echo "== #240 lower bound, kappa 1.01"; "$PY" -I "$C/lower_bound_240.py" "$H49" lower_bound_240_kappa1.01.json 1.01 | tee lower_bound_240_kappa1.01.log
echo "== lower bound self-test"; "$PY" -I "$C/lower_bound_240.py" --selftest      | tee lower_bound_240_selftest.log
echo "== #248 DOCI";             "$PY" -I "$C/doci_2fe2s.py" "$FCI" doci_248.json  | tee doci_248.log

if [ "${2:-}" = "--dmrg" ]; then
  echo "== #239 orbital localization"; "$PY" -I "$C/localizza_er.py" "$FCI" U_er.npy | tee localizza_er.log
  echo "== #239 DMRG (SU(2) singlet, M up to 750, 4 threads)"
  DMRG_ROT="$OUT/U_er.npy" "$PY" -I "$C/dmrg_2fe2s.py" "$FCI" SU2 fiedler 4 "$OUT/scratch_dmrg" dmrg_SU2_ER_M750.json 750 \
    | tee dmrg_SU2_ER_M750.log
  echo "== #239 extrapolation";        "$PY" -I "$C/estrapola.py" dmrg_SU2_ER_M750.json | tee estrapolazioni_M750.txt
fi
echo "== done; outputs in $OUT"
