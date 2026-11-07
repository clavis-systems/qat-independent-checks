#!/bin/bash
B=<work>
MSYS_NO_PATHCONV=1 docker run --rm --name ${NAME:-r257} --cpus ${THR:-8} --memory 5g -e PYTHONPATH=/w/pylib:/w/code -e OMP_NUM_THREADS=${THR:-8} -e OPENBLAS_NUM_THREADS=${THR:-8} -v "$B/ricontrollo_257:/w" -v "$B/allegato_254/estratto:/data:ro" -v "$B/tracker/data/observable-estimations/circuit-models/su2_hadron_dynamics_lsh:/t:ro" python:3.13-slim "$@"
