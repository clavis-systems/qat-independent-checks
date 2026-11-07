#!/bin/bash
# usage: run.sh name cpus mem script args...
n=$1; c=$2; m=$3; shift 3
MSYS_NO_PATHCONV=1 docker run --rm --name $n --cpus $c --memory $m -v "<work>/ricontrollo_256:/w" -v "<work>/tracker/data/observable-estimations/circuit-models/operator_loschmidt_echo:/data:ro" -w /w/code ricontrollo256:py313 "$@"
