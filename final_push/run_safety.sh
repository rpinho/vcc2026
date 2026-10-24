#!/bin/sh
# Local safety-recipe build and independent package readback. Never submits.
set -eu
cd "$(dirname "$0")/.."
NAME=${1:-example_safety}
PY=${VCC_PYTHON:-.venv/bin/python}
PERTS=data/final/bundle/pert_counts.csv
GENES=data/gene_names.csv
CTX=D=data/final/bundle/context_D.h5ad,E=data/final/bundle/context_E.h5ad,F=data/final/bundle/context_F.h5ad
BANKS=model_billmB/final/banks
RESP=model_billmB/final/response_banks
WEIGHTS=model_billm/checkpoints
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1
mkdir -p final_push/logs
nice -n 10 "$PY" -B model_billmB/panel_responses.py --banks "$BANKS" --out "$RESP"
nice -n 10 "$PY" -B model_billmB/donor_coverage.py --targets "$PERTS" --data "$RESP" --out-prefix model_billmB/final/family_coverage
nice -n 10 "$PY" -B model_billmB/panel_pipeline.py --name "$NAME" --targets "$PERTS" --contexts "$CTX" --genes "$GENES" --banks "$BANKS" --response-banks "$RESP" --weights "$WEIGHTS" --probes
nice -n 10 "$PY" -B model_billmB/panel_pipeline.py --name "$NAME" --targets "$PERTS" --contexts "$CTX" --genes "$GENES" --banks "$BANKS" --response-banks "$RESP" --weights "$WEIGHTS"
nice -n 10 "$PY" -B model_billmB/panel_prep.py --name "$NAME" --targets "$PERTS" --contexts D,E,F --genes "$GENES"
