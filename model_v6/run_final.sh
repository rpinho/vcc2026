#!/bin/zsh
# Oct 22 runbook (also the mock-final dry run): panel banks -> build -> vcc prep -> decode-verify. Never submits.
# usage (ABSOLUTE paths): [USE_PIE=1] ./run_final.sh NAME PERTS_CSV "D=/path/context_D.h5ad,E=...,F=..." [ARM=v5b] [GENES_CSV]
#   PERTS_CSV: the controls bundle's pert_counts.csv (target_gene column, optionally per context)
# Writes out/NAME.vcc, artifacts/NAME_{build,prep}.json, logs/NAME_*.log, logs/NAME_timeline.txt
set -u
cd "$(dirname "$0")"
NAME=$1; PERTS=$2; CTX=$3; ARM=${4:-v5b}; GENES=${5:-../data/gene_names.csv}
PY=${PYTHON:-../.venv/bin/python}
BANKS=artifacts/banks/$NAME; mkdir -p $BANKS logs
TL=logs/${NAME}_timeline.txt; : > $TL
stamp() { echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) $(date +%s) $*" | tee -a $TL; }
stamp start banks
pids=()
for s in K562_GWPS_raw HCT116 HEK293T CD4 K562_GWPS_bulk,H1,promoters; do
  nice -n 19 $PY banks.py --targets $PERTS --out $BANKS --sources $s > logs/${NAME}_banks_${s//,/_}.log 2>&1 &
  pids+=$!
done
# optional PIE (USE_PIE=1): 4-fold ensemble per context on MPS, concurrent with the fetches (~16-20 min/context)
if [[ ${USE_PIE:-0} == 1 ]]; then
  for kv in ${(s:,:)CTX}; do
    ( cd pie && nice -n 19 .venv/bin/python pie_tools.py run --h5ad ${kv#*=} --targets $PERTS --label ${kv%%=*} ) > logs/${NAME}_pie_${kv%%=*}.log 2>&1
  done &
  piepid=$!
fi
fail=0; for p in $pids; do wait $p || fail=1; done
[[ -n ${piepid:-} ]] && { wait $piepid; stamp end pie; }
stamp end banks fail=$fail
# a failed remote source is not fatal: model.load_banks uses whatever banks exist (K562_GWPS_bulk replaces raw)
CONTEXTS=$(echo $CTX | tr ',' '\n' | cut -d= -f1 | paste -sd, -)
stamp start build $ARM
nice -n 19 $PY pipeline.py --arm $ARM --targets $PERTS --contexts $CTX --genes $GENES --banks $BANKS \
  --promoters $BANKS/promoter_pairs.csv --out $NAME $([[ ${USE_PIE:-0} == 1 ]] && echo --pie-dir pie/out --pie-prefix ens4) > logs/${NAME}_build.log 2>&1 || { stamp build FAILED; exit 1; }
stamp end build
stamp start prep
nice -n 19 $PY prep.py $NAME --perts $PERTS --contexts $CONTEXTS --genes $GENES --delete-input > logs/${NAME}_prep.log 2>&1
code=$?
stamp end prep exit=$code
exit $code
