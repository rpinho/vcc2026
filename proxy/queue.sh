#!/bin/zsh
# Sequential proxy scoring (one scorer at a time: two concurrent scorers pushed the Mini into 13.6 GB swap).
# Phase 1: v4a,v4b,v5a,v5b on every destination task present; phase 2: legacy v3b,v1,v2b. Cached results are skipped.
cd "$(dirname "$0")"
PY=${PYTHON:-../.venv/bin/python}
while pgrep -f "run.py --dest" >/dev/null; do sleep 30; done  # never two scorers at once
for phase in v4a,v4b,v5a,v5b v3b,v1,v2b; do
  for d in HEK293T K562 H1 HCT116 CD4; do
    [[ -f tasks/$d.h5ad || $d == H1 ]] || continue
    echo "$(date -u +%FT%TZ) start $d $phase" >> logs/queue.log
    nice -n 19 $PY run.py --dest $d --recipes $phase >> logs/score_$d.log 2>&1
    echo "$(date -u +%FT%TZ) end $d $phase exit=$?" >> logs/queue.log
    $PY calibrate.py > logs/calibrate_latest.txt 2>&1
  done
done
echo "$(date -u +%FT%TZ) QUEUE DONE" >> logs/queue.log
