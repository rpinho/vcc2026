#!/bin/zsh
# One scorer at a time. Jobs: "DEST RECIPES". v5b_pie jobs wait for that destination's PIE ensemble file.
cd "$(dirname "$0")"
PY=${PYTHON:-../.venv/bin/python}; PIE=../model_v6/pie/out
while pgrep -f "run.py --dest" >/dev/null; do sleep 30; done
jobs=("K562 v4a,v4b,v5a,v5b" "HEK293T v5b_pie" "H1 v4a,v4b,v5a,v5b" "HCT116 v4a,v4b,v5a,v5b" "H1 v5b_pie" "CD4 v4a,v4b,v5a,v5b"
      "HCT116 v5b_pie" "CD4 v5b_pie" "HEK293T v3b,v1,v2b" "K562 v3b,v1,v2b" "H1 v3b,v1,v2b" "HCT116 v3b,v1,v2b" "CD4 v3b,v1,v2b")
for j in $jobs; do
  d=${j%% *}; r=${j#* }
  if [[ $r == v5b_pie ]]; then
    for i in {1..120}; do [[ -f $PIE/ens4_$d.npz ]] && break; sleep 60; done
    [[ -f $PIE/ens4_$d.npz ]] || { echo "$(date -u +%FT%TZ) skip $d v5b_pie (no PIE file)" >> logs/queue.log; continue; }
  fi
  echo "$(date -u +%FT%TZ) start $d $r" >> logs/queue.log
  nice -n 19 $PY run.py --dest $d --recipes $r >> logs/score_$d.log 2>&1
  echo "$(date -u +%FT%TZ) end $d $r exit=$?" >> logs/queue.log
  $PY calibrate.py > logs/calibrate_latest.txt 2>&1
done
echo "$(date -u +%FT%TZ) QUEUE DONE" >> logs/queue.log
