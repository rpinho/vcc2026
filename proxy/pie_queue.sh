#!/bin/zsh
# PIE 4-fold ensembles for the proxy destinations whose lines no PIE fold trained on (H1, HCT116, CD4), sequential on MPS.
cd "$(dirname "$0")"
PY=${PYTHON:-../.venv/bin/python}
while pgrep -f "pie_tools.py run" >/dev/null; do sleep 30; done
for d in H1 HCT116 CD4; do
  [[ -f ../model_v6/pie/out/ens4_$d.npz ]] && continue
  $PY - <<PYEOF
import sys; sys.path.insert(0, '.')
import run, anndata as ad
a, c, panel = run.load_task('$d')
c.write_h5ad('tasks/${d}_controls.h5ad', compression='gzip')
import pandas as pd; pd.DataFrame({'target_gene': panel}).to_csv('tasks/${d}_pie_panel.csv', index=False)
PYEOF
  (cd ../model_v6/pie && nice -n 19 .venv/bin/python pie_tools.py run --h5ad $PWD/../../proxy/tasks/${d}_controls.h5ad --targets $PWD/../../proxy/tasks/${d}_pie_panel.csv --label $d) > logs/pie_$d.log 2>&1
  echo "$(date -u +%FT%TZ) pie $d exit=$?" >> logs/queue.log
done
