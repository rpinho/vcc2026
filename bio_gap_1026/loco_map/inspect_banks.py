import runtime as rt
rt.enforce(); rt.guard(True)
import numpy as np
manifest=rt.json.loads((rt.REPO/'model_billmB/data/manifest.json').read_text())
rows=[]; shared=None
for name,info in sorted(manifest.items()):
    if name=='H1_cached': continue
    with np.load(rt.REPO/f'model_billmB/data/{name}.npz') as z:
        genes=z['genes']; measured=z['measured']; targets=z['targets']; ctrl=z['control']
        gs=set(genes[measured]); shared=gs if shared is None else shared&gs
        row=dict(name=name,family=info['family'],targets=len(targets),genes=len(genes),measured=int(measured.sum()),
                 keys=z.files,control_shape=list(ctrl.shape),control_range=float(np.max(np.abs(ctrl-ctrl.mean(0)))))
        rows.append(row); print(row,flush=True)
print('ALL BANKS INTERSECTION',len(shared),flush=True)
rt.save('results/inspection.json',dict(banks=rows,shared_genes=len(shared),free_gib=rt.guard()))
