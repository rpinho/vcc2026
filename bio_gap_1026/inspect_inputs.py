import runtime as rt
rt.enforce();rt.guard(True)
import numpy as np,h5py
for name in ['manifold_1026/inputs/data/H1_full.npz','gatecall_1026/inputs/prediction/inference_H1.npz','model_billmB/data/K562.npz']:
    with np.load(rt.REPO/name) as z:
        print(name,[(k,z[k].shape,str(z[k].dtype)) for k in z.files],flush=True)
for p in (rt.REPO/'model_sidechain').rglob('*.h5ad'):
    print('h5ad',p)
for p in (rt.REPO/'model_sidechain/approx/sources/H1').iterdir(): print('H1 file',p.name,p.stat().st_size)
for p in (rt.REPO/'model_billmB/checkpoints/H1').glob('*.json'):
    d=rt.json.loads(p.read_text());print('checkpoint',p.name,{k:v for k,v in d.items() if k not in ['log','data_sha256']})
for p in (rt.REPO/'data/expansion/references').glob('*metadata.json'):
    print('metadata',p.name,p.read_text()[:1000])
