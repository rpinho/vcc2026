"""On the official panel/context A with model_v4 banks and no exclusions, v6 heads + v5b gate must equal
model_v5's (byte-identical arrays). Runs model_v5 code read-only in a subprocess-free import."""
import runtime, sys, json, importlib.util
import numpy as np, pandas as pd, anndata as ad
sys.path.insert(0, str(runtime.ROOT))
import model as m6, v5 as v5_6
from pipeline import read_genes
genes = read_genes(); targets = pd.read_csv(runtime.VCC / 'data/pert_counts.csv').target_gene.astype(str).to_numpy()
a = ad.read_h5ad(runtime.VCC / 'data/context_A.h5ad')
# load model_v5/model.py under another name; it imports `runtime` and `vendor` relative to cwd, so alias them
import types
spec = importlib.util.spec_from_file_location('model_v5_model', runtime.VCC / 'model_v5/model.py')
m5 = importlib.util.module_from_spec(spec)
saved = sys.modules.get('vendor'); spec.loader.exec_module(m5)
m6.configure()
q6, d6, mean6, bulk6 = m6.template(a, genes, pool_k=4, seed=20260910)
q5, d5, mean5, bulk5 = m5.template(a, genes, pool_k=4, seed=20260910)
assert np.array_equal(q5, q6) and np.array_equal(d5, d6)
out = {}
(p6, b6), c6 = m6.heads(targets, genes, mean6, bulk6, 'anchor_single_mean')
(p5, b5), c5 = m5.heads(targets, genes, mean5, bulk5, 'anchor_single_mean')
out['v4a_heads_equal'] = bool(np.array_equal(p5, p6) and np.array_equal(b5, b6))
src6 = [m6.aligned_bank(n, targets, genes, 'log2fc') for n in m6.load_banks(None, True)] + [m6.cd4_family(targets, genes, 'log2fc')]
src5 = [m5.aligned_bank(n, targets, genes, 'log2fc') for n in m5.load_banks(None, True)] + [m5.cd4_family(targets, genes, 'log2fc')]
out['gate_sources_equal'] = all(np.array_equal(x, y) for s6, s5 in zip(src6, src5) for x, y in zip(s6, s5))
out['banks'] = m6.load_banks(None, True)
# fill must be a no-op on a fully covered panel
info = {}
(pf, bf), cf = m6.heads(targets, genes, mean6, bulk6, 'anchor_single_mean', fill=dict(embedding=None), info=info)
out['fill_noop_on_covered'] = bool(np.array_equal(pf, p6) and np.array_equal(bf, b6)); out['fill_info'] = info.get('coverage_fill')
print(json.dumps(out, indent=1)); runtime.save(runtime.ROOT / 'artifacts/verify_equivalence.json', out)
assert out['v4a_heads_equal'] and out['gate_sources_equal']
