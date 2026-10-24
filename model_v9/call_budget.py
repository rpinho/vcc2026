"""Borrow v8A's observed per-group call budgets, without new amplitude fitting."""
import runtime
import numpy as np
R=runtime.ROOT
_rows=runtime.json.loads((R.parent/'model_v8/artifacts/candA_gates.json').read_text())['records']
assert len(_rows)==900
CALLS={(r['context'],r['target']):int(r['n_calls']) for r in _rows}
MEDIAN={c:int(round(float(np.median([r['n_calls'] for r in _rows if r['context']==c])))) for c in 'ABC'}
NCONF={'A':694,'B':815,'C':840}
def budget(context,target=None):
 return CALLS[context,target] if target is not None else MEDIAN[context]
def supported(lfc,positive,negative,minimum=2):
 direction=np.sign(lfc)
 agree=np.where(direction>0,positive,negative);oppose=np.where(direction>0,negative,positive)
 return np.isfinite(direction)&(direction!=0)&(agree>=minimum)&(oppose==0)
