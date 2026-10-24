"""billmA learned direction for a held-out family (LOCO checkpoints exist for H1 and CD4 only; the production
checkpoint saw every family and must not be used on a public fold). Read-only use of model_billm."""
import sys
from pathlib import Path
import numpy as np

BILLM = Path(__file__).resolve().parents[1] / 'model_billm'
_CACHE = {}


def billm_learned(fold, targets, control, genes, coefficient=None):
    """Returns the ensemble's learned response [T, G] on the canonical axis (NOT yet blended)."""
    assert fold in ('H1', 'CD4'), 'only the LOCO folds H1/CD4 are honest on public contexts'
    if str(BILLM) not in sys.path:
        sys.path.insert(0, str(BILLM))
    import torch
    torch.set_num_threads(2)
    from infer import Ensemble
    if fold not in _CACHE:
        _CACHE[fold] = Ensemble(fold)
    ens = _CACHE[fold]
    assert np.array_equal(ens.src.genes, genes), 'gene axis mismatch'
    assert fold not in {b['family'] for b in ens.src.banks}
    out = ens.predict(list(targets), np.asarray(control, np.float64))
    return out['learned'].astype(np.float32)
