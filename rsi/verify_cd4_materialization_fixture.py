"""Network-free corruption tests for the independent CD4 materialization audit."""
from datetime import datetime, timezone
from pathlib import Path
import tempfile

import runtime as rt

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from verify_cd4_materialization import audit_matrix, audit_selection

lanes = [f'lane{i:02}' for i in range(24)]
groups = [('x', 'perturbation', 4), ('y', 'perturbation', 4),
          ('non-targeting', 'input_control', 80), ('non-targeting', 'reference_control', 32)]
rows = [(target, role, lane, f'guide-{target}') for target, role, n in groups for lane in lanes for _ in range(n)]
selected = pd.DataFrame(rows, columns=['target_gene', 'evaluation_role', 'batch', 'guide_id'])
selected['source_row'] = np.arange(len(selected))
blueprint = dict(targets=['x', 'y'], n_real=96, audit=dict(target_batch_quotas=dict.fromkeys(lanes, 4),
    input_control_batch_quotas=dict.fromkeys(lanes, 80), control_batch_quotas=dict.fromkeys(lanes, 32)))
genes = np.array(['gene-b', 'gene-a', 'gene-c'])
counts = (np.arange(len(selected)*3).reshape(-1, 3) % 7 + 1).astype(np.int32)
a = ad.AnnData(sparse.csr_matrix(counts), obs=selected.copy(), var=pd.DataFrame(index=genes))
checks = []
assert audit_selection(selected, blueprint) == dict(input_control=1920, perturbation=192, reference_control=768)
with tempfile.TemporaryDirectory(dir=rt.ROOT/'artifacts/tmp') as temporary:
    path = Path(temporary)/'fixture.h5ad'
    a.write_h5ad(path)
    backed = ad.read_h5ad(path, backed='r')
    try:
        assert audit_matrix(backed, selected, genes) == int(counts.sum())
    finally:
        backed.file.close()
checks.append('Positive backed-HDF5 audit: 24 lanes, exact quotas, disjoint real rows, ordered genes and independent UMI sum')


def reject(label, fn):
    try:
        fn()
    except (AssertionError, ValueError, KeyError):
        checks.append(label)
    else:
        raise AssertionError('Corruption accepted: '+label)


bad = selected.copy()
bad.loc[1, 'source_row'] = 0
reject('Duplicate source cells rejected', lambda: audit_selection(bad, blueprint))
bad = selected.copy()
bad.loc[0, 'batch'] = lanes[1]
reject('Per-target lane imbalance rejected at unchanged cell count', lambda: audit_selection(bad, blueprint))
bad = selected.copy()
bad.loc[bad.evaluation_role == 'input_control', 'target_gene'] = 'x'
reject('Response-labelled input controls rejected', lambda: audit_selection(bad, blueprint))
reject('Incorrect control arm quota rejected', lambda: audit_selection(selected.iloc[:-1], blueprint))
reject('Frozen observation order change rejected', lambda: audit_matrix(a[::-1].copy(), selected, genes))
reject('Gene-order change rejected', lambda: audit_matrix(a[:, ::-1].copy(), selected, genes))
reject('Missing measured gene rejected', lambda: audit_matrix(a[:, :2].copy(), selected, genes))
bad_a = a.copy()
bad_a.X = bad_a.X.astype(np.float32)
bad_a.X.data[0] = 0.5
reject('Fractional count dtype rejected', lambda: audit_matrix(bad_a, selected, genes))
bad_a = a.copy()
bad_a.X.data[0] = -1
reject('Negative counts rejected', lambda: audit_matrix(bad_a, selected, genes))
bad_a = a.copy()
bad_a.X.data[:3] = 0
reject('Zero-depth cell rejected', lambda: audit_matrix(bad_a, selected, genes))
rt.save('artifacts/CD4_materialization_audit_fixture.json', dict(status='PASS',
    utc=datetime.now(timezone.utc).isoformat(), checks=checks, network='none',
    verifier_sha256=rt.sha('verify_cd4_materialization.py'), fixture_sha256=rt.sha(__file__),
    limitation='Synthetic corruption tests; live task remains UNVERIFIED until successful acquisition and full audit'))
print('PASS', len(checks), 'CD4 materialization audit checks')
