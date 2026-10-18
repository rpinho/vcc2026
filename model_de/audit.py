"""Official prediction-side DE against public controls; no treated labels required.
This audits achieved calls/ranking, not accuracy. Builders can import prediction_de.
The official gate and per-target BH must be applied globally, not per gene chunk.
"""
def prediction_de(x, controls, target, *, threads=1):
    import anndata as ad
    import numpy as np,pandas as pd
    from scipy import sparse
    from dataclasses import replace
    from cell_eval2 import EvalConfig
    from cell_eval2.run import _compute_de_side
    if threads not in (1,2): raise ValueError('DE job cap is two threads')
    x=sparse.csr_matrix(x)
    if x.shape[1]!=controls.n_vars: raise ValueError('gene axes differ')
    if 'target_gene' in controls.obs and not (controls.obs.target_gene.astype(str)=='non-targeting').all():
        raise ValueError('controls must contain only non-targeting cells')
    cfg=EvalConfig.from_preset('vcc2026')
    cfg=replace(cfg,pert_col='target_gene',device='cpu',num_threads=threads,gather_threads=1,
                de=replace(cfg.de,backend='pdex'))
    a=ad.AnnData(sparse.vstack([x,sparse.csr_matrix(controls.X)],format='csr'),
        obs=pd.DataFrame({'target_gene':[str(target)]*x.shape[0]+['non-targeting']*controls.n_obs},
                         index=[f'de{i}' for i in range(x.shape[0]+controls.n_obs)]),var=controls.var.copy())
    frame=_compute_de_side(a,cfg=cfg,fp=None,store=None,side='pred')
    frame=frame.filter((frame['target']==str(target)) & (frame['feature']!=str(target)))
    return frame.sort(['p_adj','p_value','abs_log2_fold_change','feature'],descending=[False,False,True,False])
