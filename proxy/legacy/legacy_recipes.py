"""v1 / v2b / v3b LOCO ports: a copy of rsi/recipes.py (predict(), recipes v1/v2b/v3b only) with every path pointed
at rsi's archived source banks READ-ONLY and the module root at proxy/legacy.
Numerical kernels come from rsi/vendor/legacy_*.py; v3 uses its actual frozen submission source priors.
the destination family is excluded. This feedback calibration retains the original public source roster;
rsi's additional reserved-family/target restrictions are not part of the seven submitted recipes."""
import sys, json, functools, re
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
from scipy import sparse
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
RSI = HERE.parent.parent / 'rsi'; VCC = HERE.parent.parent
import legacy_pipeline as pipe
sys.modules['pipeline'] = pipe; pipe.ROOT = HERE
import legacy_v2 as v2
import legacy_v3 as v3
import legacy_profiles as profiles
profiles.ROOT = HERE
NT = 'non-targeting'; SEALED = set()
# Hold out the destination family only. Removing RPE1 or arbitrary H1 rows would
# change v3's common-response prior and would no longer reconstruct the official arm.
H1_SEALED = set()


def allowed(name, destination):
    return name.split('_')[0] not in SEALED and not (name == destination or name.startswith(destination + '_'))


@functools.lru_cache(maxsize=6)
def source(name, profile='base'):
    folder = 'cd4_sources' if profile == 'CD4' else 'sources'
    d = dict(np.load(RSI / f'artifacts/{folder}/{name}.npz')); d['name'] = name; d['family'] = 'K562' if name.startswith('K562') else name
    d['gene_ids'] = np.array([pipe.annotation_aliases().get(g, '') for g in d['genes']]); d['target_ids'] = np.array([pipe.annotation_aliases().get(t, '') for t in d['targets']]); return d


def fam_allowed(name, destination):
    fam = 'K562' if name.startswith('K562') else name
    return allowed(name, destination) and fam != destination


def sources(destination):
    out = []
    for name in ['H1', 'HepG2', 'Jurkat', 'K562', 'RPE1', 'K562_GWPS_bulk']:
        if not fam_allowed(name, destination):
            continue
        d = source(name, 'base').copy(); keep = ~np.isin(d['targets'], list(H1_SEALED)) if name == 'H1' else np.ones(len(d['targets']), bool)
        for k in ['targets', 'target_ids', 'lfc', 'fine', 'n_cells']:
            d[k] = d[k][keep]
        d['lfc'] = d['lfc'].astype(np.float32) * float(d['lfc_scale']); out.append(d)
    return out


def control_summary(ctrl):
    x = sparse.csr_matrix(ctrl.X, dtype=float); dep = np.asarray(x.sum(1)).ravel(); cm = np.asarray(x.mean(0)).ravel(); norm = x.multiply((1e4 / dep)[:, None]); avg = np.asarray(norm.mean(0)).ravel()
    var = np.maximum(np.asarray(norm.multiply(norm).mean(0)).ravel() - avg**2, 0); theta = np.clip(avg**2 / np.maximum(var - avg * np.mean(1e4 / dep), 1e-4), .1, 1000)
    genes = ctrl.var_names.to_numpy(str); return dict(genes=genes, gene_ids=np.array([pipe.annotation_aliases().get(g, '') for g in genes]), mean=cm, libs=dep, theta=theta)


def generic_sources(destination, targets):
    # These are fitted public-source parameters of the submitted v3, not a new
    # proxy-panel mean shift. Re-fitting from rsi's truncated source universe
    # would change its all-eligible fallback and erase the original failure mode.
    frozen=np.load(VCC/'model_v3/artifacts/source_shifts.npz')
    description=json.loads((VCC/'model_v3/artifacts/source_selection.json').read_text())
    out = {}
    for name in ['H1', 'HepG2', 'Jurkat', 'K562', 'RPE1', 'K562_GWPS_bulk']:
        if not fam_allowed(name, destination):
            continue
        out[name]={field:frozen[f'{name}_{field}'] for field in ('genes','base','effect')}
        out[name].update(description[name])
    frozen.close()
    return out


@functools.lru_cache(maxsize=2)
def robust_bulk(targets, genes, generic_genes, generic_effect):
    a = ad.read_h5ad(VCC / 'public_data/K562_GWPS_bulk.h5ad'); X = np.asarray(a.X, dtype=float); n = a.obs.num_cells_filtered.fillna(0).to_numpy(float)
    labels = np.array([NT if NT in s else re.match(r'^[^_]+_(.+)_[^_]+_([^_]+)$', s).group(1) for s in a.obs_names]); valid = n >= 30; cc = (labels == NT) & valid & a.obs.core_control.to_numpy(bool)
    rawgenes = a.var.gene_name.to_numpy(str); unique = np.array(sorted(set(rawgenes))); proj = sparse.csr_matrix((np.ones(len(rawgenes)), (np.arange(len(rawgenes)), pd.Index(unique).get_indexer(rawgenes))), shape=(len(rawgenes), len(unique)))
    X = np.asarray(X @ proj); rawgenes = unique; assert tuple(rawgenes) == generic_genes
    cpm = X / X.sum(1, keepdims=True) * 1e6; rawbase = np.average(X[cc], weights=n[cc], axis=0); base = rawbase / rawbase.sum() * 1e6
    sd = np.maximum(np.sqrt(np.mean((np.log2((cpm[cc] + 5) / (base + 5)) * np.sqrt(n[cc] / 170)[:, None])**2, axis=0)), 1e-3); res = np.zeros((len(targets), len(genes)), np.float32)
    ix = profiles.align(rawgenes, np.array(genes)); ok = ix >= 0; generic = np.array(generic_effect)
    for i, t in enumerate(targets):
        sel = (labels == t) & valid
        if not sel.any():
            continue
        mu = np.average(X[sel], weights=n[sel], axis=0); mu = mu / mu.sum() * 1e6; lfc = np.log2((mu + 5) / (base + 5)) - generic; z = lfc / (sd * np.sqrt(170 / n[sel].sum()))
        constructs = np.log2((cpm[sel] + 5) / (base + 5)) - generic; signok = np.mean(np.sign(constructs) == np.sign(lfc)[None, :], axis=0) >= .75
        good = (np.abs(z) >= 4) & (np.abs(lfc) >= .25) & (base >= 5) & signok
        if sel.sum() == 1:
            good &= (n[sel].sum() >= 100) & (np.abs(z) >= 5)
        mapped = np.zeros(len(genes), bool); mapped[ok] = good[ix[ok]]; mapped &= np.array(genes) != t; res[i, mapped] = .3 * np.clip(lfc[ix[mapped]], -1.5, 1.5)
    return res


def predict(ctrl, targets, destination, recipe, seed=42):
    s = control_summary(ctrl); s['panel'] = np.array(targets); rng = np.random.default_rng(seed); X = sparse.csr_matrix(ctrl.X); genes = s['genes']; n = 400
    if recipe in ['v1', 'v2b']:
        src = sources(destination)
        v2._FINE = {d['name']: dict(lfc=d['fine'], genes=d['genes'], targets=d['targets'], baseline=d['baseline'], n_cells=d['n_cells'], n_controls=d['n_controls']) for d in src}; v2._POST_STATE = None
        for t in targets:
            if recipe == 'v1':
                mu = pipe.predict_mean(s, t, src, 'transfer'); x = pipe.sample_counts(mu, s['libs'], s['theta'], n, rng)
            else:
                mu = v2.predict(s, t, src, 'postfine100'); x = v2.posterior_transport(mu, X, s['mean'], s['theta'], n, rng)
            yield t, sparse.csr_matrix(x)
    elif recipe == 'v3b':
        common = generic_sources(destination, targets); shift, _ = profiles.common_effect(common, genes, s['mean'], exclude=destination); res = np.zeros((len(targets), len(genes)))
        if 'K562_GWPS_bulk' in common:
            z = common['K562_GWPS_bulk']; res = robust_bulk(tuple(targets), tuple(genes), tuple(z['genes']), tuple(z['effect']))
        for i, t in enumerate(targets):
            yield t, v3.emit(X, s['mean'], v3.mean_vector(s['mean'], shift, genes, t, res[i]), 0 if seed == 42 else seed, i)
        # The sequential calibration runs modern arms next; release large legacy
        # banks/posterior state after their final emitted matrix has been consumed.
        source.cache_clear();robust_bulk.cache_clear()
        v2._FINE={};v2._POST_STATE=None
    else:
        raise ValueError(recipe)
