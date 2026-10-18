"""Panel source banks for ANY target list (the Oct 22 blocker: v4 banks hold only 437 targets).

`python banks.py --targets panel.csv --out artifacts/banks/<panel> [--sources ...]`

For every source the panel bank = the shared model_v4 bank rows (all 437, kept so source centring is
unchanged for known targets) UNION freshly built rows for panel targets the shared bank lacks:
  K562_GWPS_raw   anonymous byte ranges of Replogle GWPS raw dense HDF5 (figshare), 200 cells/target,
                  50 NTC per gem group (model_v4/prepare_k562_raw.py, parameterised)
  K562_GWPS_bulk  local public_data/K562_GWPS_bulk.h5ad (genome-wide pseudobulk; fallback for raw)
  H1              local public_data/H1.h5ad (model_v4/prepare_local.py)
  HCT116/HEK293T  X-Atlas-Orion Lance on HF, anonymous ranges (model_v4/probe.py + prepare_orion.py)
  CD4_*           Zhu/Marson GWCD4i pseudobulk + publisher DE stats, anonymous S3 ranges (prepare_cd4.py)
  promoters       GENCODE v47 TSS pairs + public CD4 guide metadata (prepare_promoters.py)
Every fetched cell is integrity-checked exactly as in model_v4. Network: public reads only.
"""
import runtime, argparse, time, threading, re, json, gzip, sys
import numpy as np, pandas as pd, anndata as ad, h5py, requests
from pathlib import Path
from scipy import sparse
from concurrent.futures import ThreadPoolExecutor

rt = runtime.ROOT
SHARED = runtime.VCC / 'public_data/model_v4_signature_banks'
FIGSHARE_GWPS = 'https://ndownloader.figshare.com/files/35775507'
CD4_BASE = 'https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/'
LOG = print


def log(*a):
    LOG(*a, flush=True)


def savez(path, **arrays):
    """Atomic: parallel per-source processes may read the directory while another writes."""
    path = Path(path); tmp = path.with_name(path.stem + '.partial.npz')
    np.savez_compressed(tmp, **arrays); tmp.replace(path)


# ----------------------------------------------------------------------------------------------
# moment accumulation (identical formulas to model_v4 prepare_k562_raw / prepare_orion)
# ----------------------------------------------------------------------------------------------
class Acc:
    def __init__(self, targets, batches, ng):
        self.targets = np.asarray(targets); self.batches = list(batches); nt = len(targets); nb = len(batches)
        self.nt, self.nb = nt, nb
        self.sums = np.zeros((nt + nb, ng)); self.qp = np.zeros_like(self.sums); self.q2 = np.zeros_like(self.sums)
        self.deps = np.zeros(nt + nb); self.n = np.zeros(nt + nb); self.bn = np.zeros((nt, nb)); self.bl = np.zeros((nt, nb))

    def add(self, ti, bi, cols, vals, depth):
        """One cell: ti target index or -1 for control; bi batch index; sparse (cols, vals) on bank genes;
        depth = mapped library size used for CPM (v4: counts on bank axis)."""
        g = ti if ti >= 0 else self.nt + bi
        q = vals / depth
        np.add.at(self.sums[g], cols, vals); np.add.at(self.qp[g], cols, q); np.add.at(self.q2[g], cols, q * q)
        self.deps[g] += depth; self.n[g] += 1
        if ti >= 0:
            self.bn[ti, bi] += 1; self.bl[ti, bi] += depth

    def add_block(self, groups, tis, bis, X, depth):
        """Vectorised add for a csr block X (cells x genes)."""
        X = sparse.csr_matrix(X)
        G = len(self.deps)
        M = sparse.csr_matrix((np.ones(len(groups)), (groups, np.arange(len(groups)))), shape=(G, len(groups)))
        P = sparse.diags(1.0 / depth) @ X
        self.sums += (M @ X).toarray(); self.qp += (M @ P).toarray(); self.q2 += (M @ P.multiply(P)).toarray()
        np.add.at(self.deps, groups, depth); np.add.at(self.n, groups, 1)
        ok = tis >= 0
        np.add.at(self.bn, (tis[ok], bis[ok]), 1); np.add.at(self.bl, (tis[ok], bis[ok]), depth[ok])

    def dump(self, path, genes):
        np.savez(path, targets=self.targets, batches=np.array(self.batches), genes=np.asarray(genes), sums=self.sums, qp=self.qp, q2=self.q2,
                 deps=self.deps, n=self.n, bn=self.bn, bl=self.bl)

    @classmethod
    def load(cls, path):
        d = np.load(path); a = cls(d['targets'], list(d['batches']), d['sums'].shape[1])
        for k in ['sums', 'qp', 'q2', 'deps', 'n', 'bn', 'bl']:
            setattr(a, k, d[k])
        return a, d['genes']

    def bank(self, name, genes, measured=None):
        nt = self.nt; cr = np.arange(nt, nt + self.nb); cr = cr[self.n[cr] > 0]
        assert len(cr) > 0, 'no controls'
        gp = self.sums[cr].sum(0) / self.deps[cr].sum(); gm = self.qp[cr].sum(0) / self.n[cr].sum() * 1e6
        tc = self.sums[:nt]; tn = self.n[:nt].astype(np.int64)
        tm = np.divide(self.qp[:nt] * 1e6, tn[:, None], out=np.zeros_like(tc), where=tn[:, None] > 0)
        nb = self.nb; allc = np.arange(nt, nt + nb)
        cn = self.n[allc]
        gvar = np.maximum(self.q2[cr].sum(0) / self.n[cr].sum() - (gm / 1e6)**2, 0) / max(self.n[cr].sum() - 1, 1)
        cp = np.divide(self.sums[allc], self.deps[allc, None], out=np.tile(gp, (nb, 1)), where=self.deps[allc, None] > 0)
        cm = np.divide(self.qp[allc], cn[:, None], out=np.tile(gm / 1e6, (nb, 1)), where=cn[:, None] > 0)
        cv = np.maximum(np.divide(self.q2[allc], cn[:, None], out=np.zeros_like(cm), where=cn[:, None] > 0) - cm * cm, 0) / np.maximum(cn[:, None] - 1, 1)
        cv[cn == 0] = gvar
        cw = np.divide(self.bn, tn[:, None], out=np.zeros_like(self.bn), where=tn[:, None] > 0)
        lw = np.divide(self.bl, self.bl.sum(1, keepdims=True), out=np.zeros_like(self.bl), where=self.bl.sum(1, keepdims=True) > 0)
        mp = lw @ cp; mm = cw @ cm * 1e6
        pv = np.maximum(np.divide(self.q2[:nt], tn[:, None], out=np.zeros_like(tc), where=tn[:, None] > 0) - (tm / 1e6)**2, 0) / np.maximum(tn[:, None] - 1, 1) * 1e12
        vc = (cw * cw) @ cv * 1e12
        return dict(source=np.array(name), targets=self.targets.astype(str), genes=np.asarray(genes).astype(str), n_cells=tn,
                    target_count_sums=tc.astype('f4'), target_mean_cpm=tm.astype('f4'), matched_control_probability=mp.astype('f4'),
                    matched_control_mean_cpm=mm.astype('f4'), global_control_probability=gp.astype('f4'), global_control_mean_cpm=gm.astype('f4'),
                    measured_genes=np.ones(len(genes), bool) if measured is None else measured, delta_cpm_variance=(pv + vc).astype('f4'),
                    target_mean_cpm_variance=pv.astype('f4'), matched_control_mean_cpm_variance=vc.astype('f4'), arithmetic_cpm_observed=np.array(True))


def union_bank(name, new, out_dir, shared=SHARED):
    """Shared v4 rows (all of them) + new rows for targets the shared bank lacks; same gene axis required."""
    path = shared / f'{name}.npz'; out = Path(out_dir) / f'{name}.npz'
    if not path.exists():
        if new is not None:
            savez(out, **new)
        return out
    with np.load(path) as d:
        old = {k: d[k] for k in d.files}
    if new is None or len(new['targets']) == 0:
        savez(out, **old); return out
    assert np.array_equal(old['genes'], new['genes']), f'{name}: gene axis differs between shared and new rows'
    keep = ~np.isin(new['targets'], old['targets'])
    merged = dict(old)
    for k, v in old.items():
        if isinstance(v, np.ndarray) and v.ndim >= 1 and v.shape[0] == len(old['targets']) and k not in ('genes', 'measured_genes', 'global_control_probability', 'global_control_mean_cpm'):
            merged[k] = np.concatenate([v, new[k][keep].astype(v.dtype)]) if k in new else v
    savez(out, **merged)
    log(name, 'union bank', len(old['targets']), '+', int(keep.sum()), '->', out)
    return out


def shared_targets(name):
    p = SHARED / f'{name}.npz'
    if not p.exists():
        return np.array([], str)
    with np.load(p) as d:
        return d['targets'].astype(str)


# ----------------------------------------------------------------------------------------------
# K562 GWPS raw (figshare dense HDF5)
# ----------------------------------------------------------------------------------------------
class Figshare:
    def __init__(self, url=FIGSHARE_GWPS):
        self.url = url; self.lock = threading.Lock(); self.link = {'until': 0., 'url': None, 'etag': None}
        self.s = requests.Session(); self.s.trust_env = False; self.bytes = 0

    def public_link(self):
        with self.lock:
            if time.monotonic() >= self.link['until']:
                with self.s.get(self.url + '?fresh=' + str(time.time_ns()), headers={'Range': 'bytes=0-0'}, stream=True, timeout=60) as r:
                    assert r.status_code == 206; self.link['url'] = r.url; self.link['etag'] = r.headers.get('ETag'); assert len(r.raw.read(2)) == 1
                self.link['until'] = time.monotonic() + 3
            return self.link['url'], self.link['etag']


def gwps_meta():
    obs = pd.read_parquet(rt / 'artifacts/K562_GWPS_raw.obs.parquet'); var = pd.read_parquet(rt / 'artifacts/K562_GWPS_raw.var.parquet')
    return obs, var


def gwps_layout():
    from remote_h5 import Remote
    r = Remote(FIGSHARE_GWPS, cap=50_000_000, block=65536)
    with h5py.File(r, 'r') as h:
        ds = h['X']; assert ds.chunks is None and ds.compression is None and ds.dtype == np.float32
        return ds.id.get_offset(), ds.shape


def gwps_fetch(rows, ncol, offset, cap_bytes=6_000_000_000, workers=12):
    """Yield (k, dense float64 row) for sorted source rows, in bounded 2-MB coalesced ranges."""
    fs = Figshare(); rowbytes = ncol * 4; plans = []
    for k, row in enumerate(rows):
        a = int(offset + row * rowbytes); b = a + rowbytes
        if plans and a <= plans[-1][1] + rowbytes and b - plans[-1][0] <= 2 * 1024**2:
            plans[-1][1] = b; plans[-1][2].append((k, a))
        else:
            plans.append([a, b, [(k, a)]])
    amount = sum(b - a for a, b, _ in plans); assert amount < cap_bytes, amount
    log('GWPS ranges', len(plans), 'planned GB', round(amount / 1e9, 3))
    local = threading.local()

    def fetch(plan):
        a, b, members = plan
        if not hasattr(local, 's'):
            local.s = requests.Session(); local.s.trust_env = False
        for attempt in range(6):
            try:
                url, etag = fs.public_link()
                with local.s.get(url, headers={'Range': f'bytes={a}-{b - 1}', 'If-Match': etag, 'Accept-Encoding': 'identity'}, stream=True, timeout=(10, 90)) as resp:
                    assert resp.status_code == 206, resp.status_code; buf = resp.raw.read(b - a + 1); assert len(buf) == b - a
                with fs.lock:
                    fs.bytes += b - a
                return [(k, np.frombuffer(buf, dtype=np.float32, count=ncol, offset=off - a)) for k, off in members]
            except Exception as err:
                if attempt == 5:
                    raise
                with fs.lock:
                    fs.link['until'] = 0.
                time.sleep(1 + attempt)
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for j, res in enumerate(ex.map(fetch, plans)):
            if j % 500 == 0:
                runtime.guard()
                log('GWPS', j, '/', len(plans), 'GB', round(fs.bytes / 1e9, 3), 's', round(time.time() - t0))
            for k, x in res:
                yield k, x
    gwps_fetch.last_bytes = fs.bytes


def gwps_raw_bank(targets, per_target=200, ntc_per_batch=50, seed=20261016):
    obs, var = gwps_meta()
    labels = obs.gene.astype(str).to_numpy(); batches = obs.gem_group.astype(str).to_numpy()
    rng = np.random.default_rng(seed); sel = []
    targets = [t for t in targets if (labels == t).any()]
    for t in targets:
        r = np.flatnonzero(labels == t); sel.extend(rng.choice(r, min(per_target, len(r)), replace=False))
    for b in sorted(set(batches)):
        r = np.flatnonzero((labels == 'non-targeting') & (batches == b))
        if len(r):
            sel.extend(rng.choice(r, min(ntc_per_batch, len(r)), replace=False))
    rows = np.array(sorted(sel), np.int64)
    genes0 = var.gene_name.to_numpy(str); genes = np.array(list(dict.fromkeys(genes0))); gmap = pd.Index(genes).get_indexer(genes0)
    bnames = sorted(set(batches[rows])); bix = pd.Index(bnames).get_indexer(batches[rows]); tix = pd.Index(targets).get_indexer(labels[rows])
    acc = Acc(targets, bnames, len(genes))
    proj = sparse.csr_matrix((np.ones(len(genes0)), (np.arange(len(genes0)), gmap)), shape=(len(genes0), len(genes)))
    offset, shape = gwps_layout(); assert shape == (len(obs), len(genes0))
    t0 = time.time(); buf = []; idx = []

    def flush():
        if not buf:
            return
        X = np.stack(buf); k = np.array(idx)
        Xg = (sparse.csr_matrix(X.astype(np.float64)) @ proj).toarray()
        assert np.all(Xg == np.floor(Xg)) and np.all(Xg >= 0)
        d = Xg.sum(1); assert np.all(d > 0)
        groups = np.where(tix[k] >= 0, tix[k], len(targets) + bix[k])
        acc.add_block(groups, tix[k], bix[k], sparse.csr_matrix(Xg), d)
        buf.clear(); idx.clear()
    for k, x in gwps_fetch(rows, len(genes0), offset):
        buf.append(x); idx.append(k)
        if len(buf) >= 512:
            flush()
    flush()
    acc.dump(rt / 'artifacts/tmp/acc_K562_GWPS_raw.npz', genes)
    b = acc.bank('K562_GWPS_raw', genes)
    log('GWPS raw bank', len(targets), 'targets', len(rows), 'cells', round(time.time() - t0), 's', round(getattr(gwps_fetch, 'last_bytes', 0) / 1e9, 2), 'GB')
    return b, dict(cells=len(rows), seconds=time.time() - t0, bytes=getattr(gwps_fetch, 'last_bytes', 0), targets=len(targets))


# ----------------------------------------------------------------------------------------------
# local h5ad sources (H1, K562_GWPS_bulk, Replogle lines): model_v4/prepare_local.py, parameterised
# ----------------------------------------------------------------------------------------------
def local_bank(name, targets):
    a = ad.read_h5ad(runtime.VCC / 'public_data' / f'{name}.h5ad', backed='r'); start = time.time()
    tset = set(targets)
    genes = a.var.gene_name.astype(str).to_numpy() if 'gene_name' in a.var else a.var_names.to_numpy(str)
    bulk = name.endswith('bulk')
    if bulk:
        lab = np.array(['non-targeting' if 'non-targeting' in s else re.match(r'^[^_]+_(.+)_[^_]+_([^_]+)$', s).group(1) for s in a.obs_names])
        n0 = a.obs.num_cells_filtered.fillna(0).to_numpy(float); valid = n0 >= 20; valid &= (lab != 'non-targeting') | a.obs.core_control.to_numpy(bool)
        batches = np.repeat('unmatched', a.n_obs)
    else:
        col = 'gene' if 'gene' in a.obs else 'target_gene'; lab = a.obs[col].astype(str).to_numpy()
        lab = np.array(['non-targeting' if s.lower() in ['non-targeting', 'non-targeting control', 'ntc', 'control', 'ctrl', 'unperturbed'] or s.startswith('non-targeting') else s for s in lab])
        n0 = np.ones(a.n_obs); valid = np.ones(a.n_obs, bool); bcol = next((c for c in ['gem_group', 'batch', 'sample'] if c in a.obs), None)
        batches = a.obs[bcol].astype(str).to_numpy() if bcol else np.repeat('global', a.n_obs)
    valid &= np.isin(lab, list(tset | {'non-targeting'})); rows = np.flatnonzero(valid)
    keys = list(dict.fromkeys(zip(batches[rows], lab[rows]))); lookup = {k: i for i, k in enumerate(keys)}
    group = np.array([lookup[k] for k in zip(batches[rows], lab[rows])]); n = np.bincount(group, weights=n0[rows], minlength=len(keys))
    shape = (len(keys), len(genes)); sums = np.zeros(shape); qp = np.zeros(shape); q2 = np.zeros(shape); deps = np.zeros(len(keys))
    for lo in range(0, len(rows), 512):
        rr = rows[lo:lo + 512]; x = sparse.csr_matrix(a.X[rr], dtype=float); gg = group[lo:lo + 512]
        if bulk:
            x = x.multiply(n0[rr, None]).tocsr()
        dep = np.asarray(x.sum(1)).ravel(); assert np.all(dep > 0)
        member = sparse.csr_matrix((np.ones(len(rr)), (gg, np.arange(len(rr)))), shape=(len(keys), len(rr)))
        sums += (member @ x).toarray(); proportion = x.multiply((1 / dep)[:, None]).tocsr()
        if bulk:
            proportion = proportion.multiply(n0[rr, None]).tocsr()
        qp += (member @ proportion).toarray()
        if not bulk:
            q2 += (member @ proportion.multiply(proportion)).toarray()
        np.add.at(deps, gg, dep)
    if len(set(genes)) != len(genes):
        unique = np.array(list(dict.fromkeys(genes)))
        proj = sparse.csr_matrix((np.ones(len(genes)), (np.arange(len(genes)), pd.Index(unique).get_indexer(genes))), shape=(len(genes), len(unique)))
        sums = sums @ proj; qp = qp @ proj; q2 = q2 @ proj; genes = unique
    controls = {b: i for i, (b, t) in enumerate(keys) if t == 'non-targeting'}; cr = np.array(list(controls.values()), dtype=int)
    gp = sums[cr].sum(0) / deps[cr].sum(); gm = qp[cr].sum(0) / n[cr].sum() * 1e6
    targets = np.array([t for t in targets if any(l == t for _, l in keys)])
    nt = len(targets); ng = len(genes); tc = np.zeros((nt, ng)); tm = np.zeros_like(tc); mp = np.zeros_like(tc); mm = np.zeros_like(tc); dv = np.zeros_like(tc); tn = np.zeros(nt, np.int64)
    for ti, t in enumerate(targets):
        rr = [i for i, (b, l) in enumerate(keys) if l == t]
        tn[ti] = round(n[rr].sum()); tc[ti] = sums[rr].sum(0); tm[ti] = qp[rr].sum(0) / n[rr].sum() * 1e6
        for i in rr:
            ci = controls.get(keys[i][0]); lw = deps[i] / deps[rr].sum(); cw = n[i] / n[rr].sum()
            cp = gp if ci is None else sums[ci] / deps[ci]; cm = gm if ci is None else qp[ci] / n[ci] * 1e6; mp[ti] += lw * cp; mm[ti] += cw * cm
            if ci is not None and not bulk:
                vp = np.maximum(q2[i] / n[i] - (qp[i] / n[i])**2, 0) / max(n[i] - 1, 1); vc = np.maximum(q2[ci] / n[ci] - (qp[ci] / n[ci])**2, 0) / max(n[ci] - 1, 1); dv[ti] += cw * cw * (vp + vc) * 1e12
            else:
                dv[ti] += cw * cw * (tm[ti] / max(tn[ti], 1) + cm / max(n[cr].sum(), 1))
    a.file.close()
    log(name, 'local bank', nt, 'targets', round(time.time() - start), 's')
    return dict(source=np.array(name), targets=targets.astype(str), genes=genes.astype(str), n_cells=tn, target_count_sums=tc.astype('f4'), target_mean_cpm=tm.astype('f4'),
                matched_control_probability=mp.astype('f4'), matched_control_mean_cpm=mm.astype('f4'), global_control_probability=gp.astype('f4'),
                global_control_mean_cpm=gm.astype('f4'), measured_genes=np.ones(ng, bool), delta_cpm_variance=dv.astype('f4'), arithmetic_cpm_observed=np.array(not bulk))


# ----------------------------------------------------------------------------------------------
# Orion (HCT116 / HEK293T) via Lance
# ----------------------------------------------------------------------------------------------
def orion_cells(name, labels):
    import lance
    d = lance.dataset(f'hf://datasets/slaf-project/X-Atlas-Orion/data/{name}/cells.lance', storage_options={'hf_token': ''})
    quoted = ','.join("'" + t.replace("'", "") + "'" for t in labels)
    return d.to_table(columns=['cell_integer_id', 'cell_start_index', 'n_genes_by_counts', 'total_counts', 'sample', 'gene_target', 'pass_guide_filter'],
                      filter=f'gene_target IN ({quoted})').to_pandas()


def orion_open(name):
    import lance
    from public_store import PublicStore
    from orion_offsets import offset_mapping
    p = PublicStore()
    meta = lance.dataset(f'hf://datasets/slaf-project/X-Atlas-Orion/data/{name}/expression.lance', storage_options={'hf_token': ''})
    e = lance.dataset(f's3://public/data/{name}/expression.lance', version=meta.version, storage_options=p.options(), block_size=4096)
    original, logical = offset_mapping(name, e, p)
    return p, e, original, logical, meta.version


def orion_take(e, original, logical, rr):
    """Fetch cells (DataFrame rows) -> (local cell index, gene_integer_id, value), integrity-checked."""
    length = rr.n_genes_by_counts.to_numpy(np.int64); gg = rr.cell_integer_id.to_numpy(np.int64) // 5000
    starts = rr.cell_start_index.to_numpy(np.int64) - original[gg] + logical[gg]
    ix = np.concatenate([np.arange(s, s + k, dtype=np.uint64) for s, k in zip(starts, length)])
    for attempt in range(5):
        try:
            tab = e.take(ix, columns=['cell_integer_id', 'gene_integer_id', 'value']); break
        except OSError:
            if attempt == 4:
                raise
            time.sleep(15)
    local = np.repeat(np.arange(len(rr)), length)
    assert np.array_equal(tab['cell_integer_id'].to_numpy(), np.repeat(rr.cell_integer_id.to_numpy(), length))
    val = tab['value'].to_numpy().astype(float); assert np.all(val == np.floor(val)) and np.all(val >= 0)
    np.testing.assert_array_equal(np.bincount(local, weights=val, minlength=len(rr)), rr.total_counts.to_numpy(float))
    return local, tab['gene_integer_id'].to_numpy(), val


def orion_bank(name, targets, ntc_per_sample=250, cells=None):
    t0 = time.time()
    if cells is None:
        cells = orion_cells(name, list(targets) + ['Non-Targeting'])
    cells = cells[cells.pass_guide_filter | (cells.gene_target == 'Non-Targeting')]
    v = pd.read_parquet(rt / f'artifacts/{name}_genes.parquet')
    h1 = ad.read_h5ad(runtime.VCC / 'public_data/H1.h5ad', backed='r'); h1genes = h1.var_names.to_numpy(str); h1.file.close()
    genes = np.array(list(dict.fromkeys(pd.read_csv(runtime.VCC / 'data/gene_names.csv').gene_name.tolist() + list(h1genes))))
    mapping = pd.Index(genes).get_indexer(v.gene_id.to_numpy(str))
    ntc = cells[cells.gene_target == 'Non-Targeting'].sort_values(['sample', 'cell_integer_id']).groupby('sample', observed=True).head(ntc_per_sample)
    c = pd.concat([cells[cells.gene_target.isin(targets)], ntc]).sort_values('cell_integer_id').reset_index(drop=True)
    p, e, original, logical, version = orion_open(name)
    c = c.loc[logical[c.cell_integer_id.to_numpy(np.int64) // 5000] >= 0].reset_index(drop=True)
    targets = np.array([t for t in targets if (c.gene_target == t).any()])
    samples = sorted(c['sample'].astype(str).unique()); bi = pd.Index(samples).get_indexer(c['sample'].astype(str)); ti = pd.Index(targets).get_indexer(c.gene_target.astype(str))
    acc = Acc(targets, samples, len(genes))
    log(name, 'Orion cells to fetch', len(c), 'targets', len(targets))
    try:
        for lo in range(0, len(c), 256):
            runtime.guard()
            rr = c.iloc[lo:lo + 256]; local, gid, val = orion_take(e, original, logical, rr)
            col = mapping[gid]; keep = col >= 0
            X = sparse.csr_matrix((val[keep], (local[keep], col[keep])), shape=(len(rr), len(genes)))
            d = np.asarray(X.sum(1)).ravel()
            good = d > 0
            groups = np.where(ti[lo:lo + 256] >= 0, ti[lo:lo + 256], len(targets) + bi[lo:lo + 256])
            acc.add_block(groups[good], ti[lo:lo + 256][good], bi[lo:lo + 256][good], X[good], d[good])
            if lo % 5120 == 0:
                log(name, lo, len(c), 'GB', round(p.bytes / 1e9, 3), 's', round(time.time() - t0))
    finally:
        p.close()
    acc.dump(rt / f'artifacts/tmp/acc_{name}.npz', genes)
    b = acc.bank(name, genes, measured=np.isin(np.arange(len(genes)), mapping[mapping >= 0]))
    return b, dict(cells=len(c), seconds=time.time() - t0, bytes=p.bytes, targets=len(targets), lance_version=version)


# ----------------------------------------------------------------------------------------------
# CD4 (GWCD4i pseudobulk + DE stats), model_v4/prepare_cd4.py parameterised
# ----------------------------------------------------------------------------------------------
def cd4_obs(kind):
    path = rt / f'artifacts/GWCD4i.{kind}.h5ad.obs.parquet'
    if not path.exists():
        from remote_h5 import Remote
        r = Remote(CD4_BASE + f'GWCD4i.{kind}.h5ad', cap=300_000_000)
        with h5py.File(r, 'r') as h:
            ad.io.read_elem(h['obs']).to_parquet(path)
        log('CD4 obs', kind, 'MB', r.bytes // 10**6)
    return pd.read_parquet(path)


def cd4_banks(targets):
    """Returns {name: bank dict} for CD4_<condition> pseudobulk banks and CD4_DE_statistics."""
    from remote_h5 import Remote
    t0 = time.time(); tset = set(targets); out = {}
    obs = cd4_obs('pseudobulk_merged'); var = pd.read_parquet(rt / 'artifacts/GWCD4i.pseudobulk_merged.h5ad.var.parquet')
    sel = obs.perturbed_gene_name.astype(str).isin(tset) & obs.keep_for_DE.astype(bool)
    ctrl = obs[obs.guide_type.astype(str) == 'non-targeting'].sort_index().groupby(['10xrun_id', 'donor_id', 'culture_condition'], observed=True).head(20)
    sel |= obs.index.isin(ctrl.index); rows = np.flatnonzero(sel); chosen = obs.iloc[rows].copy()
    r = Remote(CD4_BASE + 'GWCD4i.pseudobulk_merged.h5ad', cap=4_000_000_000, block=65536)
    with h5py.File(r, 'r') as h:
        g = h['X']; ptr = g['indptr'][:].astype(np.int64); lengths = ptr[rows + 1] - ptr[rows]; op = np.r_[0, np.cumsum(lengths)]; arrays = {}
        for kind in ['data', 'indices']:
            d = g[kind]; chunk = d.chunks[0]; needed = set()
            for row in rows:
                a, b = ptr[row:row + 2]
                if b > a:
                    needed.update(range(int(a) // chunk, (int(b) - 1) // chunk + 1))
            locations = {}

            def visit(info):
                number = int(info.chunk_offset[0]) // chunk
                if number in needed:
                    locations[number] = (info.byte_offset, info.size)
            d.id.chunk_iter(visit); plans = []
            for k, row in enumerate(rows):
                a = int(ptr[row]); b = int(ptr[row + 1]); dest = int(op[k])
                while a < b:
                    number = a // chunk; take = min(b, (number + 1) * chunk) - a; off = locations[number][0] + (a - number * chunk) * d.dtype.itemsize
                    plans.append((off, take * d.dtype.itemsize, dest, take)); a += take; dest += take
            groups = []
            for plan in sorted(plans):
                off, size, *_ = plan
                if groups and off <= groups[-1][1] + 16384 and max(groups[-1][1], off + size) - groups[-1][0] <= 4 * 1024**2:
                    groups[-1][1] = max(groups[-1][1], off + size); groups[-1][2].append(plan)
                else:
                    groups.append([off, off + size, [plan]])
            arr = np.empty(int(op[-1]), dtype=d.dtype); lock = threading.Lock()

            def fetch(group):
                a, b, members = group; s = requests.Session(); s.trust_env = False
                for attempt in range(5):
                    try:
                        with s.get(r.url, headers={'Range': f'bytes={a}-{b - 1}', 'If-Match': r.etag, 'Accept-Encoding': 'identity'}, timeout=(10, 90), stream=True) as resp:
                            assert resp.status_code == 206; buf = resp.raw.read(b - a + 1); assert len(buf) == b - a
                        with lock:
                            r.bytes += b - a
                        for off, size, dest, count in members:
                            arr[dest:dest + count] = np.frombuffer(buf, dtype=d.dtype, count=count, offset=off - a)
                        return
                    except Exception:
                        if attempt == 4:
                            raise
                        time.sleep(2)
            with ThreadPoolExecutor(max_workers=4) as pool:
                list(pool.map(fetch, groups))
            arrays[kind] = arr
        x = sparse.csr_matrix((arrays['data'], arrays['indices'], op), shape=(len(rows), len(var)))
    assert np.all(x.data >= 0) and np.all(x.data == np.floor(x.data)); actual = np.asarray(x.sum(1)).ravel()
    np.testing.assert_allclose(actual, chosen.total_counts.to_numpy(), rtol=1e-6)
    genes = var.gene_name.to_numpy(str)
    ptargets = np.array([t for t in targets if (chosen.perturbed_gene_name.astype(str) == t).any()])
    for cond in sorted(chosen.culture_condition.astype(str).unique()):
        c = chosen[chosen.culture_condition.astype(str) == cond]; keys = list(dict.fromkeys(zip(c['10xrun_id'].astype(str), c.donor_id.astype(str))))
        nt = len(ptargets); ng = len(genes); counts = np.zeros((nt, ng)); mp = np.zeros_like(counts); cm = np.zeros_like(counts); tm = np.zeros_like(counts); nc = np.zeros(nt, np.int64); dv = np.zeros_like(counts)
        indices = np.flatnonzero(chosen.culture_condition.astype(str).to_numpy() == cond)
        ctrl_rows = indices[chosen.iloc[indices].guide_type.astype(str).to_numpy() == 'non-targeting']; gcs = np.asarray(x[ctrl_rows].sum(0)).ravel(); gp = gcs / gcs.sum(); gm = gp * 1e6
        for ti, t in enumerate(ptargets):
            rr = indices[chosen.iloc[indices].perturbed_gene_name.astype(str).to_numpy() == t]
            if not len(rr):
                continue
            n = chosen.iloc[rr].n_cells.to_numpy(float); nc[ti] = int(n.sum()); counts[ti] = np.asarray(x[rr].sum(0)).ravel(); depth = actual[rr]; tm[ti] = counts[ti] / max(counts[ti].sum(), 1) * 1e6
            for run, donor in keys:
                pr = rr[(chosen.iloc[rr]['10xrun_id'].astype(str).to_numpy() == run) & (chosen.iloc[rr].donor_id.astype(str).to_numpy() == donor)]
                if not len(pr):
                    continue
                cr = ctrl_rows[(chosen.iloc[ctrl_rows]['10xrun_id'].astype(str).to_numpy() == run) & (chosen.iloc[ctrl_rows].donor_id.astype(str).to_numpy() == donor)]
                if not len(cr):
                    continue
                cp = np.asarray(x[cr].sum(0)).ravel(); cp /= cp.sum(); lw = actual[pr].sum() / depth.sum(); cw = chosen.iloc[pr].n_cells.sum() / n.sum(); mp[ti] += lw * cp; cm[ti] += cw * cp * 1e6
                pp = x[pr].multiply((1e6 / actual[pr])[:, None]).toarray(); cc = x[cr].multiply((1e6 / actual[cr])[:, None]).toarray()
                if len(pr) > 1:
                    dv[ti] += cw * cw * np.var(pp, axis=0, ddof=1) / len(pr)
                if len(cr) > 1:
                    dv[ti] += cw * cw * np.var(cc, axis=0, ddof=1) / len(cr)
        out[f'CD4_{cond}'] = dict(source=np.array('CD4_' + cond), targets=ptargets, genes=genes, n_cells=nc, target_count_sums=counts.astype('f4'), target_mean_cpm=tm.astype('f4'),
                                  matched_control_probability=mp.astype('f4'), matched_control_mean_cpm=cm.astype('f4'), global_control_probability=gp.astype('f4'),
                                  global_control_mean_cpm=gm.astype('f4'), measured_genes=np.ones(ng, bool), delta_cpm_variance=dv.astype('f4'), arithmetic_cpm_observed=np.array(False))
    dobs = cd4_obs('DE_stats'); dvar = pd.read_parquet(rt / 'artifacts/GWCD4i.DE_stats.h5ad.var.parquet')
    dr = np.flatnonzero(dobs.target_contrast_gene_name.astype(str).isin(tset))
    de = Remote(CD4_BASE + 'GWCD4i.DE_stats.h5ad', cap=2_000_000_000, block=65536)
    with h5py.File(de, 'r') as h:
        chosen_de = dobs.iloc[dr]; layers = {}
        for k in ['log_fc', 'lfcSE', 'adj_p_value']:
            ds = h['layers'][k]; assert ds.chunks is None
            layers[k] = np.stack([ds[int(row), :] for row in dr]) if len(dr) else np.zeros((0, ds.shape[1]), 'f4')
    dtargets = np.array([t for t in targets if (chosen_de.target_contrast_gene_name.astype(str) == t).any()])
    conditions = sorted(dobs.culture_condition.astype(str).unique()); shape = (len(conditions), len(dtargets), len(dvar))
    effects = np.zeros(shape, 'f4'); ses = np.zeros(shape, 'f4'); q = np.ones(shape, 'f4'); avail = np.zeros(shape[:2], bool); qc = np.zeros_like(avail); n = np.zeros(shape[:2], np.int64)
    for ci, cond in enumerate(conditions):
        for ti, t in enumerate(dtargets):
            rr = np.flatnonzero((chosen_de.culture_condition.astype(str).to_numpy() == cond) & (chosen_de.target_contrast_gene_name.astype(str).to_numpy() == t))
            if not len(rr):
                continue
            m = chosen_de.iloc[rr]; effects[ci, ti] = layers['log_fc'][rr].mean(0); ses[ci, ti] = layers['lfcSE'][rr].mean(0); q[ci, ti] = layers['adj_p_value'][rr].max(0)
            avail[ci, ti] = True; n[ci, ti] = int(m.n_cells_target.min())
            qc[ci, ti] = bool(((m.n_guides >= 2) & ~m.single_guide_estimate.astype(bool) & m.ontarget_significant.astype(bool) & ~m.distal_offtarget_flag.astype(bool) & ~m.low_target_gex.astype(bool)).all())
    out['CD4_DE_statistics'] = dict(targets=dtargets, genes=dvar.gene_name.to_numpy(str) if 'gene_name' in dvar else dvar.index.to_numpy(str), conditions=np.array(conditions),
                                    log2fc=effects, lfcSE=ses, adjusted_p=q, available=avail, quality_pass=qc, n_cells=n)
    return out, dict(seconds=time.time() - t0, bytes=r.bytes + de.bytes, pseudobulk_rows=len(rows), targets=len(ptargets), de_targets=len(dtargets))


def union_cd4_de(new, out_dir):
    with np.load(SHARED / 'CD4_DE_statistics.npz') as d:
        old = {k: d[k] for k in d.files}
    keep = ~np.isin(new['targets'], old['targets'])
    assert list(old['conditions']) == list(new['conditions']) and np.array_equal(old['genes'], new['genes'])
    m = dict(old); m['targets'] = np.concatenate([old['targets'], new['targets'][keep]])
    for k in ['log2fc', 'lfcSE', 'adjusted_p']:
        m[k] = np.concatenate([old[k], new[k][:, keep].astype(old[k].dtype)], axis=1)
    for k in ['available', 'quality_pass', 'n_cells']:
        m[k] = np.concatenate([old[k], new[k][:, keep].astype(old[k].dtype)], axis=1)
    savez(Path(out_dir) / 'CD4_DE_statistics.npz', **m)
    log('CD4_DE_statistics union', len(old['targets']), '+', int(keep.sum()))


# ----------------------------------------------------------------------------------------------
# promoter pairs (model_v4/prepare_promoters.py, parameterised)
# ----------------------------------------------------------------------------------------------
def promoter_pairs(targets, out_csv):
    p = rt / 'artifacts/CD4_guide_metadata.csv'
    if not p.exists():
        s = requests.Session(); s.trust_env = False
        with s.get(CD4_BASE + 'suppl_tables/sgrna_library_metadata.suppl_table.csv', stream=True, timeout=60) as r:
            r.raise_for_status(); value = r.raw.read(15_000_001); assert len(value) < 15_000_001; p.write_bytes(value)
    guides = pd.read_csv(p); rows = []
    with gzip.open(runtime.VCC / 'public_data/gencode.v47.basic.annotation.gtf.gz', 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            z = line.rstrip().split('\t')
            if z[2] != 'gene':
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]+)"', z[8])); name = attrs.get('gene_name')
            if name:
                rows.append(dict(gene=name, chromosome=z[0], tss=int(z[3] if z[6] == '+' else z[4]), strand=z[6]))
    a = pd.DataFrame(rows); a = a[~a.gene.duplicated(keep=False)].set_index('gene'); records = []
    old = pd.read_csv(rt / 'artifacts/promoter_pairs_v4.csv')
    for target in targets:
        if target not in a.index:
            continue
        t = a.loc[target]; near = a[(a.chromosome == t.chromosome) & ((a.tss - t.tss).abs() <= 5000)]
        for neighbor, b in near.iterrows():
            if neighbor == target:
                continue
            divergent = t.strand != b.strand and ((t.strand == '+' and t.tss > b.tss) or (t.strand == '-' and t.tss < b.tss))
            supported = guides[(guides.target_gene_name == target) & (guides.nearest_within2kb_nontarget_gene_name == neighbor)]
            distances = pd.to_numeric(supported.nearest_within2kb_nontarget_gene_dist, errors='coerce').abs().dropna(); n_guides = supported.sgRNA.nunique() if 'sgRNA' in supported else len(supported)
            records.append(dict(target=target, neighbor=neighbor, distance=abs(int(t.tss) - int(b.tss)), divergent=bool(divergent), guide_supported=bool(n_guides >= 2 and len(distances) >= 2),
                                guide_distance=float(distances.median()) if len(distances) else float('nan'), guide_count=int(n_guides), guide_provenance='public CD4 source library; challenge guides unavailable'))
    new = pd.DataFrame(records)
    pairs = pd.concat([old, new[~new.target.isin(old.target)]], ignore_index=True) if len(new) else old
    pairs.to_csv(out_csv, index=False)
    log('promoter pairs', len(pairs), 'rows,', pairs.target.nunique(), 'targets')
    return out_csv


# ----------------------------------------------------------------------------------------------
def build_panel_banks(targets, out_dir, sources=('K562_GWPS_bulk', 'H1', 'K562_GWPS_raw', 'HCT116', 'HEK293T', 'CD4', 'promoters')):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True); receipt = {}
    targets = [str(t) for t in targets]
    for s in sources:
        t0 = time.time()
        if s == 'promoters':
            promoter_pairs(targets, out_dir / 'promoter_pairs.csv'); receipt[s] = dict(seconds=time.time() - t0); continue
        if s == 'CD4':
            have = shared_targets('CD4_Rest'); need = [t for t in targets if t not in set(have)]
            if need:
                banks, info = cd4_banks(need)
                for name, b in banks.items():
                    if name == 'CD4_DE_statistics':
                        union_cd4_de(b, out_dir)
                    else:
                        union_bank(name, b, out_dir)
            else:
                info = dict(skipped='all targets in shared bank')
                for name in ['CD4_Rest', 'CD4_Stim8hr', 'CD4_Stim48hr', 'CD4_DE_statistics']:
                    tmp = out_dir / f'{name}.partial.npz'; tmp.write_bytes((SHARED / f'{name}.npz').read_bytes()); tmp.replace(out_dir / f'{name}.npz')
            receipt[s] = dict(info, need=len(need), seconds=time.time() - t0); continue
        have = set(shared_targets(s)); need = [t for t in targets if t not in have]
        if s in ('K562_GWPS_bulk', 'H1') and not (SHARED / f'{s}.npz').exists():
            b = local_bank(s, targets); savez(out_dir / f'{s}.npz', **b); info = dict(targets=len(b['targets']))
        elif not need:
            union_bank(s, None, out_dir); info = dict(skipped='all targets in shared bank')
        elif s == 'H1':
            b = local_bank('H1', need); info = dict(targets=len(b['targets']))
            union_bank('H1', b, out_dir)
        elif s == 'K562_GWPS_raw':
            b, info = gwps_raw_bank(need); union_bank(s, b, out_dir)
        elif s in ('HCT116', 'HEK293T'):
            b, info = orion_bank(s, need); union_bank(s, b, out_dir)
        else:
            raise ValueError(s)
        receipt[s] = dict(info, need=len(need), seconds=round(time.time() - t0, 1))
        log('BANK', s, json.dumps(receipt[s], default=str))
        runtime.save(out_dir / f'banks_receipt_{"_".join(sources)}.json', receipt)
    cov = {}; seen = set()
    for f in sorted(out_dir.glob('*.npz')):
        if f.stem.endswith('.partial'):
            continue
        try:
            with np.load(f) as d:
                tt = set(d['targets'].astype(str))
        except Exception:
            continue
        cov[f.stem] = int(sum(t in tt for t in targets)); seen |= tt
    receipt['panel_coverage'] = cov
    receipt['targets_without_any_source'] = [t for t in targets if t not in seen]
    runtime.save(out_dir / f'banks_receipt_{"_".join(sources)}.json', receipt)
    log('panel coverage', json.dumps(cov), 'no source:', len(receipt['targets_without_any_source']))
    return receipt


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--targets', required=True); ap.add_argument('--out', required=True)
    ap.add_argument('--sources', default='K562_GWPS_bulk,H1,K562_GWPS_raw,HCT116,HEK293T,CD4,promoters')
    a = ap.parse_args()
    t = pd.read_csv(a.targets); col = 'target_gene' if 'target_gene' in t else t.columns[0]
    build_panel_banks(sorted(set(t[col].astype(str)) - {'non-targeting'}), a.out, a.sources.split(','))
