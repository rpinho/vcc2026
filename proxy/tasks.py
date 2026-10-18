"""LOCO-300 proxy tasks (PIE_AND_EVIDENCE_PLAN fix #2), built outside rsi/ (rsi files are read, never written).

Panel per destination (response-blind; fixed before any model is scored):
  eligible = the 437-target v4 bank universe (official validation 300 + 137 H1-2025 targets) minus every
  Replogle-essential target (K562/RPE1/HepG2/Jurkat union), measured with >= MIN_CELLS real cells in the
  destination; ordered by sha256('proxy:<dest>:<target>'); first 300 kept. Only cell counts, essentiality and the
  hash are used. The official validation panel is itself the VCC-matched draw (0 essential, ~55% responders), so the
  panel inherits that distribution; responder share is REPORTED from K562 GWPS energy tests (not used to select).
Reference cells: all real cells up to 400 per target (seeded hash sample). Controls: NTC cells allocated to batches
  (gem group / sample / lane) in proportion to the panel's own target cells, 6,000 total. As on the server, the SAME
  control cells are the model input and the scorer's reference (the official manifest counts the 18,400 released
  controls inside the ground truth).
Depth: if the destination's median UMI/cell exceeds 20k, every cell (targets and controls) is binomially thinned by
  20k/median. Gene axis: measured genes intersected with the official 18,533 (official order), never zero-filled.

python tasks.py K562|HCT116|HEK293T|CD4
"""
import sys, os, json, hashlib, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'model_v6'))
import runtime  # model_v6 runtime (thread caps, nice, guard)
import numpy as np, pandas as pd, anndata as ad, h5py, requests, threading
from pathlib import Path
from scipy import sparse
import banks

P = Path(__file__).resolve().parent; V = runtime.VCC
MIN_CELLS = 64; MAX_CELLS = 400; N_CTRL = 6000; DEPTH = 20000
OFFICIAL = pd.read_csv(V / 'data/gene_names.csv').gene_name.astype(str).to_numpy()


def hkey(s):
    return hashlib.sha256(s.encode()).hexdigest()


def seed(s):
    return int.from_bytes(hashlib.sha256(s.encode()).digest()[:8], 'little')


def essential():
    f = P / 'artifacts/essential_union.json'
    if not f.exists():
        ess = set()
        for n in ['K562', 'RPE1', 'HepG2', 'Jurkat']:
            a = ad.read_h5ad(V / f'public_data/{n}.h5ad', backed='r'); col = 'gene' if 'gene' in a.obs else 'target_gene'
            ess |= set(a.obs[col].astype(str)); a.file.close()
        runtime.save(f, sorted(ess))
    return set(json.loads(f.read_text()))


def universe():
    return np.load(V / 'public_data/model_v4_signature_banks/HCT116.npz')['targets'].astype(str)


def gwps_responders():
    import re
    b = ad.read_h5ad(V / 'public_data/K562_GWPS_bulk.h5ad', backed='r'); o = b.obs.copy(); b.file.close()
    o['t'] = [re.match(r'^[^_]+_(.+)_[^_]+_([^_]+)$', s).group(1) if 'non-targeting' not in s else 'non-targeting' for s in o.index]
    return o.groupby('t').energy_test_p_value.min()


def select_panel(dest, counts):
    """counts: Series target -> available real cells in destination."""
    ess = essential(); official = set(pd.read_csv(V / 'data/pert_counts.csv').target_gene.astype(str))
    elig = [t for t in universe() if t not in ess and counts.get(t, 0) >= MIN_CELLS]
    panel = sorted(elig, key=lambda t: hkey(f'proxy:{dest}:{t}'))[:300]
    r = gwps_responders()
    resp = [t for t in panel if t in r.index]
    info = dict(destination=dest, eligible=len(elig), n=len(panel), official_val_targets=int(sum(t in official for t in panel)),
                essential=0, gwps_measured=len(resp), gwps_responder_share=float(np.mean([r[t] < .05 for t in resp])) if resp else None,
                real_cells_median=float(np.median([min(counts[t], MAX_CELLS) for t in panel])), rule=__doc__.split('Panel per destination')[1].split('Reference cells')[0].strip())
    return panel, info


def allocate_controls(target_batches, ctrl_df, batch_col, n_total=N_CTRL, rng=None):
    """NTC rows proportional to the panel's target cells per batch (largest remainder), capped by availability."""
    share = target_batches.value_counts(normalize=True)
    want = (share * n_total).round().astype(int)
    rows = []
    for b, k in want.items():
        pool = ctrl_df[ctrl_df[batch_col].astype(str) == str(b)]
        if len(pool) and k:
            rows.append(pool.iloc[rng.choice(len(pool), min(k, len(pool)), replace=False)])
    return pd.concat(rows) if rows else ctrl_df.iloc[:0]


def finish(dest, X, obs, genes_src, panel, info):
    """Map to measured official axis, thin depth, write task h5ad + receipt."""
    sparse.save_npz(P / 'tasks' / f'_{dest}_raw_fetch.npz', sparse.csr_matrix(X)); obs.to_parquet(P / 'tasks' / f'_{dest}_raw_obs.parquet')
    genes_src = np.asarray(genes_src).astype(str)
    present = set(genes_src); genes = np.array([g for g in OFFICIAL if g in present])
    if len(set(genes_src)) == len(genes_src):
        ix = pd.Index(genes_src).get_indexer(genes); keep = np.ones(len(genes), bool)
        X = sparse.csr_matrix(X)[:, ix]
    else:
        mapped = pd.Index(genes).get_indexer(genes_src); valid = mapped >= 0
        X = sparse.csr_matrix(X) @ sparse.csr_matrix((np.ones(valid.sum()), (np.flatnonzero(valid), mapped[valid])), shape=(len(genes_src), len(genes)))
    X = sparse.csr_matrix(X); X.data = np.rint(X.data); X = X.astype(np.int64)
    depth = np.asarray(X.sum(1)).ravel(); med = float(np.median(depth))
    thin = 1.0
    if med > DEPTH:
        thin = DEPTH / med; rng = np.random.default_rng(seed('thin:' + dest))
        X.data = rng.binomial(X.data, thin)
    X.eliminate_zeros(); X = X.astype(np.int32)
    depth = np.asarray(X.sum(1)).ravel(); ok = depth > 0
    X = X[ok]; obs = obs.iloc[np.flatnonzero(ok)].reset_index(drop=True)
    obs.index = [f'{dest}_{i}' for i in range(len(obs))]
    a = ad.AnnData(X, obs=obs, var=pd.DataFrame(index=genes)); a.uns['evaluation_only'] = True
    path = P / 'tasks' / f'{dest}.h5ad'; a.write_h5ad(path, compression='gzip', compression_opts=1)
    lab = obs.target_gene.astype(str)
    info.update(genes=int(len(genes)), native_median_umi=med, thinning_factor=thin, median_umi=float(np.median(depth[ok])),
                cells_total=int(a.n_obs), controls=int((lab == 'non-targeting').sum()),
                cells_per_target=lab[lab != 'non-targeting'].value_counts().describe().to_dict(), sha256=runtime.sha(path))
    pd.DataFrame({'target_gene': panel}).to_csv(P / 'tasks' / f'{dest}_panel.csv', index=False)
    runtime.save(P / 'tasks' / f'{dest}.json', info)
    print('TASK', dest, json.dumps({k: v for k, v in info.items() if k != 'rule'}), flush=True)


def build_gwps():
    dest = 'K562'; obs, var = banks.gwps_meta(); lab = obs.gene.astype(str)
    counts = lab.value_counts(); panel, info = select_panel(dest, counts)
    rows = []
    for t in panel:
        r = np.flatnonzero(lab.to_numpy() == t); rng = np.random.default_rng(seed(f'cells:{dest}:{t}'))
        rows.extend(rng.choice(r, min(MAX_CELLS, len(r)), replace=False))
    tsel = obs.iloc[np.sort(rows)]
    ctrl = allocate_controls(tsel.gem_group, obs[lab == 'non-targeting'], 'gem_group', rng=np.random.default_rng(seed('ctrl:' + dest)))
    sel = pd.concat([tsel, ctrl]).sort_values('source_row')
    offset, shape = banks.gwps_layout(); genes0 = var.gene_name.to_numpy(str)
    rr = sel.source_row.to_numpy(np.int64); blocks = []; buf = []
    for k, x in banks.gwps_fetch(rr, len(genes0), offset, cap_bytes=4_000_000_000):
        buf.append(x)
        if len(buf) == 1024:
            blocks.append(sparse.csr_matrix(np.stack(buf))); buf = []
    if buf:
        blocks.append(sparse.csr_matrix(np.stack(buf)))
    X = sparse.vstack(blocks, format='csr'); assert X.shape[0] == len(sel)
    assert np.all(X.data == np.floor(X.data))
    o = pd.DataFrame({'target_gene': sel.gene.astype(str).to_numpy(), 'batch': sel.gem_group.astype(str).to_numpy(),
                      'guide_id': np.where(sel.gene.astype(str) == 'non-targeting', sel.sgID_AB.astype(str), sel.sgID_AB.astype(str)), 'source_row': rr})
    info.update(source='Replogle 2022 K562 GWPS raw (figshare 35775507)', bytes=getattr(banks.gwps_fetch, 'last_bytes', None))
    finish(dest, X, o, genes0, panel, info)


def build_orion(dest):
    c = pd.read_parquet(V / f'model_v4/artifacts/{dest}_cells.parquet')
    p, e, original, logical, version = banks.orion_open(dest)
    c = c[logical[c.cell_integer_id.to_numpy(np.int64) // 5000] >= 0]
    tg = c[c.pass_guide_filter & (c.gene_target != 'Non-Targeting')]
    counts = tg.gene_target.astype(str).value_counts(); panel, info = select_panel(dest, counts)
    parts = []
    for t in panel:
        g = tg[tg.gene_target.astype(str) == t]; rng = np.random.default_rng(seed(f'cells:{dest}:{t}'))
        parts.append(g.iloc[rng.choice(len(g), min(MAX_CELLS, len(g)), replace=False)])
    tsel = pd.concat(parts)
    ctrl = allocate_controls(tsel['sample'], c[c.gene_target == 'Non-Targeting'], 'sample', rng=np.random.default_rng(seed('ctrl:' + dest)))
    sel = pd.concat([tsel, ctrl]).sort_values('cell_integer_id').reset_index(drop=True)
    v = pd.read_parquet(V / f'model_v4/artifacts/{dest}_genes.parquet'); genes0 = v.gene_id.to_numpy(str)
    blocks = []; t0 = time.time()
    try:
        for lo in range(0, len(sel), 256):
            runtime.guard()
            rr = sel.iloc[lo:lo + 256]; local, gid, val = banks.orion_take(e, original, logical, rr)
            blocks.append(sparse.csr_matrix((val, (local, gid)), shape=(len(rr), len(genes0))))
            if lo % 5120 == 0:
                print(dest, lo, len(sel), 'GB', round(p.bytes / 1e9, 3), 's', round(time.time() - t0), flush=True)
    finally:
        p.close()
    X = sparse.vstack(blocks, format='csr')
    lab = np.where(sel.gene_target.astype(str) == 'Non-Targeting', 'non-targeting', sel.gene_target.astype(str))
    o = pd.DataFrame({'target_gene': lab, 'batch': sel['sample'].astype(str).to_numpy(), 'source_row': sel.cell_integer_id.to_numpy()})
    info.update(source='X-Atlas-Orion (HF slaf-project), lance version %d' % version, bytes=p.bytes)
    finish(dest, X, o, genes0, panel, info)


def h5_csr_rows(url, rows, cap=6_000_000_000, workers=12):
    """Exact public byte ranges for selected CSR rows of a remote uncompressed-chunk h5ad (prepare_cd4 pattern)."""
    from remote_h5 import Remote
    r = Remote(url, cap=cap, block=65536)
    with h5py.File(r, 'r') as h:
        var = ad.io.read_elem(h['var']); g = h['X']; ptr = g['indptr'][:].astype(np.int64)
        rows = np.asarray(rows, np.int64); lengths = ptr[rows + 1] - ptr[rows]; op = np.r_[0, np.cumsum(lengths)]; arrays = {}
        for kind in ['data', 'indices']:
            d = g[kind]; assert d.compression is None; chunk = d.chunks[0]; needed = set()
            for row in rows:
                a, b = ptr[row:row + 2]
                if b > a:
                    needed.update(range(int(a) // chunk, (int(b) - 1) // chunk + 1))
            loc = {}

            def visit(info):
                num = int(info.chunk_offset[0]) // chunk
                if num in needed:
                    loc[num] = (info.byte_offset, info.size)
            d.id.chunk_iter(visit); plans = []
            for k, row in enumerate(rows):
                a = int(ptr[row]); b = int(ptr[row + 1]); dest = int(op[k])
                while a < b:
                    num = a // chunk; take = min(b, (num + 1) * chunk) - a; off = loc[num][0] + (a - num * chunk) * d.dtype.itemsize
                    plans.append((off, take * d.dtype.itemsize, dest, take)); a += take; dest += take
            groups = []
            for plan in sorted(plans):
                off, size, *_ = plan
                if groups and off <= groups[-1][1] + 16384 and max(groups[-1][1], off + size) - groups[-1][0] <= 2 * 1024**2:
                    groups[-1][1] = max(groups[-1][1], off + size); groups[-1][2].append(plan)
                else:
                    groups.append([off, off + size, [plan]])
            arr = np.empty(int(op[-1]), dtype=d.dtype); lock = threading.Lock(); local = threading.local()
            print(kind, 'ranges', len(groups), 'GB', round(sum(b - a for a, b, _ in groups) / 1e9, 3), flush=True)

            def fetch(group):
                a, b, members = group
                if not hasattr(local, 's'):
                    local.s = requests.Session(); local.s.trust_env = False
                for attempt in range(5):
                    try:
                        with local.s.get(r.url, headers={'Range': f'bytes={a}-{b - 1}', 'If-Match': r.etag, 'Accept-Encoding': 'identity'}, timeout=(10, 90), stream=True) as resp:
                            assert resp.status_code == 206; buf = resp.raw.read(b - a + 1); assert len(buf) == b - a
                        with lock:
                            r.bytes += b - a
                        for off, size, dest, count in members:
                            arr[dest:dest + count] = np.frombuffer(buf, dtype=d.dtype, count=count, offset=off - a)
                        return
                    except Exception:
                        if attempt == 4:
                            raise
                        time.sleep(2 + attempt)
            from concurrent.futures import ThreadPoolExecutor
            with ThreadPoolExecutor(max_workers=workers) as pool:
                for i, _ in enumerate(pool.map(fetch, groups)):
                    if i % 2000 == 0:
                        print(kind, i, len(groups), 'GB', round(r.bytes / 1e9, 3), flush=True)
            arrays[kind] = arr
    X = sparse.csr_matrix((arrays['data'], arrays['indices'], op), shape=(len(rows), len(var)))
    assert np.all(X.data >= 0) and np.all(X.data == np.floor(X.data))
    names = var.gene_name.astype(str).to_numpy() if 'gene_name' in var else var.index.to_numpy(str)
    return X, names, r.bytes


def build_cd4():
    dest = 'CD4'; url = 'https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/D4_Stim8hr.assigned_guide.h5ad'
    c = pd.read_parquet(V / 'rsi/artifacts/CD4_cell_census.parquet'); c = c[~c.low_quality.astype(bool)]
    tg = c[c.guide_type.astype(str) == 'targeting']
    counts = tg.target_gene.astype(str).value_counts(); panel, info = select_panel(dest, counts)
    parts = []
    for t in panel:
        g = tg[tg.target_gene.astype(str) == t]; rng = np.random.default_rng(seed(f'cells:{dest}:{t}'))
        parts.append(g.iloc[rng.choice(len(g), min(MAX_CELLS, len(g)), replace=False)])
    tsel = pd.concat(parts)
    ctrl = allocate_controls(tsel.lane_id, c[c.guide_type.astype(str) == 'non-targeting'], 'lane_id', rng=np.random.default_rng(seed('ctrl:' + dest)))
    sel = pd.concat([tsel, ctrl]).sort_values('source_row').reset_index(drop=True)
    X, names, nbytes = h5_csr_rows(url, sel.source_row.to_numpy(np.int64))
    tot = np.asarray(X.sum(1)).ravel(); np.testing.assert_allclose(tot, sel.total_counts.to_numpy(float), rtol=1e-6)
    lab = np.where(sel.guide_type.astype(str) == 'non-targeting', 'non-targeting', sel.target_gene.astype(str))
    o = pd.DataFrame({'target_gene': lab, 'batch': sel.lane_id.astype(str).to_numpy(), 'guide_id': sel.guide_id.astype(str).to_numpy(), 'source_row': sel.source_row.to_numpy()})
    info.update(source=url, bytes=nbytes)
    finish(dest, X, o, names, panel, info)


if __name__ == '__main__':
    d = sys.argv[1]
    (P / 'tasks').mkdir(exist_ok=True)
    {'K562': build_gwps, 'HCT116': lambda: build_orion('HCT116'), 'HEK293T': lambda: build_orion('HEK293T'), 'CD4': build_cd4}[d]()
