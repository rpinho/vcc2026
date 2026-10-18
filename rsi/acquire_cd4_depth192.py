"""Materialize the frozen79-target depth192 CD4 task, retaining parent cells.

The new1.2GB allocation is separate from both immutable earlier CD4 journals.
Metadata and raw payload checkpoints survive failures. No manifest activation.
Adapted from acquire_cd4_increment.py; the historical reader is unchanged.
"""
import runtime as rt
from transfer_budget import TransferBudget
from download_checkpoint import DownloadCheckpoint
from tree_prefetch import PrefetchRemote as Remote, read_uncompressed_vector
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib
import threading
import time
import requests
import h5py
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

NAME = 'CD4_depth192_search'
BASE = 'CD4_depth120_search'
URL = 'https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/D4_Stim8hr.assigned_guide.h5ad'


def main():
    folder = Path('artifacts/tasks')
    blueprint = folder/f'{NAME}_blueprint.json'
    plan = rt.json.loads(blueprint.read_text())
    audit = rt.json.loads(Path('artifacts/CD4_depth192_design_audit.json').read_text())
    assert audit['status'] == 'PASS independent count-only depth192 design audit'
    assert audit['blueprint_sha256'] == rt.sha(blueprint)
    fixture = rt.json.loads(Path('artifacts/cd4_depth192_fixture.json').read_text())
    assert fixture['status'] == 'PASS' and fixture['reader_sha256'] == rt.sha(__file__)
    path = folder/f'{NAME}.h5ad'
    checkpoint_root = Path('artifacts/CD4_depth192_payload_checkpoint')
    if path.exists() and not (checkpoint_root/'state.json').exists():
        raise SystemExit('CD4 depth extension already exists')
    selection_path = folder/f'{NAME}_selection.parquet'
    increment_path = folder/f'{NAME}_increment.parquet'
    assert rt.sha(selection_path) == plan['selection_sha256']
    assert rt.sha(increment_path) == plan['increment_sha256']
    selection = pd.read_parquet(selection_path).sort_values('source_row').reset_index(drop=True)
    extra = pd.read_parquet(increment_path).sort_values('source_row').reset_index(drop=True)
    assert len(extra) == audit['additional_cells'] and extra.source_row.is_unique
    assert rt.sha(folder/f'{BASE}.h5ad') == plan['parent_task_sha256']
    original_total = 0
    for journal, expected_hash in plan['previous_charge_journals'].items():
        assert rt.sha(journal) == expected_hash
        previous = [rt.json.loads(line) for line in Path(journal).read_text().splitlines()]
        original_total += sum(row['bytes'] for row in previous if row['event'] == 'charge')
    allocation = rt.json.loads(Path('artifacts/download_allocations.json').read_text())
    assert allocation['new_allocations']['CD4_depth192'] == plan['increment_ceiling_bytes']
    assert sum(allocation['new_allocations'].values()) <= allocation['limit_bytes'] == 15_000_000_000
    fixture = rt.json.loads(Path('artifacts/tree_prefetch_fixture.json').read_text())
    assert fixture['status'] == 'PASS bitwise index addresses and data'
    assert fixture['implementation_sha256'] == rt.sha('tree_prefetch.py')
    expected = rt.json.loads(Path('artifacts/CD4_cell_feasibility.json').read_text())
    source = dict(url=URL, etag=expected['etag'], size=expected['source_size'])
    budget = TransferBudget('artifacts/CD4_depth192_charges.jsonl', plan['increment_ceiling_bytes'])
    inputs = dict(blueprint=rt.sha(blueprint), increment=rt.sha(increment_path),
                  selection=rt.sha(selection_path), implementation=rt.sha(__file__),
                  census=rt.sha('artifacts/CD4_cell_census.parquet'),
                  base=plan['parent_task_sha256'], previous_journals=plan['previous_charge_journals'],
                  tree_prefetch=rt.sha('tree_prefetch.py'), source=source)
    budget.bind('frozen_inputs', inputs)
    start = time.time()
    cache_path = Path('artifacts/CD4_depth192_index.npz')
    cache_receipt = cache_path.with_suffix('.json')
    rows = extra.source_row.to_numpy(np.int64)
    if cache_receipt.exists():
        cache = rt.json.loads(cache_receipt.read_text())
        assert cache['inputs'] == inputs and cache['sha256'] == rt.sha(cache_path)
        metadata_bytes = cache['metadata_read_bytes']
        print('Reuse committed sparse index; no metadata request', flush=True)
    else:
        # All failed metadata requests remain charged inside this invocation reserve.
        budget.charge(251_000_000, 'metadata invocation upper bound including retries and probe')
        remote = Remote(URL, cap=250_000_000, block=65536)
        assert remote.etag == source['etag'] and remote.size == source['size']
        locations = {}
        with h5py.File(remote, 'r') as handle:
            var = ad.io.read_elem(handle['var'])
            genes = var.gene_name.astype(str).to_numpy() if 'gene_name' in var else var.index.to_numpy(str)
            group = handle['X']
            ptr = read_uncompressed_vector(group['indptr'], remote).astype(np.int64)
            remote.enable_tree_prefetch()
            packed = dict(genes=genes.astype(str), starts=ptr[rows], ends=ptr[rows+1])
            for kind in ['data', 'indices']:
                ds = group[kind]
                assert ds.compression is None and not ds.shuffle
                chunk = ds.chunks[0]
                needed = set()
                for a, b in zip(ptr[rows], ptr[rows+1]):
                    if b > a:
                        needed.update(range(int(a)//chunk, (int(b)-1)//chunk+1))
                loc = {}
                def visit(info):
                    number = int(info.chunk_offset[0])//chunk
                    if number in needed:
                        loc[number] = (info.byte_offset, info.size)
                ds.id.chunk_iter(visit)
                assert needed <= loc.keys()
                packed[kind+'_locations'] = np.array([[n, *loc[n]] for n in sorted(loc)], dtype=np.int64)
                packed[kind+'_dtype'] = np.array(ds.dtype.str)
                packed[kind+'_chunk'] = np.array(chunk)
                print('Indexed incremental', kind, 'metadata MB', round(remote.bytes/1e6, 2), flush=True)
        metadata_bytes = remote.bytes
        remote.http.close()
        temporary = cache_path.with_name(cache_path.stem+'.partial.npz')
        np.savez_compressed(temporary, **packed)
        temporary.replace(cache_path)
        rt.save(cache_receipt, dict(inputs=inputs, sha256=rt.sha(cache_path), metadata_read_bytes=metadata_bytes))
    budget.bind('index_cache', dict(sha256=rt.sha(cache_path), receipt_sha256=rt.sha(cache_receipt)))
    with np.load(cache_path) as cached:
        genes = cached['genes'].copy()
        starts, ends = cached['starts'].copy(), cached['ends'].copy()
        assert len(starts) == len(rows) and np.all(ends >= starts)
        outptr = np.r_[0, np.cumsum(ends-starts)]
        dtypes = {k: np.dtype(str(cached[k+'_dtype'])) for k in ['data', 'indices']}
        plans = []
        for kind in ['data', 'indices']:
            chunk = int(cached[kind+'_chunk'])
            locations = {int(n): (int(off), int(size)) for n, off, size in cached[kind+'_locations']}
            pieces = []
            for i, (begin, end) in enumerate(zip(starts, ends)):
                a, b, dest = int(begin), int(end), int(outptr[i])
                while a < b:
                    number = a//chunk
                    take = min(b, (number+1)*chunk)-a
                    offset = locations[number][0]+(a-number*chunk)*dtypes[kind].itemsize
                    pieces.append((offset, take*dtypes[kind].itemsize, dest, take, kind))
                    a += take
                    dest += take
            for piece in sorted(pieces):
                offset, size, *_ = piece
                if plans and plans[-1][3] == kind and offset <= plans[-1][1]+8192 and max(offset+size, plans[-1][1])-plans[-1][0] <= 2*1024**2:
                    plans[-1][1] = max(offset+size, plans[-1][1])
                    plans[-1][2].append(piece)
                else:
                    plans.append([offset, offset+size, [piece], kind])
    identity = dict(inputs=inputs, index_sha256=rt.sha(cache_path),
                    ranges_sha256=hashlib.sha256(np.asarray([(a, b) for a, b, *_ in plans], dtype=np.int64).tobytes()).hexdigest(),
                    outptr_sha256=hashlib.sha256(outptr.tobytes()).hexdigest())
    budget.bind('selection_and_plan', identity)
    rt.guard(int(outptr[-1])*sum(v.itemsize for v in dtypes.values()))
    checkpoint = DownloadCheckpoint(checkpoint_root, {k: ((int(outptr[-1]),), v) for k, v in dtypes.items()}, identity)
    arrays = checkpoint.arrays
    resume = checkpoint.next_plan
    assert resume <= len(plans)
    remaining = sum(b-a for a, b, *_ in plans[resume:])
    assert budget.total+remaining <= budget.limit, 'Increment cannot finish inside separate depth192 allocation'
    local = threading.local()
    def fetch(plan_range):
        a, b, pieces, kind = plan_range
        if not hasattr(local, 'http'):
            local.http = requests.Session()
            local.http.trust_env = False
        for attempt in range(3):
            rt.guard()
            budget.charge(b-a, 'increment payload range attempt')
            try:
                with local.http.get(URL, headers={'Range': f'bytes={a}-{b-1}', 'Accept-Encoding': 'identity', 'If-Match': source['etag']}, stream=True, timeout=(10, 60)) as response:
                    assert response.status_code == 206
                    buffer = response.raw.read(b-a+1)
                    assert len(buffer) == b-a
                for offset, size, dest, take, key in pieces:
                    arrays[key][dest:dest+take] = np.frombuffer(buffer, dtype=arrays[key].dtype, count=take, offset=offset-a)
                return
            except Exception:
                if attempt == 2:
                    raise
    receipt = Path('artifacts/CD4_depth192_transfer.json')
    print('Increment', len(rows), 'cells,', len(plans), 'ranges, resume', resume,
          'remaining bytes', remaining, 'ceiling', budget.limit, flush=True)
    try:
        with ThreadPoolExecutor(max_workers=3) as executor:
            for lo in range(resume, len(plans), 64):
                hi = min(lo+64, len(plans))
                for _ in executor.map(fetch, plans[lo:hi]):
                    pass
                checkpoint.commit(hi)
                if lo//64 % 8 == 0 or hi == len(plans):
                    print('Increment ranges', hi, '/', len(plans), 'charged', budget.total, flush=True)
    finally:
        rt.save(receipt, dict(bytes=budget.total, original_charged_bytes=original_total,
                             cumulative_cd4_bytes=original_total+budget.total,
                             ceiling_bytes=2_600_000_000+plan['increment_ceiling_bytes'], source=source,
                             metadata_read_bytes=metadata_bytes, resumed_at_plan=resume,
                             seconds=time.time()-start, charge_journal=str(budget.path),
                             status='ranges ended; task publication/audit separate'))
    matrix = sparse.csr_matrix((arrays['data'], arrays['indices'], outptr), shape=(len(rows), len(genes)))
    del arrays
    for lo in range(0, len(matrix.data), 1_000_000):
        values = matrix.data[lo:lo+1_000_000]
        assert np.isfinite(values).all() and np.all(values >= 0) and np.equal(values, np.floor(values)).all()
    matrix.data = matrix.data.astype(np.int32)
    official = pd.read_csv(rt.ROOT.parent/'data/gene_names.csv').iloc[:, 0].astype(str).to_numpy()
    measured = set(genes)
    readouts = np.array([gene for gene in official if gene in measured])
    mapped = pd.Index(readouts).get_indexer(genes)
    valid = mapped >= 0
    projection = sparse.csr_matrix((np.ones(valid.sum(), np.int32), (np.flatnonzero(valid), mapped[valid])), shape=(len(genes), len(readouts)))
    matrix = (matrix@projection).tocsr()
    matrix.sum_duplicates()
    matrix.sort_indices()
    assert np.all(np.asarray(matrix.sum(axis=1)).ravel() > 0)
    base = ad.read_h5ad(folder/f'{BASE}.h5ad')
    assert np.array_equal(base.var_names.to_numpy(str), readouts)
    base = base[base.obs.source_row.isin(selection.source_row)].copy()
    columns = ['target_gene', 'batch', 'guide_id', 'evaluation_role', 'source_row']
    joined = pd.concat([base.obs[columns].reset_index(drop=True), extra[columns]], ignore_index=True)
    assert joined.source_row.is_unique
    order = pd.Index(joined.source_row).get_indexer(selection.source_row)
    assert np.all(order >= 0) and len(joined) == len(selection)
    for column in columns:
        assert np.array_equal(joined.iloc[order][column].astype(str).to_numpy(), selection[column].astype(str).to_numpy()), column
    full = sparse.vstack([base.X, matrix], format='csr')[order]
    result = ad.AnnData(full, obs=selection[columns].copy(), var=pd.DataFrame(index=readouts))
    result.uns['evaluation_only'] = True
    rt.guard(2*1024**3)
    temporary = path.with_suffix('.partial.h5ad')
    result.write_h5ad(temporary, compression='gzip', compression_opts=1)
    temporary.replace(path)
    rt.save(path.with_suffix('.json'), dict(sha256=rt.sha(path), source=source,
                                          bytes=budget.total, cumulative_cd4_bytes=original_total+budget.total,
                                          parent_task_sha256=plan['parent_task_sha256'],
                                          selection_sha256=plan['selection_sha256'], blueprint_sha256=rt.sha(blueprint)))
    # Keep the frozen manifest intact until an independent live count audit passes.
    rt.save(folder/f'{NAME}_manifest_extension.json', {NAME: plan})
    checkpoint.remove()
    budget.close()
    print('DONE', NAME, result.shape, 'cumulative CD4 bytes', original_total+budget.total, flush=True)


if __name__ == '__main__':
    main()
