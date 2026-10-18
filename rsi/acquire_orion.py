"""Restore frozen Orion rows with durable range charges and payload reuse."""
import runtime as rt
from pathlib import Path
import sys
import time

import anndata as ad
import lance
import numpy as np
import pandas as pd
from scipy import sparse

from download_checkpoint import DownloadCheckpoint
from durable_public_store import PublicStore
from transfer_budget import TransferBudget


def bind_inputs(budget, inputs, checkpoint_root):
    """Accept only the observed pre-payload metadata-decoding repair.

    Original frozen bindings and failed-request charges remain append-only.
    No source or payload identity may be migrated by this exception.
    """
    old = budget.bindings.get('frozen_inputs')
    if old is None or old == inputs:
        budget.bind('frozen_inputs', inputs)
        return
    key = 'metadata_decode_fix_v1'
    migration = dict(prior_inputs=old, inputs=inputs,
                     reason='Bounded HTTP metadata content decoding; no payload existed before repair')
    if key not in budget.bindings:
        expected = dict(inputs,
            implementation='be119c8701c38551f9621a7969a84d0037ca56b3ee757413045f2dc0bcc3faf2',
            bridge='a5d48a092335c7198461c73a5940d7cdc5f7d3a343497cd13d0cf376ce010396')
        assert old == expected, 'Unrecognized pre-repair implementation or changed inputs'
        assert not checkpoint_root.exists(), 'Cannot migrate an existing payload checkpoint'
        assert set(budget.bindings) == {'limit_bytes', 'frozen_inputs'}, 'Source already bound'
        charges = [rt.json.loads(line) for line in budget.path.read_text().splitlines()
                   if rt.json.loads(line)['event'] == 'charge']
        assert charges == [dict(event='charge', bytes=1_000_000, total=1_000_000,
                              reason='repository revision metadata upper bound')]
    budget.bind(key, migration)


def main(name):
    folder = Path('artifacts/tasks')
    manifest = rt.json.loads((folder / 'manifest.json').read_text())
    task = next(k for k, value in manifest.items() if value['context'] == name)
    path = folder / (task + '.h5ad')
    checkpoint_root = Path(f'artifacts/{name}_payload_checkpoint')
    if path.exists() and not (checkpoint_root / 'state.json').exists():
        raise SystemExit('Orion destination exists')
    selection = folder / (task + '_selection.parquet')
    sel = pd.read_parquet(selection).sort_values('cell_integer_id').reset_index(drop=True)
    assert sel.cell_integer_id.is_unique
    assert rt.sha(selection) == manifest[task]['selection_sha256']
    gene_path = Path(f'artifacts/{name}_genes.parquet')
    offsets_path = Path(f'artifacts/{name}_offset_mapping.npz')
    var = pd.read_parquet(gene_path)
    off = np.load(offsets_path)
    version = int(off['version'])
    genes = var.gene_id.to_numpy(str)
    unique = np.array(list(dict.fromkeys(genes)))
    mapping = pd.Index(unique).get_indexer(genes)
    lengths = sel.n_genes_by_counts.to_numpy(np.int64)
    assert (lengths > 0).all()
    ptr = np.r_[0, np.cumsum(lengths)]
    receipt = Path(f'artifacts/{name}_raw_transfer.json')
    journal = Path(f'artifacts/{name}_charges.jsonl')
    # Never mistake the retired job's receipt for this live allocation.
    assert journal.exists() or not receipt.exists(), 'Archive prior-job transfer receipt before a fresh allocation'
    allocation = rt.json.loads(Path('artifacts/download_allocations.json').read_text())
    ceiling = allocation['new_allocations'][name]
    assert ceiling == 3_000_000_000
    budget = TransferBudget(journal, ceiling)
    store = None
    checkpoint = None
    start = time.time()
    try:
        inputs = dict(selection=rt.sha(selection), genes=rt.sha(gene_path), offsets=rt.sha(offsets_path),
                      version=version, implementation=rt.sha(__file__),
                      bridge=rt.sha(Path(__file__).with_name('durable_public_store.py')))
        bind_inputs(budget, inputs, checkpoint_root)
        store = PublicStore(budget)
        identity = dict(inputs=inputs, repository_revision=store.revision)
        rt.guard(int(ptr[-1]) * 8 + 1024**3)
        checkpoint = DownloadCheckpoint(checkpoint_root,
            dict(data=((int(ptr[-1]),), np.int32), indices=((int(ptr[-1]),), np.int32)), identity)
        resume = checkpoint.next_plan
        assert 0 <= resume <= len(sel) and (resume % 128 == 0 or resume == len(sel))
        # A completely materialized checkpoint can be published without reconnecting Lance.
        if resume < len(sel):
            expression = lance.dataset(f's3://public/data/{name}/expression.lance', version=version,
                                       storage_options=store.options(), block_size=4096)
            assert expression.version == version
            for lo in range(resume, len(sel), 128):
                rt.guard(1024**3)
                rows = sel.iloc[lo:lo + 128]
                length = rows.n_genes_by_counts.to_numpy(np.int64)
                ids = rows.cell_integer_id.to_numpy(np.int64)
                group = ids // 5000
                assert np.all(off['logical_start'][group] >= 0)
                starts = rows.cell_start_index.to_numpy(np.int64) - off['original_start'][group] + off['logical_start'][group]
                ix = np.concatenate([np.arange(s, s + n, dtype=np.uint64) for s, n in zip(starts, length)])
                tab = expression.take(ix, columns=['cell_integer_id', 'gene_integer_id', 'value'])
                np.testing.assert_array_equal(tab['cell_integer_id'].to_numpy(), np.repeat(ids, length))
                values = tab['value'].to_numpy()
                assert np.isfinite(values).all() and np.equal(values, np.floor(values)).all()
                assert (values >= 0).all() and (values <= np.iinfo(np.int32).max).all()
                local = np.repeat(np.arange(len(rows)), length)
                np.testing.assert_array_equal(np.bincount(local, weights=values, minlength=len(rows)),
                                              rows.total_counts.to_numpy())
                gene_ids = tab['gene_integer_id'].to_numpy()
                assert np.all((gene_ids >= 0) & (gene_ids < len(mapping)))
                stop = min(lo + 128, len(sel))
                segment = slice(int(ptr[lo]), int(ptr[stop]))
                checkpoint.arrays['data'][segment] = values.astype(np.int32)
                checkpoint.arrays['indices'][segment] = mapping[gene_ids]
                checkpoint.commit(stop)
                rt.save(receipt, dict(bytes=budget.total, cells_done=stop, version=version,
                    repository_revision=store.revision, status='in progress', charge_journal=str(journal)))
                if lo % 1280 == 0 or stop == len(sel):
                    print(name, stop, len(sel), 'GB', round(budget.total / 1e9, 3),
                          'sec', round(time.time() - start), 'resume', resume, flush=True)
        revision = store.revision
        store.close()
        store = None
        # CSR canonicalization mutates its buffers; keep the durable raw payload
        # untouched so failed publication can safely be retried with original ptr.
        x = sparse.csr_matrix((np.array(checkpoint.arrays['data'], copy=True),
                               np.array(checkpoint.arrays['indices'], copy=True), ptr),
                              shape=(len(sel), len(unique)))
        x.sum_duplicates()
        x.sort_indices()
        np.testing.assert_array_equal(np.asarray(x.sum(axis=1)).ravel(), sel.total_counts.to_numpy())
        obs = sel[['target_gene', 'batch', 'evaluation_role', 'source_row', 'cell_integer_id']].copy()
        obs.index = [f'{name}_{i}' for i in sel.cell_integer_id]
        a = ad.AnnData(x, obs=obs, var=pd.DataFrame(index=unique))
        a.uns['evaluation_only'] = True
        rt.guard(2 * 1024**3)
        temporary = path.with_suffix('.partial.h5ad')
        a.write_h5ad(temporary, compression='gzip', compression_opts=1)
        temporary.replace(path)
        rt.save(path.with_suffix('.json'), dict(sha256=rt.sha(path), version=version,
            source='https://huggingface.co/datasets/slaf-project/X-Atlas-Orion',
            repository_revision=revision, selection_sha256=rt.sha(selection), bytes=budget.total,
            charge_journal=str(journal)))
        rt.save(receipt, dict(bytes=budget.total, cells_done=len(sel), version=version,
            repository_revision=revision, status='complete', charge_journal=str(journal)))
        checkpoint.remove()
        print('DONE', name, a.shape, flush=True)
    finally:
        if store is not None:
            store.close()
        budget.close()


if __name__ == '__main__':
    main(sys.argv[1])
