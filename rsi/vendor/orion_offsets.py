"""Restore publisher cell offsets across out-of-order Lance ingestion chunks.
Every 5,000-cell ingestion block is independently checked against cells.lance.
"""
import runtime
import numpy as np

def offset_mapping(name,e,bridge):
 import lance
 path=runtime.ROOT/f'artifacts/{name}_offset_mapping.npz'
 if path.exists():
  d=np.load(path);assert int(d['version'])==e.version;return d['original_start'],d['logical_start']
 fragments=e.get_fragments();length=np.array([f.metadata.physical_rows for f in fragments],np.int64);prefix=np.r_[0,np.cumsum(length)];ids=np.array([f.fragment_id for f in fragments],np.uint64)
 # Physical row addresses are stable even when source ingestion chunks were appended out of order.
 boundary=runtime.ROOT/f'artifacts/{name}_fragment_boundaries.npz'
 if boundary.exists():
  d=np.load(boundary);assert np.array_equal(ids,d['fragment_ids']);first=d['first'];last=d['last']
 else:
  first=e._take_rows((ids<<np.uint64(32)).tolist(),columns=['cell_integer_id'])['cell_integer_id'].to_numpy().astype(np.int64)
  last=e._take_rows(((ids<<np.uint64(32))+(length.astype(np.uint64)-1)).tolist(),columns=['cell_integer_id'])['cell_integer_id'].to_numpy().astype(np.int64)
  np.savez(boundary,first=first,last=last,fragment_ids=ids,physical_rows=length,version=e.version)

 groups=first//5000;assert np.all(last//5000==groups),'A fragment crosses an ingestion block'
 ng=int(groups.max())+1;logical=np.full(ng,-1,np.int64);block_length=np.zeros(ng,np.int64)
 for g in range(ng):
  ii=np.flatnonzero(groups==g)
  if not len(ii):continue
  assert np.array_equal(ii,np.arange(ii[0],ii[-1]+1));assert first[ii[0]]==g*5000;assert np.all(np.diff(first[ii])>=0);logical[g]=prefix[ii[0]];block_length[g]=length[ii].sum()
 cell_meta=lance.dataset(f'hf://datasets/slaf-project/X-Atlas-Orion/data/{name}/cells.lance',storage_options={'hf_token':''})
 cells=lance.dataset(f's3://public/data/{name}/cells.lance',version=cell_meta.version,storage_options=bridge.options())
 starts=cells.take(np.arange(ng,dtype=np.uint64)*5000,columns=['cell_integer_id','cell_start_index']);assert np.array_equal(starts['cell_integer_id'].to_numpy(),np.arange(ng)*5000);original=starts['cell_start_index'].to_numpy().astype(np.int64)
 final=cells.take([cells.count_rows()-1],columns=['cell_start_index','n_genes_by_counts']);total=int(final['cell_start_index'][0].as_py())+int(final['n_genes_by_counts'][0].as_py());expected=np.diff(np.r_[original,total]);present=logical>=0;assert np.array_equal(expected[present],block_length[present]),'Block row counts disagree with source cell metadata';assert int(block_length.sum())==e.count_rows();coverage=e.count_rows()/total;assert coverage>=.99,'Source expression coverage below released minimum';missing=np.flatnonzero(~present)
 np.savez(path,original_start=original,logical_start=logical,block_length=block_length,first_cell=first,last_cell=last,version=e.version)
 runtime.save(runtime.ROOT/f'artifacts/{name}_offset_verification.json',dict(version=e.version,fragments=len(fragments),blocks=ng,total_rows=total,moved_blocks=int(np.sum(original!=logical)),max_row_displacement=int(np.max(abs(original[present]-logical[present]))),missing_blocks=missing.tolist(),missing_expression_records=int(expected[missing].sum()),expression_coverage=coverage,read_bytes=bridge.bytes,verification='Exact physical block lengths, contiguous fragment membership, source cell IDs and original cumulative offsets; every fetched cell subsequently checked again'))
 print(name,'offset map verified; moved blocks',int(np.sum(original[present]!=logical[present])),'missing blocks',missing.tolist(),'coverage',coverage,flush=True)
 return original,logical
