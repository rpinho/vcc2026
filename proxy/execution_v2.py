"""Operational scheduling and exact-result memoization; no metric math changes."""
import collections
import ast
import hashlib
import importlib
import json
import os
from pathlib import Path
import subprocess
import time
import weakref
import numpy as np

_activity=[0.,0.]


def kids_cpu():
    return 0.0


def install_guard(runtime):
    original=runtime.guard
    def guard(extra=0,floor=runtime.FLOOR_GIB,yield_kids=True):
        # Honor the existing Mini policy: a running family release/draft pauses
        # science even when it is temporarily idle. CPU use is not permission.
        original(extra,floor=floor,yield_kids=yield_kids)
    runtime.guard=guard


def numeric_adapter_sha256():
    names={'array_hash','install_sorted_mwu','install_mwu_cache','install_stat_cache'}
    tree=ast.parse(Path(__file__).read_text())
    parts=[ast.dump(n,include_attributes=False) for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
    assert len(parts)==4
    return hashlib.sha256('\n'.join(parts).encode()).hexdigest()


def array_hash(values):
    h=hashlib.sha256()
    for value in values:
        a=np.ascontiguousarray(value);h.update(str(a.dtype).encode());h.update(str(a.shape).encode());h.update(memoryview(a).cast('B'))
    return h.digest()


def install_sorted_mwu():
    """Keep the installed official MWU; the unlicensed optimization is excluded."""
    return dict(status='official installed backend; sorted-reference optimization excluded')


def install_mwu_cache(guard):
    """Cache bit-identical MWU inputs (notably all 400-cell tiled baseline groups)."""
    pdex=importlib.import_module('pdex');original=pdex.mwu
    references={};cache=collections.OrderedDict();counters=dict(hits=0,misses=0)
    def memo(x,y):
        guard()
        # pdex's reference path always supplies SparseColumnIndex as y.
        if not hasattr(y,'col_indptr') or not hasattr(x,'indptr'):
            return original(x,y)
        key=id(y.data);record=references.get(key)
        if record is None or record[0]() is not y.data:
            digest=array_hash((y.data,y.col_indptr,y.col_order,np.array([y.n_rows,y.n_cols])))
            references.clear();references[key]=(weakref.ref(y.data),digest)
        else:digest=record[1]
        group=array_hash((x.data,x.indices,x.indptr,np.array(x.shape)))
        token=(group,digest)
        if token in cache:
            counters['hits']+=1;cache.move_to_end(token);return cache[token]
        counters['misses']+=1;value=original(x,y);cache[token]=value
        while len(cache)>8:cache.popitem(last=False)
        return value
    pdex.mwu=memo
    return counters


def install_stat_cache(driver):
    """Exact sufficient-statistic memoization, including identical tiled baseline groups."""
    original=driver.stat;cache=collections.OrderedDict();counts=dict(hits=0,misses=0)
    def stat(a,cfg):
        if a.obs[cfg.pert_col].astype(str).nunique()!=1:return original(a,cfg)
        x=a.X.tocsr();sizes=np.diff(x.indptr)
        token=None
        if x.dtype.kind=='f' and len(sizes) and sizes[0]>0 and np.all(sizes==sizes[0]):
            indices=x.indices.reshape(x.shape[0],int(sizes[0]));data=x.data.reshape(x.shape[0],int(sizes[0]))
            if np.all(indices==indices[0]) and np.all(data==data[0]):
                token=('tile',x.shape,cfg.bulk_target_sum,array_hash((data[0],indices[0])))
        if token is None:token=('csr',x.shape,cfg.bulk_target_sum,array_hash((x.data,x.indices,x.indptr)))
        if token in cache:
            counts['hits']+=1;cache.move_to_end(token);return cache[token]
        counts['misses']+=1;value=original(a,cfg);cache[token]=value
        while len(cache)>800:cache.popitem(last=False)
        return value
    driver.stat=stat
    return counts
