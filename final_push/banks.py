"""Source-bank access for the local selector (billm format, canonical 18,533-gene axis; read-only).

model_billm/data/*.npz: response = log1p(5e4*p_pert) - log1p(5e4*p_ctrl_matched), control = matched control probability,
measured mask, targets, n_cells, family. model_billmB/data/{K562_GWPS_full,H1_full}_*.npz shards (Codex job owns that dir;
we only read). Label-free: public screens only.
"""
import hashlib
from pathlib import Path
import numpy as np

VCC = Path(__file__).resolve().parents[1]
BILLM = VCC / 'model_billm' / 'data'
BILLMB = VCC / 'model_billmB' / 'data'
CACHE = VCC / 'final_push' / 'artifacts' / 'cache'
CACHE.mkdir(parents=True, exist_ok=True)
TS = 50000.0


def family_of(name):
    for pre in ('CD4', 'K562', 'H1'):
        if name.startswith(pre):
            return pre
    return name


def source_of(name):
    if name.startswith('K562_GWPS'):
        return 'K562_GWPS'
    if name.startswith('H1'):
        return 'H1'
    return name


class Bank:
    """One source screen; rows gathered lazily by target from one or more npz shards."""

    def __init__(self, source, paths):
        self.source = source; self.family = family_of(source); self.paths = [Path(p) for p in paths]
        self.index = {}; self.n_cells = {}; self.genes = None; self.measured = None
        for k, p in enumerate(self.paths):
            with np.load(p, allow_pickle=False) as z:
                t = z['targets'].astype(str); n = z['n_cells']; m = z['measured'].astype(bool)
                if self.genes is None:
                    self.genes = z['genes'].astype(str); self.measured = m
                else:
                    assert np.array_equal(self.genes, z['genes'].astype(str))
                    self.measured = self.measured & m
                for i, tt in enumerate(t):
                    if tt not in self.index:   # first occurrence wins: batch-matched raw rows precede global-control shards
                        self.index[tt] = (k, i); self.n_cells[tt] = int(n[i])
        self.targets = np.array(sorted(self.index)); self._mean = None

    def __repr__(self):
        return f'Bank({self.source}, {len(self.targets)} targets, {len(self.paths)} shards)'

    def has(self, targets):
        return np.array([t in self.index for t in targets])

    def rows(self, targets, field='response'):
        out = np.full((len(targets), len(self.genes)), np.nan, np.float32); want = {}
        for j, t in enumerate(targets):
            if t in self.index:
                k, i = self.index[t]; want.setdefault(k, []).append((i, j))
        for k, pairs in want.items():
            with np.load(self.paths[k], allow_pickle=False) as z:
                arr = z[field]; ii = np.array([i for i, _ in pairs]); jj = np.array([j for _, j in pairs])
                o = np.argsort(ii); out[jj[o]] = arr[ii[o]]
        return out

    def source_mean(self):
        """Source-wide centre: mean response over ALL rows (own-target entries left out); cached on disk."""
        if self._mean is not None:
            return self._mean
        key = hashlib.sha256('|'.join(f'{p}:{p.stat().st_size}' for p in self.paths).encode()).hexdigest()[:16]
        f = CACHE / f'mean_{self.source}_{key}.npz'
        if f.exists():
            with np.load(f) as z:
                self._mean = dict(mean=z['mean'], control=z['control'], n=int(z['n']))
            return self._mean
        gi = {g: i for i, g in enumerate(self.genes)}; G = len(self.genes)
        s = np.zeros(G); cnt = np.zeros(G); c = np.zeros(G); n = 0; seen = set()
        for p in self.paths:
            with np.load(p, allow_pickle=False) as z:
                r = z['response']; ctrl = z['control']; t = z['targets'].astype(str)
                for i in range(0, len(t), 256):
                    tb = t[i:i + 256]; keep = np.array([tt not in seen for tt in tb]); seen.update(tb)
                    if not keep.any():
                        continue
                    blk = r[i:i + 256][keep].astype(np.float64); w = np.ones_like(blk)
                    for row, tt in enumerate(tb[keep]):
                        j = gi.get(tt)
                        if j is not None:
                            w[row, j] = 0
                    s += (blk * w).sum(0); cnt += w.sum(0); c += ctrl[i:i + 256][keep].astype(np.float64).sum(0); n += int(keep.sum())
        mean = np.divide(s, cnt, out=np.zeros_like(s), where=cnt > 0) * self.measured
        self._mean = dict(mean=mean.astype(np.float32), control=(c / max(n, 1)).astype(np.float32), n=n)
        np.savez(f, mean=self._mean['mean'], control=self._mean['control'], n=n)
        return self._mean


def discover(include_shards=True):
    groups = {}
    for p in sorted(BILLM.glob('*.npz')):
        groups.setdefault(source_of(p.stem), []).append(p)
    if include_shards and BILLMB.exists():
        for stem, src in [('K562_GWPS_full', 'K562_GWPS'), ('H1_full', 'H1')]:
            shards = sorted(BILLMB.glob(f'{stem}_*.npz'))
            if shards:
                groups.setdefault(src, []).extend(shards)
    return {s: Bank(s, paths) for s, paths in groups.items()}
