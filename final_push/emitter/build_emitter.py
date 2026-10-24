"""Panel- and context-agnostic dual-moment build with independently switchable emitter knobs.

Inherits v7's numerics unchanged (model.py heads, v5.py budgets, feasible_bulk.project_bulk, AtlasShift
dual_moment_counts, verify_counts.moment_certificate) and adds, each switchable (emitter.Knobs):
  per-target amplitude allocation, control-CPM silence gate (with a floor), expression-compatibility mask,
  claim at the control-400 cap with a HARD B/claim gate, independent CPM-head amplitude, call budget and
  p-value ordering hooks.

Usage (never submits):
  python build_emitter.py --name NAME --targets perts.csv --genes genes.csv \
      --contexts "D=/abs/context_D.h5ad,E=...,F=..." --banks DIR [--knobs knobs.json] [--fixture N]

--fixture N emits only N probe groups per context and writes receipts; a full run writes out/NAME.h5ad.
"""
import runtime                                   # v7 runtime: one worker, bounded disk, local TMPDIR
import argparse, json, gc, time, sys
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
from scipy import sparse
import model, v5
from feasible_bulk import project_bulk
from verify_counts import moment_certificate, audit_h5
from emitter import Knobs, allocation, silence_factors, compat_mask, ratio_gate

rt = runtime.ROOT


class Writer:
    def __init__(self, path, genes):
        import h5py
        self.n = 0; self.nnz = 0; self.obs = []; self.genes = genes; self.f = h5py.File(path, 'w')
        self.f.attrs.update({'encoding-type': 'anndata', 'encoding-version': '0.1.0'})
        self.g = self.f.create_group('X'); self.g.attrs.update({'encoding-type': 'csr_matrix', 'encoding-version': '0.1.0'})
        for k, dt in [('data', 'int32'), ('indices', 'int32'), ('indptr', 'int64')]:
            self.g.create_dataset(k, shape=(1 if k == 'indptr' else 0,), maxshape=(None,), dtype=dt,
                                  chunks=(4096 if k == 'indptr' else 262144,), compression='gzip', compression_opts=1, shuffle=True)
        self.g['indptr'][0] = 0

    def append(self, x, c, t):
        runtime.guard(); x = sparse.csr_matrix(x)
        for k, a in [('data', x.data.astype(np.int32)), ('indices', x.indices.astype(np.int32)),
                     ('indptr', x.indptr[1:].astype(np.int64) + self.nnz)]:
            z = self.g[k]; j = len(z); z.resize((j + len(a),)); z[j:] = a
        self.n += x.shape[0]; self.nnz += x.nnz; self.obs.extend([(c, t)] * x.shape[0])

    def close(self):
        self.g.attrs['shape'] = (self.n, len(self.genes))
        o = pd.DataFrame(self.obs, columns=['context', 'target_gene'], index=[f'pred_{i}' for i in range(self.n)])
        for k in o:
            o[k] = pd.Categorical(o[k])
        ad.io.write_elem(self.f, 'obs', o); ad.io.write_elem(self.f, 'var', pd.DataFrame(index=pd.Index(self.genes)))
        for k in ['uns', 'obsm', 'varm', 'obsp', 'varp', 'layers']:
            ad.io.write_elem(self.f, k, {})
        self.f.close()


class Context:
    """One destination context: heads, knobs, claim calibration, emission."""

    def __init__(self, a, targets, genes, label, ci, knobs, log=print):
        self.targets = np.asarray(targets); self.genes = np.asarray(genes); self.label = label; self.ci = ci
        self.k = knobs; self.log = log
        self.q4, self.d, self.mean, self.bulk = model.template(a, genes, pool_k=4, seed=20260910 + ci)
        # CPM head at its own amplitude, bulk head on the v5 direction (both from the same source fusion)
        (self.p, b), _ = model.heads(targets, genes, self.mean, self.bulk, 'anchor_single_mean', de_amplitude=knobs.cpm_amplitude)
        V, s = v5.control_axes(a.X, k=50, seed=ci)
        self.Lc = v5.logbulk(self.bulk)[None, :]
        delta = v5.deflate(v5.logbulk(b) - self.Lc, V)
        info = dict(context=label, control_cells=int(a.n_obs), targets=len(targets), knobs=knobs.asdict(),
                    bulk_deflation_k=50, cpm_amplitude=knobs.cpm_amplitude)
        # --- knob: expression compatibility (destination controls only)
        delta, cinfo = compat_mask(delta, self.bulk, knobs); info['compat'] = cinfo
        # --- knob: cross-source confidence, used both for allocation and (later) for call ordering
        self.confidence, ainfo = self._confidence(delta)
        info['confidence'] = ainfo
        # --- knob: per-target amplitude allocation
        f_alloc, alloc_info = allocation(delta, self.confidence, knobs); info['allocation'] = alloc_info
        # --- knob: per-target silence gate (destination controls only)
        f_sil, sinfo = silence_factors(targets, genes, self.bulk, knobs); info['silence'] = sinfo
        self.row_factor = f_alloc * f_sil
        self.L = self.Lc + delta * self.row_factor[:, None]
        # --- knob: per-cell scatter. The 4-cell pooled template under-disperses, so our own jackknife claim is
        # about half the real-cell cap; less pooling raises the claim toward it WITHOUT moving the CPM mean head
        # (the two heads are fitted independently by dual_moment_counts). 'selective' is v7b's rule: real
        # single-cell scatter on genes with no cross-source support, 4-cell pooling on supported genes, so the
        # genes that carry the DE calls keep their low within-group variance.
        self.q1 = None; self.support = None
        if knobs.scatter != 'pool4':
            self.q1, _, _, _ = model.template(a, genes, pool_k=1, seed=20260910 + ci)
            if knobs.scatter == 'selective':
                self.support = self.gene_confidence >= 1e6          # net cross-source agreement >= 1
                info['scatter'] = dict(mode='selective', supported_genes_mean=float(self.support.sum(1).mean()))
            elif knobs.scatter == 'pool1':
                info['scatter'] = dict(mode='pool1')
            else:
                raise ValueError(knobs.scatter)
        else:
            info['scatter'] = dict(mode='pool4')
        self.w = v5._row_weights(len(targets), len(genes), v5.excl_cols(targets, genes))
        # control-400 jackknife cap at this context's own depth (label-free; the server's treated cap is unknown)
        rng = np.random.default_rng(73591 + ci); X = sparse.csr_matrix(a.X)
        self.ctrl_jk = np.array([self._jk(X[rng.choice(a.n_obs, 400, replace=False)]) for _ in range(8)])
        self.cap_jk = float(self.ctrl_jk.mean())
        info['control400_jk_samples'] = self.ctrl_jk.tolist(); info['control400_claim'] = float(self.cap_jk * self.w.sum())
        self.probe = np.sort(np.random.default_rng(ci).choice(len(targets), min(8, len(targets)), replace=False))
        self.alpha = 1.0; self.info = info

    def _confidence(self, delta):
        """Net cross-source sign agreement (v5b's key, variance-weighted), per target: mean over its called genes."""
        names = model.load_banks(None, True)
        src = [model.aligned_bank(n, self.targets, self.genes, 'log2fc') for n in names]
        cd = model.cd4_family(self.targets, self.genes, 'log2fc')
        if cd is not None:
            src.append(cd)
        if not src:
            return np.full(len(self.targets), 0.5), dict(sources=[], mode='none')
        sgn = np.sign(delta); net = np.zeros(delta.shape, np.float32); zsum = np.zeros(delta.shape)
        for e, m, v in src:
            z = np.where(m & np.isfinite(v), e / np.sqrt(np.maximum(v, 1e-12)), 0.0); sig = np.abs(z) > 1.96
            net += (sig & (np.sign(z) == sgn) & (sgn != 0)).astype(np.float32)
            net -= (sig & (np.sign(z) == -sgn) & (sgn != 0)).astype(np.float32)
            zsum += np.minimum(np.abs(z), 100) * (np.sign(z) == sgn)
        self.gene_confidence = net * 1e6 + zsum                     # v9's call-ordering key
        support = (net >= 1).sum(1).astype(np.float64)
        conf = support / max(support.max(), 1.0)
        del src, cd; gc.collect()
        return conf, dict(sources=names, mode='net_source_agreement', supported_genes_mean=float(support.mean()))

    def _jk(self, x):
        return float(v5.jackknife_correction(sparse.csr_matrix(x), np.zeros(x.shape[0], np.intp), 1, v5.TS, chunk=32)[0])

    def projected(self, alpha):
        raw = v5.from_logbulk(v5.rescale(self.L, alpha, self.Lc)); z = self.d / self.d.mean()
        return project_bulk(raw, self.p, float(z.min()), float(z.max()))

    def solve_alpha(self, target):
        lo, hi = 0.0, 8.0
        B = lambda al: v5.budget(self.projected(al)[0], self.targets, self.genes)
        if B(hi) < target:
            return hi, B(hi)
        assert B(lo) <= target, 'feasible DE moments already exceed the requested bulk budget'
        for _ in range(38):
            mid = .5 * (lo + hi)
            if B(mid) < target:
                lo = mid
            else:
                hi = mid
        al = .5 * (lo + hi); return al, B(al)

    def template_for(self, i):
        if self.q1 is None:
            return self.q4
        if self.support is None:
            return self.q1
        q = self.q1.copy(); m = self.support[i]; q[:, m] = self.q4[:, m]
        q /= q.sum(1, keepdims=True); return q

    def emit(self, i, b):
        runtime.guard()
        x = model.atlas.dual_moment_counts(self.template_for(i), self.p[i], b, depths=self.d,
                                           seed=model.atlas.seed_for(self.label + ':' + self.targets[i], 20260910), iterations=400)
        cert = moment_certificate(x, self.d, self.p[i], b)
        if cert['mean_L1_upper'] > .015 or cert['bulk_L1_upper'] > .0015:
            raise RuntimeError(f'dual moment violation {self.label}/{self.targets[i]} {cert}')
        return x, cert

    def calibrate(self):
        """Iterate alpha so B = claim_rho x capped claim, then HARD-gate the realised ratio."""
        hist = []
        for it in range(5):
            bf, _ = self.projected(self.alpha)
            jks = np.array([self._jk(self.emit(int(i), bf[int(i)])[0]) for i in self.probe])
            if self.k.claim_source == 'probe':
                claim = float(jks.mean() * self.w.sum())
            elif self.k.claim_source == 'control400':
                claim = float(self.cap_jk * self.w.sum())
            else:
                claim = float(np.mean(np.minimum(jks, self.cap_jk)) * self.w.sum())
            alpha, B = self.solve_alpha(self.k.claim_rho * claim)
            hist.append(dict(iteration=it, alpha_measured=self.alpha, claim=claim, alpha_next=alpha, B_next=B,
                             probe_claim=float(jks.mean() * self.w.sum()), control400_claim=float(self.cap_jk * self.w.sum())))
            self.log(self.label, 'calibration', hist[-1], flush=True)
            prev = self.alpha; self.alpha = alpha
            if abs(alpha / max(prev, 1e-12) - 1) < .005:
                break
        self.bf, proj = self.projected(self.alpha)
        assert proj.max() <= .03, 'source bulk exceeds the public feasible-projection limit'
        B = v5.budget(self.bf, self.targets, self.genes)
        self.claim = hist[-1]['claim']
        self.gate = ratio_gate(B, self.claim, self.k, label=f'{self.label}:requested')
        self.info.update(alpha=self.alpha, calibration=hist, B_requested=B, claim_capped=self.claim,
                         ratio_gate_requested=self.gate, feasible_projection_L1_max=float(proj.max()))


def main(a):
    knobs = Knobs(**json.loads(Path(a.knobs).read_text())) if a.knobs else Knobs()
    # v7's model.py resolves its default bank dir relative to its own parent; this script lives one level deeper.
    model.BANKS = Path(a.banks) if a.banks else (rt.parents[1] / 'public_data/model_v4_signature_banks')
    assert (model.BANKS / 'K562_GWPS_raw.npz').exists() or (model.BANKS / 'K562_GWPS_bulk.npz').exists(), f'no source banks in {model.BANKS}'
    assert (rt / 'artifacts/promoter_pairs.csv').exists(), 'promoter_pairs.csv must sit in artifacts/'
    genes = pd.read_csv(a.genes).iloc[:, 0].astype(str).to_numpy()
    if genes[0] in ('gene_name', 'gene'):
        genes = genes[1:]
    tf = pd.read_csv(a.targets); col = 'target_gene' if 'target_gene' in tf else tf.columns[0]
    targets = tf.loc[tf[col].astype(str) != 'non-targeting', col].astype(str).to_numpy()
    contexts = dict(s.split('=', 1) for s in a.contexts.split(','))
    start = time.time(); recs = []; moments = []
    out = rt / 'out' / f'{a.name}.h5ad'
    writer = None if a.fixture else Writer(out, genes)
    try:
        for ci, (c, fp) in enumerate(contexts.items()):
            runtime.guard(); adata = ad.read_h5ad(fp); assert np.array_equal(adata.var_names, genes)
            ctx = Context(adata, targets, genes, c, ci, knobs); del adata; gc.collect()
            ctx.calibrate()
            ids = ctx.probe[:a.fixture] if a.fixture else range(len(targets))
            P = []; jks = []
            for i in ids:
                i = int(i); x, cert = ctx.emit(i, ctx.bf[i]); jks.append(ctx._jk(x)); P.append(x.sum(0, dtype=np.int64))
                moments.append(dict(context=c, target=targets[i], **cert))
                if writer is not None:
                    writer.append(x, c, targets[i])
                if len(jks) % 20 == 1:
                    print(a.name, c, i, 'elapsed_s', round(time.time() - start), flush=True)
            if not a.fixture:
                Pn = np.array(P) / np.array(P).sum(1, keepdims=True)
                Be = v5.budget(Pn, targets, genes)
                claim_emitted = float(ctx.w @ np.minimum(np.array(jks), ctx.cap_jk))
                ctx.info.update(B_emitted=Be, claim_pred_emitted=float(ctx.w @ np.array(jks)), claim_capped_emitted=claim_emitted,
                                ratio_gate_emitted=ratio_gate(Be, claim_emitted, knobs, label=f'{c}:emitted'))
            recs.append(ctx.info)
            print('CONTEXT DONE', c, json.dumps({k: v for k, v in ctx.info.items() if k != 'calibration'}, default=str), flush=True)
            del ctx; gc.collect()
    finally:
        if writer is not None:
            writer.close()
    rec = dict(status='prediction-side checks only; official score UNVERIFIED', name=a.name, knobs=knobs.asdict(),
               contexts=recs, moments=moments, seconds=round(time.time() - start, 1),
               source_sha256={n: runtime.sha(rt / n) for n in ['build_emitter.py', 'emitter.py', 'model.py', 'v5.py', 'feasible_bulk.py', 'vendor/atlas_model.py']})
    runtime.save(rt / 'artifacts' / f'{a.name}_{"fixture" if a.fixture else "build"}.json', rec)
    if writer is not None:
        cert = audit_h5(out, genes, targets, contexts=tuple(contexts))
        runtime.save(rt / 'artifacts' / f'{a.name}_structure.json', cert)
        print('BUILD VERIFIED', a.name, cert['shape'], cert['nnz'], flush=True)
    print('DONE', a.name, round(time.time() - start, 1), 's', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--name', required=True); p.add_argument('--targets', required=True); p.add_argument('--genes', required=True)
    p.add_argument('--contexts', required=True); p.add_argument('--banks', default=None); p.add_argument('--knobs', default=None)
    p.add_argument('--fixture', type=int, default=0)
    main(p.parse_args())
