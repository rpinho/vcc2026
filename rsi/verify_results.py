"""Independent arithmetic/protocol audit; never reimplements biological metrics."""
import runtime as rt
from pathlib import Path
from fractions import Fraction as F
import hashlib,time,itertools
import numpy as np,pandas as pd
from scipy.stats import spearmanr,pearsonr
SIX=['pds_cosine','expr_mse_unbiased_capped_norm','de_wilcoxon_lfc_nmae','de_wilcoxon_direction_fidelity_yield_raw','de_wilcoxon_direction_reach_raw','de_wilcoxon_sig_jaccard'];SHORT=dict(zip(SIX,['pds','mse','nmae','fid','reach','jac']))

def check(record,anchors):
 out={}
 for metric in SIX:
  value,base,rep=map(lambda x:F(float(x)),[record['raw'][metric],anchors['baseline'][metric],anchors['replicate'][metric]]);s=(value-base)/(rep-base)
  if metric=='expr_mse_unbiased_capped_norm':s=min(F(1),max(F(0),s))
  if metric=='de_wilcoxon_lfc_nmae':s=max(F(-6),s)
  assert abs(float(s)-record['scaled'][metric])<1e-12,(metric,s,record['scaled'][metric]);out[metric]=s
 overall=sum(out.values(),F(0))/6;assert abs(float(overall)-record['overall'])<1e-12
 return dict(exact_overall_on_binary64_inputs=str(overall),max_scaling_error=max(abs(float(out[k])-record['scaled'][k]) for k in out),status='PASS arithmetic; biological validity unverified')

def official(text=None):
 text=(rt.ROOT.parent/'OFFICIAL_SCORES.md').read_text() if text is None else text;rows={}
 for line in text.splitlines():
  if not line.startswith('|'):continue
  cols=[c.strip().replace('**','').replace('−','-') for c in line.strip('|').split('|')]
  for key in ['v1','v2b','v3b','v4a','v4b']:
   if cols and (cols[0]==key or cols[0].startswith(key+' ')):
    try:rows[key]=float(cols[1])
    except (ValueError,IndexError):pass
 return rows

def correlations(local,numbers):
 keys=sorted(set(local)&set(numbers));xx=np.array([local[k] for k in keys]);yy=np.array([numbers[k] for k in keys]);n=len(keys)
 if n<3 or np.std(xx)==0:return None
 rho=float(spearmanr(xx,yy).statistic);pear=pearsonr(xx,yy);order=all((local[a]-local[b])*(numbers[a]-numbers[b])>0 for i,a in enumerate(keys) for b in keys[i+1:]);permutations=list(itertools.permutations(yy)) if n<=8 else []
 p=sum(abs(float(spearmanr(xx,y).statistic))>=abs(rho)-1e-12 for y in permutations)/len(permutations) if permutations else None
 return dict(n=n,recipes=keys,spearman=rho,pearson=float(pear.statistic),pearson_p=float(pear.pvalue),spearman_exact_two_sided_permutation_p=p,all_pairs_rank_consistent=order,local={k:float(local[k]) for k in keys},official={k:numbers[k] for k in keys},v1_v2_order_matches=bool((local['v1']-local['v2b'])*(numbers['v1']-numbers['v2b'])>0) if {'v1','v2b'}<=set(keys) else None)

def calibrate():
 official_bytes=(rt.ROOT.parent/'OFFICIAL_SCORES.md').read_bytes();numbers=official(official_bytes.decode());records=[];audit={};protocols={}
 for folder in Path('artifacts/scores').glob('*'):
  if not(folder/'anchors.json').exists():continue
  anchors=rt.json.loads((folder/'anchors.json').read_text())
  protocols[folder.name]=anchors.get('baseline_construction','UNVERIFIED')
  for path in folder.glob('*_s[0-9]*.json'):
   r=rt.json.loads(path.read_text());audit[str(path)]=check(r,anchors);records.append(r)
 rt.save('artifacts/arithmetic_audit.json',audit)
 if not records:return None
 df=pd.DataFrame([dict(task=r['task'],recipe=r['recipe'],seed=r['seed'],overall=r['overall'],exact_control_match=r['exact_target_matched'],**{SHORT[k]:v for k,v in r['scaled'].items()}) for r in records]);df.to_csv('artifacts/calibration_scores.csv',index=False);corr={}
 canonical=df[df.seed==42]
 for task in [t for t in canonical.task.unique() if 'holdout' not in t]:
  sub=canonical[(canonical.task==task)&canonical.recipe.isin(numbers)].groupby('recipe').overall.mean();result=correlations(sub.to_dict(),numbers)
  if result:corr[task]=result
 # Common native SEARCH tasks only; a depth sensitivity never counts as another lineage.
 roles=rt.json.loads(Path('artifacts/benchmark_roles.json').read_text());primary=roles['primary_search'];native=canonical[canonical.task.isin(primary)&canonical.recipe.isin(numbers)];complete=[task for task,g in native.groupby('task') if set(numbers)<=set(g.recipe)];macro=None
 if complete:
  means=native[native.task.isin(complete)].groupby('recipe').overall.mean().to_dict();macro=correlations(means,numbers)
  if macro:macro.update(tasks=complete,n_contexts=len(complete))
 rank_gate=bool(macro and macro['n_contexts']>=2 and macro['all_pairs_rank_consistent']);manifest=rt.json.loads(Path('artifacts/tasks/manifest.json').read_text());matched=bool(complete and all(manifest[t]['exact_matched'] for t in complete));documented_baseline=bool(complete and all(isinstance(protocols.get(t),dict) for t in complete));gate=rank_gate and matched and documented_baseline;rt.save('artifacts/calibration_gate.json',dict(status='PRELIMINARY RANK CONSISTENT' if gate else 'CLOSED',phase2_allowed_by_observed_rank=rank_gate,phase2_allowed=gate,matched_control_requirement=matched,documented_baseline_requirement=documented_baseline,reason='Need all-pair rank consistency on at least two common native SEARCH contexts with target-matched controls and Arc-documented mean-response baseline. Small n remains; HOLDOUT is mandatory before promotion. A pooled-control rank correlation alone does not satisfy the benchmark design.',macro=macro))
 out=dict(official=numbers,official_file_sha256=hashlib.sha256(official_bytes).hexdigest(),read_utc=time.strftime('%Y-%m-%d %H:%M:%S UTC',time.gmtime()),baseline_protocols=protocols,per_task=corr,common_native_macro=macro,note='Small, dependent calibration set. External AtlasShift/Sidechain labels excluded without exact reproduction. All HOLDOUT tasks excluded. Perfect Spearman with n=3 has exact two-sided permutation p=1/3; it does not validate a precise score forecast. Control-mean scale is sensitivity only: Arc primary documentation specifies mean perturbation response; its hidden construct-filter identity remains unverified.')
 rt.save('artifacts/calibration.json',out);return out
if __name__=='__main__':print(calibrate())
