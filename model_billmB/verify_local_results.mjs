// Independent exact arithmetic on stored six-metric vectors, not a new metric kernel.
import fs from 'node:fs';import path from 'node:path';import crypto from 'node:crypto';
const P=import.meta.dirname,R=path.dirname(P),J=path.join(R,'model_billmB');
const six=['pds_cosine','expr_mse_unbiased_capped_norm','de_wilcoxon_lfc_nmae','de_wilcoxon_direction_fidelity_yield_raw','de_wilcoxon_direction_reach_raw','de_wilcoxon_sig_jaccard'];
const read=p=>JSON.parse(fs.readFileSync(p,'utf8'));
const hash=p=>crypto.createHash('sha256').update(fs.readFileSync(p)).digest('hex');
function gcd(a,b){a=a<0n?-a:a;while(b){[a,b]=[b,a%b];}return a;}
function rat(n,d=1n){if(d<0n){n=-n;d=-d;}const g=gcd(n,d);return [n/g,d/g];}
function binary(x){
 if(!Number.isFinite(x))throw Error('Nonfinite scalar');
 const b=new ArrayBuffer(8),v=new DataView(b);v.setFloat64(0,x);
 const bits=v.getBigUint64(0),e=Number((bits>>52n)&2047n),f=bits&((1n<<52n)-1n);
 let n=(e?1n<<52n:0n)+f;if(bits>>63n)n=-n;
 const power=e?e-1023-52:-1074;return power>=0?rat(n<<BigInt(power)):rat(n,1n<<BigInt(-power));
}
const add=(a,b)=>rat(a[0]*b[1]+b[0]*a[1],a[1]*b[1]);
const sub=(a,b)=>add(a,[-b[0],b[1]]);
const abs=a=>[a[0]<0n?-a[0]:a[0],a[1]];
const lte=(a,b)=>a[0]*b[1]<=b[0]*a[1];
const mean=xs=>{const s=xs.reduce(add,[0n,1n]);return rat(s[0],s[1]*BigInt(xs.length));};
const text=a=>`${a[0]}/${a[1]}`;
function checkMean(record){
 if(Object.keys(record.scaled).sort().join()!==[...six].sort().join())throw Error('Metric membership mismatch');
 const values=six.map(k=>record.scaled[k]);const m=mean(values.map(binary));
 const error=abs(sub(binary(record.overall),m));
 const limit=binary(4*Number.EPSILON*Math.max(1,...values.map(Math.abs)));
 if(!lte(error,limit))throw Error('Reported mean disagrees with exact rational mean');
 return {exact:m,error};
}
const spec=read(path.join(R,'rsi/artifacts/v5b_loco_spec.json'));
const canonical=fs.readFileSync(path.join(R,'data/gene_names.csv'),'utf8').trim().split('\n').slice(1).map(x=>x.trim());
const records=[];
for(const task of ['H1_rich_search','CD4_depth192_search']){
 const scorePath=path.join(J,`artifacts/eval/${task}/score.json`);if(!fs.existsSync(scorePath))continue;
 const baselinePath=path.join(R,`rsi/artifacts/baseline_scores/v5b_available_loco_r33/${task}/s20260910.json`);
 const a=read(scorePath),b=read(baselinePath),anchorPath=path.join(R,`rsi/artifacts/scores/${task}/anchors.json`),anchor=read(anchorPath);
 if(a.task!==task||b.task!==task||a.reference_sha256!==b.reference_sha256||a.reference_sha256!==anchor.real_sha256||b.anchors_sha256!==hash(anchorPath))throw Error('Reference/anchor mismatch');
 const aa=checkMean(a),bb=checkMean(b),delta=sub(aa.exact,bb.exact);
 const axisPath=path.join(J,`artifacts/eval/${task}/gene_axis.json`),axis=read(axisPath),genes=spec.tasks[task].genes;
 // Linear name search is independent of the evaluator's dictionary mapping.
 const expected=genes.map(g=>canonical.indexOf(g));
 if(expected.some(i=>i<0)||new Set(expected).size!==genes.length||JSON.stringify(expected)!==JSON.stringify(axis.canonical_indices))throw Error('Independent gene mapping mismatch');
 const geneHash=crypto.createHash('sha256').update(genes.join('\n')).digest('hex');
 if(geneHash!==axis.task_genes_sha256)throw Error('Native gene-axis identity mismatch');
 records.push({task,status:'PASS exact stored-vector aggregation/comparison and independent gene mapping',score_sha256:hash(scorePath),baseline_sha256:hash(baselinePath),anchors_sha256:hash(anchorPath),mapping_sha256:hash(axisPath),bill_exact_mean:text(aa.exact),baseline_exact_mean:text(bb.exact),exact_local_delta:text(delta),delta_sign:delta[0]>0n?1:delta[0]<0n?-1:0,bill_reported_mean_error:text(aa.error),baseline_reported_mean_error:text(bb.error),native_genes:genes.length,canonical_genes:canonical.length});
}
let macro=null;const localPath=path.join(J,'artifacts/local_scores.json');
if(records.length===2&&fs.existsSync(localPath)){
 const local=read(localPath);const expected=mean(['H1_rich_search','CD4_depth192_search'].map(t=>binary(local.tasks[t].overall)));
 if(!lte(abs(sub(binary(local.equal_context_macro),expected)),binary(4*Number.EPSILON)))throw Error('Macro mean mismatch');
 macro={status:'PASS exact equal-context aggregate check',exact_mean_of_reported_overalls:text(expected),sha256:hash(localPath)};
}
const out={utc:new Date().toISOString(),status:records.length===2&&macro?'PASS both SEARCH receipts':'PARTIAL completed SEARCH receipts only',scope:'Exact rational arithmetic on stored binary64 scaled metric aggregates and matched references/anchors; independent gene mapping. Does not rederive raw metric kernels, certify hidden performance, or verify checkpoints.',records,macro};
fs.writeFileSync(path.join(P,'local_verification.json'),JSON.stringify(out,null,2)+'\n');
console.log(out.status,records.map(r=>`${r.task}: exact local delta sign ${r.delta_sign}`).join('; '));
