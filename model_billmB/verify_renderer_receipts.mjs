// Separate exact receipt recheck. Raw counts are independently audited by
// model_v9/audit_worker.py and bound to emitted/decoded cells by package_checked.py.
import fs from 'node:fs';import path from 'node:path';import crypto from 'node:crypto';
const P=import.meta.dirname,full=process.argv.includes('--full');
const file=path.join(P,full?'full_gates.json':'renderer_preflight.json');
const bytes=fs.readFileSync(file),g=JSON.parse(bytes),expected=full?900:6;
if(g.records.length!==expected||g.status!==(full?'PASS all900':'PASS probes only; full900 still required')||(full&&!g.full900))throw Error('Incomplete or failed gates');
function bits(x){const v=new DataView(new ArrayBuffer(8));v.setFloat64(0,x);return v.getBigUint64(0);}
function fromBits(x){const v=new DataView(new ArrayBuffer(8));v.setBigUint64(0,x);return v.getFloat64(0);}
function rational(x){
 if(!Number.isFinite(x)||x<0)throw Error('Invalid nonnegative scalar');
 const b=bits(x),e=Number((b>>52n)&2047n),m=(b&((1n<<52n)-1n))+(e?1n<<52n:0n),p=e?e-1023-52:-1074;
 return p>=0?[m<<BigInt(p),1n]:[m,1n<<BigInt(-p)];
}
const fraction=x=>{const [n,d]=rational(x);return `${n}/${d}`;};
const above=(a,b)=>{const [an,ad]=rational(a),[bn,bd]=rational(b);return an*bd>=bn*ad;};
function enclosedBelow(x,n,d){
 // The integer verifier converts Fraction -> binary64 using nearest rounding.
 // One representable step upward therefore encloses the true rational bound,
 // including the case where that conversion rounded downward.
 rational(x);const upper=fromBits(bits(x)+1n),[un,ud]=rational(upper);
 if(!(un*d<n*ud))throw Error('Outward-enclosed moment bound fails exact threshold');
 return upper;
}
let maxMean=0,maxBulk=0,minJK=Infinity;
for(let k=0;k<expected;k++){
 const r=g.records[k],context='ABC'[Math.floor(k/(full?300:2))],index=full?k%300:[0,4][k%2];
 if(r.context!==context||r.index!==index||!r.all_gates_pass||!r.JK_dominance||!['totals_exact','depths_exact','nonDE_exact','own_target_exact'].every(n=>r[n]===true))throw Error('Group identity/invariant failure');
 if(!above(r.jk,r.v5b_jk))throw Error('Exact binary64 JK dominance failure');
 if(!/^[0-9a-f]{64}$/.test(r.counts_sha256)||!/^[0-9a-f]{64}$/.test(r.baseline_counts_sha256))throw Error('Missing count identities');
 const mean=enclosedBelow(r.mean_certificate.mean_L1_upper,15n,1000n),bulk=enclosedBelow(r.mean_certificate.bulk_L1_upper,15n,10000n);
 maxMean=Math.max(maxMean,mean);maxBulk=Math.max(maxBulk,bulk);minJK=Math.min(minJK,r.jk/r.v5b_jk);
}
const receipt={utc:new Date().toISOString(),status:`PASS exact ${expected}-group receipt recheck with outward-enclosed moment bounds`,scope:'Receipt arithmetic/identity check only; original independent count/moment/official-kernel audit and full emitted/decoded count hashes remain required.',groups:expected,full900:full,gates_sha256:crypto.createHash('sha256').update(bytes).digest('hex'),thresholds:{mean:'15/1000',bulk:'15/10000'},max_mean_upper_outward:maxMean,max_mean_upper_outward_exact:fraction(maxMean),max_bulk_upper_outward:maxBulk,max_bulk_upper_outward_exact:fraction(maxBulk),min_JK_ratio:minJK};
fs.writeFileSync(path.join(P,full?'full_receipt_verification.json':'probe_receipt_verification.json'),JSON.stringify(receipt,null,2)+'\n');
console.log(receipt.status);
