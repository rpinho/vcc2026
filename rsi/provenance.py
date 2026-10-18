"""Content identities, with stat-keyed memoization only for read-only large banks."""
import runtime as rt
from pathlib import Path
import hashlib,json

def identity(task=None):
 manifest=rt.json.loads(Path('artifacts/tasks/manifest.json').read_text());cd4=bool(task and manifest[task]['context']=='CD4')
 paths=[Path(p) for p in ['runtime.py','provenance.py','recipes.py','bounded_score.py','official_scoring.py','bounded_de.py','baseline_protocol.py','task_io.py','stat_cache.py','benchmark.py','raw_calibration.py','prepare_sources.py','artifacts/gene_aliases.json','artifacts/promoter_pairs.csv']]
 paths+=sorted(Path('vendor').glob('*.py'))
 paths+=[p for p in sorted(Path('artifacts/cd4_sources' if cd4 else 'artifacts/sources').glob('*.npz')) if p.stem!='RPE1']
 shared=rt.ROOT.parent/'public_data/model_v4_signature_banks';paths+=[p for p in sorted(shared.glob('*.npz')) if (p.stem=='HCT116' if cd4 else p.stem not in ['HEK293T','H1','RPE1'])]
 if cd4:paths += [Path('artifacts/cd4_banks/K562_GWPS_raw.npz')]
 bulk=rt.ROOT.parent/'public_data/K562_GWPS_bulk.h5ad';paths.append(bulk)
 file=Path('artifacts/identity_cache')/((task or 'all')+'.json');old=rt.json.loads(file.read_text()) if file.exists() else {'files':{}};entries={}
 for p in paths:
  name='shared_sources/'+p.name if p.parent==shared else 'shared_inputs/K562_GWPS_bulk.h5ad' if p==bulk else str(p);stat=p.stat();sig=[stat.st_size,stat.st_mtime_ns];prior=old['files'].get(name)
  if prior and prior['stat']==sig:entries[name]=prior
  else:entries[name]=dict(stat=sig,sha256=rt.sha(p))
 inputs=dict(sealed_H1_targets=sorted(manifest['H1_holdout']['targets']),task=task,contract=manifest[task] if task else manifest,task_file_sha256=rt.json.loads(Path('artifacts/tasks',task+'.json').read_text())['sha256'] if task else None)
 hashes={k:v['sha256'] for k,v in entries.items()};hashes['scientific_environment']=__import__('official_scoring').provenance();hashes['scientific_inputs']=inputs;digest=hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest();receipt=dict(files=entries,inputs=inputs,sha256=digest,policy='HEK293T/RPE1 excluded globally; H1 responses excluded or sealed target rows removed; control-only prediction inputs allowed');rt.save(file,receipt);snapshot=Path('artifacts/identities')/(digest+'.json')
 if not snapshot.exists():rt.save(snapshot,receipt)
 return digest
