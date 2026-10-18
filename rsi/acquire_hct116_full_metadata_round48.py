"""Pin and acquire complete HCT116 cell metadata; no expression endpoints."""
from datetime import datetime, timezone
import hashlib
import os
from pathlib import Path
import re

import runtime as rt
import requests
from durable_public_store import PublicStore
from transfer_budget import TransferBudget

ROOT = Path('artifacts/HCT116_full_cells_metadata_round48')
REVISION = '598aa5442a4ee9a8f1a383a843ab5ce142904d04'
PREFIX = 'data/HCT116/cells.lance/'
KEY = 'HCT116_full_cells_metadata_round48'
JOURNAL = Path('artifacts/'+KEY+'_charges.jsonl')
REPO = 'slaf-project/X-Atlas-Orion'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def verify(path, entry):
    data = path.read_bytes()
    assert len(data) == entry['size'], 'File size mismatch'
    if entry.get('lfs'):
        assert digest(data) == entry['lfs']['oid'], 'LFS content identity mismatch'
    else:
        assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest() == entry['oid'], 'Git blob mismatch'
    return digest(data)


def fetch_file(entry, store, budget, destination):
    """Full attempts are charged before GET; failed tmp files never count as done."""
    if destination.exists():
        return verify(destination, entry)
    url, identity = store.resolve(entry['path'])
    assert identity['size'] == entry['size']
    amount = entry['size']
    budget.charge(amount, 'Full HCT116 cells metadata GET attempt '+entry['path'])
    tmp = destination.with_suffix(destination.suffix+'.partial')
    destination.parent.mkdir(parents=True, exist_ok=True)
    with store.http.get(url, headers={'Accept-Encoding':'identity','If-Match':identity['etag']},
                        stream=True, allow_redirects=False, timeout=(10,90)) as response:
        response.raise_for_status()
        assert response.status_code == 200
        assert int(response.headers['Content-Length']) == amount
        assert response.headers['ETag'] == identity['etag']
        assert response.headers.get('Content-Encoding','identity') == 'identity'
        with tmp.open('wb') as f:
            received = 0
            while received < amount:
                block = response.raw.read(min(1024*1024, amount-received))
                assert block, 'Truncated metadata response'
                received += len(block)
                f.write(block)
            f.flush();os.fsync(f.fileno())
    result = verify(tmp, entry)
    tmp.replace(destination)
    return result


def main():
    ROOT.mkdir(exist_ok=True)
    allocation = rt.json.loads(Path('artifacts/download_allocations.json').read_text())
    assert allocation['new_allocations'][KEY] == 180_000_000
    scope = rt.json.loads(Path('artifacts/HCT116_census_scope_audit_round48.json').read_text())
    assert scope['status'].startswith('PASS')
    assert scope['excluded_name_counts']['cd4_and_full_raw_gwps_ge200_not_queried'] > 0
    budget = TransferBudget(JOURNAL,180_000_000)
    budget.bind('repository_revision',REVISION)
    budget.bind('reader_sha256',rt.sha(__file__))
    budget.bind('scope_audit_sha256',rt.sha('artifacts/HCT116_census_scope_audit_round48.json'))
    store = None
    try:
        inventory = ROOT/'inventory.json'
        if not inventory.exists():
            cap = 512_000
            budget.charge(cap,'Bounded complete pinned cells-directory inventory GET attempt')
            session=requests.Session();session.trust_env=False
            try:
                url=f'https://huggingface.co/api/datasets/{REPO}/tree/{REVISION}/data/HCT116/cells.lance?recursive=true&limit=1000'
                with session.get(url,headers={'Accept-Encoding':'identity'},stream=True,
                                 allow_redirects=False,timeout=(10,60)) as response:
                    response.raise_for_status()
                    assert response.status_code==200 and not response.headers.get('Link'), 'Paginated inventory requires separate plan'
                    body=response.raw.read(cap,decode_content=True)
                    assert len(body)<cap
                    entries=rt.json.loads(body)
                    assert isinstance(entries,list) and len(entries)<1000
                    for e in entries:
                        assert e['path'].startswith(PREFIX) and '..' not in e['path']
                rt.save(inventory,dict(revision=REVISION,complete=True,entries=entries,
                    url=url,utc=datetime.now(timezone.utc).isoformat()))
            finally:
                session.close()
        inv=rt.json.loads(inventory.read_text())
        assert inv['revision']==REVISION and inv['complete'] is True
        budget.bind('inventory_sha256',rt.sha(inventory))
        files=[e for e in inv['entries'] if e['type']=='file' and '/_transactions/' not in e['path']]
        assert files and any('/_versions/' in e['path'] for e in files)
        for e in files:
            assert e['path'].startswith(PREFIX) and '..' not in e['path']
            assert re.fullmatch('[0-9a-f]{40}',e['oid']) and e['size']>0
        outstanding=sum(e['size']+64_000 for e in files if not (ROOT/'release/cells.lance'/e['path'][len(PREFIX):]).exists())
        assert budget.total+outstanding<=budget.limit, 'Complete metadata exceeds remaining attempt budget'
        rt.guard(sum(e['size'] for e in files)+512*1024**2)
        store=PublicStore(budget)
        receipts={}
        for e in files:
            destination=ROOT/'release/cells.lance'/e['path'][len(PREFIX):]
            h=fetch_file(e,store,budget,destination)
            receipts[e['path']]=dict(size=e['size'],sha256=h,git_oid=e['oid'],lfs=e.get('lfs'),path=str(destination))
            rt.save(ROOT/'progress.json',dict(revision=REVISION,inventory_sha256=rt.sha(inventory),files=receipts,charged_bytes=budget.total))
            print('VERIFIED',e['path'],e['size'],'charged',budget.total,flush=True)
        rt.save(ROOT/'receipt.json',dict(status='ACQUIRED immutable metadata files; independent census audit pending',
            utc=datetime.now(timezone.utc).isoformat(),revision=REVISION,inventory_sha256=rt.sha(inventory),
            reader_sha256=rt.sha(__file__),files=receipts,charged_bytes=budget.total,
            journal=str(JOURNAL),expression_payloads_read=False))
    finally:
        if store:store.close()
        budget.close()


if __name__=='__main__':
    main()
