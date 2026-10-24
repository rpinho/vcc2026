"""Local unpack of the operator-downloaded final bundle; record input/output hashes."""
import common as job
import io_guard;io_guard.install()
import argparse,zipfile,tarfile,shutil
from pathlib import Path

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--archive',required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
    dest=Path(a.out).resolve();assert dest.is_relative_to(job.ROOT);dest.mkdir(parents=True,exist_ok=True)
    copied={}
    def copy(name,src):
        base=Path(name).name
        if not (base.endswith('.h5ad') or base in ['pert_counts.csv','gene_names.csv']):return
        assert base not in copied,'Repeated bundle basename: '+base
        target=dest/base
        with target.open('wb') as out:shutil.copyfileobj(src,out,8*1024**2)
        copied[base]=dict(bytes=target.stat().st_size,sha256=job.sha(target))
    if zipfile.is_zipfile(a.archive):
        with zipfile.ZipFile(a.archive) as archive:
            for member in archive.infolist():
                if not member.is_dir():
                    with archive.open(member) as f:copy(member.filename,f)
    else:
        with tarfile.open(a.archive) as archive:
            for member in archive:
                if member.isfile():
                    with archive.extractfile(member) as f:copy(member.name,f)
    assert 'pert_counts.csv' in copied and any(n.endswith('.h5ad') for n in copied)
    job.save(dest/'download_unpack_receipt.json',dict(archive_sha256=job.sha(a.archive),files=copied))
    print('LOCAL UNPACK',copied,flush=True)
if __name__=='__main__':main()
