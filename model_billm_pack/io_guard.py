"""Fail closed on writes outside job directories and reads of sealed task files."""
import sys,os
from pathlib import Path
REPO=Path(__file__).resolve().parent.parent
ROOTS=[REPO/'model_billm',REPO/'model_billm_pack']
def allowed(p):
 if isinstance(p,int):return True
 p=Path(os.fsdecode(p)).resolve()
 return any(p.is_relative_to(r) for r in ROOTS) or str(p)=='/dev/null'
def check(p):
 if p is not None and not allowed(p):raise PermissionError(f'Job write outside authorized directories: {p}')
def audit(event,args):
 if event=='open':
  p,mode,flags=args
  if isinstance(p,(str,bytes)):
   s=os.fsdecode(p).lower()
   if '/artifacts/tasks/' in s and ('holdout' in s or 'sealed' in s):raise PermissionError(f'Sealed HOLDOUT read forbidden: {p}')
  if (mode and any(c in mode for c in 'wax+')) or flags & (os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC):check(p)
 elif event in ['os.mkdir','os.remove','os.rmdir','os.chmod','os.utime']:check(args[0])
 elif event in ['os.rename','os.link','os.symlink']:
  check(args[0]);check(args[1])
def install():sys.addaudithook(audit)
