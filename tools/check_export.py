"""Offline source/privacy checks. Does not open data or credential stores."""
from pathlib import Path
import ast
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
allowed_suffixes = {'.py', '.sh', '.md', '.txt', '.sha256'}
blocked_dirs = {'_codex', '__pycache__', 'artifacts', 'logs', 'ledger', 'locks',
                'journals', 'raw', 'data', 'public_data', 'cache', 'v2_cache', 'out'}
patterns = {
    'personal absolute path': re.compile(r'(?<![A-Za-z0-9.])/(?:Users|home)/[A-Za-z0-9_.-]+'),
    'email': re.compile(r'(?<![\w.+-])[\w.+-]+@[\w.-]+\.[A-Za-z]{2,}(?![\w.-])'),
    'private key': re.compile(r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----'),
    'credential literal': re.compile(r'''(?i)\b(?:api[_-]?key|client[_-]?secret|access[_-]?token|password|passwd|vcc[_-]?token)\b\s*[:=]\s*["']([^"'\n]{8,})["']'''),
}

def main():
    findings=[];count=0;python_count=0;size=0
    for p in sorted(ROOT.rglob('*')):
        rel=p.relative_to(ROOT)
        if '.git' in rel.parts:continue
        if p.is_symlink():
            findings.append((str(rel),'symlink'));continue
        if p.is_dir():
            if p.name in blocked_dirs and str(rel) != 'rsi/vendor/sidechain_release/src/sidechain/data':
                findings.append((str(rel),'runtime/data directory'))
            continue
        if not p.is_file():continue
        count+=1;size+=p.stat().st_size
        if p.suffix not in allowed_suffixes and p.name not in {'.gitignore','LICENSE','ATLAS_LICENSE','SIDECHAIN_LICENSE','CELL_EVAL2_LICENSE'}:
            findings.append((str(rel),'unexpected file type'))
        try:text=p.read_text()
        except UnicodeDecodeError:
            findings.append((str(rel),'binary file'));continue
        if p.suffix=='.py':
            python_count+=1
            try:ast.parse(text,filename=str(rel))
            except SyntaxError as e:findings.append((str(rel),f'syntax line {e.lineno}'))
        for label,pattern in patterns.items():
            for match in pattern.finditer(text):
                # Log location/type only; never echo suspected credentials.
                line=text.count('\n',0,match.start())+1
                findings.append((str(rel),f'{label} line {line}'))
    print(f'Files scanned: {count}; Python sources parsed: {python_count}; bytes: {size}')
    for path,kind in findings:print(f'FINDING {path}: {kind}')
    print(f'Privacy/structure findings: {len(findings)}')
    print('PASS' if not findings else 'FAIL')
    return bool(findings)

if __name__=='__main__':sys.exit(main())
