"""Independent static evidence and isolation audit. Never execute fetched code."""
import argparse
import ast
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path

ROOT = Path('artifacts/responder_sources_round47')
PLAN = ROOT / 'audit_plan.json'
OUT = ROOT / 'audit.json'
COMMIT = '1d6d48a409afbbaead3e757203bccf6c883d81ab'
PIE = Path('../research/sources/ArcInstitute__pie')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def blob(data):
    return hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()


def save_new(path, value):
    with path.open('x') as f:
        f.write(json.dumps(value, indent=2) + '\n')


def freeze():
    paths = [Path(__file__), Path('acquire_responder_sources_round47.py'),
        ROOT/'ASSESSMENT.md', ROOT/'commit.json', ROOT/'tree.json', ROOT/'receipt_17.json',
        Path('artifacts/Responder_primary_code_round47_charges.jsonl'),
        Path('artifacts/download_allocations.json'),
        Path('artifacts/calibration_gate.json'), Path('artifacts/tasks/manifest.json'),
        Path('artifacts/benchmark_roles.json'), Path('artifacts/frozen_metadata_index.json'),
        Path('../research/sources/pie_tree.json'), Path('../research/sources/pie_readme.md'),
        Path('../research/sources/slim_readme.md')]
    paths += sorted(p for p in (ROOT/'release').rglob('*') if p.is_file())
    paths += sorted(p for p in PIE.rglob('*') if p.is_file())
    for task in ['H1_rich_search', 'CD4_depth192_search']:
        for recipe in ['v1', 'v2b', 'v3b', 'v4a', 'v4b']:
            paths.append(Path('artifacts/scores')/task/f'{recipe}_s42.json')
    save_new(PLAN, dict(utc=datetime.now(timezone.utc).isoformat(),
        scope='Static primary-source usability/isolation certificate only; no predictions or model metrics',
        sha256={str(p): sha(p) for p in paths}))
    print('FROZEN', PLAN)


def main():
    plan = read(PLAN)
    for path, digest in plan['sha256'].items():
        assert sha(path) == digest, path
    tree = read(ROOT/'tree.json')
    assert tree['sha'] == read(ROOT/'commit.json')['sha'] == COMMIT
    assert tree['truncated'] is False
    entries = {x['path']: x for x in tree['tree'] if x['type'] == 'blob'}
    assert len(entries) == 33
    license_names = [p for p in entries if Path(p).name.lower().startswith(('license','copying'))]
    assert not license_names
    receipt = read(ROOT/'receipt_17.json')
    assert receipt['commit'] == COMMIT and len(receipt['files']) == 17
    assert receipt['tree_sha256'] == sha(ROOT/'tree.json')
    for name, meta in receipt['files'].items():
        data = (ROOT/'release'/name).read_bytes()
        assert len(data) == entries[name]['size'] == meta['bytes']
        assert blob(data) == entries[name]['sha'] == meta['git_blob']
        assert sha(ROOT/'release'/name) == meta['sha256']
    scripts = {}
    for name in receipt['files']:
        if name.endswith('.py'):
            source = (ROOT/'release'/name).read_text()
            scripts[name] = ast.parse(source)  # parse only, never import
    source = (ROOT/'release/scripts/04_robustness_continuous_magnitude.py').read_text()
    assignments = {}
    for node in scripts['scripts/04_robustness_continuous_magnitude.py'].body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            if node.targets[0].id in ['pred_bio', 'pred_tech']:
                assignments[node.targets[0].id] = ast.literal_eval(node.value)
    assert assignments == dict(pred_bio=['max_paralog_identity_filled','log_string_degree_400','LOEUF'],
        pred_tech=['kd_strength','log_ncells','log_baseMean'])
    for literal in ['smf.ols(', 'smf.logit(', 'mc.groupby("culture_condition")',
        's.std(ddof=0)', 'mc["n_total_de_genes"] <= 2', 'mc["n_total_de_genes"] <= 100',
        'Zt[rows_self, self_idx[has_self]] = 0.0']:
        assert literal in source, literal
    tables = {}
    terms = set(x+'_z' for values in assignments.values() for x in values)
    for name in ['robust_continuous_magnitude_fit.csv', 'robust_detectability_logistic.csv']:
        with (ROOT/'release/results'/name).open() as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == 18
        assert {r['condition'] for r in rows} == {'Rest','Stim8hr','Stim48hr'}
        for cond in ['Rest','Stim8hr','Stim48hr']:
            assert {r['term'] for r in rows if r['condition'] == cond} == terms
        assert all(math.isfinite(float(r['beta'])) for r in rows)
        assert not any(r['term'].lower() in ['intercept','const'] for r in rows)
        tables[name] = dict(rows=len(rows), conditions=3, slopes_per_condition=6, intercept_rows=0)
    p_tree = read('../research/sources/pie_tree.json')
    assert p_tree['sha'] == '21a526add9ff3a81b4f962d9b16bbae570052c63' and not p_tree['truncated']
    p_entries = {x['path']: x for x in p_tree['tree'] if x['type'] == 'blob'}
    p_files = {'README.md': Path('../research/sources/pie_readme.md')}
    p_files.update({p.relative_to(PIE).as_posix():p for p in PIE.rglob('*') if p.is_file()})
    for name,p in p_files.items():
        assert blob(p.read_bytes()) == p_entries[name]['sha'], name
    assert len(p_files) == 7
    assert not any(p.endswith(('.ckpt','.safetensors','.pt','.pth')) for p in p_entries)
    assert 'CC BY-NC-SA 4.0' in p_files['README.md'].read_text()
    assert 'PREDICTION_COLUMNS: tuple[str, ...] = ("p_de", "lfc_pred", "delta_p_pred")' in p_files['src/pie/predict.py'].read_text()
    assert 'torch.load(ckpt_path' in p_files['src/pie/predict.py'].read_text()
    assert sha('../research/sources/slim_readme.md') == '6e8d10d4c0b7e37e3a67a26a127c4cb2a21ebec701aad537a5b5b670453f8e75'
    gate = read('artifacts/calibration_gate.json')
    assert gate['status'] == 'CLOSED' and gate['phase2_allowed'] is False
    assert gate['macro']['spearman'] == 0.8999999999999998
    frozen = read('artifacts/frozen_metadata_index.json')
    for path,entry in frozen['entries'].items():
        assert sha(path) == entry['sha256'], path
    assert sha(frozen['archive']) == frozen['archive_sha256']
    allocations = read('artifacts/download_allocations.json')
    assert sum(allocations['new_allocations'].values()) == allocations['total_allocated_upper_bound_bytes'] <= 15_000_000_000
    journal = [json.loads(s) for s in Path('artifacts/Responder_primary_code_round47_charges.jsonl').read_text().splitlines()]
    charged = 0
    for row in journal:
        if row['event'] == 'charge':
            charged += row['bytes']
            assert row['total'] == charged <= 5_000_000
    assert charged == receipt['charged_bytes']
    steps = [json.loads(s) for s in Path('ledger/steps.jsonl').read_text().splitlines()]
    for name in ['inspect_responder_index_round47','inspect_responder_release_round47']:
        ends = [s for s in steps if s['event']=='end' and s['name']==name]
        assert ends and ends[-1]['rc']==0
    save_new(OUT, dict(utc=datetime.now(timezone.utc).isoformat(), status='PASS static source identity, fitted-feature and isolation certificate',
        plan_sha256=sha(PLAN), assessment_sha256=sha(ROOT/'ASSESSMENT.md'),
        primary_commit=COMMIT, primary_tree_blobs=len(entries), primary_verified_files=17,
        license_file_names=license_names, parsed_scripts=len(scripts), coefficient_tables=tables,
        cached_pie_verified_files=len(p_files), slim_immutable_release_identity='UNVERIFIED cached README only',
        charged_bytes=charged, total_allocated_upper_bound_bytes=allocations['total_allocated_upper_bound_bytes'],
        frozen_members_verified=len(frozen['entries']), model_predictions=0, biological_scores=0,
        destination_family_fit_permitted_for_cd4=False, scientific_gate='CLOSED unchanged',
        verdict='No complete released permissible responder-gating inference contract established; no predictor implemented',
        limitations='Does not verify study coefficients against raw responses, external model availability, scientific generalization or every possible license statement outside pinned sources'))
    print('PASS', OUT)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze', action='store_true')
    if parser.parse_args().freeze:
        freeze()
    else:
        main()
