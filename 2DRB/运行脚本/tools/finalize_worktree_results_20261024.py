"""在20261023已校验汇集基础上，只增不覆盖地取回最终结果和新增文件。"""
from pathlib import Path
from datetime import datetime
from collections import Counter
import json
import shutil
import subprocess
import consolidate_worktree_results_20261023 as prior

ROOT = prior.ROOT
OUT = ROOT / '文档/工作树成果合并_20261024'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT / 'incremental_manifest.json').exists():
        raise SystemExit('最终汇集已完成；先核对原清单，不覆盖。')
    old = json.loads((prior.OUT / 'manifest.json').read_text(encoding='utf-8'))
    old_records = {r['source']: r for r in old['files']}
    protected = {p: prior.digest(ROOT / p) for p in old['protected_before']}
    status = subprocess.check_output(['git', '-c', 'core.quotepath=false', 'status', '--short'], cwd=ROOT.parent, encoding='utf-8')
    (OUT / 'premerge_status.json').write_text(json.dumps({'status': status, 'protected': protected}, ensure_ascii=False, indent=2), encoding='utf-8')
    # 备份既有差异，避免后续Git操作掩盖用户已有内容。
    for name, extra in [('premerge_working.patch', []), ('premerge_staged.patch', ['--cached'])]:
        (OUT / name).write_bytes(subprocess.check_output(['git', 'diff', '--binary', *extra], cwd=ROOT.parent))
    candidates = list(prior.scan())
    for tag, tree in prior.TREES.items():
        latest = tree / '文档/两组长算结果_20261024'
        if latest.exists():
            for source, skip in prior.files_under(latest):
                candidates.append((tag, source, ROOT / source.relative_to(tree), skip))
    counts = Counter()
    delta = []
    for i, (tag, source, target, skip) in enumerate(candidates):
        if skip:
            counts['excluded_compiled_or_link'] += 1
            continue
        saved = old_records.get(str(source))
        if saved and not saved['action'].startswith('skipped'):
            archived = Path(saved['destination'])
            if archived.is_file():
                a, b = source.stat(), archived.stat()
                if (a.st_size, a.st_mtime_ns) == (b.st_size, b.st_mtime_ns):
                    counts['previously_verified_unchanged_metadata'] += 1
                    continue
        source_hash = prior.digest(source)
        if saved and saved.get('source_sha256') == source_hash and Path(saved['destination']).is_file():
            expected = saved.get('destination_sha256', source_hash)
            if prior.digest(Path(saved['destination'])) == expected:
                counts['previously_verified_hash_match'] += 1
                continue
        record = {'origin': tag, 'source': str(source), 'source_sha256': source_hash, 'bytes': source.stat().st_size}
        if target.exists() and prior.digest(target) == source_hash:
            counts['already_identical'] += 1
            continue
        if target.exists():
            if prior.normalized_equal(source, target):
                counts['already_identical_normalized'] += 1
                continue
            target = OUT / '同名差异副本' / tag / source.relative_to(prior.TREES[tag])
            record['action'] = 'copied_alternative'
        else:
            record['action'] = 'copied'
        if not target.resolve().is_relative_to(ROOT.resolve()):
            raise RuntimeError(target)
        if target.exists():
            raise FileExistsError(target)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        record.update(destination=str(target), destination_sha256=prior.digest(target))
        assert record['destination_sha256'] == source_hash
        delta.append(record)
        counts[record['action']] += 1
        if i % 1000 == 0:
            print(json.dumps({'progress': i, 'total': len(candidates)}), flush=True)
    after = {p: prior.digest(ROOT / p) for p in protected}
    assert protected == after, '保护的主程序发生变化'
    result = {'date_Asia_Shanghai': '2026-10-24', 'base_manifest': str(prior.OUT / 'manifest.json'),
              'counts': dict(counts), 'protected_before': protected, 'protected_after': after,
              'files': delta, 'source_trees_preserved': True,
              'verification': '新增文件逐项SHA-256；已归档文件按原校验清单和未变元数据复用，元数据改变则重新校验。'}
    (OUT / 'incremental_manifest.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({'complete': True, 'counts': dict(counts), 'new_bytes': sum(x['bytes'] for x in delta)}, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
