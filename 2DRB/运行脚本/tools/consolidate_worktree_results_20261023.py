"""只增不覆盖地汇集两个工作树的实验成果；默认仅预览，--apply 才复制。"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / '文档/工作树成果合并_20261023'
TREES = {
    '换算法': Path('C:/Users/User/.codex/worktrees/7d1d/热对流/2DRB'),
    '从粗网格稳态找结果': Path('C:/Users/User/.codex/worktrees/0f19/热对流/2DRB'),
}
SKIP_DIRS = {'node_modules', '__pycache__', '.git', '.venv'}
SKIP_EXT = {'.exe', '.o', '.obj', '.mod', '.smod', '.pyc', '.pyo', '.dll', '.pyd', '.node'}
PROTECTED = [
    '均匀网格/2DRBOpenmp.F90', '均匀网格/3DRBOpenmp.F90',
    '均匀网格/3DRBOpenacc.F90', '均匀网格/3DRBOpenaccMpi.F90',
    '多块网格/2DRBOpenaccMultiblock.F90',
    '多块网格/2DRBOpenaccMultiblockPostprocess.F90',
]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def files_under(base):
    for directory, dirs, names in os.walk(base, followlinks=False):
        dirs[:] = sorted(d for d in dirs if d not in SKIP_DIRS and not (Path(directory)/d).is_symlink())
        for name in sorted(names):
            p = Path(directory) / name
            yield p, ('symlink' if p.is_symlink() else 'compiled_output' if p.suffix.lower() in SKIP_EXT else None)


def normalized_equal(a, b):
    # 同名继承文件若只差换行或 BOM，保留主目录原有字节。
    if max(a.stat().st_size, b.stat().st_size) > 16 * 1024 * 1024:
        return False
    try:
        return a.read_text(encoding='utf-8-sig') == b.read_text(encoding='utf-8-sig')
    except (UnicodeError, ValueError):
        return False


def scan():
    for tag, tree in TREES.items():
        bases = [tree / '运行脚本/run_cases']
        docs = tree / '文档/多块网格测试总表_20261023'
        if docs.exists():
            bases.append(docs)
        for base in bases:
            for source, skip in files_under(base):
                rel = source.relative_to(tree)
                yield tag, source, ROOT / rel, skip


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    candidates = list(scan())
    counts = Counter()
    for tag, source, target, skip in candidates:
        counts['skipped' if skip else 'existing' if target.exists() else 'new'] += 1
        if not skip and not target.exists():
            counts['new_bytes'] += source.stat().st_size
    print(json.dumps({'apply': args.apply, 'plan': dict(counts)}, ensure_ascii=False), flush=True)
    if not args.apply:
        return
    if (OUT / 'manifest.json').exists():
        raise SystemExit('已有完整合并清单；请先核对，禁止覆盖原清单。')
    protected_before = {p: digest(ROOT / p) for p in PROTECTED if (ROOT / p).exists()}
    OUT.mkdir(parents=True, exist_ok=True)
    records = []
    for index, (tag, source, target, skip) in enumerate(candidates):
        record = {'origin': tag, 'source': str(source), 'destination': str(target)}
        if skip:
            record['action'] = 'skipped_' + skip
        else:
            before = source.stat()
            source_hash = digest(source)
            record.update(bytes=before.st_size, source_sha256=source_hash)
            if target.exists():
                target_hash = digest(target)
                if target_hash == source_hash:
                    record['action'] = 'already_identical'
                elif normalized_equal(source, target):
                    record['action'] = 'already_identical_normalized'
                    record['destination_sha256'] = target_hash
                else:
                    target = OUT / '同名差异副本' / tag / source.relative_to(TREES[tag])
                    record['destination'] = str(target)
                    record['action'] = 'copied_alternative'
            else:
                record['action'] = 'copied'
            if record['action'].startswith('copied'):
                if not target.resolve().is_relative_to(ROOT.resolve()):
                    raise RuntimeError(f'目标超出当前工作区: {target}')
                if target.exists():
                    raise FileExistsError(target)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
                target_hash = digest(target)
                if target_hash != source_hash:
                    raise RuntimeError(f'复制后校验失败: {source}')
                record['destination_sha256'] = target_hash
            after = source.stat()
            if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                raise RuntimeError(f'复制期间源文件发生变化: {source}')
        records.append(record)
        if index % 500 == 0:
            print(json.dumps({'checked': index + 1, 'total': len(candidates)}, ensure_ascii=False), flush=True)
    protected_after = {p: digest(ROOT / p) for p in protected_before}
    if protected_before != protected_after:
        raise RuntimeError('主程序保护校验失败')
    result = {
        'created_local': datetime.now().astimezone().isoformat(),
        'scope': '实验副本、原始证据和测试报告的文件合并；无 Git 提交、删除、求解器修改或远端操作',
        'trees': {k: str(v) for k, v in TREES.items()},
        'counts': dict(Counter(r['action'] for r in records)),
        'protected_before': protected_before, 'protected_after': protected_after,
        'excluded_directory_names': sorted(SKIP_DIRS), 'excluded_extensions': sorted(SKIP_EXT),
        'files': records,
    }
    (OUT / 'manifest.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({'finished': True, 'counts': result['counts'], 'manifest': str(OUT / 'manifest.json')}, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
