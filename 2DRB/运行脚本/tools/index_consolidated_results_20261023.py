"""建立合并成果入口，保留各报告原有结论和来源，不重新推算数值。"""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / '文档/工作树成果合并_20261023'


def read_csv(rel):
    with (ROOT / rel).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def main():
    manifest = json.loads((OUT / 'manifest.json').read_text(encoding='utf-8'))
    paths = {
        '当前主目录': '运行脚本/run_cases/test_summary_20261023/全部数值测试.csv',
        '换算法': '文档/多块网格测试总表_20261023/全部测试清单.csv',
        '从粗网格稳态找结果': '运行脚本/run_cases/Multiblock_all_tests_summary_20261023/all_tests.json',
    }
    groups = []
    for name, rel in paths.items():
        if rel.endswith('.csv'):
            records = read_csv(rel)
            extra = {}
        else:
            data = json.loads((ROOT / rel).read_text(encoding='utf-8'))
            records = data['records']
            extra = {k: v for k, v in data.items() if k != 'records'}
        groups.append({'来源': name, '本地清单': str(ROOT / rel), '记录数': len(records), '原始记录': records, '配套检查': extra})
    index = {
        '说明': '三份原始清单按来源合并；同一算例可能重复出现，记录总数不代表独立算例数。保留各源报告时间截点，不把短时完成当作长时完成，不重算或拼接统计。',
        '来源路径迁移': manifest['trees'], '当前根目录': str(ROOT),
        '清单': groups,
    }
    (OUT / '统一测试索引.json').write_text(json.dumps(index, ensure_ascii=False, indent=2), encoding='utf-8')
    report = ROOT / '运行脚本/run_cases/Multiblock_all_tests_summary_20261023/outputs/20261023_all_tests_0f19/全部测试结果与误差.md'
    workbook = report.with_name('全部测试结果.xlsx')
    counts = manifest['counts']
    text = f'''# 两个工作树的实验成果合并

已将“换算法”和“从粗网格稳态找结果”的实验副本、原始诊断证据及测试报告汇集到当前项目目录。原始工作树继续保留。

本次是文件成果合并：未修改求解器、物理参数、数值算法或服务器任务；未创建 Git 提交，也未删除或归档工作树。两边主程序统一换行后与当前主目录内容一致，故保留主目录原文件及既有未提交修改。

## 查看结果

- [当前主聊天：全部数值测试]({(ROOT / paths['当前主目录']).as_posix()})
- [“换算法”：全部测试结果]({(ROOT / '文档/多块网格测试总表_20261023/全部测试结果.md').as_posix()})
- [“从粗网格稳态找结果”：全部测试结果与误差]({report.as_posix()})
- [“从粗网格稳态找结果”：原有 Excel 总表]({workbook.as_posix()})
- [三份清单统一索引]({(OUT / '统一测试索引.json').as_posix()})
- [逐文件来源与 SHA-256 校验清单]({(OUT / 'manifest.json').as_posix()})

## 合并边界

实验目录仍位于 `运行脚本/run_cases/` 下，按原目录名保留，算法之间没有相互替换。测试报告保持原始统计口径和结论；未把不同算法、初值、时间窗口的结果合算。

三份清单分别有 {groups[0]['记录数']}、{groups[1]['记录数']}、{groups[2]['记录数']} 条记录，存在交叉引用与重复算例，不能相加作为独立测试数量。制造场、脉冲和工程检查保留在各原始报告及统一索引的配套检查中。

本次复制 {counts.get('copied', 0)} 个文件；相同文件沿用主目录版本；同名实质差异另存 {counts.get('copied_alternative', 0)} 个文件，见校验清单。原始场、诊断二进制、历史和证据压缩包均保留；跳过依赖缓存和编译产物。各受保护主程序合并前后 SHA-256 一致。

原始快照和脚本中的历史绝对路径没有改写，避免破坏来源证据；路径迁移关系见统一索引。若重用含原工作树绝对路径的脚本，应先确认工作目录。远端目录和运行中作业不受此次本地合并影响，后续新结果仍需取回。

`run_cases` 继续遵循原有 Git 忽略规则；文件已在当前项目磁盘上汇集，并不意味着已写入 Git 历史。本次未重新查询服务器，因此报告中的运行状态仍属于各原报告的时间截点。
'''
    (OUT / 'README.md').write_text(text, encoding='utf-8')
    print(json.dumps({'records': {g['来源']: g['记录数'] for g in groups}, 'entry': str(OUT / 'README.md')}, ensure_ascii=False))


if __name__ == '__main__':
    main()
