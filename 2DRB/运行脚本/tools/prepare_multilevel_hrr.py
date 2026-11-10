"""创建独立的三级 HRR 副本；不写入已有求解器或历史算例。"""
from pathlib import Path
import difflib
import hashlib
import json
import re
import shutil

ROOT = Path(__file__).resolve().parents[2]
NAME = 'Multilevel_HRR0995_Ra9_Ra10_N2048_20261103'
OUT = ROOT / '运行脚本/run_cases' / NAME
GRID_PARENT = ROOT / '多块网格/2DRBOpenaccMultilevel.F90'
HRR_PARENT = ROOT / ('运行脚本/run_cases/'
    'Multiblock_Ra1e9_N2048_HRR0995_T1000_node07_20261027/hrr995/source/solver.F90')
NEW = ROOT / '多块网格/2DRBOpenaccMultilevelHRR.F90'
POST = ROOT / '多块网格/2DRBOpenaccMultilevelPostprocess.F90'
REMOTE = '/data2/XLLi/Multiblock/' + NAME


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding='utf-8', newline='\n')


def save(path, data):
    write(path, json.dumps(data, ensure_ascii=False, indent=2) + '\n')


def proc(text, name):
    matches = re.findall(r'^  subroutine ' + name + r'\b.*?^  end subroutine ' + name + r'\b',
                         text, flags=re.M | re.S)
    assert matches, name
    return matches[-1]


def assign(text, name, value):
    text, count = re.subn(r'(\b' + name + r'\s*=\s*)[\d.]+(?:[dDeE][+-]?\d+)?',
                          lambda m: m[1] + str(value), text, count=1)
    assert count == 1, name
    return text


def once(text, old, new):
    assert text.count(old) == 1, old[:120]
    return text.replace(old, new)


def main():
    assert not NEW.exists(), 'Do not overwrite an existing solver'
    assert not OUT.exists(), 'Do not overwrite a prepared or submitted campaign'
    OUT.mkdir(parents=True)
    protected_paths = list((ROOT / '多块网格').glob('*.F90'))
    protected_paths += [ROOT / '均匀网格' / n for n in [
        '2DRBOpenmp.F90', '3DRBOpenmp.F90', '3DRBOpenacc.F90', '3DRBOpenaccMpi.F90']]
    protected_paths += [HRR_PARENT,
        ROOT / '运行脚本/run_cases/Multiblock_Ra1e10_HRR0995_N2048_N4096_node05_20261027/N2048/source/solver.F90',
        ROOT / '运行脚本/run_cases/Multiblock_Ra1e10_HRR0995_N2048_N4096_node05_20261027/N4096/source/solver.F90']
    protected = {str(p.relative_to(ROOT)): sha(p) for p in protected_paths}
    save(OUT / 'protected_originals_before.json', protected)
    assert sha(HRR_PARENT) == '1d7b4cccde8f268aac3a62c9a631c2e1f1d83127676725faa7d5e83e24f3012a'
    grid = GRID_PARENT.read_text(encoding='utf-8-sig')
    hrr = HRR_PARENT.read_text(encoding='utf-8-sig')
    s = grid
    s = once(s, '!!! 不使用 RR、HRR、魔法参数反算或额外滤波。',
        '!!! 独立 HRR 副本：所有级别的全流场使用 H-RR，sigma=0.995；温度仍为原 Luo D2Q5。\n'
        '!!! 网格管理母版：2DRBOpenaccMultilevel.F90；流场算法移自已完成 Ra9 的双级 HRR 副本。\n'
        '!!! 只在本新文件内修改；已有源码与历史输出保持不变。')
    s = once(s, '!   自定义宏，一些选项的开关',
        '!   自定义宏，一些选项的开关\n#define REG_HRR\n'
        '#if defined(REG_CONTROL) || defined(REG_PR) || defined(REG_H90)\n'
        '#error "This independent solver requires H-RR only"\n#endif')
    s = once(s, '!#define SideHeatedHa', '!#define SideHeatedHa\n#ifdef SideHeatedHa\n'
        '#error "The inherited H-RR collision has no magnetic-force implementation"\n#endif')
    s = s.replace('默认三级：壁面1024级 -> 中间512级 -> 中心256级。',
                  '默认三级：壁面2048级 -> 中间1024级 -> 中心512级。')
    s = s.replace('64层细格作为初始选择', '128层细格作为本次选择')
    s = s.replace('63.5与959.5', '127.5与1919.5')
    s = s.replace('默认第2级32格对应64个最细格距；第3级32格对应128个最细格距。',
                  '本例第2级64格对应128个最细格距；最粗第3级为中心块，无环带厚度参数。')
    for key, value in dict(nx=2048, ny=2048, nGridLevels=3,
            finestLayerCellsLeftRight=128, finestLayerCellsBottomTop=128,
            innerLayerCellsLeftRight=64, innerLayerCellsBottomTop=64,
            Rayleigh='1.0d9').items():
        s = assign(s, key, value)
    s = once(s, 'real(kind=8), parameter :: pi=acos(-1.0d0)',
        'real(kind=8), parameter :: hrrMixingSigma=0.995d0    ! 全流场HRR二阶应力混合系数\n'
        '        real(kind=8), parameter :: pi=acos(-1.0d0)')
    s = once(s, 'Snu=1.0d0/tauf, Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)',
                  'Snu=1.0d0/tauf, Sq=Snu    ! 与双级HRR一致：接口三阶矩按同一松弛率缩放')
    s = once(s, 'real(kind=8), parameter :: outputSnapshotInterval=10.0d0',
        'real(kind=8), parameter :: NuReSampleInterval=10.0d0\n'
        '        real(kind=8), parameter :: outputSnapshotInterval=10.0d0')
    s = once(s, 'real(kind=8), parameter :: outputSnapshotInterval=0.5d0, unsteadyRunDuration=1000.0d0',
        'real(kind=8), parameter :: NuReSampleInterval=0.5d0    ! Nu/Re与有限性检查间隔，单位TFF\n'
        '        real(kind=8), parameter :: outputSnapshotInterval=5.0d0    ! 全场快照独立保存\n'
        '        real(kind=8), parameter :: unsteadyRunDuration=1000.0d0')
    s = once(s, 'integer(kind=4), parameter :: outputSnapshotIntervalItc= &',
        'integer(kind=4), parameter :: NuReSampleIntervalItc= &\n'
        '            max(1,ceiling(NuReSampleInterval*timeUnit/dble(coarsestStep)))*coarsestStep\n'
        '        integer(kind=4), parameter :: outputSnapshotIntervalItc= &')
    s = once(s, 'dble(outputSnapshotIntervalItc))*outputSnapshotIntervalItc',
                  'dble(NuReSampleIntervalItc))*NuReSampleIntervalItc')
    s = s.replace("settingsFile='SimulationSettings2DOpenaccMultilevel.txt'",
                  "settingsFile='SimulationSettings2DOpenaccMultilevelHRR.txt'")
    s = s.replace("reloadFilePrefix='reloadFile2DOpenaccMultilevel'",
                  "reloadFilePrefix='reloadFile2DOpenaccMultilevelHRR'")
    s = once(s, 'min(outputSnapshotInterval,reloadFileInterval,outputPltFileInterval)',
                  'min(NuReSampleInterval,outputSnapshotInterval,reloadFileInterval,outputPltFileInterval)')
    s = once(s, "if(nGridLevels<1) error stop 'nGridLevels must be positive'",
        "if(nGridLevels<1) error stop 'nGridLevels must be positive'\n"
        "        if(.not.ieee_is_finite(hrrMixingSigma).or.hrrMixingSigma<0.0d0.or.hrrMixingSigma>1.0d0) &\n"
        "            error stop 'HRR sigma must be finite and within [0,1]'")
    old = '''            if(mod(reloadItcOffset+itc,outputSnapshotIntervalItc)==0) then
                call calNuRe()
                if(outputSnapshotFile==1) then
                    call output_SnapshotFile()
                    lastSnapshotItc=itc
                endif
            endif'''
    new = '''            if(mod(reloadItcOffset+itc,NuReSampleIntervalItc)==0) call calNuRe()
            if(outputSnapshotFile==1.and.mod(reloadItcOffset+itc,outputSnapshotIntervalItc)==0) then
                call output_SnapshotFile()
                lastSnapshotItc=itc
            endif'''
    s = once(s, old, new)
    s = once(s, "'Original D2Q9 + Luo D2Q5; no regularization; raw coarse/fine moment scaling'",
        "'H-RR D2Q9 on every block/level + original Luo D2Q5; raw coarse/fine moment scaling'\n"
        "        write(k,'(a,es24.16e3)') 'HRR mixing sigma: ',hrrMixingSigma\n"
        "        write(k,'(a)') 'HRR: current-force velocity, local second-order gradients, third-order recursive reconstruction'\n"
        "        write(k,'(a)') 'Grid parent SHA256: " + sha(GRID_PARENT) + "'\n"
        "        write(k,'(a)') 'HRR parent SHA256: " + sha(HRR_PARENT) + "'")
    s = once(s, "write(k,'(a,4(i0,1x))') 'Sampling,restart,Tecplot,convergence intervals (fine steps): ', &\n"
        '            outputSnapshotIntervalItc,reloadFileIntervalItc,outputPltFileIntervalItc,steadyCheckIntervalItc',
        "write(k,'(a,5(i0,1x))') 'Sampling,snapshot,restart,Tecplot,convergence intervals (fine steps): ', &\n"
        '            NuReSampleIntervalItc,outputSnapshotIntervalItc,reloadFileIntervalItc, &\n'
        '            outputPltFileIntervalItc,steadyCheckIntervalItc')
    # 相同的三阶平衡态用于碰撞、pack、unpack。保留原逐矩非平衡缩放。
    for name in ['pack_interface_values', 'unpack_interface_values']:
        body = proc(s, name)
        replacement, count = re.subn(r'    meq\(0\)=rh\n.*?    meq\(8\)=rh\*ux\*uy',
            '    call regularized_equilibrium(rh,ux,uy,meq)', body, flags=re.S)
        assert count == 1
        s = once(s, body, replacement)
    s = once(s, proc(s, 'collision'), proc(hrr, 'collision') + '\n\n' +
        proc(hrr, 'regularized_node') + '\n\n' + proc(hrr, 'regularized_equilibrium'))
    # 使用独立重启版本，并在头部写入sigma；不能误读未正则化/不同sigma的检查点。
    s = s.replace('ML2DRESTART0001', 'MLHRRESTART0001')
    s = once(s, 'total,NuReSampleCount,errorU,errorT,nBlocks\n',
                  'total,NuReSampleCount,errorU,errorT,nBlocks,hrrMixingSigma\n')
    s = once(s, 'real(kind=8) :: physics(12),position(3),savedU,savedT',
                  'real(kind=8) :: physics(12),position(3),savedU,savedT,savedSigma')
    s = once(s, 'magic,config,physics,total,samples,savedU,savedT,nb\n',
                  'magic,config,physics,total,samples,savedU,savedT,nb,savedSigma\n')
    s = once(s, "if(magic/='MLHRRESTART0001') error stop 'Checkpoint is not the multilevel format'",
        "if(magic/='MLHRRESTART0001') error stop 'Checkpoint is not the HRR multilevel format'\n"
        "        if(.not.ieee_is_finite(savedSigma)) error stop 'Nonfinite checkpoint HRR sigma'\n"
        "        if(savedSigma/=hrrMixingSigma) error stop 'Checkpoint HRR sigma mismatch'")
    write(NEW, s)
    write(OUT / 'grid_parent_to_hrr.diff', ''.join(difflib.unified_diff(
        grid.splitlines(True), s.splitlines(True), fromfile=GRID_PARENT.name, tofile=NEW.name)))
    for name in ['streaming','bounceback','macro','collisionT','streamingT','bouncebackT','macroT']:
        assert proc(s, name) == proc(grid, name), name
    for name in ['collision','regularized_node','regularized_equilibrium']:
        assert proc(s, name) == proc(hrr, name), name
    reference = HRR_PARENT.parents[2]
    for src, dest in [(reference/'source/verify.F90', OUT/'verification/verify.F90'),
            (reference/'verification_inputs/hrr995/input.dat', OUT/'verification/input.dat'),
            (reference/'verification_inputs/hrr995/expected.dat', OUT/'verification/expected.dat')]:
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dest)
    cases=[]
    for exponent in [9,10]:
        name=f'Ra1e{exponent}_N2048_L3_HRR0995'
        data=assign(s,'Rayleigh',f'1.0d{exponent}')
        source=OUT/name/'source/solver.F90'
        write(source,data)
        shutil.copyfile(POST, source.parent/'postprocess.F90')
        cases.append(dict(case=name,remote=REMOTE+'/'+name,Ra=10**exponent,
            equivalent_grids=[2048,1024,512],nGridLevels=3,adjacent_ratio=2,nBlocks=9,
            finest_band_cells=128,middle_band_local_cells=64,
            Prandtl=.7,Mach=.1,hrr_sigma=.995,from_zero=True,target_tff=1000,
            NuRe_interval_tff=.5,snapshot_interval_tff=5,restart_tecplot_interval_tff=100,
            source_sha256=sha(source),postprocess_sha256=sha(source.parent/'postprocess.F90')))
    manifest=dict(name=NAME,remote=REMOTE,new_solver=str(NEW.relative_to(ROOT)),
        new_solver_sha256=sha(NEW),grid_parent=str(GRID_PARENT.relative_to(ROOT)),
        grid_parent_sha256=sha(GRID_PARENT),hrr_parent=str(HRR_PARENT.relative_to(ROOT)),
        hrr_parent_sha256=sha(HRR_PARENT),cases=cases,
        method='H-RR on all flow blocks/levels; unchanged original Luo D2Q5 temperature',
        grid='Static nested fine ring / middle ring / coarse center; nine compact blocks',
        checkpoint='Independent MLHRRESTART0001 header includes sigma',
        notes='Local block gradients inherited from two-level HRR; no new cross-level strain transfer')
    save(OUT/'manifest.json',manifest)
    assert protected == {str(p.relative_to(ROOT)):sha(p) for p in protected_paths}
    save(OUT/'protected_originals_after.json',protected)
    print(json.dumps(dict(new_solver=str(NEW),campaign=str(OUT),cases=cases,
        protected_originals_unchanged=True),ensure_ascii=False))


if __name__ == '__main__':
    main()
