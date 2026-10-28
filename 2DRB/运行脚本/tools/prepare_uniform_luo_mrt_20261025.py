"""Prepare isolated uniform-grid Luo MRT tests; never edit protected solvers."""
from pathlib import Path
import difflib
import hashlib
import json
import math
import tarfile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / '运行脚本/run_cases/Uniform_LuoMRT_Ra9_Ra10_T1000_20261025'
REMOTE = '/data2/XLLi/LBMCDE/Uniform_LuoMRT_Ra9_Ra10_T1000_20261025'
MASTER = ROOT / 'Xs/D2Q9TRT_D2Q5LuoTRT/2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90'
PROTECTED = [MASTER, ROOT/'多块网格/2DRBOpenaccMultiblock.F90',
             ROOT/'均匀网格/2DRBOpenmp.F90', ROOT/'均匀网格/2DRBOpenacc.F90']
raw = MASTER.read_bytes()
original = raw.decode('utf-8-sig').replace('\r\n', '\n')
assert not OUT.exists(), 'Refuse to replace an existing campaign'
OUT.mkdir(parents=True)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding='utf-8', newline='\n')


def once(text, old, new):
    assert text.count(old) == 1, (old, text.count(old))
    return text.replace(old, new)


protected = {str(p.relative_to(ROOT)): sha(p.read_bytes()) for p in PROTECTED}
(OUT/'source_master').mkdir()
(OUT/'source_master'/MASTER.name).write_bytes(raw)
base = once(original, '!!!    流场：D2Q9-TRT LBM-CDE；温度场：Luo D2Q5-TRT',
            '!!!    均匀网格独立测试：流场 Luo D2Q9-MRT；温度场 Luo D2Q5-TRT')
base = once(base, 'real(kind=8), parameter :: chi_nu=0.5d0',
            'real(kind=8), parameter :: chi_nu=0.0d0')
base = once(base, '        real(kind=8), parameter :: bulkViscosity=(tauf-0.5d0)*(1.0d0-chi_b)*cs2',
            '        ! Contrino/Luo 2014 §5.3：独立设置流场能量矩和四阶矩，不再是 TRT。\n'
            '        real(kind=8), parameter :: flowEnergyRate=1.2d0  ! s_e；不是温度 Qe\n'
            '        real(kind=8), parameter :: flowEpsilonRate=1.3d0 ! s_epsilon；四阶矩\n'
            '        ! chi_nu=chi_b=0：无 CDE 修正；体黏性诊断随 s_e 更新，剪切黏性不变。\n'
            '        real(kind=8), parameter :: bulkViscosity=(1.0d0/flowEnergyRate-0.5d0)*cs2')
base = once(base, '            s(1) = Snu', '            s(1) = flowEnergyRate')
base = once(base, '            s(2) = Snu', '            s(2) = flowEpsilonRate')
base = once(base, '    write(00,*) "Flow odd magic policy = original tau_0 scale"',
            '    write(00,*) "Flow collision = Luo MRT; sq follows paper Eq. (7b)"\n'
            "    write(00,*) 'flow s_e / s_epsilon =',flowEnergyRate,flowEpsilonRate\n"
            "    write(00,*) 'shear/odd product =',(1.0d0/Snu-0.5d0)*(1.0d0/Sq-0.5d0)\n"
            '    ! 保留 3/16 的参数关系不等于证明 MRT 具有全部 TRT 边界性质。')
base = once(base, '    write(00,*) "thermalScheme = Luo D2Q5-TRT  time-difference correction"',
            '    write(00,*) "thermalScheme = Luo D2Q5-TRT; no time-difference correction"')
base = once(base, '    ! Initialization\n    call initial()',
            '    ! 本诊断副本只允许无 CDE 修正的论文 MRT 配置。\n'
            '    if(chi_nu.NE.0.0d0.OR.chi_b.NE.0.0d0) error stop "MRT test requires chi=0"\n'
            '    ! Initialization\n    call initial()')
assert 's(7) = Snu' in base and 's(8) = Snu' in base
assert '(1.0d0-0.5d0*s(alpha))*mSource(alpha)' in base
assert '#define RayleighBenardCell\n' in base and '#define unsteadyFlow\n' in base
assert 'loadInitField=0' in base

configs = {}
for name, n, ra, nu_ref, re_ref in [('Ra1e9',1025,1e9,53.51,15101),
                                     ('Ra1e10',2049,1e10,93.88,61797)]:
    text = once(base, 'nx=2048, ny=2048', f'nx={n}, ny={n}')
    if ra != 1e10:
        text = once(text, 'Rayleigh=1.0d10', 'Rayleigh=1.0d9')
    write(OUT/name/'source/solver.F90', text)
    (OUT/name/'results').mkdir()
    write(OUT/name/'source/master_to_case.diff', ''.join(difflib.unified_diff(
        original.splitlines(True),text.splitlines(True),fromfile=str(MASTER.relative_to(ROOT)),
        tofile=name+'/source/solver.F90')))
    nu = .1*n*math.sqrt(.7/(3*ra))
    tau = .5+3*nu
    sq = 1/(.5+(3/16)/(tau-.5))
    configs[name] = dict(nx=n,ny=n,Ra=ra,Pr=.7,Ma=.1,chi_nu=0,chi_b=0,
        s_e=1.2,s_epsilon=1.3,s_nu=1/tau,s_q=sq,shear_odd_product=3/16,
        viscosity=nu,diffusivity=nu/.7,bulk_viscosity=(1/1.2-.5)/3,
        target_tff=1000,history_interval_tff=1,nonfinite_interval_tff=1,
        odd_diagnostic_interval_tff=.1,from_scratch=True,
        initial_condition='unchanged quiescent conductive profile; no added perturbation',
        thermal='unchanged Luo D2Q5 TRT',Nu_reference=nu_ref,Re_reference=re_ref,
        source_sha256=sha(text.encode('utf-8')))

# Small OpenACC-host smoke case checks execution only; it is not a high-Ra result.
smoke = once(base, 'nx=2048, ny=2048', 'nx=32, ny=32')
smoke = once(smoke, 'Rayleigh=1.0d10', 'Rayleigh=1.0d5')
smoke = once(smoke, 'unsteadyRunDuration=1000.0d0', 'unsteadyRunDuration=2.0d0')
write(OUT/'smoke/source/solver.F90',smoke)
(OUT/'smoke/results').mkdir()
write(OUT/'audit_flow_history.py',(ROOT/'后处理/audit_flow_history.py').read_text(encoding='utf-8-sig'))
write(OUT/'campaign.json',json.dumps(dict(remote_dir=REMOTE,master=str(MASTER),
    master_sha256=sha(raw),node='node07',order=['Ra1e9','Ra1e10'],
    reference='Zhang et al. JFM 814 (2017), Table 1; DOI 10.1017/jfm.2017.19',
    mrt_reference='Contrino et al. JCP 275 (2014), §5.3, p.267; DOI 10.1016/j.jcp.2014.06.047',
    cases=configs),ensure_ascii=False,indent=2)+'\n')
write(OUT/'protected_sources.json',json.dumps(protected,ensure_ascii=False,indent=2)+'\n')
assert {str(p.relative_to(ROOT)): sha(p.read_bytes()) for p in PROTECTED} == protected
print(json.dumps(dict(local=str(OUT),remote=REMOTE,cases=configs),ensure_ascii=False,indent=2))
