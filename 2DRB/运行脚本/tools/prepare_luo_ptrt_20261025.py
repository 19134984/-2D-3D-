"""Create the isolated Luo/P-TRT solver; protected sources are read only."""
from pathlib import Path
import difflib
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[2]
PARENT = ROOT / 'Xs/D2Q9TRT_D2Q5LuoTRT/2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90'
TARGET = ROOT / 'PTRT/2DRBOpenaccLuoPTRT.F90'
OUT = ROOT / '运行脚本/run_cases/LuoPTRT_Ra9_Ra10_20261025'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write(p, s):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(s, encoding='utf-8', newline='\n')


def once(s, old, new):
    assert s.count(old) == 1, (old[:100], s.count(old))
    return s.replace(old, new)


assert not TARGET.exists(), 'Refuse to replace an existing solver'
assert not OUT.exists(), 'Refuse to replace an existing campaign'
protected_paths = [PARENT, *ROOT.glob('均匀网格/*.F90'), ROOT/'多块网格/2DRBOpenaccMultiblock.F90']
protected = {p.relative_to(ROOT).as_posix(): sha(p) for p in protected_paths}
original = PARENT.read_text(encoding='utf-8-sig')
s = once(original, '!!!    流场：D2Q9-TRT LBM-CDE；温度场：Luo D2Q5-TRT',
         '!!!    Luo 二维独立副本：D2Q9 P-TRT；温度场保持 Luo D2Q5-TRT')
s = once(s, '        real(kind=8), parameter :: chi_nu=0.5d0',
         '        real(kind=8), parameter :: chi_nu=0.0d0')
s = once(s, '        integer(kind=4), parameter :: nx=2048, ny=2048',
         '        integer(kind=4), parameter :: nx=1024, ny=1024')
s = once(s, '        real(kind=8), parameter :: Rayleigh=1.0d10',
         '        real(kind=8), parameter :: Rayleigh=1.0d9')
s = once(s, '        real(kind=8), parameter :: Snu=1.0d0/tauf',
'''        ! 本副本专用：关闭 enablePTRT 可回退原碰撞；P1/P2 只改变 Sq。
        logical, parameter :: enablePTRT=.true.
        integer(kind=4), parameter :: ptrtSqStrategy=1 ! 1=P1 魔术参数；2=P2 固定 1.6
        real(kind=8), parameter :: Snu=1.0d0/tauf
        real(kind=8), parameter :: ptrtSqMagic=8.0d0*(2.0d0-Snu)/(8.0d0-Snu)
        real(kind=8), parameter :: ptrtSq=merge(ptrtSqMagic,1.6d0,ptrtSqStrategy.EQ.1)''')
# Keep the four legacy macro paths available in this copy.
s = re.sub(r'parameter :: Sq=([^\n]+)', r'parameter :: legacySq=\1', s)
s = once(s, '        real(kind=8), parameter :: Qk=3.0d0-dsqrt(3.0d0)',
         '        real(kind=8), parameter :: Sq=merge(ptrtSq,legacySq,enablePTRT)\n\n'
         '        real(kind=8), parameter :: Qk=3.0d0-dsqrt(3.0d0)')
s = once(s, '    ! Initialization\n    call initial()',
'''    if(enablePTRT) then
        if(chi_nu.NE.0.0d0.OR.chi_b.NE.0.0d0) error stop "P-TRT requires chi_nu=chi_b=0"
        if(ptrtSqStrategy.NE.1.AND.ptrtSqStrategy.NE.2) error stop "Invalid P-TRT Sq strategy"
        if(.not.(Snu.GT.0.0d0.AND.Snu.LT.2.0d0.AND.Sq.GT.0.0d0.AND.Sq.LT.2.0d0)) &
            error stop "Invalid P-TRT relaxation rate"
    endif
    ! Initialization
    call initial()''')
s = once(s, '    write(00,*) "thermalScheme = Luo D2Q5-TRT  time-difference correction"',
'''    if(enablePTRT) then
        write(00,*) 'Flow collision = Luo P-TRT; ptrtSqStrategy =',ptrtSqStrategy
        write(00,*) 's_e/s_nu/s_q =',Snu,Snu,Sq
        write(00,*) 'tau_nu/tau_q/Lambda_s =',1.0d0/Snu,1.0d0/Sq, &
            (1.0d0/Snu-0.5d0)*(1.0d0/Sq-0.5d0)
        write(00,*) 'ghost rule: m_post(2)=-m_post(0)-m_post(1); no independent s_epsilon'
        write(00,*) 'Force: one lattice-step momentum density; u=(j+F/2)/rho0'
    endif
    write(00,*) "thermalScheme = unchanged Luo D2Q5-TRT; no time-difference correction"''')
s = once(s, '        call collision()\n', '''        call collision()
        ! 必须在迁移覆盖 f 之前检查本次碰撞；只在原奇矩诊断间隔做额外归约。
        if(enablePTRT.AND.MOD(restartItcOffset+itc,flowOddMomentDiagnosticIntervalItc).EQ.0) &
            call monitor_ptrt_collision()
''')
# Refresh the physical velocity with the force at T^n before equilibrium/strain calculation.
s = once(s, '            densityFluctuation = rho(i,j)-rho0\n', '''            if(enablePTRT) then
                ! f 存 delta_rho；用当前 T^n 的单步力重建速度，不累计旧半步力。
                rho(i,j) = rho0+m(0)
                u(i,j) = (m(3)+0.5d0*Fx(i,j))/rho0
                v(i,j) = (m(5)+0.5d0*Fy(i,j))/rho0
            endif
            densityFluctuation = rho(i,j)-rho0
''')
old = '''            do alpha = 0, 8
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+ &
                    (1.0d0-0.5d0*s(alpha))*mSource(alpha)
            enddo'''
s = once(s, old, '''            if(enablePTRT) then
                ! 新分支直接给出完整 m*；下方旧 mSource 不能再次叠加。
                call ptrt_moment_collision(m,rho0,u(i,j),v(i,j),Fx(i,j),Fy(i,j),Snu,Sq,m_post)
            else
''' + old + '''
            endif''')
# Diagnostics must use the same cubic equilibrium as the new collision.
oldq = '            qxNeq = m4+rho0*u(i,j)\n            qyNeq = m6+rho0*v(i,j)'
assert s.count(oldq) == 2
s = s.replace(oldq, oldq + '''
            if(enablePTRT) then
                qxNeq = qxNeq-3.0d0*rho0*u(i,j)*v(i,j)*v(i,j)
                qyNeq = qyNeq-3.0d0*rho0*v(i,j)*u(i,j)*u(i,j)
            endif''')
# Keep the baseline phase order: flow then thermal. Reconstruct F^(n+1) only after T advances.
s = once(s, '        call macroT()\n', '''        call macroT()
        if(enablePTRT) call refresh_ptrt_velocity()
''')
s = once(s, '    call enter_data_2d_openacc()\n',
         '    call enter_data_2d_openacc()\n    if(enablePTRT) call refresh_ptrt_velocity()\n')
s = once(s, '    firstUnsteadyHistoryWrite = .false.\n',
         '    if(enablePTRT) call monitor_ptrt_wall_flux()\n    firstUnsteadyHistoryWrite = .false.\n')
extra = (ROOT/'运行脚本/tools/luo_ptrt_routines.inc').read_text(encoding='utf-8')
s += '\n' + extra
s = re.sub(r'(?m)^( +)(equilibriumValue = omega\(alpha\).*\n +0\.5d0\*eu\*eu/\(cs2\*cs2\)-0\.5d0\*uu/cs2\)\))',
           r'\1\2\n\1if(enablePTRT) equilibriumValue=equilibriumValue+ &\n\1    4.5d0*omega(alpha)*rho0*(eu*eu*eu-eu*uu)',s)
write(TARGET, s)
OUT.mkdir(parents=True)
write(OUT/'protected_sources.json', json.dumps(protected,ensure_ascii=False,indent=2)+'\n')
write(OUT/'master_to_ptrt.diff', ''.join(difflib.unified_diff(original.splitlines(True),s.splitlines(True),
    fromfile=PARENT.relative_to(ROOT).as_posix(),tofile=TARGET.relative_to(ROOT).as_posix())))
write(OUT/'source_master'/PARENT.name, original)
assert {p.relative_to(ROOT).as_posix(): sha(p) for p in protected_paths} == protected
print(json.dumps(dict(source=str(TARGET),sha256=sha(TARGET),protected_sources_unchanged=True),ensure_ascii=False))
