"""Small, local numerical checks. All variants/builds/results go to a temporary ASCII path.

Run: python verify_multiblock.py
Requires gfortran with OpenACC host support and numpy. No repository solver is modified.
"""
from pathlib import Path
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import time
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = HERE / '2DRBOpenaccMultiblock.F90'
PARENT = ROOT / '均匀网格/2DRBOpenacc.F90'
COMPILER = shutil.which('gfortran') or r'C:\msys64\ucrt64\bin\gfortran.exe'
BUILD = Path(tempfile.mkdtemp(prefix='mb2d-verify-'))
ENV = os.environ.copy()
ENV['PATH'] = str(Path(COMPILER).parent) + os.pathsep + ENV.get('PATH', '')
ENV['ACC_DEVICE_TYPE'] = 'host'
REPORT = {'compiler': COMPILER, 'device': 'OpenACC host', 'build_directory': str(BUILD), 'checks': []}
START = time.monotonic()

def run(args, cwd, timeout=90):
    p = subprocess.run([str(a) for a in args], cwd=cwd, env=ENV, capture_output=True, text=True,
                       encoding='utf-8', errors='replace', timeout=timeout)
    if p.returncode:
        raise RuntimeError(f"Command failed: {args}\n{p.stdout}\n{p.stderr}")
    return p.stdout

def variant(s, *, ratio=2, legacy=False, side=False, steady=False, ha=False, restart=False):
    s = s.replace('refineRatio=2', f'refineRatio={ratio}')
    if legacy:
        s = s.replace('#define EnableUseG\n', '!#define EnableUseG\n')
        s = s.replace('!#define EnableLegacyThermalScheme', '#define EnableLegacyThermalScheme')
    if side:
        for key in ['RayleighBenardCell', 'HorizontalWallsConstT', 'VerticalWallsAdiabatic']:
            s = re.sub(r'^#define '+key+r'\s*$', '!#define '+key, s, flags=re.M)
        for key in ['SideHeatedCell', 'HorizontalWallsAdiabatic', 'VerticalWallsConstT']:
            s = re.sub(r'^!#define '+key+r'\s*$', '#define '+key, s, flags=re.M)
    if steady:
        s = re.sub(r'^#define unsteadyFlow\s*$', '!#define unsteadyFlow', s, flags=re.M)
        s = re.sub(r'^!#define steadyFlow\s*$', '#define steadyFlow', s, flags=re.M)
    if ha:
        s = re.sub(r'^!#define SideHeatedHa\s*$', '#define SideHeatedHa', s, flags=re.M)
    if restart:
        s = s.replace('loadInitField=0', 'loadInitField=1')
    return s

def compile_source(name, src, driver=None, n=96, ra=1000, syntax=False):
    folder = BUILD / name
    folder.mkdir(exist_ok=True)
    if driver:
        src = re.sub(r'^\s*program main\b.*?^\s*end program main\b', driver, src, flags=re.S|re.M|re.I)
    file = folder / 'solver.F90'
    file.write_text(src, encoding='utf-8')
    exe = folder / 'solver.exe'
    flags = ['-cpp', '-fopenacc', '-ffree-line-length-none', '-O1', '-fcheck=all',
             '-fbacktrace', f'-DNX_OVERRIDE={n}', f'-DNY_OVERRIDE={n}', f'-DRAYLEIGH_OVERRIDE={ra}']
    if syntax:
        flags += ['-fsyntax-only']
    else:
        flags += ['-o', exe]
    run([COMPILER, *flags, file], folder)
    return folder, exe

def state_writer(parent=False):
    names=['f','g','u','v','T','rho','Fx','Fy','Bx_prev','By_prev']
    arrays=','.join(names if parent else ['blocks(1)%'+v for v in names])
    return f"""
    open(unit=77,file='state.bin',form='unformatted',access='stream',status='replace')
    write(77) {arrays}
    close(77)
"""

def uniform_driver(parent=False):
    steps='''
        call collision()
        call streaming()
        call bounceback()
        call macro()
        call collisionT()
        call streamingT()
        call bouncebackT()
        call macroT()
''' if parent else '        call advance_multiblock()\n'
    sync='''
    call update_host_snapshot_2d_openacc()
    call update_host_reload_2d_openacc()
    !$acc update self(Fx,Fy)
''' if parent else '    call update_host_all(.true.)\n'
    return f"""program main
    use commondata
    use openacc
    implicit none
    integer :: step
    call acc_init(acc_device_host)
    call initial()
    call enter_data_2d_openacc()
    do step=1,40
{steps}
    enddo
{sync}
{state_writer(parent)}
    call exit_data_2d_openacc()
end program main"""

CONDUCTION = '''program main
    use commondata
    use openacc
    implicit none
    integer :: b,i,j,step,a
    real(8) :: m(0:4),gv(0:4),x,y,exact,err,masserr
    call acc_init(acc_device_host)
    call initial()
    do b=1,nBlocks
        blocks(b)%gb=0.0d0
        do j=1,blocks(b)%nj
            y=blocks(b)%y0+(dble(j)-0.5d0)*blocks(b)%h
            do i=1,blocks(b)%ni
                x=blocks(b)%x0+(dble(i)-0.5d0)*blocks(b)%h
#ifdef SideHeatedCell
                exact=Thot+(Tcold-Thot)*x/dble(nx)
#else
                exact=Thot+(Tcold-Thot)*y/dble(ny)
#endif
                blocks(b)%T(i,j)=exact
                m=[exact,0.0d0,0.0d0,thermalA*exact,0.0d0]
                ! 精确线性导热的迁移前非平衡热流矩；适用于两种原 D2Q5 分支。
#ifdef SideHeatedCell
                m(1)=-(thermalA+4.0d0)/10.0d0*blocks(b)%h/blocks(b)%qk*(Tcold-Thot)/dble(nx)
#else
                m(2)=-(thermalA+4.0d0)/10.0d0*blocks(b)%h/blocks(b)%qk*(Tcold-Thot)/dble(ny)
#endif
                call thermal_populations(m,gv)
                blocks(b)%g(i,j,:)=gv
            enddo
        enddo
    enddo
    call enter_data_2d_openacc()
    do step=1,300
        call advance_multiblock()
    enddo
    call update_host_all(.true.)
    err=0.0d0; masserr=0.0d0
    do b=1,nBlocks
        do j=blocks(b)%jlo,blocks(b)%jhi
            do i=blocks(b)%ilo,blocks(b)%ihi
#ifdef SideHeatedCell
                exact=Thot+(Tcold-Thot)*(blocks(b)%x0+(dble(i)-0.5d0)*blocks(b)%h)/dble(nx)
#else
                exact=Thot+(Tcold-Thot)*(blocks(b)%y0+(dble(j)-0.5d0)*blocks(b)%h)/dble(ny)
#endif
                err=max(err,abs(blocks(b)%T(i,j)-exact))
                masserr=max(masserr,abs(blocks(b)%rho(i,j)-1.0d0))
            enddo
        enddo
    enddo
    call calNuRe()
    write(*,*) 'CONDUCTION_ERROR',err,masserr
    if (err>1.0d-10 .or. masserr>1.0d-10) error stop 'Linear conduction not preserved'
    call exit_data_2d_openacc()
end program main'''

RESTART_DRIVER = '''program main
    use commondata
    use openacc
    implicit none
    integer :: step,b,finish
    character(16) :: arg
    call get_command_argument(1,arg)
    read(arg,*) finish
    call acc_init(acc_device_host)
    call initial()
    call enter_data_2d_openacc()
    do while (itc<finish)
        call advance_multiblock()
    enddo
    call update_host_all(.true.)
    call output_ReloadFile()
    open(unit=77,file='allstate.bin',form='unformatted',access='stream',status='replace')
    do b=1,nBlocks
        write(77) blocks(b)%f,blocks(b)%g,blocks(b)%u,blocks(b)%v,blocks(b)%T,blocks(b)%rho, &
            blocks(b)%Fx,blocks(b)%Fy,blocks(b)%Bx_prev,blocks(b)%By_prev,blocks(b)%p
    enddo
    close(77)
    call exit_data_2d_openacc()
end program main'''

def main():
    source=SOURCE.read_text(encoding='utf-8-sig')
    parent=PARENT.read_text(encoding='utf-8-sig')
    before=hashlib.sha256(PARENT.read_bytes()).hexdigest()
    REPORT['parent_sha256']=before
    REPORT['source_sha256']=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    cases=[('rb_useg',{}),('rb_legacy_steady',{'legacy':True,'steady':True}),
           ('side_useg_steady',{'side':True,'steady':True}),('side_legacy',{'side':True,'legacy':True}),
           ('side_ha',{'side':True,'ha':True})]
    for name,kw in cases:
        compile_source('syntax_'+name,variant(source,**kw),syntax=True)
    REPORT['checks'].append({'syntax_openacc_cases':[name for name,_ in cases]})
    print('OpenACC syntax matrix passed',flush=True)

    for name,kw in [('rb_useg',{}),('side_legacy',{'side':True,'legacy':True})]:
        pdir,pexe=compile_source('parent_'+name,variant(parent,**kw),uniform_driver(True),n=32)
        mdir,mexe=compile_source('uniform_'+name,variant(source,ratio=1,**kw),uniform_driver(),n=32)
        run([pexe],pdir); run([mexe],mdir)
        expected=np.fromfile(pdir/'state.bin',dtype='<f8'); actual=np.fromfile(mdir/'state.bin',dtype='<f8')
        err=float(np.max(np.abs(actual-expected)))
        if not np.all(np.isfinite(actual)) or err>3e-13:
            raise AssertionError(f'Uniform regression {name}: {err}')
        REPORT['checks'].append({'uniform_regression':name,'fine_steps':40,'max_absolute_error':err})
        print('Uniform regression',name,err,flush=True)

    for name,kw in [('rb_useg',{}),('side_useg',{'side':True}),('side_legacy',{'side':True,'legacy':True})]:
        folder,exe=compile_source('conduction_'+name,variant(source,**kw),CONDUCTION,ra=10000 if kw.get('legacy') else 1000)
        stdout=run([exe],folder,timeout=150)
        line=next(line for line in stdout.splitlines() if 'CONDUCTION_ERROR' in line)
        result=[float(x) for x in line.split()[1:]]
        history=np.loadtxt(folder/'NuRe_2DOpenaccMultiblock.dat')
        if max(abs(history[[1,3,4,5]]-1))>1e-10:
            raise AssertionError(f'Conduction Nu diagnostics: {history}')
        REPORT['checks'].append({'conduction':name,'coarse_steps':300,'max_T_error':result[0],
                                 'max_rho_error':result[1],'Nu_volume_hot_cold_middle':history[[1,3,4,5]].tolist()})
        print('Conduction',name,result,flush=True)

    full,fexe=compile_source('restart_full',source,RESTART_DRIVER)
    split,sexe=compile_source('restart_split',source,RESTART_DRIVER)
    resume,rexe=compile_source('restart_resume',variant(source,restart=True),RESTART_DRIVER)
    run([fexe,'40'],full); run([sexe,'20'],split)
    shutil.copy2(rexe,split/'resume.exe')
    run([split/'resume.exe','40'],split)
    expected=np.fromfile(full/'allstate.bin',dtype='<f8'); actual=np.fromfile(split/'allstate.bin',dtype='<f8')
    err=float(np.max(np.abs(actual-expected)))
    if not np.all(np.isfinite(actual)) or err != 0.0:
        raise AssertionError(f'Restart differs from uninterrupted run: {err}')
    REPORT['checks'].append({'restart_exact':True,'fine_steps':40,'split_at':20,'max_absolute_error':err})
    print('Exact multiblock restart passed',flush=True)

    # Exercise the actual program, output clocks, binary snapshots and history-backed restart.
    smoke=source.replace('unsteadyRunDuration=1000.0d0','unsteadyRunDuration=0.1d0')
    smoke=smoke.replace('outputSnapshotInterval=0.5d0','outputSnapshotInterval=0.01d0')
    smoke=smoke.replace('reloadFileInterval=100.0d0','reloadFileInterval=0.03d0')
    smoke=smoke.replace('outputPltFileInterval=100.0d0','outputPltFileInterval=0.05d0')
    whole,wexe=compile_source('main_full',smoke)
    short,texe=compile_source('main_short',smoke.replace('unsteadyRunDuration=0.1d0','unsteadyRunDuration=0.05d0'))
    resumed,rexe=compile_source('main_resumed',variant(smoke,restart=True))
    run([wexe],whole); run([texe],short)
    shutil.copy2(rexe,short/'resume.exe'); run([short/'resume.exe'],short)
    hfull=np.loadtxt(whole/'NuRe_2DOpenaccMultiblock.dat')
    hsplit=np.loadtxt(short/'NuRe_2DOpenaccMultiblock.dat')
    if not np.array_equal(hfull,hsplit):
        raise AssertionError('Actual main program restart changes history')
    stats=(whole/'NuReStatistics_2DOpenaccMultiblock.dat').read_text()
    if 'INCOMPLETE' in stats or 'relative half-window difference:' not in stats:
        raise AssertionError('Main statistics window is not covered')
    REPORT['checks'].append({'main_smoke_and_history_restart':True,'samples':len(hfull),'target_t_ff':0.1})
    print('Actual main program, outputs and restart passed',flush=True)
    if hashlib.sha256(PARENT.read_bytes()).hexdigest()!=before:
        raise AssertionError('Parent source changed')
    REPORT['elapsed_seconds']=time.monotonic()-START
    REPORT['status']='passed'
    result=HERE/'verification_results.json'
    result.write_text(json.dumps(REPORT,indent=2),encoding='utf-8')
    print('Report:',result,flush=True)

if __name__=='__main__':
    try:
        main()
    except Exception:
        REPORT['status']='failed'
        (BUILD/'failure_report.json').write_text(json.dumps(REPORT,indent=2),encoding='utf-8')
        print('Build evidence:',BUILD,flush=True)
        raise
