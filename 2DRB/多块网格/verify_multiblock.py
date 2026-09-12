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

def compile_source(name, src, driver=None, n=96, ra=1000, syntax=False, ny=None):
    folder = BUILD / name
    folder.mkdir(exist_ok=True)
    if driver:
        src = re.sub(r'^\s*program main\b.*?^\s*end program main\b', driver, src, flags=re.S|re.M|re.I)
    file = folder / 'solver.F90'
    file.write_text(src, encoding='utf-8')
    exe = folder / 'solver.exe'
    flags = ['-cpp', '-fopenacc', '-ffree-line-length-none', '-O1', '-fcheck=all',
             '-fbacktrace', f'-DNX_OVERRIDE={n}', f'-DNY_OVERRIDE={n if ny is None else ny}', f'-DRAYLEIGH_OVERRIDE={ra}']
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

GEOMETRY_DRIVER = '''program main
    use commondata
    implicit none
    integer :: b,i,j,l,directCoarse,directFine,interpolatedFine
    real(8) :: x,y,area,xmoment,ymoment,w
    call initial()
    area=0.0d0; xmoment=0.0d0; ymoment=0.0d0
    do b=1,nBlocks
        do j=1,blocks(b)%nj
            y=blocks(b)%y0+(dble(j)-0.5d0)*blocks(b)%h
            do i=1,blocks(b)%ni
                x=blocks(b)%x0+(dble(i)-0.5d0)*blocks(b)%h
                if (abs(x-0.5d0-dble(nint(x-0.5d0)))>1.0d-12 .or. &
                    abs(y-0.5d0-dble(nint(y-0.5d0)))>1.0d-12) error stop 'Node is off the fine lattice'
                w=blocks(b)%dxWeight(i)*blocks(b)%dyWeight(j)
                if (w<0.0d0) error stop 'Negative integration weight'
                area=area+w; xmoment=xmoment+w*x; ymoment=ymoment+w*y
            enddo
        enddo
        if (blocks(b)%wall(1) .and. blocks(b)%x0+0.5d0*blocks(b)%h/=0.5d0) error stop 'Left wall moved'
        if (blocks(b)%wall(2) .and. &
            blocks(b)%x0+(blocks(b)%ni-0.5d0)*blocks(b)%h/=nx-0.5d0) error stop 'Right wall moved'
        if (blocks(b)%wall(3) .and. blocks(b)%y0+0.5d0*blocks(b)%h/=0.5d0) error stop 'Bottom wall moved'
        if (blocks(b)%wall(4) .and. &
            blocks(b)%y0+(blocks(b)%nj-0.5d0)*blocks(b)%h/=ny-0.5d0) error stop 'Top wall moved'
    enddo
    if (abs(area-dble(nx)*ny)>1.0d-9) error stop 'Area is double counted or missing'
    if (abs(xmoment-0.5d0*dble(nx)**2*ny)>1.0d-8) error stop 'Wrong first x moment'
    if (abs(ymoment-0.5d0*dble(ny)**2*nx)>1.0d-8) error stop 'Wrong first y moment'
    directCoarse=0; directFine=0; interpolatedFine=0
    do l=1,nLinks
        if (links(l)%receiver==1) then
            if (.not.all(links(l)%coincident)) error stop 'Coarse receiver should coincide with fine nodes'
            directCoarse=directCoarse+links(l)%count
        elseif (links(l)%donor==1) then
            directFine=directFine+count(links(l)%coincident)
            interpolatedFine=interpolatedFine+count(.not.links(l)%coincident)
        endif
    enddo
    if (min(directCoarse,directFine,interpolatedFine)<=0) error stop 'Missing direct or interpolation interface path'
    write(*,*) 'ALIGNED_GEOMETRY',directCoarse,directFine,interpolatedFine,area
end program main'''

PACKET_DRIVER = '''program main
    use commondata, only: packetSize,lagrange_weights,interpolate_packets
    use openacc
    implicit none
    real(8) :: p(8,8,packetSize,0:2),wx(4,2),wy(4,2),val(packetSize,2),wt(0:2),err,expected
    integer :: si(2),sj(2),i,j,k,a
    logical :: same(2)
    call acc_init(acc_device_host)
    do k=0,2
        do a=1,packetSize
            do j=1,8
                do i=1,8
                    p(i,j,a,k)=a*(dble(i)**2+dble(j)**3)*(dble(k-1)**2+2.0d0*(k-1)+4.0d0)
                enddo
            enddo
        enddo
    enddo
    si=[8,2]; sj=[8,3]; same=[.true.,.false.]
    wx(:,1)=[1.0d0,0.0d0,0.0d0,0.0d0]; wy(:,1)=wx(:,1)
    call lagrange_weights(1.5d0,wx(:,2)); call lagrange_weights(1.5d0,wy(:,2))
    wt=[-0.125d0,0.75d0,0.375d0]
    !$acc enter data copyin(p,si,sj,wx,wy,same) create(val)
    call interpolate_packets(8,8,2,p,2,si,sj,wx,wy,same,wt,val)
    !$acc wait(1)
    !$acc update self(val)
    err=0.0d0
    do a=1,packetSize
        expected=a*(8.0d0**2+8.0d0**3)*5.25d0
        err=max(err,abs(val(a,1)-expected))
        expected=a*(3.5d0**2+4.5d0**3)*5.25d0
        err=max(err,abs(val(a,2)-expected))
    enddo
    if (err>1.0d-10) error stop 'Direct/cubic-space/quadratic-time packet transfer failed'
    !$acc exit data delete(p,si,sj,wx,wy,same,val)
    write(*,*) 'PACKET_ERROR',err
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

    for nx,ny in [(96,96),(128,96)]:
        folder,exe=compile_source(f'geometry_{nx}_{ny}',source,GEOMETRY_DRIVER,n=nx,ny=ny)
        stdout=run([exe],folder)
        values=next(line for line in stdout.splitlines() if 'ALIGNED_GEOMETRY' in line).split()[1:]
        REPORT['checks'].append({'aligned_geometry':[nx,ny],'coarse_direct_nodes':int(values[0]),
            'fine_direct_from_coarse':int(values[1]),'fine_interpolated_from_coarse':int(values[2]),
            'area':float(values[3]),'linear_integrals_and_physical_walls':'passed'})
    folder,exe=compile_source('packet_transfer',source,PACKET_DRIVER)
    stdout=run([exe],folder)
    err=float(next(line for line in stdout.splitlines() if 'PACKET_ERROR' in line).split()[1])
    REPORT['checks'].append({'direct_and_interpolated_packet_transfer_error':err})
    print('Node alignment, quadrature and interface transfer passed',flush=True)

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
    # Old, staggered-grid snapshots must fail before reading state arrays.
    latest=(split/'reloadFile2DOpenaccMultiblock-latest.meta').read_text().strip()
    old=split/'old-layout.bin'
    state=(split/latest).read_bytes()
    old.write_bytes(b'MB2DRESTART0001  '[:16]+state[16:])
    (split/'reloadFile2DOpenaccMultiblock-latest.meta').write_text('old-layout.bin\n')
    bad=subprocess.run([str(split/'resume.exe'),'40'],cwd=split,env=ENV,capture_output=True)
    if bad.returncode==0 or b'Wrong checkpoint format' not in bad.stderr:
        raise AssertionError('Old staggered checkpoint was not rejected')
    REPORT['checks'].append({'old_layout_restart_rejected':True})

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
    with next(whole.glob('*Snapshot-*.bin')).open('rb') as f:
        assert f.read(16)==b'MB2DSNAPSHOT0002'
        nb,nx,ny,_=np.fromfile(f,dtype='<i4',count=4)
        np.fromfile(f,dtype='<f8',count=2)
        area=0.0
        for b in range(nb):
            ni,nj=np.fromfile(f,dtype='<i4',count=2)
            first_x,first_y,h,*box=np.fromfile(f,dtype='<f8',count=7)
            dx=np.fromfile(f,dtype='<f8',count=ni); dy=np.fromfile(f,dtype='<f8',count=nj)
            assert np.all(np.mod(first_x+np.arange(ni)*h-.5,1)==0)
            assert np.all(np.mod(first_y+np.arange(nj)*h-.5,1)==0)
            area+=dx.sum()*dy.sum()
            assert np.fromfile(f,dtype='<f8',count=4*ni*nj).size==4*ni*nj
        assert area==nx*ny and not f.read(1)
    REPORT['checks'].append({'aligned_snapshot_v2_area_and_coordinates':'passed'})
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
