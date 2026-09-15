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
    s = re.sub(r'\brefineRatio\s*=\s*\d+', f'refineRatio={ratio}', s)
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
        s = s.replace('loadInitField = 0', 'loadInitField = 1')
    return s


def manual_reload_number(src, number):
    # Change the user setting only; preserve the runtime reset for a fresh run.
    src, count = re.subn(r'(::\s*reloadFileNum\s*=\s*)\d+',
                         lambda match: match[1] + str(number), src)
    if count != 1:
        raise AssertionError('Manual reload file number declaration not found')
    return src

def compile_source(name, src, driver=None, n=96, ra=1000, syntax=False, ny=None, walls=None):
    folder = BUILD / name
    folder.mkdir(exist_ok=True)
    if driver:
        src = re.sub(r'^\s*program main\b.*?^\s*end program main\b', driver, src, flags=re.S|re.M|re.I)
    # Multiblock parameters are edited only in this temporary source; the parent still uses -D overrides.
    if '#define NX_OVERRIDE' not in src:
        height = n if ny is None else ny
        if walls is None:
            left, bottom = n // 8, height // 8
            left, right, bottom, top = left, left+1, bottom, bottom+1
        elif len(walls) == 2:
            left, bottom = walls
            left, right, bottom, top = left, left+1, bottom, bottom+1
        else:
            left, right, bottom, top = walls
        src, count = re.subn(r'\bparameter :: nx\s*=\s*\d+, ny\s*=\s*\d+',
                            f'parameter :: nx={n}, ny={height}', src)
        if count != 1:
            raise AssertionError('Manual grid parameters not found')
        for side, value in [('Left',left),('Right',right),('Bottom',bottom),('Top',top)]:
            src, count = re.subn(r'\bparameter :: fineLayerCells'+side+r'\s*=\s*\d+',
                                f'parameter :: fineLayerCells{side}={value}', src)
            if count != 1:
                raise AssertionError(f'Manual {side} wall-layer parameter not found')
        ra_literal = format(float(ra), '.16e').replace('e', 'd')
        src, count = re.subn(r'\bparameter :: Rayleigh\s*=\s*[^,\n]+',
                            f'parameter :: Rayleigh={ra_literal}', src)
        if count != 1:
            raise AssertionError('Manual Rayleigh parameter not found')
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
    arrays=','.join(names)
    select='' if parent else '    call select_block(1)\n'
    return f"""
{select}    open(unit=77,file='state.bin',form='unformatted',access='stream',status='replace')
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
        call select_block(b)
        blockGb(b)=0.0d0
        do j=1,blockNj(b)
            y=blockY0(b)+(dble(j)-0.5d0)*blockH(b)
            do i=1,blockNi(b)
                x=blockX0(b)+(dble(i)-0.5d0)*blockH(b)
#ifdef SideHeatedCell
                exact=Thot+(Tcold-Thot)*x/dble(nx)
#else
                exact=Thot+(Tcold-Thot)*y/dble(ny)
#endif
                T(i,j)=exact
                m=[exact,0.0d0,0.0d0,thermalA*exact,0.0d0]
                ! 精确线性导热的迁移前非平衡热流矩；适用于两种原 D2Q5 分支。
#ifdef SideHeatedCell
                m(1)=-(thermalA+4.0d0)/10.0d0*blockH(b)/blockQk(b)*(Tcold-Thot)/dble(nx)
#else
                m(2)=-(thermalA+4.0d0)/10.0d0*blockH(b)/blockQk(b)*(Tcold-Thot)/dble(ny)
#endif
                call thermal_populations(m,gv)
                g(i,j,:)=gv
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
        call select_block(b)
        do j=blockJlo(b),blockJhi(b)
            do i=blockIlo(b),blockIhi(b)
#ifdef SideHeatedCell
                exact=Thot+(Tcold-Thot)*(blockX0(b)+(dble(i)-0.5d0)*blockH(b))/dble(nx)
#else
                exact=Thot+(Tcold-Thot)*(blockY0(b)+(dble(j)-0.5d0)*blockH(b))/dble(ny)
#endif
                err=max(err,abs(T(i,j)-exact))
                masserr=max(masserr,abs(rho(i,j)-1.0d0))
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
    integer :: b,c,i,j,l,directCoarse,directFine,interpolatedFine
    real(8) :: x,y,area,xmoment,ymoment,w,extraX,extraY,overlap,xLeft,xRight,yBottom,yTop,intersection
    call initial()
    if (blockH(1)/=dble(refineRatio)) error stop 'Wrong coarse spacing'
    xLeft=dble(fineLayerCellsLeft)-0.5d0; xRight=dble(nx-fineLayerCellsRight)+0.5d0
    yBottom=dble(fineLayerCellsBottom)-0.5d0; yTop=dble(ny-fineLayerCellsTop)+0.5d0
    if (any(blockOwnedBox(:, 1)/=blockBaseBox(:, 1))) error stop 'Coarse integration must end at base nodes'
    if (any(blockBaseBox(:, 2)/=[0.0d0,dble(nx),0.0d0,yBottom]) .or. &
        any(blockBaseBox(:, 3)/=[0.0d0,dble(nx),yTop,dble(ny)]) .or. &
        any(blockBaseBox(:, 4)/=[0.0d0,xLeft,yBottom,yTop]) .or. &
        any(blockBaseBox(:, 5)/=[xRight,dble(nx),yBottom,yTop])) error stop 'Fine bases must use requested indices'
    if (blockX0(1)+0.5d0*blockH(1)/=dble(fineLayerCellsLeft)-0.5d0-overlapCells*refineRatio) &
        error stop 'Wrong coarse first node relative to interface'
    if (blockX0(4)+(blockNi(4)-0.5d0)*blockH(4)/= &
        dble(fineLayerCellsLeft)-0.5d0+overlapCells*refineRatio) error stop 'Wrong fine last node'
    overlap=dble(overlapCells*refineRatio)
    extraX=blockBaseBox(2, 1)-xRight
    extraY=blockBaseBox(4, 1)-yTop
    if (extraX/=0.0d0 .or. extraY/=0.0d0) &
        error stop 'Requested interfaces must not move'
    if (blockBaseBox(1, 1)/=xLeft .or. blockBaseBox(3, 1)/=yBottom) error stop 'Coarse anchor moved'
    if (modulo(blockBaseBox(2, 1)-blockBaseBox(1, 1),blockH(1))/=0.0d0 .or. &
        modulo(blockBaseBox(4, 1)-blockBaseBox(3, 1),blockH(1))/=0.0d0) &
        error stop 'Coarse base span is not an integer number of coarse spacings'
    x=blockX0(1)+(blockNi(1)-0.5d0)*blockH(1)
    y=blockY0(1)+(blockNj(1)-0.5d0)*blockH(1)
    if (x/=blockBaseBox(2, 1)+overlap .or. y/=blockBaseBox(4, 1)+overlap) &
        error stop 'Overlap must extend from the aligned coarse base'
    if (blockX0(5)+0.5d0*blockH(5)/=xRight-overlap .or. &
        blockY0(3)+0.5d0*blockH(3)/=yTop-overlap) &
        error stop 'Right/top fine blocks moved with the coarse base'
    if (blockOwnedBox(1, 5)/=blockOwnedBox(2, 1) .or. &
        blockOwnedBox(3, 3)/=blockOwnedBox(4, 1) .or. &
        blockOwnedBox(4, 4)/=blockOwnedBox(3, 3) .or. &
        blockOwnedBox(4, 5)/=blockOwnedBox(3, 3)) error stop 'Statistics split is inconsistent at corners'
    call select_block(1)
    if (blockX0(1)+(blockIlo(1)-0.5d0)*blockH(1)/=blockOwnedBox(1, 1) .or. &
        blockX0(1)+(blockIhi(1)-0.5d0)*blockH(1)/=blockOwnedBox(2, 1) .or. &
        blockY0(1)+(blockJlo(1)-0.5d0)*blockH(1)/=blockOwnedBox(3, 1) .or. &
        blockY0(1)+(blockJhi(1)-0.5d0)*blockH(1)/=blockOwnedBox(4, 1)) error stop 'Integration endpoint is not a coarse node'
    if (any(dxWeight([blockIlo(1),blockIhi(1)])/=blockH(1)/2) .or. &
        any(dyWeight([blockJlo(1),blockJhi(1)])/=blockH(1)/2) .or. &
        any(dxWeight(blockIlo(1)+1:blockIhi(1)-1)/=blockH(1)) .or. &
        any(dyWeight(blockJlo(1)+1:blockJhi(1)-1)/=blockH(1))) error stop 'Coarse weights must be composite trapezoidal'
    if (nx==1024 .and. ny==1024 .and. refineRatio==2 .and. &
        fineLayerCellsLeft==128 .and. fineLayerCellsRight==129 .and. &
        fineLayerCellsBottom==128 .and. fineLayerCellsTop==129 .and. overlapCells==2) then
        if (any(blockBaseBox(:, 1)/=[127.5d0,895.5d0,127.5d0,895.5d0])) error stop 'Wrong default coarse base'
        if (x/=899.5d0 .or. y/=899.5d0) error stop 'Wrong default coarse end nodes'
        if (blockNi(1)/=389 .or. blockNj(1)/=389) error stop 'Wrong default coarse dimensions'
    endif
    area=0.0d0; xmoment=0.0d0; ymoment=0.0d0
    do b=1,nBlocks
        call select_block(b)
        do c=b+1,nBlocks
            intersection=max(0.0d0,min(blockOwnedBox(2, b),blockOwnedBox(2, c))- &
                max(blockOwnedBox(1, b),blockOwnedBox(1, c)))* &
                max(0.0d0,min(blockOwnedBox(4, b),blockOwnedBox(4, c))- &
                max(blockOwnedBox(3, b),blockOwnedBox(3, c)))
            if (intersection>0.0d0) error stop 'Statistics rectangles overlap'
        enddo
        do j=1,blockNj(b)
            y=blockY0(b)+(dble(j)-0.5d0)*blockH(b)
            do i=1,blockNi(b)
                x=blockX0(b)+(dble(i)-0.5d0)*blockH(b)
                if (abs(x-0.5d0-dble(nint(x-0.5d0)))>1.0d-12 .or. &
                    abs(y-0.5d0-dble(nint(y-0.5d0)))>1.0d-12) error stop 'Node is off the fine lattice'
                w=dxWeight(i)*dyWeight(j)
                if (w<0.0d0) error stop 'Negative integration weight'
                area=area+w; xmoment=xmoment+w*x; ymoment=ymoment+w*y
            enddo
        enddo
        if (abs(blockH(b)*(1.0d0/blockSn(b)-0.5d0)-(1.0d0/Snu-0.5d0))>1.0d-13 .or. &
            abs(blockH(b)*(1.0d0/blockSq(b)-0.5d0)-(1.0d0/Sq-0.5d0))>1.0d-13 .or. &
            abs(blockH(b)*(1.0d0/blockQk(b)-0.5d0)-(1.0d0/Qk-0.5d0))>1.0d-13 .or. &
            abs(blockH(b)*(1.0d0/blockQn(b)-0.5d0)-(1.0d0/Qnu-0.5d0))>1.0d-13) &
            error stop 'Relaxation scaling changed physical transport coefficients'
        if (blockGb(b)/=blockH(b)*gBeta) error stop 'Wrong force scaling'
        if (blockWall(1, b) .and. blockX0(b)+0.5d0*blockH(b)/=0.5d0) error stop 'Left wall moved'
        if (blockWall(2, b) .and. &
            blockX0(b)+(blockNi(b)-0.5d0)*blockH(b)/=nx-0.5d0) error stop 'Right wall moved'
        if (blockWall(3, b) .and. blockY0(b)+0.5d0*blockH(b)/=0.5d0) error stop 'Bottom wall moved'
        if (blockWall(4, b) .and. &
            blockY0(b)+(blockNj(b)-0.5d0)*blockH(b)/=ny-0.5d0) error stop 'Top wall moved'
    enddo
    if (abs(area-dble(nx)*ny)>1.0d-9) error stop 'Area is double counted or missing'
    if (abs(xmoment-0.5d0*dble(nx)**2*ny)>1.0d-8) error stop 'Wrong first x moment'
    if (abs(ymoment-0.5d0*dble(ny)**2*nx)>1.0d-8) error stop 'Wrong first y moment'
    directCoarse=0; directFine=0; interpolatedFine=0
    do l=1,nLinks
        call select_link(l)
        if (linkReceiver(l)==1) then
            if (.not.all(linkSame)) error stop 'Coarse receiver should coincide with fine nodes'
            directCoarse=directCoarse+linkCount(l)
        elseif (linkDonor(l)==1) then
            directFine=directFine+count(linkSame)
            interpolatedFine=interpolatedFine+count(.not.linkSame)
        endif
    enddo
    if (min(directCoarse,directFine,interpolatedFine)<=0) error stop 'Missing direct or interpolation interface path'
    write(*,*) 'ALIGNED_GEOMETRY',directCoarse,directFine,interpolatedFine,area
end program main'''

PACKET_DRIVER = '''program main
    use commondata, only: packetSize,itc,refineRatio
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
    itc=0
    call coarse_time_weights(0.5d0,wt)
    if (any(wt/=[0.0d0,0.5d0,0.5d0])) error stop 'Wrong startup time weights'
    itc=refineRatio
    call coarse_time_weights(0.0d0,wt)
    if (any(wt/=[0.0d0,1.0d0,0.0d0])) error stop 'First fine collision must use current time'
    call coarse_time_weights(0.5d0,wt)
    if (any(wt/=[-0.125d0,0.75d0,0.375d0])) error stop 'Wrong midpoint time weights'
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
        call select_block(b)
        write(77) f,g,u,v,T,rho, &
            Fx,Fy,Bx_prev,By_prev,p
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

    for nx,ny in [(96,96),(128,96),(1024,1024)]:
        folder,exe=compile_source(f'geometry_{nx}_{ny}',source,GEOMETRY_DRIVER,n=nx,ny=ny)
        stdout=run([exe],folder)
        values=next(line for line in stdout.splitlines() if 'ALIGNED_GEOMETRY' in line).split()[1:]
        REPORT['checks'].append({'aligned_geometry':[nx,ny],'coarse_direct_nodes':int(values[0]),
            'fine_direct_from_coarse':int(values[1]),'fine_interpolated_from_coarse':int(values[2]),
            'area':float(values[3]),'linear_integrals_and_physical_walls':'passed',
            'four_shared_interfaces_without_rounding':'passed','fine_interfaces_unchanged':'passed',
            'statistics_at_coincident_nodes_trapezoidal_weights_and_disjoint_rectangles':'passed'})
    folder,exe=compile_source('packet_transfer',source,PACKET_DRIVER)
    stdout=run([exe],folder)
    err=float(next(line for line in stdout.splitlines() if 'PACKET_ERROR' in line).split()[1])
    REPORT['checks'].append({'direct_and_interpolated_packet_transfer_error':err,
                            'pre_collision_startup_and_substep_time_weights':'passed'})
    print('Node alignment, quadrature and interface transfer passed',flush=True)

    # Include an odd ratio and rectangular domains, retaining the supported ratio paths.
    for ratio in (3,4,8):
        dims = dict(n=24*ratio, ny=28*ratio, walls=(5*ratio,6*ratio))
        extended = variant(source,ratio=ratio)
        folder,exe=compile_source(f'geometry_r{ratio}',extended,GEOMETRY_DRIVER,**dims)
        run([exe],folder)
        for name,kw in [('rb_useg',{}),('side_useg',{'side':True}),
                        ('side_legacy',{'side':True,'legacy':True})]:
            driver=CONDUCTION.replace('step=1,300','step=1,60').replace(
                '    call update_host_all(.true.)',
                '    if (itc/=60*refineRatio) error stop "Wrong synchronized time"\n'
                '    call update_host_all(.true.)')
            folder,exe=compile_source(f'conduction_r{ratio}_{name}',variant(source,ratio=ratio,**kw),
                                      driver,ra=10000,**dims)
            stdout=run([exe],folder,timeout=150)
            result=[float(x) for x in next(line for line in stdout.splitlines()
                                          if 'CONDUCTION_ERROR' in line).split()[1:]]
            history=np.loadtxt(folder/'NuRe_2DOpenaccMultiblock.dat')
            if not np.all(np.isfinite(history)) or max(abs(history[[1,3,4,5]]-1))>1e-10:
                raise AssertionError(f'Ratio {ratio} conduction diagnostics: {history}')
            REPORT['checks'].append({'ratio':ratio,'conduction':name,'coarse_steps':60,
                'fine_steps':60*ratio,'max_T_error':result[0],'max_rho_error':result[1],
                'geometry_and_transport_scaling':'passed',
                'Nu_volume_hot_cold_middle':history[[1,3,4,5]].tolist()})
        # Exercise moving-flow restart at an intermediate synchronized time for each new ratio.
        full,fexe=compile_source(f'restart_r{ratio}_full',extended,RESTART_DRIVER,**dims)
        split,sexe=compile_source(f'restart_r{ratio}_split',extended,RESTART_DRIVER,**dims)
        resume,rexe=compile_source(f'restart_r{ratio}_resume',variant(source,ratio=ratio,restart=True),
                                  RESTART_DRIVER,**dims)
        run([fexe,str(12*ratio)],full); run([sexe,str(6*ratio)],split)
        shutil.copy2(rexe,split/'resume.exe')
        run([split/'resume.exe',str(12*ratio)],split)
        expected=np.fromfile(full/'allstate.bin',dtype='<f8')
        actual=np.fromfile(split/'allstate.bin',dtype='<f8')
        if not np.all(np.isfinite(actual)) or not np.array_equal(actual,expected):
            raise AssertionError(f'Ratio {ratio} restart mismatch or nonfinite state')
        REPORT['checks'].append({'ratio':ratio,'restart_exact':True,'fine_steps':12*ratio})
        print('Extended ratio',ratio,'geometry, conduction and exact restart passed',flush=True)

    for name,ratio,dims,message in [
        ('negative_ratio',-2,{},'refineRatio must be a positive integer'),
        ('unaligned_size',4,dict(n=98,walls=(24,24)),'must be multiples of refineRatio'),
        ('unaligned_height',4,dict(ny=98,walls=(24,24)),'must be multiples of refineRatio'),
        ('thin_wall',4,dict(walls=(8,24)),'Refined wall layer must exceed'),
        ('empty_core',3,dict(n=63,walls=(32,20)),'coarse core must have positive'),
        ('negative_core',4,dict(n=64,walls=(20,36)),'coarse core must have positive'),
        ('unaligned_core_width',4,dict(walls=(25,26,24,25)),'Central width and height must be multiples'),
        ('unaligned_core_height',4,dict(walls=(24,25,25,26)),'Central width and height must be multiples'),
        ('thin_right',2,dict(walls=(12,1,12,13)),'Refined wall layer must exceed'),
        ('thin_top',2,dict(walls=(12,13,12,1)),'Refined wall layer must exceed')]:
        folder,exe=compile_source('invalid_'+name,variant(source,ratio=ratio),GEOMETRY_DRIVER,**dims)
        bad=subprocess.run([str(exe)],cwd=folder,env=ENV,capture_output=True,text=True)
        if bad.returncode==0 or message not in bad.stderr:
            raise AssertionError(f'Invalid grid was not rejected correctly: {name}: {bad.stderr}')
    REPORT['checks'].append({'invalid_manual_grid_settings_rejected':True})

    # A small core can be valid even when fine-to-fine donors cover all fine interface nodes.
    small_geometry=GEOMETRY_DRIVER.replace(
        "    if (min(directCoarse,directFine,interpolatedFine)<=0) error stop 'Missing direct or interpolation interface path'",
        "    if (directCoarse<=0) error stop 'Missing fine-to-coarse transfer'")
    for ratio in (2,4,8):
        dims=dict(n=16*ratio,ny=18*ratio,walls=(7*ratio,8*ratio))  # Both central spans are exactly 2*ratio; there is no rounding.
        src=variant(source,ratio=ratio)
        folder,exe=compile_source(f'two_cell_core_geometry_r{ratio}',src,small_geometry,**dims)
        stdout=run([exe],folder)
        counts=[int(x) for x in next(line for line in stdout.splitlines()
                                     if 'ALIGNED_GEOMETRY' in line).split()[1:4]]
        for name,kw in [('rb_useg',{}),('side_legacy',{'side':True,'legacy':True})]:
            folder,exe=compile_source(f'two_cell_core_conduction_r{ratio}_{name}',variant(source,ratio=ratio,**kw),
                                      CONDUCTION.replace('step=1,300','step=1,60'),ra=10000,**dims)
            stdout=run([exe],folder,timeout=150)
            result=[float(x) for x in next(line for line in stdout.splitlines()
                                          if 'CONDUCTION_ERROR' in line).split()[1:]]
            history=np.loadtxt(folder/'NuRe_2DOpenaccMultiblock.dat')
            if not np.all(np.isfinite(history)) or max(abs(history[[1,3,4,5]]-1))>1e-10:
                raise AssertionError(f'Two-cell core conduction diagnostics: {history}')
            REPORT['checks'].append({'small_core_ratio':ratio,'requested_core_width_fine_spacings':2*ratio,
                'integration_core_width_fine_spacings':2*ratio,'conduction':name,
                'coarse_direct_and_fine_direct_interpolated_counts':counts,
                'max_T_error':result[0],'max_rho_error':result[1],
                'Nu_volume_hot_cold_middle':history[[1,3,4,5]].tolist()})
        full,fexe=compile_source(f'two_cell_core_restart_r{ratio}_full',src,RESTART_DRIVER,**dims)
        split,sexe=compile_source(f'two_cell_core_restart_r{ratio}_split',src,RESTART_DRIVER,**dims)
        _,rexe=compile_source(f'two_cell_core_restart_r{ratio}_resume',variant(source,ratio=ratio,restart=True),
                             RESTART_DRIVER,**dims)
        run([fexe,str(12*ratio)],full); run([sexe,str(6*ratio)],split)
        shutil.copy2(rexe,split/'resume.exe'); run([split/'resume.exe',str(12*ratio)],split)
        expected=np.fromfile(full/'allstate.bin',dtype='<f8'); actual=np.fromfile(split/'allstate.bin',dtype='<f8')
        if not np.all(np.isfinite(actual)) or not np.array_equal(actual,expected):
            raise AssertionError(f'Two-cell core restart mismatch: ratio {ratio}')
        REPORT['checks'].append({'small_core_ratio':ratio,'restart_exact':True})
        print('Small node-interface core',ratio,'geometry, conduction and restart passed',flush=True)

    # Individual indices need not divide by the ratio; paired indices must give whole central spans.
    for ratio in (2,3,4,8):
        for offset in range(1,ratio):
            dims=dict(n=24*ratio,ny=28*ratio,walls=(4*ratio+offset,4*ratio-offset+1,5*ratio+offset,5*ratio-offset+1))
            src=variant(source,ratio=ratio)
            folder,exe=compile_source(f'free_wall_geometry_r{ratio}_{offset}',src,GEOMETRY_DRIVER,**dims)
            run([exe],folder)
        dims=dict(n=24*ratio,ny=28*ratio,walls=(4*ratio+1,4*ratio,5*ratio+1,5*ratio))
        for name,kw in [('rb_useg',{}),('side_legacy',{'side':True,'legacy':True})]:
            folder,exe=compile_source(f'free_wall_conduction_r{ratio}_{name}',variant(source,ratio=ratio,**kw),
                                      CONDUCTION.replace('step=1,300','step=1,60'),ra=10000,**dims)
            stdout=run([exe],folder,timeout=150)
            result=[float(x) for x in next(line for line in stdout.splitlines()
                                          if 'CONDUCTION_ERROR' in line).split()[1:]]
            history=np.loadtxt(folder/'NuRe_2DOpenaccMultiblock.dat')
            if not np.all(np.isfinite(history)) or max(abs(history[[1,3,4,5]]-1))>1e-10:
                raise AssertionError(f'Nondivisible wall conduction diagnostics: {history}')
            REPORT['checks'].append({'nondivisible_wall_ratio':ratio,'conduction':name,
                'wall_thicknesses':dims['walls'],'max_T_error':result[0],'max_rho_error':result[1],
                'all_wall_residues_geometry':'passed','Nu_volume_hot_cold_middle':history[[1,3,4,5]].tolist()})
        full,fexe=compile_source(f'free_wall_restart_r{ratio}_full',src,RESTART_DRIVER,**dims)
        split,sexe=compile_source(f'free_wall_restart_r{ratio}_split',src,RESTART_DRIVER,**dims)
        _,rexe=compile_source(f'free_wall_restart_r{ratio}_resume',variant(source,ratio=ratio,restart=True),
                             RESTART_DRIVER,**dims)
        run([fexe,str(12*ratio)],full); run([sexe,str(6*ratio)],split)
        shutil.copy2(rexe,split/'resume.exe'); run([split/'resume.exe',str(12*ratio)],split)
        expected=np.fromfile(full/'allstate.bin',dtype='<f8'); actual=np.fromfile(split/'allstate.bin',dtype='<f8')
        if not np.all(np.isfinite(actual)) or not np.array_equal(actual,expected):
            raise AssertionError(f'Nondivisible wall restart mismatch: ratio {ratio}')
        REPORT['checks'].append({'nondivisible_wall_ratio':ratio,'restart_exact':True})
        print('Nondivisible wall thicknesses',ratio,'geometry, conduction and restart passed',flush=True)

    narrow=source.replace('overlapCells = 2','overlapCells = 1')
    folder,exe=compile_source('invalid_overlap',narrow,GEOMETRY_DRIVER)
    bad=subprocess.run([str(exe)],cwd=folder,env=ENV,capture_output=True,text=True)
    if bad.returncode==0 or 'Overlap is too narrow' not in bad.stderr:
        raise AssertionError('Invalid overlap was not rejected')
    REPORT['checks'].append({'overlapCells_1_rejected':True,'default_overlapCells':2})

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
    # Match the parent interface: manual numbered restart only when latest.meta is absent.
    manual_src=manual_reload_number(variant(source,restart=True),1)
    manual,mexe=compile_source('restart_manual_number',manual_src,RESTART_DRIVER)
    first_name='reloadFile2DOpenaccMultiblock-000000000001.bin'
    second_name='reloadFile2DOpenaccMultiblock-000000000002.bin'
    assert (split/first_name).exists() and (split/second_name).exists()
    shutil.copy2(split/first_name,manual/first_name)
    shutil.copy2(split/'NuRe_2DOpenaccMultiblock.dat',manual/'NuRe_2DOpenaccMultiblock.dat')
    first_hash=hashlib.sha256((manual/first_name).read_bytes()).hexdigest()
    run([mexe,'40'],manual)
    assert np.array_equal(np.fromfile(manual/'allstate.bin',dtype='<f8'),expected)
    assert (manual/'reloadFile2DOpenaccMultiblock-latest.meta').read_text().strip()==second_name
    assert hashlib.sha256((manual/first_name).read_bytes()).hexdigest()==first_hash
    # An existing meta pointer overrides even an invalid manually supplied file number.
    _,pexe=compile_source('restart_meta_priority',manual_reload_number(variant(source,restart=True),999),RESTART_DRIVER)
    run([pexe,'40'],manual)
    assert np.array_equal(np.fromfile(manual/'allstate.bin',dtype='<f8'),expected)
    missing,missing_exe=compile_source('restart_missing_selection',variant(source,restart=True),RESTART_DRIVER)
    bad=subprocess.run([str(missing_exe),'40'],cwd=missing,env=ENV,capture_output=True)
    assert bad.returncode!=0 and b'set reloadFileNum' in bad.stderr
    REPORT['checks'].append({'manual_number_restart_exact':True,'meta_priority':True,
        'unsteady_checkpoint_counter_continues':True,'earlier_checkpoint_preserved':True})
    fresh,fresh_exe=compile_source('fresh_ignores_manual_number',manual_reload_number(source,999),
                                   RESTART_DRIVER)
    run([fresh_exe,'40'],fresh)
    assert (fresh/first_name).exists()
    assert np.array_equal(np.fromfile(fresh/'allstate.bin',dtype='<f8'),expected)
    steady,steady_exe=compile_source('steady_restart_full',variant(source,steady=True),RESTART_DRIVER)
    steady_split,steady_split_exe=compile_source('steady_restart_split',variant(source,steady=True),RESTART_DRIVER)
    _,steady_resume_exe=compile_source('steady_restart_resume',manual_reload_number(variant(source,steady=True,restart=True),20),RESTART_DRIVER)
    run([steady_exe,'40'],steady); run([steady_split_exe,'20'],steady_split)
    assert (steady_split/'reloadFile2DOpenaccMultiblock-000000000020.bin').exists()
    (steady_split/'reloadFile2DOpenaccMultiblock-latest.meta').unlink()
    run([steady_resume_exe,'40'],steady_split)
    assert (steady_split/'reloadFile2DOpenaccMultiblock-000000000040.bin').exists()
    assert np.array_equal(np.fromfile(steady/'allstate.bin',dtype='<f8'),
                          np.fromfile(steady_split/'allstate.bin',dtype='<f8'))
    REPORT['checks'].append({'fresh_run_resets_reload_counter':True,'steady_step_number_manual_restart_exact':True})
    # Shift the right/top interface by one coarse spacing: valid meshes, but incompatible checkpoints.
    for label, walls in [('right',(12,15,12,13)),('top',(12,13,12,15))]:
        target,texe=compile_source('restart_changed_'+label,variant(source,restart=True),RESTART_DRIVER,walls=walls)
        shutil.copy2(split/first_name,target/first_name)
        (target/'reloadFile2DOpenaccMultiblock-latest.meta').write_text(first_name+'\n')
        shutil.copy2(split/'NuRe_2DOpenaccMultiblock.dat',target/'NuRe_2DOpenaccMultiblock.dat')
        bad=subprocess.run([str(texe),'40'],cwd=target,env=ENV,capture_output=True)
        assert bad.returncode!=0 and b'Restart mesh/refinement mismatch' in bad.stderr
    REPORT['checks'].append({'independent_right_top_restart_mismatch_rejected':True})
    # Previous layouts/statistics must fail before reading arrays (v5 used the old integration partition).
    latest=(split/'reloadFile2DOpenaccMultiblock-latest.meta').read_text().strip()
    old=split/'old-layout.bin'
    state=(split/latest).read_bytes()
    for version in (3,4,5,6,7):
        old.write_bytes(f'MB2DRESTART{version:04d}'.encode().ljust(16,b' ')+state[16:])
        (split/'reloadFile2DOpenaccMultiblock-latest.meta').write_text('old-layout.bin\n')
        bad=subprocess.run([str(split/'resume.exe'),'40'],cwd=split,env=ENV,capture_output=True)
        if bad.returncode==0 or b'Wrong checkpoint format' not in bad.stderr:
            raise AssertionError(f'Old checkpoint v{version} was not rejected')
    REPORT['checks'].append({'old_layout_restart_rejected':[3,4,5,6,7]})

    # Exercise the actual program, output clocks, binary snapshots and history-backed restart.
    smoke=source.replace('unsteadyRunDuration = 1000.0d0','unsteadyRunDuration = 0.1d0')
    smoke=smoke.replace('outputSnapshotInterval = 0.5d0','outputSnapshotInterval = 0.01d0')
    smoke=smoke.replace('reloadFileInterval = 100.0d0','reloadFileInterval = 0.03d0')
    smoke=smoke.replace('outputPltFileInterval = 100.0d0','outputPltFileInterval = 0.05d0')
    whole,wexe=compile_source('main_full',smoke)
    short,texe=compile_source('main_short',smoke.replace('unsteadyRunDuration = 0.1d0','unsteadyRunDuration = 0.05d0'))
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
    for flag,pattern in [('outputSnapshotFile','*Snapshot-*.bin'),('outputPltFile','*Tecplot-*.dat'),
                         ('outputReloadFile','reloadFile2DOpenaccMultiblock-*')]:
        folder,exe=compile_source('disabled_'+flag,smoke.replace(flag+' = 1',flag+' = 0'))
        run([exe],folder)
        assert not list(folder.glob(pattern)),flag
        assert np.array_equal(np.loadtxt(folder/'NuRe_2DOpenaccMultiblock.dat'),hfull),flag
    REPORT['checks'].append({'independent_output_switches_and_unchanged_sampling':'passed'})
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
        if re.search(r'maxBlocks\s*=\s*2\b', SOURCE.read_text(encoding='utf-8-sig')):
            from verify_ring import main as ring_main
            ring_main()
        else:
            main()
    except Exception:
        REPORT['status']='failed'
        (BUILD/'failure_report.json').write_text(json.dumps(REPORT,indent=2),encoding='utf-8')
        print('Build evidence:',BUILD,flush=True)
        raise
