!=============================================================
! 二维均匀网格耦合 DUGKS：D2Q9 流场 + D2Q9 温度，单 GPU OpenACC。
! Wang, Tao & Guo, Computers & Fluids 120 (2015) 70-81.
! 框架来自 均匀网格/2DRBOpenmp.F90；GPU 风格参考 2DRBOpenacc.F90。
! DUGKS 是有限体积方法：f/g 保存 tilde 分布，不是 LBM 迁移分布。
!=============================================================
#define steadyFlow
!#define unsteadyFlow
#define SideHeatedCell
!#define RayleighBenardCell
!#define VerticalWallsPeriodical   ! 仅 RB 可选，速度和温度同时周期
#if defined(steadyFlow) == defined(unsteadyFlow)
#error "Choose exactly one flow mode"
#endif
#if defined(SideHeatedCell) == defined(RayleighBenardCell)
#error "Choose exactly one convection case"
#endif
#if defined(SideHeatedCell) && defined(VerticalWallsPeriodical)
#error "Side-heated cavity cannot use periodic x boundaries"
#endif

module commondata
    use, intrinsic :: ieee_arithmetic
    implicit none
    integer, parameter :: nx=256, ny=256               ! 直接修改二维均匀网格单元数
    integer, parameter :: loadInitField=0              ! 0：新算；1：读取 latest.meta
    real(8), parameter :: dx=1.0d0, dy=1.0d0          ! 均匀网格；物理壁面在最外层单元面
    real(8), parameter :: Rayleigh=1.0d6, Prandtl=0.71d0, Mach=0.01d0
    real(8), parameter :: Thot=0.5d0, Tcold=-0.5d0, Tref=0.5d0*(Thot+Tcold)
    real(8), parameter :: temperatureDifference=Thot-Tcold, rhoRef=1.0d0
    real(8), parameter :: cs2=1.0d0/3.0d0             ! 固定格子常数，不是可调参数；D2Q9 速度为0、正负1
    real(8), parameter :: Courant=0.5d0               ! dt=CFL*min(dx,dy)/sqrt(2)
    real(8), parameter :: dt=Courant*min(dx,dy)/sqrt(6.0d0*cs2), h=0.5d0*dt
#ifdef SideHeatedCell
    real(8), parameter :: lengthUnit=dble(nx)*dx       ! 侧壁差温的特征长度
    integer, parameter :: caseCode=1
#else
    real(8), parameter :: lengthUnit=dble(ny)*dy       ! RB 的特征长度为高度
    integer, parameter :: caseCode=2
#endif
    real(8), parameter :: velocityUnit=Mach*sqrt(cs2)
    real(8), parameter :: viscosity=velocityUnit*lengthUnit*sqrt(Prandtl/Rayleigh)
    real(8), parameter :: diffusivity=viscosity/Prandtl
    real(8), parameter :: tauv=viscosity/cs2, tauc=diffusivity/cs2 ! DUGKS 没有 LBM 的减半步修正
    real(8), parameter :: gBeta=velocityUnit**2/(lengthUnit*temperatureDifference)
    real(8), parameter :: timeUnit=lengthUnit/velocityUnit ! 自由落体时间；时间推进须乘 dt
    real(8), parameter :: velocityScaleCompare=lengthUnit/diffusivity
    real(8), parameter :: pi=acos(-1.0d0), rbInitialPerturbation=1.0d-3
    real(8), parameter :: epsU=1.0d-7, epsT=1.0d-7
    integer, parameter :: steadyCheckIntervalItc=2000  ! 每次比较相隔此数量的 DUGKS 步
    integer, parameter :: steadyConsecutiveChecks=3   ! 连续通过后才判稳态收敛
    integer, parameter :: outputSnapshotFile=1, outputPltFile=1, outputReloadFile=1
    real(8), parameter :: reloadFileInterval=100.0d0, outputPltFileInterval=100.0d0
#ifdef steadyFlow
    integer, parameter :: flowCode=1, itc_max=20000000
    real(8), parameter :: outputSnapshotInterval=10.0d0
    real(8), parameter :: unsteadyRunDuration=1000.0d0
#else
    integer, parameter :: flowCode=2
    real(8), parameter :: outputSnapshotInterval=0.5d0, unsteadyRunDuration=1000.0d0
    integer, parameter :: itc_max=ceiling(unsteadyRunDuration*timeUnit/dt)
#endif
    real(8), parameter :: unsteadyAverageStartTf=0.5d0*unsteadyRunDuration
    real(8), parameter :: unsteadyAverageEndTf=unsteadyRunDuration
#ifdef VerticalWallsPeriodical
    integer, parameter :: periodicCode=1
#else
    integer, parameter :: periodicCode=0
#endif
    integer, parameter :: historyCapacity=ceiling(dble(itc_max)*dt/(outputSnapshotInterval*timeUnit))+3
    integer, parameter :: historyColumns=10
    character(len=32), parameter :: reloadMagic='DUGKS2D_D2Q9_D2Q9_V1'
    character(len=*), parameter :: reloadFilePrefix='reloadFile2DOpenaccDUGKS'
    character(len=*), parameter :: settingsFile='SimulationSettings2DOpenaccDUGKS.txt'
    character(len=*), parameter :: historyFile='NuRe_2DOpenaccDUGKS.dat'
    integer :: itc=0, restartItcOffset=0, dimensionlessTime=0
    integer :: snapshotFileNum=0, pltFileNum=0, reloadFileNum=0, convergedChecks=0
    integer :: lastSampleItc=-1, lastFluxItc=-1
    real(8) :: errorU=huge(1.0d0), errorT=huge(1.0d0)
    real(8) :: NuVolAvg, ReVolAvg, Nu_global, Nu_hot, Nu_cold, Nu_middle
    real(8) :: history(historyColumns,historyCapacity)
    real(8) :: xp(nx), yp(ny)
    integer :: ex(0:8)=[0,1,0,-1,0,1,-1,-1,1]
    integer :: ey(0:8)=[0,0,1,0,-1,1,1,-1,-1]
    integer :: opposite(0:8)=[0,3,4,1,2,7,8,5,6]
    real(8) :: omega(0:8)=[4d0/9,1d0/9,1d0/9,1d0/9,1d0/9,1d0/36,1d0/36,1d0/36,1d0/36]
    real(8), allocatable :: f(:,:,:),g(:,:,:),f_barplus(:,:,:),g_barplus(:,:,:)
    real(8), allocatable :: fsx(:,:,:),fsy(:,:,:),gsx(:,:,:),gsy(:,:,:)
    real(8), allocatable :: fluxFx(:,:,:),fluxFy(:,:,:),fluxGx(:,:,:),fluxGy(:,:,:)
    real(8), allocatable :: u(:,:),v(:,:),T(:,:),rho(:,:),up(:,:),vp(:,:),Tp(:,:)
    real(8) :: NuHotProfile(max(nx,ny)), NuColdProfile(max(nx,ny))
contains
    ! 设备顺序函数只接收标量权重，避免设备例程隐式捕获权重数组。
    real(8) function equilibrium(weight,cx,cy,scalar,ux,uy) result(eq)
        !$acc routine seq
        real(8), intent(in) :: weight,cx,cy,scalar,ux,uy
        real(8) :: cu
        cu=cx*ux+cy*uy
        eq=weight*scalar*(1.0d0+cu/cs2+0.5d0*cu*cu/(cs2*cs2)-0.5d0*(ux*ux+uy*uy)/cs2)
    end function equilibrium

    subroutine initial()
        integer :: i,j,k,ios,unit
        real(8) :: feq,source,cx,cy,ay
        if(nx<4.or.ny<4.or.Courant<=0.or.Courant>1.or.Rayleigh<=0.or.Prandtl<=0.or.Mach<=0) &
            error stop 'Invalid mesh or physical parameters'
        if(temperatureDifference<=0.or.dx<=0.or.dy<=0.or.steadyCheckIntervalItc<1) &
            error stop 'Invalid temperature, spacing, or check interval'
        if(abs(cs2-1d0/3d0)>epsilon(cs2)) error stop 'Fixed D2Q9 velocities require cs2=1/3'
        if(outputSnapshotInterval<=0.or.reloadFileInterval<=0.or.outputPltFileInterval<=0) &
            error stop 'Output intervals must be positive'
        allocate(f(nx,ny,0:8),g(nx,ny,0:8),f_barplus(nx,ny,0:8),g_barplus(nx,ny,0:8))
        allocate(fsx(nx,ny,0:8),fsy(nx,ny,0:8),gsx(nx,ny,0:8),gsy(nx,ny,0:8))
        allocate(fluxFx(0:nx,ny,0:8),fluxGx(0:nx,ny,0:8))
        allocate(fluxFy(nx,0:ny,0:8),fluxGy(nx,0:ny,0:8))
        allocate(u(nx,ny),v(nx,ny),T(nx,ny),rho(nx,ny),up(nx,ny),vp(nx,ny),Tp(nx,ny))
        do i=1,nx
            xp(i)=(dble(i)-0.5d0)*dx/lengthUnit
        enddo
        do j=1,ny
            yp(j)=(dble(j)-0.5d0)*dy/lengthUnit
        enddo
        u=0; v=0; rho=rhoRef; history=0
        do j=1,ny
            do i=1,nx
#ifdef SideHeatedCell
                T(i,j)=Thot-temperatureDifference*(dble(i)-0.5d0)/dble(nx)
#else
                T(i,j)=Thot-temperatureDifference*(dble(j)-0.5d0)/dble(ny) &
                    +rbInitialPerturbation*temperatureDifference*sin(pi*(dble(j)-0.5d0)/dble(ny)) &
                    *cos(2*pi*(dble(i)-0.5d0)/dble(nx))
#endif
                ay=gBeta*(T(i,j)-Tref)
                do k=0,8
                    cx=dble(ex(k)); cy=dble(ey(k))
                    feq=equilibrium(omega(k),cx,cy,rhoRef,0d0,0d0)
                    source=cy*ay/cs2*feq
                    f(i,j,k)=feq-0.5d0*dt*source ! 从物理平衡分布变换，物理初速严格为0
                    g(i,j,k)=equilibrium(omega(k),cx,cy,T(i,j),0d0,0d0)
                enddo
            enddo
        enddo
        up=u; vp=v; Tp=T
        if(loadInitField==1) call read_reload()
        if(loadInitField==1) call reconcile_convergence_history()
        ! 历史直接由检查点中的缓存重建，断电后多出的尾部不会重复进入统计。
        open(newunit=unit,file=historyFile,status='replace',action='write',iostat=ios)
        if(ios/=0) error stop 'Cannot open history file'
        write(unit,'(A)') '# itc t_ff NuVolAvg ReVolAvg Nu_global Nu_hot Nu_cold Nu_middle mass Tmin'
        do i=1,dimensionlessTime
            write(unit,'(10ES25.16E3)') history(:,i)
        enddo
        close(unit)
        if(loadInitField==0) then
            open(newunit=unit,file='Convergence_2DOpenaccDUGKS.dat',status='replace',action='write')
            write(unit,'(A)') '# itc t_ff errorU errorT'
            close(unit)
        endif
        if(loadInitField==0) then
            open(newunit=unit,file=settingsFile,status='replace',action='write')
        else
            open(newunit=unit,file=settingsFile,status='unknown',position='append',action='write')
        endif
        write(unit,*) 'DUGKS Wang/Tao/Guo (2015), D2Q9 + D2Q9; uniform finite volumes'
        write(unit,*) 'nx ny case flow periodic = ',nx,ny,caseCode,flowCode,periodicCode
        write(unit,*) 'Ra Pr Mach CFL dt dx dy = ',Rayleigh,Prandtl,Mach,Courant,dt,dx,dy
        write(unit,*) 'nu kappa tauv tauc Uff tff = ',viscosity,diffusivity,tauv,tauc,velocityUnit,timeUnit
        write(unit,*) 'Restart at step = ',itc
        close(unit)
    end subroutine initial

    subroutine enter_data_2d_openacc()
        !$acc enter data copyin(ex,ey,opposite,omega,f,g,u,v,T,rho)
        !$acc enter data create(f_barplus,g_barplus,fsx,fsy,gsx,gsy,fluxFx,fluxFy,fluxGx,fluxGy)
    end subroutine enter_data_2d_openacc

    subroutine update_host_2d_openacc()
        !$acc wait(1)
        !$acc update self(u,v,T,rho)
    end subroutine update_host_2d_openacc

    subroutine collision()
        ! Eq.(24),(38)：两套 bar+ 使用同一个 t^n 的物理宏观场。
        integer :: i,j,k
        real(8) :: feq,geq,source,cx,cy,ay
        !$acc parallel loop collapse(2) present(f,g,f_barplus,g_barplus,u,v,T,rho,ex,ey,omega) &
        !$acc& private(k,feq,geq,source,cx,cy,ay) async(1)
        do j=1,ny
            do i=1,nx
                ay=gBeta*(T(i,j)-Tref)
                do k=0,8
                    cx=dble(ex(k)); cy=dble(ey(k))
                    feq=equilibrium(omega(k),cx,cy,rho(i,j),u(i,j),v(i,j))
                    geq=equilibrium(omega(k),cx,cy,T(i,j),u(i,j),v(i,j))
                    source=(cy-v(i,j))*ay/cs2*feq
                    f_barplus(i,j,k)=((2*tauv-h)*f(i,j,k)+3*h*feq+3*h*tauv*source)/(2*tauv+dt)
                    g_barplus(i,j,k)=((2*tauc-h)*g(i,j,k)+3*h*geq)/(2*tauc+dt)
                enddo
            enddo
        enddo
    end subroutine collision

    subroutine reconstruction_slopes()
        ! 均匀网格线性二阶重构；内部中心差分，物理壁面邻层单边斜率。
        integer :: i,j,k,im,ip,jm,jp
        real(8) :: sx,sy
        !$acc parallel loop collapse(3) present(f_barplus,g_barplus,fsx,fsy,gsx,gsy) &
        !$acc& private(im,ip,jm,jp,sx,sy) async(1)
        do k=0,8
            do j=1,ny
                do i=1,nx
                    im=max(1,i-1); ip=min(nx,i+1); sx=dble(ip-im)*dx
#ifdef VerticalWallsPeriodical
                    im=modulo(i-2,nx)+1; ip=modulo(i,nx)+1; sx=2*dx
#endif
                    jm=max(1,j-1); jp=min(ny,j+1); sy=dble(jp-jm)*dy
                    fsx(i,j,k)=(f_barplus(ip,j,k)-f_barplus(im,j,k))/sx
                    fsy(i,j,k)=(f_barplus(i,jp,k)-f_barplus(i,jm,k))/sy
                    gsx(i,j,k)=(g_barplus(ip,j,k)-g_barplus(im,j,k))/sx
                    gsy(i,j,k)=(g_barplus(i,jp,k)-g_barplus(i,jm,k))/sy
                enddo
            enddo
        enddo
    end subroutine reconstruction_slopes

    subroutine face_flux_x()
        integer :: i,j,k,ii,il,ir,opp
        real(8) :: fb(0:8),gb(0:8),cx,cy,rx,ry,rr,tt,ux,uy,ay,tw,feq,geq,source,uwy,rhow,bsum
        logical :: wall,constantT,incoming
        !$acc parallel loop collapse(2) present(f_barplus,g_barplus,fsx,fsy,gsx,gsy,ex,ey,opposite,omega, &
        !$acc& fluxFx,fluxGx) private(k,ii,il,ir,opp,fb,gb,cx,cy,rx,ry,rr,tt,ux,uy,ay,tw,feq,geq, &
        !$acc& source,uwy,rhow,bsum,wall,constantT,incoming) async(1)
        do j=1,ny
            do i=0,nx
                wall=(i==0.or.i==nx)
#ifdef VerticalWallsPeriodical
                wall=.false.
#endif
                do k=0,8
                    cx=dble(ex(k)); cy=dble(ey(k))
#ifdef VerticalWallsPeriodical
                    il=modulo(i-1,nx)+1; ir=modulo(i,nx)+1
#else
                    il=max(1,i); ir=min(nx,i+1)
#endif
                    ry=-h*cy
                    if(cx>0) then
                        ii=il; rx=0.5d0*dx-h*cx
                        if(i==0.and.wall) rx=-0.5d0*dx-h*cx
                    elseif(cx<0) then
                        ii=ir; rx=-0.5d0*dx-h*cx
                        if(i==nx.and.wall) rx=0.5d0*dx-h*cx
                    else
                        ! 切向粒子从两侧线性重构的平均值取面值。
                        if(wall) then
                            ii=il; rx=0.5d0*dx
                            if(i==0) rx=-0.5d0*dx
                            fb(k)=f_barplus(ii,j,k)+rx*fsx(ii,j,k)+ry*fsy(ii,j,k)
                            gb(k)=g_barplus(ii,j,k)+rx*gsx(ii,j,k)+ry*gsy(ii,j,k)
                            cycle
                        endif
                        rx=0.5d0*dx
                        if(i==0.and.wall) rx=-0.5d0*dx
                        fb(k)=0.5d0*(f_barplus(il,j,k)+rx*fsx(il,j,k)+ry*fsy(il,j,k) &
                            +f_barplus(ir,j,k)-rx*fsx(ir,j,k)+ry*fsy(ir,j,k))
                        gb(k)=0.5d0*(g_barplus(il,j,k)+rx*gsx(il,j,k)+ry*gsy(il,j,k) &
                            +g_barplus(ir,j,k)-rx*gsx(ir,j,k)+ry*gsy(ir,j,k))
                        cycle
                    endif
                    fb(k)=f_barplus(ii,j,k)+rx*fsx(ii,j,k)+ry*fsy(ii,j,k)
                    gb(k)=g_barplus(ii,j,k)+rx*gsx(ii,j,k)+ry*gsy(ii,j,k)
                enddo
                constantT=.false.; tw=0d0
#ifdef SideHeatedCell
                constantT=wall
                tw=Thot
                if(i==nx) tw=Tcold
#endif
                if(wall) then
                    do k=0,8
                        incoming=(i==0.and.ex(k)>0).or.(i==nx.and.ex(k)<0)
                        if(incoming) then
                            opp=opposite(k)
                            if(constantT) then
                                gb(k)=-gb(opp)+2*omega(k)*tw ! Eq.(43)，静止恒温壁
                            else
                                gb(k)=gb(opp)               ! Eq.(44)，绝热壁
                            endif
                        endif
                    enddo
                endif
                tt=sum(gb)
                if(constantT) tt=tw
                ay=gBeta*(tt-Tref)
                if(wall) then
                    uwy=-0.5d0*h*ay ! Eq.(42)：bar 变量壁速 uw'=uw-h*a/2
                    ! 用边界关系本身解 rho_w，不采用论文近不可压缩时的 rho_w≈rhoRef 简化。
                    ! 这样 Eq.(23) 恢复后法向质量通量在舍入精度内严格为0。
                    rhow=0d0; bsum=0d0
                    do k=0,8
                        incoming=(i==0.and.ex(k)>0).or.(i==nx.and.ex(k)<0)
                        if(incoming) then
                            rhow=rhow+fb(opposite(k))
                            bsum=bsum+2*omega(k)*dble(ey(k))*uwy/cs2
                        else
                            rhow=rhow+fb(k)
                        endif
                    enddo
                    rhow=rhow/(1-bsum)
                    do k=0,8
                        incoming=(i==0.and.ex(k)>0).or.(i==nx.and.ex(k)<0)
                        if(incoming) fb(k)=fb(opposite(k))+2*rhow*omega(k)*dble(ey(k))*uwy/cs2
                    enddo
                endif
                rr=sum(fb); ux=0; uy=0
                do k=0,8
                    ux=ux+dble(ex(k))*fb(k); uy=uy+dble(ey(k))*fb(k)
                enddo
                ux=ux/rr; uy=uy/rr+0.5d0*h*ay
                if(wall) then
                    ux=0; uy=0
                endif
                do k=0,8
                    cx=dble(ex(k)); cy=dble(ey(k))
                    feq=equilibrium(omega(k),cx,cy,rr,ux,uy)
                    geq=equilibrium(omega(k),cx,cy,tt,ux,uy)
                    source=(cy-uy)*ay/cs2*feq
                    ! Eq.(23),(37)：真正半步物理分布，含碰撞与受力耦合。
                    fluxFx(i,j,k)=cx*(2*tauv*fb(k)+h*feq+tauv*h*source)/(2*tauv+h)
                    fluxGx(i,j,k)=cx*(2*tauc*gb(k)+h*geq)/(2*tauc+h)
                enddo
            enddo
        enddo
    end subroutine face_flux_x

    subroutine face_flux_y()
        integer :: i,j,k,jj,jb,jt,opp
        real(8) :: fb(0:8),gb(0:8),cx,cy,rx,ry,rr,tt,ux,uy,ay,tw,feq,geq,source,uwy,rhow,bsum
        logical :: wall,constantT,incoming
        !$acc parallel loop collapse(2) present(f_barplus,g_barplus,fsx,fsy,gsx,gsy,ex,ey,opposite,omega, &
        !$acc& fluxFy,fluxGy) private(k,jj,jb,jt,opp,fb,gb,cx,cy,rx,ry,rr,tt,ux,uy,ay,tw,feq,geq, &
        !$acc& source,uwy,rhow,bsum,wall,constantT,incoming) async(1)
        do j=0,ny
            do i=1,nx
                wall=(j==0.or.j==ny)
                do k=0,8
                    cx=dble(ex(k)); cy=dble(ey(k)); jb=max(1,j); jt=min(ny,j+1)
                    rx=-h*cx
                    if(cy>0) then
                        jj=jb; ry=0.5d0*dy-h*cy
                        if(j==0) ry=-0.5d0*dy-h*cy
                    elseif(cy<0) then
                        jj=jt; ry=-0.5d0*dy-h*cy
                        if(j==ny) ry=0.5d0*dy-h*cy
                    else
                        if(wall) then
                            jj=jb; ry=0.5d0*dy
                            if(j==0) ry=-0.5d0*dy
                            fb(k)=f_barplus(i,jj,k)+rx*fsx(i,jj,k)+ry*fsy(i,jj,k)
                            gb(k)=g_barplus(i,jj,k)+rx*gsx(i,jj,k)+ry*gsy(i,jj,k)
                            cycle
                        endif
                        ry=0.5d0*dy
                        if(j==0) ry=-0.5d0*dy
                        fb(k)=0.5d0*(f_barplus(i,jb,k)+rx*fsx(i,jb,k)+ry*fsy(i,jb,k) &
                            +f_barplus(i,jt,k)+rx*fsx(i,jt,k)-ry*fsy(i,jt,k))
                        gb(k)=0.5d0*(g_barplus(i,jb,k)+rx*gsx(i,jb,k)+ry*gsy(i,jb,k) &
                            +g_barplus(i,jt,k)+rx*gsx(i,jt,k)-ry*gsy(i,jt,k))
                        cycle
                    endif
                    fb(k)=f_barplus(i,jj,k)+rx*fsx(i,jj,k)+ry*fsy(i,jj,k)
                    gb(k)=g_barplus(i,jj,k)+rx*gsx(i,jj,k)+ry*gsy(i,jj,k)
                enddo
                constantT=.false.; tw=0d0
#ifdef RayleighBenardCell
                constantT=wall; tw=Thot
                if(j==ny) tw=Tcold
#endif
                if(wall) then
                    do k=0,8
                        incoming=(j==0.and.ey(k)>0).or.(j==ny.and.ey(k)<0)
                        if(incoming) then
                            opp=opposite(k)
                            if(constantT) then
                                gb(k)=-gb(opp)+2*omega(k)*tw
                            else
                                gb(k)=gb(opp)
                            endif
                        endif
                    enddo
                endif
                tt=sum(gb)
                if(constantT) tt=tw
                ay=gBeta*(tt-Tref)
                if(wall) then
                    uwy=-0.5d0*h*ay
                    rhow=0d0; bsum=0d0
                    do k=0,8
                        incoming=(j==0.and.ey(k)>0).or.(j==ny.and.ey(k)<0)
                        if(incoming) then
                            rhow=rhow+fb(opposite(k))
                            bsum=bsum+2*omega(k)*dble(ey(k))*uwy/cs2
                        else
                            rhow=rhow+fb(k)
                        endif
                    enddo
                    rhow=rhow/(1-bsum)
                    do k=0,8
                        incoming=(j==0.and.ey(k)>0).or.(j==ny.and.ey(k)<0)
                        if(incoming) fb(k)=fb(opposite(k))+2*rhow*omega(k)*dble(ey(k))*uwy/cs2
                    enddo
                endif
                rr=sum(fb); ux=0; uy=0
                do k=0,8
                    ux=ux+dble(ex(k))*fb(k); uy=uy+dble(ey(k))*fb(k)
                enddo
                ux=ux/rr; uy=uy/rr+0.5d0*h*ay
                if(wall) then
                    ux=0; uy=0
                endif
                do k=0,8
                    cx=dble(ex(k)); cy=dble(ey(k))
                    feq=equilibrium(omega(k),cx,cy,rr,ux,uy)
                    geq=equilibrium(omega(k),cx,cy,tt,ux,uy)
                    source=(cy-uy)*ay/cs2*feq
                    fluxFy(i,j,k)=cy*(2*tauv*fb(k)+h*feq+tauv*h*source)/(2*tauv+h)
                    fluxGy(i,j,k)=cy*(2*tauc*gb(k)+h*geq)/(2*tauc+h)
                enddo
            enddo
        enddo
    end subroutine face_flux_y

    subroutine streaming()
        ! Eq.(13),(25),(30),(39)：各面通量唯一共享，保证相邻单元严格守恒。
        integer :: i,j,k
        !$acc parallel loop collapse(3) present(f,g,f_barplus,g_barplus,fluxFx,fluxFy,fluxGx,fluxGy) async(1)
        do k=0,8
            do j=1,ny
                do i=1,nx
                    f(i,j,k)=(4*f_barplus(i,j,k)-f(i,j,k))/3 &
                        -dt/dx*(fluxFx(i,j,k)-fluxFx(i-1,j,k))-dt/dy*(fluxFy(i,j,k)-fluxFy(i,j-1,k))
                    g(i,j,k)=(4*g_barplus(i,j,k)-g(i,j,k))/3 &
                        -dt/dx*(fluxGx(i,j,k)-fluxGx(i-1,j,k))-dt/dy*(fluxGy(i,j,k)-fluxGy(i,j-1,k))
                enddo
            enddo
        enddo
        lastFluxItc=itc ! 这些面通量对应 (itc-1/2)*dt
    end subroutine streaming

    subroutine macro()
        integer :: i,j,k
        real(8) :: rr,tt,mx,my
        !$acc parallel loop collapse(2) present(f,g,rho,u,v,T,ex,ey) private(k,rr,tt,mx,my) async(1)
        do j=1,ny
            do i=1,nx
                rr=0; tt=0; mx=0; my=0
                do k=0,8
                    rr=rr+f(i,j,k); tt=tt+g(i,j,k)
                    mx=mx+dble(ex(k))*f(i,j,k); my=my+dble(ey(k))*f(i,j,k)
                enddo
                rho(i,j)=rr; T(i,j)=tt
                u(i,j)=mx/rr; v(i,j)=my/rr+0.5d0*dt*gBeta*(tt-Tref) ! Eq.(46)
            enddo
        enddo
    end subroutine macro

    subroutine check_fields()
        call update_host_2d_openacc()
        if(any(.not.ieee_is_finite(u)).or.any(.not.ieee_is_finite(v)).or. &
            any(.not.ieee_is_finite(T)).or.any(.not.ieee_is_finite(rho))) error stop 'Nonfinite physical field'
        if(minval(rho)<=0) error stop 'Nonpositive density'
    end subroutine check_fields

    subroutine check()
        integer :: unit
        call check_fields()
        errorU=sqrt(sum((u-up)**2+(v-vp)**2)/max(sum(u*u+v*v),1d-60))
        errorT=sqrt(sum((T-Tp)**2)/max(sum((T-Tref)**2),1d-60))
        up=u; vp=v; Tp=T
        if(errorU<epsU.and.errorT<epsT) then
            convergedChecks=convergedChecks+1
        else
            convergedChecks=0
        endif
        open(newunit=unit,file='Convergence_2DOpenaccDUGKS.dat',status='unknown',position='append')
        write(unit,'(I12,3ES25.16E3)') itc,itc*dt/timeUnit,errorU,errorT
        close(unit)
        write(*,'(A,I10,A,F12.6,A,2ES12.4)') 'itc=',itc,' t_ff=',itc*dt/timeUnit,' errorU/errorT=',errorU,errorT
    end subroutine check

    subroutine calNuRe(saveSample)
        logical, intent(in) :: saveSample
        integer :: i,j,unit,mid,nwall
        real(8) :: coef,weight,sumq,middle1,middle2
        call check_fields()
        coef=lengthUnit/(diffusivity*temperatureDifference)
#ifdef SideHeatedCell
        NuVolAvg=1+coef*sum(u*(T-Tref))/dble(nx*ny)
#else
        NuVolAvg=1+coef*sum(v*(T-Tref))/dble(nx*ny)
#endif
        ReVolAvg=lengthUnit/viscosity*sqrt(sum(u*u+v*v)/dble(nx*ny))
        ! 宏观量采样于整步 t^n；由该状态重构的面通量是 t^n+h 的半步预测值。
        ! Nu_global/壁 Nu 的预测时间单独说明；重构不改变 f/g 或宏观场。
        call collision(); call reconstruction_slopes(); call face_flux_x(); call face_flux_y()
        !$acc wait(1)
        !$acc update self(fluxGx,fluxGy)
        sumq=0; Nu_hot=0; Nu_cold=0; Nu_middle=0
#ifdef SideHeatedCell
        nwall=ny; mid=nx/2
        do j=1,ny
            NuHotProfile(j)=coef*sum(fluxGx(0,j,:))
            NuColdProfile(j)=coef*sum(fluxGx(nx,j,:))
            do i=0,nx
                weight=1d0
                if(i==0.or.i==nx) weight=0.5d0
                sumq=sumq+weight*sum(fluxGx(i,j,:))
            enddo
            middle1=sum(fluxGx(mid,j,:)); middle2=middle1
            if(mod(nx,2)==1) middle2=sum(fluxGx(mid+1,j,:))
            Nu_middle=Nu_middle+coef*0.5d0*(middle1+middle2)/dble(ny)
        enddo
        Nu_global=coef*sumq/dble(nx*ny)
#else
        nwall=nx; mid=ny/2
        do i=1,nx
            NuHotProfile(i)=coef*sum(fluxGy(i,0,:))
            NuColdProfile(i)=coef*sum(fluxGy(i,ny,:))
            do j=0,ny
                weight=1d0
                if(j==0.or.j==ny) weight=0.5d0
                sumq=sumq+weight*sum(fluxGy(i,j,:))
            enddo
            middle1=sum(fluxGy(i,mid,:)); middle2=middle1
            if(mod(ny,2)==1) middle2=sum(fluxGy(i,mid+1,:))
            Nu_middle=Nu_middle+coef*0.5d0*(middle1+middle2)/dble(nx)
        enddo
        Nu_global=coef*sumq/dble(nx*ny)
#endif
        Nu_hot=sum(NuHotProfile(1:nwall))/dble(nwall)
        Nu_cold=sum(NuColdProfile(1:nwall))/dble(nwall)
        if(.not.ieee_is_finite(Nu_global).or..not.ieee_is_finite(ReVolAvg)) error stop 'Nonfinite diagnostics'
        if(saveSample.and.lastSampleItc/=itc) then
            dimensionlessTime=dimensionlessTime+1
            if(dimensionlessTime>historyCapacity) error stop 'History capacity exceeded'
            history(:,dimensionlessTime)=[dble(itc),itc*dt/timeUnit,NuVolAvg,ReVolAvg,Nu_global, &
                Nu_hot,Nu_cold,Nu_middle,sum(rho)*dx*dy,minval(T)]
            open(newunit=unit,file=historyFile,status='old',position='append')
            write(unit,'(10ES25.16E3)') history(:,dimensionlessTime)
            close(unit)
            lastSampleItc=itc
        endif
    end subroutine calNuRe

    subroutine output_SnapshotFile()
        integer :: unit
        character(len=160) :: filename
        call update_host_2d_openacc()
        snapshotFileNum=snapshotFileNum+1
        write(filename,'(A,I8.8,A)') 'buoyancyCavity2DOpenaccDUGKSSnapshot-',snapshotFileNum,'.bin'
        open(newunit=unit,file=trim(filename),status='replace',access='stream',form='unformatted')
        write(unit) 'DUGKS2D_SNAPSHOT_V1',nx,ny,itc,itc*dt/timeUnit,xp,yp,u,v,T,rho
        close(unit)
    end subroutine output_SnapshotFile

    subroutine output_Tecplot()
        integer :: i,j,unit
        real(8) :: vort,psi(nx,ny),dvdx,dudy
        character(len=160) :: filename
        call update_host_2d_openacc()
        pltFileNum=pltFileNum+1
        write(filename,'(A,I8.8,A)') 'buoyancyCavity2DOpenaccDUGKSTecplot-',pltFileNum,'.plt'
        psi=0d0
        do i=1,nx
            psi(i,1)=0.25d0*dy*u(i,1) ! 壁面u=0到首个单元中心的梯形积分
            do j=2,ny
                psi(i,j)=psi(i,j-1)+0.5d0*dy*(u(i,j-1)+u(i,j))
            enddo
        enddo
        open(newunit=unit,file=trim(filename),status='replace',action='write')
        write(unit,'(A)') 'VARIABLES="x/L","y/L","u/Uff","v/Uff","T","rho","psi/(Uff L)","vort L/Uff"'
        write(unit,'(A,I0,A,I0,A,ES24.16)') 'ZONE I=',nx,', J=',ny,', DATAPACKING=POINT, SOLUTIONTIME=',itc*dt/timeUnit
        do j=1,ny
            do i=1,nx
                dvdx=(v(min(nx,i+1),j)-v(max(1,i-1),j))/(dble(min(nx,i+1)-max(1,i-1))*dx)
                dudy=(u(i,min(ny,j+1))-u(i,max(1,j-1)))/(dble(min(ny,j+1)-max(1,j-1))*dy)
                if(i==1) dvdx=(-3*v(1,j)+4*v(2,j)-v(3,j))/(2*dx)
                if(i==nx) dvdx=(3*v(nx,j)-4*v(nx-1,j)+v(nx-2,j))/(2*dx)
                if(j==1) dudy=(-3*u(i,1)+4*u(i,2)-u(i,3))/(2*dy)
                if(j==ny) dudy=(3*u(i,ny)-4*u(i,ny-1)+u(i,ny-2))/(2*dy)
#ifdef VerticalWallsPeriodical
                dvdx=(v(modulo(i,nx)+1,j)-v(modulo(i-2,nx)+1,j))/(2*dx)
#endif
                vort=(dvdx-dudy)*lengthUnit/velocityUnit
                write(unit,'(8ES25.16E3)') xp(i),yp(j),u(i,j)/velocityUnit,v(i,j)/velocityUnit,T(i,j),rho(i,j), &
                    psi(i,j)/(velocityUnit*lengthUnit),vort
            enddo
        enddo
        close(unit)
    end subroutine output_Tecplot

    subroutine final_diagnostics()
        integer :: unit,i,j,iu,iv,loc(1),nwall
        real(8) :: umid(ny),vmid(nx),uxmax,vymax,posu,posv
        call calNuRe(.false.) ! 段末非采样点只写最终诊断，不能改变固定时间序列的采样相位
        iu=nx/2; iv=ny/2
        if(mod(nx,2)==0) then
            umid=0.5d0*(u(iu,:)+u(iu+1,:))*velocityScaleCompare
        else
            umid=u(iu+1,:)*velocityScaleCompare
        endif
        if(mod(ny,2)==0) then
            vmid=0.5d0*(v(:,iv)+v(:,iv+1))*velocityScaleCompare
        else
            vmid=v(:,iv+1)*velocityScaleCompare
        endif
        loc=maxloc(umid); uxmax=umid(loc(1)); posu=yp(loc(1))
        loc=maxloc(vmid); vymax=vmid(loc(1)); posv=xp(loc(1))
        open(newunit=unit,file='FinalDiagnostics_2DOpenaccDUGKS.dat',status='replace')
        write(unit,*) 'itc t_ff errorU errorT consecutive_convergence_checks'
        write(unit,*) itc,itc*dt/timeUnit,errorU,errorT,convergedChecks
        write(unit,*) 'NuVolAvg ReVolAvg Nu_global Nu_hot Nu_cold Nu_middle'
        write(unit,'(6ES25.16E3)') NuVolAvg,ReVolAvg,Nu_global,Nu_hot,Nu_cold,Nu_middle
        write(unit,*) 'umid_max(UL/kappa) y/L vmid_max(VL/kappa) x/L (cell-center locations)'
        write(unit,'(4ES25.16E3)') uxmax,posu,vymax,posv
        write(unit,*) 'Density mass = ',sum(rho)*dx*dy,' ; Tmin/Tmax = ',minval(T),maxval(T)
        write(unit,*) 'Flux diagnostics use the t^n+h DUGKS face predictor from the t^n state.'
        close(unit)
#ifdef SideHeatedCell
        nwall=ny
#else
        nwall=nx
#endif
        open(newunit=unit,file='WallNuExtrema_2DOpenaccDUGKS.dat',status='replace')
        write(unit,'(A)') '# wall maxNu position/L minNu position/L (cell-center samples)'
        do j=1,2
            if(j==1) then
                loc=maxloc(NuHotProfile(1:nwall)); i=loc(1); uxmax=NuHotProfile(i)
#ifdef SideHeatedCell
                posu=yp(i)
#else
                posu=xp(i)
#endif
                loc=minloc(NuHotProfile(1:nwall)); i=loc(1); vymax=NuHotProfile(i)
            else
                loc=maxloc(NuColdProfile(1:nwall)); i=loc(1); uxmax=NuColdProfile(i)
#ifdef SideHeatedCell
                posu=yp(i)
#else
                posu=xp(i)
#endif
                loc=minloc(NuColdProfile(1:nwall)); i=loc(1); vymax=NuColdProfile(i)
            endif
#ifdef SideHeatedCell
            posv=yp(i)
#else
            posv=xp(i)
#endif
            write(unit,'(I3,4ES25.16E3)') j,uxmax,posu,vymax,posv
        enddo
        close(unit)
        open(newunit=unit,file='WallNu_2DOpenaccDUGKS.dat',status='replace')
        write(unit,'(A)') '# wall-coordinate/L Nu_hot Nu_cold; half-step predictor flux'
        do i=1,nwall
#ifdef SideHeatedCell
            write(unit,'(3ES25.16E3)') yp(i),NuHotProfile(i),NuColdProfile(i)
#else
            write(unit,'(3ES25.16E3)') xp(i),NuHotProfile(i),NuColdProfile(i)
#endif
        enddo
        close(unit)
        open(newunit=unit,file='Centerline_2DOpenaccDUGKS.dat',status='replace')
        write(unit,'(A)') '# vertical centerline: y/L uL/kappa; horizontal centerline follows'
        do j=1,ny
            write(unit,'(2ES25.16E3)') yp(j),umid(j)
        enddo
        write(unit,'(A)') '# horizontal centerline: x/L vL/kappa'
        do i=1,nx
            write(unit,'(2ES25.16E3)') xp(i),vmid(i)
        enddo
        close(unit)
    end subroutine final_diagnostics

    subroutine output_unsteady_NuRe_postprocess()
        ! 对真实时间戳作分段线性积分；未覆盖完整窗口时只报告覆盖情况。
        integer :: unit,n,k,part
        real(8) :: a,b,t0,t1,z0,z1,za,zb,integral(2),square(2),width,mean(2),sd(2),start,endt
        open(newunit=unit,file='UnsteadyStatistics_2DOpenaccDUGKS.dat',status='replace')
        write(unit,'(A)') '# time weighted NuVolAvg/ReVolAvg statistics; full, first half, second half'
        do part=1,3
            start=unsteadyAverageStartTf; endt=unsteadyAverageEndTf
            if(part==2) endt=0.5d0*(start+endt)
            if(part==3) start=0.5d0*(start+endt)
            integral=0; square=0; width=0
            do n=2,dimensionlessTime
                t0=history(2,n-1); t1=history(2,n)
                a=max(start,t0); b=min(endt,t1)
                if(b<=a.or.t1<=t0) cycle
                do k=1,2
                    z0=history(k+2,n-1); z1=history(k+2,n)
                    za=z0+(z1-z0)*(a-t0)/(t1-t0); zb=z0+(z1-z0)*(b-t0)/(t1-t0)
                    integral(k)=integral(k)+0.5d0*(b-a)*(za+zb)
                    square(k)=square(k)+(b-a)*(za*za+za*zb+zb*zb)/3d0
                enddo
                width=width+b-a
            enddo
            write(unit,*) 'part, requested start/end, covered duration = ',part,start,endt,width
            if(width>0) then
                mean=integral/width; sd=sqrt(max(0d0,square/width-mean*mean))
                write(unit,*) 'Nu_mean Re_mean Nu_temporal_std Re_temporal_std = ',mean,sd
            endif
            if(abs(width-(endt-start))>2*dt/timeUnit) write(unit,*) 'INCOMPLETE_WINDOW: no convergence claim'
        enddo
        close(unit)
    end subroutine output_unsteady_NuRe_postprocess

    subroutine restart_parameters(p)
        real(8), intent(out) :: p(20)
        p=[Rayleigh,Prandtl,Mach,Thot,Tcold,dx,dy,cs2,Courant,dt,tauv,tauc,gBeta,timeUnit, &
            rbInitialPerturbation,dble(steadyCheckIntervalItc),outputSnapshotInterval,reloadFileInterval, &
            outputPltFileInterval,dble(steadyConsecutiveChecks)]
    end subroutine restart_parameters

    subroutine output_ReloadFile()
        integer :: unit,ios
        real(8) :: p(20)
        character(len=160) :: filename
        call check_fields()
        !$acc update self(f,g)
        if(any(.not.ieee_is_finite(f)).or.any(.not.ieee_is_finite(g))) error stop 'Nonfinite transformed distribution'
        reloadFileNum=reloadFileNum+1
        write(filename,'(A,A,I8.8,A)') reloadFilePrefix,'-',reloadFileNum,'.bin'
        call restart_parameters(p)
        open(newunit=unit,file=trim(filename),access='stream',form='unformatted',status='replace',iostat=ios)
        if(ios/=0) error stop 'Cannot write checkpoint'
        write(unit,iostat=ios) reloadMagic,nx,ny,caseCode,flowCode,periodicCode,p
        if(ios/=0) error stop 'Checkpoint header write failed'
        write(unit,iostat=ios) itc,dimensionlessTime,snapshotFileNum,pltFileNum,reloadFileNum, &
            convergedChecks,lastSampleItc,errorU,errorT
        if(ios/=0) error stop 'Checkpoint state write failed'
        write(unit,iostat=ios) f,g,up,vp,Tp,history(:,1:dimensionlessTime)
        if(ios/=0) error stop 'Checkpoint distribution write failed'
        write(unit,iostat=ios) reloadMagic
        if(ios/=0) error stop 'Checkpoint footer write failed'
        close(unit,iostat=ios)
        if(ios/=0) error stop 'Checkpoint close failed'
        ! 先关闭完整编号二进制，再更新小型 latest 指针；旧编号检查点始终保留。
        open(newunit=unit,file=reloadFilePrefix//'-latest.meta',status='replace',action='write',iostat=ios)
        if(ios/=0) error stop 'Cannot write checkpoint metadata'
        write(unit,'(A)') reloadMagic
        write(unit,'(A)') trim(filename)
        write(unit,*) itc,itc*dt/timeUnit
        close(unit)
    end subroutine output_ReloadFile

    subroutine read_reload()
        integer :: unit,ios,nx0,ny0,c0,m0,per0,nstep,nrow,nsnap,nplt,nreload,nconv,last
        real(8) :: p(20),p0(20),eu,et
        character(len=32) :: magic,footer
        character(len=160) :: filename
        call restart_parameters(p)
        open(newunit=unit,file=reloadFilePrefix//'-latest.meta',status='old',action='read',iostat=ios)
        if(ios/=0) error stop 'Missing latest metadata; set loadInitField=0 for a fresh run'
        read(unit,'(A)',iostat=ios) magic
        if(ios/=0.or.magic/=reloadMagic) error stop 'Invalid checkpoint metadata signature'
        read(unit,'(A)',iostat=ios) filename
        close(unit)
        if(ios/=0) error stop 'Invalid checkpoint filename'
        open(newunit=unit,file=trim(filename),status='old',access='stream',form='unformatted',iostat=ios)
        if(ios/=0) error stop 'Missing checkpoint binary'
        read(unit,iostat=ios) magic,nx0,ny0,c0,m0,per0,p0
        if(ios/=0.or.magic/=reloadMagic) error stop 'Invalid checkpoint header'
        if(nx0/=nx.or.ny0/=ny.or.c0/=caseCode.or.m0/=flowCode.or.per0/=periodicCode) &
            error stop 'Incompatible checkpoint mesh, case, or boundary conditions'
        if(any(abs(p-p0)>1d-12*max(abs(p),1d-30))) error stop 'Incompatible checkpoint physical or output parameters'
        read(unit,iostat=ios) nstep,nrow,nsnap,nplt,nreload,nconv,last,eu,et
        if(ios/=0.or.nstep<0.or.nrow<0.or.nrow>historyCapacity) error stop 'Invalid checkpoint counters'
        read(unit,iostat=ios) f,g,up,vp,Tp,history(:,1:nrow)
        if(ios/=0) error stop 'Truncated checkpoint payload'
        read(unit,iostat=ios) footer
        close(unit)
        if(ios/=0.or.footer/=reloadMagic) error stop 'Incomplete checkpoint footer'
        if(any(.not.ieee_is_finite(f)).or.any(.not.ieee_is_finite(g)).or. &
            any(.not.ieee_is_finite(up)).or.any(.not.ieee_is_finite(vp)).or.any(.not.ieee_is_finite(Tp))) &
            error stop 'Checkpoint contains nonfinite fields'
        if(any(.not.ieee_is_finite(history(:,1:nrow)))) error stop 'Checkpoint contains nonfinite history'
        if(nrow>0) then
            if(nint(history(1,nrow))/=last.or.last>nstep) error stop 'Checkpoint history boundary mismatch'
        endif
        itc=nstep; restartItcOffset=nstep; dimensionlessTime=nrow; snapshotFileNum=nsnap; pltFileNum=nplt
        reloadFileNum=nreload; convergedChecks=nconv; lastSampleItc=last; errorU=eu; errorT=et
        ! 由 tilde 精确恢复物理宏观场，不把它误当作原始分布。
        block
            integer :: i,j,k
            real(8) :: mx,my
            do j=1,ny
                do i=1,nx
                    rho(i,j)=sum(f(i,j,:)); T(i,j)=sum(g(i,j,:)); mx=0; my=0
                    do k=0,8
                        mx=mx+ex(k)*f(i,j,k); my=my+ey(k)*f(i,j,k)
                    enddo
                    u(i,j)=mx/rho(i,j); v(i,j)=my/rho(i,j)+0.5d0*dt*gBeta*(T(i,j)-Tref)
                enddo
            enddo
        end block
    end subroutine read_reload

    subroutine reconcile_convergence_history()
        ! 收敛曲线不是推进必需状态；保留检查点以前的有效行，裁掉未提交尾部。
        integer :: oldunit,tempunit,newunit,ios,parseios,step,laststep
        real(8) :: tf,eu,et
        character(len=256) :: line
        logical :: exists
        inquire(file='Convergence_2DOpenaccDUGKS.dat',exist=exists)
        open(newunit=tempunit,status='scratch',action='readwrite',form='formatted')
        write(tempunit,'(A)') '# itc t_ff errorU errorT'
        laststep=-1
        if(exists) then
            open(newunit=oldunit,file='Convergence_2DOpenaccDUGKS.dat',status='old',action='read',iostat=ios)
            if(ios/=0) error stop 'Cannot read existing convergence history'
            do
                read(oldunit,'(A)',iostat=ios) line
                if(ios/=0) exit
                if(len_trim(line)==0) cycle
                if(line(1:1)=='#') cycle
                read(line,*,iostat=parseios) step,tf,eu,et
                if(parseios/=0) cycle ! 断电时未写完整的最后一行不参与恢复
                if(step<=itc.and.step>laststep) then
                    write(tempunit,'(A)') trim(line)
                    laststep=step
                endif
            enddo
            close(oldunit)
        else
            write(*,*) 'Restart: previous convergence curve is missing; a new curve starts from this checkpoint.'
        endif
        rewind(tempunit)
        open(newunit=newunit,file='Convergence_2DOpenaccDUGKS.dat',status='replace',action='write')
        do
            read(tempunit,'(A)',iostat=ios) line
            if(ios/=0) exit
            write(newunit,'(A)') trim(line)
        enddo
        close(newunit)
        close(tempunit)
    end subroutine reconcile_convergence_history

    subroutine exit_data_2d_openacc()
        !$acc wait(1)
        !$acc exit data delete(ex,ey,opposite,omega,f,g,u,v,T,rho)
        !$acc exit data delete(f_barplus,g_barplus,fsx,fsy,gsx,gsy,fluxFx,fluxFy,fluxGx,fluxGy)
        deallocate(f,g,f_barplus,g_barplus,fsx,fsy,gsx,gsy,fluxFx,fluxFy,fluxGx,fluxGy,u,v,T,rho,up,vp,Tp)
    end subroutine exit_data_2d_openacc
end module commondata

program main
#ifdef _OPENACC
    use openacc
#endif
    use commondata
    implicit none
    integer :: nextSampleItc,reloadIntervalItc,pltIntervalItc,unit,clockRate,clockStart,clockEnd
    logical :: converged
#ifdef _OPENACC
    call acc_init(acc_device_default)
#endif
    call initial()
    reloadIntervalItc=max(1,nint(reloadFileInterval*timeUnit/dt))
    pltIntervalItc=max(1,nint(outputPltFileInterval*timeUnit/dt))
    call enter_data_2d_openacc()
    call system_clock(clockStart,clockRate)
    if(itc==0) call calNuRe(.true.)
    converged=.false.
    do while(itc<itc_max)
#ifdef steadyFlow
        if(convergedChecks>=steadyConsecutiveChecks) then
            converged=.true.
            exit
        endif
#endif
        itc=itc+1
        call collision()
        call reconstruction_slopes()
        call face_flux_x()
        call face_flux_y()
        call streaming()
        call macro()
        if(mod(itc,steadyCheckIntervalItc)==0) then
#ifdef steadyFlow
            call check()
#else
            call check_fields()
#endif
        endif
        ! 第k个目标始终由绝对k*Delta_tff换算，不能先舍入步长再累计以免采样漂移。
        nextSampleItc=max(1,nint(dble(dimensionlessTime)*outputSnapshotInterval*timeUnit/dt))
        if(itc>=nextSampleItc) then
            call calNuRe(.true.)
            if(outputSnapshotFile==1) call output_SnapshotFile()
        endif
        if(outputPltFile==1.and.mod(itc,pltIntervalItc)==0) call output_Tecplot()
        if(outputReloadFile==1.and.mod(itc,reloadIntervalItc)==0) call output_ReloadFile()
    enddo
    call final_diagnostics()
    if(outputPltFile==1) call output_Tecplot()
    if(outputSnapshotFile==1) call output_SnapshotFile()
#ifdef unsteadyFlow
    call output_unsteady_NuRe_postprocess()
#endif
    if(outputReloadFile==1) call output_ReloadFile()
    call system_clock(clockEnd)
    open(newunit=unit,file=settingsFile,status='unknown',position='append')
    write(unit,*) 'Final step and t_ff = ',itc,itc*dt/timeUnit
#ifdef steadyFlow
    write(unit,*) 'Steady convergence attained = ',convergedChecks>=steadyConsecutiveChecks
#else
    write(unit,*) 'Unsteady horizon attained = ',itc>=itc_max
#endif
    write(unit,*) 'Wall seconds = ',dble(clockEnd-clockStart)/dble(max(1,clockRate))
    write(unit,*) 'Host compile/smoke checks do not establish GPU performance or benchmark agreement.'
    close(unit)
    write(*,*) 'Completed step ',itc,' t_ff=',itc*dt/timeUnit,' Nu=',NuVolAvg,' Re=',ReVolAvg
    call exit_data_2d_openacc()
end program main
