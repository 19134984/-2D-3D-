!=============================================================
!!! 二维浮力驱动自然对流 OpenACC 静态多级多块网格版本
!!! 来源：2DRBOpenaccMultiblock.F90 的独立副本；原文件保持不变。
!!! 原文件 SHA256：84ca3661eaa2544cca9bb19403b1eaa7d97e82aa4f09ccbbd7cf657de915aa2c
!!! 块内原始 D2Q9 流场、D2Q5 Luo 温度算法及逐矩非平衡缩放直接沿用。
!!! 第1级贴壁且最细；向中心每增加一级，dx、dt均乘2。
!!! 不使用 RR、HRR、魔法参数反算或额外滤波。
!!! 坐标以最细格距为单位：x=xOffset+(i-0.5)*dx，y方向同理。
!!! 同级四个紧凑矩形拼成环带；最粗级只保留中心矩形，不存整块空心数组。
!=============================================================
!   自定义宏，一些选项的开关
!#define steadyFlow
#define unsteadyFlow

!   流动模式宏的选择，两个都开、两个都关都会报错；只有二选一才通过。
#if defined(steadyFlow) && defined(unsteadyFlow)
#error "Choose only one flow mode: steadyFlow or unsteadyFlow"
#endif
#if !defined(steadyFlow) && !defined(unsteadyFlow)
#error "Define one flow mode: steadyFlow or unsteadyFlow"
#endif

!   速度边界，包括水平垂直边界无滑移，还有垂直边界速度周期
#define HorizontalWallsNoslip
#define VerticalWallsNoslip
!#define VerticalWallsPeriodicalU

!   垂直速度边界宏的选择
#if defined(VerticalWallsNoslip) && defined(VerticalWallsPeriodicalU)
#error "Choose only one vertical velocity BC: VerticalWallsNoslip or VerticalWallsPeriodicalU"
#endif
#if !defined(VerticalWallsNoslip) && !defined(VerticalWallsPeriodicalU)
#error "Define one vertical velocity BC: VerticalWallsNoslip or VerticalWallsPeriodicalU"
#endif

!   温度边界(for Rayleigh Benard Cell)，包括水平边界恒温，垂直边界温度不可穿透以及周期
#define RayleighBenardCell
#define HorizontalWallsConstT
#define VerticalWallsAdiabatic
!#define VerticalWallsPeriodicalT

!   温度边界(for Side Heated Cell)，包括水平边界温度不可穿透，垂直边界恒温,侧壁加热加磁场
!#define SideHeatedCell
!#define HorizontalWallsAdiabatic
!#define VerticalWallsConstT
!#define SideHeatedHa
!~~temperature B.C.~~

!   对流算例宏的选择
#if defined(RayleighBenardCell) && defined(SideHeatedCell)
#error "Choose only one convection case: RayleighBenardCell or SideHeatedCell"
#endif
#if !defined(RayleighBenardCell) && !defined(SideHeatedCell)
#error "Define one convection case: RayleighBenardCell or SideHeatedCell"
#endif

! 温度边界各选一种，周期交换后不能再用同一边的恒温或绝热条件覆盖。
! 该宏已启用，结果就是 1；未启用，结果就是 0。
! != 1：选择的数量不是一个，就进入下面的报错。
#if (defined(VerticalWallsConstT) + defined(VerticalWallsAdiabatic) + defined(VerticalWallsPeriodicalT)) != 1
#error "Choose one vertical temperature BC"
#endif
#if (defined(HorizontalWallsConstT) + defined(HorizontalWallsAdiabatic)) != 1
#error "Choose one horizontal temperature BC"
#endif

! 温度算法使用原 D2Q5 TRT 旧算法，平衡矩系数为 paraA。

!   自定义宏结束
!=============================================================

!=============================================================
    module commondata
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        implicit none

        !===============================================================================================
        ! 可直接修改的网格参数；nGridLevels表示总共几种分辨率，包含中心最粗级。
        ! 默认三级：壁面1024级 -> 中间512级 -> 中心256级。
        ! 只改nGridLevels即可切换为2级或4级；1级是均匀网格回归分支。
        ! nx、ny以最细格距计，必须能被2**(nGridLevels-1)整除。
        !===============================================================================================
        integer(kind=4), parameter :: nx=1024, ny=1024
        integer(kind=4), parameter :: nGridLevels=3
        integer(kind=4), parameter :: levelRefineRatio=2    ! 相邻级固定为2倍，避免直接跨越多级
        ! 固定相邻比2，位移等价于2**max(0,nGridLevels-1)，绕过NVHPC 24.3的幂常量折叠问题。
        integer(kind=4), parameter :: coarsestStep=ishft(1,max(0,nGridLevels-1))
        integer(kind=4), parameter :: interfaceSkin=2    ! 沿用原方案的两层接收节点

        ! 第1级的环带厚度按最细格子数计；64层细格作为初始选择，不代表已完成网格独立验证。
        integer(kind=4), parameter :: finestLayerCellsLeftRight=64
        integer(kind=4), parameter :: finestLayerCellsBottomTop=64
        ! 沿用原共址界面约定：默认左右分区线为63.5与959.5，两侧厚度相差一个最细格距。
        ! 这与原代码Right=Left+1的布置一致，使中心宽度可被粗格距整除。
        ! 其余非中心级的环带厚度按各自本级格子数计，不能当作最细格子数。
        ! 默认第2级32格对应64个最细格距；第3级32格对应128个最细格距。
        integer(kind=4), parameter :: innerLayerCellsLeftRight=32
        integer(kind=4), parameter :: innerLayerCellsBottomTop=32
        ! 各级中心剩余宽高还须能被下一级格距整除；初始化时逐级检查，不静默截断。

        !===============================================================================================
        ! 是否从本版本的旧算例重启；旧两级版本的重启文件格式不同，不能混用。
        !===============================================================================================
        integer(kind=4), parameter :: loadInitField=0
        character(len=256) :: reloadFileName=''
        integer(kind=4) :: reloadItcOffset=0

        !===============================================================================================
        ! 最细网格保持原始Luo参数；其余各级仅用原始粗细网格关系换算。
        !===============================================================================================
#ifdef SideHeatedCell
        real(kind=8), parameter :: lengthUnit=dble(nx)
#else
        real(kind=8), parameter :: lengthUnit=dble(ny)
#endif
        real(kind=8), parameter :: pi=acos(-1.0d0)
        real(kind=8), parameter :: Rayleigh=1.0d7, Prandtl=0.7d0, Mach=0.1d0
        real(kind=8), parameter :: Thot=0.5d0, Tcold=-0.5d0, Tref=0.5d0*(Thot+Tcold)
        real(kind=8), parameter :: tauf=0.5d0+Mach*lengthUnit*sqrt(3.0d0*Prandtl/Rayleigh)
        real(kind=8), parameter :: viscosity=(tauf-0.5d0)/3.0d0, diffusivity=viscosity/Prandtl
        real(kind=8), parameter :: paraA=20.0d0*sqrt(3.0d0)*diffusivity-4.0d0
        real(kind=8), parameter :: Snu=1.0d0/tauf, Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)
        real(kind=8), parameter :: Qk=3.0d0-sqrt(3.0d0), Qnu=4.0d0*sqrt(3.0d0)-6.0d0
        real(kind=8), parameter :: thermalGeqCoeff=10.0d0/(4.0d0+paraA)
        real(kind=8), parameter :: gBeta1=Rayleigh*viscosity*diffusivity/lengthUnit
        real(kind=8), parameter :: gBeta=gBeta1/lengthUnit/lengthUnit
        real(kind=8), parameter :: timeUnit=sqrt(lengthUnit/gBeta), velocityUnit=sqrt(gBeta*lengthUnit)
#ifdef SideHeatedHa
        real(kind=8), parameter :: Ha=20.0d0, phi=0.0d0*pi/180.0d0
        real(kind=8), parameter :: B2sigemarho=Ha**2*viscosity/lengthUnit**2
#endif

        !===============================================================================================
        ! 输出及停止参数；步数全部对齐到最粗级的同步时刻。
        !===============================================================================================
        real(kind=8), parameter :: epsU=1.0d-7, epsT=1.0d-7
#ifdef steadyFlow
        real(kind=8), parameter :: outputSnapshotInterval=10.0d0
        integer(kind=4), parameter :: itc_max=ceiling(20000000.0d0/dble(coarsestStep))*coarsestStep
#else
        real(kind=8), parameter :: outputSnapshotInterval=0.5d0, unsteadyRunDuration=1000.0d0
#endif
        real(kind=8), parameter :: reloadFileInterval=100.0d0, outputPltFileInterval=100.0d0
        integer(kind=4), parameter :: outputSnapshotFile=1, outputPltFile=1, outputReloadFile=1
        integer(kind=4), parameter :: steadyCheckIntervalItc=ceiling(2000.0d0/dble(coarsestStep))*coarsestStep
        integer(kind=4), parameter :: outputSnapshotIntervalItc= &
            max(1,ceiling(outputSnapshotInterval*timeUnit/dble(coarsestStep)))*coarsestStep
        integer(kind=4), parameter :: reloadFileIntervalItc= &
            max(1,ceiling(reloadFileInterval*timeUnit/dble(coarsestStep)))*coarsestStep
        integer(kind=4), parameter :: outputPltFileIntervalItc= &
            max(1,ceiling(outputPltFileInterval*timeUnit/dble(coarsestStep)))*coarsestStep
#ifdef unsteadyFlow
        integer(kind=4), parameter :: itc_max=ceiling(unsteadyRunDuration*timeUnit/ &
            dble(outputSnapshotIntervalItc))*outputSnapshotIntervalItc
#endif
        character(len=*), parameter :: settingsFile='SimulationSettings2DOpenaccMultilevel.txt'
        character(len=*), parameter :: NuReHistoryFile='NuRe_2DOpenaccMultiblock.dat'
        character(len=*), parameter :: convergenceFile='Convergence_2DOpenaccMultiblock.dat'
        character(len=*), parameter :: snapshotFilePrefix='buoyancyCavity2DOpenaccMultiblockSnapshot'
        character(len=*), parameter :: pltFolderPrefix='buoyancyCavity2DOpenaccMultiblockTecplot'
        character(len=*), parameter :: reloadFilePrefix='reloadFile2DOpenaccMultilevel'

        integer(kind=4) :: ex(0:8)=[0,1,0,-1,0,1,-1,-1,1], ey(0:8)=[0,0,1,0,-1,1,1,-1,-1]
        real(kind=8) :: omega(0:8), omegaT(0:4)
        integer(kind=4) :: itc=0, NuReSampleCount=0, nBlocks=0
        real(kind=8) :: errorU=100.0d0, errorT=100.0d0
        real(kind=8) :: SnuLevel(nGridLevels), SqLevel(nGridLevels), QkLevel(nGridLevels), QnuLevel(nGridLevels)

        !===============================================================================================
        ! 下列类型只在CPU上组织数组；GPU核函数仍显式传入平直数组，避免设备深复制。
        !===============================================================================================
        type grid_block
            integer :: level, ni, nj, globalI, globalJ
            real(kind=8) :: dx, xOffset, yOffset
            logical :: wall(4)
            real(kind=8), allocatable :: f(:,:,:),g(:,:,:),f_post(:,:,:),g_post(:,:,:)
            real(kind=8), allocatable :: rho(:,:),u(:,:),v(:,:),T(:,:),Fx(:,:),Fy(:,:),area(:,:)
#ifdef steadyFlow
            real(kind=8), allocatable :: up(:,:),vp(:,:),Tp(:,:)
#endif
        end type grid_block
        type grid_level
            integer :: firstBlock,lastBlock,ni,nj
            real(kind=8) :: dx,xOffset,yOffset,outerBox(4),innerBox(4)
        end type grid_level
        type node_list
            integer, allocatable :: map(:,:)    ! 第1、2行为本地i,j，第3行为缓冲区节点编号
        end type node_list
        type level_interface
            integer :: nSource=0,nFine=0,nCoarse=0,cycles=0
            integer, allocatable :: stencil(:,:)
            real(kind=8), allocatable :: weight(:,:),packet(:,:),history(:,:,:),restriction(:,:)
            type(node_list), allocatable :: coarseSource(:),fineTarget(:),fineSource(:),coarseTarget(:)
        end type level_interface
        type halo_link
            integer :: source,target,periodic
            integer, allocatable :: map(:,:)    ! 目标i,j和来源i,j；同级交换不插值
        end type halo_link
        type(grid_block), allocatable :: block(:)
        type(grid_level), allocatable :: grid(:)
        type(level_interface), allocatable :: coupling(:)
        type(halo_link), allocatable :: halo(:)
    end module commondata

!===================================================================================================
! 多级管理：只扩展网格组织、接口与同步；块内碰撞、迁移和边界子程序保留在文件末尾。
!===================================================================================================
    module multilevel_solver
        use commondata
        implicit none
    contains

!===================================================================================================
! 初始化：逐级检查整除与环带厚度，生成紧凑网格，并建立不重复积分面积。
!===================================================================================================
    subroutine initial_geometry()
        integer :: l,b,i,j,tx,ty,ix1,ix2,iy1,iy2,nxLocal,nyLocal
        real(kind=8) :: dx,x,y,totalArea,widthX,widthY
        if(nGridLevels<1) error stop 'nGridLevels must be positive'
        if(levelRefineRatio/=2) error stop 'This nested hierarchy uses adjacent refinement ratio two'
        if(min(nx,ny)<4) error stop 'nx,ny must be at least four'
        if(mod(nx,coarsestStep)/=0.or.mod(ny,coarsestStep)/=0) &
            error stop 'nx,ny must be divisible by 2**(nGridLevels-1)'
        if(interfaceSkin<2) error stop 'At least two interface receiver layers are required'
        if(loadInitField/=0.and.loadInitField/=1) error stop 'loadInitField must be 0 or 1'
        if(paraA<=-4.0d0.or.paraA>=1.0d0) error stop 'Finest Luo paraA must be in (-4,1)'
        if(min(outputSnapshotInterval,reloadFileInterval,outputPltFileInterval)<=0.0d0) &
            error stop 'Output intervals must be positive'
        nBlocks=4*(nGridLevels-1)+1
        allocate(block(nBlocks),grid(nGridLevels),coupling(nGridLevels-1))
        omega(0)=4.0d0/9.0d0
        omega(1:4)=1.0d0/9.0d0
        omega(5:8)=1.0d0/36.0d0
        omegaT(0)=(1.0d0-paraA)/5.0d0
        omegaT(1:4)=(paraA+4.0d0)/20.0d0
        b=0
        do l=1,nGridLevels
            dx=dble(levelRefineRatio**(l-1))
            grid(l)%dx=dx
            ! s_l=1/[1/2+(1/s_f-1/2)/dx_l]；paraA和所有平衡态系数保持最细级取值。
            SnuLevel(l)=1.0d0/(0.5d0+(1.0d0/Snu-0.5d0)/dx)
            SqLevel(l)=1.0d0/(0.5d0+(1.0d0/Sq-0.5d0)/dx)
            QkLevel(l)=1.0d0/(0.5d0+(1.0d0/Qk-0.5d0)/dx)
            QnuLevel(l)=1.0d0/(0.5d0+(1.0d0/Qnu-0.5d0)/dx)
            if(l==1) then
                SnuLevel(l)=Snu
                SqLevel(l)=Sq
                QkLevel(l)=Qk
                QnuLevel(l)=Qnu
            endif
            if(l==1) then
                grid(l)%outerBox=[0.0d0,dble(nx),0.0d0,dble(ny)]
                grid(l)%xOffset=0.0d0
                grid(l)%yOffset=0.0d0
                nxLocal=nx
                nyLocal=ny
            else
                grid(l)%outerBox=grid(l-1)%innerBox
                widthX=grid(l)%outerBox(2)-grid(l)%outerBox(1)
                widthY=grid(l)%outerBox(4)-grid(l)%outerBox(3)
                if(mod(nint(widthX),nint(dx))/=0.or.mod(nint(widthY),nint(dx))/=0) &
                    error stop 'A level width/height is not divisible by its spacing'
                nxLocal=nint(widthX/dx)+1+2*interfaceSkin
                nyLocal=nint(widthY/dx)+1+2*interfaceSkin
                grid(l)%xOffset=grid(l)%outerBox(1)-(interfaceSkin+0.5d0)*dx
                grid(l)%yOffset=grid(l)%outerBox(3)-(interfaceSkin+0.5d0)*dx
                if(min(nxLocal,nyLocal)-2*interfaceSkin<4) &
                    error stop 'Central level needs at least four valid nodes per direction'
            endif
            grid(l)%ni=nxLocal
            grid(l)%nj=nyLocal
            grid(l)%firstBlock=b+1
            if(l==nGridLevels) then
                grid(l)%innerBox=0.0d0
                call add_block(b,l,1,nxLocal,1,nyLocal)
            else
                if(l==1) then
                    tx=finestLayerCellsLeftRight
                    ty=finestLayerCellsBottomTop
                    ix1=tx
                    ix2=nxLocal-tx
                    iy1=ty
                    iy2=nyLocal-ty
                else
                    tx=innerLayerCellsLeftRight
                    ty=innerLayerCellsBottomTop
                    ix1=interfaceSkin+1+tx
                    ix2=nxLocal-interfaceSkin-tx
                    iy1=interfaceSkin+1+ty
                    iy2=nyLocal-interfaceSkin-ty
                endif
                if(min(tx,ty)<2*interfaceSkin+4) &
                    error stop 'Ring thickness must be at least 2*interfaceSkin+4 local cells'
                if(min(ix2-ix1,iy2-iy1)<levelRefineRatio*3) &
                    error stop 'Too many levels or rings too thick: no usable central region'
                if(mod(ix2-ix1,levelRefineRatio)/=0.or.mod(iy2-iy1,levelRefineRatio)/=0) &
                    error stop 'Inner width/height must be divisible by the adjacent refinement ratio'
                grid(l)%innerBox=[grid(l)%xOffset+(ix1-0.5d0)*dx, &
                    grid(l)%xOffset+(ix2-0.5d0)*dx,grid(l)%yOffset+(iy1-0.5d0)*dx, &
                    grid(l)%yOffset+(iy2-0.5d0)*dx]
                ! 左右条带贯穿本级存储高度；上下条带只存中间缺口，角点不重复存储。
                call add_block(b,l,1,ix1+interfaceSkin,1,nyLocal)
                call add_block(b,l,ix2-interfaceSkin,nxLocal,1,nyLocal)
                call add_block(b,l,ix1+interfaceSkin+1,ix2-interfaceSkin-1,1,iy1+interfaceSkin)
                call add_block(b,l,ix1+interfaceSkin+1,ix2-interfaceSkin-1,iy2-interfaceSkin,nyLocal)
            endif
            grid(l)%lastBlock=b
        enddo
        totalArea=0.0d0
        do b=1,nBlocks
            call allocate_block(block(b))
            do j=1,block(b)%nj
                y=block(b)%yOffset+(j-0.5d0)*block(b)%dx
                do i=1,block(b)%ni
                    x=block(b)%xOffset+(i-0.5d0)*block(b)%dx
                    block(b)%area(i,j)=owned_cell_area(block(b)%level,x,y)
                enddo
            enddo
            totalArea=totalArea+sum(block(b)%area)
        enddo
        if(abs(totalArea-dble(nx)*dble(ny))>1.0d-10*dble(nx)*dble(ny)) &
            error stop 'Owned integration areas do not cover the physical domain exactly once'
        do l=1,nGridLevels-1
            call build_interface(l)
        enddo
        call build_halos()
    end subroutine initial_geometry

    subroutine add_block(b,l,i1,i2,j1,j2)
        integer,intent(inout) :: b
        integer,intent(in) :: l,i1,i2,j1,j2
        b=b+1
        if(min(i2-i1+1,j2-j1+1)<1) error stop 'Nonpositive compact block dimension'
        block(b)%level=l
        block(b)%globalI=i1
        block(b)%globalJ=j1
        block(b)%ni=i2-i1+1
        block(b)%nj=j2-j1+1
        block(b)%dx=grid(l)%dx
        block(b)%xOffset=grid(l)%xOffset+(i1-1)*grid(l)%dx
        block(b)%yOffset=grid(l)%yOffset+(j1-1)*grid(l)%dx
        block(b)%wall=.false.
        if(l==1) block(b)%wall=[i1==1,i2==grid(l)%ni,j1==1,j2==grid(l)%nj]
    end subroutine add_block

    subroutine allocate_block(b)
        type(grid_block),intent(inout) :: b
        integer :: ni,nj
        ni=b%ni
        nj=b%nj
        allocate(b%f(ni,nj,0:8),b%g(ni,nj,0:4))
        allocate(b%f_post(0:ni+1,0:nj+1,0:8),b%g_post(0:ni+1,0:nj+1,0:4))
        allocate(b%rho(ni,nj),b%u(ni,nj),b%v(ni,nj),b%T(ni,nj),b%Fx(ni,nj),b%Fy(ni,nj),b%area(ni,nj))
#ifdef steadyFlow
        allocate(b%up(ni,nj),b%vp(ni,nj),b%Tp(ni,nj))
#endif
    end subroutine allocate_block

    real(kind=8) function rectangle_area(rect,box) result(area)
        real(kind=8),intent(in) :: rect(4),box(4)
        area=max(0.0d0,min(rect(2),box(2))-max(rect(1),box(1)))* &
            max(0.0d0,min(rect(4),box(4))-max(rect(3),box(3)))
    end function rectangle_area

    real(kind=8) function owned_cell_area(l,x,y) result(area)
        integer,intent(in) :: l
        real(kind=8),intent(in) :: x,y
        real(kind=8) :: d,rect(4)
        d=0.5d0*grid(l)%dx
        rect=[x-d,x+d,y-d,y+d]
        area=rectangle_area(rect,grid(l)%outerBox)
        if(l<nGridLevels) area=area-rectangle_area(rect,grid(l)%innerBox)
        area=max(0.0d0,area)
    end function owned_cell_area

    logical function inside_box(x,y,box)
        real(kind=8),intent(in) :: x,y,box(4)
        inside_box=x>box(1)+1.0d-10.and.x<box(2)-1.0d-10.and. &
            y>box(3)+1.0d-10.and.y<box(4)-1.0d-10
    end function inside_box

    logical function outside_box(x,y,box)
        real(kind=8),intent(in) :: x,y,box(4)
        outside_box=x<box(1)-1.0d-10.or.x>box(2)+1.0d-10.or. &
            y<box(3)-1.0d-10.or.y>box(4)+1.0d-10
    end function outside_box

    subroutine locate_level_node(l,gi,gj,b,i,j,owned)
        integer,intent(in) :: l,gi,gj
        integer,intent(out) :: b,i,j
        logical,intent(in) :: owned
        integer :: k
        b=0
        i=0
        j=0
        do k=grid(l)%firstBlock,grid(l)%lastBlock
            if(gi<block(k)%globalI.or.gi>=block(k)%globalI+block(k)%ni) cycle
            if(gj<block(k)%globalJ.or.gj>=block(k)%globalJ+block(k)%nj) cycle
            b=k
            i=gi-block(k)%globalI+1
            j=gj-block(k)%globalJ+1
            if(owned) then
                if(block(k)%area(i,j)<=0.0d0) error stop 'Interface stencil enters a receiver: increase ring thickness'
            endif
            return
        enddo
    end subroutine locate_level_node

!===================================================================================================
! 接口只在初始化时定位：粗到细仍为共址取值或4x4拉格朗日插值，细到粗仍为共址取值。
! 模板可以跨同级条带，不能把存储拼接线当作物理边界或另加单侧插值。
!===================================================================================================
    subroutine build_interface(l)
        integer,intent(in) :: l
        integer :: b,i,j,p,q,a,ix,iy,firstX,firstY,width,sb,si,sj,gi,gj,cap,ns,nf,nc
        integer,allocatable :: fb(:),fi(:),fj(:),cb(:),ci(:),cj(:),sourceMap(:,:)
        integer,allocatable :: sourceB(:),sourceI(:),sourceJ(:),donorB(:),donorI(:),donorJ(:)
        real(kind=8) :: x,y,qx,qy,wx(4),wy(4)
        nf=0
        nc=0
        do b=grid(l)%firstBlock,grid(l)%lastBlock
            do j=1,block(b)%nj
                y=block(b)%yOffset+(j-0.5d0)*grid(l)%dx
                do i=1,block(b)%ni
                    x=block(b)%xOffset+(i-0.5d0)*grid(l)%dx
                    if(inside_box(x,y,grid(l)%innerBox)) nf=nf+1
                enddo
            enddo
        enddo
        do b=grid(l+1)%firstBlock,grid(l+1)%lastBlock
            do j=1,block(b)%nj
                y=block(b)%yOffset+(j-0.5d0)*grid(l+1)%dx
                do i=1,block(b)%ni
                    x=block(b)%xOffset+(i-0.5d0)*grid(l+1)%dx
                    if(outside_box(x,y,grid(l+1)%outerBox)) nc=nc+1
                enddo
            enddo
        enddo
        coupling(l)%nFine=nf
        coupling(l)%nCoarse=nc
        allocate(fb(nf),fi(nf),fj(nf),cb(nc),ci(nc),cj(nc),donorB(nc),donorI(nc),donorJ(nc))
        allocate(coupling(l)%stencil(16,nf),coupling(l)%weight(16,nf))
        coupling(l)%stencil=1
        coupling(l)%weight=0.0d0
        allocate(sourceMap(grid(l+1)%ni,grid(l+1)%nj))
        sourceMap=0
        cap=16*nf
        allocate(sourceB(cap),sourceI(cap),sourceJ(cap))
        ns=0
        p=0
        do b=grid(l)%firstBlock,grid(l)%lastBlock
            do j=1,block(b)%nj
                y=block(b)%yOffset+(j-0.5d0)*grid(l)%dx
                do i=1,block(b)%ni
                    x=block(b)%xOffset+(i-0.5d0)*grid(l)%dx
                    if(.not.inside_box(x,y,grid(l)%innerBox)) cycle
                    p=p+1
                    fb(p)=b
                    fi(p)=i
                    fj(p)=j
                    qx=(x-grid(l+1)%xOffset)/grid(l+1)%dx+0.5d0
                    qy=(y-grid(l+1)%yOffset)/grid(l+1)%dx+0.5d0
                    call coarse_source_stencil(qx,qy,grid(l+1)%ni,grid(l+1)%nj,interfaceSkin, &
                        firstX,firstY,width,wx,wy)
                    a=0
                    do iy=1,width
                        do ix=1,width
                            a=a+1
                            gi=firstX+ix-1
                            gj=firstY+iy-1
                            if(gi<1.or.gi>size(sourceMap,1).or.gj<1.or.gj>size(sourceMap,2)) &
                                error stop 'Coarse interpolation stencil outside allocated level'
                            q=sourceMap(gi,gj)
                            if(q==0) then
                                call locate_level_node(l+1,gi,gj,sb,si,sj,.true.)
                                if(sb==0) error stop 'Coarse interpolation stencil crosses the compact storage hole'
                                ns=ns+1
                                q=ns
                                sourceMap(gi,gj)=q
                                sourceB(q)=sb
                                sourceI(q)=si
                                sourceJ(q)=sj
                            endif
                            coupling(l)%stencil(a,p)=q
                            coupling(l)%weight(a,p)=wx(ix)*wy(iy)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        p=0
        do b=grid(l+1)%firstBlock,grid(l+1)%lastBlock
            do j=1,block(b)%nj
                y=block(b)%yOffset+(j-0.5d0)*grid(l+1)%dx
                do i=1,block(b)%ni
                    x=block(b)%xOffset+(i-0.5d0)*grid(l+1)%dx
                    if(.not.outside_box(x,y,grid(l+1)%outerBox)) cycle
                    p=p+1
                    cb(p)=b
                    ci(p)=i
                    cj(p)=j
                    qx=(x-grid(l)%xOffset)/grid(l)%dx+0.5d0
                    qy=(y-grid(l)%yOffset)/grid(l)%dx+0.5d0
                    gi=nint(qx)
                    gj=nint(qy)
                    if(abs(qx-gi)>1.0d-10.or.abs(qy-gj)>1.0d-10) &
                        error stop 'Fine-to-coarse receiver is not colocated with the fine lattice'
                    call locate_level_node(l,gi,gj,sb,si,sj,.true.)
                    if(sb==0) error stop 'Fine donor missing for a coarse receiver'
                    donorB(p)=sb
                    donorI(p)=si
                    donorJ(p)=sj
                enddo
            enddo
        enddo
        coupling(l)%nSource=ns
        allocate(coupling(l)%packet(0:19,ns),coupling(l)%history(0:19,nf,0:2),coupling(l)%restriction(0:19,nc))
        coupling(l)%packet=0.0d0
        coupling(l)%history=0.0d0
        coupling(l)%restriction=0.0d0
        call group_nodes(ns,sourceB(:ns),sourceI(:ns),sourceJ(:ns),coupling(l)%coarseSource)
        call group_nodes(nf,fb,fi,fj,coupling(l)%fineTarget)
        call group_nodes(nc,donorB,donorI,donorJ,coupling(l)%fineSource)
        call group_nodes(nc,cb,ci,cj,coupling(l)%coarseTarget)
    end subroutine build_interface

    subroutine group_nodes(n,owner,ii,jj,lists)
        integer,intent(in) :: n,owner(n),ii(n),jj(n)
        type(node_list),allocatable,intent(out) :: lists(:)
        integer :: count(nBlocks),b,p,k
        allocate(lists(nBlocks))
        count=0
        do p=1,n
            count(owner(p))=count(owner(p))+1
        enddo
        do b=1,nBlocks
            allocate(lists(b)%map(3,count(b)))
        enddo
        count=0
        do p=1,n
            b=owner(p)
            count(b)=count(b)+1
            k=count(b)
            lists(b)%map(:,k)=[ii(p),jj(p),p]
        enddo
    end subroutine group_nodes

!===================================================================================================
! 同级碰撞后迁移外圈：预先记录共址节点，运行时只复制f_post或g_post。
!===================================================================================================
    subroutine build_halos()
        integer :: cap,b,i,j,gi,gj,sb,si,sj,n,p,k,mode,ng
        integer,allocatable :: src(:),dst(:),kind(:),indices(:,:),counts(:,:,:),lookup(:,:,:),filled(:)
        logical :: periodicAny
        periodicAny=.false.
#if defined(VerticalWallsPeriodicalU) || defined(VerticalWallsPeriodicalT)
        periodicAny=.true.
#endif
        cap=0
        do b=1,nBlocks
            cap=cap+2*block(b)%ni+2*block(b)%nj+4
        enddo
        allocate(src(cap),dst(cap),kind(cap),indices(4,cap),counts(0:1,nBlocks,nBlocks),lookup(0:1,nBlocks,nBlocks))
        counts=0
        lookup=0
        n=0
        do b=1,nBlocks
            do j=0,block(b)%nj+1
                do i=0,block(b)%ni+1
                    if(i/=0.and.i/=block(b)%ni+1.and.j/=0.and.j/=block(b)%nj+1) cycle
                    gi=block(b)%globalI+i-1
                    gj=block(b)%globalJ+j-1
                    mode=0
                    if(block(b)%level==1.and.periodicAny) then
                        if(gi==0.or.gi==nx+1) then
                            gi=modulo(gi-1,nx)+1
                            mode=1
                        endif
                    endif
                    call locate_level_node(block(b)%level,gi,gj,sb,si,sj,.false.)
                    if(sb==0) cycle
                    n=n+1
                    src(n)=sb
                    dst(n)=b
                    kind(n)=mode
                    indices(:,n)=[i,j,si,sj]
                    counts(mode,sb,b)=counts(mode,sb,b)+1
                enddo
            enddo
        enddo
        ng=count(counts>0)
        allocate(halo(ng),filled(ng))
        filled=0
        k=0
        do b=1,nBlocks
            do sb=1,nBlocks
                do mode=0,1
                    if(counts(mode,sb,b)==0) cycle
                    k=k+1
                    lookup(mode,sb,b)=k
                    halo(k)%source=sb
                    halo(k)%target=b
                    halo(k)%periodic=mode
                    allocate(halo(k)%map(4,counts(mode,sb,b)))
                enddo
            enddo
        enddo
        do p=1,n
            k=lookup(kind(p),src(p),dst(p))
            filled(k)=filled(k)+1
            halo(k)%map(:,filled(k))=indices(:,p)
        enddo
    end subroutine build_halos

!===================================================================================================
! 初始场和OpenACC数据寿命。所有派生类型均留在主机，核函数只接收显式形状的普通数组。
!===================================================================================================
    subroutine initial_fields()
        integer :: b,i,j,a
        real(kind=8) :: x,y
        do b=1,nBlocks
            block(b)%rho=1.0d0
            block(b)%u=0.0d0
            block(b)%v=0.0d0
            block(b)%Fx=0.0d0
            block(b)%Fy=0.0d0
            block(b)%f_post=0.0d0
            block(b)%g_post=0.0d0
            do j=1,block(b)%nj
                y=block(b)%yOffset+(j-0.5d0)*block(b)%dx
                do i=1,block(b)%ni
                    x=block(b)%xOffset+(i-0.5d0)*block(b)%dx
#ifdef SideHeatedCell
                    block(b)%T(i,j)=Thot+x/dble(nx)*(Tcold-Thot)
#else
                    block(b)%T(i,j)=Thot+y/dble(ny)*(Tcold-Thot)
#endif
                    do a=0,8
                        block(b)%f(i,j,a)=omega(a)
                    enddo
                    do a=0,4
                        block(b)%g(i,j,a)=omegaT(a)*block(b)%T(i,j)
                    enddo
                enddo
            enddo
#ifdef steadyFlow
            block(b)%up=block(b)%u
            block(b)%vp=block(b)%v
            block(b)%Tp=block(b)%T
#endif
        enddo
    end subroutine initial_fields

    subroutine enter_data_2d_openacc()
        integer :: b,l,k
        !$acc enter data copyin(ex,ey,omega,omegaT)
        do b=1,nBlocks
            call grid_data(block(b)%ni,block(b)%nj,1,block(b)%f,block(b)%g,block(b)%f_post,block(b)%g_post, &
                block(b)%rho,block(b)%u,block(b)%v,block(b)%T,block(b)%Fx,block(b)%Fy)
        enddo
        do l=1,nGridLevels-1
            call interface_data(l,1)
            do b=1,nBlocks
                call list_data(3,size(coupling(l)%coarseSource(b)%map,2),coupling(l)%coarseSource(b)%map,1)
                call list_data(3,size(coupling(l)%fineTarget(b)%map,2),coupling(l)%fineTarget(b)%map,1)
                call list_data(3,size(coupling(l)%fineSource(b)%map,2),coupling(l)%fineSource(b)%map,1)
                call list_data(3,size(coupling(l)%coarseTarget(b)%map,2),coupling(l)%coarseTarget(b)%map,1)
            enddo
        enddo
        do k=1,size(halo)
            call list_data(4,size(halo(k)%map,2),halo(k)%map,1)
        enddo
    end subroutine enter_data_2d_openacc

    subroutine grid_data(ni,nj,action,f,g,fp,gp,rho,u,v,T,Fx,Fy)
        integer,intent(in) :: ni,nj,action
        real(kind=8) :: f(ni,nj,0:8),g(ni,nj,0:4),fp(0:ni+1,0:nj+1,0:8),gp(0:ni+1,0:nj+1,0:4)
        real(kind=8) :: rho(ni,nj),u(ni,nj),v(ni,nj),T(ni,nj),Fx(ni,nj),Fy(ni,nj)
        if(action==1) then
            !$acc enter data copyin(f,g,fp,gp,rho,u,v,T,Fx,Fy)
        elseif(action==2) then
            !$acc update self(rho,u,v,T) async(1)
        elseif(action==3) then
            !$acc update self(f,g,rho,u,v,T,Fx,Fy) async(1)
        else
            !$acc exit data delete(f,g,fp,gp,rho,u,v,T,Fx,Fy)
        endif
    end subroutine grid_data

    subroutine list_data(rows,n,map,action)
        integer,intent(in) :: rows,n,action
        integer :: map(rows,n)
        if(n==0) return
        if(action==1) then
            !$acc enter data copyin(map)
        else
            !$acc exit data delete(map)
        endif
    end subroutine list_data

    subroutine interface_data(l,action)
        integer,intent(in) :: l,action
        call interface_arrays_data(coupling(l)%nSource,coupling(l)%nFine,coupling(l)%nCoarse,action, &
            coupling(l)%stencil,coupling(l)%weight,coupling(l)%packet,coupling(l)%history,coupling(l)%restriction)
    end subroutine interface_data

    subroutine interface_arrays_data(ns,nf,nc,action,stencil,weight,packet,history,restriction)
        integer,intent(in) :: ns,nf,nc,action
        integer :: stencil(16,nf)
        real(kind=8) :: weight(16,nf),packet(0:19,ns),history(0:19,nf,0:2),restriction(0:19,nc)
        if(action==1) then
            !$acc enter data copyin(stencil,weight,packet,history,restriction)
        elseif(action==2) then
            !$acc update self(history) async(1)
        else
            !$acc exit data delete(stencil,weight,packet,history,restriction)
        endif
    end subroutine interface_arrays_data

    subroutine update_host_all(full)
        logical,intent(in) :: full
        integer :: b,action
        action=2
        if(full) action=3
        do b=1,nBlocks
            call grid_data(block(b)%ni,block(b)%nj,action,block(b)%f,block(b)%g,block(b)%f_post,block(b)%g_post, &
                block(b)%rho,block(b)%u,block(b)%v,block(b)%T,block(b)%Fx,block(b)%Fy)
        enddo
        !$acc wait(1)
    end subroutine update_host_all

    subroutine exit_data_2d_openacc()
        integer :: b,l,k
        !$acc wait(1)
        do k=1,size(halo)
            call list_data(4,size(halo(k)%map,2),halo(k)%map,-1)
        enddo
        do l=1,nGridLevels-1
            do b=1,nBlocks
                call list_data(3,size(coupling(l)%coarseSource(b)%map,2),coupling(l)%coarseSource(b)%map,-1)
                call list_data(3,size(coupling(l)%fineTarget(b)%map,2),coupling(l)%fineTarget(b)%map,-1)
                call list_data(3,size(coupling(l)%fineSource(b)%map,2),coupling(l)%fineSource(b)%map,-1)
                call list_data(3,size(coupling(l)%coarseTarget(b)%map,2),coupling(l)%coarseTarget(b)%map,-1)
            enddo
            call interface_data(l,-1)
        enddo
        do b=1,nBlocks
            call grid_data(block(b)%ni,block(b)%nj,-1,block(b)%f,block(b)%g,block(b)%f_post,block(b)%g_post, &
                block(b)%rho,block(b)%u,block(b)%v,block(b)%T,block(b)%Fx,block(b)%Fy)
        enddo
        !$acc exit data delete(ex,ey,omega,omegaT)
        deallocate(block,grid,coupling,halo)
    end subroutine exit_data_2d_openacc

!===================================================================================================
! 交换量仍为原来的20项，含单位时间力与逐矩缩放；不增加正则化或滤波。
!===================================================================================================
    subroutine pack_nodes(ni,nj,n,np,map,packet,dt,sn,sq_,qk_,qn,f,g,rho,u,v,T,Fx,Fy)
        integer,intent(in) :: ni,nj,n,np,map(3,n)
        real(kind=8),intent(in) :: dt,sn,sq_,qk_,qn
        real(kind=8),intent(in) :: f(ni,nj,0:8),g(ni,nj,0:4),rho(ni,nj),u(ni,nj),v(ni,nj),T(ni,nj),Fx(ni,nj),Fy(ni,nj)
        real(kind=8),intent(inout) :: packet(0:19,np)
        integer :: p,i,j,k,a
        real(kind=8) :: fv(0:8),gv(0:4),buffer(0:19)
        if(n==0) return
        !$acc parallel loop gang vector present(map,packet,f,g,rho,u,v,T,Fx,Fy) async(1) &
        !$acc& firstprivate(dt,sn,sq_,qk_,qn) private(i,j,k,a,fv,gv,buffer)
        do p=1,n
            i=map(1,p)
            j=map(2,p)
            k=map(3,p)
            do a=0,8
                fv(a)=f(i,j,a)
            enddo
            do a=0,4
                gv(a)=g(i,j,a)
            enddo
            call pack_interface_values(dt,sn,sq_,qk_,qn,fv,gv,rho(i,j),u(i,j),v(i,j),T(i,j),Fx(i,j),Fy(i,j),buffer)
            do a=0,19
                packet(a,k)=buffer(a)
            enddo
        enddo
    end subroutine pack_nodes

    subroutine save_coarse_history(l,slot)
        integer,intent(in) :: l,slot
        integer :: b,lev
        lev=l+1
        do b=grid(lev)%firstBlock,grid(lev)%lastBlock
            call pack_nodes(block(b)%ni,block(b)%nj,size(coupling(l)%coarseSource(b)%map,2),coupling(l)%nSource, &
                coupling(l)%coarseSource(b)%map,coupling(l)%packet,grid(lev)%dx,SnuLevel(lev),SqLevel(lev), &
                QkLevel(lev),QnuLevel(lev),block(b)%f,block(b)%g,block(b)%rho,block(b)%u,block(b)%v,block(b)%T, &
                block(b)%Fx,block(b)%Fy)
        enddo
        call spatial_interpolation(coupling(l)%nSource,coupling(l)%nFine,slot,coupling(l)%stencil, &
            coupling(l)%weight,coupling(l)%packet,coupling(l)%history)
    end subroutine save_coarse_history

    subroutine spatial_interpolation(ns,nf,slot,stencil,weight,packet,history)
        integer,intent(in) :: ns,nf,slot,stencil(16,nf)
        real(kind=8),intent(in) :: weight(16,nf),packet(0:19,ns)
        real(kind=8),intent(inout) :: history(0:19,nf,0:2)
        integer :: p,a,k
        real(kind=8) :: value
        !$acc parallel loop collapse(2) present(stencil,weight,packet,history) firstprivate(slot) async(1) private(k,value)
        do p=1,nf
            do a=0,19
                value=0.0d0
                do k=1,16
                    value=value+weight(k,p)*packet(a,stencil(k,p))
                enddo
                history(a,p,slot)=value
            enddo
        enddo
    end subroutine spatial_interpolation

    subroutine interpolate_packets(ni,nj,n,nf,map,history,w0,w1,w2,dt,sn,sq_,qk_,qn,f,g,rho,u,v,T,Fx,Fy)
        integer,intent(in) :: ni,nj,n,nf,map(3,n)
        real(kind=8),intent(in) :: history(0:19,nf,0:2),w0,w1,w2,dt,sn,sq_,qk_,qn
        real(kind=8),intent(inout) :: f(ni,nj,0:8),g(ni,nj,0:4),rho(ni,nj),u(ni,nj),v(ni,nj),T(ni,nj),Fx(ni,nj),Fy(ni,nj)
        integer :: p,i,j,k,a
        real(kind=8) :: fv(0:8),gv(0:4),buffer(0:19)
        if(n==0) return
        ! 时间权重必须逐个标量捕获，保留原P100修复，避免firstprivate数组形参。
        !$acc parallel loop gang vector present(map,history,f,g,rho,u,v,T,Fx,Fy) async(1) &
        !$acc& firstprivate(w0,w1,w2,dt,sn,sq_,qk_,qn) private(i,j,k,a,fv,gv,buffer)
        do p=1,n
            i=map(1,p)
            j=map(2,p)
            k=map(3,p)
            do a=0,19
                buffer(a)=w0*history(a,k,0)+w1*history(a,k,1)+w2*history(a,k,2)
            enddo
            call unpack_interface_values(dt,sn,sq_,qk_,qn,buffer,fv,gv,rho(i,j),u(i,j),v(i,j),T(i,j),Fx(i,j),Fy(i,j))
            do a=0,8
                f(i,j,a)=fv(a)
            enddo
            do a=0,4
                g(i,j,a)=gv(a)
            enddo
        enddo
    end subroutine interpolate_packets

    subroutine coarse_to_fine(l,w)
        integer,intent(in) :: l
        real(kind=8),intent(in) :: w(0:2)
        integer :: b
        do b=grid(l)%firstBlock,grid(l)%lastBlock
            call interpolate_packets(block(b)%ni,block(b)%nj,size(coupling(l)%fineTarget(b)%map,2),coupling(l)%nFine, &
                coupling(l)%fineTarget(b)%map,coupling(l)%history,w(0),w(1),w(2),grid(l)%dx,SnuLevel(l),SqLevel(l), &
                QkLevel(l),QnuLevel(l),block(b)%f,block(b)%g,block(b)%rho,block(b)%u,block(b)%v,block(b)%T, &
                block(b)%Fx,block(b)%Fy)
        enddo
    end subroutine coarse_to_fine

    subroutine fine_to_coarse(l)
        integer,intent(in) :: l
        integer :: b,lev
        do b=grid(l)%firstBlock,grid(l)%lastBlock
            call pack_nodes(block(b)%ni,block(b)%nj,size(coupling(l)%fineSource(b)%map,2),coupling(l)%nCoarse, &
                coupling(l)%fineSource(b)%map,coupling(l)%restriction,grid(l)%dx,SnuLevel(l),SqLevel(l), &
                QkLevel(l),QnuLevel(l),block(b)%f,block(b)%g,block(b)%rho,block(b)%u,block(b)%v,block(b)%T, &
                block(b)%Fx,block(b)%Fy)
        enddo
        lev=l+1
        do b=grid(lev)%firstBlock,grid(lev)%lastBlock
            call unpack_nodes(block(b)%ni,block(b)%nj,size(coupling(l)%coarseTarget(b)%map,2),coupling(l)%nCoarse, &
                coupling(l)%coarseTarget(b)%map,coupling(l)%restriction,grid(lev)%dx,SnuLevel(lev),SqLevel(lev), &
                QkLevel(lev),QnuLevel(lev),block(b)%f,block(b)%g,block(b)%rho,block(b)%u,block(b)%v,block(b)%T, &
                block(b)%Fx,block(b)%Fy)
        enddo
    end subroutine fine_to_coarse

    subroutine unpack_nodes(ni,nj,n,np,map,packet,dt,sn,sq_,qk_,qn,f,g,rho,u,v,T,Fx,Fy)
        integer,intent(in) :: ni,nj,n,np,map(3,n)
        real(kind=8),intent(in) :: packet(0:19,np),dt,sn,sq_,qk_,qn
        real(kind=8),intent(inout) :: f(ni,nj,0:8),g(ni,nj,0:4),rho(ni,nj),u(ni,nj),v(ni,nj),T(ni,nj),Fx(ni,nj),Fy(ni,nj)
        integer :: p,i,j,k,a
        real(kind=8) :: fv(0:8),gv(0:4),buffer(0:19)
        if(n==0) return
        !$acc parallel loop gang vector present(map,packet,f,g,rho,u,v,T,Fx,Fy) async(1) &
        !$acc& firstprivate(dt,sn,sq_,qk_,qn) private(i,j,k,a,fv,gv,buffer)
        do p=1,n
            i=map(1,p)
            j=map(2,p)
            k=map(3,p)
            do a=0,19
                buffer(a)=packet(a,k)
            enddo
            call unpack_interface_values(dt,sn,sq_,qk_,qn,buffer,fv,gv,rho(i,j),u(i,j),v(i,j),T(i,j),Fx(i,j),Fy(i,j))
            do a=0,8
                f(i,j,a)=fv(a)
            enddo
            do a=0,4
                g(i,j,a)=gv(a)
            enddo
        enddo
    end subroutine unpack_nodes

    subroutine history_copy(n,history,initial)
        integer,intent(in) :: n
        logical,intent(in) :: initial
        real(kind=8),intent(inout) :: history(0:19,n,0:2)
        integer :: p,a
        !$acc parallel loop collapse(2) present(history) firstprivate(initial) async(1)
        do p=1,n
            do a=0,19
                if(initial) then
                    history(a,p,0)=history(a,p,1)
                    history(a,p,2)=history(a,p,1)
                else
                    history(a,p,0)=history(a,p,1)
                    history(a,p,1)=history(a,p,2)
                endif
            enddo
        enddo
    end subroutine history_copy

    subroutine initial_interfaces()
        integer :: l
        do l=1,nGridLevels-1
            call fine_to_coarse(l)
        enddo
        do l=1,nGridLevels-1
            call save_coarse_history(l,1)
            call history_copy(coupling(l)%nFine,coupling(l)%history,.true.)
            call coarse_to_fine(l,[0.0d0,1.0d0,0.0d0])
        enddo
        !$acc wait(1)
    end subroutine initial_interfaces

!===================================================================================================
! 递归同步：先预测本级终点，再让外侧更细级走两个子步，最后双向补齐接收层。
! 每个接口独立记录历史是否就绪；不能用只在最粗同步点变化的全局itc判断中间级启动。
!===================================================================================================
    recursive subroutine advance_multilevel(l)
        integer,intent(in) :: l
        integer :: k
        real(kind=8) :: theta,w(0:2)
        call advance_level(l)
        if(l==1) return
        call save_coarse_history(l-1,2)
        do k=1,levelRefineRatio
            theta=dble(k-1)/dble(levelRefineRatio)
            if(coupling(l-1)%cycles==0) then
                w=[0.0d0,1.0d0-theta,theta]
            else
                w=[0.5d0*theta*(theta-1.0d0),1.0d0-theta**2,0.5d0*theta*(theta+1.0d0)]
            endif
            call coarse_to_fine(l-1,w)
            call advance_multilevel(l-1)
        enddo
        call fine_to_coarse(l-1)
        call coarse_to_fine(l-1,[0.0d0,0.0d0,1.0d0])
        call history_copy(coupling(l-1)%nFine,coupling(l-1)%history,.false.)
        coupling(l-1)%cycles=coupling(l-1)%cycles+1
    end subroutine advance_multilevel

    subroutine advance_level(l)
        integer,intent(in) :: l
        integer :: b
        do b=grid(l)%firstBlock,grid(l)%lastBlock
            call collision(block(b)%ni,block(b)%nj,grid(l)%dx,block(b)%f,block(b)%f_post,block(b)%rho, &
                block(b)%u,block(b)%v,block(b)%Fx,block(b)%Fy,block(b)%T,SnuLevel(l),SqLevel(l),gBeta &
#ifdef SideHeatedHa
                ,B2sigemarho &
#endif
                )
        enddo
        call exchange_level_halos(l,.false.)
        do b=grid(l)%firstBlock,grid(l)%lastBlock
            call streaming(block(b)%ni,block(b)%nj,block(b)%f,block(b)%f_post)
            call bounceback(block(b)%ni,block(b)%nj,block(b)%f,block(b)%f_post, &
                block(b)%wall(1),block(b)%wall(2),block(b)%wall(3),block(b)%wall(4))
            call macro(block(b)%ni,block(b)%nj,grid(l)%dx,block(b)%f,block(b)%rho,block(b)%u,block(b)%v, &
                block(b)%Fx,block(b)%Fy)
        enddo
        do b=grid(l)%firstBlock,grid(l)%lastBlock
            call collisionT(block(b)%ni,block(b)%nj,block(b)%g,block(b)%g_post,block(b)%u,block(b)%v, &
                block(b)%T,QkLevel(l),QnuLevel(l))
        enddo
        call exchange_level_halos(l,.true.)
        do b=grid(l)%firstBlock,grid(l)%lastBlock
            call streamingT(block(b)%ni,block(b)%nj,block(b)%g,block(b)%g_post)
            call bouncebackT(block(b)%ni,block(b)%nj,block(b)%g,block(b)%g_post, &
                block(b)%wall(1),block(b)%wall(2),block(b)%wall(3),block(b)%wall(4))
            call macroT(block(b)%ni,block(b)%nj,block(b)%g,block(b)%T)
        enddo
    end subroutine advance_level

    subroutine exchange_level_halos(l,thermal)
        integer,intent(in) :: l
        logical,intent(in) :: thermal
        integer :: k,s,b
        logical :: periodicEnabled
        periodicEnabled=.false.
#ifdef VerticalWallsPeriodicalU
        if(.not.thermal) periodicEnabled=.true.
#endif
#ifdef VerticalWallsPeriodicalT
        if(thermal) periodicEnabled=.true.
#endif
        do k=1,size(halo)
            s=halo(k)%source
            b=halo(k)%target
            if(block(b)%level/=l) cycle
            if(halo(k)%periodic==1.and..not.periodicEnabled) cycle
            if(thermal) then
                call copy_halo(4,block(b)%ni,block(b)%nj,block(s)%ni,block(s)%nj,size(halo(k)%map,2), &
                    halo(k)%map,block(b)%g_post,block(s)%g_post)
            else
                call copy_halo(8,block(b)%ni,block(b)%nj,block(s)%ni,block(s)%nj,size(halo(k)%map,2), &
                    halo(k)%map,block(b)%f_post,block(s)%f_post)
            endif
        enddo
    end subroutine exchange_level_halos

    subroutine copy_halo(nq,ni,nj,mi,mj,n,map,target,source)
        integer,intent(in) :: nq,ni,nj,mi,mj,n,map(4,n)
        real(kind=8),intent(inout) :: target(0:ni+1,0:nj+1,0:nq)
        real(kind=8),intent(in) :: source(0:mi+1,0:mj+1,0:nq)
        integer :: p,a
        !$acc parallel loop collapse(2) present(map,target,source) async(1)
        do p=1,n
            do a=0,nq
                target(map(1,p),map(2,p),a)=source(map(3,p),map(4,p),a)
            enddo
        enddo
    end subroutine copy_halo

!===================================================================================================
! 统计仍使用原来的体积Nu、空间RMS Re与壁面温度梯度。每级按outerBox减innerBox积分。
!===================================================================================================
    subroutine calculate_NuRe(nu,re,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax)
        real(kind=8),intent(out) :: nu,re,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax
        real(kind=8) :: conv,vel2,area,d,weight
        integer :: b,i,j
        call update_host_all(.false.)
        conv=0.0d0
        vel2=0.0d0
        mass=0.0d0
        meanT=0.0d0
        hot=0.0d0
        cold=0.0d0
        tmin=huge(1.0d0)
        tmax=-tmin
        rmin=tmin
        rmax=-tmin
        do b=1,nBlocks
            d=block(b)%dx
            do j=1,block(b)%nj
                do i=1,block(b)%ni
                    area=block(b)%area(i,j)
                    if(area<=0.0d0) cycle
                    if(.not.ieee_is_finite(block(b)%rho(i,j)).or.block(b)%rho(i,j)<=0.0d0.or. &
                        .not.ieee_is_finite(block(b)%T(i,j)).or..not.ieee_is_finite(block(b)%u(i,j)).or. &
                        .not.ieee_is_finite(block(b)%v(i,j))) then
                        write(*,*) 'Invalid owned cell: step,block,level,i,j=',reloadItcOffset+itc,b,block(b)%level,i,j
                        error stop 'Nonfinite field or nonpositive density'
                    endif
#ifdef SideHeatedCell
                    conv=conv+area*block(b)%u(i,j)*block(b)%T(i,j)
#else
                    conv=conv+area*block(b)%v(i,j)*block(b)%T(i,j)
#endif
                    vel2=vel2+area*(block(b)%u(i,j)**2+block(b)%v(i,j)**2)
                    mass=mass+area*block(b)%rho(i,j)
                    meanT=meanT+area*block(b)%T(i,j)
                    tmin=min(tmin,block(b)%T(i,j))
                    tmax=max(tmax,block(b)%T(i,j))
                    rmin=min(rmin,block(b)%rho(i,j))
                    rmax=max(rmax,block(b)%rho(i,j))
                enddo
            enddo
#ifdef SideHeatedCell
            if(block(b)%wall(1).or.block(b)%wall(2)) then
                do j=1,block(b)%nj
                    weight=d/dble(ny)
                    if(block(b)%wall(1)) hot=hot+weight*(8.0d0*Thot-9.0d0*block(b)%T(1,j)+block(b)%T(2,j))/(3.0d0*d)
                    i=block(b)%ni
                    if(block(b)%wall(2)) cold=cold+weight*(-8.0d0*Tcold+9.0d0*block(b)%T(i,j)-block(b)%T(i-1,j))/(3.0d0*d)
                enddo
            endif
#else
            if(block(b)%wall(3).or.block(b)%wall(4)) then
                do i=1,block(b)%ni
                    weight=d/dble(nx)
                    if(block(b)%wall(3)) hot=hot+weight*(8.0d0*Thot-9.0d0*block(b)%T(i,1)+block(b)%T(i,2))/(3.0d0*d)
                    j=block(b)%nj
                    if(block(b)%wall(4)) cold=cold+weight*(-8.0d0*Tcold+9.0d0*block(b)%T(i,j)-block(b)%T(i,j-1))/(3.0d0*d)
                enddo
            endif
#endif
        enddo
        nu=1.0d0+conv/(dble(nx)*dble(ny))*lengthUnit/((Thot-Tcold)*diffusivity)
        re=sqrt(vel2/(dble(nx)*dble(ny)))*lengthUnit/viscosity
        meanT=meanT/(dble(nx)*dble(ny))
        hot=hot*lengthUnit/(Thot-Tcold)
        cold=cold*lengthUnit/(Thot-Tcold)
        call calculate_middle_flux(middle)
    end subroutine calculate_NuRe

    real(kind=8) function section_owned_length(l,q,position,axis) result(weight)
        integer,intent(in) :: l,axis
        real(kind=8),intent(in) :: q,position
        real(kind=8) :: lo,hi,d
        d=0.5d0*grid(l)%dx
        weight=0.0d0
        if(axis==1) then
            if(position<grid(l)%outerBox(1).or.position>=grid(l)%outerBox(2)) return
            lo=max(q-d,grid(l)%outerBox(3))
            hi=min(q+d,grid(l)%outerBox(4))
            weight=max(0.0d0,hi-lo)
            if(l<nGridLevels) then
                if(position>=grid(l)%innerBox(1).and.position<grid(l)%innerBox(2)) &
                    weight=weight-max(0.0d0,min(hi,grid(l)%innerBox(4))-max(lo,grid(l)%innerBox(3)))
            endif
        else
            if(position<grid(l)%outerBox(3).or.position>=grid(l)%outerBox(4)) return
            lo=max(q-d,grid(l)%outerBox(1))
            hi=min(q+d,grid(l)%outerBox(2))
            weight=max(0.0d0,hi-lo)
            if(l<nGridLevels) then
                if(position>=grid(l)%innerBox(3).and.position<grid(l)%innerBox(4)) &
                    weight=weight-max(0.0d0,min(hi,grid(l)%innerBox(2))-max(lo,grid(l)%innerBox(1)))
            endif
        endif
        weight=max(0.0d0,weight)
    end function section_owned_length

    subroutine sample_section(l,fixed,position,axis,tm,um,gradient)
        integer,intent(in) :: l,fixed,axis
        real(kind=8),intent(in) :: position
        real(kind=8),intent(out) :: tm,um,gradient
        integer :: first,k,b,i,j,skin
        real(kind=8) :: q,w(4),dw(4)
        skin=0
        if(l>1) skin=interfaceSkin
        if(axis==1) then
            q=(position-grid(l)%xOffset)/grid(l)%dx+0.5d0
            call section_weights(q-skin,grid(l)%ni-2*skin,first,w,dw)
        else
            q=(position-grid(l)%yOffset)/grid(l)%dx+0.5d0
            call section_weights(q-skin,grid(l)%nj-2*skin,first,w,dw)
        endif
        first=first+skin
        tm=0.0d0
        um=0.0d0
        gradient=0.0d0
        do k=1,4
            if(axis==1) then
                call locate_level_node(l,first+k-1,fixed,b,i,j,.false.)
            else
                call locate_level_node(l,fixed,first+k-1,b,i,j,.false.)
            endif
            if(b==0) error stop 'Middle-section stencil crosses an unstored node'
            tm=tm+w(k)*block(b)%T(i,j)
            if(axis==1) then
                um=um+w(k)*block(b)%u(i,j)
            else
                um=um+w(k)*block(b)%v(i,j)
            endif
            gradient=gradient+dw(k)*block(b)%T(i,j)/grid(l)%dx
        enddo
    end subroutine sample_section

    subroutine calculate_middle_flux(middle)
        real(kind=8),intent(out) :: middle
        integer :: b,l,k,fixed
        real(kind=8) :: position,q,weight,tm,um,gradient,totalLength
        middle=0.0d0
        totalLength=0.0d0
        do b=1,nBlocks
            l=block(b)%level
#ifdef SideHeatedCell
            position=0.5d0*nx
            if(position<block(b)%xOffset.or.position>=block(b)%xOffset+block(b)%ni*block(b)%dx) cycle
            do k=1,block(b)%nj
                q=block(b)%yOffset+(k-0.5d0)*block(b)%dx
                weight=section_owned_length(l,q,position,1)
                if(weight<=0.0d0) cycle
                fixed=block(b)%globalJ+k-1
                call sample_section(l,fixed,position,1,tm,um,gradient)
                middle=middle+(um*tm/diffusivity-gradient)*weight/dble(ny)
#else
            position=0.5d0*ny
            if(position<block(b)%yOffset.or.position>=block(b)%yOffset+block(b)%nj*block(b)%dx) cycle
            do k=1,block(b)%ni
                q=block(b)%xOffset+(k-0.5d0)*block(b)%dx
                weight=section_owned_length(l,q,position,2)
                if(weight<=0.0d0) cycle
                fixed=block(b)%globalI+k-1
                call sample_section(l,fixed,position,2,tm,um,gradient)
                middle=middle+(um*tm/diffusivity-gradient)*weight/dble(nx)
#endif
                totalLength=totalLength+weight
            enddo
        enddo
#ifdef SideHeatedCell
        if(abs(totalLength-ny)>1.0d-10*ny) error stop 'Middle section does not cover the domain once'
#else
        if(abs(totalLength-nx)>1.0d-10*nx) error stop 'Middle section does not cover the domain once'
#endif
        middle=middle*lengthUnit/(Thot-Tcold)
    end subroutine calculate_middle_flux

    subroutine calNuRe()
        integer :: k
        real(kind=8) :: nu,re,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax
        call calculate_NuRe(nu,re,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax)
        open(newunit=k,file=NuReHistoryFile,status='old',position='append')
        write(k,'(2(I0,1X),12(ES24.16E3,1X))') NuReSampleCount+1,reloadItcOffset+itc, &
            dble(reloadItcOffset+itc)/timeUnit,nu,re,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax
        close(k)
        NuReSampleCount=NuReSampleCount+1
        write(*,'(a,f12.5,2(a,es13.5))') 't_ff=',dble(reloadItcOffset+itc)/timeUnit,' NuVolAvg=',nu,' ReVolRMS=',re
    end subroutine calNuRe

    subroutine output_final_NuRe()
        integer :: k
        real(kind=8) :: nu,re,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax
        call calculate_NuRe(nu,re,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax)
        open(newunit=k,file='FinalNuRe_2DOpenaccMultiblock.dat',status='replace')
        write(k,'(a,i0)') '# itc_total = ',reloadItcOffset+itc
        ! 保留稳态最终诊断格式，区分终场时刻与最近一次误差检查时刻。
#ifdef steadyFlow
        write(k,'(a,l1)') '# converged = ',errorU<=epsU.and.errorT<=epsT
        write(k,'(a,i0)') '# error_check_itc = ', &
            ((reloadItcOffset+itc)/steadyCheckIntervalItc)*steadyCheckIntervalItc
        write(k,'(a,2(es24.16e3,1x))') '# errorU errorT = ',errorU,errorT
#else
        write(k,'(a)') '# flow_mode = unsteady; instantaneous final field, not a time average'
#endif
        write(k,'(a)') '# t_ff NuVolAvg ReVolRMS Nu_hot Nu_cold Nu_middle mass meanT Tmin Tmax rhoMin rhoMax'
        write(k,'(12(ES24.16E3,1X))') dble(reloadItcOffset+itc)/timeUnit,nu,re,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax
        close(k)
        open(newunit=k,file=settingsFile,status='old',position='append')
        write(k,*) 'Final diagnostics at cumulative fine step:',reloadItcOffset+itc
#ifdef steadyFlow
        write(k,*) 'Steady converged:',errorU<=epsU.and.errorT<=epsT
        write(k,*) 'errorU, errorT:',errorU,errorT
#endif
        write(k,*) 'Nu_global =',nu
        write(k,*) 'Nu_hot =',hot
        write(k,*) 'Nu_cold =',cold
        write(k,*) 'Nu_middle =',middle
        write(k,*) 'Re_rms =',re
        close(k)
        call record_output_file('FinalNuRe_2DOpenaccMultiblock.dat')
    end subroutine output_final_NuRe

    subroutine check()
#ifdef steadyFlow
        integer :: b,k
        real(kind=8) :: du,uu,dt_,tt
        call update_host_all(.false.)
        du=0.0d0
        uu=0.0d0
        dt_=0.0d0
        tt=0.0d0
        do b=1,nBlocks
            du=du+sum(block(b)%area*((block(b)%u-block(b)%up)**2+(block(b)%v-block(b)%vp)**2))
            uu=uu+sum(block(b)%area*(block(b)%u**2+block(b)%v**2))
            dt_=dt_+sum(block(b)%area*abs(block(b)%T-block(b)%Tp))
            tt=tt+sum(block(b)%area*abs(block(b)%T))
            block(b)%up=block(b)%u
            block(b)%vp=block(b)%v
            block(b)%Tp=block(b)%T
        enddo
        errorU=sqrt(du/max(uu,1.0d-300))
        errorT=dt_
        if(tt>1.0d-30) errorT=dt_/tt
        open(newunit=k,file=convergenceFile,status='unknown',position='append')
        write(k,'(i0,1x,3(es24.16e3,1x))') reloadItcOffset+itc,dble(reloadItcOffset+itc)/timeUnit,errorU,errorT
        close(k)
#endif
    end subroutine check

!===================================================================================================
! 后处理相关文本、索引和前缀与原Multiblock版本一致，保留原有画图和文本统计方式。
! 多级快照由配套MultilevelPostprocess读取；重启仍使用Multilevel独立前缀。
! 新快照包含每块的级别、格距、偏移、积分区域和完整u/v/T/rho数组。
!===================================================================================================
    function output_name(prefix,suffix) result(name)
        character(len=*),intent(in) :: prefix,suffix
        character(len=256) :: name
        character(len=48) :: timeText,stepText
        write(timeText,'(f20.8)') dble(reloadItcOffset+itc)/timeUnit
        write(stepText,'(i0)') reloadItcOffset+itc
        name=prefix//'-'//trim(adjustl(timeText))//'tff-'//trim(stepText)//'itc'//suffix
    end function output_name

    subroutine record_output_file(name)
        character(len=*),intent(in) :: name
        integer :: k
        open(newunit=k,file='OutputFiles_2DOpenaccMultiblock.dat',status='unknown',position='append')
        write(k,'(i0,1x,a)') reloadItcOffset+itc,trim(name)
        close(k)
    end subroutine record_output_file

    subroutine output_Tecplot()
        integer :: k,b,i,j,l
        character(len=256) :: name
        name=output_name(pltFolderPrefix,'.dat')
        call update_host_all(.false.)
        open(newunit=k,file=trim(name),status='replace')
        write(k,'(a)') 'VARIABLES="x/L","y/L","u","v","T","rho"'
        do b=1,nBlocks
            l=block(b)%level
            write(k,'(a,i0,a,i0,a,i0,a,i0,a,es24.16e3)') 'ZONE T="level',l,'_block',b,'", I=',block(b)%ni, &
                ', J=',block(b)%nj,', F=POINT, SOLUTIONTIME=',dble(reloadItcOffset+itc)/timeUnit
            do j=1,block(b)%nj
                do i=1,block(b)%ni
                    write(k,'(6(es24.16e3,1x))') &
                        (block(b)%xOffset+(i-0.5d0)*block(b)%dx)/lengthUnit, &
                        (block(b)%yOffset+(j-0.5d0)*block(b)%dx)/lengthUnit,block(b)%u(i,j),block(b)%v(i,j), &
                        block(b)%T(i,j),block(b)%rho(i,j)
                enddo
            enddo
        enddo
        close(k)
        call record_output_file(name)
    end subroutine output_Tecplot

    subroutine output_SnapshotFile()
        integer :: k,b,l
        character(len=16) :: magic='ML2DSNAPSHOT0001'
        character(len=256) :: name
        name=output_name(snapshotFilePrefix,'.bin')
        call update_host_all(.false.)
        open(newunit=k,file=trim(name),access='stream',form='unformatted',status='replace')
        write(k) magic,current_config(),current_physics(),reloadItcOffset+itc,nBlocks,timeUnit
        do b=1,nBlocks
            l=block(b)%level
            write(k) l,block(b)%ni,block(b)%nj,block(b)%dx,block(b)%xOffset,block(b)%yOffset
            write(k) grid(l)%outerBox,grid(l)%innerBox
            write(k) block(b)%u,block(b)%v,block(b)%T,block(b)%rho
        enddo
        close(k)
        call record_output_file(name)
    end subroutine output_SnapshotFile

    function current_config() result(config)
        integer :: config(16)
        config=[nx,ny,nGridLevels,interfaceSkin,finestLayerCellsLeftRight,finestLayerCellsBottomTop, &
            innerLayerCellsLeftRight,innerLayerCellsBottomTop,coarsestStep,0,0,0,0,0,0,0]
#ifdef SideHeatedCell
        config(10)=1
#else
        config(10)=2
#endif
#ifdef VerticalWallsPeriodicalU
        config(11)=1
#endif
#ifdef VerticalWallsPeriodicalT
        config(12)=1
#endif
#ifdef HorizontalWallsConstT
        config(13)=1
#endif
#ifdef VerticalWallsConstT
        config(14)=1
#endif
#ifdef SideHeatedHa
        config(15)=1
#endif
#ifdef steadyFlow
        config(16)=1
#else
        config(16)=2
#endif
    end function current_config

    function current_physics() result(physics)
        real(kind=8) :: physics(12)
        physics=[Rayleigh,Prandtl,Mach,Thot,Tcold,paraA,Snu,Sq,Qk,Qnu,0.0d0,0.0d0]
#ifdef SideHeatedHa
        physics(11)=B2sigemarho
        physics(12)=phi
#endif
    end function current_physics

    subroutine write_settings()
        integer :: k,l,b
        logical :: exists
        real(kind=8) :: owned,updates,stored
        inquire(file=settingsFile,exist=exists)
        if(loadInitField==1.and.exists) then
            open(newunit=k,file=settingsFile,status='old',position='append')
        else
            open(newunit=k,file=settingsFile,status='replace')
        endif
        write(k,'(a)') 'Original D2Q9 + Luo D2Q5; no regularization; raw coarse/fine moment scaling'
        write(k,'(a)') 'Parent SHA256: 84ca3661eaa2544cca9bb19403b1eaa7d97e82aa4f09ccbbd7cf657de915aa2c'
        write(k,'(a,16(i0,1x))') 'Configuration: ',current_config()
        write(k,'(a,12(es24.16e3,1x))') 'Physics: ',current_physics()
        write(k,'(a,i0)') 'nBlocks=',nBlocks
        write(k,'(a,3(es24.16e3,1x))') 'Fine viscosity, diffusivity, timeUnit: ',viscosity,diffusivity,timeUnit
        write(k,'(a,i0,a,es24.16e3)') 'Target fine steps: ',itc_max,'; t_ff: ',dble(itc_max)/timeUnit
        write(k,'(a,4(i0,1x))') 'Sampling,restart,Tecplot,convergence intervals (fine steps): ', &
            outputSnapshotIntervalItc,reloadFileIntervalItc,outputPltFileIntervalItc,steadyCheckIntervalItc
        stored=0.0d0
        updates=0.0d0
        do l=1,nGridLevels
            owned=0.0d0
            do b=grid(l)%firstBlock,grid(l)%lastBlock
                owned=owned+sum(block(b)%area)
                stored=stored+dble(block(b)%ni)*block(b)%nj
                updates=updates+dble(block(b)%ni)*block(b)%nj/grid(l)%dx
                write(k,'(a,6(i0,1x),3(es18.10,1x))') 'Block: ',b,l,block(b)%ni,block(b)%nj, &
                    block(b)%globalI,block(b)%globalJ,block(b)%dx,block(b)%xOffset,block(b)%yOffset
            enddo
            write(k,'(a,i0,a,7(es24.16e3,1x))') 'Level ',l,' dx,Snu,Sq,Qk,Qnu,area,tau: ', &
                grid(l)%dx,SnuLevel(l),SqLevel(l),QkLevel(l),QnuLevel(l),owned,1.0d0/SnuLevel(l)
            write(k,'(a,4(es24.16e3,1x))') 'outerBox: ',grid(l)%outerBox
            write(k,'(a,4(es24.16e3,1x))') 'innerBox: ',grid(l)%innerBox
        enddo
        write(k,'(a,es24.16)') 'Stored nodes / full fine nodes = ',stored/(dble(nx)*ny)
        write(k,'(a,es24.16)') 'Node updates / full fine updates (excluding transfers) = ',updates/(dble(nx)*ny)
        write(k,'(a,i0)') 'Restart step = ',reloadItcOffset
        close(k)
    end subroutine write_settings

!===================================================================================================
! 重启记录同步点完整分布函数、宏观量和每个接口的三个时间层；配置不匹配时拒绝读取。
!===================================================================================================
    subroutine output_ReloadFile()
        integer :: k,b,l,total
        character(len=16) :: magic='ML2DRESTART0001'
        character(len=256) :: name
        total=reloadItcOffset+itc
        name=output_name(reloadFilePrefix,'.bin')
        call update_host_all(.true.)
        do l=1,nGridLevels-1
            call interface_data(l,2)
        enddo
        !$acc wait(1)
        open(newunit=k,file=trim(name),access='stream',form='unformatted',status='replace')
        write(k) magic,current_config(),current_physics(),total,NuReSampleCount,errorU,errorT,nBlocks
        do b=1,nBlocks
            write(k) block(b)%level,block(b)%ni,block(b)%nj,block(b)%globalI,block(b)%globalJ
            write(k) block(b)%dx,block(b)%xOffset,block(b)%yOffset
            write(k) block(b)%f,block(b)%g,block(b)%rho,block(b)%u,block(b)%v,block(b)%T,block(b)%Fx,block(b)%Fy
#ifdef steadyFlow
            write(k) block(b)%up,block(b)%vp,block(b)%Tp
#endif
        enddo
        do l=1,nGridLevels-1
            write(k) coupling(l)%nFine,coupling(l)%cycles,coupling(l)%history
        enddo
        close(k)
        open(newunit=k,file=reloadFilePrefix//'-latest.meta',status='replace')
        write(k,'(a)') 'ML2DRESTART0001'
        write(k,'(a)') trim(name)
        write(k,'(a)') '# total_fine_steps NuReSampleCount t_ff'
        write(k,'(2(i0,1x),es24.16e3)') total,NuReSampleCount,dble(total)/timeUnit
        close(k)
        call record_output_file(name)
    end subroutine output_ReloadFile

    subroutine read_restart()
        integer :: k,b,l,ios,config(16),total,samples,nb,geometry(5),nf,cycles
        real(kind=8) :: physics(12),position(3),savedU,savedT
        character(len=16) :: magic
        character(len=256) :: name
        character(len=512) :: line
        logical :: exists
        inquire(file=reloadFilePrefix//'-latest.meta',exist=exists)
        name=reloadFileName
        if(exists) then
            open(newunit=k,file=reloadFilePrefix//'-latest.meta',status='old')
            read(k,'(a)',iostat=ios) line
            if(ios/=0.or.trim(line)/='ML2DRESTART0001') error stop 'Invalid multilevel latest metadata'
            read(k,'(a)',iostat=ios) name
            close(k)
            if(ios/=0) error stop 'Missing checkpoint name in latest metadata'
        endif
        if(len_trim(name)==0) error stop 'Set reloadFileName or provide multilevel latest metadata'
        open(newunit=k,file=trim(name),access='stream',form='unformatted',status='old',iostat=ios)
        if(ios/=0) error stop 'Cannot open multilevel checkpoint'
        read(k,iostat=ios) magic,config,physics,total,samples,savedU,savedT,nb
        if(ios/=0) error stop 'Truncated multilevel checkpoint header'
        if(magic/='ML2DRESTART0001') error stop 'Checkpoint is not the multilevel format'
        if(any(config/=current_config()).or.nb/=nBlocks) error stop 'Checkpoint mesh/BC configuration mismatch'
        if(any(abs(physics-current_physics())>1.0d-12*max(1.0d0,abs(current_physics())))) &
            error stop 'Checkpoint physical parameters mismatch'
        if(total<0.or.mod(total,coarsestStep)/=0.or.samples<0) error stop 'Checkpoint is not a valid synchronized step'
        if(.not.ieee_is_finite(savedU).or..not.ieee_is_finite(savedT)) error stop 'Nonfinite saved convergence errors'
        do b=1,nBlocks
            read(k,iostat=ios) geometry,position
            if(ios/=0) error stop 'Truncated checkpoint geometry'
            if(any(geometry/=[block(b)%level,block(b)%ni,block(b)%nj,block(b)%globalI,block(b)%globalJ])) &
                error stop 'Checkpoint block geometry mismatch'
            if(any(abs(position-[block(b)%dx,block(b)%xOffset,block(b)%yOffset])>1.0d-12)) &
                error stop 'Checkpoint block coordinates mismatch'
            read(k,iostat=ios) block(b)%f,block(b)%g,block(b)%rho,block(b)%u,block(b)%v,block(b)%T,block(b)%Fx,block(b)%Fy
            if(ios/=0) error stop 'Truncated checkpoint fields'
            if(any(.not.ieee_is_finite(block(b)%f)).or.any(.not.ieee_is_finite(block(b)%g)).or. &
                any(.not.ieee_is_finite(block(b)%rho)).or.any(.not.ieee_is_finite(block(b)%u)).or. &
                any(.not.ieee_is_finite(block(b)%v)).or.any(.not.ieee_is_finite(block(b)%T)).or. &
                any(.not.ieee_is_finite(block(b)%Fx)).or.any(.not.ieee_is_finite(block(b)%Fy))) &
                error stop 'Nonfinite checkpoint fields'
            if(any(block(b)%rho<=0.0d0.and.block(b)%area>0.0d0)) &
                error stop 'Nonpositive density in an owned checkpoint cell'
#ifdef steadyFlow
            read(k,iostat=ios) block(b)%up,block(b)%vp,block(b)%Tp
            if(ios/=0) error stop 'Truncated checkpoint convergence fields'
            if(any(.not.ieee_is_finite(block(b)%up)).or.any(.not.ieee_is_finite(block(b)%vp)).or. &
                any(.not.ieee_is_finite(block(b)%Tp))) error stop 'Nonfinite checkpoint convergence fields'
#endif
        enddo
        do l=1,nGridLevels-1
            read(k,iostat=ios) nf,cycles
            if(ios/=0) error stop 'Truncated checkpoint interface header'
            if(nf/=coupling(l)%nFine.or.cycles/=total/levelRefineRatio**l) &
                error stop 'Checkpoint interface history mismatch'
            read(k,iostat=ios) coupling(l)%history
            if(ios/=0) error stop 'Truncated checkpoint interface history'
            if(any(.not.ieee_is_finite(coupling(l)%history))) error stop 'Nonfinite checkpoint interface history'
            coupling(l)%cycles=cycles
        enddo
        read(k,iostat=ios) magic
        if(ios>=0) error stop 'Unexpected trailing checkpoint data'
        close(k)
        reloadItcOffset=total
        NuReSampleCount=samples
        errorU=savedU
        errorT=savedT
        call truncate_table(NuReHistoryFile,2,samples)
#ifdef steadyFlow
        call truncate_table(convergenceFile,1,-1)
#endif
        call truncate_table('OutputFiles_2DOpenaccMultiblock.dat',1,-1)
    end subroutine read_restart

    subroutine truncate_table(name,stepColumn,expectedRows)
        character(len=*),intent(in) :: name
        integer,intent(in) :: stepColumn,expectedRows
        integer :: input,temp,output,ios,ios2,a,step,rows,previousStep
        real(kind=8) :: values(12)
        logical :: exists
        character(len=2048) :: line
        inquire(file=name,exist=exists)
        if(.not.exists) then
            if(expectedRows>0) error stop 'Checkpoint history file is missing'
            return
        endif
        open(newunit=input,file=name,status='old')
        open(newunit=temp,file=name//'.restart.tmp',status='replace',action='readwrite')
        rows=0
        previousStep=-1
        do
            read(input,'(a)',iostat=ios) line
            if(ios/=0) exit
            if(len_trim(line)==0) cycle
            if(line(1:1)=='#') then
                write(temp,'(a)') trim(line)
                cycle
            endif
            if(stepColumn==2) then
                read(line,*,iostat=ios2) a,step
            else
                read(line,*,iostat=ios2) step
            endif
            if(ios2/=0) exit
            if(step>reloadItcOffset) exit
            if(step<0.or.step<previousStep) error stop 'Retained history steps decrease'
            if(stepColumn==2.and.step==previousStep) error stop 'Repeated retained Nu/Re history step'
            if(stepColumn==2) then
                read(line,*,iostat=ios2) a,step,values
                if(ios2/=0) error stop 'Incomplete retained Nu/Re history row'
                if(a/=rows+1.or.any(.not.ieee_is_finite(values))) error stop 'Invalid retained Nu/Re history row'
                if(abs(values(1)-dble(step)/timeUnit)>1.0d-12*max(1.0d0,abs(values(1)))) &
                    error stop 'Retained Nu/Re history time/step mismatch'
            endif
            previousStep=step
            rows=rows+1
            write(temp,'(a)') trim(line)
        enddo
        close(input)
        if(expectedRows>=0.and.rows/=expectedRows) error stop 'Checkpoint sample count and retained history disagree'
        rewind(temp)
        open(newunit=output,file=name,status='replace')
        do
            read(temp,'(a)',iostat=ios) line
            if(ios/=0) exit
            write(output,'(a)') trim(line)
        enddo
        close(output)
        close(temp,status='delete')
    end subroutine truncate_table

    subroutine initial_output_files()
        integer :: k
        open(newunit=k,file=NuReHistoryFile,status='replace')
        write(k,'(a)') '# sample itc_total t_ff NuVolAvg ReVolRMS Nu_hot Nu_cold Nu_middle mass meanT Tmin Tmax rhoMin rhoMax'
        close(k)
        open(newunit=k,file=convergenceFile,status='replace')
        write(k,'(a)') '# itc_total t_ff errorU errorT'
        close(k)
        open(newunit=k,file='OutputFiles_2DOpenaccMultiblock.dat',status='replace')
        close(k)
    end subroutine initial_output_files
    end module multilevel_solver

!=============================================================
! 主程序：设备 -> 初始化 -> 粗细同步推进 -> 统计与输出 -> 释放
!=============================================================
    program main
        use openacc
        use multilevel_solver
        implicit none
        integer :: numAccDevices,b,lastSnapshotItc,lastPltItc,lastReloadItc
        integer(kind=8) :: wallStart,wallEnd,wallRate
        real(kind=8) :: seconds,latticeUpdateCount

        numAccDevices=acc_get_num_devices(acc_device_nvidia)
        if(numAccDevices>0) then
            call acc_set_device_type(acc_device_nvidia)
            call acc_set_device_num(0,acc_device_nvidia)
        else
            call acc_set_device_type(acc_device_host)
        endif
        ! 设置设备后由首个数据区完成初始化，避免部分运行时重复acc_init报错。
        call initial_geometry()
        call initial_fields()
        if(loadInitField==1) then
            call read_restart()
        else
            reloadItcOffset=0
            call initial_output_files()
        endif
        call enter_data_2d_openacc()
        if(loadInitField==0) call initial_interfaces()
        call write_settings()
        call system_clock(wallStart,wallRate)
        lastSnapshotItc=-1
        lastPltItc=-1
        lastReloadItc=-1
        itc=0

        ! 只有所有级别到达同一物理时刻后才更新itc、采样、写快照或保存重启。
        do while(reloadItcOffset+itc<itc_max)
#ifdef steadyFlow
            if(errorU<=epsU.and.errorT<=epsT) exit
#endif
            call advance_multilevel(nGridLevels)
            itc=itc+coarsestStep
            if(mod(reloadItcOffset+itc,steadyCheckIntervalItc)==0) call check()
            if(mod(reloadItcOffset+itc,outputSnapshotIntervalItc)==0) then
                call calNuRe()
                if(outputSnapshotFile==1) then
                    call output_SnapshotFile()
                    lastSnapshotItc=itc
                endif
            endif
            if(outputPltFile==1.and.mod(reloadItcOffset+itc,outputPltFileIntervalItc)==0) then
                call output_Tecplot()
                lastPltItc=itc
            endif
            if(outputReloadFile==1.and.mod(reloadItcOffset+itc,reloadFileIntervalItc)==0) then
                call output_ReloadFile()
                lastReloadItc=itc
            endif
        enddo
        !$acc wait(1)
        call output_final_NuRe()
        if(outputSnapshotFile==1.and.lastSnapshotItc/=itc) call output_SnapshotFile()
        if(outputPltFile==1.and.lastPltItc/=itc) call output_Tecplot()
        if(outputReloadFile==1.and.lastReloadItc/=itc) call output_ReloadFile()
        call system_clock(wallEnd)
        seconds=dble(wallEnd-wallStart)/dble(wallRate)
        latticeUpdateCount=0.0d0
        do b=1,nBlocks
            latticeUpdateCount=latticeUpdateCount+dble(block(b)%ni)*block(b)%nj*dble(itc)/block(b)%dx
        enddo
        write(*,'(a,i0,a,i0)') 'Completed total fine steps: ',reloadItcOffset+itc,'; levels: ',nGridLevels
        write(*,'(a,f12.4,a,f12.4)') 'Wall seconds: ',seconds,'; MLUPS: ',latticeUpdateCount/max(seconds,1.0d-30)/1.0d6
        call exit_data_2d_openacc()
        call acc_shutdown(acc_device_default)
    end program main

!===================================================================================================
  subroutine coarse_source_stencil(qx,qy,ni,nj,skin,sourceX,sourceY,width,wx,wy)
    implicit none
    !$acc routine seq
    real(kind=8), intent(in) :: qx,qy
    integer(kind=4), intent(in) :: ni,nj,skin
    integer(kind=4), intent(out) :: sourceX,sourceY,width
    real(kind=8), intent(out) :: wx(4),wy(4)
    
    ! qx/qy 是粗网格连续下标，不是全局物理坐标；sourceX/sourceY 是来源粗数组的局部整数下标。
    ! nint 是取最近的整数。
    if(abs(qx-nint(qx))<1.0d-12 .AND. abs(qy-nint(qy))<1.0d-12) then  !两个方向都是整数粗下标，才能直接取粗节点数值
        sourceX=nint(qx)
        sourceY=nint(qy)
        width=1
        wx=0.0d0
        wy=0.0d0
        wx(1)=1.0d0
        wy(1)=1.0d0
    else                !至少一个方向的连续下标不是整数，floor 表示向下取整
        sourceX=max(1+skin,min(floor(qx)-1,ni-skin-3))
        sourceY=max(1+skin,min(floor(qy)-1,nj-skin-3))
        width=4
        call lagrange_weights(qx-sourceX,wx)   !得到x方向权重
        call lagrange_weights(qy-sourceY,wy)   !得到y方向权重
    endif
  end subroutine coarse_source_stencil

!===================================================================================================
  subroutine lagrange_weights(q,w)

    implicit none
    !$acc routine seq

    ! q 为目标点相对于模板首节点的距离，以来源网格格距为单位。
    ! 模板的相对坐标是 0、1、2、3；w(1:4) 是三次多项式的四点插值权重。

    real(kind=8), intent(in) :: q
    real(kind=8), intent(out) :: w(4)
    integer(kind=4) :: a, k

    w=1.0d0   !下面通过连续相乘来计算权重，所以初值必须是 1,如果初值为 0，结果都会是 0
    do a=0,3
        do k=0,3
            if(a/=k) w(a+1)=w(a+1)*(q-dble(k))/dble(a-k)
        enddo
    enddo
  end subroutine lagrange_weights

!===================================================================================================
  subroutine pack_interface_values(dt,Snu,Sq,Qk,Qnu,fv,gv,rh,ux,uy,temp,fx,fy,buffer)
    use commondata,only: paraA
    implicit none
    !$acc routine seq
    real(kind=8), intent(in) :: dt, Snu, Sq, Qk, Qnu, fv(0:8), gv(0:4), rh, ux, uy, temp, fx, fy
    real(kind=8), intent(out) :: buffer(0:19)
    real(kind=8) :: m(0:8), meq(0:8), fm(0:8), n(0:4), neq(0:4), s(0:8), q(0:4)
    integer(kind=4) :: a

    ! 每次只处理一个来源节点
    ! 原两级示例（相邻粗细比2）：左细接收节点(129,129)使用粗模板(3:6,3:6)，第一次传入fv/gv为粗节点(3,3)的数据。
    ! 多级时dt与四个松弛率均取当前来源级；第1级dt=1，第2级dt=2，依此类推。

    ! 1. 将当前D2Q9分布函数转换为9个流场矩：m(0)为密度，m(3)、m(5)为动量矩。
    m(0)=fv(0)+fv(1)+fv(2)+fv(3)+fv(4)+fv(5)+fv(6)+fv(7)+fv(8)
    m(1)=-4.0d0*fv(0)-fv(1)-fv(2)-fv(3)-fv(4)+2.0d0*(fv(5)+fv(6)+fv(7)+fv(8))
    m(2)=4.0d0*fv(0)-2.0d0*(fv(1)+fv(2)+fv(3)+fv(4))+fv(5)+fv(6)+fv(7)+fv(8)
    m(3)=fv(1)-fv(3)+fv(5)-fv(6)-fv(7)+fv(8)
    m(4)=-2.0d0*fv(1)+2.0d0*fv(3)+fv(5)-fv(6)-fv(7)+fv(8)
    m(5)=fv(2)-fv(4)+fv(5)+fv(6)-fv(7)-fv(8)
    m(6)=-2.0d0*fv(2)+2.0d0*fv(4)+fv(5)+fv(6)-fv(7)-fv(8)
    m(7)=fv(1)-fv(2)+fv(3)-fv(4)
    m(8)=fv(5)-fv(6)+fv(7)-fv(8)

    ! 2. 将当前D2Q5分布函数转换为5个温度矩：n(0)为温度，n(1)、n(2)为通量相关矩。
    n(0)=gv(0)+gv(1)+gv(2)+gv(3)+gv(4)
    n(1)=gv(1)-gv(3)
    n(2)=gv(2)-gv(4)
    n(3)=-4.0d0*gv(0)+gv(1)+gv(2)+gv(3)+gv(4)
    n(4)=gv(1)-gv(2)+gv(3)-gv(4)

    ! 3. 根据同一来源节点的rho、u、v、T计算平衡矩；neq表示温度平衡矩。
    meq(0)=rh
    meq(1)=rh*(-2.0d0+3.0d0*(ux*ux+uy*uy))
    meq(2)=rh*(1.0d0-3.0d0*(ux*ux+uy*uy))
    meq(3)=rh*ux
    meq(4)=-rh*ux
    meq(5)=rh*uy
    meq(6)=-rh*uy
    meq(7)=rh*(ux*ux-uy*uy)
    meq(8)=rh*ux*uy
    neq=[temp,temp*ux,temp*uy,paraA*temp,0.0d0]

    ! 4. fm为单位时间的力矩F_m，不含dt，也不含碰撞系数(1-s/2)。
    ! fx/fy就是单位时间的F_x/F_y，直接构造F_m，不乘除dt。
    ! dt显式传入且等于当前级格距；第1、2、3级分别为1、2、4。
    fm=[0.0d0,6.0d0*(ux*fx+uy*fy),-6.0d0*(ux*fx+uy*fy),fx,-fx,fy,-fy,&
        2.0d0*(ux*fx-uy*fy),ux*fy+uy*fx]

    ! 5. 各矩对应的松弛率；密度、动量和温度守恒矩的松弛率置零，由单独保存的宏观量重建。
    s=[0.0d0,Snu,Snu,0.0d0,Sq,0.0d0,Sq,Snu,Snu]
    q=[0.0d0,Qk,Qk,Qnu,Qnu]

    ! 6. 前6项依次保存rho、u、v、T、F_x、F_y；所有网格直接保存单位时间的力。
    buffer(0:5)=[rh,ux,uy,temp,fx,fy]

    ! 7. buffer(6+a)=(s_a/dt)*(m_a-meq_a+dt*(F_m,a)/2)，dt为来源网格时间步。
    ! 碰撞源项为dt*(1-s_a/2)*F_m,a。
    ! 以第2级来源dt=2为例，括号内为m_a-meq_a+(F_m,a)；细块dt=1时为m_a-meq_a+(F_m,a)/2。
    ! 第2级来源的缩放系数为s(a)/2；buffer(6)、buffer(9)、buffer(11)对应守恒矩，有限场下为零。
    do a=0,8
        buffer(6+a)=s(a)/dt*(m(a)-meq(a)+0.5d0*dt*fm(a))
    enddo

    ! 8. buffer(15:19)保存5个缩放后的温度非平衡矩；第2级来源缩放系数为q(a)/2；一般级为q(a)/dt。
    ! buffer(15)对应温度守恒矩，有限场下为零；返回后外层再对这20项做空间插值。
    do a=0,4
        buffer(15+a)=q(a)/dt*(n(a)-neq(a))
    enddo
  end subroutine pack_interface_values

!===================================================================================================
  subroutine unpack_interface_values(dt,Snu,Sq,Qk,Qnu,buffer,fv,gv,rh,ux,uy,temp,fx,fy)
    use commondata,only: paraA
    implicit none    
    !可以从设备代码中调用
    !$acc routine seq
    real(kind=8), intent(in) :: dt, Snu, Sq, Qk, Qnu, buffer(0:19)
    real(kind=8), intent(out) :: fv(0:8), gv(0:4), rh, ux, uy, temp, fx, fy
    real(kind=8) :: m(0:8), meq(0:8), fm(0:8), n(0:4), neq(0:4), s(0:8), q(0:4)
    integer(kind=4) :: a

    ! 1. 恢复接收节点的宏观量；由调用处写入接收块的rho/u/v/T等数组。
    rh=buffer(0)
    ux=buffer(1)
    uy=buffer(2)
    temp=buffer(3)
    ! 直接恢复单位时间的F_x/F_y。
    fx=buffer(4)
    fy=buffer(5)
    ! 2. 直接计算流场和温度平衡矩；neq 是温度平衡矩，不是非平衡矩。
    meq(0)=rh
    meq(1)=rh*(-2.0d0+3.0d0*(ux*ux+uy*uy))
    meq(2)=rh*(1.0d0-3.0d0*(ux*ux+uy*uy))
    meq(3)=rh*ux
    meq(4)=-rh*ux
    meq(5)=rh*uy
    meq(6)=-rh*uy
    meq(7)=rh*(ux*ux-uy*uy)
    meq(8)=rh*ux*uy
    neq=[temp,temp*ux,temp*uy,paraA*temp,0.0d0]

    ! 3. buffer(4:5)存单位时间的F_x/F_y，按式(8a)得到单位时间的F_m。
    ! fm、fx、fy都不含dt；下方恢复矩时显式乘接收时间步。
    fm=[0.0d0,6.0d0*(ux*fx+uy*fy),-6.0d0*(ux*fx+uy*fy),fx,-fx,fy,-fy,&
        2.0d0*(ux*fx-uy*fy),ux*fy+uy*fx]

    ! 4. 用接收网格的松弛率和时间步恢复完整矩；第1级到第2级时dt=2，第2级到第3级时dt=4。
    s=[0.0d0,Snu,Snu,0.0d0,Sq,0.0d0,Sq,Snu,Snu]
    q=[0.0d0,Qk,Qk,Qnu,Qnu]
    m=meq
    n=neq
    do a=0,8
        ! 式(19)反解：m_a=meq_a+(dt/s_a)*buffer(6+a)-dt*F_m,a/2。
        if(s(a)>0.0d0) m(a)=meq(a)+dt/s(a)*buffer(6+a)-0.5d0*dt*fm(a)
    enddo
    ! 守恒矩按式(21)重建，不除以零松弛率；半步动量修正为dt*F/2。
    m(0)=rh
    m(3)=rh*ux-0.5d0*dt*fx
    m(5)=rh*uy-0.5d0*dt*fy
    ! 温度非守恒矩由 buffer(16:19) 恢复，温度守恒矩单独赋值。
    do a=1,4
        n(a)=neq(a)+dt/q(a)*buffer(15+a)
    enddo
    n(0)=temp

    ! 5. 将换算后的流场矩乘以 M 逆，得到接收网格的 D2Q9 分布 f。
    ! 初始化、细子步及同步时均返回 f/g；碰撞时再统一由分布函数求矩。
    fv(0)=m(0)/9.0d0-m(1)/9.0d0+m(2)/9.0d0    !这边是乘以M逆
    fv(1)=m(0)/9.0d0-m(1)/36.0d0-m(2)/18.0d0+m(3)/6.0d0-m(4)/6.0d0 &
        +m(7)/4.0d0
    fv(2)=m(0)/9.0d0-m(1)/36.0d0-m(2)/18.0d0 &
        +m(5)/6.0d0-m(6)/6.0d0-m(7)/4.0d0
    fv(3)=m(0)/9.0d0-m(1)/36.0d0-m(2)/18.0d0-m(3)/6.0d0+m(4)/6.0d0 &
        +m(7)/4.0d0
    fv(4)=m(0)/9.0d0-m(1)/36.0d0-m(2)/18.0d0 &
        -m(5)/6.0d0+m(6)/6.0d0-m(7)/4.0d0
    fv(5)=m(0)/9.0d0+m(1)/18.0d0+m(2)/36.0d0+m(3)/6.0d0+m(4)/12.0d0 &
        +m(5)/6.0d0+m(6)/12.0d0+m(8)/4.0d0
    fv(6)=m(0)/9.0d0+m(1)/18.0d0+m(2)/36.0d0-m(3)/6.0d0-m(4)/12.0d0 &
        +m(5)/6.0d0+m(6)/12.0d0-m(8)/4.0d0
    fv(7)=m(0)/9.0d0+m(1)/18.0d0+m(2)/36.0d0-m(3)/6.0d0-m(4)/12.0d0 &
        -m(5)/6.0d0-m(6)/12.0d0+m(8)/4.0d0
    fv(8)=m(0)/9.0d0+m(1)/18.0d0+m(2)/36.0d0+m(3)/6.0d0+m(4)/12.0d0 &
        -m(5)/6.0d0-m(6)/12.0d0-m(8)/4.0d0

    ! 6. 温度矩乘以 N 逆得到 D2Q5 的 g；温度碰撞仍在流场更新速度之后进行。
    gv(0)=0.2d0*n(0)-0.2d0*n(3)
    gv(1)=0.2d0*n(0)+0.5d0*n(1)+0.05d0*n(3)+0.25d0*n(4)
    gv(2)=0.2d0*n(0)+0.5d0*n(2)+0.05d0*n(3)-0.25d0*n(4)
    gv(3)=0.2d0*n(0)-0.5d0*n(1)+0.05d0*n(3)+0.25d0*n(4)
    gv(4)=0.2d0*n(0)-0.5d0*n(2)+0.05d0*n(3)-0.25d0*n(4)
  end subroutine unpack_interface_values

!===================================================================================================
  subroutine collision(nx,ny,dt,f,f_post,rho,u,v,Fx,Fy,T,Snu,Sq,gBeta &
#ifdef SideHeatedHa
      ,B2sigemarho &
#endif
      )

    use commondata,only: Tref
#ifdef SideHeatedHa
    use commondata,only: phi
#endif
    implicit none

    ! D2Q9 流场碰撞；nx/ny、场数组和松弛参数均对应当前区域。
    ! 本子程序只生成 f_post，不负责同级细区交换或粗细插值。
    ! Fx/Fy保存单位时间的力；fSource=(I-S/2)*F_m也不含dt，碰撞时显式乘dt。

    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), intent(inout) :: f(nx,ny,0:8), f_post(0:nx+1,0:ny+1,0:8)
    real(kind=8), intent(inout) :: u(nx,ny), v(nx,ny), T(nx,ny), rho(nx,ny), Fx(nx,ny), Fy(nx,ny)
    real(kind=8), intent(in) :: dt, Snu, Sq, gBeta
#ifdef SideHeatedHa
    real(kind=8), intent(in) :: B2sigemarho
#endif

    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: m(0:8), m_post(0:8), meq(0:8)
    real(kind=8) :: s(0:8)
    real(kind=8) :: fSource(0:8)

    !$acc parallel loop gang vector collapse(2) present(f, f_post, rho, u, v, Fx, Fy, T) async(1) &
    !$acc& firstprivate(dt) private(alpha, s, m, m_post, meq, fSource)
    do j=1,ny
        do i=1,nx

            ! 流场所有节点（含接收层）统一由分布 f 计算碰撞前矩 m。
            m(0)=f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,&
                7)+f(i,j,8)
            m(1)=-4.0d0*f(i,j,0)-f(i,j,1)-f(i,j,2)-f(i,j,3)-f(i,j,4)+2.0d0*(f(i,j,&
                5)+f(i,j,&
                6)+f(i,j,7)+f(i,j,8))
            m(2)=4.0d0*f(i,j,0)-2.0d0*(f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4))+f(i,j,&
                5)+f(i,j,&
                6)+f(i,j,7)+f(i,j,8)
            m(3)=f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)
            m(4)=-2.0d0*f(i,j,1)+2.0d0*f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)
            m(5)=f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)
            m(6)=-2.0d0*f(i,j,2)+2.0d0*f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)
            m(7)=f(i,j,1)-f(i,j,2)+f(i,j,3)-f(i,j,4)
            m(8)=f(i,j,5)-f(i,j,6)+f(i,j,7)-f(i,j,8)

            meq(0)=rho(i,j)
            meq(1)=rho(i,j)*(-2.0d0+3.0d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j)))
            meq(2)=rho(i,j)*(1.0d0-3.0d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j)))
            meq(3)=rho(i,j)*u(i,j)
            meq(4)=-rho(i,j)*u(i,j)
            meq(5)=rho(i,j)*v(i,j)
            meq(6)=-rho(i,j)*v(i,j)
            meq(7)=rho(i,j)*(u(i,j)*u(i,j)-v(i,j)*v(i,j))
            meq(8)=rho(i,j)*(u(i,j)*v(i,j))

            s(0)=0.0d0    !!s_{\rho}
            s(1)=Snu    !!s_{e}
            s(2)=Snu    !!s_{\epsilon}
            s(3)=0.0d0    !!s_{j}
            s(4)=Sq    !!s_{q}
            s(5)=0.0d0    !!s_{j}
            s(6)=Sq    !!s_{q}
            s(7)=Snu    !!s_{\nu}
            s(8)=Snu    !!s_{\nu}

            Fx(i,j)=0.0d0
            Fy(i,j)=rho(i,j)*gBeta*(T(i,j)-Tref)    !动量方程上的源项，即浮力项

#ifdef    SideHeatedHa
            Fx(i,j)=0.0d0+B2sigemarho*(v(i,j)*sin(phi)*cos(phi)-u(i,j)*sin(phi)*sin(phi))
            Fy(i,j)=rho(i,j)*gBeta*(T(i,j)-Tref)+rho(i,j)*B2sigemarho*(u(i,j)*sin(phi)*cos(phi)&
                -v(i,j)*cos(phi)*cos(phi))    !动量方程上的源项，即浮力项加磁场
#endif

            ! Huang式(5)：fSource(a)=(1-s(a)/2)*F_m,a，表示单位时间的碰撞源项。
            ! Fx/Fy和fSource都不含dt；系数6-3*s=6*(1-s/2)，时间步在m_post中显式相乘。
            ! 接口式(19)使用未乘(1-s/2)的F_m，不能直接代入这里的fSource。
            fSource(0)=0.0d0
            fSource(1)=(6.0d0-3.0d0*s(1))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
            fSource(2)=-(6.0d0-3.0d0*s(2))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
            fSource(3)=(1.0d0-0.5d0*s(3))*Fx(i,j)
            fSource(4)=-(1.0d0-0.5d0*s(4))*Fx(i,j)
            fSource(5)=(1.0d0-0.5d0*s(5))*Fy(i,j)
            fSource(6)=-(1.0d0-0.5d0*s(6))*Fy(i,j)
            fSource(7)=(2.0d0-s(7))*(u(i,j)*Fx(i,j)-v(i,j)*Fy(i,j))
            fSource(8)=(1.0d0-0.5d0*s(8))*(u(i,j)*Fy(i,j)+v(i,j)*Fx(i,j))    !这边是乘以M变到矩空间，然后再乘以1-1/2S修正

            do alpha=0,8
                m_post(alpha)=m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+dt*fSource(alpha)    !显式加入dt*(I-S/2)*F_m
            enddo

            f_post(i,j,0)=m_post(0)/9.0d0-m_post(1)/9.0d0+m_post(2)/9.0d0    !这边是乘以M逆
            f_post(i,j,&
                1)=m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+m_post(3)/6.0d0-&
                m_post(4)/6.0d0 &
                +m_post(7)/4.0d0
            f_post(i,j,2)=m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                +m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)/4.0d0
            f_post(i,j,&
                3)=m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0-m_post(3)/6.0d0+&
                m_post(4)/6.0d0 &
                +m_post(7)/4.0d0
            f_post(i,j,4)=m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                -m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)/4.0d0
            f_post(i,j,&
                5)=m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+&
                m_post(4)/12.0d0 &
                +m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)/4.0d0
            f_post(i,j,&
                6)=m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-&
                m_post(4)/12.0d0 &
                +m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)/4.0d0
            f_post(i,j,&
                7)=m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-&
                m_post(4)/12.0d0 &
                -m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)/4.0d0
            f_post(i,j,&
                8)=m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+&
                m_post(4)/12.0d0 &
                -m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)/4.0d0

        enddo
    enddo

    return
  end subroutine collision

!===================================================================================================
  subroutine streaming(nx,ny,f,f_post)    !先迁移，再边界处理

    use commondata,only: ex,ey
    implicit none

    ! 从上游节点读取 f_post；相邻细区的数据在调用前已复制到迁移外圈。
    ! 触及真实墙面的入射分布，随后由 bounceback 修正。

    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), intent(inout) :: f(nx,ny,0:8), f_post(0:nx+1,0:ny+1,0:8)

    integer(kind=4) :: i, j
    integer(kind=4) :: ip, jp
    integer(kind=4) :: alpha

    !$acc parallel loop gang vector collapse(2) present(f, f_post, ex, ey) async(1) private(alpha, ip, jp)
    do j=1,ny
        do i=1,nx
            do alpha=0,8    !上游格点索引：fα(i,j) <- f_postα(i-exα, j-eyα)
                ip=i-ex(alpha)    !边界附近 (ip/jp 可能为 0 或 nx+1/ny+1)，需在 bounceback/周期边界处理中覆盖修正边界分布
                jp=j-ey(alpha)    !ghost 层在初始化中为 0，保证不会出现未初始化垃圾值

                f(i,j,alpha)=f_post(ip,jp,alpha)
            enddo
        enddo
    enddo

    return
  end subroutine streaming

!===================================================================================================
  subroutine bounceback(nx,ny,f,f_post,leftWall,rightWall,bottomWall,topWall)

    implicit none

    ! 四个墙面标志指定本区实际接触的物理边界，仅在对应边执行边界处理。
    ! 细区拼接处由直接迁移连接，粗细人工边缘由接口重建，均不在这里反弹。

    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), intent(inout) :: f(nx,ny,0:8), f_post(0:nx+1,0:ny+1,0:8)
    logical, intent(in) :: leftWall, rightWall, bottomWall, topWall

    integer(kind=4) :: i, j
    ! integer(kind=4) :: alpha


#ifdef VerticalWallsNoslip

    !$acc parallel loop gang vector present(f, f_post) async(1)
    do j=1,ny    !速度边界垂直边界静止壁无滑移，直接反弹，方向相反
        !Left side (i=1)
        if(leftWall) f(1,j,1)=f_post(1,j,3)
        if(leftWall) f(1,j,5)=f_post(1,j,7)
        if(leftWall) f(1,j,8)=f_post(1,j,6)

        !Right side (i=nx)
        if(rightWall) f(nx,j,3)=f_post(nx,j,1)
        if(rightWall) f(nx,j,6)=f_post(nx,j,8)
        if(rightWall) f(nx,j,7)=f_post(nx,j,5)
    enddo
#endif

#ifdef HorizontalWallsNoslip

    !$acc parallel loop gang vector present(f, f_post) async(1)
    do i=1,nx    !速度边界水平边界无滑移，直接反弹，方向相反
        !Bottom side (j=1)
        if(bottomWall) f(i,1,2)=f_post(i,1,4)
        if(bottomWall) f(i,1,5)=f_post(i,1,7)
        if(bottomWall) f(i,1,6)=f_post(i,1,8)

        !Top side (j=ny)
        if(topWall) f(i,ny,4)=f_post(i,ny,2)
        if(topWall) f(i,ny,7)=f_post(i,ny,5)
        if(topWall) f(i,ny,8)=f_post(i,ny,6)
    enddo
#endif

    return
  end subroutine bounceback

!===================================================================================================
  subroutine macro(nx,ny,dt,f,rho,u,v,Fx,Fy)

    implicit none

    ! Fx/Fy为单位时间的力，速度恢复显式包含dt*F/2。

    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), intent(in) :: dt
    real(kind=8), intent(in) :: f(nx,ny,0:8), Fx(nx,ny), Fy(nx,ny)
    real(kind=8), intent(inout) :: rho(nx,ny), u(nx,ny), v(nx,ny)

    integer(kind=4) :: i, j

    !$acc parallel loop gang vector collapse(2) present(f, rho, u, v, Fx, Fy) firstprivate(dt) async(1)
    do j=1,ny
        do i=1,nx
            rho(i,j)=f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,&
                6)+f(i,j,&
                7)+f(i,j,8)
            u(i,j)=(f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)+0.5d0*dt*Fx(i,&
                j))/rho(i,j)    !含力LBM的半步动量修正：rho*u = Σ f e + 0.5*dt*F，对应Guo forcing的二阶定义
            v(i,j)=(f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)+0.5d0*dt*Fy(i,&
                j))/rho(i,j)
        enddo
    enddo

    return
  end subroutine macro

!===================================================================================================
  subroutine collisionT(nx,ny,g,g_post,u,v,T,Qk,Qnu)

    use commondata,only: paraA
    implicit none

    ! Qk/Qnu 已按当前区域格距设置，温度碰撞使用本步更新后的速度。

    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), intent(inout) :: g(nx,ny,0:4), g_post(0:nx+1,0:ny+1,0:4)
    real(kind=8), intent(in) :: u(nx,ny), v(nx,ny), T(nx,ny), Qk, Qnu

    integer(kind=4) :: i, j, a
    real(kind=8) :: n(0:4), n_post(0:4), neq(0:4)
    real(kind=8) :: q(0:4)

    !$acc parallel loop gang vector collapse(2) present(g, g_post, u, v, T) async(1) &
    !$acc& private(a, n, neq, q, n_post)
    do j=1,ny
        do i=1,nx

            ! 温度所有节点（含接收层）统一由分布 g 计算碰撞前矩 n。
            n(0)=g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
            n(1)=g(i,j,1)-g(i,j,3)
            n(2)=g(i,j,2)-g(i,j,4)
            n(3)=-4.0d0*g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
            n(4)=g(i,j,1)-g(i,j,2)+g(i,j,3)-g(i,j,4)

            neq(0)=T(i,j)
            neq(1)=T(i,j)*u(i,j)
            neq(2)=T(i,j)*v(i,j)
            neq(3)=T(i,j)*paraA
            neq(4)=0.0d0

            q(0)=0.0d0
            q(1)=Qk
            q(2)=Qk
            q(3)=Qnu
            q(4)=Qnu

            n_post(0)=n(0)-q(0)*(n(0)-neq(0))
            n_post(1)=n(1)-q(1)*(n(1)-neq(1))
            n_post(2)=n(2)-q(2)*(n(2)-neq(2))
            n_post(3)=n(3)-q(3)*(n(3)-neq(3))
            n_post(4)=n(4)-q(4)*(n(4)-neq(4))

            g_post(i,j,0)=0.2d0*n_post(0)-0.2d0*n_post(3)
            g_post(i,j,1)=0.2d0*n_post(0)+0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
            g_post(i,j,2)=0.2d0*n_post(0)+0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4)
            g_post(i,j,3)=0.2d0*n_post(0)-0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
            g_post(i,j,4)=0.2d0*n_post(0)-0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4)
        enddo
    enddo

    return
  end subroutine collisionT

!===================================================================================================
  subroutine streamingT(nx,ny,g,g_post)

    use commondata,only: ex,ey
    implicit none

    ! 从上游节点读取 g_post；相邻细区的数据在调用前已复制到迁移外圈。

    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), intent(inout) :: g(nx,ny,0:4), g_post(0:nx+1,0:ny+1,0:4)

    integer(kind=4) :: i, j
    integer(kind=4) :: ip, jp
    integer(kind=4) :: alpha

    !$acc parallel loop gang vector collapse(2) present(g, g_post, ex, ey) async(1) private(alpha, ip, jp)
    do j=1,ny
        do i=1,nx
            do alpha=0,4
                ip=i-ex(alpha)
                jp=j-ey(alpha)

                g(i,j,alpha)=g_post(ip,jp,alpha)
            enddo
        enddo
    enddo

    return
  end subroutine streamingT

!===================================================================================================
  subroutine bouncebackT(nx,ny,g,g_post,leftWall,rightWall,bottomWall,topWall)

    use commondata,only: Thot,Tcold,paraA,omegaT
    implicit none

    ! 按温度边界宏处理恒温、绝热或周期边界，仅作用于四个标志指定的物理边界。

    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), intent(inout) :: g(nx,ny,0:4), g_post(0:nx+1,0:ny+1,0:4)
    logical, intent(in) :: leftWall, rightWall, bottomWall, topWall

    integer(kind=4) :: i, j
    !integer(kind=4) :: alpha


#ifdef VerticalWallsConstT

    !$acc parallel loop gang vector present(g, g_post, omegaT) async(1)
    do j=1,ny
        !Left boundary
        if(leftWall) g(1,j,1)=-g_post(1,j,3)+(4.0d0+paraA)/10.0d0*Thot
        !Right boundary
        if(rightWall) g(nx,j,3)=-g_post(nx,j,1)+(4.0d0+paraA)/10.0d0*Tcold
    enddo
#endif

#ifdef VerticalWallsAdiabatic

    !$acc parallel loop gang vector present(g, g_post) async(1)
    do j=1,ny
        !Left boundary
        if(leftWall) g(1,j,1)=g_post(1,j,3)

        !Right boundary
        if(rightWall) g(nx,j,3)=g_post(nx,j,1)
    enddo
#endif

#ifdef HorizontalWallsAdiabatic

    !$acc parallel loop gang vector present(g, g_post) async(1)
    do i=1,nx
        !Bottom side
        if(bottomWall) g(i,1,2)=g_post(i,1,4)

        !Top side
        if(topWall) g(i,ny,4)=g_post(i,ny,2)
    enddo
#endif

#ifdef HorizontalWallsConstT

    !$acc parallel loop gang vector present(g, g_post, omegaT) async(1)
    do i=1,nx
        !Bottom side
        if(bottomWall) g(i,1,2)=-g_post(i,1,4)+(4.0d0+paraA)/10.0d0*Thot
        !Top side
        if(topWall) g(i,ny,4)=-g_post(i,ny,2)+(4.0d0+paraA)/10.0d0*Tcold
    enddo
#endif

    return
  end subroutine bouncebackT

!===================================================================================================
  subroutine macroT(nx,ny,g,T)

    implicit none

    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), intent(in) :: g(nx,ny,0:4)
    real(kind=8), intent(inout) :: T(nx,ny)

    integer(kind=4) :: i, j

    !$acc parallel loop gang vector collapse(2) present(g, T) async(1)
    do j=1,ny
        do i=1,nx
            T(i,j)=g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
        enddo
    enddo

    return
  end subroutine macroT

!===================================================================================================
  subroutine section_weights(q,n,first,w,dw)

    implicit none

    ! q 是调用方编号体系中的连续下标：粗区使用粗网格连续下标，fine_section_x/y 使用全局细网格连续下标。
    ! n 为该编号体系的节点数；first 是同一编号体系中四点模板的首个整数下标。
    ! w 用于函数值，dw 是对连续下标的导数权重；换成对全局物理坐标的导数时还要除以该网格 dx。

    real(kind=8), intent(in) :: q
    integer(kind=4), intent(in) :: n
    integer(kind=4), intent(out) :: first
    real(kind=8), intent(out) :: w(4), dw(4)
    real(kind=8) :: z, term
    integer(kind=4) :: a, k, l

    first=max(1,min(floor(q)-1,n-3))  !floor 向下取整
    z=q-dble(first) !把目标位置换算成相对于模板首节点的位置

    call lagrange_weights(z,w) !四点拉格朗日插值系数
    dw=0.0d0   !温度梯度中，四个插值系数对位置的导数
    ! 对拉格朗日基函数的求导
    do a=0,3
        do k=0,3
            if(k==a) cycle
            term=1.0d0/dble(a-k) !拉格朗日基函数中第 k 个公式的求导，分子会出现 1
            do l=0,3             !乘上其余没有求导的因子
                if(l/=a .AND. l/=k) term=term*(z-dble(l))/dble(a-l)
            enddo

            dw(a+1)=dw(a+1)+term !把乘积求导产生的三项加起来
        enddo
    enddo
  end subroutine section_weights
