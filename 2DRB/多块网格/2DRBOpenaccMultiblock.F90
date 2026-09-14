!=============================================================
!!!    注释区，代码描述
!!!    二维浮力驱动自然对流 OpenACC 静态多块网格版本
!!!    块内 D2Q9/D2Q5 算法保持不变
!!!    LBM方法
!!!    MRT-LBE
!=============================================================

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

!算法切换
!启用 M1G 修正；注释掉则不使用 useG 相关修正
#define EnableUseG
!启用旧温度算法
!#define EnableLegacyThermalScheme

!   温度算法宏的选择
#if defined(EnableUseG) && defined(EnableLegacyThermalScheme)
#error "Choose only one thermal scheme: EnableUseG or EnableLegacyThermalScheme"
#endif
#if !defined(EnableUseG) && !defined(EnableLegacyThermalScheme)
#error "Define one thermal scheme: EnableUseG or EnableLegacyThermalScheme"
#endif

!   自定义宏结束
!=============================================================

!=============================================================
!   全局模块
    module commondata
        ! ieee_arithmetic 是 Fortran 标准内置模块，提供 IEEE 浮点数运算相关功能。
        ! intrinsic 明确指定使用编译器提供的内置模块，避免与同名的非内置模块混淆。
        ! only 只引入 ieee_is_finite：有限值返回 .true.，NaN 或正负无穷大返回 .false.。
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        implicit none

        !===============================================================================================
        ! 网格参数：节点编号、粗细格距与接口重叠
        !===============================================================================================
        ! nx,ny 为按最细格距计的全域长度
        ! 物理壁面在 x=0、nx 和 y=0、ny；最细块 dx=dt=1，中心粗块 dx=dt=refineRatio。
        integer(kind=4), parameter :: nx=1024, ny=1024    ! 手动设置最细网格等效分辨率，修改后重新编译
        integer(kind=4), parameter :: refineRatio=2    ! 粗/细格距及时间步之比；逐级二分取2、4、8等2的幂，1为单块
        integer(kind=4), parameter :: fineLayerCellsLeftRight=128, fineLayerCellsBottomTop=128    ! 左右/上下细块的基准边界：从对应墙向内数的细节点编号
        ! 编号128表示距墙127.5；默认左/下细块基准为127.5，右/上为896.5，不是最终的积分分界。
        ! 多块模式 nx、ny 须被 refineRatio 整除；上述节点编号为正整数，无须整除粗细比。
        ! 基准距墙为编号减0.5，至少 (overlapCells+1)*refineRatio；中心宽高须为正，插值模板按实际节点检查。
        ! 粗块以左/下基准为锚点，右/上基准向外补齐完整粗格距；默认补到897.5，细块基准不动。
        ! 积分分界取粗块共址基准：默认中心为[127.5,897.5]^2，相邻细块及角部接到同一分界，面积不重复。
        integer(kind=4), parameter :: overlapCells=2    ! 各块从自身基准边界向外延伸的粗格距数；总重叠宽度还包括基准边界差
        integer(kind=4), parameter :: interfaceSkin=2    ! 原流场->温度场分步推进需要两层人工边界；这些节点在接收后参与碰撞

        !===============================================================================================
        ! 是否从本程序检查点精确续算
        !===============================================================================================
        integer(kind=4), parameter :: loadInitField=0    ! 0: 重新初始化；1: 按 latest.meta 指向的检查点续算

        !===============================================================================================
        ! 无量纲参数及物理壁面温度
        !===============================================================================================
        real(kind=8), parameter :: Rayleigh=1.0d7    ! 手动设置瑞利数，修改后重新编译
        real(kind=8), parameter :: Prandtl=0.7d0
        real(kind=8), parameter :: Mach=0.1d0
        real(kind=8), parameter :: Thot=0.5d0
        real(kind=8), parameter :: Tcold=-0.5d0
        real(kind=8), parameter :: Tref=0.5d0*(Thot+Tcold)
        real(kind=8), parameter :: pi=acos(-1.0d0)
#ifdef SideHeatedCell
        real(kind=8), parameter :: lengthUnit=dble(nx)
#else
        real(kind=8), parameter :: lengthUnit=dble(ny)
#endif

        !===============================================================================================
        ! 输运系数：沿用原程序的最细格子单位
        !===============================================================================================
        ! 以下参数全部沿用原文件的最细格子单位，不按块重新定义 Ra、Pr、Ma 或 paraA。
        real(kind=8), parameter :: tauf=0.5d0+Mach*lengthUnit*sqrt(3.0d0*Prandtl/Rayleigh)
        real(kind=8), parameter :: viscosity=(tauf-0.5d0)/3.0d0
        real(kind=8), parameter :: diffusivity=viscosity/Prandtl

        ! 浮力项及无量纲时间、速度
        real(kind=8), parameter :: gBeta1=Rayleigh*viscosity*diffusivity/lengthUnit
        real(kind=8), parameter :: gBeta=gBeta1/lengthUnit/lengthUnit
        real(kind=8), parameter :: timeUnit=sqrt(lengthUnit/gBeta)
        real(kind=8), parameter :: velocityUnit=sqrt(gBeta*lengthUnit)
        real(kind=8), parameter :: velocityScaleCompare=lengthUnit/diffusivity

        ! 流场和温度场的多松弛系数
        real(kind=8), parameter :: Snu=1.0d0/tauf
        real(kind=8), parameter :: Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)
        real(kind=8), parameter :: paraA=20.0d0*sqrt(3.0d0)*diffusivity-4.0d0
#ifdef EnableLegacyThermalScheme
        real(kind=8), parameter :: Qk=3.0d0-sqrt(3.0d0)
        real(kind=8), parameter :: Qnu=4.0d0*sqrt(3.0d0)-6.0d0
        real(kind=8), parameter :: thermalGeqCoeff=10.0d0/(4.0d0+paraA)
        real(kind=8), parameter :: thermalA=paraA
#else
        real(kind=8), parameter :: taug=0.5d0+(tauf-0.5d0)/Prandtl
        real(kind=8), parameter :: Qk=1.0d0/taug
        real(kind=8), parameter :: Qnu=1.0d0
        real(kind=8), parameter :: thermalGeqCoeff=3.0d0
        real(kind=8), parameter :: thermalA=-2.0d0/3.0d0
#endif
#ifdef SideHeatedHa
        real(kind=8), parameter :: Ha=20.0d0
        real(kind=8), parameter :: phi=0.0d0*pi/180.0d0
        real(kind=8), parameter :: B2sigemarho=Ha**2*viscosity/lengthUnit**2
#endif

        !===============================================================================================
        ! 收敛阈值、计算时间及输出间隔
        !===============================================================================================
        real(kind=8), parameter :: epsU=1.0d-7
        real(kind=8), parameter :: epsT=1.0d-7
#ifdef steadyFlow
        real(kind=8), parameter :: outputSnapshotInterval=10.0d0
        integer(kind=4), parameter :: itc_max=20000000
#else
        real(kind=8), parameter :: outputSnapshotInterval=0.5d0
        real(kind=8), parameter :: unsteadyRunDuration=1000.0d0
        real(kind=8), parameter :: unsteadyAverageStartTf=0.5d0*unsteadyRunDuration
        real(kind=8), parameter :: unsteadyAverageEndTf=unsteadyRunDuration
        real(kind=8), parameter :: unsteadyAverageMidTf=0.5d0*(unsteadyAverageStartTf+unsteadyAverageEndTf)
        integer(kind=4), parameter :: itc_max=max(1,ceiling(unsteadyRunDuration*timeUnit))
#endif
        real(kind=8), parameter :: reloadFileInterval=100.0d0
        real(kind=8), parameter :: outputPltFileInterval=100.0d0
        integer(kind=4), parameter :: outputSnapshotFile=1
        integer(kind=4), parameter :: outputPltFile=1
        integer(kind=4), parameter :: outputReloadFile=1

        ! 输出文件命名与格式版本
        character(*), parameter :: settingsFile='SimulationSettings2DOpenaccMultiblock.txt'
        character(*), parameter :: snapshotFilePrefix='buoyancyCavity2DOpenaccMultiblockSnapshot'
        character(*), parameter :: pltFolderPrefix='buoyancyCavity2DOpenaccMultiblockTecplot'
        character(*), parameter :: reloadFilePrefix='reloadFile2DOpenaccMultiblock'
        character(*), parameter :: historyFile='NuRe_2DOpenaccMultiblock.dat'
        character(16), parameter :: restartMagic='MB2DRESTART0006'
        character(16), parameter :: snapshotMagic='MB2DSNAPSHOT0002'

        !===============================================================================================
        ! 格子方向、块数据与接口交换数据
        !===============================================================================================
        integer(kind=4), parameter :: packetSize=22
        integer(kind=4), parameter :: maxBlocks=5
        integer(kind=4) :: ex(0:8)=[0,1,0,-1,0,1,-1,-1,1], ey(0:8)=[0,0,1,0,-1,1,1,-1,-1]
        real(kind=8) :: omega(0:8), omegaT(0:4)
        integer(kind=4) :: nBlocks, itc=0, snapshotFileNum=0, pltFileNum=0
        integer(kind=4) :: nextSample=1, nextReload=1, nextPlt=1
        real(kind=8) :: errorU=100.0d0, errorT=100.0d0

        type :: gridBlock
            integer(kind=4) :: ni,nj,ilo,ihi,jlo,jhi,nh
            real(kind=8) :: h,x0,y0,sn,sq,qk,qn,gb
            ! x0,y0 是首节点减去本块半格距的虚拟面，不再等于统计分区边界。
            real(kind=8) :: ownedBox(4)    ! 不重叠物理积分分区：xmin,xmax,ymin,ymax
            real(kind=8) :: baseBox(4)    ! 添加重叠层前的基准边界；粗块右/上已向外对齐，物理墙仍记墙面坐标
            real(kind=8), allocatable :: dxWeight(:),dyWeight(:)    ! 节点控制区与 ownedBox 的交集长度
            logical :: wall(4)    ! 左、右、下、上；仅真实物理壁面施加原来的 BB/ABB
            real(kind=8), allocatable :: f(:,:,:),f_post(:,:,:),g(:,:,:),g_post(:,:,:)
            real(kind=8), allocatable :: rho(:,:),u(:,:),v(:,:),T(:,:),Fx(:,:),Fy(:,:),Bx_prev(:,:),By_prev(:,:)
            ! packet: rho,u,v,T,Fx/dt,Fy/dt,Kf(0:8),Kg(0:4),dBx/dt,dBy/dt。
            ! K=S/dt*(m-meq+F_lattice/2)。细块仅存当前面向接口的状态，粗块保留三个时间层。
            real(kind=8), allocatable :: p(:,:,:,:)
#ifdef steadyFlow
            real(kind=8), allocatable :: up(:,:),vp(:,:),Tp(:,:)
#endif
        end type

        type :: blockLink
            integer(kind=4) :: receiver,donor,count
            integer(kind=4), allocatable :: ti(:),tj(:),si(:),sj(:)
            logical, allocatable :: coincident(:)    ! 同坐标节点免除空间插值，仍作矩重标定
            real(kind=8), allocatable :: wx(:,:),wy(:,:),values(:,:)
        end type

        type(gridBlock) :: blocks(maxBlocks)
        type(blockLink) :: links(maxBlocks*maxBlocks)
        integer(kind=4) :: nLinks=0

    contains

    !===============================================================================================
    ! 初始化：构造计算块、积分分区及接口连接
    subroutine initial()
        integer(kind=4) :: b,k,overlap
        real(kind=8) :: totalArea,xLeft,xRight,yBottom,yTop,xSplit,ySplit

        if (refineRatio<1) error stop 'refineRatio must be a positive integer'
        if (min(nx,ny)<8) error stop 'At least 8 cells per direction are required'
        if (loadInitField/=0 .and. loadInitField/=1) error stop 'loadInitField must be 0 or 1'
        if (min(outputSnapshotInterval,reloadFileInterval,outputPltFileInterval)*timeUnit<dble(refineRatio)) &
            error stop 'Output intervals must be at least one synchronized coarse step'
#ifdef EnableLegacyThermalScheme
        if (paraA<=-4.0d0 .or. paraA>=1.0d0) error stop 'Legacy paraA must be in (-4,1)'
#endif
        omega(0) = 4.0d0/9.0d0
        omega(1:4) = 1.0d0/9.0d0
        omega(5:8) = 1.0d0/36.0d0
        omegaT(0) = (1.0d0-thermalA)/5.0d0
        omegaT(1:4) = (thermalA+4.0d0)/20.0d0
        if (refineRatio==1) then
            nBlocks = 1
            call make_block(blocks(1),0.0d0,dble(nx),0.0d0,dble(ny),1,0,0)
        else
#if defined(VerticalWallsPeriodicalU) || defined(VerticalWallsPeriodicalT)
            error stop 'Multiblock periodic sides require a periodic block topology; use wall BCs or refineRatio=1'
#endif
            overlap = overlapCells*refineRatio
            if (interfaceSkin<2) error stop 'Split flow/thermal advance requires at least two interface layers'
            ! 目标细块最内侧缓冲点必须位于来源粗块可用区内；四点模板另由 donor_stencil 检查。
            ! 2*overlap >= interfaceSkin*refineRatio + interfaceSkin-1；默认比值2、两层缓冲时 overlapCells>=2。
            if (2*overlap<interfaceSkin*refineRatio+interfaceSkin-1) &
                error stop 'Overlap is too narrow for the pre-collision interface layers'
            if (mod(nx,refineRatio)/=0 .or. mod(ny,refineRatio)/=0) &
                error stop 'nx,ny must be multiples of refineRatio'
            xLeft = dble(fineLayerCellsLeftRight)-0.5d0
            xRight = dble(nx)-xLeft
            yBottom = dble(fineLayerCellsBottomTop)-0.5d0
            yTop = dble(ny)-yBottom
            if (min(xLeft,yBottom)<dble(overlap+refineRatio)) &
                error stop 'Refined wall layer must exceed the overlap width by at least one coarse cell'
            if (xRight<=xLeft .or. yTop<=yBottom) &
                error stop 'The coarse core must have positive width and height'
            nBlocks = 5
            ! 计算块：中心、下、上、左、右。细块基准保持用户指定的接口细节点坐标。
            ! make_block 先确定各块 baseBox，再向邻块延伸 overlap；统计分区随后单独设置。
            ! 最后一个 2 表示粗块历史下标 0:2（三个时间层），与粗细比无关。
            call make_block(blocks(1),xLeft,xRight,yBottom,yTop,refineRatio,overlap,2)
            call make_block(blocks(2),0.0d0,dble(nx),0.0d0,yBottom,1,overlap,0)
            call make_block(blocks(3),0.0d0,dble(nx),yTop,dble(ny),1,overlap,0)
            call make_block(blocks(4),0.0d0,xLeft,yBottom,yTop,1,overlap,0)
            call make_block(blocks(5),xRight,dble(nx),yBottom,yTop,1,overlap,0)
            ! 统计分界取粗细共址的粗块基准边界，默认右/上为897.5，而细块基准仍为896.5。
            ! 上细块及左右细块共用 ySplit，右细块与中心共用 xSplit，角部不重复也不漏算。
            xSplit = blocks(1)%baseBox(2)
            ySplit = blocks(1)%baseBox(4)
            blocks(1)%ownedBox = blocks(1)%baseBox
            blocks(3)%ownedBox = [0.0d0,dble(nx),ySplit,dble(ny)]
            blocks(4)%ownedBox = [0.0d0,xLeft,yBottom,ySplit]
            blocks(5)%ownedBox = [xSplit,dble(nx),yBottom,ySplit]
        endif
        totalArea = 0.0d0
        do b = 1,nBlocks
            ! 以最终统计分区统一生成权重；中心两端均为粗节点，自动得到端点半权重的梯形积分。
            associate(bl=>blocks(b))
                call integration_weights(bl%ni,bl%x0,bl%h,bl%ownedBox(1),bl%ownedBox(2),bl%dxWeight,bl%ilo,bl%ihi)
                call integration_weights(bl%nj,bl%y0,bl%h,bl%ownedBox(3),bl%ownedBox(4),bl%dyWeight,bl%jlo,bl%jhi)
            end associate
            call initial_block(blocks(b))
            totalArea = totalArea+sum(blocks(b)%dxWeight)*sum(blocks(b)%dyWeight)
        enddo
        if (abs(totalArea-dble(nx)*ny)>1.0d-8) error stop 'Block ownership does not tile the physical domain'
        if (loadInitField==1) call read_restart()
        call build_links()
        open(newunit = k,file = settingsFile,status = 'replace')
        write(k,*) 'Parent: uniform-grid 2DRBOpenacc.F90; D2Q9 MRT / D2Q5 MRT unchanged inside blocks'
        write(k,*) 'Huang (2014): force-aware MRT scaling and four-point Lagrange; original D2Q5 adaptation.'
        write(k,*) 'Yu (2002): aligned two-way coupling and acoustic subcycling; cubic spline is not used.'
        write(k,*) 'Chen (2016): pre-collision buffer states; received f,g participate in the next collision.'
        write(k,*) 'Fine-equivalent nx,ny; refinement ratio:',nx,ny,refineRatio
        write(k,*) 'Rayleigh, Prandtl, Mach:',Rayleigh,Prandtl,Mach
        write(k,*) 'Fine tauf, Snu, Sq, Qk, Qnu, thermalA:',tauf,Snu,Sq,Qk,Qnu,thermalA
        write(k,*) 'Fine viscosity, diffusivity, gBeta, timeUnit:',viscosity,diffusivity,gBeta,timeUnit
        write(k,*) 'Interface fine-node indices from each wall; overlap in coarse spacings:', &
            fineLayerCellsLeftRight,fineLayerCellsBottomTop,overlapCells
        write(k,*) 'Owned physical area:',totalArea
        write(k,*) 'Geometry: align coarse right/top base boundaries outward, then extend each base by the overlap.'
        write(k,*) 'Fine bases stay at requested fine-node interfaces; integration partitions meet at coarse base nodes.'
        write(k,*) 'Node alignment: coarse nodes remain a subset of fine nodes; the integer phase depends on block origin.'
        write(k,*) 'Integration: clipped nodal control areas; shared coordinates do not duplicate physical area.'
        write(k,*) 'All nonconserved relaxation times satisfy h*(1/s-1/2)=constant.'
        write(k,*) 'Time: coarse prediction; prepare fine buffers BEFORE each substep; endpoint synchronization.'
        write(k,*) 'Buffers: two local layers for original split flow/thermal sequence; donors exclude these layers.'
        write(k,*) 'First coarse interval uses linear startup; subsequent intervals use three-time Lagrange interpolation.'
        write(k,*) 'Output clocks are absolute and rounded UP to a synchronized coarse step; time columns contain actual times.'
        write(k,*) 'Restart stores both coarse history and thermal B_prev; uniform-grid restart files are incompatible.'
#ifdef EnableUseG
        write(k,*) 'Thermal scheme: original D2Q5 EnableUseG, including transferred B_prev history.'
#else
        write(k,*) 'Thermal scheme: original legacy D2Q5; paraA is fixed across blocks.'
#endif
        do b = 1,nBlocks
            write(k,*) 'block,ni,nj,h,x0,y0,owned ilo,ihi,jlo,jhi:',b,blocks(b)%ni,blocks(b)%nj, &
                blocks(b)%h,blocks(b)%x0,blocks(b)%y0,blocks(b)%ilo,blocks(b)%ihi,blocks(b)%jlo,blocks(b)%jhi
            write(k,*) 'Snu,Sq,Qk,Qnu,gBeta:',blocks(b)%sn,blocks(b)%sq,blocks(b)%qk,blocks(b)%qn,blocks(b)%gb
            write(k,*) 'Owned physical rectangle:',blocks(b)%ownedBox
            write(k,*) 'Base rectangle before overlap:',blocks(b)%baseBox
            write(k,*) 'Computed first/last node x,y:',blocks(b)%x0+0.5d0*blocks(b)%h, &
                blocks(b)%x0+(dble(blocks(b)%ni)-0.5d0)*blocks(b)%h, &
                blocks(b)%y0+0.5d0*blocks(b)%h,blocks(b)%y0+(dble(blocks(b)%nj)-0.5d0)*blocks(b)%h
        enddo
        do b = 1,nLinks
            write(k,*) 'receiver, donor, direct nodes, interpolated nodes:',links(b)%receiver,links(b)%donor, &
                count(links(b)%coincident),links(b)%count-count(links(b)%coincident)
        enddo
        close(k)
        if (loadInitField==0) then
            open(newunit = k,file = historyFile,status = 'replace')
            write(k,'(a)') '# t_ff NuVolAvg ReVolRMS Nu_hot Nu_cold Nu_middle mass meanT Tmin Tmax rhoMin rhoMax'
            close(k)
        else
            call check_history()
        endif
    end subroutine initial

    !===============================================================================================
    ! 计算块：先对齐粗块基准边界，再向外延伸重叠层
    subroutine make_block(b,xlo,xhi,ylo,yhi,spacing,overlap,nh)
        type(gridBlock), intent(out) :: b
        real(kind=8), intent(in) :: xlo,xhi,ylo,yhi
        integer(kind=4), intent(in) :: spacing,overlap,nh
        real(kind=8) :: xb,xe,yb,ye

        b%h = dble(spacing)
        b%ownedBox = [xlo,xhi,ylo,yhi]
        b%baseBox = b%ownedBox
        if (spacing>1) then
            ! 先补齐粗块基准跨度，再添加重叠层；细块基准不动，统计分区由 initial 单独确定。
            ! 默认 xlo=127.5,xhi=896.5,h=2：粗右基准为897.5，加4后末节点为901.5；y同理。
            b%baseBox(2) = xlo+dble(ceiling((xhi-xlo)/b%h))*b%h
            b%baseBox(4) = ylo+dble(ceiling((yhi-ylo)/b%h))*b%h
            if (min(xlo,ylo)-dble(overlap)<=0.5d0 .or. &
                b%baseBox(2)+dble(overlap)>=dble(nx)-0.5d0 .or. &
                b%baseBox(4)+dble(overlap)>=dble(ny)-0.5d0) &
                error stop 'Aligned coarse block must remain inside the fine wall layers'
        endif
        ! 各块从自身基准边界延伸；只有细块接触物理墙，在墙处保留半个细格距。
        xb = max(0.5d0,b%baseBox(1)-dble(overlap))
        xe = min(dble(nx)-0.5d0,b%baseBox(2)+dble(overlap))
        yb = max(0.5d0,b%baseBox(3)-dble(overlap))
        ye = min(dble(ny)-0.5d0,b%baseBox(4)+dble(overlap))
        ! xb、yb 已是首节点坐标，不再加 0.5；整数粗细比保证粗节点属于细节点子集。
        b%x0 = xb-0.5d0*b%h
        b%y0 = yb-0.5d0*b%h
        ! 基准跨度已补齐，overlap 也是粗格距整数倍；末节点必须恰好到达 xe、ye。
        b%ni = nint((xe-xb)/b%h)+1
        b%nj = nint((ye-yb)/b%h)+1
        b%nh = nh
        if (abs(xb+dble(b%ni-1)*b%h-xe)>1.0d-10 .or. abs(yb+dble(b%nj-1)*b%h-ye)>1.0d-10) &
            error stop 'Computed block span must contain whole grid spacings'
        allocate(b%dxWeight(b%ni),b%dyWeight(b%nj))
        b%wall = [xb==0.5d0,xe==dble(nx)-0.5d0,yb==0.5d0,ye==dble(ny)-0.5d0]
        ! Eq. (20): h_d*(tau_d-1/2)=h_s*(tau_s-1/2)，包括非水动力矩。
        b%sn = 1.0d0/(0.5d0+(1.0d0/Snu-0.5d0)/b%h)
        b%sq = 1.0d0/(0.5d0+(1.0d0/Sq-0.5d0)/b%h)
        b%qk = 1.0d0/(0.5d0+(1.0d0/Qk-0.5d0)/b%h)
        b%qn = 1.0d0/(0.5d0+(1.0d0/Qnu-0.5d0)/b%h)
        b%gb = b%h*gBeta
        allocate(b%u(b%ni,b%nj),b%v(b%ni,b%nj),b%T(b%ni,b%nj),b%rho(b%ni,b%nj))
        allocate(b%Fx(b%ni,b%nj),b%Fy(b%ni,b%nj),b%Bx_prev(b%ni,b%nj),b%By_prev(b%ni,b%nj))
        allocate(b%f(b%ni,b%nj,0:8),b%g(b%ni,b%nj,0:4))
        allocate(b%f_post(0:b%ni+1,0:b%nj+1,0:8),b%g_post(0:b%ni+1,0:b%nj+1,0:4))
        allocate(b%p(b%ni,b%nj,packetSize,0:nh))
#ifdef steadyFlow
        allocate(b%up(b%ni,b%nj),b%vp(b%ni,b%nj),b%Tp(b%ni,b%nj))
#endif
    end subroutine make_block

    !===============================================================================================
    ! 积分权重：保证分区面积与线性坐标积分准确
    subroutine integration_weights(n,origin,h,lo,hi,w,first,last)
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: origin,h,lo,hi
        real(kind=8), intent(out) :: w(n)
        integer(kind=4), intent(out) :: first,last
        integer(kind=4) :: i
        real(kind=8) :: firstMoment,x,shift

        first = n+1
        last = 0
        firstMoment = 0.0d0
        do i = 1,n
            w(i) = max(0.0d0,min(hi,origin+dble(i)*h)-max(lo,origin+dble(i-1)*h))
            if (w(i)<=0.0d0) cycle
            first = min(first,i)
            last = i
            x = origin+(dble(i)-0.5d0)*h
            firstMoment = firstMoment+w(i)*x
        enddo
        if (abs(sum(w)-(hi-lo))>1.0d-10) error stop 'Integration weights do not cover owned interval'
        ! 中心粗块在共址节点之间积分，无需一阶矩修正；靠墙细块的裁剪仍可能不对称。
        ! 在末两个积分节点之间转移权重，补齐一阶矩，保持总长度。
        ! 两节点相距 h，转移 shift 后一阶矩增加 h*shift；无需修正的分区保持 shift=0。
        if (last>first) then
            shift = (0.5d0*(hi**2-lo**2)-firstMoment)/h
            w(last-1) = w(last-1)-shift
            w(last) = w(last)+shift
        endif
        if (any(w<0.0d0)) error stop 'Negative corrected integration weight'
        firstMoment = 0.0d0
        do i = first,last
            firstMoment = firstMoment+w(i)*(origin+(dble(i)-0.5d0)*h)
        enddo
        if (abs(firstMoment-0.5d0*(hi**2-lo**2))>1.0d-9*max(1.0d0,abs(firstMoment))) &
            error stop 'Integration weights do not integrate a linear coordinate exactly'
    end subroutine integration_weights

    !===============================================================================================
    ! 初始化各块的流场、温度场和分布函数
    subroutine initial_block(b)
        type(gridBlock), intent(inout) :: b
        integer(kind=4) :: i,j,a
        real(kind=8) :: x,y

        b%u = 0.0d0
        b%v = 0.0d0
        b%rho = 1.0d0
        b%Fx = 0.0d0
        b%Fy = 0.0d0
        b%Bx_prev = 0.0d0
        b%By_prev = 0.0d0
        b%f_post = 0.0d0
        b%g_post = 0.0d0
        b%p = 0.0d0
        do j = 1,b%nj
            y = (b%y0+(dble(j)-0.5d0)*b%h)/lengthUnit
            do i = 1,b%ni
                x = (b%x0+(dble(i)-0.5d0)*b%h)/lengthUnit
#ifdef SideHeatedCell
                b%T(i,j) = Thot+x/(dble(nx)/lengthUnit)*(Tcold-Thot)
#else
                b%T(i,j) = Thot+y/(dble(ny)/lengthUnit)*(Tcold-Thot)
                b%T(i,j) = b%T(i,j)+1.0d-3*(Thot-Tcold)*sin(2.0d0*pi*x/(dble(nx)/lengthUnit))* &
                    sin(pi*y/(dble(ny)/lengthUnit))
#endif
                do a = 0,8
                    b%f(i,j,a) = omega(a)
                enddo
                do a = 0,4
                    b%g(i,j,a) = omegaT(a)*b%T(i,j)
                enddo
            enddo
        enddo
#ifdef steadyFlow
        b%up = b%u
        b%vp = b%v
        b%Tp = b%T
#endif
    end subroutine initial_block

    !===============================================================================================
    ! OpenACC 数据区：传入数组并准备首次碰撞的接口状态
    subroutine enter_data_2d_openacc()
        integer(kind=4) :: b,l
        !$acc enter data copyin(ex,ey,omega,omegaT)
        do b = 1,nBlocks
            call block_device_data(blocks(b),.true.)
        enddo
        do l = 1,nLinks
            call link_device_data(links(l),.true.)
        enddo
        if (loadInitField==0) then
            do b = 1,nBlocks
                call pack_block(blocks(b),min(1,blocks(b)%nh))
            enddo
            if (nBlocks>1) then
                ! 首次碰撞前补齐两级接口。两次交换均读取同一份初始快照，不依赖块处理顺序。
                call exchange_interfaces(.true.,[0.0d0,1.0d0,0.0d0])
                call exchange_interfaces(.false.,[0.0d0,1.0d0,0.0d0])
                do b = 1,nBlocks
                    call pack_block(blocks(b),min(1,blocks(b)%nh))
                enddo
            endif
        endif
    end subroutine enter_data_2d_openacc

    !===============================================================================================
    subroutine block_device_data(b,entering)
        type(gridBlock), intent(inout) :: b
        logical, intent(in) :: entering
        ! 仅映射实际数组，不将含 allocatable 成员的宿主派生类型传给 GPU 内核。
        associate(f=>b%f,g=>b%g,fp=>b%f_post,gp=>b%g_post,u=>b%u,v=>b%v,T=>b%T,rho=>b%rho, &
            Fx=>b%Fx,Fy=>b%Fy,Bx=>b%Bx_prev,By=>b%By_prev,p=>b%p)
            if (entering) then
                !$acc enter data copyin(f,g,fp,gp,u,v,T,rho,Fx,Fy,Bx,By,p)
            else
                !$acc exit data delete(f,g,fp,gp,u,v,T,rho,Fx,Fy,Bx,By,p)
            endif
        end associate
    end subroutine block_device_data

    !===============================================================================================
    subroutine link_device_data(l,entering)
        type(blockLink), intent(inout) :: l
        logical, intent(in) :: entering
        associate(ti=>l%ti,tj=>l%tj,si=>l%si,sj=>l%sj,wx=>l%wx,wy=>l%wy,val=>l%values,same=>l%coincident)
            if (entering) then
                !$acc enter data copyin(ti,tj,si,sj,wx,wy,same) create(val)
            else
                !$acc exit data delete(ti,tj,si,sj,wx,wy,same,val)
            endif
        end associate
    end subroutine link_device_data

    !===============================================================================================
    subroutine update_host_block(b,full)
        type(gridBlock), intent(inout) :: b
        logical, intent(in) :: full
        associate(f=>b%f,g=>b%g,u=>b%u,v=>b%v,T=>b%T,rho=>b%rho,Fx=>b%Fx,Fy=>b%Fy, &
            Bx=>b%Bx_prev,By=>b%By_prev,p=>b%p)
            !$acc update self(u,v,T,rho) async(1)
            if (full) then
                !$acc update self(f,g,Fx,Fy,Bx,By,p) async(1)
            endif
        end associate
    end subroutine update_host_block

    !===============================================================================================
    subroutine update_host_all(full)
        logical, intent(in) :: full
        integer(kind=4) :: b

        do b = 1,nBlocks
            call update_host_block(blocks(b),full)
        enddo
        !$acc wait(1)
    end subroutine update_host_all

    !===============================================================================================
    subroutine exit_data_2d_openacc()
        integer(kind=4) :: b,l
        !$acc wait(1)
        do l = 1,nLinks
            call link_device_data(links(l),.false.)
        enddo
        do b = 1,nBlocks
            call block_device_data(blocks(b),.false.)
        enddo
        !$acc exit data delete(ex,ey,omega,omegaT)
    end subroutine exit_data_2d_openacc

    !===============================================================================================
    ! 块内推进：保持原流场与温度场的执行次序
    subroutine advance_block(b)
        type(gridBlock), intent(inout) :: b
        ! 原文件的执行次序保持不变。粗块的 dt 已吸收到松弛率及力增量中。
        ! 入口：人工边界已重建为本时刻的碰撞前状态；缓冲节点与内部节点一样参与碰撞。
        ! 出口：最外两层只作为待重建缓冲，不允许充当 donor；原算法的内部结果保留。
        call collision(b%ni,b%nj,b%f,b%f_post,b%rho,b%u,b%v,b%Fx,b%Fy,b%T,b%sn,b%sq,b%gb &
#ifdef SideHeatedHa
            ,b%h*B2sigemarho &
#endif
            )
        call streaming(b%ni,b%nj,b%f,b%f_post)
        call bounceback(b%ni,b%nj,b%f,b%f_post,b%wall(1),b%wall(2),b%wall(3),b%wall(4))
        call macro(b%ni,b%nj,b%f,b%rho,b%u,b%v,b%Fx,b%Fy)
        call collisionT(b%ni,b%nj,b%g,b%g_post,b%u,b%v,b%T,b%Bx_prev,b%By_prev,b%qk,b%qn)
        call streamingT(b%ni,b%nj,b%g,b%g_post)
        call bouncebackT(b%ni,b%nj,b%g,b%g_post,b%wall(1),b%wall(2),b%wall(3),b%wall(4))
        call macroT(b%ni,b%nj,b%g,b%T)
    end subroutine advance_block

    !===============================================================================================
    ! 多块推进：粗块预测、细块子步及两级同步
    subroutine advance_multiblock()
        integer(kind=4) :: b,k
        real(kind=8) :: theta,wt(0:2)

        if (nBlocks==1) then
            call advance_block(blocks(1))
            itc = itc+1
            return
        endif
        ! 1. 入口是粗细同步的 t 时刻。粗块缓冲已在初始交换或上次同步时补齐。
        ! 2. 粗块先预测 t+dt_c；p(:,:,:,0:2) 分别代表 t-dt_c、t、t+dt_c。
        call advance_block(blocks(1))
        call pack_block(blocks(1),2)
        do k = 1,refineRatio
            ! 3. 第 k 个细子步在 t+(k-1)*dt_f 开始：先补齐该时刻的碰撞前缓冲，再推进。
            ! 同级 donor 使用上一个子步结束的快照；粗级 donor 按三个粗时间层插值。
            theta = dble(k-1)/dble(refineRatio)
            call coarse_time_weights(theta,wt)
            call exchange_interfaces(.false.,wt)
            do b = 2,nBlocks
                call advance_block(blocks(b))
                call pack_block(blocks(b),0)
            enddo
        enddo
        ! 4. 两级均到 t+dt_c：细->粗修复粗缓冲，随后补齐细缓冲，供下次碰撞和同步输出使用。
        call exchange_interfaces(.true.,[0.0d0,0.0d0,1.0d0])
        call pack_block(blocks(1),2)
        call exchange_interfaces(.false.,[0.0d0,0.0d0,1.0d0])
        ! 5. 只在同步点滚动历史和整数时钟，重启文件必须保存完整历史。
        call rotate_coarse_history(blocks(1)%ni,blocks(1)%nj,blocks(1)%p)
        itc = itc+refineRatio
    end subroutine advance_multiblock

    !===============================================================================================
    subroutine coarse_time_weights(theta,wt)
        real(kind=8), intent(in) :: theta
        real(kind=8), intent(out) :: wt(0:2)
        ! theta 是当前细子步起点相对粗步起点的时间，不能使用子步终点时间。
        if (itc==0) then
            wt = [0.0d0,1.0d0-theta,theta]    ! 缺少负时间层时线性启动
        else
            wt = [0.5d0*theta*(theta-1.0d0),1.0d0-theta**2,0.5d0*theta*(theta+1.0d0)]
        endif
    end subroutine coarse_time_weights

    !===============================================================================================
    subroutine rotate_coarse_history(ni,nj,p)
        integer(kind=4), intent(in) :: ni,nj
        real(kind=8), intent(inout) :: p(ni,nj,packetSize,0:2)
        integer(kind=4) :: i,j,k
        !$acc parallel loop collapse(3) present(p) async(1)
        do k = 1,packetSize
            do j = 1,nj
                do i = 1,ni
                    p(i,j,k,0) = p(i,j,k,1)
                    p(i,j,k,1) = p(i,j,k,2)
                enddo
            enddo
        enddo
    end subroutine rotate_coarse_history

    !===============================================================================================
    logical function skin_node(b,i,j)
        type(gridBlock), intent(in) :: b
        integer(kind=4), intent(in) :: i,j

        skin_node = (.not.b%wall(1) .and. i<=interfaceSkin) .or. &
            (.not.b%wall(2) .and. i>b%ni-interfaceSkin) .or. &
            (.not.b%wall(3) .and. j<=interfaceSkin) .or. &
            (.not.b%wall(4) .and. j>b%nj-interfaceSkin)
    end function skin_node

    !===============================================================================================
    ! 接口来源：排除人工边界缓冲层，检查实际插值模板
    subroutine donor_stencil(receiver,x,y,donor,si,sj,wx,wy,coincident)
        integer(kind=4), intent(in) :: receiver
        real(kind=8), intent(in) :: x,y
        integer(kind=4), intent(out) :: donor,si,sj
        real(kind=8), intent(out) :: wx(4),wy(4)
        logical, intent(out) :: coincident
        integer(kind=4) :: d,il,ih,jl,jh,is,js
        real(kind=8) :: qx,qy,score,best

        best = -huge(1.0d0)
        donor = 0
        do d = 1,nBlocks
            if (d==receiver) cycle
            ! 来源模板严格避开刚推进后尚未修复的两层人工边界。
            il = 1
            ih = blocks(d)%ni
            jl = 1
            jh = blocks(d)%nj
            if (.not.blocks(d)%wall(1)) il = 1+interfaceSkin
            if (.not.blocks(d)%wall(2)) ih = ih-interfaceSkin
            if (.not.blocks(d)%wall(3)) jl = 1+interfaceSkin
            if (.not.blocks(d)%wall(4)) jh = jh-interfaceSkin
            qx = (x-blocks(d)%x0)/blocks(d)%h+0.5d0
            qy = (y-blocks(d)%y0)/blocks(d)%h+0.5d0
            if (qx<dble(il) .or. qx>dble(ih) .or. qy<dble(jl) .or. qy>dble(jh)) cycle
            if (ih-il<3 .or. jh-jl<3) cycle
            is = max(il,min(floor(qx)-1,ih-3))
            js = max(jl,min(floor(qy)-1,jh-3))
            score = min(qx-il,ih-qx,qy-jl,jh-qy)*blocks(d)%h
            ! 可用时优先同级直接交换；粗块始终从细块获取边界。
            if (blocks(d)%h==blocks(receiver)%h) score = score+1.0d6
            if (score<=best) cycle
            best = score
            donor = d
            si = is
            sj = js
            coincident = abs(qx-dble(nint(qx)))<1.0d-12 .and. abs(qy-dble(nint(qy)))<1.0d-12
            if (coincident) then
                si = nint(qx)
                sj = nint(qy)
                wx = [1.0d0,0.0d0,0.0d0,0.0d0]
                wy = wx
            else
                call lagrange_weights(qx-dble(is),wx)
                call lagrange_weights(qy-dble(js),wy)
            endif
        enddo
        if (donor==0) then
            write(*,*) 'No interior four-point donor stencil:',receiver,x,y
            error stop 'Invalid overlap geometry'
        endif
    end subroutine donor_stencil

    !===============================================================================================
    subroutine lagrange_weights(q,w)
        real(kind=8), intent(in) :: q
        real(kind=8), intent(out) :: w(4)
        integer(kind=4) :: a,k

        w = 1.0d0
        do a = 0,3
            do k = 0,3
                if (a/=k) w(a+1) = w(a+1)*(q-dble(k))/dble(a-k)
            enddo
        enddo
    end subroutine lagrange_weights

    !===============================================================================================
    subroutine section_weights(q,n,first,w,dw)
        real(kind=8), intent(in) :: q
        integer(kind=4), intent(in) :: n
        integer(kind=4), intent(out) :: first
        real(kind=8), intent(out) :: w(4),dw(4)
        real(kind=8) :: z,term
        integer(kind=4) :: a,k,l

        first = max(1,min(floor(q)-1,n-3))
        z = q-dble(first)
        call lagrange_weights(z,w)
        dw = 0.0d0
        do a = 0,3
            do k = 0,3
                if (k==a) cycle
                term = 1.0d0/dble(a-k)
                do l = 0,3
                    if (l/=a .and. l/=k) term = term*(z-dble(l))/dble(a-l)
                enddo
                dw(a+1) = dw(a+1)+term
            enddo
        enddo
    end subroutine section_weights

    !===============================================================================================
    ! 建立各块之间的接收节点、来源节点和空间权重
    subroutine build_links()
        integer(kind=4) :: b,d,i,j,l,n,si,sj,counts(maxBlocks,maxBlocks),idx(maxBlocks,maxBlocks)
        real(kind=8) :: x,y,wx(4),wy(4)
        logical :: coincident

        counts = 0
        idx = 0
        nLinks = 0
        do b = 1,nBlocks
            do j = 1,blocks(b)%nj
                do i = 1,blocks(b)%ni
                    if (.not.skin_node(blocks(b),i,j)) cycle
                    x = blocks(b)%x0+(dble(i)-0.5d0)*blocks(b)%h
                    y = blocks(b)%y0+(dble(j)-0.5d0)*blocks(b)%h
                    call donor_stencil(b,x,y,d,si,sj,wx,wy,coincident)
                    counts(b,d) = counts(b,d)+1
                enddo
            enddo
        enddo
        do b = 1,nBlocks
            do d = 1,nBlocks
                if (counts(b,d)==0) cycle
                nLinks = nLinks+1
                l = nLinks
                idx(b,d) = l
                n = counts(b,d)
                links(l)%receiver = b
                links(l)%donor = d
                links(l)%count = n
                allocate(links(l)%ti(n),links(l)%tj(n),links(l)%si(n),links(l)%sj(n))
                allocate(links(l)%coincident(n))
                allocate(links(l)%wx(4,n),links(l)%wy(4,n),links(l)%values(packetSize,n))
            enddo
        enddo
        counts = 0
        do b = 1,nBlocks
            do j = 1,blocks(b)%nj
                do i = 1,blocks(b)%ni
                    if (.not.skin_node(blocks(b),i,j)) cycle
                    x = blocks(b)%x0+(dble(i)-0.5d0)*blocks(b)%h
                    y = blocks(b)%y0+(dble(j)-0.5d0)*blocks(b)%h
                    call donor_stencil(b,x,y,d,si,sj,wx,wy,coincident)
                    counts(b,d) = counts(b,d)+1
                    n = counts(b,d)
                    l = idx(b,d)
                    links(l)%ti(n) = i
                    links(l)%tj(n) = j
                    links(l)%si(n) = si
                    links(l)%sj(n) = sj
                    links(l)%wx(:,n) = wx
                    links(l)%wy(:,n) = wy
                    links(l)%coincident(n) = coincident
                    if (blocks(b)%h>blocks(d)%h .and. .not.coincident) &
                        error stop 'Coarse interface node is not aligned with a fine node'
                enddo
            enddo
        enddo
    end subroutine build_links

    !===============================================================================================
    subroutine equilibrium_moments(rh,ux,uy,temp,meq,neq)
        !$acc routine seq
        real(kind=8), intent(in) :: rh,ux,uy,temp
        real(kind=8), intent(out) :: meq(0:8),neq(0:4)

        meq(0) = rh
        meq(1) = rh*(-2.0d0+3.0d0*(ux*ux+uy*uy))
        meq(2) = rh*(1.0d0-3.0d0*(ux*ux+uy*uy))
        meq(3) = rh*ux
        meq(4) = -rh*ux
        meq(5) = rh*uy
        meq(6) = -rh*uy
        meq(7) = rh*(ux*ux-uy*uy)
        meq(8) = rh*ux*uy
        neq = [temp,temp*ux,temp*uy,thermalA*temp,0.0d0]
    end subroutine equilibrium_moments

    !===============================================================================================
    subroutine force_moments(ux,uy,fx,fy,fm)
        !$acc routine seq
        real(kind=8), intent(in) :: ux,uy,fx,fy
        real(kind=8), intent(out) :: fm(0:8)

        fm = [0.0d0,6.0d0*(ux*fx+uy*fy),-6.0d0*(ux*fx+uy*fy),fx,-fx,fy,-fy, &
            2.0d0*(ux*fx-uy*fy),ux*fy+uy*fx]
    end subroutine force_moments

    !===============================================================================================
    subroutine pack_block(b,slot)
        type(gridBlock), intent(inout) :: b
        integer(kind=4), intent(in) :: slot

        call encode_packets(b%ni,b%nj,b%nh,slot,b%h,b%sn,b%sq,b%qk,b%qn, &
            b%f,b%g,b%rho,b%u,b%v,b%T,b%Fx,b%Fy,b%Bx_prev,b%By_prev,b%p)
    end subroutine pack_block

    !===============================================================================================
    ! 将宏观量及归一化非平衡矩编码为接口交换量
    subroutine encode_packets(ni,nj,nh,slot,h,sn,sq,qk,qn,f,g,rho,u,v,T,Fx,Fy,Bx_prev,By_prev,p)
        integer(kind=4), intent(in) :: ni,nj,nh,slot
        real(kind=8), intent(in) :: h,sn,sq,qk,qn,f(ni,nj,0:8),g(ni,nj,0:4)
        real(kind=8), intent(in) :: rho(ni,nj),u(ni,nj),v(ni,nj),T(ni,nj),Fx(ni,nj),Fy(ni,nj)
        real(kind=8), intent(in) :: Bx_prev(ni,nj),By_prev(ni,nj)
        real(kind=8), intent(inout) :: p(ni,nj,packetSize,0:nh)
        integer(kind=4) :: i,j,a
        real(kind=8) :: m(0:8),meq(0:8),fm(0:8),n(0:4),neq(0:4),s(0:8),q(0:4),src(0:4),fv(0:8),gv(0:4)
        !$acc parallel loop collapse(2) present(f,g,rho,u,v,T,Fx,Fy,Bx_prev,By_prev,p) async(1) &
        !$acc& private(a,m,meq,fm,n,neq,s,q,src,fv,gv)
        do j = 1,nj
            do i = 1,ni
                do a = 0,8
                    fv(a) = f(i,j,a)
                enddo
                do a = 0,4
                    gv(a) = g(i,j,a)
                enddo
                call flow_moments(fv,m)
                call thermal_moments(gv,n)
                call equilibrium_moments(rho(i,j),u(i,j),v(i,j),T(i,j),meq,neq)
                call force_moments(u(i,j),v(i,j),Fx(i,j),Fy(i,j),fm)
                s = [0.0d0,sn,sn,0.0d0,sq,0.0d0,sq,sn,sn]
                q = [0.0d0,qk,qk,qn,qn]
                src = 0.0d0
#ifdef EnableUseG
                ! 原 D2Q5 离散修正源。以当前 B 与上次 collisionT 保存的 B_prev 形成时间差。
                src(1) = u(i,j)*T(i,j)-Bx_prev(i,j)
                src(2) = v(i,j)*T(i,j)-By_prev(i,j)
#endif
                p(i,j,1,slot) = rho(i,j)
                p(i,j,2,slot) = u(i,j)
                p(i,j,3,slot) = v(i,j)
                p(i,j,4,slot) = T(i,j)
                p(i,j,5,slot) = Fx(i,j)/h
                p(i,j,6,slot) = Fy(i,j)/h
                do a = 0,8
                    p(i,j,7+a,slot) = s(a)/h*(m(a)-meq(a)+0.5d0*fm(a))
                enddo
                do a = 0,4
                    p(i,j,16+a,slot) = q(a)/h*(n(a)-neq(a)+0.5d0*src(a))
                enddo
                p(i,j,21,slot) = src(1)/h
                p(i,j,22,slot) = src(2)/h
            enddo
        enddo
    end subroutine encode_packets

    !===============================================================================================
    ! 接口交换：从来源块取值，再重建接收块的碰撞前状态
    subroutine exchange_interfaces(coarse_receiver,wt)
        logical, intent(in) :: coarse_receiver
        real(kind=8), intent(in) :: wt(0:2)
        integer(kind=4) :: l,b,d

        do l = 1,nLinks
            b = links(l)%receiver
            d = links(l)%donor
            if ((b==1) .neqv. coarse_receiver) cycle
            call interpolate_packets(blocks(d)%ni,blocks(d)%nj,blocks(d)%nh,blocks(d)%p,links(l)%count, &
                links(l)%si,links(l)%sj,links(l)%wx,links(l)%wy,links(l)%coincident,wt,links(l)%values)
        enddo
        do l = 1,nLinks
            b = links(l)%receiver
            if ((b==1) .neqv. coarse_receiver) cycle
            call apply_packets(blocks(b)%ni,blocks(b)%nj,blocks(b)%h,blocks(b)%sn,blocks(b)%sq,blocks(b)%qk,blocks(b)%qn, &
                blocks(b)%f,blocks(b)%g,blocks(b)%rho,blocks(b)%u,blocks(b)%v,blocks(b)%T, &
                blocks(b)%Fx,blocks(b)%Fy,blocks(b)%Bx_prev,blocks(b)%By_prev, &
                links(l)%count,links(l)%ti,links(l)%tj,links(l)%values)
        enddo
    end subroutine exchange_interfaces

    !===============================================================================================
    subroutine interpolate_packets(ni,nj,nh,p,count,si,sj,wx,wy,coincident,wt,val)
        integer(kind=4), intent(in) :: ni,nj,nh,count,si(count),sj(count)
        real(kind=8), intent(in) :: p(ni,nj,packetSize,0:nh),wx(4,count),wy(4,count),wt(0:2)
        logical, intent(in) :: coincident(count)
        real(kind=8), intent(out) :: val(packetSize,count)
        integer(kind=4) :: c,a,ix,iy,k
        real(kind=8) :: value,wk
        !$acc parallel loop collapse(2) present(p,si,sj,wx,wy,coincident,val) firstprivate(wt) async(1) private(ix,iy,k,value,wk)
        do c = 1,count
            do a = 1,packetSize
                value = 0.0d0
                do k = 0,nh
                    wk = 1.0d0
                    if (nh==2) wk = wt(k)
                    if (coincident(c)) then
                        value = value+wk*p(si(c),sj(c),a,k)
                    else
                        do iy = 1,4
                            do ix = 1,4
                                value = value+wk*wx(ix,c)*wy(iy,c)*p(si(c)+ix-1,sj(c)+iy-1,a,k)
                            enddo
                        enddo
                    endif
                enddo
                val(a,c) = value
            enddo
        enddo
    end subroutine interpolate_packets

    !===============================================================================================
    subroutine apply_packets(ni,nj,h,sn,sq,qk,qn,f,g,rho,u,v,T,Fx,Fy,Bx_prev,By_prev,count,ti,tj,val)
        integer(kind=4), intent(in) :: ni,nj,count,ti(count),tj(count)
        real(kind=8), intent(in) :: h,sn,sq,qk,qn,val(packetSize,count)
        real(kind=8), intent(inout) :: f(ni,nj,0:8),g(ni,nj,0:4),rho(ni,nj),u(ni,nj),v(ni,nj),T(ni,nj)
        real(kind=8), intent(inout) :: Fx(ni,nj),Fy(ni,nj),Bx_prev(ni,nj),By_prev(ni,nj)
        integer(kind=4) :: c,i,j,a
        real(kind=8) :: m(0:8),meq(0:8),fm(0:8),n(0:4),neq(0:4),s(0:8),q(0:4),fv(0:8),gv(0:4),src(0:4)
        !$acc parallel loop present(f,g,rho,u,v,T,Fx,Fy,Bx_prev,By_prev,ti,tj,val) async(1) &
        !$acc& private(i,j,a,m,meq,fm,n,neq,s,q,fv,gv,src)
        do c = 1,count
            i = ti(c)
            j = tj(c)
            rho(i,j) = val(1,c)
            u(i,j) = val(2,c)
            v(i,j) = val(3,c)
            T(i,j) = val(4,c)
            Fx(i,j) = h*val(5,c)
            Fy(i,j) = h*val(6,c)
            call equilibrium_moments(rho(i,j),u(i,j),v(i,j),T(i,j),meq,neq)
            call force_moments(u(i,j),v(i,j),Fx(i,j),Fy(i,j),fm)
            s = [0.0d0,sn,sn,0.0d0,sq,0.0d0,sq,sn,sn]
            q = [0.0d0,qk,qk,qn,qn]
            src = 0.0d0
            src(1) = h*val(21,c)
            src(2) = h*val(22,c)
            m = meq
            n = neq
            do a = 0,8
                if (s(a)>0.0d0) m(a) = meq(a)+h/s(a)*val(7+a,c)-0.5d0*fm(a)
            enddo
            ! Eq. (21)：守恒矩直接用宏观量和半步力重建，不除以零松弛率。
            m(0) = rho(i,j)
            m(3) = rho(i,j)*u(i,j)-0.5d0*Fx(i,j)
            m(5) = rho(i,j)*v(i,j)-0.5d0*Fy(i,j)
            do a = 1,4
                n(a) = neq(a)+h/q(a)*val(16+a,c)-0.5d0*src(a)
            enddo
            n(0) = T(i,j)
            call flow_populations(m,fv)
            call thermal_populations(n,gv)
            do a = 0,8
                f(i,j,a) = fv(a)
            enddo
            do a = 0,4
                g(i,j,a) = gv(a)
            enddo
#ifdef EnableUseG
            Bx_prev(i,j) = u(i,j)*T(i,j)-src(1)
            By_prev(i,j) = v(i,j)*T(i,j)-src(2)
#else
            Bx_prev(i,j) = 0.0d0
            By_prev(i,j) = 0.0d0
#endif
        enddo
    end subroutine apply_packets

    !===============================================================================================
    ! 以下八个块内子程序取自原文件。修改限于：显式块数组/局部尺寸/局部松弛率参数、物理壁面标记。
    !===============================================================================================

    !===============================================================================================
    subroutine collision(nx,ny,f,f_post,rho,u,v,Fx,Fy,T,Snu,Sq,gBeta &
#ifdef SideHeatedHa
        ,B2sigemarho &
#endif
        )

    implicit none
    integer(kind=4), intent(in) :: nx,ny
    real(kind=8), intent(inout) :: f(nx,ny,0:8),f_post(0:nx+1,0:ny+1,0:8)
    real(kind=8), intent(inout) :: u(nx,ny),v(nx,ny),T(nx,ny),rho(nx,ny),Fx(nx,ny),Fy(nx,ny)
    real(kind=8), intent(in) :: Snu,Sq,gBeta
#ifdef SideHeatedHa
    real(kind=8), intent(in) :: B2sigemarho
#endif

    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: m(0:8), m_post(0:8), meq(0:8)
    real(kind=8) :: s(0:8)
    real(kind=8) :: fSource(0:8)

    !$acc parallel loop gang vector collapse(2) present(f,f_post,rho,u,v,Fx,Fy,T) async(1) &
    !$acc& private(alpha,s,m,m_post,meq,fSource)
    do j = 1, ny
        do i = 1, nx

            m(0) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
            m(1) = -4.0d0*f(i,j,0)-f(i,j,1)-f(i,j,2)-f(i,j,3)-f(i,j,4)+2.0d0*(f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8))
            m(2) = 4.0d0*f(i,j,0)-2.0d0*(f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4))+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
            m(3) = f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)
            m(4) = -2.0d0*f(i,j,1)+2.0d0*f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)
            m(5) = f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)
            m(6) = -2.0d0*f(i,j,2)+2.0d0*f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)
            m(7) = f(i,j,1)-f(i,j,2)+f(i,j,3)-f(i,j,4)
            m(8) = f(i,j,5)-f(i,j,6)+f(i,j,7)-f(i,j,8)

            meq(0) = rho(i,j)
            meq(1) = rho(i,j)*( -2.0d0+3.0d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j)) )
            meq(2) = rho(i,j)*( 1.0d0-3.0d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j)) )
            meq(3) = rho(i,j)*u(i,j)
            meq(4) = -rho(i,j)*u(i,j)
            meq(5) = rho(i,j)*v(i,j)
            meq(6) = -rho(i,j)*v(i,j)
            meq(7) = rho(i,j)*( u(i,j)*u(i,j)-v(i,j)*v(i,j) )
            meq(8) = rho(i,j)*( u(i,j)*v(i,j) )

            s(0) = 0.0d0    !!s_{\rho}
            s(1) = Snu    !!s_{e}
            s(2) = Snu    !!s_{\epsilon}
            s(3) = 0.0d0    !!s_{j}
            s(4) = Sq    !!s_{q}
            s(5) = 0.0d0    !!s_{j}
            s(6) = Sq    !!s_{q}
            s(7) = Snu    !!s_{\nu}
            s(8) = Snu    !!s_{\nu}

            Fx(i,j) = 0.0d0
            Fy(i,j) = rho(i,j)*gBeta*(T(i,j)-Tref)    !动量方程上的源项，即浮力项

#ifdef    SideHeatedHa
            Fx(i,j) = 0.0d0+B2sigemarho*(v(i,j)*sin(phi)*cos(phi)-u(i,j)*sin(phi)*sin(phi))
            Fy(i,j) = rho(i,j)*gBeta*(T(i,j)-Tref)+ rho(i,j)*B2sigemarho*(u(i,j)*sin(phi)*cos(phi)&
                -v(i,j)*cos(phi)*cos(phi))    !动量方程上的源项，即浮力项加磁场
#endif

            fSource(0) = 0.0d0    !将源项F对应的贡献投影到各个矩中，并做半步修正
            fSource(1) = (6.0d0-3.0d0*s(1))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
            fSource(2) = -(6.0d0-3.0d0*s(2))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
            fSource(3) = (1.0d0-0.5d0*s(3))*Fx(i,j)
            fSource(4) = -(1.0d0-0.5d0*s(4))*Fx(i,j)
            fSource(5) = (1.0d0-0.5d0*s(5))*Fy(i,j)
            fSource(6) = -(1.0d0-0.5d0*s(6))*Fy(i,j)
            fSource(7) = (2.0d0-s(7))*(u(i,j)*Fx(i,j)-v(i,j)*Fy(i,j))
            fSource(8) = (1.0d0-0.5d0*s(8))*(u(i,j)*Fy(i,j)+v(i,j)*Fx(i,j))    !这边是乘以M变到矩空间，然后再乘以1-1/2S修正

            do alpha = 0, 8
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+fSource(alpha)    !矩空间碰撞
            enddo

            f_post(i,j,0) = m_post(0)/9.0d0-m_post(1)/9.0d0+m_post(2)/9.0d0    !这边是乘以M逆
            f_post(i,j,1) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+m_post(3)/6.0d0-m_post(4)/6.0d0 &
                +m_post(7)/4.0d0
            f_post(i,j,2) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                +m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)/4.0d0
            f_post(i,j,3) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0-m_post(3)/6.0d0+m_post(4)/6.0d0 &
                +m_post(7)/4.0d0
            f_post(i,j,4) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                -m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)/4.0d0
            f_post(i,j,5) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                +m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)/4.0d0
            f_post(i,j,6) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                +m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)/4.0d0
            f_post(i,j,7) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                -m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)/4.0d0
            f_post(i,j,8) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                -m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)/4.0d0

        enddo
    enddo
    return
    end subroutine collision

    !===============================================================================================
    subroutine streaming(nx,ny,f,f_post)    !先迁移，再边界处理

        implicit none
        integer(kind=4), intent(in) :: nx,ny
        real(kind=8), intent(inout) :: f(nx,ny,0:8),f_post(0:nx+1,0:ny+1,0:8)

        integer(kind=4) :: i, j
        integer(kind=4) :: ip, jp
        integer(kind=4) :: alpha

        !$acc parallel loop gang vector collapse(2) present(f,f_post,ex,ey) async(1) private(alpha,ip,jp)
        do j = 1, ny
            do i = 1, nx
                do alpha = 0, 8    !上游格点索引：fα(i,j) <- f_postα(i-exα, j-eyα)
                    ip = i-ex(alpha)    !边界附近 (ip/jp 可能为 0 或 nx+1/ny+1)，需在 bounceback/周期边界处理中覆盖修正边界分布
                    jp = j-ey(alpha)    !ghost 层在初始化中为 0，保证不会出现未初始化垃圾值

                    f(i,j,alpha) = f_post(ip,jp,alpha)
                enddo
            enddo
        enddo
        return
    end subroutine streaming

    !===============================================================================================
    subroutine bounceback(nx,ny,f,f_post,leftWall,rightWall,bottomWall,topWall)

        implicit none
        integer(kind=4), intent(in) :: nx,ny
        real(kind=8), intent(inout) :: f(nx,ny,0:8),f_post(0:nx+1,0:ny+1,0:8)
        logical, intent(in) :: leftWall,rightWall,bottomWall,topWall

        integer(kind=4) :: i, j
        ! integer(kind=4) :: alpha

#ifdef VerticalWallsPeriodicalU
        !$acc parallel loop gang vector present(f,f_post) async(1)
        do j = 1, ny    !速度边界垂直边界周期，直接方向相同，跨边界的入射分布
            !Left side (i=1)
            if (leftWall) f(1,j,1) = f_post(nx,j,1)
            if (leftWall) f(1,j,5) = f_post(nx,j,5)
            if (leftWall) f(1,j,8) = f_post(nx,j,8)

            !Right side (i=nx)
            if (rightWall) f(nx,j,3) = f_post(1,j,3)
            if (rightWall) f(nx,j,6) = f_post(1,j,6)
            if (rightWall) f(nx,j,7) = f_post(1,j,7)
        enddo
#endif

#ifdef VerticalWallsNoslip
        !$acc parallel loop gang vector present(f,f_post) async(1)
        do j = 1, ny    !速度边界垂直边界静止壁无滑移，直接反弹，方向相反
            !Left side (i=1)
            if (leftWall) f(1,j,1) = f_post(1,j,3)
            if (leftWall) f(1,j,5) = f_post(1,j,7)
            if (leftWall) f(1,j,8) = f_post(1,j,6)

            !Right side (i=nx)
            if (rightWall) f(nx,j,3) = f_post(nx,j,1)
            if (rightWall) f(nx,j,6) = f_post(nx,j,8)
            if (rightWall) f(nx,j,7) = f_post(nx,j,5)
        enddo
#endif

#ifdef HorizontalWallsNoslip
        !$acc parallel loop gang vector present(f,f_post) async(1)
        do i = 1, nx    !速度边界水平边界无滑移，直接反弹，方向相反
            !Bottom side (j=1)
            if (bottomWall) f(i,1,2) = f_post(i,1,4)
            if (bottomWall) f(i,1,5) = f_post(i,1,7)
            if (bottomWall) f(i,1,6) = f_post(i,1,8)

            !Top side (j=ny)
            if (topWall) f(i,ny,4) = f_post(i,ny,2)
            if (topWall) f(i,ny,7) = f_post(i,ny,5)
            if (topWall) f(i,ny,8) = f_post(i,ny,6)
        enddo
#endif

        return
    end subroutine bounceback

    !===============================================================================================
    subroutine macro(nx,ny,f,rho,u,v,Fx,Fy)

        implicit none
        integer(kind=4), intent(in) :: nx,ny
        real(kind=8), intent(in) :: f(nx,ny,0:8),Fx(nx,ny),Fy(nx,ny)
        real(kind=8), intent(inout) :: rho(nx,ny),u(nx,ny),v(nx,ny)

        integer(kind=4) :: i, j

        !$acc parallel loop gang vector collapse(2) present(f,rho,u,v,Fx,Fy) async(1)
        do j = 1, ny
            do i = 1, nx
                rho(i,j) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
                u(i,j) = ( f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)+0.5d0*Fx(i,j) )/rho(i,j)    !含力LBM的半步动量修正：rho*u = Σ f e + 0.5*F，对应Guo forcing的二阶定义
                v(i,j) = ( f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)+0.5d0*Fy(i,j) )/rho(i,j)
            enddo
        enddo
        return
    end subroutine macro

    !===============================================================================================
    subroutine collisionT(nx,ny,g,g_post,u,v,T,Bx_prev,By_prev,Qk,Qnu)

        implicit none
        integer(kind=4), intent(in) :: nx,ny
        real(kind=8), intent(inout) :: g(nx,ny,0:4),g_post(0:nx+1,0:ny+1,0:4)
        real(kind=8), intent(in) :: u(nx,ny),v(nx,ny),T(nx,ny),Qk,Qnu
        real(kind=8), intent(inout) :: Bx_prev(nx,ny),By_prev(nx,ny)

        integer(kind=4) :: i, j
        integer(kind=4) :: alpha
        real(kind=8) :: n(0:4), n_post(0:4), neq(0:4)
        real(kind=8) :: q(0:4)
        real(kind=8) :: Bx, By
        real(kind=8) :: dBx, dBy
        real(kind=8) :: SG

        SG = 1.0d0 - 0.5d0*Qk
        !$acc parallel loop gang vector collapse(2) present(g,g_post,u,v,T,Bx_prev,By_prev) async(1) &
        !$acc& private(alpha,n,neq,q,n_post,Bx,By,dBx,dBy)
        do j = 1, ny
            do i = 1, nx

                Bx = u(i,j) * T(i,j)
                By = v(i,j) * T(i,j)

#ifdef EnableUseG
                dBx = Bx - Bx_prev(i,j)
                dBy = By - By_prev(i,j)
#else
                dBx = 0.0d0
                dBy = 0.0d0
#endif

#ifdef EnableUseG
                Bx_prev(i,j) = Bx
                By_prev(i,j) = By
#endif

                n(0) = g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
                n(1) = g(i,j,1)-g(i,j,3)
                n(2) = g(i,j,2)-g(i,j,4)
                n(3) = -4.0d0*g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
                n(4) = g(i,j,1)-g(i,j,2)+g(i,j,3)-g(i,j,4)

                neq(0) = T(i,j)
                neq(1) = T(i,j)*u(i,j)
                neq(2) = T(i,j)*v(i,j)
#ifdef EnableLegacyThermalScheme
                neq(3) = T(i,j)*paraA
#else
                neq(3) = T(i,j)*(-2.0d0/3.0d0)
#endif
                neq(4) = 0.0d0

                q(0) = 0.0d0
                q(1) = Qk
                q(2) = Qk
                q(3) = Qnu
                q(4) = Qnu

                n_post(0) = n(0)-q(0)*(n(0)-neq(0))
                n_post(1) = n(1)-q(1)*(n(1)-neq(1))+ SG*dBx
                n_post(2) = n(2)-q(2)*(n(2)-neq(2))+ SG*dBy
                n_post(3) = n(3)-q(3)*(n(3)-neq(3))
                n_post(4) = n(4)-q(4)*(n(4)-neq(4))

                g_post(i,j,0) = 0.2d0*n_post(0)-0.2d0*n_post(3)
                g_post(i,j,1) = 0.2d0*n_post(0)+0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
                g_post(i,j,2) = 0.2d0*n_post(0)+0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4)
                g_post(i,j,3) = 0.2d0*n_post(0)-0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
                g_post(i,j,4) = 0.2d0*n_post(0)-0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4)
            enddo
        enddo
        return
    end subroutine collisionT

    !===============================================================================================
    subroutine streamingT(nx,ny,g,g_post)

        implicit none
        integer(kind=4), intent(in) :: nx,ny
        real(kind=8), intent(inout) :: g(nx,ny,0:4),g_post(0:nx+1,0:ny+1,0:4)

        integer(kind=4) :: i, j
        integer(kind=4) :: ip, jp
        integer(kind=4) :: alpha

        !$acc parallel loop gang vector collapse(2) present(g,g_post,ex,ey) async(1) private(alpha,ip,jp)
        do j = 1, ny
            do i = 1, nx
                do alpha = 0, 4
                    ip = i-ex(alpha)
                    jp = j-ey(alpha)

                    g(i,j,alpha) = g_post(ip,jp,alpha)
                enddo
            enddo
        enddo
        return
    end subroutine streamingT

    !===============================================================================================
    subroutine bouncebackT(nx,ny,g,g_post,leftWall,rightWall,bottomWall,topWall)

        implicit none
        integer(kind=4), intent(in) :: nx,ny
        real(kind=8), intent(inout) :: g(nx,ny,0:4),g_post(0:nx+1,0:ny+1,0:4)
        logical, intent(in) :: leftWall,rightWall,bottomWall,topWall

        integer(kind=4) :: i, j
        !integer(kind=4) :: alpha

#ifdef VerticalWallsPeriodicalT
        !$acc parallel loop gang vector present(g,g_post) async(1)
        do j = 1, ny
            !Left boundary
            if (leftWall) g(1,j,1) = g_post(nx,j,1)

            !Right boundary
            if (rightWall) g(nx,j,3) = g_post(1,j,3)
        enddo
#endif

#ifdef VerticalWallsConstT
        !$acc parallel loop gang vector present(g,g_post,omegaT) async(1)
        do j = 1, ny
            !Left boundary
#ifdef EnableLegacyThermalScheme
            if (leftWall) g(1,j,1) = -g_post(1,j,3)+(4.0d0+paraA)/10.0d0*Thot
#else
            if (leftWall) g(1,j,1) = -g_post(1,j,3)+2.0d0*omegaT(3)*Thot
#endif
            !Right boundary
#ifdef EnableLegacyThermalScheme
            if (rightWall) g(nx,j,3) = -g_post(nx,j,1)+(4.0d0+paraA)/10.0d0*Tcold
#else
            if (rightWall) g(nx,j,3) = -g_post(nx,j,1)+2.0d0*omegaT(1)*Tcold
#endif
        enddo
#endif

#ifdef VerticalWallsAdiabatic
        !$acc parallel loop gang vector present(g,g_post) async(1)
        do j = 1, ny
            !Left boundary
            if (leftWall) g(1,j,1) = g_post(1,j,3)

            !Right boundary
            if (rightWall) g(nx,j,3) = g_post(nx,j,1)
        enddo
#endif

#ifdef HorizontalWallsAdiabatic
        !$acc parallel loop gang vector present(g,g_post) async(1)
        do i = 1, nx
            !Bottom side
            if (bottomWall) g(i,1,2) = g_post(i,1,4)

            !Top side
            if (topWall) g(i,ny,4) = g_post(i,ny,2)
        enddo
#endif

#ifdef HorizontalWallsConstT
        !$acc parallel loop gang vector present(g,g_post,omegaT) async(1)
        do i = 1, nx
            !Bottom side
#ifdef EnableLegacyThermalScheme
            if (bottomWall) g(i,1,2) = -g_post(i,1,4)+(4.0d0+paraA)/10.0d0*Thot
#else
            if (bottomWall) g(i,1,2) = -g_post(i,1,4)+2.0d0*omegaT(4)*Thot
#endif
            !Top side
#ifdef EnableLegacyThermalScheme
            if (topWall) g(i,ny,4) = -g_post(i,ny,2)+(4.0d0+paraA)/10.0d0*Tcold
#else
            if (topWall) g(i,ny,4) = -g_post(i,ny,2)+2.0d0*omegaT(2)*Tcold
#endif
        enddo
#endif

        return
    end subroutine bouncebackT

    !===============================================================================================
    subroutine macroT(nx,ny,g,T)

        implicit none
        integer(kind=4), intent(in) :: nx,ny
        real(kind=8), intent(in) :: g(nx,ny,0:4)
        real(kind=8), intent(inout) :: T(nx,ny)

        integer(kind=4) :: i, j

        !$acc parallel loop gang vector collapse(2) present(g,T) async(1)
        do j = 1, ny
            do i = 1, nx
                T(i,j) = g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
            enddo
        enddo
        return
    end subroutine macroT

    !===============================================================================================
    subroutine flow_moments(fv,m)
        !$acc routine seq
        real(kind=8), intent(in) :: fv(0:8)
        real(kind=8), intent(out) :: m(0:8)

        m(0) = fv(0)+fv(1)+fv(2)+fv(3)+fv(4)+fv(5)+fv(6)+fv(7)+fv(8)
        m(1) = -4.0d0*fv(0)-fv(1)-fv(2)-fv(3)-fv(4)+2.0d0*(fv(5)+fv(6)+fv(7)+fv(8))
        m(2) = 4.0d0*fv(0)-2.0d0*(fv(1)+fv(2)+fv(3)+fv(4))+fv(5)+fv(6)+fv(7)+fv(8)
        m(3) = fv(1)-fv(3)+fv(5)-fv(6)-fv(7)+fv(8)
        m(4) = -2.0d0*fv(1)+2.0d0*fv(3)+fv(5)-fv(6)-fv(7)+fv(8)
        m(5) = fv(2)-fv(4)+fv(5)+fv(6)-fv(7)-fv(8)
        m(6) = -2.0d0*fv(2)+2.0d0*fv(4)+fv(5)+fv(6)-fv(7)-fv(8)
        m(7) = fv(1)-fv(2)+fv(3)-fv(4)
        m(8) = fv(5)-fv(6)+fv(7)-fv(8)

    end subroutine flow_moments

    !===============================================================================================
    subroutine flow_populations(m,fv)
        !$acc routine seq
        real(kind=8), intent(in) :: m(0:8)
        real(kind=8), intent(out) :: fv(0:8)

        fv(0) = m(0)/9.0d0-m(1)/9.0d0+m(2)/9.0d0    !这边是乘以M逆
        fv(1) = m(0)/9.0d0-m(1)/36.0d0-m(2)/18.0d0+m(3)/6.0d0-m(4)/6.0d0 &
            +m(7)/4.0d0
        fv(2) = m(0)/9.0d0-m(1)/36.0d0-m(2)/18.0d0 &
            +m(5)/6.0d0-m(6)/6.0d0-m(7)/4.0d0
        fv(3) = m(0)/9.0d0-m(1)/36.0d0-m(2)/18.0d0-m(3)/6.0d0+m(4)/6.0d0 &
            +m(7)/4.0d0
        fv(4) = m(0)/9.0d0-m(1)/36.0d0-m(2)/18.0d0 &
            -m(5)/6.0d0+m(6)/6.0d0-m(7)/4.0d0
        fv(5) = m(0)/9.0d0+m(1)/18.0d0+m(2)/36.0d0+m(3)/6.0d0+m(4)/12.0d0 &
            +m(5)/6.0d0+m(6)/12.0d0+m(8)/4.0d0
        fv(6) = m(0)/9.0d0+m(1)/18.0d0+m(2)/36.0d0-m(3)/6.0d0-m(4)/12.0d0 &
            +m(5)/6.0d0+m(6)/12.0d0-m(8)/4.0d0
        fv(7) = m(0)/9.0d0+m(1)/18.0d0+m(2)/36.0d0-m(3)/6.0d0-m(4)/12.0d0 &
            -m(5)/6.0d0-m(6)/12.0d0+m(8)/4.0d0
        fv(8) = m(0)/9.0d0+m(1)/18.0d0+m(2)/36.0d0+m(3)/6.0d0+m(4)/12.0d0 &
            -m(5)/6.0d0-m(6)/12.0d0-m(8)/4.0d0

    end subroutine flow_populations

    !===============================================================================================
    subroutine thermal_moments(gv,n)
        !$acc routine seq
        real(kind=8), intent(in) :: gv(0:4)
        real(kind=8), intent(out) :: n(0:4)

        n(0) = gv(0)+gv(1)+gv(2)+gv(3)+gv(4)
        n(1) = gv(1)-gv(3)
        n(2) = gv(2)-gv(4)
        n(3) = -4.0d0*gv(0)+gv(1)+gv(2)+gv(3)+gv(4)
        n(4) = gv(1)-gv(2)+gv(3)-gv(4)

    end subroutine thermal_moments

    !===============================================================================================
    subroutine thermal_populations(n,gv)
        !$acc routine seq
        real(kind=8), intent(in) :: n(0:4)
        real(kind=8), intent(out) :: gv(0:4)

        gv(0) = 0.2d0*n(0)-0.2d0*n(3)
        gv(1) = 0.2d0*n(0)+0.5d0*n(1)+0.05d0*n(3)+0.25d0*n(4)
        gv(2) = 0.2d0*n(0)+0.5d0*n(2)+0.05d0*n(3)-0.25d0*n(4)
        gv(3) = 0.2d0*n(0)-0.5d0*n(1)+0.05d0*n(3)+0.25d0*n(4)
        gv(4) = 0.2d0*n(0)-0.5d0*n(2)+0.05d0*n(3)-0.25d0*n(4)

    end subroutine thermal_populations

    !===============================================================================================
    ! 多块输出按不重叠物理分区积分，边缘节点使用裁剪权重；内部权重仍为 h^2。
    !===============================================================================================

    !===============================================================================================
    subroutine calNuRe()
        integer(kind=4) :: b,i,j,k,jm,im
        real(kind=8) :: area,conv,vel2,mass,meanT,nu,re,hot,cold,middle,dTdx,dTdy,tm,um,vm,cellArea
        real(kind=8) :: w(4),dw(4)
        real(kind=8) :: tmin,tmax,rmin,rmax,scale,xmid,ymid,xlo,xhi,ylo,yhi,h

        call update_host_all(.false.)
        area = dble(nx)*dble(ny)
        conv = 0.0d0
        vel2 = 0.0d0
        mass = 0.0d0
        meanT = 0.0d0
        hot = 0.0d0
        cold = 0.0d0
        middle = 0.0d0
        tmin = huge(1.0d0)
        tmax = -tmin
        rmin = tmin
        rmax = -tmin
        scale = lengthUnit/(Thot-Tcold)
        xmid = 0.5d0*nx
        ymid = 0.5d0*ny
        do b = 1,nBlocks
            associate(bl=>blocks(b))
                h = bl%h
                do j = bl%jlo,bl%jhi
                    do i = bl%ilo,bl%ihi
                        if (.not.ieee_is_finite(bl%T(i,j)) .or. .not.ieee_is_finite(bl%rho(i,j)) .or. &
                            .not.ieee_is_finite(bl%u(i,j)) .or. .not.ieee_is_finite(bl%v(i,j)) .or. bl%rho(i,j)<=0.0d0) then
                            write(*,*) 'Invalid state: coarse clock, block, i,j:',itc,b,i,j
                            error stop 'Nonfinite or nonpositive density in owned cells'
                        endif
                        cellArea = bl%dxWeight(i)*bl%dyWeight(j)
#ifdef SideHeatedCell
                        conv = conv+bl%u(i,j)*bl%T(i,j)*cellArea
#else
                        conv = conv+bl%v(i,j)*bl%T(i,j)*cellArea
#endif
                        vel2 = vel2+(bl%u(i,j)**2+bl%v(i,j)**2)*cellArea
                        mass = mass+bl%rho(i,j)*cellArea
                        meanT = meanT+bl%T(i,j)*cellArea
                        tmin = min(tmin,bl%T(i,j))
                        tmax = max(tmax,bl%T(i,j))
                        rmin = min(rmin,bl%rho(i,j))
                        rmax = max(rmax,bl%rho(i,j))
                    enddo
                enddo
                xlo = bl%ownedBox(1)
                xhi = bl%ownedBox(2)
                ylo = bl%ownedBox(3)
                yhi = bl%ownedBox(4)
#ifdef SideHeatedCell
                if (xlo==0.0d0) then
                    do j = bl%jlo,bl%jhi
                        hot = hot+(8.0d0*Thot-9.0d0*bl%T(1,j)+bl%T(2,j))/(3.0d0*h)*bl%dyWeight(j)/dble(ny)
                    enddo
                endif
                if (xhi==dble(nx)) then
                    do j = bl%jlo,bl%jhi
                        cold = cold+(-8.0d0*Tcold+9.0d0*bl%T(bl%ni,j)-bl%T(bl%ni-1,j))/(3.0d0*h)*bl%dyWeight(j)/dble(ny)
                    enddo
                endif
                if (xmid>=xlo .and. xmid<xhi) then
                    call section_weights((xmid-bl%x0)/h+0.5d0,bl%ni,im,w,dw)
                    do j = bl%jlo,bl%jhi
                        tm = sum(w*bl%T(im:im+3,j))
                        um = sum(w*bl%u(im:im+3,j))
                        dTdx = sum(dw*bl%T(im:im+3,j))/h
                        middle = middle+(um*tm/diffusivity-dTdx)*bl%dyWeight(j)/dble(ny)
                    enddo
                endif
#else
                if (ylo==0.0d0) then
                    do i = bl%ilo,bl%ihi
                        hot = hot+(8.0d0*Thot-9.0d0*bl%T(i,1)+bl%T(i,2))/(3.0d0*h)*bl%dxWeight(i)/dble(nx)
                    enddo
                endif
                if (yhi==dble(ny)) then
                    do i = bl%ilo,bl%ihi
                        cold = cold+(-8.0d0*Tcold+9.0d0*bl%T(i,bl%nj)-bl%T(i,bl%nj-1))/(3.0d0*h)*bl%dxWeight(i)/dble(nx)
                    enddo
                endif
                if (ymid>=ylo .and. ymid<yhi) then
                    call section_weights((ymid-bl%y0)/h+0.5d0,bl%nj,jm,w,dw)
                    do i = bl%ilo,bl%ihi
                        tm = sum(w*bl%T(i,jm:jm+3))
                        vm = sum(w*bl%v(i,jm:jm+3))
                        dTdy = sum(dw*bl%T(i,jm:jm+3))/h
                        middle = middle+(vm*tm/diffusivity-dTdy)*bl%dxWeight(i)/dble(nx)
                    enddo
                endif
#endif
            end associate
        enddo
        nu = 1.0d0+conv/area*scale/diffusivity
        re = sqrt(vel2/area)*lengthUnit/viscosity
        open(newunit = k,file = historyFile,status = 'old',position = 'append')
        write(k,'(12(ES24.16E3,1X))') dble(itc)/timeUnit,nu,re,scale*hot,scale*cold,scale*middle, &
            mass,meanT/area,tmin,tmax,rmin,rmax
        close(k)
        write(*,'(a,f12.5,a,es13.5,a,es13.5)') 't_ff=',dble(itc)/timeUnit,' NuVolAvg=',nu,' ReVolRMS=',re
    end subroutine calNuRe

    !===============================================================================================
    subroutine check()
#ifdef steadyFlow
        integer(kind=4) :: b,k,i,j
        real(kind=8) :: du,uu,dt,tt,cellArea

        du = 0.0d0
        uu = 0.0d0
        dt = 0.0d0
        tt = 0.0d0
        call update_host_all(.false.)
        do b = 1,nBlocks
            associate(bl=>blocks(b))
                do j = bl%jlo,bl%jhi
                    do i = bl%ilo,bl%ihi
                        cellArea = bl%dxWeight(i)*bl%dyWeight(j)
                        du = du+cellArea*((bl%u(i,j)-bl%up(i,j))**2+(bl%v(i,j)-bl%vp(i,j))**2)
                        uu = uu+cellArea*(bl%u(i,j)**2+bl%v(i,j)**2)
                        dt = dt+cellArea*(bl%T(i,j)-bl%Tp(i,j))**2
                        tt = tt+cellArea*bl%T(i,j)**2
                    enddo
                enddo
                bl%up = bl%u
                bl%vp = bl%v
                bl%Tp = bl%T
            end associate
        enddo
        errorU = sqrt(du/max(uu,1.0d-300))
        errorT = sqrt(dt/max(tt,1.0d-300))
        open(newunit = k,file = 'Convergence_2DOpenaccMultiblock.dat',status = 'unknown',position = 'append')
        write(k,'(I12,2(1X,ES24.16E3))') itc,errorU,errorT
        close(k)
        write(*,*) 'errorU,errorT:',errorU,errorT
#endif
    end subroutine check

    !===============================================================================================
    subroutine output_Tecplot()
        integer(kind=4) :: k,b,i,j
        character(16) :: num

        pltFileNum = pltFileNum+1
        write(num,'(I10.10)') pltFileNum
        call update_host_all(.false.)
        open(newunit = k,file = pltFolderPrefix//'-'//trim(num)//'.dat',status = 'replace')
        write(k,'(a)') 'VARIABLES="x/L","y/L","u","v","T","rho","h/L","integration_area/L^2"'
        do b = 1,nBlocks
            associate(bl=>blocks(b))
                write(k,'(a,I0,a,I0,a,I0,a,ES24.16E3)') 'ZONE T="block ',b,'", I=',bl%ihi-bl%ilo+1, &
                    ', J=',bl%jhi-bl%jlo+1,', F=POINT, SOLUTIONTIME=',dble(itc)/timeUnit
                do j = bl%jlo,bl%jhi
                    do i = bl%ilo,bl%ihi
                        write(k,'(8(ES24.16E3,1X))') (bl%x0+(dble(i)-0.5d0)*bl%h)/lengthUnit, &
                            (bl%y0+(dble(j)-0.5d0)*bl%h)/lengthUnit,bl%u(i,j),bl%v(i,j),bl%T(i,j),bl%rho(i,j),bl%h/lengthUnit, &
                            bl%dxWeight(i)*bl%dyWeight(j)/lengthUnit**2
                    enddo
                enddo
            end associate
        enddo
        close(k)
    end subroutine output_Tecplot

    !===============================================================================================
    subroutine output_SnapshotFile()
        integer(kind=4) :: k,b
        character(16) :: num

        snapshotFileNum = snapshotFileNum+1
        write(num,'(I10.10)') snapshotFileNum
        call update_host_all(.false.)
        open(newunit = k,file = snapshotFilePrefix//'-'//trim(num)//'.bin',form = 'unformatted',access = 'stream',status = 'replace')
        write(k) snapshotMagic,nBlocks,nx,ny,itc,dble(itc)/timeUnit,lengthUnit
        do b = 1,nBlocks
            associate(bl=>blocks(b))
                write(k) bl%ihi-bl%ilo+1,bl%jhi-bl%jlo+1, &
                    bl%x0+(dble(bl%ilo)-0.5d0)*bl%h,bl%y0+(dble(bl%jlo)-0.5d0)*bl%h,bl%h,bl%ownedBox
                write(k) bl%dxWeight(bl%ilo:bl%ihi),bl%dyWeight(bl%jlo:bl%jhi)
                write(k) bl%u(bl%ilo:bl%ihi,bl%jlo:bl%jhi),bl%v(bl%ilo:bl%ihi,bl%jlo:bl%jhi), &
                    bl%T(bl%ilo:bl%ihi,bl%jlo:bl%jhi),bl%rho(bl%ilo:bl%ihi,bl%jlo:bl%jhi)
            end associate
        enddo
        close(k)
    end subroutine output_SnapshotFile

    !===============================================================================================
    function model_signature() result(sig)
        integer(kind=4) :: sig(12)

        sig = 0
#ifdef steadyFlow
        sig(1) = 1
#endif
#ifdef EnableUseG
        sig(2) = 1
#endif
#ifdef SideHeatedCell
        sig(3) = 1
#endif
#ifdef HorizontalWallsNoslip
        sig(4) = 1
#endif
#ifdef VerticalWallsNoslip
        sig(5) = 1
#endif
#ifdef VerticalWallsPeriodicalU
        sig(6) = 1
#endif
#ifdef HorizontalWallsConstT
        sig(7) = 1
#endif
#ifdef HorizontalWallsAdiabatic
        sig(8) = 1
#endif
#ifdef VerticalWallsConstT
        sig(9) = 1
#endif
#ifdef VerticalWallsAdiabatic
        sig(10) = 1
#endif
#ifdef VerticalWallsPeriodicalT
        sig(11) = 1
#endif
#ifdef SideHeatedHa
        sig(12) = 1
#endif
    end function model_signature

    !===============================================================================================
    function physical_signature() result(sig)
        real(kind=8) :: sig(16)

        sig = [Rayleigh,Prandtl,Mach,Thot,Tcold,Snu,Sq,Qk,Qnu,thermalA,outputSnapshotInterval, &
            reloadFileInterval,outputPltFileInterval,0.0d0,0.0d0,0.0d0]
#ifdef SideHeatedHa
        sig(14) = Ha
        sig(15) = phi
#endif
    end function physical_signature

    !===============================================================================================
    ! 完整检查点：仅在粗细同步时写入
    subroutine output_ReloadFile()
        integer(kind=4) :: k,b
        character(16) :: num
        character(256) :: name

        call update_host_all(.true.)
        write(num,'(I12.12)') itc
        name = reloadFilePrefix//'-'//trim(num)//'.bin'
        open(newunit = k,file = trim(name),access = 'stream',form = 'unformatted',status = 'replace')
        write(k) restartMagic,nx,ny,refineRatio,fineLayerCellsLeftRight,fineLayerCellsBottomTop,overlapCells,nBlocks
        write(k) model_signature(),physical_signature()
        write(k) itc,nextSample,nextReload,nextPlt,snapshotFileNum,pltFileNum,errorU,errorT
        do b = 1,nBlocks
            associate(bl=>blocks(b))
                write(k) bl%ni,bl%nj,bl%ilo,bl%ihi,bl%jlo,bl%jhi,bl%nh,bl%x0,bl%y0,bl%h,bl%ownedBox
                write(k) bl%f,bl%g,bl%u,bl%v,bl%T,bl%rho,bl%Fx,bl%Fy,bl%Bx_prev,bl%By_prev,bl%p
#ifdef steadyFlow
                write(k) bl%up,bl%vp,bl%Tp
#endif
            end associate
        enddo
        close(k)
        ! 完整写完独立编号的 checkpoint 后再更新 latest 指针，旧 checkpoint 仍保留。
        open(newunit = k,file = reloadFilePrefix//'-latest.meta',status = 'replace')
        write(k,'(a)') trim(name)
        close(k)
    end subroutine output_ReloadFile

    !===============================================================================================
    ! 精确续算：检查网格、物理参数、统计分区与历史状态
    subroutine read_restart()
        integer(kind=4) :: k,b,ios,head(7),geom(7),sig(12)
        real(kind=8) :: phys(16),coord(7)
        character(16) :: magic
        character(256) :: name

        open(newunit = k,file = reloadFilePrefix//'-latest.meta',status = 'old',iostat = ios)
        if (ios/=0) error stop 'Missing multiblock latest.meta'
        read(k,'(a)',iostat = ios) name
        close(k)
        if (ios/=0) error stop 'Invalid latest.meta'
        open(newunit = k,file = trim(name),status = 'old',access = 'stream',form = 'unformatted',iostat = ios)
        if (ios/=0) error stop 'Missing multiblock checkpoint'
        read(k,iostat = ios) magic,head
        if (ios/=0) error stop 'Truncated multiblock checkpoint header'
        if (magic/=restartMagic) error stop 'Wrong checkpoint format'
        if (any(head/=[nx,ny,refineRatio,fineLayerCellsLeftRight,fineLayerCellsBottomTop,overlapCells,nBlocks])) &
            error stop 'Restart mesh/refinement mismatch'
        read(k) sig,phys
        if (any(sig/=model_signature()) .or. any(phys/=physical_signature())) &
            error stop 'Restart model/physics/output-cadence mismatch'
        read(k) itc,nextSample,nextReload,nextPlt,snapshotFileNum,pltFileNum,errorU,errorT
        if (itc<0 .or. mod(itc,refineRatio)/=0) error stop 'Restart is not at a synchronized time'
        do b = 1,nBlocks
            associate(bl=>blocks(b))
                read(k) geom,coord
                if (any(geom/=[bl%ni,bl%nj,bl%ilo,bl%ihi,bl%jlo,bl%jhi,bl%nh]) .or. &
                    any(coord/=[bl%x0,bl%y0,bl%h,bl%ownedBox])) error stop 'Restart block layout mismatch'
                read(k,iostat = ios) bl%f,bl%g,bl%u,bl%v,bl%T,bl%rho,bl%Fx,bl%Fy,bl%Bx_prev,bl%By_prev,bl%p
                if (ios/=0) error stop 'Incomplete checkpoint state/history'
#ifdef steadyFlow
                read(k) bl%up,bl%vp,bl%Tp
#endif
            end associate
        enddo
        close(k)
    end subroutine read_restart

    !===============================================================================================
    integer(kind = 4) function scheduled_step(index,interval) result(step)
        integer(kind=4), intent(in) :: index
        real(kind=8), intent(in) :: interval

        step = max(refineRatio,ceiling(dble(index)*interval*timeUnit/dble(refineRatio))*refineRatio)
    end function scheduled_step

    !===============================================================================================
    subroutine check_history()
        integer(kind=4) :: k,ios,n
        real(kind=8) :: values(12),lastTime,expected
        character(512) :: line

        lastTime = -1.0d0
        n = 0
        open(newunit = k,file = historyFile,status = 'old',iostat = ios)
        if (ios/=0) error stop 'Restart requires the matching Nu/Re history'
        do
            read(k,'(a)',iostat = ios) line
            if (ios<0) exit
            if (ios>0) error stop 'Nu/Re history read error'
            if (len_trim(line)==0 .or. line(1:1)=='#') cycle
            read(line,*,iostat = ios) values
            if (ios/=0 .or. .not.all(ieee_is_finite(values))) error stop 'Invalid Nu/Re history row'
            if (values(1)<=lastTime) error stop 'Nonmonotone Nu/Re history'
            n = n+1
            lastTime = values(1)
        enddo
        close(k)
        if (n/=nextSample-1) error stop 'History sample count differs from checkpoint; use the matching history'
        if (n>0) then
            expected = dble(scheduled_step(n,outputSnapshotInterval))/timeUnit
            if (abs(lastTime-expected)>1.0d-10*max(1.0d0,expected)) error stop 'Checkpoint/history time mismatch'
        endif
    end subroutine check_history

    !===============================================================================================
    subroutine average_window(t0,t1,result,coverage)
        real(kind=8), intent(in) :: t0,t1
        real(kind=8), intent(out) :: result(5),coverage
        integer(kind=4) :: k,ios
        real(kind=8) :: prev(12),row(12),aa,bb,dt,ra(5),rb(5),left(5),right(5)
        logical :: hasPrev
        character(512) :: line

        result = 0.0d0
        coverage = 0.0d0
        hasPrev = .false.
        open(newunit = k,file = historyFile,status = 'old')
        do
            read(k,'(a)',iostat = ios) line
            if (ios<0) exit
            if (ios>0) error stop 'Nu/Re history read failure'
            if (len_trim(line)==0 .or. line(1:1)=='#') cycle
            read(line,*,iostat = ios) row
            if (ios/=0 .or. .not.all(ieee_is_finite(row))) error stop 'Invalid history during averaging'
            if (hasPrev) then
                if (row(1)<=prev(1)) error stop 'Nonmonotone statistics time'
                aa = max(t0,prev(1))
                bb = min(t1,row(1))
                dt = row(1)-prev(1)
                if (bb>aa) then
                    ! 对 Re^2 积分再开根号，保持 sqrt(<u^2+v^2>_{V,t}) 的定义。
                    left = prev(2:6)
                    right = row(2:6)
                    left(2) = left(2)**2
                    right(2) = right(2)**2
                    ra = left+(right-left)*(aa-prev(1))/dt
                    rb = left+(right-left)*(bb-prev(1))/dt
                    result = result+0.5d0*(ra+rb)*(bb-aa)
                    coverage = coverage+bb-aa
                endif
            endif
            prev = row
            hasPrev = .true.
        enddo
        close(k)
        if (coverage>0.0d0) then
            result = result/coverage
            result(2) = sqrt(max(0.0d0,result(2)))
        endif
    end subroutine average_window

    !===============================================================================================
    subroutine output_unsteady_NuRe_postprocess()
#ifdef unsteadyFlow
        integer(kind=4) :: k
        real(kind=8) :: allMean(5),firstMean(5),lastMean(5),c0,c1,c2,relative(5)

        call average_window(unsteadyAverageStartTf,unsteadyAverageEndTf,allMean,c0)
        call average_window(unsteadyAverageStartTf,unsteadyAverageMidTf,firstMean,c1)
        call average_window(unsteadyAverageMidTf,unsteadyAverageEndTf,lastMean,c2)
        open(newunit = k,file = 'NuReStatistics_2DOpenaccMultiblock.dat',status = 'replace')
        write(k,*) 'Window t_ff:',unsteadyAverageStartTf,unsteadyAverageEndTf
        write(k,*) 'Covered durations, full/first/last:',c0,c1,c2
        if (abs(c0-(unsteadyAverageEndTf-unsteadyAverageStartTf))>1.0d-8 .or. c1<=0.0d0 .or. c2<=0.0d0) then
            write(k,*) 'INCOMPLETE: requested statistics window is not fully covered; no final mean is reported.'
        else
            relative = abs(lastMean-firstMean)/max(abs(allMean),1.0d-30)
            write(k,*) 'Columns: NuVolAvg Re_rms_space_time Nu_hot Nu_cold Nu_middle'
            write(k,'(a,5ES24.16E3)') 'whole:',allMean
            write(k,'(a,5ES24.16E3)') 'first:',firstMean
            write(k,'(a,5ES24.16E3)') 'last :',lastMean
            write(k,'(a,5ES24.16E3)') 'relative half-window difference:',relative
        endif
        close(k)
#endif
    end subroutine output_unsteady_NuRe_postprocess

    end module commondata

    !===============================================================================================
    ! 主程序：初始化、时间推进、输出及释放设备数据
    !===============================================================================================
    program main
        use openacc
        use commondata
        implicit none
        integer(kind=4) :: finalStep
        integer(kind=8) :: clockStart,clockEnd,clockRate

        call acc_init(acc_device_default)
        write(*,*) 'Visible OpenACC devices:',acc_get_num_devices(acc_device_default)
        call initial()
        call enter_data_2d_openacc()
        call system_clock(clockStart,clockRate)
        finalStep = ((itc_max+refineRatio-1)/refineRatio)*refineRatio
        do while (itc<finalStep)
#ifdef steadyFlow
            if (errorU<=epsU .and. errorT<=epsT) exit
#endif
            call advance_multiblock()
#ifdef steadyFlow
            if (mod(itc,2000)==0) call check()
#endif
            if (itc>=scheduled_step(nextSample,outputSnapshotInterval)) then
                call calNuRe()
                nextSample = nextSample+1
                if (outputSnapshotFile==1) call output_SnapshotFile()
            endif
            if (itc>=scheduled_step(nextPlt,outputPltFileInterval)) then
                nextPlt = nextPlt+1
                if (outputPltFile==1) call output_Tecplot()
            endif
            if (itc>=scheduled_step(nextReload,reloadFileInterval)) then
                nextReload = nextReload+1
                if (outputReloadFile==1) call output_ReloadFile()
            endif
        enddo
        !$acc wait(1)
        call system_clock(clockEnd)
        write(*,*) 'Elapsed seconds:',dble(clockEnd-clockStart)/dble(clockRate)
        call output_unsteady_NuRe_postprocess()
        if (outputPltFile==1) call output_Tecplot()
        if (outputReloadFile==1) call output_ReloadFile()
        call exit_data_2d_openacc()
    end program main
