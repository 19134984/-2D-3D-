!=============================================================
!!!    注释区，代码描述
!!!    二维浮力驱动自然对流 MPI+OpenACC 并行版本
!!!    LBM方法
!!!    MRT-LBE
!=============================================================


!=============================================================
!   自定义宏，一些选项的开关
#define steadyFlow    
!#define unsteadyFlow

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
!#define RayleighBenardCell
!#define HorizontalWallsConstT
!#define VerticalWallsAdiabatic
!#define VerticalWallsPeriodicalT



!   温度边界(for Side Heated Cell)，包括水平边界温度不可穿透，垂直边界恒温,侧壁加热加磁场
#define SideHeatedCell
#define HorizontalWallsAdiabatic
#define VerticalWallsConstT
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
!#define EnableUseG
!启用旧温度算法
#define EnableLegacyThermalScheme

!   温度算法宏的选择
#if defined(EnableUseG) && defined(EnableLegacyThermalScheme)
#error "Choose only one thermal scheme: EnableUseG or EnableLegacyThermalScheme"
#endif
#if !defined(EnableUseG) && !defined(EnableLegacyThermalScheme)
#error "Define one thermal scheme: EnableUseG or EnableLegacyThermalScheme"
#endif

!   MPI 笛卡尔分解宏：一维 y-slab 或二维 x/y block 二选一
!#define MpiDecomp1D
#define MpiDecomp2D

#if defined(MpiDecomp1D) && defined(MpiDecomp2D)
#error "Choose only one MPI decomposition mode: MpiDecomp1D or MpiDecomp2D"
#endif
#if !defined(MpiDecomp1D) && !defined(MpiDecomp2D)
#error "Define one MPI decomposition mode: MpiDecomp1D or MpiDecomp2D"
#endif

!   OpenACC GPU 分配宏：单节点测试按全局 NPROC 推出节点内 rank 数；多节点集群按 NodeCount 推出每节点 rank 数
!#define AccGpuSingleNode
#define AccGpuMultiNode

#if defined(AccGpuSingleNode) && defined(AccGpuMultiNode)
#error "Choose only one OpenACC GPU mode: AccGpuSingleNode or AccGpuMultiNode"
#endif
#if !defined(AccGpuSingleNode) && !defined(AccGpuMultiNode)
#error "Define one OpenACC GPU mode: AccGpuSingleNode or AccGpuMultiNode"
#endif

!   MPI 周期通信要求同一方向上的速度/温度周期性一致，否则 halo 语义会混乱。
#if defined(VerticalWallsPeriodicalU) && !defined(VerticalWallsPeriodicalT)
#error "MPI x-periodic decomposition requires VerticalWallsPeriodicalU and VerticalWallsPeriodicalT together"
#endif
#if !defined(VerticalWallsPeriodicalU) && defined(VerticalWallsPeriodicalT)
#error "MPI x-periodic decomposition requires VerticalWallsPeriodicalU and VerticalWallsPeriodicalT together"
#endif



!   自定义宏结束
!=============================================================


!=============================================================
!   全局模块
    module commondata
        use mpi
        implicit none
        !===============================================================================================
        ! 是否在计算前从旧算例重启
        integer(kind=4), parameter :: loadInitField=0   ! 0: 不重启；1: 按 latest .meta 自动续算

        ! 正常断电续算只需要设置 loadInitField=1；
        ! 代码只读取 <reloadFilePrefix>-latest.meta，并从里面找到最新的 .bin。
        ! 正常续算不用改 reloadFileNum；只有 latest .meta 缺失时，
        ! 才手动设置 reloadFileNum 作为保守推断编号。
        integer(kind=4) :: reloadFileNum=0              ! latest .meta 存在时会被覆盖；meta 缺失时作为手工兜底编号
        !===============================================================================================
        real(kind=8) :: reloadDimensionlessTime=0.0d0   ! 续算前已累计的 t_ff；优先从 latest .meta 读取，meta 缺失时由代码推断
        integer(kind=4) :: restartItcOffset=0           ! 续算前已累计的格子步数；优先从 latest .meta 读取，meta 缺失时由代码推断
        logical :: reloadMetadataLoaded=.false.         ! 标记是否成功读取 reload 元数据文件
        !===============================================================================================

        !===============================================================================================
        ! MPI+OpenACC 并行配置
        integer(kind=4), parameter :: NodeCount=2       ! 参与当前算例的节点数；单节点模式应保持为 1
        integer(kind=4), parameter :: GPUS_PER_NODE=1   ! 每个节点分配给本程序的 GPU 数；当前版本按一 GPU 一 MPI rank 运行
        integer(kind=4) :: MYID, NPROC
        integer(kind=4) :: NPROC_NODE                   ! 按配置期望的每节点 MPI rank 数；一 GPU 一 rank 时等于 GPUS_PER_NODE
        integer(kind=4) :: IERR
        integer(kind=4) :: COMM2D=MPI_COMM_NULL
        integer(kind=4) :: dims(2)
        integer(kind=4) :: MY_COORD(2)
        integer(kind=4) :: left, right
        integer(kind=4) :: down, top
#ifdef VerticalWallsPeriodicalU
        logical :: periods(2)=(/.true.,.false./)    ! x 方向周期；y 方向仍是上下物理壁面
#else
        logical :: periods(2)=(/.false.,.false./)   ! 默认 x/y 都不是 MPI 周期方向
#endif
        logical :: isRoot
        logical :: hasLeftBoundary, hasRightBoundary
        logical :: hasBottomBoundary, hasTopBoundary
        logical :: accDataResident=.false.                 ! 标记局部场是否已经进入 OpenACC 设备端，避免重复 exit data
        integer(kind=4) :: xLocalCount, yLocalCount
        integer(kind=4) :: xStartGlobal, xEndGlobal
        integer(kind=4) :: yStartGlobal, yEndGlobal
        integer(kind=4), allocatable :: XLocalCountAll(:), YLocalCountAll(:)
        integer(kind=4), allocatable :: XStartGlobalAll(:), YStartGlobalAll(:)
        !===============================================================================================

        !===============================================================================================
        ! 无量纲参数
        integer(kind=4), parameter :: nx=256, ny=256     !格子网格

        ! 本文件采用erf非均匀网格。
        ! erf网格拉伸强度。数值越大, 节点越向两侧物理壁面聚集。
        real(kind=8), parameter :: ISLBM_StretchA=1.5d0

        ! raw坐标中第一个内部流体节点的位置rawX(1)/rawY(1), 也就是近壁的1个lu。
        real(kind=8), parameter :: ISLBM_DxMinRaw=0.5d0*(1.0d0 + &
            erf(ISLBM_StretchA*(1.0d0/dble(nx+1)-0.5d0))/erf(0.5d0*ISLBM_StretchA))
        real(kind=8), parameter :: ISLBM_DyMinRaw=0.5d0*(1.0d0 + &
            erf(ISLBM_StretchA*(1.0d0/dble(ny+1)-0.5d0))/erf(0.5d0*ISLBM_StretchA))
#ifdef SideHeatedCell
        ! half-way壁面放在0.5*rawX(1)或0.5*rawY(1), 因此有效长度为1-ISLBM_Dx/DyMinRaw。
        ! ISLBM有效长度含多少个近壁lu: lengthUnit=(rightWall-leftWall)/ISLBM_DxMinRaw。
        real(kind=8), parameter :: lengthUnit=(1.0d0-ISLBM_DxMinRaw)/ISLBM_DxMinRaw
#else
        ! ISLBM有效长度含多少个近壁lu: lengthUnit=(topWall-bottomWall)/ISLBM_DyMinRaw。
        real(kind=8), parameter :: lengthUnit=(1.0d0-ISLBM_DyMinRaw)/ISLBM_DyMinRaw
#endif
        real(kind=8), parameter :: pi = acos(-1.0d0)

        real(kind=8), parameter :: Rayleigh=1.0d6        
        real(kind=8), parameter :: Prandtl=0.71d0       
        real(kind=8), parameter :: Mach=0.1d0
        real(kind=8), parameter :: Thot=0.5d0, Tcold=-0.5d0
        real(kind=8), parameter :: Tref=0.5d0*(Thot+Tcold)
        real(kind=8), parameter :: tauf=0.5d0+Mach*lengthUnit*dsqrt(3.0d0*Prandtl/Rayleigh) 
        real(kind=8), parameter :: viscosity=(tauf-0.5d0)/3.0d0
        real(kind=8), parameter :: diffusivity=viscosity/Prandtl
        

        ! velocityScaleCompare is used only in velocity-related post-processing to convert lattice velocity
        ! to the nondimensional velocity scale adopted by the reference paper being compared.
        real(kind=8), parameter :: velocityScaleCompare=lengthUnit/diffusivity
        ! 默认采用热扩散标度 UL/kappa；若要按自由落体标度比较，可改为 1.0d0/velocityUnit
        
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1


#ifdef  SideHeatedHa
        real(kind=8), parameter :: Ha=20.0d0                           !磁场强度
        real(kind=8), parameter :: phi=(0.0d0)*(pi/180.0d0)            !磁场角度，以水平向右为0，修改0.0d0即可
        real(kind=8), parameter :: B2sigemarho=(Ha**2*viscosity)/(lengthUnit*lengthUnit)  !动量方程上的源项系数
#endif

        ! 高阶矩参数修正
        real(kind=8), parameter :: paraA=20.0d0*dsqrt(3.0d0)*diffusivity-4.0d0



        ! 浮力项参数
        real(kind=8), parameter :: gBeta1=Rayleigh*viscosity*diffusivity/lengthUnit
        real(kind=8), parameter :: gBeta=gBeta1/lengthUnit/lengthUnit             !gbetaΔT
        
        real(kind=8), parameter :: timeUnit=dsqrt(lengthUnit/gBeta)      !无量纲时间
        real(kind=8), parameter :: velocityUnit=dsqrt(gBeta*lengthUnit)  !无量纲速度
    
        real(kind=8), parameter :: Snu=1.0d0/tauf, Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)  !动量的多松弛系数

#ifdef EnableLegacyThermalScheme
        real(kind=8), parameter :: Qk=3.0d0-dsqrt(3.0d0), Qnu=4.0d0*dsqrt(3.0d0)-6.0d0             !旧温度算法的多松弛系数
        real(kind=8), parameter :: thermalGeqCoeff=10.0d0/(4.0d0+paraA)
#else
        real(kind=8), parameter :: taug = 0.5d0 + (tauf - 0.5d0)/Prandtl
        real(kind=8), parameter :: Qnu = 1.0d0, Qk = 1.0d0/taug
        real(kind=8), parameter :: thermalGeqCoeff=3.0d0
#endif
        !===============================================================================================          
        
        !===============================================================================================
        real(kind=8), parameter :: epsU=1.0d-7, epsT=1.0d-7    ! 稳态收敛阈值   

#ifdef steadyFlow
        real(kind=8), parameter :: outputSnapshotInterval=10.0d0   ! 快照和 Nu/Re 时间序列采样间隔（单位：t_ff）
        real(kind=8), parameter :: reloadFileInterval=100.0d0  ! f/g 重启文件输出间隔（单位：t_ff）
        real(kind=8), parameter :: outputPltFileInterval=100.0d0  ! Tecplot 文件输出间隔（单位：t_ff）
        integer(kind=4), parameter :: dimensionlessTimeMax=int(12000.0d0/outputSnapshotInterval)
        integer(kind=4), parameter :: outputSnapshotFile=1   ! 是否输出后处理快照文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputPltFile=1   ! 是否输出 plt 文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputReloadFile=1 ! 是否周期输出 f/g 重启文件：0=不输出，1=输出
        integer(kind=4), parameter :: itc_max=20000000  ! 稳态：最大格子步，实际仍由 errorU/errorT 提前停止
#endif

#ifdef unsteadyFlow
        real(kind=8), parameter :: outputSnapshotInterval=0.5d0   ! uvTrho 快照和 Nu/Re 时间序列采样间隔（单位：t_ff）
        real(kind=8), parameter :: reloadFileInterval=100.0d0  ! f/g 重启文件输出间隔（单位：t_ff）
        real(kind=8), parameter :: outputPltFileInterval=100.0d0  ! Tecplot 文件周期输出间隔（单位：t_ff）
        real(kind=8), parameter :: unsteadyRunDuration=1000.0d0  ! 非稳态总目标时长，续算时只补足到该 t_ff
        ! 以下三个参数只控制非稳态结束后的 Nu/Re 统计平均窗口，不改变推进时长或采样频率。
        ! 时间以整个非稳态算例的绝对 t_ff 计，续算统计会包含旧文件中已有的历史数据。
        real(kind=8), parameter :: unsteadyAverageStartTf=0.5d0*unsteadyRunDuration  ! 平均窗口起点
        real(kind=8), parameter :: unsteadyAverageEndTf=unsteadyRunDuration          ! 平均窗口终点
        real(kind=8), parameter :: unsteadyAverageMidTf=0.5d0*(unsteadyAverageStartTf+unsteadyAverageEndTf) ! 前/后半分界
        integer(kind=4), parameter :: unsteadySampleCount=max(1, &
          int(unsteadyRunDuration/outputSnapshotInterval+0.5d0))
        ! 计数器，输出多少次快照
        integer(kind=4), parameter :: dimensionlessTimeMax=unsteadySampleCount
        integer(kind=4), parameter :: outputSnapshotFile=1   ! 是否输出后处理快照文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputPltFile=1   ! 是否输出 plt 文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputReloadFile=1 ! 是否周期输出 f/g 重启文件：0=不输出，1=输出
        integer(kind=4), parameter :: itc_max=max(1, int(unsteadyRunDuration*timeUnit+0.5d0)) ! 非稳态：由总 t_ff 自动换算格子步
#endif

        integer(kind=4) :: snapshotFileNum, pltFileNum  ! 快照/plt 输出文件的计数器
        ! 每次调用对应输出子程序时递增（用于文件名编号）

        integer(kind=4) :: dimensionlessTime
        integer(kind=4) :: outputSnapshotIntervalItc
        integer(kind=4) :: reloadFileIntervalItc, outputPltFileIntervalItc
        ! 统计/输出时间点编号（与 outputSnapshotInterval 对应）：
        ! 每调用一次 calNuRe() 就 dimensionlessTime = dimensionlessTime + 1
        ! 用于索引 NuVolAvg/ReVolAvg 数组，并用于输出的时间轴：t = reloadDimensionlessTime + dimensionlessTime*outputSnapshotInterval（单位：t_ff）

        real(kind=8) :: NuVolAvg(0:dimensionlessTimeMax), ReVolAvg(0:dimensionlessTimeMax)
        ! 体平均 Nu 和 Re 的时间序列缓存
        ! 只有在启用并调用 calNuRe() 的情况下这些数组才会被真正填充
  
        character(len=100) :: snapshotFilePrefix="buoyancyCavity2DOpenaccMpiSnapshot"
        ! 快照输出文件前缀（实际文件名形如：<snapshotFilePrefix>-<编号>.bin）

        character(len=100) :: pltFolderPrefix="buoyancyCavity2DOpenaccMpiTecplot"
        ! plt 输出文件前缀（实际文件名形如：<pltFolderPrefix>-<编号>.plt）

        character(len=100) :: reloadFilePrefix="reloadFile2DOpenaccMpi"
        ! 重启读取文件的前缀；latest meta 模式实际读取 meta 中记录的 <reloadFilePrefix>-<编号>.bin
        character(len=100) :: settingsFile="SimulationSettings2DOpenaccMpi.txt"
        !===============================================================================================

        !===============================================================================================
        !计算中需要的相关参数
        real(kind=8) :: errorU, errorT
        
        real(kind=8) :: xp(0:nx+1), yp(0:ny+1)      !无量纲的全局坐标数组，包括物理边界
        real(kind=8) :: quadWidthX(1:nx), quadWidthY(1:ny)          !非均匀网格 midpoint-rule 积分权重
        real(kind=8) :: quadSumX, quadSumY, quadSumArea
        ! 归一化坐标中的1个lattice unit; 迁移时用xp(i)-ex(alpha)*ISLBM_LatticeUnit找上游点。
        real(kind=8), parameter :: ISLBM_LatticeUnit=1.0d0/lengthUnit
        real(kind=8), allocatable :: u(:,:), v(:,:), T(:,:), rho(:,:)

#ifdef steadyFlow
        real(kind=8), allocatable :: up(:,:), vp(:,:), Tp(:,:)   !存储之前的数据，用来算收敛判据
#endif
        real(kind=8), allocatable :: f(:,:,:), f_post(:,:,:)
        real(kind=8), allocatable :: g(:,:,:), g_post(:,:,:)
        real(kind=8), allocatable :: Fx(:,:), Fy(:,:)

        real(kind=8), allocatable :: Bx_prev(:,:), By_prev(:,:)
        real(kind=8), allocatable :: u_all(:,:), v_all(:,:), T_all(:,:), rho_all(:,:)
        real(kind=8), allocatable :: Bx_prev_all(:,:), By_prev_all(:,:)
        real(kind=8), allocatable :: f_all(:,:,:), g_all(:,:,:)
        ! halo 通信缓冲区按当前 rank 的局部块尺寸分配。
        ! 后面 allocate_halo_buffers_2d_openacc_mpi() 按局部尺寸 allocate 并 enter data create 到 GPU。
        ! pointer 的作用：
        !   1) 像 allocatable 一样，可以等 xLocalCount/yLocalCount 确定后再按局部尺寸 allocate；
        !   2) 当前 OpenACC 编译器允许 pointer 出现在 host_data use_device(...) 中，
        !      这样 GPU-aware MPI 可以拿到这些 halo buffer 的设备端地址。
        ! 这里 pointer 只当作“可动态分配且可交给 use_device”的数组句柄。
        ! =>null() 表示初始时这个 pointer 不指向任何数组；后面 allocate(...) 后才变成有效数组。
        ! 这样 associated(fHaloSendDown) 的判断是可靠的，避免未初始化 pointer 状态不确定。
        real(kind=8), pointer :: fHaloSendDown(:)=>null(), fHaloSendUp(:)=>null()
        real(kind=8), pointer :: fHaloRecvDown(:)=>null(), fHaloRecvUp(:)=>null()
        real(kind=8), pointer :: fHaloSendLeft(:)=>null(), fHaloSendRight(:)=>null()
        real(kind=8), pointer :: fHaloRecvLeft(:)=>null(), fHaloRecvRight(:)=>null()
        real(kind=8), pointer :: gHaloSendDown(:)=>null(), gHaloSendUp(:)=>null()
        real(kind=8), pointer :: gHaloRecvDown(:)=>null(), gHaloRecvUp(:)=>null()
        real(kind=8), pointer :: gHaloSendLeft(:)=>null(), gHaloSendRight(:)=>null()
        real(kind=8), pointer :: gHaloRecvLeft(:)=>null(), gHaloRecvRight(:)=>null()
        integer(kind=4), allocatable :: streamInterpIndexX(:,:,:), streamInterpIndexY(:,:,:)
        real(kind=8), allocatable :: streamInterpWeightX(:,:,:), streamInterpWeightY(:,:,:)
        logical, allocatable :: streamInterpValidX(:,:), streamInterpValidY(:,:)
        integer(kind=4) :: streamStencilFallbackLocal, streamStencilFallbackGlobal
        ! OpenACC halo buffer 说明：
        !   这些数组不是用全局 nx/ny 静态开大数组，而是在每个 rank 得到局部网格后按局部尺寸分配。
        !   后续 GPU kernel 可以用 present(fHaloSendDown,...) 直接访问这些设备端缓冲区。
        !   每次 halo 交换前，代码会在 GPU 上把 f_post/g_post 的边界层 pack 到 fHalo*/gHalo*，
        !   MPI_SENDRECV 收到的数据放到对应 Recv 缓冲区后，再在 GPU 上 unpack 回 halo 层。
        !   Down/Up/Left/Right 明确表示通信方向；host_data use_device(...) 会把 GPU 地址交给 GPU-aware MPI。

        integer(kind=4) :: itc
#ifdef EnableUseG
        logical, parameter :: useG = .true.            !M1G 开关
#else
        logical, parameter :: useG = .false.           !M1G 开关
#endif

#ifdef EnableLegacyThermalScheme
        logical, parameter :: useLegacyThermalScheme = .true.            !旧算法的逻辑变量
#else
        logical, parameter :: useLegacyThermalScheme = .false.           
#endif
        real(kind=8) :: Nu_global, Nu_hot, Nu_cold, Nu_middle    !平均Nu，全场，侧壁以及中线
        real(kind=8) :: Nu_hot_max, Nu_hot_min, Nu_hot_max_position, Nu_hot_min_position    !左侧壁面的最大最小Nu，以及对应的位置
        
        
        !格子离散速度和权重
        integer(kind=4) :: ex(0:8), ey(0:8)
        data ex/0, 1, 0, -1,  0, 1, -1, -1,  1/
        data ey/0, 0, 1,  0, -1, 1,  1, -1, -1/
        real(kind=8) :: omega(0:8), omegaT(0:4)
        !===============================================================================================

    end module commondata

!   全局模块结束
!=============================================================


!=============================================================

    program main

    use openacc
    use mpi
    use commondata
    implicit none
    real(kind=8) :: timeStart, timeEnd
    character(len=24) :: ctime, string
    INTEGER(kind=4) :: time
    real(kind=8) :: timeStart2, timeEnd2
    integer(kind=4) :: numAccDevices, accDeviceId
    integer(kind=4) :: nodeComm, nodeLocalRank, nodeLocalSize
#ifdef unsteadyFlow
    integer(kind=4) :: nextSampleItc
    integer(kind=4) :: nextSampleAbsItc
    integer(kind=4) :: unsteadyItcRemaining
#endif
    

    !===============================================================================================
    !设置 MPI 笛卡尔分解和每个 rank 对应的 OpenACC 设备
    call MPI_INIT(IERR)         !MPI 初始化
    call init_mpi_cartesian()   !建立笛卡尔通信器，确定当前 rank 的局部网格范围、物理边界标记和邻居 rank
    ! GPU-aware MPI 需要每个 MPI rank 先绑定到正确 GPU；这里参考 side_acc_mono.F90 的 local_rank 绑定方式。
    ! 单节点和多节点都不能直接用全局 rank；应在同一计算节点内取 node-local rank。
    ! MPI_COMM_SPLIT_TYPE 参数含义：
    !   COMM2D               : 输入通信器；只把二维笛卡尔通信器里的 rank 拿来分组。
    !   MPI_COMM_TYPE_SHARED : 分组类型；按共享内存节点分组，通常同一计算节点上的 rank 会进入同一个 nodeComm。
    !   MYID                 : key 参数；决定 nodeComm 内部 rank 的排序，这里让节点内 rank 顺序跟 COMM2D 的 MYID 顺序一致。
    !   MPI_INFO_NULL        : info 参数；不向 MPI 提供额外提示，使用默认分组行为。
    !   nodeComm             : 输出的新通信器；每个节点得到一个只包含本节点 rank 的通信器。
    !   IERR                 : MPI 错误码。
    ! 建立单个节点通信器
    call MPI_COMM_SPLIT_TYPE(COMM2D, MPI_COMM_TYPE_SHARED, MYID, MPI_INFO_NULL, nodeComm, IERR)
    call MPI_COMM_RANK(nodeComm, nodeLocalRank, IERR)  ! 当前 rank 在本节点内的编号，用来选择 GPU
    call MPI_COMM_SIZE(nodeComm, nodeLocalSize, IERR)  ! 本节点内共有多少个 rank，用来复核一 GPU 一 rank
    call MPI_COMM_FREE(nodeComm, IERR)                 ! 已取得本节点信息，释放临时节点通信器
    if(nodeLocalSize.NE.NPROC_NODE) then
        if(nodeLocalRank.EQ.0) then
#ifdef AccGpuSingleNode
            write(*,*) "Error: AccGpuSingleNode mode expects all MPI ranks on one node."
#endif
#ifdef AccGpuMultiNode
            write(*,*) "Error: actual local ranks do not match GPUS_PER_NODE in AccGpuMultiNode mode."
#endif
            write(*,*) "Expected ranks per node =", NPROC_NODE, "; actual local ranks =", nodeLocalSize
        endif
        call MPI_ABORT(MPI_COMM_WORLD, 1006, IERR)
    endif
    ! acc_get_num_devices 查询的是“当前 MPI rank 可见”的 OpenACC 设备数。
    ! 它通常等于当前 rank 所在节点上可见的 GPU 数，不是所有节点 GPU 数之和。
    ! 若运行脚本设置了 CUDA_VISIBLE_DEVICES 等可见性限制，这里返回的是限制后的可见设备数。
    numAccDevices = acc_get_num_devices(acc_device_default)
    if(numAccDevices.LE.0) then
        if(nodeLocalRank.EQ.0) then
            write(*,*) "Error: no visible OpenACC device on this node."
        endif
        call MPI_ABORT(MPI_COMM_WORLD, 1004, IERR)
    endif
    if(nodeLocalSize.GT.numAccDevices) then
        if(nodeLocalRank.EQ.0) then
            write(*,*) "Error: one-rank-per-GPU mode requires local MPI ranks <= visible GPUs."
            write(*,*) "Local MPI ranks =", nodeLocalSize, "; visible OpenACC devices =", numAccDevices
        endif
        call MPI_ABORT(MPI_COMM_WORLD, 1005, IERR)
    endif
    ! 这里强制一 GPU 一 rank；前面已保证 nodeLocalSize <= numAccDevices，
    ! 所以 nodeLocalRank 本身就是要绑定的 GPU 编号。
    accDeviceId = nodeLocalRank
    ! acc_set_device_num 参数含义：
    !   accDeviceId        : 当前 MPI rank 要使用的 GPU 编号；这里等于 nodeLocalRank。
    !   acc_device_default : OpenACC 默认设备类型，由编译器/运行环境决定具体后端。
    ! 作用：先为当前 rank 选定后续 OpenACC kernel 和数据区使用的设备。
    call acc_set_device_num(accDeviceId, acc_device_default)
    ! acc_init 参数含义：
    !   acc_device_default : 初始化 OpenACC 默认设备类型对应的运行时和设备上下文。
    ! 作用：在已选定 accDeviceId 后初始化 OpenACC，避免多个 rank 默认先落到同一张 GPU 上。
    call acc_init(acc_device_default)
    
    call allocate_halo_buffers_2d_openacc_mpi()  ! xLocalCount/yLocalCount 已确定，按局部尺寸分配 GPU-aware MPI halo buffer
    if(loadInitField.EQ.1) then
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,*) " "
        write(00,*) "================ Restart continuation begins ================"
    else
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')   !新算例清旧日志已在 init_mpi_cartesian 完成
    endif
    string = ctime( time() )                      !ctime把 time() 返回的时间戳转换成可读的字符串
    write(00,*) 'Start: ', string                 !什么时候开始计算
    write(00,*) "Starting MPI + OpenACC >>>>>>"
    write(00,*) "MPI ranks from mpiexec:", NPROC
    write(00,*) "MPI dims:", dims(1), dims(2)
    write(00,*) "MPI rank/MY_COORD/localRange/neighbors:", MYID, MY_COORD(1), MY_COORD(2), &
        xStartGlobal, xEndGlobal, yStartGlobal, yEndGlobal, &
        left, right, down, top
#ifdef AccGpuSingleNode
    write(00,*) "OpenACC GPU mode: single-node"
#endif
#ifdef AccGpuMultiNode
    write(00,*) "OpenACC GPU mode: multi-node"
    write(00,*) "Configured NodeCount:", NodeCount
#endif
    write(00,*) "Configured GPUs per node:", GPUS_PER_NODE
    write(00,*) "Expected/actual MPI ranks per node:", NPROC_NODE, nodeLocalSize
    write(00,*) "Node-local MPI rank:", nodeLocalRank
    write(00,*) "Visible OpenACC devices on this node:", numAccDevices
    write(00,*) "OpenACC device id for this rank:", accDeviceId
    close(00)
    !===============================================================================================


    !===============================================================================================
    ! Initialization
    call initial()
    call enter_data_2d_openacc()
#ifdef unsteadyFlow
    ! 非稳态的 itc_max 是整个算例的总目标步数；
    ! 续算时 restartItcOffset 是旧算例已经完成的步数，只推进剩余步数。
    unsteadyItcRemaining = max(0, itc_max - restartItcOffset)
#endif

    !===============================================================================================
    !-----------------------------------------------------------------------------------------------

    call MPI_BARRIER(COMM2D, IERR)   !MPI 屏障同步
    call CPU_TIME(timeStart)         !当前进程累计消耗的 CPU 时间,包括并行
    timeStart2 = MPI_WTIME()         !MPI 墙钟时间
#ifdef steadyFlow
    ! 只要 (errorU > epsU 或 errorT > epsT) 且 itc <= itc_max，就继续循环
        do while( ((errorU.GT.epsU).OR.(errorT.GT.epsT)).AND. &
          (itc.LE.itc_max) )
                                                                              !换成if，就是 errorU > epsU .and. errorT > epsT
#endif
#ifdef unsteadyFlow
    do while( itc.LT.unsteadyItcRemaining )       !非稳态：续算时只推进到 unsteadyRunDuration 对应的总格子步
#endif

        itc = itc+1
        
        call collision()

        call exchange_f_post_halo_mpi()

        call streaming()

        call bounceback()

        call macro()

        call collisionT()

        call exchange_g_post_halo_mpi()

        call streamingT()

        call bouncebackT()
        
        call macroT()

#ifdef steadyFlow
        ! 周期输出按累计格子步判断，续算时才能接回不断电运行应有的输出节奏。
        if(MOD(restartItcOffset+itc,2000).EQ.0) call check()
        if( (outputPltFile.EQ.1).AND.(MOD(restartItcOffset+itc, outputPltFileIntervalItc).EQ.0) ) then
            call output_Tecplot()  !稳态模式下的可选周期 Tecplot 输出
        endif
        if( (outputSnapshotFile.EQ.1).AND.(MOD(restartItcOffset+itc, outputSnapshotIntervalItc).EQ.0) ) then
            call output_SnapshotFile()  !稳态模式下的可选周期 uvTrho 快照输出
        endif
        if( (outputReloadFile.EQ.1).AND.(MOD(restartItcOffset+itc, reloadFileIntervalItc).EQ.0) ) then
            call output_ReloadFile()      !稳态模式下的可选周期 f/g 重启文件输出
        endif
#endif

#ifdef unsteadyFlow
        do while( (reloadDimensionlessTime + real(dimensionlessTime,kind=8)*outputSnapshotInterval) &
            .LT. unsteadyRunDuration )   !判断本次还要不要继续采样
            ! 每个目标采样时刻都按绝对 t_ff 换算到本次运行段的 itc，续算时不会重复旧样本。
            nextSampleAbsItc = max(1, int((reloadDimensionlessTime + &
                real(dimensionlessTime+1,kind=8)*outputSnapshotInterval)*timeUnit+0.5d0)) !算下一个采样点对应的绝对格子步
            nextSampleItc = max(1, nextSampleAbsItc - restartItcOffset) !换成本次续算运行里的局部步数
            if(itc.LT.nextSampleItc) exit !如果本次运行还没走到这个采样步，就先退出采样循环，继续推进 LBM；走到了就调用 calNuRe()
            call calNuRe()
            if(outputSnapshotFile.EQ.1) then
                call output_SnapshotFile()          !每 0.5 t_ff 输出一次后处理 uvTrho 快照
            endif
        enddo
        if( (outputPltFile.EQ.1).AND.(MOD(restartItcOffset+itc, outputPltFileIntervalItc).EQ.0) ) then
            call output_Tecplot()  !非稳态模式下的可选周期 Tecplot 输出
        endif
        if( (outputReloadFile.EQ.1).AND.(MOD(restartItcOffset+itc, reloadFileIntervalItc).EQ.0) ) then
            call output_ReloadFile()      !非稳态模式下的可选周期 f/g 重启文件输出
        endif
#endif
     enddo

    call MPI_BARRIER(COMM2D, IERR)
    call CPU_TIME(timeEnd)         !当前进程累计消耗的 CPU 时间,包括并行
    timeEnd2 = MPI_WTIME()         !MPI 墙钟时间

#ifdef steadyFlow
    call output_Tecplot()          !输出最后一步的plt结果
    call output_SnapshotFile()              !输出最后一步的uvTrho数据
#endif

#ifdef unsteadyFlow
    if(isRoot) call output_unsteady_NuRe_postprocess()
#endif

    !===============================================================================================



    !===============================================================================================

    
#ifdef steadyFlow
    call calNuRe()
    call gather_output_fields_mpi()   ! 汇总各 rank 的局部 u/v/T/rho 到 root 的全局数组 u_all/v_all/T_all/rho_all，供稳态最终诊断使用
    call exit_data_2d_openacc()        ! 稳态最终诊断在 host 端 root 上直接使用 u_all/v_all/T_all/rho_all；这里先释放设备端局部映射
    if(isRoot) then
! 稳态最终标量诊断：
! 1) 只在 steadyFlow 收敛后调用；非稳态统计不要直接套用这一组最终诊断。
! 2) SideHeatedCell 的主热流方向为 x，用 u 和 dT/dx；RayleighBenardCell 的主热流方向为 y，用 v 和 dT/dy。
! 3) 壁面 Nu 和角点扩展默认采用半步长边界：流体节点距离物理边界 dx/2 或 dy/2。
! 4) Nu 极值、中心线速度极值使用五点最小二乘抛物线插值；中心线在偶数网格时用两侧流体节点线性插值。
! 5) 如果后续改成周期速度/温度边界，或改变半步长边界布置，这里和各后处理子程序都需要重新检查。
#ifdef SideHeatedCell                        
        call SideHeatedcalc_Nu_global()    ! x方向全场平均Nu
        call SideHeatedcalc_Nu_wall_avg()  ! 左/右壁, x中线平均Nu, 热壁Numax/Numin及位置
    
        call SideHeatedcalc_umid_max()     ! x中线上的u最大值及位置
        call SideHeatedcalc_vmid_max()     ! y中线上的v最大值及位置
#endif

#ifdef RayleighBenardCell
        call RBcalc_Nu_global()            ! y方向全场平均Nu
        call RBcalc_Nu_wall_avg()          ! 下/上壁, y中线平均Nu, 热壁Numax/Numin及位置
    
        call RBcalc_umid_max()             ! x中线上的u最大值及位置
        call RBcalc_vmid_max()             ! y中线上的v最大值及位置
#endif

#ifdef VerticalWallsNoslip
        ! psi/vort 后处理默认封闭腔体：四周无滑移，psi 在物理边界取同一常数。
        ! 若垂直边界改为周期速度边界，流函数边界补点和涡量单边差分需要另写周期版本。
        call calc_psi_vort_and_output()  ! 输出中心abs(psi), max(abs(psi))及位置；max位置用细网格样条插值
#endif
    endif
#endif





     

    open(unit=00,file=trim(settingsFile),status='unknown',position='append')        !在这个txt文件后面继续写（追加模式）
    write(00,*) "======================================================================"
    write(00,*) "Time (CPU) = ", real(timeEnd-timeStart,kind=8), "s"                             !当前进程累计消耗的 CPU 时间,包括并行
    write(00,*) "MLUPS = ", real( dble(nx)*dble(ny)*dble(itc)/(timeEnd-timeStart)/1.0d6,kind=8 )   !百万格点更新/秒
    write(00,*) "Time (MPI wall) = ", real(timeEnd2-timeStart2,kind=8), "s"                      !墙钟时间
    write(00,*) "MLUPS (MPI wall) = ", real( dble(nx)*dble(ny)*dble(itc)/(timeEnd2-timeStart2)/1.0d6,kind=8 )   !百万格点更新/秒
#ifdef steadyFlow
    write(00,'(a,1x,ES24.16E3)') "Nu_global =", real(Nu_global,kind=8)
    write(00,*) "Nu_hot    =", Nu_hot
    write(00,*) "Nu_cold   =", Nu_cold
#endif
    write(00,*) "useG =", useG



    write(00,*) "Dellocate Array......"
    call exit_data_2d_openacc()
    call deallocate_halo_buffers_2d_openacc_mpi()
    if(allocated(streamInterpIndexX)) then
        deallocate(streamInterpIndexX, streamInterpIndexY, streamInterpWeightX, streamInterpWeightY, &
          streamInterpValidX, streamInterpValidY)
    endif
    deallocate(f)
    deallocate(g)
    deallocate(f_post)
    deallocate(g_post)
    deallocate(u)
    deallocate(v)
    deallocate(T)
#ifdef steadyFlow
    deallocate(up)
    deallocate(vp)
    deallocate(Tp)
#endif
    deallocate(rho)
    deallocate(Fx)
    deallocate(Fy)
    deallocate(Bx_prev, By_prev)
    deallocate(u_all, v_all, T_all, rho_all)
    deallocate(f_all, g_all)
    deallocate(Bx_prev_all, By_prev_all)
    deallocate(XLocalCountAll, YLocalCountAll, XStartGlobalAll, YStartGlobalAll)
    write(00,*) "    "
    
    write(00,*) "Successfully: DNS completed!"
    
    string = ctime( time() )
    write(00,*) 'End:   ', string           !什么时候算完
    close(00)

    call MPI_FINALIZE(IERR)

    
    end program main

!===========================================================================================================================


!===========================================================================================================================
!===========================================================================================================================

!===========================================================================================================================
! 子程序: init_mpi_cartesian
! 作用: 初始化 MPI 进程网格、局部子域范围和邻居 rank。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===========================================================================================================================
  subroutine init_mpi_cartesian()
    use commondata
    implicit none
    integer(kind=4) :: px, py
    integer(kind=4) :: xBase, yBase, xRemainder, yRemainder
    integer(kind=4) :: xLocalMin, xLocalMax, yLocalMin, yLocalMax
    real(kind=8) :: aspectScore, bestAspectScore
    real(kind=8) :: commScore, bestCommScore
    real(kind=8), parameter :: scoreTolerance = 1.0d-12

    ! 从启动命令 mpiexec -n 读取实际 MPI 进程数；NPROC 后面用来自动分配 x/y 方向进程网格。
    call MPI_COMM_SIZE(MPI_COMM_WORLD, NPROC, IERR)
    ! 读取当前进程在 MPI_COMM_WORLD 中的编号；创建笛卡尔通信器后 MYID 可能会被重排更新。
    call MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)

#ifdef AccGpuSingleNode
    ! 单节点测试：NodeCount 应为 1；一 GPU 一 rank 时，总 MPI rank 数必须等于本节点分配的 GPU 数。
    if(NodeCount.NE.1) then
        if(MYID.EQ.0) then
            write(*,*) "Error: NodeCount must be 1 in AccGpuSingleNode mode."
            write(*,*) "NodeCount =", NodeCount
        endif
        call MPI_ABORT(MPI_COMM_WORLD, 1003, IERR)
    endif
    if(NPROC.NE.GPUS_PER_NODE) then
        if(MYID.EQ.0) then
            write(*,*) "Error: one-rank-per-GPU single-node mode requires NPROC = GPUS_PER_NODE."
            write(*,*) "NPROC =", NPROC, "; GPUS_PER_NODE =", GPUS_PER_NODE
        endif
        call MPI_ABORT(MPI_COMM_WORLD, 1003, IERR)
    endif
    NPROC_NODE = GPUS_PER_NODE
#endif
#ifdef AccGpuMultiNode
    ! 多节点集群：NodeCount 是本算例使用的节点数，GPUS_PER_NODE 是每个节点分配给本程序的 GPU 数。
    ! 当前只支持一 GPU 一 MPI rank，所以 mpiexec -n 得到的 NPROC 必须等于 NodeCount*GPUS_PER_NODE。
    if(NPROC.NE.NodeCount*GPUS_PER_NODE) then
        if(MYID.EQ.0) then
            write(*,*) "Error: one-rank-per-GPU multi-node mode requires NPROC = NodeCount*GPUS_PER_NODE."
            write(*,*) "NPROC =", NPROC, "; NodeCount =", NodeCount, "; GPUS_PER_NODE =", GPUS_PER_NODE
        endif
        call MPI_ABORT(MPI_COMM_WORLD, 1003, IERR)
    endif
    NPROC_NODE = GPUS_PER_NODE
#endif

    ! MPI 派生状态统一在初始化子程序里给默认值；真正的拓扑、边界和局部范围在下面覆盖。
    MY_COORD = 0
    left = MPI_PROC_NULL
    right = MPI_PROC_NULL
    down = MPI_PROC_NULL
    top = MPI_PROC_NULL
    isRoot = .false.
    hasLeftBoundary = .false.
    hasRightBoundary = .false.
    hasBottomBoundary = .false.
    hasTopBoundary = .false.
    xLocalCount = 0
    yLocalCount = 0
    xStartGlobal = 0
    xEndGlobal = 0
    yStartGlobal = 0
    yEndGlobal = 0

    dims = 0
#ifdef MpiDecomp1D
    ! 一维 MPI 分解：x 方向不切分，所有进程只沿 y 方向排成 1 x NPROC 的条带。
    ! 这种模式 halo 只需要上下交换
    dims(1) = 1
    dims(2) = NPROC
#endif
#ifdef MpiDecomp2D
    ! 二维 MPI 分解：遍历所有 px*py=NPROC 的候选。
    bestAspectScore = huge(1.0d0)
    bestCommScore = huge(1.0d0)
    do px = 1, NPROC
        ! px 是 x 方向的进程数；只有能整除 NPROC 时，py=NPROC/px 才是整数候选。
        if(mod(NPROC, px).NE.0) cycle
        py = NPROC / px
        ! 每个方向的进程数不能超过该方向的网格数，否则会出现空的局部块。
        if((px.GT.nx).OR.(py.GT.ny)) cycle

        ! 下面估算 nx 按 px 切分后各 rank 可能拿到的最小/最大 x 局部长度。
        ! 使用“余数前置”：前 xRemainder 个块多 1 个点。
        xBase = nx / px
        xRemainder = nx - xBase*px
        xLocalMin = xBase
        xLocalMax = xBase
        if(xRemainder.GT.0) xLocalMax = xBase + 1

        ! y 方向同理，得到候选分解下最小/最大的 y 局部长度。
        yBase = ny / py
        yRemainder = ny - yBase*py
        yLocalMin = yBase
        yLocalMax = yBase
        if(yRemainder.GT.0) yLocalMax = yBase + 1

        ! aspectScore 衡量最差局部块的长宽比偏离 1 的程度。
        ! 用 log(ratio) 是为了让 2:1 和 1:2 得到相同惩罚；
        ! 取两个极端组合，是为了覆盖所有 rank 里可能出现的最坏情况。
        aspectScore = max(abs(log(dble(xLocalMin)/dble(yLocalMax))), &
            abs(log(dble(xLocalMax)/dble(yLocalMin))))
        ! commScore 估算所有 rank-pair 内部 halo 边界的总长度。
        ! 对一条 x 分界线来说，单个 rank 发送的是 yLocalMin/yLocalMax 这样的局部段；
        ! 但这一整条分界线被 py 个 y 向子块拼起来，总长度仍然是全局 ny。
        ! y 分界线同理，总长度是全局 nx。
        ! 它只在长宽比评分几乎相同的时候作为并列判据，倾向于更少总通信边界。
        commScore = dble(px-1)*dble(ny) + dble(py-1)*dble(nx)

        ! 选择规则：
        ! 1) 优先让最差局部块更接近正方形；
        ! 2) 长宽比几乎相同时，选择内部通信边界更短的分解；
        ! 3) 仍然并列时选择更小的 px，保证同样输入下结果可复现。
        if((aspectScore.LT.bestAspectScore-scoreTolerance).OR. &
            ((abs(aspectScore-bestAspectScore).LE.scoreTolerance).AND. &
            ((commScore.LT.bestCommScore-scoreTolerance).OR. &
            ((abs(commScore-bestCommScore).LE.scoreTolerance).AND.(px.LT.dims(1)))))) then
            dims(1) = px
            dims(2) = py
            bestAspectScore = aspectScore
            bestCommScore = commScore
        endif
    enddo
#endif

    if((dims(1).LT.1).OR.(dims(2).LT.1)) then
        if(MYID.EQ.0) then
            write(*,*) "Error: no valid MPI Cartesian decomposition was found."
            write(*,*) "NPROC =", NPROC, "; mesh =", nx, ny
        endif
        call MPI_ABORT(MPI_COMM_WORLD, 1001, IERR)
    endif

    if((dims(1).GT.nx).OR.(dims(2).GT.ny)) then
        if(MYID.EQ.0) then
            write(*,*) "Error: too many MPI ranks for the grid."
            write(*,*) "dims =", dims(1), dims(2), "; mesh =", nx, ny
        endif
        ! MPI 程序中只 stop 某一个 rank 容易让其它 rank 卡在通信里；
        ! MPI_ABORT 会让 MPI_COMM_WORLD 中所有进程一起退出，1002 是这里自定义的错误码。
        call MPI_ABORT(MPI_COMM_WORLD, 1002, IERR)
    endif

    ! 按 dims(1) x dims(2) 建立二维笛卡尔通信器。
    ! 参数含义：
    !   MPI_COMM_WORLD : 原始全局通信器，包含 mpiexec 启动的所有 rank。
    !   2              : 笛卡尔拓扑维度数；这里是二维 x/y 进程网格。
    !   dims           : 每个方向上的进程数，dims(1) 对应 x，dims(2) 对应 y。
    !   periods        : 每个方向是否周期；periods(1) 对应 x，periods(2) 对应 y。
    !   .true.         : 允许 MPI 为拓扑通信重排 rank 编号。
    !   COMM2D         : 输出的新二维笛卡尔通信器。
    !   IERR           : MPI 错误码。
    call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dims, periods, .true., COMM2D, IERR)
    ! 进入 COMM2D 后重新读取本进程编号；若发生 rank 重排，MYID 会以 COMM2D 为准。
    call MPI_COMM_RANK(COMM2D, MYID, IERR)
    ! 把当前 rank 编号转换为笛卡尔坐标 MY_COORD=(x方向进程坐标, y方向进程坐标)。
    ! 参数含义：
    !   COMM2D   : 已创建好的二维笛卡尔通信器。
    !   MYID     : 当前进程在 COMM2D 中的 rank 编号。
    !   2        : 坐标维度数；这里返回二维坐标。
    !   MY_COORD : 输出坐标数组，MY_COORD(1) 对应 x 方向，MY_COORD(2) 对应 y 方向。
    !   IERR     : MPI 错误码。
    call MPI_CART_COORDS(COMM2D, MYID, 2, MY_COORD, IERR)
    ! 第 0 维对应 x 方向：找到当前子域左、右相邻 rank；非周期外边界会返回 MPI_PROC_NULL。
    ! 参数含义：
    !   COMM2D : 二维笛卡尔通信器。
    !   0      : 查询第 0 个拓扑维度，也就是 x 方向。
    !   1      : 位移步长，表示只查相邻一层进程。
    !   left   : 输出负方向邻居 rank；x 方向负方向就是左邻居。
    !   right  : 输出正方向邻居 rank；x 方向正方向就是右邻居。
    !   IERR   : MPI 错误码。
    call MPI_CART_SHIFT(COMM2D, 0, 1, left, right, IERR)
    ! 第 1 维对应 y 方向：找到当前子域下、上相邻 rank；非周期外边界会返回 MPI_PROC_NULL。
    call MPI_CART_SHIFT(COMM2D, 1, 1, down, top, IERR)

    ! 根据当前 rank 的笛卡尔坐标，计算它在全局网格中的拥有区间：
    !   xStartGlobal/yStartGlobal 是当前 rank 局部块在全局网格中的起始编号；
    !   xLocalCount/yLocalCount 是当前 rank 在 x/y 方向实际拥有的网格数量。
    call mesh_StartGlobal_LocalCount(nx, dims(1), MY_COORD(1), xStartGlobal, xLocalCount)
    call mesh_StartGlobal_LocalCount(ny, dims(2), MY_COORD(2), yStartGlobal, yLocalCount)
    ! 由起点和局部数量推出当前 rank 的全局终止编号。
    xEndGlobal = xStartGlobal + xLocalCount - 1
    yEndGlobal = yStartGlobal + yLocalCount - 1

    ! 只把全局外侧 rank 标记为物理边界；内部分区边界只做 halo 通信，不做壁面边界条件。
    hasLeftBoundary = (.not.periods(1)).AND.(MY_COORD(1).EQ.0)           ! x 非周期且位于最左侧进程列
    hasRightBoundary = (.not.periods(1)).AND.(MY_COORD(1).EQ.dims(1)-1) ! x 非周期且位于最右侧进程列
    hasBottomBoundary = MY_COORD(2).EQ.0                                ! 位于最下侧进程行
    hasTopBoundary = MY_COORD(2).EQ.dims(2)-1                           ! 位于最上侧进程行
    isRoot = MYID.EQ.0                                                  ! COMM2D 中的 0 号 rank 负责输出和汇总

    allocate(XLocalCountAll(0:NPROC-1), YLocalCountAll(0:NPROC-1))
    allocate(XStartGlobalAll(0:NPROC-1), YStartGlobalAll(0:NPROC-1))
    ! 汇总所有 rank 的局部块尺寸和全局起点；root 后续按这些信息把局部场拼回全局数组。
    ! MPI_ALLGATHER 参数含义：
    !   第 1 个变量       : 当前 rank 要发送的本地标量，例如 xLocalCount。
    !   第 2/3 个参数     : 当前 rank 发送 1 个 MPI_INTEGER。
    !   第 4 个变量       : 接收数组，通信完成后保存所有 rank 发来的对应标量。
    !   第 5/6 个参数     : 每个 rank 接收 1 个 MPI_INTEGER。
    !   COMM2D/IERR       : 笛卡尔通信器和 MPI 错误码。
    call MPI_ALLGATHER(xLocalCount, 1, MPI_INTEGER, XLocalCountAll, 1, MPI_INTEGER, COMM2D, IERR)     ! 每个 rank 的 x 方向局部数量
    call MPI_ALLGATHER(yLocalCount, 1, MPI_INTEGER, YLocalCountAll, 1, MPI_INTEGER, COMM2D, IERR)     ! 每个 rank 的 y 方向局部数量
    call MPI_ALLGATHER(xStartGlobal, 1, MPI_INTEGER, XStartGlobalAll, 1, MPI_INTEGER, COMM2D, IERR)      ! 每个 rank 的全局 x 起点
    call MPI_ALLGATHER(yStartGlobal, 1, MPI_INTEGER, YStartGlobalAll, 1, MPI_INTEGER, COMM2D, IERR)      ! 每个 rank 的全局 y 起点
    
    !root rank 继续用默认的 settingsFile，
    !非 root rank 改用带 rank 编号的日志文件，避免多个 rank 同时写同一个日志文件。
    if(.not.isRoot) then
        write(settingsFile,'("SimulationSettings2DOpenaccMpi-rank",I6.6,".txt")') MYID
    endif

    if(loadInitField.EQ.1) then
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
    else
        open(unit=00,file=trim(settingsFile),status='replace')
    endif
    write(00,*) "MPI Cartesian dims from init_mpi_cartesian (x,y):", dims(1), dims(2)
    close(00)

    return
  end subroutine init_mpi_cartesian
!===================================================================================================
! init_mpi_cartesian 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===========================================================================================================================
! 子程序: mesh_StartGlobal_LocalCount
! 作用: 对一个方向做余数前置的负载均衡切分。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===========================================================================================================================
  subroutine mesh_StartGlobal_LocalCount(TotalMeshCount, MpiCount, MpiCoord, StartGlobal, LocalCount)
    implicit none
    integer(kind=4), intent(in) :: TotalMeshCount      ! 该方向全局网格数，例如 nx 或 ny
    integer(kind=4), intent(in) :: MpiCount            ! 该方向 MPI 进程数，例如 dims(1) 或 dims(2)
    integer(kind=4), intent(in) :: MpiCoord            ! 当前 rank 在该方向上的进程坐标，例如 MY_COORD(1) 或 MY_COORD(2)
    integer(kind=4), intent(out) :: StartGlobal        ! 当前 rank 在该方向上的全局起始编号
    integer(kind=4), intent(out) :: LocalCount         ! 当前 rank 在该方向上分到的局部网格数量
    integer(kind=4) :: BaseMeshCount, RemainderMeshCount

    BaseMeshCount = TotalMeshCount / MpiCount
    RemainderMeshCount = TotalMeshCount - BaseMeshCount*MpiCount
    if(MpiCoord.LT.RemainderMeshCount) then
        LocalCount = BaseMeshCount + 1
        StartGlobal = 1 + MpiCoord*(BaseMeshCount + 1)
    else
        LocalCount = BaseMeshCount
        StartGlobal = 1 + RemainderMeshCount*(BaseMeshCount + 1) + &
            (MpiCoord-RemainderMeshCount)*BaseMeshCount
    endif

    return
  end subroutine mesh_StartGlobal_LocalCount
!===================================================================================================
! mesh_StartGlobal_LocalCount 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===========================================================================================================================
!===================================================================================================
! 子程序: initial
! 作用: 初始化网格坐标、场变量、分布函数、输出文件和重启信息。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===========================================================================================================================
  subroutine initial()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: globalI, globalJ
    integer(kind=4) :: alpha
    real(kind=8) :: un(0:8)
    real(kind=8) :: us2, Bx, By
    real(kind=8) :: xLen, yLen, rbInitPerturbAmp
    character(len=100) :: reloadFileName
    

    itc = 0
    errorU = 100.0d0
    errorT = 100.0d0 
    snapshotFileNum = 0
    pltFileNum = 0
    restartItcOffset = 0
    reloadMetadataLoaded = .false.
    outputSnapshotIntervalItc = max(1, int(outputSnapshotInterval*timeUnit+0.5d0))
    reloadFileIntervalItc = max(1, int(reloadFileInterval*timeUnit+0.5d0))
    outputPltFileIntervalItc = max(1, int(outputPltFileInterval*timeUnit+0.5d0))


    !-----------------------------------------------------------------------------------------------
    !记录各种信息在日志文件中
    open(unit=00,file=trim(settingsFile),status='unknown',position='append')  !在这个txt文件后面继续写（追加模式）
    
    if((outputSnapshotFile.EQ.1).AND.isRoot) then
        open(unit=01,file=trim(snapshotFilePrefix)//"-"//"readme",status="unknown")    !trim去掉字符串尾部空格，换了存储路径，可自己更改
        write(01,*) "snapshot file prefix exists!"
        write(01,*) "records: U_nd, V_nd, T, rho on ISLBM nonuniform nodes"
        write(01,*) "coordinates are xp(1:nx), yp(1:ny) in ", trim(snapshotFilePrefix)//"-mesh.dat"
        close(01)
        write(00,*) "Snapshot data will be stored in ", snapshotFilePrefix
    endif
    if((outputPltFile.EQ.1).AND.isRoot) then
        open(unit=01,file=trim(pltFolderPrefix)//"-"//"readme",status="unknown")     !读取路径pltFolderPrefix="../pltFile/buoyancyCavity
        write(01,*) "pltFile folder exist!"
        close(01)
        write(00,*) "Data will be stored in ", pltFolderPrefix
    endif
    if((outputReloadFile.EQ.1).AND.isRoot) then
        open(unit=01,file=trim(reloadFilePrefix)//"-"//"readme",status="unknown")
        write(01,*) "reloadFile prefix exists!"
        close(01)
        write(00,*) "Reload data will be stored in ", reloadFilePrefix
    endif
    
#ifdef EnableLegacyThermalScheme
    if( (paraA.GE.1.0d0).OR.(paraA.LE.-4.0d0) ) then                           !只有在[-4,1]才可以，要不然预警退出
        write(00,*) "----------------------------------"
        write(00,*) "paraA=", paraA
        write(00,*) "Error: condition not meet for the legacy thermal algorithm"
        write(00,*) "Ref: Luo2013, CMA"
        write(00,*) "Please try to reduce Mach number"
        write(00,*) "----------------------------------"
        stop
    endif
#endif

    write(00,*)"-------------------------------------------------------------------------------"
    write(00,*) 'Global mesh:',nx,ny
    write(00,*) 'Local mesh/range:', xLocalCount, yLocalCount, xStartGlobal, xEndGlobal, yStartGlobal, yEndGlobal
    write(00,*) 'Rayleigh=',real(Rayleigh,kind=8), '; Prandtl =',real(Prandtl,kind=8), '; Mach =',real(Mach,kind=8)
    write(00,*) "Length unit: L0 =", real(lengthUnit,kind=8)
    write(00,*) "Time unit: Sqrt(L0/(gBeta*DeltaT)) =", real(timeUnit,kind=8)
    write(00,*) "Velocity unit: Sqrt(gBeta*L0*DeltaT) =", real(velocityUnit,kind=8)
    write(00,*) "   "
    write(00,*) 'tauf=',real(tauf,kind=8)
#ifdef EnableLegacyThermalScheme
    write(00,*) "thermalScheme = legacy D2Q5 (Qk/Qnu)"
    write(00,*) 'Qk=',real(Qk,kind=8), '; Qnu=',real(Qnu,kind=8), '; paraA=',real(paraA,kind=8)
#else
    write(00,*) "thermalScheme = current D2Q5 (s_j/s_e/s_q)"
    write(00,*) 'taug=',real(taug,kind=8)
#endif
    write(00,*) "viscosity =",real(viscosity,kind=8), "; diffusivity =",real(diffusivity,kind=8)
    write(00,*) "outputSnapshotFile =", outputSnapshotFile
    write(00,*) "outputSnapshotInterval =", real(outputSnapshotInterval,kind=8), "free-fall time units"
    write(00,*) "outputSnapshotIntervalItc =", outputSnapshotIntervalItc, "in itc units"
    write(00,*) "outputPltFile =", outputPltFile
    write(00,*) "outputPltFileInterval =", real(outputPltFileInterval,kind=8), "free-fall time units"
    write(00,*) "outputPltFileIntervalItc =", outputPltFileIntervalItc, "in itc units"
    write(00,*) "outputReloadFile =", outputReloadFile
    write(00,*) "reloadFileInterval =", real(reloadFileInterval,kind=8), "free-fall time units"
    write(00,*) "reloadFileIntervalItc =", reloadFileIntervalItc, "in itc units"
#ifdef unsteadyFlow
    write(00,*) "unsteadyRunDuration =", real(unsteadyRunDuration,kind=8), "free-fall time units"
    write(00,*) "unsteadySampleCount =", unsteadySampleCount     !输出次数计数器
#endif
    if(loadInitField.EQ.1) then
        write(00,*) "Restart offsets will be read from reload metadata when available."
    endif
    write(00,*) "itc_max =",itc_max
    write(00,*) "default epsU =", real(epsU,kind=8),"; epsT =", real(epsT,kind=8)
    write(00,*) "useG =", useG
    write(00,*) "    "

#ifdef RayleighBenardCell
    write(00,*) "I am Rayleigh Benard Cell"
#endif
#ifdef SideHeatedCell
    write(00,*) "I am Side Heated Cell"
#endif
    
#ifdef steadyFlow
    write(00,*) "I am steadyFlow"
#endif
#ifdef unsteadyFlow
    write(00,*) "I am unsteadyFlow"
#endif
    !-----------------------------------------------------------------------------------------------



    !-----------------------------------------------------------------------------------------------
    ! ISLBM 节点坐标与积分权重。MPI 版本中 xp/yp/quadWidthX/quadWidthY 仍然保留为全局一维表，
    ! 每个 rank 都有完整副本；后续局部 streaming 模板再用 xStartGlobal/yStartGlobal 映射。
    call build_islbm_mesh()
    call build_islbm_quadrature()
    call allocate_islbm_streaming_stencils_mpi()
    call build_islbm_streaming_stencils_mpi()
    if((outputSnapshotFile.EQ.1).AND.isRoot) then
        call output_SnapshotMeshFile()
        write(00,*) "Snapshot mesh coordinates stored in ", trim(snapshotFilePrefix)//"-mesh.dat"
    endif
    if(isRoot) then
    write(00,*) "ISLBM mesh = erf; stretchA =", real(ISLBM_StretchA,kind=8)
        write(00,*) "ISLBM effective lengthUnit L0 =", real(lengthUnit,kind=8)
        write(00,*) "ISLBM lattice unit in normalized coordinates =", real(ISLBM_LatticeUnit,kind=8)
        write(00,*) "ISLBM quadrature sums:", real(quadSumX,kind=8), real(quadSumY,kind=8), real(quadSumArea,kind=8)
        write(00,*) "ISLBM MPI streaming fallback stencil entries:", streamStencilFallbackGlobal
        write(00,*) "Internal rank interfaces use one-layer halo values for ISLBM Lagrange streaming."
        write(00,*) "Fallback entries are physical-boundary/out-of-halo cases handled by boundary conditions."
    endif

    allocate (u(xLocalCount,yLocalCount))
    allocate (v(xLocalCount,yLocalCount))
    allocate (T(xLocalCount,yLocalCount))
    allocate (rho(xLocalCount,yLocalCount))

#ifdef steadyFlow
    allocate (up(xLocalCount,yLocalCount))
    allocate (vp(xLocalCount,yLocalCount))
    allocate (Tp(xLocalCount,yLocalCount))
#endif

    ! OpenACC 版采用 (x, y, alpha) 存储分布函数，使 GPU 上相邻线程沿 x 方向访问时更连续。
    allocate (f(xLocalCount,yLocalCount,0:8))
    allocate (f_post(0:xLocalCount+1,0:yLocalCount+1,0:8))
    allocate (g(xLocalCount,yLocalCount,0:4))
    allocate (g_post(0:xLocalCount+1,0:yLocalCount+1,0:4))

    allocate (Fx(xLocalCount,yLocalCount))
    allocate (Fy(xLocalCount,yLocalCount))

    allocate (Bx_prev(xLocalCount,yLocalCount), By_prev(xLocalCount,yLocalCount)) 
    ! 存储全场数据
    allocate (u_all(nx,ny), v_all(nx,ny), T_all(nx,ny), rho_all(nx,ny))
    allocate (f_all(nx,ny,0:8), g_all(nx,ny,0:4))
    allocate (Bx_prev_all(nx,ny), By_prev_all(nx,ny))

    !-----------------------------------------------------------------------------------------------



    !-----------------------------------------------------------------------------------------------
    !初始化
    rho = 1.0d0                                     !密度rho=1
    
    omega(0) = 4.0d0/9.0d0
    do alpha=1,4
        omega(alpha) = 1.0d0/9.0d0
    enddo
    do alpha=5,8
        omega(alpha) = 1.0d0/36.0d0
    enddo
    
#ifdef EnableLegacyThermalScheme
    omegaT(0) = (1.0d0-paraA)/5.0d0
    do alpha=1,4
        omegaT(alpha) = (paraA+4.0d0)/20.0d0
    enddo
#else
    omegaT(0) = 1.0d0/3.0d0
    do alpha=1,4
        omegaT(alpha) = 1.0d0/6.0d0
    enddo
#endif

    if(loadInitField.EQ.0) then                    !在不加载文件的情况下，都是零场为初值
    
        u = 0.0d0
        v = 0.0d0
        T = 0.0d0
        Bx= 0.0d0
        By= 0.0d0
        Bx_prev= 0.0d0
        By_prev= 0.0d0
        
        write(00,*) "Initial field is set exactly"
        if(reloadDimensionlessTime.NE.0.0d0) then        !在不加载文件的情况下，reloadDimensionlessTime必须是零
            write(00,*) "Error: since loadInitField.EQ.0, reloadDimensionlessTime should also be 0"
            stop
        endif
        
#ifdef VerticalWallsNoslip
        write(00,*) "Velocity B.C. for vertical walls are: ===No-slip wall==="
#endif
#ifdef VerticalWallsPeriodicalU
        write(00,*) "Velocity B.C. for vertical walls are: ===Periodical==="
#endif
#ifdef HorizontalWallsNoslip
        write(00,*) "Velocity B.C. for horizontal walls are: ===No-slip wall==="
#endif

#ifdef VerticalWallsConstT
    do j = 1, yLocalCount                              !MPI 局部块：用全局 i 坐标保持原来的线性初始温度
        do i = 1, xLocalCount
            globalI = xStartGlobal + i - 1
            T(i,j) = Thot + (xp(globalI)-xp(0)) / (xp(nx+1)-xp(0)) * (Tcold-Thot)
        enddo
    enddo
    write(00,*) "Temperature B.C. for vertical walls are:===Hot/cold wall==="
#endif

#ifdef HorizontalWallsConstT
    do i = 1, xLocalCount
        do j = 1, yLocalCount
            globalJ = yStartGlobal + j - 1
            T(i,j) = Thot + (yp(globalJ)-yp(0)) / (yp(ny+1)-yp(0)) * (Tcold-Thot)
        enddo
    enddo
#ifdef RayleighBenardCell
    if (Rayleigh.LT.1.0d4) then
        xLen = xp(nx+1)
        yLen = yp(ny+1)
        rbInitPerturbAmp = 1.0d-3*(Thot-Tcold)
        do i = 1, xLocalCount
            do j = 1, yLocalCount
                globalI = xStartGlobal + i - 1
                globalJ = yStartGlobal + j - 1
                T(i,j) = T(i,j) + rbInitPerturbAmp * dsin(2.0d0*pi*xp(globalI)/xLen) * dsin(pi*yp(globalJ)/yLen)
            enddo
        enddo
        write(00,'(a,1x,es12.4)') "RB initial T perturbation amplitude =", rbInitPerturbAmp
    else
        write(00,*) "RB initial T perturbation skipped because Rayleigh > 1.0d4"
    endif
#endif
    write(00,*) "Temperature B.C. for horizontal walls are:===Hot/cold wall==="
#endif

#ifdef VerticalWallsAdiabatic
    write(00,*) "Temperature B.C. for vertical walls are:===Adiabatic wall==="
#endif

#ifdef HorizontalWallsAdiabatic
    write(00,*) "Temperature B.C. for horizontal walls are:===Adiabatic wall==="
#endif

        f = 0.0d0
        g = 0.0d0
        do j = 1,yLocalCount
            do i = 1,xLocalCount
                us2 = u(i,j)*u(i,j)+v(i,j)*v(i,j)
                do alpha = 0, 8
                    un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                    f(i,j,alpha) = rho(i,j)*omega(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)  !D2Q9标准feq
                enddo
                do alpha = 0, 4
                    un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha) 
                    g(i,j,alpha) = omegaT(alpha)*T(i,j)*(1.0d0+thermalGeqCoeff*un(alpha))
                enddo
            enddo
        enddo
#ifdef EnableUseG
        do j = 1,yLocalCount
            do i = 1,xLocalCount
              Bx = u(i,j) * T(i,j)
              By = v(i,j) * T(i,j)
              Bx_prev(i,j) = Bx
              By_prev(i,j) = By
            enddo
        enddo
#endif



    elseif(loadInitField.EQ.1) then                               !续算分支：从旧的严格重启文件恢复 f/g 和输出计数
    ! 正常断电续算时，先读取 <reloadFilePrefix>-latest.meta；
    ! meta 会告诉代码实际要读哪个 <reloadFilePrefix>-*.bin，以及旧算例已经累计到哪里。
    !这里可以让 reloadDimensionlessTime 等于0，反正最后都可以通过meta文件读取到新的 reloadDimensionlessTime
        write(00,*) "Read reload metadata before choosing the restart .bin file."
        write(reloadFileName,'(i12.12)') reloadFileNum             !latest .meta 缺失时才依赖这个手工编号
        reloadFileName = adjustl(reloadFileName)                  !adjustl把字符串左对齐，把前导空格移到字符串末尾
        !latest .meta 存在时：精确读取账本并更新 reloadFileName；不存在时保守推断。
        call read_reload_metadata(reloadFileName)
        write(00,*) "Load initial field from previous simulation: ", &
            trim(reloadFilePrefix), "-", trim(reloadFileName), ".bin"
        if(.not.reloadMetadataLoaded) then
            write(00,*) "WARNING: no reload metadata file found; restart offsets were inferred."
            write(00,*) "         For exact continuation, use reload files written after this patch."
        endif
        call read_reload_fields_mpi(reloadFileName)   ! root 读取全局 restart 文件，并把对应数据分发/同步到各 rank 的局部 f/g
        call reconstruct_macro_from_fg()              ! 根据读入的 f/g 重新计算 rho/u/v/T，恢复续算所需的宏观场
        write(00,*) "Raw data is loaded from the file: ", trim(reloadFilePrefix), "-", trim(reloadFileName), ".bin"
        write(00,*) "Restart offset itc =", restartItcOffset
        write(00,*) "Restart offset time_tf =", real(reloadDimensionlessTime,kind=8)
        write(00,*) "Continue output counters: snapshot/plt/reload =", snapshotFileNum, pltFileNum, reloadFileNum
    else
        write(00,*) "Error: initial field is not properly set"   !如果 loadInitField 不是 0/1 或逻辑不一致，直接停止
        stop
    endif
    
    write(00,*)"-------------------------------------------------------------------------------"
close(00)
    
    f_post = 0.0d0
    g_post = 0.0d0
        
    if(loadInitField.EQ.0) then
        snapshotFileNum = 0
        pltFileNum = 0
        reloadFileNum = 0
        restartItcOffset = 0
        reloadDimensionlessTime = 0.0d0
#ifdef steadyFlow
        ! 新算例第一段收敛误差应从初始场开始比较。
        up = u
        vp = v
        Tp = T
#endif
    else
#ifdef steadyFlow
        ! 重启后第一段收敛误差应从载入场继续比较。
        up = u
        vp = v
        Tp = T
#endif
    endif
    dimensionlessTime = 0
    ! 新算例：清零，开始记录新历史。
    ! 续算：也清零，但不是丢旧历史；旧历史在 .dat 文件里，新数组只记录本次续算段。
    ! 写出时间轴时会叠加 reloadDimensionlessTime，所以不会从 0 t_ff 重新编号。
    NuVolAvg = 0.0d0
    ReVolAvg = 0.0d0
    
    return
  end subroutine initial
!===================================================================================================
! initial 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
!===================================================================================================
!初始化结束
!===================================================================================================


!===================================================================================================
! 子程序: enter_data_2d_openacc
! 作用: 把局部拥有区数组放到 OpenACC 设备端，后续 GPU kernel 默认直接使用 device resident 数据。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine enter_data_2d_openacc()
    use openacc
    use commondata
    implicit none

    !$acc enter data copyin(xp,yp,quadWidthX,quadWidthY,quadSumX,quadSumY,quadSumArea,ex,ey,omega,omegaT)
    !$acc enter data copyin(streamInterpIndexX,streamInterpIndexY,streamInterpWeightX, &
    !$acc& streamInterpWeightY,streamInterpValidX,streamInterpValidY)
    !$acc enter data copyin(u,v,T,rho,f,g,Fx,Fy,Bx_prev,By_prev)
    !$acc enter data create(f_post,g_post)
#ifdef steadyFlow
    !$acc enter data copyin(up,vp,Tp)
#endif
    accDataResident = .true.
    return
  end subroutine enter_data_2d_openacc
!===================================================================================================
! enter_data_2d_openacc 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: update_host_snapshot_2d_openacc
! 作用: MPI gather 或 snapshot 输出前，把设备端局部 u/v/T/rho 同步回 host。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine update_host_snapshot_2d_openacc()
    use openacc
    use commondata
    implicit none

    !$acc wait(1)
    !$acc update self(u,v,T,rho)
    return
  end subroutine update_host_snapshot_2d_openacc
!===================================================================================================
! update_host_snapshot_2d_openacc 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: update_host_tecplot_2d_openacc
! 作用: Tecplot 输出前，把设备端局部 u/v/T 同步回 host。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine update_host_tecplot_2d_openacc()
    use openacc
    use commondata
    implicit none

    !$acc wait(1)
    !$acc update self(u,v,T)
    return
  end subroutine update_host_tecplot_2d_openacc
!===================================================================================================
! update_host_tecplot_2d_openacc 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: update_host_reload_2d_openacc
! 作用: restart 输出前，把设备端局部 f/g 和 UseG 历史量同步回 host。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine update_host_reload_2d_openacc()
    use openacc
    use commondata
    implicit none

    !$acc wait(1)
    !$acc update self(f,g)
#ifdef EnableUseG
    !$acc update self(Bx_prev,By_prev)
#endif
    return
  end subroutine update_host_reload_2d_openacc
!===================================================================================================
! update_host_reload_2d_openacc 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: allocate_halo_buffers_2d_openacc_mpi
! 作用: 按当前 MPI rank 的局部块尺寸分配 halo 通信缓冲区，并在 GPU 上创建对应设备端副本。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine allocate_halo_buffers_2d_openacc_mpi()
    use commondata
    implicit none

    ! 若 fHaloSendDown 已经 associated，说明 halo buffer 已分配过；
    ! 先释放旧的 host 数组和 GPU 设备端副本，再按当前局部尺寸重新分配。
    ! 若它仍是 =>null() 的空 pointer，则 associated(...) 为 .false.，直接执行下面的 allocate。
    if(associated(fHaloSendDown)) call deallocate_halo_buffers_2d_openacc_mpi()

    ! f: 每条边只交换 3 个跨界 D2Q9 方向；x 向列交换包含上下 halo，所以长度为 3*(yLocalCount+2)。
    allocate(fHaloSendDown(3*xLocalCount), fHaloSendUp(3*xLocalCount))
    allocate(fHaloRecvDown(3*xLocalCount), fHaloRecvUp(3*xLocalCount))
    allocate(fHaloSendLeft(3*(yLocalCount+2)), fHaloSendRight(3*(yLocalCount+2)))
    allocate(fHaloRecvLeft(3*(yLocalCount+2)), fHaloRecvRight(3*(yLocalCount+2)))

    ! g: 每条边只交换 1 个跨界 D2Q5 方向；D2Q5 没有对角方向，所以列交换不需要上下 halo。
    allocate(gHaloSendDown(xLocalCount), gHaloSendUp(xLocalCount))
    allocate(gHaloRecvDown(xLocalCount), gHaloRecvUp(xLocalCount))
    allocate(gHaloSendLeft(yLocalCount), gHaloSendRight(yLocalCount))
    allocate(gHaloRecvLeft(yLocalCount), gHaloRecvRight(yLocalCount))

    !$acc enter data create(...) 的作用：
    !   在 GPU 上为这些动态分配的 halo buffer 创建设备端内存，并建立 host pointer 到 device pointer 的映射。
    !   这里不会把 CPU 端数据 copyin 到 GPU；这些 buffer 后面会直接在 GPU kernel 中 pack/clear/unpack。
    !   有了这个映射，后续 present(...) 可以确认 buffer 已在设备端，
    !   host_data use_device(...) 也能把 GPU 地址交给支持 GPU-aware 的 MPI_SENDRECV。
    !$acc enter data create(fHaloSendDown,fHaloSendUp,fHaloRecvDown,fHaloRecvUp)
    !$acc enter data create(fHaloSendLeft,fHaloSendRight,fHaloRecvLeft,fHaloRecvRight)
    !$acc enter data create(gHaloSendDown,gHaloSendUp,gHaloRecvDown,gHaloRecvUp)
    !$acc enter data create(gHaloSendLeft,gHaloSendRight,gHaloRecvLeft,gHaloRecvRight)
    return
  end subroutine allocate_halo_buffers_2d_openacc_mpi
!===================================================================================================
! allocate_halo_buffers_2d_openacc_mpi 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: deallocate_halo_buffers_2d_openacc_mpi
! 作用: 释放 halo 通信缓冲区的 GPU 设备端副本和 host 端数组。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine deallocate_halo_buffers_2d_openacc_mpi()
    use commondata
    implicit none

    if(.not.associated(fHaloSendDown)) return
    !$acc exit data delete(...) 的作用：
    !   释放这些 halo buffer 在 GPU 上的设备端内存，并删除 OpenACC 维护的 host/device 映射。
    !   delete 不会把 GPU 数据 copyout 回 CPU；这些通信缓冲区只是临时工作区，不需要保留其内容。
    !   下面的 Fortran deallocate(...) 才负责释放 CPU/host 端 pointer 数组。
    !$acc exit data delete(fHaloSendDown,fHaloSendUp,fHaloRecvDown,fHaloRecvUp)
    !$acc exit data delete(fHaloSendLeft,fHaloSendRight,fHaloRecvLeft,fHaloRecvRight)
    !$acc exit data delete(gHaloSendDown,gHaloSendUp,gHaloRecvDown,gHaloRecvUp)
    !$acc exit data delete(gHaloSendLeft,gHaloSendRight,gHaloRecvLeft,gHaloRecvRight)
    deallocate(fHaloSendDown, fHaloSendUp, fHaloRecvDown, fHaloRecvUp)
    deallocate(fHaloSendLeft, fHaloSendRight, fHaloRecvLeft, fHaloRecvRight)
    deallocate(gHaloSendDown, gHaloSendUp, gHaloRecvDown, gHaloRecvUp)
    deallocate(gHaloSendLeft, gHaloSendRight, gHaloRecvLeft, gHaloRecvRight)
    ! deallocate(...) 释放 pointer 当前指向的 host 数组内存；
    ! nullify(...) 再把 pointer 本身重置为空指针状态，使 associated(...) 重新变为 .false.。
    ! 这样后续如果再次进入分配子程序，不会把已经释放过的 pointer 误判为仍然有效。
    nullify(fHaloSendDown, fHaloSendUp, fHaloRecvDown, fHaloRecvUp)
    nullify(fHaloSendLeft, fHaloSendRight, fHaloRecvLeft, fHaloRecvRight)
    nullify(gHaloSendDown, gHaloSendUp, gHaloRecvDown, gHaloRecvUp)
    nullify(gHaloSendLeft, gHaloSendRight, gHaloRecvLeft, gHaloRecvRight)
    return
  end subroutine deallocate_halo_buffers_2d_openacc_mpi
!===================================================================================================
! deallocate_halo_buffers_2d_openacc_mpi 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: build_islbm_mesh
! 作用: 构造论文 ISLBM 使用的 half-way erf 非均匀全局坐标。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine build_islbm_mesh()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: rawX(0:nx+1), rawY(0:ny+1)
    real(kind=8) :: erfNorm, leftWall, rightWall, bottomWall, topWall, lengthX, lengthY

    ! 第一步: 生成原始erf拉伸坐标rawX/rawY。此时坐标还没有按half-way物理壁面修正。
    erfNorm = erf(0.5d0*ISLBM_StretchA)
    do i = 0, nx+1
        rawX(i) = 0.5d0*(1.0d0 + erf(ISLBM_StretchA*(dble(i)/dble(nx+1)-0.5d0))/erfNorm)
    enddo
    do j = 0, ny+1
        rawY(j) = 0.5d0*(1.0d0 + erf(ISLBM_StretchA*(dble(j)/dble(ny+1)-0.5d0))/erfNorm)
    enddo

    ! 第二步: 采用half-way壁面。物理壁面位于参考端点与第一个内部流体节点之间。
    leftWall = 0.5d0*rawX(1)
    rightWall = 1.0d0 - 0.5d0*rawX(1)
    bottomWall = 0.5d0*rawY(1)
    topWall = 1.0d0 - 0.5d0*rawY(1)

    lengthX = rightWall - leftWall
    lengthY = topWall - bottomWall
    xp(0) = 0.0d0
    xp(nx+1) = 1.0d0
    do i = 1, nx
        xp(i) = (rawX(i)-leftWall)/lengthX
    enddo
    yp(0) = 0.0d0
    yp(ny+1) = 1.0d0
    do j = 1, ny
        yp(j) = (rawY(j)-bottomWall)/lengthY
    enddo

    return
  end subroutine build_islbm_mesh
!===================================================================================================
! build_islbm_mesh 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: build_islbm_quadrature
! 作用: 基于全局非均匀坐标构造 midpoint-rule 积分权重。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine build_islbm_quadrature()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: leftGhostX, rightGhostX, bottomGhostY, topGhostY

    leftGhostX = -xp(1)
    rightGhostX = 1.0d0 + (1.0d0 - xp(nx))
    bottomGhostY = -yp(1)
    topGhostY = 1.0d0 + (1.0d0 - yp(ny))

    quadWidthX(1) = 0.5d0*(xp(2)-leftGhostX)
    do i = 2, nx-1
        quadWidthX(i) = 0.5d0*(xp(i+1)-xp(i-1))
    enddo
    quadWidthX(nx) = 0.5d0*(rightGhostX-xp(nx-1))

    quadWidthY(1) = 0.5d0*(yp(2)-bottomGhostY)
    do j = 2, ny-1
        quadWidthY(j) = 0.5d0*(yp(j+1)-yp(j-1))
    enddo
    quadWidthY(ny) = 0.5d0*(topGhostY-yp(ny-1))

    quadSumX = 0.0d0
    do i = 1, nx
        quadSumX = quadSumX + quadWidthX(i)
    enddo
    quadSumY = 0.0d0
    do j = 1, ny
        quadSumY = quadSumY + quadWidthY(j)
    enddo
    quadSumArea = quadSumX*quadSumY

    return
  end subroutine build_islbm_quadrature
!===================================================================================================
! build_islbm_quadrature 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: allocate_islbm_streaming_stencils_mpi
! 作用: 按当前 rank 的局部块尺寸分配二次 Lagrange streaming 模板。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine allocate_islbm_streaming_stencils_mpi()
    use commondata
    implicit none

    if(allocated(streamInterpIndexX)) then
        deallocate(streamInterpIndexX, streamInterpIndexY, streamInterpWeightX, streamInterpWeightY, &
          streamInterpValidX, streamInterpValidY)
    endif
    allocate(streamInterpIndexX(0:8,xLocalCount,3), streamInterpIndexY(0:8,yLocalCount,3))
    allocate(streamInterpWeightX(0:8,xLocalCount,3), streamInterpWeightY(0:8,yLocalCount,3))
    allocate(streamInterpValidX(0:8,xLocalCount), streamInterpValidY(0:8,yLocalCount))

    return
  end subroutine allocate_islbm_streaming_stencils_mpi
!===================================================================================================
! allocate_islbm_streaming_stencils_mpi 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: build_islbm_streaming_stencils_mpi
! 作用: 将全局三点模板映射到 rank-local/halo 索引；内部rank界面允许使用一层halo。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine build_islbm_streaming_stencils_mpi()
    use commondata
    implicit none
    integer(kind=4) :: alpha, i, j, globalI, globalJ
    integer(kind=4) :: idxGlobal(3), idxLocal(3)
    real(kind=8) :: w(3), target
    logical :: ok

    streamInterpIndexX = 1
    streamInterpIndexY = 1
    streamInterpWeightX = 0.0d0
    streamInterpWeightY = 0.0d0
    streamInterpValidX = .false.
    streamInterpValidY = .false.
    streamStencilFallbackLocal = 0
    streamStencilFallbackGlobal = 0

    do alpha = 0, 8
        do i = 1, xLocalCount
            globalI = xStartGlobal + i - 1
            target = xp(globalI) - dble(ex(alpha))*ISLBM_LatticeUnit
            call build_streaming_stencil_1d(nx, xp(1:nx), globalI, target, idxGlobal, w, ok)
            idxLocal = idxGlobal - xStartGlobal + 1
            if(ok.AND.all(idxLocal.GE.0).AND.all(idxLocal.LE.xLocalCount+1)) then
                streamInterpValidX(alpha,i) = .true.
                streamInterpIndexX(alpha,i,:) = idxLocal
                streamInterpWeightX(alpha,i,:) = w
            elseif(ok) then
                streamStencilFallbackLocal = streamStencilFallbackLocal + 1
            endif
        enddo
        do j = 1, yLocalCount
            globalJ = yStartGlobal + j - 1
            target = yp(globalJ) - dble(ey(alpha))*ISLBM_LatticeUnit
            call build_streaming_stencil_1d(ny, yp(1:ny), globalJ, target, idxGlobal, w, ok)
            idxLocal = idxGlobal - yStartGlobal + 1
            if(ok.AND.all(idxLocal.GE.0).AND.all(idxLocal.LE.yLocalCount+1)) then
                streamInterpValidY(alpha,j) = .true.
                streamInterpIndexY(alpha,j,:) = idxLocal
                streamInterpWeightY(alpha,j,:) = w
            elseif(ok) then
                streamStencilFallbackLocal = streamStencilFallbackLocal + 1
            endif
        enddo
    enddo

    call MPI_ALLREDUCE(streamStencilFallbackLocal, streamStencilFallbackGlobal, 1, &
        MPI_INTEGER, MPI_SUM, COMM2D, IERR)

    return
  end subroutine build_islbm_streaming_stencils_mpi
!===================================================================================================
! build_islbm_streaming_stencils_mpi 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: build_streaming_stencil_1d
! 作用: 对 off-lattice streaming 目标点, 以当前到达节点为中心选择三点模板。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine build_streaming_stencil_1d(n, xnodes, nodeIndex, target, idx, w, ok)
    implicit none
    integer(kind=4), intent(in) :: n, nodeIndex
    real(kind=8), intent(in) :: xnodes(n), target
    integer(kind=4), intent(out) :: idx(3)
    real(kind=8), intent(out) :: w(3)
    logical, intent(out) :: ok
    real(kind=8) :: xloc(3)
    real(kind=8), parameter :: tol = 1.0d-12

    idx = (/1, 1, 1/)
    w = 0.0d0
    ok = .false.
    if(n.LT.3) return
    if((target.LT.xnodes(1)-tol).OR.(target.GT.xnodes(n)+tol)) return

    if(nodeIndex.LE.1) then
        idx = (/1, 2, 3/)
    elseif(nodeIndex.GE.n) then
        idx = (/n-2, n-1, n/)
    else
        idx = (/nodeIndex-1, nodeIndex, nodeIndex+1/)
    endif

    xloc(1) = xnodes(idx(1))
    xloc(2) = xnodes(idx(2))
    xloc(3) = xnodes(idx(3))
    call build_lagrange_weights_3(xloc, target, w)
    ok = .true.

    return
  end subroutine build_streaming_stencil_1d
!===================================================================================================
! build_streaming_stencil_1d 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: build_lagrange_stencil_1d
! 作用: 对一个目标点选择三点二次 Lagrange 插值模板。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine build_lagrange_stencil_1d(n, xnodes, target, idx, w, ok)
    implicit none
    integer(kind=4), intent(in) :: n
    real(kind=8), intent(in) :: xnodes(n), target
    integer(kind=4), intent(out) :: idx(3)
    real(kind=8), intent(out) :: w(3)
    logical, intent(out) :: ok
    integer(kind=4) :: mid
    real(kind=8) :: xloc(3)
    real(kind=8), parameter :: tol = 1.0d-12

    idx = (/1, 1, 1/)
    w = 0.0d0
    ok = .false.
    if(n.LT.3) return
    if((target.LT.xnodes(1)-tol).OR.(target.GT.xnodes(n)+tol)) return

    if(target.LE.xnodes(2)) then
        idx = (/1, 2, 3/)
    elseif(target.GE.xnodes(n-1)) then
        idx = (/n-2, n-1, n/)
    else
        mid = 2
        do while((mid.LT.n-1).AND.(xnodes(mid+1).LT.target))
            mid = mid + 1
        enddo
        idx = (/mid-1, mid, mid+1/)
    endif

    xloc(1) = xnodes(idx(1))
    xloc(2) = xnodes(idx(2))
    xloc(3) = xnodes(idx(3))
    call build_lagrange_weights_3(xloc, target, w)
    ok = .true.

    return
  end subroutine build_lagrange_stencil_1d
!===================================================================================================
! build_lagrange_stencil_1d 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: build_lagrange_weights_3
! 作用: 计算三点二次 Lagrange 插值权重。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine build_lagrange_weights_3(xnode, x0, w)
    implicit none
    real(kind=8), intent(in) :: xnode(3), x0
    real(kind=8), intent(out) :: w(3)

    w(1) = ((x0-xnode(2))*(x0-xnode(3)))/((xnode(1)-xnode(2))*(xnode(1)-xnode(3)))
    w(2) = ((x0-xnode(1))*(x0-xnode(3)))/((xnode(2)-xnode(1))*(xnode(2)-xnode(3)))
    w(3) = ((x0-xnode(1))*(x0-xnode(2)))/((xnode(3)-xnode(1))*(xnode(3)-xnode(2)))

    return
  end subroutine build_lagrange_weights_3
!===================================================================================================
! build_lagrange_weights_3 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: lagrange_derivative_3
! 作用: 执行本子程序对应的初始化、迁移、碰撞、边界、通信或后处理步骤。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine lagrange_derivative_3(xnode, fnode, x0, derivativeValue)
    implicit none
    real(kind=8), intent(in) :: xnode(3), fnode(3), x0
    real(kind=8), intent(out) :: derivativeValue

    derivativeValue = &
        fnode(1)*(2.0d0*x0-xnode(2)-xnode(3))/((xnode(1)-xnode(2))*(xnode(1)-xnode(3))) + &
        fnode(2)*(2.0d0*x0-xnode(1)-xnode(3))/((xnode(2)-xnode(1))*(xnode(2)-xnode(3))) + &
        fnode(3)*(2.0d0*x0-xnode(1)-xnode(2))/((xnode(3)-xnode(1))*(xnode(3)-xnode(2)))

    return
  end subroutine lagrange_derivative_3
!===================================================================================================
! lagrange_derivative_3 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: integrate_lagrange_3_segment
! 作用: 执行本子程序对应的初始化、迁移、碰撞、边界、通信或后处理步骤。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine integrate_lagrange_3_segment(xnode, fnode, xLeft, xRight, integralValue)
    implicit none
    real(kind=8), intent(in) :: xnode(3), fnode(3), xLeft, xRight
    real(kind=8), intent(out) :: integralValue
    integer(kind=4) :: k, m, n
    real(kind=8) :: denom, nodeSum, nodeProduct, basisIntegral

    integralValue = 0.0d0
    do k = 1, 3
      if(k.EQ.1) then
        m = 2
        n = 3
      elseif(k.EQ.2) then
        m = 1
        n = 3
      else
        m = 1
        n = 2
      endif
      denom = (xnode(k)-xnode(m))*(xnode(k)-xnode(n))
      nodeSum = xnode(m) + xnode(n)
      nodeProduct = xnode(m)*xnode(n)
      basisIntegral = ((xRight**3-xLeft**3)/3.0d0 - &
          0.5d0*nodeSum*(xRight**2-xLeft**2) + &
          nodeProduct*(xRight-xLeft))/denom
      integralValue = integralValue + fnode(k)*basisIntegral
    enddo

    return
  end subroutine integrate_lagrange_3_segment
!===================================================================================================
! integrate_lagrange_3_segment 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: exit_data_2d_openacc
! 作用: 计算结束后释放 OpenACC 设备端驻留数据。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine exit_data_2d_openacc()
    use openacc
    use commondata
    implicit none

    if(.not.accDataResident) return
#ifdef steadyFlow
    !$acc exit data delete(up,vp,Tp)
#endif
    !$acc exit data delete(f_post,g_post,u,v,T,rho,f,g,Fx,Fy,Bx_prev,By_prev)
    !$acc exit data delete(streamInterpIndexX,streamInterpIndexY,streamInterpWeightX, &
    !$acc& streamInterpWeightY,streamInterpValidX,streamInterpValidY)
    !$acc exit data delete(xp,yp,quadWidthX,quadWidthY,quadSumX,quadSumY,quadSumArea,ex,ey,omega,omegaT)
    accDataResident = .false.
    return
  end subroutine exit_data_2d_openacc
!===================================================================================================
! exit_data_2d_openacc 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================



!===================================================================================================
! 子程序: exchange_f_post_halo_mpi
! 作用: 交换流场碰撞后分布函数 f_post 的一层 halo。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine exchange_f_post_halo_mpi()
    use commondata
    implicit none
    integer(kind=4) :: i, j, idx, rowCount, colCount
    integer :: status(MPI_STATUS_SIZE)

    !$acc wait(1)
    rowCount = 3*xLocalCount
    ! GPU-aware MPI 路径：在 GPU 上 pack 连续缓冲区，再用 host_data use_device 把设备端地址传给 MPI_SENDRECV。
    ! pull streaming 只会从 ghost 层读取跨界迁移回 owned 区域的方向，所以每条边只需要 3 个 D2Q9 方向。
    ! y 向发送：
    !   下边界 j=1 发给 down 的方向为 4,7,8；上边界 j=yLocalCount 发给 top 的方向为 2,5,6。
    ! 这里统一清零接收缓冲区，避免保留旧设备值。
    !$acc parallel loop gang vector present(fHaloRecvDown,fHaloRecvUp)
    do idx = 1, rowCount
        fHaloRecvDown(idx) = 0.0d0
        fHaloRecvUp(idx) = 0.0d0
    enddo
    !$acc parallel loop gang vector present(f_post,fHaloSendDown,fHaloSendUp) private(idx)
    do i = 1, xLocalCount
        idx = 3*(i-1)
        fHaloSendDown(idx+1) = f_post(i,1,4)
        fHaloSendDown(idx+2) = f_post(i,1,7)
        fHaloSendDown(idx+3) = f_post(i,1,8)
        fHaloSendUp(idx+1) = f_post(i,yLocalCount,2)
        fHaloSendUp(idx+2) = f_post(i,yLocalCount,5)
        fHaloSendUp(idx+3) = f_post(i,yLocalCount,6)
    enddo
    ! MPI_SENDRECV 参数含义：
    !   fHaloSendDown             : 发送缓冲区，当前 rank 的下边界拥有行 j=1。
    !   rowCount                   : 发送/接收元素数，3 个跨 y 边界方向乘以 x 方向 owned 长度。
    !   MPI_DOUBLE_PRECISION       : 数据类型，对应 real(kind=8)。
    !   down, 101                  : 发送目标 rank 和发送 tag；这里发给下邻居。
    !   fHaloRecvUp               : 接收缓冲区，接收到的数据后面解包到当前 rank 的上 halo 行。
    !   top, 101                   : 接收来源 rank 和接收 tag；这里从上邻居接收。
    !   COMM2D, status, IERR       : 笛卡尔通信器、通信状态数组、MPI 错误码。
    !$acc host_data use_device(fHaloSendDown,fHaloRecvUp,fHaloSendUp,fHaloRecvDown)
    call MPI_SENDRECV(fHaloSendDown, rowCount, MPI_DOUBLE_PRECISION, down, 101, &
        fHaloRecvUp, rowCount, MPI_DOUBLE_PRECISION, top, 101, COMM2D, status, IERR)
    call MPI_SENDRECV(fHaloSendUp, rowCount, MPI_DOUBLE_PRECISION, top, 102, &
        fHaloRecvDown, rowCount, MPI_DOUBLE_PRECISION, down, 102, COMM2D, status, IERR)
    !$acc end host_data
    !$acc parallel loop gang vector present(f_post,fHaloRecvDown,fHaloRecvUp) private(idx)
    do i = 1, xLocalCount
        idx = 3*(i-1)
        f_post(i,0,2) = fHaloRecvDown(idx+1)
        f_post(i,0,5) = fHaloRecvDown(idx+2)
        f_post(i,0,6) = fHaloRecvDown(idx+3)
        f_post(i,yLocalCount+1,4) = fHaloRecvUp(idx+1)
        f_post(i,yLocalCount+1,7) = fHaloRecvUp(idx+2)
        f_post(i,yLocalCount+1,8) = fHaloRecvUp(idx+3)
    enddo

    ! x 方向固定 i 的列在 Fortran 内存里不是连续块，所以先在 GPU 上 pack 到一维连续缓冲区。
    ! 这里 j 覆盖 0:yLocalCount+1，把刚刚 y 方向交换得到的上下 halo 也带上，用于补齐角点 halo。
    ! x 向发送：
    !   左边界 i=1 发给 left 的方向为 3,6,7；右边界 i=xLocalCount 发给 right 的方向为 1,5,8。
    colCount = 3*(yLocalCount+2)
    !$acc parallel loop gang vector present(fHaloRecvLeft,fHaloRecvRight)
    do idx = 1, colCount
        fHaloRecvLeft(idx) = 0.0d0
        fHaloRecvRight(idx) = 0.0d0
    enddo
    !$acc parallel loop gang vector present(f_post,fHaloSendLeft,fHaloSendRight) private(idx)
    do j = 0, yLocalCount+1
        idx = 3*j
        fHaloSendLeft(idx+1) = f_post(1,j,3)
        fHaloSendLeft(idx+2) = f_post(1,j,6)
        fHaloSendLeft(idx+3) = f_post(1,j,7)
        fHaloSendRight(idx+1) = f_post(xLocalCount,j,1)
        fHaloSendRight(idx+2) = f_post(xLocalCount,j,5)
        fHaloSendRight(idx+3) = f_post(xLocalCount,j,8)
    enddo

    !$acc host_data use_device(fHaloSendRight,fHaloRecvLeft,fHaloSendLeft,fHaloRecvRight)
    call MPI_SENDRECV(fHaloSendRight, colCount, MPI_DOUBLE_PRECISION, right, 103, &
        fHaloRecvLeft, colCount, MPI_DOUBLE_PRECISION, left, 103, COMM2D, status, IERR)
    call MPI_SENDRECV(fHaloSendLeft, colCount, MPI_DOUBLE_PRECISION, left, 104, &
        fHaloRecvRight, colCount, MPI_DOUBLE_PRECISION, right, 104, COMM2D, status, IERR)
    !$acc end host_data

    !$acc parallel loop gang vector present(f_post,fHaloRecvLeft,fHaloRecvRight) private(idx)
    do j = 0, yLocalCount+1
        idx = 3*j
        f_post(0,j,1) = fHaloRecvLeft(idx+1)
        f_post(0,j,5) = fHaloRecvLeft(idx+2)
        f_post(0,j,8) = fHaloRecvLeft(idx+3)
        f_post(xLocalCount+1,j,3) = fHaloRecvRight(idx+1)
        f_post(xLocalCount+1,j,6) = fHaloRecvRight(idx+2)
        f_post(xLocalCount+1,j,7) = fHaloRecvRight(idx+3)
    enddo

    return
  end subroutine exchange_f_post_halo_mpi
!===================================================================================================
! exchange_f_post_halo_mpi 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: exchange_g_post_halo_mpi
! 作用: 交换温度场碰撞后分布函数 g_post 的一层 halo。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine exchange_g_post_halo_mpi()
    use commondata
    implicit none
    integer(kind=4) :: i, j, rowCount, colCount
    integer :: status(MPI_STATUS_SIZE)

    !$acc wait(1)
    rowCount = xLocalCount
    ! 温度分布函数同样采用 GPU-aware MPI：先在 GPU 上 pack 到连续缓冲区，再传设备地址给 MPI。
    ! D2Q5 温度场每条边只有 1 个方向跨界：
    !   下边界发 4，上边界发 2；收到后分别写入下 halo 的 2 和上 halo 的 4。
    !$acc parallel loop gang vector present(gHaloRecvDown,gHaloRecvUp)
    do i = 1, rowCount
        gHaloRecvDown(i) = 0.0d0
        gHaloRecvUp(i) = 0.0d0
    enddo
    !$acc parallel loop gang vector present(g_post,gHaloSendDown,gHaloSendUp)
    do i = 1, xLocalCount
        gHaloSendDown(i) = g_post(i,1,4)
        gHaloSendUp(i) = g_post(i,yLocalCount,2)
    enddo
    !$acc host_data use_device(gHaloSendDown,gHaloRecvUp,gHaloSendUp,gHaloRecvDown)
    call MPI_SENDRECV(gHaloSendDown, rowCount, MPI_DOUBLE_PRECISION, down, 201, &
        gHaloRecvUp, rowCount, MPI_DOUBLE_PRECISION, top, 201, COMM2D, status, IERR)
    call MPI_SENDRECV(gHaloSendUp, rowCount, MPI_DOUBLE_PRECISION, top, 202, &
        gHaloRecvDown, rowCount, MPI_DOUBLE_PRECISION, down, 202, COMM2D, status, IERR)
    !$acc end host_data
    !$acc parallel loop gang vector present(g_post,gHaloRecvDown,gHaloRecvUp)
    do i = 1, xLocalCount
        g_post(i,0,2) = gHaloRecvDown(i)
        g_post(i,yLocalCount+1,4) = gHaloRecvUp(i)
    enddo

    ! x 向同理：左边界发 3，右边界发 1；D2Q5 没有对角方向，所以不需要交换角点 halo。
    colCount = yLocalCount
    !$acc parallel loop gang vector present(gHaloRecvLeft,gHaloRecvRight)
    do j = 1, colCount
        gHaloRecvLeft(j) = 0.0d0
        gHaloRecvRight(j) = 0.0d0
    enddo
    !$acc parallel loop gang vector present(g_post,gHaloSendLeft,gHaloSendRight)
    do j = 1, yLocalCount
        gHaloSendLeft(j) = g_post(1,j,3)
        gHaloSendRight(j) = g_post(xLocalCount,j,1)
    enddo

    !$acc host_data use_device(gHaloSendRight,gHaloRecvLeft,gHaloSendLeft,gHaloRecvRight)
    call MPI_SENDRECV(gHaloSendRight, colCount, MPI_DOUBLE_PRECISION, right, 203, &
        gHaloRecvLeft, colCount, MPI_DOUBLE_PRECISION, left, 203, COMM2D, status, IERR)
    call MPI_SENDRECV(gHaloSendLeft, colCount, MPI_DOUBLE_PRECISION, left, 204, &
        gHaloRecvRight, colCount, MPI_DOUBLE_PRECISION, right, 204, COMM2D, status, IERR)
    !$acc end host_data

    !$acc parallel loop gang vector present(g_post,gHaloRecvLeft,gHaloRecvRight)
    do j = 1, yLocalCount
        g_post(0,j,1) = gHaloRecvLeft(j)
        g_post(xLocalCount+1,j,3) = gHaloRecvRight(j)
    enddo

    return
  end subroutine exchange_g_post_halo_mpi
!===================================================================================================
! exchange_g_post_halo_mpi 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: gather_output_fields_mpi
! 作用: 汇总宏观场到 root，供输出和最终后处理使用。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine gather_output_fields_mpi()
    use commondata
    implicit none
    call update_host_snapshot_2d_openacc()   ! 先把设备端局部 u/v/T/rho 更新到 host，再执行 MPI 汇总到 u_all/v_all/T_all/rho_all
    call gather_scalar_field_mpi(u, u_all, 301)
    call gather_scalar_field_mpi(v, v_all, 401)
    call gather_scalar_field_mpi(T, T_all, 501)
    call gather_scalar_field_mpi(rho, rho_all, 601)
    return
  end subroutine gather_output_fields_mpi
!===================================================================================================
! gather_output_fields_mpi 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: gather_restart_fields_mpi
! 作用: 汇总严格重启所需的分布函数和历史项到 root。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine gather_restart_fields_mpi()
    use commondata
    implicit none
    call update_host_reload_2d_openacc()     ! restart 写出前需要 host 端拿到设备上的 f/g/Bx_prev/By_prev 最新值
    call gather_distribution_field_mpi(f, f_all, 8, 701)
    call gather_distribution_field_mpi(g, g_all, 4, 801)
#ifdef EnableUseG
    call gather_scalar_field_mpi(Bx_prev, Bx_prev_all, 901)
    call gather_scalar_field_mpi(By_prev, By_prev_all, 1001)
#endif
    return
  end subroutine gather_restart_fields_mpi
!===================================================================================================
! gather_restart_fields_mpi 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: gather_scalar_field_mpi
! 作用: 将局部二维标量场按全局坐标拼到 root 的 globalField。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine gather_scalar_field_mpi(localField, globalField, tagBase)
    use commondata
    implicit none
    real(kind=8), intent(in) :: localField(xLocalCount,yLocalCount)
    real(kind=8), intent(inout) :: globalField(nx,ny)
    integer(kind=4), intent(in) :: tagBase
    integer(kind=4) :: sourceRank, i, j, idx, count
    integer :: status(MPI_STATUS_SIZE)
    real(kind=8), allocatable :: buffer(:)

    count = xLocalCount*yLocalCount
    allocate(buffer(count))
    idx = 0
    do j = 1, yLocalCount                  !每个 rank 把自己的二维局部块 pack 成一维 buffer
        do i = 1, xLocalCount
            idx = idx + 1
            buffer(idx) = localField(i,j)
        enddo
    enddo

    if(isRoot) then
        globalField = 0.0d0
        !按全局起点解包二维标量场
        call unpack_scalar_buffer_mpi(buffer, xLocalCount, yLocalCount, xStartGlobal, yStartGlobal, globalField) 
        ! root 先解包自己的局部块，然后按 rank 顺序阻塞接收其它 rank 的 buffer。
        ! MPI_RECV 没收到匹配的 source/tag 前会停在这里；收到后才继续调用 unpack 写入全局数组。
        do sourceRank = 1, NPROC-1
            count = XLocalCountAll(sourceRank)*YLocalCountAll(sourceRank)
            if(size(buffer).NE.count) then
                deallocate(buffer)
                allocate(buffer(count))
            endif
            ! MPI_RECV 参数含义：
            !   buffer       : 接收缓冲区，用来存放 sourceRank 发来的局部二维标量场。
            !   count        : 接收元素数量，等于该 sourceRank 的 x/y 局部网格数乘积。
            !   MPI_DOUBLE_PRECISION : 数据类型，对应 real(kind=8)。
            !   sourceRank   : 指定从哪个 MPI rank 接收。
            !   tagBase+sourceRank : 消息 tag；tagBase 区分变量，sourceRank 区分发送进程。
            !   COMM2D/status/IERR : 笛卡尔通信器、通信状态数组、MPI 错误码。
            call MPI_RECV(buffer, count, MPI_DOUBLE_PRECISION, sourceRank, tagBase+sourceRank, COMM2D, status, IERR)
            call unpack_scalar_buffer_mpi(buffer, XLocalCountAll(sourceRank), YLocalCountAll(sourceRank), &
                XStartGlobalAll(sourceRank), YStartGlobalAll(sourceRank), globalField)
        enddo
    else
        call MPI_SEND(buffer, count, MPI_DOUBLE_PRECISION, 0, tagBase+MYID, COMM2D, IERR)
    endif

    deallocate(buffer)
    return
  end subroutine gather_scalar_field_mpi
!===================================================================================================
! gather_scalar_field_mpi 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: unpack_scalar_buffer_mpi
! 作用: root 端按全局起点解包二维标量场。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine unpack_scalar_buffer_mpi(buffer, blockNx, blockNy, blockXStart, blockYStart, globalField)
    use commondata, only: nx, ny
    implicit none
    real(kind=8), intent(in) :: buffer(*)
    integer(kind=4), intent(in) :: blockNx, blockNy, blockXStart, blockYStart
    real(kind=8), intent(inout) :: globalField(nx,ny)
    integer(kind=4) :: i, j, idx

    idx = 0
    do j = 1, blockNy
        do i = 1, blockNx
            idx = idx + 1
            globalField(blockXStart+i-1, blockYStart+j-1) = buffer(idx)
        enddo
    enddo
    return
  end subroutine unpack_scalar_buffer_mpi
!===================================================================================================
! unpack_scalar_buffer_mpi 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: gather_distribution_field_mpi
! 作用: 将局部分布函数按全局坐标拼到 root 的 globalField。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine gather_distribution_field_mpi(localField, globalField, qMax, tagBase)
    use commondata
    implicit none
    integer(kind=4), intent(in) :: qMax, tagBase
    real(kind=8), intent(in) :: localField(xLocalCount,yLocalCount,0:qMax)
    real(kind=8), intent(inout) :: globalField(nx,ny,0:qMax)
    integer(kind=4) :: sourceRank, i, j, alpha, idx, count
    integer :: status(MPI_STATUS_SIZE)
    real(kind=8), allocatable :: buffer(:)

    count = (qMax+1)*xLocalCount*yLocalCount
    allocate(buffer(count))
    idx = 0
    do alpha = 0, qMax
        do j = 1, yLocalCount
            do i = 1, xLocalCount
                idx = idx + 1
                buffer(idx) = localField(i,j,alpha)
            enddo
        enddo
    enddo

    if(isRoot) then
        globalField = 0.0d0
        call unpack_distribution_buffer_mpi(buffer, xLocalCount, yLocalCount, xStartGlobal, yStartGlobal, qMax, globalField)
        do sourceRank = 1, NPROC-1
            count = (qMax+1)*XLocalCountAll(sourceRank)*YLocalCountAll(sourceRank)
            if(size(buffer).NE.count) then
                deallocate(buffer)
                allocate(buffer(count))
            endif
            call MPI_RECV(buffer, count, MPI_DOUBLE_PRECISION, sourceRank, tagBase+sourceRank, COMM2D, status, IERR)
            call unpack_distribution_buffer_mpi(buffer, XLocalCountAll(sourceRank), YLocalCountAll(sourceRank), &
                XStartGlobalAll(sourceRank), YStartGlobalAll(sourceRank), qMax, globalField)
        enddo
    else
        call MPI_SEND(buffer, count, MPI_DOUBLE_PRECISION, 0, tagBase+MYID, COMM2D, IERR)
    endif

    deallocate(buffer)
    return
  end subroutine gather_distribution_field_mpi
!===================================================================================================
! gather_distribution_field_mpi 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: unpack_distribution_buffer_mpi
! 作用: root 端按全局起点解包分布函数。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine unpack_distribution_buffer_mpi(buffer, blockNx, blockNy, blockXStart, blockYStart, qMax, globalField)
    use commondata, only: nx, ny
    implicit none
    real(kind=8), intent(in) :: buffer(*)
    integer(kind=4), intent(in) :: blockNx, blockNy, blockXStart, blockYStart, qMax
    real(kind=8), intent(inout) :: globalField(nx,ny,0:qMax)
    integer(kind=4) :: i, j, alpha, idx

    idx = 0
    do alpha = 0, qMax
        do j = 1, blockNy
            do i = 1, blockNx
                idx = idx + 1
                globalField(blockXStart+i-1, blockYStart+j-1, alpha) = buffer(idx)
            enddo
        enddo
    enddo
    return
  end subroutine unpack_distribution_buffer_mpi
!===================================================================================================
! unpack_distribution_buffer_mpi 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: read_reload_fields_mpi
! 作用: root 读全局重启文件，然后按当前 MPI 分解把局部块发送给各 rank。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine read_reload_fields_mpi(reloadFileName)
    use commondata
    implicit none
    character(len=*), intent(in) :: reloadFileName
    integer(kind=4) :: i, j, alpha
    integer(kind=4) :: sourceRank, blockNx, blockNy, blockXStart, blockYStart
    integer(kind=4) :: countF, countG, countScalar
    integer :: status(MPI_STATUS_SIZE)
    real(kind=8), allocatable :: bufferF(:), bufferG(:), bufferScalar(:)

    if(isRoot) then
        f_all = 0.0d0
        g_all = 0.0d0
        Bx_prev_all = 0.0d0
        By_prev_all = 0.0d0
        open(unit=01,file=trim(reloadFilePrefix)//"-"//trim(reloadFileName)//".bin",form="unformatted", &
        access="sequential",status='old')
        write(00,*) "Reloading f, g and optional UseG history from file"
        ! f_all/g_all 的数组形状是 (i,j,alpha)。Fortran 列主序下第一维 i 连续，
        ! 所以 restart 文件读写让 i 放在最内层，按内存连续方向依次读取。
        read(01) (((f_all(i,j,alpha), i=1,nx), j=1,ny), alpha=0,8)
        read(01) (((g_all(i,j,alpha), i=1,nx), j=1,ny), alpha=0,4)
#ifdef EnableUseG
        read(01) ((Bx_prev_all(i,j), i=1,nx), j=1,ny)
        read(01) ((By_prev_all(i,j), i=1,nx), j=1,ny)
#endif
        close(01)
    endif

    if(isRoot) then
        do sourceRank = 0, NPROC-1
            blockNx = XLocalCountAll(sourceRank)
            blockNy = YLocalCountAll(sourceRank)
            blockXStart = XStartGlobalAll(sourceRank)
            blockYStart = YStartGlobalAll(sourceRank)

            countF = 9*blockNx*blockNy
            countG = 5*blockNx*blockNy
            countScalar = blockNx*blockNy
            allocate(bufferF(countF), bufferG(countG), bufferScalar(countScalar))

            call pack_distribution_global_block_mpi(f_all, 8, blockNx, blockNy, blockXStart, blockYStart, bufferF)
            call pack_distribution_global_block_mpi(g_all, 4, blockNx, blockNy, blockXStart, blockYStart, bufferG)

            if(sourceRank.EQ.0) then
                call unpack_distribution_local_buffer_mpi(bufferF, 8, f)
                call unpack_distribution_local_buffer_mpi(bufferG, 4, g)
            else
                ! MPI_SEND 参数含义：
                !   bufferF/countF/MPI_DOUBLE_PRECISION : 要发送给 sourceRank 的局部 f 数据及其数量和类型。
                !   sourceRank, 1100+sourceRank         : 目标 rank 和 tag；tag 区分 f/g/Bx/By 与目标 rank。
                !   COMM2D/IERR                         : 笛卡尔通信器和 MPI 错误码。
                call MPI_SEND(bufferF, countF, MPI_DOUBLE_PRECISION, sourceRank, 1100+sourceRank, COMM2D, IERR)
                call MPI_SEND(bufferG, countG, MPI_DOUBLE_PRECISION, sourceRank, 1200+sourceRank, COMM2D, IERR)
            endif

#ifdef EnableUseG
            call pack_scalar_global_block_mpi(Bx_prev_all, blockNx, blockNy, blockXStart, blockYStart, bufferScalar)
            if(sourceRank.EQ.0) then
                call unpack_scalar_local_buffer_mpi(bufferScalar, Bx_prev)
            else
                call MPI_SEND(bufferScalar, countScalar, MPI_DOUBLE_PRECISION, sourceRank, 1300+sourceRank, COMM2D, IERR)
            endif
            call pack_scalar_global_block_mpi(By_prev_all, blockNx, blockNy, blockXStart, blockYStart, bufferScalar)
            if(sourceRank.EQ.0) then
                call unpack_scalar_local_buffer_mpi(bufferScalar, By_prev)
            else
                call MPI_SEND(bufferScalar, countScalar, MPI_DOUBLE_PRECISION, sourceRank, 1400+sourceRank, COMM2D, IERR)
            endif
#endif
            deallocate(bufferF, bufferG, bufferScalar)
        enddo
    else
        countF = 9*xLocalCount*yLocalCount
        countG = 5*xLocalCount*yLocalCount
        countScalar = xLocalCount*yLocalCount
        allocate(bufferF(countF), bufferG(countG), bufferScalar(countScalar))

        ! MPI_RECV 是阻塞接收：当前 rank 会等到 root 发来匹配 tag 的局部块后才继续解包。
        call MPI_RECV(bufferF, countF, MPI_DOUBLE_PRECISION, 0, 1100+MYID, COMM2D, status, IERR)
        call MPI_RECV(bufferG, countG, MPI_DOUBLE_PRECISION, 0, 1200+MYID, COMM2D, status, IERR)
        call unpack_distribution_local_buffer_mpi(bufferF, 8, f)
        call unpack_distribution_local_buffer_mpi(bufferG, 4, g)

#ifdef EnableUseG
        call MPI_RECV(bufferScalar, countScalar, MPI_DOUBLE_PRECISION, 0, 1300+MYID, COMM2D, status, IERR)
        call unpack_scalar_local_buffer_mpi(bufferScalar, Bx_prev)
        call MPI_RECV(bufferScalar, countScalar, MPI_DOUBLE_PRECISION, 0, 1400+MYID, COMM2D, status, IERR)
        call unpack_scalar_local_buffer_mpi(bufferScalar, By_prev)
#else
        Bx_prev = 0.0d0
        By_prev = 0.0d0
#endif
        deallocate(bufferF, bufferG, bufferScalar)
    endif

#ifndef EnableUseG
    Bx_prev = 0.0d0
    By_prev = 0.0d0
#endif
    return
  end subroutine read_reload_fields_mpi
!===================================================================================================
! read_reload_fields_mpi 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: pack_scalar_global_block_mpi
! 作用: root 按某个 rank 的全局起点和局部尺寸，把全局二维标量场打包成一维 buffer。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine pack_scalar_global_block_mpi(globalField, blockNx, blockNy, blockXStart, blockYStart, buffer)
    use commondata, only: nx, ny
    implicit none
    real(kind=8), intent(in) :: globalField(nx,ny)
    integer(kind=4), intent(in) :: blockNx, blockNy, blockXStart, blockYStart
    real(kind=8), intent(out) :: buffer(*)
    integer(kind=4) :: i, j, idx

    idx = 0
    do j = 1, blockNy
        do i = 1, blockNx
            idx = idx + 1
            buffer(idx) = globalField(blockXStart+i-1, blockYStart+j-1)
        enddo
    enddo
    return
  end subroutine pack_scalar_global_block_mpi
!===================================================================================================
! pack_scalar_global_block_mpi 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: unpack_scalar_local_buffer_mpi
! 作用: 当前 rank 把 root 发来的标量 buffer 解包成本地二维数组。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine unpack_scalar_local_buffer_mpi(buffer, localField)
    use commondata
    implicit none
    real(kind=8), intent(in) :: buffer(*)
    real(kind=8), intent(out) :: localField(xLocalCount,yLocalCount)
    integer(kind=4) :: i, j, idx

    idx = 0
    do j = 1, yLocalCount
        do i = 1, xLocalCount
            idx = idx + 1
            localField(i,j) = buffer(idx)
        enddo
    enddo
    return
  end subroutine unpack_scalar_local_buffer_mpi
!===================================================================================================
! unpack_scalar_local_buffer_mpi 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: pack_distribution_global_block_mpi
! 作用: root 按某个 rank 的全局起点和局部尺寸，把全局分布函数打包成一维 buffer。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine pack_distribution_global_block_mpi(globalField, qMax, blockNx, blockNy, blockXStart, blockYStart, buffer)
    use commondata, only: nx, ny
    implicit none
    integer(kind=4), intent(in) :: qMax, blockNx, blockNy, blockXStart, blockYStart
    real(kind=8), intent(in) :: globalField(nx,ny,0:qMax)
    real(kind=8), intent(out) :: buffer(*)
    integer(kind=4) :: i, j, alpha, idx

    idx = 0
    do alpha = 0, qMax
        do j = 1, blockNy
            do i = 1, blockNx
                idx = idx + 1
                buffer(idx) = globalField(blockXStart+i-1, blockYStart+j-1, alpha)
            enddo
        enddo
    enddo
    return
  end subroutine pack_distribution_global_block_mpi
!===================================================================================================
! pack_distribution_global_block_mpi 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: unpack_distribution_local_buffer_mpi
! 作用: 当前 rank 把 root 发来的分布函数 buffer 解包成本地 f/g 数组。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine unpack_distribution_local_buffer_mpi(buffer, qMax, localField)
    use commondata
    implicit none
    integer(kind=4), intent(in) :: qMax
    real(kind=8), intent(in) :: buffer(*)
    real(kind=8), intent(out) :: localField(xLocalCount,yLocalCount,0:qMax)
    integer(kind=4) :: i, j, alpha, idx

    idx = 0
    do alpha = 0, qMax
        do j = 1, yLocalCount
            do i = 1, xLocalCount
                idx = idx + 1
                localField(i,j,alpha) = buffer(idx)
            enddo
        enddo
    enddo
    return
  end subroutine unpack_distribution_local_buffer_mpi
!===================================================================================================
! unpack_distribution_local_buffer_mpi 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: collision
! 作用: 完成流场分布函数 f 的碰撞更新，并处理体力项离散修正。
! 用途: 在主程序时间推进循环中调用，位于 streaming 之前。
!===================================================================================================
  subroutine collision()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: m(0:8), m_post(0:8), meq(0:8)
    real(kind=8) :: s(0:8)
    real(kind=8) :: fSource(0:8)

    !$acc parallel loop gang vector collapse(2) present(f,f_post,rho,u,v,Fx,Fy,T) async(1) &
    !$acc& private(alpha,s,m,m_post,meq,fSource)
    do j = 1, yLocalCount
        do i = 1, xLocalCount

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

          s(0) = 0.0d0      !!s_{\rho}
          s(1) = Snu !!s_{e}
          s(2) = Snu !!s_{\epsilon}
          s(3) = 0.0d0      !!s_{j} 
          s(4) = Sq !!s_{q}
          s(5) = 0.0d0      !!s_{j}
          s(6) = Sq       !!s_{q}
          s(7) = Snu !!s_{\nu}
          s(8) = Snu       !!s_{\nu}

          Fx(i,j) = 0.0d0
          Fy(i,j) = rho(i,j)*gBeta*(T(i,j)-Tref)        !动量方程上的源项，即浮力项


#ifdef    SideHeatedHa
          Fx(i,j) = 0.0d0+B2sigemarho*(v(i,j)*sin(phi)*cos(phi)-u(i,j)*sin(phi)*sin(phi))
          Fy(i,j) = rho(i,j)*gBeta*(T(i,j)-Tref)+ rho(i,j)*B2sigemarho*(u(i,j)*sin(phi)*cos(phi)&
          -v(i,j)*cos(phi)*cos(phi))                    !动量方程上的源项，即浮力项加磁场
#endif


          fSource(0) = 0.0d0                                                       !将源项F对应的贡献投影到各个矩中，并做半步修正
          fSource(1) = (6.0d0-3.0d0*s(1))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
          fSource(2) = -(6.0d0-3.0d0*s(2))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
          fSource(3) = (1.0d0-0.5d0*s(3))*Fx(i,j)
          fSource(4) = -(1.0d0-0.5d0*s(4))*Fx(i,j)
          fSource(5) = (1.0d0-0.5d0*s(5))*Fy(i,j)
          fSource(6) = -(1.0d0-0.5d0*s(6))*Fy(i,j)
          fSource(7) = (2.0d0-s(7))*(u(i,j)*Fx(i,j)-v(i,j)*Fy(i,j))
          fSource(8) = (1.0d0-0.5d0*s(8))*(u(i,j)*Fy(i,j)+v(i,j)*Fx(i,j))     !这边是乘以M变到矩空间，然后再乘以1-1/2S修正

          do alpha = 0, 8
            m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+fSource(alpha)     !矩空间碰撞
          enddo

          f_post(i,j,0) = m_post(0)/9.0d0-m_post(1)/9.0d0+m_post(2)/9.0d0                                         !这边是乘以M逆
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
!===================================================================================================
! collision 结束: 流场碰撞步骤完成。
!===================================================================================================


!===================================================================================================
! 子程序: streaming
! 作用: 完成流场分布函数 f 的迁移，把碰撞后的信息传播到相邻格点。
! 用途: 在主程序时间推进循环中调用，位于 collision 之后、bounceback 之前。
!===================================================================================================
  subroutine streaming()                                    !先迁移，再边界处理
    use commondata                                            !ISLBM迁移：内部rank界面允许使用一层halo做二次插值
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: ip, jp
    integer(kind=4) :: alpha, ii, jj
    real(kind=8) :: value
    logical :: useInterp
    
    !$acc parallel loop gang vector collapse(2) &
    !$acc& present(f,f_post,ex,ey,streamInterpValidX,streamInterpValidY) &
    !$acc& present(streamInterpIndexX,streamInterpIndexY,streamInterpWeightX,streamInterpWeightY) &
    !$acc& async(1) private(ip,jp,alpha,ii,jj,value,useInterp)
    do j = 1, yLocalCount
        do i = 1, xLocalCount
            do alpha = 0, 8
                useInterp = .false.
                if(alpha.EQ.0) then
                    f(i,j,alpha) = f_post(i,j,alpha)
                elseif(ey(alpha).EQ.0) then
                    useInterp = streamInterpValidX(alpha,i)
                    if(useInterp) then
                        value = 0.0d0
                        do ii = 1, 3
                            value = value + streamInterpWeightX(alpha,i,ii)*f_post(streamInterpIndexX(alpha,i,ii),j,alpha)
                        enddo
                        f(i,j,alpha) = value
                    endif
                elseif(ex(alpha).EQ.0) then
                    useInterp = streamInterpValidY(alpha,j)
                    if(useInterp) then
                        value = 0.0d0
                        do jj = 1, 3
                            value = value + streamInterpWeightY(alpha,j,jj)*f_post(i,streamInterpIndexY(alpha,j,jj),alpha)
                        enddo
                        f(i,j,alpha) = value
                    endif
                else
                    useInterp = streamInterpValidX(alpha,i).AND.streamInterpValidY(alpha,j)
                    if(useInterp) then
                        value = 0.0d0
                        do jj = 1, 3
                            do ii = 1, 3
                                value = value + streamInterpWeightX(alpha,i,ii)*streamInterpWeightY(alpha,j,jj)* &
                                    f_post(streamInterpIndexX(alpha,i,ii),streamInterpIndexY(alpha,j,jj),alpha)
                            enddo
                        enddo
                        f(i,j,alpha) = value
                    endif
                endif
                if((.not.useInterp).AND.(alpha.NE.0)) then
                    ip = i-ex(alpha)
                    jp = j-ey(alpha)
                    f(i,j,alpha) = f_post(ip,jp,alpha)
                endif
            enddo
        enddo
    enddo

    return
  end subroutine streaming
!===================================================================================================
! streaming 结束: 完成流场分布函数 f 的迁移，把碰撞后的信息传播到相邻格点。
!===================================================================================================




!===================================================================================================
! 子程序: bounceback
! 作用: 处理流场边界条件，包括无滑移壁面和相关反弹格式。
! 用途: 在主程序时间推进循环中调用，位于 streaming 之后、macro 之前。
!===================================================================================================
  subroutine bounceback()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    ! integer(kind=4) :: alpha

#ifdef VerticalWallsPeriodicalU      
    ! MPI 周期边界由笛卡尔周期 neighbor 的 halo 交换提供，streaming 后无需额外反弹。
#endif

#ifdef VerticalWallsNoslip
    if(hasLeftBoundary) then
        !$acc parallel loop gang vector present(f,f_post) async(1)
        do j = 1, yLocalCount                                         !只在真正拥有左物理侧壁的 rank 上做无滑移反弹
            f(1,j,1) = f_post(1,j,3)
            f(1,j,5) = f_post(1,j,7)
            f(1,j,8) = f_post(1,j,6)
        enddo
    endif
    if(hasRightBoundary) then
        !$acc parallel loop gang vector present(f,f_post) async(1)
        do j = 1, yLocalCount                                         !只在真正拥有右物理侧壁的 rank 上做无滑移反弹
            f(xLocalCount,j,3) = f_post(xLocalCount,j,1)
            f(xLocalCount,j,6) = f_post(xLocalCount,j,8)
            f(xLocalCount,j,7) = f_post(xLocalCount,j,5)
        enddo
    endif
#endif

#ifdef HorizontalWallsNoslip
    if(hasBottomBoundary) then
        !$acc parallel loop gang vector present(f,f_post) async(1)
        do i = 1, xLocalCount                                         !只在真正拥有下物理壁面的 rank 上做无滑移反弹
            f(i,1,2) = f_post(i,1,4)
            f(i,1,5) = f_post(i,1,7)
            f(i,1,6) = f_post(i,1,8)
        enddo
    endif
    if(hasTopBoundary) then
        !$acc parallel loop gang vector present(f,f_post) async(1)
        do i = 1, xLocalCount                                         !只在真正拥有上物理壁面的 rank 上做无滑移反弹
            f(i,yLocalCount,4) = f_post(i,yLocalCount,2)
            f(i,yLocalCount,7) = f_post(i,yLocalCount,5)
            f(i,yLocalCount,8) = f_post(i,yLocalCount,6)
        enddo
    endif
#endif

    return
  end subroutine bounceback
!===================================================================================================
! bounceback 结束: 处理流场边界条件，包括无滑移壁面和相关反弹格式。
!===================================================================================================




!===================================================================================================
! 子程序: macro
! 作用: 由流场分布函数恢复 rho、u、v，并加入半步力项速度修正。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine macro()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !$acc parallel loop gang vector collapse(2) present(f,rho,u,v,Fx,Fy) async(1)
    do j = 1, yLocalCount
        do i = 1, xLocalCount
            rho(i,j) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
            ! 含力LBM的半步动量修正: rho*u = sum(f e) + 0.5*F, 对应Guo forcing的二阶定义。
            u(i,j) = (f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8) + &
              0.5d0*Fx(i,j))/rho(i,j)
            v(i,j) = ( f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)+0.5d0*Fy(i,j) )/rho(i,j)
        enddo
    enddo
    
    return
  end subroutine macro
!===================================================================================================
! macro 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
!===================================================================================================
!宏观量计算结束
!===================================================================================================



!===================================================================================================
! 子程序: collisionT
! 作用: 完成温度分布函数 g 的碰撞更新，并加入热流修正项。
! 用途: 在主程序时间推进循环中调用，位于流场 macro 之后。
!===================================================================================================
    subroutine collisionT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: n(0:4), n_post(0:4), neq(0:4)
    real(kind=8) :: q(0:4)
    real(kind=8) :: Bx, By
    real(kind=8) :: dBx, dBy
    real(kind=8), parameter :: SG = 1.0d0 - 0.5d0*Qk




    !$acc parallel loop gang vector collapse(2) present(g,g_post,u,v,T,Bx_prev,By_prev) async(1) &
    !$acc& private(alpha,n,neq,q,n_post,Bx,By,dBx,dBy)
    do j = 1, yLocalCount
        do i = 1, xLocalCount

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
!===================================================================================================
! collisionT 结束: 完成温度分布函数 g 的碰撞更新，并加入热流修正项。
!===================================================================================================




!===================================================================================================
! 子程序: streamingT
! 作用: 完成温度分布函数 g 的迁移，把碰撞后的温度信息传播到相邻格点。
! 用途: 在主程序时间推进循环中调用，位于 collisionT 之后、bouncebackT 之前。
!===================================================================================================
    subroutine streamingT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: ip, jp
    integer(kind=4) :: alpha, ii, jj
    real(kind=8) :: value
    logical :: useInterp
    
    !$acc parallel loop gang vector collapse(2) &
    !$acc& present(g,g_post,ex,ey,streamInterpValidX,streamInterpValidY) &
    !$acc& present(streamInterpIndexX,streamInterpIndexY,streamInterpWeightX,streamInterpWeightY) &
    !$acc& async(1) private(ip,jp,alpha,ii,jj,value,useInterp)
    do j = 1, yLocalCount
        do i = 1, xLocalCount
            do alpha = 0, 4
                useInterp = .false.
                if(alpha.EQ.0) then
                    g(i,j,alpha) = g_post(i,j,alpha)
                elseif(ey(alpha).EQ.0) then
                    useInterp = streamInterpValidX(alpha,i)
                    if(useInterp) then
                        value = 0.0d0
                        do ii = 1, 3
                            value = value + streamInterpWeightX(alpha,i,ii)*g_post(streamInterpIndexX(alpha,i,ii),j,alpha)
                        enddo
                        g(i,j,alpha) = value
                    endif
                elseif(ex(alpha).EQ.0) then
                    useInterp = streamInterpValidY(alpha,j)
                    if(useInterp) then
                        value = 0.0d0
                        do jj = 1, 3
                            value = value + streamInterpWeightY(alpha,j,jj)*g_post(i,streamInterpIndexY(alpha,j,jj),alpha)
                        enddo
                        g(i,j,alpha) = value
                    endif
                else
                    useInterp = streamInterpValidX(alpha,i).AND.streamInterpValidY(alpha,j)
                    if(useInterp) then
                        value = 0.0d0
                        do jj = 1, 3
                            do ii = 1, 3
                                value = value + streamInterpWeightX(alpha,i,ii)*streamInterpWeightY(alpha,j,jj)* &
                                    g_post(streamInterpIndexX(alpha,i,ii),streamInterpIndexY(alpha,j,jj),alpha)
                            enddo
                        enddo
                        g(i,j,alpha) = value
                    endif
                endif
                if((.not.useInterp).AND.(alpha.NE.0)) then
                    ip = i-ex(alpha)
                    jp = j-ey(alpha)
                    g(i,j,alpha) = g_post(ip,jp,alpha)
                endif
            enddo
        enddo
    enddo

    return
    end subroutine streamingT
!===================================================================================================
! streamingT 结束: 完成温度分布函数 g 的迁移，把碰撞后的温度信息传播到相邻格点。
!===================================================================================================



!===================================================================================================
! 子程序: bouncebackT
! 作用: 处理温度边界条件，包括恒温、绝热和周期边界。
! 用途: 在主程序时间推进循环中调用，位于 streamingT 之后、macroT 之前。
!===================================================================================================
    subroutine bouncebackT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    !integer(kind=4) :: alpha

#ifdef VerticalWallsPeriodicalT 
    ! MPI 周期边界由笛卡尔周期 neighbor 的 halo 交换提供，streamingT 后无需额外覆盖。
#endif

#ifdef VerticalWallsConstT
    if(hasLeftBoundary) then
        !$acc parallel loop gang vector present(g,g_post,omegaT) async(1)
        do j = 1, yLocalCount
#ifdef EnableLegacyThermalScheme
            g(1,j,1) = -g_post(1,j,3)+(4.0d0+paraA)/10.0d0*Thot
#else
            g(1,j,1) = -g_post(1,j,3)+2.0d0*omegaT(3)*Thot
#endif
        enddo
    endif
    if(hasRightBoundary) then
        !$acc parallel loop gang vector present(g,g_post,omegaT) async(1)
        do j = 1, yLocalCount
#ifdef EnableLegacyThermalScheme
            g(xLocalCount,j,3) = -g_post(xLocalCount,j,1)+(4.0d0+paraA)/10.0d0*Tcold
#else
            g(xLocalCount,j,3) = -g_post(xLocalCount,j,1)+2.0d0*omegaT(1)*Tcold
#endif
        enddo
    endif
#endif

#ifdef VerticalWallsAdiabatic
    if(hasLeftBoundary) then
        !$acc parallel loop gang vector present(g,g_post) async(1)
        do j = 1, yLocalCount
            g(1,j,1) = g_post(1,j,3)
        enddo
    endif
    if(hasRightBoundary) then
        !$acc parallel loop gang vector present(g,g_post) async(1)
        do j = 1, yLocalCount
            g(xLocalCount,j,3) = g_post(xLocalCount,j,1)
        enddo
    endif
#endif

#ifdef HorizontalWallsAdiabatic
    if(hasBottomBoundary) then
        !$acc parallel loop gang vector present(g,g_post) async(1)
        do i = 1, xLocalCount
            g(i,1,2) = g_post(i,1,4)
        enddo
    endif
    if(hasTopBoundary) then
        !$acc parallel loop gang vector present(g,g_post) async(1)
        do i = 1, xLocalCount
            g(i,yLocalCount,4) = g_post(i,yLocalCount,2)
        enddo
    endif
#endif

#ifdef HorizontalWallsConstT
    if(hasBottomBoundary) then
        !$acc parallel loop gang vector present(g,g_post,omegaT) async(1)
        do i = 1, xLocalCount
#ifdef EnableLegacyThermalScheme
            g(i,1,2) = -g_post(i,1,4)+(4.0d0+paraA)/10.0d0*Thot
#else
            g(i,1,2) = -g_post(i,1,4)+2.0d0*omegaT(4)*Thot
#endif
        enddo
    endif
    if(hasTopBoundary) then
        !$acc parallel loop gang vector present(g,g_post,omegaT) async(1)
        do i = 1, xLocalCount
#ifdef EnableLegacyThermalScheme
            g(i,yLocalCount,4) = -g_post(i,yLocalCount,2)+(4.0d0+paraA)/10.0d0*Tcold
#else
            g(i,yLocalCount,4) = -g_post(i,yLocalCount,2)+2.0d0*omegaT(2)*Tcold
#endif
        enddo
    endif
#endif

    return
    end subroutine bouncebackT
!===================================================================================================
! bouncebackT 结束: 处理温度边界条件，包括恒温、绝热和周期边界。
!===================================================================================================




!===================================================================================================
! 子程序: macroT
! 作用: 由温度分布函数恢复宏观温度场 T。
! 用途: 在主程序时间推进循环中调用，作为温度更新链条的最后一步。
!===================================================================================================
    subroutine macroT()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !$acc parallel loop gang vector collapse(2) present(g,T) async(1)
    do j = 1, yLocalCount
        do i = 1, xLocalCount
            T(i,j) = g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
        enddo
    enddo
    
    return
    end subroutine macroT
!===================================================================================================
! macroT 结束: 由温度分布函数恢复宏观温度场 T。
!===================================================================================================



!===================================================================================================
! 子程序: reconstruct_macro_from_fg
! 作用: 从重启读回的 f/g 重新恢复 rho、u、v 和 T。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
    subroutine reconstruct_macro_from_fg()
    use commondata
    implicit none
    integer(kind=4) :: i, j, iter
    real(kind=8) :: momx, momy
    logical :: rho_bad

    ! 重启重建发生在 enter_data_2d_openacc() 之前，这里走 host 端串行恢复，避免调用设备端 macroT。
    do j = 1, yLocalCount
        do i = 1, xLocalCount
            T(i,j) = g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
        enddo
    enddo

    rho_bad = .false.
    do j = 1, yLocalCount
        do i = 1, xLocalCount
            rho(i,j) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
            momx = f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)
            momy = f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)

            if (rho(i,j).GT.0.0d0) then
                u(i,j) = momx/rho(i,j)
                v(i,j) = momy/rho(i,j)

                do iter = 1, 3
                    Fx(i,j) = 0.0d0
                    Fy(i,j) = rho(i,j)*gBeta*(T(i,j)-Tref)

#ifdef    SideHeatedHa
                    Fx(i,j) = B2sigemarho*(v(i,j)*sin(phi)*cos(phi)-u(i,j)*sin(phi)*sin(phi))
                    Fy(i,j) = rho(i,j)*gBeta*(T(i,j)-Tref)+rho(i,j)*B2sigemarho*(u(i,j)*sin(phi)*cos(phi)&
                    -v(i,j)*cos(phi)*cos(phi))
#endif

                    u(i,j) = (momx + 0.5d0*Fx(i,j))/rho(i,j)
                    v(i,j) = (momy + 0.5d0*Fy(i,j))/rho(i,j)
                enddo
            else
                rho_bad = .true.
                Fx(i,j) = 0.0d0
                Fy(i,j) = 0.0d0
                u(i,j) = 0.0d0
                v(i,j) = 0.0d0
            endif
        enddo
    enddo

    if (rho_bad) then
        write(*,*) "Warning: non-positive rho found during restart reconstruction."
        stop
    endif

    return
    end subroutine reconstruct_macro_from_fg
!===================================================================================================
! reconstruct_macro_from_fg 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
!===================================================================================================
! reconstruct_macro_from_fg end: current macro state is rebuilt; EnableUseG history is read separately
!===================================================================================================



#ifdef steadyFlow
!===================================================================================================
! 子程序: check
! 作用: 计算稳态收敛误差并写出收敛历史。
! 用途: 在 steadyFlow 模式下由主程序定期调用。
!===================================================================================================
    subroutine check()
    use commondata
    implicit none
    integer(kind=4) :: i, j, globalI, globalJ
    real(kind=8) :: error1, error2, error5, error6
    real(kind=8) :: error1Global, error2Global, error5Global, error6Global
    real(kind=8) :: areaWeight
    character(len=64) :: caseTag



    error1 = 0.0d0
    error2 = 0.0d0

    error5 = 0.0d0
    error6 = 0.0d0
    
    !$acc wait(1)
    !$acc parallel loop collapse(2) present(u,up,v,vp,T,Tp,quadWidthX,quadWidthY) &
    !$acc& private(globalI,globalJ,areaWeight) reduction(+:error1,error2,error5,error6)
    do j = 1, yLocalCount
        do i = 1, xLocalCount
            globalI = xStartGlobal + i - 1
            globalJ = yStartGlobal + j - 1
            areaWeight = quadWidthX(globalI)*quadWidthY(globalJ)
            error1 = error1+areaWeight*((u(i,j)-up(i,j))*(u(i,j)-up(i,j)) + &
                (v(i,j)-vp(i,j))*(v(i,j)-vp(i,j)))
            error2 = error2+areaWeight*(u(i,j)*u(i,j)+v(i,j)*v(i,j))
                
            error5 = error5+areaWeight*dABS( T(i,j)-Tp(i,j) )
            error6 = error6+areaWeight*dABS( T(i,j) )
                
            up(i,j) = u(i,j)
            vp(i,j) = v(i,j)
            Tp(i,j) = T(i,j)
        enddo
    enddo
    
    ! MPI_ALLREDUCE 参数含义：
    !   error1               : 当前 rank 的局部误差累加值。
    !   error1Global         : 输出的全局误差累加值，所有 rank 都会得到同一份结果。
    !   1                    : 每个 rank 只归约 1 个标量。
    !   MPI_DOUBLE_PRECISION : 数据类型，对应 real(kind=8)。
    !   MPI_SUM              : 归约操作，对所有 rank 的局部值求和。
    !   COMM2D/IERR          : 笛卡尔通信器和 MPI 错误码。
    call MPI_ALLREDUCE(error1, error1Global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, COMM2D, IERR)
    call MPI_ALLREDUCE(error2, error2Global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, COMM2D, IERR)
    call MPI_ALLREDUCE(error5, error5Global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, COMM2D, IERR)
    call MPI_ALLREDUCE(error6, error6Global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, COMM2D, IERR)

    errorU = dsqrt(error1Global)/dsqrt(error2Global)      !速度场相对L2误差：||u^n-u^{n-1}||_2 / ||u^n||_2
    errorT = error5Global/error6Global                    !温度场相对L1误差：||T^n-T^{n-1}||_1 / ||T^n||_1

  

    if(isRoot) then
        call append_convergence_tecplot('convergence.plt', restartItcOffset+itc, errorU, errorT)

    
        write(caseTag,'("Ra=",ES10.3E2,",nx=",I0,",ny=",I0,",useG=",L1,",old=",L1)') Rayleigh, nx, ny, useG,&
        &useLegacyThermalScheme  !输出收敛曲线的对比
        call append_convergence_master_tecplot('convergence_all.plt', caseTag, restartItcOffset+itc, errorU, errorT)

        write(*,'(I12,1X,ES24.16,1X,ES24.16)') restartItcOffset+itc, errorU, errorT
    endif


    return
    end subroutine check
!===================================================================================================
! check 结束: 计算稳态收敛误差并写出收敛历史。
!===================================================================================================
#endif



!===================================================================================================
! 子程序: append_convergence_tecplot
! 作用: 向单个 Tecplot 收敛文件追加一条误差记录。
! 用途: 在 check 中调用，用于输出当前算例的收敛曲线。
!===================================================================================================
subroutine append_convergence_tecplot(filename, itc, errorU, errorT)
  use commondata, only: loadInitField
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4), intent(in)  :: itc
  real(kind=8),    intent(in)  :: errorU, errorT

  integer :: u
  logical :: ex
  logical, save :: first_write = .true.      !这个局部变量在子程序返回后不会被销毁，下一次再进这个子程序时还保留上一次的值

  if (first_write) then
    inquire(file=trim(filename), exist=ex)
    if ((loadInitField.EQ.1).AND.ex) then
      ! 续算：旧收敛曲线继续追加，横坐标已经传入累计 itc。
      open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
    elseif(loadInitField.EQ.1) then
      ! 续算要求旧收敛文件必须存在；否则会丢掉断电前的收敛历史。
      write(*,*) 'Error: restart requested but convergence file is missing: ', trim(filename)
      stop
    else
      ! 新算例：清掉旧历史，避免不同算例的数据混在一起。
      open(newunit=u, file=trim(filename), status='replace', action='write', form='formatted')
      write(u,'(A)') 'VARIABLES = "itc" "errorU" "errorT"'
      write(u,'(A)') 'ZONE T="conv", F=POINT'
    endif
    write(u,'(I12,1X,ES24.16E3,1X,ES24.16E3)') itc, errorU, errorT
    close(u)

    first_write = .false.
  else
    ! 同一次运行的后续调用：追加数据行
    open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
    write(u,'(I12,1X,ES24.16E3,1X,ES24.16E3)') itc, errorU, errorT
    close(u)
  end if

end subroutine append_convergence_tecplot
!===================================================================================================
! append_convergence_tecplot 结束: 向单个 Tecplot 收敛文件追加一条误差记录。
!===================================================================================================



!===================================================================================================
! 子程序: append_convergence_master_tecplot
! 作用: 向总收敛文件追加一条带算例标签的记录。
! 用途: 预留给 check 调用，用于汇总多组算例的收敛信息。
!===================================================================================================
subroutine append_convergence_master_tecplot(filename, zoneName, itc, errorU, errorT)
  use commondata, only: loadInitField
  implicit none
  character(len=*), intent(in) :: filename, zoneName
  integer(kind=4), intent(in)  :: itc
  real(kind=8),    intent(in)  :: errorU, errorT

  logical :: ex
  integer :: u
  logical, save :: zone_started = .false.

  ! 本次运行第一次写：新算例写入新 ZONE；续算则接到旧 ZONE 后面
  if (.not. zone_started) then
    inquire(file=trim(filename), exist=ex)

    if (.not. ex) then
      if(loadInitField.EQ.1) then
        write(*,*) 'Error: restart requested but master convergence file is missing: ', trim(filename)
        stop
      endif
      open(newunit=u, file=trim(filename), status='new', action='write', form='formatted')
      write(u,'(A)') 'TITLE = "Convergence comparison"'
      write(u,'(A)') 'VARIABLES = "itc" "errorU" "errorT"'
      close(u)
    end if

    if(loadInitField.EQ.0) then
      open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
      write(u,'(A)') 'ZONE T="'//trim(zoneName)//'", F=POINT'
      close(u)
    endif

    zone_started = .true.
  end if

  ! 追加一个数据点
  open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
  write(u,'(I12,1X,ES24.16E3,1X,ES24.16E3)') itc, errorU, errorT
  close(u)
end subroutine append_convergence_master_tecplot
!===================================================================================================
! append_convergence_master_tecplot 结束: 向总收敛文件追加一条带算例标签的记录。
!===================================================================================================





!===================================================================================================
! 子程序: output_SnapshotFile
! 作用: 输出 u、v、T、rho 的二进制快照文件，供后处理分析使用。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine output_SnapshotFile()                                   !输出 uvTrho 二进制快照
    use commondata                                                   !用于后处理快照；重启读入时必须按 u,v,T,rho 顺序读取
    implicit none
    integer(kind=4) :: i, j
    character(len=100) :: filename
    ! This snapshot is for post-processing only; u/v are written after nondimensionalization.
    ! For strict restart, keep using output_ReloadFile(), which preserves the lattice-state variables.
    
#ifdef steadyFlow
    write(filename,'(i12.12)') restartItcOffset+itc                  !steadyFlow：续算时使用累计格子步编号
#endif

#ifdef unsteadyFlow
    snapshotFileNum = snapshotFileNum+1
    write(filename,'(i12.12)') snapshotFileNum   !unsteadyFlow：快照文件按调用次数编号，与 reloadFileNum 分离
#endif

    filename = adjustl(filename)

    call gather_output_fields_mpi()   ! 输出 snapshot 前，先把各 rank 的局部 u/v/T/rho 拼回 root 的全局数组
    if(.not.isRoot) return

    open(unit=03,file=trim(snapshotFilePrefix)//"-"//trim(filename)//'.bin',form="unformatted",access="sequential")    !二进制
    ! Post-processing snapshot only: write nondimensionalized u/v together with T and rho.
    ! Do not use this file for strict restart; output_ReloadFile() keeps lattice velocities for that purpose.
    write(03) ((real(velocityScaleCompare*u_all(i,j),kind=8),i=1,nx),j=1,ny)
    write(03) ((real(velocityScaleCompare*v_all(i,j),kind=8),i=1,nx),j=1,ny)
    write(03) ((real(T_all(i,j),kind=8),i=1,nx),j=1,ny)
    write(03) ((real(rho_all(i,j),kind=8), i=1,nx), j=1,ny)
    close(03)

    return
  end subroutine output_SnapshotFile
!===================================================================================================
! output_SnapshotFile 结束: 输出 u、v、T、rho 的二进制快照文件。
!===================================================================================================

!===================================================================================================
! 子程序: output_SnapshotMeshFile
! 作用: 输出所有快照共用的非均匀网格坐标。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine output_SnapshotMeshFile()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    if(.not.isRoot) return
    open(unit=07,file=trim(snapshotFilePrefix)//"-mesh.dat",status="replace",action="write")
    write(07,'(A)') "# ISLBM nonuniform mesh coordinates for snapshot binary files"
    write(07,'(A,1X,I0,1X,I0)') "nx_ny", nx, ny
    write(07,'(A)') "# xp"
    do i = 1, nx
        write(07,'(I8,1X,ES24.16E3)') i, real(xp(i),kind=8)
    enddo
    write(07,'(A)') "# yp"
    do j = 1, ny
        write(07,'(I8,1X,ES24.16E3)') j, real(yp(j),kind=8)
    enddo
    close(07)

    return
  end subroutine output_SnapshotMeshFile
!===================================================================================================
! output_SnapshotMeshFile 结束: 已输出所有快照共用的非均匀坐标。
!===================================================================================================


    

!===================================================================================================
! 子程序: output_ReloadFile
! 作用: 输出重启备份文件；基础记录为 f/g，EnableUseG 还包含 Bx_prev/By_prev 历史项。
! 用途: 在运行过程中定期调用，也在程序结束前调用。
!===================================================================================================
  subroutine output_ReloadFile()                                  !输出重启文件，名字由 reloadFilePrefix 控制
    use commondata                                                !用于严格重启；rho/u/v/T 由 f/g 重建
    implicit none
    integer(kind=4) :: i, j, alpha
    character(len=100) :: filename

#ifdef steadyFlow
    reloadFileNum = restartItcOffset+itc
    write(filename,'(i12.12)') reloadFileNum                 !steadyFlow：重启文件名使用累计格子步
#endif

#ifdef unsteadyFlow
    reloadFileNum = reloadFileNum + 1
    write(filename,'(i12.12)') reloadFileNum     !unsteadyFlow：reload 文件使用独立编号，不依赖快照输出是否开启
#endif

    filename = adjustl(filename)

    call gather_restart_fields_mpi()
    if(.not.isRoot) return

    open(unit=05,file=trim(reloadFilePrefix)//"-"//trim(filename)//'.bin',form="unformatted",access="sequential")   !二进制
    ! Strict restart files store f/g.  With EnableUseG, Bx_prev/By_prev must also be saved;
    ! otherwise the first post-reload M1G correction would lose its previous-step history.
    ! 与 read_reload_fields_mpi 保持一致：数组仍是 (i,j,alpha)，但写文件时按 Fortran 连续内存方向让 i 最内层。
    write(05) (((real(f_all(i,j,alpha),kind=8), i=1,nx), j=1,ny), alpha=0,8)
    write(05) (((real(g_all(i,j,alpha),kind=8), i=1,nx), j=1,ny), alpha=0,4)
#ifdef EnableUseG
    write(05) ((real(Bx_prev_all(i,j),kind=8), i=1,nx), j=1,ny)
    write(05) ((real(By_prev_all(i,j),kind=8), i=1,nx), j=1,ny)
#endif
    close(05)
    call write_reload_metadata(trim(filename))
    
    return
  end subroutine output_ReloadFile
!===================================================================================================
! output_ReloadFile 结束: 输出严格重启备份文件。
!===================================================================================================


!===================================================================================================
! 子程序: write_reload_metadata
! 作用: 覆盖写出最新 reload 续算账本，恢复累计步数、t_ff、输出编号和最新 .bin 文件名。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine write_reload_metadata(filename)
    use commondata
    implicit none
    character(len=*), intent(in) :: filename
    integer(kind=4) :: metaUnit, totalItc
    real(kind=8) :: totalTf

    totalItc = restartItcOffset + itc
    totalTf = real(totalItc,kind=8) / timeUnit

    open(newunit=metaUnit, file=trim(reloadFilePrefix)//'-latest.meta', &
         status='replace', action='write', form='formatted')
    write(metaUnit,'(A,1X,I0)') 'reload_meta_version', 2
#ifdef steadyFlow
    write(metaUnit,'(A,1X,A)') 'flowMode', 'steadyFlow'
#endif
#ifdef unsteadyFlow
    write(metaUnit,'(A,1X,A)') 'flowMode', 'unsteadyFlow'
#endif
    write(metaUnit,'(A,1X,I0)') 'nx', nx
    write(metaUnit,'(A,1X,I0)') 'ny', ny
    write(metaUnit,'(A,1X,A)') 'reloadFileName', trim(filename)
    write(metaUnit,'(A,1X,I0)') 'itc_total', totalItc
    write(metaUnit,'(A,1X,ES24.16E3)') 'time_tf', totalTf
    write(metaUnit,'(A,1X,I0)') 'snapshotFileNum', snapshotFileNum
    write(metaUnit,'(A,1X,I0)') 'pltFileNum', pltFileNum
    write(metaUnit,'(A,1X,I0)') 'reloadFileNum', reloadFileNum
    close(metaUnit)

    return
  end subroutine write_reload_metadata
!===================================================================================================
! write_reload_metadata 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
!===================================================================================================


!===================================================================================================
! 子程序: read_reload_metadata
! 作用: 优先读取 latest .meta；若没有，则根据手工编号做保守推断。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine read_reload_metadata(reloadFileName)
    use commondata
    implicit none
    character(len=*), intent(inout) :: reloadFileName
    character(len=64) :: label
    character(len=32) :: metaFlowMode, currentFlowMode
    character(len=100) :: metaReloadFileName
    character(len=256) :: metaFile
    integer(kind=4) :: metaUnit, ios
    integer(kind=4) :: metaVersion, metaNx, metaNy
    integer(kind=4) :: metaItc, metaSnapshotFileNum, metaPltFileNum, metaReloadFileNum
    real(kind=8) :: metaTf
    logical :: metaExists

    reloadMetadataLoaded = .false.
    metaFile = trim(reloadFilePrefix)//'-latest.meta'
    inquire(file=trim(metaFile), exist=metaExists)                 !优先检查最新账本

    if(.not.metaExists) then                                       !latest meta 不存在时，只能保守推断
        call infer_reload_offsets_without_metadata()
        return                                                     !退出当前这个子程序
    endif

    open(newunit=metaUnit, file=trim(metaFile), status='old', action='read', &
         form='formatted', iostat=ios)                             !ios == 0表示读写成功,不等于0，则代表失败
    if(ios.NE.0) then
        write(*,*) 'Error: failed to open reload metadata: ', trim(metaFile)
        stop
    endif

    read(metaUnit,*,iostat=ios) label, metaVersion
    if((ios.NE.0).OR.(trim(label).NE.'reload_meta_version').OR.(metaVersion.NE.2)) then
        write(*,*) 'Error: invalid reload metadata version in ', trim(metaFile)
        stop    !如果读取失败，或者标签不是 reload_meta_version，或者版本号不是 2，就说明文件格式不对，停止。
    endif

    read(metaUnit,*,iostat=ios) label, metaFlowMode
    if((ios.NE.0).OR.(trim(label).NE.'flowMode')) then
        write(*,*) 'Error: invalid flowMode entry in ', trim(metaFile)
        stop
    endif

    read(metaUnit,*,iostat=ios) label, metaNx
    if((ios.NE.0).OR.(trim(label).NE.'nx')) then
        write(*,*) 'Error: invalid nx entry in ', trim(metaFile)
        stop
    endif

    read(metaUnit,*,iostat=ios) label, metaNy
    if((ios.NE.0).OR.(trim(label).NE.'ny')) then
        write(*,*) 'Error: invalid ny entry in ', trim(metaFile)
        stop
    endif

    read(metaUnit,*,iostat=ios) label, metaReloadFileName
    if((ios.NE.0).OR.(trim(label).NE.'reloadFileName')) then
        write(*,*) 'Error: invalid reloadFileName entry in ', trim(metaFile)
        stop
    endif
    metaReloadFileName = adjustl(metaReloadFileName)   ! 处理左边空格，把内容左对齐
    
    read(metaUnit,*,iostat=ios) label, metaItc
    if((ios.NE.0).OR.(trim(label).NE.'itc_total')) then
        write(*,*) 'Error: invalid itc_total entry in ', trim(metaFile)
        stop
    endif

    read(metaUnit,*,iostat=ios) label, metaTf
    if((ios.NE.0).OR.(trim(label).NE.'time_tf')) then
        write(*,*) 'Error: invalid time_tf entry in ', trim(metaFile)
        stop
    endif

    read(metaUnit,*,iostat=ios) label, metaSnapshotFileNum
    if((ios.NE.0).OR.(trim(label).NE.'snapshotFileNum')) then
        write(*,*) 'Error: invalid snapshotFileNum entry in ', trim(metaFile)
        stop
    endif

    read(metaUnit,*,iostat=ios) label, metaPltFileNum
    if((ios.NE.0).OR.(trim(label).NE.'pltFileNum')) then
        write(*,*) 'Error: invalid pltFileNum entry in ', trim(metaFile)
        stop
    endif

    read(metaUnit,*,iostat=ios) label, metaReloadFileNum
    if((ios.NE.0).OR.(trim(label).NE.'reloadFileNum')) then
        write(*,*) 'Error: invalid reloadFileNum entry in ', trim(metaFile)
        stop
    endif
    close(metaUnit)

    currentFlowMode = 'unknown'
#ifdef steadyFlow
    currentFlowMode = 'steadyFlow'
#endif
#ifdef unsteadyFlow
    currentFlowMode = 'unsteadyFlow'
#endif

    if(trim(metaFlowMode).NE.trim(currentFlowMode)) then
        write(*,*) 'Error: reload metadata flowMode differs: ', trim(metaFlowMode), trim(currentFlowMode)
        stop
    endif
    if((metaNx.NE.nx).OR.(metaNy.NE.ny)) then
        write(*,*) 'Error: reload metadata mesh mismatch: ', metaNx, metaNy, nx, ny
        stop
    endif

    restartItcOffset = metaItc
    reloadDimensionlessTime = metaTf
    snapshotFileNum = metaSnapshotFileNum
    pltFileNum = metaPltFileNum
    ! reloadFileNum 是整数计数器，给后续 output_ReloadFile() 继续编号，避免覆盖旧 reload 文件。
    ! reloadFileName 是字符串文件名，本次续算马上用它打开 <reloadFilePrefix>-<reloadFileName>.bin。
    reloadFileNum = metaReloadFileNum
    reloadFileName = trim(metaReloadFileName)     ! 处理右边空格，删掉尾部空格
    reloadMetadataLoaded = .true.

    return
  end subroutine read_reload_metadata
!===================================================================================================
! read_reload_metadata 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
!===================================================================================================


!===================================================================================================
! 子程序: infer_reload_offsets_without_metadata
! 作用: 没有 latest .meta 时，只能根据文件编号和当前手工参数推断。
! 根据文件名编号和当前参数“猜一个合理值”，保证续算的时间/步数是连续的。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine infer_reload_offsets_without_metadata()
    use commondata
    implicit none

    restartItcOffset = 0
#ifdef steadyFlow
    restartItcOffset = max(0, reloadFileNum)  !稳态的 reload 文件名本来就是用 itc 写的
    if(reloadDimensionlessTime.EQ.0.0d0) then !如果没有手工给 reloadDimensionlessTime,计算一下
        reloadDimensionlessTime = real(restartItcOffset,kind=8) / timeUnit
    endif
#endif
#ifdef unsteadyFlow
    if(reloadDimensionlessTime.EQ.0.0d0) then   !非稳态的 reload 文件名不是 itc，而是第几个 reload 文件
        reloadDimensionlessTime = real(max(0,reloadFileNum),kind=8) * reloadFileInterval
    endif
    restartItcOffset = max(0, int(reloadDimensionlessTime*timeUnit+0.5d0))  !再反推格子步
    snapshotFileNum = max(0, int(reloadDimensionlessTime/outputSnapshotInterval+0.5d0)) !推断输出编号
    pltFileNum = max(0, int(reloadDimensionlessTime/outputPltFileInterval+0.5d0))
#endif

    return
  end subroutine infer_reload_offsets_without_metadata
!===================================================================================================
! infer_reload_offsets_without_metadata 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
!===================================================================================================


!===================================================================================================
! 子程序: output_Tecplot
! 作用: 输出主场变量到 Tecplot 文件，便于后处理和可视化。
! 用途: 在运行过程中按需调用，也在程序结束时调用。
!===================================================================================================
  subroutine output_Tecplot()                        !输出二进制文件
    use commondata
    implicit none
    ! Here u and v are exported as nondimensional post-processing velocities using velocityScaleCompare.
    ! Restart files should still come from output_ReloadFile(), which preserves the lattice-state information.
    integer(kind=4) :: i, j, k
    REAL(kind=4) :: zoneMarker, eohMarker   !Tecplot 二进制格式里用的两个“标记值”（299 和 357），用于告诉 Tecplot：这里开始是 zone 描述 / header 结束。
    character(len=40) :: title              !文件 Title 字符串
    character(len=40) :: V1,V2,V3,V4,V5     !变量名字符串（X,Y,U,V,T）
    integer(kind=4), parameter :: kmax=1    !二维数据也按 3D 的 IJK 写，K=1
    character(len=40) :: zoneName           !zone 名称
    character(len=100) :: filename          !输出文件名字符串
    ! MPI+OpenACC 版本：输出前会先从设备端同步局部场，再由 root 汇总全局数组。

#ifdef steadyFlow
    write(filename,'(i12.12)') restartItcOffset+itc
#endif

#ifdef unsteadyFlow
    pltFileNum = pltFileNum+1               !plt 文件按调用次数编号
    write(filename,'(i12.12)') pltFileNum
#endif

    filename = adjustl(filename)            !存储路径 pltFolderPrefix="./pltFile/buoyancyCavity000000000034.plt

    call gather_output_fields_mpi()   ! 输出 Tecplot 前，先把各 rank 的局部 u/v/T/rho 拼回 root 的全局数组
    if(.not.isRoot) return

    open(41,file=trim(pltFolderPrefix)//"-"//trim(filename)//'.plt', access='stream', form='unformatted')    !stream：字节流

    !---------------------------------------------
    zoneMarker= 299.0                                     !固定的，不需要修改
    eohMarker = 357.0

    !I. HEAD SECTION--------------------------------------
    !c--Magic number, Version number
    write(41) "#!TDV101"                                  !Tecplot 识别二进制 .plt 的“魔数/版本号”字符串,不需要修改

    !c--Integer value of 1
    write(41) 1                                           !Tecplot 规范里紧跟着一个整型值（通常表示字节序/文件类型版本等控制字段）,不需要修改

    Title="MyFirst"                                       !Tecplot 的二进制格式里字符串不是直接 write，而是逐字符写 ASCII 码，再以 0 结尾
    call dumpstring(title)                                !dumpstring() 就是干这个的

    !c-- Number of variables in this data file
    write(41) 5                                           !有 5 个变量：X, Y, U, V, T,如果需要修改，这个有变化需要修改

    !c-- Variable names.                                  !变量名依次写入,有变化需要修改
    V1='X'
    call dumpstring(V1)
    V2='Y'
    call dumpstring(V2)
    V3='U_nd'
    call dumpstring(V3)
    V4='V_nd'
    call dumpstring(V4)
    V5='T'
    call dumpstring(V5)

    !c-----Zones-----------------------------

    !c--------Zone marker. Value = 299.0                   !写入 float 299.0，告诉 Tecplot：一个 zone 的描述开始了。
    write(41) zoneMarker

    !--------Zone name.
    zoneName="ZONE 001"
    call dumpstring(zoneName)                              !zone 名称，显示在 Tecplot 里（比如 Zones 面板里）

    !---------Zone Color                                   !-1 表示使用默认配色
    write(41) -1

    !---------ZoneType                                     !0 通常表示 ORDERED（结构网格 IJK）
    write(41) 0

    !---------DataPacking 0=Block, 1=Point                 !0 = Block（先把整个 X 写完，再写整个 Y…）,1 = Point（每个点依次写 X,Y,U,V,T）

    write(41) 1

    !---------Specify Var Location. 0 = Do not specify, all data    !0：不指定，默认所有变量都在节点上（nodal）,1：会跟着一串变量位置列表（cell-centered/nodal）
    !---------is located at the nodes. 1 = Specify
    write(41) 0

    !---------Number of user defined face neighbor connections   !自定义邻接连接数；结构网格一般不需要，写 0
    ! (value >= 0)
    write(41) 0

    !---------IMax,JMax,KMax                          !结构网格尺寸：I=nx, J=ny, K=1, 有变化需要修改
    write(41) nx
    write(41) ny
    write(41) kmax

    !-----------1=Auxiliary name/value pair to follow          !0：没有辅助信息（Aux data）,然后写 357.0，告诉 Tecplot：头部结束，接下来是数据区描述/数据
    !-----------0=No more Auxiliar name/value pairs.
    write(41) 0
    write(41) eohMarker

    !----zone ------------------------------------------------------------
    ! 再写一次 299.0: Tecplot 规范里 zone header 之后的数据描述块会以 zone marker 开始。
    write(41) zoneMarker

    !--------variable data format: 1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
    ! 这里所有变量按双精度输出
    write(41) 2
    write(41) 2
    write(41) 2
    write(41) 2
    write(41) 2

    !--------Has variable sharing 0 = no, 1 = yes.            !0：不共享（每个变量都在这个文件里独立存）,1：共享（例如多个 zone 共享同一份 X,Y 坐标）
    write(41) 0

    !----------Zone number to share connectivity list with (-1 = no    !-1：不共享连接表（对 ordered zone 一般无连接表）
    ! sharing).
    write(41) -1

    !---------------------------------------------------------------------   !真正写数据（按 Point packing）
    do k=1,kmax
        do j=1,ny
            do i=1,nx
                write(41) real(xp(i),kind=8)
                write(41) real(yp(j),kind=8)
                write(41) real(velocityScaleCompare*u_all(i,j),kind=8)
                write(41) real(velocityScaleCompare*v_all(i,j),kind=8)
                write(41) real(T_all(i,j),kind=8)
            end do
        end do
    enddo
    close(41)
    !---------------------------------------------------------------------

    return
  end subroutine output_Tecplot
!===================================================================================================
! output_Tecplot 结束: 输出主场变量到 Tecplot 文件，便于后处理和可视化。
!===================================================================================================




!===================================================================================================
! 子程序: dumpstring
! 作用: 把字符串按 Tecplot 二进制格式写入文件。
! 用途: 在 output_Tecplot 和 output_Tecplot_psi_vort 中作为辅助写出工具调用。
!===================================================================================================
  subroutine dumpstring(instring)
    implicit none
    character(len=40) instring
    integer(kind=4) :: stringLength   !有效长度（去掉尾部空格）
    integer(kind=4) :: ii             !字符索引
    integer(kind=4) :: I              !字符的 ASCII 码整数

    stringLength=LEN_TRIM(instring)   !LEN_TRIM: 得到去掉尾部空格后的长度
    do ii=1,stringLength
        I=ICHAR(instring(ii:ii))      !ICHAR 把字符转成 ASCII 编码整数
        write(41) I                   !把这个整数写入文件（Tecplot 要求以整数序列写字符串）
    end do
    write(41) 0                       !最后写一个 0 作为字符串结束符

    return
  end subroutine dumpstring
!===================================================================================================
! dumpstring 结束: 把字符串按 Tecplot 二进制格式写入文件。
!===================================================================================================



!===================================================================================================
! 子程序: calNuRe
! 作用: 计算体平均 Nu / Re，并把时间序列缓存到数组中。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine calNuRe()
    use commondata
    implicit none
    integer(kind=4) :: i, j, globalI, globalJ
    real(kind=8) :: NuVolAvg_temp, ReVolAvg_temp, areaWeight
    real(kind=8) :: NuReLocal(2), NuReGlobal(2)
    real(kind=8) :: sampleTime
    logical :: exNu, exRe
    logical, save :: first_nure_write = .true.
    

    if (dimensionlessTime.GE.dimensionlessTimeMax) then
        write(*,*) "Error: dimensionlessTime exceeds dimensionlessTimeMax, please enlarge dimensionlessTimeMax"
        open(unit=00,file=trim(settingsFile),status="unknown",position="append")
        write(00,*) "Error: dimensionlessTime exceeds dimensionlessTimeMax, please enlarge dimensionlessTimeMax"
        close(00)
        stop
    endif

    dimensionlessTime = dimensionlessTime+1   !每隔 outputSnapshotInterval 个自由落体时间调用一次calNuRe
#ifdef steadyFlow
    sampleTime = real(restartItcOffset+itc,kind=8)
#endif
#ifdef unsteadyFlow
    sampleTime = reloadDimensionlessTime + real(dimensionlessTime,kind=8)*outputSnapshotInterval
#endif

    if((first_nure_write).AND.(loadInitField.EQ.1).AND.isRoot) then
        inquire(file="Nu_VolAvg.dat", exist=exNu)
        inquire(file="Re_VolAvg.dat", exist=exRe)
        if((.not.exNu).OR.(.not.exRe)) then
            write(*,*) "Error: restart requested but old Nu/Re time-series files are missing."
            open(unit=00,file=trim(settingsFile),status="unknown",position="append")
            write(00,*) "Error: restart requested but old Nu/Re time-series files are missing."
            write(00,*) "Nu_VolAvg.dat exists =", exNu
            write(00,*) "Re_VolAvg.dat exists =", exRe
            close(00)
            call MPI_ABORT(COMM2D, 2001, IERR)
        endif
    endif


    
    NuVolAvg_temp = 0.0d0    
#ifdef SideHeatedCell  
    !$acc wait(1)
    !$acc parallel loop collapse(2) present(u,T,quadWidthX,quadWidthY) &
    !$acc& private(globalI,globalJ,areaWeight) reduction(+:NuVolAvg_temp)
    do j = 1, yLocalCount
        do i = 1, xLocalCount
            globalI = xStartGlobal + i - 1
            globalJ = yStartGlobal + j - 1
            areaWeight = quadWidthX(globalI)*quadWidthY(globalJ)
            NuVolAvg_temp = NuVolAvg_temp+areaWeight*u(i,j)*(T(i,j)-Tref)     !非均匀网格 midpoint-rule 加权对流热通量
        enddo
    enddo
#endif

#ifdef RayleighBenardCell  
    !$acc wait(1)
    !$acc parallel loop collapse(2) present(v,T,quadWidthX,quadWidthY) &
    !$acc& private(globalI,globalJ,areaWeight) reduction(+:NuVolAvg_temp)
    do j = 1, yLocalCount
        do i = 1, xLocalCount
            globalI = xStartGlobal + i - 1
            globalJ = yStartGlobal + j - 1
            areaWeight = quadWidthX(globalI)*quadWidthY(globalJ)
            NuVolAvg_temp = NuVolAvg_temp+areaWeight*v(i,j)*(T(i,j)-Tref)     !非均匀网格 midpoint-rule 加权对流热通量
        enddo
    enddo
#endif

    ReVolAvg_temp = 0.0d0
    !$acc wait(1)
    !$acc parallel loop collapse(2) present(u,v,quadWidthX,quadWidthY) &
    !$acc& private(globalI,globalJ,areaWeight) reduction(+:ReVolAvg_temp)
    do j = 1, yLocalCount
        do i = 1, xLocalCount 
            globalI = xStartGlobal + i - 1
            globalJ = yStartGlobal + j - 1
            areaWeight = quadWidthX(globalI)*quadWidthY(globalJ)
            ReVolAvg_temp = ReVolAvg_temp+areaWeight*(u(i,j)*u(i,j)+v(i,j)*v(i,j))
        enddo
    enddo

    NuReLocal = (/ NuVolAvg_temp, ReVolAvg_temp /)
    call MPI_ALLREDUCE(NuReLocal, NuReGlobal, 2, MPI_DOUBLE_PRECISION, MPI_SUM, COMM2D, IERR)
    NuVolAvg(dimensionlessTime) = NuReGlobal(1)/quadSumArea*lengthUnit/diffusivity+1.0d0    !!ISLBM体平均 Nusselt 数：非均匀网格体积分加权
    ReVolAvg(dimensionlessTime) = dsqrt(NuReGlobal(2)/quadSumArea)*lengthUnit/viscosity

    if(isRoot) then
        if((first_nure_write).AND.(loadInitField.EQ.0)) then
            open(unit=01,file="Nu_VolAvg.dat",status='replace',action='write')
        else
            open(unit=01,file="Nu_VolAvg.dat",status='unknown',position='append',action='write')
        endif
        write(01,'(ES24.16E3,1X,ES24.16E3)') &
            real(sampleTime,kind=8), &
            real(NuVolAvg(dimensionlessTime),kind=8)   !以格子步数或者自由落体时间来写入
        close(01)
    endif

    if(isRoot) then
        if((first_nure_write).AND.(loadInitField.EQ.0)) then
            open(unit=02,file="Re_VolAvg.dat",status='replace',action='write')
        else
            open(unit=02,file="Re_VolAvg.dat",status='unknown',position='append',action='write')
        endif
        write(02,'(ES24.16E3,1X,ES24.16E3)') &
            real(sampleTime,kind=8), &
            real(ReVolAvg(dimensionlessTime),kind=8)
        close(02)
    endif
    first_nure_write = .false.
    if(isRoot) then
        write(*,'(a,1x,ES24.16E3)') "NuVolAvg =", real(NuVolAvg(dimensionlessTime),kind=8)
        write(*,'(a,1x,ES24.16E3)') "ReVolAvg =", real(ReVolAvg(dimensionlessTime),kind=8)
    endif
    return
  end subroutine calNuRe
!===================================================================================================
! calNuRe 结束: 计算体平均 Nu 和 Re 的时间历程统计量。
!===================================================================================================


#ifdef unsteadyFlow
!===================================================================================================
! 子程序: output_unsteady_NuRe_postprocess
! 作用: 从完整 Nu/Re 历史文件重建非稳态时间序列、运行平均值和分段窗口平均值。
! 用途: 非稳态算例结束后用于输出后处理统计文件。
!===================================================================================================
  subroutine output_unsteady_NuRe_postprocess()
    use commondata
    implicit none

    integer(kind=4) :: k, history_count
    integer(kind=4) :: whole_count, first_count, second_count
    integer(kind=4) :: iosNu, iosRe
    integer(kind=4) :: nuUnit, reUnit, seriesUnit, runningUnit
    real(kind=8) :: timeLoc, timeNu, timeRe, NuVal, ReVal, Nu_Accum, Re_Accum
    real(kind=8) :: startTf, midTf, endTf
    real(kind=8) :: Nu_WholeSum, Re_WholeSum, Nu_FirstSum, Re_FirstSum, Nu_SecondSum, Re_SecondSum
    real(kind=8) :: Nu_WholeAvg, Re_WholeAvg, Nu_FirstAvg, Re_FirstAvg, Nu_SecondAvg, Re_SecondAvg
    real(kind=8) :: Nu_FirstRelErr, Re_FirstRelErr, Nu_SecondRelErr, Re_SecondRelErr
    logical :: exNu, exRe

    inquire(file='Nu_VolAvg.dat', exist=exNu)
    inquire(file='Re_VolAvg.dat', exist=exRe)
    if((.not.exNu).or.(.not.exRe)) then
        write(*,'(A)') 'Error: Nu/Re history files are missing before postprocessing.'
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,'(A)') 'Error: Nu/Re history files are missing before postprocessing.'
        close(00)
        error stop 1
    endif

    open(newunit=nuUnit, file='Nu_VolAvg.dat', status='old', action='read', form='formatted')
    open(newunit=reUnit, file='Re_VolAvg.dat', status='old', action='read', form='formatted')

    ! These files are derived views of the full .dat history, so rebuild one continuous ZONE.
    open(newunit=seriesUnit, file='NuRe_VolAvg_2DOpenaccMpi.plt', status='replace', action='write', form='formatted')
    write(seriesUnit,'(A)') 'TITLE = "2D OpenACC MPI Nu/Re volume averages"'
    write(seriesUnit,'(A)') 'VARIABLES = "time" "NuVolAvg" "ReVolAvg"'
    write(seriesUnit,'(A)') 'ZONE T="NuReVolAvg", F=POINT'

    open(newunit=runningUnit, file='NuRe_VolAvg_runningMean_2DOpenaccMpi.plt', status='replace', action='write', form='formatted')
    write(runningUnit,'(A)') 'TITLE = "2D OpenACC MPI Nu/Re running means"'
    write(runningUnit,'(A)') 'VARIABLES = "time" "NuVolAvgMean" "ReVolAvgMean"'
    write(runningUnit,'(A)') 'ZONE T="NuReRunningMean", F=POINT'

    startTf = unsteadyAverageStartTf
    midTf = unsteadyAverageMidTf
    endTf = unsteadyAverageEndTf

    Nu_Accum = 0.0d0
    Re_Accum = 0.0d0
    Nu_WholeSum = 0.0d0
    Re_WholeSum = 0.0d0
    Nu_FirstSum = 0.0d0
    Re_FirstSum = 0.0d0
    Nu_SecondSum = 0.0d0
    Re_SecondSum = 0.0d0
    whole_count = 0
    first_count = 0
    second_count = 0
    history_count = 0

    do k = 1, unsteadySampleCount
        read(nuUnit,*,iostat=iosNu) timeNu, NuVal
        read(reUnit,*,iostat=iosRe) timeRe, ReVal
        ! iostat: 0=成功读到一行；小于0=到达文件末尾；大于0=格式或读入错误。
        ! 循环内只要不是 0，就说明文件短于 unsteadySampleCount，或者某一行格式坏了。
        if((iosNu.NE.0).OR.(iosRe.NE.0)) then
            write(*,'(A)') 'Error: Nu/Re history files are shorter than unsteadySampleCount or contain invalid rows.'
            open(unit=00,file=trim(settingsFile),status='unknown',position='append')
            write(00,'(A)') 'Error: Nu/Re history files are shorter than unsteadySampleCount or contain invalid rows.'
            close(00)
            error stop 1
        endif
        if(abs(timeNu-timeRe).GT.1.0d-10*max(1.0d0,abs(timeNu))) then     !确保 Nu 和 Re 是同一个时间采样点的数据，不是错行配对的数据
            write(*,'(A)') 'Error: Nu/Re history time columns do not match.'
            open(unit=00,file=trim(settingsFile),status='unknown',position='append')
            write(00,'(A)') 'Error: Nu/Re history time columns do not match.'
            close(00)
            error stop 1
        endif

        timeLoc = timeNu
        history_count = k
        Nu_Accum = Nu_Accum + NuVal
        Re_Accum = Re_Accum + ReVal
        write(seriesUnit,'(ES24.16E3,1X,ES24.16E3,1X,ES24.16E3)') timeLoc, NuVal, ReVal
        write(runningUnit,'(ES24.16E3,1X,ES24.16E3,1X,ES24.16E3)') &
            timeLoc, Nu_Accum/dble(history_count), Re_Accum/dble(history_count)

        if ((timeLoc >= startTf) .and. (timeLoc <= endTf)) then
            Nu_WholeSum = Nu_WholeSum + NuVal
            Re_WholeSum = Re_WholeSum + ReVal
            whole_count = whole_count + 1
        endif
        if ((timeLoc >= startTf) .and. (timeLoc < midTf)) then
            Nu_FirstSum = Nu_FirstSum + NuVal
            Re_FirstSum = Re_FirstSum + ReVal
            first_count = first_count + 1
        endif
        if ((timeLoc >= midTf) .and. (timeLoc <= endTf)) then
            Nu_SecondSum = Nu_SecondSum + NuVal
            Re_SecondSum = Re_SecondSum + ReVal
            second_count = second_count + 1
        endif
    enddo

    ! 上面的循环已经读完预期的 unsteadySampleCount 行；这里再试读一行，
    ! 不是为了继续计算，而是确认 Nu/Re 两个历史文件后面没有多余数据。
    read(nuUnit,*,iostat=iosNu) timeNu, NuVal
    read(reUnit,*,iostat=iosRe) timeRe, ReVal
    ! 任意一个 iostat 等于 0，都表示至少一个文件还成功读到了额外一行。
    ! 如果放过这种情况，后处理会静默丢掉超过 unsteadySampleCount 的尾部样本。
    if((iosNu.EQ.0).OR.(iosRe.EQ.0)) then
        write(*,'(A)') 'Error: Nu/Re history files contain more rows than unsteadySampleCount.'
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,'(A)') 'Error: Nu/Re history files contain more rows than unsteadySampleCount.'
        close(00)
        error stop 1
    endif
    ! 正常结尾必须是两个文件都到达 EOF，也就是两个 iostat 都小于 0。
    ! 其他组合说明一个文件尾部异常，或 Nu/Re 文件长度不一致。
    if(.not.((iosNu.LT.0).AND.(iosRe.LT.0))) then
        write(*,'(A)') 'Error: Nu/Re history files have inconsistent trailing rows.'
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,'(A)') 'Error: Nu/Re history files have inconsistent trailing rows.'
        close(00)
        error stop 1
    endif
    close(nuUnit)
    close(reUnit)
    close(seriesUnit)
    close(runningUnit)

    if (history_count <= 0) then
        write(*,'(A)') 'Error: no Nu/Re history samples were found before postprocessing.'
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,'(A)') 'Error: no Nu/Re history samples were found before postprocessing.'
        close(00)
        error stop 1
    endif

    if ((whole_count <= 0) .or. (first_count <= 0) .or. (second_count <= 0)) then
        write(*,'(A)') 'Error: no complete unsteady average window was found for Nu/Re postprocessing.'
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,'(A)') 'Error: no complete unsteady average window was found for Nu/Re postprocessing.'
        close(00)
        error stop 1
    endif

    Nu_WholeAvg = Nu_WholeSum / dble(whole_count)
    Re_WholeAvg = Re_WholeSum / dble(whole_count)
    Nu_FirstAvg = Nu_FirstSum / dble(first_count)
    Re_FirstAvg = Re_FirstSum / dble(first_count)
    Nu_SecondAvg = Nu_SecondSum / dble(second_count)
    Re_SecondAvg = Re_SecondSum / dble(second_count)
    ! tiny(1.0d0) 是双精度最小正规正数；这里用 max 避免整体平均值为 0 时相对误差除以 0。
    Nu_FirstRelErr = abs(Nu_FirstAvg - Nu_WholeAvg) / max(abs(Nu_WholeAvg), tiny(1.0d0))
    Re_FirstRelErr = abs(Re_FirstAvg - Re_WholeAvg) / max(abs(Re_WholeAvg), tiny(1.0d0))
    Nu_SecondRelErr = abs(Nu_SecondAvg - Nu_WholeAvg) / max(abs(Nu_WholeAvg), tiny(1.0d0))
    Re_SecondRelErr = abs(Re_SecondAvg - Re_WholeAvg) / max(abs(Re_WholeAvg), tiny(1.0d0))

    open(unit=33, file='NuRe_TimeAverage_2DOpenaccMpi.txt', status='replace', action='write', form='formatted')
    write(33,'(A)') '# 2D OpenACC MPI Nu/Re statistical-convergence window averages'
    write(33,'(A)') '# start_tf mid_tf end_tf whole_count first_count second_count ' // &
        'Nu_whole Re_whole Nu_first Re_first Nu_second Re_second ' // &
        'Nu_first_relerr Re_first_relerr Nu_second_relerr Re_second_relerr'
    write(33,'(ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,I0,1X,I0,1X,I0,1X,' // &
        'ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,' // &
        'ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3)') &
        startTf, midTf, endTf, whole_count, first_count, second_count, &
        Nu_WholeAvg, Re_WholeAvg, Nu_FirstAvg, Re_FirstAvg, Nu_SecondAvg, Re_SecondAvg, &
        Nu_FirstRelErr, Re_FirstRelErr, Nu_SecondRelErr, Re_SecondRelErr
    close(33)

    write(*,'(A)') 'Unsteady Nu/Re statistical postprocessing:'
    write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'Average window:', startTf, 'to', endTf
    write(*,'(A,1X,I0,1X,A,1X,I0,1X,A,1X,I0)') &
        'Samples whole/first/second:', whole_count, '/', first_count, '/', second_count
    write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'Whole Nu/Re average:', Nu_WholeAvg, '/', Re_WholeAvg
    write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'First-half Nu/Re average:', Nu_FirstAvg, '/', Re_FirstAvg
    write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'Second-half Nu/Re average:', Nu_SecondAvg, '/', Re_SecondAvg
    write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'First-half Nu/Re relative error:', Nu_FirstRelErr, '/', Re_FirstRelErr
    write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'Second-half Nu/Re relative error:', Nu_SecondRelErr, '/', Re_SecondRelErr

  end subroutine output_unsteady_NuRe_postprocess
!===================================================================================================
! output_unsteady_NuRe_postprocess 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
!===================================================================================================
#endif


!===================================================================================================
! 子程序: SideHeatedcalc_Nu_global
! 作用: 计算侧壁差温工况下的全场平均 Nusselt 数。
! 用途: 在 SideHeatedCell 工况结束后的后处理中调用。
!===================================================================================================
  subroutine SideHeatedcalc_Nu_global()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: dTdx, qx, sum_qx, areaWeight
    real(kind=8) :: xnode(3), fnode(3)
    real(kind=8) :: deltaT, coef

    deltaT = Thot - Tcold
    coef   = velocityScaleCompare

    sum_qx = 0.0d0

    !$omp parallel do default(none) shared(u_all,T_all,xp,quadWidthX,quadWidthY,coef) &
    !$omp private(i,j,dTdx,qx,areaWeight,xnode,fnode) reduction(+:sum_qx)
    do j = 1, ny
      do i = 1, nx
        if(i.EQ.1) then
          xnode = (/ xp(1), xp(2), xp(3) /)
          fnode = (/ T_all(1,j), T_all(2,j), T_all(3,j) /)
        elseif(i.EQ.nx) then
          xnode = (/ xp(nx-2), xp(nx-1), xp(nx) /)
          fnode = (/ T_all(nx-2,j), T_all(nx-1,j), T_all(nx,j) /)
        else
          xnode = (/ xp(i-1), xp(i), xp(i+1) /)
          fnode = (/ T_all(i-1,j), T_all(i,j), T_all(i+1,j) /)
        endif
        call lagrange_derivative_3(xnode, fnode, xp(i), dTdx)
        qx = coef*u_all(i,j)*(T_all(i,j)-Tref) - dTdx
        areaWeight = quadWidthX(i)*quadWidthY(j)
        sum_qx = sum_qx + areaWeight*qx

      enddo
    enddo
    !$omp end parallel do

    Nu_global = (sum_qx / quadSumArea) / deltaT

    ! 屏幕输出
    write(*,'(a,1x,ES24.16E3)') "Nu_global =", real(Nu_global,kind=8)


    ! 同步写入日志
    open(unit=00,file=trim(settingsFile),status="unknown",position="append")
    write(00,'(a,1x,ES24.16E3)') "Nu_global =", real(Nu_global,kind=8)
    close(00)

    return
  end subroutine SideHeatedcalc_Nu_global
!===================================================================================================
! SideHeatedcalc_Nu_global 结束: 计算侧壁差温工况下的全场平均 Nusselt 数。
!===================================================================================================



!===================================================================================================
! 子程序: SideHeatedcalc_Nu_wall_avg
! 作用: 计算侧壁差温工况下热壁、冷壁和中线的 Nusselt 数及其极值。
! 用途: 在 SideHeatedCell 工况结束后的后处理中调用。
!===================================================================================================
  subroutine SideHeatedcalc_Nu_wall_avg()
    use commondata
    implicit none
    integer(kind=4) :: j
    integer(kind=4) :: jmax, jmin
    integer(kind=4) :: midIndex(3)
    real(kind=8) :: deltaT,coef
    real(kind=8) :: qx_wall, sum_hot, sum_cold, sum_mid
    real(kind=8) :: u_mid, T_mid, dTdx_mid
    real(kind=8) :: xnode(3), fnode(3)
    real(kind=8) :: xmidNode(3), xmidWeight(3)
    logical :: midStencilValid
    real(kind=8) :: Nu_left(1:ny) 
    !------------------------------------------------------------
    ! 5-point least-squares parabola fit (general, can be one-sided)
    !------------------------------------------------------------
    integer(kind=4) :: k
    integer(kind=4) :: jj(5)
    real(kind=8) :: fk(5), yk(5)
    real(kind=8) :: fstar, ystar
    real(kind=8) :: Nu_left_ext(0:ny+1)


    deltaT = Thot - Tcold
    coef   = velocityScaleCompare

   

    !-----------------------------
    ! (1) 左侧热壁平均 Nu_hot，同时记录 Numax/Numin 及其 y 位置
    sum_hot = 0.0d0
    !$omp parallel do default(none) shared(T_all,Nu_left,quadWidthY,deltaT,xp) &
    !$omp private(j,qx_wall,xnode,fnode) reduction(+:sum_hot)
    do j = 1, ny
      ! 壁面导热通量：qx(x=0,j)
      xnode = (/ 0.0d0, xp(1), xp(2) /)
      fnode = (/ Thot, T_all(1,j), T_all(2,j) /)
      call lagrange_derivative_3(xnode, fnode, 0.0d0, qx_wall)
      qx_wall = -qx_wall
      Nu_left(j)= qx_wall / deltaT
      sum_hot   = sum_hot + Nu_left(j)*quadWidthY(j)
    enddo
    !$omp end parallel do
    Nu_hot = sum_hot / quadSumY

    ! 角点没有唯一的壁面法向Nu定义; 这里只把相邻壁面Nu复制到扩展数组,
    ! 用于后面的极值搜索和五点拟合, 避免重新引入均匀网格角点外推公式。
    Nu_left_ext(1:ny) = Nu_left(1:ny)
    Nu_left_ext(0) = Nu_left(1)
    Nu_left_ext(ny+1) = Nu_left(ny)

    ! 网格上先找最大/最小
    jmax = 0
    jmin = 0
    Nu_hot_max = Nu_left_ext(0)
    Nu_hot_min = Nu_left_ext(0)

    do j = 1, ny+1
      if (Nu_left_ext(j) > Nu_hot_max) then
        Nu_hot_max = Nu_left_ext(j)
        jmax = j
      endif
      if (Nu_left_ext(j) < Nu_hot_min) then
        Nu_hot_min = Nu_left_ext(j)
        jmin = j
      endif
    enddo


    !==========================
    ! 五点插值 least-squares parabola fit, else one-sided 5pts
    !==========================
        ! ----------- 选 Numax 的 5 个点 -----------
    if (jmax <= 2) then
      jj = (/ 0, 1, 2, 3, 4 /)
    elseif (jmax >= ny-1) then
      jj = (/ ny-3, ny-2, ny-1, ny, ny+1 /)
    else
      jj = (/ jmax-2, jmax-1, jmax, jmax+1, jmax+2 /)
    endif

    do k = 1, 5
      yk(k) = yp(jj(k))            ! yp(0)=0, yp(ny+1)=ny/lengthUnit 已经存在
      fk(k) = Nu_left_ext(jj(k))
    enddo

    call fit_parabola_ls5(yk, fk, +1, fstar, ystar)  !最小二乘法，用抛物线，拟合五个点

    Nu_hot_max = fstar
    Nu_hot_max_position = ystar



    ! ----------- 选 Numin 的 5 个点：单侧（before the minimum）-----------
    if (jmin >= 4) then
      jj = (/ jmin-4, jmin-3, jmin-2, jmin-1, jmin /)
    else
      jj = (/ 0, 1, 2, 3, 4 /)
    endif

    do k = 1, 5
      yk(k) = yp(jj(k))
      fk(k) = Nu_left_ext(jj(k))
    enddo

    call fit_parabola_ls5(yk, fk, -1, fstar, ystar)

    Nu_hot_min = fstar
    Nu_hot_min_position = ystar


    !-----------------------------
    ! (2) 右侧冷壁平均 Nu_cold
    sum_cold = 0.0d0
    !$omp parallel do default(none) shared(T_all,quadWidthY,deltaT,xp) &
    !$omp private(j,qx_wall,xnode,fnode) reduction(+:sum_cold)
    do j = 1, ny
      xnode = (/ xp(nx-1), xp(nx), 1.0d0 /)
      fnode = (/ T_all(nx-1,j), T_all(nx,j), Tcold /)
      call lagrange_derivative_3(xnode, fnode, 1.0d0, qx_wall)
      qx_wall = -qx_wall
      sum_cold = sum_cold + qx_wall/deltaT*quadWidthY(j)
    enddo
    !$omp end parallel do
    Nu_cold = sum_cold / quadSumY

    !-----------------------------
    ! (3) 竖直中线 x=1/2 的平均 Nu_middle
    sum_mid = 0.0d0

    call build_lagrange_stencil_1d(nx, xp(1:nx), 0.5d0, midIndex, xmidWeight, midStencilValid)
    if(.not.midStencilValid) then
      write(*,*) "Error: x=0.5 is outside ISLBM x nodes in SideHeatedcalc_Nu_wall_avg"
      stop
    endif
    xmidNode = (/ xp(midIndex(1)), xp(midIndex(2)), xp(midIndex(3)) /)

    !$omp parallel do default(none) shared(u_all,T_all,quadWidthY,deltaT,coef,midIndex,xmidNode,xmidWeight) &
    !$omp private(j,u_mid,T_mid,dTdx_mid,fnode) reduction(+:sum_mid)
    do j = 1, ny
      u_mid = xmidWeight(1)*u_all(midIndex(1),j) + xmidWeight(2)*u_all(midIndex(2),j) + &
        xmidWeight(3)*u_all(midIndex(3),j)
      T_mid = xmidWeight(1)*T_all(midIndex(1),j) + xmidWeight(2)*T_all(midIndex(2),j) + &
        xmidWeight(3)*T_all(midIndex(3),j)
      fnode = (/ T_all(midIndex(1),j), T_all(midIndex(2),j), T_all(midIndex(3),j) /)
      call lagrange_derivative_3(xmidNode, fnode, 0.5d0, dTdx_mid)
      sum_mid = sum_mid + (coef*u_mid*(T_mid-Tref) - dTdx_mid)/deltaT*quadWidthY(j)
    enddo
    !$omp end parallel do

    Nu_middle = sum_mid / quadSumY

    !-----------------------------
    ! 输出：屏幕 + 日志
    write(*,'(a,1x,es16.8)') "Nu_hot    =", Nu_hot
    write(*,'(a,1x,es16.8)') "Nu_cold   =", Nu_cold
    write(*,'(a,1x,es16.8)') "Nu_middle =", Nu_middle
    write(*,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_max =", Nu_hot_max, "y_max =", Nu_hot_max_position
    write(*,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_min =", Nu_hot_min, "y_min =", Nu_hot_min_position

    !open(unit=12,file="Nu_wall_avg.txt",status="unknown",position="append")
    ! write(12,...) 可用于输出稳态 Nu 极值调试信息。
    !close(12)

    open(unit=00,file=trim(settingsFile),status="unknown",position="append")
    write(00,'(a,1x,es16.8)') "Nu_hot    =", Nu_hot
    write(00,'(a,1x,es16.8)') "Nu_cold   =", Nu_cold
    write(00,'(a,1x,es16.8)') "Nu_middle =", Nu_middle
    write(00,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_max =", Nu_hot_max, "y_max =", Nu_hot_max_position
    write(00,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_min =", Nu_hot_min, "y_min =", Nu_hot_min_position
    close(00)

    return
  end subroutine SideHeatedcalc_Nu_wall_avg
!===================================================================================================
! SideHeatedcalc_Nu_wall_avg 结束: 计算侧壁差温工况下热壁、冷壁和中线的 Nusselt 数及其极值。
!===================================================================================================



!===================================================================================================
! 子程序: fit_adiabatic_wall_T4
! 作用: 用四点拟合估计绝热壁面的壁温。
! 用途: 在 SideHeated 和 RB 壁面 Nusselt 后处理中调用。
!===================================================================================================
  subroutine fit_adiabatic_wall_T4(y0, y, tt, T_wall)
    implicit none
    real(kind=8), intent(in)  :: y0
    real(kind=8), intent(in)  :: y(4), tt(4)
    real(kind=8), intent(out) :: T_wall
    real(kind=8) :: s(4)
    real(kind=8) :: S0, S1, S2, B0, B1, D
    integer(kind=4) :: k

    do k = 1, 4
      s(k) = (y(k) - y0)*(y(k) - y0)
    enddo

    S0 = 4.0d0
    S1 = 0.0d0
    S2 = 0.0d0
    B0 = 0.0d0
    B1 = 0.0d0
    do k = 1, 4
      S1 = S1 + s(k)
      S2 = S2 + s(k)*s(k)
      B0 = B0 + tt(k)
      B1 = B1 + tt(k)*s(k)
    enddo

    D = S0*S2 - S1*S1

    ! T_wall = c0
    T_wall = (B0*S2 - B1*S1) / D

    return
  end subroutine fit_adiabatic_wall_T4
!===================================================================================================
! fit_adiabatic_wall_T4 结束: 用四点拟合估计绝热壁面的壁温。
!===================================================================================================


!===================================================================================================
! 子程序: fit_parabola_ls5
! 作用: 用五点最小二乘抛物线拟合局部极值和对应位置。
! 用途: 在 Nu 极值和中心线速度极值的后处理中重复调用。
!===================================================================================================
  subroutine fit_parabola_ls5(y, f, mode, fstar, ystar)
    implicit none
    real(kind=8), intent(in)  :: y(5), f(5)
    real(kind=8), intent(out) :: fstar, ystar
    real(kind=8) :: S0, S1, S2, S3, S4
    real(kind=8) :: F0, F1, F2
    real(kind=8) :: D, DA, DB, DC
    real(kind=8) :: A, B, C
    real(kind=8) :: ymin, ymax
    integer(kind=4) :: k, kbest
    real(kind=8), parameter :: epsD=1.0d-20, epsA=1.0d-14
    integer(kind=4), intent(in) :: mode   ! +1 => max, -1 => min


  ! ----- fallback: pick max/min among the 5 samples -----
  kbest = 1
  do k = 2, 5
    if (mode == 1) then
      if (f(k) > f(kbest)) kbest = k
    else
      if (f(k) < f(kbest)) kbest = k
    endif
  enddo

    S0=0.0d0; S1=0.0d0; S2=0.0d0; S3=0.0d0; S4=0.0d0
    F0=0.0d0; F1=0.0d0; F2=0.0d0
    do k = 1, 5
      S0 = S0 + 1.0d0
      S1 = S1 + y(k)
      S2 = S2 + y(k)*y(k)
      S3 = S3 + y(k)*y(k)*y(k)
      S4 = S4 + y(k)*y(k)*y(k)*y(k)

      F0 = F0 + f(k)
      F1 = F1 + y(k)*f(k)
      F2 = F2 + y(k)*y(k)*f(k)
    enddo

    ! Solve:
    ! [S4 S3 S2][A] = [F2]
    ! [S3 S2 S1][B]   [F1]
    ! [S2 S1 S0][C]   [F0]
    D  =  S4*(S2*S0 - S1*S1) - S3*(S3*S0 - S1*S2) + S2*(S3*S1 - S2*S2)
    DA =  F2*(S2*S0 - S1*S1) - S3*(F1*S0 - S1*F0) + S2*(F1*S1 - S2*F0)
    DB =  S4*(F1*S0 - S1*F0) - F2*(S3*S0 - S1*S2) + S2*(S3*F0 - F1*S2)
    DC =  S4*(S2*F0 - F1*S1) - S3*(S3*F0 - F1*S2) + F2*(S3*S1 - S2*S2)

    if (dabs(D) > epsD) then
      A = DA / D
      B = DB / D
      C = DC / D

      if (dabs(A) > epsA) then
        ystar = -B / (2.0d0*A)
        ymin = minval(y);  ymax = maxval(y)

        if (ystar >= ymin .and. ystar <= ymax) then
          fstar = C - B*B/(4.0d0*A)
        else
          ! 顶点落在拟合区间外：退回到 5 点中最大的/最小的那个并不严格，
          ! 这里为了保持简单，通常会先选极值点
          ystar = y(kbest)
          fstar = f(kbest)
        endif
      else
        ystar = y(kbest)
        fstar = f(kbest)
      endif
    else
      ystar = y(kbest)
      fstar = f(kbest)
    endif

    return
  end subroutine fit_parabola_ls5
!===================================================================================================
! fit_parabola_ls5 结束: 用五点最小二乘抛物线拟合局部极值和对应位置。
!===================================================================================================



!===================================================================================================
! 子程序: SideHeatedcalc_umid_max
! 作用: 计算侧壁差温工况下中心线水平速度的最大值及位置。
! 用途: 在 SideHeatedCell 工况结束后的后处理中调用。
!===================================================================================================
  subroutine SideHeatedcalc_umid_max()
    use commondata
    implicit none
    integer(kind=4) :: j, k
    integer(kind=4) :: midIndex(3)
    integer(kind=4) :: j0
    integer(kind=4) :: jj(5)
    real(kind=8) :: uline(1:ny)
    real(kind=8) :: s(5), fu(5)
    real(kind=8) :: umax_grid, umax_fit, y_fit
    real(kind=8) :: xmid, xmidWeight(3)
    character(len=24) :: ctime, string
    integer(kind=4) :: time
    real(kind=8) :: coef
    logical :: midStencilValid

    coef = velocityScaleCompare

    ! ---- (1) 构造中线剖面 u(x=1/2, y_j) ----
    xmid = 0.5d0
    call build_lagrange_stencil_1d(nx, xp(1:nx), xmid, midIndex, xmidWeight, midStencilValid)
    if(.not.midStencilValid) then
      write(*,*) "Error: x=0.5 is outside ISLBM x nodes in SideHeatedcalc_umid_max"
      stop
    endif
    do j = 1, ny
      uline(j) = xmidWeight(1)*u_all(midIndex(1),j) + xmidWeight(2)*u_all(midIndex(2),j) + &
        xmidWeight(3)*u_all(midIndex(3),j)
    enddo

    ! ---- (2) 先找网格最大值点 j0 ----
    j0 = 1
    umax_grid = uline(1)
    do j = 2, ny
      if (uline(j) > umax_grid) then
        umax_grid = uline(j)
        j0 = j
      endif
    enddo

    ! ---- (3) 取 5 点（尽量对称；靠近端点时偏侧） ----
    if (j0 <= 2) then
      jj = (/ 1, 2, 3, 4, 5 /)
    elseif (j0 >= ny-1) then
      jj = (/ ny-4, ny-3, ny-2, ny-1, ny /)
    else
      jj = (/ j0-2, j0-1, j0, j0+1, j0+2 /)
    endif

    do k = 1, 5
      s(k) = yp(jj(k))       ! 自变量：y
     fu(k) = uline(jj(k))    ! 拟合量：u
    enddo

    call fit_parabola_ls5(s, fu, +1, umax_fit, y_fit)

    ! ---- 输出 ----
    write(*,'(A,1X,F12.6,1X,A,1X,F12.6,1X,A,1X,F12.6)') &
         'u_mid_max =', umax_fit*coef, 'at y =', y_fit, 'on x_mid =', xmid
         

    open(unit=00,file=trim(settingsFile),status='unknown',position='append')
    string = ctime( time() )
    write(00,*) '--- calc_umid_max --- ', string
    write(00,*) 'x_mid =', xmid
    write(00,*) 'u_mid_max =', umax_fit*coef, ' y_pos =', y_fit, ' (grid j0=', j0, ')'
    close(00)

    return
  end subroutine SideHeatedcalc_umid_max
!===================================================================================================
! SideHeatedcalc_umid_max 结束: 计算侧壁差温工况下中心线水平速度的最大值及位置。
!===================================================================================================




!===================================================================================================
! 子程序: SideHeatedcalc_vmid_max
! 作用: 计算侧壁差温工况下中心线竖直速度的最大值及位置。
! 用途: 在 SideHeatedCell 工况结束后的后处理中调用。
!===================================================================================================
  subroutine SideHeatedcalc_vmid_max()
    use commondata
    implicit none
    integer(kind=4) :: i, k
    integer(kind=4) :: midIndex(3)
    integer(kind=4) :: i0
    integer(kind=4) :: ii(5)
    real(kind=8) :: vline(1:nx)
    real(kind=8) :: s(5), fv(5)
    real(kind=8) :: vmax_grid, vmax_fit, x_fit
    real(kind=8) :: ymid, ymidWeight(3)
    character(len=24) :: ctime, string
    integer(kind=4) :: time
    real(kind=8) :: coef
    logical :: midStencilValid

    coef = velocityScaleCompare

    ! ---- (1) 构造中线剖面 v(x_i, y=1/2) ----
    ymid = 0.5d0
    call build_lagrange_stencil_1d(ny, yp(1:ny), ymid, midIndex, ymidWeight, midStencilValid)
    if(.not.midStencilValid) then
      write(*,*) "Error: y=0.5 is outside ISLBM y nodes in SideHeatedcalc_vmid_max"
      stop
    endif
    do i = 1, nx
      vline(i) = ymidWeight(1)*v_all(i,midIndex(1)) + ymidWeight(2)*v_all(i,midIndex(2)) + &
        ymidWeight(3)*v_all(i,midIndex(3))
    enddo

    ! ---- (2) 先找网格最大值点 i0 ----
    i0 = 1
    vmax_grid = vline(1)
    do i = 2, nx
      if (vline(i) > vmax_grid) then
        vmax_grid = vline(i)
        i0 = i
      endif
    enddo

    ! ---- (3) 取 5 点（尽量对称；靠近端点时偏侧） ----
    if (i0 <= 2) then
      ii = (/ 1, 2, 3, 4, 5 /)
    elseif (i0 >= nx-1) then
      ii = (/ nx-4, nx-3, nx-2, nx-1, nx /)
    else
      ii = (/ i0-2, i0-1, i0, i0+1, i0+2 /)
    endif

    do k = 1, 5
      s(k) = xp(ii(k))       ! 自变量：x
     fv(k) = vline(ii(k))    ! 拟合量：v
    enddo

    call fit_parabola_ls5(s, fv, +1, vmax_fit, x_fit)

    ! ---- 输出 ----
    write(*,'(A,1X,F12.6,1X,A,1X,F12.6,1X,A,1X,F12.6)') &
         'v_mid_max =', vmax_fit*coef, 'at x =', x_fit, 'on y_mid =', ymid


    open(unit=00,file=trim(settingsFile),status='unknown',position='append')
    string = ctime( time() )
    write(00,*) '--- calc_vmid_max --- ', string
    write(00,*) 'y_mid =', ymid
    write(00,*) 'v_mid_max =', vmax_fit*coef, ' x_pos =', x_fit, ' (grid i0=', i0, ')'
    close(00)

    return
  end subroutine SideHeatedcalc_vmid_max
!===================================================================================================
! SideHeatedcalc_vmid_max 结束: 计算侧壁差温工况下中心线竖直速度的最大值及位置。
!===================================================================================================




!===================================================================================================
! 子程序: RBcalc_Nu_global
! 作用: 计算 Rayleigh-Benard 工况下的全场平均 Nusselt 数。
! 用途: 在 RayleighBenardCell 工况结束后的后处理中调用。
!===================================================================================================
subroutine RBcalc_Nu_global()
  use commondata
  implicit none
  integer(kind=4) :: i, j
  real(kind=8) :: dTdy, qy, sum_qy, areaWeight
  real(kind=8) :: ynode(3), fnode(3)
  real(kind=8) :: deltaT, coef

  deltaT = Thot - Tcold
  coef   = velocityScaleCompare

  sum_qy = 0.0d0

  !$omp parallel do default(none) &
  !$omp shared(v_all,T_all,yp,quadWidthX,quadWidthY,coef) &
  !$omp private(i,j,dTdy,qy,areaWeight,ynode,fnode) reduction(+:sum_qy)
  do j = 1, ny
    do i = 1, nx

      if(j.EQ.1) then
        ynode = (/ yp(1), yp(2), yp(3) /)
        fnode = (/ T_all(i,1), T_all(i,2), T_all(i,3) /)
      elseif(j.EQ.ny) then
        ynode = (/ yp(ny-2), yp(ny-1), yp(ny) /)
        fnode = (/ T_all(i,ny-2), T_all(i,ny-1), T_all(i,ny) /)
      else
        ynode = (/ yp(j-1), yp(j), yp(j+1) /)
        fnode = (/ T_all(i,j-1), T_all(i,j), T_all(i,j+1) /)
      endif
      call lagrange_derivative_3(ynode, fnode, yp(j), dTdy)

      qy = coef * v_all(i,j) * (T_all(i,j) - Tref) - dTdy
      areaWeight = quadWidthX(i)*quadWidthY(j)
      sum_qy = sum_qy + areaWeight*qy

    enddo
  enddo
  !$omp end parallel do

  Nu_global = (sum_qy / quadSumArea) / deltaT

  write(*,'(a,1x,ES24.16E3)') "Nu_global =", real(Nu_global,kind=8)
  open(unit=00,file=trim(settingsFile),status="unknown",position="append")
  write(00,'(a,1x,ES24.16E3)') "Nu_global =", real(Nu_global,kind=8)
  close(00)

  return
end subroutine RBcalc_Nu_global
!===================================================================================================
! RBcalc_Nu_global 结束: 计算 Rayleigh-Benard 工况下的全场平均 Nusselt 数。
!===================================================================================================


!===================================================================================================
! 子程序: RBcalc_Nu_wall_avg
! 作用: 计算 Rayleigh-Benard 工况下热壁、冷壁和中线的 Nusselt 数及其极值。
! 用途: 在 RayleighBenardCell 工况结束后的后处理中调用。
!===================================================================================================
subroutine RBcalc_Nu_wall_avg()
  use commondata
  implicit none
  integer(kind=4) :: i, k
  integer(kind=4) :: imax, imin
  integer(kind=4) :: ii(5)
  integer(kind=4) :: midIndex(3)
  real(kind=8) :: deltaT
  real(kind=8) :: qy_wall, sum_hot, sum_cold, sum_mid, coef
  real(kind=8), dimension(1:nx) :: Nu_bot
  real(kind=8), dimension(0:nx+1) :: Nu_bot_ext
  real(kind=8) :: xfit(4), Tfit(4), T_wl, T_wr, T_wl2, T_wr2
  real(kind=8) :: xk(5), fk(5)
  real(kind=8) :: fstar, xstar
  real(kind=8) :: v_mid, T_mid, dTdy_mid
  real(kind=8) :: ynode(3), fnode(3), ymidNode(3), ymidWeight(3)
  logical :: midStencilValid


  deltaT = Thot - Tcold
  coef   = velocityScaleCompare

  !-----------------------------
  ! (1) 底部热壁平均 Nu_hot（不含角点）
  sum_hot = 0.0d0
  !$omp parallel do default(none) shared(T_all,Nu_bot,quadWidthX,deltaT,yp) &
  !$omp private(i,qy_wall,ynode,fnode) reduction(+:sum_hot)
  do i = 1, nx
    ynode = (/ 0.0d0, yp(1), yp(2) /)
    fnode = (/ Thot, T_all(i,1), T_all(i,2) /)
    call lagrange_derivative_3(ynode, fnode, 0.0d0, qy_wall)
    qy_wall   = -qy_wall
    Nu_bot(i)= qy_wall / deltaT
    sum_hot  = sum_hot + Nu_bot(i)*quadWidthX(i)
  enddo
  !$omp end parallel do
  Nu_hot = sum_hot / quadSumX

  !-----------------------------
  ! (1.1) 角点扩展：用侧壁绝热（Neumann）4点拟合得到 x=0 与 x=1 处近壁温度
  ! 左下角附近：i=1..4, j=1
  xfit(1)=xp(1);  Tfit(1)=T_all(1,1)
  xfit(2)=xp(2);  Tfit(2)=T_all(2,1)
  xfit(3)=xp(3);  Tfit(3)=T_all(3,1)
  xfit(4)=xp(4);  Tfit(4)=T_all(4,1)
  call fit_adiabatic_wall_T4(0.0d0, xfit, Tfit, T_wl)   ! 估计 T(x=0, y=yp(1))
  Tfit(1)=T_all(1,2)
  Tfit(2)=T_all(2,2)
  Tfit(3)=T_all(3,2)
  Tfit(4)=T_all(4,2)
  call fit_adiabatic_wall_T4(0.0d0, xfit, Tfit, T_wl2)  ! 估计 T(x=0, y=yp(2))

  ! 右下角附近：i=nx-3..nx, j=1
  xfit(1)=xp(nx-3);  Tfit(1)=T_all(nx-3,1)
  xfit(2)=xp(nx-2);  Tfit(2)=T_all(nx-2,1)
  xfit(3)=xp(nx-1);  Tfit(3)=T_all(nx-1,1)
  xfit(4)=xp(nx  );  Tfit(4)=T_all(nx  ,1)
  call fit_adiabatic_wall_T4(xp(nx+1), xfit, Tfit, T_wr)   ! 估计 T(x=xp(nx+1), y=yp(1))
  Tfit(1)=T_all(nx-3,2)
  Tfit(2)=T_all(nx-2,2)
  Tfit(3)=T_all(nx-1,2)
  Tfit(4)=T_all(nx  ,2)
  call fit_adiabatic_wall_T4(xp(nx+1), xfit, Tfit, T_wr2)  ! 估计 T(x=xp(nx+1), y=yp(2))

  ! 组装扩展数组：角点只用于找 max/min 与拟合
  Nu_bot_ext(1:nx) = Nu_bot(1:nx)
  ynode = (/ 0.0d0, yp(1), yp(2) /)
  fnode = (/ Thot, T_wl, T_wl2 /)
  call lagrange_derivative_3(ynode, fnode, 0.0d0, qy_wall)
  Nu_bot_ext(0) = -qy_wall / deltaT
  fnode = (/ Thot, T_wr, T_wr2 /)
  call lagrange_derivative_3(ynode, fnode, 0.0d0, qy_wall)
  Nu_bot_ext(nx+1) = -qy_wall / deltaT

  !-----------------------------
  ! (1.2) 网格上找 Numax/Numin（含角点 0 与 nx+1）
  imax = 0
  imin = 0
  Nu_hot_max = Nu_bot_ext(0)
  Nu_hot_min = Nu_bot_ext(0)

  do i = 1, nx+1
    if (Nu_bot_ext(i) > Nu_hot_max) then
      Nu_hot_max = Nu_bot_ext(i)
      imax = i
    endif
    if (Nu_bot_ext(i) < Nu_hot_min) then
      Nu_hot_min = Nu_bot_ext(i)
      imin = i
    endif
  enddo

  !-----------------------------
  ! (1.3) 对 Numax：取 5 点（尽量对称；靠近端点时单侧）并做 LS 抛物线拟合
  if (imax <= 2) then
    ii = (/ 0, 1, 2, 3, 4 /)
  elseif (imax >= nx-1) then
    ii = (/ nx-3, nx-2, nx-1, nx, nx+1 /)
  else
    ii = (/ imax-2, imax-1, imax, imax+1, imax+2 /)
  endif

  do k = 1, 5
    xk(k) = xp(ii(k))            ! xp(0)=0, xp(nx+1)=nx/lengthUnit
    fk(k) = Nu_bot_ext(ii(k))
  enddo

  call fit_parabola_ls5(xk, fk, +1, fstar, xstar)
  Nu_hot_max = fstar
  Nu_hot_max_position = xstar     ! RB: 这里是 x 位置

  !-----------------------------
  ! (1.4) 对 Numin：同样取 5 点并拟合
  if (imin <= 2) then
    ii = (/ 0, 1, 2, 3, 4 /)
  elseif (imin >= nx-1) then
    ii = (/ nx-3, nx-2, nx-1, nx, nx+1 /)
  else
    ii = (/ imin-2, imin-1, imin, imin+1, imin+2 /)
  endif

  do k = 1, 5
    xk(k) = xp(ii(k))
    fk(k) = Nu_bot_ext(ii(k))
  enddo

  call fit_parabola_ls5(xk, fk, -1, fstar, xstar)
  Nu_hot_min = fstar
  Nu_hot_min_position = xstar     ! RB: 这里是 x 位置

  !-----------------------------
  ! (2) 顶部冷壁平均 Nu_cold（不含角点）
  sum_cold = 0.0d0
  !$omp parallel do default(none) shared(T_all,quadWidthX,deltaT,yp) &
  !$omp private(i,qy_wall,ynode,fnode) reduction(+:sum_cold)
  do i = 1, nx
    ynode = (/ yp(ny-1), yp(ny), 1.0d0 /)
    fnode = (/ T_all(i,ny-1), T_all(i,ny), Tcold /)
    call lagrange_derivative_3(ynode, fnode, 1.0d0, qy_wall)
    qy_wall = -qy_wall
    sum_cold = sum_cold + qy_wall/deltaT*quadWidthX(i)
  enddo
  !$omp end parallel do
  Nu_cold = sum_cold / quadSumX

  !-----------------------------
  ! 中线的 Nusselt 数
  sum_mid = 0.0d0

  call build_lagrange_stencil_1d(ny, yp(1:ny), 0.5d0, midIndex, ymidWeight, midStencilValid)
  if(.not.midStencilValid) then
    write(*,*) "Error: y=0.5 is outside ISLBM y nodes in RBcalc_Nu_wall_avg"
    stop
  endif
  ymidNode = (/ yp(midIndex(1)), yp(midIndex(2)), yp(midIndex(3)) /)

  !$omp parallel do default(none) shared(v_all,T_all,quadWidthX,deltaT,coef,midIndex,ymidNode,ymidWeight) &
  !$omp private(i,v_mid,T_mid,dTdy_mid,fnode) reduction(+:sum_mid)
  do i = 1, nx
    v_mid = ymidWeight(1)*v_all(i,midIndex(1)) + ymidWeight(2)*v_all(i,midIndex(2)) + &
      ymidWeight(3)*v_all(i,midIndex(3))
    T_mid = ymidWeight(1)*T_all(i,midIndex(1)) + ymidWeight(2)*T_all(i,midIndex(2)) + &
      ymidWeight(3)*T_all(i,midIndex(3))
    fnode = (/ T_all(i,midIndex(1)), T_all(i,midIndex(2)), T_all(i,midIndex(3)) /)
    call lagrange_derivative_3(ymidNode, fnode, 0.5d0, dTdy_mid)
    sum_mid = sum_mid + (coef*v_mid*(T_mid-Tref) - dTdy_mid)/deltaT*quadWidthX(i)
  enddo
  !$omp end parallel do

  Nu_middle = sum_mid / quadSumX

  !-----------------------------
  ! 输出：屏幕 + 日志
  write(*,'(a,1x,es16.8)') "Nu_hot(bottom) =", Nu_hot
  write(*,'(a,1x,es16.8)') "Nu_cold(top)   =", Nu_cold
  write(*,'(a,1x,es16.8)') "Nu_middle      =", Nu_middle
  write(*,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_max =", Nu_hot_max, "x_max =", Nu_hot_max_position
  write(*,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_min =", Nu_hot_min, "x_min =", Nu_hot_min_position

  open(unit=00,file=trim(settingsFile),status="unknown",position="append")
  write(00,'(a,1x,es16.8)') "Nu_hot(bottom) =", Nu_hot
  write(00,'(a,1x,es16.8)') "Nu_cold(top)   =", Nu_cold
  write(00,'(a,1x,es16.8)') "Nu_middle      =", Nu_middle
  write(00,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_max =", Nu_hot_max, "x_max =", Nu_hot_max_position
  write(00,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_min =", Nu_hot_min, "x_min =", Nu_hot_min_position
  close(00)

  return
end subroutine RBcalc_Nu_wall_avg
!===================================================================================================
! RBcalc_Nu_wall_avg 结束: 计算 Rayleigh-Benard 工况下热壁、冷壁和中线的 Nusselt 数及其极值。
!===================================================================================================


!===================================================================================================
! 子程序: RBcalc_umid_max
! 作用: 计算 Rayleigh-Benard 工况下 x=1/2 处最大水平速度 umax 及其 y 位置。
! 用途: 在 RayleighBenardCell 工况结束后的后处理中调用。
!===================================================================================================
subroutine RBcalc_umid_max()
  use commondata
  implicit none
  integer(kind=4) :: j, k
  integer(kind=4) :: midIndex(3)
  integer(kind=4) :: j0
  integer(kind=4) :: jj(5)
  real(kind=8) :: uline(1:ny)
  real(kind=8) :: yk(5), fk(5)
  real(kind=8) :: umax_fit, y_fit
  real(kind=8) :: coef, xmid, targetX, xmidWeight(3)
  real(kind=8) :: umax_grid, yLen
  logical :: midStencilValid

  coef = velocityScaleCompare

  ! ---- (1) 取论文定义的 x=1/2 竖线剖面 u(x=0.5, y_j) ----
  targetX = 0.5d0
  xmid = targetX
  call build_lagrange_stencil_1d(nx, xp(1:nx), targetX, midIndex, xmidWeight, midStencilValid)
  if(.not.midStencilValid) then
    write(*,*) "Error: targetX is outside ISLBM x nodes in RBcalc_umid_max"
    stop
  endif
  do j = 1, ny
    uline(j) = xmidWeight(1)*u_all(midIndex(1),j) + xmidWeight(2)*u_all(midIndex(2),j) + &
      xmidWeight(3)*u_all(midIndex(3),j)
  enddo

  ! ---- (2) 找峰值所在网格点 j0 ----
  j0 = 1
  umax_grid = uline(1)
  do j = 2, ny
    if (uline(j) > umax_grid) then
      umax_grid = uline(j)
      j0 = j
    endif
  enddo

  ! ---- (3) 取 5 点（尽量对称；靠近端点时偏侧）----
  if (j0 <= 2) then
    jj = (/ 1, 2, 3, 4, 5 /)
  elseif (j0 >= ny-1) then
    jj = (/ ny-4, ny-3, ny-2, ny-1, ny /)
  else
    jj = (/ j0-2, j0-1, j0, j0+1, j0+2 /)
  endif

  do k = 1, 5
    yk(k) = yp(jj(k))
    fk(k) = uline(jj(k))
  enddo

  call fit_parabola_ls5(yk, fk, +1, umax_fit, y_fit)

  ! 对称支路下，上下半区的 y 位置等价；统一折回下半区便于与 benchmark 对照。
  yLen = yp(ny+1)
  if (y_fit > 0.5d0*yLen) y_fit = yLen - y_fit
  

  write(*,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'u_mid_max* =', umax_fit*coef, 'y =', y_fit, 'x_mid =', xmid

  open(unit=00,file=trim(settingsFile),status="unknown",position="append")
  write(00,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'u_mid_max* =', umax_fit*coef, 'y =', y_fit, 'x_mid =', xmid
  close(00)

  return
end subroutine RBcalc_umid_max
!===================================================================================================
! RBcalc_umid_max 结束: 计算 Rayleigh-Benard 工况下 x=1/2 处最大水平速度及其位置。
!===================================================================================================


!===================================================================================================
! 子程序: RBcalc_vmid_max
! 作用: 计算 Rayleigh-Benard 工况下 y=1/2 处最大垂直速度 vmax 及其 x 位置。
! 用途: 在 RayleighBenardCell 工况结束后的后处理中调用。
!===================================================================================================
subroutine RBcalc_vmid_max()
  use commondata
  implicit none
  integer(kind=4) :: i, k
  integer(kind=4) :: midIndex(3)
  integer(kind=4) :: i0
  integer(kind=4) :: ii(5)
  real(kind=8) :: vline(1:nx)
  real(kind=8) :: xk(5), fk(5)
  real(kind=8) :: vmax_fit, x_fit
  real(kind=8) :: coef, ymid, ymidWeight(3)
  real(kind=8) :: vmax_grid
  logical :: midStencilValid

  coef = velocityScaleCompare

  ! ---- (1) 取 y=1/2 中线剖面 v(x_i, y=1/2) ----
  ymid = 0.5d0
  call build_lagrange_stencil_1d(ny, yp(1:ny), ymid, midIndex, ymidWeight, midStencilValid)
  if(.not.midStencilValid) then
    write(*,*) "Error: ymid is outside ISLBM y nodes in RBcalc_vmid_max"
    stop
  endif
  do i = 1, nx
    vline(i) = ymidWeight(1)*v_all(i,midIndex(1)) + ymidWeight(2)*v_all(i,midIndex(2)) + &
      ymidWeight(3)*v_all(i,midIndex(3))
  enddo

  ! ---- (2) 找峰值所在网格点 i0 ----
  i0 = 1
  vmax_grid = vline(1)
  do i = 2, nx
    if (vline(i) > vmax_grid) then
      vmax_grid = vline(i)
      i0 = i
    endif
  enddo


  ! ---- (3) 取 5 点并拟合 ----
  if (i0 <= 2) then
    ii = (/ 1, 2, 3, 4, 5 /)
  elseif (i0 >= nx-1) then
    ii = (/ nx-4, nx-3, nx-2, nx-1, nx /)
  else
    ii = (/ i0-2, i0-1, i0, i0+1, i0+2 /)
  endif

  do k = 1, 5
    xk(k) = xp(ii(k))
    fk(k) = vline(ii(k))
  enddo

  call fit_parabola_ls5(xk, fk, +1, vmax_fit, x_fit)

  write(*,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'v_mid_max* =', vmax_fit*coef, 'x =', x_fit, 'y_mid =', ymid

  open(unit=00,file=trim(settingsFile),status="unknown",position="append")
  write(00,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'v_mid_max* =', vmax_fit*coef, 'x =', x_fit, 'y_mid =', ymid
  close(00)

  return
end subroutine RBcalc_vmid_max
!===================================================================================================
! RBcalc_vmid_max 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
!===================================================================================================

!===================================================================================================






!===================================================================================================
! 子程序: calc_psi_vort_and_output
! 作用: 计算流函数、涡量，并输出相关诊断量。
! 用途: 在主程序结束阶段调用，作为统一后处理的一部分。
!===================================================================================================
subroutine calc_psi_vort_and_output()
  use commondata
  implicit none

  integer(kind=4) :: i, j
  real(kind=8) :: coef
  real(kind=8) :: segmentIntegral
  real(kind=8) :: xNode(3), vNode(3)
  real(kind=8) :: yNode(3), uNode(3)
  real(kind=8) :: dv_dx, du_dy
  real(kind=8) :: psi(nx,ny), vort(nx,ny)
  real(kind=8) :: psiTopAbsMax


  ! for fine-grid max(|psi|)  10001*10001
  real(kind=8) :: psi_abs_max, x_at_max, y_at_max, psi_center_abs_fine

  coef = velocityScaleCompare

  !=========================================================
  ! (A) 计算流函数 psi。
  !     SideHeatedCell 按zzhao后处理从左壁积分: psi(x,y)=-∫_0^x v(mu,y)dmu。
  !     RayleighBenardCell 保留从底壁积分: psi(x,y)=∫_0^y u(x,mu)dmu。
  !     ISLBM非均匀网格下不能使用固定dx/dy的Simpson公式。
  !     每一小段用三点二次Lagrange多项式做解析积分。
  !=========================================================

#ifdef SideHeatedCell
  do j = 1, ny
    xNode = (/ 0.0d0, xp(1), xp(2) /)
    vNode = (/ 0.0d0, v_all(1,j), v_all(2,j) /)
    call integrate_lagrange_3_segment(xNode, vNode, 0.0d0, xp(1), segmentIntegral)
    psi(1,j) = -segmentIntegral*coef

    do i = 2, nx
      if (i == 2) then
        xNode = (/ 0.0d0, xp(1), xp(2) /)
        vNode = (/ 0.0d0, v_all(1,j), v_all(2,j) /)
      else
        xNode = (/ xp(i-2), xp(i-1), xp(i) /)
        vNode = (/ v_all(i-2,j), v_all(i-1,j), v_all(i,j) /)
      endif
      call integrate_lagrange_3_segment(xNode, vNode, xp(i-1), xp(i), segmentIntegral)
      psi(i,j) = psi(i-1,j) - segmentIntegral*coef
    end do

  end do
#endif

#ifdef RayleighBenardCell
  do i = 1, nx
    yNode = (/ 0.0d0, yp(1), yp(2) /)
    uNode = (/ 0.0d0, u_all(i,1), u_all(i,2) /)
    call integrate_lagrange_3_segment(yNode, uNode, 0.0d0, yp(1), segmentIntegral)
    psi(i,1) = segmentIntegral*coef

    do j = 2, ny
      if (j == 2) then
        yNode = (/ 0.0d0, yp(1), yp(2) /)
        uNode = (/ 0.0d0, u_all(i,1), u_all(i,2) /)
      else
        yNode = (/ yp(j-2), yp(j-1), yp(j) /)
        uNode = (/ u_all(i,j-2), u_all(i,j-1), u_all(i,j) /)
      endif
      call integrate_lagrange_3_segment(yNode, uNode, yp(j-1), yp(j), segmentIntegral)
      psi(i,j) = psi(i,j-1) + segmentIntegral*coef
    end do

  end do
#endif

  psiTopAbsMax = 0.0d0
#ifdef SideHeatedCell
  do j = 1, ny
    psiTopAbsMax = max(psiTopAbsMax, dabs(psi(nx,j)))
  enddo
  write(*,'(a,1x,es16.8)') "max(|psi_right_internal|) =", psiTopAbsMax
  open(unit=00,file=trim(settingsFile),status="unknown",position="append")
  write(00,'(a,1x,es16.8)') "max(|psi_right_internal|) =", psiTopAbsMax
  close(00)
#endif
#ifdef RayleighBenardCell
  do i = 1, nx
    psiTopAbsMax = max(psiTopAbsMax, dabs(psi(i,ny)))
  enddo
  write(*,'(a,1x,es16.8)') "max(|psi_top_internal|) =", psiTopAbsMax
  open(unit=00,file=trim(settingsFile),status="unknown",position="append")
  write(00,'(a,1x,es16.8)') "max(|psi_top_internal|) =", psiTopAbsMax
  close(00)
#endif

  call output_psi_center_abs(psi)     ! 基于粗网格局部四点插值得到中心点处的 abs(psi)

  !=========================================================
  ! (B) 计算涡量 vort = dv/dx - du/dy（2D）
  !     ISLBM非均匀网格下用三点Lagrange导数, 不再使用固定dx/dy中心差分。
  !=========================================================

  do j = 1, ny
    do i = 1, nx

      if(i.EQ.1) then
        xNode = (/ xp(1), xp(2), xp(3) /)
        vNode = (/ v_all(1,j), v_all(2,j), v_all(3,j) /)
      elseif(i.EQ.nx) then
        xNode = (/ xp(nx-2), xp(nx-1), xp(nx) /)
        vNode = (/ v_all(nx-2,j), v_all(nx-1,j), v_all(nx,j) /)
      else
        xNode = (/ xp(i-1), xp(i), xp(i+1) /)
        vNode = (/ v_all(i-1,j), v_all(i,j), v_all(i+1,j) /)
      endif
      call lagrange_derivative_3(xNode, vNode, xp(i), dv_dx)

      if(j.EQ.1) then
        yNode = (/ yp(1), yp(2), yp(3) /)
        uNode = (/ u_all(i,1), u_all(i,2), u_all(i,3) /)
      elseif(j.EQ.ny) then
        yNode = (/ yp(ny-2), yp(ny-1), yp(ny) /)
        uNode = (/ u_all(i,ny-2), u_all(i,ny-1), u_all(i,ny) /)
      else
        yNode = (/ yp(j-1), yp(j), yp(j+1) /)
        uNode = (/ u_all(i,j-1), u_all(i,j), u_all(i,j+1) /)
      endif
      call lagrange_derivative_3(yNode, uNode, yp(j), du_dy)

      ! psi 用了L/kappa；vort 也应在同一标度下
      ! vort = d(v*coef)/dx - d(u*coef)/dy = coef*(dv/dx - du/dy)
      vort(i,j) = coef * (dv_dx - du_dy)

    end do
  end do

  !=========================================================
  ! (C) 输出 Tecplot：X,Y,psi,vort
  !=========================================================
  call output_Tecplot_psi_vort(psi, vort)

  !=========================================================
  ! (D) 细网格 10001×10001：用三次样条插值寻找 max(|psi|)
  !=========================================================
  call calc_psi_absmax_fine_spline(psi, psi_abs_max, x_at_max, y_at_max, psi_center_abs_fine)

  write(*,'(a,1x,es16.8)') "abs(psi_center_fine) =", psi_center_abs_fine

  write(*,'(a,1x,es16.8,2x,a,1x,es16.8,2x,a,1x,es16.8)') &
       "max(|psi|) =", psi_abs_max, "x* =", x_at_max, "y* =", y_at_max

  open(unit=00,file=trim(settingsFile),status="unknown",position="append")
  write(00,'(a,1x,es16.8)') "abs(psi_center_fine) =", psi_center_abs_fine
  write(00,'(a,1x,es16.8,2x,a,1x,es16.8,2x,a,1x,es16.8)') &
       "max(|psi|) =", psi_abs_max, "x* =", x_at_max, "y* =", y_at_max
  close(00)

  return
end subroutine calc_psi_vort_and_output
!===================================================================================================
! calc_psi_vort_and_output 结束: 计算流函数、涡量，并输出相关诊断量。
!===================================================================================================



!===================================================================================================
! 子程序: output_Tecplot_psi_vort
! 作用: 把流函数和涡量场写出为 Tecplot 文件。
! 用途: 在 calc_psi_vort_and_output 中调用。
!===================================================================================================
subroutine output_Tecplot_psi_vort(psi, vort)
  use commondata
  implicit none
  real(kind=8), intent(in) :: psi(nx,ny), vort(nx,ny)

  integer(kind=4) :: i, j, k
  real(kind=4) :: zoneMarker, eohMarker
  character(len=40) :: title
  character(len=40) :: V1,V2,V3,V4
  integer(kind=4), parameter :: kmax=1
  character(len=40) :: zoneName
  character(len=100) :: filename

#ifdef steadyFlow
  write(filename,'(i12.12)') restartItcOffset+itc
#endif
#ifdef unsteadyFlow
  pltFileNum = pltFileNum+1
  write(filename,'(i12.12)') pltFileNum
#endif
  filename = adjustl(filename)

  open(41,file=trim(pltFolderPrefix)//"-psiVort-"//trim(filename)//'.plt', access='stream', form='unformatted')

  zoneMarker= 299.0
  eohMarker = 357.0

  write(41) "#!TDV101"
  write(41) 1

  title="psi-vort"
  call dumpstring(title)

  write(41) 4
  V1='X';     call dumpstring(V1)
  V2='Y';     call dumpstring(V2)
  V3='PSI';   call dumpstring(V3)
  V4='VORT';  call dumpstring(V4)

  write(41) zoneMarker
  zoneName="ZONE 001"
  call dumpstring(zoneName)

  write(41) -1
  write(41) 0
  write(41) 1
  write(41) 0
  write(41) 0

  write(41) nx
  write(41) ny
  write(41) kmax

  write(41) 0
  write(41) eohMarker

  write(41) zoneMarker

  ! 2 = Double
  write(41) 2
  write(41) 2
  write(41) 2
  write(41) 2

  write(41) 0
  write(41) -1

  do k=1,kmax
    do j=1,ny
      do i=1,nx
        write(41) real(xp(i),kind=8)
        write(41) real(yp(j),kind=8)
        write(41) real(psi(i,j),kind=8)
        write(41) real(vort(i,j),kind=8)
      end do
    end do
  end do

  close(41)
  return
end subroutine output_Tecplot_psi_vort
!===================================================================================================
! output_Tecplot_psi_vort 结束: 把流函数和涡量场写出为 Tecplot 文件。
!===================================================================================================


!===================================================================================================
! 子程序: calc_psi_absmax_fine_spline
! 作用: 用细网格样条插值搜索 abs(psi) 的最大值及位置。
! 用途: 在 calc_psi_vort_and_output 中调用。
!===================================================================================================
subroutine calc_psi_absmax_fine_spline(psi, psi_abs_max, x_at_max, y_at_max, psi_center_abs_fine)
  use commondata
  implicit none

  real(kind=8), intent(in)  :: psi(nx,ny)
  real(kind=8), intent(out) :: psi_abs_max, x_at_max, y_at_max, psi_center_abs_fine

  integer(kind=4), parameter :: nFinePerUnit = 10000
  integer(kind=4) :: i, j, k, l
  integer(kind=4) :: nFineX, nFineY
  integer(kind=4) :: nXExt, nYExt

  real(kind=8), allocatable :: xFine(:), yFine(:)
  real(kind=8), allocatable :: xExt(:), yExt(:)
  real(kind=8), allocatable :: row(:), y2row(:)
  real(kind=8), allocatable :: col(:), y2col(:)
  real(kind=8), allocatable :: psi_center_x(:)
  real(kind=8), allocatable :: psi_xfine(:,:)

  real(kind=8) :: val, xq, yq
  real(kind=8) :: xLen, yLen, xCenter, yCenter

  ! ---- 细网格坐标：0..1 等距
  ! Use the actual nondimensional domain size so the fine mesh also works for non-square cavities.
  xLen = xp(nx+1)
  yLen = yp(ny+1)
  xCenter = 0.5d0 * xLen
  yCenter = 0.5d0 * yLen
  nFineX = max(2, nint(xLen * dble(nFinePerUnit)) + 1)
  nFineY = max(2, nint(yLen * dble(nFinePerUnit)) + 1)

  allocate(xFine(nFineX), yFine(nFineY))
  do k = 1, nFineX
    xFine(k) = xLen * dble(k-1) / dble(nFineX-1)
  end do
  do l = 1, nFineY
    yFine(l) = yLen * dble(l-1) / dble(nFineY-1)
  end do

  ! ---- 为了能在 x=0/1 与 y=0/1 上插值，物理边界 psi=常数（可取0）补两个端点
  nXExt = nx + 2
  nYExt = ny + 2

  allocate(xExt(nXExt), yExt(nYExt))
  xExt(1)    = 0.0d0
  xExt(nXExt)= xp(nx+1)
  do i = 1, nx
    xExt(i+1) = xp(i)         !xExt = [0, xp(1),...,xp(nx), xLen]
  end do

  yExt(1)    = 0.0d0
  yExt(nYExt)= yp(ny+1)
  do j = 1, ny
    yExt(j+1) = yp(j)        !yExt = [0, yp(1),...,yp(ny), yLen]
  end do

  allocate(row(nXExt), y2row(nXExt))
  allocate(psi_center_x(ny))
  allocate(psi_xfine(nFineX, ny))

  ! ---- (1) 先对每个固定 y=yp(j) 的剖面做 x 方向三次样条，得到 psi(xFine, yp(j))
  do j = 1, ny
    row(1)  = 0.0d0
    row(nXExt) = 0.0d0
    do i = 1, nx
      row(i+1) = psi(i,j)     !row=[0, psi(1),...,psi(ny), 0]
    end do

    call spline_natural(nXExt, xExt, row, y2row)   !在row各节点处的二阶导

    call splint(nXExt, xExt, row, y2row, xCenter, psi_center_x(j))
    do k = 1, nFineX
      xq = xFine(k)
      call splint(nXExt, xExt, row, y2row, xq, val)   !在 10001 个细网格 x 点上采样
      psi_xfine(k,j) = val                            !在每个粗 y 层上，psi已经沿 x 被细化到 10001 个点。
    end do
  end do

  ! ---- (2) 再对每个固定 x=xFine(k) 的剖面做 y 方向三次样条，在 10001 个 yFine 上扫 max(|psi|)
  allocate(col(nYExt), y2col(nYExt))
  psi_abs_max = -1.0d0
  x_at_max    = 0.0d0
  y_at_max    = 0.0d0

  col(1)  = 0.0d0
  col(nYExt) = 0.0d0
  do j = 1, ny
    col(j+1) = psi_center_x(j)
  end do
  call spline_natural(nYExt, yExt, col, y2col)
  call splint(nYExt, yExt, col, y2col, yCenter, val)
  psi_center_abs_fine = dabs(val)

  do k = 1, nFineX
    col(1)  = 0.0d0
    col(nYExt) = 0.0d0
    do j = 1, ny
      col(j+1) = psi_xfine(k,j)        !对每一个固定的细网格 x=xFine(k)，有一条沿 y 的离散剖面数据 col
    end do

    call spline_natural(nYExt, yExt, col, y2col)

    do l = 1, nFineY
      yq = yFine(l)
      call splint(nYExt, yExt, col, y2col, yq, val)

      if (dabs(val) > psi_abs_max) then     !寻找最大的abs(psi)以及位置
        psi_abs_max = dabs(val)
        x_at_max = xFine(k)
        y_at_max = yFine(l)
      end if
    end do
  end do

  deallocate(xFine, yFine, xExt, yExt, row, y2row, col, y2col, psi_center_x, psi_xfine)
  return
end subroutine calc_psi_absmax_fine_spline
!===================================================================================================
! calc_psi_absmax_fine_spline 结束: 用细网格样条插值搜索 abs(psi) 的最大值及位置。
!===================================================================================================


!===================================================================================================
! 子程序: spline_natural
! 作用: 构造自然三次样条所需的二阶导数。
! 用途: 在 calc_psi_absmax_fine_spline 中作为样条预处理调用。
!===================================================================================================
subroutine spline_natural(n, x, y, y2)
  implicit none
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in)    :: x(*), y(*)
  real(kind=8), intent(out)   :: y2(*)

  integer(kind=4) :: i, k
  real(kind=8), allocatable :: u(:)
  real(kind=8) :: sig, p

  allocate(u(n))

  ! natural spline: y2(1)=0, y2(n)=0
  y2(1) = 0.0d0
  u(1)  = 0.0d0

  do i = 2, n-1
    sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
    p = sig*y2(i-1) + 2.0d0
    y2(i) = (sig - 1.0d0) / p
    u(i) = ( 6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))/(x(i)-x(i-1))) / (x(i+1)-x(i-1)) - sig*u(i-1) ) / p
  end do

  y2(n) = 0.0d0

  do k = n-1, 1, -1
    y2(k) = y2(k)*y2(k+1) + u(k)
  end do

  deallocate(u)
  return
end subroutine spline_natural
!===================================================================================================
! spline_natural 结束: 构造自然三次样条所需的二阶导数。
!===================================================================================================


!===================================================================================================
! 子程序: splint
! 作用: 根据样条系数在查询点上进行插值。
! 用途: 在 calc_psi_absmax_fine_spline 中重复调用。
!===================================================================================================
subroutine splint(n, xa, ya, y2a, x, y)
  implicit none
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in)    :: xa(*), ya(*), y2a(*), x
  real(kind=8), intent(out)   :: y

  integer(kind=4) :: klo, khi, k
  real(kind=8) :: h, a, b

  klo = 1
  khi = n

  do while (khi - klo > 1)
    k = (khi + klo)/2
    if (xa(k) > x) then
      khi = k
    else
      klo = k
    end if
  end do

  h = xa(khi) - xa(klo)
  if (h == 0.0d0) then
    y = ya(klo)
    return
  end if

  a = (xa(khi) - x)/h
  b = (x - xa(klo))/h

  y = a*ya(klo) + b*ya(khi) + ( (a*a*a - a)*y2a(klo) + (b*b*b - b)*y2a(khi) ) * (h*h)/6.0d0
  return
end subroutine splint
!===================================================================================================
! splint 结束: 根据样条系数在查询点上进行插值。
!===================================================================================================


!===================================================================================================
! 子程序: output_psi_center_abs
! 作用: 输出腔体中心位置的 abs(psi) 诊断结果。
! 用途: 在 calc_psi_vort_and_output 中调用。
!===================================================================================================
subroutine output_psi_center_abs(psi)
  use commondata
  implicit none
  real(kind=8), intent(in) :: psi(nx,ny)

  integer(kind=4) :: i0, j0, p, q
  integer(kind=4) :: ii(4), jj(4)
  real(kind=8) :: x0, y0, psi_center, psi_center_abs
  real(kind=8) :: x4(4), y4(4), f4(4), gx(4)

  x0 = 0.5d0 * xp(nx+1)
  y0 = 0.5d0 * yp(ny+1)

  if (nx < 4 .or. ny < 4) then
    psi_center = psi((nx+1)/2, (ny+1)/2)
  else
    i0 = 1
    do while (i0 < nx .and. xp(i0+1) <= x0)
      i0 = i0 + 1
    end do
    i0 = max(1, min(i0-1, nx-3))
    do p = 1, 4
      ii(p) = i0 + p - 1
      x4(p) = xp(ii(p))
    end do

    j0 = 1
    do while (j0 < ny .and. yp(j0+1) <= y0)
      j0 = j0 + 1
    end do
    j0 = max(1, min(j0-1, ny-3))
    do q = 1, 4
      jj(q) = j0 + q - 1
      y4(q) = yp(jj(q))
    end do

    do q = 1, 4
      do p = 1, 4
        f4(p) = psi(ii(p), jj(q))
      end do
      call interp_lagrange_4(x0, x4, f4, gx(q))
    end do

    call interp_lagrange_4(y0, y4, gx, psi_center)
  endif

  psi_center_abs = dabs(psi_center)

  ! Screen output
  write(*,'(a,1x,es16.8)') "abs(psi_center_coarse) =", psi_center_abs

  ! Log output
  open(unit=00,file=trim(settingsFile),status="unknown",position="append")
  write(00,'(a,1x,es16.8)') "abs(psi_center_coarse) =", psi_center_abs
  close(00)

  return
contains
  ! 子程序: interp_lagrange_4
  ! 作用: 对给定的四个节点执行四点 Lagrange 插值。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
  subroutine interp_lagrange_4(xq, xk, fk, fq)
    implicit none
    real(kind=8), intent(in)  :: xq
    real(kind=8), intent(in)  :: xk(4), fk(4)
    real(kind=8), intent(out) :: fq
    integer(kind=4) :: a, b
    real(kind=8) :: basis

    fq = 0.0d0
    do a = 1, 4
      basis = 1.0d0
      do b = 1, 4
        if (b /= a) basis = basis * (xq - xk(b)) / (xk(a) - xk(b))
      end do
      fq = fq + fk(a) * basis
    end do
  end subroutine interp_lagrange_4
end subroutine output_psi_center_abs
!===================================================================================================
! output_psi_center_abs 结束: 输出腔体中心位置的 abs(psi) 诊断结果。
!===================================================================================================

