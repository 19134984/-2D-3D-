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

! 温度算法使用原 D2Q5 TRT 旧算法，平衡矩系数为 paraA。

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
        integer(kind=4), parameter :: nx = 1024, ny = 1024    ! 手动设置最细网格等效分辨率，需要被refineRatio整除得到粗网格间距
        integer(kind=4), parameter :: refineRatio = 2    ! 粗/细格距及时间步之比；逐级二分取2、4、8等2的幂，1为单块
        integer(kind=4), parameter :: fineLayerCellsLeft = 128      ! 从左墙向内数的细节点编号，x = Left-0.5
        integer(kind=4), parameter :: fineLayerCellsRight = 129     ! 从右墙向内数的细节点编号，x = nx-Right+0.5
        integer(kind=4), parameter :: fineLayerCellsBottom = 128    ! 从下墙向内数的细节点编号，y = Bottom-0.5
        integer(kind=4), parameter :: fineLayerCellsTop = 129       ! 从上墙向内数的细节点编号，y = ny-Top+0.5
        ! 默认交界面为 x/y=127.5、895.5；粗细块共用这些基准边界，积分也在此分区。
        ! 中心宽度 nx-Left-Right+1、高度 ny-Bottom-Top+1 须为正，且均能被 refineRatio 整除。
        ! 右=左+1、上=下+1 时，中心宽高为 nx-2*Left、ny-2*Bottom；较大粗细比仍须检查整除。
        ! 默认 128/129 对 nx=ny=1024、粗细比 2/4/8 都满足条件；不再自动移动右/上交界面。
        ! 多块模式 nx、ny 仍须整除粗细比；四个细网格编号fineLayerCellsLeft本身无须整除，但距墙均至少 (overlapCells+1)*refineRatio，得保证重叠层足够。
        integer(kind=4), parameter :: overlapCells = 2    ! 每块向交界面外延伸的粗格距数；总重叠跨度为两倍粗格距数
        integer(kind=4), parameter :: interfaceSkin = 2    ! 粗细人工边缘上每次重建的本地节点层数；细区之间直接迁移，不使用此层数。

        !===============================================================================================
        ! 是否从本程序检查点精确续算
        !===============================================================================================
        integer(kind=4), parameter :: loadInitField = 0    ! 0: 重新初始化；1: 按 latest.meta 指向的检查点续算
        ! 正常续算只设置 loadInitField=1；优先读取 <reloadFilePrefix>-latest.meta。
        ! 只有 latest.meta 缺失时，才手动填写 reloadFileNum，选择同前缀的编号 .bin 文件。
        ! 时间、输出计数及粗块历史都从 .bin 精确恢复，不按文件编号推断累计时间。
        integer(kind=4) :: reloadFileNum = 0    ! 非稳态为独立输出编号；稳态为累计细格子步数

        !===============================================================================================
        ! 无量纲参数及物理壁面温度
        !===============================================================================================
        real(kind=8), parameter :: Rayleigh = 1.0d7    ! 手动设置瑞利数，修改后重新编译
        real(kind=8), parameter :: Prandtl = 0.7d0
        real(kind=8), parameter :: Mach = 0.1d0
        real(kind=8), parameter :: Thot = 0.5d0
        real(kind=8), parameter :: Tcold = -0.5d0
        real(kind=8), parameter :: Tref = 0.5d0*(Thot+Tcold)
        real(kind=8), parameter :: pi = acos(-1.0d0)
#ifdef SideHeatedCell
        real(kind=8), parameter :: lengthUnit = dble(nx)
#else
        real(kind=8), parameter :: lengthUnit = dble(ny)
#endif

        !===============================================================================================
        ! 输运系数：沿用最细格子单位
        !===============================================================================================
        ! 以下参数全部沿用原文件的最细格子单位，不按块重新定义 Ra、Pr、Ma 或 paraA。
        real(kind=8), parameter :: tauf = 0.5d0+Mach*lengthUnit*sqrt(3.0d0*Prandtl/Rayleigh)
        real(kind=8), parameter :: viscosity = (tauf-0.5d0)/3.0d0
        real(kind=8), parameter :: diffusivity = viscosity/Prandtl

        ! 浮力项及无量纲时间、速度
        real(kind=8), parameter :: gBeta1 = Rayleigh*viscosity*diffusivity/lengthUnit
        real(kind=8), parameter :: gBeta = gBeta1/lengthUnit/lengthUnit
        real(kind=8), parameter :: timeUnit = sqrt(lengthUnit/gBeta)
        real(kind=8), parameter :: velocityUnit = sqrt(gBeta*lengthUnit)
        real(kind=8), parameter :: velocityScaleCompare = lengthUnit/diffusivity

        ! 流场和温度场的多松弛系数
        real(kind=8), parameter :: Snu = 1.0d0/tauf
        real(kind=8), parameter :: Sq = 8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)
        real(kind=8), parameter :: paraA = 20.0d0*sqrt(3.0d0)*diffusivity-4.0d0
        real(kind=8), parameter :: Qk = 3.0d0-sqrt(3.0d0)
        real(kind=8), parameter :: Qnu = 4.0d0*sqrt(3.0d0)-6.0d0
        real(kind=8), parameter :: thermalGeqCoeff = 10.0d0/(4.0d0+paraA)
        real(kind=8), parameter :: thermalA = paraA
#ifdef SideHeatedHa
        real(kind=8), parameter :: Ha = 20.0d0
        real(kind=8), parameter :: phi = 0.0d0*pi/180.0d0
        real(kind=8), parameter :: B2sigemarho = Ha**2*viscosity/lengthUnit**2
#endif

        !===============================================================================================
        ! 收敛阈值、计算时间及输出间隔
        !===============================================================================================
        real(kind=8), parameter :: epsU = 1.0d-7
        real(kind=8), parameter :: epsT = 1.0d-7
#ifdef steadyFlow
        real(kind=8), parameter :: outputSnapshotInterval = 10.0d0    ! 快照与 Nu/Re 采样间隔，单位 t_ff
        real(kind=8), parameter :: reloadFileInterval = 100.0d0      ! 完整重启文件输出间隔，单位 t_ff
        real(kind=8), parameter :: outputPltFileInterval = 100.0d0   ! Tecplot 输出间隔，单位 t_ff
        integer(kind=4), parameter :: outputSnapshotFile = 1    ! 0: 不输出快照；1: 输出，Nu/Re 仍独立采样
        integer(kind=4), parameter :: outputPltFile = 1         ! 0: 不输出 Tecplot；1: 输出
        integer(kind=4), parameter :: outputReloadFile = 1      ! 0: 不输出重启文件；1: 输出
        integer(kind=4), parameter :: itc_max = 20000000        ! 最大细格子步数，收敛后可提前停止
#endif

#ifdef unsteadyFlow
        real(kind=8), parameter :: outputSnapshotInterval = 0.5d0    ! 快照与 Nu/Re 采样间隔，单位 t_ff
        real(kind=8), parameter :: reloadFileInterval = 100.0d0      ! 完整重启文件输出间隔，单位 t_ff
        real(kind=8), parameter :: outputPltFileInterval = 100.0d0   ! Tecplot 输出间隔，单位 t_ff
        real(kind=8), parameter :: unsteadyRunDuration = 1000.0d0    ! 绝对总目标 t_ff；续算只补足剩余时间
        ! 以下三个参数只控制结束后的统计窗口，不改变推进时长或采样频率；包含匹配的旧历史数据。
        real(kind=8), parameter :: unsteadyAverageStartTf = 0.5d0*unsteadyRunDuration
        real(kind=8), parameter :: unsteadyAverageEndTf = unsteadyRunDuration
        real(kind=8), parameter :: unsteadyAverageMidTf = 0.5d0*(unsteadyAverageStartTf+unsteadyAverageEndTf)
        integer(kind=4), parameter :: outputSnapshotFile = 1    ! 0: 不输出快照；1: 输出，Nu/Re 仍独立采样
        integer(kind=4), parameter :: outputPltFile = 1         ! 0: 不输出 Tecplot；1: 输出
        integer(kind=4), parameter :: outputReloadFile = 1      ! 0: 不输出重启文件；1: 输出
        integer(kind=4), parameter :: itc_max = max(1, &
            ceiling(unsteadyRunDuration*timeUnit))  !ceiling 是向上取整函数，返回不小于输入值的最小整数
#endif
        ! 多块输出和最终时刻向上对齐到粗细同步步；各间隔至少为一个粗步。
        ! nextSample/nextReload/nextPlt 是输出时钟，和各文件编号分开；禁用某类文件不影响其他输出。

        ! 输出文件命名与格式版本
        character(*), parameter :: settingsFile = 'SimulationSettings2DOpenaccMultiblock.txt'
        character(*), parameter :: snapshotFilePrefix = 'buoyancyCavity2DOpenaccMultiblockSnapshot'
        character(*), parameter :: pltFolderPrefix = 'buoyancyCavity2DOpenaccMultiblockTecplot'
        character(*), parameter :: reloadFilePrefix = 'reloadFile2DOpenaccMultiblock'
        character(*), parameter :: NuReHistoryFile = 'NuRe_2DOpenaccMultiblock.dat'    ! 按采样时间记录 Nu、Re 及场量统计，续算时追加
        ! Magic 为二进制文件开头的格式标识；末尾数字表示格式版本，不是时间步或输出编号。
        ! 续算文件保存恢复计算所需的状态；读取时检查此标识，不匹配则停止，防止按错误格式读取数组。
        character(16), parameter :: restartMagic = 'MB2DRESTART0011'    ! 续算文件格式 v11
        ! 场快照用于绘图和后处理；读取程序可据此识别格式，不能把快照当作完整续算文件。
        character(16), parameter :: snapshotMagic = 'MB2DSNAPSHOT0003'    ! 场快照格式 v3

        !===============================================================================================
        ! 格子方向、块数据与接口交换数据
        !===============================================================================================
        integer(kind=4) :: ex(0:8) = [0, 1, 0, -1, 0, 1, -1, -1, 1]
        integer(kind=4) :: ey(0:8) = [0, 0, 1, 0, -1, 1, 1, -1, -1]
        real(kind=8) :: omega(0:8), omegaT(0:4)
        ! itc 为累计细时间步数，新计算从 0 开始，续算从检查点恢复；粗细比为 r 时，每个粗步累计增加 r。
        ! snapshotFileNum、pltFileNum 分别为快照和 Tecplot 文件编号，仅在实际输出相应文件时增加 1。
        integer(kind=4) :: itc = 0
        integer(kind=4) :: snapshotFileNum = 0, pltFileNum = 0
        ! nextSample、nextReload、nextPlt 分别安排下一次 Nu/Re 采样、续算文件和 Tecplot 输出。
        ! 目标时间为对应序号乘输出间隔，再向上对齐到粗细同步步；初值 1 表示第一次计划输出。
        ! 这些序号与文件编号分开：关闭快照仍按时采样，nextSample 继续增加，snapshotFileNum 不增加。
        ! 续算时，各计划序号和文件编号均从检查点恢复，不重新从 1 或 0 开始。
        integer(kind=4) :: nextSample = 1, nextReload = 1, nextPlt = 1
        real(kind=8) :: errorU = 100.0d0, errorT = 100.0d0

        ! 五套独立数组：中心粗网格，以及左、右、下、上细网格；不使用块编号。
        ! 左右细区贯穿全高，上下细区只填中间部分；角点不重复存储，中心空区不分配。
        ! 各场量的 i,j 从 1 开始，坐标为 xOffset+(i-0.5)*dx、yOffset+(j-0.5)*dx。
        ! f_post/g_post 另留 0 和 n+1 外围位置，用于相邻细区的迁移数据交换。
        real(kind=8) :: dxCoarse, SnuCoarse, SqCoarse, QkCoarse, QnuCoarse, gBetaCoarse
        real(kind=8) :: SnuFine, SqFine, QkFine, QnuFine
        real(kind=8), parameter :: dxFine = 1.0d0
        integer(kind=4) :: historyLastCoarse
        real(kind=8) :: centerBox(4)    ! 中心积分区 xmin,xmax,ymin,ymax；与重叠计算区分开。

        ! 中心粗网格：尺寸、坐标偏移、积分范围和物理墙面。
        integer(kind=4) :: nxCoarse, nyCoarse, iFirstCoarse, iLastCoarse, jFirstCoarse, jLastCoarse
        real(kind=8) :: xOffsetCoarse, yOffsetCoarse, ownedBoxCoarse(4)
        logical :: wallCoarse(4)    ! 左、右、下、上是否为物理墙面。
        real(kind=8), allocatable :: f_coarse(:, :, :)
        real(kind=8), allocatable :: g_coarse(:, :, :)
        real(kind=8), allocatable :: f_post_coarse(:, :, :)
        real(kind=8), allocatable :: g_post_coarse(:, :, :)
        real(kind=8), allocatable :: rho_coarse(:, :)
        real(kind=8), allocatable :: u_coarse(:, :)
        real(kind=8), allocatable :: v_coarse(:, :)
        real(kind=8), allocatable :: T_coarse(:, :)
        real(kind=8), allocatable :: Fx_coarse(:, :)
        real(kind=8), allocatable :: Fy_coarse(:, :)
        real(kind=8), allocatable :: rhoHistory_coarse(:, :, :)
        real(kind=8), allocatable :: uHistory_coarse(:, :, :)
        real(kind=8), allocatable :: vHistory_coarse(:, :, :)
        real(kind=8), allocatable :: THistory_coarse(:, :, :)
        real(kind=8), allocatable :: FxHistory_coarse(:, :, :)
        real(kind=8), allocatable :: FyHistory_coarse(:, :, :)
        real(kind=8), allocatable :: flowNeqHistory_coarse(:, :, :, :)
        real(kind=8), allocatable :: thermalNeqHistory_coarse(:, :, :, :)
        real(kind=8), allocatable :: quadWidthX_coarse(:)
        real(kind=8), allocatable :: quadWidthY_coarse(:)
#ifdef steadyFlow
        real(kind=8), allocatable :: up_coarse(:, :)
        real(kind=8), allocatable :: vp_coarse(:, :)
        real(kind=8), allocatable :: Tp_coarse(:, :)
#endif

        ! 左侧细网格：尺寸、坐标偏移、积分范围和物理墙面。
        integer(kind=4) :: nxLeft, nyLeft, iFirstLeft, iLastLeft, jFirstLeft, jLastLeft
        real(kind=8) :: xOffsetLeft, yOffsetLeft, ownedBoxLeft(4)
        logical :: wallLeft(4)    ! 左、右、下、上是否为物理墙面。
        real(kind=8), allocatable :: f_left(:, :, :)
        real(kind=8), allocatable :: g_left(:, :, :)
        real(kind=8), allocatable :: f_post_left(:, :, :)
        real(kind=8), allocatable :: g_post_left(:, :, :)
        real(kind=8), allocatable :: rho_left(:, :)
        real(kind=8), allocatable :: u_left(:, :)
        real(kind=8), allocatable :: v_left(:, :)
        real(kind=8), allocatable :: T_left(:, :)
        real(kind=8), allocatable :: Fx_left(:, :)
        real(kind=8), allocatable :: Fy_left(:, :)
        real(kind=8), allocatable :: rhoHistory_left(:, :, :)
        real(kind=8), allocatable :: uHistory_left(:, :, :)
        real(kind=8), allocatable :: vHistory_left(:, :, :)
        real(kind=8), allocatable :: THistory_left(:, :, :)
        real(kind=8), allocatable :: FxHistory_left(:, :, :)
        real(kind=8), allocatable :: FyHistory_left(:, :, :)
        real(kind=8), allocatable :: flowNeqHistory_left(:, :, :, :)
        real(kind=8), allocatable :: thermalNeqHistory_left(:, :, :, :)
        real(kind=8), allocatable :: quadWidthX_left(:)
        real(kind=8), allocatable :: quadWidthY_left(:)
#ifdef steadyFlow
        real(kind=8), allocatable :: up_left(:, :)
        real(kind=8), allocatable :: vp_left(:, :)
        real(kind=8), allocatable :: Tp_left(:, :)
#endif

        ! 右侧细网格：尺寸、坐标偏移、积分范围和物理墙面。
        integer(kind=4) :: nxRight, nyRight, iFirstRight, iLastRight, jFirstRight, jLastRight
        real(kind=8) :: xOffsetRight, yOffsetRight, ownedBoxRight(4)
        logical :: wallRight(4)    ! 左、右、下、上是否为物理墙面。
        real(kind=8), allocatable :: f_right(:, :, :)
        real(kind=8), allocatable :: g_right(:, :, :)
        real(kind=8), allocatable :: f_post_right(:, :, :)
        real(kind=8), allocatable :: g_post_right(:, :, :)
        real(kind=8), allocatable :: rho_right(:, :)
        real(kind=8), allocatable :: u_right(:, :)
        real(kind=8), allocatable :: v_right(:, :)
        real(kind=8), allocatable :: T_right(:, :)
        real(kind=8), allocatable :: Fx_right(:, :)
        real(kind=8), allocatable :: Fy_right(:, :)
        real(kind=8), allocatable :: rhoHistory_right(:, :, :)
        real(kind=8), allocatable :: uHistory_right(:, :, :)
        real(kind=8), allocatable :: vHistory_right(:, :, :)
        real(kind=8), allocatable :: THistory_right(:, :, :)
        real(kind=8), allocatable :: FxHistory_right(:, :, :)
        real(kind=8), allocatable :: FyHistory_right(:, :, :)
        real(kind=8), allocatable :: flowNeqHistory_right(:, :, :, :)
        real(kind=8), allocatable :: thermalNeqHistory_right(:, :, :, :)
        real(kind=8), allocatable :: quadWidthX_right(:)
        real(kind=8), allocatable :: quadWidthY_right(:)
#ifdef steadyFlow
        real(kind=8), allocatable :: up_right(:, :)
        real(kind=8), allocatable :: vp_right(:, :)
        real(kind=8), allocatable :: Tp_right(:, :)
#endif

        ! 下侧细网格：尺寸、坐标偏移、积分范围和物理墙面。
        integer(kind=4) :: nxBottom, nyBottom, iFirstBottom, iLastBottom, jFirstBottom, jLastBottom
        real(kind=8) :: xOffsetBottom, yOffsetBottom, ownedBoxBottom(4)
        logical :: wallBottom(4)    ! 左、右、下、上是否为物理墙面。
        real(kind=8), allocatable :: f_bottom(:, :, :)
        real(kind=8), allocatable :: g_bottom(:, :, :)
        real(kind=8), allocatable :: f_post_bottom(:, :, :)
        real(kind=8), allocatable :: g_post_bottom(:, :, :)
        real(kind=8), allocatable :: rho_bottom(:, :)
        real(kind=8), allocatable :: u_bottom(:, :)
        real(kind=8), allocatable :: v_bottom(:, :)
        real(kind=8), allocatable :: T_bottom(:, :)
        real(kind=8), allocatable :: Fx_bottom(:, :)
        real(kind=8), allocatable :: Fy_bottom(:, :)
        real(kind=8), allocatable :: rhoHistory_bottom(:, :, :)
        real(kind=8), allocatable :: uHistory_bottom(:, :, :)
        real(kind=8), allocatable :: vHistory_bottom(:, :, :)
        real(kind=8), allocatable :: THistory_bottom(:, :, :)
        real(kind=8), allocatable :: FxHistory_bottom(:, :, :)
        real(kind=8), allocatable :: FyHistory_bottom(:, :, :)
        real(kind=8), allocatable :: flowNeqHistory_bottom(:, :, :, :)
        real(kind=8), allocatable :: thermalNeqHistory_bottom(:, :, :, :)
        real(kind=8), allocatable :: quadWidthX_bottom(:)
        real(kind=8), allocatable :: quadWidthY_bottom(:)
#ifdef steadyFlow
        real(kind=8), allocatable :: up_bottom(:, :)
        real(kind=8), allocatable :: vp_bottom(:, :)
        real(kind=8), allocatable :: Tp_bottom(:, :)
#endif

        ! 上侧细网格：尺寸、坐标偏移、积分范围和物理墙面。
        integer(kind=4) :: nxTop, nyTop, iFirstTop, iLastTop, jFirstTop, jLastTop
        real(kind=8) :: xOffsetTop, yOffsetTop, ownedBoxTop(4)
        logical :: wallTop(4)    ! 左、右、下、上是否为物理墙面。
        real(kind=8), allocatable :: f_top(:, :, :)
        real(kind=8), allocatable :: g_top(:, :, :)
        real(kind=8), allocatable :: f_post_top(:, :, :)
        real(kind=8), allocatable :: g_post_top(:, :, :)
        real(kind=8), allocatable :: rho_top(:, :)
        real(kind=8), allocatable :: u_top(:, :)
        real(kind=8), allocatable :: v_top(:, :)
        real(kind=8), allocatable :: T_top(:, :)
        real(kind=8), allocatable :: Fx_top(:, :)
        real(kind=8), allocatable :: Fy_top(:, :)
        real(kind=8), allocatable :: rhoHistory_top(:, :, :)
        real(kind=8), allocatable :: uHistory_top(:, :, :)
        real(kind=8), allocatable :: vHistory_top(:, :, :)
        real(kind=8), allocatable :: THistory_top(:, :, :)
        real(kind=8), allocatable :: FxHistory_top(:, :, :)
        real(kind=8), allocatable :: FyHistory_top(:, :, :)
        real(kind=8), allocatable :: flowNeqHistory_top(:, :, :, :)
        real(kind=8), allocatable :: thermalNeqHistory_top(:, :, :, :)
        real(kind=8), allocatable :: quadWidthX_top(:)
        real(kind=8), allocatable :: quadWidthY_top(:)
#ifdef steadyFlow
        real(kind=8), allocatable :: up_top(:, :)
        real(kind=8), allocatable :: vp_top(:, :)
        real(kind=8), allocatable :: Tp_top(:, :)
#endif

        ! History 保存接口交换所需的宏观量、单位时间力和缩放非平衡矩。
        ! 粗网格末下标 0:2 为过去、当前、预测时间层；四套细网格均只保存当前层 0。
        ! 力历史为 Fx/dx、Fy/dx：本程序 dt=dx，接收端再乘自身 dx 恢复每步力增量。
        ! 各条连接直接以方向命名；Ti/Tj 为接收节点，Si/Sj 为来源模板，不再保存块编号。
        integer(kind=4) :: coarseToLeftCount
        integer(kind=4), allocatable :: coarseToLeftTi(:), coarseToLeftTj(:), coarseToLeftSi(:), coarseToLeftSj(:)
        logical, allocatable :: coarseToLeftSame(:)
        real(kind=8), allocatable :: coarseToLeftWx(:, :), coarseToLeftWy(:, :)
        integer(kind=4) :: leftToCoarseCount
        integer(kind=4), allocatable :: leftToCoarseTi(:), leftToCoarseTj(:), leftToCoarseSi(:), leftToCoarseSj(:)
        logical, allocatable :: leftToCoarseSame(:)
        real(kind=8), allocatable :: leftToCoarseWx(:, :), leftToCoarseWy(:, :)
        integer(kind=4) :: coarseToRightCount
        integer(kind=4), allocatable :: coarseToRightTi(:), coarseToRightTj(:), coarseToRightSi(:), &
            coarseToRightSj(:)
        logical, allocatable :: coarseToRightSame(:)
        real(kind=8), allocatable :: coarseToRightWx(:, :), coarseToRightWy(:, :)
        integer(kind=4) :: rightToCoarseCount
        integer(kind=4), allocatable :: rightToCoarseTi(:), rightToCoarseTj(:), rightToCoarseSi(:), &
            rightToCoarseSj(:)
        logical, allocatable :: rightToCoarseSame(:)
        real(kind=8), allocatable :: rightToCoarseWx(:, :), rightToCoarseWy(:, :)
        integer(kind=4) :: coarseToBottomCount
        integer(kind=4), allocatable :: coarseToBottomTi(:), coarseToBottomTj(:), coarseToBottomSi(:), &
            coarseToBottomSj(:)
        logical, allocatable :: coarseToBottomSame(:)
        real(kind=8), allocatable :: coarseToBottomWx(:, :), coarseToBottomWy(:, :)
        integer(kind=4) :: bottomToCoarseCount
        integer(kind=4), allocatable :: bottomToCoarseTi(:), bottomToCoarseTj(:), bottomToCoarseSi(:), &
            bottomToCoarseSj(:)
        logical, allocatable :: bottomToCoarseSame(:)
        real(kind=8), allocatable :: bottomToCoarseWx(:, :), bottomToCoarseWy(:, :)
        integer(kind=4) :: coarseToTopCount
        integer(kind=4), allocatable :: coarseToTopTi(:), coarseToTopTj(:), coarseToTopSi(:), coarseToTopSj(:)
        logical, allocatable :: coarseToTopSame(:)
        real(kind=8), allocatable :: coarseToTopWx(:, :), coarseToTopWy(:, :)
        integer(kind=4) :: topToCoarseCount
        integer(kind=4), allocatable :: topToCoarseTi(:), topToCoarseTj(:), topToCoarseSi(:), topToCoarseSj(:)
        logical, allocatable :: topToCoarseSame(:)
        real(kind=8), allocatable :: topToCoarseWx(:, :), topToCoarseWy(:, :)
        ! 各方向依次交换，共用这一组接收临时数组；大小只取最长接口的节点数。
        real(kind=8), allocatable :: rhoReceive(:), uReceive(:), vReceive(:), TReceive(:)
        real(kind=8), allocatable :: FxReceive(:), FyReceive(:)
        real(kind=8), allocatable :: flowNeqReceive(:, :), thermalNeqReceive(:, :)
        real(kind=8), external :: owned_cell_area, section_owned_weight
        logical, external :: coarse_skin, fine_skin
        integer(kind=4), external :: scheduled_step
    end module commondata


    program main

        use openacc
        use commondata
        implicit none

        integer(kind=4) :: finalStep
        integer(kind=8) :: clockStart, clockEnd, clockRate

        ! 初始化设备、网格及场量；续算时在 initial 中恢复检查点。
        call acc_init(acc_device_default)
        write(*, *) 'Visible OpenACC devices:', acc_get_num_devices(acc_device_default)
        call initial()
        call enter_data_2d_openacc()
        call system_clock(clockStart, clockRate)
        ! itc 按最细时间步累计，停止时刻向上对齐到粗细块同步点。
        finalStep = ((itc_max+refineRatio-1)/refineRatio)*refineRatio
        do while (itc < finalStep)
#ifdef steadyFlow
            if (errorU <= epsU .and. errorT <= epsT) exit
#endif
            ! 一次循环推进一个粗步，内部执行 refineRatio 个细步及接口交换。
            call advance_multiblock()
#ifdef steadyFlow
            ! 稳态模式：定期检查速度和温度误差，两者均达标后停止。
            if (mod(itc, 2000) == 0) call check()
#endif
            ! Nu/Re 采样、Tecplot 和续算文件各用独立时钟，关闭某类文件不停止采样。
            if (itc >= scheduled_step(nextSample, outputSnapshotInterval)) then
                call calNuRe()
                nextSample = nextSample+1
                if (outputSnapshotFile == 1) call output_SnapshotFile()
            endif
            if (itc >= scheduled_step(nextPlt, outputPltFileInterval)) then
                nextPlt = nextPlt+1
                if (outputPltFile == 1) call output_Tecplot()
            endif
            if (itc >= scheduled_step(nextReload, reloadFileInterval)) then
                nextReload = nextReload+1
                if (outputReloadFile == 1) call output_ReloadFile()
            endif
        enddo
        !$acc wait(1)
        call system_clock(clockEnd)
        write(*, *) 'Elapsed seconds:', dble(clockEnd-clockStart)/dble(clockRate)
        ! 非稳态模式：按指定窗口统计 Nu/Re，并比较前后半窗口；稳态时该调用为空。
        call output_unsteady_NuRe_postprocess()
        if (outputPltFile == 1) call output_Tecplot()
        if (outputReloadFile == 1) call output_ReloadFile()
        call exit_data_2d_openacc()
    end program main


    ! 检查参数，建立五个具名网格的尺寸与偏移，初始化或读取续算数据。
    subroutine initial()

        use commondata
        implicit none

        integer(kind=4) :: k, overlap, i, j
        real(kind=8) :: totalArea, xLeft, xRight, yBottom, yTop

        if (refineRatio < 1) error stop 'refineRatio must be a positive integer'
        if (min(nx, ny) < 8) error stop 'At least 8 cells per direction are required'
        if (loadInitField /= 0 .and. loadInitField /= 1) error stop 'loadInitField must be 0 or 1'
        if (loadInitField == 0) reloadFileNum = 0
        if (min(outputSnapshotInterval, reloadFileInterval, outputPltFileInterval)*timeUnit < dble(refineRatio)) &
            error stop 'Output intervals must be at least one synchronized coarse step'
        if (paraA <= -4.0d0 .or. paraA >= 1.0d0) error stop 'Legacy paraA must be in (-4,1)'
        omega(0) = 4.0d0/9.0d0
        omega(1:4) = 1.0d0/9.0d0
        omega(5:8) = 1.0d0/36.0d0
        omegaT(0) = (1.0d0-thermalA)/5.0d0
        omegaT(1:4) = (thermalA+4.0d0)/20.0d0
        dxCoarse = dble(refineRatio)
        historyLastCoarse = 0
        if (refineRatio == 1) then
            centerBox = [0.0d0, dble(nx), 0.0d0, dble(ny)]
            nxCoarse = nx
            nyCoarse = ny
            xOffsetCoarse = 0.0d0
            yOffsetCoarse = 0.0d0
        else
#if defined(VerticalWallsPeriodicalU) || defined(VerticalWallsPeriodicalT)
            error stop 'Multiblock periodic sides require a periodic block topology; use wall BCs or refineRatio=1'
#endif
            overlap = overlapCells*refineRatio
            if (interfaceSkin < 2) error stop 'Split flow/thermal advance requires at least two interface layers'
            ! 目标细块最内侧缓冲点必须位于来源粗块可用区内；四点模板另由 donor_stencil 检查。
            ! 2*overlap >= interfaceSkin*refineRatio + interfaceSkin-1；默认比值2、两层缓冲时 overlapCells>=2。
            if (2*overlap < interfaceSkin*refineRatio+interfaceSkin-1) &
                error stop 'Overlap is too narrow for the pre-collision interface layers'
            if (mod(nx, refineRatio) /= 0 .or. mod(ny, refineRatio) /= 0) &
                error stop 'nx,ny must be multiples of refineRatio'
            xLeft = dble(fineLayerCellsLeft)-0.5d0
            xRight = dble(nx-fineLayerCellsRight)+0.5d0
            yBottom = dble(fineLayerCellsBottom)-0.5d0
            yTop = dble(ny-fineLayerCellsTop)+0.5d0
            if (min(xLeft, dble(nx)-xRight, yBottom, dble(ny)-yTop) < dble(overlap+refineRatio)) &
                error stop 'Refined wall layer must exceed the overlap width by at least one coarse cell'
            if (xRight <= xLeft .or. yTop <= yBottom) &
                error stop 'The coarse core must have positive width and height'
            if (mod(nx-fineLayerCellsLeft-fineLayerCellsRight+1, refineRatio) /= 0 .or. &
                mod(ny-fineLayerCellsBottom-fineLayerCellsTop+1, refineRatio) /= 0) &
                error stop 'Central width and height must be multiples of refineRatio; adjust the four interface indices'

            historyLastCoarse = 2
            centerBox = [xLeft, xRight, yBottom, yTop]
            nxCoarse = nint((xRight-xLeft)/dxCoarse)+2*overlapCells+1
            nyCoarse = nint((yTop-yBottom)/dxCoarse)+2*overlapCells+1
            xOffsetCoarse = xLeft-overlap-0.5d0*dxCoarse
            yOffsetCoarse = yBottom-overlap-0.5d0*dxCoarse

            ! 四个细矩形恰好拼成原来的有效细环；拼接处不增加插值接口。
            nxLeft = fineLayerCellsLeft+overlap
            nyLeft = ny
            xOffsetLeft = 0.0d0
            yOffsetLeft = 0.0d0
            nxRight = fineLayerCellsRight+overlap
            nyRight = ny
            xOffsetRight = dble(nx-nxRight)
            yOffsetRight = 0.0d0
            nxBottom = nx-nxLeft-nxRight
            nyBottom = fineLayerCellsBottom+overlap
            xOffsetBottom = dble(nxLeft)
            yOffsetBottom = 0.0d0
            nxTop = nxBottom
            nyTop = fineLayerCellsTop+overlap
            xOffsetTop = dble(nxLeft)
            yOffsetTop = dble(ny-nyTop)
            if (min(nxBottom, ny-nyBottom-nyTop) < 1) &
                error stop 'Compact fine ring requires a positive interior storage hole'
        endif
        ownedBoxCoarse = centerBox
        wallCoarse = refineRatio == 1
        SnuCoarse = 1.0d0/(0.5d0+(1.0d0/Snu-0.5d0)/dxCoarse)
        SqCoarse = 1.0d0/(0.5d0+(1.0d0/Sq-0.5d0)/dxCoarse)
        QkCoarse = 1.0d0/(0.5d0+(1.0d0/Qk-0.5d0)/dxCoarse)
        QnuCoarse = 1.0d0/(0.5d0+(1.0d0/Qnu-0.5d0)/dxCoarse)
        gBetaCoarse = dxCoarse*gBeta
        ! 细网格也沿用同一缩放表达式，保持重排存储前后的浮点运算结果一致。
        SnuFine = 1.0d0/(0.5d0+(1.0d0/Snu-0.5d0)/dxFine)
        SqFine = 1.0d0/(0.5d0+(1.0d0/Sq-0.5d0)/dxFine)
        QkFine = 1.0d0/(0.5d0+(1.0d0/Qk-0.5d0)/dxFine)
        QnuFine = 1.0d0/(0.5d0+(1.0d0/Qnu-0.5d0)/dxFine)
        if (refineRatio > 1) then
            ownedBoxLeft = [0.0d0, dble(nxLeft), 0.0d0, dble(ny)]
            ownedBoxRight = [xOffsetRight, dble(nx), 0.0d0, dble(ny)]
            ownedBoxBottom = [xOffsetBottom, xOffsetRight, 0.0d0, dble(nyBottom)]
            ownedBoxTop = [xOffsetTop, xOffsetRight, yOffsetTop, dble(ny)]
            wallLeft = [.true., .false., .true., .true.]
            wallRight = [.false., .true., .true., .true.]
            wallBottom = [.false., .false., .true., .false.]
            wallTop = [.false., .false., .false., .true.]
        endif
        call allocate_grid_arrays()
        totalArea = 0.0d0
        call initial_grid(nxCoarse, nyCoarse, historyLastCoarse, dxCoarse, xOffsetCoarse, yOffsetCoarse, &
            iFirstCoarse, iLastCoarse, jFirstCoarse, jLastCoarse, ownedBoxCoarse, f_coarse, g_coarse, &
            f_post_coarse, g_post_coarse, rho_coarse, u_coarse, v_coarse, T_coarse, Fx_coarse, Fy_coarse, &
            rhoHistory_coarse, uHistory_coarse, vHistory_coarse, THistory_coarse, FxHistory_coarse, &
            FyHistory_coarse, flowNeqHistory_coarse, thermalNeqHistory_coarse, quadWidthX_coarse, &
            quadWidthY_coarse &
#ifdef steadyFlow
            , up_coarse, vp_coarse, Tp_coarse &
#endif
        )
        do j = 1, nyCoarse
            do i = 1, nxCoarse
                totalArea = totalArea+owned_cell_area(xOffsetCoarse+(i-0.5d0)*dxCoarse, &
                    yOffsetCoarse+(j-0.5d0)*dxCoarse, dxCoarse, ownedBoxCoarse)
            enddo
        enddo
        if (refineRatio > 1) then
            call initial_grid(nxLeft, nyLeft, 0, dxFine, xOffsetLeft, yOffsetLeft, iFirstLeft, iLastLeft, &
                jFirstLeft, jLastLeft, ownedBoxLeft, f_left, g_left, f_post_left, g_post_left, rho_left, &
                u_left, v_left, T_left, Fx_left, Fy_left, rhoHistory_left, uHistory_left, vHistory_left, &
                THistory_left, FxHistory_left, FyHistory_left, flowNeqHistory_left, thermalNeqHistory_left, &
                quadWidthX_left, quadWidthY_left &
#ifdef steadyFlow
                , up_left, vp_left, Tp_left &
#endif
            )
            do j = 1, nyLeft
                do i = 1, nxLeft
                    totalArea = totalArea+owned_cell_area(xOffsetLeft+(i-0.5d0)*dxFine, &
                        yOffsetLeft+(j-0.5d0)*dxFine, dxFine, ownedBoxLeft)
                enddo
            enddo
            call initial_grid(nxRight, nyRight, 0, dxFine, xOffsetRight, yOffsetRight, iFirstRight, iLastRight, &
                jFirstRight, jLastRight, ownedBoxRight, f_right, g_right, f_post_right, g_post_right, &
                rho_right, u_right, v_right, T_right, Fx_right, Fy_right, rhoHistory_right, uHistory_right, &
                vHistory_right, THistory_right, FxHistory_right, FyHistory_right, flowNeqHistory_right, &
                thermalNeqHistory_right, quadWidthX_right, quadWidthY_right &
#ifdef steadyFlow
                , up_right, vp_right, Tp_right &
#endif
            )
            do j = 1, nyRight
                do i = 1, nxRight
                    totalArea = totalArea+owned_cell_area(xOffsetRight+(i-0.5d0)*dxFine, &
                        yOffsetRight+(j-0.5d0)*dxFine, dxFine, ownedBoxRight)
                enddo
            enddo
            call initial_grid(nxBottom, nyBottom, 0, dxFine, xOffsetBottom, yOffsetBottom, iFirstBottom, &
                iLastBottom, jFirstBottom, jLastBottom, ownedBoxBottom, f_bottom, g_bottom, f_post_bottom, &
                g_post_bottom, rho_bottom, u_bottom, v_bottom, T_bottom, Fx_bottom, Fy_bottom, &
                rhoHistory_bottom, uHistory_bottom, vHistory_bottom, THistory_bottom, FxHistory_bottom, &
                FyHistory_bottom, flowNeqHistory_bottom, thermalNeqHistory_bottom, quadWidthX_bottom, &
                quadWidthY_bottom &
#ifdef steadyFlow
                , up_bottom, vp_bottom, Tp_bottom &
#endif
            )
            do j = 1, nyBottom
                do i = 1, nxBottom
                    totalArea = totalArea+owned_cell_area(xOffsetBottom+(i-0.5d0)*dxFine, &
                        yOffsetBottom+(j-0.5d0)*dxFine, dxFine, ownedBoxBottom)
                enddo
            enddo
            call initial_grid(nxTop, nyTop, 0, dxFine, xOffsetTop, yOffsetTop, iFirstTop, iLastTop, jFirstTop, &
                jLastTop, ownedBoxTop, f_top, g_top, f_post_top, g_post_top, rho_top, u_top, v_top, T_top, &
                Fx_top, Fy_top, rhoHistory_top, uHistory_top, vHistory_top, THistory_top, FxHistory_top, &
                FyHistory_top, flowNeqHistory_top, thermalNeqHistory_top, quadWidthX_top, quadWidthY_top &
#ifdef steadyFlow
                , up_top, vp_top, Tp_top &
#endif
            )
            do j = 1, nyTop
                do i = 1, nxTop
                    totalArea = totalArea+owned_cell_area(xOffsetTop+(i-0.5d0)*dxFine, &
                        yOffsetTop+(j-0.5d0)*dxFine, dxFine, ownedBoxTop)
                enddo
            enddo
        endif
        if (abs(totalArea-dble(nx)*ny) > 1.0d-8) error stop 'Grid ownership does not tile the physical domain'
        if (loadInitField == 1) call read_restart()
        call build_interfaces()
        open(newunit = k, file = settingsFile, status = 'replace')
        write(k, *) 'Parent: uniform-grid 2DRBOpenacc.F90; D2Q9 MRT / D2Q5 MRT unchanged inside blocks'
        write(k, *) 'Huang (2014): force-aware MRT scaling and four-point Lagrange; original D2Q5 adaptation.'
        write(k, *) 'Yu (2002): aligned two-way coupling and acoustic subcycling; cubic spline is not used.'
        write(k, *) 'Chen (2016): pre-collision buffer states; received f,g participate in the next collision.'
        write(k, *) 'Fine-equivalent nx,ny; refinement ratio:', nx, ny, refineRatio
        write(k, *) 'Rayleigh, Prandtl, Mach:', Rayleigh, Prandtl, Mach
        write(k, *) 'Fine tauf, Snu, Sq, Qk, Qnu, thermalA:', tauf, Snu, Sq, Qk, Qnu, thermalA
        write(k, *) 'Fine viscosity, diffusivity, gBeta, timeUnit:', viscosity, diffusivity, gBeta, timeUnit
        write(k, *) 'Fine-node indices left/right/bottom/top; overlap in coarse spacings:', &
            fineLayerCellsLeft, fineLayerCellsRight, fineLayerCellsBottom, fineLayerCellsTop, overlapCells
        write(k, *) 'Owned physical area:', totalArea
        write(k, *) 'Geometry: central coarse rectangle and four compact fine rectangles.'
        write(k, *) 'Fine seams exchange post-collision populations directly, including diagonal links.'
        write(k, *) 'Coordinates: x=xOffset+(i-0.5)*dx, y=yOffset+(j-0.5)*dx; local i,j start at 1.'
        write(k, *) 'Interfaces preserve the connected-ring time interpolation and MRT scaling; dt=dx.'
        write(k, *) 'Restart v11 stores named coarse,left,right,bottom,top arrays; older formats need conversion.'
        write(k, *) 'coarse: nx,ny,dx,xOffset,yOffset:', nxCoarse, nyCoarse, dxCoarse, xOffsetCoarse, yOffsetCoarse
        if (refineRatio > 1) then
            write(k, *) 'left: nx,ny,dx,xOffset,yOffset:', nxLeft, nyLeft, dxFine, xOffsetLeft, yOffsetLeft
            write(k, *) 'right: nx,ny,dx,xOffset,yOffset:', nxRight, nyRight, dxFine, xOffsetRight, yOffsetRight
            write(k, *) 'bottom: nx,ny,dx,xOffset,yOffset:', nxBottom, nyBottom, dxFine, xOffsetBottom, &
            yOffsetBottom
            write(k, *) 'top: nx,ny,dx,xOffset,yOffset:', nxTop, nyTop, dxFine, xOffsetTop, yOffsetTop
        endif
        close(k)
        if (loadInitField == 0) then
            open(newunit = k, file = NuReHistoryFile, status = 'replace')
            write(k, '(a)') '# t_ff NuVolAvg ReVolRMS Nu_hot Nu_cold Nu_middle mass meanT Tmin Tmax rhoMin rhoMax'
            close(k)
        else
            call check_history()
        endif
    end subroutine initial


    !===============================================================================================
    ! 按各矩形实际尺寸分配普通数组；不分配完整 nx*ny 细网格。
    subroutine allocate_grid_arrays()
        use commondata
        implicit none

        allocate(f_coarse(nxCoarse, nyCoarse, 0:8))
        allocate(g_coarse(nxCoarse, nyCoarse, 0:4))
        allocate(f_post_coarse(0:nxCoarse+1, 0:nyCoarse+1, 0:8))
        allocate(g_post_coarse(0:nxCoarse+1, 0:nyCoarse+1, 0:4))
        allocate(rho_coarse(nxCoarse, nyCoarse))
        allocate(u_coarse(nxCoarse, nyCoarse))
        allocate(v_coarse(nxCoarse, nyCoarse))
        allocate(T_coarse(nxCoarse, nyCoarse))
        allocate(Fx_coarse(nxCoarse, nyCoarse))
        allocate(Fy_coarse(nxCoarse, nyCoarse))
        allocate(rhoHistory_coarse(nxCoarse, nyCoarse, 0:historyLastCoarse))
        allocate(uHistory_coarse(nxCoarse, nyCoarse, 0:historyLastCoarse))
        allocate(vHistory_coarse(nxCoarse, nyCoarse, 0:historyLastCoarse))
        allocate(THistory_coarse(nxCoarse, nyCoarse, 0:historyLastCoarse))
        allocate(FxHistory_coarse(nxCoarse, nyCoarse, 0:historyLastCoarse))
        allocate(FyHistory_coarse(nxCoarse, nyCoarse, 0:historyLastCoarse))
        allocate(flowNeqHistory_coarse(nxCoarse, nyCoarse, 0:8, 0:historyLastCoarse))
        allocate(thermalNeqHistory_coarse(nxCoarse, nyCoarse, 0:4, 0:historyLastCoarse))
        allocate(quadWidthX_coarse(nxCoarse))
        allocate(quadWidthY_coarse(nyCoarse))
#ifdef steadyFlow
        allocate(up_coarse(nxCoarse, nyCoarse))
        allocate(vp_coarse(nxCoarse, nyCoarse))
        allocate(Tp_coarse(nxCoarse, nyCoarse))
#endif
        if (refineRatio > 1) then
            allocate(f_left(nxLeft, nyLeft, 0:8))
            allocate(g_left(nxLeft, nyLeft, 0:4))
            allocate(f_post_left(0:nxLeft+1, 0:nyLeft+1, 0:8))
            allocate(g_post_left(0:nxLeft+1, 0:nyLeft+1, 0:4))
            allocate(rho_left(nxLeft, nyLeft))
            allocate(u_left(nxLeft, nyLeft))
            allocate(v_left(nxLeft, nyLeft))
            allocate(T_left(nxLeft, nyLeft))
            allocate(Fx_left(nxLeft, nyLeft))
            allocate(Fy_left(nxLeft, nyLeft))
            allocate(rhoHistory_left(nxLeft, nyLeft, 0:0))
            allocate(uHistory_left(nxLeft, nyLeft, 0:0))
            allocate(vHistory_left(nxLeft, nyLeft, 0:0))
            allocate(THistory_left(nxLeft, nyLeft, 0:0))
            allocate(FxHistory_left(nxLeft, nyLeft, 0:0))
            allocate(FyHistory_left(nxLeft, nyLeft, 0:0))
            allocate(flowNeqHistory_left(nxLeft, nyLeft, 0:8, 0:0))
            allocate(thermalNeqHistory_left(nxLeft, nyLeft, 0:4, 0:0))
            allocate(quadWidthX_left(nxLeft))
            allocate(quadWidthY_left(nyLeft))
#ifdef steadyFlow
            allocate(up_left(nxLeft, nyLeft))
            allocate(vp_left(nxLeft, nyLeft))
            allocate(Tp_left(nxLeft, nyLeft))
#endif
            allocate(f_right(nxRight, nyRight, 0:8))
            allocate(g_right(nxRight, nyRight, 0:4))
            allocate(f_post_right(0:nxRight+1, 0:nyRight+1, 0:8))
            allocate(g_post_right(0:nxRight+1, 0:nyRight+1, 0:4))
            allocate(rho_right(nxRight, nyRight))
            allocate(u_right(nxRight, nyRight))
            allocate(v_right(nxRight, nyRight))
            allocate(T_right(nxRight, nyRight))
            allocate(Fx_right(nxRight, nyRight))
            allocate(Fy_right(nxRight, nyRight))
            allocate(rhoHistory_right(nxRight, nyRight, 0:0))
            allocate(uHistory_right(nxRight, nyRight, 0:0))
            allocate(vHistory_right(nxRight, nyRight, 0:0))
            allocate(THistory_right(nxRight, nyRight, 0:0))
            allocate(FxHistory_right(nxRight, nyRight, 0:0))
            allocate(FyHistory_right(nxRight, nyRight, 0:0))
            allocate(flowNeqHistory_right(nxRight, nyRight, 0:8, 0:0))
            allocate(thermalNeqHistory_right(nxRight, nyRight, 0:4, 0:0))
            allocate(quadWidthX_right(nxRight))
            allocate(quadWidthY_right(nyRight))
#ifdef steadyFlow
            allocate(up_right(nxRight, nyRight))
            allocate(vp_right(nxRight, nyRight))
            allocate(Tp_right(nxRight, nyRight))
#endif
            allocate(f_bottom(nxBottom, nyBottom, 0:8))
            allocate(g_bottom(nxBottom, nyBottom, 0:4))
            allocate(f_post_bottom(0:nxBottom+1, 0:nyBottom+1, 0:8))
            allocate(g_post_bottom(0:nxBottom+1, 0:nyBottom+1, 0:4))
            allocate(rho_bottom(nxBottom, nyBottom))
            allocate(u_bottom(nxBottom, nyBottom))
            allocate(v_bottom(nxBottom, nyBottom))
            allocate(T_bottom(nxBottom, nyBottom))
            allocate(Fx_bottom(nxBottom, nyBottom))
            allocate(Fy_bottom(nxBottom, nyBottom))
            allocate(rhoHistory_bottom(nxBottom, nyBottom, 0:0))
            allocate(uHistory_bottom(nxBottom, nyBottom, 0:0))
            allocate(vHistory_bottom(nxBottom, nyBottom, 0:0))
            allocate(THistory_bottom(nxBottom, nyBottom, 0:0))
            allocate(FxHistory_bottom(nxBottom, nyBottom, 0:0))
            allocate(FyHistory_bottom(nxBottom, nyBottom, 0:0))
            allocate(flowNeqHistory_bottom(nxBottom, nyBottom, 0:8, 0:0))
            allocate(thermalNeqHistory_bottom(nxBottom, nyBottom, 0:4, 0:0))
            allocate(quadWidthX_bottom(nxBottom))
            allocate(quadWidthY_bottom(nyBottom))
#ifdef steadyFlow
            allocate(up_bottom(nxBottom, nyBottom))
            allocate(vp_bottom(nxBottom, nyBottom))
            allocate(Tp_bottom(nxBottom, nyBottom))
#endif
            allocate(f_top(nxTop, nyTop, 0:8))
            allocate(g_top(nxTop, nyTop, 0:4))
            allocate(f_post_top(0:nxTop+1, 0:nyTop+1, 0:8))
            allocate(g_post_top(0:nxTop+1, 0:nyTop+1, 0:4))
            allocate(rho_top(nxTop, nyTop))
            allocate(u_top(nxTop, nyTop))
            allocate(v_top(nxTop, nyTop))
            allocate(T_top(nxTop, nyTop))
            allocate(Fx_top(nxTop, nyTop))
            allocate(Fy_top(nxTop, nyTop))
            allocate(rhoHistory_top(nxTop, nyTop, 0:0))
            allocate(uHistory_top(nxTop, nyTop, 0:0))
            allocate(vHistory_top(nxTop, nyTop, 0:0))
            allocate(THistory_top(nxTop, nyTop, 0:0))
            allocate(FxHistory_top(nxTop, nyTop, 0:0))
            allocate(FyHistory_top(nxTop, nyTop, 0:0))
            allocate(flowNeqHistory_top(nxTop, nyTop, 0:8, 0:0))
            allocate(thermalNeqHistory_top(nxTop, nyTop, 0:4, 0:0))
            allocate(quadWidthX_top(nxTop))
            allocate(quadWidthY_top(nyTop))
#ifdef steadyFlow
            allocate(up_top(nxTop, nyTop))
            allocate(vp_top(nxTop, nyTop))
            allocate(Tp_top(nxTop, nyTop))
#endif
        endif
    end subroutine allocate_grid_arrays


    subroutine integration_weights(n, origin, dx, lo, hi, w, first, last)

        implicit none

        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: origin, dx, lo, hi
        real(kind=8), intent(out) :: w(n)
        integer(kind=4), intent(out) :: first, last
        integer(kind=4) :: i
        real(kind=8) :: firstMoment, x, shift

        first = n+1
        last = 0
        firstMoment = 0.0d0
        do i = 1, n
            w(i) = max(0.0d0, min(hi, origin+dble(i)*dx)-max(lo, origin+dble(i-1)*dx))
            if (w(i) <= 0.0d0) cycle
            first = min(first, i)
            last = i
            x = origin+(dble(i)-0.5d0)*dx
            firstMoment = firstMoment+w(i)*x
        enddo
        if (abs(sum(w)-(hi-lo)) > 1.0d-10) error stop 'Integration weights do not cover owned interval'
        ! 中心粗块在共址节点之间积分，无需一阶矩修正；靠墙细块的裁剪仍可能不对称。
        ! 在末两个积分节点之间转移权重，补齐一阶矩，保持总长度。
        ! 两节点相距 dx，转移 shift 后一阶矩增加 dx*shift；无需修正的分区保持 shift=0。
        if (last > first) then
            shift = (0.5d0*(hi**2-lo**2)-firstMoment)/dx
            w(last-1) = w(last-1)-shift
            w(last) = w(last)+shift
        endif
        if (any(w < 0.0d0)) error stop 'Negative corrected integration weight'
        firstMoment = 0.0d0
        do i = first, last
            firstMoment = firstMoment+w(i)*(origin+(dble(i)-0.5d0)*dx)
        enddo
        if (abs(firstMoment-0.5d0*(hi**2-lo**2)) > 1.0d-9*max(1.0d0, abs(firstMoment))) &
            error stop 'Integration weights do not integrate a linear coordinate exactly'
    end subroutine integration_weights


    ! 直接初始化传入的数组；ni,nj 为本地尺寸，dx 和偏移决定节点坐标。
    subroutine initial_grid(ni, nj, nh, dx, xOffset, yOffset, iFirst, iLast, jFirst, jLast, ownedBox, f, g, &
            f_post, g_post, rho, u, v, T, Fx, Fy, rhoHistory, uHistory, &
        vHistory, THistory, FxHistory, FyHistory, flowNeqHistory, thermalNeqHistory, quadWidthX, quadWidthY &
#ifdef steadyFlow
            , up, vp, Tp &
#endif
        )

        use commondata, only: nx, ny, Thot, Tcold, pi, lengthUnit, omega, omegaT
        implicit none
        integer(kind=4) :: ni, nj, nh, iFirst, iLast, jFirst, jLast
        real(kind=8) :: dx, xOffset, yOffset, ownedBox(4)

        real(kind=8) :: f(ni, nj, 0:8)
        real(kind=8) :: g(ni, nj, 0:4)
        real(kind=8) :: f_post(0:ni+1, 0:nj+1, 0:8)
        real(kind=8) :: g_post(0:ni+1, 0:nj+1, 0:4)
        real(kind=8) :: rho(ni, nj)
        real(kind=8) :: u(ni, nj)
        real(kind=8) :: v(ni, nj)
        real(kind=8) :: T(ni, nj)
        real(kind=8) :: Fx(ni, nj)
        real(kind=8) :: Fy(ni, nj)
        real(kind=8) :: rhoHistory(ni, nj, 0:nh)
        real(kind=8) :: uHistory(ni, nj, 0:nh)
        real(kind=8) :: vHistory(ni, nj, 0:nh)
        real(kind=8) :: THistory(ni, nj, 0:nh)
        real(kind=8) :: FxHistory(ni, nj, 0:nh)
        real(kind=8) :: FyHistory(ni, nj, 0:nh)
        real(kind=8) :: flowNeqHistory(ni, nj, 0:8, 0:nh)
        real(kind=8) :: thermalNeqHistory(ni, nj, 0:4, 0:nh)
        real(kind=8) :: quadWidthX(ni)
        real(kind=8) :: quadWidthY(nj)
#ifdef steadyFlow
        real(kind=8) :: up(ni, nj)
        real(kind=8) :: vp(ni, nj)
        real(kind=8) :: Tp(ni, nj)
#endif
        integer(kind=4) :: i, j, a
        real(kind=8) :: x, y

        call integration_weights(ni, xOffset, dx, &
            ownedBox(1), ownedBox(2), quadWidthX, iFirst, iLast)
        call integration_weights(nj, yOffset, dx, &
            ownedBox(3), ownedBox(4), quadWidthY, jFirst, jLast)
        u = 0.0d0
        v = 0.0d0
        rho = 1.0d0
        Fx = 0.0d0
        Fy = 0.0d0
        f_post = 0.0d0
        g_post = 0.0d0
        rhoHistory = 0.0d0
        uHistory = 0.0d0
        vHistory = 0.0d0
        THistory = 0.0d0
        FxHistory = 0.0d0
        FyHistory = 0.0d0
        flowNeqHistory = 0.0d0
        thermalNeqHistory = 0.0d0
        do j = 1, nj
            y = (yOffset+(dble(j)-0.5d0)*dx)/lengthUnit
            do i = 1, ni
                x = (xOffset+(dble(i)-0.5d0)*dx)/lengthUnit
#ifdef SideHeatedCell
                T(i, j) = Thot+x/(dble(nx)/lengthUnit)*(Tcold-Thot)
#else
                T(i, j) = Thot+y/(dble(ny)/lengthUnit)*(Tcold-Thot)
                T(i, j) = T(i, j)+1.0d-3*(Thot-Tcold)*sin(2.0d0*pi*x/(dble(nx)/lengthUnit))* &
                    sin(pi*y/(dble(ny)/lengthUnit))
#endif
                do a = 0, 8
                    f(i, j, a) = omega(a)
                enddo
                do a = 0, 4
                    g(i, j, a) = omegaT(a)*T(i, j)
                enddo
            enddo
        enddo
#ifdef steadyFlow
        up = u
        vp = v
        Tp = T
#endif
    end subroutine initial_grid


    !===============================================================================================
    ! 五套场量与历史数组一次进入设备，结束时统一释放。
    subroutine grid_device_data(entering)
        use commondata
        implicit none
        logical, intent(in) :: entering

        if (entering) then
        !$acc enter data copyin(f_coarse, g_coarse, f_post_coarse, g_post_coarse, rho_coarse, u_coarse,  &
        !$acc& v_coarse, T_coarse, Fx_coarse, Fy_coarse, rhoHistory_coarse, uHistory_coarse,  &
        !$acc& vHistory_coarse, THistory_coarse, FxHistory_coarse, FyHistory_coarse, flowNeqHistory_coarse,  &
        !$acc& thermalNeqHistory_coarse)
        if (refineRatio > 1) then
            !$acc enter data copyin(f_left, g_left, f_post_left, g_post_left, rho_left, u_left, v_left, T_left,  &
            !$acc& Fx_left, Fy_left, rhoHistory_left, uHistory_left, vHistory_left, THistory_left,  &
            !$acc& FxHistory_left, FyHistory_left, flowNeqHistory_left, thermalNeqHistory_left)
            !$acc enter data copyin(f_right, g_right, f_post_right, g_post_right, rho_right, u_right, v_right,  &
            !$acc& T_right, Fx_right, Fy_right, rhoHistory_right, uHistory_right, vHistory_right,  &
            !$acc& THistory_right, FxHistory_right, FyHistory_right, flowNeqHistory_right,  &
            !$acc& thermalNeqHistory_right)
            !$acc enter data copyin(f_bottom, g_bottom, f_post_bottom, g_post_bottom, rho_bottom, u_bottom,  &
            !$acc& v_bottom, T_bottom, Fx_bottom, Fy_bottom, rhoHistory_bottom, uHistory_bottom,  &
            !$acc& vHistory_bottom, THistory_bottom, FxHistory_bottom, FyHistory_bottom, flowNeqHistory_bottom,  &
            !$acc& thermalNeqHistory_bottom)
            !$acc enter data copyin(f_top, g_top, f_post_top, g_post_top, rho_top, u_top, v_top, T_top, Fx_top,  &
            !$acc& Fy_top, rhoHistory_top, uHistory_top, vHistory_top, THistory_top, FxHistory_top,  &
            !$acc& FyHistory_top, flowNeqHistory_top, thermalNeqHistory_top)
        endif
        else
        !$acc exit data delete(f_coarse, g_coarse, f_post_coarse, g_post_coarse, rho_coarse, u_coarse,  &
        !$acc& v_coarse, T_coarse, Fx_coarse, Fy_coarse, rhoHistory_coarse, uHistory_coarse,  &
        !$acc& vHistory_coarse, THistory_coarse, FxHistory_coarse, FyHistory_coarse, flowNeqHistory_coarse,  &
        !$acc& thermalNeqHistory_coarse)
        if (refineRatio > 1) then
            !$acc exit data delete(f_left, g_left, f_post_left, g_post_left, rho_left, u_left, v_left, T_left,  &
            !$acc& Fx_left, Fy_left, rhoHistory_left, uHistory_left, vHistory_left, THistory_left,  &
            !$acc& FxHistory_left, FyHistory_left, flowNeqHistory_left, thermalNeqHistory_left)
            !$acc exit data delete(f_right, g_right, f_post_right, g_post_right, rho_right, u_right, v_right,  &
            !$acc& T_right, Fx_right, Fy_right, rhoHistory_right, uHistory_right, vHistory_right,  &
            !$acc& THistory_right, FxHistory_right, FyHistory_right, flowNeqHistory_right,  &
            !$acc& thermalNeqHistory_right)
            !$acc exit data delete(f_bottom, g_bottom, f_post_bottom, g_post_bottom, rho_bottom, u_bottom,  &
            !$acc& v_bottom, T_bottom, Fx_bottom, Fy_bottom, rhoHistory_bottom, uHistory_bottom,  &
            !$acc& vHistory_bottom, THistory_bottom, FxHistory_bottom, FyHistory_bottom, flowNeqHistory_bottom,  &
            !$acc& thermalNeqHistory_bottom)
            !$acc exit data delete(f_top, g_top, f_post_top, g_post_top, rho_top, u_top, v_top, T_top, Fx_top,  &
            !$acc& Fy_top, rhoHistory_top, uHistory_top, vHistory_top, THistory_top, FxHistory_top,  &
            !$acc& FyHistory_top, flowNeqHistory_top, thermalNeqHistory_top)
        endif
        endif
    end subroutine grid_device_data


    !===============================================================================================
    ! 建立设备数据区，并准备首次碰撞之前的粗细接口状态。
    subroutine enter_data_2d_openacc()
        use commondata
        implicit none

        !$acc enter data copyin(ex, ey, omega, omegaT)
        call grid_device_data(.true.)
        call interface_device_data(.true.)
        if (loadInitField == 0) then
        call save_exchange_history(nxCoarse, nyCoarse, historyLastCoarse, min(1, historyLastCoarse), &
            dxCoarse, SnuCoarse, SqCoarse, QkCoarse, QnuCoarse, f_coarse, g_coarse, rho_coarse, u_coarse, &
            v_coarse, T_coarse, Fx_coarse, Fy_coarse, rhoHistory_coarse, uHistory_coarse, vHistory_coarse, &
            THistory_coarse, FxHistory_coarse, FyHistory_coarse, flowNeqHistory_coarse, thermalNeqHistory_coarse)
            if (refineRatio > 1) then
        call save_exchange_history(nxLeft, nyLeft, 0, 0, dxFine, SnuFine, SqFine, QkFine, QnuFine, f_left, &
            g_left, rho_left, u_left, v_left, T_left, Fx_left, Fy_left, rhoHistory_left, uHistory_left, &
            vHistory_left, THistory_left, FxHistory_left, FyHistory_left, flowNeqHistory_left, &
            thermalNeqHistory_left)
        call save_exchange_history(nxRight, nyRight, 0, 0, dxFine, SnuFine, SqFine, QkFine, QnuFine, &
            f_right, g_right, rho_right, u_right, v_right, T_right, Fx_right, Fy_right, rhoHistory_right, &
            uHistory_right, vHistory_right, THistory_right, FxHistory_right, FyHistory_right, &
            flowNeqHistory_right, thermalNeqHistory_right)
        call save_exchange_history(nxBottom, nyBottom, 0, 0, dxFine, SnuFine, SqFine, QkFine, QnuFine, &
            f_bottom, g_bottom, rho_bottom, u_bottom, v_bottom, T_bottom, Fx_bottom, Fy_bottom, &
            rhoHistory_bottom, uHistory_bottom, vHistory_bottom, THistory_bottom, FxHistory_bottom, &
            FyHistory_bottom, flowNeqHistory_bottom, thermalNeqHistory_bottom)
        call save_exchange_history(nxTop, nyTop, 0, 0, dxFine, SnuFine, SqFine, QkFine, QnuFine, f_top, &
            g_top, rho_top, u_top, v_top, T_top, Fx_top, Fy_top, rhoHistory_top, uHistory_top, &
            vHistory_top, THistory_top, FxHistory_top, FyHistory_top, flowNeqHistory_top, thermalNeqHistory_top)
                call fine_to_coarse()
                call coarse_to_fine([0.0d0, 1.0d0, 0.0d0])
        call save_exchange_history(nxCoarse, nyCoarse, historyLastCoarse, 1, dxCoarse, SnuCoarse, SqCoarse, &
            QkCoarse, QnuCoarse, f_coarse, g_coarse, rho_coarse, u_coarse, v_coarse, T_coarse, Fx_coarse, &
            Fy_coarse, rhoHistory_coarse, uHistory_coarse, vHistory_coarse, THistory_coarse, &
            FxHistory_coarse, FyHistory_coarse, flowNeqHistory_coarse, thermalNeqHistory_coarse)
        call save_exchange_history(nxLeft, nyLeft, 0, 0, dxFine, SnuFine, SqFine, QkFine, QnuFine, f_left, &
            g_left, rho_left, u_left, v_left, T_left, Fx_left, Fy_left, rhoHistory_left, uHistory_left, &
            vHistory_left, THistory_left, FxHistory_left, FyHistory_left, flowNeqHistory_left, &
            thermalNeqHistory_left)
        call save_exchange_history(nxRight, nyRight, 0, 0, dxFine, SnuFine, SqFine, QkFine, QnuFine, &
            f_right, g_right, rho_right, u_right, v_right, T_right, Fx_right, Fy_right, rhoHistory_right, &
            uHistory_right, vHistory_right, THistory_right, FxHistory_right, FyHistory_right, &
            flowNeqHistory_right, thermalNeqHistory_right)
        call save_exchange_history(nxBottom, nyBottom, 0, 0, dxFine, SnuFine, SqFine, QkFine, QnuFine, &
            f_bottom, g_bottom, rho_bottom, u_bottom, v_bottom, T_bottom, Fx_bottom, Fy_bottom, &
            rhoHistory_bottom, uHistory_bottom, vHistory_bottom, THistory_bottom, FxHistory_bottom, &
            FyHistory_bottom, flowNeqHistory_bottom, thermalNeqHistory_bottom)
        call save_exchange_history(nxTop, nyTop, 0, 0, dxFine, SnuFine, SqFine, QkFine, QnuFine, f_top, &
            g_top, rho_top, u_top, v_top, T_top, Fx_top, Fy_top, rhoHistory_top, uHistory_top, &
            vHistory_top, THistory_top, FxHistory_top, FyHistory_top, flowNeqHistory_top, thermalNeqHistory_top)
            endif
        endif
    end subroutine enter_data_2d_openacc


    ! 将传入数组的宏观量或完整状态从设备复制到主机。
    subroutine update_host_grid(ni, nj, nh, full, f, g, rho, u, v, T, Fx, Fy, rhoHistory, uHistory, vHistory, &
        THistory, FxHistory, FyHistory, flowNeqHistory, thermalNeqHistory)


        implicit none
        integer(kind=4) :: ni, nj, nh

        real(kind=8) :: f(ni, nj, 0:8)
        real(kind=8) :: g(ni, nj, 0:4)
        real(kind=8) :: rho(ni, nj)
        real(kind=8) :: u(ni, nj)
        real(kind=8) :: v(ni, nj)
        real(kind=8) :: T(ni, nj)
        real(kind=8) :: Fx(ni, nj)
        real(kind=8) :: Fy(ni, nj)
        real(kind=8) :: rhoHistory(ni, nj, 0:nh)
        real(kind=8) :: uHistory(ni, nj, 0:nh)
        real(kind=8) :: vHistory(ni, nj, 0:nh)
        real(kind=8) :: THistory(ni, nj, 0:nh)
        real(kind=8) :: FxHistory(ni, nj, 0:nh)
        real(kind=8) :: FyHistory(ni, nj, 0:nh)
        real(kind=8) :: flowNeqHistory(ni, nj, 0:8, 0:nh)
        real(kind=8) :: thermalNeqHistory(ni, nj, 0:4, 0:nh)
        logical, intent(in) :: full

        !$acc update self(u, v, T, rho) async(1)
        if (full) then
            !$acc update self(f, g, Fx, Fy, rhoHistory, uHistory, vHistory, THistory, FxHistory, FyHistory, &
            !$acc& flowNeqHistory, thermalNeqHistory) async(1)
        endif
    end subroutine update_host_grid


    !===============================================================================================
    ! 显式同步中心、左、右、下、上数组，供积分与文件输出使用。
    subroutine update_host_all(full)
        use commondata
        implicit none
        logical, intent(in) :: full

        call update_host_grid(nxCoarse, nyCoarse, historyLastCoarse, full, f_coarse, g_coarse, rho_coarse, &
            u_coarse, v_coarse, T_coarse, Fx_coarse, Fy_coarse, rhoHistory_coarse, uHistory_coarse, &
            vHistory_coarse, THistory_coarse, FxHistory_coarse, FyHistory_coarse, flowNeqHistory_coarse, &
            thermalNeqHistory_coarse)
        if (refineRatio > 1) then
            call update_host_grid(nxLeft, nyLeft, 0, full, f_left, g_left, rho_left, u_left, v_left, T_left, &
                Fx_left, Fy_left, rhoHistory_left, uHistory_left, vHistory_left, THistory_left, FxHistory_left, &
                FyHistory_left, flowNeqHistory_left, thermalNeqHistory_left)
            call update_host_grid(nxRight, nyRight, 0, full, f_right, g_right, rho_right, u_right, v_right, &
                T_right, Fx_right, Fy_right, rhoHistory_right, uHistory_right, vHistory_right, THistory_right, &
                FxHistory_right, FyHistory_right, flowNeqHistory_right, thermalNeqHistory_right)
            call update_host_grid(nxBottom, nyBottom, 0, full, f_bottom, g_bottom, rho_bottom, u_bottom, &
                v_bottom, T_bottom, Fx_bottom, Fy_bottom, rhoHistory_bottom, uHistory_bottom, vHistory_bottom, &
                THistory_bottom, FxHistory_bottom, FyHistory_bottom, flowNeqHistory_bottom, &
            thermalNeqHistory_bottom)
            call update_host_grid(nxTop, nyTop, 0, full, f_top, g_top, rho_top, u_top, v_top, T_top, Fx_top, &
                Fy_top, rhoHistory_top, uHistory_top, vHistory_top, THistory_top, FxHistory_top, FyHistory_top, &
                flowNeqHistory_top, thermalNeqHistory_top)
        endif
        !$acc wait(1)
    end subroutine update_host_all


    !===============================================================================================
    subroutine exit_data_2d_openacc()
        use commondata
        implicit none

        !$acc wait(1)
        call interface_device_data(.false.)
        call grid_device_data(.false.)
        !$acc exit data delete(ex, ey, omega, omegaT)
    end subroutine exit_data_2d_openacc


    ! 按碰撞、迁移、墙面处理、宏观恢复的次序推进传入的整套数组。
    subroutine advance_grid(ni, nj, dx, sn, sq, qk, qn, gb, wall, f, g, f_post, g_post, rho, u, v, T, Fx, Fy)

#ifdef SideHeatedHa
        use commondata, only: B2sigemarho
#endif
        implicit none
        integer(kind=4) :: ni, nj
        real(kind=8) :: dx, sn, sq, qk, qn, gb
        logical :: wall(4)

        real(kind=8) :: f(ni, nj, 0:8)
        real(kind=8) :: g(ni, nj, 0:4)
        real(kind=8) :: f_post(0:ni+1, 0:nj+1, 0:8)
        real(kind=8) :: g_post(0:ni+1, 0:nj+1, 0:4)
        real(kind=8) :: rho(ni, nj)
        real(kind=8) :: u(ni, nj)
        real(kind=8) :: v(ni, nj)
        real(kind=8) :: T(ni, nj)
        real(kind=8) :: Fx(ni, nj)
        real(kind=8) :: Fy(ni, nj)
        ! 原文件的执行次序保持不变。粗块的 dt 已吸收到松弛率及力增量中。
        ! 入口：人工边界已重建为本时刻的碰撞前状态；缓冲节点与内部节点一样参与碰撞。
        ! 出口：最外两层只作为待重建缓冲，不允许充当 donor；原算法的内部结果保留。
        call collision(ni, nj, f, f_post, rho, u, v, Fx, Fy, T, sn, sq, &
            gb &
#ifdef SideHeatedHa
            , dx*B2sigemarho &
#endif
            )
        call streaming(ni, nj, f, f_post)
        call bounceback(ni, nj, f, f_post, wall(1), wall(2), &
            wall(3), wall(4))
        call macro(ni, nj, f, rho, u, v, Fx, Fy)
        call collisionT(ni, nj, g, g_post, u, v, T, qk, qn)
        call streamingT(ni, nj, g, g_post)
        call bouncebackT(ni, nj, g, g_post, wall(1), wall(2), &
            wall(3), wall(4))
        call macroT(ni, nj, g, T)
    end subroutine advance_grid


    !===============================================================================================
    ! 传入中心数组与粗网格松弛率，推进一个粗时间步。
    subroutine advance_coarse()
        use commondata
        implicit none

        call advance_grid(nxCoarse, nyCoarse, dxCoarse, SnuCoarse, SqCoarse, QkCoarse, QnuCoarse, &
            gBetaCoarse, wallCoarse, f_coarse, g_coarse, f_post_coarse, g_post_coarse, rho_coarse, &
            u_coarse, v_coarse, T_coarse, Fx_coarse, Fy_coarse)
    end subroutine advance_coarse


    !===============================================================================================
    ! 四个细区分阶段同步推进；每次迁移前交换相邻细区的碰撞后分布。
    subroutine advance_fine()
        use commondata
        implicit none

        ! 四个细区先全部碰撞，再交换迁移外圈；避免读取另一细区尚未更新的数据。
        call collision(nxLeft, nyLeft, f_left, f_post_left, rho_left, u_left, v_left, Fx_left, Fy_left, &
            T_left, SnuFine, SqFine, gBeta &
#ifdef SideHeatedHa
            , B2sigemarho &
#endif
        )
        call collision(nxRight, nyRight, f_right, f_post_right, rho_right, u_right, v_right, Fx_right, &
            Fy_right, T_right, SnuFine, SqFine, gBeta &
#ifdef SideHeatedHa
            , B2sigemarho &
#endif
        )
        call collision(nxBottom, nyBottom, f_bottom, f_post_bottom, rho_bottom, u_bottom, v_bottom, &
            Fx_bottom, Fy_bottom, T_bottom, SnuFine, SqFine, gBeta &
#ifdef SideHeatedHa
            , B2sigemarho &
#endif
        )
        call collision(nxTop, nyTop, f_top, f_post_top, rho_top, u_top, v_top, Fx_top, Fy_top, T_top, &
            SnuFine, SqFine, gBeta &
#ifdef SideHeatedHa
            , B2sigemarho &
#endif
        )
        call exchange_fine_flow_halos()
        call streaming(nxLeft, nyLeft, f_left, f_post_left)
        call bounceback(nxLeft, nyLeft, f_left, f_post_left, wallLeft(1), wallLeft(2), wallLeft(3), wallLeft(4))
        call macro(nxLeft, nyLeft, f_left, rho_left, u_left, v_left, Fx_left, Fy_left)
        call streaming(nxRight, nyRight, f_right, f_post_right)
        call bounceback(nxRight, nyRight, f_right, f_post_right, wallRight(1), wallRight(2), wallRight(3), &
            wallRight(4))
        call macro(nxRight, nyRight, f_right, rho_right, u_right, v_right, Fx_right, Fy_right)
        call streaming(nxBottom, nyBottom, f_bottom, f_post_bottom)
        call bounceback(nxBottom, nyBottom, f_bottom, f_post_bottom, wallBottom(1), wallBottom(2), &
            wallBottom(3), wallBottom(4))
        call macro(nxBottom, nyBottom, f_bottom, rho_bottom, u_bottom, v_bottom, Fx_bottom, Fy_bottom)
        call streaming(nxTop, nyTop, f_top, f_post_top)
        call bounceback(nxTop, nyTop, f_top, f_post_top, wallTop(1), wallTop(2), wallTop(3), wallTop(4))
        call macro(nxTop, nyTop, f_top, rho_top, u_top, v_top, Fx_top, Fy_top)
        call collisionT(nxLeft, nyLeft, g_left, g_post_left, u_left, v_left, T_left, QkFine, QnuFine)
        call collisionT(nxRight, nyRight, g_right, g_post_right, u_right, v_right, T_right, QkFine, QnuFine)
        call collisionT(nxBottom, nyBottom, g_bottom, g_post_bottom, u_bottom, v_bottom, T_bottom, QkFine, QnuFine)
        call collisionT(nxTop, nyTop, g_top, g_post_top, u_top, v_top, T_top, QkFine, QnuFine)
        call exchange_fine_thermal_halos()
        call streamingT(nxLeft, nyLeft, g_left, g_post_left)
        call bouncebackT(nxLeft, nyLeft, g_left, g_post_left, wallLeft(1), wallLeft(2), wallLeft(3), wallLeft(4))
        call macroT(nxLeft, nyLeft, g_left, T_left)
        call streamingT(nxRight, nyRight, g_right, g_post_right)
        call bouncebackT(nxRight, nyRight, g_right, g_post_right, wallRight(1), wallRight(2), wallRight(3), &
            wallRight(4))
        call macroT(nxRight, nyRight, g_right, T_right)
        call streamingT(nxBottom, nyBottom, g_bottom, g_post_bottom)
        call bouncebackT(nxBottom, nyBottom, g_bottom, g_post_bottom, wallBottom(1), wallBottom(2), &
            wallBottom(3), wallBottom(4))
        call macroT(nxBottom, nyBottom, g_bottom, T_bottom)
        call streamingT(nxTop, nyTop, g_top, g_post_top)
        call bouncebackT(nxTop, nyTop, g_top, g_post_top, wallTop(1), wallTop(2), wallTop(3), wallTop(4))
        call macroT(nxTop, nyTop, g_top, T_top)
    end subroutine advance_fine


    !===============================================================================================
    ! 在左右细区与上下细区之间双向复制 f_post 的迁移外圈。
    subroutine exchange_fine_flow_halos()
        use commondata
        implicit none

        call copy_fine_halo(9, nxLeft, nyLeft, nint(xOffsetLeft), nint(yOffsetLeft), f_post_left, nxBottom, &
            nyBottom, nint(xOffsetBottom), nint(yOffsetBottom), f_post_bottom)
        call copy_fine_halo(9, nxBottom, nyBottom, nint(xOffsetBottom), nint(yOffsetBottom), f_post_bottom, &
            nxLeft, nyLeft, nint(xOffsetLeft), nint(yOffsetLeft), f_post_left)
        call copy_fine_halo(9, nxLeft, nyLeft, nint(xOffsetLeft), nint(yOffsetLeft), f_post_left, nxTop, &
            nyTop, nint(xOffsetTop), nint(yOffsetTop), f_post_top)
        call copy_fine_halo(9, nxTop, nyTop, nint(xOffsetTop), nint(yOffsetTop), f_post_top, nxLeft, &
            nyLeft, nint(xOffsetLeft), nint(yOffsetLeft), f_post_left)
        call copy_fine_halo(9, nxRight, nyRight, nint(xOffsetRight), nint(yOffsetRight), f_post_right, &
            nxBottom, nyBottom, nint(xOffsetBottom), nint(yOffsetBottom), f_post_bottom)
        call copy_fine_halo(9, nxBottom, nyBottom, nint(xOffsetBottom), nint(yOffsetBottom), f_post_bottom, &
            nxRight, nyRight, nint(xOffsetRight), nint(yOffsetRight), f_post_right)
        call copy_fine_halo(9, nxRight, nyRight, nint(xOffsetRight), nint(yOffsetRight), f_post_right, &
            nxTop, nyTop, nint(xOffsetTop), nint(yOffsetTop), f_post_top)
        call copy_fine_halo(9, nxTop, nyTop, nint(xOffsetTop), nint(yOffsetTop), f_post_top, nxRight, &
            nyRight, nint(xOffsetRight), nint(yOffsetRight), f_post_right)
    end subroutine exchange_fine_flow_halos


    !===============================================================================================
    ! 在左右细区与上下细区之间双向复制 g_post 的迁移外圈。
    subroutine exchange_fine_thermal_halos()
        use commondata
        implicit none

        call copy_fine_halo(5, nxLeft, nyLeft, nint(xOffsetLeft), nint(yOffsetLeft), g_post_left, nxBottom, &
            nyBottom, nint(xOffsetBottom), nint(yOffsetBottom), g_post_bottom)
        call copy_fine_halo(5, nxBottom, nyBottom, nint(xOffsetBottom), nint(yOffsetBottom), g_post_bottom, &
            nxLeft, nyLeft, nint(xOffsetLeft), nint(yOffsetLeft), g_post_left)
        call copy_fine_halo(5, nxLeft, nyLeft, nint(xOffsetLeft), nint(yOffsetLeft), g_post_left, nxTop, &
            nyTop, nint(xOffsetTop), nint(yOffsetTop), g_post_top)
        call copy_fine_halo(5, nxTop, nyTop, nint(xOffsetTop), nint(yOffsetTop), g_post_top, nxLeft, &
            nyLeft, nint(xOffsetLeft), nint(yOffsetLeft), g_post_left)
        call copy_fine_halo(5, nxRight, nyRight, nint(xOffsetRight), nint(yOffsetRight), g_post_right, &
            nxBottom, nyBottom, nint(xOffsetBottom), nint(yOffsetBottom), g_post_bottom)
        call copy_fine_halo(5, nxBottom, nyBottom, nint(xOffsetBottom), nint(yOffsetBottom), g_post_bottom, &
            nxRight, nyRight, nint(xOffsetRight), nint(yOffsetRight), g_post_right)
        call copy_fine_halo(5, nxRight, nyRight, nint(xOffsetRight), nint(yOffsetRight), g_post_right, &
            nxTop, nyTop, nint(xOffsetTop), nint(yOffsetTop), g_post_top)
        call copy_fine_halo(5, nxTop, nyTop, nint(xOffsetTop), nint(yOffsetTop), g_post_top, nxRight, &
            nyRight, nint(xOffsetRight), nint(yOffsetRight), g_post_right)
    end subroutine exchange_fine_thermal_halos


    ! 同级细区只复制真实相邻节点的碰撞后值，含 D2Q9 对角方向；不做插值或尺度变换。
    subroutine copy_fine_halo(nq, ni, nj, ix, iy, receiver, mi, mj, jx, jy, donor)
        implicit none
        integer, intent(in) :: nq, ni, nj, ix, iy, mi, mj, jx, jy
        real(8) :: receiver(0:ni+1,0:nj+1,0:nq-1), donor(0:mi+1,0:mj+1,0:nq-1)
        integer :: i,j,a,il,ih,jl,jh
        ! 全局细节点编号 = 偏移量 + 本地编号。只取 receiver 外圈与 donor 内部的交集。
        il=max(0,jx+1-ix); ih=min(ni+1,jx+mi-ix)
        jl=max(0,jy+1-iy); jh=min(nj+1,jy+mj-iy)
        !$acc parallel loop collapse(3) present(receiver,donor) async(1)
        do a=0,nq-1
            do j=jl,jh
                do i=il,ih
                    receiver(i,j,a)=donor(ix+i-jx,iy+j-jy,a)
                enddo
            enddo
        enddo
    end subroutine copy_fine_halo


    !===============================================================================================
    ! 预测粗网格，执行 refineRatio 个细步，同步接口并滚动时间历史。
    subroutine advance_multiblock()
        use commondata
        implicit none
        integer(kind=4) :: k
        real(kind=8) :: theta, wt(0:2)

        if (refineRatio == 1) then
            call advance_coarse()
            itc = itc+1
            return
        endif
        ! 粗网格先预测一个粗步；History 的 0/1/2 层保存过去、当前、预测时刻。
        call advance_coarse()
        call save_exchange_history(nxCoarse, nyCoarse, historyLastCoarse, 2, dxCoarse, SnuCoarse, SqCoarse, &
            QkCoarse, QnuCoarse, f_coarse, g_coarse, rho_coarse, u_coarse, v_coarse, T_coarse, Fx_coarse, &
            Fy_coarse, rhoHistory_coarse, uHistory_coarse, vHistory_coarse, THistory_coarse, &
            FxHistory_coarse, FyHistory_coarse, flowNeqHistory_coarse, thermalNeqHistory_coarse)
        do k = 1, refineRatio
            theta = dble(k-1)/dble(refineRatio)
            call coarse_time_weights(theta, wt)
            call coarse_to_fine(wt)
            call advance_fine()
        call save_exchange_history(nxLeft, nyLeft, 0, 0, dxFine, SnuFine, SqFine, QkFine, QnuFine, f_left, &
            g_left, rho_left, u_left, v_left, T_left, Fx_left, Fy_left, rhoHistory_left, uHistory_left, &
            vHistory_left, THistory_left, FxHistory_left, FyHistory_left, flowNeqHistory_left, &
            thermalNeqHistory_left)
        call save_exchange_history(nxRight, nyRight, 0, 0, dxFine, SnuFine, SqFine, QkFine, QnuFine, &
            f_right, g_right, rho_right, u_right, v_right, T_right, Fx_right, Fy_right, rhoHistory_right, &
            uHistory_right, vHistory_right, THistory_right, FxHistory_right, FyHistory_right, &
            flowNeqHistory_right, thermalNeqHistory_right)
        call save_exchange_history(nxBottom, nyBottom, 0, 0, dxFine, SnuFine, SqFine, QkFine, QnuFine, &
            f_bottom, g_bottom, rho_bottom, u_bottom, v_bottom, T_bottom, Fx_bottom, Fy_bottom, &
            rhoHistory_bottom, uHistory_bottom, vHistory_bottom, THistory_bottom, FxHistory_bottom, &
            FyHistory_bottom, flowNeqHistory_bottom, thermalNeqHistory_bottom)
        call save_exchange_history(nxTop, nyTop, 0, 0, dxFine, SnuFine, SqFine, QkFine, QnuFine, f_top, &
            g_top, rho_top, u_top, v_top, T_top, Fx_top, Fy_top, rhoHistory_top, uHistory_top, &
            vHistory_top, THistory_top, FxHistory_top, FyHistory_top, flowNeqHistory_top, thermalNeqHistory_top)
        enddo
        ! 同步时先细到粗，再用修正后的粗状态补细接口，最后滚动粗时间历史。
        call fine_to_coarse()
        call save_exchange_history(nxCoarse, nyCoarse, historyLastCoarse, 2, dxCoarse, SnuCoarse, SqCoarse, &
            QkCoarse, QnuCoarse, f_coarse, g_coarse, rho_coarse, u_coarse, v_coarse, T_coarse, Fx_coarse, &
            Fy_coarse, rhoHistory_coarse, uHistory_coarse, vHistory_coarse, THistory_coarse, &
            FxHistory_coarse, FyHistory_coarse, flowNeqHistory_coarse, thermalNeqHistory_coarse)
        call coarse_to_fine([0.0d0, 0.0d0, 1.0d0])
        call rotate_coarse_history(nxCoarse, nyCoarse, rhoHistory_coarse, uHistory_coarse, vHistory_coarse, &
            THistory_coarse, FxHistory_coarse, FyHistory_coarse, flowNeqHistory_coarse, thermalNeqHistory_coarse)
        itc = itc+refineRatio
    end subroutine advance_multiblock


    subroutine coarse_time_weights(theta, wt)

        use commondata, only: itc
        implicit none

        real(kind=8), intent(in) :: theta
        real(kind=8), intent(out) :: wt(0:2)
        ! theta 是当前细子步起点相对粗步起点的时间，不能使用子步终点时间。
        if (itc == 0) then
            wt = [0.0d0, 1.0d0-theta, theta]    ! 缺少负时间层时线性启动
        else
            wt = [0.5d0*theta*(theta-1.0d0), 1.0d0-theta**2, 0.5d0*theta*(theta+1.0d0)]
        endif
    end subroutine coarse_time_weights


    subroutine rotate_coarse_history(ni, nj, rhoHistory, uHistory, vHistory, THistory, FxHistory, FyHistory, &
        flowNeqHistory, thermalNeqHistory)
        implicit none
        integer(kind=4), intent(in) :: ni, nj
        real(kind=8), intent(inout) :: rhoHistory(ni, nj, 0:2)
        real(kind=8), intent(inout) :: uHistory(ni, nj, 0:2)
        real(kind=8), intent(inout) :: vHistory(ni, nj, 0:2)
        real(kind=8), intent(inout) :: THistory(ni, nj, 0:2)
        real(kind=8), intent(inout) :: FxHistory(ni, nj, 0:2)
        real(kind=8), intent(inout) :: FyHistory(ni, nj, 0:2)
        real(kind=8), intent(inout) :: flowNeqHistory(ni, nj, 0:8, 0:2)
        real(kind=8), intent(inout) :: thermalNeqHistory(ni, nj, 0:4, 0:2)
        integer(kind=4) :: i, j, a
        ! 预测层转为当前层，原当前层转为过去层；下一个粗步会重新计算第 2 层。
        !$acc parallel loop collapse(2) present(rhoHistory, uHistory, vHistory, THistory, FxHistory, FyHistory, &
        !$acc& flowNeqHistory, thermalNeqHistory) async(1) private(a)
        do j = 1, nj
            do i = 1, ni
                rhoHistory(i,j,0) = rhoHistory(i,j,1)
                rhoHistory(i,j,1) = rhoHistory(i,j,2)
                uHistory(i,j,0) = uHistory(i,j,1)
                uHistory(i,j,1) = uHistory(i,j,2)
                vHistory(i,j,0) = vHistory(i,j,1)
                vHistory(i,j,1) = vHistory(i,j,2)
                THistory(i,j,0) = THistory(i,j,1)
                THistory(i,j,1) = THistory(i,j,2)
                FxHistory(i,j,0) = FxHistory(i,j,1)
                FxHistory(i,j,1) = FxHistory(i,j,2)
                FyHistory(i,j,0) = FyHistory(i,j,1)
                FyHistory(i,j,1) = FyHistory(i,j,2)
                do a = 0, 8
                    flowNeqHistory(i,j,a,0) = flowNeqHistory(i,j,a,1)
                    flowNeqHistory(i,j,a,1) = flowNeqHistory(i,j,a,2)
                enddo
                do a = 0, 4
                    thermalNeqHistory(i,j,a,0) = thermalNeqHistory(i,j,a,1)
                    thermalNeqHistory(i,j,a,1) = thermalNeqHistory(i,j,a,2)
                enddo
            enddo
        enddo
    end subroutine rotate_coarse_history


    subroutine lagrange_weights(q, w)

        implicit none

        real(kind=8), intent(in) :: q
        real(kind=8), intent(out) :: w(4)
        integer(kind=4) :: a, k

        w = 1.0d0
        do a = 0, 3
            do k = 0, 3
                if (a /= k) w(a+1) = w(a+1)*(q-dble(k))/dble(a-k)
            enddo
        enddo
    end subroutine lagrange_weights


    subroutine section_weights(q, n, first, w, dw)

        implicit none

        real(kind=8), intent(in) :: q
        integer(kind=4), intent(in) :: n
        integer(kind=4), intent(out) :: first
        real(kind=8), intent(out) :: w(4), dw(4)
        real(kind=8) :: z, term
        integer(kind=4) :: a, k, l

        first = max(1, min(floor(q)-1, n-3))
        z = q-dble(first)
        call lagrange_weights(z, w)
        dw = 0.0d0
        do a = 0, 3
            do k = 0, 3
                if (k == a) cycle
                term = 1.0d0/dble(a-k)
                do l = 0, 3
                    if (l /= a .and. l /= k) term = term*(z-dble(l))/dble(a-l)
                enddo
                dw(a+1) = dw(a+1)+term
            enddo
        enddo
    end subroutine section_weights


    subroutine equilibrium_moments(rh, ux, uy, temp, meq, neq)

        use commondata, only: thermalA
        implicit none

        !$acc routine seq
        real(kind=8), intent(in) :: rh, ux, uy, temp
        real(kind=8), intent(out) :: meq(0:8), neq(0:4)

        meq(0) = rh
        meq(1) = rh*(-2.0d0+3.0d0*(ux*ux+uy*uy))
        meq(2) = rh*(1.0d0-3.0d0*(ux*ux+uy*uy))
        meq(3) = rh*ux
        meq(4) = -rh*ux
        meq(5) = rh*uy
        meq(6) = -rh*uy
        meq(7) = rh*(ux*ux-uy*uy)
        meq(8) = rh*ux*uy
        neq = [temp, temp*ux, temp*uy, thermalA*temp, 0.0d0]
    end subroutine equilibrium_moments


    subroutine force_moments(ux, uy, fx, fy, fm)

        implicit none

        !$acc routine seq
        real(kind=8), intent(in) :: ux, uy, fx, fy
        real(kind=8), intent(out) :: fm(0:8)

        fm = [0.0d0, 6.0d0*(ux*fx+uy*fy), -6.0d0*(ux*fx+uy*fy), fx, -fx, fy, -fy, &
            2.0d0*(ux*fx-uy*fy), ux*fy+uy*fx]
    end subroutine force_moments


    ! 保存宏观量、Fx/dx、Fy/dx 和缩放非平衡矩；dt=dx。
    subroutine save_exchange_history(ni, nj, nh, slot, dx, sn, sq, qk, qn, f, g, rho, u, v, T, Fx, Fy, &
        rhoHistory, uHistory, vHistory, THistory, FxHistory, FyHistory, flowNeqHistory, thermalNeqHistory)

        implicit none

        integer(kind=4), intent(in) :: ni, nj, nh, slot
        real(kind=8), intent(in) :: dx, sn, sq, qk, qn, f(ni, nj, 0:8), g(ni, nj, 0:4)
        real(kind=8), intent(in) :: rho(ni, nj), u(ni, nj), v(ni, nj), T(ni, nj), Fx(ni, nj), Fy(ni, nj)
        real(kind=8), intent(inout) :: rhoHistory(ni, nj, 0:nh)
        real(kind=8), intent(inout) :: uHistory(ni, nj, 0:nh)
        real(kind=8), intent(inout) :: vHistory(ni, nj, 0:nh)
        real(kind=8), intent(inout) :: THistory(ni, nj, 0:nh)
        real(kind=8), intent(inout) :: FxHistory(ni, nj, 0:nh)
        real(kind=8), intent(inout) :: FyHistory(ni, nj, 0:nh)
        real(kind=8), intent(inout) :: flowNeqHistory(ni, nj, 0:8, 0:nh)
        real(kind=8), intent(inout) :: thermalNeqHistory(ni, nj, 0:4, 0:nh)
        integer(kind=4) :: i, j, a
        real(kind=8) :: m(0:8), meq(0:8), fm(0:8), n(0:4), neq(0:4), s(0:8), q(0:4), fv(0:8), gv(0:4)
        !$acc parallel loop collapse(2) present(f, g, rho, u, v, T, Fx, Fy, rhoHistory, uHistory, vHistory, &
        !$acc& THistory, FxHistory, FyHistory, flowNeqHistory, thermalNeqHistory) async(1) &
        !$acc& private(a, m, meq, fm, n, neq, s, q, fv, gv)
        do j = 1, nj
            do i = 1, ni
                do a = 0, 8
                    fv(a) = f(i, j, a)
                enddo
                do a = 0, 4
                    gv(a) = g(i, j, a)
                enddo
                call flow_moments(fv, m)
                call thermal_moments(gv, n)
                call equilibrium_moments(rho(i, j), u(i, j), v(i, j), T(i, j), meq, neq)
                call force_moments(u(i, j), v(i, j), Fx(i, j), Fy(i, j), fm)
                s = [0.0d0, sn, sn, 0.0d0, sq, 0.0d0, sq, sn, sn]
                q = [0.0d0, qk, qk, qn, qn]
                rhoHistory(i, j, slot) = rho(i, j)
                uHistory(i, j, slot) = u(i, j)
                vHistory(i, j, slot) = v(i, j)
                THistory(i, j, slot) = T(i, j)
                FxHistory(i, j, slot) = Fx(i, j)/dx
                FyHistory(i, j, slot) = Fy(i, j)/dx
                do a = 0, 8
                    flowNeqHistory(i, j, a, slot) = s(a)/dx*(m(a)-meq(a)+0.5d0*fm(a))
                enddo
                do a = 0, 4
                    thermalNeqHistory(i, j, a, slot) = q(a)/dx*(n(a)-neq(a))
                enddo
            enddo
        enddo
    end subroutine save_exchange_history


    subroutine interpolate_exchange_history(ni, nj, nh, rhoHistory, uHistory, vHistory, THistory, FxHistory, &
        FyHistory, flowNeqHistory, thermalNeqHistory, count, si, sj, wx, wy, coincident, wt, rhoReceive, &
        uReceive, vReceive, TReceive, FxReceive, FyReceive, flowNeqReceive, thermalNeqReceive)
        implicit none
        integer(kind=4), intent(in) :: ni, nj, nh, count, si(count), sj(count)
        real(kind=8), intent(in) :: wx(4,count), wy(4,count), wt(0:2)
        logical, intent(in) :: coincident(count)
        real(kind=8), intent(in) :: rhoHistory(ni, nj, 0:nh)
        real(kind=8), intent(in) :: uHistory(ni, nj, 0:nh)
        real(kind=8), intent(in) :: vHistory(ni, nj, 0:nh)
        real(kind=8), intent(in) :: THistory(ni, nj, 0:nh)
        real(kind=8), intent(in) :: FxHistory(ni, nj, 0:nh)
        real(kind=8), intent(in) :: FyHistory(ni, nj, 0:nh)
        real(kind=8), intent(in) :: flowNeqHistory(ni, nj, 0:8, 0:nh)
        real(kind=8), intent(in) :: thermalNeqHistory(ni, nj, 0:4, 0:nh)
        real(kind=8), intent(out) :: rhoReceive(count)
        real(kind=8), intent(out) :: uReceive(count)
        real(kind=8), intent(out) :: vReceive(count)
        real(kind=8), intent(out) :: TReceive(count)
        real(kind=8), intent(out) :: FxReceive(count)
        real(kind=8), intent(out) :: FyReceive(count)
        real(kind=8), intent(out) :: flowNeqReceive(0:8, count)
        real(kind=8), intent(out) :: thermalNeqReceive(0:4, count)
        ! 每个宏观量分别插值；矩数组的额外下标只表示物理矩编号。
        call interpolate_scalar_history(ni, nj, nh, rhoHistory, count, si, sj, wx, wy, coincident, wt, &
            rhoReceive)
        call interpolate_scalar_history(ni, nj, nh, uHistory, count, si, sj, wx, wy, coincident, wt, uReceive)
        call interpolate_scalar_history(ni, nj, nh, vHistory, count, si, sj, wx, wy, coincident, wt, vReceive)
        call interpolate_scalar_history(ni, nj, nh, THistory, count, si, sj, wx, wy, coincident, wt, TReceive)
        call interpolate_scalar_history(ni, nj, nh, FxHistory, count, si, sj, wx, wy, coincident, wt, FxReceive)
        call interpolate_scalar_history(ni, nj, nh, FyHistory, count, si, sj, wx, wy, coincident, wt, FyReceive)
        call interpolate_moment_history(ni, nj, nh, 9, flowNeqHistory, count, si, sj, wx, wy, coincident, wt, &
            flowNeqReceive)
        call interpolate_moment_history(ni, nj, nh, 5, thermalNeqHistory, count, si, sj, wx, wy, coincident, wt, &
            thermalNeqReceive)
    end subroutine interpolate_exchange_history


    subroutine interpolate_scalar_history(ni, nj, nh, fieldHistory, count, si, sj, wx, wy, coincident, wt, &
            fieldReceive)

        implicit none

        integer(kind=4), intent(in) :: ni, nj, nh, count, si(count), sj(count)
        real(kind=8), intent(in) :: fieldHistory(ni, nj, 0:nh), wx(4, count), wy(4, count), wt(0:2)
        logical, intent(in) :: coincident(count)
        real(kind=8), intent(out) :: fieldReceive(count)
        integer(kind=4) :: c, ix, iy, k
        real(kind=8) :: value, wk, wt0, wt1, wt2
        ! NVHPC 24.3/P100 上数组形参 firstprivate(wt) 会在设备读取 wt(k) 时非法访问。
        ! 在主机端取出三个时间权重，按标量值捕获；保持原三时间层插值及 async(1) 次序。
        wt0 = wt(0)
        wt1 = wt(1)
        wt2 = wt(2)
        !$acc parallel loop present(fieldHistory, si, sj, wx, wy, coincident, &
        !$acc& fieldReceive) firstprivate(wt0, wt1, wt2) async(1) private(ix, iy, k, value, wk)
        do c = 1, count
            value = 0.0d0
            do k = 0, nh
                wk = 1.0d0
                if (nh == 2) then
                    select case (k)
                    case (0)
                        wk = wt0
                    case (1)
                        wk = wt1
                    case (2)
                        wk = wt2
                    end select
                endif
                if (coincident(c)) then
                    value = value+wk*fieldHistory(si(c), sj(c), k)
                else
                    do iy = 1, 4
                        do ix = 1, 4
                            value = value+wk*wx(ix, c)*wy(iy, c)*fieldHistory(si(c)+ix-1, sj(c)+iy-1, k)
                        enddo
                    enddo
                endif
            enddo
            fieldReceive(c) = value
        enddo
    end subroutine interpolate_scalar_history


    subroutine interpolate_moment_history(ni, nj, nh, momentCount, fieldHistory, count, si, sj, wx, wy, &
        coincident, wt, fieldReceive)

        implicit none

        integer(kind=4), intent(in) :: ni, nj, nh, momentCount, count, si(count), sj(count)
        real(kind=8), intent(in) :: fieldHistory(ni, nj, 0:momentCount-1, 0:nh), wx(4, count), wy(4, &
            count), wt(0:2)
        logical, intent(in) :: coincident(count)
        real(kind=8), intent(out) :: fieldReceive(0:momentCount-1, count)
        integer(kind=4) :: c, a, ix, iy, k
        real(kind=8) :: value, wk, wt0, wt1, wt2
        ! NVHPC 24.3/P100 上数组形参 firstprivate(wt) 会在设备读取 wt(k) 时非法访问。
        ! 在主机端取出三个时间权重，按标量值捕获；保持原三时间层插值及 async(1) 次序。
        wt0 = wt(0)
        wt1 = wt(1)
        wt2 = wt(2)
        !$acc parallel loop collapse(2) present(fieldHistory, si, sj, wx, wy, coincident, &
        !$acc& fieldReceive) firstprivate(wt0, wt1, wt2) async(1) private(ix, iy, k, value, wk)
        do c = 1, count
            do a = 0, momentCount-1
                value = 0.0d0
                do k = 0, nh
                    wk = 1.0d0
                    if (nh == 2) then
                        select case (k)
                        case (0)
                            wk = wt0
                        case (1)
                            wk = wt1
                        case (2)
                            wk = wt2
                        end select
                    endif
                    if (coincident(c)) then
                        value = value+wk*fieldHistory(si(c), sj(c), a, k)
                    else
                        do iy = 1, 4
                            do ix = 1, 4
                                value = value+wk*wx(ix, c)*wy(iy, c)*fieldHistory(si(c)+ix-1, sj(c)+iy-1, a, k)
                            enddo
                        enddo
                    endif
                enddo
                fieldReceive(a, c) = value
            enddo
        enddo
    end subroutine interpolate_moment_history


    ! 按接收网格 dx 和松弛率还原力与矩，再重建接收节点 f/g。
    subroutine apply_interface_data(ni, nj, dx, sn, sq, qk, qn, f, g, rho, u, v, T, Fx, Fy, count, ti, &
        tj, rhoReceive, uReceive, vReceive, TReceive, FxReceive, FyReceive, flowNeqReceive, thermalNeqReceive)

        implicit none

        integer(kind=4), intent(in) :: ni, nj, count, ti(count), tj(count)
        real(kind=8), intent(in) :: dx, sn, sq, qk, qn
real(kind=8), intent(in) :: rhoReceive(count)
        real(kind=8), intent(in) :: uReceive(count)
        real(kind=8), intent(in) :: vReceive(count)
        real(kind=8), intent(in) :: TReceive(count)
        real(kind=8), intent(in) :: FxReceive(count)
        real(kind=8), intent(in) :: FyReceive(count)
        real(kind=8), intent(in) :: flowNeqReceive(0:8, count)
        real(kind=8), intent(in) :: thermalNeqReceive(0:4, count)
        real(kind=8), intent(inout) :: f(ni, nj, 0:8), g(ni, nj, 0:4), rho(ni, nj), u(ni, nj), v(ni, nj), T(ni, nj)
        real(kind=8), intent(inout) :: Fx(ni, nj), Fy(ni, nj)
        integer(kind=4) :: c, i, j, a
        real(kind=8) :: m(0:8), meq(0:8), fm(0:8), n(0:4), neq(0:4), s(0:8), q(0:4), fv(0:8), gv(0:4)
        !$acc parallel loop present(f, g, rho, u, v, T, Fx, Fy, ti, tj, rhoReceive, uReceive, vReceive, &
        !$acc& TReceive, FxReceive, FyReceive, flowNeqReceive, thermalNeqReceive) async(1) &
        !$acc& private(i, j, a, m, meq, fm, n, neq, s, q, fv, gv)
        do c = 1, count
            i = ti(c)
            j = tj(c)
            rho(i, j) = rhoReceive(c)
            u(i, j) = uReceive(c)
            v(i, j) = vReceive(c)
            T(i, j) = TReceive(c)
            Fx(i, j) = dx*FxReceive(c)
            Fy(i, j) = dx*FyReceive(c)
            call equilibrium_moments(rho(i, j), u(i, j), v(i, j), T(i, j), meq, neq)
            call force_moments(u(i, j), v(i, j), Fx(i, j), Fy(i, j), fm)
            s = [0.0d0, sn, sn, 0.0d0, sq, 0.0d0, sq, sn, sn]
            q = [0.0d0, qk, qk, qn, qn]
            m = meq
            n = neq
            do a = 0, 8
                if (s(a) > 0.0d0) m(a) = meq(a)+dx/s(a)*flowNeqReceive(a,c)-0.5d0*fm(a)
            enddo
            ! Eq. (21)：守恒矩直接用宏观量和半步力重建，不除以零松弛率。
            m(0) = rho(i, j)
            m(3) = rho(i, j)*u(i, j)-0.5d0*Fx(i, j)
            m(5) = rho(i, j)*v(i, j)-0.5d0*Fy(i, j)
            do a = 1, 4
                n(a) = neq(a)+dx/q(a)*thermalNeqReceive(a,c)
            enddo
            n(0) = T(i, j)
            call flow_populations(m, fv)
            call thermal_populations(n, gv)
            do a = 0, 8
                f(i, j, a) = fv(a)
            enddo
            do a = 0, 4
                g(i, j, a) = gv(a)
            enddo
        enddo
    end subroutine apply_interface_data


    !===============================================================================================
    ! 分别建立 coarseToLeft/Right/Bottom/Top 及反向接口的节点和权重。
    subroutine build_interfaces()
        use commondata
        implicit none
        integer :: n

        if (refineRatio == 1) return
        call count_fine_interface(nxLeft, nyLeft, xOffsetLeft, yOffsetLeft, coarseToLeftCount)
        allocate(coarseToLeftTi(coarseToLeftCount), coarseToLeftTj(coarseToLeftCount), &
            coarseToLeftSi(coarseToLeftCount), coarseToLeftSj(coarseToLeftCount))
        allocate(coarseToLeftSame(coarseToLeftCount))
        allocate(coarseToLeftWx(4, coarseToLeftCount), coarseToLeftWy(4, coarseToLeftCount))
        call fill_fine_interface(nxLeft, nyLeft, xOffsetLeft, yOffsetLeft, coarseToLeftCount, &
            coarseToLeftTi, coarseToLeftTj, coarseToLeftSi, coarseToLeftSj, coarseToLeftWx, coarseToLeftWy, &
            coarseToLeftSame)
        call count_coarse_interface(nxLeft, nyLeft, xOffsetLeft, yOffsetLeft, leftToCoarseCount)
        allocate(leftToCoarseTi(leftToCoarseCount), leftToCoarseTj(leftToCoarseCount), &
            leftToCoarseSi(leftToCoarseCount), leftToCoarseSj(leftToCoarseCount))
        allocate(leftToCoarseSame(leftToCoarseCount))
        allocate(leftToCoarseWx(4, leftToCoarseCount), leftToCoarseWy(4, leftToCoarseCount))
        call fill_coarse_interface(nxLeft, nyLeft, xOffsetLeft, yOffsetLeft, leftToCoarseCount, &
            leftToCoarseTi, leftToCoarseTj, leftToCoarseSi, leftToCoarseSj, leftToCoarseWx, leftToCoarseWy, &
            leftToCoarseSame)
        call count_fine_interface(nxRight, nyRight, xOffsetRight, yOffsetRight, coarseToRightCount)
        allocate(coarseToRightTi(coarseToRightCount), coarseToRightTj(coarseToRightCount), &
            coarseToRightSi(coarseToRightCount), coarseToRightSj(coarseToRightCount))
        allocate(coarseToRightSame(coarseToRightCount))
        allocate(coarseToRightWx(4, coarseToRightCount), coarseToRightWy(4, coarseToRightCount))
        call fill_fine_interface(nxRight, nyRight, xOffsetRight, yOffsetRight, coarseToRightCount, &
            coarseToRightTi, coarseToRightTj, coarseToRightSi, coarseToRightSj, coarseToRightWx, &
            coarseToRightWy, coarseToRightSame)
        call count_coarse_interface(nxRight, nyRight, xOffsetRight, yOffsetRight, rightToCoarseCount)
        allocate(rightToCoarseTi(rightToCoarseCount), rightToCoarseTj(rightToCoarseCount), &
            rightToCoarseSi(rightToCoarseCount), rightToCoarseSj(rightToCoarseCount))
        allocate(rightToCoarseSame(rightToCoarseCount))
        allocate(rightToCoarseWx(4, rightToCoarseCount), rightToCoarseWy(4, rightToCoarseCount))
        call fill_coarse_interface(nxRight, nyRight, xOffsetRight, yOffsetRight, rightToCoarseCount, &
            rightToCoarseTi, rightToCoarseTj, rightToCoarseSi, rightToCoarseSj, rightToCoarseWx, &
            rightToCoarseWy, rightToCoarseSame)
        call count_fine_interface(nxBottom, nyBottom, xOffsetBottom, yOffsetBottom, coarseToBottomCount)
        allocate(coarseToBottomTi(coarseToBottomCount), coarseToBottomTj(coarseToBottomCount), &
            coarseToBottomSi(coarseToBottomCount), coarseToBottomSj(coarseToBottomCount))
        allocate(coarseToBottomSame(coarseToBottomCount))
        allocate(coarseToBottomWx(4, coarseToBottomCount), coarseToBottomWy(4, coarseToBottomCount))
        call fill_fine_interface(nxBottom, nyBottom, xOffsetBottom, yOffsetBottom, coarseToBottomCount, &
            coarseToBottomTi, coarseToBottomTj, coarseToBottomSi, coarseToBottomSj, coarseToBottomWx, &
            coarseToBottomWy, coarseToBottomSame)
        call count_coarse_interface(nxBottom, nyBottom, xOffsetBottom, yOffsetBottom, bottomToCoarseCount)
        allocate(bottomToCoarseTi(bottomToCoarseCount), bottomToCoarseTj(bottomToCoarseCount), &
            bottomToCoarseSi(bottomToCoarseCount), bottomToCoarseSj(bottomToCoarseCount))
        allocate(bottomToCoarseSame(bottomToCoarseCount))
        allocate(bottomToCoarseWx(4, bottomToCoarseCount), bottomToCoarseWy(4, bottomToCoarseCount))
        call fill_coarse_interface(nxBottom, nyBottom, xOffsetBottom, yOffsetBottom, bottomToCoarseCount, &
            bottomToCoarseTi, bottomToCoarseTj, bottomToCoarseSi, bottomToCoarseSj, bottomToCoarseWx, &
            bottomToCoarseWy, bottomToCoarseSame)
        call count_fine_interface(nxTop, nyTop, xOffsetTop, yOffsetTop, coarseToTopCount)
        allocate(coarseToTopTi(coarseToTopCount), coarseToTopTj(coarseToTopCount), &
            coarseToTopSi(coarseToTopCount), coarseToTopSj(coarseToTopCount))
        allocate(coarseToTopSame(coarseToTopCount))
        allocate(coarseToTopWx(4, coarseToTopCount), coarseToTopWy(4, coarseToTopCount))
        call fill_fine_interface(nxTop, nyTop, xOffsetTop, yOffsetTop, coarseToTopCount, coarseToTopTi, &
            coarseToTopTj, coarseToTopSi, coarseToTopSj, coarseToTopWx, coarseToTopWy, coarseToTopSame)
        call count_coarse_interface(nxTop, nyTop, xOffsetTop, yOffsetTop, topToCoarseCount)
        allocate(topToCoarseTi(topToCoarseCount), topToCoarseTj(topToCoarseCount), &
            topToCoarseSi(topToCoarseCount), topToCoarseSj(topToCoarseCount))
        allocate(topToCoarseSame(topToCoarseCount))
        allocate(topToCoarseWx(4, topToCoarseCount), topToCoarseWy(4, topToCoarseCount))
        call fill_coarse_interface(nxTop, nyTop, xOffsetTop, yOffsetTop, topToCoarseCount, topToCoarseTi, &
            topToCoarseTj, topToCoarseSi, topToCoarseSj, topToCoarseWx, topToCoarseWy, topToCoarseSame)
        if (leftToCoarseCount+rightToCoarseCount+bottomToCoarseCount+topToCoarseCount /= &
            nxCoarse*nyCoarse-(nxCoarse-2*interfaceSkin)*(nyCoarse-2*interfaceSkin)) &
            error stop 'Four fine arrays do not cover the complete coarse interface'
        n = max(coarseToLeftCount, leftToCoarseCount, coarseToRightCount, rightToCoarseCount, &
            coarseToBottomCount, bottomToCoarseCount, coarseToTopCount, topToCoarseCount)
        allocate(rhoReceive(n), uReceive(n), vReceive(n), TReceive(n), FxReceive(n), FyReceive(n))
        allocate(flowNeqReceive(0:8,n), thermalNeqReceive(0:4,n))
    end subroutine build_interfaces


    !===============================================================================================
    ! 接口节点与权重常驻设备；各方向按同一异步队列依次使用接收临时数组。
    subroutine interface_device_data(entering)
        use commondata
        implicit none
        logical, intent(in) :: entering

        if (refineRatio == 1) return
        if (entering) then
        !$acc enter data copyin(coarseToLeftTi, coarseToLeftTj, coarseToLeftSi, coarseToLeftSj,  &
        !$acc& coarseToLeftWx, coarseToLeftWy, coarseToLeftSame)
        !$acc enter data copyin(leftToCoarseTi, leftToCoarseTj, leftToCoarseSi, leftToCoarseSj,  &
        !$acc& leftToCoarseWx, leftToCoarseWy, leftToCoarseSame)
        !$acc enter data copyin(coarseToRightTi, coarseToRightTj, coarseToRightSi, coarseToRightSj,  &
        !$acc& coarseToRightWx, coarseToRightWy, coarseToRightSame)
        !$acc enter data copyin(rightToCoarseTi, rightToCoarseTj, rightToCoarseSi, rightToCoarseSj,  &
        !$acc& rightToCoarseWx, rightToCoarseWy, rightToCoarseSame)
        !$acc enter data copyin(coarseToBottomTi, coarseToBottomTj, coarseToBottomSi, coarseToBottomSj,  &
        !$acc& coarseToBottomWx, coarseToBottomWy, coarseToBottomSame)
        !$acc enter data copyin(bottomToCoarseTi, bottomToCoarseTj, bottomToCoarseSi, bottomToCoarseSj,  &
        !$acc& bottomToCoarseWx, bottomToCoarseWy, bottomToCoarseSame)
        !$acc enter data copyin(coarseToTopTi, coarseToTopTj, coarseToTopSi, coarseToTopSj, coarseToTopWx,  &
        !$acc& coarseToTopWy, coarseToTopSame)
        !$acc enter data copyin(topToCoarseTi, topToCoarseTj, topToCoarseSi, topToCoarseSj, topToCoarseWx,  &
        !$acc& topToCoarseWy, topToCoarseSame)
        !$acc enter data create(rhoReceive, uReceive, vReceive, TReceive, FxReceive, FyReceive,  &
        !$acc& flowNeqReceive, thermalNeqReceive)
        else
        !$acc exit data delete(coarseToLeftTi, coarseToLeftTj, coarseToLeftSi, coarseToLeftSj,  &
        !$acc& coarseToLeftWx, coarseToLeftWy, coarseToLeftSame)
        !$acc exit data delete(leftToCoarseTi, leftToCoarseTj, leftToCoarseSi, leftToCoarseSj,  &
        !$acc& leftToCoarseWx, leftToCoarseWy, leftToCoarseSame)
        !$acc exit data delete(coarseToRightTi, coarseToRightTj, coarseToRightSi, coarseToRightSj,  &
        !$acc& coarseToRightWx, coarseToRightWy, coarseToRightSame)
        !$acc exit data delete(rightToCoarseTi, rightToCoarseTj, rightToCoarseSi, rightToCoarseSj,  &
        !$acc& rightToCoarseWx, rightToCoarseWy, rightToCoarseSame)
        !$acc exit data delete(coarseToBottomTi, coarseToBottomTj, coarseToBottomSi, coarseToBottomSj,  &
        !$acc& coarseToBottomWx, coarseToBottomWy, coarseToBottomSame)
        !$acc exit data delete(bottomToCoarseTi, bottomToCoarseTj, bottomToCoarseSi, bottomToCoarseSj,  &
        !$acc& bottomToCoarseWx, bottomToCoarseWy, bottomToCoarseSame)
        !$acc exit data delete(coarseToTopTi, coarseToTopTj, coarseToTopSi, coarseToTopSj, coarseToTopWx,  &
        !$acc& coarseToTopWy, coarseToTopSame)
        !$acc exit data delete(topToCoarseTi, topToCoarseTj, topToCoarseSi, topToCoarseSj, topToCoarseWx,  &
        !$acc& topToCoarseWy, topToCoarseSame)
        !$acc exit data delete(rhoReceive, uReceive, vReceive, TReceive, FxReceive, FyReceive,  &
        !$acc& flowNeqReceive, thermalNeqReceive)
        endif
    end subroutine interface_device_data


    !===============================================================================================
    ! 用粗时间历史插值，分别重建左、右、下、上细网格的粗细人工边缘。
    subroutine coarse_to_fine(wt)
        use commondata
        implicit none
        real(kind=8), intent(in) :: wt(0:2)

        call interpolate_exchange_history(nxCoarse, nyCoarse, historyLastCoarse, rhoHistory_coarse, &
            uHistory_coarse, vHistory_coarse, THistory_coarse, FxHistory_coarse, FyHistory_coarse, &
            flowNeqHistory_coarse, thermalNeqHistory_coarse, coarseToLeftCount, coarseToLeftSi, &
            coarseToLeftSj, coarseToLeftWx, coarseToLeftWy, coarseToLeftSame, wt, rhoReceive, uReceive, &
            vReceive, TReceive, FxReceive, FyReceive, flowNeqReceive, thermalNeqReceive)
        call apply_interface_data(nxLeft, nyLeft, dxFine, SnuFine, SqFine, QkFine, QnuFine, f_left, g_left, &
            rho_left, u_left, v_left, T_left, Fx_left, Fy_left, coarseToLeftCount, coarseToLeftTi, &
            coarseToLeftTj, rhoReceive, uReceive, vReceive, TReceive, FxReceive, FyReceive, flowNeqReceive, &
            thermalNeqReceive)
        call interpolate_exchange_history(nxCoarse, nyCoarse, historyLastCoarse, rhoHistory_coarse, &
            uHistory_coarse, vHistory_coarse, THistory_coarse, FxHistory_coarse, FyHistory_coarse, &
            flowNeqHistory_coarse, thermalNeqHistory_coarse, coarseToRightCount, coarseToRightSi, &
            coarseToRightSj, coarseToRightWx, coarseToRightWy, coarseToRightSame, wt, rhoReceive, uReceive, &
            vReceive, TReceive, FxReceive, FyReceive, flowNeqReceive, thermalNeqReceive)
        call apply_interface_data(nxRight, nyRight, dxFine, SnuFine, SqFine, QkFine, QnuFine, f_right, &
            g_right, rho_right, u_right, v_right, T_right, Fx_right, Fy_right, coarseToRightCount, &
            coarseToRightTi, coarseToRightTj, rhoReceive, uReceive, vReceive, TReceive, FxReceive, &
            FyReceive, flowNeqReceive, thermalNeqReceive)
        call interpolate_exchange_history(nxCoarse, nyCoarse, historyLastCoarse, rhoHistory_coarse, &
            uHistory_coarse, vHistory_coarse, THistory_coarse, FxHistory_coarse, FyHistory_coarse, &
            flowNeqHistory_coarse, thermalNeqHistory_coarse, coarseToBottomCount, coarseToBottomSi, &
            coarseToBottomSj, coarseToBottomWx, coarseToBottomWy, coarseToBottomSame, wt, rhoReceive, &
            uReceive, vReceive, TReceive, FxReceive, FyReceive, flowNeqReceive, thermalNeqReceive)
        call apply_interface_data(nxBottom, nyBottom, dxFine, SnuFine, SqFine, QkFine, QnuFine, f_bottom, &
            g_bottom, rho_bottom, u_bottom, v_bottom, T_bottom, Fx_bottom, Fy_bottom, coarseToBottomCount, &
            coarseToBottomTi, coarseToBottomTj, rhoReceive, uReceive, vReceive, TReceive, FxReceive, &
            FyReceive, flowNeqReceive, thermalNeqReceive)
        call interpolate_exchange_history(nxCoarse, nyCoarse, historyLastCoarse, rhoHistory_coarse, &
            uHistory_coarse, vHistory_coarse, THistory_coarse, FxHistory_coarse, FyHistory_coarse, &
            flowNeqHistory_coarse, thermalNeqHistory_coarse, coarseToTopCount, coarseToTopSi, &
            coarseToTopSj, coarseToTopWx, coarseToTopWy, coarseToTopSame, wt, rhoReceive, uReceive, &
            vReceive, TReceive, FxReceive, FyReceive, flowNeqReceive, thermalNeqReceive)
        call apply_interface_data(nxTop, nyTop, dxFine, SnuFine, SqFine, QkFine, QnuFine, f_top, g_top, &
            rho_top, u_top, v_top, T_top, Fx_top, Fy_top, coarseToTopCount, coarseToTopTi, coarseToTopTj, &
            rhoReceive, uReceive, vReceive, TReceive, FxReceive, FyReceive, flowNeqReceive, thermalNeqReceive)
    end subroutine coarse_to_fine


    !===============================================================================================
    ! 在同步时刻，用四套细网格的当前交换数据重建中心粗网格边缘。
    subroutine fine_to_coarse()
        use commondata
        implicit none

        call interpolate_exchange_history(nxLeft, nyLeft, 0, rhoHistory_left, uHistory_left, vHistory_left, &
            THistory_left, FxHistory_left, FyHistory_left, flowNeqHistory_left, thermalNeqHistory_left, &
            leftToCoarseCount, leftToCoarseSi, leftToCoarseSj, leftToCoarseWx, leftToCoarseWy, &
            leftToCoarseSame, [1.0d0, 0.0d0, 0.0d0], rhoReceive, uReceive, vReceive, TReceive, FxReceive, &
            FyReceive, flowNeqReceive, thermalNeqReceive)
        call apply_interface_data(nxCoarse, nyCoarse, dxCoarse, SnuCoarse, SqCoarse, QkCoarse, QnuCoarse, &
            f_coarse, g_coarse, rho_coarse, u_coarse, v_coarse, T_coarse, Fx_coarse, Fy_coarse, &
            leftToCoarseCount, leftToCoarseTi, leftToCoarseTj, rhoReceive, uReceive, vReceive, TReceive, &
            FxReceive, FyReceive, flowNeqReceive, thermalNeqReceive)
        call interpolate_exchange_history(nxRight, nyRight, 0, rhoHistory_right, uHistory_right, &
            vHistory_right, THistory_right, FxHistory_right, FyHistory_right, flowNeqHistory_right, &
            thermalNeqHistory_right, rightToCoarseCount, rightToCoarseSi, rightToCoarseSj, rightToCoarseWx, &
            rightToCoarseWy, rightToCoarseSame, [1.0d0, 0.0d0, 0.0d0], rhoReceive, uReceive, vReceive, &
            TReceive, FxReceive, FyReceive, flowNeqReceive, thermalNeqReceive)
        call apply_interface_data(nxCoarse, nyCoarse, dxCoarse, SnuCoarse, SqCoarse, QkCoarse, QnuCoarse, &
            f_coarse, g_coarse, rho_coarse, u_coarse, v_coarse, T_coarse, Fx_coarse, Fy_coarse, &
            rightToCoarseCount, rightToCoarseTi, rightToCoarseTj, rhoReceive, uReceive, vReceive, TReceive, &
            FxReceive, FyReceive, flowNeqReceive, thermalNeqReceive)
        call interpolate_exchange_history(nxBottom, nyBottom, 0, rhoHistory_bottom, uHistory_bottom, &
            vHistory_bottom, THistory_bottom, FxHistory_bottom, FyHistory_bottom, flowNeqHistory_bottom, &
            thermalNeqHistory_bottom, bottomToCoarseCount, bottomToCoarseSi, bottomToCoarseSj, &
            bottomToCoarseWx, bottomToCoarseWy, bottomToCoarseSame, [1.0d0, 0.0d0, 0.0d0], rhoReceive, &
            uReceive, vReceive, TReceive, FxReceive, FyReceive, flowNeqReceive, thermalNeqReceive)
        call apply_interface_data(nxCoarse, nyCoarse, dxCoarse, SnuCoarse, SqCoarse, QkCoarse, QnuCoarse, &
            f_coarse, g_coarse, rho_coarse, u_coarse, v_coarse, T_coarse, Fx_coarse, Fy_coarse, &
            bottomToCoarseCount, bottomToCoarseTi, bottomToCoarseTj, rhoReceive, uReceive, vReceive, &
            TReceive, FxReceive, FyReceive, flowNeqReceive, thermalNeqReceive)
        call interpolate_exchange_history(nxTop, nyTop, 0, rhoHistory_top, uHistory_top, vHistory_top, &
            THistory_top, FxHistory_top, FyHistory_top, flowNeqHistory_top, thermalNeqHistory_top, &
            topToCoarseCount, topToCoarseSi, topToCoarseSj, topToCoarseWx, topToCoarseWy, topToCoarseSame, &
            [1.0d0, 0.0d0, 0.0d0], rhoReceive, uReceive, vReceive, TReceive, FxReceive, FyReceive, &
            flowNeqReceive, thermalNeqReceive)
        call apply_interface_data(nxCoarse, nyCoarse, dxCoarse, SnuCoarse, SqCoarse, QkCoarse, QnuCoarse, &
            f_coarse, g_coarse, rho_coarse, u_coarse, v_coarse, T_coarse, Fx_coarse, Fy_coarse, &
            topToCoarseCount, topToCoarseTi, topToCoarseTj, rhoReceive, uReceive, vReceive, TReceive, &
            FxReceive, FyReceive, flowNeqReceive, thermalNeqReceive)
    end subroutine fine_to_coarse


    ! 粗人工边缘的 interfaceSkin 层由细网格重建；物理墙面不在这里处理。
    logical function coarse_skin(i,j)
        use commondata, only: nxCoarse,nyCoarse,interfaceSkin,refineRatio
        implicit none
        integer, intent(in) :: i,j
        coarse_skin = refineRatio>1 .and. (i<=interfaceSkin .or. i>nxCoarse-interfaceSkin .or. &
            j<=interfaceSkin .or. j>nyCoarse-interfaceSkin)
    end function coarse_skin

    ! 用全局坐标识别细环内缘；左右上下数组的拼接线不是粗细人工边缘。
    logical function fine_skin(x,y)
        use commondata, only: centerBox,overlapCells,refineRatio,interfaceSkin
        implicit none
        real(8), intent(in) :: x,y
        real(8) :: xl,xr,yb,yt,o
        o=dble(overlapCells*refineRatio)
        xl=centerBox(1)+o; xr=centerBox(2)-o
        yb=centerBox(3)+o; yt=centerBox(4)-o
        fine_skin=(x<=xl .or. x>=xr .or. y<=yb .or. y>=yt) .and. &
            x>xl-interfaceSkin .and. x<xr+interfaceSkin .and. &
            y>yb-interfaceSkin .and. y<yt+interfaceSkin
    end function fine_skin

    ! 数出当前细矩形中需要由粗网格提供数据的节点。
    subroutine count_fine_interface(ni,nj,xOffset,yOffset,count)
        use commondata, only: fine_skin
        implicit none
        integer, intent(in) :: ni,nj
        real(8), intent(in) :: xOffset,yOffset
        integer, intent(out) :: count
        integer :: i,j
        count=0
        do j=1,nj
            do i=1,ni
                if (fine_skin(xOffset+i-.5d0,yOffset+j-.5d0)) count=count+1
            enddo
        enddo
    end subroutine count_fine_interface

    ! 记录细接收节点及对应的粗四点模板，供每步交换复用。
    subroutine fill_fine_interface(ni,nj,xOffset,yOffset,count,ti,tj,si,sj,wx,wy,same)
        use commondata, only: fine_skin
        implicit none
        integer, intent(in) :: ni,nj,count
        real(8), intent(in) :: xOffset,yOffset
        integer, intent(out) :: ti(count),tj(count),si(count),sj(count)
        real(8), intent(out) :: wx(4,count),wy(4,count)
        logical, intent(out) :: same(count)
        integer :: i,j,c
        real(8) :: x,y
        c=0
        do j=1,nj
            y=yOffset+j-.5d0
            do i=1,ni
                x=xOffset+i-.5d0
                if (.not.fine_skin(x,y)) cycle
                c=c+1; ti(c)=i; tj(c)=j
                call coarse_donor_stencil(x,y,si(c),sj(c),wx(:,c),wy(:,c),same(c))
            enddo
        enddo
        if(c/=count) error stop 'Fine interface count mismatch'
    end subroutine fill_fine_interface

    ! 根据物理坐标定位粗节点；模板必须避开粗人工缓冲层。
    subroutine coarse_donor_stencil(x,y,si,sj,wx,wy,same)
        use commondata, only: nxCoarse,nyCoarse,xOffsetCoarse,yOffsetCoarse,dxCoarse,interfaceSkin
        implicit none
        real(8), intent(in) :: x,y
        integer, intent(out) :: si,sj
        real(8), intent(out) :: wx(4),wy(4)
        logical, intent(out) :: same
        integer :: il,ih,jl,jh
        real(8) :: qx,qy
        il=1+interfaceSkin; ih=nxCoarse-interfaceSkin
        jl=1+interfaceSkin; jh=nyCoarse-interfaceSkin
        qx=(x-xOffsetCoarse)/dxCoarse+.5d0
        qy=(y-yOffsetCoarse)/dxCoarse+.5d0
        if(qx<il .or. qx>ih .or. qy<jl .or. qy>jh .or. min(ih-il,jh-jl)<3) &
            error stop 'Invalid overlap: no interior four-point coarse donor stencil'
        same=abs(qx-nint(qx))<1d-12 .and. abs(qy-nint(qy))<1d-12
        if(same) then
            si=nint(qx); sj=nint(qy)
            wx=[1d0,0d0,0d0,0d0]; wy=wx
        else
            si=max(il,min(floor(qx)-1,ih-3)); sj=max(jl,min(floor(qy)-1,jh-3))
            call lagrange_weights(qx-si,wx)
            call lagrange_weights(qy-sj,wy)
        endif
    end subroutine coarse_donor_stencil

    ! 数出粗人工边缘中坐标落入当前细矩形的节点。
    subroutine count_coarse_interface(ni,nj,xOffset,yOffset,count)
        use commondata, only: nxCoarse,nyCoarse,xOffsetCoarse,yOffsetCoarse,dxCoarse,coarse_skin
        implicit none
        integer, intent(in) :: ni,nj
        real(8), intent(in) :: xOffset,yOffset
        integer, intent(out) :: count
        integer :: i,j,si,sj
        real(8) :: x,y
        count=0
        do j=1,nyCoarse
            y=yOffsetCoarse+(j-.5d0)*dxCoarse
            do i=1,nxCoarse
                if(.not.coarse_skin(i,j)) cycle
                x=xOffsetCoarse+(i-.5d0)*dxCoarse
                si=nint(x-xOffset+.5d0); sj=nint(y-yOffset+.5d0)
                if(si<1 .or. si>ni .or. sj<1 .or. sj>nj) cycle
                count=count+1
            enddo
        enddo
    end subroutine count_coarse_interface

    ! 记录粗接收节点与共址细来源节点，禁止从细人工缓冲层取值。
    subroutine fill_coarse_interface(ni,nj,xOffset,yOffset,count,ti,tj,si,sj,wx,wy,same)
        use commondata, only: nxCoarse,nyCoarse,xOffsetCoarse,yOffsetCoarse,dxCoarse,coarse_skin,fine_skin
        implicit none
        integer, intent(in) :: ni,nj,count
        real(8), intent(in) :: xOffset,yOffset
        integer, intent(out) :: ti(count),tj(count),si(count),sj(count)
        real(8), intent(out) :: wx(4,count),wy(4,count)
        logical, intent(out) :: same(count)
        integer :: i,j,c,is,js
        real(8) :: x,y,qx,qy
        c=0
        do j=1,nyCoarse
            y=yOffsetCoarse+(j-.5d0)*dxCoarse
            do i=1,nxCoarse
                if(.not.coarse_skin(i,j)) cycle
                x=xOffsetCoarse+(i-.5d0)*dxCoarse
                qx=x-xOffset+.5d0; qy=y-yOffset+.5d0
                is=nint(qx); js=nint(qy)
                if(is<1 .or. is>ni .or. js<1 .or. js>nj) cycle
                if(abs(qx-is)>1d-12 .or. abs(qy-js)>1d-12 .or. fine_skin(x,y)) &
                    error stop 'Coarse receiver must coincide with an interior fine donor'
                c=c+1; ti(c)=i; tj(c)=j; si(c)=is; sj(c)=js
                wx(:,c)=[1d0,0d0,0d0,0d0]; wy(:,c)=wx(:,c); same(c)=.true.
            enddo
        enddo
        if(c/=count) error stop 'Coarse interface count mismatch'
    end subroutine fill_coarse_interface


    subroutine collision(nx, ny, f, f_post, rho, u, v, Fx, Fy, T, Snu, Sq, gBeta &
#ifdef SideHeatedHa
        , B2sigemarho &
#endif
        )

        use commondata, only: Tref
#ifdef SideHeatedHa
        use commondata, only: phi
#endif
        implicit none

        integer(kind=4), intent(in) :: nx, ny
        real(kind=8), intent(inout) :: f(nx, ny, 0:8), f_post(0:nx+1, 0:ny+1, 0:8)
        real(kind=8), intent(inout) :: u(nx, ny), v(nx, ny), T(nx, ny), rho(nx, ny), Fx(nx, ny), Fy(nx, ny)
        real(kind=8), intent(in) :: Snu, Sq, gBeta
#ifdef SideHeatedHa
        real(kind=8), intent(in) :: B2sigemarho
#endif

        integer(kind=4) :: i, j
        integer(kind=4) :: alpha
        real(kind=8) :: m(0:8), m_post(0:8), meq(0:8)
        real(kind=8) :: s(0:8)
        real(kind=8) :: fSource(0:8)

        !$acc parallel loop gang vector collapse(2) present(f, f_post, rho, u, v, Fx, Fy, T) async(1) &
        !$acc& private(alpha, s, m, m_post, meq, fSource)
        do j = 1, ny
            do i = 1, nx

                m(0) = f(i, j, 0)+f(i, j, 1)+f(i, j, 2)+f(i, j, 3)+f(i, j, 4)+f(i, j, 5)+f(i, j, 6)+f(i, j, &
                    7)+f(i, j, 8)
                m(1) = -4.0d0*f(i, j, 0)-f(i, j, 1)-f(i, j, 2)-f(i, j, 3)-f(i, j, 4)+2.0d0*(f(i, j, 5)+f(i, j, &
                    6)+f(i, j, 7)+f(i, j, 8))
                m(2) = 4.0d0*f(i, j, 0)-2.0d0*(f(i, j, 1)+f(i, j, 2)+f(i, j, 3)+f(i, j, 4))+f(i, j, 5)+f(i, j, &
                    6)+f(i, j, 7)+f(i, j, 8)
                m(3) = f(i, j, 1)-f(i, j, 3)+f(i, j, 5)-f(i, j, 6)-f(i, j, 7)+f(i, j, 8)
                m(4) = -2.0d0*f(i, j, 1)+2.0d0*f(i, j, 3)+f(i, j, 5)-f(i, j, 6)-f(i, j, 7)+f(i, j, 8)
                m(5) = f(i, j, 2)-f(i, j, 4)+f(i, j, 5)+f(i, j, 6)-f(i, j, 7)-f(i, j, 8)
                m(6) = -2.0d0*f(i, j, 2)+2.0d0*f(i, j, 4)+f(i, j, 5)+f(i, j, 6)-f(i, j, 7)-f(i, j, 8)
                m(7) = f(i, j, 1)-f(i, j, 2)+f(i, j, 3)-f(i, j, 4)
                m(8) = f(i, j, 5)-f(i, j, 6)+f(i, j, 7)-f(i, j, 8)

                meq(0) = rho(i, j)
                meq(1) = rho(i, j)*( -2.0d0+3.0d0*(u(i, j)*u(i, j)+v(i, j)*v(i, j)) )
                meq(2) = rho(i, j)*( 1.0d0-3.0d0*(u(i, j)*u(i, j)+v(i, j)*v(i, j)) )
                meq(3) = rho(i, j)*u(i, j)
                meq(4) = -rho(i, j)*u(i, j)
                meq(5) = rho(i, j)*v(i, j)
                meq(6) = -rho(i, j)*v(i, j)
                meq(7) = rho(i, j)*( u(i, j)*u(i, j)-v(i, j)*v(i, j) )
                meq(8) = rho(i, j)*( u(i, j)*v(i, j) )

                s(0) = 0.0d0    !!s_{\rho}
                s(1) = Snu    !!s_{e}
                s(2) = Snu    !!s_{\epsilon}
                s(3) = 0.0d0    !!s_{j}
                s(4) = Sq    !!s_{q}
                s(5) = 0.0d0    !!s_{j}
                s(6) = Sq    !!s_{q}
                s(7) = Snu    !!s_{\nu}
                s(8) = Snu    !!s_{\nu}

                Fx(i, j) = 0.0d0
                Fy(i, j) = rho(i, j)*gBeta*(T(i, j)-Tref)    !动量方程上的源项，即浮力项

#ifdef    SideHeatedHa
                Fx(i, j) = 0.0d0+B2sigemarho*(v(i, j)*sin(phi)*cos(phi)-u(i, j)*sin(phi)*sin(phi))
                Fy(i, j) = rho(i, j)*gBeta*(T(i, j)-Tref)+ rho(i, j)*B2sigemarho*(u(i, j)*sin(phi)*cos(phi)&
                    -v(i, j)*cos(phi)*cos(phi))    !动量方程上的源项，即浮力项加磁场
#endif

                fSource(0) = 0.0d0    !将源项F对应的贡献投影到各个矩中，并做半步修正
                fSource(1) = (6.0d0-3.0d0*s(1))*(u(i, j)*Fx(i, j)+v(i, j)*Fy(i, j))
                fSource(2) = -(6.0d0-3.0d0*s(2))*(u(i, j)*Fx(i, j)+v(i, j)*Fy(i, j))
                fSource(3) = (1.0d0-0.5d0*s(3))*Fx(i, j)
                fSource(4) = -(1.0d0-0.5d0*s(4))*Fx(i, j)
                fSource(5) = (1.0d0-0.5d0*s(5))*Fy(i, j)
                fSource(6) = -(1.0d0-0.5d0*s(6))*Fy(i, j)
                fSource(7) = (2.0d0-s(7))*(u(i, j)*Fx(i, j)-v(i, j)*Fy(i, j))
                fSource(8) = (1.0d0-0.5d0*s(8))*(u(i, j)*Fy(i, j)+v(i, j)*Fx(i, j))    !这边是乘以M变到矩空间，然后再乘以1-1/2S修正

                do alpha = 0, 8
                    m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+fSource(alpha)    !矩空间碰撞
                enddo

                f_post(i, j, 0) = m_post(0)/9.0d0-m_post(1)/9.0d0+m_post(2)/9.0d0    !这边是乘以M逆
                f_post(i, j, &
            1) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+m_post(3)/6.0d0-m_post(4)/6.0d0 &
                    +m_post(7)/4.0d0
                f_post(i, j, 2) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    +m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)/4.0d0
                f_post(i, j, &
            3) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0-m_post(3)/6.0d0+m_post(4)/6.0d0 &
                    +m_post(7)/4.0d0
                f_post(i, j, 4) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    -m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)/4.0d0
                f_post(i, j, &
            5) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)/4.0d0
                f_post(i, j, &
            6) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)/4.0d0
                f_post(i, j, &
            7) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)/4.0d0
                f_post(i, j, &
            8) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)/4.0d0

            enddo
        enddo
        return
    end subroutine collision


    subroutine streaming(nx, ny, f, f_post)    !先迁移，再边界处理

        use commondata, only: ex, ey
        implicit none

        integer(kind=4), intent(in) :: nx, ny
        real(kind=8), intent(inout) :: f(nx, ny, 0:8), f_post(0:nx+1, 0:ny+1, 0:8)

        integer(kind=4) :: i, j
        integer(kind=4) :: ip, jp
        integer(kind=4) :: alpha

        !$acc parallel loop gang vector collapse(2) present(f, f_post, ex, ey) async(1) private(alpha, ip, jp)
        do j = 1, ny
            do i = 1, nx
                do alpha = 0, 8    !上游格点索引：fα(i,j) <- f_postα(i-exα, j-eyα)
                    ip = i-ex(alpha)    !边界附近 (ip/jp 可能为 0 或 nx+1/ny+1)，需在 bounceback/周期边界处理中覆盖修正边界分布
                    jp = j-ey(alpha)    !ghost 层在初始化中为 0，保证不会出现未初始化垃圾值

                    f(i, j, alpha) = f_post(ip, jp, alpha)
                enddo
            enddo
        enddo
        return
    end subroutine streaming


    subroutine bounceback(nx, ny, f, f_post, leftWall, rightWall, bottomWall, topWall)

        implicit none

        integer(kind=4), intent(in) :: nx, ny
        real(kind=8), intent(inout) :: f(nx, ny, 0:8), f_post(0:nx+1, 0:ny+1, 0:8)
        logical, intent(in) :: leftWall, rightWall, bottomWall, topWall

        integer(kind=4) :: i, j
        ! integer(kind=4) :: alpha

#ifdef VerticalWallsPeriodicalU
        !$acc parallel loop gang vector present(f, f_post) async(1)
        do j = 1, ny    !速度边界垂直边界周期，直接方向相同，跨边界的入射分布
            !Left side (i=1)
            if (leftWall) f(1, j, 1) = f_post(nx, j, 1)
            if (leftWall) f(1, j, 5) = f_post(nx, j, 5)
            if (leftWall) f(1, j, 8) = f_post(nx, j, 8)

            !Right side (i=nx)
            if (rightWall) f(nx, j, 3) = f_post(1, j, 3)
            if (rightWall) f(nx, j, 6) = f_post(1, j, 6)
            if (rightWall) f(nx, j, 7) = f_post(1, j, 7)
        enddo
#endif

#ifdef VerticalWallsNoslip
        !$acc parallel loop gang vector present(f, f_post) async(1)
        do j = 1, ny    !速度边界垂直边界静止壁无滑移，直接反弹，方向相反
            !Left side (i=1)
            if (leftWall) f(1, j, 1) = f_post(1, j, 3)
            if (leftWall) f(1, j, 5) = f_post(1, j, 7)
            if (leftWall) f(1, j, 8) = f_post(1, j, 6)

            !Right side (i=nx)
            if (rightWall) f(nx, j, 3) = f_post(nx, j, 1)
            if (rightWall) f(nx, j, 6) = f_post(nx, j, 8)
            if (rightWall) f(nx, j, 7) = f_post(nx, j, 5)
        enddo
#endif

#ifdef HorizontalWallsNoslip
        !$acc parallel loop gang vector present(f, f_post) async(1)
        do i = 1, nx    !速度边界水平边界无滑移，直接反弹，方向相反
            !Bottom side (j=1)
            if (bottomWall) f(i, 1, 2) = f_post(i, 1, 4)
            if (bottomWall) f(i, 1, 5) = f_post(i, 1, 7)
            if (bottomWall) f(i, 1, 6) = f_post(i, 1, 8)

            !Top side (j=ny)
            if (topWall) f(i, ny, 4) = f_post(i, ny, 2)
            if (topWall) f(i, ny, 7) = f_post(i, ny, 5)
            if (topWall) f(i, ny, 8) = f_post(i, ny, 6)
        enddo
#endif

        return
    end subroutine bounceback


    subroutine macro(nx, ny, f, rho, u, v, Fx, Fy)

        implicit none

        integer(kind=4), intent(in) :: nx, ny
        real(kind=8), intent(in) :: f(nx, ny, 0:8), Fx(nx, ny), Fy(nx, ny)
        real(kind=8), intent(inout) :: rho(nx, ny), u(nx, ny), v(nx, ny)

        integer(kind=4) :: i, j

        !$acc parallel loop gang vector collapse(2) present(f, rho, u, v, Fx, Fy) async(1)
        do j = 1, ny
            do i = 1, nx
                rho(i, j) = f(i, j, 0)+f(i, j, 1)+f(i, j, 2)+f(i, j, 3)+f(i, j, 4)+f(i, j, 5)+f(i, j, 6)+f(i, j, &
                    7)+f(i, j, 8)
                u(i, j) = ( f(i, j, 1)-f(i, j, 3)+f(i, j, 5)-f(i, j, 6)-f(i, j, 7)+f(i, j, 8)+0.5d0*Fx(i, &
                    j) )/rho(i, j)    !含力LBM的半步动量修正：rho*u = Σ f e + 0.5*F，对应Guo forcing的二阶定义
                v(i, j) = ( f(i, j, 2)-f(i, j, 4)+f(i, j, 5)+f(i, j, 6)-f(i, j, 7)-f(i, j, 8)+0.5d0*Fy(i, &
            j) )/rho(i, j)
            enddo
        enddo
        return
    end subroutine macro


    subroutine collisionT(nx, ny, g, g_post, u, v, T, Qk, Qnu)

        use commondata, only: paraA
        implicit none

        integer(kind=4), intent(in) :: nx, ny
        real(kind=8), intent(inout) :: g(nx, ny, 0:4), g_post(0:nx+1, 0:ny+1, 0:4)
        real(kind=8), intent(in) :: u(nx, ny), v(nx, ny), T(nx, ny), Qk, Qnu

        integer(kind=4) :: i, j
        integer(kind=4) :: alpha
        real(kind=8) :: n(0:4), n_post(0:4), neq(0:4)
        real(kind=8) :: q(0:4)
        !$acc parallel loop gang vector collapse(2) present(g, g_post, u, v, T) async(1) &
        !$acc& private(alpha, n, neq, q, n_post)
        do j = 1, ny
            do i = 1, nx

                n(0) = g(i, j, 0)+g(i, j, 1)+g(i, j, 2)+g(i, j, 3)+g(i, j, 4)
                n(1) = g(i, j, 1)-g(i, j, 3)
                n(2) = g(i, j, 2)-g(i, j, 4)
                n(3) = -4.0d0*g(i, j, 0)+g(i, j, 1)+g(i, j, 2)+g(i, j, 3)+g(i, j, 4)
                n(4) = g(i, j, 1)-g(i, j, 2)+g(i, j, 3)-g(i, j, 4)

                neq(0) = T(i, j)
                neq(1) = T(i, j)*u(i, j)
                neq(2) = T(i, j)*v(i, j)
                neq(3) = T(i, j)*paraA
                neq(4) = 0.0d0

                q(0) = 0.0d0
                q(1) = Qk
                q(2) = Qk
                q(3) = Qnu
                q(4) = Qnu

                n_post(0) = n(0)-q(0)*(n(0)-neq(0))
                n_post(1) = n(1)-q(1)*(n(1)-neq(1))
                n_post(2) = n(2)-q(2)*(n(2)-neq(2))
                n_post(3) = n(3)-q(3)*(n(3)-neq(3))
                n_post(4) = n(4)-q(4)*(n(4)-neq(4))

                g_post(i, j, 0) = 0.2d0*n_post(0)-0.2d0*n_post(3)
                g_post(i, j, 1) = 0.2d0*n_post(0)+0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
                g_post(i, j, 2) = 0.2d0*n_post(0)+0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4)
                g_post(i, j, 3) = 0.2d0*n_post(0)-0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
                g_post(i, j, 4) = 0.2d0*n_post(0)-0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4)
            enddo
        enddo
        return
    end subroutine collisionT


    subroutine streamingT(nx, ny, g, g_post)

        use commondata, only: ex, ey
        implicit none

        integer(kind=4), intent(in) :: nx, ny
        real(kind=8), intent(inout) :: g(nx, ny, 0:4), g_post(0:nx+1, 0:ny+1, 0:4)

        integer(kind=4) :: i, j
        integer(kind=4) :: ip, jp
        integer(kind=4) :: alpha

        !$acc parallel loop gang vector collapse(2) present(g, g_post, ex, ey) async(1) private(alpha, ip, jp)
        do j = 1, ny
            do i = 1, nx
                do alpha = 0, 4
                    ip = i-ex(alpha)
                    jp = j-ey(alpha)

                    g(i, j, alpha) = g_post(ip, jp, alpha)
                enddo
            enddo
        enddo
        return
    end subroutine streamingT


    subroutine bouncebackT(nx, ny, g, g_post, leftWall, rightWall, bottomWall, topWall)

        use commondata, only: Thot, Tcold, paraA, omegaT
        implicit none

        integer(kind=4), intent(in) :: nx, ny
        real(kind=8), intent(inout) :: g(nx, ny, 0:4), g_post(0:nx+1, 0:ny+1, 0:4)
        logical, intent(in) :: leftWall, rightWall, bottomWall, topWall

        integer(kind=4) :: i, j
        !integer(kind=4) :: alpha

#ifdef VerticalWallsPeriodicalT
        !$acc parallel loop gang vector present(g, g_post) async(1)
        do j = 1, ny
            !Left boundary
            if (leftWall) g(1, j, 1) = g_post(nx, j, 1)

            !Right boundary
            if (rightWall) g(nx, j, 3) = g_post(1, j, 3)
        enddo
#endif

#ifdef VerticalWallsConstT
        !$acc parallel loop gang vector present(g, g_post, omegaT) async(1)
        do j = 1, ny
            !Left boundary
            if (leftWall) g(1, j, 1) = -g_post(1, j, 3)+(4.0d0+paraA)/10.0d0*Thot
            !Right boundary
            if (rightWall) g(nx, j, 3) = -g_post(nx, j, 1)+(4.0d0+paraA)/10.0d0*Tcold
        enddo
#endif

#ifdef VerticalWallsAdiabatic
        !$acc parallel loop gang vector present(g, g_post) async(1)
        do j = 1, ny
            !Left boundary
            if (leftWall) g(1, j, 1) = g_post(1, j, 3)

            !Right boundary
            if (rightWall) g(nx, j, 3) = g_post(nx, j, 1)
        enddo
#endif

#ifdef HorizontalWallsAdiabatic
        !$acc parallel loop gang vector present(g, g_post) async(1)
        do i = 1, nx
            !Bottom side
            if (bottomWall) g(i, 1, 2) = g_post(i, 1, 4)

            !Top side
            if (topWall) g(i, ny, 4) = g_post(i, ny, 2)
        enddo
#endif

#ifdef HorizontalWallsConstT
        !$acc parallel loop gang vector present(g, g_post, omegaT) async(1)
        do i = 1, nx
            !Bottom side
            if (bottomWall) g(i, 1, 2) = -g_post(i, 1, 4)+(4.0d0+paraA)/10.0d0*Thot
            !Top side
            if (topWall) g(i, ny, 4) = -g_post(i, ny, 2)+(4.0d0+paraA)/10.0d0*Tcold
        enddo
#endif

        return
    end subroutine bouncebackT


    subroutine macroT(nx, ny, g, T)

        implicit none

        integer(kind=4), intent(in) :: nx, ny
        real(kind=8), intent(in) :: g(nx, ny, 0:4)
        real(kind=8), intent(inout) :: T(nx, ny)

        integer(kind=4) :: i, j

        !$acc parallel loop gang vector collapse(2) present(g, T) async(1)
        do j = 1, ny
            do i = 1, nx
                T(i, j) = g(i, j, 0)+g(i, j, 1)+g(i, j, 2)+g(i, j, 3)+g(i, j, 4)
            enddo
        enddo
        return
    end subroutine macroT


    subroutine flow_moments(fv, m)

        implicit none

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


    subroutine flow_populations(m, fv)

        implicit none

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


    subroutine thermal_moments(gv, n)

        implicit none

        !$acc routine seq
        real(kind=8), intent(in) :: gv(0:4)
        real(kind=8), intent(out) :: n(0:4)

        n(0) = gv(0)+gv(1)+gv(2)+gv(3)+gv(4)
        n(1) = gv(1)-gv(3)
        n(2) = gv(2)-gv(4)
        n(3) = -4.0d0*gv(0)+gv(1)+gv(2)+gv(3)+gv(4)
        n(4) = gv(1)-gv(2)+gv(3)-gv(4)

    end subroutine thermal_moments


    subroutine thermal_populations(n, gv)

        implicit none

        !$acc routine seq
        real(kind=8), intent(in) :: n(0:4)
        real(kind=8), intent(out) :: gv(0:4)

        gv(0) = 0.2d0*n(0)-0.2d0*n(3)
        gv(1) = 0.2d0*n(0)+0.5d0*n(1)+0.05d0*n(3)+0.25d0*n(4)
        gv(2) = 0.2d0*n(0)+0.5d0*n(2)+0.05d0*n(3)-0.25d0*n(4)
        gv(3) = 0.2d0*n(0)-0.5d0*n(1)+0.05d0*n(3)+0.25d0*n(4)
        gv(4) = 0.2d0*n(0)-0.5d0*n(2)+0.05d0*n(3)-0.25d0*n(4)

    end subroutine thermal_populations


    ! 同步主机场量，依次累加五个积分分区的 Nu、Re 和诊断量。
    subroutine calNuRe()

        use commondata
        implicit none

        integer(kind=4) :: i, j, k, jm, im
        real(kind=8) :: area, conv, vel2, mass, meanT, nu, re, hot, cold, middle, dTdx, dTdy, tm, um, vm, cellArea
        real(kind=8) :: w(4), dw(4)
        real(kind=8) :: tmin, tmax, rmin, rmax, scale, xmid, ymid, xlo, xhi, ylo, yhi, dx

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
        call calNuRe_grid(nxCoarse, nyCoarse, dxCoarse, xOffsetCoarse, yOffsetCoarse, iFirstCoarse, &
            iLastCoarse, jFirstCoarse, jLastCoarse, ownedBoxCoarse, conv, vel2, mass, meanT, hot, cold, &
            middle, tmin, tmax, rmin, rmax, xmid, ymid, rho_coarse, u_coarse, v_coarse, T_coarse, &
            quadWidthX_coarse, quadWidthY_coarse)
        if (refineRatio > 1) then
            call calNuRe_grid(nxLeft, nyLeft, dxFine, xOffsetLeft, yOffsetLeft, iFirstLeft, iLastLeft, &
                jFirstLeft, jLastLeft, ownedBoxLeft, conv, vel2, mass, meanT, hot, cold, middle, tmin, tmax, &
                rmin, rmax, xmid, ymid, rho_left, u_left, v_left, T_left, quadWidthX_left, quadWidthY_left)
            call calNuRe_grid(nxRight, nyRight, dxFine, xOffsetRight, yOffsetRight, iFirstRight, iLastRight, &
                jFirstRight, jLastRight, ownedBoxRight, conv, vel2, mass, meanT, hot, cold, middle, tmin, tmax, &
                rmin, rmax, xmid, ymid, rho_right, u_right, v_right, T_right, quadWidthX_right, quadWidthY_right)
            call calNuRe_grid(nxBottom, nyBottom, dxFine, xOffsetBottom, yOffsetBottom, iFirstBottom, &
                iLastBottom, jFirstBottom, jLastBottom, ownedBoxBottom, conv, vel2, mass, meanT, hot, cold, &
                middle, tmin, tmax, rmin, rmax, xmid, ymid, rho_bottom, u_bottom, v_bottom, T_bottom, &
                quadWidthX_bottom, quadWidthY_bottom)
            call calNuRe_grid(nxTop, nyTop, dxFine, xOffsetTop, yOffsetTop, iFirstTop, iLastTop, jFirstTop, &
                jLastTop, ownedBoxTop, conv, vel2, mass, meanT, hot, cold, middle, tmin, tmax, rmin, rmax, &
                xmid, ymid, rho_top, u_top, v_top, T_top, quadWidthX_top, quadWidthY_top)
        endif

        ! 空间积分只覆盖 ownedBox；重叠节点按各自分区面积计权，不重复计算整块面积。
        ! Re 使用速度平方的面积平均再开方，非稳态时间统计也采用同一 RMS 定义。
        nu = 1.0d0+conv/area*scale/diffusivity
        re = sqrt(vel2/area)*lengthUnit/viscosity
        open(newunit = k, file = NuReHistoryFile, status = 'old', position = 'append')
        write(k, '(12(ES24.16E3,1X))') dble(itc)/timeUnit, nu, re, scale*hot, scale*cold, scale*middle, &
            mass, meanT/area, tmin, tmax, rmin, rmax
        close(k)
        write(*, '(a,f12.5,a,es13.5,a,es13.5)') 't_ff=', dble(itc)/timeUnit, ' NuVolAvg=', nu, ' ReVolRMS=', re
    end subroutine calNuRe


    ! 对传入数组积分；细区扣除中心面积，中线模板可以跨越细数组接缝。
    subroutine calNuRe_grid(ni, nj, dx, xOffset, yOffset, iFirst, iLast, jFirst, jLast, ownedBox, conv, &
            vel2, mass, meanT, hot, cold, middle, tmin, tmax, rmin, rmax, xmid, ymid, &
            rho, u, v, T, quadWidthX, quadWidthY)

        use commondata, only: nx, ny, Thot, Tcold, diffusivity, itc, refineRatio, ieee_is_finite, &
            owned_cell_area, section_owned_weight
        implicit none
        integer(kind=4) :: ni, nj, iFirst, iLast, jFirst, jLast
        real(kind=8) :: dx, xOffset, yOffset, ownedBox(4)

        integer(kind=4) :: i, j, jm, im
        real(kind=8) :: rho(ni, nj)
        real(kind=8) :: u(ni, nj)
        real(kind=8) :: v(ni, nj)
        real(kind=8) :: T(ni, nj)
        real(kind=8) :: quadWidthX(ni)
        real(kind=8) :: quadWidthY(nj)
        real(kind=8) :: conv, vel2, mass, meanT, hot, cold, middle, dTdx, dTdy, tm, um, vm, cellArea
        real(kind=8) :: w(4), dw(4)
        real(kind=8) :: tmin, tmax, rmin, rmax, xmid, ymid, xlo, xhi, ylo, yhi

        do j = jFirst, jLast
            do i = iFirst, iLast
                cellArea=owned_cell_area(xOffset+(i-0.5d0)*dx, yOffset+(j-0.5d0)*dx, dx, ownedBox)
                if (cellArea<=0d0) cycle
                if (.not.ieee_is_finite(T(i, j)) .or. .not.ieee_is_finite(rho(i, j)) .or. &
                    .not.ieee_is_finite(u(i, j)) .or. .not.ieee_is_finite(v(i, j)) .or. rho(i, j) <= 0.0d0) then
                    write(*, *) 'Invalid state: time, xOffset,yOffset, i,j:', itc, xOffset, yOffset, i, j
                    error stop 'Nonfinite or nonpositive density in owned cells'
                endif
                cellArea = owned_cell_area(xOffset+(i-0.5d0)*dx, yOffset+(j-0.5d0)*dx, dx, ownedBox)
#ifdef SideHeatedCell
                conv = conv+u(i, j)*T(i, j)*cellArea
#else
                conv = conv+v(i, j)*T(i, j)*cellArea
#endif
                vel2 = vel2+(u(i, j)**2+v(i, j)**2)*cellArea
                mass = mass+rho(i, j)*cellArea
                meanT = meanT+T(i, j)*cellArea
                tmin = min(tmin, T(i, j))
                tmax = max(tmax, T(i, j))
                rmin = min(rmin, rho(i, j))
                rmax = max(rmax, rho(i, j))
            enddo
        enddo
        xlo = ownedBox(1)
        xhi = ownedBox(2)
        ylo = ownedBox(3)
        yhi = ownedBox(4)
#ifdef SideHeatedCell
        if (xlo == 0.0d0) then
            do j = jFirst, jLast
                hot = hot+(8.0d0*Thot-9.0d0*T(1, j)+T(2, j))/(3.0d0*dx)*quadWidthY(j)/dble(ny)
            enddo
        endif
        if (xhi == dble(nx)) then
            do j = jFirst, jLast
                cold = cold+(-8.0d0*Tcold+9.0d0*T(ni, j)-T(ni-1, &
                    j))/(3.0d0*dx)*quadWidthY(j)/dble(ny)
            enddo
        endif
        if (xmid >= xlo .and. xmid < xhi) then
            call section_weights((xmid-xOffset)/dx+0.5d0, ni, im, w, dw)
            do j = jFirst, jLast
                if (section_owned_weight(xOffset, yOffset, dx, ownedBox,j,1,xmid)<=0d0) cycle
                if (dx==1.0d0 .and. refineRatio>1) then
                    call fine_section_x(xmid, nint(yOffset)+j, tm, um, dTdx)
                else
                    tm = sum(w*T(im:im+3, j))
                    um = sum(w*u(im:im+3, j))
                    dTdx = sum(dw*T(im:im+3, j))/dx
                endif
                middle = middle+(um*tm/diffusivity-dTdx)*section_owned_weight(xOffset, yOffset, dx, &
            ownedBox,j,1,xmid)/dble(ny)
            enddo
        endif
#else
        if (ylo == 0.0d0) then
            do i = iFirst, iLast
                hot = hot+(8.0d0*Thot-9.0d0*T(i, 1)+T(i, 2))/(3.0d0*dx)*quadWidthX(i)/dble(nx)
            enddo
        endif
        if (yhi == dble(ny)) then
            do i = iFirst, iLast
                cold = cold+(-8.0d0*Tcold+9.0d0*T(i, nj)-T(i, &
                    nj-1))/(3.0d0*dx)*quadWidthX(i)/dble(nx)
            enddo
        endif
        if (ymid >= ylo .and. ymid < yhi) then
            call section_weights((ymid-yOffset)/dx+0.5d0, nj, jm, w, dw)
            do i = iFirst, iLast
                if (section_owned_weight(xOffset, yOffset, dx, ownedBox,i,2,ymid)<=0d0) cycle
                if (dx==1.0d0 .and. refineRatio>1) then
                    call fine_section_y(nint(xOffset)+i, ymid, tm, vm, dTdy)
                else
                    tm = sum(w*T(i, jm:jm+3))
                    vm = sum(w*v(i, jm:jm+3))
                    dTdy = sum(dw*T(i, jm:jm+3))/dx
                endif
                middle = middle+(vm*tm/diffusivity-dTdy)*section_owned_weight(xOffset, yOffset, dx, &
            ownedBox,i,2,ymid)/dble(nx)
            enddo
        endif
#endif
    end subroutine calNuRe_grid


    ! 稳态模式下检查速度和温度相对于上一次检查的变化。
    subroutine check()
#ifdef steadyFlow

        use commondata
        implicit none

        integer(kind=4) :: k, i, j
        real(kind=8) :: du, uu, dt, tt, cellArea

        du = 0.0d0
        uu = 0.0d0
        dt = 0.0d0
        tt = 0.0d0
        call update_host_all(.false.)
        call check_grid(nxCoarse, nyCoarse, dxCoarse, xOffsetCoarse, yOffsetCoarse, iFirstCoarse, &
            iLastCoarse, jFirstCoarse, jLastCoarse, ownedBoxCoarse, du, uu, dt, tt, u_coarse, v_coarse, T_coarse &
#ifdef steadyFlow
            , up_coarse, vp_coarse, Tp_coarse &
#endif
        )
        if (refineRatio > 1) then
            call check_grid(nxLeft, nyLeft, dxFine, xOffsetLeft, yOffsetLeft, iFirstLeft, iLastLeft, &
                jFirstLeft, jLastLeft, ownedBoxLeft, du, uu, dt, tt, u_left, v_left, T_left &
#ifdef steadyFlow
                , up_left, vp_left, Tp_left &
#endif
            )
            call check_grid(nxRight, nyRight, dxFine, xOffsetRight, yOffsetRight, iFirstRight, iLastRight, &
                jFirstRight, jLastRight, ownedBoxRight, du, uu, dt, tt, u_right, v_right, T_right &
#ifdef steadyFlow
                , up_right, vp_right, Tp_right &
#endif
            )
            call check_grid(nxBottom, nyBottom, dxFine, xOffsetBottom, yOffsetBottom, iFirstBottom, &
                iLastBottom, jFirstBottom, jLastBottom, ownedBoxBottom, du, uu, dt, tt, u_bottom, v_bottom, &
            T_bottom &
#ifdef steadyFlow
                , up_bottom, vp_bottom, Tp_bottom &
#endif
            )
            call check_grid(nxTop, nyTop, dxFine, xOffsetTop, yOffsetTop, iFirstTop, iLastTop, jFirstTop, &
                jLastTop, ownedBoxTop, du, uu, dt, tt, u_top, v_top, T_top &
#ifdef steadyFlow
                , up_top, vp_top, Tp_top &
#endif
            )
        endif

        errorU = sqrt(du/max(uu, 1.0d-300))
        errorT = sqrt(dt/max(tt, 1.0d-300))
        open(newunit = k, file = 'Convergence_2DOpenaccMultiblock.dat', status = 'unknown', position = 'append')
        write(k, '(I12,2(1X,ES24.16E3))') itc, errorU, errorT
        close(k)
        write(*, *) 'errorU,errorT:', errorU, errorT
#endif
    end subroutine check


#ifdef steadyFlow

    ! 按实际节点面积累计收敛误差，并保存本次检查时的场。
    subroutine check_grid(ni, nj, dx, xOffset, yOffset, iFirst, iLast, jFirst, jLast, ownedBox, du, uu, dt, &
            tt, u, v, T &
#ifdef steadyFlow
            , up, vp, Tp &
#endif
        )

        use commondata, only: owned_cell_area
        implicit none
        integer(kind=4) :: ni, nj, iFirst, iLast, jFirst, jLast
        real(kind=8) :: dx, xOffset, yOffset, ownedBox(4)

        integer(kind=4) :: i, j
        real(kind=8) :: u(ni, nj)
        real(kind=8) :: v(ni, nj)
        real(kind=8) :: T(ni, nj)
#ifdef steadyFlow
        real(kind=8) :: up(ni, nj)
        real(kind=8) :: vp(ni, nj)
        real(kind=8) :: Tp(ni, nj)
#endif
        real(kind=8) :: du, uu, dt, tt, cellArea

        do j = jFirst, jLast
            do i = iFirst, iLast
                cellArea = owned_cell_area(xOffset+(i-0.5d0)*dx, yOffset+(j-0.5d0)*dx, dx, ownedBox)
                du = du+cellArea*((u(i, j)-up(i, j))**2+(v(i, j)-vp(i, j))**2)
                uu = uu+cellArea*(u(i, j)**2+v(i, j)**2)
                dt = dt+cellArea*(T(i, j)-Tp(i, j))**2
                tt = tt+cellArea*T(i, j)**2
            enddo
        enddo
        up = u
        vp = v
        Tp = T
    end subroutine check_grid


#endif

    ! 输出具名 coarse/left/right/bottom/top 区域及节点积分面积。
    subroutine output_Tecplot()

        use commondata
        implicit none

        integer(kind=4) :: k, i, j
        character(16) :: num

        pltFileNum = pltFileNum+1
        write(num, '(I10.10)') pltFileNum
        call update_host_all(.false.)
        open(newunit = k, file = pltFolderPrefix//'-'//trim(num)//'.dat', status = 'replace')
        write(k, '(a)') 'VARIABLES="x/L","y/L","u","v","T","rho","dx/L","integration_area/L^2"'
        call output_Tecplot_grid(nxCoarse, nyCoarse, dxCoarse, xOffsetCoarse, yOffsetCoarse, iFirstCoarse, &
            iLastCoarse, jFirstCoarse, jLastCoarse, ownedBoxCoarse, 'coarse', k, rho_coarse, u_coarse, &
            v_coarse, T_coarse)
        if (refineRatio > 1) then
            call output_Tecplot_grid(nxLeft, nyLeft, dxFine, xOffsetLeft, yOffsetLeft, iFirstLeft, iLastLeft, &
                jFirstLeft, jLastLeft, ownedBoxLeft, 'left', k, rho_left, u_left, v_left, T_left)
            call output_Tecplot_grid(nxRight, nyRight, dxFine, xOffsetRight, yOffsetRight, iFirstRight, &
                iLastRight, jFirstRight, jLastRight, ownedBoxRight, 'right', k, rho_right, u_right, &
            v_right, T_right)
            call output_Tecplot_grid(nxBottom, nyBottom, dxFine, xOffsetBottom, yOffsetBottom, iFirstBottom, &
                iLastBottom, jFirstBottom, jLastBottom, ownedBoxBottom, 'bottom', k, rho_bottom, u_bottom, &
                v_bottom, T_bottom)
            call output_Tecplot_grid(nxTop, nyTop, dxFine, xOffsetTop, yOffsetTop, iFirstTop, iLastTop, &
                jFirstTop, jLastTop, ownedBoxTop, 'top', k, rho_top, u_top, v_top, T_top)
        endif

        close(k)
    end subroutine output_Tecplot


    subroutine output_Tecplot_grid(ni, nj, dx, xOffset, yOffset, iFirst, iLast, jFirst, jLast, ownedBox, &
            gridName, k, rho, u, v, T)

        use commondata, only: lengthUnit, timeUnit, itc, owned_cell_area
        implicit none
        integer(kind=4) :: ni, nj, iFirst, iLast, jFirst, jLast
        real(kind=8) :: dx, xOffset, yOffset, ownedBox(4)
        character(*), intent(in) :: gridName

        integer(kind=4) :: k, i, j
        real(kind=8) :: rho(ni, nj)
        real(kind=8) :: u(ni, nj)
        real(kind=8) :: v(ni, nj)
        real(kind=8) :: T(ni, nj)

        write(k, '(a,a,a,I0,a,I0,a,ES24.16E3)') 'ZONE T="', trim(gridName), '", I=', iLast-iFirst+1, &
            ', J=', jLast-jFirst+1, ', F=POINT, SOLUTIONTIME=', dble(itc)/timeUnit
        do j = jFirst, jLast
            do i = iFirst, iLast
                write(k, '(8(ES24.16E3,1X))') ((xOffset+(dble(i)-0.5d0)*dx))/lengthUnit, &
                    ((yOffset+(dble(j)-0.5d0)*dx))/lengthUnit, u(i, j), v(i, j), T(i, j), rho(i, j), &
                    dx/lengthUnit, &
                    owned_cell_area(xOffset+(i-0.5d0)*dx, yOffset+(j-0.5d0)*dx, dx, ownedBox)/lengthUnit**2
            enddo
        enddo
    end subroutine output_Tecplot_grid


    ! 输出 v3 快照；区域数量为 1 或 5，保留显式二维积分面积。
    subroutine output_SnapshotFile()

        use commondata
        implicit none

        integer(kind=4) :: k, i,j
        character(16) :: num

        snapshotFileNum = snapshotFileNum+1
        write(num, '(I10.10)') snapshotFileNum
        call update_host_all(.false.)
        open(newunit = k, file = snapshotFilePrefix//'-'//trim(num)//'.bin', form = 'unformatted', &
            access = 'stream', status = 'replace')
        write(k) snapshotMagic, merge(5, 1, refineRatio > 1), nx, ny, itc, dble(itc)/timeUnit, lengthUnit
        call output_SnapshotFile_grid(nxCoarse, nyCoarse, dxCoarse, xOffsetCoarse, yOffsetCoarse, &
            iFirstCoarse, iLastCoarse, jFirstCoarse, jLastCoarse, ownedBoxCoarse, k, rho_coarse, u_coarse, &
            v_coarse, T_coarse, quadWidthX_coarse, quadWidthY_coarse)
        if (refineRatio > 1) then
            call output_SnapshotFile_grid(nxLeft, nyLeft, dxFine, xOffsetLeft, yOffsetLeft, iFirstLeft, &
                iLastLeft, jFirstLeft, jLastLeft, ownedBoxLeft, k, rho_left, u_left, v_left, T_left, &
                quadWidthX_left, quadWidthY_left)
            call output_SnapshotFile_grid(nxRight, nyRight, dxFine, xOffsetRight, yOffsetRight, iFirstRight, &
                iLastRight, jFirstRight, jLastRight, ownedBoxRight, k, rho_right, u_right, v_right, T_right, &
                quadWidthX_right, quadWidthY_right)
            call output_SnapshotFile_grid(nxBottom, nyBottom, dxFine, xOffsetBottom, yOffsetBottom, &
                iFirstBottom, iLastBottom, jFirstBottom, jLastBottom, ownedBoxBottom, k, rho_bottom, u_bottom, &
                v_bottom, T_bottom, quadWidthX_bottom, quadWidthY_bottom)
            call output_SnapshotFile_grid(nxTop, nyTop, dxFine, xOffsetTop, yOffsetTop, iFirstTop, iLastTop, &
                jFirstTop, jLastTop, ownedBoxTop, k, rho_top, u_top, v_top, T_top, quadWidthX_top, quadWidthY_top)
        endif

        close(k)
    end subroutine output_SnapshotFile


    subroutine output_SnapshotFile_grid(ni, nj, dx, xOffset, yOffset, iFirst, iLast, jFirst, jLast, &
            ownedBox, k, rho, u, v, T, quadWidthX, quadWidthY)

        use commondata, only: owned_cell_area
        implicit none
        integer(kind=4) :: ni, nj, iFirst, iLast, jFirst, jLast
        real(kind=8) :: dx, xOffset, yOffset, ownedBox(4)

        integer(kind=4) :: k, i, j
        real(kind=8) :: rho(ni, nj)
        real(kind=8) :: u(ni, nj)
        real(kind=8) :: v(ni, nj)
        real(kind=8) :: T(ni, nj)
        real(kind=8) :: quadWidthX(ni)
        real(kind=8) :: quadWidthY(nj)

        write(k) iLast-iFirst+1, jLast-jFirst+1, &
            xOffset+(dble(iFirst)-0.5d0)*dx, yOffset+(dble(jFirst)-0.5d0)*dx, &
            dx, ownedBox
        write(k) quadWidthX(iFirst:iLast), quadWidthY(jFirst:jLast)
        ! Snapshot v3 appends explicit 2D area after the separable coordinate weights.
        write(k) ((owned_cell_area(xOffset+(i-0.5d0)*dx, yOffset+(j-0.5d0)*dx, dx, &
            ownedBox),i=iFirst,iLast),j=jFirst,jLast)
        write(k) u(iFirst:iLast, jFirst:jLast), v(iFirst:iLast, jFirst:jLast), &
            T(iFirst:iLast, jFirst:jLast), rho(iFirst:iLast, jFirst:jLast)
    end subroutine output_SnapshotFile_grid


    ! 保存 v11 完整状态、时间历史和输出时钟，再更新 latest.meta。
    subroutine output_ReloadFile()

        use commondata
        implicit none

        integer(kind=4) :: k, currentModel(12)
        real(kind=8) :: currentPhysics(16)
        character(16) :: num
        character(256) :: name

        call update_host_all(.true.)
#ifdef steadyFlow
        reloadFileNum = itc
#else
        reloadFileNum = reloadFileNum+1
#endif
        write(num, '(I12.12)') reloadFileNum
        name = reloadFilePrefix//'-'//trim(num)//'.bin'
        open(newunit = k, file = trim(name), access = 'stream', form = 'unformatted', status = 'replace')
        write(k) restartMagic, nx, ny, refineRatio, fineLayerCellsLeft, fineLayerCellsRight, &
            fineLayerCellsBottom, fineLayerCellsTop, overlapCells, merge(5, 1, refineRatio > 1)
        call model_signature(currentModel)
        call physical_signature(currentPhysics)
        write(k) currentModel, currentPhysics
        write(k) itc, nextSample, nextReload, nextPlt, snapshotFileNum, pltFileNum, reloadFileNum, errorU, errorT
        call output_ReloadFile_grid(nxCoarse, nyCoarse, historyLastCoarse, dxCoarse, xOffsetCoarse, &
            yOffsetCoarse, iFirstCoarse, iLastCoarse, jFirstCoarse, jLastCoarse, ownedBoxCoarse, k, &
            f_coarse, g_coarse, u_coarse, v_coarse, T_coarse, rho_coarse, Fx_coarse, Fy_coarse, &
            rhoHistory_coarse, uHistory_coarse, vHistory_coarse, THistory_coarse, FxHistory_coarse, &
            FyHistory_coarse, flowNeqHistory_coarse, thermalNeqHistory_coarse &
#ifdef steadyFlow
            , up_coarse, vp_coarse, Tp_coarse &
#endif
        )
        if (refineRatio > 1) then
            call output_ReloadFile_grid(nxLeft, nyLeft, 0, dxFine, xOffsetLeft, yOffsetLeft, iFirstLeft, &
                iLastLeft, jFirstLeft, jLastLeft, ownedBoxLeft, k, f_left, g_left, u_left, v_left, T_left, &
                rho_left, Fx_left, Fy_left, rhoHistory_left, uHistory_left, vHistory_left, THistory_left, &
                FxHistory_left, FyHistory_left, flowNeqHistory_left, thermalNeqHistory_left &
#ifdef steadyFlow
                , up_left, vp_left, Tp_left &
#endif
            )
            call output_ReloadFile_grid(nxRight, nyRight, 0, dxFine, xOffsetRight, yOffsetRight, iFirstRight, &
                iLastRight, jFirstRight, jLastRight, ownedBoxRight, k, f_right, g_right, u_right, v_right, &
                T_right, rho_right, Fx_right, Fy_right, rhoHistory_right, uHistory_right, vHistory_right, &
                THistory_right, FxHistory_right, FyHistory_right, flowNeqHistory_right, thermalNeqHistory_right &
#ifdef steadyFlow
                , up_right, vp_right, Tp_right &
#endif
            )
            call output_ReloadFile_grid(nxBottom, nyBottom, 0, dxFine, xOffsetBottom, yOffsetBottom, &
                iFirstBottom, iLastBottom, jFirstBottom, jLastBottom, ownedBoxBottom, k, f_bottom, g_bottom, &
                u_bottom, v_bottom, T_bottom, rho_bottom, Fx_bottom, Fy_bottom, rhoHistory_bottom, &
                uHistory_bottom, vHistory_bottom, THistory_bottom, FxHistory_bottom, FyHistory_bottom, &
                flowNeqHistory_bottom, thermalNeqHistory_bottom &
#ifdef steadyFlow
                , up_bottom, vp_bottom, Tp_bottom &
#endif
            )
            call output_ReloadFile_grid(nxTop, nyTop, 0, dxFine, xOffsetTop, yOffsetTop, iFirstTop, iLastTop, &
                jFirstTop, jLastTop, ownedBoxTop, k, f_top, g_top, u_top, v_top, T_top, rho_top, Fx_top, &
                Fy_top, rhoHistory_top, uHistory_top, vHistory_top, THistory_top, FxHistory_top, FyHistory_top, &
                flowNeqHistory_top, thermalNeqHistory_top &
#ifdef steadyFlow
                , up_top, vp_top, Tp_top &
#endif
            )
        endif

        close(k)
        ! 完整写完独立编号的 checkpoint 后再更新 latest 指针，旧 checkpoint 仍保留。
        open(newunit = k, file = reloadFilePrefix//'-latest.meta', status = 'replace')
        write(k, '(a)') trim(name)
        close(k)
    end subroutine output_ReloadFile


    subroutine output_ReloadFile_grid(ni, nj, nh, dx, xOffset, yOffset, iFirst, iLast, jFirst, jLast, &
            ownedBox, k, f, g, u, v, T, rho, Fx, Fy, rhoHistory, uHistory, vHistory, &
        THistory, FxHistory, FyHistory, flowNeqHistory, thermalNeqHistory &
#ifdef steadyFlow
            , up, vp, Tp &
#endif
        )

        implicit none
        integer(kind=4) :: ni, nj, nh, iFirst, iLast, jFirst, jLast
        real(kind=8) :: dx, xOffset, yOffset, ownedBox(4)

        integer(kind=4) :: historyIndex

        integer(kind=4) :: k
        real(kind=8) :: f(ni, nj, 0:8)
        real(kind=8) :: g(ni, nj, 0:4)
        real(kind=8) :: u(ni, nj)
        real(kind=8) :: v(ni, nj)
        real(kind=8) :: T(ni, nj)
        real(kind=8) :: rho(ni, nj)
        real(kind=8) :: Fx(ni, nj)
        real(kind=8) :: Fy(ni, nj)
        real(kind=8) :: rhoHistory(ni, nj, 0:nh)
        real(kind=8) :: uHistory(ni, nj, 0:nh)
        real(kind=8) :: vHistory(ni, nj, 0:nh)
        real(kind=8) :: THistory(ni, nj, 0:nh)
        real(kind=8) :: FxHistory(ni, nj, 0:nh)
        real(kind=8) :: FyHistory(ni, nj, 0:nh)
        real(kind=8) :: flowNeqHistory(ni, nj, 0:8, 0:nh)
        real(kind=8) :: thermalNeqHistory(ni, nj, 0:4, 0:nh)
#ifdef steadyFlow
        real(kind=8) :: up(ni, nj)
        real(kind=8) :: vp(ni, nj)
        real(kind=8) :: Tp(ni, nj)
#endif

        write(k) ni, nj, iFirst, iLast, jFirst, jLast, nh, &
            xOffset, yOffset, dx, ownedBox
        ! 按 v11 的“每时间层：宏观量、归一化力、流场矩、温度矩”顺序写入，不另建打包数组。
        write(k) f, g, u, v, T, rho, Fx, Fy, &
            (rhoHistory(:, :, historyIndex), uHistory(:, :, historyIndex), vHistory(:, :, historyIndex), &
                THistory(:, :, historyIndex), FxHistory(:, :, historyIndex), FyHistory(:, :, historyIndex), &
                flowNeqHistory(:, :, :, historyIndex), thermalNeqHistory(:, :, :, historyIndex), historyIndex=0, &
                nh)
#ifdef steadyFlow
        write(k) up, vp, Tp
#endif
    end subroutine output_ReloadFile_grid


    ! 核对 v11 格式、网格、物理参数及输出间隔，再恢复各套数组。
    subroutine read_restart()

        use commondata
        implicit none

        integer(kind=4) :: k, ios, head(9), geom(7), sig(12), currentModel(12)
        real(kind=8) :: phys(16), coord(7), currentPhysics(16)
        character(16) :: magic
        character(16) :: num
        character(256) :: name
        logical :: metaExists

        inquire(file = reloadFilePrefix//'-latest.meta', exist = metaExists)
        if (metaExists) then
            open(newunit = k, file = reloadFilePrefix//'-latest.meta', status = 'old', iostat = ios)
            if (ios /= 0) error stop 'Cannot open multiblock latest.meta'
            read(k, '(a)', iostat = ios) name
            close(k)
            if (ios /= 0 .or. len_trim(name) == 0) error stop 'Invalid latest.meta'
        else
            if (reloadFileNum <= 0) error stop 'Missing latest.meta: set reloadFileNum to a saved checkpoint number'
            write(num, '(I12.12)') reloadFileNum
            name = reloadFilePrefix//'-'//trim(num)//'.bin'
        endif
        open(newunit = k, file = trim(name), status = 'old', access = 'stream', form = 'unformatted', iostat = ios)
        if (ios /= 0) error stop 'Missing multiblock checkpoint'
        read(k, iostat = ios) magic, head
        if (ios /= 0) error stop 'Truncated multiblock checkpoint header'
        if (magic /= restartMagic) error stop 'Expected compact v11 checkpoint; convert older ring checkpoints first'
        if (any(head /= [nx, ny, refineRatio, fineLayerCellsLeft, fineLayerCellsRight, &
            fineLayerCellsBottom, fineLayerCellsTop, overlapCells, merge(5, 1, refineRatio > 1)])) &
            error stop 'Restart mesh/refinement mismatch'
        read(k) sig, phys
        call model_signature(currentModel)
        call physical_signature(currentPhysics)
        if (any(sig /= currentModel) .or. any(phys /= currentPhysics)) &
            error stop 'Restart model/physics/output-cadence mismatch'
        read(k) itc, nextSample, nextReload, nextPlt, snapshotFileNum, pltFileNum, reloadFileNum, errorU, errorT
        if (itc < 0 .or. mod(itc, refineRatio) /= 0) error stop 'Restart is not at a synchronized time'
        call read_restart_grid(nxCoarse, nyCoarse, historyLastCoarse, dxCoarse, xOffsetCoarse, &
            yOffsetCoarse, iFirstCoarse, iLastCoarse, jFirstCoarse, jLastCoarse, ownedBoxCoarse, k, &
            f_coarse, g_coarse, u_coarse, v_coarse, T_coarse, rho_coarse, Fx_coarse, Fy_coarse, &
            rhoHistory_coarse, uHistory_coarse, vHistory_coarse, THistory_coarse, FxHistory_coarse, &
            FyHistory_coarse, flowNeqHistory_coarse, thermalNeqHistory_coarse &
#ifdef steadyFlow
            , up_coarse, vp_coarse, Tp_coarse &
#endif
        )
        if (refineRatio > 1) then
            call read_restart_grid(nxLeft, nyLeft, 0, dxFine, xOffsetLeft, yOffsetLeft, iFirstLeft, iLastLeft, &
                jFirstLeft, jLastLeft, ownedBoxLeft, k, f_left, g_left, u_left, v_left, T_left, rho_left, &
                Fx_left, Fy_left, rhoHistory_left, uHistory_left, vHistory_left, THistory_left, FxHistory_left, &
                FyHistory_left, flowNeqHistory_left, thermalNeqHistory_left &
#ifdef steadyFlow
                , up_left, vp_left, Tp_left &
#endif
            )
            call read_restart_grid(nxRight, nyRight, 0, dxFine, xOffsetRight, yOffsetRight, iFirstRight, &
                iLastRight, jFirstRight, jLastRight, ownedBoxRight, k, f_right, g_right, u_right, v_right, &
                T_right, rho_right, Fx_right, Fy_right, rhoHistory_right, uHistory_right, vHistory_right, &
                THistory_right, FxHistory_right, FyHistory_right, flowNeqHistory_right, thermalNeqHistory_right &
#ifdef steadyFlow
                , up_right, vp_right, Tp_right &
#endif
            )
            call read_restart_grid(nxBottom, nyBottom, 0, dxFine, xOffsetBottom, yOffsetBottom, iFirstBottom, &
                iLastBottom, jFirstBottom, jLastBottom, ownedBoxBottom, k, f_bottom, g_bottom, u_bottom, &
                v_bottom, T_bottom, rho_bottom, Fx_bottom, Fy_bottom, rhoHistory_bottom, uHistory_bottom, &
                vHistory_bottom, THistory_bottom, FxHistory_bottom, FyHistory_bottom, flowNeqHistory_bottom, &
                thermalNeqHistory_bottom &
#ifdef steadyFlow
                , up_bottom, vp_bottom, Tp_bottom &
#endif
            )
            call read_restart_grid(nxTop, nyTop, 0, dxFine, xOffsetTop, yOffsetTop, iFirstTop, iLastTop, &
                jFirstTop, jLastTop, ownedBoxTop, k, f_top, g_top, u_top, v_top, T_top, rho_top, Fx_top, &
                Fy_top, rhoHistory_top, uHistory_top, vHistory_top, THistory_top, FxHistory_top, FyHistory_top, &
                flowNeqHistory_top, thermalNeqHistory_top &
#ifdef steadyFlow
                , up_top, vp_top, Tp_top &
#endif
            )
        endif

        close(k)
    end subroutine read_restart


    subroutine read_restart_grid(ni, nj, nh, dx, xOffset, yOffset, iFirst, iLast, jFirst, jLast, ownedBox, &
            k, f, g, u, v, T, rho, Fx, Fy, rhoHistory, uHistory, vHistory, THistory, &
        FxHistory, FyHistory, flowNeqHistory, thermalNeqHistory &
#ifdef steadyFlow
            , up, vp, Tp &
#endif
        )

        implicit none
        integer(kind=4) :: ni, nj, nh, iFirst, iLast, jFirst, jLast
        real(kind=8) :: dx, xOffset, yOffset, ownedBox(4)

        integer(kind=4) :: historyIndex

        integer(kind=4) :: k, ios, geom(7)
        real(kind=8) :: f(ni, nj, 0:8)
        real(kind=8) :: g(ni, nj, 0:4)
        real(kind=8) :: u(ni, nj)
        real(kind=8) :: v(ni, nj)
        real(kind=8) :: T(ni, nj)
        real(kind=8) :: rho(ni, nj)
        real(kind=8) :: Fx(ni, nj)
        real(kind=8) :: Fy(ni, nj)
        real(kind=8) :: rhoHistory(ni, nj, 0:nh)
        real(kind=8) :: uHistory(ni, nj, 0:nh)
        real(kind=8) :: vHistory(ni, nj, 0:nh)
        real(kind=8) :: THistory(ni, nj, 0:nh)
        real(kind=8) :: FxHistory(ni, nj, 0:nh)
        real(kind=8) :: FyHistory(ni, nj, 0:nh)
        real(kind=8) :: flowNeqHistory(ni, nj, 0:8, 0:nh)
        real(kind=8) :: thermalNeqHistory(ni, nj, 0:4, 0:nh)
#ifdef steadyFlow
        real(kind=8) :: up(ni, nj)
        real(kind=8) :: vp(ni, nj)
        real(kind=8) :: Tp(ni, nj)
#endif
        real(kind=8) :: coord(7)
        logical :: metaExists

        read(k) geom, coord
        if (any(geom /= [ni, nj, iFirst, iLast, jFirst, jLast, nh]) .or. &
            any(coord /= [xOffset, yOffset, dx, ownedBox])) error stop 'Restart block layout mismatch'
        read(k, iostat = ios) f, g, u, v, T, rho, Fx, Fy, &
            (rhoHistory(:, :, historyIndex), uHistory(:, :, historyIndex), vHistory(:, :, historyIndex), &
                THistory(:, :, historyIndex), FxHistory(:, :, historyIndex), FyHistory(:, :, historyIndex), &
                flowNeqHistory(:, :, :, historyIndex), thermalNeqHistory(:, :, :, historyIndex), historyIndex=0, &
                nh)
        if (ios /= 0) error stop 'Incomplete checkpoint state/history'
#ifdef steadyFlow
        read(k) up, vp, Tp
#endif
    end subroutine read_restart_grid


    subroutine model_signature(sig)

        implicit none

        integer(kind=4), intent(out) :: sig(12)

        sig = 0
#ifdef steadyFlow
        sig(1) = 1
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
    end subroutine model_signature


    subroutine physical_signature(sig)

        use commondata, only: Rayleigh, Prandtl, Mach, Thot, Tcold, Snu, Sq, Qk, Qnu, thermalA, &
            outputSnapshotInterval, reloadFileInterval, outputPltFileInterval
#ifdef SideHeatedHa
        use commondata, only: Ha, phi
#endif
        implicit none

        real(kind=8), intent(out) :: sig(16)

        sig = [Rayleigh, Prandtl, Mach, Thot, Tcold, Snu, Sq, Qk, Qnu, thermalA, outputSnapshotInterval, &
            reloadFileInterval, outputPltFileInterval, 0.0d0, 0.0d0, 0.0d0]
#ifdef SideHeatedHa
        sig(14) = Ha
        sig(15) = phi
#endif
    end subroutine physical_signature


    subroutine check_history()

        use commondata, only: timeUnit, outputSnapshotInterval, NuReHistoryFile, nextSample, scheduled_step, &
            ieee_is_finite
        implicit none

        integer(kind=4) :: k, ios, n
        real(kind=8) :: values(12), lastTime, expected
        character(512) :: line

        lastTime = -1.0d0
        n = 0
        open(newunit = k, file = NuReHistoryFile, status = 'old', iostat = ios)
        if (ios /= 0) error stop 'Restart requires the matching Nu/Re history'
        do
            read(k, '(a)', iostat = ios) line
            if (ios < 0) exit
            if (ios > 0) error stop 'Nu/Re history read error'
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            read(line, *, iostat = ios) values
            if (ios /= 0 .or. .not.all(ieee_is_finite(values))) error stop 'Invalid Nu/Re history row'
            if (values(1) <= lastTime) error stop 'Nonmonotone Nu/Re history'
            n = n+1
            lastTime = values(1)
        enddo
        close(k)
        if (n /= nextSample-1) error stop 'History sample count differs from checkpoint; use the matching history'
        if (n > 0) then
            expected = dble(scheduled_step(n, outputSnapshotInterval))/timeUnit
            if (abs(lastTime-expected) > 1.0d-10*max(1.0d0, &
            expected)) error stop 'Checkpoint/history time mismatch'
        endif
    end subroutine check_history


    subroutine average_window(t0, t1, result, coverage)

        use commondata, only: NuReHistoryFile, ieee_is_finite
        implicit none

        real(kind=8), intent(in) :: t0, t1
        real(kind=8), intent(out) :: result(5), coverage
        integer(kind=4) :: k, ios
        real(kind=8) :: prev(12), row(12), aa, bb, dt, ra(5), rb(5), left(5), right(5)
        logical :: hasPrev
        character(512) :: line

        result = 0.0d0
        coverage = 0.0d0
        hasPrev = .false.
        open(newunit = k, file = NuReHistoryFile, status = 'old')
        do
            read(k, '(a)', iostat = ios) line
            if (ios < 0) exit
            if (ios > 0) error stop 'Nu/Re history read failure'
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            read(line, *, iostat = ios) row
            if (ios /= 0 .or. .not.all(ieee_is_finite(row))) error stop 'Invalid history during averaging'
            if (hasPrev) then
                if (row(1) <= prev(1)) error stop 'Nonmonotone statistics time'
                aa = max(t0, prev(1))
                bb = min(t1, row(1))
                dt = row(1)-prev(1)
                if (bb > aa) then
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
        if (coverage > 0.0d0) then
            result = result/coverage
            result(2) = sqrt(max(0.0d0, result(2)))
        endif
    end subroutine average_window


    subroutine output_unsteady_NuRe_postprocess()
#ifdef unsteadyFlow

        use commondata, only: unsteadyAverageStartTf, unsteadyAverageEndTf, unsteadyAverageMidTf
        implicit none

        integer(kind=4) :: k
        real(kind=8) :: allMean(5), firstMean(5), lastMean(5), c0, c1, c2, relative(5)

        call average_window(unsteadyAverageStartTf, unsteadyAverageEndTf, allMean, c0)
        call average_window(unsteadyAverageStartTf, unsteadyAverageMidTf, firstMean, c1)
        call average_window(unsteadyAverageMidTf, unsteadyAverageEndTf, lastMean, c2)
        open(newunit = k, file = 'NuReStatistics_2DOpenaccMultiblock.dat', status = 'replace')
        write(k, *) 'Window t_ff:', unsteadyAverageStartTf, unsteadyAverageEndTf
        write(k, *) 'Covered durations, full/first/last:', c0, c1, c2
        if (abs(c0-(unsteadyAverageEndTf-unsteadyAverageStartTf)) > 1.0d-8 .or. c1 <= 0.0d0 .or. c2 <= 0.0d0) then
            write(k, *) 'INCOMPLETE: requested statistics window is not fully covered; no final mean is reported.'
        else
            relative = abs(lastMean-firstMean)/max(abs(allMean), 1.0d-30)
            write(k, *) 'Columns: NuVolAvg Re_rms_space_time Nu_hot Nu_cold Nu_middle'
            write(k, '(a,5ES24.16E3)') 'whole:', allMean
            write(k, '(a,5ES24.16E3)') 'first:', firstMean
            write(k, '(a,5ES24.16E3)') 'last :', lastMean
            write(k, '(a,5ES24.16E3)') 'relative half-window difference:', relative
        endif
        close(k)
#endif
    end subroutine output_unsteady_NuRe_postprocess


    integer(kind=4) function scheduled_step(index, interval) result(step)

        use commondata, only: refineRatio, timeUnit
        implicit none

        integer(kind=4), intent(in) :: index
        real(kind=8), intent(in) :: interval

        step = max(refineRatio, ceiling(dble(index)*interval*timeUnit/dble(refineRatio))*refineRatio)
    end function scheduled_step


    ! 节点控制面积按真实积分分区裁剪。细区再扣除中心区域，重叠计算不重复积分。
    real(8) function owned_cell_area(x,y,dx,ownedBox) result(area)
        use commondata, only: centerBox, refineRatio
        implicit none
        real(8), intent(in) :: x,y,dx,ownedBox(4)
        real(8) :: wx,wy
        wx=max(0d0,min(x+dx/2,ownedBox(2))-max(x-dx/2,ownedBox(1)))
        wy=max(0d0,min(y+dx/2,ownedBox(4))-max(y-dx/2,ownedBox(3)))
        area=wx*wy
        if (refineRatio>1 .and. dx==1d0) then
            wx=max(0d0,min(x+dx/2,ownedBox(2),centerBox(2))-max(x-dx/2,ownedBox(1),centerBox(1)))
            wy=max(0d0,min(y+dx/2,ownedBox(4),centerBox(4))-max(y-dx/2,ownedBox(3),centerBox(3)))
            area=area-wx*wy
        endif
    end function owned_cell_area

    real(8) function section_owned_weight(xOffset,yOffset,dx,ownedBox,k,axis,position) result(w)
        use commondata, only: centerBox, refineRatio
        implicit none
        real(8), intent(in) :: xOffset,yOffset,dx,ownedBox(4),position
        integer, intent(in) :: k,axis
        real(8) :: q,lo,hi,clo,chi
        if (axis==1) then
            q=yOffset+(k-.5d0)*dx; lo=ownedBox(3); hi=ownedBox(4)
            clo=centerBox(3); chi=centerBox(4)
        else
            q=xOffset+(k-.5d0)*dx; lo=ownedBox(1); hi=ownedBox(2)
            clo=centerBox(1); chi=centerBox(2)
        endif
        w=max(0d0,min(q+dx/2,hi)-max(q-dx/2,lo))
        if (refineRatio>1 .and. dx==1d0) then
            if (axis==1) then
                if(position<centerBox(1) .or. position>=centerBox(2)) return
            else
                if(position<centerBox(3) .or. position>=centerBox(4)) return
            endif
            w=w-max(0d0,min(q+dx/2,hi,chi)-max(q-dx/2,lo,clo))
        endif
    end function section_owned_weight

    ! 中线诊断沿用全域细网格的四点模板，即使模板跨越两个细数组也不改成单侧模板。
    subroutine fine_section_x(x,j,tm,um,dTdx)
        use commondata
        implicit none
        real(8), intent(in) :: x
        integer, intent(in) :: j
        real(8), intent(out) :: tm,um,dTdx
        real(8) :: w(4),dw(4),temp(4),velocity(4),fine_value
        integer :: first,a
        external :: fine_value
        call section_weights(x+0.5d0,nx,first,w,dw)
        do a=1,4
            temp(a)=fine_value(first+a-1,j,T_left,T_right,T_bottom,T_top)
            velocity(a)=fine_value(first+a-1,j,u_left,u_right,u_bottom,u_top)
        enddo
        tm=sum(w*temp); um=sum(w*velocity); dTdx=sum(dw*temp)
    end subroutine fine_section_x

    subroutine fine_section_y(i,y,tm,vm,dTdy)
        use commondata
        implicit none
        integer, intent(in) :: i
        real(8), intent(in) :: y
        real(8), intent(out) :: tm,vm,dTdy
        real(8) :: w(4),dw(4),temp(4),velocity(4),fine_value
        integer :: first,a
        external :: fine_value
        call section_weights(y+0.5d0,ny,first,w,dw)
        do a=1,4
            temp(a)=fine_value(i,first+a-1,T_left,T_right,T_bottom,T_top)
            velocity(a)=fine_value(i,first+a-1,v_left,v_right,v_bottom,v_top)
        enddo
        tm=sum(w*temp); vm=sum(w*velocity); dTdy=sum(dw*temp)
    end subroutine fine_section_y

    ! 仅用于主机端中线积分。全局细节点编号减去偏移量，就是相应数组的本地下标。
    real(8) function fine_value(i,j,left,right,bottom,top) result(value)
        use commondata, only: nxLeft,nyLeft,nxRight,nyRight,nxBottom,nyBottom,nxTop,nyTop, &
            xOffsetRight,xOffsetBottom,xOffsetTop,yOffsetTop
        implicit none
        integer, intent(in) :: i,j
        real(8), intent(in) :: left(nxLeft,nyLeft),right(nxRight,nyRight), &
            bottom(nxBottom,nyBottom),top(nxTop,nyTop)
        if(i<=nxLeft) then
            value=left(i,j)
        else if(i>nint(xOffsetRight)) then
            value=right(i-nint(xOffsetRight),j)
        else if(j<=nyBottom) then
            value=bottom(i-nint(xOffsetBottom),j)
        else if(j>nint(yOffsetTop)) then
            value=top(i-nint(xOffsetTop),j-nint(yOffsetTop))
        else
            error stop 'Fine section requested an inactive central node'
        endif
    end function fine_value
