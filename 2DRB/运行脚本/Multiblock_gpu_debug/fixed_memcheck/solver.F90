!=============================================================
!!!    注释区，代码描述
!!!    二维浮力驱动自然对流 OpenACC 静态多块网格版本
!!!    块内 D2Q9/D2Q5 算法保持不变
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

! 温度算法：保留原 D2Q5 MRT 旧算法，平衡矩系数为 paraA，不使用热流历史修正。
#define EnableLegacyThermalScheme

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
        integer(kind=4), parameter :: nx = 256, ny = 256    ! 手动设置最细网格等效分辨率，修改后重新编译
        integer(kind=4), parameter :: refineRatio=2    ! 粗/细格距及时间步之比；逐级二分取2、4、8等2的幂，1为单块
        integer(kind=4), parameter :: fineLayerCellsLeft = 32      ! 从左墙向内数的细节点编号，x = Left-0.5
        integer(kind=4), parameter :: fineLayerCellsRight = 33     ! 从右墙向内数的细节点编号，x = nx-Right+0.5
        integer(kind=4), parameter :: fineLayerCellsBottom = 32    ! 从下墙向内数的细节点编号，y = Bottom-0.5
        integer(kind=4), parameter :: fineLayerCellsTop = 33       ! 从上墙向内数的细节点编号，y = ny-Top+0.5
        ! 默认交界面为 x/y=127.5、895.5；粗细块共用这些基准边界，积分也在此分区。
        ! 中心宽度 nx-Left-Right+1、高度 ny-Bottom-Top+1 须为正，且均能被 refineRatio 整除。
        ! 右=左+1、上=下+1 时，中心宽高为 nx-2*Left、ny-2*Bottom；较大粗细比仍须检查整除。
        ! 默认 128/129 对 nx=ny=1024、粗细比 2/4/8 都满足条件；不再自动移动右/上交界面。
        ! 多块模式 nx、ny 仍须整除粗细比；四个细网格编号本身无须整除，但距墙均至少 (overlapCells+1)*refineRatio。
        integer(kind=4), parameter :: overlapCells = 2    ! 每块向交界面外延伸的粗格距数；总重叠跨度为两倍
        integer(kind=4), parameter :: interfaceSkin = 2    ! 原流场->温度场分步推进需要两层人工边界；这些节点在接收后参与碰撞

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
        real(kind=8), parameter :: Rayleigh = 1.0d6    ! 手动设置瑞利数，修改后重新编译
        real(kind=8), parameter :: Prandtl = 0.71d0
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
        integer(kind=4), parameter :: outputSnapshotFile = 0    ! 0: 不输出快照；1: 输出，Nu/Re 仍独立采样
        integer(kind=4), parameter :: outputPltFile = 0         ! 0: 不输出 Tecplot；1: 输出
        integer(kind=4), parameter :: outputReloadFile = 0      ! 0: 不输出重启文件；1: 输出
        integer(kind=4), parameter :: itc_max = 2        ! 最大细格子步数，收敛后可提前停止
#endif

#ifdef unsteadyFlow
        real(kind=8), parameter :: outputSnapshotInterval = 1.0d0    ! 快照与 Nu/Re 采样间隔，单位 t_ff
        real(kind=8), parameter :: reloadFileInterval = 100.0d0      ! 完整重启文件输出间隔，单位 t_ff
        real(kind=8), parameter :: outputPltFileInterval = 100.0d0   ! Tecplot 输出间隔，单位 t_ff
        real(kind=8), parameter :: unsteadyRunDuration = 1000.0d0    ! 绝对总目标 t_ff；续算只补足剩余时间
        ! 以下三个参数只控制结束后的统计窗口，不改变推进时长或采样频率；包含匹配的旧历史数据。
        real(kind=8), parameter :: unsteadyAverageStartTf = 0.5d0*unsteadyRunDuration
        real(kind=8), parameter :: unsteadyAverageEndTf = unsteadyRunDuration
        real(kind=8), parameter :: unsteadyAverageMidTf = 0.5d0*(unsteadyAverageStartTf+unsteadyAverageEndTf)
        integer(kind=4), parameter :: outputSnapshotFile = 0    ! 0: 不输出快照；1: 输出，Nu/Re 仍独立采样
        integer(kind=4), parameter :: outputPltFile = 0         ! 0: 不输出 Tecplot；1: 输出
        integer(kind=4), parameter :: outputReloadFile = 0      ! 0: 不输出重启文件；1: 输出
        integer(kind=4), parameter :: itc_max = max(1, ceiling(unsteadyRunDuration*timeUnit))
#endif
        ! 多块输出和最终时刻向上对齐到粗细同步步；各间隔至少为一个粗步。
        ! nextSample/nextReload/nextPlt 是输出时钟，和各文件编号分开；禁用某类文件不影响其他输出。

        ! 输出文件命名与格式版本
        character(*), parameter :: settingsFile = 'SimulationSettings2DOpenaccMultiblock.txt'
        character(*), parameter :: snapshotFilePrefix = 'buoyancyCavity2DOpenaccMultiblockSnapshot'
        character(*), parameter :: pltFolderPrefix = 'buoyancyCavity2DOpenaccMultiblockTecplot'
        character(*), parameter :: reloadFilePrefix = 'reloadFile2DOpenaccMultiblock'
        character(*), parameter :: historyFile = 'NuRe_2DOpenaccMultiblock.dat'
        character(16), parameter :: restartMagic = 'MB2DRESTART0010'
        character(16), parameter :: snapshotMagic = 'MB2DSNAPSHOT0003'

        !===============================================================================================
        ! 格子方向、块数据与接口交换数据
        !===============================================================================================
        ! 每个节点、每个时间层交换 20 个数值：4 个宏观量 + 2 个力分量 + 9 个流场矩 + 5 个温度矩。
        ! p 的位置：1:4 为 rho,u,v,T；5:6 为 Fx/h,Fy/h；7:15 为归一化流场非平衡矩；16:20 为归一化温度非平衡矩。
        ! h 也表示本块时间步；流场矩包含力修正，温度场采用无热流历史修正的旧 MRT 算法。
        ! 20 是统一存储长度，守恒矩对应项仍保留为零；接口不再传递热流历史差分。
        integer(kind=4), parameter :: packetSize = 20
        ! 最多保存两块数据：1 为中心粗块，2 为外围连通细环；不是粗细比，也不是四周细区的数量。
        ! nBlocks 是实际使用的块数：多块模式为 2，refineRatio=2 的单块模式为 1；仅修改此上限不会自动增加网格块。
        integer(kind=4), parameter :: maxBlocks = 2
        integer(kind=4) :: ex(0:8) = [0, 1, 0, -1, 0, 1, -1, -1, 1]
        integer(kind=4) :: ey(0:8) = [0, 0, 1, 0, -1, 1, 1, -1, -1]
        real(kind=8) :: omega(0:8), omegaT(0:4)
        ! nBlocks 为实际块数，初始化时设为 2（中心粗块+细环）或 1（单块）。
        ! itc 为累计细时间步数，新计算从 0 开始，续算从检查点恢复；粗细比为 r 时，每个粗步累计增加 r。
        ! snapshotFileNum、pltFileNum 分别为快照和 Tecplot 文件编号，仅在实际输出相应文件时增加 1。
        integer(kind=4) :: nBlocks
        integer(kind=4) :: itc = 0
        integer(kind=4) :: snapshotFileNum = 0, pltFileNum = 0
        ! nextSample、nextReload、nextPlt 分别安排下一次 Nu/Re 采样、续算文件和 Tecplot 输出。
        ! 目标时间为对应序号乘输出间隔，再向上对齐到粗细同步步；初值 1 表示第一次计划输出。
        ! 这些序号与文件编号分开：关闭快照仍按时采样，nextSample 继续增加，snapshotFileNum 不增加。
        ! 续算时，各计划序号和文件编号均从检查点恢复，不重新从 1 或 0 开始。
        integer(kind=4) :: nextSample = 1, nextReload = 1, nextPlt = 1
        real(kind=8) :: errorU = 100.0d0, errorT = 100.0d0

        ! 块编号：1 为中心粗块，2 为连通细网格环；单块模式只使用 1。
        ! 细环采用全域二维存储，中心空区不推进
        ! 几何信息用普通数组保存，最后一个下标均为块编号。
        ! blockNi/blockNj 为本块含重叠层的节点数；blockNh 为历史时间层的最大下标。
        ! blockIlo:Ihi、blockJlo:Jhi 是x/y方向积分范围的起止下标，实际面积仍由积分权重确定。lo 表示下限，hi 表示上限
        integer(kind=4) :: blockNi(maxBlocks), blockNj(maxBlocks), blockNh(maxBlocks)
        integer(kind=4) :: blockIlo(maxBlocks), blockIhi(maxBlocks), blockJlo(maxBlocks), blockJhi(maxBlocks)
        real(kind=8) :: blockH(maxBlocks), blockX0(maxBlocks), blockY0(maxBlocks)     !中心粗块 blockH=2、blockX0=122.5
        ! blockH 同时表示本块格距和时间步；blockSn/Sq/Qk/Qn 为按该格距缩放后的松弛率, 本块对应的松弛率
        real(kind=8) :: blockSn(maxBlocks), blockSq(maxBlocks), blockQk(maxBlocks), blockQn(maxBlocks), blockGb(maxBlocks)
        ! 下面的参数四个位置依次是：xmin、xmax、ymin、ymax
        real(kind=8) :: blockOwnedBox(4, maxBlocks)    ! 积分外框：粗块直接使用；细环须扣除粗块内框面积
        real(kind=8) :: blockBaseBox(4, maxBlocks)     ! 加重叠层前的基准边界，粗细块共用指定交界面
        logical :: blockWall(4, maxBlocks)           ! 左、右、下、上是否为真实物理壁面，中心粗块：false、false、false、false，
                                                     ! 外围细环：true、true、true、true
        ! blockX0/blockY0 为首节点减去本块半格距的虚拟面，不等于积分分界。

        ! 各块按顺序存入连续数组，避免把所有块都扩充到最大宽高而浪费显存。
        ! 下列偏移从 0 开始；nodeOffset 用于宏观量，haloOffset 用于含迁移外圈的分布函数。
        integer(kind=4) :: nodeOffset(0:maxBlocks), haloOffset(0:maxBlocks), packetOffset(0:maxBlocks)
        integer(kind=4) :: xWeightOffset(0:maxBlocks), yWeightOffset(0:maxBlocks)
        real(kind=8), allocatable, target :: fStorage(:), gStorage(:), f_postStorage(:), g_postStorage(:)
        real(kind=8), allocatable, target :: rhoStorage(:), uStorage(:), vStorage(:), TStorage(:), FxStorage(:), &
            FyStorage(:)
        real(kind=8), allocatable, target :: pStorage(:), dxWeightStorage(:), dyWeightStorage(:)
#ifdef steadyFlow
        real(kind=8), allocatable, target :: upStorage(:), vpStorage(:), TpStorage(:)
#endif

        ! select_block(b) 仅将下面的数组指向第 b 块，不复制数据。
        ! 这样块内仍按 u(i,j)、f(i,j,k) 读写；切换块前须等当前主机调用返回。
        ! GPU 内核通过显式尺寸的形参接收本块数组，不读取这些可切换的主机指针。
        real(kind=8), pointer, contiguous :: f(:, :, :), g(:, :, :), f_post(:, :, :), g_post(:, :, :)
        real(kind=8), pointer, contiguous :: rho(:, :), u(:, :), v(:, :), T(:, :), Fx(:, :), Fy(:, :)
        real(kind=8), pointer, contiguous :: dxWeight(:), dyWeight(:)
        ! p 的 20 个量：rho,u,v,T,Fx/dt,Fy/dt,Kf(0:8),Kg(0:4)。
        ! K=S/dt*(m-meq+F_lattice/2)；粗块存三个时间层，细块只存当前层。
        real(kind=8), pointer, contiguous :: p(:, :, :, :)
#ifdef steadyFlow
        real(kind=8), pointer, contiguous :: up(:, :), vp(:, :), Tp(:, :)
#endif

        ! 接口连接：receiver 接收、donor 提供数据；每条连接只保存实际接收节点。
        integer(kind=4) :: nLinks = 0
        integer(kind=4) :: linkReceiver(maxBlocks*maxBlocks), linkDonor(maxBlocks*maxBlocks)
        integer(kind=4) :: linkCount(maxBlocks*maxBlocks), linkOffset(0:maxBlocks*maxBlocks)
        integer(kind=4), allocatable, target :: tiStorage(:), tjStorage(:), siStorage(:), sjStorage(:)
        logical, allocatable, target :: sameStorage(:)    ! 共址节点不做空间插值，仍重标定非平衡矩
        real(kind=8), allocatable, target :: wxStorage(:), wyStorage(:), valuesStorage(:)
        ! select_link(l) 选择一条连接；wx/wy 为四点拉格朗日权重，values 为交换暂存量。
        integer(kind=4), pointer, contiguous :: linkTi(:), linkTj(:), linkSi(:), linkSj(:)
        logical, pointer, contiguous :: linkSame(:)
        real(kind=8), pointer, contiguous :: linkWx(:, :), linkWy(:, :), linkValues(:, :)

        ! 外部标量函数的返回类型；具体函数排在主程序之后。
        integer(kind=4), external :: scheduled_step
        logical, external :: skin_node, fine_active
        real(kind=8), external :: owned_cell_area, section_owned_weight

    end module commondata

    !===============================================================================================
    ! 主程序：初始化、时间推进、输出及释放设备数据
    !===============================================================================================
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


    !===============================================================================================
    ! 初始化：构造计算块、积分分区及接口连接

    !===============================================================================================
    ! 子程序: initial
    ! 作用: 初始化网格、积分分区与接口连接
    !===============================================================================================
    subroutine initial()

        use commondata, only: nx, ny, refineRatio, fineLayerCellsLeft, fineLayerCellsRight, &
            fineLayerCellsBottom, fineLayerCellsTop, &
            overlapCells, interfaceSkin, loadInitField, reloadFileNum, Rayleigh, Prandtl, Mach, tauf, &
            viscosity, diffusivity, gBeta, timeUnit, Snu, Sq, paraA, Qk, Qnu, thermalA, &
            outputSnapshotInterval, reloadFileInterval, outputPltFileInterval, settingsFile, historyFile, &
            omega, omegaT, nBlocks, blockNi, blockNj, blockIlo, blockIhi, blockJlo, blockJhi, blockH, &
            blockX0, blockY0, blockSn, blockSq, blockQk, blockQn, blockGb, blockOwnedBox, blockBaseBox, &
            dxWeight, dyWeight, nLinks, linkReceiver, linkDonor, linkCount, linkSame, owned_cell_area
        implicit none

        integer(kind=4) :: b, k, overlap, i, j
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
        if (refineRatio == 1) then
            nBlocks = 1
            call make_block(1, 0.0d0, dble(nx), 0.0d0, dble(ny), 1, 0, 0)
        else
#if defined(VerticalWallsPeriodicalU) || defined(VerticalWallsPeriodicalT)
            error stop 'Multiblock periodic sides require a periodic block topology; use wall BCs or refineRatio=2'
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

            nBlocks = 2
            ! 两级各有一套状态：中心粗矩形与连通细环。
            call make_block(1, xLeft, xRight, yBottom, yTop, refineRatio, overlap, 2)
            call make_block(2, 0.0d0, dble(nx), 0.0d0, dble(ny), 1, 0, 0)
        endif
        ! 先确定所有块的大小，再分配连续存储；积分权重和场量随后逐块初始化。
        call allocate_block_arrays()
        totalArea = 0.0d0
        do b = 1, nBlocks
            ! 以最终统计分区统一生成权重；中心两端均为粗节点，自动得到端点半权重的梯形积分。
            call select_block(b)
            call integration_weights(blockNi(b), blockX0(b), blockH(b), blockOwnedBox(1, b), blockOwnedBox(2, b), dxWeight, &
                blockIlo(b), blockIhi(b))
            call integration_weights(blockNj(b), blockY0(b), blockH(b), blockOwnedBox(3, b), blockOwnedBox(4, b), dyWeight, &
                blockJlo(b), blockJhi(b))
            call initial_block(b)
            do j = 1, blockNj(b)
                do i = 1, blockNi(b)
                    totalArea = totalArea+owned_cell_area(b,i,j)
                enddo
            enddo
        enddo
        if (abs(totalArea-dble(nx)*ny) > 1.0d-8) error stop 'Block ownership does not tile the physical domain'
        if (loadInitField == 1) call read_restart()
        call build_links()
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
        write(k, *) 'Geometry: connected fine ring (block 2) and central coarse block (1); no fine/fine links.'
        write(k, *) 'Fine storage is rectangular; the central inactive hole is excluded from advance and integration.'
        write(k, *) 'Central spans must contain whole coarse spacings; interfaces and integration partitions are not shifted.'
        write(k, *) 'Node alignment: coarse nodes remain a subset of fine nodes; the integer phase depends on block origin.'
        write(k, *) 'Integration: clipped nodal control areas; shared coordinates do not duplicate physical area.'
        write(k, *) 'All nonconserved relaxation times satisfy h*(1/s-1/2)=constant.'
        write(k, *) 'Time: coarse prediction; prepare fine buffers BEFORE each substep; endpoint synchronization.'
        write(k, *) 'Buffers: two local layers for original split flow/thermal sequence; donors exclude these layers.'
        write(k, *) 'First coarse interval uses linear startup; subsequent intervals use three-time Lagrange interpolation.'
        write(k, *) 'Output clocks are absolute and rounded UP to a synchronized coarse step; time columns contain actual times.'
        write(k, *) 'Restart stores coarse time history; uniform-grid and earlier multiblock restart files are incompatible.'
        write(k, *) 'Thermal scheme: original legacy D2Q5; paraA is fixed across blocks.'
        do b = 1, nBlocks
            write(k, *) 'block,ni,nj,h,x0,y0,owned ilo,ihi,jlo,jhi:', b, blockNi(b), blockNj(b), &
                blockH(b), blockX0(b), blockY0(b), blockIlo(b), blockIhi(b), blockJlo(b), blockJhi(b)
            write(k, *) 'Snu,Sq,Qk,Qnu,gBeta:', blockSn(b), blockSq(b), blockQk(b), blockQn(b), blockGb(b)
            write(k, *) 'Owned physical rectangle:', blockOwnedBox(:, b)
            write(k, *) 'Base rectangle before overlap:', blockBaseBox(:, b)
            write(k, *) 'Computed first/last node x,y:', blockX0(b)+0.5d0*blockH(b), &
                blockX0(b)+(dble(blockNi(b))-0.5d0)*blockH(b), &
                blockY0(b)+0.5d0*blockH(b), blockY0(b)+(dble(blockNj(b))-0.5d0)*blockH(b)
        enddo
        do b = 1, nLinks
            call select_link(b)
            write(k, *) 'receiver, donor, direct nodes, interpolated nodes:', linkReceiver(b), linkDonor(b), &
                count(linkSame), linkCount(b)-count(linkSame)
        enddo
        close(k)
        if (loadInitField == 0) then
            open(newunit = k, file = historyFile, status = 'replace')
            write(k, '(a)') '# t_ff NuVolAvg ReVolRMS Nu_hot Nu_cold Nu_middle mass meanT Tmin Tmax rhoMin rhoMax'
            close(k)
        else
            call check_history()
        endif
    end subroutine initial
    !===============================================================================================


    !===============================================================================================
    ! 计算块：保持指定的粗细共址交界面，再向外延伸重叠层

    !===============================================================================================
    ! 子程序: make_block
    ! 作用: 确定块的基准边界、实际节点数和本块输运系数
    !===============================================================================================
    subroutine make_block(b, xlo, xhi, ylo, yhi, spacing, overlap, nh)

        use commondata, only: nx, ny, gBeta, Snu, Sq, Qk, Qnu, blockNi, blockNj, blockNh, blockH, blockX0, &
            blockY0, blockSn, blockSq, blockQk, blockQn, blockGb, blockOwnedBox, blockBaseBox, blockWall
        implicit none

        integer(kind=4), intent(in) :: b
        real(kind=8), intent(in) :: xlo, xhi, ylo, yhi
        integer(kind=4), intent(in) :: spacing, overlap, nh
        real(kind=8) :: xb, xe, yb, ye

        blockH(b) = dble(spacing)
        blockOwnedBox(:, b) = [xlo, xhi, ylo, yhi]
        blockBaseBox(:, b) = blockOwnedBox(:, b)
        if (spacing > 1) then
            ! 中心跨度已在 initial 中按整数检查；不取 ceiling，也不改变右/上交界面。
            ! 默认基准为 [127.5,895.5]，各侧再延伸 4，得到计算节点范围 [123.5,899.5]。
            if (min(xlo, ylo)-dble(overlap) <= 0.5d0 .or. &
                blockBaseBox(2, b)+dble(overlap) >= dble(nx)-0.5d0 .or. &
                blockBaseBox(4, b)+dble(overlap) >= dble(ny)-0.5d0) &
                error stop 'Aligned coarse block must remain inside the fine wall layers'
        endif
        ! 各块从自身基准边界延伸；只有细块接触物理墙，在墙处保留半个细格距。
        xb = max(0.5d0, blockBaseBox(1, b)-dble(overlap))
        xe = min(dble(nx)-0.5d0, blockBaseBox(2, b)+dble(overlap))
        yb = max(0.5d0, blockBaseBox(3, b)-dble(overlap))
        ye = min(dble(ny)-0.5d0, blockBaseBox(4, b)+dble(overlap))
        ! xb、yb 已是首节点坐标，不再加 0.5；整数粗细比保证粗节点属于细节点子集。
        blockX0(b) = xb-0.5d0*blockH(b)
        blockY0(b) = yb-0.5d0*blockH(b)
        ! 基准跨度和 overlap 均为粗格距整数倍；末节点必须恰好到达 xe、ye。
        blockNi(b) = nint((xe-xb)/blockH(b))+1
        blockNj(b) = nint((ye-yb)/blockH(b))+1
        blockNh(b) = nh
        if (abs(xb+dble(blockNi(b)-1)*blockH(b)-xe) > 1.0d-10 .or. abs(yb+dble(blockNj(b)-1)*blockH(b)-ye) > 1.0d-10) &
            error stop 'Computed block span must contain whole grid spacings'
        blockWall(:, b) = [xb == 0.5d0, xe == dble(nx)-0.5d0, yb == 0.5d0, ye == dble(ny)-0.5d0]
        ! Eq. (20): h_d*(tau_d-1/2)=h_s*(tau_s-1/2)，包括非水动力矩。
        blockSn(b) = 1.0d0/(0.5d0+(1.0d0/Snu-0.5d0)/blockH(b))
        blockSq(b) = 1.0d0/(0.5d0+(1.0d0/Sq-0.5d0)/blockH(b))
        blockQk(b) = 1.0d0/(0.5d0+(1.0d0/Qk-0.5d0)/blockH(b))
        blockQn(b) = 1.0d0/(0.5d0+(1.0d0/Qnu-0.5d0)/blockH(b))
        blockGb(b) = blockH(b)*gBeta
    end subroutine make_block
    !===============================================================================================


    !===============================================================================================
    ! 子程序: allocate_block_arrays
    ! 作用: 根据每块实际节点数一次性分配场量、迁移外圈和历史数组
    !===============================================================================================
    subroutine allocate_block_arrays()

        use commondata, only: packetSize, nBlocks, blockNi, blockNj, blockNh, nodeOffset, haloOffset, &
            packetOffset, xWeightOffset, yWeightOffset, fStorage, gStorage, f_postStorage, g_postStorage, &
            rhoStorage, uStorage, vStorage, TStorage, FxStorage, FyStorage, &
            pStorage, dxWeightStorage, dyWeightStorage
#ifdef steadyFlow
        use commondata, only: upStorage, vpStorage, TpStorage
#endif
        implicit none

        integer(kind=4) :: b, n, nhalo

        nodeOffset(0) = 0
        haloOffset(0) = 0
        packetOffset(0) = 0
        xWeightOffset(0) = 0
        yWeightOffset(0) = 0
        do b = 1, nBlocks
            n = blockNi(b)*blockNj(b)
            nodeOffset(b) = nodeOffset(b-1)+n
            haloOffset(b) = haloOffset(b-1)+(blockNi(b)+2)*(blockNj(b)+2)
            packetOffset(b) = packetOffset(b-1)+n*packetSize*(blockNh(b)+1)
            xWeightOffset(b) = xWeightOffset(b-1)+blockNi(b)
            yWeightOffset(b) = yWeightOffset(b-1)+blockNj(b)
        enddo
        n = nodeOffset(nBlocks)
        nhalo = haloOffset(nBlocks)
        allocate(fStorage(9*n), gStorage(5*n), f_postStorage(9*nhalo), g_postStorage(5*nhalo))
        allocate(rhoStorage(n), uStorage(n), vStorage(n), TStorage(n))
        allocate(FxStorage(n), FyStorage(n))
        allocate(pStorage(packetOffset(nBlocks)))
        allocate(dxWeightStorage(xWeightOffset(nBlocks)), dyWeightStorage(yWeightOffset(nBlocks)))
#ifdef steadyFlow
        allocate(upStorage(n), vpStorage(n), TpStorage(n))
#endif
    end subroutine allocate_block_arrays
    !===============================================================================================

    !===============================================================================================
    ! 子程序: select_block
    ! 作用: 让块内数组指向第 b 块的连续存储，保留二维和分布函数下标
    !===============================================================================================
    subroutine select_block(b)

        use commondata, only: packetSize, blockNi, blockNj, blockNh, nodeOffset, haloOffset, packetOffset, &
            xWeightOffset, yWeightOffset, fStorage, gStorage, f_postStorage, g_postStorage, rhoStorage, &
            uStorage, vStorage, TStorage, FxStorage, FyStorage, pStorage, &
            dxWeightStorage, dyWeightStorage, f, g, f_post, g_post, rho, u, v, T, Fx, Fy, &
            dxWeight, dyWeight, p
#ifdef steadyFlow
        use commondata, only: upStorage, vpStorage, TpStorage, up, vp, Tp
#endif
        implicit none

        integer(kind=4), intent(in) :: b
        integer(kind=4) :: ni, nj, first, last

        ni = blockNi(b)
        nj = blockNj(b)
        first = nodeOffset(b-1)+1
        last = nodeOffset(b)
        ! 指针重设下标只改变访问方式；不会分配数组，也不会复制流场。
        rho(1:ni, 1:nj) => rhoStorage(first:last)
        u(1:ni, 1:nj) => uStorage(first:last)
        v(1:ni, 1:nj) => vStorage(first:last)
        T(1:ni, 1:nj) => TStorage(first:last)
        Fx(1:ni, 1:nj) => FxStorage(first:last)
        Fy(1:ni, 1:nj) => FyStorage(first:last)
#ifdef steadyFlow
        up(1:ni, 1:nj) => upStorage(first:last)
        vp(1:ni, 1:nj) => vpStorage(first:last)
        Tp(1:ni, 1:nj) => TpStorage(first:last)
#endif
        f(1:ni, 1:nj, 0:8) => fStorage(9*nodeOffset(b-1)+1:9*last)
        g(1:ni, 1:nj, 0:4) => gStorage(5*nodeOffset(b-1)+1:5*last)
        f_post(0:ni+1, 0:nj+1, 0:8) => f_postStorage(9*haloOffset(b-1)+1:9*haloOffset(b))
        g_post(0:ni+1, 0:nj+1, 0:4) => g_postStorage(5*haloOffset(b-1)+1:5*haloOffset(b))
        p(1:ni, 1:nj, 1:packetSize, 0:blockNh(b)) => pStorage(packetOffset(b-1)+1:packetOffset(b))
        dxWeight => dxWeightStorage(xWeightOffset(b-1)+1:xWeightOffset(b))
        dyWeight => dyWeightStorage(yWeightOffset(b-1)+1:yWeightOffset(b))
    end subroutine select_block
    !===============================================================================================

    !===============================================================================================
    ! 子程序: allocate_link_arrays
    ! 作用: 按每条连接的实际节点数分配接口索引、插值权重和交换缓存
    !===============================================================================================
    subroutine allocate_link_arrays()

        use commondata, only: packetSize, nLinks, linkCount, linkOffset, tiStorage, tjStorage, siStorage, &
            sjStorage, sameStorage, wxStorage, wyStorage, valuesStorage
        implicit none

        integer(kind=4) :: l, n

        linkOffset(0) = 0
        do l = 1, nLinks
            linkOffset(l) = linkOffset(l-1)+linkCount(l)
        enddo
        n = linkOffset(nLinks)
        allocate(tiStorage(n), tjStorage(n), siStorage(n), sjStorage(n), sameStorage(n))
        allocate(wxStorage(4*n), wyStorage(4*n), valuesStorage(packetSize*n))
        valuesStorage = 0.0d0
    end subroutine allocate_link_arrays
    !===============================================================================================

    !===============================================================================================
    ! 子程序: select_link
    ! 作用: 选择第 l 条接口连接，不复制索引、权重或交换数据
    !===============================================================================================
    subroutine select_link(l)

        use commondata, only: packetSize, linkCount, linkOffset, tiStorage, tjStorage, siStorage, sjStorage, &
            sameStorage, wxStorage, wyStorage, valuesStorage, linkTi, linkTj, linkSi, linkSj, linkSame, &
            linkWx, linkWy, linkValues
        implicit none

        integer(kind=4), intent(in) :: l
        integer(kind=4) :: first, last, n

        first = linkOffset(l-1)+1
        last = linkOffset(l)
        n = linkCount(l)
        linkTi => tiStorage(first:last)
        linkTj => tjStorage(first:last)
        linkSi => siStorage(first:last)
        linkSj => sjStorage(first:last)
        linkSame => sameStorage(first:last)
        linkWx(1:4, 1:n) => wxStorage(4*(first-1)+1:4*last)
        linkWy(1:4, 1:n) => wyStorage(4*(first-1)+1:4*last)
        linkValues(1:packetSize, 1:n) => valuesStorage(packetSize*(first-1)+1:packetSize*last)
    end subroutine select_link
    !===============================================================================================

    !===============================================================================================
    ! 积分权重：保证分区面积与线性坐标积分准确

    !===============================================================================================
    ! 子程序: integration_weights
    ! 作用: 构造一维积分权重并检查总长度和一阶矩
    !===============================================================================================
    subroutine integration_weights(n, origin, h, lo, hi, w, first, last)

        implicit none

        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: origin, h, lo, hi
        real(kind=8), intent(out) :: w(n)
        integer(kind=4), intent(out) :: first, last
        integer(kind=4) :: i
        real(kind=8) :: firstMoment, x, shift

        first = n+1
        last = 0
        firstMoment = 0.0d0
        do i = 1, n
            w(i) = max(0.0d0, min(hi, origin+dble(i)*h)-max(lo, origin+dble(i-1)*h))
            if (w(i) <= 0.0d0) cycle
            first = min(first, i)
            last = i
            x = origin+(dble(i)-0.5d0)*h
            firstMoment = firstMoment+w(i)*x
        enddo
        if (abs(sum(w)-(hi-lo)) > 1.0d-10) error stop 'Integration weights do not cover owned interval'
        ! 中心粗块在共址节点之间积分，无需一阶矩修正；靠墙细块的裁剪仍可能不对称。
        ! 在末两个积分节点之间转移权重，补齐一阶矩，保持总长度。
        ! 两节点相距 h，转移 shift 后一阶矩增加 h*shift；无需修正的分区保持 shift=0。
        if (last > first) then
            shift = (0.5d0*(hi**2-lo**2)-firstMoment)/h
            w(last-1) = w(last-1)-shift
            w(last) = w(last)+shift
        endif
        if (any(w < 0.0d0)) error stop 'Negative corrected integration weight'
        firstMoment = 0.0d0
        do i = first, last
            firstMoment = firstMoment+w(i)*(origin+(dble(i)-0.5d0)*h)
        enddo
        if (abs(firstMoment-0.5d0*(hi**2-lo**2)) > 1.0d-9*max(1.0d0, abs(firstMoment))) &
            error stop 'Integration weights do not integrate a linear coordinate exactly'
    end subroutine integration_weights
    !===============================================================================================


    !===============================================================================================
    ! 初始化各块的流场、温度场和分布函数

    !===============================================================================================
    ! 子程序: initial_block
    ! 作用: 初始化各块的宏观量与分布函数
    !===============================================================================================
    subroutine initial_block(b)

        use commondata, only: nx, ny, Thot, Tcold, pi, lengthUnit, omega, omegaT, blockNi, blockNj, blockH, &
            blockX0, blockY0, f, g, f_post, g_post, rho, u, v, T, Fx, Fy, p
#ifdef steadyFlow
        use commondata, only: up, vp, Tp
#endif
        implicit none

        integer(kind=4), intent(in) :: b
        integer(kind=4) :: i, j, a
        real(kind=8) :: x, y

        call select_block(b)
        u = 0.0d0
        v = 0.0d0
        rho = 1.0d0
        Fx = 0.0d0
        Fy = 0.0d0
        f_post = 0.0d0
        g_post = 0.0d0
        p = 0.0d0
        do j = 1, blockNj(b)
            y = (blockY0(b)+(dble(j)-0.5d0)*blockH(b))/lengthUnit
            do i = 1, blockNi(b)
                x = (blockX0(b)+(dble(i)-0.5d0)*blockH(b))/lengthUnit
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
    end subroutine initial_block
    !===============================================================================================


    !===============================================================================================
    ! OpenACC 数据区：传入数组并准备首次碰撞的接口状态

    !===============================================================================================
    ! 子程序: enter_data_2d_openacc
    ! 作用: 建立设备数据区并准备初始接口状态
    !===============================================================================================
    subroutine enter_data_2d_openacc()

        use commondata, only: loadInitField, ex, ey, omega, omegaT, nBlocks, blockNh
        implicit none

        integer(kind=4) :: b, l
        !$acc enter data copyin(ex, ey, omega, omegaT)
        call block_device_data(.true.)
        call link_device_data(.true.)
        if (loadInitField == 0) then
            do b = 1, nBlocks
                call pack_block(b, min(1, blockNh(b)))
            enddo
            if (nBlocks > 1) then
                ! 首次碰撞前补齐两级接口。两次交换均读取同一份初始快照，不依赖块处理顺序。
                call exchange_interfaces(.true., [0.0d0, 1.0d0, 0.0d0])
                call exchange_interfaces(.false., [0.0d0, 1.0d0, 0.0d0])
                do b = 1, nBlocks
                    call pack_block(b, min(1, blockNh(b)))
                enddo
            endif
        endif
    end subroutine enter_data_2d_openacc
    !===============================================================================================


    !===============================================================================================
    ! 子程序: block_device_data
    ! 作用: 整体建立或释放所有计算块的设备存储
    !===============================================================================================
    subroutine block_device_data(entering)

        use commondata, only: fStorage, gStorage, f_postStorage, g_postStorage, rhoStorage, uStorage, &
            vStorage, TStorage, FxStorage, FyStorage, pStorage
        implicit none

        logical, intent(in) :: entering

        ! 整体映射连续存储；各块形参只是其中的连续片段，不在设备端切换指针。
        if (entering) then
            !$acc enter data copyin(fStorage, gStorage, f_postStorage, g_postStorage, rhoStorage, uStorage, &
            !$acc& vStorage, TStorage, FxStorage, FyStorage, pStorage)
        else
            !$acc exit data delete(fStorage, gStorage, f_postStorage, g_postStorage, rhoStorage, uStorage, &
            !$acc& vStorage, TStorage, FxStorage, FyStorage, pStorage)
        endif
    end subroutine block_device_data
    !===============================================================================================


    !===============================================================================================
    ! 子程序: link_device_data
    ! 作用: 管理接口连接的设备数组
    !===============================================================================================
    subroutine link_device_data(entering)

        use commondata, only: nLinks, tiStorage, tjStorage, siStorage, sjStorage, sameStorage, wxStorage, &
            wyStorage, valuesStorage
        implicit none

        logical, intent(in) :: entering

        if (nLinks == 0) return
        if (entering) then
            !$acc enter data copyin(tiStorage, tjStorage, siStorage, sjStorage, wxStorage, wyStorage, &
            !$acc& sameStorage, valuesStorage)
        else
            !$acc exit data delete(tiStorage, tjStorage, siStorage, sjStorage, wxStorage, wyStorage, sameStorage, valuesStorage)
        endif
    end subroutine link_device_data
    !===============================================================================================


    !===============================================================================================
    ! 子程序: update_host_block
    ! 作用: 将指定块的设备状态更新到主机
    !===============================================================================================
    subroutine update_host_block(b, full)

        use commondata, only: f, g, rho, u, v, T, Fx, Fy, p
        implicit none

        integer(kind=4), intent(in) :: b
        logical, intent(in) :: full

        call select_block(b)
        !$acc update self(u, v, T, rho) async(1)
        if (full) then
            !$acc update self(f, g, Fx, Fy, p) async(1)
        endif
    end subroutine update_host_block
    !===============================================================================================


    !===============================================================================================
    ! 子程序: update_host_all
    ! 作用: 将全部计算块的设备状态更新到主机
    !===============================================================================================
    subroutine update_host_all(full)

        use commondata, only: nBlocks
        implicit none

        logical, intent(in) :: full
        integer(kind=4) :: b

        do b = 1, nBlocks
            call update_host_block(b, full)
        enddo
        !$acc wait(1)
    end subroutine update_host_all
    !===============================================================================================


    !===============================================================================================
    ! 子程序: exit_data_2d_openacc
    ! 作用: 释放接口与计算块的设备数组
    !===============================================================================================
    subroutine exit_data_2d_openacc()

        use commondata, only: ex, ey, omega, omegaT
        implicit none

        integer(kind=4) :: b, l
        !$acc wait(1)
        call link_device_data(.false.)
        call block_device_data(.false.)
        !$acc exit data delete(ex, ey, omega, omegaT)
    end subroutine exit_data_2d_openacc
    !===============================================================================================


    !===============================================================================================
    ! 块内推进：保持原流场与温度场的执行次序

    !===============================================================================================
    ! 子程序: advance_block
    ! 作用: 按原算法次序推进单个计算块
    !===============================================================================================
    subroutine advance_block(b)

        use commondata, only: blockNi, blockNj, blockH, blockSn, blockSq, blockQk, blockQn, blockGb, &
            blockWall, f, g, f_post, g_post, rho, u, v, T, Fx, Fy
#ifdef SideHeatedHa
        use commondata, only: B2sigemarho
#endif
        implicit none

        integer(kind=4), intent(in) :: b
        call select_block(b)
        ! 原文件的执行次序保持不变。粗块的 dt 已吸收到松弛率及力增量中。
        ! 入口：人工边界已重建为本时刻的碰撞前状态；缓冲节点与内部节点一样参与碰撞。
        ! 出口：最外两层只作为待重建缓冲，不允许充当 donor；原算法的内部结果保留。
        call collision(blockNi(b), blockNj(b), f, f_post, rho, u, v, Fx, Fy, T, blockSn(b), blockSq(b), blockGb(b) &
#ifdef SideHeatedHa
            , blockH(b)*B2sigemarho &
#endif
            )
        call streaming(blockNi(b), blockNj(b), f, f_post)
        call bounceback(blockNi(b), blockNj(b), f, f_post, blockWall(1, b), blockWall(2, b), blockWall(3, b), blockWall(4, b))
        call macro(blockNi(b), blockNj(b), f, rho, u, v, Fx, Fy)
        call collisionT(blockNi(b), blockNj(b), g, g_post, u, v, T, blockQk(b), blockQn(b))
        call streamingT(blockNi(b), blockNj(b), g, g_post)
        call bouncebackT(blockNi(b), blockNj(b), g, g_post, blockWall(1, b), blockWall(2, b), blockWall(3, b), blockWall(4, b))
        call macroT(blockNi(b), blockNj(b), g, T)
    end subroutine advance_block
    !===============================================================================================


    !===============================================================================================
    ! 多块推进：粗块预测、细块子步及两级同步

    !===============================================================================================
    ! 子程序: advance_multiblock
    ! 作用: 执行粗块预测、细步推进和同步交换
    !===============================================================================================
    subroutine advance_multiblock()

        use commondata, only: refineRatio, nBlocks, itc, blockNi, blockNj, p
        implicit none

        integer(kind=4) :: b, k
        real(kind=8) :: theta, wt(0:2)

        if (nBlocks == 1) then
            call advance_block(1)
            itc = itc+1
            return
        endif
        ! 1. 入口是粗细同步的 t 时刻。粗块缓冲已在初始交换或上次同步时补齐。
        ! 2. 粗块先预测 t+dt_c；p(:,:,:,0:2) 分别代表 t-dt_c、t、t+dt_c。
        call advance_block(1)
        call pack_block(1, 2)
        do k = 1, refineRatio
            ! 3. 第 k 个细子步在 t+(k-1)*dt_f 开始：先补齐该时刻的碰撞前缓冲，再推进。
            ! 唯一细环的内边缘由粗级历史插值；细环内部不交换数据。
            theta = dble(k-1)/dble(refineRatio)
            call coarse_time_weights(theta, wt)
            call exchange_interfaces(.false., wt)
            do b = 2, nBlocks
                call advance_block(b)
                call pack_block(b, 0)
            enddo
        enddo
        ! 4. 两级均到 t+dt_c：细->粗修复粗缓冲，随后补齐细缓冲，供下次碰撞和同步输出使用。
        call exchange_interfaces(.true., [0.0d0, 0.0d0, 1.0d0])
        call pack_block(1, 2)
        call exchange_interfaces(.false., [0.0d0, 0.0d0, 1.0d0])
        ! 5. 只在同步点滚动历史和整数时钟，重启文件必须保存完整历史。
        call select_block(1)
        call rotate_coarse_history(blockNi(1), blockNj(1), p)
        itc = itc+refineRatio
    end subroutine advance_multiblock
    !===============================================================================================


    !===============================================================================================
    ! 子程序: coarse_time_weights
    ! 作用: 计算粗时间层到细子步起点的插值权重
    !===============================================================================================
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
    !===============================================================================================


    !===============================================================================================
    ! 子程序: rotate_coarse_history
    ! 作用: 滚动保存粗块的三个时间层
    !===============================================================================================
    subroutine rotate_coarse_history(ni, nj, p)

        use commondata, only: packetSize
        implicit none

        integer(kind=4), intent(in) :: ni, nj
        real(kind=8), intent(inout) :: p(ni, nj, packetSize, 0:2)
        integer(kind=4) :: i, j, k
        !$acc parallel loop collapse(3) present(p) async(1)
        do k = 1, packetSize
            do j = 1, nj
                do i = 1, ni
                    p(i, j, k, 0) = p(i, j, k, 1)
                    p(i, j, k, 1) = p(i, j, k, 2)
                enddo
            enddo
        enddo
    end subroutine rotate_coarse_history
    !===============================================================================================


    !===============================================================================================
    ! 函数: skin_node
    ! 作用: 判断节点是否属于人工边界缓冲层
    !===============================================================================================
    logical function skin_node(b, i, j)

        use commondata, only: interfaceSkin, blockNi, blockNj, blockWall, nBlocks, blockH, fine_active, &
            fineLayerCellsLeft, fineLayerCellsRight, fineLayerCellsBottom, fineLayerCellsTop, nx, ny, overlapCells, refineRatio
        implicit none

        integer(kind=4), intent(in) :: b
        integer(kind=4), intent(in) :: i, j

        real(8) :: x,y,xl,xr,yb,yt
        if (nBlocks>1 .and. b==2) then
            x=dble(i)-0.5d0; y=dble(j)-0.5d0
            xl=fineLayerCellsLeft-0.5d0+overlapCells*refineRatio
            xr=nx-fineLayerCellsRight+0.5d0-overlapCells*refineRatio
            yb=fineLayerCellsBottom-0.5d0+overlapCells*refineRatio
            yt=ny-fineLayerCellsTop+0.5d0-overlapCells*refineRatio
            skin_node=fine_active(i,j) .and. x>xl-interfaceSkin .and. x<xr+interfaceSkin .and. &
                y>yb-interfaceSkin .and. y<yt+interfaceSkin
            return
        endif
        skin_node = (.not.blockWall(1, b) .and. i <= interfaceSkin) .or. &
            (.not.blockWall(2, b) .and. i > blockNi(b)-interfaceSkin) .or. &
            (.not.blockWall(3, b) .and. j <= interfaceSkin) .or. &
            (.not.blockWall(4, b) .and. j > blockNj(b)-interfaceSkin)
    end function skin_node
    !===============================================================================================


    !===============================================================================================
    ! 接口来源：排除人工边界缓冲层，检查实际插值模板

    !===============================================================================================
    ! 子程序: donor_stencil
    ! 作用: 选择合法来源块及空间插值模板
    !===============================================================================================
    subroutine donor_stencil(receiver, x, y, donor, si, sj, wx, wy, coincident)

        use commondata, only: interfaceSkin, nBlocks, blockNi, blockNj, blockH, blockX0, blockY0, blockWall, fine_active, skin_node
        implicit none

        integer(kind=4), intent(in) :: receiver
        real(kind=8), intent(in) :: x, y
        integer(kind=4), intent(out) :: donor, si, sj
        real(kind=8), intent(out) :: wx(4), wy(4)
        logical, intent(out) :: coincident
        integer(kind=4) :: d, il, ih, jl, jh, is, js
        real(kind=8) :: qx, qy, score, best

        best = -huge(1.0d0)
        donor = 0
        do d = 1, nBlocks
            if (d == receiver) cycle
            ! 来源模板严格避开刚推进后尚未修复的两层人工边界。
            il = 1
            ih = blockNi(d)
            jl = 1
            jh = blockNj(d)
            if (.not.blockWall(1, d)) il = 1+interfaceSkin
            if (.not.blockWall(2, d)) ih = ih-interfaceSkin
            if (.not.blockWall(3, d)) jl = 1+interfaceSkin
            if (.not.blockWall(4, d)) jh = jh-interfaceSkin
            qx = (x-blockX0(d))/blockH(d)+0.5d0
            qy = (y-blockY0(d))/blockH(d)+0.5d0
            if (qx < dble(il) .or. qx > dble(ih) .or. qy < dble(jl) .or. qy > dble(jh)) cycle
            if (d==2) then
                if (abs(qx-nint(qx))>1d-12 .or. abs(qy-nint(qy))>1d-12) cycle
                if (.not.fine_active(nint(qx),nint(qy))) cycle
                if (skin_node(d,nint(qx),nint(qy))) cycle
            endif
            if (ih-il < 3 .or. jh-jl < 3) cycle
            is = max(il, min(floor(qx)-1, ih-3))
            js = max(jl, min(floor(qy)-1, jh-3))
            score = min(qx-il, ih-qx, qy-jl, jh-qy)*blockH(d)
            ! 只存在粗->细或细->粗连接；细环内部直接迁移。
            if (score <= best) cycle
            best = score
            donor = d
            si = is
            sj = js
            coincident = abs(qx-dble(nint(qx))) < 1.0d-12 .and. abs(qy-dble(nint(qy))) < 1.0d-12
            if (coincident) then
                si = nint(qx)
                sj = nint(qy)
                wx = [1.0d0, 0.0d0, 0.0d0, 0.0d0]
                wy = wx
            else
                call lagrange_weights(qx-dble(is), wx)
                call lagrange_weights(qy-dble(js), wy)
            endif
        enddo
        if (donor == 0) then
            write(*, *) 'No interior four-point donor stencil:', receiver, x, y
            error stop 'Invalid overlap geometry'
        endif
    end subroutine donor_stencil
    !===============================================================================================


    !===============================================================================================
    ! 子程序: lagrange_weights
    ! 作用: 计算四点拉格朗日插值权重
    !===============================================================================================
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
    !===============================================================================================


    !===============================================================================================
    ! 子程序: section_weights
    ! 作用: 计算指定截面的插值及导数权重
    !===============================================================================================
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
    !===============================================================================================


    !===============================================================================================
    ! 建立各块之间的接收节点、来源节点和空间权重

    !===============================================================================================
    ! 子程序: build_links
    ! 作用: 建立接收节点与来源节点之间的连接
    !===============================================================================================
    subroutine build_links()

        use commondata, only: maxBlocks, nBlocks, blockNi, blockNj, blockH, blockX0, blockY0, nLinks, &
            linkReceiver, linkDonor, linkCount, linkTi, linkTj, linkSi, linkSj, linkSame, linkWx, linkWy, &
            skin_node
        implicit none

        integer(kind=4) :: b, d, i, j, l, n, si, sj, counts(maxBlocks, maxBlocks), idx(maxBlocks, maxBlocks)
        real(kind=8) :: x, y, wx(4), wy(4)
        logical :: coincident

        counts = 0
        idx = 0
        nLinks = 0
        do b = 1, nBlocks
            do j = 1, blockNj(b)
                do i = 1, blockNi(b)
                    if (.not.skin_node(b, i, j)) cycle
                    x = blockX0(b)+(dble(i)-0.5d0)*blockH(b)
                    y = blockY0(b)+(dble(j)-0.5d0)*blockH(b)
                    call donor_stencil(b, x, y, d, si, sj, wx, wy, coincident)
                    counts(b, d) = counts(b, d)+1
                enddo
            enddo
        enddo
        do b = 1, nBlocks
            do d = 1, nBlocks
                if (counts(b, d) == 0) cycle
                nLinks = nLinks+1
                l = nLinks
                idx(b, d) = l
                n = counts(b, d)
                linkReceiver(l) = b
                linkDonor(l) = d
                linkCount(l) = n
            enddo
        enddo
        ! 第一遍只计数，按实际连接长度分配；第二遍填入接收点和四点来源模板。
        call allocate_link_arrays()
        counts = 0
        do b = 1, nBlocks
            do j = 1, blockNj(b)
                do i = 1, blockNi(b)
                    if (.not.skin_node(b, i, j)) cycle
                    x = blockX0(b)+(dble(i)-0.5d0)*blockH(b)
                    y = blockY0(b)+(dble(j)-0.5d0)*blockH(b)
                    call donor_stencil(b, x, y, d, si, sj, wx, wy, coincident)
                    counts(b, d) = counts(b, d)+1
                    n = counts(b, d)
                    l = idx(b, d)
                    call select_link(l)
                    linkTi(n) = i
                    linkTj(n) = j
                    linkSi(n) = si
                    linkSj(n) = sj
                    linkWx(:, n) = wx
                    linkWy(:, n) = wy
                    linkSame(n) = coincident
                    if (blockH(b) > blockH(d) .and. .not.coincident) &
                        error stop 'Coarse interface node is not aligned with a fine node'
                enddo
            enddo
        enddo
    end subroutine build_links
    !===============================================================================================


    !===============================================================================================
    ! 子程序: equilibrium_moments
    ! 作用: 计算流场与温度场的平衡矩
    !===============================================================================================
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
    !===============================================================================================


    !===============================================================================================
    ! 子程序: force_moments
    ! 作用: 计算含力修正所需的矩空间源项
    !===============================================================================================
    subroutine force_moments(ux, uy, fx, fy, fm)

        implicit none

        !$acc routine seq
        real(kind=8), intent(in) :: ux, uy, fx, fy
        real(kind=8), intent(out) :: fm(0:8)

        fm = [0.0d0, 6.0d0*(ux*fx+uy*fy), -6.0d0*(ux*fx+uy*fy), fx, -fx, fy, -fy, &
            2.0d0*(ux*fx-uy*fy), ux*fy+uy*fx]
    end subroutine force_moments
    !===============================================================================================


    !===============================================================================================
    ! 子程序: pack_block
    ! 作用: 将指定块的状态编码到交换时间层
    !===============================================================================================
    subroutine pack_block(b, slot)

        use commondata, only: blockNi, blockNj, blockNh, blockH, blockSn, blockSq, blockQk, blockQn, f, g, &
            rho, u, v, T, Fx, Fy, p
        implicit none

        integer(kind=4), intent(in) :: b
        integer(kind=4), intent(in) :: slot

        call select_block(b)
        call encode_packets(blockNi(b), blockNj(b), blockNh(b), slot, blockH(b), blockSn(b), blockSq(b), &
            blockQk(b), blockQn(b), &
            f, g, rho, u, v, T, Fx, Fy, p)
    end subroutine pack_block
    !===============================================================================================


    !===============================================================================================
    ! 将宏观量及归一化非平衡矩编码为接口交换量

    !===============================================================================================
    ! 子程序: encode_packets
    ! 作用: 编码宏观量与归一化非平衡矩
    !===============================================================================================
    subroutine encode_packets(ni, nj, nh, slot, h, sn, sq, qk, qn, f, g, rho, u, v, T, Fx, Fy, p)

        use commondata, only: packetSize
        implicit none

        integer(kind=4), intent(in) :: ni, nj, nh, slot
        real(kind=8), intent(in) :: h, sn, sq, qk, qn, f(ni, nj, 0:8), g(ni, nj, 0:4)
        real(kind=8), intent(in) :: rho(ni, nj), u(ni, nj), v(ni, nj), T(ni, nj), Fx(ni, nj), Fy(ni, nj)
        real(kind=8), intent(inout) :: p(ni, nj, packetSize, 0:nh)
        integer(kind=4) :: i, j, a
        real(kind=8) :: m(0:8), meq(0:8), fm(0:8), n(0:4), neq(0:4), s(0:8), q(0:4), fv(0:8), gv(0:4)
        !$acc parallel loop collapse(2) present(f, g, rho, u, v, T, Fx, Fy, p) async(1) &
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
                p(i, j, 1, slot) = rho(i, j)
                p(i, j, 2, slot) = u(i, j)
                p(i, j, 3, slot) = v(i, j)
                p(i, j, 4, slot) = T(i, j)
                p(i, j, 5, slot) = Fx(i, j)/h
                p(i, j, 6, slot) = Fy(i, j)/h
                do a = 0, 8
                    p(i, j, 7+a, slot) = s(a)/h*(m(a)-meq(a)+0.5d0*fm(a))
                enddo
                do a = 0, 4
                    p(i, j, 16+a, slot) = q(a)/h*(n(a)-neq(a))
                enddo
            enddo
        enddo
    end subroutine encode_packets
    !===============================================================================================


    !===============================================================================================
    ! 接口交换：从来源块取值，再重建接收块的碰撞前状态

    !===============================================================================================
    ! 子程序: exchange_interfaces
    ! 作用: 完成指定接收层级的接口交换
    !===============================================================================================
    subroutine exchange_interfaces(coarse_receiver, wt)

        use commondata, only: blockNi, blockNj, blockNh, blockH, blockSn, blockSq, blockQk, blockQn, f, g, &
            rho, u, v, T, Fx, Fy, p, nLinks, linkReceiver, linkDonor, linkCount, linkTi, &
            linkTj, linkSi, linkSj, linkSame, linkWx, linkWy, linkValues
        implicit none

        logical, intent(in) :: coarse_receiver
        real(kind=8), intent(in) :: wt(0:2)
        integer(kind=4) :: l, b, d

        do l = 1, nLinks
            b = linkReceiver(l)
            d = linkDonor(l)
            if ((b == 1) .neqv. coarse_receiver) cycle
            call select_link(l)
            call select_block(d)
            call interpolate_packets(blockNi(d), blockNj(d), blockNh(d), p, linkCount(l), &
                linkSi, linkSj, linkWx, linkWy, linkSame, wt, linkValues)
        enddo
        do l = 1, nLinks
            b = linkReceiver(l)
            if ((b == 1) .neqv. coarse_receiver) cycle
            call select_link(l)
            call select_block(b)
            call apply_packets(blockNi(b), blockNj(b), blockH(b), blockSn(b), blockSq(b), blockQk(b), blockQn(b), &
                f, g, rho, u, v, T, &
                Fx, Fy, &
                linkCount(l), linkTi, linkTj, linkValues)
        enddo
    end subroutine exchange_interfaces
    !===============================================================================================


    !===============================================================================================
    ! 子程序: interpolate_packets
    ! 作用: 执行共址取值或空间、时间插值
    !===============================================================================================
    subroutine interpolate_packets(ni, nj, nh, p, count, si, sj, wx, wy, coincident, wt, val)

        use commondata, only: packetSize
        implicit none

        integer(kind=4), intent(in) :: ni, nj, nh, count, si(count), sj(count)
        real(kind=8), intent(in) :: p(ni, nj, packetSize, 0:nh), wx(4, count), wy(4, count), wt(0:2)
        logical, intent(in) :: coincident(count)
        real(kind=8), intent(out) :: val(packetSize, count)
        integer(kind=4) :: c, a, ix, iy, k
        real(kind=8) :: value, wk, wt0, wt1, wt2
        ! NVHPC 24.3/P100 上数组形参 firstprivate(wt) 会在设备读取 wt(k) 时非法访问。
        ! 在主机端取出三个时间权重，按标量值捕获；保持原三时间层插值及 async(1) 次序。
        wt0 = wt(0)
        wt1 = wt(1)
        wt2 = wt(2)
        !$acc parallel loop collapse(2) present(p, si, sj, wx, wy, coincident, &
        !$acc& val) firstprivate(wt0, wt1, wt2) async(1) private(ix, iy, k, value, wk)
        do c = 1, count
            do a = 1, packetSize
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
                        value = value+wk*p(si(c), sj(c), a, k)
                    else
                        do iy = 1, 4
                            do ix = 1, 4
                                value = value+wk*wx(ix, c)*wy(iy, c)*p(si(c)+ix-1, sj(c)+iy-1, a, k)
                            enddo
                        enddo
                    endif
                enddo
                val(a, c) = value
            enddo
        enddo
    end subroutine interpolate_packets
    !===============================================================================================


    !===============================================================================================
    ! 子程序: apply_packets
    ! 作用: 按接收块参数重建碰撞前状态
    !===============================================================================================
    subroutine apply_packets(ni, nj, h, sn, sq, qk, qn, f, g, rho, u, v, T, Fx, Fy, count, ti, &
        tj, val)

        use commondata, only: packetSize
        implicit none

        integer(kind=4), intent(in) :: ni, nj, count, ti(count), tj(count)
        real(kind=8), intent(in) :: h, sn, sq, qk, qn, val(packetSize, count)
        real(kind=8), intent(inout) :: f(ni, nj, 0:8), g(ni, nj, 0:4), rho(ni, nj), u(ni, nj), v(ni, nj), T(ni, nj)
        real(kind=8), intent(inout) :: Fx(ni, nj), Fy(ni, nj)
        integer(kind=4) :: c, i, j, a
        real(kind=8) :: m(0:8), meq(0:8), fm(0:8), n(0:4), neq(0:4), s(0:8), q(0:4), fv(0:8), gv(0:4)
        !$acc parallel loop present(f, g, rho, u, v, T, Fx, Fy, ti, tj, val) async(1) &
        !$acc& private(i, j, a, m, meq, fm, n, neq, s, q, fv, gv)
        do c = 1, count
            i = ti(c)
            j = tj(c)
            rho(i, j) = val(1, c)
            u(i, j) = val(2, c)
            v(i, j) = val(3, c)
            T(i, j) = val(4, c)
            Fx(i, j) = h*val(5, c)
            Fy(i, j) = h*val(6, c)
            call equilibrium_moments(rho(i, j), u(i, j), v(i, j), T(i, j), meq, neq)
            call force_moments(u(i, j), v(i, j), Fx(i, j), Fy(i, j), fm)
            s = [0.0d0, sn, sn, 0.0d0, sq, 0.0d0, sq, sn, sn]
            q = [0.0d0, qk, qk, qn, qn]
            m = meq
            n = neq
            do a = 0, 8
                if (s(a) > 0.0d0) m(a) = meq(a)+h/s(a)*val(7+a, c)-0.5d0*fm(a)
            enddo
            ! Eq. (21)：守恒矩直接用宏观量和半步力重建，不除以零松弛率。
            m(0) = rho(i, j)
            m(3) = rho(i, j)*u(i, j)-0.5d0*Fx(i, j)
            m(5) = rho(i, j)*v(i, j)-0.5d0*Fy(i, j)
            do a = 1, 4
                n(a) = neq(a)+h/q(a)*val(16+a, c)
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
    end subroutine apply_packets
    !===============================================================================================


    !===============================================================================================
    ! 以下八个块内子程序取自原文件。修改限于：显式块数组/局部尺寸/局部松弛率参数、物理壁面标记。
    !===============================================================================================

    !===============================================================================================
    ! 子程序: collision
    ! 作用: 流场多松弛碰撞及含力修正
    !===============================================================================================
    subroutine collision(nx, ny, f, f_post, rho, u, v, Fx, Fy, T, Snu, Sq, gBeta &
#ifdef SideHeatedHa
        , B2sigemarho &
#endif
        )

        use commondata, only: Tref
#ifdef SideHeatedHa
        use commondata, only: phi
#endif
        use commondata, only: globalNx=>nx, refineRatio, fine_active
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
                if (refineRatio>1 .and. nx==globalNx) then
                    if (.not.fine_active(i,j)) cycle
                endif

                m(0) = f(i, j, 0)+f(i, j, 1)+f(i, j, 2)+f(i, j, 3)+f(i, j, 4)+f(i, j, 5)+f(i, j, 6)+f(i, j, 7)+f(i, j, 8)
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
                f_post(i, j, 1) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+m_post(3)/6.0d0-m_post(4)/6.0d0 &
                    +m_post(7)/4.0d0
                f_post(i, j, 2) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    +m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)/4.0d0
                f_post(i, j, 3) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0-m_post(3)/6.0d0+m_post(4)/6.0d0 &
                    +m_post(7)/4.0d0
                f_post(i, j, 4) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    -m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)/4.0d0
                f_post(i, j, 5) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)/4.0d0
                f_post(i, j, 6) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)/4.0d0
                f_post(i, j, 7) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)/4.0d0
                f_post(i, j, 8) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)/4.0d0

            enddo
        enddo
        return
    end subroutine collision
    !===============================================================================================


    !===============================================================================================
    ! 子程序: streaming
    ! 作用: 流场分布函数迁移
    !===============================================================================================
    subroutine streaming(nx, ny, f, f_post)    !先迁移，再边界处理

        use commondata, only: ex, ey
        use commondata, only: globalNx=>nx, refineRatio, fine_active
        implicit none

        integer(kind=4), intent(in) :: nx, ny
        real(kind=8), intent(inout) :: f(nx, ny, 0:8), f_post(0:nx+1, 0:ny+1, 0:8)

        integer(kind=4) :: i, j
        integer(kind=4) :: ip, jp
        integer(kind=4) :: alpha

        !$acc parallel loop gang vector collapse(2) present(f, f_post, ex, ey) async(1) private(alpha, ip, jp)
        do j = 1, ny
            do i = 1, nx
                if (refineRatio>1 .and. nx==globalNx) then
                    if (.not.fine_active(i,j)) cycle
                endif
                do alpha = 0, 8    !上游格点索引：fα(i,j) <- f_postα(i-exα, j-eyα)
                    ip = i-ex(alpha)    !边界附近 (ip/jp 可能为 0 或 nx+1/ny+1)，需在 bounceback/周期边界处理中覆盖修正边界分布
                    jp = j-ey(alpha)    !ghost 层在初始化中为 0，保证不会出现未初始化垃圾值

                    f(i, j, alpha) = f_post(ip, jp, alpha)
                enddo
            enddo
        enddo
        return
    end subroutine streaming
    !===============================================================================================


    !===============================================================================================
    ! 子程序: bounceback
    ! 作用: 流场物理壁面及周期边界处理
    !===============================================================================================
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
    !===============================================================================================


    !===============================================================================================
    ! 子程序: macro
    ! 作用: 由流场分布函数计算密度和速度
    !===============================================================================================
    subroutine macro(nx, ny, f, rho, u, v, Fx, Fy)

        use commondata, only: globalNx=>nx, refineRatio, fine_active
        implicit none

        integer(kind=4), intent(in) :: nx, ny
        real(kind=8), intent(in) :: f(nx, ny, 0:8), Fx(nx, ny), Fy(nx, ny)
        real(kind=8), intent(inout) :: rho(nx, ny), u(nx, ny), v(nx, ny)

        integer(kind=4) :: i, j

        !$acc parallel loop gang vector collapse(2) present(f, rho, u, v, Fx, Fy) async(1)
        do j = 1, ny
            do i = 1, nx
                if (refineRatio>1 .and. nx==globalNx) then
                    if (.not.fine_active(i,j)) cycle
                endif
                rho(i, j) = f(i, j, 0)+f(i, j, 1)+f(i, j, 2)+f(i, j, 3)+f(i, j, 4)+f(i, j, 5)+f(i, j, 6)+f(i, j, 7)+f(i, j, 8)
                u(i, j) = ( f(i, j, 1)-f(i, j, 3)+f(i, j, 5)-f(i, j, 6)-f(i, j, 7)+f(i, j, 8)+0.5d0*Fx(i, j) )/rho(i, j)    !含力LBM的半步动量修正：rho*u = Σ f e + 0.5*F，对应Guo forcing的二阶定义
                v(i, j) = ( f(i, j, 2)-f(i, j, 4)+f(i, j, 5)+f(i, j, 6)-f(i, j, 7)-f(i, j, 8)+0.5d0*Fy(i, j) )/rho(i, j)
            enddo
        enddo
        return
    end subroutine macro
    !===============================================================================================


    !===============================================================================================
    ! 子程序: collisionT
    ! 作用: 原 D2Q5 MRT 温度场碰撞（无热流历史修正）
    !===============================================================================================
    subroutine collisionT(nx, ny, g, g_post, u, v, T, Qk, Qnu)

        use commondata, only: paraA
        use commondata, only: globalNx=>nx, refineRatio, fine_active
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
                if (refineRatio>1 .and. nx==globalNx) then
                    if (.not.fine_active(i,j)) cycle
                endif

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
    !===============================================================================================


    !===============================================================================================
    ! 子程序: streamingT
    ! 作用: 温度分布函数迁移
    !===============================================================================================
    subroutine streamingT(nx, ny, g, g_post)

        use commondata, only: ex, ey
        use commondata, only: globalNx=>nx, refineRatio, fine_active
        implicit none

        integer(kind=4), intent(in) :: nx, ny
        real(kind=8), intent(inout) :: g(nx, ny, 0:4), g_post(0:nx+1, 0:ny+1, 0:4)

        integer(kind=4) :: i, j
        integer(kind=4) :: ip, jp
        integer(kind=4) :: alpha

        !$acc parallel loop gang vector collapse(2) present(g, g_post, ex, ey) async(1) private(alpha, ip, jp)
        do j = 1, ny
            do i = 1, nx
                if (refineRatio>1 .and. nx==globalNx) then
                    if (.not.fine_active(i,j)) cycle
                endif
                do alpha = 0, 4
                    ip = i-ex(alpha)
                    jp = j-ey(alpha)

                    g(i, j, alpha) = g_post(ip, jp, alpha)
                enddo
            enddo
        enddo
        return
    end subroutine streamingT
    !===============================================================================================


    !===============================================================================================
    ! 子程序: bouncebackT
    ! 作用: 温度场恒温、绝热及周期边界处理
    !===============================================================================================
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
    !===============================================================================================


    !===============================================================================================
    ! 子程序: macroT
    ! 作用: 由温度分布函数计算温度
    !===============================================================================================
    subroutine macroT(nx, ny, g, T)

        use commondata, only: globalNx=>nx, refineRatio, fine_active
        implicit none

        integer(kind=4), intent(in) :: nx, ny
        real(kind=8), intent(in) :: g(nx, ny, 0:4)
        real(kind=8), intent(inout) :: T(nx, ny)

        integer(kind=4) :: i, j

        !$acc parallel loop gang vector collapse(2) present(g, T) async(1)
        do j = 1, ny
            do i = 1, nx
                if (refineRatio>1 .and. nx==globalNx) then
                    if (.not.fine_active(i,j)) cycle
                endif
                T(i, j) = g(i, j, 0)+g(i, j, 1)+g(i, j, 2)+g(i, j, 3)+g(i, j, 4)
            enddo
        enddo
        return
    end subroutine macroT
    !===============================================================================================


    !===============================================================================================
    ! 子程序: flow_moments
    ! 作用: 将流场分布函数变换到矩空间
    !===============================================================================================
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
    !===============================================================================================


    !===============================================================================================
    ! 子程序: flow_populations
    ! 作用: 由流场矩重建分布函数
    !===============================================================================================
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
    !===============================================================================================


    !===============================================================================================
    ! 子程序: thermal_moments
    ! 作用: 将温度分布函数变换到矩空间
    !===============================================================================================
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
    !===============================================================================================


    !===============================================================================================
    ! 子程序: thermal_populations
    ! 作用: 由温度矩重建分布函数
    !===============================================================================================
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
    !===============================================================================================


    !===============================================================================================
    ! 多块输出按不重叠物理分区积分，边缘节点使用裁剪权重；内部权重仍为 h^2。
    !===============================================================================================

    !===============================================================================================
    ! 子程序: calNuRe
    ! 作用: 计算面积加权的 Nu、Re 及场量统计
    !===============================================================================================
    subroutine calNuRe()

        use commondata, only: nx, ny, Thot, Tcold, lengthUnit, viscosity, diffusivity, timeUnit, historyFile, &
            nBlocks, itc, blockNi, blockNj, blockIlo, blockIhi, blockJlo, blockJhi, blockH, blockX0, blockY0, &
            blockOwnedBox, rho, u, v, T, dxWeight, dyWeight, ieee_is_finite
        use commondata, only: owned_cell_area, section_owned_weight
        implicit none

        integer(kind=4) :: b, i, j, k, jm, im
        real(kind=8) :: area, conv, vel2, mass, meanT, nu, re, hot, cold, middle, dTdx, dTdy, tm, um, vm, cellArea
        real(kind=8) :: w(4), dw(4)
        real(kind=8) :: tmin, tmax, rmin, rmax, scale, xmid, ymid, xlo, xhi, ylo, yhi, h

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
        do b = 1, nBlocks
            call select_block(b)
            h = blockH(b)
            do j = blockJlo(b), blockJhi(b)
                do i = blockIlo(b), blockIhi(b)
                    cellArea=owned_cell_area(b,i,j)
                    if (cellArea<=0d0) cycle
                    if (.not.ieee_is_finite(T(i, j)) .or. .not.ieee_is_finite(rho(i, j)) .or. &
                        .not.ieee_is_finite(u(i, j)) .or. .not.ieee_is_finite(v(i, j)) .or. rho(i, j) <= 0.0d0) then
                        write(*, *) 'Invalid state: coarse clock, block, i,j:', itc, b, i, j
                        error stop 'Nonfinite or nonpositive density in owned cells'
                    endif
                    cellArea = owned_cell_area(b,i,j)
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
            xlo = blockOwnedBox(1, b)
            xhi = blockOwnedBox(2, b)
            ylo = blockOwnedBox(3, b)
            yhi = blockOwnedBox(4, b)
#ifdef SideHeatedCell
            if (xlo == 0.0d0) then
                do j = blockJlo(b), blockJhi(b)
                    hot = hot+(8.0d0*Thot-9.0d0*T(1, j)+T(2, j))/(3.0d0*h)*dyWeight(j)/dble(ny)
                enddo
            endif
            if (xhi == dble(nx)) then
                do j = blockJlo(b), blockJhi(b)
                    cold = cold+(-8.0d0*Tcold+9.0d0*T(blockNi(b), j)-T(blockNi(b)-1, j))/(3.0d0*h)*dyWeight(j)/dble(ny)
                enddo
            endif
            if (xmid >= xlo .and. xmid < xhi) then
                call section_weights((xmid-blockX0(b))/h+0.5d0, blockNi(b), im, w, dw)
                do j = blockJlo(b), blockJhi(b)
                    if (section_owned_weight(b,j,1,xmid)<=0d0) cycle
                    tm = sum(w*T(im:im+3, j))
                    um = sum(w*u(im:im+3, j))
                    dTdx = sum(dw*T(im:im+3, j))/h
                    middle = middle+(um*tm/diffusivity-dTdx)*section_owned_weight(b,j,1,xmid)/dble(ny)
                enddo
            endif
#else
            if (ylo == 0.0d0) then
                do i = blockIlo(b), blockIhi(b)
                    hot = hot+(8.0d0*Thot-9.0d0*T(i, 1)+T(i, 2))/(3.0d0*h)*dxWeight(i)/dble(nx)
                enddo
            endif
            if (yhi == dble(ny)) then
                do i = blockIlo(b), blockIhi(b)
                    cold = cold+(-8.0d0*Tcold+9.0d0*T(i, blockNj(b))-T(i, blockNj(b)-1))/(3.0d0*h)*dxWeight(i)/dble(nx)
                enddo
            endif
            if (ymid >= ylo .and. ymid < yhi) then
                call section_weights((ymid-blockY0(b))/h+0.5d0, blockNj(b), jm, w, dw)
                do i = blockIlo(b), blockIhi(b)
                    if (section_owned_weight(b,i,2,ymid)<=0d0) cycle
                    tm = sum(w*T(i, jm:jm+3))
                    vm = sum(w*v(i, jm:jm+3))
                    dTdy = sum(dw*T(i, jm:jm+3))/h
                    middle = middle+(vm*tm/diffusivity-dTdy)*section_owned_weight(b,i,2,ymid)/dble(nx)
                enddo
            endif
#endif
        enddo
        ! 空间积分只覆盖 ownedBox；重叠节点按各自分区面积计权，不重复计算整块面积。
        ! Re 使用速度平方的面积平均再开方，非稳态时间统计也采用同一 RMS 定义。
        nu = 1.0d0+conv/area*scale/diffusivity
        re = sqrt(vel2/area)*lengthUnit/viscosity
        open(newunit = k, file = historyFile, status = 'old', position = 'append')
        write(k, '(12(ES24.16E3,1X))') dble(itc)/timeUnit, nu, re, scale*hot, scale*cold, scale*middle, &
            mass, meanT/area, tmin, tmax, rmin, rmax
        close(k)
        write(*, '(a,f12.5,a,es13.5,a,es13.5)') 't_ff=', dble(itc)/timeUnit, ' NuVolAvg=', nu, ' ReVolRMS=', re
    end subroutine calNuRe
    !===============================================================================================


    !===============================================================================================
    ! 子程序: check
    ! 作用: 计算稳态速度与温度的收敛误差
    !===============================================================================================
    subroutine check()
#ifdef steadyFlow

        use commondata, only: nBlocks, itc, errorU, errorT, blockIlo, blockIhi, blockJlo, blockJhi, u, v, T, &
            dxWeight, dyWeight
#ifdef steadyFlow
        use commondata, only: up, vp, Tp
#endif
        use commondata, only: owned_cell_area, section_owned_weight
        implicit none

        integer(kind=4) :: b, k, i, j
        real(kind=8) :: du, uu, dt, tt, cellArea

        du = 0.0d0
        uu = 0.0d0
        dt = 0.0d0
        tt = 0.0d0
        call update_host_all(.false.)
        do b = 1, nBlocks
            call select_block(b)
            do j = blockJlo(b), blockJhi(b)
                do i = blockIlo(b), blockIhi(b)
                    cellArea = owned_cell_area(b,i,j)
                    du = du+cellArea*((u(i, j)-up(i, j))**2+(v(i, j)-vp(i, j))**2)
                    uu = uu+cellArea*(u(i, j)**2+v(i, j)**2)
                    dt = dt+cellArea*(T(i, j)-Tp(i, j))**2
                    tt = tt+cellArea*T(i, j)**2
                enddo
            enddo
            up = u
            vp = v
            Tp = T
        enddo
        errorU = sqrt(du/max(uu, 1.0d-300))
        errorT = sqrt(dt/max(tt, 1.0d-300))
        open(newunit = k, file = 'Convergence_2DOpenaccMultiblock.dat', status = 'unknown', position = 'append')
        write(k, '(I12,2(1X,ES24.16E3))') itc, errorU, errorT
        close(k)
        write(*, *) 'errorU,errorT:', errorU, errorT
#endif
    end subroutine check
    !===============================================================================================


    !===============================================================================================
    ! 子程序: output_Tecplot
    ! 作用: 输出各统计分区的 Tecplot 数据
    !===============================================================================================
    subroutine output_Tecplot()

        use commondata, only: lengthUnit, timeUnit, pltFolderPrefix, nBlocks, itc, pltFileNum, blockIlo, &
            blockIhi, blockJlo, blockJhi, blockH, blockX0, blockY0, rho, u, v, T, dxWeight, dyWeight
        use commondata, only: owned_cell_area, section_owned_weight
        implicit none

        integer(kind=4) :: k, b, i, j
        character(16) :: num

        pltFileNum = pltFileNum+1
        write(num, '(I10.10)') pltFileNum
        call update_host_all(.false.)
        open(newunit = k, file = pltFolderPrefix//'-'//trim(num)//'.dat', status = 'replace')
        write(k, '(a)') 'VARIABLES="x/L","y/L","u","v","T","rho","h/L","integration_area/L^2"'
        do b = 1, nBlocks
            call select_block(b)
            write(k, '(a,I0,a,I0,a,I0,a,ES24.16E3)') 'ZONE T="block ', b, '", I=', blockIhi(b)-blockIlo(b)+1, &
                ', J=', blockJhi(b)-blockJlo(b)+1, ', F=POINT, SOLUTIONTIME=', dble(itc)/timeUnit
            do j = blockJlo(b), blockJhi(b)
                do i = blockIlo(b), blockIhi(b)
                    write(k, '(8(ES24.16E3,1X))') (blockX0(b)+(dble(i)-0.5d0)*blockH(b))/lengthUnit, &
                        (blockY0(b)+(dble(j)-0.5d0)*blockH(b))/lengthUnit, u(i, j), v(i, j), T(i, j), rho(i, j), &
                        blockH(b)/lengthUnit, &
                        owned_cell_area(b,i,j)/lengthUnit**2
                enddo
            enddo
        enddo
        close(k)
    end subroutine output_Tecplot
    !===============================================================================================


    !===============================================================================================
    ! 子程序: output_SnapshotFile
    ! 作用: 输出带积分权重的多块快照
    !===============================================================================================
    subroutine output_SnapshotFile()

        use commondata, only: nx, ny, lengthUnit, timeUnit, snapshotFilePrefix, snapshotMagic, nBlocks, itc, &
            snapshotFileNum, blockIlo, blockIhi, blockJlo, blockJhi, blockH, blockX0, blockY0, blockOwnedBox, &
            rho, u, v, T, dxWeight, dyWeight
        use commondata, only: owned_cell_area, section_owned_weight
        implicit none

        integer(kind=4) :: k, b, i,j
        character(16) :: num

        snapshotFileNum = snapshotFileNum+1
        write(num, '(I10.10)') snapshotFileNum
        call update_host_all(.false.)
        open(newunit = k, file = snapshotFilePrefix//'-'//trim(num)//'.bin', form = 'unformatted', &
            access = 'stream', status = 'replace')
        write(k) snapshotMagic, nBlocks, nx, ny, itc, dble(itc)/timeUnit, lengthUnit
        do b = 1, nBlocks
            call select_block(b)
            write(k) blockIhi(b)-blockIlo(b)+1, blockJhi(b)-blockJlo(b)+1, &
                blockX0(b)+(dble(blockIlo(b))-0.5d0)*blockH(b), blockY0(b)+(dble(blockJlo(b))-0.5d0)*blockH(b), &
                blockH(b), blockOwnedBox(:, b)
            write(k) dxWeight(blockIlo(b):blockIhi(b)), dyWeight(blockJlo(b):blockJhi(b))
            ! Snapshot v3 appends explicit 2D area after the separable coordinate weights.
            write(k) ((owned_cell_area(b,i,j),i=blockIlo(b),blockIhi(b)),j=blockJlo(b),blockJhi(b))
            write(k) u(blockIlo(b):blockIhi(b), blockJlo(b):blockJhi(b)), v(blockIlo(b):blockIhi(b), blockJlo(b):blockJhi(b)), &
                T(blockIlo(b):blockIhi(b), blockJlo(b):blockJhi(b)), rho(blockIlo(b):blockIhi(b), blockJlo(b):blockJhi(b))
        enddo
        close(k)
    end subroutine output_SnapshotFile
    !===============================================================================================


    !===============================================================================================
    ! 子程序: model_signature
    ! 作用: 生成重启模型配置标识
    !===============================================================================================
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
    !===============================================================================================


    !===============================================================================================
    ! 子程序: physical_signature
    ! 作用: 生成重启物理参数与输出间隔标识
    !===============================================================================================
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
    !===============================================================================================


    !===============================================================================================
    ! 完整检查点：仅在粗细同步时写入

    !===============================================================================================
    ! 子程序: output_ReloadFile
    ! 作用: 在同步时刻保存完整检查点
    !===============================================================================================
    subroutine output_ReloadFile()

        use commondata, only: nx, ny, refineRatio, fineLayerCellsLeft, fineLayerCellsRight, &
            fineLayerCellsBottom, fineLayerCellsTop, &
            overlapCells, reloadFileNum, reloadFilePrefix, restartMagic, nBlocks, itc, snapshotFileNum, &
            pltFileNum, nextSample, nextReload, nextPlt, errorU, errorT, blockNi, blockNj, blockNh, blockIlo, &
            blockIhi, blockJlo, blockJhi, blockH, blockX0, blockY0, blockOwnedBox, f, g, rho, u, v, T, Fx, &
            Fy, p
#ifdef steadyFlow
        use commondata, only: up, vp, Tp
#endif
        implicit none

        integer(kind=4) :: k, b, currentModel(12)
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
            fineLayerCellsBottom, fineLayerCellsTop, overlapCells, nBlocks
        call model_signature(currentModel)
        call physical_signature(currentPhysics)
        write(k) currentModel, currentPhysics
        write(k) itc, nextSample, nextReload, nextPlt, snapshotFileNum, pltFileNum, reloadFileNum, errorU, errorT
        do b = 1, nBlocks
            call select_block(b)
            write(k) blockNi(b), blockNj(b), blockIlo(b), blockIhi(b), blockJlo(b), blockJhi(b), blockNh(b), &
                blockX0(b), blockY0(b), blockH(b), blockOwnedBox(:, b)
            write(k) f, g, u, v, T, rho, Fx, Fy, p
#ifdef steadyFlow
            write(k) up, vp, Tp
#endif
        enddo
        close(k)
        ! 完整写完独立编号的 checkpoint 后再更新 latest 指针，旧 checkpoint 仍保留。
        open(newunit = k, file = reloadFilePrefix//'-latest.meta', status = 'replace')
        write(k, '(a)') trim(name)
        close(k)
    end subroutine output_ReloadFile
    !===============================================================================================


    !===============================================================================================
    ! 精确续算：检查网格、物理参数、统计分区与历史状态

    !===============================================================================================
    ! 子程序: read_restart
    ! 作用: 核对配置并恢复完整检查点状态
    !===============================================================================================
    subroutine read_restart()

        use commondata, only: nx, ny, refineRatio, fineLayerCellsLeft, fineLayerCellsRight, &
            fineLayerCellsBottom, fineLayerCellsTop, &
            overlapCells, reloadFileNum, reloadFilePrefix, restartMagic, nBlocks, itc, snapshotFileNum, &
            pltFileNum, nextSample, nextReload, nextPlt, errorU, errorT, blockNi, blockNj, blockNh, blockIlo, &
            blockIhi, blockJlo, blockJhi, blockH, blockX0, blockY0, blockOwnedBox, f, g, rho, u, v, T, Fx, &
            Fy, p
#ifdef steadyFlow
        use commondata, only: up, vp, Tp
#endif
        implicit none

        integer(kind=4) :: k, b, ios, head(9), geom(7), sig(12), currentModel(12)
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
        if (magic /= restartMagic) error stop 'Wrong checkpoint format'
        if (any(head /= [nx, ny, refineRatio, fineLayerCellsLeft, fineLayerCellsRight, &
            fineLayerCellsBottom, fineLayerCellsTop, overlapCells, nBlocks])) &
            error stop 'Restart mesh/refinement mismatch'
        read(k) sig, phys
        call model_signature(currentModel)
        call physical_signature(currentPhysics)
        if (any(sig /= currentModel) .or. any(phys /= currentPhysics)) &
            error stop 'Restart model/physics/output-cadence mismatch'
        read(k) itc, nextSample, nextReload, nextPlt, snapshotFileNum, pltFileNum, reloadFileNum, errorU, errorT
        if (itc < 0 .or. mod(itc, refineRatio) /= 0) error stop 'Restart is not at a synchronized time'
        do b = 1, nBlocks
            call select_block(b)
            read(k) geom, coord
            if (any(geom /= [blockNi(b), blockNj(b), blockIlo(b), blockIhi(b), blockJlo(b), blockJhi(b), blockNh(b)]) .or. &
                any(coord /= [blockX0(b), blockY0(b), blockH(b), blockOwnedBox(:, &
                b)])) error stop 'Restart block layout mismatch'
            read(k, iostat = ios) f, g, u, v, T, rho, Fx, Fy, p
            if (ios /= 0) error stop 'Incomplete checkpoint state/history'
#ifdef steadyFlow
            read(k) up, vp, Tp
#endif
        enddo
        close(k)
    end subroutine read_restart
    !===============================================================================================


    !===============================================================================================
    ! 函数: scheduled_step
    ! 作用: 将输出时间向上对齐到粗细同步步
    !===============================================================================================
    integer(kind=4) function scheduled_step(index, interval) result(step)

        use commondata, only: refineRatio, timeUnit
        implicit none

        integer(kind=4), intent(in) :: index
        real(kind=8), intent(in) :: interval

        step = max(refineRatio, ceiling(dble(index)*interval*timeUnit/dble(refineRatio))*refineRatio)
    end function scheduled_step
    !===============================================================================================


    !===============================================================================================
    ! 子程序: check_history
    ! 作用: 核对检查点对应的统计历史文件
    !===============================================================================================
    subroutine check_history()

        use commondata, only: timeUnit, outputSnapshotInterval, historyFile, nextSample, scheduled_step, &
            ieee_is_finite
        implicit none

        integer(kind=4) :: k, ios, n
        real(kind=8) :: values(12), lastTime, expected
        character(512) :: line

        lastTime = -1.0d0
        n = 0
        open(newunit = k, file = historyFile, status = 'old', iostat = ios)
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
            if (abs(lastTime-expected) > 1.0d-10*max(1.0d0, expected)) error stop 'Checkpoint/history time mismatch'
        endif
    end subroutine check_history
    !===============================================================================================


    !===============================================================================================
    ! 子程序: average_window
    ! 作用: 对指定时间窗积分并计算平均量
    !===============================================================================================
    subroutine average_window(t0, t1, result, coverage)

        use commondata, only: historyFile, ieee_is_finite
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
        open(newunit = k, file = historyFile, status = 'old')
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
    !===============================================================================================


    !===============================================================================================
    ! 子程序: output_unsteady_NuRe_postprocess
    ! 作用: 输出完整时间窗与前后半窗的统计结果
    !===============================================================================================
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
    !===============================================================================================

    ! 细环只保留中心空区之外的节点；四个原细区的拼接处没有特殊处理。
    logical function fine_active(i,j)
        use commondata, only: nx,ny,refineRatio,overlapCells,fineLayerCellsLeft,fineLayerCellsRight, &
            fineLayerCellsBottom,fineLayerCellsTop
        implicit none
        !$acc routine seq
        integer,intent(in) :: i,j
        real(8) :: x,y,o
        x=dble(i)-.5d0; y=dble(j)-.5d0; o=dble(overlapCells*refineRatio)
        fine_active=refineRatio==1 .or. x<=fineLayerCellsLeft-.5d0+o .or. &
            x>=nx-fineLayerCellsRight+.5d0-o .or. y<=fineLayerCellsBottom-.5d0+o .or. &
            y>=ny-fineLayerCellsTop+.5d0-o
    end function fine_active

    real(8) function owned_cell_area(b,i,j) result(a)
        use commondata, only: nBlocks,blockOwnedBox,blockH,blockX0,blockY0
        implicit none
        integer,intent(in) :: b,i,j
        real(8) :: x,y,h,wx,wy
        h=blockH(b); x=blockX0(b)+(i-.5d0)*h; y=blockY0(b)+(j-.5d0)*h
        wx=max(0d0,min(x+h/2,blockOwnedBox(2,b))-max(x-h/2,blockOwnedBox(1,b)))
        wy=max(0d0,min(y+h/2,blockOwnedBox(4,b))-max(y-h/2,blockOwnedBox(3,b)))
        a=wx*wy
        if(nBlocks>1 .and. b==2) then
            wx=max(0d0,min(x+h/2,blockOwnedBox(2,1))-max(x-h/2,blockOwnedBox(1,1)))
            wy=max(0d0,min(y+h/2,blockOwnedBox(4,1))-max(y-h/2,blockOwnedBox(3,1)))
            a=a-wx*wy
        endif
    end function owned_cell_area

    real(8) function section_owned_weight(b,k,axis,position) result(w)
        use commondata, only: nBlocks,blockOwnedBox,blockH,blockX0,blockY0
        implicit none
        integer,intent(in) :: b,k,axis
        real(8),intent(in) :: position
        real(8) :: q,h,lo,hi
        h=blockH(b)
        if(axis==1) then
            q=blockY0(b)+(k-.5d0)*h; lo=blockOwnedBox(3,b);hi=blockOwnedBox(4,b)
        else
            q=blockX0(b)+(k-.5d0)*h; lo=blockOwnedBox(1,b);hi=blockOwnedBox(2,b)
        endif
        w=max(0d0,min(q+h/2,hi)-max(q-h/2,lo))
        if(nBlocks>1 .and. b==2) then
            if(axis==1) then
                if(position<blockOwnedBox(1,1) .or. position>=blockOwnedBox(2,1)) return
                lo=blockOwnedBox(3,1);hi=blockOwnedBox(4,1)
            else
                if(position<blockOwnedBox(3,1) .or. position>=blockOwnedBox(4,1)) return
                lo=blockOwnedBox(1,1);hi=blockOwnedBox(2,1)
            endif
            w=w-max(0d0,min(q+h/2,hi)-max(q-h/2,lo))
        endif
    end function section_owned_weight
