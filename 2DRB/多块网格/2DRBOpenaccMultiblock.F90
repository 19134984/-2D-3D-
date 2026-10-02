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
        ! ieee_arithmetic 是 Fortran 标准内置模块。
        ! intrinsic 明确指定使用编译器提供的内置模块。
        ! 有限值返回 .true.，NaN 或正负无穷大返回 .false.。
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        implicit none

        !===============================================================================================
        ! 网格参数
        !===============================================================================================
        ! nx,ny 为按细格距计算的全域格子长度
        ! 物理壁面在 x=0、nx 和 y=0、ny；最细块 dx=dt=1，中心粗块 dx=dt=refineRatio。
        ! 默认交界面为 x/y=127.5、895.5；粗细块共用这些基准边界，积分也用这些基准边界。
        ! 中心宽高 nx-Left-Right+1、ny-Bottom-Top+1 必须能被 refineRatio 整除。
        ! 右=左+1、上=下+1 时，中心宽高为 nx-2*Left、ny-2*Bottom；初始化时检查整除。
        ! 粗细比和四个分区位置
        integer(kind=4), parameter :: nx=1024, ny=1024    ! 全域格子长度，必须能被 refineRatio 整除
        integer(kind=4), parameter :: refineRatio=2    ! 粗/细格距之比，取2、4、8等2的幂
        integer(kind=4), parameter :: fineLayerCellsLeft=128    ! 从左墙向内的细节点编号，x = Left-0.5
        integer(kind=4), parameter :: fineLayerCellsRight=129    ! 从右墙向内的细节点编号，x = nx-Right+0.5
        integer(kind=4), parameter :: fineLayerCellsBottom=128    ! 从下墙向内的细节点编号，y = Bottom-0.5
        integer(kind=4), parameter :: fineLayerCellsTop=129    ! 从上墙向内的细节点编号，y = ny-Top+0.5
        !===============================================================================================

        !===============================================================================================
        ! 是否在计算前从旧算例重启
        !===============================================================================================
        integer(kind=4), parameter :: loadInitField=0   ! 0: 不重启；1: 按 latest .meta 自动续算

        ! 正常断电续算只需要设置 loadInitField=1；
        ! 代码只读取 <reloadFilePrefix>-latest.meta，并从里面找到最新的 .bin。
        ! 正常续算不用改 reloadFileName；只有 latest .meta 缺失时，
        ! 才手动填写已有检查点的完整文件名，旧编号文件也可以直接填写。
        character(len=256) :: reloadFileName=""        ! latest .meta 缺失时使用，包含 .bin 后缀

        ! 多块检查点的步数和时间信息保存在 .bin 中；.meta 只记录文件名，无需按编号推断。
        real(kind=8) :: reloadDimensionlessTime=0.0d0   ! 续算前已累计的 t_ff；由检查点累计步数换算
        integer(kind=4) :: reloadItcOffset=0           ! 续算前已累计的格子步数；从检查点读取
        logical :: reloadMetadataLoaded=.false.         ! 标记是否成功读取 reload 元数据
        !===============================================================================================

        !===============================================================================================
        ! 无量纲参数
        !===============================================================================================
#ifdef SideHeatedCell
        real(kind=8), parameter :: lengthUnit=dble(nx)    !侧壁差温：特征长度取 x 方向长度
#else
        real(kind=8), parameter :: lengthUnit=dble(ny)    !上下差温：特征长度取 y 方向长度
#endif
        real(kind=8), parameter :: pi=acos(-1.0d0)

        real(kind=8), parameter :: Rayleigh=1.0d7
        real(kind=8), parameter :: Prandtl=0.7d0
        real(kind=8), parameter :: Mach=0.1d0
        real(kind=8), parameter :: Thot=0.5d0, Tcold=-0.5d0
        real(kind=8), parameter :: Tref=0.5d0*(Thot+Tcold)
        ! 多块网格：以下输运参数采用细格子单位。
        real(kind=8), parameter :: tauf=0.5d0+Mach*lengthUnit*sqrt(3.0d0*Prandtl/Rayleigh)
        real(kind=8), parameter :: viscosity=(tauf-0.5d0)/3.0d0
        real(kind=8), parameter :: diffusivity=viscosity/Prandtl

        ! 多块快照当前保存格子速度，尚未使用下面的后处理比较标度。
        real(kind=8), parameter :: velocityScaleCompare=lengthUnit/diffusivity    ! 默认采用热扩散标度 UL/kappa；若要按自由落体标度比较，可改为 1.0d0/velocityUnit

#ifdef  SideHeatedHa
        real(kind=8), parameter :: Ha=20.0d0    !磁场强度
        real(kind=8), parameter :: phi=0.0d0*pi/180.0d0    !磁场角度，以水平向右为0，修改0.0d0即可
        real(kind=8), parameter :: B2sigemarho=Ha**2*viscosity/lengthUnit**2    !动量方程上的源项系数
#endif

        ! 高阶矩参数修正
        real(kind=8), parameter :: paraA=20.0d0*sqrt(3.0d0)*diffusivity-4.0d0

        ! 浮力项参数，采用细格子单位。
        real(kind=8), parameter :: gBeta1=Rayleigh*viscosity*diffusivity/lengthUnit
        real(kind=8), parameter :: gBeta=gBeta1/lengthUnit/lengthUnit    !gbetaΔT

        real(kind=8), parameter :: timeUnit=sqrt(lengthUnit/gBeta)    !无量纲时间
        real(kind=8), parameter :: velocityUnit=sqrt(gBeta*lengthUnit)    !无量纲速度

        real(kind=8), parameter :: Snu=1.0d0/tauf, Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)    !动量的多松弛系数

        ! 多块网格使用Luo的温度算法。
        real(kind=8), parameter :: Qk=3.0d0-sqrt(3.0d0), Qnu=4.0d0*sqrt(3.0d0)-6.0d0    !温度算法的多松弛系数
        real(kind=8), parameter :: thermalGeqCoeff=10.0d0/(4.0d0+paraA)
        !===============================================================================================

        !===============================================================================================
        ! 输出/备份相关设置（以自由落体时间 t_ff 为单位）
        !===============================================================================================
        real(kind=8), parameter :: epsU=1.0d-7, epsT=1.0d-7    ! 稳态收敛阈值

#ifdef steadyFlow
        real(kind=8), parameter :: outputSnapshotInterval=10.0d0    ! 快照和 Nu/Re 时间序列采样间隔（单位：t_ff）
        real(kind=8), parameter :: reloadFileInterval=100.0d0    ! f/g 重启文件输出间隔（单位：t_ff）
        real(kind=8), parameter :: outputPltFileInterval=100.0d0    ! Tecplot 文件输出间隔（单位：t_ff）
        integer(kind=4), parameter :: outputSnapshotFile=1    ! 是否输出后处理快照文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputPltFile=1    ! 是否输出 plt 文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputReloadFile=1    ! 是否周期输出 f/g 重启文件：0=不输出，1=输出
        integer(kind=4), parameter :: itc_max=20000000    ! 稳态：最大格子步，实际仍由 errorU/errorT 提前停止
#endif

#ifdef unsteadyFlow
        real(kind=8), parameter :: outputSnapshotInterval=0.5d0    ! uvTrho 快照和 Nu/Re 时间序列采样间隔（单位：t_ff）
        real(kind=8), parameter :: reloadFileInterval=100.0d0    ! f/g 重启文件输出间隔（单位：t_ff）
        real(kind=8), parameter :: outputPltFileInterval=100.0d0    ! Tecplot 文件周期输出间隔（单位：t_ff）
        real(kind=8), parameter :: unsteadyRunDuration=1000.0d0    ! 非稳态总目标时长，续算时只补足到该 t_ff
        ! 以下三个参数只控制非稳态结束后的 Nu/Re 统计平均窗口，不改变推进时长或采样频率。
        ! 时间以整个非稳态算例的绝对 t_ff 计，续算统计会包含旧文件中已有的历史数据。
        real(kind=8), parameter :: unsteadyAverageStartTf=0.5d0*unsteadyRunDuration    ! 平均窗口起点
        real(kind=8), parameter :: unsteadyAverageEndTf=unsteadyRunDuration    ! 平均窗口终点
        real(kind=8), parameter :: unsteadyAverageMidTf=0.5d0*(unsteadyAverageStartTf+unsteadyAverageEndTf)    ! 前/后半分界
        integer(kind=4), parameter :: outputSnapshotFile=1    ! 是否输出后处理快照文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputPltFile=1    ! 是否输出 plt 文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputReloadFile=1    ! 是否周期输出 f/g 重启文件：0=不输出，1=输出
        integer(kind=4), parameter :: itc_max=max(1,ceiling(unsteadyRunDuration*timeUnit))    ! 非稳态：由总 t_ff 自动换算格子步
#endif

        ! 最终输出时刻以粗步为准，只有粗步下才有完整的全场数据，各间隔至少为一个粗步。

        character(len=*), parameter :: snapshotFilePrefix="buoyancyCavity2DOpenaccMultiblockSnapshot"
        ! 快照输出文件前缀（实际文件名形如：<snapshotFilePrefix>-<累计t_ff>tff.bin）

        character(len=*), parameter :: pltFolderPrefix="buoyancyCavity2DOpenaccMultiblockTecplot"
        ! Tecplot 输出文件前缀（实际文件名形如：<pltFolderPrefix>-<累计t_ff>tff.dat）

        character(len=*), parameter :: reloadFilePrefix="reloadFile2DOpenaccMultiblock"
        ! 重启文件前缀；latest meta 模式读取 meta 中记录的文件名，新文件以累计 t_ff 命名

        character(len=*), parameter :: settingsFile="SimulationSettings2DOpenaccMultiblock.txt"
        character(*), parameter :: NuReHistoryFile='NuRe_2DOpenaccMultiblock.dat'    ! 按采样时间记录 Nu、Re 及场量统计，续算时追加
        !===============================================================================================

        !===============================================================================================
        ! 格子方向
        !===============================================================================================
        ! 接收层取2层。
        ! 粗、细网格各延伸 interfaceSkin 个本地格距。
        ! 换成细格距单位时，粗网格延伸宽度为 interfaceSkin*refineRatio，细网格为 interfaceSkin。
        integer(kind=4), parameter :: interfaceSkin=2

        integer(kind=4) :: ex(0:8)=[0,1,0,-1,0,1,-1,-1,1]
        integer(kind=4) :: ey(0:8)=[0,0,1,0,-1,1,1,-1,-1]
        real(kind=8) :: omega(0:8), omegaT(0:4)

        ! itc 为本次运行已完成的细时间步数，每次启动从 0 开始。
        ! reloadItcOffset 为续算前累计细步数；总步数为 reloadItcOffset+itc。
        ! 粗细比为 r 时，每个粗步使 itc 增加 r。
        integer(kind=4) :: itc=0

        ! NuReSampleCount 是已完成的 Nu/Re 采样次数，不是物理时间。
        ! 初值为 0，每次 calNuRe 增加 1；续算时恢复累计次数，实际 t_ff 为 (reloadItcOffset+itc)/timeUnit。
        integer(kind=4) :: NuReSampleCount=0

        ! 以下三个量分别为 Nu/Re 与快照、续算文件、Tecplot 的固定输出间隔，单位为细格子步。
        ! 间隔须向上对齐到 refineRatio 的整数倍，保证输出时粗细网格的时间是同步的。
        integer(kind=4) :: outputSnapshotIntervalItc
        integer(kind=4) :: reloadFileIntervalItc, outputPltFileIntervalItc
#ifdef steadyFlow
        integer(kind=4) :: steadyCheckIntervalItc    ! 稳态收敛检查及监测间隔，单位为细格子步；将 2000 步向上对齐到粗细同步步。
#endif
        real(kind=8) :: errorU=100.0d0, errorT=100.0d0
        !===============================================================================================

        !===============================================================================================
        ! 粗网格相应的信息。
        !===============================================================================================
        ! 左右细区贯穿全高，包含四角区域；上下细区仅存左右细区之间的节点。
        ! 各场量的 i,j 索引编号从 1 开始，物理坐标为 xOffset+(i-0.5)*dx、yOffset+(j-0.5)*dx。
        ! f_post/g_post 另留 0 和 n+1 外围位置，用于相邻细区的迁移数据交换。
        real(kind=8) :: dxCoarse
        real(kind=8), parameter :: dxFine=1.0d0
        real(kind=8) :: SnuCoarse, SqCoarse, QkCoarse, QnuCoarse

        ! 粗网格时间插值的最大时间层下标，多块模式设为 2，对应 0、1、2 三层。
        ! 一个粗步内：0=上一粗时刻，1=本轮粗步起点，2=本轮粗步预测终点。
        ! 各层保存对应时刻的流场交换量；首轮没有上一粗时刻的数据，第 0 层权重取 0。
        ! refineRatio=1 时不做粗细时间插值，设为 0，数组仅分配 0:0 一层。
        integer(kind=4) :: timeInterpMaxIndexCoarse

        ! 中心粗网格负责积分的矩形边界，依次为 xmin,xmax,ymin,ymax。
        ! 默认其值为 [127.5,895.5,127.5,895.5]。不包含重叠层
        real(kind=8) :: centerIntegrationBox(4)

        ! 积分严格从分区线 x=127.5 到 895.5、y=127.5 到 895.5（默认参数）。不包含重叠层
        ! 默认 dxCoarse=2，首末粗节点正好在分区线上：它们的积分权重各为 1，内部节点为 2。
        ! 积分下标直接去掉两侧 interfaceSkin 层，单网格则使用全部节点；不再单独保存起止编号。

        ! 粗网格坐标偏移
        real(kind=8) :: xOffsetCoarse, yOffsetCoarse

        logical :: coarseIsWall(4)    ! 粗网格左、右、下、上是否执行物理域边界条件；周期边界也为真。
        real(kind=8), allocatable :: f_coarse(:,:,:)    ! 粗网格 D2Q9 流场分布函数；第三维 0:8 为速度方向。
        real(kind=8), allocatable :: g_coarse(:,:,:)    ! 粗网格 D2Q5 温度分布函数；第三维 0:4 为速度方向。
        real(kind=8), allocatable :: f_post_coarse(:,:,:)    ! 碰撞后的流场分布函数；外围 0、n+1 供迁移。
        real(kind=8), allocatable :: g_post_coarse(:,:,:)    ! 碰撞后的温度分布函数；外围 0、n+1 供迁移。

        real(kind=8), allocatable :: rho_coarse(:,:)    ! 粗网格当前密度。
        real(kind=8), allocatable :: u_coarse(:,:)    ! 粗网格当前 x 方向速度。
        real(kind=8), allocatable :: v_coarse(:,:)    ! 粗网格当前 y 方向速度。
        real(kind=8), allocatable :: T_coarse(:,:)    ! 粗网格当前温度。
        real(kind=8), allocatable :: Fx_coarse(:,:)    ! Fx_coarse(i,j) = 0.0d0。
        real(kind=8), allocatable :: Fy_coarse(:,:)    ! Fy_coarse(i,j) = rho_coarse(i,j)*(dxCoarse * gBeta)*(T_coarse(i,j)-Tref)。

#ifdef steadyFlow
        real(kind=8), allocatable :: up_coarse(:,:)    ! 上一次稳态检查时的 x 方向速度。
        real(kind=8), allocatable :: vp_coarse(:,:)    ! 上一次稳态检查时的 y 方向速度。
        real(kind=8), allocatable :: Tp_coarse(:,:)    ! 上一次稳态检查时的温度。
#endif
        !===============================================================================================

        !===============================================================================================
        ! 左边细网格相应的信息。
        !===============================================================================================
        ! 四个细区按各自数组尺寸遍历；Offset 保存坐标偏移。
        ! 细区积分范围由中心矩形四边确定，节点面积在积分时计算。

        real(kind=8) :: xOffsetLeft, yOffsetLeft    ! 本地数组的坐标偏移。
        logical :: leftIsWall(4)    ! 左、右、下、上是否为物理墙面。
        real(kind=8), allocatable :: f_left(:,:,:)    ! 流场分布函数，第三维为 D2Q9 方向编号 0:8。
        real(kind=8), allocatable :: g_left(:,:,:)    ! 温度分布函数，第三维为 D2Q5 方向编号 0:4。
        real(kind=8), allocatable :: f_post_left(:,:,:)    ! 碰撞后的流场分布函数；前两维另留 0、n+1 外围位置供细区迁移交换。
        real(kind=8), allocatable :: g_post_left(:,:,:)    ! 碰撞后的温度分布函数；前两维另留 0、n+1 外围位置供细区迁移交换。

        real(kind=8), allocatable :: rho_left(:,:)    ! 当前密度。
        real(kind=8), allocatable :: u_left(:,:)    ! 当前 x 方向速度。
        real(kind=8), allocatable :: v_left(:,:)    ! 当前 y 方向速度。
        real(kind=8), allocatable :: T_left(:,:)    ! 当前温度。
        real(kind=8), allocatable :: Fx_left(:,:)    ! x 方向力源项；为 0。
        real(kind=8), allocatable :: Fy_left(:,:)    ! y 方向力源项；为 rho_left*gBeta*(T_left-Tref)。

#ifdef steadyFlow
        real(kind=8), allocatable :: up_left(:,:)    ! 上一次稳态检查时的 x 方向速度，用于计算收敛误差。
        real(kind=8), allocatable :: vp_left(:,:)    ! 上一次稳态检查时的 y 方向速度，用于计算收敛误差。
        real(kind=8), allocatable :: Tp_left(:,:)    ! 上一次稳态检查时的温度，用于计算收敛误差。
#endif
        !===============================================================================================

        !===============================================================================================
        ! 右边细网格相应的信息。
        !===============================================================================================

        real(kind=8) :: xOffsetRight, yOffsetRight    ! 本地数组的坐标偏移。
        logical :: rightIsWall(4)    ! 左、右、下、上是否为物理墙面。
        real(kind=8), allocatable :: f_right(:,:,:)    ! 流场分布函数，第三维为 D2Q9 方向编号 0:8。
        real(kind=8), allocatable :: g_right(:,:,:)    ! 温度分布函数，第三维为 D2Q5 方向编号 0:4。
        real(kind=8), allocatable :: f_post_right(:,:,:)    ! 碰撞后的流场分布函数；前两维另留 0、n+1 外围位置供细区迁移交换。
        real(kind=8), allocatable :: g_post_right(:,:,:)    ! 碰撞后的温度分布函数；前两维另留 0、n+1 外围位置供细区迁移交换。

        real(kind=8), allocatable :: rho_right(:,:)    ! 当前密度。
        real(kind=8), allocatable :: u_right(:,:)    ! 当前 x 方向速度。
        real(kind=8), allocatable :: v_right(:,:)    ! 当前 y 方向速度。
        real(kind=8), allocatable :: T_right(:,:)    ! 当前温度。
        real(kind=8), allocatable :: Fx_right(:,:)    ! x 方向力源项；为 0。
        real(kind=8), allocatable :: Fy_right(:,:)    ! y 方向力源项；为 rho_right*gBeta*(T_right-Tref)。

#ifdef steadyFlow
        real(kind=8), allocatable :: up_right(:,:)    ! 上一次稳态检查时的 x 方向速度，用于计算收敛误差。
        real(kind=8), allocatable :: vp_right(:,:)    ! 上一次稳态检查时的 y 方向速度，用于计算收敛误差。
        real(kind=8), allocatable :: Tp_right(:,:)    ! 上一次稳态检查时的温度，用于计算收敛误差。
#endif
        !===============================================================================================

        !===============================================================================================
        ! 下边细网格相应的信息。
        !===============================================================================================

        real(kind=8) :: xOffsetBottom, yOffsetBottom    ! 本地数组的坐标偏移。
        logical :: bottomIsWall(4)    ! 左、右、下、上是否为物理墙面。
        real(kind=8), allocatable :: f_bottom(:,:,:)    ! 流场分布函数，第三维为 D2Q9 方向编号 0:8。
        real(kind=8), allocatable :: g_bottom(:,:,:)    ! 温度分布函数，第三维为 D2Q5 方向编号 0:4。
        real(kind=8), allocatable :: f_post_bottom(:,:,:)    ! 碰撞后的流场分布函数；前两维另留 0、n+1 外围位置供细区迁移交换。
        real(kind=8), allocatable :: g_post_bottom(:,:,:)    ! 碰撞后的温度分布函数；前两维另留 0、n+1 外围位置供细区迁移交换。

        real(kind=8), allocatable :: rho_bottom(:,:)    ! 当前密度。
        real(kind=8), allocatable :: u_bottom(:,:)    ! 当前 x 方向速度。
        real(kind=8), allocatable :: v_bottom(:,:)    ! 当前 y 方向速度。
        real(kind=8), allocatable :: T_bottom(:,:)    ! 当前温度。
        real(kind=8), allocatable :: Fx_bottom(:,:)    ! x 方向力源项；为 0。
        real(kind=8), allocatable :: Fy_bottom(:,:)    ! y 方向力源项；为 rho_bottom*gBeta*(T_bottom-Tref)。

#ifdef steadyFlow
        real(kind=8), allocatable :: up_bottom(:,:)    ! 上一次稳态检查时的 x 方向速度，用于计算收敛误差。
        real(kind=8), allocatable :: vp_bottom(:,:)    ! 上一次稳态检查时的 y 方向速度，用于计算收敛误差。
        real(kind=8), allocatable :: Tp_bottom(:,:)    ! 上一次稳态检查时的温度，用于计算收敛误差。
#endif
        !===============================================================================================

        !===============================================================================================
        ! 上边细网格相应的信息。
        !===============================================================================================

        real(kind=8) :: xOffsetTop, yOffsetTop    ! 本地数组的坐标偏移。
        logical :: topIsWall(4)    ! 左、右、下、上是否为物理墙面。
        real(kind=8), allocatable :: f_top(:,:,:)    ! 流场分布函数，第三维为 D2Q9 方向编号 0:8。
        real(kind=8), allocatable :: g_top(:,:,:)    ! 温度分布函数，第三维为 D2Q5 方向编号 0:4。
        real(kind=8), allocatable :: f_post_top(:,:,:)    ! 碰撞后的流场分布函数；前两维另留 0、n+1 外围位置供细区迁移交换。
        real(kind=8), allocatable :: g_post_top(:,:,:)    ! 碰撞后的温度分布函数；前两维另留 0、n+1 外围位置供细区迁移交换。

        real(kind=8), allocatable :: rho_top(:,:)    ! 当前密度。
        real(kind=8), allocatable :: u_top(:,:)    ! 当前 x 方向速度。
        real(kind=8), allocatable :: v_top(:,:)    ! 当前 y 方向速度。
        real(kind=8), allocatable :: T_top(:,:)    ! 当前温度。
        real(kind=8), allocatable :: Fx_top(:,:)    ! x 方向力源项；为 0。
        real(kind=8), allocatable :: Fy_top(:,:)    ! y 方向力源项；为 rho_top*gBeta*(T_top-Tref)。

#ifdef steadyFlow
        real(kind=8), allocatable :: up_top(:,:)    ! 上一次稳态检查时的 x 方向速度，用于计算收敛误差。
        real(kind=8), allocatable :: vp_top(:,:)    ! 上一次稳态检查时的 y 方向速度，用于计算收敛误差。
        real(kind=8), allocatable :: Tp_top(:,:)    ! 上一次稳态检查时的温度，用于计算收敛误差。
#endif
        !===============================================================================================

        !===============================================================================================
        ! 重叠层保存数组的信息。
        !===============================================================================================
        ! 只保存粗网格时间插值需要的历史，细网格直接读取当前场。
        ! 第三维 0:5 为 rho,u,v,T,Fx/dx,Fy/dx；6:14 为缩放流场非平衡矩；15:19 为缩放温度非平衡矩。
        ! 最后一维 0/1/2 为过去、当前、预测时刻，单网格仅保留第 0 层用于兼容检查点。
        real(kind=8), allocatable :: coarseBuffer(:,:,:,:)
        !===============================================================================================

        !===============================================================================================
        ! 重叠层交换数组的信息。
        !===============================================================================================
        ! 交换时直接读取来源数组。每个线程使用 buffer(0:19) 保存局部插值结果，随即重建目标节点。
        ! 不再分配接口下标表、权重表、全局接收数组或细网格历史副本。
    end module commondata

!===================================================================================================

!===================================================================================================

!=============================================================
!   主程序
!=============================================================

    program main

    use openacc
    use commondata
    implicit none
    real(kind=8) :: timeStart, timeEnd
    real(kind=8) :: timeStart2, timeEnd2
    character(len=24) :: ctime, string
    INTEGER(kind=4) :: time
    integer(kind=4) :: numAccDevices
#ifdef unsteadyFlow
    integer(kind=4) :: unsteadyItcRemaining
#endif
    integer(kind=8) :: wallClockStart, wallClockEnd, wallClockRate
    integer(kind=4) :: lastSnapshotItc, lastPltItc, lastReloadItc    ! 本次最近输出的步数，防止结束时重复写同一步。
    real(kind=8) :: latticeUpdates    ! 统计实际节点更新次数。

    !===============================================================================================
    ! 初始化 OpenACC 设备
    if(loadInitField.EQ.1) then
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,*) " "
        write(00,*) "================ Restart continuation begins ================"
    else
        open(unit=00,file=trim(settingsFile),status='replace')
    endif

    string=ctime(time())    !ctime把 time() 返回的时间戳转换成可读的字符串
    write(00,*) 'Start: ',string    !什么时候开始计算
    write(00,*) "Starting OpenACC >>>>>>"
    call acc_init(acc_device_default)

    numAccDevices=acc_get_num_devices(acc_device_default)
    write(00,*) "Visible OpenACC devices:",numAccDevices
    close(00)
    !===============================================================================================

    !===============================================================================================
    ! Initialization
    ! 新算例初始化网格和场量；续算时在 initial 中恢复检查点。
    call initial()
    call enter_data_2d_openacc()
#ifdef unsteadyFlow
    ! 本次剩余步数；总目标先向上对齐到粗细网格同步时刻。
    unsteadyItcRemaining=max(0,((itc_max+refineRatio-1)/refineRatio)*refineRatio-reloadItcOffset)
#endif
    lastSnapshotItc=-1
    lastPltItc=-1
    lastReloadItc=-1
    !===============================================================================================

    !===============================================================================================
    ! 开始计时
    call CPU_TIME(timeStart)    !当前进程累计消耗的 CPU 时间,包括并行
    ! system_clock 返回墙钟计数器和每秒计数率；
    ! 下面用 counter/rate 把它换算成实际经过的秒数，GPU 计算耗时以墙钟时间为准。
    call system_clock(wallClockStart,wallClockRate)
    timeStart2=dble(wallClockStart)/dble(max(wallClockRate,1_8))

    !===============================================================================================

    !===============================================================================================
    ! 时间推进
    ! itc 在 advance_multiblock 中增加，终止步数向上对齐到粗细同步时刻。
#ifdef steadyFlow
    ! 稳态 itc_max 限制本次运行的步数，速度和温度误差均达标时提前停止。
    do while(((errorU .GT. epsU) .OR. (errorT .GT. epsT)) .AND. &
        (itc .LT. ((itc_max+refineRatio-1)/refineRatio)*refineRatio))
#endif
#ifdef unsteadyFlow
    ! 非稳态 itc_max 为整个算例的总目标，续算只补足剩余步数。
    do while(itc.LT.unsteadyItcRemaining)    !非稳态：续算时只推进到 unsteadyRunDuration 对应的总格子步
#endif

        ! 一次循环推进一个粗时间步，内部完成 refineRatio 个细步和接口交换。
        call advance_multiblock()

#ifdef steadyFlow
        ! 检查间隔由原来的 2000 个细步向上对齐；粗细比不能整除 2000 时也不会漏掉预定检查。
        if(MOD(reloadItcOffset+itc,steadyCheckIntervalItc) .EQ. 0) then
            call check()
            call output_steady_monitor()    ! 输出当前 Nu/Re，不增加正式采样次数。
        endif
#endif

#ifdef steadyFlow
        ! 稳态仍保留按固定间隔记录的 Nu/Re 历史，与上面的收敛监测分开。
#endif
#ifdef unsteadyFlow
        ! 非稳态按 reloadItcOffset+itc 的累计细步采样，接续原来的输出时刻。
#endif
        ! 两种模式都按固定格子步间隔采样；关闭快照文件不停止 Nu/Re 采样。
        if(MOD(reloadItcOffset+itc,outputSnapshotIntervalItc) .EQ. 0) then

            call calNuRe()

            if(outputSnapshotFile.EQ.1) then
                call output_SnapshotFile()
                lastSnapshotItc=itc
            endif

        endif

        ! 输出子程序内部将所需场量从设备同步到主机。
        if((outputPltFile.EQ.1).AND.(MOD(reloadItcOffset+itc,outputPltFileIntervalItc).EQ.0)) then

            call output_Tecplot()
            lastPltItc=itc

        endif

        ! 续算文件按自己的格子步间隔输出，保存粗细同步时刻的完整状态。
        if((outputReloadFile.EQ.1).AND.(MOD(reloadItcOffset+itc,reloadFileIntervalItc).EQ.0)) then

            call output_ReloadFile()
            lastReloadItc=itc

        endif

    enddo
    !===============================================================================================

    !===============================================================================================
    ! 等待设备计算完成，记录实际经过的时间
    !$acc wait(1)

    call CPU_TIME(timeEnd)    !当前进程累计消耗的 CPU 时间,包括并行
    ! 等待 GPU 完成后再读墙钟结束计数器；计时区间不包含下面的最终后处理。
    call system_clock(wallClockEnd,wallClockRate)
    timeEnd2=dble(wallClockEnd)/dble(max(wallClockRate,1_8))

    ! 一个同步周期中粗网格推进一次，四个细区各推进 refineRatio 次；续算只计本次新增步数。
    latticeUpdates=dble(size(f_coarse,1))*dble(size(f_coarse,2))*dble(itc/refineRatio)
    if(refineRatio > 1) then
        latticeUpdates=latticeUpdates+dble(itc)*&
            (dble(size(f_left,1))*dble(size(f_left,2))+dble(size(f_right,1))*dble(size(f_right,2))+&
             dble(size(f_bottom,1))*dble(size(f_bottom,2))+dble(size(f_top,1))*dble(size(f_top,2)))
    endif

    open(unit=00,file=trim(settingsFile),status='unknown',position='append')
    write(00,*) "======================================================================"
    write(00,*) "Time (CPU) = ",timeEnd-timeStart,"s"
    if(timeEnd > timeStart) then
        write(00,*) "MLUPS = ",latticeUpdates/(timeEnd-timeStart)/1.0d6
    else
        write(00,*) "MLUPS = unavailable (CPU timer resolution)"
    endif
    write(00,*) "Time (ACC) = ",timeEnd2-timeStart2,"s"
    if(timeEnd2 > timeStart2) then
        write(00,*) "MLUPS (ACC) = ",latticeUpdates/(timeEnd2-timeStart2)/1.0d6
    else
        write(00,*) "MLUPS (ACC) = unavailable (wall timer resolution)"
    endif
    write(00,*) "Fine steps in this run = ",itc
    write(00,*) "Lattice updates in this run = ",latticeUpdates
    close(00)

    write(*,*) 'Elapsed seconds:',timeEnd2-timeStart2
    !===============================================================================================

    !===============================================================================================
    ! 结束后的统计与最终输出
#ifdef steadyFlow
    ! 最后一步可能不在采样时刻；单独输出最终诊断，不追加 NuRe 历史，也不增加 NuReSampleCount。
    call output_final_NuRe()
#endif
#ifdef unsteadyFlow
    ! 非稳态按指定窗口统计 Nu/Re，并比较前后半窗口；不把最后一步额外算作正式采样。
    call output_unsteady_NuRe_postprocess()
#endif

    ! 保留各输出开关；若最后一步已做过周期输出，就不再重复写同一时刻。
    if((outputPltFile .EQ. 1) .AND. (lastPltItc/=itc)) call output_Tecplot()

    if((outputSnapshotFile .EQ. 1) .AND. (lastSnapshotItc/=itc)) call output_SnapshotFile()

    ! 续算文件最后写，保存最终快照和 Tecplot 输出之后的文件编号。
    ! 若刚新增了最终快照/图形文件，也需更新检查点，避免续算时重复使用文件编号。
    if(outputReloadFile .EQ. 1) then
        if(lastReloadItc/=itc .OR. &
            ((outputPltFile .EQ. 1) .AND. (lastPltItc/=itc)) .OR. &
            ((outputSnapshotFile .EQ. 1) .AND. (lastSnapshotItc/=itc))) call output_ReloadFile()
    endif
    !===============================================================================================

    !===============================================================================================
    ! 先释放设备数组，再释放主机上的五套网格、时间层和接口数组。
    open(unit=00,file=trim(settingsFile),status='unknown',position='append')
    write(00,*) "Deallocate Array......"
    call exit_data_2d_openacc()
    call deallocate_grid_arrays()

#ifdef steadyFlow
    if(errorU .LE. epsU .AND. errorT .LE. epsT) then
        write(00,*) "Steady calculation converged."
    else
        write(00,*) "Maximum synchronized step reached; steady convergence not achieved."
    endif
#endif
#ifdef unsteadyFlow
    write(00,*) "Unsteady target duration completed at a synchronized step."
#endif
    write(00,*) "Final cumulative fine step = ",reloadItcOffset+itc
    write(00,*) "Final t_ff = ",dble(reloadItcOffset+itc)/timeUnit
    string=ctime(time())
    write(00,*) 'End:   ',string    !什么时候算完
    close(00)
    !===============================================================================================

    end program main
!===================================================================================================

!===================================================================================================
! 子程序: initial
! 作用: 检查参数，初始化网格坐标、场变量、分布函数、输出文件和重启信息。
!===================================================================================================
  subroutine initial()

    use commondata
    implicit none

    ! 先确定五个区域的尺寸、坐标偏移和积分范围，再分配数组。
    ! nx/ny 是按最细格距计的全域长度；各区域的本地下标都从 1 开始。
    ! 节点坐标为 xOffset+(i-0.5)*dx；centerIntegrationBox 是积分分界，不是粗网格最外层计算节点。

    integer(kind=4) :: k, i, j
    real(kind=8) :: totalArea, cellArea, xLeft, xRight, yBottom, yTop

    itc=0
    reloadItcOffset=0
    reloadDimensionlessTime=0.0d0
    reloadMetadataLoaded=.false.

    ! 检查输入参数。
    if(refineRatio < 1) error stop 'refineRatio must be a positive integer'

    if(min(nx,ny) < 8) error stop 'At least 8 cells per direction are required'

    if(loadInitField/=0 .AND. loadInitField/=1) error stop 'loadInitField must be 0 or 1'

    if(min(outputSnapshotInterval,reloadFileInterval,outputPltFileInterval)*timeUnit < dble(refineRatio)) &
        error stop 'Output intervals must be at least one synchronized coarse step'    !要求最短的间隔也至少覆盖一个粗细同步周期,输出间隔对应的细步数 ≥ refineRatio

    if(paraA <=-4.0d0 .OR. paraA >=1.0d0) error stop 'Legacy paraA must be in (-4,1)'

    ! 将 t_ff 间隔换算为细格子步数，再向上对齐到粗细同步步。
    outputSnapshotIntervalItc=max(1,ceiling(outputSnapshotInterval*timeUnit/dble(refineRatio)))*refineRatio
    reloadFileIntervalItc=max(1,ceiling(reloadFileInterval*timeUnit/dble(refineRatio)))*refineRatio
    outputPltFileIntervalItc=max(1,ceiling(outputPltFileInterval*timeUnit/dble(refineRatio)))*refineRatio

#ifdef steadyFlow
    steadyCheckIntervalItc=ceiling(2000.0d0/dble(refineRatio))*refineRatio
#endif

    ! 设置流场和温度场的格子权重。
    omega(0)=4.0d0/9.0d0
    omega(1:4)=1.0d0/9.0d0
    omega(5:8)=1.0d0/36.0d0
    omegaT(0)=(1.0d0-paraA)/5.0d0
    omegaT(1:4)=(paraA+4.0d0)/20.0d0

    ! 确定粗格距、粗时间步及需要保存的历史层数。
    dxCoarse=dble(refineRatio)
    timeInterpMaxIndexCoarse=0
    if(refineRatio==1) then
        centerIntegrationBox=[0.0d0,dble(nx),0.0d0,dble(ny)]

        xOffsetCoarse=0.0d0
        yOffsetCoarse=0.0d0
    else
        ! || 表示“或者”。
        ! 左右速度或温度任一周期宏开启时，保留下面的报错语句；都未开启时将其去掉。
        ! 此处位于多块网格分支，目前尚未实现左右周期连接及数据交换，因此初始化时停止。
#if defined(VerticalWallsPeriodicalU) || defined(VerticalWallsPeriodicalT)
        error stop 'Multiblock periodic sides require a periodic block topology; use wall BCs or refineRatio=1'
#endif

        if(interfaceSkin < 2) error stop 'Split flow/thermal advance requires at least two interface layers'

        !检查全域整除，并确定四条分区线
        if(mod(nx,refineRatio)/=0 .OR. mod(ny,refineRatio)/=0) &
            error stop 'nx,ny must be multiples of refineRatio'
        xLeft=dble(fineLayerCellsLeft)-0.5d0    !细节点编号换成坐标
        xRight=dble(nx-fineLayerCellsRight)+0.5d0
        yBottom=dble(fineLayerCellsBottom)-0.5d0
        yTop=dble(ny-fineLayerCellsTop)+0.5d0

        ! 检查细块与粗块的重叠区域
        if(min(xLeft,dble(nx)-xRight,yBottom,dble(ny)-yTop) < dble(interfaceSkin*refineRatio+refineRatio)) &
            error stop 'Refined wall layer must exceed the overlap width by at least one coarse cell'

        if(xRight <=xLeft .OR. yTop <=yBottom) &    !保证中心区域有正宽度和高度
            error stop 'The coarse core must have positive width and height'

        if(mod(nx-fineLayerCellsLeft-fineLayerCellsRight+1,refineRatio)/=0 .OR. &    !中心宽高必须能被 refineRatio 整除
            mod(ny-fineLayerCellsBottom-fineLayerCellsTop+1,refineRatio)/=0) &
            error stop 'Central width and height must be multiples of refineRatio; adjust the four interface indices'

        timeInterpMaxIndexCoarse=2    !表示使用 0、1、2 三个时间层，以及记录中心积分区域
        centerIntegrationBox=[xLeft,xRight,yBottom,yTop]

        ! 将坐标偏移换算为细格距单位
        xOffsetCoarse=xLeft-interfaceSkin*refineRatio-0.5d0*dxCoarse
        yOffsetCoarse=yBottom-interfaceSkin*refineRatio-0.5d0*dxCoarse

        ! 四个细矩形拼成细环；仅向中心延伸 interfaceSkin 个细格距，细矩形拼接处直接迁移。

        ! 左侧细区的坐标偏移。
        xOffsetLeft=0.0d0
        yOffsetLeft=0.0d0

        ! 右侧细区的坐标偏移，包含向中心延伸的 interfaceSkin 层。
        xOffsetRight=dble(nx-(fineLayerCellsRight+interfaceSkin))
        yOffsetRight=0.0d0

        ! 下侧细区从左侧存储区之后开始。
        xOffsetBottom=dble((fineLayerCellsLeft+interfaceSkin))
        yOffsetBottom=0.0d0

        ! 上侧细区从左侧存储区之后开始，并向下延伸 interfaceSkin 层。
        xOffsetTop=dble((fineLayerCellsLeft+interfaceSkin))
        yOffsetTop=dble(ny-(fineLayerCellsTop+interfaceSkin))

        ! 检查细环中间还有空间
        if(min((nx-fineLayerCellsLeft-fineLayerCellsRight-2*interfaceSkin),ny-(fineLayerCellsBottom+interfaceSkin)-(fineLayerCellsTop+interfaceSkin)) < 1) &
            error stop 'Compact fine ring requires a positive interior storage hole'
    endif

    ! 单块模式的四边执行选定的物理域边界条件（壁面或周期）；多块模式的粗网格边界是重叠层边缘。
    coarseIsWall=refineRatio==1

    ! 把细网格的输运参数换算为粗网格参数。
    SnuCoarse=1.0d0/(0.5d0+(1.0d0/Snu-0.5d0)/dxCoarse)
    SqCoarse=1.0d0/(0.5d0+(1.0d0/Sq-0.5d0)/dxCoarse)
    QkCoarse=1.0d0/(0.5d0+(1.0d0/Qk-0.5d0)/dxCoarse)
    QnuCoarse=1.0d0/(0.5d0+(1.0d0/Qnu-0.5d0)/dxCoarse)

    ! 浮力源项的细网格表达式为 Fy = rho*gBeta*(T-Tref)。
    ! 粗时间步相当于 dxCoarse 个细时间步，调用时直接传 dxCoarse*gBeta。

    ! 细网格 dx=1，直接使用基准 Snu、Sq、Qk、Qnu。

    if(refineRatio > 1) then

        ! 左右细区接触三面外墙，上下细区只接触各自的水平外墙。数组顺序为左、右、下、上
        leftIsWall=[.true.,.false.,.true.,.true.]
        rightIsWall=[.false.,.true.,.true.,.true.]
        bottomIsWall=[.false.,.false.,.true.,.false.]
        topIsWall=[.false.,.false.,.false.,.true.]
    endif

    ! 单网格使用全域；多块粗区只遍历分区线之间的节点（含端点）。

    call allocate_grid_arrays()

    ! 初始化各区域，并核对全部实际积分面积之和。
    totalArea=0.0d0

    call initial_grid(size(f_coarse,1),size(f_coarse,2),dxCoarse,xOffsetCoarse,yOffsetCoarse,&
        f_coarse,g_coarse,&
        f_post_coarse,g_post_coarse,rho_coarse,u_coarse,v_coarse,T_coarse,Fx_coarse,Fy_coarse &
#ifdef steadyFlow
        ,up_coarse,vp_coarse,Tp_coarse &
#endif
        )
    do j=1,size(f_coarse,2)
        do i=1,size(f_coarse,1)
            call owned_cell_area(xOffsetCoarse+(i-0.5d0)*dxCoarse,&
                yOffsetCoarse+(j-0.5d0)*dxCoarse,dxCoarse,centerIntegrationBox,cellArea)
            totalArea=totalArea+cellArea
        enddo
    enddo

    if(refineRatio > 1) then

        call initial_grid(size(f_left,1),size(f_left,2),dxFine,xOffsetLeft,yOffsetLeft,&
            f_left,g_left,f_post_left,&
            g_post_left,rho_left,&
            u_left,v_left,T_left,Fx_left,Fy_left &
#ifdef steadyFlow
            ,up_left,vp_left,Tp_left &
#endif
            )
        do j=1,size(f_left,2)
            do i=1,size(f_left,1)
                call owned_cell_area(xOffsetLeft+(i-0.5d0)*dxFine,&
                    yOffsetLeft+(j-0.5d0)*dxFine,dxFine,[xOffsetLeft,xOffsetLeft+size(f_left,1),yOffsetLeft,&
                    yOffsetLeft+size(f_left,2)],cellArea)
                totalArea=totalArea+cellArea
            enddo
        enddo

        call initial_grid(size(f_right,1),size(f_right,2),dxFine,xOffsetRight,yOffsetRight,&
            f_right,g_right,f_post_right,g_post_right,&
            rho_right,u_right,v_right,T_right,Fx_right,Fy_right &
#ifdef steadyFlow
            ,up_right,vp_right,Tp_right &
#endif
            )
        do j=1,size(f_right,2)
            do i=1,size(f_right,1)
                call owned_cell_area(xOffsetRight+(i-0.5d0)*dxFine,&
                    yOffsetRight+(j-0.5d0)*dxFine,dxFine,[xOffsetRight,xOffsetRight+size(f_right,1),yOffsetRight,&
                    yOffsetRight+size(f_right,2)],cellArea)
                totalArea=totalArea+cellArea
            enddo
        enddo

        call initial_grid(size(f_bottom,1),size(f_bottom,2),dxFine,xOffsetBottom,yOffsetBottom,&
            f_bottom,&
            g_bottom,f_post_bottom,&
            g_post_bottom,rho_bottom,u_bottom,v_bottom,T_bottom,Fx_bottom,Fy_bottom &
#ifdef steadyFlow
            ,up_bottom,vp_bottom,Tp_bottom &
#endif
            )
        do j=1,size(f_bottom,2)
            do i=1,size(f_bottom,1)
                call owned_cell_area(xOffsetBottom+(i-0.5d0)*dxFine,&
                    yOffsetBottom+(j-0.5d0)*dxFine,dxFine,[xOffsetBottom,xOffsetBottom+size(f_bottom,1),&
                    yOffsetBottom,yOffsetBottom+size(f_bottom,2)],cellArea)
                totalArea=totalArea+cellArea
            enddo
        enddo

        call initial_grid(size(f_top,1),size(f_top,2),dxFine,xOffsetTop,yOffsetTop,&
            f_top,g_top,f_post_top,g_post_top,rho_top,u_top,v_top,T_top,&
            Fx_top,Fy_top &
#ifdef steadyFlow
            ,up_top,vp_top,Tp_top &
#endif
            )
        do j=1,size(f_top,2)
            do i=1,size(f_top,1)
                call owned_cell_area(xOffsetTop+(i-0.5d0)*dxFine,&
                    yOffsetTop+(j-0.5d0)*dxFine,dxFine,[xOffsetTop,xOffsetTop+size(f_top,1),yOffsetTop,&
                    yOffsetTop+size(f_top,2)],cellArea)
                totalArea=totalArea+cellArea
            enddo
        enddo
    endif

    ! 总面积检查通过后，再恢复续算数据和建立接口。
    if(abs(totalArea-dble(nx)*ny) > 1.0d-8) error stop 'Grid ownership does not tile the physical domain'
    if(loadInitField==1) call read_restart()

    call check_interfaces()

    ! 写出本次算例的参数与各区域几何信息。
    ! 主程序已记录开始时间和设备信息，此处追加参数，续算时保留此前日志。
    open(newunit=k,file=settingsFile,status='unknown',position='append')
    write(k,*) 'Parent: uniform-grid 2DRBOpenacc.F90; D2Q9 MRT / D2Q5 MRT unchanged inside blocks'
    write(k,*) 'Huang (2014): force-aware MRT scaling and four-point Lagrange; original D2Q5 adaptation.'
    write(k,*) 'Yu (2002): aligned two-way coupling and acoustic subcycling; cubic spline is not used.'
    write(k,*) 'Chen (2016): pre-collision buffer states; received f,g participate in the next collision.'
    write(k,*) 'Fine-equivalent nx,ny; refinement ratio:',nx,ny,refineRatio
    write(k,*) 'Rayleigh, Prandtl, Mach:',Rayleigh,Prandtl,Mach
    write(k,*) 'Fine tauf, Snu, Sq, Qk, Qnu, paraA:',tauf,Snu,Sq,Qk,Qnu,paraA
    write(k,*) 'Fine viscosity, diffusivity, gBeta, timeUnit:',viscosity,diffusivity,gBeta,timeUnit
    write(k,*) 'Restart offset itc =',reloadItcOffset
    write(k,*) 'Reload dimensionless time =',reloadDimensionlessTime
    write(k,*) 'Reload metadata loaded =',reloadMetadataLoaded
    write(k,*) 'outputSnapshotIntervalItc =',outputSnapshotIntervalItc
    write(k,*) 'reloadFileIntervalItc =',reloadFileIntervalItc
    write(k,*) 'outputPltFileIntervalItc =',outputPltFileIntervalItc
#ifdef steadyFlow
    write(k,*) 'steadyCheckIntervalItc =',steadyCheckIntervalItc
#endif
    write(k,*) 'Fine-node indices left/right/bottom/top:',&
        fineLayerCellsLeft,fineLayerCellsRight,fineLayerCellsBottom,fineLayerCellsTop
    write(k,*) 'Coarse extension (coarse spacings); fine extension (fine spacings); interface layers:',&
        interfaceSkin,interfaceSkin,interfaceSkin
    write(k,*) 'Owned physical area:',totalArea
    write(k,*) 'Geometry: central coarse rectangle and four compact fine rectangles.'
    write(k,*) 'Fine seams exchange post-collision populations directly, including diagonal links.'
    write(k,*) 'Coordinates: x=xOffset+(i-0.5)*dx, y=yOffset+(j-0.5)*dx; local i,j start at 1.'
    write(k,*) 'Interfaces preserve the connected-ring time interpolation and MRT scaling; dt=dx.'
    write(k,*) 'Restart stores coarse,left,right,bottom,top arrays without a text format header.'
    write(k,*) 'coarse: nx,ny,dx,xOffset,yOffset:',size(f_coarse,1),size(f_coarse,2),dxCoarse,xOffsetCoarse,yOffsetCoarse
    if(refineRatio > 1) then
        write(k,*) 'left: nx,ny,dx,xOffset,yOffset:',size(f_left,1),size(f_left,2),dxFine,xOffsetLeft,yOffsetLeft
        write(k,*) 'right: nx,ny,dx,xOffset,yOffset:',size(f_right,1),size(f_right,2),dxFine,xOffsetRight,yOffsetRight
        write(k,*) 'bottom: nx,ny,dx,xOffset,yOffset:',size(f_bottom,1),size(f_bottom,2),dxFine,xOffsetBottom,&
            yOffsetBottom
        write(k,*) 'top: nx,ny,dx,xOffset,yOffset:',size(f_top,1),size(f_top,2),dxFine,xOffsetTop,yOffsetTop
    endif

    close(k)
    if(loadInitField==0) then
        open(newunit=k,file=NuReHistoryFile,status='replace')
        write(k,'(a)') '# t_ff NuVolAvg ReVolRMS Nu_hot Nu_cold Nu_middle mass meanT Tmin Tmax rhoMin rhoMax'
        close(k)
#ifdef steadyFlow
        ! 新算例清空旧收敛记录；续算保留原文件，并从恢复的累计步数继续追加。
        open(newunit=k,file='Convergence_2DOpenaccMultiblock.dat',status='replace')
        close(k)
        open(newunit=k,file='SteadyMonitor_2DOpenaccMultiblock.dat',status='replace')
        write(k,'(a)') '# itc t_ff errorU errorT NuVolAvg ReVolRMS Nu_hot Nu_cold Nu_middle'
        close(k)
#endif
    else

        call check_history()
    endif
  end subroutine initial
!===================================================================================================

!===================================================================================================
! 子程序: allocate_grid_arrays
! 作用: 按中心、左、右、下、上五个区域的实际尺寸分配数组。
!===================================================================================================
  subroutine allocate_grid_arrays()
    use commondata
    implicit none

    ! 分别为 coarse、left、right、bottom、top 分配数组，尺寸直接由分区位置和重叠厚度计算。
    ! 先分配 f，其余场量用 size(f,1/2) 取得同一区域的尺寸，不另存尺寸变量。
    ! 场量空间范围是 1:ni、1:nj；只有 f_post/g_post 为迁移额外保留 0 和 n+1。
    ! 粗网格缓冲数组保存 0:2 三层；up/vp/Tp 保存上次稳态检查的场量。

    ! 分配中心粗网格的计算数组。
    if(refineRatio==1) then
        allocate(f_coarse(nx,ny,0:8))
    else
        allocate(f_coarse((nx-fineLayerCellsLeft-fineLayerCellsRight+1)/refineRatio+1+2*interfaceSkin,&
            (ny-fineLayerCellsBottom-fineLayerCellsTop+1)/refineRatio+1+2*interfaceSkin,0:8))
    endif
    allocate(g_coarse(size(f_coarse,1),size(f_coarse,2),0:4))
    allocate(f_post_coarse(0:size(f_coarse,1)+1,0:size(f_coarse,2)+1,0:8))
    allocate(g_post_coarse(0:size(f_coarse,1)+1,0:size(f_coarse,2)+1,0:4))
    allocate(rho_coarse(size(f_coarse,1),size(f_coarse,2)))
    allocate(u_coarse(size(f_coarse,1),size(f_coarse,2)))
    allocate(v_coarse(size(f_coarse,1),size(f_coarse,2)))
    allocate(T_coarse(size(f_coarse,1),size(f_coarse,2)))
    allocate(Fx_coarse(size(f_coarse,1),size(f_coarse,2)))
    allocate(Fy_coarse(size(f_coarse,1),size(f_coarse,2)))

    ! 粗时间层合并为一个数组；分量次序保持原检查点中的存储顺序。
    allocate(coarseBuffer(size(f_coarse,1),size(f_coarse,2),0:19,0:timeInterpMaxIndexCoarse))
    coarseBuffer=0.0d0

#ifdef steadyFlow
    allocate(up_coarse(size(f_coarse,1),size(f_coarse,2)))
    allocate(vp_coarse(size(f_coarse,1),size(f_coarse,2)))
    allocate(Tp_coarse(size(f_coarse,1),size(f_coarse,2)))
#endif
    if(refineRatio > 1) then

        ! 分配左侧细网格的计算数组。
        allocate(f_left((fineLayerCellsLeft+interfaceSkin),ny,0:8))
        allocate(g_left(size(f_left,1),size(f_left,2),0:4))
        allocate(f_post_left(0:size(f_left,1)+1,0:size(f_left,2)+1,0:8))
        allocate(g_post_left(0:size(f_left,1)+1,0:size(f_left,2)+1,0:4))
        allocate(rho_left(size(f_left,1),size(f_left,2)))
        allocate(u_left(size(f_left,1),size(f_left,2)))
        allocate(v_left(size(f_left,1),size(f_left,2)))
        allocate(T_left(size(f_left,1),size(f_left,2)))
        allocate(Fx_left(size(f_left,1),size(f_left,2)))
        allocate(Fy_left(size(f_left,1),size(f_left,2)))

#ifdef steadyFlow
        allocate(up_left(size(f_left,1),size(f_left,2)))
        allocate(vp_left(size(f_left,1),size(f_left,2)))
        allocate(Tp_left(size(f_left,1),size(f_left,2)))
#endif

        ! 分配右侧细网格的计算数组。
        allocate(f_right((fineLayerCellsRight+interfaceSkin),ny,0:8))
        allocate(g_right(size(f_right,1),size(f_right,2),0:4))
        allocate(f_post_right(0:size(f_right,1)+1,0:size(f_right,2)+1,0:8))
        allocate(g_post_right(0:size(f_right,1)+1,0:size(f_right,2)+1,0:4))
        allocate(rho_right(size(f_right,1),size(f_right,2)))
        allocate(u_right(size(f_right,1),size(f_right,2)))
        allocate(v_right(size(f_right,1),size(f_right,2)))
        allocate(T_right(size(f_right,1),size(f_right,2)))
        allocate(Fx_right(size(f_right,1),size(f_right,2)))
        allocate(Fy_right(size(f_right,1),size(f_right,2)))

#ifdef steadyFlow
        allocate(up_right(size(f_right,1),size(f_right,2)))
        allocate(vp_right(size(f_right,1),size(f_right,2)))
        allocate(Tp_right(size(f_right,1),size(f_right,2)))
#endif

        ! 分配下侧细网格的计算数组。
        allocate(f_bottom((nx-fineLayerCellsLeft-fineLayerCellsRight-2*interfaceSkin),(fineLayerCellsBottom+interfaceSkin),0:8))
        allocate(g_bottom(size(f_bottom,1),size(f_bottom,2),0:4))
        allocate(f_post_bottom(0:size(f_bottom,1)+1,0:size(f_bottom,2)+1,0:8))
        allocate(g_post_bottom(0:size(f_bottom,1)+1,0:size(f_bottom,2)+1,0:4))
        allocate(rho_bottom(size(f_bottom,1),size(f_bottom,2)))
        allocate(u_bottom(size(f_bottom,1),size(f_bottom,2)))
        allocate(v_bottom(size(f_bottom,1),size(f_bottom,2)))
        allocate(T_bottom(size(f_bottom,1),size(f_bottom,2)))
        allocate(Fx_bottom(size(f_bottom,1),size(f_bottom,2)))
        allocate(Fy_bottom(size(f_bottom,1),size(f_bottom,2)))

#ifdef steadyFlow
        allocate(up_bottom(size(f_bottom,1),size(f_bottom,2)))
        allocate(vp_bottom(size(f_bottom,1),size(f_bottom,2)))
        allocate(Tp_bottom(size(f_bottom,1),size(f_bottom,2)))
#endif

        ! 分配上侧细网格的计算数组。
        allocate(f_top((nx-fineLayerCellsLeft-fineLayerCellsRight-2*interfaceSkin),(fineLayerCellsTop+interfaceSkin),0:8))
        allocate(g_top(size(f_top,1),size(f_top,2),0:4))
        allocate(f_post_top(0:size(f_top,1)+1,0:size(f_top,2)+1,0:8))
        allocate(g_post_top(0:size(f_top,1)+1,0:size(f_top,2)+1,0:4))
        allocate(rho_top(size(f_top,1),size(f_top,2)))
        allocate(u_top(size(f_top,1),size(f_top,2)))
        allocate(v_top(size(f_top,1),size(f_top,2)))
        allocate(T_top(size(f_top,1),size(f_top,2)))
        allocate(Fx_top(size(f_top,1),size(f_top,2)))
        allocate(Fy_top(size(f_top,1),size(f_top,2)))

#ifdef steadyFlow
        allocate(up_top(size(f_top,1),size(f_top,2)))
        allocate(vp_top(size(f_top,1),size(f_top,2)))
        allocate(Tp_top(size(f_top,1),size(f_top,2)))
#endif
    endif
  end subroutine allocate_grid_arrays
!===================================================================================================

!===================================================================================================
! 子程序: initial_grid
! 作用: 初始化传入区域的宏观量、分布函数和交换历史。
!===================================================================================================
  subroutine initial_grid(ni,nj,dx,xOffset,yOffset,&
      f,g,&
      f_post,g_post,rho,u,v,T,Fx,Fy &
#ifdef steadyFlow
      ,up,vp,Tp &
#endif
      )

    use commondata,only: nx,ny,Thot,Tcold,pi,lengthUnit,omega,omegaT
    implicit none

    ! ni/nj 是传入数组的本地节点数。
    ! dx、xOffset、yOffset 把本地下标换成物理位置；交换层也需初始化。

    integer(kind=4) :: ni, nj
    real(kind=8) :: dx, xOffset, yOffset

    real(kind=8) :: f(ni,nj,0:8)
    real(kind=8) :: g(ni,nj,0:4)
    real(kind=8) :: f_post(0:ni+1,0:nj+1,0:8)
    real(kind=8) :: g_post(0:ni+1,0:nj+1,0:4)
    real(kind=8) :: rho(ni,nj)
    real(kind=8) :: u(ni,nj)
    real(kind=8) :: v(ni,nj)
    real(kind=8) :: T(ni,nj)
    real(kind=8) :: Fx(ni,nj)
    real(kind=8) :: Fy(ni,nj)

#ifdef steadyFlow
    real(kind=8) :: up(ni,nj)
    real(kind=8) :: vp(ni,nj)
    real(kind=8) :: Tp(ni,nj)
#endif
    integer(kind=4) :: i, j, a
    real(kind=8) :: x, y

    ! 初始化宏观量和碰撞后数组。
    u=0.0d0
    v=0.0d0
    rho=1.0d0
    Fx=0.0d0
    Fy=0.0d0
    f_post=0.0d0
    g_post=0.0d0

    ! 根据节点物理坐标构造初始温度，并设置平衡分布函数。
    do j=1,nj
        y=(yOffset+(dble(j)-0.5d0)*dx)/lengthUnit
        do i=1,ni
            x=(xOffset+(dble(i)-0.5d0)*dx)/lengthUnit
#ifdef SideHeatedCell
            T(i,j)=Thot+x/(dble(nx)/lengthUnit)*(Tcold-Thot)
#else
            T(i,j)=Thot+y/(dble(ny)/lengthUnit)*(Tcold-Thot)
            T(i,j)=T(i,j)+1.0d-3*(Thot-Tcold)*sin(2.0d0*pi*x/(dble(nx)/lengthUnit))*&
                sin(pi*y/(dble(ny)/lengthUnit))
#endif
            do a=0,8
                f(i,j,a)=omega(a)
            enddo

            do a=0,4
                g(i,j,a)=omegaT(a)*T(i,j)
            enddo
        enddo
    enddo
#ifdef steadyFlow

    ! 稳态检查的初始参考场。
    up=u
    vp=v
    Tp=T
#endif
  end subroutine initial_grid
!===================================================================================================

!===================================================================================================
! 子程序: grid_device_data
! 作用: 统一建立或释放五套网格场量和交换历史的 OpenACC 设备存储。
!===================================================================================================
  subroutine grid_device_data(entering)
    use commondata
    implicit none
    logical, intent(in) :: entering

    ! 场量和粗时间层常驻设备；全部交换在同一设备内直接读取来源数组。
    if(entering) then
        !$acc enter data copyin(f_coarse, g_coarse, f_post_coarse, g_post_coarse, rho_coarse, u_coarse, &
        !$acc& v_coarse, T_coarse, Fx_coarse, Fy_coarse, coarseBuffer)
        if(refineRatio > 1) then
            !$acc enter data copyin(f_left, g_left, f_post_left, g_post_left, rho_left, u_left, &
            !$acc& v_left, T_left, Fx_left, Fy_left)
            !$acc enter data copyin(f_right, g_right, f_post_right, g_post_right, rho_right, u_right, &
            !$acc& v_right, T_right, Fx_right, Fy_right)
            !$acc enter data copyin(f_bottom, g_bottom, f_post_bottom, g_post_bottom, rho_bottom, u_bottom, &
            !$acc& v_bottom, T_bottom, Fx_bottom, Fy_bottom)
            !$acc enter data copyin(f_top, g_top, f_post_top, g_post_top, rho_top, u_top, &
            !$acc& v_top, T_top, Fx_top, Fy_top)
        endif
    else
        !$acc exit data delete(f_coarse, g_coarse, f_post_coarse, g_post_coarse, rho_coarse, u_coarse, &
        !$acc& v_coarse, T_coarse, Fx_coarse, Fy_coarse, coarseBuffer)
        if(refineRatio > 1) then
            !$acc exit data delete(f_left, g_left, f_post_left, g_post_left, rho_left, u_left, &
            !$acc& v_left, T_left, Fx_left, Fy_left)
            !$acc exit data delete(f_right, g_right, f_post_right, g_post_right, rho_right, u_right, &
            !$acc& v_right, T_right, Fx_right, Fy_right)
            !$acc exit data delete(f_bottom, g_bottom, f_post_bottom, g_post_bottom, rho_bottom, u_bottom, &
            !$acc& v_bottom, T_bottom, Fx_bottom, Fy_bottom)
            !$acc exit data delete(f_top, g_top, f_post_top, g_post_top, rho_top, u_top, &
            !$acc& v_top, T_top, Fx_top, Fy_top)
        endif
    endif
  end subroutine grid_device_data
!===================================================================================================

!===================================================================================================
! 子程序: enter_data_2d_openacc
! 作用: 建立设备数据区，并补齐首次碰撞前的粗细接口状态。
!===================================================================================================
  subroutine enter_data_2d_openacc()
    use commondata
    implicit none

    !$acc enter data copyin(ex, ey, omega, omegaT)
    call grid_device_data(.true.)
    ! 新算例先保存起点，再双向补齐接口；续算直接使用读入的粗时间层。
    if(loadInitField==0) then
        call save_coarse_buffer(min(1,timeInterpMaxIndexCoarse))
        if(refineRatio > 1) then
            call fine_to_coarse()
            call coarse_to_fine([0.0d0,1.0d0,0.0d0])
            call save_coarse_buffer(1)
        endif
    endif
  end subroutine enter_data_2d_openacc
!===================================================================================================

!===================================================================================================
! 子程序: update_host_grid
! 作用: 将指定区域的宏观量或完整状态从设备同步到主机。
!===================================================================================================
  subroutine update_host_grid(ni,nj,full,f,g,rho,u,v,T,Fx,Fy)

    implicit none

    ! full 为假只同步 rho/u/v/T，供统计和快照使用；为真还同步 f/g 和力；粗历史在 update_host_all 中同步。
    ! 本例使用 async(1) 排队，同步完成由 update_host_all 中的 wait 保证。

    integer(kind=4) :: ni, nj

    real(kind=8) :: f(ni,nj,0:8)
    real(kind=8) :: g(ni,nj,0:4)
    real(kind=8) :: rho(ni,nj)
    real(kind=8) :: u(ni,nj)
    real(kind=8) :: v(ni,nj)
    real(kind=8) :: T(ni,nj)
    real(kind=8) :: Fx(ni,nj)
    real(kind=8) :: Fy(ni,nj)

    logical, intent(in) :: full

    !$acc update self(u, v, T, rho) async(1)
    if(full) then
        !$acc update self(f, g, Fx, Fy) async(1)
    endif
  end subroutine update_host_grid
!===================================================================================================

!===================================================================================================
! 子程序: update_host_all
! 作用: 依次同步五套具名数组，并等待同步完成。
!===================================================================================================
  subroutine update_host_all(full)
    use commondata
    implicit none
    logical, intent(in) :: full

    call update_host_grid(size(f_coarse,1),size(f_coarse,2),full,f_coarse,g_coarse,rho_coarse,&
        u_coarse,v_coarse,T_coarse,Fx_coarse,Fy_coarse)
    if(refineRatio > 1) then
        call update_host_grid(size(f_left,1),size(f_left,2),full,f_left,g_left,rho_left,&
            u_left,v_left,T_left,Fx_left,Fy_left)
        call update_host_grid(size(f_right,1),size(f_right,2),full,f_right,g_right,rho_right,&
            u_right,v_right,T_right,Fx_right,Fy_right)
        call update_host_grid(size(f_bottom,1),size(f_bottom,2),full,f_bottom,g_bottom,rho_bottom,&
            u_bottom,v_bottom,T_bottom,Fx_bottom,Fy_bottom)
        call update_host_grid(size(f_top,1),size(f_top,2),full,f_top,g_top,rho_top,&
            u_top,v_top,T_top,Fx_top,Fy_top)
    endif
    if(full) then
        !$acc update self(coarseBuffer) async(1)
    endif
    !$acc wait(1)
  end subroutine update_host_all
!===================================================================================================

!===================================================================================================
! 子程序: exit_data_2d_openacc
! 作用: 等待设备任务结束，释放场量、粗时间层及格子常数的设备存储。
!===================================================================================================
  subroutine exit_data_2d_openacc()
    use commondata
    implicit none

    !$acc wait(1)

    call grid_device_data(.false.)
    !$acc exit data delete(ex, ey, omega, omegaT)
  end subroutine exit_data_2d_openacc
!===================================================================================================

!===================================================================================================
! 子程序: deallocate_grid_arrays
! 作用: 释放主机端五套网格、交换时间层和接口数组；须在设备数组释放之后调用。
!===================================================================================================
  subroutine deallocate_grid_arrays()
    use commondata
    implicit none

    deallocate(coarseBuffer)
    deallocate(f_coarse,g_coarse,f_post_coarse,g_post_coarse,rho_coarse,u_coarse,v_coarse,T_coarse,Fx_coarse,Fy_coarse)
#ifdef steadyFlow
    deallocate(up_coarse,vp_coarse,Tp_coarse)
#endif
    if(refineRatio > 1) then
        deallocate(f_left,g_left,f_post_left,g_post_left,rho_left,u_left,v_left,T_left,Fx_left,Fy_left)
#ifdef steadyFlow
        deallocate(up_left,vp_left,Tp_left)
#endif
        deallocate(f_right,g_right,f_post_right,g_post_right,rho_right,u_right,v_right,T_right,Fx_right,Fy_right)
#ifdef steadyFlow
        deallocate(up_right,vp_right,Tp_right)
#endif
        deallocate(f_bottom,g_bottom,f_post_bottom,g_post_bottom,rho_bottom,u_bottom,v_bottom,T_bottom,Fx_bottom,Fy_bottom)
#ifdef steadyFlow
        deallocate(up_bottom,vp_bottom,Tp_bottom)
#endif
        deallocate(f_top,g_top,f_post_top,g_post_top,rho_top,u_top,v_top,T_top,Fx_top,Fy_top)
#ifdef steadyFlow
        deallocate(up_top,vp_top,Tp_top)
#endif
    endif
  end subroutine deallocate_grid_arrays
!===================================================================================================

!===================================================================================================
! 子程序: advance_grid
! 作用: 按原均匀网格算法次序，推进传入数组的流场和温度场。
!===================================================================================================
  subroutine advance_grid(ni,nj,dx,Snu,Sq,Qk,Qnu,gBeta,wall,f,g,f_post,g_post,rho,u,v,T,Fx,Fy)

#ifdef SideHeatedHa
    use commondata,only: B2sigemarho
#endif
    implicit none

    ! 先推进流场，再用更新后的速度推进温度场；数组和松弛率由调用处传入。
    ! wall(1:4) 依次表示左、右、下、上是否为真实墙面；粗细人工接口不能当墙面反弹。

    integer(kind=4) :: ni, nj
    real(kind=8) :: dx, Snu, Sq, Qk, Qnu, gBeta
    logical :: wall(4)

    real(kind=8) :: f(ni,nj,0:8)
    real(kind=8) :: g(ni,nj,0:4)
    real(kind=8) :: f_post(0:ni+1,0:nj+1,0:8)
    real(kind=8) :: g_post(0:ni+1,0:nj+1,0:4)
    real(kind=8) :: rho(ni,nj)
    real(kind=8) :: u(ni,nj)
    real(kind=8) :: v(ni,nj)
    real(kind=8) :: T(ni,nj)
    real(kind=8) :: Fx(ni,nj)
    real(kind=8) :: Fy(ni,nj)

    ! 粗网格的时间步长已计入松弛率和力增量。
    ! 推进前，人工边界已补成本时刻的碰撞前状态，接收节点也参加碰撞。
    ! 推进后，最外两层等待下次接口重建，不能作为反向插值的数据来源。
    call collision(ni,nj,f,f_post,rho,u,v,Fx,Fy,T,Snu,Sq,&
        gBeta &
#ifdef SideHeatedHa
        ,dx*B2sigemarho &
#endif
        )

    call streaming(ni,nj,f,f_post)

    call bounceback(ni,nj,f,f_post,wall(1),wall(2),&
        wall(3),wall(4))

    call macro(ni,nj,f,rho,u,v,Fx,Fy)

    call collisionT(ni,nj,g,g_post,u,v,T,Qk,Qnu)

    call streamingT(ni,nj,g,g_post)

    call bouncebackT(ni,nj,g,g_post,wall(1),wall(2),&
        wall(3),wall(4))

    call macroT(ni,nj,g,T)
  end subroutine advance_grid
!===================================================================================================

!===================================================================================================
! 子程序: advance_coarse
! 作用: 使用中心粗网格数组及粗网格参数，推进一个粗时间步。
!===================================================================================================
  subroutine advance_coarse()
    use commondata
    implicit none

    call advance_grid(size(f_coarse,1),size(f_coarse,2),dxCoarse,SnuCoarse,SqCoarse,QkCoarse,QnuCoarse,&
        (dxCoarse*gBeta),coarseIsWall,f_coarse,g_coarse,f_post_coarse,g_post_coarse,rho_coarse,&
        u_coarse,v_coarse,T_coarse,Fx_coarse,Fy_coarse)
  end subroutine advance_coarse
!===================================================================================================

!===================================================================================================
! 子程序: advance_fine
! 作用: 分阶段推进四个细区域，并在迁移前交换同级相邻节点数据。
!===================================================================================================
  subroutine advance_fine()
    use commondata
    implicit none

    ! 四个细区先全部碰撞，再交换相邻区的碰撞后分布，最后迁移。
    ! 流场和温度场各有一次同级交换，保证温度碰撞使用已经更新的速度。
    ! 细区拼接线只改变存储位置，不引入物理边界、空间插值或时间插值。

    ! 四个细区先全部碰撞，再交换迁移外圈；避免读取另一细区尚未更新的数据。
    call collision(size(f_left,1),size(f_left,2),f_left,f_post_left,rho_left,u_left,v_left,Fx_left,Fy_left,&
        T_left,Snu,Sq,gBeta &
#ifdef SideHeatedHa
        ,B2sigemarho &
#endif
        )

    call collision(size(f_right,1),size(f_right,2),f_right,f_post_right,rho_right,u_right,v_right,Fx_right,&
        Fy_right,T_right,Snu,Sq,gBeta &
#ifdef SideHeatedHa
        ,B2sigemarho &
#endif
        )

    call collision(size(f_bottom,1),size(f_bottom,2),f_bottom,f_post_bottom,rho_bottom,u_bottom,v_bottom,&
        Fx_bottom,Fy_bottom,T_bottom,Snu,Sq,gBeta &
#ifdef SideHeatedHa
        ,B2sigemarho &
#endif
        )

    call collision(size(f_top,1),size(f_top,2),f_top,f_post_top,rho_top,u_top,v_top,Fx_top,Fy_top,T_top,&
        Snu,Sq,gBeta &
#ifdef SideHeatedHa
        ,B2sigemarho &
#endif
        )

    ! 四区流场碰撞完成后，复制相邻细区的 f_post，再开始迁移。
    call exchange_fine_flow_halos()

    call streaming(size(f_left,1),size(f_left,2),f_left,f_post_left)

    call bounceback(size(f_left,1),size(f_left,2),f_left,f_post_left,leftIsWall(1),leftIsWall(2),leftIsWall(3),leftIsWall(4))

    call macro(size(f_left,1),size(f_left,2),f_left,rho_left,u_left,v_left,Fx_left,Fy_left)

    call streaming(size(f_right,1),size(f_right,2),f_right,f_post_right)

    call bounceback(size(f_right,1),size(f_right,2),f_right,f_post_right,rightIsWall(1),rightIsWall(2),rightIsWall(3),&
        rightIsWall(4))

    call macro(size(f_right,1),size(f_right,2),f_right,rho_right,u_right,v_right,Fx_right,Fy_right)

    call streaming(size(f_bottom,1),size(f_bottom,2),f_bottom,f_post_bottom)

    call bounceback(size(f_bottom,1),size(f_bottom,2),f_bottom,f_post_bottom,bottomIsWall(1),bottomIsWall(2),&
        bottomIsWall(3),bottomIsWall(4))

    call macro(size(f_bottom,1),size(f_bottom,2),f_bottom,rho_bottom,u_bottom,v_bottom,Fx_bottom,Fy_bottom)

    call streaming(size(f_top,1),size(f_top,2),f_top,f_post_top)

    call bounceback(size(f_top,1),size(f_top,2),f_top,f_post_top,topIsWall(1),topIsWall(2),topIsWall(3),topIsWall(4))

    call macro(size(f_top,1),size(f_top,2),f_top,rho_top,u_top,v_top,Fx_top,Fy_top)

    ! 四区速度已经更新，开始温度碰撞；此处使用更新后的 u/v。
    call collisionT(size(f_left,1),size(f_left,2),g_left,g_post_left,u_left,v_left,T_left,Qk,Qnu)

    call collisionT(size(f_right,1),size(f_right,2),g_right,g_post_right,u_right,v_right,T_right,Qk,Qnu)

    call collisionT(size(f_bottom,1),size(f_bottom,2),g_bottom,g_post_bottom,u_bottom,v_bottom,T_bottom,Qk,Qnu)

    call collisionT(size(f_top,1),size(f_top,2),g_top,g_post_top,u_top,v_top,T_top,Qk,Qnu)

    ! 四区温度碰撞完成后，复制相邻细区的 g_post，再开始温度迁移。
    call exchange_fine_thermal_halos()

    call streamingT(size(f_left,1),size(f_left,2),g_left,g_post_left)

    call bouncebackT(size(f_left,1),size(f_left,2),g_left,g_post_left,leftIsWall(1),leftIsWall(2),leftIsWall(3),leftIsWall(4))

    call macroT(size(f_left,1),size(f_left,2),g_left,T_left)

    call streamingT(size(f_right,1),size(f_right,2),g_right,g_post_right)

    call bouncebackT(size(f_right,1),size(f_right,2),g_right,g_post_right,rightIsWall(1),rightIsWall(2),rightIsWall(3),&
        rightIsWall(4))

    call macroT(size(f_right,1),size(f_right,2),g_right,T_right)

    call streamingT(size(f_bottom,1),size(f_bottom,2),g_bottom,g_post_bottom)

    call bouncebackT(size(f_bottom,1),size(f_bottom,2),g_bottom,g_post_bottom,bottomIsWall(1),bottomIsWall(2),&
        bottomIsWall(3),bottomIsWall(4))

    call macroT(size(f_bottom,1),size(f_bottom,2),g_bottom,T_bottom)

    call streamingT(size(f_top,1),size(f_top,2),g_top,g_post_top)

    call bouncebackT(size(f_top,1),size(f_top,2),g_top,g_post_top,topIsWall(1),topIsWall(2),topIsWall(3),topIsWall(4))

    call macroT(size(f_top,1),size(f_top,2),g_top,T_top)
  end subroutine advance_fine
!===================================================================================================

!===================================================================================================
! 子程序: exchange_fine_flow_halos
! 作用: 在左右细区与上下细区之间双向复制流场碰撞后的迁移外圈。
!===================================================================================================
  subroutine exchange_fine_flow_halos()
    use commondata
    implicit none

    call copy_fine_halo(9,size(f_left,1),size(f_left,2),nint(xOffsetLeft),nint(yOffsetLeft),f_post_left,size(f_bottom,1),&
        size(f_bottom,2),nint(xOffsetBottom),nint(yOffsetBottom),f_post_bottom)

    call copy_fine_halo(9,size(f_bottom,1),size(f_bottom,2),nint(xOffsetBottom),nint(yOffsetBottom),f_post_bottom,&
        size(f_left,1),size(f_left,2),nint(xOffsetLeft),nint(yOffsetLeft),f_post_left)

    call copy_fine_halo(9,size(f_left,1),size(f_left,2),nint(xOffsetLeft),nint(yOffsetLeft),f_post_left,size(f_top,1),&
        size(f_top,2),nint(xOffsetTop),nint(yOffsetTop),f_post_top)

    call copy_fine_halo(9,size(f_top,1),size(f_top,2),nint(xOffsetTop),nint(yOffsetTop),f_post_top,size(f_left,1),&
        size(f_left,2),nint(xOffsetLeft),nint(yOffsetLeft),f_post_left)

    call copy_fine_halo(9,size(f_right,1),size(f_right,2),nint(xOffsetRight),nint(yOffsetRight),f_post_right,&
        size(f_bottom,1),size(f_bottom,2),nint(xOffsetBottom),nint(yOffsetBottom),f_post_bottom)

    call copy_fine_halo(9,size(f_bottom,1),size(f_bottom,2),nint(xOffsetBottom),nint(yOffsetBottom),f_post_bottom,&
        size(f_right,1),size(f_right,2),nint(xOffsetRight),nint(yOffsetRight),f_post_right)

    call copy_fine_halo(9,size(f_right,1),size(f_right,2),nint(xOffsetRight),nint(yOffsetRight),f_post_right,&
        size(f_top,1),size(f_top,2),nint(xOffsetTop),nint(yOffsetTop),f_post_top)

    call copy_fine_halo(9,size(f_top,1),size(f_top,2),nint(xOffsetTop),nint(yOffsetTop),f_post_top,size(f_right,1),&
        size(f_right,2),nint(xOffsetRight),nint(yOffsetRight),f_post_right)
  end subroutine exchange_fine_flow_halos
!===================================================================================================

!===================================================================================================
! 子程序: exchange_fine_thermal_halos
! 作用: 在左右细区与上下细区之间双向复制温度碰撞后的迁移外圈。
!===================================================================================================
  subroutine exchange_fine_thermal_halos()
    use commondata
    implicit none

    call copy_fine_halo(5,size(f_left,1),size(f_left,2),nint(xOffsetLeft),nint(yOffsetLeft),g_post_left,size(f_bottom,1),&
        size(f_bottom,2),nint(xOffsetBottom),nint(yOffsetBottom),g_post_bottom)

    call copy_fine_halo(5,size(f_bottom,1),size(f_bottom,2),nint(xOffsetBottom),nint(yOffsetBottom),g_post_bottom,&
        size(f_left,1),size(f_left,2),nint(xOffsetLeft),nint(yOffsetLeft),g_post_left)

    call copy_fine_halo(5,size(f_left,1),size(f_left,2),nint(xOffsetLeft),nint(yOffsetLeft),g_post_left,size(f_top,1),&
        size(f_top,2),nint(xOffsetTop),nint(yOffsetTop),g_post_top)

    call copy_fine_halo(5,size(f_top,1),size(f_top,2),nint(xOffsetTop),nint(yOffsetTop),g_post_top,size(f_left,1),&
        size(f_left,2),nint(xOffsetLeft),nint(yOffsetLeft),g_post_left)

    call copy_fine_halo(5,size(f_right,1),size(f_right,2),nint(xOffsetRight),nint(yOffsetRight),g_post_right,&
        size(f_bottom,1),size(f_bottom,2),nint(xOffsetBottom),nint(yOffsetBottom),g_post_bottom)

    call copy_fine_halo(5,size(f_bottom,1),size(f_bottom,2),nint(xOffsetBottom),nint(yOffsetBottom),g_post_bottom,&
        size(f_right,1),size(f_right,2),nint(xOffsetRight),nint(yOffsetRight),g_post_right)

    call copy_fine_halo(5,size(f_right,1),size(f_right,2),nint(xOffsetRight),nint(yOffsetRight),g_post_right,&
        size(f_top,1),size(f_top,2),nint(xOffsetTop),nint(yOffsetTop),g_post_top)

    call copy_fine_halo(5,size(f_top,1),size(f_top,2),nint(xOffsetTop),nint(yOffsetTop),g_post_top,size(f_right,1),&
        size(f_right,2),nint(xOffsetRight),nint(yOffsetRight),g_post_right)
  end subroutine exchange_fine_thermal_halos
!===================================================================================================

!===================================================================================================
! 子程序: copy_fine_halo
! 作用: 根据坐标偏移，复制相邻细区域内部与接收区域迁移外圈的交集。
!===================================================================================================
  subroutine copy_fine_halo(nq,ni,nj,targetIndexX,targetIndexY,receiver,mi,mj,sourceIndexX,sourceIndexY,donor)
    implicit none

    ! nq 为方向数：流场 9，温度场 5；ni/nj 属于接收区，mi/mj 属于来源区。
    ! targetIndexX/targetIndexY、sourceIndexX/sourceIndexY 是以细格距计的两区整数偏移；全局细节点编号=偏移+本地编号。
    ! 仅复制 donor 的内部节点；receiver 的交集在 0 或 n+1 外圈，包含对角迁移需要的位置。

    integer(kind=4), intent(in) :: nq, ni, nj, targetIndexX, targetIndexY, mi, mj, sourceIndexX, sourceIndexY
    real(kind=8) :: receiver(0:ni+1,0:nj+1,0:nq-1), donor(0:mi+1,0:mj+1,0:nq-1)
    integer(kind=4) :: i, j, a

    ! 全局细节点编号 = 偏移量 + 本地编号。只取 receiver 外圈与 donor 内部的交集。

    !$acc parallel loop collapse(3) present(receiver,donor) async(1)
    do a=0,nq-1
        do j=max(0,sourceIndexY+1-targetIndexY),min(nj+1,sourceIndexY+mj-targetIndexY)
            do i=max(0,sourceIndexX+1-targetIndexX),min(ni+1,sourceIndexX+mi-targetIndexX)
                receiver(i,j,a)=donor(targetIndexX+i-sourceIndexX,targetIndexY+j-sourceIndexY,a)
            enddo
        enddo
    enddo
  end subroutine copy_fine_halo
!===================================================================================================

!===================================================================================================
! 子程序: advance_multiblock
! 作用: 完成粗网格预测、细网格子步、粗细同步及时间历史滚动。
!===================================================================================================
  subroutine advance_multiblock()
    use commondata
    implicit none
    integer(kind=4) :: k
    real(kind=8) :: timeInterpWeight(0:2)

    if(refineRatio==1) then
        call advance_coarse()
        itc=itc+1
        return
    endif

    ! 粗网格预测终点；细子步读取过去、起点和预测终点三个时间层。
    call advance_coarse()
    call save_coarse_buffer(2)
    do k=1,refineRatio
        call coarse_time_weights(dble(k-1)/dble(refineRatio),timeInterpWeight)
        call coarse_to_fine(timeInterpWeight)
        call advance_fine()
    enddo

    ! 同步时直接读细网格修正粗接收层，再补细接收层，最后滚动粗时间层。
    call fine_to_coarse()
    call save_coarse_buffer(2)
    call coarse_to_fine([0.0d0,0.0d0,1.0d0])
    call rotate_coarse_history()
    itc=itc+refineRatio
  end subroutine advance_multiblock
!===================================================================================================

!===================================================================================================
! 子程序: coarse_time_weights
! 作用: 计算细子步起始时刻在粗网格时间历史中的插值权重。
!===================================================================================================
  subroutine coarse_time_weights(theta,timeInterpWeight)

    use commondata,only: itc,reloadItcOffset
    implicit none

    ! theta=(细子步起点-t)/dtCoarse，范围从 0 到 (refineRatio-1)/refineRatio。
    ! timeInterpWeight(0:2) 分别作用于 t-dtCoarse、t、t+dtCoarse；启动时没有过去层，只用当前和预测层。

    real(kind=8), intent(in) :: theta
    real(kind=8), intent(out) :: timeInterpWeight(0:2)

    ! theta 是当前细子步起点相对粗步起点的时间，不能使用子步终点时间。
    if(reloadItcOffset+itc==0) then
        timeInterpWeight=[0.0d0,1.0d0-theta,theta]    ! 缺少负时间层时线性启动
    else
        timeInterpWeight=[0.5d0*theta*(theta-1.0d0),1.0d0-theta**2,0.5d0*theta*(theta+1.0d0)]
    endif
  end subroutine coarse_time_weights
!===================================================================================================

!===================================================================================================
! 子程序: rotate_coarse_history
! 作用: 在粗细同步后滚动粗历史，准备下一轮时间插值。
!===================================================================================================
  subroutine rotate_coarse_history()
    use commondata
    implicit none
    integer(kind=4) :: i, j, a

    !$acc parallel loop collapse(3) present(coarseBuffer) async(1)
    do a=0,19
        do j=1,size(coarseBuffer,2)
            do i=1,size(coarseBuffer,1)
                coarseBuffer(i,j,a,0)=coarseBuffer(i,j,a,1)
                coarseBuffer(i,j,a,1)=coarseBuffer(i,j,a,2)
            enddo
        enddo
    enddo
  end subroutine rotate_coarse_history
!===================================================================================================

!===================================================================================================
! 子程序: lagrange_weights
! 作用: 计算相邻四个等距节点的拉格朗日插值权重。
!===================================================================================================
  subroutine lagrange_weights(q,w)

    implicit none
    !$acc routine seq

    ! q 为目标点相对于模板首节点的距离，以来源网格格距为单位。
    ! 模板的相对坐标是 0、1、2、3；w(1:4) 是三次多项式的四点插值权重。

    real(kind=8), intent(in) :: q
    real(kind=8), intent(out) :: w(4)
    integer(kind=4) :: a, k

    w=1.0d0
    do a=0,3
        do k=0,3
            if(a/=k) w(a+1)=w(a+1)*(q-dble(k))/dble(a-k)
        enddo
    enddo
  end subroutine lagrange_weights
!===================================================================================================

!===================================================================================================
! 子程序: section_weights
! 作用: 定位中线诊断的四点模板，计算函数值及导数的插值权重。
!===================================================================================================
  subroutine section_weights(q,n,first,w,dw)

    implicit none

    ! q 为连续节点下标，n 为该方向节点数；first 是四点模板的首下标。
    ! w 用于函数值，dw 是对节点下标的导数权重，换成物理导数时还要除以 dx。

    real(kind=8), intent(in) :: q
    integer(kind=4), intent(in) :: n
    integer(kind=4), intent(out) :: first
    real(kind=8), intent(out) :: w(4), dw(4)
    real(kind=8) :: z, term
    integer(kind=4) :: a, k, l

    first=max(1,min(floor(q)-1,n-3))
    z=q-dble(first)

    call lagrange_weights(z,w)
    dw=0.0d0
    do a=0,3
        do k=0,3
            if(k==a) cycle
            term=1.0d0/dble(a-k)
            do l=0,3
                if(l/=a .AND. l/=k) term=term*(z-dble(l))/dble(a-l)
            enddo

            dw(a+1)=dw(a+1)+term
        enddo
    enddo
  end subroutine section_weights
!===================================================================================================

!===================================================================================================
! 子程序: equilibrium_moments
! 作用: 根据密度、速度和温度，构造接口重建所需的平衡矩。
!===================================================================================================
  subroutine equilibrium_moments(rh,ux,uy,temp,meq,neq)

    use commondata,only: paraA
    implicit none

    !$acc routine seq
    real(kind=8), intent(in) :: rh, ux, uy, temp
    real(kind=8), intent(out) :: meq(0:8), neq(0:4)

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
  end subroutine equilibrium_moments
!===================================================================================================

!===================================================================================================
! 子程序: force_moments
! 作用: 将宏观力转换成流场矩空间中的力项。
!===================================================================================================
  subroutine force_moments(ux,uy,fx,fy,fm)

    implicit none

    !$acc routine seq
    real(kind=8), intent(in) :: ux, uy, fx, fy
    real(kind=8), intent(out) :: fm(0:8)

    fm=[0.0d0,6.0d0*(ux*fx+uy*fy),-6.0d0*(ux*fx+uy*fy),fx,-fx,fy,-fy,&
        2.0d0*(ux*fx-uy*fy),ux*fy+uy*fx]
  end subroutine force_moments
!===================================================================================================

!===================================================================================================
! 子程序: save_coarse_buffer
! 作用: 保存指定时刻的宏观量、统一时间尺度的力和缩放非平衡矩。
!===================================================================================================
  subroutine save_coarse_buffer(slot)
    use commondata
    implicit none
    integer(kind=4), intent(in) :: slot
    integer(kind=4) :: i, j, a
    real(kind=8) :: fv(0:8), gv(0:4), buffer(0:19)

    !$acc parallel loop collapse(2) present(f_coarse, g_coarse, rho_coarse, u_coarse, v_coarse, &
    !$acc& T_coarse, Fx_coarse, Fy_coarse, coarseBuffer) async(1) private(a, fv, gv, buffer) &
    !$acc& firstprivate(dxCoarse, SnuCoarse, SqCoarse, QkCoarse, QnuCoarse)
    do j=1,size(f_coarse,2)
        do i=1,size(f_coarse,1)
            do a=0,8
                fv(a)=f_coarse(i,j,a)
            enddo
            do a=0,4
                gv(a)=g_coarse(i,j,a)
            enddo
            call pack_interface_values(dxCoarse,SnuCoarse,SqCoarse,QkCoarse,QnuCoarse,fv,gv,&
                rho_coarse(i,j),u_coarse(i,j),v_coarse(i,j),T_coarse(i,j),&
                Fx_coarse(i,j),Fy_coarse(i,j),buffer)
            do a=0,19
                coarseBuffer(i,j,a,slot)=buffer(a)
            enddo
        enddo
    enddo
  end subroutine save_coarse_buffer
!===================================================================================================

!===================================================================================================
! 子程序: pack_interface_values
! 作用: 把一个来源节点的宏观量、单位时间力和缩放非平衡矩放入局部缓冲数组。
!===================================================================================================
  subroutine pack_interface_values(dx,Snu,Sq,Qk,Qnu,fv,gv,rh,ux,uy,temp,fx,fy,buffer)
    implicit none
    !$acc routine seq
    real(kind=8), intent(in) :: dx, Snu, Sq, Qk, Qnu, fv(0:8), gv(0:4), rh, ux, uy, temp, fx, fy
    real(kind=8), intent(out) :: buffer(0:19)
    real(kind=8) :: m(0:8), meq(0:8), fm(0:8), n(0:4), neq(0:4), s(0:8), q(0:4)
    integer(kind=4) :: a

    call flow_moments(fv,m)
    call thermal_moments(gv,n)
    call equilibrium_moments(rh,ux,uy,temp,meq,neq)
    call force_moments(ux,uy,fx,fy,fm)
    s=[0.0d0,Snu,Snu,0.0d0,Sq,0.0d0,Sq,Snu,Snu]
    q=[0.0d0,Qk,Qk,Qnu,Qnu]
    buffer(0:5)=[rh,ux,uy,temp,fx/dx,fy/dx]
    do a=0,8
        buffer(6+a)=s(a)/dx*(m(a)-meq(a)+0.5d0*fm(a))
    enddo
    do a=0,4
        buffer(15+a)=q(a)/dx*(n(a)-neq(a))
    enddo
  end subroutine pack_interface_values
!===================================================================================================

!===================================================================================================
! 子程序: unpack_interface_values
! 作用: 用接收网格的松弛率和格距，把局部缓冲值还原成分布函数及宏观量。
!===================================================================================================
  subroutine unpack_interface_values(dx,Snu,Sq,Qk,Qnu,buffer,fv,gv,rh,ux,uy,temp,fx,fy)
    implicit none
    !$acc routine seq
    real(kind=8), intent(in) :: dx, Snu, Sq, Qk, Qnu, buffer(0:19)
    real(kind=8), intent(out) :: fv(0:8), gv(0:4), rh, ux, uy, temp, fx, fy
    real(kind=8) :: m(0:8), meq(0:8), fm(0:8), n(0:4), neq(0:4), s(0:8), q(0:4)
    integer(kind=4) :: a

    rh=buffer(0)
    ux=buffer(1)
    uy=buffer(2)
    temp=buffer(3)
    fx=dx*buffer(4)
    fy=dx*buffer(5)
    call equilibrium_moments(rh,ux,uy,temp,meq,neq)
    call force_moments(ux,uy,fx,fy,fm)
    s=[0.0d0,Snu,Snu,0.0d0,Sq,0.0d0,Sq,Snu,Snu]
    q=[0.0d0,Qk,Qk,Qnu,Qnu]
    m=meq
    n=neq
    do a=0,8
        if(s(a)>0.0d0) m(a)=meq(a)+dx/s(a)*buffer(6+a)-0.5d0*fm(a)
    enddo
    ! 守恒矩由宏观量和半步力重建，不除以零松弛率。
    m(0)=rh
    m(3)=rh*ux-0.5d0*fx
    m(5)=rh*uy-0.5d0*fy
    do a=1,4
        n(a)=neq(a)+dx/q(a)*buffer(15+a)
    enddo
    n(0)=temp
    call flow_populations(m,fv)
    call thermal_populations(n,gv)
  end subroutine unpack_interface_values
!===================================================================================================

!===================================================================================================
! 子程序: coarse_to_fine
! 作用: 四个细区直接从粗时间层插值，结果写入各自的接收节点。
!===================================================================================================
  subroutine coarse_to_fine(timeInterpWeight)
    use commondata
    implicit none
    real(kind=8), intent(in) :: timeInterpWeight(0:2)

    call interpolate_coarse_to_fine(size(f_left,1),size(f_left,2),xOffsetLeft,yOffsetLeft,&
        f_left,g_left,rho_left,u_left,v_left,T_left,Fx_left,Fy_left,timeInterpWeight)
    call interpolate_coarse_to_fine(size(f_right,1),size(f_right,2),xOffsetRight,yOffsetRight,&
        f_right,g_right,rho_right,u_right,v_right,T_right,Fx_right,Fy_right,timeInterpWeight)
    call interpolate_coarse_to_fine(size(f_bottom,1),size(f_bottom,2),xOffsetBottom,yOffsetBottom,&
        f_bottom,g_bottom,rho_bottom,u_bottom,v_bottom,T_bottom,Fx_bottom,Fy_bottom,timeInterpWeight)
    call interpolate_coarse_to_fine(size(f_top,1),size(f_top,2),xOffsetTop,yOffsetTop,&
        f_top,g_top,rho_top,u_top,v_top,T_top,Fx_top,Fy_top,timeInterpWeight)
  end subroutine coarse_to_fine
!===================================================================================================

!===================================================================================================
! 子程序: interpolate_coarse_to_fine
! 作用: 在实际细接收层内逐点计算来源下标和权重，无需保存接口列表或插值权重表。
!===================================================================================================
  subroutine interpolate_coarse_to_fine(ni,nj,xOffset,yOffset,f,g,rho,u,v,T,Fx,Fy,timeInterpWeight)
    use commondata,only: coarseBuffer,dxCoarse,xOffsetCoarse,yOffsetCoarse,centerIntegrationBox,&
        interfaceSkin,Snu,Sq,Qk,Qnu
    implicit none
    integer(kind=4), intent(in) :: ni, nj
    real(kind=8), intent(in) :: xOffset, yOffset, timeInterpWeight(0:2)
    real(kind=8), intent(inout) :: f(ni,nj,0:8), g(ni,nj,0:4), rho(ni,nj), u(ni,nj), v(ni,nj)
    real(kind=8), intent(inout) :: T(ni,nj), Fx(ni,nj), Fy(ni,nj)
    integer(kind=4) :: i, j, a, k, ix, iy, sourceX, sourceY, iStart, iEnd, jStart, jEnd, mi, mj
    real(kind=8) :: qx, qy, wx(4), wy(4), wk, buffer(0:19), fv(0:8), gv(0:4)
    real(kind=8) :: timeWeight0, timeWeight1, timeWeight2, x0, y0, dx
    logical :: samePosition

    ! 四个紧凑细矩形在中心物理矩形内的部分，恰好是要重建的接收层。
    iStart=max(1,floor(centerIntegrationBox(1)-xOffset+0.5d0)+1)
    iEnd=min(ni,ceiling(centerIntegrationBox(2)-xOffset+0.5d0)-1)
    jStart=max(1,floor(centerIntegrationBox(3)-yOffset+0.5d0)+1)
    jEnd=min(nj,ceiling(centerIntegrationBox(4)-yOffset+0.5d0)-1)
    mi=size(coarseBuffer,1)
    mj=size(coarseBuffer,2)
    x0=xOffsetCoarse
    y0=yOffsetCoarse
    dx=dxCoarse
    ! 保留 NVHPC 所需的标量捕获，不能改为 firstprivate(timeInterpWeight) 数组形参。
    timeWeight0=timeInterpWeight(0)
    timeWeight1=timeInterpWeight(1)
    timeWeight2=timeInterpWeight(2)
    !$acc parallel loop collapse(2) present(coarseBuffer, f, g, rho, u, v, T, Fx, Fy) async(1) &
    !$acc& firstprivate(timeWeight0, timeWeight1, timeWeight2, x0, y0, dx, mi, mj) &
    !$acc& private(a, k, ix, iy, sourceX, sourceY, qx, qy, wx, wy, wk, buffer, fv, gv, samePosition)
    do j=jStart,jEnd
        do i=iStart,iEnd
            qx=(xOffset+(i-0.5d0)-x0)/dx+0.5d0
            qy=(yOffset+(j-0.5d0)-y0)/dx+0.5d0
            samePosition=abs(qx-nint(qx))<1.0d-12 .AND. abs(qy-nint(qy))<1.0d-12
            if(samePosition) then
                sourceX=nint(qx)
                sourceY=nint(qy)
            else
                sourceX=max(1+interfaceSkin,min(floor(qx)-1,mi-interfaceSkin-3))
                sourceY=max(1+interfaceSkin,min(floor(qy)-1,mj-interfaceSkin-3))
                call lagrange_weights(qx-sourceX,wx)
                call lagrange_weights(qy-sourceY,wy)
            endif
            do a=0,19
                buffer(a)=0.0d0
                do k=0,2
                    select case(k)
                    case(0)
                        wk=timeWeight0
                    case(1)
                        wk=timeWeight1
                    case(2)
                        wk=timeWeight2
                    end select
                    if(samePosition) then
                        buffer(a)=buffer(a)+wk*coarseBuffer(sourceX,sourceY,a,k)
                    else
                        do iy=1,4
                            do ix=1,4
                                buffer(a)=buffer(a)+wk*wx(ix)*wy(iy)*&
                                    coarseBuffer(sourceX+ix-1,sourceY+iy-1,a,k)
                            enddo
                        enddo
                    endif
                enddo
            enddo
            call unpack_interface_values(1.0d0,Snu,Sq,Qk,Qnu,buffer,fv,gv,&
                rho(i,j),u(i,j),v(i,j),T(i,j),Fx(i,j),Fy(i,j))
            do a=0,8
                f(i,j,a)=fv(a)
            enddo
            do a=0,4
                g(i,j,a)=gv(a)
            enddo
        enddo
    enddo
  end subroutine interpolate_coarse_to_fine
!===================================================================================================

!===================================================================================================
! 子程序: fine_to_coarse
! 作用: 直接读取四个细区的当前场，重建粗网格接收层。
!===================================================================================================
  subroutine fine_to_coarse()
    use commondata
    implicit none

    call copy_fine_to_coarse(size(f_left,1),size(f_left,2),xOffsetLeft,yOffsetLeft,&
        f_left,g_left,rho_left,u_left,v_left,T_left,Fx_left,Fy_left)
    call copy_fine_to_coarse(size(f_right,1),size(f_right,2),xOffsetRight,yOffsetRight,&
        f_right,g_right,rho_right,u_right,v_right,T_right,Fx_right,Fy_right)
    call copy_fine_to_coarse(size(f_bottom,1),size(f_bottom,2),xOffsetBottom,yOffsetBottom,&
        f_bottom,g_bottom,rho_bottom,u_bottom,v_bottom,T_bottom,Fx_bottom,Fy_bottom)
    call copy_fine_to_coarse(size(f_top,1),size(f_top,2),xOffsetTop,yOffsetTop,&
        f_top,g_top,rho_top,u_top,v_top,T_top,Fx_top,Fy_top)
  end subroutine fine_to_coarse
!===================================================================================================

!===================================================================================================
! 子程序: coarse_interface_node
! 作用: 把粗外圈的循环编号换成下标，顺序为下边、上边及左右两边，不保存节点表。
!===================================================================================================
  subroutine coarse_interface_node(c,ni,nj,i,j)
    use commondata,only: interfaceSkin
    implicit none
    !$acc routine seq
    integer(kind=4), intent(in) :: c, ni, nj
    integer(kind=4), intent(out) :: i, j
    integer(kind=4) :: n

    if(c<=ni*interfaceSkin) then
        i=mod(c-1,ni)+1
        j=(c-1)/ni+1
    else if(c<=2*ni*interfaceSkin) then
        n=c-ni*interfaceSkin-1
        i=mod(n,ni)+1
        j=nj-interfaceSkin+n/ni+1
    else
        n=c-2*ni*interfaceSkin-1
        i=mod(n,2*interfaceSkin)+1
        if(i>interfaceSkin) i=i+ni-2*interfaceSkin
        j=interfaceSkin+n/(2*interfaceSkin)+1
    endif
  end subroutine coarse_interface_node
!===================================================================================================

!===================================================================================================
! 子程序: copy_fine_to_coarse
! 作用: 粗接收点与细节点共址，直接从对应细节点构造局部缓冲并写入粗数组。
!===================================================================================================
  subroutine copy_fine_to_coarse(ni,nj,xOffset,yOffset,f,g,rho,u,v,T,Fx,Fy)
    use commondata,only: f_coarse,g_coarse,rho_coarse,u_coarse,v_coarse,T_coarse,Fx_coarse,Fy_coarse,&
        dxCoarse,xOffsetCoarse,yOffsetCoarse,interfaceSkin,Snu,Sq,Qk,Qnu,SnuCoarse,SqCoarse,QkCoarse,QnuCoarse
    implicit none
    integer(kind=4), intent(in) :: ni, nj
    real(kind=8), intent(in) :: xOffset, yOffset, f(ni,nj,0:8), g(ni,nj,0:4)
    real(kind=8), intent(in) :: rho(ni,nj), u(ni,nj), v(ni,nj), T(ni,nj), Fx(ni,nj), Fy(ni,nj)
    integer(kind=4) :: c, i, j, a, sourceX, sourceY, mi, mj
    real(kind=8) :: fv(0:8), gv(0:4), buffer(0:19), x0, y0, dx

    mi=size(f_coarse,1)
    mj=size(f_coarse,2)
    x0=xOffsetCoarse
    y0=yOffsetCoarse
    dx=dxCoarse
    !$acc parallel loop present(f, g, rho, u, v, T, Fx, Fy, f_coarse, g_coarse, rho_coarse, &
    !$acc& u_coarse, v_coarse, T_coarse, Fx_coarse, Fy_coarse) async(1) &
    !$acc& firstprivate(mi, mj, x0, y0, dx, SnuCoarse, SqCoarse, QkCoarse, QnuCoarse) &
    !$acc& private(i, j, a, sourceX, sourceY, fv, gv, buffer)
    do c=1,2*interfaceSkin*(mi+mj-2*interfaceSkin)
        call coarse_interface_node(c,mi,mj,i,j)
        sourceX=nint((x0+(i-0.5d0)*dx)-xOffset+0.5d0)
        sourceY=nint((y0+(j-0.5d0)*dx)-yOffset+0.5d0)
        if(sourceX<1 .OR. sourceX>ni .OR. sourceY<1 .OR. sourceY>nj) cycle
        do a=0,8
            fv(a)=f(sourceX,sourceY,a)
        enddo
        do a=0,4
            gv(a)=g(sourceX,sourceY,a)
        enddo
        call pack_interface_values(1.0d0,Snu,Sq,Qk,Qnu,fv,gv,rho(sourceX,sourceY),u(sourceX,sourceY),&
            v(sourceX,sourceY),T(sourceX,sourceY),Fx(sourceX,sourceY),Fy(sourceX,sourceY),buffer)
        call unpack_interface_values(dx,SnuCoarse,SqCoarse,QkCoarse,QnuCoarse,buffer,fv,gv,&
            rho_coarse(i,j),u_coarse(i,j),v_coarse(i,j),T_coarse(i,j),Fx_coarse(i,j),Fy_coarse(i,j))
        do a=0,8
            f_coarse(i,j,a)=fv(a)
        enddo
        do a=0,4
            g_coarse(i,j,a)=gv(a)
        enddo
    enddo
  end subroutine copy_fine_to_coarse
!===================================================================================================

!===================================================================================================
! 子程序: check_interfaces
! 作用: 初始化时核对直接坐标映射、来源模板及粗接收层覆盖，不分配接口列表。
!===================================================================================================
  subroutine check_interfaces()
    use commondata
    implicit none
    integer(kind=4) :: count

    if(refineRatio==1) return
    count=0
    call check_interface_grid(size(f_left,1),size(f_left,2),xOffsetLeft,yOffsetLeft,count)
    call check_interface_grid(size(f_right,1),size(f_right,2),xOffsetRight,yOffsetRight,count)
    call check_interface_grid(size(f_bottom,1),size(f_bottom,2),xOffsetBottom,yOffsetBottom,count)
    call check_interface_grid(size(f_top,1),size(f_top,2),xOffsetTop,yOffsetTop,count)
    if(count/=2*interfaceSkin*(size(f_coarse,1)+size(f_coarse,2)-2*interfaceSkin)) &
        error stop 'Four fine arrays do not cover the complete coarse interface'
  end subroutine check_interfaces
!===================================================================================================

!===================================================================================================
! 子程序: check_interface_grid
! 作用: 检查一个细区的插值范围及其提供的粗接收节点。
!===================================================================================================
  subroutine check_interface_grid(ni,nj,xOffset,yOffset,count)
    use commondata,only: f_coarse,dxCoarse,xOffsetCoarse,yOffsetCoarse,interfaceSkin,centerIntegrationBox
    implicit none
    integer(kind=4), intent(in) :: ni, nj
    real(kind=8), intent(in) :: xOffset, yOffset
    integer(kind=4), intent(inout) :: count
    integer(kind=4) :: i, j, c, sourceX, sourceY
    real(kind=8) :: x, y, qx, qy, wx(4), wy(4)
    logical :: isFineSkin, samePosition

    do j=max(1,floor(centerIntegrationBox(3)-yOffset+0.5d0)+1),&
         min(nj,ceiling(centerIntegrationBox(4)-yOffset+0.5d0)-1)
        do i=max(1,floor(centerIntegrationBox(1)-xOffset+0.5d0)+1),&
             min(ni,ceiling(centerIntegrationBox(2)-xOffset+0.5d0)-1)
            x=xOffset+i-0.5d0
            y=yOffset+j-0.5d0
            call fine_skin(x,y,isFineSkin)
            if(.NOT.isFineSkin) error stop 'Fine receiver rectangle exceeds interface layers'
            call coarse_donor_stencil(x,y,sourceX,sourceY,wx,wy,samePosition)
        enddo
    enddo
    do c=1,2*interfaceSkin*(size(f_coarse,1)+size(f_coarse,2)-2*interfaceSkin)
        call coarse_interface_node(c,size(f_coarse,1),size(f_coarse,2),i,j)
        x=xOffsetCoarse+(i-0.5d0)*dxCoarse
        y=yOffsetCoarse+(j-0.5d0)*dxCoarse
        qx=x-xOffset+0.5d0
        qy=y-yOffset+0.5d0
        sourceX=nint(qx)
        sourceY=nint(qy)
        if(sourceX<1 .OR. sourceX>ni .OR. sourceY<1 .OR. sourceY>nj) cycle
        call fine_skin(x,y,isFineSkin)
        if(abs(qx-sourceX)>1.0d-12 .OR. abs(qy-sourceY)>1.0d-12 .OR. isFineSkin) &
            error stop 'Coarse receiver must coincide with an interior fine donor'
        count=count+1
    enddo
  end subroutine check_interface_grid
!===================================================================================================

!===================================================================================================
! 子程序: fine_skin
! 作用: 按全局坐标判断细节点是否处于细环内侧的人工边缘层。
!===================================================================================================
  subroutine fine_skin(x,y,isInterfaceNode)
    use commondata,only: centerIntegrationBox,interfaceSkin
    implicit none

    logical, intent(out) :: isInterfaceNode    ! 输出：该细节点是否属于需要接收粗网格数据的人工边缘层。

    ! 这里检查的是整个细环的内缘；四个细矩形之间的存储接缝不属于粗细人工边界。
    ! 细接收节点位于中心分区线内侧的 interfaceSkin 层，按细格距定位。

    real(kind=8), intent(in) :: x, y

    isInterfaceNode=(x <=centerIntegrationBox(1)+dble(interfaceSkin) .OR. &
        x >=centerIntegrationBox(2)-dble(interfaceSkin) .OR. &
        y <=centerIntegrationBox(3)+dble(interfaceSkin) .OR. &
        y >=centerIntegrationBox(4)-dble(interfaceSkin)) .AND. &
        x > centerIntegrationBox(1) .AND. &
        x < centerIntegrationBox(2) .AND. &
        y > centerIntegrationBox(3) .AND. &
        y < centerIntegrationBox(4)
  end subroutine fine_skin
!===================================================================================================

!===================================================================================================
! 子程序: coarse_donor_stencil
! 作用: 根据目标坐标选择粗来源节点或四点模板，并检查来源范围。
!===================================================================================================
  subroutine coarse_donor_stencil(x,y,sourceIndexX,sourceIndexY,spaceInterpWeightX,spaceInterpWeightY,samePosition)
    use commondata,only: f_coarse,xOffsetCoarse,yOffsetCoarse,dxCoarse,interfaceSkin
    implicit none

    ! continuousIndexX/Y = (坐标-Offset)/dxCoarse+0.5，是坐标对应的连续数组下标。
    ! 例如默认 x=127.5 对应粗下标 3；x=128.5 对应 3.5，需要在附近粗节点间插值。
    ! 来源范围为 1+interfaceSkin 到 nCoarse-interfaceSkin，排除四周接收层。
    ! 共址时直接记录整数节点；非共址时记录四点首下标与两个方向权重。

    real(kind=8), intent(in) :: x, y
    integer(kind=4), intent(out) :: sourceIndexX, sourceIndexY
    real(kind=8), intent(out) :: spaceInterpWeightX(4), spaceInterpWeightY(4)
    logical, intent(out) :: samePosition
    real(kind=8) :: continuousIndexX, continuousIndexY

    continuousIndexX=(x-xOffsetCoarse)/dxCoarse+.5d0
    continuousIndexY=(y-yOffsetCoarse)/dxCoarse+.5d0
    if(continuousIndexX < (1+interfaceSkin) .OR. continuousIndexX > (size(f_coarse,1)-interfaceSkin) .OR. &
        continuousIndexY < (1+interfaceSkin) .OR. continuousIndexY > (size(f_coarse,2)-interfaceSkin) .OR. &
        min((size(f_coarse,1)-interfaceSkin)-(1+interfaceSkin),(size(f_coarse,2)-interfaceSkin)-(1+interfaceSkin)) < 3) &
        error stop 'Invalid overlap: no interior four-point coarse donor stencil'
    samePosition=abs(continuousIndexX-nint(continuousIndexX)) < 1d-12 .AND. abs(continuousIndexY-&
        nint(continuousIndexY)) < 1d-12
    if(samePosition) then
        sourceIndexX=nint(continuousIndexX)
        sourceIndexY=nint(continuousIndexY)
        spaceInterpWeightX=[1d0,0d0,0d0,0d0]
        spaceInterpWeightY=spaceInterpWeightX
    else
        sourceIndexX=max((1+interfaceSkin),min(floor(continuousIndexX)-1,(size(f_coarse,1)-interfaceSkin)-3))
        sourceIndexY=max((1+interfaceSkin),min(floor(continuousIndexY)-1,(size(f_coarse,2)-interfaceSkin)-3))

        call lagrange_weights(continuousIndexX-sourceIndexX,spaceInterpWeightX)

        call lagrange_weights(continuousIndexY-sourceIndexY,spaceInterpWeightY)
    endif
  end subroutine coarse_donor_stencil
!===================================================================================================

!===================================================================================================
! 子程序: collision
! 作用: 执行 D2Q9 流场矩空间碰撞，生成碰撞后分布函数。
!===================================================================================================
  subroutine collision(nx,ny,f,f_post,rho,u,v,Fx,Fy,T,Snu,Sq,gBeta &
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

    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), intent(inout) :: f(nx,ny,0:8), f_post(0:nx+1,0:ny+1,0:8)
    real(kind=8), intent(inout) :: u(nx,ny), v(nx,ny), T(nx,ny), rho(nx,ny), Fx(nx,ny), Fy(nx,ny)
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
    do j=1,ny
        do i=1,nx

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

            fSource(0)=0.0d0    !将源项F对应的贡献投影到各个矩中，并做半步修正
            fSource(1)=(6.0d0-3.0d0*s(1))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
            fSource(2)=-(6.0d0-3.0d0*s(2))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
            fSource(3)=(1.0d0-0.5d0*s(3))*Fx(i,j)
            fSource(4)=-(1.0d0-0.5d0*s(4))*Fx(i,j)
            fSource(5)=(1.0d0-0.5d0*s(5))*Fy(i,j)
            fSource(6)=-(1.0d0-0.5d0*s(6))*Fy(i,j)
            fSource(7)=(2.0d0-s(7))*(u(i,j)*Fx(i,j)-v(i,j)*Fy(i,j))
            fSource(8)=(1.0d0-0.5d0*s(8))*(u(i,j)*Fy(i,j)+v(i,j)*Fx(i,j))    !这边是乘以M变到矩空间，然后再乘以1-1/2S修正

            do alpha=0,8
                m_post(alpha)=m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+fSource(alpha)    !矩空间碰撞
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

!===================================================================================================
! 子程序: streaming
! 作用: 按照格子方向，从上游节点迁移流场分布函数。
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

!===================================================================================================
! 子程序: bounceback
! 作用: 对接触真实物理墙面的边界执行流场边界处理。
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

#ifdef VerticalWallsPeriodicalU

    !$acc parallel loop gang vector present(f, f_post) async(1)
    do j=1,ny    !速度边界垂直边界周期，直接方向相同，跨边界的入射分布
        !Left side (i=1)
        if(leftWall) f(1,j,1)=f_post(nx,j,1)
        if(leftWall) f(1,j,5)=f_post(nx,j,5)
        if(leftWall) f(1,j,8)=f_post(nx,j,8)

        !Right side (i=nx)
        if(rightWall) f(nx,j,3)=f_post(1,j,3)
        if(rightWall) f(nx,j,6)=f_post(1,j,6)
        if(rightWall) f(nx,j,7)=f_post(1,j,7)
    enddo
#endif

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

!===================================================================================================
! 子程序: macro
! 作用: 由流场分布函数及半步力恢复密度和速度。
!===================================================================================================
  subroutine macro(nx,ny,f,rho,u,v,Fx,Fy)

    implicit none

    ! 由分布函数求密度和速度，速度包含半步力修正。

    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), intent(in) :: f(nx,ny,0:8), Fx(nx,ny), Fy(nx,ny)
    real(kind=8), intent(inout) :: rho(nx,ny), u(nx,ny), v(nx,ny)

    integer(kind=4) :: i, j

    !$acc parallel loop gang vector collapse(2) present(f, rho, u, v, Fx, Fy) async(1)
    do j=1,ny
        do i=1,nx
            rho(i,j)=f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,&
                6)+f(i,j,&
                7)+f(i,j,8)
            u(i,j)=(f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)+0.5d0*Fx(i,&
                j))/rho(i,j)    !含力LBM的半步动量修正：rho*u = Σ f e + 0.5*F，对应Guo forcing的二阶定义
            v(i,j)=(f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)+0.5d0*Fy(i,&
                j))/rho(i,j)
        enddo
    enddo

    return
  end subroutine macro
!===================================================================================================

!===================================================================================================
! 子程序: collisionT
! 作用: 执行原 D2Q5 温度碰撞，生成碰撞后温度分布函数。
!===================================================================================================
  subroutine collisionT(nx,ny,g,g_post,u,v,T,Qk,Qnu)

    use commondata,only: paraA
    implicit none

    ! Qk/Qnu 已按当前区域格距设置，温度碰撞使用本步更新后的速度。

    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), intent(inout) :: g(nx,ny,0:4), g_post(0:nx+1,0:ny+1,0:4)
    real(kind=8), intent(in) :: u(nx,ny), v(nx,ny), T(nx,ny), Qk, Qnu

    integer(kind=4) :: i, j
    real(kind=8) :: n(0:4), n_post(0:4), neq(0:4)
    real(kind=8) :: q(0:4)

    !$acc parallel loop gang vector collapse(2) present(g, g_post, u, v, T) async(1) &
    !$acc& private(n, neq, q, n_post)
    do j=1,ny
        do i=1,nx

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

!===================================================================================================
! 子程序: streamingT
! 作用: 按照格子方向，从上游节点迁移温度分布函数。
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

!===================================================================================================
! 子程序: bouncebackT
! 作用: 对接触真实物理墙面的边界执行温度边界处理。
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

#ifdef VerticalWallsPeriodicalT

    !$acc parallel loop gang vector present(g, g_post) async(1)
    do j=1,ny
        !Left boundary
        if(leftWall) g(1,j,1)=g_post(nx,j,1)

        !Right boundary
        if(rightWall) g(nx,j,3)=g_post(1,j,3)
    enddo
#endif

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

!===================================================================================================
! 子程序: macroT
! 作用: 由温度分布函数恢复宏观温度。
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

!===================================================================================================
! 子程序: flow_moments
! 作用: 将节点的九个流场分布函数转换为九个矩分量。
!===================================================================================================
  subroutine flow_moments(fv,m)

    implicit none

    !$acc routine seq
    real(kind=8), intent(in) :: fv(0:8)
    real(kind=8), intent(out) :: m(0:8)

    m(0)=fv(0)+fv(1)+fv(2)+fv(3)+fv(4)+fv(5)+fv(6)+fv(7)+fv(8)
    m(1)=-4.0d0*fv(0)-fv(1)-fv(2)-fv(3)-fv(4)+2.0d0*(fv(5)+fv(6)+fv(7)+fv(8))
    m(2)=4.0d0*fv(0)-2.0d0*(fv(1)+fv(2)+fv(3)+fv(4))+fv(5)+fv(6)+fv(7)+fv(8)
    m(3)=fv(1)-fv(3)+fv(5)-fv(6)-fv(7)+fv(8)
    m(4)=-2.0d0*fv(1)+2.0d0*fv(3)+fv(5)-fv(6)-fv(7)+fv(8)
    m(5)=fv(2)-fv(4)+fv(5)+fv(6)-fv(7)-fv(8)
    m(6)=-2.0d0*fv(2)+2.0d0*fv(4)+fv(5)+fv(6)-fv(7)-fv(8)
    m(7)=fv(1)-fv(2)+fv(3)-fv(4)
    m(8)=fv(5)-fv(6)+fv(7)-fv(8)

  end subroutine flow_moments
!===================================================================================================

!===================================================================================================
! 子程序: flow_populations
! 作用: 将流场矩分量逆变换为九个方向的分布函数。
!===================================================================================================
  subroutine flow_populations(m,fv)

    implicit none

    !$acc routine seq
    real(kind=8), intent(in) :: m(0:8)
    real(kind=8), intent(out) :: fv(0:8)

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

  end subroutine flow_populations
!===================================================================================================

!===================================================================================================
! 子程序: thermal_moments
! 作用: 将节点的五个温度分布函数转换为五个矩分量。
!===================================================================================================
  subroutine thermal_moments(gv,n)

    implicit none

    !$acc routine seq
    real(kind=8), intent(in) :: gv(0:4)
    real(kind=8), intent(out) :: n(0:4)

    n(0)=gv(0)+gv(1)+gv(2)+gv(3)+gv(4)
    n(1)=gv(1)-gv(3)
    n(2)=gv(2)-gv(4)
    n(3)=-4.0d0*gv(0)+gv(1)+gv(2)+gv(3)+gv(4)
    n(4)=gv(1)-gv(2)+gv(3)-gv(4)

  end subroutine thermal_moments
!===================================================================================================

!===================================================================================================
! 子程序: thermal_populations
! 作用: 将温度矩分量逆变换为五个方向的分布函数。
!===================================================================================================
  subroutine thermal_populations(n,gv)

    implicit none

    !$acc routine seq
    real(kind=8), intent(in) :: n(0:4)
    real(kind=8), intent(out) :: gv(0:4)

    gv(0)=0.2d0*n(0)-0.2d0*n(3)
    gv(1)=0.2d0*n(0)+0.5d0*n(1)+0.05d0*n(3)+0.25d0*n(4)
    gv(2)=0.2d0*n(0)+0.5d0*n(2)+0.05d0*n(3)-0.25d0*n(4)
    gv(3)=0.2d0*n(0)-0.5d0*n(1)+0.05d0*n(3)+0.25d0*n(4)
    gv(4)=0.2d0*n(0)-0.5d0*n(2)+0.05d0*n(3)-0.25d0*n(4)

  end subroutine thermal_populations
!===================================================================================================

!===================================================================================================
! 子程序: output_steady_monitor
! 作用: 在稳态误差检查时记录当前 Nu/Re，不改变正式采样计数和续算历史。
!===================================================================================================
  subroutine output_steady_monitor()
#ifdef steadyFlow
    use commondata,only: itc,reloadItcOffset,timeUnit,errorU,errorT
    implicit none

    integer(kind=4) :: k
    real(kind=8) :: nu, re, hot, cold, middle, mass, meanT, tmin, tmax, rmin, rmax

    call calculate_NuRe(nu,re,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax)

    open(newunit=k,file='SteadyMonitor_2DOpenaccMultiblock.dat',status='unknown',position='append')
    write(k,'(I12,8(1X,ES24.16E3))') reloadItcOffset+itc,dble(reloadItcOffset+itc)/timeUnit,errorU,errorT,nu,re,hot,cold,middle
    close(k)
#endif
  end subroutine output_steady_monitor
!===================================================================================================

!===================================================================================================
! 子程序: output_final_NuRe
! 作用: 独立输出稳态计算最后一步的诊断；达到最大步数不代表已收敛。
!===================================================================================================
  subroutine output_final_NuRe()
#ifdef steadyFlow
    use commondata,only: itc,reloadItcOffset,timeUnit,settingsFile,errorU,errorT,epsU,epsT
    implicit none

    integer(kind=4) :: k
    real(kind=8) :: nu, re, hot, cold, middle, mass, meanT, tmin, tmax, rmin, rmax

    ! 输出最后一步的诊断，不追加采样历史，也不增加 NuReSampleCount。
    call calculate_NuRe(nu,re,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax)

    open(newunit=k,file='FinalNuRe_2DOpenaccMultiblock.dat',status='replace')
    write(k,'(a,L1)') '# converged = ',errorU <=epsU .AND. errorT <=epsT
    write(k,'(a)') '# t_ff NuVolAvg ReVolRMS Nu_hot Nu_cold Nu_middle mass meanT Tmin Tmax rhoMin rhoMax'
    write(k,'(12(ES24.16E3,1X))') dble(reloadItcOffset+itc)/timeUnit,nu,re,hot,cold,middle,&
        mass,meanT,tmin,tmax,rmin,rmax
    close(k)

    open(newunit=k,file=trim(settingsFile),status='unknown',position='append')
    write(k,*) 'Final diagnostics at cumulative fine step:',reloadItcOffset+itc
    write(k,*) 'Steady converged:',errorU <=epsU .AND. errorT <=epsT
    write(k,*) 'errorU, errorT:',errorU,errorT
    write(k,*) 'Nu_global =',nu
    write(k,*) 'Nu_hot =',hot
    write(k,*) 'Nu_cold =',cold
    write(k,*) 'Nu_middle =',middle
    write(k,*) 'Re_rms =',re
    close(k)
    write(*,*) 'Final Nu_global, Re_rms:',nu,re
#endif
  end subroutine output_final_NuRe
!===================================================================================================

!===================================================================================================
! 子程序: calNuRe
! 作用: 按规定的采样时刻写入 Nu/Re 历史，并将正式采样次数增加 1。
!===================================================================================================
  subroutine calNuRe()

    use commondata,only: itc,reloadItcOffset,timeUnit,NuReHistoryFile,NuReSampleCount
    implicit none

    integer(kind=4) :: k
    real(kind=8) :: nu, re, hot, cold, middle, mass, meanT, tmin, tmax, rmin, rmax

    call calculate_NuRe(nu,re,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax)

    open(newunit=k,file=NuReHistoryFile,status='old',position='append')
    write(k,'(12(ES24.16E3,1X))') dble(reloadItcOffset+itc)/timeUnit,nu,re,hot,cold,middle,&
        mass,meanT,tmin,tmax,rmin,rmax
    close(k)
    NuReSampleCount=NuReSampleCount+1
    write(*,'(a,f12.5,a,es13.5,a,es13.5)') 't_ff=',dble(reloadItcOffset+itc)/timeUnit,' NuVolAvg=',nu,' ReVolRMS=',re
  end subroutine calNuRe
!===================================================================================================

!===================================================================================================
! 子程序: calculate_NuRe
! 作用: 只计算当前 Nu/Re 及场量统计，不写采样历史，不改变采样次数。
!===================================================================================================
  subroutine calculate_NuRe(nu,re,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax)

    use commondata
    implicit none

    ! 依次累加五个物理分区，节点实际贡献的面积由 owned_cell_area 计算。
    ! Re 使用速度平方的面积平均再开方；不能直接平均各网格各自计算的 Re。

    real(kind=8), intent(out) :: nu, re, hot, cold, middle, mass, meanT, tmin, tmax, rmin, rmax
    real(kind=8) :: conv, vel2

    call update_host_all(.false.)

    conv=0.0d0
    vel2=0.0d0
    mass=0.0d0
    meanT=0.0d0
    hot=0.0d0
    cold=0.0d0
    middle=0.0d0
    tmin=huge(1.0d0)
    tmax=-tmin
    rmin=tmin
    rmax=-tmin

    call calNuRe_grid(size(f_coarse,1),size(f_coarse,2),dxCoarse,xOffsetCoarse,yOffsetCoarse,(1+merge(interfaceSkin,0,refineRatio>1)),&
        (size(f_coarse,1)-merge(interfaceSkin,0,refineRatio>1)),(1+merge(interfaceSkin,0,refineRatio>1)),(size(f_coarse,2)-merge(interfaceSkin,0,refineRatio>1)),centerIntegrationBox,&
        conv,vel2,mass,meanT,hot,cold,&
        middle,tmin,tmax,rmin,rmax,rho_coarse,u_coarse,v_coarse,T_coarse)
    if(refineRatio > 1) then

        call calNuRe_grid(size(f_left,1),size(f_left,2),dxFine,xOffsetLeft,yOffsetLeft,1,size(f_left,1),&
            1,size(f_left,2),[xOffsetLeft,xOffsetLeft+size(f_left,1),yOffsetLeft,yOffsetLeft+size(f_left,2)],conv,vel2,mass,meanT,hot,cold,&
            middle,tmin,tmax,&
            rmin,rmax,rho_left,u_left,v_left,T_left)

        call calNuRe_grid(size(f_right,1),size(f_right,2),dxFine,xOffsetRight,yOffsetRight,1,size(f_right,1),&
            1,size(f_right,2),[xOffsetRight,xOffsetRight+size(f_right,1),yOffsetRight,yOffsetRight+size(f_right,2)],conv,vel2,mass,&
            meanT,hot,cold,&
            middle,tmin,tmax,&
            rmin,rmax,rho_right,u_right,v_right,T_right)

        call calNuRe_grid(size(f_bottom,1),size(f_bottom,2),dxFine,xOffsetBottom,yOffsetBottom,1,&
            size(f_bottom,1),1,size(f_bottom,2),[xOffsetBottom,xOffsetBottom+size(f_bottom,1),yOffsetBottom,yOffsetBottom+size(f_bottom,2)],conv,&
            vel2,mass,meanT,hot,cold,&
            middle,tmin,tmax,rmin,rmax,rho_bottom,u_bottom,v_bottom,T_bottom)

        call calNuRe_grid(size(f_top,1),size(f_top,2),dxFine,xOffsetTop,yOffsetTop,1,size(f_top,1),&
            1,&
            size(f_top,2),[xOffsetTop,xOffsetTop+size(f_top,1),yOffsetTop,yOffsetTop+size(f_top,2)],conv,vel2,mass,meanT,hot,cold,&
            middle,tmin,tmax,rmin,rmax,&
            rho_top,u_top,v_top,T_top)
    endif

    ! 空间积分等于中心、左、右、下、上五个物理分区之和；接口端点按实际面积计权。
    ! Re 使用速度平方的面积平均再开方，非稳态时间统计也采用同一 RMS 定义。
    nu=1.0d0+conv/(dble(nx)*dble(ny))*(lengthUnit/(Thot-Tcold))/diffusivity
    re=sqrt(vel2/(dble(nx)*dble(ny)))*lengthUnit/viscosity
    hot=(lengthUnit/(Thot-Tcold))*hot
    cold=(lengthUnit/(Thot-Tcold))*cold
    middle=(lengthUnit/(Thot-Tcold))*middle
    meanT=meanT/(dble(nx)*dble(ny))
  end subroutine calculate_NuRe
!===================================================================================================

!===================================================================================================
! 子程序: calNuRe_grid
! 作用: 累计传入区域的体积、壁面及中线统计贡献。
!===================================================================================================
  subroutine calNuRe_grid(ni,nj,dx,xOffset,yOffset,xStartIntegration,xEndIntegration,yStartIntegration,&
      yEndIntegration,integrationBox,conv,&
      vel2,mass,meanT,hot,cold,middle,tmin,tmax,rmin,rmax,&
      rho,u,v,T)

    use commondata,only: nx,ny,Thot,Tcold,diffusivity,itc,reloadItcOffset,refineRatio,ieee_is_finite
    implicit none

    ! integrationBox 为细数组外框或粗区物理矩形；零面积交换节点不参与统计。
    ! 细中线模板跨存储接缝时由 fine_section_x/y 按全局坐标取值，不能截成单侧模板。

    integer(kind=4) :: ni, nj, xStartIntegration, xEndIntegration, yStartIntegration, yEndIntegration
    real(kind=8) :: dx, xOffset, yOffset, integrationBox(4)

    integer(kind=4) :: i, j, jm, im
    real(kind=8) :: rho(ni,nj)
    real(kind=8) :: u(ni,nj)
    real(kind=8) :: v(ni,nj)
    real(kind=8) :: T(ni,nj)
    real(kind=8) :: conv, vel2, mass, meanT, hot, cold, middle, dTdx, dTdy, tm, um, vm, cellArea, sectionWeight
    real(kind=8) :: w(4), dw(4)
    real(kind=8) :: tmin, tmax, rmin, rmax

    do j=yStartIntegration,yEndIntegration
        do i=xStartIntegration,xEndIntegration
            call owned_cell_area(xOffset+(i-0.5d0)*dx,yOffset+(j-0.5d0)*dx,dx,integrationBox,cellArea)
            if(cellArea <=0d0) cycle
            if(.NOT. ieee_is_finite(T(i,j)) .OR. .NOT. ieee_is_finite(rho(i,j)) .OR. &
                .NOT. ieee_is_finite(u(i,j)) .OR. .NOT. ieee_is_finite(v(i,j)) .OR. rho(i,j) <=0.0d0) then
                write(*,*) 'Invalid state: time, xOffset,yOffset, i,j:',reloadItcOffset+itc,xOffset,yOffset,i,j
                error stop 'Nonfinite or nonpositive density in owned cells'
            endif
#ifdef SideHeatedCell
            conv=conv+u(i,j)*T(i,j)*cellArea
#else
            conv=conv+v(i,j)*T(i,j)*cellArea
#endif
            vel2=vel2+(u(i,j)**2+v(i,j)**2)*cellArea
            mass=mass+rho(i,j)*cellArea
            meanT=meanT+T(i,j)*cellArea
            tmin=min(tmin,T(i,j))
            tmax=max(tmax,T(i,j))
            rmin=min(rmin,rho(i,j))
            rmax=max(rmax,rho(i,j))
        enddo
    enddo

#ifdef SideHeatedCell
    if(integrationBox(1)==0.0d0) then
        do j=yStartIntegration,yEndIntegration
            hot=hot+(8.0d0*Thot-9.0d0*T(1,j)+T(2,j))/(3.0d0*dx)*dx/dble(ny)
        enddo
    endif

    if(integrationBox(2)==dble(nx)) then
        do j=yStartIntegration,yEndIntegration
            cold=cold+(-8.0d0*Tcold+9.0d0*T(ni,j)-T(ni-1,&
                j))/(3.0d0*dx)*dx/dble(ny)
        enddo
    endif

    if((0.5d0*nx) >=integrationBox(1) .AND. (0.5d0*nx) < integrationBox(2)) then

        call section_weights(((0.5d0*nx)-xOffset)/dx+0.5d0,ni,im,w,dw)
        do j=yStartIntegration,yEndIntegration
            call section_owned_weight(xOffset,yOffset,dx,integrationBox,j,1,(0.5d0*nx),sectionWeight)
            if(sectionWeight <=0d0) cycle
            if(dx==1.0d0 .AND. refineRatio > 1) then

                call fine_section_x((0.5d0*nx),nint(yOffset)+j,tm,um,dTdx)
            else
                tm=sum(w*T(im:im+3,j))
                um=sum(w*u(im:im+3,j))
                dTdx=sum(dw*T(im:im+3,j))/dx
            endif

            middle=middle+(um*tm/diffusivity-dTdx)*sectionWeight/dble(ny)
        enddo
    endif
#else

    if(integrationBox(3)==0.0d0) then
        do i=xStartIntegration,xEndIntegration
            hot=hot+(8.0d0*Thot-9.0d0*T(i,1)+T(i,2))/(3.0d0*dx)*dx/dble(nx)
        enddo
    endif

    if(integrationBox(4)==dble(ny)) then
        do i=xStartIntegration,xEndIntegration
            cold=cold+(-8.0d0*Tcold+9.0d0*T(i,nj)-T(i,&
                nj-1))/(3.0d0*dx)*dx/dble(nx)
        enddo
    endif

    if((0.5d0*ny) >=integrationBox(3) .AND. (0.5d0*ny) < integrationBox(4)) then

        call section_weights(((0.5d0*ny)-yOffset)/dx+0.5d0,nj,jm,w,dw)
        do i=xStartIntegration,xEndIntegration
            call section_owned_weight(xOffset,yOffset,dx,integrationBox,i,2,(0.5d0*ny),sectionWeight)
            if(sectionWeight <=0d0) cycle
            if(dx==1.0d0 .AND. refineRatio > 1) then

                call fine_section_y(nint(xOffset)+i,(0.5d0*ny),tm,vm,dTdy)
            else
                tm=sum(w*T(i,jm:jm+3))
                vm=sum(w*v(i,jm:jm+3))
                dTdy=sum(dw*T(i,jm:jm+3))/dx
            endif

            middle=middle+(vm*tm/diffusivity-dTdy)*sectionWeight/dble(nx)
        enddo
    endif
#endif
  end subroutine calNuRe_grid
!===================================================================================================

!===================================================================================================
! 子程序: check
! 作用: 稳态模式下汇总速度和温度误差，并输出收敛记录。
!===================================================================================================
  subroutine check()
#ifdef steadyFlow

    use commondata
    implicit none

    integer(kind=4) :: k
    real(kind=8) :: du, uu, dt, tt

    du=0.0d0
    uu=0.0d0
    dt=0.0d0
    tt=0.0d0

    call update_host_all(.false.)

    call check_grid(size(f_coarse,1),size(f_coarse,2),dxCoarse,xOffsetCoarse,yOffsetCoarse,(1+merge(interfaceSkin,0,refineRatio>1)),&
        (size(f_coarse,1)-merge(interfaceSkin,0,refineRatio>1)),(1+merge(interfaceSkin,0,refineRatio>1)),(size(f_coarse,2)-merge(interfaceSkin,0,refineRatio>1)),centerIntegrationBox,du,uu,dt,tt,&
        u_coarse,v_coarse,T_coarse &
#ifdef steadyFlow
        ,up_coarse,vp_coarse,Tp_coarse &
#endif
        )
    if(refineRatio > 1) then

        call check_grid(size(f_left,1),size(f_left,2),dxFine,xOffsetLeft,yOffsetLeft,1,size(f_left,1),&
            1,size(f_left,2),[xOffsetLeft,xOffsetLeft+size(f_left,1),yOffsetLeft,yOffsetLeft+size(f_left,2)],du,uu,dt,tt,u_left,v_left,T_left &
#ifdef steadyFlow
            ,up_left,vp_left,Tp_left &
#endif
            )

        call check_grid(size(f_right,1),size(f_right,2),dxFine,xOffsetRight,yOffsetRight,1,size(f_right,1),&
            1,size(f_right,2),[xOffsetRight,xOffsetRight+size(f_right,1),yOffsetRight,yOffsetRight+size(f_right,2)],du,uu,dt,tt,&
            u_right,v_right,T_right &
#ifdef steadyFlow
            ,up_right,vp_right,Tp_right &
#endif
            )

        call check_grid(size(f_bottom,1),size(f_bottom,2),dxFine,xOffsetBottom,yOffsetBottom,1,&
            size(f_bottom,1),1,size(f_bottom,2),[xOffsetBottom,xOffsetBottom+size(f_bottom,1),yOffsetBottom,yOffsetBottom+size(f_bottom,2)],du,uu,&
            dt,tt,u_bottom,v_bottom,&
            T_bottom &
#ifdef steadyFlow
            ,up_bottom,vp_bottom,Tp_bottom &
#endif
            )

        call check_grid(size(f_top,1),size(f_top,2),dxFine,xOffsetTop,yOffsetTop,1,size(f_top,1),&
            1,&
            size(f_top,2),[xOffsetTop,xOffsetTop+size(f_top,1),yOffsetTop,yOffsetTop+size(f_top,2)],du,uu,dt,tt,u_top,v_top,T_top &
#ifdef steadyFlow
            ,up_top,vp_top,Tp_top &
#endif
            )
    endif

    errorU=sqrt(du/max(uu,1.0d-300))
    errorT=sqrt(dt/max(tt,1.0d-300))
    open(newunit=k,file='Convergence_2DOpenaccMultiblock.dat',status='unknown',position='append')
    write(k,'(I12,2(1X,ES24.16E3))') reloadItcOffset+itc,errorU,errorT
    close(k)
    write(*,*) 'errorU,errorT:',errorU,errorT
#endif
  end subroutine check
!===================================================================================================

!===================================================================================================

#ifdef steadyFlow

!===================================================================================================
! 子程序: check_grid
! 作用: 累计一个区域的收敛误差，并保存本次检查时的宏观场。
!===================================================================================================
  subroutine check_grid(ni,nj,dx,xOffset,yOffset,xStartIntegration,xEndIntegration,yStartIntegration,&
      yEndIntegration,integrationBox,du,uu,dt,&
      tt,u,v,T &
#ifdef steadyFlow
      ,up,vp,Tp &
#endif
      )

    implicit none

    ! up/vp/Tp 是上一次收敛检查的场，不是上一时间步，也不是接口插值历史。
    ! 先按节点实际面积累计误差，再把当前 u/v/T 保存为下一次检查的参考。

    integer(kind=4) :: ni, nj, xStartIntegration, xEndIntegration, yStartIntegration, yEndIntegration
    real(kind=8) :: dx, xOffset, yOffset, integrationBox(4)

    integer(kind=4) :: i, j
    real(kind=8) :: u(ni,nj)
    real(kind=8) :: v(ni,nj)
    real(kind=8) :: T(ni,nj)
#ifdef steadyFlow
    real(kind=8) :: up(ni,nj)
    real(kind=8) :: vp(ni,nj)
    real(kind=8) :: Tp(ni,nj)
#endif
    real(kind=8) :: du, uu, dt, tt, cellArea

    do j=yStartIntegration,yEndIntegration
        do i=xStartIntegration,xEndIntegration
            call owned_cell_area(xOffset+(i-0.5d0)*dx,yOffset+(j-0.5d0)*dx,dx,integrationBox,cellArea)
            du=du+cellArea*((u(i,j)-up(i,j))**2+(v(i,j)-vp(i,j))**2)
            uu=uu+cellArea*(u(i,j)**2+v(i,j)**2)
            dt=dt+cellArea*(T(i,j)-Tp(i,j))**2
            tt=tt+cellArea*T(i,j)**2
        enddo
    enddo

    up=u
    vp=v
    Tp=T
  end subroutine check_grid
!===================================================================================================

!===================================================================================================

#endif

!===================================================================================================
! 子程序: output_time_string
! 作用: 将实际累计 t_ff 写成文件名中的时间，去掉小数末尾多余的零。
!===================================================================================================
  subroutine output_time_string(timeString)

    use commondata
    implicit none

    character(len=40), intent(out) :: timeString
    integer(kind=4) :: i

    ! 保留十二位小数；最小输出间隔为一个粗步，检查精度以避免不同时刻重名。
    if(dble(refineRatio)/timeUnit<=1.0d-12) error stop 'Increase tff filename precision for this timeUnit'
    write(timeString,'(F40.12)') dble(reloadItcOffset+itc)/timeUnit
    if(index(timeString,'*')/=0) error stop 'tff exceeds filename field width'
    timeString=adjustl(timeString)
    i=len_trim(timeString)
    do while(timeString(i:i)=='0')
        timeString(i:i)=' '
        i=i-1
    enddo
    if(timeString(i:i)=='.') timeString(i:i)=' '
  end subroutine output_time_string
!===================================================================================================

!===================================================================================================
! 子程序: output_Tecplot
! 作用: 输出各具名区域的 Tecplot 数据及节点积分面积。
!===================================================================================================
  subroutine output_Tecplot()

    use commondata
    implicit none

    integer(kind=4) :: k
    character(len=40) :: timeString

    call output_time_string(timeString)

    call update_host_all(.false.)
    open(newunit=k,file=pltFolderPrefix//'-'//trim(timeString)//'tff.dat',status='replace')
    write(k,'(a)') 'VARIABLES="x/L","y/L","u","v","T","rho","dx/L","integration_area/L^2"'

    call output_Tecplot_grid(size(f_coarse,1),size(f_coarse,2),dxCoarse,xOffsetCoarse,yOffsetCoarse,(1+merge(interfaceSkin,0,refineRatio>1)),&
        (size(f_coarse,1)-merge(interfaceSkin,0,refineRatio>1)),(1+merge(interfaceSkin,0,refineRatio>1)),(size(f_coarse,2)-merge(interfaceSkin,0,refineRatio>1)),centerIntegrationBox,&
        'coarse',k,rho_coarse,u_coarse,&
        v_coarse,T_coarse)
    if(refineRatio > 1) then

        call output_Tecplot_grid(size(f_left,1),size(f_left,2),dxFine,xOffsetLeft,yOffsetLeft,1,size(f_left,1),&
            1,size(f_left,2),[xOffsetLeft,xOffsetLeft+size(f_left,1),yOffsetLeft,yOffsetLeft+size(f_left,2)],'left',k,rho_left,u_left,v_left,T_left)

        call output_Tecplot_grid(size(f_right,1),size(f_right,2),dxFine,xOffsetRight,yOffsetRight,1,&
            size(f_right,1),1,size(f_right,2),[xOffsetRight,xOffsetRight+size(f_right,1),yOffsetRight,yOffsetRight+size(f_right,2)],'right',k,&
            rho_right,u_right,&
            v_right,T_right)

        call output_Tecplot_grid(size(f_bottom,1),size(f_bottom,2),dxFine,xOffsetBottom,yOffsetBottom,1,&
            size(f_bottom,1),1,size(f_bottom,2),[xOffsetBottom,xOffsetBottom+size(f_bottom,1),yOffsetBottom,yOffsetBottom+size(f_bottom,2)],'bottom',&
            k,rho_bottom,u_bottom,&
            v_bottom,T_bottom)

        call output_Tecplot_grid(size(f_top,1),size(f_top,2),dxFine,xOffsetTop,yOffsetTop,1,size(f_top,1),&
            1,size(f_top,2),[xOffsetTop,xOffsetTop+size(f_top,1),yOffsetTop,yOffsetTop+size(f_top,2)],'top',k,rho_top,u_top,v_top,T_top)
    endif

    close(k)
  end subroutine output_Tecplot
!===================================================================================================

!===================================================================================================
! 子程序: output_Tecplot_grid
! 作用: 向已打开的 Tecplot 文件写入一个区域的坐标及宏观场。
!===================================================================================================
  subroutine output_Tecplot_grid(ni,nj,dx,xOffset,yOffset,xStartIntegration,xEndIntegration,&
      yStartIntegration,yEndIntegration,integrationBox,&
      gridName,k,rho,u,v,T)

    use commondata,only: lengthUnit,timeUnit,itc,reloadItcOffset
    implicit none
    integer(kind=4) :: ni, nj, xStartIntegration, xEndIntegration, yStartIntegration, yEndIntegration
    real(kind=8) :: dx, xOffset, yOffset, integrationBox(4)
    character(*), intent(in) :: gridName

    integer(kind=4) :: k, i, j
    real(kind=8) :: rho(ni,nj)
    real(kind=8) :: u(ni,nj)
    real(kind=8) :: v(ni,nj)
    real(kind=8) :: T(ni,nj)
    real(kind=8) :: cellArea

    write(k,'(a,a,a,I0,a,I0,a,ES24.16E3)') 'ZONE T="',trim(gridName),'", I=',xEndIntegration-xStartIntegration+1,&
        ', J=',yEndIntegration-yStartIntegration+1,', F=POINT, SOLUTIONTIME=',dble(reloadItcOffset+itc)/timeUnit
    do j=yStartIntegration,yEndIntegration
        do i=xStartIntegration,xEndIntegration
            call owned_cell_area(xOffset+(i-0.5d0)*dx,yOffset+(j-0.5d0)*dx,dx,integrationBox,cellArea)
            write(k,'(8(ES24.16E3,1X))')((xOffset+(dble(i)-0.5d0)*dx))/lengthUnit,&
                ((yOffset+(dble(j)-0.5d0)*dx))/lengthUnit,u(i,j),v(i,j),T(i,j),rho(i,j),&
                dx/lengthUnit,&
                cellArea/lengthUnit**2
        enddo
    enddo
  end subroutine output_Tecplot_grid
!===================================================================================================

!===================================================================================================
! 子程序: output_SnapshotFile
! 作用: 输出包含几何、积分权重及宏观场的快照文件。
!===================================================================================================
  subroutine output_SnapshotFile()

    use commondata
    implicit none

    ! 快照记录区域数量、各区几何和节点积分面积；单网格为一区，多块网格为五区。
    ! 快照用于后处理；完整续算还需要 f/g 和粗历史，应读取续算检查点。

    integer(kind=4) :: k
    character(len=40) :: timeString

    call output_time_string(timeString)

    call update_host_all(.false.)
    open(newunit=k,file=snapshotFilePrefix//'-'//trim(timeString)//'tff.bin',form='unformatted',&
        access='stream',status='replace')
    write(k) merge(5,1,refineRatio > 1),nx,ny,reloadItcOffset+itc,dble(reloadItcOffset+itc)/timeUnit,lengthUnit

    call output_SnapshotFile_grid(size(f_coarse,1),size(f_coarse,2),dxCoarse,xOffsetCoarse,yOffsetCoarse,&
        (1+merge(interfaceSkin,0,refineRatio>1)),(size(f_coarse,1)-merge(interfaceSkin,0,refineRatio>1)),(1+merge(interfaceSkin,0,refineRatio>1)),(size(f_coarse,2)-merge(interfaceSkin,0,refineRatio>1)),&
        centerIntegrationBox,k,rho_coarse,u_coarse,&
        v_coarse,T_coarse)
    if(refineRatio > 1) then

        call output_SnapshotFile_grid(size(f_left,1),size(f_left,2),dxFine,xOffsetLeft,yOffsetLeft,1,&
            size(f_left,1),1,size(f_left,2),[xOffsetLeft,xOffsetLeft+size(f_left,1),yOffsetLeft,yOffsetLeft+size(f_left,2)],k,rho_left,&
            u_left,v_left,T_left)

        call output_SnapshotFile_grid(size(f_right,1),size(f_right,2),dxFine,xOffsetRight,yOffsetRight,1,&
            size(f_right,1),1,size(f_right,2),[xOffsetRight,xOffsetRight+size(f_right,1),yOffsetRight,yOffsetRight+size(f_right,2)],k,rho_right,&
            u_right,v_right,T_right)

        call output_SnapshotFile_grid(size(f_bottom,1),size(f_bottom,2),dxFine,xOffsetBottom,yOffsetBottom,&
            1,size(f_bottom,1),1,size(f_bottom,2),&
            [xOffsetBottom,xOffsetBottom+size(f_bottom,1),yOffsetBottom,yOffsetBottom+size(f_bottom,2)],k,rho_bottom,u_bottom,&
            v_bottom,T_bottom)

        call output_SnapshotFile_grid(size(f_top,1),size(f_top,2),dxFine,xOffsetTop,yOffsetTop,1,size(f_top,1),&
            1,size(f_top,2),[xOffsetTop,xOffsetTop+size(f_top,1),yOffsetTop,yOffsetTop+size(f_top,2)],k,rho_top,u_top,v_top,T_top)
    endif

    close(k)
  end subroutine output_SnapshotFile
!===================================================================================================

!===================================================================================================
! 子程序: output_SnapshotFile_grid
! 作用: 向快照文件写入一个区域的几何、权重和宏观量。
!===================================================================================================
  subroutine output_SnapshotFile_grid(ni,nj,dx,xOffset,yOffset,xStartIntegration,xEndIntegration,&
      yStartIntegration,yEndIntegration,&
      integrationBox,k,rho,u,v,T)

    implicit none
    integer(kind=4) :: ni, nj, xStartIntegration, xEndIntegration, yStartIntegration, yEndIntegration
    real(kind=8) :: dx, xOffset, yOffset, integrationBox(4)

    integer(kind=4) :: k, i, j
    real(kind=8) :: rho(ni,nj)
    real(kind=8) :: u(ni,nj)
    real(kind=8) :: v(ni,nj)
    real(kind=8) :: T(ni,nj)
    real(kind=8) :: cellAreaRow(xStartIntegration:xEndIntegration)    ! 暂存一行节点的积分面积，按原 i 快、j 慢顺序写入流文件。

    write(k) xEndIntegration-xStartIntegration+1,yEndIntegration-yStartIntegration+1,&
        xOffset+(dble(xStartIntegration)-0.5d0)*dx,yOffset+(dble(yStartIntegration)-0.5d0)*dx,&
        dx,integrationBox
    ! 保持已有快照格式；一维宽度仅在写文件时计算，不在 commondata 重复保存。
    write(k) (max(0.0d0,min(integrationBox(2),xOffset+i*dx)-&
        max(integrationBox(1),xOffset+(i-1)*dx)),i=xStartIntegration,xEndIntegration),&
        (max(0.0d0,min(integrationBox(4),yOffset+j*dx)-&
        max(integrationBox(3),yOffset+(j-1)*dx)),j=yStartIntegration,yEndIntegration)
    ! 一维积分权重之后保存显式的二维积分面积。
    do j=yStartIntegration,yEndIntegration
        do i=xStartIntegration,xEndIntegration
            call owned_cell_area(xOffset+(i-0.5d0)*dx,yOffset+(j-0.5d0)*dx,dx,&
                integrationBox,cellAreaRow(i))
        enddo
        write(k) cellAreaRow
    enddo
    write(k) u(xStartIntegration:xEndIntegration,yStartIntegration:yEndIntegration),&
        v(xStartIntegration:xEndIntegration,yStartIntegration:yEndIntegration),&
        T(xStartIntegration:xEndIntegration,yStartIntegration:yEndIntegration),&
        rho(xStartIntegration:xEndIntegration,yStartIntegration:yEndIntegration)
  end subroutine output_SnapshotFile_grid
!===================================================================================================

!===================================================================================================
! 子程序: output_ReloadFile
! 作用: 保存完整检查点，并在写入完成后更新 latest.meta。
!===================================================================================================
  subroutine output_ReloadFile()

    use commondata
    implicit none

    ! 按 coarse、left、right、bottom、top 的固定顺序保存各区数据。
    ! 粗历史和输出时钟也必须写入，否则恢复后时间插值及采样次序无法保持连续。

    integer(kind=4) :: k, currentModel(12)
    real(kind=8) :: currentPhysics(16)
    character(len=40) :: timeString
    character(len=256) :: name

    call update_host_all(.true.)
    call output_time_string(timeString)
    name=reloadFilePrefix//'-'//trim(timeString)//'tff.bin'
    open(newunit=k,file=trim(name),access='stream',form='unformatted',status='replace')
    write(k) nx,ny,refineRatio,fineLayerCellsLeft,fineLayerCellsRight,&
        fineLayerCellsBottom,fineLayerCellsTop,interfaceSkin,merge(5,1,refineRatio > 1)

    ! 保存模型开关和物理参数，使续算能够检查配置是否匹配。
    call model_signature(currentModel)

    call physical_signature(currentPhysics)
    write(k) currentModel,currentPhysics
    ! 三个输出序号依次为采样、续算文件和 Tecplot 的下一序号，用于续算时核对输出间隔。
    ! 旧文件的三个命名计数位置保留并写零，保持检查点布局兼容。
    write(k) reloadItcOffset+itc,NuReSampleCount+1,&
        (reloadItcOffset+itc)/reloadFileIntervalItc+1,&
        (reloadItcOffset+itc)/outputPltFileIntervalItc+1,&
        0,0,0,errorU,errorT

    call output_ReloadFile_grid(size(f_coarse,1),size(f_coarse,2),timeInterpMaxIndexCoarse,dxCoarse,xOffsetCoarse,&
        yOffsetCoarse,(1+merge(interfaceSkin,0,refineRatio>1)),(size(f_coarse,1)-merge(interfaceSkin,0,refineRatio>1)),(1+merge(interfaceSkin,0,refineRatio>1)),&
        (size(f_coarse,2)-merge(interfaceSkin,0,refineRatio>1)),centerIntegrationBox,k,&
        f_coarse,g_coarse,u_coarse,v_coarse,T_coarse,rho_coarse,Fx_coarse,Fy_coarse &
#ifdef steadyFlow
        ,up_coarse,vp_coarse,Tp_coarse &
#endif
        )
    if(refineRatio > 1) then

        call output_ReloadFile_grid(size(f_left,1),size(f_left,2),0,dxFine,xOffsetLeft,yOffsetLeft,1,&
            size(f_left,1),1,size(f_left,2),[xOffsetLeft,xOffsetLeft+size(f_left,1),yOffsetLeft,yOffsetLeft+size(f_left,2)],k,f_left,g_left,&
            u_left,v_left,T_left,&
            rho_left,Fx_left,Fy_left &
#ifdef steadyFlow
            ,up_left,vp_left,Tp_left &
#endif
            )

        call output_ReloadFile_grid(size(f_right,1),size(f_right,2),0,dxFine,xOffsetRight,yOffsetRight,1,&
            size(f_right,1),1,size(f_right,2),[xOffsetRight,xOffsetRight+size(f_right,1),yOffsetRight,yOffsetRight+size(f_right,2)],k,f_right,&
            g_right,u_right,v_right,&
            T_right,rho_right,Fx_right,Fy_right &
#ifdef steadyFlow
            ,up_right,vp_right,Tp_right &
#endif
            )

        call output_ReloadFile_grid(size(f_bottom,1),size(f_bottom,2),0,dxFine,xOffsetBottom,yOffsetBottom,&
            1,size(f_bottom,1),1,size(f_bottom,2),&
            [xOffsetBottom,xOffsetBottom+size(f_bottom,1),yOffsetBottom,yOffsetBottom+size(f_bottom,2)],k,f_bottom,g_bottom,&
            u_bottom,v_bottom,T_bottom,rho_bottom,Fx_bottom,Fy_bottom &
#ifdef steadyFlow
            ,up_bottom,vp_bottom,Tp_bottom &
#endif
            )

        call output_ReloadFile_grid(size(f_top,1),size(f_top,2),0,dxFine,xOffsetTop,yOffsetTop,1,size(f_top,1),&
            1,size(f_top,2),[xOffsetTop,xOffsetTop+size(f_top,1),yOffsetTop,yOffsetTop+size(f_top,2)],k,f_top,g_top,u_top,v_top,T_top,&
            rho_top,Fx_top,&
            Fy_top &
#ifdef steadyFlow
            ,up_top,vp_top,Tp_top &
#endif
            )
    endif

    close(k)
    ! 完整写完当前时刻的 checkpoint 后再更新 latest 指针，其他时刻的 checkpoint 仍保留。
    open(newunit=k,file=reloadFilePrefix//'-latest.meta',status='replace')
    write(k,'(a)') trim(name)
    close(k)
  end subroutine output_ReloadFile
!===================================================================================================

!===================================================================================================
! 子程序: output_ReloadFile_grid
! 作用: 向检查点写入一个区域的几何、分布函数、宏观量和历史。
!===================================================================================================
  subroutine output_ReloadFile_grid(ni,nj,nh,dx,xOffset,yOffset,xStartIntegration,xEndIntegration,&
      yStartIntegration,yEndIntegration,&
      integrationBox,k,f,g,u,v,T,rho,Fx,Fy &
#ifdef steadyFlow
      ,up,vp,Tp &
#endif
      )

    use commondata,only: coarseBuffer,dxCoarse
    implicit none

    ! 每个时间层按宏观量、归一化力、流场矩、温度矩的固定顺序写入。
    ! 读写顺序与 read_restart_grid 完全对应；稳态模式另存 up/vp/Tp。

    integer(kind=4) :: ni, nj, nh, xStartIntegration, xEndIntegration, yStartIntegration, yEndIntegration
    real(kind=8) :: dx, xOffset, yOffset, integrationBox(4)

    integer(kind=4) :: a
    real(kind=8) :: unusedHistory(ni,nj)

    integer(kind=4) :: k
    real(kind=8) :: f(ni,nj,0:8)
    real(kind=8) :: g(ni,nj,0:4)
    real(kind=8) :: u(ni,nj)
    real(kind=8) :: v(ni,nj)
    real(kind=8) :: T(ni,nj)
    real(kind=8) :: rho(ni,nj)
    real(kind=8) :: Fx(ni,nj)
    real(kind=8) :: Fy(ni,nj)

#ifdef steadyFlow
    real(kind=8) :: up(ni,nj)
    real(kind=8) :: vp(ni,nj)
    real(kind=8) :: Tp(ni,nj)
#endif

    write(k) ni,nj,xStartIntegration,xEndIntegration,yStartIntegration,yEndIntegration,nh,&
        xOffset,yOffset,dx,integrationBox
    ! 每个时间层按宏观量、归一化力、流场矩、温度矩的顺序写入。
    write(k) f,g,u,v,T,rho,Fx,Fy
    if(dx==dxCoarse) then
        write(k) coarseBuffer
    else
        ! 仅在主机写文件时复用一个平面，保留旧细历史在检查点中的位置。
        unusedHistory=0.0d0
        do a=0,19
            write(k) unusedHistory
        enddo
    endif
#ifdef steadyFlow
    write(k) up,vp,Tp
#endif
  end subroutine output_ReloadFile_grid
!===================================================================================================

!===================================================================================================
! 子程序: read_restart
! 作用: 检查网格及物理配置，恢复时间计数和五套网格状态。
!===================================================================================================
  subroutine read_restart()

    use commondata
    implicit none

    ! 文件直接从网格参数开始；核对网格、物理参数及输出间隔后才恢复场量。
    ! 旧 v10 连通细环文件需先用 convert_restart_v10.py 转换，不能按五套数组直接读取。
    ! interfaceSkin 决定各细数组尺寸和偏移，由 read_restart_grid 核对；旧对称布局不可直接续算。
    ! 历史对称布局需使用匹配该布局的源码，不能只修改当前的接收层数。

    integer(kind=4) :: k, ios, head(9), sig(12), currentModel(12)
    integer(kind=4) :: savedOutputIndex(3), unusedFileNumbers(3)
    real(kind=8) :: phys(16), currentPhysics(16)
    character(len=256) :: name
    logical :: metaExists

    inquire(file=reloadFilePrefix//'-latest.meta',exist=metaExists)

    ! 优先从 latest.meta 找到检查点；meta 不存在时使用手动填写的文件名。
    if(metaExists) then
        open(newunit=k,file=reloadFilePrefix//'-latest.meta',status='old',iostat=ios)
        if(ios/=0) error stop 'Cannot open multiblock latest.meta'
        read(k,'(a)',iostat=ios) name
        close(k)
        if(ios/=0 .OR. len_trim(name)==0) error stop 'Invalid latest.meta'
    else
        if(len_trim(reloadFileName)==0) error stop 'Missing latest.meta: set reloadFileName to a saved checkpoint filename'
        name=reloadFileName
    endif

    open(newunit=k,file=trim(name),status='old',access='stream',form='unformatted',iostat=ios)
    if(ios/=0) error stop 'Missing multiblock checkpoint'
    read(k,iostat=ios) head
    if(ios/=0) error stop 'Truncated multiblock checkpoint header'
    if(any(head/=[nx,ny,refineRatio,fineLayerCellsLeft,fineLayerCellsRight,&
        fineLayerCellsBottom,fineLayerCellsTop,interfaceSkin,merge(5,1,refineRatio > 1)])) &
        error stop 'Restart mesh/refinement mismatch'
    read(k) sig,phys

    call model_signature(currentModel)

    call physical_signature(currentPhysics)
    if(any(sig/=currentModel) .OR. any(phys/=currentPhysics)) &
        error stop 'Restart model/physics/output-cadence mismatch'
    ! 文件中的三个序号依次对应采样、续算文件和 Tecplot；仅用于恢复采样次数及检查间隔一致性。
    ! 兼容读取旧命名计数；输出文件名只由累计步数换算的 t_ff 决定。
    read(k) reloadItcOffset,savedOutputIndex,unusedFileNumbers,errorU,errorT
    if(reloadItcOffset < 0 .OR. mod(reloadItcOffset,refineRatio)/=0) error stop 'Restart is not at a synchronized time'
    if(any(savedOutputIndex/=[reloadItcOffset/outputSnapshotIntervalItc+1,&
        reloadItcOffset/reloadFileIntervalItc+1,reloadItcOffset/outputPltFileIntervalItc+1])) &
        error stop 'Restart output schedule differs from the fixed step intervals'
    itc=0
    NuReSampleCount=savedOutputIndex(1)-1

    call read_restart_grid(size(f_coarse,1),size(f_coarse,2),timeInterpMaxIndexCoarse,dxCoarse,xOffsetCoarse,&
        yOffsetCoarse,(1+merge(interfaceSkin,0,refineRatio>1)),(size(f_coarse,1)-merge(interfaceSkin,0,refineRatio>1)),(1+merge(interfaceSkin,0,refineRatio>1)),&
        (size(f_coarse,2)-merge(interfaceSkin,0,refineRatio>1)),centerIntegrationBox,k,&
        f_coarse,g_coarse,u_coarse,v_coarse,T_coarse,rho_coarse,Fx_coarse,Fy_coarse &
#ifdef steadyFlow
        ,up_coarse,vp_coarse,Tp_coarse &
#endif
        )
    if(refineRatio > 1) then

        call read_restart_grid(size(f_left,1),size(f_left,2),0,dxFine,xOffsetLeft,yOffsetLeft,1,&
            size(f_left,1),&
            1,size(f_left,2),[xOffsetLeft,xOffsetLeft+size(f_left,1),yOffsetLeft,yOffsetLeft+size(f_left,2)],k,f_left,g_left,u_left,v_left,&
            T_left,rho_left,&
            Fx_left,Fy_left &
#ifdef steadyFlow
            ,up_left,vp_left,Tp_left &
#endif
            )

        call read_restart_grid(size(f_right,1),size(f_right,2),0,dxFine,xOffsetRight,yOffsetRight,1,&
            size(f_right,1),1,size(f_right,2),[xOffsetRight,xOffsetRight+size(f_right,1),yOffsetRight,yOffsetRight+size(f_right,2)],k,f_right,&
            g_right,u_right,v_right,&
            T_right,rho_right,Fx_right,Fy_right &
#ifdef steadyFlow
            ,up_right,vp_right,Tp_right &
#endif
            )

        call read_restart_grid(size(f_bottom,1),size(f_bottom,2),0,dxFine,xOffsetBottom,yOffsetBottom,1,&
            size(f_bottom,1),1,size(f_bottom,2),[xOffsetBottom,xOffsetBottom+size(f_bottom,1),yOffsetBottom,yOffsetBottom+size(f_bottom,2)],k,&
            f_bottom,g_bottom,u_bottom,&
            v_bottom,T_bottom,rho_bottom,Fx_bottom,Fy_bottom &
#ifdef steadyFlow
            ,up_bottom,vp_bottom,Tp_bottom &
#endif
            )

        call read_restart_grid(size(f_top,1),size(f_top,2),0,dxFine,xOffsetTop,yOffsetTop,1,size(f_top,1),&
            1,size(f_top,2),[xOffsetTop,xOffsetTop+size(f_top,1),yOffsetTop,yOffsetTop+size(f_top,2)],k,f_top,g_top,u_top,v_top,T_top,&
            rho_top,Fx_top,&
            Fy_top &
#ifdef steadyFlow
            ,up_top,vp_top,Tp_top &
#endif
            )
    endif

    close(k)
    reloadDimensionlessTime=dble(reloadItcOffset)/timeUnit
    reloadMetadataLoaded=.true.
  end subroutine read_restart
!===================================================================================================

!===================================================================================================
! 子程序: read_restart_grid
! 作用: 核对一个区域的尺寸和坐标，并读取其完整状态与历史。
!===================================================================================================
  subroutine read_restart_grid(ni,nj,nh,dx,xOffset,yOffset,xStartIntegration,xEndIntegration,&
      yStartIntegration,yEndIntegration,integrationBox,&
      k,f,g,u,v,T,rho,Fx,Fy &
#ifdef steadyFlow
      ,up,vp,Tp &
#endif
      )

    use commondata,only: coarseBuffer,dxCoarse
    implicit none
    integer(kind=4) :: ni, nj, nh, xStartIntegration, xEndIntegration, yStartIntegration, yEndIntegration
    real(kind=8) :: dx, xOffset, yOffset, integrationBox(4)

    integer(kind=4) :: a
    real(kind=8) :: unusedHistory(ni,nj)

    integer(kind=4) :: k, ios, geom(7)
    real(kind=8) :: f(ni,nj,0:8)
    real(kind=8) :: g(ni,nj,0:4)
    real(kind=8) :: u(ni,nj)
    real(kind=8) :: v(ni,nj)
    real(kind=8) :: T(ni,nj)
    real(kind=8) :: rho(ni,nj)
    real(kind=8) :: Fx(ni,nj)
    real(kind=8) :: Fy(ni,nj)

#ifdef steadyFlow
    real(kind=8) :: up(ni,nj)
    real(kind=8) :: vp(ni,nj)
    real(kind=8) :: Tp(ni,nj)
#endif
    real(kind=8) :: coord(7)

    read(k) geom,coord
    if(any(geom/=[ni,nj,xStartIntegration,xEndIntegration,yStartIntegration,yEndIntegration,nh]) .OR. &
        any(coord/=[xOffset,yOffset,dx,integrationBox])) error stop 'Restart block layout mismatch'
    read(k,iostat=ios) f,g,u,v,T,rho,Fx,Fy
    if(ios/=0) error stop 'Incomplete checkpoint state'
    if(dx==dxCoarse) then
        read(k,iostat=ios) coarseBuffer
        if(ios/=0) error stop 'Incomplete coarse history'
    else
        ! 读完旧细历史占用的位置；推进直接使用细网格当前场。
        do a=0,19
            read(k,iostat=ios) unusedHistory
            if(ios/=0) error stop 'Incomplete legacy fine history'
        enddo
    endif
#ifdef steadyFlow
    read(k) up,vp,Tp
#endif
  end subroutine read_restart_grid
!===================================================================================================

!===================================================================================================
! 子程序: model_signature
! 作用: 将当前物理模型和边界宏配置写入检查点配置标识。
!===================================================================================================
  subroutine model_signature(sig)

    implicit none

    integer(kind=4), intent(out) :: sig(12)

    sig=0
#ifdef steadyFlow
    sig(1)=1
#endif
#ifdef SideHeatedCell
    sig(3)=1
#endif
#ifdef HorizontalWallsNoslip
    sig(4)=1
#endif
#ifdef VerticalWallsNoslip
    sig(5)=1
#endif
#ifdef VerticalWallsPeriodicalU
    sig(6)=1
#endif
#ifdef HorizontalWallsConstT
    sig(7)=1
#endif
#ifdef HorizontalWallsAdiabatic
    sig(8)=1
#endif
#ifdef VerticalWallsConstT
    sig(9)=1
#endif
#ifdef VerticalWallsAdiabatic
    sig(10)=1
#endif
#ifdef VerticalWallsPeriodicalT
    sig(11)=1
#endif
#ifdef SideHeatedHa
    sig(12)=1
#endif
  end subroutine model_signature
!===================================================================================================

!===================================================================================================
! 子程序: physical_signature
! 作用: 记录物理参数和输出间隔，供续算时核对兼容性。
!===================================================================================================
  subroutine physical_signature(sig)

    use commondata,only: Rayleigh,Prandtl,Mach,Thot,Tcold,Snu,Sq,Qk,Qnu,paraA,&
        outputSnapshotInterval,reloadFileInterval,outputPltFileInterval
#ifdef SideHeatedHa
    use commondata,only: Ha,phi
#endif
    implicit none

    real(kind=8), intent(out) :: sig(16)

    sig=[Rayleigh,Prandtl,Mach,Thot,Tcold,Snu,Sq,Qk,Qnu,paraA,outputSnapshotInterval,&
        reloadFileInterval,outputPltFileInterval,0.0d0,0.0d0,0.0d0]
#ifdef SideHeatedHa
    sig(14)=Ha
    sig(15)=phi
#endif
  end subroutine physical_signature
!===================================================================================================

!===================================================================================================
! 子程序: check_history
! 作用: 检查已有 NuRe 历史是否与检查点的采样时钟一致。
!===================================================================================================
  subroutine check_history()

    use commondata,only: timeUnit,outputSnapshotIntervalItc,NuReHistoryFile,NuReSampleCount,&
        ieee_is_finite
    implicit none

    ! 续算只追加新的 NuRe 样本；已有文件的行数必须等于已完成采样次数 NuReSampleCount。

    integer(kind=4) :: k, ios, n
    real(kind=8) :: values(12), lastTime
    character(len=512) :: line

    lastTime=-1.0d0
    n=0
    open(newunit=k,file=NuReHistoryFile,status='old',iostat=ios)
    if(ios/=0) error stop 'Restart requires the matching Nu/Re history'
    do
        read(k,'(a)',iostat=ios) line
        if(ios < 0) exit
        if(ios > 0) error stop 'Nu/Re history read error'
        if(len_trim(line)==0 .OR. line(1:1)=='#') cycle
        read(line,*,iostat=ios) values
        if(ios/=0 .OR. .NOT. all(ieee_is_finite(values))) error stop 'Invalid Nu/Re history row'
        if(values(1) <=lastTime) error stop 'Nonmonotone Nu/Re history'
        n=n+1
        lastTime=values(1)
    enddo

    close(k)
    if(n/=NuReSampleCount) error stop 'History sample count differs from checkpoint; use the matching history'
    if(n > 0) then

        if(abs(lastTime-(dble(n)*dble(outputSnapshotIntervalItc)/timeUnit)) > 1.0d-10*max(1.0d0,&
            (dble(n)*dble(outputSnapshotIntervalItc)/timeUnit))) error stop 'Checkpoint/history time mismatch'
    endif
  end subroutine check_history
!===================================================================================================

!===================================================================================================
! 子程序: average_window
! 作用: 对指定时间窗口内的历史数据积分，并检查窗口覆盖情况。
!===================================================================================================
  subroutine average_window(t0,t1,result,coverage)

    use commondata,only: NuReHistoryFile,ieee_is_finite
    implicit none

    real(kind=8), intent(in) :: t0, t1
    real(kind=8), intent(out) :: result(5), coverage
    integer(kind=4) :: k, ios
    real(kind=8) :: prev(12), row(12), aa, bb, ra(5), rb(5), left(5), right(5)
    logical :: hasPrev
    character(len=512) :: line

    result=0.0d0
    coverage=0.0d0
    hasPrev=.false.
    open(newunit=k,file=NuReHistoryFile,status='old')
    do
        read(k,'(a)',iostat=ios) line
        if(ios < 0) exit
        if(ios > 0) error stop 'Nu/Re history read failure'
        if(len_trim(line)==0 .OR. line(1:1)=='#') cycle
        read(line,*,iostat=ios) row
        if(ios/=0 .OR. .NOT. all(ieee_is_finite(row))) error stop 'Invalid history during averaging'
        if(hasPrev) then
            if(row(1) <=prev(1)) error stop 'Nonmonotone statistics time'
            aa=max(t0,prev(1))
            bb=min(t1,row(1))

            if(bb > aa) then
                ! 对 Re^2 积分再开根号，保持 sqrt(<u^2+v^2>_{V,t}) 的定义。
                left=prev(2:6)
                right=row(2:6)
                left(2)=left(2)**2
                right(2)=right(2)**2
                ra=left+(right-left)*(aa-prev(1))/(row(1)-prev(1))
                rb=left+(right-left)*(bb-prev(1))/(row(1)-prev(1))
                result=result+0.5d0*(ra+rb)*(bb-aa)
                coverage=coverage+bb-aa
            endif
        endif

        prev=row
        hasPrev=.true.
    enddo

    close(k)
    if(coverage > 0.0d0) then
        result=result/coverage
        result(2)=sqrt(max(0.0d0,result(2)))
    endif
  end subroutine average_window
!===================================================================================================

!===================================================================================================
! 子程序: output_unsteady_NuRe_postprocess
! 作用: 输出非稳态窗口平均及前后半窗口对比结果。
!===================================================================================================
  subroutine output_unsteady_NuRe_postprocess()
#ifdef unsteadyFlow

    use commondata,only: unsteadyAverageStartTf,unsteadyAverageEndTf,unsteadyAverageMidTf
    implicit none

    integer(kind=4) :: k
    real(kind=8) :: allMean(5), firstMean(5), lastMean(5), c0, c1, c2

    call average_window(unsteadyAverageStartTf,unsteadyAverageEndTf,allMean,c0)

    call average_window(unsteadyAverageStartTf,unsteadyAverageMidTf,firstMean,c1)

    call average_window(unsteadyAverageMidTf,unsteadyAverageEndTf,lastMean,c2)
    open(newunit=k,file='NuReStatistics_2DOpenaccMultiblock.dat',status='replace')
    write(k,*) 'Window t_ff:',unsteadyAverageStartTf,unsteadyAverageEndTf
    write(k,*) 'Covered durations, full/first/last:',c0,c1,c2
    if(abs(c0-(unsteadyAverageEndTf-unsteadyAverageStartTf)) > 1.0d-8 .OR. c1 <=0.0d0 .OR. c2 <=0.0d0) then
        write(k,*) 'INCOMPLETE: requested statistics window is not fully covered; no final mean is reported.'
    else

        write(k,*) 'Columns: NuVolAvg Re_rms_space_time Nu_hot Nu_cold Nu_middle'
        write(k,'(a,5ES24.16E3)') 'whole:',allMean
        write(k,'(a,5ES24.16E3)') 'first:',firstMean
        write(k,'(a,5ES24.16E3)') 'last :',lastMean
        write(k,'(a,5ES24.16E3)') 'relative half-window difference:',&
            abs(lastMean-firstMean)/max(abs(allMean),1.0d-30)
    endif

    close(k)
#endif
  end subroutine output_unsteady_NuRe_postprocess
!===================================================================================================

!===================================================================================================
! 子程序: owned_cell_area
! 作用: 累计节点在互不重叠物理分区中的面积。
!===================================================================================================
  subroutine owned_cell_area(x,y,dx,integrationBox,area)
    use commondata,only: centerIntegrationBox,refineRatio
    implicit none

    real(kind=8), intent(in) :: x, y, dx, integrationBox(4)
    real(kind=8), intent(out) :: area
    real(kind=8) :: xlo, xhi, ylo, yhi, middleWidth

    ! integrationBox 为当前存储区外框（粗区则为中心物理矩形）。
    ! 节点控制面积只向所属物理分区贡献；交换用的多存节点不另算一个区域。
    xlo=max(x-dx/2,integrationBox(1))
    xhi=min(x+dx/2,integrationBox(2))
    ylo=max(y-dx/2,integrationBox(3))
    yhi=min(y+dx/2,integrationBox(4))
    area=0.0d0
    if(xhi <=xlo .OR. yhi <=ylo) return

    if(refineRatio==1 .OR. dx > 1.0d0) then
        area=(xhi-xlo)*(yhi-ylo)    ! 单网格全域，或中心粗矩形。
        return
    endif

    ! 四个互不重叠的细物理矩形：左、右贯穿全高，下、上位于两者之间。
    ! 分别相加；没有“整个存储区积分再扣重叠区”的步骤。
    area=max(0.0d0,min(xhi,centerIntegrationBox(1))-xlo)*(yhi-ylo)    ! 左
    area=area+max(0.0d0,xhi-max(xlo,centerIntegrationBox(2)))*(yhi-ylo)    ! 右
    middleWidth=max(0.0d0,min(xhi,centerIntegrationBox(2))-max(xlo,centerIntegrationBox(1)))
    area=area+middleWidth*max(0.0d0,min(yhi,centerIntegrationBox(3))-ylo)    ! 下
    area=area+middleWidth*max(0.0d0,yhi-max(ylo,centerIntegrationBox(4)))    ! 上
  end subroutine owned_cell_area
!===================================================================================================

!===================================================================================================
! 子程序: section_owned_weight
! 作用: 将中线穿过的物理分区线段相加。
!===================================================================================================
  subroutine section_owned_weight(xOffset,yOffset,dx,integrationBox,k,axis,position,w)
    use commondata,only: centerIntegrationBox,refineRatio
    implicit none

    real(kind=8), intent(out) :: w    ! 输出：本节点在指定中线上实际负责的积分长度。

    ! axis=1 表示固定 x 的竖直中线，axis=2 表示固定 y 的水平中线。
    ! 先定位本节点控制区间；细区处于中心投影内时只累加两侧物理线段。

    real(kind=8), intent(in) :: xOffset, yOffset, dx, integrationBox(4), position
    integer(kind=4), intent(in) :: k, axis
    real(kind=8) :: q, lo, hi, clo, chi

    if(axis==1) then
        q=yOffset+(k-.5d0)*dx
        lo=integrationBox(3)
        hi=integrationBox(4)
        clo=centerIntegrationBox(3)
        chi=centerIntegrationBox(4)
    else
        q=xOffset+(k-.5d0)*dx
        lo=integrationBox(1)
        hi=integrationBox(2)
        clo=centerIntegrationBox(1)
        chi=centerIntegrationBox(2)
    endif

    w=max(0d0,min(q+dx/2,hi)-max(q-dx/2,lo))
    if(refineRatio > 1 .AND. dx==1d0) then
        if(axis==1) then
            if(position < centerIntegrationBox(1) .OR. position >=centerIntegrationBox(2)) return
        else
            if(position < centerIntegrationBox(3) .OR. position >=centerIntegrationBox(4)) return
        endif

        ! 中心投影内的细中线由两段物理区间相加，不计算交换层长度。
        w=max(0d0,min(q+dx/2,hi,clo)-max(q-dx/2,lo)) &
            +max(0d0,min(q+dx/2,hi)-max(q-dx/2,lo,chi))
    endif
  end subroutine section_owned_weight
!===================================================================================================

!===================================================================================================
! 子程序: fine_section_x
! 作用: 跨具名细数组读取 x 方向四点模板，计算中线温度、速度和温度导数。
!===================================================================================================
  subroutine fine_section_x(x,j,tm,um,dTdx)
    use commondata
    implicit none

    ! 模板位置按全域 nx 定位；四个温度和速度样本分别从相应具名细数组读取。

    real(kind=8), intent(in) :: x
    integer(kind=4), intent(in) :: j
    real(kind=8), intent(out) :: tm, um, dTdx
    real(kind=8) :: w(4), dw(4), temp(4), velocity(4)
    integer(kind=4) :: first, a

    call section_weights(x+0.5d0,nx,first,w,dw)
    do a=1,4
        call fine_value(first+a-1,j,T_left,T_right,T_bottom,T_top,temp(a))
        call fine_value(first+a-1,j,u_left,u_right,u_bottom,u_top,velocity(a))
    enddo

    tm=sum(w*temp)
    um=sum(w*velocity)
    dTdx=sum(dw*temp)
  end subroutine fine_section_x
!===================================================================================================

!===================================================================================================
! 子程序: fine_section_y
! 作用: 跨具名细数组读取 y 方向四点模板，计算中线温度、速度和温度导数。
!===================================================================================================
  subroutine fine_section_y(i,y,tm,vm,dTdy)
    use commondata
    implicit none

    ! 模板位置按全域 ny 定位；四个温度和速度样本分别从相应具名细数组读取。

    integer(kind=4), intent(in) :: i
    real(kind=8), intent(in) :: y
    real(kind=8), intent(out) :: tm, vm, dTdy
    real(kind=8) :: w(4), dw(4), temp(4), velocity(4)
    integer(kind=4) :: first, a

    call section_weights(y+0.5d0,ny,first,w,dw)
    do a=1,4
        call fine_value(i,first+a-1,T_left,T_right,T_bottom,T_top,temp(a))
        call fine_value(i,first+a-1,v_left,v_right,v_bottom,v_top,velocity(a))
    enddo

    tm=sum(w*temp)
    vm=sum(w*velocity)
    dTdy=sum(dw*temp)
  end subroutine fine_section_y
!===================================================================================================

!===================================================================================================
! 子程序: fine_value
! 作用: 将全局细节点编号转换成本地下标，读取对应细数组中的标量值。
!===================================================================================================
  subroutine fine_value(i,j,left,right,bottom,top,value)
    use commondata,only: f_left,f_right,f_bottom,f_top,&
        xOffsetRight,xOffsetBottom,xOffsetTop,yOffsetTop
    implicit none

    real(kind=8), intent(out) :: value    ! 输出：指定全局细节点对应的温度或速度等标量值。

    ! i/j 是全局细节点编号；减去所选具名数组的整数偏移后才是本地下标。
    ! 只用于主机端中线后处理；若请求落入未分配的中心空区，立即报错。

    integer(kind=4), intent(in) :: i, j
    real(kind=8), intent(in) :: left(size(f_left,1),size(f_left,2)), right(size(f_right,1),size(f_right,2)), &
        bottom(size(f_bottom,1),size(f_bottom,2)),top(size(f_top,1),size(f_top,2))

    if(i <=size(f_left,1)) then
        value=left(i,j)
    else if(i > nint(xOffsetRight)) then
        value=right(i-nint(xOffsetRight),j)
    else if(j <=size(f_bottom,2)) then
        value=bottom(i-nint(xOffsetBottom),j)
    else if(j > nint(yOffsetTop)) then
        value=top(i-nint(xOffsetTop),j-nint(yOffsetTop))
    else
        error stop 'Fine section requested an inactive central node'
    endif
  end subroutine fine_value
!===================================================================================================

!===================================================================================================
