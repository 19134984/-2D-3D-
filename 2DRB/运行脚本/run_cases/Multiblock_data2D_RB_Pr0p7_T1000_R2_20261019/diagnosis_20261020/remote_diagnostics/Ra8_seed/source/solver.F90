!=============================================================
!!!    注释区，代码描述
!!!    二维浮力驱动自然对流 OpenACC 静态多块网格版本
!!!    块内 D2Q9/D2Q5 算法保持不变
!!!    LBM方法
!!!    MRT-LBE
!
! 坐标与下标约定（全文统一使用以下名称）：
! 1. 局部下标 i,j：当前网格块数组中的整数编号；各场数组从 1 开始，f_post/g_post 另有外圈。
! 2. 全局物理坐标 x,y：节点在整个计算域中的位置，统一以细格距为长度单位，用实数表示。
!    x=xOffset+(i-0.5)*dx，y=yOffset+(j-0.5)*dx；这里的坐标尚未除以 lengthUnit。
!    xOffset/yOffset 是该网格块在全局物理坐标中的偏移量。
! 3. 全局细网格下标 I,J：统一细网格编号体系中的整数编号，对应 x=I-0.5、y=J-0.5。
!    当前网格的共址节点满足 I=x+0.5、J=y+0.5；全局细网格下标不等于任意块的局部下标。
! 默认例子：粗网格 dx=2、xOffset=122.5，局部下标 i=1 -> 全局物理坐标 x=123.5 -> I=124。
!           粗网格局部下标 i=3 与左细区局部下标 i=128 都位于 x=127.5，对应 I=128。
! 左细区偏移为零且 dx=1，局部下标恰好等于全局细网格下标；其他块需要按偏移和格距换算。
! 插值另用连续下标 q=(x-xOffset)/dx+0.5；例如 x=128.5 对应粗网格 q=3.5。
! 连续下标用于定位节点之间的位置，可能不是整数，不能直接作数组下标；y 方向同理。
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
        integer(kind=4), parameter :: nx=512, ny=512    ! 全域格子长度，必须能被 refineRatio 整除
        integer(kind=4), parameter :: refineRatio=2    ! 粗/细格距之比，取2、4、8等2的幂
        ! 接收层取2层。
        ! 粗、细网格各延伸 interfaceSkin 个本地格距。
        ! 换成细格距单位时，粗网格延伸宽度为 interfaceSkin*refineRatio，细网格为 interfaceSkin。
        integer(kind=4), parameter :: interfaceSkin=2

        integer(kind=4), parameter :: fineLayerCellsLeft=64    ! 从左墙向内的细节点编号，x = Left-0.5
        integer(kind=4), parameter :: fineLayerCellsRight=65    ! 从右墙向内的细节点编号，x = nx-Right+0.5
        integer(kind=4), parameter :: fineLayerCellsBottom=64    ! 从下墙向内的细节点编号，y = Bottom-0.5
        integer(kind=4), parameter :: fineLayerCellsTop=65    ! 从上墙向内的细节点编号，y = ny-Top+0.5
        !===============================================================================================

        !===============================================================================================
        ! 是否在计算前从旧算例重启
        !===============================================================================================
        integer(kind=4), parameter :: loadInitField=0   ! 0: 重新计算；1: 续算，优先读取latest.meta

        ! 正常断电续算只需要设置 loadInitField=1；
        ! 优先读取<reloadFilePrefix>-latest.meta，其中直接保存最新bin文件名、配置和累计计数。
        ! 只有latest.meta不存在时，才使用手工填写的reloadFileName和reloadItcOffset；步数需与所选文件一致。
        character(len=256) :: reloadFileName=""        ! 手工续算文件名，例如reloadFile2DOpenaccMultiblock-100tff-100000itc.bin
        integer(kind=4) :: reloadItcOffset=0           ! 续算前累计细格子步数；有latest.meta时读取其中的值，没有时在此手工填写

        ! 每次保存bin后覆盖同一个latest.meta，不另存各时刻的meta；bin只存数组。
        real(kind=8) :: reloadDimensionlessTime=0.0d0   ! 续算前累计t_ff，由reloadItcOffset/timeUnit计算，不需要手工填写
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

        real(kind=8), parameter :: Rayleigh=1.0d8
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
        real(kind=8), parameter :: dxFine=1.0d0
        ! 粗网格dt=dxCoarse；保持当前松弛率换算，细网格参数仍按基准定义。
        real(kind=8), parameter :: dxCoarse=dble(refineRatio)
        real(kind=8), parameter :: SnuCoarse=1.0d0/(0.5d0+(1.0d0/Snu-0.5d0)/dxCoarse)
        real(kind=8), parameter :: SqCoarse=1.0d0/(0.5d0+(1.0d0/Sq-0.5d0)/dxCoarse)
        real(kind=8), parameter :: QkCoarse=1.0d0/(0.5d0+(1.0d0/Qk-0.5d0)/dxCoarse)
        real(kind=8), parameter :: QnuCoarse=1.0d0/(0.5d0+(1.0d0/Qnu-0.5d0)/dxCoarse)
        !===============================================================================================

        !===============================================================================================
        ! 输出/备份相关设置（以自由落体时间 t_ff 为单位）
        !===============================================================================================
        real(kind=8), parameter :: epsU=1.0d-7, epsT=1.0d-7    ! 稳态收敛阈值

#ifdef steadyFlow
        real(kind=8), parameter :: outputSnapshotInterval=10.0d0    ! 快照和 Nu/Re 时间序列采样间隔（单位：t_ff）
        real(kind=8), parameter :: reloadFileInterval=100.0d0    ! f/g 重启文件输出间隔（单位：t_ff）
        real(kind=8), parameter :: outputPltFileInterval=100.0d0    ! Tecplot 文件输出间隔（单位：t_ff）
        integer(kind=4), parameter :: outputSnapshotFile=0    ! 是否输出后处理快照文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputPltFile=0    ! 是否输出 plt 文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputReloadFile=1    ! 是否周期输出 f/g 重启文件：0=不输出，1=输出
        ! ceiling 向上取整，保证终止步数为 refineRatio 的整数倍。
        ! 稳态：最大格子步，实际仍由 errorU/errorT 提前停止
        integer(kind=4), parameter :: itc_max=ceiling(20000000.0d0/dble(refineRatio))*refineRatio
        integer(kind=4), parameter :: steadyCheckIntervalItc=ceiling(2000.0d0/dble(refineRatio))*refineRatio
#endif

#ifdef unsteadyFlow
        real(kind=8), parameter :: outputSnapshotInterval=0.5d0    ! uvTrho 快照和 Nu/Re 时间序列采样间隔（单位：t_ff）
        real(kind=8), parameter :: reloadFileInterval=100.0d0    ! f/g 重启文件输出间隔（单位：t_ff）
        real(kind=8), parameter :: outputPltFileInterval=100.0d0    ! Tecplot 文件周期输出间隔（单位：t_ff）
        real(kind=8), parameter :: unsteadyRunDuration=30.0d0    ! 非稳态总目标时长，续算时只补足到该 t_ff
        integer(kind=4), parameter :: outputSnapshotFile=0    ! 是否输出后处理快照文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputPltFile=0    ! 是否输出 plt 文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputReloadFile=1    ! 是否周期输出 f/g 重启文件：0=不输出，1=输出
#endif

        ! ceiling向上取整：t_ff先换成粗步并向上取整，再换回细步，所有间隔均落在同步时刻。
        integer(kind=4), parameter :: outputSnapshotIntervalItc= &
            ceiling(outputSnapshotInterval*timeUnit/dble(refineRatio))*refineRatio
        integer(kind=4), parameter :: reloadFileIntervalItc= &
            ceiling(reloadFileInterval*timeUnit/dble(refineRatio))*refineRatio
        integer(kind=4), parameter :: outputPltFileIntervalItc= &
            ceiling(outputPltFileInterval*timeUnit/dble(refineRatio))*refineRatio
#ifdef unsteadyFlow
        ! 累计目标向上对齐到采样步；可在续算时增大，不属于物理兼容性约束。
        integer(kind=4), parameter :: itc_max=ceiling(unsteadyRunDuration*timeUnit/ &
            dble(max(1,outputSnapshotIntervalItc)))*max(1,outputSnapshotIntervalItc)
#endif

        character(len=*), parameter :: snapshotFilePrefix="buoyancyCavity2DOpenaccMultiblockSnapshot"
        ! 快照输出文件前缀（实际文件名形如：<snapshotFilePrefix>-<累计t_ff>tff-<累计细步数>itc.bin）

        character(len=*), parameter :: pltFolderPrefix="buoyancyCavity2DOpenaccMultiblockTecplot"
        ! Tecplot 输出文件前缀（实际文件名形如：<pltFolderPrefix>-<累计t_ff>tff-<累计细步数>itc.dat）

        character(len=*), parameter :: reloadFilePrefix="reloadFile2DOpenaccMultiblock"
        ! 重启文件前缀；文件名为<prefix>-<累计t_ff>tff-<累计细步数>itc.bin，latest.meta记录最新文件名

        character(len=*), parameter :: settingsFile="SimulationSettings2DOpenaccMultiblock.txt"
        character(*), parameter :: NuReHistoryFile='NuRe_2DOpenaccMultiblock.dat'    ! 按采样时间记录 Nu、Re 及场量统计，续算时追加
        !===============================================================================================

        !===============================================================================================
        ! 格子方向
        !===============================================================================================
        integer(kind=4) :: ex(0:8)=[0,1,0,-1,0,1,-1,-1,1]
        integer(kind=4) :: ey(0:8)=[0,0,1,0,-1,1,1,-1,-1]
        real(kind=8) :: omega(0:8), omegaT(0:4)

        ! itc 为本次运行已完成的细时间步数，每次启动从 0 开始。
        ! reloadItcOffset 为续算前累计细步数；总步数为 reloadItcOffset+itc。
        ! 粗细比为 r 时，每个粗步使 itc 增加 r。
        integer(kind=4) :: itc=0
#ifdef unsteadyFlow
        integer(kind=4) :: unsteadyItcRemaining    ! 本次尚需推进的细格子步数；初始化时扣除续算已有步数。
#endif

        ! NuReSampleCount 是已完成的 Nu/Re 采样次数，不是物理时间。
        ! 初值为 0，每次 calNuRe 增加 1；续算时恢复累计次数，实际 t_ff 为 (reloadItcOffset+itc)/timeUnit。
        integer(kind=4) :: NuReSampleCount=0

        real(kind=8) :: errorU=100.0d0, errorT=100.0d0
        !===============================================================================================

        !===============================================================================================
        ! 粗网格相应的信息。
        !===============================================================================================
        ! 左右细区贯穿全高，包含四角区域；上下细区仅存左右细区之间的节点。
        ! 各场量的局部下标 i,j 从 1 开始；对应的全局物理坐标为 xOffset+(i-0.5)*dx、yOffset+(j-0.5)*dx。
        ! f_post/g_post 另留 0 和 n+1 外围位置，用于相邻细区的迁移数据交换。


        ! 中心粗网格负责积分的矩形边界，依次为 xmin,xmax,ymin,ymax。
        ! 默认其值为 [127.5,895.5,127.5,895.5]。不包含重叠层
        ! 默认 dxCoarse=2，首末粗节点正好在分区线上：它们的积分权重各为 1，内部节点为 2。
        real(kind=8) :: centerIntegrationBox(4)



        ! 粗网格在全局物理坐标中的偏移量；局部下标 i=1 的节点位置还需加 0.5*dxCoarse。
        real(kind=8) :: xOffsetCoarse, yOffsetCoarse

        logical :: coarseIsWall(4)    ! 粗网格左、右、下、上是否执行物理域边界条件；周期边界也为真。
        ! 普通节点和接口接收层始终保存 f/g；collision/collisionT 统一由分布函数计算矩。
        ! 各自迁移后恢复分布函数；初始化、同步点、输出及续算时五区 f/g 均为分布函数。
        real(kind=8), allocatable :: f_coarse(:,:,:)    ! 粗网格 D2Q9 流场分布函数；第三维 0:8 为速度方向。
        real(kind=8), allocatable :: g_coarse(:,:,:)    ! 粗网格 D2Q5 温度分布函数；第三维 0:4 为速度方向。
        real(kind=8), allocatable :: f_post_coarse(:,:,:)    ! 碰撞后的流场分布函数；外围 0、n+1 供迁移。
        real(kind=8), allocatable :: g_post_coarse(:,:,:)    ! 碰撞后的温度分布函数；外围 0、n+1 供迁移。

        real(kind=8), allocatable :: rho_coarse(:,:)    ! 粗网格当前密度。
        real(kind=8), allocatable :: u_coarse(:,:)    ! 粗网格当前 x 方向速度。
        real(kind=8), allocatable :: v_coarse(:,:)    ! 粗网格当前 y 方向速度。
        real(kind=8), allocatable :: T_coarse(:,:)    ! 粗网格当前温度。
        ! 五个区域统一保存单位时间的F_x、F_y，不含dt或碰撞系数(1-s/2)。
        ! 时间步在碰撞、速度恢复和接口矩转换中显式相乘；粗块dt=dxCoarse，细块dt=dxFine=1。
        real(kind=8), allocatable :: Fx_coarse(:,:)    ! 无磁场时为0；SideHeatedHa模式另含磁场贡献。
        real(kind=8), allocatable :: Fy_coarse(:,:)    ! 单位时间的y方向力；无磁场时为rho_coarse*gBeta*(T_coarse-Tref)。

#ifdef steadyFlow
        real(kind=8), allocatable :: up_coarse(:,:)    ! 上一次稳态检查时的 x 方向速度。
        real(kind=8), allocatable :: vp_coarse(:,:)    ! 上一次稳态检查时的 y 方向速度。
        real(kind=8), allocatable :: Tp_coarse(:,:)    ! 上一次稳态检查时的温度。
#endif
        !===============================================================================================

        !===============================================================================================
        ! 左边细网格相应的信息。
        !===============================================================================================
        ! 四个细区按局部下标遍历各自数组；Offset 保存各区在全局物理坐标中的偏移量。
        ! 细区积分范围由中心矩形四边确定，节点面积在积分时计算。

        real(kind=8) :: xOffsetLeft, yOffsetLeft    ! 当前细区在全局物理坐标中的偏移量。
        logical :: leftIsWall(4)    ! 左、右、下、上是否为物理墙面。
        real(kind=8), allocatable :: f_left(:,:,:)    ! 流场分布函数，第三维为 D2Q9 方向编号 0:8。
        real(kind=8), allocatable :: g_left(:,:,:)    ! 温度分布函数，第三维为 D2Q5 方向编号 0:4。
        real(kind=8), allocatable :: f_post_left(:,:,:)    ! 碰撞后的流场分布函数；前两维另留 0、n+1 外围位置供细区迁移交换。
        real(kind=8), allocatable :: g_post_left(:,:,:)    ! 碰撞后的温度分布函数；前两维另留 0、n+1 外围位置供细区迁移交换。

        real(kind=8), allocatable :: rho_left(:,:)    ! 当前密度。
        real(kind=8), allocatable :: u_left(:,:)    ! 当前 x 方向速度。
        real(kind=8), allocatable :: v_left(:,:)    ! 当前 y 方向速度。
        real(kind=8), allocatable :: T_left(:,:)    ! 当前温度。
        real(kind=8), allocatable :: Fx_left(:,:)    ! 单位时间的x方向力；无磁场时为0。
        real(kind=8), allocatable :: Fy_left(:,:)    ! 单位时间的y方向力；无磁场时为rho_left*gBeta*(T_left-Tref)。

#ifdef steadyFlow
        real(kind=8), allocatable :: up_left(:,:)    ! 上一次稳态检查时的 x 方向速度，用于计算收敛误差。
        real(kind=8), allocatable :: vp_left(:,:)    ! 上一次稳态检查时的 y 方向速度，用于计算收敛误差。
        real(kind=8), allocatable :: Tp_left(:,:)    ! 上一次稳态检查时的温度，用于计算收敛误差。
#endif
        !===============================================================================================

        !===============================================================================================
        ! 右边细网格相应的信息。
        !===============================================================================================

        real(kind=8) :: xOffsetRight, yOffsetRight    ! 当前细区在全局物理坐标中的偏移量。
        logical :: rightIsWall(4)    ! 左、右、下、上是否为物理墙面。
        real(kind=8), allocatable :: f_right(:,:,:)    ! 流场分布函数，第三维为 D2Q9 方向编号 0:8。
        real(kind=8), allocatable :: g_right(:,:,:)    ! 温度分布函数，第三维为 D2Q5 方向编号 0:4。
        real(kind=8), allocatable :: f_post_right(:,:,:)    ! 碰撞后的流场分布函数；前两维另留 0、n+1 外围位置供细区迁移交换。
        real(kind=8), allocatable :: g_post_right(:,:,:)    ! 碰撞后的温度分布函数；前两维另留 0、n+1 外围位置供细区迁移交换。

        real(kind=8), allocatable :: rho_right(:,:)    ! 当前密度。
        real(kind=8), allocatable :: u_right(:,:)    ! 当前 x 方向速度。
        real(kind=8), allocatable :: v_right(:,:)    ! 当前 y 方向速度。
        real(kind=8), allocatable :: T_right(:,:)    ! 当前温度。
        real(kind=8), allocatable :: Fx_right(:,:)    ! 单位时间的x方向力；无磁场时为0。
        real(kind=8), allocatable :: Fy_right(:,:)    ! 单位时间的y方向力；无磁场时为rho_right*gBeta*(T_right-Tref)。

#ifdef steadyFlow
        real(kind=8), allocatable :: up_right(:,:)    ! 上一次稳态检查时的 x 方向速度，用于计算收敛误差。
        real(kind=8), allocatable :: vp_right(:,:)    ! 上一次稳态检查时的 y 方向速度，用于计算收敛误差。
        real(kind=8), allocatable :: Tp_right(:,:)    ! 上一次稳态检查时的温度，用于计算收敛误差。
#endif
        !===============================================================================================

        !===============================================================================================
        ! 下边细网格相应的信息。
        !===============================================================================================

        real(kind=8) :: xOffsetBottom, yOffsetBottom    ! 当前细区在全局物理坐标中的偏移量。
        logical :: bottomIsWall(4)    ! 左、右、下、上是否为物理墙面。
        real(kind=8), allocatable :: f_bottom(:,:,:)    ! 流场分布函数，第三维为 D2Q9 方向编号 0:8。
        real(kind=8), allocatable :: g_bottom(:,:,:)    ! 温度分布函数，第三维为 D2Q5 方向编号 0:4。
        real(kind=8), allocatable :: f_post_bottom(:,:,:)    ! 碰撞后的流场分布函数；前两维另留 0、n+1 外围位置供细区迁移交换。
        real(kind=8), allocatable :: g_post_bottom(:,:,:)    ! 碰撞后的温度分布函数；前两维另留 0、n+1 外围位置供细区迁移交换。

        real(kind=8), allocatable :: rho_bottom(:,:)    ! 当前密度。
        real(kind=8), allocatable :: u_bottom(:,:)    ! 当前 x 方向速度。
        real(kind=8), allocatable :: v_bottom(:,:)    ! 当前 y 方向速度。
        real(kind=8), allocatable :: T_bottom(:,:)    ! 当前温度。
        real(kind=8), allocatable :: Fx_bottom(:,:)    ! 单位时间的x方向力；无磁场时为0。
        real(kind=8), allocatable :: Fy_bottom(:,:)    ! 单位时间的y方向力；无磁场时为rho_bottom*gBeta*(T_bottom-Tref)。

#ifdef steadyFlow
        real(kind=8), allocatable :: up_bottom(:,:)    ! 上一次稳态检查时的 x 方向速度，用于计算收敛误差。
        real(kind=8), allocatable :: vp_bottom(:,:)    ! 上一次稳态检查时的 y 方向速度，用于计算收敛误差。
        real(kind=8), allocatable :: Tp_bottom(:,:)    ! 上一次稳态检查时的温度，用于计算收敛误差。
#endif
        !===============================================================================================

        !===============================================================================================
        ! 上边细网格相应的信息。
        !===============================================================================================

        real(kind=8) :: xOffsetTop, yOffsetTop    ! 当前细区在全局物理坐标中的偏移量。
        logical :: topIsWall(4)    ! 左、右、下、上是否为物理墙面。
        real(kind=8), allocatable :: f_top(:,:,:)    ! 流场分布函数，第三维为 D2Q9 方向编号 0:8。
        real(kind=8), allocatable :: g_top(:,:,:)    ! 温度分布函数，第三维为 D2Q5 方向编号 0:4。
        real(kind=8), allocatable :: f_post_top(:,:,:)    ! 碰撞后的流场分布函数；前两维另留 0、n+1 外围位置供细区迁移交换。
        real(kind=8), allocatable :: g_post_top(:,:,:)    ! 碰撞后的温度分布函数；前两维另留 0、n+1 外围位置供细区迁移交换。

        real(kind=8), allocatable :: rho_top(:,:)    ! 当前密度。
        real(kind=8), allocatable :: u_top(:,:)    ! 当前 x 方向速度。
        real(kind=8), allocatable :: v_top(:,:)    ! 当前 y 方向速度。
        real(kind=8), allocatable :: T_top(:,:)    ! 当前温度。
        real(kind=8), allocatable :: Fx_top(:,:)    ! 单位时间的x方向力；无磁场时为0。
        real(kind=8), allocatable :: Fy_top(:,:)    ! 单位时间的y方向力；无磁场时为rho_top*gBeta*(T_top-Tref)。

#ifdef steadyFlow
        real(kind=8), allocatable :: up_top(:,:)    ! 上一次稳态检查时的 x 方向速度，用于计算收敛误差。
        real(kind=8), allocatable :: vp_top(:,:)    ! 上一次稳态检查时的 y 方向速度，用于计算收敛误差。
        real(kind=8), allocatable :: Tp_top(:,:)    ! 上一次稳态检查时的温度，用于计算收敛误差。
#endif
        !===============================================================================================

        !===============================================================================================
        ! 重叠层保存数组的信息。
        !===============================================================================================
        ! 只保存细接收节点所需的三个粗时刻；第二维为左、右、下、上四区，各区编号都从1开始。
        ! 每节点20个量：rho,u,v,T,Fx,Fy（单位时间力）、9个缩放流场非平衡矩、5个缩放温度非平衡矩。
        ! 空间插值直接读取粗场；时间插值在接收节点内完成，不另设中转数组。
        real(kind=8), allocatable :: coarseHistoryPrevious(:,:)    ! 前一粗时刻
        real(kind=8), allocatable :: coarseHistoryCurrent(:,:)     ! 当前粗时刻
        real(kind=8), allocatable :: coarseHistoryNext(:,:)        ! 下一粗时刻预测数据
        !===============================================================================================
    end module commondata

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
    logical :: outputExists
    character(len=40) :: outputTime
    character(len=256) :: outputName
    integer(kind=8) :: wallClockStart, wallClockEnd, wallClockRate
    integer(kind=4) :: lastSnapshotItc, lastPltItc, lastReloadItc    ! 本次最近输出的步数，防止结束时重复写同一步。
    real(kind=8) :: latticeUpdateCount    ! 统计实际节点更新次数。
    integer(kind=4) :: k    ! 一个粗时间步内的细网格子步编号。
    real(kind=8) :: timeInterpWeight(0:2)    ! 粗网格三个时间层的插值权重。

    !===============================================================================================
    ! 初始化 OpenACC 设备
    if(loadInitField.EQ.1) then
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')  !append,继续往后保存
        write(00,*) " "
        write(00,*) "================ Restart continuation begins ================"
    else
        open(unit=00,file=trim(settingsFile),status='replace')  !replace,覆盖原有内容
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

    ! 重新计算首先需要知道粗细网格的接收区的数据
    ! 多块计算：先从粗场保存 t=0 的接口数据，再填充双方的接收层。
    ! [0,1,0] 只取 coarseHistoryCurrent，不涉及中间时刻的时间插值。
    ! 续算直接使用检查点恢复的场量和三个时间层。
    if(loadInitField==0 .AND. refineRatio > 1) then
        call save_coarse_history(coarseHistoryCurrent) !四个细网格重叠层区域从粗网格获取的数据，被提前算好并保存起来，宏观量对，非平衡还没转换好。
        call fine_to_coarse()  !从细网格当前的数据中取值，经过粗细网格转换，填入中心粗网格的接收层，这里返回的是 fg。
        call coarse_to_fine([0.0d0,1.0d0,0.0d0]) !从粗网格当前的数据中取值，经过粗细网格转换，填入四个细网格的接收层 fg。
    endif

    lastSnapshotItc=-1    !记录本次最近一次输出的位置
    lastPltItc=-1         !初始设为 -1，表示本次程序启动后，还没有输出过这一类文件。
    lastReloadItc=-1

    if(loadInitField==1) then   !续算时
        ! 如果无需再推进，避免把恢复时刻已有的文件重复记入输出索引。
        ! 重启和 plt 的实际输出间隔如果不相同，得做这个检测。
        call output_time_string(outputTime)  !返回累计 tff 的有效字符串，100，100.5这种，outputTime是字符串
        write(outputName,'(a,a,a,i0,a)') snapshotFilePrefix//'-',trim(outputTime),'tff-',reloadItcOffset+itc,'itc.bin'
        inquire(file=trim(outputName),exist=outputExists) !检查这个时刻的快照文件是否已经存在
        if(outputExists) lastSnapshotItc=0 !本次运行的起点 itc=0，已经有对应的快照文件，可以当作这个时刻已经输出过了。
        write(outputName,'(a,a,a,i0,a)') pltFolderPrefix//'-',trim(outputTime),'tff-',reloadItcOffset+itc,'itc.dat'
        inquire(file=trim(outputName),exist=outputExists) !检查这个时刻的 plt 文件是否已经存在
        if(outputExists) lastPltItc=0 !本次运行的起点 itc=0，已经有对应的 plt 文件，可以当作这个时刻已经输出过了。
        lastReloadItc=0
    endif
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
    ! 每轮粗细同步后更新 itc；各个步数已在参数设置和初始化时确定。
#ifdef steadyFlow
    ! 稳态 itc_max 限制本次运行的步数，速度和温度误差均达标时提前停止。
    do while(((errorU .GT. epsU) .OR. (errorT .GT. epsT)) .AND. &
        (itc .LT. itc_max))
#endif
#ifdef unsteadyFlow
    ! 非稳态 itc_max 为整个算例的总目标，续算只补足剩余步数。
    do while(itc.LT.unsteadyItcRemaining)    !非稳态：续算时只推进到 unsteadyRunDuration 对应的总格子步
#endif

        ! 一次循环推进一个粗时间步，完成 refineRatio 个细步和接口交换。
        if(refineRatio==1) then  !单网格情况直接推进
            call advance_coarse() !中心区域推进
            itc=itc+refineRatio
        else
            ! 粗网格预测终点（以下 t+2 的注释以 refineRatio=2 为例）；细子步读取过去、起点和预测终点三个时间层。
            call advance_coarse() !粗网格：t+2, 细网格：t
            call save_coarse_history(coarseHistoryNext) !保存粗网格在 t+2 的接口数据
            do k=1,refineRatio
                call coarse_time_weights(dble(k-1)/dble(refineRatio),timeInterpWeight) !计算子步起点的时间权重
                call coarse_to_fine(timeInterpWeight) !准备当前细子步所需的接收层 f/g
                call advance_fine()  !细网格推进一步
            enddo   !此时粗网格：t+2, 细网格：t+2

            ! 同步时更新双方接收层的 f/g，供诊断、续算和下一轮粗网格直接推进。
            call fine_to_coarse()   !细网格已经到达 t+2，因此从细网格读取当前数据，重新填充粗网格接收层
            call coarse_to_fine([0.0d0,0.0d0,1.0d0]) !补齐细接收层的完整 f、g
            call rotate_coarse_history() !把时间历史向前移动一层
            itc=itc+refineRatio  !更新本次完成的细步数
        endif
        ! 此时粗细网格已经到达同一时刻，接收层已恢复为正常 f、g，供后面的诊断和文件输出。

#ifdef steadyFlow
        ! 到达检查间隔时计算收敛误差。
        if(MOD(reloadItcOffset+itc,steadyCheckIntervalItc) .EQ. 0) then
            call check()
        endif
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

        ! 到达续算文件的输出间隔时保存当前状态。
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
    latticeUpdateCount=dble(size(f_coarse,1))*dble(size(f_coarse,2))*dble(itc/refineRatio) !统计本次粗网格推进了多少个节点
    if(refineRatio > 1) then
        latticeUpdateCount=latticeUpdateCount+dble(itc)*&
            (dble(size(f_left,1))*dble(size(f_left,2))+dble(size(f_right,1))*dble(size(f_right,2))+&
             dble(size(f_bottom,1))*dble(size(f_bottom,2))+dble(size(f_top,1))*dble(size(f_top,2)))
    endif   !统计本次实际推进了多少个节点

    open(unit=00,file=trim(settingsFile),status='unknown',position='append')  !append在已有内容后追加
    write(00,*) "======================================================================"
    write(00,*) "Time (CPU) = ",timeEnd-timeStart,"s"
    if(timeEnd > timeStart) then
        write(00,*) "MLUPS = ",latticeUpdateCount/(timeEnd-timeStart)/1.0d6
    else
        write(00,*) "MLUPS = unavailable (CPU timer resolution)"
    endif
    write(00,*) "Time (ACC) = ",timeEnd2-timeStart2,"s"
    if(timeEnd2 > timeStart2) then
        write(00,*) "MLUPS (ACC) = ",latticeUpdateCount/(timeEnd2-timeStart2)/1.0d6
    else
        write(00,*) "MLUPS (ACC) = unavailable (wall timer resolution)"
    endif
    write(00,*) "Fine steps in this run = ",itc   !本次推进的细步数和节点更新次数
    write(00,*) "Lattice updates in this run = ",latticeUpdateCount
    close(00)

    write(*,*) 'Elapsed seconds:',timeEnd2-timeStart2
    !===============================================================================================

    !===============================================================================================
    ! 结束后的统计与最终输出
#ifdef steadyFlow
    ! 最后一步可能不在采样时刻；单独输出最终诊断，不追加 NuRe 历史，也不增加 NuReSampleCount。
    call output_final_NuRe()
#endif

    ! 保留各输出开关；若最后一步已做过周期输出，就不再重复写同一时刻。
    ! 如果最后一步已经在循环中输出过，就不再重复输出；否则补一份终止时刻的数据。
    if((outputPltFile .EQ. 1) .AND. (lastPltItc/=itc)) call output_Tecplot()

    if((outputSnapshotFile .EQ. 1) .AND. (lastSnapshotItc/=itc)) call output_SnapshotFile()

    ! 本次尚未保存终止步时，补写最后一份续算文件。
    if((outputReloadFile .EQ. 1) .AND. (lastReloadItc/=itc)) call output_ReloadFile()
    !===============================================================================================

    !===============================================================================================
    ! 先释放设备数组，再释放主机上的五套网格、时间层和接口数组。
    open(unit=00,file=trim(settingsFile),status='unknown',position='append')
    write(00,*) "Deallocate Array......"
    call exit_data_2d_openacc()  !删除设备数据
    deallocate(coarseHistoryPrevious,coarseHistoryCurrent,coarseHistoryNext) !释放主机上的动态数组
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

#ifdef steadyFlow
    if(errorU .LE. epsU .AND. errorT .LE. epsT) then
        write(00,*) "Steady calculation converged."
    else   !已经达到最大同步步数，但尚未达到稳态收敛要求。
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
    ! nx/ny 是按细格距计的全域长度；各区域场数组的局部下标都从 1 开始。
    ! 局部下标 i 对应的全局物理坐标为 xOffset+(i-0.5)*dx；centerIntegrationBox 保存全局物理坐标分界。

    integer(kind=4) :: k, i, j, alpha, bufferLength, ios, n
    real(kind=8) :: xLeft, xRight, yBottom, yTop
    real(kind=8) :: x, y, us2, un(0:8)
    real(kind=8) :: values(12)
    real(kind=8), allocatable :: history(:,:)
#ifdef steadyFlow
    real(kind=8), allocatable :: convergence(:,:)
#endif
    integer(kind=4) :: sample, step
    character(len=512) :: line

    ! 1. 集中检查输入参数。
    ! 粗细比不能小于1
    if(refineRatio < 1) error stop 'refineRatio must be a positive integer'

    ! 单网格中线诊断使用四点模板；多块分支在下面检查中心节点数和各数组宽度。
    if(refineRatio==1 .AND. min(nx,ny)<4) error stop 'At least four nodes per direction are required'

    ! 初始场加载只能为0或1
    if(loadInitField/=0 .AND. loadInitField/=1) error stop 'loadInitField must be 0 or 1'

    ! 正时间间隔由模块中的ceiling保证至少一个粗步。
    if(min(outputSnapshotInterval,reloadFileInterval,outputPltFileInterval) <= 0.0d0) &
        error stop 'Output intervals must be positive'

    ! 参数a范围检查
    if(paraA <=-4.0d0 .OR. paraA >=1.0d0) error stop 'Legacy paraA must be in (-4,1)'

    if(refineRatio > 1) then
        ! 检查界面层数
        if(interfaceSkin < 2) error stop 'Split flow/thermal advance requires at least two interface layers'

        !检查全域网格整除
        if(mod(nx,refineRatio)/=0 .OR. mod(ny,refineRatio)/=0) &
            error stop 'nx,ny must be multiples of refineRatio'

        ! 最外粗接收节点对应的全局细网格下标必须>=1，例如Left-interfaceSkin*refineRatio>=1。
        if(min(fineLayerCellsLeft,fineLayerCellsRight,fineLayerCellsBottom,fineLayerCellsTop) <= &
            interfaceSkin*refineRatio) &
            error stop 'Refined layer must contain the outermost coarse receiver node'

        ! 检查粗网格还有空间，确保重叠层有足够的空间
        if(min((nx-fineLayerCellsLeft-fineLayerCellsRight-2*interfaceSkin), &
            ny-(fineLayerCellsBottom+interfaceSkin)-(fineLayerCellsTop+interfaceSkin)) < 1) &
            error stop 'Compact fine ring requires a positive interior storage hole'

        !中心宽高必须能被 refineRatio 整除
        if(mod(nx-fineLayerCellsLeft-fineLayerCellsRight+1,refineRatio)/=0 .OR. &
            mod(ny-fineLayerCellsBottom-fineLayerCellsTop+1,refineRatio)/=0) &
            error stop 'Central width and height must be multiples of refineRatio; adjust the four interface indices'

        ! 中心粗网格在x、y两个方向上，是否都有至少4个可用于空间插值的节点。
        if(min(nx-fineLayerCellsLeft-fineLayerCellsRight+1, &
            ny-fineLayerCellsBottom-fineLayerCellsTop+1)/refineRatio < 3) &
            error stop 'Central coarse region needs at least four nodes per direction'
    endif

    ! 2. 设置本次计数、格子权重和各区坐标；固定输出步数已在模块中定义。
    itc=0
    ! 重新计算时清零累计步数；续算时保留手工设置，若有latest.meta再由其中的值覆盖。
    if(loadInitField==0) then
        reloadItcOffset=0
    endif
    reloadDimensionlessTime=0.0d0
    reloadMetadataLoaded=.false.


    ! 设置流场和温度场的格子权重。
    omega(0)=4.0d0/9.0d0
    omega(1:4)=1.0d0/9.0d0
    omega(5:8)=1.0d0/36.0d0
    omegaT(0)=(1.0d0-paraA)/5.0d0
    omegaT(1:4)=(paraA+4.0d0)/20.0d0

    ! 粗格距和粗时间步均为细网格的 refineRatio 倍。

    if(refineRatio==1) then    ! 单块网格情况,  centerIntegrationBox 为全域, 不存在细网格了
        centerIntegrationBox=[0.0d0,dble(nx),0.0d0,dble(ny)]
        xOffsetCoarse=0.0d0
        yOffsetCoarse=0.0d0
    else       !下面是多块网格的情况
        ! 确定四条分区线
        xLeft=dble(fineLayerCellsLeft)-0.5d0    !左分区线的全局细网格下标 I=128，对应全局物理坐标 xLeft=127.5
        xRight=dble(nx-fineLayerCellsRight)+0.5d0    !右分区线的全局细网格下标 I=896，对应全局物理坐标 xRight=895.5
        yBottom=dble(fineLayerCellsBottom)-0.5d0    !下分区线的全局细网格下标 J=128，对应全局物理坐标 yBottom=127.5
        yTop=dble(ny-fineLayerCellsTop)+0.5d0    !上分区线的全局细网格下标 J=896，对应全局物理坐标 yTop=895.5

        centerIntegrationBox=[xLeft,xRight,yBottom,yTop]    ! [127.5,895.5,127.5,895.5]

        ! 各区偏移和格距均使用细格距长度单位；以下节点位置均为全局物理坐标。
        !xOffsetCoarse = 122.5
        !粗区局部下标 i：  1       2       3       4
        !全局物理坐标 x： 123.5   125.5   127.5   129.5
        !全局细网格下标 I：124     126     128     130
        !用途：        接收层   接收层   分区线   中心内部
        !yOffsetCoarse = 122.5
        !粗区局部下标 j：  1       2       3       4
        !全局物理坐标 y： 123.5   125.5   127.5   129.5
        !全局细网格下标 J：124     126     128     130
        !用途：         接收层   接收层   分区线   中心内部
        xOffsetCoarse=xLeft-interfaceSkin*refineRatio-0.5d0*dxCoarse  !!xOffsetCoarse = 122.5
        yOffsetCoarse=yBottom-interfaceSkin*refineRatio-0.5d0*dxCoarse    !yOffsetCoarse = 122.5

        ! 四个细矩形拼成细环；仅向中心延伸 interfaceSkin 个细格距，细矩形和细矩形拼接处直接迁移。

        ! 左侧细区的坐标偏移。
        xOffsetLeft=0.0d0
        yOffsetLeft=0.0d0

        ! 右侧细区的坐标偏移，包含向中心延伸的 interfaceSkin 层。
        xOffsetRight=dble(nx-(fineLayerCellsRight+interfaceSkin))    !xOffsetRight = 893.0
        yOffsetRight=0.0d0

        ! 下侧细区从左侧重叠区之后开始。
        xOffsetBottom=dble((fineLayerCellsLeft+interfaceSkin))    !xOffsetBottom = 130.0
        yOffsetBottom=0.0d0

        ! 上侧细区从左侧重叠区之后开始，并向下延伸 interfaceSkin 层。
        xOffsetTop=dble((fineLayerCellsLeft+interfaceSkin))    !xOffsetTop = 130.0
        yOffsetTop=dble(ny-(fineLayerCellsTop+interfaceSkin))    !yOffsetTop = 893.0

    endif

    ! 单网格四边为物理边界；多块粗区四边均为接收层。
    if(refineRatio==1) then
        coarseIsWall=.true.
    else
        coarseIsWall=.false.
    endif

    if(refineRatio > 1) then

        ! 左右细区接触三面外墙，上下细区只接触各自的水平外墙。数组顺序为左、右、下、上
        leftIsWall=[.true.,.false.,.true.,.true.]
        rightIsWall=[.false.,.true.,.true.,.true.]
        bottomIsWall=[.false.,.false.,.true.,.false.]
        topIsWall=[.false.,.false.,.false.,.true.]
    endif


    ! 3. 按各区尺寸分配分布函数、宏观场和迁移临时数组；多块模式还分配四个细区，稳态模式另分配检查用的旧场。
    ! 同时分配并清零前一、当前、下一粗时刻的接口历史
    ! 分别为 coarse、left、right、bottom、top 分配数组，尺寸直接由分区位置和重叠厚度计算。
    ! 先分配 f，其余场量用 size(f,1/2) 取得同一区域的尺寸，不另存尺寸变量。
    ! 场量空间范围是 1:ni、1:nj；左右细区在x方向加 interfaceSkin，上下细区在y方向加 interfaceSkin。
    ! 壁面的全局物理坐标为0，相邻细网格首个流体节点的全局物理坐标为0.5
    ! f_post/g_post 为迁移额外保留 0 和 n+1，延续基准代码的半格壁面处理。
    ! 三个小数组保存接口时间历史；up/vp/Tp 保存上次稳态检查的场量。

    ! 分配中心粗网格的计算数组。
    if(refineRatio==1) then  ! 单块网格情况,  centerIntegrationBox 为全域, 不存在细网格了
        allocate(f_coarse(nx,ny,0:8))
    else     !多块模式下，x 方向的节点数，带有重叠层的节点数
        ! f_coarse(1:389,1:389,0:8)
        allocate(f_coarse((nx-fineLayerCellsLeft-fineLayerCellsRight+1)/refineRatio+1+2*interfaceSkin,&
            (ny-fineLayerCellsBottom-fineLayerCellsTop+1)/refineRatio+1+2*interfaceSkin,0:8))
    endif      !size 返回的是数某个维度的元素个数
    allocate(g_coarse(size(f_coarse,1),size(f_coarse,2),0:4))
    allocate(f_post_coarse(0:size(f_coarse,1)+1,0:size(f_coarse,2)+1,0:8))
    allocate(g_post_coarse(0:size(f_coarse,1)+1,0:size(f_coarse,2)+1,0:4))
    allocate(rho_coarse(size(f_coarse,1),size(f_coarse,2)))
    allocate(u_coarse(size(f_coarse,1),size(f_coarse,2)))
    allocate(v_coarse(size(f_coarse,1),size(f_coarse,2)))
    allocate(T_coarse(size(f_coarse,1),size(f_coarse,2)))
    allocate(Fx_coarse(size(f_coarse,1),size(f_coarse,2)))
    allocate(Fy_coarse(size(f_coarse,1),size(f_coarse,2)))

#ifdef steadyFlow
    allocate(up_coarse(size(f_coarse,1),size(f_coarse,2)))
    allocate(vp_coarse(size(f_coarse,1),size(f_coarse,2)))
    allocate(Tp_coarse(size(f_coarse,1),size(f_coarse,2)))
#endif


    if(refineRatio > 1) then
        ! 分配左侧细网格的计算数组。带有重叠层的节点数
        ! f_left(1:130,1:1024,0:8)
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

        ! 分配右侧细网格的计算数组。带有重叠层的节点数
        ! f_right(1:131,1:1024,0:8)
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

        ! 分配下侧细网格的计算数组。y方向带有重叠层的节点数， x 方向只占左右细数组之间剩下的空隙，不再向左右增加重叠层
        ! f_bottom(1:763,1:130,0:8)
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

        ! 分配上侧细网格的计算数组。y方向带有重叠层的节点数， x 方向只占左右细数组之间剩下的空隙，不再向左右增加重叠层
        ! f_top(1:763,1:131,0:8)
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

    ! 左右接收区各 interfaceSkin*(ny-Bottom-Top) 个节点；上下各 interfaceSkin*size(f_bottom,1)。
    ! 默认左右1534、上下1526个节点，每节点20个量，最长30680；四区分别占一列，不拼接编号。
    if(refineRatio>1) then
        bufferLength=20*interfaceSkin*max(ny-fineLayerCellsBottom-fineLayerCellsTop,size(f_bottom,1))
        allocate(coarseHistoryPrevious(bufferLength,4),coarseHistoryCurrent(bufferLength,4), &
            coarseHistoryNext(bufferLength,4))
    else
        ! 单网格不交换，各留一个位置，统一设备分配和释放。
        allocate(coarseHistoryPrevious(1,1),coarseHistoryCurrent(1,1),coarseHistoryNext(1,1))
    endif
    coarseHistoryPrevious=0.0d0
    coarseHistoryCurrent=0.0d0
    coarseHistoryNext=0.0d0

    ! 4. 设置五个区域的初始场，并累计不重复的物理积分面积。

    ! 中心粗网格数组，包括它的重叠接收层的初始化，还有初始温度和初始fg的设置
    ! 初始化宏观量和碰撞后数组。
    u_coarse=0.0d0
    v_coarse=0.0d0
    rho_coarse=1.0d0
    Fx_coarse=0.0d0
    Fy_coarse=0.0d0
    f_post_coarse=0.0d0
    g_post_coarse=0.0d0

    ! 根据节点的全局物理坐标构造初始温度，并设置平衡分布函数。
    do j=1,size(f_coarse,2)
        y=(yOffsetCoarse+(dble(j)-0.5d0)*dxCoarse)/lengthUnit
        do i=1,size(f_coarse,1)
            x=(xOffsetCoarse+(dble(i)-0.5d0)*dxCoarse)/lengthUnit
#ifdef SideHeatedCell
            T_coarse(i,j)=Thot+x/(dble(nx)/lengthUnit)*(Tcold-Thot)
#else
            T_coarse(i,j)=Thot+y/(dble(ny)/lengthUnit)*(Tcold-Thot)+1d-3*sin(2d0*pi*x)*sin(pi*y)
#endif
            us2=u_coarse(i,j)*u_coarse(i,j)+v_coarse(i,j)*v_coarse(i,j)
            do alpha=0,8
                un(alpha)=u_coarse(i,j)*ex(alpha)+v_coarse(i,j)*ey(alpha)
                f_coarse(i,j,alpha)=rho_coarse(i,j)*omega(alpha)*&
                    (1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)  !D2Q9标准feq
            enddo

            do alpha=0,4
                un(alpha)=u_coarse(i,j)*ex(alpha)+v_coarse(i,j)*ey(alpha)
                g_coarse(i,j,alpha)=omegaT(alpha)*T_coarse(i,j)*(1.0d0+thermalGeqCoeff*un(alpha))
            enddo
        enddo
    enddo

#ifdef steadyFlow
    ! 稳态检查的初始参考场。
    up_coarse=u_coarse
    vp_coarse=v_coarse
    Tp_coarse=T_coarse
#endif

    !多块细网格
    if(refineRatio > 1) then  
        ! 初始化左侧细网格，包括重叠层，并设置初始温度和初始fg
        ! 初始化宏观量和碰撞后数组。
        u_left=0.0d0
        v_left=0.0d0
        rho_left=1.0d0
        Fx_left=0.0d0
        Fy_left=0.0d0
        f_post_left=0.0d0
        g_post_left=0.0d0

        ! 根据节点的全局物理坐标构造初始温度，并设置平衡分布函数。
        do j=1,size(f_left,2)
            y=(yOffsetLeft+(dble(j)-0.5d0)*dxFine)/lengthUnit
            do i=1,size(f_left,1)
                x=(xOffsetLeft+(dble(i)-0.5d0)*dxFine)/lengthUnit
#ifdef SideHeatedCell
                T_left(i,j)=Thot+x/(dble(nx)/lengthUnit)*(Tcold-Thot)
#else
                T_left(i,j)=Thot+y/(dble(ny)/lengthUnit)*(Tcold-Thot)+1d-3*sin(2d0*pi*x)*sin(pi*y)
#endif
                us2=u_left(i,j)*u_left(i,j)+v_left(i,j)*v_left(i,j)
                do alpha=0,8
                    un(alpha)=u_left(i,j)*ex(alpha)+v_left(i,j)*ey(alpha)
                    f_left(i,j,alpha)=rho_left(i,j)*omega(alpha)*&
                        (1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)  !D2Q9标准feq
                enddo

                do alpha=0,4
                    un(alpha)=u_left(i,j)*ex(alpha)+v_left(i,j)*ey(alpha)
                    g_left(i,j,alpha)=omegaT(alpha)*T_left(i,j)*(1.0d0+thermalGeqCoeff*un(alpha))
                enddo
            enddo
        enddo

#ifdef steadyFlow
        ! 稳态检查的初始参考场。
        up_left=u_left
        vp_left=v_left
        Tp_left=T_left
#endif

        ! 初始化右侧细网格，包括重叠层，并设置初始温度和初始fg
        ! 初始化宏观量和碰撞后数组。
        u_right=0.0d0
        v_right=0.0d0
        rho_right=1.0d0
        Fx_right=0.0d0
        Fy_right=0.0d0
        f_post_right=0.0d0
        g_post_right=0.0d0

        ! 根据节点的全局物理坐标构造初始温度，并设置平衡分布函数。
        do j=1,size(f_right,2)
            y=(yOffsetRight+(dble(j)-0.5d0)*dxFine)/lengthUnit
            do i=1,size(f_right,1)
                x=(xOffsetRight+(dble(i)-0.5d0)*dxFine)/lengthUnit
#ifdef SideHeatedCell
                T_right(i,j)=Thot+x/(dble(nx)/lengthUnit)*(Tcold-Thot)
#else
                T_right(i,j)=Thot+y/(dble(ny)/lengthUnit)*(Tcold-Thot)+1d-3*sin(2d0*pi*x)*sin(pi*y)
#endif
                us2=u_right(i,j)*u_right(i,j)+v_right(i,j)*v_right(i,j)
                do alpha=0,8
                    un(alpha)=u_right(i,j)*ex(alpha)+v_right(i,j)*ey(alpha)
                    f_right(i,j,alpha)=rho_right(i,j)*omega(alpha)*&
                        (1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)  !D2Q9标准feq
                enddo

                do alpha=0,4
                    un(alpha)=u_right(i,j)*ex(alpha)+v_right(i,j)*ey(alpha)
                    g_right(i,j,alpha)=omegaT(alpha)*T_right(i,j)*(1.0d0+thermalGeqCoeff*un(alpha))
                enddo
            enddo
        enddo

#ifdef steadyFlow
        ! 稳态检查的初始参考场。
        up_right=u_right
        vp_right=v_right
        Tp_right=T_right
#endif

        ! 初始化下侧细网格，包括y方向重叠层，并设置初始温度和初始fg
        ! 初始化宏观量和碰撞后数组。
        u_bottom=0.0d0
        v_bottom=0.0d0
        rho_bottom=1.0d0
        Fx_bottom=0.0d0
        Fy_bottom=0.0d0
        f_post_bottom=0.0d0
        g_post_bottom=0.0d0

        ! 根据节点的全局物理坐标构造初始温度，并设置平衡分布函数。
        do j=1,size(f_bottom,2)
            y=(yOffsetBottom+(dble(j)-0.5d0)*dxFine)/lengthUnit
            do i=1,size(f_bottom,1)
                x=(xOffsetBottom+(dble(i)-0.5d0)*dxFine)/lengthUnit
#ifdef SideHeatedCell
                T_bottom(i,j)=Thot+x/(dble(nx)/lengthUnit)*(Tcold-Thot)
#else
                T_bottom(i,j)=Thot+y/(dble(ny)/lengthUnit)*(Tcold-Thot)+1d-3*sin(2d0*pi*x)*sin(pi*y)
#endif
                us2=u_bottom(i,j)*u_bottom(i,j)+v_bottom(i,j)*v_bottom(i,j)
                do alpha=0,8
                    un(alpha)=u_bottom(i,j)*ex(alpha)+v_bottom(i,j)*ey(alpha)
                    f_bottom(i,j,alpha)=rho_bottom(i,j)*omega(alpha)*&
                        (1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)  !D2Q9标准feq
                enddo

                do alpha=0,4
                    un(alpha)=u_bottom(i,j)*ex(alpha)+v_bottom(i,j)*ey(alpha)
                    g_bottom(i,j,alpha)=omegaT(alpha)*T_bottom(i,j)*(1.0d0+thermalGeqCoeff*un(alpha))
                enddo
            enddo
        enddo

#ifdef steadyFlow
        ! 稳态检查的初始参考场。
        up_bottom=u_bottom
        vp_bottom=v_bottom
        Tp_bottom=T_bottom
#endif

        ! 初始化上侧细网格，包括重叠层
        ! 初始化宏观量和碰撞后数组。
        u_top=0.0d0
        v_top=0.0d0
        rho_top=1.0d0
        Fx_top=0.0d0
        Fy_top=0.0d0
        f_post_top=0.0d0
        g_post_top=0.0d0

        ! 根据节点的全局物理坐标构造初始温度，并设置平衡分布函数。
        do j=1,size(f_top,2)
            y=(yOffsetTop+(dble(j)-0.5d0)*dxFine)/lengthUnit
            do i=1,size(f_top,1)
                x=(xOffsetTop+(dble(i)-0.5d0)*dxFine)/lengthUnit
#ifdef SideHeatedCell
                T_top(i,j)=Thot+x/(dble(nx)/lengthUnit)*(Tcold-Thot)
#else
                T_top(i,j)=Thot+y/(dble(ny)/lengthUnit)*(Tcold-Thot)+1d-3*sin(2d0*pi*x)*sin(pi*y)
#endif
                us2=u_top(i,j)*u_top(i,j)+v_top(i,j)*v_top(i,j)
                do alpha=0,8
                    un(alpha)=u_top(i,j)*ex(alpha)+v_top(i,j)*ey(alpha)
                    f_top(i,j,alpha)=rho_top(i,j)*omega(alpha)*&
                        (1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)  !D2Q9标准feq
                enddo

                do alpha=0,4
                    un(alpha)=u_top(i,j)*ex(alpha)+v_top(i,j)*ey(alpha)
                    g_top(i,j,alpha)=omegaT(alpha)*T_top(i,j)*(1.0d0+thermalGeqCoeff*un(alpha))
                enddo
            enddo
        enddo

#ifdef steadyFlow
        ! 稳态检查的初始参考场。
        up_top=u_top
        vp_top=v_top
        Tp_top=T_top
#endif

    endif

    ! 5. 检查已构建的网格；续算时在这里读取并检查文件。

    ! 续算下，恢复各个区域的数据
    if(loadInitField==1) then
        call read_restart()    ! 恢复五个区域、三个接口时间层以及累计计数。

        ! 只保留检查点已经完成的采样。先完整读入并核对前缀，再覆盖原文件；旧尾部直接删除。
        ! 先把检查点之前需要保留的数据读入内存，确认正确后，再覆盖原文件，从而删除检查点之后的旧记录。
        allocate(history(12,NuReSampleCount))
        open(newunit=k,file=NuReHistoryFile,status='old',iostat=ios)   ! Nu/Re 历史
        if(ios/=0) error stop 'Restart requires the matching Nu/Re history'

        n=0
        do while(n<NuReSampleCount)  !循环只读取检查点对应的次数，不要求把整个文件读完
            read(k,'(a)',iostat=ios) line !# sample itc_total t_ff NuVolAvg ReVolRMS Nu_hot Nu_cold Nu_middle
                                          ! mass meanT Tmin Tmax rhoMin rhoMax
            if(ios/=0) error stop 'History is shorter than checkpoint sample count'
            if(len_trim(line)==0 .OR. line(1:1)=='#') cycle  !跳过，继续读下一行,不会增加 n，所以标题不算一次采样
            read(line,*,iostat=ios) sample,step,values !读取的是字符串 line，不是再次从文件读取下一行
            if(ios/=0 .OR. .NOT.all(ieee_is_finite(values))) error stop 'Invalid retained history row'

            n=n+1  ! 增加有效记录行数，核对采样编号和步数
            if(sample/=n .OR. step/=n*outputSnapshotIntervalItc) error stop 'History sample/step mismatch'
            if(abs(values(1)-dble(step)/timeUnit)>1d-12*max(1d0,abs(values(1)))) &
                error stop 'History time/step mismatch'  ! 再核对这一行的时间
            history(:,n)=values
        enddo
        close(k)

#ifdef steadyFlow
        ! 收敛检查频率可以不同；按累计检查步数保留对应的误差记录。
        allocate(convergence(3,reloadItcOffset/steadyCheckIntervalItc))
        open(newunit=k,file='Convergence_2DOpenaccMultiblock.dat',status='old',iostat=ios)
        if(ios/=0) error stop 'Restart requires convergence history'

        n=0
        do while(n<size(convergence,2))
            read(k,'(a)',iostat=ios) line  !itc_total  t_ff  errorU  errorT
            if(ios/=0) error stop 'Convergence history is shorter than checkpoint'
            if(len_trim(line)==0 .OR. line(1:1)=='#') cycle
            read(line,*,iostat=ios) step,values(1:3)
            if(ios/=0 .OR. .NOT.all(ieee_is_finite(values(1:3)))) error stop 'Invalid convergence row'

            n=n+1
            if(step/=n*steadyCheckIntervalItc) error stop 'Convergence step mismatch'
            if(abs(values(1)-dble(step)/timeUnit)>1d-12*max(1d0,abs(values(1)))) &
                error stop 'Convergence time mismatch'
            convergence(:,n)=values(1:3)
        enddo
        close(k)
#endif
        ! 两类前缀均检查成功才覆盖；次数与步数按整数写入。
        open(newunit=k,file=NuReHistoryFile,status='replace') !替换原文件
        write(k,'(a)') '# sample itc_total t_ff NuVolAvg ReVolRMS Nu_hot Nu_cold Nu_middle '// &
            'mass meanT Tmin Tmax rhoMin rhoMax' !抬头
        do n=1,NuReSampleCount
            write(k,'(2(I0,1X),12(ES24.16E3,1X))') n,n*outputSnapshotIntervalItc,history(:,n)
        enddo
        close(k)
        deallocate(history)

#ifdef steadyFlow
        open(newunit=k,file='Convergence_2DOpenaccMultiblock.dat',status='replace')
        write(k,'(a)') '# itc_total t_ff errorU errorT'
        do n=1,size(convergence,2)
            write(k,'(I0,1X,3(ES24.16E3,1X))') n*steadyCheckIntervalItc,convergence(:,n)
        enddo
        close(k)
        deallocate(convergence)
#endif

        call truncate_output_files() !截断快照，plt 以及FinalNuRe这三个输出文件到重启时刻

    endif  ! 续算处理结束

#ifdef unsteadyFlow
    ! 总目标 itc_max 已经对齐，且 read_restart 已检查 reloadItcOffset，因此相减后仍为整粗步。
    unsteadyItcRemaining=max(0,itc_max-reloadItcOffset)  !非稳态本次还需要跑多少步
#endif

    ! 6. 写出运行设置；全新计算建立输出文件，续算保留已有历史。
    ! 主程序已写入本次开始时间和设备信息，此处继续追加，避免覆盖启动日志。
    open(newunit=k,file=settingsFile,status='unknown',position='append')
    if(loadInitField==0) then
        ! 首次计算记录完整参数及各区域几何信息。
        write(k,*) 'Fine-equivalent nx,ny; refinement ratio:',nx,ny,refineRatio
        write(k,*) 'Rayleigh, Prandtl, Mach:',Rayleigh,Prandtl,Mach
        write(k,*) 'Fine tauf, Snu, Sq, Qk, Qnu, paraA:',tauf,Snu,Sq,Qk,Qnu,paraA
        write(k,*) 'Fine viscosity, diffusivity, gBeta, timeUnit:',viscosity,diffusivity,gBeta,timeUnit
        write(k,*) 'Target itc, actual t_ff =',itc_max,dble(itc_max)/timeUnit
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
        write(k,*) 'Geometry: central coarse rectangle and four compact fine rectangles.'
        write(k,*) 'Fine seams exchange post-collision populations directly, including diagonal links.'
        write(k,*) 'Coordinates: x=xOffset+(i-0.5)*dx, y=yOffset+(j-0.5)*dx; local i,j start at 1.'
        write(k,*) 'Interfaces preserve the connected-ring time interpolation and MRT scaling; dt=dx.'
        write(k,*) 'Restart: paired meta (configuration/counters) and bin (coarse,left,right,bottom,top arrays).'
        write(k,*) 'coarse: nx,ny,dx,xOffset,yOffset:',size(f_coarse,1),size(f_coarse,2),dxCoarse,xOffsetCoarse,yOffsetCoarse
        if(refineRatio > 1) then
            write(k,*) 'left: nx,ny,dx,xOffset,yOffset:',size(f_left,1),size(f_left,2),dxFine,xOffsetLeft,yOffsetLeft
            write(k,*) 'right: nx,ny,dx,xOffset,yOffset:',size(f_right,1),size(f_right,2),dxFine,xOffsetRight,yOffsetRight
            write(k,*) 'bottom: nx,ny,dx,xOffset,yOffset:',size(f_bottom,1),size(f_bottom,2),dxFine,xOffsetBottom,&
                yOffsetBottom
            write(k,*) 'top: nx,ny,dx,xOffset,yOffset:',size(f_top,1),size(f_top,2),dxFine,xOffsetTop,yOffsetTop
        endif
    else
        ! 续算只记录本次恢复的累计步数和tff，不重复写入完整设置。
        write(k,*) 'Reload offset itc =',reloadItcOffset
        write(k,*) 'Reload dimensionless time =',reloadDimensionlessTime
#ifdef unsteadyFlow
        write(k,*) 'Cumulative target itc, t_ff =',itc_max,dble(itc_max)/timeUnit
        write(k,*) 'Remaining fine steps =',unsteadyItcRemaining
#endif
    endif
    close(k)

    if(loadInitField==0) then         ! 新算例建立输出索引和历史标题
        open(newunit=k,file='OutputFiles_2DOpenaccMultiblock.dat',status='replace')
        close(k)
        open(newunit=k,file=NuReHistoryFile,status='replace')
        write(k,'(a)') '# sample itc_total t_ff NuVolAvg ReVolRMS Nu_hot Nu_cold Nu_middle mass meanT Tmin Tmax rhoMin rhoMax'
        close(k)   
#ifdef steadyFlow
        open(newunit=k,file='Convergence_2DOpenaccMultiblock.dat',status='replace')
        write(k,'(a)') '# itc_total t_ff errorU errorT'
        close(k)
#endif
    endif
  end subroutine initial
!===================================================================================================

!===================================================================================================
! 子程序: read_restart
! 作用: 检查网格及物理配置，恢复时间计数和五套网格状态。
!===================================================================================================
  subroutine read_restart()
    use commondata
    implicit none
    integer(kind=4) :: k, m, ios, head(9), sig(12), currentModel(12), geom(6)
    real(kind=8) :: phys(15), currentPhysics(15), coord(7), savedTime
#ifdef steadyFlow
    integer :: savedCheckInterval
#endif
    character(len=256) :: name, metaName, line
    logical :: exists

    ! latest存在时必须读取并核对；只有不存在时才使用手工设置，读取失败不能转入手工分支。
    metaName=reloadFilePrefix//'-latest.meta'
    inquire(file=trim(metaName),exist=exists)  !看文件是否存在
    if(exists) then   !文件存在
        open(newunit=m,file=trim(metaName),status='old',iostat=ios)
        if(ios/=0) error stop 'Cannot open latest.meta' !检查文件是否损坏

        read(m,'(a)',iostat=ios) line     !元数据格式
        if(ios/=0 .OR. trim(line)/='MB2D_RESTART 1') error stop 'Unsupported checkpoint metadata'

        read(m,'(a)',iostat=ios) name   !例如reloadFile2DOpenaccMultiblock-100tff-100000itc.bin
        if(ios/=0 .OR. len_trim(name)==0) error stop 'Missing checkpoint array filename'

        read(m,'(a)') line  !'(a)' 表示按字符读取     !# nx ny refineRatio Left Right Bottom Top interfaceSkin blockCount
        read(m,*,iostat=ios) head   !9个物理参数
        if(ios/=0) error stop 'Incomplete checkpoint mesh metadata'
        if(any(head/=[nx,ny,refineRatio,fineLayerCellsLeft,fineLayerCellsRight, & !any只要有一项不一致，就报错。
            fineLayerCellsBottom,fineLayerCellsTop,interfaceSkin,merge(5,1,refineRatio>1)])) &
            error stop 'Restart mesh/refinement mismatch'

        read(m,'(a)') line !自定义宏
        read(m,*,iostat=ios) sig
        if(ios/=0) error stop 'Incomplete checkpoint model metadata'

        read(m,'(a)') line !# Rayleigh Prandtl Mach Thot Tcold Snu Sq Qk Qnu paraA snapshotInterval reloadInterval pltInterval Ha phi
        read(m,*,iostat=ios) phys
        if(ios/=0 .OR. .NOT.all(ieee_is_finite(phys))) error stop 'Invalid checkpoint physics metadata'

        currentModel=0     !将数组的12个元素全部置零
#ifdef steadyFlow
        currentModel(1)=1     !根据编译开关，把对应元素设成1
#endif
#ifdef SideHeatedCell
        currentModel(3)=1
#endif
#ifdef HorizontalWallsNoslip
        currentModel(4)=1
#endif
#ifdef VerticalWallsNoslip
        currentModel(5)=1
#endif
#ifdef VerticalWallsPeriodicalU
        currentModel(6)=1
#endif
#ifdef HorizontalWallsConstT
        currentModel(7)=1
#endif
#ifdef HorizontalWallsAdiabatic
        currentModel(8)=1
#endif
#ifdef VerticalWallsConstT
        currentModel(9)=1
#endif
#ifdef VerticalWallsAdiabatic
        currentModel(10)=1
#endif
#ifdef VerticalWallsPeriodicalT
        currentModel(11)=1
#endif
#ifdef SideHeatedHa
        currentModel(12)=1
#endif
        ! 前13项为公共物理参数和输出间隔，第14、15项为Ha、phi；非磁场模式保持0。
        currentPhysics=[Rayleigh,Prandtl,Mach,Thot,Tcold,Snu,Sq,Qk,Qnu,paraA,outputSnapshotInterval,&
            reloadFileInterval,outputPltFileInterval,0.0d0,0.0d0]  !用数组构造物理信息
#ifdef SideHeatedHa
        currentPhysics(14)=Ha
        currentPhysics(15)=phi
#endif

        if(any(sig/=currentModel) .OR. any(phys/=currentPhysics)) &
            error stop 'Restart model/physics/output-cadence mismatch'

        read(m,'(a)') line !# itc_total NuReSampleCount t_ff errorU errorT
        read(m,*,iostat=ios) reloadItcOffset,NuReSampleCount,savedTime,errorU,errorT
        if(ios/=0) error stop 'Incomplete checkpoint counters'
        if(reloadItcOffset<0 .OR. mod(reloadItcOffset,refineRatio)/=0) &  !累计步数必须非负，而且处在粗细同步时刻
            error stop 'Restart is not at a synchronized time'
        !采样次数应该与累计步数对应
        if(NuReSampleCount/=reloadItcOffset/outputSnapshotIntervalItc) error stop 'Invalid restart sample count'
        !由步数重新计算时间

        reloadDimensionlessTime=dble(reloadItcOffset)/timeUnit
        !与文件中的 savedTime 比较，确认时间和步数一致
        if(.NOT.all(ieee_is_finite([savedTime,errorU,errorT]))) error stop 'Nonfinite restart metadata'
        if(abs(savedTime-reloadDimensionlessTime)>1d-12*max(1d0,abs(savedTime))) error stop 'Restart time mismatch'

#ifdef steadyFlow
        read(m,'(a)') line !# steadyCheckIntervalItc
        read(m,*,iostat=ios) savedCheckInterval
        if(ios/=0) error stop 'Missing steady check interval'
        if(savedCheckInterval/=steadyCheckIntervalItc) error stop 'Steady check interval mismatch'
#endif
        ! 各区几何只在meta中核对；bin只含数组，i最快、j其次、分布方向最后。
        read(m,'(a)') line  ! 粗网格：ni nj iStart iEnd jStart jEnd
        read(m,*,iostat=ios) geom
        if(ios/=0) error stop 'Incomplete coarse geometry'

        read(m,*,iostat=ios) coord ![xOffset, yOffset, dx, xmin, xmax, ymin, ymax]
        if(ios/=0) error stop 'Incomplete coarse coordinates'
        if(any(geom/=[size(f_coarse,1),size(f_coarse,2), &
            1+merge(interfaceSkin,0,refineRatio>1),size(f_coarse,1)-merge(interfaceSkin,0,refineRatio>1), &
            1+merge(interfaceSkin,0,refineRatio>1),size(f_coarse,2)-merge(interfaceSkin,0,refineRatio>1)]) .OR. &
            any(coord/=[xOffsetCoarse,yOffsetCoarse,dxCoarse, &
            centerIntegrationBox])) error stop 'Restart coarse layout mismatch'


        if(refineRatio>1) then

        read(m,'(a)') line  ! left: ni nj iStart iEnd jStart jEnd
        read(m,*,iostat=ios) geom
        if(ios/=0) error stop 'Incomplete left geometry'

        read(m,*,iostat=ios) coord  ![xOffset, yOffset, dx, xmin, xmax, ymin, ymax]
        if(ios/=0) error stop 'Incomplete left coordinates'
        if(any(geom/=[size(f_left,1),size(f_left,2), &
            1,size(f_left,1),1,size(f_left,2)]) .OR. &
            any(coord/=[xOffsetLeft,yOffsetLeft,dxFine, &
            xOffsetLeft,xOffsetLeft+size(f_left,1),yOffsetLeft,yOffsetLeft+size(f_left,2)])) error stop 'Restart left layout mismatch'

        read(m,'(a)') line  ! right: ni nj iStart iEnd jStart jEnd；随后为xOffset yOffset dx及积分外框。
        read(m,*,iostat=ios) geom
        if(ios/=0) error stop 'Incomplete right geometry'

        read(m,*,iostat=ios) coord  ![xOffset, yOffset, dx, xmin, xmax, ymin, ymax]
        if(ios/=0) error stop 'Incomplete right coordinates'
        if(any(geom/=[size(f_right,1),size(f_right,2), &
            1,size(f_right,1),1,size(f_right,2)]) .OR. &
            any(coord/=[xOffsetRight,yOffsetRight,dxFine, &
            xOffsetRight,xOffsetRight+size(f_right,1),yOffsetRight,yOffsetRight+size(f_right,2)])) error stop 'Restart right layout mismatch'
        
        read(m,'(a)') line  ! bottom: ni nj iStart iEnd jStart jEnd；随后为xOffset yOffset dx及积分外框。
        read(m,*,iostat=ios) geom 
        if(ios/=0) error stop 'Incomplete bottom geometry'
        read(m,*,iostat=ios) coord  ![xOffset, yOffset, dx, xmin, xmax, ymin, ymax]
        if(ios/=0) error stop 'Incomplete bottom coordinates'
        if(any(geom/=[size(f_bottom,1),size(f_bottom,2), &
            1,size(f_bottom,1),1,size(f_bottom,2)]) .OR. &
            any(coord/=[xOffsetBottom,yOffsetBottom,dxFine, &
            xOffsetBottom,xOffsetBottom+size(f_bottom,1),yOffsetBottom,yOffsetBottom+size(f_bottom,2)])) error stop 'Restart bottom layout mismatch'
        
        read(m,'(a)') line  ! top: ni nj iStart iEnd jStart jEnd；随后为xOffset yOffset dx及积分外框。
        read(m,*,iostat=ios) geom
        if(ios/=0) error stop 'Incomplete top geometry'
        read(m,*,iostat=ios) coord  ![xOffset, yOffset, dx, xmin, xmax, ymin, ymax]
        if(ios/=0) error stop 'Incomplete top coordinates'
        if(any(geom/=[size(f_top,1),size(f_top,2), &
            1,size(f_top,1),1,size(f_top,2)]) .OR. &
            any(coord/=[xOffsetTop,yOffsetTop,dxFine, &
            xOffsetTop,xOffsetTop+size(f_top,1),yOffsetTop,yOffsetTop+size(f_top,2)])) error stop 'Restart top layout mismatch'
        endif  ! 细网格下的网格判断
        close(m)
        reloadMetadataLoaded=.true.
    else     !meta文件不存在
        ! bin只存数组：没有meta时，直接使用手工累计步数，不从时间或文件名反推。
        name=trim(reloadFileName)
        if(len_trim(name)==0) error stop 'No latest.meta; set reloadFileName and reloadItcOffset' !文件名不为空
        if(reloadItcOffset<0 .OR. mod(reloadItcOffset,refineRatio)/=0) &   ! 累计步数必须非负，而且处在粗细同步时刻
            error stop 'Manual restart step must be nonnegative and synchronized'

        NuReSampleCount=reloadItcOffset/outputSnapshotIntervalItc !反算计数
        reloadDimensionlessTime=dble(reloadItcOffset)/timeUnit  !反算无量纲时间 tff
        ! 未读到meta，保留初始errorU/errorT；稳态旧场从bin恢复，下一次check重新计算误差。
    endif



    if(len_trim(name)<5) error stop 'Invalid restart filename'  !bin文件名长度检查
    if(name(len_trim(name)-3:len_trim(name))/='.bin') error stop 'Restart filename must end in .bin' !取文件名最后四个字符，确认它们是 .bin。
    
    ! 恢复完整场量与接口时间历史；手工分支要求使用同一网格、模型和输出间隔生成的bin。
    ! 'stream'	按连续数据流访问,'unformatted'	读取二进制数据
    open(newunit=k,file=trim(name),access='stream',form='unformatted',status='old',iostat=ios)
    if(ios/=0) error stop 'Cannot open checkpoint arrays'
    ! 粗网格数据
    read(k,iostat=ios) f_coarse,g_coarse,u_coarse,v_coarse,T_coarse,rho_coarse,Fx_coarse,Fy_coarse
    if(ios/=0) error stop 'Incomplete coarse state'
    ! 读取粗网格时间历史，保存的是接口插值所需的数据，不是三个完整粗网格场。
    read(k,iostat=ios) coarseHistoryPrevious,coarseHistoryCurrent,coarseHistoryNext
    if(ios/=0) error stop 'Incomplete coarse time history'

    ! 稳态下恢复上一次检查的旧场
#ifdef steadyFlow
    read(k,iostat=ios) up_coarse,vp_coarse,Tp_coarse
    if(ios/=0) error stop 'Incomplete coarse convergence reference'
#endif

    if(refineRatio>1) then  !依次恢复四个细区域
    read(k,iostat=ios) f_left,g_left,u_left,v_left,T_left,rho_left,Fx_left,Fy_left
    if(ios/=0) error stop 'Incomplete left state'

#ifdef steadyFlow
    read(k,iostat=ios) up_left,vp_left,Tp_left
    if(ios/=0) error stop 'Incomplete left convergence reference'
#endif

    read(k,iostat=ios) f_right,g_right,u_right,v_right,T_right,rho_right,Fx_right,Fy_right
    if(ios/=0) error stop 'Incomplete right state'

#ifdef steadyFlow
    read(k,iostat=ios) up_right,vp_right,Tp_right
    if(ios/=0) error stop 'Incomplete right convergence reference'
#endif

    read(k,iostat=ios) f_bottom,g_bottom,u_bottom,v_bottom,T_bottom,rho_bottom,Fx_bottom,Fy_bottom
    if(ios/=0) error stop 'Incomplete bottom state'

#ifdef steadyFlow
    read(k,iostat=ios) up_bottom,vp_bottom,Tp_bottom
    if(ios/=0) error stop 'Incomplete bottom convergence reference'
#endif

    read(k,iostat=ios) f_top,g_top,u_top,v_top,T_top,rho_top,Fx_top,Fy_top
    if(ios/=0) error stop 'Incomplete top state'

#ifdef steadyFlow
    read(k,iostat=ios) up_top,vp_top,Tp_top
    if(ios/=0) error stop 'Incomplete top convergence reference'
#endif

    endif
    close(k)
    !有 .meta 时，这个名字来自元数据；没有 .meta 时，它就是你手工填写的名字。
    reloadFileName=trim(name)
  end subroutine read_restart
!===================================================================================================

!===================================================================================================
! 子程序: truncate_output_files
! 作用: 续算回到某个检查点后，删除这个检查点之后的旧快照、plt 等输出，并同步整理输出索引
!===================================================================================================
  subroutine truncate_output_files()
    use commondata,only: reloadItcOffset,snapshotFilePrefix,pltFolderPrefix
    implicit none
    integer :: k, u, ios, n, i, step
    integer, allocatable :: steps(:)
    character(len=256), allocatable :: names(:)
    character(len=512) :: line
    character(len=256) :: name
    logical :: exists

    open(newunit=k,file='OutputFiles_2DOpenaccMultiblock.dat',status='old',iostat=ios) !打开统计清单
    if(ios/=0) error stop 'Missing output file index'

    n=0
    do
        read(k,'(a)',iostat=ios) line
        if(ios<0) exit !到文件末尾：退出循环。
        if(ios>0) error stop 'Output index read error' !读取出错：停止程序。
        if(len_trim(line)==0) cycle !空行：跳过，不计数。
        n=n+1   !计数有效行
    enddo

    allocate(steps(n),names(n)) !分配数组
    rewind(k) !把文件的读取位置移回开头

    i=0
    do while(i<n)
        read(k,'(a)',iostat=ios) line
        if(ios/=0) error stop 'Incomplete output index'
        if(len_trim(line)==0) cycle

        read(line,*,iostat=ios) step,name !在 line 里面阅读信息
        if(ios/=0 .OR. step<0) error stop 'Invalid output index row'

        ! index(字符串,要查找的内容) 的作用是查找子字符串：
        !- 找到：返回第一次出现的位置，从 1 开始计数。
        !- 找不到：返回 0。
        !achar(整数) 则把 ASCII 编码转换成对应字符,achar(92) 指的是\，反斜杠 
        if(index(name,'/')>0 .OR. index(name,achar(92))>0 .OR. index(name,':')>0 .OR. &
            index(name,'..')>0) error stop 'Output index contains a path'  !如果文件名包含 /、\、: 或 .. 中的任何一种，就报错。
        !不是以快照前缀开头；不是以 plt 前缀开头；也不是最终 Nu/Re 结果文件，那就报错
        if(index(name,snapshotFilePrefix//'-')/=1 .AND. index(name,pltFolderPrefix//'-')/=1 .AND. &
            trim(name)/='FinalNuRe_2DOpenaccMultiblock.dat') error stop 'Unexpected output filename' 
        i=i+1
        steps(i)=step
        names(i)=name
    enddo
    close(k)

    ! 前缀已完整读入，删除检查点之后的输出，不保存备份。
    do i=1,n
        if(steps(i)<=reloadItcOffset) cycle  !直接跳过，文件保留
        inquire(file=trim(names(i)),exist=exists) !文件是否存在
        if(exists) then
            open(newunit=u,file=trim(names(i)),status='old',iostat=ios)
            if(ios/=0) error stop 'Cannot open stale output for deletion' !打开失败就停止，手动解决
            close(u,status='delete',iostat=ios) !关闭并删除这个文件
            if(ios/=0) error stop 'Cannot delete stale output' !检查能否关闭并删除
        endif
    enddo

    open(newunit=k,file='OutputFiles_2DOpenaccMultiblock.dat',status='replace') !'replace' 覆盖原内容
    do i=1,n
        if(steps(i)<=reloadItcOffset) write(k,'(I0,1X,A)') steps(i),trim(names(i)) !只把需要保留的记录写回
    enddo
    close(k)
    deallocate(steps,names)
  end subroutine truncate_output_files
!===================================================================================================

!===================================================================================================
! 子程序: enter_data_2d_openacc
! 作用: 将主机上初始化或恢复的数据复制到设备，建立常驻数据区。
!===================================================================================================
  subroutine enter_data_2d_openacc()
    use commondata
    implicit none

    ! 将主机上已经初始化或恢复的数据复制到设备，建立常驻存储。
    !$acc enter data copyin(ex, ey, omega, omegaT)
    !$acc enter data copyin(f_coarse, g_coarse, f_post_coarse, g_post_coarse, rho_coarse, u_coarse, &
    !$acc& v_coarse, T_coarse, Fx_coarse, Fy_coarse, coarseHistoryPrevious,coarseHistoryCurrent,coarseHistoryNext)
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

  end subroutine enter_data_2d_openacc
!===================================================================================================

!===================================================================================================
! 子程序: save_coarse_history
! 作用: 从当前粗场完成空间插值，分别保存四区所需的时间历史，各区编号都从 1 开始。
!===================================================================================================
  subroutine save_coarse_history(history)
    use commondata
    implicit none
    real(kind=8), intent(inout) :: history(size(coarseHistoryCurrent,1),size(coarseHistoryCurrent,2))

    if(refineRatio==1) return  !直接返回

    ! 接收范围由分区设置和 interfaceSkin 得到，传入的是各细区的局部下标范围。
    ! 默认局部下标：左i=129:130,j=129:895；右i=1:2,j=129:895；下i=1:763,j=129:130；上i=1:763,j=1:2。
    call save_coarse_history_grid(history,1,fineLayerCellsLeft+1,fineLayerCellsLeft+interfaceSkin, &
        fineLayerCellsBottom+1,ny-fineLayerCellsTop,xOffsetLeft,yOffsetLeft)
    call save_coarse_history_grid(history,2,1,interfaceSkin, &
        fineLayerCellsBottom+1,ny-fineLayerCellsTop,xOffsetRight,yOffsetRight)
    call save_coarse_history_grid(history,3,1,size(f_bottom,1), &
        fineLayerCellsBottom+1,fineLayerCellsBottom+interfaceSkin,xOffsetBottom,yOffsetBottom)
    call save_coarse_history_grid(history,4,1,size(f_top,1), &
        1,interfaceSkin,xOffsetTop,yOffsetTop)
  end subroutine save_coarse_history
!===================================================================================================

!===================================================================================================
! 子程序: save_coarse_history_grid
! 作用: 直接从粗场取共址值或空间插值，只保存给定细接收矩形所需的20个交换量。
!===================================================================================================
  subroutine save_coarse_history_grid(history,region,iStart,iEnd,jStart,jEnd,xOffset,yOffset)
    use commondata
    implicit none
    real(kind=8), intent(inout) :: history(size(coarseHistoryCurrent,1),size(coarseHistoryCurrent,2))
    integer(kind=4), intent(in) :: region,iStart,iEnd,jStart,jEnd
    real(kind=8), intent(in) :: xOffset,yOffset
    integer(kind=4) :: i,j,a,ix,iy,ii,jj,sourceX,sourceY,width,p,ni,nj
    real(kind=8) :: qx,qy,wx(4),wy(4),fv(0:8),gv(0:4),buffer(0:19),value(0:19)

    ni=size(f_coarse,1)  !ni=389
    nj=size(f_coarse,2)  !nj=389
    !$acc parallel loop collapse(2) present(f_coarse,g_coarse,rho_coarse,u_coarse,v_coarse, &
    !$acc& T_coarse,Fx_coarse,Fy_coarse,history) async(1) &
    !$acc& firstprivate(region,iStart,iEnd,jStart,xOffset,yOffset,ni,nj, &
    !$acc& xOffsetCoarse,yOffsetCoarse) &
    !$acc& private(a,ix,iy,ii,jj,sourceX,sourceY,width,p,qx,qy,wx,wy,fv,gv,buffer,value)
    do j=jStart,jEnd
        do i=iStart,iEnd
            qx=(xOffset+(i-0.5d0)-xOffsetCoarse)/dxCoarse+0.5d0   !细区局部下标 i -> 全局物理坐标 x -> 粗网格连续下标 qx
            qy=(yOffset+(j-0.5d0)-yOffsetCoarse)/dxCoarse+0.5d0   !细区局部下标 j -> 全局物理坐标 y -> 粗网格连续下标 qy
            !确定取哪个粗节点，或者用哪些粗节点插值，还有权重
            call coarse_source_stencil(qx,qy,ni,nj,interfaceSkin,sourceX,sourceY,width,wx,wy)
            buffer=0.0d0         !清空当前目标节点的累加数组
            do iy=1,width        !width指的是插值节点数
                jj=sourceY+iy-1  !来源粗数组的局部下标；例如 sourceY=3 时 jj=3,4,5,6
                do ix=1,width
                    ii=sourceX+ix-1
                    do a=0,8     !对每个来源粗节点，先读出分布函数
                        fv(a)=f_coarse(ii,jj,a)
                    enddo
                    do a=0,4
                        gv(a)=g_coarse(ii,jj,a)
                    enddo
                    call pack_interface_values(dxCoarse,SnuCoarse,SqCoarse,QkCoarse,QnuCoarse,fv,gv, &
                        rho_coarse(ii,jj),u_coarse(ii,jj),v_coarse(ii,jj),T_coarse(ii,jj), &
                        Fx_coarse(ii,jj),Fy_coarse(ii,jj),value)  !将这个来源粗节点整理成 20 个交换量
                    do a=0,19
                        buffer(a)=buffer(a)+wx(ix)*wy(iy)*value(a) !四点拉格朗日插值在xy方向
                    enddo
                enddo
            enddo  !当前目标细接收节点的空间插值结束
            ! 每个节点连续放 20 个分量；区内 i 先变，j 后变。
            p=20*((j-jStart)*(iEnd-iStart+1)+i-iStart) !当前节点前面的节点总数，再乘以每个节点的 20 个交换量
            do a=0,19
                history(p+a+1,region)=buffer(a)
            enddo
        enddo
    enddo
  end subroutine save_coarse_history_grid
!===================================================================================================

!===================================================================================================
! 子程序: coarse_source_stencil
! 作用: 由目标的粗网格连续下标 qx/qy，找到共址粗节点或四乘四空间插值模板。
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

    w=1.0d0   !下面通过连续相乘来计算权重，所以初值必须是 1,如果初值为 0，结果都会是 0
    do a=0,3
        do k=0,3
            if(a/=k) w(a+1)=w(a+1)*(q-dble(k))/dble(a-k)
        enddo
    enddo
  end subroutine lagrange_weights
!===================================================================================================

!===================================================================================================
! 子程序: pack_interface_values
! 作用: 把一个来源节点的宏观量、单位时间力和缩放非平衡矩放入局部缓冲数组。
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
    ! 当前refineRatio=2：左细接收节点(129,129)使用粗模板(3:6,3:6)，第一次传入fv/gv为粗节点(3,3)的数据。
    ! 此次dt=dxCoarse=2，Snu/Sq/Qk/Qnu均为对应粗松弛率；细到粗调用则传入细网格参数。

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
    ! dt显式传入：粗块dt=dxCoarse，细块dt=dxFine=1。
    fm=[0.0d0,6.0d0*(ux*fx+uy*fy),-6.0d0*(ux*fx+uy*fy),fx,-fx,fy,-fy,&
        2.0d0*(ux*fx-uy*fy),ux*fy+uy*fx]

    ! 5. 各矩对应的松弛率；密度、动量和温度守恒矩的松弛率置零，由单独保存的宏观量重建。
    s=[0.0d0,Snu,Snu,0.0d0,Sq,0.0d0,Sq,Snu,Snu]
    q=[0.0d0,Qk,Qk,Qnu,Qnu]

    ! 6. 前6项依次保存rho、u、v、T、F_x、F_y；所有网格直接保存单位时间的力。
    buffer(0:5)=[rh,ux,uy,temp,fx,fy]

    ! 7. buffer(6+a)=(s_a/dt)*(m_a-meq_a+dt*(F_m,a)/2)，dt为来源网格时间步。
    ! 碰撞源项为dt*(1-s_a/2)*F_m,a。
    ! 本例粗块dt=2，括号内为m_a-meq_a+(F_m,a)；细块dt=1时为m_a-meq_a+(F_m,a)/2。
    ! 粗来源的缩放系数为s(a)/2；buffer(6)、buffer(9)、buffer(11)对应守恒矩，有限场下为零。
    do a=0,8
        buffer(6+a)=s(a)/dt*(m(a)-meq(a)+0.5d0*dt*fm(a))
    enddo

    ! 8. buffer(15:19)保存5个缩放后的温度非平衡矩；粗来源缩放系数为q(a)/2。
    ! buffer(15)对应温度守恒矩，有限场下为零；返回后外层再对这20项做空间插值。
    do a=0,4
        buffer(15+a)=q(a)/dt*(n(a)-neq(a))
    enddo
  end subroutine pack_interface_values
!===================================================================================================

!===================================================================================================
! 子程序: fine_to_coarse
! 作用: 直接读取四个细区的当前场，重建粗网格接收层。
!===================================================================================================
  subroutine fine_to_coarse()
    use commondata
    implicit none

    ! 接收层始终保存完整 f/g；初始化及每轮同步后重建，下一轮粗网格直接使用。
    ! 仅在粗细同步时传输：左右各 interfaceSkin 列贯穿全高，上下行去掉已由左右覆盖的角点。
    ! 默认粗区局部下标：左i=1:2，右i=388:389；下j=1:2，上j=388:389，二者i=3:387。
    call copy_fine_to_coarse(1,interfaceSkin,1,size(f_coarse,2))
    call copy_fine_to_coarse(size(f_coarse,1)-interfaceSkin+1,size(f_coarse,1),1,size(f_coarse,2))
    call copy_fine_to_coarse(interfaceSkin+1,size(f_coarse,1)-interfaceSkin,1,interfaceSkin)
    call copy_fine_to_coarse(interfaceSkin+1,size(f_coarse,1)-interfaceSkin, &
        size(f_coarse,2)-interfaceSkin+1,size(f_coarse,2))
  end subroutine fine_to_coarse
!===================================================================================================

!===================================================================================================
! 子程序: copy_fine_to_coarse
! 作用: 粗接收点与细节点共址，直接从对应细节点构造局部缓冲并写入粗数组。
!===================================================================================================
  subroutine copy_fine_to_coarse(iStart,iEnd,jStart,jEnd)  !接收方——粗区局部下标范围
    use commondata
    implicit none
    integer(kind=4), intent(in) :: iStart,iEnd,jStart,jEnd
    integer(kind=4) :: i,j,a,sourceX,sourceY
    real(kind=8) :: fv(0:8),gv(0:4),buffer(0:19),rh,ux,uy,temp,fx,fy

    ! 由粗区局部下标得到全局细网格下标；默认粗(1,1) -> 全局物理坐标(123.5,123.5) -> (I,J)=(124,124)。
    !$acc parallel loop collapse(2) present(f_coarse,g_coarse,rho_coarse,u_coarse,v_coarse,T_coarse, &
    !$acc& Fx_coarse,Fy_coarse, &
    !$acc& f_left,g_left,rho_left,u_left,v_left,T_left,Fx_left,Fy_left, &
    !$acc& f_right,g_right,rho_right,u_right,v_right,T_right,Fx_right,Fy_right, &
    !$acc& f_bottom,g_bottom,rho_bottom,u_bottom,v_bottom,T_bottom,Fx_bottom,Fy_bottom, &
    !$acc& f_top,g_top,rho_top,u_top,v_top,T_top,Fx_top,Fy_top) async(1) &
    !$acc& private(a,sourceX,sourceY,fv,gv,buffer,rh,ux,uy,temp,fx,fy)
    do j=jStart,jEnd      !粗区局部下标
        do i=iStart,iEnd  !粗区局部下标(1,1) ← 全局细网格下标(124,124)
            sourceX=fineLayerCellsLeft+(i-1-interfaceSkin)*refineRatio  !全局细网格下标
            sourceY=fineLayerCellsBottom+(j-1-interfaceSkin)*refineRatio !全局细网格下标
            if(sourceX<=fineLayerCellsLeft+interfaceSkin) then  !判断这个细节点存在哪个细区数组里
                do a=0,8                                      !从左细网格读取
                    fv(a)=f_left(sourceX,sourceY,a)
                enddo
                do a=0,4
                    gv(a)=g_left(sourceX,sourceY,a)
                enddo
                rh=rho_left(sourceX,sourceY)
                ux=u_left(sourceX,sourceY)
                uy=v_left(sourceX,sourceY)
                temp=T_left(sourceX,sourceY)
                fx=Fx_left(sourceX,sourceY)
                fy=Fy_left(sourceX,sourceY)      !来源细节点的数据
            else if(sourceX>nx-fineLayerCellsRight-interfaceSkin) then  !右细区
                sourceX=sourceX-(nx-fineLayerCellsRight-interfaceSkin)  !全局细网格下标换成右细区局部下标
                do a=0,8
                    fv(a)=f_right(sourceX,sourceY,a)
                enddo
                do a=0,4
                    gv(a)=g_right(sourceX,sourceY,a)
                enddo
                rh=rho_right(sourceX,sourceY)
                ux=u_right(sourceX,sourceY)
                uy=v_right(sourceX,sourceY)
                temp=T_right(sourceX,sourceY)
                fx=Fx_right(sourceX,sourceY)
                fy=Fy_right(sourceX,sourceY)
            else if(sourceY<=fineLayerCellsBottom+interfaceSkin) then  ! 下细区,前两项不满足，且 sourceY<=130
                sourceX=sourceX-(fineLayerCellsLeft+interfaceSkin)     ! 全局细网格下标换成下细区局部下标
                do a=0,8
                    fv(a)=f_bottom(sourceX,sourceY,a)
                enddo
                do a=0,4
                    gv(a)=g_bottom(sourceX,sourceY,a)
                enddo
                rh=rho_bottom(sourceX,sourceY)
                ux=u_bottom(sourceX,sourceY)
                uy=v_bottom(sourceX,sourceY)
                temp=T_bottom(sourceX,sourceY)
                fx=Fx_bottom(sourceX,sourceY)
                fy=Fy_bottom(sourceX,sourceY)
            else                                                     !上细区
                sourceX=sourceX-(fineLayerCellsLeft+interfaceSkin)   !全局细网格下标换成上细区局部下标
                sourceY=sourceY-(ny-fineLayerCellsTop-interfaceSkin) !全局细网格下标换成上细区局部下标
                do a=0,8
                    fv(a)=f_top(sourceX,sourceY,a)
                enddo
                do a=0,4
                    gv(a)=g_top(sourceX,sourceY,a)
                enddo
                rh=rho_top(sourceX,sourceY)
                ux=u_top(sourceX,sourceY)
                uy=v_top(sourceX,sourceY)
                temp=T_top(sourceX,sourceY)
                fx=Fx_top(sourceX,sourceY)
                fy=Fy_top(sourceX,sourceY)
            endif
            !传入细网格的时间步和松弛率,用细网格参数打包
            call pack_interface_values(1.0d0,Snu,Sq,Qk,Qnu,fv,gv,rh,ux,uy,temp,fx,fy,buffer)
            !恢复宏观量及单位时间力，直接写入当前粗接收节点
            !用粗网格参数重建矩；推进前保留矩，初始化及同步后才还原分布函数
            call unpack_interface_values(dxCoarse,SnuCoarse,SqCoarse,QkCoarse,QnuCoarse,buffer,fv,gv, &
                rho_coarse(i,j),u_coarse(i,j),v_coarse(i,j),T_coarse(i,j),Fx_coarse(i,j),Fy_coarse(i,j))
            do a=0,8       !将换算后的流场分布 f 写入粗接收层
                f_coarse(i,j,a)=fv(a)
            enddo
            do a=0,4
                g_coarse(i,j,a)=gv(a)
            enddo
        enddo
    enddo
  end subroutine copy_fine_to_coarse
!===================================================================================================

!===================================================================================================
! 子程序: unpack_interface_values
! 作用: 用接收网格的松弛率和dt恢复完整矩及宏观量；推进前保留矩，同步点还原f/g。
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

    ! 1. 恢复接收节点的宏观量；本次细到粗调用写入 rho_coarse(i,j) 等数组。
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

    ! 4. 用接收网格的松弛率和时间步恢复完整矩；细到粗时 dt=dxCoarse=2。
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

!===================================================================================================
! 子程序: coarse_to_fine
! 作用: 四个细区按打包顺序读取粗时间层，结果写入各自的接收节点。
!===================================================================================================
  subroutine coarse_to_fine(timeInterpWeight)
    use commondata
    implicit none
    real(kind=8), intent(in) :: timeInterpWeight(0:2)

    ! 初始化、细子步开始及同步点均重建接收层 f/g，数组内容始终为分布函数。
    ! 读取时间历史 → 时间取值或插值 → 写入细接收层
    ! 读取历史数组第 1 列，
    ! 给左细区局部下标 i=129:130、j=129:895 的接收节点恢复数据。
    call interpolate_coarse_to_fine(1,size(f_left,1),size(f_left,2), &
        fineLayerCellsLeft+1,fineLayerCellsLeft+interfaceSkin,fineLayerCellsBottom+1,ny-fineLayerCellsTop, &
        f_left,g_left,rho_left,u_left,v_left,T_left,Fx_left,Fy_left,timeInterpWeight)
    call interpolate_coarse_to_fine(2,size(f_right,1),size(f_right,2), &
        1,interfaceSkin,fineLayerCellsBottom+1,ny-fineLayerCellsTop, &
        f_right,g_right,rho_right,u_right,v_right,T_right,Fx_right,Fy_right,timeInterpWeight)
    call interpolate_coarse_to_fine(3,size(f_bottom,1),size(f_bottom,2), &
        1,size(f_bottom,1),fineLayerCellsBottom+1,fineLayerCellsBottom+interfaceSkin, &
        f_bottom,g_bottom,rho_bottom,u_bottom,v_bottom,T_bottom,Fx_bottom,Fy_bottom,timeInterpWeight)
    call interpolate_coarse_to_fine(4,size(f_top,1),size(f_top,2), &
        1,size(f_top,1),1,interfaceSkin, &
        f_top,g_top,rho_top,u_top,v_top,T_top,Fx_top,Fy_top,timeInterpWeight)
  end subroutine coarse_to_fine
!===================================================================================================

!===================================================================================================
! 子程序: interpolate_coarse_to_fine
! 作用: 对指定接收矩形做时间插值，直接转换并写入细网格。
!===================================================================================================
  subroutine interpolate_coarse_to_fine(region,ni,nj,iStart,iEnd,jStart,jEnd, &
      f,g,rho,u,v,T,Fx,Fy,timeInterpWeight)
    use commondata,only: coarseHistoryPrevious,coarseHistoryCurrent,coarseHistoryNext,Snu,Sq,Qk,Qnu
    implicit none
    integer(kind=4), intent(in) :: region,ni,nj,iStart,iEnd,jStart,jEnd
    real(kind=8), intent(in) :: timeInterpWeight(0:2)
    real(kind=8), intent(inout) :: f(ni,nj,0:8),g(ni,nj,0:4),rho(ni,nj),u(ni,nj),v(ni,nj)
    real(kind=8), intent(inout) :: T(ni,nj),Fx(ni,nj),Fy(ni,nj)
    integer(kind=4) :: i,j,a,p
    real(kind=8) :: buffer(0:19),fv(0:8),gv(0:4)
    real(kind=8) :: timeWeight0,timeWeight1,timeWeight2

    ! 保留 NVHPC 所需的标量捕获，不能改为 firstprivate(timeInterpWeight) 数组形参。
    timeWeight0=timeInterpWeight(0)
    timeWeight1=timeInterpWeight(1)
    timeWeight2=timeInterpWeight(2)
    ! 时间插值使用局部 buffer；按细网格参数恢复矩后转成 f/g，写入接收层。
    !$acc parallel loop collapse(2) present(coarseHistoryPrevious,coarseHistoryCurrent,coarseHistoryNext, &
    !$acc& f,g,rho,u,v,T,Fx,Fy) async(1) &
    !$acc& firstprivate(region,iStart,iEnd,jStart,timeWeight0,timeWeight1,timeWeight2) private(a,p,buffer,fv,gv)
    do j=jStart,jEnd      !找到当前细节点在历史数组中的位置
        do i=iStart,iEnd  !被插值的细网格范围
            p=20*((j-jStart)*(iEnd-iStart+1)+i-iStart)
            do a=0,19    !对当前节点的 20 个交换量做时间取值或插值
                buffer(a)=timeWeight0*coarseHistoryPrevious(p+a+1,region)+ &
                    timeWeight1*coarseHistoryCurrent(p+a+1,region)+timeWeight2*coarseHistoryNext(p+a+1,region)
            enddo
            call unpack_interface_values(1.0d0,Snu,Sq,Qk,Qnu,buffer,fv,gv, &  !接收方是细网格 fg
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
! 子程序: output_time_string
! 作用: 将实际累计 t_ff 写成文件名中的时间，去掉小数末尾多余的零。
!===================================================================================================
  subroutine output_time_string(timeString)

    use commondata
    implicit none

    character(len=40), intent(out) :: timeString
    integer(kind=4) :: i

    ! t_ff保留十二位小数供查看；文件名另含累计整数步数，不同保存步不会因时间舍入而重名。
    write(timeString,'(F40.12)') dble(reloadItcOffset+itc)/timeUnit  !累计 tff
    if(index(timeString,'*')/=0) error stop 'tff exceeds filename field width' !检查数字是否超过可容纳的宽度。
    timeString=adjustl(timeString) !adjustl 将字符串左对齐，前面的空格移到末尾
    i=len_trim(timeString) !len_trim 返回忽略末尾空格后的长度，因而 i 指向最后一个有效字符。
    do while(timeString(i:i)=='0')
        timeString(i:i)=' '  !如果它是 '0'，就改成空格，再向前检查一个字符
        i=i-1
    enddo
    if(timeString(i:i)=='.') timeString(i:i)=' ' !如果只剩一个小数点，也去掉它
  end subroutine output_time_string
!===================================================================================================

!===================================================================================================
! 子程序: advance_coarse
! 作用: 使用中心粗网格数组及粗网格参数，推进一个粗时间步。
!===================================================================================================
  subroutine advance_coarse()
    use commondata
    implicit none

    ! 粗接收层 f/g 已由初始化、续算读取或上一轮同步准备好，此处直接推进，不再重复交换。
    ! 单网格时同样使用 f/g 推进，所有节点都在碰撞内部转换为矩。
    call advance_grid(size(f_coarse,1),size(f_coarse,2),dxCoarse,SnuCoarse,SqCoarse,QkCoarse,QnuCoarse,&
        gBeta,coarseIsWall,f_coarse,g_coarse,f_post_coarse,g_post_coarse,rho_coarse,&
        u_coarse,v_coarse,T_coarse,Fx_coarse,Fy_coarse)  !把中心粗网格推进一个粗时间步。
  end subroutine advance_coarse
!===================================================================================================

!===================================================================================================
! 子程序: advance_grid
! 作用: 按原均匀网格算法次序，推进传入数组的流场和温度场。
!===================================================================================================
  subroutine advance_grid(ni,nj,dt,Snu,Sq,Qk,Qnu,gBeta,wall,f,g,f_post,g_post,rho,u,v,T,Fx,Fy)

#ifdef SideHeatedHa
    use commondata,only: B2sigemarho
#endif
    implicit none

    ! 先推进流场，再用更新后的速度推进温度场；数组和松弛率由调用处传入。
    ! wall(1:4) 依次表示左、右、下、上是否为真实墙面；粗细人工接口不能当墙面反弹。

    integer(kind=4) :: ni, nj, j, a
    real(kind=8) :: dt, Snu, Sq, Qk, Qnu, gBeta
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

    ! dt显式传入：粗块dt=dxCoarse，细块dt=dxFine=1。
    ! gBeta和磁场系数不预乘dt；Fx/Fy保存F，碰撞与速度恢复各自显式乘dt。
    ! 普通节点和粗接收层均保存 f/g，碰撞统一执行 f->m 和 g->n。
    ! 推进后，最外两层等待下次接口重建，不能作为反向插值的数据来源。
    call collision(ni,nj,dt,f,f_post,rho,u,v,Fx,Fy,T,Snu,Sq,&
        gBeta &
#ifdef SideHeatedHa
        ,B2sigemarho &
#endif
        )

#ifdef VerticalWallsPeriodicalU
    ! 单网格自己左右连接；多块粗区没有物理外边界。
    if(wall(1) .AND. wall(2)) then
        !$acc parallel loop collapse(2) present(f_post) async(1)
        do a=0,8
            do j=1,nj
                f_post(0,j,a)=f_post(ni,j,a)
                f_post(ni+1,j,a)=f_post(1,j,a)
            enddo
        enddo
    endif
#endif
    call streaming(ni,nj,f,f_post)

    call bounceback(ni,nj,f,f_post,wall(1),wall(2),&
        wall(3),wall(4))

    call macro(ni,nj,dt,f,rho,u,v,Fx,Fy)

    call collisionT(ni,nj,g,g_post,u,v,T,Qk,Qnu)

#ifdef VerticalWallsPeriodicalT
    ! 单网格自己左右连接；多块粗区没有物理外边界。
    if(wall(1) .AND. wall(2)) then
        !$acc parallel loop collapse(2) present(g_post) async(1)
        do a=0,4
            do j=1,nj
                g_post(0,j,a)=g_post(ni,j,a)
                g_post(ni+1,j,a)=g_post(1,j,a)
            enddo
        enddo
    endif
#endif
    call streamingT(ni,nj,g,g_post)

    call bouncebackT(ni,nj,g,g_post,wall(1),wall(2),&
        wall(3),wall(4))

    call macroT(ni,nj,g,T)
  end subroutine advance_grid
!===================================================================================================

!===================================================================================================
! 子程序: collision
! 作用: 执行 D2Q9 流场矩空间碰撞，生成碰撞后分布函数。
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
    if(reloadItcOffset+itc==0) then  !累计步数为零，计算刚开始，没有真实的负时间历史层。
        timeInterpWeight=[0.0d0,1.0d0-theta,theta]    ! 缺少负时间层时线性启动
    else   !三个等间隔时间点的二次拉格朗日插值权重
        timeInterpWeight=[0.5d0*theta*(theta-1.0d0),1.0d0-theta**2,0.5d0*theta*(theta+1.0d0)]
    endif
  end subroutine coarse_time_weights
!===================================================================================================

!===================================================================================================
! 子程序: advance_fine
! 作用: 分阶段推进四个细区域，并在迁移前交换同级相邻节点数据。
!===================================================================================================
  subroutine advance_fine()
    use commondata
    implicit none

    ! 调用前 coarse_to_fine 已按当前子步时刻更新细接收层的 f/g。
    ! 所有节点统一由 f/g 求矩，再进行矩空间碰撞；每个节点只碰撞一次。
    ! 先完成流场推进并更新速度，再由原有 g 求温度矩并进行温度碰撞。
    ! 四个细区先全部碰撞，再交换相邻区的碰撞后分布，最后迁移。
    ! 流场和温度场各有一次同级交换，保证温度碰撞使用已经更新的速度。
    ! 细区拼接线只改变存储位置，不引入物理边界、空间插值或时间插值。

    ! 四个细区先全部碰撞，再交换迁移外圈；避免读取另一细区尚未更新的数据。
    call collision(size(f_left,1),size(f_left,2),dxFine,f_left,f_post_left,rho_left,u_left,v_left,Fx_left,Fy_left,&
        T_left,Snu,Sq,gBeta &
#ifdef SideHeatedHa
        ,B2sigemarho &
#endif
        )

    call collision(size(f_right,1),size(f_right,2),dxFine,f_right,f_post_right,rho_right,u_right,v_right,Fx_right,&
        Fy_right,T_right,Snu,Sq,gBeta &
#ifdef SideHeatedHa
        ,B2sigemarho &
#endif
        )

    call collision(size(f_bottom,1),size(f_bottom,2),dxFine,f_bottom,f_post_bottom,rho_bottom,u_bottom,v_bottom,&
        Fx_bottom,Fy_bottom,T_bottom,Snu,Sq,gBeta &
#ifdef SideHeatedHa
        ,B2sigemarho &
#endif
        )

    call collision(size(f_top,1),size(f_top,2),dxFine,f_top,f_post_top,rho_top,u_top,v_top,Fx_top,Fy_top,T_top,&
        Snu,Sq,gBeta &
#ifdef SideHeatedHa
        ,B2sigemarho &
#endif
        )

    ! 四区流场碰撞完成后，复制相邻细区的 f_post，再开始迁移。
    call exchange_fine_flow_halos()

    call streaming(size(f_left,1),size(f_left,2),f_left,f_post_left)

    call bounceback(size(f_left,1),size(f_left,2),f_left,f_post_left,leftIsWall(1),leftIsWall(2),leftIsWall(3),leftIsWall(4))

    call macro(size(f_left,1),size(f_left,2),dxFine,f_left,rho_left,u_left,v_left,Fx_left,Fy_left)

    call streaming(size(f_right,1),size(f_right,2),f_right,f_post_right)

    call bounceback(size(f_right,1),size(f_right,2),f_right,f_post_right,rightIsWall(1),rightIsWall(2),rightIsWall(3),&
        rightIsWall(4))

    call macro(size(f_right,1),size(f_right,2),dxFine,f_right,rho_right,u_right,v_right,Fx_right,Fy_right)

    call streaming(size(f_bottom,1),size(f_bottom,2),f_bottom,f_post_bottom)

    call bounceback(size(f_bottom,1),size(f_bottom,2),f_bottom,f_post_bottom,bottomIsWall(1),bottomIsWall(2),&
        bottomIsWall(3),bottomIsWall(4))

    call macro(size(f_bottom,1),size(f_bottom,2),dxFine,f_bottom,rho_bottom,u_bottom,v_bottom,Fx_bottom,Fy_bottom)

    call streaming(size(f_top,1),size(f_top,2),f_top,f_post_top)

    call bounceback(size(f_top,1),size(f_top,2),f_top,f_post_top,topIsWall(1),topIsWall(2),topIsWall(3),topIsWall(4))

    call macro(size(f_top,1),size(f_top,2),dxFine,f_top,rho_top,u_top,v_top,Fx_top,Fy_top)

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
    integer(kind=4) :: j, a
    !参数中先出现的是接收方，后出现的是来源方,nint表示把实数偏移量取最近整数，交换的是同一级细网格之间的碰撞后分布函数
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
#ifdef VerticalWallsPeriodicalU
    ! 左右细区贯穿全高；只交换外侧迁移列，不涉及中心粗细接收层。
    ! 斜方向的j偏移由streaming完成；越过上下墙面的方向随后由壁面条件覆盖。
    !$acc parallel loop collapse(2) present(f_post_left,f_post_right) async(1)
    do a=0,8
        do j=1,ny
            f_post_left(0,j,a)=f_post_right(size(f_post_right,1)-2,j,a)
            f_post_right(size(f_post_right,1)-1,j,a)=f_post_left(1,j,a)
        enddo
    enddo
#endif
  end subroutine exchange_fine_flow_halos
!===================================================================================================

!===================================================================================================
! 子程序: exchange_fine_thermal_halos
! 作用: 在左右细区与上下细区之间双向复制温度碰撞后的迁移外圈。
!===================================================================================================
  subroutine exchange_fine_thermal_halos()
    use commondata
    implicit none
    integer(kind=4) :: j, a

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
#ifdef VerticalWallsPeriodicalT
    ! 左右细区贯穿全高；只交换外侧迁移列，不涉及中心粗细接收层。
    ! 斜方向的j偏移由streaming完成；越过上下墙面的方向随后由壁面条件覆盖。
    !$acc parallel loop collapse(2) present(g_post_left,g_post_right) async(1)
    do a=0,4
        do j=1,ny
            g_post_left(0,j,a)=g_post_right(size(g_post_right,1)-2,j,a)
            g_post_right(size(g_post_right,1)-1,j,a)=g_post_left(1,j,a)
        enddo
    enddo
#endif
  end subroutine exchange_fine_thermal_halos
!===================================================================================================

!===================================================================================================
! 子程序: copy_fine_halo
! 作用: 根据坐标偏移，复制相邻细区域内部与接收区域迁移外圈的交集。
!===================================================================================================
  subroutine copy_fine_halo(nq,ni,nj,targetIndexX,targetIndexY,receiver,mi,mj,sourceIndexX,sourceIndexY,donor)
    implicit none

    ! nq 为方向数：流场 9，温度场 5；ni/nj 属于接收区，mi/mj 属于来源区。
    ! targetIndexX/targetIndexY、sourceIndexX/sourceIndexY 是两细区的整数下标偏移；dxFine=1 时与坐标偏移数值相同。
    ! 仅复制 donor 的内部节点；receiver 的交集在 0 或 n+1 外圈，包含对角迁移需要的位置。

    integer(kind=4), intent(in) :: nq, ni, nj, targetIndexX, targetIndexY, mi, mj, sourceIndexX, sourceIndexY
    real(kind=8) :: receiver(0:ni+1,0:nj+1,0:nq-1), donor(0:mi+1,0:mj+1,0:nq-1)
    integer(kind=4) :: i, j, a

    ! 全局细网格下标 = 整数下标偏移 + 局部下标。只取 receiver 外圈与 donor 内部的交集。

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
! 子程序: rotate_coarse_history
! 作用: 在粗细同步后滚动粗历史，准备下一轮时间插值。
!===================================================================================================
  subroutine rotate_coarse_history()
    use commondata
    implicit none
    integer(kind=4) :: p,region

    !$acc parallel loop collapse(2) present(coarseHistoryPrevious,coarseHistoryCurrent,coarseHistoryNext) async(1)
    do region=1,4
        do p=1,size(coarseHistoryCurrent,1)
            coarseHistoryPrevious(p,region)=coarseHistoryCurrent(p,region)
            coarseHistoryCurrent(p,region)=coarseHistoryNext(p,region)
        enddo
    enddo
  end subroutine rotate_coarse_history
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

    du=0.0d0 !速度变化的平方积分
    uu=0.0d0 !当前速度的平方积分
    dt=0.0d0 !温度变化绝对值的积分
    tt=0.0d0 !当前温度绝对值的积分

    call update_host_all(.false.) !传输数据到 CPU

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
    ! 温度沿用基准代码的相对 L1 判据；分母接近零时使用绝对变化量。
    if(tt .GT. 1.0d-30) then
        errorT=dt/tt
    else
        errorT=dt
    endif
    open(newunit=k,file='Convergence_2DOpenaccMultiblock.dat',status='unknown',position='append')
    write(k,'(I0,1X,3(ES24.16E3,1X))') reloadItcOffset+itc,dble(reloadItcOffset+itc)/timeUnit,errorU,errorT
    close(k)
    write(*,*) 'itc_total:',reloadItcOffset+itc,' errorU,errorT:',errorU,errorT
#endif
  end subroutine check
!===================================================================================================

!===================================================================================================
! 子程序: update_host_all
! 作用: 依次同步五套具名数组，并等待同步完成。
!===================================================================================================
  subroutine update_host_all(full)
    use commondata
    implicit none
    logical, intent(in) :: full
    ! 同一同步步只传一次宏观场；完整检查点需要时再补传f/g、力和接口时间历史。
    integer(kind=4), save :: macroItc=-1, stateItc=-1 !save 表示：退出子程序后仍然保留变量的值，下次调用继续使用。

    if(full) then
        if(stateItc==itc) return !判断是否跳过已经完成的同步。
    else
        if(macroItc==itc) return
    endif
    call update_host_grid(size(f_coarse,1),size(f_coarse,2),full,macroItc/=itc,f_coarse,g_coarse,rho_coarse,&
        u_coarse,v_coarse,T_coarse,Fx_coarse,Fy_coarse) !处理中心区域,传输哪些数据
    if(refineRatio > 1) then
        call update_host_grid(size(f_left,1),size(f_left,2),full,macroItc/=itc,f_left,g_left,rho_left,&
            u_left,v_left,T_left,Fx_left,Fy_left) !处理左区域,传输哪些数据
        call update_host_grid(size(f_right,1),size(f_right,2),full,macroItc/=itc,f_right,g_right,rho_right,&
            u_right,v_right,T_right,Fx_right,Fy_right) !处理右区域,传输哪些数据
        call update_host_grid(size(f_bottom,1),size(f_bottom,2),full,macroItc/=itc,f_bottom,g_bottom,rho_bottom,&
            u_bottom,v_bottom,T_bottom,Fx_bottom,Fy_bottom) !处理下区域,传输哪些数据
        call update_host_grid(size(f_top,1),size(f_top,2),full,macroItc/=itc,f_top,g_top,rho_top,&
            u_top,v_top,T_top,Fx_top,Fy_top) !处理上区域,传输哪些数据
    endif
    if(full) then !重启之后，细网格还需要利用这些历史继续做时间插值
        !$acc update self(coarseHistoryPrevious,coarseHistoryCurrent,coarseHistoryNext) async(1)
    endif
    !$acc wait(1)
    macroItc=itc
    if(full) stateItc=itc
  end subroutine update_host_all
!===================================================================================================

!===================================================================================================
! 子程序: update_host_grid
! 作用: 将指定区域的宏观量或完整状态从设备同步到主机。
!===================================================================================================
  subroutine update_host_grid(ni,nj,full,macros,f,g,rho,u,v,T,Fx,Fy)

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

    logical, intent(in) :: full, macros

    if(macros) then !需要传宏观场
        !$acc update self(u, v, T, rho) async(1)
    endif
    if(full) then !需要完整状态
        !$acc update self(f, g, Fx, Fy) async(1)
    endif
  end subroutine update_host_grid
!===================================================================================================

!===================================================================================================
! 子程序: check_grid
! 作用: 累计一个区域的收敛误差，并保存本次检查时的宏观场。
!===================================================================================================
! 仅稳态模式编译本子程序；非稳态不使用上一次检查场 up/vp/Tp。
#ifdef steadyFlow
  subroutine check_grid(ni,nj,dx,xOffset,yOffset,xStartIntegration,xEndIntegration,yStartIntegration,&
      yEndIntegration,integrationBox,du,uu,dt,&
      tt,u,v,T,up,vp,Tp)

    use commondata,only: centerIntegrationBox,refineRatio
    implicit none

    ! up/vp/Tp 是上一次收敛检查的场。
    ! 先按节点实际面积累计误差，再把当前 u/v/T 保存为下一次检查的参考。

    integer(kind=4) :: ni, nj, xStartIntegration, xEndIntegration, yStartIntegration, yEndIntegration
    real(kind=8) :: dx, xOffset, yOffset, integrationBox(4)

    real(kind=8) :: x, y, xlo, xhi, ylo, yhi, middleWidth
    integer(kind=4) :: i, j
    real(kind=8) :: u(ni,nj)
    real(kind=8) :: v(ni,nj)
    real(kind=8) :: T(ni,nj)
    real(kind=8) :: up(ni,nj)
    real(kind=8) :: vp(ni,nj)
    real(kind=8) :: Tp(ni,nj)
    real(kind=8) :: du, uu, dt, tt, cellArea

    do j=yStartIntegration,yEndIntegration !先遍历传入的节点范围
        do i=xStartIntegration,xEndIntegration
            ! 节点控制面积先裁剪到本区范围；x/y 是以细格距为单位的全局物理坐标。
            x=xOffset+(i-0.5d0)*dx
            y=yOffset+(j-0.5d0)*dx
            !一个节点对应的完整控制区域为：
            ![x-\Delta x/2,\ x+\Delta x/2][y-\Delta x/2,\ y+\Delta x/2].
            !完整粗节点通常按面积 4 加权；
            !完整细节点按面积 1 加权
            xlo=max(x-dx/2,integrationBox(1))  !integrationBox取交集
            xhi=min(x+dx/2,integrationBox(2))
            ylo=max(y-dx/2,integrationBox(3))
            yhi=min(y+dx/2,integrationBox(4))
            cellArea=0.0d0
            if(xhi > xlo .AND. yhi > ylo) then
                if(refineRatio==1 .OR. dx > 1.0d0) then
                    cellArea=(xhi-xlo)*(yhi-ylo)    ! 单网格全域或中心粗物理矩形。
                else
                    ! 分别累加四个细物理矩形的面积，排除中心粗物理区域。
                    cellArea=max(0.0d0,min(xhi,centerIntegrationBox(1))-xlo)*(yhi-ylo)    ! 左
                    cellArea=cellArea+max(0.0d0,xhi-max(xlo,centerIntegrationBox(2)))*(yhi-ylo)    ! 右
                    middleWidth=max(0.0d0,min(xhi,centerIntegrationBox(2))-max(xlo,centerIntegrationBox(1)))
                    cellArea=cellArea+middleWidth*max(0.0d0,min(yhi,centerIntegrationBox(3))-ylo)    ! 下
                    cellArea=cellArea+middleWidth*max(0.0d0,yhi-max(ylo,centerIntegrationBox(4)))    ! 上
                endif
            endif
            du=du+cellArea*((u(i,j)-up(i,j))**2+(v(i,j)-vp(i,j))**2) !面积乘节点误差，再累计
            uu=uu+cellArea*(u(i,j)**2+v(i,j)**2)
            ! 温度误差按实际面积累加绝对值。
            dt=dt+cellArea*dABS( T(i,j)-Tp(i,j) )
            tt=tt+cellArea*dABS( T(i,j) )
        enddo
    enddo

    up=u   !最后保存当前场
    vp=v
    Tp=T
  end subroutine check_grid
#endif
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
    !这里计算的是当前时刻的空间统计结果，不进行时间平均。
    !mass	密度的全域面积积分,meanT	全域面积平均温度
    call calculate_NuRe(nu,re,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax) !计算当前全域统计量

    !新算例初始化时已经建立文件并写好标题；续算时前面的逻辑已经保留重启时刻之前的有效历史，所以这里直接追加即可。
    open(newunit=k,file=NuReHistoryFile,status='old',position='append')
    write(k,'(2(I0,1X),12(ES24.16E3,1X))') NuReSampleCount+1,reloadItcOffset+itc, &
        dble(reloadItcOffset+itc)/timeUnit,nu,re,hot,cold,middle,&
        mass,meanT,tmin,tmax,rmin,rmax
    close(k)
    NuReSampleCount=NuReSampleCount+1
    write(*,'(a,f12.5,2(a,es13.5))') 't_ff=',dble(reloadItcOffset+itc)/timeUnit,' NuVolAvg=',nu,' ReVolRMS=',re
  end subroutine calNuRe
!===================================================================================================

!===================================================================================================
! 子程序: calculate_NuRe
! 作用: 只计算当前 Nu/Re 及场量统计，不写采样历史，不改变采样次数。
!===================================================================================================
  subroutine calculate_NuRe(nu,re,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax)

    use commondata
    implicit none

    ! 依次累加五个物理分区，节点实际贡献的面积在各区统计循环内直接计算。
    ! RMS Re 使用速度平方的面积平均再开方。
    ! 先累加全域积分，不能直接平均各网格各自计算的 Re。

    real(kind=8), intent(out) :: nu, re, hot, cold, middle, mass, meanT, tmin, tmax, rmin, rmax
    real(kind=8) :: conv, vel2

    call update_host_all(.false.) !确保 CPU 上的 rho/u/v/T 是当前场

    conv=0.0d0 !沿传热方向的速度与温度乘积的面积积分
    vel2=0.0d0 !速度大小平方的面积积分。
    mass=0.0d0
    meanT=0.0d0
    hot=0.0d0
    cold=0.0d0
    middle=0.0d0
    tmin=huge(1.0d0) !huge(1.0d0) 返回该双精度实数类型能够表示的最大有限值
    tmax=-tmin
    rmin=tmin
    rmax=-tmin

    call calNuRe_grid(size(f_coarse,1),size(f_coarse,2),dxCoarse,xOffsetCoarse,yOffsetCoarse,(1+merge(interfaceSkin,0,refineRatio>1)),&
        (size(f_coarse,1)-merge(interfaceSkin,0,refineRatio>1)),(1+merge(interfaceSkin,0,refineRatio>1)),(size(f_coarse,2)-merge(interfaceSkin,0,refineRatio>1)),centerIntegrationBox,&
        conv,vel2,mass,meanT,hot,cold,&
        middle,tmin,tmax,rmin,rmax,rho_coarse,u_coarse,v_coarse,T_coarse) !累计中心粗网格的贡献
    
    if(refineRatio > 1) then   !多块情况下，再累计四个细区

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
    ! 原 RMS Re 使用速度平方的面积平均再开方，非稳态时间统计保持同一 RMS 定义。
    nu=1.0d0+conv/(dble(nx)*dble(ny))*(lengthUnit/(Thot-Tcold))/diffusivity
    re=sqrt(vel2/(dble(nx)*dble(ny)))*lengthUnit/viscosity
    hot=(lengthUnit/(Thot-Tcold))*hot
    cold=(lengthUnit/(Thot-Tcold))*cold
    middle=(lengthUnit/(Thot-Tcold))*middle
    meanT=meanT/(dble(nx)*dble(ny))  !平均温度
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

    use commondata,only: nx,ny,Thot,Tcold,diffusivity,itc,reloadItcOffset,refineRatio,ieee_is_finite,centerIntegrationBox
    implicit none

    ! integrationBox 为细数组外框或粗区物理矩形；零面积交换节点不参与统计。
    ! fine_section_x/y 将目标全局物理坐标换成全局细网格连续下标，再跨细数组读取模板，不能截成单侧模板。

    integer(kind=4) :: ni, nj, xStartIntegration, xEndIntegration, yStartIntegration, yEndIntegration
    real(kind=8) :: dx, xOffset, yOffset, integrationBox(4)

    real(kind=8) :: x, y, xlo, xhi, ylo, yhi, middleWidth
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
            ! 节点控制面积先裁剪到本区范围；x/y 是以细格距为单位的全局物理坐标。
            x=xOffset+(i-0.5d0)*dx
            y=yOffset+(j-0.5d0)*dx
            xlo=max(x-dx/2,integrationBox(1))
            xhi=min(x+dx/2,integrationBox(2))
            ylo=max(y-dx/2,integrationBox(3))
            yhi=min(y+dx/2,integrationBox(4))
            cellArea=0.0d0
            if(xhi > xlo .AND. yhi > ylo) then
                if(refineRatio==1 .OR. dx > 1.0d0) then
                    cellArea=(xhi-xlo)*(yhi-ylo)    ! 单网格全域或中心粗物理矩形。
                else
                    ! 分别累加四个细物理矩形的面积，排除中心粗物理区域。
                    cellArea=max(0.0d0,min(xhi,centerIntegrationBox(1))-xlo)*(yhi-ylo)    ! 左
                    cellArea=cellArea+max(0.0d0,xhi-max(xlo,centerIntegrationBox(2)))*(yhi-ylo)    ! 右
                    middleWidth=max(0.0d0,min(xhi,centerIntegrationBox(2))-max(xlo,centerIntegrationBox(1)))
                    cellArea=cellArea+middleWidth*max(0.0d0,min(yhi,centerIntegrationBox(3))-ylo)    ! 下
                    cellArea=cellArea+middleWidth*max(0.0d0,yhi-max(ylo,centerIntegrationBox(4)))    ! 上
                endif
            endif
            if(cellArea <=0d0) cycle !跳过没有积分贡献的节点,跳到当前 i 循环的下一次迭代，不是退出整个子程序。
            if(.NOT. ieee_is_finite(T(i,j)) .OR. .NOT. ieee_is_finite(rho(i,j)) .OR. &
                .NOT. ieee_is_finite(u(i,j)) .OR. .NOT. ieee_is_finite(v(i,j)) .OR. rho(i,j) <=0.0d0) then
                write(*,*) 'Invalid state: time, xOffset,yOffset, i,j:',reloadItcOffset+itc,xOffset,yOffset,i,j
                error stop 'Nonfinite or nonpositive density in owned cells' !检查参与统计的节点是否出现异常。
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
    if(integrationBox(1)==0.0d0) then  !判断当前区域是否接触左壁面
        do j=yStartIntegration,yEndIntegration
            hot=hot+(8.0d0*Thot-9.0d0*T(1,j)+T(2,j))/(3.0d0*dx)*dx/dble(ny)
        enddo
    endif

    if(integrationBox(2)==dble(nx)) then  !判断当前区域是否接触右壁面
        do j=yStartIntegration,yEndIntegration
            cold=cold+(-8.0d0*Tcold+9.0d0*T(ni,j)-T(ni-1,&
                j))/(3.0d0*dx)*dx/dble(ny)
        enddo
    endif

    if((0.5d0*nx) >=integrationBox(1) .AND. (0.5d0*nx) < integrationBox(2)) then !判断中线是否经过当前区域
    !细网格传入的是整个细数组对应的矩形范围，包含伸入中心粗区的接收层。粗网格传入的是中心粗区的实际积分范围。
        call section_weights(((0.5d0*nx)-xOffset)/dx+0.5d0,ni,im,w,dw) !中线未必恰好落在某一行节点上,找到节点，计算插值系数和和拉格朗日基函数的系数的导数
        do j=yStartIntegration,yEndIntegration
            call section_owned_weight(xOffset,yOffset,dx,integrationBox,j,1,(0.5d0*nx),sectionWeight) !当前节点在这条中线上实际负责多长的一段
            if(sectionWeight <=0d0) cycle

            if(dx==1.0d0 .AND. refineRatio > 1) then
                call fine_section_x((0.5d0*nx),nint(yOffset)+j,tm,um,dTdx) !从细数组读取模板数据，算出目标位置的温度、速度和温度导数
            else
                tm=sum(w*T(im:im+3,j))
                um=sum(w*u(im:im+3,j))
                dTdx=sum(dw*T(im:im+3,j))/dx
            endif

            middle=middle+(um*tm/diffusivity-dTdx)*sectionWeight/dble(ny)
        enddo
    endif
#else

    if(integrationBox(3)==0.0d0) then  !判断当前区域是否接触下壁面
        do i=xStartIntegration,xEndIntegration
            hot=hot+(8.0d0*Thot-9.0d0*T(i,1)+T(i,2))/(3.0d0*dx)*dx/dble(nx)
        enddo
    endif

    if(integrationBox(4)==dble(ny)) then  !判断当前区域是否接触上壁面
        do i=xStartIntegration,xEndIntegration
            cold=cold+(-8.0d0*Tcold+9.0d0*T(i,nj)-T(i,&
                nj-1))/(3.0d0*dx)*dx/dble(nx)
        enddo
    endif

    if((0.5d0*ny) >=integrationBox(3) .AND. (0.5d0*ny) < integrationBox(4)) then !判断中线是否经过当前区域

        call section_weights(((0.5d0*ny)-yOffset)/dx+0.5d0,nj,jm,w,dw) !中线未必恰好落在某一行节点上,找到节点，计算插值系数和和拉格朗日基函数的系数的导数
        do i=xStartIntegration,xEndIntegration
            call section_owned_weight(xOffset,yOffset,dx,integrationBox,i,2,(0.5d0*ny),sectionWeight) !当前节点在这条中线上实际负责多长的一段
            if(sectionWeight <=0d0) cycle !贡献为零就跳过

            if(dx==1.0d0 .AND. refineRatio > 1) then !多块细网格
                call fine_section_y(nint(xOffset)+i,(0.5d0*ny),tm,vm,dTdy) !从细数组读取模板数据，算出目标位置的温度、速度和温度导数
            else
                tm=sum(w*T(i,jm:jm+3)) !中线温度，粗网格
                vm=sum(w*v(i,jm:jm+3)) !中线竖直速度，粗网格
                dTdy=sum(dw*T(i,jm:jm+3))/dx !中线温度的 y 方向导数，粗网格
            endif

            middle=middle+(vm*tm/diffusivity-dTdy)*sectionWeight/dble(nx) !计算局部热通量，乘以积分长度并累加
        enddo
    endif
#endif
  end subroutine calNuRe_grid
!===================================================================================================

!===================================================================================================
! 子程序: section_weights
! 作用: 定位中线诊断的四点模板，计算函数值及导数的插值权重。
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
    ! k 是当前区域沿积分方向的局部下标；q 和 position 都是全局物理坐标，不是连续下标。
    ! 先定位本节点控制区间；细区处于中心投影内时只累加两侧物理线段。

    real(kind=8), intent(in) :: xOffset, yOffset, dx, integrationBox(4), position
    integer(kind=4), intent(in) :: k, axis
    real(kind=8) :: q, lo, hi, clo, chi

    if(axis==1) then  !固定 x 的竖直中线
        q=yOffset+(k-.5d0)*dx
        lo=integrationBox(3)
        hi=integrationBox(4)
        clo=centerIntegrationBox(3)
        chi=centerIntegrationBox(4)
    else            !固定 y 的水平中线
        q=xOffset+(k-.5d0)*dx  !当前节点的全局物理坐标 x
        lo=integrationBox(1)   !当前网格块传入的 x 范围
        hi=integrationBox(2)   !当前网格块传入的 x 范围
        clo=centerIntegrationBox(1)  !中心粗区的 x 范围，即 127.5、895.5
        chi=centerIntegrationBox(2)  !中心粗区的 x 范围，即 127.5、895.5
    endif

    w=max(0d0,min(q+dx/2,hi)-max(q-dx/2,lo)) !计算节点控制区间[q-\Delta x/2,\ q+\Delta x/2]与当前积分范围的交集长度
    !对于粗网格，到这里就算完了。细网格还没处理完，因为它传入的范围包含接收层
    if(refineRatio > 1 .AND. dx==1d0) then !多块下的细网格
        if(axis==1) then
            if(position < centerIntegrationBox(1) .OR. position >=centerIntegrationBox(2)) return
        else
            if(position < centerIntegrationBox(3) .OR. position >=centerIntegrationBox(4)) return !判断这条截面是否穿过中心粗区
        endif  !如果截面在中心粗区的下方或上方，就不用扣除中心部分，保留前面计算的长度直接返回

        ! 如果经过，就从当前节点的控制线段中排除中心粗区负责的部分，只保留该节点落在两侧细物理区的长度
        w=max(0d0,min(q+dx/2,hi,clo)-max(q-dx/2,lo)) &    !排除重叠层的影响，只计算自己真正区域的长度
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

    ! x 为目标全局物理坐标，j 为全局细网格下标；目标 y=j-0.5。
    ! x+0.5 为全局细网格连续下标；first 为模板首个全局细网格下标，再由 fine_value 换成各区局部下标。

    real(kind=8), intent(in) :: x
    integer(kind=4), intent(in) :: j
    real(kind=8), intent(out) :: tm, um, dTdx
    real(kind=8) :: w(4), dw(4), temp(4), velocity(4)
    integer(kind=4) :: first, a

    call section_weights(x+0.5d0,nx,first,w,dw) !找到节点，计算插值系数和和拉格朗日基函数的系数的导数
    do a=1,4
        call fine_value(first+a-1,j,T_left,T_right,T_bottom,T_top,temp(a)) !按全局细网格下标选取细数组，再换成该区局部下标读取。
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

    ! i 为全局细网格下标，y 为目标全局物理坐标；目标 x=i-0.5。
    ! y+0.5 为全局细网格连续下标；first 为模板首个全局细网格下标，再由 fine_value 换成各区局部下标。

    integer(kind=4), intent(in) :: i
    real(kind=8), intent(in) :: y
    real(kind=8), intent(out) :: tm, vm, dTdy
    real(kind=8) :: w(4), dw(4), temp(4), velocity(4)
    integer(kind=4) :: first, a

    call section_weights(y+0.5d0,ny,first,w,dw) !找到节点，计算插值系数和和拉格朗日基函数的系数的导数
    do a=1,4
        call fine_value(i,first+a-1,T_left,T_right,T_bottom,T_top,temp(a)) !把四个节点的温度和速度读出来
        call fine_value(i,first+a-1,v_left,v_right,v_bottom,v_top,velocity(a)) !把四个节点的温度和速度读出来
    enddo

    tm=sum(w*temp)   !用四个样本计算目标位置的数值
    vm=sum(w*velocity)
    dTdy=sum(dw*temp)
  end subroutine fine_section_y
!===================================================================================================

!===================================================================================================
! 子程序: fine_value
! 作用: 将全局细网格下标转换成对应细区的局部下标，读取细数组中的标量值。
!===================================================================================================
  subroutine fine_value(i,j,left,right,bottom,top,value)
    use commondata,only: f_left,f_right,f_bottom,f_top,&
        xOffsetRight,xOffsetBottom,xOffsetTop,yOffsetTop
    implicit none

    real(kind=8), intent(out) :: value    ! 输出：指定全局细网格下标对应节点的温度或速度等标量值。

    ! 这里的 i/j 是全局细网格下标，对应全局物理坐标 x=i-0.5、y=j-0.5。
    ! 减去所选细区的整数偏移后才是局部下标；如全局 i=950 对应右细区局部下标 950-893=57。
    ! 只用于主机端中线后处理；若请求落入未分配的中心空区，立即报错。

    integer(kind=4), intent(in) :: i, j
    real(kind=8), intent(in) :: left(size(f_left,1),size(f_left,2)), right(size(f_right,1),size(f_right,2)), &
        bottom(size(f_bottom,1),size(f_bottom,2)),top(size(f_top,1),size(f_top,2))

    if(i <=size(f_left,1)) then
        value=left(i,j) 
    else if(i > nint(xOffsetRight)) then  !nint 取最近整数
        value=right(i-nint(xOffsetRight),j)
    else if(j <=size(f_bottom,2)) then
        value=bottom(i-nint(xOffsetBottom),j)
    else if(j > nint(yOffsetTop)) then
        value=top(i-nint(xOffsetTop),j-nint(yOffsetTop))
    else   !处于四个细数组中间的空缺区域，没有分配细网格数组来保存它的数据
        error stop 'Fine section requested an inactive central node'
    endif
  end subroutine fine_value
!===================================================================================================

!===================================================================================================
! 子程序: output_final_NuRe
! 作用: 独立输出稳态计算最后一步的诊断；达到最大步数不代表已收敛。
!===================================================================================================
  subroutine output_final_NuRe()
#ifdef steadyFlow
    use commondata,only: itc,reloadItcOffset,timeUnit,settingsFile,errorU,errorT,epsU,epsT,steadyCheckIntervalItc
    implicit none

    integer(kind=4) :: k
    real(kind=8) :: nu, re, hot, cold, middle, mass, meanT, tmin, tmax, rmin, rmax

    ! 输出最后一步的诊断，不追加采样历史，也不增加 NuReSampleCount。
    call calculate_NuRe(nu,re,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax)

    open(newunit=k,file='FinalNuRe_2DOpenaccMultiblock.dat',status='replace')
    write(k,'(a,I0)') '# itc_total = ',reloadItcOffset+itc
    write(k,'(a,L1)') '# converged = ',errorU <=epsU .AND. errorT <=epsT
    write(k,'(a,I0)') '# error_check_itc = ', &
        ((reloadItcOffset+itc)/steadyCheckIntervalItc)*steadyCheckIntervalItc
    write(k,'(a,2(ES24.16E3,1X))') '# errorU errorT = ',errorU,errorT
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
    call record_output_file('FinalNuRe_2DOpenaccMultiblock.dat')
    write(*,*) 'Final Nu_global, Re_rms:',nu,re
#endif
  end subroutine output_final_NuRe
!===================================================================================================

!===================================================================================================
! 子程序: output_Tecplot
! 作用: 输出各具名区域的坐标和宏观场，用于 Tecplot 可视化。
!===================================================================================================
  subroutine output_Tecplot()

    use commondata
    implicit none

    integer(kind=4) :: k
    character(len=40) :: timeString
    character(len=256) :: name

    call output_time_string(timeString) !!返回累计 tff 的有效字符串，100，100.5这种，timeString是字符串
    write(name,'(a,a,a,i0,a)') pltFolderPrefix//'-',trim(timeString),'tff-',reloadItcOffset+itc,'itc.dat'

    call update_host_all(.false.)  !传输数据到 CPU
    open(newunit=k,file=trim(name),status='replace')
    write(k,'(a)') 'VARIABLES="x/L","y/L","u","v","T","rho"'

    call output_Tecplot_grid(size(f_coarse,1),size(f_coarse,2),dxCoarse,xOffsetCoarse,yOffsetCoarse,(1+merge(interfaceSkin,0,refineRatio>1)),&
        (size(f_coarse,1)-merge(interfaceSkin,0,refineRatio>1)),(1+merge(interfaceSkin,0,refineRatio>1)),(size(f_coarse,2)-merge(interfaceSkin,0,refineRatio>1)),&
        'coarse',k,rho_coarse,u_coarse,&
        v_coarse,T_coarse)
    if(refineRatio > 1) then

        call output_Tecplot_grid(size(f_left,1),size(f_left,2),dxFine,xOffsetLeft,yOffsetLeft,1,size(f_left,1),&
            1,size(f_left,2),'left',k,rho_left,u_left,v_left,T_left)

        call output_Tecplot_grid(size(f_right,1),size(f_right,2),dxFine,xOffsetRight,yOffsetRight,1,&
            size(f_right,1),1,size(f_right,2),'right',k,&
            rho_right,u_right,&
            v_right,T_right)

        call output_Tecplot_grid(size(f_bottom,1),size(f_bottom,2),dxFine,xOffsetBottom,yOffsetBottom,1,&
            size(f_bottom,1),1,size(f_bottom,2),'bottom',&
            k,rho_bottom,u_bottom,&
            v_bottom,T_bottom)

        call output_Tecplot_grid(size(f_top,1),size(f_top,2),dxFine,xOffsetTop,yOffsetTop,1,size(f_top,1),&
            1,size(f_top,2),'top',k,rho_top,u_top,v_top,T_top)
    endif

    close(k)
    call record_output_file(name)    !OutputFiles输出
  end subroutine output_Tecplot
!===================================================================================================

!===================================================================================================
! 子程序: output_Tecplot_grid
! 作用: 向已打开的 Tecplot 文件写入一个区域的坐标及宏观场。
!===================================================================================================
  subroutine output_Tecplot_grid(ni,nj,dx,xOffset,yOffset,xStartIntegration,xEndIntegration,&
      yStartIntegration,yEndIntegration,&
      gridName,k,rho,u,v,T)

    use commondata,only: lengthUnit,timeUnit,itc,reloadItcOffset
    implicit none
    integer(kind=4) :: ni, nj, xStartIntegration, xEndIntegration, yStartIntegration, yEndIntegration
    real(kind=8) :: dx, xOffset, yOffset
    character(*), intent(in) :: gridName

    integer(kind=4) :: k, i, j
    real(kind=8) :: rho(ni,nj)
    real(kind=8) :: u(ni,nj)
    real(kind=8) :: v(ni,nj)
    real(kind=8) :: T(ni,nj)

    write(k,'(a,a,a,I0,a,I0,a,ES24.16E3)') 'ZONE T="',trim(gridName),'", I=',xEndIntegration-xStartIntegration+1,&
        ', J=',yEndIntegration-yStartIntegration+1,', F=POINT, SOLUTIONTIME=',dble(reloadItcOffset+itc)/timeUnit
    do j=yStartIntegration,yEndIntegration
        do i=xStartIntegration,xEndIntegration
            ! 局部下标 i/j 经偏移和 dx 换成全局物理坐标，再除以 lengthUnit；另写 u/v/T/rho，共六列。
            write(k,'(6(ES24.16E3,1X))')((xOffset+(dble(i)-0.5d0)*dx))/lengthUnit,&
                ((yOffset+(dble(j)-0.5d0)*dx))/lengthUnit,u(i,j),v(i,j),T(i,j),rho(i,j)
        enddo
    enddo
  end subroutine output_Tecplot_grid
!===================================================================================================

!===================================================================================================
! 子程序: output_SnapshotFile
! 作用: 将五区完整宏观场拼成四个一维数组，保留原始局部下标和坐标偏移，不存积分权重。
!===================================================================================================
  subroutine output_SnapshotFile()    !四个细数组全部保存，包括部分积分权重为零的接收节点
    use commondata
    implicit none
    integer(kind=4) :: k, first, last, count
    real(kind=8), allocatable :: uall(:), vall(:), Tall(:), rhoall(:)
    character(len=40) :: timeString
    character(len=256) :: name

    call update_host_all(.false.)     !传输数据到 CPU
    call output_time_string(timeString) !返回累计 tff 的有效字符串，100，100.5这种，timeString是字符串
    write(name,'(a,a,a,i0,a)') snapshotFilePrefix//'-',trim(timeString),'tff-',reloadItcOffset+itc,'itc.bin'
    ! 五区均保存完整数组，包含接收层；后处理局部下标与主程序一致，积分时再按实际物理范围计权。
    ! 默认粗区完整尺寸389*389，xOffsetCoarse/yOffsetCoarse仍为122.5，局部下标1对应全局物理坐标123.5。
    count=size(T_coarse)   !中心粗网格的完整节点数
    if(refineRatio>1) count=count+size(T_left)+size(T_right)+size(T_bottom)+size(T_top) ! 四周细网格
    allocate(uall(count),vall(count),Tall(count),rhoall(count))
    !按照“粗、左、右、下、上”的顺序拼接
    first=1
    last=size(T_coarse) !Fortran 默认第一个下标变化最快,先遍历 i，再增加 j
    uall(first:last)=reshape(u_coarse,[last-first+1]) !reshape：把二维数组按规定顺序展开成一维
    vall(first:last)=reshape(v_coarse,[last-first+1]) ![last-first+1]：指定展开后的一维数组长度
    Tall(first:last)=reshape(T_coarse,[last-first+1])
    rhoall(first:last)=reshape(rho_coarse,[last-first+1])

    if(refineRatio>1) then
        first=last+1              !左细区
        last=last+size(T_left)
        uall(first:last)=reshape(u_left,[last-first+1])
        vall(first:last)=reshape(v_left,[last-first+1])
        Tall(first:last)=reshape(T_left,[last-first+1])
        rhoall(first:last)=reshape(rho_left,[last-first+1])

        first=last+1              !右细区
        last=last+size(T_right)
        uall(first:last)=reshape(u_right,[last-first+1])
        vall(first:last)=reshape(v_right,[last-first+1])
        Tall(first:last)=reshape(T_right,[last-first+1])
        rhoall(first:last)=reshape(rho_right,[last-first+1])

        first=last+1              !下细区
        last=last+size(T_bottom)
        uall(first:last)=reshape(u_bottom,[last-first+1])
        vall(first:last)=reshape(v_bottom,[last-first+1])
        Tall(first:last)=reshape(T_bottom,[last-first+1])
        rhoall(first:last)=reshape(rho_bottom,[last-first+1])

        first=last+1              !上细区
        last=last+size(T_top)
        uall(first:last)=reshape(u_top,[last-first+1])
        vall(first:last)=reshape(v_top,[last-first+1])
        Tall(first:last)=reshape(T_top,[last-first+1])
        rhoall(first:last)=reshape(rho_top,[last-first+1])
    endif

    ! 每段内i最快、j其次，段顺序coarse/left/right/bottom/top；不插值成全域均匀网格。
    ! 只存四个原始场与少量几何/物理信息；后处理由这些信息重建面积。
    !'unformatted'：按二进制数据写入，不转换成文本数字。
    !'stream'：连续写入数据，不添加顺序非格式化文件的记录标记
    open(newunit=k,file=trim(name),form='unformatted',access='stream',status='replace')
    ! 0002格式：完整粗数组、原始偏移和实际interfaceSkin；只记录累计步数和累计tff，不分存续算计数。
    write(k) 'MB2DSNAPSHOT0002',merge(5,1,refineRatio>1),nx,ny,reloadItcOffset+itc
#ifdef SideHeatedCell
    write(k) 1
#else
    write(k) 0
#endif
    write(k) refineRatio,interfaceSkin,fineLayerCellsLeft,fineLayerCellsRight,fineLayerCellsBottom,fineLayerCellsTop
    write(k) dble(reloadItcOffset+itc)/timeUnit,lengthUnit,viscosity,diffusivity,Thot,Tcold,Tref,timeUnit
    write(k) centerIntegrationBox
    ! 每区直接保存完整尺寸和原始坐标偏移，不再用首节点坐标重新编号。
    write(k) size(T_coarse,1),size(T_coarse,2),xOffsetCoarse,yOffsetCoarse,dxCoarse
    if(refineRatio>1) then
        write(k) size(T_left,1),size(T_left,2),xOffsetLeft,yOffsetLeft,dxFine
        write(k) size(T_right,1),size(T_right,2),xOffsetRight,yOffsetRight,dxFine
        write(k) size(T_bottom,1),size(T_bottom,2),xOffsetBottom,yOffsetBottom,dxFine
        write(k) size(T_top,1),size(T_top,2),xOffsetTop,yOffsetTop,dxFine
    endif
    write(k) uall,vall,Tall,rhoall
    close(k)
    deallocate(uall,vall,Tall,rhoall)
    call record_output_file(name)  !OutputFiles输出
  end subroutine output_SnapshotFile
!===================================================================================================


!===================================================================================================
! 子程序: output_ReloadFile
! 作用: 保存完整检查点，并在写入完成后更新 latest.meta。
!===================================================================================================
  subroutine output_ReloadFile()
    use commondata
    implicit none
    integer(kind=4) :: k, m, currentModel(12)
    real(kind=8) :: currentPhysics(15)
    character(len=40) :: timeString
    character(len=256) :: name, metaName

    call update_host_all(.true.)  !传输数据到 CPU
    call output_time_string(timeString)   !返回累计 tff 的有效字符串，100，100.5这种，timeString是字符串
    ! 文件名同时记录实际t_ff和累计细格子步数，手工续算直接填写这个整数步数。
    write(name,'(a,a,a,i0,a)') reloadFilePrefix//'-',trim(timeString),'tff-',reloadItcOffset+itc,'itc.bin'
    metaName=reloadFilePrefix//'-latest.meta'
    open(newunit=k,file=trim(name),access='stream',form='unformatted',status='replace')
    ! 五区顺序固定：coarse、left、right、bottom、top；不在bin中夹入标量或几何信息。
    write(k) f_coarse,g_coarse,u_coarse,v_coarse,T_coarse,rho_coarse,Fx_coarse,Fy_coarse
    write(k) coarseHistoryPrevious,coarseHistoryCurrent,coarseHistoryNext
#ifdef steadyFlow
    write(k) up_coarse,vp_coarse,Tp_coarse
#endif
    if(refineRatio>1) then
    write(k) f_left,g_left,u_left,v_left,T_left,rho_left,Fx_left,Fy_left
#ifdef steadyFlow
    write(k) up_left,vp_left,Tp_left
#endif
    write(k) f_right,g_right,u_right,v_right,T_right,rho_right,Fx_right,Fy_right
#ifdef steadyFlow
    write(k) up_right,vp_right,Tp_right
#endif
    write(k) f_bottom,g_bottom,u_bottom,v_bottom,T_bottom,rho_bottom,Fx_bottom,Fy_bottom
#ifdef steadyFlow
    write(k) up_bottom,vp_bottom,Tp_bottom
#endif
    write(k) f_top,g_top,u_top,v_top,T_top,rho_top,Fx_top,Fy_top
#ifdef steadyFlow
    write(k) up_top,vp_top,Tp_top
#endif
    endif
    close(k)
    ! 完整数组落盘后，覆盖唯一的latest.meta，直接记录本次检查点信息。
    open(newunit=m,file=trim(metaName),status='replace')
    write(m,'(a)') 'MB2D_RESTART 1'
    write(m,'(a)') trim(name)
    write(m,'(a)') '# nx ny refineRatio Left Right Bottom Top interfaceSkin blockCount'
    write(m,*) nx,ny,refineRatio,fineLayerCellsLeft,fineLayerCellsRight, &
        fineLayerCellsBottom,fineLayerCellsTop,interfaceSkin,merge(5,1,refineRatio>1)
    currentModel=0     !将数组的12个元素全部置零
#ifdef steadyFlow
    currentModel(1)=1     !根据编译开关，把对应元素设成1
#endif
#ifdef SideHeatedCell
    currentModel(3)=1
#endif
#ifdef HorizontalWallsNoslip
    currentModel(4)=1
#endif
#ifdef VerticalWallsNoslip
    currentModel(5)=1
#endif
#ifdef VerticalWallsPeriodicalU
    currentModel(6)=1
#endif
#ifdef HorizontalWallsConstT
    currentModel(7)=1
#endif
#ifdef HorizontalWallsAdiabatic
    currentModel(8)=1
#endif
#ifdef VerticalWallsConstT
    currentModel(9)=1
#endif
#ifdef VerticalWallsAdiabatic
    currentModel(10)=1
#endif
#ifdef VerticalWallsPeriodicalT
    currentModel(11)=1
#endif
#ifdef SideHeatedHa
    currentModel(12)=1
#endif
    ! 前13项为公共物理参数和输出间隔，第14、15项为Ha、phi；非磁场模式保持0。
    currentPhysics=[Rayleigh,Prandtl,Mach,Thot,Tcold,Snu,Sq,Qk,Qnu,paraA,outputSnapshotInterval,&
        reloadFileInterval,outputPltFileInterval,0.0d0,0.0d0]  !用数组构造物理信息
#ifdef SideHeatedHa
    currentPhysics(14)=Ha
    currentPhysics(15)=phi
#endif
    write(m,'(a)') '# steady unused SideHeated UhorizontalWall UverticalWall Uperiodic '// &
        'ThorizontalConst ThorizontalAdiabatic TverticalConst TverticalAdiabatic Tperiodic Ha'
    write(m,*) currentModel
    write(m,'(a)') '# Rayleigh Prandtl Mach Thot Tcold Snu Sq Qk Qnu paraA '// &
        'snapshotInterval reloadInterval pltInterval Ha phi'
    write(m,'(15(ES24.16E3,1X))') currentPhysics
    write(m,'(a)') '# itc_total NuReSampleCount t_ff errorU errorT'
    write(m,'(2(I0,1X),3(ES24.16E3,1X))') reloadItcOffset+itc,NuReSampleCount, &
        dble(reloadItcOffset+itc)/timeUnit,errorU,errorT
#ifdef steadyFlow
    write(m,'(a)') '# steadyCheckIntervalItc'
    write(m,*) steadyCheckIntervalItc
#endif
    write(m,'(a)') '# coarse: ni nj iStart iEnd jStart jEnd / xOffset yOffset dx xmin xmax ymin ymax'
    write(m,*) size(f_coarse,1),size(f_coarse,2), &
        1+merge(interfaceSkin,0,refineRatio>1),size(f_coarse,1)-merge(interfaceSkin,0,refineRatio>1), &
        1+merge(interfaceSkin,0,refineRatio>1),size(f_coarse,2)-merge(interfaceSkin,0,refineRatio>1)
    write(m,'(7(ES24.16E3,1X))') xOffsetCoarse,yOffsetCoarse,dxCoarse, &
        centerIntegrationBox
    if(refineRatio>1) then
    write(m,'(a)') '# left: ni nj iStart iEnd jStart jEnd / xOffset yOffset dx xmin xmax ymin ymax'
    write(m,*) size(f_left,1),size(f_left,2), &
        1,size(f_left,1),1,size(f_left,2)
    write(m,'(7(ES24.16E3,1X))') xOffsetLeft,yOffsetLeft,dxFine, &
        xOffsetLeft,xOffsetLeft+size(f_left,1),yOffsetLeft,yOffsetLeft+size(f_left,2)
    write(m,'(a)') '# right: ni nj iStart iEnd jStart jEnd / xOffset yOffset dx xmin xmax ymin ymax'
    write(m,*) size(f_right,1),size(f_right,2), &
        1,size(f_right,1),1,size(f_right,2)
    write(m,'(7(ES24.16E3,1X))') xOffsetRight,yOffsetRight,dxFine, &
        xOffsetRight,xOffsetRight+size(f_right,1),yOffsetRight,yOffsetRight+size(f_right,2)
    write(m,'(a)') '# bottom: ni nj iStart iEnd jStart jEnd / xOffset yOffset dx xmin xmax ymin ymax'
    write(m,*) size(f_bottom,1),size(f_bottom,2), &
        1,size(f_bottom,1),1,size(f_bottom,2)
    write(m,'(7(ES24.16E3,1X))') xOffsetBottom,yOffsetBottom,dxFine, &
        xOffsetBottom,xOffsetBottom+size(f_bottom,1),yOffsetBottom,yOffsetBottom+size(f_bottom,2)
    write(m,'(a)') '# top: ni nj iStart iEnd jStart jEnd / xOffset yOffset dx xmin xmax ymin ymax'
    write(m,*) size(f_top,1),size(f_top,2), &
        1,size(f_top,1),1,size(f_top,2)
    write(m,'(7(ES24.16E3,1X))') xOffsetTop,yOffsetTop,dxFine, &
        xOffsetTop,xOffsetTop+size(f_top,1),yOffsetTop,yOffsetTop+size(f_top,2)
    endif
    ! 仅记录目标，不把目标作为续算兼容性约束。
    write(m,'(a)') '# requested_target_itc requested_target_tff'
    write(m,'(I0,1X,ES24.16E3)') itc_max,dble(itc_max)/timeUnit
    close(m)
  end subroutine output_ReloadFile
!===================================================================================================

!===================================================================================================

!===================================================================================================
! 输出索引只记录本程序生成的文件和准确整数步数，不依赖目录扫描或文件名时间的舍入。
! 同时供续算删除未来输出、独立后处理选择有效快照使用。
!===================================================================================================
  subroutine record_output_file(name)
    use commondata,only: reloadItcOffset,itc
    implicit none
    character(*), intent(in) :: name
    integer :: k

    open(newunit=k,file='OutputFiles_2DOpenaccMultiblock.dat',status='unknown',position='append')
    write(k,'(I0,1X,A)') reloadItcOffset+itc,trim(name)
    close(k)
  end subroutine record_output_file
!===================================================================================================

!===================================================================================================
! 子程序: exit_data_2d_openacc
! 作用: 等待设备任务结束，释放场量、粗时间层及格子常数的设备存储。
!===================================================================================================
  subroutine exit_data_2d_openacc()
    use commondata
    implicit none

    !$acc wait(1)

    !$acc exit data delete(f_coarse, g_coarse, f_post_coarse, g_post_coarse, rho_coarse, u_coarse, &
    !$acc& v_coarse, T_coarse, Fx_coarse, Fy_coarse, coarseHistoryPrevious,coarseHistoryCurrent,coarseHistoryNext)
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
    !$acc exit data delete(ex, ey, omega, omegaT)
  end subroutine exit_data_2d_openacc
!===================================================================================================
