!=============================================================
!!!    注释区，代码描述
!!!    三维浮力驱动自然对流
!!!    D3Q19 流场 + D3Q7 温度场
! Coordinate convention: z/k/w is buoyancy/vertical; y/j/v is the lateral direction.
!=============================================================

!=============================================================
!自定义宏，一些选项的开关
#define steadyFlow
!#define unsteadyFlow

! 流动模式宏的选择：稳态用残差提前停止，非稳态按指定自由落体时间运行。
! 两个都开、两个都关都会报错；只有二选一才通过。
#if defined(steadyFlow) && defined(unsteadyFlow)
#error "Choose only one flow mode: steadyFlow or unsteadyFlow"
#endif
#if !defined(steadyFlow) && !defined(unsteadyFlow)
#error "Define one flow mode: steadyFlow or unsteadyFlow"
#endif

! 速度边界：上下壁为 z-normal，左右壁为 y-normal，前后壁为 x-normal
! HorizontalWalls*: z/k 上下壁；VerticalWalls*: y/j 左右壁；SpanwiseWalls*: x/i 前后壁
#define HorizontalWallsNoslip
#define VerticalWallsNoslip
#define SpanwiseWallsNoslip
!#define HorizontalWallsPeriodicalU
!#define VerticalWallsPeriodicalU
!#define SpanwiseWallsPeriodicalU

! 速度边界宏的选择：每个方向在无滑移和周期之间二选一。
#if defined(HorizontalWallsNoslip) && defined(HorizontalWallsPeriodicalU)
#error "Choose only one horizontal velocity BC: HorizontalWallsNoslip or HorizontalWallsPeriodicalU"
#endif
#if !defined(HorizontalWallsNoslip) && !defined(HorizontalWallsPeriodicalU)
#error "Define one horizontal velocity BC: HorizontalWallsNoslip or HorizontalWallsPeriodicalU"
#endif
#if defined(VerticalWallsNoslip) && defined(VerticalWallsPeriodicalU)
#error "Choose only one y-normal/left-right velocity BC: VerticalWallsNoslip or VerticalWallsPeriodicalU"
#endif
#if !defined(VerticalWallsNoslip) && !defined(VerticalWallsPeriodicalU)
#error "Define one y-normal/left-right velocity BC: VerticalWallsNoslip or VerticalWallsPeriodicalU"
#endif
#if defined(SpanwiseWallsNoslip) && defined(SpanwiseWallsPeriodicalU)
#error "Choose only one x-normal/front-back velocity BC: SpanwiseWallsNoslip or SpanwiseWallsPeriodicalU"
#endif
#if !defined(SpanwiseWallsNoslip) && !defined(SpanwiseWallsPeriodicalU)
#error "Define one x-normal/front-back velocity BC: SpanwiseWallsNoslip or SpanwiseWallsPeriodicalU"
#endif

!温度边界(for Rayleigh Benard Cell)：上下 z/k 壁恒温，左右 y/j 和前后 x/i 壁绝热或周期
!#define RayleighBenardCell
!#define HorizontalWallsConstT
!#define VerticalWallsAdiabatic
!#define VerticalWallsPeriodicalT
!#define SpanwiseWallsAdiabatic
!#define SpanwiseWallsPeriodicalT



!温度边界(for Side Heated Cell)：左右 y/j 壁恒温，上下 z/k 和前后 x/i 壁绝热
#define SideHeatedCell
#define VerticalWallsConstT
#define HorizontalWallsAdiabatic
#define SpanwiseWallsAdiabatic
!~~temperature B.C.~~



! 对流算例宏的选择：SideHeatedCell 和 RayleighBenardCell 的热流方向、边界和后处理不同，不能混开。
#if defined(RayleighBenardCell) && defined(SideHeatedCell)
#error "Choose only one case macro: RayleighBenardCell or SideHeatedCell"
#endif
#if !defined(RayleighBenardCell) && !defined(SideHeatedCell)
#error "Define one case macro: RayleighBenardCell or SideHeatedCell"
#endif
#if defined(RayleighBenardCell) && !defined(HorizontalWallsConstT)
#error "RayleighBenardCell uses z-direction hot/cold walls: define HorizontalWallsConstT"
#endif
#if defined(HorizontalWallsConstT) && (defined(HorizontalWallsAdiabatic) || defined(HorizontalWallsPeriodicalT))
#error "Choose only one horizontal temperature BC"
#endif
#if defined(VerticalWallsConstT) && (defined(VerticalWallsAdiabatic) || defined(VerticalWallsPeriodicalT))
#error "Choose only one y-normal/left-right temperature BC"
#endif
#if defined(SpanwiseWallsAdiabatic) && defined(SpanwiseWallsPeriodicalT)
#error "Choose only one x-normal/front-back temperature BC: SpanwiseWallsAdiabatic or SpanwiseWallsPeriodicalT"
#endif



!算法切换
!启用 M1G 修正；注释掉则不使用 useG 相关修正
!#define EnableUseG
!启用旧算法
#define EnableLegacyThermalScheme

! 温度算法宏的选择：legacy D3Q7 和 EnableUseG 历史通量修正属于两套温度模型，不能混用。
#if defined(EnableUseG) && defined(EnableLegacyThermalScheme)
#error "Choose only one thermal branch: EnableUseG or EnableLegacyThermalScheme"
#endif
#if !defined(EnableUseG) && !defined(EnableLegacyThermalScheme)
#error "Define one thermal branch: EnableUseG or EnableLegacyThermalScheme"
#endif


!   自定义宏结束
!=============================================================


!=============================================================
!   全局模块
module commondata3d
  implicit none
  !===============================================================================================
  ! 格子离散速度数：流场 D3Q19，温度场 D3Q7。
  integer(kind=4), parameter :: qf=19, qt=7
  !===============================================================================================

  !===============================================================================================
  ! 是否在计算前从旧算例重启
  integer(kind=4), parameter :: loadInitField=0   ! 0: 不重启；1: 从 reloadFilePrefix-*.bin 读取初值

  ! 正常断电续算只需要设置 loadInitField=1；
  ! 代码会读取 <reloadFilePrefix>-latest.meta，并从里面找到最新的 .bin。
  ! 正常续算不用改 reloadFileNum；只有 latest .meta 缺失时，
  ! 才手动设置 reloadFileNum 作为保守推断编号。
  integer(kind=4) :: reloadFileNum=0              ! latest .meta 存在时会被覆盖；meta 缺失时作为手工兜底编号
  !===============================================================================================
  real(kind=8) :: reloadDimensionlessTime=0.0d0   ! 续算前已累计的 t_ff；优先从 latest .meta 读取，meta 缺失时由代码推断
  integer(kind=4) :: restartItcOffset=0           ! 续算前已累计的格子步数；优先从 latest .meta 读取，meta 缺失时由代码推断
  logical :: reloadMetadataLoaded=.false.         ! 标记是否成功读取 reload 元数据文件
  !===============================================================================================

  !===============================================================================================
  ! 无量纲参数
  integer(kind=4), parameter :: nx=80, ny=80, nz=80     ! 三个方向的流体节点数；物理壁面在 0 和 N 处
#ifdef SideHeatedCell
  integer(kind=4), parameter :: meshModeUniform=0, meshModeErf=1
  integer(kind=4), parameter :: meshMode=meshModeErf
  real(kind=8), parameter :: islbmStretchA=1.5d0
  real(kind=8), parameter :: islbmDxMinRaw=0.5d0*(1.0d0+erf(islbmStretchA*(1.0d0/dble(nx+1)-0.5d0))/erf(0.5d0*islbmStretchA))
  real(kind=8), parameter :: islbmDyMinRaw=0.5d0*(1.0d0+erf(islbmStretchA*(1.0d0/dble(ny+1)-0.5d0))/erf(0.5d0*islbmStretchA))
  real(kind=8), parameter :: islbmDzMinRaw=0.5d0*(1.0d0+erf(islbmStretchA*(1.0d0/dble(nz+1)-0.5d0))/erf(0.5d0*islbmStretchA))
  real(kind=8), parameter :: lengthUnit=(1.0d0-islbmDyMinRaw)/islbmDyMinRaw     ! 侧壁差温：特征长度取 y 方向 half-way 有效距离
#else
  integer(kind=4), parameter :: meshModeUniform=0, meshModeErf=1
  integer(kind=4), parameter :: meshMode=meshModeErf
  real(kind=8), parameter :: islbmStretchA=1.5d0
  real(kind=8), parameter :: islbmDxMinRaw=0.5d0*(1.0d0+erf(islbmStretchA*(1.0d0/dble(nx+1)-0.5d0))/erf(0.5d0*islbmStretchA))
  real(kind=8), parameter :: islbmDyMinRaw=0.5d0*(1.0d0+erf(islbmStretchA*(1.0d0/dble(ny+1)-0.5d0))/erf(0.5d0*islbmStretchA))
  real(kind=8), parameter :: islbmDzMinRaw=0.5d0*(1.0d0+erf(islbmStretchA*(1.0d0/dble(nz+1)-0.5d0))/erf(0.5d0*islbmStretchA))
  real(kind=8), parameter :: lengthUnit=(1.0d0-islbmDzMinRaw)/islbmDzMinRaw     ! RB 上下差温：特征长度取 z 方向 half-way 有效距离
#endif
  real(kind=8), parameter :: pi=acos(-1.0d0)

  real(kind=8), parameter :: Rayleigh=1.0d5
  real(kind=8), parameter :: Prandtl=0.71d0
  real(kind=8), parameter :: Mach=0.1d0
  real(kind=8), parameter :: Thot=0.5d0, Tcold=-0.5d0
  real(kind=8), parameter :: Tref=0.5d0*(Thot+Tcold)       ! 参考温度，计算热膨胀系数和无量纲温度用
  real(kind=8), parameter :: tauf=0.5d0+Mach*lengthUnit*dsqrt(3.0d0*Prandtl/Rayleigh)
  real(kind=8), parameter :: viscosity=(tauf-0.5d0)/3.0d0  ! 动量扩散率 nu
  real(kind=8), parameter :: diffusivity=viscosity/Prandtl ! 热扩散率 kappa

  real(kind=8), parameter :: cs2T=0.25d0                   ! UseG D3Q7 温度格子的 cs_T^2

  ! 高阶矩参数修正 aT，legacy D3Q7 温度算法使用。
  real(kind=8), parameter :: paraA=42.0d0*dsqrt(3.0d0)*diffusivity-6.0d0



  ! 速度后处理的比较标度，只影响输出诊断，不参与主时间推进。
  real(kind=8), parameter :: velocityScaleCompare=lengthUnit/diffusivity

  ! 浮力项参数；浮力沿 z 方向施加，输出变量 W_nd 是竖直速度
  real(kind=8), parameter :: gBeta1=Rayleigh*viscosity*diffusivity/lengthUnit
  real(kind=8), parameter :: gBeta=gBeta1/lengthUnit/lengthUnit
  real(kind=8), parameter :: timeUnit=dsqrt(lengthUnit/gBeta)      ! 1 个自由落体时间对应的格子步数
  real(kind=8), parameter :: velocityUnit=dsqrt(gBeta*lengthUnit)  ! 自由落体速度标度

  ! 动量方程的多松弛系数
  real(kind=8), parameter :: Se=1.0d0/tauf, Seps=1.0d0/tauf
  real(kind=8), parameter :: Snu=1.0d0/tauf, Spi=1.0d0/tauf
  real(kind=8), parameter :: Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)
  real(kind=8), parameter :: Sm=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)

  ! 温度方程的多松弛系数：legacy 使用 Xu D3Q7 参数，UseG 使用 Chai 型通量修正。
#ifdef EnableLegacyThermalScheme
  real(kind=8), parameter :: Qk=3.0d0-dsqrt(3.0d0), Qnu=4.0d0*dsqrt(3.0d0)-6.0d0
  real(kind=8), parameter :: thermalGeqCoeff=7.0d0/((6.0d0+paraA)*cs2T)
#else
  real(kind=8), parameter :: taug=0.5d0+diffusivity/cs2T
  real(kind=8), parameter :: Qnu=1.0d0, Qk=1.0d0/taug
  real(kind=8), parameter :: thermalGeqCoeff=1.0d0/cs2T
#endif
  !===============================================================================================

  !===============================================================================================
  ! 输出/备份相关设置（以自由落体时间 t_ff 为单位）
  real(kind=8), parameter :: epsU=1.0d-7, epsT=1.0d-7    ! 稳态收敛阈值

#ifdef steadyFlow
  real(kind=8), parameter :: outputSnapshotInterval=10.0d0   ! 快照和 Nu/Re 时间序列采样间隔（单位：t_ff）
  real(kind=8), parameter :: reloadFileInterval=100.0d0  ! f/g 重启文件输出间隔（单位：t_ff）
  real(kind=8), parameter :: outputPltFileInterval=100.0d0  ! Tecplot 文件输出间隔（单位：t_ff）
  integer(kind=4), parameter :: dimensionlessTimeMax=int(12000.0d0/outputSnapshotInterval)
  integer(kind=4), parameter :: outputSnapshotFile=1  ! 是否输出后处理快照文件：0=不输出，1=输出
  integer(kind=4), parameter :: outputPltFile=1       ! 是否输出 Tecplot 文件：0=不输出，1=输出
  integer(kind=4), parameter :: outputReloadFile=1    ! 是否周期输出 f/g 重启文件：0=不输出，1=输出
  integer(kind=4), parameter :: itc_max=20000000      ! 稳态最大格子步，实际可由 errorU/errorT 提前停止
#endif

#ifdef unsteadyFlow
  real(kind=8), parameter :: outputSnapshotInterval=0.5d0   ! uvwTrho 快照和 Nu/Re 时间序列采样间隔（单位：t_ff）
  real(kind=8), parameter :: reloadFileInterval=100.0d0  ! f/g 重启文件输出间隔（单位：t_ff）
  real(kind=8), parameter :: outputPltFileInterval=100.0d0  ! Tecplot 文件周期输出间隔（单位：t_ff）
  real(kind=8), parameter :: unsteadyRunDuration=1000.0d0  ! 非稳态阶段固定运行的自由落体时间长度
  ! 以下三个参数只控制非稳态结束后的 Nu/Re 统计平均窗口，不改变推进时长或采样频率。
  ! 时间窗口按完整算例的绝对 t_ff 计；续算后处理会从完整 .dat 历史重建。
  real(kind=8), parameter :: unsteadyAverageStartTf=0.5d0*unsteadyRunDuration  ! 平均窗口起点
  real(kind=8), parameter :: unsteadyAverageEndTf=unsteadyRunDuration          ! 平均窗口终点
  real(kind=8), parameter :: unsteadyAverageMidTf=0.5d0*(unsteadyAverageStartTf+unsteadyAverageEndTf) ! 前/后半分界
  integer(kind=4), parameter :: unsteadySampleCount=max(1, int(unsteadyRunDuration/outputSnapshotInterval+0.5d0)) ! 输出次数计数器
  integer(kind=4), parameter :: dimensionlessTimeMax=unsteadySampleCount
  integer(kind=4), parameter :: outputSnapshotFile=1  ! 是否输出后处理快照文件：0=不输出，1=输出
  integer(kind=4), parameter :: outputPltFile=1       ! 是否输出 Tecplot 文件：0=不输出，1=输出
  integer(kind=4), parameter :: outputReloadFile=1    ! 是否周期输出 f/g 重启文件：0=不输出，1=输出
  integer(kind=4), parameter :: itc_max=max(1, int(unsteadyRunDuration*timeUnit+0.5d0)) ! 非稳态总时长换算成格子步
#endif

  integer(kind=4) :: snapshotFileNum, pltFileNum  ! 快照/plt 输出文件计数器；reload 使用 reloadFileNum 独立编号
  integer(kind=4) :: dimensionlessTime
  integer(kind=4) :: outputSnapshotIntervalItc
  integer(kind=4) :: reloadFileIntervalItc, outputPltFileIntervalItc
  ! dimensionlessTime 与 outputSnapshotInterval 对应：
  ! 每调用一次 calNuRe() 就加 1，用于索引 NuVolAvg/ReVolAvg 和写出 t_ff 时间轴。

  real(kind=8) :: NuVolAvg(0:dimensionlessTimeMax), ReVolAvg(0:dimensionlessTimeMax)
  ! 体平均 Nu 和 Re 的时间序列缓存
  ! 只有在启用并调用 calNuRe() 的情况下这些数组才会被真正填充。

  character(len=100) :: snapshotFilePrefix="buoyancyCavity3DOpenmpSnapshot"
  ! 快照输出文件前缀，实际文件名形如：<snapshotFilePrefix>-<编号>.bin。
  character(len=100) :: pltFolderPrefix="buoyancyCavity3DOpenmpTecplot"
  ! Tecplot 输出文件前缀，实际文件名形如：<pltFolderPrefix>-<类型>-<编号>.plt。
  character(len=100) :: reloadFilePrefix="reloadFile3DOpenmp"
  ! 严格重启文件前缀，实际读取/写出：<reloadFilePrefix>-<编号>.bin。
  character(len=100) :: settingsFile="SimulationSettings3DOpenmp.txt"
  !===============================================================================================

  !===============================================================================================
  ! 计算中需要的相关参数
  real(kind=8) :: errorU, errorT

  real(kind=8) :: xp(0:nx+1), yp(0:ny+1), zp(0:nz+1)   ! 无量纲坐标数组，包括边界壁面位置
  real(kind=8) :: quadWx(1:nx), quadWy(1:ny), quadWz(1:nz), quadSumX, quadSumY, quadSumZ, quadSumVolume
  real(kind=8), parameter :: islbmShift=1.0d0/lengthUnit
  integer(kind=4) :: stream_ix(0:qf-1,1:nx,3), stream_iy(0:qf-1,1:ny,3), stream_iz(0:qf-1,1:nz,3)
  real(kind=8) :: stream_wx(0:qf-1,1:nx,3), stream_wy(0:qf-1,1:ny,3), stream_wz(0:qf-1,1:nz,3)
  logical :: stream_x_valid(0:qf-1,1:nx), stream_y_valid(0:qf-1,1:ny), stream_z_valid(0:qf-1,1:nz)
  integer(kind=4) :: stream_ixT(0:qt-1,1:nx,3), stream_iyT(0:qt-1,1:ny,3), stream_izT(0:qt-1,1:nz,3)
  real(kind=8) :: stream_wxT(0:qt-1,1:nx,3), stream_wyT(0:qt-1,1:ny,3), stream_wzT(0:qt-1,1:nz,3)
  logical :: stream_x_validT(0:qt-1,1:nx), stream_y_validT(0:qt-1,1:ny), stream_z_validT(0:qt-1,1:nz)
  real(kind=8), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:), T(:,:,:), rho(:,:,:) ! 宏观变量

#ifdef steadyFlow
  real(kind=8), allocatable :: up(:,:,:), vp(:,:,:), wp(:,:,:), Tp(:,:,:) ! 上一次 check 的场，用于稳态残差
#endif
  real(kind=8), allocatable :: f(:,:,:,:), f_post(:,:,:,:)      ! 流场分布函数及碰撞后缓冲
  real(kind=8), allocatable :: g(:,:,:,:), g_post(:,:,:,:)      ! 温度分布函数及碰撞后缓冲
  real(kind=8), allocatable :: Bx_prev(:,:,:), By_prev(:,:,:), Bz_prev(:,:,:) ! EnableUseG 的上一时间步 B=uT

  integer(kind=4) :: itc
#ifdef EnableUseG
  logical, parameter :: useG=.true.            ! 仅用于日志和收敛曲线标记
#else
  logical, parameter :: useG=.false.           ! 仅用于日志和收敛曲线标记
#endif

#ifdef EnableLegacyThermalScheme
  logical, parameter :: useLegacyThermalScheme=.true.  ! 仅用于日志和收敛曲线标记
#else
  logical, parameter :: useLegacyThermalScheme=.false. ! 仅用于日志和收敛曲线标记
#endif

#ifdef steadyFlow
  real(kind=8) :: Nu_global, Nu_hot, Nu_cold, Nu_middle ! 稳态最终诊断：全场、热壁、冷壁和中面 Nu
  real(kind=8) :: Nu_hot_max, Nu_hot_min, Nu_hot_max_position, Nu_hot_min_position ! 稳态最终诊断：热壁局部 Nu 极值及位置
#endif

  ! 格子离散速度、反向索引和权重
  integer(kind=4) :: ex(0:qf-1), ey(0:qf-1), ez(0:qf-1), opp(0:qf-1)
  integer(kind=4) :: exT(0:qt-1), eyT(0:qt-1), ezT(0:qt-1), oppT(0:qt-1)
  real(kind=8) :: omega(0:qf-1), omegaT(0:qt-1)

  !===============================================================================================
end module commondata3d


!   全局模块结束
!=============================================================


!=============================================================
!   主程序


program main3d
  use omp_lib
  use commondata3d
  implicit none

  real(kind=8) :: timeStart, timeEnd
  real(kind=8) :: timeStart2, timeEnd2
  character(len=24) :: ctime
  character(len=24) :: string
  integer(kind=4) :: time
  integer(kind=4) :: myMaxThreads
#ifdef unsteadyFlow
  integer(kind=4) :: nextSampleItc, nextSampleAbsItc, unsteadyItcRemaining
#endif


  if (loadInitField .EQ. 1) then
    open(unit=00, file=trim(settingsFile), status='unknown', position='append')
    write(00,*) '================ Restart continuation begins ================'
  else
    open(unit=00, file=trim(settingsFile), status='replace')
  endif
  string = ctime(time())
  write(00,*) 'Start: ', string
  write(00,*) 'Starting OpenMP >>>>>>'
  call OMP_set_num_threads(24)
  myMaxThreads = OMP_get_max_threads()
  write(00,*) 'Max Running threads:', myMaxThreads
  close(00)

  call initial()
#ifdef unsteadyFlow
  ! 非稳态的 itc_max 是整个算例的总目标步数；
  ! 续算时 restartItcOffset 是旧算例已经完成的步数，本次只推进剩余步数。
  unsteadyItcRemaining = max(0, itc_max - restartItcOffset)
#endif

  call CPU_TIME(timeStart)
  timeStart2 = OMP_get_wtime()

#ifdef steadyFlow
  do while( ((errorU.GT.epsU).OR.(errorT.GT.epsT)).AND.(itc.LE.itc_max) )
#endif
#ifdef unsteadyFlow
  do while( itc.LT.unsteadyItcRemaining )
#endif
  
    itc = itc+1

    call collision()

    call streaming()

    call bounceback()

    call macro()

    call collisionT()

    call streamingT()

    call bouncebackT()

    call macroT()

#ifdef steadyFlow
    ! 周期输出按累计格子步判断，续算时才能接回不断电运行应有的输出节奏。
    if(MOD(restartItcOffset+itc,2000).EQ.0) call check()
    if ((outputPltFile .EQ. 1) .AND. (mod(restartItcOffset+itc, outputPltFileIntervalItc) .EQ. 0)) then
      call output_Tecplot()
    endif
    if ((outputSnapshotFile .EQ. 1) .AND. (mod(restartItcOffset+itc, outputSnapshotIntervalItc) .EQ. 0)) then
      call output_SnapshotFile()  !稳态模式下的可选周期 uvwTrho 快照输出
    endif
    if ((outputReloadFile .EQ. 1) .AND. (mod(restartItcOffset+itc, reloadFileIntervalItc) .EQ. 0)) then
      call output_ReloadFile()
    endif
#endif

#ifdef unsteadyFlow
    do while( (reloadDimensionlessTime + real(dimensionlessTime,kind=8)*outputSnapshotInterval) &
         .LT. unsteadyRunDuration )
      ! 每个目标采样时刻都按绝对 t_ff 换算到本次运行段的 itc，续算时不会重复旧样本。
      nextSampleAbsItc = max(1, int((reloadDimensionlessTime + &
           real(dimensionlessTime+1,kind=8)*outputSnapshotInterval)*timeUnit+0.5d0))
      nextSampleItc = max(1, nextSampleAbsItc - restartItcOffset)
      if(itc.LT.nextSampleItc) exit
      call calNuRe()
      if(outputSnapshotFile.EQ.1) then
        call output_SnapshotFile()     !每 0.5 t_ff 输出一次 u、v、w、T、rho 的二进制快照文件
      endif
    enddo
    if ((outputPltFile .EQ. 1) .AND. (mod(restartItcOffset+itc, outputPltFileIntervalItc) .EQ. 0)) then
      call output_Tecplot()
    endif
    if ((outputReloadFile .EQ. 1) .AND. (mod(restartItcOffset+itc, reloadFileIntervalItc) .EQ. 0)) then
      call output_ReloadFile()
    endif
#endif
  enddo

  call CPU_TIME(timeEnd)
  timeEnd2 = OMP_get_wtime()

#ifdef steadyFlow
  call output_Tecplot()
  call output_SnapshotFile()    !输出 u、v、w、T、rho 的二进制快照文件
#endif 

#ifdef unsteadyFlow
  call output_unsteady_NuRe_postprocess()
#endif

    !===============================================================================================



    !===============================================================================================

    
#ifdef steadyFlow
!侧壁加热和RB对流的计算Nu不一样；这些最终标量诊断只用于稳态收敛后评估
#ifdef SideHeatedCell
  call SideHeatedcalc_Nu_global()
  call SideHeatedcalc_Nu_wall_avg()
  call SideHeatedcalc_Nu_zmid_wall_mean()
  call SideHeatedcalc_umid_max()
  call SideHeatedcalc_vmid_max()
  call SideHeatedcalc_wmid_max()
  call SideHeatedcalc_centerline_uv_max()
  call SideHeatedcalc_kinetic_energy_avg()
#endif

#ifdef RayleighBenardCell
  call RBcalc_Nu_global()
  call RBcalc_Nu_wall_avg()
  call RBcalc_umid_max()
  call RBcalc_vmid_max()
  call RBcalc_wmid_max()
#endif

  call calNuRe()
#endif

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) '======================================================================'
  write(00,*) 'Time (CPU) = ', real(timeEnd - timeStart, kind=8), 's'
  write(00,*) 'MLUPS = ', &
       real(dble(nx) * dble(ny) * dble(nz) * dble(itc) / &
       & max(timeEnd - timeStart, 1.0d-12) / 1.0d6, kind=8)
  write(00,*) 'Time (OMP) = ', real(timeEnd2 - timeStart2, kind=8), 's'
  write(00,*) 'MLUPS (OMP) = ', &
       real(dble(nx) * dble(ny) * dble(nz) * dble(itc) / &
       & max(timeEnd2 - timeStart2, 1.0d-12) / 1.0d6, kind=8)
#ifdef steadyFlow
  write(00,*) 'Nu_global =', Nu_global
  write(00,*) 'Nu_hot    =', Nu_hot
  write(00,*) 'Nu_cold   =', Nu_cold
  write(00,*) 'Nu_middle =', Nu_middle
#endif
  write(00,*) 'Deallocate Array......'


  deallocate(f, f_post, g, g_post)
  deallocate(u, v, w, T, rho)
#ifdef steadyFlow
  deallocate(up, vp, wp, Tp)
#endif
  deallocate(Bx_prev, By_prev, Bz_prev)


  write(00,*) 'Successfully: DNS completed!'
  string = ctime(time())
  write(00,*) 'End:   ', string
  close(00)

end program main3d

!   主程序结束
!=============================================================



!===========================================================================================================================
! 子程序: initial
! 作用: 初始化网格坐标、场变量、分布函数、输出文件和重启信息。
! 用途: 在主程序进入时间推进前调用，完成三维算例的全部启动准备。
!===========================================================================================================================
subroutine initial()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k, alpha
  real(kind=8) :: eu, u2Loc
  real(kind=8) :: xLen, zLen, rbInitPerturbAmp
  character(len=100) :: reloadFileName

  itc = 0
  errorU = 100.0d0
  errorT = 100.0d0
  snapshotFileNum = 0
  pltFileNum = 0
  restartItcOffset = 0
  reloadMetadataLoaded = .false.

  ! 把按自由落体时间给出的快照/备份间隔换算成格子步数 itc
  outputSnapshotIntervalItc = max(1, int(outputSnapshotInterval * timeUnit + 0.5d0))
  reloadFileIntervalItc = max(1, int(reloadFileInterval * timeUnit + 0.5d0))
  outputPltFileIntervalItc = max(1, int(outputPltFileInterval * timeUnit + 0.5d0))

  !-----------------------------------------------------------------------------------------------
  ! 记录各种信息在日志文件中
  !-----------------------------------------------------------------------------------------------
  open(unit=00, file=trim(settingsFile), status='unknown', position='append')

  if (outputSnapshotFile .EQ. 1) then
    open(unit=01, file=trim(snapshotFilePrefix)//'-'//'readme', status='unknown')
    write(01,*) 'snapshotFile folder exist!'
    close(01)
    write(00,*) 'Snapshot data will be stored in ', trim(snapshotFilePrefix)
  endif

  if (outputPltFile .EQ. 1) then
    open(unit=01, file=trim(pltFolderPrefix)//'-'//'readme', status='unknown')
    write(01,*) 'pltFile folder exist!'
    close(01)
    write(00,*) 'Data will be stored in ', trim(pltFolderPrefix)
  endif

  if (outputReloadFile .EQ. 1) then
    open(unit=01, file=trim(reloadFilePrefix)//'-'//'readme', status='unknown')
    write(01,*) 'reloadFile prefix exists!'
    close(01)
    write(00,*) 'Reload data will be stored in ', trim(reloadFilePrefix)
  endif

#ifdef EnableLegacyThermalScheme
  if ((paraA .GE. 1.0d0) .OR. (paraA .LE. -6.0d0)) then
    write(00,*) '----------------------------------'
    write(00,*) 'paraA =', paraA
    write(00,*) 'Error: condition not met for the legacy 3D thermal algorithm'
    write(00,*) 'Ref: Lattice Boltzmann simulations of three-dimensional thermal convective flows at high Rayleigh number'
    write(00,*) 'Please try to reduce Mach number'
    write(00,*) '----------------------------------'
    close(00)
    stop
  endif
#endif

  write(00,*) '-------------------------------------------------------------------------------'
  write(00,*) 'Mesh:', nx, ny, nz
  write(00,*) 'Rayleigh=', real(Rayleigh,kind=8), '; Prandtl =', real(Prandtl,kind=8), '; Mach =', real(Mach,kind=8)
  write(00,*) 'Length unit: L0 =', real(lengthUnit,kind=8)
  write(00,*) 'gravityDirection = z ; vertical velocity component = W_nd'
#ifdef SideHeatedCell
  write(00,*) 'thermalWallDistanceDirection = y'
#endif
#ifdef RayleighBenardCell
  write(00,*) 'thermalWallDistanceDirection = z'
#endif
  write(00,*) 'Time unit: Sqrt(L0/(gBeta*DeltaT)) =', real(timeUnit,kind=8)
  write(00,*) 'Velocity unit: Sqrt(gBeta*L0*DeltaT) =', real(velocityUnit,kind=8)
  write(00,*) '   '
  write(00,*) 'tauf =', real(tauf,kind=8)
#ifdef EnableLegacyThermalScheme
  write(00,*) 'thermalScheme = legacy D3Q7 '
  write(00,*) 'Qk =', real(Qk,kind=8), '; Qnu =', real(Qnu,kind=8), '; paraA =', real(paraA,kind=8)
#else
  write(00,*) 'thermalScheme = current D3Q7 (EnableUseG branch)'
  write(00,*) 'taug =', real(taug,kind=8), '; cs2T =', real(cs2T,kind=8)
  write(00,*) 'Qk =', real(Qk,kind=8), '; Qnu =', real(Qnu,kind=8)
#endif
  write(00,*) 'thermalGeqCoeff =', real(thermalGeqCoeff,kind=8)
  write(00,*) 'viscosity =', real(viscosity,kind=8), '; diffusivity =', real(diffusivity,kind=8)
  write(00,*) 'outputSnapshotInterval =', real(outputSnapshotInterval,kind=8), ' free-fall time units'
  write(00,*) 'outputSnapshotIntervalItc =', outputSnapshotIntervalItc, ' in itc units'
  write(00,*) 'outputPltFile =', outputPltFile
  write(00,*) 'outputPltFileInterval =', real(outputPltFileInterval,kind=8), ' free-fall time units'
  write(00,*) 'outputPltFileIntervalItc =', outputPltFileIntervalItc, ' in itc units'
  write(00,*) 'outputReloadFile =', outputReloadFile
  write(00,*) 'reloadFileInterval =', real(reloadFileInterval,kind=8), ' free-fall time units'
  write(00,*) 'reloadFileIntervalItc =', reloadFileIntervalItc, ' in itc units'
#ifdef unsteadyFlow
  write(00,*) 'unsteadyRunDuration =', real(unsteadyRunDuration,kind=8), ' free-fall time units'
  write(00,*) 'unsteadySampleCount =', unsteadySampleCount
#endif
  if (loadInitField .EQ. 1) then
    write(00,*) 'Restart offsets will be read from reload metadata when available.'
  endif
  write(00,*) 'itc_max =', itc_max
  write(00,*) 'default epsU =', real(epsU,kind=8), '; epsT =', real(epsT,kind=8)
  write(00,*) 'useG =', useG
  write(00,*) 'useLegacyThermalScheme =', useLegacyThermalScheme
  write(00,*) '    '

#ifdef RayleighBenardCell
  write(00,*) 'I am 3D periodic Rayleigh-Benard Cell'
#endif
#ifdef SideHeatedCell
  write(00,*) 'I am 3D Side Heated Closed Cavity'
#endif

#ifdef steadyFlow
  write(00,*) 'I am steadyFlow'
#endif
#ifdef unsteadyFlow
  write(00,*) 'I am unsteadyFlow'
#endif
  write(00,*) 'OpenMP only; MPI/OpenACC are not included in this file'

  !-----------------------------------------------------------------------------------------------
  
  
  
  !-----------------------------------------------------------------------------------------------
  ! ISLBM 节点坐标、积分权重和三维插值迁移模板
  !-----------------------------------------------------------------------------------------------
  call init_lattice_constants()
  call build_islbm_mesh_3d()
  call build_islbm_quadrature_3d()
  call prepare_islbm_streaming_stencils_3d()
  write(00,*) "ISLBM meshMode =", meshMode, "; stretchA =", real(islbmStretchA,kind=8)
  write(00,*) "ISLBM effective lengthUnit L0 =", real(lengthUnit,kind=8)
  write(00,*) "ISLBM streaming shift =", real(islbmShift,kind=8)
  write(00,*) "ISLBM quadrature volume =", real(quadSumVolume,kind=8)



  allocate(u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), T(nx,ny,nz), rho(nx,ny,nz))
#ifdef steadyFlow
  allocate(up(nx,ny,nz), vp(nx,ny,nz), wp(nx,ny,nz), Tp(nx,ny,nz))
#endif
  allocate(f(0:qf-1,nx,ny,nz), f_post(0:qf-1,0:nx+1,0:ny+1,0:nz+1))
  allocate(g(0:qt-1,nx,ny,nz), g_post(0:qt-1,0:nx+1,0:ny+1,0:nz+1))
  allocate(Bx_prev(nx,ny,nz), By_prev(nx,ny,nz), Bz_prev(nx,ny,nz))

  !-----------------------------------------------------------------------------------------------
  
  
  
  !-----------------------------------------------------------------------------------------------
  ! 初始化
  !-----------------------------------------------------------------------------------------------
  call init_lattice_constants()  !速度集和权重

  rho = 1.0d0
  f = 0.0d0
  g = 0.0d0
  f_post = 0.0d0
  g_post = 0.0d0
  Bx_prev = 0.0d0
  By_prev = 0.0d0
  Bz_prev = 0.0d0

  if (loadInitField .EQ. 0) then
    write(00,*) 'Initial field is set exactly'
    if (reloadDimensionlessTime .NE. 0.0d0) then
      write(00,*) 'Error: since loadInitField .EQ. 0, reloadDimensionlessTime should also be 0'
      close(00)
      stop
    endif

#ifdef VerticalWallsNoslip
    write(00,*) 'Velocity B.C. for vertical y/j left-right walls are: ===No-slip wall==='
#endif
#ifdef VerticalWallsPeriodicalU
    write(00,*) 'Velocity B.C. for vertical y/j left-right walls are: ===Periodical==='
#endif
#ifdef HorizontalWallsNoslip
    write(00,*) 'Velocity B.C. for horizontal z/k top-bottom walls are: ===No-slip wall==='
#endif
#ifdef HorizontalWallsPeriodicalU
    write(00,*) 'Velocity B.C. for horizontal z/k top-bottom walls are: ===Periodical==='
#endif
#ifdef SpanwiseWallsNoslip
    write(00,*) 'Velocity B.C. for spanwise x/i front-back walls are: ===No-slip wall==='
#endif
#ifdef SpanwiseWallsPeriodicalU
    write(00,*) 'Velocity B.C. for spanwise x/i front-back walls are: ===Periodical==='
#endif

    u = 0.0d0
    v = 0.0d0
    w = 0.0d0
    T = 0.0d0

#ifdef VerticalWallsConstT
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          T(i,j,k) = Thot + (yp(j) - yp(0)) / (yp(ny+1) - yp(0)) * (Tcold - Thot)
        enddo
      enddo
    enddo
    write(00,*) 'Temperature B.C. for vertical y/j left-right walls are: ===Hot/cold wall==='
#endif

#ifdef VerticalWallsAdiabatic
    write(00,*) 'Temperature B.C. for vertical y/j left-right walls are: ===Adiabatic wall==='
#endif

#ifdef VerticalWallsPeriodicalT
    write(00,*) 'Temperature B.C. for vertical y/j left-right walls are: ===Periodical==='
#endif

#ifdef HorizontalWallsAdiabatic
    write(00,*) 'Temperature B.C. for horizontal z/k top-bottom walls are: ===Adiabatic wall==='
#endif

#ifdef HorizontalWallsPeriodicalT
    write(00,*) 'Temperature B.C. for horizontal z/k top-bottom walls are: ===Periodical==='
#endif

#ifdef HorizontalWallsConstT
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          T(i,j,k) = Thot + (zp(k) - zp(0)) / (zp(nz+1) - zp(0)) * (Tcold - Thot)
        enddo
      enddo
    enddo
    write(00,*) 'Temperature B.C. for horizontal z/k top-bottom walls are: ===Hot/cold wall==='
#ifdef RayleighBenardCell
    if (Rayleigh .LE. 1.0d4) then
      xLen = xp(nx+1)
      zLen = zp(nz+1)
      rbInitPerturbAmp = 1.0d-3 * (Thot - Tcold)
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            T(i,j,k) = T(i,j,k) + rbInitPerturbAmp * dcos(2.0d0 * pi * xp(i) / xLen) * dsin(pi * zp(k) / zLen)
          enddo
        enddo
      enddo
      write(00,'(a,1x,ES24.16E3)') '3D RB initial T perturbation amplitude =', rbInitPerturbAmp
    else
      write(00,*) '3D RB initial T perturbation skipped because Rayleigh > 1.0d4'
    endif
#endif
#endif

#ifdef SpanwiseWallsAdiabatic
    write(00,*) 'Temperature B.C. for spanwise x/i front-back walls are: ===Adiabatic wall==='
#endif

#ifdef SpanwiseWallsPeriodicalT
    write(00,*) 'Temperature B.C. for spanwise x/i front-back walls are: ===Periodical==='
#endif

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          u2Loc = u(i,j,k) * u(i,j,k) + v(i,j,k) * v(i,j,k) + w(i,j,k) * w(i,j,k)
          do alpha = 0, qf-1
            eu = dble(ex(alpha)) * u(i,j,k) + dble(ey(alpha)) * v(i,j,k) + dble(ez(alpha)) * w(i,j,k)
            f(alpha,i,j,k) = omega(alpha) * rho(i,j,k) * (1.0d0 + 3.0d0 * eu + 4.5d0 * eu * eu - 1.5d0 * u2Loc)
          enddo
          do alpha = 0, qt-1
            eu = dble(exT(alpha)) * u(i,j,k) + dble(eyT(alpha)) * v(i,j,k) + dble(ezT(alpha)) * w(i,j,k)
            g(alpha,i,j,k) = omegaT(alpha) * T(i,j,k) * (1.0d0 + thermalGeqCoeff * eu)
          enddo
#ifdef EnableUseG
            Bx_prev(i,j,k) = u(i,j,k) * T(i,j,k)
            By_prev(i,j,k) = v(i,j,k) * T(i,j,k)
            Bz_prev(i,j,k) = w(i,j,k) * T(i,j,k)
#endif
        enddo
      enddo
    enddo

  elseif (loadInitField .EQ. 1) then
    ! 正常断电续算时，先读取 <reloadFilePrefix>-latest.meta；
    ! meta 会告诉代码实际要读哪个 <reloadFilePrefix>-*.bin，以及旧算例已经累计到哪里。
    write(00,*) 'Read reload metadata before choosing the restart .bin file.'
    write(reloadFileName,'(i12.12)') reloadFileNum             ! latest .meta 缺失时才依赖这个手工编号
    reloadFileName = adjustl(reloadFileName)                  ! adjustl 把前导空格移到字符串末尾
    call read_reload_metadata(reloadFileName)
    write(00,*) 'Load initial field from previous simulation: ', &
         trim(reloadFilePrefix), '-', trim(reloadFileName), '.bin'
    if (.not. reloadMetadataLoaded) then
      write(00,*) 'WARNING: no reload metadata file found; restart offsets were inferred.'
      write(00,*) '         For exact continuation, use reload files written after this patch.'
    endif
    open(unit=01, file=trim(reloadFilePrefix)//'-'//trim(reloadFileName)//'.bin', &
         form='unformatted', access='sequential', status='old')
    ! Strict restart files store f and g. With EnableUseG, they also store
    ! Bx_prev/By_prev/Bz_prev because the M1G correction needs previous heat-flux history.
    write(00,*) 'Reloading f and g from file'
    read(01) ((((f(alpha,i,j,k), i=1,nx), j=1,ny), k=1,nz), alpha=0,qf-1)
    read(01) ((((g(alpha,i,j,k), i=1,nx), j=1,ny), k=1,nz), alpha=0,qt-1)
#ifdef EnableUseG
    write(00,*) 'Reloading Bx_prev, By_prev and Bz_prev from file'
    read(01) (((Bx_prev(i,j,k), i=1,nx), j=1,ny), k=1,nz)
    read(01) (((By_prev(i,j,k), i=1,nx), j=1,ny), k=1,nz)
    read(01) (((Bz_prev(i,j,k), i=1,nx), j=1,ny), k=1,nz)
#endif
    close(01)
    call reconstruct_macro_from_fg()
    write(00,*) 'Raw data is loaded from the file: ', trim(reloadFilePrefix), '-', trim(reloadFileName), '.bin'
    write(00,*) 'Restart offset itc =', restartItcOffset
    write(00,*) 'Restart offset time_tf =', real(reloadDimensionlessTime,kind=8)
    write(00,*) 'Continue output counters: snapshot/plt/reload =', snapshotFileNum, pltFileNum, reloadFileNum
  else
    write(00,*) 'Error: initial field is not properly set'
    close(00)
    stop
  endif

  write(00,*) '-------------------------------------------------------------------------------'
  close(00)

  f_post = 0.0d0
  g_post = 0.0d0
  if (loadInitField .EQ. 0) then
    snapshotFileNum = 0
    pltFileNum = 0
    reloadFileNum = 0
    restartItcOffset = 0
    reloadDimensionlessTime = 0.0d0
#ifdef steadyFlow
    ! 新算例第一段收敛误差应从初始场开始比较。
    up = u
    vp = v
    wp = w
    Tp = T
#endif
  else
#ifdef steadyFlow
    ! 重启后第一段收敛误差应从载入场继续比较。
    up = u
    vp = v
    wp = w
    Tp = T
#endif
  endif
  dimensionlessTime = 0
  ! 新算例：清零，开始记录新历史。
  ! 续算：也清零，但不是丢旧历史；旧历史在 .dat 文件里，新数组只记录本次续算段。
  ! 写出时间轴时会叠加 reloadDimensionlessTime，所以不会从 0 t_ff 重新编号。
  NuVolAvg = 0.0d0
  ReVolAvg = 0.0d0

end subroutine initial


!===========================================================================================================================
! 子程序: init_lattice_constants
! 作用: 初始化 D3Q19 / D3Q7 的离散速度、反向索引和权重。
!===========================================================================================================================
subroutine init_lattice_constants()
  use commondata3d
  implicit none

  ! D3Q19 顺序与与 D3Q7 的编号同步（前 7 个速度一样）
  ex  = (/ 0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0 /)
  ey  = (/ 0,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1 /)
  ez  = (/ 0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1 /)
  opp = (/ 0,  2,  1,  4,  3,  6,  5, 10,  9,  8,  7, 14, 13, 12, 11, 18, 17, 16, 15 /)

  omega(0) = 1.0d0 / 3.0d0
  omega(1:6) = 1.0d0 / 18.0d0
  omega(7:18) = 1.0d0 / 36.0d0

  ! D3Q7 温度场
  exT  = (/ 0,  1, -1,  0,  0,  0,  0 /)
  eyT  = (/ 0,  0,  0,  1, -1,  0,  0 /)
  ezT  = (/ 0,  0,  0,  0,  0,  1, -1 /)
  oppT = (/ 0,  2,  1,  4,  3,  6,  5 /)

  ! 温度权重按所选热算法分支设置
#ifdef EnableLegacyThermalScheme
  omegaT(0) = (1.0d0 - paraA) / 7.0d0
  omegaT(1:6) = (6.0d0 + paraA) / 42.0d0
#else
  omegaT(0) = 1.0d0 / 4.0d0
  omegaT(1:6) = 1.0d0 / 8.0d0
#endif

end subroutine init_lattice_constants




!===========================================================================================================================
subroutine build_islbm_mesh_3d()
  use commondata3d
  implicit none
  integer(kind=4) :: i, j, k
  real(kind=8) :: rawX(0:nx+1), rawY(0:ny+1), rawZ(0:nz+1)
  real(kind=8) :: erfNorm, wall0, wall1, len

  if(meshMode.EQ.meshModeErf) then
    erfNorm = erf(0.5d0*islbmStretchA)
    do i = 0, nx+1
      rawX(i) = 0.5d0*(1.0d0 + erf(islbmStretchA*(dble(i)/dble(nx+1)-0.5d0))/erfNorm)
    enddo
    do j = 0, ny+1
      rawY(j) = 0.5d0*(1.0d0 + erf(islbmStretchA*(dble(j)/dble(ny+1)-0.5d0))/erfNorm)
    enddo
    do k = 0, nz+1
      rawZ(k) = 0.5d0*(1.0d0 + erf(islbmStretchA*(dble(k)/dble(nz+1)-0.5d0))/erfNorm)
    enddo
  else
    rawX(0) = 0.0d0; rawX(nx+1) = 1.0d0
    do i = 1, nx
      rawX(i) = (dble(i)-0.5d0)/dble(nx)
    enddo
    rawY(0) = 0.0d0; rawY(ny+1) = 1.0d0
    do j = 1, ny
      rawY(j) = (dble(j)-0.5d0)/dble(ny)
    enddo
    rawZ(0) = 0.0d0; rawZ(nz+1) = 1.0d0
    do k = 1, nz
      rawZ(k) = (dble(k)-0.5d0)/dble(nz)
    enddo
  endif

  wall0 = merge(0.5d0*rawX(1), 0.0d0, meshMode.EQ.meshModeErf)
  wall1 = merge(1.0d0-0.5d0*rawX(1), 1.0d0, meshMode.EQ.meshModeErf)
  len = wall1 - wall0
  xp(0) = 0.0d0; xp(nx+1) = 1.0d0
  do i = 1, nx
    xp(i) = (rawX(i)-wall0)/len
  enddo

  wall0 = merge(0.5d0*rawY(1), 0.0d0, meshMode.EQ.meshModeErf)
  wall1 = merge(1.0d0-0.5d0*rawY(1), 1.0d0, meshMode.EQ.meshModeErf)
  len = wall1 - wall0
  yp(0) = 0.0d0; yp(ny+1) = 1.0d0
  do j = 1, ny
    yp(j) = (rawY(j)-wall0)/len
  enddo

  wall0 = merge(0.5d0*rawZ(1), 0.0d0, meshMode.EQ.meshModeErf)
  wall1 = merge(1.0d0-0.5d0*rawZ(1), 1.0d0, meshMode.EQ.meshModeErf)
  len = wall1 - wall0
  zp(0) = 0.0d0; zp(nz+1) = 1.0d0
  do k = 1, nz
    zp(k) = (rawZ(k)-wall0)/len
  enddo

  return
end subroutine build_islbm_mesh_3d

subroutine build_islbm_quadrature_3d()
  use commondata3d
  implicit none
  integer(kind=4) :: i, j, k

  do i = 1, nx
    quadWx(i) = 0.5d0*(xp(i+1)-xp(i-1))
  enddo
  do j = 1, ny
    quadWy(j) = 0.5d0*(yp(j+1)-yp(j-1))
  enddo
  do k = 1, nz
    quadWz(k) = 0.5d0*(zp(k+1)-zp(k-1))
  enddo
  quadSumX = sum(quadWx)
  quadSumY = sum(quadWy)
  quadSumZ = sum(quadWz)
  quadSumVolume = quadSumX*quadSumY*quadSumZ

  return
end subroutine build_islbm_quadrature_3d

subroutine prepare_islbm_streaming_stencils_3d()
  use commondata3d
  implicit none
  integer(kind=4) :: alpha, i, j, k, idx(3)
  real(kind=8) :: ww(3), target
  logical :: ok

  stream_ix = 1; stream_iy = 1; stream_iz = 1
  stream_wx = 0.0d0; stream_wy = 0.0d0; stream_wz = 0.0d0
  stream_x_valid = .false.; stream_y_valid = .false.; stream_z_valid = .false.
  do alpha = 0, qf-1
    do i = 1, nx
      target = xp(i) - dble(ex(alpha))*islbmShift
      call build_lagrange_stencil_1d(nx, xp(1:nx), target, idx, ww, ok)
      stream_x_valid(alpha,i) = ok; stream_ix(alpha,i,:) = idx; stream_wx(alpha,i,:) = ww
    enddo
    do j = 1, ny
      target = yp(j) - dble(ey(alpha))*islbmShift
      call build_lagrange_stencil_1d(ny, yp(1:ny), target, idx, ww, ok)
      stream_y_valid(alpha,j) = ok; stream_iy(alpha,j,:) = idx; stream_wy(alpha,j,:) = ww
    enddo
    do k = 1, nz
      target = zp(k) - dble(ez(alpha))*islbmShift
      call build_lagrange_stencil_1d(nz, zp(1:nz), target, idx, ww, ok)
      stream_z_valid(alpha,k) = ok; stream_iz(alpha,k,:) = idx; stream_wz(alpha,k,:) = ww
    enddo
  enddo

  stream_ixT = 1; stream_iyT = 1; stream_izT = 1
  stream_wxT = 0.0d0; stream_wyT = 0.0d0; stream_wzT = 0.0d0
  stream_x_validT = .false.; stream_y_validT = .false.; stream_z_validT = .false.
  do alpha = 0, qt-1
    do i = 1, nx
      target = xp(i) - dble(exT(alpha))*islbmShift
      call build_lagrange_stencil_1d(nx, xp(1:nx), target, idx, ww, ok)
      stream_x_validT(alpha,i) = ok; stream_ixT(alpha,i,:) = idx; stream_wxT(alpha,i,:) = ww
    enddo
    do j = 1, ny
      target = yp(j) - dble(eyT(alpha))*islbmShift
      call build_lagrange_stencil_1d(ny, yp(1:ny), target, idx, ww, ok)
      stream_y_validT(alpha,j) = ok; stream_iyT(alpha,j,:) = idx; stream_wyT(alpha,j,:) = ww
    enddo
    do k = 1, nz
      target = zp(k) - dble(ezT(alpha))*islbmShift
      call build_lagrange_stencil_1d(nz, zp(1:nz), target, idx, ww, ok)
      stream_z_validT(alpha,k) = ok; stream_izT(alpha,k,:) = idx; stream_wzT(alpha,k,:) = ww
    enddo
  enddo

  return
end subroutine prepare_islbm_streaming_stencils_3d

subroutine build_lagrange_stencil_1d(n, xnodes, target, idx, ww, ok)
  implicit none
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in) :: xnodes(n), target
  integer(kind=4), intent(out) :: idx(3)
  real(kind=8), intent(out) :: ww(3)
  logical, intent(out) :: ok
  integer(kind=4) :: mid
  real(kind=8) :: xloc(3)
  real(kind=8), parameter :: tol=1.0d-12

  idx = (/1, 1, 1/)
  ww = 0.0d0
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
  xloc = (/xnodes(idx(1)), xnodes(idx(2)), xnodes(idx(3))/)
  call lagrange_weights_3(xloc, target, ww)
  ok = .true.

  return
end subroutine build_lagrange_stencil_1d

subroutine lagrange_weights_3(xnode, x0, ww)
  implicit none
  real(kind=8), intent(in) :: xnode(3), x0
  real(kind=8), intent(out) :: ww(3)

  ww(1) = ((x0-xnode(2))*(x0-xnode(3)))/((xnode(1)-xnode(2))*(xnode(1)-xnode(3)))
  ww(2) = ((x0-xnode(1))*(x0-xnode(3)))/((xnode(2)-xnode(1))*(xnode(2)-xnode(3)))
  ww(3) = ((x0-xnode(1))*(x0-xnode(2)))/((xnode(3)-xnode(1))*(xnode(3)-xnode(2)))

  return
end subroutine lagrange_weights_3

! 子程序: collision
! 作用: 流场碰撞步骤，在矩空间完成松弛并加入浮力源项修正。
! 用途: 在主程序时间推进循环中调用，位于 streaming 之前。
!===========================================================================================================================
subroutine collision()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k, alpha
  real(kind=8) :: m(0:qf-1), meq(0:qf-1), m_post(0:qf-1)
  real(kind=8) :: s(0:qf-1), fSource(0:qf-1)
  real(kind=8) :: rhoLoc, uLoc, vLoc, wLoc, u2, uDotF
  real(kind=8) :: FxLoc, FyLoc, FzLoc

  ! 流场采用 D3Q19 MRT。
  !$omp parallel do collapse(3) default(none) shared(f,f_post,rho,u,v,w,T) &
  !$omp private(i,j,k,alpha,m,meq,m_post,s,fSource,rhoLoc,uLoc,vLoc,wLoc,u2,uDotF,FxLoc,FyLoc,FzLoc)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        rhoLoc = rho(i,j,k)
        uLoc = u(i,j,k)
        vLoc = v(i,j,k)
        wLoc = w(i,j,k)
        u2 = uLoc * uLoc + vLoc * vLoc + wLoc * wLoc

        !-----------------------------------------------------------------------------------------
        ! D3Q19 正变换: m = M * f
        !-----------------------------------------------------------------------------------------
        m(0) = f(0,i,j,k) + f(1,i,j,k) + f(2,i,j,k) + f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + &
             f(6,i,j,k) + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + &
             f(12,i,j,k) + f(13,i,j,k) + f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k)

        m(1) = -30.0d0 * f(0,i,j,k) - 11.0d0 * (f(1,i,j,k) + f(2,i,j,k) + f(3,i,j,k) + f(4,i,j,k) + &
             f(5,i,j,k) + f(6,i,j,k)) + 8.0d0 * (f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + &
             f(11,i,j,k) + f(12,i,j,k) + f(13,i,j,k) + f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) + &
             f(17,i,j,k) + f(18,i,j,k))

        m(2) = 12.0d0 * f(0,i,j,k) - 4.0d0 * (f(1,i,j,k) + f(2,i,j,k) + f(3,i,j,k) + f(4,i,j,k) + &
             f(5,i,j,k) + f(6,i,j,k)) + (f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + &
             f(11,i,j,k) + f(12,i,j,k) + f(13,i,j,k) + f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) + &
             f(17,i,j,k) + f(18,i,j,k))

        m(3) =  f(1,i,j,k) - f(2,i,j,k) + f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) - f(10,i,j,k) + &
             f(11,i,j,k) - f(12,i,j,k) + f(13,i,j,k) - f(14,i,j,k)

        m(4) = -4.0d0 * f(1,i,j,k) + 4.0d0 * f(2,i,j,k) + f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) - &
             f(10,i,j,k) + f(11,i,j,k) - f(12,i,j,k) + f(13,i,j,k) - f(14,i,j,k)

        m(5) =  f(3,i,j,k) - f(4,i,j,k) + f(7,i,j,k) + f(8,i,j,k) - f(9,i,j,k) - f(10,i,j,k) + &
             f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) - f(18,i,j,k)
             
        m(6) = -4.0d0 * f(3,i,j,k) + 4.0d0 * f(4,i,j,k) + f(7,i,j,k) + f(8,i,j,k) - f(9,i,j,k) - &
             f(10,i,j,k) + f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) - f(18,i,j,k)

        m(7) =  f(5,i,j,k) - f(6,i,j,k) + f(11,i,j,k) + f(12,i,j,k) - f(13,i,j,k) - f(14,i,j,k) + &
             f(15,i,j,k) + f(16,i,j,k) - f(17,i,j,k) - f(18,i,j,k)

        m(8) = -4.0d0 * f(5,i,j,k) + 4.0d0 * f(6,i,j,k) + f(11,i,j,k) + f(12,i,j,k) - f(13,i,j,k) - &
             f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) - f(17,i,j,k) - f(18,i,j,k)

        m(9) =  2.0d0 * (f(1,i,j,k) + f(2,i,j,k)) - (f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + f(6,i,j,k)) + &
             (f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + f(12,i,j,k) + &
             f(13,i,j,k) + f(14,i,j,k)) - 2.0d0 * (f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k))

        m(10) = -4.0d0 * (f(1,i,j,k) + f(2,i,j,k)) + 2.0d0 * (f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + f(6,i,j,k)) + &
             (f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + f(12,i,j,k) + &
             f(13,i,j,k) + f(14,i,j,k)) - 2.0d0 * (f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k))

        m(11) =  (f(3,i,j,k) + f(4,i,j,k)) - (f(5,i,j,k) + f(6,i,j,k)) + &
             (f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k)) - &
             (f(11,i,j,k) + f(12,i,j,k) + f(13,i,j,k) + f(14,i,j,k))

        m(12) = -2.0d0 * (f(3,i,j,k) + f(4,i,j,k)) + 2.0d0 * (f(5,i,j,k) + f(6,i,j,k)) + &
             (f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k)) - &
             (f(11,i,j,k) + f(12,i,j,k) + f(13,i,j,k) + f(14,i,j,k))

        m(13) =  f(7,i,j,k) - f(8,i,j,k) - f(9,i,j,k) + f(10,i,j,k)

        m(14) =  f(15,i,j,k) - f(16,i,j,k) - f(17,i,j,k) + f(18,i,j,k)

        m(15) =  f(11,i,j,k) - f(12,i,j,k) - f(13,i,j,k) + f(14,i,j,k)

        m(16) =  f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) - f(10,i,j,k) - f(11,i,j,k) + f(12,i,j,k) - &
             f(13,i,j,k) + f(14,i,j,k)

        m(17) = -f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(15,i,j,k) - f(16,i,j,k) + &
             f(17,i,j,k) - f(18,i,j,k)

        m(18) =  f(11,i,j,k) + f(12,i,j,k) - f(13,i,j,k) - f(14,i,j,k) - f(15,i,j,k) - f(16,i,j,k) + &
             f(17,i,j,k) + f(18,i,j,k)

        !-----------------------------------------------------------------------------------------
        ! 平衡矩
        !-----------------------------------------------------------------------------------------
        meq(0)  = rhoLoc
        meq(1)  = rhoLoc * (-11.0d0 + 19.0d0 * u2)
        meq(2)  = rhoLoc * (3.0d0 - 11.0d0 * u2 / 2.0d0)
        meq(3)  = rhoLoc * uLoc
        meq(4)  = -2.0d0 * rhoLoc * uLoc / 3.0d0
        meq(5)  = rhoLoc * vLoc
        meq(6)  = -2.0d0 * rhoLoc * vLoc / 3.0d0
        meq(7)  = rhoLoc * wLoc
        meq(8)  = -2.0d0 * rhoLoc * wLoc / 3.0d0
        meq(9)  = rhoLoc * (2.0d0 * uLoc * uLoc - vLoc * vLoc - wLoc * wLoc)
        meq(10) = -0.5d0 * meq(9)
        meq(11) = rhoLoc * (vLoc * vLoc - wLoc * wLoc)
        meq(12) = -0.5d0 * meq(11)
        meq(13) = rhoLoc * uLoc * vLoc
        meq(14) = rhoLoc * vLoc * wLoc
        meq(15) = rhoLoc * uLoc * wLoc
        meq(16) = 0.0d0
        meq(17) = 0.0d0
        meq(18) = 0.0d0

        !-----------------------------------------------------------------------------------------
        ! 松弛系数
        !-----------------------------------------------------------------------------------------
        s(0)  = 0.0d0   !! s_rho
        s(1)  = Se      !! s_e
        s(2)  = Seps    !! s_eps
        s(3)  = 0.0d0   !! s_jx
        s(4)  = Sq      !! s_qx
        s(5)  = 0.0d0   !! s_jy
        s(6)  = Sq      !! s_qy
        s(7)  = 0.0d0   !! s_jz
        s(8)  = Sq      !! s_qz
        s(9)  = Snu     !! s_pxx
        s(10) = Spi     !! s_pixx
        s(11) = Snu     !! s_pww
        s(12) = Spi     !! s_piww
        s(13) = Snu     !! s_pxy
        s(14) = Snu     !! s_pyz
        s(15) = Snu     !! s_pxz
        s(16) = Sm      !! s_mx
        s(17) = Sm      !! s_my
        s(18) = Sm      !! s_mz

        !-----------------------------------------------------------------------------------------
        ! 体力项
        !-----------------------------------------------------------------------------------------
        FxLoc = 0.0d0
        FyLoc = 0.0d0
        FzLoc = rhoLoc * gBeta * (T(i,j,k) - Tref)
        uDotF = uLoc * FxLoc + vLoc * FyLoc + wLoc * FzLoc

        fSource(0)  = 0.0d0
        fSource(1)  = (1.0d0 - 0.5d0 * s(1))  * 38.0d0 * uDotF
        fSource(2)  = (1.0d0 - 0.5d0 * s(2))  * (-11.0d0) * uDotF
        fSource(3)  = (1.0d0 - 0.5d0 * s(3))  * FxLoc
        fSource(4)  = (1.0d0 - 0.5d0 * s(4))  * (-2.0d0 / 3.0d0) * FxLoc
        fSource(5)  = (1.0d0 - 0.5d0 * s(5))  * FyLoc
        fSource(6)  = (1.0d0 - 0.5d0 * s(6))  * (-2.0d0 / 3.0d0) * FyLoc
        fSource(7)  = (1.0d0 - 0.5d0 * s(7))  * FzLoc
        fSource(8)  = (1.0d0 - 0.5d0 * s(8))  * (-2.0d0 / 3.0d0) * FzLoc
        fSource(9)  = (1.0d0 - 0.5d0 * s(9))  * (4.0d0 * uLoc * FxLoc - 2.0d0 * vLoc * FyLoc - 2.0d0 * wLoc * FzLoc)
        fSource(10) = (1.0d0 - 0.5d0 * s(10)) * (-2.0d0 * uLoc * FxLoc + vLoc * FyLoc + wLoc * FzLoc)
        fSource(11) = (1.0d0 - 0.5d0 * s(11)) * (2.0d0 * vLoc * FyLoc - 2.0d0 * wLoc * FzLoc)
        fSource(12) = (1.0d0 - 0.5d0 * s(12)) * (-vLoc * FyLoc + wLoc * FzLoc)
        fSource(13) = (1.0d0 - 0.5d0 * s(13)) * (uLoc * FyLoc + vLoc * FxLoc)
        fSource(14) = (1.0d0 - 0.5d0 * s(14)) * (vLoc * FzLoc + wLoc * FyLoc)
        fSource(15) = (1.0d0 - 0.5d0 * s(15)) * (uLoc * FzLoc + wLoc * FxLoc)
        fSource(16) = 0.0d0
        fSource(17) = 0.0d0
        fSource(18) = 0.0d0

        !-----------------------------------------------------------------------------------------
        ! 矩空间碰撞
        !-----------------------------------------------------------------------------------------
        do alpha = 0, qf-1
          m_post(alpha) = m(alpha) - s(alpha) * (m(alpha) - meq(alpha)) + fSource(alpha)
        enddo

        !-----------------------------------------------------------------------------------------
        ! D3Q19 逆变换: f_post = M^{-1} * m_post
        !-----------------------------------------------------------------------------------------
        f_post(0,i,j,k) =  m_post(0) / 19.0d0 - 5.0d0 * m_post(1) / 399.0d0 + m_post(2) / 21.0d0

        f_post(1,i,j,k) =  m_post(0) / 19.0d0 - 11.0d0 * m_post(1) / 2394.0d0 - m_post(2) / 63.0d0 + &
             m_post(3) / 10.0d0 - m_post(4) / 10.0d0 + m_post(9) / 18.0d0 - m_post(10) / 18.0d0

        f_post(2,i,j,k) =  m_post(0) / 19.0d0 - 11.0d0 * m_post(1) / 2394.0d0 - m_post(2) / 63.0d0 - &
             m_post(3) / 10.0d0 + m_post(4) / 10.0d0 + m_post(9) / 18.0d0 - m_post(10) / 18.0d0

        f_post(3,i,j,k) =  m_post(0) / 19.0d0 - 11.0d0 * m_post(1) / 2394.0d0 - m_post(2) / 63.0d0 + &
             m_post(5) / 10.0d0 - m_post(6) / 10.0d0 - m_post(9) / 36.0d0 + m_post(10) / 36.0d0 + &
             m_post(11) / 12.0d0 - m_post(12) / 12.0d0

        f_post(4,i,j,k) =  m_post(0) / 19.0d0 - 11.0d0 * m_post(1) / 2394.0d0 - m_post(2) / 63.0d0 - &
             m_post(5) / 10.0d0 + m_post(6) / 10.0d0 - m_post(9) / 36.0d0 + m_post(10) / 36.0d0 + &
             m_post(11) / 12.0d0 - m_post(12) / 12.0d0

        f_post(5,i,j,k) =  m_post(0) / 19.0d0 - 11.0d0 * m_post(1) / 2394.0d0 - m_post(2) / 63.0d0 + &
             m_post(7) / 10.0d0 - m_post(8) / 10.0d0 - m_post(9) / 36.0d0 + m_post(10) / 36.0d0 - &
             m_post(11) / 12.0d0 + m_post(12) / 12.0d0

        f_post(6,i,j,k) =  m_post(0) / 19.0d0 - 11.0d0 * m_post(1) / 2394.0d0 - m_post(2) / 63.0d0 - &
             m_post(7) / 10.0d0 + m_post(8) / 10.0d0 - m_post(9) / 36.0d0 + m_post(10) / 36.0d0 - &
             m_post(11) / 12.0d0 + m_post(12) / 12.0d0

        f_post(7,i,j,k) =  m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 + &
             m_post(3) / 10.0d0 + m_post(4) / 40.0d0 + m_post(5) / 10.0d0 + m_post(6) / 40.0d0 + &
             m_post(9) / 36.0d0 + m_post(10) / 72.0d0 + m_post(11) / 12.0d0 + m_post(12) / 24.0d0 + &
             m_post(13) / 4.0d0 + m_post(16) / 8.0d0 - m_post(17) / 8.0d0

        f_post(8,i,j,k) =  m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 - &
             m_post(3) / 10.0d0 - m_post(4) / 40.0d0 + m_post(5) / 10.0d0 + m_post(6) / 40.0d0 + &
             m_post(9) / 36.0d0 + m_post(10) / 72.0d0 + m_post(11) / 12.0d0 + m_post(12) / 24.0d0 - &
             m_post(13) / 4.0d0 - m_post(16) / 8.0d0 - m_post(17) / 8.0d0

        f_post(9,i,j,k) =  m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 + &
             m_post(3) / 10.0d0 + m_post(4) / 40.0d0 - m_post(5) / 10.0d0 - m_post(6) / 40.0d0 + &
             m_post(9) / 36.0d0 + m_post(10) / 72.0d0 + m_post(11) / 12.0d0 + m_post(12) / 24.0d0 - &
             m_post(13) / 4.0d0 + m_post(16) / 8.0d0 + m_post(17) / 8.0d0

        f_post(10,i,j,k) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 - &
             m_post(3) / 10.0d0 - m_post(4) / 40.0d0 - m_post(5) / 10.0d0 - m_post(6) / 40.0d0 + &
             m_post(9) / 36.0d0 + m_post(10) / 72.0d0 + m_post(11) / 12.0d0 + m_post(12) / 24.0d0 + &
             m_post(13) / 4.0d0 - m_post(16) / 8.0d0 + m_post(17) / 8.0d0

        f_post(11,i,j,k) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 + &
             m_post(3) / 10.0d0 + m_post(4) / 40.0d0 + m_post(7) / 10.0d0 + m_post(8) / 40.0d0 + &
             m_post(9) / 36.0d0 + m_post(10) / 72.0d0 - m_post(11) / 12.0d0 - m_post(12) / 24.0d0 + &
             m_post(15) / 4.0d0 - m_post(16) / 8.0d0 + m_post(18) / 8.0d0

        f_post(12,i,j,k) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 - &
             m_post(3) / 10.0d0 - m_post(4) / 40.0d0 + m_post(7) / 10.0d0 + m_post(8) / 40.0d0 + &
             m_post(9) / 36.0d0 + m_post(10) / 72.0d0 - m_post(11) / 12.0d0 - m_post(12) / 24.0d0 - &
             m_post(15) / 4.0d0 + m_post(16) / 8.0d0 + m_post(18) / 8.0d0

        f_post(13,i,j,k) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 + &
             m_post(3) / 10.0d0 + m_post(4) / 40.0d0 - m_post(7) / 10.0d0 - m_post(8) / 40.0d0 + &
             m_post(9) / 36.0d0 + m_post(10) / 72.0d0 - m_post(11) / 12.0d0 - m_post(12) / 24.0d0 - &
             m_post(15) / 4.0d0 - m_post(16) / 8.0d0 - m_post(18) / 8.0d0

        f_post(14,i,j,k) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 - &
             m_post(3) / 10.0d0 - m_post(4) / 40.0d0 - m_post(7) / 10.0d0 - m_post(8) / 40.0d0 + &
             m_post(9) / 36.0d0 + m_post(10) / 72.0d0 - m_post(11) / 12.0d0 - m_post(12) / 24.0d0 + &
             m_post(15) / 4.0d0 + m_post(16) / 8.0d0 - m_post(18) / 8.0d0

        f_post(15,i,j,k) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 + &
             m_post(5) / 10.0d0 + m_post(6) / 40.0d0 + m_post(7) / 10.0d0 + m_post(8) / 40.0d0 - &
             m_post(9) / 18.0d0 - m_post(10) / 36.0d0 + m_post(14) / 4.0d0 + m_post(17) / 8.0d0 - &
             m_post(18) / 8.0d0

        f_post(16,i,j,k) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 - &
             m_post(5) / 10.0d0 - m_post(6) / 40.0d0 + m_post(7) / 10.0d0 + m_post(8) / 40.0d0 - &
             m_post(9) / 18.0d0 - m_post(10) / 36.0d0 - m_post(14) / 4.0d0 - m_post(17) / 8.0d0 - &
             m_post(18) / 8.0d0

        f_post(17,i,j,k) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 + &
             m_post(5) / 10.0d0 + m_post(6) / 40.0d0 - m_post(7) / 10.0d0 - m_post(8) / 40.0d0 - &
             m_post(9) / 18.0d0 - m_post(10) / 36.0d0 - m_post(14) / 4.0d0 + m_post(17) / 8.0d0 + &
             m_post(18) / 8.0d0

        f_post(18,i,j,k) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 - &
             m_post(5) / 10.0d0 - m_post(6) / 40.0d0 - m_post(7) / 10.0d0 - m_post(8) / 40.0d0 - &
             m_post(9) / 18.0d0 - m_post(10) / 36.0d0 + m_post(14) / 4.0d0 - m_post(17) / 8.0d0 + &
             m_post(18) / 8.0d0
      enddo
    enddo
  enddo
  !$omp end parallel do

end subroutine collision


!===========================================================================================================================
! 子程序: streaming
! 作用: 对流场分布函数执行三维 pull streaming。
! 用途: 在主程序时间推进循环中调用，位于 collision 之后、bounceback 之前。
!===========================================================================================================================
subroutine streaming()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k, ii, jj, kk, alpha
  real(kind=8) :: value

  f = 0.0d0
  !$omp parallel do collapse(3) default(none) shared(f,f_post,stream_x_valid,stream_y_valid,stream_z_valid,stream_ix,stream_iy,stream_iz,stream_wx,stream_wy,stream_wz) private(i,j,k,ii,jj,kk,alpha,value)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        do alpha = 0, qf-1
          if(stream_x_valid(alpha,i).AND.stream_y_valid(alpha,j).AND.stream_z_valid(alpha,k)) then
            value = 0.0d0
            do kk = 1, 3
              do jj = 1, 3
                do ii = 1, 3
                  value = value + stream_wx(alpha,i,ii)*stream_wy(alpha,j,jj)*stream_wz(alpha,k,kk) * &
                    f_post(alpha,stream_ix(alpha,i,ii),stream_iy(alpha,j,jj),stream_iz(alpha,k,kk))
                enddo
              enddo
            enddo
            f(alpha,i,j,k) = value
          endif
        enddo
      enddo
    enddo
  enddo
  !$omp end parallel do

end subroutine streaming


!===========================================================================================================================
! 子程序: bounceback
! 作用: 施加流场边界条件，包括无滑移壁面和周期边界配套处理。
! 用途: 在主程序时间推进循环中调用，位于 streaming 之后、macro 之前。
!===========================================================================================================================
subroutine bounceback()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k

#ifdef SpanwiseWallsPeriodicalU
  !$omp parallel do collapse(2) default(none) shared(f,f_post) private(j,k)
  do k = 1, nz
    do j = 1, ny
      ! Front boundary (i=1): incoming populations with ex=+1
      f(1, 1,j,k) = f_post(1, nx,j,k)
      f(7, 1,j,k) = f_post(7, nx,j,k)
      f(9, 1,j,k) = f_post(9, nx,j,k)
      f(11,1,j,k) = f_post(11,nx,j,k)
      f(13,1,j,k) = f_post(13,nx,j,k)

      ! Back boundary (i=nx): incoming populations with ex=-1
      f(2, nx,j,k) = f_post(2, 1,j,k)
      f(8, nx,j,k) = f_post(8, 1,j,k)
      f(10,nx,j,k) = f_post(10,1,j,k)
      f(12,nx,j,k) = f_post(12,1,j,k)
      f(14,nx,j,k) = f_post(14,1,j,k)
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef SpanwiseWallsNoslip
  !$omp parallel do collapse(2) default(none) shared(f,f_post) private(j,k)
  do k = 1, nz
    do j = 1, ny
      ! Front no-slip wall (i=1): incoming populations with ex=+1
      f(1, 1,j,k) = f_post(2, 1,j,k)
      f(7, 1,j,k) = f_post(10,1,j,k)
      f(9, 1,j,k) = f_post(8, 1,j,k)
      f(11,1,j,k) = f_post(14,1,j,k)
      f(13,1,j,k) = f_post(12,1,j,k)

      ! Back no-slip wall (i=nx): incoming populations with ex=-1
      f(2, nx,j,k) = f_post(1, nx,j,k)
      f(8, nx,j,k) = f_post(9, nx,j,k)
      f(10,nx,j,k) = f_post(7, nx,j,k)
      f(12,nx,j,k) = f_post(13,nx,j,k)
      f(14,nx,j,k) = f_post(11,nx,j,k)
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef VerticalWallsPeriodicalU
  !$omp parallel do collapse(2) default(none) shared(f,f_post) private(i,k)
  do k = 1, nz
    do i = 1, nx
      ! Left boundary (j=1): incoming populations with ey=+1
      f(3, i,1,k) = f_post(3, i,ny,k)
      f(7, i,1,k) = f_post(7, i,ny,k)
      f(8, i,1,k) = f_post(8, i,ny,k)
      f(15,i,1,k) = f_post(15,i,ny,k)
      f(17,i,1,k) = f_post(17,i,ny,k)

      ! Right boundary (j=ny): incoming populations with ey=-1
      f(4, i,ny,k) = f_post(4, i,1,k)
      f(9, i,ny,k) = f_post(9, i,1,k)
      f(10,i,ny,k) = f_post(10,i,1,k)
      f(16,i,ny,k) = f_post(16,i,1,k)
      f(18,i,ny,k) = f_post(18,i,1,k)
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef VerticalWallsNoslip
  !$omp parallel do collapse(2) default(none) shared(f,f_post) private(i,k)
  do k = 1, nz
    do i = 1, nx
      ! Left no-slip wall (j=1): incoming populations with ey=+1
      f(3, i,1,k) = f_post(4, i,1,k)
      f(7, i,1,k) = f_post(10,i,1,k)
      f(8, i,1,k) = f_post(9, i,1,k)
      f(15,i,1,k) = f_post(18,i,1,k)
      f(17,i,1,k) = f_post(16,i,1,k)

      ! Right no-slip wall (j=ny): incoming populations with ey=-1
      f(4, i,ny,k) = f_post(3, i,ny,k)
      f(9, i,ny,k) = f_post(8, i,ny,k)
      f(10,i,ny,k) = f_post(7, i,ny,k)
      f(16,i,ny,k) = f_post(17,i,ny,k)
      f(18,i,ny,k) = f_post(15,i,ny,k)
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef HorizontalWallsPeriodicalU
  !$omp parallel do collapse(2) default(none) shared(f,f_post) private(i,j)
  do j = 1, ny
    do i = 1, nx
      ! Bottom boundary (k=1): incoming populations with ez=+1
      f(5, i,j,1) = f_post(5, i,j,nz)
      f(11,i,j,1) = f_post(11,i,j,nz)
      f(12,i,j,1) = f_post(12,i,j,nz)
      f(15,i,j,1) = f_post(15,i,j,nz)
      f(16,i,j,1) = f_post(16,i,j,nz)

      ! Top boundary (k=nz): incoming populations with ez=-1
      f(6, i,j,nz) = f_post(6, i,j,1)
      f(13,i,j,nz) = f_post(13,i,j,1)
      f(14,i,j,nz) = f_post(14,i,j,1)
      f(17,i,j,nz) = f_post(17,i,j,1)
      f(18,i,j,nz) = f_post(18,i,j,1)
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef HorizontalWallsNoslip
  !$omp parallel do collapse(2) default(none) shared(f,f_post) private(i,j)
  do j = 1, ny
    do i = 1, nx
      ! Bottom no-slip wall (k=1): incoming populations with ez=+1
      f(5, i,j,1) = f_post(6, i,j,1)
      f(11,i,j,1) = f_post(14,i,j,1)
      f(12,i,j,1) = f_post(13,i,j,1)
      f(15,i,j,1) = f_post(18,i,j,1)
      f(16,i,j,1) = f_post(17,i,j,1)

      ! Top no-slip wall (k=nz): incoming populations with ez=-1
      f(6, i,j,nz) = f_post(5, i,j,nz)
      f(13,i,j,nz) = f_post(12,i,j,nz)
      f(14,i,j,nz) = f_post(11,i,j,nz)
      f(17,i,j,nz) = f_post(16,i,j,nz)
      f(18,i,j,nz) = f_post(15,i,j,nz)
    enddo
  enddo
  !$omp end parallel do
#endif

end subroutine bounceback


!===========================================================================================================================
! 子程序: macro
! 作用: 由流场分布函数恢复 rho、u、v、w 以及浮力项。
! 用途: 在主程序时间推进循环中调用，作为流场更新链条的最后一步。
!===========================================================================================================================
subroutine macro()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k
  real(kind=8) :: FzLoc

  !$omp parallel do collapse(3) default(none) &
  !$omp& shared(f,rho,u,v,w,T) &
  !$omp& private(i,j,k,FzLoc)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        rho(i,j,k) = f(0,i,j,k) + f(1,i,j,k) + f(2,i,j,k) + f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + &
             f(6,i,j,k) + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + &
             f(12,i,j,k) + f(13,i,j,k) + f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k)

        FzLoc = rho(i,j,k) * gBeta * (T(i,j,k) - Tref)

        u(i,j,k) = ( f(1,i,j,k) - f(2,i,j,k) + f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) - f(10,i,j,k) + &
             f(11,i,j,k) - f(12,i,j,k) + f(13,i,j,k) - f(14,i,j,k) ) / rho(i,j,k)

        v(i,j,k) = ( f(3,i,j,k) - f(4,i,j,k) + f(7,i,j,k) + f(8,i,j,k) - f(9,i,j,k) - f(10,i,j,k) + &
             f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) - f(18,i,j,k) ) / rho(i,j,k)

        w(i,j,k) = ( f(5,i,j,k) - f(6,i,j,k) + f(11,i,j,k) + f(12,i,j,k) - f(13,i,j,k) - f(14,i,j,k) + &
             f(15,i,j,k) + f(16,i,j,k) - f(17,i,j,k) - f(18,i,j,k) + 0.5d0 * FzLoc ) / rho(i,j,k)
      enddo
    enddo
  enddo
  !$omp end parallel do

end subroutine macro


!===========================================================================================================================
! 子程序: collisionT
! 作用: 温度场碰撞步骤
! 用途: 在主程序时间推进循环中调用，位于流场 macro 之后。
!===========================================================================================================================
subroutine collisionT()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k
  real(kind=8) :: n(0:qt-1), neq(0:qt-1), q(0:qt-1), n_post(0:qt-1)
  real(kind=8) :: Bx, By, Bz, dBx, dBy, dBz
  real(kind=8), parameter :: SG = 1.0d0 - 0.5d0 * Qk

  ! 温度场采用 D3Q7 MRT
  !$omp parallel do collapse(3) default(none) shared(g,g_post,u,v,w,T,Bx_prev,By_prev,Bz_prev) &
  !$omp private(i,j,k,n,neq,q,n_post,Bx,By,Bz,dBx,dBy,dBz)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        Bx = u(i,j,k) * T(i,j,k)
        By = v(i,j,k) * T(i,j,k)
        Bz = w(i,j,k) * T(i,j,k)

#ifdef EnableUseG
          dBx = Bx - Bx_prev(i,j,k)
          dBy = By - By_prev(i,j,k)
          dBz = Bz - Bz_prev(i,j,k)
          Bx_prev(i,j,k) = Bx
          By_prev(i,j,k) = By
          Bz_prev(i,j,k) = Bz
#else
          dBx = 0.0d0
          dBy = 0.0d0
          dBz = 0.0d0
#endif

        !-----------------------------------------------------------------------------------------
        ! D3Q7 正变换: n = M * g
        !-----------------------------------------------------------------------------------------
        n(0) = g(0,i,j,k) + g(1,i,j,k) + g(2,i,j,k) + g(3,i,j,k) + g(4,i,j,k) + g(5,i,j,k) + g(6,i,j,k)
        n(1) = g(1,i,j,k) - g(2,i,j,k)
        n(2) = g(3,i,j,k) - g(4,i,j,k)
        n(3) = g(5,i,j,k) - g(6,i,j,k)
        n(4) = -6.0d0 * g(0,i,j,k) + g(1,i,j,k) + g(2,i,j,k) + g(3,i,j,k) + g(4,i,j,k) + g(5,i,j,k) + g(6,i,j,k)
        n(5) = 2.0d0 * g(1,i,j,k) + 2.0d0 * g(2,i,j,k) - g(3,i,j,k) - g(4,i,j,k) - g(5,i,j,k) - g(6,i,j,k)
        n(6) = g(3,i,j,k) + g(4,i,j,k) - g(5,i,j,k) - g(6,i,j,k)

        ! 平衡矩
        neq(0) = T(i,j,k)
        neq(1) = Bx
        neq(2) = By
        neq(3) = Bz
#ifdef EnableLegacyThermalScheme
        neq(4) = paraA * T(i,j,k)
#else
        neq(4) = -0.75d0 * T(i,j,k)
#endif
        neq(5) = 0.0d0
        neq(6) = 0.0d0

        q(0) = 0.0d0
        q(1) = Qk
        q(2) = Qk
        q(3) = Qk
        q(4) = Qnu
        q(5) = Qnu
        q(6) = Qnu

        n_post(0) = n(0) - q(0) * (n(0) - neq(0))
        n_post(1) = n(1) - q(1) * (n(1) - neq(1)) + SG * dBx
        n_post(2) = n(2) - q(2) * (n(2) - neq(2)) + SG * dBy
        n_post(3) = n(3) - q(3) * (n(3) - neq(3)) + SG * dBz
        n_post(4) = n(4) - q(4) * (n(4) - neq(4))
        n_post(5) = n(5) - q(5) * (n(5) - neq(5))
        n_post(6) = n(6) - q(6) * (n(6) - neq(6))

        !-----------------------------------------------------------------------------------------
        ! D3Q7 逆变换: g_post = M^{-1} * n_post
        !-----------------------------------------------------------------------------------------
        g_post(0,i,j,k) = n_post(0) / 7.0d0 - n_post(4) / 7.0d0
        g_post(1,i,j,k) = n_post(0) / 7.0d0 + n_post(1) / 2.0d0 + n_post(4) / 42.0d0 + n_post(5) / 6.0d0
        g_post(2,i,j,k) = n_post(0) / 7.0d0 - n_post(1) / 2.0d0 + n_post(4) / 42.0d0 + n_post(5) / 6.0d0
        g_post(3,i,j,k) = n_post(0) / 7.0d0 + n_post(2) / 2.0d0 + n_post(4) / 42.0d0 - n_post(5) / 12.0d0 + &
             n_post(6) / 4.0d0
        g_post(4,i,j,k) = n_post(0) / 7.0d0 - n_post(2) / 2.0d0 + n_post(4) / 42.0d0 - n_post(5) / 12.0d0 + &
             n_post(6) / 4.0d0
        g_post(5,i,j,k) = n_post(0) / 7.0d0 + n_post(3) / 2.0d0 + n_post(4) / 42.0d0 - n_post(5) / 12.0d0 - &
             n_post(6) / 4.0d0
        g_post(6,i,j,k) = n_post(0) / 7.0d0 - n_post(3) / 2.0d0 + n_post(4) / 42.0d0 - n_post(5) / 12.0d0 - &
             n_post(6) / 4.0d0
      enddo
    enddo
  enddo
  !$omp end parallel do

  return
end subroutine collisionT
!===========================================================================================================================
! collisionT 结束: 完成温度分布函数 g 的碰撞更新
!===========================================================================================================================


!===========================================================================================================================
! 子程序: streamingT
! 作用: 对温度分布函数执行三维 pull streaming。
! 用途: 在主程序时间推进循环中调用，位于 collisionT 之后、bouncebackT 之前。
!===========================================================================================================================
subroutine streamingT()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k, ii, jj, kk, alpha
  real(kind=8) :: value

  g = 0.0d0
  !$omp parallel do collapse(3) default(none) shared(g,g_post,stream_x_validT,stream_y_validT,stream_z_validT,stream_ixT,stream_iyT,stream_izT,stream_wxT,stream_wyT,stream_wzT) private(i,j,k,ii,jj,kk,alpha,value)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        do alpha = 0, qt-1
          if(stream_x_validT(alpha,i).AND.stream_y_validT(alpha,j).AND.stream_z_validT(alpha,k)) then
            value = 0.0d0
            do kk = 1, 3
              do jj = 1, 3
                do ii = 1, 3
                  value = value + stream_wxT(alpha,i,ii)*stream_wyT(alpha,j,jj)*stream_wzT(alpha,k,kk) * &
                    g_post(alpha,stream_ixT(alpha,i,ii),stream_iyT(alpha,j,jj),stream_izT(alpha,k,kk))
                enddo
              enddo
            enddo
            g(alpha,i,j,k) = value
          endif
        enddo
      enddo
    enddo
  enddo
  !$omp end parallel do

  return
end subroutine streamingT
!===========================================================================================================================
! streamingT 结束: 完成温度分布函数 g 的迁移，把碰撞后的温度信息传播到相邻格点。
!===========================================================================================================================


!===========================================================================================================================
! 子程序: bouncebackT
! 作用: 施加温度边界条件，包括恒温、绝热和周期边界。
! 用途: 在主程序时间推进循环中调用，位于 streamingT 之后、macroT 之前。
!===========================================================================================================================
subroutine bouncebackT()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k

#ifdef SpanwiseWallsPeriodicalT
  !$omp parallel do collapse(2) default(none) shared(g,g_post) private(j,k)
  do k = 1, nz
    do j = 1, ny
      g(1,1,j,k)  = g_post(1,nx,j,k)
      g(2,nx,j,k) = g_post(2,1,j,k)
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef SpanwiseWallsAdiabatic
  !$omp parallel do collapse(2) default(none) shared(g,g_post) private(j,k)
  do k = 1, nz
    do j = 1, ny
      g(1,1,j,k)  = g_post(2,1,j,k)
      g(2,nx,j,k) = g_post(1,nx,j,k)
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef VerticalWallsPeriodicalT
  !$omp parallel do collapse(2) default(none) shared(g,g_post) private(i,k)
  do k = 1, nz
    do i = 1, nx
      g(3,i,1,k)  = g_post(3,i,ny,k)
      g(4,i,ny,k) = g_post(4,i,1,k)
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef VerticalWallsConstT
  !$omp parallel do collapse(2) default(none) shared(g,g_post,omegaT) private(i,k)
  do k = 1, nz
    do i = 1, nx
#ifdef EnableLegacyThermalScheme
      g(3,i,1,k)  = -g_post(4,i,1,k)  + (6.0d0 + paraA) / 21.0d0 * Thot
      g(4,i,ny,k) = -g_post(3,i,ny,k) + (6.0d0 + paraA) / 21.0d0 * Tcold
#else
      g(3,i,1,k)  = -g_post(4,i,1,k)  + 2.0d0 * omegaT(3) * Thot
      g(4,i,ny,k) = -g_post(3,i,ny,k) + 2.0d0 * omegaT(4) * Tcold
#endif
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef VerticalWallsAdiabatic
  !$omp parallel do collapse(2) default(none) shared(g,g_post) private(i,k)
  do k = 1, nz
    do i = 1, nx
      g(3,i,1,k)  = g_post(4,i,1,k)
      g(4,i,ny,k) = g_post(3,i,ny,k)
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef HorizontalWallsConstT
  !$omp parallel do collapse(2) default(none) shared(g,g_post,omegaT) private(i,j)
  do j = 1, ny
    do i = 1, nx
#ifdef EnableLegacyThermalScheme
      g(5,i,j,1)  = -g_post(6,i,j,1)  + (6.0d0 + paraA) / 21.0d0 * Thot
      g(6,i,j,nz) = -g_post(5,i,j,nz) + (6.0d0 + paraA) / 21.0d0 * Tcold
#else
      g(5,i,j,1)  = -g_post(6,i,j,1)  + 2.0d0 * omegaT(5) * Thot
      g(6,i,j,nz) = -g_post(5,i,j,nz) + 2.0d0 * omegaT(6) * Tcold
#endif
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef HorizontalWallsPeriodicalT
  !$omp parallel do collapse(2) default(none) shared(g,g_post) private(i,j)
  do j = 1, ny
    do i = 1, nx
      g(5,i,j,1)  = g_post(5,i,j,nz)
      g(6,i,j,nz) = g_post(6,i,j,1)
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef HorizontalWallsAdiabatic
  !$omp parallel do collapse(2) default(none) shared(g,g_post) private(i,j)
  do j = 1, ny
    do i = 1, nx
      g(5,i,j,1)  = g_post(6,i,j,1)
      g(6,i,j,nz) = g_post(5,i,j,nz)
    enddo
  enddo
  !$omp end parallel do
#endif

  return
end subroutine bouncebackT
!===========================================================================================================================
! bouncebackT 结束: 处理温度边界条件，包括恒温、绝热和周期边界。
!===========================================================================================================================


!===========================================================================================================================
! 子程序: macroT
! 作用: 由温度分布函数恢复温度场，并更新历史热流项。
! 用途: 在主程序时间推进循环中调用，作为温度更新链条的最后一步。
!===========================================================================================================================
subroutine macroT()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k

  ! 温度恢复就是对 7 个方向的 g 求和
  !$omp parallel do collapse(3) default(none) shared(g,T) private(i,j,k)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        T(i,j,k) = g(0,i,j,k) + g(1,i,j,k) + g(2,i,j,k) + g(3,i,j,k) + &
                   g(4,i,j,k) + g(5,i,j,k) + g(6,i,j,k)
      enddo
    enddo
  enddo
  !$omp end parallel do

  return
end subroutine macroT
!===========================================================================================================================
! macroT 结束: 由温度分布函数恢复宏观温度场 T。
!===========================================================================================================================


!===========================================================================================================================
! 子程序: reconstruct_macro_from_fg
! 作用: 从重启读回的 f/g 重新恢复宏观场，避免备份文件格式过重。
! 用途: 在 loadInitField=1 的重启路径中调用，用于从严格重启文件恢复宏观量。
!===========================================================================================================================
subroutine reconstruct_macro_from_fg()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k, iter
  real(kind=8) :: momx, momy, momz
  real(kind=8) :: FxLoc, FyLoc, FzLoc
  logical :: rho_bad

  ! 重启时从 f/g 重构 T、rho、u、v、w；EnableUseG 的 B 历史项由严格重启文件单独读回。
  call macroT()
  rho_bad = .false.

  !$omp parallel do collapse(3) default(none) shared(f,rho,u,v,w,T) &
  !$omp private(i,j,k,iter,momx,momy,momz,FxLoc,FyLoc,FzLoc) reduction(.or.:rho_bad)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        rho(i,j,k) = f(0,i,j,k) + f(1,i,j,k) + f(2,i,j,k) + f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + &
             f(6,i,j,k) + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + &
             f(12,i,j,k) + f(13,i,j,k) + f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k)

        momx = f(1,i,j,k) - f(2,i,j,k) + f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) - f(10,i,j,k) + &
             f(11,i,j,k) - f(12,i,j,k) + f(13,i,j,k) - f(14,i,j,k)

        momy = f(3,i,j,k) - f(4,i,j,k) + f(7,i,j,k) + f(8,i,j,k) - f(9,i,j,k) - f(10,i,j,k) + &
             f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) - f(18,i,j,k)

        momz = f(5,i,j,k) - f(6,i,j,k) + f(11,i,j,k) + f(12,i,j,k) - f(13,i,j,k) - f(14,i,j,k) + &
             f(15,i,j,k) + f(16,i,j,k) - f(17,i,j,k) - f(18,i,j,k)

        if (rho(i,j,k) .GT. 0.0d0) then
          u(i,j,k) = momx / rho(i,j,k)
          v(i,j,k) = momy / rho(i,j,k)
          w(i,j,k) = momz / rho(i,j,k)
          ! 力项里含有 rho 和 T，因此这里做几次简单迭代来恢复带半步修正的速度
          do iter = 1, 3
            FxLoc = 0.0d0
            FyLoc = 0.0d0
            FzLoc = rho(i,j,k) * gBeta * (T(i,j,k) - Tref)
            u(i,j,k) = (momx + 0.5d0 * FxLoc) / rho(i,j,k)
            v(i,j,k) = (momy + 0.5d0 * FyLoc) / rho(i,j,k)
            w(i,j,k) = (momz + 0.5d0 * FzLoc) / rho(i,j,k)
          enddo
        else
          rho_bad = .true.
          u(i,j,k) = 0.0d0
          v(i,j,k) = 0.0d0
          w(i,j,k) = 0.0d0
        endif
      enddo
    enddo
  enddo
  !$omp end parallel do

  if (rho_bad) then
    write(*,*) 'Warning: non-positive rho found during restart reconstruction.'
    stop
  endif

  return
end subroutine reconstruct_macro_from_fg
!===========================================================================================================================
! reconstruct_macro_from_fg end: current macro state is rebuilt; EnableUseG history is read separately
!===========================================================================================================================


#ifdef steadyFlow
!===========================================================================================================================
! 子程序: check
! 作用: 计算稳态收敛误差，并按需写入收敛历史。
! 用途: 在 steadyFlow 模式下由主程序定期调用。
!===========================================================================================================================
subroutine check()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k
  real(kind=8) :: error1, error2, error5, error6
  character(len=80) :: caseTag

  ! 误差定义：errorU 用速度场相对 L2 误差
  ! errorT 用温度场相对 L1 误差
  error1 = 0.0d0
  error2 = 0.0d0
  error5 = 0.0d0
  error6 = 0.0d0

  !$omp parallel do collapse(3) default(none) &
  !$omp& shared(u,up,v,vp,w,wp,T,Tp) private(i,j,k) &
  !$omp& reduction(+:error1,error2,error5,error6)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        error1 = error1 + (u(i,j,k)-up(i,j,k))*(u(i,j,k)-up(i,j,k)) &
                        + (v(i,j,k)-vp(i,j,k))*(v(i,j,k)-vp(i,j,k)) &
                        + (w(i,j,k)-wp(i,j,k))*(w(i,j,k)-wp(i,j,k))
        error2 = error2 + u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)
        error5 = error5 + dabs(T(i,j,k)-Tp(i,j,k))
        error6 = error6 + dabs(T(i,j,k))
        up(i,j,k) = u(i,j,k)
        vp(i,j,k) = v(i,j,k)
        wp(i,j,k) = w(i,j,k)
        Tp(i,j,k) = T(i,j,k)
      enddo
    enddo
  enddo
  !$omp end parallel do

  if (error2 .GT. 1.0d-30) then
    errorU = dsqrt(error1) / dsqrt(error2)
  else
    errorU = dsqrt(error1)
  endif
  if (error6 .GT. 1.0d-30) then
    errorT = error5 / error6
  else
    errorT = error5
  endif

  call append_convergence_tecplot('convergence3D.plt', restartItcOffset+itc, errorU, errorT)
  write(caseTag,'("Ra=",ES10.3E2,",nx=",I0,",ny=",I0,",nz=",I0,",useG=",L1,",old=",L1)') &
       Rayleigh, nx, ny, nz, useG, useLegacyThermalScheme
  call append_convergence_master_tecplot('convergence_all_3D.plt', caseTag, restartItcOffset+itc, errorU, errorT)
  write(*,'(I12,1X,ES24.16E3,1X,ES24.16E3)') restartItcOffset+itc, errorU, errorT

end subroutine check
#endif


!===========================================================================================================================
! 子程序: append_convergence_tecplot
! 作用: 向单个收敛历史文件追加一条误差记录。
!===========================================================================================================================
subroutine append_convergence_tecplot(filename, itcLoc, errorULoc, errorTLoc)
  use commondata3d, only: loadInitField
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4), intent(in) :: itcLoc
  real(kind=8), intent(in) :: errorULoc, errorTLoc
  integer(kind=4) :: u
  logical :: ex
  logical, save :: first_write = .true.

  if (first_write) then
    inquire(file=trim(filename), exist=ex)
    if ((loadInitField .EQ. 1) .AND. ex) then
      ! 续算：旧收敛曲线继续追加，横坐标已经传入累计 itc。
      open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
    elseif (loadInitField .EQ. 1) then
      ! 续算要求旧收敛文件必须存在；否则会丢掉断电前的收敛历史。
      write(*,*) 'Error: restart requested but convergence file is missing: ', trim(filename)
      stop
    else
      ! 新算例：清掉旧历史，避免不同算例的数据混在一起。
      open(newunit=u, file=trim(filename), status='replace', action='write', form='formatted')
      write(u,'(A)') 'VARIABLES = "itc" "errorU" "errorT"'
      write(u,'(A)') 'ZONE T="conv3d", F=POINT'
    endif
    write(u,'(I12,1X,ES24.16E3,1X,ES24.16E3)') itcLoc, errorULoc, errorTLoc
    close(u)
    first_write = .false.
  else
    ! 同一次运行的后续调用：追加数据行
    open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
    write(u,'(I12,1X,ES24.16E3,1X,ES24.16E3)') itcLoc, errorULoc, errorTLoc
    close(u)
  endif

end subroutine append_convergence_tecplot


!===========================================================================================================================
! 子程序: append_convergence_master_tecplot
! 作用: 向带 zone 名称的收敛历史文件追加一条记录。
!===========================================================================================================================
subroutine append_convergence_master_tecplot(filename, zoneName, itcLoc, errorULoc, errorTLoc)
  use commondata3d, only: loadInitField
  implicit none
  character(len=*), intent(in) :: filename, zoneName
  integer(kind=4), intent(in) :: itcLoc
  real(kind=8), intent(in) :: errorULoc, errorTLoc
  logical :: ex
  integer(kind=4) :: u
  logical, save :: zone_started = .false.

  if (.not. zone_started) then
    inquire(file=trim(filename), exist=ex)
    if (.not. ex) then
      if (loadInitField .EQ. 1) then
        write(*,*) 'Error: restart requested but master convergence file is missing: ', trim(filename)
        stop
      endif
      open(newunit=u, file=trim(filename), status='new', action='write', form='formatted')
      write(u,'(A)') 'TITLE = "Convergence comparison 3D"'
      write(u,'(A)') 'VARIABLES = "itc" "errorU" "errorT"'
      close(u)
    endif
    if (loadInitField .EQ. 0) then
      open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
      write(u,'(A)') 'ZONE T="'//trim(zoneName)//'", F=POINT'
      close(u)
    endif
    zone_started = .true.
  endif

  open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
  write(u,'(I12,1X,ES24.16E3,1X,ES24.16E3)') itcLoc, errorULoc, errorTLoc
  close(u)

end subroutine append_convergence_master_tecplot


!===========================================================================================================================
! 子程序: output_SnapshotFile
! 作用: 输出三维快照二进制文件，供后处理或继续分析使用。
! 用途: 在运行过程中按需调用，也在程序结束时调用。
!===========================================================================================================================
subroutine output_SnapshotFile()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k
  character(len=100) :: filename

  ! This snapshot is for post-processing only; u/v/w are written after nondimensionalization.
  ! For strict restart, keep using output_ReloadFile(), which preserves the lattice-state variables.
#ifdef steadyFlow
  write(filename,'(i12.12)') restartItcOffset+itc
#endif
#ifdef unsteadyFlow
  snapshotFileNum = snapshotFileNum + 1
  write(filename,'(i12.12)') snapshotFileNum
#endif

  filename = adjustl(filename)
  open(unit=03, file=trim(snapshotFilePrefix)//'-'//trim(filename)//'.bin', form='unformatted', access='sequential')
  write(03) (((real(velocityScaleCompare*u(i,j,k), kind=8), i=1,nx), j=1,ny), k=1,nz)
  write(03) (((real(velocityScaleCompare*v(i,j,k), kind=8), i=1,nx), j=1,ny), k=1,nz)
  write(03) (((real(velocityScaleCompare*w(i,j,k), kind=8), i=1,nx), j=1,ny), k=1,nz)
  write(03) (((real(T(i,j,k), kind=8), i=1,nx), j=1,ny), k=1,nz)
  write(03) (((real(rho(i,j,k), kind=8), i=1,nx), j=1,ny), k=1,nz)
  close(03)

  return
end subroutine output_SnapshotFile
!===========================================================================================================================
! output_SnapshotFile 结束: 输出 u、v、w、T、rho 的二进制快照文件。
!===========================================================================================================================


!===========================================================================================================================
! 子程序: output_ReloadFile
! 作用: 输出重启备份文件；基础记录为 f/g，EnableUseG 还包含 Bx_prev/By_prev/Bz_prev 历史项。
! 用途: 在运行过程中定期调用，也在程序结束前调用。
!===========================================================================================================================
subroutine output_ReloadFile()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k, alpha
  character(len=100) :: filename

#ifdef steadyFlow
  reloadFileNum = restartItcOffset+itc
  write(filename,'(i12.12)') reloadFileNum
#endif
#ifdef unsteadyFlow
  reloadFileNum = reloadFileNum + 1
  write(filename,'(i12.12)') reloadFileNum
#endif

  filename = adjustl(filename)
  open(unit=05, file=trim(reloadFilePrefix)//'-'//trim(filename)//'.bin', form='unformatted', access='sequential')
  ! Strict restart files store f/g. With EnableUseG, Bx_prev/By_prev/Bz_prev must also be saved;
  ! otherwise the first post-reload M1G correction would lose its previous-step history.
  write(05) ((((real(f(alpha,i,j,k), kind=8), i=1,nx), j=1,ny), k=1,nz), alpha=0,qf-1)
  write(05) ((((real(g(alpha,i,j,k), kind=8), i=1,nx), j=1,ny), k=1,nz), alpha=0,qt-1)
#ifdef EnableUseG
  write(05) (((real(Bx_prev(i,j,k), kind=8), i=1,nx), j=1,ny), k=1,nz)
  write(05) (((real(By_prev(i,j,k), kind=8), i=1,nx), j=1,ny), k=1,nz)
  write(05) (((real(Bz_prev(i,j,k), kind=8), i=1,nx), j=1,ny), k=1,nz)
#endif
  close(05)
  call write_reload_metadata(trim(filename))

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) 'Backup f/g restart state to: ', trim(reloadFilePrefix), '-', trim(filename), '.bin'
  write(00,*) 'Backup restart metadata to: ', trim(reloadFilePrefix), '-latest.meta'
  close(00)

  return
end subroutine output_ReloadFile
!===========================================================================================================================
! output_ReloadFile 结束: 输出严格重启备份文件。
!===========================================================================================================================


!===========================================================================================================================
! 子程序: write_reload_metadata
! 作用: 覆盖写出最新 reload 续算账本，恢复累计步数、t_ff、输出编号和最新 .bin 文件名。
!===========================================================================================================================
subroutine write_reload_metadata(filename)
  use commondata3d
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4) :: metaUnit, totalItc
  real(kind=8) :: totalTf

  totalItc = restartItcOffset + itc
  totalTf = real(totalItc,kind=8) / timeUnit

  open(newunit=metaUnit, file=trim(reloadFilePrefix)//'-latest.meta', &
       status='replace', action='write', form='formatted')
  write(metaUnit,'(A,1X,I0)') 'reload_meta_version', 3
#ifdef steadyFlow
  write(metaUnit,'(A,1X,A)') 'flowMode', 'steadyFlow'
#endif
#ifdef unsteadyFlow
  write(metaUnit,'(A,1X,A)') 'flowMode', 'unsteadyFlow'
#endif
  write(metaUnit,'(A,1X,I0)') 'nx', nx
  write(metaUnit,'(A,1X,I0)') 'ny', ny
  write(metaUnit,'(A,1X,I0)') 'nz', nz
  write(metaUnit,'(A,1X,A)') 'reloadFileName', trim(filename)
  write(metaUnit,'(A,1X,I0)') 'itc_total', totalItc
  write(metaUnit,'(A,1X,ES24.16E3)') 'time_tf', totalTf
  write(metaUnit,'(A,1X,I0)') 'snapshotFileNum', snapshotFileNum
  write(metaUnit,'(A,1X,I0)') 'pltFileNum', pltFileNum
  write(metaUnit,'(A,1X,I0)') 'reloadFileNum', reloadFileNum
  close(metaUnit)

  return
end subroutine write_reload_metadata
!===========================================================================================================================


!===========================================================================================================================
! 子程序: read_reload_metadata
! 作用: 优先读取 latest .meta；若没有，则根据手工编号做保守推断。
!===========================================================================================================================
subroutine read_reload_metadata(reloadFileName)
  use commondata3d
  implicit none
  character(len=*), intent(inout) :: reloadFileName
  character(len=64) :: label
  character(len=32) :: metaFlowMode, currentFlowMode
  character(len=100) :: metaReloadFileName
  character(len=256) :: metaFile
  integer(kind=4) :: metaUnit, ios
  integer(kind=4) :: metaVersion, metaNx, metaNy, metaNz
  integer(kind=4) :: metaItc, metaSnapshotFileNum, metaPltFileNum, metaReloadFileNum
  real(kind=8) :: metaTf
  logical :: metaExists

  reloadMetadataLoaded = .false.
  metaFile = trim(reloadFilePrefix)//'-latest.meta'
  inquire(file=trim(metaFile), exist=metaExists)                 ! 优先检查最新账本

  if (.not. metaExists) then                                     ! latest meta 不存在时，只能保守推断
    call infer_reload_offsets_without_metadata()
    return
  endif

  open(newunit=metaUnit, file=trim(metaFile), status='old', action='read', &
       form='formatted', iostat=ios)                             ! ios==0 表示成功，非 0 表示打开失败
  if (ios .NE. 0) then
    write(*,*) 'Error: failed to open reload metadata: ', trim(metaFile)
    stop
  endif

  read(metaUnit,*,iostat=ios) label, metaVersion
  if ((ios .NE. 0) .OR. (trim(label) .NE. 'reload_meta_version') .OR. (metaVersion .NE. 3)) then
    write(*,*) 'Error: invalid reload metadata version in ', trim(metaFile)
    stop
  endif

  read(metaUnit,*,iostat=ios) label, metaFlowMode
  if ((ios .NE. 0) .OR. (trim(label) .NE. 'flowMode')) then
    write(*,*) 'Error: invalid flowMode entry in ', trim(metaFile)
    stop
  endif

  read(metaUnit,*,iostat=ios) label, metaNx
  if ((ios .NE. 0) .OR. (trim(label) .NE. 'nx')) then
    write(*,*) 'Error: invalid nx entry in ', trim(metaFile)
    stop
  endif

  read(metaUnit,*,iostat=ios) label, metaNy
  if ((ios .NE. 0) .OR. (trim(label) .NE. 'ny')) then
    write(*,*) 'Error: invalid ny entry in ', trim(metaFile)
    stop
  endif

  read(metaUnit,*,iostat=ios) label, metaNz
  if ((ios .NE. 0) .OR. (trim(label) .NE. 'nz')) then
    write(*,*) 'Error: invalid nz entry in ', trim(metaFile)
    stop
  endif

  read(metaUnit,*,iostat=ios) label, metaReloadFileName
  if ((ios .NE. 0) .OR. (trim(label) .NE. 'reloadFileName')) then
    write(*,*) 'Error: invalid reloadFileName entry in ', trim(metaFile)
    stop
  endif
  metaReloadFileName = adjustl(metaReloadFileName)

  read(metaUnit,*,iostat=ios) label, metaItc
  if ((ios .NE. 0) .OR. (trim(label) .NE. 'itc_total')) then
    write(*,*) 'Error: invalid itc_total entry in ', trim(metaFile)
    stop
  endif

  read(metaUnit,*,iostat=ios) label, metaTf
  if ((ios .NE. 0) .OR. (trim(label) .NE. 'time_tf')) then
    write(*,*) 'Error: invalid time_tf entry in ', trim(metaFile)
    stop
  endif

  read(metaUnit,*,iostat=ios) label, metaSnapshotFileNum
  if ((ios .NE. 0) .OR. (trim(label) .NE. 'snapshotFileNum')) then
    write(*,*) 'Error: invalid snapshotFileNum entry in ', trim(metaFile)
    stop
  endif

  read(metaUnit,*,iostat=ios) label, metaPltFileNum
  if ((ios .NE. 0) .OR. (trim(label) .NE. 'pltFileNum')) then
    write(*,*) 'Error: invalid pltFileNum entry in ', trim(metaFile)
    stop
  endif

  read(metaUnit,*,iostat=ios) label, metaReloadFileNum
  if ((ios .NE. 0) .OR. (trim(label) .NE. 'reloadFileNum')) then
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

  if (trim(metaFlowMode) .NE. trim(currentFlowMode)) then
    write(*,*) 'Error: reload metadata flowMode differs: ', trim(metaFlowMode), trim(currentFlowMode)
    stop
  endif
  if ((metaNx .NE. nx) .OR. (metaNy .NE. ny) .OR. (metaNz .NE. nz)) then
    write(*,*) 'Error: reload metadata mesh mismatch: ', metaNx, metaNy, metaNz, nx, ny, nz
    stop
  endif

  restartItcOffset = metaItc
  reloadDimensionlessTime = metaTf
  snapshotFileNum = metaSnapshotFileNum
  pltFileNum = metaPltFileNum
  ! reloadFileNum 是整数计数器，给后续 output_ReloadFile() 继续编号，避免覆盖旧 reload 文件。
  ! reloadFileName 是字符串文件名，本次续算马上用它打开 <reloadFilePrefix>-<reloadFileName>.bin。
  reloadFileNum = metaReloadFileNum
  reloadFileName = trim(metaReloadFileName)
  reloadMetadataLoaded = .true.

  return
end subroutine read_reload_metadata
!===========================================================================================================================


!===========================================================================================================================
! 子程序: infer_reload_offsets_without_metadata
! 作用: 没有 latest .meta 时，只能根据文件编号和当前手工参数推断。
! 根据文件名编号和当前参数“猜一个合理值”，保证续算的时间/步数尽量连续。
!===========================================================================================================================
subroutine infer_reload_offsets_without_metadata()
  use commondata3d
  implicit none

  restartItcOffset = 0
#ifdef steadyFlow
  restartItcOffset = max(0, reloadFileNum)  ! 稳态的 reload 文件名本来就是用 itc 写的
  if (reloadDimensionlessTime .EQ. 0.0d0) then
    reloadDimensionlessTime = real(restartItcOffset,kind=8) / timeUnit
  endif
#endif
#ifdef unsteadyFlow
  if (reloadDimensionlessTime .EQ. 0.0d0) then
    reloadDimensionlessTime = real(max(0,reloadFileNum),kind=8) * reloadFileInterval
  endif
  restartItcOffset = max(0, int(reloadDimensionlessTime*timeUnit+0.5d0))
  snapshotFileNum = max(0, int(reloadDimensionlessTime/outputSnapshotInterval+0.5d0))
  pltFileNum = max(0, int(reloadDimensionlessTime/outputPltFileInterval+0.5d0))
#endif

  return
end subroutine infer_reload_offsets_without_metadata
!===========================================================================================================================


!===========================================================================================================================
! 子程序: output_Tecplot
! 作用: 输出全场体数据以及 x/y/z 三个中面切片，便于快速查看三维流场结构。
! 用途: 在运行过程中按需调用，也在程序结束时调用。
!===========================================================================================================================
subroutine output_Tecplot()
  use commondata3d
  implicit none

  character(len=100) :: filename

  ! 3D 这里同时输出整体 Tecplot 体数据和三个中面切片；
#ifdef steadyFlow
  write(filename,'(i12.12)') restartItcOffset+itc
#endif
#ifdef unsteadyFlow
  pltFileNum = pltFileNum + 1
  write(filename,'(i12.12)') pltFileNum
#endif

  filename = adjustl(filename)
  call write_full_fields_plt(trim(pltFolderPrefix)//'-fullField-'//trim(filename)//'.plt')
  call write_midplane_x(trim(pltFolderPrefix)//'-midX-'//trim(filename)//'.plt')
  call write_midplane_y(trim(pltFolderPrefix)//'-midY-'//trim(filename)//'.plt')
  call write_midplane_z(trim(pltFolderPrefix)//'-midZ-'//trim(filename)//'.plt')
  call write_midplane_stream_x(trim(pltFolderPrefix)//'-midXpsiVort-'//trim(filename)//'.plt')
  call write_midplane_stream_y(trim(pltFolderPrefix)//'-midYpsiVort-'//trim(filename)//'.plt')
  call write_midplane_stream_z(trim(pltFolderPrefix)//'-midZpsiVort-'//trim(filename)//'.plt')

  return
end subroutine output_Tecplot
!===========================================================================================================================
! output_Tecplot 结束: 输出主场变量到 Tecplot 切片文件，便于后处理和可视化。
!===========================================================================================================================


!===========================================================================================================================
! 子程序: calNuRe
! 作用: 计算体平均 Nu / Re，并把时间序列缓存到数组中。
! 用途: 在主程序时间推进过程中按输出间隔调用，也在程序结束阶段补记一次。
!===========================================================================================================================
subroutine calNuRe()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k
  real(kind=8) :: NuVolAvg_temp, ReVolAvg_temp
  real(kind=8) :: volumeWeight
  real(kind=8) :: sampleTime
  logical :: exNu, exRe
  logical, save :: first_nure_write = .true.

  ! 这里记录的是时间序列版本的体平均 Nu / Re：
  ! NuVolAvg : 体平均对流热通量对应的 Nu
  ! ReVolAvg : 全域速度模体平均对应的 Reynolds 数
  if (dimensionlessTime .GE. dimensionlessTimeMax) then
    write(*,*) 'Error: dimensionlessTime exceeds dimensionlessTimeMax, please enlarge dimensionlessTimeMax'
    open(unit=00, file=trim(settingsFile), status='unknown', position='append')
    write(00,*) 'Error: dimensionlessTime exceeds dimensionlessTimeMax, please enlarge dimensionlessTimeMax'
    close(00)
    stop
  endif
  dimensionlessTime = dimensionlessTime + 1
#ifdef steadyFlow
  sampleTime = real(restartItcOffset+itc,kind=8)
#endif
#ifdef unsteadyFlow
  sampleTime = reloadDimensionlessTime + real(dimensionlessTime,kind=8)*outputSnapshotInterval
#endif

  if ((first_nure_write) .AND. (loadInitField .EQ. 1)) then
    inquire(file='Nu_VolAvg_3D.dat', exist=exNu)
    inquire(file='Re_VolAvg_3D.dat', exist=exRe)
    if ((.not. exNu) .OR. (.not. exRe)) then
      write(*,*) 'Error: restart requested but old Nu/Re time-series files are missing.'
      open(unit=00, file=trim(settingsFile), status='unknown', position='append')
      write(00,*) 'Error: restart requested but old Nu/Re time-series files are missing.'
      write(00,*) 'Nu_VolAvg_3D.dat exists =', exNu
      write(00,*) 'Re_VolAvg_3D.dat exists =', exRe
      close(00)
      stop
    endif
  endif

  NuVolAvg_temp = 0.0d0
#ifdef SideHeatedCell
  !$omp parallel do collapse(3) default(none) shared(v,T,quadWx,quadWy,quadWz) private(i,j,k,volumeWeight) reduction(+:NuVolAvg_temp)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        volumeWeight = quadWx(i)*quadWy(j)*quadWz(k)
        NuVolAvg_temp = NuVolAvg_temp + volumeWeight*v(i,j,k)*(T(i,j,k)-Tref)
      enddo
    enddo
  enddo
  !$omp end parallel do
#else
  !$omp parallel do collapse(3) default(none) shared(w,T,quadWx,quadWy,quadWz) private(i,j,k,volumeWeight) reduction(+:NuVolAvg_temp)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        volumeWeight = quadWx(i)*quadWy(j)*quadWz(k)
        NuVolAvg_temp = NuVolAvg_temp + volumeWeight*w(i,j,k)*(T(i,j,k)-Tref)
      enddo
    enddo
  enddo
  !$omp end parallel do
#endif

  NuVolAvg(dimensionlessTime) = NuVolAvg_temp/quadSumVolume*lengthUnit/diffusivity + 1.0d0
  if ((first_nure_write) .AND. (loadInitField .EQ. 0)) then
    open(unit=01, file='Nu_VolAvg_3D.dat', status='replace', action='write')
  else
    open(unit=01, file='Nu_VolAvg_3D.dat', status='unknown', position='append', action='write')
  endif
  write(01,'(ES24.16E3,1X,ES24.16E3)') &
       real(sampleTime, kind=8), NuVolAvg(dimensionlessTime)
  close(01)

  ReVolAvg_temp = 0.0d0
  !$omp parallel do collapse(3) default(none) shared(u,v,w,quadWx,quadWy,quadWz) private(i,j,k,volumeWeight) reduction(+:ReVolAvg_temp)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        volumeWeight = quadWx(i)*quadWy(j)*quadWz(k)
        ReVolAvg_temp = ReVolAvg_temp + volumeWeight*dsqrt(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k))
      enddo
    enddo
  enddo
  !$omp end parallel do
  ReVolAvg(dimensionlessTime) = ReVolAvg_temp/quadSumVolume*lengthUnit/viscosity
  if ((first_nure_write) .AND. (loadInitField .EQ. 0)) then
    open(unit=02, file='Re_VolAvg_3D.dat', status='replace', action='write')
  else
    open(unit=02, file='Re_VolAvg_3D.dat', status='unknown', position='append', action='write')
  endif
  write(02,'(ES24.16E3,1X,ES24.16E3)') &
       real(sampleTime, kind=8), ReVolAvg(dimensionlessTime)
  close(02)
  first_nure_write = .false.

  write(*,'(a,1x,ES24.16E3)') 'NuVolAvg =', NuVolAvg(dimensionlessTime)
  write(*,'(a,1x,ES24.16E3)') 'ReVolAvg =', ReVolAvg(dimensionlessTime)

  return
end subroutine calNuRe
!===========================================================================================================================
! calNuRe 结束: 计算体平均 Nu 和 Re 的时间历程统计量。
!===========================================================================================================================


#ifdef unsteadyFlow
!===========================================================================================================================
! Subroutine: output_unsteady_NuRe_postprocess
! Purpose: rebuild unsteady Nu/Re series, running means, and window averages from full .dat history.
!===========================================================================================================================
subroutine output_unsteady_NuRe_postprocess()
  use commondata3d
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

  inquire(file='Nu_VolAvg_3D.dat', exist=exNu)
  inquire(file='Re_VolAvg_3D.dat', exist=exRe)
  if ((.not. exNu) .OR. (.not. exRe)) then
    write(*,'(A)') 'Error: Nu/Re history files are missing before postprocessing.'
    open(unit=00, file=trim(settingsFile), status='unknown', position='append')
    write(00,'(A)') 'Error: Nu/Re history files are missing before postprocessing.'
    close(00)
    error stop 1
  endif

  open(newunit=nuUnit, file='Nu_VolAvg_3D.dat', status='old', action='read', form='formatted')
  open(newunit=reUnit, file='Re_VolAvg_3D.dat', status='old', action='read', form='formatted')

  ! These files are derived views of the full .dat history, so rebuild one continuous ZONE.
  open(newunit=seriesUnit, file='NuRe_VolAvg_3DOpenmp.plt', status='replace', action='write', form='formatted')
  write(seriesUnit,'(A)') 'TITLE = "3D OpenMP Nu/Re volume averages"'
  write(seriesUnit,'(A)') 'VARIABLES = "time" "NuVolAvg" "ReVolAvg"'
  write(seriesUnit,'(A)') 'ZONE T="NuReVolAvg", F=POINT'

  open(newunit=runningUnit, file='NuRe_VolAvg_runningMean_3DOpenmp.plt', status='replace', action='write', &
       form='formatted')
  write(runningUnit,'(A)') 'TITLE = "3D OpenMP Nu/Re running means"'
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
    if ((iosNu .NE. 0) .OR. (iosRe .NE. 0)) then
      write(*,'(A)') 'Error: Nu/Re history files are shorter than unsteadySampleCount or contain invalid rows.'
      open(unit=00, file=trim(settingsFile), status='unknown', position='append')
      write(00,'(A)') 'Error: Nu/Re history files are shorter than unsteadySampleCount or contain invalid rows.'
      close(00)
      error stop 1
    endif
    ! 确保 Nu 和 Re 是同一个时间采样点的数据，不是错行配对的数据。
    if (abs(timeNu-timeRe) .GT. 1.0d-10*max(1.0d0,abs(timeNu))) then
      write(*,'(A)') 'Error: Nu/Re history time columns do not match.'
      open(unit=00, file=trim(settingsFile), status='unknown', position='append')
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
  if ((iosNu .EQ. 0) .OR. (iosRe .EQ. 0)) then
    write(*,'(A)') 'Error: Nu/Re history files contain more rows than unsteadySampleCount.'
    open(unit=00, file=trim(settingsFile), status='unknown', position='append')
    write(00,'(A)') 'Error: Nu/Re history files contain more rows than unsteadySampleCount.'
    close(00)
    error stop 1
  endif
  ! 正常结尾必须是两个文件都到达 EOF，也就是两个 iostat 都小于 0。
  ! 其他组合说明一个文件尾部异常，或 Nu/Re 文件长度不一致。
  if (.not. ((iosNu .LT. 0) .AND. (iosRe .LT. 0))) then
    write(*,'(A)') 'Error: Nu/Re history files have inconsistent trailing rows.'
    open(unit=00, file=trim(settingsFile), status='unknown', position='append')
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
    open(unit=00, file=trim(settingsFile), status='unknown', position='append')
    write(00,'(A)') 'Error: no Nu/Re history samples were found before postprocessing.'
    close(00)
    error stop 1
  endif

  if ((whole_count <= 0) .or. (first_count <= 0) .or. (second_count <= 0)) then
    write(*,'(A)') 'Error: no complete unsteady average window was found for Nu/Re postprocessing.'
    open(unit=00, file=trim(settingsFile), status='unknown', position='append')
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

  open(unit=33, file='NuRe_TimeAverage_3DOpenmp.txt', status='replace', action='write', form='formatted')
  write(33,'(A)') '# 3D OpenMP Nu/Re statistical-convergence window averages'
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
!===========================================================================================================================
#endif


#ifdef steadyFlow
!===========================================================================================================================
! 子程序: SideHeatedcalc_Nu_global
! 作用: 计算侧壁差温工况下的全场平均 Nusselt 数。
! 用途: 在 SideHeatedCell 工况结束后的后处理中调用。
!===========================================================================================================================
subroutine SideHeatedcalc_Nu_global()
  use commondata3d
  implicit none
  integer(kind=4) :: i, j, k
  real(kind=8) :: dy, dTdy, qy, sum_qy
  real(kind=8) :: deltaT, coef

  dy = 1.0d0 / lengthUnit
  deltaT = Thot - Tcold
  coef = velocityScaleCompare
  sum_qy = 0.0d0

  !$omp parallel do collapse(3) default(none) shared(v,T,dy,coef) private(i,j,k,dTdy,qy) reduction(+:sum_qy)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        if (j .EQ. 1) then
          dTdy = (-3.0d0*T(i,1,k) - T(i,2,k) + 4.0d0*Thot) / (3.0d0*dy)
        elseif (j .EQ. ny) then
          dTdy = (-4.0d0*Tcold + 3.0d0*T(i,ny,k) + T(i,ny-1,k)) / (3.0d0*dy)
        else
          dTdy = (T(i,j-1,k) - T(i,j+1,k)) / (2.0d0*dy)
        endif

        qy = coef * v(i,j,k) * (T(i,j,k) - Tref) + dTdy
        sum_qy = sum_qy + qy
      enddo
    enddo
  enddo
  !$omp end parallel do

  Nu_global = (sum_qy / dble(nx * ny * nz)) / deltaT

  write(*,'(a,1x,ES24.16E3)') 'Nu_global =', Nu_global
  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,'(a,1x,ES24.16E3)') 'Nu_global =', Nu_global
  close(00)

  return
end subroutine SideHeatedcalc_Nu_global
!===========================================================================================================================
! SideHeatedcalc_Nu_global 结束: 计算侧壁差温工况下的全场平均 Nusselt 数。
!===========================================================================================================================
#endif


#ifdef steadyFlow
!===========================================================================================================================
! 子程序: RBcalc_Nu_global
! 作用: 计算 Rayleigh-Benard 工况下的全场平均 Nusselt 数。
! 用途: 在 RayleighBenardCell 工况结束后的后处理中调用。
!===========================================================================================================================
subroutine RBcalc_Nu_global()
  use commondata3d
  implicit none
  integer(kind=4) :: i, j, k
  real(kind=8) :: dz, dTdz, qz, sum_qz
  real(kind=8) :: deltaT, coef

  dz = 1.0d0 / lengthUnit
  deltaT = Thot - Tcold
  coef = velocityScaleCompare
  sum_qz = 0.0d0

  !$omp parallel do collapse(3) default(none) shared(w,T,dz,coef) private(i,j,k,dTdz,qz) reduction(+:sum_qz)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        if (k .EQ. 1) then
          dTdz = (3.0d0*T(i,j,1) + T(i,j,2) - 4.0d0*Thot) / (3.0d0*dz)
        elseif (k .EQ. nz) then
          dTdz = (4.0d0*Tcold - 3.0d0*T(i,j,nz) - T(i,j,nz-1)) / (3.0d0*dz)
        else
          dTdz = (T(i,j,k+1) - T(i,j,k-1)) / (2.0d0*dz)
        endif

        qz = coef * w(i,j,k) * (T(i,j,k) - Tref) - dTdz
        sum_qz = sum_qz + qz
      enddo
    enddo
  enddo
  !$omp end parallel do

  Nu_global = (sum_qz / dble(nx * ny * nz)) / deltaT

  write(*,'(a,1x,ES24.16E3)') 'Nu_global =', Nu_global
  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,'(a,1x,ES24.16E3)') 'Nu_global =', Nu_global
  close(00)

  return
end subroutine RBcalc_Nu_global
!===========================================================================================================================
! RBcalc_Nu_global 结束: 计算 Rayleigh-Benard 工况下的全场平均 Nusselt 数。
!===========================================================================================================================
#endif


#ifdef steadyFlow
!===========================================================================================================================
! 子程序: SideHeatedcalc_Nu_wall_avg
! 作用: 计算侧壁差温工况下热壁、冷壁和中面平均 Nusselt 数及其极值。沿 z/k 平均后的热壁局部 Nu 在 x/i 方向上的最大位置
! 用途: 在 SideHeatedCell 工况结束后的后处理中调用。
!===========================================================================================================================
subroutine SideHeatedcalc_Nu_wall_avg()
  use commondata3d
  implicit none
  integer(kind=4) :: i, imax, imin, jL, jR, jMid, k, m
  integer(kind=4) :: ii(5)
  real(kind=8) :: dy, deltaT, coef
  real(kind=8) :: qy_wall, sum_hot, sum_cold, sum_mid
  real(kind=8) :: T_wf, T_wb
  real(kind=8) :: xfit(4), Tfit(4)
  real(kind=8) :: xk(5), fk(5), fstar, xstar
  real(kind=8) :: Nu_hot_line(1:nx), Nu_hot_ext(0:nx+1), T_hot_avg(1:nx)

  dy = 1.0d0 / lengthUnit
  deltaT = Thot - Tcold
  coef = velocityScaleCompare

  !$omp parallel do default(none) shared(T,T_hot_avg) private(i,k)
  do i = 1, nx
    T_hot_avg(i) = 0.0d0
    do k = 1, nz
      T_hot_avg(i) = T_hot_avg(i) + T(i,1,k)
    enddo
    T_hot_avg(i) = T_hot_avg(i) / dble(nz)
  enddo
  !$omp end parallel do

  sum_hot = 0.0d0
  !$omp parallel do default(none) shared(T,Nu_hot_line,dy,deltaT) private(i,k,qy_wall) reduction(+:sum_hot)
  do i = 1, nx
    Nu_hot_line(i) = 0.0d0
    do k = 1, nz
      qy_wall = (8.0d0*Thot - 9.0d0*T(i,1,k) + T(i,2,k)) / (3.0d0*dy)
      Nu_hot_line(i) = Nu_hot_line(i) + qy_wall / deltaT
    enddo
    Nu_hot_line(i) = Nu_hot_line(i) / dble(nz)
    sum_hot = sum_hot + Nu_hot_line(i)
  enddo
  !$omp end parallel do
  Nu_hot = sum_hot / dble(nx)

  Nu_hot_ext(1:nx) = Nu_hot_line(1:nx)
  xfit(1) = xp(1);  Tfit(1) = T_hot_avg(1)
  xfit(2) = xp(2);  Tfit(2) = T_hot_avg(2)
  xfit(3) = xp(3);  Tfit(3) = T_hot_avg(3)
  xfit(4) = xp(4);  Tfit(4) = T_hot_avg(4)
  call fit_adiabatic_wall_T4(0.0d0, xfit, Tfit, T_wf)
  Nu_hot_ext(0) = (2.0d0 * (Thot - T_wf) / dy) / deltaT

  xfit(1) = xp(nx-3);  Tfit(1) = T_hot_avg(nx-3)
  xfit(2) = xp(nx-2);  Tfit(2) = T_hot_avg(nx-2)
  xfit(3) = xp(nx-1);  Tfit(3) = T_hot_avg(nx-1)
  xfit(4) = xp(nx  );  Tfit(4) = T_hot_avg(nx  )
  call fit_adiabatic_wall_T4(xp(nx+1), xfit, Tfit, T_wb)
  Nu_hot_ext(nx+1) = (2.0d0 * (Thot - T_wb) / dy) / deltaT

  imax = 0
  imin = 0
  Nu_hot_max = Nu_hot_ext(0)
  Nu_hot_min = Nu_hot_ext(0)
  do i = 1, nx+1
    if (Nu_hot_ext(i) .GT. Nu_hot_max) then
      Nu_hot_max = Nu_hot_ext(i)
      imax = i
    endif
    if (Nu_hot_ext(i) .LT. Nu_hot_min) then
      Nu_hot_min = Nu_hot_ext(i)
      imin = i
    endif
  enddo

  if (imax .LE. 2) then
    ii = (/ 0, 1, 2, 3, 4 /)
  elseif (imax .GE. nx-1) then
    ii = (/ nx-3, nx-2, nx-1, nx, nx+1 /)
  else
    ii = (/ imax-2, imax-1, imax, imax+1, imax+2 /)
  endif
  do m = 1, 5
    xk(m) = xp(ii(m))
    fk(m) = Nu_hot_ext(ii(m))
  enddo
  call fit_parabola_ls5(xk, fk, +1, fstar, xstar)
  Nu_hot_max = fstar
  Nu_hot_max_position = xstar

  if (imin .LE. 2) then
    ii = (/ 0, 1, 2, 3, 4 /)
  elseif (imin .GE. nx-1) then
    ii = (/ nx-3, nx-2, nx-1, nx, nx+1 /)
  else
    ii = (/ imin-2, imin-1, imin, imin+1, imin+2 /)
  endif
  do m = 1, 5
    xk(m) = xp(ii(m))
    fk(m) = Nu_hot_ext(ii(m))
  enddo
  call fit_parabola_ls5(xk, fk, -1, fstar, xstar)
  Nu_hot_min = fstar
  Nu_hot_min_position = xstar

  sum_cold = 0.0d0
  !$omp parallel do collapse(2) default(none) shared(T,dy,deltaT) private(i,k,qy_wall) reduction(+:sum_cold)
  do k = 1, nz
    do i = 1, nx
      qy_wall = (-8.0d0*Tcold + 9.0d0*T(i,ny,k) - T(i,ny-1,k)) / (3.0d0*dy)
      sum_cold = sum_cold + qy_wall / deltaT
    enddo
  enddo
  !$omp end parallel do
  Nu_cold = sum_cold / dble(nx * nz)

  sum_mid = 0.0d0
  if (mod(ny,2) .EQ. 1) then
    jMid = (ny + 1) / 2
    !$omp parallel do collapse(2) default(none) shared(v,T,jMid,dy,deltaT,coef) private(i,k) reduction(+:sum_mid)
    do k = 1, nz
      do i = 1, nx
        sum_mid = sum_mid + (coef * v(i,jMid,k) * (T(i,jMid,k) - Tref) + &
             (T(i,jMid-1,k) - T(i,jMid+1,k)) / (2.0d0*dy)) / deltaT
      enddo
    enddo
    !$omp end parallel do
  else
    jL = ny / 2
    jR = jL + 1
    !$omp parallel do collapse(2) default(none) shared(v,T,jL,jR,dy,deltaT,coef) private(i,k) reduction(+:sum_mid)
    do k = 1, nz
      do i = 1, nx
        sum_mid = sum_mid + (coef * 0.5d0 * (v(i,jL,k) * (T(i,jL,k) - Tref) + &
             v(i,jR,k) * (T(i,jR,k) - Tref)) + (T(i,jL,k) - T(i,jR,k)) / dy) / deltaT
      enddo
    enddo
    !$omp end parallel do
  endif
  Nu_middle = sum_mid / dble(nx * nz)

  write(*,'(a,1x,ES24.16E3)') 'Nu_hot(left)  =', Nu_hot
  write(*,'(a,1x,ES24.16E3)') 'Nu_cold(right)=', Nu_cold
  write(*,'(a,1x,ES24.16E3)') 'Nu_middle     =', Nu_middle
  write(*,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') 'Nu_hot_max =', Nu_hot_max, 'x_max =', Nu_hot_max_position
  write(*,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') 'Nu_hot_min =', Nu_hot_min, 'x_min =', Nu_hot_min_position

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,'(a,1x,ES24.16E3)') 'Nu_hot(left)  =', Nu_hot
  write(00,'(a,1x,ES24.16E3)') 'Nu_cold(right)=', Nu_cold
  write(00,'(a,1x,ES24.16E3)') 'Nu_middle     =', Nu_middle
  write(00,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') 'Nu_hot_max =', Nu_hot_max, 'x_max =', Nu_hot_max_position
  write(00,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') 'Nu_hot_min =', Nu_hot_min, 'x_min =', Nu_hot_min_position
  close(00)

  return
end subroutine SideHeatedcalc_Nu_wall_avg
!===========================================================================================================================
! SideHeatedcalc_Nu_wall_avg 结束: 计算侧壁差温工况下热壁、冷壁和中面的 Nusselt 数及其极值。
!===========================================================================================================================
#endif


!===========================================================================================================================
! 子程序: SideHeatedcalc_Nu_zmid_wall_mean
! 作用: 计算侧壁差温工况下 z=Lz/2 中平面上 y=0 和 y=1 两条壁线各自的平均 Nusselt 数。
! 用途: 保持与主侧壁 Nu 统计相同的壁面导数离散方式，只在 z 方向插值到中平面后沿 x 求平均。
!===========================================================================================================================
subroutine SideHeatedcalc_Nu_zmid_wall_mean()
  use commondata3d
  implicit none
  integer(kind=4) :: i, kL, kR
  real(kind=8) :: targetZ, weight
  real(kind=8) :: dy, deltaT
  real(kind=8) :: Tleft1, Tleft2, Tright1, Tright2
  real(kind=8) :: qy_hot, qy_cold, Nu_left_mean, Nu_right_mean

  targetZ = 0.5d0 * zp(nz+1)
  call find_bracketing_index(zp, nz, targetZ, kL, kR, weight)

  dy = 1.0d0 / lengthUnit
  deltaT = Thot - Tcold
  Nu_left_mean = 0.0d0
  Nu_right_mean = 0.0d0

  !$omp parallel do default(none) shared(T,kL,kR,weight,dy,deltaT) private(i,Tleft1,Tleft2,Tright1,Tright2,qy_hot,qy_cold) &
  !$omp reduction(+:Nu_left_mean,Nu_right_mean)
  do i = 1, nx
    call interp_scalar_z(kL, kR, weight, i, 1,    T, Tleft1)
    call interp_scalar_z(kL, kR, weight, i, 2,    T, Tleft2)
    call interp_scalar_z(kL, kR, weight, i, ny,   T, Tright1)
    call interp_scalar_z(kL, kR, weight, i, ny-1, T, Tright2)

    qy_hot  = ( 8.0d0 * Thot  - 9.0d0 * Tleft1  + Tleft2 ) / (3.0d0 * dy)
    qy_cold = (-8.0d0 * Tcold + 9.0d0 * Tright1 - Tright2) / (3.0d0 * dy)

    Nu_left_mean  = Nu_left_mean  + qy_hot  / deltaT
    Nu_right_mean = Nu_right_mean + qy_cold / deltaT
  enddo
  !$omp end parallel do

  Nu_left_mean  = Nu_left_mean  / dble(nx)
  Nu_right_mean = Nu_right_mean / dble(nx)
  write(*,'(a,1x,ES24.16E3)') 'Nu_zmid_left   =', Nu_left_mean
  write(*,'(a,1x,ES24.16E3)') 'Nu_zmid_right  =', Nu_right_mean
  write(*,'(a,1x,ES24.16E3)') 'z_mid          =', targetZ

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,'(a,1x,ES24.16E3)') 'Nu_zmid_left   =', Nu_left_mean
  write(00,'(a,1x,ES24.16E3)') 'Nu_zmid_right  =', Nu_right_mean
  write(00,'(a,1x,ES24.16E3)') 'z_mid          =', targetZ
  close(00)

  return
end subroutine SideHeatedcalc_Nu_zmid_wall_mean
!===========================================================================================================================
! SideHeatedcalc_Nu_zmid_wall_mean 结束: 计算 z=Lz/2 截面上左右 y 壁线各自的平均 Nu。
!===========================================================================================================================


#ifdef steadyFlow
!===========================================================================================================================
! 子程序: RBcalc_Nu_wall_avg
! 作用: 计算 Rayleigh-Benard 工况下热壁、冷壁和中面的 Nusselt 数及其极值。沿 y/j 平均后的热壁局部 Nu 在 x/i 方向上的最大位置
! 用途: 在 RayleighBenardCell 工况结束后的后处理中调用。
!===========================================================================================================================
subroutine RBcalc_Nu_wall_avg()
  use commondata3d
  implicit none
  integer(kind=4) :: i, imax, imin, j, kB, kT, kMid, m
  integer(kind=4) :: ii(5)
  real(kind=8) :: dx, dz, deltaT, coef
  real(kind=8) :: qz_wall, sum_hot, sum_cold, sum_mid
  real(kind=8) :: T_wl, T_wr
  real(kind=8) :: xfit(4), Tfit(4)
  real(kind=8) :: xk(5), fk(5), fstar, xstar
  real(kind=8) :: Nu_bot(1:nx), Nu_bot_ext(0:nx+1), T_bot_avg(1:nx)

  dx = 1.0d0 / lengthUnit
  dz = 1.0d0 / lengthUnit
  deltaT = Thot - Tcold
  coef = velocityScaleCompare

  !$omp parallel do default(none) shared(T,T_bot_avg) private(i,j)
  do i = 1, nx
    T_bot_avg(i) = 0.0d0
    do j = 1, ny
      T_bot_avg(i) = T_bot_avg(i) + T(i,j,1)
    enddo
    T_bot_avg(i) = T_bot_avg(i) / dble(ny)
  enddo
  !$omp end parallel do

  sum_hot = 0.0d0
  !$omp parallel do default(none) shared(T,Nu_bot,dz,deltaT) private(i,j,qz_wall) reduction(+:sum_hot)
  do i = 1, nx
    Nu_bot(i) = 0.0d0
    do j = 1, ny
      qz_wall = (8.0d0*Thot - 9.0d0*T(i,j,1) + T(i,j,2)) / (3.0d0*dz)
      Nu_bot(i) = Nu_bot(i) + qz_wall / deltaT
    enddo
    Nu_bot(i) = Nu_bot(i) / dble(ny)
    sum_hot = sum_hot + Nu_bot(i)
  enddo
  !$omp end parallel do
  Nu_hot = sum_hot / dble(nx)

  Nu_bot_ext(1:nx) = Nu_bot(1:nx)
  xfit(1) = xp(1);  Tfit(1) = T_bot_avg(1)
  xfit(2) = xp(2);  Tfit(2) = T_bot_avg(2)
  xfit(3) = xp(3);  Tfit(3) = T_bot_avg(3)
  xfit(4) = xp(4);  Tfit(4) = T_bot_avg(4)
  call fit_adiabatic_wall_T4(0.0d0, xfit, Tfit, T_wl)
  Nu_bot_ext(0) = (2.0d0 * (Thot - T_wl) / dz) / deltaT

  xfit(1) = xp(nx-3);  Tfit(1) = T_bot_avg(nx-3)
  xfit(2) = xp(nx-2);  Tfit(2) = T_bot_avg(nx-2)
  xfit(3) = xp(nx-1);  Tfit(3) = T_bot_avg(nx-1)
  xfit(4) = xp(nx  );  Tfit(4) = T_bot_avg(nx  )
  call fit_adiabatic_wall_T4(xp(nx+1), xfit, Tfit, T_wr)
  Nu_bot_ext(nx+1) = (2.0d0 * (Thot - T_wr) / dz) / deltaT

  imax = 0
  imin = 0
  Nu_hot_max = Nu_bot_ext(0)
  Nu_hot_min = Nu_bot_ext(0)
  do i = 1, nx+1
    if (Nu_bot_ext(i) .GT. Nu_hot_max) then
      Nu_hot_max = Nu_bot_ext(i)
      imax = i
    endif
    if (Nu_bot_ext(i) .LT. Nu_hot_min) then
      Nu_hot_min = Nu_bot_ext(i)
      imin = i
    endif
  enddo

  if (imax .LE. 2) then
    ii = (/ 0, 1, 2, 3, 4 /)
  elseif (imax .GE. nx-1) then
    ii = (/ nx-3, nx-2, nx-1, nx, nx+1 /)
  else
    ii = (/ imax-2, imax-1, imax, imax+1, imax+2 /)
  endif
  do m = 1, 5
    xk(m) = xp(ii(m))
    fk(m) = Nu_bot_ext(ii(m))
  enddo
  call fit_parabola_ls5(xk, fk, +1, fstar, xstar)
  Nu_hot_max = fstar
  Nu_hot_max_position = xstar

  if (imin .LE. 2) then
    ii = (/ 0, 1, 2, 3, 4 /)
  elseif (imin .GE. nx-1) then
    ii = (/ nx-3, nx-2, nx-1, nx, nx+1 /)
  else
    ii = (/ imin-2, imin-1, imin, imin+1, imin+2 /)
  endif
  do m = 1, 5
    xk(m) = xp(ii(m))
    fk(m) = Nu_bot_ext(ii(m))
  enddo
  call fit_parabola_ls5(xk, fk, -1, fstar, xstar)
  Nu_hot_min = fstar
  Nu_hot_min_position = xstar

  sum_cold = 0.0d0
  !$omp parallel do collapse(2) default(none) shared(T,dz,deltaT) private(i,j,qz_wall) reduction(+:sum_cold)
  do j = 1, ny
    do i = 1, nx
      qz_wall = (-8.0d0*Tcold + 9.0d0*T(i,j,nz) - T(i,j,nz-1)) / (3.0d0*dz)
      sum_cold = sum_cold + qz_wall / deltaT
    enddo
  enddo
  !$omp end parallel do
  Nu_cold = sum_cold / dble(nx * ny)

  sum_mid = 0.0d0
  if (mod(nz,2) .EQ. 1) then
    kMid = (nz + 1) / 2
    !$omp parallel do collapse(2) default(none) shared(w,T,kMid,dz,deltaT,coef) private(i,j) reduction(+:sum_mid)
    do j = 1, ny
      do i = 1, nx
        sum_mid = sum_mid + (coef * w(i,j,kMid) * (T(i,j,kMid) - Tref) - &
             (T(i,j,kMid+1) - T(i,j,kMid-1)) / (2.0d0*dz)) / deltaT
      enddo
    enddo
    !$omp end parallel do
  else
    kB = nz / 2
    kT = kB + 1
    !$omp parallel do collapse(2) default(none) shared(w,T,kB,kT,dz,deltaT,coef) private(i,j) reduction(+:sum_mid)
    do j = 1, ny
      do i = 1, nx
        sum_mid = sum_mid + (coef * 0.5d0 * (w(i,j,kB) * (T(i,j,kB) - Tref) + &
             w(i,j,kT) * (T(i,j,kT) - Tref)) + (T(i,j,kB) - T(i,j,kT)) / dz) / deltaT
      enddo
    enddo
    !$omp end parallel do
  endif
  Nu_middle = sum_mid / dble(nx * ny)

  write(*,'(a,1x,ES24.16E3)') 'Nu_hot(bottom) =', Nu_hot
  write(*,'(a,1x,ES24.16E3)') 'Nu_cold(top)   =', Nu_cold
  write(*,'(a,1x,ES24.16E3)') 'Nu_middle      =', Nu_middle
  write(*,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') 'Nu_hot_max =', Nu_hot_max, 'x_max =', Nu_hot_max_position
  write(*,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') 'Nu_hot_min =', Nu_hot_min, 'x_min =', Nu_hot_min_position

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,'(a,1x,ES24.16E3)') 'Nu_hot(bottom) =', Nu_hot
  write(00,'(a,1x,ES24.16E3)') 'Nu_cold(top)   =', Nu_cold
  write(00,'(a,1x,ES24.16E3)') 'Nu_middle      =', Nu_middle
  write(00,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') 'Nu_hot_max =', Nu_hot_max, 'x_max =', Nu_hot_max_position
  write(00,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') 'Nu_hot_min =', Nu_hot_min, 'x_min =', Nu_hot_min_position
  close(00)

  return
end subroutine RBcalc_Nu_wall_avg
!===========================================================================================================================
! RBcalc_Nu_wall_avg 结束: 计算 Rayleigh-Benard 工况下热壁、冷壁和中面的 Nusselt 数及其极值。
!===========================================================================================================================
#endif


!===========================================================================================================================
! 子程序: fit_adiabatic_wall_T4
! 作用: 用四点拟合估计绝热壁面的壁温。
! 用途: 在 SideHeated 和 RB 壁面 Nusselt 后处理中调用。
!===========================================================================================================================
subroutine fit_adiabatic_wall_T4(y0, y, tt, T_wall)
  implicit none
  real(kind=8), intent(in)  :: y0
  real(kind=8), intent(in)  :: y(4), tt(4)
  real(kind=8), intent(out) :: T_wall
  real(kind=8) :: s(4)
  real(kind=8) :: S0, S1, S2, B0, B1, D
  integer(kind=4) :: k

  do k = 1, 4
    s(k) = (y(k) - y0) * (y(k) - y0)
  enddo

  S0 = 4.0d0
  S1 = 0.0d0
  S2 = 0.0d0
  B0 = 0.0d0
  B1 = 0.0d0
  do k = 1, 4
    S1 = S1 + s(k)
    S2 = S2 + s(k) * s(k)
    B0 = B0 + tt(k)
    B1 = B1 + tt(k) * s(k)
  enddo

  D = S0 * S2 - S1 * S1
  T_wall = (B0 * S2 - B1 * S1) / D

  return
end subroutine fit_adiabatic_wall_T4
!===========================================================================================================================
! fit_adiabatic_wall_T4 结束: 用四点拟合估计绝热壁面的壁温。
!===========================================================================================================================


!===========================================================================================================================
! 子程序: fit_parabola_ls5
! 作用: 用五点最小二乘抛物线拟合局部极值和对应位置。
! 用途: 在 Nu 极值和中心面速度极值的后处理中重复调用。
!===========================================================================================================================
subroutine fit_parabola_ls5(y, f, mode, fstar, ystar)
  implicit none
  real(kind=8), intent(in)  :: y(5), f(5)
  integer(kind=4), intent(in) :: mode
  real(kind=8), intent(out) :: fstar, ystar
  real(kind=8) :: S0, S1, S2, S3, S4
  real(kind=8) :: F0, F1, F2
  real(kind=8) :: D, DA, DB, DC
  real(kind=8) :: A, B, C
  real(kind=8) :: ymin, ymax
  integer(kind=4) :: k, kbest
  real(kind=8), parameter :: epsD = 1.0d-20, epsA = 1.0d-14

  kbest = 1
  do k = 2, 5
    if (mode .EQ. 1) then
      if (f(k) .GT. f(kbest)) kbest = k
    else
      if (f(k) .LT. f(kbest)) kbest = k
    endif
  enddo

  S0 = 0.0d0
  S1 = 0.0d0
  S2 = 0.0d0
  S3 = 0.0d0
  S4 = 0.0d0
  F0 = 0.0d0
  F1 = 0.0d0
  F2 = 0.0d0
  do k = 1, 5
    S0 = S0 + 1.0d0
    S1 = S1 + y(k)
    S2 = S2 + y(k) * y(k)
    S3 = S3 + y(k) * y(k) * y(k)
    S4 = S4 + y(k) * y(k) * y(k) * y(k)

    F0 = F0 + f(k)
    F1 = F1 + f(k) * y(k)
    F2 = F2 + f(k) * y(k) * y(k)
  enddo

  D  = S0*(S2*S4 - S3*S3) - S1*(S1*S4 - S2*S3) + S2*(S1*S3 - S2*S2)
  DA = F0*(S2*S4 - S3*S3) - S1*(F1*S4 - S3*F2) + S2*(F1*S3 - S2*F2)
  DB = S0*(F1*S4 - S3*F2) - F0*(S1*S4 - S2*S3) + S2*(S1*F2 - F1*S2)
  DC = S0*(S2*F2 - F1*S3) - S1*(S1*F2 - F1*S2) + F0*(S1*S3 - S2*S2)

  if (dabs(D) .LE. epsD) then
    ystar = y(kbest)
    fstar = f(kbest)
    return
  endif

  A = DA / D
  B = DB / D
  C = DC / D

  ymin = minval(y)
  ymax = maxval(y)

  if (dabs(C) .LE. epsA) then
    ystar = y(kbest)
    fstar = f(kbest)
    return
  endif

  if ((mode .EQ. 1 .AND. C .GE. 0.0d0) .OR. (mode .EQ. -1 .AND. C .LE. 0.0d0)) then
    ystar = y(kbest)
    fstar = f(kbest)
    return
  endif

  ystar = -B / (2.0d0 * C)
  if ((ystar .LT. ymin) .OR. (ystar .GT. ymax)) then
    ystar = y(kbest)
    fstar = f(kbest)
    return
  endif

  fstar = A + B * ystar + C * ystar * ystar

  return
end subroutine fit_parabola_ls5
!===========================================================================================================================
! fit_parabola_ls5 结束: 用五点最小二乘抛物线拟合局部极值和对应位置。
!===========================================================================================================================


!===========================================================================================================================
! 子程序: SideHeatedcalc_umid_max
! 作用: 计算侧壁差温工况下 x=Lx/2 中面上的 u 最大值及其位置。
! 用途: 在 SideHeatedCell 工况结束后的后处理中调用。
!===========================================================================================================================
subroutine SideHeatedcalc_umid_max()
  call calc_umid_max_common('SideHeatedcalc_umid_max')
end subroutine SideHeatedcalc_umid_max


!===========================================================================================================================
! 子程序: SideHeatedcalc_vmid_max
! 作用: 计算侧壁差温工况下 y=Ly/2 中面上的 v 最大值及其位置。
!===========================================================================================================================
subroutine SideHeatedcalc_vmid_max()
  call calc_vmid_max_common('SideHeatedcalc_vmid_max')
end subroutine SideHeatedcalc_vmid_max


!===========================================================================================================================
! 子程序: SideHeatedcalc_wmid_max
! 作用: 计算侧壁差温工况下 z=Lz/2 中面上的 w 最大值及其位置。
!===========================================================================================================================
subroutine SideHeatedcalc_wmid_max()
  call calc_wmid_max_common('SideHeatedcalc_wmid_max')
end subroutine SideHeatedcalc_wmid_max


!===========================================================================================================================
! 子程序: SideHeatedcalc_centerline_uv_max
! 作用: 计算侧壁差温工况下 z=Lz/2 平面内两条中心线上的极大速度及其位置。
! 用途: 统计 x=Lx/2、z=Lz/2 垂直中心线上的 u 最大值及对应 y，
!       以及 y=Ly/2、z=Lz/2 水平中心线上的 v 最大值及对应 x。
!===========================================================================================================================
subroutine SideHeatedcalc_centerline_uv_max()
  use commondata3d
  implicit none
  integer(kind=4) :: i, j, iL, iR, jL, jR, kL, kR, iBest, jBest
  real(kind=8) :: targetX, targetY, targetZ, wx, wy, wz
  real(kind=8) :: val, valL, valR, valB, valT
  real(kind=8) :: umax, vmax, yAtU, xAtV

  targetX = 0.5d0 * xp(nx+1)
  targetY = 0.5d0 * yp(ny+1)
  targetZ = 0.5d0 * zp(nz+1)
  call find_bracketing_index(xp, nx, targetX, iL, iR, wx)
  call find_bracketing_index(yp, ny, targetY, jL, jR, wy)
  call find_bracketing_index(zp, nz, targetZ, kL, kR, wz)

  umax = -huge(1.0d0)
  jBest = 1
  do j = 1, ny
    ! 在 z=Lz/2 上先做 z 向插值，再在 x=Lx/2 上做 x 向插值，得到中心竖线上的 u。
    if (kL .EQ. kR) then
      valL = u(iL,j,kL)
      valR = u(iR,j,kL)
    else
      valL = (1.0d0 - wz) * u(iL,j,kL) + wz * u(iL,j,kR)
      valR = (1.0d0 - wz) * u(iR,j,kL) + wz * u(iR,j,kR)
    endif

    if (iL .EQ. iR) then
      val = valL
    else
      val = (1.0d0 - wx) * valL + wx * valR
    endif

    if (val .GT. umax) then
      umax = val
      jBest = j
    endif
  enddo
  yAtU = yp(jBest)

  vmax = -huge(1.0d0)
  iBest = 1
  do i = 1, nx
    ! 在 z=Lz/2 上先做 z 向插值，再在 y=Ly/2 上做 y 向插值，得到中心横线上的 v。
    if (kL .EQ. kR) then
      valB = v(i,jL,kL)
      valT = v(i,jR,kL)
    else
      valB = (1.0d0 - wz) * v(i,jL,kL) + wz * v(i,jL,kR)
      valT = (1.0d0 - wz) * v(i,jR,kL) + wz * v(i,jR,kR)
    endif

    if (jL .EQ. jR) then
      val = valB
    else
      val = (1.0d0 - wy) * valB + wy * valT
    endif

    if (val .GT. vmax) then
      vmax = val
      iBest = i
    endif
  enddo
  xAtV = xp(iBest)

  write(*,'(A,1X,ES24.16E3,2X,A,1X,ES24.16E3,2X,A,1X,ES24.16E3,2X,A,1X,ES24.16E3)') &
       'u_centerline_max =', umax*velocityScaleCompare, 'at y =', yAtU, 'on x =', targetX, 'z =', targetZ
  write(*,'(A,1X,ES24.16E3,2X,A,1X,ES24.16E3,2X,A,1X,ES24.16E3,2X,A,1X,ES24.16E3)') &
       'v_centerline_max =', vmax*velocityScaleCompare, 'at x =', xAtV, 'on y =', targetY, 'z =', targetZ

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) '--- SideHeatedcalc_centerline_uv_max ---'
  write(00,'(a,1x,ES24.16E3)') 'z_mid =', targetZ
  write(00,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
       'u_centerline_max =', umax*velocityScaleCompare, 'x_mid =', targetX, 'y_pos =', yAtU
  write(00,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
       'v_centerline_max =', vmax*velocityScaleCompare, 'y_mid =', targetY, 'x_pos =', xAtV
  close(00)

end subroutine SideHeatedcalc_centerline_uv_max


!===========================================================================================================================
! 子程序: SideHeatedcalc_kinetic_energy_avg
! 作用: 计算侧壁差温工况下全域平均动能 E = 0.5 * <|u|^2>。
! 用途: 使用当前速度后处理缩放后的无量纲速度，统计整个计算域的平均动能。
!===========================================================================================================================
subroutine SideHeatedcalc_kinetic_energy_avg()
  use commondata3d
  implicit none
  integer(kind=4) :: i, j, k
  real(kind=8) :: coef, energyAvg

  coef = velocityScaleCompare
  energyAvg = 0.0d0

  !$omp parallel do collapse(3) default(none) shared(u,v,w,coef) private(i,j,k) reduction(+:energyAvg)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        energyAvg = energyAvg + (coef * u(i,j,k))**2 + (coef * v(i,j,k))**2 + (coef * w(i,j,k))**2
      enddo
    enddo
  enddo
  !$omp end parallel do

  energyAvg = 0.5d0 * energyAvg / dble(nx * ny * nz)

  write(*,'(A,1X,ES24.16E3)') 'KineticEnergyAvg =', energyAvg

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) '--- SideHeatedcalc_kinetic_energy_avg ---'
  write(00,'(a,1x,ES24.16E3)') 'KineticEnergyAvg =', energyAvg
  close(00)

end subroutine SideHeatedcalc_kinetic_energy_avg


!===========================================================================================================================
! 子程序: RBcalc_umid_max
! 作用: 计算 Rayleigh-Benard 工况下 x=Lx/2 中面上的 u 最大值及其位置。
!===========================================================================================================================
subroutine RBcalc_umid_max()
  call calc_umid_max_common('RBcalc_umid_max')
end subroutine RBcalc_umid_max


!===========================================================================================================================
! 子程序: RBcalc_vmid_max
! 作用: 计算 Rayleigh-Benard 工况下 y=Ly/2 中面上的 v 最大值及其位置。
!===========================================================================================================================
subroutine RBcalc_vmid_max()
  call calc_vmid_max_common('RBcalc_vmid_max')
end subroutine RBcalc_vmid_max


!===========================================================================================================================
! 子程序: RBcalc_wmid_max
! 作用: 计算 Rayleigh-Benard 工况下 z=Lz/2 中面上的 w 最大值及其位置。
!===========================================================================================================================
subroutine RBcalc_wmid_max()
  call calc_wmid_max_common('RBcalc_wmid_max')
end subroutine RBcalc_wmid_max


!===========================================================================================================================
! 子程序: calc_umid_max_common
! 作用: 供 3D 后处理复用，统计 x=Lx/2 中面上的 u 最大值及其位置。
!===========================================================================================================================
subroutine calc_umid_max_common(logTag)
  use commondata3d
  implicit none
  character(len=*), intent(in) :: logTag
  integer(kind=4) :: j, k, iL, iR, jBest, kBest
  real(kind=8) :: targetX, weight, val, umax, yAtU, zAtU

  targetX = 0.5d0 * xp(nx+1)
  call find_bracketing_index(xp, nx, targetX, iL, iR, weight)

  umax = -huge(1.0d0)
  jBest = 1
  kBest = 1
  do k = 1, nz
    do j = 1, ny
      call interp_scalar_x(iL, iR, weight, j, k, u, val)
      if (val .GT. umax) then
        umax = val
        jBest = j
        kBest = k
      endif
    enddo
  enddo
  yAtU = yp(jBest)
  zAtU = zp(kBest)

  write(*,'(A,1X,ES24.16E3,2X,A,1X,ES24.16E3,2X,A,1X,ES24.16E3,2X,A,1X,ES24.16E3)') &
       'u_mid_max =', umax*velocityScaleCompare, 'at y =', yAtU, 'z =', zAtU, 'on x_mid =', targetX

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) '--- ', trim(logTag), ' ---'
  write(00,'(a,1x,ES24.16E3)') 'x_mid =', targetX
  write(00,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
       'u_mid_max =', umax*velocityScaleCompare, 'y_pos =', yAtU, 'z_pos =', zAtU
  close(00)

end subroutine calc_umid_max_common


!===========================================================================================================================
! 子程序: calc_vmid_max_common
! 作用: 供 3D 后处理复用，统计 y=Ly/2 中面上的 v 最大值及其位置。
!===========================================================================================================================
subroutine calc_vmid_max_common(logTag)
  use commondata3d
  implicit none
  character(len=*), intent(in) :: logTag
  integer(kind=4) :: i, k, jL, jR, iBest, kBest
  real(kind=8) :: targetY, weight, val, vmax, xAtV, zAtV

  targetY = 0.5d0 * yp(ny+1)
  call find_bracketing_index(yp, ny, targetY, jL, jR, weight)

  vmax = -huge(1.0d0)
  iBest = 1
  kBest = 1
  do k = 1, nz
    do i = 1, nx
      call interp_scalar_y(jL, jR, weight, i, k, v, val)
      if (val .GT. vmax) then
        vmax = val
        iBest = i
        kBest = k
      endif
    enddo
  enddo
  xAtV = xp(iBest)
  zAtV = zp(kBest)

  write(*,'(A,1X,ES24.16E3,2X,A,1X,ES24.16E3,2X,A,1X,ES24.16E3,2X,A,1X,ES24.16E3)') &
       'v_mid_max =', vmax*velocityScaleCompare, 'at x =', xAtV, 'z =', zAtV, 'on y_mid =', targetY

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) '--- ', trim(logTag), ' ---'
  write(00,'(a,1x,ES24.16E3)') 'y_mid =', targetY
  write(00,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
       'v_mid_max =', vmax*velocityScaleCompare, 'x_pos =', xAtV, 'z_pos =', zAtV
  close(00)

end subroutine calc_vmid_max_common


!===========================================================================================================================
! 子程序: calc_wmid_max_common
! 作用: 供 3D 后处理复用，统计 z=Lz/2 中面上的 w 最大值及其位置。
!===========================================================================================================================
subroutine calc_wmid_max_common(logTag)
  use commondata3d
  implicit none
  character(len=*), intent(in) :: logTag
  integer(kind=4) :: i, j, kL, kR, iBest, jBest
  real(kind=8) :: targetZ, weight, val, wmax, xAtW, yAtW

  ! W is the buoyancy/vertical component after the z-direction convention change.
  targetZ = 0.5d0 * zp(nz+1)
  call find_bracketing_index(zp, nz, targetZ, kL, kR, weight)

  wmax = -huge(1.0d0)
  iBest = 1
  jBest = 1
  do j = 1, ny
    do i = 1, nx
      call interp_scalar_z(kL, kR, weight, i, j, w, val)
      if (val .GT. wmax) then
        wmax = val
        iBest = i
        jBest = j
      endif
    enddo
  enddo
  xAtW = xp(iBest)
  yAtW = yp(jBest)

  write(*,'(A,1X,ES24.16E3,2X,A,1X,ES24.16E3,2X,A,1X,ES24.16E3,2X,A,1X,ES24.16E3)') &
       'w_mid_max =', wmax*velocityScaleCompare, 'at x =', xAtW, 'y =', yAtW, 'on z_mid =', targetZ

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) '--- ', trim(logTag), ' ---'
  write(00,'(a,1x,ES24.16E3)') 'z_mid =', targetZ
  write(00,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
       'w_mid_max =', wmax*velocityScaleCompare, 'x_pos =', xAtW, 'y_pos =', yAtW
  close(00)

end subroutine calc_wmid_max_common


!===========================================================================================================================
! 子程序: find_bracketing_index
! 作用: 在一维坐标数组中寻找包围目标点的左右索引及插值权重。
!===========================================================================================================================
subroutine find_bracketing_index(coord, n, target, iL, iR, weight)
  implicit none
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in) :: coord(0:n+1), target
  integer(kind=4), intent(out) :: iL, iR
  real(kind=8), intent(out) :: weight

  ! 给定目标位置 target，找到左右两个包围点以及线性插值权重
  iL = 1
  do while ((iL .LT. n) .AND. (coord(iL+1) .LT. target))
    iL = iL + 1
  enddo
  iR = min(iL + 1, n)

  if (dabs(coord(iR) - coord(iL)) .LE. 1.0d-14) then
    weight = 0.0d0
  else
    weight = (target - coord(iL)) / (coord(iR) - coord(iL))
  endif
  weight = max(0.0d0, min(1.0d0, weight))

end subroutine find_bracketing_index


!===========================================================================================================================
! 子程序: interp_scalar_x
! 作用: 在 x 方向对标量场做线性插值。
!===========================================================================================================================
subroutine interp_scalar_x(iL, iR, weight, j, k, field, val)
  use commondata3d
  implicit none
  integer(kind=4), intent(in) :: iL, iR, j, k
  real(kind=8), intent(in) :: weight
  real(kind=8), intent(in) :: field(nx,ny,nz)
  real(kind=8), intent(out) :: val

  ! 在 x 方向做一维线性插值；y、z 固定
  if (iL .EQ. iR) then
    val = field(iL,j,k)
  else
    val = (1.0d0 - weight) * field(iL,j,k) + weight * field(iR,j,k)
  endif

end subroutine interp_scalar_x


!===========================================================================================================================
! 子程序: interp_scalar_y
! 作用: 在 y 方向对标量场做线性插值。
!===========================================================================================================================
subroutine interp_scalar_y(jL, jR, weight, i, k, field, val)
  use commondata3d
  implicit none
  integer(kind=4), intent(in) :: jL, jR, i, k
  real(kind=8), intent(in) :: weight
  real(kind=8), intent(in) :: field(nx,ny,nz)
  real(kind=8), intent(out) :: val

  ! 在 y 方向做一维线性插值；x、z 固定
  if (jL .EQ. jR) then
    val = field(i,jL,k)
  else
    val = (1.0d0 - weight) * field(i,jL,k) + weight * field(i,jR,k)
  endif

end subroutine interp_scalar_y


!===========================================================================================================================
! 子程序: interp_scalar_z
! 作用: 在 z 方向对标量场做线性插值。
!===========================================================================================================================
subroutine interp_scalar_z(kL, kR, weight, i, j, field, val)
  use commondata3d
  implicit none
  integer(kind=4), intent(in) :: kL, kR, i, j
  real(kind=8), intent(in) :: weight
  real(kind=8), intent(in) :: field(nx,ny,nz)
  real(kind=8), intent(out) :: val

  ! 在 z 方向做一维线性插值；x、y 固定
  if (kL .EQ. kR) then
    val = field(i,j,kL)
  else
    val = (1.0d0 - weight) * field(i,j,kL) + weight * field(i,j,kR)
  endif

end subroutine interp_scalar_z


!===========================================================================================================================
! 子程序: write_midplane_x
! 作用: 输出 x=Lx/2 中面的 Tecplot 切片文件。
!===========================================================================================================================
subroutine write_midplane_x(filename)
  use commondata3d
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4) :: j, k, iL, iR
  real(kind=8) :: targetX, weight, valU, valV, valW, valT
  real(kind=8) :: uSlice(ny,nz), vSlice(ny,nz), wSlice(ny,nz), tSlice(ny,nz)

  ! 输出 x=Lx/2 的 y-z 截面；W_nd 是 z 方向竖直速度
  targetX = 0.5d0 * xp(nx+1)
  call find_bracketing_index(xp, nx, targetX, iL, iR, weight)
  do k = 1, nz
    do j = 1, ny
      call interp_scalar_x(iL, iR, weight, j, k, u, valU)
      call interp_scalar_x(iL, iR, weight, j, k, v, valV)
      call interp_scalar_x(iL, iR, weight, j, k, w, valW)
      call interp_scalar_x(iL, iR, weight, j, k, T, valT)
      uSlice(j,k) = velocityScaleCompare * valU
      vSlice(j,k) = velocityScaleCompare * valV
      wSlice(j,k) = velocityScaleCompare * valW
      tSlice(j,k) = valT
    enddo
  enddo

  call write_slice_fields_plt(filename, 'midX-fields', 'Y', 'Z', 'U_nd', 'V_nd', 'W_nd', 'T', &
       ny, nz, yp(1:ny), zp(1:nz), uSlice, vSlice, wSlice, tSlice)

end subroutine write_midplane_x


!===========================================================================================================================
! 子程序: write_midplane_y
! 作用: 输出 y=Ly/2 中面的 Tecplot 切片文件。
!===========================================================================================================================
subroutine write_midplane_y(filename)
  use commondata3d
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4) :: i, k, jL, jR
  real(kind=8) :: targetY, weight, valU, valV, valW, valT
  real(kind=8) :: uSlice(nx,nz), vSlice(nx,nz), wSlice(nx,nz), tSlice(nx,nz)

  ! 输出 y=Ly/2 的 x-z 截面，最适合直接看 RB 主循环胞和竖直速度 W_nd
  targetY = 0.5d0 * yp(ny+1)
  call find_bracketing_index(yp, ny, targetY, jL, jR, weight)
  do k = 1, nz
    do i = 1, nx
      call interp_scalar_y(jL, jR, weight, i, k, u, valU)
      call interp_scalar_y(jL, jR, weight, i, k, v, valV)
      call interp_scalar_y(jL, jR, weight, i, k, w, valW)
      call interp_scalar_y(jL, jR, weight, i, k, T, valT)
      uSlice(i,k) = velocityScaleCompare * valU
      vSlice(i,k) = velocityScaleCompare * valV
      wSlice(i,k) = velocityScaleCompare * valW
      tSlice(i,k) = valT
    enddo
  enddo

  call write_slice_fields_plt(filename, 'midY-fields', 'X', 'Z', 'U_nd', 'V_nd', 'W_nd', 'T', &
       nx, nz, xp(1:nx), zp(1:nz), uSlice, vSlice, wSlice, tSlice)

end subroutine write_midplane_y


!===========================================================================================================================
! 子程序: write_midplane_z
! 作用: 输出 z=Lz/2 中面的 Tecplot 切片文件。
!===========================================================================================================================
subroutine write_midplane_z(filename)
  use commondata3d
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4) :: i, j, kL, kR
  real(kind=8) :: targetZ, weight, valU, valV, valW, valT
  real(kind=8) :: uSlice(nx,ny), vSlice(nx,ny), wSlice(nx,ny), tSlice(nx,ny)

  ! 输出 z=Lz/2 的 x-y 水平中截面
  ! W is the buoyancy/vertical component after the z-direction convention change.
  targetZ = 0.5d0 * zp(nz+1)
  call find_bracketing_index(zp, nz, targetZ, kL, kR, weight)
  do j = 1, ny
    do i = 1, nx
      call interp_scalar_z(kL, kR, weight, i, j, u, valU)
      call interp_scalar_z(kL, kR, weight, i, j, v, valV)
      call interp_scalar_z(kL, kR, weight, i, j, w, valW)
      call interp_scalar_z(kL, kR, weight, i, j, T, valT)
      uSlice(i,j) = velocityScaleCompare * valU
      vSlice(i,j) = velocityScaleCompare * valV
      wSlice(i,j) = velocityScaleCompare * valW
      tSlice(i,j) = valT
    enddo
  enddo

  call write_slice_fields_plt(filename, 'midZ-fields', 'X', 'Y', 'U_nd', 'V_nd', 'W_nd', 'T', &
       nx, ny, xp(1:nx), yp(1:ny), uSlice, vSlice, wSlice, tSlice)

end subroutine write_midplane_z


!===========================================================================================================================
! 子程序: write_midplane_stream_x
! 作用: 输出 x=Lx/2 中面的流函数/涡量诊断切片。
!===========================================================================================================================
subroutine write_midplane_stream_x(filename)
  use commondata3d
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4) :: j, k, m, uout, iL, iR
  real(kind=8) :: targetX, weight
  real(kind=8) :: dy, dz, coef
  real(kind=8) :: v1, v2, vMid, inc
  real(kind=8) :: dw_dy, dv_dz
  real(kind=8) :: vSlice(ny,nz), wSlice(ny,nz)
  real(kind=8) :: psiSlice(ny,nz), vortSlice(ny,nz)

  targetX = 0.5d0 * xp(nx+1)
  call find_bracketing_index(xp, nx, targetX, iL, iR, weight)
  dy = 1.0d0 / lengthUnit
  dz = 1.0d0 / lengthUnit
  coef = velocityScaleCompare

  do k = 1, nz
    do j = 1, ny
      call interp_scalar_x(iL, iR, weight, j, k, v, vSlice(j,k))
      call interp_scalar_x(iL, iR, weight, j, k, w, wSlice(j,k))
    enddo
  enddo

  if (nz .EQ. 1) then
    psiSlice(:,1) = 0.0d0
  else
    do j = 1, ny
      v1 = vSlice(j,1) * coef
      v2 = vSlice(j,2) * coef
      psiSlice(j,1) = dz * (21.0d0 * v1 - v2) / 72.0d0
      do k = 2, nz
        m = k - 1
        if (nz .EQ. 2) then
          vMid = 0.5d0 * (vSlice(j,1) + vSlice(j,2)) * coef
        elseif (m .EQ. 1) then
          vMid = (3.0d0/8.0d0 * vSlice(j,1) + 3.0d0/4.0d0 * vSlice(j,2) - 1.0d0/8.0d0 * vSlice(j,3)) * coef
        else
          vMid = (-1.0d0/8.0d0 * vSlice(j,m-1) + 3.0d0/4.0d0 * vSlice(j,m) + 3.0d0/8.0d0 * vSlice(j,m+1)) * coef
        endif
        inc = dz / 6.0d0 * (vSlice(j,k-1) * coef + 4.0d0 * vMid + vSlice(j,k) * coef)
        psiSlice(j,k) = psiSlice(j,k-1) + inc
      enddo
    enddo
  endif

  do k = 1, nz
    do j = 1, ny
      if (ny .EQ. 1) then
        dw_dy = 0.0d0
      elseif (j .EQ. 1) then
        if (ny .EQ. 2) then
          dw_dy = (wSlice(2,k) - wSlice(1,k)) / dy
        else
          dw_dy = (-3.0d0 * wSlice(1,k) + 4.0d0 * wSlice(2,k) - wSlice(3,k)) / (2.0d0 * dy)
        endif
      elseif (j .EQ. ny) then
        if (ny .EQ. 2) then
          dw_dy = (wSlice(ny,k) - wSlice(ny-1,k)) / dy
        else
          dw_dy = (3.0d0 * wSlice(ny,k) - 4.0d0 * wSlice(ny-1,k) + wSlice(ny-2,k)) / (2.0d0 * dy)
        endif
      else
        dw_dy = (wSlice(j+1,k) - wSlice(j-1,k)) / (2.0d0 * dy)
      endif

      if (nz .EQ. 1) then
        dv_dz = 0.0d0
      elseif (k .EQ. 1) then
        if (nz .EQ. 2) then
          dv_dz = (vSlice(j,2) - vSlice(j,1)) / dz
        else
          dv_dz = (-3.0d0 * vSlice(j,1) + 4.0d0 * vSlice(j,2) - vSlice(j,3)) / (2.0d0 * dz)
        endif
      elseif (k .EQ. nz) then
        if (nz .EQ. 2) then
          dv_dz = (vSlice(j,nz) - vSlice(j,nz-1)) / dz
        else
          dv_dz = (3.0d0 * vSlice(j,nz) - 4.0d0 * vSlice(j,nz-1) + vSlice(j,nz-2)) / (2.0d0 * dz)
        endif
      else
        dv_dz = (vSlice(j,k+1) - vSlice(j,k-1)) / (2.0d0 * dz)
      endif

      vortSlice(j,k) = coef * (dw_dy - dv_dz)
    enddo
  enddo

  call write_slice_psi_vort_plt(filename, 'midX-psi-vort', 'Y', 'Z', 'PSI_YZ', 'VORT_X', &
       ny, nz, yp(1:ny), zp(1:nz), psiSlice, vortSlice)

end subroutine write_midplane_stream_x


!===========================================================================================================================
! 子程序: write_midplane_stream_y
! 作用: 输出 y=Ly/2 中面的流函数/涡量诊断切片。
!===========================================================================================================================
subroutine write_midplane_stream_y(filename)
  use commondata3d
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4) :: i, k, m, uout, jL, jR
  real(kind=8) :: targetY, weight
  real(kind=8) :: dx, dz, coef
  real(kind=8) :: u1, u2, uMid, inc
  real(kind=8) :: du_dz, dw_dx
  real(kind=8) :: uSlice(nx,nz), wSlice(nx,nz)
  real(kind=8) :: psiSlice(nx,nz), vortSlice(nx,nz)

  targetY = 0.5d0 * yp(ny+1)
  call find_bracketing_index(yp, ny, targetY, jL, jR, weight)
  dx = 1.0d0 / lengthUnit
  dz = 1.0d0 / lengthUnit
  coef = velocityScaleCompare

  do k = 1, nz
    do i = 1, nx
      call interp_scalar_y(jL, jR, weight, i, k, u, uSlice(i,k))
      call interp_scalar_y(jL, jR, weight, i, k, w, wSlice(i,k))
    enddo
  enddo

  if (nz .EQ. 1) then
    psiSlice(:,1) = 0.0d0
  else
    do i = 1, nx
      u1 = uSlice(i,1) * coef
      u2 = uSlice(i,2) * coef
      psiSlice(i,1) = dz * (21.0d0 * u1 - u2) / 72.0d0
      do k = 2, nz
        m = k - 1
        if (nz .EQ. 2) then
          uMid = 0.5d0 * (uSlice(i,1) + uSlice(i,2)) * coef
        elseif (m .EQ. 1) then
          uMid = (3.0d0/8.0d0 * uSlice(i,1) + 3.0d0/4.0d0 * uSlice(i,2) - 1.0d0/8.0d0 * uSlice(i,3)) * coef
        else
          uMid = (-1.0d0/8.0d0 * uSlice(i,m-1) + 3.0d0/4.0d0 * uSlice(i,m) + 3.0d0/8.0d0 * uSlice(i,m+1)) * coef
        endif
        inc = dz / 6.0d0 * (uSlice(i,k-1) * coef + 4.0d0 * uMid + uSlice(i,k) * coef)
        psiSlice(i,k) = psiSlice(i,k-1) + inc
      enddo
    enddo
  endif

  do k = 1, nz
    do i = 1, nx
      if (nz .EQ. 1) then
        du_dz = 0.0d0
      elseif (k .EQ. 1) then
        if (nz .EQ. 2) then
          du_dz = (uSlice(i,2) - uSlice(i,1)) / dz
        else
          du_dz = (-3.0d0 * uSlice(i,1) + 4.0d0 * uSlice(i,2) - uSlice(i,3)) / (2.0d0 * dz)
        endif
      elseif (k .EQ. nz) then
        if (nz .EQ. 2) then
          du_dz = (uSlice(i,nz) - uSlice(i,nz-1)) / dz
        else
          du_dz = (3.0d0 * uSlice(i,nz) - 4.0d0 * uSlice(i,nz-1) + uSlice(i,nz-2)) / (2.0d0 * dz)
        endif
      else
        du_dz = (uSlice(i,k+1) - uSlice(i,k-1)) / (2.0d0 * dz)
      endif

      if (nx .EQ. 1) then
        dw_dx = 0.0d0
      elseif (i .EQ. 1) then
        if (nx .EQ. 2) then
          dw_dx = (wSlice(2,k) - wSlice(1,k)) / dx
        else
          dw_dx = (-3.0d0 * wSlice(1,k) + 4.0d0 * wSlice(2,k) - wSlice(3,k)) / (2.0d0 * dx)
        endif
      elseif (i .EQ. nx) then
        if (nx .EQ. 2) then
          dw_dx = (wSlice(nx,k) - wSlice(nx-1,k)) / dx
        else
          dw_dx = (3.0d0 * wSlice(nx,k) - 4.0d0 * wSlice(nx-1,k) + wSlice(nx-2,k)) / (2.0d0 * dx)
        endif
      else
        dw_dx = (wSlice(i+1,k) - wSlice(i-1,k)) / (2.0d0 * dx)
      endif

      vortSlice(i,k) = coef * (du_dz - dw_dx)
    enddo
  enddo

  call write_slice_psi_vort_plt(filename, 'midY-psi-vort', 'X', 'Z', 'PSI_XZ', 'VORT_Y', &
       nx, nz, xp(1:nx), zp(1:nz), psiSlice, vortSlice)

end subroutine write_midplane_stream_y


!===========================================================================================================================
! 子程序: write_midplane_stream_z
! 作用: 输出 z=Lz/2 中面的流函数/涡量诊断切片。
!===========================================================================================================================
subroutine write_midplane_stream_z(filename)
  use commondata3d
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4) :: i, j, m, uout, kL, kR
  real(kind=8) :: targetZ, weight
  real(kind=8) :: dx, dy, coef
  real(kind=8) :: u1, u2, uMid, inc
  real(kind=8) :: dv_dx, du_dy
  real(kind=8) :: uSlice(nx,ny), vSlice(nx,ny)
  real(kind=8) :: psiSlice(nx,ny), vortSlice(nx,ny)

  targetZ = 0.5d0 * zp(nz+1)
  call find_bracketing_index(zp, nz, targetZ, kL, kR, weight)
  dx = 1.0d0 / lengthUnit
  dy = 1.0d0 / lengthUnit
  coef = velocityScaleCompare

  do j = 1, ny
    do i = 1, nx
      call interp_scalar_z(kL, kR, weight, i, j, u, uSlice(i,j))
      call interp_scalar_z(kL, kR, weight, i, j, v, vSlice(i,j))
    enddo
  enddo

  if (ny .EQ. 1) then
    psiSlice(:,1) = 0.0d0
  else
    do i = 1, nx
      u1 = uSlice(i,1) * coef
      u2 = uSlice(i,2) * coef
      psiSlice(i,1) = dy * (21.0d0 * u1 - u2) / 72.0d0
      do j = 2, ny
        m = j - 1
        if (ny .EQ. 2) then
          uMid = 0.5d0 * (uSlice(i,1) + uSlice(i,2)) * coef
        elseif (m .EQ. 1) then
          uMid = (3.0d0/8.0d0 * uSlice(i,1) + 3.0d0/4.0d0 * uSlice(i,2) - 1.0d0/8.0d0 * uSlice(i,3)) * coef
        else
          uMid = (-1.0d0/8.0d0 * uSlice(i,m-1) + 3.0d0/4.0d0 * uSlice(i,m) + 3.0d0/8.0d0 * uSlice(i,m+1)) * coef
        endif
        inc = dy / 6.0d0 * (uSlice(i,j-1) * coef + 4.0d0 * uMid + uSlice(i,j) * coef)
        psiSlice(i,j) = psiSlice(i,j-1) + inc
      enddo
    enddo
  endif

  do j = 1, ny
    do i = 1, nx
      if (nx .EQ. 1) then
        dv_dx = 0.0d0
      elseif (i .EQ. 1) then
        if (nx .EQ. 2) then
          dv_dx = (vSlice(2,j) - vSlice(1,j)) / dx
        else
          dv_dx = (-3.0d0 * vSlice(1,j) + 4.0d0 * vSlice(2,j) - vSlice(3,j)) / (2.0d0 * dx)
        endif
      elseif (i .EQ. nx) then
        if (nx .EQ. 2) then
          dv_dx = (vSlice(nx,j) - vSlice(nx-1,j)) / dx
        else
          dv_dx = (3.0d0 * vSlice(nx,j) - 4.0d0 * vSlice(nx-1,j) + vSlice(nx-2,j)) / (2.0d0 * dx)
        endif
      else
        dv_dx = (vSlice(i+1,j) - vSlice(i-1,j)) / (2.0d0 * dx)
      endif

      if (ny .EQ. 1) then
        du_dy = 0.0d0
      elseif (j .EQ. 1) then
        if (ny .EQ. 2) then
          du_dy = (uSlice(i,2) - uSlice(i,1)) / dy
        else
          du_dy = (-3.0d0 * uSlice(i,1) + 4.0d0 * uSlice(i,2) - uSlice(i,3)) / (2.0d0 * dy)
        endif
      elseif (j .EQ. ny) then
        if (ny .EQ. 2) then
          du_dy = (uSlice(i,ny) - uSlice(i,ny-1)) / dy
        else
          du_dy = (3.0d0 * uSlice(i,ny) - 4.0d0 * uSlice(i,ny-1) + uSlice(i,ny-2)) / (2.0d0 * dy)
        endif
      else
        du_dy = (uSlice(i,j+1) - uSlice(i,j-1)) / (2.0d0 * dy)
      endif

      vortSlice(i,j) = coef * (dv_dx - du_dy)
    enddo
  enddo

  call write_slice_psi_vort_plt(filename, 'midZ-psi-vort', 'X', 'Y', 'PSI_XY', 'VORT_Z', &
       nx, ny, xp(1:nx), yp(1:ny), psiSlice, vortSlice)

end subroutine write_midplane_stream_z


!===========================================================================================================================
! 子程序: write_full_fields_plt
! 作用: 以 Tecplot 二进制 plt 格式输出三维全场的 X/Y/Z/U/V/W/T，变量按双精度写出。
!===========================================================================================================================
subroutine write_full_fields_plt(filename)
  use commondata3d
  implicit none

  character(len=*), intent(in) :: filename
  integer(kind=4) :: i, j, k, uout
  real(kind=4) :: zoneMarker, eohMarker ! Tecplot binary markers are single-precision control records.
  character(len=40) :: zoneName

  !用“流式访问”打开,意思是不按“记录/行”来组织，而是像字节流一样连续写，适合写 .plt 这种二进制格式文件。
  open(newunit=uout, file=trim(filename), access='stream', form='unformatted', status='replace')

  zoneMarker = 299.0
  eohMarker = 357.0

  write(uout) '#!TDV101'
  write(uout) 1

  call dumpstring_to_unit(uout, 'full-field-3d')

  write(uout) 7
  call dumpstring_to_unit(uout, 'X')
  call dumpstring_to_unit(uout, 'Y')
  call dumpstring_to_unit(uout, 'Z')
  call dumpstring_to_unit(uout, 'U')
  call dumpstring_to_unit(uout, 'V')
  call dumpstring_to_unit(uout, 'W')
  call dumpstring_to_unit(uout, 'T')

  write(uout) zoneMarker
  zoneName = 'ZONE 001'
  call dumpstring_to_unit(uout, zoneName)

  write(uout) -1
  write(uout) 0
  write(uout) 1
  write(uout) 0
  write(uout) 0

  write(uout) nx
  write(uout) ny
  write(uout) nz

  write(uout) 0
  write(uout) eohMarker

  write(uout) zoneMarker

  ! 2 = Double, so Tecplot reads all seven variables in double precision.
  write(uout) 2
  write(uout) 2
  write(uout) 2
  write(uout) 2
  write(uout) 2
  write(uout) 2
  write(uout) 2

  write(uout) 0
  write(uout) -1

  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        write(uout) xp(i)
        write(uout) yp(j)
        write(uout) zp(k)
        write(uout) u(i,j,k)
        write(uout) v(i,j,k)
        write(uout) w(i,j,k)
        write(uout) T(i,j,k)
      enddo
    enddo
  enddo

  close(uout)
  return
end subroutine write_full_fields_plt


!===========================================================================================================================
! 子程序: write_slice_fields_plt
! 作用: 以 Tecplot 二进制 plt 格式输出二维切片上的 4 个场变量，变量按双精度写出。
!===========================================================================================================================
subroutine write_slice_fields_plt(filename, title, var1Name, var2Name, field1Name, field2Name, field3Name, field4Name, &
     ni, nj, coord1, coord2, field1, field2, field3, field4)
  implicit none
  character(len=*), intent(in) :: filename, title, var1Name, var2Name
  character(len=*), intent(in) :: field1Name, field2Name, field3Name, field4Name
  integer(kind=4), intent(in) :: ni, nj
  real(kind=8), intent(in) :: coord1(ni), coord2(nj)
  real(kind=8), intent(in) :: field1(ni,nj), field2(ni,nj), field3(ni,nj), field4(ni,nj)
  integer(kind=4) :: i, j, uout
  real(kind=4) :: zoneMarker, eohMarker ! Tecplot binary markers are single-precision control records.
  character(len=40) :: zoneName

  open(newunit=uout, file=trim(filename), access='stream', form='unformatted', status='replace')

  zoneMarker = 299.0
  eohMarker = 357.0

  write(uout) '#!TDV101'
  write(uout) 1

  call dumpstring_to_unit(uout, title)

  write(uout) 6
  call dumpstring_to_unit(uout, var1Name)
  call dumpstring_to_unit(uout, var2Name)
  call dumpstring_to_unit(uout, field1Name)
  call dumpstring_to_unit(uout, field2Name)
  call dumpstring_to_unit(uout, field3Name)
  call dumpstring_to_unit(uout, field4Name)

  write(uout) zoneMarker
  zoneName = 'ZONE 001'
  call dumpstring_to_unit(uout, zoneName)

  write(uout) -1
  write(uout) 0
  write(uout) 1
  write(uout) 0
  write(uout) 0

  write(uout) ni
  write(uout) nj
  write(uout) 1

  write(uout) 0
  write(uout) eohMarker

  write(uout) zoneMarker

  ! 2 = Double, so Tecplot reads all six variables in double precision.
  write(uout) 2
  write(uout) 2
  write(uout) 2
  write(uout) 2
  write(uout) 2
  write(uout) 2

  write(uout) 0
  write(uout) -1

  do j = 1, nj
    do i = 1, ni
      write(uout) coord1(i)
      write(uout) coord2(j)
      write(uout) field1(i,j)
      write(uout) field2(i,j)
      write(uout) field3(i,j)
      write(uout) field4(i,j)
    enddo
  enddo

  close(uout)
  return
end subroutine write_slice_fields_plt


!===========================================================================================================================
! 子程序: write_slice_psi_vort_plt
! 作用: 以 Tecplot 二进制 plt 格式输出二维切面的流函数/涡量数据，变量按双精度写出。
!===========================================================================================================================
subroutine write_slice_psi_vort_plt(filename, title, var1Name, var2Name, psiName, vortName, ni, nj, coord1, coord2, psi, vort)
  implicit none
  character(len=*), intent(in) :: filename, title, var1Name, var2Name, psiName, vortName
  integer(kind=4), intent(in) :: ni, nj
  real(kind=8), intent(in) :: coord1(ni), coord2(nj), psi(ni,nj), vort(ni,nj)
  integer(kind=4) :: i, j, uout
  real(kind=4) :: zoneMarker, eohMarker ! Tecplot binary markers are single-precision control records.
  character(len=40) :: zoneName

  open(newunit=uout, file=trim(filename), access='stream', form='unformatted', status='replace')

  zoneMarker = 299.0
  eohMarker = 357.0

  write(uout) '#!TDV101'
  write(uout) 1

  call dumpstring_to_unit(uout, title)

  write(uout) 4
  call dumpstring_to_unit(uout, var1Name)
  call dumpstring_to_unit(uout, var2Name)
  call dumpstring_to_unit(uout, psiName)
  call dumpstring_to_unit(uout, vortName)

  write(uout) zoneMarker
  zoneName = 'ZONE 001'
  call dumpstring_to_unit(uout, zoneName)

  write(uout) -1
  write(uout) 0
  write(uout) 1
  write(uout) 0
  write(uout) 0

  write(uout) ni
  write(uout) nj
  write(uout) 1

  write(uout) 0
  write(uout) eohMarker

  write(uout) zoneMarker

  ! 2 = Double, so Tecplot reads all four variables in double precision.
  write(uout) 2
  write(uout) 2
  write(uout) 2
  write(uout) 2

  write(uout) 0
  write(uout) -1

  do j = 1, nj
    do i = 1, ni
      write(uout) coord1(i)
      write(uout) coord2(j)
      write(uout) psi(i,j)
      write(uout) vort(i,j)
    enddo
  enddo

  close(uout)
  return
end subroutine write_slice_psi_vort_plt


!===========================================================================================================================
! 子程序: dumpstring_to_unit
! 作用: 按 Tecplot 二进制字符串格式向指定文件单元写入字符串。
!===========================================================================================================================
subroutine dumpstring_to_unit(iunit, instring)
  implicit none
  integer(kind=4), intent(in) :: iunit
  character(len=*), intent(in) :: instring
  integer(kind=4) :: stringLength, ii, ich

  stringLength = len_trim(instring)
  do ii = 1, stringLength
    ich = ichar(instring(ii:ii))
    write(iunit) ich
  enddo
  write(iunit) 0

  return
end subroutine dumpstring_to_unit
