!=============================================================
! 坐标约定：z/k/w 是浮力和竖直方向；y/j/v 是侧向方向。
!!!    注释区，代码描述
!!!    三维浮力驱动自然对流
!!!    D3Q19 流场 + D3Q7 温度场
!=============================================================

!=============================================================
!   自定义宏，一些选项的开关
#define steadyFlow
!#define unsteadyFlow

!   流动模式宏的选择
#if defined(steadyFlow) && defined(unsteadyFlow)
#error "Choose only one flow mode: steadyFlow or unsteadyFlow"
#endif
#if !defined(steadyFlow) && !defined(unsteadyFlow)
#error "Define one flow mode: steadyFlow or unsteadyFlow"
#endif

! 速度边界：上下壁为 z-normal，左右壁为 y-normal，前后壁为 x-normal。
! HorizontalWalls*: z/k 上下壁；VerticalWalls*: y/j 左右壁；SpanwiseWalls*: x/i 前后壁。
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
#error "Choose only one vertical velocity BC: VerticalWallsNoslip or VerticalWallsPeriodicalU"
#endif
#if !defined(VerticalWallsNoslip) && !defined(VerticalWallsPeriodicalU)
#error "Define one vertical velocity BC: VerticalWallsNoslip or VerticalWallsPeriodicalU"
#endif
#if defined(SpanwiseWallsNoslip) && defined(SpanwiseWallsPeriodicalU)
#error "Choose only one spanwise velocity BC: SpanwiseWallsNoslip or SpanwiseWallsPeriodicalU"
#endif
#if !defined(SpanwiseWallsNoslip) && !defined(SpanwiseWallsPeriodicalU)
#error "Define one spanwise velocity BC: SpanwiseWallsNoslip or SpanwiseWallsPeriodicalU"
#endif



! Rayleigh-Benard 温度边界：上下 z/k 壁恒温，左右 y/j 和前后 x/i 壁绝热或周期。
!#define RayleighBenardCell
!#define HorizontalWallsConstT
!#define VerticalWallsAdiabatic
!#define VerticalWallsPeriodicalT
!#define SpanwiseWallsAdiabatic
!#define SpanwiseWallsPeriodicalT




! Side Heated Cell温度边界：左右 y/j 壁恒温，上下 z/k 和前后 x/i 壁绝热。
#define SideHeatedCell
#define VerticalWallsConstT
#define HorizontalWallsAdiabatic
#define SpanwiseWallsAdiabatic
!~~温度边界条件~~

!   物理算例宏的选择
#if defined(RayleighBenardCell) && defined(SideHeatedCell)
#error "Choose only one thermal case: RayleighBenardCell or SideHeatedCell"
#endif
#if !defined(RayleighBenardCell) && !defined(SideHeatedCell)
#error "Define one thermal case: RayleighBenardCell or SideHeatedCell"
#endif
#if defined(RayleighBenardCell) && !defined(HorizontalWallsConstT)
#error "RayleighBenardCell uses z-direction hot/cold walls: define HorizontalWallsConstT"
#endif
#if defined(HorizontalWallsConstT) && (defined(HorizontalWallsAdiabatic) || defined(HorizontalWallsPeriodicalT))
#error "Choose only one horizontal temperature BC"
#endif
#if defined(VerticalWallsConstT) && (defined(VerticalWallsAdiabatic) || defined(VerticalWallsPeriodicalT))
#error "Choose only one vertical temperature BC"
#endif
#if defined(SpanwiseWallsAdiabatic) && defined(SpanwiseWallsPeriodicalT)
#error "Choose only one x-normal/front-back temperature BC: SpanwiseWallsAdiabatic or SpanwiseWallsPeriodicalT"
#endif



!算法切换
!启用 M1G 修正；注释掉则不使用 useG 相关修正
!#define EnableUseG
!启用旧算法
#define EnableLegacyThermalScheme

!   温度算法宏的选择
#if defined(EnableUseG) && defined(EnableLegacyThermalScheme)
#error "Choose only one thermal algorithm: EnableUseG or EnableLegacyThermalScheme"
#endif
#if !defined(EnableUseG) && !defined(EnableLegacyThermalScheme)
#error "Define one thermal algorithm: EnableUseG or EnableLegacyThermalScheme"
#endif


!   自定义宏结束
!=============================================================


!=============================================================
!   全局模块
module commondata3dOpenacc
  implicit none

  !===============================================================================================
  ! 格子离散速度数
  integer(kind=4), parameter :: qf=19, qt=7
  !===============================================================================================

  !===============================================================================================
  ! 是否在计算前从旧算例重启
  integer(kind=4), parameter :: loadInitField=0

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
  integer(kind=4), parameter :: nx=120, ny=120, nz=120
#ifdef SideHeatedCell
  real(kind=8), parameter :: lengthUnit=dble(ny)     ! 侧壁差温：冷热壁距离沿 y 方向。
#else
  real(kind=8), parameter :: lengthUnit=dble(nz)     ! RB：冷热壁距离沿 z 方向。
#endif
  real(kind=8), parameter :: pi=acos(-1.0d0)

  real(kind=8), parameter :: Rayleigh=1.0d7
  real(kind=8), parameter :: Prandtl=0.71d0
  real(kind=8), parameter :: Mach=0.1d0
  real(kind=8), parameter :: Thot=0.5d0, Tcold=-0.5d0
  real(kind=8), parameter :: Tref=0.5d0*(Thot+Tcold)
  real(kind=8), parameter :: tauf=0.5d0+Mach*lengthUnit*dsqrt(3.0d0*Prandtl/Rayleigh)
  real(kind=8), parameter :: viscosity=(tauf-0.5d0)/3.0d0
  real(kind=8), parameter :: diffusivity=viscosity/Prandtl

  real(kind=8), parameter :: cs2T=0.25d0

  ! 高阶矩参数修正 aT
  real(kind=8), parameter :: paraA=42.0d0*dsqrt(3.0d0)*diffusivity-6.0d0

  ! velocityScaleCompare 把格子速度换算到热扩散标度 U*L/alpha。
  ! 它用于速度输出，也用于 Nu 诊断中的对流热通量部分。
  real(kind=8), parameter :: velocityScaleCompare=lengthUnit/diffusivity

  ! 浮力项参数
  ! 浮力沿 z 方向施加；输出文件中的 W_nd 是竖直速度。
  real(kind=8), parameter :: gBeta1=Rayleigh*viscosity*diffusivity/lengthUnit
  real(kind=8), parameter :: gBeta=gBeta1/lengthUnit/lengthUnit
  real(kind=8), parameter :: timeUnit=dsqrt(lengthUnit/gBeta)
  real(kind=8), parameter :: velocityUnit=dsqrt(gBeta*lengthUnit)

  ! 动量方程的多松弛系数
  real(kind=8), parameter :: Se=1.0d0/tauf, Seps=1.0d0/tauf
  real(kind=8), parameter :: Snu=1.0d0/tauf, Spi=1.0d0/tauf
  real(kind=8), parameter :: Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)
  real(kind=8), parameter :: Sm=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)

  ! 温度方程的多松弛系数
#ifdef EnableLegacyThermalScheme
  real(kind=8), parameter :: Qk=3.0d0-dsqrt(3.0d0), Qnu=4.0d0*dsqrt(3.0d0)-6.0d0
  real(kind=8), parameter :: thermalGeqCoeff=7.0d0/((6.0d0+paraA)*cs2T)
#else
  real(kind=8), parameter :: taug=0.5d0+diffusivity/cs2T
  real(kind=8), parameter :: Qnu=1.0d0, Qk=1.0d0/taug
  real(kind=8), parameter :: thermalGeqCoeff=1.0d0/cs2T
#endif

  !===============================================================================================
  ! 输出/备份相关设置（以自由落体时间 t_ff 为单位）
  real(kind=8), parameter :: epsU=1.0d-7, epsT=1.0d-7

#ifdef steadyFlow
  real(kind=8), parameter :: outputSnapshotInterval=10.0d0   ! Nu/Re 时间序列采样间隔（单位：t_ff）
  real(kind=8), parameter :: reloadFileInterval=100.0d0  ! f/g 重启文件输出间隔（单位：t_ff）
  real(kind=8), parameter :: outputPltFileInterval=100.0d0  ! Tecplot 文件输出间隔（单位：t_ff）
  integer(kind=4), parameter :: dimensionlessTimeMax=int(12000.0d0/outputSnapshotInterval)
  integer(kind=4), parameter :: outputSnapshotFile=1
  integer(kind=4), parameter :: outputPltFile=1
  integer(kind=4), parameter :: outputReloadFile=1
  integer(kind=4), parameter :: itc_max=20000000
#endif

#ifdef unsteadyFlow
  real(kind=8), parameter :: outputSnapshotInterval=0.5d0   ! uvwTrho 和 Nu/Re 时间序列采样间隔（单位：t_ff）
  real(kind=8), parameter :: reloadFileInterval=100.0d0  ! f/g 重启文件输出间隔（单位：t_ff）
  real(kind=8), parameter :: outputPltFileInterval=100.0d0  ! Tecplot 文件周期输出间隔（单位：t_ff）
  real(kind=8), parameter :: unsteadyRunDuration=1000.0d0
  ! 以下三个参数只控制非稳态结束后的 Nu/Re 统计平均窗口，不改变推进时长或采样频率。
  ! 时间窗口按完整算例的绝对 t_ff 计；续算后处理会从完整 .dat 历史重建。
  real(kind=8), parameter :: unsteadyAverageStartTf=0.5d0*unsteadyRunDuration  ! 平均窗口起点
  real(kind=8), parameter :: unsteadyAverageEndTf=unsteadyRunDuration          ! 平均窗口终点
  real(kind=8), parameter :: unsteadyAverageMidTf=0.5d0*(unsteadyAverageStartTf+unsteadyAverageEndTf) ! 前/后半分界
  integer(kind=4), parameter :: unsteadySampleCount=max(1, int(unsteadyRunDuration/outputSnapshotInterval+0.5d0))
  integer(kind=4), parameter :: dimensionlessTimeMax=unsteadySampleCount
  integer(kind=4), parameter :: outputSnapshotFile=1
  integer(kind=4), parameter :: outputPltFile=1
  integer(kind=4), parameter :: outputReloadFile=1
  integer(kind=4), parameter :: itc_max=max(1, int(unsteadyRunDuration*timeUnit+0.5d0))
#endif

  integer(kind=4) :: snapshotFileNum, pltFileNum
  integer(kind=4) :: dimensionlessTime
  integer(kind=4) :: outputSnapshotIntervalItc
  integer(kind=4) :: reloadFileIntervalItc, outputPltFileIntervalItc

  real(kind=8) :: NuVolAvg(0:dimensionlessTimeMax), ReVolAvg(0:dimensionlessTimeMax)

  character(len=100) :: snapshotFilePrefix="buoyancyCavity3DOpenaccSnapshot"
  character(len=100) :: pltFolderPrefix="buoyancyCavity3DOpenaccTecplot"
  character(len=100) :: reloadFilePrefix="reloadFile3DOpenacc"
  character(len=100) :: settingsFile="SimulationSettings3DOpenacc.txt"
  !===============================================================================================

  real(kind=8) :: errorU, errorT

  real(kind=8) :: xp(0:nx+1), yp(0:ny+1), zp(0:nz+1)
  real(kind=8), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:), T(:,:,:), rho(:,:,:)

#ifdef steadyFlow
  ! 稳态误差判据需要保存上一次输出时刻的场
  real(kind=8), allocatable :: up(:,:,:), vp(:,:,:), wp(:,:,:), Tp(:,:,:)
#endif

  real(kind=8), allocatable :: f(:,:,:,:), f_post(:,:,:,:)
  real(kind=8), allocatable :: g(:,:,:,:), g_post(:,:,:,:)
  real(kind=8), allocatable :: Bx_prev(:,:,:), By_prev(:,:,:), Bz_prev(:,:,:)

  integer(kind=4) :: itc

#ifdef EnableUseG
  logical, parameter :: useG = .true.
#else
  logical, parameter :: useG = .false.
#endif

#ifdef EnableLegacyThermalScheme
  logical, parameter :: useLegacyThermalScheme=.true.
#else
  logical, parameter :: useLegacyThermalScheme=.false.
#endif

#ifdef steadyFlow
  real(kind=8) :: Nu_global, Nu_hot, Nu_cold, Nu_middle
  real(kind=8) :: Nu_hot_max, Nu_hot_min, Nu_hot_max_position, Nu_hot_min_position
#endif

  integer(kind=4) :: ex(0:qf-1), ey(0:qf-1), ez(0:qf-1), opp(0:qf-1)
  real(kind=8)    :: omega(0:qf-1)

  integer(kind=4) :: exT(0:qt-1), eyT(0:qt-1), ezT(0:qt-1), oppT(0:qt-1)
  real(kind=8)    :: omegaT(0:qt-1)
  !===============================================================================================
end module commondata3dOpenacc


!   全局模块结束
!=============================================================


!=============================================================
!   主程序


program main3dOpenacc
  use openacc
  use commondata3dOpenacc
  implicit none

  real(kind=8) :: timeStart, timeEnd
  real(kind=8) :: timeStart2, timeEnd2
  character(len=24) :: ctime
  character(len=24) :: string
  integer(kind=4) :: time
  integer(kind=4) :: numAccDevices
  integer(kind=8) :: wallClockStart, wallClockEnd, wallClockRate
#ifdef unsteadyFlow
  integer(kind=4) :: nextSampleItc, nextSampleAbsItc, unsteadyItcRemaining
#endif

  ! 先写日志头，并初始化 OpenACC 设备
  if (loadInitField .EQ. 1) then
    open(unit=00, file=trim(settingsFile), status='unknown', position='append')
    write(00,*) '================ Restart continuation begins ================'
  else
    open(unit=00, file=trim(settingsFile), status='replace')
  endif
  string = ctime(time())
  write(00,*) 'Start: ', string
  write(00,*) 'Starting OpenACC >>>>>>'
  call acc_init(acc_device_default)
  numAccDevices = acc_get_num_devices(acc_device_default)
  write(00,*) 'Visible OpenACC devices:', numAccDevices
  close(00)

  call initial3d()
#ifdef unsteadyFlow
  ! 非稳态的 itc_max 是整个算例的总目标步数；
  ! 续算时 restartItcOffset 是旧算例已经完成的步数，本次只推进剩余步数。
  unsteadyItcRemaining = max(0, itc_max - restartItcOffset)
#endif
  call enter_data_3d_openacc()     !把主要数组和常量映射到 OpenACC 设备端

  call CPU_TIME(timeStart)
  call system_clock(wallClockStart, wallClockRate)   !wallClockStart 保存“当前这一刻”的时钟计数值, tick 数
  timeStart2 = dble(wallClockStart) / dble(max(wallClockRate,1_8))  !把“tick 数”换成“秒数”

#ifdef steadyFlow
  do while (((errorU .GT. epsU) .OR. (errorT .GT. epsT)) .AND. (itc .LE. itc_max))
#endif
#ifdef unsteadyFlow
  do while (itc .LT. unsteadyItcRemaining)
#endif
    itc = itc + 1

    call collision3d()
    call streaming3d()
    call bounceback3d()
    call macro3d()

    call collisionT3d()
    call streamingT3d()
    call bouncebackT3d()
    call macroT3d()

#ifdef steadyFlow
    ! 周期输出按累计格子步判断，续算时才能接回不断电运行应有的输出节奏。
    if (mod(restartItcOffset+itc, 2000) .EQ. 0) call check3d()
    if ((outputPltFile .EQ. 1) .AND. (mod(restartItcOffset+itc, outputPltFileIntervalItc) .EQ. 0)) then
      call update_host_tecplot_3d_openacc()
      call output_Tecplot3d()
    endif
    if ((outputSnapshotFile .EQ. 1) .AND. (mod(restartItcOffset+itc, outputSnapshotIntervalItc) .EQ. 0)) then
      call update_host_snapshot_3d_openacc()
      call output_SnapshotFile3d()
    endif
    if ((outputReloadFile .EQ. 1) .AND. (mod(restartItcOffset+itc, reloadFileIntervalItc) .EQ. 0)) then
      call update_host_reload_3d_openacc()
      call output_ReloadFile3d()
    endif
#endif

#ifdef unsteadyFlow
    do while( (reloadDimensionlessTime + real(dimensionlessTime,kind=8)*outputSnapshotInterval) &
         .LT. unsteadyRunDuration )
      ! 每个目标采样时刻都按绝对 t_ff 换算到本次运行段的 itc，续算时不会重复旧样本。
      nextSampleAbsItc = max(1, int((reloadDimensionlessTime + &
           real(dimensionlessTime+1,kind=8)*outputSnapshotInterval)*timeUnit+0.5d0))
      nextSampleItc = max(1, nextSampleAbsItc - restartItcOffset)
      if(itc .LT. nextSampleItc) exit
      call calNuRe3d()
      if (outputSnapshotFile .EQ. 1) then
        call update_host_snapshot_3d_openacc()
        call output_SnapshotFile3d()     !每 0.5 t_ff 输出一次 u、v、w、T、rho 的二进制快照文件
      endif
    enddo
    if ((outputPltFile .EQ. 1) .AND. (mod(restartItcOffset+itc, outputPltFileIntervalItc) .EQ. 0)) then
      call update_host_tecplot_3d_openacc()
      call output_Tecplot3d()
    endif
    if ((outputReloadFile .EQ. 1) .AND. (mod(restartItcOffset+itc, reloadFileIntervalItc) .EQ. 0)) then
      call update_host_reload_3d_openacc()
      call output_ReloadFile3d()
    endif
#endif
  enddo

  !$acc wait(1)
  call CPU_TIME(timeEnd)
  call system_clock(wallClockEnd, wallClockRate)
  timeEnd2 = dble(wallClockEnd) / dble(max(wallClockRate,1_8))

#ifdef steadyFlow
  call update_host_snapshot_3d_openacc()
  call output_Tecplot3d()
  call output_SnapshotFile3d()
#endif

#ifdef unsteadyFlow
  call output_unsteady_NuRe_postprocess3d()
#endif


#ifdef steadyFlow
! 这些最终标量诊断只用于稳态收敛后评估
#ifdef SideHeatedCell
  call SideHeatedcalc_Nu_global3d()
  call SideHeatedcalc_Nu_wall_avg3d()
  call SideHeatedcalc_umid_max3d()
  call SideHeatedcalc_vmid_max3d()
  call SideHeatedcalc_wmid_max3d()
#endif

#ifdef RayleighBenardCell
  call RBcalc_Nu_global3d()
  call RBcalc_Nu_wall_avg3d()
  call RBcalc_umid_max3d()
  call RBcalc_vmid_max3d()
  call RBcalc_wmid_max3d()
#endif

  call calNuRe3d()
#endif


  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) '======================================================================'
  write(00,*) 'Time (CPU) = ', real(timeEnd - timeStart, kind=8), 's'
  write(00,*) 'MLUPS = ', &
       real(dble(nx) * dble(ny) * dble(nz) * dble(itc) / &
       & max(timeEnd - timeStart, 1.0d-12) / 1.0d6, kind=8)
  write(00,*) 'Time (ACC) = ', real(timeEnd2 - timeStart2, kind=8), 's'
  write(00,*) 'MLUPS (ACC) = ', &
       real(dble(nx) * dble(ny) * dble(nz) * dble(itc) / &
       & max(timeEnd2 - timeStart2, 1.0d-12) / 1.0d6, kind=8)
#ifdef steadyFlow
  write(00,'(a,1x,ES24.16E3)') 'Nu_global =', real(Nu_global,kind=8)
  write(00,'(a,1x,ES24.16E3)') 'Nu_hot    =', real(Nu_hot,kind=8)
  write(00,'(a,1x,ES24.16E3)') 'Nu_cold   =', real(Nu_cold,kind=8)
  write(00,'(a,1x,ES24.16E3)') 'Nu_middle =', real(Nu_middle,kind=8)
#endif
  write(00,*) 'useG =', useG
  write(00,*) 'useLegacyThermalScheme =', useLegacyThermalScheme
  write(00,*) 'Deallocate Array......'
  close(00)

  call exit_data_3d_openacc()

  deallocate(f, f_post, g, g_post)
  deallocate(u, v, w, T, rho)
#ifdef steadyFlow
  deallocate(up, vp, wp, Tp)
#endif
  deallocate(Bx_prev, By_prev, Bz_prev)

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) 'Successfully: DNS completed!'
  string = ctime(time())
  write(00,*) 'End:   ', string
  close(00)

end program main3dOpenacc

!   主程序结束
!=============================================================



!===========================================================================================================================
! 子程序: initial3d
! 作用: 初始化网格坐标、场变量、分布函数、输出文件和重启信息。
!===========================================================================================================================
subroutine initial3d()
  use commondata3dOpenacc
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

  ! 把按自由落体时间给出的输出/备份间隔换算成格子步数 itc
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

  if(outputReloadFile.EQ.1) then
        open(unit=01,file=trim(reloadFilePrefix)//"-"//"readme",status="unknown")
        write(01,*) "reloadFile prefix exists!"
        close(01)
        write(00,*) "Reload data will be stored in ", reloadFilePrefix
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
  write(00,*) 'outputSnapshotFile =', outputSnapshotFile
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
  write(00,*) 'OpenACC GPU version; MPI/OpenMP are not included in this file'

  !-----------------------------------------------------------------------------------------------



  !-----------------------------------------------------------------------------------------------
  ! 节点坐标数组
  !-----------------------------------------------------------------------------------------------
  xp(0) = 0.0d0
  xp(nx+1) = dble(nx)
  do i = 1, nx
    xp(i) = dble(i) - 0.5d0
  enddo
  xp = xp / lengthUnit

  yp(0) = 0.0d0
  yp(ny+1) = dble(ny)
  do j = 1, ny
    yp(j) = dble(j) - 0.5d0
  enddo
  yp = yp / lengthUnit

  zp(0) = 0.0d0
  zp(nz+1) = dble(nz)
  do k = 1, nz
    zp(k) = dble(k) - 0.5d0
  enddo
  zp = zp / lengthUnit



  allocate(u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), T(nx,ny,nz), rho(nx,ny,nz))
#ifdef steadyFlow
  allocate(up(nx,ny,nz), vp(nx,ny,nz), wp(nx,ny,nz), Tp(nx,ny,nz))
#endif
  allocate(f(nx,ny,nz,0:qf-1), f_post(0:nx+1,0:ny+1,0:nz+1,0:qf-1))
  allocate(g(nx,ny,nz,0:qt-1), g_post(0:nx+1,0:ny+1,0:nz+1,0:qt-1))
  allocate(Bx_prev(nx,ny,nz), By_prev(nx,ny,nz), Bz_prev(nx,ny,nz))

  !-----------------------------------------------------------------------------------------------



  !-----------------------------------------------------------------------------------------------
  ! 初始化
  !-----------------------------------------------------------------------------------------------
  call init_lattice_constants_3d()  !速度集和权重

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
    if (Rayleigh .LT. 1.0d4) then
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
      write(00,'(a,1x,ES24.16E3)') '3D RB initial T perturbation amplitude =', real(rbInitPerturbAmp,kind=8)
    else
      write(00,*) '3D RB initial T perturbation skipped because Rayleigh >= 1.0d4'
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
            f(i,j,k,alpha) = omega(alpha) * rho(i,j,k) * (1.0d0 + 3.0d0 * eu + 4.5d0 * eu * eu - 1.5d0 * u2Loc)
          enddo
          do alpha = 0, qt-1
            eu = dble(exT(alpha)) * u(i,j,k) + dble(eyT(alpha)) * v(i,j,k) + dble(ezT(alpha)) * w(i,j,k)
            g(i,j,k,alpha) = omegaT(alpha) * T(i,j,k) * (1.0d0 + thermalGeqCoeff * eu)
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
    ! 这一步在 enter_data_3d_openacc() 之前完成，后续 copyin 会把恢复后的主机端状态送到 GPU。
    write(00,*) 'Read reload metadata before choosing the restart .bin file.'
    write(reloadFileName,'(i12.12)') reloadFileNum             ! latest .meta 缺失时才依赖这个手工编号
    reloadFileName = adjustl(reloadFileName)                  ! adjustl 把前导空格移到字符串末尾
    call read_reload_metadata3d(reloadFileName)
    write(00,*) 'Load initial field from previous simulation: ', &
         trim(reloadFilePrefix), '-', trim(reloadFileName), '.bin'
    if (.not. reloadMetadataLoaded) then
      write(00,*) 'WARNING: no reload metadata file found; restart offsets were inferred.'
      write(00,*) '         For exact continuation, use reload files written after this patch.'
    endif
    open(unit=01, file=trim(reloadFilePrefix)//'-'//trim(reloadFileName)//'.bin', &
         form='unformatted', access='sequential', status='old')
    write(00,*) 'Reloading strict restart state from file'
    read(01) ((((f(i,j,k,alpha), i=1,nx), j=1,ny), k=1,nz), alpha=0,qf-1)
    read(01) ((((g(i,j,k,alpha), i=1,nx), j=1,ny), k=1,nz), alpha=0,qt-1)
#ifdef EnableUseG
    read(01) (((Bx_prev(i,j,k), i=1,nx), j=1,ny), k=1,nz)
    read(01) (((By_prev(i,j,k), i=1,nx), j=1,ny), k=1,nz)
    read(01) (((Bz_prev(i,j,k), i=1,nx), j=1,ny), k=1,nz)
#endif
    close(01)
    call reconstruct_macro_from_fg3d()
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

end subroutine initial3d
!===========================================================================================================================
! 子程序: enter_data_3d_openacc
! 作用: 在主时间推进前把主要数组和常量映射到 OpenACC 设备端。
!===========================================================================================================================
subroutine enter_data_3d_openacc()
  use openacc
  use commondata3dOpenacc
  implicit none

  ! copyin: 把主机端已有初值复制到设备端，并让这些数组在整个时间推进期间常驻 GPU。
  !$acc enter data copyin(xp,yp,zp,ex,ey,ez,opp,omega,exT,eyT,ezT,oppT,omegaT)
  !$acc enter data copyin(u,v,w,T,rho,f,g,Bx_prev,By_prev,Bz_prev)
  ! create: 只在设备端分配中间缓冲，不从主机复制初值；f_post/g_post 会在 GPU kernel 内被完全覆盖。
  !$acc enter data create(f_post,g_post)
#ifdef steadyFlow
  ! 稳态误差判据需要保留上一时刻场，因此也放到设备端常驻。
  !$acc enter data copyin(up,vp,wp,Tp)
#endif
end subroutine enter_data_3d_openacc


!===========================================================================================================================
! 子程序: update_host_snapshot_3d_openacc
! 作用: 输出 u/v/w/T/rho 快照前，把对应宏观场同步回主机端。
!===========================================================================================================================
subroutine update_host_snapshot_3d_openacc()
  use commondata3dOpenacc
  implicit none

  !$acc wait(1)
  !$acc update self(u,v,w,T,rho)
end subroutine update_host_snapshot_3d_openacc


!===========================================================================================================================
! 子程序: update_host_tecplot_3d_openacc
! 作用: 输出 Tecplot 主场和 CPU 后处理切片前，同步可视化需要的宏观场。
!===========================================================================================================================
subroutine update_host_tecplot_3d_openacc()
  use commondata3dOpenacc
  implicit none

  !$acc wait(1)
  !$acc update self(u,v,w,T)
end subroutine update_host_tecplot_3d_openacc


!===========================================================================================================================
! 子程序: update_host_reload_3d_openacc
! 作用: 写严格重启文件前，只同步分布函数和 EnableUseG 的历史热流。
!===========================================================================================================================
subroutine update_host_reload_3d_openacc()
  use commondata3dOpenacc
  implicit none

  !$acc wait(1)
  !$acc update self(f,g)
#ifdef EnableUseG
  !$acc update self(Bx_prev,By_prev,Bz_prev)
#endif
end subroutine update_host_reload_3d_openacc


!===========================================================================================================================
! 子程序: exit_data_3d_openacc
! 作用: 在计算结束后释放设备端数据映射。
!===========================================================================================================================
subroutine exit_data_3d_openacc()
  use commondata3dOpenacc
  implicit none

#ifdef steadyFlow
  ! delete: 释放设备端映射关系；不会删除主机端的 Fortran 数组，只是清掉 GPU 常驻数据。
  !$acc exit data delete(up,vp,wp,Tp)
#endif
  !$acc exit data delete(f_post,g_post,u,v,w,T,rho,f,g,Bx_prev,By_prev,Bz_prev)
  !$acc exit data delete(xp,yp,zp,ex,ey,ez,opp,omega,exT,eyT,ezT,oppT,omegaT)
end subroutine exit_data_3d_openacc


!===========================================================================================================================
! 子程序: init_lattice_constants_3d
! 作用: 初始化 D3Q19 / D3Q7 的离散速度、反向索引和权重。
!===========================================================================================================================
subroutine init_lattice_constants_3d()
  use commondata3dOpenacc
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

end subroutine init_lattice_constants_3d



!===========================================================================================================================
! 子程序: collision3d
! 作用: 流场碰撞步骤，在矩空间完成松弛并加入浮力源项修正。
!===========================================================================================================================
subroutine collision3d()
  use commondata3dOpenacc
  implicit none

  integer(kind=4) :: i, j, k, alpha
  real(kind=8) :: m(0:qf-1), meq(0:qf-1), m_post(0:qf-1)
  real(kind=8) :: s(0:qf-1), fSource(0:qf-1)
  real(kind=8) :: rhoLoc, uLoc, vLoc, wLoc, u2, uDotF
  real(kind=8) :: FxLoc, FyLoc, FzLoc

  ! 流场采用 D3Q19 MRT。这里把正变换和逆变换都显式展开
  ! parallel loop : 为整个三重循环生成一个 GPU kernel。
  ! gang/vector   : 指示 OpenACC 把迭代映射到粗粒度/细粒度并行层次。
  ! collapse(3)   : 把 i-j-k 三层循环展平，扩大可并行的迭代空间。
  ! default(none) : 强制显式写出数据属性，避免变量被编译器偷偷按默认规则处理。
  ! present(...)  : 这些数组已经由 enter data 放到 GPU，这里直接使用，不再重复拷贝。
  ! async(1)      : 投递到编号 1 的异步队列；同一队列中的 kernel 按顺序执行，但 CPU 不必原地等待。
  ! private(...)  : 每个并行迭代拥有自己的临时标量/小数组，避免线程之间相互覆盖。
  !$acc parallel loop gang vector collapse(3) default(none) present(f,f_post,rho,u,v,w,T) async(1) &
  !$acc& private(alpha,m,meq,m_post,s,fSource,rhoLoc,uLoc,vLoc,wLoc,u2,uDotF,FxLoc,FyLoc,FzLoc)
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
        m(0) = f(i,j,k,0) + f(i,j,k,1) + f(i,j,k,2) + f(i,j,k,3) + f(i,j,k,4) + f(i,j,k,5) + &
             f(i,j,k,6) + f(i,j,k,7) + f(i,j,k,8) + f(i,j,k,9) + f(i,j,k,10) + f(i,j,k,11) + &
             f(i,j,k,12) + f(i,j,k,13) + f(i,j,k,14) + f(i,j,k,15) + f(i,j,k,16) + f(i,j,k,17) + f(i,j,k,18)

        m(1) = -30.0d0 * f(i,j,k,0) - 11.0d0 * (f(i,j,k,1) + f(i,j,k,2) + f(i,j,k,3) + f(i,j,k,4) + &
             f(i,j,k,5) + f(i,j,k,6)) + 8.0d0 * (f(i,j,k,7) + f(i,j,k,8) + f(i,j,k,9) + f(i,j,k,10) + &
             f(i,j,k,11) + f(i,j,k,12) + f(i,j,k,13) + f(i,j,k,14) + f(i,j,k,15) + f(i,j,k,16) + &
             f(i,j,k,17) + f(i,j,k,18))

        m(2) = 12.0d0 * f(i,j,k,0) - 4.0d0 * (f(i,j,k,1) + f(i,j,k,2) + f(i,j,k,3) + f(i,j,k,4) + &
             f(i,j,k,5) + f(i,j,k,6)) + (f(i,j,k,7) + f(i,j,k,8) + f(i,j,k,9) + f(i,j,k,10) + &
             f(i,j,k,11) + f(i,j,k,12) + f(i,j,k,13) + f(i,j,k,14) + f(i,j,k,15) + f(i,j,k,16) + &
             f(i,j,k,17) + f(i,j,k,18))

        m(3) =  f(i,j,k,1) - f(i,j,k,2) + f(i,j,k,7) - f(i,j,k,8) + f(i,j,k,9) - f(i,j,k,10) + &
             f(i,j,k,11) - f(i,j,k,12) + f(i,j,k,13) - f(i,j,k,14)

        m(4) = -4.0d0 * f(i,j,k,1) + 4.0d0 * f(i,j,k,2) + f(i,j,k,7) - f(i,j,k,8) + f(i,j,k,9) - &
             f(i,j,k,10) + f(i,j,k,11) - f(i,j,k,12) + f(i,j,k,13) - f(i,j,k,14)

        m(5) =  f(i,j,k,3) - f(i,j,k,4) + f(i,j,k,7) + f(i,j,k,8) - f(i,j,k,9) - f(i,j,k,10) + &
             f(i,j,k,15) - f(i,j,k,16) + f(i,j,k,17) - f(i,j,k,18)

        m(6) = -4.0d0 * f(i,j,k,3) + 4.0d0 * f(i,j,k,4) + f(i,j,k,7) + f(i,j,k,8) - f(i,j,k,9) - &
             f(i,j,k,10) + f(i,j,k,15) - f(i,j,k,16) + f(i,j,k,17) - f(i,j,k,18)

        m(7) =  f(i,j,k,5) - f(i,j,k,6) + f(i,j,k,11) + f(i,j,k,12) - f(i,j,k,13) - f(i,j,k,14) + &
             f(i,j,k,15) + f(i,j,k,16) - f(i,j,k,17) - f(i,j,k,18)

        m(8) = -4.0d0 * f(i,j,k,5) + 4.0d0 * f(i,j,k,6) + f(i,j,k,11) + f(i,j,k,12) - f(i,j,k,13) - &
             f(i,j,k,14) + f(i,j,k,15) + f(i,j,k,16) - f(i,j,k,17) - f(i,j,k,18)

        m(9) =  2.0d0 * (f(i,j,k,1) + f(i,j,k,2)) - (f(i,j,k,3) + f(i,j,k,4) + f(i,j,k,5) + f(i,j,k,6)) + &
             (f(i,j,k,7) + f(i,j,k,8) + f(i,j,k,9) + f(i,j,k,10) + f(i,j,k,11) + f(i,j,k,12) + &
             f(i,j,k,13) + f(i,j,k,14)) - 2.0d0 * (f(i,j,k,15) + f(i,j,k,16) + f(i,j,k,17) + f(i,j,k,18))

        m(10) = -4.0d0 * (f(i,j,k,1) + f(i,j,k,2)) + 2.0d0 * (f(i,j,k,3) + f(i,j,k,4) + f(i,j,k,5) + f(i,j,k,6)) + &
             (f(i,j,k,7) + f(i,j,k,8) + f(i,j,k,9) + f(i,j,k,10) + f(i,j,k,11) + f(i,j,k,12) + &
             f(i,j,k,13) + f(i,j,k,14)) - 2.0d0 * (f(i,j,k,15) + f(i,j,k,16) + f(i,j,k,17) + f(i,j,k,18))

        m(11) =  (f(i,j,k,3) + f(i,j,k,4)) - (f(i,j,k,5) + f(i,j,k,6)) + &
             (f(i,j,k,7) + f(i,j,k,8) + f(i,j,k,9) + f(i,j,k,10)) - &
             (f(i,j,k,11) + f(i,j,k,12) + f(i,j,k,13) + f(i,j,k,14))

        m(12) = -2.0d0 * (f(i,j,k,3) + f(i,j,k,4)) + 2.0d0 * (f(i,j,k,5) + f(i,j,k,6)) + &
             (f(i,j,k,7) + f(i,j,k,8) + f(i,j,k,9) + f(i,j,k,10)) - &
             (f(i,j,k,11) + f(i,j,k,12) + f(i,j,k,13) + f(i,j,k,14))

        m(13) =  f(i,j,k,7) - f(i,j,k,8) - f(i,j,k,9) + f(i,j,k,10)

        m(14) =  f(i,j,k,15) - f(i,j,k,16) - f(i,j,k,17) + f(i,j,k,18)

        m(15) =  f(i,j,k,11) - f(i,j,k,12) - f(i,j,k,13) + f(i,j,k,14)

        m(16) =  f(i,j,k,7) - f(i,j,k,8) + f(i,j,k,9) - f(i,j,k,10) - f(i,j,k,11) + f(i,j,k,12) - &
             f(i,j,k,13) + f(i,j,k,14)

        m(17) = -f(i,j,k,7) - f(i,j,k,8) + f(i,j,k,9) + f(i,j,k,10) + f(i,j,k,15) - f(i,j,k,16) + &
             f(i,j,k,17) - f(i,j,k,18)

        m(18) =  f(i,j,k,11) + f(i,j,k,12) - f(i,j,k,13) - f(i,j,k,14) - f(i,j,k,15) - f(i,j,k,16) + &
             f(i,j,k,17) + f(i,j,k,18)

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
        s(0)  = 0.0d0
        s(1)  = Se
        s(2)  = Seps
        s(3)  = 0.0d0
        s(4)  = Sq
        s(5)  = 0.0d0
        s(6)  = Sq
        s(7)  = 0.0d0
        s(8)  = Sq
        s(9)  = Snu
        s(10) = Spi
        s(11) = Snu
        s(12) = Spi
        s(13) = Snu
        s(14) = Snu
        s(15) = Snu
        s(16) = Sm
        s(17) = Sm
        s(18) = Sm

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
        f_post(i,j,k,0) =  m_post(0) / 19.0d0 - 5.0d0 * m_post(1) / 399.0d0 + m_post(2) / 21.0d0

        f_post(i,j,k,1) =  m_post(0) / 19.0d0 - 11.0d0 * m_post(1) / 2394.0d0 - m_post(2) / 63.0d0 + &
             m_post(3) / 10.0d0 - m_post(4) / 10.0d0 + m_post(9) / 18.0d0 - m_post(10) / 18.0d0

        f_post(i,j,k,2) =  m_post(0) / 19.0d0 - 11.0d0 * m_post(1) / 2394.0d0 - m_post(2) / 63.0d0 - &
             m_post(3) / 10.0d0 + m_post(4) / 10.0d0 + m_post(9) / 18.0d0 - m_post(10) / 18.0d0

        f_post(i,j,k,3) =  m_post(0) / 19.0d0 - 11.0d0 * m_post(1) / 2394.0d0 - m_post(2) / 63.0d0 + &
             m_post(5) / 10.0d0 - m_post(6) / 10.0d0 - m_post(9) / 36.0d0 + m_post(10) / 36.0d0 + &
             m_post(11) / 12.0d0 - m_post(12) / 12.0d0

        f_post(i,j,k,4) =  m_post(0) / 19.0d0 - 11.0d0 * m_post(1) / 2394.0d0 - m_post(2) / 63.0d0 - &
             m_post(5) / 10.0d0 + m_post(6) / 10.0d0 - m_post(9) / 36.0d0 + m_post(10) / 36.0d0 + &
             m_post(11) / 12.0d0 - m_post(12) / 12.0d0

        f_post(i,j,k,5) =  m_post(0) / 19.0d0 - 11.0d0 * m_post(1) / 2394.0d0 - m_post(2) / 63.0d0 + &
             m_post(7) / 10.0d0 - m_post(8) / 10.0d0 - m_post(9) / 36.0d0 + m_post(10) / 36.0d0 - &
             m_post(11) / 12.0d0 + m_post(12) / 12.0d0

        f_post(i,j,k,6) =  m_post(0) / 19.0d0 - 11.0d0 * m_post(1) / 2394.0d0 - m_post(2) / 63.0d0 - &
             m_post(7) / 10.0d0 + m_post(8) / 10.0d0 - m_post(9) / 36.0d0 + m_post(10) / 36.0d0 - &
             m_post(11) / 12.0d0 + m_post(12) / 12.0d0

        f_post(i,j,k,7) =  m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 + &
             m_post(3) / 10.0d0 + m_post(4) / 40.0d0 + m_post(5) / 10.0d0 + m_post(6) / 40.0d0 + &
             m_post(9) / 36.0d0 + m_post(10) / 72.0d0 + m_post(11) / 12.0d0 + m_post(12) / 24.0d0 + &
             m_post(13) / 4.0d0 + m_post(16) / 8.0d0 - m_post(17) / 8.0d0

        f_post(i,j,k,8) =  m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 - &
             m_post(3) / 10.0d0 - m_post(4) / 40.0d0 + m_post(5) / 10.0d0 + m_post(6) / 40.0d0 + &
             m_post(9) / 36.0d0 + m_post(10) / 72.0d0 + m_post(11) / 12.0d0 + m_post(12) / 24.0d0 - &
             m_post(13) / 4.0d0 - m_post(16) / 8.0d0 - m_post(17) / 8.0d0

        f_post(i,j,k,9) =  m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 + &
             m_post(3) / 10.0d0 + m_post(4) / 40.0d0 - m_post(5) / 10.0d0 - m_post(6) / 40.0d0 + &
             m_post(9) / 36.0d0 + m_post(10) / 72.0d0 + m_post(11) / 12.0d0 + m_post(12) / 24.0d0 - &
             m_post(13) / 4.0d0 + m_post(16) / 8.0d0 + m_post(17) / 8.0d0

        f_post(i,j,k,10) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 - &
             m_post(3) / 10.0d0 - m_post(4) / 40.0d0 - m_post(5) / 10.0d0 - m_post(6) / 40.0d0 + &
             m_post(9) / 36.0d0 + m_post(10) / 72.0d0 + m_post(11) / 12.0d0 + m_post(12) / 24.0d0 + &
             m_post(13) / 4.0d0 - m_post(16) / 8.0d0 + m_post(17) / 8.0d0

        f_post(i,j,k,11) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 + &
             m_post(3) / 10.0d0 + m_post(4) / 40.0d0 + m_post(7) / 10.0d0 + m_post(8) / 40.0d0 + &
             m_post(9) / 36.0d0 + m_post(10) / 72.0d0 - m_post(11) / 12.0d0 - m_post(12) / 24.0d0 + &
             m_post(15) / 4.0d0 - m_post(16) / 8.0d0 + m_post(18) / 8.0d0

        f_post(i,j,k,12) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 - &
             m_post(3) / 10.0d0 - m_post(4) / 40.0d0 + m_post(7) / 10.0d0 + m_post(8) / 40.0d0 + &
             m_post(9) / 36.0d0 + m_post(10) / 72.0d0 - m_post(11) / 12.0d0 - m_post(12) / 24.0d0 - &
             m_post(15) / 4.0d0 + m_post(16) / 8.0d0 + m_post(18) / 8.0d0

        f_post(i,j,k,13) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 + &
             m_post(3) / 10.0d0 + m_post(4) / 40.0d0 - m_post(7) / 10.0d0 - m_post(8) / 40.0d0 + &
             m_post(9) / 36.0d0 + m_post(10) / 72.0d0 - m_post(11) / 12.0d0 - m_post(12) / 24.0d0 - &
             m_post(15) / 4.0d0 - m_post(16) / 8.0d0 - m_post(18) / 8.0d0

        f_post(i,j,k,14) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 - &
             m_post(3) / 10.0d0 - m_post(4) / 40.0d0 - m_post(7) / 10.0d0 - m_post(8) / 40.0d0 + &
             m_post(9) / 36.0d0 + m_post(10) / 72.0d0 - m_post(11) / 12.0d0 - m_post(12) / 24.0d0 + &
             m_post(15) / 4.0d0 + m_post(16) / 8.0d0 - m_post(18) / 8.0d0

        f_post(i,j,k,15) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 + &
             m_post(5) / 10.0d0 + m_post(6) / 40.0d0 + m_post(7) / 10.0d0 + m_post(8) / 40.0d0 - &
             m_post(9) / 18.0d0 - m_post(10) / 36.0d0 + m_post(14) / 4.0d0 + m_post(17) / 8.0d0 - &
             m_post(18) / 8.0d0

        f_post(i,j,k,16) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 - &
             m_post(5) / 10.0d0 - m_post(6) / 40.0d0 + m_post(7) / 10.0d0 + m_post(8) / 40.0d0 - &
             m_post(9) / 18.0d0 - m_post(10) / 36.0d0 - m_post(14) / 4.0d0 - m_post(17) / 8.0d0 - &
             m_post(18) / 8.0d0

        f_post(i,j,k,17) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 + &
             m_post(5) / 10.0d0 + m_post(6) / 40.0d0 - m_post(7) / 10.0d0 - m_post(8) / 40.0d0 - &
             m_post(9) / 18.0d0 - m_post(10) / 36.0d0 - m_post(14) / 4.0d0 + m_post(17) / 8.0d0 + &
             m_post(18) / 8.0d0

        f_post(i,j,k,18) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 - &
             m_post(5) / 10.0d0 - m_post(6) / 40.0d0 - m_post(7) / 10.0d0 - m_post(8) / 40.0d0 - &
             m_post(9) / 18.0d0 - m_post(10) / 36.0d0 + m_post(14) / 4.0d0 - m_post(17) / 8.0d0 + &
             m_post(18) / 8.0d0
      enddo
    enddo
  enddo

end subroutine collision3d


!===========================================================================================================================
! 子程序: streaming3d
! 作用: 对流场分布函数执行三维 pull streaming。
!===========================================================================================================================
subroutine streaming3d()
  use commondata3dOpenacc
  implicit none

  integer(kind=4) :: i, j, k, ip, jp, kp, alpha

  ! pull streaming：当前格点 (i,j,k) 从上游格点 (i-ex, j-ey, k-ez) 拉取分布函数
 !$acc parallel loop gang vector collapse(3) default(none) present(f,f_post,ex,ey,ez) async(1) private(alpha,ip,jp,kp)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        do alpha = 0, qf-1
          ip = i - ex(alpha)
          jp = j - ey(alpha)
          kp = k - ez(alpha)
          f(i,j,k,alpha) = f_post(ip,jp,kp,alpha)
        enddo
      enddo
    enddo
  enddo

end subroutine streaming3d


!===========================================================================================================================
! 子程序: bounceback3d
! 作用: 施加流场边界条件，包括无滑移壁面和周期边界配套处理。
!===========================================================================================================================
subroutine bounceback3d()
  use commondata3dOpenacc
  implicit none

  integer(kind=4) :: i, j, k

#ifdef SpanwiseWallsPeriodicalU
  !$acc parallel loop gang vector collapse(2) default(none) present(f,f_post) async(1)
  do k = 1, nz
    do j = 1, ny
      ! 前边界 (i=1)：处理 ex=+1 的入射分布函数
      f(1,j,k,1) = f_post(nx,j,k,1)
      f(1,j,k,7) = f_post(nx,j,k,7)
      f(1,j,k,9) = f_post(nx,j,k,9)
      f(1,j,k,11) = f_post(nx,j,k,11)
      f(1,j,k,13) = f_post(nx,j,k,13)

      ! 后边界 (i=nx)：处理 ex=-1 的入射分布函数
      f(nx,j,k,2) = f_post(1,j,k,2)
      f(nx,j,k,8) = f_post(1,j,k,8)
      f(nx,j,k,10) = f_post(1,j,k,10)
      f(nx,j,k,12) = f_post(1,j,k,12)
      f(nx,j,k,14) = f_post(1,j,k,14)
    enddo
  enddo
#endif

#ifdef SpanwiseWallsNoslip
 !$acc parallel loop gang vector collapse(2) default(none) present(f,f_post) async(1)
  do k = 1, nz
    do j = 1, ny
      ! 前侧无滑移壁面 (i=1)：处理 ex=+1 的入射分布函数
      f(1,j,k,1) = f_post(1,j,k,2)
      f(1,j,k,7) = f_post(1,j,k,10)
      f(1,j,k,9) = f_post(1,j,k,8)
      f(1,j,k,11) = f_post(1,j,k,14)
      f(1,j,k,13) = f_post(1,j,k,12)

      ! 后侧无滑移壁面 (i=nx)：处理 ex=-1 的入射分布函数
      f(nx,j,k,2) = f_post(nx,j,k,1)
      f(nx,j,k,8) = f_post(nx,j,k,9)
      f(nx,j,k,10) = f_post(nx,j,k,7)
      f(nx,j,k,12) = f_post(nx,j,k,13)
      f(nx,j,k,14) = f_post(nx,j,k,11)
    enddo
  enddo
#endif

#ifdef VerticalWallsPeriodicalU
 !$acc parallel loop gang vector collapse(2) default(none) present(f,f_post) async(1)
  do k = 1, nz
    do i = 1, nx
      ! 左边界 (j=1)：处理 ey=+1 的入射分布函数
      f(i,1,k,3) = f_post(i,ny,k,3)
      f(i,1,k,7) = f_post(i,ny,k,7)
      f(i,1,k,8) = f_post(i,ny,k,8)
      f(i,1,k,15) = f_post(i,ny,k,15)
      f(i,1,k,17) = f_post(i,ny,k,17)

      ! 右边界 (j=ny)：处理 ey=-1 的入射分布函数
      f(i,ny,k,4) = f_post(i,1,k,4)
      f(i,ny,k,9) = f_post(i,1,k,9)
      f(i,ny,k,10) = f_post(i,1,k,10)
      f(i,ny,k,16) = f_post(i,1,k,16)
      f(i,ny,k,18) = f_post(i,1,k,18)
    enddo
  enddo
#endif

#ifdef VerticalWallsNoslip
 !$acc parallel loop gang vector collapse(2) default(none) present(f,f_post) async(1)
  do k = 1, nz
    do i = 1, nx
      ! 左侧无滑移壁面 (j=1)：处理 ey=+1 的入射分布函数
      f(i,1,k,3) = f_post(i,1,k,4)
      f(i,1,k,7) = f_post(i,1,k,10)
      f(i,1,k,8) = f_post(i,1,k,9)
      f(i,1,k,15) = f_post(i,1,k,18)
      f(i,1,k,17) = f_post(i,1,k,16)

      ! 右侧无滑移壁面 (j=ny)：处理 ey=-1 的入射分布函数
      f(i,ny,k,4) = f_post(i,ny,k,3)
      f(i,ny,k,9) = f_post(i,ny,k,8)
      f(i,ny,k,10) = f_post(i,ny,k,7)
      f(i,ny,k,16) = f_post(i,ny,k,15)
      f(i,ny,k,18) = f_post(i,ny,k,17)
    enddo
  enddo
#endif

#ifdef HorizontalWallsPeriodicalU
 !$acc parallel loop gang vector collapse(2) default(none) present(f,f_post) async(1)
  do j = 1, ny
    do i = 1, nx
      ! 下边界 (k=1)：处理 ez=+1 的入射分布函数
      f(i,j,1,5) = f_post(i,j,nz,5)
      f(i,j,1,11) = f_post(i,j,nz,11)
      f(i,j,1,12) = f_post(i,j,nz,12)
      f(i,j,1,15) = f_post(i,j,nz,15)
      f(i,j,1,16) = f_post(i,j,nz,16)

      ! 上边界 (k=nz)：处理 ez=-1 的入射分布函数
      f(i,j,nz,6) = f_post(i,j,1,6)
      f(i,j,nz,13) = f_post(i,j,1,13)
      f(i,j,nz,14) = f_post(i,j,1,14)
      f(i,j,nz,17) = f_post(i,j,1,17)
      f(i,j,nz,18) = f_post(i,j,1,18)
    enddo
  enddo
#endif

#ifdef HorizontalWallsNoslip
 !$acc parallel loop gang vector collapse(2) default(none) present(f,f_post) async(1)
  do j = 1, ny
    do i = 1, nx
      ! 下侧无滑移壁面 (k=1)：处理 ez=+1 的入射分布函数
      f(i,j,1,5) = f_post(i,j,1,6)
      f(i,j,1,11) = f_post(i,j,1,14)
      f(i,j,1,12) = f_post(i,j,1,13)
      f(i,j,1,15) = f_post(i,j,1,18)
      f(i,j,1,16) = f_post(i,j,1,17)

      ! 上侧无滑移壁面 (k=nz)：处理 ez=-1 的入射分布函数
      f(i,j,nz,6) = f_post(i,j,nz,5)
      f(i,j,nz,13) = f_post(i,j,nz,12)
      f(i,j,nz,14) = f_post(i,j,nz,11)
      f(i,j,nz,17) = f_post(i,j,nz,16)
      f(i,j,nz,18) = f_post(i,j,nz,15)
    enddo
  enddo
#endif

end subroutine bounceback3d


!===========================================================================================================================
! 子程序: macro3d
! 作用: 由流场分布函数恢复 rho、u、v、w 以及浮力项。
!===========================================================================================================================
subroutine macro3d()
  use commondata3dOpenacc
  implicit none

  integer(kind=4) :: i, j, k
  real(kind=8) :: FzLoc


 !$acc parallel loop gang vector collapse(3) default(none) present(f,rho,u,v,w,T) async(1) private(FzLoc)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        rho(i,j,k) = f(i,j,k,0) + f(i,j,k,1) + f(i,j,k,2) + f(i,j,k,3) + &
             f(i,j,k,4) + f(i,j,k,5) + f(i,j,k,6) + f(i,j,k,7) + &
             f(i,j,k,8) + f(i,j,k,9) + f(i,j,k,10) + f(i,j,k,11) + &
             f(i,j,k,12) + f(i,j,k,13) + f(i,j,k,14) + f(i,j,k,15) + &
             f(i,j,k,16) + f(i,j,k,17) + f(i,j,k,18)

        FzLoc = rho(i,j,k) * gBeta * (T(i,j,k) - Tref)

        u(i,j,k) = ( f(i,j,k,1) - f(i,j,k,2) + f(i,j,k,7) - f(i,j,k,8) + f(i,j,k,9) - f(i,j,k,10) + &
             f(i,j,k,11) - f(i,j,k,12) + f(i,j,k,13) - f(i,j,k,14) ) / rho(i,j,k)

        v(i,j,k) = ( f(i,j,k,3) - f(i,j,k,4) + f(i,j,k,7) + f(i,j,k,8) - f(i,j,k,9) - f(i,j,k,10) + &
             f(i,j,k,15) - f(i,j,k,16) + f(i,j,k,17) - f(i,j,k,18) ) / rho(i,j,k)

        w(i,j,k) = ( f(i,j,k,5) - f(i,j,k,6) + f(i,j,k,11) + f(i,j,k,12) - f(i,j,k,13) - f(i,j,k,14) + &
             f(i,j,k,15) + f(i,j,k,16) - f(i,j,k,17) - f(i,j,k,18) + 0.5d0 * FzLoc ) / rho(i,j,k)
      enddo
    enddo
  enddo

end subroutine macro3d


!===========================================================================================================================
! 子程序: collisionT3d
! 作用: 温度场碰撞步骤，算法口径尽量保持与 2DRB 的对流扩散处理一致。
!===========================================================================================================================
subroutine collisionT3d()
  use commondata3dOpenacc
  implicit none

  integer(kind=4) :: i, j, k
  real(kind=8) :: n(0:qt-1), neq(0:qt-1), q(0:qt-1), n_post(0:qt-1)
  real(kind=8) :: Bx, By, Bz, dBx, dBy, dBz
  real(kind=8), parameter :: SG = 1.0d0 - 0.5d0 * Qk

  ! 温度场采用 D3Q7 MRT
  !$acc parallel loop gang vector collapse(3) default(none) present(g,g_post,u,v,w,T,Bx_prev,By_prev,Bz_prev) async(1) &
  !$acc& private(n,neq,q,n_post,Bx,By,Bz,dBx,dBy,dBz)
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
        n(0) = g(i,j,k,0) + g(i,j,k,1) + g(i,j,k,2) + g(i,j,k,3) + g(i,j,k,4) + g(i,j,k,5) + g(i,j,k,6)
        n(1) = g(i,j,k,1) - g(i,j,k,2)
        n(2) = g(i,j,k,3) - g(i,j,k,4)
        n(3) = g(i,j,k,5) - g(i,j,k,6)
        n(4) = -6.0d0 * g(i,j,k,0) + g(i,j,k,1) + g(i,j,k,2) + g(i,j,k,3) + g(i,j,k,4) + g(i,j,k,5) + g(i,j,k,6)
        n(5) = 2.0d0 * g(i,j,k,1) + 2.0d0 * g(i,j,k,2) - g(i,j,k,3) - g(i,j,k,4) - g(i,j,k,5) - g(i,j,k,6)
        n(6) = g(i,j,k,3) + g(i,j,k,4) - g(i,j,k,5) - g(i,j,k,6)

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
        g_post(i,j,k,0) = n_post(0) / 7.0d0 - n_post(4) / 7.0d0
        g_post(i,j,k,1) = n_post(0) / 7.0d0 + n_post(1) / 2.0d0 + n_post(4) / 42.0d0 + n_post(5) / 6.0d0
        g_post(i,j,k,2) = n_post(0) / 7.0d0 - n_post(1) / 2.0d0 + n_post(4) / 42.0d0 + n_post(5) / 6.0d0
        g_post(i,j,k,3) = n_post(0) / 7.0d0 + n_post(2) / 2.0d0 + n_post(4) / 42.0d0 - n_post(5) / 12.0d0 + &
             n_post(6) / 4.0d0
        g_post(i,j,k,4) = n_post(0) / 7.0d0 - n_post(2) / 2.0d0 + n_post(4) / 42.0d0 - n_post(5) / 12.0d0 + &
             n_post(6) / 4.0d0
        g_post(i,j,k,5) = n_post(0) / 7.0d0 + n_post(3) / 2.0d0 + n_post(4) / 42.0d0 - n_post(5) / 12.0d0 - &
             n_post(6) / 4.0d0
        g_post(i,j,k,6) = n_post(0) / 7.0d0 - n_post(3) / 2.0d0 + n_post(4) / 42.0d0 - n_post(5) / 12.0d0 - &
             n_post(6) / 4.0d0
      enddo
    enddo
  enddo

end subroutine collisionT3d


!===========================================================================================================================
! 子程序: streamingT3d
! 作用: 对温度分布函数执行三维 pull streaming。
!===========================================================================================================================
subroutine streamingT3d()
  use commondata3dOpenacc
  implicit none

  integer(kind=4) :: i, j, k, ip, jp, kp, alpha

  ! 温度场同样采用 pull streaming
 !$acc parallel loop gang vector collapse(3) default(none) present(g,g_post,exT,eyT,ezT) async(1) private(alpha,ip,jp,kp)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        do alpha = 0, qt-1
          ip = i - exT(alpha)
          jp = j - eyT(alpha)
          kp = k - ezT(alpha)
          g(i,j,k,alpha) = g_post(ip,jp,kp,alpha)
        enddo
      enddo
    enddo
  enddo

  return
end subroutine streamingT3d


!===========================================================================================================================
! 子程序: bouncebackT3d
! 作用: 施加温度边界条件，包括恒温、绝热和周期边界。
!===========================================================================================================================
subroutine bouncebackT3d()
  use commondata3dOpenacc
  implicit none

  integer(kind=4) :: i, j, k

#ifdef SpanwiseWallsPeriodicalT
  !$acc parallel loop gang vector collapse(2) default(none) present(g,g_post) async(1)
  do k = 1, nz
    do j = 1, ny
      g(1,j,k,1)  = g_post(nx,j,k,1)
      g(nx,j,k,2) = g_post(1,j,k,2)
    enddo
  enddo
#endif

#ifdef VerticalWallsPeriodicalT
 !$acc parallel loop gang vector collapse(2) default(none) present(g,g_post) async(1)
  do k = 1, nz
    do i = 1, nx
      g(i,1,k,3)  = g_post(i,ny,k,3)
      g(i,ny,k,4) = g_post(i,1,k,4)
    enddo
  enddo
#endif

#ifdef VerticalWallsConstT
  !$acc parallel loop gang vector collapse(2) default(none) present(g,g_post,omegaT) async(1)
  do k = 1, nz
    do i = 1, nx
#ifdef EnableLegacyThermalScheme
      g(i,1,k,3)  = -g_post(i,1,k,4)  + (6.0d0 + paraA) / 21.0d0 * Thot
      g(i,ny,k,4) = -g_post(i,ny,k,3) + (6.0d0 + paraA) / 21.0d0 * Tcold
#else
      g(i,1,k,3)  = -g_post(i,1,k,4)  + 2.0d0 * omegaT(3) * Thot
      g(i,ny,k,4) = -g_post(i,ny,k,3) + 2.0d0 * omegaT(4) * Tcold
#endif
    enddo
  enddo
#endif

#ifdef SpanwiseWallsAdiabatic
 !$acc parallel loop gang vector collapse(2) default(none) present(g,g_post) async(1)
  do k = 1, nz
    do j = 1, ny
      g(1,j,k,1)  = g_post(1,j,k,2)
      g(nx,j,k,2) = g_post(nx,j,k,1)
    enddo
  enddo
#endif

#ifdef VerticalWallsAdiabatic
 !$acc parallel loop gang vector collapse(2) default(none) present(g,g_post) async(1)
  do k = 1, nz
    do i = 1, nx
      g(i,1,k,3)  = g_post(i,1,k,4)
      g(i,ny,k,4) = g_post(i,ny,k,3)
    enddo
  enddo
#endif

#ifdef HorizontalWallsConstT
 !$acc parallel loop gang vector collapse(2) default(none) present(g,g_post,omegaT) async(1)
  do j = 1, ny
    do i = 1, nx
#ifdef EnableLegacyThermalScheme
      g(i,j,1,5)  = -g_post(i,j,1,6)  + (6.0d0 + paraA) / 21.0d0 * Thot
      g(i,j,nz,6) = -g_post(i,j,nz,5) + (6.0d0 + paraA) / 21.0d0 * Tcold
#else
      g(i,j,1,5)  = -g_post(i,j,1,6)  + 2.0d0 * omegaT(5) * Thot
      g(i,j,nz,6) = -g_post(i,j,nz,5) + 2.0d0 * omegaT(6) * Tcold
#endif
    enddo
  enddo
#endif

#ifdef HorizontalWallsPeriodicalT
 !$acc parallel loop gang vector collapse(2) default(none) present(g,g_post) async(1)
  do j = 1, ny
    do i = 1, nx
      g(i,j,1,5)  = g_post(i,j,nz,5)
      g(i,j,nz,6) = g_post(i,j,1,6)
    enddo
  enddo
#endif

#ifdef HorizontalWallsAdiabatic
 !$acc parallel loop gang vector collapse(2) default(none) present(g,g_post) async(1)
  do j = 1, ny
    do i = 1, nx
      g(i,j,1,5)  = g_post(i,j,1,6)
      g(i,j,nz,6) = g_post(i,j,nz,5)
    enddo
  enddo
#endif

  return
end subroutine bouncebackT3d


!===========================================================================================================================
! 子程序: macroT3d
! 作用: 由温度分布函数恢复温度场，并更新历史热流项。
!===========================================================================================================================
subroutine macroT3d()
  use commondata3dOpenacc
  implicit none

  integer(kind=4) :: i, j, k

  ! 温度恢复就是对 7 个方向的 g 求和
 !$acc parallel loop gang vector collapse(3) default(none) present(g,T) async(1)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        T(i,j,k) = g(i,j,k,0) + g(i,j,k,1) + g(i,j,k,2) + g(i,j,k,3) + &
                   g(i,j,k,4) + g(i,j,k,5) + g(i,j,k,6)
      enddo
    enddo
  enddo

  return
end subroutine macroT3d


!===========================================================================================================================
! 子程序: reconstruct_macro_from_fg3d
! 作用: 从重启读回的 f/g 重新恢复宏观场。
!===========================================================================================================================
subroutine reconstruct_macro_from_fg3d()
  use commondata3dOpenacc
  implicit none

  integer(kind=4) :: i, j, k, iter
  real(kind=8) :: momx, momy, momz
  real(kind=8) :: FxLoc, FyLoc, FzLoc
  logical :: rho_bad

  ! 严格重启文件会保存 EnableUseG 的历史热流；这里只从 f/g 恢复当前宏观场。
  rho_bad = .false.
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        T(i,j,k) = g(i,j,k,0) + g(i,j,k,1) + g(i,j,k,2) + g(i,j,k,3) + &
                   g(i,j,k,4) + g(i,j,k,5) + g(i,j,k,6)

        rho(i,j,k) = f(i,j,k,0) + f(i,j,k,1) + f(i,j,k,2) + f(i,j,k,3) + &
             f(i,j,k,4) + f(i,j,k,5) + f(i,j,k,6) + f(i,j,k,7) + &
             f(i,j,k,8) + f(i,j,k,9) + f(i,j,k,10) + f(i,j,k,11) + &
             f(i,j,k,12) + f(i,j,k,13) + f(i,j,k,14) + f(i,j,k,15) + &
             f(i,j,k,16) + f(i,j,k,17) + f(i,j,k,18)

        momx = f(i,j,k,1) - f(i,j,k,2) + f(i,j,k,7) - f(i,j,k,8) + f(i,j,k,9) - f(i,j,k,10) + &
             f(i,j,k,11) - f(i,j,k,12) + f(i,j,k,13) - f(i,j,k,14)

        momy = f(i,j,k,3) - f(i,j,k,4) + f(i,j,k,7) + f(i,j,k,8) - f(i,j,k,9) - f(i,j,k,10) + &
             f(i,j,k,15) - f(i,j,k,16) + f(i,j,k,17) - f(i,j,k,18)

        momz = f(i,j,k,5) - f(i,j,k,6) + f(i,j,k,11) + f(i,j,k,12) - f(i,j,k,13) - f(i,j,k,14) + &
             f(i,j,k,15) + f(i,j,k,16) - f(i,j,k,17) - f(i,j,k,18)

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

  if (rho_bad) then
    write(*,*) 'Warning: non-positive rho found during restart reconstruction.'
    stop
  endif
end subroutine reconstruct_macro_from_fg3d


#ifdef steadyFlow
!===========================================================================================================================
! 子程序: check3d
! 作用: 计算稳态收敛误差，并按需写入收敛历史。
!===========================================================================================================================
subroutine check3d()
  use commondata3dOpenacc
  implicit none

  integer(kind=4) :: i, j, k
  real(kind=8) :: error1, error2, error5, error6
  character(len=80) :: caseTag

  ! 误差定义：errorU 用速度场相对 L2 误差
  ! errorT 用温度场相对 L1 误差
  !$acc wait(1)
  error1 = 0.0d0
  error2 = 0.0d0
  error5 = 0.0d0
  error6 = 0.0d0

  ! reduction(+:...) : 每个线程先做局部累加，kernel 结束后再安全归并成 error1/error2/error5/error6。
 !$acc parallel loop collapse(3) default(none) &
 !$acc& present(u,up,v,vp,w,wp,T,Tp) &
 !$acc& reduction(+:error1,error2,error5,error6)
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

  call append_convergence_tecplot3d('convergence3D.plt', restartItcOffset+itc, errorU, errorT)

  write(caseTag,'("Ra=",ES10.3E2,",nx=",I0,",ny=",I0,",nz=",I0,",useG=",L1,",old=",L1)') &
       real(Rayleigh,kind=8), nx, ny, nz, useG, useLegacyThermalScheme
  call append_convergence_master_tecplot3d('convergence_all_3D.plt', caseTag, restartItcOffset+itc, errorU, errorT)
  write(*,'(I12,1X,ES24.16E3,1X,ES24.16E3)') restartItcOffset+itc, real(errorU,kind=8), real(errorT,kind=8)

end subroutine check3d
#endif


!===========================================================================================================================
! 子程序: append_convergence_tecplot3d
! 作用: 向单个收敛历史文件追加一条误差记录。
!===========================================================================================================================
subroutine append_convergence_tecplot3d(filename, itcLoc, errorULoc, errorTLoc)
  use commondata3dOpenacc, only: loadInitField
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
    write(u,'(I12,1X,ES24.16E3,1X,ES24.16E3)') itcLoc, real(errorULoc,kind=8), real(errorTLoc,kind=8)
    close(u)
    first_write = .false.
  else
    ! 同一次运行的后续调用：追加数据行
    open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
    write(u,'(I12,1X,ES24.16E3,1X,ES24.16E3)') itcLoc, real(errorULoc,kind=8), real(errorTLoc,kind=8)
    close(u)
  endif

end subroutine append_convergence_tecplot3d


!===========================================================================================================================
! 子程序: append_convergence_master_tecplot3d
! 作用: 向带 zone 名称的收敛历史文件追加一条记录。
!===========================================================================================================================
subroutine append_convergence_master_tecplot3d(filename, zoneName, itcLoc, errorULoc, errorTLoc)
  use commondata3dOpenacc, only: loadInitField
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
  write(u,'(I12,1X,ES24.16E3,1X,ES24.16E3)') itcLoc, real(errorULoc,kind=8), real(errorTLoc,kind=8)
  close(u)

end subroutine append_convergence_master_tecplot3d


!===========================================================================================================================
! 子程序: output_SnapshotFile3d
! 作用: 输出三维快照二进制文件，供后处理或继续分析使用。
!===========================================================================================================================
subroutine output_SnapshotFile3d()
  use commondata3dOpenacc
  implicit none

  integer(kind=4) :: i, j, k
  character(len=100) :: filename


#ifdef steadyFlow
  write(filename,'(i12.12)') restartItcOffset+itc
#endif
#ifdef unsteadyFlow
  snapshotFileNum = snapshotFileNum + 1
  write(filename,'(i12.12)') snapshotFileNum
#endif

  filename = adjustl(filename)
  open(unit=03, file=trim(snapshotFilePrefix)//'-'//trim(filename)//'.bin', form='unformatted', access='sequential')
  write(03) (((real(velocityScaleCompare*u(i,j,k),kind=8), i=1,nx), j=1,ny), k=1,nz)
  write(03) (((real(velocityScaleCompare*v(i,j,k),kind=8), i=1,nx), j=1,ny), k=1,nz)
  write(03) (((real(velocityScaleCompare*w(i,j,k),kind=8), i=1,nx), j=1,ny), k=1,nz)
  write(03) (((real(T(i,j,k),kind=8), i=1,nx), j=1,ny), k=1,nz)
  write(03) (((real(rho(i,j,k),kind=8), i=1,nx), j=1,ny), k=1,nz)
  close(03)

  return
end subroutine output_SnapshotFile3d


!===========================================================================================================================
! 子程序: output_ReloadFile3d
! 作用: 备份严格重启状态，供后续重启继续计算。
!===========================================================================================================================
subroutine output_ReloadFile3d()
  use commondata3dOpenacc
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
  write(05) ((((real(f(i,j,k,alpha),kind=8), i=1,nx), j=1,ny), k=1,nz), alpha=0,qf-1)
  write(05) ((((real(g(i,j,k,alpha),kind=8), i=1,nx), j=1,ny), k=1,nz), alpha=0,qt-1)
#ifdef EnableUseG
  write(05) (((real(Bx_prev(i,j,k),kind=8), i=1,nx), j=1,ny), k=1,nz)
  write(05) (((real(By_prev(i,j,k),kind=8), i=1,nx), j=1,ny), k=1,nz)
  write(05) (((real(Bz_prev(i,j,k),kind=8), i=1,nx), j=1,ny), k=1,nz)
#endif
  close(05)
  call write_reload_metadata3d(trim(filename))

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) 'Backup f/g restart state to: ', trim(reloadFilePrefix), '-', trim(filename), '.bin'
  write(00,*) 'Backup restart metadata to: ', trim(reloadFilePrefix), '-latest.meta'
  close(00)

  return
end subroutine output_ReloadFile3d


!===========================================================================================================================
! 子程序: write_reload_metadata3d
! 作用: 覆盖写出最新 reload 续算账本，恢复累计步数、t_ff、输出编号和最新 .bin 文件名。
!===========================================================================================================================
subroutine write_reload_metadata3d(filename)
  use commondata3dOpenacc
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4) :: metaUnit, totalItc
  real(kind=8) :: totalTf

  totalItc = restartItcOffset + itc
  totalTf = real(totalItc,kind=8) / timeUnit

  open(newunit=metaUnit, file=trim(reloadFilePrefix)//'-latest.meta', &
       status='replace', action='write', form='formatted')
  write(metaUnit,'(A,1X,I0)') 'reload_meta_version', 4
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
end subroutine write_reload_metadata3d


!===========================================================================================================================
! 子程序: read_reload_metadata3d
! 作用: 优先读取 latest .meta；若没有，则根据手工编号做保守推断。
!===========================================================================================================================
subroutine read_reload_metadata3d(reloadFileName)
  use commondata3dOpenacc
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
    call infer_reload_offsets_without_metadata3d()
    return
  endif

  open(newunit=metaUnit, file=trim(metaFile), status='old', action='read', &
       form='formatted', iostat=ios)                             ! ios==0 表示成功，非 0 表示打开失败
  if (ios .NE. 0) then
    write(*,*) 'Error: failed to open reload metadata: ', trim(metaFile)
    stop
  endif

  read(metaUnit,*,iostat=ios) label, metaVersion
  if ((ios .NE. 0) .OR. (trim(label) .NE. 'reload_meta_version') .OR. (metaVersion .NE. 4)) then
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
  ! reloadFileNum 是整数计数器，给后续 output_ReloadFile3d() 继续编号，避免覆盖旧 reload 文件。
  ! reloadFileName 是字符串文件名，本次续算马上用它打开 <reloadFilePrefix>-<reloadFileName>.bin。
  reloadFileNum = metaReloadFileNum
  reloadFileName = trim(metaReloadFileName)
  reloadMetadataLoaded = .true.

  return
end subroutine read_reload_metadata3d


!===========================================================================================================================
! 子程序: infer_reload_offsets_without_metadata3d
! 作用: 没有 latest .meta 时，只能根据文件编号和当前手工参数推断。
! 根据文件名编号和当前参数“猜一个合理值”，保证续算的时间/步数尽量连续。
!===========================================================================================================================
subroutine infer_reload_offsets_without_metadata3d()
  use commondata3dOpenacc
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
end subroutine infer_reload_offsets_without_metadata3d


!===========================================================================================================================
! 子程序: output_Tecplot3d
! 作用: 输出全场体数据以及 x/y/z 三个中面切片，便于快速查看三维流场结构。
!===========================================================================================================================
subroutine output_Tecplot3d()
  use commondata3dOpenacc
  implicit none

  character(len=100) :: filename

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
end subroutine output_Tecplot3d


!===========================================================================================================================
! 子程序: calNuRe3d
! 作用: 计算体平均 Nu / Re，并把时间序列缓存到数组中。
!===========================================================================================================================
subroutine calNuRe3d()
  use commondata3dOpenacc
  implicit none

  integer(kind=4) :: i, j, k
  real(kind=8) :: NuVolAvg_temp, ReVolAvg_temp
  real(kind=8) :: sampleTime
  logical :: exNu, exRe
  logical, save :: first_nure_write = .true.

  !$acc wait(1)
  ! 时间序列的体平均 Nu / Re：
  ! NuVolAvg : 体平均对流热通量对应的 Nu
  ! ReVolAvg : 全域 RMS 速度对应的 Reynolds 数
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
 !$acc parallel loop collapse(3) default(none) present(v,T) reduction(+:NuVolAvg_temp)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        NuVolAvg_temp = NuVolAvg_temp + v(i,j,k) * T(i,j,k)
      enddo
    enddo
  enddo
#else
 !$acc parallel loop collapse(3) default(none) present(w,T) reduction(+:NuVolAvg_temp)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        NuVolAvg_temp = NuVolAvg_temp + w(i,j,k) * T(i,j,k)
      enddo
    enddo
  enddo
#endif

  NuVolAvg(dimensionlessTime) = NuVolAvg_temp / dble(nx * ny * nz) * velocityScaleCompare + 1.0d0
  if ((first_nure_write) .AND. (loadInitField .EQ. 0)) then
    open(unit=01, file='Nu_VolAvg_3D.dat', status='replace', action='write')
  else
    open(unit=01, file='Nu_VolAvg_3D.dat', status='unknown', position='append', action='write')
  endif
  write(01,'(2(ES24.16E3,1X))') &
       real(sampleTime, kind=8), &
       real(NuVolAvg(dimensionlessTime), kind=8)
  close(01)

  ReVolAvg_temp = 0.0d0
 !$acc parallel loop collapse(3) default(none) present(u,v,w) reduction(+:ReVolAvg_temp)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        ReVolAvg_temp = ReVolAvg_temp + u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)
      enddo
    enddo
  enddo
  ReVolAvg(dimensionlessTime) = dsqrt(ReVolAvg_temp / dble(nx * ny * nz)) * lengthUnit / viscosity
  if ((first_nure_write) .AND. (loadInitField .EQ. 0)) then
    open(unit=02, file='Re_VolAvg_3D.dat', status='replace', action='write')
  else
    open(unit=02, file='Re_VolAvg_3D.dat', status='unknown', position='append', action='write')
  endif
  write(02,'(2(ES24.16E3,1X))') &
       real(sampleTime, kind=8), &
       real(ReVolAvg(dimensionlessTime), kind=8)
  close(02)
  first_nure_write = .false.

  write(*,'(a,1x,ES24.16E3)') 'NuVolAvg =', real(NuVolAvg(dimensionlessTime),kind=8)
  write(*,'(a,1x,ES24.16E3)') 'ReVolAvg =', real(ReVolAvg(dimensionlessTime),kind=8)

  return
end subroutine calNuRe3d


#ifdef unsteadyFlow
!===========================================================================================================================
! 子程序: output_unsteady_NuRe_postprocess3d
! 作用: 从完整 .dat 历史重建非稳态 Nu/Re 序列、运行平均和窗口平均。
!===========================================================================================================================
subroutine output_unsteady_NuRe_postprocess3d()
  use commondata3dOpenacc
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

  ! 这些文件是完整 .dat 历史的派生视图，因此重建为一个连续的 ZONE。
  open(newunit=seriesUnit, file='NuRe_VolAvg_3DOpenacc.plt', status='replace', action='write', form='formatted')
  write(seriesUnit,'(A)') 'TITLE = "3D OpenACC Nu/Re volume averages"'
  write(seriesUnit,'(A)') 'VARIABLES = "time" "NuVolAvg" "ReVolAvg"'
  write(seriesUnit,'(A)') 'ZONE T="NuReVolAvg", F=POINT'

  open(newunit=runningUnit, file='NuRe_VolAvg_runningMean_3DOpenacc.plt', status='replace', action='write', &
       form='formatted')
  write(runningUnit,'(A)') 'TITLE = "3D OpenACC Nu/Re running means"'
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

  open(unit=33, file='NuRe_TimeAverage_3DOpenacc.txt', status='replace', action='write', form='formatted')
  write(33,'(A)') '# 3D OpenACC Nu/Re statistical-convergence window averages'
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

end subroutine output_unsteady_NuRe_postprocess3d
!===========================================================================================================================
#endif


#ifdef steadyFlow
!===========================================================================================================================
! 子程序: RBcalc_Nu_global3d
! 作用: 计算三维算例的全局平均 Nusselt 数。
!===========================================================================================================================
subroutine RBcalc_Nu_global3d()
  use commondata3dOpenacc
  implicit none
  integer(kind=4) :: i, j, k
  real(kind=8) :: dz, dTdz, qz, sum_qz
  real(kind=8) :: deltaT, coef

  dz = 1.0d0 / lengthUnit
  deltaT = Thot - Tcold
  coef = velocityScaleCompare
  sum_qz = 0.0d0

 !$acc parallel loop collapse(3) default(none) present(w,T) reduction(+:sum_qz) private(dTdz,qz)
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

  Nu_global = (sum_qz / dble(nx * ny * nz)) / deltaT

  write(*,'(a,1x,ES24.16E3)') 'Nu_global =', real(Nu_global,kind=8)
  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,'(a,1x,ES24.16E3)') 'Nu_global =', real(Nu_global,kind=8)
  close(00)

  return
end subroutine RBcalc_Nu_global3d


!===========================================================================================================================
! 子程序: RBcalc_Nu_wall_avg3d
! 作用: 计算壁面平均 Nu、中心面 Nu 以及相关后处理量。
!===========================================================================================================================
subroutine RBcalc_Nu_wall_avg3d()
  use commondata3dOpenacc
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

 !$acc parallel loop default(none) present(T) copyout(T_bot_avg)
  do i = 1, nx
    T_bot_avg(i) = 0.0d0
    do j = 1, ny
      T_bot_avg(i) = T_bot_avg(i) + T(i,j,1)
    enddo
    T_bot_avg(i) = T_bot_avg(i) / dble(ny)
  enddo

  sum_hot = 0.0d0
 !$acc parallel loop default(none) present(T) copyout(Nu_bot) reduction(+:sum_hot) private(qz_wall)
  do i = 1, nx
    Nu_bot(i) = 0.0d0
    do j = 1, ny
      qz_wall = (8.0d0*Thot - 9.0d0*T(i,j,1) + T(i,j,2)) / (3.0d0*dz)
      Nu_bot(i) = Nu_bot(i) + qz_wall / deltaT
    enddo
    Nu_bot(i) = Nu_bot(i) / dble(ny)
    sum_hot = sum_hot + Nu_bot(i)
  enddo
  Nu_hot = sum_hot / dble(nx)

  Nu_bot_ext(1:nx) = Nu_bot(1:nx)
  xfit(1) = xp(1);  Tfit(1) = T_bot_avg(1)
  xfit(2) = xp(2);  Tfit(2) = T_bot_avg(2)
  xfit(3) = xp(3);  Tfit(3) = T_bot_avg(3)
  xfit(4) = xp(4);  Tfit(4) = T_bot_avg(4)
  call fit_adiabatic_wall_T4_3d(0.0d0, xfit, Tfit, T_wl)
  Nu_bot_ext(0) = (2.0d0 * (Thot - T_wl) / dz) / deltaT

  xfit(1) = xp(nx-3);  Tfit(1) = T_bot_avg(nx-3)
  xfit(2) = xp(nx-2);  Tfit(2) = T_bot_avg(nx-2)
  xfit(3) = xp(nx-1);  Tfit(3) = T_bot_avg(nx-1)
  xfit(4) = xp(nx  );  Tfit(4) = T_bot_avg(nx  )
  call fit_adiabatic_wall_T4_3d(xp(nx+1), xfit, Tfit, T_wr)
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
  call fit_parabola_ls5_3d(xk, fk, +1, fstar, xstar)
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
  call fit_parabola_ls5_3d(xk, fk, -1, fstar, xstar)
  Nu_hot_min = fstar
  Nu_hot_min_position = xstar

  sum_cold = 0.0d0
 !$acc parallel loop collapse(2) default(none) present(T) reduction(+:sum_cold) private(qz_wall)
  do j = 1, ny
    do i = 1, nx
      qz_wall = (-8.0d0*Tcold + 9.0d0*T(i,j,nz) - T(i,j,nz-1)) / (3.0d0*dz)
      sum_cold = sum_cold + qz_wall / deltaT
    enddo
  enddo
  Nu_cold = sum_cold / dble(nx * ny)

  sum_mid = 0.0d0
  if (mod(nz,2) .EQ. 1) then
    kMid = (nz + 1) / 2
 !$acc parallel loop collapse(2) default(none) present(w,T) reduction(+:sum_mid)
    do j = 1, ny
      do i = 1, nx
        sum_mid = sum_mid + (coef * w(i,j,kMid) * (T(i,j,kMid) - Tref) - &
             (T(i,j,kMid+1) - T(i,j,kMid-1)) / (2.0d0 * dz)) / deltaT
      enddo
    enddo
  else
    kB = nz / 2
    kT = kB + 1
 !$acc parallel loop collapse(2) default(none) present(w,T) reduction(+:sum_mid)
    do j = 1, ny
      do i = 1, nx
        sum_mid = sum_mid + (coef * 0.5d0 * (w(i,j,kB) * (T(i,j,kB) - Tref) + &
             w(i,j,kT) * (T(i,j,kT) - Tref)) + (T(i,j,kB) - T(i,j,kT)) / dz) / deltaT
      enddo
    enddo
  endif
  Nu_middle = sum_mid / dble(nx * ny)

  write(*,'(a,1x,ES24.16E3)') 'Nu_hot(bottom) =', real(Nu_hot,kind=8)
  write(*,'(a,1x,ES24.16E3)') 'Nu_cold(top)   =', real(Nu_cold,kind=8)
  write(*,'(a,1x,ES24.16E3)') 'Nu_middle      =', real(Nu_middle,kind=8)
  write(*,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
       'Nu_hot_max =', real(Nu_hot_max,kind=8), 'x_max =', real(Nu_hot_max_position,kind=8)
  write(*,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
       'Nu_hot_min =', real(Nu_hot_min,kind=8), 'x_min =', real(Nu_hot_min_position,kind=8)
  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,'(a,1x,ES24.16E3)') 'Nu_hot(bottom) =', real(Nu_hot,kind=8)
  write(00,'(a,1x,ES24.16E3)') 'Nu_cold(top)   =', real(Nu_cold,kind=8)
  write(00,'(a,1x,ES24.16E3)') 'Nu_middle      =', real(Nu_middle,kind=8)
  write(00,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
       'Nu_hot_max =', real(Nu_hot_max,kind=8), 'x_max =', real(Nu_hot_max_position,kind=8)
  write(00,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
       'Nu_hot_min =', real(Nu_hot_min,kind=8), 'x_min =', real(Nu_hot_min_position,kind=8)
  close(00)

  return
end subroutine RBcalc_Nu_wall_avg3d
#endif


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
  use commondata3dOpenacc
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
  use commondata3dOpenacc
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
  use commondata3dOpenacc
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
  use commondata3dOpenacc
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4) :: j, k, iL, iR
  real(kind=8) :: targetX, weight, valU, valV, valW, valT
  real(kind=8) :: uSlice(ny,nz), vSlice(ny,nz), wSlice(ny,nz), tSlice(ny,nz)

  ! midX 是 y-z 截面；浮力沿 z 方向作用，所以 W_nd 是竖直速度。
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
  use commondata3dOpenacc
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4) :: i, k, jL, jR
  real(kind=8) :: targetY, weight, valU, valV, valW, valT
  real(kind=8) :: uSlice(nx,nz), vSlice(nx,nz), wSlice(nx,nz), tSlice(nx,nz)

  ! midY 是 x-z 截面；浮力沿 z 方向作用，所以 W_nd 是竖直速度。
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
  use commondata3dOpenacc
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4) :: i, j, kL, kR
  real(kind=8) :: targetZ, weight, valU, valV, valW, valT
  real(kind=8) :: uSlice(nx,ny), vSlice(nx,ny), wSlice(nx,ny), tSlice(nx,ny)

  ! midZ 是 x-y 水平中截面。
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


#ifdef steadyFlow
!===========================================================================================================================
! 子程序: SideHeatedcalc_Nu_global3d
! 作用: 计算侧壁差温工况下的全场平均 Nusselt 数。
!===========================================================================================================================
subroutine SideHeatedcalc_Nu_global3d()
  use commondata3dOpenacc
  implicit none
  integer(kind=4) :: i, j, k
  real(kind=8) :: dy, dTdy, qy, sum_qy
  real(kind=8) :: deltaT, coef

  dy = 1.0d0 / lengthUnit
  deltaT = Thot - Tcold
  coef = velocityScaleCompare
  sum_qy = 0.0d0

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

  Nu_global = (sum_qy / dble(nx * ny * nz)) / deltaT

  write(*,'(a,1x,ES24.16E3)') 'Nu_global =', real(Nu_global,kind=8)
  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,'(a,1x,ES24.16E3)') 'Nu_global =', real(Nu_global,kind=8)
  close(00)

  return
end subroutine SideHeatedcalc_Nu_global3d


!===========================================================================================================================
! 子程序: SideHeatedcalc_Nu_wall_avg3d
! 作用: 计算侧壁差温工况下热壁、冷壁和中面的 Nusselt 数及其极值。
!===========================================================================================================================
subroutine SideHeatedcalc_Nu_wall_avg3d()
  use commondata3dOpenacc
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

  do i = 1, nx
    T_hot_avg(i) = 0.0d0
    do k = 1, nz
      T_hot_avg(i) = T_hot_avg(i) + T(i,1,k)
    enddo
    T_hot_avg(i) = T_hot_avg(i) / dble(nz)
  enddo

  sum_hot = 0.0d0
  do i = 1, nx
    Nu_hot_line(i) = 0.0d0
    do k = 1, nz
      qy_wall = (8.0d0*Thot - 9.0d0*T(i,1,k) + T(i,2,k)) / (3.0d0*dy)
      Nu_hot_line(i) = Nu_hot_line(i) + qy_wall / deltaT
    enddo
    Nu_hot_line(i) = Nu_hot_line(i) / dble(nz)
    sum_hot = sum_hot + Nu_hot_line(i)
  enddo
  Nu_hot = sum_hot / dble(nx)

  Nu_hot_ext(1:nx) = Nu_hot_line(1:nx)
  xfit(1) = xp(1);  Tfit(1) = T_hot_avg(1)
  xfit(2) = xp(2);  Tfit(2) = T_hot_avg(2)
  xfit(3) = xp(3);  Tfit(3) = T_hot_avg(3)
  xfit(4) = xp(4);  Tfit(4) = T_hot_avg(4)
  call fit_adiabatic_wall_T4_3d(0.0d0, xfit, Tfit, T_wf)
  Nu_hot_ext(0) = (2.0d0 * (Thot - T_wf) / dy) / deltaT

  xfit(1) = xp(nx-3);  Tfit(1) = T_hot_avg(nx-3)
  xfit(2) = xp(nx-2);  Tfit(2) = T_hot_avg(nx-2)
  xfit(3) = xp(nx-1);  Tfit(3) = T_hot_avg(nx-1)
  xfit(4) = xp(nx  );  Tfit(4) = T_hot_avg(nx  )
  call fit_adiabatic_wall_T4_3d(xp(nx+1), xfit, Tfit, T_wb)
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
  call fit_parabola_ls5_3d(xk, fk, +1, fstar, xstar)
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
  call fit_parabola_ls5_3d(xk, fk, -1, fstar, xstar)
  Nu_hot_min = fstar
  Nu_hot_min_position = xstar

  sum_cold = 0.0d0
  do k = 1, nz
    do i = 1, nx
      qy_wall = (-8.0d0*Tcold + 9.0d0*T(i,ny,k) - T(i,ny-1,k)) / (3.0d0*dy)
      sum_cold = sum_cold + qy_wall / deltaT
    enddo
  enddo
  Nu_cold = sum_cold / dble(nx * nz)

  sum_mid = 0.0d0
  if (mod(ny,2) .EQ. 1) then
    jMid = (ny + 1) / 2
    do k = 1, nz
      do i = 1, nx
        sum_mid = sum_mid + (coef * v(i,jMid,k) * (T(i,jMid,k) - Tref) + &
             (T(i,jMid-1,k) - T(i,jMid+1,k)) / (2.0d0*dy)) / deltaT
      enddo
    enddo
  else
    jL = ny / 2
    jR = jL + 1
    do k = 1, nz
      do i = 1, nx
        sum_mid = sum_mid + (coef * 0.5d0 * (v(i,jL,k) * (T(i,jL,k) - Tref) + &
             v(i,jR,k) * (T(i,jR,k) - Tref)) + (T(i,jL,k) - T(i,jR,k)) / dy) / deltaT
      enddo
    enddo
  endif
  Nu_middle = sum_mid / dble(nx * nz)

  write(*,'(a,1x,ES24.16E3)') 'Nu_hot(left)  =', real(Nu_hot,kind=8)
  write(*,'(a,1x,ES24.16E3)') 'Nu_cold(right)=', real(Nu_cold,kind=8)
  write(*,'(a,1x,ES24.16E3)') 'Nu_middle     =', real(Nu_middle,kind=8)
  write(*,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
       'Nu_hot_max =', real(Nu_hot_max,kind=8), 'x_max =', real(Nu_hot_max_position,kind=8)
  write(*,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
       'Nu_hot_min =', real(Nu_hot_min,kind=8), 'x_min =', real(Nu_hot_min_position,kind=8)

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,'(a,1x,ES24.16E3)') 'Nu_hot(left)  =', real(Nu_hot,kind=8)
  write(00,'(a,1x,ES24.16E3)') 'Nu_cold(right)=', real(Nu_cold,kind=8)
  write(00,'(a,1x,ES24.16E3)') 'Nu_middle     =', real(Nu_middle,kind=8)
  write(00,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
       'Nu_hot_max =', real(Nu_hot_max,kind=8), 'x_max =', real(Nu_hot_max_position,kind=8)
  write(00,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
       'Nu_hot_min =', real(Nu_hot_min,kind=8), 'x_min =', real(Nu_hot_min_position,kind=8)
  close(00)

  return
end subroutine SideHeatedcalc_Nu_wall_avg3d
#endif


!===========================================================================================================================
! 子程序: fit_adiabatic_wall_T4_3d
! 作用: 用四点拟合估计绝热壁面的壁温。
!===========================================================================================================================
subroutine fit_adiabatic_wall_T4_3d(y0, y, tt, T_wall)
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
end subroutine fit_adiabatic_wall_T4_3d


!===========================================================================================================================
! 子程序: fit_parabola_ls5_3d
! 作用: 用五点最小二乘抛物线拟合局部极值和对应位置。
!===========================================================================================================================
subroutine fit_parabola_ls5_3d(y, f, mode, fstar, ystar)
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
end subroutine fit_parabola_ls5_3d


!===========================================================================================================================
! 子程序: SideHeatedcalc_umid_max3d
! 作用: 侧壁差温工况的 u 中面最大值诊断入口。
!===========================================================================================================================
subroutine SideHeatedcalc_umid_max3d()
  call calc_umid_max_common3d('SideHeatedcalc_umid_max')
end subroutine SideHeatedcalc_umid_max3d

!===========================================================================================================================
! 子程序: SideHeatedcalc_vmid_max3d
! 作用: 侧壁差温工况的 v 中面最大值诊断入口。
!===========================================================================================================================
subroutine SideHeatedcalc_vmid_max3d()
  call calc_vmid_max_common3d('SideHeatedcalc_vmid_max')
end subroutine SideHeatedcalc_vmid_max3d

!===========================================================================================================================
! 子程序: SideHeatedcalc_wmid_max3d
! 作用: 侧壁差温工况的 w 中面最大值诊断入口。
!===========================================================================================================================
subroutine SideHeatedcalc_wmid_max3d()
  call calc_wmid_max_common3d('SideHeatedcalc_wmid_max')
end subroutine SideHeatedcalc_wmid_max3d

!===========================================================================================================================
! 子程序: RBcalc_umid_max3d
! 作用: Rayleigh-Benard 工况的 u 中面最大值诊断入口。
!===========================================================================================================================
subroutine RBcalc_umid_max3d()
  call calc_umid_max_common3d('RBcalc_umid_max')
end subroutine RBcalc_umid_max3d

!===========================================================================================================================
! 子程序: RBcalc_vmid_max3d
! 作用: Rayleigh-Benard 工况的 v 中面最大值诊断入口。
!===========================================================================================================================
subroutine RBcalc_vmid_max3d()
  call calc_vmid_max_common3d('RBcalc_vmid_max')
end subroutine RBcalc_vmid_max3d

!===========================================================================================================================
! 子程序: RBcalc_wmid_max3d
! 作用: Rayleigh-Benard 工况的 w 中面最大值诊断入口。
!===========================================================================================================================
subroutine RBcalc_wmid_max3d()
  call calc_wmid_max_common3d('RBcalc_wmid_max')
end subroutine RBcalc_wmid_max3d


!===========================================================================================================================
! 子程序: calc_umid_max_common3d
! 作用: 在 x 中面插值搜索 u 的最大值，并输出对应的 y/z 位置。
!===========================================================================================================================
subroutine calc_umid_max_common3d(logTag)
  use commondata3dOpenacc
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
       'u_mid_max =', real(umax*velocityScaleCompare,kind=8), 'at y =', real(yAtU,kind=8), &
       'z =', real(zAtU,kind=8), 'on x_mid =', real(targetX,kind=8)

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) '--- ', trim(logTag), ' ---'
  write(00,*) 'x_mid =', targetX
  write(00,*) 'u_mid_max =', umax*velocityScaleCompare, ' y_pos =', yAtU, ' z_pos =', zAtU
  close(00)

end subroutine calc_umid_max_common3d


!===========================================================================================================================
! 子程序: calc_vmid_max_common3d
! 作用: 在 y 中面插值搜索 v 的最大值，并输出对应的 x/z 位置。
!===========================================================================================================================
subroutine calc_vmid_max_common3d(logTag)
  use commondata3dOpenacc
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
       'v_mid_max =', real(vmax*velocityScaleCompare,kind=8), 'at x =', real(xAtV,kind=8), &
       'z =', real(zAtV,kind=8), 'on y_mid =', real(targetY,kind=8)

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) '--- ', trim(logTag), ' ---'
  write(00,*) 'y_mid =', targetY
  write(00,*) 'v_mid_max =', vmax*velocityScaleCompare, ' x_pos =', xAtV, ' z_pos =', zAtV
  close(00)

end subroutine calc_vmid_max_common3d


!===========================================================================================================================
! 子程序: calc_wmid_max_common3d
! 作用: 在 z 中面插值搜索 w 的最大值，并输出对应的 x/y 位置。
!===========================================================================================================================
subroutine calc_wmid_max_common3d(logTag)
  use commondata3dOpenacc
  implicit none
  character(len=*), intent(in) :: logTag
  integer(kind=4) :: i, j, kL, kR, iBest, jBest
  real(kind=8) :: targetZ, weight, val, wmax, xAtW, yAtW

  ! 按当前 z 方向约定，W 是浮力方向/竖直方向速度分量。
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
       'w_mid_max =', real(wmax*velocityScaleCompare,kind=8), 'at x =', real(xAtW,kind=8), &
       'y =', real(yAtW,kind=8), 'on z_mid =', real(targetZ,kind=8)

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) '--- ', trim(logTag), ' ---'
  write(00,*) 'z_mid =', targetZ
  write(00,*) 'w_mid_max =', wmax*velocityScaleCompare, ' x_pos =', xAtW, ' y_pos =', yAtW
  close(00)

end subroutine calc_wmid_max_common3d


!===========================================================================================================================
! 子程序: write_midplane_stream_x
! 作用: 输出 x=Lx/2 中面的流函数/涡量诊断切片。
!===========================================================================================================================
subroutine write_midplane_stream_x(filename)
  use commondata3dOpenacc
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4) :: j, k, m, iL, iR
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
  use commondata3dOpenacc
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4) :: i, k, m, jL, jR
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
  use commondata3dOpenacc
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4) :: i, j, m, kL, kR
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
  use commondata3dOpenacc
  implicit none

  character(len=*), intent(in) :: filename
  integer(kind=4) :: i, j, k, uout
  real(kind=4) :: zoneMarker, eohMarker
  character(len=40) :: zoneName

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
        write(uout) real(xp(i),kind=8)
        write(uout) real(yp(j),kind=8)
        write(uout) real(zp(k),kind=8)
        write(uout) real(u(i,j,k),kind=8)
        write(uout) real(v(i,j,k),kind=8)
        write(uout) real(w(i,j,k),kind=8)
        write(uout) real(T(i,j,k),kind=8)
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
  real(kind=4) :: zoneMarker, eohMarker
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
      write(uout) real(coord1(i),kind=8)
      write(uout) real(coord2(j),kind=8)
      write(uout) real(field1(i,j),kind=8)
      write(uout) real(field2(i,j),kind=8)
      write(uout) real(field3(i,j),kind=8)
      write(uout) real(field4(i,j),kind=8)
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
  real(kind=4) :: zoneMarker, eohMarker
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
  write(uout) 2
  write(uout) 2
  write(uout) 2
  write(uout) 2

  write(uout) 0
  write(uout) -1

  do j = 1, nj
    do i = 1, ni
      write(uout) real(coord1(i),kind=8)
      write(uout) real(coord2(j),kind=8)
      write(uout) real(psi(i,j),kind=8)
      write(uout) real(vort(i,j),kind=8)
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
