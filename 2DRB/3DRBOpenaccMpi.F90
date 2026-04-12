!=============================================================
!!!    注释区，代码描述
!!!    三维浮力驱动自然对流
!!!    D3Q19 流场 + D3Q7 温度场
!!!    主循环按 collision/streaming/boundary/macro 展开
!!!    参考 3DRB.F90 的主流程与后处理结构，并保留 MPI/OpenACC 必需部分
!=============================================================

!=============================================================
!   自定义宏，一些选项的开关
!steadyFlow / unsteadyFlow : 稳态或非稳态运行模式
#define steadyFlow
!#define unsteadyFlow

!速度边界，包括水平垂直展向边界无滑移，还有垂直展向边界速度周期
!spanwise 表示 z 方向前后展向壁面
#define HorizontalWallsNoslip
#define VerticalWallsNoslip
#define SpanwiseWallsNoslip
!#define VerticalWallsPeriodicalU
!#define SpanwiseWallsPeriodicalU

!温度边界(for Rayleigh Benard Cell)，包括水平边界恒温，垂直/展向边界温度不可穿透以及周期
!#define RayleighBenardCell
!#define HorizontalWallsConstT
!#define VerticalWallsAdiabatic
!#define VerticalWallsPeriodicalT
!#define SpanwiseWallsAdiabatic
!#define SpanwiseWallsPeriodicalT

!温度边界(for Side Heated Cell)，包括水平/展向边界温度不可穿透，垂直边界恒温
#define SideHeatedCell
#define VerticalWallsConstT
#define HorizontalWallsAdiabatic
#define SpanwiseWallsAdiabatic

!算法切换
! RB 小 Ra 初始温度扰动：只对 RayleighBenardCell 且 Rayleigh<=1e4 生效
#define EnableRBInitPerturbation3D
#define EnableUseG
!#define EnableLegacyThermalScheme


!   自定义宏结束
!=============================================================


!=============================================================
!   全局模块
module commondata3dOpenaccMpi
  implicit none

  ! qf / qt 分别表示流场和温度场的离散速度数
  integer(kind=4), parameter :: qf = 19
  integer(kind=4), parameter :: qt = 7

  ! 重启控制参数：
  ! loadInitField=0 表示从初值启动
  ! loadInitField=1 表示从当前 rank 对应的 backupFile3D-rankXXXX-*.bin 重启
  integer(kind=4), parameter :: loadInitField = 0
  integer(kind=4), parameter :: reloadDimensionlessTime = 0
  integer(kind=4), parameter :: reloadbinFileNum = 0


  integer(kind=4), parameter :: nx = 120, ny = 120, nzGlobal = 120
#ifdef SideHeatedCell
  real(kind=8), parameter :: lengthUnit = dble(nx)
#else
  real(kind=8), parameter :: lengthUnit = dble(ny)
#endif
  real(kind=8), parameter :: pi = acos(-1.0d0)
  integer(kind=4), parameter :: itc_max = 20000000
  real(kind=8),    parameter :: outputFrequency = 100.0d0
  ! 主要控制参数：Rayleigh / Prandtl / Mach / 热壁温度
  real(kind=8), parameter :: Rayleigh = 1.0d7
  real(kind=8), parameter :: Prandtl = 0.71d0
  real(kind=8), parameter :: Mach = 0.1d0
  real(kind=8), parameter :: Thot = 0.5d0
  real(kind=8), parameter :: Tcold = -0.5d0
  real(kind=8), parameter :: Tref = 0.5d0 * (Thot + Tcold)

  ! 松弛时间和输运系数
  real(kind=8), parameter :: tauf = 0.5d0 + Mach * lengthUnit * dsqrt(3.0d0 * Prandtl / Rayleigh)
  real(kind=8), parameter :: viscosity = (tauf - 0.5d0) / 3.0d0
  real(kind=8), parameter :: diffusivity = viscosity / Prandtl
  real(kind=8), parameter :: cs2T = 0.25d0
  real(kind=8), parameter :: paraA = 42.0d0 * dsqrt(3.0d0) * diffusivity - 6.0d0

  ! 动量方程的多松弛系数
  real(kind=8), parameter :: Se = 1.0d0 / tauf, Seps = 1.0d0 / tauf
  real(kind=8), parameter :: Snu = 1.0d0 / tauf, Spi = 1.0d0 / tauf
  real(kind=8), parameter :: Sq = 8.0d0 * (2.0d0 * tauf - 1.0d0) / (8.0d0 * tauf - 1.0d0)
  real(kind=8), parameter :: Sm = 8.0d0 * (2.0d0 * tauf - 1.0d0) / (8.0d0 * tauf - 1.0d0)

  ! 温度方程的多松弛系数
#ifdef EnableLegacyThermalScheme
  real(kind=8), parameter :: Qk = 3.0d0 - dsqrt(3.0d0), Qnu = 4.0d0 * dsqrt(3.0d0) - 6.0d0
  real(kind=8), parameter :: thermalGeqCoeff = 7.0d0 / ((6.0d0 + paraA) * cs2T)
#else
  real(kind=8), parameter :: taug = 0.5d0 + diffusivity / cs2T
  real(kind=8), parameter :: Qnu = 1.0d0, Qk = 1.0d0 / taug
  real(kind=8), parameter :: thermalGeqCoeff = 1.0d0 / cs2T
#endif


  ! heatFluxScale         : Nu 相关量的系数，等于 L/kappa
  ! velocityScaleCompare  : 输出速度时使用的无量纲速度标度
  real(kind=8), parameter :: heatFluxScale = lengthUnit / diffusivity
  real(kind=8), parameter :: velocityScaleCompare = lengthUnit / diffusivity

  ! 浮力项系数：
  ! gBeta 对应 Fy = rho * gBeta * (T - Tref)
  real(kind=8), parameter :: gBeta1 = Rayleigh * viscosity * diffusivity / lengthUnit
  real(kind=8), parameter :: gBeta = gBeta1 / lengthUnit / lengthUnit
  real(kind=8), parameter :: timeUnit = dsqrt(lengthUnit / gBeta)
  real(kind=8), parameter :: velocityUnit = dsqrt(gBeta * lengthUnit)

  ! 输出相关控制量
  integer(kind=4), parameter :: dimensionlessTimeMax = max(1, int(12000.0d0 / outputFrequency))
  integer(kind=4), parameter :: backupInterval = 1000      ! 备份间隔

  ! 稳态判据
  real(kind=8), parameter :: epsU = 1.0d-7
  real(kind=8), parameter :: epsT = 1.0d-7

  integer(kind=4), parameter :: outputBinFile = 0
  integer(kind=4), parameter :: outputPltFile = 0

  integer(kind=4) :: binFileNum, pltFileNum
  integer(kind=4) :: dimensionlessTime
  integer(kind=4) :: outputIntervalItc, backupIntervalItc
  integer(kind=4), parameter :: mpiRoot = 0
  integer(kind=4) :: nz, mpiRank, mpiSize, mpiErr, mpiLeft, mpiRight
  integer(kind=4) :: zStartGlobal, zEndGlobal, accDeviceId, numAccDevices
  logical :: isRoot, isFirstZRank, isLastZRank
  character(len=16) :: rankTag

  ! 体平均 Nu / Re 的时间序列缓存
  real(kind=8) :: NuVolAvg(0:dimensionlessTimeMax), ReVolAvg(0:dimensionlessTimeMax)

  ! 输出文件命名保留 MPI rank 后缀，但基名尽量与 3DRB 对齐
  character(len=100) :: binFilePrefix = "buoyancyCavity3DbinFile"
  character(len=100) :: pltFilePrefix = "buoyancyCavity3DTecplot"
  character(len=100) :: reloadFilePrefix = "backupFile3D"
  character(len=100) :: settingsFile = "SimulationSettings3D.txt"

  real(kind=8) :: errorU, errorT

  ! 几何坐标数组：包含物理边界点 0 和 nx+1/ny+1/nz+1
  real(kind=8) :: xp(0:nx+1), yp(0:ny+1)
  real(kind=8), allocatable :: zp(:)
  ! 宏观场：u,v,w,T,rho
  real(kind=8), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:), T(:,:,:), rho(:,:,:)

#ifdef steadyFlow
  ! 稳态误差判据需要保存上一次输出时刻的场
  real(kind=8), allocatable :: up(:,:,:), vp(:,:,:), wp(:,:,:), Tp(:,:,:)
#endif

  ! f/g 是当前分布函数，f_post/g_post 是碰撞后、迁移前的分布函数，post 数组带 ghost 层便于直接 pull streaming
  real(kind=8), allocatable :: f(:,:,:,:), f_post(:,:,:,:)
  real(kind=8), allocatable :: g(:,:,:,:), g_post(:,:,:,:)
  ! 温度方程中的历史热流项
  real(kind=8), allocatable :: Bx_prev(:,:,:), By_prev(:,:,:), Bz_prev(:,:,:)
  real(kind=8), allocatable :: fSendLower(:), fSendUpper(:), fRecvLower(:), fRecvUpper(:)
  real(kind=8), allocatable :: gSendLower(:), gSendUpper(:), gRecvLower(:), gRecvUpper(:)

  integer(kind=4) :: itc

#ifdef EnableUseG
  logical, parameter :: useG = .true.
#else
  logical, parameter :: useG = .false.
#endif

#ifdef EnableLegacyThermalScheme
  logical, parameter :: useLegacyThermalScheme = .true.
#else
  logical, parameter :: useLegacyThermalScheme = .false.
#endif

  real(kind=8) :: Nu_global, Nu_hot, Nu_cold, Nu_middle
  real(kind=8) :: Nu_hot_max, Nu_hot_min, Nu_hot_max_position, Nu_hot_min_position

  ! D3Q19 / D3Q7 的离散速度、反向索引和权重
  integer(kind=4) :: ex(0:qf-1), ey(0:qf-1), ez(0:qf-1), opp(0:qf-1)
  real(kind=8)    :: omega(0:qf-1)

  integer(kind=4) :: exT(0:qt-1), eyT(0:qt-1), ezT(0:qt-1), oppT(0:qt-1)
  real(kind=8)    :: omegaT(0:qt-1)

end module commondata3dOpenaccMpi


!   全局模块结束
!=============================================================


!=============================================================
!   主程序


program main3dOpenaccMpi
  use mpi
  use openacc
  use commondata3dOpenaccMpi
  implicit none

  real(kind=8) :: timeStart, timeEnd
  real(kind=8) :: timeStart2, timeEnd2
  character(len=24) :: ctime
  character(len=24) :: string
  integer(kind=4) :: time
  integer(kind=8) :: wallClockStart, wallClockEnd, wallClockRate

  call MPI_Init(mpiErr)
  call MPI_Comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)
  call MPI_Comm_size(MPI_COMM_WORLD, mpiSize, mpiErr)
  isRoot = (mpiRank .EQ. mpiRoot)

  call setup_mpi_decomposition_3d()
  call compose_ranked_file_prefixes_3d()

  ! 每个 MPI rank 独立记录自己的日志和输出前缀，避免并发写同一个文件。
  open(unit=00, file=trim(settingsFile), status='replace')
  string = ctime(time())
  write(00,*) 'Start: ', string
  write(00,*) 'Starting OpenACC+MPI >>>>>>', ' rank =', mpiRank, '; size =', mpiSize
  write(00,*) 'Local nz =', nz, '; Global z range =', zStartGlobal, zEndGlobal, '; nzGlobal =', nzGlobal

  call acc_init(acc_device_default)
  numAccDevices = acc_get_num_devices(acc_device_default)
  if (numAccDevices .GT. 0) then
    accDeviceId = mod(mpiRank, numAccDevices)
    call acc_set_device_num(accDeviceId, acc_device_default)
  else
    accDeviceId = -1
  endif
  write(00,*) 'Visible OpenACC devices:', numAccDevices
  write(00,*) 'Assigned OpenACC device:', accDeviceId
  close(00)

  call initial3d()
  call enter_data_3d_openacc()

  call CPU_TIME(timeStart)
  call system_clock(wallClockStart, wallClockRate)
  timeStart2 = dble(wallClockStart) / dble(max(wallClockRate,1_8))

  do while (((errorU .GT. epsU) .OR. (errorT .GT. epsT)) .AND. (itc .LE. itc_max))
    itc = itc + 1

    call collision3d()
    call fill_periodic_ghosts_f_post()
    call streaming3d()
    call bounceback3d()
    call macro3d()

    call collisionT3d()
    call fill_periodic_ghosts_g_post()
    call streamingT3d()
    call bouncebackT3d()
    call macroT3d()

#ifdef steadyFlow
    if (mod(itc, 2000) .EQ. 0) call check3d()
#endif

    if (mod(itc, outputIntervalItc) .EQ. 0) then
      !call calNuRe3d()

#ifdef steadyFlow
      if ((outputPltFile .EQ. 1) .AND. (mod(itc, backupIntervalItc) .EQ. 0)) then
        call update_host_all_3d_openacc()
        call output_Tecplot3d()
      endif
#endif

#ifdef unsteadyFlow
      if (outputBinFile .EQ. 1) then
        call update_host_all_3d_openacc()
        call output_binary3d()
        if (mod(itc, backupIntervalItc) .EQ. 0) call backupData3d()
      endif
      if (outputPltFile .EQ. 1) then
        call update_host_all_3d_openacc()
        call output_Tecplot3d()
      endif
#endif
    endif
  enddo

  call CPU_TIME(timeEnd)
  call system_clock(wallClockEnd, wallClockRate)
  timeEnd2 = dble(wallClockEnd) / dble(max(wallClockRate,1_8))

  call update_host_all_3d_openacc()

#ifdef steadyFlow
  if (outputPltFile .EQ. 1) call output_Tecplot3d()
  if (outputBinFile .EQ. 1) call output_binary3d()
#endif

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
  call output_Tecplot3d()

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) '======================================================================'
  write(00,*) 'Time (CPU) = ', real(timeEnd - timeStart, kind=8), 's'
  write(00,*) 'MLUPS = ', &
       real(dble(nx) * dble(ny) * dble(nzGlobal) * dble(itc) / &
       & max(timeEnd - timeStart, 1.0d-12) / 1.0d6, kind=8)
  write(00,*) 'Time (Wall) = ', real(timeEnd2 - timeStart2, kind=8), 's'
  write(00,*) 'MLUPS (Wall) = ', &
       real(dble(nx) * dble(ny) * dble(nzGlobal) * dble(itc) / &
       & max(timeEnd2 - timeStart2, 1.0d-12) / 1.0d6, kind=8)
  write(00,*) 'Nu_global =', Nu_global
  write(00,*) 'Nu_hot    =', Nu_hot
  write(00,*) 'Nu_cold   =', Nu_cold
  write(00,*) 'Nu_middle =', Nu_middle
  write(00,*) 'useG =', useG
  write(00,*) 'useLegacyThermalScheme =', useLegacyThermalScheme
  if (outputBinFile .EQ. 1) then
    call backupData3d()
  endif
  write(00,*) 'Deallocate Array......'
  close(00)

  call exit_data_3d_openacc()

  deallocate(f, f_post, g, g_post)
  deallocate(u, v, w, T, rho)
  deallocate(zp)
#ifdef steadyFlow
  deallocate(up, vp, wp, Tp)
#endif
  deallocate(Bx_prev, By_prev, Bz_prev)
  deallocate(fSendLower, fSendUpper, fRecvLower, fRecvUpper)
  deallocate(gSendLower, gSendUpper, gRecvLower, gRecvUpper)

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) 'Successfully: DNS completed!'
  string = ctime(time())
  write(00,*) 'End:   ', string
  close(00)

  call acc_shutdown(acc_device_default)
  call MPI_Finalize(mpiErr)

end program main3dOpenaccMpi

!   主程序结束
!=============================================================


!===========================================================================================================================
! 子程序: setup_mpi_decomposition_3d
! 作用: 采用 z 方向一维 slab 分解，为 MPI 多 GPU 计算准备本地子域信息。
!===========================================================================================================================
subroutine setup_mpi_decomposition_3d()
  use mpi
  use commondata3dOpenaccMpi
  implicit none
  integer(kind=4) :: baseNz, remainderNz

  if (mpiSize .GT. nzGlobal) then
    if (mpiRank .EQ. mpiRoot) then
      write(*,*) 'Error: mpiSize must not be larger than nzGlobal in 1D z decomposition.'
    endif
    call MPI_Abort(MPI_COMM_WORLD, 1, mpiErr)
  endif

  baseNz = nzGlobal / mpiSize
  remainderNz = mod(nzGlobal, mpiSize)

  if (mpiRank .LT. remainderNz) then
    nz = baseNz + 1
    zStartGlobal = mpiRank * (baseNz + 1) + 1
  else
    nz = baseNz
    zStartGlobal = remainderNz * (baseNz + 1) + (mpiRank - remainderNz) * baseNz + 1
  endif
  zEndGlobal = zStartGlobal + nz - 1

  if (nz .LE. 0) then
    write(*,*) 'Error: local nz <= 0 on rank ', mpiRank
    call MPI_Abort(MPI_COMM_WORLD, 2, mpiErr)
  endif

  isFirstZRank = (zStartGlobal .EQ. 1)
  isLastZRank = (zEndGlobal .EQ. nzGlobal)

#ifdef SpanwiseWallsPeriodicalU
  mpiLeft = mod(mpiRank - 1 + mpiSize, mpiSize)
  mpiRight = mod(mpiRank + 1, mpiSize)
#else
  if (isFirstZRank) then
    mpiLeft = MPI_PROC_NULL
  else
    mpiLeft = mpiRank - 1
  endif
  if (isLastZRank) then
    mpiRight = MPI_PROC_NULL
  else
    mpiRight = mpiRank + 1
  endif
#endif
end subroutine setup_mpi_decomposition_3d


!===========================================================================================================================
! 子程序: compose_ranked_file_prefixes_3d
! 作用: 为每个 MPI rank 生成独立的日志、重启和输出文件名前缀。
!===========================================================================================================================
subroutine compose_ranked_file_prefixes_3d()
  use commondata3dOpenaccMpi
  implicit none
  write(rankTag,'("rank",I4.4)') mpiRank
  binFilePrefix = 'buoyancyCavity3DbinFile-' // trim(rankTag)
  pltFilePrefix = 'buoyancyCavity3DTecplot-' // trim(rankTag)
  reloadFilePrefix = 'backupFile3D-' // trim(rankTag)
  settingsFile = 'SimulationSettings3D-' // trim(rankTag) // '.txt'
end subroutine compose_ranked_file_prefixes_3d




!===========================================================================================================================
! 子程序: initial3d
! 作用: 初始化网格坐标、场变量、分布函数、输出文件和重启信息。
!===========================================================================================================================
subroutine initial3d()
  use commondata3dOpenaccMpi
  implicit none

  integer(kind=4) :: i, j, k, alpha
  real(kind=8) :: eu, u2Loc
  real(kind=8) :: xLen, yLen, rbInitPerturbAmp
  character(len=100) :: reloadFileName

  itc = 0
  errorU = 100.0d0
  errorT = 100.0d0

  outputIntervalItc = max(1, int(outputFrequency * timeUnit))
  backupIntervalItc = max(1, int(backupInterval * timeUnit))

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  if (outputBinFile .EQ. 1) write(00,*) 'Binary snapshots will be stored with prefix: ', trim(binFilePrefix)
  if (outputPltFile .EQ. 1) write(00,*) 'Tecplot slices will be stored with prefix: ', trim(pltFilePrefix)

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
  write(00,*) 'Mesh local :', nx, ny, nz
  write(00,*) 'Mesh global:', nx, ny, nzGlobal
  write(00,*) 'MPI rank/size =', mpiRank, mpiSize, '; global z range =', zStartGlobal, zEndGlobal
  write(00,*) 'Rayleigh=', real(Rayleigh,kind=8), '; Prandtl =', real(Prandtl,kind=8), '; Mach =', real(Mach,kind=8)
  write(00,*) 'Length unit: L0 =', real(lengthUnit,kind=8)
  write(00,*) 'Time unit: Sqrt(L0/(gBeta*DeltaT)) =', real(timeUnit,kind=8)
  write(00,*) 'Velocity unit: Sqrt(gBeta*L0*DeltaT) =', real(velocityUnit,kind=8)
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
  write(00,*) 'outputFrequency =', real(outputFrequency,kind=8), ' free-fall time units'
  write(00,*) '......................  or ', outputIntervalItc, ' in itc units'
  write(00,*) 'backupInterval =', backupInterval, ' free-fall time units'
  write(00,*) '.................... or ', backupIntervalItc, ' in itc units'
  if (loadInitField .EQ. 1) write(00,*) 'reloadDimensionlessTime =', reloadDimensionlessTime
  write(00,*) 'itc_max =', itc_max
  write(00,*) 'default epsU =', real(epsU,kind=8), '; epsT =', real(epsT,kind=8)
  write(00,*) 'useG =', useG
  write(00,*) 'useLegacyThermalScheme =', useLegacyThermalScheme

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
  write(00,*) 'OpenACC+MPI version; z-slab domain decomposition in MPI, kernels on OpenACC devices'

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

  allocate(zp(0:nz+1))
  zp(0) = dble(zStartGlobal - 1)
  zp(nz+1) = dble(zEndGlobal)
  do k = 1, nz
    zp(k) = dble(zStartGlobal + k - 1) - 0.5d0
  enddo
  zp = zp / lengthUnit

  allocate(u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), T(nx,ny,nz), rho(nx,ny,nz))
#ifdef steadyFlow
  allocate(up(nx,ny,nz), vp(nx,ny,nz), wp(nx,ny,nz), Tp(nx,ny,nz))
#endif
  allocate(f(0:qf-1,nx,ny,nz), f_post(0:qf-1,0:nx+1,0:ny+1,0:nz+1))
  allocate(g(0:qt-1,nx,ny,nz), g_post(0:qt-1,0:nx+1,0:ny+1,0:nz+1))
  allocate(Bx_prev(nx,ny,nz), By_prev(nx,ny,nz), Bz_prev(nx,ny,nz))
  allocate(fSendLower(qf*(nx+2)*(ny+2)), fSendUpper(qf*(nx+2)*(ny+2)))
  allocate(fRecvLower(qf*(nx+2)*(ny+2)), fRecvUpper(qf*(nx+2)*(ny+2)))
  allocate(gSendLower(qt*(nx+2)*(ny+2)), gSendUpper(qt*(nx+2)*(ny+2)))
  allocate(gRecvLower(qt*(nx+2)*(ny+2)), gRecvUpper(qt*(nx+2)*(ny+2)))

  call init_lattice_constants_3d()

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
    if (reloadDimensionlessTime .NE. 0) then
      write(00,*) 'Error: since loadInitField .EQ. 0, reloadDimensionlessTime should also be 0'
      close(00)
      stop
    endif

#ifdef VerticalWallsNoslip
    write(00,*) 'Velocity B.C. for vertical walls are: ===No-slip wall==='
#endif
#ifdef VerticalWallsPeriodicalU
    write(00,*) 'Velocity B.C. for vertical walls are: ===Periodical==='
#endif
#ifdef HorizontalWallsNoslip
    write(00,*) 'Velocity B.C. for horizontal walls are: ===No-slip wall==='
#endif
#ifdef SpanwiseWallsNoslip
    write(00,*) 'Velocity B.C. for spanwise walls are: ===No-slip wall==='
#endif
#ifdef SpanwiseWallsPeriodicalU
    write(00,*) 'Velocity B.C. for spanwise walls are: ===Periodical==='
#endif

    u = 0.0d0
    v = 0.0d0
    w = 0.0d0
    T = 0.0d0

#ifdef VerticalWallsConstT
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          T(i,j,k) = Thot + (xp(i) - xp(0)) / (xp(nx+1) - xp(0)) * (Tcold - Thot)
        enddo
      enddo
    enddo
    write(00,*) 'Temperature B.C. for vertical walls are: ===Hot/cold wall==='
#endif

#ifdef HorizontalWallsConstT
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          T(i,j,k) = Thot + (yp(j) - yp(0)) / (yp(ny+1) - yp(0)) * (Tcold - Thot)
        enddo
      enddo
    enddo
    write(00,*) 'Temperature B.C. for horizontal walls are: ===Hot/cold wall==='
#ifdef RayleighBenardCell
#ifdef EnableRBInitPerturbation3D
    if (Rayleigh .LE. 1.0d4) then
      xLen = xp(nx+1)
      yLen = yp(ny+1)
      rbInitPerturbAmp = 1.0d-3 * (Thot - Tcold)
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            T(i,j,k) = T(i,j,k) + rbInitPerturbAmp * dcos(2.0d0 * pi * xp(i) / xLen) * dsin(pi * yp(j) / yLen)
          enddo
        enddo
      enddo
      write(00,'(a,1x,es12.4)') '3D RB initial T perturbation amplitude =', rbInitPerturbAmp
    else
      write(00,*) '3D RB initial T perturbation skipped because Rayleigh > 1.0d4'
    endif
#endif
#endif
#endif

#ifdef VerticalWallsAdiabatic
    write(00,*) 'Temperature B.C. for vertical walls are: ===Adiabatic wall==='
#endif
#ifdef VerticalWallsPeriodicalT
    write(00,*) 'Temperature B.C. for vertical walls are: ===Periodical==='
#endif
#ifdef HorizontalWallsAdiabatic
    write(00,*) 'Temperature B.C. for horizontal walls are: ===Adiabatic wall==='
#endif
#ifdef SpanwiseWallsAdiabatic
    write(00,*) 'Temperature B.C. for spanwise walls are: ===Adiabatic wall==='
#endif
#ifdef SpanwiseWallsPeriodicalT
    write(00,*) 'Temperature B.C. for spanwise walls are: ===Periodical==='
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
    if (reloadDimensionlessTime .EQ. 0) then
      write(00,*) 'WARNING: since loadInitField .EQ. 1, please confirm reloadDimensionlessTime', reloadDimensionlessTime
      close(00)
      stop
    endif
    write(00,*) 'Load initial field from previous simulation: ', trim(reloadFilePrefix), '- >>>'
    write(reloadFileName,'(i0)') reloadbinFileNum
    open(unit=01, file=trim(reloadFilePrefix)//'-'//trim(adjustl(reloadFileName))//'.bin', &
         form='unformatted', access='sequential', status='old')
    write(00,*) 'Reloading f and g from file'
    read(01) ((((f(alpha,i,j,k), i=1,nx), j=1,ny), k=1,nz), alpha=0,qf-1)
    read(01) ((((g(alpha,i,j,k), i=1,nx), j=1,ny), k=1,nz), alpha=0,qt-1)
    close(01)
    call reconstruct_macro_from_fg3d()
    write(00,*) 'Raw data is loaded from the file: ', trim(reloadFilePrefix), '-', trim(adjustl(reloadFileName)), '.bin'
  else
    write(00,*) 'Error: initial field is not properly set'
    close(00)
    stop
  endif

  write(00,*) '-------------------------------------------------------------------------------'
  close(00)

#ifdef steadyFlow
  up = 0.0d0
  vp = 0.0d0
  wp = 0.0d0
  Tp = 0.0d0
#endif

  f_post = 0.0d0
  g_post = 0.0d0
  binFileNum = 0
  pltFileNum = 0
  dimensionlessTime = 0
  NuVolAvg = 0.0d0
  ReVolAvg = 0.0d0

end subroutine initial3d
!===========================================================================================================================
! 子程序: enter_data_3d_openacc
! 作用: 在主时间推进前把主要数组和常量映射到 OpenACC 设备端。
!===========================================================================================================================
subroutine enter_data_3d_openacc()
  use openacc
  use commondata3dOpenaccMpi
  implicit none

  !$acc enter data copyin(xp,yp,zp,ex,ey,ez,opp,omega,exT,eyT,ezT,oppT,omegaT)
  !$acc enter data copyin(u,v,w,T,rho,f,g,Bx_prev,By_prev,Bz_prev)
  !$acc enter data create(f_post,g_post)
#ifdef steadyFlow
  !$acc enter data copyin(up,vp,wp,Tp)
#endif
end subroutine enter_data_3d_openacc


!===========================================================================================================================
! 子程序: update_host_all_3d_openacc
! 作用: 在需要输出或写文件前，把设备端场变量同步回主机端。
!===========================================================================================================================
subroutine update_host_all_3d_openacc()
  use commondata3dOpenaccMpi
  implicit none

  !$acc update self(u,v,w,T,rho,f,g,Bx_prev,By_prev,Bz_prev)
#ifdef steadyFlow
  !$acc update self(up,vp,wp,Tp)
#endif
end subroutine update_host_all_3d_openacc


!===========================================================================================================================
! 子程序: exit_data_3d_openacc
! 作用: 在计算结束后释放设备端数据映射。
!===========================================================================================================================
subroutine exit_data_3d_openacc()
  use commondata3dOpenaccMpi
  implicit none

#ifdef steadyFlow
  !$acc exit data delete(up,vp,wp,Tp)
#endif
  !$acc exit data delete(f_post,g_post,u,v,w,T,rho,f,g,Bx_prev,By_prev,Bz_prev)
  !$acc exit data delete(xp,yp,zp,ex,ey,ez,opp,omega,exT,eyT,ezT,oppT,omegaT)
end subroutine exit_data_3d_openacc

!===========================================================================================================================
! 子程序: exchange_f_post_halo_mpi
! 作用: 交换流场分布函数 f_post 的 z 向 halo 层，供多 GPU pull streaming 使用。
!===========================================================================================================================
subroutine exchange_f_post_halo_mpi()
  use mpi
  use commondata3dOpenaccMpi
  implicit none
  integer(kind=4) :: i, j, alpha, idx, nBuf

  if (mpiSize .EQ. 1) then
#ifdef SpanwiseWallsPeriodicalU
    !$acc parallel loop collapse(2) present(f_post)
    do j = 1, ny
      do i = 0, nx+1
        do alpha = 0, qf-1
          f_post(alpha,i,j,0) = f_post(alpha,i,j,nz)
          f_post(alpha,i,j,nz+1) = f_post(alpha,i,j,1)
        enddo
      enddo
    enddo
#endif
    return
  endif

  nBuf = qf * (nx + 2) * (ny + 2)
  !$acc update self(f_post(:,:,:,1), f_post(:,:,:,nz))

  idx = 0
  do j = 0, ny+1
    do i = 0, nx+1
      do alpha = 0, qf-1
        idx = idx + 1
        fSendLower(idx) = f_post(alpha,i,j,1)
        fSendUpper(idx) = f_post(alpha,i,j,nz)
      enddo
    enddo
  enddo

  call MPI_Sendrecv(fSendLower, nBuf, MPI_DOUBLE_PRECISION, mpiLeft, 101, &
       fRecvUpper, nBuf, MPI_DOUBLE_PRECISION, mpiRight, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpiErr)
  call MPI_Sendrecv(fSendUpper, nBuf, MPI_DOUBLE_PRECISION, mpiRight, 102, &
       fRecvLower, nBuf, MPI_DOUBLE_PRECISION, mpiLeft, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpiErr)

  if (mpiLeft .NE. MPI_PROC_NULL) then
    idx = 0
    do j = 0, ny+1
      do i = 0, nx+1
        do alpha = 0, qf-1
          idx = idx + 1
          f_post(alpha,i,j,0) = fRecvLower(idx)
        enddo
      enddo
    enddo
  endif

  if (mpiRight .NE. MPI_PROC_NULL) then
    idx = 0
    do j = 0, ny+1
      do i = 0, nx+1
        do alpha = 0, qf-1
          idx = idx + 1
          f_post(alpha,i,j,nz+1) = fRecvUpper(idx)
        enddo
      enddo
    enddo
  endif

  !$acc update device(f_post(:,:,:,0), f_post(:,:,:,nz+1))
end subroutine exchange_f_post_halo_mpi


!===========================================================================================================================
! 子程序: exchange_g_post_halo_mpi
! 作用: 交换温度分布函数 g_post 的 z 向 halo 层，供多 GPU pull streaming 使用。
!===========================================================================================================================
subroutine exchange_g_post_halo_mpi()
  use mpi
  use commondata3dOpenaccMpi
  implicit none
  integer(kind=4) :: i, j, alpha, idx, nBuf

  if (mpiSize .EQ. 1) then
#ifdef SpanwiseWallsPeriodicalT
    !$acc parallel loop collapse(2) present(g_post)
    do j = 1, ny
      do i = 0, nx+1
        do alpha = 0, qt-1
          g_post(alpha,i,j,0) = g_post(alpha,i,j,nz)
          g_post(alpha,i,j,nz+1) = g_post(alpha,i,j,1)
        enddo
      enddo
    enddo
#endif
    return
  endif

  nBuf = qt * (nx + 2) * (ny + 2)
  !$acc update self(g_post(:,:,:,1), g_post(:,:,:,nz))

  idx = 0
  do j = 0, ny+1
    do i = 0, nx+1
      do alpha = 0, qt-1
        idx = idx + 1
        gSendLower(idx) = g_post(alpha,i,j,1)
        gSendUpper(idx) = g_post(alpha,i,j,nz)
      enddo
    enddo
  enddo

  call MPI_Sendrecv(gSendLower, nBuf, MPI_DOUBLE_PRECISION, mpiLeft, 201, &
       gRecvUpper, nBuf, MPI_DOUBLE_PRECISION, mpiRight, 201, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpiErr)
  call MPI_Sendrecv(gSendUpper, nBuf, MPI_DOUBLE_PRECISION, mpiRight, 202, &
       gRecvLower, nBuf, MPI_DOUBLE_PRECISION, mpiLeft, 202, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpiErr)

  if (mpiLeft .NE. MPI_PROC_NULL) then
    idx = 0
    do j = 0, ny+1
      do i = 0, nx+1
        do alpha = 0, qt-1
          idx = idx + 1
          g_post(alpha,i,j,0) = gRecvLower(idx)
        enddo
      enddo
    enddo
  endif

  if (mpiRight .NE. MPI_PROC_NULL) then
    idx = 0
    do j = 0, ny+1
      do i = 0, nx+1
        do alpha = 0, qt-1
          idx = idx + 1
          g_post(alpha,i,j,nz+1) = gRecvUpper(idx)
        enddo
      enddo
    enddo
  endif

  !$acc update device(g_post(:,:,:,0), g_post(:,:,:,nz+1))
end subroutine exchange_g_post_halo_mpi


!===========================================================================================================================
! 子程序: init_lattice_constants_3d
! 作用: 初始化 D3Q19 / D3Q7 的离散速度、反向索引和权重。
!===========================================================================================================================
subroutine init_lattice_constants_3d()
  use commondata3dOpenaccMpi
  implicit none

  ! D3Q19 顺序与 D3Q7 的编号同步（前 7 个速度一样）
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
  use commondata3dOpenaccMpi
  implicit none

  integer(kind=4) :: i, j, k, alpha
  real(kind=8) :: m(0:qf-1), meq(0:qf-1), m_post(0:qf-1)
  real(kind=8) :: s(0:qf-1), fSource(0:qf-1)
  real(kind=8) :: rhoLoc, uLoc, vLoc, wLoc, u2, uDotF
  real(kind=8) :: FxLoc, FyLoc, FzLoc

  !$acc parallel loop collapse(3) present(f,f_post,rho,u,v,w,T)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        rhoLoc = rho(i,j,k)
        uLoc = u(i,j,k)
        vLoc = v(i,j,k)
        wLoc = w(i,j,k)
        u2 = uLoc * uLoc + vLoc * vLoc + wLoc * wLoc

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

        m(3) = f(1,i,j,k) - f(2,i,j,k) + f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) - f(10,i,j,k) + &
             f(11,i,j,k) - f(12,i,j,k) + f(13,i,j,k) - f(14,i,j,k)

        m(4) = -4.0d0 * f(1,i,j,k) + 4.0d0 * f(2,i,j,k) + f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) - &
             f(10,i,j,k) + f(11,i,j,k) - f(12,i,j,k) + f(13,i,j,k) - f(14,i,j,k)

        m(5) = f(3,i,j,k) - f(4,i,j,k) + f(7,i,j,k) + f(8,i,j,k) - f(9,i,j,k) - f(10,i,j,k) + &
             f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) - f(18,i,j,k)

        m(6) = -4.0d0 * f(3,i,j,k) + 4.0d0 * f(4,i,j,k) + f(7,i,j,k) + f(8,i,j,k) - f(9,i,j,k) - &
             f(10,i,j,k) + f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) - f(18,i,j,k)

        m(7) = f(5,i,j,k) - f(6,i,j,k) + f(11,i,j,k) + f(12,i,j,k) - f(13,i,j,k) - f(14,i,j,k) + &
             f(15,i,j,k) + f(16,i,j,k) - f(17,i,j,k) - f(18,i,j,k)

        m(8) = -4.0d0 * f(5,i,j,k) + 4.0d0 * f(6,i,j,k) + f(11,i,j,k) + f(12,i,j,k) - f(13,i,j,k) - &
             f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) - f(17,i,j,k) - f(18,i,j,k)

        m(9) = 2.0d0 * (f(1,i,j,k) + f(2,i,j,k)) - (f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + f(6,i,j,k)) + &
             (f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + f(12,i,j,k) + &
             f(13,i,j,k) + f(14,i,j,k)) - 2.0d0 * (f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k))

        m(10) = -4.0d0 * (f(1,i,j,k) + f(2,i,j,k)) + 2.0d0 * (f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + f(6,i,j,k)) + &
             (f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + f(12,i,j,k) + &
             f(13,i,j,k) + f(14,i,j,k)) - 2.0d0 * (f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k))

        m(11) = (f(3,i,j,k) + f(4,i,j,k)) - (f(5,i,j,k) + f(6,i,j,k)) + &
             (f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k)) - &
             (f(11,i,j,k) + f(12,i,j,k) + f(13,i,j,k) + f(14,i,j,k))

        m(12) = -2.0d0 * (f(3,i,j,k) + f(4,i,j,k)) + 2.0d0 * (f(5,i,j,k) + f(6,i,j,k)) + &
             (f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k)) - &
             (f(11,i,j,k) + f(12,i,j,k) + f(13,i,j,k) + f(14,i,j,k))

        m(13) = f(7,i,j,k) - f(8,i,j,k) - f(9,i,j,k) + f(10,i,j,k)
        m(14) = f(15,i,j,k) - f(16,i,j,k) - f(17,i,j,k) + f(18,i,j,k)
        m(15) = f(11,i,j,k) - f(12,i,j,k) - f(13,i,j,k) + f(14,i,j,k)
        m(16) = f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) - f(10,i,j,k) - f(11,i,j,k) + f(12,i,j,k) - &
             f(13,i,j,k) + f(14,i,j,k)
        m(17) = -f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(15,i,j,k) - f(16,i,j,k) + &
             f(17,i,j,k) - f(18,i,j,k)
        m(18) = f(11,i,j,k) + f(12,i,j,k) - f(13,i,j,k) - f(14,i,j,k) - f(15,i,j,k) - f(16,i,j,k) + &
             f(17,i,j,k) + f(18,i,j,k)

        meq(0) = rhoLoc
        meq(1) = rhoLoc * (-11.0d0 + 19.0d0 * u2)
        meq(2) = rhoLoc * (3.0d0 - 11.0d0 * u2 / 2.0d0)
        meq(3) = rhoLoc * uLoc
        meq(4) = -2.0d0 * rhoLoc * uLoc / 3.0d0
        meq(5) = rhoLoc * vLoc
        meq(6) = -2.0d0 * rhoLoc * vLoc / 3.0d0
        meq(7) = rhoLoc * wLoc
        meq(8) = -2.0d0 * rhoLoc * wLoc / 3.0d0
        meq(9) = rhoLoc * (2.0d0 * uLoc * uLoc - vLoc * vLoc - wLoc * wLoc)
        meq(10) = -0.5d0 * meq(9)
        meq(11) = rhoLoc * (vLoc * vLoc - wLoc * wLoc)
        meq(12) = -0.5d0 * meq(11)
        meq(13) = rhoLoc * uLoc * vLoc
        meq(14) = rhoLoc * vLoc * wLoc
        meq(15) = rhoLoc * uLoc * wLoc
        meq(16) = 0.0d0
        meq(17) = 0.0d0
        meq(18) = 0.0d0

        s(0) = 0.0d0
        s(1) = Se
        s(2) = Seps
        s(3) = 0.0d0
        s(4) = Sq
        s(5) = 0.0d0
        s(6) = Sq
        s(7) = 0.0d0
        s(8) = Sq
        s(9) = Snu
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
        FyLoc = rhoLoc * gBeta * (T(i,j,k) - Tref)
        FzLoc = 0.0d0
        uDotF = uLoc * FxLoc + vLoc * FyLoc + wLoc * FzLoc

        fSource(0) = 0.0d0
        fSource(1) = (1.0d0 - 0.5d0 * s(1)) * 38.0d0 * uDotF
        fSource(2) = (1.0d0 - 0.5d0 * s(2)) * (-11.0d0) * uDotF
        fSource(3) = (1.0d0 - 0.5d0 * s(3)) * FxLoc
        fSource(4) = (1.0d0 - 0.5d0 * s(4)) * (-2.0d0 / 3.0d0) * FxLoc
        fSource(5) = (1.0d0 - 0.5d0 * s(5)) * FyLoc
        fSource(6) = (1.0d0 - 0.5d0 * s(6)) * (-2.0d0 / 3.0d0) * FyLoc
        fSource(7) = (1.0d0 - 0.5d0 * s(7)) * FzLoc
        fSource(8) = (1.0d0 - 0.5d0 * s(8)) * (-2.0d0 / 3.0d0) * FzLoc
        fSource(9) = (1.0d0 - 0.5d0 * s(9)) * (4.0d0 * uLoc * FxLoc - 2.0d0 * vLoc * FyLoc - 2.0d0 * wLoc * FzLoc)
        fSource(10) = (1.0d0 - 0.5d0 * s(10)) * (-2.0d0 * uLoc * FxLoc + vLoc * FyLoc + wLoc * FzLoc)
        fSource(11) = (1.0d0 - 0.5d0 * s(11)) * (2.0d0 * vLoc * FyLoc - 2.0d0 * wLoc * FzLoc)
        fSource(12) = (1.0d0 - 0.5d0 * s(12)) * (-vLoc * FyLoc + wLoc * FzLoc)
        fSource(13) = (1.0d0 - 0.5d0 * s(13)) * (uLoc * FyLoc + vLoc * FxLoc)
        fSource(14) = (1.0d0 - 0.5d0 * s(14)) * (vLoc * FzLoc + wLoc * FyLoc)
        fSource(15) = (1.0d0 - 0.5d0 * s(15)) * (uLoc * FzLoc + wLoc * FxLoc)
        fSource(16) = 0.0d0
        fSource(17) = 0.0d0
        fSource(18) = 0.0d0

        do alpha = 0, qf-1
          m_post(alpha) = m(alpha) - s(alpha) * (m(alpha) - meq(alpha)) + fSource(alpha)
        enddo

        f_post(0,i,j,k) = m_post(0) / 19.0d0 - 5.0d0 * m_post(1) / 399.0d0 + m_post(2) / 21.0d0

        f_post(1,i,j,k) = m_post(0) / 19.0d0 - 11.0d0 * m_post(1) / 2394.0d0 - m_post(2) / 63.0d0 + &
             m_post(3) / 10.0d0 - m_post(4) / 10.0d0 + m_post(9) / 18.0d0 - m_post(10) / 18.0d0

        f_post(2,i,j,k) = m_post(0) / 19.0d0 - 11.0d0 * m_post(1) / 2394.0d0 - m_post(2) / 63.0d0 - &
             m_post(3) / 10.0d0 + m_post(4) / 10.0d0 + m_post(9) / 18.0d0 - m_post(10) / 18.0d0

        f_post(3,i,j,k) = m_post(0) / 19.0d0 - 11.0d0 * m_post(1) / 2394.0d0 - m_post(2) / 63.0d0 + &
             m_post(5) / 10.0d0 - m_post(6) / 10.0d0 - m_post(9) / 36.0d0 + m_post(10) / 36.0d0 + &
             m_post(11) / 12.0d0 - m_post(12) / 12.0d0

        f_post(4,i,j,k) = m_post(0) / 19.0d0 - 11.0d0 * m_post(1) / 2394.0d0 - m_post(2) / 63.0d0 - &
             m_post(5) / 10.0d0 + m_post(6) / 10.0d0 - m_post(9) / 36.0d0 + m_post(10) / 36.0d0 + &
             m_post(11) / 12.0d0 - m_post(12) / 12.0d0

        f_post(5,i,j,k) = m_post(0) / 19.0d0 - 11.0d0 * m_post(1) / 2394.0d0 - m_post(2) / 63.0d0 + &
             m_post(7) / 10.0d0 - m_post(8) / 10.0d0 - m_post(9) / 36.0d0 + m_post(10) / 36.0d0 - &
             m_post(11) / 12.0d0 + m_post(12) / 12.0d0

        f_post(6,i,j,k) = m_post(0) / 19.0d0 - 11.0d0 * m_post(1) / 2394.0d0 - m_post(2) / 63.0d0 - &
             m_post(7) / 10.0d0 + m_post(8) / 10.0d0 - m_post(9) / 36.0d0 + m_post(10) / 36.0d0 - &
             m_post(11) / 12.0d0 + m_post(12) / 12.0d0

        f_post(7,i,j,k) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 + &
             m_post(3) / 10.0d0 + m_post(4) / 40.0d0 + m_post(5) / 10.0d0 + m_post(6) / 40.0d0 + &
             m_post(9) / 36.0d0 + m_post(10) / 72.0d0 + m_post(11) / 12.0d0 + m_post(12) / 24.0d0 + &
             m_post(13) / 4.0d0 + m_post(16) / 8.0d0 - m_post(17) / 8.0d0

        f_post(8,i,j,k) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 - &
             m_post(3) / 10.0d0 - m_post(4) / 40.0d0 + m_post(5) / 10.0d0 + m_post(6) / 40.0d0 + &
             m_post(9) / 36.0d0 + m_post(10) / 72.0d0 + m_post(11) / 12.0d0 + m_post(12) / 24.0d0 - &
             m_post(13) / 4.0d0 - m_post(16) / 8.0d0 - m_post(17) / 8.0d0

        f_post(9,i,j,k) = m_post(0) / 19.0d0 + 4.0d0 * m_post(1) / 1197.0d0 + m_post(2) / 252.0d0 + &
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

end subroutine collision3d


!===========================================================================================================================
! 子程序: fill_periodic_ghosts_f_post
! 作用: 为流场分布函数补齐周期 ghost 层，便于统一 pull streaming。
!===========================================================================================================================
subroutine fill_periodic_ghosts_f_post()
  use commondata3dOpenaccMpi
  implicit none

  integer(kind=4) :: i, j, k, alpha

#ifdef VerticalWallsPeriodicalU
  ! x 方向周期仍然在各 rank 内局部补 ghost 层。
 !$acc parallel loop collapse(2) present(f_post)
  do k = 1, nz
    do j = 1, ny
      do alpha = 0, qf-1
        f_post(alpha,0,j,k) = f_post(alpha,nx,j,k)
        f_post(alpha,nx+1,j,k) = f_post(alpha,1,j,k)
      enddo
    enddo
  enddo
#endif

  ! z 方向由 MPI halo 交换补齐；当只有一个 rank 且 z 周期时，会退化成本地复制。
  call exchange_f_post_halo_mpi()

end subroutine fill_periodic_ghosts_f_post


!===========================================================================================================================
! 子程序: streaming3d
! 作用: 对流场分布函数执行三维 pull streaming。
!===========================================================================================================================
subroutine streaming3d()
  use commondata3dOpenaccMpi
  implicit none

  integer(kind=4) :: i, j, k, ip, jp, kp, alpha

  ! pull streaming：当前格点 (i,j,k) 从上游格点 (i-ex, j-ey, k-ez) 拉取分布函数
 !$acc parallel loop collapse(3) present(f,f_post,ex,ey,ez)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        do alpha = 0, qf-1
          ip = i - ex(alpha)
          jp = j - ey(alpha)
          kp = k - ez(alpha)
          f(alpha,i,j,k) = f_post(alpha,ip,jp,kp)
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
  use commondata3dOpenaccMpi
  implicit none

  integer(kind=4) :: i, j, k, alpha

#ifdef VerticalWallsPeriodicalU
 !$acc parallel loop collapse(2) present(f,f_post,ex)
  do k = 1, nz
    do j = 1, ny
      do alpha = 0, qf-1
        if (ex(alpha) .EQ. 1)  f(alpha,1,j,k)  = f_post(alpha,nx,j,k)
        if (ex(alpha) .EQ. -1) f(alpha,nx,j,k) = f_post(alpha,1,j,k)
      enddo
    enddo
  enddo
#endif

#ifdef VerticalWallsNoslip
 !$acc parallel loop collapse(2) present(f,f_post,ex,opp)
  do k = 1, nz
    do j = 1, ny
      do alpha = 0, qf-1
        if (ex(alpha) .EQ. 1)  f(alpha,1,j,k)  = f_post(opp(alpha),1,j,k)
        if (ex(alpha) .EQ. -1) f(alpha,nx,j,k) = f_post(opp(alpha),nx,j,k)
      enddo
    enddo
  enddo
#endif

#ifdef HorizontalWallsNoslip
 !$acc parallel loop collapse(2) present(f,f_post,ey,opp)
  do k = 1, nz
    do i = 1, nx
      do alpha = 0, qf-1
        if (ey(alpha) .EQ. 1)  f(alpha,i,1,k)  = f_post(opp(alpha),i,1,k)
        if (ey(alpha) .EQ. -1) f(alpha,i,ny,k) = f_post(opp(alpha),i,ny,k)
      enddo
    enddo
  enddo
#endif

#ifdef SpanwiseWallsPeriodicalU
  if (isFirstZRank) then
 !$acc parallel loop collapse(2) present(f,f_post,ez)
    do j = 1, ny
      do i = 1, nx
        do alpha = 0, qf-1
          if (ez(alpha) .EQ. 1)  f(alpha,i,j,1)  = f_post(alpha,i,j,0)
        enddo
      enddo
    enddo
  endif
  if (isLastZRank) then
 !$acc parallel loop collapse(2) present(f,f_post,ez)
    do j = 1, ny
      do i = 1, nx
        do alpha = 0, qf-1
          if (ez(alpha) .EQ. -1) f(alpha,i,j,nz) = f_post(alpha,i,j,nz+1)
        enddo
      enddo
    enddo
  endif
#endif

#ifdef SpanwiseWallsNoslip
  if (isFirstZRank) then
 !$acc parallel loop collapse(2) present(f,f_post,ez,opp)
    do j = 1, ny
      do i = 1, nx
        do alpha = 0, qf-1
          if (ez(alpha) .EQ. 1)  f(alpha,i,j,1)  = f_post(opp(alpha),i,j,1)
        enddo
      enddo
    enddo
  endif
  if (isLastZRank) then
 !$acc parallel loop collapse(2) present(f,f_post,ez,opp)
    do j = 1, ny
      do i = 1, nx
        do alpha = 0, qf-1
          if (ez(alpha) .EQ. -1) f(alpha,i,j,nz) = f_post(opp(alpha),i,j,nz)
        enddo
      enddo
    enddo
  endif
#endif

end subroutine bounceback3d


!===========================================================================================================================
! 子程序: macro3d
! 作用: 由流场分布函数恢复 rho、u、v、w 以及浮力项。
!===========================================================================================================================
subroutine macro3d()
  use commondata3dOpenaccMpi
  implicit none

  integer(kind=4) :: i, j, k
  real(kind=8) :: FyLoc

  !$acc parallel loop collapse(3) present(f,rho,u,v,w,T)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        rho(i,j,k) = f(0,i,j,k) + f(1,i,j,k) + f(2,i,j,k) + f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + &
             f(6,i,j,k) + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + &
             f(12,i,j,k) + f(13,i,j,k) + f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k)

        FyLoc = rho(i,j,k) * gBeta * (T(i,j,k) - Tref)

        u(i,j,k) = (f(1,i,j,k) - f(2,i,j,k) + f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) - f(10,i,j,k) + &
             f(11,i,j,k) - f(12,i,j,k) + f(13,i,j,k) - f(14,i,j,k)) / rho(i,j,k)

        v(i,j,k) = (f(3,i,j,k) - f(4,i,j,k) + f(7,i,j,k) + f(8,i,j,k) - f(9,i,j,k) - f(10,i,j,k) + &
             f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) - f(18,i,j,k) + 0.5d0 * FyLoc) / rho(i,j,k)

        w(i,j,k) = (f(5,i,j,k) - f(6,i,j,k) + f(11,i,j,k) + f(12,i,j,k) - f(13,i,j,k) - f(14,i,j,k) + &
             f(15,i,j,k) + f(16,i,j,k) - f(17,i,j,k) - f(18,i,j,k)) / rho(i,j,k)
      enddo
    enddo
  enddo

end subroutine macro3d


!===========================================================================================================================
! 子程序: collisionT3d
! 作用: 温度场碰撞步骤，算法口径尽量保持与 2DRB 的对流扩散处理一致。
!===========================================================================================================================
subroutine collisionT3d()
  use commondata3dOpenaccMpi
  implicit none

  integer(kind=4) :: i, j, k
  real(kind=8) :: n(0:qt-1), neq(0:qt-1), n_post(0:qt-1), q(0:qt-1)
  real(kind=8) :: Bx, By, Bz, dBx, dBy, dBz
  real(kind=8), parameter :: SG = 1.0d0 - 0.5d0 * Qk

  !$acc parallel loop collapse(3) present(g,g_post,u,v,w,T,Bx_prev,By_prev,Bz_prev)
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

        n(0) = g(0,i,j,k) + g(1,i,j,k) + g(2,i,j,k) + g(3,i,j,k) + g(4,i,j,k) + g(5,i,j,k) + g(6,i,j,k)
        n(1) = g(1,i,j,k) - g(2,i,j,k)
        n(2) = g(3,i,j,k) - g(4,i,j,k)
        n(3) = g(5,i,j,k) - g(6,i,j,k)
        n(4) = -6.0d0 * g(0,i,j,k) + g(1,i,j,k) + g(2,i,j,k) + g(3,i,j,k) + g(4,i,j,k) + g(5,i,j,k) + g(6,i,j,k)
        n(5) = 2.0d0 * g(1,i,j,k) + 2.0d0 * g(2,i,j,k) - g(3,i,j,k) - g(4,i,j,k) - g(5,i,j,k) - g(6,i,j,k)
        n(6) = g(3,i,j,k) + g(4,i,j,k) - g(5,i,j,k) - g(6,i,j,k)

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

end subroutine collisionT3d


!===========================================================================================================================
!===========================================================================================================================
! 子程序: fill_periodic_ghosts_g_post
! 作用: 为温度分布函数补齐周期 ghost 层，便于统一 pull streaming。
!===========================================================================================================================
subroutine fill_periodic_ghosts_g_post()
  use commondata3dOpenaccMpi
  implicit none

  integer(kind=4) :: i, j, k, alpha

#ifdef VerticalWallsPeriodicalT
 !$acc parallel loop collapse(2) present(g_post)
  do k = 1, nz
    do j = 1, ny
      do alpha = 0, qt-1
        g_post(alpha,0,j,k) = g_post(alpha,nx,j,k)
        g_post(alpha,nx+1,j,k) = g_post(alpha,1,j,k)
      enddo
    enddo
  enddo
#endif

  call exchange_g_post_halo_mpi()

end subroutine fill_periodic_ghosts_g_post


!===========================================================================================================================
! 子程序: streamingT3d
! 作用: 对温度分布函数执行三维 pull streaming。
!===========================================================================================================================
subroutine streamingT3d()
  use commondata3dOpenaccMpi
  implicit none

  integer(kind=4) :: i, j, k, ip, jp, kp, alpha

 !$acc parallel loop collapse(3) present(g,g_post,exT,eyT,ezT)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        do alpha = 0, qt-1
          ip = i - exT(alpha)
          jp = j - eyT(alpha)
          kp = k - ezT(alpha)
          g(alpha,i,j,k) = g_post(alpha,ip,jp,kp)
        enddo
      enddo
    enddo
  enddo

end subroutine streamingT3d


!===========================================================================================================================
! 子程序: bouncebackT3d
! 作用: 施加温度边界条件，包括恒温、绝热和周期边界。
!===========================================================================================================================
subroutine bouncebackT3d()
  use commondata3dOpenaccMpi
  implicit none

  integer(kind=4) :: i, j, k, alpha

#ifdef VerticalWallsPeriodicalT
 !$acc parallel loop collapse(2) present(g,g_post,exT)
  do k = 1, nz
    do j = 1, ny
      do alpha = 0, qt-1
        if (exT(alpha) .EQ. 1)  g(alpha,1,j,k)  = g_post(alpha,nx,j,k)
        if (exT(alpha) .EQ. -1) g(alpha,nx,j,k) = g_post(alpha,1,j,k)
      enddo
    enddo
  enddo
#endif

#ifdef VerticalWallsConstT
 !$acc parallel loop collapse(2) present(g,g_post,exT,oppT,omegaT)
  do k = 1, nz
    do j = 1, ny
      do alpha = 0, qt-1
#ifdef EnableLegacyThermalScheme
        if (exT(alpha) .EQ. 1)  g(alpha,1,j,k)  = -g_post(oppT(alpha),1,j,k)  + (6.0d0 + paraA) / 21.0d0 * Thot
        if (exT(alpha) .EQ. -1) g(alpha,nx,j,k) = -g_post(oppT(alpha),nx,j,k) + (6.0d0 + paraA) / 21.0d0 * Tcold
#else
        if (exT(alpha) .EQ. 1)  g(alpha,1,j,k)  = -g_post(oppT(alpha),1,j,k)  + 2.0d0 * omegaT(alpha) * Thot
        if (exT(alpha) .EQ. -1) g(alpha,nx,j,k) = -g_post(oppT(alpha),nx,j,k) + 2.0d0 * omegaT(alpha) * Tcold
#endif
      enddo
    enddo
  enddo
#endif

#ifdef VerticalWallsAdiabatic
 !$acc parallel loop collapse(2) present(g,g_post,exT,oppT)
  do k = 1, nz
    do j = 1, ny
      do alpha = 0, qt-1
        if (exT(alpha) .EQ. 1)  g(alpha,1,j,k)  = g_post(oppT(alpha),1,j,k)
        if (exT(alpha) .EQ. -1) g(alpha,nx,j,k) = g_post(oppT(alpha),nx,j,k)
      enddo
    enddo
  enddo
#endif

#ifdef HorizontalWallsAdiabatic
 !$acc parallel loop collapse(2) present(g,g_post,eyT,oppT)
  do k = 1, nz
    do i = 1, nx
      do alpha = 0, qt-1
        if (eyT(alpha) .EQ. 1)  g(alpha,i,1,k)  = g_post(oppT(alpha),i,1,k)
        if (eyT(alpha) .EQ. -1) g(alpha,i,ny,k) = g_post(oppT(alpha),i,ny,k)
      enddo
    enddo
  enddo
#endif

#ifdef HorizontalWallsConstT
 !$acc parallel loop collapse(2) present(g,g_post,eyT,oppT,omegaT)
  do k = 1, nz
    do i = 1, nx
      do alpha = 0, qt-1
#ifdef EnableLegacyThermalScheme
        if (eyT(alpha) .EQ. 1)  g(alpha,i,1,k)  = -g_post(oppT(alpha),i,1,k)  + (6.0d0 + paraA) / 21.0d0 * Thot
        if (eyT(alpha) .EQ. -1) g(alpha,i,ny,k) = -g_post(oppT(alpha),i,ny,k) + (6.0d0 + paraA) / 21.0d0 * Tcold
#else
        if (eyT(alpha) .EQ. 1)  g(alpha,i,1,k)  = -g_post(oppT(alpha),i,1,k)  + 2.0d0 * omegaT(alpha) * Thot
        if (eyT(alpha) .EQ. -1) g(alpha,i,ny,k) = -g_post(oppT(alpha),i,ny,k) + 2.0d0 * omegaT(alpha) * Tcold
#endif
      enddo
    enddo
  enddo
#endif

#ifdef SpanwiseWallsPeriodicalT
  if (isFirstZRank) then
 !$acc parallel loop collapse(2) present(g,g_post,ezT)
    do j = 1, ny
      do i = 1, nx
        do alpha = 0, qt-1
          if (ezT(alpha) .EQ. 1)  g(alpha,i,j,1)  = g_post(alpha,i,j,0)
        enddo
      enddo
    enddo
  endif
  if (isLastZRank) then
 !$acc parallel loop collapse(2) present(g,g_post,ezT)
    do j = 1, ny
      do i = 1, nx
        do alpha = 0, qt-1
          if (ezT(alpha) .EQ. -1) g(alpha,i,j,nz) = g_post(alpha,i,j,nz+1)
        enddo
      enddo
    enddo
  endif
#endif

#ifdef SpanwiseWallsAdiabatic
  if (isFirstZRank) then
 !$acc parallel loop collapse(2) present(g,g_post,ezT,oppT)
    do j = 1, ny
      do i = 1, nx
        do alpha = 0, qt-1
          if (ezT(alpha) .EQ. 1)  g(alpha,i,j,1)  = g_post(oppT(alpha),i,j,1)
        enddo
      enddo
    enddo
  endif
  if (isLastZRank) then
 !$acc parallel loop collapse(2) present(g,g_post,ezT,oppT)
    do j = 1, ny
      do i = 1, nx
        do alpha = 0, qt-1
          if (ezT(alpha) .EQ. -1) g(alpha,i,j,nz) = g_post(oppT(alpha),i,j,nz)
        enddo
      enddo
    enddo
  endif
#endif

  return
end subroutine bouncebackT3d


!===========================================================================================================================
! 子程序: macroT3d
! 作用: 由温度分布函数恢复温度场，并更新历史热流项。
!===========================================================================================================================
subroutine macroT3d()
  use commondata3dOpenaccMpi
  implicit none

  integer(kind=4) :: i, j, k

 !$acc parallel loop collapse(3) present(g,T)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        T(i,j,k) = g(0,i,j,k) + g(1,i,j,k) + g(2,i,j,k) + g(3,i,j,k) + &
                   g(4,i,j,k) + g(5,i,j,k) + g(6,i,j,k)
      enddo
    enddo
  enddo

  return
end subroutine macroT3d

!===========================================================================================================================
! 子程序: reconstruct_macro_from_fg3d
! 作用: 从重启读回的 f/g 重新恢复宏观场，避免备份文件格式过重。
!===========================================================================================================================
subroutine reconstruct_macro_from_fg3d()
  use commondata3dOpenaccMpi
  implicit none

  integer(kind=4) :: i, j, k, alpha, iter
  real(kind=8) :: momx, momy, momz
  real(kind=8) :: FxLoc, FyLoc, FzLoc
  logical :: rho_bad

  ! 重启时只存了 f 和 g，所以这里统一把 T、rho、u、v、w 以及历史热流都重构回来
  ! 这一步发生在 enter_data_3d_openacc() 之前，因此保持主机端重构更稳妥。
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        T(i,j,k) = 0.0d0
        do alpha = 0, qt-1
          T(i,j,k) = T(i,j,k) + g(alpha,i,j,k)
        enddo
      enddo
    enddo
  enddo

  rho_bad = .false.
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        rho(i,j,k) = 0.0d0
        momx = 0.0d0
        momy = 0.0d0
        momz = 0.0d0
        do alpha = 0, qf-1
          rho(i,j,k) = rho(i,j,k) + f(alpha,i,j,k)
          momx = momx + f(alpha,i,j,k) * dble(ex(alpha))
          momy = momy + f(alpha,i,j,k) * dble(ey(alpha))
          momz = momz + f(alpha,i,j,k) * dble(ez(alpha))
        enddo

        if (rho(i,j,k) .GT. 0.0d0) then
          u(i,j,k) = momx / rho(i,j,k)
          v(i,j,k) = momy / rho(i,j,k)
          w(i,j,k) = momz / rho(i,j,k)
          do iter = 1, 3
            FxLoc = 0.0d0
            FyLoc = rho(i,j,k) * gBeta * (T(i,j,k) - Tref)
            FzLoc = 0.0d0
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

#ifdef EnableUseG
        Bx_prev(i,j,k) = u(i,j,k) * T(i,j,k)
        By_prev(i,j,k) = v(i,j,k) * T(i,j,k)
        Bz_prev(i,j,k) = w(i,j,k) * T(i,j,k)
#else
        Bx_prev(i,j,k) = 0.0d0
        By_prev(i,j,k) = 0.0d0
        Bz_prev(i,j,k) = 0.0d0
#endif
      enddo
    enddo
  enddo

  if (rho_bad) then
    write(*,*) 'Warning: non-positive rho found during restart reconstruction.'
    stop
  endif

  return
end subroutine reconstruct_macro_from_fg3d


#ifdef steadyFlow
!===========================================================================================================================
! 子程序: check3d
! 作用: 计算稳态收敛误差，并按需写入收敛历史。
!===========================================================================================================================
subroutine check3d()
  use mpi
  use commondata3dOpenaccMpi
  implicit none

  integer(kind=4) :: i, j, k
  real(kind=8) :: error1, error2, error5, error6
  character(len=80) :: caseTag

  error1 = 0.0d0
  error2 = 0.0d0
  error5 = 0.0d0
  error6 = 0.0d0

 !$acc parallel loop collapse(3) &
 !$acc& present(u,up,v,vp,w,wp,T,Tp) &
 !$acc& reduction(+:error1,error2,error5,error6)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        error1 = error1 + (u(i,j,k)-up(i,j,k))*(u(i,j,k)-up(i,j,k)) + &
                        (v(i,j,k)-vp(i,j,k))*(v(i,j,k)-vp(i,j,k)) + &
                        (w(i,j,k)-wp(i,j,k))*(w(i,j,k)-wp(i,j,k))
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

  call MPI_Allreduce(MPI_IN_PLACE, error1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  call MPI_Allreduce(MPI_IN_PLACE, error2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  call MPI_Allreduce(MPI_IN_PLACE, error5, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  call MPI_Allreduce(MPI_IN_PLACE, error6, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpiErr)

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

  if (isRoot) then
    call append_convergence_tecplot3d('convergence3D_mpi.plt', itc, errorU, errorT)
    write(caseTag,'("Ra=",ES10.3E2,",nx=",I0,",ny=",I0,",nz=",I0,",useG=",L1,",old=",L1,",mpi=",I0)') &
         Rayleigh, nx, ny, nzGlobal, useG, useLegacyThermalScheme, mpiSize
    call append_convergence_master_tecplot3d('convergence_all_3D_mpi.plt', caseTag, itc, errorU, errorT)
    write(*,'(I12,1X,ES24.16,1X,ES24.16)') itc, errorU, errorT
  endif

  return
end subroutine check3d
#endif


!===========================================================================================================================
! 子程序: append_convergence_tecplot3d
! 作用: 向单个收敛历史文件追加一条误差记录。
!===========================================================================================================================
subroutine append_convergence_tecplot3d(filename, itcLoc, errorULoc, errorTLoc)
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4), intent(in) :: itcLoc
  real(kind=8), intent(in) :: errorULoc, errorTLoc
  integer(kind=4) :: u
  logical, save :: first_write = .true.

  if (first_write) then
    open(newunit=u, file=trim(filename), status='replace', action='write', form='formatted')
    write(u,'(A)') 'VARIABLES = "itc" "errorU" "errorT"'
    write(u,'(A)') 'ZONE T="conv3d", F=POINT'
    write(u,'(I12,1X,ES24.16,1X,ES24.16)') itcLoc, errorULoc, errorTLoc
    close(u)
    first_write = .false.
  else
    open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
    write(u,'(I12,1X,ES24.16,1X,ES24.16)') itcLoc, errorULoc, errorTLoc
    close(u)
  endif

end subroutine append_convergence_tecplot3d


!===========================================================================================================================
! 子程序: append_convergence_master_tecplot3d
! 作用: 向带 zone 名称的收敛历史文件追加一条记录。
!===========================================================================================================================
subroutine append_convergence_master_tecplot3d(filename, zoneName, itcLoc, errorULoc, errorTLoc)
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
      open(newunit=u, file=trim(filename), status='replace', action='write', form='formatted')
      write(u,'(A)') 'TITLE = "Convergence comparison 3D"'
      write(u,'(A)') 'VARIABLES = "itc" "errorU" "errorT"'
      close(u)
    endif
    open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
    write(u,'(A)') 'ZONE T="'//trim(zoneName)//'", F=POINT'
    close(u)
    zone_started = .true.
  endif

  open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
  write(u,'(I12,1X,ES24.16,1X,ES24.16)') itcLoc, errorULoc, errorTLoc
  close(u)

end subroutine append_convergence_master_tecplot3d


!===========================================================================================================================
! 子程序: output_binary3d
! 作用: 输出三维快照二进制文件，供后处理或继续分析使用。
!===========================================================================================================================
subroutine output_binary3d()
  use commondata3dOpenaccMpi
  implicit none

  integer(kind=4) :: i, j, k
  character(len=100) :: filename

  ! 这是给后处理看的快照文件
  ! 输出的是已经乘上 velocityScaleCompare 的无量纲速度场
#ifdef steadyFlow
  write(filename,'(i12.12)') itc
#endif
#ifdef unsteadyFlow
  binFileNum = binFileNum + 1
  if (loadInitField .EQ. 0) write(filename,'(i12.12)') binFileNum
  if (loadInitField .EQ. 1) write(filename,'(i12.12)') binFileNum + reloadbinFileNum
#endif

  filename = adjustl(filename)
  open(unit=03, file=trim(binFilePrefix)//'-'//trim(filename)//'.bin', form='unformatted', access='sequential')
  write(03) (((velocityScaleCompare*u(i,j,k), i=1,nx), j=1,ny), k=1,nz)
  write(03) (((velocityScaleCompare*v(i,j,k), i=1,nx), j=1,ny), k=1,nz)
  write(03) (((velocityScaleCompare*w(i,j,k), i=1,nx), j=1,ny), k=1,nz)
  write(03) (((T(i,j,k), i=1,nx), j=1,ny), k=1,nz)
  write(03) (((rho(i,j,k), i=1,nx), j=1,ny), k=1,nz)
  close(03)

end subroutine output_binary3d


!===========================================================================================================================
! 子程序: backupData3d
! 作用: 备份 f/g 分布函数，供后续重启继续计算。
!===========================================================================================================================
subroutine backupData3d()
  use commondata3dOpenaccMpi
  implicit none

  integer(kind=4) :: i, j, k, alpha
  character(len=100) :: filename

  ! 这是严格重启文件
  ! 只写 f 和 g，后续由 reconstruct_macro_from_fg3d() 恢复宏观量
#ifdef steadyFlow
  write(filename,'(i0)') itc
#endif
#ifdef unsteadyFlow
  if (loadInitField .EQ. 0) write(filename,'(i0)') binFileNum
  if (loadInitField .EQ. 1) write(filename,'(i0)') binFileNum + reloadbinFileNum
#endif

  filename = adjustl(filename)
  open(unit=05, file=trim(reloadFilePrefix)//'-'//trim(filename)//'.bin', form='unformatted', access='sequential')
  write(05) ((((f(alpha,i,j,k), i=1,nx), j=1,ny), k=1,nz), alpha=0,qf-1)
  write(05) ((((g(alpha,i,j,k), i=1,nx), j=1,ny), k=1,nz), alpha=0,qt-1)
  close(05)

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) 'Backup f and g to the file: ', trim(reloadFilePrefix), '-', trim(filename), '.bin'
  close(00)

end subroutine backupData3d


!===========================================================================================================================
! 子程序: output_Tecplot3d
! 作用: MPI 版本当前仍只输出中面切片，这里保留和 3DRB 一致的调用入口。
!===========================================================================================================================
subroutine output_Tecplot3d()
  call output_midplanes_tecplot3d()
  return
end subroutine output_Tecplot3d


!===========================================================================================================================
! 子程序: output_midplanes_tecplot3d
! 作用: 输出 x/y/z 三个中面切片，便于快速查看三维流场结构。
!===========================================================================================================================
subroutine output_midplanes_tecplot3d()
  use commondata3dOpenaccMpi
  implicit none

  character(len=100) :: tag

  ! 3D 第一版不输出整体 Tecplot 体数据，只输出三个中面切片便于观察主流型和温度分布
#ifdef steadyFlow
  write(tag,'(i12.12)') itc
#endif
#ifdef unsteadyFlow
  pltFileNum = pltFileNum + 1
  write(tag,'(i12.12)') pltFileNum
#endif

  tag = adjustl(tag)
  call write_midplane_x(trim(pltFilePrefix)//'-midX-'//trim(tag)//'.dat')
  call write_midplane_y(trim(pltFilePrefix)//'-midY-'//trim(tag)//'.dat')
  call write_midplane_z(trim(pltFilePrefix)//'-midZ-'//trim(tag)//'.dat')

end subroutine output_midplanes_tecplot3d
subroutine calNuRe3d()
  use mpi
  use commondata3dOpenaccMpi
  implicit none

  integer(kind=4) :: i, j, k
  real(kind=8) :: NuVolAvg_temp, ReVolAvg_temp

  if (dimensionlessTime .GE. dimensionlessTimeMax) return
  dimensionlessTime = dimensionlessTime + 1

  NuVolAvg_temp = 0.0d0
#ifdef SideHeatedCell
 !$acc parallel loop collapse(3) present(u,T) reduction(+:NuVolAvg_temp)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        NuVolAvg_temp = NuVolAvg_temp + u(i,j,k) * T(i,j,k)
      enddo
    enddo
  enddo
#else
 !$acc parallel loop collapse(3) present(v,T) reduction(+:NuVolAvg_temp)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        NuVolAvg_temp = NuVolAvg_temp + v(i,j,k) * T(i,j,k)
      enddo
    enddo
  enddo
#endif
  call MPI_Allreduce(MPI_IN_PLACE, NuVolAvg_temp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  NuVolAvg(dimensionlessTime) = NuVolAvg_temp / dble(nx * ny * nzGlobal) * lengthUnit / diffusivity + 1.0d0

  ReVolAvg_temp = 0.0d0
 !$acc parallel loop collapse(3) present(u,v,w) reduction(+:ReVolAvg_temp)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        ReVolAvg_temp = ReVolAvg_temp + u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)
      enddo
    enddo
  enddo
  call MPI_Allreduce(MPI_IN_PLACE, ReVolAvg_temp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  ReVolAvg(dimensionlessTime) = dsqrt(ReVolAvg_temp / dble(nx * ny * nzGlobal)) * lengthUnit / viscosity

  if (isRoot) then
    open(unit=01, file='Nu_VolAvg_3D_MPI.dat', status='unknown', position='append')
    write(01,*) real(reloadDimensionlessTime + dimensionlessTime * outputFrequency, kind=8), NuVolAvg(dimensionlessTime)
    close(01)
    open(unit=02, file='Re_VolAvg_3D_MPI.dat', status='unknown', position='append')
    write(02,*) real(reloadDimensionlessTime + dimensionlessTime * outputFrequency, kind=8), ReVolAvg(dimensionlessTime)
    close(02)
    write(*,'(a,1x,es16.8)') 'NuVolAvg =', NuVolAvg(dimensionlessTime)
    write(*,'(a,1x,es16.8)') 'ReVolAvg =', ReVolAvg(dimensionlessTime)
  endif

end subroutine calNuRe3d


subroutine SideHeatedcalc_Nu_global3d()
  call RBcalc_Nu_global3d()
end subroutine SideHeatedcalc_Nu_global3d

subroutine SideHeatedcalc_Nu_wall_avg3d()
  call RBcalc_Nu_wall_avg3d()
end subroutine SideHeatedcalc_Nu_wall_avg3d


subroutine RBcalc_Nu_global3d()
  use mpi
  use commondata3dOpenaccMpi
  implicit none
  integer(kind=4) :: i, j, k
  real(kind=8) :: Nu_sum, dx, deltaT, coef, dT, qx

#ifdef SideHeatedCell
  dx = 1.0d0 / lengthUnit
  deltaT = Thot - Tcold
  coef = heatFluxScale
  Nu_sum = 0.0d0
 !$acc parallel loop collapse(3) present(u,T) reduction(+:Nu_sum)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        if (i .EQ. 1) then
          dT = (-3.0d0*T(1,j,k) - T(2,j,k) + 4.0d0*Thot) / (3.0d0*dx)
        elseif (i .EQ. nx) then
          dT = (-4.0d0*Tcold + 3.0d0*T(nx,j,k) + T(nx-1,j,k)) / (3.0d0*dx)
        else
          dT = (T(i-1,j,k) - T(i+1,j,k)) / (2.0d0*dx)
        endif
        qx = coef * u(i,j,k) * (T(i,j,k) - Tref) + dT
        Nu_sum = Nu_sum + qx
      enddo
    enddo
  enddo
  call MPI_Allreduce(MPI_IN_PLACE, Nu_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  Nu_global = (Nu_sum / dble(nx * ny * nzGlobal)) / deltaT
#else
  Nu_sum = 0.0d0
 !$acc parallel loop collapse(3) present(v,T) reduction(+:Nu_sum)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        Nu_sum = Nu_sum + v(i,j,k) * T(i,j,k)
      enddo
    enddo
  enddo
  call MPI_Allreduce(MPI_IN_PLACE, Nu_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  Nu_global = 1.0d0 + heatFluxScale * Nu_sum / dble(nx * ny * nzGlobal)
#endif

  if (isRoot) then
    write(*,'(a,1x,es16.8)') 'Nu_global =', Nu_global
    open(unit=00, file=trim(settingsFile), status='unknown', position='append')
    write(00,'(a,1x,es16.8)') 'Nu_global =', Nu_global
    close(00)
  endif

end subroutine RBcalc_Nu_global3d


subroutine RBcalc_Nu_wall_avg3d()
  use mpi
  use commondata3dOpenaccMpi
  implicit none
  integer(kind=4) :: i, j, k, iL, iR, iMid, jB, jT
  real(kind=8) :: dx, dy, deltaT, coef
  real(kind=8) :: sum_hot, sum_cold, sum_mid, q_wall

#ifdef SideHeatedCell
  dx = 1.0d0 / lengthUnit
  deltaT = Thot - Tcold
  coef = heatFluxScale

  sum_hot = 0.0d0
 !$acc parallel loop collapse(2) present(T) reduction(+:sum_hot)
  do k = 1, nz
    do j = 1, ny
      q_wall = 2.0d0 * (Thot - T(1,j,k)) / dx
      sum_hot = sum_hot + q_wall / deltaT
    enddo
  enddo
  call MPI_Allreduce(MPI_IN_PLACE, sum_hot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  Nu_hot = sum_hot / dble(ny * nzGlobal)

  sum_cold = 0.0d0
 !$acc parallel loop collapse(2) present(T) reduction(+:sum_cold)
  do k = 1, nz
    do j = 1, ny
      q_wall = 2.0d0 * (T(nx,j,k) - Tcold) / dx
      sum_cold = sum_cold + q_wall / deltaT
    enddo
  enddo
  call MPI_Allreduce(MPI_IN_PLACE, sum_cold, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  Nu_cold = sum_cold / dble(ny * nzGlobal)

  sum_mid = 0.0d0
  if (mod(nx,2) .EQ. 1) then
    iMid = (nx + 1) / 2
 !$acc parallel loop collapse(2) present(u,T) reduction(+:sum_mid)
    do k = 1, nz
      do j = 1, ny
        sum_mid = sum_mid + (coef * u(iMid,j,k) * (T(iMid,j,k) - Tref) + &
             (T(iMid-1,j,k) - T(iMid+1,j,k)) / (2.0d0 * dx)) / deltaT
      enddo
    enddo
  else
    iL = nx / 2
    iR = iL + 1
 !$acc parallel loop collapse(2) present(u,T) reduction(+:sum_mid)
    do k = 1, nz
      do j = 1, ny
        sum_mid = sum_mid + (coef * 0.5d0 * (u(iL,j,k) * (T(iL,j,k) - Tref) + &
             u(iR,j,k) * (T(iR,j,k) - Tref)) + (T(iL,j,k) - T(iR,j,k)) / dx) / deltaT
      enddo
    enddo
  endif
  call MPI_Allreduce(MPI_IN_PLACE, sum_mid, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  Nu_middle = sum_mid / dble(ny * nzGlobal)

  if (isRoot) then
    write(*,'(a,1x,es16.8)') 'Nu_hot(left)  =', Nu_hot
    write(*,'(a,1x,es16.8)') 'Nu_cold(right)=', Nu_cold
    write(*,'(a,1x,es16.8)') 'Nu_middle     =', Nu_middle
    open(unit=00, file=trim(settingsFile), status='unknown', position='append')
    write(00,'(a,1x,es16.8)') 'Nu_hot(left)  =', Nu_hot
    write(00,'(a,1x,es16.8)') 'Nu_cold(right)=', Nu_cold
    write(00,'(a,1x,es16.8)') 'Nu_middle     =', Nu_middle
    close(00)
  endif
#else
  dy = 1.0d0 / lengthUnit
  deltaT = Thot - Tcold
  coef = heatFluxScale

  sum_hot = 0.0d0
 !$acc parallel loop collapse(2) present(T) reduction(+:sum_hot)
  do k = 1, nz
    do i = 1, nx
      q_wall = 2.0d0 * (Thot - T(i,1,k)) / dy
      sum_hot = sum_hot + q_wall / deltaT
    enddo
  enddo
  call MPI_Allreduce(MPI_IN_PLACE, sum_hot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  Nu_hot = sum_hot / dble(nx * nzGlobal)

  sum_cold = 0.0d0
 !$acc parallel loop collapse(2) present(T) reduction(+:sum_cold)
  do k = 1, nz
    do i = 1, nx
      q_wall = 2.0d0 * (T(i,ny,k) - Tcold) / dy
      sum_cold = sum_cold + q_wall / deltaT
    enddo
  enddo
  call MPI_Allreduce(MPI_IN_PLACE, sum_cold, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  Nu_cold = sum_cold / dble(nx * nzGlobal)

  sum_mid = 0.0d0
  if (mod(ny,2) .EQ. 1) then
    jB = (ny + 1) / 2
 !$acc parallel loop collapse(2) present(v,T) reduction(+:sum_mid)
    do k = 1, nz
      do i = 1, nx
        sum_mid = sum_mid + (coef * v(i,jB,k) * (T(i,jB,k) - Tref) - &
             (T(i,jB+1,k) - T(i,jB-1,k)) / (2.0d0 * dy)) / deltaT
      enddo
    enddo
  else
    jB = ny / 2
    jT = jB + 1
 !$acc parallel loop collapse(2) present(v,T) reduction(+:sum_mid)
    do k = 1, nz
      do i = 1, nx
        sum_mid = sum_mid + (coef * 0.5d0 * (v(i,jB,k) * (T(i,jB,k) - Tref) + &
             v(i,jT,k) * (T(i,jT,k) - Tref)) + (T(i,jB,k) - T(i,jT,k)) / dy) / deltaT
      enddo
    enddo
  endif
  call MPI_Allreduce(MPI_IN_PLACE, sum_mid, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  Nu_middle = sum_mid / dble(nx * nzGlobal)

  if (isRoot) then
    write(*,'(a,1x,es16.8)') 'Nu_hot(bottom) =', Nu_hot
    write(*,'(a,1x,es16.8)') 'Nu_cold(top)   =', Nu_cold
    write(*,'(a,1x,es16.8)') 'Nu_middle      =', Nu_middle
    open(unit=00, file=trim(settingsFile), status='unknown', position='append')
    write(00,'(a,1x,es16.8)') 'Nu_hot(bottom) =', Nu_hot
    write(00,'(a,1x,es16.8)') 'Nu_cold(top)   =', Nu_cold
    write(00,'(a,1x,es16.8)') 'Nu_middle      =', Nu_middle
    close(00)
  endif
#endif

end subroutine RBcalc_Nu_wall_avg3d


subroutine SideHeatedcalc_umid_max3d()
  call calc_umid_max_common3d('SideHeatedcalc_umid_max')
end subroutine SideHeatedcalc_umid_max3d

subroutine SideHeatedcalc_vmid_max3d()
  call calc_vmid_max_common3d('SideHeatedcalc_vmid_max')
end subroutine SideHeatedcalc_vmid_max3d

subroutine SideHeatedcalc_wmid_max3d()
  call calc_wmid_max_common3d('SideHeatedcalc_wmid_max')
end subroutine SideHeatedcalc_wmid_max3d

subroutine RBcalc_umid_max3d()
  call calc_umid_max_common3d('RBcalc_umid_max')
end subroutine RBcalc_umid_max3d

subroutine RBcalc_vmid_max3d()
  call calc_vmid_max_common3d('RBcalc_vmid_max')
end subroutine RBcalc_vmid_max3d

subroutine RBcalc_wmid_max3d()
  call calc_wmid_max_common3d('RBcalc_wmid_max')
end subroutine RBcalc_wmid_max3d

subroutine calc_umid_max_common3d(logTag)
  use mpi
  use commondata3dOpenaccMpi
  implicit none
  character(len=*), intent(in) :: logTag
  integer(kind=4) :: j, k, iL, iR, jBest, kBest, r, bestIdx
  real(kind=8) :: targetX, weight, val, umax, yAtU, zAtU
  real(kind=8) :: localU(3)
  real(kind=8), allocatable :: gatherU(:)

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
  localU = (/ umax, yAtU, zAtU /)

  allocate(gatherU(3*mpiSize))
  call MPI_Gather(localU, 3, MPI_DOUBLE_PRECISION, gatherU, 3, MPI_DOUBLE_PRECISION, mpiRoot, MPI_COMM_WORLD, mpiErr)

  if (isRoot) then
    bestIdx = 1
    do r = 1, mpiSize-1
      if (gatherU(3*r+1) .GT. gatherU(3*bestIdx-2)) bestIdx = r + 1
    enddo
    umax = gatherU(3*bestIdx-2)
    yAtU = gatherU(3*bestIdx-1)
    zAtU = gatherU(3*bestIdx)

    write(*,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
         'u_mid_max =', umax*velocityScaleCompare, 'at y =', yAtU, 'z =', zAtU, 'on x_mid =', targetX

    open(unit=00, file=trim(settingsFile), status='unknown', position='append')
    write(00,*) '--- ', trim(logTag), ' ---'
    write(00,*) 'x_mid =', targetX
    write(00,*) 'u_mid_max =', umax*velocityScaleCompare, ' y_pos =', yAtU, ' z_pos =', zAtU
    close(00)
  endif

  deallocate(gatherU)
end subroutine calc_umid_max_common3d

subroutine calc_vmid_max_common3d(logTag)
  use mpi
  use commondata3dOpenaccMpi
  implicit none
  character(len=*), intent(in) :: logTag
  integer(kind=4) :: i, k, jL, jR, iBest, kBest, r, bestIdx
  real(kind=8) :: targetY, weight, val, vmax, xAtV, zAtV
  real(kind=8) :: localV(3)
  real(kind=8), allocatable :: gatherV(:)

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
  localV = (/ vmax, xAtV, zAtV /)

  allocate(gatherV(3*mpiSize))
  call MPI_Gather(localV, 3, MPI_DOUBLE_PRECISION, gatherV, 3, MPI_DOUBLE_PRECISION, mpiRoot, MPI_COMM_WORLD, mpiErr)

  if (isRoot) then
    bestIdx = 1
    do r = 1, mpiSize-1
      if (gatherV(3*r+1) .GT. gatherV(3*bestIdx-2)) bestIdx = r + 1
    enddo
    vmax = gatherV(3*bestIdx-2)
    xAtV = gatherV(3*bestIdx-1)
    zAtV = gatherV(3*bestIdx)

    write(*,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
         'v_mid_max =', vmax*velocityScaleCompare, 'at x =', xAtV, 'z =', zAtV, 'on y_mid =', targetY

    open(unit=00, file=trim(settingsFile), status='unknown', position='append')
    write(00,*) '--- ', trim(logTag), ' ---'
    write(00,*) 'y_mid =', targetY
    write(00,*) 'v_mid_max =', vmax*velocityScaleCompare, ' x_pos =', xAtV, ' z_pos =', zAtV
    close(00)
  endif

  deallocate(gatherV)
end subroutine calc_vmid_max_common3d

subroutine calc_wmid_max_common3d(logTag)
  use mpi
  use commondata3dOpenaccMpi
  implicit none
  character(len=*), intent(in) :: logTag
  integer(kind=4) :: i, j, kL, kR, iBest, jBest, r, bestIdx
  real(kind=8) :: targetZ, weight, val, wmax, xAtW, yAtW
  real(kind=8) :: localW(3)
  real(kind=8), allocatable :: gatherW(:)

  targetZ = 0.5d0 * dble(nzGlobal) / lengthUnit
  if ((targetZ .GT. zp(0)) .AND. (targetZ .LE. zp(nz+1))) then
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
    localW = (/ wmax, xAtW, yAtW /)
  else
    localW = (/ -huge(1.0d0), 0.0d0, 0.0d0 /)
  endif

  allocate(gatherW(3*mpiSize))
  call MPI_Gather(localW, 3, MPI_DOUBLE_PRECISION, gatherW, 3, MPI_DOUBLE_PRECISION, mpiRoot, MPI_COMM_WORLD, mpiErr)

  if (isRoot) then
    bestIdx = 1
    do r = 1, mpiSize-1
      if (gatherW(3*r+1) .GT. gatherW(3*bestIdx-2)) bestIdx = r + 1
    enddo
    wmax = gatherW(3*bestIdx-2)
    xAtW = gatherW(3*bestIdx-1)
    yAtW = gatherW(3*bestIdx)

    write(*,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
         'w_mid_max =', wmax*velocityScaleCompare, 'at x =', xAtW, 'y =', yAtW, 'on z_mid =', targetZ

    open(unit=00, file=trim(settingsFile), status='unknown', position='append')
    write(00,*) '--- ', trim(logTag), ' ---'
    write(00,*) 'z_mid =', targetZ
    write(00,*) 'w_mid_max =', wmax*velocityScaleCompare, ' x_pos =', xAtW, ' y_pos =', yAtW
    close(00)
  endif

  deallocate(gatherW)
end subroutine calc_wmid_max_common3d


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
  use commondata3dOpenaccMpi
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
  use commondata3dOpenaccMpi
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
  use commondata3dOpenaccMpi
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
  use commondata3dOpenaccMpi
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4) :: j, k, uout, iL, iR
  real(kind=8) :: targetX, weight, valU, valV, valW, valT

  ! 输出 x=Lx/2 这张中面，便于看 y-z 平面流型
  targetX = 0.5d0 * xp(nx+1)
  call find_bracketing_index(xp, nx, targetX, iL, iR, weight)
  open(newunit=uout, file=trim(filename), status='replace', action='write', form='formatted')
  write(uout,'(A)') 'VARIABLES = "Y" "Z" "U_nd" "V_nd" "W_nd" "T"'
  write(uout,'(A,I0,A,I0,A)') 'ZONE I=', ny, ', J=', nz, ', F=POINT'
  do k = 1, nz
    do j = 1, ny
      call interp_scalar_x(iL, iR, weight, j, k, u, valU)
      call interp_scalar_x(iL, iR, weight, j, k, v, valV)
      call interp_scalar_x(iL, iR, weight, j, k, w, valW)
      call interp_scalar_x(iL, iR, weight, j, k, T, valT)
      write(uout,'(6(1X,ES24.16))') yp(j), zp(k), velocityScaleCompare*valU, &
           velocityScaleCompare*valV, velocityScaleCompare*valW, valT
    enddo
  enddo
  close(uout)

end subroutine write_midplane_x


!===========================================================================================================================
! 子程序: write_midplane_y
! 作用: 输出 y=Ly/2 中面的 Tecplot 切片文件。
!===========================================================================================================================
subroutine write_midplane_y(filename)
  use commondata3dOpenaccMpi
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4) :: i, k, uout, jL, jR
  real(kind=8) :: targetY, weight, valU, valV, valW, valT

  ! 输出 y=Ly/2 中面，最适合直接看 RB 主循环胞
  targetY = 0.5d0 * yp(ny+1)
  call find_bracketing_index(yp, ny, targetY, jL, jR, weight)
  open(newunit=uout, file=trim(filename), status='replace', action='write', form='formatted')
  write(uout,'(A)') 'VARIABLES = "X" "Z" "U_nd" "V_nd" "W_nd" "T"'
  write(uout,'(A,I0,A,I0,A)') 'ZONE I=', nx, ', J=', nz, ', F=POINT'
  do k = 1, nz
    do i = 1, nx
      call interp_scalar_y(jL, jR, weight, i, k, u, valU)
      call interp_scalar_y(jL, jR, weight, i, k, v, valV)
      call interp_scalar_y(jL, jR, weight, i, k, w, valW)
      call interp_scalar_y(jL, jR, weight, i, k, T, valT)
      write(uout,'(6(1X,ES24.16))') xp(i), zp(k), velocityScaleCompare*valU, &
           velocityScaleCompare*valV, velocityScaleCompare*valW, valT
    enddo
  enddo
  close(uout)

end subroutine write_midplane_y


!===========================================================================================================================
! 子程序: write_midplane_z
! 作用: 输出 z=Lz/2 中面的 Tecplot 切片文件。
!===========================================================================================================================
subroutine write_midplane_z(filename)
  use commondata3dOpenaccMpi
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4) :: i, j, uout, kL, kR
  real(kind=8) :: targetZ, weight, valU, valV, valW, valT

  ! 输出 z=Lz/2 中面，便于观察第三方向上的结构
  ! MPI 按 z 分块后，只让拥有该中面的 rank 写这一张切片，避免多 rank 重复写。
  targetZ = 0.5d0 * dble(nzGlobal) / lengthUnit
  if ((targetZ .LE. zp(0)) .OR. (targetZ .GT. zp(nz+1))) return

  call find_bracketing_index(zp, nz, targetZ, kL, kR, weight)
  open(newunit=uout, file=trim(filename), status='replace', action='write', form='formatted')
  write(uout,'(A)') 'VARIABLES = "X" "Y" "U_nd" "V_nd" "W_nd" "T"'
  write(uout,'(A,I0,A,I0,A)') 'ZONE I=', nx, ', J=', ny, ', F=POINT'
  do j = 1, ny
    do i = 1, nx
      call interp_scalar_z(kL, kR, weight, i, j, u, valU)
      call interp_scalar_z(kL, kR, weight, i, j, v, valV)
      call interp_scalar_z(kL, kR, weight, i, j, w, valW)
      call interp_scalar_z(kL, kR, weight, i, j, T, valT)
      write(uout,'(6(1X,ES24.16))') xp(i), yp(j), velocityScaleCompare*valU, &
           velocityScaleCompare*valV, velocityScaleCompare*valW, valT
    enddo
  enddo
  close(uout)

end subroutine write_midplane_z







