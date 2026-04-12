!=============================================================
!!!    注释区，代码描述
!!!    三维浮力驱动自然对流
!!!    D3Q19 流场 + D3Q7 温度场
!!!    参考 3DRB.F90 的主流程与后处理结构
!=============================================================

!=============================================================
!   自定义宏，一些选项的开关
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
!~~temperature B.C.~~



!算法切换
!启用 M1G 修正；注释掉则不使用 useG 相关修正
#define EnableUseG
!启用旧算法
!#define EnableLegacyThermalScheme


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

  ! 在 loadInitField=1 的前提下：
  integer(kind=4), parameter :: reloadDimensionlessTime=0
  integer(kind=4), parameter :: reloadbinFileNum=0
  !===============================================================================================

  !===============================================================================================
  ! 无量纲参数
  integer(kind=4), parameter :: nx=80, ny=80, nz=80
#ifdef SideHeatedCell
  real(kind=8), parameter :: lengthUnit=dble(nx)
#else
  real(kind=8), parameter :: lengthUnit=dble(ny)
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

  ! heatFluxScale is used in Nu/heat-flux post-processing and should stay consistent with the Nu definition.
  real(kind=8), parameter :: heatFluxScale=lengthUnit/diffusivity

  ! velocityScaleCompare is used only in velocity-related post-processing to convert lattice velocity
  ! to the nondimensional velocity scale adopted by the reference paper being compared.
  real(kind=8), parameter :: velocityScaleCompare=lengthUnit/diffusivity

  ! 浮力项参数
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
  integer(kind=4), parameter :: itc_max=20000000
  real(kind=8), parameter :: outputFrequency=100.0d0
  integer(kind=4), parameter :: dimensionlessTimeMax=int(12000.0d0/outputFrequency)
  integer(kind=4), parameter :: backupInterval=1000

  real(kind=8), parameter :: epsU=1.0d-7, epsT=1.0d-7

  integer(kind=4), parameter :: outputBinFile=0
  integer(kind=4), parameter :: outputPltFile=0

  integer(kind=4) :: binFileNum, pltFileNum
  integer(kind=4) :: dimensionlessTime
  integer(kind=4) :: outputIntervalItc, backupIntervalItc

  real(kind=8) :: NuVolAvg(0:dimensionlessTimeMax), ReVolAvg(0:dimensionlessTimeMax)

  character(len=100) :: binFolderPrefix="buoyancyCavity3DOpenaccbinFile"
  character(len=100) :: pltFolderPrefix="buoyancyCavity3DOpenaccTecplot"
  character(len=100) :: reloadFilePrefix="backupFile3DOpenacc"
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

  real(kind=8) :: Nu_global, Nu_hot, Nu_cold, Nu_middle
  real(kind=8) :: Nu_hot_max, Nu_hot_min, Nu_hot_max_position, Nu_hot_min_position

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

  ! 先写日志头，并初始化 OpenACC 设备
  open(unit=00, file=trim(settingsFile), status='replace')
  string = ctime(time())
  write(00,*) 'Start: ', string
  write(00,*) 'Starting OpenACC >>>>>>'
  call acc_init(acc_device_default)
  numAccDevices = acc_get_num_devices(acc_device_default)
  write(00,*) 'Visible OpenACC devices:', numAccDevices
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
  call output_Tecplot3d()
  call output_binary3d()
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
       real(dble(nx) * dble(ny) * dble(nz) * dble(itc) / &
       & max(timeEnd - timeStart, 1.0d-12) / 1.0d6, kind=8)
  write(00,*) 'Time (ACC) = ', real(timeEnd2 - timeStart2, kind=8), 's'
  write(00,*) 'MLUPS (ACC) = ', &
       real(dble(nx) * dble(ny) * dble(nz) * dble(itc) / &
       & max(timeEnd2 - timeStart2, 1.0d-12) / 1.0d6, kind=8)
  write(00,*) 'Nu_global =', Nu_global
  write(00,*) 'Nu_hot    =', Nu_hot
  write(00,*) 'Nu_cold   =', Nu_cold
  write(00,*) 'Nu_middle =', Nu_middle
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
  real(kind=8) :: xLen, yLen, rbInitPerturbAmp
  character(len=100) :: reloadFileName

  itc = 0
  errorU = 100.0d0
  errorT = 100.0d0

  ! 把按自由落体时间给出的输出/备份间隔换算成格子步数 itc
  outputIntervalItc = max(1, int(outputFrequency * timeUnit))
  backupIntervalItc = max(1, int(backupInterval * timeUnit))

  !-----------------------------------------------------------------------------------------------
  ! 记录各种信息在日志文件中
  !-----------------------------------------------------------------------------------------------
  open(unit=00, file=trim(settingsFile), status='unknown', position='append')

  if (outputBinFile .EQ. 1) then
    open(unit=01, file=trim(binFolderPrefix)//'-'//'readme', status='unknown')
    write(01,*) 'binFile folder exist!'
    close(01)
    write(00,*) 'Data will be stored in ', trim(binFolderPrefix)
  endif

  if (outputPltFile .EQ. 1) then
    open(unit=01, file=trim(pltFolderPrefix)//'-'//'readme', status='unknown')
    write(01,*) 'pltFile folder exist!'
    close(01)
    write(00,*) 'Data will be stored in ', trim(pltFolderPrefix)
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
  write(00,*) 'outputFrequency =', real(outputFrequency,kind=8), ' free-fall time units'
  write(00,*) '......................  or ', outputIntervalItc, ' in itc units'
  write(00,*) 'backupInterval =', backupInterval, ' free-fall time units'
  write(00,*) '.................... or ', backupIntervalItc, ' in itc units'
  if (loadInitField .EQ. 1) then
    write(00,*) 'reloadDimensionlessTime =', reloadDimensionlessTime
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
  allocate(f(0:qf-1,nx,ny,nz), f_post(0:qf-1,0:nx+1,0:ny+1,0:nz+1))
  allocate(g(0:qt-1,nx,ny,nz), g_post(0:qt-1,0:nx+1,0:ny+1,0:nz+1))
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
    write(reloadFileName,*) reloadbinFileNum
    reloadFileName = adjustl(reloadFileName)
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
  use commondata3dOpenacc
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
  use commondata3dOpenacc
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
  use commondata3dOpenacc
  implicit none

#ifdef steadyFlow
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
 !$acc parallel loop collapse(3) present(f,f_post,rho,u,v,w,T)
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
        FyLoc = rhoLoc * gBeta * (T(i,j,k) - Tref)
        FzLoc = 0.0d0
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

end subroutine collision3d


!===========================================================================================================================
! 子程序: fill_periodic_ghosts_f_post
! 作用: 为流场分布函数补齐周期 ghost 层，便于统一 pull streaming。
!===========================================================================================================================
subroutine fill_periodic_ghosts_f_post()
  use commondata3dOpenacc
  implicit none

  integer(kind=4) :: i, j, k, alpha

#ifdef VerticalWallsPeriodicalU
  ! x 方向是周期边界时，先把左右 ghost 层补齐
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

#ifdef SpanwiseWallsPeriodicalU
  ! z 方向是周期边界时，再把前后 ghost 层补齐
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

end subroutine fill_periodic_ghosts_f_post


!===========================================================================================================================
! 子程序: streaming3d
! 作用: 对流场分布函数执行三维 pull streaming。
!===========================================================================================================================
subroutine streaming3d()
  use commondata3dOpenacc
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
  use commondata3dOpenacc
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
 !$acc parallel loop collapse(2) present(f,f_post,ez)
  do j = 1, ny
    do i = 1, nx
      do alpha = 0, qf-1
        if (ez(alpha) .EQ. 1)  f(alpha,i,j,1)  = f_post(alpha,i,j,nz)
        if (ez(alpha) .EQ. -1) f(alpha,i,j,nz) = f_post(alpha,i,j,1)
      enddo
    enddo
  enddo
#endif

#ifdef SpanwiseWallsNoslip
 !$acc parallel loop collapse(2) present(f,f_post,ez,opp)
  do j = 1, ny
    do i = 1, nx
      do alpha = 0, qf-1
        if (ez(alpha) .EQ. 1)  f(alpha,i,j,1)  = f_post(opp(alpha),i,j,1)
        if (ez(alpha) .EQ. -1) f(alpha,i,j,nz) = f_post(opp(alpha),i,j,nz)
      enddo
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
  real(kind=8) :: FyLoc

 !$acc parallel loop collapse(3) present(f,rho,u,v,w,T)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        rho(i,j,k) = f(0,i,j,k) + f(1,i,j,k) + f(2,i,j,k) + f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + &
             f(6,i,j,k) + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + &
             f(12,i,j,k) + f(13,i,j,k) + f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k)

        FyLoc = rho(i,j,k) * gBeta * (T(i,j,k) - Tref)

        u(i,j,k) = ( f(1,i,j,k) - f(2,i,j,k) + f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) - f(10,i,j,k) + &
             f(11,i,j,k) - f(12,i,j,k) + f(13,i,j,k) - f(14,i,j,k) ) / rho(i,j,k)

        v(i,j,k) = ( f(3,i,j,k) - f(4,i,j,k) + f(7,i,j,k) + f(8,i,j,k) - f(9,i,j,k) - f(10,i,j,k) + &
             f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) - f(18,i,j,k) + 0.5d0 * FyLoc ) / rho(i,j,k)

        w(i,j,k) = ( f(5,i,j,k) - f(6,i,j,k) + f(11,i,j,k) + f(12,i,j,k) - f(13,i,j,k) - f(14,i,j,k) + &
             f(15,i,j,k) + f(16,i,j,k) - f(17,i,j,k) - f(18,i,j,k) ) / rho(i,j,k)
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

end subroutine collisionT3d


!===========================================================================================================================
! 子程序: fill_periodic_ghosts_g_post
! 作用: 为温度分布函数补齐周期 ghost 层，便于统一 pull streaming。
!===========================================================================================================================
subroutine fill_periodic_ghosts_g_post()
  use commondata3dOpenacc
  implicit none

  integer(kind=4) :: i, j, k, alpha

#ifdef VerticalWallsPeriodicalT
  ! x 方向是周期边界时，先把左右 ghost 层补齐
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

#ifdef SpanwiseWallsPeriodicalT
  ! z 方向是周期边界时，再把前后 ghost 层补齐
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

end subroutine fill_periodic_ghosts_g_post


!===========================================================================================================================
! 子程序: streamingT3d
! 作用: 对温度分布函数执行三维 pull streaming。
!===========================================================================================================================
subroutine streamingT3d()
  use commondata3dOpenacc
  implicit none

  integer(kind=4) :: i, j, k, ip, jp, kp, alpha

  ! 温度场同样采用 pull streaming
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

  return
end subroutine streamingT3d


!===========================================================================================================================
! 子程序: bouncebackT3d
! 作用: 施加温度边界条件，包括恒温、绝热和周期边界。
!===========================================================================================================================
subroutine bouncebackT3d()
  use commondata3dOpenacc
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
 !$acc parallel loop collapse(2) present(g,g_post,ezT)
  do j = 1, ny
    do i = 1, nx
      do alpha = 0, qt-1
        if (ezT(alpha) .EQ. 1)  g(alpha,i,j,1)  = g_post(alpha,i,j,nz)
        if (ezT(alpha) .EQ. -1) g(alpha,i,j,nz) = g_post(alpha,i,j,1)
      enddo
    enddo
  enddo
#endif

#ifdef SpanwiseWallsAdiabatic
 !$acc parallel loop collapse(2) present(g,g_post,ezT,oppT)
  do j = 1, ny
    do i = 1, nx
      do alpha = 0, qt-1
        if (ezT(alpha) .EQ. 1)  g(alpha,i,j,1)  = g_post(oppT(alpha),i,j,1)
        if (ezT(alpha) .EQ. -1) g(alpha,i,j,nz) = g_post(oppT(alpha),i,j,nz)
      enddo
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
  use commondata3dOpenacc
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
          ! 力项里含有 rho 和 T，因此这里做几次简单迭代来恢复带半步修正的速度
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

  call append_convergence_tecplot3d('convergence3D.plt', itc, errorU, errorT)
  write(caseTag,'("Ra=",ES10.3E2,",nx=",I0,",ny=",I0,",nz=",I0,",useG=",L1,",old=",L1)') &
       Rayleigh, nx, ny, nz, useG, useLegacyThermalScheme
  call append_convergence_master_tecplot3d('convergence_all_3D.plt', caseTag, itc, errorU, errorT)
  write(*,'(I12,1X,ES24.16,1X,ES24.16)') itc, errorU, errorT

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
  use commondata3dOpenacc
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
  open(unit=03, file=trim(binFolderPrefix)//'-'//trim(filename)//'.bin', form='unformatted', access='sequential')
  write(03) (((velocityScaleCompare*u(i,j,k), i=1,nx), j=1,ny), k=1,nz)
  write(03) (((velocityScaleCompare*v(i,j,k), i=1,nx), j=1,ny), k=1,nz)
  write(03) (((velocityScaleCompare*w(i,j,k), i=1,nx), j=1,ny), k=1,nz)
  write(03) (((T(i,j,k), i=1,nx), j=1,ny), k=1,nz)
  write(03) (((rho(i,j,k), i=1,nx), j=1,ny), k=1,nz)
  close(03)

  return
end subroutine output_binary3d


!===========================================================================================================================
! 子程序: backupData3d
! 作用: 备份 f/g 分布函数，供后续重启继续计算。
!===========================================================================================================================
subroutine backupData3d()
  use commondata3dOpenacc
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

  return
end subroutine backupData3d


!===========================================================================================================================
! 子程序: output_Tecplot3d
! 作用: 输出全场体数据以及 x/y/z 三个中面切片，便于快速查看三维流场结构。
!===========================================================================================================================
subroutine output_Tecplot3d()
  use commondata3dOpenacc
  implicit none

  character(len=100) :: filename

#ifdef steadyFlow
  write(filename,'(i12.12)') itc
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

  ! 这里记录的是时间序列版本的体平均 Nu / Re：
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

  NuVolAvg(dimensionlessTime) = NuVolAvg_temp / dble(nx * ny * nz) * lengthUnit / diffusivity + 1.0d0
  open(unit=01, file='Nu_VolAvg_3D.dat', status='unknown', position='append')
  write(01,*) real(reloadDimensionlessTime + dimensionlessTime * outputFrequency, kind=8), NuVolAvg(dimensionlessTime)
  close(01)

  ReVolAvg_temp = 0.0d0
 !$acc parallel loop collapse(3) present(u,v,w) reduction(+:ReVolAvg_temp)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        ReVolAvg_temp = ReVolAvg_temp + u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)
      enddo
    enddo
  enddo
  ReVolAvg(dimensionlessTime) = dsqrt(ReVolAvg_temp / dble(nx * ny * nz)) * lengthUnit / viscosity
  open(unit=02, file='Re_VolAvg_3D.dat', status='unknown', position='append')
  write(02,*) real(reloadDimensionlessTime + dimensionlessTime * outputFrequency, kind=8), ReVolAvg(dimensionlessTime)
  close(02)

  write(*,'(a,1x,es16.8)') 'NuVolAvg =', NuVolAvg(dimensionlessTime)
  write(*,'(a,1x,es16.8)') 'ReVolAvg =', ReVolAvg(dimensionlessTime)

  return
end subroutine calNuRe3d


!===========================================================================================================================
! 子程序: RBcalc_Nu_global3d
! 作用: 计算三维算例的全局平均 Nusselt 数。
!===========================================================================================================================
subroutine RBcalc_Nu_global3d()
  use commondata3dOpenacc
  implicit none
  integer(kind=4) :: i, j, k
  real(kind=8) :: dy, dTdy, qy, sum_qy
  real(kind=8) :: deltaT, coef

  dy = 1.0d0 / lengthUnit
  deltaT = Thot - Tcold
  coef = velocityScaleCompare
  sum_qy = 0.0d0

 !$acc parallel loop collapse(3) present(v,T) reduction(+:sum_qy)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        if (j .EQ. 1) then
          dTdy = (3.0d0*T(i,1,k) + T(i,2,k) - 4.0d0*Thot) / (3.0d0*dy)
        elseif (j .EQ. ny) then
          dTdy = (4.0d0*Tcold - 3.0d0*T(i,ny,k) - T(i,ny-1,k)) / (3.0d0*dy)
        else
          dTdy = (T(i,j+1,k) - T(i,j-1,k)) / (2.0d0*dy)
        endif

        qy = coef * v(i,j,k) * (T(i,j,k) - Tref) - dTdy
        sum_qy = sum_qy + qy
      enddo
    enddo
  enddo

  Nu_global = (sum_qy / dble(nx * ny * nz)) / deltaT

  write(*,'(a,1x,es16.8)') 'Nu_global =', Nu_global
  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,'(a,1x,es16.8)') 'Nu_global =', Nu_global
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
  integer(kind=4) :: i, imax, imin, jB, jT, jMid, k, m
  integer(kind=4) :: ii(5)
  real(kind=8) :: dx, dy, deltaT, coef
  real(kind=8) :: qy_wall, sum_hot, sum_cold, sum_mid
  real(kind=8) :: T_wl, T_wr
  real(kind=8) :: xfit(4), Tfit(4)
  real(kind=8) :: xk(5), fk(5), fstar, xstar
  real(kind=8) :: Nu_bot(1:nx), Nu_bot_ext(0:nx+1), T_bot_avg(1:nx)

  dx = 1.0d0 / lengthUnit
  dy = 1.0d0 / lengthUnit
  deltaT = Thot - Tcold
  coef = heatFluxScale

 !$acc parallel loop present(T,T_bot_avg)
  do i = 1, nx
    T_bot_avg(i) = 0.0d0
    do k = 1, nz
      T_bot_avg(i) = T_bot_avg(i) + T(i,1,k)
    enddo
    T_bot_avg(i) = T_bot_avg(i) / dble(nz)
  enddo

  sum_hot = 0.0d0
 !$acc parallel loop present(T,Nu_bot) reduction(+:sum_hot)
  do i = 1, nx
    Nu_bot(i) = 0.0d0
    do k = 1, nz
      qy_wall = (8.0d0*Thot - 9.0d0*T(i,1,k) + T(i,2,k)) / (3.0d0*dy)
      Nu_bot(i) = Nu_bot(i) + qy_wall / deltaT
    enddo
    Nu_bot(i) = Nu_bot(i) / dble(nz)
    sum_hot = sum_hot + Nu_bot(i)
  enddo
  Nu_hot = sum_hot / dble(nx)

  Nu_bot_ext(1:nx) = Nu_bot(1:nx)
  xfit(1) = xp(1);  Tfit(1) = T_bot_avg(1)
  xfit(2) = xp(2);  Tfit(2) = T_bot_avg(2)
  xfit(3) = xp(3);  Tfit(3) = T_bot_avg(3)
  xfit(4) = xp(4);  Tfit(4) = T_bot_avg(4)
  call fit_adiabatic_wall_T4_3d(0.0d0, xfit, Tfit, T_wl)
  Nu_bot_ext(0) = (2.0d0 * (Thot - T_wl) / dy) / deltaT

  xfit(1) = xp(nx-3);  Tfit(1) = T_bot_avg(nx-3)
  xfit(2) = xp(nx-2);  Tfit(2) = T_bot_avg(nx-2)
  xfit(3) = xp(nx-1);  Tfit(3) = T_bot_avg(nx-1)
  xfit(4) = xp(nx  );  Tfit(4) = T_bot_avg(nx  )
  call fit_adiabatic_wall_T4_3d(xp(nx+1), xfit, Tfit, T_wr)
  Nu_bot_ext(nx+1) = (2.0d0 * (Thot - T_wr) / dy) / deltaT

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
 !$acc parallel loop collapse(2) present(T) reduction(+:sum_cold)
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
 !$acc parallel loop collapse(2) present(v,T) reduction(+:sum_mid)
    do k = 1, nz
      do i = 1, nx
        sum_mid = sum_mid + (coef * v(i,jMid,k) * (T(i,jMid,k) - Tref) - &
             (T(i,jMid+1,k) - T(i,jMid-1,k)) / (2.0d0 * dy)) / deltaT
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
  Nu_middle = sum_mid / dble(nx * nz)

  write(*,'(a,1x,es16.8)') 'Nu_hot(bottom) =', Nu_hot
  write(*,'(a,1x,es16.8)') 'Nu_cold(top)   =', Nu_cold
  write(*,'(a,1x,es16.8)') 'Nu_middle      =', Nu_middle
  write(*,'(a,1x,es16.8,2x,a,1x,es16.8)') 'Nu_hot_max =', Nu_hot_max, 'x_max =', Nu_hot_max_position
  write(*,'(a,1x,es16.8,2x,a,1x,es16.8)') 'Nu_hot_min =', Nu_hot_min, 'x_min =', Nu_hot_min_position
  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,'(a,1x,es16.8)') 'Nu_hot(bottom) =', Nu_hot
  write(00,'(a,1x,es16.8)') 'Nu_cold(top)   =', Nu_cold
  write(00,'(a,1x,es16.8)') 'Nu_middle      =', Nu_middle
  write(00,'(a,1x,es16.8,2x,a,1x,es16.8)') 'Nu_hot_max =', Nu_hot_max, 'x_max =', Nu_hot_max_position
  write(00,'(a,1x,es16.8,2x,a,1x,es16.8)') 'Nu_hot_min =', Nu_hot_min, 'x_min =', Nu_hot_min_position
  close(00)

  return
end subroutine RBcalc_Nu_wall_avg3d


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
! 子程序: SideHeatedcalc_Nu_global3d
! 作用: 计算侧壁差温工况下的全场平均 Nusselt 数。
!===========================================================================================================================
subroutine SideHeatedcalc_Nu_global3d()
  use commondata3dOpenacc
  implicit none
  integer(kind=4) :: i, j, k
  real(kind=8) :: dx, dTdx, qx, sum_qx
  real(kind=8) :: deltaT, coef

  dx = 1.0d0 / lengthUnit
  deltaT = Thot - Tcold
  coef = velocityScaleCompare
  sum_qx = 0.0d0

  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        if (i .EQ. 1) then
          dTdx = (-3.0d0*T(1,j,k) - T(2,j,k) + 4.0d0*Thot) / (3.0d0*dx)
        elseif (i .EQ. nx) then
          dTdx = (-4.0d0*Tcold + 3.0d0*T(nx,j,k) + T(nx-1,j,k)) / (3.0d0*dx)
        else
          dTdx = (T(i-1,j,k) - T(i+1,j,k)) / (2.0d0*dx)
        endif

        qx = coef * u(i,j,k) * (T(i,j,k) - Tref) + dTdx
        sum_qx = sum_qx + qx
      enddo
    enddo
  enddo

  Nu_global = (sum_qx / dble(nx * ny * nz)) / deltaT

  write(*,'(a,1x,es16.8)') 'Nu_global =', Nu_global
  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,'(a,1x,es16.8)') 'Nu_global =', Nu_global
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
  integer(kind=4) :: iL, iR, iMid, j, jmax, jmin, k, m
  integer(kind=4) :: jj(5)
  real(kind=8) :: dx, deltaT, coef
  real(kind=8) :: qx_wall, sum_hot, sum_cold, sum_mid
  real(kind=8) :: T_wb, T_wt
  real(kind=8) :: yfit(4), Tfit(4)
  real(kind=8) :: yk(5), fk(5), fstar, ystar
  real(kind=8) :: Nu_left(1:ny), Nu_left_ext(0:ny+1), T_left_avg(1:ny)

  dx = 1.0d0 / lengthUnit
  deltaT = Thot - Tcold
  coef = heatFluxScale

  do j = 1, ny
    T_left_avg(j) = 0.0d0
    do k = 1, nz
      T_left_avg(j) = T_left_avg(j) + T(1,j,k)
    enddo
    T_left_avg(j) = T_left_avg(j) / dble(nz)
  enddo

  sum_hot = 0.0d0
  do j = 1, ny
    Nu_left(j) = 0.0d0
    do k = 1, nz
      qx_wall = (8.0d0*Thot - 9.0d0*T(1,j,k) + T(2,j,k)) / (3.0d0*dx)
      Nu_left(j) = Nu_left(j) + qx_wall / deltaT
    enddo
    Nu_left(j) = Nu_left(j) / dble(nz)
    sum_hot = sum_hot + Nu_left(j)
  enddo
  Nu_hot = sum_hot / dble(ny)

  Nu_left_ext(1:ny) = Nu_left(1:ny)
  yfit(1) = yp(1);  Tfit(1) = T_left_avg(1)
  yfit(2) = yp(2);  Tfit(2) = T_left_avg(2)
  yfit(3) = yp(3);  Tfit(3) = T_left_avg(3)
  yfit(4) = yp(4);  Tfit(4) = T_left_avg(4)
  call fit_adiabatic_wall_T4_3d(0.0d0, yfit, Tfit, T_wb)
  Nu_left_ext(0) = (2.0d0 * (Thot - T_wb) / dx) / deltaT

  yfit(1) = yp(ny-3);  Tfit(1) = T_left_avg(ny-3)
  yfit(2) = yp(ny-2);  Tfit(2) = T_left_avg(ny-2)
  yfit(3) = yp(ny-1);  Tfit(3) = T_left_avg(ny-1)
  yfit(4) = yp(ny  );  Tfit(4) = T_left_avg(ny  )
  call fit_adiabatic_wall_T4_3d(yp(ny+1), yfit, Tfit, T_wt)
  Nu_left_ext(ny+1) = (2.0d0 * (Thot - T_wt) / dx) / deltaT

  jmax = 0
  jmin = 0
  Nu_hot_max = Nu_left_ext(0)
  Nu_hot_min = Nu_left_ext(0)
  do j = 1, ny+1
    if (Nu_left_ext(j) .GT. Nu_hot_max) then
      Nu_hot_max = Nu_left_ext(j)
      jmax = j
    endif
    if (Nu_left_ext(j) .LT. Nu_hot_min) then
      Nu_hot_min = Nu_left_ext(j)
      jmin = j
    endif
  enddo

  if (jmax .LE. 2) then
    jj = (/ 0, 1, 2, 3, 4 /)
  elseif (jmax .GE. ny-1) then
    jj = (/ ny-3, ny-2, ny-1, ny, ny+1 /)
  else
    jj = (/ jmax-2, jmax-1, jmax, jmax+1, jmax+2 /)
  endif
  do m = 1, 5
    yk(m) = yp(jj(m))
    fk(m) = Nu_left_ext(jj(m))
  enddo
  call fit_parabola_ls5_3d(yk, fk, +1, fstar, ystar)
  Nu_hot_max = fstar
  Nu_hot_max_position = ystar

  if (jmin .GE. 4) then
    jj = (/ jmin-4, jmin-3, jmin-2, jmin-1, jmin /)
  else
    jj = (/ 0, 1, 2, 3, 4 /)
  endif
  do m = 1, 5
    yk(m) = yp(jj(m))
    fk(m) = Nu_left_ext(jj(m))
  enddo
  call fit_parabola_ls5_3d(yk, fk, -1, fstar, ystar)
  Nu_hot_min = fstar
  Nu_hot_min_position = ystar

  sum_cold = 0.0d0
  do j = 1, ny
    do k = 1, nz
      qx_wall = (-8.0d0*Tcold + 9.0d0*T(nx,j,k) - T(nx-1,j,k)) / (3.0d0*dx)
      sum_cold = sum_cold + qx_wall / deltaT
    enddo
  enddo
  Nu_cold = sum_cold / dble(ny * nz)

  sum_mid = 0.0d0
  if (mod(nx,2) .EQ. 1) then
    iMid = (nx + 1) / 2
    do k = 1, nz
      do j = 1, ny
        sum_mid = sum_mid + (coef * u(iMid,j,k) * (T(iMid,j,k) - Tref) + &
             (T(iMid-1,j,k) - T(iMid+1,j,k)) / (2.0d0*dx)) / deltaT
      enddo
    enddo
  else
    iL = nx / 2
    iR = iL + 1
    do k = 1, nz
      do j = 1, ny
        sum_mid = sum_mid + (coef * 0.5d0 * (u(iL,j,k) * (T(iL,j,k) - Tref) + &
             u(iR,j,k) * (T(iR,j,k) - Tref)) + (T(iL,j,k) - T(iR,j,k)) / dx) / deltaT
      enddo
    enddo
  endif
  Nu_middle = sum_mid / dble(ny * nz)

  write(*,'(a,1x,es16.8)') 'Nu_hot(left)  =', Nu_hot
  write(*,'(a,1x,es16.8)') 'Nu_cold(right)=', Nu_cold
  write(*,'(a,1x,es16.8)') 'Nu_middle     =', Nu_middle
  write(*,'(a,1x,es16.8,2x,a,1x,es16.8)') 'Nu_hot_max =', Nu_hot_max, 'y_max =', Nu_hot_max_position
  write(*,'(a,1x,es16.8,2x,a,1x,es16.8)') 'Nu_hot_min =', Nu_hot_min, 'y_min =', Nu_hot_min_position

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,'(a,1x,es16.8)') 'Nu_hot(left)  =', Nu_hot
  write(00,'(a,1x,es16.8)') 'Nu_cold(right)=', Nu_cold
  write(00,'(a,1x,es16.8)') 'Nu_middle     =', Nu_middle
  write(00,'(a,1x,es16.8,2x,a,1x,es16.8)') 'Nu_hot_max =', Nu_hot_max, 'y_max =', Nu_hot_max_position
  write(00,'(a,1x,es16.8,2x,a,1x,es16.8)') 'Nu_hot_min =', Nu_hot_min, 'y_min =', Nu_hot_min_position
  close(00)

  return
end subroutine SideHeatedcalc_Nu_wall_avg3d


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
!===========================================================================================================================
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


!===========================================================================================================================
! 子程序: calc_umid_max_common3d
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

  write(*,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'u_mid_max =', umax*velocityScaleCompare, 'at y =', yAtU, 'z =', zAtU, 'on x_mid =', targetX

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) '--- ', trim(logTag), ' ---'
  write(00,*) 'x_mid =', targetX
  write(00,*) 'u_mid_max =', umax*velocityScaleCompare, ' y_pos =', yAtU, ' z_pos =', zAtU
  close(00)

end subroutine calc_umid_max_common3d


!===========================================================================================================================
! 子程序: calc_vmid_max_common3d
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

  write(*,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'v_mid_max =', vmax*velocityScaleCompare, 'at x =', xAtV, 'z =', zAtV, 'on y_mid =', targetY

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) '--- ', trim(logTag), ' ---'
  write(00,*) 'y_mid =', targetY
  write(00,*) 'v_mid_max =', vmax*velocityScaleCompare, ' x_pos =', xAtV, ' z_pos =', zAtV
  close(00)

end subroutine calc_vmid_max_common3d


!===========================================================================================================================
! 子程序: calc_wmid_max_common3d
!===========================================================================================================================
subroutine calc_wmid_max_common3d(logTag)
  use commondata3dOpenacc
  implicit none
  character(len=*), intent(in) :: logTag
  integer(kind=4) :: i, j, kL, kR, iBest, jBest
  real(kind=8) :: targetZ, weight, val, wmax, xAtW, yAtW

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

  write(*,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'w_mid_max =', wmax*velocityScaleCompare, 'at x =', xAtW, 'y =', yAtW, 'on z_mid =', targetZ

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



