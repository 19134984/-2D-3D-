!=============================================================
!!!    注释区，代码描述
!!!    三维浮力驱动自然对流
!!!    D3Q19 流场 + D3Q7 温度场
!=============================================================

!=============================================================
!自定义宏，一些选项的开关
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
module commondata3d
  implicit none
  !===============================================================================================
  ! 格子离散速度数
  integer(kind=4), parameter :: qf=19, qt=7
  !===============================================================================================

  !===============================================================================================
  ! 是否在计算前从旧算例重启
  integer(kind=4), parameter :: loadInitField=0   ! 0: 不重启；1: 从 backupFile3D-*.bin 读取初值

  ! 在 loadInitField=1 的前提下：
  integer(kind=4), parameter :: reloadDimensionlessTime=0  ! 旧算例累计的无量纲时间
  integer(kind=4), parameter :: reloadbinFileNum=0         ! 读取的备份文件编号：backupFile3D-<reloadbinFileNum>.bin
  !===============================================================================================

  !===============================================================================================
  ! 无量纲参数
  integer(kind=4), parameter :: nx=40, ny=40, nz=40
#ifdef SideHeatedCell
  real(kind=8), parameter :: lengthUnit=dble(nx)     ! 侧壁差温：特征长度取 x 方向长度
#else
  real(kind=8), parameter :: lengthUnit=dble(ny)     ! 上下差温：特征长度取 y 方向长度
#endif
  real(kind=8), parameter :: pi=acos(-1.0d0)

  real(kind=8), parameter :: Rayleigh=2.0d3
  real(kind=8), parameter :: Prandtl=0.71d0
  real(kind=8), parameter :: Mach=0.1d0
  real(kind=8), parameter :: Thot=0.5d0, Tcold=-0.5d0
  real(kind=8), parameter :: Tref=0.5d0*(Thot+Tcold)
  real(kind=8), parameter :: tauf=0.5d0+Mach*lengthUnit*dsqrt(3.0d0*Prandtl/Rayleigh)
  real(kind=8), parameter :: viscosity=(tauf-0.5d0)/3.0d0
  real(kind=8), parameter :: diffusivity=viscosity/Prandtl

  real(kind=8), parameter :: cs2T=0.25d0

  ! 高阶矩参数修正aT
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

  !温度方程的多松弛系数
#ifdef EnableLegacyThermalScheme
  real(kind=8), parameter :: Qk=3.0d0-dsqrt(3.0d0), Qnu=4.0d0*dsqrt(3.0d0)-6.0d0
  real(kind=8), parameter :: thermalGeqCoeff=21.0d0/(6.0d0+paraA)
#else
  real(kind=8), parameter :: taug=0.5d0+diffusivity/cs2T
  real(kind=8), parameter :: Qnu=1.0d0, Qk=1.0d0/taug
  real(kind=8), parameter :: thermalGeqCoeff=1.0d0/cs2T
#endif
  !===============================================================================================

  !===============================================================================================
  ! 输出/备份相关设置（以自由落体时间 t_ff 为单位）
  integer(kind=4), parameter :: itc_max=20000000
  real(kind=8), parameter :: outputFrequency=100.0d0   ! 每隔 outputFrequency 个自由落体时间输出/统计一次
  integer(kind=4), parameter :: dimensionlessTimeMax=max(1,int(12000.0d0/outputFrequency))
  integer(kind=4), parameter :: backupInterval=1000    ! 备份间隔（自由落体时间单位）

  real(kind=8), parameter :: epsU=1.0d-7, epsT=1.0d-7

  integer(kind=4), parameter :: outputBinFile=0
  integer(kind=4), parameter :: outputPltFile=0

  integer(kind=4) :: binFileNum, pltFileNum
  integer(kind=4) :: dimensionlessTime
  integer(kind=4) :: outputIntervalItc, backupIntervalItc

  real(kind=8) :: NuVolAvg(0:dimensionlessTimeMax), ReVolAvg(0:dimensionlessTimeMax)
  ! 体平均 Nu 和 Re 的时间序列缓存

  character(len=100) :: binFolderPrefix="buoyancyCavity3DbinFile"
  character(len=100) :: pltFolderPrefix="buoyancyCavity3DTecplot"
  character(len=100) :: reloadFilePrefix="backupFile3D"
  character(len=100) :: settingsFile="SimulationSettings3D.txt"
  !===============================================================================================

  !===============================================================================================
  ! 计算中需要的相关参数
  real(kind=8) :: errorU, errorT

  real(kind=8) :: xp(0:nx+1), yp(0:ny+1), zp(0:nz+1)   ! 无量纲坐标数组，包括边界
  real(kind=8), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:), T(:,:,:), rho(:,:,:)

#ifdef steadyFlow
  real(kind=8), allocatable :: up(:,:,:), vp(:,:,:), wp(:,:,:), Tp(:,:,:)
#endif
  real(kind=8), allocatable :: f(:,:,:,:), f_post(:,:,:,:)
  real(kind=8), allocatable :: g(:,:,:,:), g_post(:,:,:,:)
  real(kind=8), allocatable :: Fx(:,:,:), Fy(:,:,:), Fz(:,:,:)
  real(kind=8), allocatable :: Bx_prev(:,:,:), By_prev(:,:,:), Bz_prev(:,:,:)

  integer(kind=4) :: itc
#ifdef EnableUseG
  logical, parameter :: useG=.true.            
#else
  logical, parameter :: useG=.false.           
#endif

#ifdef EnableLegacyThermalScheme
  logical, parameter :: useLegacyThermalScheme=.true.
#else
  logical, parameter :: useLegacyThermalScheme=.false.
#endif

  real(kind=8) :: Nu_global, Nu_hot, Nu_cold, Nu_middle
  real(kind=8) :: Nu_hot_max, Nu_hot_min, Nu_hot_max_position, Nu_hot_min_position

  ! 格子离散速度、反向索引和权重
  integer(kind=4) :: ex(0:qf-1), ey(0:qf-1), ez(0:qf-1), opp(0:qf-1)
  integer(kind=4) :: exT(0:qt-1), eyT(0:qt-1), ezT(0:qt-1), oppT(0:qt-1)
  real(kind=8) :: omega(0:qf-1), omegaT(0:qt-1)

  ! MRT 变换矩阵及逆矩阵
  real(kind=8) :: M19(qf,qf), Minv19(qf,qf)
  real(kind=8) :: M7(qt,qt), Minv7(qt,qt)
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


  open(unit=00, file=trim(settingsFile), status='unknown')
  string = ctime(time())
  write(00,*) 'Start: ', string
  write(00,*) 'Starting OpenMP >>>>>>'
  call OMP_set_num_threads(24)
  myMaxThreads = OMP_get_max_threads()
  write(00,*) 'Max Running threads:', myMaxThreads
  close(00)

  call initial()

  call CPU_TIME(timeStart)
  timeStart2 = OMP_get_wtime()

  do while( ((errorU.GT.epsU).OR.(errorT.GT.epsT)).AND.(itc.LE.itc_max) )
  
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
    if(MOD(itc,2000).EQ.0) call check()
#endif

    if( MOD(itc, int(outputFrequency*timeUnit)).EQ.0 ) then
      call calNuRe()

#ifdef steadyFlow
      if( (outputPltFile.EQ.1).AND.(MOD(itc, backupInterval*int(timeUnit)).EQ.0) ) call output_Tecplot()
#endif

#ifdef unsteadyFlow
      if(outputBinFile.EQ.1) then
        call output_binary()
        if(MOD(itc, int(backupInterval/outputFrequency)*int(outputFrequency*timeUnit)).EQ.0) call backupData()
      endif
      if(outputPltFile.EQ.1) call output_Tecplot()
#endif
    endif
  enddo

  call CPU_TIME(timeEnd)
  timeEnd2 = OMP_get_wtime()

#ifdef steadyFlow
  call output_Tecplot()
  call output_binary()
#endif

#ifdef SideHeatedCell
  call SideHeatedcalc_Nu_global()
  call SideHeatedcalc_Nu_wall_avg()
  call SideHeatedcalc_umid_max()
  call SideHeatedcalc_vmid_max()
  call SideHeatedcalc_wmid_max()
#endif

#ifdef RayleighBenardCell
  call RBcalc_Nu_global()
  call RBcalc_Nu_wall_avg()
  call RBcalc_umid_max()
  call RBcalc_vmid_max()
  call RBcalc_wmid_max()
#endif

  call calNuRe()

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
  write(00,*) 'Nu_global =', Nu_global
  write(00,*) 'Nu_hot    =', Nu_hot
  write(00,*) 'Nu_cold   =', Nu_cold
  write(00,*) 'Nu_middle =', Nu_middle
  write(00,*) 'Deallocate Array......'
  close(00)

  deallocate(f, f_post, g, g_post)
  deallocate(u, v, w, T, rho)
#ifdef steadyFlow
  deallocate(up, vp, wp, Tp)
#endif
  deallocate(Fx, Fy, Fz, Bx_prev, By_prev, Bz_prev)

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
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
  real(kind=8) :: feq(qf), geq(qt)
  real(kind=8) :: xLen, yLen, rbInitPerturbAmp
  character(len=100) :: reloadFileName

  ! 初始时把误差设大，这样主循环一定能先进入一轮
  itc = 0
  errorU = 100.0d0
  errorT = 100.0d0

  ! 把按自由落体时间给出的输出/备份间隔换算成格子步数 itc
  outputIntervalItc = max(1, int(outputFrequency * timeUnit))
  backupIntervalItc = max(1, int(backupInterval * timeUnit))

  !-----------------------------------------------------------------------------------------------
  ! 节点坐标数组与内存分配
  !-----------------------------------------------------------------------------------------------

  ! 网格点采用 cell-center 布置，内点坐标取 i-0.5、j-0.5、k-0.5，再统一除以 lengthUnit 做无量纲
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
  allocate(Fx(nx,ny,nz), Fy(nx,ny,nz), Fz(nx,ny,nz))
  allocate(Bx_prev(nx,ny,nz), By_prev(nx,ny,nz), Bz_prev(nx,ny,nz))

  !-----------------------------------------------------------------------------------------------
  ! 初始化离散速度、权重，以及 MRT 变换矩阵
  !-----------------------------------------------------------------------------------------------
  call init_lattice_constants()
  call init_mrt_matrices()

#ifdef EnableLegacyThermalScheme
  if ((paraA .GE. 1.0d0) .OR. (paraA .LE. -6.0d0)) then
    open(unit=00, file=trim(settingsFile), status='unknown', position='append')
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

  !-----------------------------------------------------------------------------------------------
  ! 记录各种信息在日志文件中
  !-----------------------------------------------------------------------------------------------
  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) '-------------------------------------------------------------------------------'
#ifdef SideHeatedCell
  write(00,*) 'I am 3D Side Heated Closed Cavity'
#else
  write(00,*) 'I am 3D periodic Rayleigh-Benard Cell'
#endif
  write(00,*) 'Mesh:', nx, ny, nz
  write(00,*) 'Rayleigh=', real(Rayleigh,kind=8), '; Prandtl =', real(Prandtl,kind=8), '; Mach =', real(Mach,kind=8)
  write(00,*) 'Length unit: L0 =', real(lengthUnit,kind=8)
  write(00,*) 'Time unit: Sqrt(L0/(gBeta*DeltaT)) =', real(timeUnit,kind=8)
  write(00,*) 'Velocity unit: Sqrt(gBeta*L0*DeltaT) =', real(velocityUnit,kind=8)
  write(00,*) 'tauf =', real(tauf,kind=8)
#ifdef EnableLegacyThermalScheme
  write(00,*) 'thermalScheme = legacy D3Q7 (paper-based MRT branch)'
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
#ifdef steadyFlow
  write(00,*) 'I am steadyFlow'
#endif
#ifdef unsteadyFlow
  write(00,*) 'I am unsteadyFlow'
#endif
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
#ifdef VerticalWallsConstT
  write(00,*) 'Temperature B.C. for vertical walls are: ===Hot/cold wall==='
#endif
#ifdef VerticalWallsAdiabatic
  write(00,*) 'Temperature B.C. for vertical walls are: ===Adiabatic wall==='
#endif
#ifdef VerticalWallsPeriodicalT
  write(00,*) 'Temperature B.C. for vertical walls are: ===Periodical==='
#endif
#ifdef HorizontalWallsConstT
  write(00,*) 'Temperature B.C. for horizontal walls are: ===Hot/cold wall==='
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
  write(00,*) 'OpenMP only; MPI/OpenACC are not included in this file'

  !-----------------------------------------------------------------------------------------------
  ! 初始化
  !-----------------------------------------------------------------------------------------------
  rho = 1.0d0
  Fx = 0.0d0
  Fy = 0.0d0
  Fz = 0.0d0
  f = 0.0d0
  g = 0.0d0
  f_post = 0.0d0
  g_post = 0.0d0

  if (loadInitField .EQ. 0) then
    ! 解析初值：速度先全部置零，温度先给导热基态
    write(00,*) 'Initial field is set exactly'
    if (reloadDimensionlessTime .NE. 0) then
      write(00,*) 'Error: since loadInitField .EQ. 0, reloadDimensionlessTime should also be 0'
      close(00)
      stop
    endif

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
#endif

#ifdef HorizontalWallsConstT
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          T(i,j,k) = Thot + (yp(j) - yp(0)) / (yp(ny+1) - yp(0)) * (Tcold - Thot)
        enddo
      enddo
    enddo
#ifdef RayleighBenardCell
    if (Rayleigh .LE. 1.0d4) then
      ! 小 Ra 时导热解附近的增长很慢，给一个很小的初始温度扰动更容易触发对流模态
      ! 这里扰动只沿 x-y 变化，z 方向保持均匀，先对应最简单的卷胞结构
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

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ! 用当前宏观场直接生成 feq / geq，作为初始分布函数
          call compute_feq(rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), feq)
          call compute_geq(T(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), geq)
          do alpha = 0, qf-1
            f(alpha,i,j,k) = feq(alpha+1)
          enddo
          do alpha = 0, qt-1
            g(alpha,i,j,k) = geq(alpha+1)
          enddo
          Bx_prev(i,j,k) = u(i,j,k) * T(i,j,k)
          By_prev(i,j,k) = v(i,j,k) * T(i,j,k)
          Bz_prev(i,j,k) = w(i,j,k) * T(i,j,k)
        enddo
      enddo
    enddo

  elseif (loadInitField .EQ. 1) then
    ! 重启时只读 f 和 g，宏观量后面统一重构，避免文件格式反复变化
    write(00,*) 'Load initial field from previous simulation: ', trim(reloadFilePrefix), '- >>>'
    write(reloadFileName,'(i0)') reloadbinFileNum
    open(unit=01, file=trim(reloadFilePrefix)//'-'//trim(adjustl(reloadFileName))//'.bin', &
         form='unformatted', access='sequential', status='old')
    write(00,*) 'Reloading f and g from file'
    read(01) ((((f(alpha,i,j,k), i=1,nx), j=1,ny), k=1,nz), alpha=0,qf-1)
    read(01) ((((g(alpha,i,j,k), i=1,nx), j=1,ny), k=1,nz), alpha=0,qt-1)
    close(01)
    call reconstruct_macro_from_fg()
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

  binFileNum = 0
  pltFileNum = 0
  dimensionlessTime = 0
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

  ! D3Q19 顺序按参考文献 Eq. (6)：先是 6 个轴向速度，再依次是 xy / xz / yz 三组面对角速度
  ! 注意这与一些 LBM 教材里常见的编号不同，opp(alpha) 也要和这套编号同步修改
  ex  = (/ 0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0 /)
  ey  = (/ 0,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1 /)
  ez  = (/ 0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1 /)
  opp = (/ 0,  2,  1,  4,  3,  6,  5, 10,  9,  8,  7, 14, 13, 12, 11, 18, 17, 16, 15 /)

  omega(0) = 1.0d0 / 3.0d0
  omega(1:6) = 1.0d0 / 18.0d0
  omega(7:18) = 1.0d0 / 36.0d0

  ! D3Q7 温度场：顺序与文献 Eq. (15) 一致
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
! 子程序: init_mrt_matrices
! 作用: 初始化 D3Q19 / D3Q7 的变换矩阵及其逆矩阵。
!===========================================================================================================================
subroutine init_mrt_matrices()
  use commondata3d
  implicit none

  ! M / Minv 在程序启动时一次性构造完，后续碰撞步骤直接拿来做矩空间变换
  ! 这里不再调用通用求逆，而是利用文献这套正交基底的已知模长，直接给出逆矩阵
  call build_basis_matrix_f(M19)
  call init_inverse_matrix_f()

  call build_basis_matrix_g(M7)
  call init_inverse_matrix_g()

end subroutine init_mrt_matrices


!===========================================================================================================================
! 子程序: init_inverse_matrix_f
! 作用: 根据当前正交基底直接构造 D3Q19 逆矩阵。
!===========================================================================================================================
subroutine init_inverse_matrix_f()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j
  real(kind=8), parameter :: invRowNorm2(qf) = (/ &
       1.0d0/19.0d0,   1.0d0/2394.0d0, 1.0d0/252.0d0, &
       1.0d0/10.0d0,   1.0d0/40.0d0,   1.0d0/10.0d0,   1.0d0/40.0d0, &
       1.0d0/10.0d0,   1.0d0/40.0d0,   1.0d0/36.0d0,   1.0d0/72.0d0, &
       1.0d0/12.0d0,   1.0d0/24.0d0,   1.0d0/4.0d0,    1.0d0/4.0d0, &
       1.0d0/4.0d0,    1.0d0/8.0d0,    1.0d0/8.0d0,    1.0d0/8.0d0 /)

  ! 对正交基底有 Minv = transpose(M) * diag(1/||row_i||^2)
  ! 这部分是 3D 相比 2DRB 额外多出来的“矩阵准备”，本质上仍然是在给后面的 MRT 碰撞提供 M 和 Minv。
  do i = 1, qf
    do j = 1, qf
      Minv19(j,i) = M19(i,j) * invRowNorm2(i)
    enddo
  enddo

end subroutine init_inverse_matrix_f


!===========================================================================================================================
! 子程序: init_inverse_matrix_g
! 作用: 根据当前正交基底直接构造 D3Q7 逆矩阵。
!===========================================================================================================================
subroutine init_inverse_matrix_g()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j
  real(kind=8), parameter :: invRowNorm2(qt) = (/ &
       1.0d0/7.0d0, 1.0d0/2.0d0, 1.0d0/2.0d0, 1.0d0/2.0d0, &
       1.0d0/42.0d0, 1.0d0/12.0d0, 1.0d0/4.0d0 /)

  do i = 1, qt
    do j = 1, qt
      Minv7(j,i) = M7(i,j) * invRowNorm2(i)
    enddo
  enddo

end subroutine init_inverse_matrix_g


!===========================================================================================================================
! 子程序: build_basis_matrix_f
! 作用: 按既定矩基顺序填写 D3Q19 变换矩阵。
!===========================================================================================================================
subroutine build_basis_matrix_f(M)
  use commondata3d
  implicit none

  real(kind=8), intent(out) :: M(qf,qf)
  integer(kind=4) :: alpha
  real(kind=8) :: exa, eya, eza, ex2, ey2, ez2, e2, e4

  ! 这里直接按文献给出的 D3Q19 正交基底来生成 M19，和文献中的速度编号、矩空间定义保持一一对应。
  ! 读代码时可以把它理解成 3D 版的“矩空间变换准备”，角色和 2DRB 里的常量定义类似。
  M = 0.0d0
  do alpha = 0, qf-1
    exa = dble(ex(alpha))
    eya = dble(ey(alpha))
    eza = dble(ez(alpha))
    ex2 = exa * exa
    ey2 = eya * eya
    ez2 = eza * eza
    e2 = ex2 + ey2 + ez2
    e4 = e2 * e2

    M( 1,alpha+1) = 1.0d0
    M( 2,alpha+1) = 19.0d0 * e2 - 30.0d0
    M( 3,alpha+1) = 0.5d0 * (21.0d0 * e4 - 53.0d0 * e2 + 24.0d0)
    M( 4,alpha+1) = exa
    M( 5,alpha+1) = (5.0d0 * e2 - 9.0d0) * exa
    M( 6,alpha+1) = eya
    M( 7,alpha+1) = (5.0d0 * e2 - 9.0d0) * eya
    M( 8,alpha+1) = eza
    M( 9,alpha+1) = (5.0d0 * e2 - 9.0d0) * eza
    M(10,alpha+1) = 3.0d0 * ex2 - e2
    M(11,alpha+1) = (3.0d0 * e2 - 5.0d0) * (3.0d0 * ex2 - e2)
    M(12,alpha+1) = ey2 - ez2
    M(13,alpha+1) = (3.0d0 * e2 - 5.0d0) * (ey2 - ez2)
    M(14,alpha+1) = exa * eya
    M(15,alpha+1) = eya * eza
    M(16,alpha+1) = exa * eza
    M(17,alpha+1) = (ey2 - ez2) * exa
    M(18,alpha+1) = (ez2 - ex2) * eya
    M(19,alpha+1) = (ex2 - ey2) * eza
  enddo

end subroutine build_basis_matrix_f


!===========================================================================================================================
! 子程序: build_basis_matrix_g
! 作用: 按既定矩基顺序填写 D3Q7 变换矩阵。
!===========================================================================================================================
subroutine build_basis_matrix_g(M)
  use commondata3d
  implicit none

  real(kind=8), intent(out) :: M(qt,qt)
  integer(kind=4) :: alpha
  real(kind=8) :: exa, eya, eza, ex2, ey2, ez2, e2

  ! 这里按文献 Eq. (16) 的 7 个基底逐行生成 M7。
  ! legacy/current 两个温度分支共用这一套矩空间基底，只是平衡矩和松弛参数不同。
  M = 0.0d0
  do alpha = 0, qt-1
    exa = dble(exT(alpha))
    eya = dble(eyT(alpha))
    eza = dble(ezT(alpha))
    ex2 = exa * exa
    ey2 = eya * eya
    ez2 = eza * eza
    e2 = ex2 + ey2 + ez2

    M(1,alpha+1) = 1.0d0
    M(2,alpha+1) = exa
    M(3,alpha+1) = eya
    M(4,alpha+1) = eza
    M(5,alpha+1) = -6.0d0 + 7.0d0 * e2
    M(6,alpha+1) = 3.0d0 * ex2 - e2
    M(7,alpha+1) = ey2 - ez2
  enddo

end subroutine build_basis_matrix_g











!===========================================================================================================================
! 子程序: compute_feq
! 作用: 根据宏观量计算 D3Q19 流场平衡分布函数。
!===========================================================================================================================
subroutine compute_feq(rhoLoc, uLoc, vLoc, wLoc, feq)
  use commondata3d
  implicit none

  real(kind=8), intent(in) :: rhoLoc, uLoc, vLoc, wLoc
  real(kind=8), intent(out) :: feq(qf)

  integer(kind=4) :: alpha
  real(kind=8) :: eu, u2

  ! 标准二阶平衡分布：保留到速度平方项
  u2 = uLoc * uLoc + vLoc * vLoc + wLoc * wLoc
  do alpha = 0, qf-1
    eu = dble(ex(alpha)) * uLoc + dble(ey(alpha)) * vLoc + dble(ez(alpha)) * wLoc
    feq(alpha+1) = omega(alpha) * rhoLoc * (1.0d0 + 3.0d0 * eu + 4.5d0 * eu * eu - 1.5d0 * u2)
  enddo

end subroutine compute_feq


!===========================================================================================================================
! 子程序: compute_geq
! 作用: 根据温度和速度计算 D3Q7 温度平衡分布函数。
!===========================================================================================================================
subroutine compute_geq(TLoc, uLoc, vLoc, wLoc, geq)
  use commondata3d
  implicit none

  real(kind=8), intent(in) :: TLoc, uLoc, vLoc, wLoc
  real(kind=8), intent(out) :: geq(qt)

  integer(kind=4) :: alpha
  real(kind=8) :: eu

  ! 温度场采用线性平衡分布；不同热算法分支体现在权重和线性系数上
  do alpha = 0, qt-1
    eu = dble(exT(alpha)) * uLoc + dble(eyT(alpha)) * vLoc + dble(ezT(alpha)) * wLoc
    geq(alpha+1) = omegaT(alpha) * TLoc * (1.0d0 + thermalGeqCoeff * eu)
  enddo

end subroutine compute_geq




!===========================================================================================================================
! 子程序: collision
! 作用: 流场碰撞步骤，在矩空间完成松弛并加入浮力源项修正。
! 用途: 在主程序时间推进循环中调用，位于 streaming 之前。
!===========================================================================================================================
subroutine collision()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k, alpha
  real(kind=8) :: fLoc(qf), m(qf), meq(qf), mPost(qf), fPostLoc(qf)
  real(kind=8) :: rhoLoc, uLoc, vLoc, wLoc, u2, uDotF
  real(kind=8) :: FxLoc, FyLoc, FzLoc
  ! 流场仍采用 D3Q19 的矩空间碰撞：
  ! 先做 m = M*f，再按文献的平衡矩 meq 和 commondata3d 中的松弛率逐个碰撞，
  ! 体力项也先投到矩空间，再乘以 (I-S/2) 做半步修正。

  !$omp parallel do collapse(3) default(none) shared(f,f_post,rho,u,v,w,T,Fx,Fy,Fz,M19,Minv19) &
  !$omp private(i,j,k,alpha,fLoc,m,meq,mPost,fPostLoc,rhoLoc,uLoc,vLoc,wLoc,u2,uDotF,FxLoc,FyLoc,FzLoc)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        ! ---- Step 1: 读取当前格点的分布函数与宏观量 ----
        do alpha = 0, qf-1
          fLoc(alpha+1) = f(alpha,i,j,k)
        enddo

        rhoLoc = rho(i,j,k)
        uLoc = u(i,j,k)
        vLoc = v(i,j,k)
        wLoc = w(i,j,k)
        u2 = uLoc * uLoc + vLoc * vLoc + wLoc * wLoc

        ! ---- Step 2: 计算体力项 ----
        FxLoc = 0.0d0
        FyLoc = rhoLoc * gBeta * (T(i,j,k) - Tref)
        FzLoc = 0.0d0
        Fx(i,j,k) = FxLoc
        Fy(i,j,k) = FyLoc
        Fz(i,j,k) = FzLoc
        uDotF = uLoc * FxLoc + vLoc * FyLoc + wLoc * FzLoc

        ! x / z 方向无体力，RB 浮力只在 y 方向

        ! ---- Step 3: 把分布函数变到矩空间，并构造平衡矩 ----
        m = matmul(M19, fLoc)

        meq(1)  = rhoLoc
        meq(2)  = rhoLoc * (-11.0d0 + 19.0d0 * u2)
        meq(3)  = rhoLoc * (3.0d0 - 11.0d0 * u2 / 2.0d0)
        meq(4)  = rhoLoc * uLoc
        meq(5)  = -2.0d0 * rhoLoc * uLoc / 3.0d0
        meq(6)  = rhoLoc * vLoc
        meq(7)  = -2.0d0 * rhoLoc * vLoc / 3.0d0
        meq(8)  = rhoLoc * wLoc
        meq(9)  = -2.0d0 * rhoLoc * wLoc / 3.0d0
        meq(10) = rhoLoc * (2.0d0 * uLoc * uLoc - vLoc * vLoc - wLoc * wLoc)
        meq(11) = -0.5d0 * meq(10)
        meq(12) = rhoLoc * (vLoc * vLoc - wLoc * wLoc)
        meq(13) = -0.5d0 * meq(12)
        meq(14) = rhoLoc * uLoc * vLoc
        meq(15) = rhoLoc * vLoc * wLoc
        meq(16) = rhoLoc * uLoc * wLoc
        meq(17) = 0.0d0
        meq(18) = 0.0d0
        meq(19) = 0.0d0

        ! ---- Step 4: 在矩空间完成碰撞，并加入体力项修正 ----
        mPost(1)  = m(1)
        mPost(2)  = m(2)  - Se    * (m(2)  - meq(2))  + (1.0d0 - 0.5d0 * Se   ) * 38.0d0 * uDotF
        mPost(3)  = m(3)  - Seps  * (m(3)  - meq(3))  + (1.0d0 - 0.5d0 * Seps ) * (-11.0d0) * uDotF
        mPost(4)  = m(4)  + FxLoc
        mPost(5)  = m(5)  - Sq    * (m(5)  - meq(5))  + (1.0d0 - 0.5d0 * Sq   ) * (-2.0d0 / 3.0d0) * FxLoc
        mPost(6)  = m(6)  + FyLoc
        mPost(7)  = m(7)  - Sq    * (m(7)  - meq(7))  + (1.0d0 - 0.5d0 * Sq   ) * (-2.0d0 / 3.0d0) * FyLoc
        mPost(8)  = m(8)  + FzLoc
        mPost(9)  = m(9)  - Sq    * (m(9)  - meq(9))  + (1.0d0 - 0.5d0 * Sq   ) * (-2.0d0 / 3.0d0) * FzLoc
        mPost(10) = m(10) - Snu   * (m(10) - meq(10)) + &
             (1.0d0 - 0.5d0 * Snu) * &
             (4.0d0 * uLoc * FxLoc - 2.0d0 * vLoc * FyLoc - 2.0d0 * wLoc * FzLoc)
        mPost(11) = m(11) - Spi   * (m(11) - meq(11)) + &
             (1.0d0 - 0.5d0 * Spi) * &
             (-2.0d0 * uLoc * FxLoc + vLoc * FyLoc + wLoc * FzLoc)
        mPost(12) = m(12) - Snu   * (m(12) - meq(12)) + &
             (1.0d0 - 0.5d0 * Snu) * &
             (2.0d0 * vLoc * FyLoc - 2.0d0 * wLoc * FzLoc)
        mPost(13) = m(13) - Spi   * (m(13) - meq(13)) + (1.0d0 - 0.5d0 * Spi) * (-vLoc * FyLoc + wLoc * FzLoc)
        mPost(14) = m(14) - Snu   * (m(14) - meq(14)) + (1.0d0 - 0.5d0 * Snu) * (uLoc * FyLoc + vLoc * FxLoc)
        mPost(15) = m(15) - Snu   * (m(15) - meq(15)) + (1.0d0 - 0.5d0 * Snu) * (vLoc * FzLoc + wLoc * FyLoc)
        mPost(16) = m(16) - Snu   * (m(16) - meq(16)) + (1.0d0 - 0.5d0 * Snu) * (uLoc * FzLoc + wLoc * FxLoc)
        mPost(17) = m(17) - Sm    * m(17)
        mPost(18) = m(18) - Sm    * m(18)
        mPost(19) = m(19) - Sm    * m(19)

        ! ---- Step 5: 由逆矩阵回到速度空间 ----
        fPostLoc = matmul(Minv19, mPost)

        do alpha = 0, qf-1
          f_post(alpha,i,j,k) = fPostLoc(alpha+1)
        enddo
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

  integer(kind=4) :: i, j, k, ip, jp, kp, alpha

#ifdef VerticalWallsPeriodicalU
  !$omp parallel do collapse(2) default(none) shared(f_post) private(j,k,alpha)
  do k = 1, nz
    do j = 1, ny
      do alpha = 0, qf-1
        f_post(alpha,0,j,k) = f_post(alpha,nx,j,k)
        f_post(alpha,nx+1,j,k) = f_post(alpha,1,j,k)
      enddo
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef SpanwiseWallsPeriodicalU
  !$omp parallel do collapse(2) default(none) shared(f_post) private(i,j,alpha)
  do j = 1, ny
    do i = 0, nx+1
      do alpha = 0, qf-1
        f_post(alpha,i,j,0) = f_post(alpha,i,j,nz)
        f_post(alpha,i,j,nz+1) = f_post(alpha,i,j,1)
      enddo
    enddo
  enddo
  !$omp end parallel do
#endif

  ! pull streaming：当前格点 (i,j,k) 从上游格点 (i-ex, j-ey, k-ez) 拉取分布函数
  !$omp parallel do collapse(3) default(none) shared(f,f_post,ex,ey,ez) private(i,j,k,ip,jp,kp,alpha)
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

  integer(kind=4) :: i, j, k, alpha

#ifdef VerticalWallsNoslip
  !$omp parallel do collapse(2) default(none) shared(f,f_post,ex,opp) private(j,k,alpha)
  do k = 1, nz
    do j = 1, ny
      do alpha = 0, qf-1
        if (ex(alpha) .EQ. 1)  f(alpha,1,j,k)  = f_post(opp(alpha),1,j,k)
        if (ex(alpha) .EQ. -1) f(alpha,nx,j,k) = f_post(opp(alpha),nx,j,k)
      enddo
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef HorizontalWallsNoslip
  !$omp parallel do collapse(2) default(none) shared(f,f_post,ey,opp) private(i,k,alpha)
  do k = 1, nz
    do i = 1, nx
      do alpha = 0, qf-1
        if (ey(alpha) .EQ. 1)  f(alpha,i,1,k)  = f_post(opp(alpha),i,1,k)
        if (ey(alpha) .EQ. -1) f(alpha,i,ny,k) = f_post(opp(alpha),i,ny,k)
      enddo
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef SpanwiseWallsNoslip
  !$omp parallel do collapse(2) default(none) shared(f,f_post,ez,opp) private(i,j,alpha)
  do j = 1, ny
    do i = 1, nx
      do alpha = 0, qf-1
        if (ez(alpha) .EQ. 1)  f(alpha,i,j,1)  = f_post(opp(alpha),i,j,1)
        if (ez(alpha) .EQ. -1) f(alpha,i,j,nz) = f_post(opp(alpha),i,j,nz)
      enddo
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

  integer(kind=4) :: i, j, k, alpha
  real(kind=8) :: momx, momy, momz, rhoLoc

  ! ---- Step 1: 由分布函数恢复密度与动量 ----
  ! ---- Step 2: 用 0.5F 做半步修正，恢复三分量速度 ----
  !$omp parallel do collapse(3) default(none) &
  !$omp& shared(f,rho,u,v,w,Fx,Fy,Fz,ex,ey,ez) &
  !$omp& private(i,j,k,alpha,momx,momy,momz,rhoLoc)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        rhoLoc = 0.0d0
        momx = 0.0d0
        momy = 0.0d0
        momz = 0.0d0
        do alpha = 0, qf-1
          rhoLoc = rhoLoc + f(alpha,i,j,k)
          momx = momx + f(alpha,i,j,k) * dble(ex(alpha))
          momy = momy + f(alpha,i,j,k) * dble(ey(alpha))
          momz = momz + f(alpha,i,j,k) * dble(ez(alpha))
        enddo
        rho(i,j,k) = rhoLoc
        u(i,j,k) = (momx + 0.5d0 * Fx(i,j,k)) / rhoLoc
        v(i,j,k) = (momy + 0.5d0 * Fy(i,j,k)) / rhoLoc
        w(i,j,k) = (momz + 0.5d0 * Fz(i,j,k)) / rhoLoc
      enddo
    enddo
  enddo
  !$omp end parallel do

end subroutine macro


!===========================================================================================================================
! 子程序: collisionT
! 作用: 温度场碰撞步骤；legacy 分支走文献 D3Q7-MRT，current 分支走现有 EnableUseG 方案。
! 用途: 在主程序时间推进循环中调用，位于流场 macro 之后。
!===========================================================================================================================
subroutine collisionT()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k, alpha
  real(kind=8) :: gLoc(qt), n(qt), neq(qt), q(qt), nPost(qt), gPostLoc(qt)
  real(kind=8) :: Bx, By, Bz, dBx, dBy, dBz
  real(kind=8), parameter :: SG = 1.0d0 - 0.5d0 * Qk

  ! Step 1: 读取当前格点的温度分布函数，并构造热流 B = uT。
  ! Step 2: 若启用 useG，则用当前时刻与上一时刻的热流差 dB 做修正。
  ! Step 3: 把 g 变换到矩空间，构造平衡矩与松弛率。
  ! Step 4: 在矩空间完成碰撞，再逆变换回速度空间写入 g_post。
  !$omp parallel do collapse(3) default(none) shared(g,g_post,u,v,w,T,Bx_prev,By_prev,Bz_prev,M7,Minv7) &
  !$omp private(i,j,k,alpha,gLoc,n,neq,q,nPost,gPostLoc,Bx,By,Bz,dBx,dBy,dBz)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        do alpha = 0, qt-1
          gLoc(alpha+1) = g(alpha,i,j,k)
        enddo

        Bx = u(i,j,k) * T(i,j,k)
        By = v(i,j,k) * T(i,j,k)
        Bz = w(i,j,k) * T(i,j,k)

        if (useG) then
          dBx = Bx - Bx_prev(i,j,k)
          dBy = By - By_prev(i,j,k)
          dBz = Bz - Bz_prev(i,j,k)
          Bx_prev(i,j,k) = Bx
          By_prev(i,j,k) = By
          Bz_prev(i,j,k) = Bz
        else
          dBx = 0.0d0
          dBy = 0.0d0
          dBz = 0.0d0
        endif

        ! 先把温度分布函数变到矩空间
        n = matmul(M7, gLoc)

        ! 平衡矩按所选热算法分支写入：
        ! legacy 分支对应文献中的 [T, uT, vT, wT, aT*T, 0, 0]
        ! current 分支对应现有 EnableUseG 方案
        neq(1) = T(i,j,k)
        neq(2) = Bx
        neq(3) = By
        neq(4) = Bz
#ifdef EnableLegacyThermalScheme
        neq(5) = paraA * T(i,j,k)
#else
        neq(5) = -0.75d0 * T(i,j,k)
#endif
        neq(6) = 0.0d0
        neq(7) = 0.0d0

        q(1) = 0.0d0
        q(2:4) = Qk
        q(5:7) = Qnu

        nPost(1) = n(1) - q(1) * (n(1) - neq(1))
        nPost(2:4) = n(2:4) - q(2:4) * (n(2:4) - neq(2:4))
        nPost(5:7) = n(5:7) - q(5:7) * (n(5:7) - neq(5:7))
        nPost(2) = nPost(2) + SG * dBx
        nPost(3) = nPost(3) + SG * dBy
        nPost(4) = nPost(4) + SG * dBz

        gPostLoc = matmul(Minv7, nPost)

        do alpha = 0, qt-1
          g_post(alpha,i,j,k) = gPostLoc(alpha+1)
        enddo
      enddo
    enddo
  enddo
  !$omp end parallel do

  return
end subroutine collisionT
!===========================================================================================================================
! collisionT 结束: 完成温度分布函数 g 的碰撞更新，并加入热流修正项。
!===========================================================================================================================


!===========================================================================================================================
! 子程序: streamingT
! 作用: 对温度分布函数执行三维 pull streaming。
! 用途: 在主程序时间推进循环中调用，位于 collisionT 之后、bouncebackT 之前。
!===========================================================================================================================
subroutine streamingT()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k, ip, jp, kp, alpha

  ! Step 1: 若某些方向是周期边界，则先补齐 g_post 的 ghost layer，
  !         这样后续的 pull streaming 就能沿用和内部点完全一致的写法。
#ifdef VerticalWallsPeriodicalT
  !$omp parallel do collapse(2) default(none) shared(g_post) private(j,k,alpha)
  do k = 1, nz
    do j = 1, ny
      do alpha = 0, qt-1
        g_post(alpha,0,j,k) = g_post(alpha,nx,j,k)
        g_post(alpha,nx+1,j,k) = g_post(alpha,1,j,k)
      enddo
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef SpanwiseWallsPeriodicalT
  !$omp parallel do collapse(2) default(none) shared(g_post) private(i,j,alpha)
  do j = 1, ny
    do i = 0, nx+1
      do alpha = 0, qt-1
        g_post(alpha,i,j,0) = g_post(alpha,i,j,nz)
        g_post(alpha,i,j,nz+1) = g_post(alpha,i,j,1)
      enddo
    enddo
  enddo
  !$omp end parallel do
#endif

  ! Step 2: 对内部节点执行 pull streaming，把碰撞后的温度信息拉回当前格点。
  ! 温度场同样采用 pull streaming
  !$omp parallel do collapse(3) default(none) shared(g,g_post,exT,eyT,ezT) private(i,j,k,ip,jp,kp,alpha)
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

  integer(kind=4) :: i, j, k, alpha

  ! 周期边界已经在 streamingT 中通过 ghost layer 处理完毕，
  ! 这里继续按 2DRB 的思路处理恒温和绝热边界。
#ifdef VerticalWallsConstT
  !$omp parallel do collapse(2) default(none) shared(g,g_post,exT,oppT,omegaT) private(j,k,alpha)
  do k = 1, nz
    do j = 1, ny
      do alpha = 0, qt-1
        if (exT(alpha) .EQ. 1)  g(alpha,1,j,k)  = -g_post(oppT(alpha),1,j,k)  + 2.0d0 * omegaT(alpha) * Thot
        if (exT(alpha) .EQ. -1) g(alpha,nx,j,k) = -g_post(oppT(alpha),nx,j,k) + 2.0d0 * omegaT(alpha) * Tcold
      enddo
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef VerticalWallsAdiabatic
  !$omp parallel do collapse(2) default(none) shared(g,g_post,exT,oppT) private(j,k,alpha)
  do k = 1, nz
    do j = 1, ny
      do alpha = 0, qt-1
        if (exT(alpha) .EQ. 1)  g(alpha,1,j,k)  = g_post(oppT(alpha),1,j,k)
        if (exT(alpha) .EQ. -1) g(alpha,nx,j,k) = g_post(oppT(alpha),nx,j,k)
      enddo
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef HorizontalWallsAdiabatic
  !$omp parallel do collapse(2) default(none) shared(g,g_post,eyT,oppT) private(i,k,alpha)
  do k = 1, nz
    do i = 1, nx
      do alpha = 0, qt-1
        if (eyT(alpha) .EQ. 1)  g(alpha,i,1,k)  = g_post(oppT(alpha),i,1,k)
        if (eyT(alpha) .EQ. -1) g(alpha,i,ny,k) = g_post(oppT(alpha),i,ny,k)
      enddo
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef HorizontalWallsConstT
  !$omp parallel do collapse(2) default(none) shared(g,g_post,eyT,oppT,omegaT) private(i,k,alpha)
  do k = 1, nz
    do i = 1, nx
      do alpha = 0, qt-1
        if (eyT(alpha) .EQ. 1)  g(alpha,i,1,k)  = -g_post(oppT(alpha),i,1,k)  + 2.0d0 * omegaT(alpha) * Thot
        if (eyT(alpha) .EQ. -1) g(alpha,i,ny,k) = -g_post(oppT(alpha),i,ny,k) + 2.0d0 * omegaT(alpha) * Tcold
      enddo
    enddo
  enddo
  !$omp end parallel do
#endif

#ifdef SpanwiseWallsAdiabatic
  !$omp parallel do collapse(2) default(none) shared(g,g_post,ezT,oppT) private(i,j,alpha)
  do j = 1, ny
    do i = 1, nx
      do alpha = 0, qt-1
        if (ezT(alpha) .EQ. 1)  g(alpha,i,j,1)  = g_post(oppT(alpha),i,j,1)
        if (ezT(alpha) .EQ. -1) g(alpha,i,j,nz) = g_post(oppT(alpha),i,j,nz)
      enddo
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

  integer(kind=4) :: i, j, k, alpha, iter
  real(kind=8) :: momx, momy, momz
  logical :: rho_bad

  ! 重启时只存了 f 和 g，所以这里统一把 T、rho、u、v、w 以及历史热流都重构回来
  call macroT()
  rho_bad = .false.

  !$omp parallel do collapse(3) default(none) shared(f,rho,u,v,w,T,Fx,Fy,Fz,Bx_prev,By_prev,Bz_prev,ex,ey,ez) &
  !$omp private(i,j,k,alpha,iter,momx,momy,momz) reduction(.or.:rho_bad)
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
            Fx(i,j,k) = 0.0d0
            Fy(i,j,k) = rho(i,j,k) * gBeta * (T(i,j,k) - Tref)
            Fz(i,j,k) = 0.0d0
            u(i,j,k) = (momx + 0.5d0 * Fx(i,j,k)) / rho(i,j,k)
            v(i,j,k) = (momy + 0.5d0 * Fy(i,j,k)) / rho(i,j,k)
            w(i,j,k) = (momz + 0.5d0 * Fz(i,j,k)) / rho(i,j,k)
          enddo
        else
          rho_bad = .true.
          u(i,j,k) = 0.0d0
          v(i,j,k) = 0.0d0
          w(i,j,k) = 0.0d0
          Fx(i,j,k) = 0.0d0
          Fy(i,j,k) = 0.0d0
          Fz(i,j,k) = 0.0d0
        endif

        if (useG) then
          Bx_prev(i,j,k) = u(i,j,k) * T(i,j,k)
          By_prev(i,j,k) = v(i,j,k) * T(i,j,k)
          Bz_prev(i,j,k) = w(i,j,k) * T(i,j,k)
        else
          Bx_prev(i,j,k) = 0.0d0
          By_prev(i,j,k) = 0.0d0
          Bz_prev(i,j,k) = 0.0d0
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
! reconstruct_macro_from_fg end: restart state is fully rebuilt from the reloaded distributions
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

  call append_convergence_tecplot('convergence3D.plt', itc, errorU, errorT)
  write(caseTag,'("Ra=",ES10.3E2,",nx=",I0,",ny=",I0,",nz=",I0,",useG=",L1,",useLegacyThermalScheme=",L1)') &
       Rayleigh, nx, ny, nz, useG, useLegacyThermalScheme
  call append_convergence_master_tecplot('convergence_all_3D.plt', caseTag, itc, errorU, errorT)
  write(*,'(I12,1X,ES24.16,1X,ES24.16)') itc, errorU, errorT

end subroutine check
#endif


!===========================================================================================================================
! 子程序: append_convergence_tecplot
! 作用: 向单个收敛历史文件追加一条误差记录。
!===========================================================================================================================
subroutine append_convergence_tecplot(filename, itcLoc, errorULoc, errorTLoc)
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

end subroutine append_convergence_tecplot


!===========================================================================================================================
! 子程序: append_convergence_master_tecplot
! 作用: 向带 zone 名称的收敛历史文件追加一条记录。
!===========================================================================================================================
subroutine append_convergence_master_tecplot(filename, zoneName, itcLoc, errorULoc, errorTLoc)
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

end subroutine append_convergence_master_tecplot


!===========================================================================================================================
! 子程序: output_binary
! 作用: 输出三维快照二进制文件，供后处理或继续分析使用。
! 用途: 在运行过程中按需调用，也在程序结束时调用。
!===========================================================================================================================
subroutine output_binary()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k
  character(len=100) :: filename

  ! This snapshot is for post-processing only; u/v/w are written after nondimensionalization.
  ! For strict restart, keep using backupData(), which preserves the lattice-state variables.
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
end subroutine output_binary
!===========================================================================================================================
! output_binary 结束: 输出 u、v、w、T、rho 的二进制快照文件。
!===========================================================================================================================


!===========================================================================================================================
! 子程序: backupData
! 作用: 备份 f/g 分布函数，供后续重启继续计算。
! 用途: 在运行过程中定期调用，也在程序结束前调用。
!===========================================================================================================================
subroutine backupData()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k, alpha
  character(len=100) :: filename

  ! Strict restart snapshots store only f and g; rho/u/v/w/T are reconstructed after reload.
#ifdef steadyFlow
  write(filename,'(i0)') itc
#endif
#ifdef unsteadyFlow
  if (loadInitField .EQ. 0) write(filename,'(i0)') binFileNum
  if (loadInitField .EQ. 1) write(filename,'(i0)') binFileNum + reloadbinFileNum
#endif

  filename = adjustl(filename)
  open(unit=05, file='backupFile3D-'//trim(filename)//'.bin', form='unformatted', access='sequential')
  write(05) ((((f(alpha,i,j,k), i=1,nx), j=1,ny), k=1,nz), alpha=0,qf-1)
  write(05) ((((g(alpha,i,j,k), i=1,nx), j=1,ny), k=1,nz), alpha=0,qt-1)
  close(05)

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) 'Backup f and g to the file: backupFile3D-', trim(filename), '.bin'
  close(00)

  return
end subroutine backupData
!===========================================================================================================================
! backupData 结束: 输出包含 f、g 的重启备份文件。
!===========================================================================================================================


!===========================================================================================================================
! 子程序: output_Tecplot
! 作用: 输出 x/y/z 三个中面切片，便于快速查看三维流场结构。
! 用途: 在运行过程中按需调用，也在程序结束时调用。
!===========================================================================================================================
subroutine output_Tecplot()
  use commondata3d
  implicit none

  character(len=100) :: filename

  ! 3D 这里不直接输出整体 Tecplot 体数据，而是输出三个中面切片；
  ! 入口名称仍保持 output_Tecplot()，这样和 2DRB 主程序调用保持一致。
#ifdef steadyFlow
  write(filename,'(i12.12)') itc
#endif
#ifdef unsteadyFlow
  pltFileNum = pltFileNum + 1
  write(filename,'(i12.12)') pltFileNum
#endif

  filename = adjustl(filename)
  call write_midplane_x(trim(pltFolderPrefix)//'-midX-'//trim(filename)//'.dat')
  call write_midplane_y(trim(pltFolderPrefix)//'-midY-'//trim(filename)//'.dat')
  call write_midplane_z(trim(pltFolderPrefix)//'-midZ-'//trim(filename)//'.dat')

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
  !$omp parallel do collapse(3) default(none) shared(u,T) private(i,j,k) reduction(+:NuVolAvg_temp)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        NuVolAvg_temp = NuVolAvg_temp + u(i,j,k) * T(i,j,k)
      enddo
    enddo
  enddo
  !$omp end parallel do
#else
  !$omp parallel do collapse(3) default(none) shared(v,T) private(i,j,k) reduction(+:NuVolAvg_temp)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        NuVolAvg_temp = NuVolAvg_temp + v(i,j,k) * T(i,j,k)
      enddo
    enddo
  enddo
  !$omp end parallel do
#endif

  NuVolAvg(dimensionlessTime) = NuVolAvg_temp / dble(nx * ny * nz) * lengthUnit / diffusivity + 1.0d0
  open(unit=01, file='Nu_VolAvg_3D.dat', status='unknown', position='append')
  write(01,*) real(reloadDimensionlessTime + dimensionlessTime * outputFrequency, kind=8), NuVolAvg(dimensionlessTime)
  close(01)

  ReVolAvg_temp = 0.0d0
  !$omp parallel do collapse(3) default(none) shared(u,v,w) private(i,j,k) reduction(+:ReVolAvg_temp)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        ReVolAvg_temp = ReVolAvg_temp + u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)
      enddo
    enddo
  enddo
  !$omp end parallel do
  ReVolAvg(dimensionlessTime) = dsqrt(ReVolAvg_temp / dble(nx * ny * nz)) * lengthUnit / viscosity
  open(unit=02, file='Re_VolAvg_3D.dat', status='unknown', position='append')
  write(02,*) real(reloadDimensionlessTime + dimensionlessTime * outputFrequency, kind=8), ReVolAvg(dimensionlessTime)
  close(02)

  write(*,'(a,1x,es16.8)') 'NuVolAvg =', NuVolAvg(dimensionlessTime)
  write(*,'(a,1x,es16.8)') 'ReVolAvg =', ReVolAvg(dimensionlessTime)

  return
end subroutine calNuRe
!===========================================================================================================================
! calNuRe 结束: 计算体平均 Nu 和 Re 的时间历程统计量。
!===========================================================================================================================


!===========================================================================================================================
! 子程序: SideHeatedcalc_Nu_global
! 作用: 计算侧壁差温工况下的全场平均 Nusselt 数。
! 用途: 在 SideHeatedCell 工况结束后的后处理中调用。
!===========================================================================================================================
subroutine SideHeatedcalc_Nu_global()
  use commondata3d
  implicit none
  integer(kind=4) :: i, j, k
  real(kind=8) :: dx, dTdx, qx, sum_qx
  real(kind=8) :: deltaT, coef

  dx = 1.0d0 / lengthUnit
  deltaT = Thot - Tcold
  coef = heatFluxScale
  sum_qx = 0.0d0

  !$omp parallel do collapse(3) default(none) shared(u,T,dx,coef) private(i,j,k,dTdx,qx) reduction(+:sum_qx)
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
  !$omp end parallel do

  Nu_global = (sum_qx / dble(nx * ny * nz)) / deltaT

  write(*,'(a,1x,es16.8)') 'Nu_global =', Nu_global
  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,'(a,1x,es16.8)') 'Nu_global =', Nu_global
  close(00)

  return
end subroutine SideHeatedcalc_Nu_global
!===========================================================================================================================
! SideHeatedcalc_Nu_global 结束: 计算侧壁差温工况下的全场平均 Nusselt 数。
!===========================================================================================================================


!===========================================================================================================================
! 子程序: RBcalc_Nu_global
! 作用: 计算 Rayleigh-Benard 工况下的全场平均 Nusselt 数。
! 用途: 在 RayleighBenardCell 工况结束后的后处理中调用。
!===========================================================================================================================
subroutine RBcalc_Nu_global()
  use commondata3d
  implicit none
  integer(kind=4) :: i, j, k
  real(kind=8) :: dy, dTdy, qy, sum_qy
  real(kind=8) :: deltaT, coef

  dy = 1.0d0 / lengthUnit
  deltaT = Thot - Tcold
  coef = heatFluxScale
  sum_qy = 0.0d0

  !$omp parallel do collapse(3) default(none) shared(v,T,dy,coef) private(i,j,k,dTdy,qy) reduction(+:sum_qy)
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
  !$omp end parallel do

  Nu_global = (sum_qy / dble(nx * ny * nz)) / deltaT

  write(*,'(a,1x,es16.8)') 'Nu_global =', Nu_global
  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,'(a,1x,es16.8)') 'Nu_global =', Nu_global
  close(00)

  return
end subroutine RBcalc_Nu_global
!===========================================================================================================================
! RBcalc_Nu_global 结束: 计算 Rayleigh-Benard 工况下的全场平均 Nusselt 数。
!===========================================================================================================================


!===========================================================================================================================
! 子程序: SideHeatedcalc_Nu_wall_avg
! 作用: 计算侧壁差温工况下热壁、冷壁和中面平均 Nusselt 数及其极值。
! 用途: 在 SideHeatedCell 工况结束后的后处理中调用。
!===========================================================================================================================
subroutine SideHeatedcalc_Nu_wall_avg()
  use commondata3d
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

  !$omp parallel do default(none) shared(T,T_left_avg) private(j,k)
  do j = 1, ny
    T_left_avg(j) = 0.0d0
    do k = 1, nz
      T_left_avg(j) = T_left_avg(j) + T(1,j,k)
    enddo
    T_left_avg(j) = T_left_avg(j) / dble(nz)
  enddo
  !$omp end parallel do

  sum_hot = 0.0d0
  !$omp parallel do default(none) shared(T,Nu_left,dx,deltaT) private(j,k,qx_wall) reduction(+:sum_hot)
  do j = 1, ny
    Nu_left(j) = 0.0d0
    do k = 1, nz
      qx_wall = (8.0d0*Thot - 9.0d0*T(1,j,k) + T(2,j,k)) / (3.0d0*dx)
      Nu_left(j) = Nu_left(j) + qx_wall / deltaT
    enddo
    Nu_left(j) = Nu_left(j) / dble(nz)
    sum_hot = sum_hot + Nu_left(j)
  enddo
  !$omp end parallel do
  Nu_hot = sum_hot / dble(ny)

  Nu_left_ext(1:ny) = Nu_left(1:ny)
  yfit(1) = yp(1);  Tfit(1) = T_left_avg(1)
  yfit(2) = yp(2);  Tfit(2) = T_left_avg(2)
  yfit(3) = yp(3);  Tfit(3) = T_left_avg(3)
  yfit(4) = yp(4);  Tfit(4) = T_left_avg(4)
  call fit_adiabatic_wall_T4(0.0d0, yfit, Tfit, T_wb)
  Nu_left_ext(0) = (2.0d0 * (Thot - T_wb) / dx) / deltaT

  yfit(1) = yp(ny-3);  Tfit(1) = T_left_avg(ny-3)
  yfit(2) = yp(ny-2);  Tfit(2) = T_left_avg(ny-2)
  yfit(3) = yp(ny-1);  Tfit(3) = T_left_avg(ny-1)
  yfit(4) = yp(ny  );  Tfit(4) = T_left_avg(ny  )
  call fit_adiabatic_wall_T4(yp(ny+1), yfit, Tfit, T_wt)
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
  call fit_parabola_ls5(yk, fk, +1, fstar, ystar)
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
  call fit_parabola_ls5(yk, fk, -1, fstar, ystar)
  Nu_hot_min = fstar
  Nu_hot_min_position = ystar

  sum_cold = 0.0d0
  !$omp parallel do default(none) shared(T,dx,deltaT) private(j,k,qx_wall) reduction(+:sum_cold)
  do j = 1, ny
    do k = 1, nz
      qx_wall = (-8.0d0*Tcold + 9.0d0*T(nx,j,k) - T(nx-1,j,k)) / (3.0d0*dx)
      sum_cold = sum_cold + qx_wall / deltaT
    enddo
  enddo
  !$omp end parallel do
  Nu_cold = sum_cold / dble(ny * nz)

  sum_mid = 0.0d0
  if (mod(nx,2) .EQ. 1) then
    iMid = (nx + 1) / 2
    !$omp parallel do collapse(2) default(none) shared(u,T,iMid,dx,deltaT,coef) private(j,k) reduction(+:sum_mid)
    do k = 1, nz
      do j = 1, ny
        sum_mid = sum_mid + (coef * u(iMid,j,k) * (T(iMid,j,k) - Tref) + &
             (T(iMid-1,j,k) - T(iMid+1,j,k)) / (2.0d0*dx)) / deltaT
      enddo
    enddo
    !$omp end parallel do
  else
    iL = nx / 2
    iR = iL + 1
    !$omp parallel do collapse(2) default(none) shared(u,T,iL,iR,dx,deltaT,coef) private(j,k) reduction(+:sum_mid)
    do k = 1, nz
      do j = 1, ny
        sum_mid = sum_mid + (coef * 0.5d0 * (u(iL,j,k) * (T(iL,j,k) - Tref) + &
             u(iR,j,k) * (T(iR,j,k) - Tref)) + (T(iL,j,k) - T(iR,j,k)) / dx) / deltaT
      enddo
    enddo
    !$omp end parallel do
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
end subroutine SideHeatedcalc_Nu_wall_avg
!===========================================================================================================================
! SideHeatedcalc_Nu_wall_avg 结束: 计算侧壁差温工况下热壁、冷壁和中面的 Nusselt 数及其极值。
!===========================================================================================================================


!===========================================================================================================================
! 子程序: RBcalc_Nu_wall_avg
! 作用: 计算 Rayleigh-Benard 工况下热壁、冷壁和中面的 Nusselt 数及其极值。
! 用途: 在 RayleighBenardCell 工况结束后的后处理中调用。
!===========================================================================================================================
subroutine RBcalc_Nu_wall_avg()
  use commondata3d
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

  !$omp parallel do default(none) shared(T,T_bot_avg) private(i,k)
  do i = 1, nx
    T_bot_avg(i) = 0.0d0
    do k = 1, nz
      T_bot_avg(i) = T_bot_avg(i) + T(i,1,k)
    enddo
    T_bot_avg(i) = T_bot_avg(i) / dble(nz)
  enddo
  !$omp end parallel do

  sum_hot = 0.0d0
  !$omp parallel do default(none) shared(T,Nu_bot,dy,deltaT) private(i,k,qy_wall) reduction(+:sum_hot)
  do i = 1, nx
    Nu_bot(i) = 0.0d0
    do k = 1, nz
      qy_wall = (8.0d0*Thot - 9.0d0*T(i,1,k) + T(i,2,k)) / (3.0d0*dy)
      Nu_bot(i) = Nu_bot(i) + qy_wall / deltaT
    enddo
    Nu_bot(i) = Nu_bot(i) / dble(nz)
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
  Nu_bot_ext(0) = (2.0d0 * (Thot - T_wl) / dy) / deltaT

  xfit(1) = xp(nx-3);  Tfit(1) = T_bot_avg(nx-3)
  xfit(2) = xp(nx-2);  Tfit(2) = T_bot_avg(nx-2)
  xfit(3) = xp(nx-1);  Tfit(3) = T_bot_avg(nx-1)
  xfit(4) = xp(nx  );  Tfit(4) = T_bot_avg(nx  )
  call fit_adiabatic_wall_T4(xp(nx+1), xfit, Tfit, T_wr)
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
  !$omp parallel do default(none) shared(T,dy,deltaT) private(i,k,qy_wall) reduction(+:sum_cold)
  do i = 1, nx
    do k = 1, nz
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
        sum_mid = sum_mid + (coef * v(i,jMid,k) * (T(i,jMid,k) - Tref) - &
             (T(i,jMid+1,k) - T(i,jMid-1,k)) / (2.0d0*dy)) / deltaT
      enddo
    enddo
    !$omp end parallel do
  else
    jB = ny / 2
    jT = jB + 1
    !$omp parallel do collapse(2) default(none) shared(v,T,jB,jT,dy,deltaT,coef) private(i,k) reduction(+:sum_mid)
    do k = 1, nz
      do i = 1, nx
        sum_mid = sum_mid + (coef * 0.5d0 * (v(i,jB,k) * (T(i,jB,k) - Tref) + &
             v(i,jT,k) * (T(i,jT,k) - Tref)) + (T(i,jB,k) - T(i,jT,k)) / dy) / deltaT
      enddo
    enddo
    !$omp end parallel do
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
end subroutine RBcalc_Nu_wall_avg
!===========================================================================================================================
! RBcalc_Nu_wall_avg 结束: 计算 Rayleigh-Benard 工况下热壁、冷壁和中面的 Nusselt 数及其极值。
!===========================================================================================================================


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

  write(*,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'u_mid_max =', umax*velocityScaleCompare, 'at y =', yAtU, 'z =', zAtU, 'on x_mid =', targetX

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) '--- ', trim(logTag), ' ---'
  write(00,*) 'x_mid =', targetX
  write(00,*) 'u_mid_max =', umax*velocityScaleCompare, ' y_pos =', yAtU, ' z_pos =', zAtU
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

  write(*,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'v_mid_max =', vmax*velocityScaleCompare, 'at x =', xAtV, 'z =', zAtV, 'on y_mid =', targetY

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) '--- ', trim(logTag), ' ---'
  write(00,*) 'y_mid =', targetY
  write(00,*) 'v_mid_max =', vmax*velocityScaleCompare, ' x_pos =', xAtV, ' z_pos =', zAtV
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
  use commondata3d
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
  use commondata3d
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4) :: i, j, uout, kL, kR
  real(kind=8) :: targetZ, weight, valU, valV, valW, valT

  ! 输出 z=Lz/2 中面，便于观察第三方向上的结构
  targetZ = 0.5d0 * zp(nz+1)
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



