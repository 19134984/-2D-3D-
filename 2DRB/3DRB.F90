!=============================================================
!!!    3D buoyancy-driven natural convection
!!!    D3Q19 flow + D3Q7 temperature
!!!    主循环按 collision/streaming/boundary/macro 展开
!!!    后处理也放在同一文件里，便于直接查阅
!=============================================================

! 宏开关：
! steadyFlow / unsteadyFlow : 稳态或非稳态运行模式
#define steadyFlow
!#define unsteadyFlow

! EnableRBInitPerturbation3D : 小 Ra 时给 3D RB 一个初始温度扰动，帮助脱离导热支路
! EnableUseG                : 是否启用温度方程中的偏差修正
#define EnableRBInitPerturbation3D
#define EnableUseG


module commondata3d
  implicit none

  ! qf / qt 分别表示流场和温度场的离散速度数
  integer(kind=4), parameter :: qf = 19
  integer(kind=4), parameter :: qt = 7

  ! 重启控制参数：
  ! loadInitField=0 表示从初值启动
  ! loadInitField=1 表示从 backupFile3D-*.bin 重启
  integer(kind=4), parameter :: loadInitField = 0
  integer(kind=4), parameter :: reloadDimensionlessTime = 0
  integer(kind=4), parameter :: reloadbinFileNum = 0


  ! 网格：2:1:1 的 3D periodic RB 几何设成 40x20x20
  integer(kind=4), parameter :: nx = 40, ny = 20, nz = 20
  real(kind=8), parameter :: lengthUnit = dble(ny)
  real(kind=8), parameter :: pi = acos(-1.0d0)
  integer(kind=4), parameter :: itc_max = 20000000
  real(kind=8),    parameter :: outputFrequency = 100.0d0
  ! 主要控制参数：Rayleigh / Prandtl / Mach / 热壁温度
  real(kind=8), parameter :: Rayleigh = 2.0d3
  real(kind=8), parameter :: Prandtl = 0.71d0
  real(kind=8), parameter :: Mach = 0.1d0
  real(kind=8), parameter :: Thot = 0.5d0
  real(kind=8), parameter :: Tcold = -0.5d0
  real(kind=8), parameter :: Tref = 0.5d0 * (Thot + Tcold)

  ! 松弛时间和输运系数：
  ! tauf 仍沿用流场的黏性定义
  ! 温度场这里保留 2DRB 的算法口径：
  ! 1) D3Q7 的变换矩阵 N 采用参考文献的正交基底
  ! 2) 但热扩散率仍用固定权重 + taug 的方式给定
  ! 3) 也就是“基底学文献，碰撞参数沿用 2DRB 思路”
  real(kind=8), parameter :: tauf = 0.5d0 + Mach * lengthUnit * dsqrt(3.0d0 * Prandtl / Rayleigh)
  real(kind=8), parameter :: viscosity = (tauf - 0.5d0) / 3.0d0
  real(kind=8), parameter :: diffusivity = viscosity / Prandtl
  real(kind=8), parameter :: cs2T = 0.25d0
  real(kind=8), parameter :: taug = 0.5d0 + diffusivity / cs2T
  real(kind=8), parameter :: s_eT = 1.0d0
  real(kind=8), parameter :: s_qT = 1.0d0
  real(kind=8), parameter :: s_jT = 1.0d0 / taug


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

  ! 体平均 Nu / Re 的时间序列缓存
  real(kind=8) :: NuVolAvg(0:dimensionlessTimeMax), ReVolAvg(0:dimensionlessTimeMax)

  ! 输出文件命名与 2D 版分开，避免共存时互相覆盖
  character(len=100) :: binFilePrefix = "buoyancyCavity3D"
  character(len=100) :: pltFilePrefix = "buoyancyCavity3D"
  character(len=100) :: reloadFilePrefix = "backupFile3D"
  character(len=100) :: settingsFile = "SimulationSettings3D.txt"

  real(kind=8) :: errorU, errorT

  ! 几何坐标数组：包含物理边界点 0 和 nx+1/ny+1/nz+1
  real(kind=8) :: xp(0:nx+1), yp(0:ny+1), zp(0:nz+1)
  ! 宏观场：u,v,w,T,rho
  real(kind=8), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:), T(:,:,:), rho(:,:,:)

#ifdef steadyFlow
  ! 稳态误差判据需要保存上一次输出时刻的场
  real(kind=8), allocatable :: up(:,:,:), vp(:,:,:), wp(:,:,:), Tp(:,:,:)
#endif

  ! f/g 是当前分布函数，f_post/g_post 是碰撞后、迁移前的分布函数，post 数组带 ghost 层便于直接 pull streaming
  real(kind=8), allocatable :: f(:,:,:,:), f_post(:,:,:,:)
  real(kind=8), allocatable :: g(:,:,:,:), g_post(:,:,:,:)
  ! 流场力项，以及温度方程中的历史热流项
  real(kind=8), allocatable :: Fx(:,:,:), Fy(:,:,:), Fz(:,:,:)
  real(kind=8), allocatable :: Bx_prev(:,:,:), By_prev(:,:,:), Bz_prev(:,:,:)

  integer(kind=4) :: itc

#ifdef EnableUseG
  logical, parameter :: useG = .true.
#else
  logical, parameter :: useG = .false.
#endif

  real(kind=8) :: Nu_global, Nu_hot, Nu_cold, Nu_middle

  ! D3Q19 / D3Q7 的离散速度、反向索引和权重
  integer(kind=4) :: ex(0:qf-1), ey(0:qf-1), ez(0:qf-1), opp(0:qf-1)
  real(kind=8)    :: omega(0:qf-1)

  integer(kind=4) :: exT(0:qt-1), eyT(0:qt-1), ezT(0:qt-1), oppT(0:qt-1)
  real(kind=8)    :: omegaT(0:qt-1)

  real(kind=8) :: M19(qf,qf), Minv19(qf,qf)
  real(kind=8) :: M7(qt,qt),  Minv7(qt,qt)
  ! M 和 Minv 在初始化时按当前正交基底构造，既保持文献编号一致，也方便后续自己改基底做试验
end module commondata3d


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

  ! 先写日志头，并固定 OpenMP 线程数
  ! 这里故意不引入 MPI / OpenACC，保持实现路径单纯
  open(unit=00, file=trim(settingsFile), status='replace')
  string = ctime(time())
  write(00,*) 'Start: ', string
  write(00,*) 'Starting OpenMP >>>>>>'
  call OMP_set_num_threads(24)
  myMaxThreads = OMP_get_max_threads()
  write(00,*) 'Max Running threads:', myMaxThreads
  close(00)

  call initial3d()

  call CPU_TIME(timeStart)
  timeStart2 = OMP_get_wtime()

  ! 主推进：
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
    ! 稳态模式下每隔一段步数检查一次相对误差，避免每步都做全场比较
    if (mod(itc, 2000) .EQ. 0) call check3d()
#endif

    ! 每到一个物理输出间隔，就统计一次 Nu/Re，并按需输出切面或重启文件
    if (mod(itc, outputIntervalItc) .EQ. 0) then
      call calNuRe3d()

#ifdef steadyFlow
      if ((outputPltFile .EQ. 1) .AND. (mod(itc, backupIntervalItc) .EQ. 0)) then
        call output_midplanes_tecplot3d()
      endif
#endif

#ifdef unsteadyFlow
      if (outputBinFile .EQ. 1) then
        call output_binary3d()
        if (mod(itc, backupIntervalItc) .EQ. 0) call backupData3d()
      endif
      if (outputPltFile .EQ. 1) call output_midplanes_tecplot3d()
#endif
    endif
  enddo

  call CPU_TIME(timeEnd)
  timeEnd2 = OMP_get_wtime()

#ifdef steadyFlow
  ! 稳态模式退出循环后，补一次最终场输出
  if (outputPltFile .EQ. 1) call output_midplanes_tecplot3d()
  if (outputBinFile .EQ. 1) call output_binary3d()
#endif

  ! 最终诊断：先输出全局 Nu，再输出壁面 Nu 和三中面速度极值
  call RBcalc_Nu_global3d()
  call RBcalc_Nu_wall_avg3d()
  call RBcalc_midplane_velocity_max3d()
  call calNuRe3d()

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
  write(00,*) 'useG =', useG
  if (outputBinFile .EQ. 1) then
    call backupData3d()
  endif
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


subroutine initial3d()
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

  ! 先初始化离散速度、权重，以及 MRT 变换矩阵
  call init_lattice_constants_3d()
  call init_mrt_matrices_3d()

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

  ! 把这次算例的主要信息写入日志文件，方便后续查设置
  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) '-------------------------------------------------------------------------------'
  write(00,*) 'I am 3D periodic Rayleigh-Benard Cell'
  write(00,*) 'Mesh:', nx, ny, nz
  write(00,*) 'Rayleigh=', real(Rayleigh,kind=8), '; Prandtl =', real(Prandtl,kind=8), '; Mach =', real(Mach,kind=8)
  write(00,*) 'Length unit: L0 =', real(lengthUnit,kind=8)
  write(00,*) 'Time unit: Sqrt(L0/(gBeta*DeltaT)) =', real(timeUnit,kind=8)
  write(00,*) 'Velocity unit: Sqrt(gBeta*L0*DeltaT) =', real(velocityUnit,kind=8)
  write(00,*) 'tauf =', real(tauf,kind=8), '; taug =', real(taug,kind=8)
  write(00,*) 'cs2T =', real(cs2T,kind=8), '; s_jT =', real(s_jT,kind=8), &
       '; s_eT =', real(s_eT,kind=8), '; s_qT =', real(s_qT,kind=8)
  write(00,*) 'viscosity =', real(viscosity,kind=8), '; diffusivity =', real(diffusivity,kind=8)
  write(00,*) 'outputFrequency =', real(outputFrequency,kind=8), ' free-fall time units'
  write(00,*) '......................  or ', outputIntervalItc, ' in itc units'
  write(00,*) 'backupInterval =', backupInterval, ' free-fall time units'
  write(00,*) '.................... or ', backupIntervalItc, ' in itc units'
  write(00,*) 'itc_max =', itc_max
  write(00,*) 'default epsU =', real(epsU,kind=8), '; epsT =', real(epsT,kind=8)
  write(00,*) 'OpenMP only; MPI/OpenACC are not included in this file'
  close(00)

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
    u = 0.0d0
    v = 0.0d0
    w = 0.0d0
    T = 0.0d0

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          T(i,j,k) = Thot + (yp(j) - yp(0)) / (yp(ny+1) - yp(0)) * (Tcold - Thot)
        enddo
      enddo
    enddo


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
      open(unit=00, file=trim(settingsFile), status='unknown', position='append')
      write(00,'(a,1x,es12.4)') '3D RB initial T perturbation amplitude =', rbInitPerturbAmp
      close(00)
    endif


    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ! 用当前宏观场直接生成 feq / geq，作为初始分布函数
          call compute_feq_d3q19(rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), feq)
          call compute_geq_d3q7(T(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), geq)
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
    write(reloadFileName,'(i0)') reloadbinFileNum
    open(unit=01, file=trim(reloadFilePrefix)//'-'//trim(adjustl(reloadFileName))//'.bin', &
         form='unformatted', access='sequential', status='old')
    read(01) ((((f(alpha,i,j,k), i=1,nx), j=1,ny), k=1,nz), alpha=0,qf-1)
    read(01) ((((g(alpha,i,j,k), i=1,nx), j=1,ny), k=1,nz), alpha=0,qt-1)
    close(01)
    call reconstruct_macro_from_fg3d()
  else
    open(unit=00, file=trim(settingsFile), status='unknown', position='append')
    write(00,*) 'Error: initial field is not properly set'
    close(00)
    stop
  endif

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

end subroutine initial3d


subroutine init_lattice_constants_3d()
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

  ! 这里温度权重按 2DRB 思路固定下来：w0 = 1/4，w1~w6 = 1/8，对应 cs2T = 1/4
  omegaT(0) = 1.0d0 / 4.0d0
  omegaT(1:6) = 1.0d0 / 8.0d0

end subroutine init_lattice_constants_3d


subroutine init_mrt_matrices_3d()
  use commondata3d
  implicit none

  ! M / Minv 在程序启动时一次性构造完，后续碰撞步骤直接拿来做矩空间变换
  ! 这里不再调用通用求逆，而是利用文献这套正交基底的已知模长，直接给出逆矩阵
  call build_basis_matrix_d3q19(M19)
  call init_inverse_matrix_d3q19()

  call build_basis_matrix_d3q7(M7)
  call init_inverse_matrix_d3q7()

end subroutine init_mrt_matrices_3d


subroutine init_inverse_matrix_d3q19()
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
  do i = 1, qf
    do j = 1, qf
      Minv19(j,i) = M19(i,j) * invRowNorm2(i)
    enddo
  enddo

end subroutine init_inverse_matrix_d3q19


subroutine init_inverse_matrix_d3q7()
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

end subroutine init_inverse_matrix_d3q7


subroutine build_basis_matrix_d3q19(M)
  use commondata3d
  implicit none

  real(kind=8), intent(out) :: M(qf,qf)
  integer(kind=4) :: alpha
  real(kind=8) :: exa, eya, eza, ex2, ey2, ez2, e2, e4

  ! 这里直接按文献给出的 D3Q19 正交基底来生成 M19，和文献中的速度编号、矩空间定义保持一一对应
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

end subroutine build_basis_matrix_d3q19


subroutine build_basis_matrix_d3q7(M)
  use commondata3d
  implicit none

  real(kind=8), intent(out) :: M(qt,qt)
  integer(kind=4) :: alpha
  real(kind=8) :: exa, eya, eza, ex2, ey2, ez2, e2

  ! 这里按文献 Eq. (16) 的 7 个基底逐行生成 M7
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

end subroutine build_basis_matrix_d3q7











subroutine compute_feq_d3q19(rhoLoc, uLoc, vLoc, wLoc, feq)
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

end subroutine compute_feq_d3q19


subroutine compute_geq_d3q7(TLoc, uLoc, vLoc, wLoc, geq)
  use commondata3d
  implicit none

  real(kind=8), intent(in) :: TLoc, uLoc, vLoc, wLoc
  real(kind=8), intent(out) :: geq(qt)

  integer(kind=4) :: alpha
  real(kind=8) :: eu

  ! 温度场采用线性平衡分布，结构和 2D 版保持一致
  do alpha = 0, qt-1
    eu = dble(exT(alpha)) * uLoc + dble(eyT(alpha)) * vLoc + dble(ezT(alpha)) * wLoc
    geq(alpha+1) = omegaT(alpha) * TLoc * (1.0d0 + eu / cs2T)
  enddo

end subroutine compute_geq_d3q7




subroutine collision3d()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k, alpha
  real(kind=8) :: fLoc(qf), m(qf), meq(qf), mPost(qf), fPostLoc(qf)
  real(kind=8) :: rhoLoc, uLoc, vLoc, wLoc, u2, uDotF
  real(kind=8) :: FxLoc, FyLoc, FzLoc
  real(kind=8), parameter :: s_e  = 1.0d0 / tauf
  real(kind=8), parameter :: s_eps = 1.0d0 / tauf
  real(kind=8), parameter :: s_nu = 1.0d0 / tauf
  real(kind=8), parameter :: s_pi = 1.0d0 / tauf
  real(kind=8), parameter :: s_q  = 8.0d0 * (2.0d0 * tauf - 1.0d0) / (8.0d0 * tauf - 1.0d0)
  real(kind=8), parameter :: s_m  = 8.0d0 * (2.0d0 * tauf - 1.0d0) / (8.0d0 * tauf - 1.0d0)

  ! 这里改回和 2DRB 同类型的矩空间碰撞：
  ! 先做 m = M*f，再按文献的平衡矩 meq 和松弛率逐个碰撞，
  ! 体力项也先投到矩空间，再乘以 (I-S/2) 做半步修正。

  !$omp parallel do collapse(3) default(none) shared(f,f_post,rho,u,v,w,T,Fx,Fy,Fz,M19,Minv19) &
  !$omp private(i,j,k,alpha,fLoc,m,meq,mPost,fPostLoc,rhoLoc,uLoc,vLoc,wLoc,u2,uDotF,FxLoc,FyLoc,FzLoc)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        do alpha = 0, qf-1
          fLoc(alpha+1) = f(alpha,i,j,k)
        enddo

        rhoLoc = rho(i,j,k)
        uLoc = u(i,j,k)
        vLoc = v(i,j,k)
        wLoc = w(i,j,k)
        u2 = uLoc * uLoc + vLoc * vLoc + wLoc * wLoc

        FxLoc = 0.0d0
        FyLoc = rhoLoc * gBeta * (T(i,j,k) - Tref)
        FzLoc = 0.0d0
        Fx(i,j,k) = FxLoc
        Fy(i,j,k) = FyLoc
        Fz(i,j,k) = FzLoc
        uDotF = uLoc * FxLoc + vLoc * FyLoc + wLoc * FzLoc

        ! x / z 方向无体力，RB 浮力只在 y 方向

        ! 先把分布函数变到矩空间，再做松弛
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

        ! 矩空间碰撞加体力修正：结构上和 2DRB 一样，只是这里使用 D3Q19 的矩定义
        mPost(1)  = m(1)
        mPost(2)  = m(2)  - s_e   * (m(2)  - meq(2))  + (1.0d0 - 0.5d0 * s_e  ) * 38.0d0 * uDotF
        mPost(3)  = m(3)  - s_eps * (m(3)  - meq(3))  + (1.0d0 - 0.5d0 * s_eps) * (-11.0d0) * uDotF
        mPost(4)  = m(4)  + FxLoc
        mPost(5)  = m(5)  - s_q   * (m(5)  - meq(5))  + (1.0d0 - 0.5d0 * s_q  ) * (-2.0d0 / 3.0d0) * FxLoc
        mPost(6)  = m(6)  + FyLoc
        mPost(7)  = m(7)  - s_q   * (m(7)  - meq(7))  + (1.0d0 - 0.5d0 * s_q  ) * (-2.0d0 / 3.0d0) * FyLoc
        mPost(8)  = m(8)  + FzLoc
        mPost(9)  = m(9)  - s_q   * (m(9)  - meq(9))  + (1.0d0 - 0.5d0 * s_q  ) * (-2.0d0 / 3.0d0) * FzLoc
        mPost(10) = m(10) - s_nu  * (m(10) - meq(10)) + &
             (1.0d0 - 0.5d0 * s_nu ) * &
             (4.0d0 * uLoc * FxLoc - 2.0d0 * vLoc * FyLoc - 2.0d0 * wLoc * FzLoc)
        mPost(11) = m(11) - s_pi  * (m(11) - meq(11)) + &
             (1.0d0 - 0.5d0 * s_pi ) * &
             (-2.0d0 * uLoc * FxLoc + vLoc * FyLoc + wLoc * FzLoc)
        mPost(12) = m(12) - s_nu  * (m(12) - meq(12)) + &
             (1.0d0 - 0.5d0 * s_nu ) * &
             (2.0d0 * vLoc * FyLoc - 2.0d0 * wLoc * FzLoc)
        mPost(13) = m(13) - s_pi  * (m(13) - meq(13)) + (1.0d0 - 0.5d0 * s_pi ) * (-vLoc * FyLoc + wLoc * FzLoc)
        mPost(14) = m(14) - s_nu  * (m(14) - meq(14)) + (1.0d0 - 0.5d0 * s_nu ) * (uLoc * FyLoc + vLoc * FxLoc)
        mPost(15) = m(15) - s_nu  * (m(15) - meq(15)) + (1.0d0 - 0.5d0 * s_nu ) * (vLoc * FzLoc + wLoc * FyLoc)
        mPost(16) = m(16) - s_nu  * (m(16) - meq(16)) + (1.0d0 - 0.5d0 * s_nu ) * (uLoc * FzLoc + wLoc * FxLoc)
        mPost(17) = m(17) - s_m   * m(17)
        mPost(18) = m(18) - s_m   * m(18)
        mPost(19) = m(19) - s_m   * m(19)

        ! 再由逆矩阵回到速度空间
        fPostLoc = matmul(Minv19, mPost)

        do alpha = 0, qf-1
          f_post(alpha,i,j,k) = fPostLoc(alpha+1)
        enddo
      enddo
    enddo
  enddo
  !$omp end parallel do

end subroutine collision3d


subroutine fill_periodic_ghosts_f_post()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k, alpha

  ! x 和 z 方向都是周期边界，这里先把碰撞后分布写入 ghost 层，随后 streaming3d 直接 pull 即可
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

end subroutine fill_periodic_ghosts_f_post


subroutine streaming3d()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k, ip, jp, kp, alpha

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

end subroutine streaming3d


subroutine bounceback3d()
  use commondata3d
  implicit none

  integer(kind=4) :: i, k, alpha

  ! 这里只处理 y=1 和 y=ny 两个水平壁面的无滑移反弹
  ! x/z 方向由于本来就是周期边界，所以不在这里处理
  !$omp parallel do collapse(2) default(none) shared(f,f_post,ey,opp) private(i,k,alpha)
  do k = 1, nz
    do i = 1, nx
      do alpha = 0, qf-1
        if (ey(alpha) .EQ. 1) f(alpha,i,1,k) = f_post(opp(alpha),i,1,k)
        if (ey(alpha) .EQ. -1) f(alpha,i,ny,k) = f_post(opp(alpha),i,ny,k)
      enddo
    enddo
  enddo
  !$omp end parallel do

end subroutine bounceback3d


subroutine macro3d()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k, alpha
  real(kind=8) :: momx, momy, momz, rhoLoc

  ! 由分布函数恢复宏观密度和三分量速度
  ! 速度恢复里带 0.5F 的半步修正，对应前面的体力项半步离散
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

end subroutine macro3d


subroutine collisionT3d()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k, alpha
  real(kind=8) :: gLoc(qt), n(qt), neq(qt), nPost(qt), gPostLoc(qt)
  real(kind=8) :: Bx, By, Bz, dBx, dBy, dBz
  real(kind=8), parameter :: SG = 1.0d0 - 0.5d0 * s_jT

  ! 温度场的思路和 2D 一样：
  ! 1) 先构造当前热通量 B = uT
  ! 2) 若启用 useG，则用相邻时刻差分 dB 做修正；3) 在矩空间里做温度碰撞
  !$omp parallel do collapse(3) default(none) shared(g,g_post,u,v,w,T,Bx_prev,By_prev,Bz_prev,M7,Minv7) &
  !$omp private(i,j,k,alpha,gLoc,n,neq,nPost,gPostLoc,Bx,By,Bz,dBx,dBy,dBz)
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
        else
          dBx = 0.0d0
          dBy = 0.0d0
          dBz = 0.0d0
        endif

        Bx_prev(i,j,k) = Bx
        By_prev(i,j,k) = By
        Bz_prev(i,j,k) = Bz

        ! 先把温度分布函数变到矩空间
        n = matmul(M7, gLoc)

        ! 这里按 2DRB 的思路直接写平衡矩，而不是跟文献那套 aT/qkappa 走
        ! [T, uT, vT, wT, -3/4*T, 0, 0]
        neq(1) = T(i,j,k)
        neq(2) = Bx
        neq(3) = By
        neq(4) = Bz
        neq(5) = -0.75d0 * T(i,j,k)
        neq(6) = 0.0d0
        neq(7) = 0.0d0

        nPost(1) = n(1)
        nPost(2:4) = n(2:4) - s_jT * (n(2:4) - neq(2:4))
        nPost(5) = n(5) - s_eT * (n(5) - neq(5))
        nPost(6:7) = n(6:7) - s_qT * (n(6:7) - neq(6:7))
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

end subroutine collisionT3d


subroutine fill_periodic_ghosts_g_post()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k, alpha

  ! 温度场和流场一样，x/z 方向先填周期 ghost 层
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

end subroutine fill_periodic_ghosts_g_post


subroutine streamingT3d()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k, ip, jp, kp, alpha

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

end subroutine streamingT3d


subroutine bouncebackT3d()
  use commondata3d
  implicit none

  integer(kind=4) :: i, k, alpha

  ! 温度边界只在上下水平壁处理：
  ! y=1   为热壁 Thot
  ! y=ny  为冷壁 Tcold
  ! 写法沿用 2D 版“反弹 + 壁温修正”思路
  !$omp parallel do collapse(2) default(none) shared(g,g_post,eyT,oppT,omegaT) private(i,k,alpha)
  do k = 1, nz
    do i = 1, nx
      do alpha = 0, qt-1
        if (eyT(alpha) .EQ. 1) g(alpha,i,1,k) = -g_post(oppT(alpha),i,1,k) + 2.0d0 * omegaT(alpha) * Thot
        if (eyT(alpha) .EQ. -1) g(alpha,i,ny,k) = -g_post(oppT(alpha),i,ny,k) + 2.0d0 * omegaT(alpha) * Tcold
      enddo
    enddo
  enddo
  !$omp end parallel do

end subroutine bouncebackT3d


subroutine macroT3d()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k, alpha
  real(kind=8) :: TLoc

  ! 温度恢复就是对 7 个方向的 g 求和
  !$omp parallel do collapse(3) default(none) shared(g,T) private(i,j,k,alpha,TLoc)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        TLoc = 0.0d0
        do alpha = 0, qt-1
          TLoc = TLoc + g(alpha,i,j,k)
        enddo
        T(i,j,k) = TLoc
      enddo
    enddo
  enddo
  !$omp end parallel do

end subroutine macroT3d


subroutine reconstruct_macro_from_fg3d()
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k, alpha, iter
  real(kind=8) :: momx, momy, momz
  logical :: rho_bad

  ! 重启时只存了 f 和 g，所以这里统一把 T、rho、u、v、w 以及历史热流都重构回来
  call macroT3d()
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

        Bx_prev(i,j,k) = u(i,j,k) * T(i,j,k)
        By_prev(i,j,k) = v(i,j,k) * T(i,j,k)
        Bz_prev(i,j,k) = w(i,j,k) * T(i,j,k)
      enddo
    enddo
  enddo
  !$omp end parallel do

  if (rho_bad) stop 'Warning: non-positive rho found during restart reconstruction.'

end subroutine reconstruct_macro_from_fg3d


#ifdef steadyFlow
subroutine check3d()
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

  call append_convergence_tecplot3d('convergence3D.plt', itc, errorU, errorT)
  write(caseTag,'("Ra=",ES10.3E2,",nx=",I0,",ny=",I0,",nz=",I0,",useG=",L1)') Rayleigh, nx, ny, nz, useG
  call append_convergence_master_tecplot3d('convergence_all_3D.plt', caseTag, itc, errorU, errorT)
  write(*,'(I12,1X,ES24.16,1X,ES24.16)') itc, errorU, errorT

end subroutine check3d
#endif


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


subroutine output_binary3d()
  use commondata3d
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


subroutine backupData3d()
  use commondata3d
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
  open(unit=05, file='backupFile3D-'//trim(filename)//'.bin', form='unformatted', access='sequential')
  write(05) ((((f(alpha,i,j,k), i=1,nx), j=1,ny), k=1,nz), alpha=0,qf-1)
  write(05) ((((g(alpha,i,j,k), i=1,nx), j=1,ny), k=1,nz), alpha=0,qt-1)
  close(05)

  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,*) 'Backup f and g to the file: backupFile3D-', trim(filename), '.bin'
  close(00)

end subroutine backupData3d


subroutine output_midplanes_tecplot3d()
  use commondata3d
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
  use commondata3d
  implicit none

  integer(kind=4) :: i, j, k
  real(kind=8) :: NuVolAvg_temp, ReVolAvg_temp

  ! 这里记录的是时间序列版本的体平均 Nu / Re：
  ! NuVolAvg : 1 + (L/kappa) * <vT>
  ! ReVolAvg : 全域 RMS 速度对应的 Reynolds 数
  if (dimensionlessTime .GE. dimensionlessTimeMax) return
  dimensionlessTime = dimensionlessTime + 1

  NuVolAvg_temp = 0.0d0
  !$omp parallel do collapse(3) default(none) shared(v,T) private(i,j,k) reduction(+:NuVolAvg_temp)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        NuVolAvg_temp = NuVolAvg_temp + v(i,j,k) * T(i,j,k)
      enddo
    enddo
  enddo
  !$omp end parallel do
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

end subroutine calNuRe3d


subroutine RBcalc_Nu_global3d()
  use commondata3d
  implicit none
  integer(kind=4) :: i, j, k
  real(kind=8) :: Nu_sum

  ! 最终全局 Nu：对整个体积做平均
  Nu_sum = 0.0d0
  !$omp parallel do collapse(3) default(none) shared(v,T) private(i,j,k) reduction(+:Nu_sum)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        Nu_sum = Nu_sum + v(i,j,k) * T(i,j,k)
      enddo
    enddo
  enddo
  !$omp end parallel do

  Nu_global = 1.0d0 + heatFluxScale * Nu_sum / dble(nx * ny * nz)
  write(*,'(a,1x,es16.8)') 'Nu_global =', Nu_global
  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,'(a,1x,es16.8)') 'Nu_global =', Nu_global
  close(00)

end subroutine RBcalc_Nu_global3d


subroutine RBcalc_Nu_wall_avg3d()
  use commondata3d
  implicit none
  integer(kind=4) :: i, jB, jT, k
  real(kind=8) :: dy, deltaT, coef
  real(kind=8) :: sum_hot, sum_cold, sum_mid, qy_wall

  ! 这一版先做“平面平均”的 3D RB 后处理：
  ! Nu_hot    : 底壁 x-z 平面平均
  ! Nu_cold   : 顶壁 x-z 平面平均
  ! Nu_middle : y=1/2 中平面平均
  dy = 1.0d0 / lengthUnit
  deltaT = Thot - Tcold
  coef = heatFluxScale

  sum_hot = 0.0d0
  !$omp parallel do collapse(2) default(none) shared(T,dy,deltaT) private(i,k,qy_wall) reduction(+:sum_hot)
  do k = 1, nz
    do i = 1, nx
      qy_wall = 2.0d0 * (Thot - T(i,1,k)) / dy
      sum_hot = sum_hot + qy_wall / deltaT
    enddo
  enddo
  !$omp end parallel do
  Nu_hot = sum_hot / dble(nx * nz)

  sum_cold = 0.0d0
  !$omp parallel do collapse(2) default(none) shared(T,dy,deltaT) private(i,k,qy_wall) reduction(+:sum_cold)
  do k = 1, nz
    do i = 1, nx
      qy_wall = 2.0d0 * (T(i,ny,k) - Tcold) / dy
      sum_cold = sum_cold + qy_wall / deltaT
    enddo
  enddo
  !$omp end parallel do
  Nu_cold = sum_cold / dble(nx * nz)

  sum_mid = 0.0d0
  if (mod(ny,2) .EQ. 1) then
    jB = (ny + 1) / 2
    !$omp parallel do collapse(2) default(none) shared(v,T,jB,dy,deltaT,coef) private(i,k) reduction(+:sum_mid)
    do k = 1, nz
      do i = 1, nx
        sum_mid = sum_mid + &
             (coef * v(i,jB,k) * (T(i,jB,k) - Tref) - &
             (T(i,jB+1,k) - T(i,jB-1,k)) / (2.0d0 * dy)) / deltaT
      enddo
    enddo
    !$omp end parallel do
  else
    jB = ny / 2
    jT = jB + 1
    !$omp parallel do collapse(2) default(none) shared(v,T,jB,jT,dy,deltaT,coef) private(i,k) reduction(+:sum_mid)
    do k = 1, nz
      do i = 1, nx
        sum_mid = sum_mid + (coef * 0.5d0 * (v(i,jB,k)*(T(i,jB,k)-Tref) + v(i,jT,k)*(T(i,jT,k)-Tref)) &
                 + (T(i,jB,k) - T(i,jT,k)) / dy) / deltaT
      enddo
    enddo
    !$omp end parallel do
  endif
  Nu_middle = sum_mid / dble(nx * nz)

  write(*,'(a,1x,es16.8)') 'Nu_hot(bottom) =', Nu_hot
  write(*,'(a,1x,es16.8)') 'Nu_cold(top)   =', Nu_cold
  write(*,'(a,1x,es16.8)') 'Nu_middle      =', Nu_middle
  open(unit=00, file=trim(settingsFile), status='unknown', position='append')
  write(00,'(a,1x,es16.8)') 'Nu_hot(bottom) =', Nu_hot
  write(00,'(a,1x,es16.8)') 'Nu_cold(top)   =', Nu_cold
  write(00,'(a,1x,es16.8)') 'Nu_middle      =', Nu_middle
  close(00)

end subroutine RBcalc_Nu_wall_avg3d


subroutine RBcalc_midplane_velocity_max3d()
  use commondata3d
  implicit none
  integer(kind=4) :: i, j, k, iL, iR, jL, jR, kL, kR, jBest, kBest, iBest
  real(kind=8) :: targetX, targetY, targetZ, wx, wy, wz, val, umax, vmax, wmax
  real(kind=8) :: yAtU, zAtU, xAtV, zAtV, xAtW, yAtW

  ! u_max : 在 x = Lx/2 的中面上找最大 u，输出 (y,z)
  ! v_max : 在 y = Ly/2 的中面上找最大 v，输出 (x,z)
  ! w_max : 在 z = Lz/2 的中面上找最大 w，输出 (x,y)
  targetX = 0.5d0 * xp(nx+1)
  targetY = 0.5d0 * yp(ny+1)
  targetZ = 0.5d0 * zp(nz+1)
  call find_bracketing_index(xp, nx, targetX, iL, iR, wx)
  call find_bracketing_index(yp, ny, targetY, jL, jR, wy)
  call find_bracketing_index(zp, nz, targetZ, kL, kR, wz)

  umax = -huge(1.0d0); jBest = 1; kBest = 1
  do k = 1, nz
    do j = 1, ny
      call interp_scalar_x(iL, iR, wx, j, k, u, val)
      if (val .GT. umax) then
        umax = val; jBest = j; kBest = k
      endif
    enddo
  enddo
  yAtU = yp(jBest); zAtU = zp(kBest)

  vmax = -huge(1.0d0); iBest = 1; kBest = 1
  do k = 1, nz
    do i = 1, nx
      call interp_scalar_y(jL, jR, wy, i, k, v, val)
      if (val .GT. vmax) then
        vmax = val; iBest = i; kBest = k
      endif
    enddo
  enddo
  xAtV = xp(iBest); zAtV = zp(kBest)

  wmax = -huge(1.0d0); iBest = 1; jBest = 1
  do j = 1, ny
    do i = 1, nx
      call interp_scalar_z(kL, kR, wz, i, j, w, val)
      if (val .GT. wmax) then
        wmax = val; iBest = i; jBest = j
      endif
    enddo
  enddo
  xAtW = xp(iBest); yAtW = yp(jBest)

  write(*,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'u_max* =', umax*velocityScaleCompare, 'y =', yAtU, 'z =', zAtU, &
       'x_mid =', targetX
  write(*,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'v_max* =', vmax*velocityScaleCompare, 'x =', xAtV, 'z =', zAtV, &
       'y_mid =', targetY
  write(*,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'w_max* =', wmax*velocityScaleCompare, 'x =', xAtW, 'y =', yAtW, &
       'z_mid =', targetZ

end subroutine RBcalc_midplane_velocity_max3d


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


