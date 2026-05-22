module commondata
        use mpi
        implicit none

        ! Grid parameters
        integer(kind=4), parameter :: nx=513, ny=nx
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1

        real(kind=8), parameter :: Ra=1e5          !Rayleigh number
        real(kind=8), parameter :: Ma=0.1d0         !Mach number
        real(kind=8), parameter :: Pr=0.71d0        !Prandtl number
        real(kind=8), parameter :: Thot=1.0d0, Tcold=0.0d0, Tref=0.0d0
        real(kind=8) :: kappa                        !Thermal expansion coefficient
        real(kind=8) :: gbeta                        !Volume expansion coefficient * Gravitational acceleration

        ! Parameters in 0 system
        real(kind=8), parameter :: rho0 = 1.0d0
        real(kind=8) :: length0                   ! Characteristic length
        real(kind=8) :: viscosity0                ! Kinematic viscosity
        real(kind=8) :: dt0                       ! Time step
        real(kind=8) :: dx0                       ! Grid spacing

        ! Parameters in LB system
        real(kind=8) :: length_LB                 ! Characteristic length
        real(kind=8) :: viscosity_LB              ! Kinematic viscosity
        real(kind=8) :: tauf, taut                ! Relaxation time
        real(kind=8) :: dt                       ! Time step

        ! Iteration control
        integer(kind=4) :: itc = 0            ! Current iteration count
        integer(kind=4), parameter :: itc_max = 20000  ! Maximum iterations

        ! Convergence criteria
        real(kind=8) :: errorU, errorT              ! Current error
        real(kind=8), parameter :: epsU=1e-6, epsT=1e-6  ! Convergence threshold

        ! Grid coordinates
        real(kind=8) :: xGrid(0:nx+1), y(0:ny+1)
        real(kind=8), allocatable :: yGrid(:)

        ! Flow variables
        real(kind=8), allocatable :: u(:,:), v(:,:)      ! Velocity components
        real(kind=8), allocatable :: temp(:,:)           ! temperature
        real(kind=8), allocatable :: rho(:,:)            ! Density
        real(kind=8), allocatable :: up(:,:), vp(:,:)    ! Previous velocity components for error checking
        real(kind=8), allocatable :: utemp(:,:)          ! Previous temperature for error checking
        real(kind=8), allocatable :: Fx(:,:), Fy(:,:)    ! force components
        real(8), allocatable :: u_all(:,:), v_all(:,:), rho_all(:,:), T_all(:,:)

        ! Distribution functions
        real(kind=8), allocatable :: f(:,:,:), f_post(:,:,:)  ! Current and post-collision distributions
        real(kind=8),allocatable :: g(:,:,:), g_post(:,:,:)

        ! MRT relaxation parameters
        real(kind=8) :: omega_U(0:8), omega_T(0:4)  ! Relaxation rates for MRT

        ! Lattice directions
        integer(kind=4) :: ex(0:8), ey(0:8)  ! Lattice velocity directions
        data ex/0, 1, 0, -1,  0, 1, -1, -1,  1/
        data ey/0, 0, 1,  0, -1, 1,  1, -1, -1/

        ! Additional MRT parameters
        real(kind=8) :: Snu, Sq, sig_k

        integer(kind=4), allocatable :: inter_x(:,:), inter_y(:,:)
!----------------------------MPI---------------------------------------------------
        integer(kind=4) :: NPROC, MYID, IERR,ISTATUS(MPI_STATUS_SIZE)
        integer(kind=4) :: nyLocal
        integer(kind=4) :: upid, downid
        integer(kind=4), allocatable :: start1d(:), end1d(:), count1d(:), displ1d(:)
        integer(kind=4), allocatable :: start2d(:), end2d(:), count2d(:), displ2d(:)
end module commondata

subroutine mesh()
    use commondata
    implicit none

    integer(kind=4) :: i, j
    real(kind=8) :: dx(nx+1), dy(ny+1)
    real(kind=8) :: constA
    allocate(yGrid(0:nyLocal+1))

    if(MYID == 0)then
        ! Compute grid coordinates
        constA = 3.2d0
        do i = 0, nx+1
        !原始 erf 范围：[-erf(0.5A), +erf(0.5A)]
        !除以 erf(0.5A) 后：[-1, +1]
        !加 1 后：[0, 2]
        !乘 0.5 后：[0, 1]
        !中间坐标跳得大 -> 网格间距大 -> 中间较疏
        !靠近两端坐标跳得小 -> 网格间距小 -> 边界附近较密
            xGrid(i) = 0.5d0 * (erf(constA  * (dble(i) / dble(nx+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
        end do
        do j = 0, ny+1
            y(j) = 0.5d0 * (erf(constA  * (dble(j) / dble(ny+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
        end do

        ! 由相邻坐标点相减得到非均匀网格的局部间距。
        ! dx(k) = xGrid(k)-xGrid(k-1)，dy(k) = y(k)-y(k-1)。
        dx(1:nx+1) = xGrid(1:nx+1) - xGrid(0:nx)
        dy(1:ny+1) = y(1:ny+1) - y(0:ny)
        ! 取第一段 x 向网格间距作为 0-system 中的参考空间步长。
        dx0 = dx(1)
        ! 这里把 0-system 的参考时间步长设为 dx0。
        dt0 = dx0
        write(*,*) "nx =", nx, ", ny =", ny
        write(*,*) "Ra =", real(Ra)
        write(*,*) "  "
        write(*,*) "---in 0 system---"
        write(*,*) "deltaX = ", dx0
        write(*,*) "deltaT = ", dt0

        ! 这里的 xGrid(0) 和 xGrid(nx+1) 是归一化几何辅助端点 0 和 1。
        ! 本参考代码没有把完整几何长度 1 直接作为后续 Ma/Ra 公式中的特征长度，
        ! 而是采用 length0 = 1-dx0 的有效长度口径。
        ! 这个口径通常和近壁半步/有效边界距离有关；若按物理边界位于
        ! x=0.5*dx0 与 x=1-0.5*dx0 理解，则两侧半步合计正好扣掉一个 dx0。
        ! 由于 xGrid(1)=dx0、xGrid(nx+1)=1，所以 length0 = xGrid(nx+1)-xGrid(1) = 1-dx0。
        length0 = xGrid(nx+1)-xGrid(1)
        write(*,*) "length0 = ", real(length0)
        write(*,*) "  "

        ! Calculate viscosity based on LB unit
        length_LB = 1.0d0 / dx(1)
        dt = dt0 * length_LB

        ! 下面这两行把内部坐标整体减去半个近壁参考步长 dx0/2，
        ! 可以理解为把左/下物理边界作为距离坐标的零点。
        ! 这样第一层内部点从 xGrid(1)=dx0 变为 dx0/2，换算到 LB 单位后约为 0.5，
        ! 表示第一层计算点距离物理边界半个 lattice unit。
        ! 该坐标会被后续 interpolate() 和输出使用，既是存储/输出坐标，也是插值迁移的几何坐标。
        ! 注意这里统一使用近壁参考步长 dx0/2 做整体偏移，而不是每一段非均匀局部间距的半值。
        ! y 方向使用同样的非均匀映射和参考步长，因此做相同的半步平移。
        xGrid(1:nx) = xGrid(1:nx)-dx0/2.0d0
        y(1:ny) = y(1:ny)-dx0/2.0d0
        ! 最后把 0-system 坐标放大到 LB 单位；此时 dx0 对应 1 个 lattice unit。
        xGrid=xGrid*length_LB
        y=y*length_LB
    endif
      !0 system:
      !腔体坐标约在 [0,1]
      !dx0 是归一化几何坐标下的参考网格间距
      !dt0 = dx0

      !LB system:
      !用 length_LB = 1/dx0 把坐标放大
      !使第一段网格间距变成 1
      !dt = dt0 * length_LB，通常也变成 1


    ! 广播 LB 时间步 dt：
    ! dt 是发送/接收缓冲区；1 表示发送 1 个 real(kind=8) 标量；
    ! MPI_REAL8 是数据类型；0 表示由 rank 0 作为 root 发出；
    ! MPI_COMM_WORLD 表示所有 rank 都参与；IERR 接收 MPI 调用错误码。
    call MPI_BCAST(dt,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    ! length_LB 是从 0-system 转到 LB system 的长度缩放因子，所有 rank 必须一致。
    call MPI_BCAST(length_LB,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    ! x 方向没有被 MPI 切分，所以把完整的 xGrid(1:nx) 广播给所有 rank。
    call MPI_BCAST(xGrid(1),nx,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    ! y 方向被一维切分：rank 0 按 count1d/displ1d 把全局 y 坐标分发到各 rank 的 yGrid(1:nyLocal)。
    ! 参数含义：y(1)=root 端发送缓冲区起点，count1d=每个 rank 接收的 y 坐标个数，
    !           displ1d=每个 rank 在 y(1:ny) 中的发送偏移量，MPI_REAL8=发送数据类型，
    !           yGrid(1)=本 rank 接收缓冲区起点，count1d(MYID)=本 rank 实际接收个数，
    !           MPI_REAL8=接收数据类型，0=root rank，MPI_COMM_WORLD=通信域，IERR=错误码。
    call MPI_SCATTERV(y(1),count1d,displ1d,MPI_REAL8,yGrid(1),count1d(MYID),MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    ! 这里只交换一维坐标数组 yGrid 的 halo，因此每次只发送 1 个坐标值，不是整条流场边界。
    ! 把本 rank 的第一层 yGrid(1) 发给下方邻居，同时从上方邻居接收其边界到 yGrid(nyLocal+1)。
    ! 真正的分布函数边界交换在后面的 update() 中完成，那里会发送 9*nx 或 5*nx 个数据。
    ! 参数含义：yGrid(1)=发送缓冲区起点，1=发送 1 个坐标，MPI_REAL8=发送类型，
    !           downid=发送目标 rank，0=发送标签，yGrid(nyLocal+1)=接收缓冲区起点，
    !           1=接收 1 个坐标，MPI_REAL8=接收类型，upid=接收来源 rank，0=接收标签，
    !           MPI_COMM_WORLD=通信域，ISTATUS=接收状态，IERR=错误码。
    call MPI_SENDRECV(yGrid(1),1,MPI_REAL8,downid,0,yGrid(nyLocal+1),1,MPI_REAL8,upid,0&
    ,MPI_COMM_WORLD,ISTATUS,IERR)

    ! 交换 yGrid 的下侧 halo：把本 rank 的最后一层 yGrid(nyLocal) 发给上方邻居，同时从下方邻居接收其边界到 yGrid(0)。
    call MPI_SENDRECV(yGrid(nyLocal),1,MPI_REAL8,upid,1,yGrid(0),1,MPI_REAL8,downid,1&
    ,MPI_COMM_WORLD,ISTATUS,IERR)

    return
end subroutine

subroutine initial()
    use commondata
    implicit none

    integer(kind=4) :: i, j, alpha
    real(kind=8) :: un(0:8)
    real(kind=8) :: velocitySquare
    allocate(inter_x(nx,3))
    allocate(inter_y(nyLocal,3))

    ! Initialize iteration count and error
    itc = 0
    errorU = 100.0d0
    errorT = 100.0d0

    ! Calculate viscosity based on LB unit
    ! 这里的 Ma-L-Ra-Pr 关系可理解为 viscosity_LB = Ma * L * sqrt(Pr) / sqrt(3*Ra)。
    ! 核心是特征长度 L 的换算：mesh() 中 length0 = xGrid(nx+1)-xGrid(1)。
    ! 这是半步边界处理下的有效长度：两侧各扣半个近壁步长，总共扣掉 dx(1)。
    ! 因为 xGrid(nx+1)=1，xGrid(1)=dx(1)，所以 length0 = 1-dx(1)。
    ! 转到 LB 单位后：length0/dx(1) = (1-dx(1))/dx(1) = 1/dx(1)-1 = length_LB-1。
    ! 因此这里用 length_LB-1.0d0 作为该参考代码的有效内部特征长度。
    viscosity_LB = (Ma*(length_LB-1.0d0)*dsqrt(Pr))/dsqrt(3.0d0*Ra)
    kappa = viscosity_LB/Pr
    gbeta = Ra*viscosity_LB*kappa/((length_LB-1.0d0)**3)

    ! Calculate relaxation time
    tauf = viscosity_LB * 3.0d0 + 0.5d0

    ! Calculate MRT relaxation parameters
    Snu = 1.0d0/tauf
    Sq = 8.0d0*(2.0d0-Snu)/(8.0d0-Snu)
    ! sig_k 是温度 D2Q5 MRT 模型中热扩散相关矩量的松弛率，可理解为 1/tau_g。
    ! 由 Snu=1/tauf 且 tauf-0.5=3*viscosity_LB，可得
    ! sig_k = 1/[0.5 + 4*viscosity_LB/Pr] = 1/[0.5 + 4*kappa]。
    ! 即 D2Q5 温度模型采用 kappa=(tau_g-0.5)/4，因此 sig_k=1/tau_g。
    sig_k = 1.0d0/(0.5d0+4.0d0*(1.0d0/Snu-0.5d0)/(3.0d0*Pr))

    if(MYID == 0)then
        ! Output initial parameters for verification
        write(*,*) "---in system LB"
        write(*,*) "deltaX = ", xGrid(1)-xGrid(0)
        write(*,*) "deltaT = ", dt
        write(*,*) "characteristic length   =", real(length_LB), "l.u."
        write(*,*) "viscosity_LB =", real(viscosity_LB), "l.u.^2/t.s."
        ! 这个比值约等于 dt_uniform/dt_nonuniform；非均匀网格近壁加密后 dx(1) 更小，
        ! 因此用 dx(1) 定义的时间步更小。例如 length_LB/nx=1.6 时，
        ! 当前非均匀网格时间步约为均匀网格时间步的 1/1.6。换句话说，非均匀网格为了稳定/尺度一致，需要更多时间步。
        write(*,*) "timeStep ratio for (uniform) / (non-uniform) : ", real(length_LB / dble(nx))
        write(*,*) "tauf =", real(tauf)
        write(*,*) "    "
    end if

    ! Allocate flow variables
    allocate (u(nx,nyLocal))
    allocate (v(nx,nyLocal))
    allocate (rho(nx,nyLocal))
    allocate (up(nx,nyLocal))
    allocate (vp(nx,nyLocal))
    allocate (temp(nx,nyLocal))
    allocate (utemp(nx,nyLocal))
    allocate (Fx(nx,nyLocal))
    allocate (Fy(nx,nyLocal))

    if(MYID == 0)then            !root保存全场数据
        allocate (u_all(nx,ny))
        allocate (v_all(nx,ny))
        allocate (rho_all(nx,ny))
        allocate (T_all(nx,ny))
    end if

    allocate (f(0:8,nx,nyLocal))
    allocate (f_post(0:8,nx,0:nyLocal+1))
    allocate (g(0:4,nx,nyLocal))
    allocate (g_post(0:4,nx,0:nyLocal+1))

    ! Initialize flow variables
    rho = rho0
    temp = 0.0d0
    utemp=0.0d0
    u = 0.0d0
    v = 0.0d0
    up = 0.0d0
    vp = 0.0d0

    do j=1,nyLocal
        do i=1,nx
            temp(i,j) = dble(i-1)/dble(nx-1)*(Tcold-Thot)+Thot
        enddo
    enddo

    omega_U(0) = 4.0d0/9.0d0
    do alpha=1,4
        omega_U(alpha) = 1.0d0/9.0d0
    enddo
    do alpha=5,8
        omega_U(alpha) = 1.0d0/36.0d0
    enddo

    omega_T(0) = 1.0d0/2.0d0
    do alpha=1,4
        omega_T(alpha) = 1.0d0/8.0d0
    enddo

    do j = 1, nyLocal
        do i = 1, nx
            velocitySquare = u(i,j)*u(i,j)+v(i,j)*v(i,j)
            do alpha = 0,8
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                f(alpha, i, j) = rho(i,j)*omega_U(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*velocitySquare)
            enddo

            do alpha=0,4
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                g(alpha,i,j)=omega_T(alpha)*temp(i,j)*(1.0d0+4.0d0*un(alpha))
            end do
        enddo
    enddo

    ! 为 x 方向插值预先记录每个 i 对应的三个参考点下标。
    ! 内部点使用中心模板 (i-1, i, i+1)；边界点改用单侧模板，避免访问 0 或 nx+1。
    do i = 1, nx
        if(i == 1)then
            ! 左边界没有 i-1，因此使用右侧的合法点 (i+1, i, i+2)。
            inter_x(i,:) = (/i+1, i, i+2/)

        elseif(i == nx)then
            ! 右边界没有 i+1，因此使用左侧的合法点 (i-1, i, i-2)。
            inter_x(i,:) = (/i-1, i, i-2/)

        else
            ! 内部点使用左右各一个邻点的三点插值模板。
            inter_x(i,:) = (/i-1, i, i+1/)
        end if
    enddo

    if(MYID == 0)then
        do j = 1, nyLocal
            if(j == 1)then
                inter_y(j,:) = (/j+1, j, j+2/)

            else
                inter_y(j,:) = (/j-1, j, j+1/)
            end if
        end do

    elseif(MYID == NPROC-1)then
        do j = 1, nyLocal
            if(j == nyLocal)then
                inter_y(j,:) = (/j-1, j, j-2/)

            else
                inter_y(j,:) = (/j-1, j, j+1/)
            end if
        end do

    else
        do j = 1, nyLocal
            inter_y(j,:) = (/j-1, j, j+1/)
        end do
    end if

    return
end subroutine initial

program main
    use commondata
    implicit none
    real(8) :: timestart, timeEnd
    !----------------------------------------------------------
    ! 初始化 MPI 运行环境；后续 MPI 调用都必须在它之后执行。
    call MPI_INIT(IERR)
    ! 获取 MPI_COMM_WORLD 通信域中的总进程数，保存到 NPROC。
    call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,IERR)
    ! 获取当前进程在 MPI_COMM_WORLD 中的编号 MYID，编号范围为 0 到 NPROC-1。
    call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)

    ! 本程序沿 y 方向做一维区域切分：MYID+1 是上方邻居，MYID-1 是下方邻居。
    upid=MYID+1
    downid=MYID-1
    ! 最下方 rank 没有下方邻居，用 MPI_PROC_NULL 表示该方向不通信。
    if(MYID == 0)downid=MPI_PROC_NULL
    ! 最上方 rank 没有上方邻居，用 MPI_PROC_NULL 表示该方向不通信。
    if(MYID == NPROC-1)upid=MPI_PROC_NULL
    !----------------------------------------------------------
    call StartEnd(1, ny)
    ! 从切分结果中取出当前 rank 负责的 y 方向本地网格数。
    nyLocal = count1d(MYID)

    ! 计时前先同步所有 rank，避免有的进程还没准备好就开始计时。
    call MPI_BARRIER(MPI_COMM_WORLD, IERR)
    ! 记录 MPI 墙钟时间，后面用来统计总运行时间和 MLUPS。
    timestart = MPI_WTIME()

    ! 生成非均匀网格坐标，并将完整 x 坐标、局部 y 坐标及 yGrid 坐标 halo 分发给各 rank。
    call mesh()

    call initial()

    do while(((errorU > epsU).or.(errorT > epsT)).AND.(itc < itc_max))

        itc = itc+1

        call collision_U()

        call collision_T()

        call update()

        call interpolate()

        call bounceback_u()

        call bounceback_T()

        call macro_u()

        call macro_t()

        if(MOD(itc,2000).EQ.0) then
            call check()
        endif
    enddo

    call MPI_BARRIER(MPI_COMM_WORLD, IERR)
    timeEnd = MPI_WTIME()

    if(MYID == 0)then
        write(*,*)"Time=", timeEnd-timestart, "s"
        write(*,*) "MLUPS = ", real(dble(nx*ny)/1e6*dble(itc)/(timeEnd-timeStart))
    end if

    call gather_data()

    if(MYID == 0)then
        call output_ASCII()
        write(*,*) "Successfully: DNS completed!"
    end if

    deallocate(u)
    deallocate(v)
    deallocate(rho)
    deallocate(up)
    deallocate(vp)
    deallocate(f)
    deallocate(f_post)
    deallocate(g)
    deallocate(g_post)
    deallocate(temp)
    deallocate(utemp)
    deallocate(yGrid)
    deallocate(inter_y)

    if(MYID == 0)then
        deallocate(u_all)
        deallocate(v_all)
        deallocate(rho_all)
        deallocate(T_all)
    end if

    deallocate(start1d)
    deallocate(end1d)
    deallocate(count1d)
    deallocate(displ1d)
    deallocate(start2d)
    deallocate(end2d)
    deallocate(count2d)
    deallocate(displ2d)

    ! 所有 MPI 通信和并行计算完成后，关闭 MPI 运行环境；此后不再调用 MPI 通信例程。
    call MPI_FINALIZE(IERR)

    stop
end program main

subroutine StartEnd(iS1, iS2)
    use commondata
    implicit none
    integer(kind=4):: leng, iBlock
    integer(kind=4) :: ir
    integer(kind=4) :: iS1, iS2
    integer(kind=4) :: i
    allocate (start1d(0:NPROC-1))
    allocate (end1d(0:NPROC-1))
    allocate (count1d(0:NPROC-1))
    allocate (displ1d(0:NPROC-1))

    allocate (start2d(0:NPROC-1))
    allocate (end2d(0:NPROC-1))
    allocate (count2d(0:NPROC-1))
    allocate (displ2d(0:NPROC-1))

    ! 待切分区间的总长度；这里调用 StartEnd(1, ny)，所以 leng = ny。
    leng = iS2-iS1+1
    ! 每个 rank 至少分到的基础网格数，整数除法会自动舍去余数。
    iBlock = leng/NPROC
    ! 不能被 NPROC 整除的剩余网格数；前 ir 个 rank 每个多分 1 个。
    ir= leng-iBlock*NPROC

    do i = 0, NPROC-1

        if(i.LT.ir) then
            ! 前 ir 个 rank 负责多出来的余数部分，因此本地 y 方向长度为 iBlock+1。
            count1d(i) = iBlock+1
            ! 当前 rank 在全局 y 区间中的起始下标。
            start1d(i) = iS1+i*(iBlock+1)
            ! 当前 rank 在全局 y 区间中的结束下标。
            end1d(i) = start1d(i)+count1d(i)-1
            !-----------------------------------------------------------
            ! 二维场按完整 x 方向乘以本地 y 长度来统计元素数，用于 MPI_GATHERV。
            count2d(i) = (iBlock+1)*nx
            ! 当前 rank 的二维场在全局二维数组一维连续存储视角下的起始位置。
            start2d(i) = iS1+i*(iBlock+1)*nx
            ! 当前 rank 的二维场在全局二维数组一维连续存储视角下的结束位置。
            end2d(i) = start2d(i)+count2d(i)-1

        else
            ! 其余 rank 没有额外余数，只负责基础长度 iBlock。
            count1d(i) = iBlock
            ! 起始下标需要加上 ir，表示前 ir 个 rank 已经各自多拿了 1 个网格。
            start1d(i) = iS1+i*iBlock+ir
            ! 当前 rank 在全局 y 区间中的结束下标。
            end1d(i) = start1d(i)+count1d(i)-1
            !-----------------------------------------------------------
            ! 二维场元素数仍然等于完整 x 方向长度乘以本地 y 长度。
            count2d(i) = iBlock*nx
            ! 二维连续存储视角下的起始位置，同样要跳过前 ir 个 rank 多拿的 nx 个元素。
            start2d(i) = iS1+i*iBlock*nx+ir*nx
            ! 当前 rank 的二维场在全局二维数组一维连续存储视角下的结束位置。
            end2d(i) = start2d(i)+count2d(i)-1
        endif

        ! ScatterV/GatherV 使用的是相对发送缓冲区起点的偏移量，而不是全局下标。
        !从发送/接收大数组的起点开始，跳过多少个元素，才轮到第 i 个 rank 的数据
        displ1d(i) = start1d(i)-iS1
        ! 二维场按一维连续内存收集时使用的偏移量。
        displ2d(i) = start2d(i)-iS1

    enddo
    return
end subroutine StartEnd

subroutine collision_U()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: s(0:8)
    real(kind=8) :: m(0:8)
    real(kind=8) :: m_post(0:8)
    real(kind=8) :: meq(0:8)
    real(kind=8) :: fSource(0:8)

    do j=1,nyLocal
        do i=1,nx

    m(0) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
    m(1) = -4.0d0*f(0,i,j)-f(1,i,j)-f(2,i,j)-f(3,i,j)-f(4,i,j)+2.0d0*(f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j))
    m(2) = 4.0d0*f(0,i,j)-2.0d0*(f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j))+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
    m(3) = f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
    m(4) = -2.0d0*f(1,i,j)+2.0d0*f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
    m(5) = f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)
    m(6) = -2.0d0*f(2,i,j)+2.0d0*f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)
    m(7) = f(1,i,j)-f(2,i,j)+f(3,i,j)-f(4,i,j)
    m(8) = f(5,i,j)-f(6,i,j)+f(7,i,j)-f(8,i,j)

            meq(0) = rho(i,j)
            meq(1) = rho(i,j)*( -2.0d0+3.0d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j)) )
            meq(2) = rho(i,j)*( 1.0d0-3.0d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j)) )
            meq(3) = rho(i,j)*u(i,j)
            meq(4) = -rho(i,j)*u(i,j)
            meq(5) = rho(i,j)*v(i,j)
            meq(6) = -rho(i,j)*v(i,j)
            meq(7) = rho(i,j)*(u(i,j)*u(i,j)-v(i,j)*v(i,j))
            meq(8) = rho(i,j)*(u(i,j)*v(i,j))

            s(0) = 0.0d0
            s(1) = Snu
            s(2) = Snu
            s(3) = 0.0d0
            s(4) = Sq
            s(5) = 0.0d0
            s(6) = Sq
            s(7) = Snu
            s(8) = Snu

            Fx(i,j) = 0.0d0
            Fy(i,j) = gbeta*temp(i,j)*rho(i,j)

            fSource(0) = 0.0d0
            fSource(1) = (6.0d0-3.0d0*s(1))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
            fSource(2) = -(6.0d0-3.0d0*s(2))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
            fSource(3) = Fx(i,j)
            fSource(4) = -(1.0d0-0.5d0*s(4))*Fx(i,j)
            fSource(5) = Fy(i,j)
            fSource(6) = -(1.0d0-0.5d0*s(6))*Fy(i,j)
            fSource(7) = (2.0d0-s(7))*(u(i,j)*Fx(i,j)-v(i,j)*Fy(i,j))
            fSource(8) = (1-0.5d0*s(8))*(v(i,j)*Fx(i,j)+u(i,j)*Fy(i,j))


            do alpha=0,8
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+fSource(alpha)*dt
            enddo

    f_post(0,i,j) = ( m_post(0)-m_post(1)+m_post(2) )/9.0d0
    f_post(1,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+m_post(3)/6.0d0-m_post(4)/6.0d0 &
                    +m_post(7)*0.25d0
    f_post(2,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    +m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)*0.25d0
    f_post(3,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0-m_post(3)/6.0d0+m_post(4)/6.0d0 &
                    +m_post(7)*0.25d0
    f_post(4,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    -m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)*0.25d0
    f_post(5,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)*0.25d0
    f_post(6,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)*0.25d0
    f_post(7,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)*0.25d0
    f_post(8,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)*0.25d0
        enddo
    enddo

    do j = 1, nyLocal
        do i = 1, nx
            f(0,i,j) = f_post(0,i,j)
        enddo
    enddo

    return
end subroutine collision_U

subroutine collision_T()
    use commondata
    implicit none
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: n(0:4), n_post(0:4), neq(0:4)
    real(kind=8) :: Q(0:4)

    do j=1,nyLocal
        do i=1,nx
            n(0) = g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
            n(1) = g(1,i,j)-g(3,i,j)
            n(2) = g(2,i,j)-g(4,i,j)
            n(3) = g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
            n(4) = g(1,i,j)-g(2,i,j)+g(3,i,j)-g(4,i,j)

            Q(0) = 1.0d0
            Q(1) = sig_k
            Q(2) = sig_k
            Q(3) = 1.9d0
            Q(4) = 1.9d0

            neq(0) = temp(i,j)
            neq(1) = temp(i,j)*u(i,j)
            neq(2) = temp(i,j)*v(i,j)
            neq(3) = temp(i,j)*0.5d0
            neq(4) = 0.0d0

            do alpha=0,4
                n_post(alpha)=n(alpha)-Q(alpha)*(n(alpha)-neq(alpha))
            enddo

            g_post(0,i,j) = n_post(0)-n_post(3)
            g_post(1,i,j) = n_post(1)/2.0d0+n_post(3)/4.0d0+n_post(4)/4.0d0
            g_post(2,i,j) = n_post(2)/2.0d0+n_post(3)/4.0d0-n_post(4)/4.0d0
            g_post(3,i,j) = -n_post(1)/2.0d0+n_post(3)/4.0d0+n_post(4)/4.0d0
            g_post(4,i,j) = -n_post(2)/2.0d0+n_post(3)/4.0d0-n_post(4)/4.0d0

        enddo
    enddo

    do j = 1, nyLocal
        do i = 1, nx
            g(0,i,j) = g_post(0,i,j)
        enddo
    enddo

    return
end subroutine collision_T


subroutine update()
    use commondata
    implicit none

    ! 交换 f_post 的上侧 halo：发送本 rank 的下边界 j=1 给 downid，
    ! 同时从 upid 接收其边界层，填入本 rank 的上侧 halo j=nyLocal+1。
    ! f_post 是 D2Q9，整条 y 边界包含 nx 个格点、每点 9 个分布函数，所以数量为 9*nx。
    call MPI_SENDRECV(f_post(0,1,1),9*nx,MPI_REAL8,downid,0,f_post(0,1,nyLocal+1),9*nx,MPI_REAL8,upid,0&
    ,MPI_COMM_WORLD,ISTATUS,IERR)

    ! 交换 f_post 的下侧 halo：发送本 rank 的上边界 j=nyLocal 给 upid，
    ! 同时从 downid 接收其边界层，填入本 rank 的下侧 halo j=0。
    call MPI_SENDRECV(f_post(0,1,nyLocal),9*nx,MPI_REAL8,upid,1,f_post(0,1,0),9*nx,MPI_REAL8,downid,1&
    ,MPI_COMM_WORLD,ISTATUS,IERR)

    ! 对温度分布函数 g_post 做同样的上侧 halo 交换。
    ! g_post 是 D2Q5，整条 y 边界包含 nx 个格点、每点 5 个分布函数，所以数量为 5*nx。
    call MPI_SENDRECV(g_post(0,1,1),5*nx,MPI_REAL8,downid,0,g_post(0,1,nyLocal+1),5*nx,MPI_REAL8,upid,0&
    ,MPI_COMM_WORLD,ISTATUS,IERR)

    ! 对温度分布函数 g_post 做同样的下侧 halo 交换。
    ! 这些 halo 会在后续 interpolate() 中被边界附近的三点插值访问。
    call MPI_SENDRECV(g_post(0,1,nyLocal),5*nx,MPI_REAL8,upid,1,g_post(0,1,0),5*nx,MPI_REAL8,downid,1&
    ,MPI_COMM_WORLD,ISTATUS,IERR)

    return
end subroutine

subroutine interpolate()
    use commondata
    implicit none
    real(kind=8) :: interpolateF, delta_x, delta_y
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: f0, f1, f2, g0, g1, g2

        ! 在非均匀网格上完成分布函数迁移 streaming。
        ! f_post/g_post 是碰撞后、仍存放在非均匀网格节点上的分布函数。
        ! 迁移时并不是另建一个均匀网格，而是仍按标准 LBM 离散速度 e_alpha*dt, 就是按照最小的dt来移动的。
        ! 预测这些 f_post/g_post 从原节点出发后的新位置。
        ! 由于 xGrid/yGrid 间距不均匀，迁移后的新位置通常不会正好落在非均匀网格节点上；
        ! 因此若想得到当前非均匀节点 (i,j) 在下一时刻的 f/g 值，就要用附近 3x3 个
        ! 已按 e_alpha*dt 迁移后的点做二维三点 Lagrange 插值。
        ! 这里选择最小网格间距作为 LB 格子步长，使一个时间步内的迁移距离受限；
        ! 对内部点而言，出发/到达位置会落在所选三点模板覆盖范围内，主要做内插而不是外插。
        ! update() 已经提前补好了 f_post/g_post 的 y 向 halo，边界附近插值可访问 j=0 或 j=nyLocal+1。
        do j = 1, nyLocal
            do i = 1, nx
                do alpha = 1, 8
                    ! 当前离散速度方向在一个时间步内对应的位移。
                    ! 可把周围 f_post(alpha,*,*) 看成先从原非均匀节点按 e_alpha*dt 迁移，
                    ! 即原位置 (xGrid(k),yGrid(l)) 变成 (xGrid(k)+delta_x,yGrid(l)+delta_y)。
                    ! 这些迁移后的点携带下一时刻的数据，但一般不在非均匀网格节点上，
                    ! 所以下面要把它们插值回目标节点 (xGrid(i),yGrid(j))。
                    delta_x=dble(ex(alpha))*dt
                    delta_y=dble(ey(alpha))*dt

            ! 先固定三个 x 参考点，分别沿 y 方向做三点插值，得到 f0/f1/f2。
            ! f0/f1/f2 不是 alpha=0/1/2 的方向分量，而是二维插值的中间值：
            !   f0：固定在第 1 个 x 参考列 inter_x(i,1)，用三个 y 参考点插值得到；
            !   f1：固定在第 2 个 x 参考列 inter_x(i,2)，用三个 y 参考点插值得到；
            !   f2：固定在第 3 个 x 参考列 inter_x(i,3)，用三个 y 参考点插值得到。
            ! 每个 y 插值点的坐标都加 delta_y，表示这些 f_post 已经沿 alpha 方向迁移过；
            ! 插值目标 yGrid(inter_y(j,2)) 通常就是当前目标节点的 yGrid(j)。
            f0 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), f_post(alpha,inter_x(i,1), inter_y(j,1)), f_post(alpha,inter_x(i,1), inter_y(j,2))&
                , f_post(alpha,inter_x(i,1), inter_y(j,3)))

            f1 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), f_post(alpha,inter_x(i,2), inter_y(j,1)), f_post(alpha,inter_x(i,2), inter_y(j,2))&
                , f_post(alpha,inter_x(i,2), inter_y(j,3)))

            f2 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), f_post(alpha,inter_x(i,3), inter_y(j,1)), f_post(alpha,inter_x(i,3), inter_y(j,2))&
                , f_post(alpha,inter_x(i,3), inter_y(j,3)))

            ! 再用 f0/f1/f2 沿 x 方向插值，得到当前格点处迁移后的 f(alpha,i,j)。
            ! 这里的 xGrid(inter_x(i,*))+delta_x 同样表示三个 x 参考列迁移后的坐标；
            ! 目标 xGrid(inter_x(i,2)) 通常就是当前目标节点的 xGrid(i)。
            ! 因此最终得到的是“迁移并投影回非均匀网格后”的下一时刻 f(alpha,i,j)。
            f(alpha, i, j) = interpolateF(xGrid(inter_x(i,1))+delta_x, xGrid(inter_x(i,2))+delta_x, &
                            xGrid(inter_x(i,3))+delta_x, xGrid(inter_x(i,2)), f0, f1, f2)

                end do
            enddo
        enddo

        ! 温度分布函数 g_post 使用同样的二维插值迁移；D2Q5 只需要处理 alpha=1:4。
        do j = 1, nyLocal
            do i = 1, nx
                do alpha = 1, 4
                    ! 温度分布函数同样按 e_alpha*dt 预测迁移后的位置，再插值回非均匀网格格点。
                    delta_x=dble(ex(alpha))*dt
                    delta_y=dble(ey(alpha))*dt

            ! 先沿 y 方向插值得到三个 x 参考位置上的 g0/g1/g2。
            ! g0/g1/g2 与 f0/f1/f2 的含义相同，只是这里处理的是 D2Q5 温度分布函数 g_post。
            g0 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), g_post(alpha,inter_x(i,1), inter_y(j,1)), g_post(alpha,inter_x(i,1), inter_y(j,2))&
                , g_post(alpha,inter_x(i,1), inter_y(j,3)))

            g1 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), g_post(alpha,inter_x(i,2), inter_y(j,1)), g_post(alpha,inter_x(i,2), inter_y(j,2))&
                , g_post(alpha,inter_x(i,2), inter_y(j,3)))

            g2 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), g_post(alpha,inter_x(i,3), inter_y(j,1)), g_post(alpha,inter_x(i,3), inter_y(j,2))&
                , g_post(alpha,inter_x(i,3), inter_y(j,3)))

            ! 再沿 x 方向插值，得到迁移并投影回非均匀网格后的 g(alpha,i,j)。
            g(alpha, i, j) = interpolateF(xGrid(inter_x(i,1))+delta_x, xGrid(inter_x(i,2))+delta_x, &
                            xGrid(inter_x(i,3))+delta_x, xGrid(inter_x(i,2)), g0, g1, g2)
                enddo
            end do
        end do
end subroutine

!!NOTE: consider using compiler-specific directives to suggest inlining if necessary.
pure function interpolateF(x0, x1, x2, x, f0, f1, f2) result(f_interp)
    implicit none
    real(kind=8), intent(in) :: x0, x1, x2, x, f0, f1, f2
    real(kind=8) :: f_interp

    ! 三点 Lagrange 插值公式：用 (x0,f0)、(x1,f1)、(x2,f2) 重构目标位置 x 处的函数值。
    f_interp = ((x - x1) * (x - x2)) / ((x0 - x1) * (x0 - x2)) * f0 + &
               ((x - x0) * (x - x2)) / ((x1 - x0) * (x1 - x2)) * f1 + &
               ((x - x0) * (x - x1)) / ((x2 - x0) * (x2 - x1)) * f2

    return
end function interpolateF

subroutine bounceback_u()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    do j=1,nyLocal
        !Left side
        f(1,1,j) = f_post(3,1,j)
        f(5,1,j) = f_post(7,1,j)
        f(8,1,j) = f_post(6,1,j)

        !Right side
        f(3,nx,j) = f_post(1,nx,j)
        f(6,nx,j) = f_post(8,nx,j)
        f(7,nx,j) = f_post(5,nx,j)
    enddo

    if(MYID==NPROC-1)then
        do i=1,nx
            !Top side
            f(4,i,nyLocal) = f_post(2,i,nyLocal)
            f(7,i,nyLocal) = f_post(5,i,nyLocal)
            f(8,i,nyLocal) = f_post(6,i,nyLocal)
        end do

    elseif(MYID==0)then
        do i=1,nx
            !Bottom side
            f(2,i,1) = f_post(4,i,1)
            f(5,i,1) = f_post(7,i,1)
            f(6,i,1) = f_post(8,i,1)
        end do
    end if
return
end subroutine bounceback_u

subroutine bounceback_T()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    do j=1, nyLocal
        !left
        g(1,1,j) = -g_post(3,1,j) + 2.0d0 * omega_T(1) * Thot
        !!right
        g(3,nx,j) = -g_post(1,nx,j) + 2.0d0 * omega_T(3) * Tcold
    end do

    if(MYID==NPROC-1)then
        do i=1,nx
            !Top side
            g(4,i,nyLocal) = g_post(2,i,nyLocal)
        end do

    elseif(MYID==0)then
        do i=1,nx
            !Bottom side
            g(2,i,1) = g_post(4,i,1)
        end do
    end if

    return
end subroutine

subroutine macro_u()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    do j=1, nyLocal
        do i=1, nx
            rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
            u(i,j) = ( f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)+0.5d0*dt*Fx(i,j))/rho(i,j)
            v(i,j) = ( f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)+0.5d0*dt*Fy(i,j))/rho(i,j)
        enddo
    enddo
    return
end subroutine macro_u

subroutine macro_t()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    do j=1, nyLocal
        do i=1, nx
            temp(i,j) = g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
        end do
    end do

    return
end subroutine macro_t

subroutine check()
    use commondata
    implicit none
    integer :: i, j
    real(kind=8) :: error1, error2, error1_all, error2_all
    real(kind=8) :: error3, error4, error3_all, error4_all

    error1 = 0.0d0
    error2 = 0.0d0
    error1_all = 0.0d0
    error2_all = 0.0d0
    error3 = 0.0d0
    error4 = 0.0d0
    error3_all = 0.0d0
    error4_all = 0.0d0

    do j=1,nyLocal
        do i=1,nx
            error1 = error1+(u(i,j)-up(i,j))*(u(i,j)-up(i,j))+(v(i,j)-vp(i,j))*(v(i,j)-vp(i,j))
            error2 = error2+u(i,j)*u(i,j)+v(i,j)*v(i,j)

            up(i,j) = u(i,j)
            vp(i,j) = v(i,j)
        enddo
    enddo

    ! 将各 rank 上的局部速度误差平方和 error1 做全局求和，并把结果返回给所有 rank。
    ! 参数含义：error1=本 rank 发送的局部值，error1_all=接收全局归约结果，
    !           1=归约 1 个标量，MPI_REAL8=数据类型，MPI_SUM=求和操作，
    !           MPI_COMM_WORLD=参与归约的通信域，IERR=错误码。
    call MPI_ALLREDUCE(error1,error1_all,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
    call MPI_ALLREDUCE(error2,error2_all,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)

    errorU = sqrt(error1_all)/sqrt(error2_all)

    do j=1,nyLocal
        do i=1,nx
            error3 = error3+(temp(i,j)-utemp(i,j))**2
            error4 = error4+temp(i,j)**2

            utemp(i,j) = temp(i,j)
        enddo
    enddo

    call MPI_ALLREDUCE(error3,error3_all,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
    call MPI_ALLREDUCE(error4,error4_all,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)

    errorT = sqrt(error3_all)/sqrt(error4_all)

    if(MYID==0)then
       write(*,*) itc,' ',errorU,' ',errorT
    end if

    return
end subroutine check

subroutine gather_data()
    use commondata
    implicit none

    ! 将每个 rank 的局部二维速度场 u(nx,nyLocal) 收集到 root rank 的全局数组 u_all(nx,ny)。
    ! 参数含义：u=本 rank 发送缓冲区起点，count2d(MYID)=本 rank 发送元素数，
    !           MPI_REAL8=发送数据类型，u_all=root 端接收全局数组，
    !           count2d=每个 rank 接收元素数数组，displ2d=每个 rank 数据在 u_all 中的接收偏移量，
    !           MPI_REAL8=接收数据类型，0=root rank，MPI_COMM_WORLD=通信域，IERR=错误码。
    call MPI_GATHERV(u,count2d(MYID),MPI_REAL8,u_all,count2d,displ2d,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
    call MPI_GATHERV(v,count2d(MYID),MPI_REAL8,v_all,count2d,displ2d,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
    call MPI_GATHERV(rho,count2d(MYID),MPI_REAL8,rho_all,count2d,displ2d,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
    call MPI_GATHERV(temp,count2d(MYID),MPI_REAL8,T_all,count2d,displ2d,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
    return
end subroutine

subroutine output_ASCII()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    character(len=100) :: filename

    write(filename,*) itc
    filename = adjustl(filename)


    open(unit=02,file='MRTcavity-'//trim(filename)//'.dat',status='unknown')

    write(02,*) 'TITLE="thermal convective flows"'
    write(02,*) 'VARIABLES="X" "Y" "U" "V" "T" "rho"'
    write(02,101) nx, ny
    do j=1,ny
        do i=1,nx
            write(02,100) xGrid(i), y(j), u_all(i,j), v_all(i,j), T_all(i,j), rho_all(i,j)
        enddo
    enddo
100 format(1x,2(e11.4,' '),10(e13.6,' '))
101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
    close(02)

    return
end subroutine output_ASCII
