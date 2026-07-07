!!!    This program sloves Buoyancy Driven Cavity Flow problem using Lattice Boltzmann Method
!!!    Lattice Boltzmann Equation with MRT-LBE model

!!!~~velocity B.C.~~
#define HorizontalWallsNoslip
#define VerticalWallsNoslip
!!!~~velocity B.C.~~

!!!~~temperature B.C. (for Side Heated Cell)~~
#define SideHeatedCell
#define HorizontalWallsAdiabatic
#define VerticalWallsConstT
!!!~~temperature B.C.~~

!!!~~output~~
! #define __OUT_PUT__

module commondata
    implicit none

    integer(kind=4), parameter :: total_nx=8193, total_ny=total_nx     !----Section 1----!
    integer :: nx, ny, i_start_global, j_start_global
    real(kind=8), parameter :: lengthUnit=dble(total_ny)     !----Section 1----!
    
    real(kind=8), parameter :: Rayleigh=1e8        !----Section 2----!
    real(kind=8), parameter :: Prandtl=0.71d0       !----Section 2----!
    real(kind=8), parameter :: Mach=0.1d0           !----Section 2----!
    
    real(kind=8), parameter :: outputFrequency=1.0d0 !~unit free fall time                            !----Section 3----!
    
    integer(kind=4), parameter :: dimensionlessTimeMax=int(3000/outputFrequency)  !----Section 3----!

    !----Section 3----!
    real(kind=8), parameter :: epsU=1e-6                     !----Section 3----!
    real(kind=8), parameter :: epsT=1e-6                     !----Section 3----!
    
    integer(kind=4), parameter :: outputBinFile=0, outputPltFile=1                 !----Section *----!
            
    real(kind=8), parameter :: Pi=4.0d0*datan(1.0d0)
    !-----------------------------------------------------------------------------------------------        
    !----Section 1----!
    integer(kind=4), parameter :: nxHalf=(total_nx-1)/2+1, nyHalf=(total_ny-1)/2+1
    
    !-----------------------------------------------------------------------------------------------
    !----Section 2----!
    real(kind=8), parameter :: rho0=1.0d0 !~!~m.u.
    real(kind=8), parameter :: Thot=1.0d0, Tcold=0.0d0, Tref=0.0d0
    real(kind=8), parameter :: tauf=0.5d0+Mach*lengthUnit*DSQRT(3.0d0*Prandtl/Rayleigh)
    real(kind=8), parameter :: viscosity=(tauf-0.5d0)/3.0d0
    real(kind=8), parameter :: diffusivity=viscosity/Prandtl
    
    real(kind=8), parameter :: paraA=20.0d0*dsqrt(3.0d0)*diffusivity-4.0d0
    real(kind=8), parameter :: gBeta1=Rayleigh*viscosity*diffusivity/lengthUnit
    real(kind=8), parameter :: gBeta=gBeta1/lengthUnit/lengthUnit
    
    real(kind=8), parameter :: timeUnit=dsqrt(lengthUnit/gBeta)  !!dble(ny*ny)/diffusivity


    real(kind=8), parameter :: Snu=1.0d0/tauf, Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)
    real(kind=8), parameter :: Qd=3.0d0-dsqrt(3.0d0), Qnu=4.0d0*dsqrt(3.0d0)-6.0d0
    !-----------------------------------------------------------------------------------------------
    !----Section 3----!
    real(kind=8) :: errorU, errorT

    !-----------------------------------------------------------------------------------------------
    !----Section 5----!
    real(kind=8) :: xp(0:total_nx+1), yp(0:total_ny+1)
    real(kind=8), allocatable :: u(:,:), v(:,:), T(:,:), rho(:,:)

    real(kind=8), allocatable :: up(:,:), vp(:,:), Tp(:,:)

    real(kind=8), allocatable :: f(:,:,:), f_post(:,:,:)
    real(kind=8), allocatable :: g(:,:,:), g_post(:,:,:)
    real(kind=8), allocatable :: Fx(:,:), Fy(:,:)
    
    integer :: index
    integer(kind=4), parameter :: ex(0:8) = (/ 0, 1, 0, -1,  0, 1, -1, -1,  1 /) 
    integer(kind=4), parameter :: ey(0:8) = (/ 0, 0, 1,  0, -1, 1,  1, -1, -1 /)
    integer(kind=4), parameter :: r(1:8) = (/ 3, 4, 1, 2, 7, 8, 5, 6 /)
    real(kind=8), parameter :: omega(0:8) = (/ 4.0d0/9.0d0, &
                                        (1.0d0/9.0d0, index = 1, 4), &
                                        (1.0d0/36.0d0, index = 5, 8) /) 
    real(kind=8), parameter :: omegaT(0:4) = (/ (1.0d0-paraA)/5.0d0, &
                                            ((paraA+4.0d0)/20.0d0, index = 1,4) /)
    !-----------------------------------------------------------------------------------------------
    !----Section 6----!
    integer(kind=4) :: itc
    integer(kind=4), parameter :: itc_max=dimensionlessTimeMax*int(outputFrequency*timeUnit)


    ! 在GPU端为这些模块变量建立设备副本，后续OpenACC核函数可直接访问。
    ! declare create只负责让OpenACC运行时维护对应的设备端存储/映射关系。
    !$acc declare create(u, v, T, rho, up, vp, Tp, f, f_post, g, g_post, Fx, Fy) &
    !$acc create(nx, ny)
end module commondata
        
module mpi_data
    implicit none
    integer :: rc, rank, num_process
    ! MPI二维笛卡尔拓扑用的一维数组：
    ! dims表示x/y方向的进程数，coords表示当前rank在x/y方向的坐标。
    integer :: dims(0:1) = (/0, 0/), coords(0:1)
    logical :: periods(0:1)
    data periods/2*.false./
    integer :: comm2d, rank2d
    integer :: nbr_left, nbr_right, nbr_top, nbr_bottom
    real(8), allocatable :: send_pos(:), recv_pos(:), send_neg(:), recv_neg(:)
    real(8), allocatable :: g_send_pos_x(:), g_recv_pos_x(:), g_send_neg_x(:), g_recv_neg_x(:)
    integer, parameter :: f_tag_x_pos(0:2) = (/ 1, 5, 8 /)
    integer, parameter :: f_tag_x_neg(0:2) = (/ 3, 6, 7 /)
    integer, parameter :: f_tag_y_pos(0:2) = (/ 2, 5, 6 /)
    integer, parameter :: f_tag_y_neg(0:2) = (/ 4, 7, 8 /)

    ! 在GPU端为MPI通信相关变量建立设备副本：
    ! send/recv数组用于边界数据交换，dims/coords用于GPU端识别当前MPI子区域位置。
    ! 这里只建立OpenACC设备端映射；通信缓冲区本身仍在allocate_all中分配。
    !$acc declare create(send_pos(:), recv_pos(:), send_neg(:), recv_neg(:)) &
    !$acc create(g_send_pos_x(:), g_recv_pos_x(:), g_send_neg_x(:), g_recv_neg_x(:)) &
    !$acc create(dims, coords)
end module mpi_data
    

program main
    use mpi
    use mpi_data    
    use commondata
    implicit none
    real(kind=8) :: timeStart, timeEnd
    real(8) :: start_time, end_time
    integer :: i, j, idx, tmp

    call mpi_starts()

    ! mpi_starts在CPU端确定了当前rank的局部网格尺寸nx/ny和拓扑信息dims/coords。
    ! 由于这些变量前面用declare create建立了GPU端副本，这里需要显式同步到设备端。
    !$acc update device(nx, ny, dims, coords)
    
    call allocate_all()

    call initial()

#ifdef __OUT_PUT__
    call output()
#endif
                
    call CPU_TIME(timeStart)

    call MPI_Barrier(MPI_COMM_WORLD, rc)
    start_time = MPI_Wtime()

    do while( ((errorU.GT.epsU).OR.(errorT.GT.epsT)).AND.(itc.LE.itc_max))

        itc = itc+1
        
        call flow_update()

        if(MOD(itc,2000).EQ.0) call check()

        if(MOD(itc,12000).EQ.0) exit

    enddo

    call MPI_Barrier(MPI_COMM_WORLD, rc)
    end_time = MPI_Wtime()
        
    call CPU_TIME(timeEnd)
    if(rank == 0) then
        write(*,*) "Time (CPU) = ", real(timeEnd-timeStart), "s"
        write(*,*) "Time (MPI) = ", real(end_time - start_time), "s"
    endif

#ifdef __OUT_PUT__
    call output()
#endif
    
   call free_all()

   if(rank == 0) then
    write(*,*) "Successfully: DNS completed!"
    endif

    call MPI_Finalize(rc)

end program main


subroutine mpi_starts()
    use openacc 
    use mpi
    use mpi_data
    use commondata
    implicit none
    integer :: num_gpus, gpu_id
    integer :: local_comm, local_rank
    integer :: name_len, tmp
    ! MPI_MAX_PROCESSOR_NAME是MPI规定的处理器/节点名最大长度。
    ! processor_name用于保存当前rank所在节点名，name_len保存实际名称长度。
    character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name

    call MPI_Init(rc)

    call MPI_Comm_size(MPI_COMM_WORLD, num_process, rc)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, rc)
    ! 获取当前MPI进程所在节点名：processor_name为输出名称，name_len为名称长度，rc为返回状态码。
    call MPI_Get_processor_name(processor_name, name_len, rc)

    !!! decomposition the domain 
    ! call MPI_Dims_create(num_process, 2, dims, rc)
    dims(0) = 1
    ! 根据总MPI进程数、全局网格尺寸和周期性设置，选择二维进程拓扑dims(0:1)。
    ! 这里dims(0)=1已提前固定，所以该函数主要确定y方向需要分成多少个MPI子块。
    call MPI_Dims_create_2d(num_process, dims, total_nx, total_ny, periods, rc)
    
    ! 创建二维MPI笛卡尔通信器：
    ! MPI_COMM_WORLD为原始全局通信器；2表示二维拓扑；
    ! dims给出x/y方向进程数；periods给出x/y方向是否周期；
    ! .true.允许MPI重新排列rank以优化拓扑；comm2d为输出的新通信器；rc为返回状态码。
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, .true., comm2d, rc)
    if(rank == 0) then
        write(*,*) "dimens is x*y = ", dims(0), "x", dims(1)
    endif

    ! get my new rank in decomposition，得到新的rank
    call MPI_Comm_rank(comm2d, rank2d, rc)
    ! write(*,*) "process ", rank2d, " of total ", num_process, "is alive."

    ! determain sub-domain size
    ! 从二维笛卡尔通信器中取回拓扑信息：
    ! dims/periods分别返回各方向进程数和周期性，coords返回当前rank的二维坐标。
    ! 这里主要使用coords(0:1)来计算当前rank负责的全局子区域起点和局部尺寸。
    call MPI_Cart_get(comm2d, 2, dims, periods, coords, rc)
    
    ! 分别沿x/y方向做一维网格切分（余数前置）：
    ! 输入全局网格数total_nx/total_ny、当前rank拓扑坐标coords和方向进程数dims；
    ! 输出当前rank的局部网格尺寸nx/ny，以及在全局网格中的起点i_start_global/j_start_global。
    call decompose_1d(total_nx, nx, coords(0), dims(0), i_start_global)
    call decompose_1d(total_ny, ny, coords(1), dims(1), j_start_global)
    if(rank == 0) then
        write(*,*) "Total nx * ny =", total_nx, " * ", total_ny
        write(*,*) "local nx * ny =", nx, " * ", ny
    endif
    ! write(*,*) "coords = ", coords(1), coords(2)
    ! write(*,*) "nx*ny = ", nx, ny

    ! 查询当前rank在二维笛卡尔拓扑中的相邻rank，用于后续MPI边界/halo交换。
    ! 参数含义：comm2d为二维通信器；第二个参数0/1分别表示x/y方向；
    ! 第三个参数1表示沿该方向偏移一个子块；后两个输出分别是负方向和正方向邻居。
    ! x方向得到左/右邻居，y方向得到下/上邻居；若边界无邻居且非周期，邻居为MPI_PROC_NULL。
    call MPI_Cart_shift(comm2d, 0, 1, nbr_left, nbr_right, rc)
    call MPI_Cart_shift(comm2d, 1, 1, nbr_bottom, nbr_top, rc)


    local_rank = -1
    ! MPI_Comm_split_type参数含义：
    ! comm2d为输入通信器；MPI_COMM_TYPE_SHARED表示按共享内存节点分组；
    ! rank2d作为key，用于决定新通信器local_comm中的rank排序；
    ! MPI_INFO_NULL表示不提供额外分组提示；local_comm为输出的节点内通信器；rc为返回状态码。
    call MPI_Comm_split_type(comm2d, MPI_COMM_TYPE_SHARED, rank2d, MPI_INFO_NULL, local_comm, rc)
    ! 获取当前rank在节点内的编号local_rank
    call MPI_Comm_rank(local_comm, local_rank, rc)
    ! 释放临时通信器local_comm
    call MPI_Comm_free(local_comm, rc)
    
    ! 查询当前节点上OpenACC可见的NVIDIA GPU数量。
    num_gpus = 0
    num_gpus = acc_get_num_devices(acc_device_nvidia)
    ! 每个节点只让local_rank=0的进程打印一次GPU数量，避免同节点重复输出。
    if (local_rank .eq. 0) then
        write(*,*) "I from rank", rank2d, "we have", num_gpus, "gpus"
    endif
    ! 如果没有可用GPU，则终止MPI作业。
    if (num_gpus .le. 0) then
        if (rank2d .eq. 0) then
            write(*,*) 'No NVIDIA GPUs available'
            call MPI_Abort(MPI_COMM_WORLD, 1, rc)
        endif
    else
            ! 用节点内local_rank对GPU数量取模，为当前MPI rank选择GPU。
            ! 例如4张GPU时，local_rank=0/1/2/3分别绑定GPU 0/1/2/3；
            ! 若本节点rank数多于GPU数，则多个rank会循环共享GPU。
            gpu_id = mod(local_rank, num_gpus)
            write(*,*) "i'm local rank", local_rank, "rank", rank2d, "using gpu ", gpu_id
            ! 设置当前MPI rank后续OpenACC计算使用的NVIDIA GPU编号。
            call acc_set_device_num(gpu_id, acc_device_nvidia)
    endif

    !! Call acc_init after acc_set_device_num to avoid multiple contexts on device 0 in multi GPU systems
    call acc_init(acc_device_nvidia)
contains
    subroutine decompose_1d(total_n, local_n, rank, num_process, i_start_global)
        implicit none
        integer, intent(in) :: total_n, rank, num_process
        integer, intent(out) :: local_n, i_start_global

        local_n = total_n / num_process

        if (rank < MOD(total_n, num_process)) then
            local_n = local_n + 1
        endif

        if (local_n > total_n / num_process) then ! --- 5 5 '5' 4 4 4
            i_start_global = local_n * rank
        else                    ! --- 5 5 5 4 '4' 4
            i_start_global = local_n * rank + mod(total_n, num_process)
        endif

    end subroutine decompose_1d

    subroutine MPI_Dims_create_2d(num_process, dims, total_nx, total_ny, periods, rc)
        integer, intent(in) :: num_process, total_nx, total_ny, rc
        logical, intent(in) :: periods(0:1)
        integer, intent(inout) :: dims(0:1)
        integer :: i, j, k, is(0:1), ie(0:1)
        real :: message, diff, lx, ly, lz
        integer :: restrict(0:1)

        ! determine the dimensions of cartesian topologies to minimize the message exchange
        ! under user provided restrictions if dims(0:2) != 0

        ! 保存调用者预先给出的dims限制：
        ! restrict(k)=0表示第k个方向可自由搜索，非0表示该方向的进程数被固定。
        restrict = dims

        ! diff保存当前找到的最小通信量估计值。
        ! 先给一个足够大的初值，保证第一组合法划分可以更新它。
        diff = dble(total_nx) * dble(total_ny) * 10.0d0

        ! 默认两个方向都从1搜索到num_process。
        is = 1
        ie = num_process
        ! 如果某个方向已被restrict固定，就把该方向的搜索起点和终点都设为固定值。
        do i = 0, 1
            if (restrict(i) .NE. 0) then
                ! 例如restrict(0)=1时，x方向只允许i=1，不再枚举其他x方向分块数。
                is(i) = restrict(i)
                ie(i) = restrict(i)
            endif
        enddo
        

        do i = is(0), ie(0)
            do j = is(1), ie(1)
                if (i * j == num_process) then
                    message = 0.0d0
                    ! 当前候选划分为 i*j 个MPI子块：
                    ! lx/ly是单个rank负责的局部子区域在x/y方向上的近似长度。
                    lx = dble(total_nx) / dble(i)
                    ly = dble(total_ny) / dble(j)

                    ! 估计“单个rank最多”需要交换的边界长度，而不是全局所有切分线总长度。
                    if (i > 1) then
                        ! x方向被切分时，每个rank至少可能和一侧邻居交换一条竖直边界，长度约为ly。
                        message = message + ly
                        if (i > 2 .OR. periods(0)) then
                            ! i>2时中间rank有左右两个邻居；若x方向周期，即使i=2也有两侧通信。
                            ! 因此这里最多只再加一条ly，不是按全局(i-1)条切分线累加。
                            message = message + ly 
                        endif
                    endif
                    if (j > 1) then
                        ! y方向被切分时，每个rank至少可能和一侧邻居交换一条水平边界，长度约为lx。
                        message = message + lx
                        if (j > 2 .OR. periods(1)) then
                            ! j>2时中间rank有上下两个邻居；若y方向周期，即使j=2也有两侧通信。
                            message = message + lx 
                        endif
                    endif


                    if (message < diff) then
                        ! 保存当前找到的最小“单rank最大通信边界长度”对应的二维进程划分。
                        diff = message
                        dims(0) = i
                        dims(1) = j
                    endif
                endif
            enddo
        enddo

    end subroutine MPI_Dims_create_2d

end subroutine mpi_starts



subroutine allocate_all()
    use commondata
    use mpi_data
    implicit none
    integer :: max_length

    allocate (u(nx,ny))
    allocate (v(nx,ny))
    allocate (T(nx,ny))
    allocate (rho(nx,ny))
    

    allocate (up(nx,ny))
    allocate (vp(nx,ny))
    allocate (Tp(nx,ny))

    
    allocate (f(nx,ny,0:8))
    allocate (f_post(0:nx+1,0:ny+1,0:8))
    allocate (g(nx,ny,0:4))
    allocate (g_post(0:nx+1,0:ny+1,0:4))
    
    allocate (Fx(nx,ny))
    allocate (Fy(nx,ny))

    ! allocate buffer layer
    ! max_length取x/y两个方向通信边界长度的较大值，即max(nx,ny)+2。
    ! +2包含两侧halo/ghost索引0和nx+1或ny+1，保证x/y方向MPI通信缓冲区都够用。
    max_length = nx+2
    if (ny > nx) then
        max_length = ny+2
    endif

    allocate(send_pos(1 : 3*max_length))
    allocate(recv_pos(1 : 3*max_length))
    allocate(send_neg(1 : 3*max_length))
    allocate(recv_neg(1 : 3*max_length))

    allocate(g_send_pos_x(max_length))
    allocate(g_recv_pos_x(max_length))
    allocate(g_send_neg_x(max_length))
    allocate(g_recv_neg_x(max_length))

end subroutine allocate_all


subroutine free_all()
    use commondata
    use mpi_data
    implicit none

    deallocate(f)
    deallocate(g)
    deallocate(f_post)
    deallocate(g_post)

    deallocate(u)
    deallocate(v)
    deallocate(T)
    deallocate(rho)

    deallocate(up)
    deallocate(vp)
    deallocate(Tp)

    deallocate(Fx)
    deallocate(Fy)

    ! free buffer layer
    deallocate(send_pos)
    deallocate(recv_pos)
    deallocate(send_neg)
    deallocate(recv_neg)

    deallocate(g_send_pos_x)
    deallocate(g_recv_pos_x)
    deallocate(g_send_neg_x)
    deallocate(g_recv_neg_x)

end subroutine free_all
 

subroutine initial()
    use mpi_data
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: un(0:8)
    real(kind=8) :: us2

    
    itc = 0
    errorU = 100.0d0
    errorT = 100.0d0
    !格子长度
    xp(0) = 0.0d0
    xp(total_nx+1) = dble(total_nx)
    do i=1,total_nx
        xp(i) = dble(i)-0.5d0
    enddo
    yp(0) = 0.0d0
    yp(total_ny+1) = dble(total_ny)
    do j=1,total_ny
        yp(j) = dble(j)-0.5d0
    enddo
    

    ! OpenACC并行循环指令：
    ! parallel表示在GPU上开启并行计算区域；loop表示紧跟的循环要并行执行；
    ! collapse(2)表示把后面j/i两层嵌套循环合并成一个二维格点任务集合并行。
    !$acc parallel loop collapse(2)
    do j = 1, ny
        do i = 1, nx
            rho(i, j) = rho0
            u(i, j) = 0.0d0
            v(i, j) = 0.0d0
            T(i, j) = 0.0d0

            up(i, j) = 0.0d0
            vp(i, j) = 0.0d0
            Tp(i, j) = 0.0d0
        enddo
    enddo

    ! 开启OpenACC kernels区域，由编译器分析下面的循环并生成GPU kernel执行。
    ! 这里用于在GPU端初始化温度场T；具体循环并行性由后面的!$acc loop说明。
    !$acc kernels
#ifdef VerticalWallsConstT
    ! 这里不用写parallel，因为外层!$acc kernels已经开启了OpenACC计算区域。
    ! loop只负责修饰下面的循环；collapse(2)把j/i两层循环合并后交给该kernels区域并行执行。
    !$acc loop collapse(2)
    do j=1,ny
        do i=1,nx
            T(i, j) = dble(i_start_global + i - 1) / dble(total_nx - 1) * (Tcold - Thot) + Thot
        enddo
    enddo
#endif
    !$acc end kernels

    ! f = 0.0d0
    ! g = 0.0d0
    
    ! 并行初始化f/g分布函数；collapse(2)并行展开j/i格点循环。
    ! private表示每个并行格点任务都有自己的alpha、us2和un临时变量副本，
    ! 避免不同GPU线程同时写同一临时变量造成数据竞争。
    !$acc parallel loop collapse(2) private(alpha, us2, un)
    do j=1,ny
        do i=1,nx
            us2 = u(i,j)*u(i,j)+v(i,j)*v(i,j)
            do alpha=0,8
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                f(i,j,alpha) = rho(i,j)*omega(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)
            enddo
            do alpha=0,4
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                g(i,j,alpha) = T(i,j)*omegaT(alpha)*(1.0d0+10.0d0/(4.0d0+paraA)*un(alpha))
            enddo
        enddo
    enddo
     
    ! f_post = 0.0d0
    ! g_post = 0.0d0
    
    return
end subroutine initial

    
subroutine collision(i_start, i_end, j_start, j_end)
    ! routine表示为该子程序生成OpenACC设备端例程，使其可在GPU端执行/调用。
    ! gang表示该例程内部包含gang级并行工作，下面的i/j循环会作为较粗粒度GPU并行循环。
    ! nohost表示只生成设备端版本，不额外生成OpenACC主机端例程版本。
    !$acc routine gang nohost
    use commondata
    implicit none
    integer, intent(in) :: i_start, i_end, j_start, j_end
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: m(0:8), m_post(0:8), meq(0:8)
    real(kind=8) :: s(0:8)
    real(kind=8) :: fSource(0:8)

    ! loop independent表示程序员声明各个(i,j)格点迭代彼此无数据依赖，可安全并行。
    ! collapse(2)把j/i两层循环合并成一个二维格点任务集合。
    ! private为每个并行格点任务提供独立的m/meq/m_post/s/fSource临时数组副本。
    !$acc loop independent collapse(2) &
    !$acc private(m, meq, m_post, s, fSource)
    do j = j_start, j_end
        do i = i_start, i_end
            m(0) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
            m(1) = -4.0d0*f(i,j,0)-f(i,j,1)-f(i,j,2)-f(i,j,3)-f(i,j,4)+2.0d0*(f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8))
            m(2) = 4.0d0*f(i,j,0)-2.0d0*(f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4))+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
            m(3) = f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)
            m(4) = -2.0d0*f(i,j,1)+2.0d0*f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)
            m(5) = f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)
            m(6) = -2.0d0*f(i,j,2)+2.0d0*f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)
            m(7) = f(i,j,1)-f(i,j,2)+f(i,j,3)-f(i,j,4)
            m(8) = f(i,j,5)-f(i,j,6)+f(i,j,7)-f(i,j,8)

            meq(0) = rho(i,j)
            meq(1) = rho(i,j)*( -2.0d0+3.0d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j)) )
            meq(2) = rho(i,j)*( 1.0d0-3.0d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j)) )
            meq(3) = rho(i,j)*u(i,j)
            meq(4) = -rho(i,j)*u(i,j)
            meq(5) = rho(i,j)*v(i,j)
            meq(6) = -rho(i,j)*v(i,j)
            meq(7) = rho(i,j)*( u(i,j)*u(i,j)-v(i,j)*v(i,j) )
            meq(8) = rho(i,j)*( u(i,j)*v(i,j) ) 

            s(0) = 0.0d0      !!s_{\rho}
            s(1) = Snu !!s_{e}
            s(2) = Snu !!s_{\epsilon}
            s(3) = 0.0d0      !!s_{j} 
            s(4) = Sq !!s_{q}
            s(5) = 0.0d0      !!s_{j}
            s(6) = Sq       !!s_{q}
            s(7) = Snu !!s_{\nu}
            s(8) = Snu       !!s_{\nu}

            Fx(i,j) = 0.0d0
            Fy(i,j) = rho(i,j)*gBeta*(T(i,j)-Tref)

            fSource(0) = 0.0d0
            fSource(1) = (6.0d0-3.0d0*s(1))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
            fSource(2) = -(6.0d0-3.0d0*s(2))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
            fSource(3) = (1.0d0-0.5d0*s(3))*Fx(i,j)
            fSource(4) = -(1.0d0-0.5d0*s(4))*Fx(i,j)
            fSource(5) = (1.0d0-0.5d0*s(5))*Fy(i,j)
            fSource(6) = -(1.0d0-0.5d0*s(6))*Fy(i,j)
            fSource(7) = (2.0d0-s(7))*(u(i,j)*Fx(i,j)-v(i,j)*Fy(i,j))
            fSource(8) = (1.0d0-0.5d0*s(8))*(u(i,j)*Fy(i,j)+v(i,j)*Fx(i,j))

            ! 当前外层(i,j)格点循环已经并行；这里alpha=0..8只是单个格点内部的小循环。
            ! loop seq表示该内层循环在当前并行任务内顺序执行，不再拆成额外GPU并行循环。
            !$acc loop seq
            do alpha=0,8
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+fSource(alpha)
            enddo

            f_post(i,j,0) = ( m_post(0)-m_post(1)+m_post(2) )/9.0d0
            f_post(i,j,1) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+m_post(3)/6.0d0-m_post(4)/6.0d0 &
                            +m_post(7)/4.0d0
            f_post(i,j,2) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                            +m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)/4.0d0
            f_post(i,j,3) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0-m_post(3)/6.0d0+m_post(4)/6.0d0 &
                            +m_post(7)/4.0d0
            f_post(i,j,4) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                            -m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)/4.0d0
            f_post(i,j,5) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                            +m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)/4.0d0
            f_post(i,j,6) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                            +m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)/4.0d0
            f_post(i,j,7) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                            -m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)/4.0d0
            f_post(i,j,8) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                            -m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)/4.0d0

        enddo
    enddo
    
    return
end subroutine collision


subroutine streaming()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: ip, jp
    integer(kind=4) :: alpha
    
    !$acc parallel loop independent collapse(2)
    do j=1,ny
        do i=1,nx
            !$acc loop seq
            do alpha=0,8
                ip = i-ex(alpha)
                jp = j-ey(alpha)
            
                f(i,j,alpha) = f_post(ip,jp,alpha)
            enddo
        enddo
    enddo
    
    return
end subroutine streaming


subroutine bounceback()
    use mpi_data
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !$acc kernels
#ifdef VerticalWallsNoslip
    !Left side (i=1)
    if (coords(0) == 0) then
        !$acc loop independent
        do j=1,ny 
            f(1,j,1) = f_post(1,j,3)
            f(1,j,5) = f_post(1,j,7)
            f(1,j,8) = f_post(1,j,6)
        enddo
    endif

    !Right side (i=nx)
    if (coords(0) == dims(0)-1) then    
        !$acc loop independent
        do j=1,ny
            f(nx,j,3) = f_post(nx,j,1)
            f(nx,j,6) = f_post(nx,j,8)
            f(nx,j,7) = f_post(nx,j,5)
        enddo
    endif
#endif

#ifdef HorizontalWallsNoslip
    !Bottom side (j=1)
    if (coords(1) == 0) then
        !$acc loop independent
        do i=1,nx 
            f(i,1,2) = f_post(i,1,4)
            f(i,1,5) = f_post(i,1,7)
            f(i,1,6) = f_post(i,1,8)
        enddo
    endif

    !Top side (j=ny)
    if (coords(1) == dims(1) - 1) then
        !$acc loop independent
        do i=1,nx
            f(i,ny,4) = f_post(i,ny,2)
            f(i,ny,7) = f_post(i,ny,5)
            f(i,ny,8) = f_post(i,ny,6)
        enddo
    endif
#endif
    !$acc end kernels

end subroutine bounceback


subroutine macro()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !$acc parallel loop independent collapse(2)
    do j=1,ny
        do i=1,nx
            rho(i,j) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
            u(i,j) = ( f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)+0.5d0*Fx(i,j) )/rho(i,j)
            v(i,j) = ( f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)+0.5d0*Fy(i,j) )/rho(i,j)
        enddo
    enddo

    return
end subroutine macro
    

subroutine collisionT(i_start, i_end, j_start, j_end)
    !$acc routine gang nohost
    use commondata
    implicit none
    integer, intent(in) :: i_start, i_end, j_start, j_end
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    !------------------------
    real(kind=8) :: n(0:4), n_post(0:4), neq(0:4)
    real(kind=8) :: q(0:4)

    !$acc loop independent collapse(2) &
    !$acc private(n, n_post, neq, q)
    do j = j_start, j_end
        do i = i_start, i_end 
            n(0) = g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
            n(1) = g(i,j,1)-g(i,j,3)
            n(2) = g(i,j,2)-g(i,j,4)
            n(3) = -4.0d0*g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
            n(4) = g(i,j,1)-g(i,j,2)+g(i,j,3)-g(i,j,4)
            
            neq(0) = T(i,j)
            neq(1) = T(i,j)*u(i,j)
            neq(2) = T(i,j)*v(i,j)
            neq(3) = T(i,j)*paraA
            neq(4) = 0.0d0
        
            q(0) = 0.0d0
            q(1) = Qd
            q(2) = Qd
            q(3) = Qnu
            q(4) = Qnu
            
            !$acc loop seq
            do alpha=0,4
                n_post(alpha) = n(alpha)-q(alpha)*(n(alpha)-neq(alpha))
            enddo
            
            g_post(i,j,0) = 0.2d0*n_post(0)-0.2d0*n_post(3)
            g_post(i,j,1) = 0.2d0*n_post(0)+0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
            g_post(i,j,2) = 0.2d0*n_post(0)+0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4)
            g_post(i,j,3) = 0.2d0*n_post(0)-0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
            g_post(i,j,4) = 0.2d0*n_post(0)-0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4) 
        enddo
    enddo
    

    return
end subroutine collisionT


subroutine streamingT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: ip, jp
    integer(kind=4) :: alpha
    
    !$acc parallel loop independent collapse(2) &
    !$acc private(ip, jp)
    do j = 1, ny
        do i = 1, nx
            !$acc loop seq
            do alpha = 0, 4
                ip = i-ex(alpha)
                jp = j-ey(alpha)
            
                g(i,j,alpha) = g_post(ip,jp,alpha)
            enddo
        enddo
    enddo
    
    return
end subroutine streamingT


subroutine bouncebackT()
    use mpi_data
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !$acc kernels
#ifdef HorizontalWallsAdiabatic
    !Bottom side
    if (coords(1) == 0) then
        !$acc loop independent
        do i = 1, nx 
            g(i, 1, 2) = g_post(i, 1, 4)
        enddo
    endif

    ! Top side
    if (coords(1) == dims(1) - 1) then
        !$acc loop independent
        do i = 1, nx 
            g(i, ny, 4) = g_post(i, ny, 2)
        enddo
    endif
#endif

#ifdef VerticalWallsConstT
    !Left side
    if (coords(0) == 0) then
        !$acc loop independent
        do j = 1, ny 
            g(1,j,1) = -g_post(1,j,3)+(4.0d0+paraA)/10.0d0*Thot
        enddo
    endif

    !Right side
    if (coords(0) == dims(0) - 1) then
        !$acc loop independent
        do j = 1, ny 
            g(nx,j,3) = -g_post(nx,j,1)+(4.0d0+paraA)/10.0d0*Tcold
        enddo
    endif
#endif
    !$acc end kernels

end subroutine bouncebackT


subroutine macroT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    
    !$acc parallel loop independent collapse(2)
    do j = 1, ny
        do i = 1, nx
            T(i,j) = g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
        enddo
    enddo

    return
end subroutine macroT
    
subroutine flow_update()
    use commondata
    use mpi_data
    use mpi
    implicit none

    !$acc kernels
    call collision(1, nx, 1, ny)
    !$acc end kernels

    !$acc kernels
    call collisionT(1, nx, 1, ny)
    !$acc end kernels

    call message_passing_f()

    call message_passing_g()

    call streaming()

    call bounceback()

    call streamingT()

    call bouncebackT()

    call macro()
    
    call macroT()

end subroutine


subroutine check()
    use mpi
    use mpi_data
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: error1, error2, error5, error6
    real(8) :: total_error1, total_error2, total_error5, total_error6

    error1 = 0.0d0
    error2 = 0.0d0

    error5 = 0.0d0
    error6 = 0.0d0
    
    !$acc kernels 
    !$acc loop independent reduction(+:error1, error2, error5, error6)
    do j=1,ny
        do i=1,nx
            error1 = error1+(u(i,j)-up(i,j))*(u(i,j)-up(i,j))+(v(i,j)-vp(i,j))*(v(i,j)-vp(i,j))
            error2 = error2+u(i,j)*u(i,j)+v(i,j)*v(i,j)
            
            error5 = error5+dABS( T(i,j)-Tp(i,j) )
            error6 = error6+dABS( T(i,j) )
            
            up(i,j) = u(i,j)
            vp(i,j) = v(i,j)
            Tp(i,j) = T(i,j)
        enddo
    enddo
    !$acc end kernels

    call MPI_Barrier(comm2d, rc)

    call MPI_ALLreduce(error1, total_error1, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(error2, total_error2, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(error5, total_error5, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(error6, total_error6, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    
    errorU = dsqrt(total_error1)/dsqrt(total_error2)
    errorT = total_error5 / total_error6

    if (rank == 0) then
        write(*,*) itc,' ',errorU,' ',errorT
    endif

    return
end subroutine check


subroutine message_passing_f()
    use mpi
    use mpi_data
    use commondata
    implicit none
    integer :: i, j, idx, tmp

    ! ------------ exchange message along x ----------------
    if (dims(0) > 1) then
        ! pack the data
        !$acc parallel loop independent collapse(2) private(tmp)
        do idx = 0, 2
            do j = 0, ny+1
                ! f_post(nx,j,dir)/f_post(1,j,dir)是竖直边界数据，在Fortran内存中随j变化为跨步访问。
                ! tmp把这些非连续边界点按idx分段打包成连续的一维send_pos/send_neg缓冲区，供MPI_Sendrecv使用。
                tmp = (ny+2)*idx + j + 1
                send_pos(tmp) = f_post(nx, j, f_tag_x_pos(idx))
                send_neg(tmp) = f_post(1, j, f_tag_x_neg(idx))
            enddo
        enddo

        ! host_data表示下面的MPI调用仍在CPU端执行；use_device表示传入这些数组的GPU设备端地址。
        ! send_pos/recv_pos/send_neg/recv_neg分别是正/负方向的发送和接收缓冲区。
        ! 若MPI支持GPU-aware通信，MPI_Sendrecv可直接使用GPU缓冲区，避免GPU-CPU来回拷贝。
        !$acc host_data use_device(send_pos, recv_pos, send_neg, recv_neg)
        ! message passing to right(i++)
        call MPI_Sendrecv(send_pos, 3*(ny+2), MPI_DOUBLE_PRECISION, nbr_right, 158, &
            recv_pos, 3*(ny+2), MPI_DOUBLE_PRECISION, nbr_left, 158, &
            comm2d, MPI_STATUS_IGNORE, rc)

        ! message passing to left(i--)
        call MPI_Sendrecv(send_neg, 3*(ny+2), MPI_DOUBLE_PRECISION, nbr_left, 367, &
            recv_neg, 3*(ny+2), MPI_DOUBLE_PRECISION, nbr_right, 367, &
            comm2d, MPI_STATUS_IGNORE, rc)   
        !$acc end host_data  
                
        ! unpack the data 
        !$acc parallel loop independent collapse(2) private(tmp)
        do idx = 0, 2
            do j = 0, ny+1
                tmp = (ny+2)*idx + j + 1
                f_post(0, j, f_tag_x_pos(idx)) = recv_pos(tmp)
                f_post(nx+1, j, f_tag_x_neg(idx)) = recv_neg(tmp)
            enddo
        enddo      
    endif

    ! ------------ exchange message along y ----------------
    if(dims(1) > 1) then
        ! pack the data
        !$acc parallel loop independent collapse(2) private(tmp)
        do idx = 0, 2
            do i = 0, nx+1
                tmp = (nx+2)*idx + i + 1
                send_pos(tmp) = f_post(i, ny, f_tag_y_pos(idx))
                send_neg(tmp) = f_post(i, 1, f_tag_y_neg(idx))
            enddo
        enddo

        !$acc host_data use_device(send_pos, recv_pos, send_neg, recv_neg)
        ! message passing to top(j++)
        call MPI_Sendrecv(send_pos, 3*(nx+2), MPI_DOUBLE_PRECISION, nbr_top, 256, &
            recv_pos, 3*(nx+2), MPI_DOUBLE_PRECISION, nbr_bottom, 256, &
            comm2d, MPI_STATUS_IGNORE, rc)
        ! message passing to bottom(j--)
            call MPI_Sendrecv(send_neg, 3*(nx+2), MPI_DOUBLE_PRECISION, nbr_bottom, 478, &
            recv_neg, 3*(nx+2), MPI_DOUBLE_PRECISION, nbr_top, 478, &
            comm2d, MPI_STATUS_IGNORE, rc)
        !$acc end host_data

        ! unpack the data
        !$acc parallel loop independent collapse(2) private(tmp)
        do idx = 0, 2
            do i = 0, nx+1
                tmp = (nx+2)*idx + i + 1
                f_post(i, 0, f_tag_y_pos(idx)) = recv_pos(tmp)
                f_post(i, ny+1, f_tag_y_neg(idx)) = recv_neg(tmp)
            enddo
        enddo    

    endif

    ! message exchange with corners implicitly happens 
                
end subroutine message_passing_f


subroutine message_passing_g()
    use mpi
    use mpi_data
    use commondata
    implicit none
    integer :: i, j

    ! ------------ exchange message along x ----------------
    if (dims(0) > 1) then
        ! pack the data send to right
        !$acc parallel loop independent
        do j = 0, ny+1
            g_send_pos_x(j+1) = g_post(nx, j, 1)
            g_send_neg_x(j+1) = g_post(1, j, 3)
        enddo

        !$acc host_data use_device(g_send_pos_x, g_recv_pos_x, g_send_neg_x, g_recv_neg_x)
        ! message passing to right(i++)
        call MPI_Sendrecv(g_send_pos_x, ny+2, MPI_DOUBLE_PRECISION, nbr_right, 1, &
            g_recv_pos_x, ny+2, MPI_DOUBLE_PRECISION, nbr_left, 1, &
            comm2d, MPI_STATUS_IGNORE, rc)

        ! message passing to left(i--)
        call MPI_Sendrecv(g_send_neg_x, ny+2, MPI_DOUBLE_PRECISION, nbr_left, 3, &
            g_recv_neg_x, ny+2, MPI_DOUBLE_PRECISION, nbr_right, 3, &
            comm2d, MPI_STATUS_IGNORE, rc)   
        !$acc end host_data  
                    
        ! depack the data from right
        !$acc parallel loop independent
        do j = 0, ny+1
            g_post(0, j, 1) = g_recv_pos_x(j+1)
            g_post(nx+1, j, 3) = g_recv_neg_x(j+1)
        enddo   
    endif

    ! ------------ exchange message along y ----------------
    if (dims(1) > 0) then
        !$acc host_data use_device(g_post)
        ! message passing to top(j++)
        call MPI_Sendrecv(g_post(0, ny, 2), nx+2, MPI_DOUBLE_PRECISION, nbr_top, 2, &
            g_post(0, 0, 2), nx+2, MPI_DOUBLE_PRECISION, nbr_bottom, 2, &
            comm2d, MPI_STATUS_IGNORE, rc)

        ! message passing to bottom(j--)
        call MPI_Sendrecv(g_post(0, 1, 4), nx+2, MPI_DOUBLE_PRECISION, nbr_bottom, 4, &
            g_post(0, ny+1, 4), nx+2, MPI_DOUBLE_PRECISION, nbr_top, 4, &
            comm2d, MPI_STATUS_IGNORE, rc)
        !$acc end host_data    
    endif

    ! message exchange with corners implicitly happens 

end subroutine message_passing_g


#ifdef __OUT_PUT__
    subroutine output()
        use mpi
        use commondata
        use mpi_data
        integer :: i, j
        integer :: p_rank, num(0:3) ,dx = 0, dy = 0, new_coords(0:1)
        real(8), allocatable :: total_u(:, :), total_v(:, :), total_rho(:, :), total_T(:, :)
        real(8), allocatable :: tmp_u(:, :), tmp_v(:, :), tmp_rho(:, :), tmp_T(:, :)

        ! update self表示把设备端(GPU)的u/v/rho/T同步回主机端(CPU)，供后续MPI发送和文件输出使用。
        ! if_present表示只有这些变量已存在于OpenACC设备端时才更新，避免变量不在设备端时报错。
        !$acc update if_present self(u,v,rho,T)

        if (rank2d > 0) then
            ! rank != 0 send data
            num(0) = nx
            num(1) = ny
            num(2) = i_start_global
            num(3) = j_start_global
            ! send to rank 0
            ! MPI_Send参数含义：num为发送缓冲区；4表示发送4个整数；
            ! MPI_INTEGER为数据类型；第一个0为目标rank 0；第二个0为消息tag；
            ! comm2d为二维笛卡尔通信器；rc为返回状态码。
            ! num(0:3)依次保存局部块尺寸nx/ny和全局起点i_start_global/j_start_global。
            call MPI_Send(num, 4, MPI_INTEGER, 0, 0, comm2d, rc)    ! block size and origion
            call MPI_Send(u, nx*ny, MPI_REAL8, 0, 1, comm2d, rc)
            call MPI_Send(v, nx*ny, MPI_REAL8, 0, 2, comm2d, rc)
            call MPI_Send(rho, nx*ny, MPI_REAL8, 0, 3, comm2d, rc)
            call MPI_Send(T, nx*ny, MPI_REAL8, 0, 4, comm2d, rc)
        else
            ! rank 0 collect data
            allocate(total_u(total_nx, total_ny))
            allocate(total_v(total_nx, total_ny))
            allocate(total_rho(total_nx, total_ny))
            allocate(total_T(total_nx, total_ny))

            dx = i_start_global
            dy = j_start_global

            ! collect data from rank 0
            do j = 1, ny
                do i = 1, nx
                    total_u(dx + i, dy + j) = u(i, j)
                    total_v(dx + i, dy + j) = v(i, j)
                    total_rho(dx + i, dy + j) = rho(i, j)
                    total_T(dx + i, dy + j) = T(i, j)
                enddo
            enddo

            ! collect data from all other processors
            do p_rank = 1, dims(0) * dims(1) - 1
                ! receive the block size and origion
                call MPI_Recv(num, 4, MPI_INTEGER, p_rank, 0, comm2d, MPI_STATUS_IGNORE, rc)
                ! creat buffer
                allocate(tmp_u(num(0), num(1)))
                allocate(tmp_v(num(0), num(1)))
                allocate(tmp_rho(num(0), num(1)))
                allocate(tmp_T(num(0), num(1)))
                
                ! receive data
                call MPI_Recv(tmp_u, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 1, comm2d, MPI_STATUS_IGNORE, rc)
                call MPI_Recv(tmp_v, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 2, comm2d, MPI_STATUS_IGNORE, rc)
                call MPI_Recv(tmp_rho, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 3, comm2d, MPI_STATUS_IGNORE, rc)
                call MPI_Recv(tmp_T, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 4, comm2d, MPI_STATUS_IGNORE, rc)

                ! determine the origin
                dx = num(2)
                dy = num(3)

                ! assign data
                do j = 1, num(1)
                    do i = 1, num(0)
                        total_u(dx + i, dy + j) = tmp_u(i, j)
                        total_v(dx + i, dy + j) = tmp_v(i, j)
                        total_rho(dx + i, dy + j) = tmp_rho(i, j)
                        total_T(dx + i, dy + j) = tmp_T(i, j)
                    enddo
                enddo

                deallocate(tmp_u)
                deallocate(tmp_v)
                deallocate(tmp_rho)
                deallocate(tmp_T)
            enddo
            
            call output_Tecplot(xp, yp, total_u, total_v, total_rho, total_T, total_nx, total_ny, itc)
            call output_binary(total_u, total_v, total_rho, total_T, total_nx, total_ny, itc)
            ! call output_ASCII(total_u, total_v, total_rho, total_T, total_nx, total_ny, itc)

            deallocate(total_u)
            deallocate(total_v)
            deallocate(total_rho)
            deallocate(total_T)

        endif

    end subroutine output

    subroutine output_binary(u, v, rho, T, nx, ny, itc)
        implicit none
        integer, intent(in) :: nx, ny, itc
        real(8), intent(in) :: u(nx, ny), v(nx, ny), rho(nx, ny), T(nx, ny)
        integer(kind=4) :: i, j
        character(len=100) :: filename


        write(filename,*) itc
        filename = adjustl(filename)

        open(unit=03,file="buoyancyCavity-"//trim(filename)//'.bin',form="unformatted",access="sequential")
        write(03) ((u(i,j),i=1,nx),j=1,ny)
        write(03) ((v(i,j),i=1,nx),j=1,ny)
        write(03) ((T(i,j),i=1,nx),j=1,ny)
        close(03)

        return
    end subroutine output_binary

    subroutine output_Tecplot(xp, yp, u, v, rho, T, nx, ny, itc)
        implicit none
        integer, intent(in) :: nx, ny, itc
        real(8), intent(in) :: xp(0:nx+1), yp(0:ny+1)
        real(8), intent(in) :: u(nx, ny), v(nx, ny), rho(nx, ny), T(nx, ny)
        integer(kind=4) :: i, j, k
        character(len=9) :: B2
        REAL(kind=4) :: zoneMarker, eohMarker
        character(len=40) :: title
        character(len=40) :: V1,V2,V3,V4,V5,V6,V7,V8
        integer(kind=4), parameter :: kmax=1
        character(len=40) :: zoneName
        character(len=100) :: filename


        write(filename,'(i12.12)') itc
        filename = adjustl(filename)
        
        open(unit=41,file="buoyancyCavity-"//trim(filename)//'.plt', access='stream', form='unformatted')

        !---------------------------------------------
        zoneMarker= 299.0
        eohMarker = 357.0

        !I. HEAD SECTION--------------------------------------
        !c--Magic number, Version number
        write(41) '#!TDV101'

        !c--Integer value of 1
        write(41) 1

        Title='MyFirst'
        call dumpstring(title)

        !c-- Number of variables in this data file

        write(41) 6

        !c-- Variable names.
        V1='X'
        call dumpstring(V1)
        V2='Y'
        call dumpstring(V2)
        V3='U'
        call dumpstring(V3)
        V4='V'
        call dumpstring(V4)
        V5='T'
        call dumpstring(V5)
        V6='rho'
        call dumpstring(V6)



        !c-----Zones-----------------------------

        !c--------Zone marker. Value = 299.0
        write(41) zoneMarker

        !--------Zone name.
        zoneName='ZONE 001'
        call dumpstring(zoneName)

        !---------Zone Color
        write(41) -1

        !---------ZoneType
        write(41) 0

        !---------DataPacking 0=Block, 1=Point
        write(41) 1

        !---------Specify Var Location. 0 = Do not specify, all data
        !---------is located at the nodes. 1 = Specify
        write(41) 0

        !---------Number of user defined face neighbor connections
        ! (value >= 0)
        write(41) 0

        !---------IMax,JMax,KMax
        write(41) nx
        write(41) ny
        write(41) kmax

        !-----------1=Auxiliary name/value pair to follow
        !-----------0=No more Auxiliar name/value pairs.
        write(41) 0
        write(41) eohMarker

        !----zone ------------------------------------------------------------
        write(41) zoneMarker

        !--------variable data format, 1=Float, 2=Double, 3=LongInt,4=ShortInt, 5=Byte, 6=Bit
        write(41) 1
        write(41) 1
        write(41) 1
        write(41) 1
        write(41) 1
        write(41) 1


        !--------Has variable sharing 0 = no, 1 = yes.
        write(41) 0

        !----------Zone number to share connectivity list with (-1 = no
        ! sharing).
        write(41) -1

        !---------------------------------------------------------------------
        do k=1,kmax
            do j=1,ny
                do i=1,nx
                    write(41) real(xp(i))
                    write(41) real(yp(j))
                    write(41) real(u(i,j))
                    write(41) real(v(i,j))
                    write(41) real(T(i,j))
                    write(41) real(rho(i,j))
                end do
            end do
        enddo
        close(41)
        !---------------------------------------------------------------------

        return
    end subroutine output_Tecplot

    subroutine dumpstring(instring)
        implicit none
        character(len=40) instring
        integer(kind=4) :: stringLength
        integer(kind=4) :: ii
        integer(kind=4) :: I

        stringLength=LEN_TRIM(instring)
        do ii=1,stringLength
            I=ICHAR(instring(ii:ii))
            write(41) I
        end do
        write(41) 0

        return
    end subroutine dumpstring

#endif
