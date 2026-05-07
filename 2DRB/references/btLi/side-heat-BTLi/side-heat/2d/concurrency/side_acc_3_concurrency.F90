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


    !$acc declare create(u, v, T, rho, up, vp, Tp, f, f_post, g, g_post, Fx, Fy) &
    !$acc create(nx, ny)
end module commondata
        
module mpi_data
    implicit none
    integer :: rc, rank, num_process
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
    character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name

    call MPI_Init(rc)

    call MPI_Comm_size(MPI_COMM_WORLD, num_process, rc)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, rc)
    call MPI_Get_processor_name(processor_name, name_len, rc)

    !!! decomposition the domain 
    ! call MPI_Dims_create(num_process, 2, dims, rc)
    ! dims(0) = 1
    call MPI_Dims_create_2d(num_process, dims, total_nx, total_ny, periods, rc)
    
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, .true., comm2d, rc)
    if(rank == 0) then
        write(*,*) "dimens is x*y = ", dims(0), "x", dims(1)
    endif

    ! get my new rank in decomposition
    call MPI_Comm_rank(comm2d, rank2d, rc)
    ! write(*,*) "process ", rank2d, " of total ", num_process, "is alive."

    ! determain sub-domain size
    call MPI_Cart_get(comm2d, 2, dims, periods, coords, rc)
    call decompose_1d(total_nx, nx, coords(0), dims(0), i_start_global)
    call decompose_1d(total_ny, ny, coords(1), dims(1), j_start_global)
    if(rank == 0) then
        write(*,*) "Total nx * ny =", total_nx, " * ", total_ny
        write(*,*) "local nx * ny =", nx, " * ", ny
    endif
    ! write(*,*) "coords = ", coords(1), coords(2)
    ! write(*,*) "nx*ny = ", nx, ny

    ! get the neighbors
    call MPI_Cart_shift(comm2d, 0, 1, nbr_left, nbr_right, rc)
    call MPI_Cart_shift(comm2d, 1, 1, nbr_bottom, nbr_top, rc)

    local_rank = -1
    call MPI_Comm_split_type(comm2d, MPI_COMM_TYPE_SHARED, rank2d, MPI_INFO_NULL, local_comm, rc)
    call MPI_Comm_rank(local_comm, local_rank, rc)
    call MPI_Comm_free(local_comm, rc)
    
    num_gpus = 0
    num_gpus = acc_get_num_devices(acc_device_nvidia)
    if (local_rank .eq. 0) then
        write(*,*) "I from rank", rank2d, "we have", num_gpus, "gpus"
    endif
    if (num_gpus .le. 0) then
        if (rank2d .eq. 0) then
            write(*,*) 'No NVIDIA GPUs available'
            call MPI_Abort(MPI_COMM_WORLD, 1, rc)
        endif
        else
            gpu_id = mod(local_rank, num_gpus)
            write(*,*) "i'm local rank", local_rank, "rank", rank2d, "using gpu ", gpu_id
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

        restrict = dims

        diff = dble(total_nx) * dble(total_ny) * 10.0d0

        is = 1
        ie = num_process
        ! if user set restrictions
        do i = 0, 1
            if (restrict(i) .NE. 0) then
                is(i) = restrict(i)
                ie(i) = restrict(i)
            endif
        enddo
        

        do i = is(0), ie(0)
            do j = is(1), ie(1)
                if (i * j == num_process) then
                    message = 0.0d0
                    ! local nx, local ny, local nz
                    lx = dble(total_nx) / dble(i)
                    ly = dble(total_ny) / dble(j)

                    ! maximum message need to exchange for a process
                    if (i > 1) then
                        ! if divide at this dimision
                        message = message + ly
                        if (i > 2 .OR. periods(0)) then
                            ! if dims() > 2 or is periodic
                            message = message + ly 
                        endif
                    endif
                    if (j > 1) then
                        message = message + lx
                        if (j > 2 .OR. periods(1)) then
                            message = message + lx 
                        endif
                    endif


                    if (message < diff) then
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

    !$acc kernels
#ifdef VerticalWallsConstT
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
    !$acc routine gang nohost
    use commondata
    implicit none
    integer, intent(in) :: i_start, i_end, j_start, j_end
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: m(0:8), m_post(0:8), meq(0:8)
    real(kind=8) :: s(0:8)
    real(kind=8) :: fSource(0:8)

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
    
    !$acc parallel loop independent collapse(2) async(2)
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

    !$acc kernels async(2)
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

    !$acc parallel loop independent collapse(2) async(2)
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
    !$acc private(ip, jp) async(1)
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

    !$acc kernels async(1)
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
    real(kind=8) :: a
    
    !$acc parallel loop independent collapse(2) async(1)
    do j = 1, ny
        do i = 1, nx
            T(i,j) = g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
        enddo
    enddo

    return
end subroutine macroT


subroutine compute_boundary_f_yz()
    use commondata
    use mpi_data
    use mpi
    implicit none
    integer :: i, j, idx, tmp

    !$acc kernels async(3)
    call collision(1, nx, 1, 1)
    !$acc end kernels

    !$acc kernels async(4)
    call collision(1, nx, ny, ny)
    !$acc end kernels

end subroutine compute_boundary_f_yz

    
subroutine flow_update()
    use commondata
    use mpi_data
    use mpi
    implicit none

    call compute_boundary_f_yz()
    !$acc wait(3,4)

    ! async(3)
    call pack_f_y()
        !$acc kernels async(1)
        call collisionT(1, nx, 1, ny)
        !$acc end kernels
        !$acc wait(3)
        call mpi_f_y()
    ! async(3)
    call unpack_f_y()
    !$acc wait

    ! async(3)
    call pack_g_x()
    ! async(4)
    call pack_g_y()
        !$acc kernels async(2)
        call collision(1, nx, 2, ny-1)
        !$acc end kernels
        !$acc wait(3, 4)
        call mpi_g_x()
        call mpi_g_y()
    ! async(3)
    call unpack_g_x()
    ! async(4)
    call unpack_g_y()
    !$acc wait

    ! async(3)
    call pack_f_x()
        ! kernels async(1)
        call streamingT()
        call bouncebackT()
        call macroT()
        ! end kernels
        !$acc wait(3)
        call mpi_f_x()
    ! async(2)
    call unpack_f_x()

    ! async(2)
    call streaming()
    call bounceback()
    call macro()
    
    !$acc wait

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


subroutine mpi_f_y()
    use mpi
    use mpi_data
    use commondata
    implicit none

    ! ------------ exchange message along y ----------------
    if(dims(1) > 1) then
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
    endif

endsubroutine 

subroutine pack_f_y
    use mpi_data
    use commondata
    implicit none
    integer :: i, j, tmp, idx

    if(dims(1) > 1) then
        !$acc parallel loop independent collapse(2) private(tmp) async(3)
        do idx = 0, 2
            do i = 0, nx+1
                tmp = (nx+2)*idx + i + 1
                send_pos(tmp) = f_post(i, ny, f_tag_y_pos(idx))
                send_neg(tmp) = f_post(i, 1, f_tag_y_neg(idx))
            enddo
        enddo        
    endif

endsubroutine

subroutine unpack_f_y()
    use mpi_data
    use commondata
    implicit none
    integer :: i, j, tmp, idx

    if(dims(1) > 1) then
        !$acc parallel loop independent collapse(2) private(tmp) async(3)
        do idx = 0, 2
            do i = 0, nx+1
                tmp = (nx+2)*idx + i + 1
                f_post(i, 0, f_tag_y_pos(idx)) = recv_pos(tmp)
                f_post(i, ny+1, f_tag_y_neg(idx)) = recv_neg(tmp)
            enddo
        enddo      
    endif

endsubroutine
 

subroutine mpi_f_x()
    use mpi
    use mpi_data
    use commondata
    implicit none

    ! ------------ exchange message along x ----------------
    if (dims(0) > 1) then
    ! message passing to right(i++)
    !$acc host_data use_device(send_pos, recv_pos, send_neg, recv_neg)
    call MPI_Sendrecv(send_pos, 3*(ny+2), MPI_DOUBLE_PRECISION, nbr_right, 158, &
        recv_pos, 3*(ny+2), MPI_DOUBLE_PRECISION, nbr_left, 158, &
        comm2d, MPI_STATUS_IGNORE, rc)
    
    ! message passing to left(i--)
    call MPI_Sendrecv(send_neg, 3*(ny+2), MPI_DOUBLE_PRECISION, nbr_left, 367, &
        recv_neg, 3*(ny+2), MPI_DOUBLE_PRECISION, nbr_right, 367, &
        comm2d, MPI_STATUS_IGNORE, rc)   
    !$acc end host_data  
    endif

endsubroutine 

subroutine pack_f_x()
    use mpi
    use mpi_data
    use commondata
    implicit none
    integer :: i, j, idx, tmp

    if (dims(0) > 1) then
        !$acc parallel loop independent collapse(2) private(tmp) async(3)
        do idx = 0, 2
            do j = 0, ny+1
                tmp = (ny+2)*idx + j + 1
                send_pos(tmp) = f_post(nx, j, f_tag_x_pos(idx))
                send_neg(tmp) = f_post(1, j, f_tag_x_neg(idx))
            enddo
        enddo
    endif
    

endsubroutine 

subroutine unpack_f_x()
    use mpi
    use mpi_data
    use commondata
    implicit none
    integer :: i, j, idx, tmp

    if (dims(0) > 1) then
        !$acc parallel loop independent collapse(2) private(tmp) async(2)
        do idx = 0, 2
            do j = 0, ny+1
                tmp = (ny+2)*idx + j + 1
                f_post(0, j, f_tag_x_pos(idx)) = recv_pos(tmp)
                f_post(nx+1, j, f_tag_x_neg(idx)) = recv_neg(tmp)
            enddo
        enddo      
    endif

endsubroutine 


subroutine mpi_g_y()
    use mpi
    use mpi_data
    use commondata
    implicit none

    ! ------------ exchange message along y ----------------
    if (dims(1) > 1) then
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

endsubroutine 

subroutine pack_g_y()
    return
endsubroutine 

subroutine unpack_g_y()
    return
endsubroutine 


subroutine mpi_g_x()
    use mpi
    use mpi_data
    use commondata
    implicit none

    if (dims(0) > 1) then
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
    endif

endsubroutine 

subroutine pack_g_x()
    use mpi
    use mpi_data
    use commondata
    implicit none
    integer :: i, j

    if (dims(0) > 1) then
        !$acc parallel loop independent async(3)
        do j = 0, ny+1
            g_send_pos_x(j+1) = g_post(nx, j, 1)
            g_send_neg_x(j+1) = g_post(1, j, 3)
        enddo
    endif

endsubroutine 

subroutine unpack_g_x()
    use mpi
    use mpi_data
    use commondata
    implicit none
    integer :: i, j

    if (dims(0) > 1) then
        !$acc parallel loop independent async(3)
        do j = 0, ny+1
            g_post(0, j, 1) = g_recv_pos_x(j+1)
            g_post(nx+1, j, 3) = g_recv_neg_x(j+1)
        enddo   
    endif
  

endsubroutine 


#ifdef __OUT_PUT__
    subroutine output()
        use mpi
        use commondata
        use mpi_data
        integer :: i, j
        integer :: p_rank, num(0:3) ,dx = 0, dy = 0, new_coords(0:1)
        real(8), allocatable :: total_u(:, :), total_v(:, :), total_rho(:, :), total_T(:, :)
        real(8), allocatable :: tmp_u(:, :), tmp_v(:, :), tmp_rho(:, :), tmp_T(:, :)

        !$acc update if_present self(u,v,rho,T)

        if (rank2d > 0) then
            ! rank != 0 send data
            num(0) = nx
            num(1) = ny
            num(2) = i_start_global
            num(3) = j_start_global
            ! send to rank 0
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