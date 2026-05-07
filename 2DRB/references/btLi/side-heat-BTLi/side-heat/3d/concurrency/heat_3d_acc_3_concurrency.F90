!!!    This program sloves Buoyancy Driven Cavity Flow problem using Lattice Boltzmann Method
!!!    Lattice Boltzmann Equation with MRT-LBE model

!!==velocity B.C.==
#define noslipWalls
!!==velocity B.C.==

!~~temperature B.C. (for cavity flow benchmark)~~
#define benchmarkCavity
#define BackFrontWallsAdiabatic
#define LeftRightWallsConstT
#define TopBottomPlatesAdiabatic
!~~temperature B.C.~~

module commondata
    implicit none
    
    integer, parameter :: loadInitField=0
    integer, parameter :: flowReversal=0  !! value 0: do NOT simulate reversal; value 1: do simulate reversal

    integer, parameter :: total_nx = 385, total_ny = total_nx, total_nz = total_nx
    integer :: nx, ny, nz, i_start_global, j_start_global, k_start_global
    integer, parameter :: nxHalf = (total_nx-1)/2+1, nyHalf = (total_ny-1)/2+1, nzHalf = (total_nz-1)/2+1
    
    real(kind=8), parameter :: Rayleigh=1e7
    real(kind=8), parameter :: Prandtl=0.71d0
    real(kind=8), parameter :: Mach=0.1d0

    real(kind=8), parameter :: tauf=0.5d0+Mach*dble(total_nz)*DSQRT(3.0d0*Prandtl/Rayleigh)
    real(kind=8), parameter :: Thot=1.0d0, Tcold=0.0d0, Tref=0.0d0
    real(kind=8), parameter :: viscosity=(tauf-0.5d0)/3.0d0
    real(kind=8), parameter :: diffusivity=viscosity/Prandtl
    
    real(kind=8), parameter :: Ekman=0.001d0
    real(kind=8), parameter :: omegaRatating=viscosity/2.0d0/Ekman/dble(total_nz*total_nz)
    
    real(kind=8), parameter :: paraA=42.0d0*dsqrt(3.0d0)*diffusivity-6.0d0
    real(kind=8), parameter :: gBeta1=Rayleigh*viscosity*diffusivity/dble(total_nz)
    real(kind=8), parameter :: gBeta=gBeta1/dble(total_nz*total_nz)
    
    real(kind=8), parameter :: timeUnit=dsqrt(dble(total_nz)/gBeta)  !!dble(ny*ny)/diffusivity
    integer, parameter :: dimensionlessTimeMax=1000
    integer, parameter :: flowReversalTime=20000
    integer :: itc
    integer, parameter :: itc_max=INT(dimensionlessTimeMax*timeUnit)
    
    real(kind=8), parameter :: epsU=1e-8
    real(kind=8), parameter :: epsT=1e-8
    real(kind=8) :: errorU, errorT
    
    real(kind=8) :: xp(0:total_nx+1), yp(0:total_ny+1), zp(0:total_nz+1)
    real(kind=8), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:), T(:,:,:)
    real(kind=8), allocatable :: rho(:,:,:)
    real(kind=8), allocatable :: up(:,:,:), vp(:,:,:), wp(:,:,:), Tp(:,:,:)

    real(kind=8), allocatable :: f(:,:,:,:), f_post(:,:,:,:)
    real(kind=8), allocatable :: g(:,:,:,:), g_post(:,:,:,:)
    
    integer, parameter :: ex(0:18) = (/ 0,  &
                                        1, -1,  0,  0,  0,  0, &
                                        1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0 /)
    integer, parameter :: ey(0:18) = (/ 0, &
                                        0,  0,  1, -1,  0,  0, &
                                        1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1 /)  
    integer, parameter :: ez(0:18) = (/ 0, &
                                        0,  0,  0,  0,  1, -1, &
                                        0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1 /)
    
    real(kind=8), allocatable :: Fx(:,:,:), Fy(:,:,:), Fz(:,:,:)
    
    real(kind=8), parameter :: Snu=1.0d0/tauf, Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)
    real(kind=8), parameter :: Qd=3.0d0-dsqrt(3.0d0), Qnu=4.0d0*dsqrt(3.0d0)-6.0d0    

    !$acc declare create(u, v, w, T, rho, up, vp, wp, Tp, f, f_post, g, g_post, nx, ny, nz, Fx, Fy, Fz)

end module commondata


module mpi_data
    implicit none
    integer :: rc, rank, num_process 
    integer :: dims(0:2) = (/ 0, 0, 0 /), coords(0:2)
    logical :: periods(0:2) = (/ .false., .false., .false. /)
    integer :: comm3d, rank3d
    integer :: nbr_surface(1:6) 
    integer :: surface_length_x, surface_length_y, surface_length_z
    ! buffer layer 
    real(8), allocatable :: send_pos(:), recv_pos(:), send_neg(:), recv_neg(:)
    real(8), allocatable :: g_send_pos_x(:), g_recv_pos_x(:), g_send_neg_x(:), g_recv_neg_x(:)
    real(8), allocatable :: g_send_pos_y(:), g_recv_pos_y(:), g_send_neg_y(:), g_recv_neg_y(:)
    ! real(8), allocatable :: g_send_pos_z(:), g_recv_pos_z(:), g_send_neg_z(:), g_recv_neg_z(:)
    ! f tags at different directions
    integer, parameter :: f_tag_x_pos(0:4) = (/ 1, 7, 9, 11, 13 /)
    integer, parameter :: f_tag_x_neg(0:4) = (/ 2, 8, 10 ,12, 14 /)
    integer, parameter :: f_tag_y_pos(0:4) = (/ 3, 7, 8, 15, 17 /)
    integer, parameter :: f_tag_y_neg(0:4) = (/ 4, 9, 10, 16, 18 /)
    integer, parameter :: f_tag_z_pos(0:4) = (/ 5, 11, 12, 15, 16 /)
    integer, parameter :: f_tag_z_neg(0:4) = (/ 6, 13, 14, 17, 18 /)

    !$acc declare create(send_pos(:), recv_pos(:), send_neg(:), recv_neg(:)) &
    !$acc create(g_send_pos_x(:), g_recv_pos_x(:), g_send_neg_x(:), g_recv_neg_x(:)) &
    !$acc create(g_send_pos_y(:), g_recv_pos_y(:), g_send_neg_y(:), g_recv_neg_y(:)) &
    !$acc create(surface_length_x, surface_length_y, surface_length_z, dims, coords)
end module mpi_data

! #define __OUT_PUT__

program main
    use mpi    
    use mpi_data
    use commondata
    implicit none
    real(kind=8) :: start, finish
    real(8) :: start_time, end_time

    call mpi_starts()

    !$acc update device(surface_length_x, surface_length_y, surface_length_z, dims, coords, nx, ny, nz)

    call allocate_all()

    call initial()

#ifdef __OUT_PUT__
    call output()
#endif

    call CPU_TIME(start)

    call MPI_Barrier(MPI_COMM_WORLD, rc)
    start_time = MPI_Wtime()

    do while( ((errorU.GT.epsU).OR.(errorT.GT.epsT)).AND.(itc.LE.itc_max) )

        itc = itc+1

        call flow_update()

        if(MOD(itc,2000).EQ.0) then
            call check()
        endif

        ! timer test
        if (mod(itc, 6000) == 0) then
            exit
        endif
    enddo

    call MPI_Barrier(MPI_COMM_WORLD, rc)
    end_time = MPI_Wtime()
    
    call CPU_TIME(finish)

    if (rank == 0) then
        write(*,*) "---------------------------------------------"
        write(*,*) "Time (CPU) = ", real(finish-start), "s"
        write(*,*) "Time (MPI) = ", real(end_time - start_time), "s"
        write(*,*) "---------------------------------------------"
    endif

#ifdef __OUT_PUT__
    call output()
#endif

    if (rank == 0) then
        write(*,*) "Deallocate Array..."
    endif

    call free_all()

    if (rank == 0) then
        write(*,*) "    "
        write(*,*) "Successfully: DNS completed!"
    endif

    call MPI_Finalize(rc)

end program main

subroutine mpi_starts()
    use commondata
    use mpi_data
    use mpi
    use openacc
    implicit none
    integer :: num_gpus, gpu_id
    integer :: local_comm, local_rank
    integer :: name_len, tmp
    character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name

    call MPI_Init(rc)

    call MPI_Comm_size(MPI_COMM_WORLD, num_process, rc)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, rc)
    call MPI_Get_processor_name(processor_name, name_len, rc)

    !!! ---- Decomposition the domain
    ! call MPI_Dims_create(num_process, 3, dims, rc)
    call MPI_Dims_create_3d(num_process, dims, total_nx, total_ny, total_nz, periods, rc)


    call MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, .true., comm3d, rc)
    if(rank == 0) then
        write(*,*) "Using ", num_process, "processers."
        write(*,*) "dimens is x*y*z = ", dims(0), "x", dims(1), "x", dims(2)
    endif

    ! get my new rank in decomposition
    call MPI_Comm_rank(comm3d, rank3d, rc)
    ! write(*,*) "process ", rank3d, " of total ", num_process, "is alive."

    ! determain sub-domain size
    call MPI_Cart_get(comm3d, 3, dims, periods, coords, rc)
    call decompose_1d(total_nx, nx, coords(0), dims(0), i_start_global)
    call decompose_1d(total_ny, ny, coords(1), dims(1), j_start_global)
    call decompose_1d(total_nz, nz, coords(2), dims(2), k_start_global)
    surface_length_x = (ny+2) * (nz+2)
    surface_length_y = (nx+2) * (nz+2)
    surface_length_z = (nx+2) * (ny+2)
    ! write(*,*) "coords = ", coords(1), coords(2), coords(3), "nx*ny = ", nx, ny, nz
    if(rank == 0) then
        write(*,*) "Total nx * ny *nz =", total_nx, " * ", total_ny, " * ", total_nz
        write(*,*) "local nx * ny *nz =", nx, " * ", ny, " * ", nz
    endif

    ! get the neighbors
    call MPI_Cart_shift(comm3d, 0, 1, nbr_surface(2), nbr_surface(1), rc)
    call MPI_Cart_shift(comm3d, 1, 1, nbr_surface(4), nbr_surface(3), rc)
    call MPI_Cart_shift(comm3d, 2, 1, nbr_surface(6), nbr_surface(5), rc)
    ! write(*,*) "I'm process ", rank3d, "My neighbor surfaces are", nbr_surface(1), nbr_surface(2), nbr_surface(3), nbr_surface(4), nbr_surface(5), nbr_surface(6)

    local_rank = -1
    call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank3d, MPI_INFO_NULL, local_comm, rc)
    call MPI_Comm_rank(local_comm, local_rank, rc)
    call MPI_Comm_free(local_comm, rc)
    
    num_gpus = 0
    num_gpus = acc_get_num_devices(acc_device_nvidia)
    if (local_rank .eq. 0) then
        write(*,*) "I from rank", rank3d, "we have", num_gpus, "gpus"
    endif
    if (num_gpus .le. 0) then
        if (rank3d .eq. 0) then
            write(*,'(a)'), 'No NVIDIA GPUs available'
            call MPI_Abort(MPI_COMM_WORLD, 1, rc)
        endif
        else
            gpu_id = mod(local_rank, num_gpus)
            write(*,*) "i'm local rank", local_rank, "rank", rank3d, "using gpu ", gpu_id
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

    subroutine MPI_Dims_create_3d(num_process, dims, total_nx, total_ny, total_nz, periods, rc)
        integer, intent(in) :: num_process, total_nx, total_ny, total_nz, rc
        logical, intent(in) :: periods(0:2)
        integer, intent(inout) :: dims(0:2)
        integer :: i, j, k, is(0:2), ie(0:2)
        real :: message, diff, lx, ly, lz
        integer :: restrict(0:2)

        ! determine the dimensions of cartesian topologies to minimize the message exchange
        ! under user provided restrictions if dims(0:2) != 0

        restrict = dims

        diff = dble(total_nx) * dble(total_ny) * dble(total_nz)

        is = 1
        ie = num_process
        ! if user set restrictions
        do i = 0, 2
            if (restrict(i) .NE. 0) then
                is(i) = restrict(i)
                ie(i) = restrict(i)
            endif
        enddo
        

        do i = is(0), ie(0)
            do j = is(1), ie(1)
                do k = is(2), ie(2)
                    if (i * j * k == num_process) then
                        message = 0.0d0
                        ! local nx, local ny, local nz
                        lx = dble(total_nx) / dble(i)
                        ly = dble(total_ny) / dble(j)
                        lz = dble(total_nz) / dble(k)

                        ! maximum message need to exchange for a process
                        if (i > 1) then
                            ! if divide at this dimision
                            message = message + ly * lz
                            if (i > 2 .OR. periods(0)) then
                                ! if dims() > 2 or is periodic
                                message = message + ly * lz
                            endif
                        endif
                        if (j > 1) then
                            message = message + lx * lz
                            if (j > 2 .OR. periods(1)) then
                                message = message + lx * lz
                            endif
                        endif
                        if (k > 1) then
                            message = message + lx * ly
                            if (k > 2 .OR. periods(2)) then
                                message = message + lx * ly
                            endif
                        endif

                        if (message < diff) then
                            diff = message
                            dims(0) = i
                            dims(1) = j
                            dims(2) = k
                        endif
                    endif
                enddo
            enddo
        enddo

    end subroutine MPI_Dims_create_3d

end subroutine mpi_starts


subroutine allocate_all()
    use commondata
    use mpi_data
    integer :: max_length

    allocate (u(nx,ny,nz))
    allocate (v(nx,ny,nz))
    allocate (w(nx,ny,nz))
    allocate (T(nx,ny,nz))
    allocate (rho(nx,ny,nz))
    allocate (up(nx,ny,nz))
    allocate (vp(nx,ny,nz))
    allocate (wp(nx,ny,nz))
    allocate (Tp(nx,ny,nz))
    
    allocate (f(nx,ny,nz,0:18))
    allocate (f_post(0:nx+1,0:ny+1,0:nz+1,0:18))
    allocate (g(nx,ny,nz,0:6))
    allocate (g_post(0:nx+1,0:ny+1,0:nz+1,0:6))
    
    allocate (Fx(nx,ny,nz))
    allocate (Fy(nx,ny,nz))
    allocate (Fz(nx,ny,nz))

    max_length = surface_length_x
    if (max_length < surface_length_y) then
        max_length = surface_length_y
    endif
    if (max_length < surface_length_z) then
        max_length = surface_length_z
    endif

    allocate(send_pos(1 : 5*max_length))
    allocate(recv_pos(1 : 5*max_length))
    allocate(send_neg(1 : 5*max_length))
    allocate(recv_neg(1 : 5*max_length))

    allocate(g_send_pos_x(max_length))
    allocate(g_recv_pos_x(max_length))
    allocate(g_send_neg_x(max_length))
    allocate(g_recv_neg_x(max_length))

    allocate(g_send_pos_y(max_length))
    allocate(g_recv_pos_y(max_length))
    allocate(g_send_neg_y(max_length))
    allocate(g_recv_neg_y(max_length))

    ! allocate(g_send_pos_z(max_length))
    ! allocate(g_recv_pos_z(max_length))
    ! allocate(g_send_neg_z(max_length))
    ! allocate(g_recv_neg_z(max_length))

end subroutine

subroutine free_all()
    use commondata
    use mpi_data

    deallocate(f)
    deallocate(g)
    deallocate(f_post)
    deallocate(g_post)
    deallocate(u)
    deallocate(v)
    deallocate(w)
    deallocate(T)
    deallocate(up)
    deallocate(vp)
    deallocate(wp)
    deallocate(Tp)
    deallocate(rho)
    deallocate(Fx)
    deallocate(Fy)
    deallocate(Fz)

    deallocate(send_pos)
    deallocate(recv_pos)
    deallocate(send_neg)
    deallocate(recv_neg)

    deallocate(g_send_pos_x)
    deallocate(g_recv_pos_x)
    deallocate(g_send_neg_x)
    deallocate(g_recv_neg_x)

    deallocate(g_send_pos_y)
    deallocate(g_recv_pos_y)
    deallocate(g_send_neg_y)
    deallocate(g_recv_neg_y)

    ! deallocate(g_send_pos_z)
    ! deallocate(g_recv_pos_z)
    ! deallocate(g_send_neg_z)
    ! deallocate(g_recv_neg_z)

end subroutine


subroutine initial()
    use mpi_data
    use commondata
    implicit none
    integer :: i, j, k
    integer :: alpha
    real(kind=8) :: un(0:18), unT(0:6)
    real(kind=8) :: us2
    real(kind=8) :: omega(0:18), omegaT(0:6)

    itc = 0
    errorU = 100.0d0
    errorT = 100.0d0
    
    if (rank == 0) then
        write(*,*) 'Mesh:',nx,ny,nz
        write(*,*) 'Rayleigh=',real(Rayleigh), ', Prandtl=',real(Prandtl), ', Mach=',real(Mach)
        write(*,*) "Ekman=",real(Ekman)
        write(*,*) "   "
        write(*,*) 'tauf=',real(tauf)
        write(*,*) "paraA=",real(paraA)
        write(*,*) "viscosity=",real(viscosity), ", diffusivity=",real(diffusivity)
        write(*,*) "omegaRatating=", real(omegaRatating)
        write(*,*) "itc_max=",itc_max
        write(*,*) "Ouput will begin at", int(timeUnit)
        write(*,*) "Output interval is", int(timeUnit)
        write(*,*) "Time unit: Sqrt(L0/(gBeta*DeltaT))=",dsqrt(dble(total_ny)/gBeta)
        write(*,*) "Velocity unit: Sqrt(gBeta*L0*DeltaT)=",dsqrt(gBeta*dble(total_ny))
        write(*,*) "    "
    endif

    xp(0) = 0.0d0
    xp(total_nx+1) = dble(total_nx)
    do i = 1, total_nx
        xp(i) = dble(i)-0.5d0
    enddo
    yp(0) = 0.0d0
    yp(total_ny+1) = dble(total_ny)
    do j = 1, total_ny
        yp(j) = dble(j)-0.5d0
    enddo
    zp(0) = 0.0d0
    zp(total_nz+1) = dble(nz)
    do k = 1, total_nz
        zp(k) = dble(k)-0.5d0
    enddo

    omega(0) = 1.0d0/3.0d0
    do alpha=1,6
        omega(alpha) = 1.0d0/18.0d0
    enddo
    do alpha=7,18
        omega(alpha) = 1.0d0/36.0d0
    enddo
    
    omegaT(0) = (1.0d0-paraA)/7.0d0
    do alpha=1,6
        omegaT(alpha) = (paraA+6.0d0)/42.0d0
    enddo


    !$acc parallel loop collapse(3)
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                rho(i, j, k) = 1.0d0
                u(i, j, k) = 0.0d0
                v(i, j, k) = 0.0d0
                w(i, j, k) = 0.0d0
                T(i, j, k) = 0.0d0
                up(i, j, k) = 0.0d0
                vp(i, j, k) = 0.0d0
                wp(i, j, k) = 0.0d0
                Tp(i, j, k) = 0.0d0
            enddo
        enddo
    enddo
    

    
    if(loadInitField.EQ.0) then 
    !$acc kernels
#ifdef LeftRightWallsConstT
    ! left wall y=1
    if (coords(1) == 0) then
        !$acc loop collapse(2)
        do k=1,nz
            do i=1,nx
                T(i,1,k) = Thot
            enddo
        enddo
    endif

    ! right wall y = ny
    if (coords(1) == dims(1) - 1) then
        !$acc loop collapse(2)
        do k=1,nz
            do i=1,nx
                T(i,ny,k) = Tcold
            enddo
        enddo
    endif

#endif
#ifdef TopBottomPlatesConstT
    ! bottom wall k = 1
    if (coords(2) == 0) then
        !$acc loop collapse(2)
        do j=1,ny
            do i=1,nx
                T(i,j,1) = Thot
            enddo
        enddo
    endif


    ! top wall k = nz
    if (coords(2) == dims(2) - 1) then
        !$acc loop collapse(2)
        do j=1,ny
            do i=1,nx
                T(i,j,nz) = Tcold
            enddo
        enddo
    endif
#endif
    !$acc end kernels

        !$acc parallel loop collapse(3) private(alpha, us2, un, unT)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    us2 = u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k)
                    do alpha=0,18
                        un(alpha) = u(i,j,k)*ex(alpha)+v(i,j,k)*ey(alpha)+w(i,j,k)*ez(alpha)
                        f(i,j,k,alpha) = rho(i,j,k)*omega(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)
                    enddo
                    do alpha=0,6
                        unT(alpha) = u(i,j,k)*ex(alpha)+v(i,j,k)*ey(alpha)+w(i,j,k)*ez(alpha)
                        g(i,j,k,alpha) = omegaT(alpha)*T(i,j,k)*(1.0d0+21.0d0/(6.0d0+paraA)*unT(alpha))
                    enddo
                enddo
            enddo
        enddo
    else
        write(*,*) "Error: initial field is not properly set"
    endif
    

    return
end subroutine initial


subroutine collision(i_start, i_end, j_start, j_end, k_start, k_end)
    !$acc routine gang nohost
    use commondata
    implicit none
    integer, intent(in) :: i_start, i_end, j_start, j_end, k_start, k_end
    integer :: i, j, k
    integer :: alpha
    real(kind=8) :: m(0:18), m_post(0:18), meq(0:18)
    real(kind=8) :: s(0:18)
    real(kind=8) :: fSource(0:18)

    !$acc loop independent collapse(3) &
    !$acc private(m, s, meq, m_post)
    do k = k_start, k_end
        do j = j_start, j_end
            do i = i_start, i_end
    !--------------------------------------------------------------------------------------------------------------------
    !---m0    
    m(0) =f(i,j,k,0)+f(i,j,k,1)+f(i,j,k,2)+f(i,j,k,3)+f(i,j,k,4)+f(i,j,k,5)+f(i,j,k,6) &
            +f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,11)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,18)
    !---m1
    m(1) = -30.0d0*f(i,j,k,0)-11.0d0*( f(i,j,k,1)+f(i,j,k,2)+f(i,j,k,3)+f(i,j,k,4)+f(i,j,k,5)+f(i,j,k,6) ) &
    +8.0d0*( f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,11)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,18) )
    !---m2
    m(2) = 12.0d0*f(i,j,k,0)-4.0d0*( f(i,j,k,1)+f(i,j,k,2)+f(i,j,k,3)+f(i,j,k,4)+f(i,j,k,5)+f(i,j,k,6) ) &
            +f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,11)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,18)
    !---m3
    m(3) = f(i,j,k,1)-f(i,j,k,2)+f(i,j,k,7)-f(i,j,k,8)+f(i,j,k,9)-f(i,j,k,10)+f(i,j,k,11)-f(i,j,k,12)+f(i,j,k,13)-f(i,j,k,14)
    !---m4
    m(4) = -4.0d0*(f(i,j,k,1)-f(i,j,k,2))+f(i,j,k,7)-f(i,j,k,8)+f(i,j,k,9)-f(i,j,k,10) &
                    +f(i,j,k,11)-f(i,j,k,12)+f(i,j,k,13)-f(i,j,k,14)
    !---m5
    m(5) = f(i,j,k,3)-f(i,j,k,4)+f(i,j,k,7)+f(i,j,k,8)-f(i,j,k,9)-f(i,j,k,10)+f(i,j,k,15)-f(i,j,k,16)+f(i,j,k,17)-f(i,j,k,18)
    !---m6
    m(6) = -4.0d0*(f(i,j,k,3)-f(i,j,k,4))+f(i,j,k,7)+f(i,j,k,8)-f(i,j,k,9)-f(i,j,k,10) &
                +f(i,j,k,15)-f(i,j,k,16)+f(i,j,k,17)-f(i,j,k,18)
    !---m7
    m(7) = f(i,j,k,5)-f(i,j,k,6)+f(i,j,k,11)+f(i,j,k,12)-f(i,j,k,13)-f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,16)-f(i,j,k,17)-f(i,j,k,18)
    !---m8
    m(8) = -4.0d0*(f(i,j,k,5)-f(i,j,k,6))+f(i,j,k,11)+f(i,j,k,12)-f(i,j,k,13)-f(i,j,k,14) &
            +f(i,j,k,15)+f(i,j,k,16)-f(i,j,k,17)-f(i,j,k,18)
    !---m9
    m(9) = 2.0d0*(f(i,j,k,1)+f(i,j,k,2))-f(i,j,k,3)-f(i,j,k,4)-f(i,j,k,5)-f(i,j,k,6) &
            +f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,11)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14)-2.0d0*( f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,18) )
    !---m10
     m(10) = -4.0d0*(f(i,j,k,1)+f(i,j,k,2))+2.0d0*(f(i,j,k,3)+f(i,j,k,4)+f(i,j,k,5)+f(i,j,k,6)) &
            +f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,11)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14)-2.0d0*( f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,18) )
    !---m11
    m(11) = f(i,j,k,3)+f(i,j,k,4)-f(i,j,k,5)-f(i,j,k,6)+f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)-( f(i,j,k,11)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14) )
    !---m12
     m(12) = -2.0d0*(f(i,j,k,3)+f(i,j,k,4)-f(i,j,k,5)-f(i,j,k,6))+( f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10) )-( f(i,j,k,11)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14) )
    !---m13
    m(13) = f(i,j,k,7)-f(i,j,k,8)-f(i,j,k,9)+f(i,j,k,10)
    !---m14
    m(14) = f(i,j,k,15)-f(i,j,k,16)-f(i,j,k,17)+f(i,j,k,18)
    !---m15
    m(15) = f(i,j,k,11)-f(i,j,k,12)-f(i,j,k,13)+f(i,j,k,14)
    !---m16
    m(16) = f(i,j,k,7)-f(i,j,k,8)+f(i,j,k,9)-f(i,j,k,10)-f(i,j,k,11)+f(i,j,k,12)-f(i,j,k,13)+f(i,j,k,14)
    !---m17
    m(17) = -f(i,j,k,7)-f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,15)-f(i,j,k,16)+f(i,j,k,17)-f(i,j,k,18)
    !---m18
    m(18) = f(i,j,k,11)+f(i,j,k,12)-f(i,j,k,13)-f(i,j,k,14)-f(i,j,k,15)-f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,18)
    !--------------------------------------------------------------------------------------------------------------------

    meq(0) = rho(i,j,k)
    meq(1) = -11.0d0*rho(i,j,k)+19.0d0*rho(i,j,k)*( u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k) )
    meq(2) = 3.0d0*rho(i,j,k)-11.0d0/2.0d0*rho(i,j,k)*( u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k) )
    meq(3) = rho(i,j,k)*u(i,j,k)
    meq(4) = -2.0d0/3.0d0*meq(3)
    meq(5) = rho(i,j,k)*v(i,j,k)
    meq(6) = -2.0d0/3.0d0*meq(5)
    meq(7) = rho(i,j,k)*w(i,j,k)
    meq(8) = -2.0d0/3.0d0*meq(7)
    meq(9) = rho(i,j,k)*(2.0d0*u(i,j,k)*u(i,j,k)-v(i,j,k)*v(i,j,k)-w(i,j,k)*w(i,j,k))
    meq(10) = -0.5d0*meq(9)
    meq(11) = rho(i,j,k)*(v(i,j,k)*v(i,j,k)-w(i,j,k)*w(i,j,k))
    meq(12) = -0.5d0*meq(11)
    meq(13) = rho(i,j,k)*(u(i,j,k)*v(i,j,k))
    meq(14) = rho(i,j,k)*(v(i,j,k)*w(i,j,k))
    meq(15) = rho(i,j,k)*(w(i,j,k)*u(i,j,k))
    meq(16) = 0.0d0
    meq(17) = 0.0d0
    meq(18) = 0.0d0

            s(0) = 0.0d0 
            s(1) = Snu  !!!s_{e}
            s(2) = Snu   !!! s_{epsilon}
            s(3) = 0.0d0 
            s(4) = Sq   !!! s_{q}
            s(5) = 0.0d0 
            s(6) = Sq   !!! s_{q}
            s(7) = 0.0d0 
            s(8) = Sq   !!! s_{q}
            s(9) = Snu !!! s_{nu}
            s(10) = Snu   !!! s_{pi}
            s(11) = Snu   !!! s_{nu}
            s(12) = Snu !!! s_{pi}
            s(13) = Snu !!! s_{nu}
            s(14) = Snu   !!! s_{nu}
            s(15) = Snu   !!! s_{nu}
            s(16) = Sq   !!! s_{m}
            s(17) = Sq   !!! s_{m}
            s(18) = Sq   !!! s_{m}

            Fx(i,j,k) = 0.0d0
            Fy(i,j,k) = 0.0d0
            Fz(i,j,k) = rho(i,j,k)*gBeta*(T(i,j,k)-Tref)

            fSource(0) = 0.0d0
            fSource(1) = 38.0d0*(u(i,j,k)*Fx(i,j,k)+v(i,j,k)*Fy(i,j,k)+w(i,j,k)*Fz(i,j,k))
            fSource(2) = -11.0d0*(u(i,j,k)*Fx(i,j,k)+v(i,j,k)*Fy(i,j,k)+w(i,j,k)*Fz(i,j,k))
            fSource(3) = Fx(i,j,k)
            fSource(4) = -2.0d0/3.0d0*Fx(i,j,k)
            fSource(5) = Fy(i,j,k)
            fSource(6) = -2.0d0/3.0d0*Fy(i,j,k)
            fSource(7) = Fz(i,j,k)
            fSource(8) = -2.0d0/3.0d0*Fz(i,j,k)
            fSource(9) = 4.0d0*u(i,j,k)*Fx(i,j,k)-2.0d0*v(i,j,k)*Fy(i,j,k)-2.0d0*w(i,j,k)*Fz(i,j,k)
            fSource(10) = -2.0d0*u(i,j,k)*Fx(i,j,k)+v(i,j,k)*Fy(i,j,k)+w(i,j,k)*Fz(i,j,k) 
            fSource(11) = 2.0d0*v(i,j,k)*Fy(i,j,k)-2.0d0*w(i,j,k)*Fz(i,j,k)
            fSource(12) = -v(i,j,k)*Fy(i,j,k)+w(i,j,k)*Fz(i,j,k)
            fSource(13) = u(i,j,k)*Fy(i,j,k)+v(i,j,k)*Fx(i,j,k)
            fSource(14) = v(i,j,k)*Fz(i,j,k)+w(i,j,k)*Fy(i,j,k)
            fSource(15) = u(i,j,k)*Fz(i,j,k)+w(i,j,k)*Fx(i,j,k)
            fSource(16) = 0.0d0
            fSource(17) = 0.0d0
            fSource(18) = 0.0d0

            !$acc loop seq
            do alpha=0,18
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+(1.0d0-0.5d0*s(alpha))*fSource(alpha)
            enddo

    f_post(i,j,k,0) = m_post(0)/19.0d0-5.0d0/399.0d0*m_post(1)+m_post(2)/21.0d0

    f_post(i,j,k,1) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
                +( m_post(3)-m_post(4) )*0.1d0+( m_post(9)-m_post(10) )/18.0d0

    f_post(i,j,k,2) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
                -( m_post(3)-m_post(4) )*0.1d0+( m_post(9)-m_post(10) )/18.0d0

    f_post(i,j,k,3) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
                +( m_post(5)-m_post(6) )*0.1d0-( m_post(9)-m_post(10) )/36.0d0 &
                +( m_post(11)-m_post(12) )/12.0d0

    f_post(i,j,k,4) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
                -( m_post(5)-m_post(6) )*0.1d0-( m_post(9)-m_post(10) )/36.0d0 &
                +( m_post(11)-m_post(12) )/12.0d0

    f_post(i,j,k,5) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
                +( m_post(7)-m_post(8) )*0.1d0-( m_post(9)-m_post(10) )/36.0d0 &
                -( m_post(11)-m_post(12) )/12.0d0

    f_post(i,j,k,6) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
                -( m_post(7)-m_post(8) )*0.1d0-( m_post(9)-m_post(10) )/36.0d0 &
                -( m_post(11)-m_post(12) )/12.0d0

    f_post(i,j,k,7) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
                +0.025d0*( 4.0d0*m_post(3)+m_post(4)+4.0d0*m_post(5)+m_post(6) ) &
                +m_post(9)/36.0d0+m_post(10)/72.0d0+m_post(11)/12.0d0+m_post(12)/24.0d0 &
                +m_post(13)*0.25d0+( m_post(16)-m_post(17) )*0.125d0

    f_post(i,j,k,8) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
                -0.025d0*( 4.0d0*m_post(3)+m_post(4)-4.0d0*m_post(5)-m_post(6) ) &
                +m_post(9)/36.0d0+m_post(10)/72.0d0+m_post(11)/12.0d0+m_post(12)/24.0d0 &
                -m_post(13)*0.25d0-( m_post(16)+m_post(17) )*0.125d0

    f_post(i,j,k,9) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
                +0.025d0*( 4.0d0*m_post(3)+m_post(4)-4.0d0*m_post(5)-m_post(6) ) &
                +m_post(9)/36.0d0+m_post(10)/72.0d0+m_post(11)/12.0d0+m_post(12)/24.0d0 &
                -m_post(13)*0.25d0+( m_post(16)+m_post(17) )*0.125d0

    f_post(i,j,k,10) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
                -0.025d0*( 4.0d0*m_post(3)+m_post(4)+4.0d0*m_post(5)+m_post(6) ) &
                +m_post(9)/36.0d0+m_post(10)/72.0d0+m_post(11)/12.0d0+m_post(12)/24.0d0 &
                +m_post(13)*0.25d0-( m_post(16)-m_post(17) )*0.125d0

    f_post(i,j,k,11) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
                +0.025d0*( 4.0d0*m_post(3)+m_post(4)+4.0d0*m_post(7)+m_post(8) ) &
                +m_post(9)/36.0d0+m_post(10)/72.0d0-m_post(11)/12.0d0-m_post(12)/24.0d0 &
                +0.25d0*m_post(15)-0.1250d0*( m_post(16)-m_post(18) )

    f_post(i,j,k,12) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
                -0.025d0*( 4.0d0*m_post(3)+m_post(4)-4.0d0*m_post(7)-m_post(8) ) &
                +m_post(9)/36.0d0+m_post(10)/72.0d0-m_post(11)/12.0d0-m_post(12)/24.0d0 &
                -0.25d0*m_post(15)+0.125d0*( m_post(16)+m_post(18) )

    f_post(i,j,k,13) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
                +0.025d0*( 4.0d0*m_post(3)+m_post(4)-4.0d0*m_post(7)-m_post(8) ) &
                +m_post(9)/36.0d0+m_post(10)/72.0d0-m_post(11)/12.0d0-m_post(12)/24.0d0 &
                -0.25d0*m_post(15)-0.125d0*( m_post(16)+m_post(18) )

    f_post(i,j,k,14) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
                -0.025d0*( 4.0d0*m_post(3)+m_post(4)+4.0d0*m_post(7)+m_post(8) ) &
                +m_post(9)/36.0d0+m_post(10)/72.0d0-m_post(11)/12.0d0-m_post(12)/24.0d0 &
                +0.25d0*m_post(15)+0.125d0*( m_post(16)-m_post(18) )

    f_post(i,j,k,15) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
                +( 4.0d0*m_post(5)+m_post(6)+4.0d0*m_post(7)+m_post(8) )*0.025d0 &
                -( m_post(9)+m_post(10)*0.5d0 )/18.0d0 &
                +0.25d0*m_post(14)+0.125d0*( m_post(17)-m_post(18) )

    f_post(i,j,k,16) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
                -( 4.0d0*m_post(5)+m_post(6)-4.0d0*m_post(7)-m_post(8) )*0.025d0 &
                -( m_post(9)+m_post(10)*0.5d0 )/18.0d0 &
                -0.25d0*m_post(14)-0.125d0*( m_post(17)+m_post(18) )

    f_post(i,j,k,17) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
                +( 4.0d0*m_post(5)+m_post(6)-4.0d0*m_post(7)-m_post(8) )*0.025d0 &
                -( m_post(9)+m_post(10)*0.5d0 )/18.0d0 &
                -0.25d0*m_post(14)+0.125d0*( m_post(17)+m_post(18) )

    f_post(i,j,k,18) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
                -( 4.0d0*m_post(5)+m_post(6)+4.0d0*m_post(7)+m_post(8) )*0.025d0 &
                -( m_post(9)+m_post(10)*0.5d0 )/18.0d0 &
                +0.25d0*m_post(14)-0.125d0*( m_post(17)-m_post(18) )

            enddo
        enddo
    enddo
    
    return
end subroutine collision

subroutine streaming()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: ip, jp, kp
    integer :: alpha
    
    !$acc parallel loop independent collapse(3) async(2)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                !$acc loop seq
                do alpha=0,18
                    ip = i-ex(alpha)
                    jp = j-ey(alpha)
                    kp = k-ez(alpha)
                    
                    f(i,j,k,alpha) = f_post(ip,jp,kp,alpha)
                
                enddo
            enddo
        enddo
    enddo
    
    return
end subroutine streaming

subroutine bounceback()
    use mpi_data
    use commondata
    implicit none
    integer :: i, j, k

    !$acc kernels async(2)
#ifdef noslipWalls
    !Back plane (i = 1)
    if (coords(0) == 0) then
        !$acc loop independent collapse(2) 
        do k=1,nz
            do j=1,ny
                f(1,j,k,1) = f_post(1,j,k,2)
                f(1,j,k,7) = f_post(1,j,k,10)
                f(1,j,k,9) = f_post(1,j,k,8)
                f(1,j,k,11) = f_post(1,j,k,14)
                f(1,j,k,13) = f_post(1,j,k,12)
            enddo
        enddo
    endif

    ! !Front plane (i=nx)
    if (coords(0) == dims(0) - 1) then
        !$acc loop independent collapse(2)
        do k=1,nz
            do j=1,ny
                f(nx,j,k,2) = f_post(nx,j,k,1)
                f(nx,j,k,8) = f_post(nx,j,k,9)
                f(nx,j,k,10) = f_post(nx,j,k,7)
                f(nx,j,k,12) = f_post(nx,j,k,13)
                f(nx,j,k,14) = f_post(nx,j,k,11)
            enddo
        enddo
    endif

    ! Left plane (j=1)
    if (coords(1) == 0) then
        !$acc loop independent collapse(2)
        do k=1,nz
            do i=1,nx
                f(i,1,k,3) = f_post(i,1,k,4)
                f(i,1,k,7) = f_post(i,1,k,10)
                f(i,1,k,8) = f_post(i,1,k,9)
                f(i,1,k,15) = f_post(i,1,k,18)
                f(i,1,k,17) = f_post(i,1,k,16)
            enddo
        enddo
    endif

    ! Right plane (j=ny)
    if  (coords(1) == dims(1) - 1) then
        !$acc loop independent collapse(2)
        do k=1,nz
            do i=1,nx
                !Right plane (j=ny)
                f(i,ny,k,4) = f_post(i,ny,k,3)
                f(i,ny,k,9) = f_post(i,ny,k,8)
                f(i,ny,k,10) = f_post(i,ny,k,7)
                f(i,ny,k,16) = f_post(i,ny,k,17)
                f(i,ny,k,18) = f_post(i,ny,k,15)
            enddo
        enddo
    endif
#endif

    ! Bottom side (k=1)
    if (coords(2) == 0) then
        !$acc loop independent collapse(2)
        do j=1,ny
            do i=1,nx
                f(i,j,1,5) = f_post(i,j,1,6)
                f(i,j,1,11) = f_post(i,j,1,14)
                f(i,j,1,12) = f_post(i,j,1,13)
                f(i,j,1,15) = f_post(i,j,1,18)
                f(i,j,1,16) = f_post(i,j,1,17)
            enddo
        enddo
    endif


    !Top side (k=nz)
    if (coords(2) == dims(2) - 1) then
        !$acc loop independent collapse(2)
        do j=1,ny
            do i=1,nx
                f(i,j,nz,6) = f_post(i,j,nz,5)
                f(i,j,nz,13) = f_post(i,j,nz,12)
                f(i,j,nz,14) = f_post(i,j,nz,11)
                f(i,j,nz,17) = f_post(i,j,nz,16)
                f(i,j,nz,18) = f_post(i,j,nz,15)
            enddo
        enddo
    endif
    !$acc end kernels

    return
end subroutine bounceback
    
subroutine macro()
    use commondata
    implicit none
    integer :: i, j, k

    !$acc parallel loop independent collapse(3) async(2)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                rho(i,j,k) = f(i,j,k,0)+f(i,j,k,1)+f(i,j,k,2)+f(i,j,k,3)+f(i,j,k,4)+f(i,j,k,5)+f(i,j,k,6) &
            +f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,11)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,18)
            
                u(i,j,k) = ( f(i,j,k,1)-f(i,j,k,2)+f(i,j,k,7)-f(i,j,k,8)+f(i,j,k,9)-f(i,j,k,10)+f(i,j,k,11)-f(i,j,k,12)+f(i,j,k,13)-f(i,j,k,14)+0.5d0*Fx(i,j,k) )/rho(i,j,k)
            
                v(i,j,k) = ( f(i,j,k,3)-f(i,j,k,4)+f(i,j,k,7)+f(i,j,k,8)-f(i,j,k,9)-f(i,j,k,10)+f(i,j,k,15)-f(i,j,k,16)+f(i,j,k,17)-f(i,j,k,18)+0.5d0*Fy(i,j,k) )/rho(i,j,k)
            
                w(i,j,k) = ( f(i,j,k,5)-f(i,j,k,6)+f(i,j,k,11)+f(i,j,k,12)-f(i,j,k,13)-f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,16)-f(i,j,k,17)-f(i,j,k,18)+0.5d0*Fz(i,j,k) )/rho(i,j,k)
            
            enddo
        enddo
    enddo

    return
end subroutine macro

subroutine collisionT(i_start, i_end, j_start, j_end, k_start, k_end)
    !$acc routine gang nohost
    use commondata
    implicit none
    integer, intent(in) :: i_start, i_end, j_start, j_end, k_start, k_end
    integer :: i, j, k
    integer :: alpha
    !------------------------
    real(kind=8) :: n(0:6), n_post(0:6), neq(0:6)
    real(kind=8) :: q(0:6)

    !$acc loop independent collapse(3) &
    !$acc private(n, q, neq, n_post)
    do k = k_start, k_end
        do j = j_start, j_end
            do i = i_start, i_end
            
    n(0) = g(i,j,k,0)+g(i,j,k,1)+g(i,j,k,2)+g(i,j,k,3)+g(i,j,k,4)+g(i,j,k,5)+g(i,j,k,6)
    n(1) = g(i,j,k,1)-g(i,j,k,2)
    n(2) = g(i,j,k,3)-g(i,j,k,4)
    n(3) = g(i,j,k,5)-g(i,j,k,6)
    n(4) = -6.0d0*g(i,j,k,0)+g(i,j,k,1)+g(i,j,k,2)+g(i,j,k,3)+g(i,j,k,4)+g(i,j,k,5)+g(i,j,k,6)
    n(5) = 2.0d0*g(i,j,k,1)+2.0d0*g(i,j,k,2)-g(i,j,k,3)-g(i,j,k,4)-g(i,j,k,5)-g(i,j,k,6)
    n(6) = g(i,j,k,3)+g(i,j,k,4)-g(i,j,k,5)-g(i,j,k,6)
        
            neq(0) = T(i,j,k)
            neq(1) = T(i,j,k)*u(i,j,k)
            neq(2) = T(i,j,k)*v(i,j,k)
            neq(3) = T(i,j,k)*w(i,j,k)
            neq(4) = T(i,j,k)*paraA
            neq(5) = 0.0d0
            neq(6) = 0.0d0
        
            q(0) = 0.0d0
            q(1) = Qd
            q(2) = Qd
            q(3) = Qd
            q(4) = Qnu
            q(5) = Qnu
            q(6) = Qnu
        
            !$acc loop seq
            do alpha=0,6
                n_post(alpha) = n(alpha)-q(alpha)*(n(alpha)-neq(alpha))
            enddo
        
    g_post(i,j,k,0) = n_post(0)/7.0d0-n_post(4)/7.0d0
    g_post(i,j,k,1) = n_post(0)/7.0d0+0.5d0*n_post(1)+n_post(4)/42.0d0+n_post(5)/6.0d0
    g_post(i,j,k,2) = n_post(0)/7.0d0-0.5d0*n_post(1)+n_post(4)/42.0d0+n_post(5)/6.0d0
    g_post(i,j,k,3) = n_post(0)/7.0d0+0.5d0*n_post(2)+n_post(4)/42.0d0-n_post(5)/12.0d0+0.25d0*n_post(6)
    g_post(i,j,k,4) = n_post(0)/7.0d0-0.5d0*n_post(2)+n_post(4)/42.0d0-n_post(5)/12.0d0+0.25d0*n_post(6)
    g_post(i,j,k,5) = n_post(0)/7.0d0+0.5d0*n_post(3)+n_post(4)/42.0d0-n_post(5)/12.0d0-0.25d0*n_post(6)
    g_post(i,j,k,6) = n_post(0)/7.0d0-0.5d0*n_post(3)+n_post(4)/42.0d0-n_post(5)/12.0d0-0.25d0*n_post(6)
        
            enddo
        enddo
    enddo
    
    return
end subroutine collisionT

subroutine streamingT()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: ip, jp, kp
    integer :: alpha
    
    !$acc parallel loop independent collapse(3) async(1)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                !$acc loop seq
                do alpha=0,6
                    ip = i-ex(alpha)
                    jp = j-ey(alpha)
                    kp = k-ez(alpha)
                    
                    g(i,j,k,alpha) = g_post(ip,jp,kp,alpha)
                    
                enddo
            enddo
        enddo
    enddo

    return
end subroutine streamingT

subroutine bouncebackT()
    use mpi_data
    use commondata
    implicit none
    integer :: i, j, k

    !$acc kernels async(1)
#ifdef TopBottomPlatesAdiabatic
    ! Bottom side (k=1)
    if (coords(2) == 0) then
        !$acc loop independent collapse(2)
        do j=1,ny
            do i=1,nx
                g(i,j,1,5) = g_post(i,j,1,6)
            enddo
        enddo
    endif

    ! Top side (k=nz)
    if (coords(2) == dims(2) - 1) then
        !$acc loop independent collapse(2)
        do j=1,ny
            do i=1,nx
                g(i,j,nz,6) = g_post(i,j,nz,5)
            enddo
        enddo
    endif
#endif
#ifdef TopBottomPlatesConstT
    ! Bottom side (k=1)
    if (coords(2) == 0) then
        !$acc loop independent collapse(2)
        do j=1,ny
            do i=1,nx
                g(i,j,1,5) = -g_post(i,j,1,6)+(6.0d0+paraA)/21.0d0*Thot
            enddo
        enddo
    endif

    ! Top side (k=nz)
    if (coords(2) == dims(2) - 1) then
        !$acc loop independent collapse(2)
        do j=1,ny
            do i=1,nx
                g(i,j,nz,6) = -g_post(i,j,nz,5)+(6.0d0+paraA)/21.0d0*Tcold
            enddo
        enddo
    endif
#endif
#ifdef LeftRightWallsConstT
    ! Left side (j=1)
    if (coords(1) == 0) then
        !$acc loop independent collapse(2)
        do k=1,nz
            do i=1,nx
                g(i,1,k,3) = -g_post(i,1,k,4)+(6.0d0+paraA)/21.0d0*Thot
            enddo
        enddo
    endif

    ! Right side (j=ny)
    if (coords(1) == dims(1) - 1) then
        !$acc loop independent collapse(2)
        do k=1,nz
            do i=1,nx
                g(i,ny,k,4) = -g_post(i,ny,k,3)+(6.0d0+paraA)/21.0d0*Tcold
            enddo
        enddo
    endif
#endif
#ifdef LeftRightWallsAdiabatic
    ! Left side (j=1)
    if (coords(1) == 0) then
        !$acc loop independent collapse(2)
        do k=1,nz
            do i=1,nx
                g(i,1,k,3) = g_post(i,1,k,4)
            enddo
        enddo
    endif

    ! Right side (j=ny)
    if (coords(1) == dims(1) - 1) then
        !$acc loop independent collapse(2)
        do k=1,nz
            do i=1,nx
                g(i,ny,k,4) = g_post(i,ny,k,3)
            enddo
        enddo
    endif
#endif   
#ifdef BackFrontWallsAdiabatic
    ! Back side (i=1)
    if (coords(0) == 0) then
        !$acc loop independent collapse(2)
        do k=1,nz
            do j=1,ny
                g(1,j,k,1) = g_post(1,j,k,2)
            enddo
        enddo
    endif

    !Front side (i=nx)
    if (coords(0) == dims(0) - 1) then
        !$acc loop independent collapse(2)
        do k=1,nz
            do j=1,ny          
                g(nx,j,k,2) = g_post(nx,j,k,1)
            enddo
        enddo
    endif
#endif
    !$acc end kernels

    return
end subroutine bouncebackT

subroutine macroT()
    use commondata
    implicit none
    integer :: i, j, k

    !$acc parallel loop independent collapse(3) async(1)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                T(i,j,k) = g(i,j,k,0)+g(i,j,k,1)+g(i,j,k,2)+g(i,j,k,3)+g(i,j,k,4)+g(i,j,k,5)+g(i,j,k,6)
            enddo
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
    call collision(1, nx, 1, ny, 1, 1)
    !$acc end kernels

    !$acc kernels async(4)
    call collision(1, nx, 1, ny, nz, nz)
    !$acc end kernels

    !$acc kernels async(5)
    call collision(1, nx, 1, 1, 1, nz)
    !$acc end kernels

    !$acc kernels async(6)
    call collision(1, nx, ny, ny, 1, nz)
    !$acc end kernels

end subroutine compute_boundary_f_yz


subroutine flow_update()
    use commondata
    use mpi_data
    use mpi
    implicit none

    ! async(3, 4, 5, 6)
    call compute_boundary_f_yz()
    !$acc wait

    ! async(3)
    call pack_f_z()
        !$acc kernels async(1)
        call collisionT(1, nx, 1, ny, 1, nz)
        !$acc end kernels
        !$acc wait(3)
        call mpi_f_z()
    ! async(3)
    call unpack_f_z()
    !$acc wait

    ! async(3,4,5)
    call pack_f_y_g()
        !$acc kernels async(2)
        call collision(1, nx, 2, ny-1, 2, nz-1)
        !$acc end kernels
        !$acc wait(3,4,5)
        call mpi_f_y_g()
    ! async(3,4,5)
    call unpack_f_y_g()
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
    ! async(2)
    call bounceback()
    ! async(2)
    call macro()

    !$acc wait


end subroutine


subroutine mpi_f_y_g()
    use mpi
    use mpi_data
    use commondata
    implicit none

    call mpi_f_y()

    call mpi_g_y()
    call mpi_g_x()
    call mpi_g_z()

contains
    subroutine mpi_f_y()
        use mpi
        use mpi_data
        use commondata
        implicit none

        ! ------------ exchange message along y ----------------
        ! message passing to (j++)
        !$acc host_data use_device(send_pos, recv_pos, send_neg, recv_neg)
        call MPI_Sendrecv(send_pos, 5 * surface_length_y, MPI_DOUBLE_PRECISION, nbr_surface(3), 44, &
            recv_pos, 5 * surface_length_y, MPI_DOUBLE_PRECISION, nbr_surface(4), 44, &
            comm3d, MPI_STATUS_IGNORE, rc)

        ! message passing to (j--)
        call MPI_Sendrecv(send_neg, 5*surface_length_y, MPI_DOUBLE_PRECISION, nbr_surface(4), 55, &
            recv_neg, 5*surface_length_y, MPI_DOUBLE_PRECISION, nbr_surface(3), 55, &
            comm3d, MPI_STATUS_IGNORE, rc)
        !$acc end host_data

    end subroutine

    subroutine mpi_g_x()
        use mpi
        use mpi_data
        use commondata
        implicit none
    
        ! ------------ exchange message along x ----------------
        ! message passing to (i++)
        !$acc host_data use_device(g_send_pos_x, g_recv_pos_x, g_send_neg_x, g_recv_neg_x)
        call MPI_Sendrecv(g_send_pos_x, surface_length_x, MPI_DOUBLE_PRECISION, nbr_surface(1), 1, &
            g_recv_pos_x, surface_length_x, MPI_DOUBLE_PRECISION, nbr_surface(2), 1, &
            comm3d, MPI_STATUS_IGNORE, rc)

        ! message passing to (i--)
        call MPI_Sendrecv(g_send_neg_x, surface_length_x, MPI_DOUBLE_PRECISION, nbr_surface(2), 2, &
            g_recv_neg_x, surface_length_x, MPI_DOUBLE_PRECISION, nbr_surface(1), 2, &
            comm3d, MPI_STATUS_IGNORE, rc)
        !$acc end host_data
    
    end subroutine

    subroutine mpi_g_y()
        use mpi
        use mpi_data
        use commondata
        implicit none
    
        ! ------------ exchange message along y ----------------
        ! message passing to (j++)
        !$acc host_data use_device(g_send_pos_y, g_recv_pos_y, g_send_neg_y, g_recv_neg_y)
        call MPI_Sendrecv(g_send_pos_y, surface_length_y, MPI_DOUBLE_PRECISION, nbr_surface(3), 3, &
            g_recv_pos_y, surface_length_y, MPI_DOUBLE_PRECISION, nbr_surface(4), 3, &
            comm3d, MPI_STATUS_IGNORE, rc)

        ! message passing to (j--)
        call MPI_Sendrecv(g_send_neg_y, surface_length_y, MPI_DOUBLE_PRECISION, nbr_surface(4), 4, &
            g_recv_neg_y, surface_length_y, MPI_DOUBLE_PRECISION, nbr_surface(3), 4, &
            comm3d, MPI_STATUS_IGNORE, rc)
        !$acc end host_data
    
    end subroutine

    subroutine mpi_g_z()
        use mpi
        use mpi_data
        use commondata
        implicit none

        ! ! ------------ exchange message along z ----------------
        ! ! message passing to (k++)
        ! !$acc host_data use_device(g_post)  
        ! call MPI_Sendrecv(g_send_pos_z, surface_length_z, MPI_DOUBLE_PRECISION, nbr_surface(5), 5, &
        !     g_recv_pos_z, surface_length_z, MPI_DOUBLE_PRECISION, nbr_surface(6), 5, &
        !     comm3d, MPI_STATUS_IGNORE, rc)
    
        ! ! message passing to (k--)
        ! call MPI_Sendrecv(g_send_neg_z, surface_length_z, MPI_DOUBLE_PRECISION, nbr_surface(6), 6, &
        !     g_recv_neg_z, surface_length_z, MPI_DOUBLE_PRECISION, nbr_surface(5), 6, &
        !     comm3d, MPI_STATUS_IGNORE, rc)
        ! !$acc end host_data

        
        !$acc host_data use_device(g_post)  
        ! message passing to (k++)
        call MPI_Sendrecv(g_post(0, 0, nz, 5), surface_length_z, MPI_DOUBLE_PRECISION, nbr_surface(5), 5, &
            g_post(0, 0, 0, 5), surface_length_z, MPI_DOUBLE_PRECISION, nbr_surface(6), 5, &
            comm3d, MPI_STATUS_IGNORE, rc)     

        ! message passing to (k--)
        call MPI_Sendrecv(g_post(0, 0, 1, 6), surface_length_z, MPI_DOUBLE_PRECISION, nbr_surface(6), 6, &
            g_post(0, 0, nz+1, 6), surface_length_z, MPI_DOUBLE_PRECISION, nbr_surface(5), 6, &
            comm3d, MPI_STATUS_IGNORE, rc)     
        !$acc end host_data



    end subroutine

end subroutine

subroutine pack_f_y_g()
    use mpi_data
    use commondata

    ! async(3)
    call pack_f_y()
    ! async(4,5)
    call pack_g()

contains
    subroutine pack_f_y()
        use mpi_data
        use commondata
        implicit none
        integer :: i, j, k, idx, tmp

        ! pack the message to send
        !$acc parallel loop independent collapse(3) private(tmp) async(3)
        do idx = 0, 4
            do k = 0, nz+1
                do i = 0, nx+1
                    tmp = i + (nx+2) * k + surface_length_y * idx + 1
                    send_pos(tmp) = f_post(i, ny, k, f_tag_y_pos(idx))
                    send_neg(tmp) = f_post(i, 1, k, f_tag_y_neg(idx))
                enddo
            enddo
        enddo

    end subroutine

    subroutine pack_g()
        use mpi_data
        use commondata
        implicit none
        integer :: i, j, k, tmp
            
        ! pack the message to send along x
        !$acc parallel loop independent collapse(2) private(tmp) async(4)
        do k = 0, nz+1
            do j = 0, ny+1
                tmp = j + (ny+2) * k + 1
                g_send_pos_x(tmp) = g_post(nx, j, k, 1)
                g_send_neg_x(tmp) = g_post(1, j, k, 2)
            enddo
        enddo

        ! pack the message to send along y
        !$acc parallel loop independent collapse(2) private(tmp) async(5)
        do k = 0, nz+1
            do i = 0, nx+1
                tmp = i + (nx+2) * k + 1
                g_send_pos_y(tmp) = g_post(i, ny, k, 3)
                g_send_neg_y(tmp) = g_post(i, 1, k, 4)
            enddo
        enddo

        ! ! pack the message to send along z
        ! !$acc parallel loop independent collapse(2) private(tmp)
        ! do j = 0, ny+1
        !     do i = 0, nx+1
        !         tmp = i + (nx+2) * j + 1
        !         g_send_pos_z(tmp) = g_post(i, j, nz, 5)
        !         g_send_neg_z(tmp) = g_post(i, j, 1, 6)
        !     enddo
        ! enddo
    
    end subroutine

end subroutine

subroutine unpack_f_y_g()
    use mpi_data
    use commondata

    ! async(3)
    call unpack_f_y()
    ! async(4,5)
    call unpack_g()

contains

    subroutine unpack_f_y()
        use mpi_data
        use commondata
        implicit none
        integer :: i, j, k, idx, tmp

        ! depack the recvied message
        !$acc parallel loop independent collapse(3) private(tmp) async(3)
        do idx = 0, 4
            do k = 0, nz+1
                do i = 0, nx+1
                    tmp = i + (nx+2) * k + surface_length_y * idx + 1
                    f_post(i, 0, k, f_tag_y_pos(idx)) = recv_pos(tmp)
                    f_post(i, ny+1, k, f_tag_y_neg(idx)) = recv_neg(tmp)
                enddo
            enddo
        enddo

    end subroutine

    subroutine unpack_g()
        use mpi_data
        use commondata
        implicit none
        integer :: i, j, k, tmp
    
        ! unpack the recvied message along x
        !$acc parallel loop independent collapse(2) private(tmp) async(4)
        do k = 0, nz+1
            do j = 0, ny+1
                tmp = j + (ny+2) * k + 1
                g_post(0, j, k, 1) = g_recv_pos_x(tmp) 
                g_post(nx+1, j, k, 2) = g_recv_neg_x(tmp) 
            enddo
        enddo

        ! unpack the recvied message along y
        !$acc parallel loop independent collapse(2) private(tmp) async(5)
        do k = 0, nz+1
            do i = 0, nx+1
                tmp = i + (nx+2) * k + 1
                g_post(i, ny+1, k, 4) = g_recv_neg_y(tmp)
                g_post(i, 0, k, 3) = g_recv_pos_y(tmp)
            enddo
        enddo

        ! ! unpack the recvied message along z
        ! !$acc parallel loop independent collapse(2) private(tmp)
        ! do j = 0, ny+1
        !     do i = 0, nx+1
        !         tmp = i + (nx+2) * j + 1
        !         g_post(i, j, 0, 5) = g_recv_pos_z(tmp) 
        !         g_post(i, j, nz+1, 6) = g_recv_neg_z(tmp) 
        !     enddo
        ! enddo
    
    end subroutine

end subroutine

subroutine check()
    use mpi
    use commondata
    use mpi_data
    implicit none
    integer :: i, j, k
    real(kind=8) :: error1, error2, error5, error6
    real(8) :: total_error1, total_error2, total_error5, total_error6

    error1 = 0.0d0
    error2 = 0.0d0

    error5 = 0.0d0
    error6 = 0.0d0
    
    !$acc kernels

    !$acc loop independent collapse(3) &
    !$acc reduction(+:error1, error2, error5, error6)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                error1 = error1+(u(i,j,k)-up(i,j,k))*(u(i,j,k)-up(i,j,k))+(v(i,j,k)-vp(i,j,k))*(v(i,j,k)-vp(i,j,k))+(w(i,j,k)-wp(i,j,k))*(w(i,j,k)-wp(i,j,k))
                error2 = error2+u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k)
                
                error5 = error5+dABS( T(i,j,k)-Tp(i,j,k) )
                error6 = error6+dABS( T(i,j,k) )
                
                up(i,j,k) = u(i,j,k)
                vp(i,j,k) = v(i,j,k)
                wp(i,j,k) = w(i,j,k)
                Tp(i,j,k) = T(i,j,k)
            enddo
        enddo
    enddo

    !$acc end kernels

    call MPI_Barrier(comm3d, rc)

    call MPI_ALLreduce(error1, total_error1, 1, MPI_REAL8, MPI_SUM, comm3d, rc)
    call MPI_ALLreduce(error2, total_error2, 1, MPI_REAL8, MPI_SUM, comm3d, rc)
    call MPI_ALLreduce(error5, total_error5, 1, MPI_REAL8, MPI_SUM, comm3d, rc)
    call MPI_ALLreduce(error6, total_error6, 1, MPI_REAL8, MPI_SUM, comm3d, rc)
    
    errorU = dsqrt(total_error1)/dsqrt(total_error2)
    errorT = total_error5 / total_error6

    if (rank == 0) then
        write(*,*) itc,' ',errorU,' ',errorT
    endif

    return
end subroutine check

subroutine pack_f_z()
    use mpi_data
    use commondata
    implicit none
    integer :: i, j, k, idx, tmp

    ! pack the message to send
    !$acc parallel loop independent collapse(3) private(tmp) async(3)
    do idx = 0, 4
        do j = 0, ny+1
            do i = 0, nx+1
                tmp = i + (nx+2) * j + surface_length_z * idx + 1
                send_pos(tmp) = f_post(i, j, nz, f_tag_z_pos(idx))
                send_neg(tmp) = f_post(i, j, 1, f_tag_z_neg(idx))
            enddo
        enddo
    enddo
end subroutine

subroutine mpi_f_z()
    use mpi
    use mpi_data
    use commondata
    implicit none
    integer :: idx, req(4)

    ! ------------ exchange message along z ----------------
    ! message passing to (k++)
    !$acc host_data use_device(send_pos, recv_pos, send_neg, recv_neg)
    call MPI_Sendrecv(send_pos, 5*surface_length_z, MPI_DOUBLE_PRECISION, nbr_surface(5), 55, &
        recv_pos, 5*surface_length_z, MPI_DOUBLE_PRECISION, nbr_surface(6), 55, &
        comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (k--)
    call MPI_Sendrecv(send_neg, 5*surface_length_z, MPI_DOUBLE_PRECISION, nbr_surface(6), 66, &
        recv_neg, 5*surface_length_z, MPI_DOUBLE_PRECISION, nbr_surface(5), 66, &
        comm3d, MPI_STATUS_IGNORE, rc)
    !$acc end host_data

    call MPI_Waitall(4, req, MPI_STATUSES_IGNORE, rc)

end subroutine

subroutine unpack_f_z()
    use mpi_data
    use commondata
    implicit none
    integer :: i, j, k, idx, tmp
    ! unpack the recvied message
    !$acc parallel loop independent collapse(3) private(tmp) async(3)
    do idx = 0, 4
        do j = 0, ny+1
            do i = 0, nx+1
                tmp = i + (nx+2) * j + surface_length_z * idx + 1
                f_post(i, j, 0, f_tag_z_pos(idx)) = recv_pos(tmp) 
                f_post(i, j, nz+1, f_tag_z_neg(idx)) = recv_neg(tmp) 
            enddo
        enddo
    enddo
endsubroutine

subroutine pack_f_x()
    use mpi_data
    use commondata
    implicit none
    integer :: i, j, k, idx, tmp

    ! pack the message to send
    !$acc parallel loop independent collapse(3) private(tmp) async(3)
    do idx = 0, 4
        do k = 0, nz+1
            do j = 0, ny+1
                tmp = j + (ny+2) * k + surface_length_x * idx + 1
                send_pos(tmp) = f_post(nx, j, k, f_tag_x_pos(idx))
                send_neg(tmp) = f_post(1, j, k, f_tag_x_neg(idx))
            enddo
        enddo
    enddo

end subroutine

subroutine mpi_f_x()
    use mpi
    use mpi_data
    use commondata
    implicit none
    integer :: req(4)

    ! ------------ exchange message along x ----------------
    ! message passing to (i++)
    !$acc host_data use_device(send_pos, recv_pos, send_neg, recv_neg)
    call MPI_Sendrecv(send_pos, 5*surface_length_x, MPI_DOUBLE_PRECISION, nbr_surface(1), 22, &
        recv_pos, 5*surface_length_x, MPI_DOUBLE_PRECISION, nbr_surface(2), 22, &
        comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i--)
    call MPI_Sendrecv(send_neg, 5*surface_length_x, MPI_DOUBLE_PRECISION, nbr_surface(2), 33, &
        recv_neg, 5*surface_length_x, MPI_DOUBLE_PRECISION, nbr_surface(1), 33, &
        comm3d, MPI_STATUS_IGNORE, rc)
    !$acc end host_data

    call MPI_Waitall(4, req, MPI_STATUSES_IGNORE, rc)

end subroutine

subroutine unpack_f_x()
    use mpi_data
    use commondata
    implicit none
    integer :: i, j, k, idx, tmp
    !$acc parallel loop independent collapse(3) private(tmp) async(2)
    do idx = 0, 4
        do k = 0, nz+1
            do j = 0, ny+1
                tmp = j + (ny+2) * k + surface_length_x * idx + 1
                f_post(0, j, k, f_tag_x_pos(idx)) = recv_pos(tmp) 
                f_post(nx+1, j, k, f_tag_x_neg(idx)) = recv_neg(tmp) 
            enddo
        enddo
    enddo

end subroutine


#ifdef __OUT_PUT__
    subroutine output()
        use mpi
        use commondata
        implicit none
        integer :: i, j, k
        integer :: p_rank, num(0:5), new_coords(0:2), dx = 0, dy = 0, dz = 0
        real(8), allocatable :: total_u(:, :, :), total_v(:, :, :), total_w(:, :, :), total_rho(:, :, :), total_T(:, :, :)
        real(8), allocatable :: tmp_u(:, :, :), tmp_v(:, :, :), tmp_w(:, :, :), tmp_rho(:, :, :), tmp_T(:, :, :)
        
        ! use rank 0 to receive data and output the results

        if (rank3d > 0) then  !!! ----  rank != 0 send data
            ! collect the rank information
            num(0) = nx
            num(1) = ny
            num(2) = nz
            num(3) = i_start_global
            num(4) = j_start_global
            num(5) = k_start_global
            ! send to rank 0
            call MPI_Send(num, 6, MPI_INTEGER, 0, 0, comm3d, rc)    ! rank information
            call MPI_Send(u, nx*ny*nz, MPI_REAL8, 0, 1, comm3d, rc)
            call MPI_Send(v, nx*ny*nz, MPI_REAL8, 0, 2, comm3d, rc)
            call MPI_Send(w, nx*ny*nz, MPI_REAL8, 0, 3, comm3d, rc)
            call MPI_Send(rho, nx*ny*nz, MPI_REAL8, 0, 4, comm3d, rc)
            call MPI_Send(T, nx*ny*nz, MPI_REAL8, 0, 5, comm3d, rc)
        else    
            !!! ---- rank 0 collect data
            ! allocate array
            allocate(total_u(total_nx, total_ny, total_nz))
            allocate(total_v(total_nx, total_ny, total_nz))
            allocate(total_w(total_nx, total_ny, total_nz))
            allocate(total_rho(total_nx, total_ny, total_nz))
            allocate(total_T(total_nx, total_ny, total_nz))

            ! determine the origin
            dx = i_start_global
            dy = j_start_global
            dz = k_start_global

            ! collect data from rank 0
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        total_u(dx + i, dy + j, dz + k) = u(i, j, k)
                        total_v(dx + i, dy + j, dz + k) = v(i, j, k)
                        total_w(dx + i, dy + j, dz + k) = w(i, j, k)
                        total_rho(dx + i, dy + j, dz + k) = rho(i, j, k)
                        total_T(dx + i, dy + j, dz + k) = T(i, j, k)
                    enddo
                enddo
            enddo

            ! collect data from all other processors
            do p_rank = 1, dims(0) * dims(1) * dims(2) - 1

                call MPI_Cart_coords(comm3d, p_rank, 3, new_coords, rc)

                ! receive the block size and origion
                call MPI_Recv(num, 6, MPI_INTEGER, p_rank, 0, comm3d, MPI_STATUS_IGNORE, rc)

                ! creat buffer
                allocate(tmp_u(num(0), num(1), num(2)))
                allocate(tmp_v(num(0), num(1), num(2)))
                allocate(tmp_w(num(0), num(1), num(2)))
                allocate(tmp_rho(num(0), num(1), num(2)))
                allocate(tmp_T(num(0), num(1), num(2)))

                ! receive data
                call MPI_Recv(tmp_u, num(0) * num(1) * num(2), MPI_REAL8, p_rank, 1, comm3d, MPI_STATUS_IGNORE, rc)
                call MPI_Recv(tmp_v, num(0) * num(1) * num(2), MPI_REAL8, p_rank, 2, comm3d, MPI_STATUS_IGNORE, rc)
                call MPI_Recv(tmp_w, num(0) * num(1) * num(2), MPI_REAL8, p_rank, 3, comm3d, MPI_STATUS_IGNORE, rc)
                call MPI_Recv(tmp_rho, num(0) * num(1) * num(2), MPI_REAL8, p_rank, 4, comm3d, MPI_STATUS_IGNORE, rc)
                call MPI_Recv(tmp_T, num(0) * num(1) * num(2), MPI_REAL8, p_rank, 5, comm3d, MPI_STATUS_IGNORE, rc)


                ! determine the origin
                dx = num(3)
                dy = num(4)
                dz = num(5)

                ! assign data
                do k = 1, num(2)
                    do j = 1, num(1)
                        do i = 1, num(0)
                            total_u(i + dx, j + dy, k + dz) = tmp_u(i, j, k)
                            total_v(i + dx, j + dy, k + dz) = tmp_v(i, j, k)
                            total_w(i + dx, j + dy, k + dz) = tmp_w(i, j, k)
                            total_rho(i + dx, j + dy, k + dz) = tmp_rho(i, j, k)
                            total_T(i + dx, j + dy, k + dz) = tmp_T(i, j, k)
                        enddo
                    enddo
                enddo

                ! de-allocate buffer arrays
                deallocate(tmp_u)
                deallocate(tmp_v)
                deallocate(tmp_w)
                deallocate(tmp_rho)
                deallocate(tmp_T)
            enddo

            ! after collect total_* data, then output
            ! call output_ASCII(xp, yp, zp, total_u, total_v, total_w, total_rho, total_T, total_nx, total_ny, total_nz, itc)
            call output_binary(total_u, total_v, total_w, total_rho, total_T, total_nx, total_ny, total_nz, itc)
            call output_Tecplot(xp, yp, zp, total_u, total_v, total_w, total_rho, total_T, total_nx, total_ny, total_nz, itc)
            ! call out_Velocity_Nu(total_u, total_v, total_w, total_T, total_nx, total_ny, total_nz, diffusivity, lengthUnit)

            ! de-allocate total arrays
            deallocate(total_u)
            deallocate(total_v)
            deallocate(total_w)
            deallocate(total_rho)
            deallocate(total_T)
        endif

    end subroutine output



    subroutine output_binary(u, v, w, rho, T, nx, ny, nz, itc)
        implicit none
        integer, intent(in) :: nx, ny, nz, itc
        real(8), intent(in) :: u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz), rho(nx, ny, nz), T(nx, ny, nz)
        integer :: i, j, k
        character(len=100) :: filename


    #ifdef benchmarkCavity
        write(filename,*) itc
        filename = adjustl(filename)
    #endif

        open(unit=01,file='buoyancyCavity-'//trim(filename)//'.bin',form="unformatted",access="sequential")
        write(01) (((u(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        write(01) (((v(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        write(01) (((w(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        write(01) (((T(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        close(01)

        return
    end subroutine output_binary




    subroutine output_Tecplot(xp, yp, zp, u, v, w, rho, T, nx, ny, nz, itc)
        implicit none
        integer, intent(in) :: nx, ny, nz, itc
        real(8), intent(in) :: xp(0:nx+1), yp(0:ny+1), zp(0:nz+1)
        real(8), intent(in) :: u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz), rho(nx, ny, nz), T(nx, ny, nz)
        integer :: i, j, k
        character(len=9) :: B2
        REAL(4) :: zoneMarker, eohMarker
        character(len=40) :: title
        character(len=40) :: V1,V2,V3,V4,V5,V6,V7
        character(len=40) :: zoneName
        character(len=100) :: filename
    

    #ifdef benchmarkCavity
        write(filename,*) itc
        filename = adjustl(filename)
    #endif
        open(41,file='buoyancyCavity-'//trim(filename)//'.plt', access='stream', form='unformatted')

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

        !c-- Number of variables in this data file (here 5 variables)
        write(41) 7

        !c-- Variable names.
        V1='X'
        call dumpstring(V1)
        V2='Y'
        call dumpstring(V2)
        V3='Z'
        call dumpstring(V3)
        
        V4='U'
        call dumpstring(V4)
        V5='V'
        call dumpstring(V5)
        V6='W'
        call dumpstring(V6)
        
        V7='T'
        call dumpstring(V7)

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
        write(41) nz

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
        write(41) 1

        !--------Has variable sharing 0 = no, 1 = yes.
        write(41) 0

        !----------Zone number to share connectivity list with (-1 = no
        ! sharing).
        write(41) -1

        !---------------------------------------------------------------------
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    write(41) real(xp(i))
                    write(41) real(yp(j))
                    write(41) real(zp(k))
                    write(41) real(u(i,j,k))
                    write(41) real(v(i,j,k))
                    write(41) real(w(i,j,k))
                    write(41) real(T(i,j,k))
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
        integer :: stringLength
        integer :: ii
        integer :: I

        stringLength=LEN_TRIM(instring)
        do ii=1,stringLength
            I=ICHAR(instring(ii:ii))
            write(41) I
        end do
        write(41) 0

        return
    end subroutine dumpstring

#endif
