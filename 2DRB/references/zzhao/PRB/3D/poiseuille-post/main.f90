module commondata
        implicit none
        ! Grid parameters
        integer(kind=4), parameter :: nz=41, nx=41, ny=nz*2+1
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1

        real(kind=8), parameter :: Reynolds=18.3d0
        real(kind=8), parameter :: u_max=0.1d0
        real(kind=8) :: viscosity_LB
        real(kind=8) :: length_LB
        real(kind=8) :: lengthScale
        real(kind=8) :: H
        real(kind=8) :: viscosity

        ! Grid coordinates
        real(kind=8) :: xGrid(nx), yGrid(ny), zGrid(0:nz+1)

        ! Flow variables
        real(kind=8), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)      ! Velocity components
        real(kind=8), allocatable ::  vana(:)
        real(kind=8) :: dt0, dx0
        integer(kind=4), parameter :: n=5
        real(kind=8) :: time(5)=(/0.1546, 0.3091, 0.6182, 1.2365, 2.4729/)

end module commondata

program main
    use commondata
    implicit none

    call initial()
    call loadData()
    !call calV()
    call solution()
    call calErr()

    deallocate (u)
    deallocate (v)
    deallocate (w)
    deallocate (vana)
    return
end program main

subroutine initial()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k
    real(kind=8) :: constA
    real(kind=8) :: dz(nz+1)

    allocate (u(nx,ny,nz))
    allocate (v(nx,ny,nz))
    allocate (w(nx,ny,nz))
    allocate (vana(nz))

    constA = 3.2d0
    do i = 1, nx
        xGrid(i) = (dble(i) - 1.0d0) / (dble(nx) - 1.0d0)
    end do
    do j = 1, ny
        yGrid(j) = (dble(j) - 1.0d0) / (dble(ny) - 1.0d0)
    end do
    do k = 0, nz+1
        zGrid(k) = 0.5d0 * (erf(constA  * (dble(k) / dble(nz+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
    end do

    ! Compute grid spacing using array slicing
    dz(1:nz+1) = zGrid(1:nz+1) - zGrid(0:nz)
    dx0 = dz(1)
    dt0 = dx0
    ! Calculate viscosity based on LB unit
    length_LB = 1.0d0 / dz(1)

    ! Compute grid spacing in system LB
    zGrid(1:nz) = zGrid(1:nz)-dx0/2.0d0
    zGrid(nz+1) = zGrid(nz+1)-dx0

    H = zGrid(nz+1)/2.0d0
    viscosity = 2.0d0*u_max*H/Reynolds

end subroutine

subroutine loadData()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k
    integer(kind=4) :: ios
    real(kind=8) :: value1, value2, value3, value4, value5, value6, value7
    character(len=100) :: line

    open(unit=01,file='Poiseuille3D-5.dat',form="formatted", access="sequential", status='old', iostat=ios)

    !Check if the file has been successfully opened
    if (ios /= 0) then
        print *, 'Error opening file!'
        stop
    end if

    do i = 1, 3
        read(01, '(A)') line
    end do

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                read(01,*) value1, value2, value3, value4, value5, value6, value7
                u(i, j, k) = value4
                v(i, j, k) = value5
                w(i, j, k) = value6
            end do
        end do
    end do

    close(01)

    return
end subroutine

!subroutine calV()
!    use commondata
!    implicit none
!    integer(kind=4) :: i, j, k
!
!    open(unit=m,file='cal solution- 1.dat',status="unknown")
!    do k = 1, nz
!        write(01,*)  v(nxHalf, nyHalf, k)/u_max, zGrid(k)/zGrid(nz+1)
!    end do
!    close(01)
!
!end subroutine

subroutine solution()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k
    real(kind=8) :: pi=4.0d0*atan(1.0d0)
    real(kind=8) :: tn
    real(kind=8) :: term_sec

    tn = time(5)
    open(unit=j,file='Theoretical solution-5.dat',status="unknown")
    do i=1,nz
        term_sec=0.0d0
        do k=0,300
            term_sec=term_sec+((-1)**k*4.0d0)*cos((k+0.5d0)*pi*(zGrid(i)/H-1.0d0))&
                    *exp(-(k+0.5d0)**2*pi**2*tn)/((k+0.5d0)*pi)**3
        end do
        vana(i)=(1.0d0-(zGrid(i)/H-1.0d0)**2)-term_sec
        write(j,*)  vana(i), zGrid(i)/zGrid(nz+1)
    enddo
    close(j)

end subroutine

subroutine calErr()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k
    real(kind=8) :: errv(nz)

    do k = 1, nz
        errv(k) = 0.0d0
        !do j = 1, ny
           ! do i = 1, nx
                errv(k) = v(nxHalf,nyHalf,k)/u_max-vana(k)
           ! end do
        !end do
    end do

    open(unit=100,file='error v-5.dat',status="unknown")
    do k = 1, nz
        write(100,*)  errv(k), zGrid(k)/zGrid(nz+1)
    end do
    close(100)

end subroutine
