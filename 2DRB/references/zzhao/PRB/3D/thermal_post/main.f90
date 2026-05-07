module commondata
        implicit none
        ! Grid parameters
        integer(kind=4), parameter :: nz = 31, ny = 61, nx = 31
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1, nzHalf=(nz-1)/2+1

        real(kind=8), parameter :: A = 2.0d0      !axial aspect ratios of the channel
        real(kind=8), parameter :: B = 1.0d0     !transverse aspect ratio of the channel

        real(kind=8), parameter :: Ra=100.0d0          !Rayleigh number
        real(kind=8), parameter :: Reynolds=10.0d0      !Reynolds number
        real(kind=8), parameter :: Pr=0.71d0        !Prandtl number
        real(kind=8) :: dt                        !Time step
        real(kind=8) :: length_LB, lengthScale
        real(kind=8) :: U0   !characteristic velocity
        real(kind=8) :: viscosity_LB, tauf
        ! Grid coordinates
        real(kind=8) :: xGrid(nx), yGrid(ny), zGrid(0:nz+1)

        ! Flow variables
        real(kind=8), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)     ! Velocity components
        real(kind=8), allocatable :: T(:,:,:)           ! temperature
        real(kind=8) :: vAna(nz), TAna(nz)
end module commondata

program main
    use commondata
    implicit none

    call initial()

    !call loadData()
    !call getU()
    !call getT()
    call getAna()
    !call getError()
    write(*,*)"End"

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
    allocate (T(nx,ny,nz))

    constA = 3.2d0
     do i = 1, nx
        xGrid(i) = (dble(i) - 1.0d0) / (dble(nx)-1.0d0)
    end do

    do j = 1, ny
        yGrid(j) = (dble(j) - 1.0d0) / (dble(ny)-1.0d0)
    end do

    do k = 0, nz+1
        zGrid(k) = 0.5d0 * (erf(constA  * (dble(k) / dble(nz+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
    end do

    ! Compute grid spacing using array slicing
    dz(1:nz+1) = zGrid(1:nz+1) - zGrid(0:nz)
    !dy(1) = yGrid(2) - yGrid(1)
    ! Compute grid spacing in system LB
    length_LB = 1.0d0 / dz(1)
    zGrid = zGrid-dz(1)/2.0d0
    zGrid(nz+1) = zGrid(nz+1)-dz(1)/2.0d0
    lengthScale = length_LB-1.0d0

    zGrid=zGrid*length_LB
    xGrid=xGrid*lengthScale*B
    yGrid=yGrid*lengthScale*A

    viscosity_LB = (1.0d0/0.8d0-0.5d0)/3.0d0
    U0 = viscosity_LB*Reynolds/(length_LB-1.0d0)
    write(*,*)"u0=",  U0

    return
end subroutine

subroutine loadData()
     use commondata
    implicit none
    integer(kind=4) :: i, j, k
    integer(kind=4) :: ios
    real(kind=8) :: value1, value2, value3, value4, value5, value6, value7
    character(len=100) :: line
    character(len=100) :: filename

    write(filename,*) ny
    filename = adjustl(filename)

    open(unit=01,file='Poiseuille-Benard3D-31.dat',form="formatted", access="sequential", status='old', iostat=ios)

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
                T(i, j, k) = value7
            end do
        end do
    end do

    close(01)

    return
end subroutine

subroutine getU()
    use commondata
    implicit none

    integer(kind=4) :: k

    open(unit=03,file='v.dat',status='unknown')
        do k = 1, nz
            write(03,*) v(nxHalf, nyHalf, k), zGrid(k)/lengthScale
        end do
    close(03)

    return
end subroutine

subroutine getT()
    use commondata
    implicit none

    integer(kind=4) :: k

    open(unit=04,file='T.dat',status='unknown')
        do k = 1, nz
            write(04,*) T(nxHalf, nyHalf, k), zGrid(k)/lengthScale
        end do
    close(04)

    return
end subroutine

subroutine getAna()
    use commondata
    implicit none

    integer(kind=4) :: k

    open(unit=05,file='T_ana.dat',status='unknown')
        do k = 1, nz
            write(05,*) (exp(Pr*Reynolds*zGrid(k)/lengthScale)-1.0d0)/(exp(Pr*Reynolds)-1.0d0), zGrid(k)/lengthScale
        end do
    close(05)

    open(unit=06,file='v_ana.dat',status='unknown')
        do k = 1, nz
            write(06,*) U0*(exp(Reynolds*zGrid(k)/lengthScale)-1.0d0)/(exp(Reynolds)-1.0d0), zGrid(k)/lengthScale
        end do
    close(06)

    return
end subroutine

!subroutine getError()
!    use commondata
!    implicit none
!
!    integer(kind=4) :: i, j
!    real(kind=8) :: errU, errT
!    real(kind=8) :: errx_u, errx_T
!    real(kind=8) :: err2U, err2T
!    real(kind=8) :: integrate
!
!    do j = 1, ny
!        uAna(j) = U0*(exp(Reynolds*yGrid(j)/lengthScale)-1.0d0)/(exp(Reynolds)-1.0d0)
!        TAna(j) = (exp(Pr*Reynolds*yGrid(j)/lengthScale)-1.0d0)/(exp(Pr*Reynolds)-1.0d0)
!    end do
!
!    errU = 0.0d0
!    errT = 0.0d0
!    errx_u = 0.0d0
!    errx_T = 0.0d0
!    err2U = 0.0d0
!    err2T = 0.0d0
!
!    do j = 1, ny-1
!        errU = errU+(yGrid(j+1)-yGrid(j))*((uAna(j)-u(nxHalf,j))**2+(uAna(j+1)-u(nxHalf,j+1))**2)/2.0d0
!        errT = errT+(yGrid(j+1)-yGrid(j))*((TAna(j)-T(nxHalf,j))**2+(TAna(j+1)-T(nxHalf,j+1))**2)/2.0d0
!        errx_u = errx_u+(yGrid(j+1)-yGrid(j))*((uAna(j))**2+(uAna(j+1))**2)/2.0d0
!        errx_T = errx_T+(yGrid(j+1)-yGrid(j))*((TAna(j))**2+(TAna(j+1))**2)/2.0d0
!    end do
!    err2U = dsqrt(errU)/dsqrt(errx_u)
!    err2T = dsqrt(errT)/dsqrt(errx_T)
!    write(*,*)"err2U=", err2U, "err2T=", err2T
!
!    open(unit=07,file='errorU.dat',position='append')
!        write(07,*) (2.0d0*yGrid(1))/lengthScale, err2U
!    close(07)
!
!    open(unit=08,file='errorT.dat',position='append')
!        write(08,*) (2.0d0*yGrid(1))/lengthScale, err2T
!    close(08)
!
!    open(unit=02,file='y=2x.dat',status='unknown')
!        do j = 1, 7
!            write(02,*) yGrid(j)/lengthScale, 10**(2.0d0*log10(yGrid(j)/lengthScale)+3.3d0)
!        end do
!    close(02)
!    return
!end subroutine
