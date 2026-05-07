module commondata
        implicit none
        ! Grid parameters
        integer(kind=4), parameter :: nz = 41, ny = 601, nx = 121
        integer(kind=4), parameter :: y10=(ny-1)/25+(ny-1)/5+1, y30=(ny-1)/25+3*(ny-1)/5+1
        integer(kind=4), parameter :: nxHalf = (nx-1)/2+1
        real(kind=8), parameter :: A = 50.0d0      !axial aspect ratios of the channel
        real(kind=8), parameter :: B = 10.0d0     !transverse aspect ratio of the channel
        real(kind=8), parameter :: Ae = 2.0d0
        real(kind=8) :: viscosity_LB
        real(kind=8) :: dt                        !Time step
        real(kind=8) :: length_LB
        real(kind=8) :: lengthScale

        ! Grid coordinates
        real(kind=8) :: xGrid(0:nx+1), yGrid(ny), zGrid(0:nz+1)

        ! Flow variables
        real(kind=8), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)      ! Velocity components
        real(kind=8), allocatable :: T(:,:,:), P(:,:,:), rho(:,:,:)
        real(8) :: U0, Umean
        real(kind=8) :: coeff_integral(ny)                                       ! Integral coefficient
end module commondata

program main
    use commondata
    implicit none

    call initial()
    call loadData()
    call get_velocity()
    call get_temperture()
    call energy()
    call temperature()
    call pressure()
    call get_Nu()

    deallocate (u)
    deallocate (v)
    deallocate (w)
    deallocate (T)
    deallocate (P)
    deallocate (rho)
    write(*,*)"end"
    return
end program main

subroutine initial()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k, n
    real(kind=8) :: constA, coeff_v
    real(kind=8) :: dz(nz+1)
    real(kind=8) :: pi = 4.0d0*atan(1.0d0)
    allocate (u(nx,ny,nz))
    allocate (v(nx,ny,nz))
    allocate (w(nx,ny,nz))
    allocate (T(nx,ny,nz))
    allocate (P(nx,ny,nz))
    allocate (rho(nx,ny,nz))
    constA = 3.2d0
    do i = 1, nx
        xGrid(i) = (dble(i) - 0.5d0) / dble(nx)
    end do
    xGrid(nx+1) = 1.0d0
    do j = 1, ny
        yGrid(j) = (dble(j) - 1.0d0) / (dble(ny)-1.0d0)
    end do

    do k = 0, nz+1
        zGrid(k) = 0.5d0 * (erf(constA  * (dble(k) / dble(nz+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
    end do

    ! Compute grid spacing using array slicing
    dz(1:nz+1) = zGrid(1:nz+1) - zGrid(0:nz)

    ! Compute grid spacing in system LB
    length_LB = 1.0d0 / dz(1)

    lengthScale = length_LB-1.0d0

   ! viscosity_LB = U0*lengthScale/Reynolds

    ! Compute grid spacing in system LB
    zGrid(1:nz) = zGrid(1:nz)-dz(1)/2.0d0
    zGrid=zGrid*length_LB
    xGrid=xGrid*lengthScale*B
    yGrid=yGrid*lengthScale*A
    write(*,*)"length_z=", zGrid(nz)+zGrid(1), "length_x=", xGrid(nx+1), "length_y=", yGrid(ny)

     U0 = 0.02d0
    ! Calculate relaxation time
    coeff_v = 0.0d0
    do n = 0, 25
        coeff_v = coeff_v+(-1)**(n+1)*sin((2.0d0*dble(n)+1.0d0)*pi/2.0d0)*(1.0d0-exp(-2.0d0*(2.0d0*dble(n)+1.0d0)&
                *pi*B/2.0d0))/((2.0d0*dble(n)+1.0d0)**5*B*(1.0d0+exp(-2.0d0*(2.0d0*dble(n)+1.0d0)*pi*B/2.0d0)))
    end do
    Umean = (1.0d0+192.0d0/pi**5*coeff_v)*U0
    write(*,*)"Umean=",Umean
    return
end subroutine

subroutine loadData()
     use commondata
    implicit none
    integer(kind=4) :: i, j, k
    integer(kind=4) :: ios
    real(kind=8) :: value1, value2, value3, value4, value5, value6, value7, value8, value9
    character(len=100) :: line
    character(len=100) :: filename

    write(filename,*) nz
    filename = adjustl(filename)

    open(unit=01,file='Poiseuille-Benard3DCBC.dat',form="formatted", access="sequential", status='old', iostat=ios)

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
                read(01,*) value1, value2, value3, value4, value5, value6, value7,value8, value9
                u(i, j, k) = value4
                v(i, j, k) = value5
                w(i, j, k) = value6
                T(i, j, k) = value7
                P(i, j, k) = value8
                rho(i,j,k) = value9
            end do
        end do
    end do

    close(01)

    return
end subroutine

subroutine get_velocity()
    use commondata
    implicit none
    integer(kind=4) :: i
!----------------------------------------------------Flow velocity-----------------------------------------
    open(unit=03,file='v(x=10,z=0.5).dat',status='unknown')
        do i = 1, nxHalf
            write(03,*) xGrid(i)/lengthScale, v(i,(ny-1)/25+(ny-1)/5+1,(nz-1)/2+1)/Umean
        end do
    close(03)


    open(unit=02,file='v(x=30,z=0.5).dat',status='unknown')
        do i = 1, nxHalf
            write(02,*) xGrid(i)/lengthScale, v(i,(ny-1)/25+3*(ny-1)/5+1,(nz-1)/2+1)/Umean
        end do
    close(02)


!    open(unit=05,file='v(x=10,z=0.2).dat',status='unknown')
!        do i = 1, nxHalf
!            write(05,*) xGrid(i)/lengthScale, v(i,(ny-1)/25+(ny-1)/5+1,(nz-1)/5+1)/Umean
!        end do
!    close(05)
!
!
!    open(unit=07,file='v(x=30,z=0.2).dat',status='unknown')
!        do i = 1, nxHalf
!            write(07,*) xGrid(i)/lengthScale, v(i,(ny-1)/25+3*(ny-1)/5+1,(nz-1)/5+1)/Umean
!        end do
!    close(07)
!---------------------------------------vertical velocity---------------------------------------------------
    open(unit=13,file='w(x=10,z=0.5).dat',status='unknown')
        do i = 1, nxHalf
            write(13,*) xGrid(i)/lengthScale, w(i,(ny-1)/25+(ny-1)/5+1,(nz-1)/2+1)/Umean
        end do
    close(13)


    open(unit=12,file='w(x=30,z=0.5).dat',status='unknown')
        do i = 1, nxHalf
            write(12,*) xGrid(i)/lengthScale, w(i,(ny-1)/25+3*(ny-1)/5+1,(nz-1)/2+1)/Umean
        end do
    close(12)


!    open(unit=15,file='w(x=10,z=0.2).dat',status='unknown')
!        do i = 1, nxHalf
!            write(15,*) xGrid(i)/lengthScale, w(i,(ny-1)/25+(ny-1)/5+1,(nz-1)/5+1)/Umean
!        end do
!    close(15)
!
!
!    open(unit=17,file='w(x=30,z=0.2).dat',status='unknown')
!        do i = 1, nxHalf
!            write(17,*) xGrid(i)/lengthScale, w(i,(ny-1)/25+3*(ny-1)/5+1,(nz-1)/5+1)/Umean
!        end do
!    close(17)

!------------------------------------------Exhibition speed------------------------------
    open(unit=23,file='u(x=10,z=0.5).dat',status='unknown')
        do i = 1, nxHalf
            write(23,*) xGrid(i)/lengthScale, u(i,(ny-1)/25+(ny-1)/5+1,(nz-1)/2+1)/Umean
        end do
    close(23)


    open(unit=22,file='u(x=30,z=0.5).dat',status='unknown')
        do i = 1, nxHalf
            write(22,*) xGrid(i)/lengthScale, u(i,(ny-1)/25+3*(ny-1)/5+1,(nz-1)/2+1)/Umean
        end do
    close(22)


!    open(unit=25,file='u(x=10,z=0.2).dat',status='unknown')
!        do i = 1, nxHalf
!            write(25,*) xGrid(i)/lengthScale, u(i,(ny-1)/25+(ny-1)/5+1,(nz-1)/5+1)/Umean
!        end do
!    close(25)
!
!
!    open(unit=27,file='u(x=30,z=0.2).dat',status='unknown')
!        do i = 1, nxHalf
!            write(27,*) xGrid(i)/lengthScale, u(i,(ny-1)/25+3*(ny-1)/5+1,(nz-1)/5+1)/Umean
!        end do
!    close(27)
end subroutine

subroutine get_temperture()
    use commondata
    implicit none

     integer(kind=4) :: i

    open(unit=33,file='T(x=10,z=0.5).dat',status='unknown')
        do i = 1, nxHalf
            write(33,*) xGrid(i)/lengthScale, T(i,(ny-1)/25+(ny-1)/5+1,(nz-1)/2+1)
        end do
    close(33)


    open(unit=32,file='T(x=30,z=0.5).dat',status='unknown')
        do i = 1, nxHalf
            write(32,*) xGrid(i)/lengthScale, T(i,(ny-1)/25+3*(ny-1)/5+1,(nz-1)/2+1)
        end do
    close(32)


!    open(unit=35,file='T(x=10,z=0.2).dat',status='unknown')
!        do i = 1, nxHalf
!            write(35,*) xGrid(i)/lengthScale, T(i,(ny-1)/25+(ny-1)/5+1,(nz-1)/5)
!        end do
!    close(35)

!
!    open(unit=37,file='T(x=30,z=0.2).dat',status='unknown')
!        do i = 1, nxHalf
!            write(37,*) xGrid(i)/lengthScale, T(i,(ny-1)/25+3*(ny-1)/5+1,(nz-1)/5)
!        end do
!    close(37)
end subroutine

subroutine get_Nu()
    use commondata
    implicit none

    real(kind=8) :: Nut_10(nxHalf), Nut_30(nxHalf), Nub_10(nxHalf), Nub_30(nxHalf)
    integer(kind=8) :: i
    real(kind=8) :: delta1, delta2

    delta1 = zGrid(2)/lengthScale-zGrid(1)/lengthScale
    delta2 = zGrid(3)/lengthScale-zGrid(1)/lengthScale
    do i = 1, nxHalf
        Nub_10(i) = delta2*(T(i,y10,2)-delta1**2*T(i,y10,3)/delta2**2+(delta1**2/delta2**2-1.0d0)*T(i,y10,1))&
                    /(delta1*delta2-delta1**2)
        Nub_30(i) = delta2*(T(i,y30,2)-delta1**2*T(i,y30,3)/delta2**2+(delta1**2/delta2**2-1.0d0)*T(i,y30,1))&
            /(delta1*delta2-delta1**2)
    end do

    delta1 = zGrid(nz)/lengthScale-zGrid(nz-1)/lengthScale
    delta2 = zGrid(nz)/lengthScale-zGrid(nz-2)/lengthScale
    do i = 1, nxHalf
        Nut_10(i) = -delta2*(T(i,y10,nz-1)-delta1**2*T(i,y10,nz-2)/delta2**2+(delta1**2/delta2**2-1.0d0)*T(i,y10,nz))&
            /(delta1*delta2-delta1**2)
        Nut_30(i) = -delta2*(T(i,y30,nz-1)-delta1**2*T(i,y30,nz-2)/delta2**2+(delta1**2/delta2**2-1.0d0)*T(i,y30,nz))&
            /(delta1*delta2-delta1**2)
    end do

    open(unit=43,file='Nub_10.dat',status='unknown')
        do i = 1, nxHalf
            write(43,*) xGrid(i)/lengthScale, -Nub_10(i)
        end do
    close(43)


    open(unit=42,file='Nub_30.dat',status='unknown')
        do i = 1, nxHalf
            write(42,*) xGrid(i)/lengthScale,  -Nub_30(i)
        end do
    close(42)


    open(unit=45,file='Nut_10.dat',status='unknown')
        do i = 1, nxHalf
            write(45,*) xGrid(i)/lengthScale, -Nut_10(i)
        end do
    close(45)


    open(unit=47,file='Nut_30.dat',status='unknown')
        do i = 1, nxHalf
            write(47,*) xGrid(i)/lengthScale, -Nut_30(i)
        end do
    close(47)

    return
end subroutine

subroutine energy()
    use commondata
    implicit none

    integer(kind=4) :: i, j, k
    real(kind=8) :: velocitysquare(nx,ny,nz)
    real(kind=8) :: inter_z(nx,ny), inter_y(nx)
    real(kind=8) :: kenergy

    do k = 1, nz
        do j = (ny-1)/25, ny
            do i = 1, nx
                velocitysquare(i,j,k) = (u(i,j,k)/Umean)**2+(v(i,j,k)/Umean)**2+(w(i,j,k)/Umean)**2
            end do
        end do
    end do

    inter_z = 0.0d0
    do k = 1, nz-1
        do j = (ny-1)/25, ny
            do i = 1, nx
                inter_z(i,j) = inter_z(i,j)+(zGrid(k+1)-zGrid(k))*(velocitysquare(i,j,k)+velocitysquare(i,j,k+1))/2.0d0
            end do
        end do
    end do

    inter_y = 0.0d0
    do j = (ny-1)/25, ny-1
        do i = 1, nx
            inter_y(i) = inter_y(i)+(yGrid(j+1)-yGrid(j))*(inter_z(i,j)+inter_z(i,j+1))/2.0d0
        end do
    end do

    kenergy = 0.0d0
    do i = 1, nx-1
        kenergy = kenergy+(xGrid(i+1)-xGrid(i))*(inter_y(i)+inter_y(i+1))/2.0d0
    end do

    kenergy = kenergy/(lengthScale*xGrid(nx+1)*(yGrid(ny)-yGrid((ny-1)/25)))
    write(*,*)"energy=", kenergy

    return
end subroutine

subroutine temperature()
    use commondata
    implicit none

    integer(kind=4) :: i, j, k
    real(kind=8) :: inter_z(nx,ny), inter_y(nx)
    real(kind=8) :: Tmean

    inter_z = 0.0d0
    do k = 1, nz-1
        do j = (ny-1)/25, ny
            do i = 1, nx
                inter_z(i,j) = inter_z(i,j)+(zGrid(k+1)-zGrid(k))*(T(i,j,k)+T(i,j,k+1))/2.0d0
            end do
        end do
    end do

    inter_y = 0.0d0
    do j = (ny-1)/25, ny-1
        do i = 1, nx
            inter_y(i) = inter_y(i)+(yGrid(j+1)-yGrid(j))*(inter_z(i,j)+inter_z(i,j+1))/2.0d0
        end do
    end do

    Tmean = 0.0d0
    do i = 1, nx-1
        Tmean = Tmean+(xGrid(i+1)-xGrid(i))*(inter_y(i)+inter_y(i+1))/2.0d0
    end do

    Tmean = Tmean/(lengthScale*xGrid(nx+1)*(yGrid(ny)-yGrid((ny-1)/25)))

    write(*,*)"Tmean=", Tmean

    return
end subroutine

subroutine pressure()
    use commondata
    implicit none

    integer(kind=4) :: i, j, k
    real(kind=8) :: P0, P1 ,diff_P
    real(kind=8) :: inter_z(nx)

    inter_z = 0.0d0
    do k = 1, nz-1
        do i = 1, nx
            inter_z(i) = inter_z(i)+(zGrid(k+1)-zGrid(k))*(P(i,(ny-1)/25,k)/(Umean**2)&
                        +P(i,(ny-1)/25,k+1)/(Umean**2))/2.0d0
        end do
    end do

    P0 = 0.0d0
    do i = 1, nx-1
        P0 = P0+(xGrid(i+1)-xGrid(i))*(inter_z(i)+inter_z(i+1))/2.0d0
    end do

    inter_z = 0.0d0
    do k = 1, nz-1
        do i = 1, nx
            inter_z(i) = inter_z(i)+(zGrid(k+1)-zGrid(k))*(P(i,ny,k)/(Umean**2)&
                        +P(i,ny,k+1)/(Umean**2))/2.0d0
        end do
    end do

    P1 = 0.0d0
    do i = 1, nx-1
        P1 = P1+(xGrid(i+1)-xGrid(i))*(inter_z(i)+inter_z(i+1))/2.0d0
    end do

    diff_P = (P0-P1)/(lengthScale*xGrid(nx+1))

    write(*,*)"diff_P=",diff_P
    return
end subroutine
