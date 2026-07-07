module commondata
        implicit none
        ! Grid parameters
        integer(kind=4), parameter :: ny=65, nx=ny*20+1
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1

        integer(kind=4), parameter :: B = 20          !aspect ratio

        real(kind=8) :: dt                        !Time step
        real(kind=8) :: length_LB
        real(kind=8), parameter :: Ur = 0.01d0   !characteristic velocity

        ! Grid coordinates
        real(kind=8) :: xGrid(0:nx+1), yGrid(0:ny+1)

        ! Flow variables
        real(kind=8), allocatable :: u(:,:), v(:,:)      ! Velocity components
        real(kind=8), allocatable :: T(:,:)           ! temperature
        real(kind=8), allocatable :: P(:,:)            !pressure
        real(kind=8) :: vorticity(nx, ny)
        real(kind=8) :: Nu_Avg=0.0d0, Nu_space

end module commondata

program main
    use commondata
    implicit none

    call initial()

    call loadData()
    call getspeed_extreme()
    call calNu()
    call get_x5_Info()
    call get_x10_Info()
    call getVorticity()
    call streaming_function()
   ! Nu_Avg = Nu_Avg+Nu_space

    deallocate (u)
    deallocate (v)
    deallocate (T)
    deallocate (P)

    return
end program main

subroutine initial()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: constA
    real(kind=8) :: dy(ny+1)

    allocate (u(nx,ny))
    allocate (v(nx,ny))
    allocate (T(nx,ny))
    allocate (P(nx,ny))

    constA = 3.2d0
    do i = 1, nx
        xGrid(i) = (dble(i) - 1.0d0) / (dble(nx) - 1.0d0)
    end do
    do j = 0, ny+1
        yGrid(j) = 0.5d0 * (erf(constA  * (dble(j) / dble(ny+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
    end do

    ! Compute grid spacing using array slicing
    dy(1:ny+1) = yGrid(1:ny+1) - yGrid(0:ny)

    ! Compute grid spacing in system LB
    length_LB = 1.0d0 / dy(1)
    yGrid(1:ny+1) = yGrid(1:ny+1)-dy(1)/2.0d0
    yGrid(ny+1) = yGrid(ny+1)-yGrid(1)
    xGrid=xGrid*(length_LB-1.0d0)*dble(B)
    yGrid=yGrid*length_LB

    xGrid=xGrid/(yGrid(ny)+yGrid(1))
    yGrid=yGrid/(yGrid(ny)+yGrid(1))
    write(*,*)"y1=",yGrid(1),"y(ny)=",yGrid(ny+1)
    write(*,*)"x1=",xGrid(1),"x(nx)=",xGrid(nx)

    return
end subroutine

subroutine loadData()
     use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: ios
    real(kind=8) :: value1, value2, value3, value4, value5, value6
    character(len=100) :: filename
    character(len=100) :: line

    write(filename,*) dt/((yGrid(ny)+yGrid(1))/Ur)
    filename = adjustl(filename)
    open(unit=01,file='MRTcavity-5982440.dat',form="formatted", access="sequential", status='old', iostat=ios)
    !open(unit=01, file=trim(filename)//'.dat', form="formatted", access="sequential", status='old', iostat=ios)

    !Check if the file has been successfully opened
    if (ios /= 0) then
        print *, 'Error opening file!'
        stop
    end if

    do i = 1, 3
        read(01, '(A)') line
    end do

    do j = 1, ny
        do i = 1, nx
            read(01,*) value1, value2, value3, value4, value5, value6
            u(i, j) = value3 / Ur
            v(i, j) = value4 / Ur
            T(i, j) = value5
            P(i, j) = value6
        end do
    end do

    close(01)

    return
end subroutine

subroutine get_x5_Info()
    use commondata
    implicit none
    integer(kind=4) :: j
    integer(kind=4) :: localx = (nxHalf-1)/2+1
    real(kind=8) :: x1, x2

    open(unit=03, file='T_x5.dat', status='unknown', action='write')
        do j = 1, ny
            write(03, *) yGrid(j), T(localx,j)
        end do
    close(03)

    open(unit=04, file='U_x5.dat', status='unknown', action='write')
        do j = 1, ny
            write(04, *) yGrid(j), u(localx,j)
        end do
    close(04)

    open(unit=05, file='V_x5.dat', status='unknown', action='write')
        do j = 1, ny
            write(05, *) yGrid(j), v(localx,j)
        end do
    close(05)

    open(unit=09, file='P_x5.dat', status='unknown', action='write')
        do j = 1, ny
            write(09, *) yGrid(j), (P(localx,j)-1.0d0/3.0d0)/(Ur)**2
        end do
    close(09)

!--------------------------------------derivatives---------------------------------------------
    x1 = xGrid(localx)-xGrid(localx-1)
    x2 = xGrid(localx+1)-xGrid(localx)

    open(unit=11, file='dudx_5.dat', status='unknown', action='write')
        do j = 1, ny
            write(11, *) yGrid(j), x2*(u(localx-1,j)-x1**2*u(localx+1,j)/x2**2+&
            (x1**2/x2**2-1.0d0)*u(localx,j))/(-x1*x2-x1**2)
        end do
    close(11)


    open(unit=12, file='dvdx_5.dat', status='unknown', action='write')
        do j = 1, ny
            write(12, *) yGrid(j), x2*(v(localx-1,j)-x1**2*v(localx+1,j)/x2**2+&
            (x1**2/x2**2-1.0d0)*v(localx,j))/(-x1*x2-x1**2)
        end do
    close(12)

    open(unit=13, file='dTdx_5.dat', status='unknown', action='write')
        do j = 1, ny
            write(13, *) yGrid(j), x2*(T(localx-1,j)-x1**2*T(localx+1,j)/x2**2+&
            (x1**2/x2**2-1.0d0)*T(localx,j))/(-x1*x2-x1**2)
        end do
    close(13)

    return
end subroutine get_x5_Info

subroutine get_x10_Info()
    use commondata
    implicit none
    integer(kind=4) :: j
    real(kind=8) :: x1, x2

    open(unit=06, file='T_x10.dat', status='unknown', action='write')
        do j = 1, ny
            write(06, *) yGrid(j), T(nxHalf,j)
        end do
    close(06)

    open(unit=07, file='U_x10.dat', status='unknown', action='write')
        do j = 1, ny
            write(07, *) yGrid(j), u(nxHalf,j)
        end do
    close(07)

    open(unit=08, file='V_x10.dat', status='unknown', action='write')
        do j = 1, ny
            write(08, *) yGrid(j), v(nxHalf,j)
        end do
    close(08)

    open(unit=10, file='P_x10.dat', status='unknown', action='write')
        do j = 1, ny
            write(10, *) yGrid(j),  (P(nxHalf,j)-1.0d0/3.0d0)/(Ur)**2
        end do
    close(10)
!--------------------------------------derivatives---------------------------------------------
    x1 = xGrid(nxHalf)-xGrid(nxHalf-1)
    x2 = xGrid(nxHalf+1)-xGrid(nxHalf)

    open(unit=21, file='dudx_10.dat', status='unknown', action='write')
        do j = 1, ny
            write(21, *) yGrid(j), x2*(u(nxHalf-1,j)-x1**2*u(nxHalf+1,j)/x2**2+&
            (x1**2/x2**2-1.0d0)*u(nxHalf,j))/(-x1*x2-x1**2)
        end do
    close(21)

    open(unit=22, file='dvdx_10.dat', status='unknown', action='write')
        do j = 1, ny
            write(22, *) yGrid(j), x2*(v(nxHalf-1,j)-x1**2*v(nxHalf+1,j)/x2**2+&
            (x1**2/x2**2-1.0d0)*v(nxHalf,j))/(-x1*x2-x1**2)
        end do
    close(22)

    open(unit=23, file='dTdx_10.dat', status='unknown', action='write')
        do j = 1, ny
            write(23, *) yGrid(j), x2*(T(nxHalf-1,j)-x1**2*T(nxHalf+1,j)/x2**2+&
            (x1**2/x2**2-1.0d0)*T(nxHalf,j))/(-x1*x2-x1**2)
        end do
    close(23)

    return
end subroutine get_x10_Info

subroutine getspeed_extreme()
    use commondata
    implicit none

    integer(kind=4) :: i, j
    real(kind=8) :: u_max, u_min, v_max, v_min

    u_max = maxval(u(1:nxHalf,:))
    u_min = minval(u(1:nxHalf,:))
    v_max = maxval(v(1:nxHalf,:))
    v_min = minval(v(1:nxHalf,:))
    write(*,*)"u_max=",u_max,"u_min=",u_min
    write(*,*)"v_max=",v_max,"v_min=",v_min
end subroutine

subroutine calNu()
    use commondata
    implicit none

    integer(kind=4) :: i, j
    real(kind=8) :: delta_y1, delta_y2
    real(kind=8) :: Nu(nxHalf,ny)


    !bottom
    do i = 1, nxHalf
        j = 1
        delta_y1 = yGrid(j+1)-yGrid(j)
        delta_y2 = yGrid(j+2)-yGrid(j)

        Nu(i,j) = delta_y2*(T(i,j+1)-delta_y1**2*T(i,j+2)/delta_y2**2+(delta_y1**2/delta_y2**2-1.0d0)*T(i,j))&
            /(delta_y1*delta_y2-delta_y1**2)
    end do

    !top
    do i = 1, nxHalf
        j = ny
        delta_y1 = yGrid(j)-yGrid(j-1)
        delta_y2 = yGrid(j)-yGrid(j-2)

        Nu(i,j) = -delta_y2*(T(i,j-1)-delta_y1**2*T(i,j-2)/delta_y2**2+(delta_y1**2/delta_y2**2-1.0d0)*T(i,j))&
            /(delta_y1*delta_y2-delta_y1**2)
    end do

    Nu_space = 0.0d0
    do i = 1, nxHalf
        Nu_space = Nu_space-(Nu(i,1)+Nu(i,ny))/2.0d0
    end do

    Nu_space = Nu_space / (dble(nxHalf)-1.0d0)
    write(*,*)"Nu_space=",Nu_space
end subroutine

subroutine getVorticity()
    use commondata
    implicit none

    integer(kind=4) :: i, j
    real(kind=8) :: dvdx(nx,ny), dudy(nx,ny)
    real(kind=8) :: delta1, delta2
    character(len=100) :: filename
    integer(kind=4) :: localx = (nxHalf-1)/2+1

    do j = 1, ny
        do i = 1, nx
            if(i == 1)then
                delta1 = xGrid(i+1)-xGrid(i)
                delta2 = xGrid(i+2)-xGrid(i)
                dvdx(i,j) = delta2*(v(i+1,j)-delta1**2.0d0*v(i+2,j)/delta2**2.0d0+(delta1**2.0d0/delta2**2.0d0&
                            -1.0d0)*v(i,j))/(delta1*delta2-delta1**2.0d0)

            elseif(i == nx)then
                delta1 = xGrid(i)-xGrid(i-1)
                delta2 = xGrid(i)-xGrid(i-2)
                dvdx(i,j) = -delta2*(v(i-1,j)-delta1**2.0d0*v(i-2,j)/delta2**2.0d0+(delta1**2.0d0/delta2**2.0d0&
                            -1.0d0)*v(i,j))/(delta1*delta2-delta1**2.0d0)

            else
                delta1 = xGrid(i)-xGrid(i-1)
                delta2 = xGrid(i+1)-xGrid(i)
                dvdx(i,j) = delta2*(v(i-1,j)-delta1**2.0d0*v(i+1,j)/delta2**2.0d0+(delta1**2.0d0/delta2**2.0d0&
                            -1.0d0)*v(i,j))/(-delta1*delta2-delta1**2.0d0)
            end if

            if(j == 1)then
                delta1 = yGrid(j+1)-yGrid(j)
                delta2 = yGrid(j+2)-yGrid(j)
                dudy(i,j) = delta2*(u(i,j+1)-delta1**2.0d0*u(i,j+2)/delta2**2.0d0+(delta1**2.0d0/delta2**2.0d0&
                            -1.0d0)*u(i,j))/(delta1*delta2-delta1**2.0d0)

            elseif(j == ny)then
                delta1 = yGrid(j)-yGrid(j-1)
                delta2 = yGrid(j)-yGrid(j-2)
                dudy(i,j) = -delta2*(u(i,j-1)-delta1**2.0d0*u(i,j-2)/delta2**2.0d0+(delta1**2.0d0/delta2**2.0d0&
                            -1.0d0)*u(i,j))/(delta1*delta2-delta1**2.0d0)

            else
                delta1 = yGrid(j)-yGrid(j-1)
                delta2 = yGrid(j+1)-yGrid(j)
                dudy(i,j) = delta2*(u(i,j-1)-delta1**2.0d0*u(i,j+1)/delta2**2.0d0+(delta1**2.0d0/delta2**2.0d0&
                            -1.0d0)*u(i,j))/(-delta1*delta2-delta1**2.0d0)
            end if

            vorticity(i,j) = dvdx(i,j)-dudy(i,j)
        end do
    end do

    open(unit=03, file='vorticity_5.dat', status='unknown', action='write')
        do j = 1, ny
            write(03, *) yGrid(j), vorticity(localx,j)
        end do
    close(03)


    open(unit=04, file='vorticity_10.dat', status='unknown', action='write')
        do j = 1, ny
            write(04, *) yGrid(j), vorticity(nxHalf,j)
        end do
    close(04)


    open(unit=02,file='vorticity.dat',status='unknown')

    write(02,*) 'TITLE="Flow field vorticity"'
    write(02,*) 'VARIABLES="X" "Y" "vorticity" '
    write(02,101) nx, ny
    do j=1,ny
        do i=1,nx
            write(02,100) xGrid(i), yGrid(j), vorticity(i,j)
        enddo
    enddo
100 format(1x,2(e11.4,' '),(e13.6,' '))
101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
    close(02)

end subroutine

subroutine streaming_function()
    use commondata
    implicit none

    integer(kind=4) :: i, j
    real(kind=8) :: phi(nx,ny), coeff_integral(ny)

    do j = 1, ny
        coeff_integral(j) = (yGrid(j+1)-yGrid(j-1))/2.0d0
    end do

    phi = 0.0d0
    do j = 1, ny
        do i = 1, nx
           phi(i,j) = u(i,j)*coeff_integral(j)+phi(i,j-1)
        end do
    end do


 open(unit=01,file='streaming_function.dat',status='unknown')

    write(01,*) 'TITLE="Flow field streaming_function"'
    write(01,*) 'VARIABLES="X" "Y" "streaming_function" '
    write(01,101) nx, ny
    do j=1,ny
        do i=1,nx
            write(01,100) xGrid(i), yGrid(j), phi(i,j)
        enddo
    enddo
100 format(1x,2(e11.4,' '),(e13.6,' '))
101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
    close(01)
    return
end subroutine
