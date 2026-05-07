module commondata
        implicit none
        ! Grid parameters
        integer(kind=4), parameter :: ny=33, nx=ny*20+1
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1

        real(kind=8), parameter :: Ra=1e4          !Rayleigh number
        real(kind=8), parameter :: Pr=2.0d0/3.0d0        !Prandtl number
        real(kind=8), parameter :: Reynolds=10.0d0
        real(kind=8), parameter :: Thot=1.0d0, Tcold=0.0d0, Tref=0.0d0
        real(kind=8) :: kappa                        !Thermal expansion coefficient
        real(kind=8) :: gbeta                        !Volume expansion coefficient * Gravitational acceleration
        integer(kind=4), parameter :: B = 20          !aspect ratio

        ! Parameters in 0 system
        real(kind=8), parameter :: rho0 = 1.0d0
        real(kind=8) :: dt0                       ! Time step
        real(kind=8) :: dx0                       ! Grid spacing
        real(kind=8) :: viscosity0
        real(kind=8) :: kappa0
        ! Parameters in LB system
        real(kind=8) :: length_LB                 ! Characteristic length
        real(kind=8) :: viscosity_LB              ! Kinematic viscosity
        real(kind=8) :: tauf, taut                ! Relaxation time
        real(kind=8) :: dt                        !Time step
        real(kind=8) :: ax
        real(kind=8), parameter :: Ur = 0.01d0   !characteristic velocity
        real(kind=8) :: u_ave                     !Average entrance speed
        real(kind=8) :: lemda                     !CBC coefficient
        ! Iteration control
        integer(kind=4) :: itc = 0            ! Current iteration count
        integer(kind=4) :: itc_max              ! Maximum iterations

        ! Grid coordinates
        real(kind=8) :: xGrid(nx), yGrid(0:ny+1)

        ! Flow variables
        real(kind=8), allocatable :: u(:,:), v(:,:)      ! Velocity components
        real(kind=8), allocatable :: temp(:,:)           ! temperature
        real(kind=8), allocatable :: rho(:,:)            ! Density
        real(kind=8), allocatable :: Fx(:,:), Fy(:,:)    ! force components

        ! Distribution functions
        real(kind=8), allocatable :: f(:,:,:), f_post(:,:,:)  ! Current and post-collision distributions
        real(kind=8), allocatable :: g(:,:,:), g_post(:,:,:)
        real(kind=8) :: f_last(ny,0:8)

        ! MRT relaxation parameters
        real(kind=8) :: omega_U(0:8), omega_T(0:4)  ! Relaxation rates for MRT

        ! Lattice directions
        integer(kind=4) :: ex(0:8), ey(0:8)  ! Lattice velocity directions
        data ex/0, 1, 0, -1,  0, 1, -1, -1,  1/
        data ey/0, 0, 1,  0, -1, 1,  1, -1, -1/

        ! Additional MRT parameters
        real(kind=8) :: Snu, Sq, sig_k

        integer(kind=4) :: inter_x(nx,3), inter_y(ny,3)

        real(kind=8), parameter :: outputFrequency=1.0d0
        integer(kind=4), parameter :: dimensionlessTimeMax=int(12000/outputFrequency)
        real(kind=8) :: timeUnit
end module commondata

subroutine initial()
    use commondata
    implicit none
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: un(0:8)
    real(kind=8) :: velocitySquare
    real(kind=8) :: dy(ny+1)
    real(kind=8) :: constA

    ! Allocate flow variables
    allocate (u(nx,ny))
    allocate (v(nx,ny))
    allocate (rho(nx,ny))
    allocate (temp(nx,ny))
    allocate (Fx(nx,ny))
    allocate (Fy(nx,ny))

    allocate (f(nx,ny,0:8))
    allocate (f_post(nx,ny,0:8))
    allocate (g(nx,ny,0:4))
    allocate (g_post(nx,ny,0:4))

    ! Compute grid coordinates
    constA = 3.2d0
    do i = 1, nx
        xGrid(i) = (dble(i) - 1.0d0) / (dble(nx) - 1.0d0)
    end do
    do j = 0, ny+1
        yGrid(j) = 0.5d0 * (erf(constA  * (dble(j) / dble(ny+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
    end do

    ! Compute grid spacing using array slicing
    dy(1:ny+1) = yGrid(1:ny+1) - yGrid(0:ny)
    dx0 = dy(1)
    dt0 = dx0
    write(*,*) "---in 0 system---"
    write(*,*) "deltaX = ", dx0
    write(*,*) "deltaT = ", dt0
    viscosity0 =  Ur*(yGrid(ny+1)-yGrid(1))/Reynolds
    kappa0 = viscosity0/Pr
    write(*,*) "viscosity0 = ", real(viscosity0)
    write(*,*) "kappa0 = ", real(kappa0)

    ! Compute grid spacing in system LB
    length_LB = 1.0d0 / dy(1)
    dt = dt0 * length_LB
    yGrid(1:ny) = yGrid(1:ny)-dx0/2.0d0
    yGrid(ny+1) = yGrid(ny)+dx0/2.0d0
    xGrid=xGrid*(length_LB-1.0d0)*dble(B)
    yGrid=yGrid*length_LB
    ! Initialize flow variables
    rho = rho0
    temp = 0.0d0
    u = 0.0d0
    v = 0.0d0

    do j = 1, ny
        do i = 1, nx
            temp(i,j) = 1.0d0-yGrid(j)/(yGrid(ny)+yGrid(1))
            u(i,j) = 6.0d0*yGrid(j)/(yGrid(ny)+yGrid(1))*(1.0d0-yGrid(j)/(yGrid(ny)+yGrid(1)))*Ur
        enddo
    enddo

!    u_ave = 0.0d0
!    do j = 0, ny
!        u_ave = u_ave+(u(1,j)+u(1,j+1))*(yGrid(j+1)-yGrid(j))/2.0d0
!    end do
!    u_ave = u_ave/(length_LB-1.0d0)

    !lemda = Ur*dt/(xGrid(nx)-xGrid(nx-1))

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

    ! Output initial parameters for verification
    write(*,*) "nx =", nx, ", ny =", ny
    write(*,*) "Ra =", real(Ra)

    ! Initialize iteration count and error
    itc = 0

    ! Calculate viscosity based on LB unit
    viscosity_LB = Ur*(length_LB-1.0d0)/Reynolds
    kappa = viscosity_LB/Pr
    gbeta = Ra*viscosity_LB*kappa/((length_LB-1.0d0)**3)
    write(*,*) "---in LB unit---"
    write(*,*) "characteristic length   =", real(length_LB), "l.u."
    write(*,*) "viscosity_LB =", real(viscosity_LB), "l.u.^2/t.s."
    write(*,*) "timeStep ratio for (uniform) / (non-uniform) : ", real(length_LB / dble(nx))
    write(*,*) "    "

    ! Calculate relaxation time
    tauf = viscosity_LB * 3.0d0 + 0.5d0
    write(*,*) "tauf =", real(tauf)
    write(*,*)"_______________________________", 616250 *dt/((yGrid(ny)+yGrid(1))/Ur)
    write(*,*)"x(nxHalf)=",xGrid(nxHalf)
    ! Calculate MRT relaxation parameters
    Snu = 1.0d0/tauf
    Sq = 8.0d0*(2.0d0-Snu)/(8.0d0-Snu)
    sig_k = 1.0d0/(0.5d0+4.0d0*(1.0d0/Snu-0.5d0)/(3.0d0*Pr))

    timeUnit = dsqrt((length_LB-1.0d0) / gbeta)

    itc_max = 100000

    write(*,*)"itc_max=",itc_max

    do j = 1, ny
        do i = 1, nx
            velocitySquare = u(i,j)*u(i,j)+v(i,j)*v(i,j)
            do alpha = 0,8
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                f(i, j, alpha) = rho(i,j)*omega_U(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*velocitySquare)
            enddo

            do alpha=0,4
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                g(i,j,alpha) = omega_T(alpha)*temp(i,j)*(1.0d0+4.0d0*un(alpha))
            end do
        enddo
    enddo

    do j = 1, ny
        do alpha = 0, 8
            f_last(j,alpha) = f(nx,j,alpha)
        end do
    end do

    do i = 1, nx
        if(i == 1)then
            inter_x(i,:) = (/i+1, i, i+2/)

        elseif(i == nx)then
            inter_x(i,:) = (/i-1, i, i-2/)

        else
            inter_x(i,:) = (/i-1, i, i+1/)
        end if
    enddo

    do j = 1, ny
        if(j == 1)then
            inter_y(j,:) = (/j+1, j, j+2/)

        elseif(j == ny)then
            inter_y(j,:) = (/j-1, j, j-2/)

        else
            inter_y(j,:) = (/j-1, j, j+1/)
        end if
    enddo

    return
end subroutine initial

program main
    use commondata
    implicit none
    real(8) :: timestart, timeEnd
    integer(kind=4) :: i, num, ios
    integer, dimension(25) :: itc1
    call cpu_time(timestart)

    call initial()

!    open(unit=01,file='itc.dat',form="formatted", access="sequential", status='old', iostat=ios)
!
!    if (ios /= 0) then
!        print *, 'Error opening file!'
!        stop
!    end if
!
!    do i = 1, 25
!        read(01,*) num
!        itc1(i) = num
!    end do
!
!    close(01)
!    write(*,*)itc1
    !$acc data copy(u,v,rho,temp) copyin(xGrid,yGrid,ex,ey,f,g,inter_x,inter_y,f_last)&
    !$acc create(g_post,f_post,Fx,Fy)
    !do while((errorU > epsU).AND.(itc < itc_max))
    do while(itc < itc_max)
        itc = itc+1

        call collision_U()

        call collision_T()

        call interpolate()

        call bounceback_u()

        call bounceback_T()

        call macro_u()

        call macro_t()

        if((MOD(itc,10).EQ.0).and.itc > 3000) then
            call getT()
        end if

        if(itc == itc_max) then
        !if((MOD(itc,100000).EQ.0).and.itc > 300000) then
            call output_ASCII()
            !call getT_x5()
            !call getU_x5()
            !call getT_x10()
            !call getU_x10()
           ! write(*,*)itc
        end if
    enddo
    !$acc end data

    call cpu_time(timeEnd)
    write(*,*)"Time=", timeEnd-timestart, "s"
    write(*,*) "MLUPS = ", real(dble(nx*ny)/1e6*dble(itc)/(timeEnd-timeStart))
    !call output_ASCII()

    deallocate(u)
    deallocate(v)
    deallocate(rho)
    deallocate(f)
    deallocate(f_post)
    deallocate(g)
    deallocate(g_post)
    deallocate(temp)

    stop
end program main

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

!$acc parallel loop private(m,m_post,s,meq,fSource) gang vector collapse(2)
    do j=1,ny
        do i=1,nx

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
            Fy(i,j) = gbeta*(temp(i,j)-Tref)*rho(i,j)

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

            f_post(i,j,0) = ( m_post(0)-m_post(1)+m_post(2) )/9.0d0
            f_post(i,j,1) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+m_post(3)/6.0d0-m_post(4)/6.0d0 &
                            +m_post(7)*0.25d0
            f_post(i,j,2) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                            +m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)*0.25d0
            f_post(i,j,3) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0-m_post(3)/6.0d0+m_post(4)/6.0d0 &
                            +m_post(7)*0.25d0
            f_post(i,j,4) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                            -m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)*0.25d0
            f_post(i,j,5) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                            +m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)*0.25d0
            f_post(i,j,6) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                            +m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)*0.25d0
            f_post(i,j,7) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                            -m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)*0.25d0
            f_post(i,j,8) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                            -m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)*0.25d0
        enddo
    enddo
!$acc end parallel loop

!$acc parallel loop default(none) present(f,f_post)  gang vector collapse(2)
    do j = 1, ny
        do i = 1, nx
            f(i,j,0) = f_post(i,j,0)
        enddo
    enddo
!$acc end parallel loop

    return
end subroutine collision_U

subroutine collision_T()
    use commondata
    implicit none
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: n(0:4), n_post(0:4), neq(0:4)
    real(kind=8) :: Q(0:4)

!$acc parallel loop private(n,n_post,Q,neq) gang vector collapse(2)
    do j=1,ny
        do i=1,nx
            n(0) = g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
            n(1) = g(i,j,1)-g(i,j,3)
            n(2) = g(i,j,2)-g(i,j,4)
            n(3) = g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
            n(4) = g(i,j,1)-g(i,j,2)+g(i,j,3)-g(i,j,4)

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

            g_post(i,j,0) = n_post(0)-n_post(3)
            g_post(i,j,1) = n_post(1)/2.0d0+n_post(3)/4.0d0+n_post(4)/4.0d0
            g_post(i,j,2) = n_post(2)/2.0d0+n_post(3)/4.0d0-n_post(4)/4.0d0
            g_post(i,j,3) = -n_post(1)/2.0d0+n_post(3)/4.0d0+n_post(4)/4.0d0
            g_post(i,j,4) = -n_post(2)/2.0d0+n_post(3)/4.0d0-n_post(4)/4.0d0

        enddo
    enddo
!$acc end parallel loop

!$acc parallel loop default(none) present(g,g_post)  gang vector collapse(2)
    do j = 1, ny
        do i = 1, nx
            g(i,j,0) = g_post(i,j,0)
        enddo
    enddo
!$acc end parallel loop

    return
end subroutine collision_T

subroutine interpolate()
    use commondata
    implicit none
    real(kind=8) :: interpolateF, delta_x, delta_y
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: f0, f1, f2, g0, g1, g2
!$acc routine (interpolateF) seq
!$acc parallel loop present(f,f_post,ex,ey,xGrid,yGrid) gang vector collapse(2)
        do j = 1, ny
            do i = 1, nx
                do alpha = 1, 8
                    delta_x=dble(ex(alpha))*dt
                    delta_y=dble(ey(alpha))*dt

            f0 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), f_post(inter_x(i,1), inter_y(j,1), alpha), f_post(inter_x(i,1), inter_y(j,2)&
                , alpha), f_post(inter_x(i,1), inter_y(j,3), alpha))

            f1 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), f_post(inter_x(i,2), inter_y(j,1), alpha), f_post(inter_x(i,2), inter_y(j,2)&
                , alpha), f_post(inter_x(i,2), inter_y(j,3), alpha))

            f2 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), f_post(inter_x(i,3), inter_y(j,1), alpha), f_post(inter_x(i,3), inter_y(j,2)&
                , alpha), f_post(inter_x(i,3), inter_y(j,3), alpha))

            f(i, j, alpha) = interpolateF(xGrid(inter_x(i,1))+delta_x, xGrid(inter_x(i,2))+delta_x, &
                            xGrid(inter_x(i,3))+delta_x, xGrid(inter_x(i,2)), f0, f1, f2)

                end do
            enddo
        enddo
!$acc end parallel loop
!$acc parallel loop present(g,g_post,ex,ey,xGrid,yGrid) gang vector collapse(2)
        do j = 1, ny
            do i = 1, nx
                do alpha = 1, 4
                    delta_x=dble(ex(alpha))*dt
                    delta_y=dble(ey(alpha))*dt

            g0 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), g_post(inter_x(i,1), inter_y(j,1), alpha), g_post(inter_x(i,1), inter_y(j,2)&
                , alpha), g_post(inter_x(i,1), inter_y(j,3), alpha))

            g1 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), g_post(inter_x(i,2), inter_y(j,1), alpha), g_post(inter_x(i,2), inter_y(j,2)&
                , alpha), g_post(inter_x(i,2), inter_y(j,3), alpha))

            g2 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), g_post(inter_x(i,3), inter_y(j,1), alpha), g_post(inter_x(i,3), inter_y(j,2)&
                , alpha), g_post(inter_x(i,3), inter_y(j,3), alpha))

            g(i, j, alpha) = interpolateF(xGrid(inter_x(i,1))+delta_x, xGrid(inter_x(i,2))+delta_x, &
                            xGrid(inter_x(i,3))+delta_x, xGrid(inter_x(i,2)), g0, g1, g2)
                enddo
            end do
        end do
!$acc end parallel loop
end subroutine

!!NOTE: consider using compiler-specific directives to suggest inlining if necessary.
pure function interpolateF(x0, x1, x2, x, f0, f1, f2) result(f_interp)
    implicit none
    !$acc routine (interpolateF) seq
    real(kind=8), intent(in) :: x0, x1, x2, x, f0, f1, f2
    real(kind=8) :: f_interp

    ! Interpolation formula
    f_interp = ((x - x1) * (x - x2)) / ((x0 - x1) * (x0 - x2)) * f0 + &
               ((x - x0) * (x - x2)) / ((x1 - x0) * (x1 - x2)) * f1 + &
               ((x - x0) * (x - x1)) / ((x2 - x0) * (x2 - x1)) * f2

    return
end function interpolateF

subroutine bounceback_u()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    !$acc parallel loop present(f,f_post)
    do i=1, nx
        !Bottom side
        f(i,1,2) = f_post(i,1,4)
        f(i,1,5) = f_post(i,1,7)
        f(i,1,6) = f_post(i,1,8)

        !Top side
        f(i,ny,4) = f_post(i,ny,2)
        f(i,ny,7) = f_post(i,ny,5)
        f(i,ny,8) = f_post(i,ny,6)
    enddo
    !$acc end parallel

    !$acc parallel loop present(f,f_post)
    do j = 1, ny
        !Left side
        rho(1,j) = (f(1,j,0)+f(1,j,2)+f(1,j,4)+2.0d0*(f(1,j,3)+f(1,j,6)+f(1,j,7)))/(1.0d0-u(1,j))
        f(1,j,1) = f(1,j,3)+2.0d0*rho(1,j)*u(1,j)/3.0d0
        f(1,j,5) = f(1,j,7)-(f(1,j,2)-f(1,j,4))/2.0d0-dt*Fy(1,j)/4.0d0+rho(1,j)*u(1,j)/6.0d0
        f(1,j,8) = f(1,j,6)+(f(1,j,2)-f(1,j,4))/2.0d0+dt*Fy(1,j)/4.0d0+rho(1,j)*u(1,j)/6.0d0

        !Right side
        rho(nx-1,j) = f(nx-1,j,0)+f(nx-1,j,1)+f(nx-1,j,2)+f(nx-1,j,3)+f(nx-1,j,4)+f(nx-1,j,5)+f(nx-1,j,6)+f(nx-1,j,7)+f(nx-1,j,8)
        u(nx-1,j) = (f(nx-1,j,1)-f(nx-1,j,3)+f(nx-1,j,5)-f(nx-1,j,6)-f(nx-1,j,7)+f(nx-1,j,8))/rho(nx-1,j)
        v(nx-1,j) = (f(nx-1,j,2)-f(nx-1,j,4)+f(nx-1,j,5)+f(nx-1,j,6)-f(nx-1,j,7)-f(nx-1,j,8)+0.5d0*dt*Fy(nx-1,j))/rho(nx-1,j)
        lemda = dsqrt(*dt/(xGrid(nx)-xGrid(nx-1))
        f(nx,j,3) = (f_last(j,3)+lemda*f(nx-1,j,3))/(1.0d0+lemda)
        f(nx,j,6) = (f_last(j,6)+lemda*f(nx-1,j,6))/(1.0d0+lemda)
        f(nx,j,7) = (f_last(j,7)+lemda*f(nx-1,j,7))/(1.0d0+lemda)

        f_last(j,3) = f(nx,j,3)
        f_last(j,6) = f(nx,j,6)
        f_last(j,7) = f(nx,j,7)
    enddo
    !$acc end parallel

    return
end subroutine

subroutine bounceback_T()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: geq
    !$acc routine (geq) seq
    !$acc parallel loop present(g,g_post)
    do i=1, nx
        !top
        g(i,ny,4) = -g_post(i,ny,2) + Tcold / 4.0d0
        !bottom
        g(i,1,2) = -g_post(i,1,4) + Thot / 4.0d0
    end do
    !$acc end parallel

    !$acc parallel loop present(g,g_post)
        do j = 1, ny
            !left
            temp(2,j) = g(2,j,0)+g(2,j,1)+g(2,j,2)+g(2,j,3)+g(2,j,4)
            rho(2,j) = f(2,j,0)+f(2,j,1)+f(2,j,2)+f(2,j,3)+f(2,j,4)+f(2,j,5)+f(2,j,6)+f(2,j,7)+f(2,j,8)
            u(2,j) = (f(2,j,1)-f(2,j,3)+f(2,j,5)-f(2,j,6)-f(2,j,7)+f(2,j,8)+0.5d0*dt*Fx(2,j))/rho(2,j)
            v(2,j) = (f(2,j,2)-f(2,j,4)+f(2,j,5)+f(2,j,6)-f(2,j,7)-f(2,j,8)+0.5d0*dt*Fy(2,j))/rho(2,j)

            g(1,j,1) = geq(1,temp(1,j),u(1,j),v(1,j))+g(2,j,1)-geq(1,temp(2,j),u(2,j),v(2,j))
            !right
            g(nx,j,3) = g(nx-1,j,3)
        end do
    !$acc end parallel
return
end subroutine bounceback_T

subroutine macro_u()
    use commondata
    implicit none
    integer(kind=4) :: i, j

!$acc parallel loop default(none) present(rho,u,v,f,Fx,Fy) gang vector collapse(2)
    do j=1, ny
        do i=3, nx
            rho(i,j) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
            u(i,j) = (f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)+0.5d0*dt*Fx(i,j))/rho(i,j)
            v(i,j) = (f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)+0.5d0*dt*Fy(i,j))/rho(i,j)
        enddo
    enddo
!$acc end parallel loop
    return
end subroutine macro_u

subroutine macro_t()
    use commondata
    implicit none
    integer(kind=4) :: i, j

!$acc parallel loop default(none) present(temp,g) gang vector collapse(2)
    do j=1, ny
        do i=3, nx
            temp(i,j) = g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
        end do
    end do
!$acc end parallel loop
    return
end subroutine macro_t

function feq(alpha,rho,u,v)
    !$acc routine (feq) seq
    implicit none
    real(8) :: feq
    integer :: alpha
    real(8) :: rho, u, v
    real(8) :: us2, un
    real(8) :: omega
    integer :: ex(0:8), ey(0:8)
    data ex/0, 1, 0, -1,  0, 1, -1, -1,  1/
    data ey/0, 0, 1,  0, -1, 1,  1, -1, -1/

    if (alpha.EQ.0) then
        omega = 4.0d0/9.0d0
    elseif( (alpha.GE.1).AND.(alpha.LE.4) ) then
        omega = 1.0d0/9.0d0
    elseif( (alpha.GE.5).AND.(alpha.LE.8) ) then
        omega = 1.0d0/36.0d0
    else
        write(*,*) "error in function feq"
        stop
    endif

    us2 = u*u+v*v
    un = u*ex(alpha)+v*ey(alpha)
    feq = omega*rho*(1.0d0+3.0d0*un+4.5d0*un*un-1.5d0*us2)

    return
end function feq

function geq(alpha,temp,u,v)
    !$acc routine (geq) seq
    implicit none
    real(8) :: geq
    integer :: alpha
    real(8) :: temp, u, v
    real(8) :: un
    real(8) :: omega
    integer :: ex(0:4), ey(0:4)
    data ex/0, 1, 0, -1,  0/
    data ey/0, 0, 1,  0, -1/

    if (alpha.EQ.0) then
        omega = 1.0d0/2.0d0
    elseif( (alpha.GE.1).AND.(alpha.LE.4) ) then
        omega = 1.0d0/8.0d0
    else
        write(*,*) "error in function geq"
        stop
    endif

    un = u*ex(alpha)+v*ey(alpha)
    geq = omega*temp*(1.0d0+4.0d0*un)
    return
end function geq

subroutine output_ASCII()
    use commondata
    implicit none
    integer :: i, j
    character(len=100) :: filename

    write(filename,*) itc
    filename = adjustl(filename)

    !$acc update self(u, v, temp, rho)
    open(unit=01,file='MRTcavity-'//trim(filename)//'.dat',status='unknown')

    write(01,*) 'TITLE="thermal convective flows"'
    write(01,*) 'VARIABLES="X" "Y" "U" "V" "T" "P"'
    write(01,101) nx, ny
    do j=1,ny
        do i=1,nx
            write(01,100) xGrid(i), yGrid(j), u(i,j), v(i,j), temp(i,j), rho(i,j)/3.0d0
        enddo
    enddo
100 format(1x,2(e11.4,' '),10(e13.6,' '))
101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
    close(01)

    return
end subroutine output_ASCII

subroutine getT()
    use commondata
    implicit none
    integer(kind=4) :: localx = (nxHalf-1)/2+1
    !$acc update self(temp)
    open(unit=02, file='T.dat', status='unknown', position='append', action='write')
    write(02, *) itc, temp(localx,nyHalf)
    close(02)

    return
end subroutine getT

subroutine getT_x5()
    use commondata
    implicit none
    integer(kind=4) :: j
    integer(kind=4) :: localx = (nxHalf-1)/2+1
    character(len=100) :: filename

    write(filename,*) itc
    filename = adjustl(filename)
    !$acc update self(temp)
    open(unit=03, file='T_x5-'//trim(filename)//'.dat', status='unknown', action='write')
        do j = 1, ny
            write(03, *) yGrid(j)/(yGrid(ny)+yGrid(1)), temp(localx,j)
        end do
    close(03)

    return
end subroutine getT_x5

subroutine getU_x5()
    use commondata
    implicit none
    integer(kind=4) :: j
    integer(kind=4) :: localx = (nxHalf-1)/2+1
    character(len=100) :: filename

    write(filename,*) itc
    filename = adjustl(filename)
    !$acc update self(u,v)
    open(unit=04, file='U_x5-'//trim(filename)//'.dat', status='unknown', action='write')
        do j = 1, ny
            write(04, *) yGrid(j)/(yGrid(ny)+yGrid(1)), u(localx,j)/Ur
        end do
    close(04)

    open(unit=05, file='v_x5-'//trim(filename)//'.dat', status='unknown', action='write')
        do j = 1, ny
            write(05, *) yGrid(j)/(yGrid(ny)+yGrid(1)), v(localx,j)/Ur
        end do
    close(05)

    return
end subroutine getU_x5

subroutine getT_x10()
    use commondata
    implicit none
    integer(kind=4) :: j
    character(len=100) :: filename

    write(filename,*) itc
    filename = adjustl(filename)

    !$acc update self(temp)
    open(unit=06,file='T_x10-'//trim(filename)//'.dat',status='unknown', action='write')
        do j = 1, ny
            write(06, *) yGrid(j)/(yGrid(ny)+yGrid(1)), temp(nxHalf,j)
        end do
    close(06)

    return
end subroutine getT_x10

subroutine getU_x10()
    use commondata
    implicit none
    integer(kind=4) :: j
    character(len=100) :: filename

    write(filename,*) itc
    filename = adjustl(filename)
    !$acc update self(u,v)
    open(unit=07, file='U_x10-'//trim(filename)//'.dat', status='unknown', action='write')
        do j = 1, ny
            write(07, *) yGrid(j)/(yGrid(ny)+yGrid(1)), u(nxHalf,j)/Ur
        end do
    close(07)

    open(unit=08, file='v_x10'//trim(filename)//'.dat', status='unknown', action='write')
        do j = 1, ny
            write(08, *) yGrid(j)/(yGrid(ny)+yGrid(1)), v(nxHalf,j)/Ur
        end do
    close(08)

    return
end subroutine getU_x10
