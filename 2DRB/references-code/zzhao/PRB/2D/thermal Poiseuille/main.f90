module commondata
        implicit none

        ! Grid parametersi
        integer(kind=4), parameter :: ny=33,nx=2*ny+1
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1

        real(kind=8), parameter :: Ra=1000          !Rayleigh number
        real(kind=8), parameter :: Reynolds=15      !Reynolds number
        real(kind=8), parameter :: Pr=0.71d0        !Prandtl number
        real(kind=8), parameter :: Thot=1.0d0, Tcold=0.0d0, Tref=0.0d0
        real(kind=8) :: U0, vIn, vOut
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
        integer(kind=4), parameter :: itc_max = 50000000 ! Maximum iterations

        ! Convergence criteria
        real(kind=8) :: errorU, errorT              ! Current error
        real(kind=8), parameter :: epsU=1e-9, epsT=1e-9  ! Convergence threshold

        ! Grid coordinates
        real(kind=8) :: xGrid(nx), yGrid(0:ny+1)

        ! Flow variables
        real(kind=8), allocatable :: u(:,:), v(:,:)      ! Velocity components
        real(kind=8), allocatable :: temp(:,:)           ! temperature
        real(kind=8), allocatable :: rho(:,:)            ! Density
        real(kind=8), allocatable :: up(:,:), vp(:,:)    ! Previous velocity components for error checking
        real(kind=8), allocatable :: utemp(:,:)          ! Previous temperature for error checking
        real(kind=8), allocatable :: Fx(:,:), Fy(:,:)    ! force components

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

        integer(kind=4) :: inter_x(nx,3), inter_y(ny,3)

end module commondata

subroutine initial()
    use commondata
    implicit none
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: un(0:8)
    real(kind=8) :: velocitySquare
    real(kind=8) :: dy(ny+1)
    real(kind=8) :: constA

    ! Output initial parameters for verification
    write(*,*) "nx =", nx, ", ny =", ny
    write(*,*) "Ra =", real(Ra)

    ! Initialize iteration count and error
    itc = 0
    errorU = 100.0d0
    errorT = 100.0d0

    ! Compute grid coordinates
    constA = 3.2d0
    do i = 1, nx
        xGrid(i) = (dble(i) - 1.0d0) / (dble(nx)-1.0d0)
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

    Snu = 0.8d0
    Sq = 8.0d0*(2.0d0-Snu)/(8.0d0-Snu)
    tauf = 1.0d0/Snu
    sig_k = 1.0d0/(0.5d0+4.0d0*(1.0d0/Snu-0.5d0)/(3.0d0*Pr))
    ! Calculate viscosity based on LB unit
    length_LB = 1.0d0 / dy(1)
    dt = dt0 * length_LB
    !viscosity_LB = U0*(length_LB-1.0d0)/Reynolds
    viscosity_LB = (tauf-0.5d0)/3.0d0
    U0 = viscosity_LB*Reynolds/(length_LB-1.0d0)
    vIn = U0
    vOut = U0
    kappa = viscosity_LB/Pr
    gbeta = Ra*viscosity_LB*kappa/((length_LB-1.0d0)**3)
    write(*,*) "---in LB unit---"
    write(*,*) "U0=", U0
    write(*,*) "characteristic length   =", real(length_LB), "l.u."
    write(*,*) "viscosity_LB =", real(viscosity_LB), "l.u.^2/t.s."
    write(*,*) "timeStep ratio for (uniform) / (non-uniform) : ", real(length_LB / dble(nx))
    write(*,*)  gbeta
    write(*,*) "    "

    ! Compute grid spacing in system LB
    yGrid = yGrid-dy(1)/2.0d0
    yGrid = yGrid*length_LB
    xGrid = 2.0d0*xGrid*(length_LB-1.0d0)
    write(*,*)"x=",xGrid(nx), "y=",yGrid(ny)+yGrid(1)
    ! Calculate relaxation time
    !tauf = viscosity_LB * 3.0d0 + 0.5d0
    write(*,*) "tauf =", real(tauf)

    ! Calculate MRT relaxation parameters
!    Snu = 1.0d0/tauf

    ! Allocate flow variables
    allocate (u(nx,ny))
    allocate (v(nx,ny))
    allocate (rho(nx,ny))
    allocate (up(nx,ny))
    allocate (vp(nx,ny))
    allocate (temp(nx,ny))
    allocate (utemp(nx,ny))
    allocate (Fx(nx,ny))
    allocate (Fy(nx,ny))

    allocate (f(nx,ny,0:8))
    allocate (f_post(nx,ny,0:8))
    allocate (g(nx,ny,0:4))
    allocate (g_post(nx,ny,0:4))

    ! Initialize flow variables
    rho = rho0
    temp = 0.0d0
    utemp=0.0d0
    u = 0.0d0
    v = 0.0d0
    up = 0.0d0
    vp = 0.0d0

    do j=1,ny
        do i=1,nx
            temp(i,j) = dble(j-1)/dble(ny-1)*(Thot-Tcold)+Thot
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

    do j = 1, ny
        do i = 1, nx
            velocitySquare = u(i,j)*u(i,j)+v(i,j)*v(i,j)
            do alpha = 0,8
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                f(i, j, alpha) = rho(i,j)*omega_U(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*velocitySquare)
            enddo

            do alpha=0,4
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                g(i,j,alpha)=omega_T(alpha)*temp(i,j)*(1.0d0+4.0d0*un(alpha))
            end do
        enddo
    enddo

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

    call cpu_time(timestart)

    call initial()

    !$acc data copy(u,v,rho,temp) copyin(xGrid,yGrid,ex,ey,f,g,inter_x,inter_y)&
    !$acc create(g_post,f_post,up,vp,utemp,Fx,Fy)
    do while(((errorU > epsU).or.(errorT > epsT)).AND.(itc < itc_max))

        itc = itc+1

        call collision_U()

        call collision_T()

        call interpolate()

        call bounceback_u()

        call macro_u()

        call bounceback_T()

        call macro_t()

        if(MOD(itc,2000).EQ.0) then
            call check()
        endif

    enddo
    !$acc end data

    call cpu_time(timeEnd)
    write(*,*)"Time=", timeEnd-timestart, "s"
    write(*,*) "MLUPS = ", real(dble(nx*ny)/1e6*dble(itc)/(timeEnd-timeStart))
    call output_ASCII()

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

    !$acc parallel loop default(none) present(f,f_post)
    do j=1, ny
        !Left side
        f(1,j,1) = f(nx,j,1)
        f(1,j,5) = f(nx,j,5)
        f(1,j,8) = f(nx,j,8)

        !Right side
        f(nx,j,3) = f(1,j,3)
        f(nx,j,6) = f(1,j,6)
        f(nx,j,7) = f(1,j,7)
    enddo
    !$acc end parallel
    !$acc parallel loop default(none) present(rho,f,f_post)
    do i=1, nx
        !Bottom side
        f(i,1,2) = f_post(i,1,4)+2.0d0*rho(i,1)*vIn/3.0d0
        f(i,1,5) = f_post(i,1,7)+rho(i,1)*vIn/6.0d0
        f(i,1,6) = f_post(i,1,8)+rho(i,1)*vIn/6.0d0

        !Top side
        f(i,ny,4) = f_post(i,ny,2)-2.0d0*rho(i,ny)*vOut/3.0d0
        f(i,ny,7) = f_post(i,ny,5)-rho(i,ny)*(U0+vOut)/6.0d0
        f(i,ny,8) = f_post(i,ny,6)-rho(i,ny)*(-U0+vOut)/6.0d0
    enddo
    !$acc end parallel
return
end subroutine bounceback_u

subroutine bounceback_T()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: geq
    !$acc routine (geq) seq
!$acc parallel loop present(g,g_post)
    do j=1, ny
        !left
        g(1,j,1) = g(nx,j,1)
        !!right
        g(nx,j,3) = g(1,j,3)
    end do
!$acc end parallel
!$acc parallel loop present(g,g_post)
    do i=1, nx
        !top
        g(i,ny,4) = -g_post(i,ny,2)+2.0d0*geq(2,Thot,U0,0.0d0)

        !bottom
        g(i,1,2) = -g_post(i,1,4)+2.0d0*geq(4,Tcold,0.0d0,0.0d0)

    end do
!$acc end parallel
    return
end subroutine

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

subroutine macro_u()
    use commondata
    implicit none
    integer(kind=4) :: i, j

!$acc parallel loop default(none) present(rho,u,v,f,Fx,Fy) gang vector collapse(2)
    do j=1, ny
        do i=1, nx
            rho(i,j) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
            u(i,j) = (f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)+0.5d0*dt*Fx(i,j) )/rho(i,j)
            v(i,j) = (f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)+0.5d0*dt*Fy(i,j) )/rho(i,j)
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
        do i=1, nx
            temp(i,j) = g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
        end do
    end do
!$acc end parallel loop
    return
end subroutine macro_t

subroutine check()
    use commondata
    implicit none
    integer :: i, j
    real(kind=8) :: error1, error2,error3, error4
    error1 = 0.0d0
    error2 = 0.0d0
    error3 = 0.0d0
    error4 = 0.0d0

 !$acc parallel loop default(none) reduction(+:error1,error2) present(u,v,up,vp)
    do j=1,ny
        do i=1,nx
            error1 = error1+(u(i,j)-up(i,j))**2.0d0+(v(i,j)-vp(i,j))**2.0d0
            error2 = error2+u(i,j)**2.0d0+v(i,j)**2.0d0
            up(i,j) = u(i,j)
            vp(i,j) = v(i,j)
        enddo
    enddo
!$acc end parallel loop

    errorU = dsqrt(error1)/dsqrt(error2)

!$acc parallel loop default(none) reduction(+:error3,error4) present(temp,utemp)
    do j=1,ny
        do i=1,nx
            error3 = error3+(temp(i,j)-utemp(i,j))**2
            error4 = error4+temp(i,j)*temp(i,j)

            utemp(i,j) = temp(i,j)
        end do
    end do
!$acc end parallel
    errorT = dsqrt(error3)/dsqrt(error4)


    write(*,*) itc,' ',errorU,' ',errorT

    return
end subroutine check

subroutine output_ASCII()
    use commondata
    implicit none
    integer :: i, j
    character(len=100) :: filename

    write(filename,*) ny
    filename = adjustl(filename)


    open(unit=02,file='MRTcavity-'//trim(filename)//'.dat',status='unknown')

    write(02,*) 'TITLE="thermal convective flows"'
    write(02,*) 'VARIABLES="X" "Y" "U" "V" "T" "rho"'
    write(02,101) nx, ny
    do j=1,ny
        do i=1,nx
            write(02,100) xGrid(i), yGrid(j), u(i,j), v(i,j), temp(i,j), rho(i,j)
        enddo
    enddo
100 format(1x,2(e11.4,' '),10(e13.6,' '))
101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
    close(02)

    return
end subroutine output_ASCII
