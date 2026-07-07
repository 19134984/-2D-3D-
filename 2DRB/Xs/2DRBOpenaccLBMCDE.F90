!=============================================================
!!!    2D side-heated natural convection
!!!    LBM-CDE reproduction draft based on pdf/LBM-CDE.pdf
!!!    D2Q9 BGK for both flow and temperature fields
!!!    OpenACC single-GPU port of 2DRBOpenmpLBMCDE.F90
!=============================================================


!=============================================================
!   Case switches
#define steadyFlow
!#define unsteadyFlow

#define HorizontalWallsNoslip
#define VerticalWallsNoslip

#define SideHeatedCell
#define HorizontalWallsAdiabatic
#define VerticalWallsConstT

!   The 2D cases in LBM-CDE.pdf use D2Q9 for both velocity and scalar fields.
#define TemperatureD2Q9

#ifndef ITC_MAX_OVERRIDE
#define ITC_MAX_OVERRIDE 20000000
#endif
#ifndef RAYLEIGH_OVERRIDE
#define RAYLEIGH_OVERRIDE 1.0d5
#endif
!=============================================================


!=============================================================
!   Global data
    module commondata
        implicit none

        integer(kind=4), parameter :: nx=128, ny=128
        real(kind=8), parameter :: rho0=1.0d0
        real(kind=8), parameter :: cs2=1.0d0/3.0d0
        real(kind=8), parameter :: pi=acos(-1.0d0)

#ifdef SideHeatedCell
        real(kind=8), parameter :: lengthUnit=dble(nx)
#else
        real(kind=8), parameter :: lengthUnit=dble(ny)
#endif

        real(kind=8), parameter :: Rayleigh=RAYLEIGH_OVERRIDE
        real(kind=8), parameter :: Prandtl=0.71d0
        real(kind=8), parameter :: Mach=0.1d0
        real(kind=8), parameter :: Thot=0.5d0, Tcold=-0.5d0
        real(kind=8), parameter :: Tref=0.5d0*(Thot+Tcold)
        real(kind=8), parameter :: deltaT=Thot-Tcold

        ! LBM-CDE free parameters.  chi=0 recovers the default BGK diffusivity/viscosity.
        real(kind=8), parameter :: chi=0.8d0
        real(kind=8), parameter :: chi_s=chi
        real(kind=8), parameter :: chi_b=0.0d0
        real(kind=8), parameter :: chi_kappa=chi
        real(kind=8), parameter :: heatSourceQ=0.0d0

        ! Target physical transport coefficients in lattice units.
        real(kind=8), parameter :: viscosity=Mach*lengthUnit*dsqrt(Prandtl/(3.0d0*Rayleigh))
        real(kind=8), parameter :: diffusivity=viscosity/Prandtl

        ! LBM-CDE relaxation times: nu=(tau_fL-0.5)*(1-chi_s)*cs2.
        real(kind=8), parameter :: taufL=0.5d0+viscosity/(cs2*(1.0d0-chi_s))
        real(kind=8), parameter :: taugL=0.5d0+diffusivity/(cs2*(1.0d0-chi_kappa))
        real(kind=8), parameter :: bulkViscosity=(taufL-0.5d0)*(1.0d0-chi_b)*cs2
        real(kind=8), parameter :: omegaF=1.0d0/taufL
        real(kind=8), parameter :: omegaG=1.0d0/taugL

        ! Boussinesq force in the weakly compressible LBM-CDE momentum equation.
        real(kind=8), parameter :: gBeta=Rayleigh*viscosity*diffusivity/(deltaT*lengthUnit**3)
        real(kind=8), parameter :: timeUnit=dsqrt(lengthUnit/(gBeta*deltaT))
        real(kind=8), parameter :: velocityUnit=dsqrt(gBeta*deltaT*lengthUnit)

#ifdef steadyFlow
        real(kind=8), parameter :: epsU=1.0d-7, epsT=1.0d-7
        integer(kind=4), parameter :: itc_max=ITC_MAX_OVERRIDE
        integer(kind=4), parameter :: checkInterval=1000
#endif
#ifdef unsteadyFlow
        integer(kind=4), parameter :: itc_max=ITC_MAX_OVERRIDE
        integer(kind=4), parameter :: checkInterval=1000
#endif

        integer(kind=4), parameter :: ex(0:8)=(/0, 1, 0, -1, 0, 1, -1, -1, 1/)
        integer(kind=4), parameter :: ey(0:8)=(/0, 0, 1, 0, -1, 1, 1, -1, -1/)
        integer(kind=4), parameter :: opp(0:8)=(/0, 3, 4, 1, 2, 7, 8, 5, 6/)

        real(kind=8), parameter :: omega(0:8)=(/ &
            4.0d0/9.0d0, &
            1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, &
            1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0 /)

        ! lambdaT supplies the Tp/(rho0*cs2) term in Eq. (16):
        ! sum(lambdaT)=0 and sum(lambdaT*e_i_alpha*e_i_beta)=cs2*delta_alpha_beta.
        real(kind=8), parameter :: lambdaT(0:8)=(/ &
            -5.0d0/9.0d0, &
            1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, &
            1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0 /)

        real(kind=8), allocatable :: u(:,:), v(:,:), T(:,:), rho(:,:), pressure(:,:)
        real(kind=8), allocatable :: up(:,:), vp(:,:), Tp(:,:)
        real(kind=8), allocatable :: Fx(:,:), Fy(:,:)
        real(kind=8), allocatable :: gradTx(:,:), gradTy(:,:)
        real(kind=8), allocatable :: Sxx(:,:), Sxy(:,:), Syy(:,:), Sdiv(:,:)
        real(kind=8), allocatable :: f(:,:,:), f_post(:,:,:)
        real(kind=8), allocatable :: g(:,:,:), g_post(:,:,:)
        real(kind=8) :: xp(0:nx+1), yp(0:ny+1)

        real(kind=8) :: errorU, errorT
        real(kind=8) :: Nu_global, Nu_hot, Nu_cold, Nu_middle
        real(kind=8) :: ReVolAvg, umax_mid, vmax_mid, Numax_mid
        real(kind=8) :: y_umax_mid, x_vmax_mid, y_Numax_mid
        integer(kind=4) :: itc

        character(len=100) :: settingsFile="SimulationSettings2DOpenaccLBMCDE.txt"
        character(len=100) :: pltFile="buoyancyCavity2DOpenaccLBMCDE.plt"
        character(len=100) :: historyFile="NuRe_2DOpenaccLBMCDE.dat"

    contains

        pure integer(kind=4) function qex(alpha)
            !$acc routine seq
            implicit none
            integer(kind=4), intent(in) :: alpha

            select case(alpha)
            case(1,5,8)
                qex = 1
            case(3,6,7)
                qex = -1
            case default
                qex = 0
            end select
        end function qex

        pure integer(kind=4) function qey(alpha)
            !$acc routine seq
            implicit none
            integer(kind=4), intent(in) :: alpha

            select case(alpha)
            case(2,5,6)
                qey = 1
            case(4,7,8)
                qey = -1
            case default
                qey = 0
            end select
        end function qey

        pure real(kind=8) function qomega(alpha)
            !$acc routine seq
            implicit none
            integer(kind=4), intent(in) :: alpha

            select case(alpha)
            case(0)
                qomega = 4.0d0/9.0d0
            case(1,2,3,4)
                qomega = 1.0d0/9.0d0
            case default
                qomega = 1.0d0/36.0d0
            end select
        end function qomega

        pure real(kind=8) function qlambdaT(alpha)
            !$acc routine seq
            implicit none
            integer(kind=4), intent(in) :: alpha

            select case(alpha)
            case(0)
                qlambdaT = -5.0d0/9.0d0
            case(1,2,3,4)
                qlambdaT = 1.0d0/9.0d0
            case default
                qlambdaT = 1.0d0/36.0d0
            end select
        end function qlambdaT

        pure real(kind=8) function calc_feq(alpha, rhoLoc, uLoc, vLoc)
            !$acc routine seq
            implicit none
            integer(kind=4), intent(in) :: alpha
            real(kind=8), intent(in) :: rhoLoc, uLoc, vLoc
            real(kind=8) :: eu, uu, deltaRhoLoc

            eu = dble(qex(alpha))*uLoc+dble(qey(alpha))*vLoc
            uu = uLoc*uLoc+vLoc*vLoc
            deltaRhoLoc = rhoLoc-rho0
            calc_feq = qomega(alpha)*(deltaRhoLoc + &
                rho0*(eu/cs2+0.5d0*eu*eu/(cs2*cs2)-0.5d0*uu/cs2))
        end function calc_feq

        pure real(kind=8) function calc_geq(alpha, TLoc, uLoc, vLoc, rhoLoc)
            !$acc routine seq
            implicit none
            integer(kind=4), intent(in) :: alpha
            real(kind=8), intent(in) :: TLoc, uLoc, vLoc, rhoLoc
            real(kind=8) :: eu, uu, pressureLoc

            eu = dble(qex(alpha))*uLoc+dble(qey(alpha))*vLoc
            uu = uLoc*uLoc+vLoc*vLoc
            pressureLoc = cs2*(rhoLoc-rho0)
            calc_geq = qomega(alpha)*TLoc*(1.0d0+eu/cs2+0.5d0*eu*eu/(cs2*cs2)-0.5d0*uu/cs2) + &
                qlambdaT(alpha)*TLoc*pressureLoc/(rho0*cs2)
        end function calc_geq

    end module commondata
!=============================================================


!=============================================================
    program main
        use openacc
        use commondata
        implicit none
        real(kind=8) :: timeStart, timeEnd, wallStart, wallEnd
        integer(kind=8) :: wallClockStart, wallClockEnd, wallClockRate

        call acc_init(acc_device_default)
        call initial()
        call enter_data_2d_openacc()

        call cpu_time(timeStart)
        call system_clock(wallClockStart, wallClockRate)

#ifdef steadyFlow
        do while(((errorU.GT.epsU).OR.(errorT.GT.epsT)).AND.(itc.LT.itc_max))
#endif
#ifdef unsteadyFlow
        do while(itc.LT.itc_max)
#endif
            itc = itc+1

            call compute_force()
            call compute_strain_rate()
            call collision()
            call streaming()
            call bounceback()
            call macro()

            call compute_force()
            call compute_temperature_gradient()
            call collisionT()
            call streamingT()
            call bouncebackT()
            call macroT()

            if(mod(itc, checkInterval).EQ.0) then
                call compute_force()
                call compute_temperature_gradient()
                call update_host_monitor_2d_openacc()
                call check()
                call calNuRe()
                call append_history()
            endif
        enddo

        call compute_force()
        call compute_temperature_gradient()
        call update_host_monitor_2d_openacc()
        call calNuRe()
        call output_Tecplot()

        call cpu_time(timeEnd)
        call system_clock(wallClockEnd)
        wallStart = 0.0d0
        wallEnd = dble(wallClockEnd-wallClockStart)/dble(wallClockRate)
        call final_log(timeEnd-timeStart, wallEnd-wallStart)
        call exit_data_2d_openacc()
        call finalize_arrays()
    end program main
!=============================================================


!=============================================================
    subroutine initial()
        use commondata
        implicit none
        integer(kind=4) :: i, j, alpha

        itc = 0
        errorU = 100.0d0
        errorT = 100.0d0
        Nu_global = 0.0d0
        Nu_hot = 0.0d0
        Nu_cold = 0.0d0
        Nu_middle = 0.0d0
        ReVolAvg = 0.0d0
        umax_mid = 0.0d0
        vmax_mid = 0.0d0
        Numax_mid = 0.0d0
        y_umax_mid = 0.0d0
        x_vmax_mid = 0.0d0
        y_Numax_mid = 0.0d0

        allocate(u(nx,ny), v(nx,ny), T(nx,ny), rho(nx,ny), pressure(nx,ny))
        allocate(up(nx,ny), vp(nx,ny), Tp(nx,ny))
        allocate(Fx(nx,ny), Fy(nx,ny))
        allocate(gradTx(nx,ny), gradTy(nx,ny))
        allocate(Sxx(nx,ny), Sxy(nx,ny), Syy(nx,ny), Sdiv(nx,ny))
        allocate(f(0:8,nx,ny), f_post(0:8,0:nx+1,0:ny+1))
        allocate(g(0:8,nx,ny), g_post(0:8,0:nx+1,0:ny+1))

        do i = 0, nx+1
            xp(i) = dble(i)-0.5d0
        enddo
        xp(0) = 0.0d0
        xp(nx+1) = dble(nx)
        xp = xp/lengthUnit

        do j = 0, ny+1
            yp(j) = dble(j)-0.5d0
        enddo
        yp(0) = 0.0d0
        yp(ny+1) = dble(ny)
        yp = yp/lengthUnit

        rho = rho0
        u = 0.0d0
        v = 0.0d0
        pressure = 0.0d0
        Fx = 0.0d0
        Fy = 0.0d0
        gradTx = 0.0d0
        gradTy = 0.0d0
        Sxx = 0.0d0
        Sxy = 0.0d0
        Syy = 0.0d0
        Sdiv = 0.0d0

        do j = 1, ny
            do i = 1, nx
                T(i,j) = Thot + (xp(i)-xp(0))/(xp(nx+1)-xp(0))*(Tcold-Thot)
            enddo
        enddo

        f = 0.0d0
        g = 0.0d0
        f_post = 0.0d0
        g_post = 0.0d0
        do j = 1, ny
            do i = 1, nx
                do alpha = 0, 8
                    f(alpha,i,j) = calc_feq(alpha, rho(i,j), u(i,j), v(i,j))
                    g(alpha,i,j) = calc_geq(alpha, T(i,j), u(i,j), v(i,j), rho(i,j))
                enddo
            enddo
        enddo

        up = u
        vp = v
        Tp = T

        open(unit=00, file=trim(settingsFile), status='replace', action='write')
        write(00,*) "2D OpenACC LBM-CDE reproduction"
        write(00,*) "Mesh =", nx, ny
        write(00,*) "Device model = OpenACC single GPU"
        write(00,*) "Case = SideHeatedCell"
        write(00,*) "Flow lattice = D2Q9 BGK"
        write(00,*) "Temperature lattice = D2Q9 BGK"
        write(00,*) "Rayleigh =", real(Rayleigh,kind=8), "Prandtl =", real(Prandtl,kind=8)
        write(00,*) "Mach =", real(Mach,kind=8), "Length unit =", real(lengthUnit,kind=8)
        write(00,*) "chi_s =", real(chi_s,kind=8), "chi_kappa =", real(chi_kappa,kind=8), "chi_b =", real(chi_b,kind=8)
        write(00,*) "viscosity =", real(viscosity,kind=8), "diffusivity =", real(diffusivity,kind=8)
        write(00,*) "taufL =", real(taufL,kind=8), "taugL =", real(taugL,kind=8)
        write(00,*) "gBeta =", real(gBeta,kind=8)
        write(00,*) "timeUnit =", real(timeUnit,kind=8), "velocityUnit =", real(velocityUnit,kind=8)
        write(00,*) "itc_max =", itc_max, "checkInterval =", checkInterval
        close(00)

        open(unit=01, file=trim(historyFile), status='replace', action='write')
        write(01,'(A)') "# itc time_ff errorU errorT Nu_hot Nu_cold Nu_middle ReVolAvg"
        close(01)
    end subroutine initial
!=============================================================


!=============================================================
    subroutine enter_data_2d_openacc()
        use commondata
        implicit none

        !$acc enter data copyin(u,v,T,rho,pressure,Fx,Fy,gradTx,gradTy)
        !$acc enter data copyin(Sxx,Sxy,Syy,Sdiv,f,f_post,g,g_post)
    end subroutine enter_data_2d_openacc
!=============================================================


!=============================================================
    subroutine update_host_monitor_2d_openacc()
        use commondata
        implicit none

        !$acc update self(u,v,T,rho,pressure,gradTx,gradTy)
    end subroutine update_host_monitor_2d_openacc
!=============================================================


!=============================================================
    subroutine exit_data_2d_openacc()
        use commondata
        implicit none

        !$acc exit data delete(u,v,T,rho,pressure,Fx,Fy,gradTx,gradTy)
        !$acc exit data delete(Sxx,Sxy,Syy,Sdiv,f,f_post,g,g_post)
    end subroutine exit_data_2d_openacc
!=============================================================


!=============================================================
    subroutine compute_force()
        use commondata
        implicit none
        integer(kind=4) :: i, j

        !$acc parallel loop collapse(2) present(Fx,Fy,T,rho,pressure) private(i,j)
        do j = 1, ny
            do i = 1, nx
                pressure(i,j) = cs2*(rho(i,j)-rho0)
                Fx(i,j) = 0.0d0
                Fy(i,j) = rho0*gBeta*(T(i,j)-Tref)
            enddo
        enddo
    end subroutine compute_force
!=============================================================


!=============================================================
    subroutine compute_strain_rate()
        use commondata
        implicit none
        integer(kind=4) :: i, j, alpha
        real(kind=8) :: feqLoc, neqxx, neqxy, neqyy
        real(kind=8) :: denomDiag, denomShear, denomDiv, coeffTrace

        denomDiag = 2.0d0*viscosity + rho0*cs2
        denomShear = 4.0d0*viscosity + 2.0d0*rho0*cs2
        denomDiv = 2.0d0*bulkViscosity + rho0*cs2
        coeffTrace = (2.0d0*viscosity-2.0d0*bulkViscosity)/(2.0d0*denomDiag)

        !$acc parallel loop collapse(2) present(f,rho,u,v,Fx,Fy,Sxx,Sxy,Syy,Sdiv) &
        !$acc& private(i,j,alpha,feqLoc,neqxx,neqxy,neqyy)
        do j = 1, ny
            do i = 1, nx
                neqxx = 0.0d0
                neqxy = 0.0d0
                neqyy = 0.0d0
                do alpha = 0, 8
                    feqLoc = calc_feq(alpha, rho(i,j), u(i,j), v(i,j))
                    neqxx = neqxx + dble(qex(alpha)*qex(alpha))*(f(alpha,i,j)-feqLoc)
                    neqxy = neqxy + dble(qex(alpha)*qey(alpha))*(f(alpha,i,j)-feqLoc)
                    neqyy = neqyy + dble(qey(alpha)*qey(alpha))*(f(alpha,i,j)-feqLoc)
                enddo
                Sdiv(i,j) = -(neqxx+neqyy+u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))/denomDiv
                Sxx(i,j) = -(neqxx+u(i,j)*Fx(i,j))/denomDiag - coeffTrace*Sdiv(i,j)
                Syy(i,j) = -(neqyy+v(i,j)*Fy(i,j))/denomDiag - coeffTrace*Sdiv(i,j)
                Sxy(i,j) = -(2.0d0*neqxy+u(i,j)*Fy(i,j)+v(i,j)*Fx(i,j))/denomShear
            enddo
        enddo
    end subroutine compute_strain_rate
!=============================================================


!=============================================================
    subroutine compute_temperature_gradient()
        use commondata
        implicit none
        integer(kind=4) :: i, j, alpha
        real(kind=8) :: geqLoc, neqx, neqy, denom

        !$acc parallel loop collapse(2) present(g,T,u,v,rho,pressure,Fx,Fy,gradTx,gradTy) &
        !$acc& private(i,j,alpha,geqLoc,neqx,neqy,denom)
        do j = 1, ny
            do i = 1, nx
                neqx = 0.0d0
                neqy = 0.0d0
                do alpha = 0, 8
                    geqLoc = calc_geq(alpha, T(i,j), u(i,j), v(i,j), rho(i,j))
                    neqx = neqx + dble(qex(alpha))*(g(alpha,i,j)-geqLoc)
                    neqy = neqy + dble(qey(alpha))*(g(alpha,i,j)-geqLoc)
                enddo
                denom = cs2*(2.0d0*taugL*(1.0d0-chi_kappa)+chi_kappa) + pressure(i,j)/rho0
                denom = sign(max(abs(denom), tiny(1.0d0)), denom)
                gradTx(i,j) = -(2.0d0*neqx + T(i,j)*Fx(i,j)/rho0)/denom
                gradTy(i,j) = -(2.0d0*neqy + T(i,j)*Fy(i,j)/rho0)/denom
            enddo
        enddo
    end subroutine compute_temperature_gradient
!=============================================================


!=============================================================
    subroutine collision()
        use commondata
        implicit none
        integer(kind=4) :: i, j, alpha
        real(kind=8) :: eu, eF, uF, feqLoc
        real(kind=8) :: Axx, Axy, Ayy, hermite2, phi
        real(kind=8), parameter :: sourcePref=1.0d0-0.5d0/taufL

        !$acc parallel loop collapse(3) present(f_post) private(i,j,alpha)
        do j = 0, ny+1
            do i = 0, nx+1
                do alpha = 0, 8
                    f_post(alpha,i,j) = 0.0d0
                enddo
            enddo
        enddo

        !$acc parallel loop collapse(2) present(f,f_post,rho,u,v,Fx,Fy,Sxx,Sxy,Syy,Sdiv) &
        !$acc& private(i,j,alpha,eu,eF,uF,feqLoc,Axx,Axy,Ayy,hermite2,phi)
        do j = 1, ny
            do i = 1, nx
                uF = u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j)
                Axx = chi_s*Sxx(i,j)+0.5d0*(chi_b-chi_s)*Sdiv(i,j)
                Ayy = chi_s*Syy(i,j)+0.5d0*(chi_b-chi_s)*Sdiv(i,j)
                Axy = chi_s*Sxy(i,j)
                do alpha = 0, 8
                    eu = dble(qex(alpha))*u(i,j)+dble(qey(alpha))*v(i,j)
                    eF = dble(qex(alpha))*Fx(i,j)+dble(qey(alpha))*Fy(i,j)
                    feqLoc = calc_feq(alpha, rho(i,j), u(i,j), v(i,j))
                    hermite2 = (dble(qex(alpha)*qex(alpha))/cs2-1.0d0)*Axx + &
                        (dble(qey(alpha)*qey(alpha))/cs2-1.0d0)*Ayy + &
                        2.0d0*dble(qex(alpha)*qey(alpha))/cs2*Axy
                    phi = sourcePref*qomega(alpha)*(eF/cs2 + eu*eF/(cs2*cs2) - uF/cs2 + rho0*hermite2)
                    f_post(alpha,i,j) = f(alpha,i,j)-omegaF*(f(alpha,i,j)-feqLoc)+phi
                enddo
            enddo
        enddo
    end subroutine collision
!=============================================================


!=============================================================
    subroutine collisionT()
        use commondata
        implicit none
        integer(kind=4) :: i, j, alpha
        real(kind=8) :: geqLoc, scalarSource, psi, eu
        real(kind=8), parameter :: sourcePref=1.0d0-0.5d0/taugL

        !$acc parallel loop collapse(3) present(g_post) private(i,j,alpha)
        do j = 0, ny+1
            do i = 0, nx+1
                do alpha = 0, 8
                    g_post(alpha,i,j) = 0.0d0
                enddo
            enddo
        enddo

        !$acc parallel loop collapse(2) present(g,g_post,T,u,v,rho,pressure,Fx,Fy,gradTx,gradTy) &
        !$acc& private(i,j,alpha,geqLoc,scalarSource,psi,eu)
        do j = 1, ny
            do i = 1, nx
                do alpha = 0, 8
                    geqLoc = calc_geq(alpha, T(i,j), u(i,j), v(i,j), rho(i,j))
                    eu = dble(qex(alpha))*u(i,j)+dble(qey(alpha))*v(i,j)
                    scalarSource = heatSourceQ*(1.0d0+eu/cs2) + &
                        (dble(qex(alpha))*(pressure(i,j)*gradTx(i,j)+T(i,j)*Fx(i,j)) + &
                        dble(qey(alpha))*(pressure(i,j)*gradTy(i,j)+T(i,j)*Fy(i,j)))/(rho0*cs2) + &
                        chi_kappa*(dble(qex(alpha))*gradTx(i,j)+dble(qey(alpha))*gradTy(i,j))
                    psi = sourcePref*qomega(alpha)*scalarSource
                    g_post(alpha,i,j) = g(alpha,i,j)-omegaG*(g(alpha,i,j)-geqLoc)+psi
                enddo
            enddo
        enddo
    end subroutine collisionT
!=============================================================


!=============================================================
    subroutine streaming()
        use commondata
        implicit none
        integer(kind=4) :: i, j, alpha, ip, jp

        !$acc parallel loop collapse(2) present(f,f_post) private(i,j,alpha,ip,jp)
        do j = 1, ny
            do i = 1, nx
                do alpha = 0, 8
                    ip = i-qex(alpha)
                    jp = j-qey(alpha)
                    f(alpha,i,j) = f_post(alpha,ip,jp)
                enddo
            enddo
        enddo
    end subroutine streaming
!=============================================================


!=============================================================
    subroutine streamingT()
        use commondata
        implicit none
        integer(kind=4) :: i, j, alpha, ip, jp

        !$acc parallel loop collapse(2) present(g,g_post) private(i,j,alpha,ip,jp)
        do j = 1, ny
            do i = 1, nx
                do alpha = 0, 8
                    ip = i-qex(alpha)
                    jp = j-qey(alpha)
                    g(alpha,i,j) = g_post(alpha,ip,jp)
                enddo
            enddo
        enddo
    end subroutine streamingT
!=============================================================


!=============================================================
    subroutine bounceback()
        use commondata
        implicit none
        integer(kind=4) :: i, j

#ifdef VerticalWallsNoslip
        !$acc parallel loop present(f,f_post) private(j)
        do j = 1, ny
            f(1,1,j) = f_post(3,1,j)
            f(5,1,j) = f_post(7,1,j)
            f(8,1,j) = f_post(6,1,j)

            f(3,nx,j) = f_post(1,nx,j)
            f(6,nx,j) = f_post(8,nx,j)
            f(7,nx,j) = f_post(5,nx,j)
        enddo
#endif

#ifdef HorizontalWallsNoslip
        !$acc parallel loop present(f,f_post) private(i)
        do i = 1, nx
            f(2,i,1) = f_post(4,i,1)
            f(5,i,1) = f_post(7,i,1)
            f(6,i,1) = f_post(8,i,1)

            f(4,i,ny) = f_post(2,i,ny)
            f(7,i,ny) = f_post(5,i,ny)
            f(8,i,ny) = f_post(6,i,ny)
        enddo
#endif
    end subroutine bounceback
!=============================================================


!=============================================================
    subroutine bouncebackT()
        use commondata
        implicit none
        integer(kind=4) :: i, j

#ifdef VerticalWallsConstT
        !$acc parallel loop present(g,g_post,pressure) private(j)
        do j = 1, ny
            g(1,1,j) = -g_post(3,1,j)+2.0d0*qomega(1)*Thot + &
                2.0d0*qlambdaT(1)*Thot*pressure(1,j)/(rho0*cs2)
            g(5,1,j) = -g_post(7,1,j)+2.0d0*qomega(5)*Thot + &
                2.0d0*qlambdaT(5)*Thot*pressure(1,j)/(rho0*cs2)
            g(8,1,j) = -g_post(6,1,j)+2.0d0*qomega(8)*Thot + &
                2.0d0*qlambdaT(8)*Thot*pressure(1,j)/(rho0*cs2)

            g(3,nx,j) = -g_post(1,nx,j)+2.0d0*qomega(3)*Tcold + &
                2.0d0*qlambdaT(3)*Tcold*pressure(nx,j)/(rho0*cs2)
            g(6,nx,j) = -g_post(8,nx,j)+2.0d0*qomega(6)*Tcold + &
                2.0d0*qlambdaT(6)*Tcold*pressure(nx,j)/(rho0*cs2)
            g(7,nx,j) = -g_post(5,nx,j)+2.0d0*qomega(7)*Tcold + &
                2.0d0*qlambdaT(7)*Tcold*pressure(nx,j)/(rho0*cs2)
        enddo
#endif

#ifdef HorizontalWallsAdiabatic
        !$acc parallel loop present(g,g_post) private(i)
        do i = 1, nx
            g(2,i,1) = g_post(4,i,1)
            g(5,i,1) = g_post(7,i,1)
            g(6,i,1) = g_post(8,i,1)

            g(4,i,ny) = g_post(2,i,ny)
            g(7,i,ny) = g_post(5,i,ny)
            g(8,i,ny) = g_post(6,i,ny)
        enddo
#endif

#if defined(VerticalWallsConstT) && defined(HorizontalWallsAdiabatic)
        ! D2Q9 corner diagonals belong to both a constant-temperature side wall
        ! and an adiabatic horizontal wall.  Let the side-wall Dirichlet condition
        ! set the four conflicting diagonal populations after the horizontal pass.
        !$acc serial present(g,g_post,pressure)
        g(5,1,1) = -g_post(7,1,1)+2.0d0*qomega(5)*Thot + &
            2.0d0*qlambdaT(5)*Thot*pressure(1,1)/(rho0*cs2)
        g(8,1,ny) = -g_post(6,1,ny)+2.0d0*qomega(8)*Thot + &
            2.0d0*qlambdaT(8)*Thot*pressure(1,ny)/(rho0*cs2)
        g(6,nx,1) = -g_post(8,nx,1)+2.0d0*qomega(6)*Tcold + &
            2.0d0*qlambdaT(6)*Tcold*pressure(nx,1)/(rho0*cs2)
        g(7,nx,ny) = -g_post(5,nx,ny)+2.0d0*qomega(7)*Tcold + &
            2.0d0*qlambdaT(7)*Tcold*pressure(nx,ny)/(rho0*cs2)
        !$acc end serial
#endif
    end subroutine bouncebackT
!=============================================================


!=============================================================
    subroutine macro()
        use commondata
        implicit none
        integer(kind=4) :: i, j
        real(kind=8) :: momx, momy

        !$acc parallel loop collapse(2) present(f,rho,u,v,Fx,Fy) private(i,j,momx,momy)
        do j = 1, ny
            do i = 1, nx
                rho(i,j) = rho0 + f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+ &
                    f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
                momx = f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
                momy = f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)
                u(i,j) = (momx+0.5d0*Fx(i,j))/rho0
                v(i,j) = (momy+0.5d0*Fy(i,j))/rho0
            enddo
        enddo
    end subroutine macro
!=============================================================


!=============================================================
    subroutine macroT()
        use commondata
        implicit none
        integer(kind=4) :: i, j

        !$acc parallel loop collapse(2) present(g,T) private(i,j)
        do j = 1, ny
            do i = 1, nx
                T(i,j) = g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)+ &
                    g(5,i,j)+g(6,i,j)+g(7,i,j)+g(8,i,j) + 0.5d0*heatSourceQ
            enddo
        enddo
    end subroutine macroT
!=============================================================


!=============================================================
    subroutine check()
        use commondata
        implicit none
        integer(kind=4) :: i, j
        real(kind=8) :: errUSum, refUSum, errTSum, refTSum

        errUSum = 0.0d0
        refUSum = 0.0d0
        errTSum = 0.0d0
        refTSum = 0.0d0
        do j = 1, ny
            do i = 1, nx
                errUSum = errUSum+(u(i,j)-up(i,j))**2+(v(i,j)-vp(i,j))**2
                refUSum = refUSum+u(i,j)**2+v(i,j)**2
                errTSum = errTSum+(T(i,j)-Tp(i,j))**2
                refTSum = refTSum+T(i,j)**2
            enddo
        enddo
        errorU = dsqrt(errUSum/max(refUSum,tiny(1.0d0)))
        errorT = dsqrt(errTSum/max(refTSum,tiny(1.0d0)))
        up = u
        vp = v
        Tp = T

        write(*,'(A,I10,2(A,ES14.6))') "itc=", itc, " errorU=", errorU, " errorT=", errorT
    end subroutine check
!=============================================================


!=============================================================
    subroutine calNuRe()
        use commondata
        implicit none
        integer(kind=4) :: i, j, imid, jmid
        real(kind=8) :: sumHot, sumCold, sumMiddle, sumRe, sumGlobal
        real(kind=8) :: localNu

        imid = (nx+1)/2
        jmid = (ny+1)/2
        sumHot = 0.0d0
        sumCold = 0.0d0
        sumMiddle = 0.0d0
        sumRe = 0.0d0
        sumGlobal = 0.0d0
        umax_mid = 0.0d0
        vmax_mid = 0.0d0
        Numax_mid = -huge(1.0d0)
        y_umax_mid = yp(1)
        x_vmax_mid = xp(1)
        y_Numax_mid = yp(1)

        do j = 1, ny
            localNu = -lengthUnit*gradTx(1,j)/deltaT
            sumHot = sumHot + localNu
            localNu = -lengthUnit*gradTx(nx,j)/deltaT
            sumCold = sumCold + localNu
            localNu = -lengthUnit*gradTx(imid,j)/deltaT
            sumMiddle = sumMiddle + localNu
            if(localNu.GT.Numax_mid) then
                Numax_mid = localNu
                y_Numax_mid = yp(j)
            endif
            if(abs(u(imid,j)).GT.abs(umax_mid)) then
                umax_mid = u(imid,j)
                y_umax_mid = yp(j)
            endif
        enddo

        do i = 1, nx
            if(abs(v(i,jmid)).GT.abs(vmax_mid)) then
                vmax_mid = v(i,jmid)
                x_vmax_mid = xp(i)
            endif
        enddo

        do j = 1, ny
            do i = 1, nx
                sumRe = sumRe + u(i,j)*u(i,j)+v(i,j)*v(i,j)
                sumGlobal = sumGlobal + u(i,j)*(T(i,j)-Tref)
            enddo
        enddo

        Nu_hot = sumHot/dble(ny)
        Nu_cold = sumCold/dble(ny)
        Nu_middle = sumMiddle/dble(ny)
        Nu_global = 1.0d0 + sumGlobal/dble(nx*ny)*lengthUnit/diffusivity
        ReVolAvg = dsqrt(sumRe/dble(nx*ny))*lengthUnit/viscosity
        umax_mid = umax_mid/velocityUnit
        vmax_mid = vmax_mid/velocityUnit
    end subroutine calNuRe
!=============================================================


!=============================================================
    subroutine append_history()
        use commondata
        implicit none

        open(unit=11, file=trim(historyFile), status='old', position='append', action='write')
        write(11,'(I12,1X,7ES24.16E3)') itc, dble(itc)/timeUnit, errorU, errorT, &
            Nu_hot, Nu_cold, Nu_middle, ReVolAvg
        close(11)

        open(unit=00, file=trim(settingsFile), status='old', position='append', action='write')
        write(00,'(A,I12,2(A,ES16.8))') "itc=", itc, " errorU=", errorU, " errorT=", errorT
        write(00,'(4(A,ES16.8))') " Nu_hot=", Nu_hot, " Nu_cold=", Nu_cold, &
            " Nu_middle=", Nu_middle, " Re=", ReVolAvg
        close(00)
    end subroutine append_history
!=============================================================


!=============================================================
    subroutine output_Tecplot()
        use commondata
        implicit none
        integer(kind=4) :: i, j

        open(unit=41, file=trim(pltFile), status='replace', action='write', form='formatted')
        write(41,'(A)') 'TITLE = "2D OpenACC LBM-CDE side-heated cavity"'
        write(41,'(A)') 'VARIABLES = "X" "Y" "U" "V" "T" "RHO" "GradTx" "GradTy"'
        write(41,'(A,I0,A,I0,A)') 'ZONE I=', nx, ', J=', ny, ', F=POINT'
        do j = 1, ny
            do i = 1, nx
                write(41,'(8ES24.16E3)') xp(i), yp(j), u(i,j), v(i,j), T(i,j), rho(i,j), gradTx(i,j), gradTy(i,j)
            enddo
        enddo
        close(41)
    end subroutine output_Tecplot
!=============================================================


!=============================================================
    subroutine final_log(cpuSeconds, wallSeconds)
        use commondata
        implicit none
        real(kind=8), intent(in) :: cpuSeconds, wallSeconds

        write(*,'(A,I12)') "Final itc =", itc
        write(*,'(A,ES16.8,A,ES16.8)') "Final errorU/errorT =", errorU, " / ", errorT
        write(*,'(A,ES16.8,A,ES16.8)') "Nu hot/cold =", Nu_hot, " / ", Nu_cold
        write(*,'(A,ES16.8,A,ES16.8)') "umax/vmax mid =", umax_mid, " / ", vmax_mid

        open(unit=00, file=trim(settingsFile), status='old', position='append', action='write')
        write(00,*) "Final itc =", itc
        write(00,*) "Final errorU =", real(errorU,kind=8), "errorT =", real(errorT,kind=8)
        write(00,*) "Nu_hot =", real(Nu_hot,kind=8), "Nu_cold =", real(Nu_cold,kind=8)
        write(00,*) "Nu_middle =", real(Nu_middle,kind=8), "Nu_global =", real(Nu_global,kind=8)
        write(00,*) "ReVolAvg =", real(ReVolAvg,kind=8)
        write(00,*) "umax_mid/u0 =", real(umax_mid,kind=8), "y =", real(y_umax_mid,kind=8)
        write(00,*) "vmax_mid/u0 =", real(vmax_mid,kind=8), "x =", real(x_vmax_mid,kind=8)
        write(00,*) "Numax_mid =", real(Numax_mid,kind=8), "y =", real(y_Numax_mid,kind=8)
        write(00,*) "CPU seconds =", real(cpuSeconds,kind=8)
        write(00,*) "Wall seconds =", real(wallSeconds,kind=8)
        close(00)
    end subroutine final_log
!=============================================================


!=============================================================
    subroutine finalize_arrays()
        use commondata
        implicit none

        deallocate(u, v, T, rho, pressure)
        deallocate(up, vp, Tp)
        deallocate(Fx, Fy)
        deallocate(gradTx, gradTy)
        deallocate(Sxx, Sxy, Syy, Sdiv)
        deallocate(f, f_post, g, g_post)
    end subroutine finalize_arrays
!=============================================================
