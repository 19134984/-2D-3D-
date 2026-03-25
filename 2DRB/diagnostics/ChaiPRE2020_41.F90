!=============================================================
!    ChaiPRE2020.F90
!    2D sidewall-heated natural convection
!    CPU + OpenMP
!    D2Q9 MRT for flow + D2Q5 MRT for temperature
!=============================================================

module commondata
    implicit none

    integer(kind=4), parameter :: nx = 41, ny = 41
    integer(kind=4), parameter :: qf = 9, qt = 5
    integer(kind=4), parameter :: itc_max = 5000000
    integer(kind=4), parameter :: checkInterval = 200
    integer(kind=4), parameter :: outputEveryChecks = 10
    integer(kind=4), parameter :: numThreads = 1
    integer(kind=4), parameter :: nxHalf = int(0.5d0 * (dble(nx) + 1.0d0))
    integer(kind=4), parameter :: nyHalf = int(0.5d0 * (dble(ny) + 1.0d0))

    real(kind=8), parameter :: pi = acos(-1.0d0)
    real(kind=8), parameter :: dt = 1.0d0
    real(kind=8), parameter :: cs2Flow = 1.0d0 / 3.0d0
    real(kind=8), parameter :: cs4Flow = cs2Flow * cs2Flow
    real(kind=8), parameter :: cs2Temp = 1.0d0 / 3.0d0
    real(kind=8), parameter :: lengthUnit = dble(nx)

    real(kind=8), parameter :: Rayleigh = 1.0d3
    real(kind=8), parameter :: Prandtl = 0.71d0
    real(kind=8), parameter :: Mach = 0.1d0
    real(kind=8), parameter :: Thot = 0.5d0
    real(kind=8), parameter :: Tcold = -0.5d0
    real(kind=8), parameter :: Tref = 0.5d0 * (Thot + Tcold)
    real(kind=8), parameter :: rho0 = 1.0d0
    real(kind=8), parameter :: tauf = 0.5d0 + Mach * lengthUnit * dsqrt(3.0d0 * Prandtl / Rayleigh)
    real(kind=8), parameter :: viscosity = (tauf - 0.5d0) * cs2Flow
    real(kind=8), parameter :: diffusivity = viscosity / Prandtl
    real(kind=8), parameter :: S2s = 1.0d0 / (viscosity / (cs2Flow * dt) + 0.5d0)
    real(kind=8), parameter :: S2b = 1.0d0
    real(kind=8), parameter :: SepsFlow = 1.0d0
    real(kind=8), parameter :: SqFlow = 8.0d0 * (2.0d0 * tauf - 1.0d0) / (8.0d0 * tauf - 1.0d0)
    real(kind=8), parameter :: sT = 1.0d0 / (diffusivity / (cs2Temp * dt) + 0.5d0)
    real(kind=8), parameter :: gBeta1 = Rayleigh * viscosity * diffusivity / lengthUnit
    real(kind=8), parameter :: gBeta = gBeta1 / (lengthUnit * lengthUnit)
    real(kind=8), parameter :: timeUnit = dsqrt(lengthUnit / gBeta)
    real(kind=8), parameter :: velocityUnit = dsqrt(gBeta * lengthUnit)
    real(kind=8), parameter :: epsU = 1.0d-7
    real(kind=8), parameter :: epsT = 1.0d-7

    logical, parameter :: outputPltFile = .true.
    logical, parameter :: useGjTemp = .true.
    logical, parameter :: useChaiStrict = .true.

    integer(kind=4) :: exFlow(0:qf-1), eyFlow(0:qf-1)
    integer(kind=4) :: exTemp(0:qt-1), eyTemp(0:qt-1)
    integer(kind=4) :: oppFlow(0:qf-1), oppTemp(0:qt-1)

    real(kind=8) :: wFlow(0:qf-1), wTemp(0:qt-1)
    real(kind=8) :: Mflow(0:qf-1,0:qf-1), MinvFlow(0:qf-1,0:qf-1), Sflow(0:qf-1)
    real(kind=8) :: Mtemp(0:qt-1,0:qt-1), MinvTemp(0:qt-1,0:qt-1), Stemp(0:qt-1)

    real(kind=8) :: xp(0:nx+1), yp(0:ny+1)
    real(kind=8), allocatable :: rho(:,:), u(:,:), v(:,:), T(:,:)
    real(kind=8), allocatable :: uPrev(:,:), vPrev(:,:), TPrev(:,:)
    real(kind=8), allocatable :: uStage(:,:), vStage(:,:), TStage(:,:)
    real(kind=8), allocatable :: Fx(:,:), Fy(:,:)
    real(kind=8), allocatable :: BxPrev(:,:), ByPrev(:,:)
    real(kind=8), allocatable :: f(:,:,:), fPost(:,:,:), fNext(:,:,:)
    real(kind=8), allocatable :: g(:,:,:), gPost(:,:,:), gNext(:,:,:)

    integer(kind=4) :: itc
    integer(kind=4) :: pltFileNum
    integer(kind=4) :: dimensionlessTime
    real(kind=8) :: errorU, errorT
    real(kind=8) :: Nu_global, Nu_hot, Nu_cold, Nu_middle
    real(kind=8) :: Nu_hot_max, Nu_hot_min, Nu_hot_max_position, Nu_hot_min_position
    real(kind=8), allocatable :: NuVolAvg(:), ReVolAvg(:)
    character(len=100) :: pltFolderPrefix = "./ChaiPRE2020"

contains

    pure real(kind=8) function basis_flow(mode, cx, cy)
        implicit none
        integer(kind=4), intent(in) :: mode, cx, cy
        real(kind=8) :: x, y

        x = dble(cx)
        y = dble(cy)

        select case (mode)
        case (0)
            basis_flow = 1.0d0
        case (1)
            basis_flow = x
        case (2)
            basis_flow = y
        case (3)
            basis_flow = x * x + y * y
        case (4)
            basis_flow = x * x - y * y
        case (5)
            basis_flow = x * y
        case (6)
            basis_flow = x * x * y
        case (7)
            basis_flow = x * y * y
        case default
            basis_flow = x * x * y * y
        end select
    end function basis_flow

    pure real(kind=8) function basis_temp(mode, cx, cy)
        implicit none
        integer(kind=4), intent(in) :: mode, cx, cy
        real(kind=8) :: x, y

        x = dble(cx)
        y = dble(cy)

        select case (mode)
        case (0)
            basis_temp = 1.0d0
        case (1)
            basis_temp = x
        case (2)
            basis_temp = y
        case (3)
            basis_temp = x * x + y * y
        case default
            basis_temp = x * x - y * y
        end select
    end function basis_temp

end module commondata


program main
    use omp_lib
    use commondata
    implicit none

    real(kind=8) :: timeStart, timeEnd, wallStart, wallEnd
    integer(kind=4) :: outputCounter

    call build_lattices()
    call build_mrt_matrices()
    call initial()

    open(unit=00, file="SimulationSettings.txt", status="unknown")
    write(00,'(a)') "ChaiPRE2020 sidewall-heated natural convection"
    write(00,'(a,2(1x,i8))') "Grid nx ny =", nx, ny
    write(00,'(a,1x,es16.8)') "Ra =", Rayleigh
    write(00,'(a,1x,es16.8)') "Pr =", Prandtl
    write(00,'(a,1x,es16.8)') "Mach =", Mach
    write(00,'(a,1x,es16.8)') "Thot =", Thot
    write(00,'(a,1x,es16.8)') "Tcold =", Tcold
    write(00,'(a,1x,es16.8)') "Tref =", Tref
    write(00,'(a,1x,es16.8)') "nu =", viscosity
    write(00,'(a,1x,es16.8)') "alpha =", diffusivity
    write(00,'(a,1x,es16.8)') "S2s =", S2s
    write(00,'(a,1x,es16.8)') "S2b =", S2b
    write(00,'(a,1x,es16.8)') "SepsFlow =", SepsFlow
    write(00,'(a,1x,es16.8)') "SqFlow =", SqFlow
    write(00,'(a,1x,es16.8)') "sT =", sT
    write(00,'(a,1x,es16.8)') "gBeta =", gBeta
    write(00,'(a,1x,l1)') "useChaiStrict =", useChaiStrict
    call OMP_set_num_threads(numThreads)
    write(00,'(a,1x,i8)') "OpenMP threads =", OMP_get_max_threads()
    close(00)

    outputCounter = 0
    call cpu_time(timeStart)
    wallStart = OMP_get_wtime()

    do while (((errorU .gt. epsU) .or. (errorT .gt. epsT)) .and. (itc .lt. itc_max))
        itc = itc + 1

        call update_buoyancy_force()
        call collision_flow()
        call streaming_flow()
        call bounceback_flow()
        call macro_flow()

        call collision_temp()
        call streaming_temp()
        call bounceback_temp()
        call macro_temp()

        if (mod(itc, checkInterval) .eq. 0) then
            call check()
            call calNuRe()
            outputCounter = outputCounter + 1

            write(*,'(a,1x,i10,2(1x,a,1x,es16.8))') "itc =", itc, "errU =", errorU, "errT =", errorT
            open(unit=00, file="SimulationSettings.txt", status="unknown", position="append")
            write(00,'(a,1x,i10,2(1x,a,1x,es16.8))') "itc =", itc, "errU =", errorU, "errT =", errorT
            close(00)

            if (mod(outputCounter, outputEveryChecks) .eq. 0) then
                call SideHeatedcalc_Nu_global()
                call SideHeatedcalc_Nu_wall_avg()
                call SideHeatedcalc_umid_max()
                call SideHeatedcalc_vmid_max()
                if (outputPltFile) call output_Tecplot()
            endif
        endif
    end do

    call cpu_time(timeEnd)
    wallEnd = OMP_get_wtime()

    call SideHeatedcalc_Nu_global()
    call SideHeatedcalc_Nu_wall_avg()
    call SideHeatedcalc_umid_max()
    call SideHeatedcalc_vmid_max()
    call output_centerline_profiles()
    if (outputPltFile) call output_Tecplot()

    open(unit=00, file="SimulationSettings.txt", status="unknown", position="append")
    write(00,'(a,1x,i12)') "Finished itc =", itc
    write(00,'(a,1x,es16.8)') "Final errorU =", errorU
    write(00,'(a,1x,es16.8)') "Final errorT =", errorT
    write(00,'(a,1x,es16.8)') "CPU time =", timeEnd - timeStart
    write(00,'(a,1x,es16.8)') "Wall time =", wallEnd - wallStart
    close(00)

end program main


subroutine write_timestamp(unit_id, label)
    implicit none
    integer(kind=4), intent(in) :: unit_id
    character(len=*), intent(in) :: label
    character(len=8) :: d
    character(len=10) :: t
    character(len=5) :: z

    call date_and_time(date=d, time=t, zone=z)
    write(unit_id,'(a,1x,a4,"-",a2,"-",a2,1x,a2,":",a2,":",a2,1x,a)') trim(label), &
        d(1:4), d(5:6), d(7:8), t(1:2), t(3:4), t(5:6), trim(z)
end subroutine write_timestamp


subroutine build_lattices()
    use commondata
    implicit none

    exFlow = (/ 0, 1, 0, -1, 0, 1, -1, -1, 1 /)
    eyFlow = (/ 0, 0, 1, 0, -1, 1, 1, -1, -1 /)
    oppFlow = (/ 0, 3, 4, 1, 2, 7, 8, 5, 6 /)
    wFlow = (/ 4.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, &
               1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0 /)

    exTemp = (/ 0, 1, 0, -1, 0 /)
    eyTemp = (/ 0, 0, 1, 0, -1 /)
    oppTemp = (/ 0, 3, 4, 1, 2 /)
    wTemp = (/ 1.0d0/3.0d0, 1.0d0/6.0d0, 1.0d0/6.0d0, 1.0d0/6.0d0, 1.0d0/6.0d0 /)
end subroutine build_lattices


subroutine build_mrt_matrices()
    use commondata
    implicit none
    Mflow = 0.0d0
    Mflow(0,:) = (/ 1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0 /)
    Mflow(1,:) = (/ -4.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0,  2.0d0,  2.0d0,  2.0d0,  2.0d0 /)
    Mflow(2,:) = (/ 4.0d0, -2.0d0, -2.0d0, -2.0d0, -2.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0 /)
    Mflow(3,:) = (/ 0.0d0,  1.0d0,  0.0d0, -1.0d0,  0.0d0,  1.0d0, -1.0d0, -1.0d0,  1.0d0 /)
    Mflow(4,:) = (/ 0.0d0, -2.0d0,  0.0d0,  2.0d0,  0.0d0,  1.0d0, -1.0d0, -1.0d0,  1.0d0 /)
    Mflow(5,:) = (/ 0.0d0,  0.0d0,  1.0d0,  0.0d0, -1.0d0,  1.0d0,  1.0d0, -1.0d0, -1.0d0 /)
    Mflow(6,:) = (/ 0.0d0,  0.0d0, -2.0d0,  0.0d0,  2.0d0,  1.0d0,  1.0d0, -1.0d0, -1.0d0 /)
    Mflow(7,:) = (/ 0.0d0,  1.0d0, -1.0d0,  1.0d0, -1.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0 /)
    Mflow(8,:) = (/ 0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  1.0d0, -1.0d0,  1.0d0, -1.0d0 /)

    Mtemp = 0.0d0
    Mtemp(0,:) = (/ 1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0 /)
    Mtemp(1,:) = (/ 0.0d0,  1.0d0,  0.0d0, -1.0d0,  0.0d0 /)
    Mtemp(2,:) = (/ 0.0d0,  0.0d0,  1.0d0,  0.0d0, -1.0d0 /)
    Mtemp(3,:) = (/ -4.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0 /)
    Mtemp(4,:) = (/ 0.0d0,  1.0d0, -1.0d0,  1.0d0, -1.0d0 /)

    call invert_matrix(qf, Mflow, MinvFlow)
    call invert_matrix(qt, Mtemp, MinvTemp)

    Sflow = (/ 0.0d0, S2b, SepsFlow, 0.0d0, SqFlow, 0.0d0, SqFlow, S2s, S2s /)
    Stemp = (/ 0.0d0, sT, sT, 1.0d0, 1.0d0 /)
end subroutine build_mrt_matrices


subroutine initial()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k
    real(kind=8) :: feq(0:qf-1), geq(0:qt-1)

    allocate(rho(1:nx,1:ny), u(1:nx,1:ny), v(1:nx,1:ny), T(1:nx,1:ny))
    allocate(uPrev(1:nx,1:ny), vPrev(1:nx,1:ny), TPrev(1:nx,1:ny))
    allocate(uStage(1:nx,1:ny), vStage(1:nx,1:ny), TStage(1:nx,1:ny))
    allocate(Fx(1:nx,1:ny), Fy(1:nx,1:ny))
    allocate(BxPrev(1:nx,1:ny), ByPrev(1:nx,1:ny))
    allocate(f(0:qf-1,0:nx+1,0:ny+1), fPost(0:qf-1,0:nx+1,0:ny+1), fNext(0:qf-1,0:nx+1,0:ny+1))
    allocate(g(0:qt-1,0:nx+1,0:ny+1), gPost(0:qt-1,0:nx+1,0:ny+1), gNext(0:qt-1,0:nx+1,0:ny+1))
    allocate(NuVolAvg(0:itc_max / checkInterval + 5), ReVolAvg(0:itc_max / checkInterval + 5))

    xp(0) = 0.0d0
    xp(nx+1) = dble(nx) / lengthUnit
    do i = 1, nx
        xp(i) = (dble(i) - 0.5d0) / lengthUnit
    end do

    yp(0) = 0.0d0
    yp(ny+1) = dble(ny) / lengthUnit
    do j = 1, ny
        yp(j) = (dble(j) - 0.5d0) / lengthUnit
    end do

    f = 0.0d0
    fPost = 0.0d0
    fNext = 0.0d0
    g = 0.0d0
    gPost = 0.0d0
    gNext = 0.0d0

    do j = 1, ny
        do i = 1, nx
            rho(i,j) = rho0
            u(i,j) = 0.0d0
            v(i,j) = 0.0d0
            T(i,j) = Thot + (Tcold - Thot) * xp(i)
            Fx(i,j) = 0.0d0
            Fy(i,j) = rho(i,j) * gBeta * (T(i,j) - Tref)
            BxPrev(i,j) = u(i,j) * T(i,j)
            ByPrev(i,j) = v(i,j) * T(i,j)

            call compute_feq_flow(rho(i,j), u(i,j), v(i,j), feq)
            do k = 0, qf-1
                f(k,i,j) = feq(k)
            end do

            call compute_geq_temp(T(i,j), u(i,j), v(i,j), geq)
            do k = 0, qt-1
                g(k,i,j) = geq(k)
            end do
        end do
    end do

    uPrev = u
    vPrev = v
    TPrev = T
    uStage = u
    vStage = v
    TStage = T
    errorU = 1.0d0
    errorT = 1.0d0
    itc = 0
    pltFileNum = 0
    dimensionlessTime = 0
    NuVolAvg = 0.0d0
    ReVolAvg = 0.0d0
end subroutine initial


subroutine update_buoyancy_force()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !$omp parallel do default(none) shared(Fx,Fy,T,rho) private(i,j)
    do j = 1, ny
        do i = 1, nx
            Fx(i,j) = 0.0d0
            Fy(i,j) = rho(i,j) * gBeta * (T(i,j) - Tref)
        end do
    end do
    !$omp end parallel do
end subroutine update_buoyancy_force


subroutine collision_flow()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k
    real(kind=8) :: m(0:qf-1), meq(0:qf-1), mPost(0:qf-1)
    real(kind=8) :: s(0:qf-1), fSource(0:qf-1)

    if (useChaiStrict) then
        call collision_flow_strict()
        return
    endif

    !$omp parallel do default(none) shared(f,fPost,rho,u,v,Fx,Fy,T) &
    !$omp private(i,j,k,m,meq,mPost,s,fSource)
    do j = 1, ny
        do i = 1, nx
            m(0) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
            m(1) = -4.0d0*f(0,i,j)-f(1,i,j)-f(2,i,j)-f(3,i,j)-f(4,i,j)+ &
                2.0d0*(f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j))
            m(2) = 4.0d0*f(0,i,j)-2.0d0*(f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j))+ &
                f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
            m(3) = f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
            m(4) = -2.0d0*f(1,i,j)+2.0d0*f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
            m(5) = f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)
            m(6) = -2.0d0*f(2,i,j)+2.0d0*f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)
            m(7) = f(1,i,j)-f(2,i,j)+f(3,i,j)-f(4,i,j)
            m(8) = f(5,i,j)-f(6,i,j)+f(7,i,j)-f(8,i,j)

            meq(0) = rho(i,j)
            meq(1) = rho(i,j) * (-2.0d0 + 3.0d0 * (u(i,j)*u(i,j) + v(i,j)*v(i,j)))
            meq(2) = rho(i,j) * ( 1.0d0 - 3.0d0 * (u(i,j)*u(i,j) + v(i,j)*v(i,j)))
            meq(3) = rho(i,j) * u(i,j)
            meq(4) = -rho(i,j) * u(i,j)
            meq(5) = rho(i,j) * v(i,j)
            meq(6) = -rho(i,j) * v(i,j)
            meq(7) = rho(i,j) * (u(i,j)*u(i,j) - v(i,j)*v(i,j))
            meq(8) = rho(i,j) * (u(i,j) * v(i,j))

            s(0) = 0.0d0
            s(1) = S2b
            s(2) = SepsFlow
            s(3) = 0.0d0
            s(4) = SqFlow
            s(5) = 0.0d0
            s(6) = SqFlow
            s(7) = S2s
            s(8) = S2s

            Fx(i,j) = 0.0d0
            Fy(i,j) = rho(i,j) * gBeta * (T(i,j) - Tref)

            fSource(0) = 0.0d0
            fSource(1) = (6.0d0 - 3.0d0 * s(1)) * (u(i,j) * Fx(i,j) + v(i,j) * Fy(i,j))
            fSource(2) = -(6.0d0 - 3.0d0 * s(2)) * (u(i,j) * Fx(i,j) + v(i,j) * Fy(i,j))
            fSource(3) = (1.0d0 - 0.5d0 * s(3)) * Fx(i,j)
            fSource(4) = -(1.0d0 - 0.5d0 * s(4)) * Fx(i,j)
            fSource(5) = (1.0d0 - 0.5d0 * s(5)) * Fy(i,j)
            fSource(6) = -(1.0d0 - 0.5d0 * s(6)) * Fy(i,j)
            fSource(7) = (2.0d0 - s(7)) * (u(i,j) * Fx(i,j) - v(i,j) * Fy(i,j))
            fSource(8) = (1.0d0 - 0.5d0 * s(8)) * (u(i,j) * Fy(i,j) + v(i,j) * Fx(i,j))

            do k = 0, qf-1
                mPost(k) = m(k) - s(k) * (m(k) - meq(k)) + fSource(k)
            end do

            fPost(0,i,j) = mPost(0)/9.0d0 - mPost(1)/9.0d0 + mPost(2)/9.0d0
            fPost(1,i,j) = mPost(0)/9.0d0 - mPost(1)/36.0d0 - mPost(2)/18.0d0 + mPost(3)/6.0d0 - &
                mPost(4)/6.0d0 + mPost(7)/4.0d0
            fPost(2,i,j) = mPost(0)/9.0d0 - mPost(1)/36.0d0 - mPost(2)/18.0d0 + mPost(5)/6.0d0 - &
                mPost(6)/6.0d0 - mPost(7)/4.0d0
            fPost(3,i,j) = mPost(0)/9.0d0 - mPost(1)/36.0d0 - mPost(2)/18.0d0 - mPost(3)/6.0d0 + &
                mPost(4)/6.0d0 + mPost(7)/4.0d0
            fPost(4,i,j) = mPost(0)/9.0d0 - mPost(1)/36.0d0 - mPost(2)/18.0d0 - mPost(5)/6.0d0 + &
                mPost(6)/6.0d0 - mPost(7)/4.0d0
            fPost(5,i,j) = mPost(0)/9.0d0 + mPost(1)/18.0d0 + mPost(2)/36.0d0 + mPost(3)/6.0d0 + &
                mPost(4)/12.0d0 + mPost(5)/6.0d0 + mPost(6)/12.0d0 + mPost(8)/4.0d0
            fPost(6,i,j) = mPost(0)/9.0d0 + mPost(1)/18.0d0 + mPost(2)/36.0d0 - mPost(3)/6.0d0 - &
                mPost(4)/12.0d0 + mPost(5)/6.0d0 + mPost(6)/12.0d0 - mPost(8)/4.0d0
            fPost(7,i,j) = mPost(0)/9.0d0 + mPost(1)/18.0d0 + mPost(2)/36.0d0 - mPost(3)/6.0d0 - &
                mPost(4)/12.0d0 - mPost(5)/6.0d0 - mPost(6)/12.0d0 + mPost(8)/4.0d0
            fPost(8,i,j) = mPost(0)/9.0d0 + mPost(1)/18.0d0 + mPost(2)/36.0d0 + mPost(3)/6.0d0 + &
                mPost(4)/12.0d0 - mPost(5)/6.0d0 - mPost(6)/12.0d0 - mPost(8)/4.0d0
        end do
    end do
    !$omp end parallel do
end subroutine collision_flow


subroutine collision_flow_strict()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k, l
    real(kind=8) :: feq(0:qf-1), src(0:qf-1), gaux(0:qf-1)
    real(kind=8) :: m(0:qf-1), meq(0:qf-1), mF(0:qf-1), mG(0:qf-1), mPost(0:qf-1)

    !$omp parallel do default(none) shared(f,fPost,rho,u,v,Fx,Fy,Mflow,MinvFlow,Sflow) &
    !$omp private(i,j,k,l,feq,src,gaux,m,meq,mF,mG,mPost)
    do j = 1, ny
        do i = 1, nx
            call compute_feq_flow(rho(i,j), u(i,j), v(i,j), feq)
            call compute_flow_source(u(i,j), v(i,j), Fx(i,j), Fy(i,j), src, gaux)

            do k = 0, qf-1
                m(k) = 0.0d0
                meq(k) = 0.0d0
                mF(k) = 0.0d0
                mG(k) = 0.0d0
                do l = 0, qf-1
                    m(k) = m(k) + Mflow(k,l) * f(l,i,j)
                    meq(k) = meq(k) + Mflow(k,l) * feq(l)
                    mF(k) = mF(k) + Mflow(k,l) * src(l)
                    mG(k) = mG(k) + Mflow(k,l) * gaux(l)
                end do
                mPost(k) = m(k) - Sflow(k) * (m(k) - meq(k)) + (1.0d0 - 0.5d0 * Sflow(k)) * mF(k) + mG(k)
            end do

            do k = 0, qf-1
                fPost(k,i,j) = 0.0d0
                do l = 0, qf-1
                    fPost(k,i,j) = fPost(k,i,j) + MinvFlow(k,l) * mPost(l)
                end do
            end do
        end do
    end do
    !$omp end parallel do
end subroutine collision_flow_strict


subroutine streaming_flow()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k, ii, jj

    fNext = 0.0d0

    !$omp parallel do default(none) shared(fPost,fNext,exFlow,eyFlow) private(i,j,k,ii,jj)
    do j = 1, ny
        do i = 1, nx
            do k = 0, qf-1
                ii = i + exFlow(k)
                jj = j + eyFlow(k)
                fNext(k,ii,jj) = fPost(k,i,j)
            end do
        end do
    end do
    !$omp end parallel do

    f = fNext
end subroutine streaming_flow


subroutine bounceback_flow()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !$omp parallel do default(none) shared(f,fNext) private(j)
    do j = 1, ny
        f(1,1,j) = fNext(3,0,j)
        f(5,1,j) = fNext(7,0,j)
        f(8,1,j) = fNext(6,0,j)

        f(3,nx,j) = fNext(1,nx+1,j)
        f(6,nx,j) = fNext(8,nx+1,j)
        f(7,nx,j) = fNext(5,nx+1,j)
    end do
    !$omp end parallel do

    !$omp parallel do default(none) shared(f,fNext) private(i)
    do i = 1, nx
        f(2,i,1) = fNext(4,i,0)
        f(5,i,1) = fNext(7,i,0)
        f(6,i,1) = fNext(8,i,0)

        f(4,i,ny) = fNext(2,i,ny+1)
        f(7,i,ny) = fNext(5,i,ny+1)
        f(8,i,ny) = fNext(6,i,ny+1)
    end do
    !$omp end parallel do
end subroutine bounceback_flow


subroutine macro_flow()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k
    real(kind=8) :: rhoLocal, ux, uy

    !$omp parallel do default(none) shared(f,rho,u,v,Fx,Fy,exFlow,eyFlow) private(i,j,k,rhoLocal,ux,uy)
    do j = 1, ny
        do i = 1, nx
            rhoLocal = 0.0d0
            ux = 0.0d0
            uy = 0.0d0
            do k = 0, qf-1
                rhoLocal = rhoLocal + f(k,i,j)
                ux = ux + dble(exFlow(k)) * f(k,i,j)
                uy = uy + dble(eyFlow(k)) * f(k,i,j)
            end do

            rho(i,j) = rhoLocal
            if (useChaiStrict) then
                u(i,j) = (ux + 0.5d0 * dt * Fx(i,j)) / rhoLocal
                v(i,j) = (uy + 0.5d0 * dt * Fy(i,j)) / rhoLocal
            else
                u(i,j) = (ux + 0.5d0 * dt * Fx(i,j)) / rhoLocal
                v(i,j) = (uy + 0.5d0 * dt * Fy(i,j)) / rhoLocal
            endif
        end do
    end do
    !$omp end parallel do
end subroutine macro_flow


subroutine collision_temp()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k
    real(kind=8) :: n(0:qt-1), neq(0:qt-1), q(0:qt-1), nPost(0:qt-1)
    real(kind=8) :: bx, by, dBx, dBy
    real(kind=8) :: sg

    if (useChaiStrict) then
        call collision_temp_strict()
        return
    endif

    sg = 1.0d0 - 0.5d0 * sT

    !$omp parallel do default(none) shared(g,gPost,u,v,T,BxPrev,ByPrev,sg) &
    !$omp private(i,j,k,n,neq,q,nPost,bx,by,dBx,dBy)
    do j = 1, ny
        do i = 1, nx
            bx = u(i,j) * T(i,j)
            by = v(i,j) * T(i,j)

            if (useGjTemp) then
                dBx = bx - BxPrev(i,j)
                dBy = by - ByPrev(i,j)
                BxPrev(i,j) = bx
                ByPrev(i,j) = by
            else
                dBx = 0.0d0
                dBy = 0.0d0
            endif

            n(0) = g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
            n(1) = g(1,i,j)-g(3,i,j)
            n(2) = g(2,i,j)-g(4,i,j)
            n(3) = -4.0d0*g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
            n(4) = g(1,i,j)-g(2,i,j)+g(3,i,j)-g(4,i,j)

            neq(0) = T(i,j)
            neq(1) = T(i,j) * u(i,j)
            neq(2) = T(i,j) * v(i,j)
            neq(3) = -2.0d0 * T(i,j) / 3.0d0
            neq(4) = 0.0d0

            q(0) = 0.0d0
            q(1) = sT
            q(2) = sT
            q(3) = 1.0d0
            q(4) = 1.0d0

            nPost(0) = n(0) - q(0) * (n(0) - neq(0))
            nPost(1) = n(1) - q(1) * (n(1) - neq(1)) + sg * dBx
            nPost(2) = n(2) - q(2) * (n(2) - neq(2)) + sg * dBy
            nPost(3) = n(3) - q(3) * (n(3) - neq(3))
            nPost(4) = n(4) - q(4) * (n(4) - neq(4))

            gPost(0,i,j) = 0.2d0*nPost(0) - 0.2d0*nPost(3)
            gPost(1,i,j) = 0.2d0*nPost(0) + 0.5d0*nPost(1) + 0.05d0*nPost(3) + 0.25d0*nPost(4)
            gPost(2,i,j) = 0.2d0*nPost(0) + 0.5d0*nPost(2) + 0.05d0*nPost(3) - 0.25d0*nPost(4)
            gPost(3,i,j) = 0.2d0*nPost(0) - 0.5d0*nPost(1) + 0.05d0*nPost(3) + 0.25d0*nPost(4)
            gPost(4,i,j) = 0.2d0*nPost(0) - 0.5d0*nPost(2) + 0.05d0*nPost(3) - 0.25d0*nPost(4)
        end do
    end do
    !$omp end parallel do
end subroutine collision_temp


subroutine collision_temp_strict()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k, l
    real(kind=8) :: geq(0:qt-1), gj(0:qt-1)
    real(kind=8) :: n(0:qt-1), neq(0:qt-1), nSrc(0:qt-1), nPost(0:qt-1)

    !$omp parallel do default(none) shared(g,gPost,u,v,T,BxPrev,ByPrev,Mtemp,MinvTemp,Stemp) &
    !$omp private(i,j,k,l,geq,gj,n,neq,nSrc,nPost)
    do j = 1, ny
        do i = 1, nx
            call compute_geq_temp(T(i,j), u(i,j), v(i,j), geq)
            call compute_temp_source(T(i,j), u(i,j), v(i,j), BxPrev(i,j), ByPrev(i,j), gj)

            do k = 0, qt-1
                n(k) = 0.0d0
                neq(k) = 0.0d0
                nSrc(k) = 0.0d0
                do l = 0, qt-1
                    n(k) = n(k) + Mtemp(k,l) * g(l,i,j)
                    neq(k) = neq(k) + Mtemp(k,l) * geq(l)
                    nSrc(k) = nSrc(k) + Mtemp(k,l) * gj(l)
                end do
                nPost(k) = n(k) - Stemp(k) * (n(k) - neq(k)) + nSrc(k)
            end do

            do k = 0, qt-1
                gPost(k,i,j) = 0.0d0
                do l = 0, qt-1
                    gPost(k,i,j) = gPost(k,i,j) + MinvTemp(k,l) * nPost(l)
                end do
            end do
        end do
    end do
    !$omp end parallel do
end subroutine collision_temp_strict


subroutine streaming_temp()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k, ii, jj

    gNext = 0.0d0

    !$omp parallel do default(none) shared(gPost,gNext,exTemp,eyTemp) private(i,j,k,ii,jj)
    do j = 1, ny
        do i = 1, nx
            do k = 0, qt-1
                ii = i + exTemp(k)
                jj = j + eyTemp(k)
                gNext(k,ii,jj) = gPost(k,i,j)
            end do
        end do
    end do
    !$omp end parallel do

    g = gNext
end subroutine streaming_temp


subroutine bounceback_temp()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !$omp parallel do default(none) shared(g,gNext,wTemp) private(j)
    do j = 1, ny
        g(1,1,j) = 2.0d0 * wTemp(1) * Thot - gNext(3,0,j)
        g(3,nx,j) = 2.0d0 * wTemp(3) * Tcold - gNext(1,nx+1,j)
    end do
    !$omp end parallel do

    !$omp parallel do default(none) shared(g,gNext) private(i)
    do i = 1, nx
        g(2,i,1) = gNext(4,i,0)
        g(4,i,ny) = gNext(2,i,ny+1)
    end do
    !$omp end parallel do
end subroutine bounceback_temp


subroutine macro_temp()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k
    real(kind=8) :: tempLocal

    !$omp parallel do default(none) shared(g,T,u,v,BxPrev,ByPrev) private(i,j,k,tempLocal)
    do j = 1, ny
        do i = 1, nx
            tempLocal = 0.0d0
            do k = 0, qt-1
                tempLocal = tempLocal + g(k,i,j)
            end do
            T(i,j) = tempLocal
        end do
    end do
    !$omp end parallel do

    !$omp parallel do default(none) shared(BxPrev,ByPrev,u,v,T) private(i,j)
    do j = 1, ny
        do i = 1, nx
            BxPrev(i,j) = u(i,j) * T(i,j)
            ByPrev(i,j) = v(i,j) * T(i,j)
        end do
    end do
    !$omp end parallel do
end subroutine macro_temp


subroutine check()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: numU, denU, numT, denT, du, dv, dtmp

    numU = 0.0d0
    denU = 0.0d0
    numT = 0.0d0
    denT = 0.0d0

    !$omp parallel do default(none) shared(u,v,T,uPrev,vPrev,TPrev) private(i,j,du,dv,dtmp) &
    !$omp reduction(+:numU,denU,numT,denT)
    do j = 1, ny
        do i = 1, nx
            du = u(i,j) - uPrev(i,j)
            dv = v(i,j) - vPrev(i,j)
            dtmp = T(i,j) - TPrev(i,j)

            numU = numU + du * du + dv * dv
            denU = denU + uPrev(i,j) * uPrev(i,j) + vPrev(i,j) * vPrev(i,j)
            numT = numT + dtmp * dtmp
            denT = denT + TPrev(i,j) * TPrev(i,j)

            uPrev(i,j) = u(i,j)
            vPrev(i,j) = v(i,j)
            TPrev(i,j) = T(i,j)
        end do
    end do
    !$omp end parallel do

    errorU = dsqrt(numU / max(denU, 1.0d-30))
    errorT = dsqrt(numT / max(denT, 1.0d-30))
end subroutine check


subroutine calNuRe()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: dx, dTdx, sumNu, sumRe, velMag

    dimensionlessTime = dimensionlessTime + 1
    if (dimensionlessTime .gt. ubound(NuVolAvg, 1)) then
        write(*,*) "Error: dimensionlessTime exceeds array bound in calNuRe"
        stop
    endif

    dx = 1.0d0 / lengthUnit
    sumNu = 0.0d0
    sumRe = 0.0d0

    !$omp parallel do default(none) shared(u,v,T,dx) private(i,j,dTdx,velMag) reduction(+:sumNu,sumRe)
    do j = 1, ny
        do i = 1, nx
            if (i == 1) then
                dTdx = (-3.0d0 * T(1,j) - T(2,j) + 4.0d0 * Thot) / (3.0d0 * dx)
            elseif (i == nx) then
                dTdx = (-4.0d0 * Tcold + 3.0d0 * T(nx,j) + T(nx-1,j)) / (3.0d0 * dx)
            else
                dTdx = (T(i-1,j) - T(i+1,j)) / (2.0d0 * dx)
            endif

            sumNu = sumNu + (lengthUnit / diffusivity) * u(i,j) * (T(i,j) - Tref) + dTdx
            velMag = dsqrt(u(i,j) * u(i,j) + v(i,j) * v(i,j))
            sumRe = sumRe + velMag * lengthUnit / viscosity
        end do
    end do
    !$omp end parallel do

    NuVolAvg(dimensionlessTime) = (sumNu / dble(nx * ny)) / (Thot - Tcold)
    ReVolAvg(dimensionlessTime) = sumRe / dble(nx * ny)

    open(unit=11, file="Nu_VolAvg.dat", status="unknown", position="append")
    write(11,'(i12,1x,es16.8)') itc, NuVolAvg(dimensionlessTime)
    close(11)

    open(unit=12, file="Re_VolAvg.dat", status="unknown", position="append")
    write(12,'(i12,1x,es16.8)') itc, ReVolAvg(dimensionlessTime)
    close(12)
end subroutine calNuRe


subroutine output_Tecplot()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k
    real(kind=4) :: zoneMarker, eohMarker
    character(len=40) :: title, v1, v2, v3, v4, v5
    character(len=40) :: zoneName
    character(len=100) :: filename
    integer(kind=4), parameter :: kmax = 1

    pltFileNum = pltFileNum + 1
    write(filename,'(i12.12)') pltFileNum
    filename = adjustl(filename)

    open(41, file=trim(pltFolderPrefix)//"-"//trim(filename)//".plt", &
        access="stream", form="unformatted")

    zoneMarker = 299.0
    eohMarker = 357.0

    write(41) "#!TDV101"
    write(41) 1

    title = "ChaiPRE2020"
    call dumpstring(title)

    write(41) 5
    v1 = "X"
    call dumpstring(v1)
    v2 = "Y"
    call dumpstring(v2)
    v3 = "U"
    call dumpstring(v3)
    v4 = "V"
    call dumpstring(v4)
    v5 = "T"
    call dumpstring(v5)

    write(41) zoneMarker
    zoneName = "ZONE 001"
    call dumpstring(zoneName)
    write(41) -1
    write(41) 0
    write(41) 1
    write(41) 0
    write(41) 0
    write(41) nx
    write(41) ny
    write(41) kmax
    write(41) 0
    write(41) eohMarker

    write(41) zoneMarker
    write(41) 2
    write(41) 2
    write(41) 2
    write(41) 2
    write(41) 2
    write(41) 0
    write(41) -1

    do k = 1, kmax
        do j = 1, ny
            do i = 1, nx
                write(41) real(xp(i), kind=8)
                write(41) real(yp(j), kind=8)
                write(41) real(u(i,j), kind=8)
                write(41) real(v(i,j), kind=8)
                write(41) real(T(i,j), kind=8)
            end do
        end do
    end do

    close(41)
end subroutine output_Tecplot


subroutine dumpstring(instring)
    implicit none
    character(len=*), intent(in) :: instring
    integer(kind=4) :: stringLength
    integer(kind=4) :: ii, codeInt

    stringLength = len_trim(instring)
    do ii = 1, stringLength
        codeInt = ichar(instring(ii:ii))
        write(41) codeInt
    end do
    write(41) 0
end subroutine dumpstring


subroutine SideHeatedcalc_Nu_global()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: dx, dTdx, qx, sum_qx, deltaT, coef

    dx = 1.0d0 / lengthUnit
    deltaT = Thot - Tcold
    coef = lengthUnit / diffusivity
    sum_qx = 0.0d0

    !$omp parallel do default(none) shared(u,T,dx,coef) private(i,j,dTdx,qx) reduction(+:sum_qx)
    do j = 1, ny
        do i = 1, nx
            if (i == 1) then
                dTdx = (-3.0d0 * T(1,j) - T(2,j) + 4.0d0 * Thot) / (3.0d0 * dx)
            elseif (i == nx) then
                dTdx = (-4.0d0 * Tcold + 3.0d0 * T(nx,j) + T(nx-1,j)) / (3.0d0 * dx)
            else
                dTdx = (T(i-1,j) - T(i+1,j)) / (2.0d0 * dx)
            endif

            qx = coef * u(i,j) * (T(i,j) - Tref) + dTdx
            sum_qx = sum_qx + qx
        end do
    end do
    !$omp end parallel do

    Nu_global = (sum_qx / dble(nx * ny)) / deltaT

    write(*,'(a,1x,es16.8)') "Nu_global =", Nu_global
    open(unit=00, file="SimulationSettings.txt", status="unknown", position="append")
    write(00,'(a,1x,es16.8)') "Nu_global =", Nu_global
    close(00)
end subroutine SideHeatedcalc_Nu_global


subroutine SideHeatedcalc_Nu_wall_avg()
    use commondata
    implicit none
    integer(kind=4) :: j, iMid, k, iL, iR
    integer(kind=4) :: jmax, jmin
    integer(kind=4) :: jj(5)
    real(kind=8) :: dx, deltaT, coef
    real(kind=8) :: qx_wall, sum_hot, sum_cold, sum_mid
    real(kind=8) :: Nu_left(1:ny), Nu_left_ext(0:ny+1)
    real(kind=8) :: T_wb, T_wt
    real(kind=8) :: yfit(4), Tfit(4)
    real(kind=8) :: yk(5), fk(5), fstar, ystar

    dx = 1.0d0 / lengthUnit
    deltaT = Thot - Tcold
    coef = lengthUnit / diffusivity

    sum_hot = 0.0d0
    !$omp parallel do default(none) shared(T,Nu_left,dx,deltaT) private(j,qx_wall) reduction(+:sum_hot)
    do j = 1, ny
        qx_wall = 2.0d0 * (Thot - T(1,j)) / dx
        Nu_left(j) = qx_wall / deltaT
        sum_hot = sum_hot + Nu_left(j)
    end do
    !$omp end parallel do
    Nu_hot = sum_hot / dble(ny)

    Nu_left_ext(1:ny) = Nu_left(1:ny)

    yfit(1) = yp(1); Tfit(1) = T(1,1)
    yfit(2) = yp(2); Tfit(2) = T(1,2)
    yfit(3) = yp(3); Tfit(3) = T(1,3)
    yfit(4) = yp(4); Tfit(4) = T(1,4)
    call fit_adiabatic_wall_T4(0.0d0, yfit, Tfit, T_wb)
    Nu_left_ext(0) = (2.0d0 * (Thot - T_wb) / dx) / deltaT

    yfit(1) = yp(ny-3); Tfit(1) = T(1,ny-3)
    yfit(2) = yp(ny-2); Tfit(2) = T(1,ny-2)
    yfit(3) = yp(ny-1); Tfit(3) = T(1,ny-1)
    yfit(4) = yp(ny);   Tfit(4) = T(1,ny)
    call fit_adiabatic_wall_T4(yp(ny+1), yfit, Tfit, T_wt)
    Nu_left_ext(ny+1) = (2.0d0 * (Thot - T_wt) / dx) / deltaT

    jmax = 0
    jmin = 0
    Nu_hot_max = Nu_left_ext(0)
    Nu_hot_min = Nu_left_ext(0)
    do j = 1, ny + 1
        if (Nu_left_ext(j) > Nu_hot_max) then
            Nu_hot_max = Nu_left_ext(j)
            jmax = j
        endif
        if (Nu_left_ext(j) < Nu_hot_min) then
            Nu_hot_min = Nu_left_ext(j)
            jmin = j
        endif
    end do

    if (jmax <= 2) then
        jj = (/ 0, 1, 2, 3, 4 /)
    elseif (jmax >= ny - 1) then
        jj = (/ ny-3, ny-2, ny-1, ny, ny+1 /)
    else
        jj = (/ jmax-2, jmax-1, jmax, jmax+1, jmax+2 /)
    endif
    do k = 1, 5
        yk(k) = yp(jj(k))
        fk(k) = Nu_left_ext(jj(k))
    end do
    call fit_parabola_ls5(yk, fk, +1, fstar, ystar)
    Nu_hot_max = fstar
    Nu_hot_max_position = ystar

    if (jmin >= 4) then
        jj = (/ jmin-4, jmin-3, jmin-2, jmin-1, jmin /)
    else
        jj = (/ 0, 1, 2, 3, 4 /)
    endif
    do k = 1, 5
        yk(k) = yp(jj(k))
        fk(k) = Nu_left_ext(jj(k))
    end do
    call fit_parabola_ls5(yk, fk, -1, fstar, ystar)
    Nu_hot_min = fstar
    Nu_hot_min_position = ystar

    sum_cold = 0.0d0
    !$omp parallel do default(none) shared(T,dx,deltaT) private(j,qx_wall) reduction(+:sum_cold)
    do j = 1, ny
        qx_wall = 2.0d0 * (T(nx,j) - Tcold) / dx
        sum_cold = sum_cold + qx_wall / deltaT
    end do
    !$omp end parallel do
    Nu_cold = sum_cold / dble(ny)

    sum_mid = 0.0d0
    if (mod(nx,2) == 1) then
        iMid = int(0.5d0 * (dble(nx) + 1.0d0))
        !$omp parallel do default(none) shared(u,T,iMid,dx,deltaT,coef) private(j) reduction(+:sum_mid)
        do j = 1, ny
            sum_mid = sum_mid + (coef * u(iMid,j) * (T(iMid,j) - Tref) + &
                (T(iMid-1,j) - T(iMid+1,j)) / (2.0d0 * dx)) / deltaT
        end do
        !$omp end parallel do
    else
        iL = int(0.5d0 * dble(nx))
        iR = iL + 1
        !$omp parallel do default(none) shared(u,T,iL,iR,dx,deltaT,coef) private(j) reduction(+:sum_mid)
        do j = 1, ny
            sum_mid = sum_mid + (coef * 0.5d0 * (u(iL,j) * (T(iL,j) - Tref) + u(iR,j) * (T(iR,j) - Tref)) + &
                (T(iL,j) - T(iR,j)) / dx) / deltaT
        end do
        !$omp end parallel do
    endif
    Nu_middle = sum_mid / dble(ny)

    write(*,'(a,1x,es16.8)') "Nu_hot    =", Nu_hot
    write(*,'(a,1x,es16.8)') "Nu_cold   =", Nu_cold
    write(*,'(a,1x,es16.8)') "Nu_middle =", Nu_middle
    write(*,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_max =", Nu_hot_max, "y_max =", Nu_hot_max_position
    write(*,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_min =", Nu_hot_min, "y_min =", Nu_hot_min_position

    open(unit=00, file="SimulationSettings.txt", status="unknown", position="append")
    write(00,'(a,1x,es16.8)') "Nu_hot    =", Nu_hot
    write(00,'(a,1x,es16.8)') "Nu_cold   =", Nu_cold
    write(00,'(a,1x,es16.8)') "Nu_middle =", Nu_middle
    write(00,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_max =", Nu_hot_max, "y_max =", Nu_hot_max_position
    write(00,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_min =", Nu_hot_min, "y_min =", Nu_hot_min_position
    close(00)
end subroutine SideHeatedcalc_Nu_wall_avg


subroutine fit_adiabatic_wall_T4(y0, y, tt, T_wall)
    implicit none
    real(kind=8), intent(in) :: y0
    real(kind=8), intent(in) :: y(4), tt(4)
    real(kind=8), intent(out) :: T_wall
    real(kind=8) :: s(4), s0, s1, s2, b0, b1, det
    integer(kind=4) :: k

    do k = 1, 4
        s(k) = (y(k) - y0) * (y(k) - y0)
    end do

    s0 = 4.0d0
    s1 = 0.0d0
    s2 = 0.0d0
    b0 = 0.0d0
    b1 = 0.0d0
    do k = 1, 4
        s1 = s1 + s(k)
        s2 = s2 + s(k) * s(k)
        b0 = b0 + tt(k)
        b1 = b1 + tt(k) * s(k)
    end do

    det = s0 * s2 - s1 * s1
    T_wall = (b0 * s2 - b1 * s1) / det
end subroutine fit_adiabatic_wall_T4


subroutine fit_parabola_ls5(y, f, mode, fstar, ystar)
    implicit none
    real(kind=8), intent(in) :: y(5), f(5)
    integer(kind=4), intent(in) :: mode
    real(kind=8), intent(out) :: fstar, ystar
    real(kind=8) :: s0, s1, s2, s3, s4, f0, f1, f2
    real(kind=8) :: det, da, db, dc, a, b, c, ymin, ymax
    integer(kind=4) :: k, kbest
    real(kind=8), parameter :: epsD = 1.0d-20, epsA = 1.0d-14

    kbest = 1
    do k = 2, 5
        if (mode == 1) then
            if (f(k) > f(kbest)) kbest = k
        else
            if (f(k) < f(kbest)) kbest = k
        endif
    end do

    s0 = 0.0d0; s1 = 0.0d0; s2 = 0.0d0; s3 = 0.0d0; s4 = 0.0d0
    f0 = 0.0d0; f1 = 0.0d0; f2 = 0.0d0
    do k = 1, 5
        s0 = s0 + 1.0d0
        s1 = s1 + y(k)
        s2 = s2 + y(k) * y(k)
        s3 = s3 + y(k) * y(k) * y(k)
        s4 = s4 + y(k) * y(k) * y(k) * y(k)
        f0 = f0 + f(k)
        f1 = f1 + y(k) * f(k)
        f2 = f2 + y(k) * y(k) * f(k)
    end do

    det = s4 * (s2 * s0 - s1 * s1) - s3 * (s3 * s0 - s1 * s2) + s2 * (s3 * s1 - s2 * s2)
    da = f2 * (s2 * s0 - s1 * s1) - s3 * (f1 * s0 - s1 * f0) + s2 * (f1 * s1 - s2 * f0)
    db = s4 * (f1 * s0 - s1 * f0) - f2 * (s3 * s0 - s1 * s2) + s2 * (s3 * f0 - f1 * s2)
    dc = s4 * (s2 * f0 - f1 * s1) - s3 * (s3 * f0 - f1 * s2) + f2 * (s3 * s1 - s2 * s2)

    if (dabs(det) > epsD) then
        a = da / det
        b = db / det
        c = dc / det
        if (dabs(a) > epsA) then
            ystar = -b / (2.0d0 * a)
            ymin = minval(y)
            ymax = maxval(y)
            if ((ystar >= ymin) .and. (ystar <= ymax)) then
                fstar = c - b * b / (4.0d0 * a)
            else
                ystar = y(kbest)
                fstar = f(kbest)
            endif
        else
            ystar = y(kbest)
            fstar = f(kbest)
        endif
    else
        ystar = y(kbest)
        fstar = f(kbest)
    endif
end subroutine fit_parabola_ls5


subroutine SideHeatedcalc_umid_max()
    use commondata
    implicit none
    integer(kind=4) :: j, k, iMid, iL, iR, j0, jj(5)
    real(kind=8) :: uline(1:ny), s(5), fu(5), umaxFit, yFit, xmid, coef, tmp

    coef = lengthUnit / diffusivity

    if (mod(nx,2) == 1) then
        iMid = int(0.5d0 * (dble(nx) + 1.0d0))
        xmid = xp(iMid)
        do j = 1, ny
            uline(j) = u(iMid,j)
        end do
    else
        iL = int(0.5d0 * dble(nx))
        iR = iL + 1
        xmid = 0.5d0 * (xp(iL) + xp(iR))
        do j = 1, ny
            uline(j) = 0.5d0 * (u(iL,j) + u(iR,j))
        end do
    endif

    j0 = 1
    tmp = uline(1)
    do j = 2, ny
        if (uline(j) > tmp) then
            tmp = uline(j)
            j0 = j
        endif
    end do

    if (j0 <= 2) then
        jj = (/ 1, 2, 3, 4, 5 /)
    elseif (j0 >= ny - 1) then
        jj = (/ ny-4, ny-3, ny-2, ny-1, ny /)
    else
        jj = (/ j0-2, j0-1, j0, j0+1, j0+2 /)
    endif

    do k = 1, 5
        s(k) = yp(jj(k))
        fu(k) = uline(jj(k))
    end do

    call fit_parabola_ls5(s, fu, +1, umaxFit, yFit)

    write(*,'(a,1x,es16.8,1x,a,1x,es16.8,1x,a,1x,es16.8)') &
        "u_mid_max =", umaxFit * coef, "at y =", yFit, "on x_mid =", xmid
    open(unit=00, file="SimulationSettings.txt", status="unknown", position="append")
    write(00,'(a,1x,es16.8,1x,a,1x,es16.8,1x,a,1x,es16.8)') &
        "u_mid_max =", umaxFit * coef, "at y =", yFit, "on x_mid =", xmid
    close(00)
end subroutine SideHeatedcalc_umid_max


subroutine SideHeatedcalc_vmid_max()
    use commondata
    implicit none
    integer(kind=4) :: i, k, jMid, jB, jT, i0, ii(5)
    real(kind=8) :: vline(1:nx), s(5), fv(5), vmaxFit, xFit, ymid, coef, tmp

    coef = lengthUnit / diffusivity

    if (mod(ny,2) == 1) then
        jMid = int(0.5d0 * (dble(ny) + 1.0d0))
        ymid = yp(jMid)
        do i = 1, nx
            vline(i) = v(i,jMid)
        end do
    else
        jB = int(0.5d0 * dble(ny))
        jT = jB + 1
        ymid = 0.5d0 * (yp(jB) + yp(jT))
        do i = 1, nx
            vline(i) = 0.5d0 * (v(i,jB) + v(i,jT))
        end do
    endif

    i0 = 1
    tmp = vline(1)
    do i = 2, nx
        if (vline(i) > tmp) then
            tmp = vline(i)
            i0 = i
        endif
    end do

    if (i0 <= 2) then
        ii = (/ 1, 2, 3, 4, 5 /)
    elseif (i0 >= nx - 1) then
        ii = (/ nx-4, nx-3, nx-2, nx-1, nx /)
    else
        ii = (/ i0-2, i0-1, i0, i0+1, i0+2 /)
    endif

    do k = 1, 5
        s(k) = xp(ii(k))
        fv(k) = vline(ii(k))
    end do

    call fit_parabola_ls5(s, fv, +1, vmaxFit, xFit)

    write(*,'(a,1x,es16.8,1x,a,1x,es16.8,1x,a,1x,es16.8)') &
        "v_mid_max =", vmaxFit * coef, "at x =", xFit, "on y_mid =", ymid
    open(unit=00, file="SimulationSettings.txt", status="unknown", position="append")
    write(00,'(a,1x,es16.8,1x,a,1x,es16.8,1x,a,1x,es16.8)') &
        "v_mid_max =", vmaxFit * coef, "at x =", xFit, "on y_mid =", ymid
    close(00)
end subroutine SideHeatedcalc_vmid_max


subroutine output_centerline_profiles()
    use commondata
    implicit none
    integer(kind=4) :: i, j, iMid, jMid, iL, iR, jB, jT
    real(kind=8) :: uMid, vMid

    open(unit=21, file="centerline_u_y.dat", status="unknown")
    if (mod(nx,2) == 1) then
        iMid = int(0.5d0 * (dble(nx) + 1.0d0))
        do j = 1, ny
            write(21,'(2(1x,es20.12))') yp(j), u(iMid,j) * lengthUnit / diffusivity
        end do
    else
        iL = int(0.5d0 * dble(nx))
        iR = iL + 1
        do j = 1, ny
            uMid = 0.5d0 * (u(iL,j) + u(iR,j)) * lengthUnit / diffusivity
            write(21,'(2(1x,es20.12))') yp(j), uMid
        end do
    endif
    close(21)

    open(unit=22, file="centerline_v_x.dat", status="unknown")
    if (mod(ny,2) == 1) then
        jMid = int(0.5d0 * (dble(ny) + 1.0d0))
        do i = 1, nx
            write(22,'(2(1x,es20.12))') xp(i), v(i,jMid) * lengthUnit / diffusivity
        end do
    else
        jB = int(0.5d0 * dble(ny))
        jT = jB + 1
        do i = 1, nx
            vMid = 0.5d0 * (v(i,jB) + v(i,jT)) * lengthUnit / diffusivity
            write(22,'(2(1x,es20.12))') xp(i), vMid
        end do
    endif
    close(22)
end subroutine output_centerline_profiles


subroutine compute_feq_flow(rhoLocal, ux, uy, feq)
    use commondata
    implicit none
    real(kind=8), intent(in) :: rhoLocal, ux, uy
    real(kind=8), intent(out) :: feq(0:qf-1)
    integer(kind=4) :: k
    real(kind=8) :: cu, uu

    uu = ux * ux + uy * uy
    do k = 0, qf-1
        cu = dble(exFlow(k)) * ux + dble(eyFlow(k)) * uy
        feq(k) = wFlow(k) * rhoLocal * (1.0d0 + cu / cs2Flow + &
            0.5d0 * cu * cu / cs4Flow - 0.5d0 * uu / cs2Flow)
    end do
end subroutine compute_feq_flow


subroutine compute_geq_temp(tempLocal, ux, uy, geq)
    use commondata
    implicit none
    real(kind=8), intent(in) :: tempLocal, ux, uy
    real(kind=8), intent(out) :: geq(0:qt-1)
    integer(kind=4) :: k
    real(kind=8) :: cu

    do k = 0, qt-1
        cu = dble(exTemp(k)) * ux + dble(eyTemp(k)) * uy
        geq(k) = wTemp(k) * (tempLocal + tempLocal * cu / cs2Temp)
    end do
end subroutine compute_geq_temp


subroutine compute_flow_source(ux, uy, fxBar, fyBar, src, gaux)
    use commondata
    implicit none
    real(kind=8), intent(in) :: ux, uy, fxBar, fyBar
    real(kind=8), intent(out) :: src(0:qf-1), gaux(0:qf-1)
    integer(kind=4) :: k
    real(kind=8) :: ccx, ccy, qxx, qyy, qxy
    real(kind=8) :: wxx, wyy, wxy, trw, devxx, devyy
    real(kind=8) :: mxx, myy, mxy

    do k = 0, qf-1
        src(k) = wFlow(k) * (dble(exFlow(k)) * fxBar + dble(eyFlow(k)) * fyBar) / cs2Flow
    end do

    wxx = 2.0d0 * fxBar * ux
    wyy = 2.0d0 * fyBar * uy
    wxy = fxBar * uy + ux * fyBar
    trw = wxx + wyy
    devxx = wxx - 0.5d0 * trw
    devyy = wyy - 0.5d0 * trw

    mxx = 0.5d0 * (1.0d0 - 0.5d0 * S2b) * trw + (1.0d0 - 0.5d0 * S2s) * devxx
    myy = 0.5d0 * (1.0d0 - 0.5d0 * S2b) * trw + (1.0d0 - 0.5d0 * S2s) * devyy
    mxy = (1.0d0 - 0.5d0 * S2s) * wxy

    do k = 0, qf-1
        ccx = dble(exFlow(k))
        ccy = dble(eyFlow(k))
        qxx = ccx * ccx - cs2Flow
        qyy = ccy * ccy - cs2Flow
        qxy = ccx * ccy
        gaux(k) = wFlow(k) * (mxx * qxx + 2.0d0 * mxy * qxy + myy * qyy) / &
            (2.0d0 * cs4Flow)
    end do
end subroutine compute_flow_source


subroutine compute_temp_source(tempLocal, ux, uy, bxOld, byOld, gj)
    use commondata
    implicit none
    real(kind=8), intent(in) :: tempLocal, ux, uy, bxOld, byOld
    real(kind=8), intent(out) :: gj(0:qt-1)
    integer(kind=4) :: k
    real(kind=8) :: bx, by, dBxdt, dBydt, m1gx, m1gy

    if (.not. useGjTemp) then
        gj = 0.0d0
        return
    endif

    bx = ux * tempLocal
    by = uy * tempLocal
    dBxdt = (bx - bxOld) / dt
    dBydt = (by - byOld) / dt
    m1gx = (1.0d0 - 0.5d0 * sT) * dBxdt
    m1gy = (1.0d0 - 0.5d0 * sT) * dBydt

    do k = 0, qt-1
        gj(k) = wTemp(k) * (dble(exTemp(k)) * m1gx + dble(eyTemp(k)) * m1gy) / cs2Temp
    end do
end subroutine compute_temp_source


subroutine invert_matrix(n, aIn, aOut)
    implicit none
    integer(kind=4), intent(in) :: n
    real(kind=8), intent(in) :: aIn(0:n-1,0:n-1)
    real(kind=8), intent(out) :: aOut(0:n-1,0:n-1)
    integer(kind=4) :: i, j
    integer(kind=4) :: pivotRow
    real(kind=8) :: aug(1:n,1:2*n)
    real(kind=8) :: pivotValue, factor, maxAbs

    aug = 0.0d0
    do i = 1, n
        do j = 1, n
            aug(i,j) = aIn(i-1,j-1)
        end do
        aug(i,n+i) = 1.0d0
    end do

    do i = 1, n
        pivotRow = i
        maxAbs = dabs(aug(i,i))
        do j = i + 1, n
            if (dabs(aug(j,i)) > maxAbs) then
                maxAbs = dabs(aug(j,i))
                pivotRow = j
            endif
        end do

        if (maxAbs <= 1.0d-30) then
            write(*,*) "Error: singular matrix in invert_matrix"
            stop
        endif

        if (pivotRow /= i) call swap_rows(n, 2*n, aug, i, pivotRow)

        pivotValue = aug(i,i)
        do j = 1, 2 * n
            aug(i,j) = aug(i,j) / pivotValue
        end do

        do pivotRow = 1, n
            if (pivotRow /= i) then
                factor = aug(pivotRow,i)
                do j = 1, 2 * n
                    aug(pivotRow,j) = aug(pivotRow,j) - factor * aug(i,j)
                end do
            endif
        end do
    end do

    do i = 1, n
        do j = 1, n
            aOut(i-1,j-1) = aug(i,n+j)
        end do
    end do
end subroutine invert_matrix


subroutine swap_rows(nrow, ncol, a, row1, row2)
    implicit none
    integer(kind=4), intent(in) :: nrow, ncol, row1, row2
    real(kind=8), intent(inout) :: a(1:nrow,1:ncol)
    integer(kind=4) :: j
    real(kind=8) :: temp

    do j = 1, ncol
        temp = a(row1,j)
        a(row1,j) = a(row2,j)
        a(row2,j) = temp
    end do
end subroutine swap_rows
