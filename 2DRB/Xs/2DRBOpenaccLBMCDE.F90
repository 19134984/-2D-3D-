!===================================================================================================
!!! 二维侧壁差温自然对流：LBM-CDE 数值测试程序（OpenACC 单 GPU）
!!!
!!! 论文来源：pdf/LBM-CDE.pdf
!!! 离散模型：流场 D2Q9-BGK + 温度场 D2Q9-BGK
!!! 主要论文对应：
!!!   1) 流场平衡态与源项：式 (4)、式 (26)，有效运动黏度见式 (27)；
!!!   2) 温度平衡态与源项：式 (16)、式 (28)，有效扩散率见式 (29)；
!!!   3) 局部应变率和温度梯度：式 (31)-(35)；
!!!   4) 无滑移、恒温和绝热边界：式 (36)、式 (38)、式 (39)；
!!!   5) 侧热方腔算例及 Nu 定义：式 (50)-(52)。
!!!
!!! 程序定位：用于检查论文方法和做小规模数值测试，不包含重启、MPI 或复杂后处理。
!!! 默认网格为 128x128，便于快速测试；论文中 Ra=1e5/1e6 使用 256x256，Ra=1e7 使用 512x512。
!!! 若要严格复现论文表 2，应同时匹配网格、tau_fL、统计方式和收敛状态，而不能只修改 Ra。
!===================================================================================================


!===================================================================================================
! 算例开关
! steadyFlow：达到 errorU/errorT 阈值后停止；unsteadyFlow：只按 itc_max 停止。
#define steadyFlow
!#define unsteadyFlow

#if defined(steadyFlow) && defined(unsteadyFlow)
#error "Choose only one flow mode: steadyFlow or unsteadyFlow"
#endif
#if !defined(steadyFlow) && !defined(unsteadyFlow)
#error "Define one flow mode: steadyFlow or unsteadyFlow"
#endif

! 本程序固定为封闭方腔：四壁速度均采用 halfway bounce-back 无滑移边界。
#define HorizontalWallsNoslip
#define VerticalWallsNoslip

! 温度边界固定为：左热、右冷、上下绝热。
#define SideHeatedCell
#define HorizontalWallsAdiabatic
#define VerticalWallsConstT

! 当前实现没有其它边界分支，下面的检查用于防止以后改宏时得到一个“能编译但边界不完整”的程序。
#if !defined(HorizontalWallsNoslip) || !defined(VerticalWallsNoslip)
#error "This test requires no-slip velocity boundaries on all four walls"
#endif
#if !defined(SideHeatedCell) || !defined(HorizontalWallsAdiabatic) || !defined(VerticalWallsConstT)
#error "This test requires the side-heated cavity temperature boundaries"
#endif

! 编译时可用 -DNX_OVERRIDE=256 -DNY_OVERRIDE=256 改网格，而不必直接改源文件。
#ifndef NX_OVERRIDE
#define NX_OVERRIDE 128
#endif
#ifndef NY_OVERRIDE
#define NY_OVERRIDE 128
#endif
#if NX_OVERRIDE != NY_OVERRIDE
#error "The paper test is a square cavity and requires NX_OVERRIDE == NY_OVERRIDE"
#endif
#ifndef ITC_MAX_OVERRIDE
#define ITC_MAX_OVERRIDE 20000000
#endif
#ifndef RAYLEIGH_OVERRIDE
#define RAYLEIGH_OVERRIDE 1.0d5
#endif
!===================================================================================================


!===================================================================================================
! 全局数据模块
    module commondata
        implicit none

        !-------------------------------------------------------------------------------------------
        ! 网格与 D2Q9 常量
        integer(kind=4), parameter :: nx=NX_OVERRIDE, ny=NY_OVERRIDE ! 流体节点数；壁面位于首末节点外半格
        real(kind=8), parameter :: rho0=1.0d0                        ! 弱可压缩模型的参考密度
        real(kind=8), parameter :: cs2=1.0d0/3.0d0                   ! D2Q9 格子声速平方

#ifdef SideHeatedCell
        real(kind=8), parameter :: lengthUnit=dble(nx)               ! 侧热方腔的特征长度 L
#else
        real(kind=8), parameter :: lengthUnit=dble(ny)
#endif

        !-------------------------------------------------------------------------------------------
        ! 无量纲工况参数
        real(kind=8), parameter :: Rayleigh=RAYLEIGH_OVERRIDE        ! Ra，论文式 (51)
        real(kind=8), parameter :: Prandtl=0.71d0                    ! Pr，论文表 2 使用 0.71
        real(kind=8), parameter :: Mach=0.1d0                        ! 用于确定格子黏度的数值 Mach 数
        real(kind=8), parameter :: Thot=0.5d0, Tcold=-0.5d0          ! 左热壁、右冷壁温度
        real(kind=8), parameter :: Tref=0.5d0*(Thot+Tcold)           ! Boussinesq 参考温度 T0
        real(kind=8), parameter :: deltaT=Thot-Tcold                 ! 温差 DeltaT

        !-------------------------------------------------------------------------------------------
        ! LBM-CDE 自由参数
        ! chi_s 调节剪切黏度，chi_b 调节体黏度，chi_kappa 调节温度扩散率。
        ! chi=0 退化为标准 BGK 输运系数；chi 必须小于 1，才能保持正黏度/正扩散率。
        ! 论文自然对流算例比较了 chi=0 和 chi=0.8，并在所有算例中令 chi_b=0。
        real(kind=8), parameter :: chi=0.8d0
        real(kind=8), parameter :: chi_s=chi
        real(kind=8), parameter :: chi_b=0.0d0
        real(kind=8), parameter :: chi_kappa=chi
        real(kind=8), parameter :: heatSourceQ=0.0d0                  ! 方腔无体热源；保留 Q 以实现完整式 (28)/(35)

        !-------------------------------------------------------------------------------------------
        ! 格子单位下的目标输运系数
        ! 这一写法沿用仓库均匀网格/ISLBM 程序：先由 Ra、Pr、Ma 确定 nu，再令 kappa=nu/Pr。
        real(kind=8), parameter :: viscosity=Mach*lengthUnit*dsqrt(Prandtl/(3.0d0*Rayleigh))
        real(kind=8), parameter :: diffusivity=viscosity/Prandtl

        ! 论文式 (27)/(29)，dt=1：
        !   nu    = (tau_fL-0.5)*(1-chi_s)*cs2
        !   kappa = (tau_gL-0.5)*(1-chi_kappa)*cs2
        ! 因此 chi 接近 1 时，即使 nu/kappa 很小，松弛时间仍可远离 0.5，从而改善稳定性。
        real(kind=8), parameter :: taufL=0.5d0+viscosity/(cs2*(1.0d0-chi_s))
        real(kind=8), parameter :: taugL=0.5d0+diffusivity/(cs2*(1.0d0-chi_kappa))
        real(kind=8), parameter :: bulkViscosity=(taufL-0.5d0)*(1.0d0-chi_b)*cs2 ! D=2 时的 nu_B
        real(kind=8), parameter :: omegaF=1.0d0/taufL
        real(kind=8), parameter :: omegaG=1.0d0/taugL

        !-------------------------------------------------------------------------------------------
        ! Boussinesq 浮力与自由落体标度，见论文式 (50)/(51)。
        ! Fy>0 表示热流体向 +y 上升；gBeta 在这里代表 g*beta。
        real(kind=8), parameter :: gBeta=Rayleigh*viscosity*diffusivity/(deltaT*lengthUnit**3)
        real(kind=8), parameter :: timeUnit=dsqrt(lengthUnit/(gBeta*deltaT))
        real(kind=8), parameter :: velocityUnit=dsqrt(gBeta*deltaT*lengthUnit)

        !-------------------------------------------------------------------------------------------
        ! 推进与收敛设置。误差定义见 check()，它不是论文模型的一部分。
#ifdef steadyFlow
        real(kind=8), parameter :: epsU=1.0d-7, epsT=1.0d-7
        integer(kind=4), parameter :: itc_max=ITC_MAX_OVERRIDE
        integer(kind=4), parameter :: checkInterval=1000
#endif
#ifdef unsteadyFlow
        integer(kind=4), parameter :: itc_max=ITC_MAX_OVERRIDE
        integer(kind=4), parameter :: checkInterval=1000
#endif

        !-------------------------------------------------------------------------------------------
        ! 宏观场和论文式 (31)-(35) 所需的局部辅助量
        real(kind=8), allocatable :: u(:,:), v(:,:), T(:,:), rho(:,:), pressure(:,:)
        real(kind=8), allocatable :: up(:,:), vp(:,:), Tp(:,:)       ! 上一次检查时的宏观场，只在 CPU 端使用
        real(kind=8), allocatable :: Fx(:,:), Fy(:,:)                ! 力密度 F
        real(kind=8), allocatable :: gradTx(:,:), gradTy(:,:)        ! 由非平衡矩直接得到的温度梯度
        real(kind=8), allocatable :: Sxx(:,:), Sxy(:,:), Syy(:,:), Sdiv(:,:) ! 应变率与速度散度
        real(kind=8), allocatable :: f(:,:,:), f_post(:,:,:)         ! 流场迁移后/碰撞后分布函数
        real(kind=8), allocatable :: g(:,:,:), g_post(:,:,:)         ! 温度场迁移后/碰撞后分布函数
        real(kind=8) :: xp(0:nx+1), yp(0:ny+1)                       ! 归一化坐标，含物理壁面位置

        !-------------------------------------------------------------------------------------------
        ! 收敛量和后处理量
        real(kind=8) :: errorU, errorT
        real(kind=8) :: Nu_global, Nu_hot, Nu_cold, Nu_middle
        real(kind=8) :: ReVolAvg, umax_mid, vmax_mid, Numax_mid
        real(kind=8) :: y_umax_mid, x_vmax_mid, y_Numax_mid
        integer(kind=4) :: itc

        character(len=100) :: settingsFile="SimulationSettings2DOpenaccLBMCDE.txt"
        character(len=100) :: pltFile="buoyancyCavity2DOpenaccLBMCDE.plt"
        character(len=100) :: historyFile="NuRe_2DOpenaccLBMCDE.dat"

    contains

        !-------------------------------------------------------------------------------------------
        ! OpenACC 设备端 D2Q9 查询函数
        !
        ! 这里不在 GPU 核函数中直接引用模块 parameter 数组，而使用顺序设备函数返回 e_i、omega_i、lambda_i。
        ! 这样可避免不同 OpenACC 编译器对模块常量数组设备可见性处理不一致；修改速度编号时必须同步修改下列函数。
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

        ! D2Q9 权重：静止方向 4/9，轴向 1/9，对角方向 1/36。
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

        ! 论文式 (16) 的 lambda_i：非零方向等于 omega_i，仅把零速度权重改为 -5/9。
        ! 因此 sum(lambda_i)=0、一阶矩为 0、二阶矩为 cs2*delta_ab，压力只进入二阶矩。
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

        !-------------------------------------------------------------------------------------------
        ! 流场平衡分布函数，离散形式对应论文式 (4)。
        ! f 的零阶矩是密度扰动 deltaRho=rho-rho0，所以平衡态首项不是完整 rho。
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

        !-------------------------------------------------------------------------------------------
        ! 温度平衡分布函数，对应论文式 (16)。
        ! 第一项是二阶 Hermite 平衡态；第二项 lambda_i*T*p/(rho0*cs2) 建立压力与温度方程的耦合。
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
!===================================================================================================


!===================================================================================================
! 主程序：保持“流场一步 -> 温度场一步 -> 定期同步和诊断”的推进顺序。
    program main
        use openacc
        use commondata
        implicit none
        real(kind=8) :: timeStart, timeEnd, wallSeconds
        integer(kind=8) :: wallClockStart, wallClockEnd, wallClockRate
        integer(kind=4) :: numAccDevices

        ! 初始化当前 OpenACC 默认设备。实际设备由编译器和 ACC_DEVICE_TYPE/ACC_DEVICE_NUM 等环境设置决定。
        call acc_init(acc_device_default)
        numAccDevices = acc_get_num_devices(acc_device_default)
        call initial(numAccDevices)

        ! 设备数据只在这里建立一次，并在整个时间推进期间常驻，避免每个格子步重复传输大数组。
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

            !---------------------------------------------------------------------------------------
            ! 流场推进：F/p -> 局部应变率 -> 碰撞 -> 拉式迁移 -> 无滑移边界 -> 宏观量。
            ! 应变率必须在碰撞前更新，因为论文式 (26) 的 chi_s/chi_b 修正项显式依赖 S_ab。
            call compute_force()
            call compute_strain_rate()
            call collision()
            call streaming()
            call bounceback()
            call macro()

            !---------------------------------------------------------------------------------------
            ! 温度场推进：先用新流场重算 F/p 和 grad(T)，再执行温度碰撞、迁移、边界和宏观量。
            ! 论文式 (28) 同时依赖压力、力和温度梯度，所以这些量不能沿用更早的旧值。
            call compute_force()
            call compute_temperature_gradient()
            call collisionT()
            call streamingT()
            call bouncebackT()
            call macroT()

            if(mod(itc, checkInterval).EQ.0) then
                ! 先在设备端刷新诊断量，再只同步 CPU 后处理所需的宏观数组。
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
        wallSeconds = dble(wallClockEnd-wallClockStart)/dble(max(wallClockRate,1_8))
        call final_log(timeEnd-timeStart, wallSeconds)

        ! 所有主机输出结束后，再释放设备端和主机端数组。
        call exit_data_2d_openacc()
        call finalize_arrays()
    end program main
!===================================================================================================


!===================================================================================================
! 子程序: initial
! 作用: 分配数组，建立半格点均匀网格，以静止导热解初始化宏观场和两个平衡分布函数。
! 说明: f_post/g_post 在设备端用 create 建立，第一次碰撞会先清零，因此无须主机初始化。
!===================================================================================================
    subroutine initial(numAccDevices)
        use commondata
        implicit none
        integer(kind=4), intent(in) :: numAccDevices
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

        ! 流体节点位于 (i-1/2,j-1/2)，物理壁面位于 x/y=0 和 x/y=L；
        ! 这与后面的 halfway bounce-back / anti-bounce-back 边界位置一致。
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

        ! 初始温度取纯导热线性分布：左壁 Thot，右壁 Tcold；初始速度为零。
        do j = 1, ny
            do i = 1, nx
                T(i,j) = Thot + (xp(i)-xp(0))/(xp(nx+1)-xp(0))*(Tcold-Thot)
            enddo
        enddo

        f = 0.0d0
        g = 0.0d0
        ! 由初始宏观量建立平衡分布，避免第一步人为引入额外非平衡扰动。
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
        write(00,*) "Visible OpenACC default devices =", numAccDevices
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
        write(00,*) "Paper mesh reminder: Ra=1e5/1e6 -> 256x256; Ra=1e7 -> 512x512"
        write(00,*) "This run uses the Ma-derived tau values printed above, not the paper's prescribed tau_fL table."
        close(00)

        open(unit=01, file=trim(historyFile), status='replace', action='write')
        write(01,'(A)') "# itc time_ff errorU errorT Nu_hot Nu_cold Nu_middle ReVolAvg"
        close(01)
    end subroutine initial
!===================================================================================================


!===================================================================================================
! 子程序: enter_data_2d_openacc
! 作用: 建立一次性的设备数据生命周期。
! 说明: 初值已存在的数组用 copyin；每一步都会完整覆盖的碰撞后数组用 create，避免无意义传输。
!===================================================================================================
    subroutine enter_data_2d_openacc()
        use commondata
        implicit none

        !$acc enter data copyin(u,v,T,rho,pressure,Fx,Fy,gradTx,gradTy)
        !$acc enter data copyin(Sxx,Sxy,Syy,Sdiv,f,g)
        !$acc enter data create(f_post,g_post)
    end subroutine enter_data_2d_openacc
!===================================================================================================


!===================================================================================================
! 子程序: update_host_monitor_2d_openacc
! 作用: 把收敛判断、Nu/Re 和 Tecplot 输出需要的宏观量从设备同步到主机。
! 说明: f/g、应变率和力不参与当前 CPU 后处理，所以不传回，减少检查时的数据流量。
!===================================================================================================
    subroutine update_host_monitor_2d_openacc()
        use commondata
        implicit none

        !$acc update self(u,v,T,rho,pressure,gradTx,gradTy)
    end subroutine update_host_monitor_2d_openacc
!===================================================================================================


!===================================================================================================
! 子程序: exit_data_2d_openacc
! 作用: 计算和最终主机输出完成后，释放本程序建立的全部设备数据。
!===================================================================================================
    subroutine exit_data_2d_openacc()
        use commondata
        implicit none

        !$acc exit data delete(u,v,T,rho,pressure,Fx,Fy,gradTx,gradTy)
        !$acc exit data delete(Sxx,Sxy,Syy,Sdiv,f,f_post,g,g_post)
    end subroutine exit_data_2d_openacc
!===================================================================================================


!===================================================================================================
! 子程序: compute_force
! 作用: 由当前 rho/T 计算弱可压缩压力 p 和 Boussinesq 力 F，供两个分布函数共同使用。
! 论文对应: p=cs2*(rho-rho0)，Fy=rho0*g*beta*(T-Tref)，见式 (50)。
!===================================================================================================
    subroutine compute_force()
        use commondata
        implicit none
        integer(kind=4) :: i, j

        !$acc parallel loop collapse(2) present(Fx,Fy,T,rho,pressure) private(i,j)
        do j = 1, ny
            do i = 1, nx
                pressure(i,j) = cs2*(rho(i,j)-rho0)       ! He-Luo 弱可压缩状态方程
                Fx(i,j) = 0.0d0                           ! 本算例没有水平方向体力
                Fy(i,j) = rho0*gBeta*(T(i,j)-Tref)        ! 热流体向 +y 方向上升
            enddo
        enddo
    end subroutine compute_force
!===================================================================================================


!===================================================================================================
! 子程序: compute_strain_rate
! 作用: 由流场局部非平衡二阶矩 f-feq 计算 Sxx、Sxy、Syy 和 Sdiv，不使用有限差分。
! 论文对应: 非对角分量用式 (31)，速度散度用式 (32)，对角分量用式 (33)，本程序 D=2、dt=1。
! 说明: 论文式中使用动力黏度 mu=rho0*nu，本程序显式换算后再构造各分母。
!===================================================================================================
    subroutine compute_strain_rate()
        use commondata
        implicit none
        integer(kind=4) :: i, j, alpha
        real(kind=8) :: feqLoc, neqxx, neqxy, neqyy
        real(kind=8) :: muShear, muBulk, denomDiag, denomShear, denomDiv, coeffTrace

        ! 下列分母分别来自式 (33)、式 (31) 和式 (32)。
        muShear = rho0*viscosity
        muBulk = rho0*bulkViscosity
        denomDiag = 2.0d0*muShear + rho0*cs2
        denomShear = 4.0d0*muShear + 2.0d0*rho0*cs2
        denomDiv = 2.0d0*muBulk + rho0*cs2
        coeffTrace = (2.0d0*muShear-2.0d0*muBulk)/(2.0d0*denomDiag)

        !$acc parallel loop collapse(2) present(f,rho,u,v,Fx,Fy,Sxx,Sxy,Syy,Sdiv) &
        !$acc& private(i,j,alpha,feqLoc,neqxx,neqxy,neqyy)
        do j = 1, ny
            do i = 1, nx
                neqxx = 0.0d0
                neqxy = 0.0d0
                neqyy = 0.0d0
                do alpha = 0, 8
                    feqLoc = calc_feq(alpha, rho(i,j), u(i,j), v(i,j))
                    ! neq_ab = sum_i e_ia*e_ib*(f_i-feq_i)
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
!===================================================================================================


!===================================================================================================
! 子程序: compute_temperature_gradient
! 作用: 由温度分布函数的局部非平衡一阶矩 g-geq 直接得到 grad(T)，避免非局部有限差分。
! 论文对应: 式 (35)，本程序 dt=1；分子必须同时保留 T*F/rho0 和 u*Q 两个修正项。
!===================================================================================================
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
                ! 式 (35) 的分母。pressure/rho0 很小时仍保留，以保持与论文的压力耦合一致。
                denom = cs2*(2.0d0*taugL*(1.0d0-chi_kappa)+chi_kappa) + pressure(i,j)/rho0
                ! 防止极端非物理解使分母恰好为零；正常低 Mach 方腔不会触发这一保护。
                denom = sign(max(abs(denom), tiny(1.0d0)), denom)
                gradTx(i,j) = -(2.0d0*neqx + T(i,j)*Fx(i,j)/rho0 + u(i,j)*heatSourceQ)/denom
                gradTy(i,j) = -(2.0d0*neqy + T(i,j)*Fy(i,j)/rho0 + v(i,j)*heatSourceQ)/denom
            enddo
        enddo
    end subroutine compute_temperature_gradient
!===================================================================================================


!===================================================================================================
! 子程序: collision
! 作用: 执行流场 BGK 碰撞，并加入论文式 (26) 的 Guo 力项和 chi_s/chi_b 应变率修正项。
! 公式结构: f_post = f - (f-feq)/tau_fL + Phi；sourcePref=1-1/(2*tau_fL)。
! OpenACC: 空间循环并行，每个格点内部的 9 个方向由该线程顺序完成，没有跨格点写冲突。
!===================================================================================================
    subroutine collision()
        use commondata
        implicit none
        integer(kind=4) :: i, j, alpha
        real(kind=8) :: eu, eF, uF, feqLoc
        real(kind=8) :: Axx, Axy, Ayy, hermite2, phi
        real(kind=8), parameter :: sourcePref=1.0d0-0.5d0/taufL

        ! f_post 含一层 ghost。每步先清零可保证拉式迁移不会读到上一步遗留的 ghost 数据；
        ! 真实入射分布随后由 bounceback()/bouncebackT() 覆盖。这里优先保证测试代码清楚可靠。
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
                ! A_ab = chi_s*S_ab + (chi_b-chi_s)*Sdiv*delta_ab/D，二维 D=2。
                Axx = chi_s*Sxx(i,j)+0.5d0*(chi_b-chi_s)*Sdiv(i,j)
                Ayy = chi_s*Syy(i,j)+0.5d0*(chi_b-chi_s)*Sdiv(i,j)
                Axy = chi_s*Sxy(i,j)
                do alpha = 0, 8
                    eu = dble(qex(alpha))*u(i,j)+dble(qey(alpha))*v(i,j)
                    eF = dble(qex(alpha))*Fx(i,j)+dble(qey(alpha))*Fy(i,j)
                    feqLoc = calc_feq(alpha, rho(i,j), u(i,j), v(i,j))
                    ! (e_a*e_b/cs2-delta_ab)*A_ab；交叉项 xy 与 yx 合并后带系数 2。
                    hermite2 = (dble(qex(alpha)*qex(alpha))/cs2-1.0d0)*Axx + &
                        (dble(qey(alpha)*qey(alpha))/cs2-1.0d0)*Ayy + &
                        2.0d0*dble(qex(alpha)*qey(alpha))/cs2*Axy
                    phi = sourcePref*qomega(alpha)*(eF/cs2 + eu*eF/(cs2*cs2) - uF/cs2 + rho0*hermite2)
                    f_post(alpha,i,j) = f(alpha,i,j)-omegaF*(f(alpha,i,j)-feqLoc)+phi
                enddo
            enddo
        enddo
    end subroutine collision
!===================================================================================================


!===================================================================================================
! 子程序: collisionT
! 作用: 执行温度场 BGK 碰撞，并加入论文式 (28) 的完整标量源项 Psi_i。
! 源项三部分: p*grad(T)+T*F 耦合项、体热源 Q 项、chi_kappa*grad(T) 扩散率调节项。
!===================================================================================================
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
                    ! 大括号内逐项对应论文式 (28)，外层 qomega*sourcePref 在下一行乘入。
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
!===================================================================================================


!===================================================================================================
! 子程序: streaming
! 作用: 对流场执行拉式迁移：当前格点沿 alpha 方向的分布来自上游 (i-ex,j-ey) 的碰撞后分布。
! 说明: 拉式写法使每个线程只写自己的 f(:,i,j)，适合 GPU，且边界入射方向可在随后统一覆盖。
!===================================================================================================
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
!===================================================================================================


!===================================================================================================
! 子程序: streamingT
! 作用: 温度分布函数采用与流场相同的 D2Q9 拉式迁移。
!===================================================================================================
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
!===================================================================================================


!===================================================================================================
! 子程序: bounceback
! 作用: 对四个静止壁面施加标准 halfway bounce-back，对应论文式 (36) 在 uw=0 时的形式。
! 方向编号: 1/3 为东/西，2/4 为北/南，5/7 和 6/8 为两对反向对角速度。
!===================================================================================================
    subroutine bounceback()
        use commondata
        implicit none
        integer(kind=4) :: i, j

#ifdef VerticalWallsNoslip
        ! 左壁补入向东的 1/5/8，右壁补入向西的 3/6/7。
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
        ! 下壁补入向北的 2/5/6，上壁补入向南的 4/7/8。
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
!===================================================================================================


!===================================================================================================
! 子程序: bouncebackT
! 作用: 左/右恒温壁采用 anti-bounce-back（论文式 (38)），上/下绝热壁采用 bounce-back（式 (39)）。
! 关键点: 式 (38) 不只有 2*omega_i*Twall，还必须保留 2*lambda_i*Twall*p/(rho0*cs2) 压力项。
!===================================================================================================
    subroutine bouncebackT()
        use commondata
        implicit none
        integer(kind=4) :: i, j

#ifdef VerticalWallsConstT
        ! 左热壁 Thot；右冷壁 Tcold。每个赋值只更新从壁面进入流体的未知分布。
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
        ! 零法向热通量：把撞向壁面的分布原值反弹，不添加恒温修正。
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
        ! 四个角点的对角分布同时属于恒温侧壁和绝热水平壁。
        ! 上面的水平壁循环会覆盖侧壁先写入的对角分布，因此这里最后再写一次，使恒温 Dirichlet 条件优先。
        ! 角点只有 4 个标量赋值，使用 serial 比额外启动并行 kernel 更直接。
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
!===================================================================================================


!===================================================================================================
! 子程序: macro
! 作用: 由迁移和边界处理后的 f 重构 rho/u/v。
! 论文对应: sum(f_i)=deltaRho，rho0*u=sum(e_i*f_i)+F/2；半力修正来自梯形离散变换。
!===================================================================================================
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
!===================================================================================================


!===================================================================================================
! 子程序: macroT
! 作用: 由 g 的零阶矩重构温度；存在体热源时必须加 Q/2，本方腔 Q=0。
!===================================================================================================
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
!===================================================================================================


!===================================================================================================
! 子程序: check
! 作用: 在 CPU 端计算相邻两个检查时刻之间的速度/温度相对 L2 变化量，并更新参考场。
! 说明: 这是本测试程序的稳态停止判据，不是论文 LBM-CDE 方程的一部分；检查间隔为 checkInterval。
!===================================================================================================
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
!===================================================================================================


!===================================================================================================
! 子程序: calNuRe
! 作用: 在 CPU 端计算论文表 2 相关的壁面/中心线 Nu、中心线速度极值，以及补充的体平均 Nu/Re。
! 论文对应: Nu=-L/DeltaT*dT/dx，冷壁平均见式 (52)，中心竖线最大 Nu 见式 (52) 后的定义。
! 注意: 偶数网格的 x=0.5L/y=0.5H 位于两条流体节点线之间，因此本程序对两侧节点做线性插值。
!===================================================================================================
    subroutine calNuRe()
        use commondata
        implicit none
        integer(kind=4) :: i, j, imidLeft, imidRight, jmidBottom, jmidTop
        real(kind=8) :: sumHot, sumCold, sumMiddle, sumRe, sumGlobal
        real(kind=8) :: localNu, centerGradTx, centerU, centerV

        if(mod(nx,2).EQ.0) then
            imidLeft = nx/2
            imidRight = nx/2+1
        else
            imidLeft = nint(0.5d0*dble(nx))
            imidRight = imidLeft
        endif
        if(mod(ny,2).EQ.0) then
            jmidBottom = ny/2
            jmidTop = ny/2+1
        else
            jmidBottom = nint(0.5d0*dble(ny))
            jmidTop = jmidBottom
        endif

        sumHot = 0.0d0
        sumCold = 0.0d0
        sumMiddle = 0.0d0
        sumRe = 0.0d0
        sumGlobal = 0.0d0
        umax_mid = -huge(1.0d0)
        vmax_mid = -huge(1.0d0)
        Numax_mid = -huge(1.0d0)
        y_umax_mid = yp(1)
        x_vmax_mid = xp(1)
        y_Numax_mid = yp(1)

        do j = 1, ny
            ! gradTx 由论文式 (35) 在局部节点直接得到；壁面值用相邻的第一/最后一层流体节点代表。
            localNu = -lengthUnit*gradTx(1,j)/deltaT
            sumHot = sumHot + localNu
            localNu = -lengthUnit*gradTx(nx,j)/deltaT
            sumCold = sumCold + localNu

            ! x=0.5L 中心线：偶数网格平均中心左右两列，奇数网格两个下标相同，公式自动退化为单列值。
            centerGradTx = 0.5d0*(gradTx(imidLeft,j)+gradTx(imidRight,j))
            centerU = 0.5d0*(u(imidLeft,j)+u(imidRight,j))
            localNu = -lengthUnit*centerGradTx/deltaT
            sumMiddle = sumMiddle + localNu
            if(localNu.GT.Numax_mid) then
                Numax_mid = localNu
                y_Numax_mid = yp(j)
            endif
            ! 论文表 2 的 umax 是“最大正值”，不是绝对值最大的带符号速度。
            if(centerU.GT.umax_mid) then
                umax_mid = centerU
                y_umax_mid = yp(j)
            endif
        enddo

        do i = 1, nx
            ! y=0.5H 中心线同理在上下两行之间插值。
            centerV = 0.5d0*(v(i,jmidBottom)+v(i,jmidTop))
            ! 论文表 2 的 vmax 同样取最大正值；其位置应落在靠近左侧热壁的一侧。
            if(centerV.GT.vmax_mid) then
                vmax_mid = centerV
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
        ! 补充诊断：侧热方腔主热流方向为 x，因此体平均对流项使用 u*(T-Tref)。
        Nu_global = 1.0d0 + sumGlobal/dble(nx*ny)*lengthUnit/(diffusivity*deltaT)
        ReVolAvg = dsqrt(sumRe/dble(nx*ny))*lengthUnit/viscosity
        ! 论文表 2 的速度采用自由落体速度 u0 归一化。
        umax_mid = umax_mid/velocityUnit
        vmax_mid = vmax_mid/velocityUnit
    end subroutine calNuRe
!===================================================================================================


!===================================================================================================
! 子程序: append_history
! 作用: 每次检查后追加收敛误差和主要 Nu/Re；时间列采用自由落体时间 t/t0=itc/timeUnit。
!===================================================================================================
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
!===================================================================================================


!===================================================================================================
! 子程序: output_Tecplot
! 作用: 输出当前宏观场和由论文式 (35) 得到的温度梯度，便于检查边界层和局部 Nu。
! 说明: 调用前必须先执行 update_host_monitor_2d_openacc()，文件 I/O 始终在主机端进行。
!===================================================================================================
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
!===================================================================================================


!===================================================================================================
! 子程序: final_log
! 作用: 把最终收敛量、论文表 2 相关诊断和 CPU/墙钟时间写到终端及设置文件。
!===================================================================================================
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
!===================================================================================================


!===================================================================================================
! 子程序: finalize_arrays
! 作用: 在设备数据已经删除后，释放全部主机 allocatable 数组。
!===================================================================================================
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
!===================================================================================================
