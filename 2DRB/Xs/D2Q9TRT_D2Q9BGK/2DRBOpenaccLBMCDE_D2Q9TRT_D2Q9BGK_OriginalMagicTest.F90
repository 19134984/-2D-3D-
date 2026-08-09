!===================================================================================================
!!! 二维高 Rayleigh 数自然对流：LBM-CDE（OpenACC 单 GPU）
!!!
!!! 默认工况：Rayleigh-Benard，Ra=1.0e10，Pr=0.7，论文等效 2049x2049 网格点，非稳态统计。
!!! 可选工况：侧壁加热方腔；可选推进：稳态收敛或非稳态定时统计。
!!! 算法来源：pdf/LBM-CDE.pdf；目标 RB 工况来源：
!!! pdf/Statistics of kinetic and thermal energy dissipation rates in two-dimensional turbulent
!!! Rayleigh-Benard convection.pdf。
!!! 离散模型：流场 D2Q9-TRT；温度场严格采用原 LBM-CDE 的 D2Q9-BGK。
!!! 原始 magic 对照试验：流场奇模态按基准偶松弛尺度 (tau_fL-0.5) 满足 magic 3/16。
!!! 主要论文对应：
!!!   1) 流场平衡态与源项：式 (4)、式 (26)，有效运动黏度见式 (27)；
!!!   2) 温度平衡态与源项：式 (16)、式 (28)，有效扩散率见式 (29)；
!!!   3) 局部应变率和温度梯度：式 (31)-(35)；
!!!   4) 无滑移、恒温和绝热边界：式 (36)、式 (38)、式 (39)；
!!!   5) 侧热方腔算例及 Nu 定义：式 (50)-(52)；RB 的主热流方向改为 y。
!!!
!!! 重要网格口径：目标 DNS 的 2049^2 是包含两端物理边界的节点数，故 Delta/H=1/(2049-1)。
!!! 本程序采用 halfway 壁面与 2048^2 个流体中心节点，同样有 Delta/H=1/2048；两者空间步长等效。
!!! D2Q9 的 halfway BB/ABB 与 BGK 内点并非严格二阶全局一致，正式计算仍需独立检查边界层
!!! 分辨率、耗散尺度和网格收敛，不能仅凭点数相同宣称复现 DNS。
!===================================================================================================


!===================================================================================================
! 算例开关
! steadyFlow：达到 errorU/errorT 阈值后停止；unsteadyFlow：按自由落体时间采样统计。
! 未从编译命令指定时，默认选择高 Ra 所需的 unsteadyFlow。
! 这里的“湍流”是按网格直接解析的非稳态 LBM，不含额外 LES/RANS/SGS 黏度模型。
#if !defined(steadyFlow) && !defined(unsteadyFlow)
#define unsteadyFlow
#endif

#if defined(steadyFlow) && defined(unsteadyFlow)
#error "Choose only one flow mode: steadyFlow or unsteadyFlow"
#endif
! 本程序固定为封闭方腔：四壁速度均采用 halfway bounce-back 无滑移边界。
#define HorizontalWallsNoslip
#define VerticalWallsNoslip

! 未从编译命令指定时，默认选择 Rayleigh-Benard：底热、顶冷、左右绝热。
! 侧热方腔可用 -DSideHeatedCell 选择：左热、右冷、上下绝热。
#if !defined(RayleighBenardCell) && !defined(SideHeatedCell)
#define RayleighBenardCell
#endif
#if defined(RayleighBenardCell) && defined(SideHeatedCell)
#error "Choose only one thermal case: RayleighBenardCell or SideHeatedCell"
#endif
#ifdef RayleighBenardCell
#define HorizontalWallsConstT
#define VerticalWallsAdiabatic
#endif
#ifdef SideHeatedCell
#define HorizontalWallsAdiabatic
#define VerticalWallsConstT
#endif

! 下面的检查用于防止以后改宏时得到一个“能编译但边界不完整”的程序。
#if !defined(HorizontalWallsNoslip) || !defined(VerticalWallsNoslip)
#error "This test requires no-slip velocity boundaries on all four walls"
#endif
#if defined(RayleighBenardCell) && (!defined(HorizontalWallsConstT) || !defined(VerticalWallsAdiabatic))
#error "RayleighBenardCell requires horizontal constant-T and vertical adiabatic walls"
#endif
#if defined(SideHeatedCell) && (!defined(HorizontalWallsAdiabatic) || !defined(VerticalWallsConstT))
#error "SideHeatedCell requires horizontal adiabatic and vertical constant-T walls"
#endif

! 目标 DNS 使用 2049^2 个含边界节点；halfway LBM 用 2048^2 个流体中心节点得到相同 Delta/H。
! 可用 -DNX_OVERRIDE=... -DNY_OVERRIDE=... 覆盖流体节点数以做网格筛查。
#ifndef NX_OVERRIDE
#define NX_OVERRIDE 2048
#endif
#ifndef NY_OVERRIDE
#define NY_OVERRIDE 2048
#endif
#if NX_OVERRIDE != NY_OVERRIDE
#error "The current cavity cases require NX_OVERRIDE == NY_OVERRIDE"
#endif
#ifndef ITC_MAX_OVERRIDE
#define ITC_MAX_OVERRIDE -1
#endif
#ifndef RAYLEIGH_OVERRIDE
#define RAYLEIGH_OVERRIDE 1.0d10
#endif
#ifndef PRANDTL_OVERRIDE
#define PRANDTL_OVERRIDE 0.7d0
#endif
#ifndef MACH_OVERRIDE
#define MACH_OVERRIDE 0.1d0
#endif
#ifndef BASE_TAUF_OVERRIDE
#define BASE_TAUF_OVERRIDE 0.64d0
#endif
#ifndef CHI_B_OVERRIDE
#define CHI_B_OVERRIDE 0.0d0
#endif
#ifndef BASE_TAUG_OVERRIDE
#define BASE_TAUG_OVERRIDE 0.64d0
#endif
#if defined(ThermalD2Q5) || defined(D2Q5_QUARTIC_RATES)
#error "Temperature is fixed to the original D2Q9-BGK LBM-CDE; D2Q5/Xu/Luo branches are disabled"
#endif
#define ThermalD2Q9
#ifndef FLOW_ODD_BASE_MAGIC_TEST
#define FLOW_ODD_BASE_MAGIC_TEST
#endif
#if defined(FLOW_ODD_UNIT) || defined(FLOW_BGK) || defined(FLOW_ODD_RATE_OVERRIDE)
#error "Direct odd-rate/BGK branches are disabled; this test fixes the original/base flow magic policy"
#endif
#ifndef FLOW_MAGIC_PARAMETER_OVERRIDE
#define FLOW_MAGIC_PARAMETER_OVERRIDE (3.0d0/16.0d0)
#endif
#ifndef RB_PERTURBATION_OVERRIDE
#define RB_PERTURBATION_OVERRIDE 1.0d-3
#endif
#ifndef CHECK_INTERVAL_OVERRIDE
#define CHECK_INTERVAL_OVERRIDE 1000
#endif
#ifndef UNSTEADY_TRANSIENT_TF_OVERRIDE
#define UNSTEADY_TRANSIENT_TF_OVERRIDE 500.0d0
#endif
#ifndef UNSTEADY_AVERAGE_TF_OVERRIDE
#define UNSTEADY_AVERAGE_TF_OVERRIDE 500.0d0
#endif
#ifndef UNSTEADY_SAMPLE_TF_OVERRIDE
#define UNSTEADY_SAMPLE_TF_OVERRIDE 5.0d0
#endif
#ifndef UNSTEADY_MONITOR_TF_OVERRIDE
#define UNSTEADY_MONITOR_TF_OVERRIDE 5.0d0
#endif
#ifndef UNSTEADY_SNAPSHOT_TF_OVERRIDE
#define UNSTEADY_SNAPSHOT_TF_OVERRIDE 100.0d0
#endif
#ifndef UNSTEADY_TECPLOT_TF_OVERRIDE
#define UNSTEADY_TECPLOT_TF_OVERRIDE 100.0d0
#endif
#ifndef UNSTEADY_OUTPUT_SNAPSHOT_OVERRIDE
#define UNSTEADY_OUTPUT_SNAPSHOT_OVERRIDE 1
#endif
#ifndef UNSTEADY_OUTPUT_TECPLOT_OVERRIDE
#define UNSTEADY_OUTPUT_TECPLOT_OVERRIDE 1
#endif
!===================================================================================================


!===================================================================================================
! 全局数据模块
    module commondata
        implicit none

        !-------------------------------------------------------------------------------------------
        ! 网格与格子常量
        integer(kind=4), parameter :: nx=NX_OVERRIDE, ny=NY_OVERRIDE ! 流体节点数；壁面位于首末节点外半格
        real(kind=8), parameter :: rho0=1.0d0                        ! 弱可压缩模型的参考密度
        real(kind=8), parameter :: cs2=1.0d0/3.0d0                   ! D2Q9 格子声速平方
        integer(kind=4), parameter :: thermalQ=8
        real(kind=8), parameter :: cT2=cs2

#ifdef SideHeatedCell
        real(kind=8), parameter :: lengthUnit=dble(nx)               ! 侧热方腔的特征长度 L
#else
        real(kind=8), parameter :: lengthUnit=dble(ny)
#endif

        !-------------------------------------------------------------------------------------------
        ! 无量纲工况参数
        real(kind=8), parameter :: Rayleigh=RAYLEIGH_OVERRIDE        ! Ra，论文式 (51)
        real(kind=8), parameter :: Prandtl=PRANDTL_OVERRIDE          ! 用户指定默认 Pr=0.7
        real(kind=8), parameter :: Mach=MACH_OVERRIDE                ! 低 Mach 渐近假设；默认 0.1
        real(kind=8), parameter :: Thot=0.5d0, Tcold=-0.5d0          ! 热壁/冷壁温度
        real(kind=8), parameter :: Tref=0.5d0*(Thot+Tcold)           ! Boussinesq 参考温度 T0
        real(kind=8), parameter :: deltaT=Thot-Tcold                 ! 温差 DeltaT

        !-------------------------------------------------------------------------------------------
        ! 格子单位下的目标输运系数
        ! 物理映射不把 nu/kappa/gBeta 当作独立手调输入：由 Ra、Pr、Ma、H_LB 自动确定。
        real(kind=8), parameter :: viscosity=Mach*lengthUnit*dsqrt(Prandtl/(3.0d0*Rayleigh))
        real(kind=8), parameter :: diffusivity=viscosity/Prandtl

        !-------------------------------------------------------------------------------------------
        ! LBM-CDE 数值自由度与自动映射
        ! 论文式 (27)/(29)，dt=1：
        !   nu    = (tau_fL-0.5)*(1-chi_s)*cs2
        !   kappa = (tau_gL-0.5)*(1-chi_kappa)*cT2
        ! Ra/Pr/Ma/H_LB 只能确定 nu/kappa，不能同时唯一确定 tau 与 chi。默认采用论文讨论的
        ! base tau=0.64 稳定化策略，再由目标 nu/kappa 自动反算 chi；这不是手工指定物理输运系数。
        ! 若需复现论文的固定-chi 扫描，可用 -DCHI_S_OVERRIDE=.../-DCHI_KAPPA_OVERRIDE=...，
        ! 此时相应 tau 由物理输运系数自动反算。
#ifdef CHI_S_OVERRIDE
        real(kind=8), parameter :: chi_s=CHI_S_OVERRIDE
        real(kind=8), parameter :: taufL=0.5d0+viscosity/(cs2*(1.0d0-chi_s))
#else
        real(kind=8), parameter :: taufL=BASE_TAUF_OVERRIDE
        real(kind=8), parameter :: chi_s=1.0d0-viscosity/(cs2*(taufL-0.5d0))
#endif
#ifdef CHI_KAPPA_OVERRIDE
        real(kind=8), parameter :: chi_kappa=CHI_KAPPA_OVERRIDE
        real(kind=8), parameter :: taugL=0.5d0+diffusivity/(cT2*(1.0d0-chi_kappa))
#else
        real(kind=8), parameter :: taugL=BASE_TAUG_OVERRIDE
        real(kind=8), parameter :: chi_kappa=1.0d0-diffusivity/(cT2*(taugL-0.5d0))
#endif
        real(kind=8), parameter :: chi_b=CHI_B_OVERRIDE
        real(kind=8), parameter :: rbPerturbation=RB_PERTURBATION_OVERRIDE
        real(kind=8), parameter :: heatSourceQ=0.0d0                  ! 方腔无体热源；保留 Q 以实现完整式 (28)/(35)
        real(kind=8), parameter :: bulkViscosity=(taufL-0.5d0)*(1.0d0-chi_b)*cs2 ! D=2 时的 nu_B
        real(kind=8), parameter :: omegaF=1.0d0/taufL
        real(kind=8), parameter :: omegaG=1.0d0/taugL
        real(kind=8), parameter :: sigmaFlowEffective=viscosity/cs2
        real(kind=8), parameter :: sigmaFlowBase=taufL-0.5d0
        real(kind=8), parameter :: flowMagicParameter=FLOW_MAGIC_PARAMETER_OVERRIDE
        real(kind=8), parameter :: omegaFOdd=1.0d0/(0.5d0+flowMagicParameter/sigmaFlowBase)
        real(kind=8), parameter :: taufOdd=1.0d0/omegaFOdd
        real(kind=8), parameter :: sigmaFlowOdd=1.0d0/omegaFOdd-0.5d0
        real(kind=8), parameter :: magicProductBase=sigmaFlowBase*sigmaFlowOdd
        real(kind=8), parameter :: magicProductEffective=sigmaFlowEffective*sigmaFlowOdd
        real(kind=8), parameter :: sigmaThermalEffective=diffusivity/cT2

        !-------------------------------------------------------------------------------------------
        ! Boussinesq 浮力与自由落体标度，见论文式 (50)/(51)。压力 p 始终指
        ! 弱可压缩压力扰动 p'=cs2*(rho-rho0)，而不是包含常数热力学基准压的绝对压力。
        ! Fy>0 表示热流体向 +y 上升；gBeta 在这里代表 g*beta。
        real(kind=8), parameter :: gBeta=Rayleigh*viscosity*diffusivity/(deltaT*lengthUnit**3)
        real(kind=8), parameter :: timeUnit=dsqrt(lengthUnit/(gBeta*deltaT))
        real(kind=8), parameter :: velocityUnit=dsqrt(gBeta*deltaT*lengthUnit)

        !-------------------------------------------------------------------------------------------
        ! 推进与收敛设置。误差定义见 check()，它不是论文模型的一部分。
#ifdef steadyFlow
        real(kind=8), parameter :: epsU=1.0d-7, epsT=1.0d-7
        integer(kind=4), parameter :: itc_max=merge(ITC_MAX_OVERRIDE,20000000,ITC_MAX_OVERRIDE.GT.0)
        integer(kind=4), parameter :: checkInterval=max(1,CHECK_INTERVAL_OVERRIDE)
#endif
#ifdef unsteadyFlow
        real(kind=8), parameter :: transientTimeTf=UNSTEADY_TRANSIENT_TF_OVERRIDE
        real(kind=8), parameter :: averagingTimeTf=UNSTEADY_AVERAGE_TF_OVERRIDE
        real(kind=8), parameter :: sampleTimeTf=UNSTEADY_SAMPLE_TF_OVERRIDE
        real(kind=8), parameter :: monitorTimeTf=UNSTEADY_MONITOR_TF_OVERRIDE
        real(kind=8), parameter :: snapshotTimeTf=UNSTEADY_SNAPSHOT_TF_OVERRIDE
        real(kind=8), parameter :: tecplotTimeTf=UNSTEADY_TECPLOT_TF_OVERRIDE
        integer(kind=4), parameter :: outputSnapshotFile=UNSTEADY_OUTPUT_SNAPSHOT_OVERRIDE
        integer(kind=4), parameter :: outputTecplotFile=UNSTEADY_OUTPUT_TECPLOT_OVERRIDE
        real(kind=8), parameter :: sampleTimeForParameters=merge(sampleTimeTf,1.0d0,sampleTimeTf.GT.0.0d0)
        real(kind=8), parameter :: monitorTimeForParameters=merge(monitorTimeTf,1.0d0,monitorTimeTf.GT.0.0d0)
        integer(kind=4), parameter :: averageStartItc=max(0,nint(transientTimeTf*timeUnit))
        integer(kind=4), parameter :: sampleIntervalItc=max(1,nint(sampleTimeForParameters*timeUnit))
        integer(kind=4), parameter :: targetStatisticSamples=max(1,nint(max(averagingTimeTf,0.0d0)/sampleTimeForParameters))
        integer(kind=4), parameter :: targetSnapshotFiles=max(0, &
            int(max(averagingTimeTf,0.0d0)/max(snapshotTimeTf,1.0d-12)+1.0d-12))
        integer(kind=4), parameter :: targetTecplotFiles=max(0, &
            int(max(averagingTimeTf,0.0d0)/max(tecplotTimeTf,1.0d-12)+1.0d-12))
        integer(kind=4), parameter :: defaultUnsteadyMaxItc=averageStartItc+targetStatisticSamples*sampleIntervalItc
        integer(kind=4), parameter :: itc_max=merge(ITC_MAX_OVERRIDE,defaultUnsteadyMaxItc,ITC_MAX_OVERRIDE.GT.0)
        integer(kind=4), parameter :: checkInterval=max(1,nint(monitorTimeForParameters*timeUnit))
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
        real(kind=8) :: Nu_hot_halfway, Nu_cold_halfway, Nu_middleFD
        real(kind=8) :: ReVolAvg, speedSquaredVol, umax_mid, vmax_mid, Numax_mid
        real(kind=8) :: epsKineticVol, epsThermalVol, NuKineticDiss, NuThermalDiss
        real(kind=8) :: kineticDissExactRatio, thermalDissExactRatio
        real(kind=8) :: maxMachLocal, minTemperature, maxTemperature, maxDensityDeviation
        real(kind=8) :: pos_umax_mid, pos_vmax_mid, pos_Numax_mid
        real(kind=8) :: sumNuGlobal, sumNuGlobal2, sumNuHot, sumNuHot2
        real(kind=8) :: sumNuCold, sumNuCold2, sumNuMiddle, sumNuMiddle2
        real(kind=8) :: sumSpeedSquaredVol, sumSpeedSquaredVol2
        real(kind=8) :: sumNuKineticDiss, sumNuKineticDiss2
        real(kind=8) :: sumNuThermalDiss, sumNuThermalDiss2
        real(kind=8) :: sumNuGlobalFirstHalf, sumNuGlobalSecondHalf
        real(kind=8) :: sumSpeedSquaredFirstHalf, sumSpeedSquaredSecondHalf
        real(kind=8) :: maxStatisticCFL
        real(kind=8), allocatable :: sumTemperatureProfile(:), sumTemperatureSquaredProfile(:)
        integer(kind=4) :: statisticSamples, firstHalfSamples, secondHalfSamples
        integer(kind=4) :: firstStatisticItc, lastStatisticItc
        integer(kind=4) :: itc, lastCheckItc
#ifdef unsteadyFlow
        integer(kind=4) :: snapshotFileNum, tecplotFileNum, lastSnapshotItc, lastTecplotItc
#endif

        character(len=100) :: settingsFile="SimulationSettings2DOpenaccLBMCDE.txt"
        character(len=100) :: pltFile="buoyancyCavity2DOpenaccLBMCDE.plt"
        character(len=100) :: historyFile="NuRe_2DOpenaccLBMCDE.dat"
        character(len=100) :: statisticsFile="NuReStatistics_2DOpenaccLBMCDE.dat"
#ifdef unsteadyFlow
        character(len=100) :: snapshotFilePrefix="buoyancyCavity2DOpenaccLBMCDESnapshot"
        character(len=100) :: snapshotIndexFile="buoyancyCavity2DOpenaccLBMCDESnapshot-index.dat"
        character(len=100) :: tecplotFilePrefix="buoyancyCavity2DOpenaccLBMCDETecplot"
#endif

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

        ! 温度格子严格使用 LBM-CDE 原文 D2Q9 权重，与流场 D2Q9 权重相同。
        pure real(kind=8) function tomega(alpha)
            !$acc routine seq
            implicit none
            integer(kind=4), intent(in) :: alpha

            tomega = qomega(alpha)
        end function tomega

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
        ! 压力项为 lambda_i*T*p/(rho0*cs2)，保持零/一阶压力矩为零、二阶压力矩为 T*p/rho0。
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
#ifdef unsteadyFlow
        integer(kind=4) :: nextSnapshotItc, nextTecplotItc
        logical :: monitorDue, statisticDue, snapshotDue, tecplotDue
#endif

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

#ifdef steadyFlow
            if(mod(itc, checkInterval).EQ.0) then
                ! 先在设备端刷新诊断量，再只同步 CPU 后处理所需的宏观数组。
                call compute_force()
                call compute_strain_rate()
                call compute_temperature_gradient()
                call update_host_monitor_2d_openacc()
                call check()
                call calNuRe()
                call append_history()
            endif
#endif
#ifdef unsteadyFlow
            ! 标量统计、二进制瞬时场和 Tecplot 使用彼此独立的自由落体时间间隔。
            ! 瞬时场只在统计窗口内保存；目标时刻由绝对 t_ff 逐次换算，避免整数步间隔累积漂移。
            monitorDue = mod(itc,checkInterval).EQ.0
            statisticDue = (itc.GT.averageStartItc).AND. &
                (mod(itc-averageStartItc,sampleIntervalItc).EQ.0)
            snapshotDue = .false.
            if((outputSnapshotFile.EQ.1).AND.(snapshotFileNum.LT.targetSnapshotFiles)) then
                nextSnapshotItc = max(1,nint((transientTimeTf+ &
                    dble(snapshotFileNum+1)*snapshotTimeTf)*timeUnit))
                snapshotDue = itc.GE.nextSnapshotItc
            endif
            tecplotDue = .false.
            if((outputTecplotFile.EQ.1).AND.(tecplotFileNum.LT.targetTecplotFiles)) then
                nextTecplotItc = max(1,nint((transientTimeTf+ &
                    dble(tecplotFileNum+1)*tecplotTimeTf)*timeUnit))
                tecplotDue = itc.GE.nextTecplotItc
            endif

            ! Tecplot 与标量诊断都需要温度梯度；重合时只做一次设备到主机同步。
            if(monitorDue.OR.statisticDue.OR.tecplotDue) then
                call compute_force()
                call compute_strain_rate()
                call compute_temperature_gradient()
                call update_host_monitor_2d_openacc()
                if(monitorDue.OR.statisticDue) then
                    call check()
                    call calNuRe()
                    if(statisticDue) call accumulate_statistics()
                    call append_history()
                endif
                if(snapshotDue) call output_SnapshotFile()
                if(tecplotDue) call output_periodic_Tecplot()
            elseif(snapshotDue) then
                call update_host_snapshot_2d_openacc()
                call output_SnapshotFile()
            endif
#endif
        enddo

        call compute_force()
        call compute_strain_rate()
        call compute_temperature_gradient()
        call update_host_monitor_2d_openacc()
        if(lastCheckItc.NE.itc) call check()
        call calNuRe()
#ifdef unsteadyFlow
        ! 无论最后一个格子步是否恰好命中周期，都保留最终瞬时场；避免只留下较早的周期快照。
        if((outputSnapshotFile.EQ.1).AND.(lastSnapshotItc.NE.itc)) call output_SnapshotFile()
        call write_statistics()
#endif
        call output_Tecplot(pltFile)

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
        real(kind=8) :: piValue
#ifdef RayleighBenardCell
        real(kind=8) :: yNormalized
#endif

        itc = 0
        lastCheckItc = -1
        errorU = 100.0d0
        errorT = 100.0d0
        Nu_global = 0.0d0
        Nu_hot = 0.0d0
        Nu_cold = 0.0d0
        Nu_middle = 0.0d0
        Nu_hot_halfway = 0.0d0
        Nu_cold_halfway = 0.0d0
        Nu_middleFD = 0.0d0
        ReVolAvg = 0.0d0
        speedSquaredVol = 0.0d0
        epsKineticVol = 0.0d0
        epsThermalVol = 0.0d0
        NuKineticDiss = 0.0d0
        NuThermalDiss = 0.0d0
        kineticDissExactRatio = 0.0d0
        thermalDissExactRatio = 0.0d0
        maxMachLocal = 0.0d0
        minTemperature = 0.0d0
        maxTemperature = 0.0d0
        maxDensityDeviation = 0.0d0
        umax_mid = 0.0d0
        vmax_mid = 0.0d0
        Numax_mid = 0.0d0
        pos_umax_mid = 0.0d0
        pos_vmax_mid = 0.0d0
        pos_Numax_mid = 0.0d0
        sumNuGlobal = 0.0d0
        sumNuGlobal2 = 0.0d0
        sumNuHot = 0.0d0
        sumNuHot2 = 0.0d0
        sumNuCold = 0.0d0
        sumNuCold2 = 0.0d0
        sumNuMiddle = 0.0d0
        sumNuMiddle2 = 0.0d0
        sumSpeedSquaredVol = 0.0d0
        sumSpeedSquaredVol2 = 0.0d0
        sumNuKineticDiss = 0.0d0
        sumNuKineticDiss2 = 0.0d0
        sumNuThermalDiss = 0.0d0
        sumNuThermalDiss2 = 0.0d0
        sumNuGlobalFirstHalf = 0.0d0
        sumNuGlobalSecondHalf = 0.0d0
        sumSpeedSquaredFirstHalf = 0.0d0
        sumSpeedSquaredSecondHalf = 0.0d0
        maxStatisticCFL = 0.0d0
        statisticSamples = 0
        firstHalfSamples = 0
        secondHalfSamples = 0
        firstStatisticItc = -1
        lastStatisticItc = -1
#ifdef unsteadyFlow
        snapshotFileNum = 0
        tecplotFileNum = 0
        lastSnapshotItc = -1
        lastTecplotItc = -1
#endif

        if((Rayleigh.LE.0.0d0).OR.(Prandtl.LE.0.0d0).OR.(Mach.LE.0.0d0)) then
            error stop "Ra, Pr and Mach must all be positive"
        endif
        if((chi_s.GE.1.0d0).OR.(chi_b.GE.1.0d0).OR.(chi_kappa.GE.1.0d0)) then
            error stop "chi_s, chi_b and chi_kappa must be smaller than one"
        endif
        if((taufL.LE.0.5d0).OR.(taugL.LE.0.5d0)) then
            error stop "The mapped shear/thermal relaxation times must be larger than 0.5"
        endif
        if((omegaFOdd.LE.0.0d0).OR.(omegaFOdd.GE.2.0d0)) then
            error stop "The TRT odd relaxation rate must lie strictly between zero and two"
        endif
        if(sigmaFlowBase.LE.0.0d0) then
            error stop "The original/base flow Henon shift must be positive"
        endif
        if(sigmaFlowEffective.LE.0.0d0) then
            error stop "The chi_s-corrected effective flow Hénon shift must be positive"
        endif
        if(flowMagicParameter.LE.0.0d0) then
            error stop "The TRT magic parameter must be positive"
        endif
#ifdef unsteadyFlow
        if((transientTimeTf.LT.0.0d0).OR.(averagingTimeTf.LE.0.0d0).OR. &
            (sampleTimeTf.LE.0.0d0).OR.(monitorTimeTf.LE.0.0d0)) then
            error stop "Unsteady transient/average/sample/monitor times are invalid"
        endif
        if((outputSnapshotFile.EQ.1).AND.(snapshotTimeTf.LE.0.0d0)) then
            error stop "The unsteady snapshot interval must be positive when snapshot output is enabled"
        endif
        if((outputTecplotFile.EQ.1).AND.(tecplotTimeTf.LE.0.0d0)) then
            error stop "The unsteady Tecplot interval must be positive when Tecplot output is enabled"
        endif
        if((outputSnapshotFile.NE.0).AND.(outputSnapshotFile.NE.1)) then
            error stop "UNSTEADY_OUTPUT_SNAPSHOT_OVERRIDE must be zero or one"
        endif
        if((outputTecplotFile.NE.0).AND.(outputTecplotFile.NE.1)) then
            error stop "UNSTEADY_OUTPUT_TECPLOT_OVERRIDE must be zero or one"
        endif
#endif

        allocate(u(nx,ny), v(nx,ny), T(nx,ny), rho(nx,ny), pressure(nx,ny))
        allocate(up(nx,ny), vp(nx,ny), Tp(nx,ny))
        allocate(Fx(nx,ny), Fy(nx,ny))
        allocate(gradTx(nx,ny), gradTy(nx,ny))
        allocate(Sxx(nx,ny), Sxy(nx,ny), Syy(nx,ny), Sdiv(nx,ny))
        allocate(f(0:8,nx,ny), f_post(0:8,0:nx+1,0:ny+1))
        ! 温度分布严格按原 LBM-CDE 的 D2Q9 方向 0:8 分配。
        allocate(g(0:thermalQ,nx,ny), g_post(0:thermalQ,0:nx+1,0:ny+1))
        allocate(sumTemperatureProfile(ny), sumTemperatureSquaredProfile(ny))
        sumTemperatureProfile = 0.0d0
        sumTemperatureSquaredProfile = 0.0d0

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

        ! 初始温度取相应方向的纯导热线性分布。RB 额外加入确定性、零均值的小扰动，
        ! 以打破精确导热态的离散对称性；扰动在四条物理壁面处解析为零。
        piValue = dacos(-1.0d0)
        do j = 1, ny
            do i = 1, nx
#ifdef RayleighBenardCell
                T(i,j) = Thot + (yp(j)-yp(0))/(yp(ny+1)-yp(0))*(Tcold-Thot) + &
                    rbPerturbation*deltaT*dsin(2.0d0*piValue*xp(i))*dsin(piValue*yp(j))
                ! 纯导热基态满足 dp/dy=rho0*gBeta*(T_cond-Tref)。取体平均压力为零，
                ! 可避免从均匀密度启动时产生与碰撞参数混杂的全域静水调整瞬态。
                yNormalized = (yp(j)-yp(0))/(yp(ny+1)-yp(0))
                pressure(i,j) = rho0*gBeta*lengthUnit* &
                    (0.5d0*yNormalized-0.5d0*yNormalized*yNormalized-1.0d0/12.0d0)
                rho(i,j) = rho0+pressure(i,j)/cs2
#endif
#ifdef SideHeatedCell
                T(i,j) = Thot + (xp(i)-xp(0))/(xp(nx+1)-xp(0))*(Tcold-Thot)
#endif
            enddo
        enddo

        f = 0.0d0
        g = 0.0d0
        ! 由初始宏观量建立平衡分布，避免第一步人为引入额外非平衡扰动。
        do j = 1, ny
            do i = 1, nx
                do alpha = 0, 8
                    f(alpha,i,j) = calc_feq(alpha, rho(i,j), u(i,j), v(i,j))
                enddo
                do alpha = 0, thermalQ
                    g(alpha,i,j) = calc_geq(alpha, T(i,j), u(i,j), v(i,j), rho(i,j))
                enddo
            enddo
        enddo

        up = u
        vp = v
        Tp = T

        open(unit=00, file=trim(settingsFile), status='replace', action='write')
        write(00,*) "2D OpenACC high-Ra LBM-CDE"
        write(00,*) "Halfway-LBM fluid-node mesh =", nx, ny
        write(00,*) "Equivalent boundary-node resolution =", nx+1, ny+1
        write(00,*) "Device model = OpenACC single GPU"
        write(00,*) "Visible OpenACC default devices =", numAccDevices
#ifdef RayleighBenardCell
        write(00,*) "Case = RayleighBenardCell (bottom hot, top cold, sidewalls adiabatic)"
#endif
#ifdef SideHeatedCell
        write(00,*) "Case = SideHeatedCell"
#endif
#ifdef steadyFlow
        write(00,*) "Flow mode = steadyFlow"
#endif
#ifdef unsteadyFlow
        write(00,*) "Flow mode = unsteadyFlow"
        write(00,*) "Turbulence treatment = time-resolved LBM without an SGS/RANS closure"
#endif
        write(00,*) "Flow lattice = D2Q9 TRT with LBM-CDE even stress source"
        write(00,*) "Flow odd policy = original magic parameter on base (tau_fL-0.5) even scale"
        write(00,*) "Temperature lattice = D2Q9 BGK LBM-CDE"
        write(00,*) "Temperature cT2 =", real(cT2,kind=8)
        write(00,*) "Rayleigh =", real(Rayleigh,kind=8), "Prandtl =", real(Prandtl,kind=8)
        write(00,*) "Mach =", real(Mach,kind=8), "Length unit =", real(lengthUnit,kind=8)
        write(00,*) "chi_s =", real(chi_s,kind=8), "chi_kappa =", real(chi_kappa,kind=8), "chi_b =", real(chi_b,kind=8)
#ifdef CHI_S_OVERRIDE
        write(00,*) "Flow mapping policy = explicit chi_s; tau_fL derived automatically"
#else
        write(00,*) "Flow mapping policy = base tau_fL target; chi_s derived automatically"
#endif
#ifdef CHI_KAPPA_OVERRIDE
        write(00,*) "Thermal mapping policy = explicit chi_kappa; tau_gL derived automatically"
#else
        write(00,*) "Thermal mapping policy = base tau_gL target; chi_kappa derived automatically"
#endif
        write(00,*) "viscosity =", real(viscosity,kind=8), "diffusivity =", real(diffusivity,kind=8)
        write(00,*) "taufL =", real(taufL,kind=8), "taugL =", real(taugL,kind=8)
        write(00,*) "tau margins from 0.5 =", real(taufL-0.5d0,kind=8), real(taugL-0.5d0,kind=8)
        write(00,*) "tau flow even/odd nominal =", real(taufL,kind=8), real(taufOdd,kind=8)
        write(00,*) "omega flow even/odd =", real(omegaF,kind=8), real(omegaFOdd,kind=8)
        write(00,*) "sigma flow base/effective/odd =", real(sigmaFlowBase,kind=8), &
            real(sigmaFlowEffective,kind=8), real(sigmaFlowOdd,kind=8)
        write(00,*) "TRT original/base magic target/actual =", real(flowMagicParameter,kind=8), &
            real(magicProductBase,kind=8)
        write(00,*) "TRT chi_s-corrected effective product (diagnostic only) =", &
            real(magicProductEffective,kind=8)
        write(00,*) "thermal effective sigma =", real(sigmaThermalEffective,kind=8)
        write(00,*) "gBeta =", real(gBeta,kind=8)
        write(00,*) "timeUnit =", real(timeUnit,kind=8), "velocityUnit =", real(velocityUnit,kind=8)
        write(00,*) "Automatic mapping check: Ra =", &
            real(gBeta*deltaT*lengthUnit**3/(viscosity*diffusivity),kind=8)
        write(00,*) "Automatic mapping check: Pr =", real(viscosity/diffusivity,kind=8), &
            "Mach =", real(velocityUnit/dsqrt(cs2),kind=8)
        write(00,*) "itc_max =", itc_max, "checkInterval =", checkInterval
        write(00,*) "Pressure convention: p = cs2*(rho-rho0), i.e. the weakly compressible pressure fluctuation."
        write(00,*) "Boundary warning: the selected collision with halfway BB/ABB needs independent wall-accuracy checks."
        write(00,*) "Reference DNS target: 2-D RB, Pr=0.7, Ra=1e10, 2049^2 boundary-including grid points."
        write(00,*) "nu/kappa/gBeta are derived from Ra/Pr/Ma/H_LB; the selected base-tau policy determines chi."
#ifdef RayleighBenardCell
        write(00,*) "RB initial perturbation amplitude / DeltaT =", real(rbPerturbation,kind=8)
        write(00,*) "RB pressure initialization = hydrostatic balance of the conductive profile"
#endif
#ifdef unsteadyFlow
        write(00,*) "Transient/free-fall time =", real(transientTimeTf,kind=8)
        write(00,*) "Averaging/free-fall time =", real(averagingTimeTf,kind=8)
        write(00,*) "Sample interval/free-fall time =", real(sampleTimeTf,kind=8)
        write(00,*) "Average start itc =", averageStartItc, "sample interval itc =", sampleIntervalItc
        write(00,*) "Target statistic samples =", targetStatisticSamples
        write(00,*) "Snapshot output enabled =", outputSnapshotFile
        write(00,*) "Snapshot interval/free-fall time =", real(snapshotTimeTf,kind=8)
        write(00,*) "Target snapshot files in statistics window =", targetSnapshotFiles
        write(00,*) "Tecplot output enabled =", outputTecplotFile
        write(00,*) "Tecplot interval/free-fall time =", real(tecplotTimeTf,kind=8)
        write(00,*) "Target periodic Tecplot files in statistics window =", targetTecplotFiles
#endif
        close(00)

        open(unit=01, file=trim(historyFile), status='replace', action='write')
        write(01,'(A)') "# itc time_ff errorU errorT Nu_global Nu_hot Nu_cold Nu_middle " // &
            "Nu_hot_halfway Nu_cold_halfway Nu_middle_fd ReVolAvg " // &
            "Nu_eps_u Nu_eps_T eps_u_exact_ratio eps_T_exact_ratio maxMach minT maxT maxAbsRhoDeviation samples"
        close(01)
#ifdef unsteadyFlow
        if(outputSnapshotFile.EQ.1) then
            open(unit=01, file=trim(snapshotFilePrefix)//"-readme", status='replace', action='write')
            write(01,'(A)') "Fortran sequential unformatted float64 records: u/U_ff, v/U_ff, T, rho."
            write(01,'(A,I0,A,I0)') "Each record contains a field with I=", nx, ", J=", ny
            write(01,'(A)') "Files are numbered by output order; see the index file for itc and t/t_ff."
            write(01,'(A)') "These snapshots are for post-processing, not strict f/g restart."
            close(01)
            open(unit=01, file=trim(snapshotIndexFile), status='replace', action='write')
            write(01,'(A)') "# snapshot_id itc time_ff filename"
            close(01)
        endif
#endif
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

        !$acc update self(u,v,T,rho,pressure,gradTx,gradTy,Sxx,Sxy,Syy)
    end subroutine update_host_monitor_2d_openacc
!===================================================================================================


!===================================================================================================
! 子程序: update_host_snapshot_2d_openacc
! 作用: 只同步二进制瞬时场需要的 u/v/T/rho，避免每 0.5 t_ff 传回应变率和梯度等额外数组。
!===================================================================================================
    subroutine update_host_snapshot_2d_openacc()
        use commondata
        implicit none

        !$acc update self(u,v,T,rho)
    end subroutine update_host_snapshot_2d_openacc
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
! 这里的 p 是压力扰动 p'，不是绝对热力学压力；这也统一了 geq、Psi 和 ABB 中压力项的含义。
!===================================================================================================
    subroutine compute_force()
        use commondata
        implicit none
        integer(kind=4) :: i, j

        !$acc parallel loop collapse(2) present(Fx,Fy,T,rho,pressure) private(i,j)
        do j = 1, ny
            do i = 1, nx
                pressure(i,j) = cs2*(rho(i,j)-rho0)       ! LBM-CDE 使用的弱可压缩压力扰动 p'
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
! 审稿勘误: 已发表式 (33) 的迹/散度项前误印为负号；由式 (30) 回代后应为正号。
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
                Sxx(i,j) = -(neqxx+u(i,j)*Fx(i,j))/denomDiag + coeffTrace*Sdiv(i,j)
                Syy(i,j) = -(neqyy+v(i,j)*Fy(i,j))/denomDiag + coeffTrace*Sdiv(i,j)
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
                do alpha = 0, thermalQ
                    geqLoc = calc_geq(alpha, T(i,j), u(i,j), v(i,j), rho(i,j))
                    neqx = neqx + dble(qex(alpha))*(g(alpha,i,j)-geqLoc)
                    neqy = neqy + dble(qey(alpha))*(g(alpha,i,j)-geqLoc)
                enddo
                ! 式 (35) 的分母。pressure/rho0 很小时仍保留，以保持与论文的压力耦合一致。
                denom = cT2*(2.0d0*taugL*(1.0d0-chi_kappa)+chi_kappa) + pressure(i,j)/rho0
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
! 作用: 执行流场 TRT 碰撞，并加入论文式 (26) 的 Guo 力项和 chi_s/chi_b 应变率修正项。
! 每一对反向速度分别采用偶/奇松弛率；原始源 R 也先分解成 R+/R-，再分别乘
! (1-s_even/2) 与 (1-s_odd/2)。应力修正是偶源，浮力项同时含偶、奇部分。
!===================================================================================================
    subroutine collision()
        use commondata
        implicit none
        integer(kind=4) :: i, j, alpha, beta, pair
        real(kind=8) :: euA, euB, eFA, eFB, uF, feqA, feqB
        real(kind=8) :: Axx, Axy, Ayy, hermiteA, hermiteB
        real(kind=8) :: rawA, rawB, dEven, dOdd, rEven, rOdd
        real(kind=8), parameter :: sourcePrefEven=1.0d0-0.5d0*omegaF
        real(kind=8), parameter :: sourcePrefOdd=1.0d0-0.5d0*omegaFOdd

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
        !$acc& private(i,j,alpha,beta,pair,euA,euB,eFA,eFB,uF,feqA,feqB,Axx,Axy,Ayy, &
        !$acc& hermiteA,hermiteB,rawA,rawB,dEven,dOdd,rEven,rOdd)
        do j = 1, ny
            do i = 1, nx
                uF = u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j)
                ! A_ab = chi_s*S_ab + (chi_b-chi_s)*Sdiv*delta_ab/D，二维 D=2。
                ! 审稿勘误说明：连续式 (12) 的体黏度 Hermite 系数应含 1/2；这里采用最终式 (26)
                ! 的 A_ab，并由 D2Q9 二阶矩对称收缩产生正确系数，不能再额外乘或除一个 2。
                Axx = chi_s*Sxx(i,j)+0.5d0*(chi_b-chi_s)*Sdiv(i,j)
                Ayy = chi_s*Syy(i,j)+0.5d0*(chi_b-chi_s)*Sdiv(i,j)
                Axy = chi_s*Sxy(i,j)
                ! 静止方向只含偶模。
                alpha = 0
                feqA = calc_feq(alpha, rho(i,j), u(i,j), v(i,j))
                hermiteA = -Axx-Ayy
                rawA = qomega(alpha)*(-uF/cs2+rho0*hermiteA)
                f_post(alpha,i,j) = f(alpha,i,j)-omegaF*(f(alpha,i,j)-feqA)+sourcePrefEven*rawA

                ! 四对反向速度：(1,3)、(2,4)、(5,7)、(6,8)。
                do pair = 1, 4
                    select case(pair)
                    case(1)
                        alpha = 1
                        beta = 3
                    case(2)
                        alpha = 2
                        beta = 4
                    case(3)
                        alpha = 5
                        beta = 7
                    case default
                        alpha = 6
                        beta = 8
                    end select

                    euA = dble(qex(alpha))*u(i,j)+dble(qey(alpha))*v(i,j)
                    euB = dble(qex(beta))*u(i,j)+dble(qey(beta))*v(i,j)
                    eFA = dble(qex(alpha))*Fx(i,j)+dble(qey(alpha))*Fy(i,j)
                    eFB = dble(qex(beta))*Fx(i,j)+dble(qey(beta))*Fy(i,j)
                    feqA = calc_feq(alpha, rho(i,j), u(i,j), v(i,j))
                    feqB = calc_feq(beta, rho(i,j), u(i,j), v(i,j))
                    hermiteA = (dble(qex(alpha)*qex(alpha))/cs2-1.0d0)*Axx + &
                        (dble(qey(alpha)*qey(alpha))/cs2-1.0d0)*Ayy + &
                        2.0d0*dble(qex(alpha)*qey(alpha))/cs2*Axy
                    hermiteB = (dble(qex(beta)*qex(beta))/cs2-1.0d0)*Axx + &
                        (dble(qey(beta)*qey(beta))/cs2-1.0d0)*Ayy + &
                        2.0d0*dble(qex(beta)*qey(beta))/cs2*Axy
                    rawA = qomega(alpha)*(eFA/cs2+euA*eFA/(cs2*cs2)-uF/cs2+rho0*hermiteA)
                    rawB = qomega(beta)*(eFB/cs2+euB*eFB/(cs2*cs2)-uF/cs2+rho0*hermiteB)

                    dEven = 0.5d0*((f(alpha,i,j)-feqA)+(f(beta,i,j)-feqB))
                    dOdd = 0.5d0*((f(alpha,i,j)-feqA)-(f(beta,i,j)-feqB))
                    rEven = 0.5d0*(rawA+rawB)
                    rOdd = 0.5d0*(rawA-rawB)
                    f_post(alpha,i,j) = f(alpha,i,j)-omegaF*dEven-omegaFOdd*dOdd + &
                        sourcePrefEven*rEven+sourcePrefOdd*rOdd
                    f_post(beta,i,j) = f(beta,i,j)-omegaF*dEven+omegaFOdd*dOdd + &
                        sourcePrefEven*rEven-sourcePrefOdd*rOdd
                enddo
            enddo
        enddo
    end subroutine collision
!===================================================================================================


!===================================================================================================
! 子程序: collisionT
! 作用: 严格执行原 LBM-CDE 的 D2Q9-BGK 温度碰撞，保留式 (28) 的
! Q、p*grad(T)+T*F 与 chi_kappa*grad(T) 三部分，不混入 Xu/Luo 或 D2Q5-MRT 分支。
!===================================================================================================
    subroutine collisionT()
        use commondata
        implicit none
        integer(kind=4) :: i, j, alpha
        real(kind=8) :: geqLoc, scalarSource, psi, eu
        real(kind=8), parameter :: sourcePref=1.0d0-0.5d0*omegaG

        !$acc parallel loop collapse(3) present(g_post) private(i,j,alpha)
        do j = 0, ny+1
            do i = 0, nx+1
                do alpha = 0, thermalQ
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
                    ! 大括号内逐项对应论文式 (28)，外层温度权重与梯形前因子在下一行乘入。
                    scalarSource = heatSourceQ*(1.0d0+eu/cT2) + &
                        (dble(qex(alpha))*(pressure(i,j)*gradTx(i,j)+T(i,j)*Fx(i,j)) + &
                        dble(qey(alpha))*(pressure(i,j)*gradTy(i,j)+T(i,j)*Fy(i,j)))/(rho0*cT2) + &
                        chi_kappa*(dble(qex(alpha))*gradTx(i,j)+dble(qey(alpha))*gradTy(i,j))
                    psi = sourcePref*tomega(alpha)*scalarSource
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
! 作用: 温度 D2Q9 分布函数执行拉式迁移。
!===================================================================================================
    subroutine streamingT()
        use commondata
        implicit none
        integer(kind=4) :: i, j, alpha, ip, jp

        !$acc parallel loop collapse(2) present(g,g_post) private(i,j,alpha,ip,jp)
        do j = 1, ny
            do i = 1, nx
                do alpha = 0, thermalQ
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
! 作用: 恒温壁采用 anti-bounce-back（论文式 (38)），绝热壁采用 bounce-back（式 (39)）。
! 关键点: 式 (38) 不只有 2*omega_i*Twall，还必须保留压力平衡项；分母为 rho0*cs2。
! 边界精度: 这些规则反射的是含完整源项的碰撞后分布；但 BGK 内点与 halfway BB/ABB 在壁面
! 仍有隐藏截断误差，因此本实现只保证边界语义和质量/热流方向一致，不宣称全局二阶精度。
!===================================================================================================
    subroutine bouncebackT()
        use commondata
        implicit none
        integer(kind=4) :: i, j

#ifdef VerticalWallsConstT
        ! 左热壁 Thot；右冷壁 Tcold。每个赋值只更新从壁面进入流体的未知分布。
        !$acc parallel loop present(g,g_post,pressure) private(j)
        do j = 1, ny
            g(1,1,j) = -g_post(3,1,j)+2.0d0*tomega(1)*Thot + &
                2.0d0*qlambdaT(1)*Thot*pressure(1,j)/(rho0*cT2)
#ifdef ThermalD2Q9
            g(5,1,j) = -g_post(7,1,j)+2.0d0*tomega(5)*Thot + &
                2.0d0*qlambdaT(5)*Thot*pressure(1,j)/(rho0*cT2)
            g(8,1,j) = -g_post(6,1,j)+2.0d0*tomega(8)*Thot + &
                2.0d0*qlambdaT(8)*Thot*pressure(1,j)/(rho0*cT2)
#endif

            g(3,nx,j) = -g_post(1,nx,j)+2.0d0*tomega(3)*Tcold + &
                2.0d0*qlambdaT(3)*Tcold*pressure(nx,j)/(rho0*cT2)
#ifdef ThermalD2Q9
            g(6,nx,j) = -g_post(8,nx,j)+2.0d0*tomega(6)*Tcold + &
                2.0d0*qlambdaT(6)*Tcold*pressure(nx,j)/(rho0*cT2)
            g(7,nx,j) = -g_post(5,nx,j)+2.0d0*tomega(7)*Tcold + &
                2.0d0*qlambdaT(7)*Tcold*pressure(nx,j)/(rho0*cT2)
#endif
        enddo
#endif

#ifdef VerticalWallsAdiabatic
        ! RB 左右绝热壁：零法向热通量。
        !$acc parallel loop present(g,g_post) private(j)
        do j = 1, ny
            g(1,1,j) = g_post(3,1,j)
#ifdef ThermalD2Q9
            g(5,1,j) = g_post(7,1,j)
            g(8,1,j) = g_post(6,1,j)
#endif

            g(3,nx,j) = g_post(1,nx,j)
#ifdef ThermalD2Q9
            g(6,nx,j) = g_post(8,nx,j)
            g(7,nx,j) = g_post(5,nx,j)
#endif
        enddo
#endif

#ifdef HorizontalWallsConstT
        ! RB 下热壁 Thot；上冷壁 Tcold。
        !$acc parallel loop present(g,g_post,pressure) private(i)
        do i = 1, nx
            g(2,i,1) = -g_post(4,i,1)+2.0d0*tomega(2)*Thot + &
                2.0d0*qlambdaT(2)*Thot*pressure(i,1)/(rho0*cT2)
#ifdef ThermalD2Q9
            g(5,i,1) = -g_post(7,i,1)+2.0d0*tomega(5)*Thot + &
                2.0d0*qlambdaT(5)*Thot*pressure(i,1)/(rho0*cT2)
            g(6,i,1) = -g_post(8,i,1)+2.0d0*tomega(6)*Thot + &
                2.0d0*qlambdaT(6)*Thot*pressure(i,1)/(rho0*cT2)
#endif

            g(4,i,ny) = -g_post(2,i,ny)+2.0d0*tomega(4)*Tcold + &
                2.0d0*qlambdaT(4)*Tcold*pressure(i,ny)/(rho0*cT2)
#ifdef ThermalD2Q9
            g(7,i,ny) = -g_post(5,i,ny)+2.0d0*tomega(7)*Tcold + &
                2.0d0*qlambdaT(7)*Tcold*pressure(i,ny)/(rho0*cT2)
            g(8,i,ny) = -g_post(6,i,ny)+2.0d0*tomega(8)*Tcold + &
                2.0d0*qlambdaT(8)*Tcold*pressure(i,ny)/(rho0*cT2)
#endif
        enddo
#endif

#ifdef HorizontalWallsAdiabatic
        ! 零法向热通量：把撞向壁面的分布原值反弹，不添加恒温修正。
        !$acc parallel loop present(g,g_post) private(i)
        do i = 1, nx
            g(2,i,1) = g_post(4,i,1)
#ifdef ThermalD2Q9
            g(5,i,1) = g_post(7,i,1)
            g(6,i,1) = g_post(8,i,1)
#endif

            g(4,i,ny) = g_post(2,i,ny)
#ifdef ThermalD2Q9
            g(7,i,ny) = g_post(5,i,ny)
            g(8,i,ny) = g_post(6,i,ny)
#endif
        enddo
#endif

#if defined(VerticalWallsConstT) && defined(HorizontalWallsAdiabatic) && defined(ThermalD2Q9)
        ! 四个角点的对角分布同时属于恒温侧壁和绝热水平壁。
        ! 上面的水平壁循环会覆盖侧壁先写入的对角分布，因此这里最后再写一次，使恒温 Dirichlet 条件优先。
        ! 角点只有 4 个标量赋值，使用 serial 比额外启动并行 kernel 更直接。
        !$acc serial present(g,g_post,pressure)
        g(5,1,1) = -g_post(7,1,1)+2.0d0*tomega(5)*Thot + &
            2.0d0*qlambdaT(5)*Thot*pressure(1,1)/(rho0*cT2)
        g(8,1,ny) = -g_post(6,1,ny)+2.0d0*tomega(8)*Thot + &
            2.0d0*qlambdaT(8)*Thot*pressure(1,ny)/(rho0*cT2)
        g(6,nx,1) = -g_post(8,nx,1)+2.0d0*tomega(6)*Tcold + &
            2.0d0*qlambdaT(6)*Tcold*pressure(nx,1)/(rho0*cT2)
        g(7,nx,ny) = -g_post(5,nx,ny)+2.0d0*tomega(7)*Tcold + &
            2.0d0*qlambdaT(7)*Tcold*pressure(nx,ny)/(rho0*cT2)
        !$acc end serial
#endif

#if defined(HorizontalWallsConstT) && defined(VerticalWallsAdiabatic) && defined(ThermalD2Q9)
        ! RB 角点同时属于恒温水平壁和绝热侧壁；最后重写对角入射分布，使 Dirichlet 条件优先。
        !$acc serial present(g,g_post,pressure)
        g(5,1,1) = -g_post(7,1,1)+2.0d0*tomega(5)*Thot + &
            2.0d0*qlambdaT(5)*Thot*pressure(1,1)/(rho0*cT2)
        g(6,nx,1) = -g_post(8,nx,1)+2.0d0*tomega(6)*Thot + &
            2.0d0*qlambdaT(6)*Thot*pressure(nx,1)/(rho0*cT2)
        g(8,1,ny) = -g_post(6,1,ny)+2.0d0*tomega(8)*Tcold + &
            2.0d0*qlambdaT(8)*Tcold*pressure(1,ny)/(rho0*cT2)
        g(7,nx,ny) = -g_post(5,nx,ny)+2.0d0*tomega(7)*Tcold + &
            2.0d0*qlambdaT(7)*Tcold*pressure(nx,ny)/(rho0*cT2)
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
        integer(kind=4) :: i, j, alpha
        real(kind=8) :: temperatureSum

        !$acc parallel loop collapse(2) present(g,T) private(i,j,alpha,temperatureSum)
        do j = 1, ny
            do i = 1, nx
                temperatureSum = 0.0d0
                do alpha = 0, thermalQ
                    temperatureSum = temperatureSum+g(alpha,i,j)
                enddo
                T(i,j) = temperatureSum+0.5d0*heatSourceQ
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
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        implicit none
        integer(kind=4) :: i, j, badI, badJ
        real(kind=8) :: errUSum, refUSum, errTSum, refTSum
        logical :: badState

        errUSum = 0.0d0
        refUSum = 0.0d0
        errTSum = 0.0d0
        refTSum = 0.0d0
        badState = .false.
        badI = -1
        badJ = -1
        do j = 1, ny
            do i = 1, nx
                if((.NOT.ieee_is_finite(u(i,j))).OR.(.NOT.ieee_is_finite(v(i,j))).OR. &
                    (.NOT.ieee_is_finite(T(i,j))).OR.(.NOT.ieee_is_finite(rho(i,j))).OR. &
                    (.NOT.ieee_is_finite(gradTx(i,j))).OR.(.NOT.ieee_is_finite(gradTy(i,j)))) then
                    if(.NOT.badState) then
                        badState = .true.
                        badI = i
                        badJ = j
                    endif
                    cycle
                endif
                errUSum = errUSum+(u(i,j)-up(i,j))**2+(v(i,j)-vp(i,j))**2
                refUSum = refUSum+u(i,j)**2+v(i,j)**2
                errTSum = errTSum+(T(i,j)-Tp(i,j))**2
                refTSum = refTSum+T(i,j)**2
            enddo
        enddo
        if(badState) then
            open(unit=00, file=trim(settingsFile), status='old', position='append', action='write')
            write(00,*) "NONFINITE state detected: itc/i/j =", itc, badI, badJ
            write(00,*) "u/v/T/rho/gradTx/gradTy =", u(badI,badJ), v(badI,badJ), T(badI,badJ), &
                rho(badI,badJ), gradTx(badI,badJ), gradTy(badI,badJ)
            close(00)
            error stop "Nonfinite LBM state; see SimulationSettings2DOpenaccLBMCDE.txt"
        endif
        errorU = dsqrt(errUSum/max(refUSum,tiny(1.0d0)))
        errorT = dsqrt(errTSum/max(refTSum,tiny(1.0d0)))
        up = u
        vp = v
        Tp = T
        lastCheckItc = itc

        write(*,'(A,I10,2(A,ES14.6))') "itc=", itc, " errorU=", errorU, " errorT=", errorT
    end subroutine check
!===================================================================================================


!===================================================================================================
! 子程序: calNuRe
! 作用: 在 CPU 端计算热/冷壁 Nu、中心截面总热流、体平均 Nu/rms Re、耗散恒等式诊断。
! 侧热工况主热流方向为 x；RB 主热流方向为 y。Numax_mid 仍只用导热梯度，便于与论文表 2 对照；
! Nu_middle 则包含对流和导热两部分，用于截面热流守恒检查。
! 注意: 偶数网格的中心线位于两条流体节点线之间，因此对相邻两条节点线做线性插值。
! 壁面 Nu 同时给出局部非平衡梯度值与半格壁面温差值，便于交叉检查梯度重构。
!===================================================================================================
    subroutine calNuRe()
        use commondata
        implicit none
        integer(kind=4) :: i, j, imidLeft, imidRight, jmidBottom, jmidTop
        real(kind=8) :: sumHot, sumCold, sumMiddle, sumRe, sumGlobal
        real(kind=8) :: sumHotHalfway, sumColdHalfway, sumMiddleFD
        real(kind=8) :: sumEpsKinetic, sumEpsThermal, cellCount, speedSquared
        real(kind=8) :: localNu, centerGradNormal, centerGradNormalFD, centerU, centerV, centerT

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
        sumHotHalfway = 0.0d0
        sumColdHalfway = 0.0d0
        sumMiddleFD = 0.0d0
        sumRe = 0.0d0
        sumGlobal = 0.0d0
        sumEpsKinetic = 0.0d0
        sumEpsThermal = 0.0d0
        cellCount = dble(nx)*dble(ny)
        maxMachLocal = 0.0d0
        minTemperature = huge(1.0d0)
        maxTemperature = -huge(1.0d0)
        maxDensityDeviation = 0.0d0
        umax_mid = -huge(1.0d0)
        vmax_mid = -huge(1.0d0)
        Numax_mid = -huge(1.0d0)
        pos_umax_mid = 0.0d0
        pos_vmax_mid = 0.0d0
        pos_Numax_mid = 0.0d0

#ifdef SideHeatedCell
        do j = 1, ny
            ! gradTx 由论文式 (35) 在局部节点直接得到；壁面值用相邻的第一/最后一层流体节点代表。
            localNu = -lengthUnit*gradTx(1,j)/deltaT
            sumHot = sumHot + localNu
            localNu = -lengthUnit*gradTx(nx,j)/deltaT
            sumCold = sumCold + localNu
            sumHotHalfway = sumHotHalfway+2.0d0*lengthUnit*(Thot-T(1,j))/deltaT
            sumColdHalfway = sumColdHalfway+2.0d0*lengthUnit*(T(nx,j)-Tcold)/deltaT

            ! x=0.5L 中心线：偶数网格平均中心左右两列，奇数网格两个下标相同，公式自动退化为单列值。
            centerGradNormal = 0.5d0*(gradTx(imidLeft,j)+gradTx(imidRight,j))
            if(imidLeft.NE.imidRight) then
                centerGradNormalFD = T(imidRight,j)-T(imidLeft,j)
            else
                centerGradNormalFD = 0.5d0*(T(imidLeft+1,j)-T(imidLeft-1,j))
            endif
            centerU = 0.5d0*(u(imidLeft,j)+u(imidRight,j))
            centerT = 0.5d0*(T(imidLeft,j)+T(imidRight,j))
            localNu = -lengthUnit*centerGradNormal/deltaT
            sumMiddle = sumMiddle + lengthUnit*(centerU*(centerT-Tref)-diffusivity*centerGradNormal)/ &
                (diffusivity*deltaT)
            sumMiddleFD = sumMiddleFD + lengthUnit*(centerU*(centerT-Tref)-diffusivity*centerGradNormalFD)/ &
                (diffusivity*deltaT)
            if(localNu.GT.Numax_mid) then
                Numax_mid = localNu
                pos_Numax_mid = yp(j)
            endif
            ! 论文表 2 的 umax 是“最大正值”，不是绝对值最大的带符号速度。
            if(centerU.GT.umax_mid) then
                umax_mid = centerU
                pos_umax_mid = yp(j)
            endif
        enddo

        do i = 1, nx
            ! y=0.5H 中心线同理在上下两行之间插值。
            centerV = 0.5d0*(v(i,jmidBottom)+v(i,jmidTop))
            ! 论文表 2 的 vmax 同样取最大正值；其位置应落在靠近左侧热壁的一侧。
            if(centerV.GT.vmax_mid) then
                vmax_mid = centerV
                pos_vmax_mid = xp(i)
            endif
        enddo
#endif

#ifdef RayleighBenardCell
        do i = 1, nx
            ! RB 底热/顶冷壁的无量纲法向导热通量；两处梯度均为负，因此 Nu 取正值。
            localNu = -lengthUnit*gradTy(i,1)/deltaT
            sumHot = sumHot + localNu
            localNu = -lengthUnit*gradTy(i,ny)/deltaT
            sumCold = sumCold + localNu
            sumHotHalfway = sumHotHalfway+2.0d0*lengthUnit*(Thot-T(i,1))/deltaT
            sumColdHalfway = sumColdHalfway+2.0d0*lengthUnit*(T(i,ny)-Tcold)/deltaT

            ! y=0.5H 截面总热流；Numax_mid 只记录其中的导热梯度局部最大值。
            centerGradNormal = 0.5d0*(gradTy(i,jmidBottom)+gradTy(i,jmidTop))
            if(jmidBottom.NE.jmidTop) then
                centerGradNormalFD = T(i,jmidTop)-T(i,jmidBottom)
            else
                centerGradNormalFD = 0.5d0*(T(i,jmidBottom+1)-T(i,jmidBottom-1))
            endif
            centerV = 0.5d0*(v(i,jmidBottom)+v(i,jmidTop))
            centerU = 0.5d0*(u(i,jmidBottom)+u(i,jmidTop))
            centerT = 0.5d0*(T(i,jmidBottom)+T(i,jmidTop))
            localNu = -lengthUnit*centerGradNormal/deltaT
            sumMiddle = sumMiddle + lengthUnit*(centerV*(centerT-Tref)-diffusivity*centerGradNormal)/ &
                (diffusivity*deltaT)
            sumMiddleFD = sumMiddleFD + lengthUnit*(centerV*(centerT-Tref)-diffusivity*centerGradNormalFD)/ &
                (diffusivity*deltaT)
            if(localNu.GT.Numax_mid) then
                Numax_mid = localNu
                pos_Numax_mid = xp(i)
            endif
            if(centerU.GT.umax_mid) then
                umax_mid = centerU
                pos_umax_mid = xp(i)
            endif
        enddo

        do j = 1, ny
            centerV = 0.5d0*(v(imidLeft,j)+v(imidRight,j))
            if(centerV.GT.vmax_mid) then
                vmax_mid = centerV
                pos_vmax_mid = yp(j)
            endif
        enddo
#endif

        do j = 1, ny
            do i = 1, nx
                speedSquared = u(i,j)*u(i,j)+v(i,j)*v(i,j)
                sumRe = sumRe+speedSquared
                maxMachLocal = max(maxMachLocal,dsqrt(speedSquared/cs2))
                minTemperature = min(minTemperature,T(i,j))
                maxTemperature = max(maxTemperature,T(i,j))
                maxDensityDeviation = max(maxDensityDeviation,abs(rho(i,j)-rho0))
                sumEpsKinetic = sumEpsKinetic+2.0d0*viscosity* &
                    (Sxx(i,j)*Sxx(i,j)+Syy(i,j)*Syy(i,j)+2.0d0*Sxy(i,j)*Sxy(i,j))
                sumEpsThermal = sumEpsThermal+diffusivity* &
                    (gradTx(i,j)*gradTx(i,j)+gradTy(i,j)*gradTy(i,j))
#ifdef SideHeatedCell
                sumGlobal = sumGlobal + u(i,j)*(T(i,j)-Tref)
#endif
#ifdef RayleighBenardCell
                sumGlobal = sumGlobal + v(i,j)*(T(i,j)-Tref)
#endif
            enddo
        enddo

#ifdef SideHeatedCell
        Nu_hot = sumHot/dble(ny)
        Nu_cold = sumCold/dble(ny)
        Nu_middle = sumMiddle/dble(ny)
        Nu_hot_halfway = sumHotHalfway/dble(ny)
        Nu_cold_halfway = sumColdHalfway/dble(ny)
        Nu_middleFD = sumMiddleFD/dble(ny)
#endif
#ifdef RayleighBenardCell
        Nu_hot = sumHot/dble(nx)
        Nu_cold = sumCold/dble(nx)
        Nu_middle = sumMiddle/dble(nx)
        Nu_hot_halfway = sumHotHalfway/dble(nx)
        Nu_cold_halfway = sumColdHalfway/dble(nx)
        Nu_middleFD = sumMiddleFD/dble(nx)
#endif
        ! 在线性导热基准上加入体平均对流热通量；主速度分量由工况分支在上面选择。
        Nu_global = 1.0d0+sumGlobal/cellCount*lengthUnit/(diffusivity*deltaT)
        speedSquaredVol = sumRe/cellCount
        ReVolAvg = dsqrt(speedSquaredVol)*lengthUnit/viscosity
        epsKineticVol = sumEpsKinetic/cellCount
        epsThermalVol = sumEpsThermal/cellCount
#ifdef RayleighBenardCell
        ! 2-D RB 的精确全局耗散关系：eps_T=kappa*DeltaT^2/H^2*Nu，
        ! eps_u=nu^3/H^4*(Nu-1)*Ra/Pr^2。偏差同时反映统计不足、边界和离散误差。
        NuThermalDiss = epsThermalVol*lengthUnit**2/(diffusivity*deltaT**2)
        NuKineticDiss = 1.0d0+epsKineticVol*lengthUnit**4*Prandtl**2/ &
            (viscosity**3*Rayleigh)
        if(abs(Nu_global-1.0d0).GT.1.0d-12) then
            kineticDissExactRatio = (NuKineticDiss-1.0d0)/(Nu_global-1.0d0)
        else
            kineticDissExactRatio = 0.0d0
        endif
        if(abs(Nu_global).GT.1.0d-12) then
            thermalDissExactRatio = NuThermalDiss/Nu_global
        else
            thermalDissExactRatio = 0.0d0
        endif
#endif
#ifdef SideHeatedCell
        ! 上述精确耗散关系只适用于底热顶冷 RB；侧热方腔不伪造对应的等效 Nu 或比值。
        NuKineticDiss = 0.0d0
        NuThermalDiss = 0.0d0
        kineticDissExactRatio = 0.0d0
        thermalDissExactRatio = 0.0d0
#endif
        ! 速度极值采用自由落体速度 u0 归一化；Re 使用 rms 全速度和黏性标度。
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
        write(11,'(I12,1X,19(ES24.16E3,1X),I10)') itc, dble(itc)/timeUnit, errorU, errorT, &
            Nu_global, Nu_hot, Nu_cold, Nu_middle, Nu_hot_halfway, Nu_cold_halfway, Nu_middleFD, &
            ReVolAvg, NuKineticDiss, NuThermalDiss, &
            kineticDissExactRatio, thermalDissExactRatio, maxMachLocal, minTemperature, &
            maxTemperature, maxDensityDeviation, statisticSamples
        close(11)

        open(unit=00, file=trim(settingsFile), status='old', position='append', action='write')
        write(00,'(A,I12,2(A,ES16.8))') "itc=", itc, " errorU=", errorU, " errorT=", errorT
        write(00,'(5(A,ES16.8),A,I8)') " Nu_global=", Nu_global, " Nu_hot=", Nu_hot, &
            " Nu_cold=", Nu_cold, " Nu_middle=", Nu_middle, " Re=", ReVolAvg, &
            " samples=", statisticSamples
        write(00,'(3(A,ES16.8))') " Nu_hot_halfway=", Nu_hot_halfway, &
            " Nu_cold_halfway=", Nu_cold_halfway, " Nu_middle_fd=", Nu_middleFD
        write(00,'(6(A,ES16.8))') " Nu_eps_u=", NuKineticDiss, " Nu_eps_T=", NuThermalDiss, &
            " maxMach=", maxMachLocal, " minT=", minTemperature, " maxT=", maxTemperature, &
            " maxAbsRhoDeviation=", maxDensityDeviation
        write(00,'(2(A,ES16.8))') " eps_u/exact=", kineticDissExactRatio, &
            " eps_T/exact=", thermalDissExactRatio
        close(00)
    end subroutine append_history
!===================================================================================================


!===================================================================================================
! 子程序: accumulate_statistics
! 作用: 对高 Ra 非稳态阶段的瞬时 Nu、体积均方速度、耗散量和水平温度剖面做在线累积。
! Re 必须最后由 <u^2+v^2>_(V,t) 开方，不能先对每个时刻的 Re 开方再做时间平均。
!===================================================================================================
    subroutine accumulate_statistics()
        use commondata
        implicit none
        integer(kind=4) :: i, j
        real(kind=8) :: temperatureLineSum, temperatureSquaredLineSum

        statisticSamples = statisticSamples+1
        if(firstStatisticItc.LT.0) firstStatisticItc = itc
        lastStatisticItc = itc
        sumNuGlobal = sumNuGlobal+Nu_global
        sumNuGlobal2 = sumNuGlobal2+Nu_global*Nu_global
        sumNuHot = sumNuHot+Nu_hot
        sumNuHot2 = sumNuHot2+Nu_hot*Nu_hot
        sumNuCold = sumNuCold+Nu_cold
        sumNuCold2 = sumNuCold2+Nu_cold*Nu_cold
        sumNuMiddle = sumNuMiddle+Nu_middle
        sumNuMiddle2 = sumNuMiddle2+Nu_middle*Nu_middle
        sumSpeedSquaredVol = sumSpeedSquaredVol+speedSquaredVol
        sumSpeedSquaredVol2 = sumSpeedSquaredVol2+speedSquaredVol*speedSquaredVol
        sumNuKineticDiss = sumNuKineticDiss+NuKineticDiss
        sumNuKineticDiss2 = sumNuKineticDiss2+NuKineticDiss*NuKineticDiss
        sumNuThermalDiss = sumNuThermalDiss+NuThermalDiss
        sumNuThermalDiss2 = sumNuThermalDiss2+NuThermalDiss*NuThermalDiss
        maxStatisticCFL = max(maxStatisticCFL,maxMachLocal*dsqrt(cs2))

#ifdef unsteadyFlow
        if(statisticSamples.LE.max(1,targetStatisticSamples/2)) then
            firstHalfSamples = firstHalfSamples+1
            sumNuGlobalFirstHalf = sumNuGlobalFirstHalf+Nu_global
            sumSpeedSquaredFirstHalf = sumSpeedSquaredFirstHalf+speedSquaredVol
        else
            secondHalfSamples = secondHalfSamples+1
            sumNuGlobalSecondHalf = sumNuGlobalSecondHalf+Nu_global
            sumSpeedSquaredSecondHalf = sumSpeedSquaredSecondHalf+speedSquaredVol
        endif
#endif

        ! 论文用温度 rms 峰值到壁面的距离定义热边界层厚度 delta_theta_rms。
        do j = 1, ny
            temperatureLineSum = 0.0d0
            temperatureSquaredLineSum = 0.0d0
            do i = 1, nx
                temperatureLineSum = temperatureLineSum+T(i,j)
                temperatureSquaredLineSum = temperatureSquaredLineSum+T(i,j)*T(i,j)
            enddo
            sumTemperatureProfile(j) = sumTemperatureProfile(j)+temperatureLineSum/dble(nx)
            sumTemperatureSquaredProfile(j) = sumTemperatureSquaredProfile(j)+ &
                temperatureSquaredLineSum/dble(nx)
        enddo
    end subroutine accumulate_statistics
!===================================================================================================


!===================================================================================================
! 子程序: write_statistics
! 作用: 输出非稳态样本的时间平均和总体标准差；无样本时明确写出原因，不伪造平均量。
!===================================================================================================
    subroutine write_statistics()
        use commondata
        implicit none
        real(kind=8) :: invN, meanGlobal, meanHot, meanCold, meanMiddle, meanRe
        real(kind=8) :: meanSpeedSquared, stdSpeedSquared
        real(kind=8) :: meanNuKineticDiss, meanNuThermalDiss
        real(kind=8) :: stdGlobal, stdHot, stdCold, stdMiddle
        real(kind=8) :: stdNuKineticDiss, stdNuThermalDiss
        real(kind=8) :: firstHalfNu, secondHalfNu, firstHalfRe, secondHalfRe
        real(kind=8) :: halfWindowNuRelativeDifference, halfWindowReRelativeDifference
#ifdef RayleighBenardCell
        integer(kind=4) :: j, lowerPeakIndex, upperPeakIndex, thermalBLGridPoints
        real(kind=8) :: rmsTemperature, lowerRmsMaximum, upperRmsMaximum
        real(kind=8) :: thermalBLThickness, kineticExactRatioMean, thermalExactRatioMean
        real(kind=8) :: meanEpsKinetic, meanEpsThermal, exactEpsKinetic, exactEpsThermal
        real(kind=8) :: kolmogorovScaleOverH, batchelorScaleOverH
        real(kind=8) :: gridOverKolmogorov, gridOverBatchelor, timeStepOverKolmogorovTime
#endif

        open(unit=21, file=trim(statisticsFile), status='replace', action='write')
        write(21,'(A)') "# Time statistics from the configured post-transient sampling window"
        write(21,'(A,I0)') "# samples ", statisticSamples
        write(21,'(A,2(1X,I0))') "# first_last_sample_itc", firstStatisticItc, lastStatisticItc
        write(21,'(A,2(1X,I0))') "# first_second_half_samples", firstHalfSamples, secondHalfSamples
        if(statisticSamples.GT.0) then
            invN = 1.0d0/dble(statisticSamples)
            meanGlobal = sumNuGlobal*invN
            meanHot = sumNuHot*invN
            meanCold = sumNuCold*invN
            meanMiddle = sumNuMiddle*invN
            meanSpeedSquared = sumSpeedSquaredVol*invN
            meanRe = dsqrt(max(0.0d0,meanSpeedSquared))*lengthUnit/viscosity
            meanNuKineticDiss = sumNuKineticDiss*invN
            meanNuThermalDiss = sumNuThermalDiss*invN
            stdGlobal = dsqrt(max(0.0d0,sumNuGlobal2*invN-meanGlobal*meanGlobal))
            stdHot = dsqrt(max(0.0d0,sumNuHot2*invN-meanHot*meanHot))
            stdCold = dsqrt(max(0.0d0,sumNuCold2*invN-meanCold*meanCold))
            stdMiddle = dsqrt(max(0.0d0,sumNuMiddle2*invN-meanMiddle*meanMiddle))
            stdSpeedSquared = dsqrt(max(0.0d0,sumSpeedSquaredVol2*invN-meanSpeedSquared*meanSpeedSquared))
            stdNuKineticDiss = dsqrt(max(0.0d0, &
                sumNuKineticDiss2*invN-meanNuKineticDiss*meanNuKineticDiss))
            stdNuThermalDiss = dsqrt(max(0.0d0, &
                sumNuThermalDiss2*invN-meanNuThermalDiss*meanNuThermalDiss))

            write(21,'(A)') "# metric mean population_std"
            write(21,'(A,2(1X,ES24.16E3))') "Nu_global", meanGlobal, stdGlobal
            write(21,'(A,2(1X,ES24.16E3))') "Nu_hot", meanHot, stdHot
            write(21,'(A,2(1X,ES24.16E3))') "Nu_cold", meanCold, stdCold
            write(21,'(A,2(1X,ES24.16E3))') "Nu_middle", meanMiddle, stdMiddle
            write(21,'(A,1X,ES24.16E3)') "Re_Table1_sqrt_volume_time_mean", meanRe
            write(21,'(A,2(1X,ES24.16E3))') "speed_squared_volume", meanSpeedSquared, stdSpeedSquared
            if((firstHalfSamples.GT.0).AND.(secondHalfSamples.GT.0)) then
                firstHalfNu = sumNuGlobalFirstHalf/dble(firstHalfSamples)
                secondHalfNu = sumNuGlobalSecondHalf/dble(secondHalfSamples)
                firstHalfRe = dsqrt(max(0.0d0,sumSpeedSquaredFirstHalf/dble(firstHalfSamples)))* &
                    lengthUnit/viscosity
                secondHalfRe = dsqrt(max(0.0d0,sumSpeedSquaredSecondHalf/dble(secondHalfSamples)))* &
                    lengthUnit/viscosity
                halfWindowNuRelativeDifference = abs(firstHalfNu-secondHalfNu)/ &
                    max(0.5d0*(abs(firstHalfNu)+abs(secondHalfNu)),tiny(1.0d0))
                halfWindowReRelativeDifference = abs(firstHalfRe-secondHalfRe)/ &
                    max(0.5d0*(abs(firstHalfRe)+abs(secondHalfRe)),tiny(1.0d0))
                write(21,'(A,2(1X,ES24.16E3))') "Nu_first_second_half", firstHalfNu, secondHalfNu
                write(21,'(A,2(1X,ES24.16E3))') "Re_first_second_half", firstHalfRe, secondHalfRe
                write(21,'(A,2(1X,ES24.16E3))') "half_window_relative_difference_Nu_Re", &
                    halfWindowNuRelativeDifference, halfWindowReRelativeDifference
            endif
#ifdef RayleighBenardCell
            write(21,'(A,2(1X,ES24.16E3))') "Nu_eps_u", meanNuKineticDiss, stdNuKineticDiss
            write(21,'(A,2(1X,ES24.16E3))') "Nu_eps_T", meanNuThermalDiss, stdNuThermalDiss

            if(meanGlobal.GT.1.0d0) then
                meanEpsKinetic = (meanNuKineticDiss-1.0d0)*viscosity**3*Rayleigh/ &
                    (lengthUnit**4*Prandtl**2)
                meanEpsThermal = meanNuThermalDiss*diffusivity*deltaT**2/lengthUnit**2
                exactEpsKinetic = viscosity**3/lengthUnit**4*(meanGlobal-1.0d0)*Rayleigh/Prandtl**2
                exactEpsThermal = diffusivity*deltaT**2/lengthUnit**2*meanGlobal
                kineticExactRatioMean = (meanNuKineticDiss-1.0d0)/(meanGlobal-1.0d0)
                thermalExactRatioMean = meanNuThermalDiss/meanGlobal
                kolmogorovScaleOverH = dsqrt(Prandtl)/(Rayleigh*(meanGlobal-1.0d0))**0.25d0
                batchelorScaleOverH = kolmogorovScaleOverH/dsqrt(Prandtl)
                gridOverKolmogorov = (1.0d0/lengthUnit)/kolmogorovScaleOverH
                gridOverBatchelor = (1.0d0/lengthUnit)/batchelorScaleOverH
                timeStepOverKolmogorovTime = (1.0d0/timeUnit)/dsqrt(Prandtl/(meanGlobal-1.0d0))
            else
                meanEpsKinetic = 0.0d0
                meanEpsThermal = 0.0d0
                exactEpsKinetic = 0.0d0
                exactEpsThermal = 0.0d0
                kineticExactRatioMean = 0.0d0
                thermalExactRatioMean = 0.0d0
                kolmogorovScaleOverH = huge(1.0d0)
                batchelorScaleOverH = huge(1.0d0)
                gridOverKolmogorov = 0.0d0
                gridOverBatchelor = 0.0d0
                timeStepOverKolmogorovTime = 0.0d0
            endif

            ! 论文第 3.3 节：delta_theta_rms 是温度 rms 剖面峰值到最近水平壁面的距离。
            lowerPeakIndex = 1
            upperPeakIndex = ny
            lowerRmsMaximum = -1.0d0
            upperRmsMaximum = -1.0d0
            do j = 1, max(1,ny/2)
                rmsTemperature = dsqrt(max(0.0d0,sumTemperatureSquaredProfile(j)*invN - &
                    (sumTemperatureProfile(j)*invN)**2))
                if(rmsTemperature.GT.lowerRmsMaximum) then
                    lowerRmsMaximum = rmsTemperature
                    lowerPeakIndex = j
                endif
            enddo
            do j = min(ny,ny/2+1), ny
                rmsTemperature = dsqrt(max(0.0d0,sumTemperatureSquaredProfile(j)*invN - &
                    (sumTemperatureProfile(j)*invN)**2))
                if(rmsTemperature.GT.upperRmsMaximum) then
                    upperRmsMaximum = rmsTemperature
                    upperPeakIndex = j
                endif
            enddo
            thermalBLThickness = 0.5d0*((dble(lowerPeakIndex)-0.5d0) + &
                (dble(ny-upperPeakIndex)+0.5d0))/lengthUnit
            thermalBLGridPoints = max(1,ceiling(thermalBLThickness*lengthUnit))

            write(21,'(A)') "# Zhang et al. JFM 814 (2017), Table 1 compatible columns"
            write(21,'(A)') "# Pr Ra Nx Nz Nu Re eps_u_over_exact eps_T_over_exact NBL dg_over_eta dg_over_etaB dt_over_tau_eta"
            write(21,'(2(ES24.16E3,1X),2(I8,1X),4(ES24.16E3,1X),I8,3(1X,ES24.16E3))') &
                Prandtl, Rayleigh, nx+1, ny+1, meanGlobal, meanRe, kineticExactRatioMean, &
                thermalExactRatioMean, thermalBLGridPoints, gridOverKolmogorov, &
                gridOverBatchelor, timeStepOverKolmogorovTime
            write(21,'(A,1X,ES24.16E3)') "# delta_theta_rms_over_H", thermalBLThickness
            write(21,'(A,2(1X,ES24.16E3))') "# eps_u_direct_exact", meanEpsKinetic, exactEpsKinetic
            write(21,'(A,2(1X,ES24.16E3))') "# eps_T_direct_exact", meanEpsThermal, exactEpsThermal
            write(21,'(A,1X,ES24.16E3)') "# eta_over_H", kolmogorovScaleOverH
            write(21,'(A,1X,ES24.16E3)') "# etaB_over_H", batchelorScaleOverH
            write(21,'(A,1X,ES24.16E3)') "# maximum_lattice_CFL", maxStatisticCFL
#endif
        else
            write(21,'(A)') "# no samples: the run ended before the configured averaging window/sample time"
        endif
        close(21)

        open(unit=00, file=trim(settingsFile), status='old', position='append', action='write')
        write(00,*) "Statistic samples =", statisticSamples
        if(statisticSamples.GT.0) then
            write(00,*) "Time-averaged Nu_global =", real(meanGlobal,kind=8), "std =", real(stdGlobal,kind=8)
            write(00,*) "Time-averaged Nu_hot =", real(meanHot,kind=8), "std =", real(stdHot,kind=8)
            write(00,*) "Time-averaged Nu_cold =", real(meanCold,kind=8), "std =", real(stdCold,kind=8)
            write(00,*) "Time-averaged Nu_middle =", real(meanMiddle,kind=8), "std =", real(stdMiddle,kind=8)
            write(00,*) "Table-1 Re=sqrt(<u2+v2>V,t)*H/nu =", real(meanRe,kind=8)
            write(00,*) "Time-averaged <u2+v2>V =", real(meanSpeedSquared,kind=8), &
                "std =", real(stdSpeedSquared,kind=8)
            if((firstHalfSamples.GT.0).AND.(secondHalfSamples.GT.0)) then
                write(00,*) "First/second-half Nu =", real(firstHalfNu,kind=8), real(secondHalfNu,kind=8), &
                    "relative difference =", real(halfWindowNuRelativeDifference,kind=8)
                write(00,*) "First/second-half Re =", real(firstHalfRe,kind=8), real(secondHalfRe,kind=8), &
                    "relative difference =", real(halfWindowReRelativeDifference,kind=8)
            endif
#ifdef RayleighBenardCell
            write(00,*) "Time-averaged Nu_eps_u =", real(meanNuKineticDiss,kind=8), &
                "std =", real(stdNuKineticDiss,kind=8)
            write(00,*) "Time-averaged Nu_eps_T =", real(meanNuThermalDiss,kind=8), &
                "std =", real(stdNuThermalDiss,kind=8)
            write(00,*) "Table-1 eps_u/exact =", real(kineticExactRatioMean,kind=8), &
                "eps_T/exact =", real(thermalExactRatioMean,kind=8)
            write(00,*) "Mean/direct and exact eps_u =", real(meanEpsKinetic,kind=8), &
                real(exactEpsKinetic,kind=8)
            write(00,*) "Mean/direct and exact eps_T =", real(meanEpsThermal,kind=8), &
                real(exactEpsThermal,kind=8)
            write(00,*) "Table-1 NBL/dg_eta/dg_etaB/dt_tauEta =", thermalBLGridPoints, &
                real(gridOverKolmogorov,kind=8), real(gridOverBatchelor,kind=8), &
                real(timeStepOverKolmogorovTime,kind=8)
            write(00,*) "Maximum lattice CFL in statistics window =", real(maxStatisticCFL,kind=8)
#endif
        endif
        close(00)
    end subroutine write_statistics
!===================================================================================================


!===================================================================================================
! 子程序: output_SnapshotFile
! 作用: 保存统计窗口内的无量纲瞬时 u/v、温度和密度，供耗散 PDF、剖面和流场后处理使用。
! 说明: 仅保存宏观场，不包含 f/g，因此不能替代严格的重启文件。
!===================================================================================================
#ifdef unsteadyFlow
    subroutine output_SnapshotFile()
        use commondata
        implicit none
        integer(kind=4) :: i, j
        character(len=12) :: fileNumber
        character(len=100) :: outputFile

        snapshotFileNum = snapshotFileNum+1
        write(fileNumber,'(I12.12)') snapshotFileNum
        outputFile = trim(snapshotFilePrefix)//"-"//fileNumber//".bin"

        open(unit=43, file=trim(outputFile), status='replace', action='write', &
            form='unformatted', access='sequential')
        write(43) ((real(u(i,j)/velocityUnit,kind=8),i=1,nx),j=1,ny)
        write(43) ((real(v(i,j)/velocityUnit,kind=8),i=1,nx),j=1,ny)
        write(43) ((real(T(i,j),kind=8),i=1,nx),j=1,ny)
        write(43) ((real(rho(i,j),kind=8),i=1,nx),j=1,ny)
        close(43)

        open(unit=44, file=trim(snapshotIndexFile), status='old', position='append', action='write')
        write(44,'(I8,1X,I12,1X,ES24.16E3,1X,A)') snapshotFileNum, itc, &
            dble(itc)/timeUnit, trim(outputFile)
        close(44)
        lastSnapshotItc = itc
    end subroutine output_SnapshotFile
!===================================================================================================


!===================================================================================================
! 子程序: output_periodic_Tecplot
! 作用: 为非稳态统计窗口生成带编号的 Tecplot 瞬时场，避免覆盖先前时刻。
!===================================================================================================
    subroutine output_periodic_Tecplot()
        use commondata
        implicit none
        character(len=12) :: fileNumber
        character(len=100) :: outputFile

        tecplotFileNum = tecplotFileNum+1
        write(fileNumber,'(I12.12)') tecplotFileNum
        outputFile = trim(tecplotFilePrefix)//"-"//fileNumber//".plt"
        call output_Tecplot(outputFile)
        lastTecplotItc = itc
    end subroutine output_periodic_Tecplot
#endif
!===================================================================================================


!===================================================================================================
! 子程序: output_Tecplot
! 作用: 输出当前宏观场和由论文式 (35) 得到的温度梯度，便于检查边界层和局部 Nu。
! 说明: 调用前必须先执行 update_host_monitor_2d_openacc()，文件 I/O 始终在主机端进行。
!===================================================================================================
    subroutine output_Tecplot(outputFile)
        use commondata
        implicit none
        integer(kind=4) :: i, j
        character(len=*), intent(in) :: outputFile

        open(unit=41, file=trim(outputFile), status='replace', action='write', form='formatted')
#ifdef RayleighBenardCell
        write(41,'(A)') 'TITLE = "2D OpenACC LBM-CDE Rayleigh-Benard cell"'
#endif
#ifdef SideHeatedCell
        write(41,'(A)') 'TITLE = "2D OpenACC LBM-CDE side-heated cavity"'
#endif
        write(41,'(A)') 'VARIABLES = "X" "Y" "U" "V" "T" "RHO" "GradTx" "GradTy"'
        write(41,'(A,I0,A,I0,A)') 'ZONE I=', nx, ', J=', ny, ', F=POINT'
        do j = 1, ny
            do i = 1, nx
                write(41,'(8(ES24.16E3,1X))') xp(i), yp(j), u(i,j), v(i,j), T(i,j), &
                    rho(i,j), gradTx(i,j), gradTy(i,j)
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
        write(*,'(A,ES16.8,A,ES16.8)') "Nu global/middle =", Nu_global, " / ", Nu_middle
        write(*,'(A,ES16.8,A,ES16.8,A,ES16.8)') "Nu halfway hot/cold =", Nu_hot_halfway, &
            " / ", Nu_cold_halfway, " middle FD =", Nu_middleFD
        write(*,'(A,ES16.8,A,ES16.8)') "Nu eps_u/eps_T =", NuKineticDiss, " / ", NuThermalDiss
        write(*,'(A,ES16.8,A,ES16.8)') "umax/vmax mid =", umax_mid, " / ", vmax_mid
#ifdef unsteadyFlow
        write(*,'(A,I8)') "Statistic samples =", statisticSamples
        write(*,'(A,I8,A,I8)') "Snapshot/Tecplot files =", snapshotFileNum, " / ", tecplotFileNum
#endif

        open(unit=00, file=trim(settingsFile), status='old', position='append', action='write')
        write(00,*) "Final itc =", itc
        write(00,*) "Final errorU =", real(errorU,kind=8), "errorT =", real(errorT,kind=8)
        write(00,*) "Nu_hot =", real(Nu_hot,kind=8), "Nu_cold =", real(Nu_cold,kind=8)
        write(00,*) "Nu_middle =", real(Nu_middle,kind=8), "Nu_global =", real(Nu_global,kind=8)
        write(00,*) "Nu_hot_halfway =", real(Nu_hot_halfway,kind=8), &
            "Nu_cold_halfway =", real(Nu_cold_halfway,kind=8), "Nu_middle_fd =", real(Nu_middleFD,kind=8)
        write(00,*) "ReVolAvg =", real(ReVolAvg,kind=8)
        write(00,*) "epsKineticVol =", real(epsKineticVol,kind=8), &
            "epsThermalVol =", real(epsThermalVol,kind=8)
        write(00,*) "Nu_eps_u =", real(NuKineticDiss,kind=8), &
            "Nu_eps_T =", real(NuThermalDiss,kind=8)
        write(00,*) "eps_u/exact =", real(kineticDissExactRatio,kind=8), &
            "eps_T/exact =", real(thermalDissExactRatio,kind=8)
        write(00,*) "maxMach/minT/maxT/maxAbsRhoDeviation =", real(maxMachLocal,kind=8), &
            real(minTemperature,kind=8), real(maxTemperature,kind=8), real(maxDensityDeviation,kind=8)
        write(00,*) "umax_mid/u0 =", real(umax_mid,kind=8), "position along line =", real(pos_umax_mid,kind=8)
        write(00,*) "vmax_mid/u0 =", real(vmax_mid,kind=8), "position along line =", real(pos_vmax_mid,kind=8)
        write(00,*) "Numax_mid(conductive) =", real(Numax_mid,kind=8), &
            "position along line =", real(pos_Numax_mid,kind=8)
#ifdef unsteadyFlow
        write(00,*) "Snapshot files written =", snapshotFileNum, "last snapshot itc =", lastSnapshotItc
        write(00,*) "Periodic Tecplot files written =", tecplotFileNum, "last Tecplot itc =", lastTecplotItc
#endif
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
        deallocate(sumTemperatureProfile, sumTemperatureSquaredProfile)
    end subroutine finalize_arrays
!===================================================================================================
