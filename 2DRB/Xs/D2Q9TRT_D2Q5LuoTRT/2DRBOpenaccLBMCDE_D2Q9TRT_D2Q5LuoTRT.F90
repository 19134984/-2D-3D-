!=============================================================
!!!    注释区，代码描述
!!!    二维浮力驱动自然对流 OpenACC 并行版本
!!!    流场：D2Q9-TRT LBM-CDE；温度场：Luo D2Q5-TRT
!=============================================================


!=============================================================
!   自定义宏，一些选项的开关
!#define steadyFlow
#define unsteadyFlow

!   流动模式宏的选择，两个都开、两个都关都会报错；只有二选一才通过。
#if defined(steadyFlow) && defined(unsteadyFlow)
#error "Choose only one flow mode: steadyFlow or unsteadyFlow"
#endif
#if !defined(steadyFlow) && !defined(unsteadyFlow)
#error "Define one flow mode: steadyFlow or unsteadyFlow"
#endif

!   流场碰撞策略（四选一）：三种 TRT 奇模态策略，或速度场 BGK 对照。
#define FLOW_ODD_ORIGINAL_MAGIC
!#define FLOW_ODD_EFFECTIVE_MAGIC
!#define FLOW_ODD_FIXED_SQ
!#define FLOW_BGK

#if defined(FLOW_ODD_ORIGINAL_MAGIC) && defined(FLOW_ODD_EFFECTIVE_MAGIC)
#error "Choose only one flow odd magic policy"
#endif
#if defined(FLOW_ODD_ORIGINAL_MAGIC) && defined(FLOW_ODD_FIXED_SQ)
#error "Choose only one flow odd magic policy"
#endif
#if defined(FLOW_ODD_EFFECTIVE_MAGIC) && defined(FLOW_ODD_FIXED_SQ)
#error "Choose only one flow odd magic policy"
#endif
#if defined(FLOW_ODD_ORIGINAL_MAGIC) && defined(FLOW_BGK)
#error "Choose only one flow collision policy"
#endif
#if defined(FLOW_ODD_EFFECTIVE_MAGIC) && defined(FLOW_BGK)
#error "Choose only one flow collision policy"
#endif
#if defined(FLOW_ODD_FIXED_SQ) && defined(FLOW_BGK)
#error "Choose only one flow collision policy"
#endif
#if !defined(FLOW_ODD_ORIGINAL_MAGIC) && !defined(FLOW_ODD_EFFECTIVE_MAGIC) && !defined(FLOW_ODD_FIXED_SQ) && !defined(FLOW_BGK)
#error "Define one flow collision policy"
#endif

!   速度边界，包括水平垂直边界无滑移，还有垂直边界速度周期
#define HorizontalWallsNoslip
#define VerticalWallsNoslip
!#define VerticalWallsPeriodicalU

!   垂直速度边界宏的选择
#if defined(VerticalWallsNoslip) && defined(VerticalWallsPeriodicalU)
#error "Choose only one vertical velocity BC: VerticalWallsNoslip or VerticalWallsPeriodicalU"
#endif
#if !defined(VerticalWallsNoslip) && !defined(VerticalWallsPeriodicalU)
#error "Define one vertical velocity BC: VerticalWallsNoslip or VerticalWallsPeriodicalU"
#endif

!   温度边界(for Rayleigh Benard Cell)，包括水平边界恒温，垂直边界温度不可穿透以及周期
#define RayleighBenardCell
#define HorizontalWallsConstT
#define VerticalWallsAdiabatic
!#define VerticalWallsPeriodicalT



!   温度边界(for Side Heated Cell)，包括水平边界温度不可穿透，垂直边界恒温,侧壁加热加磁场
!#define SideHeatedCell
!#define HorizontalWallsAdiabatic
!#define VerticalWallsConstT
!#define SideHeatedHa
!~~temperature B.C.~~

!   对流算例宏的选择
#if defined(RayleighBenardCell) && defined(SideHeatedCell)
#error "Choose only one convection case: RayleighBenardCell or SideHeatedCell"
#endif
#if !defined(RayleighBenardCell) && !defined(SideHeatedCell)
#error "Define one convection case: RayleighBenardCell or SideHeatedCell"
#endif


!   自定义宏结束
!=============================================================


!=============================================================
!   全局模块
    module commondata
        implicit none
        !===============================================================================================
        ! 是否在计算前从旧算例重启
        integer(kind=4), parameter :: loadInitField=0   ! 0: 不重启；1: 按 latest .meta 自动续算

        ! 正常断电续算只需要设置 loadInitField=1；
        ! 代码只读取 <reloadFilePrefix>-latest.meta，并从里面找到最新的 .bin。
        ! 正常续算不用改 reloadFileNum；只有 latest .meta 缺失时，
        ! 才手动设置 reloadFileNum 作为保守推断编号。重启文件的编号
        integer(kind=4) :: reloadFileNum=0              ! latest .meta 存在时会被覆盖；meta 缺失时作为手工兜底编号
        !===============================================================================================
        real(kind=8) :: restartTfOffset=0.0d0            ! 续算前已累计的 t/t_ff；优先从 latest .meta 读取，meta 缺失时由代码推断
        integer(kind=4) :: restartItcOffset=0           ! 续算前已累计的格子步数；优先从 latest .meta 读取，meta 缺失时由代码推断
        logical :: reloadMetadataLoaded=.false.         ! 标记是否成功读取 reload 元数据文件
        !===============================================================================================

        !===============================================================================================
        ! 无量纲参数
        integer(kind=4), parameter :: nx=2048, ny=2048     !流体节点数
        real(kind=8), parameter :: rho0=1.0d0              !弱可压缩模型参考密度
        real(kind=8), parameter :: cs2=1.0d0/3.0d0         !D2Q9 格子声速平方
#ifdef SideHeatedCell
        real(kind=8), parameter :: lengthUnit=dble(nx)     !侧壁差温：特征长度取 x 方向长度
#else
        real(kind=8), parameter :: lengthUnit=dble(ny)     !上下差温：特征长度取 y 方向长度
#endif
        real(kind=8), parameter :: pi = acos(-1.0d0)

        real(kind=8), parameter :: Rayleigh=1.0d10
        real(kind=8), parameter :: Prandtl=0.7d0
        real(kind=8), parameter :: Mach=0.1d0
        real(kind=8), parameter :: Thot=0.5d0, Tcold=-0.5d0
        real(kind=8), parameter :: Tref=0.5d0*(Thot+Tcold)
        real(kind=8), parameter :: deltaT=Thot-Tcold
        real(kind=8), parameter :: viscosity=Mach*lengthUnit*dsqrt(Prandtl/(3.0d0*Rayleigh))
        real(kind=8), parameter :: diffusivity=viscosity/Prandtl !目标热扩散率 kappa
        real(kind=8), parameter :: chi_nu=0.5d0          !手动设置的流场剪切修正参数
        real(kind=8), parameter :: chi_b=0.0d0           !手动设置的流场体黏度修正参数
        real(kind=8), parameter :: tauf=0.5d0+viscosity/(cs2*(1.0d0-chi_nu)) !基础 tau_0
        real(kind=8), parameter :: bulkViscosity=(tauf-0.5d0)*(1.0d0-chi_b)*cs2


        ! velocityScaleCompare is used only in velocity-related post-processing to convert lattice velocity
        ! to the nondimensional velocity scale adopted by the reference paper being compared.
        ! 默认采用热扩散标度 UL/kappa；若按自由落体标度比较，可改为 1.0d0/velocityUnit。
        real(kind=8), parameter :: velocityScaleCompare=lengthUnit/diffusivity

        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1


#ifdef  SideHeatedHa
        real(kind=8), parameter :: Ha=20.0d0                           !磁场强度
        real(kind=8), parameter :: phi=(0.0d0)*(pi/180.0d0)            !磁场角度，以水平向右为0，修改0.0d0即可
        real(kind=8), parameter :: B2sigemarho=(Ha**2*viscosity)/(lengthUnit*lengthUnit)  !动量方程上的源项系数
#endif

        !===============================================================================================
        ! Luo D2Q5-TRT 温度模型参数
        real(kind=8), parameter :: paraA=20.0d0*dsqrt(3.0d0)*diffusivity-4.0d0



        ! 浮力项参数
        real(kind=8), parameter :: gBeta=Rayleigh*viscosity*diffusivity/(deltaT*lengthUnit**3) !g*beta

        real(kind=8), parameter :: timeUnit=dsqrt(lengthUnit/(gBeta*deltaT)) !一个自由落体时间对应的格子步数
        real(kind=8), parameter :: velocityUnit=dsqrt(gBeta*deltaT*lengthUnit) !自由落体速度

        real(kind=8), parameter :: Snu=1.0d0/tauf
        real(kind=8), parameter :: flowMagicParameter=3.0d0/16.0d0
#ifdef FLOW_ODD_ORIGINAL_MAGIC
        real(kind=8), parameter :: Sq=1.0d0/(0.5d0+flowMagicParameter/(tauf-0.5d0))
#endif
#ifdef FLOW_ODD_EFFECTIVE_MAGIC
        real(kind=8), parameter :: Sq=1.0d0/(0.5d0+flowMagicParameter/(viscosity/cs2))
#endif
#ifdef FLOW_ODD_FIXED_SQ
        real(kind=8), parameter :: Sq=1.0d0 ! 固定奇模态松弛率对照：tau_q=1
#endif
#ifdef FLOW_BGK
        real(kind=8), parameter :: Sq=Snu ! 速度场 BGK：所有非守恒矩均使用 tau_0
#endif


        real(kind=8), parameter :: Qk=3.0d0-dsqrt(3.0d0)

        real(kind=8), parameter :: Qe=4.0d0*dsqrt(3.0d0)-6.0d0
        real(kind=8), parameter :: Qnu=4.0d0*dsqrt(3.0d0)-6.0d0
        ! 初始化平衡态 geq_i=omegaT_i*T*(1+thermalGeqCoeff*e_i.u) 中的速度系数；
        ! 该取值保证 D2Q5 平衡态的一阶矩严格恢复为 T*u 和 T*v。
        real(kind=8), parameter :: thermalGeqCoeff=10.0d0/(4.0d0+paraA)
        !===============================================================================================

        !===============================================================================================
        ! 输出/备份相关设置（以自由落体时间 t_ff 为单位）
        real(kind=8), parameter :: epsU=1.0d-7, epsT=1.0d-7    ! 稳态收敛阈值

#ifdef steadyFlow
        real(kind=8), parameter :: outputSnapshotInterval=10.0d0   ! 快照和 Nu/Re 时间序列采样间隔（单位：t_ff）
        real(kind=8), parameter :: reloadFileInterval=100.0d0  ! f/g 重启文件输出间隔（单位：t_ff）
        real(kind=8), parameter :: outputPltFileInterval=100.0d0  ! Tecplot 文件输出间隔（单位：t_ff）
        integer(kind=4), parameter :: statisticSampleCountMax= &
            int(12000.0d0/outputSnapshotInterval)
        integer(kind=4), parameter :: outputSnapshotFile=1   ! 是否输出后处理快照文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputPltFile=1   ! 是否输出 plt 文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputReloadFile=1 ! 是否周期输出 f/g 重启文件：0=不输出，1=输出
        ! 稳态最多推进 12000 t_ff；
        ! 实际计算仍可由 errorU/errorT 达到阈值而提前停止。
        integer(kind=4), parameter :: itc_max=max(1,int(12000.0d0*timeUnit+0.5d0))
#endif

#ifdef unsteadyFlow
        real(kind=8), parameter :: statisticSampleInterval=1.0d0 ! Nu/Re、耗散和温度剖面统计间隔（单位：t_ff）
        real(kind=8), parameter :: outputSnapshotInterval=20.0d0 ! uvTrho 全场快照输出间隔（单位：t_ff）
        real(kind=8), parameter :: reloadFileInterval=100.0d0  ! f/g 重启文件输出间隔（单位：t_ff）
        real(kind=8), parameter :: outputPltFileInterval=100.0d0  ! Tecplot 文件周期输出间隔（单位：t_ff）
        real(kind=8), parameter :: unsteadyRunDuration=1000.0d0  ! 非稳态总目标时长，续算时只补足到该 t_ff
        ! 以下三个参数控制非稳态统计平均窗口，不改变推进时长或采样频率。
        ! 时间以整个算例的绝对 t_ff 计；续算会重读 Nu/Re、耗散和温度剖面历史，恢复完整窗口累计量。
        ! 默认统计窗口取总时长后 1/2；Nu/Re 最终结果取该窗口后半段，即默认对应总时长最后 1/4。
        real(kind=8), parameter :: unsteadyAverageStartTf=0.5d0*unsteadyRunDuration  ! 平均窗口起点
        real(kind=8), parameter :: unsteadyAverageEndTf=unsteadyRunDuration          ! 平均窗口终点
        real(kind=8), parameter :: unsteadyAverageMidTf=0.5d0*(unsteadyAverageStartTf+unsteadyAverageEndTf) ! 前/后半分界
        ! Nu/Re、耗散、可压缩性和温度剖面只在统计窗口内采样；快照/plt/reload 仍从计算开始输出。
        ! 本次进程 Nu/Re 数组的容量；续算后只保存本次进程新产生的样本。
        integer(kind=4), parameter :: statisticSampleCountMax=max(1, &
            nint((unsteadyAverageEndTf-unsteadyAverageStartTf)/statisticSampleInterval)+1)
        integer(kind=4), parameter :: outputSnapshotFile=1   ! 是否独立输出 uvTrho 快照：0=不输出，1=输出
        integer(kind=4), parameter :: outputPltFile=1   ! 是否输出 plt 文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputReloadFile=1 ! 是否周期输出 f/g 重启文件：0=不输出，1=输出
        integer(kind=4), parameter :: itc_max=max(1, int(unsteadyRunDuration*timeUnit+0.5d0)) ! 非稳态：由总 t_ff 自动换算格子步
#endif

        integer(kind=4) :: snapshotFileNum, pltFileNum  ! 快照/plt 输出文件的计数器
        ! 每次调用对应输出子程序时递增（用于文件名编号）

        integer(kind=4) :: currentRunStatisticSampleCount ! 本次进程新计算的样本数
        integer(kind=4) :: outputSnapshotIntervalItc
        integer(kind=4) :: reloadFileIntervalItc, outputPltFileIntervalItc


        real(kind=8) :: NuVolAvg(0:statisticSampleCountMax)
        real(kind=8) :: ReVolAvg(0:statisticSampleCountMax)
        ! NuVolAvg：当前采样时刻由体平均对流热通量得到的 Nu。
        ! ReVolAvg：当前采样时刻的空间 RMS Re。

        character(len=100) :: snapshotFilePrefix="buoyancyCavity2DOpenaccLBMCDE_D2Q5Snapshot"
        ! 快照输出文件前缀（实际文件名形如：<snapshotFilePrefix>-<编号>.bin）

        character(len=100) :: pltFolderPrefix="buoyancyCavity2DOpenaccLBMCDE_D2Q5Tecplot"
        ! plt 输出文件前缀（实际文件名形如：<pltFolderPrefix>-<编号>.plt）

        character(len=100) :: reloadFilePrefix="reloadFile2DOpenaccLBMCDE_D2Q5"
        ! 重启读取文件的前缀；latest meta 模式实际读取 meta 中记录的 <reloadFilePrefix>-<编号>.bin

        character(len=100) :: settingsFile="SimulationSettings2DOpenaccLBMCDE_D2Q5.txt"
        !===============================================================================================

        !===============================================================================================
        !计算中需要的相关参数
        real(kind=8) :: errorU, errorT

        real(kind=8) :: xp(0:nx+1), yp(0:ny+1)      !无量纲的坐标数组，包括边界
        real(kind=8), allocatable :: u(:,:), v(:,:), T(:,:), rho(:,:)

#ifdef steadyFlow
        real(kind=8), allocatable :: up(:,:), vp(:,:), Tp(:,:)   !存储之前的数据，用来算收敛判据
#endif
        real(kind=8), allocatable :: f(:,:,:), f_post(:,:,:)
        real(kind=8), allocatable :: g(:,:,:), g_post(:,:,:)
        real(kind=8), allocatable :: Fx(:,:), Fy(:,:)
        real(kind=8), allocatable :: Sxx(:,:), Sxy(:,:), Syy(:,:), Sdiv(:,:) !局部应变率

        integer(kind=4) :: itc
#ifdef steadyFlow
        real(kind=8) :: Nu_global, Nu_hot, Nu_cold, Nu_middle    !平均Nu，全场，侧壁以及中线
        real(kind=8) :: Nu_hot_max, Nu_hot_min, Nu_hot_max_position, Nu_hot_min_position    !左侧壁面的最大最小Nu，以及对应的位置
#endif


        !格子离散速度和权重
        integer(kind=4) :: ex(0:8), ey(0:8)
        data ex/0, 1, 0, -1,  0, 1, -1, -1,  1/
        data ey/0, 0, 1,  0, -1, 1,  1, -1, -1/
        real(kind=8) :: omega(0:8)  !D2Q9 流场权重
        real(kind=8) :: omegaT(0:4) !D2Q5 温度权重

#ifdef unsteadyFlow
        ! 命名约定：xxxVol 是当前时刻的空间求和；xxxVolAvg 是除以 nx*ny 后的瞬时体平均。
        !             sumXxxVolAvg 是瞬时体平均量的时间求和；xxxVolTimeAvg 是最终空间-时间平均。
        real(kind=8) :: speedSquaredVolAvg        ! <u^2+v^2>_V；最终 Re=sqrt(<u^2+v^2>_V,t)*H/nu
        real(kind=8) :: epsKineticVolAvg          ! <epsilon_u>_V；由非平衡应变率计算的动能耗散
        real(kind=8) :: epsThermalVolAvg          ! <epsilon_T>_V；由温度梯度计算的热耗散
        real(kind=8) :: densityFluctuationSquaredVolAvg ! <(rho-rho0)^2>_V；开平方后才是瞬时 RMS 密度
        real(kind=8) :: velocityDivergenceSquaredVolAvg ! <(div u)^2>_V；开平方后才是瞬时 RMS 散度
        real(kind=8) :: maxMachLocal, minTemperature, maxTemperature ! 数值稳定性辅助诊断

        ! 统计窗口累计量：用于最终空间-时间平均、耗散精确关系和可压缩性检查。
        real(kind=8) :: sumNuVolAvg                  ! sum Nu_vol_avg(t)，瞬时体平均 Nu 的时间求和
        real(kind=8) :: sumSpeedSquaredVolAvg        ! sum <u^2+v^2>_V，用于文献定义的空间-时间 RMS Re
        real(kind=8) :: sumEpsKineticVolAvg, sumEpsThermalVolAvg ! 两类瞬时体平均耗散的时间求和
        real(kind=8) :: sumDensityFluctuationSquaredVolAvg ! sum <(rho-rho0)^2>_V，用于空间-时间 RMS
        real(kind=8) :: sumVelocityDivergenceSquaredVolAvg ! sum <(div u)^2>_V，用于空间-时间 RMS
        real(kind=8) :: maxStatisticCFL              ! 统计窗口采样点中的最大格子 CFL=max|u|
        ! 以下两个剖面累计量保留温度 RMS 边界层诊断，用于计算热边界层几个网格。N_BL 暂时先由 H/(2*Nu) 估计。
        real(kind=8), allocatable :: sumTemperatureXAvgProfile(:), sumTemperatureSquaredXAvgProfile(:)
        integer(kind=4) :: cumulativeStatisticSampleCount ! 完整统计窗口累计样本数（历史恢复+本次新增）
        character(len=100) :: dissipationHistoryFile="DissipationHistory_2DOpenaccLBMCDE.dat" ! 瞬时原始量历史
        character(len=100) :: temperatureProfileHistoryFile= &
            "TemperatureProfileHistory_2DOpenaccLBMCDE.bin" ! 每个统计采样的水平平均 T(y) 与 T^2(y)，体平均，供续算恢复
        character(len=100) :: temperatureRmsProfileFile= &
            "TemperatureRmsProfile_2DOpenaccLBMCDE_D2Q5.dat" ! 最终 z/H 与 theta_rms(z) 剖面，时间平均，根据T_rms峰值确定边界层厚度
        character(len=100) :: statisticsFile="NuReDissStatistics_2DOpenaccLBMCDE.dat" ! 最终统计结果
#endif
        !===============================================================================================

    end module commondata

!   全局模块结束
!=============================================================


!=============================================================

    program main

    use openacc
    use commondata
    implicit none
    real(kind=8) :: timeStart, timeEnd
    real(kind=8) :: timeStart2, timeEnd2
    character(len=24) :: ctime, string
    INTEGER(kind=4) :: time
    integer(kind=4) :: numAccDevices
    integer(kind=8) :: wallClockStart, wallClockEnd, wallClockRate


    !===============================================================================================
    ! 初始化 OpenACC 设备
    if(loadInitField.EQ.1) then
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,*) " "
        write(00,*) "================ Restart continuation begins ================"
    else
        open(unit=00,file=trim(settingsFile),status='replace')
    endif
    string = ctime( time() )                      !ctime把 time() 返回的时间戳转换成可读的字符串
    write(00,*) 'Start: ', string                 !什么时候开始计算
    write(00,*) "Starting OpenACC >>>>>>"
    call acc_init(acc_device_default)
    numAccDevices = acc_get_num_devices(acc_device_default)
    write(00,*) "Visible OpenACC devices:", numAccDevices
    close(00)
    !===============================================================================================


    !===============================================================================================
    ! Initialization
    call initial()
    call enter_data_2d_openacc()

    !===============================================================================================

    call CPU_TIME(timeStart)         !当前进程累计消耗的 CPU 时间,包括并行
    ! system_clock 返回墙钟计数器和每秒计数率；
    ! 下面用 counter/rate 把它换算成实际经过的秒数。
    call system_clock(wallClockStart, wallClockRate)
    timeStart2 = dble(wallClockStart) / dble(max(wallClockRate,1_8))

    !开始计算
#ifdef steadyFlow
    do while( ((errorU.GT.epsU).OR.(errorT.GT.epsT)).AND.(itc.LE.itc_max) )
#endif
#ifdef unsteadyFlow
    ! 非稳态按累计格子步推进；续算时自然从 restartItcOffset 继续到 itc_max。
    do while( (restartItcOffset+itc).LT.itc_max )
#endif

        itc = itc+1

        call collision()

        call streaming()

        call bounceback()

        call macro()

        call collisionT()

        call streamingT()

        call bouncebackT()

        call macroT()

#ifdef steadyFlow
        ! 输出按累计格子步判断；否则从 1050tf 续算会在 1150tf 才输出，
        ! 而不是接回不断电运行应有的 1100tf、1200tf、...
        if(MOD(restartItcOffset+itc,2000).EQ.0) call check()
        if( (outputPltFile.EQ.1).AND.(MOD(restartItcOffset+itc, outputPltFileIntervalItc).EQ.0) ) then
            call update_host_tecplot_2d_openacc()
            call output_Tecplot()  !稳态模式下的可选周期 Tecplot 输出
        endif
        if( (outputSnapshotFile.EQ.1).AND.(MOD(restartItcOffset+itc, outputSnapshotIntervalItc).EQ.0) ) then
            call update_host_snapshot_2d_openacc()
            call output_SnapshotFile()  !稳态模式下的可选周期 uvTrho 快照输出
        endif
        if( (outputReloadFile.EQ.1).AND.(MOD(restartItcOffset+itc, reloadFileIntervalItc).EQ.0) ) then
            call update_host_reload_2d_openacc()
            call output_ReloadFile()      !稳态模式下的可选周期 f/g 重启文件输出
        endif
#endif

#ifdef unsteadyFlow
        ! 统计采样，统计窗口内的 Nu/Re、耗散和温度剖面。
        ! 每个目标时刻都从绝对 t/t_ff 换算为累计格子步，当前完成样本个数为 cumulativeStatisticSampleCount。
        if((unsteadyAverageStartTf+real(cumulativeStatisticSampleCount,kind=8)* &
            statisticSampleInterval.LE.unsteadyAverageEndTf).AND. &
           ((restartItcOffset+itc).GE.max(1,nint((unsteadyAverageStartTf+ &
            real(cumulativeStatisticSampleCount,kind=8)*statisticSampleInterval)*timeUnit)))) then
            ! 一次全场遍历同时计算 Nu、Re、两类耗散和可压缩性指标，避免重复计算速度平方。
            call calculate_unsteady_sample()
            !$acc update self(T)
            call accumulate_dissipation_statistics() !所有样本都位于统计窗口内
        endif
        ! 快照按 outputSnapshotInterval 输出。
        if( (outputSnapshotFile.EQ.1).AND. &
            (MOD(restartItcOffset+itc,outputSnapshotIntervalItc).EQ.0) ) then
            call update_host_snapshot_2d_openacc()
            call output_SnapshotFile()
        endif
        if( (outputPltFile.EQ.1).AND.(MOD(restartItcOffset+itc, outputPltFileIntervalItc).EQ.0) ) then
            call update_host_tecplot_2d_openacc()
            call output_Tecplot()  !非稳态模式下的可选周期 Tecplot 输出
        endif
        if( (outputReloadFile.EQ.1).AND.(MOD(restartItcOffset+itc, reloadFileIntervalItc).EQ.0) ) then
            call update_host_reload_2d_openacc()
            call output_ReloadFile()      !非稳态模式下的可选周期 f/g 重启文件输出
        endif
#endif
     enddo

    !$acc wait(1)
    call CPU_TIME(timeEnd)         !当前进程累计消耗的 CPU 时间,包括并行
    ! 取墙钟结束计数器，用于后面输出 OpenACC 实际耗时。
    call system_clock(wallClockEnd, wallClockRate)
    timeEnd2 = dble(wallClockEnd) / dble(max(wallClockRate,1_8))
    call update_host_snapshot_2d_openacc()

#ifdef steadyFlow
    call output_Tecplot()          !输出最后一步的plt结果
    call output_SnapshotFile()              !输出最后一步的uvTrho数据
#endif

#ifdef unsteadyFlow
    call output_unsteady_NuRe_postprocess()
    call write_dissipation_statistics()
#endif

    !===============================================================================================



    !===============================================================================================


#ifdef steadyFlow
! 稳态最终标量诊断：
! 1) 只在 steadyFlow 收敛后调用
! 2) SideHeatedCell 的主热流方向为 x，用 u 和 dT/dx；RayleighBenardCell 的主热流方向为 y，用 v 和 dT/dy。
! 3) 壁面 Nu 和角点扩展默认采用半步长边界：流体节点距离物理边界 dx/2 或 dy/2。
! 4) Nu 极值、中心线速度极值使用五点最小二乘抛物线插值；中心线在偶数网格时用两侧流体节点线性插值。
! 5) 如果后续改成周期速度/温度边界，或改变半步长边界布置，这里和各后处理子程序都需要重新检查。
#ifdef SideHeatedCell
    call SideHeatedcalc_Nu_global()          ! 全场平均Nu
    call SideHeatedcalc_Nu_wall_avg()  ! 热/冷壁, 中线平均Nu,以及热壁最大Numax和Numin以及位置，都采用五点最小二乘法插值出来

    call SideHeatedcalc_umid_max()     !中心线上的最大速度及其位置，也是用五点最小二乘法插值出来
    call SideHeatedcalc_vmid_max()
#endif

#ifdef RayleighBenardCell
    call RBcalc_Nu_global()          ! 全场平均Nu
    call RBcalc_Nu_wall_avg()  ! 热/冷壁, 中线平均Nu,以及热壁最大Numax和Numin以及位置，都采用五点最小二乘法插值出来

    call RBcalc_umid_max()     !中心线上的最大速度及其位置，也是用五点最小二乘法插值出来
    call RBcalc_vmid_max()
#endif
#endif


#ifdef steadyFlow
#ifdef VerticalWallsNoslip
    ! psi/vort 后处理默认封闭腔体：四周无滑移，psi 在物理边界取同一常数。
    ! 若垂直边界改为周期速度边界，流函数边界补点和涡量单边差分需要另写周期版本。
    call calc_psi_vort_and_output()  ! 输出中心abs(psi), max(abs(psi))及位置；max位置用细网格样条插值
#endif

    call calNuRe()
#endif






    open(unit=00,file=trim(settingsFile),status='unknown',position='append')        !在这个txt文件后面继续写（追加模式）
    write(00,*) "======================================================================"
    write(00,*) "Time (CPU) = ", real(timeEnd-timeStart,kind=8), "s"                             !当前进程累计消耗的 CPU 时间,包括并行
    write(00,*) "MLUPS = ", real( dble(nx)*dble(ny)*dble(itc)/(timeEnd-timeStart)/1.0d6,kind=8 )   !百万格点更新/秒
    write(00,*) "Time (ACC) = ", real(timeEnd2-timeStart2,kind=8), "s"                           !墙钟时间
    write(00,*) "MLUPS (ACC) = ", real( dble(nx)*dble(ny)*dble(itc)/(timeEnd2-timeStart2)/1.0d6,kind=8 )   !百万格点更新/秒
#ifdef steadyFlow
    write(00,'(a,1x,ES24.16E3)') "Nu_global =", real(Nu_global,kind=8)
    write(00,'(a,1x,ES24.16E3)') "Nu_hot    =", real(Nu_hot,kind=8)
    write(00,'(a,1x,ES24.16E3)') "Nu_cold   =", real(Nu_cold,kind=8)
#endif
    write(00,*) "Dellocate Array......"
    call exit_data_2d_openacc()
    deallocate(f)
    deallocate(g)
    deallocate(f_post)
    deallocate(g_post)
    deallocate(u)
    deallocate(v)
    deallocate(T)
#ifdef steadyFlow
    deallocate(up)
    deallocate(vp)
    deallocate(Tp)
#endif
    deallocate(rho)
    deallocate(Fx)
    deallocate(Fy)
    deallocate(Sxx)
    deallocate(Sxy)
    deallocate(Syy)
    deallocate(Sdiv)
#ifdef unsteadyFlow
    deallocate(sumTemperatureXAvgProfile)
    deallocate(sumTemperatureSquaredXAvgProfile)
#endif
    write(00,*) "    "

    write(00,*) "Successfully: DNS completed!"

    string = ctime( time() )
    write(00,*) 'End:   ', string           !什么时候算完
    close(00)


    end program main

!===========================================================================================================================


!===================================================================================================
! 子程序: initial
! 作用: 初始化网格坐标、场变量、分布函数、输出文件和重启信息。
!===========================================================================================================================
  subroutine initial()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: un(0:8)
    real(kind=8) :: us2
    character(len=100) :: reloadFileName


    itc = 0
    errorU = 100.0d0
    errorT = 100.0d0
    snapshotFileNum = 0
    pltFileNum = 0
    restartItcOffset = 0
#ifdef unsteadyFlow
    cumulativeStatisticSampleCount = 0
#endif
    reloadMetadataLoaded = .false.
    outputSnapshotIntervalItc = max(1, int(outputSnapshotInterval*timeUnit+0.5d0))
    reloadFileIntervalItc = max(1, int(reloadFileInterval*timeUnit+0.5d0))
    outputPltFileIntervalItc = max(1, int(outputPltFileInterval*timeUnit+0.5d0))


    !-----------------------------------------------------------------------------------------------
    !记录各种信息在日志文件中
    open(unit=00,file=trim(settingsFile),status='unknown',position='append')  !在这个txt文件后面继续写（追加模式）

    if(outputSnapshotFile.EQ.1) then
        open(unit=01,file=trim(snapshotFilePrefix)//"-"//"readme",status="unknown")    !trim去掉字符串尾部空格，换了存储路径，可自己更改
        write(01,*) "snapshot file prefix exists!"
        close(01)
        write(00,*) "Snapshot data will be stored in ", snapshotFilePrefix
    endif
    if(outputPltFile.EQ.1) then
        open(unit=01,file=trim(pltFolderPrefix)//"-"//"readme",status="unknown")     !读取路径pltFolderPrefix="../pltFile/buoyancyCavity
        write(01,*) "pltFile folder exist!"
        close(01)
        write(00,*) "Data will be stored in ", pltFolderPrefix
    endif
    if(outputReloadFile.EQ.1) then
        open(unit=01,file=trim(reloadFilePrefix)//"-"//"readme",status="unknown")
        write(01,*) "reloadFile prefix exists!"
        close(01)
        write(00,*) "Reload data will be stored in ", reloadFilePrefix
    endif

    if( (paraA.GE.1.0d0).OR.(paraA.LE.-4.0d0) ) then                           !只有在[-4,1]才可以，要不然预警退出
        write(00,*) "----------------------------------"
        write(00,*) "paraA=", paraA
        write(00,*) "Error: condition not meet for the legacy thermal algorithm"
        write(00,*) "Ref: Luo2013, CMA"
        write(00,*) "Please try to reduce Mach number"
        write(00,*) "----------------------------------"
        stop
    endif
    if((chi_nu.GE.1.0d0).OR.(chi_b.GT.1.0d0)) then
        write(00,*) "Error: chi_nu must be smaller than one and chi_b must not exceed one"
        stop
    endif
    if((Snu.LE.0.0d0).OR.(Snu.GE.2.0d0).OR.(Sq.LE.0.0d0).OR.(Sq.GE.2.0d0)) then
        write(00,*) "Error: flow TRT relaxation rates must lie strictly between zero and two"
        stop
    endif
    if((Qk.LE.0.0d0).OR.(Qk.GE.2.0d0).OR.(Qe.LE.0.0d0).OR.(Qe.GE.2.0d0).OR. &
        (Qnu.LE.0.0d0).OR.(Qnu.GE.2.0d0)) then
        write(00,*) "Error: Luo D2Q5 relaxation rates must lie strictly between zero and two"
        stop
    endif
#ifdef unsteadyFlow
    ! 统计时刻直接从窗口起点递增；统计时长能被采样间隔整除由用户设置参数时保证。
    if((statisticSampleInterval.LE.0.0d0).OR.(unsteadyAverageStartTf.LT.0.0d0).OR. &
        (unsteadyAverageEndTf.LT.unsteadyAverageStartTf).OR. &
        (unsteadyAverageEndTf.GT.unsteadyRunDuration+1.0d-12)) then
        write(00,*) "Error: invalid unsteady statistics window or sampling interval"
        stop
    endif
#endif

    write(00,*)"-------------------------------------------------------------------------------"
    write(00,*) 'Mesh:',nx,ny
    write(00,*) 'Rayleigh=',real(Rayleigh,kind=8), '; Prandtl =',real(Prandtl,kind=8), '; Mach =',real(Mach,kind=8)
    write(00,*) "Length unit: L0 =", real(lengthUnit,kind=8)
    write(00,*) "Time unit: Sqrt(L0/(gBeta*DeltaT)) =", real(timeUnit,kind=8)
    write(00,*) "Velocity unit: Sqrt(gBeta*L0*DeltaT) =", real(velocityUnit,kind=8)
    write(00,*) "   "
    write(00,*) 'chi_nu=',real(chi_nu,kind=8), '; chi_b=',real(chi_b,kind=8)
    write(00,*) 'tau_0/tauf=',real(tauf,kind=8), '; Snu=',real(Snu,kind=8), '; Sq=',real(Sq,kind=8)
#ifdef FLOW_ODD_ORIGINAL_MAGIC
    write(00,*) "Flow odd magic policy = original tau_0 scale"
#endif
#ifdef FLOW_ODD_EFFECTIVE_MAGIC
    write(00,*) "Flow odd magic policy = chi_nu-corrected effective scale"
#endif
#ifdef FLOW_ODD_FIXED_SQ
    write(00,*) "Flow odd magic policy = fixed Sq = 1 control"
#endif
#ifdef FLOW_BGK
    write(00,*) "Flow collision policy = BGK (Sq = Snu; no flow magic parameter)"
#endif
    write(00,*) "thermalScheme = Luo D2Q5-TRT  time-difference correction"
    write(00,*) 'Qk=',real(Qk,kind=8), '; Qe=',real(Qe,kind=8), &
        '; Qnu=',real(Qnu,kind=8), '; paraA=',real(paraA,kind=8)
    write(00,*) "viscosity =",real(viscosity,kind=8), "; diffusivity =",real(diffusivity,kind=8)
    write(00,*) "outputSnapshotFile =", outputSnapshotFile
    write(00,*) "outputSnapshotInterval =", real(outputSnapshotInterval,kind=8), "free-fall time units"
    write(00,*) "outputSnapshotIntervalItc =", outputSnapshotIntervalItc, "in itc units"
    write(00,*) "outputPltFile =", outputPltFile
    write(00,*) "outputPltFileInterval =", real(outputPltFileInterval,kind=8), "free-fall time units"
    write(00,*) "outputPltFileIntervalItc =", outputPltFileIntervalItc, "in itc units"
    write(00,*) "outputReloadFile =", outputReloadFile
    write(00,*) "reloadFileInterval =", real(reloadFileInterval,kind=8), "free-fall time units"
    write(00,*) "reloadFileIntervalItc =", reloadFileIntervalItc, "in itc units"
#ifdef unsteadyFlow
    write(00,*) "statisticSampleInterval =", real(statisticSampleInterval,kind=8), &
        "free-fall time units"
    write(00,*) "unsteadyRunDuration =", real(unsteadyRunDuration,kind=8), "free-fall time units"
    write(00,*) "maximum samples in Nu/Re arrays =", statisticSampleCountMax
    write(00,*) "start/mid/end statistic sample time_tf =", &
        unsteadyAverageStartTf, unsteadyAverageMidTf, unsteadyAverageEndTf
    write(00,*) "start/mid/end statistic sample itc =", &
        max(1,nint(unsteadyAverageStartTf*timeUnit)), &
        max(1,nint(unsteadyAverageMidTf*timeUnit)), &
        max(1,nint(unsteadyAverageEndTf*timeUnit))
    write(00,*) "Statistics histories cover only the averaging window; time axis uses prescribed t/t_ff"
    write(00,*) "dissipationHistoryFile =", trim(dissipationHistoryFile)
    write(00,*) "temperatureProfileHistoryFile =", trim(temperatureProfileHistoryFile)
#endif
#ifdef steadyFlow
    write(00,*) "Steady convergence/Nu-Re history axis = cumulative lattice itc"
#endif
    if(loadInitField.EQ.1) then
        write(00,*) "Restart offsets will be read from reload metadata when available."
    endif
    write(00,*) "itc_max =",itc_max
    write(00,*) "default epsU =", real(epsU,kind=8),"; epsT =", real(epsT,kind=8)
    write(00,*) "    "

#ifdef RayleighBenardCell
    write(00,*) "I am Rayleigh Benard Cell"
#endif
#ifdef SideHeatedCell
    write(00,*) "I am Side Heated Cell"
#endif

#ifdef steadyFlow
    write(00,*) "I am steadyFlow"
#endif
#ifdef unsteadyFlow
    write(00,*) "I am unsteadyFlow"
#endif
    write(00,*) "OpenACC GPU version; final wall-fit and psi/vorticity diagnostics stay on host"
    !-----------------------------------------------------------------------------------------------



    !-----------------------------------------------------------------------------------------------
    !节点坐标数组
    xp(0) = 0.0d0
    xp(nx+1) = dble(nx)
    do i=1,nx
        xp(i) = dble(i)-0.5d0
    enddo
    xp = xp / lengthUnit

    yp(0) = 0.0d0
    yp(ny+1) = dble(ny)
    do j=1,ny
        yp(j) = dble(j)-0.5d0
    enddo
    yp = yp / lengthUnit

    allocate (u(nx,ny))
    allocate (v(nx,ny))
    allocate (T(nx,ny))
    allocate (rho(nx,ny))

#ifdef steadyFlow
    allocate (up(nx,ny))
    allocate (vp(nx,ny))
    allocate (Tp(nx,ny))
#endif

    allocate (f(nx,ny,0:8))
    allocate (f_post(0:nx+1,0:ny+1,0:8))
    allocate (g(nx,ny,0:4))
    allocate (g_post(0:nx+1,0:ny+1,0:4))

    allocate (Fx(nx,ny))
    allocate (Fy(nx,ny))
    allocate (Sxx(nx,ny), Sxy(nx,ny), Syy(nx,ny), Sdiv(nx,ny))

#ifdef unsteadyFlow
    allocate (sumTemperatureXAvgProfile(ny), sumTemperatureSquaredXAvgProfile(ny))
#endif

    !-----------------------------------------------------------------------------------------------



    !-----------------------------------------------------------------------------------------------
    !初始化
    rho = 1.0d0                                     !密度rho=1

    omega(0) = 4.0d0/9.0d0
    do alpha=1,4
        omega(alpha) = 1.0d0/9.0d0
    enddo
    do alpha=5,8
        omega(alpha) = 1.0d0/36.0d0
    enddo

    ! D2Q5 温度权重，paraA 的 -4<paraA<1 限制保证所有权重为正。
    omegaT(0) = (1.0d0-paraA)/5.0d0
    do alpha=1,4
        omegaT(alpha) = (paraA+4.0d0)/20.0d0
    enddo

    if(loadInitField.EQ.0) then                    !在不加载文件的情况下，都是零场为初值

        u = 0.0d0
        v = 0.0d0
        T = 0.0d0
        Fx = 0.0d0
        Fy = 0.0d0
        Sxx = 0.0d0
        Sxy = 0.0d0
        Syy = 0.0d0
        Sdiv = 0.0d0

        write(00,*) "Initial field is set exactly"
        if(restartTfOffset.NE.0.0d0) then !在不加载文件的情况下，restartTfOffset必须是零
            write(00,*) "Error: since loadInitField.EQ.0, restartTfOffset should also be 0"
            stop
        endif

#ifdef VerticalWallsNoslip
        write(00,*) "Velocity B.C. for vertical walls are: ===No-slip wall==="
#endif
#ifdef VerticalWallsPeriodicalU
        write(00,*) "Velocity B.C. for vertical walls are: ===Periodical==="
#endif
#ifdef HorizontalWallsNoslip
        write(00,*) "Velocity B.C. for horizontal walls are: ===No-slip wall==="
#endif

#ifdef VerticalWallsConstT
    do j = 1, ny                                   !在不加载文件的情况下，初始化温度是分层的，xp(0)和xp(nx+1)分别是左侧和右侧边界的坐标，xp(i)是内部节点的坐标，线性插值出每个节点的初始温度
        do i = 1, nx
            T(i,j) = Thot + (xp(i)-xp(0)) / (xp(nx+1)-xp(0)) * (Tcold-Thot)
        enddo
    enddo
    write(00,*) "Temperature B.C. for vertical walls are:===Hot/cold wall==="
#endif

#ifdef HorizontalWallsConstT
    do i = 1, nx
        do j = 1, ny
            T(i,j) = Thot + (yp(j)-yp(0)) / (yp(ny+1)-yp(0)) * (Tcold-Thot)
        enddo
    enddo
    write(00,*) "Temperature B.C. for horizontal walls are:===Hot/cold wall==="
#endif

#ifdef VerticalWallsAdiabatic
    write(00,*) "Temperature B.C. for vertical walls are:===Adiabatic wall==="
#endif

#ifdef HorizontalWallsAdiabatic
    write(00,*) "Temperature B.C. for horizontal walls are:===Adiabatic wall==="
#endif

        f = 0.0d0
        g = 0.0d0
        do j = 1,ny
            do i = 1,nx
                us2 = u(i,j)*u(i,j)+v(i,j)*v(i,j)
                do alpha = 0, 8
                    un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                    f(i,j,alpha) = omega(alpha)*((rho(i,j)-rho0) + &
                        rho0*(un(alpha)/cs2+0.5d0*un(alpha)*un(alpha)/(cs2*cs2)-0.5d0*us2/cs2))
                enddo
                do alpha = 0, 4
                    un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                    g(i,j,alpha) = omegaT(alpha)*T(i,j)*(1.0d0+thermalGeqCoeff*un(alpha))
                enddo
            enddo
        enddo
    elseif(loadInitField.EQ.1) then                               !续算分支：从旧的严格重启文件恢复 f/g 和输出计数
    ! 正常断电续算时，先读取 <reloadFilePrefix>-latest.meta；
    ! meta 会告诉代码实际要读哪个 <reloadFilePrefix>-*.bin，以及旧算例已经累计到哪里。
    ! 正常续算时 restartTfOffset 可保持为 0，实际值由 latest.meta 恢复。
        write(00,*) "Read reload metadata before choosing the restart .bin file."
        write(reloadFileName,'(i12.12)') reloadFileNum             !latest .meta 缺失时才依赖这个手工编号
        reloadFileName = adjustl(reloadFileName)                  !adjustl把字符串左对齐，把前导空格移到字符串末尾
        !latest .meta 存在时：精确读取账本并更新 reloadFileName；不存在时保守推断。
        call read_reload_metadata(reloadFileName)
        write(00,*) "Load initial field from previous simulation: ", &
            trim(reloadFilePrefix), "-", trim(reloadFileName), ".bin"
        if(.not.reloadMetadataLoaded) then
            write(00,*) "WARNING: no reload metadata file found; restart offsets were inferred."
            write(00,*) "         For exact continuation, use reload files written after this patch."
        endif
        open(unit=01,file=trim(reloadFilePrefix)//"-"//trim(reloadFileName)//".bin",form="unformatted", &
        access="sequential",status='old')  !unformatted是二进制,sequential：按记录顺序读写
            ! 严格重启文件只保存完整的 f、g 分布函数。
            write(00,*) "Reloading f and g from file"
            read(01) (((f(i,j,alpha), i=1,nx), j=1,ny), alpha=0,8)      !先 i，再 j，再 alpha
            read(01) (((g(i,j,alpha), i=1,nx), j=1,ny), alpha=0,4)
        close(01)
        call reconstruct_macro_from_fg()
        write(00,*) "Raw data is loaded from the file: ", trim(reloadFilePrefix), "-", trim(reloadFileName), ".bin"
        write(00,*) "Restart offset itc =", restartItcOffset
        write(00,*) "Restart offset time_tf =", real(restartTfOffset,kind=8)
#ifdef unsteadyFlow
        write(00,*) "Restart cumulative statistic samples =", cumulativeStatisticSampleCount
#endif
        write(00,*) "Continue output counters: snapshot/plt/reload =", snapshotFileNum, pltFileNum, reloadFileNum
    else
        write(00,*) "Error: initial field is not properly set"  !如果 loadInitField 不是 0/1 或逻辑不一致，直接停止
        stop
    endif

    write(00,*)"-------------------------------------------------------------------------------"
close(00)

    f_post = 0.0d0
    g_post = 0.0d0

    if(loadInitField.EQ.0) then
        snapshotFileNum = 0
        pltFileNum = 0
        reloadFileNum = 0
        restartItcOffset = 0
        restartTfOffset = 0.0d0
#ifdef steadyFlow
        ! 新算例第一段收敛误差应从初始场开始比较。
        up = u
        vp = v
        Tp = T
#endif
    else
#ifdef steadyFlow
        ! 重启后第一段收敛误差应从载入场继续比较。
        up = u
        vp = v
        Tp = T
#endif
    endif
    ! 新算例先清零；续算也先清零，随后从耗散与温度剖面历史恢复到重启时刻。
    currentRunStatisticSampleCount = 0
    NuVolAvg = 0.0d0
    ReVolAvg = 0.0d0
#ifdef unsteadyFlow
    speedSquaredVolAvg = 0.0d0
    epsKineticVolAvg = 0.0d0
    epsThermalVolAvg = 0.0d0
    densityFluctuationSquaredVolAvg = 0.0d0
    velocityDivergenceSquaredVolAvg = 0.0d0
    maxMachLocal = 0.0d0
    minTemperature = 0.0d0
    maxTemperature = 0.0d0
    sumNuVolAvg = 0.0d0
    sumSpeedSquaredVolAvg = 0.0d0
    sumEpsKineticVolAvg = 0.0d0
    sumEpsThermalVolAvg = 0.0d0
    sumDensityFluctuationSquaredVolAvg = 0.0d0
    sumVelocityDivergenceSquaredVolAvg = 0.0d0
    maxStatisticCFL = 0.0d0
    sumTemperatureXAvgProfile = 0.0d0
    sumTemperatureSquaredXAvgProfile = 0.0d0
    if(loadInitField.EQ.0) cumulativeStatisticSampleCount = 0
    if(loadInitField.EQ.1) call restore_dissipation_statistics()
#endif

    return
  end subroutine initial
!===================================================================================================


!===================================================================================================
! 子程序: enter_data_2d_openacc
! 作用: 在主时间推进前把主要数组和常量映射到 OpenACC 设备端。
!===================================================================================================
  subroutine enter_data_2d_openacc()
    use openacc
    use commondata
    implicit none

    !$acc enter data copyin(xp,yp,ex,ey,omega,omegaT)
    !$acc enter data copyin(u,v,T,rho,f,g,Fx,Fy,Sxx,Sxy,Syy,Sdiv)
    !$acc enter data create(f_post,g_post)
#ifdef steadyFlow
    !$acc enter data copyin(up,vp,Tp)
#endif
  end subroutine enter_data_2d_openacc
!===================================================================================================

!===================================================================================================
! 子程序: update_host_snapshot_2d_openacc
! 作用: 在主机端输出 uvTrho 快照或做 CPU 后处理前，同步宏观场。
!===================================================================================================
  subroutine update_host_snapshot_2d_openacc()
    use commondata
    implicit none

    !$acc wait(1)
    !$acc update self(u,v,T,rho)
  end subroutine update_host_snapshot_2d_openacc
!===================================================================================================

!===================================================================================================
! 子程序: update_host_tecplot_2d_openacc
! 作用: Tecplot 主场输出只需要 u、v、T。
!===================================================================================================
  subroutine update_host_tecplot_2d_openacc()
    use commondata
    implicit none

    !$acc wait(1)
    !$acc update self(u,v,T)
  end subroutine update_host_tecplot_2d_openacc
!===================================================================================================

!===================================================================================================
! 子程序: update_host_reload_2d_openacc
! 作用: 严格重启文件输出前，把完整的 f、g 分布函数同步回主机。
!===================================================================================================
  subroutine update_host_reload_2d_openacc()
    use commondata
    implicit none

    !$acc wait(1)
    !$acc update self(f,g)
  end subroutine update_host_reload_2d_openacc
!===================================================================================================

!===================================================================================================
! 子程序: exit_data_2d_openacc
! 作用: 在计算结束后释放设备端常驻数据。
!===================================================================================================
  subroutine exit_data_2d_openacc()
    use commondata
    implicit none

#ifdef steadyFlow
    !$acc exit data delete(up,vp,Tp)
#endif
    !$acc exit data delete(f_post,g_post,u,v,T,rho,f,g,Fx,Fy,Sxx,Sxy,Syy,Sdiv)
    !$acc exit data delete(xp,yp,ex,ey,omega,omegaT)
  end subroutine exit_data_2d_openacc
!===================================================================================================



!===================================================================================================
! 子程序: collision
! 作用: 在一个格点循环中计算力、非平衡应变率，并完成矩空间 D2Q9 TRT LBM-CDE 碰撞。
! 说明: m0/m3/m5 是密度和动量守恒矩，松弛率为零；TRT 时其余偶矩使用 Snu、奇矩使用 Sq；
!       FLOW_BGK 时 Sq=Snu，因此全部非守恒矩使用同一 tau_0。
!       A_ab=chi_nu*S_ab+(chi_b-chi_nu)*Sdiv*delta_ab/2，并作为矩空间源项的一部分加入。
! 动力黏度为 mu=rho0*nu。
!===================================================================================================
  subroutine collision()
    use commondata
    implicit none
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: m(0:8), meq(0:8), s(0:8), mSource(0:8), m_post(0:8)
    real(kind=8) :: densityFluctuation, uu, uF
    real(kind=8) :: neqTrace, neqxx, neqxy, neqyy
    real(kind=8) :: Axx, Axy, Ayy
    real(kind=8), parameter :: muShear=rho0*viscosity
    real(kind=8), parameter :: muBulk=rho0*bulkViscosity
    real(kind=8), parameter :: denomDiag=2.0d0*muShear+rho0*cs2
    real(kind=8), parameter :: denomShear=4.0d0*muShear+2.0d0*rho0*cs2
    real(kind=8), parameter :: denomDiv=2.0d0*muBulk+rho0*cs2
    real(kind=8), parameter :: coeffTrace=(2.0d0*muShear-2.0d0*muBulk)/(2.0d0*denomDiag)

    !$acc parallel loop gang vector collapse(2) default(none) &
    !$acc& present(f,f_post,rho,u,v,T,Fx,Fy,Sxx,Sxy,Syy,Sdiv) async(1) &
    !$acc& private(alpha,m,meq,s,mSource,m_post,densityFluctuation,uu,uF, &
    !$acc& neqTrace,neqxx,neqxy,neqyy,Axx,Axy,Ayy)
    do j = 1, ny
        do i = 1, nx
            ! Boussinesq 力使用参考密度 rho0。
            Fx(i,j) = 0.0d0
            Fy(i,j) = rho0*gBeta*(T(i,j)-Tref)
#ifdef SideHeatedHa
            Fx(i,j) = rho0*B2sigemarho*(v(i,j)*sin(phi)*cos(phi)-u(i,j)*sin(phi)*sin(phi))
            Fy(i,j) = Fy(i,j)+rho0*B2sigemarho*(u(i,j)*sin(phi)*cos(phi)-v(i,j)*cos(phi)*cos(phi))
#endif

            ! D2Q9 分布函数变换到矩空间。
            m(0) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+ &
                f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
            m(1) = -4.0d0*f(i,j,0)-f(i,j,1)-f(i,j,2)-f(i,j,3)-f(i,j,4)+ &
                2.0d0*(f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8))
            m(2) = 4.0d0*f(i,j,0)-2.0d0*(f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4))+ &
                f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
            m(3) = f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)
            m(4) = -2.0d0*f(i,j,1)+2.0d0*f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)
            m(5) = f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)
            m(6) = -2.0d0*f(i,j,2)+2.0d0*f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)
            m(7) = f(i,j,1)-f(i,j,2)+f(i,j,3)-f(i,j,4)
            m(8) = f(i,j,5)-f(i,j,6)+f(i,j,7)-f(i,j,8)

            densityFluctuation = rho(i,j)-rho0
            uu = u(i,j)*u(i,j)+v(i,j)*v(i,j)
            meq(0) = densityFluctuation
            meq(1) = -2.0d0*densityFluctuation+3.0d0*rho0*uu
            meq(2) = densityFluctuation-3.0d0*rho0*uu
            meq(3) = rho0*u(i,j)
            meq(4) = -rho0*u(i,j)
            meq(5) = rho0*v(i,j)
            meq(6) = -rho0*v(i,j)
            meq(7) = rho0*(u(i,j)*u(i,j)-v(i,j)*v(i,j))
            meq(8) = rho0*u(i,j)*v(i,j)

            ! 二阶非平衡矩：trace=(m1+4*m0)/3，m7=Pi_xx-Pi_yy，m8=Pi_xy。
            neqTrace = ((m(1)-meq(1))+4.0d0*(m(0)-meq(0)))/3.0d0
            neqxx = 0.5d0*(neqTrace+m(7)-meq(7))
            neqyy = 0.5d0*(neqTrace-m(7)+meq(7))
            neqxy = m(8)-meq(8)

            uF = u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j)
            Sdiv(i,j) = -(neqxx+neqyy+uF)/denomDiv
            Sxx(i,j) = -(neqxx+u(i,j)*Fx(i,j))/denomDiag+coeffTrace*Sdiv(i,j)
            Syy(i,j) = -(neqyy+v(i,j)*Fy(i,j))/denomDiag+coeffTrace*Sdiv(i,j)
            Sxy(i,j) = -(2.0d0*neqxy+u(i,j)*Fy(i,j)+v(i,j)*Fx(i,j))/denomShear

            Axx = chi_nu*Sxx(i,j)+0.5d0*(chi_b-chi_nu)*Sdiv(i,j)
            Ayy = chi_nu*Syy(i,j)+0.5d0*(chi_b-chi_nu)*Sdiv(i,j)
            Axy = chi_nu*Sxy(i,j)

            ! 原始矩源项 = Guo 力项 + LBM-CDE 二阶 Hermite 应力修正。
            mSource(0) = 0.0d0
            mSource(1) = 6.0d0*uF+2.0d0*rho0*(Axx+Ayy)
            mSource(2) = -6.0d0*uF-2.0d0*rho0*(Axx+Ayy)
            mSource(3) = Fx(i,j)
            mSource(4) = -Fx(i,j)
            mSource(5) = Fy(i,j)
            mSource(6) = -Fy(i,j)
            mSource(7) = 2.0d0*(u(i,j)*Fx(i,j)-v(i,j)*Fy(i,j))+ &
                2.0d0/3.0d0*rho0*(Axx-Ayy)
            mSource(8) = u(i,j)*Fy(i,j)+v(i,j)*Fx(i,j)+2.0d0/3.0d0*rho0*Axy

            ! m0/m3/m5 为守恒矩；TRT 时非守恒奇矩 m4/m6 用 Sq；FLOW_BGK 时 Sq=Snu。
            s(0) = 0.0d0
            s(1) = Snu
            s(2) = Snu
            s(3) = 0.0d0
            s(4) = Sq
            s(5) = 0.0d0
            s(6) = Sq
            s(7) = Snu
            s(8) = Snu

            do alpha = 0, 8
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+ &
                    (1.0d0-0.5d0*s(alpha))*mSource(alpha)
            enddo

            ! 乘 M^{-1} 回到速度空间。
            f_post(i,j,0) = m_post(0)/9.0d0-m_post(1)/9.0d0+m_post(2)/9.0d0
            f_post(i,j,1) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+ &
                m_post(3)/6.0d0-m_post(4)/6.0d0+m_post(7)/4.0d0
            f_post(i,j,2) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+ &
                m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)/4.0d0
            f_post(i,j,3) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0- &
                m_post(3)/6.0d0+m_post(4)/6.0d0+m_post(7)/4.0d0
            f_post(i,j,4) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0- &
                m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)/4.0d0
            f_post(i,j,5) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+ &
                m_post(3)/6.0d0+m_post(4)/12.0d0+m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)/4.0d0
            f_post(i,j,6) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0- &
                m_post(3)/6.0d0-m_post(4)/12.0d0+m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)/4.0d0
            f_post(i,j,7) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0- &
                m_post(3)/6.0d0-m_post(4)/12.0d0-m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)/4.0d0
            f_post(i,j,8) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+ &
                m_post(3)/6.0d0+m_post(4)/12.0d0-m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)/4.0d0
        enddo
    enddo
    return
  end subroutine collision
!===================================================================================================


!===================================================================================================
! 子程序: streaming
! 作用: 完成流场分布函数 f 的迁移，把碰撞后的信息传播到相邻格点。
! 用途: 在主程序时间推进循环中调用，位于 collision 之后、bounceback 之前。
!===================================================================================================
  subroutine streaming()                                    !先迁移，再边界处理
    use commondata                                            !迁移步骤：pull streaming，把碰撞后的 f_post 拉取到当前格点
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: ip, jp
    integer(kind=4) :: alpha

    !$acc parallel loop gang vector collapse(2) default(none) present(f,f_post,ex,ey) async(1) private(alpha,ip,jp)
    do j = 1, ny
        do i = 1, nx
            do alpha = 0, 8                        !上游格点索引：fα(i,j) <- f_postα(i-exα, j-eyα)
                ip = i-ex(alpha)                   !边界附近 (ip/jp 可能为 0 或 nx+1/ny+1)，需在 bounceback/周期边界处理中覆盖修正边界分布
                jp = j-ey(alpha)                   !ghost 层在初始化中为 0，保证不会出现未初始化垃圾值

                f(i,j,alpha) = f_post(ip,jp,alpha)
            enddo
        enddo
    enddo
    return
  end subroutine streaming
!===================================================================================================
! streaming 结束: 完成流场分布函数 f 的迁移，把碰撞后的信息传播到相邻格点。
!===================================================================================================




!===================================================================================================
! 子程序: bounceback
! 作用: 处理流场边界条件，包括无滑移壁面和相关反弹格式。
! 用途: 在主程序时间推进循环中调用，位于 streaming 之后、macro 之前。
!===================================================================================================
  subroutine bounceback()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    ! integer(kind=4) :: alpha

#ifdef VerticalWallsPeriodicalU
    !$acc parallel loop gang vector default(none) present(f,f_post) async(1)
    do j = 1, ny                                                  !速度边界垂直边界周期，直接方向相同，跨边界的入射分布
        !Left side (i=1)
        f(1,j,1) = f_post(nx,j,1)
        f(1,j,5) = f_post(nx,j,5)
        f(1,j,8) = f_post(nx,j,8)

        !Right side (i=nx)
        f(nx,j,3) = f_post(1,j,3)
        f(nx,j,6) = f_post(1,j,6)
        f(nx,j,7) = f_post(1,j,7)
    enddo
#endif

#ifdef VerticalWallsNoslip
    !$acc parallel loop gang vector default(none) present(f,f_post) async(1)
    do j = 1, ny                                                 !速度边界垂直边界静止壁无滑移，直接反弹，方向相反
        !Left side (i=1)
        f(1,j,1) = f_post(1,j,3)
        f(1,j,5) = f_post(1,j,7)
        f(1,j,8) = f_post(1,j,6)

        !Right side (i=nx)
        f(nx,j,3) = f_post(nx,j,1)
        f(nx,j,6) = f_post(nx,j,8)
        f(nx,j,7) = f_post(nx,j,5)
    enddo
#endif

#ifdef HorizontalWallsNoslip
    !$acc parallel loop gang vector default(none) present(f,f_post) async(1)
    do i = 1, nx                                                  !速度边界水平边界无滑移，直接反弹，方向相反
        !Bottom side (j=1)
        f(i,1,2) = f_post(i,1,4)
        f(i,1,5) = f_post(i,1,7)
        f(i,1,6) = f_post(i,1,8)

        !Top side (j=ny)
        f(i,ny,4) = f_post(i,ny,2)
        f(i,ny,7) = f_post(i,ny,5)
        f(i,ny,8) = f_post(i,ny,6)
    enddo
#endif

    return
  end subroutine bounceback
!===================================================================================================
! bounceback 结束: 处理流场边界条件，包括无滑移壁面和相关反弹格式。
!===================================================================================================



!===================================================================================================

!===================================================================================================
! 子程序: macro
! 作用: 由流场分布函数恢复 rho、u、v，并加入半步力项速度修正。
!===================================================================================================
  subroutine macro()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !$acc parallel loop gang vector collapse(2) default(none) present(f,rho,u,v,Fx,Fy) async(1)
    do j = 1, ny
        do i = 1, nx
            rho(i,j) = rho0+f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+ &
                f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
            u(i,j) = (f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8) + &
                0.5d0*Fx(i,j))/rho0
            v(i,j) = (f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8) + &
                0.5d0*Fy(i,j))/rho0
        enddo
    enddo
    return
  end subroutine macro
!===================================================================================================

!===================================================================================================



!===================================================================================================
! 子程序: collisionT
! 作用: 完成固定的 Luo D2Q5-TRT 温度分布函数碰撞更新。
! 用途: 在主程序时间推进循环中调用，位于流场 macro 之后。
!===================================================================================================
    subroutine collisionT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: n(0:4), n_post(0:4), neq(0:4)
    real(kind=8) :: q(0:4)




    !$acc parallel loop gang vector collapse(2) default(none) present(g,g_post,u,v,T) async(1) &
    !$acc& private(alpha,n,neq,q,n_post)
    do j = 1, ny
        do i = 1, nx
          n(0) = g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
          n(1) = g(i,j,1)-g(i,j,3)
          n(2) = g(i,j,2)-g(i,j,4)
          n(3) = -4.0d0*g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
          n(4) = g(i,j,1)-g(i,j,2)+g(i,j,3)-g(i,j,4)

          ! 平衡矩：n_eq(0) 保证温度守恒，n_eq(1:2) 输运 T*u/T*v；
          ! n_eq(3) 由 paraA 控制热扩散映射，n_eq(4)=0 不携带宏观物理量。
          neq(0) = T(i,j)
          neq(1) = T(i,j)*u(i,j)
          neq(2) = T(i,j)*v(i,j)
          neq(3) = T(i,j)*paraA
          neq(4) = 0.0d0

          ! 各矩松弛率：温度守恒矩为 0，两个通量矩共用 Qk，其余分别使用 Qe、Qnu。
          q(0) = 0.0d0
          q(1) = Qk
          q(2) = Qk
          q(3) = Qe
          q(4) = Qnu


          n_post(0) = n(0)-q(0)*(n(0)-neq(0))
          n_post(1) = n(1)-q(1)*(n(1)-neq(1))
          n_post(2) = n(2)-q(2)*(n(2)-neq(2))
          n_post(3) = n(3)-q(3)*(n(3)-neq(3))
          n_post(4) = n(4)-q(4)*(n(4)-neq(4))


          g_post(i,j,0) = 0.2d0*n_post(0)-0.2d0*n_post(3)
          g_post(i,j,1) = 0.2d0*n_post(0)+0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
          g_post(i,j,2) = 0.2d0*n_post(0)+0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4)
          g_post(i,j,3) = 0.2d0*n_post(0)-0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
          g_post(i,j,4) = 0.2d0*n_post(0)-0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4)
        enddo
    enddo
    return
    end subroutine collisionT
!===================================================================================================
! collisionT 结束: 完成固定的 Luo D2Q5-TRT 温度碰撞更新。
!===================================================================================================




!===================================================================================================
! 子程序: streamingT
! 作用: 完成温度分布函数 g 的迁移，把碰撞后的温度信息传播到相邻格点。
! 用途: 在主程序时间推进循环中调用，位于 collisionT 之后、bouncebackT 之前。
!===================================================================================================
    subroutine streamingT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: ip, jp
    integer(kind=4) :: alpha

    !$acc parallel loop gang vector collapse(2) default(none) present(g,g_post,ex,ey) async(1) private(alpha,ip,jp)
    do j = 1, ny
        do i = 1, nx
            do alpha = 0, 4
                ip = i-ex(alpha)
                jp = j-ey(alpha)

                g(i,j,alpha) = g_post(ip,jp,alpha)
            enddo
        enddo
    enddo
    return
    end subroutine streamingT
!===================================================================================================
! streamingT 结束: 完成温度分布函数 g 的迁移，把碰撞后的温度信息传播到相邻格点。
!===================================================================================================



!===================================================================================================
! 子程序: bouncebackT
! 作用: 处理温度边界条件，包括恒温、绝热和周期边界。
! 用途: 在主程序时间推进循环中调用，位于 streamingT 之后、macroT 之前。
!===================================================================================================
    subroutine bouncebackT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    !integer(kind=4) :: alpha

#ifdef VerticalWallsPeriodicalT
    !$acc parallel loop gang vector default(none) present(g,g_post) async(1)
    do j = 1, ny
        !Left boundary
        g(1,j,1) = g_post(nx,j,1)

        !Right boundary
        g(nx,j,3) = g_post(1,j,3)
    enddo
#endif

#ifdef VerticalWallsConstT
    !$acc parallel loop gang vector default(none) present(g,g_post,omegaT) async(1)
    do j = 1, ny
        !Left boundary
        g(1,j,1) = -g_post(1,j,3)+(4.0d0+paraA)/10.0d0*Thot
        !Right boundary
        g(nx,j,3) = -g_post(nx,j,1)+(4.0d0+paraA)/10.0d0*Tcold
    enddo
#endif

#ifdef VerticalWallsAdiabatic
    !$acc parallel loop gang vector default(none) present(g,g_post) async(1)
    do j = 1, ny
        !Left boundary
        g(1,j,1) = g_post(1,j,3)

        !Right boundary
        g(nx,j,3) = g_post(nx,j,1)
    enddo
#endif

#ifdef HorizontalWallsAdiabatic
    !$acc parallel loop gang vector default(none) present(g,g_post) async(1)
    do i = 1, nx
        !Bottom side
        g(i,1,2) = g_post(i,1,4)

        !Top side
        g(i,ny,4) = g_post(i,ny,2)
    enddo
#endif

#ifdef HorizontalWallsConstT
    !$acc parallel loop gang vector default(none) present(g,g_post,omegaT) async(1)
    do i = 1, nx
        !Bottom side
        g(i,1,2) = -g_post(i,1,4)+(4.0d0+paraA)/10.0d0*Thot
        !Top side
        g(i,ny,4) = -g_post(i,ny,2)+(4.0d0+paraA)/10.0d0*Tcold
    enddo
#endif

    return
    end subroutine bouncebackT
!===================================================================================================
! bouncebackT 结束: 处理温度边界条件，包括恒温、绝热和周期边界。
!===================================================================================================




!===================================================================================================
! 子程序: macroT
! 作用: 由温度分布函数恢复宏观温度场 T。
! 用途: 在主程序时间推进循环中调用，作为温度更新链条的最后一步。
!===================================================================================================
    subroutine macroT()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !$acc parallel loop gang vector collapse(2) default(none) present(g,T) async(1)
    do j = 1, ny
        do i = 1, nx
            T(i,j) = g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
        enddo
    enddo
    return
    end subroutine macroT
!===================================================================================================
! macroT 结束: 由温度分布函数恢复宏观温度场 T。
!===================================================================================================



!===================================================================================================
! 子程序: reconstruct_macro_from_fg
! 作用: 从重启读回的 f/g 重新恢复 rho、u、v 和 T。
!===================================================================================================
    subroutine reconstruct_macro_from_fg()
    use commondata
    implicit none
    integer(kind=4) :: i, j, iter
    real(kind=8) :: momx, momy
    logical :: rho_bad

    rho_bad = .false.
    do j = 1, ny
        do i = 1, nx
            T(i,j) = g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
            rho(i,j) = rho0+f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+ &
                f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
            momx = f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)
            momy = f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)

            if (rho(i,j).GT.0.0d0) then
                u(i,j) = momx/rho0
                v(i,j) = momy/rho0

                do iter = 1, 3
                    Fx(i,j) = 0.0d0
                    Fy(i,j) = rho0*gBeta*(T(i,j)-Tref)

#ifdef    SideHeatedHa
                    Fx(i,j) = rho0*B2sigemarho*(v(i,j)*sin(phi)*cos(phi)-u(i,j)*sin(phi)*sin(phi))
                    Fy(i,j) = rho0*gBeta*(T(i,j)-Tref)+rho0*B2sigemarho*(u(i,j)*sin(phi)*cos(phi)&
                    -v(i,j)*cos(phi)*cos(phi))
#endif

                    u(i,j) = (momx + 0.5d0*Fx(i,j))/rho0
                    v(i,j) = (momy + 0.5d0*Fy(i,j))/rho0
                enddo
            else
                rho_bad = .true.
                Fx(i,j) = 0.0d0
                Fy(i,j) = 0.0d0
                u(i,j) = 0.0d0
                v(i,j) = 0.0d0
            endif
        enddo
    enddo

    if (rho_bad) then
        write(*,*) "Warning: non-positive rho found during restart reconstruction."
        stop
    endif

    return
    end subroutine reconstruct_macro_from_fg
!===================================================================================================
! reconstruct_macro_from_fg end: 重新恢复 rho、u、v 和 T
!===================================================================================================



#ifdef steadyFlow
!===================================================================================================
! 子程序: check
! 作用: 计算稳态收敛误差并写出收敛历史。
! 用途: 在 steadyFlow 模式下由主程序定期调用。
!===================================================================================================
    subroutine check()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: error1, error2, error5, error6
    character(len=64) :: caseTag



    !$acc wait(1)
    error1 = 0.0d0
    error2 = 0.0d0

    error5 = 0.0d0
    error6 = 0.0d0

    !$acc parallel loop collapse(2) default(none) present(u,up,v,vp,T,Tp) reduction(+:error1,error2,error5,error6)
    do j = 1, ny
        do i = 1, nx
            error1 = error1+(u(i,j)-up(i,j))*(u(i,j)-up(i,j))+(v(i,j)-vp(i,j))*(v(i,j)-vp(i,j))
            error2 = error2+u(i,j)*u(i,j)+v(i,j)*v(i,j)

            error5 = error5+dABS( T(i,j)-Tp(i,j) )
            error6 = error6+dABS( T(i,j) )

            up(i,j) = u(i,j)
            vp(i,j) = v(i,j)
            Tp(i,j) = T(i,j)
        enddo
    enddo
    if (error2 .GT. 1.0d-30) then
        errorU = dsqrt(error1)/dsqrt(error2)                 !速度场相对L2误差：||u^n-u^{n-1}||_2 / ||u^n||_2
    else
        errorU = dsqrt(error1)
    endif
    if (error6 .GT. 1.0d-30) then
        errorT = error5/error6                               !温度场相对L1误差：||T^n-T^{n-1}||_1 / ||T^n||_1
    else
        errorT = error5
    endif



    call append_convergence_tecplot('convergence2DOpenacc.plt', restartItcOffset+itc, errorU, errorT)


    write(caseTag,'("Ra=",ES24.16E3,",nx=",I0,",ny=",I0,",thermal=LuoD2Q5TRT")') &
        Rayleigh, nx, ny  !输出收敛曲线的对比
    call append_convergence_master_tecplot('convergence_all_2DOpenacc.plt', caseTag, restartItcOffset+itc, errorU, errorT)

    write(*,'(I12,1X,ES24.16E3,1X,ES24.16E3)') restartItcOffset+itc, errorU, errorT


    return
    end subroutine check
!===================================================================================================
! check 结束: 计算稳态收敛误差并写出收敛历史。
!===================================================================================================
#endif



!===================================================================================================
! 子程序: append_convergence_tecplot
! 作用: 向单个 Tecplot 收敛文件追加一条误差记录。
! 用途: 在 check 中调用，用于输出当前算例的收敛曲线。
!===================================================================================================
subroutine append_convergence_tecplot(filename, itc, errorU, errorT)
  use commondata, only: loadInitField
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4), intent(in)  :: itc
  real(kind=8),    intent(in)  :: errorU, errorT

  integer :: u
  logical :: ex
  logical, save :: first_write = .true.      !这个局部变量在子程序返回后不会被销毁，下一次再进这个子程序时还保留上一次的值

  if (first_write) then
    inquire(file=trim(filename), exist=ex)
    if ((loadInitField.EQ.1).AND.ex) then
      ! 续算：旧收敛曲线继续追加，横坐标已经传入累计 itc。
      open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
    elseif(loadInitField.EQ.1) then
      ! 续算要求旧收敛文件必须存在；否则会丢掉断电前的收敛历史。
      write(*,*) 'Error: restart requested but convergence file is missing: ', trim(filename)
      stop
    else
      ! 新算例：清掉旧历史，避免不同算例的数据混在一起。
      open(newunit=u, file=trim(filename), status='replace', action='write', form='formatted')
      write(u,'(A)') 'VARIABLES = "itc" "errorU" "errorT"'
      write(u,'(A)') 'ZONE T="conv", F=POINT'
    endif
    write(u,'(I12,1X,ES24.16E3,1X,ES24.16E3)') itc, errorU, errorT
    close(u)

    first_write = .false.
  else
    ! 同一次运行的后续调用：追加数据行
    open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
    write(u,'(I12,1X,ES24.16E3,1X,ES24.16E3)') itc, errorU, errorT
    close(u)
  end if

end subroutine append_convergence_tecplot
!===================================================================================================
! append_convergence_tecplot 结束: 向单个 Tecplot 收敛文件追加一条误差记录。
!===================================================================================================



!===================================================================================================
! 子程序: append_convergence_master_tecplot
! 作用: 向总收敛文件追加一条带算例标签的记录。
! 用途: 预留给 check 调用，用于汇总多组算例的收敛信息。
!===================================================================================================
subroutine append_convergence_master_tecplot(filename, zoneName, itc, errorU, errorT)
  use commondata, only: loadInitField
  implicit none
  character(len=*), intent(in) :: filename, zoneName
  integer(kind=4), intent(in)  :: itc
  real(kind=8),    intent(in)  :: errorU, errorT

  logical :: ex
  integer :: u
  logical, save :: zone_started = .false.

  ! 本次运行第一次写：新算例写入新 ZONE；续算则接到旧 ZONE 后面
  if (.not. zone_started) then
    inquire(file=trim(filename), exist=ex)

    if (.not. ex) then
      if(loadInitField.EQ.1) then
        write(*,*) 'Error: restart requested but master convergence file is missing: ', trim(filename)
        stop
      endif
      open(newunit=u, file=trim(filename), status='new', action='write', form='formatted')
      write(u,'(A)') 'TITLE = "Convergence comparison"'
      write(u,'(A)') 'VARIABLES = "itc" "errorU" "errorT"'
      close(u)
    end if

    if(loadInitField.EQ.0) then
      open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
      write(u,'(A)') 'ZONE T="'//trim(zoneName)//'", F=POINT'
      close(u)
    endif

    zone_started = .true.
  end if

  ! 追加一个数据点
  open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
  write(u,'(I12,1X,ES24.16E3,1X,ES24.16E3)') itc, errorU, errorT
  close(u)
end subroutine append_convergence_master_tecplot
!===================================================================================================
! append_convergence_master_tecplot 结束: 向总收敛文件追加一条带算例标签的记录。
!===================================================================================================





!===================================================================================================
! 子程序: output_SnapshotFile
! 作用: 输出 u、v、T、rho 的二进制快照文件，供后处理分析使用。
!===================================================================================================
  subroutine output_SnapshotFile()                                   !输出 uvTrho 二进制快照
    use commondata                                                   !用于后处理快照
    implicit none
    integer(kind=4) :: i, j
    character(len=100) :: filename
    ! This snapshot is for post-processing only; u/v are written after nondimensionalization.
    ! For strict restart, keep using output_ReloadFile(), which preserves the lattice-state variables.

#ifdef steadyFlow
    write(filename,'(i12.12)') restartItcOffset+itc                  !steadyFlow：续算时使用累计格子步编号
#endif

#ifdef unsteadyFlow
    snapshotFileNum = snapshotFileNum+1
    write(filename,'(i12.12)') snapshotFileNum                       !unsteadyFlow：快照文件按调用次数编号，与 reloadFileNum 分离
#endif

    filename = adjustl(filename)

    open(unit=03,file=trim(snapshotFilePrefix)//"-"//trim(filename)//'.bin',form="unformatted",access="sequential")    !二进制
    ! Post-processing snapshot only: write nondimensionalized u/v together with T and rho.
    ! Do not use this file for strict restart; output_ReloadFile() keeps lattice velocities for that purpose.
    write(03) ((real(velocityScaleCompare*u(i,j),kind=8),i=1,nx),j=1,ny)
    write(03) ((real(velocityScaleCompare*v(i,j),kind=8),i=1,nx),j=1,ny)
    write(03) ((real(T(i,j),kind=8),i=1,nx),j=1,ny)
    write(03) ((real(rho(i,j),kind=8), i=1,nx), j=1,ny)
    close(03)

    return
  end subroutine output_SnapshotFile
!===================================================================================================
! output_SnapshotFile 结束: 输出 u、v、T、rho 的二进制快照文件。
!===================================================================================================




!===================================================================================================
! 子程序: output_ReloadFile
! 作用: 输出包含 f、g 的严格重启备份文件；rho/u/v/T 在读取后由 f/g 重建。
! 用途: 在运行过程中定期调用，也在程序结束前调用。
!===================================================================================================
  subroutine output_ReloadFile()                                  !输出fg，存储在当前路径，名字由 reloadFilePrefix 控制
    use commondata                                                !用于重启，包含 f,g；读入后调用 reconstruct_macro_from_fg() 重建宏观量
    implicit none
    integer(kind=4) :: i, j, alpha
    character(len=100) :: filename

#ifdef steadyFlow
    reloadFileNum = restartItcOffset+itc
    write(filename,'(i12.12)') reloadFileNum                 !steadyFlow：重启文件名使用累计格子步
#endif

#ifdef unsteadyFlow
    reloadFileNum = reloadFileNum + 1
    write(filename,'(i12.12)') reloadFileNum                !unsteadyFlow：reload 文件使用独立编号，不依赖快照输出是否开启
#endif

    filename = adjustl(filename)

    open(unit=05,file=trim(reloadFilePrefix)//"-"//trim(filename)//'.bin',form="unformatted",access="sequential")   !二进制
    ! 严格重启快照保存完整的 f、g 分布函数。
    write(05) (((real(f(i,j,alpha),kind=8), i=1,nx), j=1,ny), alpha=0,8)
    write(05) (((real(g(i,j,alpha),kind=8), i=1,nx), j=1,ny), alpha=0,4)
    close(05)
    call write_reload_metadata(trim(filename))

    open(unit=00,file=trim(settingsFile),status='unknown',position='append')
    write(00,*) "Backup f/g restart state to: ", trim(reloadFilePrefix), "-", trim(filename),".bin"
    write(00,*) "Backup restart metadata to: ", trim(reloadFilePrefix), "-latest.meta"
    close(00)

    return
  end subroutine output_ReloadFile
!===================================================================================================
! output_ReloadFile 结束: 输出包含 f、g 的严格重启备份文件。
!===================================================================================================


!===================================================================================================
! 子程序: write_reload_metadata
! 作用: 最新 reload 续算账本，保存累计步数、t_ff、统计样本数、输出编号和最新 .bin 文件名。
!===================================================================================================
  subroutine write_reload_metadata(filename)
    use commondata
    implicit none
    character(len=*), intent(in) :: filename
    integer(kind=4) :: metaUnit

    open(newunit=metaUnit, file=trim(reloadFilePrefix)//'-latest.meta', &
         status='replace', action='write', form='formatted')
    write(metaUnit,'(A,1X,I0)') 'reload_meta_version', 4
#ifdef steadyFlow
    write(metaUnit,'(A,1X,A)') 'flowMode', 'steadyFlow'
#endif
#ifdef unsteadyFlow
    write(metaUnit,'(A,1X,A)') 'flowMode', 'unsteadyFlow'
#endif
    write(metaUnit,'(A,1X,I0)') 'nx', nx
    write(metaUnit,'(A,1X,I0)') 'ny', ny
    write(metaUnit,'(A,1X,A)') 'reloadFileName', trim(filename)
    write(metaUnit,'(A,1X,I0)') 'itc_total', restartItcOffset+itc
    write(metaUnit,'(A,1X,ES24.16E3)') 'time_tf', &
        real(restartItcOffset+itc,kind=8)/timeUnit
#ifdef unsteadyFlow
    write(metaUnit,'(A,1X,I0)') 'cumulativeStatisticSampleCount', cumulativeStatisticSampleCount
#endif
    write(metaUnit,'(A,1X,I0)') 'snapshotFileNum', snapshotFileNum
    write(metaUnit,'(A,1X,I0)') 'pltFileNum', pltFileNum
    write(metaUnit,'(A,1X,I0)') 'reloadFileNum', reloadFileNum
    close(metaUnit)

    return
  end subroutine write_reload_metadata
!===================================================================================================


!===================================================================================================
! 子程序: read_reload_metadata
! 作用: 优先读取 latest .meta；若没有，则根据手工编号做保守推断。
!===================================================================================================
  subroutine read_reload_metadata(reloadFileName)
    use commondata
    implicit none
    character(len=*), intent(inout) :: reloadFileName
    character(len=64) :: label
    character(len=32) :: metaFlowMode, currentFlowMode
    character(len=100) :: metaReloadFileName
    character(len=256) :: metaFile
    integer(kind=4) :: metaUnit, ios
    integer(kind=4) :: metaVersion, metaNx, metaNy
    logical :: metaExists

    reloadMetadataLoaded = .false.
    metaFile = trim(reloadFilePrefix)//'-latest.meta'
    inquire(file=trim(metaFile), exist=metaExists)                 !优先检查最新账本

    if(.not.metaExists) then                                       !latest meta 不存在时，只能保守推断
        call infer_reload_offsets_without_metadata()
        return                                                     !退出当前这个子程序
    endif

    open(newunit=metaUnit, file=trim(metaFile), status='old', action='read', &
         form='formatted', iostat=ios)                             !ios == 0表示读写成功,不等于0，则代表失败
    if(ios.NE.0) then
        write(*,*) 'Error: failed to open reload metadata: ', trim(metaFile)
        stop
    endif

    read(metaUnit,*,iostat=ios) label, metaVersion
    if((ios.NE.0).OR.(trim(label).NE.'reload_meta_version').OR.(metaVersion.NE.4)) then
        write(*,*) 'Error: invalid reload metadata version in ', trim(metaFile)
        stop    !如果读取失败，或者标签/版本号不对，就说明文件格式不兼容，停止。
    endif

    read(metaUnit,*,iostat=ios) label, metaFlowMode
    if((ios.NE.0).OR.(trim(label).NE.'flowMode')) then
        write(*,*) 'Error: invalid flowMode entry in ', trim(metaFile)
        stop
    endif

    read(metaUnit,*,iostat=ios) label, metaNx
    if((ios.NE.0).OR.(trim(label).NE.'nx')) then
        write(*,*) 'Error: invalid nx entry in ', trim(metaFile)
        stop
    endif

    read(metaUnit,*,iostat=ios) label, metaNy
    if((ios.NE.0).OR.(trim(label).NE.'ny')) then
        write(*,*) 'Error: invalid ny entry in ', trim(metaFile)
        stop
    endif

    read(metaUnit,*,iostat=ios) label, metaReloadFileName
    if((ios.NE.0).OR.(trim(label).NE.'reloadFileName')) then
        write(*,*) 'Error: invalid reloadFileName entry in ', trim(metaFile)
        stop
    endif
    metaReloadFileName = adjustl(metaReloadFileName)   ! 处理左边空格，把内容左对齐

    read(metaUnit,*,iostat=ios) label, restartItcOffset
    if((ios.NE.0).OR.(trim(label).NE.'itc_total')) then
        write(*,*) 'Error: invalid itc_total entry in ', trim(metaFile)
        stop
    endif

    read(metaUnit,*,iostat=ios) label, restartTfOffset
    if((ios.NE.0).OR.(trim(label).NE.'time_tf')) then
        write(*,*) 'Error: invalid time_tf entry in ', trim(metaFile)
        stop
    endif

#ifdef unsteadyFlow
    read(metaUnit,*,iostat=ios) label, cumulativeStatisticSampleCount
    if((ios.NE.0).OR.(trim(label).NE.'cumulativeStatisticSampleCount').OR. &
       (cumulativeStatisticSampleCount.LT.0)) then
        write(*,*) 'Error: invalid cumulativeStatisticSampleCount entry in ', trim(metaFile)
        stop
    endif
#endif

    read(metaUnit,*,iostat=ios) label, snapshotFileNum
    if((ios.NE.0).OR.(trim(label).NE.'snapshotFileNum')) then
        write(*,*) 'Error: invalid snapshotFileNum entry in ', trim(metaFile)
        stop
    endif

    read(metaUnit,*,iostat=ios) label, pltFileNum
    if((ios.NE.0).OR.(trim(label).NE.'pltFileNum')) then
        write(*,*) 'Error: invalid pltFileNum entry in ', trim(metaFile)
        stop
    endif

    read(metaUnit,*,iostat=ios) label, reloadFileNum
    if((ios.NE.0).OR.(trim(label).NE.'reloadFileNum')) then
        write(*,*) 'Error: invalid reloadFileNum entry in ', trim(metaFile)
        stop
    endif
    close(metaUnit)

    currentFlowMode = 'unknown'
#ifdef steadyFlow
    currentFlowMode = 'steadyFlow'
#endif
#ifdef unsteadyFlow
    currentFlowMode = 'unsteadyFlow'
#endif

    if(trim(metaFlowMode).NE.trim(currentFlowMode)) then
        write(*,*) 'Error: reload metadata flowMode differs: ', trim(metaFlowMode), trim(currentFlowMode)
        stop
    endif
    if((metaNx.NE.nx).OR.(metaNy.NE.ny)) then
        write(*,*) 'Error: reload metadata mesh mismatch: ', metaNx, metaNy, nx, ny
        stop
    endif

    ! reloadFileNum 是整数计数器，给后续 output_ReloadFile() 继续编号，避免覆盖旧 reload 文件。
    ! reloadFileName 是字符串文件名，本次续算马上用它打开 <reloadFilePrefix>-<reloadFileName>.bin。
    reloadFileName = trim(metaReloadFileName)     ! 处理右边空格，删掉尾部空格
    reloadMetadataLoaded = .true.

    return
  end subroutine read_reload_metadata
!===================================================================================================


!===================================================================================================
! 子程序: infer_reload_offsets_without_metadata
! 作用: 没有 latest .meta 时，只能根据文件编号和当前手工参数推断。
! 根据文件名编号和当前参数“猜一个合理值”，保证续算的时间/步数是连续的。
!===================================================================================================
  subroutine infer_reload_offsets_without_metadata()
    use commondata
    implicit none

    restartItcOffset = 0
#ifdef steadyFlow
    restartItcOffset = max(0, reloadFileNum)  !稳态的 reload 文件名本来就是用 itc 写的
    if(restartTfOffset.EQ.0.0d0) then !如果没有手工给续算时间，则由累计格子步换算
        restartTfOffset = real(restartItcOffset,kind=8) / timeUnit
    endif
#endif
#ifdef unsteadyFlow
    if(restartTfOffset.EQ.0.0d0) then !非稳态的 reload 文件名不是 itc，而是第几个 reload 文件
        restartTfOffset = real(max(0,reloadFileNum),kind=8) * reloadFileInterval
    endif
    restartItcOffset = max(0, int(restartTfOffset*timeUnit+0.5d0))  !再反推格子步
    snapshotFileNum = max(0, int(restartTfOffset/outputSnapshotInterval+0.5d0)) !推断输出编号
    pltFileNum = max(0, int(restartTfOffset/outputPltFileInterval+0.5d0))
    ! 没有 latest.meta 时无法可靠知道重启场对应的统计样本数；统计窗口开始后禁止猜测。
    if(restartTfOffset.GE.unsteadyAverageStartTf) then
        write(*,*) 'Error: unsteady restart after statistics begin requires latest metadata version 4.'
        stop
    endif
#endif

    return
  end subroutine infer_reload_offsets_without_metadata
!===================================================================================================




!===================================================================================================
! 子程序: output_Tecplot
! 作用: 输出主场变量到 Tecplot 文件，便于后处理和可视化。
! 用途: 在运行过程中按需调用，也在程序结束时调用。
!===================================================================================================
  subroutine output_Tecplot()                        !输出二进制文件
    use commondata
    implicit none
    ! Here u and v are exported as nondimensional post-processing velocities using velocityScaleCompare.
    ! Restart files should still come from output_ReloadFile(), which preserves the lattice-state information.
    integer(kind=4) :: i, j, k
    REAL(kind=4) :: zoneMarker, eohMarker   !Tecplot 二进制格式里用的两个“标记值”（299 和 357），用于告诉 Tecplot：这里开始是 zone 描述 / header 结束。
    character(len=40) :: title              !文件 Title 字符串
    character(len=40) :: V1,V2,V3,V4,V5     !变量名字符串（X,Y,U,V,T）
    integer(kind=4), parameter :: kmax=1    !二维数据也按 3D 的 IJK 写，K=1
    character(len=40) :: zoneName           !zone 名称
    character(len=100) :: filename          !输出文件名字符串
#ifdef steadyFlow
    write(filename,'(i12.12)') restartItcOffset+itc
#endif

#ifdef unsteadyFlow
    pltFileNum = pltFileNum+1               !Tecplot 文件使用独立输出编号
    write(filename,'(i12.12)') pltFileNum
#endif

    filename = adjustl(filename)            !存储路径 pltFolderPrefix="./pltFile/buoyancyCavity000000000034.plt

    open(41,file=trim(pltFolderPrefix)//"-"//trim(filename)//'.plt', access='stream', form='unformatted')    !stream：字节流

    !---------------------------------------------
    zoneMarker= 299.0                                     !固定的，不需要修改
    eohMarker = 357.0

    !I. HEAD SECTION--------------------------------------
    !c--Magic number, Version number
    write(41) "#!TDV101"                                  !Tecplot 识别二进制 .plt 的“魔数/版本号”字符串,不需要修改

    !c--Integer value of 1
    write(41) 1                                           !Tecplot 规范里紧跟着一个整型值（通常表示字节序/文件类型版本等控制字段）,不需要修改

    Title="MyFirst"                                       !Tecplot 的二进制格式里字符串不是直接 write，而是逐字符写 ASCII 码，再以 0 结尾
    call dumpstring(title)                                !dumpstring() 就是干这个的

    !c-- Number of variables in this data file
    write(41) 5                                           !有 5 个变量：X, Y, U, V, T,如果需要修改，这个有变化需要修改

    !c-- Variable names.                                  !变量名依次写入,有变化需要修改
    V1='X'
    call dumpstring(V1)
    V2='Y'
    call dumpstring(V2)
    V3='U_nd'
    call dumpstring(V3)
    V4='V_nd'
    call dumpstring(V4)
    V5='T'
    call dumpstring(V5)

    !c-----Zones-----------------------------

    !c--------Zone marker. Value = 299.0                   !写入 float 299.0，告诉 Tecplot：一个 zone 的描述开始了。
    write(41) zoneMarker

    !--------Zone name.
    zoneName="ZONE 001"
    call dumpstring(zoneName)                              !zone 名称，显示在 Tecplot 里（比如 Zones 面板里）

    !---------Zone Color                                   !-1 表示使用默认配色
    write(41) -1

    !---------ZoneType                                     !0 通常表示 ORDERED（结构网格 IJK）
    write(41) 0

    !---------DataPacking 0=Block, 1=Point                 !0 = Block（先把整个 X 写完，再写整个 Y…）,1 = Point（每个点依次写 X,Y,U,V,T）

    write(41) 1

    !---------Specify Var Location. 0 = Do not specify, all data    !0：不指定，默认所有变量都在节点上（nodal）,1：会跟着一串变量位置列表（cell-centered/nodal）
    !---------is located at the nodes. 1 = Specify
    write(41) 0

    !---------Number of user defined face neighbor connections   !自定义邻接连接数；结构网格一般不需要，写 0
    ! (value >= 0)
    write(41) 0

    !---------IMax,JMax,KMax                          !结构网格尺寸：I=nx, J=ny, K=1, 有变化需要修改
    write(41) nx
    write(41) ny
    write(41) kmax

    !-----------1=Auxiliary name/value pair to follow          !0：没有辅助信息（Aux data）,然后写 357.0，告诉 Tecplot：头部结束，接下来是数据区描述/数据
    !-----------0=No more Auxiliar name/value pairs.
    write(41) 0
    write(41) eohMarker

    ! 再写一次 299.0：Tecplot 规范中，zone header 后的数据描述块以 zone marker 开始。
    write(41) zoneMarker

    !--------variable data format, 1=Float, 2=Double, 3=LongInt,4=ShortInt, 5=Byte, 6=Bit  !每个变量的数据格式（这里都是 float）,双精度就是2
    write(41) 2
    write(41) 2
    write(41) 2
    write(41) 2
    write(41) 2

    !--------Has variable sharing 0 = no, 1 = yes.            !0：不共享（每个变量都在这个文件里独立存）,1：共享（例如多个 zone 共享同一份 X,Y 坐标）
    write(41) 0

    !----------Zone number to share connectivity list with (-1 = no    !-1：不共享连接表（对 ordered zone 一般无连接表）
    ! sharing).
    write(41) -1

    !---------------------------------------------------------------------   !真正写数据（按 Point packing）
    do k=1,kmax
        do j=1,ny
            do i=1,nx
                write(41) real(xp(i),kind=8)
                write(41) real(yp(j),kind=8)
                write(41) real(velocityScaleCompare*u(i,j),kind=8)
                write(41) real(velocityScaleCompare*v(i,j),kind=8)
                write(41) real(T(i,j),kind=8)
            end do
        end do
    enddo
    close(41)
    !---------------------------------------------------------------------

    return
  end subroutine output_Tecplot
!===================================================================================================
! output_Tecplot 结束: 输出主场变量到 Tecplot 文件，便于后处理和可视化。
!===================================================================================================




!===================================================================================================
! 子程序: dumpstring
! 作用: 把字符串按 Tecplot 二进制格式写入文件。
! 用途: 在 output_Tecplot 和 output_Tecplot_psi_vort 中作为辅助写出工具调用。
!===================================================================================================
  subroutine dumpstring(instring)
    implicit none
    character(len=40) instring
    integer(kind=4) :: stringLength   !有效长度（去掉尾部空格）
    integer(kind=4) :: ii             !字符索引
    integer(kind=4) :: I              !字符的 ASCII 码整数

    stringLength=LEN_TRIM(instring)   !LEN_TRIM: 得到去掉尾部空格后的长度
    do ii=1,stringLength
        I=ICHAR(instring(ii:ii))      !ICHAR 把字符转成 ASCII 编码整数
        write(41) I                   !把这个整数写入文件（Tecplot 要求以整数序列写字符串）
    end do
    write(41) 0                       !最后写一个 0 作为字符串结束符

    return
  end subroutine dumpstring
!===================================================================================================
! dumpstring 结束: 把字符串按 Tecplot 二进制格式写入文件。
!===================================================================================================



!===================================================================================================
#ifdef steadyFlow
! 子程序: calNuRe
! 作用: 稳态计算结束时计算一次体平均 Nu/Re，并以累计格子步 itc 作为第一列。
! 说明: 非稳态采样由 calculate_unsteady_sample 完成，第一列改用无量纲自由落体时间 t/t_ff。
!===================================================================================================
  subroutine calNuRe()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: NuVolAvg_temp    !体平均 Nu
    real(kind=8) :: ReVolAvg_temp    !体平均 Re
    logical :: exNu, exRe
    logical, save :: first_nure_write = .true.


    !$acc wait(1)
    if(currentRunStatisticSampleCount.GE.statisticSampleCountMax) then
        write(*,*) "Error: current-run statistic samples exceed the Nu/Re array capacity"
        open(unit=00,file=trim(settingsFile),status="unknown",position="append")
        write(00,*) "Error: current-run statistic samples exceed the Nu/Re array capacity"
        close(00)
        stop
    endif
    currentRunStatisticSampleCount = currentRunStatisticSampleCount+1 ! 稳态结束时只调用一次

    if((first_nure_write).AND.(loadInitField.EQ.1)) then
        inquire(file="Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat", exist=exNu)
        inquire(file="Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat", exist=exRe)
        if((.not.exNu).OR.(.not.exRe)) then
            write(*,*) "Error: restart requested but old Nu/Re time-series files are missing."
            open(unit=00,file=trim(settingsFile),status="unknown",position="append")
            write(00,*) "Error: restart requested but old Nu/Re time-series files are missing."
            write(00,*) "Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat exists =", exNu
            write(00,*) "Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat exists =", exRe
            close(00)
            stop
        endif
    endif



    NuVolAvg_temp = 0.0d0
#ifdef SideHeatedCell
    !$acc parallel loop collapse(2) default(none) present(u,T) reduction(+:NuVolAvg_temp)
    do j = 1, ny
        do i = 1, nx
            NuVolAvg_temp = NuVolAvg_temp+u(i,j)*(T(i,j)-Tref)     !对流热通量
        enddo
    enddo
#endif

#ifdef RayleighBenardCell
    !$acc parallel loop collapse(2) default(none) present(v,T) reduction(+:NuVolAvg_temp)
    do j = 1, ny
        do i = 1, nx
            NuVolAvg_temp = NuVolAvg_temp+v(i,j)*(T(i,j)-Tref)     !对流热通量
        enddo
    enddo
#endif


    NuVolAvg(currentRunStatisticSampleCount) = &
        NuVolAvg_temp/dble(nx*ny)*lengthUnit/diffusivity+1.0d0

    if((first_nure_write).AND.(loadInitField.EQ.0)) then
        open(unit=01,file="Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat",status='replace',action='write')
        write(01,'(A)') '# itc Nu_volume'
    else
        open(unit=01,file="Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat",status='unknown',position='append',action='write')
    endif
    write(01,'(I12,1X,ES24.16E3)') restartItcOffset+itc, &
        real(NuVolAvg(currentRunStatisticSampleCount),kind=8)
    close(01)

    ReVolAvg_temp = 0.0d0
    !$acc parallel loop collapse(2) default(none) present(u,v) reduction(+:ReVolAvg_temp)
    do j = 1, ny
        do i = 1, nx
            ReVolAvg_temp = ReVolAvg_temp+(u(i,j)*u(i,j)+v(i,j)*v(i,j))
        enddo
    enddo
    ReVolAvg(currentRunStatisticSampleCount) = &
        dsqrt(ReVolAvg_temp/dble(nx*ny))*lengthUnit/viscosity


    if((first_nure_write).AND.(loadInitField.EQ.0)) then
        open(unit=02,file="Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat",status='replace',action='write')
        write(02,'(A)') '# itc Re_spatial_rms'
    else
        open(unit=02,file="Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat",status='unknown',position='append',action='write')
    endif
    write(02,'(I12,1X,ES24.16E3)') restartItcOffset+itc, &
        real(ReVolAvg(currentRunStatisticSampleCount),kind=8)
    close(02)
    first_nure_write = .false.
    write(*,'(a,1x,ES24.16E3)') "NuVolAvg =", real(NuVolAvg(currentRunStatisticSampleCount),kind=8)
    write(*,'(a,1x,ES24.16E3)') "ReVolAvg =", real(ReVolAvg(currentRunStatisticSampleCount),kind=8)
    return
  end subroutine calNuRe
!===================================================================================================
! calNuRe 结束: 稳态计算结束时输出一次体平均 Nu 和空间 RMS Re。
!===================================================================================================
#endif


#ifdef unsteadyFlow
!===================================================================================================
! 子程序: restore_dissipation_statistics
! 作用: 续算时从完整耗散历史和温度剖面历史恢复统计窗口的全部在线累计量。
! 说明: 耗散历史恢复标量统计；二进制温度剖面历史恢复每个高度上的 <T>_x 和 <T^2>_x。
!       两套历史必须都与重启时刻一致，防止从较旧流场续算时重复累计较新的后处理样本。
!===================================================================================================
  subroutine restore_dissipation_statistics()
    use commondata
    implicit none
    integer(kind=4) :: ios, dissipationUnit, profileUnit, k
    real(kind=8) :: currentStatisticSampleTf, nuVolAvgHistory, reVolRmsHistory
    real(kind=8) :: epsKineticVolAvgHistory, epsThermalVolAvgHistory
    real(kind=8) :: densityFluctuationVolRmsHistory, velocityDivergenceVolRmsHistory
    real(kind=8) :: maxMachHistory, minTemperatureHistory, maxTemperatureHistory
    real(kind=8) :: lastDissipationTf, profileTf, lastProfileTf, tfTolerance
    real(kind=8) :: temperatureXAvgProfileSample(ny), temperatureSquaredXAvgProfileSample(ny)
    character(len=512) :: historyHeader
    logical :: dissipationHistoryExists, profileHistoryExists

    ! latest.meta 已经保存了该重启场对应的累计统计样本数；为零时没有历史需要恢复。
    if(cumulativeStatisticSampleCount.LE.0) return
    tfTolerance = 1.0d-10*max(1.0d0,abs(restartTfOffset),abs(statisticSampleInterval))

    inquire(file=trim(dissipationHistoryFile),exist=dissipationHistoryExists)
    if(.not.dissipationHistoryExists) then
        write(*,*) 'Error: dissipation history is missing for restart statistics restoration.'
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,*) 'Error: missing ', trim(dissipationHistoryFile)
        close(00)
        error stop 1
    endif

    open(newunit=dissipationUnit,file=trim(dissipationHistoryFile),status='old', &
        action='read',form='formatted')
    read(dissipationUnit,'(A)',iostat=ios) historyHeader
    ! 重命名只改变表头，不改变十列数据顺序；续算时同时兼容旧表头和当前表头。
    if((ios.NE.0).OR.((index(historyHeader, &
        'rho_fluctuation_vol_rms div_u_vol_rms maxMach minT maxT').EQ.0).AND. &
        (index(historyHeader,'rho_fluctuation_rms div_u_rms maxMach minT maxT').EQ.0))) then
        close(dissipationUnit)
        write(*,*) 'Error: dissipation history format is incompatible with restart restoration.'
        error stop 1
    endif

    lastDissipationTf = -huge(1.0d0)
    do k = 1, cumulativeStatisticSampleCount   !读取之前的统计数据
        read(dissipationUnit,*,iostat=ios) currentStatisticSampleTf, nuVolAvgHistory, reVolRmsHistory, &
            epsKineticVolAvgHistory, epsThermalVolAvgHistory, densityFluctuationVolRmsHistory, &
            velocityDivergenceVolRmsHistory, maxMachHistory, minTemperatureHistory, maxTemperatureHistory
        if(ios.NE.0) then
            close(dissipationUnit)
            write(*,*) 'Error: dissipation history is shorter than cumulativeStatisticSampleCount or invalid.'
            error stop 1
        endif
        if((k.GT.1).AND.(currentStatisticSampleTf.LE.lastDissipationTf+tfTolerance)) then
            close(dissipationUnit)
            write(*,*) 'Error: dissipation history times are not strictly increasing.'
            error stop 1
        endif
        if(currentStatisticSampleTf.GT.restartTfOffset+tfTolerance) then
            close(dissipationUnit)
            write(*,*) 'Error: dissipation history is newer than the selected restart field.'
            error stop 1
        endif
        lastDissipationTf = currentStatisticSampleTf  !恢复之前的统计量
        sumNuVolAvg = sumNuVolAvg+nuVolAvgHistory
        sumSpeedSquaredVolAvg = sumSpeedSquaredVolAvg+ &
            (reVolRmsHistory*viscosity/lengthUnit)*(reVolRmsHistory*viscosity/lengthUnit)
        sumEpsKineticVolAvg = sumEpsKineticVolAvg+epsKineticVolAvgHistory
        sumEpsThermalVolAvg = sumEpsThermalVolAvg+epsThermalVolAvgHistory
        sumDensityFluctuationSquaredVolAvg = sumDensityFluctuationSquaredVolAvg+ &
            densityFluctuationVolRmsHistory*densityFluctuationVolRmsHistory
        sumVelocityDivergenceSquaredVolAvg = sumVelocityDivergenceSquaredVolAvg+ &
            velocityDivergenceVolRmsHistory*velocityDivergenceVolRmsHistory
        maxStatisticCFL = max(maxStatisticCFL,maxMachHistory*dsqrt(cs2))
    enddo
    read(dissipationUnit,*,iostat=ios) currentStatisticSampleTf
    close(dissipationUnit)
    if(ios.GE.0) then
        write(*,*) 'Error: dissipation history contains samples newer than the restart metadata.'
        error stop 1
    endif

    inquire(file=trim(temperatureProfileHistoryFile),exist=profileHistoryExists)
    if(.not.profileHistoryExists) then
        write(*,*) 'Error: temperature-profile history is missing for restart statistics restoration.'
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,*) 'Error: missing ', trim(temperatureProfileHistoryFile)
        close(00)
        error stop 1
    endif

    lastProfileTf = -huge(1.0d0)
    open(newunit=profileUnit,file=trim(temperatureProfileHistoryFile),status='old', &
        action='read',form='unformatted',access='stream')
    do k = 1, cumulativeStatisticSampleCount            !读取之前的统计数据
        read(profileUnit,iostat=ios) profileTf
        if(ios.NE.0) then
            close(profileUnit)
            write(*,*) 'Error: temperature-profile history is shorter than cumulativeStatisticSampleCount or invalid.'
            error stop 1
        endif
        read(profileUnit,iostat=ios) temperatureXAvgProfileSample, temperatureSquaredXAvgProfileSample
        if(ios.NE.0) then
            close(profileUnit)
            write(*,*) 'Error: incomplete record in temperature-profile history.'
            error stop 1
        endif
        if((k.GT.1).AND.(profileTf.LE.lastProfileTf+tfTolerance)) then
            close(profileUnit)
            write(*,*) 'Error: temperature-profile history times are not strictly increasing.'
            error stop 1
        endif
        if(profileTf.GT.restartTfOffset+tfTolerance) then
            close(profileUnit)
            write(*,*) 'Error: temperature-profile history is newer than the selected restart field.'
            error stop 1
        endif
        lastProfileTf = profileTf
        sumTemperatureXAvgProfile = sumTemperatureXAvgProfile+temperatureXAvgProfileSample
        sumTemperatureSquaredXAvgProfile = &
            sumTemperatureSquaredXAvgProfile+temperatureSquaredXAvgProfileSample
    enddo
    read(profileUnit,iostat=ios) profileTf
    close(profileUnit)
    if(ios.GE.0) then
        write(*,*) 'Error: temperature-profile history contains samples newer than the restart metadata.'
        error stop 1
    endif
    if(abs(lastProfileTf-lastDissipationTf).GT.tfTolerance) then
        write(*,*) 'Error: dissipation and temperature-profile restart statistics do not match.'
        error stop 1
    endif

    open(unit=00,file=trim(settingsFile),status='unknown',position='append')
    write(00,*) 'Restored dissipation/temperature-profile samples =', cumulativeStatisticSampleCount
    write(00,*) 'Restored statistics through time_tf =', lastDissipationTf
    close(00)
  end subroutine restore_dissipation_statistics
!===================================================================================================


!===================================================================================================
! 子程序: calculate_unsteady_sample
! 作用: 一次全场遍历同时计算非稳态瞬时 Nu、Re、两类耗散及可压缩性指标。
! 调用: 只在配置的统计窗口内按 statisticSampleInterval 调用；快照、plt 和 reload 不调用本程序。
! 时间轴: 给定 t/t_ff 后，用 itc=nint((t/t_ff)*timeUnit) 选择最近格子步；写文件仍以给定 t/t_ff 为横坐标。
! 说明: 此处只记录原始瞬时量，不在每个采样点反算耗散 Nu 或检查精确关系。
!       Xu/Zhang 所需的耗散 Nu 和精确关系比值在统计窗口结束后由空间-时间平均量统一计算。
!       速度耗散使用 LBM-CDE 局部应变率；温度梯度使用二阶中心/单边有限差分。
!===================================================================================================
  subroutine calculate_unsteady_sample()
    use commondata
    implicit none
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: gradX, gradY, speedSquared, eu, uu, feq
    real(kind=8) :: neqxx, neqxy, neqyy
    ! xxxVol 是当前采样时刻对全部格点的空间求和，尚未除以 nx*ny。
    real(kind=8) :: nuConvectiveFluxVol, speedSquaredVol
    real(kind=8) :: epsKineticVol, epsThermalVol
    real(kind=8) :: densityFluctuationSquaredVol
    real(kind=8) :: velocityDivergenceSquaredVol
    real(kind=8) :: maxMachWork, minTemperatureWork, maxTemperatureWork
    real(kind=8) :: currentStatisticSampleTf
    logical :: exNu, exRe
    logical, save :: firstUnsteadyHistoryWrite = .true.
    real(kind=8), parameter :: muShear=rho0*viscosity
    real(kind=8), parameter :: muBulk=rho0*bulkViscosity
    real(kind=8), parameter :: denomDiag=2.0d0*muShear+rho0*cs2
    real(kind=8), parameter :: denomShear=4.0d0*muShear+2.0d0*rho0*cs2
    real(kind=8), parameter :: denomDiv=2.0d0*muBulk+rho0*cs2
    real(kind=8), parameter :: coeffTrace=(2.0d0*muShear-2.0d0*muBulk)/(2.0d0*denomDiag)

    !$acc wait(1)
    if(currentRunStatisticSampleCount.GE.statisticSampleCountMax) then
        write(*,*) "Error: current-run statistic samples exceed the Nu/Re array capacity"
        open(unit=00,file=trim(settingsFile),status="unknown",position="append")
        write(00,*) "Error: current-run statistic samples exceed the Nu/Re array capacity"
        close(00)
        stop
    endif
    currentRunStatisticSampleCount = currentRunStatisticSampleCount+1
    currentStatisticSampleTf = unsteadyAverageStartTf+ &
        real(cumulativeStatisticSampleCount,kind=8)*statisticSampleInterval

    ! 只有重启时刻已经越过统计窗口首个样本，续算才必须接在旧 Nu/Re 历史之后。
    ! 若在统计窗口开始前重启，cumulativeStatisticSampleCount=0，首个统计样本会新建历史文件。
    if((firstUnsteadyHistoryWrite).AND.(cumulativeStatisticSampleCount.GT.0)) then
        inquire(file="Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat", exist=exNu)
        inquire(file="Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat", exist=exRe)
        if((.not.exNu).OR.(.not.exRe)) then
            write(*,*) "Error: restart requested but old Nu/Re time-series files are missing."
            open(unit=00,file=trim(settingsFile),status="unknown",position="append")
            write(00,*) "Error: restart requested but old Nu/Re time-series files are missing."
            write(00,*) "Nu history exists =", exNu
            write(00,*) "Re history exists =", exRe
            close(00)
            stop
        endif
    endif

    nuConvectiveFluxVol = 0.0d0
    speedSquaredVol = 0.0d0
    epsKineticVol = 0.0d0
    epsThermalVol = 0.0d0
    densityFluctuationSquaredVol = 0.0d0
    velocityDivergenceSquaredVol = 0.0d0
    maxMachWork = 0.0d0
    minTemperatureWork = huge(1.0d0)
    maxTemperatureWork = -huge(1.0d0)

    !$acc parallel loop gang vector collapse(2) default(none) &
    !$acc& present(f,u,v,T,rho,Fx,Fy,Sxx,Sxy,Syy,Sdiv,ex,ey,omega) &
    !$acc& reduction(+:nuConvectiveFluxVol,speedSquaredVol) &
    !$acc& reduction(+:epsKineticVol,epsThermalVol) &
    !$acc& reduction(+:densityFluctuationSquaredVol,velocityDivergenceSquaredVol) &
    !$acc& reduction(max:maxMachWork,maxTemperatureWork) &
    !$acc& reduction(min:minTemperatureWork) &
    !$acc& private(alpha,gradX,gradY,speedSquared,eu,uu,feq,neqxx,neqxy,neqyy)
    do j = 1, ny
        do i = 1, nx
            ! 温度在本步 macroT 后已经更新，因此这里同步刷新采样时刻的力和应变率。
            Fx(i,j) = 0.0d0
            Fy(i,j) = rho0*gBeta*(T(i,j)-Tref)
#ifdef SideHeatedHa
            Fx(i,j) = rho0*B2sigemarho*(v(i,j)*sin(phi)*cos(phi)-u(i,j)*sin(phi)*sin(phi))
            Fy(i,j) = Fy(i,j)+rho0*B2sigemarho*(u(i,j)*sin(phi)*cos(phi)-v(i,j)*cos(phi)*cos(phi))
#endif

            neqxx = 0.0d0
            neqxy = 0.0d0
            neqyy = 0.0d0
            uu = u(i,j)*u(i,j)+v(i,j)*v(i,j)
            do alpha = 0, 8
                eu = dble(ex(alpha))*u(i,j)+dble(ey(alpha))*v(i,j)
                feq = omega(alpha)*((rho(i,j)-rho0)+rho0*(eu/cs2 + &
                    0.5d0*eu*eu/(cs2*cs2)-0.5d0*uu/cs2))
                neqxx = neqxx+dble(ex(alpha)*ex(alpha))*(f(i,j,alpha)-feq)
                neqxy = neqxy+dble(ex(alpha)*ey(alpha))*(f(i,j,alpha)-feq)
                neqyy = neqyy+dble(ey(alpha)*ey(alpha))*(f(i,j,alpha)-feq)
            enddo
            Sdiv(i,j) = -(neqxx+neqyy+u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))/denomDiv
            Sxx(i,j) = -(neqxx+u(i,j)*Fx(i,j))/denomDiag+coeffTrace*Sdiv(i,j)
            Syy(i,j) = -(neqyy+v(i,j)*Fy(i,j))/denomDiag+coeffTrace*Sdiv(i,j)
            Sxy(i,j) = -(2.0d0*neqxy+u(i,j)*Fy(i,j)+v(i,j)*Fx(i,j))/denomShear

            if(i.EQ.1) then
                gradX = 0.5d0*(-3.0d0*T(1,j)+4.0d0*T(2,j)-T(3,j))
            elseif(i.EQ.nx) then
                gradX = 0.5d0*(3.0d0*T(nx,j)-4.0d0*T(nx-1,j)+T(nx-2,j))
            else
                gradX = 0.5d0*(T(i+1,j)-T(i-1,j))
            endif
            if(j.EQ.1) then
                gradY = 0.5d0*(-3.0d0*T(i,1)+4.0d0*T(i,2)-T(i,3))
            elseif(j.EQ.ny) then
                gradY = 0.5d0*(3.0d0*T(i,ny)-4.0d0*T(i,ny-1)+T(i,ny-2))
            else
                gradY = 0.5d0*(T(i,j+1)-T(i,j-1))
            endif

            speedSquared = u(i,j)*u(i,j)+v(i,j)*v(i,j)
#ifdef SideHeatedCell
            ! 侧壁差温方腔的主热流方向为 x，因此使用 u(T-Tref)。
            nuConvectiveFluxVol = nuConvectiveFluxVol+u(i,j)*(T(i,j)-Tref)
#endif
#ifdef RayleighBenardCell
            ! Rayleigh-Benard 的主热流方向为 y，因此使用 v(T-Tref)。
            nuConvectiveFluxVol = nuConvectiveFluxVol+v(i,j)*(T(i,j)-Tref)
#endif
            speedSquaredVol = speedSquaredVol+speedSquared
            epsKineticVol = epsKineticVol+2.0d0*viscosity* &
                (Sxx(i,j)*Sxx(i,j)+Syy(i,j)*Syy(i,j)+2.0d0*Sxy(i,j)*Sxy(i,j))
            epsThermalVol = epsThermalVol+ &
                diffusivity*(gradX*gradX+gradY*gradY)
            ! Xu 可压缩性指标的空间平方平均；最后在全场归一化并开平方。
            densityFluctuationSquaredVol = densityFluctuationSquaredVol+ &
                (rho(i,j)-rho0)*(rho(i,j)-rho0)
            velocityDivergenceSquaredVol = velocityDivergenceSquaredVol+ &
                Sdiv(i,j)*Sdiv(i,j)
            maxMachWork = max(maxMachWork,dsqrt(speedSquared/cs2))
            minTemperatureWork = min(minTemperatureWork,T(i,j))
            maxTemperatureWork = max(maxTemperatureWork,T(i,j))
        enddo
    enddo

    ! 前面是累加，现在是空间平均。
    NuVolAvg(currentRunStatisticSampleCount) = &
        1.0d0+nuConvectiveFluxVol/dble(nx*ny)*lengthUnit/diffusivity
    speedSquaredVolAvg = speedSquaredVol/dble(nx*ny)
    ReVolAvg(currentRunStatisticSampleCount) = &
        dsqrt(max(0.0d0,speedSquaredVolAvg))*lengthUnit/viscosity
    epsKineticVolAvg = epsKineticVol/dble(nx*ny)
    epsThermalVolAvg = epsThermalVol/dble(nx*ny)
    densityFluctuationSquaredVolAvg = densityFluctuationSquaredVol/dble(nx*ny)
    velocityDivergenceSquaredVolAvg = velocityDivergenceSquaredVol/dble(nx*ny)
    maxMachLocal = maxMachWork
    minTemperature = minTemperatureWork
    maxTemperature = maxTemperatureWork

    if((firstUnsteadyHistoryWrite).AND.(cumulativeStatisticSampleCount.EQ.0)) then
        open(unit=01,file="Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat",status='replace',action='write')
        open(unit=02,file="Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat",status='replace',action='write')
    else
        open(unit=01,file="Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat",status='unknown',position='append',action='write')
        open(unit=02,file="Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat",status='unknown',position='append',action='write')
    endif
    write(01,'(ES24.16E3,1X,ES24.16E3)') currentStatisticSampleTf, &
        NuVolAvg(currentRunStatisticSampleCount)
    write(02,'(ES24.16E3,1X,ES24.16E3)') currentStatisticSampleTf, &
        ReVolAvg(currentRunStatisticSampleCount)
    close(01)
    close(02)

    if((firstUnsteadyHistoryWrite).AND.(cumulativeStatisticSampleCount.EQ.0)) then
        open(unit=12,file=trim(dissipationHistoryFile),status='replace',action='write')
        write(12,'(A)') '# time_over_tff Nu_vol_avg Re_vol_rms eps_u_vol_avg eps_T_vol_avg ' // &
            'rho_fluctuation_vol_rms div_u_vol_rms maxMach minT maxT'
    else
        open(unit=12,file=trim(dissipationHistoryFile),status='unknown',position='append',action='write')
    endif
    ! 历史文件需要瞬时空间 RMS；这里直接由瞬时空间均方量开平方，不再保存重复中间变量。
    write(12,'(10(ES24.16E3,1X))') currentStatisticSampleTf, &
        NuVolAvg(currentRunStatisticSampleCount), ReVolAvg(currentRunStatisticSampleCount), &
        epsKineticVolAvg, epsThermalVolAvg, &
        dsqrt(max(0.0d0,densityFluctuationSquaredVolAvg)), &
        dsqrt(max(0.0d0,velocityDivergenceSquaredVolAvg)), &
        maxMachLocal, minTemperature, maxTemperature
    close(12)
    firstUnsteadyHistoryWrite = .false.
    write(*,'(A,1X,ES24.16E3)') 'Nu_volume_instantaneous =', &
        NuVolAvg(currentRunStatisticSampleCount)
    write(*,'(A,1X,ES24.16E3)') 'Re_instantaneous_rms =', &
        ReVolAvg(currentRunStatisticSampleCount)
  end subroutine calculate_unsteady_sample
!===================================================================================================


!===================================================================================================
! 子程序: accumulate_dissipation_statistics
! 作用: 对统计窗口内的 Nu、速度平方、两类耗散、Xu 可压缩性指标和温度剖面做在线累加。
! 调用: calculate_unsteady_sample 已经完成当前采样时刻的空间求和与体平均；
!       主程序只在 unsteadyAverageStartTf 至 unsteadyAverageEndTf 的统计窗口内调用本程序。
! 说明: 速度平方必须保留，因为文献 Re 的定义是 sqrt(<u^2+v^2>_V,t)*H/nu。
!       前后半段 Nu/Re 统一由 output_unsteady_NuRe_postprocess 重读完整历史后计算，
!       本程序不再维护一套重复的前后半段累计量。
!       温度的一、二阶矩继续保留，用于以后分析 T_rms(y) 峰值定义的边界层厚度；
!       当前 Table 1 的 N_BL 在 write_dissipation_statistics 中由 H/(2*Nu) 直接估计。
!       每个统计采样的水平温度剖面同时追加到二进制历史，供后续进程恢复完整累计量。
!       本程序只做等时间间隔样本的时间求和；除以样本数及精确关系检查在
!       write_dissipation_statistics 中统一完成。
!===================================================================================================
  subroutine accumulate_dissipation_statistics()
    use commondata
    implicit none
    integer(kind=4) :: i, j, profileUnit ! i/j 为温度剖面空间循环；profileUnit 为二进制历史文件单元号
    real(kind=8) :: currentStatisticSampleTf ! 当前统计样本对应的自由落体时间 t/t_ff
    real(kind=8) :: temperatureXSum, temperatureSquaredXSum ! 固定高度上 T 和 T^2 的水平格点求和
    real(kind=8) :: temperatureXAvgProfileSample(ny)        ! 当前时刻的水平平均温度 <T>_x(y)
    real(kind=8) :: temperatureSquaredXAvgProfileSample(ny) ! 当前时刻的水平平均温度平方 <T^2>_x(y)
    logical, save :: firstProfileHistoryWrite = .true.      ! 标记本次进程是否第一次写温度剖面历史

    currentStatisticSampleTf = unsteadyAverageStartTf+ &
        real(cumulativeStatisticSampleCount,kind=8)*statisticSampleInterval

    ! 对已经完成空间平均的瞬时量做时间求和。
    ! 固定采样间隔下，最终除以 cumulativeStatisticSampleCount 即为离散时间平均。
    sumNuVolAvg = sumNuVolAvg+NuVolAvg(currentRunStatisticSampleCount)
    ! 文献 Re 需要先对 <u^2+v^2>_V 做时间平均，最后统一开平方，不能直接平均瞬时 Re。
    sumSpeedSquaredVolAvg = sumSpeedSquaredVolAvg+speedSquaredVolAvg
    ! 两类耗散先累计瞬时体平均量，统计结束后得到 <epsilon_u>_V,t 和 <epsilon_T>_V,t。
    sumEpsKineticVolAvg = sumEpsKineticVolAvg+epsKineticVolAvg
    sumEpsThermalVolAvg = sumEpsThermalVolAvg+epsThermalVolAvg
    ! Xu 的弱可压缩性指标采用空间-时间 RMS，因此累计均方量而不是瞬时 RMS。
    sumDensityFluctuationSquaredVolAvg = &
        sumDensityFluctuationSquaredVolAvg+densityFluctuationSquaredVolAvg
    sumVelocityDivergenceSquaredVolAvg = &
        sumVelocityDivergenceSquaredVolAvg+velocityDivergenceSquaredVolAvg
    ! maxMachLocal*cs 恢复为最大格子速度；在 dx=dt=1 时也就是格子 CFL。
    maxStatisticCFL = max(maxStatisticCFL,maxMachLocal*dsqrt(cs2))

    ! 对每个高度 y_j 做水平平均：
    ! temperatureXAvgProfileSample(j)        = <T>_x(y_j,t)
    ! temperatureSquaredXAvgProfileSample(j) = <T^2>_x(y_j,t)
    ! 两者再做时间求和，最终使用 T_rms(y)=sqrt(<T^2>_x,t-<T>_x,t^2)。
    do j = 1, ny
        temperatureXSum = 0.0d0
        temperatureSquaredXSum = 0.0d0
        do i = 1, nx
            temperatureXSum = temperatureXSum+T(i,j)
            temperatureSquaredXSum = temperatureSquaredXSum+T(i,j)*T(i,j)
        enddo
        temperatureXAvgProfileSample(j) = temperatureXSum/dble(nx)
        temperatureSquaredXAvgProfileSample(j) = temperatureSquaredXSum/dble(nx)
        sumTemperatureXAvgProfile(j) = &
            sumTemperatureXAvgProfile(j)+temperatureXAvgProfileSample(j)
        sumTemperatureSquaredXAvgProfile(j) = &
            sumTemperatureSquaredXAvgProfile(j)+temperatureSquaredXAvgProfileSample(j)
    enddo

    ! 温度剖面不在标量耗散历史中，因此单独保存为无格式流二进制文件。
    ! 没有旧统计样本时，首个统计样本覆盖旧文件；恢复过旧统计样本时则继续追加。
    ! 每个样本依次保存 t/t_ff、<T>_x(:) 和 <T^2>_x(:)，供续算恢复完整时间累计量。
    if(firstProfileHistoryWrite.AND.(cumulativeStatisticSampleCount.EQ.0)) then
        open(newunit=profileUnit,file=trim(temperatureProfileHistoryFile),status='replace', &
            action='write',form='unformatted',access='stream')
    else
        open(newunit=profileUnit,file=trim(temperatureProfileHistoryFile),status='unknown', &
            action='write',form='unformatted',access='stream',position='append')
    endif
    write(profileUnit) currentStatisticSampleTf
    write(profileUnit) temperatureXAvgProfileSample, temperatureSquaredXAvgProfileSample
    close(profileUnit)
    ! save 变量在本次进程后续调用中保持为 .false.，避免再次覆盖刚写入的历史文件。
    firstProfileHistoryWrite = .false.
    cumulativeStatisticSampleCount = cumulativeStatisticSampleCount+1 ! 当前样本已完整累计
  end subroutine accumulate_dissipation_statistics
!===================================================================================================


!===================================================================================================
! 子程序: write_dissipation_statistics
! 作用: 在统计窗口结束后，用空间-时间平均量统一生成 Xu 与 Zhang 两套文献对照结果。
! 调用: 非稳态计算结束时调用一次；本程序不推进流场，也不重新遍历二维瞬时流场。
! 输入: accumulate_dissipation_statistics 在线累计的标量和温度剖面；续算时这些累计量已经
!       由 restore_dissipation_statistics 恢复，所以最终结果覆盖完整统计窗口。
! Xu 输出: 对流通量 Nu、壁面 Nu、两种耗散反算 Nu、后三者相对 Nu_vol 的差异，
!          以及 RMS Re、RMS 密度波动和 RMS 速度散度。
! Zhang 输出: 直接计算耗散、Nu 精确关系给出的耗散、两者比值，以及 N_BL 和分辨率指标。
!             N_BL 采用 H/(2*Nu) 的整体估计，同时保留温度 RMS 峰值边界层诊断。
! 文件: statisticsFile 保存汇总量；temperatureRmsProfileFile 保存 z/H-theta_rms(z) 剖面；
!       settingsFile 只追加一份便于运行结束后快速查看的摘要。
! 说明: 所有文献比较都先完成空间-时间平均，再计算比值或反算 Nu；不平均瞬时比值。
!===================================================================================================
  subroutine write_dissipation_statistics()
    use commondata
    implicit none
    integer(kind=4) :: j                         ! 壁法向网格索引
    integer(kind=4) :: lowerPeakIndex, upperPeakIndex ! 下、上半区 theta_rms 峰值索引
    integer(kind=4) :: profileUnit               ! 温度 RMS 剖面输出文件单元号
    integer(kind=4) :: thermalBLGridPointsGlobalEstimate ! H/(2*Nu) 对应的 N_BL
    integer(kind=4) :: thermalBLGridPointsRms    ! theta_rms 峰值厚度对应的 N_BL
    real(kind=8) :: invN                         ! 统计样本数的倒数，用于时间平均
    real(kind=8) :: nuVolTimeAvg                 ! 对流热通量 Nu 的时间平均
    real(kind=8) :: speedSquaredVolTimeAvg       ! <u^2+v^2>_V,t
    real(kind=8) :: reVolTimeRms                 ! sqrt(<u^2+v^2>_V,t)*H/nu
    real(kind=8) :: epsKineticVolTimeAvg, epsThermalVolTimeAvg ! 两类空间-时间平均耗散率
    real(kind=8) :: densityFluctuationSquaredVolTimeAvg ! <(rho-rho0)^2>_V,t
    real(kind=8) :: velocityDivergenceSquaredVolTimeAvg ! <(div u)^2>_V,t
    real(kind=8) :: densityFluctuationVolTimeRms, velocityDivergenceVolTimeRms ! 开平方后的 RMS
    real(kind=8) :: wallNormalSpacing            ! 无量纲壁法向网格间距 Delta_g/H
    real(kind=8) :: nuBottomWallTimeAvg, nuTopWallTimeAvg, nuWallTimeAvg ! 上下壁及两壁平均 Nu
    real(kind=8) :: nuFromKineticDissipation, nuFromThermalDissipation ! 两个耗散反算 Nu
    real(kind=8) :: kineticDissipationExact, thermalDissipationExact ! Nu 精确关系给出的耗散
    real(kind=8) :: kineticDissipationRatio, thermalDissipationRatio ! 直接耗散/精确关系耗散
    real(kind=8) :: wallNuDifferencePercent      ! Nu_wall 相对 Nu_vol 的百分差
    real(kind=8) :: kineticNuDifferencePercent   ! Nu_kin 相对 Nu_vol 的百分差
    real(kind=8) :: thermalNuDifferencePercent   ! Nu_th 相对 Nu_vol 的百分差
    real(kind=8) :: rmsTemperature               ! 当前高度的 theta_rms
    real(kind=8) :: lowerRmsMaximum, upperRmsMaximum ! 下、上半区 theta_rms 最大值
    real(kind=8) :: thermalBLThicknessGlobalEstimate ! delta_theta/H=1/(2*Nu)
    real(kind=8) :: thermalBLThicknessRms        ! 上下 theta_rms 峰值距离的平均值
    real(kind=8) :: thermalBLThicknessRelativeDifferencePercent ! 两种厚度的百分差
    real(kind=8) :: zOverH                       ! 流体节点中心的无量纲高度
    real(kind=8) :: etaOverH, etaBOverH           ! Kolmogorov、Batchelor 长度相对 H
    real(kind=8) :: gridOverEta, gridOverEtaB, timeStepOverEta ! Zhang 表一分辨率指标

    !-----------------------------------------------------------------------------------------------
    ! 第一部分：建立最终统计文件，并检查统计窗口内是否确实存在样本。
    ! 即使没有样本也输出说明，避免留下一个看似正常但内容无效的结果文件。
    open(unit=22,file=trim(statisticsFile),status='replace',action='write')
    write(22,'(A)') '# Xu/Zhang benchmark statistics over the configured window, including restored restart history'
    write(22,'(A)') '# unsteady time axis: t/t_ff; configured window times are also mapped to lattice itc'
    write(22,'(A,I0)') '# cumulative_statistics_samples ', cumulativeStatisticSampleCount
    write(22,'(A,3(1X,I0))') '# configured_start_mid_end_sample_itc', &
        max(1,nint(unsteadyAverageStartTf*timeUnit)), &
        max(1,nint(unsteadyAverageMidTf*timeUnit)), &
        max(1,nint(unsteadyAverageEndTf*timeUnit))
    if(cumulativeStatisticSampleCount.LE.0) then
        write(22,'(A)') '# no samples: run ended before the averaging window'
        close(22)
        return
    endif

    !-----------------------------------------------------------------------------------------------
    ! 第二部分：把 accumulate_dissipation_statistics 保存的时间求和量除以样本数。
    ! 采样间隔固定，因此算术平均就是当前统计窗口上的离散时间平均。
    invN = 1.0d0/dble(cumulativeStatisticSampleCount)
    nuVolTimeAvg = sumNuVolAvg*invN
    speedSquaredVolTimeAvg = sumSpeedSquaredVolAvg*invN
    ! Xu 与 Zhang 均使用 Re=sqrt(<u^2+v^2>_V,t)*H/nu，而不是 <Re_instant>_t。
    reVolTimeRms = dsqrt(max(0.0d0,speedSquaredVolTimeAvg))*lengthUnit/viscosity
    epsKineticVolTimeAvg = sumEpsKineticVolAvg*invN
    epsThermalVolTimeAvg = sumEpsThermalVolAvg*invN
    ! Xu 的两个整体可压缩性指标采用空间-时间 RMS，先平均平方再统一开平方。
    densityFluctuationSquaredVolTimeAvg = sumDensityFluctuationSquaredVolAvg*invN
    velocityDivergenceSquaredVolTimeAvg = sumVelocityDivergenceSquaredVolAvg*invN
    densityFluctuationVolTimeRms = dsqrt(max(0.0d0,densityFluctuationSquaredVolTimeAvg))
    velocityDivergenceVolTimeRms = dsqrt(max(0.0d0,velocityDivergenceSquaredVolTimeAvg))

    !-----------------------------------------------------------------------------------------------
    ! 第三部分：先初始化只对 RB 对流有意义的壁面 Nu、耗散 Nu 和精确关系结果。
    ! 侧壁加热分支不会进入下面的 RayleighBenardCell 计算，因而保持这些安全初值。
    wallNormalSpacing = 1.0d0/lengthUnit
    nuBottomWallTimeAvg = 0.0d0
    nuTopWallTimeAvg = 0.0d0
    nuWallTimeAvg = 0.0d0
    nuFromKineticDissipation = 0.0d0
    nuFromThermalDissipation = 0.0d0
    kineticDissipationExact = 0.0d0
    thermalDissipationExact = 0.0d0
    kineticDissipationRatio = 0.0d0
    thermalDissipationRatio = 0.0d0
    wallNuDifferencePercent = 0.0d0
    kineticNuDifferencePercent = 0.0d0
    thermalNuDifferencePercent = 0.0d0
#ifdef RayleighBenardCell
    !-----------------------------------------------------------------------------------------------
    ! 第四部分（Xu 表四）：从平均壁面梯度和两类平均耗散率得到三套独立 Nu，
    ! 并统一以体积对流热通量得到的 nuVolTimeAvg 为基准计算百分差。
    ! Xu 表 4 的 Nu_wall：先对上下壁面及时间平均，再取两壁平均。
    ! 恒温壁位于首/末流体节点外半个网格处，因此沿用 RBcalc_Nu_wall_avg 的
    ! 二阶半格点壁面梯度；温度梯度是线性的，可直接使用已经累计的 <T>_x,t 剖面。
    nuBottomWallTimeAvg = (8.0d0*Thot-9.0d0*sumTemperatureXAvgProfile(1)*invN+ &
        sumTemperatureXAvgProfile(2)*invN)/(3.0d0*wallNormalSpacing*deltaT)
    nuTopWallTimeAvg = (-8.0d0*Tcold+9.0d0*sumTemperatureXAvgProfile(ny)*invN- &
        sumTemperatureXAvgProfile(ny-1)*invN)/(3.0d0*wallNormalSpacing*deltaT)
    nuWallTimeAvg = 0.5d0*(nuBottomWallTimeAvg+nuTopWallTimeAvg)

    ! Xu 表 4 的 Nu_kin 和 Nu_th：把空间-时间平均耗散率反算为两个等效 Nusselt 数。
    nuFromKineticDissipation = 1.0d0+ &
        epsKineticVolTimeAvg*lengthUnit**4*Prandtl**2/(viscosity**3*Rayleigh)
    nuFromThermalDissipation = &
        epsThermalVolTimeAvg*lengthUnit**2/(diffusivity*deltaT**2)
    if(abs(nuVolTimeAvg).GT.1.0d-12) then
        wallNuDifferencePercent = &
            100.0d0*abs(nuWallTimeAvg-nuVolTimeAvg)/abs(nuVolTimeAvg)
        kineticNuDifferencePercent = &
            100.0d0*abs(nuFromKineticDissipation-nuVolTimeAvg)/abs(nuVolTimeAvg)
        thermalNuDifferencePercent = &
            100.0d0*abs(nuFromThermalDissipation-nuVolTimeAvg)/abs(nuVolTimeAvg)
    endif

    ! Zhang 风格：直接计算“应变率/温度梯度耗散”与“Nu 精确关系耗散”的比值。
    ! 分子是代码直接得到的空间-时间平均耗散；分母由同一个 nuVolTimeAvg 代入精确关系得到。
    kineticDissipationExact = viscosity**3/lengthUnit**4* &
        (nuVolTimeAvg-1.0d0)*Rayleigh/(Prandtl**2)
    thermalDissipationExact = diffusivity*deltaT**2/lengthUnit**2*nuVolTimeAvg
    if(abs(kineticDissipationExact).GT.1.0d-30) &
        kineticDissipationRatio = epsKineticVolTimeAvg/kineticDissipationExact
    if(abs(thermalDissipationExact).GT.1.0d-30) &
        thermalDissipationRatio = epsThermalVolTimeAvg/thermalDissipationExact
#endif

    !-----------------------------------------------------------------------------------------------
    ! 第五部分：初始化 Zhang 表一的热边界层和最小尺度分辨率指标。
    ! 这些量只适用于上下壁恒温的 Rayleigh-BenardCell。
    thermalBLThicknessGlobalEstimate = 0.0d0
    thermalBLThicknessRms = 0.0d0
    thermalBLThicknessRelativeDifferencePercent = 0.0d0
    thermalBLGridPointsGlobalEstimate = 0
    thermalBLGridPointsRms = 0
    etaOverH = 0.0d0
    etaBOverH = 0.0d0
    gridOverEta = 0.0d0
    gridOverEtaB = 0.0d0
    timeStepOverEta = 0.0d0
#ifdef RayleighBenardCell
    ! 在下半区和上半区分别搜索 theta_rms 最大值，避免把上下两个边界层峰值混在一起。
    ! theta_rms(j)=sqrt(<T^2>_x,t-<T>_x,t^2)，max 用于抑制舍入误差造成的微小负数。
    lowerPeakIndex = 1
    upperPeakIndex = ny
    lowerRmsMaximum = -1.0d0
    upperRmsMaximum = -1.0d0
    do j = 1, max(1,ny/2)
        rmsTemperature = dsqrt(max(0.0d0,sumTemperatureSquaredXAvgProfile(j)*invN - &
            (sumTemperatureXAvgProfile(j)*invN)**2))
        if(rmsTemperature.GT.lowerRmsMaximum) then
            lowerRmsMaximum = rmsTemperature
            lowerPeakIndex = j
        endif
    enddo
    do j = min(ny,ny/2+1), ny
        rmsTemperature = dsqrt(max(0.0d0,sumTemperatureSquaredXAvgProfile(j)*invN - &
            (sumTemperatureXAvgProfile(j)*invN)**2))
        if(rmsTemperature.GT.upperRmsMaximum) then
            upperRmsMaximum = rmsTemperature
            upperPeakIndex = j
        endif
    enddo
    ! 保留 Zhang 后文的温度 RMS 峰值定义，便于以后分析边界层结构；
    ! 当前与 Table 1 对照使用的 N_BL 则采用更直接的整体估计 delta_theta/H=1/(2*Nu)。
    thermalBLThicknessRms = 0.5d0*((dble(lowerPeakIndex)-0.5d0) + &
        (dble(ny-upperPeakIndex)+0.5d0))/lengthUnit
    ! ceiling(x) 返回不小于 x 的最小整数，例如 ceiling(10.2)=11；
    ! 这里向上取整可避免把边界层覆盖的非整数网格数向下截断，从而低估边界层分辨率。
    thermalBLGridPointsRms = max(1,ceiling(thermalBLThicknessRms*lengthUnit))
    if(nuVolTimeAvg.GT.0.0d0) then
        ! 整体估计是当前 Table 1 的 N_BL 基准：delta_theta/H=1/(2*Nu)，
        ! 再除以 Delta_g/H=1/lengthUnit，并同样用 ceiling 向上取整为完整网格点数。
        thermalBLThicknessGlobalEstimate = 1.0d0/(2.0d0*nuVolTimeAvg)
        thermalBLGridPointsGlobalEstimate = &
            max(1,ceiling(thermalBLThicknessGlobalEstimate*lengthUnit))
        thermalBLThicknessRelativeDifferencePercent = 100.0d0* &
            abs(thermalBLThicknessRms-thermalBLThicknessGlobalEstimate)/ &
            thermalBLThicknessGlobalEstimate
    endif

    ! 输出完整温度 RMS 剖面。z/H 取流体节点中心位置，首末节点距物理壁面均为半格。
    open(newunit=profileUnit,file=trim(temperatureRmsProfileFile),status='replace',action='write')
    write(profileUnit,'(A)') '# z_over_H theta_rms'
    write(profileUnit,'(A,1X,ES24.16E3)') '# delta_theta_global_estimate_over_H', &
        thermalBLThicknessGlobalEstimate
    write(profileUnit,'(A,1X,ES24.16E3)') '# delta_theta_rms_over_H', thermalBLThicknessRms
    write(profileUnit,'(A,1X,I0)') '# N_BL_global_estimate', thermalBLGridPointsGlobalEstimate
    write(profileUnit,'(A,1X,I0)') '# N_BL_from_temperature_rms_peak', thermalBLGridPointsRms
    write(profileUnit,'(A,1X,ES24.16E3,A)') '# delta_theta_relative_difference_percent', &
        thermalBLThicknessRelativeDifferencePercent, '%'
    do j = 1, ny
        zOverH = (dble(j)-0.5d0)/lengthUnit
        rmsTemperature = dsqrt(max(0.0d0,sumTemperatureSquaredXAvgProfile(j)*invN - &
            (sumTemperatureXAvgProfile(j)*invN)**2))
        write(profileUnit,'(2(ES24.16E3,1X))') zOverH, rmsTemperature
    enddo
    close(profileUnit)

    ! Zhang 表一的三个分辨率指标：Delta_g/eta、Delta_g/eta_B 和 Delta_t/tau_eta。
    ! eta、eta_B、tau_eta 均由时间平均 Nu 的 RB 精确关系得到；只有 Nu>1 时公式有效。
    ! 本代码采用 Delta_g/H=1/lengthUnit、Delta_t/t_ff=1/timeUnit。
    etaOverH = huge(1.0d0)
    etaBOverH = huge(1.0d0)
    gridOverEta = 0.0d0
    gridOverEtaB = 0.0d0
    timeStepOverEta = 0.0d0
    if(nuVolTimeAvg.GT.1.0d0) then
        etaOverH = dsqrt(Prandtl)/(Rayleigh*(nuVolTimeAvg-1.0d0))**0.25d0
        etaBOverH = etaOverH/dsqrt(Prandtl)
        gridOverEta = (1.0d0/lengthUnit)/etaOverH
        gridOverEtaB = (1.0d0/lengthUnit)/etaBOverH
        timeStepOverEta = (1.0d0/timeUnit)/dsqrt(Prandtl/(nuVolTimeAvg-1.0d0))
    endif
#endif

    !-----------------------------------------------------------------------------------------------
    ! 第六部分：把完整文献对照量写入 statisticsFile。
    ! 英文键名保持稳定，便于后续脚本按名称读取；中文物理解释保留在本程序注释中。
    write(22,'(A)') '# Xu et al. Table 4: Nu_volume is the reference for all relative differences'
    write(22,'(A,1X,ES24.16E3)') 'Nu_volume_from_heat_flux', nuVolTimeAvg
    write(22,'(A,1X,ES24.16E3)') 'Re_rms_space_time', reVolTimeRms
    write(22,'(A,1X,ES24.16E3)') 'rho_fluctuation_rms_space_time', densityFluctuationVolTimeRms
    write(22,'(A,1X,ES24.16E3)') 'velocity_divergence_rms_space_time', velocityDivergenceVolTimeRms
#ifdef RayleighBenardCell
    write(22,'(A,1X,ES24.16E3)') 'Nu_wall_from_mean_wall_gradient', nuWallTimeAvg
    write(22,'(A,1X,ES24.16E3,A)') 'Nu_wall_relative_difference_percent', wallNuDifferencePercent, '%'
    write(22,'(A,1X,ES24.16E3)') 'Nu_kinetic_from_mean_dissipation', nuFromKineticDissipation
    write(22,'(A,1X,ES24.16E3,A)') 'Nu_kinetic_relative_difference_percent', kineticNuDifferencePercent, '%'
    write(22,'(A,1X,ES24.16E3)') 'Nu_thermal_from_mean_dissipation', nuFromThermalDissipation
    write(22,'(A,1X,ES24.16E3,A)') 'Nu_thermal_relative_difference_percent', thermalNuDifferencePercent, '%'

    write(22,'(A)') '# Zhang et al.: directly calculated dissipation over exact-relation dissipation'
    write(22,'(A,1X,ES24.16E3)') 'eps_u_calculated_space_time_mean', epsKineticVolTimeAvg
    write(22,'(A,1X,ES24.16E3)') 'eps_u_exact_from_Nu', kineticDissipationExact
    write(22,'(A,1X,ES24.16E3)') 'R_u_calculated_over_exact', kineticDissipationRatio
    write(22,'(A,1X,ES24.16E3)') 'eps_T_calculated_space_time_mean', epsThermalVolTimeAvg
    write(22,'(A,1X,ES24.16E3)') 'eps_T_exact_from_Nu', thermalDissipationExact
    write(22,'(A,1X,ES24.16E3)') 'R_T_calculated_over_exact', thermalDissipationRatio
#else
    write(22,'(A)') '# RB-only dissipation Nu and exact-relation ratios are not output for side-heated flow'
    write(22,'(A,1X,ES24.16E3)') 'eps_u_calculated_space_time_mean', epsKineticVolTimeAvg
    write(22,'(A,1X,ES24.16E3)') 'eps_T_calculated_space_time_mean', epsThermalVolTimeAvg
#endif
    write(22,'(A)') '# Nu/Re half-window convergence is rebuilt once from the complete Nu/Re history'
    write(22,'(A)') '# See NuRe_TimeAverage_2DOpenaccLBMCDE_D2Q5.txt for first/second-half results'
#ifdef RayleighBenardCell
    write(22,'(A)') '# Zhang Table 1: N_BL uses the global estimate delta_theta/H=1/(2*Nu_volume)'
    write(22,'(A,1X,ES24.16E3)') 'delta_theta_global_estimate_over_H', &
        thermalBLThicknessGlobalEstimate
    write(22,'(A,1X,I0)') 'N_BL_global_estimate', thermalBLGridPointsGlobalEstimate
    write(22,'(A)') '# Retained diagnostic: thermal boundary layer from the T_rms peak position'
    write(22,'(A,1X,ES24.16E3)') 'delta_theta_rms_over_H', thermalBLThicknessRms
    write(22,'(A,1X,I0)') 'N_BL_from_temperature_rms_peak', thermalBLGridPointsRms
    write(22,'(A,1X,ES24.16E3,A)') 'delta_theta_relative_difference_percent', &
        thermalBLThicknessRelativeDifferencePercent, '%'
    write(22,'(A,1X,A)') 'temperature_rms_profile_file', trim(temperatureRmsProfileFile)
    write(22,'(A,3(1X,ES24.16E3))') 'grid_over_eta_grid_over_etaB_dt_over_tauEta', &
        gridOverEta, gridOverEtaB, timeStepOverEta
#endif
    write(22,'(A,1X,ES24.16E3)') 'maximum_lattice_CFL', maxStatisticCFL
    close(22)

    !-----------------------------------------------------------------------------------------------
    ! 第七部分：在 settingsFile 末尾追加精简摘要，便于不打开详细统计文件时快速核对结果。
    open(unit=00,file=trim(settingsFile),status='unknown',position='append',action='write')
    write(00,*) 'Cumulative dissipation/statistics samples =', cumulativeStatisticSampleCount
    write(00,*) 'Mean Nu_flux/Re_rms =', nuVolTimeAvg, reVolTimeRms
    write(00,*) 'RMS density fluctuation/divergence =', &
        densityFluctuationVolTimeRms, velocityDivergenceVolTimeRms
    write(00,*) 'Mean eps_u/eps_T =', epsKineticVolTimeAvg, epsThermalVolTimeAvg
#ifdef RayleighBenardCell
    write(00,*) 'Xu Nu_vol/Nu_wall/Nu_kin/Nu_th =', &
        nuVolTimeAvg, nuWallTimeAvg, nuFromKineticDissipation, nuFromThermalDissipation
    write(00,'(A,1X,ES16.8,A,1X,A,1X,ES16.8,A,1X,A,1X,ES16.8,A)') &
        'Xu relative differences wall/kin/th =', wallNuDifferencePercent, '%', '/', &
        kineticNuDifferencePercent, '%', '/', thermalNuDifferencePercent, '%'
    write(00,*) 'Zhang R_u/R_T =', kineticDissipationRatio, thermalDissipationRatio
#endif
#ifdef RayleighBenardCell
    write(00,*) 'Zhang N_BL from H/(2*Nu) =', thermalBLGridPointsGlobalEstimate
    write(00,*) 'Temperature-rms BL grid points (retained) =', thermalBLGridPointsRms
    write(00,*) 'Temperature-rms profile file =', trim(temperatureRmsProfileFile)
#endif
    close(00)
  end subroutine write_dissipation_statistics
!===================================================================================================
! write_dissipation_statistics 结束：只生成最终统计和剖面文件，不改变任何流场变量。
!===================================================================================================


!===================================================================================================
! 子程序: output_unsteady_NuRe_postprocess
! 作用: 从完整统计窗口历史重建瞬时 Nu/Re 曲线、运行平均和前后半段收敛结果。
! 横坐标: 所有非稳态曲线均使用无量纲自由落体时间 t/t_ff，不直接输出格子步 itc。
! 说明: Re_instant(t)=sqrt(<u^2+v^2>_V)*H/nu，只表示当前采样时刻的空间 RMS。
!       文献 Re_rms=sqrt(<u^2+v^2>_V,t)*H/nu，因此运行值必须计算 sqrt(<Re_instant^2>_t)，
!       不能再对 Re_instant 做普通算术平均。运行统计从 unsteadyAverageStartTf 开始，不含初始过渡段。
!===================================================================================================
  subroutine output_unsteady_NuRe_postprocess()
    use commondata
    implicit none

    integer(kind=4) :: k, running_count
    integer(kind=4) :: whole_count, first_count, second_count
    integer(kind=4) :: iosNu, iosRe
    integer(kind=4) :: nuUnit, reUnit, seriesUnit, runningUnit
    real(kind=8) :: currentStatisticSampleTf, nuTf, reTf, NuVal, ReVal
    real(kind=8) :: Nu_RunningSum, ReSquared_RunningSum
    real(kind=8) :: startTf, midTf, endTf
    real(kind=8) :: Nu_WholeSum, ReSquared_WholeSum
    real(kind=8) :: Nu_FirstSum, ReSquared_FirstSum, Nu_SecondSum, ReSquared_SecondSum
    real(kind=8) :: Nu_WholeAvg, Re_WholeAvg, Nu_FirstAvg, Re_FirstAvg, Nu_SecondAvg, Re_SecondAvg
    real(kind=8) :: Nu_FirstVsSecondPercent, Re_FirstVsSecondPercent
    logical :: exNu, exRe, statisticsConverged

    inquire(file='Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat', exist=exNu)
    inquire(file='Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat', exist=exRe)
    if((.not.exNu).or.(.not.exRe)) then
        write(*,'(A)') 'Error: Nu/Re history files are missing before postprocessing.'
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,'(A)') 'Error: Nu/Re history files are missing before postprocessing.'
        close(00)
        error stop 1
    endif
    if(cumulativeStatisticSampleCount.LE.0) then
        write(*,'(A)') 'Error: no Nu/Re history samples were recorded before postprocessing.'
        error stop 1
    endif

    open(newunit=nuUnit, file='Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat', status='old', action='read', form='formatted')
    open(newunit=reUnit, file='Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat', status='old', action='read', form='formatted')

    ! These files are derived views of the full .dat history, so rebuild one continuous ZONE.
    open(newunit=seriesUnit, file='NuRe_InstantaneousVolAvg_2DOpenaccLBMCDE_D2Q5.plt', &
        status='replace', action='write', form='formatted')
    write(seriesUnit,'(A)') 'TITLE = "2D OpenACC instantaneous Nu/Re in statistics window"'
    write(seriesUnit,'(A)') 'VARIABLES = "time_over_tff" "NuVolumeInstantaneous" "ReInstantaneousRms"'
    write(seriesUnit,'(A)') 'ZONE T="NuReInstantaneous", F=POINT'

    open(newunit=runningUnit, file='NuRe_RunningStatistics_2DOpenaccLBMCDE_D2Q5.plt', &
        status='replace', action='write', form='formatted')
    write(runningUnit,'(A)') 'TITLE = "2D OpenACC running statistics in averaging window"'
    write(runningUnit,'(A)') 'VARIABLES = "time_over_tff" "NuVolumeRunningMean" "ReRunningRms"'
    write(runningUnit,'(A)') 'ZONE T="NuReRunningStatistics", F=POINT'

    startTf = unsteadyAverageStartTf
    midTf = unsteadyAverageMidTf
    endTf = unsteadyAverageEndTf

    Nu_RunningSum = 0.0d0
    ReSquared_RunningSum = 0.0d0
    Nu_WholeSum = 0.0d0
    ReSquared_WholeSum = 0.0d0
    Nu_FirstSum = 0.0d0
    ReSquared_FirstSum = 0.0d0
    Nu_SecondSum = 0.0d0
    ReSquared_SecondSum = 0.0d0
    whole_count = 0
    first_count = 0
    second_count = 0
    running_count = 0

    do k = 1, cumulativeStatisticSampleCount
        read(nuUnit,*,iostat=iosNu) nuTf, NuVal
        read(reUnit,*,iostat=iosRe) reTf, ReVal
        ! iostat: 0=成功读到一行；小于0=到达文件末尾；大于0=格式或读入错误。
        ! 循环内只要不是 0，就说明历史短于 latest.meta 记录的样本总数，或者某一行格式坏了。
        if((iosNu.NE.0).OR.(iosRe.NE.0)) then
            write(*,'(A)') 'Error: Nu/Re histories are shorter than cumulativeStatisticSampleCount or invalid.'
            open(unit=00,file=trim(settingsFile),status='unknown',position='append')
            write(00,'(A)') 'Error: Nu/Re histories are shorter than cumulativeStatisticSampleCount or invalid.'
            close(00)
            error stop 1
        endif
        if(abs(nuTf-reTf).GT.1.0d-10*max(1.0d0,abs(nuTf))) then     !确保 Nu 和 Re 是同一个时间采样点的数据，不是错行配对的数据
            write(*,'(A)') 'Error: Nu/Re history time columns do not match.'
            open(unit=00,file=trim(settingsFile),status='unknown',position='append')
            write(00,'(A)') 'Error: Nu/Re history time columns do not match.'
            close(00)
            error stop 1
        endif

        currentStatisticSampleTf = nuTf
        write(seriesUnit,'(ES24.16E3,1X,ES24.16E3,1X,ES24.16E3)') &
            currentStatisticSampleTf, NuVal, ReVal

        if ((currentStatisticSampleTf >= startTf) .and. (currentStatisticSampleTf <= endTf)) then
            ! Nu 对热通量是线性的，使用普通运行平均；Re 必须先累计 Re_instant^2 再开平方。
            running_count = running_count+1
            Nu_RunningSum = Nu_RunningSum+NuVal
            ReSquared_RunningSum = ReSquared_RunningSum+ReVal*ReVal
            write(runningUnit,'(ES24.16E3,1X,ES24.16E3,1X,ES24.16E3)') &
                currentStatisticSampleTf, Nu_RunningSum/dble(running_count), &
                dsqrt(max(0.0d0,ReSquared_RunningSum/dble(running_count)))
            Nu_WholeSum = Nu_WholeSum + NuVal
            ReSquared_WholeSum = ReSquared_WholeSum + ReVal*ReVal
            whole_count = whole_count + 1
        endif

        if ((currentStatisticSampleTf >= startTf) .and. (currentStatisticSampleTf < midTf)) then
            Nu_FirstSum = Nu_FirstSum + NuVal
            ReSquared_FirstSum = ReSquared_FirstSum + ReVal*ReVal
            first_count = first_count + 1
        endif
        if ((currentStatisticSampleTf >= midTf) .and. (currentStatisticSampleTf <= endTf)) then
            Nu_SecondSum = Nu_SecondSum + NuVal
            ReSquared_SecondSum = ReSquared_SecondSum + ReVal*ReVal
            second_count = second_count + 1
        endif
    enddo

    ! 上面的循环已经读完 latest.meta 记录的 cumulativeStatisticSampleCount 行；这里再试读一行，
    ! 不是为了继续计算，而是确认 Nu/Re 两个历史文件后面没有多余数据。
    read(nuUnit,*,iostat=iosNu) nuTf, NuVal
    read(reUnit,*,iostat=iosRe) reTf, ReVal
    ! 任意一个 iostat 等于 0，都表示至少一个文件还成功读到了额外一行。
    ! 如果放过这种情况，后处理会静默丢掉超过重启元数据记录的尾部样本。
    if((iosNu.EQ.0).OR.(iosRe.EQ.0)) then
        write(*,'(A)') 'Error: Nu/Re history files contain more rows than cumulativeStatisticSampleCount.'
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,'(A)') 'Error: Nu/Re history files contain more rows than cumulativeStatisticSampleCount.'
        close(00)
        error stop 1
    endif
    ! 正常结尾必须是两个文件都到达 EOF，也就是两个 iostat 都小于 0。
    ! 其他组合说明一个文件尾部异常，或 Nu/Re 文件长度不一致。
    if(.not.((iosNu.LT.0).AND.(iosRe.LT.0))) then
        write(*,'(A)') 'Error: Nu/Re history files have inconsistent trailing rows.'
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,'(A)') 'Error: Nu/Re history files have inconsistent trailing rows.'
        close(00)
        error stop 1
    endif
    close(nuUnit)
    close(reUnit)
    close(seriesUnit)
    close(runningUnit)

    if ((whole_count <= 0) .or. (first_count <= 0) .or. (second_count <= 0)) then
        write(*,'(A)') 'Error: no complete unsteady average window was found for Nu/Re postprocessing.'
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,'(A)') 'Error: no complete unsteady average window was found for Nu/Re postprocessing.'
        close(00)
        error stop 1
    endif

    Nu_WholeAvg = Nu_WholeSum / dble(whole_count)
    Re_WholeAvg = dsqrt(max(0.0d0,ReSquared_WholeSum/dble(whole_count)))
    Nu_FirstAvg = Nu_FirstSum / dble(first_count)
    Re_FirstAvg = dsqrt(max(0.0d0,ReSquared_FirstSum/dble(first_count)))
    Nu_SecondAvg = Nu_SecondSum / dble(second_count)
    Re_SecondAvg = dsqrt(max(0.0d0,ReSquared_SecondSum/dble(second_count)))
    ! 按统计收敛判据，以后半窗最终结果为分母比较前后半窗；乘 100 后以百分数输出。
    Nu_FirstVsSecondPercent = 100.0d0*abs(Nu_FirstAvg-Nu_SecondAvg) / &
        max(abs(Nu_SecondAvg),tiny(1.0d0))
    Re_FirstVsSecondPercent = 100.0d0*abs(Re_FirstAvg-Re_SecondAvg) / &
        max(abs(Re_SecondAvg),tiny(1.0d0))
    statisticsConverged = (Nu_FirstVsSecondPercent.LT.1.0d0).AND. &
        (Re_FirstVsSecondPercent.LT.1.0d0)
    ! 该判断只输出诊断信息，不控制程序退出或后续时间推进。

    open(unit=33, file='NuRe_TimeAverage_2DOpenaccLBMCDE_D2Q5.txt', &
        status='replace', action='write', form='formatted')
    write(33,'(A)') '# 2D OpenACC Nu mean and literature-defined RMS Re in the averaging window'
    write(33,'(A)') '# final Nu/Re use the second half; relative-error columns are percentages (%), target < 1.0'
    write(33,'(A)') '# start_tf mid_tf end_tf whole_count first_count second_count ' // &
        'Nu_mean_whole Re_rms_whole Nu_mean_first Re_rms_first Nu_final_second Re_final_second ' // &
        'Nu_first_vs_final_percent Re_first_vs_final_percent'
    write(33,'(ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,I0,1X,I0,1X,I0,1X,' // &
        'ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,' // &
        'ES24.16E3,A,1X,ES24.16E3,A)') &
        startTf, midTf, endTf, whole_count, first_count, second_count, &
        Nu_WholeAvg, Re_WholeAvg, Nu_FirstAvg, Re_FirstAvg, Nu_SecondAvg, Re_SecondAvg, &
        Nu_FirstVsSecondPercent, '%', Re_FirstVsSecondPercent, '%'
    if(statisticsConverged) then
        write(33,'(A)') '# 达到统计收敛了。'
    else
        write(33,'(A)') '# 未达到统计收敛，当前计算时间不足，还需要继续往后算。'
    endif
    close(33)

    write(*,'(A)') 'Unsteady Nu/Re statistical postprocessing:'
    write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'Average window t/t_ff:', startTf, 'to', endTf
    write(*,'(A,1X,I0,1X,A,1X,I0,1X,A,1X,I0)') &
        'Samples whole/first/second:', whole_count, '/', first_count, '/', second_count
    write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'Whole Nu mean/Re rms:', Nu_WholeAvg, '/', Re_WholeAvg
    write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'First-half Nu mean/Re rms:', Nu_FirstAvg, '/', Re_FirstAvg
    write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') &
        'Final second-half Nu mean/Re rms:', Nu_SecondAvg, '/', Re_SecondAvg
    write(*,'(A,1X,ES16.8,A,1X,A,1X,ES16.8,A)') &
        'First-half vs final relative error:', Nu_FirstVsSecondPercent, '%', '/', Re_FirstVsSecondPercent, '%'
    if(statisticsConverged) then
        write(*,'(A)') '达到统计收敛了。'
    else
        write(*,'(A)') '未达到统计收敛，当前计算时间不足，还需要继续往后算。'
    endif

  end subroutine output_unsteady_NuRe_postprocess
!===================================================================================================
! output_unsteady_NuRe_postprocess 结束: 输出瞬时曲线和按文献定义累计得到的运行 Nu/Re。
!===================================================================================================
#endif


!===================================================================================================
#ifdef steadyFlow
! 子程序: SideHeatedcalc_Nu_global
! 作用: 计算侧壁差温工况下的全场平均 Nusselt 数。
! 用途: 在 SideHeatedCell 工况结束后的后处理中调用。
!===================================================================================================
  subroutine SideHeatedcalc_Nu_global()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: dx, dT, qx, sum_qx
    real(kind=8) :: temperatureDifference, coef

    ! 网格间距
    dx = 1.0d0 / lengthUnit
    temperatureDifference = Thot - Tcold
    coef   = velocityScaleCompare

    sum_qx = 0.0d0

    !$acc parallel loop collapse(2) default(none) present(u,T) reduction(+:sum_qx) private(dT,qx)
    do j = 1, ny
      do i = 1, nx

        if (i == 1) then
          ! i=1: 节点位于 x=dx/2，利用 (wall, i=1, i=2) 二次插值给出 dT/dx在x=dx/2 的二阶近似（边界特别处理）
          dT = (-3.0d0*T(1,j) - T(2,j) + 4.0d0*Thot ) / (3.0d0*dx)
        elseif (i == nx) then
          ! i=nx: 节点位于 x=L-dx/2，利用 (i=nx-1, i=nx, wall) 二次插值给出 dT/dx在x=L-dx/2 的二阶近似（边界特别处理）
          dT = ( -4.0d0*Tcold + 3.0d0*T(nx,j) + T(nx-1,j) ) / (3.0d0*dx)
        else
          ! 1<i<nx: 中心差分
          dT = ( T(i-1,j) - T(i+1,j) ) / (2.0d0*dx)
        endif

        qx = coef*u(i,j)*(T(i,j)-Tref) + dT
        sum_qx = sum_qx + qx

      enddo
    enddo
    Nu_global = (sum_qx / dble(nx*ny)) / temperatureDifference

    ! 屏幕输出
    write(*,'(a,1x,ES24.16E3)') "Nu_global =", real(Nu_global,kind=8)


    ! 同步写入日志
    open(unit=00,file=trim(settingsFile),status="unknown",position="append")
    write(00,'(a,1x,ES24.16E3)') "Nu_global =", real(Nu_global,kind=8)
    close(00)

    return
  end subroutine SideHeatedcalc_Nu_global
!===================================================================================================
! SideHeatedcalc_Nu_global 结束: 计算侧壁差温工况下的全场平均 Nusselt 数。
!===================================================================================================



!===================================================================================================
! 子程序: SideHeatedcalc_Nu_wall_avg
! 作用: 计算侧壁差温工况下热壁、冷壁和中线的 Nusselt 数及其极值。
! 用途: 在 SideHeatedCell 工况结束后的后处理中调用。
!===================================================================================================
  subroutine SideHeatedcalc_Nu_wall_avg()
    use commondata
    implicit none
    integer(kind=4) :: j, iMid
    integer(kind=4) :: jmax, jmin
    real(kind=8) :: dx, dy, temperatureDifference, coef
    real(kind=8) :: qx_wall, sum_hot, sum_cold, sum_mid
    real(kind=8) :: denom
    real(kind=8) :: Nu_left(1:ny)
    real(kind=8) :: f_m2,f_m1,f_0,f_p1,f_p2
    !------------------------------------------------------------
    ! 5-point least-squares parabola fit (general, can be one-sided)
    !------------------------------------------------------------
    integer(kind=4) :: k,iL,iR
    integer(kind=4) :: jj(5)
    real(kind=8) :: xk(5), fk(5), yk(5)
    real(kind=8) :: S0, S1, S2, S3, S4
    real(kind=8) :: F0, F1, F2
    real(kind=8) :: D, Da, Db, Dc
    real(kind=8) :: a, b, c
    real(kind=8) :: delta, fstar, ystar
    real(kind=8) :: xmin, xmax
    real(kind=8), parameter :: epsD = 1.0d-20, epsA = 1.0d-14
    real(kind=8) :: Nu_left_ext(0:ny+1)
    real(kind=8) :: T_wb, T_wt     ! wall temperature at y=0 and y=1 on i=1 vertical line
    real(kind=8) :: yfit(4), Tfit(4)


    dx = 1.0d0 / lengthUnit
    dy = 1.0d0 / lengthUnit
    temperatureDifference = Thot - Tcold
    coef   = velocityScaleCompare



    !-----------------------------
    ! (1) 左侧热壁平均 Nu_hot，同时记录 Numax/Numin 及其 y 位置
    sum_hot = 0.0d0
    do j = 1, ny
      ! 壁面导热通量：qx(x=0,j)
      qx_wall = (8.0d0*Thot - 9.0d0*T(1,j) + T(2,j)) / (3.0d0*dx)
      Nu_left(j)= qx_wall / temperatureDifference
      sum_hot   = sum_hot + Nu_left(j)
    enddo
    Nu_hot = sum_hot / dble(ny)

    ! 计算左壁面上下两个角点的热通量
    ! 先把中间 j=1..ny 的值复制过来
    Nu_left_ext(1:ny) = Nu_left(1:ny)

    ! ---------- 左下角：在 i=1 这条竖线上，用 j=1..4 拟合得到 y=0 的温度 ----------
    yfit(1) = yp(1);  Tfit(1) = T(1,1)
    yfit(2) = yp(2);  Tfit(2) = T(1,2)
    yfit(3) = yp(3);  Tfit(3) = T(1,3)
    yfit(4) = yp(4);  Tfit(4) = T(1,4)

    call fit_adiabatic_wall_T4(0.0d0, yfit, Tfit, T_wb)   ! 得到 T(y=0) = T_wb,拟合绝热壁面温度，用4个点

    Nu_left_ext(0) = (2.0d0*(Thot-T_wb)/dx)/temperatureDifference   ! 角点局部 Nu

    ! ---------- 左上角：在 i=1 这条竖线上，用 j=ny-3..ny 拟合得到 y=1 的温度 ----------
    yfit(1) = yp(ny-3);  Tfit(1) = T(1,ny-3)
    yfit(2) = yp(ny-2);  Tfit(2) = T(1,ny-2)
    yfit(3) = yp(ny-1);  Tfit(3) = T(1,ny-1)
    yfit(4) = yp(ny  );  Tfit(4) = T(1,ny  )

    call fit_adiabatic_wall_T4(yp(ny+1), yfit, Tfit, T_wt)   ! 得到顶壁温度 T(y=yp(ny+1))

    Nu_left_ext(ny+1) = (2.0d0*(Thot-T_wt)/dx)/temperatureDifference  ! 角点局部 Nu





    ! 网格上先找最大/最小
    jmax = 0
    jmin = 0
    Nu_hot_max = Nu_left_ext(0)
    Nu_hot_min = Nu_left_ext(0)

    do j = 1, ny+1
      if (Nu_left_ext(j) > Nu_hot_max) then
        Nu_hot_max = Nu_left_ext(j)
        jmax = j
      endif
      if (Nu_left_ext(j) < Nu_hot_min) then
        Nu_hot_min = Nu_left_ext(j)
        jmin = j
      endif
    enddo


    !==========================
    ! 五点插值 least-squares parabola fit, else one-sided 5pts
    !==========================
        ! ----------- 选 Numax 的 5 个点 -----------
    if (jmax <= 2) then
      jj = (/ 0, 1, 2, 3, 4 /)
    elseif (jmax >= ny-1) then
      jj = (/ ny-3, ny-2, ny-1, ny, ny+1 /)
    else
      jj = (/ jmax-2, jmax-1, jmax, jmax+1, jmax+2 /)
    endif

    do k = 1, 5
      yk(k) = yp(jj(k))            ! yp(0)=0, yp(ny+1)=ny/lengthUnit 已经存在
      fk(k) = Nu_left_ext(jj(k))
    enddo

    call fit_parabola_ls5(yk, fk, +1, fstar, ystar)  !最小二乘法，用抛物线，拟合五个点

    Nu_hot_max = fstar
    Nu_hot_max_position = ystar



    ! ----------- 选 Numin 的 5 个点：单侧（before the minimum）-----------
    if (jmin >= 4) then
      jj = (/ jmin-4, jmin-3, jmin-2, jmin-1, jmin /)
    else
      jj = (/ 0, 1, 2, 3, 4 /)
    endif

    do k = 1, 5
      yk(k) = yp(jj(k))
      fk(k) = Nu_left_ext(jj(k))
    enddo

    call fit_parabola_ls5(yk, fk, -1, fstar, ystar)

    Nu_hot_min = fstar
    Nu_hot_min_position = ystar


    !-----------------------------
    ! (2) 右侧冷壁平均 Nu_cold
    sum_cold = 0.0d0
    do j = 1, ny
      qx_wall = (-8.0d0*Tcold + 9.0d0*T(nx,j) - T(nx-1,j)) / (3.0d0*dx)
      sum_cold = sum_cold+qx_wall/temperatureDifference
    enddo
    Nu_cold = (sum_cold / dble(ny))

    !-----------------------------
    ! (3) 竖直中线 x=1/2 的平均 Nu_middle
    sum_mid = 0.0d0

    if (mod(nx,2) == 1) then
      iMid = (nx + 1)/2

      do j = 1, ny
        sum_mid = sum_mid+(coef*u(iMid,j)*(T(iMid,j)-Tref) + &
            (T(iMid-1,j)-T(iMid+1,j))/(2.0d0*dx))/temperatureDifference
      enddo

    else
      iL = nx/2
      iR = iL + 1

      do j = 1, ny
        sum_mid = sum_mid + (coef*( 0.5d0*( u(iL,j)*(T(iL,j)-Tref) + u(iR,j)*(T(iR,j)-Tref) )) &
        + (T(iL,j)-T(iR,j))/dx)/temperatureDifference
      enddo
    endif

    Nu_middle = (sum_mid / dble(ny))

    !-----------------------------
    ! 输出：屏幕 + 日志
    write(*,'(a,1x,ES24.16E3)') "Nu_hot    =", real(Nu_hot,kind=8)
    write(*,'(a,1x,ES24.16E3)') "Nu_cold   =", real(Nu_cold,kind=8)
    write(*,'(a,1x,ES24.16E3)') "Nu_middle =", real(Nu_middle,kind=8)
    write(*,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
        "Nu_hot_max =", real(Nu_hot_max,kind=8), "y_max =", real(Nu_hot_max_position,kind=8)
    write(*,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
        "Nu_hot_min =", real(Nu_hot_min,kind=8), "y_min =", real(Nu_hot_min_position,kind=8)



    open(unit=00,file=trim(settingsFile),status="unknown",position="append")
    write(00,'(a,1x,ES24.16E3)') "Nu_hot    =", real(Nu_hot,kind=8)
    write(00,'(a,1x,ES24.16E3)') "Nu_cold   =", real(Nu_cold,kind=8)
    write(00,'(a,1x,ES24.16E3)') "Nu_middle =", real(Nu_middle,kind=8)
    write(00,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
        "Nu_hot_max =", real(Nu_hot_max,kind=8), "y_max =", real(Nu_hot_max_position,kind=8)
    write(00,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
        "Nu_hot_min =", real(Nu_hot_min,kind=8), "y_min =", real(Nu_hot_min_position,kind=8)
    close(00)

    return
  end subroutine SideHeatedcalc_Nu_wall_avg
!===================================================================================================
! SideHeatedcalc_Nu_wall_avg 结束: 计算侧壁差温工况下热壁、冷壁和中线的 Nusselt 数及其极值。
!===================================================================================================
#endif



!===================================================================================================
! 子程序: fit_adiabatic_wall_T4
! 作用: 用四点拟合估计绝热壁面的壁温。
! 用途: 在 SideHeated 和 RB 壁面 Nusselt 后处理中调用。
!===================================================================================================
  subroutine fit_adiabatic_wall_T4(y0, y, tt, T_wall)
    implicit none
    real(kind=8), intent(in)  :: y0
    real(kind=8), intent(in)  :: y(4), tt(4)
    real(kind=8), intent(out) :: T_wall
    real(kind=8) :: s(4)
    real(kind=8) :: S0, S1, S2, B0, B1, D
    integer(kind=4) :: k

    do k = 1, 4
      s(k) = (y(k) - y0)*(y(k) - y0)
    enddo

    S0 = 4.0d0
    S1 = 0.0d0
    S2 = 0.0d0
    B0 = 0.0d0
    B1 = 0.0d0
    do k = 1, 4
      S1 = S1 + s(k)
      S2 = S2 + s(k)*s(k)
      B0 = B0 + tt(k)
      B1 = B1 + tt(k)*s(k)
    enddo

    D = S0*S2 - S1*S1

    ! T_wall = c0
    T_wall = (B0*S2 - B1*S1) / D

    return
  end subroutine fit_adiabatic_wall_T4
!===================================================================================================
! fit_adiabatic_wall_T4 结束: 用四点拟合估计绝热壁面的壁温。
!===================================================================================================


!===================================================================================================
! 子程序: fit_parabola_ls5
! 作用: 用五点最小二乘抛物线拟合局部极值和对应位置。
! 用途: 在 Nu 极值和中心线速度极值的后处理中重复调用。
!===================================================================================================
  subroutine fit_parabola_ls5(y, f, mode, fstar, ystar)
    implicit none
    real(kind=8), intent(in)  :: y(5), f(5)
    real(kind=8), intent(out) :: fstar, ystar
    real(kind=8) :: S0, S1, S2, S3, S4
    real(kind=8) :: F0, F1, F2
    real(kind=8) :: D, DA, DB, DC
    real(kind=8) :: A, B, C
    real(kind=8) :: ymin, ymax
    integer(kind=4) :: k, kbest
    real(kind=8), parameter :: epsD=1.0d-20, epsA=1.0d-14
    integer(kind=4), intent(in) :: mode   ! +1 => max, -1 => min


  ! ----- fallback: pick max/min among the 5 samples -----
  kbest = 1
  do k = 2, 5
    if (mode == 1) then
      if (f(k) > f(kbest)) kbest = k
    else
      if (f(k) < f(kbest)) kbest = k
    endif
  enddo

    S0=0.0d0; S1=0.0d0; S2=0.0d0; S3=0.0d0; S4=0.0d0
    F0=0.0d0; F1=0.0d0; F2=0.0d0
    do k = 1, 5
      S0 = S0 + 1.0d0
      S1 = S1 + y(k)
      S2 = S2 + y(k)*y(k)
      S3 = S3 + y(k)*y(k)*y(k)
      S4 = S4 + y(k)*y(k)*y(k)*y(k)

      F0 = F0 + f(k)
      F1 = F1 + y(k)*f(k)
      F2 = F2 + y(k)*y(k)*f(k)
    enddo

    ! Solve:
    ! [S4 S3 S2][A] = [F2]
    ! [S3 S2 S1][B]   [F1]
    ! [S2 S1 S0][C]   [F0]
    D  =  S4*(S2*S0 - S1*S1) - S3*(S3*S0 - S1*S2) + S2*(S3*S1 - S2*S2)
    DA =  F2*(S2*S0 - S1*S1) - S3*(F1*S0 - S1*F0) + S2*(F1*S1 - S2*F0)
    DB =  S4*(F1*S0 - S1*F0) - F2*(S3*S0 - S1*S2) + S2*(S3*F0 - F1*S2)
    DC =  S4*(S2*F0 - F1*S1) - S3*(S3*F0 - F1*S2) + F2*(S3*S1 - S2*S2)

    if (dabs(D) > epsD) then
      A = DA / D
      B = DB / D
      C = DC / D

      if (dabs(A) > epsA) then
        ystar = -B / (2.0d0*A)
        ymin = minval(y);  ymax = maxval(y)

        if (ystar >= ymin .and. ystar <= ymax) then
          fstar = C - B*B/(4.0d0*A)
        else
          ! 顶点落在拟合区间外：退回到 5 点中最大的/最小的那个并不严格，
          ! 这里为了保持简单，通常会先选极值点
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

    return
  end subroutine fit_parabola_ls5
!===================================================================================================
! fit_parabola_ls5 结束: 用五点最小二乘抛物线拟合局部极值和对应位置。
!===================================================================================================



!===================================================================================================
! 子程序: SideHeatedcalc_umid_max
! 作用: 计算侧壁差温工况下中心线水平速度的最大值及位置。
! 用途: 在 SideHeatedCell 工况结束后的后处理中调用。
!===================================================================================================
  subroutine SideHeatedcalc_umid_max()
    use commondata
    implicit none
    integer(kind=4) :: j, k
    integer(kind=4) :: iMid, iL, iR
    integer(kind=4) :: j0
    integer(kind=4) :: jj(5)
    real(kind=8) :: uline(1:ny)
    real(kind=8) :: s(5), fu(5)
    real(kind=8) :: umax_grid, umax_fit, y_fit
    real(kind=8) :: xmid
    character(len=24) :: ctime, string
    integer(kind=4) :: time
    real(kind=8) :: coef

    coef = velocityScaleCompare

    ! ---- (1) 构造中线剖面 u(x=1/2, y_j) ----
    if (mod(nx,2) == 1) then
      iMid = (nx + 1)/2
      xmid = xp(iMid)
      do j = 1, ny
        uline(j) = u(iMid,j)
      enddo
    else
      iL = nx/2
      iR = iL + 1
      xmid = 0.5d0*(xp(iL) + xp(iR))
      do j = 1, ny
        uline(j) = 0.5d0*(u(iL,j) + u(iR,j))
      enddo
    endif

    ! ---- (2) 先找网格最大值点 j0 ----
    j0 = 1
    umax_grid = uline(1)
    do j = 2, ny
      if (uline(j) > umax_grid) then
        umax_grid = uline(j)
        j0 = j
      endif
    enddo

    ! ---- (3) 取 5 点（尽量对称；靠近端点时偏侧） ----
    if (j0 <= 2) then
      jj = (/ 1, 2, 3, 4, 5 /)
    elseif (j0 >= ny-1) then
      jj = (/ ny-4, ny-3, ny-2, ny-1, ny /)
    else
      jj = (/ j0-2, j0-1, j0, j0+1, j0+2 /)
    endif

    do k = 1, 5
      s(k) = yp(jj(k))       ! 自变量：y
     fu(k) = uline(jj(k))    ! 拟合量：u
    enddo

    call fit_parabola_ls5(s, fu, +1, umax_fit, y_fit)

    ! ---- 输出 ----
    write(*,'(A,1X,ES24.16E3,1X,A,1X,ES24.16E3,1X,A,1X,ES24.16E3)') &
         'u_mid_max =', umax_fit*coef, 'at y =', y_fit, 'on x_mid =', xmid


    open(unit=00,file=trim(settingsFile),status='unknown',position='append')
    string = ctime( time() )
    write(00,*) '--- calc_umid_max --- ', string
    write(00,'(A,1X,ES24.16E3)') 'x_mid =', real(xmid,kind=8)
    write(00,'(A,1X,ES24.16E3,1X,A,1X,ES24.16E3,1X,A,I0,A)') &
        'u_mid_max =', real(umax_fit*coef,kind=8), 'y_pos =', real(y_fit,kind=8), ' (grid j0=', j0, ')'
    close(00)

    return
  end subroutine SideHeatedcalc_umid_max
!===================================================================================================
! SideHeatedcalc_umid_max 结束: 计算侧壁差温工况下中心线水平速度的最大值及位置。
!===================================================================================================




!===================================================================================================
! 子程序: SideHeatedcalc_vmid_max
! 作用: 计算侧壁差温工况下中心线竖直速度的最大值及位置。
! 用途: 在 SideHeatedCell 工况结束后的后处理中调用。
!===================================================================================================
  subroutine SideHeatedcalc_vmid_max()
    use commondata
    implicit none
    integer(kind=4) :: i, k
    integer(kind=4) :: jMid, jB, jT
    integer(kind=4) :: i0
    integer(kind=4) :: ii(5)
    real(kind=8) :: vline(1:nx)
    real(kind=8) :: s(5), fv(5)
    real(kind=8) :: vmax_grid, vmax_fit, x_fit
    real(kind=8) :: ymid
    character(len=24) :: ctime, string
    integer(kind=4) :: time
    real(kind=8) :: coef

    coef = velocityScaleCompare

    ! ---- (1) 构造中线剖面 v(x_i, y=1/2) ----
    if (mod(ny,2) == 1) then
      jMid = (ny + 1)/2
      ymid = yp(jMid)
      do i = 1, nx
        vline(i) = v(i,jMid)
      enddo
    else
      jB = ny/2
      jT = jB + 1
      ymid = 0.5d0*(yp(jB) + yp(jT))
      do i = 1, nx
        vline(i) = 0.5d0*(v(i,jB) + v(i,jT))
      enddo
    endif

    ! ---- (2) 先找网格最大值点 i0 ----
    i0 = 1
    vmax_grid = vline(1)
    do i = 2, nx
      if (vline(i) > vmax_grid) then
        vmax_grid = vline(i)
        i0 = i
      endif
    enddo

    ! ---- (3) 取 5 点（尽量对称；靠近端点时偏侧） ----
    if (i0 <= 2) then
      ii = (/ 1, 2, 3, 4, 5 /)
    elseif (i0 >= nx-1) then
      ii = (/ nx-4, nx-3, nx-2, nx-1, nx /)
    else
      ii = (/ i0-2, i0-1, i0, i0+1, i0+2 /)
    endif

    do k = 1, 5
      s(k) = xp(ii(k))       ! 自变量：x
     fv(k) = vline(ii(k))    ! 拟合量：v
    enddo

    call fit_parabola_ls5(s, fv, +1, vmax_fit, x_fit)

    ! ---- 输出 ----
    write(*,'(A,1X,ES24.16E3,1X,A,1X,ES24.16E3,1X,A,1X,ES24.16E3)') &
         'v_mid_max =', vmax_fit*coef, 'at x =', x_fit, 'on y_mid =', ymid


    open(unit=00,file=trim(settingsFile),status='unknown',position='append')
    string = ctime( time() )
    write(00,*) '--- calc_vmid_max --- ', string
    write(00,'(A,1X,ES24.16E3)') 'y_mid =', real(ymid,kind=8)
    write(00,'(A,1X,ES24.16E3,1X,A,1X,ES24.16E3,1X,A,I0,A)') &
        'v_mid_max =', real(vmax_fit*coef,kind=8), 'x_pos =', real(x_fit,kind=8), ' (grid i0=', i0, ')'
    close(00)

    return
  end subroutine SideHeatedcalc_vmid_max
!===================================================================================================
! SideHeatedcalc_vmid_max 结束: 计算侧壁差温工况下中心线竖直速度的最大值及位置。
!===================================================================================================




!===================================================================================================
#ifdef steadyFlow
! 子程序: RBcalc_Nu_global
! 作用: 计算 Rayleigh-Benard 工况下的全场平均 Nusselt 数。
! 用途: 在 RayleighBenardCell 工况结束后的后处理中调用。
!===================================================================================================
subroutine RBcalc_Nu_global()
  use commondata
  implicit none
  integer(kind=4) :: i, j
  real(kind=8) :: dy, dTdy, qy, sum_qy
  real(kind=8) :: temperatureDifference, coef

  dy     = 1.0d0 / lengthUnit
  temperatureDifference = Thot - Tcold
  coef   = velocityScaleCompare

  sum_qy = 0.0d0

  !$acc parallel loop collapse(2) default(none) present(v,T) reduction(+:sum_qy) private(dTdy,qy)
  do j = 1, ny
    do i = 1, nx

      if (j == 1) then
        ! y=dy/2：用 (wall, j=1, j=2) 二次插值给 dT/dy 的二阶近似
        dTdy = ( 3.0d0*T(i,1) + T(i,2) - 4.0d0*Thot ) / (3.0d0*dy)

      elseif (j == ny) then
        ! y=1-dy/2：用 (j=ny-1, j=ny, wall)
        dTdy = ( 4.0d0*Tcold - 3.0d0*T(i,ny) - T(i,ny-1) ) / (3.0d0*dy)

      else
        dTdy = ( T(i,j+1) - T(i,j-1) ) / (2.0d0*dy)
      endif

      qy = coef * v(i,j) * (T(i,j) - Tref) - dTdy
      sum_qy = sum_qy + qy

    enddo
  enddo
  Nu_global = (sum_qy / dble(nx*ny)) / temperatureDifference

  write(*,'(a,1x,ES24.16E3)') "Nu_global =", real(Nu_global,kind=8)
  open(unit=00,file=trim(settingsFile),status="unknown",position="append")
  write(00,'(a,1x,ES24.16E3)') "Nu_global =", real(Nu_global,kind=8)
  close(00)

  return
end subroutine RBcalc_Nu_global
!===================================================================================================
! RBcalc_Nu_global 结束: 计算 Rayleigh-Benard 工况下的全场平均 Nusselt 数。
!===================================================================================================


!===================================================================================================
! 子程序: RBcalc_Nu_wall_avg
! 作用: 计算 Rayleigh-Benard 工况下热壁、冷壁和中线的 Nusselt 数及其极值。
! 用途: 在 RayleighBenardCell 工况结束后的后处理中调用。
!===================================================================================================
subroutine RBcalc_Nu_wall_avg()
  use commondata
  implicit none
  integer(kind=4) :: i, k
  integer(kind=4) :: imax, imin
  integer(kind=4) :: jMid, jB, jT
  integer(kind=4) :: ii(5)
  real(kind=8) :: dx, dy, temperatureDifference
  real(kind=8) :: qy_wall, sum_hot, sum_cold, sum_mid, coef
  real(kind=8), dimension(1:nx) :: Nu_bot
  real(kind=8), dimension(0:nx+1) :: Nu_bot_ext
  real(kind=8) :: xfit(4), Tfit(4), T_wl, T_wr
  real(kind=8) :: xk(5), fk(5)
  real(kind=8) :: fstar, xstar


  dx     = 1.0d0 / lengthUnit
  dy     = 1.0d0 / lengthUnit
  temperatureDifference = Thot - Tcold
  coef   = velocityScaleCompare

  !-----------------------------
  ! (1) 底部热壁平均 Nu_hot（不含角点）
  sum_hot = 0.0d0
  do i = 1, nx
    qy_wall   = (8.0d0*Thot - 9.0d0*T(i,1) + T(i,2)) / (3.0d0*dy)
    Nu_bot(i)= qy_wall / temperatureDifference
    sum_hot  = sum_hot + Nu_bot(i)
  enddo
  Nu_hot = sum_hot / dble(nx)

  !-----------------------------
  ! (1.1) 角点扩展：用侧壁绝热（Neumann）4点拟合得到 x=0 与 x=1 处 y=dy/2 的温度
  ! 左下角附近：i=1..4, j=1
  xfit(1)=xp(1);  Tfit(1)=T(1,1)
  xfit(2)=xp(2);  Tfit(2)=T(2,1)
  xfit(3)=xp(3);  Tfit(3)=T(3,1)
  xfit(4)=xp(4);  Tfit(4)=T(4,1)
  call fit_adiabatic_wall_T4(0.0d0, xfit, Tfit, T_wl)   ! 估计 T(x=0, y=dy/2)

  ! 右下角附近：i=nx-3..nx, j=1
  xfit(1)=xp(nx-3);  Tfit(1)=T(nx-3,1)
  xfit(2)=xp(nx-2);  Tfit(2)=T(nx-2,1)
  xfit(3)=xp(nx-1);  Tfit(3)=T(nx-1,1)
  xfit(4)=xp(nx  );  Tfit(4)=T(nx  ,1)
  call fit_adiabatic_wall_T4(xp(nx+1), xfit, Tfit, T_wr)   ! 估计 T(x=xp(nx+1), y=dy/2)

  ! 组装扩展数组：角点只用于找 max/min 与拟合
  Nu_bot_ext(1:nx) = Nu_bot(1:nx)
  Nu_bot_ext(0)    = (2.0d0*(Thot-T_wl)/dy)/temperatureDifference
  Nu_bot_ext(nx+1) = (2.0d0*(Thot-T_wr)/dy)/temperatureDifference

  !-----------------------------
  ! (1.2) 网格上找 Numax/Numin（含角点 0 与 nx+1）
  imax = 0
  imin = 0
  Nu_hot_max = Nu_bot_ext(0)
  Nu_hot_min = Nu_bot_ext(0)

  do i = 1, nx+1
    if (Nu_bot_ext(i) > Nu_hot_max) then
      Nu_hot_max = Nu_bot_ext(i)
      imax = i
    endif
    if (Nu_bot_ext(i) < Nu_hot_min) then
      Nu_hot_min = Nu_bot_ext(i)
      imin = i
    endif
  enddo

  !-----------------------------
  ! (1.3) 对 Numax：取 5 点（尽量对称；靠近端点时单侧）并做 LS 抛物线拟合
  if (imax <= 2) then
    ii = (/ 0, 1, 2, 3, 4 /)
  elseif (imax >= nx-1) then
    ii = (/ nx-3, nx-2, nx-1, nx, nx+1 /)
  else
    ii = (/ imax-2, imax-1, imax, imax+1, imax+2 /)
  endif

  do k = 1, 5
    xk(k) = xp(ii(k))            ! xp(0)=0, xp(nx+1)=nx/lengthUnit
    fk(k) = Nu_bot_ext(ii(k))
  enddo

  call fit_parabola_ls5(xk, fk, +1, fstar, xstar)
  Nu_hot_max = fstar
  Nu_hot_max_position = xstar     ! RB: 这里是 x 位置

  !-----------------------------
  ! (1.4) 对 Numin：同样取 5 点并拟合
  if (imin <= 2) then
    ii = (/ 0, 1, 2, 3, 4 /)
  elseif (imin >= nx-1) then
    ii = (/ nx-3, nx-2, nx-1, nx, nx+1 /)
  else
    ii = (/ imin-2, imin-1, imin, imin+1, imin+2 /)
  endif

  do k = 1, 5
    xk(k) = xp(ii(k))
    fk(k) = Nu_bot_ext(ii(k))
  enddo

  call fit_parabola_ls5(xk, fk, -1, fstar, xstar)
  Nu_hot_min = fstar
  Nu_hot_min_position = xstar     ! RB: 这里是 x 位置

  !-----------------------------
  ! (2) 顶部冷壁平均 Nu_cold（不含角点）
  sum_cold = 0.0d0
  do i = 1, nx
    qy_wall = (-8.0d0*Tcold + 9.0d0*T(i,ny) - T(i,ny-1)) / (3.0d0*dy)
    sum_cold = sum_cold+qy_wall/temperatureDifference
  enddo
  Nu_cold = sum_cold / dble(nx)

  !-----------------------------
  ! 中线的 Nusselt 数
  sum_mid = 0.0d0

  if (mod(ny,2) == 1) then
    jMid = (ny + 1)/2

    do i = 1, nx
      sum_mid = sum_mid+(coef*v(i,jMid)*(T(i,jMid)-Tref) - &
          (T(i,jMid+1)-T(i,jMid-1))/(2.0d0*dy))/temperatureDifference
    enddo

  else
    jB = ny/2
    jT = jB + 1

    do i = 1, nx
      sum_mid = sum_mid + (coef*( 0.5d0*( v(i,jB)*(T(i,jB)-Tref) + v(i,jT)*(T(i,jT)-Tref) )) &
      + (T(i,jB)-T(i,jT))/dy)/temperatureDifference
    enddo
  endif

  Nu_middle = sum_mid / dble(nx)

  !-----------------------------
  ! 输出：屏幕 + 日志
  write(*,'(a,1x,ES24.16E3)') "Nu_hot(bottom) =", real(Nu_hot,kind=8)
  write(*,'(a,1x,ES24.16E3)') "Nu_cold(top)   =", real(Nu_cold,kind=8)
  write(*,'(a,1x,ES24.16E3)') "Nu_middle      =", real(Nu_middle,kind=8)
  write(*,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
      "Nu_hot_max =", real(Nu_hot_max,kind=8), "x_max =", real(Nu_hot_max_position,kind=8)
  write(*,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
      "Nu_hot_min =", real(Nu_hot_min,kind=8), "x_min =", real(Nu_hot_min_position,kind=8)

  open(unit=00,file=trim(settingsFile),status="unknown",position="append")
  write(00,'(a,1x,ES24.16E3)') "Nu_hot(bottom) =", real(Nu_hot,kind=8)
  write(00,'(a,1x,ES24.16E3)') "Nu_cold(top)   =", real(Nu_cold,kind=8)
  write(00,'(a,1x,ES24.16E3)') "Nu_middle      =", real(Nu_middle,kind=8)
  write(00,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
      "Nu_hot_max =", real(Nu_hot_max,kind=8), "x_max =", real(Nu_hot_max_position,kind=8)
  write(00,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
      "Nu_hot_min =", real(Nu_hot_min,kind=8), "x_min =", real(Nu_hot_min_position,kind=8)
  close(00)

  return
end subroutine RBcalc_Nu_wall_avg
!===================================================================================================
! RBcalc_Nu_wall_avg 结束: 计算 Rayleigh-Benard 工况下热壁、冷壁和中线的 Nusselt 数及其极值。
!===================================================================================================
#endif


!===================================================================================================
! 子程序: RBcalc_umid_max
! 作用: 计算 Rayleigh-Benard 工况下 x=1/2 处最大水平速度 umax 及其 y 位置。
! 用途: 在 RayleighBenardCell 工况结束后的后处理中调用。
!===================================================================================================
subroutine RBcalc_umid_max()
  use commondata
  implicit none
  integer(kind=4) :: j, k
  integer(kind=4) :: iMid, iL, iR
  integer(kind=4) :: j0
  integer(kind=4) :: jj(5)
  real(kind=8) :: uline(1:ny)
  real(kind=8) :: yk(5), fk(5)
  real(kind=8) :: umax_fit, y_fit
  real(kind=8) :: coef, xmid, targetX, w
  real(kind=8) :: umax_grid, yLen

  coef = velocityScaleCompare

  ! ---- (1) 取论文定义的 x=1/2 竖线剖面 u(x=0.5, y_j) ----
  targetX = 0.5d0
  xmid = targetX

  iL = 1
  do while (iL < nx .and. xp(iL+1) < targetX)
    iL = iL + 1
  enddo
  iR = min(iL + 1, nx)

  if (dabs(xp(iL) - targetX) <= 1.0d-14) then
    do j = 1, ny
      uline(j) = u(iL,j)
    enddo
  elseif (iR == iL) then
    do j = 1, ny
      uline(j) = u(iL,j)
    enddo
  else
    w = (targetX - xp(iL)) / (xp(iR) - xp(iL))
    do j = 1, ny
      uline(j) = (1.0d0 - w) * u(iL,j) + w * u(iR,j)
    enddo
  endif

  ! ---- (2) 找峰值所在网格点 j0 ----
  j0 = 1
  umax_grid = uline(1)
  do j = 2, ny
    if (uline(j) > umax_grid) then
      umax_grid = uline(j)
      j0 = j
    endif
  enddo

  ! ---- (3) 取 5 点（尽量对称；靠近端点时偏侧）----
  if (j0 <= 2) then
    jj = (/ 1, 2, 3, 4, 5 /)
  elseif (j0 >= ny-1) then
    jj = (/ ny-4, ny-3, ny-2, ny-1, ny /)
  else
    jj = (/ j0-2, j0-1, j0, j0+1, j0+2 /)
  endif

  do k = 1, 5
    yk(k) = yp(jj(k))
    fk(k) = uline(jj(k))
  enddo

  call fit_parabola_ls5(yk, fk, +1, umax_fit, y_fit)

  ! 对称支路下，上下半区的 y 位置等价；统一折回下半区便于与 benchmark 对照。
  yLen = yp(ny+1)
  if (y_fit > 0.5d0*yLen) y_fit = yLen - y_fit


  write(*,'(A,1X,ES24.16E3,2X,A,1X,ES24.16E3,2X,A,1X,ES24.16E3)') &
       'u_mid_max* =', umax_fit*coef, 'y =', y_fit, 'x_mid =', xmid

  open(unit=00,file=trim(settingsFile),status="unknown",position="append")
  write(00,'(A,1X,ES24.16E3,2X,A,1X,ES24.16E3,2X,A,1X,ES24.16E3)') &
       'u_mid_max* =', umax_fit*coef, 'y =', y_fit, 'x_mid =', xmid
  close(00)

  return
end subroutine RBcalc_umid_max
!===================================================================================================
! RBcalc_umid_max 结束: 计算 Rayleigh-Benard 工况下 x=1/2 处最大水平速度及其位置。
!===================================================================================================


!===================================================================================================
! 子程序: RBcalc_vmid_max
! 作用: 计算 Rayleigh-Benard 工况下 y=1/2 处最大垂直速度 vmax 及其 x 位置。
! 用途: 在 RayleighBenardCell 工况结束后的后处理中调用。
!===================================================================================================
subroutine RBcalc_vmid_max()
  use commondata
  implicit none
  integer(kind=4) :: i, k
  integer(kind=4) :: jMid, jB, jT
  integer(kind=4) :: i0
  integer(kind=4) :: ii(5)
  real(kind=8) :: vline(1:nx)
  real(kind=8) :: xk(5), fk(5)
  real(kind=8) :: vmax_fit, x_fit
  real(kind=8) :: coef, ymid
  real(kind=8) :: vmax_grid

  coef = velocityScaleCompare

  ! ---- (1) 取 y=1/2 中线剖面 v(x_i, y=1/2) ----
  if (mod(ny,2) == 1) then
    jMid = (ny + 1)/2
    ymid = yp(jMid)
    do i = 1, nx
      vline(i) = v(i,jMid)
    enddo
  else
    jB = ny/2
    jT = jB + 1
    ymid = 0.5d0*(yp(jB) + yp(jT))
    do i = 1, nx
      vline(i) = 0.5d0*(v(i,jB) + v(i,jT))
    enddo
  endif

  ! ---- (2) 找峰值所在网格点 i0 ----
  i0 = 1
  vmax_grid = vline(1)
  do i = 2, nx
    if (vline(i) > vmax_grid) then
      vmax_grid = vline(i)
      i0 = i
    endif
  enddo


  ! ---- (3) 取 5 点并拟合 ----
  if (i0 <= 2) then
    ii = (/ 1, 2, 3, 4, 5 /)
  elseif (i0 >= nx-1) then
    ii = (/ nx-4, nx-3, nx-2, nx-1, nx /)
  else
    ii = (/ i0-2, i0-1, i0, i0+1, i0+2 /)
  endif

  do k = 1, 5
    xk(k) = xp(ii(k))
    fk(k) = vline(ii(k))
  enddo

  call fit_parabola_ls5(xk, fk, +1, vmax_fit, x_fit)

  write(*,'(A,1X,ES24.16E3,2X,A,1X,ES24.16E3,2X,A,1X,ES24.16E3)') &
       'v_mid_max* =', vmax_fit*coef, 'x =', x_fit, 'y_mid =', ymid

  open(unit=00,file=trim(settingsFile),status="unknown",position="append")
  write(00,'(A,1X,ES24.16E3,2X,A,1X,ES24.16E3,2X,A,1X,ES24.16E3)') &
       'v_mid_max* =', vmax_fit*coef, 'x =', x_fit, 'y_mid =', ymid
  close(00)

  return
end subroutine RBcalc_vmid_max
!===================================================================================================







!===================================================================================================
! 子程序: calc_psi_vort_and_output
! 作用: 计算流函数、涡量，并输出相关诊断量。
! 用途: 仅在稳态计算结束且垂直边界无滑移时调用，非稳态不计算该项。
!===================================================================================================
subroutine calc_psi_vort_and_output()
  use commondata
  implicit none

  integer(kind=4) :: i, j, m
  real(kind=8) :: dx, dy, coef
  real(kind=8) :: u1, u2, u_mid, inc
  real(kind=8) :: dv_dx, du_dy
  real(kind=8) :: psi(nx,ny), vort(nx,ny)


  ! for fine-grid max(|psi|)  10001*10001
  real(kind=8) :: psi_abs_max, x_at_max, y_at_max, psi_center_abs_fine


  dx   = 1.0d0 / lengthUnit
  dy   = 1.0d0 / lengthUnit
  coef = velocityScaleCompare

  !=========================================================
  ! (A) 计算流函数 psi：  psi(x,y)=∫_0^y u(x,mu)dmu
  !     积分用“累积 Simpson”，中点 u(y=整点) 用二阶多项式插值
  !=========================================================

  do i = 1, nx

    ! ---- (A1) 第一个半格点 y=dy/2 的积分：用二次多项式并强制壁面 u=0（no-slip）
    ! 这里第一小段 [0,dy/2] 单独处理，整体阶数由它控制
    u1 = u(i,1) * coef
    u2 = u(i,2) * coef
    psi(i,1) = dy * (21.0d0*u1 - u2) / 72.0d0   !第一个点的psi

    ! ---- (A2) 从 y=(j-1/2)dy 到 y=(j+1/2)dy 的每一段长度为 dy，用 Simpson：
    ! ∫_{y_{j-1/2}}^{y_{j+1/2}} u dy ≈ dy/6 * [ u_{j-1/2} + 4 u_{j} + u_{j+1/2} ]
    do j = 2, ny
      m = j - 1   ! 这一段的中点在 y = m*dy（整点）

      ! ---- (A2.1) 计算中点速度 u(y=m*dy)：
      ! 用二阶 Lagrange 插值（三点），把半格点值插到整点
      !
      ! 对 m>=2：用 (m-1/2, m+1/2) 及更下方的 (m-3/2) 三点 => 系数(-1/8, 3/4, 3/8)
      ! 对 m=1（靠近底壁）：只能用最靠近底部的三点 => 系数( 3/8, 3/4,-1/8)
      if (m == 1) then
        u_mid = ( 3.0d0/8.0d0*u(i,1) + 3.0d0/4.0d0*u(i,2) - 1.0d0/8.0d0*u(i,3) ) * coef
      else
        u_mid = ( -1.0d0/8.0d0*u(i,m-1) + 3.0d0/4.0d0*u(i,m) + 3.0d0/8.0d0*u(i,m+1) ) * coef
      end if

      inc = dy/6.0d0 * ( (u(i,j-1)*coef) + 4.0d0*u_mid + (u(i,j)*coef) )
      psi(i,j) = psi(i,j-1) + inc
    end do

  end do


  call output_psi_center_abs(psi)     ! 基于粗网格局部四点插值得到中心点处的 abs(psi)

  !=========================================================
  ! (B) 计算涡量 vort = dv/dx - du/dy（2D）
  !     内部用中心差分；边界节点用“壁在半格距”假设下的二阶单边公式
  !=========================================================

  do j = 1, ny
    do i = 1, nx

      ! ---- dv/dx
      if (i == 1) then
        ! 在 x=dx/2 处，用 (wall, i=1, i=2) 二次拟合得到二阶近似
        dv_dx = ( 3.0d0*v(1,j) + v(2,j) - 4.0d0*0.0d0 ) / (3.0d0*dx)
      elseif (i == nx) then
        dv_dx = ( -3.0d0*v(nx,j) - v(nx-1,j) + 4.0d0*0.0d0 ) / (3.0d0*dx)
      else
        dv_dx = ( v(i+1,j) - v(i-1,j) ) / (2.0d0*dx)
      end if

      ! ---- du/dy
      if (j == 1) then
        du_dy = ( 3.0d0*u(i,1) + u(i,2) - 4.0d0*0.0d0 ) / (3.0d0*dy)
      elseif (j == ny) then
        du_dy = ( -3.0d0*u(i,ny) - u(i,ny-1) + 4.0d0*0.0d0 ) / (3.0d0*dy)
      else
        du_dy = ( u(i,j+1) - u(i,j-1) ) / (2.0d0*dy)
      end if

      ! psi 用了L/kappa；vort 也应在同一标度下
      ! vort = d(v*coef)/dx - d(u*coef)/dy = coef*(dv/dx - du/dy)
      vort(i,j) = coef * (dv_dx - du_dy)

    end do
  end do

  !=========================================================
  ! (C) 输出 Tecplot：X,Y,psi,vort
  !=========================================================
  call output_Tecplot_psi_vort(psi, vort)

  !=========================================================
  ! (D) 细网格 10001×10001：用三次样条插值寻找 max(|psi|)
  !=========================================================
  call calc_psi_absmax_fine_spline(psi, psi_abs_max, x_at_max, y_at_max, psi_center_abs_fine)

  write(*,'(a,1x,ES24.16E3)') "abs(psi_center_fine) =", real(psi_center_abs_fine,kind=8)

  write(*,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
       "max(|psi|) =", real(psi_abs_max,kind=8), "x* =", real(x_at_max,kind=8), &
       "y* =", real(y_at_max,kind=8)

  open(unit=00,file=trim(settingsFile),status="unknown",position="append")
  write(00,'(a,1x,ES24.16E3)') "abs(psi_center_fine) =", real(psi_center_abs_fine,kind=8)
  write(00,'(a,1x,ES24.16E3,2x,a,1x,ES24.16E3,2x,a,1x,ES24.16E3)') &
       "max(|psi|) =", real(psi_abs_max,kind=8), "x* =", real(x_at_max,kind=8), &
       "y* =", real(y_at_max,kind=8)
  close(00)

  return
end subroutine calc_psi_vort_and_output
!===================================================================================================
! calc_psi_vort_and_output 结束: 计算流函数、涡量，并输出相关诊断量。
!===================================================================================================



!===================================================================================================
! 子程序: output_Tecplot_psi_vort
! 作用: 把流函数和涡量场写出为 Tecplot 文件。
! 用途: 在 calc_psi_vort_and_output 中调用。
!===================================================================================================
subroutine output_Tecplot_psi_vort(psi, vort)
  use commondata
  implicit none
  real(kind=8), intent(in) :: psi(nx,ny), vort(nx,ny)

  integer(kind=4) :: i, j, k
  real(kind=4) :: zoneMarker, eohMarker
  character(len=40) :: title
  character(len=40) :: V1,V2,V3,V4
  integer(kind=4), parameter :: kmax=1
  character(len=40) :: zoneName
  character(len=100) :: filename

#ifdef steadyFlow
  write(filename,'(i12.12)') restartItcOffset+itc
#endif
#ifdef unsteadyFlow
  pltFileNum = pltFileNum+1
  write(filename,'(i12.12)') pltFileNum
#endif
  filename = adjustl(filename)

  open(41,file=trim(pltFolderPrefix)//"-psiVort-"//trim(filename)//'.plt', access='stream', form='unformatted')

  zoneMarker= 299.0
  eohMarker = 357.0

  write(41) "#!TDV101"
  write(41) 1

  title="psi-vort"
  call dumpstring(title)

  write(41) 4
  V1='X';     call dumpstring(V1)
  V2='Y';     call dumpstring(V2)
  V3='PSI';   call dumpstring(V3)
  V4='VORT';  call dumpstring(V4)

  write(41) zoneMarker
  zoneName="ZONE 001"
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

  ! 2 = Double
  write(41) 2
  write(41) 2
  write(41) 2
  write(41) 2

  write(41) 0
  write(41) -1

  do k=1,kmax
    do j=1,ny
      do i=1,nx
        write(41) real(xp(i),kind=8)
        write(41) real(yp(j),kind=8)
        write(41) real(psi(i,j),kind=8)
        write(41) real(vort(i,j),kind=8)
      end do
    end do
  end do

  close(41)
  return
end subroutine output_Tecplot_psi_vort
!===================================================================================================
! output_Tecplot_psi_vort 结束: 把流函数和涡量场写出为 Tecplot 文件。
!===================================================================================================


!===================================================================================================
! 子程序: calc_psi_absmax_fine_spline
! 作用: 用细网格样条插值搜索 abs(psi) 的最大值及位置。
! 用途: 在 calc_psi_vort_and_output 中调用。
!===================================================================================================
subroutine calc_psi_absmax_fine_spline(psi, psi_abs_max, x_at_max, y_at_max, psi_center_abs_fine)
  use commondata
  implicit none

  real(kind=8), intent(in)  :: psi(nx,ny)
  real(kind=8), intent(out) :: psi_abs_max, x_at_max, y_at_max, psi_center_abs_fine

  integer(kind=4), parameter :: nFinePerUnit = 10000
  integer(kind=4) :: i, j, k, l
  integer(kind=4) :: nFineX, nFineY
  integer(kind=4) :: nXExt, nYExt

  real(kind=8), allocatable :: xFine(:), yFine(:)
  real(kind=8), allocatable :: xExt(:), yExt(:)
  real(kind=8), allocatable :: row(:), y2row(:)
  real(kind=8), allocatable :: col(:), y2col(:)
  real(kind=8), allocatable :: psi_center_x(:)
  real(kind=8), allocatable :: psi_xfine(:,:)

  real(kind=8) :: val, xq, yq
  real(kind=8) :: xLen, yLen, xCenter, yCenter

  ! ---- 细网格坐标：0..1 等距
  ! Use the actual nondimensional domain size so the fine mesh also works for non-square cavities.
  xLen = xp(nx+1)
  yLen = yp(ny+1)
  xCenter = 0.5d0 * xLen
  yCenter = 0.5d0 * yLen
  nFineX = max(2, nint(xLen * dble(nFinePerUnit)) + 1)
  nFineY = max(2, nint(yLen * dble(nFinePerUnit)) + 1)

  allocate(xFine(nFineX), yFine(nFineY))
  do k = 1, nFineX
    xFine(k) = xLen * dble(k-1) / dble(nFineX-1)
  end do
  do l = 1, nFineY
    yFine(l) = yLen * dble(l-1) / dble(nFineY-1)
  end do

  ! ---- 为了能在 x=0/1 与 y=0/1 上插值，物理边界 psi=常数（可取0）补两个端点
  nXExt = nx + 2
  nYExt = ny + 2

  allocate(xExt(nXExt), yExt(nYExt))
  xExt(1)    = 0.0d0
  xExt(nXExt)= xp(nx+1)
  do i = 1, nx
    xExt(i+1) = xp(i)         !xExt = [0, xp(1),...,xp(nx), xLen]
  end do

  yExt(1)    = 0.0d0
  yExt(nYExt)= yp(ny+1)
  do j = 1, ny
    yExt(j+1) = yp(j)        !yExt = [0, yp(1),...,yp(ny), yLen]
  end do

  allocate(row(nXExt), y2row(nXExt))
  allocate(psi_center_x(ny))
  allocate(psi_xfine(nFineX, ny))

  ! ---- (1) 先对每个固定 y=yp(j) 的剖面做 x 方向三次样条，得到 psi(xFine, yp(j))
  do j = 1, ny
    row(1)  = 0.0d0
    row(nXExt) = 0.0d0
    do i = 1, nx
      row(i+1) = psi(i,j)     !row=[0, psi(1),...,psi(ny), 0]
    end do

    call spline_natural(nXExt, xExt, row, y2row)   !在row各节点处的二阶导

    call splint(nXExt, xExt, row, y2row, xCenter, psi_center_x(j))
    do k = 1, nFineX
      xq = xFine(k)
      call splint(nXExt, xExt, row, y2row, xq, val)   !在 10001 个细网格 x 点上采样
      psi_xfine(k,j) = val                            !在每个粗 y 层上，psi已经沿 x 被细化到 10001 个点。
    end do
  end do

  ! ---- (2) 再对每个固定 x=xFine(k) 的剖面做 y 方向三次样条，在 10001 个 yFine 上扫 max(|psi|)
  allocate(col(nYExt), y2col(nYExt))
  psi_abs_max = -1.0d0
  x_at_max    = 0.0d0
  y_at_max    = 0.0d0

  col(1)  = 0.0d0
  col(nYExt) = 0.0d0
  do j = 1, ny
    col(j+1) = psi_center_x(j)
  end do
  call spline_natural(nYExt, yExt, col, y2col)
  call splint(nYExt, yExt, col, y2col, yCenter, val)
  psi_center_abs_fine = dabs(val)

  do k = 1, nFineX
    col(1)  = 0.0d0
    col(nYExt) = 0.0d0
    do j = 1, ny
      col(j+1) = psi_xfine(k,j)        !对每一个固定的细网格 x=xFine(k)，有一条沿 y 的离散剖面数据 col
    end do

    call spline_natural(nYExt, yExt, col, y2col)

    do l = 1, nFineY
      yq = yFine(l)
      call splint(nYExt, yExt, col, y2col, yq, val)

      if (dabs(val) > psi_abs_max) then     !寻找最大的abs(psi)以及位置
        psi_abs_max = dabs(val)
        x_at_max = xFine(k)
        y_at_max = yFine(l)
      end if
    end do
  end do

  deallocate(xFine, yFine, xExt, yExt, row, y2row, col, y2col, psi_center_x, psi_xfine)
  return
end subroutine calc_psi_absmax_fine_spline
!===================================================================================================
! calc_psi_absmax_fine_spline 结束: 用细网格样条插值搜索 abs(psi) 的最大值及位置。
!===================================================================================================


!===================================================================================================
! 子程序: spline_natural
! 作用: 构造自然三次样条所需的二阶导数。
! 用途: 在 calc_psi_absmax_fine_spline 中作为样条预处理调用。
!===================================================================================================
subroutine spline_natural(n, x, y, y2)
  implicit none
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in)    :: x(*), y(*)
  real(kind=8), intent(out)   :: y2(*)

  integer(kind=4) :: i, k
  real(kind=8), allocatable :: u(:)
  real(kind=8) :: sig, p

  allocate(u(n))

  ! natural spline: y2(1)=0, y2(n)=0
  y2(1) = 0.0d0
  u(1)  = 0.0d0

  do i = 2, n-1
    sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
    p = sig*y2(i-1) + 2.0d0
    y2(i) = (sig - 1.0d0) / p
    u(i) = ( 6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))/(x(i)-x(i-1))) / (x(i+1)-x(i-1)) - sig*u(i-1) ) / p
  end do

  y2(n) = 0.0d0

  do k = n-1, 1, -1
    y2(k) = y2(k)*y2(k+1) + u(k)
  end do

  deallocate(u)
  return
end subroutine spline_natural
!===================================================================================================
! spline_natural 结束: 构造自然三次样条所需的二阶导数。
!===================================================================================================


!===================================================================================================
! 子程序: splint
! 作用: 根据样条系数在查询点上进行插值。
! 用途: 在 calc_psi_absmax_fine_spline 中重复调用。
!===================================================================================================
subroutine splint(n, xa, ya, y2a, x, y)
  implicit none
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in)    :: xa(*), ya(*), y2a(*), x
  real(kind=8), intent(out)   :: y

  integer(kind=4) :: klo, khi, k
  real(kind=8) :: h, a, b

  klo = 1
  khi = n

  do while (khi - klo > 1)
    k = (khi + klo)/2
    if (xa(k) > x) then
      khi = k
    else
      klo = k
    end if
  end do

  h = xa(khi) - xa(klo)
  if (h == 0.0d0) then
    y = ya(klo)
    return
  end if

  a = (xa(khi) - x)/h
  b = (x - xa(klo))/h

  y = a*ya(klo) + b*ya(khi) + ( (a*a*a - a)*y2a(klo) + (b*b*b - b)*y2a(khi) ) * (h*h)/6.0d0
  return
end subroutine splint
!===================================================================================================
! splint 结束: 根据样条系数在查询点上进行插值。
!===================================================================================================


!===================================================================================================
! 子程序: output_psi_center_abs
! 作用: 输出腔体中心位置的 abs(psi) 诊断结果。
! 用途: 在 calc_psi_vort_and_output 中调用。
!===================================================================================================
subroutine output_psi_center_abs(psi)
  use commondata
  implicit none
  real(kind=8), intent(in) :: psi(nx,ny)

  integer(kind=4) :: i0, j0, p, q
  integer(kind=4) :: ii(4), jj(4)
  real(kind=8) :: x0, y0, psi_center, psi_center_abs
  real(kind=8) :: x4(4), y4(4), f4(4), gx(4)

  x0 = 0.5d0 * xp(nx+1)
  y0 = 0.5d0 * yp(ny+1)

  if (nx < 4 .or. ny < 4) then
    psi_center = psi((nx+1)/2, (ny+1)/2)
  else
    i0 = 1
    do while (i0 < nx .and. xp(i0+1) <= x0)
      i0 = i0 + 1
    end do
    i0 = max(1, min(i0-1, nx-3))
    do p = 1, 4
      ii(p) = i0 + p - 1
      x4(p) = xp(ii(p))
    end do

    j0 = 1
    do while (j0 < ny .and. yp(j0+1) <= y0)
      j0 = j0 + 1
    end do
    j0 = max(1, min(j0-1, ny-3))
    do q = 1, 4
      jj(q) = j0 + q - 1
      y4(q) = yp(jj(q))
    end do

    do q = 1, 4
      do p = 1, 4
        f4(p) = psi(ii(p), jj(q))
      end do
      call interp_lagrange_4(x0, x4, f4, gx(q))
    end do

    call interp_lagrange_4(y0, y4, gx, psi_center)
  endif

  psi_center_abs = dabs(psi_center)

  ! Screen output
  write(*,'(a,1x,ES24.16E3)') "abs(psi_center_coarse) =", real(psi_center_abs,kind=8)

  ! Log output
  open(unit=00,file=trim(settingsFile),status="unknown",position="append")
  write(00,'(a,1x,ES24.16E3)') "abs(psi_center_coarse) =", real(psi_center_abs,kind=8)
  close(00)

  return
contains
  ! 子程序: interp_lagrange_4
  ! 作用: 对给定的四个节点执行四点 Lagrange 插值。
  subroutine interp_lagrange_4(xq, xk, fk, fq)
    implicit none
    real(kind=8), intent(in)  :: xq
    real(kind=8), intent(in)  :: xk(4), fk(4)
    real(kind=8), intent(out) :: fq
    integer(kind=4) :: a, b
    real(kind=8) :: basis

    fq = 0.0d0
    do a = 1, 4
      basis = 1.0d0
      do b = 1, 4
        if (b /= a) basis = basis * (xq - xk(b)) / (xk(a) - xk(b))
      end do
      fq = fq + fk(a) * basis
    end do
  end subroutine interp_lagrange_4
end subroutine output_psi_center_abs
!===================================================================================================
! output_psi_center_abs 结束: 输出腔体中心位置的 abs(psi) 诊断结果。
!===================================================================================================
