!=============================================================
!!!    二维均匀网格混合 RLB-4：高阶正则化流场 + RK4 宏观温度
!!!    流场 D2Q9，温度中心二阶空间差分，单 GPU OpenACC
!!!    依据: Nee (2025)，式 (6)--(8)、(13)--(16) 及第 2.4 节。
!!!    主框架参照 2DRBOpenmp，GPU 数据与输出参照 2DRBOpenacc。
!=============================================================

!=============================================================
!   手动选择流动模式：侧热通常用 steadyFlow，RB 通常用 unsteadyFlow。
#define steadyFlow
!#define unsteadyFlow

!   手动选择算例：两项只能启用一项；边界由算例自动配套。
#define SideHeatedCell
!#define RayleighBenardCell

#if defined(steadyFlow) && defined(unsteadyFlow)
#error "Choose only one flow mode: steadyFlow or unsteadyFlow"
#endif
#if !defined(steadyFlow) && !defined(unsteadyFlow)
#error "Define one flow mode: steadyFlow or unsteadyFlow"
#endif
#if defined(SideHeatedCell) && defined(RayleighBenardCell)
#error "Choose only one convection case: SideHeatedCell or RayleighBenardCell"
#endif
#if !defined(SideHeatedCell) && !defined(RayleighBenardCell)
#error "Define one convection case: SideHeatedCell or RayleighBenardCell"
#endif

!   沿用本项目封闭腔体：四壁无滑移，RB 侧壁绝热，不增加周期几何。
#define HorizontalWallsNoslip
#define VerticalWallsNoslip
#ifdef SideHeatedCell
#define VerticalWallsConstT
#define HorizontalWallsAdiabatic
#else
#define HorizontalWallsConstT
#define VerticalWallsAdiabatic
#endif
!   原文未完整给出温度空间模板；补齐选择和耦合精度见本目录说明。
!=============================================================

!   全局模块
    module commondata
        implicit none
        !===============================================================================================
        ! 是否在计算前从旧算例重启
        integer(kind=4), parameter :: loadInitField=0   ! 0: 不重启；1: 按 latest .meta 自动续算

        ! 正常断电续算只需要设置 loadInitField=1；
        ! 代码只读取 <reloadFilePrefix>-latest.meta，并从里面找到最新的 .bin。
        ! 检查点二进制、meta 与 Nu/Re 历史须配套保留；缺少 meta 时直接停止。
        integer(kind=4) :: reloadFileNum=0              ! 从 latest .meta 恢复的重启文件编号
        !===============================================================================================
        real(kind=8) :: reloadDimensionlessTime=0.0d0   ! 续算前已累计的 t_ff；从 latest .meta 读取
        integer(kind=4) :: restartItcOffset=0           ! 续算前已累计的格子步数；从 latest .meta 读取
        logical :: reloadMetadataLoaded=.false.         ! 标记是否成功读取 reload 元数据文件
        !===============================================================================================

        !===============================================================================================
        ! 无量纲参数
        integer(kind=4), parameter :: nx=256, ny=256     ! 直接修改二维均匀网格数
#ifdef SideHeatedCell
        real(kind=8), parameter :: lengthUnit=dble(nx)     !侧壁差温：特征长度取 x 方向长度
#else
        real(kind=8), parameter :: lengthUnit=dble(ny)     !上下差温：特征长度取 y 方向长度
#endif
        real(kind=8), parameter :: pi = acos(-1.0d0)

        real(kind=8), parameter :: Rayleigh=1.0d6       ! 直接修改瑞利数
        real(kind=8), parameter :: Prandtl=0.71d0
        real(kind=8), parameter :: Mach=dsqrt(3.0d0)*0.01d0  ! 默认 Uff=0.01；Ma=Uff/cs
        real(kind=8), parameter :: Thot=1.0d0, Tcold=0.0d0  ! 温度以冷热壁温差归一化
        real(kind=8), parameter :: Tref=0.5d0*(Thot+Tcold)
        real(kind=8), parameter :: tauf=0.5d0+Mach*lengthUnit*dsqrt(3.0d0*Prandtl/Rayleigh) 
        real(kind=8), parameter :: viscosity=(tauf-0.5d0)/3.0d0
        real(kind=8), parameter :: diffusivity=viscosity/Prandtl
        

        ! velocityScaleCompare is used only in velocity-related post-processing to convert lattice velocity
        ! to the nondimensional velocity scale adopted by the reference paper being compared.
        real(kind=8), parameter :: velocityScaleCompare=lengthUnit/diffusivity   ! 默认采用热扩散标度 UL/kappa；若要按自由落体标度比较，可改为 1.0d0/velocityUnit
        
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1



        real(kind=8), parameter :: temperatureDifference=Thot-Tcold  ! 正的冷热壁温差



        ! 浮力项参数
        real(kind=8), parameter :: gBeta1=Rayleigh*viscosity*diffusivity/lengthUnit
        real(kind=8), parameter :: gBeta=gBeta1/lengthUnit/lengthUnit/temperatureDifference             !gbetaΔT
        
        real(kind=8), parameter :: timeUnit=dsqrt(lengthUnit/(gBeta*temperatureDifference))      !无量纲时间
        real(kind=8), parameter :: velocityUnit=dsqrt(gBeta*temperatureDifference*lengthUnit)  !无量纲速度
    
        ! 流场按原文 raw 动量与显式力更新；不能混入 Guo 半步力修正。
        real(kind=8), parameter :: cs2=1.0d0/3.0d0
        real(kind=8), parameter :: rlbRelaxationRate=1.0d0/tauf
        integer(kind=4), parameter :: thermalSubsteps=1         ! 每个格子步内的温度 RK4 子步数
        integer(kind=4), parameter :: steadyCheckIntervalItc=2000 ! 稳态误差/有限性检查间隔
        real(kind=8), parameter :: rbInitialPerturbation=1.0d-3  ! RB 初始扰动，相对于冷热壁温差
        character(len=32), parameter :: reloadMagic='HYBRID_RLB4_D2Q9_T_V1'
        integer(kind=4) :: sampleCountOffset=0                 ! 已保存 Nu/Re 历史样本数
        !===============================================================================================          
        
        !===============================================================================================
        ! 输出/备份相关设置（以自由落体时间 t_ff 为单位）
        real(kind=8), parameter :: epsU=1.0d-7, epsT=1.0d-7    ! 稳态收敛阈值   

#ifdef steadyFlow
        real(kind=8), parameter :: outputSnapshotInterval=10.0d0   ! 快照和 Nu/Re 时间序列采样间隔（单位：t_ff）
        real(kind=8), parameter :: reloadFileInterval=100.0d0  ! f/T 重启文件输出间隔（单位：t_ff）
        real(kind=8), parameter :: outputPltFileInterval=100.0d0  ! Tecplot 文件输出间隔（单位：t_ff）
        integer(kind=4), parameter :: dimensionlessTimeMax=int(12000.0d0/outputSnapshotInterval)
        integer(kind=4), parameter :: outputSnapshotFile=1   ! 是否输出后处理快照文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputPltFile=1   ! 是否输出 plt 文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputReloadFile=1 ! 是否周期输出 f/T 重启文件：0=不输出，1=输出
        integer(kind=4), parameter :: itc_max=20000000  ! 稳态：最大格子步，实际仍由 errorU/errorT 提前停止
#endif

#ifdef unsteadyFlow
        real(kind=8), parameter :: outputSnapshotInterval=0.5d0   ! uvTrho 快照和 Nu/Re 时间序列采样间隔（单位：t_ff）
        real(kind=8), parameter :: reloadFileInterval=100.0d0  ! f/T 重启文件输出间隔（单位：t_ff）
        real(kind=8), parameter :: outputPltFileInterval=100.0d0  ! Tecplot 文件周期输出间隔（单位：t_ff）
        real(kind=8), parameter :: unsteadyRunDuration=1000.0d0  ! 非稳态总目标时长，续算时只补足到该 t_ff
        ! 以下三个参数只控制非稳态结束后的 Nu/Re 统计平均窗口，不改变推进时长或采样频率。
        ! 时间以整个非稳态算例的绝对 t_ff 计，续算统计会包含旧文件中已有的历史数据。
        real(kind=8), parameter :: unsteadyAverageStartTf=0.5d0*unsteadyRunDuration  ! 平均窗口起点
        real(kind=8), parameter :: unsteadyAverageEndTf=unsteadyRunDuration          ! 平均窗口终点
        real(kind=8), parameter :: unsteadyAverageMidTf=0.5d0*(unsteadyAverageStartTf+unsteadyAverageEndTf) ! 前/后半分界
        integer(kind=4), parameter :: unsteadySampleCount=max(0, int(unsteadyRunDuration/outputSnapshotInterval+1.0d-12))   !计数器，输出多少次快照
        integer(kind=4), parameter :: dimensionlessTimeMax=max(1,unsteadySampleCount)
        integer(kind=4), parameter :: outputSnapshotFile=1   ! 是否输出后处理快照文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputPltFile=1   ! 是否输出 plt 文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputReloadFile=1 ! 是否周期输出 f/T 重启文件：0=不输出，1=输出
        integer(kind=4), parameter :: itc_max=max(1, int(unsteadyRunDuration*timeUnit+0.5d0)) ! 非稳态：由总 t_ff 自动换算格子步
#endif

        integer(kind=4) :: snapshotFileNum, pltFileNum  ! 快照/plt 输出文件的计数器
        ! 每次调用对应输出子程序时递增（用于文件名编号）

        integer(kind=4) :: dimensionlessTime
        integer(kind=4) :: outputSnapshotIntervalItc
        integer(kind=4) :: reloadFileIntervalItc, outputPltFileIntervalItc
        ! 统计/输出时间点编号（与 outputSnapshotInterval 对应）：
        ! 每调用一次 calNuRe() 就 dimensionlessTime = dimensionlessTime + 1
        ! 用于索引 NuVolAvg/ReVolAvg 数组，并用于输出的时间轴：t = (sampleCountOffset+dimensionlessTime)*outputSnapshotInterval（单位：t_ff）

        real(kind=8) :: NuVolAvg(0:dimensionlessTimeMax), ReVolAvg(0:dimensionlessTimeMax)
        ! 体平均 Nu 和 Re 的时间序列缓存
        ! 只有在启用并调用 calNuRe() 的情况下这些数组才会被真正填充
  
        character(len=100) :: snapshotFilePrefix="buoyancyCavity2DOpenaccHybridRLB4Snapshot"
        ! 快照输出文件前缀（实际文件名形如：<snapshotFilePrefix>-<编号>.bin）

        character(len=100) :: pltFolderPrefix="buoyancyCavity2DOpenaccHybridRLB4Tecplot"
        ! plt 输出文件前缀（实际文件名形如：<pltFolderPrefix>-<编号>.plt）

        character(len=100) :: reloadFilePrefix="reloadFile2DOpenaccHybridRLB4"
        ! 重启读取文件的前缀；latest meta 模式实际读取 meta 中记录的 <reloadFilePrefix>-<编号>.bin
        
        character(len=100) :: settingsFile="SimulationSettings2DOpenaccHybridRLB4.txt"
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
        ! RK4 工作数组只作临时量；整步状态为 f 和 T，检查点不保存工作数组。
        real(kind=8), allocatable :: temperatureBase(:,:), temperatureStage(:,:)
        real(kind=8), allocatable :: temperatureRhs(:,:), temperatureIncrement(:,:)
        real(kind=8), allocatable :: Fx(:,:), Fy(:,:)


        integer(kind=4) :: itc

#ifdef steadyFlow
        real(kind=8) :: Nu_global, Nu_hot, Nu_cold, Nu_middle    !平均Nu，全场，侧壁以及中线
        real(kind=8) :: Nu_hot_max, Nu_hot_min, Nu_hot_max_position, Nu_hot_min_position    !左侧壁面的最大最小Nu，以及对应的位置
#endif
        
        
        !格子离散速度和权重
        integer(kind=4) :: ex(0:8), ey(0:8)
        data ex/0, 1, 0, -1,  0, 1, -1, -1,  1/
        data ey/0, 0, 1,  0, -1, 1,  1, -1, -1/
        real(kind=8) :: omega(0:8)
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
#ifdef unsteadyFlow
    integer(kind=4) :: nextSampleItc
    integer(kind=4) :: nextSampleAbsItc
    integer(kind=4) :: unsteadyItcRemaining
#endif
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
    ! 目标已经达到时不重复写同一时刻的稳态历史；续算前应提高目标总步数/时长。
    if((loadInitField.EQ.1).AND.(restartItcOffset.GE.itc_max)) then
        write(*,'(A)') 'Checkpoint already reaches the requested horizon; no new steps were taken.'
        stop
    endif
#ifdef steadyFlow
    ! 已收敛检查点不重复追加同一时刻历史；收紧 eps 后仍可继续推进。
    if((loadInitField.EQ.1).AND.(errorU.LE.epsU).AND.(errorT.LE.epsT)) then
        write(*,'(A)') 'Restart state already satisfies steady criteria; no new steps were taken.'
        stop
    endif
#endif
    call enter_data_2d_openacc()
#ifdef unsteadyFlow
    ! 非稳态的 itc_max 是整个算例的总目标步数；
    ! 续算时 restartItcOffset 是旧算例已经完成的步数，本次只推进剩余步数。
    unsteadyItcRemaining = max(0, itc_max - restartItcOffset)
#endif

    !===============================================================================================

    call CPU_TIME(timeStart)         !当前进程累计消耗的 CPU 时间,包括并行
    ! system_clock 返回墙钟计数器和每秒计数率；
    ! 下面用 counter/rate 把它换算成实际经过的秒数。
    call system_clock(wallClockStart, wallClockRate)
    timeStart2 = dble(wallClockStart) / dble(max(wallClockRate,1_8))
#ifdef steadyFlow
    do while( ((errorU.GT.epsU).OR.(errorT.GT.epsT)).AND.(restartItcOffset+itc.LT.itc_max) )
#endif
#ifdef unsteadyFlow
    do while( itc.LT.unsteadyItcRemaining )   !非稳态：续算时只推进到 unsteadyRunDuration 对应的总格子步
#endif

        itc = itc+1
        
        ! 按论文流程先推进流场，再以本步恢复的速度推进温度。
        ! RK4 四阶段冻结同一速度；这不代表整个流热耦合为四阶时间精度。
        call collision()
        call streaming()
        call bounceback()
        call macro()
        call advance_temperature_rk4(1.0d0)
        if(mod(restartItcOffset+itc,steadyCheckIntervalItc).EQ.0) call check_fields_hybrid()

#ifdef steadyFlow
        ! 周期输出按累计格子步判断；否则从 1050tf 续算会在 1150tf 才输出，
        ! 而不是接回不断电运行应有的 1100tf、1200tf、...
        if(MOD(restartItcOffset+itc,steadyCheckIntervalItc).EQ.0) call check()
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
            call output_ReloadFile()      !稳态模式下的可选周期 f/T 重启文件输出
        endif
#endif

#ifdef unsteadyFlow
        do while(sampleCountOffset+dimensionlessTime.LT.unsteadySampleCount)   !判断本次还要不要继续采样
            ! 每个目标采样时刻都按绝对 t_ff 换算到本次运行段的 itc，续算时不会重复旧样本。
            nextSampleAbsItc = max(1, int( &
                real(sampleCountOffset+dimensionlessTime+1,kind=8)*outputSnapshotInterval*timeUnit+0.5d0)) !算下一个采样点对应的绝对格子步
            nextSampleItc = max(1, nextSampleAbsItc - restartItcOffset) !换成本次续算运行里的局部步数
            if(itc.LT.nextSampleItc) exit !如果本次运行还没走到这个采样步，就先退出采样循环，继续推进 LBM；走到了就调用 calNuRe()
            call calNuRe()
            if(outputSnapshotFile.EQ.1) then
                call update_host_snapshot_2d_openacc()
                call output_SnapshotFile()          !每 0.5 t_ff 输出一次后处理 uvTrho 快照
            endif
        enddo
        if( (outputPltFile.EQ.1).AND.(MOD(restartItcOffset+itc, outputPltFileIntervalItc).EQ.0) ) then
            call update_host_tecplot_2d_openacc()
            call output_Tecplot()  !非稳态模式下的可选周期 Tecplot 输出
        endif
        if( (outputReloadFile.EQ.1).AND.(MOD(restartItcOffset+itc, reloadFileIntervalItc).EQ.0) ) then
            call update_host_reload_2d_openacc()
            call output_ReloadFile()      !非稳态模式下的可选周期 f/T 重启文件输出
        endif
#endif
     enddo

    !$acc wait(1)
    call CPU_TIME(timeEnd)         !当前进程累计消耗的 CPU 时间,包括并行
    ! 取墙钟结束计数器，用于后面输出 OpenACC 实际耗时。
    call system_clock(wallClockEnd, wallClockRate)
    timeEnd2 = dble(wallClockEnd) / dble(max(wallClockRate,1_8))
    call check_fields_hybrid()
    call update_host_snapshot_2d_openacc()
    call output_thermal_budget()

#ifdef steadyFlow
    if(outputPltFile.EQ.1) call output_Tecplot() ! 输出最后一步
    if(outputSnapshotFile.EQ.1) call output_SnapshotFile()
#endif

#ifdef unsteadyFlow
    call output_unsteady_NuRe_postprocess()
#endif

    !===============================================================================================



    !===============================================================================================

    
#ifdef steadyFlow
! 稳态最终标量诊断：
! 1) 在 steadyFlow 结束时输出；是否收敛要另看 errorU/errorT
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


#ifdef VerticalWallsNoslip
    ! psi/vort 后处理默认封闭腔体：四周无滑移，psi 在物理边界取同一常数。
    ! 若垂直边界改为周期速度边界，流函数边界补点和涡量单边差分需要另写周期版本。
    call calc_psi_vort_and_output()  ! 输出中心abs(psi), max(abs(psi))及位置；max位置用细网格样条插值
#endif

#ifdef steadyFlow
    call calNuRe()
#endif

    ! 最后一个检查点在最终宏观输出/计数之后保存，恢复时计数和 f/T 属于同一步。
    if(outputReloadFile.EQ.1) then
        call update_host_reload_2d_openacc()
        call output_ReloadFile()
    endif

     

    open(unit=00,file=trim(settingsFile),status='unknown',position='append')        !在这个txt文件后面继续写（追加模式）
    write(00,*) "======================================================================"
    write(00,*) "Time (CPU) = ", real(timeEnd-timeStart,kind=8), "s"                             !当前进程累计消耗的 CPU 时间,包括并行
    write(00,*) "MLUPS = ", real( dble(nx)*dble(ny)*dble(itc)/max(timeEnd-timeStart,epsilon(1.0d0))/1.0d6,kind=8 )   !百万格点更新/秒
    write(00,*) "Time (ACC) = ", real(timeEnd2-timeStart2,kind=8), "s"                           !墙钟时间
    write(00,*) "MLUPS (ACC) = ", real( dble(nx)*dble(ny)*dble(itc)/max(timeEnd2-timeStart2,epsilon(1.0d0))/1.0d6,kind=8 )   !百万格点更新/秒
#ifdef steadyFlow
    if((errorU.LE.epsU).AND.(errorT.LE.epsT)) then
        write(00,*) 'Steady criteria satisfied.'
    else
        write(00,*) 'Step limit reached before steady criteria; final diagnostics are not converged benchmarks.'
    endif
    write(00,'(a,1x,ES24.16E3)') "Nu_global =", real(Nu_global,kind=8)
    write(00,'(a,1x,ES24.16E3)') "Nu_hot    =", real(Nu_hot,kind=8)
    write(00,'(a,1x,ES24.16E3)') "Nu_cold   =", real(Nu_cold,kind=8)
#endif
    write(00,*) "HybridRLB4: raw-momentum RLB flow; central-space RK4 temperature; explicit sequential coupling."



    write(00,*) "Dellocate Array......"
    call exit_data_2d_openacc()
    deallocate(f)
    deallocate(f_post)
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
    deallocate(temperatureBase,temperatureStage,temperatureRhs,temperatureIncrement)
    write(00,*) "    "
    
    write(00,*) "Time integration completed; benchmark validation is separate."
    
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
    real(kind=8) :: feqLocal(0:8)
    real(kind=8) :: xLen, yLen, rbInitPerturbAmp
    character(len=100) :: reloadFileName
    character(len=32) :: binaryMagic
    integer(kind=4) :: binaryNx, binaryNy, binaryQf, binaryQg
    

    itc = 0
    sampleCountOffset = 0
    errorU = 100.0d0
    errorT = 100.0d0 
    snapshotFileNum = 0
    pltFileNum = 0
    restartItcOffset = 0
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
    

    write(00,*)"-------------------------------------------------------------------------------"
    write(00,*) 'Mesh:',nx,ny
    write(00,*) 'Rayleigh=',real(Rayleigh,kind=8), '; Prandtl =',real(Prandtl,kind=8), '; Mach =',real(Mach,kind=8)
    write(00,*) "Length unit: L0 =", real(lengthUnit,kind=8)
    write(00,*) "Time unit: Sqrt(L0/(gBeta*DeltaT)) =", real(timeUnit,kind=8)
    write(00,*) "Velocity unit: Sqrt(gBeta*L0*DeltaT) =", real(velocityUnit,kind=8)
    write(00,*) "   "
    write(00,*) 'tauf=',real(tauf,kind=8)
    write(00,*) "thermalScheme = central2-space RK4; flowScheme = high-order regularized D2Q9"
    write(00,*) 'thermalSubsteps=',thermalSubsteps
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
    write(00,*) "unsteadyRunDuration =", real(unsteadyRunDuration,kind=8), "free-fall time units"
    write(00,*) "unsteadySampleCount =", unsteadySampleCount
#endif
    if(loadInitField.EQ.1) then
        write(00,*) "Restart requires matching binary data, metadata and Nu/Re histories."
    endif
    write(00,*) "itc_max =",itc_max
    write(00,*) "default epsU =", real(epsU,kind=8),"; epsT =", real(epsT,kind=8)
    write(00,*) "HybridRLB4: raw-momentum RLB flow; central-space RK4 temperature; explicit sequential coupling."
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
    allocate(temperatureBase(nx,ny),temperatureRhs(nx,ny),temperatureIncrement(nx,ny))
    allocate(temperatureStage(0:nx+1,0:ny+1))

    allocate (Fx(nx,ny))
    allocate (Fy(nx,ny))


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
    
    if(thermalSubsteps.LT.1) error stop 'thermalSubsteps must be positive'
    if((nx.LT.6).OR.(ny.LT.6)) error stop 'Use at least 6 cells per direction for wall diagnostics'
    ! 二次恒温 ghost 的扩散谱上界高于标准五点内区；取 10*kappa 作保守界。
    ! 经典 RK4 负实轴范围约 2.785；此检查不保证带对流项时也稳定。
    if(10.0d0*diffusivity/dble(thermalSubsteps).GT.2.7d0) &
        error stop 'Increase thermalSubsteps to satisfy the RK4 diffusion stability limit' 

    if(loadInitField.EQ.0) then                    !在不加载文件的情况下，都是零场为初值
    
        u = 0.0d0
        v = 0.0d0
        T = 0.0d0
        
        write(00,*) "Initial field is set exactly"
        if(reloadDimensionlessTime.NE.0.0d0) then        !在不加载文件的情况下，reloadDimensionlessTime必须是零
            write(00,*) "Error: since loadInitField.EQ.0, reloadDimensionlessTime should also be 0"
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
#ifdef RayleighBenardCell
    ! 非稳态 RB 从完全对称的导热态启动时需要一个确定性小扰动来触发对流。
    ! 扰动幅值和模态与同组 LBMCDE Rayleigh-Benard 非稳态算例保持一致。
    xLen = xp(nx+1)
    yLen = yp(ny+1)
    rbInitPerturbAmp = rbInitialPerturbation*(Thot-Tcold)
    do i = 1, nx
        do j = 1, ny
            T(i,j) = T(i,j) + rbInitPerturbAmp * dsin(2.0d0*pi*xp(i)/xLen) * dsin(pi*yp(j)/yLen)
        enddo
    enddo
    write(00,'(a,1x,es12.4)') "RB initial T perturbation amplitude =", rbInitPerturbAmp
#endif
    write(00,*) "Temperature B.C. for horizontal walls are:===Hot/cold wall==="
#endif

#ifdef VerticalWallsAdiabatic
    write(00,*) "Temperature B.C. for vertical walls are:===Adiabatic wall==="
#endif

#ifdef HorizontalWallsAdiabatic
    write(00,*) "Temperature B.C. for horizontal walls are:===Adiabatic wall==="
#endif

        f = 0.0d0
        do j = 1,ny
            do i = 1,nx
                Fx(i,j) = 0.0d0
                Fy(i,j) = rho(i,j)*gBeta*(T(i,j)-Tref)
                call rlb_equilibrium(rho(i,j),u(i,j),v(i,j),feqLocal)
                f(i,j,0:8) = feqLocal
            enddo
        enddo

    elseif(loadInitField.EQ.1) then
        call read_reload_metadata(reloadFileName)
        open(unit=01,file=trim(reloadFilePrefix)//'-'//trim(reloadFileName)//'.bin', &
             form='unformatted',access='sequential',status='old')
        read(01) binaryMagic, binaryNx, binaryNy, binaryQf, binaryQg
        if ((binaryMagic.NE.reloadMagic).OR.(binaryNx.NE.nx).OR.(binaryNy.NE.ny) &
            .OR.(binaryQf.NE.9).OR.(binaryQg.NE.0)) error stop 'HybridRLB4 restart header mismatch'
        read(01) (((f(i,j,alpha),i=1,nx),j=1,ny),alpha=0,8)
        read(01) T
#ifdef steadyFlow
        read(01) up,vp,Tp,errorU,errorT
#endif
        close(01)
        call reconstruct_macro_from_fT()
        write(00,*) 'HybridRLB4 restart loaded; accumulated steps =',restartItcOffset
    else
        write(00,*) "Error: initial field is not properly set"  !如果 loadInitField 不是 0/1 或逻辑不一致，直接停止
        stop
    endif
    
    write(00,*)"-------------------------------------------------------------------------------"
close(00)
    
    f_post = 0.0d0
    temperatureBase=0.0d0
    temperatureStage=0.0d0
    temperatureRhs=0.0d0
    temperatureIncrement=0.0d0
        
    if(loadInitField.EQ.0) then
        snapshotFileNum = 0
        pltFileNum = 0
        reloadFileNum = 0
        restartItcOffset = 0
        reloadDimensionlessTime = 0.0d0
#ifdef steadyFlow
        ! 新算例第一段收敛误差应从初始场开始比较。
        up = u
        vp = v
        Tp = T
#endif
    endif
    dimensionlessTime = 0
    !新算例：清零，开始记录新历史。
    !续算：也清零，但不是丢旧历史；旧历史在 .dat 文件里，新数组只记录本次续算段。
    NuVolAvg = 0.0d0
    ReVolAvg = 0.0d0
    
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

    !$acc enter data copyin(xp,yp,ex,ey,omega)
    !$acc enter data create(f_post,temperatureBase,temperatureStage,temperatureRhs,temperatureIncrement)
    !$acc enter data copyin(u,v,T,rho,f,Fx,Fy)
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
! 作用: 同步整步状态 f/T，供 HybridRLB4 严格检查点保存。
!===================================================================================================
  subroutine update_host_reload_2d_openacc()
    use commondata
    implicit none

    !$acc wait(1)
    !$acc update self(f,T)
#ifdef steadyFlow
    !$acc update self(up,vp,Tp)
#endif
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
    !$acc exit data delete(f_post,u,v,T,rho,f,Fx,Fy)
    !$acc exit data delete(temperatureBase,temperatureStage,temperatureRhs,temperatureIncrement)
    !$acc exit data delete(xp,yp,ex,ey,omega)
  end subroutine exit_data_2d_openacc
!===================================================================================================



!===================================================================================================
! 子程序: collision
! 作用: 完成流场分布函数 f 的碰撞更新，并处理体力项离散修正。
! 用途: 在主程序时间推进循环中调用，位于 streaming 之前。
!===================================================================================================

subroutine collision()
    use commondata
    implicit none
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: pop(0:8), popPost(0:8)
    real(kind=8) :: rhoLocal, ux, uy, forceX, forceY

    !$acc parallel loop gang vector collapse(2) default(none) present(f,f_post,rho,u,v,T,Fx,Fy) async(1) &
    !$acc& private(alpha,pop,popPost,rhoLocal,ux,uy,forceX,forceY)
    do j=1,ny
        do i=1,nx
            rhoLocal=rho(i,j)
            ux=u(i,j)
            uy=v(i,j)
            forceX=0.0d0
            forceY=rhoLocal*gBeta*(T(i,j)-Tref)
            Fx(i,j)=forceX
            Fy(i,j)=forceY
            do alpha=0,8
                pop(alpha)=f(i,j,alpha)
            enddo
            call rlb_collide_node(pop,rhoLocal,ux,uy,forceX,forceY,popPost)
            do alpha=0,8
                f_post(i,j,alpha)=popPost(alpha)
            enddo
        enddo
    enddo
end subroutine collision


!===================================================================================================
! collision 结束: 流场碰撞步骤完成。
!===================================================================================================


!===================================================================================================
! 子程序: streaming
! 作用: 完成流场分布函数 f 的迁移，把碰撞后的信息传播到相邻格点。
! 用途: 在主程序时间推进循环中调用，位于 collision 之后、bounceback 之前。
!===================================================================================================

subroutine streaming()
    use commondata
    implicit none
    integer(kind=4) :: i,j,alpha,ip,jp
    !$acc parallel loop gang vector collapse(2) default(none) present(f,f_post,ex,ey) async(1) &
    !$acc& private(alpha,ip,jp)
    do j=1,ny
      do i=1,nx
        do alpha=0,8
          ip=i-ex(alpha)
          jp=j-ey(alpha)
          ! 域外链交给边界补足，不读取设备端未初始化的 halo。
          if ((ip.GE.1).AND.(ip.LE.nx).AND.(jp.GE.1).AND.(jp.LE.ny)) then
            f(i,j,alpha)=f_post(ip,jp,alpha)
          endif
        enddo
      enddo
    enddo
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
! 作用: 按原文 raw 动量定义恢复 rho、u、v，不加入半步力。
!===================================================================================================

subroutine macro()
    use commondata
    implicit none
    integer(kind=4) :: i,j
    real(kind=8) :: momx,momy
    !$acc parallel loop gang vector collapse(2) default(none) present(f,rho,u,v) async(1) private(momx,momy)
    do j=1,ny
      do i=1,nx
        rho(i,j)=sum(f(i,j,0:8))
        momx=f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)
        momy=f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)
        ! Nee 式 (8)：raw 动量速度。原文式 (6) 的简单力项不另加半步力。
        u(i,j)=momx/rho(i,j)
        v(i,j)=momy/rho(i,j)
      enddo
    enddo
  end subroutine macro

!===================================================================================================
! 子程序: reconstruct_macro_from_fT
! 作用: 由检查点中的 f 与独立温度 T 恢复宏观状态。
!===================================================================================================
subroutine reconstruct_macro_from_fT()
    use commondata
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer(kind=4) :: i,j
    real(kind=8) :: momx,momy
    if((.not.all(ieee_is_finite(f))).OR.(.not.all(ieee_is_finite(T)))) &
        error stop 'Nonfinite HybridRLB4 restart distributions or temperature'
    do j=1,ny
      do i=1,nx
        rho(i,j)=sum(f(i,j,0:8))
        if(rho(i,j).LE.0.0d0) error stop 'Nonpositive HybridRLB4 restart density'
        momx=f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)
        momy=f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)
        u(i,j)=momx/rho(i,j)
        v(i,j)=momy/rho(i,j)
        Fx(i,j)=0.0d0
        Fy(i,j)=rho(i,j)*gBeta*(T(i,j)-Tref)
      enddo
    enddo
#ifdef steadyFlow
    if((.not.all(ieee_is_finite(up))).OR.(.not.all(ieee_is_finite(vp))) &
        .OR.(.not.all(ieee_is_finite(Tp))).OR.(.not.ieee_is_finite(errorU)) &
        .OR.(.not.ieee_is_finite(errorT))) error stop 'Nonfinite HybridRLB4 convergence state'
#endif
  end subroutine reconstruct_macro_from_fT


!===================================================================================================
! reconstruct_macro_from_fT end: 重新恢复 rho、u、v 和 T
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
    character(len=128) :: caseTag



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
            error6 = error6+temperatureDifference
                
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
        errorT = error5/error6                               !温度场平均L1变化 / DeltaT，不依赖温标平移
    else
        errorT = error5
    endif

  

    call append_convergence_tecplot('convergence2DOpenaccHybridRLB4.plt', restartItcOffset+itc, errorU, errorT)

    
    write(caseTag,'("HybridRLB4,Ra=",ES24.16E3,",nx=",I0,",ny=",I0)') Rayleigh,nx,ny
    call append_convergence_master_tecplot('convergence_all_2DOpenaccHybridRLB4.plt', caseTag, restartItcOffset+itc, errorU, errorT)

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
! 作用: 输出包含 f、T 的严格重启备份文件；rho/u/v/T 在读取后由 f 重建 rho/u/v；T 单独保存。
! 用途: 在运行过程中定期调用，也在程序结束前调用。
!===================================================================================================
  subroutine output_ReloadFile()                                  !输出fg，存储在当前路径，名字由 reloadFilePrefix 控制
    use commondata                                                !用于重启，包含 f,T；读入后调用 reconstruct_macro_from_fT() 重建宏观量
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
    ! HybridRLB4 restart stores its own scheme/grid header and D2Q9 f and scalar T.
    write(05) reloadMagic,nx,ny,9,0
    write(05) (((real(f(i,j,alpha),kind=8), i=1,nx), j=1,ny), alpha=0,8)
    write(05) T
#ifdef steadyFlow
    write(05) up,vp,Tp,errorU,errorT
#endif
    close(05)
    call write_reload_metadata(trim(filename))
    
    open(unit=00,file=trim(settingsFile),status='unknown',position='append')
    write(00,*) "Backup f/T restart state to: ", trim(reloadFilePrefix), "-", trim(filename),".bin"
    write(00,*) "Backup restart metadata to: ", trim(reloadFilePrefix), "-latest.meta"
    close(00)
    
    return
  end subroutine output_ReloadFile
!===================================================================================================
! output_ReloadFile 结束: 输出包含 f、T 的严格重启备份文件。
!===================================================================================================


!===================================================================================================
! 子程序: write_reload_metadata
! 作用: 覆盖写出最新 reload 续算账本，恢复累计步数、t_ff、输出编号和最新 .bin 文件名。
!===================================================================================================

subroutine write_reload_metadata(filename)
    use commondata
    implicit none
    character(len=*),intent(in) :: filename
    integer(kind=4) :: metaUnit,totalItc
    real(kind=8) :: physicalParameters(14)
    character(len=32) :: flowMode,caseName
    totalItc=restartItcOffset+itc
    call hybrid_restart_parameters(physicalParameters,flowMode,caseName)
    open(newunit=metaUnit,file=trim(reloadFilePrefix)//'-latest.meta', &
         status='replace',action='write',form='formatted')
    write(metaUnit,'(A,1X,I0)') 'reload_meta_version',801
    write(metaUnit,'(A,1X,A)') 'scheme',reloadMagic
    write(metaUnit,'(A,1X,A)') 'flowMode',trim(flowMode)
    write(metaUnit,'(A,1X,A)') 'caseName',trim(caseName)
    write(metaUnit,'(A,1X,I0)') 'nx',nx
    write(metaUnit,'(A,1X,I0)') 'ny',ny
    write(metaUnit,'(A,1X,14(ES24.16E3,1X))') 'physicalParameters',physicalParameters
    write(metaUnit,'(A,1X,A)') 'reloadFileName',trim(filename)
    write(metaUnit,'(A,1X,I0)') 'itc_total',totalItc
    write(metaUnit,'(A,1X,ES24.16E3)') 'time_tf',dble(totalItc)/timeUnit
    write(metaUnit,'(A,1X,I0)') 'snapshotFileNum',snapshotFileNum
    write(metaUnit,'(A,1X,I0)') 'pltFileNum',pltFileNum
    write(metaUnit,'(A,1X,I0)') 'reloadFileNum',reloadFileNum
    write(metaUnit,'(A,1X,I0)') 'nureSampleCount',sampleCountOffset+dimensionlessTime
    close(metaUnit)
  end subroutine write_reload_metadata

!===================================================================================================


!===================================================================================================
! 子程序: read_reload_metadata
! 作用: 严格读取 latest .meta；缺失或不匹配则停止。
!===================================================================================================

subroutine read_reload_metadata(reloadFileName)
    use commondata
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    character(len=*),intent(inout) :: reloadFileName
    integer(kind=4) :: metaUnit,ios,metaVersion,metaNx,metaNy,metaItc
    integer(kind=4) :: savedSnapshot,savedPlt,savedReload,savedSamples
    real(kind=8) :: savedParameters(14),expectedParameters(14),metaTf
    character(len=64) :: label
    character(len=32) :: savedScheme,savedFlow,savedCase,currentFlow,currentCase
    logical :: exists
    inquire(file=trim(reloadFilePrefix)//'-latest.meta',exist=exists)
    if(.not.exists) error stop 'HybridRLB4 restart requires its matching latest.meta file'
    call hybrid_restart_parameters(expectedParameters,currentFlow,currentCase)
    open(newunit=metaUnit,file=trim(reloadFilePrefix)//'-latest.meta', &
         status='old',action='read',form='formatted')
    read(metaUnit,*,iostat=ios) label,metaVersion
    if(ios.NE.0) error stop 'Invalid HybridRLB4 metadata'
    if((label.NE.'reload_meta_version').OR.(metaVersion.NE.801)) error stop 'Not an HybridRLB4 metadata file'
    read(metaUnit,*,iostat=ios) label,savedScheme
    if(ios.NE.0) error stop 'Invalid HybridRLB4 scheme entry'
    if((label.NE.'scheme').OR.(savedScheme.NE.reloadMagic)) error stop 'HybridRLB4 scheme mismatch'
    read(metaUnit,*,iostat=ios) label,savedFlow
    if(ios.NE.0) error stop 'Invalid HybridRLB4 flow mode'
    if((label.NE.'flowMode').OR.(savedFlow.NE.currentFlow)) error stop 'HybridRLB4 flow mode mismatch'
    read(metaUnit,*,iostat=ios) label,savedCase
    if(ios.NE.0) error stop 'Invalid HybridRLB4 case entry'
    if((label.NE.'caseName').OR.(savedCase.NE.currentCase)) error stop 'HybridRLB4 boundary/case mismatch'
    read(metaUnit,*,iostat=ios) label,metaNx
    if(ios.NE.0) error stop 'Invalid HybridRLB4 nx'
    if((label.NE.'nx').OR.(metaNx.NE.nx)) error stop 'HybridRLB4 nx mismatch'
    read(metaUnit,*,iostat=ios) label,metaNy
    if(ios.NE.0) error stop 'Invalid HybridRLB4 ny'
    if((label.NE.'ny').OR.(metaNy.NE.ny)) error stop 'HybridRLB4 ny mismatch'
    read(metaUnit,*,iostat=ios) label,savedParameters
    if(ios.NE.0) error stop 'Invalid HybridRLB4 physical parameters'
    if(label.NE.'physicalParameters') error stop 'Missing HybridRLB4 physical parameters'
    if(.not.all(ieee_is_finite(savedParameters))) error stop 'Nonfinite HybridRLB4 restart parameters'
    if(any(abs(savedParameters-expectedParameters).GT.1.0d-12*max(1.0d0,abs(expectedParameters)))) &
        error stop 'HybridRLB4 restart physical parameters or sample interval differ'
    read(metaUnit,*,iostat=ios) label,reloadFileName
    if(ios.NE.0) error stop 'Invalid HybridRLB4 binary name'
    if(label.NE.'reloadFileName') error stop 'Missing HybridRLB4 binary name'
    read(metaUnit,*,iostat=ios) label,metaItc
    if(ios.NE.0) error stop 'Invalid HybridRLB4 step offset'
    if((label.NE.'itc_total').OR.(metaItc.LT.0)) error stop 'Invalid HybridRLB4 step offset'
    read(metaUnit,*,iostat=ios) label,metaTf
    if(ios.NE.0) error stop 'Invalid HybridRLB4 time offset'
    if((label.NE.'time_tf').OR.(.not.ieee_is_finite(metaTf))) error stop 'Invalid HybridRLB4 time offset'
    if(abs(metaTf-dble(metaItc)/timeUnit).GT.1.0d-12*max(1.0d0,abs(metaTf))) &
        error stop 'HybridRLB4 time and step offsets disagree'
    read(metaUnit,*,iostat=ios) label,savedSnapshot
    if(ios.NE.0) error stop 'Invalid HybridRLB4 snapshot counter'
    if((label.NE.'snapshotFileNum').OR.(savedSnapshot.LT.0)) error stop 'Invalid HybridRLB4 snapshot counter'
    read(metaUnit,*,iostat=ios) label,savedPlt
    if(ios.NE.0) error stop 'Invalid HybridRLB4 Tecplot counter'
    if((label.NE.'pltFileNum').OR.(savedPlt.LT.0)) error stop 'Invalid HybridRLB4 Tecplot counter'
    read(metaUnit,*,iostat=ios) label,savedReload
    if(ios.NE.0) error stop 'Invalid HybridRLB4 reload counter'
    if((label.NE.'reloadFileNum').OR.(savedReload.LT.0)) error stop 'Invalid HybridRLB4 reload counter'
    read(metaUnit,*,iostat=ios) label,savedSamples
    if(ios.NE.0) error stop 'Invalid HybridRLB4 Nu/Re counter'
    if((label.NE.'nureSampleCount').OR.(savedSamples.LT.0)) error stop 'Invalid HybridRLB4 Nu/Re counter'
    close(metaUnit)
    restartItcOffset=metaItc
    reloadDimensionlessTime=metaTf
    snapshotFileNum=savedSnapshot
    pltFileNum=savedPlt
    reloadFileNum=savedReload
    sampleCountOffset=savedSamples
    call validate_restart_nure_history()
    reloadMetadataLoaded=.true.
  end subroutine read_reload_metadata

!===================================================================================================


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

    !----zone ------------------------------------------------------------ !再写一次 299.0：Tecplot 规范里“zone header之后的数据描述块”会以 zone marker 开始。
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
! 子程序: calNuRe
! 作用: 计算体平均 Nu / Re，并把时间序列缓存到数组中。
!===================================================================================================
  subroutine calNuRe()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: NuVolAvg_temp    !体平均 Nu
    real(kind=8) :: ReVolAvg_temp    !体平均 Re
    real(kind=8) :: sampleTime
    logical :: exNu, exRe
    logical, save :: first_nure_write = .true.
    

    !$acc wait(1)
    if (dimensionlessTime.GE.dimensionlessTimeMax) then
        write(*,*) "Error: dimensionlessTime exceeds dimensionlessTimeMax, please enlarge dimensionlessTimeMax"
        open(unit=00,file=trim(settingsFile),status="unknown",position="append")
        write(00,*) "Error: dimensionlessTime exceeds dimensionlessTimeMax, please enlarge dimensionlessTimeMax"
        close(00)
        stop
    endif

    call check_fields_hybrid()
    dimensionlessTime = dimensionlessTime+1   !每隔 outputSnapshotInterval 个自由落体时间调用一次calNuRe
#ifdef steadyFlow
    sampleTime = real(restartItcOffset+itc,kind=8)
#endif
#ifdef unsteadyFlow
    sampleTime = real(sampleCountOffset+dimensionlessTime,kind=8)*outputSnapshotInterval
#endif

    if((first_nure_write).AND.(loadInitField.EQ.1).AND.(sampleCountOffset.GT.0)) then
        inquire(file="Nu_VolAvg_2DOpenaccHybridRLB4.dat", exist=exNu)
        inquire(file="Re_VolAvg_2DOpenaccHybridRLB4.dat", exist=exRe)
        if((.not.exNu).OR.(.not.exRe)) then
            write(*,*) "Error: restart requested but old Nu/Re time-series files are missing."
            open(unit=00,file=trim(settingsFile),status="unknown",position="append")
            write(00,*) "Error: restart requested but old Nu/Re time-series files are missing."
            write(00,*) "Nu_VolAvg_2DOpenaccHybridRLB4.dat exists =", exNu
            write(00,*) "Re_VolAvg_2DOpenaccHybridRLB4.dat exists =", exRe
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


    NuVolAvg(dimensionlessTime) = NuVolAvg_temp/dble(nx*ny)*lengthUnit/(diffusivity*temperatureDifference)+1.0d0    !!体平均 Nusselt 数 = 1 + (常数系数) × 体平均对流热通量

    if((first_nure_write).AND.(sampleCountOffset.EQ.0)) then
        open(unit=01,file="Nu_VolAvg_2DOpenaccHybridRLB4.dat",status='replace',action='write')
    else
        open(unit=01,file="Nu_VolAvg_2DOpenaccHybridRLB4.dat",status='unknown',position='append',action='write')
    endif
    write(01,'(ES24.16E3,1X,ES24.16E3)') &
        real(sampleTime,kind=8), &
        real(NuVolAvg(dimensionlessTime),kind=8)   !以格子步数或者自由落体时间来写入
    close(01)

    ReVolAvg_temp = 0.0d0
    !$acc parallel loop collapse(2) default(none) present(u,v) reduction(+:ReVolAvg_temp)
    do j = 1, ny
        do i = 1, nx 
            ReVolAvg_temp = ReVolAvg_temp+(u(i,j)*u(i,j)+v(i,j)*v(i,j))
        enddo
    enddo
    ReVolAvg(dimensionlessTime) = dsqrt(ReVolAvg_temp/dble(nx*ny))*lengthUnit/viscosity    !全域体平均 RMS-Reynolds 数


    if((first_nure_write).AND.(sampleCountOffset.EQ.0)) then
        open(unit=02,file="Re_VolAvg_2DOpenaccHybridRLB4.dat",status='replace',action='write')
    else
        open(unit=02,file="Re_VolAvg_2DOpenaccHybridRLB4.dat",status='unknown',position='append',action='write')
    endif
    write(02,'(ES24.16E3,1X,ES24.16E3)') &
        real(sampleTime,kind=8), &
        real(ReVolAvg(dimensionlessTime),kind=8)
    close(02)
    first_nure_write = .false.
    write(*,'(a,1x,ES24.16E3)') "NuVolAvg =", real(NuVolAvg(dimensionlessTime),kind=8)
    write(*,'(a,1x,ES24.16E3)') "ReVolAvg =", real(ReVolAvg(dimensionlessTime),kind=8)
    return
  end subroutine calNuRe
!===================================================================================================
! calNuRe 结束: 计算体平均 Nu 和 Re 的时间历程统计量。
!===================================================================================================


#ifdef unsteadyFlow
!===================================================================================================
! Subroutine: output_unsteady_NuRe_postprocess
! Purpose: rebuild unsteady Nu/Re series, running means, and window averages from full .dat history.
!===================================================================================================
  subroutine output_unsteady_NuRe_postprocess()
    use commondata
    implicit none

    integer(kind=4) :: k, history_count
    integer(kind=4) :: whole_count, first_count, second_count
    integer(kind=4) :: iosNu, iosRe
    integer(kind=4) :: nuUnit, reUnit, seriesUnit, runningUnit
    real(kind=8) :: timeLoc, timeNu, timeRe, NuVal, ReVal, Nu_Accum, Re_Accum
    real(kind=8) :: startTf, midTf, endTf
    real(kind=8) :: Nu_WholeSum, Re_WholeSum, Nu_FirstSum, Re_FirstSum, Nu_SecondSum, Re_SecondSum
    real(kind=8) :: Nu_WholeAvg, Re_WholeAvg, Nu_FirstAvg, Re_FirstAvg, Nu_SecondAvg, Re_SecondAvg
    real(kind=8) :: Nu_FirstRelErr, Re_FirstRelErr, Nu_SecondRelErr, Re_SecondRelErr
    logical :: exNu, exRe

    inquire(file='Nu_VolAvg_2DOpenaccHybridRLB4.dat', exist=exNu)
    inquire(file='Re_VolAvg_2DOpenaccHybridRLB4.dat', exist=exRe)
    if((.not.exNu).or.(.not.exRe)) then
        write(*,'(A)') 'Error: Nu/Re history files are missing before postprocessing.'
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,'(A)') 'Error: Nu/Re history files are missing before postprocessing.'
        close(00)
        error stop 1
    endif

    open(newunit=nuUnit, file='Nu_VolAvg_2DOpenaccHybridRLB4.dat', status='old', action='read', form='formatted')
    open(newunit=reUnit, file='Re_VolAvg_2DOpenaccHybridRLB4.dat', status='old', action='read', form='formatted')

    ! These files are derived views of the full .dat history, so rebuild one continuous ZONE.
    open(newunit=seriesUnit, file='NuRe_VolAvg_2DOpenaccHybridRLB4.plt', status='replace', action='write', form='formatted')
    write(seriesUnit,'(A)') 'TITLE = "2D OpenACC Nu/Re volume averages"'
    write(seriesUnit,'(A)') 'VARIABLES = "time" "NuVolAvg" "ReVolAvg"'
    write(seriesUnit,'(A)') 'ZONE T="NuReVolAvg", F=POINT'

    open(newunit=runningUnit, file='NuRe_VolAvg_runningMean_2DOpenaccHybridRLB4.plt', status='replace', action='write', form='formatted')
    write(runningUnit,'(A)') 'TITLE = "2D OpenACC Nu/Re running means"'
    write(runningUnit,'(A)') 'VARIABLES = "time" "NuVolAvgMean" "ReVolAvgMean"'
    write(runningUnit,'(A)') 'ZONE T="NuReRunningMean", F=POINT'

    startTf = unsteadyAverageStartTf
    midTf = unsteadyAverageMidTf
    endTf = unsteadyAverageEndTf

    Nu_Accum = 0.0d0
    Re_Accum = 0.0d0
    Nu_WholeSum = 0.0d0
    Re_WholeSum = 0.0d0
    Nu_FirstSum = 0.0d0
    Re_FirstSum = 0.0d0
    Nu_SecondSum = 0.0d0
    Re_SecondSum = 0.0d0
    whole_count = 0
    first_count = 0
    second_count = 0
    history_count = 0

    do k = 1, sampleCountOffset+dimensionlessTime
        read(nuUnit,*,iostat=iosNu) timeNu, NuVal
        read(reUnit,*,iostat=iosRe) timeRe, ReVal
        ! iostat: 0=成功读到一行；小于0=到达文件末尾；大于0=格式或读入错误。
        ! 循环内只要不是 0，就说明文件短于实际记录样本数，或者某一行格式坏了。
        if((iosNu.NE.0).OR.(iosRe.NE.0)) then
            write(*,'(A)') 'Error: Nu/Re history files are shorter than the recorded sample count or contain invalid rows.'
            open(unit=00,file=trim(settingsFile),status='unknown',position='append')
            write(00,'(A)') 'Error: Nu/Re history files are shorter than the recorded sample count or contain invalid rows.'
            close(00)
            error stop 1
        endif
        if(abs(timeNu-timeRe).GT.1.0d-10*max(1.0d0,abs(timeNu))) then     !确保 Nu 和 Re 是同一个时间采样点的数据，不是错行配对的数据
            write(*,'(A)') 'Error: Nu/Re history time columns do not match.'
            open(unit=00,file=trim(settingsFile),status='unknown',position='append')
            write(00,'(A)') 'Error: Nu/Re history time columns do not match.'
            close(00)
            error stop 1
        endif

        timeLoc = timeNu
        history_count = k
        Nu_Accum = Nu_Accum + NuVal
        Re_Accum = Re_Accum + ReVal
        write(seriesUnit,'(ES24.16E3,1X,ES24.16E3,1X,ES24.16E3)') timeLoc, NuVal, ReVal
        write(runningUnit,'(ES24.16E3,1X,ES24.16E3,1X,ES24.16E3)') &
            timeLoc, Nu_Accum/dble(history_count), Re_Accum/dble(history_count)

        if ((timeLoc >= startTf) .and. (timeLoc <= endTf)) then
            Nu_WholeSum = Nu_WholeSum + NuVal
            Re_WholeSum = Re_WholeSum + ReVal
            whole_count = whole_count + 1
        endif
        if ((timeLoc >= startTf) .and. (timeLoc < midTf)) then
            Nu_FirstSum = Nu_FirstSum + NuVal
            Re_FirstSum = Re_FirstSum + ReVal
            first_count = first_count + 1
        endif
        if ((timeLoc >= midTf) .and. (timeLoc <= endTf)) then
            Nu_SecondSum = Nu_SecondSum + NuVal
            Re_SecondSum = Re_SecondSum + ReVal
            second_count = second_count + 1
        endif
    enddo

    ! 上面的循环已经读完实际记录的样本行；这里再试读一行，
    ! 不是为了继续计算，而是确认 Nu/Re 两个历史文件后面没有多余数据。
    read(nuUnit,*,iostat=iosNu) timeNu, NuVal
    read(reUnit,*,iostat=iosRe) timeRe, ReVal
    ! 任意一个 iostat 等于 0，都表示至少一个文件还成功读到了额外一行。
    ! 如果放过这种情况，后处理会静默丢掉超过记录样本数的尾部数据。
    if((iosNu.EQ.0).OR.(iosRe.EQ.0)) then
        write(*,'(A)') 'Error: Nu/Re history files contain more rows than the recorded sample count.'
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,'(A)') 'Error: Nu/Re history files contain more rows than the recorded sample count.'
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

    if (history_count <= 0) then
        write(*,'(A)') 'Error: no Nu/Re history samples were found before postprocessing.'
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,'(A)') 'Error: no Nu/Re history samples were found before postprocessing.'
        close(00)
        error stop 1
    endif

    if ((whole_count <= 0) .or. (first_count <= 0) .or. (second_count <= 0)) then
        write(*,'(A)') 'Error: no complete unsteady average window was found for Nu/Re postprocessing.'
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,'(A)') 'Error: no complete unsteady average window was found for Nu/Re postprocessing.'
        close(00)
        error stop 1
    endif

    Nu_WholeAvg = Nu_WholeSum / dble(whole_count)
    Re_WholeAvg = Re_WholeSum / dble(whole_count)
    Nu_FirstAvg = Nu_FirstSum / dble(first_count)
    Re_FirstAvg = Re_FirstSum / dble(first_count)
    Nu_SecondAvg = Nu_SecondSum / dble(second_count)
    Re_SecondAvg = Re_SecondSum / dble(second_count)
    ! tiny(1.0d0) 是双精度最小正规正数；这里用 max 避免整体平均值为 0 时相对误差除以 0。
    Nu_FirstRelErr = abs(Nu_FirstAvg - Nu_WholeAvg) / max(abs(Nu_WholeAvg), tiny(1.0d0))
    Re_FirstRelErr = abs(Re_FirstAvg - Re_WholeAvg) / max(abs(Re_WholeAvg), tiny(1.0d0))
    Nu_SecondRelErr = abs(Nu_SecondAvg - Nu_WholeAvg) / max(abs(Nu_WholeAvg), tiny(1.0d0))
    Re_SecondRelErr = abs(Re_SecondAvg - Re_WholeAvg) / max(abs(Re_WholeAvg), tiny(1.0d0))

    open(unit=33, file='NuRe_TimeAverage_2DOpenaccHybridRLB4.txt', status='replace', action='write', form='formatted')
    write(33,'(A)') '# 2D OpenACC Nu/Re statistical-convergence window averages'
    write(33,'(A)') '# start_tf mid_tf end_tf whole_count first_count second_count ' // &
        'Nu_whole Re_whole Nu_first Re_first Nu_second Re_second ' // &
        'Nu_first_relerr Re_first_relerr Nu_second_relerr Re_second_relerr'
    write(33,'(ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,I0,1X,I0,1X,I0,1X,' // &
        'ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,' // &
        'ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3)') &
        startTf, midTf, endTf, whole_count, first_count, second_count, &
        Nu_WholeAvg, Re_WholeAvg, Nu_FirstAvg, Re_FirstAvg, Nu_SecondAvg, Re_SecondAvg, &
        Nu_FirstRelErr, Re_FirstRelErr, Nu_SecondRelErr, Re_SecondRelErr
    close(33)

    write(*,'(A)') 'Unsteady Nu/Re statistical postprocessing:'
    write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'Average window:', startTf, 'to', endTf
    write(*,'(A,1X,I0,1X,A,1X,I0,1X,A,1X,I0)') &
        'Samples whole/first/second:', whole_count, '/', first_count, '/', second_count
    write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'Whole Nu/Re average:', Nu_WholeAvg, '/', Re_WholeAvg
    write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'First-half Nu/Re average:', Nu_FirstAvg, '/', Re_FirstAvg
    write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'Second-half Nu/Re average:', Nu_SecondAvg, '/', Re_SecondAvg
    write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'First-half Nu/Re relative error:', Nu_FirstRelErr, '/', Re_FirstRelErr
    write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'Second-half Nu/Re relative error:', Nu_SecondRelErr, '/', Re_SecondRelErr

  end subroutine output_unsteady_NuRe_postprocess
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
    real(kind=8) :: deltaT, coef

    ! 网格间距
    dx = 1.0d0 / lengthUnit
    deltaT = Thot - Tcold
    coef   = (lengthUnit/diffusivity)

    sum_qx = 0.0d0

    !$acc parallel loop collapse(2) default(none) present(u,T) reduction(+:sum_qx) private(dT,qx) firstprivate(dx,coef)
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
    Nu_global = (sum_qx / dble(nx*ny)) / deltaT

    ! 屏幕输出
    write(*,'(a,1x,ES24.16E3)') "Nu_global =", real(Nu_global,kind=8)


    ! 同步写入日志
    open(unit=00,file=trim(settingsFile),status="unknown",position="append")
    if((errorU.LE.epsU).AND.(errorT.LE.epsT)) then
        write(00,*) 'Steady criteria satisfied.'
    else
        write(00,*) 'Step limit reached before steady criteria; final diagnostics are not converged benchmarks.'
    endif
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
    real(kind=8) :: dx, dy, deltaT,coef
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
    deltaT = Thot - Tcold
    coef   = (lengthUnit/diffusivity)

   

    !-----------------------------
    ! (1) 左侧热壁平均 Nu_hot，同时记录 Numax/Numin 及其 y 位置
    sum_hot = 0.0d0
    do j = 1, ny
      ! 壁面导热通量：qx(x=0,j)
      qx_wall = (8.0d0*Thot - 9.0d0*T(1,j) + T(2,j)) / (3.0d0*dx)
      Nu_left(j)= qx_wall / deltaT
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

    Nu_left_ext(0) = ( 2.0d0 * (Thot- T_wb) / dx ) / deltaT   ! 角点局部 Nu

    ! ---------- 左上角：在 i=1 这条竖线上，用 j=ny-3..ny 拟合得到 y=1 的温度 ----------
    yfit(1) = yp(ny-3);  Tfit(1) = T(1,ny-3)
    yfit(2) = yp(ny-2);  Tfit(2) = T(1,ny-2)
    yfit(3) = yp(ny-1);  Tfit(3) = T(1,ny-1)
    yfit(4) = yp(ny  );  Tfit(4) = T(1,ny  )

    call fit_adiabatic_wall_T4(yp(ny+1), yfit, Tfit, T_wt)   ! 得到顶壁温度 T(y=yp(ny+1))

    Nu_left_ext(ny+1) = ( 2.0d0 * (Thot-T_wt) / dx ) / deltaT  ! 角点局部 Nu





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
      sum_cold = sum_cold + qx_wall/ deltaT
    enddo
    Nu_cold = (sum_cold / dble(ny))

    !-----------------------------
    ! (3) 竖直中线 x=1/2 的平均 Nu_middle
    sum_mid = 0.0d0

    if (mod(nx,2) == 1) then
      iMid = (nx + 1)/2

      do j = 1, ny
        sum_mid = sum_mid + ( coef*u(iMid,j)*(T(iMid,j)-Tref) + (T(iMid-1,j)-T(iMid+1,j))/(2.0d0*dx) ) / deltaT
      enddo

    else
      iL = nx/2
      iR = iL + 1

      do j = 1, ny
        sum_mid = sum_mid + (coef*( 0.5d0*( u(iL,j)*(T(iL,j)-Tref) + u(iR,j)*(T(iR,j)-Tref) )) &
        + (T(iL,j)-T(iR,j))/dx )/ deltaT
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
  real(kind=8) :: deltaT, coef

  dy     = 1.0d0 / lengthUnit
  deltaT = Thot - Tcold
  coef   = (lengthUnit/diffusivity)

  sum_qy = 0.0d0

  !$acc parallel loop collapse(2) default(none) present(v,T) reduction(+:sum_qy) private(dTdy,qy) firstprivate(dy,coef)
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
  Nu_global = (sum_qy / dble(nx*ny)) / deltaT

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
  real(kind=8) :: dx, dy, deltaT
  real(kind=8) :: qy_wall, sum_hot, sum_cold, sum_mid, coef
  real(kind=8), dimension(1:nx) :: Nu_bot
  real(kind=8), dimension(0:nx+1) :: Nu_bot_ext
  real(kind=8) :: xfit(4), Tfit(4), T_wl, T_wr
  real(kind=8) :: xk(5), fk(5)
  real(kind=8) :: fstar, xstar


  dx     = 1.0d0 / lengthUnit
  dy     = 1.0d0 / lengthUnit
  deltaT = Thot - Tcold
  coef   = (lengthUnit/diffusivity)

  !-----------------------------
  ! (1) 底部热壁平均 Nu_hot（不含角点）
  sum_hot = 0.0d0
  do i = 1, nx
    qy_wall   = (8.0d0*Thot - 9.0d0*T(i,1) + T(i,2)) / (3.0d0*dy)
    Nu_bot(i)= qy_wall / deltaT
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
  Nu_bot_ext(0)    = ( 2.0d0 * (Thot - T_wl) / dy ) / deltaT
  Nu_bot_ext(nx+1) = ( 2.0d0 * (Thot - T_wr) / dy ) / deltaT

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
    sum_cold = sum_cold + qy_wall / deltaT
  enddo
  Nu_cold = sum_cold / dble(nx)

  !-----------------------------
  ! 中线的 Nusselt 数
  sum_mid = 0.0d0

  if (mod(ny,2) == 1) then
    jMid = (ny + 1)/2

    do i = 1, nx
      sum_mid = sum_mid + ( coef*v(i,jMid)*(T(i,jMid)-Tref) - (T(i,jMid+1)-T(i,jMid-1))/(2.0d0*dy) ) / deltaT
    enddo

  else
    jB = ny/2
    jT = jB + 1

    do i = 1, nx
      sum_mid = sum_mid + (coef*( 0.5d0*( v(i,jB)*(T(i,jB)-Tref) + v(i,jT)*(T(i,jT)-Tref) )) &
      + (T(i,jB)-T(i,jT))/dy ) / deltaT
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
! 用途: 在主程序结束阶段调用，作为统一后处理的一部分。
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




!===================================================================================================
! 子程序: hybrid_equilibrium
! 作用: 用完整中心矩平衡态初始化 D2Q9，避免把高阶平衡态退化成二阶截断多项式。
!===================================================================================================



!===================================================================================================

!===================================================================================================
! 子程序: hybrid_restart_parameters
! 作用: 保存影响 f/T 物理含义和历史时间轴的参数；目标总时长可增加，其余必须匹配。
!===================================================================================================
  subroutine hybrid_restart_parameters(parameters,flowMode,caseName)
    use commondata
    implicit none
    real(kind=8),intent(out) :: parameters(14)
    character(len=*),intent(out) :: flowMode,caseName
    parameters=(/Rayleigh,Prandtl,Mach,Thot,Tcold,cs2,rlbRelaxationRate,diffusivity, &
        timeUnit,outputSnapshotInterval,dble(thermalSubsteps),dble(steadyCheckIntervalItc),2.0d0,1.0d0/)
#ifdef steadyFlow
    flowMode='steadyFlow'
#else
    flowMode='unsteadyFlow'
#endif
#ifdef SideHeatedCell
    caseName='SideHeatedCell'
#else
    caseName='RayleighBenardCell'
#endif
  end subroutine hybrid_restart_parameters
!===================================================================================================
!===================================================================================================
! 子程序: validate_restart_nure_history
! 作用: 续算历史必须与检查点同一边界；禁止悄悄拼接丢失或超前的样本。
!===================================================================================================
  subroutine validate_restart_nure_history()
    use commondata
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer(kind=4) :: nuUnit,reUnit,k,iosNu,iosRe
    real(kind=8) :: tn,tr,nuv,rev,expectedTime,previousTime
    logical :: nuExists,reExists
    inquire(file='Nu_VolAvg_2DOpenaccHybridRLB4.dat',exist=nuExists)
    inquire(file='Re_VolAvg_2DOpenaccHybridRLB4.dat',exist=reExists)
    if((sampleCountOffset.EQ.0).AND.(.not.nuExists).AND.(.not.reExists)) return
    if((.not.nuExists).OR.(.not.reExists)) error stop 'HybridRLB4 restart Nu/Re histories missing'
    open(newunit=nuUnit,file='Nu_VolAvg_2DOpenaccHybridRLB4.dat',status='old',action='read')
    open(newunit=reUnit,file='Re_VolAvg_2DOpenaccHybridRLB4.dat',status='old',action='read')
    previousTime=-1.0d0
    do k=1,sampleCountOffset
      read(nuUnit,*,iostat=iosNu) tn,nuv
      read(reUnit,*,iostat=iosRe) tr,rev
      if((iosNu.NE.0).OR.(iosRe.NE.0)) error stop 'HybridRLB4 restart histories are shorter than checkpoint'
      if(.not.all(ieee_is_finite((/tn,tr,nuv,rev/)))) error stop 'Nonfinite HybridRLB4 restart histories'
#ifdef unsteadyFlow
      expectedTime=dble(k)*outputSnapshotInterval
      if(max(abs(tn-expectedTime),abs(tr-expectedTime)).GT.1.0d-10*max(1.0d0,expectedTime)) &
          error stop 'HybridRLB4 restart history time mismatch'
#else
      if(abs(tn-tr).GT.1.0d-10*max(1.0d0,abs(tn))) error stop 'HybridRLB4 steady history time mismatch'
      if((tn.LE.previousTime).OR.(tn.GT.dble(restartItcOffset))) &
          error stop 'HybridRLB4 steady history is not ordered within checkpoint'
#endif
      previousTime=tn
    enddo
    read(nuUnit,*,iostat=iosNu) tn,nuv
    read(reUnit,*,iostat=iosRe) tr,rev
    if((iosNu.GE.0).OR.(iosRe.GE.0)) error stop 'HybridRLB4 histories extend beyond checkpoint; align them first'
    close(nuUnit)
    close(reUnit)
  end subroutine validate_restart_nure_history
!===================================================================================================
! 子程序: check_fields_hybrid
! 作用: 在检查/采样边界检查有限性，避免 NaN 被误当成稳态收敛后继续输出。
!===================================================================================================
  subroutine check_fields_hybrid()
    use commondata
    implicit none
    integer(kind=4) :: i,j,badField
    badField=0
    !$acc wait(1)
    !$acc parallel loop collapse(2) default(none) present(rho,u,v,T) reduction(max:badField)
    do j=1,ny
      do i=1,nx
        if(.not.((rho(i,j).GT.0.0d0).AND.(rho(i,j).LE.huge(1.0d0)) &
            .AND.(abs(u(i,j)).LE.huge(1.0d0)).AND.(abs(v(i,j)).LE.huge(1.0d0)) &
            .AND.(abs(T(i,j)).LE.huge(1.0d0)))) badField=1
      enddo
    enddo
    if(badField.NE.0) error stop 'Nonfinite HybridRLB4 state or nonpositive density'
  end subroutine check_fields_hybrid
!===================================================================================================


!===================================================================================================
! rlb_equilibrium：式 (14) 的两个一维二阶 Hermite 因子相乘，保留全部三阶、四阶混合项。
! D2Q9 方向按本项目固定为 0,E,N,W,S,NE,NW,SW,SE；cs2 固定为标准 D2Q9 的 1/3。
! 局部常量避免设备例程隐式捕获可变模块方向、权重数组。
!===================================================================================================
subroutine rlb_equilibrium(rhoLocal,ux,uy,popEq)
    use commondata, only: cs2
    implicit none
    !$acc routine seq
    real(kind=8), intent(in) :: rhoLocal, ux, uy
    real(kind=8), intent(out) :: popEq(0:8)
    integer(kind=4), parameter :: latticeX(0:8)=(/0,1,0,-1,0,1,-1,-1,1/)
    integer(kind=4), parameter :: latticeY(0:8)=(/0,0,1,0,-1,1,1,-1,-1/)
    real(kind=8), parameter :: latticeWeight(0:8)=(/4.0d0/9.0d0, &
        1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0, &
        1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0/)
    integer(kind=4) :: alpha
    real(kind=8) :: cx, cy, hxx, hyy, factorX, factorY

    do alpha=0,8
        cx=dble(latticeX(alpha))
        cy=dble(latticeY(alpha))
        hxx=cx*cx-cs2
        hyy=cy*cy-cs2
        factorX=1.0d0+ux*cx/cs2+ux*ux*hxx/(2.0d0*cs2*cs2)
        factorY=1.0d0+uy*cy/cs2+uy*uy*hyy/(2.0d0*cs2*cs2)
        popEq(alpha)=rhoLocal*latticeWeight(alpha)*factorX*factorY
    enddo
end subroutine rlb_equilibrium

!===================================================================================================
! rlb_collide_node：式 (16) 只提取二阶非平衡张量，式 (15) 递归构造更高阶系数。
! 原文式 (15) 的 tau_ab 与式 (16) 的 xi_ab 视为同一二阶非平衡张量，未引入额外比例。
! 不读取原始分布中的三阶、四阶非平衡矩，不增加滤波或黏性修正。
!===================================================================================================
subroutine rlb_collide_node(pop,rhoLocal,ux,uy,forceX,forceY,popPost)
    use commondata, only: cs2, rlbRelaxationRate
    implicit none
    !$acc routine seq
    real(kind=8), intent(in) :: pop(0:8), rhoLocal, ux, uy, forceX, forceY
    real(kind=8), intent(out) :: popPost(0:8)
    integer(kind=4), parameter :: latticeX(0:8)=(/0,1,0,-1,0,1,-1,-1,1/)
    integer(kind=4), parameter :: latticeY(0:8)=(/0,0,1,0,-1,1,1,-1,-1/)
    real(kind=8), parameter :: latticeWeight(0:8)=(/4.0d0/9.0d0, &
        1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0, &
        1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0/)
    real(kind=8), parameter :: cs4=cs2*cs2, cs6=cs2*cs2*cs2, cs8=cs2*cs2*cs2*cs2
    integer(kind=4) :: alpha
    real(kind=8) :: popEq(0:8), difference, xiXX, xiXY, xiYY
    real(kind=8) :: a21, a12, a22, cx, cy, hxx, hyy, regularized, source

    call rlb_equilibrium(rhoLocal,ux,uy,popEq)
    xiXX=0.0d0
    xiXY=0.0d0
    xiYY=0.0d0
    do alpha=0,8
        cx=dble(latticeX(alpha))
        cy=dble(latticeY(alpha))
        difference=pop(alpha)-popEq(alpha)
        xiXX=xiXX+cx*cx*difference
        xiXY=xiXY+cx*cy*difference
        xiYY=xiYY+cy*cy*difference
    enddo

    ! 三阶与四阶 Hermite 系数：两个 1/2 和一个 1/4 来自对称指标组合，不能漏掉。
    a21=uy*xiXX+2.0d0*ux*xiXY
    a12=ux*xiYY+2.0d0*uy*xiXY
    a22=uy*uy*xiXX+4.0d0*ux*uy*xiXY+ux*ux*xiYY
    do alpha=0,8
        cx=dble(latticeX(alpha))
        cy=dble(latticeY(alpha))
        hxx=cx*cx-cs2
        hyy=cy*cy-cs2
        regularized=latticeWeight(alpha)* &
            ((xiXX*hxx+2.0d0*xiXY*cx*cy+xiYY*hyy)/(2.0d0*cs4) &
             +a21*hxx*cy/(2.0d0*cs6)+a12*cx*hyy/(2.0d0*cs6) &
             +a22*hxx*hyy/(4.0d0*cs8))
        source=latticeWeight(alpha)*(cx*forceX+cy*forceY)/cs2
        popPost(alpha)=popEq(alpha)+(1.0d0-rlbRelaxationRate)*regularized+source
    enddo
end subroutine rlb_collide_node

!===================================================================================================
! 子程序: thermal_boundary
! 作用: 每个 RK 阶段补齐半格壁面的温度 ghost；只填 halo，不覆盖流体格心。
! 恒温用二次延拓，使 ghost 数值通量与项目三点半格壁面 Nu 一致。
!===================================================================================================
  subroutine thermal_boundary(stage)
    use commondata
    implicit none
    real(kind=8),intent(inout) :: stage(0:nx+1,0:ny+1)
    integer(kind=4) :: i,j
    !$acc parallel loop default(none) present(stage) async(1)
    do j=1,ny
#ifdef SideHeatedCell
      stage(0,j)=(8.0d0*Thot-6.0d0*stage(1,j)+stage(2,j))/3.0d0
      stage(nx+1,j)=(8.0d0*Tcold-6.0d0*stage(nx,j)+stage(nx-1,j))/3.0d0
#else
      stage(0,j)=stage(1,j)
      stage(nx+1,j)=stage(nx,j)
#endif
    enddo
    !$acc parallel loop default(none) present(stage) async(1)
    do i=0,nx+1
#ifdef SideHeatedCell
      stage(i,0)=stage(i,1)
      stage(i,ny+1)=stage(i,ny)
#else
      stage(i,0)=(8.0d0*Thot-6.0d0*stage(i,1)+stage(i,2))/3.0d0
      stage(i,ny+1)=(8.0d0*Tcold-6.0d0*stage(i,ny)+stage(i,ny-1))/3.0d0
#endif
    enddo
  end subroutine thermal_boundary
!===================================================================================================

!===================================================================================================
! 子程序: thermal_rhs
! 作用: 在格子单位下离散原文式 (2)：-u*dT/dx-v*dT/dy+kappa*laplacian(T)。
! 原文未完整给出 RLB-4 空间模板；此处补齐为中心二阶，dx=dy=1。
! 常温平移不改变梯度；不增加上风耗散、滤波或人工裁剪。
!===================================================================================================
  subroutine thermal_rhs(stage,rhs)
    use commondata
    implicit none
    real(kind=8),intent(inout) :: stage(0:nx+1,0:ny+1)
    real(kind=8),intent(out) :: rhs(nx,ny)
    integer(kind=4) :: i,j
    real(kind=8) :: dtdx,dtdy,laplacian
    call thermal_boundary(stage)
    !$acc parallel loop gang vector collapse(2) default(none) present(stage,rhs,u,v) async(1) &
    !$acc& private(dtdx,dtdy,laplacian)
    do j=1,ny
      do i=1,nx
        dtdx=0.5d0*(stage(i+1,j)-stage(i-1,j))
        dtdy=0.5d0*(stage(i,j+1)-stage(i,j-1))
        laplacian=stage(i+1,j)+stage(i-1,j)+stage(i,j+1)+stage(i,j-1)-4.0d0*stage(i,j)
        rhs(i,j)=-u(i,j)*dtdx-v(i,j)*dtdy+diffusivity*laplacian
      enddo
    enddo
  end subroutine thermal_rhs
!===================================================================================================

!===================================================================================================
! 子程序: advance_temperature_rk4
! 作用: 固定本步流场速度，对温度半离散方程作经典四阶 RK 推进。
! dt 为本次温度推进的格子时间；主程序每个流场步调用 dt=1。
! 温度时间积分为 RK4；整体耦合仍受流场/显式力和顺序耦合精度限制。
!===================================================================================================
  subroutine advance_temperature_rk4(dt)
    use commondata
    implicit none
    real(kind=8),intent(in) :: dt
    real(kind=8) :: subDt
    integer(kind=4) :: i,j,substep
    subDt=dt/dble(thermalSubsteps)
    do substep=1,thermalSubsteps
      !$acc parallel loop gang vector collapse(2) default(none) &
      !$acc& present(T,temperatureBase,temperatureStage) async(1)
      do j=1,ny
        do i=1,nx
          temperatureBase(i,j)=T(i,j)
          temperatureStage(i,j)=T(i,j)
        enddo
      enddo
      call thermal_rhs(temperatureStage,temperatureRhs)
      !$acc parallel loop gang vector collapse(2) default(none) &
      !$acc& present(temperatureBase,temperatureStage,temperatureRhs,temperatureIncrement) &
      !$acc& firstprivate(subDt) async(1)
      do j=1,ny
        do i=1,nx
          temperatureIncrement(i,j)=temperatureRhs(i,j)
          temperatureStage(i,j)=temperatureBase(i,j)+0.5d0*subDt*temperatureRhs(i,j)
        enddo
      enddo
      call thermal_rhs(temperatureStage,temperatureRhs)
      !$acc parallel loop gang vector collapse(2) default(none) &
      !$acc& present(temperatureBase,temperatureStage,temperatureRhs,temperatureIncrement) &
      !$acc& firstprivate(subDt) async(1)
      do j=1,ny
        do i=1,nx
          temperatureIncrement(i,j)=temperatureIncrement(i,j)+2.0d0*temperatureRhs(i,j)
          temperatureStage(i,j)=temperatureBase(i,j)+0.5d0*subDt*temperatureRhs(i,j)
        enddo
      enddo
      call thermal_rhs(temperatureStage,temperatureRhs)
      !$acc parallel loop gang vector collapse(2) default(none) &
      !$acc& present(temperatureBase,temperatureStage,temperatureRhs,temperatureIncrement) &
      !$acc& firstprivate(subDt) async(1)
      do j=1,ny
        do i=1,nx
          temperatureIncrement(i,j)=temperatureIncrement(i,j)+2.0d0*temperatureRhs(i,j)
          temperatureStage(i,j)=temperatureBase(i,j)+subDt*temperatureRhs(i,j)
        enddo
      enddo
      call thermal_rhs(temperatureStage,temperatureRhs)
      !$acc parallel loop gang vector collapse(2) default(none) &
      !$acc& present(T,temperatureBase,temperatureRhs,temperatureIncrement) firstprivate(subDt) async(1)
      do j=1,ny
        do i=1,nx
          T(i,j)=temperatureBase(i,j)+subDt*(temperatureIncrement(i,j)+temperatureRhs(i,j))/6.0d0
        enddo
      enddo
    enddo
  end subroutine advance_temperature_rk4
!===================================================================================================

!===================================================================================================
! 子程序: output_thermal_budget
! 作用: 结束时独立报告温度半离散方程与壁面数值热流的闭合，以及离散速度散度。
! Rheat=sum(rhs)+Qout；对流形式在非严格无散速度下不保证全离散热量守恒。
! 本诊断只使用当前瞬时场，不能据此判定统计收敛或网格无关。
!===================================================================================================
  subroutine output_thermal_budget()
    use commondata
    implicit none
    integer(kind=4) :: i,j,budgetUnit
    real(kind=8) :: totalTemperature,totalRhs,outwardHeat,divSquare,divMaximum
    real(kind=8) :: uxMinus,uxPlus,vyMinus,vyPlus,divergence,heatScale,divScale
    !$acc parallel loop gang vector collapse(2) default(none) present(T,temperatureStage) async(1)
    do j=1,ny
      do i=1,nx
        temperatureStage(i,j)=T(i,j)
      enddo
    enddo
    call thermal_rhs(temperatureStage,temperatureRhs)
    !$acc wait(1)
    totalTemperature=0.0d0
    totalRhs=0.0d0
    outwardHeat=0.0d0
    divSquare=0.0d0
    divMaximum=0.0d0
    !$acc parallel loop gang vector collapse(2) default(none) &
    !$acc& present(T,u,v,temperatureStage,temperatureRhs) &
    !$acc& reduction(+:totalTemperature,totalRhs,outwardHeat,divSquare) reduction(max:divMaximum) &
    !$acc& private(uxMinus,uxPlus,vyMinus,vyPlus,divergence)
    do j=1,ny
      do i=1,nx
        totalTemperature=totalTemperature+T(i,j)
        totalRhs=totalRhs+temperatureRhs(i,j)
        if(i.EQ.1) outwardHeat=outwardHeat+diffusivity*(T(i,j)-temperatureStage(0,j))
        if(i.EQ.nx) outwardHeat=outwardHeat+diffusivity*(T(i,j)-temperatureStage(nx+1,j))
        if(j.EQ.1) outwardHeat=outwardHeat+diffusivity*(T(i,j)-temperatureStage(i,0))
        if(j.EQ.ny) outwardHeat=outwardHeat+diffusivity*(T(i,j)-temperatureStage(i,ny+1))
        ! 无滑移壁面放在半格位置，速度 ghost 取奇延拓。
        uxMinus=-u(i,j)
        uxPlus=-u(i,j)
        vyMinus=-v(i,j)
        vyPlus=-v(i,j)
        if(i.GT.1) uxMinus=u(i-1,j)
        if(i.LT.nx) uxPlus=u(i+1,j)
        if(j.GT.1) vyMinus=v(i,j-1)
        if(j.LT.ny) vyPlus=v(i,j+1)
        divergence=0.5d0*(uxPlus-uxMinus+vyPlus-vyMinus)
        divSquare=divSquare+divergence*divergence
        divMaximum=max(divMaximum,abs(divergence))
      enddo
    enddo
#ifdef SideHeatedCell
    heatScale=diffusivity*temperatureDifference*dble(ny)/lengthUnit
#else
    heatScale=diffusivity*temperatureDifference*dble(nx)/lengthUnit
#endif
    divScale=max(velocityUnit/lengthUnit,tiny(1.0d0))
    open(newunit=budgetUnit,file='ThermalBudget_2DOpenaccHybridRLB4.txt',status='replace',action='write')
    write(budgetUnit,'(A)') '# Instantaneous scalar-equation budget; lattice dx=dy=dt_flow=1'
    write(budgetUnit,'(A)') '# tf meanT sum_rhs Qout Rheat Rheat_normalized div_rms div_max div_rms_normalized'
    write(budgetUnit,'(9(ES24.16E3,1X))') dble(restartItcOffset+itc)/timeUnit, &
        totalTemperature/dble(nx*ny),totalRhs,outwardHeat,totalRhs+outwardHeat, &
        (totalRhs+outwardHeat)/max(heatScale,tiny(1.0d0)), &
        sqrt(divSquare/dble(nx*ny)),divMaximum,sqrt(divSquare/dble(nx*ny))/divScale
    close(budgetUnit)
  end subroutine output_thermal_budget
!===================================================================================================
