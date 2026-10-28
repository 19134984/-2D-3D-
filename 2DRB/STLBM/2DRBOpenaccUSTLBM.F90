!=============================================================
!!!    二维单 GPU OpenACC：D2Q9 流场 + D2Q9 内能场 USTLBM
!!!    Chen, Shu and Tan, Applied Sciences 7 (2017), 277
!!!    DOI: 10.3390/app7030277；主要实现依据：式 (39)--(47)
!!!    主体/输出风格沿用本项目二维 OpenACC 框架；旧 STLBM 文件不修改。
!!!    rho/u/v/T 是整步状态，非平衡量只是半时间、半格位置的临时缓存。
!!!    二阶校正必须用 N(i+ex,j+ey)-N(i,j)，不能沿用初版 STLBM 上游单侧和。
!!!    浮力用新旧温度的梯形平均，温度先于本步速度完成更新。
!!!    保留项目半格壁面：1:nx,1:ny 全是流体格心，壁面位于 0,nx 和 0,ny。
!!!    论文壁上节点边界改为与此几何一致的 ghost 延拓，不覆盖首末流体层。
!!!    二维内能平衡态二阶矩为 2/3，故 tauT=0.5+1.5*kappa。
!!!    原文式 (16) 的 5/3 与二维平衡态矩不一致，不能直接代入本二维实现。
!!!    时间/空间二阶来自原文推导；RB 是项目扩展，需单独做瞬态与统计验证。
!!!    早期“无条件稳定”的证明后来有修正；本程序不保证任意参数均稳定。
!=============================================================

!=============================================================
!   自定义宏：默认先验证侧壁稳态；切换 RB 时同时切换模式和温度边界。
#define steadyFlow
!#define unsteadyFlow

#if defined(steadyFlow) && defined(unsteadyFlow)
#error "Choose only one flow mode: steadyFlow or unsteadyFlow"
#endif
#if !defined(steadyFlow) && !defined(unsteadyFlow)
#error "Define steadyFlow or unsteadyFlow"
#endif

#define HorizontalWallsNoslip
#define VerticalWallsNoslip
!#define VerticalWallsPeriodicalU

!   RB：上下恒温；侧面绝热无滑移，或速度/温度同时周期。
!#define RayleighBenardCell
!#define HorizontalWallsConstT
!#define VerticalWallsAdiabatic
!#define VerticalWallsPeriodicalT

!   侧壁：左热右冷、上下绝热、四壁无滑移。
#define SideHeatedCell
#define HorizontalWallsAdiabatic
#define VerticalWallsConstT
!#define SideHeatedHa

#if defined(RayleighBenardCell) && defined(SideHeatedCell)
#error "Choose only one convection case"
#endif
#if !defined(RayleighBenardCell) && !defined(SideHeatedCell)
#error "Define RayleighBenardCell or SideHeatedCell"
#endif
#if defined(VerticalWallsNoslip) && defined(VerticalWallsPeriodicalU)
#error "Choose only one vertical velocity BC"
#endif
#if !defined(VerticalWallsNoslip) && !defined(VerticalWallsPeriodicalU)
#error "Define a vertical velocity BC"
#endif
#ifndef HorizontalWallsNoslip
#error "HorizontalWallsNoslip is required"
#endif
#ifdef SideHeatedCell
#if !defined(VerticalWallsNoslip) || !defined(HorizontalWallsAdiabatic) || !defined(VerticalWallsConstT)
#error "SideHeatedCell requires no-slip walls, horizontal adiabatic and vertical constant T"
#endif
#if defined(HorizontalWallsConstT) || defined(VerticalWallsAdiabatic) || defined(VerticalWallsPeriodicalT)
#error "Conflicting SideHeatedCell thermal BCs"
#endif
#endif
#ifdef RayleighBenardCell
#ifndef HorizontalWallsConstT
#error "RayleighBenardCell requires HorizontalWallsConstT"
#endif
#if defined(VerticalWallsAdiabatic) && defined(VerticalWallsPeriodicalT)
#error "Choose adiabatic or periodic side temperature BC"
#endif
#if !defined(VerticalWallsAdiabatic) && !defined(VerticalWallsPeriodicalT)
#error "Define adiabatic or periodic side temperature BC"
#endif
#if defined(VerticalWallsPeriodicalU) != defined(VerticalWallsPeriodicalT)
#error "Periodic RB requires both velocity and temperature periodicity"
#endif
#if defined(HorizontalWallsAdiabatic) || defined(VerticalWallsConstT) || defined(SideHeatedHa)
#error "Conflicting RayleighBenardCell options"
#endif
#endif
!   自定义宏结束
!=============================================================

!=============================================================
!   全局模块：主要参数仍在 commondata 中直接修改。
    module commondata
        implicit none
        integer(kind=4), parameter :: loadInitField=0   ! 0: 新算；1: 读取本版本 latest.meta 严格续算
        integer(kind=4) :: reloadFileNum=0             ! latest.meta 恢复的重启编号
        real(kind=8) :: reloadDimensionlessTime=0.0d0   ! 重启边界的实际累计 t_ff
        integer(kind=4) :: restartItcOffset=0
        logical :: reloadMetadataLoaded=.false.

        !===============================================================================================
        ! 默认基准：Pr=0.71、Ra=1e6 的二维侧壁自然对流；不代表已经网格无关。
        integer(kind=4), parameter :: nx=128, ny=128
#ifdef SideHeatedCell
        real(kind=8), parameter :: lengthUnit=dble(nx)
#else
        real(kind=8), parameter :: lengthUnit=dble(ny)
#endif
        real(kind=8), parameter :: pi=acos(-1.0d0)
        real(kind=8), parameter :: Rayleigh=1.0d6
        real(kind=8), parameter :: Prandtl=0.71d0
        real(kind=8), parameter :: Mach=0.1d0
        real(kind=8), parameter :: Thot=0.5d0, Tcold=-0.5d0
        real(kind=8), parameter :: Tref=0.5d0*(Thot+Tcold)
        real(kind=8), parameter :: temperatureDifference=Thot-Tcold
        real(kind=8), parameter :: tauf=0.5d0+Mach*lengthUnit*dsqrt(3.0d0*Prandtl/Rayleigh)
        real(kind=8), parameter :: viscosity=(tauf-0.5d0)/3.0d0
        real(kind=8), parameter :: diffusivity=viscosity/Prandtl
        ! D2Q9 内能模型的零速二阶矩为 (2/3)*rho*T，而不是 D2Q5 的 1/3。
        real(kind=8), parameter :: tauT=0.5d0+1.5d0*diffusivity
        real(kind=8), parameter :: gBeta=Rayleigh*viscosity*diffusivity &
            /(temperatureDifference*lengthUnit**3)
        real(kind=8), parameter :: velocityUnit=Mach/dsqrt(3.0d0)  ! 自由落体速度，等于 sqrt(gBeta*DeltaT*L)
        real(kind=8), parameter :: timeUnit=lengthUnit/velocityUnit ! 自由落体时间对应的格子步数
        real(kind=8), parameter :: velocityScaleCompare=lengthUnit/diffusivity ! 输出速度采用 UL/kappa 标度
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1
        ! RB 温度初扰动始终由此参数控制；0 表示关闭，避免高 Ra 无扰动导热态长期不启动。
        real(kind=8), parameter :: rbInitialPerturbation=1.0d-3
#ifdef SideHeatedHa
        real(kind=8), parameter :: Ha=20.0d0
        real(kind=8), parameter :: phi=0.0d0*pi/180.0d0
        real(kind=8), parameter :: B2sigemarho=Ha**2*viscosity/lengthUnit**2
#endif
        real(kind=8), parameter :: epsU=1.0d-7, epsT=1.0d-7
        integer(kind=4), parameter :: steadyCheckIntervalItc=2000 ! 收敛检查间隔
        integer(kind=4), parameter :: stateCheckIntervalItc=2000  ! 有限值和正密度检查间隔
#ifdef steadyFlow
        real(kind=8), parameter :: outputSnapshotInterval=10.0d0
        real(kind=8), parameter :: reloadFileInterval=100.0d0
        real(kind=8), parameter :: outputPltFileInterval=100.0d0
        integer(kind=4), parameter :: dimensionlessTimeMax=1200
        integer(kind=4), parameter :: itc_max=20000000
#endif
#ifdef unsteadyFlow
        real(kind=8), parameter :: outputSnapshotInterval=0.5d0
        real(kind=8), parameter :: reloadFileInterval=100.0d0
        real(kind=8), parameter :: outputPltFileInterval=100.0d0
        real(kind=8), parameter :: unsteadyRunDuration=1000.0d0
        real(kind=8), parameter :: unsteadyAverageStartTf=0.5d0*unsteadyRunDuration
        real(kind=8), parameter :: unsteadyAverageEndTf=unsteadyRunDuration
        real(kind=8), parameter :: unsteadyAverageMidTf=0.5d0*(unsteadyAverageStartTf+unsteadyAverageEndTf)
        integer(kind=4), parameter :: unsteadySampleCount=int(unsteadyRunDuration/outputSnapshotInterval+1.0d-10)
        integer(kind=4), parameter :: dimensionlessTimeMax=unsteadySampleCount
        integer(kind=4), parameter :: itc_max=max(1,int(unsteadyRunDuration*timeUnit+0.5d0))
#endif
        integer(kind=4), parameter :: outputSnapshotFile=1
        integer(kind=4), parameter :: outputPltFile=1
        integer(kind=4), parameter :: outputReloadFile=1
        integer(kind=4) :: snapshotFileNum, pltFileNum, dimensionlessTime
        integer(kind=4) :: sampleIndexOffset=0 ! 续算前已完成的绝对采样序号，不用实际重启时刻平移采样网格
        integer(kind=4) :: outputSnapshotIntervalItc, reloadFileIntervalItc, outputPltFileIntervalItc
        real(kind=8) :: NuVolAvg(0:dimensionlessTimeMax), ReVolAvg(0:dimensionlessTimeMax)
        character(len=100) :: snapshotFilePrefix='buoyancyCavity2DOpenaccUSTLBMSnapshot'
        character(len=100) :: pltFolderPrefix='buoyancyCavity2DOpenaccUSTLBMTecplot'
        character(len=100) :: reloadFilePrefix='reloadFile2DOpenaccUSTLBM'
        character(len=100) :: settingsFile='SimulationSettings2DOpenaccUSTLBM.txt'
        character(len=16), parameter :: restartSignature='USTLBM2D2017-V1'
#ifdef SideHeatedCell
#ifdef SideHeatedHa
        character(len=32), parameter :: caseConfiguration='SideHeated-noslip-Ha'
#else
        character(len=32), parameter :: caseConfiguration='SideHeated-noslip'
#endif
#else
#ifdef VerticalWallsPeriodicalU
        character(len=32), parameter :: caseConfiguration='RB-periodic'
#else
        character(len=32), parameter :: caseConfiguration='RB-noslip-adiabatic'
#endif
#endif
#ifdef SideHeatedHa
        real(kind=8), parameter :: restartPhysics(7)=[Rayleigh,Prandtl,Mach,Thot,Tcold,Ha,phi]
#else
        real(kind=8), parameter :: restartPhysics(7)=[Rayleigh,Prandtl,Mach,Thot,Tcold,0.0d0,0.0d0]
#endif
        !===============================================================================================
        real(kind=8) :: errorU, errorT
        real(kind=8) :: xp(0:nx+1), yp(0:ny+1)
        real(kind=8), allocatable :: rho(:,:), u(:,:), v(:,:), T(:,:)
#ifdef steadyFlow
        real(kind=8), allocatable :: up(:,:), vp(:,:), Tp(:,:)
        real(kind=8) :: Nu_global, Nu_hot, Nu_cold, Nu_middle
        real(kind=8) :: Nu_hot_max, Nu_hot_min, Nu_hot_max_position, Nu_hot_min_position
#endif
        ! old/star 的 ghost 只用于格子端点求值；不进入输出或体平均。
        real(kind=8), allocatable :: rho_old(:,:), u_old(:,:), v_old(:,:), T_old(:,:)
        real(kind=8), allocatable :: rho_star(:,:), u_star(:,:), v_star(:,:), T_star(:,:)
        ! N_alpha(i,j) 位于 (r_i-e_alpha/2,t_n+1/2)，并非格心值。
        real(kind=8), allocatable :: f_neq(:,:,:), g_neq(:,:,:)
        integer(kind=4) :: itc
        logical, parameter :: useUSTLBM=.true.
        integer(kind=4) :: ex(0:8), ey(0:8)
        data ex/0,1,0,-1,0,1,-1,-1,1/
        data ey/0,0,1,0,-1,1,1,-1,-1/
        real(kind=8) :: omega(0:8)
    end module commondata
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
    ! 单进程只初始化默认设备（默认设备号为 0），不做 MPI 或多设备分配。
    call acc_init(acc_device_default)
    numAccDevices = acc_get_num_devices(acc_device_default)
    write(00,*) "Visible OpenACC devices:", numAccDevices
    close(00)
    !===============================================================================================


    !===============================================================================================
    ! Initialization
    call initial()
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
        
        call ustlbm_step()
        if(MOD(restartItcOffset+itc,stateCheckIntervalItc).EQ.0) call check_ustlbm_state_acc()

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
            call output_ReloadFile()      !稳态模式下的可选周期 rho/u/v/T 重启文件输出
        endif
#endif

#ifdef unsteadyFlow
        do while( dble(sampleIndexOffset+dimensionlessTime+1)*outputSnapshotInterval &
            .LE.unsteadyRunDuration+1.0d-12*max(1.0d0,unsteadyRunDuration) )
            ! 用整数采样序号固定全局时刻；重启发生在非采样格子步时也不平移后续时间轴。
            nextSampleAbsItc=max(1,int(dble(sampleIndexOffset+dimensionlessTime+1) &
                *outputSnapshotInterval*timeUnit+0.5d0))
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
            call output_ReloadFile()      !非稳态模式下的可选周期 rho/u/v/T 重启文件输出
        endif
#endif
     enddo

    !$acc wait(1)
    call CPU_TIME(timeEnd)         !当前进程累计消耗的 CPU 时间,包括并行
    ! 取墙钟结束计数器，用于后面输出 OpenACC 实际耗时。
    call system_clock(wallClockEnd, wallClockRate)
    timeEnd2 = dble(wallClockEnd) / dble(max(wallClockRate,1_8))
    call check_ustlbm_state_acc()
    call update_host_snapshot_2d_openacc()
    if( (outputReloadFile.EQ.1).AND.(MOD(restartItcOffset+itc, reloadFileIntervalItc).NE.0) ) then
        call output_ReloadFile()      !结束步不是周期输出点时，额外保存最终 rho/u/v/T 重启状态
    endif

#ifdef steadyFlow
    call output_Tecplot()          !输出最后一步的plt结果
    call output_SnapshotFile()              !输出最后一步的uvTrho数据
#endif

#ifdef unsteadyFlow
    call output_unsteady_NuRe_postprocess()
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


#ifdef VerticalWallsNoslip
    ! psi/vort 后处理默认封闭腔体：四周无滑移，psi 在物理边界取同一常数。
    ! 若垂直边界改为周期速度边界，流函数边界补点和涡量单边差分需要另写周期版本。
    call calc_psi_vort_and_output()  ! 输出中心abs(psi), max(abs(psi))及位置；max位置用细网格样条插值
#endif

#ifdef steadyFlow
    call calNuRe()
#endif




     

    open(unit=00,file=trim(settingsFile),status='unknown',position='append')        !在这个txt文件后面继续写（追加模式）
    write(00,*) "======================================================================"
    write(00,*) "Time (CPU) = ", real(timeEnd-timeStart,kind=8), "s"                             !当前进程累计消耗的 CPU 时间,包括并行
    if(timeEnd.GT.timeStart) then
        write(00,*) 'MLUPS =',dble(nx)*dble(ny)*dble(itc)/(timeEnd-timeStart)/1.0d6
    else
        write(00,*) 'MLUPS unavailable: CPU elapsed time is below timer resolution'
    endif
    write(00,*) "Time (ACC) = ", real(timeEnd2-timeStart2,kind=8), "s"                           !墙钟时间
    if(timeEnd2.GT.timeStart2) then
        write(00,*) 'MLUPS (ACC) =',dble(nx)*dble(ny)*dble(itc)/(timeEnd2-timeStart2)/1.0d6
    else
        write(00,*) 'MLUPS (ACC) unavailable: wall elapsed time is below timer resolution'
    endif
#ifdef steadyFlow
    write(00,'(a,1x,ES24.16E3)') "Nu_global =", real(Nu_global,kind=8)
    write(00,'(a,1x,ES24.16E3)') "Nu_hot    =", real(Nu_hot,kind=8)
    write(00,'(a,1x,ES24.16E3)') "Nu_cold   =", real(Nu_cold,kind=8)
#endif
    write(00,*) 'algorithm = USTLBM2017'
    write(00,*) 'total_itc =',restartItcOffset+itc
    write(00,*) 'actual_time_tf =',dble(restartItcOffset+itc)/timeUnit
#ifdef steadyFlow
    write(00,*) 'steadyConverged =',errorU.LE.epsU.AND.errorT.LE.epsT
    if(errorU.GT.epsU.OR.errorT.GT.epsT) write(00,*) 'Iteration limit reached before steady convergence.'
#endif



    write(00,*) 'Deallocate USTLBM arrays'
    call exit_data_2d_openacc()
    deallocate(rho,u,v,T,rho_old,u_old,v_old,T_old,rho_star,u_star,v_star,T_star,f_neq,g_neq)
#ifdef steadyFlow
    deallocate(up,vp,Tp)
#endif
    write(00,*) "    "
    
    write(00,*) "USTLBM run completed; check convergence and statistical diagnostics."
    
    string = ctime( time() )
    write(00,*) 'End:   ', string           !什么时候算完
    close(00)

    
    end program main


!===================================================================================================
! 子程序: initial
! 作用: 初始化项目半格网格、宏观整步状态及输出；续算只接受本版本的完整元数据。
!===================================================================================================
  subroutine initial()
    use commondata
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer(kind=4) :: i,j,alpha,ios,storedNx,storedNy
    character(len=100) :: reloadFileName
    character(len=16) :: storedSignature

    if(nx.LT.8.OR.ny.LT.8) error stop 'At least 8 cells per direction are required by wall-fit diagnostics'
    if(Rayleigh.LE.0.0d0.OR.Prandtl.LE.0.0d0.OR.Mach.LE.0.0d0.OR.temperatureDifference.LE.0.0d0) &
        error stop 'Ra, Pr, Mach and Thot-Tcold must be positive'
    if(loadInitField.NE.0.AND.loadInitField.NE.1) error stop 'loadInitField must be 0 or 1'
    if(steadyCheckIntervalItc.LE.0.OR.stateCheckIntervalItc.LE.0) error stop 'Check intervals must be positive'
    if(outputSnapshotInterval.LE.0.0d0.OR.reloadFileInterval.LE.0.0d0.OR.outputPltFileInterval.LE.0.0d0) &
        error stop 'Output intervals must be positive'
    itc=0
    errorU=100.0d0
    errorT=100.0d0
    snapshotFileNum=0
    pltFileNum=0
    restartItcOffset=0
    sampleIndexOffset=0
    reloadMetadataLoaded=.false.
    outputSnapshotIntervalItc=max(1,int(outputSnapshotInterval*timeUnit+0.5d0))
    reloadFileIntervalItc=max(1,int(reloadFileInterval*timeUnit+0.5d0))
    outputPltFileIntervalItc=max(1,int(outputPltFileInterval*timeUnit+0.5d0))
    xp(0)=0.0d0
    xp(nx+1)=dble(nx)/lengthUnit
    yp(0)=0.0d0
    yp(ny+1)=dble(ny)/lengthUnit
    do i=1,nx
        xp(i)=(dble(i)-0.5d0)/lengthUnit
    enddo
    do j=1,ny
        yp(j)=(dble(j)-0.5d0)/lengthUnit
    enddo
    allocate(rho(nx,ny),u(nx,ny),v(nx,ny),T(nx,ny))
    allocate(rho_old(0:nx+1,0:ny+1),u_old(0:nx+1,0:ny+1),v_old(0:nx+1,0:ny+1),T_old(0:nx+1,0:ny+1))
    allocate(rho_star(0:nx+1,0:ny+1),u_star(0:nx+1,0:ny+1),v_star(0:nx+1,0:ny+1),T_star(0:nx+1,0:ny+1))
    allocate(f_neq(0:nx+1,0:ny+1,0:8),g_neq(0:nx+1,0:ny+1,0:8))
#ifdef steadyFlow
    allocate(up(nx,ny),vp(nx,ny),Tp(nx,ny))
#endif
    omega(0)=4.0d0/9.0d0
    do alpha=1,4
        omega(alpha)=1.0d0/9.0d0
    enddo
    do alpha=5,8
        omega(alpha)=1.0d0/36.0d0
    enddo
    if(loadInitField.EQ.0) then
        if(reloadDimensionlessTime.NE.0.0d0) error stop 'New run requires reloadDimensionlessTime=0'
        rho=1.0d0
        u=0.0d0
        v=0.0d0
        do j=1,ny
            do i=1,nx
#ifdef SideHeatedCell
                T(i,j)=Thot+(Tcold-Thot)*xp(i)/xp(nx+1)
#else
                T(i,j)=Thot+(Tcold-Thot)*yp(j)/yp(ny+1)
                T(i,j)=T(i,j)+rbInitialPerturbation*temperatureDifference &
                    *cos(2.0d0*pi*xp(i)/xp(nx+1))*sin(pi*yp(j)/yp(ny+1))
#endif
            enddo
        enddo
        reloadFileNum=0
    else
        call read_reload_metadata(reloadFileName)
        open(unit=01,file=trim(reloadFilePrefix)//'-'//trim(reloadFileName)//'.bin', &
            form='unformatted',access='sequential',status='old',iostat=ios)
        if(ios.NE.0) error stop 'Cannot open USTLBM restart binary'
        read(01,iostat=ios) storedSignature,storedNx,storedNy
        if(ios.NE.0) error stop 'Invalid USTLBM binary header'
        if(storedSignature.NE.restartSignature.OR.storedNx.NE.nx.OR.storedNy.NE.ny) &
            error stop 'Restart algorithm or grid mismatch'
        read(01) rho
        read(01) u
        read(01) v
        read(01) T
        close(01)
        if(.not.all(ieee_is_finite(rho)).OR..not.all(ieee_is_finite(u)).OR. &
           .not.all(ieee_is_finite(v)).OR..not.all(ieee_is_finite(T)).OR.any(rho.LE.0.0d0)) &
            error stop 'Non-finite or non-positive-density restart state'
        call align_ustlbm_restart_history() ! 先对齐 checkpoint 边界，再继续写 Nu/Re 历史
    endif
#ifdef steadyFlow
    up=u
    vp=v
    Tp=T
#endif
    dimensionlessTime=0
    NuVolAvg=0.0d0
    ReVolAvg=0.0d0
    open(unit=00,file=trim(settingsFile),status='unknown',position='append')
    write(00,*) 'algorithm = USTLBM2017, D2Q9 momentum + D2Q9 internal energy'
    write(00,*) 'Reference: Chen, Shu, Tan, Applied Sciences 7 (2017) 277, Eqs.39-47'
    write(00,*) 'caseConfiguration =',caseConfiguration
    write(00,*) 'Mesh:',nx,ny,'; cell centers; physical walls at half a lattice spacing'
    write(00,*) 'Rayleigh=',Rayleigh,'; Prandtl=',Prandtl,'; Mach=',Mach
    write(00,*) 'Thot=',Thot,'; Tcold=',Tcold,'; gBeta=',gBeta
    write(00,*) 'viscosity=',viscosity,'; diffusivity=',diffusivity
    write(00,*) 'tauf=',tauf,'; tauT=',tauT,'; thermal second moment=2/3'
    write(00,*) 'timeUnit=',timeUnit,'; velocityUnit=',velocityUnit
    write(00,*) 'velocityScaleCompare=',velocityScaleCompare
    write(00,*) 'outputSnapshotInterval=',outputSnapshotInterval,'; reloadFileInterval=',reloadFileInterval
    write(00,*) 'outputPltFileInterval=',outputPltFileInterval,'; itc_max=',itc_max
    write(00,*) 'epsU=',epsU,'; epsT=',epsT
    write(00,*) 'restartItcOffset=',restartItcOffset,'; reloadDimensionlessTime=',reloadDimensionlessTime
#ifdef unsteadyFlow
    write(00,*) 'unsteadyRunDuration=',unsteadyRunDuration,'; rbInitialPerturbation=',rbInitialPerturbation
#endif
    write(00,*) 'No unconditional-stability guarantee; convergence and benchmark accuracy are separate checks.'
    close(00)
  end subroutine initial
!===================================================================================================

!===================================================================================================
! 子程序: enter_data_2d_openacc / exit_data_2d_openacc
! 作用: 整步宏观状态和半格非平衡缓存常驻单个设备，所有推进核使用 async(1)。
!===================================================================================================
  subroutine enter_data_2d_openacc()
    use commondata
    implicit none
    !$acc enter data copyin(xp,yp,ex,ey,omega,rho,u,v,T)
    !$acc enter data create(rho_old,u_old,v_old,T_old,rho_star,u_star,v_star,T_star,f_neq,g_neq)
#ifdef steadyFlow
    !$acc enter data copyin(up,vp,Tp)
#endif
  end subroutine enter_data_2d_openacc

  subroutine exit_data_2d_openacc()
    use commondata
    implicit none
    !$acc wait(1)
#ifdef steadyFlow
    !$acc exit data delete(up,vp,Tp)
#endif
    !$acc exit data delete(rho_old,u_old,v_old,T_old,rho_star,u_star,v_star,T_star,f_neq,g_neq)
    !$acc exit data delete(xp,yp,ex,ey,omega,rho,u,v,T)
  end subroutine exit_data_2d_openacc
!===================================================================================================

!===================================================================================================
! 子程序: update_host_snapshot_2d_openacc / update_host_tecplot_2d_openacc / update_host_reload_2d_openacc
! 作用: CPU 输出或最终后处理前，只同步所需的整步宏观量。
!===================================================================================================
  subroutine update_host_snapshot_2d_openacc()
    use commondata
    implicit none
    !$acc wait(1)
    !$acc update self(u,v,T,rho)
  end subroutine update_host_snapshot_2d_openacc

  subroutine update_host_tecplot_2d_openacc()
    use commondata
    implicit none
    !$acc wait(1)
    !$acc update self(u,v,T)
  end subroutine update_host_tecplot_2d_openacc

  subroutine update_host_reload_2d_openacc()
    use commondata
    implicit none
    call check_ustlbm_state_acc() ! 写 checkpoint 前拒绝非有限场，避免生成不可用的续算状态
    !$acc wait(1)
    !$acc update self(rho,u,v,T)
  end subroutine update_host_reload_2d_openacc
!===================================================================================================

!===================================================================================================
! 子程序: ustlbm_equilibrium
! 作用: 由格子方向与宏观变量计算单方向平衡态；全部输入为标量，避免设备隐式捕获局部数组。
! 二维内能 e=T (R=1)，sum(g_eq)=rho*T、sum(e_alpha*g_eq)=rho*T*u。
! 零速 sum(ex_alpha**2*g_eq)=(2/3)*rho*T；因此 kappa=(2/3)*(tauT-0.5)。
! g_eq(0) 可以为负，这是论文内能平衡态的定义，不能逐方向做正性裁剪。
!===================================================================================================
  subroutine ustlbm_equilibrium(alpha,exAlpha,eyAlpha,weight,rho0,u0,v0,temp0,feq0,geq0)
    implicit none
    !$acc routine seq
    integer(kind=4), intent(in) :: alpha,exAlpha,eyAlpha
    real(kind=8), intent(in) :: weight,rho0,u0,v0,temp0
    real(kind=8), intent(out) :: feq0,geq0
    real(kind=8) :: eu,u2,rhoE
    eu=dble(exAlpha)*u0+dble(eyAlpha)*v0
    u2=u0*u0+v0*v0
    rhoE=rho0*temp0
    feq0=weight*rho0*(1.0d0+3.0d0*eu+4.5d0*eu*eu-1.5d0*u2)
    if(alpha.EQ.0) then
        geq0=-2.0d0*rhoE*u2/3.0d0
    elseif(alpha.LE.4) then
        geq0=rhoE*(1.0d0/6.0d0+eu/6.0d0+0.5d0*eu*eu-u2/6.0d0)
    else
        geq0=rhoE*(1.0d0/12.0d0+eu/6.0d0+0.125d0*eu*eu-u2/24.0d0)
    endif
  end subroutine ustlbm_equilibrium
!===================================================================================================

!===================================================================================================
! 子程序: build_ustlbm_old_ghosts_acc
! 作用: 延拓上一个整步的宏观量，供预测步和半格非平衡重构使用。
! old/star 的索引 0 是墙外格心；xp(0)/yp(0) 则仍是后处理物理墙坐标，二者不能混用。
! 同一 ghost 仅从真实内点读入；角点先处理 x，再处理 y，避免原地边界写读竞争。
!===================================================================================================
  subroutine build_ustlbm_old_ghosts_acc()
    use commondata
    implicit none
    integer(kind=4) :: i,j,ic,jc
    real(kind=8) :: rho0,u0,v0,temp0
    !$acc parallel loop gang vector collapse(2) default(none) &
    !$acc& present(rho,u,v,T,rho_old,u_old,v_old,T_old) async(1) private(ic,jc,rho0,u0,v0,temp0)
    do j=0,ny+1
        do i=0,nx+1
            ic=min(max(i,1),nx)
            jc=min(max(j,1),ny)
#ifdef VerticalWallsPeriodicalU
            if(i.EQ.0) ic=nx
            if(i.EQ.nx+1) ic=1
#endif
            rho0=rho(ic,jc)
            u0=u(ic,jc)
            v0=v(ic,jc)
            temp0=T(ic,jc)
#ifdef VerticalWallsNoslip
            if(i.EQ.0.OR.i.EQ.nx+1) then
                u0=-u0
                v0=-v0
            endif
#endif
#ifdef VerticalWallsConstT
            if(i.EQ.0) temp0=2.0d0*Thot-temp0
            if(i.EQ.nx+1) temp0=2.0d0*Tcold-temp0
#endif
            ! 绝热温度和密度采用偶延拓；无滑移速度采用奇延拓。
            if(j.EQ.0.OR.j.EQ.ny+1) then
                u0=-u0
                v0=-v0
            endif
#ifdef HorizontalWallsConstT
            if(j.EQ.0) temp0=2.0d0*Thot-temp0
            if(j.EQ.ny+1) temp0=2.0d0*Tcold-temp0
#endif
            rho_old(i,j)=rho0
            u_old(i,j)=u0
            v_old(i,j)=v0
            T_old(i,j)=temp0
        enddo
    enddo
  end subroutine build_ustlbm_old_ghosts_acc
!===================================================================================================

!===================================================================================================
! 子程序: build_ustlbm_star_ghosts_acc
! 作用: 对已完成预测的宏观量施加同一半格墙延拓；只写 ghost，所有来源都是真实内点。
!===================================================================================================
  subroutine build_ustlbm_star_ghosts_acc()
    use commondata
    implicit none
    integer(kind=4) :: i,j,ic,jc
    real(kind=8) :: rho0,u0,v0,temp0
    !$acc parallel loop gang vector collapse(2) default(none) present(rho_star,u_star,v_star,T_star) &
    !$acc& async(1) private(ic,jc,rho0,u0,v0,temp0)
    do j=0,ny+1
        do i=0,nx+1
            if(i.EQ.0.OR.i.EQ.nx+1.OR.j.EQ.0.OR.j.EQ.ny+1) then
                ic=min(max(i,1),nx)
                jc=min(max(j,1),ny)
#ifdef VerticalWallsPeriodicalU
                if(i.EQ.0) ic=nx
                if(i.EQ.nx+1) ic=1
#endif
                rho0=rho_star(ic,jc)
                u0=u_star(ic,jc)
                v0=v_star(ic,jc)
                temp0=T_star(ic,jc)
#ifdef VerticalWallsNoslip
                if(i.EQ.0.OR.i.EQ.nx+1) then
                    u0=-u0
                    v0=-v0
                endif
#endif
#ifdef VerticalWallsConstT
                if(i.EQ.0) temp0=2.0d0*Thot-temp0
                if(i.EQ.nx+1) temp0=2.0d0*Tcold-temp0
#endif
                if(j.EQ.0.OR.j.EQ.ny+1) then
                    u0=-u0
                    v0=-v0
                endif
#ifdef HorizontalWallsConstT
                if(j.EQ.0) temp0=2.0d0*Thot-temp0
                if(j.EQ.ny+1) temp0=2.0d0*Tcold-temp0
#endif
                rho_star(i,j)=rho0
                u_star(i,j)=u0
                v_star(i,j)=v0
                T_star(i,j)=temp0
            endif
        enddo
    enddo
  end subroutine build_ustlbm_star_ghosts_acc
!===================================================================================================

!===================================================================================================
! 子程序: ustlbm_step
! 作用: 按论文式 (39)--(47) 完成一个整步，dt=dx=1。
! 1) 预测平衡态上游求和得到 rhoStar、rhoStar*uStar、rhoStar*TStar。
! 2) N_alpha(P)=-tau*[EeqStar(P)-EeqOld(P-e_alpha)]，存于半格 P-e_alpha/2。
! 3) rho=rhoStar；动量/内能校正为 -(1-1/tau)*sum[N(P+e_alpha)-N(P)]。
! 4) 先得到 Tnew，再用 0.5*rhoStar*gBeta*(Told+Tnew-2*Tref) 加完整一步浮力。
! 本式已更新宏观动量，不能再额外叠加传统 Guo 宏观速度的半步力修正。
!===================================================================================================
  subroutine ustlbm_step()
    use commondata
    implicit none
    integer(kind=4) :: i,j,alpha,ip,jp
    real(kind=8) :: feq0,geq0,rhoSum,momxSum,momySum,energySum
    real(kind=8) :: feqStar,geqStar,deltaF,deltaG,corrX,corrY,corrEnergy,energyNew,tempNew,velx,vely
#ifdef SideHeatedHa
    real(kind=8) :: dxx,dxy,dyy,rhsx,rhsy,a11,a12,a22,det
#endif

    call build_ustlbm_old_ghosts_acc()
    !$acc parallel loop gang vector collapse(2) default(none) &
    !$acc& present(rho_old,u_old,v_old,T_old,rho_star,u_star,v_star,T_star,ex,ey,omega) &
    !$acc& async(1) private(alpha,ip,jp,feq0,geq0,rhoSum,momxSum,momySum,energySum)
    do j=1,ny
        do i=1,nx
            rhoSum=0.0d0
            momxSum=0.0d0
            momySum=0.0d0
            energySum=0.0d0
            !$acc loop seq
            do alpha=0,8
                ip=i-ex(alpha)
                jp=j-ey(alpha)
                call ustlbm_equilibrium(alpha,ex(alpha),ey(alpha),omega(alpha), &
                    rho_old(ip,jp),u_old(ip,jp),v_old(ip,jp),T_old(ip,jp),feq0,geq0)
                rhoSum=rhoSum+feq0
                momxSum=momxSum+dble(ex(alpha))*feq0
                momySum=momySum+dble(ey(alpha))*feq0
                energySum=energySum+geq0
            enddo
            rho_star(i,j)=rhoSum
            u_star(i,j)=momxSum/rhoSum
            v_star(i,j)=momySum/rhoSum
            T_star(i,j)=energySum/rhoSum
        enddo
    enddo
    call build_ustlbm_star_ghosts_acc()

    !$acc parallel loop gang vector collapse(2) default(none) &
    !$acc& present(rho_old,u_old,v_old,T_old,rho_star,u_star,v_star,T_star,f_neq,g_neq,ex,ey,omega) &
    !$acc& async(1) private(alpha,ip,jp,feq0,geq0,feqStar,geqStar)
    do j=0,ny+1
        do i=0,nx+1
            !$acc loop seq
            do alpha=0,8
                ip=i-ex(alpha)
                jp=j-ey(alpha)
                f_neq(i,j,alpha)=0.0d0
                g_neq(i,j,alpha)=0.0d0
                ! corrector 只会读取源端点有效的方向；无效槽先置零且不做越界读取。
                if(ip.GE.0.AND.ip.LE.nx+1.AND.jp.GE.0.AND.jp.LE.ny+1) then
                    call ustlbm_equilibrium(alpha,ex(alpha),ey(alpha),omega(alpha), &
                        rho_star(i,j),u_star(i,j),v_star(i,j),T_star(i,j),feqStar,geqStar)
                    call ustlbm_equilibrium(alpha,ex(alpha),ey(alpha),omega(alpha), &
                        rho_old(ip,jp),u_old(ip,jp),v_old(ip,jp),T_old(ip,jp),feq0,geq0)
                    f_neq(i,j,alpha)=-tauf*(feqStar-feq0)
                    g_neq(i,j,alpha)=-tauT*(geqStar-geq0)
                endif
            enddo
        enddo
    enddo

    !$acc parallel loop gang vector collapse(2) default(none) &
    !$acc& present(rho,u,v,T,rho_old,u_old,v_old,T_old,rho_star,u_star,v_star,T_star,f_neq,g_neq,ex,ey) &
    !$acc& async(1) &
#ifdef SideHeatedHa
    !$acc& private(dxx,dxy,dyy,rhsx,rhsy,a11,a12,a22,det) &
#endif
    !$acc& private(alpha,ip,jp,deltaF,deltaG,corrX,corrY,corrEnergy,energyNew,tempNew,velx,vely)
    do j=1,ny
        do i=1,nx
            corrX=0.0d0
            corrY=0.0d0
            corrEnergy=0.0d0
            !$acc loop seq
            do alpha=0,8
                ip=i+ex(alpha)
                jp=j+ey(alpha)
                deltaF=f_neq(ip,jp,alpha)-f_neq(i,j,alpha)
                deltaG=g_neq(ip,jp,alpha)-g_neq(i,j,alpha)
                corrX=corrX+dble(ex(alpha))*deltaF
                corrY=corrY+dble(ey(alpha))*deltaF
                corrEnergy=corrEnergy+deltaG
            enddo
            energyNew=rho_star(i,j)*T_star(i,j)-(1.0d0-1.0d0/tauT)*corrEnergy
            tempNew=energyNew/rho_star(i,j)
            velx=u_star(i,j)-(1.0d0-1.0d0/tauf)*corrX/rho_star(i,j)
            vely=v_star(i,j)-(1.0d0-1.0d0/tauf)*corrY/rho_star(i,j) &
                +0.5d0*gBeta*(T_old(i,j)+tempNew-2.0d0*Tref)
#ifdef SideHeatedHa
            ! 保留原可选线性磁阻力；用梯形积分解局部 2x2 系统，避免把二阶推进降回显式一阶。
            dxx=-B2sigemarho*sin(phi)**2
            dxy=B2sigemarho*sin(phi)*cos(phi)
            dyy=-B2sigemarho*cos(phi)**2
            rhsx=velx+0.5d0*(dxx*u_old(i,j)+dxy*v_old(i,j))
            rhsy=vely+0.5d0*(dxy*u_old(i,j)+dyy*v_old(i,j))
            a11=1.0d0-0.5d0*dxx
            a12=-0.5d0*dxy
            a22=1.0d0-0.5d0*dyy
            det=a11*a22-a12*a12
            velx=(a22*rhsx-a12*rhsy)/det
            vely=(a11*rhsy-a12*rhsx)/det
#endif
            rho(i,j)=rho_star(i,j)
            u(i,j)=velx
            v(i,j)=vely
            T(i,j)=tempNew
        enddo
    enddo
  end subroutine ustlbm_step
!===================================================================================================

!===================================================================================================
! 子程序: check_ustlbm_state_acc
! 作用: 检查非有限数和非正密度；错误时终止，不把坏密度悄悄重置为 1。
! 封闭墙 ghost 不等于逐链反弹，另记录密度范围和平均值以观察质量漂移。
!===================================================================================================
  subroutine check_ustlbm_state_acc()
    use commondata
    implicit none
    integer(kind=4) :: i,j,badCount
    real(kind=8) :: rhoSum,rhoMin,rhoMax,speedMax
    !$acc wait(1)
    badCount=0
    rhoSum=0.0d0
    rhoMin=huge(1.0d0)
    rhoMax=-huge(1.0d0)
    speedMax=0.0d0
    !$acc parallel loop collapse(2) default(none) present(rho,u,v,T) &
    !$acc& reduction(+:badCount,rhoSum) reduction(min:rhoMin) reduction(max:rhoMax,speedMax)
    do j=1,ny
        do i=1,nx
            if(.not.(rho(i,j).GT.0.0d0.AND.rho(i,j).LT.huge(1.0d0))) badCount=badCount+1
            if(.not.(abs(u(i,j)).LE.huge(1.0d0))) badCount=badCount+1
            if(.not.(abs(v(i,j)).LE.huge(1.0d0))) badCount=badCount+1
            if(.not.(abs(T(i,j)).LE.huge(1.0d0))) badCount=badCount+1
            rhoSum=rhoSum+rho(i,j)
            rhoMin=min(rhoMin,rho(i,j))
            rhoMax=max(rhoMax,rho(i,j))
            speedMax=max(speedMax,sqrt(u(i,j)**2+v(i,j)**2))
        enddo
    enddo
    if(badCount.GT.0) then
        write(*,*) 'USTLBM invalid state at total_itc=',restartItcOffset+itc,'; invalid entries=',badCount
        error stop 'Non-finite field or non-positive density'
    endif
    open(unit=00,file=trim(settingsFile),status='unknown',position='append')
    write(00,'(A,I0,4(1X,ES24.16E3))') 'state_check itc/rho_mean/rho_min/rho_max/maxMa: ', &
        restartItcOffset+itc,rhoSum/dble(nx*ny),rhoMin,rhoMax,sqrt(3.0d0)*speedMax
    close(00)
  end subroutine check_ustlbm_state_acc
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

  

    call append_convergence_tecplot('convergence2DOpenaccUSTLBM.plt', restartItcOffset+itc, errorU, errorT)

    
    write(caseTag,'("Ra=",ES10.3E2,",nx=",I0,",ny=",I0,",USTLBM2017")') Rayleigh,nx,ny
    call append_convergence_master_tecplot('convergence_all_2DOpenaccUSTLBM.plt', caseTag, restartItcOffset+itc, errorU, errorT)

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
      error stop 1
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
        error stop 1
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
! 作用: 输出 USTLBM 重启备份文件；直接记录宏观 rho/u/v/T。
! 用途: 在运行过程中定期调用，也在程序结束前调用。
!===================================================================================================
  subroutine output_ReloadFile()                                  !输出 reload 文件，名字由 reloadFilePrefix 控制
    use commondata                                                !USTLBM 严格重启直接读回 rho/u/v/T
    implicit none
    integer(kind=4) :: i, j
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
    ! USTLBM 不演化分布函数；严格重启保存宏观状态量。
    write(05) restartSignature,nx,ny
    write(05) ((real(rho(i,j),kind=8), i=1,nx), j=1,ny)
    write(05) ((real(u(i,j),kind=8), i=1,nx), j=1,ny)
    write(05) ((real(v(i,j),kind=8), i=1,nx), j=1,ny)
    write(05) ((real(T(i,j),kind=8), i=1,nx), j=1,ny)
    close(05)
    call write_reload_metadata(trim(filename))
    
    open(unit=00,file=trim(settingsFile),status='unknown',position='append')
    write(00,*) "Backup USTLBM rho/u/v/T restart state to: ", trim(reloadFilePrefix), "-", trim(filename),".bin"
    write(00,*) "Backup restart metadata to: ", trim(reloadFilePrefix), "-latest.meta"
    close(00)
    
    return
  end subroutine output_ReloadFile
!===================================================================================================
! output_ReloadFile 结束: 输出 USTLBM 严格重启备份文件。
!===================================================================================================


!===================================================================================================
! 子程序: write_reload_metadata
! 作用: 覆盖写出最新 reload 续算账本，恢复累计步数、t_ff、输出编号和最新 .bin 文件名。
!===================================================================================================
  subroutine write_reload_metadata(filename)
    use commondata
    implicit none
    character(len=*), intent(in) :: filename
    integer(kind=4) :: metaUnit, totalItc
    real(kind=8) :: totalTf

    totalItc = restartItcOffset + itc
    totalTf = real(totalItc,kind=8) / timeUnit

    open(newunit=metaUnit, file=trim(reloadFilePrefix)//'-latest.meta', &
         status='replace', action='write', form='formatted')
    write(metaUnit,'(A,1X,I0)') 'reload_meta_version', 4
    write(metaUnit,'(A,1X,A)') 'algorithm',restartSignature
    write(metaUnit,'(A,1X,A)') 'caseConfiguration',caseConfiguration
    write(metaUnit,'(A,7(1X,ES24.16E3))') 'physical_parameters',restartPhysics
#ifdef steadyFlow
    write(metaUnit,'(A,1X,A)') 'flowMode', 'steadyFlow'
#endif
#ifdef unsteadyFlow
    write(metaUnit,'(A,1X,A)') 'flowMode', 'unsteadyFlow'
#endif
    write(metaUnit,'(A,1X,I0)') 'nx', nx
    write(metaUnit,'(A,1X,I0)') 'ny', ny
    write(metaUnit,'(A,1X,A)') 'reloadFileName', trim(filename)
    write(metaUnit,'(A,1X,I0)') 'itc_total', totalItc
    write(metaUnit,'(A,1X,ES24.16E3)') 'time_tf', totalTf
    write(metaUnit,'(A,1X,I0)') 'snapshotFileNum', snapshotFileNum
    write(metaUnit,'(A,1X,I0)') 'pltFileNum', pltFileNum
    write(metaUnit,'(A,1X,I0)') 'reloadFileNum', reloadFileNum
    write(metaUnit,'(A,1X,I0)') 'sampleIndex', sampleIndexOffset+dimensionlessTime
    write(metaUnit,'(A,1X,ES24.16E3)') 'sampling_interval', outputSnapshotInterval
    close(metaUnit)

    return
  end subroutine write_reload_metadata
!===================================================================================================


!===================================================================================================
! 子程序: read_reload_metadata
! 作用: 严格读取本算法 latest.meta，核对物理参数、边界和采样网格；禁止推断旧格式。
!===================================================================================================
  subroutine read_reload_metadata(reloadFileName)
    use commondata
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    character(len=*), intent(inout) :: reloadFileName
    character(len=64) :: label
    character(len=32) :: metaFlowMode, currentFlowMode
    character(len=100) :: metaReloadFileName
    character(len=256) :: metaFile
    integer(kind=4) :: metaUnit, ios
    integer(kind=4) :: metaVersion, metaNx, metaNy
    integer(kind=4) :: metaItc, metaSnapshotFileNum, metaPltFileNum, metaReloadFileNum
    real(kind=8) :: metaTf
    logical :: metaExists
    character(len=32) :: metaAlgorithm,metaCase
    real(kind=8) :: metaPhysics(7)
    real(kind=8) :: metaSamplingInterval
    integer(kind=4) :: metaSampleIndex

    reloadMetadataLoaded = .false.
    metaFile = trim(reloadFilePrefix)//'-latest.meta'
    inquire(file=trim(metaFile), exist=metaExists)                 !优先检查最新账本

    if(.not.metaExists) then
        error stop 'USTLBM restart requires its latest.meta; inferred/old STLBM restarts are not accepted'
        return                                                     !退出当前这个子程序
    endif

    open(newunit=metaUnit, file=trim(metaFile), status='old', action='read', &
         form='formatted', iostat=ios)                             !ios == 0表示读写成功,不等于0，则代表失败
    if(ios.NE.0) then
        write(*,*) 'Error: failed to open reload metadata: ', trim(metaFile)
        error stop 1
    endif

    read(metaUnit,*,iostat=ios) label, metaVersion
    if((ios.NE.0).OR.(trim(label).NE.'reload_meta_version').OR.(metaVersion.NE.4)) then
        write(*,*) 'Error: invalid reload metadata version in ', trim(metaFile)
        error stop 'Unsupported USTLBM metadata version'
    endif

    read(metaUnit,*,iostat=ios) label,metaAlgorithm
    if(ios.NE.0.OR.trim(label).NE.'algorithm') error stop 'Missing USTLBM algorithm tag'
    if(trim(metaAlgorithm).NE.trim(restartSignature)) error stop 'Restart algorithm differs'
    read(metaUnit,*,iostat=ios) label,metaCase
    if(ios.NE.0.OR.trim(label).NE.'caseConfiguration') error stop 'Missing USTLBM BC tag'
    if(trim(metaCase).NE.trim(caseConfiguration)) error stop 'Restart boundary configuration differs'
    read(metaUnit,*,iostat=ios) label,metaPhysics
    if(ios.NE.0.OR.trim(label).NE.'physical_parameters') error stop 'Missing USTLBM physical parameters'
    if(.not.all(ieee_is_finite(metaPhysics))) error stop 'Non-finite restart parameters'
    if(any(abs(metaPhysics-restartPhysics).GT.1.0d-12*max(1.0d0,abs(restartPhysics)))) &
        error stop 'Restart physical parameters differ'

    read(metaUnit,*,iostat=ios) label, metaFlowMode
    if((ios.NE.0).OR.(trim(label).NE.'flowMode')) then
        write(*,*) 'Error: invalid flowMode entry in ', trim(metaFile)
        error stop 1
    endif

    read(metaUnit,*,iostat=ios) label, metaNx
    if((ios.NE.0).OR.(trim(label).NE.'nx')) then
        write(*,*) 'Error: invalid nx entry in ', trim(metaFile)
        error stop 1
    endif

    read(metaUnit,*,iostat=ios) label, metaNy
    if((ios.NE.0).OR.(trim(label).NE.'ny')) then
        write(*,*) 'Error: invalid ny entry in ', trim(metaFile)
        error stop 1
    endif

    read(metaUnit,*,iostat=ios) label, metaReloadFileName
    if((ios.NE.0).OR.(trim(label).NE.'reloadFileName')) then
        write(*,*) 'Error: invalid reloadFileName entry in ', trim(metaFile)
        error stop 1
    endif
    metaReloadFileName = adjustl(metaReloadFileName)   ! 处理左边空格，把内容左对齐
    
    read(metaUnit,*,iostat=ios) label, metaItc
    if((ios.NE.0).OR.(trim(label).NE.'itc_total')) then
        write(*,*) 'Error: invalid itc_total entry in ', trim(metaFile)
        error stop 1
    endif

    read(metaUnit,*,iostat=ios) label, metaTf
    if((ios.NE.0).OR.(trim(label).NE.'time_tf')) then
        write(*,*) 'Error: invalid time_tf entry in ', trim(metaFile)
        error stop 1
    endif

    read(metaUnit,*,iostat=ios) label, metaSnapshotFileNum
    if((ios.NE.0).OR.(trim(label).NE.'snapshotFileNum')) then
        write(*,*) 'Error: invalid snapshotFileNum entry in ', trim(metaFile)
        error stop 1
    endif

    read(metaUnit,*,iostat=ios) label, metaPltFileNum
    if((ios.NE.0).OR.(trim(label).NE.'pltFileNum')) then
        write(*,*) 'Error: invalid pltFileNum entry in ', trim(metaFile)
        error stop 1
    endif

    read(metaUnit,*,iostat=ios) label, metaReloadFileNum
    if((ios.NE.0).OR.(trim(label).NE.'reloadFileNum')) then
        write(*,*) 'Error: invalid reloadFileNum entry in ', trim(metaFile)
        error stop 1
    endif
    read(metaUnit,*,iostat=ios) label,metaSampleIndex
    if(ios.NE.0.OR.trim(label).NE.'sampleIndex') error stop 'Missing absolute sample index'
    read(metaUnit,*,iostat=ios) label,metaSamplingInterval
    if(ios.NE.0.OR.trim(label).NE.'sampling_interval') error stop 'Missing sampling interval'
    close(metaUnit)

    if(.not.ieee_is_finite(metaTf).OR..not.ieee_is_finite(metaSamplingInterval)) &
        error stop 'Non-finite restart time or sampling interval'
    if(min(metaItc,metaSnapshotFileNum,metaPltFileNum,metaReloadFileNum,metaSampleIndex).LT.0) &
        error stop 'Negative restart counter'
    if(abs(metaTf-dble(metaItc)/timeUnit).GT.1.0d-12*max(1.0d0,abs(metaTf))) &
        error stop 'Inconsistent restart step and physical time'
#ifdef unsteadyFlow
    if(abs(metaSamplingInterval-outputSnapshotInterval).GT.1.0d-12*outputSnapshotInterval) &
        error stop 'Unsteady restart requires the same sampling interval'
#endif

    currentFlowMode = 'unknown'
#ifdef steadyFlow
    currentFlowMode = 'steadyFlow'
#endif
#ifdef unsteadyFlow
    currentFlowMode = 'unsteadyFlow'
#endif

    if(trim(metaFlowMode).NE.trim(currentFlowMode)) then
        write(*,*) 'Error: reload metadata flowMode differs: ', trim(metaFlowMode), trim(currentFlowMode)
        error stop 1
    endif
    if((metaNx.NE.nx).OR.(metaNy.NE.ny)) then
        write(*,*) 'Error: reload metadata mesh mismatch: ', metaNx, metaNy, nx, ny
        error stop 1
    endif

    restartItcOffset = metaItc
    reloadDimensionlessTime = metaTf
    snapshotFileNum = metaSnapshotFileNum
    pltFileNum = metaPltFileNum
    ! reloadFileNum 是整数计数器，给后续 output_ReloadFile() 继续编号，避免覆盖旧 reload 文件。
    ! reloadFileName 是字符串文件名，本次续算马上用它打开 <reloadFilePrefix>-<reloadFileName>.bin。
    reloadFileNum = metaReloadFileNum
    sampleIndexOffset=metaSampleIndex
    reloadFileName = trim(metaReloadFileName)     ! 处理右边空格，删掉尾部空格
    reloadMetadataLoaded = .true.

    return
  end subroutine read_reload_metadata
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
        error stop 1
    endif

    dimensionlessTime = dimensionlessTime+1   !每隔 outputSnapshotInterval 个自由落体时间调用一次calNuRe
#ifdef steadyFlow
    sampleTime = real(restartItcOffset+itc,kind=8)
#endif
#ifdef unsteadyFlow
    sampleTime = dble(sampleIndexOffset+dimensionlessTime)*outputSnapshotInterval
#endif

#ifdef unsteadyFlow
    ! 非稳态历史在重启前必须已对齐；稳态只在结束时写 Nu/Re，不能要求它们预先存在。
    if((first_nure_write).AND.(loadInitField.EQ.1)) then
        inquire(file="Nu_VolAvg_2DOpenaccUSTLBM.dat", exist=exNu)
        inquire(file="Re_VolAvg_2DOpenaccUSTLBM.dat", exist=exRe)
        if((.not.exNu).OR.(.not.exRe)) then
            write(*,*) "Error: restart requested but old Nu/Re time-series files are missing."
            open(unit=00,file=trim(settingsFile),status="unknown",position="append")
            write(00,*) "Error: restart requested but old Nu/Re time-series files are missing."
            write(00,*) "Nu_VolAvg_2DOpenaccUSTLBM.dat exists =", exNu
            write(00,*) "Re_VolAvg_2DOpenaccUSTLBM.dat exists =", exRe
            close(00)
            error stop 1
        endif
    endif
#endif


    
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


    ! Nu 的热流标度固定为 L/(kappa*DeltaT)，与快照速度的比较标度无关。
    NuVolAvg(dimensionlessTime)=1.0d0+NuVolAvg_temp*lengthUnit &
        /(dble(nx*ny)*diffusivity*temperatureDifference)

    if((first_nure_write).AND.(loadInitField.EQ.0)) then
        open(unit=01,file="Nu_VolAvg_2DOpenaccUSTLBM.dat",status='replace',action='write')
    else
        open(unit=01,file="Nu_VolAvg_2DOpenaccUSTLBM.dat",status='unknown',position='append',action='write')
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


    if((first_nure_write).AND.(loadInitField.EQ.0)) then
        open(unit=02,file="Re_VolAvg_2DOpenaccUSTLBM.dat",status='replace',action='write')
    else
        open(unit=02,file="Re_VolAvg_2DOpenaccUSTLBM.dat",status='unknown',position='append',action='write')
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

    inquire(file='Nu_VolAvg_2DOpenaccUSTLBM.dat', exist=exNu)
    inquire(file='Re_VolAvg_2DOpenaccUSTLBM.dat', exist=exRe)
    if((.not.exNu).or.(.not.exRe)) then
        write(*,'(A)') 'Error: Nu/Re history files are missing before postprocessing.'
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,'(A)') 'Error: Nu/Re history files are missing before postprocessing.'
        close(00)
        error stop 1
    endif

    open(newunit=nuUnit, file='Nu_VolAvg_2DOpenaccUSTLBM.dat', status='old', action='read', form='formatted')
    open(newunit=reUnit, file='Re_VolAvg_2DOpenaccUSTLBM.dat', status='old', action='read', form='formatted')

    ! These files are derived views of the full .dat history, so rebuild one continuous ZONE.
    open(newunit=seriesUnit, file='NuRe_VolAvg_2DOpenaccUSTLBM.plt', status='replace', action='write', form='formatted')
    write(seriesUnit,'(A)') 'TITLE = "2D OpenACC USTLBM Nu/Re volume averages"'
    write(seriesUnit,'(A)') 'VARIABLES = "time" "NuVolAvg" "ReVolAvg"'
    write(seriesUnit,'(A)') 'ZONE T="NuReVolAvg", F=POINT'

    open(newunit=runningUnit, file='NuRe_VolAvg_runningMean_2DOpenaccUSTLBM.plt', &
        status='replace', action='write', form='formatted')
    write(runningUnit,'(A)') 'TITLE = "2D OpenACC USTLBM Nu/Re running means"'
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

    do k = 1, unsteadySampleCount
        read(nuUnit,*,iostat=iosNu) timeNu, NuVal
        read(reUnit,*,iostat=iosRe) timeRe, ReVal
        ! iostat: 0=成功读到一行；小于0=到达文件末尾；大于0=格式或读入错误。
        ! 循环内只要不是 0，就说明文件短于 unsteadySampleCount，或者某一行格式坏了。
        if((iosNu.NE.0).OR.(iosRe.NE.0)) then
            write(*,'(A)') 'Error: Nu/Re history files are shorter than unsteadySampleCount or contain invalid rows.'
            open(unit=00,file=trim(settingsFile),status='unknown',position='append')
            write(00,'(A)') 'Error: Nu/Re history files are shorter than unsteadySampleCount or contain invalid rows.'
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

    ! 上面的循环已经读完预期的 unsteadySampleCount 行；这里再试读一行，
    ! 不是为了继续计算，而是确认 Nu/Re 两个历史文件后面没有多余数据。
    read(nuUnit,*,iostat=iosNu) timeNu, NuVal
    read(reUnit,*,iostat=iosRe) timeRe, ReVal
    ! 任意一个 iostat 等于 0，都表示至少一个文件还成功读到了额外一行。
    ! 如果放过这种情况，后处理会静默丢掉超过 unsteadySampleCount 的尾部样本。
    if((iosNu.EQ.0).OR.(iosRe.EQ.0)) then
        write(*,'(A)') 'Error: Nu/Re history files contain more rows than unsteadySampleCount.'
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,'(A)') 'Error: Nu/Re history files contain more rows than unsteadySampleCount.'
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

    open(unit=33, file='NuRe_TimeAverage_2DOpenaccUSTLBM.txt', status='replace', action='write', form='formatted')
    write(33,'(A)') '# 2D OpenACC USTLBM Nu/Re statistical-convergence window averages'
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
    coef   = lengthUnit/diffusivity ! 热流采用热扩散标度，独立于输出速度比较标度

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
    coef   = lengthUnit/diffusivity ! 热流采用热扩散标度，独立于输出速度比较标度

   

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
  coef   = lengthUnit/diffusivity ! 热流采用热扩散标度，独立于输出速度比较标度

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
  coef   = lengthUnit/diffusivity ! 热流采用热扩散标度，独立于输出速度比较标度

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
! 子程序: align_ustlbm_restart_history
! 作用: 先核验两份 Nu/Re 历史在 checkpoint 之前的完整前缀，再恢复到共同采样边界。
! 超出 checkpoint 的尾部（含断电留下的坏尾行）只在原文件逐字节备份成功后截去。
! 本过程仅在主机执行；steadyFlow 不依赖此历史，保持为空操作。
!===================================================================================================
subroutine align_ustlbm_restart_history()
  use commondata
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: iostat_end, iostat_eor, int64
  implicit none
#ifdef unsteadyFlow
  type saved_history_line
    character(len=:), allocatable :: text
  end type saved_history_line
  type checked_history
    logical :: exists=.false., repair=.false.
    type(saved_history_line), allocatable :: prefix(:)
    real(kind=8), allocatable :: sampleTime(:)
  end type checked_history
  type(checked_history) :: nuHistory, reHistory
  character(len=*), parameter :: nuName='Nu_VolAvg_2DOpenaccUSTLBM.dat'
  character(len=*), parameter :: reName='Re_VolAvg_2DOpenaccUSTLBM.dat'
  integer(kind=4) :: k
  real(kind=8) :: expectedTime, timeTolerance

  if(sampleIndexOffset.LT.0) error stop 'Negative checkpoint history sample count'
  ! 两文件全部通过前缀核验前，不创建备份、不替换任一历史文件。
  call inspect_history(nuName,nuHistory)
  call inspect_history(reName,reHistory)
  do k=1,sampleIndexOffset
    expectedTime=dble(k)*outputSnapshotInterval
    timeTolerance=64.0d0*epsilon(1.0d0)*max(1.0d0,abs(expectedTime))
    if(abs(nuHistory%sampleTime(k)-reHistory%sampleTime(k)).GT.timeTolerance) &
      error stop 'Nu/Re checkpoint history times are not paired'
  enddo

  ! 所有需要修复的已有文件均先完成备份，之后才允许改写任何原文件。
  if(nuHistory%repair.AND.nuHistory%exists) call backup_history(nuName)
  if(reHistory%repair.AND.reHistory%exists) call backup_history(reName)
  if(nuHistory%repair) call replace_with_prefix(nuName,nuHistory)
  if(reHistory%repair) call replace_with_prefix(reName,reHistory)

contains

  subroutine inspect_history(filename,history)
    character(len=*), intent(in) :: filename
    type(checked_history), intent(out) :: history
    integer(kind=4) :: historyUnit, ios, row
    real(kind=8) :: rowTime, rowValue, target, tolerance
    character(len=:), allocatable :: line

    allocate(history%prefix(sampleIndexOffset),history%sampleTime(sampleIndexOffset))
    inquire(file=filename,exist=history%exists)
    if(.not.history%exists) then
      if(sampleIndexOffset.GT.0) then
        write(*,'(A,1X,A)') 'Missing required checkpoint history:',filename
        error stop 'Checkpoint history prefix is missing'
      endif
      history%repair=.true. ! 零样本 checkpoint 可以创建空文件，供后续 append 使用。
      return
    endif
    open(newunit=historyUnit,file=filename,status='old',action='read',form='formatted',iostat=ios)
    if(ios.NE.0) error stop 'Cannot open checkpoint history for validation'
    do row=1,sampleIndexOffset
      call read_history_line(historyUnit,line,ios)
      if(ios.NE.0) then
        write(*,'(A,1X,A,1X,A,I0)') 'Incomplete checkpoint history:',filename,'row=',row
        error stop 'Checkpoint history prefix is incomplete or unreadable'
      endif
      read(line,*,iostat=ios) rowTime,rowValue
      if(ios.NE.0) then
        write(*,'(A,1X,A,1X,A,I0)') 'Invalid checkpoint history:',filename,'row=',row
        error stop 'Checkpoint history prefix contains a damaged row'
      endif
      if(.not.ieee_is_finite(rowTime).OR..not.ieee_is_finite(rowValue)) &
        error stop 'Checkpoint history prefix contains a non-finite value'
      target=dble(row)*outputSnapshotInterval
      tolerance=64.0d0*epsilon(1.0d0)*max(1.0d0,abs(target))
      if(abs(rowTime-target).GT.tolerance) then
        write(*,'(A,1X,A,1X,A,I0)') 'Off-grid checkpoint history:',filename,'row=',row
        error stop 'Checkpoint history times must equal k times the sampling interval'
      endif
      history%prefix(row)%text=line
      history%sampleTime(row)=rowTime
    enddo
    ! 前缀之后任何额外记录或坏尾行都必须去除；前缀本身已完整核验。
    call read_history_line(historyUnit,line,ios)
    history%repair=(ios.NE.iostat_end)
    close(historyUnit,iostat=ios)
    if(ios.NE.0) error stop 'Cannot close checkpoint history after validation'
  end subroutine inspect_history

  subroutine read_history_line(historyUnit,line,ios)
    integer(kind=4), intent(in) :: historyUnit
    character(len=:), allocatable, intent(out) :: line
    integer(kind=4), intent(out) :: ios
    character(len=1024) :: chunk
    integer(kind=4) :: received

    line=''
    do
      received=0
      read(historyUnit,'(A)',advance='no',size=received,iostat=ios) chunk
      if(received.GT.0) line=line//chunk(:received)
      if(ios.EQ.0) cycle
      if(ios.EQ.iostat_eor) then
        ios=0
      elseif(ios.EQ.iostat_end.AND.len(line).GT.0) then
        ios=0 ! 无换行的最后一行仍须按内容核验，不能误当作缺失前缀。
      endif
      return
    enddo
  end subroutine read_history_line

  subroutine backup_history(filename)
    character(len=*), intent(in) :: filename
    integer(kind=4) :: sourceUnit, backupUnit, ios, attempt, chunkLength
    integer(kind=int64) :: fileSize, position
    character(len=65536) :: buffer
    character(len=512) :: backupName
    logical :: alreadyExists

    open(newunit=sourceUnit,file=filename,status='old',action='read', &
      access='stream',form='unformatted',iostat=ios)
    if(ios.NE.0) error stop 'Cannot open original history for exact backup'
    inquire(unit=sourceUnit,size=fileSize,iostat=ios)
    if(ios.NE.0) error stop 'Cannot obtain original history size'
    attempt=0
    do
      if(attempt.EQ.huge(attempt)) error stop 'No unique history backup name remains'
      attempt=attempt+1
      write(backupName,'(A,A,I0,A,I0,A)') filename,'.restart-',restartItcOffset,'-',attempt,'.bak'
      inquire(file=trim(backupName),exist=alreadyExists)
      if(alreadyExists) cycle
      ! status=new 在检查与打开之间也禁止覆盖已存在的备份。
      open(newunit=backupUnit,file=trim(backupName),status='new',action='write', &
        access='stream',form='unformatted',iostat=ios)
      if(ios.EQ.0) exit
      inquire(file=trim(backupName),exist=alreadyExists)
      if(.not.alreadyExists) error stop 'Cannot create a unique history backup'
    enddo
    position=1_int64
    do while(position.LE.fileSize)
      chunkLength=int(min(int(len(buffer),int64),fileSize-position+1_int64))
      read(sourceUnit,pos=position,iostat=ios) buffer(:chunkLength)
      if(ios.NE.0) error stop 'Cannot read complete original history for backup'
      write(backupUnit,iostat=ios) buffer(:chunkLength)
      if(ios.NE.0) error stop 'Cannot write complete history backup'
      position=position+int(chunkLength,int64)
    enddo
    close(sourceUnit,iostat=ios)
    if(ios.NE.0) error stop 'Cannot close original history after backup'
    close(backupUnit,iostat=ios)
    if(ios.NE.0) error stop 'Cannot finish complete history backup'
    write(*,'(A,1X,A)') 'Saved complete pre-repair history:',trim(backupName)
  end subroutine backup_history

  subroutine replace_with_prefix(filename,history)
    character(len=*), intent(in) :: filename
    type(checked_history), intent(in) :: history
    integer(kind=4) :: historyUnit, ios, row

    open(newunit=historyUnit,file=filename,status='replace',action='write',form='formatted',iostat=ios)
    if(ios.NE.0) error stop 'Cannot restore checkpoint history prefix'
    do row=1,sampleIndexOffset
      write(historyUnit,'(A)',iostat=ios) history%prefix(row)%text
      if(ios.NE.0) error stop 'Cannot write complete checkpoint history prefix'
    enddo
    close(historyUnit,iostat=ios)
    if(ios.NE.0) error stop 'Cannot finish checkpoint history prefix restoration'
    write(*,'(A,1X,A,1X,A,I0)') 'Aligned checkpoint history:',filename,'samples=',sampleIndexOffset
  end subroutine replace_with_prefix
#endif
end subroutine align_ustlbm_restart_history
!===================================================================================================

