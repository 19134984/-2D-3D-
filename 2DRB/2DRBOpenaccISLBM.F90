!=============================================================
!!!    注释区，代码描述
!!!    二维浮力驱动自然对流 OpenACC 并行版本
!!!    LBM方法
!!!    MRT-LBE
!=============================================================


!=============================================================
!   自定义宏，一些选项的开关
#define steadyFlow
!#define unsteadyFlow

!   流动模式宏的选择，两个都开、两个都关都会报错；只有二选一才通过。
#if defined(steadyFlow) && defined(unsteadyFlow)
#error "Choose only one flow mode: steadyFlow or unsteadyFlow"
#endif
#if !defined(steadyFlow) && !defined(unsteadyFlow)
#error "Define one flow mode: steadyFlow or unsteadyFlow"
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
!#define RayleighBenardCell
!#define HorizontalWallsConstT
!#define VerticalWallsAdiabatic
!#define VerticalWallsPeriodicalT



!   温度边界(for Side Heated Cell)，包括水平边界温度不可穿透，垂直边界恒温,侧壁加热加磁场
#define SideHeatedCell
#define HorizontalWallsAdiabatic
#define VerticalWallsConstT
!#define SideHeatedHa  
!~~temperature B.C.~~

!   对流算例宏的选择
#if defined(RayleighBenardCell) && defined(SideHeatedCell)
#error "Choose only one convection case: RayleighBenardCell or SideHeatedCell"
#endif
#if !defined(RayleighBenardCell) && !defined(SideHeatedCell)
#error "Define one convection case: RayleighBenardCell or SideHeatedCell"
#endif

!算法切换
!启用 M1G 修正；注释掉则不使用 useG 相关修正
!#define EnableUseG
!启用旧温度算法
#define EnableLegacyThermalScheme

!   温度算法宏的选择
#if defined(EnableUseG) && defined(EnableLegacyThermalScheme)
#error "Choose only one thermal scheme: EnableUseG or EnableLegacyThermalScheme"
#endif
#if !defined(EnableUseG) && !defined(EnableLegacyThermalScheme)
#error "Define one thermal scheme: EnableUseG or EnableLegacyThermalScheme"
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
        real(kind=8) :: reloadDimensionlessTime=0.0d0   ! 续算前已累计的 t_ff；优先从 latest .meta 读取，meta 缺失时由代码推断
        integer(kind=4) :: restartItcOffset=0           ! 续算前已累计的格子步数；优先从 latest .meta 读取，meta 缺失时由代码推断
        logical :: reloadMetadataLoaded=.false.         ! 标记是否成功读取 reload 元数据文件
        !===============================================================================================

        !===============================================================================================
        ! 无量纲参数
        integer(kind=4), parameter :: nx=513, ny=513     !格子网格

        ! 本文件采用erf非均匀网格。
        ! erf网格拉伸强度。数值越大, 节点越向两侧物理壁面聚集。
        real(kind=8), parameter :: ISLBM_StretchA=1.5d0

        ! raw坐标中第一个内部流体节点的位置rawX(1)/rawY(1), 也就是近壁的1个lu。
        real(kind=8), parameter :: ISLBM_DxMinRaw=0.5d0*(1.0d0 + &
            erf(ISLBM_StretchA*(1.0d0/dble(nx+1)-0.5d0))/erf(0.5d0*ISLBM_StretchA))
        real(kind=8), parameter :: ISLBM_DyMinRaw=0.5d0*(1.0d0 + &
            erf(ISLBM_StretchA*(1.0d0/dble(ny+1)-0.5d0))/erf(0.5d0*ISLBM_StretchA))
#ifdef SideHeatedCell
        ! half-way壁面放在0.5*rawX(1)或0.5*rawY(1), 因此有效长度为1-ISLBM_Dx/DyMinRaw。
        ! ISLBM有效长度含多少个近壁lu: lengthUnit=(rightWall-leftWall)/ISLBM_DxMinRaw。
        real(kind=8), parameter :: lengthUnit=(1.0d0-ISLBM_DxMinRaw)/ISLBM_DxMinRaw
#else
        ! ISLBM有效长度含多少个近壁lu: lengthUnit=(topWall-bottomWall)/ISLBM_DyMinRaw。
        real(kind=8), parameter :: lengthUnit=(1.0d0-ISLBM_DyMinRaw)/ISLBM_DyMinRaw
#endif
        real(kind=8), parameter :: pi = acos(-1.0d0)

        real(kind=8), parameter :: Rayleigh=1.0d8
        real(kind=8), parameter :: Prandtl=0.71d0
        real(kind=8), parameter :: Mach=0.1d0
        real(kind=8), parameter :: Thot=0.5d0, Tcold=-0.5d0
        real(kind=8), parameter :: Tref=0.5d0*(Thot+Tcold)
        real(kind=8), parameter :: tauf=0.5d0+Mach*lengthUnit*dsqrt(3.0d0*Prandtl/Rayleigh) 
        real(kind=8), parameter :: viscosity=(tauf-0.5d0)/3.0d0
        real(kind=8), parameter :: diffusivity=viscosity/Prandtl
        

        ! velocityScaleCompare is used only in velocity-related post-processing to convert lattice velocity
        ! to the nondimensional velocity scale adopted by the reference paper being compared.
        real(kind=8), parameter :: velocityScaleCompare=lengthUnit/diffusivity
        ! 默认采用热扩散标度 UL/kappa；若要按自由落体标度比较，可改为 1.0d0/velocityUnit
        
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1


#ifdef  SideHeatedHa
        real(kind=8), parameter :: Ha=20.0d0                           !磁场强度
        real(kind=8), parameter :: phi=(0.0d0)*(pi/180.0d0)            !磁场角度，以水平向右为0，修改0.0d0即可
        real(kind=8), parameter :: B2sigemarho=(Ha**2*viscosity)/(lengthUnit*lengthUnit)  !动量方程上的源项系数
#endif

        ! 高阶矩参数修正
        real(kind=8), parameter :: paraA=20.0d0*dsqrt(3.0d0)*diffusivity-4.0d0



        ! 浮力项参数
        real(kind=8), parameter :: gBeta1=Rayleigh*viscosity*diffusivity/lengthUnit
        real(kind=8), parameter :: gBeta=gBeta1/lengthUnit/lengthUnit             !gbetaΔT
        
        real(kind=8), parameter :: timeUnit=dsqrt(lengthUnit/gBeta)      !无量纲时间
        real(kind=8), parameter :: velocityUnit=dsqrt(gBeta*lengthUnit)  !无量纲速度
    
        real(kind=8), parameter :: Snu=1.0d0/tauf, Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)  !动量的多松弛系数

#ifdef EnableLegacyThermalScheme
        real(kind=8), parameter :: Qk=3.0d0-dsqrt(3.0d0), Qnu=4.0d0*dsqrt(3.0d0)-6.0d0             !旧温度算法的多松弛系数
        real(kind=8), parameter :: thermalGeqCoeff=10.0d0/(4.0d0+paraA)
#else
        real(kind=8), parameter :: taug = 0.5d0 + (tauf - 0.5d0)/Prandtl
        real(kind=8), parameter :: Qnu = 1.0d0, Qk = 1.0d0/taug
        real(kind=8), parameter :: thermalGeqCoeff=3.0d0
#endif
        !===============================================================================================          
        
        !===============================================================================================
        ! 输出/备份相关设置（以自由落体时间 t_ff 为单位）
        real(kind=8), parameter :: epsU=1.0d-7, epsT=1.0d-7    ! 稳态收敛阈值   

#ifdef steadyFlow
        real(kind=8), parameter :: outputSnapshotInterval=10.0d0   ! 快照和 Nu/Re 时间序列采样间隔（单位：t_ff）
        real(kind=8), parameter :: reloadFileInterval=100.0d0  ! f/g 重启文件输出间隔（单位：t_ff）
        real(kind=8), parameter :: outputPltFileInterval=100.0d0  ! Tecplot 文件输出间隔（单位：t_ff）
        integer(kind=4), parameter :: dimensionlessTimeMax=int(12000.0d0/outputSnapshotInterval)
        integer(kind=4), parameter :: outputSnapshotFile=1   ! 是否输出后处理快照文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputPltFile=1   ! 是否输出 plt 文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputReloadFile=1 ! 是否周期输出 f/g 重启文件：0=不输出，1=输出
        integer(kind=4), parameter :: itc_max=20000000  ! 稳态：最大格子步，实际仍由 errorU/errorT 提前停止
#endif

#ifdef unsteadyFlow
        real(kind=8), parameter :: outputSnapshotInterval=0.5d0   ! uvTrho 快照和 Nu/Re 时间序列采样间隔（单位：t_ff）
        real(kind=8), parameter :: reloadFileInterval=100.0d0  ! f/g 重启文件输出间隔（单位：t_ff）
        real(kind=8), parameter :: outputPltFileInterval=100.0d0  ! Tecplot 文件周期输出间隔（单位：t_ff）
        real(kind=8), parameter :: unsteadyRunDuration=1000.0d0  ! 非稳态总目标时长，续算时只补足到该 t_ff
        ! 以下三个参数只控制非稳态结束后的 Nu/Re 统计平均窗口，不改变推进时长或采样频率。
        ! 时间以整个非稳态算例的绝对 t_ff 计，续算统计会包含旧文件中已有的历史数据。
        real(kind=8), parameter :: unsteadyAverageStartTf=0.5d0*unsteadyRunDuration  ! 平均窗口起点
        real(kind=8), parameter :: unsteadyAverageEndTf=unsteadyRunDuration          ! 平均窗口终点
        real(kind=8), parameter :: unsteadyAverageMidTf=0.5d0*(unsteadyAverageStartTf+unsteadyAverageEndTf) ! 前/后半分界
        integer(kind=4), parameter :: unsteadySampleCount=max(1, &
          int(unsteadyRunDuration/outputSnapshotInterval+0.5d0))
        ! 计数器，输出多少次快照
        integer(kind=4), parameter :: dimensionlessTimeMax=unsteadySampleCount
        integer(kind=4), parameter :: outputSnapshotFile=1   ! 是否输出后处理快照文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputPltFile=1   ! 是否输出 plt 文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputReloadFile=1 ! 是否周期输出 f/g 重启文件：0=不输出，1=输出
        integer(kind=4), parameter :: itc_max=max(1, int(unsteadyRunDuration*timeUnit+0.5d0)) ! 非稳态：由总 t_ff 自动换算格子步
#endif

        integer(kind=4) :: snapshotFileNum, pltFileNum  ! 快照/plt 输出文件的计数器
        ! 每次调用对应输出子程序时递增（用于文件名编号）

        integer(kind=4) :: dimensionlessTime
        integer(kind=4) :: outputSnapshotIntervalItc
        integer(kind=4) :: reloadFileIntervalItc, outputPltFileIntervalItc
        ! 统计/输出时间点编号（与 outputSnapshotInterval 对应）：
        ! 每调用一次 calNuRe() 就 dimensionlessTime = dimensionlessTime + 1
        ! 用于索引 NuVolAvg/ReVolAvg 数组，并用于输出的时间轴：t = reloadDimensionlessTime + dimensionlessTime*outputSnapshotInterval（单位：t_ff）

        real(kind=8) :: NuVolAvg(0:dimensionlessTimeMax), ReVolAvg(0:dimensionlessTimeMax)
        ! 体平均 Nu 和 Re 的时间序列缓存
        ! 只有在启用并调用 calNuRe() 的情况下这些数组才会被真正填充
  
        character(len=100) :: snapshotFilePrefix="buoyancyCavity2DOpenaccSnapshot"
        ! 快照输出文件前缀（实际文件名形如：<snapshotFilePrefix>-<编号>.bin）

        character(len=100) :: pltFolderPrefix="buoyancyCavity2DOpenaccTecplot"
        ! plt 输出文件前缀（实际文件名形如：<pltFolderPrefix>-<编号>.plt）

        character(len=100) :: reloadFilePrefix="reloadFile2DOpenacc"
        ! 重启读取文件的前缀；latest meta 模式实际读取 meta 中记录的 <reloadFilePrefix>-<编号>.bin
        
        character(len=100) :: settingsFile="SimulationSettings2DOpenacc.txt"
        !===============================================================================================

        !===============================================================================================
        !计算中需要的相关参数
        real(kind=8) :: errorU, errorT
        
        real(kind=8) :: xp(0:nx+1), yp(0:ny+1)      !无量纲的坐标数组，包括边界
        real(kind=8) :: quadWidthX(1:nx), quadWidthY(1:ny), quadSumX, quadSumY, quadSumArea
        ! 归一化坐标中的1个lattice unit; 迁移时用xp(i)-ex(alpha)*ISLBM_LatticeUnit找上游点。
        real(kind=8), parameter :: ISLBM_LatticeUnit=1.0d0/lengthUnit
        integer(kind=4) :: streamInterpIndexX(0:8,1:nx,3), streamInterpIndexY(0:8,1:ny,3)
        real(kind=8) :: streamInterpWeightX(0:8,1:nx,3), streamInterpWeightY(0:8,1:ny,3)
        logical :: streamInterpValidX(0:8,1:nx), streamInterpValidY(0:8,1:ny)
        real(kind=8), allocatable :: u(:,:), v(:,:), T(:,:), rho(:,:)

#ifdef steadyFlow
        real(kind=8), allocatable :: up(:,:), vp(:,:), Tp(:,:)   !存储之前的数据，用来算收敛判据
#endif
        real(kind=8), allocatable :: f(:,:,:), f_post(:,:,:)
        real(kind=8), allocatable :: g(:,:,:), g_post(:,:,:)
        real(kind=8), allocatable :: Fx(:,:), Fy(:,:)

        real(kind=8), allocatable :: Bx_prev(:,:), By_prev(:,:)

        integer(kind=4) :: itc
#ifdef EnableUseG
        logical, parameter :: useG = .true.            !M1G 开关
#else
        logical, parameter :: useG = .false.           !M1G 开关
#endif

#ifdef EnableLegacyThermalScheme
        logical, parameter :: useLegacyThermalScheme = .true.            !旧算法的逻辑变量
#else
        logical, parameter :: useLegacyThermalScheme = .false.           
#endif
#ifdef steadyFlow
        real(kind=8) :: Nu_global, Nu_hot, Nu_cold, Nu_middle    !平均Nu，全场，侧壁以及中线
        real(kind=8) :: Nu_hot_max, Nu_hot_min, Nu_hot_max_position, Nu_hot_min_position    !左侧壁面的最大最小Nu，以及对应的位置
#endif
        
        
        !格子离散速度和权重
        integer(kind=4) :: ex(0:8), ey(0:8)
        data ex/0, 1, 0, -1,  0, 1, -1, -1,  1/
        data ey/0, 0, 1,  0, -1, 1,  1, -1, -1/
        real(kind=8) :: omega(0:8), omegaT(0:4)
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
    do while( ((errorU.GT.epsU).OR.(errorT.GT.epsT)).AND.(itc.LE.itc_max) )
#endif
#ifdef unsteadyFlow
    do while( itc.LT.unsteadyItcRemaining )   !非稳态：续算时只推进到 unsteadyRunDuration 对应的总格子步
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
        ! 周期输出按累计格子步判断；否则从 1050tf 续算会在 1150tf 才输出，
        ! 而不是接回不断电运行应有的 1100tf、1200tf、...
        if(MOD(restartItcOffset+itc,2000).EQ.0) then
            call check()
            call output_steady_monitor()
        endif
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
        do while( (reloadDimensionlessTime + real(dimensionlessTime,kind=8)*outputSnapshotInterval) &
            .LT. unsteadyRunDuration )   !判断本次还要不要继续采样
            ! 每个目标采样时刻都按绝对 t_ff 换算到本次运行段的 itc，续算时不会重复旧样本。
            nextSampleAbsItc = max(1, int((reloadDimensionlessTime + &
                real(dimensionlessTime+1,kind=8)*outputSnapshotInterval)*timeUnit+0.5d0)) !算下一个采样点对应的绝对格子步
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
    write(00,*) "MLUPS = ", real( dble(nx)*dble(ny)*dble(itc)/(timeEnd-timeStart)/1.0d6,kind=8 )   !百万格点更新/秒
    write(00,*) "Time (ACC) = ", real(timeEnd2-timeStart2,kind=8), "s"                           !墙钟时间
    write(00,*) "MLUPS (ACC) = ", real( dble(nx)*dble(ny)*dble(itc)/(timeEnd2-timeStart2)/1.0d6,kind=8 )   !百万格点更新/秒
#ifdef steadyFlow
    write(00,'(a,1x,ES24.16E3)') "Nu_global =", real(Nu_global,kind=8)
    write(00,'(a,1x,ES24.16E3)') "Nu_hot    =", real(Nu_hot,kind=8)
    write(00,'(a,1x,ES24.16E3)') "Nu_cold   =", real(Nu_cold,kind=8)
#endif
    write(00,*) "useG =", useG



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
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===========================================================================================================================
  subroutine initial()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: un(0:8)
    real(kind=8) :: us2, Bx, By
    real(kind=8) :: xLen, yLen, rbInitPerturbAmp
    character(len=100) :: reloadFileName
    

    itc = 0
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
        write(01,*) "records: U_nd, V_nd, T, rho on ISLBM nonuniform nodes"
        write(01,*) "coordinates are xp(1:nx), yp(1:ny) in ", trim(snapshotFilePrefix)//"-mesh.dat"
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
    
#ifdef EnableLegacyThermalScheme
    if( (paraA.GE.1.0d0).OR.(paraA.LE.-4.0d0) ) then                           !只有在[-4,1]才可以，要不然预警退出
        write(00,*) "----------------------------------"
        write(00,*) "paraA=", paraA
        write(00,*) "Error: condition not meet for the legacy thermal algorithm"
        write(00,*) "Ref: Luo2013, CMA"
        write(00,*) "Please try to reduce Mach number"
        write(00,*) "----------------------------------"
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
    write(00,*) 'tauf=',real(tauf,kind=8)
#ifdef EnableLegacyThermalScheme
    write(00,*) "thermalScheme = legacy D2Q5 (Qk/Qnu)"
    write(00,*) 'Qk=',real(Qk,kind=8), '; Qnu=',real(Qnu,kind=8), '; paraA=',real(paraA,kind=8)
#else
    write(00,*) "thermalScheme = current D2Q5 (s_j/s_e/s_q)"
    write(00,*) 'taug=',real(taug,kind=8)
#endif
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
        write(00,*) "Restart offsets will be read from reload metadata when available."
    endif
    write(00,*) "itc_max =",itc_max
    write(00,*) "default epsU =", real(epsU,kind=8),"; epsT =", real(epsT,kind=8)
    write(00,*) "useG =", useG
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
    write(00,*) "OpenACC GPU ISLBM version; final wall-fit and psi/vorticity diagnostics stay on host"
    !-----------------------------------------------------------------------------------------------



    !-----------------------------------------------------------------------------------------------
    ! ISLBM节点坐标与迁移模板。坐标按 half-way 物理壁面归一化到 [0,1]，
    ! 第一个流体节点位于 0.5/L0，最后一个流体节点位于 1-0.5/L0。
    call build_islbm_mesh()
    call build_islbm_quadrature()
    call build_islbm_streaming_stencils()
    if(outputSnapshotFile.EQ.1) then
        call output_SnapshotMeshFile()
        write(00,*) "Snapshot mesh coordinates stored in ", trim(snapshotFilePrefix)//"-mesh.dat"
    endif
    write(00,*) "ISLBM mesh = erf; stretchA =", real(ISLBM_StretchA,kind=8)
    write(00,*) "ISLBM effective lengthUnit L0 =", real(lengthUnit,kind=8)
    write(00,*) "ISLBM lattice unit in normalized coordinates =", real(ISLBM_LatticeUnit,kind=8)
    write(00,*) "ISLBM quadrature sums =", real(quadSumX,kind=8), real(quadSumY,kind=8), real(quadSumArea,kind=8)

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

    allocate (Bx_prev(nx,ny), By_prev(nx,ny)) 

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
    
#ifdef EnableLegacyThermalScheme
    omegaT(0) = (1.0d0-paraA)/5.0d0
    do alpha=1,4
        omegaT(alpha) = (paraA+4.0d0)/20.0d0
    enddo
#else
    omegaT(0) = 1.0d0/3.0d0
    do alpha=1,4
        omegaT(alpha) = 1.0d0/6.0d0
    enddo
#endif

    if(loadInitField.EQ.0) then                    !在不加载文件的情况下，都是零场为初值
    
        u = 0.0d0
        v = 0.0d0
        T = 0.0d0
        Bx= 0.0d0
        By= 0.0d0
        Bx_prev= 0.0d0
        By_prev= 0.0d0
        
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
    if (Rayleigh.LT.1.0d4) then
        xLen = xp(nx+1)
        yLen = yp(ny+1)
        rbInitPerturbAmp = 1.0d-3*(Thot-Tcold)
        do i = 1, nx
            do j = 1, ny
                T(i,j) = T(i,j) + rbInitPerturbAmp * dsin(2.0d0*pi*xp(i)/xLen) * dsin(pi*yp(j)/yLen)
            enddo
        enddo
        write(00,'(a,1x,es12.4)') "RB initial T perturbation amplitude =", rbInitPerturbAmp
    else
        write(00,*) "RB initial T perturbation skipped because Rayleigh > 1.0d4"
    endif
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
        g = 0.0d0
        do j = 1,ny
            do i = 1,nx
                us2 = u(i,j)*u(i,j)+v(i,j)*v(i,j)
                do alpha = 0, 8
                    un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                    f(i,j,alpha) = rho(i,j)*omega(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)  !D2Q9标准feq
                enddo
                do alpha = 0, 4
                    un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha) 
                    g(i,j,alpha) = omegaT(alpha)*T(i,j)*(1.0d0+thermalGeqCoeff*un(alpha))
                enddo
            enddo
        enddo
#ifdef EnableUseG
        do j = 1,ny
            do i = 1,nx
              Bx = u(i,j) * T(i,j)
              By = v(i,j) * T(i,j)
              Bx_prev(i,j) = Bx
              By_prev(i,j) = By
            enddo
        enddo
#endif



    elseif(loadInitField.EQ.1) then                               !续算分支：从旧的严格重启文件恢复 f/g 和输出计数
    ! 正常断电续算时，先读取 <reloadFilePrefix>-latest.meta；
    ! meta 会告诉代码实际要读哪个 <reloadFilePrefix>-*.bin，以及旧算例已经累计到哪里。
    !这里可以让 reloadDimensionlessTime 等于0，反正最后都可以通过meta文件读取到新的 reloadDimensionlessTime
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
            ! Strict restart files store f and g; EnableUseG also stores the previous heat-flux history.
            write(00,*) "Reloading f, g and optional UseG history from file"
            read(01) (((f(i,j,alpha), i=1,nx), j=1,ny), alpha=0,8)      !先 i，再 j，再 alpha
            read(01) (((g(i,j,alpha), i=1,nx), j=1,ny), alpha=0,4)
#ifdef EnableUseG
            read(01) ((Bx_prev(i,j), i=1,nx), j=1,ny)
            read(01) ((By_prev(i,j), i=1,nx), j=1,ny)
#endif
        close(01)
        call reconstruct_macro_from_fg()
        write(00,*) "Raw data is loaded from the file: ", trim(reloadFilePrefix), "-", trim(reloadFileName), ".bin"
        write(00,*) "Restart offset itc =", restartItcOffset
        write(00,*) "Restart offset time_tf =", real(reloadDimensionlessTime,kind=8)
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
        reloadDimensionlessTime = 0.0d0
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
    dimensionlessTime = 0
    !新算例：清零，开始记录新历史。
    !续算：也清零，但不是丢旧历史；旧历史在 .dat 文件里，新数组只记录本次续算段。
    NuVolAvg = 0.0d0
    ReVolAvg = 0.0d0
    
    return
  end subroutine initial
!===================================================================================================
! initial 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
!===================================================================================================
!初始化结束
!===================================================================================================



!===================================================================================================
!===================================================================================================
! 子程序: build_islbm_mesh
! 作用: 执行本子程序对应的初始化、迁移、碰撞、边界、通信或后处理步骤。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine build_islbm_mesh()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: rawX(0:nx+1), rawY(0:ny+1)
    real(kind=8) :: erfNorm, leftWall, rightWall, bottomWall, topWall, lengthX, lengthY

    ! 第一步: 生成原始erf拉伸坐标rawX/rawY。此时坐标还没有按half-way物理壁面修正。
    erfNorm = erf(0.5d0*ISLBM_StretchA)
    do i = 0, nx+1
        rawX(i) = 0.5d0*(1.0d0 + erf(ISLBM_StretchA*(dble(i)/dble(nx+1)-0.5d0))/erfNorm)
    enddo
    do j = 0, ny+1
        rawY(j) = 0.5d0*(1.0d0 + erf(ISLBM_StretchA*(dble(j)/dble(ny+1)-0.5d0))/erfNorm)
    enddo

    ! 第二步: 采用half-way壁面。物理壁面位于参考端点与第一个内部流体节点之间。
    leftWall = 0.5d0*rawX(1)
    rightWall = 1.0d0 - 0.5d0*rawX(1)
    bottomWall = 0.5d0*rawY(1)
    topWall = 1.0d0 - 0.5d0*rawY(1)

    lengthX = rightWall - leftWall
    lengthY = topWall - bottomWall
    xp(0) = 0.0d0
    xp(nx+1) = 1.0d0
    do i = 1, nx
        xp(i) = (rawX(i)-leftWall)/lengthX
    enddo
    yp(0) = 0.0d0
    yp(ny+1) = 1.0d0
    do j = 1, ny
        yp(j) = (rawY(j)-bottomWall)/lengthY
    enddo

    return
  end subroutine build_islbm_mesh
!===================================================================================================
! build_islbm_mesh 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: build_islbm_quadrature
! 作用: 执行本子程序对应的初始化、迁移、碰撞、边界、通信或后处理步骤。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine build_islbm_quadrature()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: leftGhostX, rightGhostX, bottomGhostY, topGhostY

    ! midpoint-rule在边界控制体上的写法:
    ! 每个流体节点代表一个控制体, 控制体边界取相邻节点的中点。
    ! 对边界流体节点, 控制体的一侧边界是物理壁面0/1, 另一侧是第1、2个流体节点的中点。
    ! 这等价于在注释里引入虚拟点 x_0_virtual=-xp(1), x_{nx+1}_virtual=1+(1-xp(nx)):
    !   quadWidthX(1)=0.5*(xp(2)-x_0_virtual) = (xp(1)-0) + 0.5*(xp(2)-xp(1))
    ! 内部节点仍为0.5*(xp(i+1)-xp(i-1))。这样quadSumX/quadSumY对应完整物理长度1。
    leftGhostX = -xp(1)
    rightGhostX = 1.0d0 + (1.0d0 - xp(nx))
    bottomGhostY = -yp(1)
    topGhostY = 1.0d0 + (1.0d0 - yp(ny))

    quadWidthX(1) = 0.5d0*(xp(2)-leftGhostX)
    do i = 2, nx-1
        quadWidthX(i) = 0.5d0*(xp(i+1)-xp(i-1))
    enddo
    quadWidthX(nx) = 0.5d0*(rightGhostX-xp(nx-1))

    quadWidthY(1) = 0.5d0*(yp(2)-bottomGhostY)
    do j = 2, ny-1
        quadWidthY(j) = 0.5d0*(yp(j+1)-yp(j-1))
    enddo
    quadWidthY(ny) = 0.5d0*(topGhostY-yp(ny-1))

    quadSumX = 0.0d0
    do i = 1, nx
        quadSumX = quadSumX + quadWidthX(i)
    enddo
    quadSumY = 0.0d0
    do j = 1, ny
        quadSumY = quadSumY + quadWidthY(j)
    enddo
    quadSumArea = quadSumX*quadSumY

    return
  end subroutine build_islbm_quadrature
!===================================================================================================
! build_islbm_quadrature 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: build_islbm_streaming_stencils
! 作用: 执行本子程序对应的初始化、迁移、碰撞、边界、通信或后处理步骤。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine build_islbm_streaming_stencils()
    use commondata
    implicit none
    integer(kind=4) :: alpha, i, j
    integer(kind=4) :: idx(3)
    real(kind=8) :: w(3), target
    logical :: ok

    streamInterpIndexX = 1
    streamInterpIndexY = 1
    streamInterpWeightX = 0.0d0
    streamInterpWeightY = 0.0d0
    streamInterpValidX = .false.
    streamInterpValidY = .false.

    do alpha = 0, 8
        do i = 1, nx
            ! pull形式: f_alpha(x_i,t+1) 需要上一时刻碰撞后沿-alpha方向的出发点。
            ! 对均匀格子 target 会正好落在相邻格点; 非均匀格子一般是off-lattice点,
            ! 所以用迁移后的 f_post 在当前节点附近的三点模板上插值。
            target = xp(i) - dble(ex(alpha))*ISLBM_LatticeUnit
            call build_streaming_stencil_1d(nx, xp(1:nx), i, target, idx, w, ok)
            streamInterpValidX(alpha,i) = ok
            streamInterpIndexX(alpha,i,:) = idx
            streamInterpWeightX(alpha,i,:) = w
        enddo
        do j = 1, ny
            target = yp(j) - dble(ey(alpha))*ISLBM_LatticeUnit
            call build_streaming_stencil_1d(ny, yp(1:ny), j, target, idx, w, ok)
            streamInterpValidY(alpha,j) = ok
            streamInterpIndexY(alpha,j,:) = idx
            streamInterpWeightY(alpha,j,:) = w
        enddo
    enddo

    return
  end subroutine build_islbm_streaming_stencils
!===================================================================================================
! build_islbm_streaming_stencils 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: build_streaming_stencil_1d
! 作用: 执行本子程序对应的初始化、迁移、碰撞、边界、通信或后处理步骤。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine build_streaming_stencil_1d(n, xnodes, nodeIndex, target, idx, w, ok)
    implicit none
    integer(kind=4), intent(in) :: n, nodeIndex
    real(kind=8), intent(in) :: xnodes(n), target
    integer(kind=4), intent(out) :: idx(3)
    real(kind=8), intent(out) :: w(3)
    logical, intent(out) :: ok
    real(kind=8) :: xloc(3)
    real(kind=8), parameter :: tol = 1.0d-12

    idx = (/1, 1, 1/)
    w = 0.0d0
    ok = .false.
    if(n.LT.3) return
    if((target.LT.xnodes(1)-tol).OR.(target.GT.xnodes(n)+tol)) return

    ! ISLBM streaming的插值模板以当前到达节点为中心选择。
    ! 靠近物理边界时, 指向壁面的分布函数会先由边界条件处理并直接return;
    ! 仍需要插值的方向只会使用一侧三点模板(/1,2,3/)或(/n-2,n-1,n/)。
    if(nodeIndex.LE.1) then
        idx = (/1, 2, 3/)
    elseif(nodeIndex.GE.n) then
        idx = (/n-2, n-1, n/)
    else
        idx = (/nodeIndex-1, nodeIndex, nodeIndex+1/)
    endif

    xloc(1) = xnodes(idx(1))
    xloc(2) = xnodes(idx(2))
    xloc(3) = xnodes(idx(3))
    call build_lagrange_weights_3(xloc, target, w)
    ok = .true.

    return
  end subroutine build_streaming_stencil_1d
!===================================================================================================
! build_streaming_stencil_1d 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: build_lagrange_stencil_1d
! 作用: 执行本子程序对应的初始化、迁移、碰撞、边界、通信或后处理步骤。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine build_lagrange_stencil_1d(n, xnodes, target, idx, w, ok)
    implicit none
    integer(kind=4), intent(in) :: n
    real(kind=8), intent(in) :: xnodes(n), target
    integer(kind=4), intent(out) :: idx(3)
    real(kind=8), intent(out) :: w(3)
    logical, intent(out) :: ok
    integer(kind=4) :: mid
    real(kind=8) :: xloc(3)
    real(kind=8), parameter :: tol = 1.0d-12

    idx = (/1, 1, 1/)
    w = 0.0d0
    ok = .false.
    if(n.LT.3) return
    if((target.LT.xnodes(1)-tol).OR.(target.GT.xnodes(n)+tol)) return

    if(target.LE.xnodes(2)) then
        idx = (/1, 2, 3/)
    elseif(target.GE.xnodes(n-1)) then
        idx = (/n-2, n-1, n/)
    else
        mid = 2
        do while((mid.LT.n-1).AND.(xnodes(mid+1).LT.target))
            mid = mid + 1
        enddo
        idx = (/mid-1, mid, mid+1/)
    endif

    xloc(1) = xnodes(idx(1))
    xloc(2) = xnodes(idx(2))
    xloc(3) = xnodes(idx(3))
    call build_lagrange_weights_3(xloc, target, w)
    ok = .true.

    return
  end subroutine build_lagrange_stencil_1d
!===================================================================================================
! build_lagrange_stencil_1d 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: build_lagrange_weights_3
! 作用: 执行本子程序对应的初始化、迁移、碰撞、边界、通信或后处理步骤。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine build_lagrange_weights_3(xnode, x0, w)
    implicit none
    real(kind=8), intent(in) :: xnode(3), x0
    real(kind=8), intent(out) :: w(3)

    w(1) = ((x0-xnode(2))*(x0-xnode(3)))/((xnode(1)-xnode(2))*(xnode(1)-xnode(3)))
    w(2) = ((x0-xnode(1))*(x0-xnode(3)))/((xnode(2)-xnode(1))*(xnode(2)-xnode(3)))
    w(3) = ((x0-xnode(1))*(x0-xnode(2)))/((xnode(3)-xnode(1))*(xnode(3)-xnode(2)))

    return
  end subroutine build_lagrange_weights_3
!===================================================================================================
! build_lagrange_weights_3 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: lagrange_derivative_3
! 作用: 执行本子程序对应的初始化、迁移、碰撞、边界、通信或后处理步骤。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine lagrange_derivative_3(xnode, fnode, x0, derivativeValue)
    implicit none
    real(kind=8), intent(in) :: xnode(3), fnode(3), x0
    real(kind=8), intent(out) :: derivativeValue

    derivativeValue = &
        fnode(1)*(2.0d0*x0-xnode(2)-xnode(3))/((xnode(1)-xnode(2))*(xnode(1)-xnode(3))) + &
        fnode(2)*(2.0d0*x0-xnode(1)-xnode(3))/((xnode(2)-xnode(1))*(xnode(2)-xnode(3))) + &
        fnode(3)*(2.0d0*x0-xnode(1)-xnode(2))/((xnode(3)-xnode(1))*(xnode(3)-xnode(2)))

    return
  end subroutine lagrange_derivative_3
!===================================================================================================
! lagrange_derivative_3 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================

!===================================================================================================
! 子程序: integrate_lagrange_3_segment
! 作用: 执行本子程序对应的初始化、迁移、碰撞、边界、通信或后处理步骤。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine integrate_lagrange_3_segment(xnode, fnode, xLeft, xRight, integralValue)
    implicit none
    real(kind=8), intent(in) :: xnode(3), fnode(3), xLeft, xRight
    real(kind=8), intent(out) :: integralValue
    integer(kind=4) :: k, m, n
    real(kind=8) :: denom, nodeSum, nodeProduct, basisIntegral

    integralValue = 0.0d0
    do k = 1, 3
      if(k.EQ.1) then
        m = 2
        n = 3
      elseif(k.EQ.2) then
        m = 1
        n = 3
      else
        m = 1
        n = 2
      endif
      denom = (xnode(k)-xnode(m))*(xnode(k)-xnode(n))
      nodeSum = xnode(m) + xnode(n)
      nodeProduct = xnode(m)*xnode(n)
      basisIntegral = ((xRight**3-xLeft**3)/3.0d0 - &
          0.5d0*nodeSum*(xRight**2-xLeft**2) + &
          nodeProduct*(xRight-xLeft))/denom
      integralValue = integralValue + fnode(k)*basisIntegral
    enddo

    return
  end subroutine integrate_lagrange_3_segment
!===================================================================================================
! integrate_lagrange_3_segment 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
!===================================================================================================


!===================================================================================================
! 子程序: enter_data_2d_openacc
! 作用: 在主时间推进前把主要数组和常量映射到 OpenACC 设备端。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine enter_data_2d_openacc()
    use openacc
    use commondata
    implicit none

    !$acc enter data copyin(xp,yp,quadWidthX,quadWidthY,quadSumX,quadSumY,quadSumArea,ex,ey,omega,omegaT)
    !$acc enter data copyin(streamInterpIndexX,streamInterpIndexY,streamInterpWeightX, &
    !$acc& streamInterpWeightY,streamInterpValidX,streamInterpValidY)
    !$acc enter data copyin(u,v,T,rho,f,g,Fx,Fy,Bx_prev,By_prev)
    !$acc enter data create(f_post,g_post)
#ifdef steadyFlow
    !$acc enter data copyin(up,vp,Tp)
#endif
  end subroutine enter_data_2d_openacc
!===================================================================================================
! enter_data_2d_openacc 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
!===================================================================================================

!===================================================================================================
! 子程序: update_host_snapshot_2d_openacc
! 作用: 在主机端输出 uvTrho 快照或做 CPU 后处理前，同步宏观场。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine update_host_snapshot_2d_openacc()
    use commondata
    implicit none

    !$acc wait(1)
    !$acc update self(u,v,T,rho)
  end subroutine update_host_snapshot_2d_openacc
!===================================================================================================
! update_host_snapshot_2d_openacc 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
!===================================================================================================

!===================================================================================================
! 子程序: update_host_tecplot_2d_openacc
! 作用: Tecplot 主场输出只需要 u、v、T。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine update_host_tecplot_2d_openacc()
    use commondata
    implicit none

    !$acc wait(1)
    !$acc update self(u,v,T)
  end subroutine update_host_tecplot_2d_openacc
!===================================================================================================
! update_host_tecplot_2d_openacc 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
!===================================================================================================

!===================================================================================================
! 子程序: update_host_reload_2d_openacc
! 作用: 严格重启文件需要 f/g；EnableUseG 还需要保存历史通量。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine update_host_reload_2d_openacc()
    use commondata
    implicit none

    !$acc wait(1)
    !$acc update self(f,g)
#ifdef EnableUseG
    !$acc update self(Bx_prev,By_prev)
#endif
  end subroutine update_host_reload_2d_openacc
!===================================================================================================
! update_host_reload_2d_openacc 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
!===================================================================================================

!===================================================================================================
! 子程序: exit_data_2d_openacc
! 作用: 在计算结束后释放设备端常驻数据。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine exit_data_2d_openacc()
    use commondata
    implicit none

#ifdef steadyFlow
    !$acc exit data delete(up,vp,Tp)
#endif
    !$acc exit data delete(f_post,g_post,u,v,T,rho,f,g,Fx,Fy,Bx_prev,By_prev)
    !$acc exit data delete(streamInterpIndexX,streamInterpIndexY,streamInterpWeightX, &
    !$acc& streamInterpWeightY,streamInterpValidX,streamInterpValidY)
    !$acc exit data delete(xp,yp,quadWidthX,quadWidthY,quadSumX,quadSumY,quadSumArea,ex,ey,omega,omegaT)
  end subroutine exit_data_2d_openacc
!===================================================================================================
! exit_data_2d_openacc 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
!===================================================================================================



!===================================================================================================
! 子程序: collision
! 作用: 完成流场分布函数 f 的碰撞更新，并处理体力项离散修正。
! 用途: 在主程序时间推进循环中调用，位于 streaming 之前。
!===================================================================================================
  subroutine collision()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: m(0:8), m_post(0:8), meq(0:8)
    real(kind=8) :: s(0:8)
    real(kind=8) :: fSource(0:8)

    !$acc parallel loop gang vector collapse(2) default(none) present(f,f_post,rho,u,v,Fx,Fy,T) async(1) &
    !$acc& private(alpha,s,m,m_post,meq,fSource)
    do j = 1, ny
        do i = 1, nx

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
          meq(7) = rho(i,j)*( u(i,j)*u(i,j)-v(i,j)*v(i,j) )
          meq(8) = rho(i,j)*( u(i,j)*v(i,j) ) 

          s(0) = 0.0d0      !!s_{\rho}
          s(1) = Snu !!s_{e}
          s(2) = Snu !!s_{\epsilon}
          s(3) = 0.0d0      !!s_{j} 
          s(4) = Sq !!s_{q}
          s(5) = 0.0d0      !!s_{j}
          s(6) = Sq       !!s_{q}
          s(7) = Snu !!s_{\nu}
          s(8) = Snu       !!s_{\nu}

          Fx(i,j) = 0.0d0
          Fy(i,j) = rho(i,j)*gBeta*(T(i,j)-Tref)        !动量方程上的源项，即浮力项


#ifdef    SideHeatedHa
          Fx(i,j) = 0.0d0+B2sigemarho*(v(i,j)*sin(phi)*cos(phi)-u(i,j)*sin(phi)*sin(phi))
          Fy(i,j) = rho(i,j)*gBeta*(T(i,j)-Tref)+ rho(i,j)*B2sigemarho*(u(i,j)*sin(phi)*cos(phi)&
          -v(i,j)*cos(phi)*cos(phi))                    !动量方程上的源项，即浮力项加磁场
#endif


          fSource(0) = 0.0d0                                                       !将源项F对应的贡献投影到各个矩中，并做半步修正
          fSource(1) = (6.0d0-3.0d0*s(1))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
          fSource(2) = -(6.0d0-3.0d0*s(2))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
          fSource(3) = (1.0d0-0.5d0*s(3))*Fx(i,j)
          fSource(4) = -(1.0d0-0.5d0*s(4))*Fx(i,j)
          fSource(5) = (1.0d0-0.5d0*s(5))*Fy(i,j)
          fSource(6) = -(1.0d0-0.5d0*s(6))*Fy(i,j)
          fSource(7) = (2.0d0-s(7))*(u(i,j)*Fx(i,j)-v(i,j)*Fy(i,j))
          fSource(8) = (1.0d0-0.5d0*s(8))*(u(i,j)*Fy(i,j)+v(i,j)*Fx(i,j))     !这边是乘以M变到矩空间，然后再乘以1-1/2S修正

          do alpha = 0, 8
            m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+fSource(alpha)     !矩空间碰撞
          enddo

          f_post(i,j,0) = m_post(0)/9.0d0-m_post(1)/9.0d0+m_post(2)/9.0d0                                         !这边是乘以M逆
          f_post(i,j,1) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+m_post(3)/6.0d0-m_post(4)/6.0d0 &
                    +m_post(7)/4.0d0
          f_post(i,j,2) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    +m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)/4.0d0
          f_post(i,j,3) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0-m_post(3)/6.0d0+m_post(4)/6.0d0 &
                    +m_post(7)/4.0d0
          f_post(i,j,4) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    -m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)/4.0d0
          f_post(i,j,5) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)/4.0d0
          f_post(i,j,6) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)/4.0d0
          f_post(i,j,7) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)/4.0d0
          f_post(i,j,8) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)/4.0d0

        enddo
    enddo
    return
  end subroutine collision
!===================================================================================================
! collision 结束: 流场碰撞步骤完成。
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
    integer(kind=4) :: alpha
    integer(kind=4) :: ii, jj
    real(kind=8) :: value
    
    !$acc parallel loop gang vector collapse(2) default(none) &
    !$acc& present(f,f_post,ex,ey,streamInterpValidX,streamInterpValidY) &
    !$acc& present(streamInterpIndexX,streamInterpIndexY,streamInterpWeightX,streamInterpWeightY) &
    !$acc& async(1) private(alpha,ii,jj,value)
    do j = 1, ny
        do i = 1, nx
            do alpha = 0, 8
                if((ex(alpha).EQ.0).AND.streamInterpValidY(alpha,j)) then
                    value = 0.0d0
                    do jj = 1, 3
                        value = value + streamInterpWeightY(alpha,j,jj)*f_post(i,streamInterpIndexY(alpha,j,jj),alpha)
                    enddo
                    f(i,j,alpha) = value
                elseif((ey(alpha).EQ.0).AND.streamInterpValidX(alpha,i)) then
                    value = 0.0d0
                    do ii = 1, 3
                        value = value + streamInterpWeightX(alpha,i,ii)*f_post(streamInterpIndexX(alpha,i,ii),j,alpha)
                    enddo
                    f(i,j,alpha) = value
                elseif(streamInterpValidX(alpha,i).AND.streamInterpValidY(alpha,j)) then
                    value = 0.0d0
                    do jj = 1, 3
                        do ii = 1, 3
                            value = value + streamInterpWeightX(alpha,i,ii)*streamInterpWeightY(alpha,j,jj)* &
                                f_post(streamInterpIndexX(alpha,i,ii),streamInterpIndexY(alpha,j,jj),alpha)
                        enddo
                    enddo
                    f(i,j,alpha) = value
                else
                    f(i,j,alpha) = 0.0d0
                endif
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
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine macro()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !$acc parallel loop gang vector collapse(2) default(none) present(f,rho,u,v,Fx,Fy) async(1)
    do j = 1, ny
        do i = 1, nx
            rho(i,j) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
            ! 含力LBM的半步动量修正: rho*u = sum(f e) + 0.5*F, 对应Guo forcing的二阶定义。
            u(i,j) = (f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8) + &
              0.5d0*Fx(i,j))/rho(i,j)
            v(i,j) = ( f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)+0.5d0*Fy(i,j) )/rho(i,j)
        enddo
    enddo
    return
  end subroutine macro
!===================================================================================================
! macro 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
!===================================================================================================

!===================================================================================================



!===================================================================================================
! 子程序: collisionT
! 作用: 完成温度分布函数 g 的碰撞更新，并加入热流修正项。
! 用途: 在主程序时间推进循环中调用，位于流场 macro 之后。
!===================================================================================================
    subroutine collisionT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: n(0:4), n_post(0:4), neq(0:4)
    real(kind=8) :: q(0:4)
    real(kind=8) :: Bx, By
    real(kind=8) :: dBx, dBy
    real(kind=8), parameter :: SG = 1.0d0 - 0.5d0*Qk




    !$acc parallel loop gang vector collapse(2) default(none) present(g,g_post,u,v,T,Bx_prev,By_prev) async(1) &
    !$acc& private(alpha,n,neq,q,n_post,Bx,By,dBx,dBy)
    do j = 1, ny
        do i = 1, nx

            Bx = u(i,j) * T(i,j)
            By = v(i,j) * T(i,j)

#ifdef EnableUseG
            dBx = Bx - Bx_prev(i,j)
            dBy = By - By_prev(i,j)
#else
            dBx = 0.0d0
            dBy = 0.0d0
#endif

#ifdef EnableUseG
            Bx_prev(i,j) = Bx
            By_prev(i,j) = By
#endif

          n(0) = g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
          n(1) = g(i,j,1)-g(i,j,3)
          n(2) = g(i,j,2)-g(i,j,4)
          n(3) = -4.0d0*g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
          n(4) = g(i,j,1)-g(i,j,2)+g(i,j,3)-g(i,j,4)
        
          neq(0) = T(i,j)
          neq(1) = T(i,j)*u(i,j)
          neq(2) = T(i,j)*v(i,j)
#ifdef EnableLegacyThermalScheme
          neq(3) = T(i,j)*paraA
#else
          neq(3) = T(i,j)*(-2.0d0/3.0d0)
#endif
          neq(4) = 0.0d0
        
          q(0) = 0.0d0
          q(1) = Qk
          q(2) = Qk
          q(3) = Qnu
          q(4) = Qnu
        
          
          n_post(0) = n(0)-q(0)*(n(0)-neq(0))
          n_post(1) = n(1)-q(1)*(n(1)-neq(1))+ SG*dBx
          n_post(2) = n(2)-q(2)*(n(2)-neq(2))+ SG*dBy
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
! collisionT 结束: 完成温度分布函数 g 的碰撞更新，并加入热流修正项。
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
    integer(kind=4) :: alpha
    integer(kind=4) :: ii, jj
    real(kind=8) :: value
    
    !$acc parallel loop gang vector collapse(2) default(none) &
    !$acc& present(g,g_post,ex,ey,streamInterpValidX,streamInterpValidY) &
    !$acc& present(streamInterpIndexX,streamInterpIndexY,streamInterpWeightX,streamInterpWeightY) &
    !$acc& async(1) private(alpha,ii,jj,value)
    do j = 1, ny
        do i = 1, nx
            do alpha = 0, 4
                if((ex(alpha).EQ.0).AND.streamInterpValidY(alpha,j)) then
                    value = 0.0d0
                    do jj = 1, 3
                        value = value + streamInterpWeightY(alpha,j,jj)*g_post(i,streamInterpIndexY(alpha,j,jj),alpha)
                    enddo
                    g(i,j,alpha) = value
                elseif((ey(alpha).EQ.0).AND.streamInterpValidX(alpha,i)) then
                    value = 0.0d0
                    do ii = 1, 3
                        value = value + streamInterpWeightX(alpha,i,ii)*g_post(streamInterpIndexX(alpha,i,ii),j,alpha)
                    enddo
                    g(i,j,alpha) = value
                elseif(streamInterpValidX(alpha,i).AND.streamInterpValidY(alpha,j)) then
                    value = 0.0d0
                    do jj = 1, 3
                        do ii = 1, 3
                            value = value + streamInterpWeightX(alpha,i,ii)*streamInterpWeightY(alpha,j,jj)* &
                                g_post(streamInterpIndexX(alpha,i,ii),streamInterpIndexY(alpha,j,jj),alpha)
                        enddo
                    enddo
                    g(i,j,alpha) = value
                else
                    g(i,j,alpha) = 0.0d0
                endif
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
#ifdef EnableLegacyThermalScheme
        g(1,j,1) = -g_post(1,j,3)+(4.0d0+paraA)/10.0d0*Thot
#else
        g(1,j,1) = -g_post(1,j,3)+2.0d0*omegaT(3)*Thot
#endif
        !Right boundary
#ifdef EnableLegacyThermalScheme
        g(nx,j,3) = -g_post(nx,j,1)+(4.0d0+paraA)/10.0d0*Tcold
#else
        g(nx,j,3) = -g_post(nx,j,1)+2.0d0*omegaT(1)*Tcold
#endif
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
#ifdef EnableLegacyThermalScheme
        g(i,1,2) = -g_post(i,1,4)+(4.0d0+paraA)/10.0d0*Thot
#else
        g(i,1,2) = -g_post(i,1,4)+2.0d0*omegaT(4)*Thot
#endif
        !Top side
#ifdef EnableLegacyThermalScheme
        g(i,ny,4) = -g_post(i,ny,2)+(4.0d0+paraA)/10.0d0*Tcold
#else
        g(i,ny,4) = -g_post(i,ny,2)+2.0d0*omegaT(2)*Tcold
#endif
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
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
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
            rho(i,j) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
            momx = f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)
            momy = f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)

            if (rho(i,j).GT.0.0d0) then
                u(i,j) = momx/rho(i,j)
                v(i,j) = momy/rho(i,j)

                do iter = 1, 3
                    Fx(i,j) = 0.0d0
                    Fy(i,j) = rho(i,j)*gBeta*(T(i,j)-Tref)

#ifdef    SideHeatedHa
                    Fx(i,j) = B2sigemarho*(v(i,j)*sin(phi)*cos(phi)-u(i,j)*sin(phi)*sin(phi))
                    Fy(i,j) = rho(i,j)*gBeta*(T(i,j)-Tref)+rho(i,j)*B2sigemarho*(u(i,j)*sin(phi)*cos(phi)&
                    -v(i,j)*cos(phi)*cos(phi))
#endif

                    u(i,j) = (momx + 0.5d0*Fx(i,j))/rho(i,j)
                    v(i,j) = (momy + 0.5d0*Fy(i,j))/rho(i,j)
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
! reconstruct_macro_from_fg 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
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
    
    !$acc parallel loop collapse(2) default(none) present(u,up,v,vp,T,Tp) &
    !$acc& reduction(+:error1,error2,error5,error6)
    do j = 1, ny
        do i = 1, nx
            error1 = error1+((u(i,j)-up(i,j))*(u(i,j)-up(i,j)) + &
                (v(i,j)-vp(i,j))*(v(i,j)-vp(i,j)))
            error2 = error2+(u(i,j)*u(i,j)+v(i,j)*v(i,j))
                
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

    
    write(caseTag,'("Ra=",ES24.16E3,",nx=",I0,",ny=",I0,",useG=",L1,",old=",L1)') Rayleigh, nx, ny, useG,&
    &useLegacyThermalScheme  !输出收敛曲线的对比
    call append_convergence_master_tecplot('convergence_all_2DOpenacc.plt', caseTag, restartItcOffset+itc, errorU, errorT)

    write(*,'(I12,1X,ES24.16E3,1X,ES24.16E3)') restartItcOffset+itc, errorU, errorT


    return
    end subroutine check
!===================================================================================================
! check 结束: 计算稳态收敛误差并写出收敛历史。
!===================================================================================================

!===================================================================================================
! 子程序: output_steady_monitor
! 作用: 在 steadyFlow 中按 check 的间隔输出当前误差和瞬时体平均 Nu/Re。
! 说明: 只做诊断, 不推进 dimensionlessTime, 也不写 Nu_VolAvg/Re_VolAvg 正式采样文件。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
    subroutine output_steady_monitor()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer :: monitorUnit
    real(kind=8) :: NuVolAvg_temp, ReVolAvg_temp
    real(kind=8) :: NuVolAvg_inst, ReVolAvg_inst
    real(kind=8) :: areaWeight
    logical :: monitorFileExists
    logical, save :: first_monitor_write = .true.

    !$acc wait(1)
    NuVolAvg_temp = 0.0d0
#ifdef SideHeatedCell
    !$acc parallel loop collapse(2) default(none) present(u,T,quadWidthX,quadWidthY) &
    !$acc& private(areaWeight) reduction(+:NuVolAvg_temp)
    do j = 1, ny
        do i = 1, nx
            areaWeight = quadWidthX(i)*quadWidthY(j)
            NuVolAvg_temp = NuVolAvg_temp+areaWeight*u(i,j)*(T(i,j)-Tref)
        enddo
    enddo
#endif

#ifdef RayleighBenardCell
    !$acc parallel loop collapse(2) default(none) present(v,T,quadWidthX,quadWidthY) &
    !$acc& private(areaWeight) reduction(+:NuVolAvg_temp)
    do j = 1, ny
        do i = 1, nx
            areaWeight = quadWidthX(i)*quadWidthY(j)
            NuVolAvg_temp = NuVolAvg_temp+areaWeight*v(i,j)*(T(i,j)-Tref)
        enddo
    enddo
#endif

    ReVolAvg_temp = 0.0d0
    !$acc parallel loop collapse(2) default(none) present(u,v,quadWidthX,quadWidthY) &
    !$acc& private(areaWeight) reduction(+:ReVolAvg_temp)
    do j = 1, ny
        do i = 1, nx
            areaWeight = quadWidthX(i)*quadWidthY(j)
            ReVolAvg_temp = ReVolAvg_temp+areaWeight*(u(i,j)*u(i,j)+v(i,j)*v(i,j))
        enddo
    enddo

    NuVolAvg_inst = NuVolAvg_temp/quadSumArea*lengthUnit/diffusivity+1.0d0
    ReVolAvg_inst = dsqrt(ReVolAvg_temp/quadSumArea)*lengthUnit/viscosity

    if(first_monitor_write) then
        inquire(file="steady_monitor.dat", exist=monitorFileExists)
        if((loadInitField.EQ.0).OR.(.not.monitorFileExists)) then
            open(newunit=monitorUnit,file="steady_monitor.dat",status="replace",action="write")
            write(monitorUnit,'(A)') "# itc errorU errorT NuVolAvg ReVolAvg"
        else
            open(newunit=monitorUnit,file="steady_monitor.dat",status="old",position="append",action="write")
        endif
        first_monitor_write = .false.
    else
        open(newunit=monitorUnit,file="steady_monitor.dat",status="unknown",position="append",action="write")
    endif

    write(monitorUnit,'(I12,4(1X,ES24.16E3))') restartItcOffset+itc, &
        real(errorU,kind=8), real(errorT,kind=8), real(NuVolAvg_inst,kind=8), real(ReVolAvg_inst,kind=8)
    close(monitorUnit)

    return
    end subroutine output_steady_monitor
!===================================================================================================
! output_steady_monitor 结束: 已输出 steadyFlow 诊断时间序列。
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
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
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
! 子程序: output_SnapshotMeshFile
! 作用: 输出所有快照共用的非均匀网格坐标。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine output_SnapshotMeshFile()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    open(unit=07,file=trim(snapshotFilePrefix)//"-mesh.dat",status="replace",action="write")
    write(07,'(A)') "# ISLBM nonuniform mesh coordinates for snapshot binary files"
    write(07,'(A,1X,I0,1X,I0)') "nx_ny", nx, ny
    write(07,'(A)') "# xp"
    do i = 1, nx
        write(07,'(I8,1X,ES24.16E3)') i, real(xp(i),kind=8)
    enddo
    write(07,'(A)') "# yp"
    do j = 1, ny
        write(07,'(I8,1X,ES24.16E3)') j, real(yp(j),kind=8)
    enddo
    close(07)

    return
  end subroutine output_SnapshotMeshFile
!===================================================================================================
! output_SnapshotMeshFile 结束: 已输出所有快照共用的非均匀坐标。
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
    ! Strict restart snapshots store f and g; EnableUseG also stores the previous heat-flux history.
    write(05) (((real(f(i,j,alpha),kind=8), i=1,nx), j=1,ny), alpha=0,8)
    write(05) (((real(g(i,j,alpha),kind=8), i=1,nx), j=1,ny), alpha=0,4)
#ifdef EnableUseG
    write(05) ((real(Bx_prev(i,j),kind=8), i=1,nx), j=1,ny)
    write(05) ((real(By_prev(i,j),kind=8), i=1,nx), j=1,ny)
#endif
    close(05)
    call write_reload_metadata(trim(filename))
    
    return
  end subroutine output_ReloadFile
!===================================================================================================
! output_ReloadFile 结束: 输出包含 f、g 的严格重启备份文件。
!===================================================================================================


!===================================================================================================
! 子程序: write_reload_metadata
! 作用: 覆盖写出最新 reload 续算账本，恢复累计步数、t_ff、输出编号和最新 .bin 文件名。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
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
    write(metaUnit,'(A,1X,I0)') 'reload_meta_version', 3
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
    close(metaUnit)

    return
  end subroutine write_reload_metadata
!===================================================================================================
! write_reload_metadata 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
!===================================================================================================


!===================================================================================================
! 子程序: read_reload_metadata
! 作用: 优先读取 latest .meta；若没有，则根据手工编号做保守推断。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
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
    integer(kind=4) :: metaItc, metaSnapshotFileNum, metaPltFileNum, metaReloadFileNum
    real(kind=8) :: metaTf
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
    if((ios.NE.0).OR.(trim(label).NE.'reload_meta_version').OR.(metaVersion.NE.3)) then
        write(*,*) 'Error: invalid reload metadata version in ', trim(metaFile)
        stop    !如果读取失败，或者标签不是 reload_meta_version，或者版本号不是 3，就说明文件格式不对，停止。
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
    
    read(metaUnit,*,iostat=ios) label, metaItc
    if((ios.NE.0).OR.(trim(label).NE.'itc_total')) then
        write(*,*) 'Error: invalid itc_total entry in ', trim(metaFile)
        stop
    endif

    read(metaUnit,*,iostat=ios) label, metaTf
    if((ios.NE.0).OR.(trim(label).NE.'time_tf')) then
        write(*,*) 'Error: invalid time_tf entry in ', trim(metaFile)
        stop
    endif

    read(metaUnit,*,iostat=ios) label, metaSnapshotFileNum
    if((ios.NE.0).OR.(trim(label).NE.'snapshotFileNum')) then
        write(*,*) 'Error: invalid snapshotFileNum entry in ', trim(metaFile)
        stop
    endif

    read(metaUnit,*,iostat=ios) label, metaPltFileNum
    if((ios.NE.0).OR.(trim(label).NE.'pltFileNum')) then
        write(*,*) 'Error: invalid pltFileNum entry in ', trim(metaFile)
        stop
    endif

    read(metaUnit,*,iostat=ios) label, metaReloadFileNum
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

    restartItcOffset = metaItc
    reloadDimensionlessTime = metaTf
    snapshotFileNum = metaSnapshotFileNum
    pltFileNum = metaPltFileNum
    ! reloadFileNum 是整数计数器，给后续 output_ReloadFile() 继续编号，避免覆盖旧 reload 文件。
    ! reloadFileName 是字符串文件名，本次续算马上用它打开 <reloadFilePrefix>-<reloadFileName>.bin。
    reloadFileNum = metaReloadFileNum
    reloadFileName = trim(metaReloadFileName)     ! 处理右边空格，删掉尾部空格
    reloadMetadataLoaded = .true.

    return
  end subroutine read_reload_metadata
!===================================================================================================
! read_reload_metadata 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
!===================================================================================================


!===================================================================================================
! 子程序: infer_reload_offsets_without_metadata
! 作用: 没有 latest .meta 时，只能根据文件编号和当前手工参数推断。
! 根据文件名编号和当前参数“猜一个合理值”，保证续算的时间/步数是连续的。
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine infer_reload_offsets_without_metadata()
    use commondata
    implicit none

    restartItcOffset = 0
#ifdef steadyFlow
    restartItcOffset = max(0, reloadFileNum)  !稳态的 reload 文件名本来就是用 itc 写的
    if(reloadDimensionlessTime.EQ.0.0d0) then !如果没有手工给 reloadDimensionlessTime,计算一下
        reloadDimensionlessTime = real(restartItcOffset,kind=8) / timeUnit
    endif
#endif
#ifdef unsteadyFlow
    if(reloadDimensionlessTime.EQ.0.0d0) then   !非稳态的 reload 文件名不是 itc，而是第几个 reload 文件
        reloadDimensionlessTime = real(max(0,reloadFileNum),kind=8) * reloadFileInterval
    endif
    restartItcOffset = max(0, int(reloadDimensionlessTime*timeUnit+0.5d0))  !再反推格子步
    snapshotFileNum = max(0, int(reloadDimensionlessTime/outputSnapshotInterval+0.5d0)) !推断输出编号
    pltFileNum = max(0, int(reloadDimensionlessTime/outputPltFileInterval+0.5d0))
#endif

    return
  end subroutine infer_reload_offsets_without_metadata
!===================================================================================================
! infer_reload_offsets_without_metadata 结束: 已完成本子程序对应的计算或数据处理步骤。
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

    !----zone ------------------------------------------------------------
    ! 再写一次 299.0: Tecplot 规范里 zone header 之后的数据描述块会以 zone marker 开始。
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
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
!===================================================================================================
  subroutine calNuRe()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: NuVolAvg_temp    !体平均 Nu 的对流热通量积分
    real(kind=8) :: ReVolAvg_temp    !体平均 Re 的速度模积分
    real(kind=8) :: areaWeight
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

    dimensionlessTime = dimensionlessTime+1   !每隔 outputSnapshotInterval 个自由落体时间调用一次calNuRe
#ifdef steadyFlow
    sampleTime = real(restartItcOffset+itc,kind=8)
#endif
#ifdef unsteadyFlow
    sampleTime = reloadDimensionlessTime + real(dimensionlessTime,kind=8)*outputSnapshotInterval
#endif

    if((first_nure_write).AND.(loadInitField.EQ.1)) then
        inquire(file="Nu_VolAvg_2DOpenacc.dat", exist=exNu)
        inquire(file="Re_VolAvg_2DOpenacc.dat", exist=exRe)
        if((.not.exNu).OR.(.not.exRe)) then
            write(*,*) "Error: restart requested but old Nu/Re time-series files are missing."
            open(unit=00,file=trim(settingsFile),status="unknown",position="append")
            write(00,*) "Error: restart requested but old Nu/Re time-series files are missing."
            write(00,*) "Nu_VolAvg_2DOpenacc.dat exists =", exNu
            write(00,*) "Re_VolAvg_2DOpenacc.dat exists =", exRe
            close(00)
            stop
        endif
    endif


    
    NuVolAvg_temp = 0.0d0    
#ifdef SideHeatedCell  
    !$acc parallel loop collapse(2) default(none) present(u,T,quadWidthX,quadWidthY) private(areaWeight) reduction(+:NuVolAvg_temp)
    do j = 1, ny
        do i = 1, nx
            areaWeight = quadWidthX(i)*quadWidthY(j)
            NuVolAvg_temp = NuVolAvg_temp+areaWeight*u(i,j)*(T(i,j)-Tref)     !非均匀网格 midpoint-rule 加权对流热通量
        enddo
    enddo
#endif

#ifdef RayleighBenardCell  
    !$acc parallel loop collapse(2) default(none) present(v,T,quadWidthX,quadWidthY) private(areaWeight) reduction(+:NuVolAvg_temp)
    do j = 1, ny
        do i = 1, nx
            areaWeight = quadWidthX(i)*quadWidthY(j)
            NuVolAvg_temp = NuVolAvg_temp+areaWeight*v(i,j)*(T(i,j)-Tref)     !非均匀网格 midpoint-rule 加权对流热通量
        enddo
    enddo
#endif


    NuVolAvg(dimensionlessTime) = NuVolAvg_temp/quadSumArea*lengthUnit/diffusivity+1.0d0    !!ISLBM体平均 Nusselt 数：非均匀网格体积分加权

    if((first_nure_write).AND.(loadInitField.EQ.0)) then
        open(unit=01,file="Nu_VolAvg_2DOpenacc.dat",status='replace',action='write')
    else
        open(unit=01,file="Nu_VolAvg_2DOpenacc.dat",status='unknown',position='append',action='write')
    endif
    write(01,'(ES24.16E3,1X,ES24.16E3)') &
        real(sampleTime,kind=8), &
        real(NuVolAvg(dimensionlessTime),kind=8)   !以格子步数或者自由落体时间来写入
    close(01)

    ReVolAvg_temp = 0.0d0
    !$acc parallel loop collapse(2) default(none) present(u,v,quadWidthX,quadWidthY) private(areaWeight) reduction(+:ReVolAvg_temp)
    do j = 1, ny
        do i = 1, nx 
            areaWeight = quadWidthX(i)*quadWidthY(j)
            ReVolAvg_temp = ReVolAvg_temp+areaWeight*(u(i,j)*u(i,j)+v(i,j)*v(i,j))
        enddo
    enddo
    ReVolAvg(dimensionlessTime) = dsqrt(ReVolAvg_temp/quadSumArea)*lengthUnit/viscosity    !ISLBM全域体平均 Reynolds 数：速度RMS非均匀加权平均


    if((first_nure_write).AND.(loadInitField.EQ.0)) then
        open(unit=02,file="Re_VolAvg_2DOpenacc.dat",status='replace',action='write')
    else
        open(unit=02,file="Re_VolAvg_2DOpenacc.dat",status='unknown',position='append',action='write')
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
! 子程序: output_unsteady_NuRe_postprocess
! 作用: 从完整 Nu/Re 历史文件重建非稳态时间序列、运行平均值和分段窗口平均值。
! 用途: 非稳态算例结束后用于输出后处理统计文件。
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

    inquire(file='Nu_VolAvg_2DOpenacc.dat', exist=exNu)
    inquire(file='Re_VolAvg_2DOpenacc.dat', exist=exRe)
    if((.not.exNu).or.(.not.exRe)) then
        write(*,'(A)') 'Error: Nu/Re history files are missing before postprocessing.'
        open(unit=00,file=trim(settingsFile),status='unknown',position='append')
        write(00,'(A)') 'Error: Nu/Re history files are missing before postprocessing.'
        close(00)
        error stop 1
    endif

    open(newunit=nuUnit, file='Nu_VolAvg_2DOpenacc.dat', status='old', action='read', form='formatted')
    open(newunit=reUnit, file='Re_VolAvg_2DOpenacc.dat', status='old', action='read', form='formatted')

    ! These files are derived views of the full .dat history, so rebuild one continuous ZONE.
    open(newunit=seriesUnit, file='NuRe_VolAvg_2DOpenacc.plt', status='replace', action='write', form='formatted')
    write(seriesUnit,'(A)') 'TITLE = "2D OpenACC Nu/Re volume averages"'
    write(seriesUnit,'(A)') 'VARIABLES = "time" "NuVolAvg" "ReVolAvg"'
    write(seriesUnit,'(A)') 'ZONE T="NuReVolAvg", F=POINT'

    open(newunit=runningUnit, file='NuRe_VolAvg_runningMean_2DOpenacc.plt', status='replace', action='write', form='formatted')
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

    open(unit=33, file='NuRe_TimeAverage_2DOpenacc.txt', status='replace', action='write', form='formatted')
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
! output_unsteady_NuRe_postprocess 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
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
    real(kind=8) :: dTdx, qx, sum_qx, areaWeight
    real(kind=8) :: xnode(3), fnode(3)
    real(kind=8) :: deltaT, coef

    deltaT = Thot - Tcold
    coef   = velocityScaleCompare

    sum_qx = 0.0d0

    do j = 1, ny
      do i = 1, nx

        if(i.EQ.1) then
            xnode = (/ xp(1), xp(2), xp(3) /)
            fnode = (/ T(1,j), T(2,j), T(3,j) /)
        elseif(i.EQ.nx) then
            xnode = (/ xp(nx-2), xp(nx-1), xp(nx) /)
            fnode = (/ T(nx-2,j), T(nx-1,j), T(nx,j) /)
        else
            xnode = (/ xp(i-1), xp(i), xp(i+1) /)
            fnode = (/ T(i-1,j), T(i,j), T(i+1,j) /)
        endif
        call lagrange_derivative_3(xnode, fnode, xp(i), dTdx)
        qx = coef*u(i,j)*(T(i,j)-Tref) - dTdx
        areaWeight = quadWidthX(i)*quadWidthY(j)
        sum_qx = sum_qx + areaWeight*qx

      enddo
    enddo
    Nu_global = (sum_qx / quadSumArea) / deltaT

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
    integer(kind=4) :: j
    integer(kind=4) :: jmax, jmin
    integer(kind=4) :: midIndex(3)
    real(kind=8) :: deltaT,coef
    real(kind=8) :: qx_wall, sum_hot, sum_cold, sum_mid
    real(kind=8) :: u_mid, T_mid, dTdx_mid
    real(kind=8) :: xnode(3), fnode(3)
    real(kind=8) :: xmidNode(3), xmidWeight(3)
    logical :: midStencilValid
    real(kind=8) :: Nu_left(1:ny) 
    !------------------------------------------------------------
    ! 5-point least-squares parabola fit (general, can be one-sided)
    !------------------------------------------------------------
    integer(kind=4) :: k
    integer(kind=4) :: jj(5)
    real(kind=8) :: fk(5), yk(5)
    real(kind=8) :: fstar, ystar
    real(kind=8) :: Nu_left_ext(0:ny+1)


    deltaT = Thot - Tcold
    coef   = velocityScaleCompare

   

    !-----------------------------
    ! (1) 左侧热壁平均 Nu_hot，同时记录 Numax/Numin 及其 y 位置
    sum_hot = 0.0d0
    do j = 1, ny
      ! 壁面导热通量：qx(x=0,j)
      xnode = (/ 0.0d0, xp(1), xp(2) /)
      fnode = (/ Thot, T(1,j), T(2,j) /)
      call lagrange_derivative_3(xnode, fnode, 0.0d0, qx_wall)
      qx_wall = -qx_wall
      Nu_left(j)= qx_wall / deltaT
      sum_hot   = sum_hot + Nu_left(j)*quadWidthY(j)
    enddo
    Nu_hot = sum_hot / quadSumY

    ! 角点没有唯一的壁面法向Nu定义; 这里只把相邻壁面Nu复制到扩展数组,
    ! 用于后面的极值搜索和五点拟合, 避免重新引入均匀网格角点外推公式。
    Nu_left_ext(1:ny) = Nu_left(1:ny)
    Nu_left_ext(0) = Nu_left(1)
    Nu_left_ext(ny+1) = Nu_left(ny)

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
      xnode = (/ xp(nx-1), xp(nx), 1.0d0 /)
      fnode = (/ T(nx-1,j), T(nx,j), Tcold /)
      call lagrange_derivative_3(xnode, fnode, 1.0d0, qx_wall)
      qx_wall = -qx_wall
      sum_cold = sum_cold + qx_wall/deltaT*quadWidthY(j)
    enddo
    Nu_cold = sum_cold / quadSumY

    !-----------------------------
    ! (3) 竖直中线 x=1/2 的平均 Nu_middle
    sum_mid = 0.0d0

    call build_lagrange_stencil_1d(nx, xp(1:nx), 0.5d0, midIndex, xmidWeight, midStencilValid)
    if(.not.midStencilValid) then
      write(*,*) "Error: x=0.5 is outside ISLBM x nodes in SideHeatedcalc_Nu_wall_avg"
      stop
    endif
    xmidNode = (/ xp(midIndex(1)), xp(midIndex(2)), xp(midIndex(3)) /)

    do j = 1, ny
      u_mid = xmidWeight(1)*u(midIndex(1),j) + xmidWeight(2)*u(midIndex(2),j) + &
        xmidWeight(3)*u(midIndex(3),j)
      T_mid = xmidWeight(1)*T(midIndex(1),j) + xmidWeight(2)*T(midIndex(2),j) + &
        xmidWeight(3)*T(midIndex(3),j)
      fnode = (/ T(midIndex(1),j), T(midIndex(2),j), T(midIndex(3),j) /)
      call lagrange_derivative_3(xmidNode, fnode, 0.5d0, dTdx_mid)
      sum_mid = sum_mid + (coef*u_mid*(T_mid-Tref) - dTdx_mid)/deltaT*quadWidthY(j)
    enddo

    Nu_middle = sum_mid / quadSumY

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
    integer(kind=4) :: midIndex(3)
    integer(kind=4) :: j0
    integer(kind=4) :: jj(5)
    real(kind=8) :: uline(1:ny)
    real(kind=8) :: s(5), fu(5)
    real(kind=8) :: umax_grid, umax_fit, y_fit
    real(kind=8) :: xmid, xmidWeight(3)
    logical :: midStencilValid
    character(len=24) :: ctime, string
    integer(kind=4) :: time
    real(kind=8) :: coef

    coef = velocityScaleCompare

    ! ---- (1) 构造中线剖面 u(x=1/2, y_j) ----
    xmid = 0.5d0
    call build_lagrange_stencil_1d(nx, xp(1:nx), xmid, midIndex, xmidWeight, midStencilValid)
    if(.not.midStencilValid) then
      write(*,*) "Error: xmid is outside ISLBM x nodes in SideHeatedcalc_umid_max"
      stop
    endif
    do j = 1, ny
      uline(j) = xmidWeight(1)*u(midIndex(1),j) + xmidWeight(2)*u(midIndex(2),j) + &
        xmidWeight(3)*u(midIndex(3),j)
    enddo

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
    integer(kind=4) :: midIndex(3)
    integer(kind=4) :: i0
    integer(kind=4) :: ii(5)
    real(kind=8) :: vline(1:nx)
    real(kind=8) :: s(5), fv(5)
    real(kind=8) :: vmax_grid, vmax_fit, x_fit
    real(kind=8) :: ymid, ymidWeight(3)
    logical :: midStencilValid
    character(len=24) :: ctime, string
    integer(kind=4) :: time
    real(kind=8) :: coef

    coef = velocityScaleCompare

    ! ---- (1) 构造中线剖面 v(x_i, y=1/2) ----
    ymid = 0.5d0
    call build_lagrange_stencil_1d(ny, yp(1:ny), ymid, midIndex, ymidWeight, midStencilValid)
    if(.not.midStencilValid) then
      write(*,*) "Error: ymid is outside ISLBM y nodes in SideHeatedcalc_vmid_max"
      stop
    endif
    do i = 1, nx
      vline(i) = ymidWeight(1)*v(i,midIndex(1)) + ymidWeight(2)*v(i,midIndex(2)) + &
        ymidWeight(3)*v(i,midIndex(3))
    enddo

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
  real(kind=8) :: dTdy, qy, sum_qy, areaWeight
  real(kind=8) :: ynode(3), fnode(3)
  real(kind=8) :: deltaT, coef

  deltaT = Thot - Tcold
  coef   = velocityScaleCompare

  sum_qy = 0.0d0

  do j = 1, ny
    do i = 1, nx

      if(j.EQ.1) then
        ynode = (/ yp(1), yp(2), yp(3) /)
        fnode = (/ T(i,1), T(i,2), T(i,3) /)
      elseif(j.EQ.ny) then
        ynode = (/ yp(ny-2), yp(ny-1), yp(ny) /)
        fnode = (/ T(i,ny-2), T(i,ny-1), T(i,ny) /)
      else
        ynode = (/ yp(j-1), yp(j), yp(j+1) /)
        fnode = (/ T(i,j-1), T(i,j), T(i,j+1) /)
      endif
      call lagrange_derivative_3(ynode, fnode, yp(j), dTdy)

      qy = coef * v(i,j) * (T(i,j) - Tref) - dTdy
      areaWeight = quadWidthX(i)*quadWidthY(j)
      sum_qy = sum_qy + areaWeight*qy

    enddo
  enddo
  Nu_global = (sum_qy / quadSumArea) / deltaT

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
  integer(kind=4) :: midIndex(3)
  integer(kind=4) :: ii(5)
  real(kind=8) :: deltaT
  real(kind=8) :: qy_wall, sum_hot, sum_cold, sum_mid, coef
  real(kind=8) :: v_mid, T_mid, dTdy_mid
  real(kind=8), dimension(1:nx) :: Nu_bot
  real(kind=8), dimension(0:nx+1) :: Nu_bot_ext
  real(kind=8) :: xfit(4), Tfit(4), T_wl, T_wr, T_wl2, T_wr2
  real(kind=8) :: xk(5), fk(5)
  real(kind=8) :: ynode(3), fnode(3)
  real(kind=8) :: ymidNode(3), ymidWeight(3)
  real(kind=8) :: fstar, xstar
  logical :: midStencilValid


  deltaT = Thot - Tcold
  coef   = velocityScaleCompare

  !-----------------------------
  ! (1) 底部热壁平均 Nu_hot（不含角点）
  sum_hot = 0.0d0
  do i = 1, nx
    ynode = (/ 0.0d0, yp(1), yp(2) /)
    fnode = (/ Thot, T(i,1), T(i,2) /)
    call lagrange_derivative_3(ynode, fnode, 0.0d0, qy_wall)
    qy_wall = -qy_wall
    Nu_bot(i)= qy_wall / deltaT
    sum_hot  = sum_hot + Nu_bot(i)*quadWidthX(i)
  enddo
  Nu_hot = sum_hot / quadSumX

  !-----------------------------
  ! (1.1) 角点扩展：用侧壁绝热（Neumann）4点拟合得到 x=0 与 x=1 处近壁温度
  ! 左下角附近：i=1..4, j=1 和 j=2
  xfit(1)=xp(1);  Tfit(1)=T(1,1)
  xfit(2)=xp(2);  Tfit(2)=T(2,1)
  xfit(3)=xp(3);  Tfit(3)=T(3,1)
  xfit(4)=xp(4);  Tfit(4)=T(4,1)
  call fit_adiabatic_wall_T4(0.0d0, xfit, Tfit, T_wl)
  xfit(1)=xp(1);  Tfit(1)=T(1,2)
  xfit(2)=xp(2);  Tfit(2)=T(2,2)
  xfit(3)=xp(3);  Tfit(3)=T(3,2)
  xfit(4)=xp(4);  Tfit(4)=T(4,2)
  call fit_adiabatic_wall_T4(0.0d0, xfit, Tfit, T_wl2)

  ! 右下角附近：i=nx-3..nx, j=1 和 j=2
  xfit(1)=xp(nx-3);  Tfit(1)=T(nx-3,1)
  xfit(2)=xp(nx-2);  Tfit(2)=T(nx-2,1)
  xfit(3)=xp(nx-1);  Tfit(3)=T(nx-1,1)
  xfit(4)=xp(nx  );  Tfit(4)=T(nx  ,1)
  call fit_adiabatic_wall_T4(1.0d0, xfit, Tfit, T_wr)
  xfit(1)=xp(nx-3);  Tfit(1)=T(nx-3,2)
  xfit(2)=xp(nx-2);  Tfit(2)=T(nx-2,2)
  xfit(3)=xp(nx-1);  Tfit(3)=T(nx-1,2)
  xfit(4)=xp(nx  );  Tfit(4)=T(nx  ,2)
  call fit_adiabatic_wall_T4(1.0d0, xfit, Tfit, T_wr2)

  ! 组装扩展数组：角点只用于找 max/min 与拟合
  Nu_bot_ext(1:nx) = Nu_bot(1:nx)
  ynode = (/ 0.0d0, yp(1), yp(2) /)
  fnode = (/ Thot, T_wl, T_wl2 /)
  call lagrange_derivative_3(ynode, fnode, 0.0d0, qy_wall)
  Nu_bot_ext(0) = -qy_wall / deltaT
  fnode = (/ Thot, T_wr, T_wr2 /)
  call lagrange_derivative_3(ynode, fnode, 0.0d0, qy_wall)
  Nu_bot_ext(nx+1) = -qy_wall / deltaT

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
    ynode = (/ yp(ny-1), yp(ny), 1.0d0 /)
    fnode = (/ T(i,ny-1), T(i,ny), Tcold /)
    call lagrange_derivative_3(ynode, fnode, 1.0d0, qy_wall)
    qy_wall = -qy_wall
    sum_cold = sum_cold + qy_wall/deltaT*quadWidthX(i)
  enddo
  Nu_cold = sum_cold / quadSumX

  !-----------------------------
  ! 中线的 Nusselt 数：在物理中线 y=0.5 上插值。
  sum_mid = 0.0d0

  call build_lagrange_stencil_1d(ny, yp(1:ny), 0.5d0, midIndex, ymidWeight, midStencilValid)
  if(.not.midStencilValid) then
    write(*,*) "Error: y=0.5 is outside ISLBM y nodes in RBcalc_Nu_wall_avg"
    stop
  endif
  ymidNode = (/ yp(midIndex(1)), yp(midIndex(2)), yp(midIndex(3)) /)

  do i = 1, nx
    v_mid = ymidWeight(1)*v(i,midIndex(1)) + ymidWeight(2)*v(i,midIndex(2)) + &
      ymidWeight(3)*v(i,midIndex(3))
    T_mid = ymidWeight(1)*T(i,midIndex(1)) + ymidWeight(2)*T(i,midIndex(2)) + &
      ymidWeight(3)*T(i,midIndex(3))
    fnode = (/ T(i,midIndex(1)), T(i,midIndex(2)), T(i,midIndex(3)) /)
    call lagrange_derivative_3(ymidNode, fnode, 0.5d0, dTdy_mid)
    sum_mid = sum_mid + (coef*v_mid*(T_mid-Tref) - dTdy_mid)/deltaT*quadWidthX(i)
  enddo

  Nu_middle = sum_mid / quadSumX

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
  integer(kind=4) :: midIndex(3)
  integer(kind=4) :: j0
  integer(kind=4) :: jj(5)
  real(kind=8) :: uline(1:ny)
  real(kind=8) :: yk(5), fk(5)
  real(kind=8) :: umax_fit, y_fit
  real(kind=8) :: coef, xmid, targetX, xmidWeight(3)
  real(kind=8) :: umax_grid, yLen
  logical :: midStencilValid

  coef = velocityScaleCompare

  ! ---- (1) 取论文定义的 x=1/2 竖线剖面 u(x=0.5, y_j) ----
  targetX = 0.5d0
  xmid = targetX

  call build_lagrange_stencil_1d(nx, xp(1:nx), targetX, midIndex, xmidWeight, midStencilValid)
  if(.not.midStencilValid) then
    write(*,*) "Error: targetX is outside ISLBM x nodes in RBcalc_umid_max"
    stop
  endif
  do j = 1, ny
    uline(j) = xmidWeight(1)*u(midIndex(1),j) + xmidWeight(2)*u(midIndex(2),j) + &
      xmidWeight(3)*u(midIndex(3),j)
  enddo

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
  integer(kind=4) :: midIndex(3)
  integer(kind=4) :: i0
  integer(kind=4) :: ii(5)
  real(kind=8) :: vline(1:nx)
  real(kind=8) :: xk(5), fk(5)
  real(kind=8) :: vmax_fit, x_fit
  real(kind=8) :: coef, ymid, ymidWeight(3)
  real(kind=8) :: vmax_grid
  logical :: midStencilValid

  coef = velocityScaleCompare

  ! ---- (1) 取 y=1/2 中线剖面 v(x_i, y=1/2) ----
  ymid = 0.5d0
  call build_lagrange_stencil_1d(ny, yp(1:ny), ymid, midIndex, ymidWeight, midStencilValid)
  if(.not.midStencilValid) then
    write(*,*) "Error: ymid is outside ISLBM y nodes in RBcalc_vmid_max"
    stop
  endif
  do i = 1, nx
    vline(i) = ymidWeight(1)*v(i,midIndex(1)) + ymidWeight(2)*v(i,midIndex(2)) + &
      ymidWeight(3)*v(i,midIndex(3))
  enddo

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
! RBcalc_vmid_max 结束: 已完成本子程序对应的计算或数据处理步骤。
!===================================================================================================
!===================================================================================================







!===================================================================================================
! 子程序: calc_psi_vort_and_output
! 作用: 计算流函数、涡量，并输出相关诊断量。
! 用途: 在主程序结束阶段调用，作为统一后处理的一部分。
!===================================================================================================
subroutine calc_psi_vort_and_output()
  use commondata
  implicit none

  integer(kind=4) :: i, j
  real(kind=8) :: coef
  real(kind=8) :: segmentIntegral
  real(kind=8) :: xNode(3), vNode(3)
  real(kind=8) :: yNode(3), uNode(3)
  real(kind=8) :: dv_dx, du_dy
  real(kind=8) :: psi(nx,ny), vort(nx,ny)
  real(kind=8) :: psiTopAbsMax


  ! for fine-grid max(|psi|)  10001*10001
  real(kind=8) :: psi_abs_max, x_at_max, y_at_max, psi_center_abs_fine


  coef = velocityScaleCompare

  !=========================================================
  ! (A) 计算流函数 psi。
  !     SideHeatedCell 按zzhao后处理从左壁积分: psi(x,y)=-∫_0^x v(mu,y)dmu。
  !     RayleighBenardCell 保留从底壁积分: psi(x,y)=∫_0^y u(x,mu)dmu。
  !     ISLBM非均匀网格下不能使用固定dx/dy的Simpson公式。
  !     每一小段用三点二次Lagrange多项式做解析积分。
  !=========================================================

#ifdef SideHeatedCell
  do j = 1, ny
    xNode = (/ 0.0d0, xp(1), xp(2) /)
    vNode = (/ 0.0d0, v(1,j), v(2,j) /)
    call integrate_lagrange_3_segment(xNode, vNode, 0.0d0, xp(1), segmentIntegral)
    psi(1,j) = -segmentIntegral*coef

    do i = 2, nx
      if (i == 2) then
        xNode = (/ 0.0d0, xp(1), xp(2) /)
        vNode = (/ 0.0d0, v(1,j), v(2,j) /)
      else
        xNode = (/ xp(i-2), xp(i-1), xp(i) /)
        vNode = (/ v(i-2,j), v(i-1,j), v(i,j) /)
      endif
      call integrate_lagrange_3_segment(xNode, vNode, xp(i-1), xp(i), segmentIntegral)
      psi(i,j) = psi(i-1,j) - segmentIntegral*coef
    end do

  end do
#endif

#ifdef RayleighBenardCell
  do i = 1, nx
    yNode = (/ 0.0d0, yp(1), yp(2) /)
    uNode = (/ 0.0d0, u(i,1), u(i,2) /)
    call integrate_lagrange_3_segment(yNode, uNode, 0.0d0, yp(1), segmentIntegral)
    psi(i,1) = segmentIntegral*coef

    do j = 2, ny
      if (j == 2) then
        yNode = (/ 0.0d0, yp(1), yp(2) /)
        uNode = (/ 0.0d0, u(i,1), u(i,2) /)
      else
        yNode = (/ yp(j-2), yp(j-1), yp(j) /)
        uNode = (/ u(i,j-2), u(i,j-1), u(i,j) /)
      endif
      call integrate_lagrange_3_segment(yNode, uNode, yp(j-1), yp(j), segmentIntegral)
      psi(i,j) = psi(i,j-1) + segmentIntegral*coef
    end do

  end do
#endif

  psiTopAbsMax = 0.0d0
#ifdef SideHeatedCell
  do j = 1, ny
    psiTopAbsMax = max(psiTopAbsMax, dabs(psi(nx,j)))
  enddo
  write(*,'(a,1x,es16.8)') "max(|psi_right_internal|) =", psiTopAbsMax
  open(unit=00,file=trim(settingsFile),status="unknown",position="append")
  write(00,'(a,1x,es16.8)') "max(|psi_right_internal|) =", psiTopAbsMax
  close(00)
#endif
#ifdef RayleighBenardCell
  do i = 1, nx
    psiTopAbsMax = max(psiTopAbsMax, dabs(psi(i,ny)))
  enddo
  write(*,'(a,1x,es16.8)') "max(|psi_top_internal|) =", psiTopAbsMax
  open(unit=00,file=trim(settingsFile),status="unknown",position="append")
  write(00,'(a,1x,es16.8)') "max(|psi_top_internal|) =", psiTopAbsMax
  close(00)
#endif


  call output_psi_center_abs(psi)     ! 基于粗网格局部四点插值得到中心点处的 abs(psi)

  !=========================================================
  ! (B) 计算涡量 vort = dv/dx - du/dy（2D）
  !     ISLBM非均匀网格下用三点Lagrange导数, 不再使用固定dx/dy中心差分。
  !=========================================================

  do j = 1, ny
    do i = 1, nx

      if(i.EQ.1) then
        xNode = (/ xp(1), xp(2), xp(3) /)
        vNode = (/ v(1,j), v(2,j), v(3,j) /)
      elseif(i.EQ.nx) then
        xNode = (/ xp(nx-2), xp(nx-1), xp(nx) /)
        vNode = (/ v(nx-2,j), v(nx-1,j), v(nx,j) /)
      else
        xNode = (/ xp(i-1), xp(i), xp(i+1) /)
        vNode = (/ v(i-1,j), v(i,j), v(i+1,j) /)
      endif
      call lagrange_derivative_3(xNode, vNode, xp(i), dv_dx)

      if(j.EQ.1) then
        yNode = (/ yp(1), yp(2), yp(3) /)
        uNode = (/ u(i,1), u(i,2), u(i,3) /)
      elseif(j.EQ.ny) then
        yNode = (/ yp(ny-2), yp(ny-1), yp(ny) /)
        uNode = (/ u(i,ny-2), u(i,ny-1), u(i,ny) /)
      else
        yNode = (/ yp(j-1), yp(j), yp(j+1) /)
        uNode = (/ u(i,j-1), u(i,j), u(i,j+1) /)
      endif
      call lagrange_derivative_3(yNode, uNode, yp(j), du_dy)

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
! 用途: 由主程序、时间推进或后处理流程按需调用，保持与参考 ISLBM 代码的接口风格一致。
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

