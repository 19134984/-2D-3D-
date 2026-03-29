!=============================================================
!!!    注释区，代码描述
!!!    浮力驱动自然对流（上下加热）
!!!    LBM方法
!!!    MRT-LBE
!=============================================================


!=============================================================
!   自定义宏，一些选项的开关
#define steadyFlow    
!~!!#define unsteadyFlow

!   速度边界，包括水平垂直边界无滑移，还有垂直边界速度周期
#define HorizontalWallsNoslip
#define VerticalWallsNoslip
!~#define VerticalWallsPeriodicalU

!   温度边界(for Rayleigh Benard Cell)，包括水平边界恒温，垂直边界温度不可穿透以及周期
!#define RayleighBenardCell
!#define HorizontalWallsConstT
!#define VerticalWallsAdiabatic
!~#define VerticalWallsPeriodicalT


!   温度边界(for Side Heated Cell)，包括水平边界温度不可穿透，垂直边界恒温,侧壁加热加磁场
#define SideHeatedCell
#define HorizontalWallsAdiabatic
#define VerticalWallsConstT
!#define SideHeatedHa          
!~~temperature B.C.~~

!   自定义宏结束
!=============================================================


!=============================================================
!   全局模块
    module commondata
        implicit none
        !===============================================================================================
        ! 是否在计算前从旧算例重启
        integer(kind=4), parameter :: loadInitField=0   ! 0: 不重启；1: 从 backupFile-*.bin 读取初值

        ! 在 loadInitField=1 的前提下：
        integer(kind=4), parameter :: reloadDimensionlessTime=0  ! 旧算例累计的无量纲时间（用于续写 Nu/Re 输出横坐标）
        integer(kind=4), parameter :: reloadbinFileNum=0         ! 读取的备份文件编号：backupFile-<reloadbinFileNum>.bin
        !===============================================================================================

        !===============================================================================================
        ! 无量纲参数
        integer(kind=4), parameter :: nx=22, ny=22     !格子网格
#ifdef SideHeatedCell
        real(kind=8), parameter :: lengthUnit=dble(nx)     !侧壁差温：特征长度取 x 方向长度
#else
        real(kind=8), parameter :: lengthUnit=dble(ny)     !上下差温：特征长度取 y 方向长度
#endif
        real(kind=8), parameter :: pi = acos(-1.0d0)

        real(kind=8), parameter :: Rayleigh=1.0d6        
        real(kind=8), parameter :: Prandtl=0.71d0       
        real(kind=8), parameter :: Mach=0.1d0
        real(kind=8), parameter :: Thot=0.5d0, Tcold=-0.5d0
        !real(kind=8), parameter :: Tref=0.0d0
        real(kind=8), parameter :: Tref=0.5d0*(Thot+Tcold)
        real(kind=8), parameter :: tauf=0.5d0+Mach*lengthUnit*dsqrt(3.0d0*Prandtl/Rayleigh) 
        real(kind=8), parameter :: viscosity=(tauf-0.5d0)/3.0d0
        real(kind=8), parameter :: diffusivity=viscosity/Prandtl
        
        ! heatFluxScale is used in Nu/heat-flux post-processing and should stay consistent with the Nu definition.
        real(kind=8), parameter :: heatFluxScale=lengthUnit/diffusivity          ! 固定采用 L/kappa

        ! velocityScaleCompare is used only in velocity-related post-processing to convert lattice velocity
        ! to the nondimensional velocity scale adopted by the reference paper being compared.
        real(kind=8), parameter :: velocityScaleCompare=lengthUnit/diffusivity   ! 默认采用热扩散标度 UL/kappa；若要按自由落体标度比较，可改为 1.0d0/velocityUnit
        
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
        !real(kind=8), parameter :: Qd=3.0d0-dsqrt(3.0d0), Qnu=4.0d0*dsqrt(3.0d0)-6.0d0             !温度的多松弛系数
        real(kind=8), parameter :: taug = 0.5d0 + (tauf - 0.5d0)/Prandtl
        real(kind=8), parameter :: s_e = 1.0d0, s_q = 1.0d0, s_j = 1.0d0/taug
        !===============================================================================================          
        
        !===============================================================================================
        ! 输出/备份相关设置（以自由落体时间 t_ff 为单位）
        real(kind=8), parameter :: outputFrequency=100.0d0   ! 每隔 outputFrequency 个自由落体时间 t_ff 输出/统计一次（时间间隔）

        integer(kind=4), parameter :: dimensionlessTimeMax=int(12000/outputFrequency)  ! 用于 NuVolAvg/ReVolAvg 数组的最大记录点数；每 outputFrequency*t_ff 记录一次，则最多覆盖约 12000*t_ff

        integer(kind=4), parameter :: backupInterval=1000  ! 备份间隔（自由落体时间单位），为了停电情况下，可以继续计算
        
        real(kind=8), parameter :: epsU=1.0d-7, epsT=1.0d-7    ! 稳态收敛阈值   

        integer(kind=4), parameter :: outputBinFile=0   ! 是否输出 bin 文件：0=不输出，1=输出
        integer(kind=4), parameter :: outputPltFile=0   ! 是否输出 plt 文件：0=不输出，1=输出

        integer(kind=4) :: binFileNum, pltFileNum  ! bin/plt 输出文件的计数器
        ! - unsteadyFlow：每次输出递增（用于文件名编号）
        ! - steadyFlow：bin/plt 文件名通常直接用 itc

        integer(kind=4) :: dimensionlessTime
        ! 统计/输出时间点编号（与 outputFrequency 对应）：
        ! 每调用一次 calNuRe() 就 dimensionlessTime = dimensionlessTime + 1
        ! 用于索引 NuVolAvg/ReVolAvg 数组，并用于输出的时间轴：t = reloadDimensionlessTime + dimensionlessTime*outputFrequency（单位：t_ff）

        real(kind=8) :: NuVolAvg(0:dimensionlessTimeMax), ReVolAvg(0:dimensionlessTimeMax)
        ! 体平均 Nu 和 Re 的时间序列缓存
        ! 只有在启用并调用 calNuRe() 的情况下这些数组才会被真正填充
  
        character(len=100) :: binFolderPrefix="./binFile/buoyancyCavity"
        ! bin 输出文件前缀（实际文件名形如：<binFolderPrefix>-<编号>.bin）

        character(len=100) :: pltFolderPrefix="./pltFile/buoyancyCavity"
        ! plt 输出文件前缀（实际文件名形如：<pltFolderPrefix>-<编号>.plt）

        character(len=100) :: reloadFilePrefix="./reloadFile/backupFile"
        ! 重启读取文件的前缀（实际读取：<reloadFilePrefix>-<reloadbinFileNum>.bin）
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
        real(kind=8), allocatable :: Bx_prev(:,:), By_prev(:,:)
        integer(kind=4) :: itc
        integer(kind=4), parameter :: itc_max=20000000 !格子时间步长
        logical, parameter :: useG = .true.            !M1G 开关
        real(kind=8) :: Nu_global, Nu_hot, Nu_cold, Nu_middle    !平均Nu，全场，侧壁以及中线
        real(kind=8) :: Nu_hot_max, Nu_hot_min, Nu_hot_max_position, Nu_hot_min_position    !左侧壁面的最大最小Nu，以及对应的位置
        
        
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

    use omp_lib
    use commondata
    implicit none
    real(kind=8) :: timeStart, timeEnd
    character(len=24) :: ctime, string
    INTEGER(kind=4) :: time
    real(kind=8) :: timeStart2, timeEnd2
    integer(kind=4) :: myMaxThreads
    

    !===============================================================================================
    !设置并行核数
    open(unit=00,file="SimulationSettings.txt",status='unknown')   !打开（或创建）txt文件，准备写入
    string = ctime( time() )                      !ctime把 time() 返回的时间戳转换成可读的字符串
    write(00,*) 'Start: ', string                 !什么时候开始计算
    write(00,*) "Starting OpenMP >>>>>>"
    call OMP_set_num_threads(24)                   !使用 24 个线程
    myMaxThreads = OMP_get_max_threads()           !查询最大可用线程数
    write(00,*) "Max Running threads:",myMaxThreads
    close(00)
    !===============================================================================================


    !===============================================================================================
    ! Initialization
    call initial()

    !===============================================================================================
    !-----------------------------------------------------------------------------------------------

    call CPU_TIME(timeStart)         !当前进程累计消耗的 CPU 时间,包括并行
    timeStart2 = OMP_get_wtime()     !墙钟时间(实际耗时，不包括并行)
    do while( ((errorU.GT.epsU).OR.(errorT.GT.epsT)).AND.(itc.LE.itc_max) )   !只要 (errorU > epsU 或 errorT > epsT) 且 itc ≤ itc_max，就继续循环
                                                                              !换成if，就是 errorU > epsU .and. errorT > epsT

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
        if(MOD(itc,2000).EQ.0) call check()
#endif
        
         if( MOD(itc, int(outputFrequency*timeUnit)).EQ.0 ) then  ! 达到一个输出间隔outputFrequency，就执行一次计算体平均Nu,是对流通量，不是全场平均Nu

             call calNuRe()
            
#ifdef steadyFlow
             if( (outputPltFile.EQ.1).AND.(MOD(itc, backupInterval*int(timeUnit)).EQ.0) ) call output_Tecplot()  !plt 输出备份间隔uvT
#endif
            
#ifdef unsteadyFlow
             if(outputBinFile.EQ.1) then
                     call output_binary()          !bin 输出后处理间隔uvTrho数据
                     if(MOD(itc, int(backupInterval/outputFrequency)*int(outputFrequency*timeUnit)).EQ.0) call backupData()  !输出备份间隔uvTfg
             endif   
             if(outputPltFile.EQ.1) call output_Tecplot()  !plt 输出后处理uvT快照
#endif
        endif
     enddo

    call CPU_TIME(timeEnd)         !当前进程累计消耗的 CPU 时间,包括并行
    timeEnd2 = OMP_get_wtime()     !墙钟时间(实际耗时，不包括并行)

#ifdef steadyFlow
    call output_Tecplot()          !输出最后一步的plt结果
    call output_binary()              !输出最后一步的bin结果
#endif
    !===============================================================================================



    !===============================================================================================

    
!侧壁加热和RB对流的计算Nu不一样
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


    call calc_psi_vort_and_output()  !输出腔体中心的abs(psi),以及最大的abs(psi)max以及位置（采用细网格插值出来）






     

    open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')        !在这个txt文件后面继续写（追加模式）
    write(00,*) "======================================================================"
    write(00,*) "Time (CPU) = ", real(timeEnd-timeStart,kind=8), "s"                             !当前进程累计消耗的 CPU 时间,包括并行
    write(00,*) "MLUPS = ", real( dble(nx)*dble(ny)*dble(itc)/(timeEnd-timeStart)/1.0d6,kind=8 )   !百万格点更新/秒
    write(00,*) "Time (OMP) = ", real(timeEnd2-timeStart2,kind=8), "s"                           !墙钟时间
    write(00,*) "MLUPS (OMP) = ", real( dble(nx)*dble(ny)*dble(itc)/(timeEnd2-timeStart2)/1.0d6,kind=8 )   !百万格点更新/秒
    write(00,*) "Nu_global =", Nu_global
    write(00,*) "Nu_hot    =", Nu_hot
    write(00,*) "Nu_cold   =", Nu_cold
    write(00,*) "useG =", useG



    write(00,*) "Dellocate Array......"
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


!===========================================================================================================================
!===========================================================================================================================

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
    real(kind=8) :: us2, Bx, By
    character(len=100) :: reloadFileName
    

    itc = 0
    errorU = 100.0d0
    errorT = 100.0d0 


    !-----------------------------------------------------------------------------------------------
    !记录各种信息在日志文件中
    open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')  !在这个txt文件后面继续写（追加模式）
    
    if(outputBinFile.EQ.1) then
        open(unit=01,file=trim(binFolderPrefix)//"-"//"readme",status="unknown")    !trim去掉字符串尾部空格，换了存储路径，可自己更改
        write(01,*) "binFile folder exist!"                                         !读取路径binFolderPrefix="../binFile/buoyancyCavity
        close(01)
        write(00,*) "Data will be stored in ", binFolderPrefix
    endif
    if(outputPltFile.EQ.1) then
        open(unit=01,file=trim(pltFolderPrefix)//"-"//"readme",status="unknown")     !读取路径pltFolderPrefix="../pltFile/buoyancyCavity
        write(01,*) "pltFile folder exist!"
        close(01)
        write(00,*) "Data will be stored in ", pltFolderPrefix
    endif
    
    !if( (paraA.GE.1.0d0).OR.(paraA.LE.-4.0d0) ) then                           !只有在[-4,1]才可以，要不然预警
    !    write(00,*) "----------------------------------"
    !    write(00,*) "paraA=", paraA
    !    write(00,*) "Error: condition not meet for the algorithm"
    !    write(00,*) "Ref: Luo2013, CMA"
    !    write(00,*) "Please try to reduce Mach number"
    !    write(00,*) "----------------------------------"
    !    stop
    !endif

    write(00,*)"-------------------------------------------------------------------------------"
    write(00,*) 'Mesh:',nx,ny
    write(00,*) 'Rayleigh=',real(Rayleigh,kind=8), '; Prandtl =',real(Prandtl,kind=8), '; Mach =',real(Mach,kind=8)
    write(00,*) "Length unit: L0 =", real(lengthUnit,kind=8)
    write(00,*) "Time unit: Sqrt(L0/(gBeta*DeltaT)) =", real(timeUnit,kind=8)
    write(00,*) "Velocity unit: Sqrt(gBeta*L0*DeltaT) =", real(velocityUnit,kind=8)
    write(00,*) "   "
    write(00,*) 'tauf=',real(tauf,kind=8)
    write(00,*) 'taug=',real(taug,kind=8)
    write(00,*) "viscosity =",real(viscosity,kind=8), "; diffusivity =",real(diffusivity,kind=8)
    write(00,*) "outputFrequency =", real(outputFrequency,kind=8), "free-fall time units"
    write(00,*) "......................  or ",  int(outputFrequency*timeUnit), "in itc units"
    write(00,*) "backupInterval =", backupInterval, " free-fall time units"
    write(00,*) ".................... or ", int(backupInterval/outputFrequency)*int(outputFrequency*timeUnit), "in itc units"
    if(loadInitField.EQ.1) then
        write(00,*) "reloadDimensionlessTime=", reloadDimensionlessTime
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

    allocate (f(0:8,nx,ny))
    allocate (f_post(0:8,0:nx+1,0:ny+1))
    allocate (g(0:4,nx,ny))
    allocate (g_post(0:4,0:nx+1,0:ny+1))

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
    
    !omegaT(0) = (1.0d0-paraA)/5.0d0
    !do alpha=1,4
    !    omegaT(alpha) = (paraA+4.0d0)/20.0d0
    !enddo

        omegaT(0) = 1.0d0/3.0d0
    do alpha=1,4
        omegaT(alpha) = 1.0d0/6.0d0
    enddo

    if(loadInitField.EQ.0) then                    !在不加载文件的情况下，都是零场为初值
    
        u = 0.0d0
        v = 0.0d0
        T = 0.0d0
        
        write(00,*) "Initial field is set exactly"
        if(reloadDimensionlessTime.NE.0) then        !在不加载文件的情况下，reloadDimensionlessTime必须是零
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
    do j = 1, ny                                   !在不加载文件的情况下，初始化温度是分层的，其实这个有点问题，边界才是精确的Thot和Tcold
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
                    f(alpha,i,j) = rho(i,j)*omega(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)  !D2Q9标准feq
                enddo
                do alpha = 0, 4
                    un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha) 
                    !g(alpha,i,j) = omegaT(alpha)*T(i,j)*(1.0d0+10.0d0/(4.0d0+paraA)*un(alpha))                        !D2Q5线性geq，不含有u^2项，这个是含有修正的geq
                    g(alpha,i,j) = omegaT(alpha)*T(i,j)*(1.0d0+3.0d0*un(alpha))  
                enddo
            enddo
        enddo

        do j = 1,ny
            do i = 1,nx
              Bx = u(i,j) * T(i,j)
              By = v(i,j) * T(i,j)
              Bx_prev(i,j) = Bx
              By_prev(i,j) = By
            enddo
        enddo




    elseif(loadInitField.EQ.1) then                               !在加载文件的情况下，读取路径 reloadFilePrefix="./reloadFile/backupFile
        if(reloadDimensionlessTime.EQ.0) then                     !在加载文件的情况下，reloadDimensionlessTime最好是非零
            write(00,*) "WARNING: since loadInitField.EQ.1, please confirm reloadDimensionlessTime", reloadDimensionlessTime
            stop
        endif
        write(00,*) "Load initial field from previous simulation: ../reloadFile/backupFile- >>>"
        write(reloadFileName, *) reloadbinFileNum                 !换了个变量Name
        reloadFileName = adjustl(reloadFileName)                  !adjustl把字符串左对齐，把前导空格移到字符串末尾
        open(unit=01,file=trim(reloadFilePrefix)//"-"//trim(reloadFileName)//".bin",form="unformatted", &
        access="sequential",status='old')  !unformatted是二进制,sequential：按记录顺序读写
            ! Strict restart files store only f and g; rho, u, v and T are rebuilt after reading.
            write(00,*) "Reloading f and g from file"
            read(01) (((f(alpha,i,j), i=1,nx), j=1,ny), alpha=0,8)      !先 i，再 j，再 alpha
            read(01) (((g(alpha,i,j), i=1,nx), j=1,ny), alpha=0,4)
        close(01)
        call reconstruct_macro_from_fg()
        write(00,*) "Raw data is loaded from the file: backupFile-",trim(reloadFileName),".bin"
    else
        write(00,*) "Error: initial field is not properly set"                                                  !如果 loadInitField 不是 0/1 或逻辑不一致，直接停止
        stop
    endif
    
    write(00,*)"-------------------------------------------------------------------------------"
close(00)
    
#ifdef steadyFlow
    up = 0.0d0
    vp = 0.0d0
    Tp = 0.0d0
#endif
    
    f_post = 0.0d0
    g_post = 0.0d0
        
    binFileNum = 0
    pltFileNum = 0
    dimensionlessTime = 0
    
    NuVolAvg = 0.0d0
    ReVolAvg = 0.0d0
    
    return
  end subroutine initial
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

    !$omp parallel do default(none) shared(f,f_post,rho,u,v,Fx,Fy,T) private(i,j,alpha,s,m,m_post,meq,fSource) 
    do j = 1, ny
        do i = 1, nx

          m(0) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
          m(1) = -4.0d0*f(0,i,j)-f(1,i,j)-f(2,i,j)-f(3,i,j)-f(4,i,j)+2.0d0*(f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j))
          m(2) = 4.0d0*f(0,i,j)-2.0d0*(f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j))+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
          m(3) = f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
          m(4) = -2.0d0*f(1,i,j)+2.0d0*f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
          m(5) = f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)
          m(6) = -2.0d0*f(2,i,j)+2.0d0*f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)
          m(7) = f(1,i,j)-f(2,i,j)+f(3,i,j)-f(4,i,j)
          m(8) = f(5,i,j)-f(6,i,j)+f(7,i,j)-f(8,i,j)

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

          f_post(0,i,j) = m_post(0)/9.0d0-m_post(1)/9.0d0+m_post(2)/9.0d0                                         !这边是乘以M逆
          f_post(1,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+m_post(3)/6.0d0-m_post(4)/6.0d0 &
                    +m_post(7)/4.0d0
          f_post(2,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    +m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)/4.0d0
          f_post(3,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0-m_post(3)/6.0d0+m_post(4)/6.0d0 &
                    +m_post(7)/4.0d0
          f_post(4,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    -m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)/4.0d0
          f_post(5,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)/4.0d0
          f_post(6,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)/4.0d0
          f_post(7,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)/4.0d0
          f_post(8,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)/4.0d0

        enddo
    enddo
    !$omp end parallel do
    
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
    integer(kind=4) :: ip, jp
    integer(kind=4) :: alpha
    
    !$omp parallel do default(none) shared(f,f_post,ex,ey) private(i,j,ip,jp,alpha)
    do j = 1, ny
        do i = 1, nx
            do alpha = 0, 8                        !上游格点索引：fα(i,j) <- f_postα(i-exα, j-eyα)
                ip = i-ex(alpha)                   !边界附近 (ip/jp 可能为 0 或 nx+1/ny+1)，需在 bounceback/周期边界处理中覆盖修正边界分布
                jp = j-ey(alpha)                   !ghost 层在初始化中为 0，保证不会出现未初始化垃圾值
                
                f(alpha,i,j) = f_post(alpha,ip,jp)
            enddo
        enddo
    enddo
    !$omp end parallel do

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
    !$omp parallel do default(none) shared(f, f_post) private(j)
    do j = 1, ny                                                  !速度边界垂直边界周期，直接方向相同，跨边界的入射分布      
        !Left side (i=1)
        f(1,1,j) = f_post(1,nx,j)
        f(5,1,j) = f_post(5,nx,j)
        f(8,1,j) = f_post(8,nx,j)

        !Right side (i=nx)
        f(3,nx,j) = f_post(3,1,j)
        f(6,nx,j) = f_post(6,1,j)
        f(7,nx,j) = f_post(7,1,j)
    enddo
    !$omp end parallel do
#endif

#ifdef VerticalWallsNoslip
    !$omp parallel do default(none) shared(f, f_post) private(j)
    do j = 1, ny                                                 !速度边界垂直边界静止壁无滑移，直接反弹，方向相反
        !Left side (i=1)
        f(1,1,j) = f_post(3,1,j)
        f(5,1,j) = f_post(7,1,j)
        f(8,1,j) = f_post(6,1,j)

        !Right side (i=nx)
        f(3,nx,j) = f_post(1,nx,j)
        f(6,nx,j) = f_post(8,nx,j)
        f(7,nx,j) = f_post(5,nx,j)
    enddo
    !$omp end parallel do
#endif

#ifdef HorizontalWallsNoslip
    !$omp parallel do default(none) shared(f, f_post) private(i)
    do i = 1, nx                                                  !速度边界水平边界无滑移，直接反弹，方向相反
        !Bottom side (j=1)
        f(2,i,1) = f_post(4,i,1)
        f(5,i,1) = f_post(7,i,1)
        f(6,i,1) = f_post(8,i,1)

        !Top side (j=ny)
        f(4,i,ny) = f_post(2,i,ny)
        f(7,i,ny) = f_post(5,i,ny)
        f(8,i,ny) = f_post(6,i,ny)
    enddo
    !$omp end parallel do
#endif

    return
  end subroutine bounceback
!===================================================================================================
! bounceback 结束: 处理流场边界条件，包括无滑移壁面和相关反弹格式。
!===================================================================================================



!===================================================================================================

!===================================================================================================
! 子程序: macro
!===================================================================================================
  subroutine macro()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !$omp parallel do default(none) shared(f,rho,u,v,Fx,Fy) private(i,j)
    do j = 1, ny
        do i = 1, nx
            rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
            u(i,j) = ( f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)+0.5d0*Fx(i,j) )/rho(i,j)     !含力LBM的半步动量修正：rho*u = Σ f e + 0.5*F，对应Guo forcing的二阶定义
            v(i,j) = ( f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)+0.5d0*Fy(i,j) )/rho(i,j)
        enddo
    enddo
    !$omp end parallel do
    
    return
  end subroutine macro
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
    real(kind=8), parameter :: SG = 1.0d0 - 0.5d0*s_j




    !$omp parallel do default(none) shared(g,g_post,u,v,T,Bx_prev,By_prev) private(i,j,alpha,n,neq,q,n_post,Bx,By,dBx,dBy) 
    do j = 1, ny
        do i = 1, nx

            Bx = u(i,j) * T(i,j)
            By = v(i,j) * T(i,j)

            if (useG) then
              dBx = Bx - Bx_prev(i,j)
              dBy = By - By_prev(i,j)
            else
              dBx = 0.0d0
              dBy = 0.0d0
            end if

            Bx_prev(i,j) = Bx
            By_prev(i,j) = By


          n(0) = g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
          n(1) = g(1,i,j)-g(3,i,j)
          n(2) = g(2,i,j)-g(4,i,j)
          n(3) = -4.0d0*g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
          n(4) = g(1,i,j)-g(2,i,j)+g(3,i,j)-g(4,i,j)
        
          neq(0) = T(i,j)
          neq(1) = T(i,j)*u(i,j)
          neq(2) = T(i,j)*v(i,j)
          neq(3) = T(i,j)*(-2.0d0/3.0d0)
          neq(4) = 0.0d0
        
          q(0) = 0.0d0
          q(1) = s_j
          q(2) = s_j
          q(3) = s_e
          q(4) = s_q
        
          
          n_post(0) = n(0)-q(0)*(n(0)-neq(0))
          n_post(1) = n(1)-q(1)*(n(1)-neq(1))+ SG*dBx
          n_post(2) = n(2)-q(2)*(n(2)-neq(2))+ SG*dBy
          n_post(3) = n(3)-q(3)*(n(3)-neq(3))
          n_post(4) = n(4)-q(4)*(n(4)-neq(4))
          
        
          g_post(0,i,j) = 0.2d0*n_post(0)-0.2d0*n_post(3)
          g_post(1,i,j) = 0.2d0*n_post(0)+0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
          g_post(2,i,j) = 0.2d0*n_post(0)+0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4)
          g_post(3,i,j) = 0.2d0*n_post(0)-0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
          g_post(4,i,j) = 0.2d0*n_post(0)-0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4) 
        enddo
    enddo
    !$omp end parallel do
    
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
    integer(kind=4) :: ip, jp
    integer(kind=4) :: alpha
    
    !$omp parallel do default(none) shared(g, g_post, ex, ey) private(i, j, ip, jp, alpha)
    do j = 1, ny
        do i = 1, nx
            do alpha = 0, 4
                ip = i-ex(alpha)
                jp = j-ey(alpha)
                
                g(alpha,i,j) = g_post(alpha,ip,jp)
            enddo
        enddo
    enddo
    !$omp end parallel do

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
    !$omp parallel do default(none) shared(g,g_post) private(j)
    do j = 1, ny
        !Left boundary
        g(1,1,j) = g_post(1,nx,j)

        !Right boundary
        g(3,nx,j) = g_post(3,1,j)
    enddo
    !$omp end parallel do
#endif

#ifdef VerticalWallsConstT
    !$omp parallel do default(none) shared(g,g_post,omegaT) private(j)
    do j = 1, ny
        !Left boundary
        !g(1,1,j) = -g_post(3,1,j)+(4.0d0+paraA)/10.0d0*Thot
        g(1,1,j) = -g_post(3,1,j)+2.0d0*omegaT(3)*Thot

        !Right boundary
        !g(3,nx,j) = -g_post(1,nx,j)+(4.0d0+paraA)/10.0d0*Tcold
        g(3,nx,j) = -g_post(1,nx,j)+2.0d0*omegaT(1)*Tcold
    enddo
    !$omp end parallel do
#endif

#ifdef VerticalWallsAdiabatic
    !$omp parallel do default(none) shared(g,g_post) private(j)
    do j = 1, ny
        !Left boundary
        g(1,1,j) = g_post(3,1,j)

        !Right boundary
        g(3,nx,j) = g_post(1,nx,j)
    enddo
    !$omp end parallel do
#endif

#ifdef HorizontalWallsAdiabatic
    !$omp parallel do default(none) shared(g,g_post) private(i)
    do i = 1, nx
        !Bottom side
        g(2,i,1) = g_post(4,i,1)

        !Top side
        g(4,i,ny) = g_post(2,i,ny)
    enddo
    !$omp end parallel do
#endif

#ifdef HorizontalWallsConstT
    !$omp parallel do default(none) shared(g,g_post,omegaT) private(i)
    do i = 1, nx
        !Bottom side
        !g(2,i,1) = -g_post(4,i,1)+(4.0d0+paraA)/10.0d0*Thot
        g(2,i,1) = -g_post(4,i,1)+2.0d0*omegaT(4)*Thot

        !Top side
        !g(4,i,ny) = -g_post(2,i,ny)+(4.0d0+paraA)/10.0d0*Tcold
        g(4,i,ny) = -g_post(2,i,ny)+2.0d0*omegaT(2)*Tcold
    enddo
    !$omp end parallel do
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

    !$omp parallel do default(none) shared(g, T) private(i,j)
    do j = 1, ny
        do i = 1, nx
            T(i,j) = g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
        enddo
    enddo
    !$omp end parallel do
    
    return
    end subroutine macroT
!===================================================================================================
! macroT 结束: 由温度分布函数恢复宏观温度场 T。
!===================================================================================================



!===================================================================================================
! reconstruct_macro_from_fg: rebuild rho, u, v, T and previous heat flux from reloaded f and g
!===================================================================================================
    subroutine reconstruct_macro_from_fg()
    use commondata
    implicit none
    integer(kind=4) :: i, j, iter
    real(kind=8) :: momx, momy
    logical :: rho_bad

    call macroT()

    rho_bad = .false.
    !$omp parallel do default(none) shared(f,rho,u,v,T,Fx,Fy,Bx_prev,By_prev) private(i,j,iter,momx,momy) &
    !$omp reduction(.or.:rho_bad)
    do j = 1, ny
        do i = 1, nx
            rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
            momx = f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
            momy = f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)

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

            Bx_prev(i,j) = u(i,j)*T(i,j)
            By_prev(i,j) = v(i,j)*T(i,j)
        enddo
    enddo
    !$omp end parallel do

    if (rho_bad) then
        write(*,*) "Warning: non-positive rho found during restart reconstruction."
        stop
    endif

    return
    end subroutine reconstruct_macro_from_fg
!===================================================================================================
! reconstruct_macro_from_fg end: restart state is fully rebuilt from the reloaded distributions
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



    error1 = 0.0d0
    error2 = 0.0d0

    error5 = 0.0d0
    error6 = 0.0d0
    
    !$omp parallel do default(none) shared(u,up,v,vp,T,Tp) private(i,j) reduction(+:error1,error2,error5,error6)
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
    !$omp end parallel do 
    
    errorU = dsqrt(error1)/dsqrt(error2)                 !速度场相对L2误差：||u^n-u^{n-1}||_2 / ||u^n||_2
    errorT = error5/error6                               !温度场相对L1误差：||T^n-T^{n-1}||_1 / ||T^n||_1

  

    call append_convergence_tecplot('convergence.plt', itc, errorU, errorT)

    !write(caseTag,'("nx=",I0,",ny=",I0,",useG=",L1)') nx, ny, useG         !输出收敛曲线的对比
    write(caseTag,'("Ra=",ES10.3E2,",nx=",I0,",ny=",I0,",useG=",L1)') Rayleigh, nx, ny, useG     
    call append_convergence_master_tecplot('convergence_all.plt', caseTag, itc, errorU, errorT)

    write(*,'(I12,1X,ES24.16,1X,ES24.16)') itc, errorU, errorT


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
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4), intent(in)  :: itc
  real(kind=8),    intent(in)  :: errorU, errorT

  integer :: u
  logical, save :: first_write = .true.

  if (first_write) then
    ! 每次程序新运行的第一次调用：直接覆盖旧文件
    open(newunit=u, file=trim(filename), status='replace', action='write', form='formatted')

    write(u,'(A)') 'VARIABLES = "itc" "errorU" "errorT"'
    write(u,'(A)') 'ZONE T="conv", F=POINT'
    write(u,'(I12,1X,ES24.16,1X,ES24.16)') itc, errorU, errorT
    close(u)

    first_write = .false.
  else
    ! 同一次运行的后续调用：追加数据行
    open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
    write(u,'(I12,1X,ES24.16,1X,ES24.16)') itc, errorU, errorT
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
  implicit none
  character(len=*), intent(in) :: filename, zoneName
  integer(kind=4), intent(in)  :: itc
  real(kind=8),    intent(in)  :: errorU, errorT

  logical :: ex
  integer :: u
  logical, save :: zone_started = .false.

  ! 本次运行第一次写：确保文件头存在 + 写入一个新ZONE
  if (.not. zone_started) then
    inquire(file=trim(filename), exist=ex)

    if (.not. ex) then
      open(newunit=u, file=trim(filename), status='replace', action='write', form='formatted')
      write(u,'(A)') 'TITLE = "Convergence comparison"'
      write(u,'(A)') 'VARIABLES = "itc" "errorU" "errorT"'
      close(u)
    end if

    open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
    write(u,'(A)') 'ZONE T="'//trim(zoneName)//'", F=POINT'
    close(u)

    zone_started = .true.
  end if

  ! 追加一个数据点
  open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
  write(u,'(I12,1X,ES24.16,1X,ES24.16)') itc, errorU, errorT
  close(u)
end subroutine append_convergence_master_tecplot
!===================================================================================================
! append_convergence_master_tecplot 结束: 向总收敛文件追加一条带算例标签的记录。
!===================================================================================================





!===================================================================================================
! File output helpers
!===========================================================================================================================
!===================================================================================================
!===================================================================================================
  subroutine output_binary()                                         !输出uvTrho，存储在binFolderPrefix="../binFile/buoyancyCavity-000000001234.bin
    use commondata                                                   !用于后处理快照；重启读入时必须按 u,v,T,rho 顺序读取
    implicit none
    integer(kind=4) :: i, j
    character(len=100) :: filename
    ! This snapshot is for post-processing only; u/v are written after nondimensionalization.
    ! For strict restart, keep using backupData(), which preserves the lattice-state variables.
    
#ifdef steadyFlow
    write(filename,*) itc                                            !steadyFlow：bin/plt文件名通常直接用 itc 来编写（格子时间步长） 
#endif

#ifdef unsteadyFlow
    binFileNum = binFileNum+1                 !unsteadyFlow：文件名用输出序号 binFileNum（每次输出自增，与 tff 间隔由 outputFrequency 控制）
    if(loadInitField.EQ.0) write(filename,*) binFileNum
    if(loadInitField.EQ.1) write(filename,*) binFileNum+reloadbinFileNum
#endif

    !unsteadyFlow 下：按“输出次数”编号，也就是调用一次加一

    filename = adjustl(filename)

    open(unit=03,file=trim(binFolderPrefix)//"-"//trim(filename)//'.bin',form="unformatted",access="sequential")    !二进制
    ! Post-processing snapshot only: write nondimensionalized u/v together with T and rho.
    ! Do not use this file for strict restart; backupData() keeps lattice velocities for that purpose.
    write(03) ((velocityScaleCompare*u(i,j),i=1,nx),j=1,ny)
    write(03) ((velocityScaleCompare*v(i,j),i=1,nx),j=1,ny)
    write(03) ((T(i,j),i=1,nx),j=1,ny)
    write(03) ((rho(i,j), i=1,nx), j=1,ny)
    close(03)

    return
  end subroutine output_binary
!===================================================================================================
! output_binary 结束: 输出 u、v、T、rho 的二进制快照文件。
!===================================================================================================


    

!===================================================================================================
! 子程序: backupData
! 作用: 输出包含 f、g、u、v、T 的重启备份文件。
! 用途: 在运行过程中定期调用，也在程序结束前调用。
!===================================================================================================
  subroutine backupData()                                         !输出fg，存储在当前路径，名字是backupFile-1000.bin
    use commondata                                                !用于重启，包含 f,g；读入时必须先读 f,g 再读 u,v,T（无 rho）
    implicit none
    integer(kind=4) :: i, j, alpha
    character(len=100) :: filename

#ifdef steadyFlow
    write(filename,*) itc                                    !steadyFlow：bin/plt文件名通常直接用 itc 来编写（格子时间步长）                               
#endif

#ifdef unsteadyFlow
    if(loadInitField.EQ.0) write(filename,*) binFileNum     !unsteadyFlow：文件名用输出序号 binFileNum（每次输出自增，与 tff 间隔由 outputFrequency 控制）
    if(loadInitField.EQ.1) write(filename,*) binFileNum+reloadbinFileNum
#endif

    filename = adjustl(filename)

    open(unit=05,file='backupFile-'//trim(filename)//'.bin',form="unformatted",access="sequential")   !二进制
    ! Strict restart snapshots store only f and g; rho/u/v/T are reconstructed after reload.
    write(05) (((f(alpha,i,j), i=1,nx), j=1,ny), alpha=0,8)
    write(05) (((g(alpha,i,j), i=1,nx), j=1,ny), alpha=0,4)
    close(05)
    
    open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
    write(00,*) "Backup  f and g to the file: backupFile-", trim(filename),".bin"
    close(00)
    
    return
  end subroutine backupData
!===================================================================================================
! backupData 结束: 输出包含 f、g、u、v、T 的重启备份文件。
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
    ! Restart files should still come from backupData(), which preserves the lattice-state information.
    integer(kind=4) :: i, j, k
    REAL(kind=4) :: zoneMarker, eohMarker   !Tecplot 二进制格式里用的两个“标记值”（299 和 357），用于告诉 Tecplot：这里开始是 zone 描述 / header 结束。
    character(len=40) :: title              !文件 Title 字符串
    character(len=40) :: V1,V2,V3,V4,V5     !变量名字符串（X,Y,U,V,T）
    integer(kind=4), parameter :: kmax=1    !二维数据也按 3D 的 IJK 写，K=1
    character(len=40) :: zoneName           !zone 名称
    character(len=100) :: filename          !输出文件名字符串
    
    !$acc update self(u,v,T)                !OpenACC指令，忽略

#ifdef steadyFlow
    write(filename,'(i12.12)') itc
#endif

#ifdef unsteadyFlow
    pltFileNum = pltFileNum+1               !和前面的binFileNum一样
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
!===========================================================================================================================
!===================================================================================================
! 子程序: calNuRe
!===================================================================================================
  subroutine calNuRe()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: NuVolAvg_temp    !体平均 Nu
    real(kind=8) :: ReVolAvg_temp    !体平均 Re
    

    if (dimensionlessTime.GE.dimensionlessTimeMax) then
        write(*,*) "Error: dimensionlessTime exceeds dimensionlessTimeMax, please enlarge dimensionlessTimeMax"
        open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
        write(00,*) "Error: dimensionlessTime exceeds dimensionlessTimeMax, please enlarge dimensionlessTimeMax"
        close(00)
        stop
    endif

    dimensionlessTime = dimensionlessTime+1   !每隔 outputFrequency 个自由落体时间调用一次calNuRe
    
    NuVolAvg_temp = 0.0d0    
    !$omp parallel do default(none) shared(v,T) private(i,j) reduction(+:NuVolAvg_temp)
    do j = 1, ny
        do i = 1, nx
            NuVolAvg_temp = NuVolAvg_temp+v(i,j)*T(i,j)     !对流热通量
        enddo
    enddo
    !$omp end parallel do
    NuVolAvg(dimensionlessTime) = NuVolAvg_temp/dble(nx*ny)*lengthUnit/diffusivity+1.0d0    !!体平均 Nusselt 数 = 1 + (常数系数) × 体平均对流热通量

    open(unit=01,file="Nu_VolAvg.dat",status='unknown',position='append')
    write(01,*) real(reloadDimensionlessTime+dimensionlessTime*outputFrequency,kind=8), NuVolAvg(dimensionlessTime)   !以自由落体时间来写入
    close(01)

    ReVolAvg_temp = 0.0d0
    !$omp parallel do default(none) shared(u,v) private(i,j) reduction(+:ReVolAvg_temp)
    do j = 1, ny
        do i = 1, nx 
            ReVolAvg_temp = ReVolAvg_temp+(u(i,j)*u(i,j)+v(i,j)*v(i,j))
        enddo
    enddo
    !$omp end parallel do
    ReVolAvg(dimensionlessTime) = dsqrt(ReVolAvg_temp/dble(nx*ny))*lengthUnit/viscosity    !全域体平均 RMS-Reynolds 数


    open(unit=02,file="Re_VolAvg.dat",status='unknown',position='append')
    write(02,*) real(reloadDimensionlessTime+dimensionlessTime*outputFrequency,kind=8), ReVolAvg(dimensionlessTime)  !!for print purpose only
    close(02)
    write(*,*)  ReVolAvg(dimensionlessTime)
    return
  end subroutine calNuRe
!===================================================================================================
! calNuRe 结束: 计算体平均 Nu 和 Re 的时间历程统计量。
!===================================================================================================


!===================================================================================================
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
    coef   = heatFluxScale

    sum_qx = 0.0d0

    !$omp parallel do default(none) shared(u,T,dx,coef) private(i,j,dT,qx) reduction(+:sum_qx)
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
    !$omp end parallel do

    Nu_global = (sum_qx / dble(nx*ny)) / deltaT

    ! 屏幕输出
    write(*,'(a,1x,es16.8)') "Nu_global =", Nu_global


    ! 同步写入日志
    open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
    write(00,'(a,1x,es16.8)') "Nu_global =", Nu_global
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
    coef   = heatFluxScale

   

    !-----------------------------
    ! (1) 左侧热壁平均 Nu_hot，同时记录 Numax/Numin 及其 y 位置
    sum_hot = 0.0d0
    !$omp parallel do default(none) shared(T,Nu_left,dx,deltaT) private(j,qx_wall) reduction(+:sum_hot)
    do j = 1, ny
      ! 壁面导热通量：qx(x=0,j)
      qx_wall   = 2.0d0 * (Thot - T(1,j) ) / dx
      Nu_left(j)= qx_wall / deltaT
      sum_hot   = sum_hot + Nu_left(j)
    enddo
    !$omp end parallel do
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
    !$omp parallel do default(none) shared(T,dx,deltaT) private(j,qx_wall) reduction(+:sum_cold)
    do j = 1, ny
      qx_wall  = 2.0d0 * (T(nx,j) - Tcold) / dx
      sum_cold = sum_cold + qx_wall/ deltaT
    enddo
    !$omp end parallel do
    Nu_cold = (sum_cold / dble(ny))

    !-----------------------------
    ! (3) 竖直中线 x=1/2 的平均 Nu_middle
    sum_mid = 0.0d0

    if (mod(nx,2) == 1) then
      iMid = (nx + 1)/2

      !$omp parallel do default(none) shared(u,T,iMid,dx,deltaT,coef) private(j) reduction(+:sum_mid)
      do j = 1, ny
        sum_mid = sum_mid + ( coef*u(iMid,j)*(T(iMid,j)-Tref) + (T(iMid-1,j)-T(iMid+1,j))/(2.0d0*dx) ) / deltaT
      enddo
      !$omp end parallel do

    else
      iL = nx/2
      iR = iL + 1

      !$omp parallel do default(none) shared(u,T,iL,iR,dx,deltaT,coef) private(j) reduction(+:sum_mid)
      do j = 1, ny
        sum_mid = sum_mid + (coef*( 0.5d0*( u(iL,j)*(T(iL,j)-Tref) + u(iR,j)*(T(iR,j)-Tref) )) &
        + (T(iL,j)-T(iR,j))/dx )/ deltaT
      enddo
      !$omp end parallel do
    endif

    Nu_middle = (sum_mid / dble(ny))

    !-----------------------------
    ! 输出：屏幕 + 日志
    write(*,'(a,1x,es16.8)') "Nu_hot    =", Nu_hot
    write(*,'(a,1x,es16.8)') "Nu_cold   =", Nu_cold
    write(*,'(a,1x,es16.8)') "Nu_middle =", Nu_middle
    write(*,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_max =", Nu_hot_max, "y_max =", Nu_hot_max_position
    write(*,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_min =", Nu_hot_min, "y_min =", Nu_hot_min_position

    !open(unit=12,file="Nu_wall_avg.txt",status="unknown",position="append")
    !write(12,'(i12,1x,8(1x,es16.8))') itc, Nu_hot, Nu_cold, Nu_middle, Nu_hot_max, Nu_hot_max_position, Nu_hot_min, Nu_hot_min_position
    !close(12)

    open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
    write(00,'(a,1x,es16.8)') "Nu_hot    =", Nu_hot
    write(00,'(a,1x,es16.8)') "Nu_cold   =", Nu_cold
    write(00,'(a,1x,es16.8)') "Nu_middle =", Nu_middle
    write(00,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_max =", Nu_hot_max, "y_max =", Nu_hot_max_position
    write(00,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_min =", Nu_hot_min, "y_min =", Nu_hot_min_position
    close(00)

    return
  end subroutine SideHeatedcalc_Nu_wall_avg
!===================================================================================================
! SideHeatedcalc_Nu_wall_avg 结束: 计算侧壁差温工况下热壁、冷壁和中线的 Nusselt 数及其极值。
!===================================================================================================



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
    write(*,'(A,1X,F12.6,1X,A,1X,F12.6,1X,A,1X,F12.6)') &
         'u_mid_max =', umax_fit*coef, 'at y =', y_fit, 'on x_mid =', xmid
         

    open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
    string = ctime( time() )
    write(00,*) '--- calc_umid_max --- ', string
    write(00,*) 'x_mid =', xmid
    write(00,*) 'u_mid_max =', umax_fit*coef, ' y_pos =', y_fit, ' (grid j0=', j0, ')'
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
    write(*,'(A,1X,F12.6,1X,A,1X,F12.6,1X,A,1X,F12.6)') &
         'v_mid_max =', vmax_fit*coef, 'at x =', x_fit, 'on y_mid =', ymid


    open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
    string = ctime( time() )
    write(00,*) '--- calc_vmid_max --- ', string
    write(00,*) 'y_mid =', ymid
    write(00,*) 'v_mid_max =', vmax_fit*coef, ' x_pos =', x_fit, ' (grid i0=', i0, ')'
    close(00)

    return
  end subroutine SideHeatedcalc_vmid_max
!===================================================================================================
! SideHeatedcalc_vmid_max 结束: 计算侧壁差温工况下中心线竖直速度的最大值及位置。
!===================================================================================================




!===================================================================================================
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
  coef   = heatFluxScale

  sum_qy = 0.0d0

  !$omp parallel do default(none) &
  !$omp shared(v,T,dy,coef) &
  !$omp private(i,j,dTdy,qy) reduction(+:sum_qy)
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
  !$omp end parallel do

  Nu_global = (sum_qy / dble(nx*ny)) / deltaT

  write(*,'(a,1x,es16.8)') "Nu_global =", Nu_global
  open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
  write(00,'(a,1x,es16.8)') "Nu_global =", Nu_global
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
  coef   = heatFluxScale

  !-----------------------------
  ! (1) 底部热壁平均 Nu_hot（不含角点）
  sum_hot = 0.0d0
  !$omp parallel do default(none) shared(T,Nu_bot,dy,deltaT) private(i,qy_wall) reduction(+:sum_hot)
  do i = 1, nx
    qy_wall  = 2.0d0 * (Thot - T(i,1)) / dy
    Nu_bot(i)= qy_wall / deltaT
    sum_hot  = sum_hot + Nu_bot(i)
  enddo
  !$omp end parallel do
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
  !$omp parallel do default(none) shared(T,dy,deltaT) private(i,qy_wall) reduction(+:sum_cold)
  do i = 1, nx
    qy_wall  = 2.0d0 * ( T(i,ny) - Tcold ) / dy
    sum_cold = sum_cold + qy_wall / deltaT
  enddo
  !$omp end parallel do
  Nu_cold = sum_cold / dble(nx)

  !-----------------------------
  ! 输出：屏幕 + 日志
  sum_mid = 0.0d0

  if (mod(ny,2) == 1) then
    jMid = (ny + 1)/2

    !$omp parallel do default(none) shared(v,T,jMid,dy,deltaT,coef) private(i) reduction(+:sum_mid)
    do i = 1, nx
      sum_mid = sum_mid + ( coef*v(i,jMid)*(T(i,jMid)-Tref) - (T(i,jMid+1)-T(i,jMid-1))/(2.0d0*dy) ) / deltaT
    enddo
    !$omp end parallel do

  else
    jB = ny/2
    jT = jB + 1

    !$omp parallel do default(none) shared(v,T,jB,jT,dy,deltaT,coef) private(i) reduction(+:sum_mid)
    do i = 1, nx
      sum_mid = sum_mid + (coef*( 0.5d0*( v(i,jB)*(T(i,jB)-Tref) + v(i,jT)*(T(i,jT)-Tref) )) &
      + (T(i,jB)-T(i,jT))/dy ) / deltaT
    enddo
    !$omp end parallel do
  endif

  Nu_middle = sum_mid / dble(nx)

  !-----------------------------
  ! 输出：屏幕 + 日志
  write(*,'(a,1x,es16.8)') "Nu_hot(bottom) =", Nu_hot
  write(*,'(a,1x,es16.8)') "Nu_cold(top)   =", Nu_cold
  write(*,'(a,1x,es16.8)') "Nu_middle      =", Nu_middle
  write(*,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_max =", Nu_hot_max, "x_max =", Nu_hot_max_position
  write(*,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_min =", Nu_hot_min, "x_min =", Nu_hot_min_position

  open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
  write(00,'(a,1x,es16.8)') "Nu_hot(bottom) =", Nu_hot
  write(00,'(a,1x,es16.8)') "Nu_cold(top)   =", Nu_cold
  write(00,'(a,1x,es16.8)') "Nu_middle      =", Nu_middle
  write(00,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_max =", Nu_hot_max, "x_max =", Nu_hot_max_position
  write(00,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_min =", Nu_hot_min, "x_min =", Nu_hot_min_position
  close(00)

  return
end subroutine RBcalc_Nu_wall_avg
!===================================================================================================
! RBcalc_Nu_wall_avg 结束: 计算 Rayleigh-Benard 工况下热壁、冷壁和中线的 Nusselt 数及其极值。
!===================================================================================================


!===================================================================================================
! 子程序: RBcalc_umid_max
! 作用: 计算 Rayleigh-Benard 工况下中心线速度极值及位置。
! 用途: 在 RayleighBenardCell 工况结束后的后处理中调用。
!===================================================================================================
subroutine RBcalc_umid_max()
  use commondata
  implicit none
  integer(kind=4) :: i, k
  integer(kind=4) :: jMid, jB, jT
  integer(kind=4) :: i0
  integer(kind=4) :: ii(5)
  real(kind=8) :: uline(1:nx)
  real(kind=8) :: xk(5), fk(5)
  real(kind=8) :: umax_fit, x_fit
  real(kind=8) :: coef, ymid
  real(kind=8) :: umax_grid

  coef = velocityScaleCompare

  ! ---- (1) 取 y=1/2 中线剖面 u(x_i, y=1/2) ----
  if (mod(ny,2) == 1) then
    jMid = (ny + 1)/2
    ymid = yp(jMid)
    do i = 1, nx
      uline(i) = u(i,jMid)
    enddo
  else
    jB = ny/2
    jT = jB + 1
    ymid = 0.5d0*(yp(jB) + yp(jT))
    do i = 1, nx
      uline(i) = 0.5d0*(u(i,jB) + u(i,jT))
    enddo
  endif

  ! ---- (2) 找峰值所在网格点 i0----------------
  i0 = 1
  umax_grid = uline(1)
  do i = 2, nx
    if (uline(i) > umax_grid) then
      umax_grid = uline(i)
      i0 = i
    endif
  enddo

  ! ---- (3) 取 5 点（尽量对称；靠近端点时偏侧）----
  if (i0 <= 2) then
    ii = (/ 1, 2, 3, 4, 5 /)
  elseif (i0 >= nx-1) then
    ii = (/ nx-4, nx-3, nx-2, nx-1, nx /)
  else
    ii = (/ i0-2, i0-1, i0, i0+1, i0+2 /)
  endif

  do k = 1, 5
    xk(k) = xp(ii(k))
    fk(k) = uline(ii(k))   
  enddo

  call fit_parabola_ls5(xk, fk, +1, umax_fit, x_fit)
  

  write(*,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'u_mid_max* =', umax_fit*coef, 'x =', x_fit, 'y_mid =', ymid

  open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
  write(00,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'u_mid_max* =', umax_fit*coef, 'x =', x_fit, 'y_mid =', ymid
  close(00)

  return
end subroutine RBcalc_umid_max
!===================================================================================================
! RBcalc_umid_max 结束: 计算 Rayleigh-Benard 工况下中心线速度极值及位置。
!===================================================================================================


!===================================================================================================
! 子程序: RBcalc_vmid_max
! 作用: 计算 Rayleigh-Benard 工况下另一方向中心线速度极值及位置。
! 用途: 在 RayleighBenardCell 工况结束后的后处理中调用。
!===================================================================================================
subroutine RBcalc_vmid_max()
  use commondata
  implicit none
  integer(kind=4) :: j, k
  integer(kind=4) :: iMid, iL, iR
  integer(kind=4) :: j0
  integer(kind=4) :: jj(5)
  real(kind=8) :: vline(1:ny)
  real(kind=8) :: yk(5), fk(5)
  real(kind=8) :: vmax_fit, y_fit
  real(kind=8) :: coef, xmid
  real(kind=8) :: vmax_grid

  coef = velocityScaleCompare

  ! ---- (1) 取 x=1/2 中线剖面 v(x=1/2, y_j) ----
  if (mod(nx,2) == 1) then
    iMid = (nx + 1)/2
    xmid = xp(iMid)
    do j = 1, ny
      vline(j) = v(iMid,j)
    enddo
  else
    iL = nx/2
    iR = iL + 1
    xmid = 0.5d0*(xp(iL) + xp(iR))
    do j = 1, ny
      vline(j) = 0.5d0*(v(iL,j) + v(iR,j))
    enddo
  endif

  ! ---- (2) 找峰值所在网格点 j0 ----
  j0 = 1
  vmax_grid = vline(1)
  do j = 2, ny
    if (vline(j) > vmax_grid) then
      vmax_grid = vline(j)
      j0 = j
    endif
  enddo


  ! ---- (3) 取 5 点并拟合 ----
  if (j0 <= 2) then
    jj = (/ 1, 2, 3, 4, 5 /)
  elseif (j0 >= ny-1) then
    jj = (/ ny-4, ny-3, ny-2, ny-1, ny /)
  else
    jj = (/ j0-2, j0-1, j0, j0+1, j0+2 /)
  endif

  do k = 1, 5
    yk(k) = yp(jj(k))
    fk(k) = vline(jj(k))
  enddo

  call fit_parabola_ls5(yk, fk, +1, vmax_fit, y_fit)

  write(*,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'v_mid_max* =', vmax_fit*coef, 'y =', y_fit, 'x_mid =', xmid

  open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
  write(00,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'v_mid_max* =', vmax_fit*coef, 'y =', y_fit, 'x_mid =', xmid
  close(00)

  return
end subroutine RBcalc_vmid_max
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

  write(*,'(a,1x,es16.8)') "abs(psi_center_fine) =", psi_center_abs_fine

  write(*,'(a,1x,es16.8,2x,a,1x,es16.8,2x,a,1x,es16.8)') &
       "max(|psi|) =", psi_abs_max, "x* =", x_at_max, "y* =", y_at_max

  open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
  write(00,'(a,1x,es16.8)') "abs(psi_center_fine) =", psi_center_abs_fine
  write(00,'(a,1x,es16.8,2x,a,1x,es16.8,2x,a,1x,es16.8)') &
       "max(|psi|) =", psi_abs_max, "x* =", x_at_max, "y* =", y_at_max
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
  write(filename,'(i12.12)') itc
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
        write(41) xp(i)
        write(41) yp(j)
        write(41) psi(i,j)
        write(41) vort(i,j)
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
  write(*,'(a,1x,es16.8)') "abs(psi_center_coarse) =", psi_center_abs

  ! Log output
  open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
  write(00,'(a,1x,es16.8)') "abs(psi_center_coarse) =", psi_center_abs
  close(00)

  return
contains
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

