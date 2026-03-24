!=============================================================
!!!    娉ㄩ噴鍖猴紝浠ｇ爜鎻忚堪
!!!    娴姏椹卞姩鑷劧瀵规祦锛堜笂涓嬪姞鐑級
!!!    LBM鏂规硶
!!!    MRT-LBE
!=============================================================

!   鑷畾涔夊畯锛屼竴浜涢€夐」鐨勫紑鍏?
#define steadyFlow    
!~!!#define unsteadyFlow

!   閫熷害杈圭晫锛屽寘鎷按骞冲瀭鐩磋竟鐣屾棤婊戠Щ锛岃繕鏈夊瀭鐩磋竟鐣岄€熷害鍛ㄦ湡
#define HorizontalWallsNoslip
#define VerticalWallsNoslip
!~#define VerticalWallsPeriodicalU

!   娓╁害杈圭晫(for Rayleigh Benard Cell)锛屽寘鎷按骞宠竟鐣屾亽娓╋紝鍨傜洿杈圭晫娓╁害涓嶅彲绌块€忎互鍙婂懆鏈?
!#define RayleighBenardCell
!#define HorizontalWallsConstT
!#define VerticalWallsAdiabatic
!~#define VerticalWallsPeriodicalT


!   娓╁害杈圭晫(for Side Heated Cell)锛屽寘鎷按骞宠竟鐣屾俯搴︿笉鍙┛閫忥紝鍨傜洿杈圭晫鎭掓俯,渚у鍔犵儹鍔犵鍦?
#define SideHeatedCell
#define HorizontalWallsAdiabatic
#define VerticalWallsConstT
!#define SideHeatedHa          
!~~temperature B.C.~~

!   鍏ㄥ眬妯″潡
    module commondata
        implicit none
        !===============================================================================================
        ! 鏄惁鍦ㄨ绠楀墠浠庢棫绠椾緥閲嶅惎
        integer(kind=4), parameter :: loadInitField=0   ! 0: 涓嶉噸鍚紱1: 浠?backupFile-*.bin 璇诲彇鍒濆€?

        ! 鍦?loadInitField=1 鐨勫墠鎻愪笅锛?
        integer(kind=4), parameter :: loadInitFG=0      ! 0: 涓嶈 f,g锛堢敤璇诲叆鐨?u,v,T 閲嶆柊鏋勯€犲钩琛″垎甯冿級锛?: 璇诲叆 f,g锛堜弗鏍奸噸鍚級
        integer(kind=4), parameter :: loadInitRho=0     ! 0: 涓嶈 rho锛堝厛璁句负1锛涜嫢 loadInitFG=1 鍒欏悗缁?rho=危f 浼氳鐩栵級锛?: 浠庢枃浠惰鍏?rho
        integer(kind=4), parameter :: reloadDimensionlessTime=0  ! 鏃х畻渚嬬疮璁＄殑鏃犻噺绾叉椂闂达紙鐢ㄤ簬缁啓 Nu/Re 杈撳嚭妯潗鏍囷級
        integer(kind=4), parameter :: reloadbinFileNum=0         ! 璇诲彇鐨勫浠芥枃浠剁紪鍙凤細backupFile-<reloadbinFileNum>.bin
        !===============================================================================================

        !-----------------------------------------------------------------------------------------------
#if 0
        ! 鏃犻噺绾插弬鏁?
        integer(kind=4), parameter :: nx=22, ny=22     !鏍煎瓙缃戞牸
        real(kind=8), parameter :: lengthUnit=dble(ny)     !鏃犻噺绾查暱搴?
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
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1


#ifdef  SideHeatedHa
        real(kind=8), parameter :: Ha=20.0d0                           !纾佸満寮哄害
        real(kind=8), parameter :: phi=(0.0d0)*(pi/180.0d0)            !纾佸満瑙掑害锛屼互姘村钩鍚戝彸涓?锛屼慨鏀?.0d0鍗冲彲
        real(kind=8), parameter :: B2sigemarho=(Ha**2*viscosity)/(lengthUnit*lengthUnit)  !鍔ㄩ噺鏂圭▼涓婄殑婧愰」绯绘暟
#endif

        ! 楂橀樁鐭╁弬鏁颁慨姝?
        real(kind=8), parameter :: paraA=20.0d0*dsqrt(3.0d0)*diffusivity-4.0d0



        ! 娴姏椤瑰弬鏁?
        real(kind=8), parameter :: gBeta1=Rayleigh*viscosity*diffusivity/lengthUnit
        real(kind=8), parameter :: gBeta=gBeta1/lengthUnit/lengthUnit             !gbeta
        
        real(kind=8), parameter :: timeUnit=dsqrt(lengthUnit/gBeta)      !鏃犻噺绾叉椂闂?
        real(kind=8), parameter :: velocityUnit=dsqrt(gBeta*lengthUnit)  !鏃犻噺绾查€熷害
    
        real(kind=8), parameter :: Snu=1.0d0/tauf, Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)  !鍔ㄩ噺鐨勫鏉惧紱绯绘暟
        !real(kind=8), parameter :: Qd=3.0d0-dsqrt(3.0d0), Qnu=4.0d0*dsqrt(3.0d0)-6.0d0             !娓╁害鐨勫鏉惧紱绯绘暟
        real(kind=8), parameter :: taug = 0.5d0 + (tauf - 0.5d0)/Prandtl
        real(kind=8), parameter :: s_e = 1.0d0, s_q = 1.0d0, s_j = 1.0d0/taug
#endif
        integer(kind=4), parameter :: nx=41, ny=41
        integer(kind=4), parameter :: meshModeUniform=0, meshModeErf=1
        integer(kind=4), parameter :: meshMode=meshModeErf
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1
        real(kind=8), parameter :: stretchA=1.5d0
        real(kind=8), parameter :: pi = acos(-1.0d0)

        real(kind=8), parameter :: Rayleigh=1.0d6
        real(kind=8), parameter :: Prandtl=0.71d0
        real(kind=8), parameter :: Mach=0.1d0
        real(kind=8), parameter :: Thot=0.5d0, Tcold=-0.5d0
        real(kind=8), parameter :: Tref=0.5d0*(Thot+Tcold)

        real(kind=8) :: lengthUnit
        real(kind=8) :: tauf, viscosity, diffusivity
        real(kind=8) :: paraA
        real(kind=8) :: gBeta1, gBeta
        real(kind=8) :: timeUnit, velocityUnit
        real(kind=8) :: Snu, Sq
        real(kind=8) :: taug, s_e, s_q, s_j
        real(kind=8) :: uniformShift
#ifdef  SideHeatedHa
        real(kind=8), parameter :: Ha=20.0d0
        real(kind=8), parameter :: phi=(0.0d0)*(pi/180.0d0)
        real(kind=8) :: B2sigemarho
#endif
        !-----------------------------------------------------------------------------------------------           
        
        !===============================================================================================
        ! 杈撳嚭/澶囦唤鐩稿叧璁剧疆锛堜互鑷敱钀戒綋鏃堕棿 t_ff 涓哄崟浣嶏級
        real(kind=8), parameter :: outputFrequency=100.0d0   ! 姣忛殧 outputFrequency 涓嚜鐢辫惤浣撴椂闂?t_ff 杈撳嚭/缁熻涓€娆★紙鏃堕棿闂撮殧锛?

        integer(kind=4), parameter :: dimensionlessTimeMax=int(12000/outputFrequency)  ! 鐢ㄤ簬 NuVolAvg/ReVolAvg 鏁扮粍鐨勬渶澶ц褰曠偣鏁帮紱姣?outputFrequency*t_ff 璁板綍涓€娆★紝鍒欐渶澶氳鐩栫害 12000*t_ff

        integer(kind=4), parameter :: backupInterval=1000  ! 澶囦唤闂撮殧锛堣嚜鐢辫惤浣撴椂闂村崟浣嶏級锛屼负浜嗗仠鐢垫儏鍐典笅锛屽彲浠ョ户缁绠?
        
        real(kind=8), parameter :: epsU=1.0d-9, epsT=1.0d-9    ! paper Eq. (30) convergence thresholds for steady side-heated runs

        integer(kind=4), parameter :: outputBinFile=0   ! 鏄惁杈撳嚭 bin 鏂囦欢锛?=涓嶈緭鍑猴紝1=杈撳嚭
        integer(kind=4), parameter :: outputPltFile=0   ! 鏄惁杈撳嚭 plt 鏂囦欢锛?=涓嶈緭鍑猴紝1=杈撳嚭

                integer(kind=4) :: binFileNum, pltFileNum  ! bin/plt 杈撳嚭鏂囦欢鐨勮鏁板櫒
        ! - unsteadyFlow锛氭瘡娆¤緭鍑洪€掑锛堢敤浜庢枃浠跺悕缂栧彿锛?
        ! - steadyFlow锛歜in/plt 鏂囦欢鍚嶉€氬父鐩存帴鐢?itc

        integer(kind=4) :: dimensionlessTime
        ! 缁熻/杈撳嚭鏃堕棿鐐圭紪鍙凤紙涓?outputFrequency 瀵瑰簲锛夛細
        ! 姣忚皟鐢ㄤ竴娆?calNuRe() 灏?dimensionlessTime = dimensionlessTime + 1
        ! 鐢ㄤ簬绱㈠紩 NuVolAvg/ReVolAvg 鏁扮粍锛屽苟鐢ㄤ簬杈撳嚭鐨勬椂闂磋酱锛歵 = reloadDimensionlessTime + dimensionlessTime*outputFrequency锛堝崟浣嶏細t_ff锛?

        real(kind=8) :: NuVolAvg(0:dimensionlessTimeMax), ReVolAvg(0:dimensionlessTimeMax)
        ! 浣撳钩鍧?Nu 鍜?Re 鐨勬椂闂村簭鍒楃紦瀛?
        ! 鍙湁鍦ㄥ惎鐢ㄥ苟璋冪敤 calNuRe() 鐨勬儏鍐典笅杩欎簺鏁扮粍鎵嶄細琚湡姝ｅ～鍏?
  
        character(len=100) :: binFolderPrefix="./binFile/buoyancyCavity"
        ! bin 杈撳嚭鏂囦欢鍓嶇紑锛堝疄闄呮枃浠跺悕褰㈠锛?binFolderPrefix>-<缂栧彿>.bin锛?

        character(len=100) :: pltFolderPrefix="./pltFile/buoyancyCavity"
        ! plt 杈撳嚭鏂囦欢鍓嶇紑锛堝疄闄呮枃浠跺悕褰㈠锛?pltFolderPrefix>-<缂栧彿>.plt锛?

        character(len=100) :: reloadFilePrefix="./reloadFile/backupFile"
        ! 閲嶅惎璇诲彇鏂囦欢鐨勫墠缂€锛堝疄闄呰鍙栵細<reloadFilePrefix>-<reloadbinFileNum>.bin锛?
        !===============================================================================================

        !-----------------------------------------------------------------------------------------------
        !璁＄畻涓渶瑕佺殑鐩稿叧鍙傛暟
        real(kind=8) :: errorU, errorT
        
        real(kind=8) :: xp(0:nx+1), yp(0:ny+1)      !鏃犻噺绾茬殑鍧愭爣鏁扮粍锛屽寘鎷竟鐣?
        real(kind=8) :: wx(1:nx), wy(1:ny)
        real(kind=8) :: Lx_eff, Ly_eff
        real(kind=8) :: quadSumX, quadSumY, quadSumArea
        real(kind=8) :: maxStretchRatioX, maxStretchRatioY, maxAspectRatio
        real(kind=8) :: symmetryErrorX, symmetryErrorY
        integer(kind=4) :: stream_ix(0:8,1:nx,3), stream_iy(0:8,1:ny,3)
        real(kind=8) :: stream_wx(0:8,1:nx,3), stream_wy(0:8,1:ny,3)
        logical :: stream_x_valid(0:8,1:nx), stream_y_valid(0:8,1:ny)
        real(kind=8), allocatable :: u(:,:), v(:,:), T(:,:), rho(:,:)
#ifdef steadyFlow
        real(kind=8), allocatable :: up(:,:), vp(:,:), Tp(:,:)   !瀛樺偍涔嬪墠鐨勬暟鎹紝鐢ㄦ潵绠楁敹鏁涘垽鎹?
#endif
        real(kind=8), allocatable :: f(:,:,:), f_post(:,:,:)
        real(kind=8), allocatable :: g(:,:,:), g_post(:,:,:)
        real(kind=8), allocatable :: Fx(:,:), Fy(:,:)
        real(kind=8), allocatable :: Bx_prev(:,:), By_prev(:,:)
        integer(kind=4) :: itc
        integer(kind=4), parameter :: itc_max=20000000 !鏍煎瓙鏃堕棿姝ラ暱
        logical, parameter :: useG = .true.            !M1G 寮€鍏?
        real(kind=8) :: Nu_global, Nu_hot, Nu_cold, Nu_middle    !骞冲潎Nu锛屽叏鍦猴紝渚у浠ュ強涓嚎
        real(kind=8) :: Nu_hot_max, Nu_hot_min, Nu_hot_max_position, Nu_hot_min_position    !宸︿晶澹侀潰鐨勬渶澶ф渶灏廚u锛屼互鍙婂搴旂殑浣嶇疆
        
        
        !鏍煎瓙绂绘暎閫熷害鍜屾潈閲?
        integer(kind=4) :: ex(0:8), ey(0:8)
        data ex/0, 1, 0, -1,  0, 1, -1, -1,  1/
        data ey/0, 0, 1,  0, -1, 1,  1, -1, -1/
        real(kind=8) :: omega(0:8), omegaT(0:4)
        !-----------------------------------------------------------------------------------------------

    end module commondata



    program main

    use omp_lib
    use commondata
    implicit none
    real(kind=8) :: timeStart, timeEnd
    character(len=24) :: ctime, string
    INTEGER(kind=4) :: time
    real(kind=8) :: timeStart2, timeEnd2
    integer(kind=4) :: myMaxThreads
    

    !-----------------------------------------------------------------------------------------------
    !璁剧疆骞惰鏍告暟
    open(unit=00,file="SimulationSettings.txt",status='unknown')   !鎵撳紑锛堟垨鍒涘缓锛塼xt鏂囦欢锛屽噯澶囧啓鍏?
    string = ctime( time() )                      !ctime鎶?time() 杩斿洖鐨勬椂闂存埑杞崲鎴愬彲璇荤殑瀛楃涓?
    write(00,*) 'Start: ', string                 !浠€涔堟椂鍊欏紑濮嬭绠?
    write(00,*) "Starting OpenMP >>>>>>"
    call OMP_set_num_threads(24)                   !浣跨敤 24 涓嚎绋?
    myMaxThreads = OMP_get_max_threads()           !鏌ヨ鏈€澶у彲鐢ㄧ嚎绋嬫暟
    write(00,*) "Max Running threads:",myMaxThreads
    close(00)
    !-----------------------------------------------------------------------------------------------


    call initial()

    call CPU_TIME(timeStart)         !褰撳墠杩涚▼绱娑堣€楃殑 CPU 鏃堕棿,鍖呮嫭骞惰
    timeStart2 = OMP_get_wtime()     !澧欓挓鏃堕棿(瀹為檯鑰楁椂锛屼笉鍖呮嫭骞惰)
    do while( ((errorU.GT.epsU).OR.(errorT.GT.epsT)).AND.(itc.LE.itc_max) )   !鍙 (errorU > epsU 鎴?errorT > epsT) 涓?itc 鈮?itc_max锛屽氨缁х画寰幆
                                                                              !鎹㈡垚if锛屽氨鏄?errorU > epsU .and. errorT > epsT

        itc = itc+1
        
        call collision()

        call streaming_interp()

        call bounceback()

        call macro()

        call collisionT()

        call streamingT_interp()

        call bouncebackT()
        
        call macroT()

#ifdef steadyFlow
        if(MOD(itc,2000).EQ.0) call check()
#endif
        
         if( MOD(itc, int(outputFrequency*timeUnit)).EQ.0 ) then  ! 杈惧埌涓€涓緭鍑洪棿闅攐utputFrequency锛屽氨鎵ц涓€娆¤绠椾綋骞冲潎Nu,鏄娴侀€氶噺锛屼笉鏄叏鍦哄钩鍧嘚u

            call calNuRe_islbm()
            
#ifdef steadyFlow
             !if( (outputPltFile.EQ.1).AND.(MOD(itc, backupInterval*int(timeUnit)).EQ.0) ) call output_Tecplot()  !plt 杈撳嚭澶囦唤闂撮殧uvT
#endif
            
#ifdef unsteadyFlow
             if(outputBinFile.EQ.1) then
                     call output_binary()          !bin 杈撳嚭澶囦唤闂撮殧uvTrho
                     if(MOD(itc, int(backupInterval/outputFrequency)*int(outputFrequency*timeUnit)).EQ.0) call backupData()  !杈撳嚭澶囦唤闂撮殧uvTfg
             endif   
             if(outputPltFile.EQ.1) call output_Tecplot()  !plt 杈撳嚭澶囦唤闂撮殧uvT
#endif
        endif
     enddo





    call CPU_TIME(timeEnd)         !褰撳墠杩涚▼绱娑堣€楃殑 CPU 鏃堕棿,鍖呮嫭骞惰
    timeEnd2 = OMP_get_wtime()     !澧欓挓鏃堕棿(瀹為檯鑰楁椂锛屼笉鍖呮嫭骞惰)


    call output_Tecplot() 


!渚у鍔犵儹鍜孯B瀵规祦鐨勮绠桸u涓嶄竴鏍?
#ifdef SideHeatedCell                        
    call SideHeatedcalc_Nu_global()          ! 鍏ㄥ満骞冲潎Nu
    call SideHeatedcalc_Nu_wall_avg()  ! 鐑?鍐峰, 涓嚎骞冲潎Nu,浠ュ強鐑鏈€澶umax鍜孨umin浠ュ強浣嶇疆锛岄兘閲囩敤浜旂偣鏈€灏忎簩涔樻硶鎻掑€煎嚭鏉?
    
    call SideHeatedcalc_umid_max()     !涓績绾夸笂鐨勬渶澶ч€熷害鍙婂叾浣嶇疆锛屼篃鏄敤浜旂偣鏈€灏忎簩涔樻硶鎻掑€煎嚭鏉?
    call SideHeatedcalc_vmid_max()
#endif

#ifdef RayleighBenardCell
    call RBcalc_Nu_global()          ! 鍏ㄥ満骞冲潎Nu
    call RBcalc_Nu_wall_avg()  ! 鐑?鍐峰, 涓嚎骞冲潎Nu,浠ュ強鐑鏈€澶umax鍜孨umin浠ュ強浣嶇疆锛岄兘閲囩敤浜旂偣鏈€灏忎簩涔樻硶鎻掑€煎嚭鏉?
    
    call RBcalc_umid_max()     !涓績绾夸笂鐨勬渶澶ч€熷害鍙婂叾浣嶇疆锛屼篃鏄敤浜旂偣鏈€灏忎簩涔樻硶鎻掑€煎嚭鏉?
    call RBcalc_vmid_max()
#endif


    call calc_psi_vort_and_output()  !杈撳嚭鑵斾綋涓績鐨刟bs(psi),浠ュ強鏈€澶х殑abs(psi)max浠ュ強浣嶇疆锛堥噰鐢ㄧ粏缃戞牸鎻掑€煎嚭鏉ワ級





    if(outputBinFile.EQ.1) call backupData()       !杈撳嚭澶囦唤闂撮殧uvTfg
    
    if(outputPltFile.EQ.1) call output_Tecplot()   !plt 杈撳嚭澶囦唤闂撮殧uvT
     

    open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')        !鍦ㄨ繖涓猼xt鏂囦欢鍚庨潰缁х画鍐欙紙杩藉姞妯″紡锛?
    write(00,*) "Time (CPU) = ", real(timeEnd-timeStart,kind=8), "s"                             !褰撳墠杩涚▼绱娑堣€楃殑 CPU 鏃堕棿,鍖呮嫭骞惰
    write(00,*) "MLUPS = ", real( dble(nx)*dble(ny)*dble(itc)/(timeEnd-timeStart)/1.0d6,kind=8 )   !鐧句竾鏍肩偣鏇存柊/绉?    write(00,*) "Time (OMP) = ", real(timeEnd2-timeStart2,kind=8), "s"                           !澧欓挓鏃堕棿
    write(00,*) "MLUPS (OMP) = ", real( dble(nx)*dble(ny)*dble(itc)/(timeEnd2-timeStart2)/1.0d6,kind=8 )   !鐧句竾鏍肩偣鏇存柊/绉?    write(00,*) "Nu_global =", Nu_global
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
    write(00,*) 'End:   ', string           !浠€涔堟椂鍊欑畻瀹?
    close(00)

    
    end program main


  subroutine setup_mesh_and_runtime()
    use commondata
    implicit none

    call build_mesh_coordinates()
    call setup_runtime_parameters()
    call build_quadrature_weights()
    call prepare_streaming_stencils()
    call interpolation_self_check()

    return
  end subroutine setup_mesh_and_runtime


  subroutine build_mesh_coordinates()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: erfNorm

    xp(0) = 0.0d0
    xp(nx+1) = 1.0d0
    yp(0) = 0.0d0
    yp(ny+1) = 1.0d0

    if (meshMode == meshModeErf .and. dabs(stretchA) > 1.0d-14) then
      erfNorm = erf(0.5d0*stretchA)
      do i = 1, nx
        xp(i) = 0.5d0 * ( 1.0d0 + &
          erf(stretchA*(dble(i)/dble(nx+1) - 0.5d0)) / erfNorm )
      enddo
      do j = 1, ny
        yp(j) = 0.5d0 * ( 1.0d0 + &
          erf(stretchA*(dble(j)/dble(ny+1) - 0.5d0)) / erfNorm )
      enddo
    else
      do i = 1, nx
        xp(i) = (dble(i) - 0.5d0) / dble(nx)
      enddo
      do j = 1, ny
        yp(j) = (dble(j) - 0.5d0) / dble(ny)
      enddo
    endif

    return
  end subroutine build_mesh_coordinates


  subroutine setup_runtime_parameters()
    use commondata
    implicit none

    if (meshMode == meshModeErf) then
      lengthUnit = 1.0d0 / yp(1) - 1.0d0
    else
      lengthUnit = dble(ny)
    endif

    uniformShift = 1.0d0 / lengthUnit
    tauf = 0.5d0 + Mach*lengthUnit*dsqrt(3.0d0*Prandtl/Rayleigh)
    viscosity = (tauf - 0.5d0) / 3.0d0
    diffusivity = viscosity / Prandtl
    paraA = 20.0d0*dsqrt(3.0d0)*diffusivity - 4.0d0
    gBeta1 = Rayleigh*viscosity*diffusivity/lengthUnit
    gBeta = gBeta1/(lengthUnit*lengthUnit)
    timeUnit = dsqrt(lengthUnit/gBeta)
    velocityUnit = dsqrt(gBeta*lengthUnit)
    Snu = 1.0d0 / tauf
    Sq = 8.0d0*(2.0d0*tauf - 1.0d0)/(8.0d0*tauf - 1.0d0)
    taug = 0.5d0 + (tauf - 0.5d0)/Prandtl
    s_e = 1.0d0
    s_q = 1.0d0
    s_j = 1.0d0 / taug
#ifdef SideHeatedHa
    B2sigemarho = (Ha**2*viscosity)/(lengthUnit*lengthUnit)
#endif

    return
  end subroutine setup_runtime_parameters


  subroutine build_quadrature_weights()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: ratio, dxLoc, dyLoc

    do i = 1, nx
      wx(i) = 0.5d0 * (xp(i+1) - xp(i-1))
      if (wx(i) <= 0.0d0) then
        write(*,*) "Error: non-monotonic x mesh"
        stop
      endif
    enddo
    do j = 1, ny
      wy(j) = 0.5d0 * (yp(j+1) - yp(j-1))
      if (wy(j) <= 0.0d0) then
        write(*,*) "Error: non-monotonic y mesh"
        stop
      endif
    enddo

    quadSumX = sum(wx)
    quadSumY = sum(wy)
    quadSumArea = quadSumX * quadSumY
    Lx_eff = 1.0d0 - xp(1)
    Ly_eff = 1.0d0 - yp(1)

    maxStretchRatioX = 1.0d0
    do i = 1, nx-1
      ratio = wx(i+1) / wx(i)
      if (ratio < 1.0d0) ratio = 1.0d0 / ratio
      if (ratio > maxStretchRatioX) maxStretchRatioX = ratio
    enddo

    maxStretchRatioY = 1.0d0
    do j = 1, ny-1
      ratio = wy(j+1) / wy(j)
      if (ratio < 1.0d0) ratio = 1.0d0 / ratio
      if (ratio > maxStretchRatioY) maxStretchRatioY = ratio
    enddo

    symmetryErrorX = 0.0d0
    do i = 1, nx
      symmetryErrorX = max(symmetryErrorX, dabs(xp(i) + xp(nx+1-i) - 1.0d0))
    enddo

    symmetryErrorY = 0.0d0
    do j = 1, ny
      symmetryErrorY = max(symmetryErrorY, dabs(yp(j) + yp(ny+1-j) - 1.0d0))
    enddo

    maxAspectRatio = 1.0d0
    do j = 1, ny
      dyLoc = wy(j)
      do i = 1, nx
        dxLoc = wx(i)
        ratio = dxLoc / dyLoc
        if (ratio < 1.0d0) ratio = 1.0d0 / ratio
        if (ratio > maxAspectRatio) maxAspectRatio = ratio
      enddo
    enddo

    return
  end subroutine build_quadrature_weights


  subroutine prepare_streaming_stencils()
    use commondata
    implicit none
    integer(kind=4) :: alpha, i, j
    real(kind=8) :: xsrc, ysrc
    integer(kind=4) :: idx(3)
    real(kind=8) :: w(3)
    logical :: ok

    do alpha = 0, 8
      do i = 1, nx
        if (ex(alpha) == 0) then
          stream_x_valid(alpha,i) = .true.
          stream_ix(alpha,i,1) = i
          stream_ix(alpha,i,2) = i
          stream_ix(alpha,i,3) = i
          stream_wx(alpha,i,1) = 1.0d0
          stream_wx(alpha,i,2) = 0.0d0
          stream_wx(alpha,i,3) = 0.0d0
        else
          xsrc = xp(i) - dble(ex(alpha))*uniformShift
          call build_lagrange_stencil_1d(nx, xp(1:nx), xsrc, idx, w, ok)
          stream_x_valid(alpha,i) = ok
          stream_ix(alpha,i,:) = idx
          stream_wx(alpha,i,:) = w
        endif
      enddo

      do j = 1, ny
        if (ey(alpha) == 0) then
          stream_y_valid(alpha,j) = .true.
          stream_iy(alpha,j,1) = j
          stream_iy(alpha,j,2) = j
          stream_iy(alpha,j,3) = j
          stream_wy(alpha,j,1) = 1.0d0
          stream_wy(alpha,j,2) = 0.0d0
          stream_wy(alpha,j,3) = 0.0d0
        else
          ysrc = yp(j) - dble(ey(alpha))*uniformShift
          call build_lagrange_stencil_1d(ny, yp(1:ny), ysrc, idx, w, ok)
          stream_y_valid(alpha,j) = ok
          stream_iy(alpha,j,:) = idx
          stream_wy(alpha,j,:) = w
        endif
      enddo
    enddo

    return
  end subroutine prepare_streaming_stencils


  subroutine build_lagrange_stencil_1d(n, xnodes, target, idx, w, ok)
    implicit none
    integer(kind=4), intent(in) :: n
    real(kind=8), intent(in) :: xnodes(1:n), target
    integer(kind=4), intent(out) :: idx(3)
    real(kind=8), intent(out) :: w(3)
    logical, intent(out) :: ok
    integer(kind=4) :: k, mid
    real(kind=8), parameter :: eps = 1.0d-12
    real(kind=8) :: xloc(3)

    idx = (/ 1, 1, 1 /)
    w = 0.0d0
    ok = .false.

    if (n < 3) return
    if (target < xnodes(1) - eps .or. target > xnodes(n) + eps) return

    if (target <= xnodes(2)) then
      mid = 2
    elseif (target >= xnodes(n-1)) then
      mid = n - 1
    else
      mid = 2
      do k = 2, n-1
        if (target >= xnodes(k) .and. target <= xnodes(k+1)) then
          if (dabs(target - xnodes(k)) <= dabs(xnodes(k+1) - target)) then
            mid = k
          else
            mid = k + 1
          endif
          exit
        endif
      enddo
      if (mid < 2) mid = 2
      if (mid > n-1) mid = n-1
    endif

    idx = (/ mid-1, mid, mid+1 /)
    xloc(1) = xnodes(idx(1))
    xloc(2) = xnodes(idx(2))
    xloc(3) = xnodes(idx(3))
    call lagrange_weights_3(xloc, target, w)
    ok = .true.

    return
  end subroutine build_lagrange_stencil_1d


  subroutine lagrange_weights_3(xnode, x0, w)
    implicit none
    real(kind=8), intent(in) :: xnode(3), x0
    real(kind=8), intent(out) :: w(3)

    w(1) = (x0 - xnode(2))*(x0 - xnode(3)) / ((xnode(1) - xnode(2))*(xnode(1) - xnode(3)))
    w(2) = (x0 - xnode(1))*(x0 - xnode(3)) / ((xnode(2) - xnode(1))*(xnode(2) - xnode(3)))
    w(3) = (x0 - xnode(1))*(x0 - xnode(2)) / ((xnode(3) - xnode(1))*(xnode(3) - xnode(2)))

    return
  end subroutine lagrange_weights_3


  real(kind=8) function lagrange_interp_3(xnode, fnode, x0)
    implicit none
    real(kind=8), intent(in) :: xnode(3), fnode(3), x0
    real(kind=8) :: w(3)

    call lagrange_weights_3(xnode, x0, w)
    lagrange_interp_3 = w(1)*fnode(1) + w(2)*fnode(2) + w(3)*fnode(3)

    return
  end function lagrange_interp_3


  real(kind=8) function lagrange_derivative_3(xnode, fnode, x0)
    implicit none
    real(kind=8), intent(in) :: xnode(3), fnode(3), x0

    lagrange_derivative_3 = &
      fnode(1) * (2.0d0*x0 - xnode(2) - xnode(3)) / ((xnode(1) - xnode(2))*(xnode(1) - xnode(3))) + &
      fnode(2) * (2.0d0*x0 - xnode(1) - xnode(3)) / ((xnode(2) - xnode(1))*(xnode(2) - xnode(3))) + &
      fnode(3) * (2.0d0*x0 - xnode(1) - xnode(2)) / ((xnode(3) - xnode(1))*(xnode(3) - xnode(2)))

    return
  end function lagrange_derivative_3


  subroutine interpolation_self_check()
    use commondata
    implicit none
    integer(kind=4) :: alpha, i, j, ii, jj
    real(kind=8) :: xsrc, ysrc
    real(kind=8) :: xnode(3), ynode(3)
    real(kind=8) :: errConst, errLin, errQuad, errBiquad, errUniform
    real(kind=8) :: s0, s1, s2, approx

    errConst = 0.0d0
    errLin = 0.0d0
    errQuad = 0.0d0
    errBiquad = 0.0d0
    errUniform = 0.0d0

    do alpha = 0, 8
      do i = 1, nx
        if (stream_x_valid(alpha,i)) then
          xsrc = xp(i) - dble(ex(alpha))*uniformShift
          do ii = 1, 3
            xnode(ii) = xp(stream_ix(alpha,i,ii))
          enddo
          s0 = stream_wx(alpha,i,1) + stream_wx(alpha,i,2) + stream_wx(alpha,i,3)
          s1 = stream_wx(alpha,i,1)*xnode(1) + stream_wx(alpha,i,2)*xnode(2) + stream_wx(alpha,i,3)*xnode(3)
          s2 = stream_wx(alpha,i,1)*xnode(1)*xnode(1) + stream_wx(alpha,i,2)*xnode(2)*xnode(2) + &
               stream_wx(alpha,i,3)*xnode(3)*xnode(3)
          errConst = max(errConst, dabs(s0 - 1.0d0))
          errLin = max(errLin, dabs(s1 - xsrc))
          errQuad = max(errQuad, dabs(s2 - xsrc*xsrc))
          if (meshMode == meshModeUniform .and. i-ex(alpha) >= 1 .and. i-ex(alpha) <= nx) then
            errUniform = max(errUniform, dabs(xsrc - xp(i-ex(alpha))))
          endif
        endif
      enddo

      do j = 1, ny
        if (stream_y_valid(alpha,j)) then
          ysrc = yp(j) - dble(ey(alpha))*uniformShift
          do jj = 1, 3
            ynode(jj) = yp(stream_iy(alpha,j,jj))
          enddo
          s0 = stream_wy(alpha,j,1) + stream_wy(alpha,j,2) + stream_wy(alpha,j,3)
          s1 = stream_wy(alpha,j,1)*ynode(1) + stream_wy(alpha,j,2)*ynode(2) + stream_wy(alpha,j,3)*ynode(3)
          s2 = stream_wy(alpha,j,1)*ynode(1)*ynode(1) + stream_wy(alpha,j,2)*ynode(2)*ynode(2) + &
               stream_wy(alpha,j,3)*ynode(3)*ynode(3)
          errConst = max(errConst, dabs(s0 - 1.0d0))
          errLin = max(errLin, dabs(s1 - ysrc))
          errQuad = max(errQuad, dabs(s2 - ysrc*ysrc))
          if (meshMode == meshModeUniform .and. j-ey(alpha) >= 1 .and. j-ey(alpha) <= ny) then
            errUniform = max(errUniform, dabs(ysrc - yp(j-ey(alpha))))
          endif
        endif
      enddo
    enddo

    do alpha = 0, 8
      if (ex(alpha) == 0 .or. ey(alpha) == 0) cycle
      do j = 1, ny
        if (.not. stream_y_valid(alpha,j)) cycle
        ysrc = yp(j) - dble(ey(alpha))*uniformShift
        do jj = 1, 3
          ynode(jj) = yp(stream_iy(alpha,j,jj))
        enddo
        do i = 1, nx
          if (.not. stream_x_valid(alpha,i)) cycle
          xsrc = xp(i) - dble(ex(alpha))*uniformShift
          do ii = 1, 3
            xnode(ii) = xp(stream_ix(alpha,i,ii))
          enddo
          approx = 0.0d0
          do jj = 1, 3
            do ii = 1, 3
              approx = approx + stream_wx(alpha,i,ii)*stream_wy(alpha,j,jj) * &
                (xnode(ii)*xnode(ii) + xnode(ii)*ynode(jj) + ynode(jj)*ynode(jj))
            enddo
          enddo
          errBiquad = max(errBiquad, dabs(approx - (xsrc*xsrc + xsrc*ysrc + ysrc*ysrc)))
        enddo
      enddo
    enddo

    open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
    write(00,*) "quadrature sum x/y/area =", quadSumX, quadSumY, quadSumArea
    write(00,*) "effective lengths =", Lx_eff, Ly_eff
    write(00,*) "max stretch ratio x/y =", maxStretchRatioX, maxStretchRatioY
    write(00,*) "max aspect ratio =", maxAspectRatio
    write(00,*) "symmetry error x/y =", symmetryErrorX, symmetryErrorY
    write(00,*) "interp err const/lin/quad/biquad =", errConst, errLin, errQuad, errBiquad
    if (meshMode == meshModeUniform) then
      write(00,*) "uniform backshift error =", errUniform
    endif
    close(00)

    return
  end subroutine interpolation_self_check


  real(kind=8) function interpolate_line_x(fieldLine, targetX)
    use commondata
    implicit none
    real(kind=8), intent(in) :: fieldLine(1:nx)
    real(kind=8), intent(in) :: targetX
    integer(kind=4) :: idx(3)
    real(kind=8) :: w(3), xnode(3), fnode(3)
    logical :: ok
    real(kind=8) :: lagrange_interp_3

    call build_lagrange_stencil_1d(nx, xp(1:nx), targetX, idx, w, ok)
    if (.not. ok) then
      interpolate_line_x = 0.0d0
      return
    endif

    xnode(1) = xp(idx(1))
    xnode(2) = xp(idx(2))
    xnode(3) = xp(idx(3))
    fnode(1) = fieldLine(idx(1))
    fnode(2) = fieldLine(idx(2))
    fnode(3) = fieldLine(idx(3))
    interpolate_line_x = lagrange_interp_3(xnode, fnode, targetX)

    return
  end function interpolate_line_x


  real(kind=8) function interpolate_line_y(fieldLine, targetY)
    use commondata
    implicit none
    real(kind=8), intent(in) :: fieldLine(1:ny)
    real(kind=8), intent(in) :: targetY
    integer(kind=4) :: idx(3)
    real(kind=8) :: w(3), ynode(3), fnode(3)
    logical :: ok
    real(kind=8) :: lagrange_interp_3

    call build_lagrange_stencil_1d(ny, yp(1:ny), targetY, idx, w, ok)
    if (.not. ok) then
      interpolate_line_y = 0.0d0
      return
    endif

    ynode(1) = yp(idx(1))
    ynode(2) = yp(idx(2))
    ynode(3) = yp(idx(3))
    fnode(1) = fieldLine(idx(1))
    fnode(2) = fieldLine(idx(2))
    fnode(3) = fieldLine(idx(3))
    interpolate_line_y = lagrange_interp_3(ynode, fnode, targetY)

    return
  end function interpolate_line_y


  real(kind=8) function derivative_line_x(fieldLine, targetX)
    use commondata
    implicit none
    real(kind=8), intent(in) :: fieldLine(1:nx)
    real(kind=8), intent(in) :: targetX
    integer(kind=4) :: idx(3)
    real(kind=8) :: w(3), xnode(3), fnode(3)
    logical :: ok
    real(kind=8) :: lagrange_derivative_3

    call build_lagrange_stencil_1d(nx, xp(1:nx), targetX, idx, w, ok)
    if (.not. ok) then
      derivative_line_x = 0.0d0
      return
    endif

    xnode(1) = xp(idx(1))
    xnode(2) = xp(idx(2))
    xnode(3) = xp(idx(3))
    fnode(1) = fieldLine(idx(1))
    fnode(2) = fieldLine(idx(2))
    fnode(3) = fieldLine(idx(3))
    derivative_line_x = lagrange_derivative_3(xnode, fnode, targetX)

    return
  end function derivative_line_x


  real(kind=8) function derivative_line_y(fieldLine, targetY)
    use commondata
    implicit none
    real(kind=8), intent(in) :: fieldLine(1:ny)
    real(kind=8), intent(in) :: targetY
    integer(kind=4) :: idx(3)
    real(kind=8) :: w(3), ynode(3), fnode(3)
    logical :: ok
    real(kind=8) :: lagrange_derivative_3

    call build_lagrange_stencil_1d(ny, yp(1:ny), targetY, idx, w, ok)
    if (.not. ok) then
      derivative_line_y = 0.0d0
      return
    endif

    ynode(1) = yp(idx(1))
    ynode(2) = yp(idx(2))
    ynode(3) = yp(idx(3))
    fnode(1) = fieldLine(idx(1))
    fnode(2) = fieldLine(idx(2))
    fnode(3) = fieldLine(idx(3))
    derivative_line_y = lagrange_derivative_3(ynode, fnode, targetY)

    return
  end function derivative_line_y


  real(kind=8) function wall_derivative_x_left(twall, t1, t2)
    use commondata
    implicit none
    real(kind=8), intent(in) :: twall, t1, t2
    real(kind=8) :: dx1, dx2, q

    dx1 = xp(1)
    dx2 = xp(2) - xp(1)
    q = dx2 / dx1
    wall_derivative_x_left = ( -4.0d0*q*(q+1.0d0)*twall + (2.0d0*q+1.0d0)**2*t1 - t2 ) / &
      ( q*(2.0d0*q+1.0d0)*dx1 )

    return
  end function wall_derivative_x_left


  real(kind=8) function wall_derivative_x_right(tnm1, tn, twall)
    use commondata
    implicit none
    real(kind=8), intent(in) :: tnm1, tn, twall
    real(kind=8) :: dx1, dx2, q

    dx1 = 1.0d0 - xp(nx)
    dx2 = xp(nx) - xp(nx-1)
    q = dx2 / dx1
    wall_derivative_x_right = ( 4.0d0*q*(q+1.0d0)*twall - (2.0d0*q+1.0d0)**2*tn + tnm1 ) / &
      ( q*(2.0d0*q+1.0d0)*dx1 )

    return
  end function wall_derivative_x_right


  real(kind=8) function wall_derivative_y_bottom(twall, t1, t2)
    use commondata
    implicit none
    real(kind=8), intent(in) :: twall, t1, t2
    real(kind=8) :: dy1, dy2, q

    dy1 = yp(1)
    dy2 = yp(2) - yp(1)
    q = dy2 / dy1
    wall_derivative_y_bottom = ( -4.0d0*q*(q+1.0d0)*twall + (2.0d0*q+1.0d0)**2*t1 - t2 ) / &
      ( q*(2.0d0*q+1.0d0)*dy1 )

    return
  end function wall_derivative_y_bottom


  real(kind=8) function wall_derivative_y_top(tnm1, tn, twall)
    use commondata
    implicit none
    real(kind=8), intent(in) :: tnm1, tn, twall
    real(kind=8) :: dy1, dy2, q

    dy1 = 1.0d0 - yp(ny)
    dy2 = yp(ny) - yp(ny-1)
    q = dy2 / dy1
    wall_derivative_y_top = ( 4.0d0*q*(q+1.0d0)*twall - (2.0d0*q+1.0d0)**2*tn + tnm1 ) / &
      ( q*(2.0d0*q+1.0d0)*dy1 )

    return
  end function wall_derivative_y_top


!===========================================================================================================================
!濡備笅閮芥槸瀛愮▼搴?

  subroutine initial()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: un(0:8)
    real(kind=8) :: us2, Bx, By
    character(len=100) :: reloadFileName
    real(kind=8) :: dev
    

    itc = 0
    errorU = 100.0d0
    errorT = 100.0d0 

    call setup_mesh_and_runtime()


    !-----------------------------------------------------------------------------------------------
    !璁板綍鍚勭淇℃伅鍦ㄦ棩蹇楁枃浠朵腑
    open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')  !鍦ㄨ繖涓猼xt鏂囦欢鍚庨潰缁х画鍐欙紙杩藉姞妯″紡锛?
    
    if(outputBinFile.EQ.1) then
        open(unit=01,file=trim(binFolderPrefix)//"-"//"readme",status="unknown")    !trim鍘绘帀瀛楃涓插熬閮ㄧ┖鏍硷紝鎹簡瀛樺偍璺緞锛屽彲鑷繁鏇存敼
        write(01,*) "binFile folder exist!"                                         !璇诲彇璺緞binFolderPrefix="../binFile/buoyancyCavity
        close(01)
        write(00,*) "Data will be stored in ", binFolderPrefix
    endif
    if(outputPltFile.EQ.1) then
        open(unit=01,file=trim(pltFolderPrefix)//"-"//"readme",status="unknown")     !璇诲彇璺緞pltFolderPrefix="../pltFile/buoyancyCavity
        write(01,*) "pltFile folder exist!"
        close(01)
        write(00,*) "Data will be stored in ", pltFolderPrefix
    endif
    
    !if( (paraA.GE.1.0d0).OR.(paraA.LE.-4.0d0) ) then                           !鍙湁鍦╗-4,1]鎵嶅彲浠ワ紝瑕佷笉鐒堕璀?
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
    write(00,*) 'meshMode =', meshMode, '; stretchA =', stretchA
    write(00,*) 'Rayleigh=',real(Rayleigh,kind=8), '; Prandtl =',real(Prandtl,kind=8), '; Mach =',real(Mach,kind=8)
    write(00,*) "Length unit: L0 =", real(lengthUnit,kind=8)
    write(00,*) "Time unit: Sqrt(L0/(gBeta*DeltaT)) =", real(timeUnit,kind=8)
    write(00,*) "Velocity unit: Sqrt(gBeta*L0*DeltaT) =", real(velocityUnit,kind=8)
    write(00,*) "   "
    write(00,*) 'tauf=',real(tauf,kind=8)
    write(00,*) 'taug=',real(taug,kind=8)
    write(00,*) "viscosity =",real(viscosity,kind=8), "; diffusivity =",real(diffusivity,kind=8)
    write(00,*) "outputFrequency =", real(outputFrequency,kind=8), "tf"
    write(00,*) "......................  or ",  int(outputFrequency*timeUnit), "in itc units"
    write(00,*) "backupInterval =", backupInterval, " free-fall time units"
    write(00,*) ".................... or ", int(backupInterval/outputFrequency)*int(outputFrequency*timeUnit), "itc units"
    if(loadInitField.EQ.1) then
        write(00,*) "reloadDimensionlessTime=", reloadDimensionlessTime
    endif
    write(00,*) "itc_max =",itc_max
    write(00,*) "default epsU =", real(epsU,kind=8),"; epsT =", real(epsT,kind=8)
    write(00,*) "useG =", useG
    write(00,*) "quadrature sums =", quadSumX, quadSumY, quadSumArea
    write(00,*) "mesh symmetry =", symmetryErrorX, symmetryErrorY
    write(00,*) "stretch ratios =", maxStretchRatioX, maxStretchRatioY, "; aspect =", maxAspectRatio
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
    !鑺傜偣鍧愭爣鏁扮粍
    xp(0) = 0.0d0
    xp(nx+1) = dble(nx)
    do i=1,nx
        xp(i) = dble(i)-0.5d0      ! 0 | 0.5 1.5 2.5 3.5 | 4   锛屽寘鍚竟鐣岀偣
    enddo
    xp = xp / lengthUnit           !鐢ㄥ崟浣嶉暱搴︽棤閲忕翰鍖?

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
    call build_mesh_coordinates()
    call build_quadrature_weights()
    call prepare_streaming_stencils()
    !-----------------------------------------------------------------------------------------------



    !-----------------------------------------------------------------------------------------------
    !鍒濆鍖?
    rho = 1.0d0                                     !瀵嗗害rho=1
    
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

    if(loadInitField.EQ.0) then                    !鍦ㄤ笉鍔犺浇鏂囦欢鐨勬儏鍐典笅锛岄兘鏄浂鍦轰负鍒濆€?
    
        u = 0.0d0
        v = 0.0d0
        T = 0.0d0
        
        write(00,*) "Initial field is set exactly"
        if(reloadDimensionlessTime.NE.0) then        !鍦ㄤ笉鍔犺浇鏂囦欢鐨勬儏鍐典笅锛宺eloadDimensionlessTime蹇呴』鏄浂
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
    do j = 1, ny                                   !鍦ㄤ笉鍔犺浇鏂囦欢鐨勬儏鍐典笅锛屽垵濮嬪寲娓╁害鏄垎灞傜殑锛屽叾瀹炶繖涓湁鐐归棶棰橈紝杈圭晫鎵嶆槸绮剧‘鐨凾hot鍜孴cold
        !T(1,j) = Thot
        !T(nx,j) = Tcold
        do i = 1, nx
            T(i,j) = Thot + (Tcold-Thot)*xp(i)
        enddo
    enddo
    write(00,*) "Temperature B.C. for vertical walls are:===Hot/cold wall==="
#endif

#ifdef HorizontalWallsConstT
    do i = 1, nx
        !T(i,1) = Thot
        !T(i,ny) = Tcold
        do j = 1, ny
            T(i,j) = Thot + (Tcold-Thot)*yp(j)
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
                    f(alpha,i,j) = rho(i,j)*omega(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)  !D2Q9鏍囧噯feq
                enddo
                do alpha = 0, 4
                    un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha) 
                    !g(alpha,i,j) = omegaT(alpha)*T(i,j)*(1.0d0+10.0d0/(4.0d0+paraA)*un(alpha))                        !D2Q5绾挎€eq锛屼笉鍚湁u^2椤癸紝杩欎釜鏄惈鏈変慨姝ｇ殑geq
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




    elseif(loadInitField.EQ.1) then                               !鍦ㄥ姞杞芥枃浠剁殑鎯呭喌涓嬶紝璇诲彇璺緞 reloadFilePrefix="./reloadFile/backupFile
        if(reloadDimensionlessTime.EQ.0) then                     !鍦ㄥ姞杞芥枃浠剁殑鎯呭喌涓嬶紝reloadDimensionlessTime鏈€濂芥槸闈為浂
            write(00,*) "WARNING: since loadInitField.EQ.1, please confirm reloadDimensionlessTime", reloadDimensionlessTime
            !stop
        endif
        write(00,*) "Load initial field from previous simulation: ../reloadFile/backupFile- >>>"
        write(reloadFileName, *) reloadbinFileNum                 !鎹簡涓彉閲廚ame
        reloadFileName = adjustl(reloadFileName)                  !adjustl鎶婂瓧绗︿覆宸﹀榻愶紝鎶婂墠瀵肩┖鏍肩Щ鍒板瓧绗︿覆鏈熬
        
        open(unit=01,file=trim(reloadFilePrefix)//"-"//trim(reloadFileName)//".bin",form="unformatted", &
        access="sequential",status='old')  !unformatted鏄簩杩涘埗,sequential锛氭寜璁板綍椤哄簭璇诲啓
            if(loadInitFG.EQ.1) then
                write(00,*) "Reloading f and g from file"
                read(01) (((f(alpha,i,j), i=1,nx), j=1,ny), alpha=0,8)      !鍏?i锛屽啀 j锛屽啀 alpha
                read(01) (((g(alpha,i,j), i=1,nx), j=1,ny), alpha=0,4)
            else
                write(00,*) "Not reloading f and g from file"               !濡傛灉璇荤殑鏄?backupData鐢熸垚鐨?backupFile-*.bin锛歭oadInitFG 蹇呴』涓?1锛堝惁鍒欒閿欎綅锛?
                stop               
            endif

            write(00,*) "Reloading u, v, T from file"
            read(01) ((u(i,j), i=1,nx), j=1,ny)
            read(01) ((v(i,j), i=1,nx), j=1,ny)
            read(01) ((T(i,j), i=1,nx), j=1,ny)
            
            if(loadInitRho.EQ.1) then
                write(00,*) "Reloading rho from file"                       !濡傛灉璇荤殑鏄?backupData鐢熸垚鐨?backupFile-*.bin锛歭oadInitRho 蹇呴』涓?0锛堝洜涓烘枃浠堕噷娌℃湁 rho锛?
                read(01) ((rho(i,j), i=1,nx), j=1,ny)
            else
                write(00,*) "Not reloading rho from file, rho will be set to 1.0"
                rho = 1.0d0
            endif
        close(01)

        if(loadInitFG.EQ.0) then                                            !鍦ㄥ姞杞芥枃浠剁殑鎯呭喌涓?鐢ㄥ钩琛″垎甯冮噸鏂版瀯閫?f,g 
            write(00,*) "Not reloading f and g from file, f and g will be set as equlibrium distribution functions"
            do j = 1, ny
                do i = 1, nx
                    us2 = u(i,j)*u(i,j)+v(i,j)*v(i,j)
                    do alpha = 0, 8
                        un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                        f(alpha,i,j) = rho(i,j)*omega(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)
                    enddo
                    do alpha = 0, 4
                        un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                        g(alpha,i,j) = omegaT(alpha)*T(i,j)*(1.0d0+10.0d0/(4.0d0+paraA)*un(alpha))
                    enddo
                enddo
            enddo
        endif

        dev = 0.0d0
        do j=1,ny
            do i=1,nx
                rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)     !鐢ㄥ垎甯冨嚱鏁?f 閲嶅缓瀵嗗害 rho
                dev = dev+g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)-T(i,j)                                   !璁＄畻娓╁害涓€鑷存€у亸宸?
            enddo
        enddo
        write(00,*) "RELOAD: Deviation in temperature: ", real(dev,kind=8)
        if(dev.GT.1.0d0) then
            write(00,*) "Error: too large Deviation when reload data!"                                          !杈撳嚭鍋忓樊骞跺仛鈥滅矖闃堝€尖€濆垽閿?
            stop
        endif
        write(00,*) "Raw data is loaded from the file: backupFile-",trim(reloadFileName),".bin"                 !鎵撳嵃鈥滄棫鏂囦欢鈥濅俊鎭?
    else
        write(00,*) "Error: initial field is not properly set"                                                  !濡傛灉 loadInitField 涓嶆槸 0/1 鎴栭€昏緫涓嶄竴鑷达紝鐩存帴鍋滄
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



  subroutine collision()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: m(0:8), m_post(0:8), meq(0:8)
    real(kind=8) :: s(0:8)
    real(kind=8) :: fSource(0:8)

    !$omp parallel do default(none) shared(f,f_post,rho,u,v,Fx,Fy,T,Snu,Sq,gBeta) private(i,j,alpha,s,m,m_post,meq,fSource) 
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
          Fy(i,j) = rho(i,j)*gBeta*(T(i,j)-Tref)        !鍔ㄩ噺鏂圭▼涓婄殑婧愰」锛屽嵆娴姏椤?


#ifdef    SideHeatedHa
          Fx(i,j) = 0.0d0+B2sigemarho*(v(i,j)*sin(phi)*cos(phi)-u(i,j)*sin(phi)*sin(phi))
          Fy(i,j) = rho(i,j)*gBeta*(T(i,j)-Tref)+ rho(i,j)*B2sigemarho*(u(i,j)*sin(phi)*cos(phi)&
          -v(i,j)*cos(phi)*cos(phi))                    !鍔ㄩ噺鏂圭▼涓婄殑婧愰」锛屽嵆娴姏椤瑰姞纾佸満
#endif


          fSource(0) = 0.0d0                                                       !灏嗘簮椤笷瀵瑰簲鐨勮础鐚姇褰卞埌鍚勪釜鐭╀腑锛屽苟鍋氬崐姝ヤ慨姝?
          fSource(1) = (6.0d0-3.0d0*s(1))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
          fSource(2) = -(6.0d0-3.0d0*s(2))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
          fSource(3) = (1.0d0-0.5d0*s(3))*Fx(i,j)
          fSource(4) = -(1.0d0-0.5d0*s(4))*Fx(i,j)
          fSource(5) = (1.0d0-0.5d0*s(5))*Fy(i,j)
          fSource(6) = -(1.0d0-0.5d0*s(6))*Fy(i,j)
          fSource(7) = (2.0d0-s(7))*(u(i,j)*Fx(i,j)-v(i,j)*Fy(i,j))
          fSource(8) = (1.0d0-0.5d0*s(8))*(u(i,j)*Fy(i,j)+v(i,j)*Fx(i,j))     !杩欒竟鏄箻浠鍙樺埌鐭╃┖闂达紝鐒跺悗鍐嶄箻浠?-1/2S淇

          do alpha = 0, 8
            m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+fSource(alpha)     !鐭╃┖闂寸鎾?
          enddo

          f_post(0,i,j) = m_post(0)/9.0d0-m_post(1)/9.0d0+m_post(2)/9.0d0                                         !杩欒竟鏄箻浠閫?
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


  subroutine streaming()                                    !鍏堣縼绉伙紝鍐嶈竟鐣屽鐞?
    use commondata                                            !杩佺Щ姝ラ锛歱ull streaming锛屾妸纰版挒鍚庣殑 f_post 鎷夊彇鍒板綋鍓嶆牸鐐?
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: ip, jp
    integer(kind=4) :: alpha
    
    !$omp parallel do default(none) shared(f,f_post,ex,ey) private(i,j,ip,jp,alpha)
    do j = 1, ny
        do i = 1, nx
            do alpha = 0, 8                        !涓婃父鏍肩偣绱㈠紩锛歠伪(i,j) <- f_post伪(i-ex伪, j-ey伪)
                ip = i-ex(alpha)                   !杈圭晫闄勮繎 (ip/jp 鍙兘涓?0 鎴?nx+1/ny+1)锛岄渶鍦?bounceback/鍛ㄦ湡杈圭晫澶勭悊涓鐩栦慨姝ｈ竟鐣屽垎甯?
                jp = j-ey(alpha)                   !ghost 灞傚湪鍒濆鍖栦腑涓?0锛屼繚璇佷笉浼氬嚭鐜版湭鍒濆鍖栧瀮鍦惧€?
                
                f(alpha,i,j) = f_post(alpha,ip,jp)
            enddo
        enddo
    enddo
    !$omp end parallel do

    return
  end subroutine streaming




  subroutine bounceback()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    ! integer(kind=4) :: alpha

#ifdef VerticalWallsPeriodicalU      
    !$omp parallel do default(none) shared(f, f_post) private(j)
    do j = 1, ny                                                  !閫熷害杈圭晫鍨傜洿杈圭晫鍛ㄦ湡锛岀洿鎺ユ柟鍚戠浉鍚岋紝璺ㄨ竟鐣岀殑鍏ュ皠鍒嗗竷      
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
    do j = 1, ny                                                 !閫熷害杈圭晫鍨傜洿杈圭晫闈欐澹佹棤婊戠Щ锛岀洿鎺ュ弽寮癸紝鏂瑰悜鐩稿弽
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
    do i = 1, nx                                                  !閫熷害杈圭晫姘村钩杈圭晫鏃犳粦绉伙紝鐩存帴鍙嶅脊锛屾柟鍚戠浉鍙?
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

  !streaming锛氬厛鎶?interior 鐨?f 鎷夊ソ锛堣竟鐣屽鍙兘涓嶅叏瀵癸級
  !bounceback锛氱敤 f_post 鐨勨€滃嚭灏勫垎甯冣€濆幓濉€滃叆灏勫垎甯冣€濓紙杩欏搴?half-way bounceback 鐨勫父瑙佸疄鐜版柟寮忥級


  subroutine macro()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !$omp parallel do default(none) shared(f,rho,u,v,Fx,Fy) private(i,j)
    do j = 1, ny
        do i = 1, nx
            rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
            u(i,j) = ( f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)+0.5d0*Fx(i,j) )/rho(i,j)     !鍚姏LBM鐨勫崐姝ュ姩閲忎慨姝ｏ細rho*u = 危 f e + 0.5*F锛屽搴擥uo forcing鐨勪簩闃跺畾涔?
            v(i,j) = ( f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)+0.5d0*Fy(i,j) )/rho(i,j)
        enddo
    enddo
    !$omp end parallel do
    
    return
  end subroutine macro




  subroutine collisionT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: n(0:4), n_post(0:4), neq(0:4)
    real(kind=8) :: q(0:4)
    real(kind=8) :: Bx, By
    real(kind=8) :: dBx, dBy
    real(kind=8) :: SG
    SG = 1.0d0 - 0.5d0*s_j

    !$omp parallel do default(none) shared(g,g_post,u,v,T,Bx_prev,By_prev,SG,s_j,s_e,s_q) private(i,j,alpha,n,neq,q,n_post,Bx,By,dBx,dBy) 
    do j = 1, ny
        do i = 1, nx

            Bx = u(i,j) * T(i,j)
            By = v(i,j) * T(i,j)

            if (useG) then
              dBx = Bx - Bx_prev(i,j)
              dBy = By - By_prev(i,j)
              Bx_prev(i,j) = Bx
              By_prev(i,j) = By
            else
              dBx = 0.0d0
              dBy = 0.0d0
            end if



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


#ifdef steadyFlow
    subroutine check()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: error1, error2, error5, error6



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
    
    errorU = dsqrt(error1)/dsqrt(error2)                 !閫熷害鍦虹浉瀵筁2璇樊锛殀|u^n-u^{n-1}||_2 / ||u^n||_2
    errorT = error5/error6                               !娓╁害鍦虹浉瀵筁1璇樊锛殀|T^n-T^{n-1}||_1 / ||T^n||_1

    !open(unit=01,file='convergence.plt',status='unknown',position='append')  !璁板綍鏀舵暃鍘嗗彶鏁版嵁锛坕tc, errorU, errorT锛夛紝鐢ㄤ簬鍚庡鐞嗙敾鏀舵暃鏇茬嚎锛涙湰韬笉鍋氱粯鍥?
    !write(01,*) itc,' ',errorU,' ',errorT
    !close(01)
    !write(*,*) itc,' ',errorU,' ',errorT    

    call append_convergence_tecplot('convergence.plt', itc, errorU, errorT)

    !write(caseTag,'("nx=",I0,",ny=",I0,",useG=",L1)') nx, ny, useG         !杈撳嚭鏀舵暃鏇茬嚎鐨勫姣?
    !call append_convergence_master_tecplot('convergence_all.plt', caseTag, itc, errorU, errorT)

    write(*,'(I12,1X,ES24.16,1X,ES24.16)') itc, errorU, errorT


    return
    end subroutine check
#endif



subroutine append_convergence_tecplot(filename, itc, errorU, errorT)
  implicit none
  character(len=*), intent(in) :: filename
  integer(kind=4), intent(in)  :: itc
  real(kind=8),    intent(in)  :: errorU, errorT

  integer :: u
  logical, save :: first_write = .true.

  if (first_write) then
    ! 姣忔绋嬪簭鏂拌繍琛岀殑绗竴娆¤皟鐢細鐩存帴瑕嗙洊鏃ф枃浠?
    open(newunit=u, file=trim(filename), status='replace', action='write', form='formatted')

    write(u,'(A)') 'TITLE = "Convergence history; itc='//trim(adjustl(itoa(itc)))// &
                   ', errorU='//trim(adjustl(dtoa(errorU)))// &
                   ', errorT='//trim(adjustl(dtoa(errorT)))//'"'
    write(u,'(A)') 'VARIABLES = "itc" "errorU" "errorT"'
    write(u,'(A)') 'ZONE T="conv", F=POINT'
    write(u,'(I12,1X,ES24.16,1X,ES24.16)') itc, errorU, errorT
    close(u)

    first_write = .false.
  else
    ! 鍚屼竴娆¤繍琛岀殑鍚庣画璋冪敤锛氳拷鍔犳暟鎹
    open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
    write(u,'(I12,1X,ES24.16,1X,ES24.16)') itc, errorU, errorT
    close(u)
  end if

contains
  function itoa(i) result(s)
    integer(kind=4), intent(in) :: i
    character(len=32) :: s
    write(s,'(I0)') i
  end function itoa

  function dtoa(x) result(s)
    real(kind=8), intent(in) :: x
    character(len=64) :: s
    write(s,'(ES16.8)') x
  end function dtoa
end subroutine append_convergence_tecplot



subroutine append_convergence_master_tecplot(filename, zoneName, itc, errorU, errorT)
  implicit none
  character(len=*), intent(in) :: filename, zoneName
  integer(kind=4), intent(in)  :: itc
  real(kind=8),    intent(in)  :: errorU, errorT

  logical :: ex
  integer :: u
  logical, save :: zone_started = .false.

  ! 鏈杩愯绗竴娆″啓锛氱‘淇濇枃浠跺ご瀛樺湪 + 鍐欏叆涓€涓柊ZONE
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

  ! 杩藉姞涓€涓暟鎹偣
  open(newunit=u, file=trim(filename), status='old', position='append', action='write', form='formatted')
  write(u,'(I12,1X,ES24.16,1X,ES24.16)') itc, errorU, errorT
  close(u)
end subroutine append_convergence_master_tecplot





  subroutine output_binary()                                         !杈撳嚭uvTrho锛屽瓨鍌ㄥ湪binFolderPrefix="../binFile/buoyancyCavity-000000001234.bin
    use commondata                                                   !鐢ㄤ簬鍚庡鐞嗗揩鐓э紱閲嶅惎璇诲叆鏃跺繀椤绘寜 u,v,T,rho 椤哄簭璇诲彇
    implicit none
    integer(kind=4) :: i, j
    character(len=100) :: filename
    
#ifdef steadyFlow
    write(filename,*) itc                                            !steadyFlow锛歜in/plt鏂囦欢鍚嶉€氬父鐩存帴鐢?itc 鏉ョ紪鍐欙紙鏍煎瓙鏃堕棿姝ラ暱锛?
#endif

#ifdef unsteadyFlow
    binFileNum = binFileNum+1                 !unsteadyFlow锛氭枃浠跺悕鐢ㄨ緭鍑哄簭鍙?binFileNum锛堟瘡娆¤緭鍑鸿嚜澧烇紝涓?tff 闂撮殧鐢?outputFrequency 鎺у埗锛?
    if(loadInitField.EQ.0) write(filename,*) binFileNum
    if(loadInitField.EQ.1) write(filename,*) binFileNum+reloadbinFileNum
#endif

    filename = adjustl(filename)

    open(unit=03,file=trim(binFolderPrefix)//"-"//trim(filename)//'.bin',form="unformatted",access="sequential")    !浜岃繘鍒?
    write(03) ((u(i,j),i=1,nx),j=1,ny)
    write(03) ((v(i,j),i=1,nx),j=1,ny)
    write(03) ((T(i,j),i=1,nx),j=1,ny)
    write(03) ((rho(i,j), i=1,nx), j=1,ny)
    close(03)

    return
  end subroutine output_binary


    

  subroutine backupData()                                         !杈撳嚭fguvT锛屽瓨鍌ㄥ湪褰撳墠璺緞锛屽悕瀛楁槸backupFile-1000.bin
    use commondata                                                !鐢ㄤ簬閲嶅惎锛屽寘鍚?f,g锛涜鍏ユ椂蹇呴』鍏堣 f,g 鍐嶈 u,v,T锛堟棤 rho锛?
    implicit none
    integer(kind=4) :: i, j, alpha
    character(len=100) :: filename

#ifdef steadyFlow
    write(filename,*) itc                                    !steadyFlow锛歜in/plt鏂囦欢鍚嶉€氬父鐩存帴鐢?itc 鏉ョ紪鍐欙紙鏍煎瓙鏃堕棿姝ラ暱锛?                              
#endif

#ifdef unsteadyFlow
    if(loadInitField.EQ.0) write(filename,*) binFileNum     !unsteadyFlow锛氭枃浠跺悕鐢ㄨ緭鍑哄簭鍙?binFileNum锛堟瘡娆¤緭鍑鸿嚜澧烇紝涓?tff 闂撮殧鐢?outputFrequency 鎺у埗锛?
    if(loadInitField.EQ.1) write(filename,*) binFileNum+reloadbinFileNum
#endif

    filename = adjustl(filename)

    open(unit=05,file='backupFile-'//trim(filename)//'.bin',form="unformatted",access="sequential")   !浜岃繘鍒?
    write(05) (((f(alpha,i,j), i=1,nx), j=1,ny), alpha=0,8)
    write(05) (((g(alpha,i,j), i=1,nx), j=1,ny), alpha=0,4)
    write(05) ((u(i,j), i=1,nx), j=1,ny)
    write(05) ((v(i,j), i=1,nx), j=1,ny)
    write(05) ((T(i,j), i=1,nx), j=1,ny)
    close(05)
    
    open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
    write(00,*) "Backup  f, g, u, v, T to the file: backupFile-", trim(filename),".bin"
    close(00)
    
    return
  end subroutine backupData

    
    

  subroutine output_Tecplot()                        !杈撳嚭浜岃繘鍒舵枃浠?
    use commondata
    implicit none
    integer(kind=4) :: i, j, k
    REAL(kind=4) :: zoneMarker, eohMarker   !Tecplot 浜岃繘鍒舵牸寮忛噷鐢ㄧ殑涓や釜鈥滄爣璁板€尖€濓紙299 鍜?357锛夛紝鐢ㄤ簬鍛婅瘔 Tecplot锛氳繖閲屽紑濮嬫槸 zone 鎻忚堪 / header 缁撴潫銆?
    character(len=40) :: title              !鏂囦欢 Title 瀛楃涓?
    character(len=40) :: V1,V2,V3,V4,V5     !鍙橀噺鍚嶅瓧绗︿覆锛圶,Y,U,V,T锛?
    integer(kind=4), parameter :: kmax=1    !浜岀淮鏁版嵁涔熸寜 3D 鐨?IJK 鍐欙紝K=1
    character(len=40) :: zoneName           !zone 鍚嶇О
    character(len=100) :: filename          !杈撳嚭鏂囦欢鍚嶅瓧绗︿覆
    
    !$acc update self(u,v,T)                !OpenACC鎸囦护锛屽拷鐣?

#ifdef steadyFlow
    write(filename,'(i12.12)') itc
#endif

#ifdef unsteadyFlow
    pltFileNum = pltFileNum+1               !鍜屽墠闈㈢殑binFileNum涓€鏍?
    write(filename,'(i12.12)') pltFileNum
#endif

    filename = adjustl(filename)            !瀛樺偍璺緞 pltFolderPrefix="./pltFile/buoyancyCavity000000000034.plt

    open(41,file=trim(pltFolderPrefix)//"-"//trim(filename)//'.plt', access='stream', form='unformatted')    !stream锛氬瓧鑺傛祦

    !---------------------------------------------
    zoneMarker= 299.0                                     !鍥哄畾鐨勶紝涓嶉渶瑕佷慨鏀?
    eohMarker = 357.0

    !I. HEAD SECTION--------------------------------------
    !c--Magic number, Version number
    write(41) "#!TDV101"                                  !Tecplot 璇嗗埆浜岃繘鍒?.plt 鐨勨€滈瓟鏁?鐗堟湰鍙封€濆瓧绗︿覆,涓嶉渶瑕佷慨鏀?

    !c--Integer value of 1
    write(41) 1                                           !Tecplot 瑙勮寖閲岀揣璺熺潃涓€涓暣鍨嬪€硷紙閫氬父琛ㄧず瀛楄妭搴?鏂囦欢绫诲瀷鐗堟湰绛夋帶鍒跺瓧娈碉級,涓嶉渶瑕佷慨鏀?

    Title="MyFirst"                                       !Tecplot 鐨勪簩杩涘埗鏍煎紡閲屽瓧绗︿覆涓嶆槸鐩存帴 write锛岃€屾槸閫愬瓧绗﹀啓 ASCII 鐮侊紝鍐嶄互 0 缁撳熬
    call dumpstring(title)                                !dumpstring() 灏辨槸骞茶繖涓殑

    !c-- Number of variables in this data file
    write(41) 5                                           !鏈?5 涓彉閲忥細X, Y, U, V, T,濡傛灉闇€瑕佷慨鏀癸紝杩欎釜鏈夊彉鍖栭渶瑕佷慨鏀?

    !c-- Variable names.                                  !鍙橀噺鍚嶄緷娆″啓鍏?鏈夊彉鍖栭渶瑕佷慨鏀?
    V1='X'
    call dumpstring(V1)
    V2='Y'
    call dumpstring(V2)
    V3='U'
    call dumpstring(V3)
    V4='V'
    call dumpstring(V4)
    V5='T'
    call dumpstring(V5)

    !c-----Zones-----------------------------

    !c--------Zone marker. Value = 299.0                   !鍐欏叆 float 299.0锛屽憡璇?Tecplot锛氫竴涓?zone 鐨勬弿杩板紑濮嬩簡銆?
    write(41) zoneMarker

    !--------Zone name.
    zoneName="ZONE 001"
    call dumpstring(zoneName)                              !zone 鍚嶇О锛屾樉绀哄湪 Tecplot 閲岋紙姣斿 Zones 闈㈡澘閲岋級

    !---------Zone Color                                   !-1 琛ㄧず浣跨敤榛樿閰嶈壊
    write(41) -1

    !---------ZoneType                                     !0 閫氬父琛ㄧず ORDERED锛堢粨鏋勭綉鏍?IJK锛?
    write(41) 0

    !---------DataPacking 0=Block, 1=Point                 !0 = Block锛堝厛鎶婃暣涓?X 鍐欏畬锛屽啀鍐欐暣涓?Y鈥︼級,1 = Point锛堟瘡涓偣渚濇鍐?X,Y,U,V,T锛?

    write(41) 1

    !---------Specify Var Location. 0 = Do not specify, all data    !0锛氫笉鎸囧畾锛岄粯璁ゆ墍鏈夊彉閲忛兘鍦ㄨ妭鐐逛笂锛坣odal锛?1锛氫細璺熺潃涓€涓插彉閲忎綅缃垪琛紙cell-centered/nodal锛?
    !---------is located at the nodes. 1 = Specify
    write(41) 0

    !---------Number of user defined face neighbor connections   !鑷畾涔夐偦鎺ヨ繛鎺ユ暟锛涚粨鏋勭綉鏍间竴鑸笉闇€瑕侊紝鍐?0
    ! (value >= 0)
    write(41) 0

    !---------IMax,JMax,KMax                          !缁撴瀯缃戞牸灏哄锛欼=nx, J=ny, K=1, 鏈夊彉鍖栭渶瑕佷慨鏀?
    write(41) nx
    write(41) ny
    write(41) kmax

    !-----------1=Auxiliary name/value pair to follow          !0锛氭病鏈夎緟鍔╀俊鎭紙Aux data锛?鐒跺悗鍐?357.0锛屽憡璇?Tecplot锛氬ご閮ㄧ粨鏉燂紝鎺ヤ笅鏉ユ槸鏁版嵁鍖烘弿杩?鏁版嵁
    !-----------0=No more Auxiliar name/value pairs.
    write(41) 0
    write(41) eohMarker

    !----zone ------------------------------------------------------------ !鍐嶅啓涓€娆?299.0锛歍ecplot 瑙勮寖閲屸€渮one header涔嬪悗鐨勬暟鎹弿杩板潡鈥濅細浠?zone marker 寮€濮嬨€?
    write(41) zoneMarker

    !--------variable data format, 1=Float, 2=Double, 3=LongInt,4=ShortInt, 5=Byte, 6=Bit  !姣忎釜鍙橀噺鐨勬暟鎹牸寮忥紙杩欓噷閮芥槸 float锛?鍙岀簿搴﹀氨鏄?
    write(41) 2
    write(41) 2
    write(41) 2
    write(41) 2
    write(41) 2

    !--------Has variable sharing 0 = no, 1 = yes.            !0锛氫笉鍏变韩锛堟瘡涓彉閲忛兘鍦ㄨ繖涓枃浠堕噷鐙珛瀛橈級,1锛氬叡浜紙渚嬪澶氫釜 zone 鍏变韩鍚屼竴浠?X,Y 鍧愭爣锛?
    write(41) 0

    !----------Zone number to share connectivity list with (-1 = no    !-1锛氫笉鍏变韩杩炴帴琛紙瀵?ordered zone 涓€鑸棤杩炴帴琛級
    ! sharing).
    write(41) -1

    !---------------------------------------------------------------------   !鐪熸鍐欐暟鎹紙鎸?Point packing锛?
    do k=1,kmax
        do j=1,ny
            do i=1,nx
                write(41) real(xp(i),kind=8)
                write(41) real(yp(j),kind=8)
                write(41) real(u(i,j),kind=8)
                write(41) real(v(i,j),kind=8)
                write(41) real(T(i,j),kind=8)
            end do
        end do
    enddo
    close(41)
    !---------------------------------------------------------------------

    return
  end subroutine output_Tecplot




  subroutine dumpstring(instring)
    implicit none
    character(len=40) instring
    integer(kind=4) :: stringLength   !鏈夋晥闀垮害锛堝幓鎺夊熬閮ㄧ┖鏍硷級
    integer(kind=4) :: ii             !瀛楃绱㈠紩
    integer(kind=4) :: I              !瀛楃鐨?ASCII 鐮佹暣鏁?

    stringLength=LEN_TRIM(instring)   !LEN_TRIM: 寰楀埌鍘绘帀灏鹃儴绌烘牸鍚庣殑闀垮害
    do ii=1,stringLength
        I=ICHAR(instring(ii:ii))      !ICHAR 鎶婂瓧绗﹁浆鎴?ASCII 缂栫爜鏁存暟
        write(41) I                   !鎶婅繖涓暣鏁板啓鍏ユ枃浠讹紙Tecplot 瑕佹眰浠ユ暣鏁板簭鍒楀啓瀛楃涓诧級
    end do
    write(41) 0                       !鏈€鍚庡啓涓€涓?0 浣滀负瀛楃涓茬粨鏉熺

    return
  end subroutine dumpstring


  subroutine calNuRe()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: NuVolAvg_temp    !浣撳钩鍧?Nu
    real(kind=8) :: ReVolAvg_temp    !浣撳钩鍧?Re
    
    ! 鍘熶唬鐮侊細
    ! dimensionlessTime = dimensionlessTime+1   !姣忛殧 outputFrequency 涓嚜鐢辫惤浣撴椂闂磋皟鐢ㄤ竴娆alNuRe
    if (dimensionlessTime.GE.dimensionlessTimeMax) then
        write(*,*) "Error: dimensionlessTime exceeds dimensionlessTimeMax, please enlarge dimensionlessTimeMax"
        open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
        write(00,*) "Error: dimensionlessTime exceeds dimensionlessTimeMax, please enlarge dimensionlessTimeMax"
        close(00)
        stop
    endif

    dimensionlessTime = dimensionlessTime+1   !姣忛殧 outputFrequency 涓嚜鐢辫惤浣撴椂闂磋皟鐢ㄤ竴娆alNuRe
    
    NuVolAvg_temp = 0.0d0    
    !$omp parallel do default(none) shared(v,T) private(i,j) reduction(+:NuVolAvg_temp)
    do j = 1, ny
        do i = 1, nx
            NuVolAvg_temp = NuVolAvg_temp+v(i,j)*T(i,j)     !瀵规祦鐑€氶噺
        enddo
    enddo
    !$omp end parallel do
    NuVolAvg(dimensionlessTime) = NuVolAvg_temp/dble(nx*ny)*lengthUnit/diffusivity+1.0d0    !!浣撳钩鍧?Nusselt 鏁?= 1 + (甯告暟绯绘暟) 脳 浣撳钩鍧囧娴佺儹閫氶噺

    open(unit=01,file="Nu_VolAvg.dat",status='unknown',position='append')
    write(01,*) real(reloadDimensionlessTime+dimensionlessTime*outputFrequency,kind=8), NuVolAvg(dimensionlessTime)   !浠ヨ嚜鐢辫惤浣撴椂闂存潵鍐欏叆
    close(01)

    ReVolAvg_temp = 0.0d0
    !$omp parallel do default(none) shared(u,v) private(i,j) reduction(+:ReVolAvg_temp)
    do j = 1, ny
        do i = 1, nx 
            ReVolAvg_temp = ReVolAvg_temp+(u(i,j)*u(i,j)+v(i,j)*v(i,j))
        enddo
    enddo
    !$omp end parallel do
    ReVolAvg(dimensionlessTime) = dsqrt(ReVolAvg_temp/dble(nx*ny))*lengthUnit/viscosity    !鍏ㄥ煙浣撳钩鍧?RMS-Reynolds 鏁?


    open(unit=02,file="Re_VolAvg.dat",status='unknown',position='append')
    write(02,*) real(reloadDimensionlessTime+dimensionlessTime*outputFrequency,kind=8), ReVolAvg(dimensionlessTime)  !!for print purpose only
    close(02)
    write(*,*)  ReVolAvg(dimensionlessTime)
    return
  end subroutine calNuRe
  subroutine SideHeatedcalc_Nu_global()
    use commondata
    implicit none

    call SideHeatedcalc_Nu_global_islbm()

    return
  end subroutine SideHeatedcalc_Nu_global
  subroutine SideHeatedcalc_Nu_wall_avg()
    use commondata
    implicit none

    call SideHeatedcalc_Nu_wall_avg_islbm()

    return
  end subroutine SideHeatedcalc_Nu_wall_avg
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
          ! 椤剁偣钀藉湪鎷熷悎鍖洪棿澶栵細閫€鍥炲埌 5 鐐逛腑鏈€澶х殑/鏈€灏忕殑閭ｄ釜骞朵笉涓ユ牸锛?
          ! 杩欓噷涓轰簡淇濇寔绠€鍗曪紝閫氬父浼氬厛閫夋瀬鍊肩偣
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
  subroutine SideHeatedcalc_umid_max()
    use commondata
    implicit none

    call SideHeatedcalc_umid_max_islbm()

    return
  end subroutine SideHeatedcalc_umid_max
  subroutine SideHeatedcalc_vmid_max()
    use commondata
    implicit none

    call SideHeatedcalc_vmid_max_islbm()

    return
  end subroutine SideHeatedcalc_vmid_max
subroutine RBcalc_Nu_global()
  use commondata
  implicit none

  call RBcalc_Nu_global_islbm()

  return
end subroutine RBcalc_Nu_global
subroutine RBcalc_Nu_wall_avg()
  use commondata
  implicit none

  call RBcalc_Nu_wall_avg_islbm()

  return
end subroutine RBcalc_Nu_wall_avg
subroutine RBcalc_umid_max()
  use commondata
  implicit none

  call RBcalc_umid_max_islbm()

  return
end subroutine RBcalc_umid_max
subroutine RBcalc_vmid_max()
  use commondata
  implicit none

  call RBcalc_vmid_max_islbm()

  return
end subroutine RBcalc_vmid_max
subroutine streaming_interp()
  use commondata
  implicit none
  integer(kind=4) :: i, j
  integer(kind=4) :: alpha, ii, jj
  real(kind=8) :: value

  f = 0.0d0

  !$omp parallel do default(none) shared(f,f_post,ex,ey,stream_x_valid,stream_y_valid,stream_ix,stream_iy,stream_wx,stream_wy) private(i,j,alpha,ii,jj,value)
  do j = 1, ny
    do i = 1, nx
      do alpha = 0, 8
        if (alpha == 0) then
          f(alpha,i,j) = f_post(alpha,i,j)
        elseif (ey(alpha) == 0) then
          if (stream_x_valid(alpha,i)) then
            value = 0.0d0
            do ii = 1, 3
              value = value + stream_wx(alpha,i,ii) * f_post(alpha,stream_ix(alpha,i,ii),j)
            enddo
            f(alpha,i,j) = value
          endif
        elseif (ex(alpha) == 0) then
          if (stream_y_valid(alpha,j)) then
            value = 0.0d0
            do jj = 1, 3
              value = value + stream_wy(alpha,j,jj) * f_post(alpha,i,stream_iy(alpha,j,jj))
            enddo
            f(alpha,i,j) = value
          endif
        else
          if (stream_x_valid(alpha,i) .and. stream_y_valid(alpha,j)) then
            value = 0.0d0
            do jj = 1, 3
              do ii = 1, 3
                value = value + stream_wx(alpha,i,ii) * stream_wy(alpha,j,jj) * &
                  f_post(alpha,stream_ix(alpha,i,ii),stream_iy(alpha,j,jj))
              enddo
            enddo
            f(alpha,i,j) = value
          endif
        endif
      enddo
    enddo
  enddo
  !$omp end parallel do

  return
end subroutine streaming_interp


subroutine streamingT_interp()
  use commondata
  implicit none
  integer(kind=4) :: i, j
  integer(kind=4) :: alpha, ii, jj
  real(kind=8) :: value

  g = 0.0d0

  !$omp parallel do default(none) shared(g,g_post,ex,ey,stream_x_valid,stream_y_valid,stream_ix,stream_iy,stream_wx,stream_wy) private(i,j,alpha,ii,jj,value)
  do j = 1, ny
    do i = 1, nx
      do alpha = 0, 4
        if (alpha == 0) then
          g(alpha,i,j) = g_post(alpha,i,j)
        elseif (ey(alpha) == 0) then
          if (stream_x_valid(alpha,i)) then
            value = 0.0d0
            do ii = 1, 3
              value = value + stream_wx(alpha,i,ii) * g_post(alpha,stream_ix(alpha,i,ii),j)
            enddo
            g(alpha,i,j) = value
          endif
        elseif (ex(alpha) == 0) then
          if (stream_y_valid(alpha,j)) then
            value = 0.0d0
            do jj = 1, 3
              value = value + stream_wy(alpha,j,jj) * g_post(alpha,i,stream_iy(alpha,j,jj))
            enddo
            g(alpha,i,j) = value
          endif
        else
          if (stream_x_valid(alpha,i) .and. stream_y_valid(alpha,j)) then
            value = 0.0d0
            do jj = 1, 3
              do ii = 1, 3
                value = value + stream_wx(alpha,i,ii) * stream_wy(alpha,j,jj) * &
                  g_post(alpha,stream_ix(alpha,i,ii),stream_iy(alpha,j,jj))
              enddo
            enddo
            g(alpha,i,j) = value
          endif
        endif
      enddo
    enddo
  enddo
  !$omp end parallel do

  return
end subroutine streamingT_interp


subroutine calNuRe_islbm()
  use commondata
  implicit none
  integer(kind=4) :: i, j
  real(kind=8) :: fluxVol, speedVol, areaWeight

  if (dimensionlessTime .GE. dimensionlessTimeMax) then
    write(*,*) "Error: dimensionlessTime exceeds dimensionlessTimeMax, please enlarge dimensionlessTimeMax"
    open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
    write(00,*) "Error: dimensionlessTime exceeds dimensionlessTimeMax, please enlarge dimensionlessTimeMax"
    close(00)
    stop
  endif

  dimensionlessTime = dimensionlessTime + 1
  fluxVol = 0.0d0
  speedVol = 0.0d0

  do j = 1, ny
    do i = 1, nx
      areaWeight = wx(i) * wy(j)
#ifdef SideHeatedCell
      fluxVol = fluxVol + areaWeight * u(i,j) * T(i,j)
#endif
#ifdef RayleighBenardCell
      fluxVol = fluxVol + areaWeight * v(i,j) * T(i,j)
#endif
      speedVol = speedVol + areaWeight * dsqrt(u(i,j)*u(i,j) + v(i,j)*v(i,j))
    enddo
  enddo

  NuVolAvg(dimensionlessTime) = fluxVol/quadSumArea*lengthUnit/diffusivity + 1.0d0
  ReVolAvg(dimensionlessTime) = speedVol/quadSumArea*lengthUnit/viscosity

  open(unit=01,file="Nu_VolAvg.dat",status='unknown',position='append')
  write(01,*) real(reloadDimensionlessTime+dimensionlessTime*outputFrequency,kind=8), NuVolAvg(dimensionlessTime)
  close(01)

  open(unit=02,file="Re_VolAvg.dat",status='unknown',position='append')
  write(02,*) real(reloadDimensionlessTime+dimensionlessTime*outputFrequency,kind=8), ReVolAvg(dimensionlessTime)
  close(02)

  return
end subroutine calNuRe_islbm


subroutine SideHeatedcalc_Nu_global_islbm()
  use commondata
  implicit none
  integer(kind=4) :: i, j
  real(kind=8) :: sumUT

  sumUT = 0.0d0
  do j = 1, ny
    do i = 1, nx
      sumUT = sumUT + wx(i) * wy(j) * u(i,j) * T(i,j)
    enddo
  enddo

  Nu_global = sumUT/quadSumArea*lengthUnit/diffusivity + 1.0d0

  write(*,'(a,1x,es16.8)') "Nu_global =", Nu_global
  open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
  write(00,'(a,1x,es16.8)') "Nu_global =", Nu_global
  close(00)

  return
end subroutine SideHeatedcalc_Nu_global_islbm


subroutine SideHeatedcalc_Nu_wall_avg_islbm()
  use commondata
  implicit none
  integer(kind=4) :: j, jmax, jmin, k
  integer(kind=4) :: jj(5)
  real(kind=8) :: coef, sum_hot, sum_cold, sum_mid
  real(kind=8) :: Nu_left(1:ny), Nu_right, u_mid, T_mid, dTdx_mid
  real(kind=8) :: yk(5), fk(5), fstar, ystar
  real(kind=8) :: wall_derivative_x_left, wall_derivative_x_right
  real(kind=8) :: interpolate_line_x, derivative_line_x

  coef = lengthUnit / diffusivity
  sum_hot = 0.0d0
  sum_cold = 0.0d0
  sum_mid = 0.0d0

  do j = 1, ny
    Nu_left(j) = -wall_derivative_x_left(Thot, T(1,j), T(2,j))
    Nu_right = -wall_derivative_x_right(T(nx-1,j), T(nx,j), Tcold)
    u_mid = interpolate_line_x(u(1:nx,j), 0.5d0)
    T_mid = interpolate_line_x(T(1:nx,j), 0.5d0)
    dTdx_mid = derivative_line_x(T(1:nx,j), 0.5d0)

    sum_hot = sum_hot + Nu_left(j) * wy(j)
    sum_cold = sum_cold + Nu_right * wy(j)
    sum_mid = sum_mid + (coef*u_mid*T_mid - dTdx_mid) * wy(j)
  enddo

  Nu_hot = sum_hot / quadSumY
  Nu_cold = sum_cold / quadSumY
  Nu_middle = sum_mid / quadSumY

  jmax = 1
  jmin = 1
  do j = 2, ny
    if (Nu_left(j) > Nu_left(jmax)) jmax = j
    if (Nu_left(j) < Nu_left(jmin)) jmin = j
  enddo

  if (ny >= 5) then
    if (jmax <= 2) then
      jj = (/ 1, 2, 3, 4, 5 /)
    elseif (jmax >= ny-1) then
      jj = (/ ny-4, ny-3, ny-2, ny-1, ny /)
    else
      jj = (/ jmax-2, jmax-1, jmax, jmax+1, jmax+2 /)
    endif
    do k = 1, 5
      yk(k) = yp(jj(k))
      fk(k) = Nu_left(jj(k))
    enddo
    call fit_parabola_ls5(yk, fk, +1, fstar, ystar)
    Nu_hot_max = fstar
    Nu_hot_max_position = ystar

    if (jmin <= 2) then
      jj = (/ 1, 2, 3, 4, 5 /)
    elseif (jmin >= ny-1) then
      jj = (/ ny-4, ny-3, ny-2, ny-1, ny /)
    else
      jj = (/ jmin-2, jmin-1, jmin, jmin+1, jmin+2 /)
    endif
    do k = 1, 5
      yk(k) = yp(jj(k))
      fk(k) = Nu_left(jj(k))
    enddo
    call fit_parabola_ls5(yk, fk, -1, fstar, ystar)
    Nu_hot_min = fstar
    Nu_hot_min_position = ystar
  else
    Nu_hot_max = Nu_left(jmax)
    Nu_hot_min = Nu_left(jmin)
    Nu_hot_max_position = yp(jmax)
    Nu_hot_min_position = yp(jmin)
  endif

  write(*,'(a,1x,es16.8)') "Nu_hot    =", Nu_hot
  write(*,'(a,1x,es16.8)') "Nu_cold   =", Nu_cold
  write(*,'(a,1x,es16.8)') "Nu_middle =", Nu_middle
  write(*,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_max =", Nu_hot_max, "y_max =", Nu_hot_max_position
  write(*,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_min =", Nu_hot_min, "y_min =", Nu_hot_min_position

  open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
  write(00,'(a,1x,es16.8)') "Nu_hot    =", Nu_hot
  write(00,'(a,1x,es16.8)') "Nu_cold   =", Nu_cold
  write(00,'(a,1x,es16.8)') "Nu_middle =", Nu_middle
  write(00,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_max =", Nu_hot_max, "y_max =", Nu_hot_max_position
  write(00,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_min =", Nu_hot_min, "y_min =", Nu_hot_min_position
  close(00)

  return
end subroutine SideHeatedcalc_Nu_wall_avg_islbm


subroutine SideHeatedcalc_umid_max_islbm()
  use commondata
  implicit none
  integer(kind=4) :: j, j0, k
  integer(kind=4) :: jj(5)
  real(kind=8) :: uline(1:ny)
  real(kind=8) :: s(5), fu(5)
  real(kind=8) :: umax_fit, y_fit, xmid, coef
  real(kind=8) :: interpolate_line_x

  coef = lengthUnit / diffusivity
  xmid = 0.5d0

  do j = 1, ny
    uline(j) = interpolate_line_x(u(1:nx,j), xmid)
  enddo

  j0 = 1
  do j = 2, ny
    if (uline(j) > uline(j0)) j0 = j
  enddo

  if (ny >= 5) then
    if (j0 <= 2) then
      jj = (/ 1, 2, 3, 4, 5 /)
    elseif (j0 >= ny-1) then
      jj = (/ ny-4, ny-3, ny-2, ny-1, ny /)
    else
      jj = (/ j0-2, j0-1, j0, j0+1, j0+2 /)
    endif
    do k = 1, 5
      s(k) = yp(jj(k))
      fu(k) = uline(jj(k))
    enddo
    call fit_parabola_ls5(s, fu, +1, umax_fit, y_fit)
  else
    umax_fit = uline(j0)
    y_fit = yp(j0)
  endif

  write(*,'(A,1X,F12.6,1X,A,1X,F12.6,1X,A,1X,F12.6)') &
       'u_mid_max =', umax_fit*coef, 'at y =', y_fit, 'on x_mid =', xmid

  open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
  write(00,*) 'x_mid =', xmid
  write(00,*) 'u_mid_max =', umax_fit*coef, ' y_pos =', y_fit
  close(00)

  return
end subroutine SideHeatedcalc_umid_max_islbm


subroutine SideHeatedcalc_vmid_max_islbm()
  use commondata
  implicit none
  integer(kind=4) :: i, j, i0, k
  integer(kind=4) :: ii(5)
  real(kind=8) :: vline(1:nx), col(1:ny)
  real(kind=8) :: s(5), fv(5)
  real(kind=8) :: vmax_fit, x_fit, ymid, coef
  real(kind=8) :: interpolate_line_y

  coef = lengthUnit / diffusivity
  ymid = 0.5d0

  do i = 1, nx
    do j = 1, ny
      col(j) = v(i,j)
    enddo
    vline(i) = interpolate_line_y(col, ymid)
  enddo

  i0 = 1
  do i = 2, nx
    if (vline(i) > vline(i0)) i0 = i
  enddo

  if (nx >= 5) then
    if (i0 <= 2) then
      ii = (/ 1, 2, 3, 4, 5 /)
    elseif (i0 >= nx-1) then
      ii = (/ nx-4, nx-3, nx-2, nx-1, nx /)
    else
      ii = (/ i0-2, i0-1, i0, i0+1, i0+2 /)
    endif
    do k = 1, 5
      s(k) = xp(ii(k))
      fv(k) = vline(ii(k))
    enddo
    call fit_parabola_ls5(s, fv, +1, vmax_fit, x_fit)
  else
    vmax_fit = vline(i0)
    x_fit = xp(i0)
  endif

  write(*,'(A,1X,F12.6,1X,A,1X,F12.6,1X,A,1X,F12.6)') &
       'v_mid_max =', vmax_fit*coef, 'at x =', x_fit, 'on y_mid =', ymid

  open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
  write(00,*) 'y_mid =', ymid
  write(00,*) 'v_mid_max =', vmax_fit*coef, ' x_pos =', x_fit
  close(00)

  return
end subroutine SideHeatedcalc_vmid_max_islbm


subroutine RBcalc_Nu_global_islbm()
  use commondata
  implicit none
  integer(kind=4) :: i, j
  real(kind=8) :: sumVT

  sumVT = 0.0d0
  do j = 1, ny
    do i = 1, nx
      sumVT = sumVT + wx(i) * wy(j) * v(i,j) * T(i,j)
    enddo
  enddo

  Nu_global = sumVT/quadSumArea*lengthUnit/diffusivity + 1.0d0

  write(*,'(a,1x,es16.8)') "Nu_global =", Nu_global
  open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
  write(00,'(a,1x,es16.8)') "Nu_global =", Nu_global
  close(00)

  return
end subroutine RBcalc_Nu_global_islbm


subroutine RBcalc_Nu_wall_avg_islbm()
  use commondata
  implicit none
  integer(kind=4) :: i, imax, imin, k
  integer(kind=4) :: ii(5)
  real(kind=8) :: sum_hot, sum_cold
  real(kind=8) :: Nu_bot(1:nx), Nu_top
  real(kind=8) :: xk(5), fk(5), fstar, xstar
  real(kind=8) :: wall_derivative_y_bottom, wall_derivative_y_top

  sum_hot = 0.0d0
  sum_cold = 0.0d0

  do i = 1, nx
    Nu_bot(i) = -wall_derivative_y_bottom(Thot, T(i,1), T(i,2))
    Nu_top = -wall_derivative_y_top(T(i,ny-1), T(i,ny), Tcold)
    sum_hot = sum_hot + Nu_bot(i) * wx(i)
    sum_cold = sum_cold + Nu_top * wx(i)
  enddo

  Nu_hot = sum_hot / quadSumX
  Nu_cold = sum_cold / quadSumX

  imax = 1
  imin = 1
  do i = 2, nx
    if (Nu_bot(i) > Nu_bot(imax)) imax = i
    if (Nu_bot(i) < Nu_bot(imin)) imin = i
  enddo

  if (nx >= 5) then
    if (imax <= 2) then
      ii = (/ 1, 2, 3, 4, 5 /)
    elseif (imax >= nx-1) then
      ii = (/ nx-4, nx-3, nx-2, nx-1, nx /)
    else
      ii = (/ imax-2, imax-1, imax, imax+1, imax+2 /)
    endif
    do k = 1, 5
      xk(k) = xp(ii(k))
      fk(k) = Nu_bot(ii(k))
    enddo
    call fit_parabola_ls5(xk, fk, +1, fstar, xstar)
    Nu_hot_max = fstar
    Nu_hot_max_position = xstar

    if (imin <= 2) then
      ii = (/ 1, 2, 3, 4, 5 /)
    elseif (imin >= nx-1) then
      ii = (/ nx-4, nx-3, nx-2, nx-1, nx /)
    else
      ii = (/ imin-2, imin-1, imin, imin+1, imin+2 /)
    endif
    do k = 1, 5
      xk(k) = xp(ii(k))
      fk(k) = Nu_bot(ii(k))
    enddo
    call fit_parabola_ls5(xk, fk, -1, fstar, xstar)
    Nu_hot_min = fstar
    Nu_hot_min_position = xstar
  else
    Nu_hot_max = Nu_bot(imax)
    Nu_hot_min = Nu_bot(imin)
    Nu_hot_max_position = xp(imax)
    Nu_hot_min_position = xp(imin)
  endif

  write(*,'(a,1x,es16.8)') "Nu_hot(bottom) =", Nu_hot
  write(*,'(a,1x,es16.8)') "Nu_cold(top)   =", Nu_cold
  write(*,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_max =", Nu_hot_max, "x_max =", Nu_hot_max_position
  write(*,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_min =", Nu_hot_min, "x_min =", Nu_hot_min_position

  open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
  write(00,'(a,1x,es16.8)') "Nu_hot(bottom) =", Nu_hot
  write(00,'(a,1x,es16.8)') "Nu_cold(top)   =", Nu_cold
  write(00,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_max =", Nu_hot_max, "x_max =", Nu_hot_max_position
  write(00,'(a,1x,es16.8,2x,a,1x,es16.8)') "Nu_hot_min =", Nu_hot_min, "x_min =", Nu_hot_min_position
  close(00)

  return
end subroutine RBcalc_Nu_wall_avg_islbm


subroutine RBcalc_umid_max_islbm()
  use commondata
  implicit none
  integer(kind=4) :: i, j, i0, k
  integer(kind=4) :: ii(5)
  real(kind=8) :: uline(1:nx), col(1:ny)
  real(kind=8) :: xk(5), fk(5)
  real(kind=8) :: umax_fit, x_fit, coef, ymid
  real(kind=8) :: interpolate_line_y

  coef = lengthUnit / diffusivity
  ymid = 0.5d0

  do i = 1, nx
    do j = 1, ny
      col(j) = u(i,j)
    enddo
    uline(i) = interpolate_line_y(col, ymid)
  enddo

  i0 = 1
  do i = 2, nx
    if (uline(i) > uline(i0)) i0 = i
  enddo

  if (nx >= 5) then
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
  else
    umax_fit = uline(i0)
    x_fit = xp(i0)
  endif

  write(*,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'u_mid_abs_max* =', umax_fit*coef, 'x =', x_fit, 'y_mid =', ymid

  open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
  write(00,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'u_mid_abs_max* =', umax_fit*coef, 'x =', x_fit, 'y_mid =', ymid
  close(00)

  return
end subroutine RBcalc_umid_max_islbm


subroutine RBcalc_vmid_max_islbm()
  use commondata
  implicit none
  integer(kind=4) :: j, j0, k
  integer(kind=4) :: jj(5)
  real(kind=8) :: vline(1:ny)
  real(kind=8) :: yk(5), fk(5)
  real(kind=8) :: vmax_fit, y_fit, coef, xmid
  real(kind=8) :: interpolate_line_x

  coef = lengthUnit / diffusivity
  xmid = 0.5d0

  do j = 1, ny
    vline(j) = interpolate_line_x(v(1:nx,j), xmid)
  enddo

  j0 = 1
  do j = 2, ny
    if (vline(j) > vline(j0)) j0 = j
  enddo

  if (ny >= 5) then
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
  else
    vmax_fit = vline(j0)
    y_fit = yp(j0)
  endif

  write(*,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'v_mid_abs_max* =', vmax_fit*coef, 'y =', y_fit, 'x_mid =', xmid

  open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
  write(00,'(A,1X,ES16.8,2X,A,1X,ES16.8,2X,A,1X,ES16.8)') &
       'v_mid_abs_max* =', vmax_fit*coef, 'y =', y_fit, 'x_mid =', xmid
  close(00)

  return
end subroutine RBcalc_vmid_max_islbm


subroutine output_centerline_profiles()
  use commondata
  implicit none
  integer(kind=4) :: i, j
  real(kind=8) :: targetX, targetY, coef, lineValue
  real(kind=8) :: col(1:ny)
  real(kind=8) :: interpolate_line_x, interpolate_line_y

  targetX = 0.5d0
  targetY = 0.5d0
  coef = lengthUnit / diffusivity

  open(unit=21,file="centerline_u_y.dat",status='replace')
  do j = 1, ny
    lineValue = interpolate_line_x(u(1:nx,j), targetX)
    write(21,'(2(1x,es24.16))') yp(j), lineValue*coef
  enddo
  close(21)

  open(unit=22,file="centerline_v_x.dat",status='replace')
  do i = 1, nx
    do j = 1, ny
      col(j) = v(i,j)
    enddo
    lineValue = interpolate_line_y(col, targetY)
    write(22,'(2(1x,es24.16))') xp(i), lineValue*coef
  enddo
  close(22)

  open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
  write(00,*) "Centerline profiles written at x=0.5 and y=0.5"
  close(00)

  return
end subroutine output_centerline_profiles


subroutine calc_psi_vort_and_output_safe()
  use commondata
  implicit none

  if (meshMode /= meshModeUniform) then
    write(*,*) "psi/vorticity output is skipped on non-uniform meshes."
    open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
    write(00,*) "psi/vorticity output is skipped on non-uniform meshes."
    close(00)
    return
  endif

  call calc_psi_vort_and_output()

  return
end subroutine calc_psi_vort_and_output_safe


subroutine calc_psi_vort_and_output()
  use commondata
  implicit none

  integer(kind=4) :: i, j, m
  real(kind=8) :: dx, dy, coef
  real(kind=8) :: u1, u2, u_mid, inc
  real(kind=8) :: dv_dx, du_dy
  real(kind=8) :: psi(nx,ny), vort(nx,ny)


  ! for fine-grid max(|psi|)  10001*10001
  real(kind=8) :: psi_abs_max, x_at_max, y_at_max

  call output_centerline_profiles()
  if (meshMode /= meshModeUniform) then
    write(*,*) "psi/vorticity output is skipped on non-uniform meshes."
    open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
    write(00,*) "psi/vorticity output is skipped on non-uniform meshes."
    close(00)
    return
  endif

  ! 鍘熶唬鐮侊細
  ! dx   = 1.0d0 / nx
  dx   = 1.0d0 / lengthUnit
  dy   = 1.0d0 / lengthUnit
  coef = lengthUnit / diffusivity   ! L/kappa锛屽搴旀枃鐚爣搴?
  !=========================================================
  ! (A) 璁＄畻娴佸嚱鏁?psi锛? psi(x,y)=鈭玙0^y u(x,mu)dmu
  !     绉垎鐢ㄢ€滅疮绉?Simpson鈥濓紝涓偣 u(y=鏁寸偣) 鐢ㄤ簩闃跺椤瑰紡鎻掑€?
  !=========================================================

  do i = 1, nx

    ! ---- (A1) 绗竴涓崐鏍肩偣 y=dy/2 鐨勭Н鍒嗭細鐢ㄤ簩娆″椤瑰紡骞跺己鍒跺闈?u=0锛坣o-slip锛?
    ! 杩欓噷绗竴灏忔 [0,dy/2] 鍗曠嫭澶勭悊锛屾暣浣撻樁鏁扮敱瀹冩帶鍒?
    u1 = u(i,1) * coef
    u2 = u(i,2) * coef
    psi(i,1) = dy * (21.0d0*u1 - u2) / 72.0d0   !绗竴涓偣鐨刾si

    ! ---- (A2) 浠?y=(j-1/2)dy 鍒?y=(j+1/2)dy 鐨勬瘡涓€娈甸暱搴︿负 dy锛岀敤 Simpson锛?
    ! 鈭玙{y_{j-1/2}}^{y_{j+1/2}} u dy 鈮?dy/6 * [ u_{j-1/2} + 4 u_{j} + u_{j+1/2} ]
    do j = 2, ny
      m = j - 1   ! 杩欎竴娈电殑涓偣鍦?y = m*dy锛堟暣鐐癸級

      ! ---- (A2.1) 璁＄畻涓偣閫熷害 u(y=m*dy)锛?
      ! 鐢ㄤ簩闃?Lagrange 鎻掑€硷紙涓夌偣锛夛紝鎶婂崐鏍肩偣鍊兼彃鍒版暣鐐?
      !
      ! 瀵?m>=2锛氱敤 (m-1/2, m+1/2) 鍙婃洿涓嬫柟鐨?(m-3/2) 涓夌偣 => 绯绘暟(-1/8, 3/4, 3/8)
      ! 瀵?m=1锛堥潬杩戝簳澹侊級锛氬彧鑳界敤鏈€闈犺繎搴曢儴鐨勪笁鐐?=> 绯绘暟( 3/8, 3/4,-1/8)
      if (m == 1) then
        u_mid = ( 3.0d0/8.0d0*u(i,1) + 3.0d0/4.0d0*u(i,2) - 1.0d0/8.0d0*u(i,3) ) * coef
      else
        u_mid = ( -1.0d0/8.0d0*u(i,m-1) + 3.0d0/4.0d0*u(i,m) + 3.0d0/8.0d0*u(i,m+1) ) * coef
      end if

      inc = dy/6.0d0 * ( (u(i,j-1)*coef) + 4.0d0*u_mid + (u(i,j)*coef) )
      psi(i,j) = psi(i,j-1) + inc
    end do

  end do


  call output_psi_center_abs(psi)     !杈撳嚭涓績鐨刟bs(psi)

  !=========================================================
  ! (B) 璁＄畻娑￠噺 vort = dv/dx - du/dy锛?D锛?
  !     鍐呴儴鐢ㄤ腑蹇冨樊鍒嗭紱杈圭晫鑺傜偣鐢ㄢ€滃鍦ㄥ崐鏍艰窛鈥濆亣璁句笅鐨勪簩闃跺崟杈瑰叕寮?
  !=========================================================

  do j = 1, ny
    do i = 1, nx

      ! ---- dv/dx
      if (i == 1) then
        ! 鍦?x=dx/2 澶勶紝鐢?(wall, i=1, i=2) 浜屾鎷熷悎寰楀埌浜岄樁杩戜技
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

      ! psi 鐢ㄤ簡L/kappa锛泇ort 涔熷簲鍦ㄥ悓涓€鏍囧害涓?
      ! vort = d(v*coef)/dx - d(u*coef)/dy = coef*(dv/dx - du/dy)
      vort(i,j) = coef * (dv_dx - du_dy)

    end do
  end do

  !=========================================================
  ! (C) 杈撳嚭 Tecplot锛歑,Y,psi,vort
  !=========================================================
  call output_Tecplot_psi_vort(psi, vort)

  !=========================================================
  ! (D) 缁嗙綉鏍?10001脳10001锛氱敤涓夋鏍锋潯鎻掑€煎鎵?max(|psi|)
  !=========================================================
  call calc_psi_absmax_fine_spline(psi, psi_abs_max, x_at_max, y_at_max)

  write(*,'(a,1x,es16.8,2x,a,1x,es16.8,2x,a,1x,es16.8)') &
       "max(|psi|) =", psi_abs_max, "x* =", x_at_max, "y* =", y_at_max

  open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
  write(00,'(a,1x,es16.8,2x,a,1x,es16.8,2x,a,1x,es16.8)') &
       "max(|psi|) =", psi_abs_max, "x* =", x_at_max, "y* =", y_at_max
  close(00)

  return
end subroutine calc_psi_vort_and_output



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


subroutine calc_psi_absmax_fine_spline(psi, psi_abs_max, x_at_max, y_at_max)
  use commondata
  implicit none

  real(kind=8), intent(in)  :: psi(nx,ny)
  real(kind=8), intent(out) :: psi_abs_max, x_at_max, y_at_max

  integer(kind=4), parameter :: nFine = 10001    !缁嗙綉鏍?
  integer(kind=4) :: i, j, k, l
  integer(kind=4) :: nXExt, nYExt

  real(kind=8), allocatable :: xFine(:), yFine(:)
  real(kind=8), allocatable :: xExt(:), yExt(:)
  real(kind=8), allocatable :: row(:), y2row(:)
  real(kind=8), allocatable :: col(:), y2col(:)
  real(kind=8), allocatable :: psi_xfine(:,:)   ! (nFine, ny)

  real(kind=8) :: val, xq, yq

  ! ---- 缁嗙綉鏍煎潗鏍囷細0..1 绛夎窛
  allocate(xFine(nFine), yFine(nFine))
  do k = 1, nFine
    xFine(k) = dble(k-1) / dble(nFine-1)   !绛夊垎鎴?10000 娈?
    yFine(k) = dble(k-1) / dble(nFine-1)
  end do

  ! ---- 涓轰簡鑳藉湪 x=0/1 涓?y=0/1 涓婃彃鍊硷紝鐗╃悊杈圭晫 psi=甯告暟锛堝彲鍙?锛夎ˉ涓や釜绔偣
  nXExt = nx + 2
  nYExt = ny + 2
  ! 鍘熶唬鐮侊細
  ! allocate(xExt(nYExt), yExt(nYExt))
  allocate(xExt(nXExt), yExt(nYExt))
  xExt(1)    = 0.0d0
  ! 鍘熶唬鐮侊細
  ! xExt(nX)   = 1.0d0
  xExt(nXExt)= 1.0d0
  do i = 1, nx
    xExt(i+1) = xp(i)         !xExt = [0, xp(1),...,xp(nx), 1]
  end do

  yExt(1)    = 0.0d0
  ! 鍘熶唬鐮侊細
  ! yExt(nY)   = 1.0d0
  yExt(nYExt)= 1.0d0
  do j = 1, ny
    yExt(j+1) = yp(j)        !yExt = [0, yp(1),...,yp(ny), 1]
  end do

  allocate(row(nXExt), y2row(nXExt))
  allocate(psi_xfine(nFine, ny))

  ! ---- (1) 鍏堝姣忎釜鍥哄畾 y=yp(j) 鐨勫墫闈㈠仛 x 鏂瑰悜涓夋鏍锋潯锛屽緱鍒?psi(xFine, yp(j))
  do j = 1, ny
    row(1)  = 0.0d0
    row(nXExt) = 0.0d0
    do i = 1, nx
      row(i+1) = psi(i,j)     !row=[0, psi(1),...,psi(ny), 0]
    end do

    call spline_natural(nXExt, xExt, row, y2row)   !鍦╮ow鍚勮妭鐐瑰鐨勪簩闃跺

    do k = 1, nFine
      xq = xFine(k)
      call splint(nXExt, xExt, row, y2row, xq, val)   !鍦?10001 涓粏缃戞牸 x 鐐逛笂閲囨牱
      psi_xfine(k,j) = val                            !鍦ㄦ瘡涓矖 y 灞備笂锛宲si宸茬粡娌?x 琚粏鍖栧埌 10001 涓偣銆?
    end do
  end do

  ! ---- (2) 鍐嶅姣忎釜鍥哄畾 x=xFine(k) 鐨勫墫闈㈠仛 y 鏂瑰悜涓夋鏍锋潯锛屽湪 10001 涓?yFine 涓婃壂 max(|psi|)
  allocate(col(nYExt), y2col(nYExt))
  psi_abs_max = -1.0d0
  x_at_max    = 0.0d0
  y_at_max    = 0.0d0

  do k = 1, nFine
    col(1)  = 0.0d0
    col(nYExt) = 0.0d0
    do j = 1, ny
      col(j+1) = psi_xfine(k,j)        !瀵规瘡涓€涓浐瀹氱殑缁嗙綉鏍?x=xFine(k)锛屾湁涓€鏉℃部 y 鐨勭鏁ｅ墫闈㈡暟鎹?col
    end do

    call spline_natural(nYExt, yExt, col, y2col)

    do l = 1, nFine
      yq = yFine(l)
      call splint(nYExt, yExt, col, y2col, yq, val)

      if (dabs(val) > psi_abs_max) then     !瀵绘壘鏈€澶х殑abs(psi)浠ュ強浣嶇疆
        psi_abs_max = dabs(val)
        x_at_max = xFine(k)
        y_at_max = yFine(l)
      end if
    end do
  end do

  deallocate(xFine, yFine, xExt, yExt, row, y2row, col, y2col, psi_xfine)
  return
end subroutine calc_psi_absmax_fine_spline


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
  if (dabs(h) <= 1.0d-15) then
    y = ya(klo)
    return
  end if

  a = (xa(khi) - x)/h
  b = (x - xa(klo))/h

  y = a*ya(klo) + b*ya(khi) + ( (a*a*a - a)*y2a(klo) + (b*b*b - b)*y2a(khi) ) * (h*h)/6.0d0
  return
end subroutine splint


subroutine output_psi_center_abs(psi)
  use commondata
  implicit none
  real(kind=8), intent(in) :: psi(nx,ny)

  integer(kind=4) :: iC, jC
  integer(kind=4) :: iL, iR, jB, jT
  real(kind=8) :: psi_center, psi_center_abs

  !--------------------------------------------
  ! Center value at (x,y) = (0.5,0.5) on coarse grid
  ! xp(i)=(i-0.5)/lengthUnit, yp(j)=(j-0.5)/lengthUnit
  !--------------------------------------------
  if (mod(nx,2) == 1 .and. mod(ny,2) == 1) then
    ! Odd-odd: center is exactly on a node
    iC = (nx + 1)/2
    jC = (ny + 1)/2
    psi_center = psi(iC,jC)

  elseif (mod(nx,2) == 0 .and. mod(ny,2) == 0) then
    ! Even-even: center is at the middle of a cell (4-point average)
    iL = ishft(nx,-1)
    iR = iL + 1
    jB = ishft(ny,-1)
    jT = jB + 1
    psi_center = 0.25d0*( psi(iL,jB) + psi(iR,jB) + psi(iL,jT) + psi(iR,jT) )

  elseif (mod(nx,2) == 0 .and. mod(ny,2) == 1) then
    ! Even nx, odd ny: average in x at centerline row
    iL = ishft(nx,-1)
    iR = iL + 1
    jC = (ny + 1)/2
    psi_center = 0.5d0*( psi(iL,jC) + psi(iR,jC) )

  else
    ! Odd nx, even ny: average in y at centerline column
    iC = (nx + 1)/2
    jB = ishft(ny,-1)
    jT = jB + 1
    psi_center = 0.5d0*( psi(iC,jB) + psi(iC,jT) )
  endif

  psi_center_abs = dabs(psi_center)

  ! Screen output
  write(*,'(a,1x,es16.8)') "abs(psi_center) =", psi_center_abs

  ! Log output
  open(unit=00,file="SimulationSettings.txt",status="unknown",position="append")
  write(00,'(a,1x,es16.8)') "abs(psi_center) =", psi_center_abs
  close(00)

  return
end subroutine output_psi_center_abs


