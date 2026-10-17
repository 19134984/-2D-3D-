! 独立后处理程序：CPU运行，不依赖OpenACC或求解器的commondata模块。
! 修改下面参数区即可重新统计，不必重新运行求解器。
module postdata
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    ! 0: 非稳态窗口统计；1: 稳态最终结果与误差曲线，不做时间平均。
    integer, parameter :: postprocessSteady=0
    logical, parameter :: processSnapshots=.true.
    real(kind=8), parameter :: unsteadyAverageStartTf=500.0d0
    real(kind=8), parameter :: unsteadyAverageEndTf=1000.0d0
    real(kind=8), parameter :: unsteadyAverageMidTf= &
        0.5d0*(unsteadyAverageStartTf+unsteadyAverageEndTf)
    character(*), parameter :: NuReHistoryFile='NuRe_2DOpenaccMultiblock.dat'
    character(*), parameter :: outputIndexFile='OutputFiles_2DOpenaccMultiblock.dat'
    character(*), parameter :: snapshotFilePrefix='buoyancyCavity2DOpenaccMultiblockSnapshot'
    ! 下列参数和几何由快照读入，不手工重复填写物理参数。
    integer :: nx,ny,itc,refineRatio
    integer, parameter :: reloadItcOffset=0, interfaceSkin=0
    real(kind=8), parameter :: dxFine=1d0
    real(kind=8) :: lengthUnit,timeUnit,viscosity,diffusivity,Thot,Tcold,Tref,dxCoarse
    real(kind=8) :: centerIntegrationBox(4)
    logical :: sideHeated
    real(kind=8) :: xOffsetCoarse, yOffsetCoarse
    real(kind=8), allocatable :: u_coarse(:,:),v_coarse(:,:),T_coarse(:,:),rho_coarse(:,:)
    real(kind=8) :: xOffsetLeft, yOffsetLeft
    real(kind=8), allocatable :: u_left(:,:),v_left(:,:),T_left(:,:),rho_left(:,:)
    real(kind=8) :: xOffsetRight, yOffsetRight
    real(kind=8), allocatable :: u_right(:,:),v_right(:,:),T_right(:,:),rho_right(:,:)
    real(kind=8) :: xOffsetBottom, yOffsetBottom
    real(kind=8), allocatable :: u_bottom(:,:),v_bottom(:,:),T_bottom(:,:),rho_bottom(:,:)
    real(kind=8) :: xOffsetTop, yOffsetTop
    real(kind=8), allocatable :: u_top(:,:),v_top(:,:),T_top(:,:),rho_top(:,:)
end module postdata

program postprocess
    use postdata
    implicit none
    integer :: k,u,ios,n,sample,step,prevStep,interval
    real(kind=8) :: values(12),previousTime,dt,firstTime,lastTime
    real(kind=8) :: allMean(5),firstMean(5),lastMean(5),c0,c1,c2
    character(len=512) :: line
    character(len=256) :: name

    ! 首先核对连续采样编号、整数步数与实际时间；不把缺测间隙跨过去积分。
    open(newunit=k,file=NuReHistoryFile,status='old')
    open(newunit=u,file='NuReCurve_2DOpenaccMultiblock.dat',status='replace')
    write(u,'(a)') '# sample itc_total t_ff NuVolAvg ReVolRMS Nu_hot Nu_cold Nu_middle mass meanT Tmin Tmax rhoMin rhoMax'
    n=0
    prevStep=0
    previousTime=0d0
    interval=0
    do
        read(k,'(a)',iostat=ios) line
        if(ios<0) exit
        if(ios>0) error stop 'History read error'
        if(len_trim(line)==0 .OR. line(1:1)=='#') cycle
        read(line,*,iostat=ios) sample,step,values
        if(ios/=0 .OR. .NOT.all(ieee_is_finite(values))) error stop 'Invalid Nu/Re history'
        n=n+1
        if(n==1) then
            interval=step
            dt=values(1)
            firstTime=values(1)
            if(interval<=0 .OR. dt<=0d0) error stop 'Invalid first sample'
        endif
        if(sample/=n .OR. step/=n*interval) error stop 'Missing or duplicate sample'
        if(values(1)<=previousTime .OR. abs(values(1)-dble(n)*dt)>1d-10*max(1d0,abs(values(1)))) &
            error stop 'History time does not match fixed sampling'
        previousTime=values(1)
        prevStep=step
        write(u,'(a)') trim(line)
    enddo
    close(k)
    close(u)
    if(n==0 .AND. postprocessSteady==0) error stop 'No samples for unsteady statistics'
    lastTime=previousTime

    if(postprocessSteady==1) then
        ! 每次check记录的误差直接用于画图，无需从快照重新定义误差。
        open(newunit=k,file='Convergence_2DOpenaccMultiblock.dat',status='old')
        open(newunit=u,file='ConvergenceCurve_2DOpenaccMultiblock.dat',status='replace')
        write(u,'(a)') '# itc_total t_ff errorU errorT'
        prevStep=-1
        do
            read(k,'(a)',iostat=ios) line
            if(ios<0) exit
            if(ios>0) error stop 'Convergence read error'
            if(len_trim(line)==0 .OR. line(1:1)=='#') cycle
            read(line,*,iostat=ios) step,values(1:3)
            if(ios/=0 .OR. .NOT.all(ieee_is_finite(values(1:3)))) error stop 'Invalid convergence data'
            if(step<=prevStep) error stop 'Nonmonotone convergence steps'
            prevStep=step
            write(u,'(a)') trim(line)
        enddo
        close(k)
        close(u)
        ! 原样保留最终时刻与converged标志；不把最后一个周期样本当作最终状态。
        open(newunit=k,file='FinalNuRe_2DOpenaccMultiblock.dat',status='old')
        open(newunit=u,file='SteadySummary_2DOpenaccMultiblock.dat',status='replace')
        do
            read(k,'(a)',iostat=ios) line
            if(ios<0) exit
            if(ios>0) error stop 'Final result read error'
            write(u,'(a)') trim(line)
        enddo
        close(k)
        close(u)
    else
        if(unsteadyAverageStartTf<0d0 .OR. unsteadyAverageStartTf>=unsteadyAverageMidTf .OR. &
            unsteadyAverageMidTf>=unsteadyAverageEndTf) error stop 'Invalid averaging window'
        open(newunit=u,file='NuReStatistics_2DOpenaccMultiblock.dat',status='replace')
        write(u,*) 'Window t_ff:',unsteadyAverageStartTf,unsteadyAverageEndTf
        if(firstTime>unsteadyAverageStartTf .OR. lastTime<unsteadyAverageEndTf) then
            write(u,'(a)') 'INCOMPLETE: the history does not cover the requested window; no final mean.'
        else
            call average_window(unsteadyAverageStartTf,unsteadyAverageEndTf,allMean,c0)
            call average_window(unsteadyAverageStartTf,unsteadyAverageMidTf,firstMean,c1)
            call average_window(unsteadyAverageMidTf,unsteadyAverageEndTf,lastMean,c2)
            write(u,*) 'Covered durations:',c0,c1,c2
            write(u,'(a)') '# NuVolAvg Re_rms_space_time Nu_hot Nu_cold Nu_middle'
            write(u,'(a,5ES24.16E3)') 'whole:',allMean
            write(u,'(a,5ES24.16E3)') 'first:',firstMean
            write(u,'(a,5ES24.16E3)') 'last :',lastMean
            write(u,'(a,5ES24.16E3)') 'relative half-window difference:', &
                abs(lastMean-firstMean)/max(abs(allMean),1d-30)
        endif
        close(u)
    endif

    if(processSnapshots) then
        open(newunit=k,file=outputIndexFile,status='old')
        open(newunit=u,file='SnapshotNuRe_2DOpenaccMultiblock.dat',status='replace')
        write(u,'(a)') '# itc_total t_ff NuVolAvg ReVolRMS Nu_hot Nu_cold Nu_middle mass meanT Tmin Tmax rhoMin rhoMax'
        do
            read(k,'(a)',iostat=ios) line
            if(ios<0) exit
            if(ios>0) error stop 'Output index read error'
            if(len_trim(line)==0) cycle
            read(line,*,iostat=ios) step,name
            if(ios/=0) error stop 'Invalid output index'
            if(index(name,snapshotFilePrefix//'-')/=1) cycle
            ! 索引保存完整的<tff>tff-<累计细步数>itc.bin文件名，直接打开，不从文件名反推时间或步数。
            call read_snapshot(trim(name))
            if(itc/=step) error stop 'Snapshot/index step mismatch'
            call calculate_NuRe(values(2),values(3),values(4),values(5),values(6),values(7), &
                values(8),values(9),values(10),values(11),values(12))
            values(1)=dble(itc)/timeUnit
            write(u,'(I0,1X,12(ES24.16E3,1X))') itc,values
        enddo
        close(k)
        close(u)
    endif
    write(*,'(a)') 'Postprocessing complete; solver fields and histories were not modified.'
end program postprocess

subroutine average_window(t0,t1,result,coverage)
    use postdata,only: NuReHistoryFile
    implicit none
    real(kind=8),intent(in) :: t0,t1
    real(kind=8),intent(out) :: result(5),coverage
    real(kind=8) :: prev(12),row(12),a,b,left(5),right(5),ra(5),rb(5)
    integer :: k,ios,sample,step
    logical :: hasPrev
    character(len=512) :: line
    result=0d0
    coverage=0d0
    hasPrev=.false.
    open(newunit=k,file=NuReHistoryFile,status='old')
    do
        read(k,'(a)',iostat=ios) line
        if(ios<0) exit
        if(ios>0) error stop 'Averaging read error'
        if(len_trim(line)==0 .OR. line(1:1)=='#') cycle
        read(line,*,iostat=ios) sample,step,row
        if(ios/=0) error stop 'Invalid averaging row'
        if(hasPrev) then
            a=max(t0,prev(1))
            b=min(t1,row(1))
            if(b>a) then
                left=prev(2:6)
                right=row(2:6)
                ! 保持时空均方根定义：对瞬时空间RMS的平方积分，再开方。
                left(2)=left(2)**2
                right(2)=right(2)**2
                ra=left+(right-left)*(a-prev(1))/(row(1)-prev(1))
                rb=left+(right-left)*(b-prev(1))/(row(1)-prev(1))
                result=result+0.5d0*(ra+rb)*(b-a)
                coverage=coverage+b-a
            endif
        endif
        prev=row
        hasPrev=.true.
    enddo
    close(k)
    if(abs(coverage-(t1-t0))>1d-10*max(1d0,t1-t0)) error stop 'Incomplete average window'
    result=result/coverage
    result(2)=sqrt(max(0d0,result(2)))
end subroutine average_window

subroutine read_snapshot(name)
    use postdata
    implicit none
    character(*),intent(in) :: name
    character(len=16) :: magic
    integer :: k,ios,blocks,caseMode,r,ni(5),nj(5),first,last,count
    real(kind=8) :: xfirst(5),yfirst(5),dx(5),savedTime
    real(kind=8),allocatable :: uall(:),vall(:),Tall(:),rhoall(:)
    open(newunit=k,file=name,form='unformatted',access='stream',status='old')
    read(k,iostat=ios) magic,blocks,nx,ny,itc,caseMode
    if(ios/=0 .OR. magic/='MB2DSNAPSHOT0001') error stop 'Invalid snapshot header'
    if((blocks/=1 .AND. blocks/=5) .OR. min(nx,ny)<4) error stop 'Invalid snapshot mesh'
    if(caseMode/=0 .AND. caseMode/=1) error stop 'Invalid snapshot case'
    sideHeated=caseMode==1
    read(k,iostat=ios) savedTime,lengthUnit,viscosity,diffusivity,Thot,Tcold,Tref,timeUnit,centerIntegrationBox
    if(ios/=0) error stop 'Incomplete snapshot physics'
    if(.NOT.all(ieee_is_finite(centerIntegrationBox))) error stop 'Invalid snapshot integration box'
    if(.NOT.all(ieee_is_finite([savedTime,lengthUnit,viscosity,diffusivity,timeUnit,Thot,Tcold,Tref]))) &
        error stop 'Nonfinite snapshot metadata'
    if(min(lengthUnit,viscosity,diffusivity,timeUnit)<=0d0 .OR. Thot==Tcold) error stop 'Invalid snapshot units'
    if(abs(savedTime-dble(itc)/timeUnit)>1d-12*max(1d0,abs(savedTime))) error stop 'Snapshot time mismatch'
    count=0
    do r=1,blocks
        read(k,iostat=ios) ni(r),nj(r),xfirst(r),yfirst(r),dx(r)
        if(ios/=0 .OR. min(ni(r),nj(r))<1 .OR. dx(r)<=0d0) error stop 'Invalid snapshot block'
        if(.NOT.all(ieee_is_finite([xfirst(r),yfirst(r),dx(r)]))) error stop 'Nonfinite block geometry'
        count=count+ni(r)*nj(r)
    enddo
    dxCoarse=dx(1)
    refineRatio=nint(dxCoarse)
    allocate(uall(count),vall(count),Tall(count),rhoall(count))
    read(k,iostat=ios) uall,vall,Tall,rhoall
    if(ios/=0) error stop 'Incomplete snapshot fields'
    close(k)
    first=1
    last=0
    if(allocated(T_coarse)) deallocate(u_coarse,v_coarse,T_coarse,rho_coarse)
    allocate(u_coarse(ni(1),nj(1)),v_coarse(ni(1),nj(1)), &
        T_coarse(ni(1),nj(1)),rho_coarse(ni(1),nj(1)))
    first=last+1
    last=last+ni(1)*nj(1)
    xOffsetCoarse=xfirst(1)-0.5d0*dx(1)
    yOffsetCoarse=yfirst(1)-0.5d0*dx(1)
    u_coarse=reshape(uall(first:last),[ni(1),nj(1)])
    v_coarse=reshape(vall(first:last),[ni(1),nj(1)])
    T_coarse=reshape(Tall(first:last),[ni(1),nj(1)])
    rho_coarse=reshape(rhoall(first:last),[ni(1),nj(1)])
    if(blocks==5) then
    if(allocated(T_left)) deallocate(u_left,v_left,T_left,rho_left)
    allocate(u_left(ni(2),nj(2)),v_left(ni(2),nj(2)), &
        T_left(ni(2),nj(2)),rho_left(ni(2),nj(2)))
    first=last+1
    last=last+ni(2)*nj(2)
    xOffsetLeft=xfirst(2)-0.5d0*dx(2)
    yOffsetLeft=yfirst(2)-0.5d0*dx(2)
    u_left=reshape(uall(first:last),[ni(2),nj(2)])
    v_left=reshape(vall(first:last),[ni(2),nj(2)])
    T_left=reshape(Tall(first:last),[ni(2),nj(2)])
    rho_left=reshape(rhoall(first:last),[ni(2),nj(2)])
    if(allocated(T_right)) deallocate(u_right,v_right,T_right,rho_right)
    allocate(u_right(ni(3),nj(3)),v_right(ni(3),nj(3)), &
        T_right(ni(3),nj(3)),rho_right(ni(3),nj(3)))
    first=last+1
    last=last+ni(3)*nj(3)
    xOffsetRight=xfirst(3)-0.5d0*dx(3)
    yOffsetRight=yfirst(3)-0.5d0*dx(3)
    u_right=reshape(uall(first:last),[ni(3),nj(3)])
    v_right=reshape(vall(first:last),[ni(3),nj(3)])
    T_right=reshape(Tall(first:last),[ni(3),nj(3)])
    rho_right=reshape(rhoall(first:last),[ni(3),nj(3)])
    if(allocated(T_bottom)) deallocate(u_bottom,v_bottom,T_bottom,rho_bottom)
    allocate(u_bottom(ni(4),nj(4)),v_bottom(ni(4),nj(4)), &
        T_bottom(ni(4),nj(4)),rho_bottom(ni(4),nj(4)))
    first=last+1
    last=last+ni(4)*nj(4)
    xOffsetBottom=xfirst(4)-0.5d0*dx(4)
    yOffsetBottom=yfirst(4)-0.5d0*dx(4)
    u_bottom=reshape(uall(first:last),[ni(4),nj(4)])
    v_bottom=reshape(vall(first:last),[ni(4),nj(4)])
    T_bottom=reshape(Tall(first:last),[ni(4),nj(4)])
    rho_bottom=reshape(rhoall(first:last),[ni(4),nj(4)])
    if(allocated(T_top)) deallocate(u_top,v_top,T_top,rho_top)
    allocate(u_top(ni(5),nj(5)),v_top(ni(5),nj(5)), &
        T_top(ni(5),nj(5)),rho_top(ni(5),nj(5)))
    first=last+1
    last=last+ni(5)*nj(5)
    xOffsetTop=xfirst(5)-0.5d0*dx(5)
    yOffsetTop=yfirst(5)-0.5d0*dx(5)
    u_top=reshape(uall(first:last),[ni(5),nj(5)])
    v_top=reshape(vall(first:last),[ni(5),nj(5)])
    T_top=reshape(Tall(first:last),[ni(5),nj(5)])
    rho_top=reshape(rhoall(first:last),[ni(5),nj(5)])
    endif
    deallocate(uall,vall,Tall,rhoall)
end subroutine read_snapshot

! 以下空间积分与主程序采用相同公式；面积从快照几何重建，不从文件读取面积数组。
  subroutine calculate_NuRe(nu,re,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax)

    use postdata
    implicit none

    ! 依次累加五个物理分区，节点实际贡献的面积在各区统计循环内直接计算。
    ! Re 使用速度平方的面积平均再开方；不能直接平均各网格各自计算的 Re。

    real(kind=8), intent(out) :: nu, re, hot, cold, middle, mass, meanT, tmin, tmax, rmin, rmax
    real(kind=8) :: conv, vel2


    conv=0.0d0
    vel2=0.0d0
    mass=0.0d0
    meanT=0.0d0
    hot=0.0d0
    cold=0.0d0
    middle=0.0d0
    tmin=huge(1.0d0)
    tmax=-tmin
    rmin=tmin
    rmax=-tmin

    call calNuRe_grid(size(T_coarse,1),size(T_coarse,2),dxCoarse,xOffsetCoarse,yOffsetCoarse,(1+merge(interfaceSkin,0,refineRatio>1)),&
        (size(T_coarse,1)-merge(interfaceSkin,0,refineRatio>1)),(1+merge(interfaceSkin,0,refineRatio>1)),(size(T_coarse,2)-merge(interfaceSkin,0,refineRatio>1)),centerIntegrationBox,&
        conv,vel2,mass,meanT,hot,cold,&
        middle,tmin,tmax,rmin,rmax,rho_coarse,u_coarse,v_coarse,T_coarse)
    if(refineRatio > 1) then

        call calNuRe_grid(size(T_left,1),size(T_left,2),dxFine,xOffsetLeft,yOffsetLeft,1,size(T_left,1),&
            1,size(T_left,2),[xOffsetLeft,xOffsetLeft+size(T_left,1),yOffsetLeft,yOffsetLeft+size(T_left,2)],conv,vel2,mass,meanT,hot,cold,&
            middle,tmin,tmax,&
            rmin,rmax,rho_left,u_left,v_left,T_left)

        call calNuRe_grid(size(T_right,1),size(T_right,2),dxFine,xOffsetRight,yOffsetRight,1,size(T_right,1),&
            1,size(T_right,2),[xOffsetRight,xOffsetRight+size(T_right,1),yOffsetRight,yOffsetRight+size(T_right,2)],conv,vel2,mass,&
            meanT,hot,cold,&
            middle,tmin,tmax,&
            rmin,rmax,rho_right,u_right,v_right,T_right)

        call calNuRe_grid(size(T_bottom,1),size(T_bottom,2),dxFine,xOffsetBottom,yOffsetBottom,1,&
            size(T_bottom,1),1,size(T_bottom,2),[xOffsetBottom,xOffsetBottom+size(T_bottom,1),yOffsetBottom,yOffsetBottom+size(T_bottom,2)],conv,&
            vel2,mass,meanT,hot,cold,&
            middle,tmin,tmax,rmin,rmax,rho_bottom,u_bottom,v_bottom,T_bottom)

        call calNuRe_grid(size(T_top,1),size(T_top,2),dxFine,xOffsetTop,yOffsetTop,1,size(T_top,1),&
            1,&
            size(T_top,2),[xOffsetTop,xOffsetTop+size(T_top,1),yOffsetTop,yOffsetTop+size(T_top,2)],conv,vel2,mass,meanT,hot,cold,&
            middle,tmin,tmax,rmin,rmax,&
            rho_top,u_top,v_top,T_top)
    endif

    ! 空间积分等于中心、左、右、下、上五个物理分区之和；接口端点按实际面积计权。
    ! Re 使用速度平方的面积平均再开方，非稳态时间统计也采用同一 RMS 定义。
    nu=1.0d0+conv/(dble(nx)*dble(ny))*(lengthUnit/(Thot-Tcold))/diffusivity
    re=sqrt(vel2/(dble(nx)*dble(ny)))*lengthUnit/viscosity
    hot=(lengthUnit/(Thot-Tcold))*hot
    cold=(lengthUnit/(Thot-Tcold))*cold
    middle=(lengthUnit/(Thot-Tcold))*middle
    meanT=meanT/(dble(nx)*dble(ny))
  end subroutine calculate_NuRe

  subroutine calNuRe_grid(ni,nj,dx,xOffset,yOffset,xStartIntegration,xEndIntegration,yStartIntegration,&
      yEndIntegration,integrationBox,conv,&
      vel2,mass,meanT,hot,cold,middle,tmin,tmax,rmin,rmax,&
      rho,u,v,T)

    use postdata,only: nx,ny,Thot,Tcold,diffusivity,itc,reloadItcOffset,refineRatio,ieee_is_finite,centerIntegrationBox,sideHeated
    implicit none

    ! integrationBox 为细数组外框或粗区物理矩形；零面积交换节点不参与统计。
    ! 细中线模板跨存储接缝时由 fine_section_x/y 按全局坐标取值，不能截成单侧模板。

    integer(kind=4) :: ni, nj, xStartIntegration, xEndIntegration, yStartIntegration, yEndIntegration
    real(kind=8) :: dx, xOffset, yOffset, integrationBox(4)

    real(kind=8) :: x, y, xlo, xhi, ylo, yhi, middleWidth
    integer(kind=4) :: i, j, jm, im
    real(kind=8) :: rho(ni,nj)
    real(kind=8) :: u(ni,nj)
    real(kind=8) :: v(ni,nj)
    real(kind=8) :: T(ni,nj)
    real(kind=8) :: conv, vel2, mass, meanT, hot, cold, middle, dTdx, dTdy, tm, um, vm, cellArea, sectionWeight
    real(kind=8) :: w(4), dw(4)
    real(kind=8) :: tmin, tmax, rmin, rmax

    do j=yStartIntegration,yEndIntegration
        do i=xStartIntegration,xEndIntegration
            ! 节点控制面积先裁剪到本区范围；这里的坐标仍以细格距为单位。
            x=xOffset+(i-0.5d0)*dx
            y=yOffset+(j-0.5d0)*dx
            xlo=max(x-dx/2,integrationBox(1))
            xhi=min(x+dx/2,integrationBox(2))
            ylo=max(y-dx/2,integrationBox(3))
            yhi=min(y+dx/2,integrationBox(4))
            cellArea=0.0d0
            if(xhi > xlo .AND. yhi > ylo) then
                if(refineRatio==1 .OR. dx > 1.0d0) then
                    cellArea=(xhi-xlo)*(yhi-ylo)    ! 单网格全域或中心粗物理矩形。
                else
                    ! 分别累加四个细物理矩形的面积，排除中心粗物理区域。
                    cellArea=max(0.0d0,min(xhi,centerIntegrationBox(1))-xlo)*(yhi-ylo)    ! 左
                    cellArea=cellArea+max(0.0d0,xhi-max(xlo,centerIntegrationBox(2)))*(yhi-ylo)    ! 右
                    middleWidth=max(0.0d0,min(xhi,centerIntegrationBox(2))-max(xlo,centerIntegrationBox(1)))
                    cellArea=cellArea+middleWidth*max(0.0d0,min(yhi,centerIntegrationBox(3))-ylo)    ! 下
                    cellArea=cellArea+middleWidth*max(0.0d0,yhi-max(ylo,centerIntegrationBox(4)))    ! 上
                endif
            endif
            if(cellArea <=0d0) cycle
            if(.NOT. ieee_is_finite(T(i,j)) .OR. .NOT. ieee_is_finite(rho(i,j)) .OR. &
                .NOT. ieee_is_finite(u(i,j)) .OR. .NOT. ieee_is_finite(v(i,j)) .OR. rho(i,j) <=0.0d0) then
                write(*,*) 'Invalid state: time, xOffset,yOffset, i,j:',reloadItcOffset+itc,xOffset,yOffset,i,j
                error stop 'Nonfinite or nonpositive density in owned cells'
            endif
    if(sideHeated) then
            conv=conv+u(i,j)*T(i,j)*cellArea
    else
            conv=conv+v(i,j)*T(i,j)*cellArea
    endif
            vel2=vel2+(u(i,j)**2+v(i,j)**2)*cellArea
            mass=mass+rho(i,j)*cellArea
            meanT=meanT+T(i,j)*cellArea
            tmin=min(tmin,T(i,j))
            tmax=max(tmax,T(i,j))
            rmin=min(rmin,rho(i,j))
            rmax=max(rmax,rho(i,j))
        enddo
    enddo

    if(sideHeated) then
    if(integrationBox(1)==0.0d0) then
        do j=yStartIntegration,yEndIntegration
            hot=hot+(8.0d0*Thot-9.0d0*T(1,j)+T(2,j))/(3.0d0*dx)*dx/dble(ny)
        enddo
    endif

    if(integrationBox(2)==dble(nx)) then
        do j=yStartIntegration,yEndIntegration
            cold=cold+(-8.0d0*Tcold+9.0d0*T(ni,j)-T(ni-1,&
                j))/(3.0d0*dx)*dx/dble(ny)
        enddo
    endif

    if((0.5d0*nx) >=integrationBox(1) .AND. (0.5d0*nx) < integrationBox(2)) then

        call section_weights(((0.5d0*nx)-xOffset)/dx+0.5d0,ni,im,w,dw)
        do j=yStartIntegration,yEndIntegration
            call section_owned_weight(xOffset,yOffset,dx,integrationBox,j,1,(0.5d0*nx),sectionWeight)
            if(sectionWeight <=0d0) cycle
            if(dx==1.0d0 .AND. refineRatio > 1) then

                call fine_section_x((0.5d0*nx),nint(yOffset)+j,tm,um,dTdx)
            else
                tm=sum(w*T(im:im+3,j))
                um=sum(w*u(im:im+3,j))
                dTdx=sum(dw*T(im:im+3,j))/dx
            endif

            middle=middle+(um*tm/diffusivity-dTdx)*sectionWeight/dble(ny)
        enddo
    endif
    else

    if(integrationBox(3)==0.0d0) then
        do i=xStartIntegration,xEndIntegration
            hot=hot+(8.0d0*Thot-9.0d0*T(i,1)+T(i,2))/(3.0d0*dx)*dx/dble(nx)
        enddo
    endif

    if(integrationBox(4)==dble(ny)) then
        do i=xStartIntegration,xEndIntegration
            cold=cold+(-8.0d0*Tcold+9.0d0*T(i,nj)-T(i,&
                nj-1))/(3.0d0*dx)*dx/dble(nx)
        enddo
    endif

    if((0.5d0*ny) >=integrationBox(3) .AND. (0.5d0*ny) < integrationBox(4)) then

        call section_weights(((0.5d0*ny)-yOffset)/dx+0.5d0,nj,jm,w,dw)
        do i=xStartIntegration,xEndIntegration
            call section_owned_weight(xOffset,yOffset,dx,integrationBox,i,2,(0.5d0*ny),sectionWeight)
            if(sectionWeight <=0d0) cycle
            if(dx==1.0d0 .AND. refineRatio > 1) then

                call fine_section_y(nint(xOffset)+i,(0.5d0*ny),tm,vm,dTdy)
            else
                tm=sum(w*T(i,jm:jm+3))
                vm=sum(w*v(i,jm:jm+3))
                dTdy=sum(dw*T(i,jm:jm+3))/dx
            endif

            middle=middle+(vm*tm/diffusivity-dTdy)*sectionWeight/dble(nx)
        enddo
    endif
    endif
  end subroutine calNuRe_grid

  subroutine section_weights(q,n,first,w,dw)

    implicit none

    ! q 为连续节点下标，n 为该方向节点数；first 是四点模板的首下标。
    ! w 用于函数值，dw 是对节点下标的导数权重，换成物理导数时还要除以 dx。

    real(kind=8), intent(in) :: q
    integer(kind=4), intent(in) :: n
    integer(kind=4), intent(out) :: first
    real(kind=8), intent(out) :: w(4), dw(4)
    real(kind=8) :: z, term
    integer(kind=4) :: a, k, l

    first=max(1,min(floor(q)-1,n-3))
    z=q-dble(first)

    call lagrange_weights(z,w)
    dw=0.0d0
    do a=0,3
        do k=0,3
            if(k==a) cycle
            term=1.0d0/dble(a-k)
            do l=0,3
                if(l/=a .AND. l/=k) term=term*(z-dble(l))/dble(a-l)
            enddo

            dw(a+1)=dw(a+1)+term
        enddo
    enddo
  end subroutine section_weights

  subroutine lagrange_weights(q,w)

    implicit none
    !$acc routine seq

    ! q 为目标点相对于模板首节点的距离，以来源网格格距为单位。
    ! 模板的相对坐标是 0、1、2、3；w(1:4) 是三次多项式的四点插值权重。

    real(kind=8), intent(in) :: q
    real(kind=8), intent(out) :: w(4)
    integer(kind=4) :: a, k

    w=1.0d0   !下面通过连续相乘来计算权重，所以初值必须是 1,如果初值为 0，结果都会是 0
    do a=0,3
        do k=0,3
            if(a/=k) w(a+1)=w(a+1)*(q-dble(k))/dble(a-k)
        enddo
    enddo
  end subroutine lagrange_weights

  subroutine section_owned_weight(xOffset,yOffset,dx,integrationBox,k,axis,position,w)
    use postdata,only: centerIntegrationBox,refineRatio
    implicit none

    real(kind=8), intent(out) :: w    ! 输出：本节点在指定中线上实际负责的积分长度。

    ! axis=1 表示固定 x 的竖直中线，axis=2 表示固定 y 的水平中线。
    ! 先定位本节点控制区间；细区处于中心投影内时只累加两侧物理线段。

    real(kind=8), intent(in) :: xOffset, yOffset, dx, integrationBox(4), position
    integer(kind=4), intent(in) :: k, axis
    real(kind=8) :: q, lo, hi, clo, chi

    if(axis==1) then
        q=yOffset+(k-.5d0)*dx
        lo=integrationBox(3)
        hi=integrationBox(4)
        clo=centerIntegrationBox(3)
        chi=centerIntegrationBox(4)
    else
        q=xOffset+(k-.5d0)*dx
        lo=integrationBox(1)
        hi=integrationBox(2)
        clo=centerIntegrationBox(1)
        chi=centerIntegrationBox(2)
    endif

    w=max(0d0,min(q+dx/2,hi)-max(q-dx/2,lo))
    if(refineRatio > 1 .AND. dx==1d0) then
        if(axis==1) then
            if(position < centerIntegrationBox(1) .OR. position >=centerIntegrationBox(2)) return
        else
            if(position < centerIntegrationBox(3) .OR. position >=centerIntegrationBox(4)) return
        endif

        ! 中心投影内的细中线由两段物理区间相加，不计算交换层长度。
        w=max(0d0,min(q+dx/2,hi,clo)-max(q-dx/2,lo)) &
            +max(0d0,min(q+dx/2,hi)-max(q-dx/2,lo,chi))
    endif
  end subroutine section_owned_weight

  subroutine fine_section_x(x,j,tm,um,dTdx)
    use postdata
    implicit none

    ! 模板位置按全域 nx 定位；四个温度和速度样本分别从相应具名细数组读取。

    real(kind=8), intent(in) :: x
    integer(kind=4), intent(in) :: j
    real(kind=8), intent(out) :: tm, um, dTdx
    real(kind=8) :: w(4), dw(4), temp(4), velocity(4)
    integer(kind=4) :: first, a

    call section_weights(x+0.5d0,nx,first,w,dw)
    do a=1,4
        call fine_value(first+a-1,j,T_left,T_right,T_bottom,T_top,temp(a))
        call fine_value(first+a-1,j,u_left,u_right,u_bottom,u_top,velocity(a))
    enddo

    tm=sum(w*temp)
    um=sum(w*velocity)
    dTdx=sum(dw*temp)
  end subroutine fine_section_x

  subroutine fine_section_y(i,y,tm,vm,dTdy)
    use postdata
    implicit none

    ! 模板位置按全域 ny 定位；四个温度和速度样本分别从相应具名细数组读取。

    integer(kind=4), intent(in) :: i
    real(kind=8), intent(in) :: y
    real(kind=8), intent(out) :: tm, vm, dTdy
    real(kind=8) :: w(4), dw(4), temp(4), velocity(4)
    integer(kind=4) :: first, a

    call section_weights(y+0.5d0,ny,first,w,dw)
    do a=1,4
        call fine_value(i,first+a-1,T_left,T_right,T_bottom,T_top,temp(a))
        call fine_value(i,first+a-1,v_left,v_right,v_bottom,v_top,velocity(a))
    enddo

    tm=sum(w*temp)
    vm=sum(w*velocity)
    dTdy=sum(dw*temp)
  end subroutine fine_section_y

  subroutine fine_value(i,j,left,right,bottom,top,value)
    use postdata,only: T_left,T_right,T_bottom,T_top,&
        xOffsetRight,xOffsetBottom,xOffsetTop,yOffsetTop
    implicit none

    real(kind=8), intent(out) :: value    ! 输出：指定全局细节点对应的温度或速度等标量值。

    ! i/j 是全局细节点编号；减去所选具名数组的整数偏移后才是本地下标。
    ! 只用于主机端中线后处理；若请求落入未分配的中心空区，立即报错。

    integer(kind=4), intent(in) :: i, j
    real(kind=8), intent(in) :: left(size(T_left,1),size(T_left,2)), right(size(T_right,1),size(T_right,2)), &
        bottom(size(T_bottom,1),size(T_bottom,2)),top(size(T_top,1),size(T_top,2))

    if(i <=size(T_left,1)) then
        value=left(i,j)
    else if(i > nint(xOffsetRight)) then
        value=right(i-nint(xOffsetRight),j)
    else if(j <=size(T_bottom,2)) then
        value=bottom(i-nint(xOffsetBottom),j)
    else if(j > nint(yOffsetTop)) then
        value=top(i-nint(xOffsetTop),j-nint(yOffsetTop))
    else
        error stop 'Fine section requested an inactive central node'
    endif
  end subroutine fine_value

