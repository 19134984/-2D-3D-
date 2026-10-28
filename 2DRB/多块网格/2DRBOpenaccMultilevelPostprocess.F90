! 独立后处理程序：CPU运行，不依赖OpenACC或求解器的commondata模块。
! 修改下面宏开关和参数区后重新编译本程序即可重新统计，不必重新运行求解器。
!=============================================================
! 自定义宏：与主程序一样，稳态和非稳态只能选择一个。
!#define steadyFlow
#define unsteadyFlow

#if defined(steadyFlow) && defined(unsteadyFlow)
#error "Choose only one flow mode: steadyFlow or unsteadyFlow"
#endif
#if !defined(steadyFlow) && !defined(unsteadyFlow)
#error "Define one flow mode: steadyFlow or unsteadyFlow"
#endif

!=============================================================

! 原后处理程序的独立副本：历史文件列顺序、统计窗口和输出文件名保持一致。
! processSnapshots=.true.时，可同时读取旧MB2DSNAPSHOT0002和新ML2DSNAPSHOT0001快照。
! 稳态最终Nu/Re仍由主求解器直接输出，独立后处理负责整理曲线及可选快照核对。
module postdata
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    logical, parameter :: processSnapshots=.false.    ! 是否读取快照重算瞬时空间指标：.true.=开启，.false.=关闭。
    ! 非稳态时间窗口
    real(kind=8), parameter :: unsteadyAverageStartTf=500.0d0 !非稳态统计起点
    real(kind=8), parameter :: unsteadyAverageEndTf=1000.0d0  !非稳态统计终点
    real(kind=8), parameter :: unsteadyAverageMidTf= &        !非稳态统计中点
        0.5d0*(unsteadyAverageStartTf+unsteadyAverageEndTf)
    character(*), parameter :: NuReHistoryFile='NuRe_2DOpenaccMultiblock.dat'
    character(*), parameter :: outputIndexFile='OutputFiles_2DOpenaccMultiblock.dat' !已输出文件的累计步数和文件名
    character(*), parameter :: snapshotFilePrefix='buoyancyCavity2DOpenaccMultiblockSnapshot'

    integer :: nx,ny,nGridLevels,nBlocks,interfaceSkin,totalItc,itc=0,reloadItcOffset=0
    real(kind=8) :: totalTff,lengthUnit,timeUnit,viscosity,diffusivity,Thot,Tcold,Tref
    logical :: sideHeated
    type grid_block
        integer :: level,ni,nj,globalI,globalJ
        real(kind=8) :: dx,xOffset,yOffset
        logical :: wall(4)
        real(kind=8),allocatable :: rho(:,:),u(:,:),v(:,:),T(:,:),area(:,:)
    end type grid_block
    type grid_level
        integer :: firstBlock,lastBlock,ni,nj
        real(kind=8) :: dx,xOffset,yOffset,outerBox(4),innerBox(4)
    end type grid_level
    type(grid_block),allocatable :: block(:)
    type(grid_level),allocatable :: grid(:)
contains
!===================================================================================================
! 根据文件标识选择读取器；不把多级结果重采样成均匀网格，直接按原分区积分。
!===================================================================================================
    subroutine read_snapshot(name)
        character(len=*),intent(in) :: name
        character(len=16) :: magic
        integer :: k,ios
        open(newunit=k,file=name,access='stream',form='unformatted',status='old',iostat=ios)
        if(ios/=0) error stop 'Cannot open snapshot'
        read(k,iostat=ios) magic
        if(ios/=0) error stop 'Incomplete snapshot identifier'
        rewind(k)
        if(allocated(block)) deallocate(block,grid)
        select case(magic)
        case('ML2DSNAPSHOT0001')
            call read_multilevel_snapshot(k)
        case('MB2DSNAPSHOT0002')
            call read_legacy_snapshot(k)
        case default
            error stop 'Unsupported snapshot format'
        end select
        close(k)
        reloadItcOffset=totalItc
        call prepare_geometry()
    end subroutine read_snapshot

    subroutine allocate_fields(b)
        type(grid_block),intent(inout) :: b
        allocate(b%u(b%ni,b%nj),b%v(b%ni,b%nj),b%T(b%ni,b%nj),b%rho(b%ni,b%nj),b%area(b%ni,b%nj))
    end subroutine allocate_fields

    subroutine read_multilevel_snapshot(k)
        integer,intent(in) :: k
        integer :: ios,config(16),b,l
        integer(kind=8) :: position,fileSize
        real(kind=8) :: physics(12),outer(4),inner(4)
        logical,allocatable :: seen(:)
        character(len=16) :: magic
        read(k,iostat=ios) magic,config,physics,totalItc,nBlocks,timeUnit
        if(ios/=0) error stop 'Incomplete multilevel snapshot header'
        nx=config(1)
        ny=config(2)
        nGridLevels=config(3)
        interfaceSkin=config(4)
        if(nGridLevels<1.or.min(nx,ny)<4.or.interfaceSkin<2) error stop 'Invalid multilevel snapshot dimensions'
        if(nBlocks/=4*(nGridLevels-1)+1.or.config(9)/=2**(nGridLevels-1)) &
            error stop 'Invalid multilevel snapshot hierarchy'
        if(totalItc<0.or.mod(totalItc,config(9))/=0) error stop 'Unsynchronized multilevel snapshot'
        if(config(10)/=1.and.config(10)/=2) error stop 'Invalid snapshot convection case'
        if(.not.all(ieee_is_finite(physics)).or..not.ieee_is_finite(timeUnit)) error stop 'Nonfinite snapshot physics'
        if(physics(2)<=0.0d0.or.physics(7)<=0.0d0) error stop 'Invalid snapshot transport parameters'
        sideHeated=config(10)==1
        lengthUnit=dble(ny)
        if(sideHeated) lengthUnit=dble(nx)
        Thot=physics(4)
        Tcold=physics(5)
        Tref=0.5d0*(Thot+Tcold)
        viscosity=(1.0d0/physics(7)-0.5d0)/3.0d0
        diffusivity=viscosity/physics(2)
        if(min(viscosity,diffusivity,timeUnit)<=0.0d0.or.Thot==Tcold) error stop 'Invalid snapshot units'
        totalTff=dble(totalItc)/timeUnit
        allocate(block(nBlocks),grid(nGridLevels),seen(nGridLevels))
        seen=.false.
        do b=1,nBlocks
            read(k,iostat=ios) l,block(b)%ni,block(b)%nj,block(b)%dx,block(b)%xOffset,block(b)%yOffset,outer,inner
            if(ios/=0) error stop 'Incomplete multilevel snapshot geometry'
            if(l<1.or.l>nGridLevels.or.min(block(b)%ni,block(b)%nj)<1) error stop 'Invalid snapshot block dimensions'
            block(b)%level=l
            if(.not.all(ieee_is_finite([block(b)%dx,block(b)%xOffset,block(b)%yOffset,outer,inner]))) &
                error stop 'Nonfinite snapshot geometry'
            if(block(b)%dx/=dble(2**(l-1))) error stop 'Snapshot level spacing mismatch'
            if(seen(l)) then
                if(any(abs(outer-grid(l)%outerBox)>1.0d-12).or.any(abs(inner-grid(l)%innerBox)>1.0d-12)) &
                    error stop 'Inconsistent geometry within one level'
            else
                grid(l)%outerBox=outer
                grid(l)%innerBox=inner
                grid(l)%dx=block(b)%dx
                seen(l)=.true.
            endif
            call allocate_fields(block(b))
            read(k,iostat=ios) block(b)%u,block(b)%v,block(b)%T,block(b)%rho
            if(ios/=0) error stop 'Incomplete multilevel snapshot fields'
        enddo
        if(.not.all(seen)) error stop 'A snapshot grid level is missing'
        inquire(unit=k,pos=position,size=fileSize)
        if(position-1/=fileSize) error stop 'Unexpected trailing snapshot bytes'
    end subroutine read_multilevel_snapshot

!===================================================================================================
! 原五块快照读取顺序保持原样：先读全部几何，再读拼接的四种场，最后还原各块。
!===================================================================================================
    subroutine read_legacy_snapshot(k)
        integer,intent(in) :: k
        integer :: ios,caseMode,ratio,fineLeft,fineRight,fineBottom,fineTop,b,l,count,first,last
        integer,allocatable :: ni(:),nj(:)
        real(kind=8),allocatable :: xoff(:),yoff(:),dx(:),uall(:),vall(:),Tall(:),rhoall(:)
        real(kind=8) :: center(4)
        character(len=16) :: magic
        read(k,iostat=ios) magic,nBlocks,nx,ny,totalItc,caseMode
        if(ios/=0) error stop 'Incomplete legacy snapshot header'
        if((nBlocks/=1.and.nBlocks/=5).or.min(nx,ny)<4) error stop 'Invalid legacy snapshot mesh'
        if(caseMode/=0.and.caseMode/=1) error stop 'Invalid legacy snapshot case'
        read(k,iostat=ios) ratio,interfaceSkin,fineLeft,fineRight,fineBottom,fineTop
        if(ios/=0.or.ratio<1.or.interfaceSkin<0) error stop 'Invalid legacy snapshot refinement'
        if((nBlocks==1.and.ratio/=1).or.(nBlocks==5.and.ratio<=1)) error stop 'Invalid legacy hierarchy'
        if(totalItc<0.or.mod(totalItc,ratio)/=0) error stop 'Invalid legacy snapshot step'
        read(k,iostat=ios) totalTff,lengthUnit,viscosity,diffusivity,Thot,Tcold,Tref,timeUnit,center
        if(ios/=0) error stop 'Incomplete legacy snapshot physics'
        if(.not.all(ieee_is_finite([totalTff,lengthUnit,viscosity,diffusivity,Thot,Tcold,Tref,timeUnit,center]))) &
            error stop 'Nonfinite legacy snapshot physics'
        if(min(lengthUnit,viscosity,diffusivity,timeUnit)<=0.0d0.or.Thot==Tcold) error stop 'Invalid legacy snapshot units'
        if(abs(totalTff-dble(totalItc)/timeUnit)>1.0d-12*max(1.0d0,abs(totalTff))) error stop 'Legacy snapshot time mismatch'
        sideHeated=caseMode==1
        nGridLevels=1
        if(nBlocks==5) nGridLevels=2
        allocate(block(nBlocks),grid(nGridLevels),ni(nBlocks),nj(nBlocks),xoff(nBlocks),yoff(nBlocks),dx(nBlocks))
        count=0
        do b=1,nBlocks
            read(k,iostat=ios) ni(b),nj(b),xoff(b),yoff(b),dx(b)
            if(ios/=0.or.min(ni(b),nj(b))<1.or.dx(b)<=0.0d0) error stop 'Invalid legacy block'
            if(.not.all(ieee_is_finite([dx(b),xoff(b),yoff(b)]))) error stop 'Nonfinite legacy geometry'
            count=count+ni(b)*nj(b)
        enddo
        if(dx(1)/=dble(ratio)) error stop 'Legacy coarse spacing mismatch'
        if(nBlocks==5) then
            if(any(dx(2:)/=1.0d0)) error stop 'Legacy fine spacing mismatch'
        endif
        allocate(uall(count),vall(count),Tall(count),rhoall(count))
        read(k,iostat=ios) uall,vall,Tall,rhoall
        if(ios/=0) error stop 'Incomplete legacy snapshot fields'
        grid(1)%dx=1.0d0
        grid(1)%outerBox=[0.0d0,dble(nx),0.0d0,dble(ny)]
        grid(1)%innerBox=0.0d0
        if(nBlocks==5) then
            grid(1)%innerBox=center
            grid(2)%dx=dble(ratio)
            grid(2)%outerBox=center
            grid(2)%innerBox=0.0d0
        endif
        last=0
        do b=1,nBlocks
            l=1
            if(nBlocks==5.and.b==1) l=2
            block(b)%level=l
            block(b)%ni=ni(b)
            block(b)%nj=nj(b)
            block(b)%dx=dx(b)
            block(b)%xOffset=xoff(b)
            block(b)%yOffset=yoff(b)
            call allocate_fields(block(b))
            first=last+1
            last=last+ni(b)*nj(b)
            block(b)%u=reshape(uall(first:last),[ni(b),nj(b)])
            block(b)%v=reshape(vall(first:last),[ni(b),nj(b)])
            block(b)%T=reshape(Tall(first:last),[ni(b),nj(b)])
            block(b)%rho=reshape(rhoall(first:last),[ni(b),nj(b)])
        enddo
    end subroutine read_legacy_snapshot

    subroutine prepare_geometry()
        integer :: l,b,i,j,localCount
        real(kind=8) :: d,x,y,totalArea
        do l=1,nGridLevels
            d=grid(l)%dx
            if(l==1) then
                grid(l)%xOffset=0.0d0
                grid(l)%yOffset=0.0d0
                grid(l)%ni=nx
                grid(l)%nj=ny
            else
                grid(l)%xOffset=grid(l)%outerBox(1)-(interfaceSkin+0.5d0)*d
                grid(l)%yOffset=grid(l)%outerBox(3)-(interfaceSkin+0.5d0)*d
                grid(l)%ni=nint((grid(l)%outerBox(2)-grid(l)%outerBox(1))/d)+1+2*interfaceSkin
                grid(l)%nj=nint((grid(l)%outerBox(4)-grid(l)%outerBox(3))/d)+1+2*interfaceSkin
            endif
            grid(l)%firstBlock=nBlocks+1
            grid(l)%lastBlock=0
            localCount=0
            do b=1,nBlocks
                if(block(b)%level/=l) cycle
                localCount=localCount+1
                grid(l)%firstBlock=min(grid(l)%firstBlock,b)
                grid(l)%lastBlock=max(grid(l)%lastBlock,b)
            enddo
            if(localCount/=merge(1,4,l==nGridLevels)) error stop 'Invalid number of blocks in a grid level'
        enddo
        totalArea=0.0d0
        do b=1,nBlocks
            l=block(b)%level
            d=block(b)%dx
            block(b)%globalI=nint((block(b)%xOffset-grid(l)%xOffset)/d)+1
            block(b)%globalJ=nint((block(b)%yOffset-grid(l)%yOffset)/d)+1
            block(b)%wall=.false.
            if(l==1) block(b)%wall=[abs(block(b)%xOffset)<1.0d-12, &
                abs(block(b)%xOffset+block(b)%ni*d-nx)<1.0d-12,abs(block(b)%yOffset)<1.0d-12, &
                abs(block(b)%yOffset+block(b)%nj*d-ny)<1.0d-12]
            do j=1,block(b)%nj
                y=block(b)%yOffset+(j-0.5d0)*d
                do i=1,block(b)%ni
                    x=block(b)%xOffset+(i-0.5d0)*d
                    block(b)%area(i,j)=owned_cell_area(l,x,y)
                enddo
            enddo
            totalArea=totalArea+sum(block(b)%area)
        enddo
        if(abs(totalArea-dble(nx)*ny)>1.0d-10*dble(nx)*ny) error stop 'Snapshot geometry does not cover the domain once'
    end subroutine prepare_geometry
    real(kind=8) function rectangle_area(rect,box) result(area)
        real(kind=8),intent(in) :: rect(4),box(4)
        area=max(0.0d0,min(rect(2),box(2))-max(rect(1),box(1)))* &
            max(0.0d0,min(rect(4),box(4))-max(rect(3),box(3)))
    end function rectangle_area
    real(kind=8) function owned_cell_area(l,x,y) result(area)
        integer,intent(in) :: l
        real(kind=8),intent(in) :: x,y
        real(kind=8) :: d,rect(4)
        d=0.5d0*grid(l)%dx
        rect=[x-d,x+d,y-d,y+d]
        area=rectangle_area(rect,grid(l)%outerBox)
        if(l<nGridLevels) area=area-rectangle_area(rect,grid(l)%innerBox)
        area=max(0.0d0,area)
    end function owned_cell_area

!===================================================================================================
subroutine locate_level_node(l,gi,gj,b,i,j,owned)
        integer,intent(in) :: l,gi,gj
        integer,intent(out) :: b,i,j
        logical,intent(in) :: owned
        integer :: k
        b=0
        i=0
        j=0
        do k=grid(l)%firstBlock,grid(l)%lastBlock
            if(gi<block(k)%globalI.or.gi>=block(k)%globalI+block(k)%ni) cycle
            if(gj<block(k)%globalJ.or.gj>=block(k)%globalJ+block(k)%nj) cycle
            b=k
            i=gi-block(k)%globalI+1
            j=gj-block(k)%globalJ+1
            if(owned) then
                if(block(k)%area(i,j)<=0.0d0) error stop 'Interface stencil enters a receiver: increase ring thickness'
            endif
            return
        enddo
    end subroutine locate_level_node

!===================================================================================================
subroutine calculate_NuRe(nu,re,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax)
        real(kind=8),intent(out) :: nu,re,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax
        real(kind=8) :: conv,vel2,area,d,weight
        integer :: b,i,j
        conv=0.0d0
        vel2=0.0d0
        mass=0.0d0
        meanT=0.0d0
        hot=0.0d0
        cold=0.0d0
        tmin=huge(1.0d0)
        tmax=-tmin
        rmin=tmin
        rmax=-tmin
        do b=1,nBlocks
            d=block(b)%dx
            do j=1,block(b)%nj
                do i=1,block(b)%ni
                    area=block(b)%area(i,j)
                    if(area<=0.0d0) cycle
                    if(.not.ieee_is_finite(block(b)%rho(i,j)).or.block(b)%rho(i,j)<=0.0d0.or. &
                        .not.ieee_is_finite(block(b)%T(i,j)).or..not.ieee_is_finite(block(b)%u(i,j)).or. &
                        .not.ieee_is_finite(block(b)%v(i,j))) then
                        write(*,*) 'Invalid owned cell: step,block,level,i,j=',reloadItcOffset+itc,b,block(b)%level,i,j
                        error stop 'Nonfinite field or nonpositive density'
                    endif
        if(sideHeated) then
                    conv=conv+area*block(b)%u(i,j)*block(b)%T(i,j)
        else
                    conv=conv+area*block(b)%v(i,j)*block(b)%T(i,j)
        endif
                    vel2=vel2+area*(block(b)%u(i,j)**2+block(b)%v(i,j)**2)
                    mass=mass+area*block(b)%rho(i,j)
                    meanT=meanT+area*block(b)%T(i,j)
                    tmin=min(tmin,block(b)%T(i,j))
                    tmax=max(tmax,block(b)%T(i,j))
                    rmin=min(rmin,block(b)%rho(i,j))
                    rmax=max(rmax,block(b)%rho(i,j))
                enddo
            enddo
        if(sideHeated) then
            if(block(b)%wall(1).or.block(b)%wall(2)) then
                do j=1,block(b)%nj
                    weight=d/dble(ny)
                    if(block(b)%wall(1)) hot=hot+weight*(8.0d0*Thot-9.0d0*block(b)%T(1,j)+block(b)%T(2,j))/(3.0d0*d)
                    i=block(b)%ni
                    if(block(b)%wall(2)) cold=cold+weight*(-8.0d0*Tcold+9.0d0*block(b)%T(i,j)-block(b)%T(i-1,j))/(3.0d0*d)
                enddo
            endif
        else
            if(block(b)%wall(3).or.block(b)%wall(4)) then
                do i=1,block(b)%ni
                    weight=d/dble(nx)
                    if(block(b)%wall(3)) hot=hot+weight*(8.0d0*Thot-9.0d0*block(b)%T(i,1)+block(b)%T(i,2))/(3.0d0*d)
                    j=block(b)%nj
                    if(block(b)%wall(4)) cold=cold+weight*(-8.0d0*Tcold+9.0d0*block(b)%T(i,j)-block(b)%T(i,j-1))/(3.0d0*d)
                enddo
            endif
        endif
        enddo
        nu=1.0d0+conv/(dble(nx)*dble(ny))*lengthUnit/((Thot-Tcold)*diffusivity)
        re=sqrt(vel2/(dble(nx)*dble(ny)))*lengthUnit/viscosity
        meanT=meanT/(dble(nx)*dble(ny))
        hot=hot*lengthUnit/(Thot-Tcold)
        cold=cold*lengthUnit/(Thot-Tcold)
        call calculate_middle_flux(middle)
    end subroutine calculate_NuRe

!===================================================================================================
    real(kind=8) function section_owned_length(l,q,position,axis) result(weight)
        integer,intent(in) :: l,axis
        real(kind=8),intent(in) :: q,position
        real(kind=8) :: lo,hi,d
        d=0.5d0*grid(l)%dx
        weight=0.0d0
        if(axis==1) then
            if(position<grid(l)%outerBox(1).or.position>=grid(l)%outerBox(2)) return
            lo=max(q-d,grid(l)%outerBox(3))
            hi=min(q+d,grid(l)%outerBox(4))
            weight=max(0.0d0,hi-lo)
            if(l<nGridLevels) then
                if(position>=grid(l)%innerBox(1).and.position<grid(l)%innerBox(2)) &
                    weight=weight-max(0.0d0,min(hi,grid(l)%innerBox(4))-max(lo,grid(l)%innerBox(3)))
            endif
        else
            if(position<grid(l)%outerBox(3).or.position>=grid(l)%outerBox(4)) return
            lo=max(q-d,grid(l)%outerBox(1))
            hi=min(q+d,grid(l)%outerBox(2))
            weight=max(0.0d0,hi-lo)
            if(l<nGridLevels) then
                if(position>=grid(l)%innerBox(3).and.position<grid(l)%innerBox(4)) &
                    weight=weight-max(0.0d0,min(hi,grid(l)%innerBox(2))-max(lo,grid(l)%innerBox(1)))
            endif
        endif
        weight=max(0.0d0,weight)
    end function section_owned_length

!===================================================================================================
subroutine sample_section(l,fixed,position,axis,tm,um,gradient)
        integer,intent(in) :: l,fixed,axis
        real(kind=8),intent(in) :: position
        real(kind=8),intent(out) :: tm,um,gradient
        integer :: first,k,b,i,j,skin
        real(kind=8) :: q,w(4),dw(4)
        skin=0
        if(l>1) skin=interfaceSkin
        if(axis==1) then
            q=(position-grid(l)%xOffset)/grid(l)%dx+0.5d0
            call section_weights(q-skin,grid(l)%ni-2*skin,first,w,dw)
        else
            q=(position-grid(l)%yOffset)/grid(l)%dx+0.5d0
            call section_weights(q-skin,grid(l)%nj-2*skin,first,w,dw)
        endif
        first=first+skin
        tm=0.0d0
        um=0.0d0
        gradient=0.0d0
        do k=1,4
            if(axis==1) then
                call locate_level_node(l,first+k-1,fixed,b,i,j,.false.)
            else
                call locate_level_node(l,fixed,first+k-1,b,i,j,.false.)
            endif
            if(b==0) error stop 'Middle-section stencil crosses an unstored node'
            tm=tm+w(k)*block(b)%T(i,j)
            if(axis==1) then
                um=um+w(k)*block(b)%u(i,j)
            else
                um=um+w(k)*block(b)%v(i,j)
            endif
            gradient=gradient+dw(k)*block(b)%T(i,j)/grid(l)%dx
        enddo
    end subroutine sample_section

    subroutine calculate_middle_flux(middle)
        real(kind=8),intent(out) :: middle
        integer :: b,l,k,fixed,axis,limit
        real(kind=8) :: position,q,weight,tm,um,gradient,totalLength,fullLength
        middle=0.0d0
        totalLength=0.0d0
        do b=1,nBlocks
            l=block(b)%level
            if(sideHeated) then
                axis=1
                position=0.5d0*nx
                fullLength=dble(ny)
                limit=block(b)%nj
                if(position<block(b)%xOffset.or.position>=block(b)%xOffset+block(b)%ni*block(b)%dx) cycle
            else
                axis=2
                position=0.5d0*ny
                fullLength=dble(nx)
                limit=block(b)%ni
                if(position<block(b)%yOffset.or.position>=block(b)%yOffset+block(b)%nj*block(b)%dx) cycle
            endif
            do k=1,limit
                if(sideHeated) then
                    q=block(b)%yOffset+(k-0.5d0)*block(b)%dx
                    fixed=block(b)%globalJ+k-1
                else
                    q=block(b)%xOffset+(k-0.5d0)*block(b)%dx
                    fixed=block(b)%globalI+k-1
                endif
                weight=section_owned_length(l,q,position,axis)
                if(weight<=0.0d0) cycle
                call sample_section(l,fixed,position,axis,tm,um,gradient)
                middle=middle+(um*tm/diffusivity-gradient)*weight/fullLength
                totalLength=totalLength+weight
            enddo
        enddo
        if(abs(totalLength-fullLength)>1.0d-10*fullLength) error stop 'Middle section does not cover the domain once'
        middle=middle*lengthUnit/(Thot-Tcold)
    end subroutine calculate_middle_flux

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
    open(newunit=k,file=NuReHistoryFile,status='old')   !打开主程序生成的NuRe历史文件
    open(newunit=u,file='NuReCurve_2DOpenaccMultiblock.dat',status='replace')  !创建曲线文件
    write(u,'(a)') '# sample itc_total t_ff NuVolAvg ReVolRMS Nu_hot Nu_cold Nu_middle mass meanT Tmin Tmax rhoMin rhoMax'
    n=0                  !初始化读取过程中的计数
    prevStep=0           !上一条记录的累计步数
    previousTime=0d0     !上一条记录的累计 tff
    interval=0           !两次正式采样之间相隔的细格子步数
    do                   !不断读取，直到遇到文件结尾
        read(k,'(a)',iostat=ios) line  !读取一行
        if(ios<0) exit      !到达文件结尾，退出 do 循环
        if(ios>0) error stop 'History read error'  !发生读取错误，停止程序
        if(len_trim(line)==0 .OR. line(1:1)=='#') cycle !跳过空行和说明行
        read(line,*,iostat=ios) sample,step,values  !将这一行文本解析成数值
        if(ios/=0 .OR. .NOT.all(ieee_is_finite(values))) error stop 'Invalid Nu/Re history' !检查能否按规定读取完整
        n=n+1
        if(n==1) then    !第一次读到数值记录时
            interval=step   !采样间隔，单位为细格子步。
            dt=values(1)    !采样间隔，单位为 tff
            firstTime=values(1) !第一条采样记录的时间
            if(interval<=0 .OR. dt<=0d0) error stop 'Invalid first sample'
        endif
        !能发现采样缺失、重复或步数不符合固定间隔的问题
        if(sample/=n .OR. step/=n*interval) error stop 'Missing or duplicate sample' !核对采样编号和累计步数
        if(values(1)<=previousTime .OR. abs(values(1)-dble(n)*dt)>1d-10*max(1d0,abs(values(1)))) &
            error stop 'History time does not match fixed sampling' !核对累计时间
        previousTime=values(1)
        prevStep=step
        write(u,'(a)') trim(line) !写进曲线文件
    enddo
    close(k)
    close(u)
#ifdef unsteadyFlow
    if(n==0) error stop 'No samples for unsteady statistics'
#endif
    lastTime=previousTime

 !稳态：整理收敛误差曲线，保存最终结果摘要
 !非稳态：计算指定时间窗口内的统计结果
#ifdef steadyFlow
        ! 每次check记录的误差直接用于画图，无需从快照重新定义误差。
        open(newunit=k,file='Convergence_2DOpenaccMultiblock.dat',status='old')
        open(newunit=u,file='ConvergenceCurve_2DOpenaccMultiblock.dat',status='replace')
        write(u,'(a)') '# itc_total t_ff errorU errorT'  !打开稳态误差历史，创建误差曲线文件
        prevStep=-1    !保存上一条记录的累计步数
        do
            read(k,'(a)',iostat=ios) line
            if(ios<0) exit   !遇到结尾退出
            if(ios>0) error stop 'Convergence read error' !读取错误则停止
            if(len_trim(line)==0 .OR. line(1:1)=='#') cycle !跳过空行和以 # 开头的说明行
            read(line,*,iostat=ios) step,values(1:3)
            if(ios/=0 .OR. .NOT.all(ieee_is_finite(values(1:3)))) error stop 'Invalid convergence data'
            if(step<=prevStep) error stop 'Nonmonotone convergence steps' !累计步数严格增加
            prevStep=step
            write(u,'(a)') trim(line)
        enddo
        close(k)
        close(u)
        ! 稳态最终结果直接查看主程序输出的 FinalNuRe_2DOpenaccMultiblock.dat，不再另存副本。
#else     
        if(unsteadyAverageStartTf<0d0 .OR. unsteadyAverageStartTf>=unsteadyAverageMidTf .OR. & 
            unsteadyAverageMidTf>=unsteadyAverageEndTf) error stop 'Invalid averaging window' !先检查时间窗口
        open(newunit=u,file='NuReStatistics_2DOpenaccMultiblock.dat',status='replace') !创建统计文件，检查历史是否覆盖窗口
        write(u,*) 'Window t_ff:',unsteadyAverageStartTf,unsteadyAverageEndTf
        if(firstTime>unsteadyAverageStartTf .OR. lastTime<unsteadyAverageEndTf) then !统计区间不全
            write(u,'(a)') 'INCOMPLETE: the history does not cover the requested window; no final mean.'
        else   !统计区间全了
            call average_window(unsteadyAverageStartTf,unsteadyAverageEndTf,allMean,c0)  !整个统计区间
            call average_window(unsteadyAverageStartTf,unsteadyAverageMidTf,firstMean,c1) !前半个统计区间
            call average_window(unsteadyAverageMidTf,unsteadyAverageEndTf,lastMean,c2)    !后半个统计区间
            write(u,*) 'Covered durations:',c0,c1,c2  !写入统计时长
            write(u,'(a)') '# NuVolAvg Re_rms_space_time Nu_hot Nu_cold Nu_middle'
            write(u,'(a,5ES24.16E3)') 'whole:',allMean   !整个窗口
            write(u,'(a,5ES24.16E3)') 'first:',firstMean !前半窗口
            write(u,'(a,5ES24.16E3)') 'last :',lastMean  !后半窗口
            write(u,'(a,5ES24.16E3)') 'relative half-window difference:', & !前后半窗口的差异
                abs(lastMean-firstMean)/max(abs(allMean),1d-30)
        endif
        close(u)
#endif

    if(processSnapshots) then   !核对结果：验证快照保存、读取和空间积分是否正确。
        open(newunit=k,file=outputIndexFile,status='old')
        open(newunit=u,file='SnapshotNuRe_2DOpenaccMultiblock.dat',status='replace')
        write(u,'(a)') '# itc_total t_ff NuVolAvg ReVolRMS Nu_hot Nu_cold Nu_middle mass meanT Tmin Tmax rhoMin rhoMax'
        do
            read(k,'(a)',iostat=ios) line
            if(ios<0) exit
            if(ios>0) error stop 'Output index read error'
            if(len_trim(line)==0.or.line(1:1)=='#') cycle
            read(line,*,iostat=ios) step,name
            if(ios/=0) error stop 'Invalid output index'
            if(index(name,snapshotFilePrefix//'-')/=1) cycle
            ! 索引保存完整的<tff>tff-<累计细步数>itc.bin文件名，直接打开，不从文件名反推时间或步数。
            call read_snapshot(trim(name))
            if(totalItc/=step) error stop 'Snapshot/index step mismatch'
            call calculate_NuRe(values(2),values(3),values(4),values(5),values(6),values(7), &
                values(8),values(9),values(10),values(11),values(12))
            values(1)=totalTff
            write(u,'(I0,1X,12(ES24.16E3,1X))') totalItc,values
        enddo
        close(k)
        close(u)
    endif
    write(*,'(a)') 'Postprocessing complete; solver fields and histories were not modified.'
end program postprocess

!===================================================================================================
subroutine average_window(t0,t1,result,coverage)
    use postdata,only: NuReHistoryFile
    implicit none
    real(kind=8),intent(in) :: t0,t1
    real(kind=8),intent(out) :: result(5),coverage
    real(kind=8) :: prev(12),values(12),a,b,left(5),right(5),ra(5),rb(5)
    integer :: k,ios,sample,step
    logical :: hasPrev
    character(len=512) :: line


    result=0d0  !返回五项统计结果
    coverage=0d0 !实际参与积分的时间长度
    hasPrev=.false.  !是否已经读到上一条记录
    open(newunit=k,file=NuReHistoryFile,status='old')
    do
        read(k,'(a)',iostat=ios) line
        if(ios<0) exit  !遇到文件结尾退出
        if(ios>0) error stop 'Averaging read error' 
        if(len_trim(line)==0 .OR. line(1:1)=='#') cycle !空行和说明行跳过
        read(line,*,iostat=ios) sample,step,values
        if(ios/=0) error stop 'Invalid averaging row'
        if(hasPrev) then   !两个时刻之间才能构成一段时间，因此第一条不计算积分，只保存到 prev
            a=max(t0,prev(1))  !截取落在统计窗口内的时间段,前一时刻
            b=min(t1,values(1))   !截取落在统计窗口内的时间段,当前时刻
            if(b>a) then
                left=prev(2:6)  !上一条记录的 Nu、Re、三个截面 Nu
                right=values(2:6)  !当前一条记录的对应数值
                ! 保持时空均方根定义：对瞬时空间RMS的平方积分，再开方。
                left(2)=left(2)**2
                right(2)=right(2)**2
                !a、b 不一定正好等于已有采样时刻
                ra=left+(right-left)*(a-prev(1))/(values(1)-prev(1)) !计算截取后两个端点上的值
                rb=left+(right-left)*(b-prev(1))/(values(1)-prev(1)) !计算截取后两个端点上的值
                result=result+0.5d0*(ra+rb)*(b-a) !用梯形公式累计积分和时间长度
                coverage=coverage+b-a             !用梯形公式累计积分和时间长度
            endif
        endif
        prev=values
        hasPrev=.true.
    enddo
    close(k)
    if(abs(coverage-(t1-t0))>1d-10*max(1d0,t1-t0)) error stop 'Incomplete average window' !检查实际积分时长是否等于目标窗口
    result=result/coverage        !积分除以时长，RMS Re 再开方
    result(2)=sqrt(max(0d0,result(2)))
end subroutine average_window

!===================================================================================================
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

!===================================================================================================
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
