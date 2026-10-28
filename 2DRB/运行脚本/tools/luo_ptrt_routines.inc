!===================================================================================================
! 子程序: ptrt_moment_collision
! 作用: 单节点闭合 P-TRT 碰撞；顺序沿用母版 (delta_rho,e,epsilon,jx,Qx,jy,Qy,pxx,pxy)。
! 说明: ux/uy=(j+F/2)/rho0；Fx/Fy 是 dt=1 的动量增量密度；本子程序已包含完整外力。
!===================================================================================================
subroutine ptrt_moment_collision(m,rhoReference,ux,uy,forceX,forceY,rateNu,rateQ,mp)
    implicit none
    !$acc routine seq
    real(kind=8), intent(in) :: m(0:8),rhoReference,ux,uy,forceX,forceY,rateNu,rateQ
    real(kind=8), intent(out) :: mp(0:8)
    real(kind=8) :: eqE,eqXX,eqXY,eqQX,eqQY,aNu,bQ,uDotF
    aNu=1.0d0-0.5d0*rateNu
    bQ=1.0d0-0.5d0*rateQ
    uDotF=ux*forceX+uy*forceY
    eqE=-2.0d0*m(0)+3.0d0*rhoReference*(ux*ux+uy*uy)
    eqXX=rhoReference*(ux*ux-uy*uy)
    eqXY=rhoReference*ux*uy
    eqQX=rhoReference*ux*(3.0d0*uy*uy-1.0d0)
    eqQY=rhoReference*uy*(3.0d0*ux*ux-1.0d0)
    mp(0)=m(0)
    mp(3)=m(3)+forceX
    mp(5)=m(5)+forceY
    mp(1)=m(1)-rateNu*(m(1)-eqE)+6.0d0*aNu*uDotF
    mp(7)=m(7)-rateNu*(m(7)-eqXX)+2.0d0*aNu*(ux*forceX-uy*forceY)
    mp(8)=m(8)-rateNu*(m(8)-eqXY)+aNu*(ux*forceY+uy*forceX)
    mp(4)=m(4)-rateQ*(m(4)-eqQX)-bQ*forceX
    mp(6)=m(6)-rateQ*(m(6)-eqQY)-bQ*forceY
    ! epsilon 并非独立的纯 ghost。最后恢复它，严格令 delta_rho+e*+epsilon*=0。
    mp(2)=-mp(0)-mp(1)
end subroutine ptrt_moment_collision
!===================================================================================================

!===================================================================================================
! 子程序: refresh_ptrt_velocity
! 作用: 温度更新后，用 F^(n+1) 重建物理速度；只重建，不改 f、不追加动量冲量。
!===================================================================================================
subroutine refresh_ptrt_velocity()
    use commondata
    implicit none
    integer(kind=4) :: i,j
    !$acc parallel loop gang vector collapse(2) default(none) present(f,T,u,v,Fx,Fy) async(1)
    do j=1,ny
        do i=1,nx
            Fx(i,j)=0.0d0
            Fy(i,j)=rho0*gBeta*(T(i,j)-Tref)
            u(i,j)=(f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)+0.5d0*Fx(i,j))/rho0
            v(i,j)=(f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)+0.5d0*Fy(i,j))/rho0
        enddo
    enddo
end subroutine refresh_ptrt_velocity
!===================================================================================================

!===================================================================================================
! 子程序: monitor_ptrt_collision
! 作用: 在迁移前稀疏检查本地质量、冲量、ghost；不把迁移后的 ghost 当成过滤误差。
!===================================================================================================
subroutine monitor_ptrt_collision()
    use commondata
    implicit none
    integer(kind=4) :: i,j,k,unitNo,badCount,badI,badJ
    real(kind=8) :: massError,forceError,ghostPre,ghostPost,qMax,uMax,rhoMax,tMin,tMax,thermalMax
    real(kind=8) :: r,rp,e,ep,eps,epsp,jx,jy,jxp,jyp,qx,qy,thermalMoment
    real(kind=8), parameter :: finiteLimit=huge(1.0d0)
    logical :: exists
    massError=0.0d0; forceError=0.0d0; ghostPre=0.0d0; ghostPost=0.0d0
    qMax=0.0d0; uMax=0.0d0; rhoMax=0.0d0; tMin=huge(1.0d0); tMax=-huge(1.0d0)
    thermalMax=0.0d0; badCount=0
    !$acc wait(1)
    !$acc parallel loop gang vector collapse(2) default(none) present(f,f_post,g,u,v,T,Fx,Fy) &
    !$acc& private(k,r,rp,e,ep,eps,epsp,jx,jy,jxp,jyp,qx,qy,thermalMoment) &
    !$acc& reduction(max:massError,forceError,ghostPre,ghostPost,qMax,uMax,rhoMax,tMax,thermalMax) &
    !$acc& reduction(min:tMin) reduction(+:badCount)
    do j=1,ny
        do i=1,nx
            r=0.0d0; rp=0.0d0
            do k=0,8
                r=r+f(i,j,k); rp=rp+f_post(i,j,k)
                if(.not.(abs(f_post(i,j,k)).LE.finiteLimit)) badCount=badCount+1
                if(.not.(abs(f(i,j,k)).LE.finiteLimit)) badCount=badCount+1
            enddo
            do k=0,4
                if(.not.(abs(g(i,j,k)).LE.finiteLimit)) badCount=badCount+1
            enddo
            if(.not.(abs(T(i,j)).LE.finiteLimit)) badCount=badCount+1
            e=-4*f(i,j,0)-sum(f(i,j,1:4))+2*sum(f(i,j,5:8))
            ep=-4*f_post(i,j,0)-sum(f_post(i,j,1:4))+2*sum(f_post(i,j,5:8))
            eps=4*f(i,j,0)-2*sum(f(i,j,1:4))+sum(f(i,j,5:8))
            epsp=4*f_post(i,j,0)-2*sum(f_post(i,j,1:4))+sum(f_post(i,j,5:8))
            jx=f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)
            jy=f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)
            jxp=f_post(i,j,1)-f_post(i,j,3)+f_post(i,j,5)-f_post(i,j,6)-f_post(i,j,7)+f_post(i,j,8)
            jyp=f_post(i,j,2)-f_post(i,j,4)+f_post(i,j,5)+f_post(i,j,6)-f_post(i,j,7)-f_post(i,j,8)
            qx=-2*f(i,j,1)+2*f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8) &
                -rho0*u(i,j)*(3*v(i,j)*v(i,j)-1)
            qy=-2*f(i,j,2)+2*f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8) &
                -rho0*v(i,j)*(3*u(i,j)*u(i,j)-1)
            thermalMoment=-4*g(i,j,0)+sum(g(i,j,1:4))-paraA*T(i,j)
            massError=max(massError,abs(rp-r))
            forceError=max(forceError,abs(jxp-jx-Fx(i,j)),abs(jyp-jy-Fy(i,j)))
            ghostPre=max(ghostPre,abs(r+e+eps)/9.0d0)
            ghostPost=max(ghostPost,abs(rp+ep+epsp))
            qMax=max(qMax,sqrt(qx*qx+qy*qy))
            uMax=max(uMax,sqrt(u(i,j)**2+v(i,j)**2)); rhoMax=max(rhoMax,abs(r))
            tMin=min(tMin,T(i,j)); tMax=max(tMax,T(i,j))
            thermalMax=max(thermalMax,abs(thermalMoment),abs(g(i,j,1)-g(i,j,2)+g(i,j,3)-g(i,j,4)))
        enddo
    enddo
    inquire(file='PTRTCollisionHistory.dat',exist=exists)
    open(newunit=unitNo,file='PTRTCollisionHistory.dat',status='unknown',position='append')
    if(.not.exists) write(unitNo,'(A)') '# tff step E_mass E_force E_ghost Sghost_pre Qneq_max speed_max '// &
        'delta_rho_max Tmin Tmax thermal_high_neq_max nonfinite_population_count'
    write(unitNo,'(ES24.16E3,1X,I0,1X,10(ES24.16E3,1X),I0)') &
        dble(restartItcOffset+itc)/timeUnit,restartItcOffset+itc,massError,forceError,ghostPost,ghostPre, &
        qMax,uMax,rhoMax,tMin,tMax,thermalMax,badCount
    close(unitNo)
    if(badCount.GT.0) then
        !$acc update self(f,f_post,g,T)
        badI=-1; badJ=-1
        find_bad: do j=1,ny
            do i=1,nx
                do k=0,8
                    if(.not.(abs(f_post(i,j,k)).LE.finiteLimit).OR. &
                       .not.(abs(f(i,j,k)).LE.finiteLimit)) then
                        badI=i; badJ=j
                        write(*,*) 'PTRT_NONFINITE flow: step/i/j/k=',restartItcOffset+itc,i,j,k
                        exit find_bad
                    endif
                enddo
                do k=0,4
                    if(.not.(abs(g(i,j,k)).LE.finiteLimit)) then
                        badI=i; badJ=j
                        write(*,*) 'PTRT_NONFINITE thermal: step/i/j/k=',restartItcOffset+itc,i,j,k
                        exit find_bad
                    endif
                enddo
            enddo
        enddo find_bad
        write(*,*) 'First detected sampled failure; checkpoint is not overwritten:',badI,badJ
        error stop 86
    endif
end subroutine monitor_ptrt_collision
!===================================================================================================

!===================================================================================================
! 子程序: monitor_ptrt_wall_flux
! 作用: 保留母版 RB 壁面梯度公式、正号与归一化；同原 Nu/Re 每 1 t_ff 采样。
!===================================================================================================
subroutine monitor_ptrt_wall_flux()
    use commondata
    implicit none
    integer(kind=4) :: i,unitNo
    real(kind=8) :: nuBottom,nuTop
    logical :: exists
    nuBottom=0.0d0; nuTop=0.0d0
    !$acc wait(1)
    !$acc parallel loop default(none) present(T) reduction(+:nuBottom,nuTop)
    do i=1,nx
        nuBottom=nuBottom+(8*Thot-9*T(i,1)+T(i,2))*lengthUnit/(3*deltaT)
        nuTop=nuTop+(-8*Tcold+9*T(i,ny)-T(i,ny-1))*lengthUnit/(3*deltaT)
    enddo
    inquire(file='PTRTWallHistory.dat',exist=exists)
    open(newunit=unitNo,file='PTRTWallHistory.dat',status='unknown',position='append')
    if(.not.exists) write(unitNo,'(A)') '# tff Nu_bottom Nu_top Nu_volume Re_rms'
    write(unitNo,'(5(ES24.16E3,1X))') dble(restartItcOffset+itc)/timeUnit, &
        nuBottom/dble(nx),nuTop/dble(nx),NuVolAvg(currentRunStatisticSampleCount), &
        ReVolAvg(currentRunStatisticSampleCount)
    close(unitNo)
end subroutine monitor_ptrt_wall_flux
!===================================================================================================
