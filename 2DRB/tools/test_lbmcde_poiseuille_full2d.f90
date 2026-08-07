program test_lbmcde_poiseuille_full2d
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    integer, parameter :: dp=kind(1.0d0)
    real(dp), parameter :: cs2=1.0_dp/3.0_dp
    real(dp), parameter :: tau_base=0.62_dp
    real(dp), parameter :: sigma_base=tau_base-0.5_dp
    real(dp), parameter :: omega_even=1.0_dp/tau_base
    real(dp), parameter :: chi_b=0.0_dp
    integer, parameter :: ex(0:8)=[0,1,0,-1,0,1,-1,-1,1]
    integer, parameter :: ey(0:8)=[0,0,1,0,-1,1,1,-1,-1]
    integer, parameter :: opposite(0:8)=[0,3,4,1,2,7,8,5,6]
    real(dp), parameter :: weight(0:8)=[4.0_dp/9.0_dp, &
        1.0_dp/9.0_dp,1.0_dp/9.0_dp,1.0_dp/9.0_dp,1.0_dp/9.0_dp, &
        1.0_dp/36.0_dp,1.0_dp/36.0_dp,1.0_dp/36.0_dp,1.0_dp/36.0_dp]

    integer :: nx, ny, steps, check_span, policy
    character(len=32) :: arg, policy_name
    real(dp) :: viscosity, chi_s, sigma_effective, sigma_odd, omega_odd
    real(dp) :: target_umax, accel_x, reynolds, ma_input, convergence_tolerance
    real(dp), allocatable :: f(:,:,:), post(:,:,:), next_f(:,:,:)
    real(dp), allocatable :: profile(:)
    real(dp), allocatable :: u_current(:,:),v_current(:,:)
    real(dp), allocatable :: u_previous_step(:,:),v_previous_step(:,:)
    real(dp), allocatable :: u_previous_check(:,:),v_previous_check(:,:)

    call get_command_argument(1,arg)
    if(len_trim(arg)==0) error stop &
        'usage: exe NY effective|bgk [max_steps] [check_span] [Re] [Ma] [tolerance]'
    read(arg,*) ny
    call get_command_argument(2,policy_name)
    select case(trim(policy_name))
    case('effective')
        policy=2
    case('bgk')
        policy=3
    case default
        error stop 'policy must be effective or bgk'
    end select
    steps=500000
    call get_command_argument(3,arg)
    if(len_trim(arg)>0) read(arg,*) steps
    check_span=max(1000,steps/10)
    call get_command_argument(4,arg)
    if(len_trim(arg)>0) read(arg,*) check_span
    reynolds=1000.0_dp
    call get_command_argument(5,arg)
    if(len_trim(arg)>0) read(arg,*) reynolds
    ma_input=0.1_dp
    call get_command_argument(6,arg)
    if(len_trim(arg)>0) read(arg,*) ma_input
    convergence_tolerance=1.0e-8_dp
    call get_command_argument(7,arg)
    if(len_trim(arg)>0) read(arg,*) convergence_tolerance

    nx=2*ny
    target_umax=ma_input*sqrt(cs2)
    viscosity=target_umax*real(ny,dp)/reynolds
    chi_s=1.0_dp-viscosity/(cs2*sigma_base)
    sigma_effective=viscosity/cs2
    select case(policy)
    case(2)
        sigma_odd=3.0_dp/(16.0_dp*sigma_effective)
    case(3)
        sigma_odd=sigma_base
    end select
    omega_odd=1.0_dp/(0.5_dp+sigma_odd)
    accel_x=8.0_dp*viscosity*target_umax/real(ny*ny,dp)

    allocate(f(0:8,nx,ny),post(0:8,nx,ny),next_f(0:8,nx,ny))
    allocate(profile(ny))
    allocate(u_current(nx,ny),v_current(nx,ny))
    allocate(u_previous_step(nx,ny),v_previous_step(nx,ny))
    allocate(u_previous_check(nx,ny),v_previous_check(nx,ny))
    f=0.0_dp
    u_current=0.0_dp; v_current=0.0_dp
    u_previous_step=0.0_dp; v_previous_step=0.0_dp
    u_previous_check=0.0_dp; v_previous_check=0.0_dp

    !$acc enter data copyin(f) create(post,next_f,profile,u_current,v_current, &
    !$acc u_previous_step,v_previous_step)
    call march(f,post,next_f,profile,u_current,v_current,u_previous_step,v_previous_step, &
        u_previous_check,v_previous_check,nx,ny,steps,check_span, &
        accel_x,viscosity,chi_s,omega_odd,target_umax,ma_input, &
        convergence_tolerance,policy_name)
    !$acc exit data delete(f,post,next_f,profile,u_current,v_current, &
    !$acc u_previous_step,v_previous_step)

contains

    subroutine march(f,post,next_f,profile,u_current,v_current,u_previous_step, &
        v_previous_step,u_previous_check,v_previous_check,nx,ny,steps,check_span, &
        accel_x,viscosity,chi_s,omega_odd,target_umax,ma_input, &
        convergence_tolerance,policy_name)
        integer, intent(in) :: nx,ny,steps,check_span
        real(dp), intent(inout) :: f(0:8,nx,ny),post(0:8,nx,ny),next_f(0:8,nx,ny)
        real(dp), intent(inout) :: profile(ny)
        real(dp), intent(inout) :: u_current(nx,ny),v_current(nx,ny)
        real(dp), intent(inout) :: u_previous_step(nx,ny),v_previous_step(nx,ny)
        real(dp), intent(inout) :: u_previous_check(nx,ny),v_previous_check(nx,ny)
        real(dp), intent(in) :: accel_x,viscosity,chi_s,omega_odd,target_umax,ma_input
        real(dp), intent(in) :: convergence_tolerance
        character(len=*), intent(in) :: policy_name
        integer :: iteration,checkpoint_count
        real(dp) :: relative_change,relative_change_window,relative_change_one_step
        logical :: converged

        relative_change=huge(1.0_dp)
        relative_change_window=huge(1.0_dp)
        relative_change_one_step=huge(1.0_dp)
        checkpoint_count=0
        converged=.false.
        do iteration=1,steps
            if(mod(iteration,check_span)==0) then
                ! f is the state at n-1 here.  Preserve its full velocity field so that
                ! a period-2 orbit cannot alias through an even checkpoint interval.
                call velocity_field(f,u_previous_step,v_previous_step,nx,ny,accel_x)
                !$acc update self(u_previous_step,v_previous_step)
            endif
            call collide(f,post,nx,ny,accel_x,viscosity,chi_s,omega_odd)
            call stream_bounce(post,next_f,nx,ny)
            call copy_state(next_f,f,nx,ny)
            if(mod(iteration,1000)==0) then
                call average_profile(f,profile,nx,ny,accel_x)
                !$acc update self(profile)
                if(any(.not.ieee_is_finite(profile))) then
                    write(*,'(a)') 'FAIL_NONFINITE'
                    write(*,'(a,a)') 'policy=',trim(policy_name)
                    write(*,'(a,i0)') 'nx=',nx
                    write(*,'(a,i0)') 'ny=',ny
                    write(*,'(a,i0)') 'first_detected_iteration=',iteration
                    write(*,'(a,es24.16)') 'reynolds=',reynolds
                    write(*,'(a,es24.16)') 'viscosity=',viscosity
                    write(*,'(a,es24.16)') 'omega_odd=',omega_odd
                    error stop 70
                endif
            endif
            if(mod(iteration,check_span)==0) then
                call velocity_field(f,u_current,v_current,nx,ny,accel_x)
                !$acc update self(u_current,v_current)
                if(any(.not.ieee_is_finite(u_current)) .or. &
                    any(.not.ieee_is_finite(v_current))) then
                    write(*,'(a)') 'FAIL_NONFINITE_FULL_VELOCITY'
                    write(*,'(a,a)') 'policy=',trim(policy_name)
                    write(*,'(a,i0)') 'first_detected_iteration=',iteration
                    error stop 70
                endif
                relative_change_one_step=velocity_relative_change(u_current,v_current, &
                    u_previous_step,v_previous_step,nx,ny)
                checkpoint_count=checkpoint_count+1
                if(checkpoint_count>1) then
                    relative_change_window=velocity_relative_change(u_current,v_current, &
                        u_previous_check,v_previous_check,nx,ny)
                    relative_change=max(relative_change_window,relative_change_one_step)
                    write(*,'(a,i0,a,es16.8,a,es16.8,a,es16.8)') &
                        'CHECK iteration=',iteration, &
                        ' relative_change_window=',relative_change_window, &
                        ' relative_change_one_step=',relative_change_one_step, &
                        ' relative_change_combined=',relative_change
                    if(relative_change<=convergence_tolerance) then
                        converged=.true.
                        exit
                    endif
                endif
                u_previous_check=u_current
                v_previous_check=v_current
            endif
        enddo
        call average_profile(f,profile,nx,ny,accel_x)
        !$acc update self(profile,f)
        call report_case(f,profile,nx,ny,steps,check_span,accel_x,viscosity,chi_s, &
            omega_odd,target_umax,ma_input,relative_change,relative_change_window, &
            relative_change_one_step,convergence_tolerance, &
            min(iteration,steps),converged,policy_name)
        if(.not.converged) error stop 71
    end subroutine march

    subroutine collide(f,post,nx,ny,accel_x,viscosity,chi_s,omega_odd)
        integer, intent(in) :: nx,ny
        real(dp), intent(in) :: f(0:8,nx,ny)
        real(dp), intent(out) :: post(0:8,nx,ny)
        real(dp), intent(in) :: accel_x,viscosity,chi_s,omega_odd
        integer :: i,j,a,b,p
        integer, parameter :: pa(4)=[1,2,5,6], pb(4)=[3,4,7,8]
        real(dp) :: rho_p,u,v,uu,eu,feq(0:8),neq(0:8),raw(0:8)
        real(dp) :: neq_xx,neq_xy,neq_yy,bulk_viscosity,denom_diag,denom_shear
        real(dp) :: denom_div,coeff_trace,uf,strain_div,strain_xx,strain_yy,strain_xy
        real(dp) :: corr_xx,corr_yy,corr_xy,hermite,direction_force
        real(dp) :: d_even,d_odd,r_even,r_odd,pref_even,pref_odd

        bulk_viscosity=sigma_base*(1.0_dp-chi_b)*cs2
        denom_diag=2.0_dp*viscosity+cs2
        denom_shear=4.0_dp*viscosity+2.0_dp*cs2
        denom_div=2.0_dp*bulk_viscosity+cs2
        coeff_trace=(2.0_dp*viscosity-2.0_dp*bulk_viscosity)/(2.0_dp*denom_diag)
        pref_even=1.0_dp-0.5_dp*omega_even
        pref_odd=1.0_dp-0.5_dp*omega_odd

        !$acc parallel loop collapse(2) present(f,post) private(i,j,a,b,p,rho_p,u,v,uu,eu, &
        !$acc feq,neq,raw,neq_xx,neq_xy,neq_yy,uf,strain_div,strain_xx,strain_yy,strain_xy, &
        !$acc corr_xx,corr_yy,corr_xy,hermite,direction_force,d_even,d_odd,r_even,r_odd)
        !$omp parallel do collapse(2) default(shared) private(i,j,a,b,p,rho_p,u,v,uu,eu, &
        !$omp feq,neq,raw,neq_xx,neq_xy,neq_yy,uf,strain_div,strain_xx,strain_yy,strain_xy, &
        !$omp corr_xx,corr_yy,corr_xy,hermite,direction_force,d_even,d_odd,r_even,r_odd)
        do j=1,ny
            do i=1,nx
                rho_p=sum(f(:,i,j))
                u=0.5_dp*accel_x
                v=0.0_dp
                do a=0,8
                    u=u+real(ex(a),dp)*f(a,i,j)
                    v=v+real(ey(a),dp)*f(a,i,j)
                enddo
                uu=u*u+v*v
                do a=0,8
                    eu=real(ex(a),dp)*u+real(ey(a),dp)*v
                    feq(a)=weight(a)*(rho_p+eu/cs2+0.5_dp*eu*eu/(cs2*cs2)-0.5_dp*uu/cs2)
                    neq(a)=f(a,i,j)-feq(a)
                enddo
                neq_xx=0.0_dp; neq_xy=0.0_dp; neq_yy=0.0_dp
                do a=0,8
                    neq_xx=neq_xx+real(ex(a)*ex(a),dp)*neq(a)
                    neq_xy=neq_xy+real(ex(a)*ey(a),dp)*neq(a)
                    neq_yy=neq_yy+real(ey(a)*ey(a),dp)*neq(a)
                enddo
                uf=u*accel_x
                strain_div=-(neq_xx+neq_yy+uf)/denom_div
                strain_xx=-(neq_xx+uf)/denom_diag+coeff_trace*strain_div
                strain_yy=-neq_yy/denom_diag+coeff_trace*strain_div
                strain_xy=-(2.0_dp*neq_xy+v*accel_x)/denom_shear
                corr_xx=chi_s*strain_xx+0.5_dp*(chi_b-chi_s)*strain_div
                corr_yy=chi_s*strain_yy+0.5_dp*(chi_b-chi_s)*strain_div
                corr_xy=chi_s*strain_xy
                do a=0,8
                    eu=real(ex(a),dp)*u+real(ey(a),dp)*v
                    direction_force=real(ex(a),dp)*accel_x
                    hermite=(real(ex(a)*ex(a),dp)/cs2-1.0_dp)*corr_xx + &
                        (real(ey(a)*ey(a),dp)/cs2-1.0_dp)*corr_yy + &
                        2.0_dp*real(ex(a)*ey(a),dp)*corr_xy/cs2
                    raw(a)=weight(a)*(direction_force/cs2+eu*direction_force/(cs2*cs2) &
                        -uf/cs2+hermite)
                enddo
                post(0,i,j)=f(0,i,j)-omega_even*(f(0,i,j)-feq(0))+pref_even*raw(0)
                do p=1,4
                    a=pa(p); b=pb(p)
                    d_even=0.5_dp*(f(a,i,j)-feq(a)+f(b,i,j)-feq(b))
                    d_odd=0.5_dp*(f(a,i,j)-feq(a)-f(b,i,j)+feq(b))
                    r_even=0.5_dp*(raw(a)+raw(b))
                    r_odd=0.5_dp*(raw(a)-raw(b))
                    post(a,i,j)=f(a,i,j)-omega_even*d_even-omega_odd*d_odd + &
                        pref_even*r_even+pref_odd*r_odd
                    post(b,i,j)=f(b,i,j)-omega_even*d_even+omega_odd*d_odd + &
                        pref_even*r_even-pref_odd*r_odd
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine collide

    subroutine stream_bounce(post,next_f,nx,ny)
        integer, intent(in) :: nx,ny
        real(dp), intent(in) :: post(0:8,nx,ny)
        real(dp), intent(out) :: next_f(0:8,nx,ny)
        integer :: i,j,a,ip,jp
        !$acc parallel loop collapse(3) present(post,next_f) private(i,j,a,ip,jp)
        !$omp parallel do collapse(3) default(shared) private(i,j,a,ip,jp)
        do j=1,ny
            do i=1,nx
                do a=0,8
                    ip=1+modulo(i-1-ex(a),nx)
                    jp=j-ey(a)
                    if(jp<1 .or. jp>ny) then
                        next_f(a,i,j)=post(opposite(a),i,j)
                    else
                        next_f(a,i,j)=post(a,ip,jp)
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine stream_bounce

    subroutine copy_state(source,destination,nx,ny)
        integer, intent(in) :: nx,ny
        real(dp), intent(in) :: source(0:8,nx,ny)
        real(dp), intent(out) :: destination(0:8,nx,ny)
        integer :: i,j,a
        !$acc parallel loop collapse(3) present(source,destination)
        !$omp parallel do collapse(3) default(shared) private(i,j,a)
        do j=1,ny
            do i=1,nx
                do a=0,8
                    destination(a,i,j)=source(a,i,j)
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine copy_state

    subroutine velocity_field(f,u_field,v_field,nx,ny,accel_x)
        integer, intent(in) :: nx,ny
        real(dp), intent(in) :: f(0:8,nx,ny),accel_x
        real(dp), intent(out) :: u_field(nx,ny),v_field(nx,ny)
        integer :: i,j,a
        real(dp) :: ux,uy
        !$acc parallel loop collapse(2) present(f,u_field,v_field) private(i,j,a,ux,uy)
        !$omp parallel do collapse(2) default(shared) private(i,j,a,ux,uy)
        do j=1,ny
            do i=1,nx
                ux=0.5_dp*accel_x
                uy=0.0_dp
                do a=0,8
                    ux=ux+real(ex(a),dp)*f(a,i,j)
                    uy=uy+real(ey(a),dp)*f(a,i,j)
                enddo
                u_field(i,j)=ux
                v_field(i,j)=uy
            enddo
        enddo
        !$omp end parallel do
    end subroutine velocity_field

    function velocity_relative_change(u_new,v_new,u_old,v_old,nx,ny) result(value)
        integer, intent(in) :: nx,ny
        real(dp), intent(in) :: u_new(nx,ny),v_new(nx,ny),u_old(nx,ny),v_old(nx,ny)
        real(dp) :: value,numerator,denominator
        numerator=sqrt(sum((u_new-u_old)**2+(v_new-v_old)**2))
        denominator=sqrt(sum(u_new**2+v_new**2))
        value=numerator/max(denominator,tiny(1.0_dp))
    end function velocity_relative_change

    subroutine average_profile(f,profile,nx,ny,accel_x)
        integer, intent(in) :: nx,ny
        real(dp), intent(in) :: f(0:8,nx,ny),accel_x
        real(dp), intent(out) :: profile(ny)
        integer :: i,j,a
        real(dp) :: sum_velocity
        !$acc parallel loop gang present(f,profile) private(i,j,a,sum_velocity)
        do j=1,ny
            sum_velocity=0.0_dp
            !$acc loop vector reduction(+:sum_velocity)
            do i=1,nx
                sum_velocity=sum_velocity+0.5_dp*accel_x
                do a=0,8
                    sum_velocity=sum_velocity+real(ex(a),dp)*f(a,i,j)
                enddo
            enddo
            profile(j)=sum_velocity/real(nx,dp)
        enddo
    end subroutine average_profile

    subroutine report_case(f,profile,nx,ny,steps,check_span,accel_x,viscosity,chi_s, &
        omega_odd,target_umax,ma_input,relative_change,relative_change_window, &
        relative_change_one_step,convergence_tolerance, &
        actual_steps,converged,policy_name)
        integer, intent(in) :: nx,ny,steps,check_span
        real(dp), intent(in) :: f(0:8,nx,ny),profile(ny),accel_x,viscosity,chi_s
        real(dp), intent(in) :: omega_odd,target_umax,ma_input,relative_change
        real(dp), intent(in) :: relative_change_window,relative_change_one_step
        real(dp), intent(in) :: convergence_tolerance
        integer, intent(in) :: actual_steps
        logical, intent(in) :: converged
        character(len=*), intent(in) :: policy_name
        integer :: i,j,a,i_left,i_right
        real(dp) :: y,exact,l2_num,l2_den,max_x_nonuniform,uij,rho_p,max_rho
        real(dp) :: xmean,sigma_effective_local,sigma_odd_local
        real(dp) :: u_left,u_right,u_mid,local_error,mean_error,error_variance
        real(dp) :: max_abs_error

        i_left=nx/2
        i_right=i_left+1
        l2_num=0.0_dp; l2_den=0.0_dp
        mean_error=0.0_dp; max_abs_error=0.0_dp
        do j=1,ny
            y=real(j,dp)-0.5_dp
            exact=accel_x*y*(real(ny,dp)-y)/(2.0_dp*viscosity)
            u_left=0.5_dp*accel_x
            u_right=0.5_dp*accel_x
            do a=0,8
                u_left=u_left+real(ex(a),dp)*f(a,i_left,j)
                u_right=u_right+real(ex(a),dp)*f(a,i_right,j)
            enddo
            u_mid=0.5_dp*(u_left+u_right)
            local_error=u_mid-exact
            l2_num=l2_num+local_error**2
            l2_den=l2_den+exact**2
            mean_error=mean_error+local_error
            max_abs_error=max(max_abs_error,abs(local_error))
        enddo
        mean_error=mean_error/real(ny,dp)
        error_variance=0.0_dp
        do j=1,ny
            y=real(j,dp)-0.5_dp
            exact=accel_x*y*(real(ny,dp)-y)/(2.0_dp*viscosity)
            u_left=0.5_dp*accel_x
            u_right=0.5_dp*accel_x
            do a=0,8
                u_left=u_left+real(ex(a),dp)*f(a,i_left,j)
                u_right=u_right+real(ex(a),dp)*f(a,i_right,j)
            enddo
            u_mid=0.5_dp*(u_left+u_right)
            local_error=u_mid-exact
            error_variance=error_variance+(local_error-mean_error)**2
        enddo
        error_variance=sqrt(error_variance/real(ny,dp))

        max_x_nonuniform=0.0_dp; max_rho=0.0_dp
        do j=1,ny
            xmean=profile(j)
            do i=1,nx
                uij=0.5_dp*accel_x
                rho_p=0.0_dp
                do a=0,8
                    uij=uij+real(ex(a),dp)*f(a,i,j)
                    rho_p=rho_p+f(a,i,j)
                enddo
                max_x_nonuniform=max(max_x_nonuniform,abs(uij-xmean))
                max_rho=max(max_rho,abs(rho_p))
            enddo
        enddo
        sigma_effective_local=viscosity/cs2
        sigma_odd_local=1.0_dp/omega_odd-0.5_dp
        write(*,'(a)') 'RESULT'
        write(*,'(a,a)') 'policy=',trim(policy_name)
        write(*,'(a,i0)') 'nx=',nx
        write(*,'(a,i0)') 'ny=',ny
        write(*,'(a,i0)') 'max_steps=',steps
        write(*,'(a,i0)') 'actual_steps=',actual_steps
        write(*,'(a,i0)') 'check_span=',check_span
        write(*,'(a,l1)') 'converged=',converged
        write(*,'(a,es24.16)') 'convergence_tolerance=',convergence_tolerance
        write(*,'(a,es24.16)') 'ma=',ma_input
        write(*,'(a,es24.16)') 'reynolds=',reynolds
        write(*,'(a,es24.16)') 'viscosity=',viscosity
        write(*,'(a,es24.16)') 'chi_s=',chi_s
        write(*,'(a,es24.16)') 'sigma_base=',sigma_base
        write(*,'(a,es24.16)') 'sigma_effective=',sigma_effective_local
        write(*,'(a,es24.16)') 'omega_even=',omega_even
        write(*,'(a,es24.16)') 'omega_odd=',omega_odd
        write(*,'(a,es24.16)') 'sigma_effective_times_sigma_odd=', &
            sigma_effective_local*sigma_odd_local
        write(*,'(a,es24.16)') 'target_umax=',target_umax
        write(*,'(a,es24.16)') 'acceleration_x=',accel_x
        write(*,'(a,es24.16)') 'relative_change_window=',relative_change_window
        write(*,'(a,es24.16)') 'relative_change_one_step=',relative_change_one_step
        write(*,'(a,es24.16)') 'relative_change_combined=',relative_change
        write(*,'(a,es24.16)') 'midline_relative_l2_error=',sqrt(l2_num/l2_den)
        write(*,'(a,es24.16)') 'midline_mean_local_error_over_umax=',mean_error/target_umax
        write(*,'(a,es24.16)') 'midline_max_abs_local_error_over_umax=', &
            max_abs_error/target_umax
        write(*,'(a,es24.16)') 'midline_local_error_std_over_umax=', &
            error_variance/target_umax
        write(*,'(a,es24.16)') 'max_x_nonuniformity_over_umax=',max_x_nonuniform/target_umax
        write(*,'(a,es24.16)') 'max_density_perturbation=',max_rho
        write(*,'(a)') 'MIDLINE_PROFILE_BEGIN'
        do j=1,ny
            y=real(j,dp)-0.5_dp
            exact=accel_x*y*(real(ny,dp)-y)/(2.0_dp*viscosity)
            u_left=0.5_dp*accel_x
            u_right=0.5_dp*accel_x
            do a=0,8
                u_left=u_left+real(ex(a),dp)*f(a,i_left,j)
                u_right=u_right+real(ex(a),dp)*f(a,i_right,j)
            enddo
            u_mid=0.5_dp*(u_left+u_right)
            local_error=u_mid-exact
            write(*,'(a,i0,a,es24.16,a,es24.16,a,es24.16,a,es24.16)') &
                'MIDLINE j=',j,' y_over_h=',y/real(ny,dp), &
                ' u_num_over_umax=',u_mid/target_umax, &
                ' u_exact_over_umax=',exact/target_umax, &
                ' local_error_over_umax=',local_error/target_umax
        enddo
        write(*,'(a)') 'MIDLINE_PROFILE_END'
    end subroutine report_case

end program test_lbmcde_poiseuille_full2d
