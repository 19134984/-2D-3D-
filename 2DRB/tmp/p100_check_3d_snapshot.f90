program p100_check_3d_snapshot
  implicit none
  integer, parameter :: nx=257, ny=257, nz=257
  real(8), parameter :: stretchA=1.5d0, Thot=0.5d0, Tcold=-0.5d0
  real(8), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:), T(:,:,:), rho(:,:,:)
  real(8) :: xp(1:nx), yp(1:ny), zp(1:nz)
  real(8) :: raw, erfNorm, c0, c1, c2, dTdy, NuVal
  real(8) :: urms, vrms, wrms, vrmsSlice, wrmsSlice, nuMin, nuMax, nuSum
  integer :: i, j, k, p, idx(9)

  allocate(u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), T(nx,ny,nz), rho(nx,ny,nz))

  open(10, file='buoyancyCavity3DOpenaccISLBMSnapshot-000001338000.bin', &
       form='unformatted', access='sequential', status='old')
  read(10) (((u(i,j,k), i=1,nx), j=1,ny), k=1,nz)
  read(10) (((v(i,j,k), i=1,nx), j=1,ny), k=1,nz)
  read(10) (((w(i,j,k), i=1,nx), j=1,ny), k=1,nz)
  read(10) (((T(i,j,k), i=1,nx), j=1,ny), k=1,nz)
  read(10) (((rho(i,j,k), i=1,nx), j=1,ny), k=1,nz)
  close(10)

  erfNorm = erf(0.5d0*stretchA)
  do i = 1, nx
    raw = 0.5d0*(1.0d0 + erf(stretchA*(dble(i)/dble(nx+1)-0.5d0))/erfNorm)
    xp(i) = (raw - 0.5d0*first_raw(nx, stretchA, erfNorm))/(1.0d0-first_raw(nx, stretchA, erfNorm))
  enddo
  do j = 1, ny
    raw = 0.5d0*(1.0d0 + erf(stretchA*(dble(j)/dble(ny+1)-0.5d0))/erfNorm)
    yp(j) = (raw - 0.5d0*first_raw(ny, stretchA, erfNorm))/(1.0d0-first_raw(ny, stretchA, erfNorm))
  enddo
  do k = 1, nz
    raw = 0.5d0*(1.0d0 + erf(stretchA*(dble(k)/dble(nz+1)-0.5d0))/erfNorm)
    zp(k) = (raw - 0.5d0*first_raw(nz, stretchA, erfNorm))/(1.0d0-first_raw(nz, stretchA, erfNorm))
  enddo

  urms = sqrt(sum(u*u)/dble(nx*ny*nz))
  vrms = sqrt(sum(v*v)/dble(nx*ny*nz))
  wrms = sqrt(sum(w*w)/dble(nx*ny*nz))
  write(*,'(a,3(1x,es16.8))') 'global_rms_u_v_w', urms, vrms, wrms
  write(*,'(a,4(1x,es16.8))') 'T_minmax_rho_minmax', minval(T), maxval(T), minval(rho), maxval(rho)

  c0 = (0.0d0-yp(1)-yp(2))/((0.0d0-yp(1))*(0.0d0-yp(2)))
  c1 = (0.0d0-yp(2))/((yp(1)-0.0d0)*(yp(1)-yp(2)))
  c2 = (0.0d0-yp(1))/((yp(2)-0.0d0)*(yp(2)-yp(1)))

  idx = (/1,2,3,nx/4,nx/2,3*nx/4,nx-2,nx-1,nx/)
  nuMin = huge(1.0d0)
  nuMax = -huge(1.0d0)
  nuSum = 0.0d0
  do i = 1, nx
    NuVal = 0.0d0
    do k = 1, nz
      dTdy = c0*Thot + c1*T(i,1,k) + c2*T(i,2,k)
      NuVal = NuVal - dTdy
    enddo
    NuVal = NuVal/dble(nz)
    nuMin = min(nuMin, NuVal)
    nuMax = max(nuMax, NuVal)
    nuSum = nuSum + NuVal
  enddo
  write(*,'(a,3(1x,es16.8))') 'Nu_hot_x_simple_mean_min_max', nuSum/dble(nx), nuMin, nuMax
  do p = 1, size(idx)
    i = idx(p)
    NuVal = 0.0d0
    do k = 1, nz
      dTdy = c0*Thot + c1*T(i,1,k) + c2*T(i,2,k)
      NuVal = NuVal - dTdy
    enddo
    NuVal = NuVal/dble(nz)
    vrmsSlice = sqrt(sum(v(i,:,:)*v(i,:,:))/dble(ny*nz))
    wrmsSlice = sqrt(sum(w(i,:,:)*w(i,:,:))/dble(ny*nz))
    write(*,'(a,i5,4(1x,es16.8))') 'x_sample_i_x_Nu_vrms_wrms', i, xp(i), NuVal, vrmsSlice, wrmsSlice
  enddo

contains
  real(8) function first_raw(n, a, norm)
    integer, intent(in) :: n
    real(8), intent(in) :: a, norm
    first_raw = 0.5d0*(1.0d0 + erf(a*(1.0d0/dble(n+1)-0.5d0))/norm)
  end function first_raw
end program p100_check_3d_snapshot
