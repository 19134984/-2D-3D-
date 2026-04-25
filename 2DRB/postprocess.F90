program postprocess
  use, intrinsic :: iso_fortran_env, only: error_unit, iostat_end
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)

  character(len=*), parameter :: default_nu_file = 'Nu_VolAvg_2DOpenacc.dat'
  character(len=*), parameter :: default_re_file = 'Re_VolAvg_2DOpenacc.dat'
  character(len=*), parameter :: default_raw_file = 'NuRe_VolAvg_2DOpenacc.plt'
  character(len=*), parameter :: default_mean_file = 'NuRe_VolAvg_runningMean_2DOpenacc.plt'

  character(len=512) :: nu_file, re_file, raw_file, mean_file
  real(kind=dp), allocatable :: time_nu(:), nu(:)
  real(kind=dp), allocatable :: time_re(:), re(:)
  real(kind=dp), allocatable :: nu_mean(:), re_mean(:)
  integer :: n_nu, n_re

  call parse_command_line(nu_file, re_file, raw_file, mean_file)

  call read_series(trim(nu_file), time_nu, nu, n_nu)
  call read_series(trim(re_file), time_re, re, n_re)
  call validate_series(time_nu, time_re, n_nu, n_re)

  allocate(nu_mean(n_nu), re_mean(n_nu))
  call compute_running_mean(nu, nu_mean, n_nu)
  call compute_running_mean(re, re_mean, n_nu)

  call write_raw_tecplot(trim(raw_file), time_nu, nu, re, n_nu)
  call write_mean_tecplot(trim(mean_file), time_nu, nu_mean, re_mean, n_nu)

  write(*,'(A)') 'Finished Nu/Re postprocessing.'
  write(*,'(A,1X,A)') 'Raw Tecplot file:', trim(raw_file)
  write(*,'(A,1X,A)') 'Running-mean Tecplot file:', trim(mean_file)

contains

  subroutine parse_command_line(nu_file, re_file, raw_file, mean_file)
    character(len=*), intent(out) :: nu_file, re_file, raw_file, mean_file
    integer :: nargs

    nargs = command_argument_count()

    if (nargs == 0) then
      nu_file = default_nu_file
      re_file = default_re_file
      raw_file = default_raw_file
      mean_file = default_mean_file
    else if (nargs == 4) then
      call get_command_argument(1, nu_file)
      call get_command_argument(2, re_file)
      call get_command_argument(3, raw_file)
      call get_command_argument(4, mean_file)
    else
      write(error_unit,'(A)') 'Usage: postprocess.exe [Nu.dat Re.dat raw.plt mean.plt]'
      write(error_unit,'(A)') 'Run with no arguments to use the default 2DRBOpenacc file names.'
      error stop 1
    end if
  end subroutine parse_command_line

  subroutine read_series(filename, time, value, n)
    character(len=*), intent(in) :: filename
    real(kind=dp), allocatable, intent(out) :: time(:), value(:)
    integer, intent(out) :: n

    integer :: u, ios, line_no
    real(kind=dp) :: t_tmp, v_tmp
    logical :: exists

    inquire(file=filename, exist=exists)
    if (.not. exists) then
      write(error_unit,'(A,1X,A)') 'Input file does not exist:', trim(filename)
      error stop 1
    end if

    n = 0
    line_no = 0
    open(newunit=u, file=filename, status='old', action='read', form='formatted')
    do
      read(u,*,iostat=ios) t_tmp, v_tmp
      if (ios == iostat_end) exit
      line_no = line_no + 1
      if (ios /= 0) then
        write(error_unit,'(A,1X,A,1X,A,I0)') 'Failed to read two numeric columns from', &
             trim(filename), 'at line', line_no
        close(u)
        error stop 1
      end if
      n = n + 1
    end do
    close(u)

    if (n <= 0) then
      write(error_unit,'(A,1X,A)') 'Input file is empty:', trim(filename)
      error stop 1
    end if

    allocate(time(n), value(n))

    line_no = 0
    open(newunit=u, file=filename, status='old', action='read', form='formatted')
    do line_no = 1, n
      read(u,*,iostat=ios) time(line_no), value(line_no)
      if (ios /= 0) then
        write(error_unit,'(A,1X,A,1X,A,I0)') 'Failed while rereading', trim(filename), &
             'at line', line_no
        close(u)
        error stop 1
      end if
    end do
    close(u)
  end subroutine read_series

  subroutine validate_series(time_nu, time_re, n_nu, n_re)
    real(kind=dp), intent(in) :: time_nu(:), time_re(:)
    integer, intent(in) :: n_nu, n_re

    integer :: k
    real(kind=dp) :: tol

    if (n_nu /= n_re) then
      write(error_unit,'(A,I0,A,I0)') 'Nu/Re sample counts differ: Nu=', n_nu, ', Re=', n_re
      error stop 1
    end if

    do k = 1, n_nu
      tol = 1.0d-10 * max(1.0d0, abs(time_nu(k)), abs(time_re(k)))
      if (abs(time_nu(k) - time_re(k)) > tol) then
        write(error_unit,'(A,I0)') 'Nu/Re time values differ at sample', k
        write(error_unit,'(A,1X,ES24.16)') 'Nu time:', time_nu(k)
        write(error_unit,'(A,1X,ES24.16)') 'Re time:', time_re(k)
        error stop 1
      end if
    end do
  end subroutine validate_series

  subroutine compute_running_mean(value, value_mean, n)
    real(kind=dp), intent(in) :: value(:)
    real(kind=dp), intent(out) :: value_mean(:)
    integer, intent(in) :: n

    integer :: k
    real(kind=dp) :: accum

    accum = 0.0d0
    do k = 1, n
      accum = accum + value(k)
      value_mean(k) = accum / dble(k)
    end do
  end subroutine compute_running_mean

  subroutine write_raw_tecplot(filename, time, nu, re, n)
    character(len=*), intent(in) :: filename
    real(kind=dp), intent(in) :: time(:), nu(:), re(:)
    integer, intent(in) :: n

    integer :: u, k

    open(newunit=u, file=filename, status='replace', action='write', form='formatted')
    write(u,'(A)') 'TITLE = "2D OpenACC Nu/Re volume averages"'
    write(u,'(A)') 'VARIABLES = "time" "NuVolAvg" "ReVolAvg"'
    write(u,'(A)') 'ZONE T="NuReVolAvg", F=POINT'
    do k = 1, n
      write(u,'(ES24.16,1X,ES24.16,1X,ES24.16)') time(k), nu(k), re(k)
    end do
    close(u)
  end subroutine write_raw_tecplot

  subroutine write_mean_tecplot(filename, time, nu_mean, re_mean, n)
    character(len=*), intent(in) :: filename
    real(kind=dp), intent(in) :: time(:), nu_mean(:), re_mean(:)
    integer, intent(in) :: n

    integer :: u, k

    open(newunit=u, file=filename, status='replace', action='write', form='formatted')
    write(u,'(A)') 'TITLE = "2D OpenACC Nu/Re running means"'
    write(u,'(A)') 'VARIABLES = "time" "NuVolAvgMean" "ReVolAvgMean"'
    write(u,'(A)') 'ZONE T="NuReRunningMean", F=POINT'
    do k = 1, n
      write(u,'(ES24.16,1X,ES24.16,1X,ES24.16)') time(k), nu_mean(k), re_mean(k)
    end do
    close(u)
  end subroutine write_mean_tecplot

end program postprocess
