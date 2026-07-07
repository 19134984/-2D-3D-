!===================================================================================================
! 程序: postprocess_statAvg_2DOpenacc
! 作用: 读取 2DRBOpenacc 输出的 Nu/Re 时间序列，并计算指定时间段内的统计收敛量。
!===================================================================================================
program postprocess_statAvg_2DOpenacc
  implicit none

  ! 这个程序和 postprocess.F90 分工不同：
  ! postprocess.F90 用来看时间序列和累计平均曲线，判断什么时候达到统计稳态；
  ! 本程序在给定的统计稳态起点之后，计算到 1000 t_ff 为止的整体/前半/后半平均值。

  ! 输入文件名写在源码中，对应 2DRBOpenacc.F90 中 calNuRe() 写出的两个 .dat 文件。
  character(len=*), parameter :: nu_file = 'Nu_VolAvg.dat'
  character(len=*), parameter :: re_file = 'Re_VolAvg.dat'

  ! Nu/Re 文件中已知有 2000 行数据，因此直接按固定行数读取。
  integer(kind=4), parameter :: sample_count = 2000

  ! 统计平均时间段。average_start_tf 根据 postprocess.F90 看到的统计稳态起点手动修改。
  real(kind=8), parameter :: average_start_tf = 500.0d0
  real(kind=8), parameter :: average_end_tf = 1000.0d0
  real(kind=8), parameter :: average_mid_tf = 0.5d0 * (average_start_tf + average_end_tf)

  ! 输出文件名
  character(len=*), parameter :: avg_file = 'NuRe_TimeAverage_2DOpenacc.txt'

  real(kind=8), allocatable :: time_nu(:), nu(:)
  real(kind=8), allocatable :: time_re(:), re(:)

  integer(kind=4) :: n
  integer(kind=4) :: whole_sample_count, first_sample_count, second_sample_count
  real(kind=8) :: nu_whole_average, re_whole_average
  real(kind=8) :: nu_first_average, re_first_average
  real(kind=8) :: nu_second_average, re_second_average
  real(kind=8) :: nu_first_relerr, re_first_relerr
  real(kind=8) :: nu_second_relerr, re_second_relerr

  n = sample_count
  allocate(time_nu(n), nu(n))
  allocate(time_re(n), re(n))

  call read_series(trim(nu_file), time_nu, nu, n)
  call read_series(trim(re_file), time_re, re, n)

  call validate_time_series(time_nu, time_re, n)

  call compute_window_averages(time_nu, nu, re, n, average_start_tf, average_mid_tf, average_end_tf, &
       nu_whole_average, re_whole_average, whole_sample_count, &
       nu_first_average, re_first_average, first_sample_count, &
       nu_second_average, re_second_average, second_sample_count, &
       nu_first_relerr, re_first_relerr, nu_second_relerr, re_second_relerr)

  call write_average_result(trim(avg_file), average_start_tf, average_mid_tf, average_end_tf, &
       whole_sample_count, first_sample_count, second_sample_count, &
       nu_whole_average, re_whole_average, nu_first_average, re_first_average, &
       nu_second_average, re_second_average, nu_first_relerr, re_first_relerr, &
       nu_second_relerr, re_second_relerr)

  write(*,'(A)') 'Finished Nu/Re interval time averaging.'
  write(*,'(A,1X,ES16.8)') 'Average start t_ff:', average_start_tf
  write(*,'(A,1X,ES16.8)') 'Average middle t_ff:', average_mid_tf
  write(*,'(A,1X,ES16.8)') 'Average end t_ff:', average_end_tf
  write(*,'(A,1X,I0)') 'Whole samples:', whole_sample_count
  write(*,'(A,1X,I0)') 'First-half samples:', first_sample_count
  write(*,'(A,1X,I0)') 'Second-half samples:', second_sample_count
  write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'Whole Nu/Re average:', &
       nu_whole_average, '/', re_whole_average
  write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'First-half Nu/Re average:', &
       nu_first_average, '/', re_first_average
  write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'Second-half Nu/Re average:', &
       nu_second_average, '/', re_second_average
  write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'First-half Nu/Re relative error:', &
       nu_first_relerr, '/', re_first_relerr
  write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') 'Second-half Nu/Re relative error:', &
       nu_second_relerr, '/', re_second_relerr
  write(*,'(A,1X,A)') 'Average file:', trim(avg_file)

end program postprocess_statAvg_2DOpenacc
!===================================================================================================
! postprocess_statAvg_2DOpenacc 结束。
!===================================================================================================

!===================================================================================================
! 子程序: read_series
! 作用: 读取一个两列时间序列文件，第一列为时间，第二列为统计量。
!===================================================================================================
subroutine read_series(filename, time, value, n)
  implicit none

  character(len=*), intent(in) :: filename
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(out) :: time(n), value(n)

  integer(kind=4) :: line_no
  logical :: exists

  inquire(file=filename, exist=exists)
  if (.not. exists) then
    write(*,'(A,1X,A)') 'Input file does not exist:', trim(filename)
    error stop 1
  end if

  open(unit=10, file=filename, status='old', action='read', form='formatted')
  do line_no = 1, n
    read(10,*) time(line_no), value(line_no)
  end do
  close(10)

end subroutine read_series
!===================================================================================================
! read_series 结束。
!===================================================================================================

!===================================================================================================
! 子程序: validate_time_series
! 作用: 检查 Nu/Re 两个时间序列是否能按行一一配对。
!===================================================================================================
subroutine validate_time_series(time_nu, time_re, n)
  implicit none

  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in) :: time_nu(n), time_re(n)

  integer(kind=4) :: k
  real(kind=8) :: tol

  do k = 1, n
    ! 用相对量级设置一个很小的容差，避免文本读写造成的最后几位差异误判。
    tol = 1.0d-10 * max(1.0d0, abs(time_nu(k)), abs(time_re(k)))
    if (abs(time_nu(k) - time_re(k)) > tol) then
      write(*,'(A,I0)') 'Nu/Re time values differ at sample', k
      write(*,'(A,1X,ES24.16)') 'Nu time:', time_nu(k)
      write(*,'(A,1X,ES24.16)') 'Re time:', time_re(k)
      error stop 1
    end if
  end do

end subroutine validate_time_series
!===================================================================================================
! validate_time_series 结束。
!===================================================================================================

!===================================================================================================
! 子程序: compute_window_averages
! 作用: 对整体/前半/后半窗口分别做时间平均，并计算半窗口相对整体的误差。
!===================================================================================================
subroutine compute_window_averages(time, nu, re, n, start_tf, mid_tf, end_tf, &
     nu_whole_average, re_whole_average, whole_sample_count, &
     nu_first_average, re_first_average, first_sample_count, &
     nu_second_average, re_second_average, second_sample_count, &
     nu_first_relerr, re_first_relerr, nu_second_relerr, re_second_relerr)
  implicit none

  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in) :: time(n), nu(n), re(n)
  real(kind=8), intent(in) :: start_tf, mid_tf, end_tf
  real(kind=8), intent(out) :: nu_whole_average, re_whole_average
  real(kind=8), intent(out) :: nu_first_average, re_first_average
  real(kind=8), intent(out) :: nu_second_average, re_second_average
  real(kind=8), intent(out) :: nu_first_relerr, re_first_relerr
  real(kind=8), intent(out) :: nu_second_relerr, re_second_relerr
  integer(kind=4), intent(out) :: whole_sample_count, first_sample_count, second_sample_count

  integer(kind=4) :: k
  real(kind=8) :: nu_whole_sum, re_whole_sum
  real(kind=8) :: nu_first_sum, re_first_sum
  real(kind=8) :: nu_second_sum, re_second_sum

  if ((start_tf >= mid_tf) .or. (mid_tf >= end_tf)) then
    write(*,'(A)') 'Average window times must satisfy start_tf < mid_tf < end_tf.'
    error stop 1
  end if

  nu_whole_sum = 0.0d0
  re_whole_sum = 0.0d0
  nu_first_sum = 0.0d0
  re_first_sum = 0.0d0
  nu_second_sum = 0.0d0
  re_second_sum = 0.0d0

  whole_sample_count = 0
  first_sample_count = 0
  second_sample_count = 0

  do k = 1, n
    if ((time(k) >= start_tf) .and. (time(k) <= end_tf)) then
      nu_whole_sum = nu_whole_sum + nu(k)
      re_whole_sum = re_whole_sum + re(k)
      whole_sample_count = whole_sample_count + 1
    end if

    if ((time(k) >= start_tf) .and. (time(k) < mid_tf)) then  
      nu_first_sum = nu_first_sum + nu(k)
      re_first_sum = re_first_sum + re(k)
      first_sample_count = first_sample_count + 1
    end if

    if ((time(k) >= mid_tf) .and. (time(k) <= end_tf)) then
      nu_second_sum = nu_second_sum + nu(k)
      re_second_sum = re_second_sum + re(k)
      second_sample_count = second_sample_count + 1
    end if
  end do

  if ((whole_sample_count <= 0) .or. (first_sample_count <= 0) .or. &
      (second_sample_count <= 0)) then
    write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8)') &
         'No samples found in interval:', start_tf, 'to', end_tf
    error stop 1
  end if

  if (first_sample_count + second_sample_count /= whole_sample_count) then
    write(*,'(A)') 'First-half and second-half sample counts do not match whole-window count.'
    error stop 1
  end if

  nu_whole_average = nu_whole_sum / dble(whole_sample_count)
  re_whole_average = re_whole_sum / dble(whole_sample_count)
  nu_first_average = nu_first_sum / dble(first_sample_count)
  re_first_average = re_first_sum / dble(first_sample_count)
  nu_second_average = nu_second_sum / dble(second_sample_count)
  re_second_average = re_second_sum / dble(second_sample_count)

  if ((abs(nu_whole_average) <= tiny(1.0d0)) .or. &   !tiny(1.0d0) 是 Fortran 里对应双精度实数的一个极小正数
      (abs(re_whole_average) <= tiny(1.0d0))) then
    write(*,'(A)') 'Whole-window average is zero; relative error cannot be computed.'
    error stop 1
  end if

  nu_first_relerr = abs(nu_first_average - nu_whole_average) / abs(nu_whole_average)
  re_first_relerr = abs(re_first_average - re_whole_average) / abs(re_whole_average)
  nu_second_relerr = abs(nu_second_average - nu_whole_average) / abs(nu_whole_average)
  re_second_relerr = abs(re_second_average - re_whole_average) / abs(re_whole_average)

end subroutine compute_window_averages
!===================================================================================================
! compute_window_averages 结束。
!===================================================================================================

!===================================================================================================
! 子程序: write_average_result
! 作用: 写出最终统计平均值。
!===================================================================================================
subroutine write_average_result(filename, start_tf, mid_tf, end_tf, &
     whole_sample_count, first_sample_count, second_sample_count, &
     nu_whole_average, re_whole_average, nu_first_average, re_first_average, &
     nu_second_average, re_second_average, nu_first_relerr, re_first_relerr, &
     nu_second_relerr, re_second_relerr)
  implicit none

  character(len=*), intent(in) :: filename
  real(kind=8), intent(in) :: start_tf, mid_tf, end_tf
  integer(kind=4), intent(in) :: whole_sample_count, first_sample_count, second_sample_count
  real(kind=8), intent(in) :: nu_whole_average, re_whole_average
  real(kind=8), intent(in) :: nu_first_average, re_first_average
  real(kind=8), intent(in) :: nu_second_average, re_second_average
  real(kind=8), intent(in) :: nu_first_relerr, re_first_relerr
  real(kind=8), intent(in) :: nu_second_relerr, re_second_relerr

  open(unit=20, file=filename, status='replace', action='write', form='formatted')
  write(20,'(A)') '# 2D OpenACC Nu/Re statistical-convergence window averages'
  write(20,'(A)') '# start_tf mid_tf end_tf whole_count first_count second_count ' // &
       'Nu_whole Re_whole Nu_first Re_first Nu_second Re_second ' // &
       'Nu_first_relerr Re_first_relerr Nu_second_relerr Re_second_relerr'
  write(20,'(ES24.16,1X,ES24.16,1X,ES24.16,1X,I0,1X,I0,1X,I0,1X,' // &
       'ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16,1X,' // &
       'ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16)') &
       start_tf, mid_tf, end_tf, whole_sample_count, first_sample_count, &
       second_sample_count, nu_whole_average, re_whole_average, nu_first_average, &
       re_first_average, nu_second_average, re_second_average, nu_first_relerr, &
       re_first_relerr, nu_second_relerr, re_second_relerr
  close(20)

end subroutine write_average_result
!===================================================================================================
! write_average_result 结束。
!===================================================================================================
