!===================================================================================================
! 程序: postprocess
! 作用: 读取 2DRBOpenacc 输出的 Nu/Re 体平均时间序列，并生成 Tecplot 后处理文件。
!===================================================================================================
program postprocess
  implicit none

  ! 后处理用途:
  ! 1) 读取 2DRBOpenacc.F90 输出的 NuVolAvg / ReVolAvg 时间序列。
  ! 2) 输出一个原始序列 Tecplot 文件。
  ! 3) 输出一个累计时间平均 Tecplot 文件，其中第 k 个点为 mean(1:k)。

  ! 输入文件名写在源码中，对应 2DRBOpenacc.F90 中 calNuRe() 写出的两个 .dat 文件。
  character(len=*), parameter :: nu_file = 'Nu_VolAvg.dat'
  character(len=*), parameter :: re_file = 'Re_VolAvg.dat'

  ! Nu/Re 文件中已知有 2000 行数据，因此直接按固定行数读取。
  integer(kind=4), parameter :: sample_count = 2000

  ! 输出文件名
  character(len=*), parameter :: series_file = 'NuRe_VolAvg_2DOpenacc.plt'
  character(len=*), parameter :: mean_file = 'NuRe_VolAvg_runningMean_2DOpenacc.plt'

  ! time_nu/time_re 分别来自 Nu/Re 文件的第一列；后续会检查两者是否一致。
  real(kind=8), allocatable :: time_nu(:), nu(:)
  real(kind=8), allocatable :: time_re(:), re(:)

  ! 累计平均序列：nu_mean(k)=sum(nu(1:k))/k，re_mean 同理。
  real(kind=8), allocatable :: nu_mean(:), re_mean(:)

  integer(kind=4) :: n

  n = sample_count
  allocate(time_nu(n), nu(n))
  allocate(time_re(n), re(n))
  allocate(nu_mean(n), re_mean(n))

  ! 分别读取 Nu 和 Re 的两列数据：第一列为时间，第二列为体平均量。
  call read_series(trim(nu_file), time_nu, nu, n)
  call read_series(trim(re_file), time_re, re, n)

  ! Nu/Re 必须按同一采样时刻一一对应，否则合并到同一个 Tecplot 文件会失去意义。
  call validate_series(time_nu, time_re, n)

  ! 包含当前点的累计时间平均：第 5 个点平均前 5 个样本，第 8 个点平均前 8 个样本。
  call compute_running_mean(nu, nu_mean, n)
  call compute_running_mean(re, re_mean, n)

  ! 输出两个独立的 ASCII Tecplot 文件：原始序列和累计平均序列。
  call write_series_tecplot(trim(series_file), time_nu, nu, re, n)
  call write_mean_tecplot(trim(mean_file), time_nu, nu_mean, re_mean, n)

  write(*,'(A)') 'Finished Nu/Re postprocessing.'
  write(*,'(A,1X,A)') 'Series Tecplot file:', trim(series_file)
  write(*,'(A,1X,A)') 'Running-mean Tecplot file:', trim(mean_file)

end program postprocess
!===================================================================================================
! postprocess 结束: Nu/Re 后处理流程完成。
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
    ! 文件缺失属于错误退出；error stop 后面的 1 是非零退出码，便于脚本判断运行失败。
    error stop 1
  end if

  open(unit=10, file=filename, status='old', action='read', form='formatted')
  do line_no = 1, n
    read(10,*) time(line_no), value(line_no)
  end do
  close(10)

end subroutine read_series
!===================================================================================================
! read_series 结束: 两列时间序列读取完成。
!===================================================================================================

!===================================================================================================
! 子程序: validate_series
! 作用: 检查 Nu/Re 两个时间序列是否能按行一一配对。
!===================================================================================================
subroutine validate_series(time_nu, time_re, n)
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

end subroutine validate_series
!===================================================================================================
! validate_series 结束: Nu/Re 时间列一致性检查完成。
!===================================================================================================

!===================================================================================================
! 子程序: compute_running_mean
! 作用: 计算包含当前样本的累计时间平均序列。
!===================================================================================================
subroutine compute_running_mean(value, value_mean, n)
  implicit none

  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in) :: value(n)
  real(kind=8), intent(out) :: value_mean(n)

  integer(kind=4) :: k
  real(kind=8) :: accum

  accum = 0.0d0
  do k = 1, n
    accum = accum + value(k)
    value_mean(k) = accum / dble(k)
  end do

end subroutine compute_running_mean
!===================================================================================================
! compute_running_mean 结束: 累计时间平均计算完成。
!===================================================================================================

!===================================================================================================
! 子程序: write_series_tecplot
! 作用: 写出原始 NuVolAvg/ReVolAvg 时间序列到 ASCII Tecplot 文件。
!===================================================================================================
subroutine write_series_tecplot(filename, time, nu, re, n)
  implicit none

  character(len=*), intent(in) :: filename
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in) :: time(n), nu(n), re(n)

  integer(kind=4) :: k

  open(unit=20, file=filename, status='replace', action='write', form='formatted')
  ! ASCII Tecplot 的基本头部；F=POINT 表示每一行是一组变量值。
  write(20,'(A)') 'TITLE = "2D OpenACC Nu/Re volume averages"'
  write(20,'(A)') 'VARIABLES = "time" "NuVolAvg" "ReVolAvg"'
  write(20,'(A)') 'ZONE T="NuReVolAvg", F=POINT'
  do k = 1, n
    write(20,'(ES24.16,1X,ES24.16,1X,ES24.16)') time(k), nu(k), re(k)
  end do
  close(20)

end subroutine write_series_tecplot
!===================================================================================================
! write_series_tecplot 结束: 原始 Nu/Re Tecplot 文件写出完成。
!===================================================================================================

!===================================================================================================
! 子程序: write_mean_tecplot
! 作用: 写出 NuVolAvg/ReVolAvg 的累计时间平均序列到 ASCII Tecplot 文件。
!===================================================================================================
subroutine write_mean_tecplot(filename, time, nu_mean, re_mean, n)
  implicit none

  character(len=*), intent(in) :: filename
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in) :: time(n), nu_mean(n), re_mean(n)

  integer(kind=4) :: k

  open(unit=21, file=filename, status='replace', action='write', form='formatted')
  ! 这里的时间列仍然使用原始采样时间，第二/三列是从第一个样本累计到当前样本的均值。
  write(21,'(A)') 'TITLE = "2D OpenACC Nu/Re running means"'
  write(21,'(A)') 'VARIABLES = "time" "NuVolAvgMean" "ReVolAvgMean"'
  write(21,'(A)') 'ZONE T="NuReRunningMean", F=POINT'
  do k = 1, n
    write(21,'(ES24.16,1X,ES24.16,1X,ES24.16)') time(k), nu_mean(k), re_mean(k)
  end do
  close(21)

end subroutine write_mean_tecplot
!===================================================================================================
! write_mean_tecplot 结束: 累计平均 Nu/Re Tecplot 文件写出完成。
!===================================================================================================
