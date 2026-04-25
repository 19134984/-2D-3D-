program postprocess
  implicit none

  ! 后处理用途:
  ! 1) 读取 2DRBOpenacc.F90 输出的 NuVolAvg / ReVolAvg 时间序列。
  ! 2) 输出一个原始序列 Tecplot 文件。
  ! 3) 输出一个累计时间平均 Tecplot 文件，其中第 k 个点为 mean(1:k)。

  ! 默认输入文件名，对应 2DRBOpenacc.F90 中 calNuRe() 追加写出的两个 .dat 文件。
  character(len=*), parameter :: default_nu_file = 'Nu_VolAvg_2DOpenacc.dat'
  character(len=*), parameter :: default_re_file = 'Re_VolAvg_2DOpenacc.dat'

  ! 默认输出文件名，使用 ASCII Tecplot 格式，便于直接查看和手工检查。
  character(len=*), parameter :: default_raw_file = 'NuRe_VolAvg_2DOpenacc.plt'
  character(len=*), parameter :: default_mean_file = 'NuRe_VolAvg_runningMean_2DOpenacc.plt'

  character(len=512) :: nu_file, re_file, raw_file, mean_file

  ! time_nu/time_re 分别来自 Nu/Re 文件的第一列；后续会检查两者是否一致。
  real(kind=8), allocatable :: time_nu(:), nu(:)
  real(kind=8), allocatable :: time_re(:), re(:)

  ! 累计平均序列：nu_mean(k)=sum(nu(1:k))/k，re_mean 同理。
  real(kind=8), allocatable :: nu_mean(:), re_mean(:)
  integer(kind=4) :: n_nu, n_re

  ! 若无命令行参数，则使用默认文件名；若有参数，必须完整给出 4 个文件名。
  call parse_command_line(nu_file, re_file, raw_file, mean_file)

  ! 分别读取 Nu 和 Re 的两列数据：第一列为时间，第二列为体平均量。
  call read_series(trim(nu_file), time_nu, nu, n_nu)
  call read_series(trim(re_file), time_re, re, n_re)

  ! Nu/Re 必须按同一采样时刻一一对应，否则合并到同一个 Tecplot 文件会失去意义。
  call validate_series(time_nu, time_re, n_nu, n_re)

  allocate(nu_mean(n_nu), re_mean(n_nu))

  ! 包含当前点的累计时间平均：第 5 个点平均前 5 个样本，第 8 个点平均前 8 个样本。
  call compute_running_mean(nu, nu_mean, n_nu)
  call compute_running_mean(re, re_mean, n_nu)

  ! 输出两个独立的 ASCII Tecplot 文件：原始序列和累计平均序列。
  call write_raw_tecplot(trim(raw_file), time_nu, nu, re, n_nu)
  call write_mean_tecplot(trim(mean_file), time_nu, nu_mean, re_mean, n_nu)

  write(*,'(A)') 'Finished Nu/Re postprocessing.'
  write(*,'(A,1X,A)') 'Raw Tecplot file:', trim(raw_file)
  write(*,'(A,1X,A)') 'Running-mean Tecplot file:', trim(mean_file)

contains

  ! 读取命令行参数。
  ! 用法:
  !   postprocess.exe
  !   postprocess.exe Nu.dat Re.dat raw.plt mean.plt
  subroutine parse_command_line(nu_file, re_file, raw_file, mean_file)
    character(len=*), intent(out) :: nu_file, re_file, raw_file, mean_file
    integer(kind=4) :: nargs

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
      write(*,'(A)') 'Usage: postprocess.exe [Nu.dat Re.dat raw.plt mean.plt]'
      write(*,'(A)') 'Run with no arguments to use the default 2DRBOpenacc file names.'
      error stop 1
    end if
  end subroutine parse_command_line

  ! 读取一个两列时间序列文件。
  ! 第一遍只统计有效数据行数 n，这样数组可以按真实样本数分配；
  ! 第二遍再把 time/value 读入数组。
  subroutine read_series(filename, time, value, n)
    character(len=*), intent(in) :: filename
    real(kind=8), allocatable, intent(out) :: time(:), value(:)
    integer(kind=4), intent(out) :: n

    integer(kind=4) :: ios, line_no
    real(kind=8) :: t_tmp, v_tmp
    logical :: exists

    inquire(file=filename, exist=exists)
    if (.not. exists) then
      write(*,'(A,1X,A)') 'Input file does not exist:', trim(filename)
      error stop 1
    end if

    ! 第一遍读取：检查每一行是否是两个数，并统计样本数。
    n = 0
    line_no = 0
    open(unit=10, file=filename, status='old', action='read', form='formatted')
    do
      read(10,*,iostat=ios) t_tmp, v_tmp
      if (ios < 0) exit
      line_no = line_no + 1
      if (ios /= 0) then
        write(*,'(A,1X,A,1X,A,I0)') 'Failed to read two numeric columns from', &
             trim(filename), 'at line', line_no
        close(10)
        error stop 1
      end if
      n = n + 1
    end do
    close(10)

    if (n <= 0) then
      write(*,'(A,1X,A)') 'Input file is empty:', trim(filename)
      error stop 1
    end if

    allocate(time(n), value(n))

    ! 第二遍读取：把数据真正装入数组。
    line_no = 0
    open(unit=10, file=filename, status='old', action='read', form='formatted')
    do line_no = 1, n
      read(10,*,iostat=ios) time(line_no), value(line_no)
      if (ios /= 0) then
        write(*,'(A,1X,A,1X,A,I0)') 'Failed while rereading', trim(filename), &
             'at line', line_no
        close(10)
        error stop 1
      end if
    end do
    close(10)
  end subroutine read_series

  ! 检查 Nu 和 Re 两个时间序列是否能逐行配对。
  ! 如果样本数不同，或者同一行的时间值不同，就停止并报错。
  subroutine validate_series(time_nu, time_re, n_nu, n_re)
    real(kind=8), intent(in) :: time_nu(:), time_re(:)
    integer(kind=4), intent(in) :: n_nu, n_re

    integer(kind=4) :: k
    real(kind=8) :: tol

    if (n_nu /= n_re) then
      write(*,'(A,I0,A,I0)') 'Nu/Re sample counts differ: Nu=', n_nu, ', Re=', n_re
      error stop 1
    end if

    do k = 1, n_nu
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

  ! 计算包含当前样本的累计平均:
  ! value_mean(k) = (value(1)+...+value(k))/k。
  subroutine compute_running_mean(value, value_mean, n)
    real(kind=8), intent(in) :: value(:)
    real(kind=8), intent(out) :: value_mean(:)
    integer(kind=4), intent(in) :: n

    integer(kind=4) :: k
    real(kind=8) :: accum

    accum = 0.0d0
    do k = 1, n
      accum = accum + value(k)
      value_mean(k) = accum / dble(k)
    end do
  end subroutine compute_running_mean

  ! 写出原始 NuVolAvg/ReVolAvg 序列到 ASCII Tecplot 文件。
  subroutine write_raw_tecplot(filename, time, nu, re, n)
    character(len=*), intent(in) :: filename
    real(kind=8), intent(in) :: time(:), nu(:), re(:)
    integer(kind=4), intent(in) :: n

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
  end subroutine write_raw_tecplot

  ! 写出 NuVolAvg/ReVolAvg 的累计时间平均序列到 ASCII Tecplot 文件。
  subroutine write_mean_tecplot(filename, time, nu_mean, re_mean, n)
    character(len=*), intent(in) :: filename
    real(kind=8), intent(in) :: time(:), nu_mean(:), re_mean(:)
    integer(kind=4), intent(in) :: n

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

end program postprocess
