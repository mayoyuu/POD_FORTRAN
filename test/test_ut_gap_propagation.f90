!> @file test_ut_gap_propagation.f90
!> @brief 读取 EMDAC 间隔前散点与后验协方差，用 UT 传播，与间隔后散点对比
program test_ut_gap_propagation
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_engine_module, only: pod_engine_init
    use pod_filter_ut_module, only: ut_filter
    use pod_spice, only: str2et
    implicit none

    character(len=MAX_STRING_LEN) :: scatter_before, scatter_after
    character(len=MAX_STRING_LEN) :: output_dir, config_file, arg_str
    integer :: i, n_lines
    real(DP), allocatable :: data(:, :), mean_vec(:), cov_mat(:, :)
    real(DP) :: et_before, et_after
    real(DP) :: noise_Q(6, 6), sigma_a
    real(DP) :: ut_cov(6, 6)
    type(ut_filter) :: ut

    ! 默认路径
    config_file    = 'config/config.txt'
    scatter_before = 'output/L1Halo-2_emdac_n5_before_gap.txt'
    scatter_after  = 'output/L1Halo-2_emdac_n5_after_gap.txt'
    output_dir     = 'output/L1Halo-2_emdac_n5_ut'

    ! 命令行参数
    i = 1
    do while (i <= command_argument_count())
        call get_command_argument(i, arg_str)
        select case (trim(arg_str))
            case ('-config'); call get_command_argument(i+1, arg_str); config_file    = trim(arg_str); i = i + 1
            case ('-before'); call get_command_argument(i+1, arg_str); scatter_before = trim(arg_str); i = i + 1
            case ('-after');  call get_command_argument(i+1, arg_str); scatter_after  = trim(arg_str); i = i + 1
            case ('-out');    call get_command_argument(i+1, arg_str); output_dir     = trim(arg_str); i = i + 1
        end select
        i = i + 1
    end do

    write(*,*) '=================================================='
    write(*,*) '  UT 长间隔协方差传播测试'
    write(*,*) '=================================================='
    write(*,*) '间隔前散点 : ', trim(scatter_before)
    write(*,*) '间隔后散点 : ', trim(scatter_after)
    write(*,*) '输出目录   : ', trim(output_dir)

    ! 1. 初始化引擎
    call pod_engine_init(trim(config_file))

    ! ===== 第一部分：间隔前（后验）散点 → UT 传播 =====
    write(*,*) ''
    write(*,*) '--- 间隔前散点（测后）高斯性诊断 ---'
    call read_scatter_file(scatter_before, data, n_lines)
    write(*,*) '  读取 ', n_lines, ' 个粒子'
    allocate(mean_vec(6), cov_mat(6, 6))
    call compute_mean_cov(data, n_lines, mean_vec, cov_mat)
    call gaussianity_check(data, n_lines, mean_vec, cov_mat)

    call write_cov_txt(trim(output_dir) // '/cov_before_gap.txt', cov_mat)
    write(*,*) '  已保存 cov_before_gap.txt'

    call str2et('2027-01-02T23:24:00.000', et_before)
    call str2et('2027-01-13T12:44:00.384', et_after)

    write(*,*) '--- UT 传播 ', et_after - et_before, ' 秒 ---'
    call ut%filter_init(et_before, mean_vec, cov_mat)

    sigma_a = 1.0e-11_DP
    noise_Q = 0.0_DP
    do i = 1, 3
        noise_Q(i,i)     = ((et_after - et_before)**4 / 4.0_DP) * sigma_a**2
        noise_Q(i+3,i+3) = (et_after - et_before)**2 * sigma_a**2
    end do

    call ut%time_update(et_after, noise_Q)

    call ut%get_current_cov(ut_cov)
    call write_cov_txt(trim(output_dir) // '/cov_after_gap.txt', ut_cov)
    write(*,*) '  已保存 cov_after_gap.txt'

    deallocate(data)

    ! ===== 第二部分：间隔后散点对比 =====
    write(*,*) ''
    write(*,*) '--- 间隔后散点高斯性诊断 ---'
    call read_scatter_file(scatter_after, data, n_lines)
    write(*,*) '  读取 ', n_lines, ' 个粒子'
    call compute_mean_cov(data, n_lines, mean_vec, cov_mat)
    call gaussianity_check(data, n_lines, mean_vec, cov_mat)

    call compare_cov(ut_cov, cov_mat)

    deallocate(data, mean_vec, cov_mat)

    write(*,*) '=================================================='
    write(*,*) '  完成！'
    write(*,*) '=================================================='

contains

    subroutine read_scatter_file(filename, samples, n)
        character(len=*), intent(in) :: filename
        real(DP), allocatable, intent(out) :: samples(:, :)
        integer, intent(out) :: n
        integer :: u, ios, line_num, cnt, p
        real(DP) :: v(6)
        character(len=512) :: line_buf

        open(newunit=u, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) '[ERROR] 无法打开: ', trim(filename)
            error stop
        end if
        cnt = 0
        do
            read(u, '(A)', iostat=ios) line_buf
            if (ios /= 0) exit
            if (len_trim(line_buf) > 0) cnt = cnt + 1
        end do
        rewind(u)

        n = cnt
        allocate(samples(6, n))
        do p = 1, n
            read(u, *, iostat=ios) line_num, v(1), v(2), v(3), v(4), v(5), v(6)
            if (ios /= 0) then
                write(*,*) '[ERROR] 读取行 ', p
                error stop
            end if
            samples(:, p) = v
        end do
        close(u)
    end subroutine read_scatter_file

    subroutine compute_mean_cov(samples, n, mean_vec, cov_mat)
        real(DP), intent(in) :: samples(:, :)
        integer, intent(in) :: n
        real(DP), intent(out) :: mean_vec(6), cov_mat(6, 6)
        integer :: i, j

        mean_vec = 0.0_DP
        do i = 1, n
            mean_vec = mean_vec + samples(:, i)
        end do
        mean_vec = mean_vec / real(n, DP)

        cov_mat = 0.0_DP
        do i = 1, n
            do j = 1, 6
                cov_mat(:, j) = cov_mat(:, j) &
                    + (samples(:, i) - mean_vec) * (samples(j, i) - mean_vec(j))
            end do
        end do
        cov_mat = cov_mat / real(n - 1, DP)
    end subroutine compute_mean_cov

    subroutine gaussianity_check(samples, n, mean_vec, cov_mat)
        real(DP), intent(in) :: samples(:, :), mean_vec(6), cov_mat(6, 6)
        integer, intent(in) :: n
        real(DP) :: diff(6), skew(6), ekurt(6), m2(6), m3(6), m4(6)
        real(DP) :: L(6, 6), y(6), md, md_mean, md_var, md_min, md_max
        integer :: i, j, k
        character(len=6), parameter :: labels(6) = ['x     ', 'y     ', 'z     ', &
                                                      'vx    ', 'vy    ', 'vz    ']

        write(*,*) '  维度    偏度        超额峰度'

        skew = 0.0_DP; ekurt = 0.0_DP; m2 = 0.0_DP; m3 = 0.0_DP; m4 = 0.0_DP
        do i = 1, n
            diff = samples(:, i) - mean_vec
            do j = 1, 6
                m2(j) = m2(j) + diff(j)**2
                m3(j) = m3(j) + diff(j)**3
                m4(j) = m4(j) + diff(j)**4
            end do
        end do
        do j = 1, 6
            skew(j)  = (m3(j) / n) / (m2(j) / n)**1.5_DP
            ekurt(j) = (m4(j) / n) / (m2(j) / n)**2 - 3.0_DP
        end do
        do j = 1, 6
            write(*, '(A, A, F10.4, F12.4)') '  ', labels(j), skew(j), ekurt(j)
        end do

        ! Cholesky 分解
        L = 0.0_DP
        do j = 1, 6
            do i = j, 6
                L(i, j) = cov_mat(i, j)
                do k = 1, j - 1
                    L(i, j) = L(i, j) - L(i, k) * L(j, k)
                end do
                if (i == j) then
                    if (L(j, j) <= 0.0_DP) then
                        write(*,*) '[WARN] 协方差非正定，跳过马氏距离'
                        return
                    end if
                    L(j, j) = sqrt(L(j, j))
                else
                    L(i, j) = L(i, j) / L(j, j)
                end if
            end do
        end do

        md_mean = 0.0_DP; md_var = 0.0_DP
        md_min = huge(1.0_DP); md_max = 0.0_DP
        do i = 1, n
            diff = samples(:, i) - mean_vec
            do j = 1, 6
                y(j) = diff(j)
                do k = 1, j - 1
                    y(j) = y(j) - L(j, k) * y(k)
                end do
                y(j) = y(j) / L(j, j)
            end do
            md = dot_product(y, y)
            md_mean = md_mean + md
            md_var  = md_var  + md * md
            if (md < md_min) md_min = md
            if (md > md_max) md_max = md
        end do
        md_mean = md_mean / n
        md_var  = md_var / n - md_mean**2

        write(*,*) '  马氏距离 (应 ~ chi2(6): mean=6, var=12):'
        write(*, '(A, F8.3)') '    mean  = ', md_mean
        write(*, '(A, F8.3)') '    std   = ', sqrt(md_var)
        write(*, '(A, F8.3)') '    min   = ', md_min
        write(*, '(A, F8.3)') '    max   = ', md_max
        write(*,*) ''
    end subroutine gaussianity_check

    subroutine compare_cov(ut_cov, emp_cov)
        real(DP), intent(in) :: ut_cov(6, 6), emp_cov(6, 6)
        real(DP) :: ut_std(6), emp_std(6), ratio(6), frob_diff, frob_ut, ut_trace, emp_trace
        integer :: j
        character(len=6), parameter :: labels(6) = ['x     ', 'y     ', 'z     ', &
                                                      'vx    ', 'vy    ', 'vz    ']

        write(*,*) '--- UT 传播协方差 vs 散点经验协方差 ---'
        write(*,*) '  维度   UT std        经验 std      比值(UT/经验)'
        do j = 1, 6
            ut_std(j)  = sqrt(ut_cov(j, j))
            emp_std(j) = sqrt(emp_cov(j, j))
            ratio(j) = ut_std(j) / emp_std(j)
            write(*, '(A, A, 2ES14.5, F10.4)') '  ', labels(j), ut_std(j), emp_std(j), ratio(j)
        end do

        frob_diff = sqrt(sum((ut_cov - emp_cov)**2))
        frob_ut   = sqrt(sum(ut_cov**2))
        ut_trace  = ut_cov(1,1)+ut_cov(2,2)+ut_cov(3,3)+ut_cov(4,4)+ut_cov(5,5)+ut_cov(6,6)
        emp_trace = emp_cov(1,1)+emp_cov(2,2)+emp_cov(3,3)+emp_cov(4,4)+emp_cov(5,5)+emp_cov(6,6)
        write(*,*)
        write(*, '(A, ES10.3)') '  协方差之差 Frobenius 范数 : ', frob_diff
        write(*, '(A, F8.4)')   '  相对差异 (diff / ||UT||): ', frob_diff / frob_ut
        write(*, '(A, ES14.5)') '  UT 协方差迹               : ', ut_trace
        write(*, '(A, ES14.5)') '  经验协方差迹               : ', emp_trace
        write(*,*) ''
    end subroutine compare_cov

    subroutine write_cov_txt(filename, cov_mat)
        character(len=*), intent(in) :: filename
        real(DP), intent(in) :: cov_mat(6, 6)
        integer :: u, ios, j

        open(newunit=u, file=trim(filename), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(*,*) '[ERROR] 无法创建: ', trim(filename)
            return
        end if
        do j = 1, 6
            write(u, '(6(ES25.16E3,A))') cov_mat(j, 1), ' ', cov_mat(j, 2), ' ', &
                cov_mat(j, 3), ' ', cov_mat(j, 4), ' ', cov_mat(j, 5), ' ', cov_mat(j, 6)
        end do
        close(u)
    end subroutine write_cov_txt

end program test_ut_gap_propagation
