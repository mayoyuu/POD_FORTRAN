! ============================================================================
! 测试: sample_particles_from_gmm 采样准确性验证
!       从 n5.gmms.json 解析两个 STEP 的 GMM，采样 1e5 粒子，
!       对比经验统计量与理论值，输出散点 CSV。
! ============================================================================
program test_sample_particles_from_gmm
    use pod_global, only: DP
    use pod_filter_emdac_module, only: emdac_filter
    use pod_uq_gmm_state_module, only: uq_gmm_state_type, gaussian_component
    use pod_random_module, only: init_random_seed

    implicit none

    integer, parameter :: DIM = 6
    integer, parameter :: N_PARTICLES = 100000
    integer, parameter :: MAX_STEPS = 2
    character(len=*), parameter :: OUT_DIR = 'OPM/L1Halo-1/'
    character(len=*), parameter :: JSON_FILE = &
        'OPM/L1Halo-1/L1Halo-1_supp_single_R91_1h_emdac_n5.gmms.json'

    type(emdac_filter) :: filter
    type(uq_gmm_state_type) :: gmms(MAX_STEPS)
    character(len=40) :: step_labels(MAX_STEPS)
    integer :: s

    call init_random_seed(.true.)
    call parse_gmms_json(JSON_FILE, gmms, step_labels)

    do s = 1, MAX_STEPS
        write(*,*) '>>> 测试 STEP ', s, ': ', trim(step_labels(s))
        call test_gmm_step(step_labels(s), gmms(s))
    end do

    write(*,*) '============================================================'
    write(*,*) '  测试完成！输出目录: ', OUT_DIR
    write(*,*) '============================================================'

contains

    ! ====================================================================
    ! 解析 n5.gmms.json，提取两个 STEP 的 GMM 参数
    ! 文件结构:
    !   "STEPS": [
    !       { "INDEX": 1, "LABEL": "...", "N_COMPONENTS": N, "COMPONENTS": [
    !           { "INDEX":1, "WEIGHT":w, "MEAN":[...], "COV":[[...],...] },
    !           ...
    !       ]},
    !       { "INDEX": 2, ... }
    !   ]
    ! ====================================================================
    subroutine parse_gmms_json(filename, gmms_out, labels_out)
        character(len=*), intent(in) :: filename
        type(uq_gmm_state_type), intent(out) :: gmms_out(MAX_STEPS)
        character(len=40), intent(out) :: labels_out(MAX_STEPS)

        integer :: u, ios, step_idx, comp_idx, cov_row, i
        character(len=512) :: line
        real(DP) :: w, m(DIM), c(DIM, DIM)
        integer :: ncomp
        logical :: in_components, in_cov, skip_cov_close

        open(newunit=u, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) '[ERROR] 无法打开: ', trim(filename); stop
        end if

        step_idx = 0
        comp_idx = 0
        in_components = .false.
        in_cov = .false.
        skip_cov_close = .false.
        cov_row = 0

        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)

            ! ---- 检测 STEP 层级 ----
            if (index(line, '"LABEL"') > 0 .and. .not. in_components) then
                step_idx = step_idx + 1
                labels_out(step_idx) = extract_string_value(line)
                cycle
            end if

            if (index(line, '"N_COMPONENTS"') > 0 .and. .not. in_components) then
                ncomp = extract_int_value(line)
                call gmms_out(step_idx)%allocate_components(ncomp, DIM)
                cycle
            end if

            if (index(line, '"COMPONENTS"') > 0) then
                in_components = .true.
                comp_idx = 0
                cycle
            end if

            ! ---- 检测分量层级 (在 COMPONENTS 内部) ----
            if (in_components .and. index(line, '"INDEX"') > 0) then
                comp_idx = comp_idx + 1
                ! 重置当前分量的临时缓冲区
                w = 0.0_DP; m = 0.0_DP; c = 0.0_DP
                in_cov = .false.; cov_row = 0
                cycle
            end if

            if (in_components .and. comp_idx > 0 .and. .not. in_cov) then
                if (index(line, '"WEIGHT"') > 0) then
                    w = extract_double_value(line)
                    cycle
                end if
                if (index(line, '"MEAN"') > 0) then
                    call parse_array_values(line, m, DIM)
                    cycle
                end if
                if (index(line, '"COV"') > 0) then
                    in_cov = .true.
                    cov_row = 0
                    cycle
                end if
            end if

            ! ---- 读取 COV 矩阵的行 ----
            if (in_cov .and. index(line, '[') > 0) then
                cov_row = cov_row + 1
                call parse_array_values(line, c(cov_row, :), DIM)
                if (cov_row >= DIM) then
                    in_cov = .false.
                    skip_cov_close = .true.   ! 下一行 ']' 是 COV 数组的关闭，跳过
                    cov_row = 0
                    ! 分量完整: 初始化到 GMM
                    call gmms_out(step_idx)%components(comp_idx)%init(w, m, c)
                end if
                cycle
            end if

            ! ---- 退出 COMPONENTS 数组 ----
            if (in_components .and. .not. in_cov) then
                if (is_closing_bracket(line, ']')) then
                    if (skip_cov_close) then
                        ! 这是 COV 矩阵的关闭 ']'，忽略
                        skip_cov_close = .false.
                    else
                        ! 这是 COMPONENTS 数组的关闭 ']'
                        in_components = .false.
                        comp_idx = 0
                    end if
                    cycle
                end if
            end if
        end do

        close(u)
        write(*,*) '[解析] 成功加载 ', step_idx, ' 个 STEP from ', trim(filename)
    end subroutine parse_gmms_json

    ! ====================================================================
    ! 辅助函数: 提取 "KEY": value 中的字符串值
    ! ====================================================================
    function extract_string_value(line) result(val)
        character(len=*), intent(in) :: line
        character(len=:), allocatable :: val
        integer :: p1, p2, n
        p1 = index(line, ':')
        if (p1 > 0) then
            line_scan: do n = p1+1, len_trim(line)
                if (line(n:n) == '"') then
                    p1 = n + 1
                    do p2 = p1, len_trim(line)
                        if (line(p2:p2) == '"') exit line_scan
                    end do
                end if
            end do line_scan
            val = line(p1:p2-1)
        else
            val = ''
        end if
    end function extract_string_value

    ! ====================================================================
    ! 辅助函数: 提取 "KEY": INT_VALUE 中的整数值
    ! ====================================================================
    function extract_int_value(line) result(ival)
        character(len=*), intent(in) :: line
        integer :: ival, p1
        character(len=64) :: tmp
        p1 = index(line, ':')
        if (p1 > 0) then
            tmp = adjustl(line(p1+1:))
            call remove_trailing_comma(tmp)
            read(tmp, *) ival
        else
            ival = 0
        end if
    end function extract_int_value

    ! ====================================================================
    ! 辅助函数: 提取 "KEY": DOUBLE_VALUE 中的浮点值
    ! ====================================================================
    function extract_double_value(line) result(dval)
        character(len=*), intent(in) :: line
        real(DP) :: dval
        integer :: p1
        character(len=64) :: tmp
        p1 = index(line, ':')
        if (p1 > 0) then
            tmp = adjustl(line(p1+1:))
            call remove_trailing_comma(tmp)
            read(tmp, *) dval
        else
            dval = 0.0_DP
        end if
    end function extract_double_value

    ! ====================================================================
    ! 辅助子程序: 解析 [...] 中的浮点数组，按逗号分隔
    ! ====================================================================
    subroutine parse_array_values(line, arr, n)
        character(len=*), intent(in) :: line
        real(DP), intent(out) :: arr(n)
        integer, intent(in) :: n

        character(len=512) :: work
        integer :: i, p1, p2

        ! 截取 [ 到 ] 之间的内容
        p1 = index(line, '[')
        p2 = index(line, ']', back=.true.)
        if (p1 == 0 .or. p2 == 0) return
        work = adjustl(line(p1+1:p2-1))

        do i = 1, n
            p2 = index(work, ',')
            if (p2 > 0) then
                read(work(1:p2-1), *) arr(i)
                work = adjustl(work(p2+1:))
            else
                read(work, *) arr(i)
                exit
            end if
        end do
    end subroutine parse_array_values

    ! ====================================================================
    ! 辅助子程序: 去除末尾逗号
    ! ====================================================================
    subroutine remove_trailing_comma(s)
        character(len=*), intent(inout) :: s
        integer :: n
        n = len_trim(s)
        if (n > 0 .and. s(n:n) == ',') s(n:n) = ' '
    end subroutine remove_trailing_comma

    ! ====================================================================
    ! 辅助函数: 判断某行是否只包含指定的闭括号 (忽略空白)
    ! ====================================================================
    function is_closing_bracket(line, bracket) result(yes)
        character(len=*), intent(in) :: line, bracket
        logical :: yes
        character(len=:), allocatable :: trimmed
        integer :: n
        trimmed = trim(line)
        n = len(trimmed)
        ! 允许末尾有逗号 (JSON 数组元素分隔)
        if (n >= 1) then
            if (trimmed(n:n) == ',') trimmed = trimmed(1:n-1)
        end if
        yes = (trim(trimmed) == bracket)
    end function is_closing_bracket

    ! ====================================================================
    ! 对单个 GMM 执行采样、统计对比、输出报告和散点 CSV
    ! ====================================================================
    subroutine test_gmm_step(label, gmm)
        character(len=*), intent(in) :: label
        type(uq_gmm_state_type), intent(in) :: gmm

        real(DP), allocatable :: samples(:, :)
        integer, allocatable :: num_per_comp(:)
        real(DP) :: emp_mean(DIM), emp_cov(DIM, DIM)
        real(DP) :: frob_norm, frob_sum, diff_val
        integer :: comp_idx, start_idx, end_idx, n_actual, ncomp
        integer :: remainder, max_weight_idx, p_count, i, j, k
        real(DP) :: max_weight
        character(len=256) :: fn
        character(len=20) :: dim_labels(DIM) = &
            ['x  ', 'y  ', 'z  ', 'vx ', 'vy ', 'vz ']
        integer :: u_report, u_csv, ios

        ncomp = gmm%n_components
        allocate(num_per_comp(ncomp))

        ! ---- 采样 ----
        call filter%set_n_particles(N_PARTICLES)
        filter%gmm_state = gmm
        allocate(samples(DIM, N_PARTICLES))
        call filter%sample_particles_from_gmm(samples)

        ! ---- 粒子数分配 (与 sample_particles_from_gmm 相同逻辑) ----
        p_count = 0
        max_weight = -1.0_DP
        max_weight_idx = 1
        do comp_idx = 1, ncomp
            num_per_comp(comp_idx) = &
                nint(gmm%components(comp_idx)%weight * real(N_PARTICLES, DP))
            p_count = p_count + num_per_comp(comp_idx)
            if (gmm%components(comp_idx)%weight > max_weight) then
                max_weight = gmm%components(comp_idx)%weight
                max_weight_idx = comp_idx
            end if
        end do
        remainder = N_PARTICLES - p_count
        num_per_comp(max_weight_idx) = num_per_comp(max_weight_idx) + remainder

        ! ---- 统计报告 ----
        fn = OUT_DIR // 'test_gmm_sampling_report.txt'
        if (label == 'AFTER_TIME_UPDATE') then
            open(newunit=u_report, file=trim(fn), status='replace', &
                 action='write', iostat=ios)
        else
            open(newunit=u_report, file=trim(fn), status='old', &
                 position='append', action='write', iostat=ios)
        end if

        write(u_report, '(A)') '============================================================'
        write(u_report, '(A,A)') '  GMM 采样准确性测试 - ', trim(label)
        write(u_report, '(A,I0)') '  粒子总数: ', N_PARTICLES
        write(u_report, '(A,I0)') '  分量数:   ', ncomp
        write(u_report, '(A)') '============================================================'
        write(u_report, '(A)') ''

        do comp_idx = 1, ncomp
            start_idx = 1
            do k = 1, comp_idx - 1
                start_idx = start_idx + num_per_comp(k)
            end do
            end_idx = start_idx + num_per_comp(comp_idx) - 1
            n_actual = num_per_comp(comp_idx)

            if (n_actual <= 0) then
                write(u_report, '(A,I0,A)') '[分量 ', comp_idx, '] 无粒子，跳过'
                cycle
            end if

            ! 经验均值
            emp_mean = 0.0_DP
            do j = 1, n_actual
                emp_mean(:) = emp_mean(:) + samples(:, start_idx + j - 1)
            end do
            emp_mean = emp_mean / real(n_actual, DP)

            ! 经验协方差 (无偏估计)
            emp_cov = 0.0_DP
            do j = 1, n_actual
                do k = 1, DIM
                    emp_cov(:, k) = emp_cov(:, k) &
                        + (samples(:, start_idx + j - 1) - emp_mean) &
                        * (samples(k, start_idx + j - 1) - emp_mean(k))
                end do
            end do
            emp_cov = emp_cov / real(n_actual - 1, DP)

            ! ---- 写入报告 ----
            write(u_report, '(A)') &
                '------------------------------------------------------------'
            write(u_report, '(A,I0)') '  分量 ', comp_idx
            write(u_report, '(A)') &
                '------------------------------------------------------------'
            write(u_report, '(A,ES15.6,A,ES15.6)') &
                '  理论权重: ', gmm%components(comp_idx)%weight, &
                '  实际占比: ', real(n_actual, DP) / real(N_PARTICLES, DP)

            write(u_report, '(A)') '  --- 均值对比 (维度: theo | emp | diff) ---'
            do i = 1, DIM
                diff_val = emp_mean(i) - gmm%components(comp_idx)%mean(i)
                write(u_report, '(A,A3,A,ES20.12,A,ES20.12,A,ES15.6)') &
                    '    ', dim_labels(i), ': ', &
                    gmm%components(comp_idx)%mean(i), ' | ', &
                    emp_mean(i), ' | ', diff_val
            end do

            write(u_report, '(A)') '  --- 协方差差异矩阵 (emp - theo) ---'
            frob_norm = 0.0_DP; frob_sum = 0.0_DP
            do i = 1, DIM
                write(u_report, '(A,6ES14.5)') '    ', &
                    (emp_cov(i, k) - gmm%components(comp_idx)%cov(i, k), k=1, DIM)
                do k = 1, DIM
                    frob_norm = frob_norm &
                        + (emp_cov(i, k) - gmm%components(comp_idx)%cov(i, k))**2
                    frob_sum = frob_sum &
                        + abs(gmm%components(comp_idx)%cov(i, k))
                end do
            end do
            frob_norm = sqrt(frob_norm)
            write(u_report, '(A,ES15.6)') '  Frobenius 差异: ', frob_norm
            if (frob_sum > 1.0D-30) then
                write(u_report, '(A,ES15.6)') '  相对 Frobenius: ', &
                    frob_norm / (frob_sum / real(DIM * DIM, DP))
            end if
            write(u_report, '(A)') ''
        end do
        close(u_report)
        write(*,*) '[报告] ', trim(fn)

        ! ---- 散点 CSV ----
        if (label == 'AFTER_TIME_UPDATE') then
            fn = OUT_DIR // 'test_gmm_sampling_step1_scatter.csv'
        else
            fn = OUT_DIR // 'test_gmm_sampling_step2_scatter.csv'
        end if

        open(newunit=u_csv, file=trim(fn), status='replace', &
             action='write', iostat=ios)
        write(u_csv, '(A)') 'x,y,z,vx,vy,vz,comp_id'

        do comp_idx = 1, ncomp
            start_idx = 1
            do k = 1, comp_idx - 1
                start_idx = start_idx + num_per_comp(k)
            end do
            end_idx = start_idx + num_per_comp(comp_idx) - 1
            n_actual = num_per_comp(comp_idx)
            if (n_actual <= 0) cycle

            do j = 1, n_actual
                write(u_csv, '(6(ES20.12,A),I0)') &
                    samples(1, start_idx + j - 1), ',', &
                    samples(2, start_idx + j - 1), ',', &
                    samples(3, start_idx + j - 1), ',', &
                    samples(4, start_idx + j - 1), ',', &
                    samples(5, start_idx + j - 1), ',', &
                    samples(6, start_idx + j - 1), ',', comp_idx
            end do
        end do
        close(u_csv)
        write(*,*) '[CSV]   ', trim(fn)

        deallocate(samples, num_per_comp)
    end subroutine test_gmm_step

end program test_sample_particles_from_gmm
