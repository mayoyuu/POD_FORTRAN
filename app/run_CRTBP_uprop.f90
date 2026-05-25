program run_CRTBP_uprop
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_uq_crtbp_mc_module, only: crtbp_mc_propagate
    use pod_uq_crtbp_da_module, only: crtbp_da_propagate
    use pod_dace_classes
    implicit none

    real(DP) :: mu, t_end
    real(DP) :: nominal_state(6)
    real(DP) :: cov(6,6)
    character(len=10) :: method
    integer  :: n_samples, da_order
    real(DP) :: rel_tol, abs_tol, dt_min, dt_max
    integer  :: max_steps_int
    character(len=MAX_STRING_LEN) :: output_prefix

    real(DP), allocatable :: final_samples(:,:)
    real(DP) :: final_mean(6), final_cov(6,6), propagated_ref(6)

    integer :: i, j

    ! ---- 维度转换专用变量 ----
    real(DP) :: LU, TU, t_end_dim
    real(DP) :: pp_factor, pv_factor, vv_factor

    ! =========================================================
    !  1. 硬编码配置区域 (融合自原 make_crtbp_config)
    ! =========================================================
    mu = 0.012153614091892_DP
    LU = 384400.0_DP
    TU = 375190.26189464843_DP
    t_end_dim = 86400.0_DP*4  ! 4天的传播时间 (秒)

    ! 标称状态 (已经是无量纲)
    nominal_state = [ &
        0.622521490762424_DP, 0.725842321858840_DP, 0.0_DP, &
        -0.214221454673397_DP, 1.061689685054569_DP, 0.0_DP ]

    ! 初始化有量纲的协方差矩阵 (对角线)
    cov = 0.0_DP
    cov(1,1) = 3.0_DP
    cov(2,2) = 3.0_DP
    cov(3,3) = 3.0_DP
    cov(4,4) = 3.33e-5_DP
    cov(5,5) = 3.33e-5_DP
    cov(6,6) = 3.33e-5_DP

    ! UQ 和 积分器参数
    method = 'MC'
    n_samples = 10000
    da_order = 4

    rel_tol = 1.0e-12_DP
    abs_tol = 1.0e-12_DP
    dt_min  = 1.0e-10_DP
    dt_max  = 0.1_DP
    max_steps_int = 100000

    output_prefix = './output/crtbp_uprop_MC'

    ! =========================================================
    !  2. 自动无量纲化处理
    ! =========================================================
    ! 时间无量纲化
    t_end = t_end_dim / TU

    ! 协方差矩阵转换系数计算
    pp_factor = 1.0_DP / (LU * LU)
    pv_factor = TU / (LU * LU)
    vv_factor = (TU * TU) / (LU * LU)

    ! 执行协方差矩阵的无量纲化转换
    do i = 1, 3
        do j = 1, 3
            cov(i,j) = cov(i,j) * pp_factor
        end do
    end do
    do i = 1, 3
        do j = 4, 6
            cov(i,j) = cov(i,j) * pv_factor
            cov(j,i) = cov(j,i) * pv_factor
        end do
    end do
    do i = 4, 6
        do j = 4, 6
            cov(i,j) = cov(i,j) * vv_factor
        end do
    end do

    ! ---- 校验无量纲化后的状态量级 (可选安全机制) ----
    if (maxval(abs(nominal_state(1:3))) > 100.0_DP .or. &
        maxval(abs(nominal_state(4:6))) > 10.0_DP) then
        print *, 'WARNING: Input state appears to be dimensional.'
        print *, 'CRTBP expects non-dimensional units (positions ~O(1), velocities ~O(1)).'
        print *, 'Position magnitudes found: ', maxval(abs(nominal_state(1:3)))
        print *, 'Velocity magnitudes found: ', maxval(abs(nominal_state(4:6)))
    end if

    ! ---- 打印配置汇总信息 ----
    print *, '========================================'
    print *, '  CRTBP Uncertainty Propagation (Hardcoded)'
    print *, '========================================'
    print '(A,F12.8)',  '  mu              = ', mu
    print '(A,6F10.6)', '  initial state   = ', nominal_state
    print '(A,F14.4)',  '  propagation T   = ', t_end
    print '(A,A)',      '  method          = ', trim(method)
    if (trim(method) == 'MC') then
        print '(A,I0)', '  num_samples     = ', n_samples
    else
        print '(A,I0)', '  da_order        = ', da_order
        print '(A,I0)', '  num_samples     = ', n_samples
    end if
    print '(A,ES10.3)', '  rel_tol         = ', rel_tol
    print '(A,ES10.3)', '  abs_tol         = ', abs_tol
    print '(A,A)',      '  output prefix   = ', trim(output_prefix)
    print *, '========================================'

    ! ---- 执行蒙特卡洛或微分代数传播 ----
    if (trim(method) == 'MC') then
        call crtbp_mc_propagate(nominal_state, cov, mu, t_end, &
            n_samples, rel_tol, abs_tol, dt_min, dt_max, max_steps_int, &
            final_samples, final_mean, final_cov, .true.)
    else if (trim(method) == 'DA') then
        call dace_initialize(da_order, 6)
        call crtbp_da_propagate(nominal_state, cov, mu, t_end, &
            da_order, n_samples, rel_tol, abs_tol, dt_min, dt_max, max_steps_int, &
            final_samples, final_mean, final_cov, propagated_ref, .true.)
    else
        print *, 'Error: method must be "MC" or "DA", got: ', trim(method)
        stop
    end if

    ! ---- 写入 CSV 结果文件 ----
    call write_csv_particles(trim(output_prefix) // '_particles.csv', final_samples)
    call write_csv_stats(trim(output_prefix) // '_stats.csv', final_mean, final_cov)

    print *, 'Results saved to:'
    print *, '  ', trim(output_prefix) // '_particles.csv'
    print *, '  ', trim(output_prefix) // '_stats.csv'

    deallocate(final_samples)

contains

    ! ================================================================
    !  CSV 输出工具
    ! ================================================================

    subroutine write_csv_particles(filename, samples)
        character(len=*), intent(in) :: filename
        real(DP), intent(in) :: samples(:,:)
        integer :: unit, i

        open(newunit=unit, file=trim(filename), status='replace', action='write')
        write(unit, '(A)') 'x,y,z,vx,vy,vz'
        do i = 1, size(samples, 2)
            write(unit, '(*(ES22.14, :, ","))') samples(:, i)
        end do
        close(unit)
    end subroutine write_csv_particles

    subroutine write_csv_stats(filename, mean_vec, cov_mat)
        character(len=*), intent(in) :: filename
        real(DP), intent(in) :: mean_vec(6), cov_mat(6,6)
        integer :: unit, i

        open(newunit=unit, file=trim(filename), status='replace', action='write')
        write(unit, '(A)') '# Mean'
        write(unit, '(*(ES22.14, :, ","))') mean_vec(:)
        write(unit, '(A)') '# Covariance Matrix'
        do i = 1, 6
            write(unit, '(*(ES22.14, :, ","))') cov_mat(i, :)
        end do
        close(unit)
    end subroutine write_csv_stats

end program run_CRTBP_uprop