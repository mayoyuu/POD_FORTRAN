!> @file test_l2halo1_gap_da_ut.f90
!> @brief L2Halo-1 init OPM 固定 10 天传播: DA-MC(4) + UKF(alpha=1.0), GMM n=3/5 fitting
!> UKF 使用 pod_filter_ut_module 的 time_update，alpha=1.0
!> 输出至 output/0627_L1HALO_DAMC4_UKF10_GMM/
program test_l2halo1_gap_da_ut
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_engine_module, only: pod_engine_init
    use pod_data_format_module, only: load_initial_opm, write_json_opm
    use pod_filter_ut_module, only: ut_filter
    use pod_uq_propagation, only: run_uq_propagation, METHOD_DA
    use pod_uq_state_module, only: uq_state_type
    use pod_uq_gmm_state_module, only: uq_gmm_state_type
    use pod_gmm_math_module, only: fit_gmm_to_particles
    use pod_dace_classes, only: dace_initialize
    implicit none

    character(len=MAX_STRING_LEN) :: opm_file, config_file
    character(len=MAX_STRING_LEN) :: output_dir
    real(DP) :: et0, mean0(6), cov0(6, 6), target_et, dt
    real(DP) :: ukf_mean(6), ukf_cov(6, 6), noise_Q(6, 6)
    type(uq_state_type) :: init_da, final_da
    type(uq_gmm_state_type) :: gmm3, gmm5
    type(ut_filter) :: my_filter
    integer :: i

    config_file = 'config/config.txt'
    opm_file    = 'OPM/L2Halo-1/L2Halo-1_init.opm.json'
    output_dir  = 'output/0627_L1HALO_DAMC4_UKF10_GMM'

    write(*,*) '=================================================='
    write(*,*) '  L2Halo-1 固定 10 天传播: DA-MC(4) vs UKF(a=1.0)'
    write(*,*) '  输出: ', trim(output_dir)
    write(*,*) '=================================================='

    call pod_engine_init(trim(config_file))

    call load_initial_opm(trim(opm_file), et0, mean0, cov0)
    write(*,*) '初始历元 et0 = ', et0
    write(*,*) '初始状态 = ', mean0

    dt = 10.0_DP * 86400.0_DP
    target_et = et0 + dt

    write(*,*) ''
    write(*,*) '========== 传播时间信息 =========='
    write(*,*) '初始历元: ', et0
    write(*,*) '目标历元: ', target_et
    write(*,*) '传播时间: ', dt/86400.0_DP, ' 天'
    write(*,*) '=================================='

    call system('mkdir -p ' // trim(output_dir))

    ! ===== DA-MC(4) 传播: 100k 粒子 =====
    write(*,*) ''
    write(*,*) '--- DA-MC(4) 传播 (100k 粒子, 10 天) ---'

    call dace_initialize(4, 6)

    call run_uq_propagation(mean0, cov0, et0, 0.0_DP, dt, &
                            METHOD_DA, 100000, .false., &
                            init_da, final_da, da_order=4)

    write(*,*) 'DA-MC 最终位置标准差: ', sqrt(final_da%cov(1,1)), &
                                         sqrt(final_da%cov(2,2)), sqrt(final_da%cov(3,3))
    write(*,*) 'DA-MC 最终速度标准差: ', sqrt(final_da%cov(4,4)), &
                                         sqrt(final_da%cov(5,5)), sqrt(final_da%cov(6,6))

    call write_scatter_txt(trim(output_dir) // '/L2Halo-1_DA_MC4_particles.txt', &
                           final_da%samples, size(final_da%samples, 2))
    call write_mean_cov_txt(trim(output_dir) // '/L2Halo-1_DA_MC4_mean_cov.txt', &
                            final_da%mean, final_da%cov)

    call gaussianity_check(final_da%samples, size(final_da%samples, 2), &
                           final_da%mean, final_da%cov, 'DA-MC(4)')

    call write_json_opm(trim(output_dir) // '/L2Halo-1_DA_MC4', &
                        final_da%mean, final_da%cov, &
                        rms=0.0_DP, obj_id='L2Halo-1', et_last=target_et)

    ! ===== GMM 拟合 n=3, n=5 =====
    write(*,*) ''
    write(*,*) '--- GMM 拟合 (n=3, n=5) ---'

    call gmm3%allocate_components(3, 6)
    call fit_gmm_to_particles(final_da%samples, gmm3, 50, 1.0e-4_DP)
    write(*,*) 'GMM n=3 完成'
    do i = 1, 3
        write(*,*) '  comp', i, ' weight=', gmm3%components(i)%weight
    end do
    call write_json_opm(trim(output_dir) // '/L2Halo-1_DA_MC4_GMM3', &
                        final_da%mean, final_da%cov, gmm_state=gmm3, &
                        rms=0.0_DP, obj_id='L2Halo-1', et_last=target_et)

    call gmm5%allocate_components(5, 6)
    call fit_gmm_to_particles(final_da%samples, gmm5, 50, 1.0e-4_DP)
    write(*,*) 'GMM n=5 完成'
    do i = 1, 5
        write(*,*) '  comp', i, ' weight=', gmm5%components(i)%weight
    end do
    call write_json_opm(trim(output_dir) // '/L2Halo-1_DA_MC4_GMM5', &
                        final_da%mean, final_da%cov, gmm_state=gmm5, &
                        rms=0.0_DP, obj_id='L2Halo-1', et_last=target_et)

    ! ===== UKF time_update (alpha=1.0) =====
    write(*,*) ''
    write(*,*) '--- UKF time_update (alpha=1.0, 10 天) ---'

    call my_filter%filter_init(et0, mean0, cov0, alpha=1.0_DP)
    write(*,*) 'UKF alpha = 1.0, 滤波器已初始化'

    noise_Q = 0.0_DP
    call my_filter%time_update(target_et, noise_Q)

    call my_filter%get_current_state(ukf_mean)
    call my_filter%get_current_cov(ukf_cov)

    write(*,*) 'UKF 最终位置标准差: ', sqrt(ukf_cov(1,1)), &
                                     sqrt(ukf_cov(2,2)), sqrt(ukf_cov(3,3))
    write(*,*) 'UKF 最终速度标准差: ', sqrt(ukf_cov(4,4)), &
                                     sqrt(ukf_cov(5,5)), sqrt(ukf_cov(6,6))

    call write_mean_cov_txt(trim(output_dir) // '/L2Halo-1_UKF10_mean_cov.txt', &
                            ukf_mean, ukf_cov)

    call write_json_opm(trim(output_dir) // '/L2Halo-1_UKF10', &
                        ukf_mean, ukf_cov, &
                        rms=0.0_DP, obj_id='L2Halo-1', et_last=target_et)

    ! ===== 对比 =====
    write(*,*) ''
    call compare_da_ukf(final_da%mean, final_da%cov, ukf_mean, ukf_cov)

    call init_da%deallocate_memory()
    call final_da%deallocate_memory()

    write(*,*) ''
    write(*,*) '=================================================='
    write(*,*) '  完成！结果已保存至: ', trim(output_dir)
    write(*,*) '=================================================='

contains

    subroutine write_scatter_txt(filename, samples, n)
        character(len=*), intent(in) :: filename
        real(DP), intent(in) :: samples(:, :)
        integer, intent(in) :: n
        integer :: u, p
        open(newunit=u, file=trim(filename), status='replace', action='write')
        do p = 1, n
            write(u, '(I0, 6(1X, ES25.16E3))') p, samples(:, p)
        end do
        close(u)
        write(*,*) '已保存: ', trim(filename)
    end subroutine write_scatter_txt

    subroutine write_mean_cov_txt(filename, mean_vec, cov_mat)
        character(len=*), intent(in) :: filename
        real(DP), intent(in) :: mean_vec(6), cov_mat(6, 6)
        integer :: u, j
        open(newunit=u, file=trim(filename), status='replace', action='write')
        write(u, '(A)') '# Mean'
        write(u, '(6(ES25.16E3, 1X))') mean_vec
        write(u, '(A)') '# Covariance'
        do j = 1, 6
            write(u, '(6(ES25.16E3, 1X))') cov_mat(j, :)
        end do
        close(u)
        write(*,*) '已保存: ', trim(filename)
    end subroutine write_mean_cov_txt

    subroutine gaussianity_check(samples, n, mean_vec, cov_mat, label)
        real(DP), intent(in) :: samples(:, :), mean_vec(6), cov_mat(6, 6)
        integer, intent(in) :: n
        character(len=*), intent(in) :: label
        real(DP) :: diff(6), skew(6), ekurt(6), m2(6), m3(6), m4(6)
        real(DP) :: L(6, 6), y(6), md, md_mean, md_var, md_min, md_max
        integer :: i, j, k

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

        write(*,*) '[' // trim(label) // '] 高斯性:'
        write(*, '(A, F10.4, F12.4)') '   x     ', skew(1), ekurt(1)
        write(*, '(A, F10.4, F12.4)') '   y     ', skew(2), ekurt(2)
        write(*, '(A, F10.4, F12.4)') '   z     ', skew(3), ekurt(3)
        write(*, '(A, F10.4, F12.4)') '   vx    ', skew(4), ekurt(4)
        write(*, '(A, F10.4, F12.4)') '   vy    ', skew(5), ekurt(5)
        write(*, '(A, F10.4, F12.4)') '   vz    ', skew(6), ekurt(6)

        L = 0.0_DP
        do j = 1, 6
            do i = j, 6
                L(i, j) = cov_mat(i, j)
                do k = 1, j - 1
                    L(i, j) = L(i, j) - L(i, k) * L(j, k)
                end do
                if (i == j) then
                    if (L(j, j) <= 0.0_DP) then
                        write(*,*) '[WARN] 协方差非正定'
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

        write(*,*) '  马氏距离:'
        write(*, '(A, F8.3, A, F8.3, A, F8.3, A, F8.3)') &
            '    mean=', md_mean, '  std=', sqrt(md_var), &
            '  min=', md_min, '  max=', md_max
    end subroutine gaussianity_check

    subroutine compare_da_ukf(da_mean, da_cov, ukf_mean, ukf_cov)
        real(DP), intent(in) :: da_mean(6), da_cov(6, 6), ukf_mean(6), ukf_cov(6, 6)
        real(DP) :: da_std(6), ukf_std(6), ratio(6), frob_diff, frob_da
        real(DP) :: da_trace, ukf_trace, mean_diff(6)
        integer :: j
        character(len=6), parameter :: labels(6) = &
            ['x     ', 'y     ', 'z     ', 'vx    ', 'vy    ', 'vz    ']

        write(*,*) '--- DA-MC(4) vs UKF(alpha=1.0) ---'
        write(*,*) '  维度   DA-MC std      UKF std       比值(UKF/DA)'
        do j = 1, 6
            da_std(j) = sqrt(da_cov(j, j))
            ukf_std(j) = sqrt(ukf_cov(j, j))
            ratio(j) = ukf_std(j) / da_std(j)
            write(*, '(A, A, 2ES14.5, F10.4)') '  ', labels(j), da_std(j), ukf_std(j), ratio(j)
        end do

        mean_diff = da_mean - ukf_mean
        frob_diff = sqrt(sum((da_cov - ukf_cov)**2))
        frob_da   = sqrt(sum(da_cov**2))
        da_trace   = da_cov(1,1)+da_cov(2,2)+da_cov(3,3)+da_cov(4,4)+da_cov(5,5)+da_cov(6,6)
        ukf_trace  = ukf_cov(1,1)+ukf_cov(2,2)+ukf_cov(3,3)+ukf_cov(4,4)+ukf_cov(5,5)+ukf_cov(6,6)

        write(*,*)
        write(*, '(A, 6(ES10.2, 1X))') '  均值差 (DA-UKF): ', mean_diff
        write(*, '(A, ES10.3)') '  协方差 Frobenius 差异 : ', frob_diff
        write(*, '(A, F8.4)')   '  相对差异 (diff/||DA||): ', frob_diff / frob_da
        write(*, '(A, ES14.5)') '  DA-MC 协方差迹         : ', da_trace
        write(*, '(A, ES14.5)') '  UKF   协方差迹         : ', ukf_trace
    end subroutine compare_da_ukf

end program test_l2halo1_gap_da_ut
