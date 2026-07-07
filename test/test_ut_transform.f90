!> @file test_ut_transform.f90
!> @brief 纯数学验证 pod_uq_transform_ut_module 的核心函数
!>   Test 1: Sigma 点自洽性 (3D, 非对角协方差)
!>   Test 2: 线性变换与解析解对比 (2D)
!>   Test 3: 非线性极坐标变换 + MC 参考对比 (2D)
!>   Test 4: 角度环绕修复验证 (is_angle)
program test_ut_transform
    use pod_global, only: DP
    use pod_uq_transform_ut_module, only: generate_sigma_points, reconstruct_ut_moments
    implicit none

    interface
        subroutine dpotrf(uplo, n, a, lda, info)
            integer, parameter :: DP2 = selected_real_kind(15, 307)
            character(len=1), intent(in) :: uplo
            integer, intent(in)          :: n
            real(DP2), intent(inout)     :: a(lda, *)
            integer, intent(in)          :: lda
            integer, intent(out)         :: info
        end subroutine dpotrf
    end interface

    integer :: passed, failed
    passed = 0
    failed = 0

    write(*,*) '========================================'
    write(*,*) '  UT Transform 单元测试'
    write(*,*) '========================================'

    call test_sigma_self_consistency(passed, failed)
    call test_linear_transform(passed, failed)
    call test_nonlinear_vs_mc(passed, failed)
    call test_angle_wraparound(passed, failed)

    write(*,*) '========================================'
    write(*,'(A,I0,A,I0,A)') '  结果: ', passed, ' PASSED, ', failed, ' FAILED'
    write(*,*) '========================================'

    if (failed > 0) stop 1

contains

    !> ===================================================================
    !> Test 1: Sigma 点自洽性
    !> 用恒等映射验证 generate + reconstruct 能精确恢复输入
    !> ===================================================================
    subroutine test_sigma_self_consistency(passed, failed)
        integer, intent(inout) :: passed, failed

        real(DP) :: mean_x(3), cov_x(3,3), kappa
        real(DP), allocatable :: sigmas(:,:), weights(:)
        real(DP), allocatable :: mean_y(:), P_yy(:,:), P_xy(:,:)
        integer :: i
        logical :: ok

        write(*,*) '--- Test 1: Sigma 点自洽性 (3D) ---'

        mean_x = [1.0_DP, 2.0_DP, 3.0_DP]
        cov_x = reshape([ &
            4.0_DP, 1.2_DP, 0.8_DP, &
            1.2_DP, 3.0_DP, 0.6_DP, &
            0.8_DP, 0.6_DP, 2.0_DP  &
        ], [3, 3])
        kappa = 0.5_DP

        call generate_sigma_points(mean_x, cov_x, kappa, sigmas, weights)
        call reconstruct_ut_moments(sigmas, sigmas, weights, mean_y, P_yy, P_xy)

        ok = .true.

        ! 验证均值恢复
        do i = 1, 3
            if (abs(mean_y(i) - mean_x(i)) > 1.0e-12_DP) then
                write(*,'(A,I1,A,ES15.6,A,ES15.6)') '  FAIL: mean(', i, ') = ', mean_y(i), ' expected ', mean_x(i)
                ok = .false.
            end if
        end do

        ! 验证协方差恢复
        do i = 1, 3
            if (abs(P_yy(i,i) - cov_x(i,i)) > 1.0e-12_DP) then
                write(*,'(A,I1,A,ES15.6,A,ES15.6)') '  FAIL: cov(', i, ',', i, ') = ', P_yy(i,i), ' expected ', cov_x(i,i)
                ok = .false.
            end if
        end do
        if (abs(P_yy(1,2) - cov_x(1,2)) > 1.0e-12_DP) then
            write(*,'(A,ES15.6,A,ES15.6)') '  FAIL: cov(1,2) = ', P_yy(1,2), ' expected ', cov_x(1,2)
            ok = .false.
        end if

        if (ok) then
            write(*,*) '  PASS'
            passed = passed + 1
        else
            failed = failed + 1
        end if

        deallocate(sigmas, weights, mean_y, P_yy, P_xy)
    end subroutine test_sigma_self_consistency

    !> ===================================================================
    !> Test 2: 线性变换 f(x) = A*x + b
    !> UT 结果应与线性误差传播理论值一致
    !> ===================================================================
    subroutine test_linear_transform(passed, failed)
        integer, intent(inout) :: passed, failed

        real(DP) :: mean_x(2), cov_x(2,2), kappa
        real(DP) :: A(2,2), b(2)
        real(DP), allocatable :: sigmas(:,:), weights(:), sigmas_y(:,:)
        real(DP), allocatable :: mean_y(:), P_yy(:,:), P_xy(:,:)
        real(DP) :: expected_mean(2), expected_cov(2,2)
        integer :: i
        logical :: ok

        write(*,*) '--- Test 2: 线性变换 f(x)=A*x+b (2D) ---'

        mean_x = [2.0_DP, -1.0_DP]
        cov_x = reshape([1.5_DP, 0.3_DP, 0.3_DP, 0.8_DP], [2, 2])
        A = reshape([2.0_DP, 0.0_DP, 0.0_DP, 3.0_DP], [2, 2])
        b = [1.0_DP, 1.0_DP]
        kappa = 0.5_DP

        expected_mean = matmul(A, mean_x) + b   ! [5, -2]
        expected_cov  = matmul(A, matmul(cov_x, transpose(A)))  ! [[6, 0.9], [0.9, 7.2]]

        call generate_sigma_points(mean_x, cov_x, kappa, sigmas, weights)

        allocate(sigmas_y(2, size(sigmas, 2)))
        do i = 1, size(sigmas, 2)
            sigmas_y(:, i) = matmul(A, sigmas(:, i)) + b
        end do

        call reconstruct_ut_moments(sigmas, sigmas_y, weights, mean_y, P_yy, P_xy)

        ok = .true.

        do i = 1, 2
            if (abs(mean_y(i) - expected_mean(i)) > 1.0e-12_DP) then
                write(*,'(A,I1,A,ES15.6,A,ES15.6)') '  FAIL: mean(', i, ') = ', mean_y(i), ' expected ', expected_mean(i)
                ok = .false.
            end if
        end do

        do i = 1, 2
            if (abs(P_yy(i,i) - expected_cov(i,i)) > 1.0e-12_DP) then
                write(*,'(A,I1,A,ES15.6,A,ES15.6)') '  FAIL: cov(', i, ',', i, ') = ', P_yy(i,i), ' expected ', expected_cov(i,i)
                ok = .false.
            end if
        end do
        if (abs(P_yy(1,2) - expected_cov(1,2)) > 1.0e-12_DP) then
            write(*,'(A,ES15.6,A,ES15.6)') '  FAIL: cov(1,2) = ', P_yy(1,2), ' expected ', expected_cov(1,2)
            ok = .false.
        end if

        if (ok) then
            write(*,*) '  PASS'
            passed = passed + 1
        else
            failed = failed + 1
        end if

        deallocate(sigmas, weights, sigmas_y, mean_y, P_yy, P_xy)
    end subroutine test_linear_transform

    !> ===================================================================
    !> Test 3: 非线性极坐标变换 + MC 参考对比
    !> f([x,y]) = [r, theta], 用 Monte Carlo 验证 UT 精度
    !> ===================================================================
    subroutine test_nonlinear_vs_mc(passed, failed)
        integer, intent(inout) :: passed, failed

        real(DP) :: mean_x(2), cov_x(2,2), kappa
        real(DP), allocatable :: sigmas(:,:), weights(:), sigmas_y(:,:)
        real(DP), allocatable :: mean_y(:), P_yy(:,:), P_xy(:,:)
        real(DP), parameter :: PI = 3.14159265358979323846_DP

        integer, parameter :: N_MC = 100000
        real(DP) :: mc_samples(2, N_MC), mc_r(N_MC), mc_theta(N_MC)
        real(DP) :: mc_mean(2), mc_cov(2,2), L(2,2)
        real(DP) :: r_scale
        integer :: i, info
        logical :: ok

        write(*,*) '--- Test 3: 非线性极坐标变换 vs MC (2D) ---'

        ! 输入: 远离原点的分布，避免 r 接近 0 的奇异性
        mean_x = [100.0_DP, 0.0_DP]
        cov_x = reshape([4.0_DP, 0.0_DP, 0.0_DP, 9.0_DP], [2, 2])
        kappa = 0.5_DP

        ! --- UT ---
        call generate_sigma_points(mean_x, cov_x, kappa, sigmas, weights)

        allocate(sigmas_y(2, size(sigmas, 2)))
        do i = 1, size(sigmas, 2)
            sigmas_y(1, i) = sqrt(sigmas(1,i)**2 + sigmas(2,i)**2)          ! r
            sigmas_y(2, i) = atan2(sigmas(2,i), sigmas(1,i))                ! theta
        end do

        call reconstruct_ut_moments(sigmas, sigmas_y, weights, mean_y, P_yy, P_xy)

        ! --- MC ---
        L = cov_x
        call dpotrf('L', 2, L, 2, info)
        if (info /= 0) then
            write(*,*) '  FAIL: Cholesky 分解失败 (MC 采样)'
            failed = failed + 1
            deallocate(sigmas, weights, sigmas_y, mean_y, P_yy, P_xy)
            return
        end if
        L(1, 2) = 0.0_DP

        call random_seed()
        call random_number(mc_samples)
        do i = 1, N_MC
            ! Box-Muller 生成 2 个独立 N(0,1)
            r_scale = sqrt(-2.0_DP * log(mc_samples(1,i)))
            mc_samples(1,i) = r_scale * cos(2.0_DP * PI * mc_samples(2,i))
            mc_samples(2,i) = r_scale * sin(2.0_DP * PI * mc_samples(2,i))
        end do

        ! 转换为目标分布的样本: mean_x + L * z
        do i = 1, N_MC
            mc_samples(:, i) = mean_x + matmul(L, mc_samples(:, i))
            mc_r(i)     = sqrt(mc_samples(1,i)**2 + mc_samples(2,i)**2)
            mc_theta(i) = atan2(mc_samples(2,i), mc_samples(1,i))
        end do

        mc_mean(1) = sum(mc_r)     / real(N_MC, DP)
        mc_mean(2) = sum(mc_theta) / real(N_MC, DP)

        mc_cov = 0.0_DP
        do i = 1, N_MC
            mc_cov(1,1) = mc_cov(1,1) + (mc_r(i) - mc_mean(1))**2
            mc_cov(2,2) = mc_cov(2,2) + (mc_theta(i) - mc_mean(2))**2
            mc_cov(1,2) = mc_cov(1,2) + (mc_r(i) - mc_mean(1)) * (mc_theta(i) - mc_mean(2))
        end do
        mc_cov = mc_cov / real(N_MC - 1, DP)
        mc_cov(2,1) = mc_cov(1,2)

        ok = .true.

        ! 均值: UT 与 MC 的误差应小于 MC 统计误差 (~sigma/sqrt(N))
        write(*,'(A,ES12.4,A,ES12.4)') '  UT  mean_r = ', mean_y(1), ', mean_theta = ', mean_y(2)
        write(*,'(A,ES12.4,A,ES12.4)') '  MC  mean_r = ', mc_mean(1), ', mean_theta = ', mc_mean(2)

        if (abs(mean_y(1) - mc_mean(1)) > 3.0_DP * sqrt(mc_cov(1,1)/real(N_MC, DP))) then
            write(*,*) '  FAIL: mean_r 超出 3-sigma MC 误差'
            ok = .false.
        end if

        if (abs(mean_y(2) - mc_mean(2)) > 3.0_DP * sqrt(mc_cov(2,2)/real(N_MC, DP))) then
            write(*,*) '  FAIL: mean_theta 超出 3-sigma MC 误差'
            ok = .false.
        end if

        ! 协方差: 检查对角线元素的相对误差
        write(*,'(A,ES12.4,A,ES12.4)') '  UT  std_r = ', sqrt(P_yy(1,1)), ', std_theta = ', sqrt(P_yy(2,2))
        write(*,'(A,ES12.4,A,ES12.4)') '  MC  std_r = ', sqrt(mc_cov(1,1)), ', std_theta = ', sqrt(mc_cov(2,2))

        if (abs(P_yy(1,1) - mc_cov(1,1)) / mc_cov(1,1) > 0.10_DP) then
            write(*,'(A,F6.2,A)') '  FAIL: r 方差相对误差 ', abs(P_yy(1,1)-mc_cov(1,1))/mc_cov(1,1)*100, '% > 10%'
            ok = .false.
        end if

        if (ok) then
            write(*,*) '  PASS'
            passed = passed + 1
        else
            failed = failed + 1
        end if

        deallocate(sigmas, weights, sigmas_y, mean_y, P_yy, P_xy)
    end subroutine test_nonlinear_vs_mc

    !> ===================================================================
    !> Test 4: 角度环绕修复验证
    !> 手动构造跨 0/360 边界的不对称 sigma 点集，
    !> 验证 is_angle 的 sin/cos 平均能正确恢复圆形均值
    !> ===================================================================
    subroutine test_angle_wraparound(passed, failed)
        integer, intent(inout) :: passed, failed

        real(DP), allocatable :: sigmas(:,:), weights(:), sigmas_y(:,:)
        real(DP), allocatable :: mean_no_fix(:), P_yy_nofix(:,:), P_xy_nofix(:,:)
        real(DP), allocatable :: mean_fixed(:), P_yy_fixed(:,:), P_xy_fixed(:,:)
        real(DP), parameter :: PI = 3.14159265358979323846_DP
        real(DP) :: deg2rad, correct_circular_mean_deg
        logical :: is_angle(2)
        logical :: ok

        deg2rad = PI / 180.0_DP

        write(*,*) '--- Test 4: 角度环绕修复验证 ---'

        ! 手动构造 5 个不对称 sigma 点 (2 维: r, theta)
        ! theta 取在 350°, 355°, 2°, 5°, 8° (跨 0° 边界，不对称)
        ! 权重: 中心点 350° 权重 0.5, 其余各 0.125
        ! 线性平均 ≈ 0.5*350 + 0.125*(355+2+5+8) = 175 + 46.25 = 221.25° (错误)
        ! sin/cos 平均 ≈ 354° (正确)
        allocate(sigmas(2, 5))
        allocate(weights(5))
        allocate(sigmas_y(2, 5))

        sigmas(1, 1:5) = 100.0_DP   ! radius
        sigmas(2, 1)   = 350.0_DP * deg2rad
        sigmas(2, 2)   = 355.0_DP * deg2rad
        sigmas(2, 3)   =   2.0_DP * deg2rad
        sigmas(2, 4)   =   5.0_DP * deg2rad
        sigmas(2, 5)   =   8.0_DP * deg2rad

        weights = [0.5_DP, 0.125_DP, 0.125_DP, 0.125_DP, 0.125_DP]

        sigmas_y = sigmas   ! 恒等映射

        ! 不修复角度
        call reconstruct_ut_moments(sigmas, sigmas_y, weights, mean_no_fix, P_yy_nofix, P_xy_nofix)
        ! 修复角度 (第二维标记为角度)
        is_angle = [.false., .true.]
        call reconstruct_ut_moments(sigmas, sigmas_y, weights, mean_fixed, P_yy_fixed, P_xy_fixed, is_angle)

        write(*,'(A,F8.3,A)') '  无修复均值  : ', mean_no_fix(2) * 180.0_DP / PI, ' deg'
        write(*,'(A,F8.3,A)') '  修复后均值  : ', mean_fixed(2) * 180.0_DP / PI, ' deg'

        ! 计算理论圆形均值 (使用 sin/cos 加权)
        correct_circular_mean_deg = atan2( &
            sum(weights * sin(sigmas(2, :))), &
            sum(weights * cos(sigmas(2, :))) ) * 180.0_DP / PI

        write(*,'(A,F8.3,A)') '  理论圆形均值: ', correct_circular_mean_deg, ' deg'

        ok = .true.

        ! 无修复均值应该接近 221° (线性平均的明显错误结果)
        if (abs(mean_no_fix(2) * 180.0_DP / PI - 221.25_DP) > 1.0_DP) then
            write(*,'(A,F8.3,A)') '  FAIL: 无修复均值应 ≈ 221.25°, 实际= ', mean_no_fix(2) * 180.0_DP / PI, ' deg'
            ok = .false.
        end if

        ! 修复后均值应接近理论圆形均值 (~354°)
        if (abs(mean_fixed(2) * 180.0_DP / PI - correct_circular_mean_deg) > 0.1_DP) then
            write(*,*) '  FAIL: 修复后均值与理论圆形均值偏差过大'
            ok = .false.
        end if

        ! 半径不应受 is_angle 影响
        if (abs(mean_no_fix(1) - mean_fixed(1)) > 1.0e-12_DP) then
            write(*,*) '  FAIL: is_angle 不应影响非角度维度'
            ok = .false.
        end if

        if (ok) then
            write(*,*) '  PASS'
            passed = passed + 1
        else
            failed = failed + 1
        end if

        deallocate(sigmas, weights, sigmas_y, mean_no_fix, P_yy_nofix, P_xy_nofix)
        deallocate(mean_fixed, P_yy_fixed, P_xy_fixed)
    end subroutine test_angle_wraparound

end program test_ut_transform
