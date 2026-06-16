!> CRTBP 不确定性传播综合对比测试
!>
!> 对比四种方法的传播结果和耗时:
!>   MC (Monte Carlo) — 数值积分参考解
!>   DA (Differential Algebra) — 单 DA 扩展 + 采样求值
!>   STT (State Transition Tensor) — 增广积分 + STT 多项式映射
!>   ADS (Automatic Domain Splitting) — 域分裂 + 分段 DA 求值
!>
!> 基于 test_stt_crtbp.f90 的 Test 6 模式扩展。
program test_uq_crtbp_comparison
    use pod_global, only: DP
    use pod_dace_classes, only: dace_initialize
    use pod_uq_crtbp_da_module, only: crtbp_da_propagate_deviates
    use pod_uq_crtbp_mc_module, only: crtbp_mc_propagate_deviates
    use pod_uq_prop_stt_module, only: uq_stt_propagator, stt_propagate_deviates
    use pod_uq_crtbp_ads_module, only: crtbp_ads_propagate_deviates
    implicit none

    integer :: passed, failed
    real(DP) :: mu, dro_x0(6), dro_period
    real(DP), allocatable :: deviates(:,:)
    real(DP), allocatable :: mc_samples(:,:), da_samples(:,:)
    real(DP), allocatable :: stt_samples(:,:), ads_samples(:,:)
    real(DP), allocatable :: mc_dev_subset(:,:)
    real(DP) :: mc_mean(6), mc_cov(6,6)
    real(DP) :: da_mean(6), da_cov(6,6), da_ref(6)
    real(DP) :: stt_mean(6), stt_cov(6,6)
    real(DP) :: ads_mean(6), ads_cov(6,6)
    type(uq_stt_propagator) :: stt_prop
    real(DP) :: mc_time, da_time, stt_time, ads_time
    integer(8) :: t1, t2, t_rate
    integer :: test_order, n_dev, n_mc, unit, io, k, i, j
    integer :: ads_n_patches
    real(DP) :: err_toll(6)
    character(len=64) :: fname

    passed = 0
    failed = 0

    write(*,*) '========================================'
    write(*,*) '  CRTBP Uncertainty Propagation'
    write(*,*) '  MC vs DA vs STT vs ADS Comparison'
    write(*,*) '========================================'

    ! Parameters (same as test_stt_crtbp Test 6)
    test_order = 2
    n_dev = 100000
    n_mc = 100000
    mu = 0.012153614091892_DP
    dro_x0 = [1.1309265107780351_DP, 0.0_DP, 0.0_DP, &
              0.0_DP, -0.46540743845849059_DP, 0.0_DP]
    dro_period = 2.3017923284002024_DP * 1.5_DP

    write(*,'(A,I0,A,I0,A,I0)') ' Config: order=', test_order, &
        ', N_dev=', n_dev, ', N_mc=', n_mc

    ! ---- 1. Load deviates ----
    write(*,*) 'Loading deviates from rand_list200km_0.7ms.txt...'
    allocate(deviates(6, n_dev))
    open(newunit=unit, file='../rand_list200km_0.7ms.txt', status='old', action='read')
    do k = 1, n_dev
        read(unit, *, iostat=io) deviates(1, k), deviates(2, k), deviates(3, k), &
                                  deviates(4, k), deviates(5, k), deviates(6, k)
        if (io /= 0) then
            write(*,*) 'FAIL: error reading line', k
            failed = failed + 1
            close(unit)
            go to 99
        end if
    end do
    close(unit)
    write(*,*) 'Loaded', n_dev, 'deviates'

    ! ---- 2. MC reference (n_mc subset) ----
    write(*,*) 'Running MC...'
    allocate(mc_dev_subset(6, n_mc))
    mc_dev_subset = deviates(:, 1:n_mc)
    call system_clock(t1, t_rate)
    call crtbp_mc_propagate_deviates(dro_x0, mc_dev_subset, mu, dro_period, &
        1.0e-14_DP, 1.0e-14_DP, 1.0e-10_DP, 0.01_DP, 100000, &
        mc_samples, mc_mean, mc_cov, .false.)
    call system_clock(t2, t_rate)
    mc_time = real(t2 - t1, DP) / real(t_rate, DP)

    ! ---- 3. DA ----
    write(*,*) 'Running DA...'
    call dace_initialize(test_order, 6)
    call system_clock(t1, t_rate)
    call crtbp_da_propagate_deviates(dro_x0, deviates, mu, dro_period, &
        test_order, 1.0e-14_DP, 1.0e-14_DP, 1.0e-10_DP, 0.01_DP, 100000, &
        da_samples, da_mean, da_cov, da_ref, .false.)
    call system_clock(t2, t_rate)
    da_time = real(t2 - t1, DP) / real(t_rate, DP)

    ! ---- 4. STT ----
    write(*,*) 'Running STT...'
    call stt_prop%set_stt_order(test_order)
    call stt_prop%set_stt_mu(mu)
    call stt_prop%set_stt_tolerances(abs_tol=1.0e-14_DP, rel_tol=1.0e-14_DP, &
        dt_min=1.0e-10_DP, dt_max=0.01_DP, max_steps=100000)
    call stt_prop%set_verbosity(.false.)
    call system_clock(t1, t_rate)
    call stt_propagate_deviates(stt_prop, 0.0_DP, dro_period, &
        dro_x0, deviates, stt_samples, stt_mean, stt_cov)
    call system_clock(t2, t_rate)
    stt_time = real(t2 - t1, DP) / real(t_rate, DP)

    ! ---- 5. ADS ----
    write(*,*) 'Running ADS...'
    err_toll = 1.0d-4
    call system_clock(t1, t_rate)
    call crtbp_ads_propagate_deviates(dro_x0, deviates, mu, dro_period, &
        test_order, 12, err_toll, &
        1.0e-14_DP, 1.0e-14_DP, 1.0e-10_DP, 0.01_DP, 100000, &
        ads_samples, ads_mean, ads_cov, ads_n_patches, .true.)
    call system_clock(t2, t_rate)
    ads_time = real(t2 - t1, DP) / real(t_rate, DP)

    ! ---- 6. Results comparison ----
    write(*,*) '========================================'
    write(*,*) '  Results'
    write(*,*) '========================================'

    write(*,'(A,6ES14.6)') ' MC  mean:', mc_mean
    write(*,'(A,6ES14.6)') ' DA  mean:', da_mean
    write(*,'(A,6ES14.6)') ' STT mean:', stt_mean
    write(*,'(A,6ES14.6)') ' ADS mean:', ads_mean

    write(*,'(A,ES14.6)') ' |DA - MC|_max  :', maxval(abs(da_mean - mc_mean))
    write(*,'(A,ES14.6)') ' |STT - MC|_max :', maxval(abs(stt_mean - mc_mean))
    write(*,'(A,ES14.6)') ' |ADS - MC|_max :', maxval(abs(ads_mean - mc_mean))

    ! Covariance comparison
    call compare_cov('DA  vs MC', da_cov, mc_cov, 0.10_DP, passed, failed)
    call compare_cov('STT vs MC', stt_cov, mc_cov, 0.10_DP, passed, failed)
    call compare_cov('ADS vs MC', ads_cov, mc_cov, 0.15_DP, passed, failed)

    ! ---- 7. Timing ----
    write(*,*) '========================================'
    write(*,*) '  Timing'
    write(*,*) '========================================'
    write(*,'(A,F10.3,A)') ' MC  : ', mc_time,  ' s'
    write(*,'(A,F10.3,A)') ' DA  : ', da_time,  ' s'
    write(*,'(A,F10.3,A)') ' STT : ', stt_time, ' s'
    write(*,'(A,F10.3,A,I0,A)') ' ADS : ', ads_time, ' s (', ads_n_patches, ' patches)'
    write(*,'(A,F10.1)') ' Speedup vs MC (DA/MC):  ', mc_time / max(da_time, 1.0e-6_DP)
    write(*,'(A,F10.1)') ' Speedup vs MC (STT/MC): ', mc_time / max(stt_time, 1.0e-6_DP)
    write(*,'(A,F10.1)') ' Speedup vs MC (ADS/MC): ', mc_time / max(ads_time, 1.0e-6_DP)

    ! ---- 8. Output samples ----
    call write_samples('veri_crtbp_mc.txt', mc_samples, n_mc)
    write(fname, '(A,I0,A)') 'o', test_order, '_veri_crtbp_da.txt'
    call write_samples(trim(fname), da_samples, n_dev)
    write(fname, '(A,I0,A)') 'o', test_order, '_veri_crtbp_stt.txt'
    call write_samples(trim(fname), stt_samples, n_dev)
    write(fname, '(A,I0,A)') 'o', test_order, '_veri_crtbp_ads.txt'
    call write_samples(trim(fname), ads_samples, n_dev)

    ! ---- Cleanup ----
    deallocate(deviates, mc_dev_subset)
    if (allocated(mc_samples)) deallocate(mc_samples)
    if (allocated(da_samples)) deallocate(da_samples)
    if (allocated(stt_samples)) deallocate(stt_samples)
    if (allocated(ads_samples)) deallocate(ads_samples)

99  continue
    write(*,*) '========================================'
    write(*,'(A,I0,A,I0,A)') ' Results: ', passed, ' passed, ', failed, ' failed'
    write(*,*) '========================================'
    if (failed > 0) stop 1

contains

    subroutine compare_cov(label, cov1, cov2, threshold, passed, failed)
        character(len=*), intent(in) :: label
        real(DP), intent(in) :: cov1(6,6), cov2(6,6)
        real(DP), intent(in) :: threshold
        integer, intent(inout) :: passed, failed
        real(DP) :: fn1, fn2, diff_fn, rel_err
        integer :: r, c
        fn1 = 0.0_DP; fn2 = 0.0_DP; diff_fn = 0.0_DP
        do r = 1, 6
            do c = 1, 6
                fn1 = fn1 + cov1(r,c)**2
                fn2 = fn2 + cov2(r,c)**2
                diff_fn = diff_fn + (cov1(r,c) - cov2(r,c))**2
            end do
        end do
        rel_err = sqrt(diff_fn) / max(sqrt(fn2), 1.0e-15_DP)
        write(*,'(A,A,ES14.6)') ' Cov ', label, ' rel err:', rel_err
        if (rel_err < threshold) then
            write(*,'(A,A,A)') '  PASS: ', label, ' covariance matches'
            passed = passed + 1
        else
            write(*,'(A,A,A)') '  FAIL: ', label, ' covariance mismatch'
            failed = failed + 1
        end if
    end subroutine compare_cov

    subroutine write_samples(fname, samples, n)
        character(len=*), intent(in) :: fname
        real(DP), intent(in) :: samples(:,:)
        integer, intent(in) :: n
        integer :: u, k
        open(newunit=u, file=trim(fname), status='replace', action='write')
        do k = 1, n
            write(u, '(6ES22.14)') samples(:, k)
        end do
        close(u)
    end subroutine write_samples

end program test_uq_crtbp_comparison
