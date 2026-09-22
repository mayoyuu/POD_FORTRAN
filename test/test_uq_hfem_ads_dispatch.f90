program test_uq_hfem_ads_dispatch
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_uq_propagation, only: run_uq_propagation, METHOD_ADS
    use pod_uq_state_module, only: uq_state_type
    use pod_uq_hfem_ads_module, only: hfem_ads_options_type, hfem_ads_stats_type
    use pod_da_force_model_module, only: cleanup_gravity_network
    implicit none

    real(DP) :: nominal(6), covariance(6,6)
    integer :: failures

    failures = 0
    nominal = [1000.0_DP,-2000.0_DP,3000.0_DP,1.0_DP,-2.0_DP,3.0_DP]
    covariance = 0.0_DP
    covariance(1,1) = 1.0_DP
    covariance(2,2) = 1.0_DP
    covariance(3,3) = 1.0_DP
    covariance(4,4) = 1.0e-8_DP
    covariance(5,5) = 1.0e-8_DP
    covariance(6,6) = 1.0e-8_DP
    call pod_engine_init('config/config.txt')
    call run_case(0.0_DP, 6, failures)
    call run_case(0.02_DP, 7, failures)
    call cleanup_gravity_network()
    if (failures /= 0) then
        write(*,'(A,I0)') 'FAIL: unified HFEM ADS dispatch checks: ', failures
        stop 1
    end if
    write(*,'(A)') 'PASS: unified UQ dispatch preserves HFEM ADS clouds'

contains

    subroutine run_case(srp_sigma, expected_dimension, n_fail)
        real(DP), intent(in) :: srp_sigma
        integer, intent(in) :: expected_dimension
        integer, intent(inout) :: n_fail
        type(hfem_ads_options_type) :: options
        type(hfem_ads_stats_type) :: stats
        type(uq_state_type) :: initial_state, final_state

        options%srp_sigma = srp_sigma
        call run_uq_propagation(nominal, covariance, 0.0_DP, 0.0_DP, 0.0_DP, &
            METHOD_ADS, 16, .false., initial_state, final_state, da_order=2, &
            ads_options=options, ads_stats=stats)
        call assert_true(allocated(final_state%samples), &
            'dispatcher returns final particles', n_fail)
        call assert_true(size(final_state%samples,1) == expected_dimension, &
            'dispatcher selects six or seven dimensions', n_fail)
        call assert_true(size(final_state%samples,2) == 16, &
            'dispatcher retains requested particle count', n_fail)
        call assert_true(allocated(final_state%mean) .and. allocated(final_state%cov), &
            'dispatcher computes final moments', n_fail)
        call assert_true(stats%requested_count == 16 .and. &
            stats%propagated_count == 16, 'ADS stats retain particle count', n_fail)
        call initial_state%deallocate_memory()
        call final_state%deallocate_memory()
    end subroutine run_case

    subroutine assert_true(condition, label, n_fail)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        integer, intent(inout) :: n_fail
        if (.not. condition) then
            write(*,'(A,A)') 'FAIL: ', trim(label)
            n_fail = n_fail + 1
        end if
    end subroutine assert_true

end program test_uq_hfem_ads_dispatch
