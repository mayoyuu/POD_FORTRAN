program test_hfem_ads_zero_duration_cloud
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_uq_ads_coordinates_module, only: ADS_COORD_COMPONENT, &
        ADS_COORD_WHITENED, ads_coordinate_map_type, ads_build_coordinate_map, &
        ads_unit_to_physical
    use pod_uq_hfem_ads_module, only: hfem_ads_options_type, &
        hfem_ads_stats_type, hfem_ads_propagate
    use pod_da_force_model_module, only: cleanup_gravity_network
    implicit none

    real(DP) :: nominal(6), covariance(6,6)
    integer :: failures

    failures = 0
    nominal = [1000.0_DP, -2000.0_DP, 3000.0_DP, &
        1.0_DP, -2.0_DP, 3.0_DP]
    covariance = 0.0_DP
    covariance(1,1) = 4.0_DP
    covariance(2,2) = 9.0_DP
    covariance(3,3) = 16.0_DP
    covariance(4,4) = 1.0e-8_DP
    covariance(5,5) = 4.0e-8_DP
    covariance(6,6) = 9.0e-8_DP
    covariance(2,1) = 1.0_DP
    covariance(1,2) = 1.0_DP

    call pod_engine_init('config/config.txt')
    call run_case(ADS_COORD_COMPONENT, 0.0_DP, failures)
    call run_case(ADS_COORD_WHITENED, 0.0_DP, failures)
    call run_case(ADS_COORD_COMPONENT, 0.02_DP, failures)
    call cleanup_gravity_network()

    if (failures /= 0) then
        write(*,'(A,I0)') 'FAIL: HFEM ADS zero-duration checks: ', failures
        stop 1
    end if
    write(*,'(A)') 'PASS: HFEM ADS preserves 6D/7D zero-duration clouds'

contains

    subroutine run_case(coordinate_mode, srp_sigma, n_fail)
        integer, intent(in) :: coordinate_mode
        real(DP), intent(in) :: srp_sigma
        integer, intent(inout) :: n_fail
        type(ads_coordinate_map_type) :: coordinate_map
        type(hfem_ads_options_type) :: options
        type(hfem_ads_stats_type) :: stats
        real(DP), allocatable :: input_samples(:,:), output_samples(:,:)
        real(DP) :: unit_points(7,4), deviation(6), eta
        character(len=256) :: message
        integer :: i, n_variables, status

        options%coordinate_mode = coordinate_mode
        options%domain_sigma = 3.0_DP
        options%srp_sigma = srp_sigma
        call ads_build_coordinate_map(nominal, covariance, coordinate_mode, &
            options%domain_sigma, srp_sigma, coordinate_map, status, message)
        call assert_true(status == 0, 'coordinate map builds', n_fail)
        n_variables = coordinate_map%n_variables
        allocate(input_samples(n_variables,4))
        unit_points = 0.0_DP
        unit_points(1:6,2) = [0.5_DP,-0.4_DP,0.3_DP,-0.2_DP,0.1_DP,-0.6_DP]
        unit_points(1:6,3) = [-0.8_DP,0.7_DP,-0.6_DP,0.5_DP,-0.4_DP,0.3_DP]
        unit_points(1:6,4) = [1.0_DP,-1.0_DP,0.0_DP,0.25_DP,-0.25_DP,0.5_DP]
        if (n_variables == 7) unit_points(7,:) = [0.0_DP,-0.5_DP,0.5_DP,1.0_DP]
        do i = 1, 4
            call ads_unit_to_physical(coordinate_map, unit_points(1:n_variables,i), &
                deviation, eta, status)
            input_samples(1:6,i) = nominal + deviation
            if (n_variables == 7) input_samples(7,i) = eta
        end do

        call hfem_ads_propagate(nominal, covariance, 0.0_DP, 0.0_DP, 0.0_DP, &
            input_samples, options, output_samples, stats, status, message)
        call assert_true(status == 0, 'zero-duration propagation succeeds', n_fail)
        call assert_true(size(output_samples,1) == n_variables, &
            'output dimension follows optional SRP coordinate', n_fail)
        call assert_true(size(output_samples,2) == 4, &
            'all four in-domain particles are retained', n_fail)
        call assert_true(maxval(abs(output_samples-input_samples)) < 1.0e-11_DP, &
            'zero-duration propagation preserves the complete cloud', n_fail)
        call assert_true(stats%input_count == 4 .and. stats%inside_count == 4 .and. &
            stats%outside_count == 0 .and. stats%propagated_count == 4 .and. &
            stats%written_count == 4, 'particle counters remain consistent', n_fail)
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

end program test_hfem_ads_zero_duration_cloud
