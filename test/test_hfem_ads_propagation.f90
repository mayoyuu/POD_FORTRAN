program test_hfem_ads_propagation
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config
    use pod_data_format_module, only: load_initial_opm
    use pod_uq_ads_coordinates_module, only: ADS_COORD_COMPONENT, &
        ads_coordinate_map_type, ads_build_coordinate_map, ads_unit_to_physical
    use pod_uq_hfem_ads_module, only: hfem_ads_options_type, &
        hfem_ads_stats_type, hfem_ads_propagate
    use pod_da_force_model_module, only: clear_srp_scale_uncertainty, &
        cleanup_gravity_network
    implicit none

    real(DP), parameter :: PROPAGATION_SECONDS = 604800.0_DP
    real(DP) :: epoch0, state0(6), covariance(6,6), half_width(6)
    integer :: failures

    failures = 0
    half_width = [300.0_DP,300.0_DP,300.0_DP, &
        9.0e-4_DP,9.0e-4_DP,9.0e-4_DP]
    call pod_engine_init('config/config.txt')
    call load_initial_opm('OPM/L1Halo-1/L1Halo-1_init.opm.json', &
        epoch0, state0, covariance)
    covariance = 0.0_DP
    covariance(1,1) = (half_width(1)/3.0_DP)**2
    covariance(2,2) = (half_width(2)/3.0_DP)**2
    covariance(3,3) = (half_width(3)/3.0_DP)**2
    covariance(4,4) = (half_width(4)/3.0_DP)**2
    covariance(5,5) = (half_width(5)/3.0_DP)**2
    covariance(6,6) = (half_width(6)/3.0_DP)**2
    config%use_srp = .true.

    call run_case(0.0_DP, failures)
    call run_case(0.02_DP, failures)
    call clear_srp_scale_uncertainty()
    call cleanup_gravity_network()
    if (failures /= 0) then
        write(*,'(A,I0)') 'FAIL: HFEM ADS propagation checks: ', failures
        stop 1
    end if
    write(*,'(A)') 'PASS: production HFEM ADS propagates 6D and optional 7D clouds'

contains

    subroutine run_case(srp_sigma, n_fail)
        real(DP), intent(in) :: srp_sigma
        integer, intent(inout) :: n_fail
        type(hfem_ads_options_type) :: options
        type(hfem_ads_stats_type) :: stats
        type(ads_coordinate_map_type) :: coordinate_map
        real(DP), allocatable :: input_samples(:,:), output_samples(:,:)
        real(DP) :: unit_points(7,4), deviation(6), eta, sensitivity
        character(len=256) :: message
        integer :: i, n_variables, status

        options%coordinate_mode = ADS_COORD_COMPONENT
        options%domain_sigma = 3.0_DP
        options%srp_sigma = srp_sigma
        options%da_order = 4
        options%max_split_depth = 8
        call ads_build_coordinate_map(state0, covariance, options%coordinate_mode, &
            options%domain_sigma, srp_sigma, coordinate_map, status, message)
        n_variables = coordinate_map%n_variables
        allocate(input_samples(n_variables,4))
        unit_points = 0.0_DP
        unit_points(1:6,2) = [0.25_DP,-0.25_DP,0.5_DP, &
            -0.5_DP,0.125_DP,-0.125_DP]
        if (n_variables == 7) unit_points(7,3:4) = [-1.0_DP,1.0_DP]
        do i = 1, 4
            call ads_unit_to_physical(coordinate_map, unit_points(1:n_variables,i), &
                deviation, eta, status)
            input_samples(1:6,i) = state0 + deviation
            if (n_variables == 7) input_samples(7,i) = eta
        end do

        call hfem_ads_propagate(state0, covariance, epoch0, 0.0_DP, &
            PROPAGATION_SECONDS, input_samples, options, output_samples, &
            stats, status, message)
        call assert_true(status == 0, 'production propagation status', n_fail)
        call assert_true(stats%n_variables == n_variables, &
            'DACE dimension follows optional seventh coordinate', n_fail)
        call assert_true(stats%n_patches > 0, 'ADS accepts at least one Patch', n_fail)
        call assert_true(size(output_samples,2) == 4, &
            'long propagation retains the complete input cloud', n_fail)
        call assert_true(all(ieee_is_finite(output_samples)), &
            'propagated error cloud is finite', n_fail)
        if (n_variables == 7) then
            sensitivity = maxval(abs(output_samples(1:6,4)-output_samples(1:6,3)))
            call assert_true(sensitivity > 1.0e-12_DP, &
                'global SRP scale uncertainty changes the trajectory', n_fail)
            call assert_true(maxval(abs(output_samples(7,:)-input_samples(7,:))) == 0.0_DP, &
                'SRP error coordinate is retained with each particle', n_fail)
        end if
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

end program test_hfem_ads_propagation
