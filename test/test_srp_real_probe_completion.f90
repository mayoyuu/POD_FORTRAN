! A finite checkpoint before t_end must not be reported as successful propagation.
program test_srp_real_probe_completion
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config
    use pod_data_format_module, only: load_initial_opm
    use pod_uq_hfem_ads_module, only: hfem_ads_options_type
    use pod_srp_ads_history_module, only: advance_real_probes
    implicit none
    type(hfem_ads_options_type) :: options
    real(DP) :: epoch, initial(6), covariance(6,6), states(6,1)
    integer :: status
    character(len=256) :: message

    call pod_engine_init('config/config_no_srp.txt')
    config%use_earth_nspheric = .false.
    config%use_moon_nspheric = .false.
    call load_initial_opm('OPM/DRO-2/DRO-2_init.opm.json', epoch, initial, covariance)
    states(:,1) = initial
    options%max_steps = 2
    call advance_real_probes(states, [0.0_DP], epoch, 0.0_DP, 3600.0_DP, options, status, message)
    if (status == 0) error stop 'FAIL: incomplete Real integration was accepted as t_end'
    if (index(message,'end time') == 0) error stop 'FAIL: incomplete integration reason was lost'
    if (any(states(:,1) /= initial)) error stop 'FAIL: incomplete checkpoint overwrote probe state'

    options%max_steps = 10000
    call advance_real_probes(states, [0.0_DP], epoch, 0.0_DP, 60.0_DP, options, status, message)
    if (status /= 0) then
        print *, trim(message)
        error stop 'FAIL: complete Real integration was rejected'
    end if
    initial = states(:,1)
    call advance_real_probes(states, [0.0_DP], epoch, 60.0_DP, 60.0_DP, options, status, message)
    if (status /= 0 .or. any(states(:,1) /= initial)) error stop 'FAIL: zero-duration probe changed'
    print *, 'PASS: Real probes reject incomplete endpoints and accept complete intervals'
end program test_srp_real_probe_completion
