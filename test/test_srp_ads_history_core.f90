program test_srp_ads_history_core
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_data_format_module, only: load_initial_opm
    use pod_uq_ads_coordinates_module, only: ADS_COORD_WHITENED, &
        ads_unit_to_physical
    use pod_uq_hfem_ads_module, only: hfem_ads_options_type, &
        hfem_ads_stats_type, hfem_ads_history_type, hfem_ads_propagate, &
        hfem_ads_history_init, hfem_ads_history_advance, &
        hfem_ads_history_evaluate, hfem_ads_history_destroy
    use pod_da_force_model_module, only: cleanup_gravity_network
    use pod_srp_ads_snapshot_module, only: write_ads_patch_snapshot
    use pod_srp_ads_history_module, only: advance_real_probes
    implicit none

    type(hfem_ads_options_type) :: options
    type(hfem_ads_stats_type) :: stats
    type(hfem_ads_history_type) :: history
    real(DP) :: epoch0, state0(6), covariance(6,6), unit_points(7,2)
    real(DP) :: values0(6,2), values1(6,2), values2(6,2), deviation(6), eta
    real(DP) :: physical(7,2), real_states(6,2), rebuilt(6), coefficient, monomial
    real(DP), allocatable :: one_shot(:,:)
    logical :: found(2)
    integer :: status, i, io, file_unit, patch_id, component, powers(7), k
    character(len=256) :: message, line

    call pod_engine_init('config/config.txt')
    call load_initial_opm('input/DROb_20251210_9.opm', epoch0, state0, covariance)
    options%coordinate_mode = ADS_COORD_WHITENED
    options%srp_sigma = 0.02_DP
    options%max_split_depth = 2
    options%error_tolerance = [1.0e3_DP,1.0e3_DP,1.0e3_DP, &
        1.0_DP,1.0_DP,1.0_DP]
    unit_points = 0.0_DP
    unit_points(:,2) = [0.3_DP,-0.2_DP,0.1_DP,-0.2_DP, &
        0.1_DP,-0.3_DP,0.5_DP]

    call hfem_ads_history_init(history, state0, covariance, epoch0, options, status, message)
    call require(status == 0, 'history initialization: '//trim(message))
    call require(maxval(abs(matmul(history%coordinate_map%basis6, &
        transpose(history%coordinate_map%basis6))-9.0_DP*covariance)) / &
        max(1.0_DP,maxval(abs(covariance))) < 1.0e-11_DP, &
        'whitened 3-sigma basis reproduces the full covariance')
    call hfem_ads_history_evaluate(history, unit_points, values0, found, status)
    call require(status == 0 .and. all(found), 'initial domain evaluation')
    do i = 1, 2
        call ads_unit_to_physical(history%coordinate_map, unit_points(:,i), &
            deviation, eta, status)
        call require(status == 0, 'unit-to-physical mapping')
        physical(1:6,i) = state0 + deviation
        physical(7,i) = eta
    end do
    call require(abs(physical(7,2)-0.03_DP) < 1.0e-14_DP, &
        'seventh coordinate maps to 3-sigma relative SRP error')
    call require(maxval(abs(values0-physical(1:6,:))) < 1.0e-9_DP, &
        'initial polynomial is the OPM uncertainty map')

    call hfem_ads_history_advance(history, 60.0_DP, status, message)
    call require(status == 0, 'first hourly step: '//trim(message))
    call hfem_ads_history_evaluate(history, unit_points, values1, found, status)
    call require(status == 0 .and. all(found), 'first stepped evaluation')
    call hfem_ads_history_advance(history, 120.0_DP, status, message)
    call require(status == 0, 'second hourly step: '//trim(message))
    call hfem_ads_history_evaluate(history, unit_points, values2, found, status)
    call require(status == 0 .and. all(found), 'second stepped evaluation')
    call require(history%current_time == 120.0_DP, 'history time advances')
    call require(history%stats%n_patches > 0, 'history keeps accepted patches')
    real_states=physical(1:6,:)
    call advance_real_probes(real_states,physical(7,:),epoch0,0.0_DP,120.0_DP, &
        options,status,message)
    call require(status == 0, 'independent real SRP integration: '//trim(message))
    call require(maxval(abs(real_states(1:3,:)-values2(1:3,:))) < 1.0e-4_DP, &
        'short-time real and ADS positions agree')
    call require(maxval(abs(real_states(4:6,:)-values2(4:6,:))) < 1.0e-8_DP, &
        'short-time real and ADS velocities agree')
    call require(history%previous_available .and. history%previous_time == 60.0_DP, &
        'previous checkpoint preserved')
    call require(history%previous_domain%n_patches == 1, 'test uses one previous patch')
    call write_ads_patch_snapshot(history,'/tmp/pod_srp_ads_previous',.true.,status,message)
    call require(status == 0, 'export previous checkpoint snapshot: '//trim(message))
    rebuilt=0.0_DP
    open(newunit=file_unit,file='/tmp/pod_srp_ads_previous_coeff.csv', &
        status='old',action='read',iostat=io)
    call require(io == 0, 'open previous checkpoint coefficients')
    read(file_unit,'(A)',iostat=io) line
    do
        read(file_unit,'(A)',iostat=io) line
        if(io /= 0) exit
        read(line,*,iostat=io) patch_id,component,powers,coefficient
        call require(io == 0, 'parse previous checkpoint coefficient')
        monomial=coefficient
        do k=1,7
            monomial=monomial*unit_points(k,2)**powers(k)
        end do
        rebuilt(component)=rebuilt(component)+monomial
    end do
    close(file_unit)
    call require(maxval(abs(rebuilt-values1(:,2))) < 1.0e-8_DP, &
        'previous checkpoint snapshot reconstructs ADS state')
    call hfem_ads_history_destroy(history)

    call hfem_ads_propagate(state0, covariance, epoch0, 0.0_DP, 120.0_DP, &
        physical, options, one_shot, stats, status, message)
    call require(status == 0, 'one-shot comparison: '//trim(message))
    call require(size(one_shot,2) == 2, 'one-shot retains both samples')
    call require(maxval(abs(values2-one_shot(1:6,:))) < 1.0e-4_DP, &
        'incremental and one-shot ADS agree')
    call cleanup_gravity_network()
    write(*,'(A)') 'PASS: SRP ADS history preserves the 7D domain across steps'

contains
    subroutine require(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description
        if (condition) return
        write(*,'(A)') 'FAIL: '//trim(description)
        stop 1
    end subroutine require
end program test_srp_ads_history_core
