!> @file test_halo_srp_scale_uncertainty_15day.f90
!! @brief Fifteen-day SRP-only scale uncertainty propagation and export.
!!
!! Four nominal force models share the same spacecraft and initial orbit:
!! cannonball, then Sun/Earth/Moon-pointing box-wing.  The only DA variable is
!! delta_s=0.1*q in a_SRP=(1+delta_s)*a_nominal.  The initial orbit is a DA
!! constant.  Each complete hourly DA map is evaluated by deterministic
!! quadrature for equal-variance uniform and Gaussian inputs.  Deterministic
!! Real trajectories close the loop at uniform bounds and Gaussian +/-3 sigma.
program test_halo_srp_scale_uncertainty_15day
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config, validate_config
    use pod_data_format_module, only: load_initial_opm
    use pod_force_model_module, only: &
        set_propagation_epoch_real => set_propagation_epoch, &
        set_srp_scale_error, clear_srp_scale_error
    use pod_da_force_model_module, only: &
        set_propagation_epoch_da => set_propagation_epoch, &
        init_gravity_network_da => init_gravity_network, &
        cleanup_gravity_network_da => cleanup_gravity_network, &
        set_srp_scale_uncertainty, clear_srp_scale_uncertainty, &
        set_srp_ballistic_parameters
    use pod_integrator_module, only: adaptive_step_integrate, METHOD_RKF78
    use pod_da_integrator_module, only: da_adaptive_step_integrate
    use pod_dace_classes, only: AlgebraicVector, CompiledDA, dace_initialize, &
                                dace_set_to, dace_get_to, active_da_count, &
                                assignment(=)
    use pod_frame_module, only: build_rtn_rotation, transform_vector_to_rtn, &
                                transform_covariance6_to_rtn
    use pod_uncertainty_diagnostics_module, only: &
        gauss_legendre_probability_rule, gauss_normal_probability_rule, &
        compute_weighted_moments, compute_covariance_axes, compute_effective_rank
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    integer, parameter :: N_MODELS=4, N_DISTRIBUTIONS=2, N_EPOCHS=361
    integer, parameter :: N_QUAD=9, N_VALIDATION=4
    integer, parameter :: DA_MAX_ORDER=8, N_DA_VARS=1, MAX_STEPS_PER_HOUR=1000
    integer, parameter :: STRATEGY_CONSTANT4=1, STRATEGY_SCHEDULED=2
    integer, parameter :: STRATEGY_CONSTANT8=3
    real(DP), parameter :: OUTPUT_STEP_S=3600.0_DP
    real(DP), parameter :: DURATION_S=360.0_DP*OUTPUT_STEP_S
    real(DP), parameter :: SCALE_SPAN=0.1_DP
    real(DP), parameter :: SIGMA_Q=1.0_DP/sqrt(3.0_DP)
    real(DP), parameter :: SKEW_LIMIT=0.1_DP, KURT_LIMIT=0.2_DP
    real(DP), parameter :: NONLINEAR_LIMIT=0.05_DP
    real(DP), parameter :: REL_TOL=1.0e-12_DP, ABS_TOL=1.0e-12_DP
    real(DP), parameter :: DT_MIN_S=1.0e-6_DP, DT_MAX_S=3600.0_DP
    real(DP), parameter :: POSITION_ABS_TOL_KM=1.0e-2_DP
    real(DP), parameter :: VELOCITY_ABS_TOL_KMS=1.0e-8_DP
    real(DP), parameter :: VALIDATION_REL_TOL=1.0e-3_DP
    real(DP), parameter :: SPACECRAFT_MASS_KG=1200.0_DP
    real(DP), parameter :: CANNONBALL_AREA_M2=9.0_DP
    real(DP), parameter :: CANNONBALL_CR=1.25_DP
    real(DP), parameter :: CANNONBALL_SMR=CANNONBALL_AREA_M2/SPACECRAFT_MASS_KG
    real(DP), parameter :: SOLAR_PRESSURE_1AU=1367.0_DP/299792458.0_DP
    character(len=*), parameter :: CONFIG_FILE='config/config.txt'
    character(len=*), parameter :: HALO_OPM='OPM/L1Halo-1/L1Halo-1_init.opm.json'
    character(len=*), parameter :: OUTPUT_DIR='SRP/260903_srp_uncertainty'
    character(len=32), parameter :: MODEL_NAMES(N_MODELS)=[character(len=32) :: &
        'cannonball','box_wing_sun','box_wing_earth','box_wing_moon']
    character(len=16), parameter :: DISTRIBUTION_NAMES(N_DISTRIBUTIONS)= &
        [character(len=16) :: 'uniform','Gaussian']
    character(len=16), parameter :: STRATEGY_NAMES(3)= &
        [character(len=16) :: 'constant_4','scheduled_4_6_8','constant_8']
    real(DP), parameter :: VALIDATION_ERRORS(N_VALIDATION)= &
        [-sqrt(3.0_DP)*SCALE_SPAN,-SCALE_SPAN,SCALE_SPAN, &
          sqrt(3.0_DP)*SCALE_SPAN]

    real(DP) :: epoch0, state0(6), covariance0(6,6), initial_nd(6)
    real(DP) :: uniform_nodes(N_QUAD), uniform_weights(N_QUAD)
    real(DP) :: gaussian_nodes(N_QUAD), gaussian_weights(N_QUAD)
    real(DP) :: da_nominal(6,N_EPOCHS)
    real(DP) :: da_validation(6,N_VALIDATION,N_EPOCHS)
    real(DP) :: real_nominal(6,N_EPOCHS), real_perturbed(6,N_EPOCHS)
    real(DP) :: first_gaussian_invalid(N_MODELS,N_DISTRIBUTIONS)
    real(DP) :: first_ellipsoid_invalid(N_MODELS,N_DISTRIBUTIONS)
    real(DP) :: max_skew(N_MODELS,N_DISTRIBUTIONS)
    real(DP) :: max_kurt(N_MODELS,N_DISTRIBUTIONS)
    real(DP) :: max_eta_nl(N_MODELS,N_DISTRIBUTIONS)
    real(DP) :: max_eta_r(N_MODELS,N_DISTRIBUTIONS)
    real(DP) :: max_eta_v(N_MODELS,N_DISTRIBUTIONS)
    real(DP) :: final_mean(6,N_MODELS,N_DISTRIBUTIONS)
    real(DP) :: final_pos_axis(3,N_MODELS,N_DISTRIBUTIONS)
    real(DP) :: final_vel_axis(3,N_MODELS,N_DISTRIBUTIONS)
    real(DP) :: order3_schedule_pos_diff(N_MODELS)
    real(DP) :: order3_schedule_vel_diff(N_MODELS)
    real(DP) :: scheduled_phase_seconds(3),constant4_phase_seconds(3)
    real(DP) :: constant8_phase_seconds(3)
    real(DP) :: constant4_final(6,N_VALIDATION)
    real(DP) :: constant8_final(6,N_VALIDATION)
    integer :: history_unit, validation_unit, summary_unit, timing_unit
    integer :: status, model_id, validation_id, epoch_id, exit_status
    logical :: all_validation_pass

    call pod_engine_init(CONFIG_FILE)
    call init_gravity_network_da()
    call load_initial_opm(HALO_OPM,epoch0,state0,covariance0)
    call set_propagation_epoch_real(epoch0)
    call set_propagation_epoch_da(epoch0)
    initial_nd(1:3)=state0(1:3)/config%LU
    initial_nd(4:6)=state0(4:6)/config%VU
    call set_srp_ballistic_parameters(Cr=CANNONBALL_CR,SMR=CANNONBALL_SMR, &
                                      RP=SOLAR_PRESSURE_1AU)
    call clear_srp_scale_error()
    call clear_srp_scale_uncertainty()
    call assert_true(da_order_for_output_hour(0)==4, &
                     '0 h must use fourth-order DA')
    call assert_true(da_order_for_output_hour(100)==4, &
                     '100 h must still use fourth-order DA')
    call assert_true(da_order_for_output_hour(101)==6, &
                     '101 h must use sixth-order DA')
    call assert_true(da_order_for_output_hour(250)==6, &
                     '250 h must still use sixth-order DA')
    call assert_true(da_order_for_output_hour(251)==8, &
                     '251 h must use eighth-order DA')
    call assert_true(da_order_for_output_hour(360)==8, &
                     '360 h must use eighth-order DA')

    call gauss_legendre_probability_rule(uniform_nodes,uniform_weights,status)
    call assert_true(status==0,'failed to build uniform quadrature')
    call gauss_normal_probability_rule(SIGMA_Q,gaussian_nodes, &
                                       gaussian_weights,status)
    call assert_true(status==0,'failed to build Gaussian quadrature')

    first_gaussian_invalid=-1.0_DP
    first_ellipsoid_invalid=-1.0_DP
    max_skew=0.0_DP
    max_kurt=0.0_DP
    max_eta_nl=0.0_DP
    max_eta_r=0.0_DP
    max_eta_v=0.0_DP
    final_mean=0.0_DP
    final_pos_axis=0.0_DP
    final_vel_axis=0.0_DP
    order3_schedule_pos_diff=0.0_DP
    order3_schedule_vel_diff=0.0_DP

    call execute_command_line('mkdir -p '//OUTPUT_DIR,exitstat=exit_status)
    call assert_true(exit_status==0,'failed to create SRP output directory')
    call open_outputs(history_unit,validation_unit,summary_unit,timing_unit)

    all_validation_pass=.true.
    write(*,'(a)') '15-day SRP-only scale uncertainty propagation'
    do model_id=1,N_MODELS
        call configure_srp_case(model_id)
        write(*,'(a,a)') '  DA model: ',trim(MODEL_NAMES(model_id))
        call propagate_da_history(model_id,history_unit,da_nominal, &
                                  da_validation,scheduled_phase_seconds)
        call compute_order3_schedule_diagnostic(model_id, &
                                        da_validation(:,:,N_EPOCHS), &
                                        order3_schedule_pos_diff(model_id), &
                                        order3_schedule_vel_diff(model_id))

        call propagate_da_timing(STRATEGY_CONSTANT4,constant4_final, &
                                 constant4_phase_seconds)
        call propagate_da_timing(STRATEGY_CONSTANT8,constant8_final, &
                                 constant8_phase_seconds)
        call write_timing_records(model_id,constant4_phase_seconds, &
            scheduled_phase_seconds,constant8_phase_seconds,constant4_final, &
            da_validation(:,:,N_EPOCHS),constant8_final,timing_unit)

        call propagate_real_history(0.0_DP,real_nominal)
        do validation_id=1,N_VALIDATION
            call propagate_real_history(VALIDATION_ERRORS(validation_id), &
                                        real_perturbed)
            do epoch_id=1,N_EPOCHS
                call write_validation_record(model_id,validation_id,epoch_id, &
                    da_nominal(:,epoch_id),da_validation(:,validation_id,epoch_id), &
                    real_nominal(:,epoch_id),real_perturbed(:,epoch_id), &
                    validation_unit,all_validation_pass)
            end do
        end do
        call assert_true(all_validation_pass, &
                         trim(MODEL_NAMES(model_id))//' DA/Real closure failed')
    end do

    call write_summary(summary_unit)
    close(history_unit)
    close(validation_unit)
    close(summary_unit)
    close(timing_unit)
    call clear_srp_scale_error()
    call clear_srp_scale_uncertainty()
    call cleanup_gravity_network_da()
    call assert_true(active_da_count()==0,'DA handles remain after full test')
    call assert_true(all_validation_pass,'one or more DA/Real records failed')
    write(*,'(a)') 'PASS: 15-day SRP uncertainty histories and validation CSVs written.'

contains

    !> Active DA order attached to an hourly output epoch.
    pure integer function da_order_for_output_hour(hour) result(order)
        integer, intent(in) :: hour

        if(hour<=100) then
            order=4
        else if(hour<=250) then
            order=6
        else
            order=8
        end if
    end function da_order_for_output_hour

    pure integer function da_order_for_strategy(strategy,hour) result(order)
        integer, intent(in) :: strategy,hour

        select case(strategy)
        case(STRATEGY_CONSTANT4)
            order=4
        case(STRATEGY_SCHEDULED)
            order=da_order_for_output_hour(hour)
        case(STRATEGY_CONSTANT8)
            order=8
        case default
            order=-1
        end select
    end function da_order_for_strategy

    pure integer function order_slot(order) result(slot)
        integer, intent(in) :: order

        select case(order)
        case(4)
            slot=1
        case(6)
            slot=2
        case(8)
            slot=3
        case default
            slot=0
        end select
    end function order_slot

    !> Configure one nominal SRP model without changing spacecraft properties.
    subroutine configure_srp_case(id)
        integer, intent(in) :: id

        config%use_srp=.true.
        config%srp_mass_kg=SPACECRAFT_MASS_KG
        config%srp_box_dimensions_m=[2.0_DP,3.0_DP,4.0_DP]
        config%srp_box_optical=[0.30_DP,0.40_DP,0.30_DP]
        config%srp_array_total_area_m2=24.0_DP
        config%srp_array_tracking_mode='single_axis'
        config%srp_array_hinge_axis_body=[0.0_DP,1.0_DP,0.0_DP]
        config%srp_array_reference_normal_body=[1.0_DP,0.0_DP,0.0_DP]
        config%srp_array_front_optical=[0.10_DP,0.80_DP,0.10_DP]
        config%srp_array_back_optical=[0.60_DP,0.20_DP,0.20_DP]
        config%srp_roll_reference='orbit_normal'
        config%srp_primary_axis_body=[0.0_DP,0.0_DP,1.0_DP]
        config%srp_secondary_axis_body=[0.0_DP,1.0_DP,0.0_DP]
        config%srp_pressure_1au_n_m2=SOLAR_PRESSURE_1AU
        config%srp_geometry_tolerance=1.0e-12_DP
        config%srp_scale_da_span=0.0_DP
        config%srp_attitude_bias_span_arcsec=0.0_DP
        config%srp_array_angle_span_deg=0.0_DP

        select case(id)
        case(1)
            config%srp_model='cannonball'
            config%srp_attitude_mode='sun'
        case(2)
            config%srp_model='box_wing'
            config%srp_attitude_mode='sun'
        case(3)
            config%srp_model='box_wing'
            config%srp_attitude_mode='earth'
        case(4)
            config%srp_model='box_wing'
            config%srp_attitude_mode='moon'
        case default
            error stop 'invalid SRP model id'
        end select
        call assert_true(validate_config(),'invalid configuration: '// &
                         trim(MODEL_NAMES(id)))
    end subroutine configure_srp_case

    !> Open all reproducible CSV outputs and emit their schemas.
    subroutine open_outputs(history,validation,summary,timing)
        integer, intent(out) :: history,validation,summary,timing
        integer :: ios

        open(newunit=history,file=OUTPUT_DIR//'/srp_uncertainty_history.csv', &
             status='replace',action='write',iostat=ios)
        call assert_true(ios==0,'cannot open uncertainty history CSV')
        open(newunit=validation,file=OUTPUT_DIR//'/srp_da_real_validation.csv', &
             status='replace',action='write',iostat=ios)
        call assert_true(ios==0,'cannot open DA/Real validation CSV')
        open(newunit=summary,file=OUTPUT_DIR//'/srp_uncertainty_summary.csv', &
             status='replace',action='write',iostat=ios)
        call assert_true(ios==0,'cannot open uncertainty summary CSV')
        open(newunit=timing,file=OUTPUT_DIR//'/srp_da_order_timing.csv', &
             status='replace',action='write',iostat=ios)
        call assert_true(ios==0,'cannot open DA order timing CSV')
        call write_history_header(history)
        call write_validation_header(validation)
        call write_summary_header(summary)
        call write_timing_header(timing)
    end subroutine open_outputs

    subroutine write_history_header(unit)
        integer, intent(in) :: unit
        write(unit,'(a)',advance='no') &
            'model,distribution,time_hours,da_order,mean_error_x_km,'// &
            'mean_error_y_km,'// &
            'mean_error_z_km,'
        write(unit,'(a)',advance='no') &
            'mean_error_vx_kms,mean_error_vy_kms,mean_error_vz_kms,'// &
            'mean_error_r_km,mean_error_t_km,'
        write(unit,'(a)',advance='no') &
            'mean_error_n_km,mean_error_vr_kms,mean_error_vt_kms,'// &
            'mean_error_vn_kms,'
        call write_covariance_header(unit,'p_i')
        call write_covariance_header(unit,'p_rtn')
        write(unit,'(a)',advance='no') &
            'pos_sigma_1_km,pos_sigma_2_km,pos_sigma_3_km,'
        write(unit,'(a)',advance='no') &
            'vel_sigma_1_kms,vel_sigma_2_kms,vel_sigma_3_kms,'
        call write_direction_header(unit,'pos_axis_i')
        call write_direction_header(unit,'pos_axis_rtn')
        call write_direction_header(unit,'vel_axis_i')
        call write_direction_header(unit,'vel_axis_rtn')
        write(unit,'(a)',advance='no') &
            'skew_x,skew_y,skew_z,skew_vx,skew_vy,skew_vz,'
        write(unit,'(a)',advance='no') &
            'kurt_x,kurt_y,kurt_z,kurt_vx,kurt_vy,kurt_vz,'
        write(unit,'(a)') 'pos_pc1_skew,pos_pc1_kurt,vel_pc1_skew,'// &
            'vel_pc1_kurt,eta_nl_pos,eta_nl_vel,eta_nl,eta_r,eta_v,'// &
            'effective_rank,gaussian_valid,ellipsoid_valid'
    end subroutine write_history_header

    subroutine write_covariance_header(unit,prefix)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: prefix
        integer :: i,j
        do i=1,6
            do j=i,6
                write(unit,'(a,"_",i0,i0,",")',advance='no') trim(prefix),i,j
            end do
        end do
    end subroutine write_covariance_header

    subroutine write_direction_header(unit,prefix)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: prefix
        integer :: axis,component
        character(len=1), parameter :: labels(3)=['x','y','z']
        do axis=1,3
            do component=1,3
                write(unit,'(a,"_",i0,"_",a,",")',advance='no') &
                    trim(prefix),axis,labels(component)
            end do
        end do
    end subroutine write_direction_header

    subroutine write_validation_header(unit)
        integer, intent(in) :: unit
        write(unit,'(a)') 'model,scale_error,time_hours,da_order,'// &
            'da_dx_km,da_dy_km,da_dz_km,da_dvx_kms,da_dvy_kms,da_dvz_kms,'// &
            'real_dx_km,real_dy_km,real_dz_km,real_dvx_kms,real_dvy_kms,'// &
            'real_dvz_kms,diff_x_km,diff_y_km,diff_z_km,diff_vx_kms,'// &
            'diff_vy_kms,diff_vz_kms,pos_diff_norm_km,vel_diff_norm_kms,'// &
            'pos_tolerance_km,vel_tolerance_kms,passed'
    end subroutine write_validation_header

    subroutine write_summary_header(unit)
        integer, intent(in) :: unit
        write(unit,'(a)') 'model,distribution,first_gaussian_invalid_hour,'// &
            'first_ellipsoid_invalid_hour,max_abs_skew,max_abs_excess_kurtosis,'// &
            'max_eta_nl,max_eta_r,max_eta_v,final_mean_error_x_km,'// &
            'final_mean_error_y_km,final_mean_error_z_km,'// &
            'final_mean_error_vx_kms,final_mean_error_vy_kms,'// &
            'final_mean_error_vz_kms,final_pos_sigma_1_km,'// &
            'final_pos_sigma_2_km,'// &
            'final_pos_sigma_3_km,final_vel_sigma_1_kms,final_vel_sigma_2_kms,'// &
            'final_vel_sigma_3_kms,order3_scheduled_pos_diff_km,'// &
            'order3_scheduled_vel_diff_kms'
    end subroutine write_summary_header

    subroutine write_timing_header(unit)
        integer, intent(in) :: unit
        write(unit,'(a)') 'model,strategy,engine_max_order,'// &
            'order4_seconds,order6_seconds,order8_seconds,total_seconds,'// &
            'speedup_vs_constant8,max_final_pos_diff_vs_constant8_km,'// &
            'max_final_vel_diff_vs_constant8_kms'
    end subroutine write_timing_header

    !> Propagate one scheduled-order, one-variable DA map on the hourly grid.
    subroutine propagate_da_history(id,unit,nominal_history,validation_history, &
                                    phase_seconds)
        integer, intent(in) :: id,unit
        real(DP), intent(out) :: nominal_history(6,N_EPOCHS)
        real(DP), intent(out) :: validation_history(6,N_VALIDATION,N_EPOCHS)
        real(DP), intent(out) :: phase_seconds(3)
        type(AlgebraicVector) :: current_state,next_state
        real(DP), allocatable :: times(:),nominal_steps(:,:)
        real(DP) :: start_time,end_time
        integer :: i,hour,n_steps,handles_before,active_order
        integer :: clock_start,clock_end,clock_rate

        handles_before=active_da_count()
        call assert_true(handles_before==0,'DA handles nonzero before model')
        call dace_initialize(DA_MAX_ORDER,N_DA_VARS)
        call dace_set_to(4)
        call set_srp_scale_uncertainty(1,0.0_DP,SCALE_SPAN)
        phase_seconds=0.0_DP
        call current_state%init(6)
        do i=1,6
            current_state%elements(i)=initial_nd(i)
        end do
        call process_da_epoch(id,1,0.0_DP,4,current_state,unit, &
                              nominal_history,validation_history)

        do hour=1,N_EPOCHS-1
            active_order=da_order_for_output_hour(hour)
            call dace_set_to(active_order)
            call assert_true(dace_get_to()==active_order, &
                             'DACE truncation order was not applied')
            start_time=real(hour-1,DP)*OUTPUT_STEP_S/config%TU
            end_time=real(hour,DP)*OUTPUT_STEP_S/config%TU
            call system_clock(clock_start,clock_rate)
            call da_adaptive_step_integrate(current_state,start_time,end_time, &
                METHOD_RKF78,times,nominal_steps,next_state,n_steps, &
                max_steps_in=MAX_STEPS_PER_HOUR,rel_tol_in=REL_TOL, &
                abs_tol_in=ABS_TOL,dt_min_in=DT_MIN_S,dt_max_in=DT_MAX_S)
            call system_clock(clock_end)
            phase_seconds(order_slot(active_order))= &
                phase_seconds(order_slot(active_order))+ &
                real(clock_end-clock_start,DP)/real(clock_rate,DP)
            call assert_true(n_steps>1,'DA hourly segment produced no step')
            call process_da_epoch(id,hour+1,real(hour,DP),active_order, &
                                  next_state,unit,nominal_history, &
                                  validation_history)
            current_state=next_state
            call next_state%destroy()
            if(allocated(times)) deallocate(times)
            if(allocated(nominal_steps)) deallocate(nominal_steps)
        end do

        call current_state%destroy()
        call next_state%destroy()
        call clear_srp_scale_uncertainty()
        call assert_true(active_da_count()==handles_before, &
                         'DA handle leak in hourly model propagation')
    end subroutine propagate_da_history

    !> Compile one hourly state map once, then evaluate every requested point.
    subroutine process_da_epoch(id,index,time_hours,active_order,state_da,unit, &
                                nominal_history,validation_history)
        integer, intent(in) :: id,index,active_order,unit
        real(DP), intent(in) :: time_hours
        type(AlgebraicVector), intent(in) :: state_da
        real(DP), intent(inout) :: nominal_history(6,N_EPOCHS)
        real(DP), intent(inout) :: validation_history(6,N_VALIDATION,N_EPOCHS)
        type(CompiledDA) :: compiled
        real(DP), allocatable :: evaluated_nd(:)
        real(DP) :: nominal_nd(6),nominal_physical(6),derivative_physical(6)
        real(DP) :: q
        integer :: i

        nominal_nd=state_da%cons()
        call nondimensional_to_physical(nominal_nd,nominal_physical)
        nominal_history(:,index)=nominal_physical
        do i=1,6
            derivative_physical(i)=state_da%elements(i)%get_deriv_value(1)
        end do
        derivative_physical(1:3)=derivative_physical(1:3)*config%LU
        derivative_physical(4:6)=derivative_physical(4:6)*config%VU

        compiled=state_da%compile()
        do i=1,N_VALIDATION
            q=VALIDATION_ERRORS(i)/SCALE_SPAN
            evaluated_nd=compiled%eval([q])
            call nondimensional_to_physical(evaluated_nd, &
                                            validation_history(:,i,index))
        end do
        call analyze_distribution(id,1,index,time_hours,active_order,compiled, &
            nominal_physical,derivative_physical,uniform_nodes,uniform_weights,unit)
        call analyze_distribution(id,2,index,time_hours,active_order,compiled, &
            nominal_physical,derivative_physical,gaussian_nodes,gaussian_weights,unit)
        call compiled%destroy()
        if(allocated(evaluated_nd)) deallocate(evaluated_nd)
    end subroutine process_da_epoch

    !> Evaluate one distribution and write all covariance/shape diagnostics.
    subroutine analyze_distribution(id,distribution_id,index,time_hours, &
        active_order,compiled,nominal,derivative,nodes,weights,unit)
        integer, intent(in) :: id,distribution_id,index,active_order,unit
        real(DP), intent(in) :: time_hours,nominal(6),derivative(6)
        real(DP), intent(in) :: nodes(N_QUAD),weights(N_QUAD)
        type(CompiledDA), intent(in) :: compiled

        real(DP) :: samples(6,N_QUAD),linear_samples(6,N_QUAD)
        real(DP) :: mean_value(6),covariance(6,6),covariance_rtn(6,6)
        real(DP) :: position_covariance(3,3),velocity_covariance(3,3)
        real(DP) :: skewness(6),kurtosis(6),c_rtn_i(3,3)
        real(DP) :: mean_rtn(6),pos_eigenvalues(3),vel_eigenvalues(3)
        real(DP) :: pos_vectors_i(3,3),vel_vectors_i(3,3)
        real(DP) :: pos_vectors_rtn(3,3),vel_vectors_rtn(3,3)
        real(DP) :: pos_axis(3),vel_axis(3),evaluated(6)
        real(DP) :: projected(1,N_QUAD),projected_mean(1)
        real(DP) :: projected_cov(1,1),projected_skew(1),projected_kurt(1)
        real(DP) :: pos_pc_skew,pos_pc_kurt,vel_pc_skew,vel_pc_kurt
        real(DP) :: eta_pos,eta_vel,eta_nl,eta_r,eta_v
        real(DP) :: nonlinear_pos,nonlinear_vel,pos_scale,vel_scale
        real(DP) :: max_skew_here,max_kurt_here
        real(DP) :: covariance_upper(21),covariance_rtn_upper(21)
        integer :: j,status_local,rank_value
        logical :: gaussian_valid,ellipsoid_valid
        character(len=5) :: gaussian_text,ellipsoid_text

        do j=1,N_QUAD
            evaluated=compiled%eval([nodes(j)])
            call nondimensional_to_physical(evaluated,samples(:,j))
            samples(:,j)=samples(:,j)-nominal
            linear_samples(:,j)=derivative*nodes(j)
        end do
        call compute_weighted_moments(samples,weights,mean_value,covariance, &
                                      skewness,kurtosis,status_local)
        call assert_true(status_local==0,'weighted moment failure')
        call assert_true(all(ieee_is_finite(mean_value)) .and. &
                         all(ieee_is_finite(covariance)), &
                         'non-finite uncertainty moments')
        call assert_true(maxval(abs(covariance-transpose(covariance)))<=1.0e-10_DP* &
                         max(1.0_DP,maxval(abs(covariance))), &
                         'asymmetric uncertainty covariance')

        position_covariance=covariance(1:3,1:3)
        velocity_covariance=covariance(4:6,4:6)
        call compute_covariance_axes(position_covariance,pos_eigenvalues, &
                                     pos_vectors_i,pos_axis,status_local)
        call assert_true(status_local==0,'position covariance axes failed')
        call compute_covariance_axes(velocity_covariance,vel_eigenvalues, &
                                     vel_vectors_i,vel_axis,status_local)
        call assert_true(status_local==0,'velocity covariance axes failed')
        call compute_effective_rank(covariance,1.0e-12_DP,rank_value,status_local)
        call assert_true(status_local==0,'effective covariance rank failed')

        call build_rtn_rotation(nominal(1:3),nominal(4:6),c_rtn_i, &
                                status_local)
        call assert_true(status_local==0,'nominal state has undefined RTN frame')
        call transform_vector_to_rtn(mean_value(1:3),c_rtn_i,mean_rtn(1:3))
        call transform_vector_to_rtn(mean_value(4:6),c_rtn_i,mean_rtn(4:6))
        call transform_covariance6_to_rtn(covariance,c_rtn_i,covariance_rtn)
        pos_vectors_rtn=matmul(c_rtn_i,pos_vectors_i)
        vel_vectors_rtn=matmul(c_rtn_i,vel_vectors_i)

        do j=1,N_QUAD
            projected(1,j)=dot_product(samples(1:3,j)-mean_value(1:3), &
                                       pos_vectors_i(:,1))
        end do
        call compute_weighted_moments(projected,weights,projected_mean, &
            projected_cov,projected_skew,projected_kurt,status_local)
        call assert_true(status_local==0,'position principal moment failed')
        pos_pc_skew=projected_skew(1)
        pos_pc_kurt=projected_kurt(1)
        do j=1,N_QUAD
            projected(1,j)=dot_product(samples(4:6,j)-mean_value(4:6), &
                                       vel_vectors_i(:,1))
        end do
        call compute_weighted_moments(projected,weights,projected_mean, &
            projected_cov,projected_skew,projected_kurt,status_local)
        call assert_true(status_local==0,'velocity principal moment failed')
        vel_pc_skew=projected_skew(1)
        vel_pc_kurt=projected_kurt(1)

        nonlinear_pos=0.0_DP
        nonlinear_vel=0.0_DP
        do j=1,N_QUAD
            nonlinear_pos=nonlinear_pos+weights(j)* &
                sum((samples(1:3,j)-linear_samples(1:3,j))**2)
            nonlinear_vel=nonlinear_vel+weights(j)* &
                sum((samples(4:6,j)-linear_samples(4:6,j))**2)
        end do
        pos_scale=sqrt(max(sum(pos_eigenvalues),0.0_DP))
        vel_scale=sqrt(max(sum(vel_eigenvalues),0.0_DP))
        if(pos_scale>tiny(1.0_DP)) then
            eta_pos=sqrt(max(nonlinear_pos,0.0_DP))/pos_scale
        else
            eta_pos=0.0_DP
        end if
        if(vel_scale>tiny(1.0_DP)) then
            eta_vel=sqrt(max(nonlinear_vel,0.0_DP))/vel_scale
        else
            eta_vel=0.0_DP
        end if
        eta_nl=max(eta_pos,eta_vel)
        if(pos_eigenvalues(1)>tiny(1.0_DP)) then
            eta_r=sqrt(max(pos_eigenvalues(2)+pos_eigenvalues(3),0.0_DP)/ &
                       pos_eigenvalues(1))
        else
            eta_r=0.0_DP
        end if
        if(vel_eigenvalues(1)>tiny(1.0_DP)) then
            eta_v=sqrt(max(vel_eigenvalues(2)+vel_eigenvalues(3),0.0_DP)/ &
                       vel_eigenvalues(1))
        else
            eta_v=0.0_DP
        end if

        max_skew_here=max(maxval(abs(skewness)),abs(pos_pc_skew), &
                          abs(vel_pc_skew))
        max_kurt_here=max(maxval(abs(kurtosis)),abs(pos_pc_kurt), &
                          abs(vel_pc_kurt))
        gaussian_valid=max_skew_here<=SKEW_LIMIT .and. &
                       max_kurt_here<=KURT_LIMIT
        ellipsoid_valid=eta_nl<=NONLINEAR_LIMIT .and. &
                        eta_r<=NONLINEAR_LIMIT .and. &
                        eta_v<=NONLINEAR_LIMIT
        if(distribution_id==1) then
            gaussian_text='NA'
        else if(gaussian_valid) then
            gaussian_text='true'
        else
            gaussian_text='false'
            if(first_gaussian_invalid(id,distribution_id)<0.0_DP) &
                first_gaussian_invalid(id,distribution_id)=time_hours
        end if
        if(ellipsoid_valid) then
            ellipsoid_text='true'
        else
            ellipsoid_text='false'
            if(first_ellipsoid_invalid(id,distribution_id)<0.0_DP) &
                first_ellipsoid_invalid(id,distribution_id)=time_hours
        end if

        max_skew(id,distribution_id)=max(max_skew(id,distribution_id), &
                                         max_skew_here)
        max_kurt(id,distribution_id)=max(max_kurt(id,distribution_id), &
                                         max_kurt_here)
        max_eta_nl(id,distribution_id)=max(max_eta_nl(id,distribution_id),eta_nl)
        max_eta_r(id,distribution_id)=max(max_eta_r(id,distribution_id),eta_r)
        max_eta_v(id,distribution_id)=max(max_eta_v(id,distribution_id),eta_v)
        if(index==N_EPOCHS) then
            final_mean(:,id,distribution_id)=mean_value
            final_pos_axis(:,id,distribution_id)=pos_axis
            final_vel_axis(:,id,distribution_id)=vel_axis
        end if

        covariance_upper=pack_upper6(covariance)
        covariance_rtn_upper=pack_upper6(covariance_rtn)
        write(unit,'(*(g0,:,","))') trim(MODEL_NAMES(id)), &
            trim(DISTRIBUTION_NAMES(distribution_id)),time_hours,active_order, &
            mean_value, &
            mean_rtn,covariance_upper,covariance_rtn_upper,pos_axis,vel_axis, &
            pos_vectors_i,pos_vectors_rtn,vel_vectors_i,vel_vectors_rtn, &
            skewness,kurtosis,pos_pc_skew,pos_pc_kurt,vel_pc_skew,vel_pc_kurt, &
            eta_pos,eta_vel,eta_nl,eta_r,eta_v,rank_value, &
            trim(gaussian_text),trim(ellipsoid_text)
    end subroutine analyze_distribution

    pure function pack_upper6(matrix) result(packed)
        real(DP), intent(in) :: matrix(6,6)
        real(DP) :: packed(21)
        integer :: i,j,k
        k=0
        do i=1,6
            do j=i,6
                k=k+1
                packed(k)=matrix(i,j)
            end do
        end do
    end function pack_upper6

    !> Propagate one deterministic Real trajectory on the hourly grid.
    subroutine propagate_real_history(scale_error,history)
        real(DP), intent(in) :: scale_error
        real(DP), intent(out) :: history(6,N_EPOCHS)
        real(DP) :: current_state(6),physical_state(6)
        real(DP) :: start_time,end_time
        real(DP), allocatable :: times(:),states(:,:)
        integer :: hour,n_steps

        call set_srp_scale_error(scale_error)
        current_state=initial_nd
        call nondimensional_to_physical(current_state,history(:,1))
        do hour=1,N_EPOCHS-1
            start_time=real(hour-1,DP)*OUTPUT_STEP_S/config%TU
            end_time=real(hour,DP)*OUTPUT_STEP_S/config%TU
            call adaptive_step_integrate(current_state,start_time,end_time, &
                METHOD_RKF78,times,states,n_steps, &
                max_steps_in=MAX_STEPS_PER_HOUR,rel_tol_in=REL_TOL, &
                abs_tol_in=ABS_TOL,dt_min_in=DT_MIN_S,dt_max_in=DT_MAX_S)
            call assert_true(n_steps>1,'Real hourly segment produced no step')
            current_state=states(n_steps,:)
            call nondimensional_to_physical(current_state,physical_state)
            history(:,hour+1)=physical_state
            if(allocated(times)) deallocate(times)
            if(allocated(states)) deallocate(states)
        end do
        call clear_srp_scale_error()
    end subroutine propagate_real_history

    !> Compare a saved scheduled-order DA evaluation with one Real trajectory.
    subroutine write_validation_record(id,value_id,index,nominal_da,perturbed_da, &
        nominal_real,perturbed_real,unit,all_pass)
        integer, intent(in) :: id,value_id,index,unit
        real(DP), intent(in) :: nominal_da(6),perturbed_da(6)
        real(DP), intent(in) :: nominal_real(6),perturbed_real(6)
        logical, intent(inout) :: all_pass
        real(DP) :: da_error(6),real_error(6),difference(6)
        real(DP) :: pos_diff,vel_diff,pos_tolerance,vel_tolerance
        logical :: passed
        character(len=5) :: passed_text

        da_error=perturbed_da-nominal_da
        real_error=perturbed_real-nominal_real
        difference=da_error-real_error
        pos_diff=sqrt(sum(difference(1:3)**2))
        vel_diff=sqrt(sum(difference(4:6)**2))
        pos_tolerance=max(POSITION_ABS_TOL_KM,VALIDATION_REL_TOL* &
                          sqrt(sum(real_error(1:3)**2)))
        vel_tolerance=max(VELOCITY_ABS_TOL_KMS,VALIDATION_REL_TOL* &
                          sqrt(sum(real_error(4:6)**2)))
        passed=pos_diff<=pos_tolerance .and. vel_diff<=vel_tolerance
        if(passed) then
            passed_text='true'
        else
            passed_text='false'
            all_pass=.false.
        end if
        write(unit,'(*(g0,:,","))') trim(MODEL_NAMES(id)), &
            VALIDATION_ERRORS(value_id),real(index-1,DP), &
            da_order_for_output_hour(index-1),da_error,real_error, &
            difference,pos_diff,vel_diff,pos_tolerance,vel_tolerance, &
            trim(passed_text)
    end subroutine write_validation_record

    !> Compare a one-shot third-order endpoint with the scheduled endpoint.
    subroutine compute_order3_schedule_diagnostic(id,scheduled_values, &
                                                   pos_diff,vel_diff)
        integer, intent(in) :: id
        real(DP), intent(in) :: scheduled_values(6,N_VALIDATION)
        real(DP), intent(out) :: pos_diff,vel_diff
        type(AlgebraicVector) :: initial_state,final_state
        type(CompiledDA) :: compiled
        real(DP), allocatable :: times(:),nominal_steps(:,:),evaluated(:)
        real(DP) :: order3_value(6),q
        integer :: i,j,n_steps,handles_before

        if(id<1) continue
        handles_before=active_da_count()
        call assert_true(handles_before==0,'handles before order convergence run')
        call dace_initialize(3,N_DA_VARS)
        call set_srp_scale_uncertainty(1,0.0_DP,SCALE_SPAN)
        call initial_state%init(6)
        do i=1,6
            initial_state%elements(i)=initial_nd(i)
        end do
        call da_adaptive_step_integrate(initial_state,0.0_DP,DURATION_S/config%TU, &
            METHOD_RKF78,times,nominal_steps,final_state,n_steps, &
            max_steps_in=10000,rel_tol_in=REL_TOL,abs_tol_in=ABS_TOL, &
            dt_min_in=DT_MIN_S,dt_max_in=DT_MAX_S)
        call assert_true(n_steps>1,'third-order endpoint propagation failed')
        compiled=final_state%compile()
        pos_diff=0.0_DP
        vel_diff=0.0_DP
        do j=1,N_VALIDATION
            q=VALIDATION_ERRORS(j)/SCALE_SPAN
            evaluated=compiled%eval([q])
            call nondimensional_to_physical(evaluated,order3_value)
            pos_diff=max(pos_diff,sqrt(sum((order3_value(1:3)- &
                                           scheduled_values(1:3,j))**2)))
            vel_diff=max(vel_diff,sqrt(sum((order3_value(4:6)- &
                                           scheduled_values(4:6,j))**2)))
        end do
        call compiled%destroy()
        call initial_state%destroy()
        call final_state%destroy()
        call clear_srp_scale_uncertainty()
        if(allocated(times)) deallocate(times)
        if(allocated(nominal_steps)) deallocate(nominal_steps)
        if(allocated(evaluated)) deallocate(evaluated)
        call assert_true(active_da_count()==handles_before, &
                         'DA handles leaked in order convergence run')
    end subroutine compute_order3_schedule_diagnostic

    !> Time only the hourly DA integrations for one fixed order strategy.
    subroutine propagate_da_timing(strategy,final_values,phase_seconds)
        integer, intent(in) :: strategy
        real(DP), intent(out) :: final_values(6,N_VALIDATION)
        real(DP), intent(out) :: phase_seconds(3)
        type(AlgebraicVector) :: current_state,next_state
        type(CompiledDA) :: compiled
        real(DP), allocatable :: times(:),nominal_steps(:,:),evaluated(:)
        real(DP) :: start_time,end_time,q
        integer :: i,j,hour,n_steps,handles_before,active_order,slot
        integer :: clock_start,clock_end,clock_rate

        handles_before=active_da_count()
        call assert_true(handles_before==0, &
                         'DA handles nonzero before timing trajectory')
        call assert_true(strategy>=STRATEGY_CONSTANT4 .and. &
                         strategy<=STRATEGY_CONSTANT8, &
                         'unknown DA timing strategy')
        call dace_initialize(DA_MAX_ORDER,N_DA_VARS)
        active_order=da_order_for_strategy(strategy,0)
        call dace_set_to(active_order)
        call set_srp_scale_uncertainty(1,0.0_DP,SCALE_SPAN)
        phase_seconds=0.0_DP
        call current_state%init(6)
        do i=1,6
            current_state%elements(i)=initial_nd(i)
        end do

        do hour=1,N_EPOCHS-1
            active_order=da_order_for_strategy(strategy,hour)
            slot=order_slot(active_order)
            call assert_true(slot>0,'invalid active DA order in timing run')
            call dace_set_to(active_order)
            call assert_true(dace_get_to()==active_order, &
                             'timing run failed to set DACE order')
            start_time=real(hour-1,DP)*OUTPUT_STEP_S/config%TU
            end_time=real(hour,DP)*OUTPUT_STEP_S/config%TU
            call system_clock(clock_start,clock_rate)
            call da_adaptive_step_integrate(current_state,start_time,end_time, &
                METHOD_RKF78,times,nominal_steps,next_state,n_steps, &
                max_steps_in=MAX_STEPS_PER_HOUR,rel_tol_in=REL_TOL, &
                abs_tol_in=ABS_TOL,dt_min_in=DT_MIN_S,dt_max_in=DT_MAX_S)
            call system_clock(clock_end)
            phase_seconds(slot)=phase_seconds(slot)+ &
                real(clock_end-clock_start,DP)/real(clock_rate,DP)
            call assert_true(n_steps>1,'timing DA hourly segment failed')
            current_state=next_state
            call next_state%destroy()
            if(allocated(times)) deallocate(times)
            if(allocated(nominal_steps)) deallocate(nominal_steps)
        end do

        compiled=current_state%compile()
        do j=1,N_VALIDATION
            q=VALIDATION_ERRORS(j)/SCALE_SPAN
            evaluated=compiled%eval([q])
            call nondimensional_to_physical(evaluated,final_values(:,j))
        end do
        call compiled%destroy()
        call current_state%destroy()
        call next_state%destroy()
        call clear_srp_scale_uncertainty()
        if(allocated(evaluated)) deallocate(evaluated)
        call assert_true(active_da_count()==handles_before, &
                         'DA handles leaked in timing trajectory')
    end subroutine propagate_da_timing

    !> Emit constant-4, scheduled and constant-8 timing comparisons.
    subroutine write_timing_records(id,phase4,phase_scheduled,phase8, &
                                    final4,final_scheduled,final8,unit)
        integer, intent(in) :: id,unit
        real(DP), intent(in) :: phase4(3),phase_scheduled(3),phase8(3)
        real(DP), intent(in) :: final4(6,N_VALIDATION)
        real(DP), intent(in) :: final_scheduled(6,N_VALIDATION)
        real(DP), intent(in) :: final8(6,N_VALIDATION)
        real(DP) :: full8_seconds

        full8_seconds=sum(phase8)
        call assert_true(full8_seconds>0.0_DP .and. &
                         ieee_is_finite(full8_seconds), &
                         'invalid full eighth-order timing')
        call write_one_timing_record(id,STRATEGY_CONSTANT4,phase4, &
                                     final4,final8,full8_seconds,unit)
        call write_one_timing_record(id,STRATEGY_SCHEDULED,phase_scheduled, &
                                     final_scheduled,final8,full8_seconds,unit)
        call write_one_timing_record(id,STRATEGY_CONSTANT8,phase8, &
                                     final8,final8,full8_seconds,unit)
    end subroutine write_timing_records

    subroutine write_one_timing_record(id,strategy,phase_seconds,values, &
                                       reference8,full8_seconds,unit)
        integer, intent(in) :: id,strategy,unit
        real(DP), intent(in) :: phase_seconds(3)
        real(DP), intent(in) :: values(6,N_VALIDATION)
        real(DP), intent(in) :: reference8(6,N_VALIDATION),full8_seconds
        real(DP) :: total_seconds,speedup,pos_diff,vel_diff
        integer :: j

        total_seconds=sum(phase_seconds)
        call assert_true(total_seconds>0.0_DP .and. &
                         ieee_is_finite(total_seconds), &
                         'invalid DA order timing total')
        speedup=full8_seconds/total_seconds
        pos_diff=0.0_DP
        vel_diff=0.0_DP
        do j=1,N_VALIDATION
            pos_diff=max(pos_diff,sqrt(sum((values(1:3,j)- &
                                            reference8(1:3,j))**2)))
            vel_diff=max(vel_diff,sqrt(sum((values(4:6,j)- &
                                            reference8(4:6,j))**2)))
        end do
        write(unit,'(*(g0,:,","))') trim(MODEL_NAMES(id)), &
            trim(STRATEGY_NAMES(strategy)),DA_MAX_ORDER,phase_seconds, &
            total_seconds,speedup,pos_diff,vel_diff
        write(*,'(2x,a,1x,a,": ",f8.3," s, speedup ",f6.3,"x")') &
            trim(MODEL_NAMES(id)),trim(STRATEGY_NAMES(strategy)), &
            total_seconds,speedup
    end subroutine write_one_timing_record

    !> Write one compact record for each model/distribution pair.
    subroutine write_summary(unit)
        integer, intent(in) :: unit
        integer :: id,distribution_id
        character(len=32) :: gaussian_time,ellipsoid_time

        do id=1,N_MODELS
            do distribution_id=1,N_DISTRIBUTIONS
                if(distribution_id==1) then
                    gaussian_time='NA'
                else
                    gaussian_time=first_time_text( &
                        first_gaussian_invalid(id,distribution_id))
                end if
                ellipsoid_time=first_time_text( &
                    first_ellipsoid_invalid(id,distribution_id))
                write(unit,'(*(g0,:,","))') trim(MODEL_NAMES(id)), &
                    trim(DISTRIBUTION_NAMES(distribution_id)), &
                    trim(gaussian_time),trim(ellipsoid_time), &
                    max_skew(id,distribution_id), &
                    max_kurt(id,distribution_id), &
                    max_eta_nl(id,distribution_id), &
                    max_eta_r(id,distribution_id), &
                    max_eta_v(id,distribution_id), &
                    final_mean(:,id,distribution_id), &
                    final_pos_axis(:,id,distribution_id), &
                    final_vel_axis(:,id,distribution_id), &
                    order3_schedule_pos_diff(id), &
                    order3_schedule_vel_diff(id)
            end do
        end do
    end subroutine write_summary

    !> Render an hourly threshold crossing, preserving the no-crossing state.
    function first_time_text(time_hours) result(value)
        real(DP), intent(in) :: time_hours
        character(len=32) :: value

        if(time_hours<0.0_DP) then
            value='not_exceeded'
        else
            write(value,'(f0.1)') time_hours
        end if
    end function first_time_text

    !> Convert the internal LU/TU state to km and km/s for all CSV products.
    subroutine nondimensional_to_physical(state_nd,state_physical)
        real(DP), intent(in) :: state_nd(:)
        real(DP), intent(out) :: state_physical(6)

        call assert_true(size(state_nd)==6, &
                         'state conversion requires six elements')
        state_physical(1:3)=state_nd(1:3)*config%LU
        state_physical(4:6)=state_nd(4:6)*config%VU
    end subroutine nondimensional_to_physical

    !> Local assertion helper with a useful failing-test message.
    subroutine assert_true(condition,message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if(.not.condition) then
            write(*,'(a)') 'FAIL: '//trim(message)
            error stop 1
        end if
    end subroutine assert_true

end program test_halo_srp_scale_uncertainty_15day
