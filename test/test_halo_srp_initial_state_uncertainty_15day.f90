!> @file test_halo_srp_initial_state_uncertainty_15day.f90
!! @brief Fifteen-day initial-state Gaussian uncertainty comparison of SRP models.
!!
!! The only uncertain quantities are the six J2000 initial-state components:
!! independent one-sigma values of 10 km and 0.03 m/s.  Cannonball and
!! Sun/Earth/Moon-pointing box-wing models use deterministic spacecraft/SRP
!! properties.  One shared antithetic Gaussian ensemble is evaluated through
!! scheduled fourth/sixth/eighth-order DA maps every hour.  Twelve axis-aligned
!! +/-3-sigma Real trajectories provide an independent closure check.
program test_halo_srp_initial_state_uncertainty_15day
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config,validate_config
    use pod_data_format_module, only: load_initial_opm
    use pod_force_model_module, only: &
        set_propagation_epoch_real=>set_propagation_epoch, &
        clear_srp_scale_error
    use pod_da_force_model_module, only: &
        set_propagation_epoch_da=>set_propagation_epoch, &
        init_gravity_network_da=>init_gravity_network, &
        cleanup_gravity_network_da=>cleanup_gravity_network, &
        clear_srp_scale_uncertainty,set_srp_ballistic_parameters
    use pod_integrator_module, only: adaptive_step_integrate,METHOD_RKF78
    use pod_da_integrator_module, only: da_adaptive_step_integrate
    use pod_dace_classes
    use pod_frame_module, only: build_rtn_rotation,transform_covariance6_to_rtn
    use pod_uncertainty_diagnostics_module, only: &
        generate_antithetic_standard_normal_samples,compute_weighted_moments, &
        compute_covariance_axes,compute_effective_rank
    use pod_measurement_base_module, only: observation_station
    use pod_obs_io_module, only: station_record,preload_stations,find_station_by_id
    use pod_measurement_model_module, only: compute_measurement
    use pod_measurement_da_module, only: compute_measurement_da
    use pod_basicmath_module, only: wrap_angle_rad,eigenvalue_decomposition
    use,intrinsic :: ieee_arithmetic,only: ieee_is_finite
    implicit none

    integer,parameter :: N_MODELS=4,N_EPOCHS=361,N_SAMPLES=4096
    integer,parameter :: N_DA_VARS=6,N_BOUNDARY=12,DA_MAX_ORDER=8
    integer,parameter :: MAX_STEPS_PER_HOUR=1000
    real(DP),parameter :: OUTPUT_STEP_S=3600.0_DP
    real(DP),parameter :: REL_TOL=1.0e-12_DP,ABS_TOL=1.0e-12_DP
    real(DP),parameter :: DT_MIN_S=1.0e-6_DP,DT_MAX_S=3600.0_DP
    real(DP),parameter :: POSITION_ABS_TOL_KM=1.0e-2_DP
    real(DP),parameter :: VELOCITY_ABS_TOL_KMS=1.0e-8_DP
    real(DP),parameter :: ANGLE_ABS_TOL_RAD=1.0e-10_DP
    real(DP),parameter :: VALIDATION_REL_TOL=1.0e-3_DP
    real(DP),parameter :: NONLINEAR_LIMIT=0.05_DP
    real(DP),parameter :: SKEW_LIMIT=0.10_DP,KURT_LIMIT=0.20_DP
    real(DP),parameter :: RANK_TOL=1.0e-8_DP
    real(DP),parameter :: RAD_TO_ARCSEC=206264.80624709636_DP
    real(DP),parameter :: SPACECRAFT_MASS_KG=1200.0_DP
    real(DP),parameter :: CANNONBALL_AREA_M2=9.0_DP,CANNONBALL_CR=1.25_DP
    real(DP),parameter :: CANNONBALL_SMR=CANNONBALL_AREA_M2/SPACECRAFT_MASS_KG
    real(DP),parameter :: SOLAR_PRESSURE_1AU=1367.0_DP/299792458.0_DP
    character(len=*),parameter :: CONFIG_FILE='config/config.txt'
    character(len=*),parameter :: HALO_OPM='OPM/L1Halo-1/L1Halo-1_init.opm.json'
    character(len=*),parameter :: SITE_FILE='config/site-used.json'
    character(len=*),parameter :: OUTPUT_DIR='SRP/260907_initial_state_uncertainty'
    character(len=32),parameter :: MODEL_NAMES(N_MODELS)=[character(len=32) :: &
        'cannonball','box_wing_sun','box_wing_earth','box_wing_moon']
    real(DP),parameter :: INITIAL_SIGMA(6)= &
        [10.0_DP,10.0_DP,10.0_DP,3.0e-5_DP,3.0e-5_DP,3.0e-5_DP]

    real(DP) :: epoch0,state0(6),covariance0(6,6),initial_nd(6),sigma_nd(6)
    real(DP) :: gaussian_q(N_DA_VARS,N_SAMPLES),weights(N_SAMPLES)
    real(DP) :: boundary_q(N_DA_VARS,N_BOUNDARY)
    real(DP) :: nominal_store(6,N_EPOCHS,N_MODELS)
    real(DP) :: mean_store(6,N_EPOCHS,N_MODELS)
    real(DP) :: covariance_store(6,6,N_EPOCHS,N_MODELS)
    real(DP) :: pos_axis_store(3,N_EPOCHS,N_MODELS)
    real(DP) :: tangent_cov_store(2,2,N_EPOCHS,N_MODELS)
    real(DP) :: tangent_axis_store(2,N_EPOCHS,N_MODELS)
    real(DP) :: first_distorted_hour(N_MODELS),max_eta_nl(N_MODELS)
    real(DP) :: max_eta_cov(N_MODELS),max_eta_mean(N_MODELS)
    real(DP) :: max_abs_skew(N_MODELS),max_abs_kurt(N_MODELS)
    type(station_record),allocatable :: station_list(:)
    type(observation_station) :: r91
    integer :: history_unit,validation_unit,comparison_unit,summary_unit,sampling_unit
    integer :: performance_unit
    integer :: model_id,status,exit_status,total_validation_records
    integer :: failed_validation_records,handles_at_start
    integer :: validation_failures_by_model(N_MODELS)
    real(DP) :: timing_real(N_MODELS),timing_da_order(3,N_MODELS)
    real(DP) :: timing_state_batch(N_MODELS),timing_angle_batch(N_MODELS)
    real(DP) :: timing_analysis_other(N_MODELS),timing_total(N_MODELS)
    real(DP) :: model_time_start,model_time_end

    handles_at_start=active_da_count()
    call assert_true(handles_at_start==0,'DA handles nonzero before test')
    call pod_engine_init(CONFIG_FILE)
    call init_gravity_network_da()
    call load_initial_opm(HALO_OPM,epoch0,state0,covariance0)
    call set_propagation_epoch_real(epoch0)
    call set_propagation_epoch_da(epoch0)
    call preload_stations(SITE_FILE,station_list)
    r91=find_station_by_id('R91',station_list)
    initial_nd(1:3)=state0(1:3)/config%LU
    initial_nd(4:6)=state0(4:6)/config%VU
    sigma_nd(1:3)=INITIAL_SIGMA(1:3)/config%LU
    sigma_nd(4:6)=INITIAL_SIGMA(4:6)/config%VU
    call set_srp_ballistic_parameters(Cr=CANNONBALL_CR,SMR=CANNONBALL_SMR, &
                                      RP=SOLAR_PRESSURE_1AU)
    call clear_srp_scale_error()
    call clear_srp_scale_uncertainty()

    call assert_true(da_order_for_output_hour(0)==4,'0 h order')
    call assert_true(da_order_for_output_hour(100)==4,'100 h order')
    call assert_true(da_order_for_output_hour(101)==6,'101 h order')
    call assert_true(da_order_for_output_hour(250)==6,'250 h order')
    call assert_true(da_order_for_output_hour(251)==8,'251 h order')
    call generate_antithetic_standard_normal_samples(gaussian_q,260907,status)
    call assert_true(status==0,'failed to generate Gaussian ensemble')
    weights=1.0_DP/real(N_SAMPLES,DP)
    call build_boundary_points(boundary_q)

    call execute_command_line('mkdir -p '//OUTPUT_DIR,exitstat=exit_status)
    call assert_true(exit_status==0,'failed to create output directory')
    call open_outputs()
    call verify_initial_samples()
    first_distorted_hour=-1.0_DP
    max_eta_nl=0.0_DP
    max_eta_cov=0.0_DP
    max_eta_mean=0.0_DP
    max_abs_skew=0.0_DP
    max_abs_kurt=0.0_DP
    timing_real=0.0_DP
    timing_da_order=0.0_DP
    timing_state_batch=0.0_DP
    timing_angle_batch=0.0_DP
    timing_analysis_other=0.0_DP
    timing_total=0.0_DP
    total_validation_records=0
    failed_validation_records=0
    validation_failures_by_model=0

    write(*,'(a)') '15-day initial-state uncertainty by deterministic SRP model'
    do model_id=1,N_MODELS
        call configure_srp_case(model_id)
        write(*,'(2x,a)') trim(MODEL_NAMES(model_id))
        call cpu_time(model_time_start)
        call run_one_model(model_id)
        call cpu_time(model_time_end)
        timing_total(model_id)=model_time_end-model_time_start
        call write_performance_row(model_id)
        write(*,'(4x,a,f10.3,a)') 'CPU time: ',timing_total(model_id),' s'
    end do
    call write_summary()

    close(history_unit)
    close(validation_unit)
    close(comparison_unit)
    close(summary_unit)
    close(sampling_unit)
    close(performance_unit)
    call clear_srp_scale_error()
    call clear_srp_scale_uncertainty()
    call cleanup_gravity_network_da()
    if(allocated(station_list)) deallocate(station_list)
    call assert_true(total_validation_records==N_MODELS*N_BOUNDARY*N_EPOCHS, &
                     'unexpected number of 3-sigma validation rows')
    call assert_true(failed_validation_records==0, &
                     'one or more 3-sigma DA/Real closure records failed')
    call assert_true(active_da_count()==handles_at_start,'DA handle leak in full test')
    write(*,'(a)') 'PASS: 15-day initial-state uncertainty outputs written.'

contains

    pure integer function da_order_for_output_hour(hour) result(order)
        integer,intent(in) :: hour
        if(hour<=100) then
            order=4
        else if(hour<=250) then
            order=6
        else
            order=8
        end if
    end function da_order_for_output_hour

    pure integer function da_order_slot(order) result(slot)
        integer,intent(in) :: order
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
    end function da_order_slot

    subroutine build_boundary_points(points)
        real(DP),intent(out) :: points(N_DA_VARS,N_BOUNDARY)
        integer :: i
        points=0.0_DP
        do i=1,N_DA_VARS
            points(i,2*i-1)=3.0_DP
            points(i,2*i)=-3.0_DP
        end do
    end subroutine build_boundary_points

    subroutine configure_srp_case(id)
        integer,intent(in) :: id
        config%use_srp=.true.
        config%srp_mass_kg=SPACECRAFT_MASS_KG
        config%srp_box_dimensions_m=[2.0_DP,3.0_DP,4.0_DP]
        config%srp_box_optical=[0.30_DP,0.40_DP,0.30_DP]
        config%srp_array_total_area_m2=24.0_DP
        config%srp_array_tracking_mode='single_axis'
        config%srp_array_hinge_axis_body=[0.0_DP,1.0_DP,0.0_DP]
        config%srp_array_reference_normal_body=[1.0_DP,0.0_DP,0.0_DP]
        config%srp_array_front_optical=[0.85_DP,0.08_DP,0.07_DP]
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
        call assert_true(validate_config(),'invalid SRP configuration')
    end subroutine configure_srp_case

    subroutine verify_initial_samples()
        real(DP) :: mean_q(6),cov_q(6,6),skew_q(6),kurt_q(6)
        integer :: local_status,i
        call compute_weighted_moments(gaussian_q,weights,mean_q,cov_q, &
                                      skew_q,kurt_q,local_status)
        call assert_true(local_status==0,'initial sample moment failure')
        call assert_true(maxval(abs(mean_q))<1.0e-14_DP,'initial sample mean')
        do i=1,6
            call assert_true(abs(cov_q(i,i)-1.0_DP)<1.0e-12_DP, &
                             'initial sample marginal variance')
        end do
        write(sampling_unit,'(a)') &
            'seed,n_samples,mean_q1,mean_q2,mean_q3,mean_q4,mean_q5,mean_q6,'// &
            'var_q1,var_q2,var_q3,var_q4,var_q5,var_q6,'// &
            'skew_q1,skew_q2,skew_q3,skew_q4,skew_q5,skew_q6,'// &
            'kurt_q1,kurt_q2,kurt_q3,kurt_q4,kurt_q5,kurt_q6'
        write(sampling_unit,'(*(g0,:,","))') 260907,N_SAMPLES,mean_q, &
            [(cov_q(i,i),i=1,6)],skew_q,kurt_q
    end subroutine verify_initial_samples

    subroutine run_one_model(id)
        integer,intent(in) :: id
        real(DP) :: real_states(6,0:N_BOUNDARY,N_EPOCHS)
        real(DP) :: real_angles(2,0:N_BOUNDARY,N_EPOCHS)
        real(DP) :: zero_q(6)
        real(DP) :: time_start,time_end
        integer :: boundary

        zero_q=0.0_DP
        call cpu_time(time_start)
        call propagate_real_history(zero_q,real_states(:,0,:),real_angles(:,0,:))
        do boundary=1,N_BOUNDARY
            call propagate_real_history(boundary_q(:,boundary), &
                real_states(:,boundary,:),real_angles(:,boundary,:))
        end do
        call cpu_time(time_end)
        timing_real(id)=timing_real(id)+time_end-time_start
        call propagate_da_history(id,real_states,real_angles)
    end subroutine run_one_model

    !> Propagate one deterministic initial-state sample at exact hourly endpoints.
    subroutine propagate_real_history(q,state_history,angle_history)
        real(DP),intent(in) :: q(6)
        real(DP),intent(out) :: state_history(6,N_EPOCHS)
        real(DP),intent(out) :: angle_history(2,N_EPOCHS)
        real(DP) :: current_nd(6),physical(6),start_time,end_time
        real(DP),allocatable :: times(:),states(:,:)
        integer :: hour,n_steps

        current_nd=initial_nd+sigma_nd*q
        call nondimensional_to_physical(current_nd,physical)
        state_history(:,1)=physical
        call compute_measurement(physical,epoch0,r91,'OPTICAL',angle_history(:,1))
        do hour=1,N_EPOCHS-1
            start_time=real(hour-1,DP)*OUTPUT_STEP_S/config%TU
            end_time=real(hour,DP)*OUTPUT_STEP_S/config%TU
            call adaptive_step_integrate(current_nd,start_time,end_time, &
                METHOD_RKF78,times,states,n_steps,max_steps_in=MAX_STEPS_PER_HOUR, &
                rel_tol_in=REL_TOL,abs_tol_in=ABS_TOL,dt_min_in=DT_MIN_S, &
                dt_max_in=DT_MAX_S)
            call assert_true(n_steps>1,'Real hourly segment produced no step')
            current_nd=states(n_steps,:)
            call nondimensional_to_physical(current_nd,physical)
            state_history(:,hour+1)=physical
            call compute_measurement(physical,epoch0+real(hour,DP)*OUTPUT_STEP_S, &
                                     r91,'OPTICAL',angle_history(:,hour+1))
            if(allocated(times)) deallocate(times)
            if(allocated(states)) deallocate(states)
        end do
    end subroutine propagate_real_history

    !> Propagate one six-variable DA state and process every hourly map.
    subroutine propagate_da_history(id,real_states,real_angles)
        integer,intent(in) :: id
        real(DP),intent(in) :: real_states(6,0:N_BOUNDARY,N_EPOCHS)
        real(DP),intent(in) :: real_angles(2,0:N_BOUNDARY,N_EPOCHS)
        type(AlgebraicVector) :: current_state,next_state
        type(DA) :: variable,scaled
        real(DP),allocatable :: times(:),nominal_steps(:,:)
        real(DP) :: start_time,end_time,expected,time_start,time_end
        integer :: i,j,hour,n_steps,active_order,handles_before

        handles_before=active_da_count()
        call assert_true(handles_before==0,'DA handles nonzero before model')
        call dace_initialize(DA_MAX_ORDER,N_DA_VARS)
        call dace_set_to(4)
        call current_state%init(6)
        call variable%init()
        call scaled%init()
        do i=1,6
            call variable%destroy()
            call variable%init_var(i)
            call da_mul(variable,sigma_nd(i),scaled)
            call da_add(scaled,initial_nd(i),current_state%elements(i))
        end do
        call variable%destroy()
        call scaled%destroy()
        do i=1,6
            call assert_close(current_state%elements(i)%cons(),initial_nd(i), &
                              1.0e-14_DP,'initial DA constant')
            do j=1,6
                expected=0.0_DP
                if(i==j) expected=sigma_nd(i)
                call assert_close(current_state%elements(i)%get_deriv_value(j), &
                                  expected,1.0e-14_DP,'initial DA derivative')
            end do
        end do
        call analyze_and_write_epoch(id,1,0,current_state,real_states,real_angles)

        do hour=1,N_EPOCHS-1
            active_order=da_order_for_output_hour(hour)
            call dace_set_to(active_order)
            call assert_true(dace_get_to()==active_order,'DACE order not applied')
            start_time=real(hour-1,DP)*OUTPUT_STEP_S/config%TU
            end_time=real(hour,DP)*OUTPUT_STEP_S/config%TU
            call cpu_time(time_start)
            call da_adaptive_step_integrate(current_state,start_time,end_time, &
                METHOD_RKF78,times,nominal_steps,next_state,n_steps, &
                max_steps_in=MAX_STEPS_PER_HOUR,rel_tol_in=REL_TOL, &
                abs_tol_in=ABS_TOL,dt_min_in=DT_MIN_S,dt_max_in=DT_MAX_S)
            call cpu_time(time_end)
            timing_da_order(da_order_slot(active_order),id)= &
                timing_da_order(da_order_slot(active_order),id)+time_end-time_start
            call assert_true(n_steps>1,'DA hourly segment produced no step')
            call analyze_and_write_epoch(id,hour+1,hour,next_state, &
                                         real_states,real_angles)
            current_state=next_state
            call next_state%destroy()
            if(allocated(times)) deallocate(times)
            if(allocated(nominal_steps)) deallocate(nominal_steps)
        end do
        call current_state%destroy()
        call next_state%destroy()
        call assert_true(active_da_count()==handles_before, &
                         'DA handle leak in model propagation')
    end subroutine propagate_da_history

    subroutine analyze_and_write_epoch(id,index,hour,state_da,real_states,real_angles)
        integer,intent(in) :: id,index,hour
        type(AlgebraicVector),intent(in) :: state_da
        real(DP),intent(in) :: real_states(6,0:N_BOUNDARY,N_EPOCHS)
        real(DP),intent(in) :: real_angles(2,0:N_BOUNDARY,N_EPOCHS)
        type(AlgebraicVector) :: physical_da,angle_da
        type(CompiledDA) :: compiled_state,compiled_angle
        real(DP) :: samples(6,N_SAMPLES),linear_samples(6,N_SAMPLES)
        real(DP) :: rtn_samples(6,N_SAMPLES)
        real(DP) :: radec_samples(2,N_SAMPLES),tangent_samples(2,N_SAMPLES)
        real(DP) :: linear_tangent(2,N_SAMPLES),evaluated_nd(6),evaluated(6)
        real(DP) :: evaluated_angle(2),nominal(6),nominal_nd(6),nominal_angle(2)
        real(DP) :: evaluated_nd_batch(6,N_SAMPLES)
        real(DP) :: evaluated_angle_batch(2,N_SAMPLES)
        real(DP) :: derivative(6,6),angle_derivative(2,6),rotation(3,3)
        real(DP) :: block_rotation(6,6),mean_i(6),cov_i(6,6)
        real(DP) :: mean_rtn(6),cov_rtn(6,6),skew6(6),kurt6(6)
        real(DP) :: radec_mean(2),radec_cov(2,2),tangent_mean(2),tangent_cov(2,2)
        real(DP) :: skew2(2),kurt2(2),dummy_skew6(6),dummy_kurt6(6)
        real(DP) :: linear_mean6(6),linear_cov6(6,6)
        real(DP) :: linear_mean2(2),linear_cov2(2,2)
        real(DP) :: pos_eval(3),vel_eval(3),pos_vec_i(3,3),vel_vec_i(3,3)
        real(DP) :: pos_axis_i(3),vel_axis_i(3),pos_vec_rtn(3,3)
        real(DP) :: vel_vec_rtn(3,3),pos_axis_rtn(3),vel_axis_rtn(3)
        real(DP) :: pos_eval_rtn(3),vel_eval_rtn(3)
        real(DP) :: radec_eval(2),radec_vec(2,2),radec_axis(2)
        real(DP) :: tangent_eval(2),tangent_vec(2,2),tangent_axis(2)
        real(DP) :: state_pc_skew,state_pc_kurt,tangent_pc_skew,tangent_pc_kurt
        real(DP) :: eta_nl_state,eta_cov_state,eta_mean_state
        real(DP) :: eta_nl_tangent,eta_cov_tangent,eta_mean_tangent
        real(DP) :: position_rms,velocity_rms,tangent_rms
        real(DP) :: cov_i_upper(21),cov_rtn_upper(21),radec_upper(3),tangent_upper(3)
        real(DP) :: unit_scale,analysis_start,analysis_end
        real(DP) :: state_batch_start,state_batch_end
        real(DP) :: angle_batch_start,angle_batch_end
        real(DP) :: state_batch_elapsed,angle_batch_elapsed
        integer :: i,j,local_status,rank_state,rank_tangent
        logical :: distorted
        character(len=5) :: distorted_text

        call cpu_time(analysis_start)
        nominal_nd=state_da%cons()
        call nondimensional_to_physical(nominal_nd,nominal)
        do i=1,6
            unit_scale=merge(config%LU,config%VU,i<=3)
            do j=1,6
                derivative(i,j)=unit_scale* &
                    state_da%elements(i)%get_deriv_value(j)
            end do
        end do

        call cpu_time(state_batch_start)
        compiled_state=state_da%compile()
        call compiled_state%eval_batch_into(gaussian_q,evaluated_nd_batch,local_status)
        call assert_true(local_status==0,'compiled state batch evaluation')
        call cpu_time(state_batch_end)
        state_batch_elapsed=state_batch_end-state_batch_start
        timing_state_batch(id)=timing_state_batch(id)+state_batch_elapsed
        if(index==1) then
            call compiled_state%eval_into(gaussian_q(:,1),evaluated_nd,local_status)
            call assert_true(local_status==0,'compiled state scalar cross-check')
            call assert_close(maxval(abs(evaluated_nd-evaluated_nd_batch(:,1))), &
                              0.0_DP,1.0e-14_DP,'state scalar/batch equivalence')
        end if
        do j=1,N_SAMPLES
            evaluated_nd=evaluated_nd_batch(:,j)
            call nondimensional_to_physical(evaluated_nd,evaluated)
            samples(:,j)=evaluated-nominal
            linear_samples(:,j)=matmul(derivative,gaussian_q(:,j))
        end do
        call compute_weighted_moments(samples,weights,mean_i,cov_i, &
                                      skew6,kurt6,local_status)
        call assert_true(local_status==0,'J2000 moment computation')
        call compute_weighted_moments(linear_samples,weights,linear_mean6, &
            linear_cov6,dummy_skew6,dummy_kurt6,local_status)
        call assert_true(local_status==0,'linear state moment computation')
        call validate_covariance(cov_i,'J2000 covariance')

        call build_rtn_rotation(nominal(1:3),nominal(4:6),rotation,local_status)
        call assert_true(local_status==0,'undefined nominal RTN frame')
        block_rotation=0.0_DP
        block_rotation(1:3,1:3)=rotation
        block_rotation(4:6,4:6)=rotation
        rtn_samples=matmul(block_rotation,samples)
        call compute_weighted_moments(rtn_samples,weights,mean_rtn,cov_rtn, &
                                      skew6,kurt6,local_status)
        call assert_true(local_status==0,'RTN moment computation')
        call validate_covariance(cov_rtn,'RTN covariance')
        call transform_covariance6_to_rtn(cov_i,rotation,linear_cov6)
        call assert_close(maxval(abs(cov_rtn-linear_cov6)),0.0_DP, &
            1.0e-8_DP*max(1.0_DP,maxval(abs(cov_i))),'RTN sample covariance')

        call compute_covariance_axes(cov_i(1:3,1:3),pos_eval,pos_vec_i, &
                                     pos_axis_i,local_status)
        call assert_true(local_status==0,'J2000 position axes')
        call compute_covariance_axes(cov_i(4:6,4:6),vel_eval,vel_vec_i, &
                                     vel_axis_i,local_status)
        call assert_true(local_status==0,'J2000 velocity axes')
        call compute_covariance_axes(cov_rtn(1:3,1:3),pos_eval_rtn,pos_vec_rtn, &
                                     pos_axis_rtn,local_status)
        call assert_true(local_status==0,'RTN position axes')
        call compute_covariance_axes(cov_rtn(4:6,4:6),vel_eval_rtn,vel_vec_rtn, &
                                     vel_axis_rtn,local_status)
        call assert_true(local_status==0,'RTN velocity axes')
        call assert_close(maxval(abs(pos_eval-pos_eval_rtn)),0.0_DP, &
            1.0e-8_DP*max(1.0_DP,maxval(pos_eval)),'J2000/RTN position eigenvalues')
        call assert_close(maxval(abs(vel_eval-vel_eval_rtn)),0.0_DP, &
            1.0e-8_DP*max(1.0_DP,maxval(vel_eval)),'J2000/RTN velocity eigenvalues')
        call compute_effective_rank(cov_i,RANK_TOL,rank_state,local_status)
        call assert_true(local_status==0,'state effective rank')

        call cpu_time(angle_batch_start)
        call physical_da%init(6)
        do i=1,6
            unit_scale=merge(config%LU,config%VU,i<=3)
            call da_mul(state_da%elements(i),unit_scale,physical_da%elements(i))
        end do
        call compute_measurement_da(physical_da,epoch0+real(hour,DP)*OUTPUT_STEP_S, &
                                    r91,'OPTICAL',angle_da)
        nominal_angle=angle_da%cons()
        call assert_true(sqrt(sum((nominal(1:3)-real_states(1:3,0,index))**2))<= &
                         POSITION_ABS_TOL_KM,'nominal DA/Real position mismatch')
        call assert_true(sqrt(sum((nominal(4:6)-real_states(4:6,0,index))**2))<= &
                         VELOCITY_ABS_TOL_KMS,'nominal DA/Real velocity mismatch')
        call assert_true(sqrt(wrap_angle_rad(nominal_angle(1)- &
                         real_angles(1,0,index))**2+(nominal_angle(2)- &
                         real_angles(2,0,index))**2)<=ANGLE_ABS_TOL_RAD, &
                         'nominal DA/Real R91 angle mismatch')
        do i=1,2
            do j=1,6
                angle_derivative(i,j)=angle_da%elements(i)%get_deriv_value(j)
            end do
        end do
        compiled_angle=angle_da%compile()
        call compiled_angle%eval_batch_into(gaussian_q,evaluated_angle_batch, &
                                            local_status)
        call assert_true(local_status==0,'compiled RADEC batch evaluation')
        call cpu_time(angle_batch_end)
        angle_batch_elapsed=angle_batch_end-angle_batch_start
        timing_angle_batch(id)=timing_angle_batch(id)+angle_batch_elapsed
        if(index==1) then
            call compiled_angle%eval_into(gaussian_q(:,1),evaluated_angle,local_status)
            call assert_true(local_status==0,'compiled RADEC scalar cross-check')
            call assert_close(maxval(abs(evaluated_angle-evaluated_angle_batch(:,1))), &
                              0.0_DP,1.0e-14_DP,'RADEC scalar/batch equivalence')
        end if
        do j=1,N_SAMPLES
            evaluated_angle=evaluated_angle_batch(:,j)
            radec_samples(1,j)=wrap_angle_rad(evaluated_angle(1)-nominal_angle(1))
            radec_samples(2,j)=evaluated_angle(2)-nominal_angle(2)
            tangent_samples(1,j)=radec_samples(1,j)*cos(nominal_angle(2))
            tangent_samples(2,j)=radec_samples(2,j)
            linear_tangent(1,j)=cos(nominal_angle(2))* &
                dot_product(angle_derivative(1,:),gaussian_q(:,j))
            linear_tangent(2,j)=dot_product(angle_derivative(2,:),gaussian_q(:,j))
        end do
        call compute_weighted_moments(radec_samples,weights,radec_mean,radec_cov, &
                                      skew2,kurt2,local_status)
        call assert_true(local_status==0,'RADEC moment computation')
        call compute_weighted_moments(tangent_samples,weights,tangent_mean, &
                                      tangent_cov,skew2,kurt2,local_status)
        call assert_true(local_status==0,'tangent moment computation')
        call compute_weighted_moments(linear_tangent,weights,linear_mean2, &
            linear_cov2,skew2,kurt2,local_status)
        call assert_true(local_status==0,'linear tangent moment computation')
        call validate_covariance(radec_cov,'RADEC covariance')
        call validate_covariance(tangent_cov,'tangent covariance')
        call covariance_axes2(radec_cov,radec_eval,radec_vec,radec_axis,local_status)
        call assert_true(local_status==0,'RADEC covariance axes')
        call covariance_axes2(tangent_cov,tangent_eval,tangent_vec,tangent_axis, &
                              local_status)
        call assert_true(local_status==0,'tangent covariance axes')
        call compute_effective_rank(tangent_cov,RANK_TOL,rank_tangent,local_status)
        call assert_true(local_status==0,'tangent effective rank')

        call principal_shape(samples(1:3,:),mean_i(1:3),pos_eval,pos_vec_i, &
                             state_pc_skew,state_pc_kurt)
        call principal_shape(samples(4:6,:),mean_i(4:6),vel_eval,vel_vec_i, &
                             skew2(1),kurt2(1))
        state_pc_skew=max(state_pc_skew,skew2(1))
        state_pc_kurt=max(state_pc_kurt,kurt2(1))
        call principal_shape(tangent_samples,tangent_mean,tangent_eval,tangent_vec, &
                             tangent_pc_skew,tangent_pc_kurt)

        position_rms=sqrt(max(sum(pos_eval),0.0_DP))
        velocity_rms=sqrt(max(sum(vel_eval),0.0_DP))
        tangent_rms=sqrt(max(sum(tangent_eval),0.0_DP))
        eta_nl_state=max(normalized_rms_difference(samples(1:3,:), &
            linear_samples(1:3,:),position_rms),normalized_rms_difference( &
            samples(4:6,:),linear_samples(4:6,:),velocity_rms))
        eta_cov_state=max(relative_frobenius(cov_i(1:3,1:3), &
            covariance_of(linear_samples(1:3,:))),relative_frobenius( &
            cov_i(4:6,4:6),covariance_of(linear_samples(4:6,:))))
        eta_mean_state=max(normalized_mean(mean_i(1:3),position_rms), &
                           normalized_mean(mean_i(4:6),velocity_rms))
        eta_nl_tangent=normalized_rms_difference(tangent_samples,linear_tangent, &
                                                 tangent_rms)
        eta_cov_tangent=relative_frobenius(tangent_cov,linear_cov2)
        eta_mean_tangent=normalized_mean(tangent_mean,tangent_rms)
        distorted=eta_nl_state>NONLINEAR_LIMIT .or. &
            eta_cov_state>NONLINEAR_LIMIT .or. &
            eta_nl_tangent>NONLINEAR_LIMIT .or. eta_cov_tangent>NONLINEAR_LIMIT .or. &
            state_pc_skew>SKEW_LIMIT .or. tangent_pc_skew>SKEW_LIMIT .or. &
            state_pc_kurt>KURT_LIMIT .or. tangent_pc_kurt>KURT_LIMIT
        distorted_text=merge('true ','false',distorted)

        nominal_store(:,index,id)=nominal
        mean_store(:,index,id)=mean_i
        covariance_store(:,:,index,id)=cov_i
        pos_axis_store(:,index,id)=pos_axis_i
        tangent_cov_store(:,:,index,id)=tangent_cov
        tangent_axis_store(:,index,id)=tangent_axis
        max_eta_nl(id)=max(max_eta_nl(id),eta_nl_state,eta_nl_tangent)
        max_eta_cov(id)=max(max_eta_cov(id),eta_cov_state,eta_cov_tangent)
        max_eta_mean(id)=max(max_eta_mean(id),eta_mean_state,eta_mean_tangent)
        max_abs_skew(id)=max(max_abs_skew(id),state_pc_skew,tangent_pc_skew)
        max_abs_kurt(id)=max(max_abs_kurt(id),state_pc_kurt,tangent_pc_kurt)
        if(distorted .and. first_distorted_hour(id)<0.0_DP) &
            first_distorted_hour(id)=real(hour,DP)

        cov_i_upper=pack_upper(cov_i)
        cov_rtn_upper=pack_upper(cov_rtn)
        radec_upper=pack_upper2(radec_cov)
        tangent_upper=pack_upper2(tangent_cov)
        write(history_unit,'(*(g0,:,","))') trim(MODEL_NAMES(id)),hour, &
            da_order_for_output_hour(hour),nominal,mean_i,mean_rtn,cov_i_upper, &
            cov_rtn_upper,pos_axis_i,pos_axis_rtn,vel_axis_i,vel_axis_rtn, &
            pos_vec_i,pos_vec_rtn,vel_vec_i,vel_vec_rtn,radec_mean,radec_upper, &
            radec_axis,radec_vec,radec_axis*RAD_TO_ARCSEC,tangent_mean,tangent_upper, &
            tangent_axis,tangent_vec,tangent_axis*RAD_TO_ARCSEC, &
            position_rms,velocity_rms,tangent_rms,eta_nl_state,eta_cov_state, &
            eta_mean_state,eta_nl_tangent,eta_cov_tangent,eta_mean_tangent, &
            state_pc_skew,state_pc_kurt,tangent_pc_skew,tangent_pc_kurt, &
            rank_state,rank_tangent,trim(distorted_text)

        call write_validation_rows(id,index,hour,compiled_state,compiled_angle, &
            nominal,nominal_angle,real_states,real_angles)
        if(id>1) call write_model_comparison(id,index,hour,position_rms, &
                                             velocity_rms,tangent_rms)

        call compiled_state%destroy()
        call compiled_angle%destroy()
        call angle_da%destroy()
        call physical_da%destroy()
        call cpu_time(analysis_end)
        timing_analysis_other(id)=timing_analysis_other(id)+max(0.0_DP, &
            analysis_end-analysis_start-state_batch_elapsed-angle_batch_elapsed)
    end subroutine analyze_and_write_epoch

    subroutine nondimensional_to_physical(nd,physical)
        real(DP),intent(in) :: nd(6)
        real(DP),intent(out) :: physical(6)
        physical(1:3)=nd(1:3)*config%LU
        physical(4:6)=nd(4:6)*config%VU
    end subroutine nondimensional_to_physical

    !> Compare all twelve axis-aligned three-sigma points against Real propagation.
    subroutine write_validation_rows(id,index,hour,state_map,angle_map, &
                                     nominal_da,angle_nominal_da,real_states,real_angles)
        integer,intent(in) :: id,index,hour
        type(CompiledDA),intent(in) :: state_map,angle_map
        real(DP),intent(in) :: nominal_da(6),angle_nominal_da(2)
        real(DP),intent(in) :: real_states(6,0:N_BOUNDARY,N_EPOCHS)
        real(DP),intent(in) :: real_angles(2,0:N_BOUNDARY,N_EPOCHS)
        real(DP) :: evaluated_nd(6),evaluated(6),evaluated_angle(2)
        real(DP) :: da_error(6),real_error(6),difference(6)
        real(DP) :: da_angle_error(2),real_angle_error(2),angle_difference(2)
        real(DP) :: pos_diff,vel_diff,angle_diff,pos_tol,vel_tol,angle_tol
        integer :: boundary,axis,local_status
        logical :: passed
        character(len=5) :: passed_text

        do boundary=1,N_BOUNDARY
            axis=(boundary+1)/2
            call state_map%eval_into(boundary_q(:,boundary),evaluated_nd,local_status)
            call assert_true(local_status==0,'3-sigma state map evaluation')
            call nondimensional_to_physical(evaluated_nd,evaluated)
            da_error=evaluated-nominal_da
            real_error=real_states(:,boundary,index)-real_states(:,0,index)
            difference=da_error-real_error
            call angle_map%eval_into(boundary_q(:,boundary),evaluated_angle,local_status)
            call assert_true(local_status==0,'3-sigma angle map evaluation')
            da_angle_error=[wrap_angle_rad(evaluated_angle(1)-angle_nominal_da(1)), &
                            evaluated_angle(2)-angle_nominal_da(2)]
            real_angle_error=[wrap_angle_rad(real_angles(1,boundary,index)- &
                                             real_angles(1,0,index)), &
                              real_angles(2,boundary,index)-real_angles(2,0,index)]
            angle_difference=[wrap_angle_rad(da_angle_error(1)-real_angle_error(1)), &
                              da_angle_error(2)-real_angle_error(2)]
            pos_diff=sqrt(sum(difference(1:3)**2))
            vel_diff=sqrt(sum(difference(4:6)**2))
            angle_diff=sqrt(sum(angle_difference**2))
            pos_tol=max(POSITION_ABS_TOL_KM,VALIDATION_REL_TOL* &
                        sqrt(sum(real_error(1:3)**2)))
            vel_tol=max(VELOCITY_ABS_TOL_KMS,VALIDATION_REL_TOL* &
                        sqrt(sum(real_error(4:6)**2)))
            angle_tol=max(ANGLE_ABS_TOL_RAD,VALIDATION_REL_TOL* &
                          sqrt(sum(real_angle_error**2)))
            passed=pos_diff<=pos_tol .and. vel_diff<=vel_tol .and. &
                   angle_diff<=angle_tol
            passed_text=merge('true ','false',passed)
            total_validation_records=total_validation_records+1
            if(.not.passed) then
                failed_validation_records=failed_validation_records+1
                validation_failures_by_model(id)=validation_failures_by_model(id)+1
            end if
            write(validation_unit,'(*(g0,:,","))') trim(MODEL_NAMES(id)),hour, &
                da_order_for_output_hour(hour),axis,boundary_q(axis,boundary), &
                da_error,real_error,difference,pos_diff,vel_diff,pos_tol,vel_tol, &
                da_angle_error,real_angle_error,angle_difference,angle_diff, &
                angle_tol,trim(passed_text)
        end do
    end subroutine write_validation_rows

    subroutine write_model_comparison(id,index,hour,position_rms,velocity_rms,tangent_rms)
        integer,intent(in) :: id,index,hour
        real(DP),intent(in) :: position_rms,velocity_rms,tangent_rms
        real(DP) :: ref_pos_rms,ref_vel_rms,ref_tangent_rms
        real(DP) :: nominal_difference(6),mean_difference(6)
        real(DP) :: pos_cov_difference,vel_cov_difference,tangent_cov_difference

        ref_pos_rms=sqrt(max(trace_matrix(covariance_store(1:3,1:3,index,1)),0.0_DP))
        ref_vel_rms=sqrt(max(trace_matrix(covariance_store(4:6,4:6,index,1)),0.0_DP))
        ref_tangent_rms=sqrt(max(trace_matrix(tangent_cov_store(:,:,index,1)),0.0_DP))
        nominal_difference=nominal_store(:,index,id)-nominal_store(:,index,1)
        mean_difference=mean_store(:,index,id)-mean_store(:,index,1)
        pos_cov_difference=frobenius_norm(covariance_store(1:3,1:3,index,id)- &
            covariance_store(1:3,1:3,index,1))/max(frobenius_norm( &
            covariance_store(1:3,1:3,index,1)),tiny(1.0_DP))
        vel_cov_difference=frobenius_norm(covariance_store(4:6,4:6,index,id)- &
            covariance_store(4:6,4:6,index,1))/max(frobenius_norm( &
            covariance_store(4:6,4:6,index,1)),tiny(1.0_DP))
        tangent_cov_difference=frobenius_norm(tangent_cov_store(:,:,index,id)- &
            tangent_cov_store(:,:,index,1))/max(frobenius_norm( &
            tangent_cov_store(:,:,index,1)),tiny(1.0_DP))
        write(comparison_unit,'(*(g0,:,","))') trim(MODEL_NAMES(id)),hour, &
            nominal_difference,mean_difference,position_rms/ref_pos_rms, &
            velocity_rms/ref_vel_rms,tangent_rms/ref_tangent_rms, &
            pos_axis_store(1,index,id)/pos_axis_store(1,index,1), &
            tangent_axis_store(1,index,id)/tangent_axis_store(1,index,1), &
            pos_cov_difference,vel_cov_difference,tangent_cov_difference
    end subroutine write_model_comparison

    subroutine validate_covariance(covariance,label)
        real(DP),intent(in) :: covariance(:,:)
        character(len=*),intent(in) :: label
        real(DP),allocatable :: symmetric(:,:),wr(:),wi(:),vectors(:,:)
        real(DP) :: scale
        integer :: i,n,info
        call assert_true(size(covariance,1)==size(covariance,2),trim(label)//' square')
        call assert_true(all(ieee_is_finite(covariance)),trim(label)//' finite')
        scale=max(1.0_DP,maxval(abs(covariance)))
        call assert_true(maxval(abs(covariance-transpose(covariance)))<= &
                         1.0e-10_DP*scale,trim(label)//' symmetric')
        do i=1,size(covariance,1)
            call assert_true(covariance(i,i)>=-1.0e-12_DP*scale, &
                             trim(label)//' nonnegative diagonal')
        end do
        n=size(covariance,1)
        allocate(symmetric(n,n),wr(n),wi(n),vectors(n,n))
        symmetric=0.5_DP*(covariance+transpose(covariance))
        call eigenvalue_decomposition(symmetric,wr,wi,vectors,info)
        call assert_true(info==0,trim(label)//' eigensolver')
        call assert_true(maxval(abs(wi))<=1.0e-10_DP*scale, &
                         trim(label)//' real eigenvalues')
        call assert_true(minval(wr)>=-1.0e-10_DP*scale, &
                         trim(label)//' positive semidefinite')
        deallocate(symmetric,wr,wi,vectors)
    end subroutine validate_covariance

    subroutine covariance_axes2(covariance,eigenvalues,eigenvectors,axis_sigma,local_status)
        real(DP),intent(in) :: covariance(2,2)
        real(DP),intent(out) :: eigenvalues(2),eigenvectors(2,2),axis_sigma(2)
        integer,intent(out) :: local_status
        real(DP) :: a,b,d,root,scale,vector_norm

        a=covariance(1,1)
        b=0.5_DP*(covariance(1,2)+covariance(2,1))
        d=covariance(2,2)
        root=sqrt(max((0.5_DP*(a-d))**2+b*b,0.0_DP))
        eigenvalues=[0.5_DP*(a+d)+root,0.5_DP*(a+d)-root]
        scale=max(1.0_DP,maxval(abs(covariance)))
        local_status=0
        if(eigenvalues(2)<-1.0e-12_DP*scale) then
            local_status=-1
            eigenvectors=0.0_DP
            axis_sigma=0.0_DP
            return
        end if
        eigenvalues=max(eigenvalues,0.0_DP)
        if(abs(b)>100.0_DP*epsilon(1.0_DP)*scale) then
            eigenvectors(:,1)=[b,eigenvalues(1)-a]
            vector_norm=sqrt(sum(eigenvectors(:,1)**2))
            eigenvectors(:,1)=eigenvectors(:,1)/vector_norm
        else if(a>=d) then
            eigenvectors(:,1)=[1.0_DP,0.0_DP]
        else
            eigenvectors(:,1)=[0.0_DP,1.0_DP]
        end if
        eigenvectors(:,2)=[-eigenvectors(2,1),eigenvectors(1,1)]
        axis_sigma=sqrt(eigenvalues)
    end subroutine covariance_axes2

    subroutine principal_shape(samples,mean_value,eigenvalues,eigenvectors, &
                               maximum_skew,maximum_kurtosis)
        real(DP),intent(in) :: samples(:,:),mean_value(:),eigenvalues(:)
        real(DP),intent(in) :: eigenvectors(:,:)
        real(DP),intent(out) :: maximum_skew,maximum_kurtosis
        real(DP) :: projected(1,N_SAMPLES),projected_mean(1),projected_cov(1,1)
        real(DP) :: projected_skew(1),projected_kurt(1),largest
        integer :: axis,j,local_status

        maximum_skew=0.0_DP
        maximum_kurtosis=0.0_DP
        largest=max(eigenvalues(1),tiny(1.0_DP))
        do axis=1,size(eigenvalues)
            if(eigenvalues(axis)/largest<RANK_TOL) cycle
            do j=1,N_SAMPLES
                projected(1,j)=dot_product(samples(:,j)-mean_value, &
                                           eigenvectors(:,axis))
            end do
            call compute_weighted_moments(projected,weights,projected_mean, &
                projected_cov,projected_skew,projected_kurt,local_status)
            call assert_true(local_status==0,'principal-axis moment computation')
            maximum_skew=max(maximum_skew,abs(projected_skew(1)))
            maximum_kurtosis=max(maximum_kurtosis,abs(projected_kurt(1)))
        end do
    end subroutine principal_shape

    function covariance_of(samples) result(covariance)
        real(DP),intent(in) :: samples(:,:)
        real(DP) :: covariance(size(samples,1),size(samples,1))
        real(DP) :: mean_value(size(samples,1)),deviation(size(samples,1))
        integer :: j,i
        mean_value=sum(samples,dim=2)/real(size(samples,2),DP)
        covariance=0.0_DP
        do j=1,size(samples,2)
            deviation=samples(:,j)-mean_value
            do i=1,size(samples,1)
                covariance(:,i)=covariance(:,i)+ &
                    deviation*deviation(i)/real(size(samples,2),DP)
            end do
        end do
        covariance=0.5_DP*(covariance+transpose(covariance))
    end function covariance_of

    function normalized_rms_difference(samples,linear_samples,scale) result(value)
        real(DP),intent(in) :: samples(:,:),linear_samples(:,:),scale
        real(DP) :: value,total
        integer :: j
        total=0.0_DP
        do j=1,size(samples,2)
            total=total+sum((samples(:,j)-linear_samples(:,j))**2)/ &
                        real(size(samples,2),DP)
        end do
        value=sqrt(max(total,0.0_DP))/max(scale,tiny(1.0_DP))
    end function normalized_rms_difference

    pure function relative_frobenius(actual,linear) result(value)
        real(DP),intent(in) :: actual(:,:),linear(:,:)
        real(DP) :: value
        value=frobenius_norm(actual-linear)/max(frobenius_norm(actual),tiny(1.0_DP))
    end function relative_frobenius

    pure function normalized_mean(mean_value,scale) result(value)
        real(DP),intent(in) :: mean_value(:),scale
        real(DP) :: value
        value=sqrt(sum(mean_value**2))/max(scale,tiny(1.0_DP))
    end function normalized_mean

    pure function frobenius_norm(matrix) result(value)
        real(DP),intent(in) :: matrix(:,:)
        real(DP) :: value
        value=sqrt(sum(matrix**2))
    end function frobenius_norm

    pure function trace_matrix(matrix) result(value)
        real(DP),intent(in) :: matrix(:,:)
        real(DP) :: value
        integer :: i
        value=0.0_DP
        do i=1,min(size(matrix,1),size(matrix,2))
            value=value+matrix(i,i)
        end do
    end function trace_matrix

    pure function pack_upper(matrix) result(packed)
        real(DP),intent(in) :: matrix(6,6)
        real(DP) :: packed(21)
        integer :: i,j,k
        k=0
        do i=1,6
            do j=i,6
                k=k+1
                packed(k)=matrix(i,j)
            end do
        end do
    end function pack_upper

    pure function pack_upper2(matrix) result(packed)
        real(DP),intent(in) :: matrix(2,2)
        real(DP) :: packed(3)
        packed=[matrix(1,1),matrix(1,2),matrix(2,2)]
    end function pack_upper2

    subroutine open_outputs()
        integer :: ios
        open(newunit=history_unit,file=OUTPUT_DIR//'/initial_state_uncertainty_history.csv', &
             status='replace',action='write',iostat=ios)
        call assert_true(ios==0,'cannot open history CSV')
        open(newunit=validation_unit,file=OUTPUT_DIR//'/initial_state_da_real_3sigma.csv', &
             status='replace',action='write',iostat=ios)
        call assert_true(ios==0,'cannot open validation CSV')
        open(newunit=comparison_unit,file=OUTPUT_DIR//'/initial_state_model_comparison.csv', &
             status='replace',action='write',iostat=ios)
        call assert_true(ios==0,'cannot open model comparison CSV')
        open(newunit=summary_unit,file=OUTPUT_DIR//'/initial_state_uncertainty_summary.csv', &
             status='replace',action='write',iostat=ios)
        call assert_true(ios==0,'cannot open summary CSV')
        open(newunit=sampling_unit,file=OUTPUT_DIR//'/initial_state_sampling_check.csv', &
             status='replace',action='write',iostat=ios)
        call assert_true(ios==0,'cannot open sampling check CSV')
        open(newunit=performance_unit,file=OUTPUT_DIR//'/initial_state_performance.csv', &
             status='replace',action='write',iostat=ios)
        call assert_true(ios==0,'cannot open performance CSV')
        call write_history_header()
        write(validation_unit,'(a)') &
            'model,time_hours,da_order,initial_axis,q_sigma,da_dx_km,da_dy_km,'// &
            'da_dz_km,da_dvx_kms,da_dvy_kms,da_dvz_kms,real_dx_km,real_dy_km,'// &
            'real_dz_km,real_dvx_kms,real_dvy_kms,real_dvz_kms,diff_x_km,'// &
            'diff_y_km,diff_z_km,diff_vx_kms,diff_vy_kms,diff_vz_kms,'// &
            'position_diff_norm_km,velocity_diff_norm_kms,position_tolerance_km,'// &
            'velocity_tolerance_kms,da_dra_rad,da_ddec_rad,real_dra_rad,'// &
            'real_ddec_rad,diff_ra_rad,diff_dec_rad,angle_diff_norm_rad,'// &
            'angle_tolerance_rad,passed'
        write(comparison_unit,'(a)') &
            'model_vs_cannonball,time_hours,nominal_dx_km,nominal_dy_km,'// &
            'nominal_dz_km,nominal_dvx_kms,nominal_dvy_kms,nominal_dvz_kms,'// &
            'mean_dx_km,mean_dy_km,mean_dz_km,mean_dvx_kms,mean_dvy_kms,'// &
            'mean_dvz_kms,position_rms_ratio,velocity_rms_ratio,'// &
            'tangent_rms_ratio,position_long_axis_ratio,tangent_long_axis_ratio,'// &
            'position_covariance_relative_frobenius_difference,'// &
            'velocity_covariance_relative_frobenius_difference,'// &
            'tangent_covariance_relative_frobenius_difference'
        write(summary_unit,'(a)') &
            'model,first_distorted_hour,max_eta_nl,max_eta_cov,max_eta_mean,'// &
            'max_abs_pc_skew,max_abs_pc_excess_kurtosis,'// &
            'final_position_rms_km,final_position_long_sigma_km,'// &
            'final_position_short_sigma_km,final_velocity_rms_kms,'// &
            'final_tangent_rms_rad,final_tangent_long_sigma_rad,'// &
            'final_tangent_short_sigma_rad,validation_records,validation_failures'
        write(performance_unit,'(a)') &
            'model,n_samples,n_epochs,real_propagation_cpu_s,'// &
            'da_integration_order4_cpu_s,da_integration_order6_cpu_s,'// &
            'da_integration_order8_cpu_s,state_compile_batch_eval_cpu_s,'// &
            'angle_build_compile_batch_eval_cpu_s,other_analysis_output_cpu_s,'// &
            'total_cpu_s'
    end subroutine open_outputs

    subroutine write_performance_row(id)
        integer,intent(in) :: id
        write(performance_unit,'(*(g0,:,","))') trim(MODEL_NAMES(id)), &
            N_SAMPLES,N_EPOCHS,timing_real(id),timing_da_order(:,id), &
            timing_state_batch(id),timing_angle_batch(id), &
            timing_analysis_other(id),timing_total(id)
        flush(performance_unit)
    end subroutine write_performance_row

    subroutine write_summary()
        real(DP) :: position_rms,velocity_rms,tangent_rms
        integer :: id
        do id=1,N_MODELS
            position_rms=sqrt(max(trace_matrix( &
                covariance_store(1:3,1:3,N_EPOCHS,id)),0.0_DP))
            velocity_rms=sqrt(max(trace_matrix( &
                covariance_store(4:6,4:6,N_EPOCHS,id)),0.0_DP))
            tangent_rms=sqrt(max(trace_matrix( &
                tangent_cov_store(:,:,N_EPOCHS,id)),0.0_DP))
            write(summary_unit,'(*(g0,:,","))') trim(MODEL_NAMES(id)), &
                first_distorted_hour(id),max_eta_nl(id),max_eta_cov(id), &
                max_eta_mean(id),max_abs_skew(id),max_abs_kurt(id), &
                position_rms,pos_axis_store(1,N_EPOCHS,id), &
                pos_axis_store(3,N_EPOCHS,id),velocity_rms,tangent_rms, &
                tangent_axis_store(1,N_EPOCHS,id), &
                tangent_axis_store(2,N_EPOCHS,id),N_BOUNDARY*N_EPOCHS, &
                validation_failures_by_model(id)
        end do
    end subroutine write_summary

    subroutine write_history_header()
        write(history_unit,'(a)',advance='no') 'model,time_hours,da_order,'
        call write_vector_labels(history_unit,'nominal',6)
        call write_vector_labels(history_unit,'mean_i',6)
        call write_vector_labels(history_unit,'mean_rtn',6)
        call write_upper_labels(history_unit,'cov_i',6)
        call write_upper_labels(history_unit,'cov_rtn',6)
        call write_vector_labels(history_unit,'pos_sigma_i',3)
        call write_vector_labels(history_unit,'pos_sigma_rtn',3)
        call write_vector_labels(history_unit,'vel_sigma_i',3)
        call write_vector_labels(history_unit,'vel_sigma_rtn',3)
        call write_matrix_labels(history_unit,'pos_axis_i',3,3)
        call write_matrix_labels(history_unit,'pos_axis_rtn',3,3)
        call write_matrix_labels(history_unit,'vel_axis_i',3,3)
        call write_matrix_labels(history_unit,'vel_axis_rtn',3,3)
        call write_vector_labels(history_unit,'radec_mean',2)
        call write_upper_labels(history_unit,'radec_cov',2)
        call write_vector_labels(history_unit,'radec_sigma',2)
        call write_matrix_labels(history_unit,'radec_axis',2,2)
        call write_vector_labels(history_unit,'radec_sigma_arcsec',2)
        call write_vector_labels(history_unit,'tangent_mean',2)
        call write_upper_labels(history_unit,'tangent_cov',2)
        call write_vector_labels(history_unit,'tangent_sigma',2)
        call write_matrix_labels(history_unit,'tangent_axis',2,2)
        call write_vector_labels(history_unit,'tangent_sigma_arcsec',2)
        write(history_unit,'(a)') &
            'position_rms_km,velocity_rms_kms,tangent_rms_rad,'// &
            'eta_nl_state,eta_cov_state,eta_mean_state,eta_nl_tangent,'// &
            'eta_cov_tangent,eta_mean_tangent,max_state_pc_skew,'// &
            'max_state_pc_excess_kurtosis,max_tangent_pc_skew,'// &
            'max_tangent_pc_excess_kurtosis,state_effective_rank,'// &
            'tangent_effective_rank,distorted'
    end subroutine write_history_header

    subroutine write_vector_labels(unit,prefix,n)
        integer,intent(in) :: unit,n
        character(len=*),intent(in) :: prefix
        integer :: i
        do i=1,n
            write(unit,'(a,"_",i0,",")',advance='no') trim(prefix),i
        end do
    end subroutine write_vector_labels

    subroutine write_upper_labels(unit,prefix,n)
        integer,intent(in) :: unit,n
        character(len=*),intent(in) :: prefix
        integer :: i,j
        do i=1,n
            do j=i,n
                write(unit,'(a,"_",i0,"_",i0,",")',advance='no') trim(prefix),i,j
            end do
        end do
    end subroutine write_upper_labels

    subroutine write_matrix_labels(unit,prefix,nrow,ncol)
        integer,intent(in) :: unit,nrow,ncol
        character(len=*),intent(in) :: prefix
        integer :: i,j
        do j=1,ncol
            do i=1,nrow
                write(unit,'(a,"_",i0,"_",i0,",")',advance='no') trim(prefix),i,j
            end do
        end do
    end subroutine write_matrix_labels

    subroutine assert_close(actual,expected,tolerance,label)
        real(DP),intent(in) :: actual,expected,tolerance
        character(len=*),intent(in) :: label
        if(abs(actual-expected)>tolerance) then
            write(*,'(a,3(1x,es24.16))') 'FAIL: '//trim(label), &
                actual,expected,tolerance
            error stop 1
        end if
    end subroutine assert_close

    subroutine assert_true(condition,label)
        logical,intent(in) :: condition
        character(len=*),intent(in) :: label
        if(.not.condition) then
            write(*,'(a)') 'FAIL: '//trim(label)
            error stop 1
        end if
    end subroutine assert_true

end program test_halo_srp_initial_state_uncertainty_15day
