!> @file test_halo_srp_acceleration_history_15day.f90
!! @brief Export SRP acceleration along four independently propagated Halo orbits.
!!
!! Run from the project root:
!!   fpm test --target test_halo_srp_acceleration_history_15day
!!
!! The optional command-line argument --opm <file> selects an external OPM.
!! With no argument this test uses the repository L1Halo-1 Earth-relative
!! J2000 state, propagating it independently with cannonball,
!! Sun-pointing box-wing, Earth-pointing box-wing, and Moon-pointing box-wing.
!! The analytic pointing law uses +Z of the body frame as its primary axis and
!! the orbit normal as its roll reference. Each model's acceleration is sampled
!! on ITS OWN trajectory, not on a common reference trajectory.
!!
!! This is a nominal Real experiment: no initial-state error, SRP scale error,
!! attitude bias, array-angle uncertainty, Gaussian sampling, or DA arithmetic.
!! The simplified production SRP model does not activate eclipse transitions.
!! Other force contributions are inherited from CONFIG_FILE; only SRP model
!! and spacecraft geometry are explicitly selected here.
!!
!! Integrator inputs are nondimensional. Output state is in km and km/s; SRP
!! acceleration is in km/s2. The SRP evaluator needs ABSOLUTE TDB seconds,
!! whereas RKF78 uses elapsed time divided by config%TU. Hourly segmentation
!! gives exact 0..360 h sampling endpoints without state interpolation.
!! RTN is an instantaneous vector rotation using each orbit's own r and v,
!! not a rotating-frame state derivative. The CSV contains SRP ONLY, not total
!! gravitational/non-gravitational acceleration.
program test_halo_srp_acceleration_history_15day
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config, validate_config
    use pod_data_format_module, only: load_initial_opm
    use pod_force_model_module, only: set_propagation_epoch, &
        clear_srp_scale_error, compute_solar_radiation_pressure, AU_KM
    use pod_integrator_module, only: adaptive_step_integrate, METHOD_RKF78
    use pod_frame_module, only: build_rtn_rotation, transform_vector_to_rtn
    use pod_spice, only: get_body_state
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    integer, parameter :: N_MODELS=4, N_HOURS=360, N_EPOCHS=N_HOURS+1
    integer, parameter :: MAX_STEPS_PER_HOUR=1000
    real(DP), parameter :: HOUR_S=3600.0_DP
    real(DP), parameter :: REL_TOL=1.0e-12_DP, ABS_TOL=1.0e-12_DP
    real(DP), parameter :: DT_MIN_S=1.0e-6_DP, DT_MAX_S=3600.0_DP
    real(DP), parameter :: ENDPOINT_TOL_S=1.0e-6_DP
    real(DP), parameter :: NORM_REL_TOL=1.0e-12_DP
    real(DP), parameter :: MASS_KG=1200.0_DP
    real(DP), parameter :: PRESSURE_1AU=1367.0_DP/299792458.0_DP
    character(len=*), parameter :: CONFIG_FILE='config/config.txt'
    character(len=*), parameter :: DEFAULT_OPM='OPM/L1Halo-1/L1Halo-1_init.opm.json'
    character(len=*), parameter :: OUTPUT_DIR_BASE='SRP/260915_srp_acceleration_history'
    character(len=32), parameter :: MODEL_NAMES(N_MODELS)=[character(len=32) :: &
        'cannonball','box_wing_sun','box_wing_earth','box_wing_moon']
    character(len=*), parameter :: CSV_HEADER= &
        'model,time_hours,epoch_tdb_s,x_km,y_km,z_km,vx_km_s,vy_km_s,vz_km_s,'// &
        'srp_ax_j2000_km_s2,srp_ay_j2000_km_s2,srp_az_j2000_km_s2,'// &
        'srp_norm_km_s2,srp_ar_rtn_km_s2,srp_at_rtn_km_s2,srp_an_rtn_km_s2'

    real(DP) :: epoch0,state0(6),covariance0(6,6),initial_nd(6)
    real(DP) :: cpu_start,cpu_end
    integer :: output_unit,ios,exit_status,model_id,record_count
    character(len=96) :: context='initialization'
    character(len=512) :: target_opm=DEFAULT_OPM
    character(len=512) :: output_dir=OUTPUT_DIR_BASE
    character(len=512) :: output_csv
    logical :: external_opm,opm_exists

    call cpu_time(cpu_start)
    call parse_command_line(external_opm)
    if(external_opm) output_dir=OUTPUT_DIR_BASE//'_'//trim(orbit_tag(target_opm))
    output_csv=trim(output_dir)//'/srp_acceleration_history.csv'
    inquire(file=trim(target_opm),exist=opm_exists)
    call assert_true(opm_exists,'OPM file does not exist: '//trim(target_opm))
    call pod_engine_init(CONFIG_FILE)
    call assert_true(trim(config%reference_frame)=='J2000','requires J2000 configuration')
    call assert_true(trim(config%central_body)=='EARTH','requires Earth-relative configuration')
    ! The OPM covariance is intentionally ignored: all four states are nominal.
    call load_initial_opm(trim(target_opm),epoch0,state0,covariance0)
    call assert_true(all(ieee_is_finite(state0)),'initial state must be finite')
    call set_propagation_epoch(epoch0)
    call clear_srp_scale_error()
    initial_nd(1:3)=state0(1:3)/config%LU
    initial_nd(4:6)=state0(4:6)/config%VU
    call configure_geometry()

    call execute_command_line('mkdir -p '//trim(output_dir),exitstat=exit_status)
    call assert_true(exit_status==0,'cannot create output directory')
    open(newunit=output_unit,file=trim(output_csv),status='replace',action='write',iostat=ios)
    call assert_true(ios==0,'cannot open output CSV')
    write(output_unit,'(a)',iostat=ios) CSV_HEADER
    call assert_true(ios==0,'cannot write CSV header')
    record_count=0
    write(*,'(a)') '15-day nominal SRP acceleration history: '//trim(target_opm)
    do model_id=1,N_MODELS
        context=trim(MODEL_NAMES(model_id))
        call configure_model(model_id)
        call propagate_and_export(model_id)
        flush(output_unit)
    end do
    close(output_unit,iostat=ios)
    call assert_true(ios==0,'cannot close CSV')
    context='CSV readback'
    call assert_true(record_count==N_MODELS*N_EPOCHS,'unexpected exported record count')
    call verify_written_csv()
    call clear_srp_scale_error()
    call cpu_time(cpu_end)
    write(*,'(a,f12.3,a)') 'Total CPU time: ',cpu_end-cpu_start,' s'
    write(*,'(a)') 'Output: '//trim(output_csv)
    write(*,'(a,i0,a)') 'PASS: ',record_count,' SRP acceleration records exported and read back.'

contains

    !> Parse the optional external OPM path while preserving the legacy default.
    subroutine parse_command_line(has_external_opm)
        logical,intent(out) :: has_external_opm
        character(len=512) :: argument
        integer :: argument_count,index

        has_external_opm=.false.
        argument_count=command_argument_count()
        index=1
        do while(index<=argument_count)
            call get_command_argument(index,argument)
            select case(trim(argument))
            case('-h','--help')
                call print_usage()
                stop 0
            case('-opm','--opm')
                if(index==argument_count) then
                    write(*,'(a,a)') 'ERROR: missing OPM path after ',trim(argument)
                    error stop 2
                end if
                call get_command_argument(index+1,target_opm)
                if(len_trim(target_opm)==0 .or. target_opm(1:1)=='-') then
                    write(*,'(a,a)') 'ERROR: missing OPM path after ',trim(argument)
                    error stop 2
                end if
                has_external_opm=.true.
                index=index+1
            case default
                write(*,'(a,a)') 'ERROR: unknown argument: ',trim(argument)
                call print_usage()
                error stop 2
            end select
            index=index+1
        end do
    end subroutine parse_command_line

    subroutine print_usage()
        write(*,'(a)') 'Usage: fpm test --target test_halo_srp_acceleration_history_15day -- [--opm <file>]'
        write(*,'(a)') '  --opm <file>  Earth-relative J2000 OPM; default: '//DEFAULT_OPM
    end subroutine print_usage

    !> Return a filesystem-safe orbit tag from an OPM basename.
    function orbit_tag(path) result(tag)
        character(len=*),intent(in) :: path
        character(len=256) :: tag,base
        integer :: slash,backslash,dot,i,code

        base=trim(path)
        slash=index(base,'/',back=.true.)
        backslash=index(base,achar(92),back=.true.)
        slash=max(slash,backslash)
        if(slash>0) base=base(slash+1:)
        dot=index(base,'.opm.json',back=.true.)
        if(dot>1) then
            tag=base(1:dot-1)
        else
            dot=index(base,'.json',back=.true.)
            if(dot>1) then
                tag=base(1:dot-1)
            else
                tag=base
            end if
        end if
        do i=1,len_trim(tag)
            code=iachar(tag(i:i))
            if(.not.((code>=iachar('a') .and. code<=iachar('z')) .or. &
                     (code>=iachar('A') .and. code<=iachar('Z')) .or. &
                     (code>=iachar('0') .and. code<=iachar('9')) .or. &
                     tag(i:i)=='-' .or. tag(i:i)=='_')) tag(i:i)='_'
        end do
        if(len_trim(tag)==0) tag='orbit'
    end function orbit_tag

    !> Match the latest geometry used by the existing deterministic Halo tests.
    !! Optical arrays are absorption, specular reflection, diffuse reflection.
    !! The cannonball law has one fixed reference area rather than projected
    !! box/array areas; its config-driven values are checked at the first epoch.
    subroutine configure_geometry()
        config%use_srp=.true.
        config%srp_mass_kg=MASS_KG
        config%srp_box_dimensions_m=[2.0_DP,3.0_DP,4.0_DP]
        config%srp_box_optical=[0.30_DP,0.40_DP,0.30_DP]
        config%srp_array_total_area_m2=24.0_DP
        config%srp_array_tracking_mode='single_axis'
        config%srp_array_hinge_axis_body=[0.0_DP,1.0_DP,0.0_DP]
        config%srp_array_reference_normal_body=[1.0_DP,0.0_DP,0.0_DP]
        config%srp_array_front_optical=[0.85_DP,0.08_DP,0.07_DP]
        config%srp_array_back_optical=[0.60_DP,0.20_DP,0.20_DP]
        config%srp_primary_axis_body=[0.0_DP,0.0_DP,1.0_DP]
        config%srp_secondary_axis_body=[0.0_DP,1.0_DP,0.0_DP]
        config%srp_roll_reference='orbit_normal'
        config%srp_pressure_1au_n_m2=PRESSURE_1AU
        config%srp_geometry_tolerance=1.0e-12_DP
        config%srp_scale_da_span=0.0_DP
        config%srp_attitude_bias_span_arcsec=0.0_DP
        config%srp_array_angle_span_deg=0.0_DP
    end subroutine configure_geometry

    !> Only the SRP model/pointing target varies between the four cases.
    subroutine configure_model(id)
        integer,intent(in) :: id
        config%srp_model='box_wing'
        select case(id)
        case(1)
            config%srp_model='cannonball'
            config%srp_attitude_mode='sun'
        case(2)
            config%srp_attitude_mode='sun'
        case(3)
            config%srp_attitude_mode='earth'
        case(4)
            config%srp_attitude_mode='moon'
        case default
            call assert_true(.false.,'invalid model identifier')
        end select
        call assert_true(validate_config(),'invalid SRP configuration')
    end subroutine configure_model

    !> Start afresh for each model, advance to exact hourly endpoints, and export.
    subroutine propagate_and_export(id)
        integer,intent(in) :: id
        real(DP) :: current_nd(6),physical(6),a_i(3),a_rtn(3),rotation(3,3)
        real(DP) :: absolute_epoch,start_time,end_time,a_norm,norm_min,norm_max,norm_sum
        real(DP),allocatable :: times(:),states(:,:)
        integer :: hour,n_steps,frame_status,write_status
        character(len=256) :: frame_message

        current_nd=initial_nd
        norm_min=huge(1.0_DP)
        norm_max=0.0_DP
        norm_sum=0.0_DP
        do hour=0,N_HOURS
            write(context,'(a," @ ",i0," h")') trim(MODEL_NAMES(id)),hour
            physical(1:3)=current_nd(1:3)*config%LU
            physical(4:6)=current_nd(4:6)*config%VU
            absolute_epoch=epoch0+real(hour,DP)*HOUR_S
            call assert_true(all(ieee_is_finite(physical)),'non-finite propagated state')
            ! No optional ballistic overrides: use exactly the law called by
            ! compute_acceleration inside the production Real integrator.
            call compute_solar_radiation_pressure(physical(1:3),absolute_epoch,a_i, &
                                                  velocity=physical(4:6))
            call assert_true(all(ieee_is_finite(a_i)),'non-finite SRP acceleration')
            a_norm=norm2(a_i)
            call assert_true(a_norm>0.0_DP,'unexpected zero sunlight acceleration')
            if(id==1 .and. hour==0) call verify_cannonball_configuration(physical,a_i)
            call build_rtn_rotation(physical(1:3),physical(4:6),rotation, &
                                    frame_status,frame_message)
            call assert_true(frame_status==0,'undefined RTN frame: '//trim(frame_message))
            call transform_vector_to_rtn(a_i,rotation,a_rtn)
            call assert_true(all(ieee_is_finite(a_rtn)),'non-finite RTN acceleration')
            call assert_true(abs(norm2(a_rtn)-a_norm)<=NORM_REL_TOL*a_norm, &
                             'RTN rotation changed acceleration norm')
            write(output_unit,'(*(g0,:,","))',iostat=write_status) &
                trim(MODEL_NAMES(id)),hour,absolute_epoch,physical,a_i,a_norm,a_rtn
            call assert_true(write_status==0,'cannot write acceleration row')
            record_count=record_count+1
            norm_min=min(norm_min,a_norm)
            norm_max=max(norm_max,a_norm)
            norm_sum=norm_sum+a_norm
            if(hour==N_HOURS) exit

            start_time=real(hour,DP)*HOUR_S/config%TU
            end_time=real(hour+1,DP)*HOUR_S/config%TU
            call adaptive_step_integrate(current_nd,start_time,end_time,METHOD_RKF78, &
                times,states,n_steps,max_steps_in=MAX_STEPS_PER_HOUR, &
                rel_tol_in=REL_TOL,abs_tol_in=ABS_TOL, &
                dt_min_in=DT_MIN_S,dt_max_in=DT_MAX_S)
            call assert_true(n_steps>1 .and. n_steps<=MAX_STEPS_PER_HOUR,'invalid step count')
            call assert_true(abs(times(n_steps)*config%TU-real(hour+1,DP)*HOUR_S) &
                             <=ENDPOINT_TOL_S,'integration did not reach hourly endpoint')
            current_nd=states(n_steps,:)
            deallocate(times,states)
        end do
        call assert_true(.not.allocated(times) .and. .not.allocated(states),'history cleanup failed')
        write(*,'(a)') '  '//trim(MODEL_NAMES(id))
        write(*,'(a,3(1x,es24.16))') '    final r [km]   :',physical(1:3)
        write(*,'(a,3(1x,es24.16))') '    final v [km/s] :',physical(4:6)
        write(*,'(a,3(1x,es24.16))') '    |SRP| min/max/mean [km/s2]:', &
            norm_min,norm_max,norm_sum/real(N_EPOCHS,DP)
    end subroutine propagate_and_export

    !> Ensure the recorded cannonball force uses the stated Cr and A/m defaults.
    subroutine verify_cannonball_configuration(physical,acceleration)
        real(DP),intent(in) :: physical(6),acceleration(3)
        real(DP) :: sun_position(3),sun_velocity(3),relative(3),distance,expected(3)
        call get_body_state('SUN',epoch0,'EARTH',sun_position,sun_velocity)
        relative=physical(1:3)-sun_position
        distance=norm2(relative)
        expected=config%srp_cannonball_cr* &
                 (config%srp_cannonball_area_m2/config%srp_mass_kg)*PRESSURE_1AU* &
                 (AU_KM/distance)**2*(relative/distance)*1.0e-3_DP
        call assert_true(norm2(acceleration-expected)<=NORM_REL_TOL*norm2(expected), &
                         'production cannonball defaults differ from stated parameters')
    end subroutine verify_cannonball_configuration

    !> Validate the actual written file, not just the number of write calls.
    subroutine verify_written_csv()
        character(len=2048) :: line
        character(len=32) :: name
        real(DP) :: values(14)
        integer :: unit,read_status,row,expected_model,expected_hour,hour
        open(newunit=unit,file=trim(output_csv),status='old',action='read',iostat=read_status)
        call assert_true(read_status==0,'cannot reopen CSV')
        read(unit,'(a)',iostat=read_status) line
        call assert_true(read_status==0,'cannot read CSV header')
        call assert_true(trim(line)==CSV_HEADER,'incorrect CSV header')
        row=0
        do
            read(unit,'(a)',iostat=read_status) line
            if(read_status<0) exit
            call assert_true(read_status==0,'cannot read CSV row')
            call assert_true(row<N_MODELS*N_EPOCHS,'too many CSV rows')
            read(line,*,iostat=read_status) name,hour,values
            call assert_true(read_status==0,'cannot parse CSV row')
            expected_model=1+row/N_EPOCHS
            expected_hour=mod(row,N_EPOCHS)
            call assert_true(trim(name)==trim(MODEL_NAMES(expected_model)),'incorrect model order')
            call assert_true(hour==expected_hour,'incorrect hourly sequence')
            call assert_true(all(ieee_is_finite(values)),'non-finite numeric CSV field')
            call assert_true(abs(values(1)-(epoch0+real(hour,DP)*HOUR_S))<=ENDPOINT_TOL_S, &
                             'incorrect absolute TDB epoch')
            call assert_true(values(11)>0.0_DP,'invalid stored SRP magnitude')
            call assert_true(abs(norm2(values(8:10))-values(11))<= &
                             NORM_REL_TOL*values(11),'stored J2000 norm mismatch')
            call assert_true(abs(norm2(values(12:14))-values(11))<= &
                             NORM_REL_TOL*values(11),'stored RTN norm mismatch')
            row=row+1
        end do
        close(unit)
        call assert_true(row==N_MODELS*N_EPOCHS,'incorrect CSV data-row count')
    end subroutine verify_written_csv

    subroutine assert_true(condition,label)
        logical,intent(in) :: condition
        character(len=*),intent(in) :: label
        if(.not.condition) then
            write(*,'(a)') 'FAIL ['//trim(context)//']: '//trim(label)
            error stop 1
        end if
    end subroutine assert_true
end program test_halo_srp_acceleration_history_15day
