!> P1: all SRP laws are evaluated on one deterministic Cannonball trajectory.
program run_boxwing_mode_extract
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config, validate_config
    use pod_data_format_module, only: load_initial_opm
    use pod_force_model_module, only: set_propagation_epoch, clear_srp_scale_error, &
        clear_srp_ballistic_parameters, compute_solar_radiation_pressure
    use pod_integrator_module, only: adaptive_step_integrate, METHOD_RKF78
    use pod_frame_module, only: build_rtn_rotation, transform_vector_to_rtn
    use pod_spacecraft_geometry, only: reset_srp_geometry_diagnostics, &
        get_srp_geometry_diagnostics, srp_geometry_diagnostics_type
    use pod_spice, only: spice_init, str2et
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    real(DP), parameter :: REL_TOL=1.0e-12_DP, ABS_TOL=1.0e-12_DP
    real(DP), parameter :: DT_MIN_S=1.0e-6_DP, DT_MAX_S=3600.0_DP
    integer, parameter :: MAX_STEPS=10000
    character(len=*), parameter :: CSV_HEADER= &
        'time_hours,epoch_tdb_s,ref_x_km,ref_y_km,ref_z_km,ref_vx_km_s,ref_vy_km_s,ref_vz_km_s,'// &
        'a_cb_x_km_s2,a_cb_y_km_s2,a_cb_z_km_s2,a_cb_norm_km_s2,'// &
        'a_bwe_x_km_s2,a_bwe_y_km_s2,a_bwe_z_km_s2,a_bwm_x_km_s2,a_bwm_y_km_s2,a_bwm_z_km_s2,'// &
        'eE_x,eE_y,eE_z,eE_norm,eM_x,eM_y,eM_z,eM_norm'
    character(len=1024) :: opm_file='OPM/DRO-2/DRO-2_init.opm.json'
    character(len=1024) :: config_file='config/config_box_wing_earth.txt', prefix='', csv_file, meta_file
    character(len=1024) :: working_dir, executable_file
    character(len=256) :: message
    real(DP) :: duration_hours=360.0_DP, sample_hours=1.0_DP
    real(DP) :: epoch0, state0(6), covariance(6,6), current_nd(6), physical(6)
    real(DP) :: time_hours, absolute_epoch, cb(3), bwe(3), bwm(3), cb_norm
    real(DP) :: earth_mode(3), moon_mode(3), earth_norm, moon_norm, rotation(3,3), rotated(3)
    real(DP) :: max_earth=-1.0_DP, max_moon=-1.0_DP, earth_time=0.0_DP, moon_time=0.0_DP
    real(DP) :: start_nd, end_nd, intervals
    real(DP), allocatable :: times(:), states(:,:)
    integer :: intervals_count, index, csv_unit, metadata_unit, io, n_steps, frame_status
    logical :: exists
    type(srp_geometry_diagnostics_type) :: diagnostics

    call parse_arguments()
    call get_environment_variable('PWD',working_dir,status=io)
    call require(io==0 .and. len_trim(working_dir)>0, 'PWD must identify the project directory')
    call get_command_argument(0,executable_file)
    csv_file=trim(prefix)//'.csv'
    meta_file=trim(prefix)//'_meta.json'
    inquire(file=trim(csv_file),exist=exists)
    call require(.not.exists,'output CSV already exists; choose a new staging prefix')
    inquire(file=trim(meta_file),exist=exists)
    call require(.not.exists,'output metadata already exists; choose a new staging prefix')
    inquire(file=trim(opm_file),exist=exists)
    call require(exists,'OPM file does not exist')
    inquire(file=trim(config_file),exist=exists)
    call require(exists,'configuration file does not exist')
    call pod_engine_init(trim(config_file))
    call require(trim(config%reference_frame)=='J2000','requires J2000')
    call require(trim(config%central_body)=='EARTH','requires Earth center')
    call clear_srp_scale_error()
    call clear_srp_ballistic_parameters()
    config%use_srp=.true.
    config%srp_scale_da_span=0.0_DP
    config%srp_attitude_bias_span_arcsec=0.0_DP
    config%srp_array_angle_span_deg=0.0_DP
    ! Validate the actual geometry even when the input config selects Cannonball.
    config%srp_model='box_wing'
    config%srp_attitude_mode='earth'
    call require(validate_config(),'invalid Earth-pointing Box-wing configuration')
    config%srp_attitude_mode='moon'
    call require(validate_config(),'invalid Moon-pointing Box-wing configuration')
    config%srp_model='cannonball'
    call load_initial_opm(trim(opm_file),epoch0,state0,covariance)
    call require(ieee_is_finite(epoch0) .and. all(ieee_is_finite(state0)),'initial state must be finite')
    call set_propagation_epoch(epoch0)
    call reset_srp_geometry_diagnostics()
    current_nd(1:3)=state0(1:3)/config%LU
    current_nd(4:6)=state0(4:6)/config%VU
    open(newunit=csv_unit,file=trim(csv_file),status='new',action='write',iostat=io)
    call require(io==0,'cannot create staging CSV')
    write(csv_unit,'(A)',iostat=io) CSV_HEADER
    call require(io==0,'cannot write CSV header')
    intervals_count=nint(duration_hours/sample_hours)
    do index=0,intervals_count
        time_hours=real(index,DP)*sample_hours
        absolute_epoch=epoch0+3600.0_DP*time_hours
        physical(1:3)=current_nd(1:3)*config%LU
        physical(4:6)=current_nd(4:6)*config%VU
        call require(all(ieee_is_finite(physical)),'reference state is not finite')
        ! Never advance a Box-wing trajectory: only the instantaneous SRP selector changes.
        config%srp_model='cannonball'
        call compute_solar_radiation_pressure(physical(1:3),absolute_epoch,cb,velocity=physical(4:6))
        config%srp_model='box_wing'
        config%srp_attitude_mode='earth'
        call compute_solar_radiation_pressure(physical(1:3),absolute_epoch,bwe,velocity=physical(4:6))
        config%srp_attitude_mode='moon'
        call compute_solar_radiation_pressure(physical(1:3),absolute_epoch,bwm,velocity=physical(4:6))
        config%srp_model='cannonball'
        call require(all(ieee_is_finite(cb)) .and. all(ieee_is_finite(bwe)) .and. &
                     all(ieee_is_finite(bwm)),'SRP acceleration is not finite')
        cb_norm=norm2(cb)
        call require(cb_norm>0.0_DP,'Cannonball acceleration is zero')
        earth_mode=(bwe-cb)/cb_norm
        moon_mode=(bwm-cb)/cb_norm
        earth_norm=norm2(earth_mode)
        moon_norm=norm2(moon_mode)
        call require(all(ieee_is_finite(earth_mode)) .and. all(ieee_is_finite(moon_mode)), &
                     'mode is not finite')
        call build_rtn_rotation(physical(1:3),physical(4:6),rotation,frame_status,message)
        call require(frame_status==0,'undefined RTN frame: '//trim(message))
        call transform_vector_to_rtn(earth_mode,rotation,rotated)
        call require(abs(norm2(rotated)-earth_norm)<=1.0e-12_DP*earth_norm,'Earth RTN norm changed')
        call transform_vector_to_rtn(moon_mode,rotation,rotated)
        call require(abs(norm2(rotated)-moon_norm)<=1.0e-12_DP*moon_norm,'Moon RTN norm changed')
        if(earth_norm>max_earth) then
            max_earth=earth_norm
            earth_time=time_hours
        end if
        if(moon_norm>max_moon) then
            max_moon=moon_norm
            moon_time=time_hours
        end if
        write(csv_unit,'(*(ES25.17E3,:,","))',iostat=io) &
            time_hours,absolute_epoch,physical,cb,cb_norm,bwe,bwm,earth_mode,earth_norm,moon_mode,moon_norm
        call require(io==0,'cannot write CSV data')
        if(mod(index,24)==0 .or. index==intervals_count) then
            flush(csv_unit)
            write(*,'(A,F8.2,A,2(1X,ES15.7))') 'P1 reference hour ',time_hours,'; mode norms:',earth_norm,moon_norm
            flush(6)
        end if
        if(index==intervals_count) exit
        start_nd=3600.0_DP*time_hours/config%TU
        end_nd=3600.0_DP*(time_hours+sample_hours)/config%TU
        call adaptive_step_integrate(current_nd,start_nd,end_nd,METHOD_RKF78,times,states,n_steps, &
            max_steps_in=MAX_STEPS,rel_tol_in=REL_TOL,abs_tol_in=ABS_TOL, &
            dt_min_in=DT_MIN_S,dt_max_in=DT_MAX_S)
        call require(n_steps>1 .and. n_steps<=MAX_STEPS,'invalid integration step count')
        call require(abs(times(n_steps)*config%TU-3600.0_DP*(time_hours+sample_hours))<=1.0e-6_DP, &
                     'integrator did not reach sampling endpoint')
        current_nd=states(n_steps,:)
        deallocate(times,states)
    end do
    close(csv_unit,iostat=io)
    call require(io==0,'cannot close staging CSV')
    call get_srp_geometry_diagnostics(diagnostics)
    call write_metadata()
    call clear_srp_scale_error()
    call clear_srp_ballistic_parameters()
    write(*,'(A,ES24.16,A,F8.2)') 'max ||e_E|| = ',max_earth,' at hour ',earth_time
    write(*,'(A,ES24.16,A,F8.2)') 'max ||e_M|| = ',max_moon,' at hour ',moon_time
    write(*,'(A,L1)') '30% physical gate: ',max_earth<=0.3_DP .and. max_moon<=0.3_DP
    write(*,'(A)') 'Staging data require Python audit/finalization: '//trim(csv_file)

contains

    subroutine require(condition,text)
        logical,intent(in) :: condition
        character(len=*),intent(in) :: text
        if(.not.condition) then
            call clear_srp_scale_error()
            call clear_srp_ballistic_parameters()
            write(*,'(A)') 'P1 extraction failed: '//trim(text)
            error stop 1
        end if
    end subroutine require

    subroutine parse_arguments()
        character(len=1024) :: argument,value
        integer :: arg_index,parse_status
        arg_index=1
        if(command_argument_count()==2) then
            call get_command_argument(1,argument)
            if(trim(argument)=='--utc-to-tdb') then
                call get_command_argument(2,value)
                call spice_init()
                call str2et(trim(value),epoch0)
                call require(ieee_is_finite(epoch0),'SPICE epoch is not finite')
                write(*,'(A,ES25.17E3)') 'P1_EPOCH_TDB=',epoch0
                stop 0
            end if
        end if
        do while(arg_index<=command_argument_count())
            call get_command_argument(arg_index,argument)
            if(trim(argument)=='--help' .or. trim(argument)=='-h') then
                write(*,'(A)') 'run_boxwing_mode_extract -o PREFIX [-opm FILE -cfg FILE --hours 360 --sample-hours 1]'
                write(*,'(A)') 'Writes PREFIX.csv and PREFIX_meta.json; no existing output is replaced.'
                write(*,'(A)') '--utc-to-tdb UTC_STRING performs the authoritative SPICE epoch audit.'
                stop 0
            end if
            call require(arg_index<command_argument_count(),'option requires a value: '//trim(argument))
            call get_command_argument(arg_index+1,value)
            select case(trim(argument))
            case('-opm')
                opm_file=value
            case('-cfg')
                config_file=value
            case('-o')
                prefix=value
            case('--hours','--sample-hours')
                if(trim(argument)=='--hours') then
                    read(value,*,iostat=parse_status) duration_hours
                else
                    read(value,*,iostat=parse_status) sample_hours
                end if
                call require(parse_status==0,'invalid numeric option: '//trim(argument))
            case default
                call require(.false.,'unknown option: '//trim(argument))
            end select
            arg_index=arg_index+2
        end do
        call require(len_trim(prefix)>0,'-o staging prefix is required')
        call require(ieee_is_finite(duration_hours) .and. ieee_is_finite(sample_hours),'duration/cadence must be finite')
        call require(duration_hours>0.0_DP .and. sample_hours>0.0_DP,'duration/cadence must be positive')
        intervals=duration_hours/sample_hours
        call require(intervals>=1.0_DP .and. intervals<=1.0e6_DP,'invalid sampling interval count')
        call require(abs(real(nint(intervals),DP)*sample_hours-duration_hours)<=1.0e-12_DP, &
                     'duration must contain whole sampling intervals')
    end subroutine parse_arguments

    function absolute_path(path) result(absolute)
        character(len=*),intent(in) :: path
        character(len=:),allocatable :: absolute
        if(path(1:1)=='/') then
            absolute=trim(path)
        else
            absolute=trim(working_dir)//'/'//trim(path)
        end if
    end function absolute_path

    function quoted(value) result(text)
        character(len=*),intent(in) :: value
        character(len=:),allocatable :: text
        integer :: character_index,code
        text=achar(34)
        do character_index=1,len_trim(value)
            code=iachar(value(character_index:character_index))
            select case(code)
            case(34,92)
                text=text//achar(92)//value(character_index:character_index)
            case(10)
                text=text//achar(92)//'n'
            case(13)
                text=text//achar(92)//'r'
            case(9)
                text=text//achar(92)//'t'
            case default
                call require(code>=32,'control character in metadata string')
                text=text//value(character_index:character_index)
            end select
        end do
        text=text//achar(34)
    end function quoted

    subroutine json_string(key,value,comma)
        character(len=*),intent(in) :: key,value
        logical,intent(in),optional :: comma
        character(len=1) :: suffix
        suffix=','
        if(present(comma)) then
            if(.not.comma) suffix=' '
        end if
        write(metadata_unit,'(A)',iostat=io) '  '//quoted(key)//': '//quoted(value)//suffix
        call require(io==0,'cannot write metadata string')
    end subroutine json_string

    subroutine json_real(key,value,comma)
        character(len=*),intent(in) :: key
        real(DP),intent(in) :: value
        logical,intent(in),optional :: comma
        character(len=1) :: suffix
        suffix=','
        if(present(comma)) then
            if(.not.comma) suffix=' '
        end if
        call require(ieee_is_finite(value),'metadata number is not finite: '//trim(key))
        write(metadata_unit,'(A,ES25.17E3,A)',iostat=io) '  '//quoted(key)//': ',value,suffix
        call require(io==0,'cannot write metadata number')
    end subroutine json_real

    subroutine json_vector(key,value)
        character(len=*),intent(in) :: key
        real(DP),intent(in) :: value(3)
        call require(all(ieee_is_finite(value)),'metadata vector is not finite')
        write(metadata_unit,'(A,ES25.17E3,A,ES25.17E3,A,ES25.17E3,A)',iostat=io) &
            '  '//quoted(key)//': [',value(1),',',value(2),',',value(3),'],'
        call require(io==0,'cannot write metadata vector')
    end subroutine json_vector

    subroutine write_metadata()
        open(newunit=metadata_unit,file=trim(meta_file),status='new',action='write',iostat=io)
        call require(io==0,'cannot create staging metadata')
        write(metadata_unit,'(A)') '{'
        call json_string('schema_version','boxwing_srp_modes_v1')
        write(metadata_unit,'(A)') '  "audit_complete": false,'
        call json_string('frame','J2000')
        call json_string('center','EARTH')
        call json_string('timescale','TDB')
        call json_string('illumination','full_sun')
        call json_string('reference_model','cannonball')
        call json_string('evaluation_contract','common_reference_state')
        write(metadata_unit,'(A)') '  "units": {"position":"km","velocity":"km/s",'// &
            '"acceleration":"km/s2","modes":"dimensionless"},'
        call json_string('opm_file',absolute_path(opm_file))
        call json_string('config_file',absolute_path(config_file))
        call json_string('executable_file',absolute_path(executable_file))
        call json_real('epoch0_tdb_s',epoch0)
        call json_real('duration_hours',duration_hours)
        call json_real('sample_hours',sample_hours)
        write(metadata_unit,'(A,I0,A)') '  "row_count": ',intervals_count+1,','
        call json_real('max_eE_norm',max_earth)
        call json_real('max_eE_time_hours',earth_time)
        call json_real('max_eM_norm',max_moon)
        call json_real('max_eM_time_hours',moon_time)
        write(metadata_unit,'(A)') '  "integrator": {'
        call json_string('method','RKF78')
        call json_real('rel_tol',REL_TOL)
        call json_real('abs_tol',ABS_TOL)
        call json_real('dt_min_s',DT_MIN_S)
        call json_real('dt_max_s',DT_MAX_S,.false.)
        write(metadata_unit,'(A)') '  },'
        write(metadata_unit,'(A)') '  "physical_parameters": {'
        call json_real('srp_cannonball_cr',config%srp_cannonball_cr)
        call json_real('srp_cannonball_area_m2',config%srp_cannonball_area_m2)
        call json_real('srp_mass_kg',config%srp_mass_kg)
        call json_vector('srp_box_dimensions_m',config%srp_box_dimensions_m)
        call json_vector('srp_box_optical',config%srp_box_optical)
        call json_real('srp_array_total_area_m2',config%srp_array_total_area_m2)
        call json_string('srp_array_tracking_mode',config%srp_array_tracking_mode)
        call json_vector('srp_array_hinge_axis_body',config%srp_array_hinge_axis_body)
        call json_vector('srp_array_reference_normal_body',config%srp_array_reference_normal_body)
        call json_vector('srp_array_front_optical',config%srp_array_front_optical)
        call json_vector('srp_array_back_optical',config%srp_array_back_optical)
        call json_vector('srp_primary_axis_body',config%srp_primary_axis_body)
        call json_vector('srp_secondary_axis_body',config%srp_secondary_axis_body)
        call json_string('srp_roll_reference',config%srp_roll_reference)
        call json_real('srp_pressure_1au_n_m2',config%srp_pressure_1au_n_m2)
        call json_real('srp_geometry_tolerance',config%srp_geometry_tolerance,.false.)
        write(metadata_unit,'(A)') '  },'
        write(metadata_unit,'(A)') '  "geometry_fallback_diagnostics": {'
        write(metadata_unit,'(A,I0,A)') '  "real_count": ',diagnostics%real_fallback_count,','
        call json_string('first_real_message',trim(diagnostics%first_real_message),.false.)
        write(metadata_unit,'(A)') '  }'
        write(metadata_unit,'(A)') '}'
        close(metadata_unit,iostat=io)
        call require(io==0,'cannot close metadata')
    end subroutine write_metadata
end program run_boxwing_mode_extract
