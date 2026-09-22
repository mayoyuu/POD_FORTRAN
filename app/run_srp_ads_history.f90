!> Hourly ADS error-domain history with independent real-SRP validation.
program run_srp_ads_history
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config
    use pod_data_format_module, only: load_initial_opm
    use pod_uq_ads_coordinates_module, only: ADS_COORD_WHITENED, ads_unit_to_physical
    use pod_uq_hfem_ads_module, only: hfem_ads_options_type, hfem_ads_history_type, &
        hfem_ads_history_init, hfem_ads_history_advance, hfem_ads_history_evaluate, &
        hfem_ads_history_destroy
    use pod_srp_ads_history_module, only: SHAPE_COUNT, VALIDATION_COUNT, &
        build_study_points, assess_accuracy, advance_real_probes, compute_shape_metrics, &
        validate_json_opm_fields
    use pod_srp_ads_snapshot_module, only: write_ads_patch_snapshot
    implicit none
    type(hfem_ads_options_type) :: options
    type(hfem_ads_history_type) :: history
    character(len=MAX_STRING_LEN) :: opm_file, prefix, config_file, arg, message, path
    real(DP) :: epoch0, nominal0(6), cov(6,6), eta_sigma
    real(DP) :: pos_tolerance, vel_tolerance, duration_hours, save_hours
    real(DP) :: t_seconds, prev_seconds, final_seconds, current_hour, previous_hour
    real(DP) :: max_validation_pos, max_validation_vel, worst_pos, worst_vel
    real(DP) :: lower(6), upper(6), pos_max, vel_max, pos_axes(3), vel_axes(3)
    real(DP) :: mapped(7), nominal(6), delta(6)
    real(DP), allocatable :: points(:,:), ads(:,:), truth(:,:), eta(:), truth_initial(:,:)
    integer, allocatable :: validation_ids(:)
    logical, allocatable :: found(:)
    logical :: has_opm, has_prefix, has_eta, failed, any_failed
    integer :: argc, i, j, unit_history, unit_shape, unit_initial
    integer :: status, step_index, steps, worst_id, first_worst_id, first_failed_step, previous_pass_step
    integer :: first_uncertified_step, depth_limited, last_pass_step
    character(len=16) :: hour_tag

    opm_file=''
    prefix=''
    config_file='config/config.txt'
    eta_sigma=-1.0_DP
    duration_hours=360.0_DP
    save_hours=1.0_DP
    pos_tolerance=0.1_DP
    vel_tolerance=1.0e-6_DP
    has_opm=.false.
    has_prefix=.false.
    has_eta=.false.
    argc=command_argument_count()
    i=1
    do while(i<=argc)
        call get_command_argument(i,arg)
        if(trim(arg)=='--help'.or.trim(arg)=='-h') then
            call usage()
            stop
        end if
        if(i==argc) call fatal('option requires a value: '//trim(arg))
        select case(trim(arg))
        case('-opm','--opm')
            call get_command_argument(i+1,opm_file)
            has_opm=.true.
        case('-o','--output')
            call get_command_argument(i+1,prefix)
            has_prefix=.true.
        case('--eta-sigma')
            call get_command_argument(i+1,arg)
            read(arg,*,iostat=status) eta_sigma
            if(status/=0) call fatal('invalid --eta-sigma')
            has_eta=.true.
        case('--hours')
            call get_command_argument(i+1,arg)
            read(arg,*,iostat=status) duration_hours
            if(status/=0) call fatal('invalid --hours')
        case('--days')
            call get_command_argument(i+1,arg)
            read(arg,*,iostat=status) duration_hours
            if(status/=0) call fatal('invalid --days')
            duration_hours=24.0_DP*duration_hours
        case('--save-hours')
            call get_command_argument(i+1,arg)
            read(arg,*,iostat=status) save_hours
            if(status/=0) call fatal('invalid --save-hours')
        case('--da-order')
            call get_command_argument(i+1,arg)
            read(arg,*,iostat=status) options%da_order
            if(status/=0) call fatal('invalid --da-order')
        case('--max-depth')
            call get_command_argument(i+1,arg)
            read(arg,*,iostat=status) options%max_split_depth
            if(status/=0) call fatal('invalid --max-depth')
        case('--pos-tol-km')
            call get_command_argument(i+1,arg)
            read(arg,*,iostat=status) pos_tolerance
            if(status/=0) call fatal('invalid --pos-tol-km')
        case('--vel-tol-kms')
            call get_command_argument(i+1,arg)
            read(arg,*,iostat=status) vel_tolerance
            if(status/=0) call fatal('invalid --vel-tol-kms')
        case('-cfg','--config')
            call get_command_argument(i+1,config_file)
        case default
            call fatal('unknown option: '//trim(arg))
        end select
        i=i+2
    end do
    if(.not.has_opm.or..not.has_prefix.or..not.has_eta) then
        call usage()
        call fatal('-opm, -o, and --eta-sigma are required')
    end if
    if(.not.ieee_is_finite(eta_sigma).or..not.ieee_is_finite(duration_hours).or. &
       .not.ieee_is_finite(save_hours).or..not.ieee_is_finite(pos_tolerance).or. &
       .not.ieee_is_finite(vel_tolerance)) call fatal('numeric arguments must be finite')
    if(eta_sigma<=0.0_DP.or.duration_hours<=0.0_DP.or.save_hours<=0.0_DP.or. &
       pos_tolerance<=0.0_DP.or.vel_tolerance<=0.0_DP) &
        call fatal('eta sigma, duration, cadence, and tolerances must be positive')
    if(options%da_order<1.or.options%max_split_depth<0) call fatal('invalid DA order or split depth')
    i=len_trim(opm_file)
    if(i<4) call fatal('input must be JSON .opm or .opm.json')
    if(opm_file(i-3:i)/='.opm') then
        if(i<9) call fatal('input must be JSON .opm or .opm.json')
        if(opm_file(i-8:i)/='.opm.json') call fatal('input must be JSON .opm or .opm.json')
    end if
    call validate_json_opm_fields(trim(opm_file),status,message)
    if(status/=0) call fatal('invalid JSON OPM: '//trim(message))
    call pod_engine_init(trim(config_file))
    config%use_srp=.true.
    call load_initial_opm(trim(opm_file),epoch0,nominal0,cov)
    options%coordinate_mode=ADS_COORD_WHITENED
    options%srp_sigma=eta_sigma
    options%error_tolerance=[pos_tolerance,pos_tolerance,pos_tolerance, &
        vel_tolerance,vel_tolerance,vel_tolerance]
    call hfem_ads_history_init(history,nominal0,cov,epoch0,options,status,message)
    if(status/=0) call fatal('ADS init: '//trim(message))
    call build_study_points(points,validation_ids)
    allocate(ads(6,SHAPE_COUNT),found(SHAPE_COUNT),truth(6,VALIDATION_COUNT), &
        truth_initial(6,VALIDATION_COUNT),eta(VALIDATION_COUNT))
    open(newunit=unit_initial,file=trim(prefix)//'_initial_points.csv', &
        status='replace',action='write',iostat=status)
    if(status/=0) call fatal('cannot open initial-points output')
    write(unit_initial,'(A)') 'point_id,u1,u2,u3,u4,u5,u6,u7,eta,x0_km,y0_km,z0_km,vx0_kms,vy0_kms,vz0_kms,validation'
    do j=1,SHAPE_COUNT
        call ads_unit_to_physical(history%coordinate_map,points(:,j),mapped(1:6),mapped(7),status)
        if(status/=0) call fatal('initial point mapping failed')
        mapped(1:6)=mapped(1:6)+nominal0
        write(unit_initial,'(I0,14(",",ES25.16E3),",",I0)') j,points(:,j), &
            mapped(7),mapped(1:6),merge(1,0,j<=VALIDATION_COUNT)
        if(j<=VALIDATION_COUNT) then
            truth(:,j)=mapped(1:6)
            eta(j)=mapped(7)
        end if
    end do
    close(unit_initial)
    truth_initial=truth
    open(newunit=unit_history,file=trim(prefix)//'_history.csv', &
        status='replace',action='write',iostat=status)
    if(status/=0) call fatal('cannot open history output')
    open(newunit=unit_shape,file=trim(prefix)//'_shape.csv', &
        status='replace',action='write',iostat=status)
    if(status/=0) call fatal('cannot open shape output')
    write(unit_history,'(A)') 'hour,epoch_et,nominal_x_km,nominal_y_km,nominal_z_km,'// &
        'nominal_vx_kms,nominal_vy_kms,nominal_vz_kms,'// &
        'lower_dx_km,lower_dy_km,lower_dz_km,lower_dvx_kms,lower_dvy_kms,lower_dvz_kms,'// &
        'upper_dx_km,upper_dy_km,upper_dz_km,upper_dvx_kms,upper_dvy_kms,upper_dvz_kms,'// &
        'pos_max_km,vel_max_kms,pos_axis1_km,pos_axis2_km,pos_axis3_km,'// &
        'vel_axis1_kms,vel_axis2_kms,vel_axis3_kms,patches,depth_limited,'// &
        'split_dim1,split_dim2,split_dim3,split_dim4,split_dim5,split_dim6,split_dim7,'// &
        'bfs_iterations,max_queue_size,step_elapsed_seconds,'// &
        'max_validation_pos_km,max_validation_vel_kms,worst_point_id,'// &
        'validation_failed,precision_unconfirmed'
    write(unit_shape,'(A)') 'hour,point_id,dx_km,dy_km,dz_km,dvx_kms,dvy_kms,dvz_kms'
    first_failed_step=-1
    previous_pass_step=-1
    first_uncertified_step=-1
    last_pass_step=-1
    final_seconds=duration_hours*3600.0_DP
    steps=ceiling(duration_hours/save_hours)
    prev_seconds=0.0_DP
    do step_index=0,steps
        t_seconds=min(final_seconds,real(step_index,DP)*save_hours*3600.0_DP)
        current_hour=t_seconds/3600.0_DP
        previous_hour=prev_seconds/3600.0_DP
        if(step_index>0) then
            call hfem_ads_history_advance(history,t_seconds,status,message)
            if(status/=0) call fatal('ADS advance: '//trim(message))
            call advance_real_probes(truth,eta,epoch0,prev_seconds,t_seconds, &
                options,status,message)
            if(status/=0) call fatal('real SRP advance: '//trim(message))
        end if
        call hfem_ads_history_evaluate(history,points,ads,found,status)
        if(status/=0.or..not.all(found)) call fatal('ADS point evaluation failed')
        nominal=truth(:,143)
        call compute_shape_metrics(ads,nominal,lower,upper,pos_max,vel_max, &
            pos_axes,vel_axes,status)
        if(status/=0) call fatal('shape diagnostics failed')
        call assess_accuracy(ads(:,validation_ids),truth,pos_tolerance, &
            vel_tolerance,max_validation_pos,max_validation_vel,worst_id,failed)
        depth_limited=history%stats%depth_limited_patches
        if(failed.and.first_failed_step<0) then
            first_failed_step=step_index
            previous_pass_step=last_pass_step
            first_worst_id=worst_id
            worst_pos=sqrt(sum((ads(1:3,worst_id)-truth(1:3,worst_id))**2))
            worst_vel=sqrt(sum((ads(4:6,worst_id)-truth(4:6,worst_id))**2))
        end if
        if(.not.failed) last_pass_step=step_index
        if(depth_limited>0.and.first_uncertified_step<0) first_uncertified_step=step_index
        write(unit_history,'(*(G0,:,","))') current_hour,epoch0+t_seconds,nominal, &
            lower,upper,pos_max,vel_max,pos_axes,vel_axes, &
            history%domain%n_patches,depth_limited,history%stats%split_counts, &
            history%stats%bfs_iterations,history%stats%max_queue_size, &
            history%stats%elapsed_seconds,max_validation_pos, &
            max_validation_vel,worst_id,merge(1,0,failed), &
            merge(1,0,depth_limited>0.or.first_failed_step>=0)
        do j=1,SHAPE_COUNT
            delta=ads(:,j)-nominal
            write(unit_shape,'(ES25.16E3,",",I0,6(",",ES25.16E3))') current_hour,j,delta
        end do
        flush(unit_history)
        flush(unit_shape)
        if(step_index==0.or.step_index==steps.or. &
           step_index==first_failed_step.or.step_index==first_uncertified_step) then
            call snapshot_current()
        end if
        if(step_index>0.and.(step_index==first_failed_step.or. &
            step_index==first_uncertified_step)) then
            write(hour_tag,'(I7.7)') nint(prev_seconds)
            path=trim(prefix)//'_patch_s'//trim(hour_tag)
            call write_ads_patch_snapshot(history,trim(path),.true.,status,message)
            if(status/=0) call fatal('previous ADS snapshot: '//trim(message))
        end if
        write(*,'(A,F8.2,A,I0,A,ES11.3,A,ES11.3)') 'hour ',current_hour, &
            ' patches ',history%domain%n_patches,' validation km ', &
            max_validation_pos,' km/s ',max_validation_vel
        prev_seconds=t_seconds
    end do
    close(unit_history)
    close(unit_shape)
    call write_report()
    call hfem_ads_history_destroy(history)
    if(first_failed_step>=0) then
        write(*,'(A,F8.2)') 'First observed validation failure at hour ', &
            min(duration_hours,real(first_failed_step,DP)*save_hours)
    else
        write(*,'(A)') 'No loss of accuracy observed at the validated times and points.'
    end if
    if(first_uncertified_step>=0) write(*,'(A,F8.2)') &
        'First depth-limited, tolerance-unconfirmed hour ', &
        min(duration_hours,real(first_uncertified_step,DP)*save_hours)
contains
    subroutine snapshot_current()
        write(hour_tag,'(I7.7)') nint(t_seconds)
        path=trim(prefix)//'_patch_s'//trim(hour_tag)
        call write_ads_patch_snapshot(history,trim(path),.false.,status,message)
        if(status/=0) call fatal('ADS snapshot: '//trim(message))
    end subroutine snapshot_current

    subroutine write_report()
        integer :: u, ios
        open(newunit=u,file=trim(prefix)//'_report.json',status='replace', &
            action='write',iostat=ios)
        if(ios/=0) call fatal('cannot write JSON report')
        write(u,'(A)') '{'
        write(u,'(A)') '  "method": "incremental ADS with independent real SRP probes",'
        write(u,'(A)') '  "sampled_ranges": true,'
        write(u,'(A)') '  "continuous_domain_certified": false,'
        write(u,'(A,I0,A)') '  "shape_point_count": ',SHAPE_COUNT,','
        write(u,'(A,I0,A)') '  "validation_point_count": ',VALIDATION_COUNT,','
        write(u,'(A,ES25.16E3,A)') '  "duration_hours": ',duration_hours,','
        write(u,'(A,ES25.16E3,A)') '  "save_hours": ',save_hours,','
        write(u,'(A,ES25.16E3,A)') '  "eta_sigma": ',eta_sigma,','
        write(u,'(A,ES25.16E3,A)') '  "position_tolerance_km": ',pos_tolerance,','
        write(u,'(A,ES25.16E3,A)') '  "velocity_tolerance_kms": ',vel_tolerance,','
        write(u,'(A,I0,A)') '  "da_order": ',options%da_order,','
        write(u,'(A,I0,A)') '  "max_split_depth": ',options%max_split_depth,','
        if(first_failed_step>=0) then
            write(u,'(A,ES25.16E3,A)') '  "first_failed_hour": ', &
                min(duration_hours,real(first_failed_step,DP)*save_hours),','
            if(previous_pass_step>=0) then
                write(u,'(A,ES25.16E3,A)') '  "previous_passing_hour": ', &
                    min(duration_hours,real(previous_pass_step,DP)*save_hours),','
            else
                write(u,'(A)') '  "previous_passing_hour": null,'
            end if
            write(u,'(A,I0,A)') '  "worst_point_id": ',first_worst_id,','
            write(u,'(A,ES25.16E3,A)') '  "worst_point_position_error_km": ',worst_pos,','
            write(u,'(A,ES25.16E3,A)') '  "worst_point_velocity_error_kms": ',worst_vel,','
            write(u,'(A)') '  "accuracy_status": "observed_failure",'
        else
            write(u,'(A)') '  "first_failed_hour": null,'
            write(u,'(A)') '  "previous_passing_hour": null,'
            write(u,'(A)') '  "worst_point_id": null,'
            write(u,'(A)') '  "worst_point_position_error_km": null,'
            write(u,'(A)') '  "worst_point_velocity_error_kms": null,'
            write(u,'(A)') '  "accuracy_status": "no_failure_observed_at_validated_times_and_points",'
        end if
        if(first_uncertified_step>=0) then
            write(u,'(A,ES25.16E3)') '  "first_unconfirmed_hour": ', &
                min(duration_hours,real(first_uncertified_step,DP)*save_hours)
        else
            write(u,'(A)') '  "first_unconfirmed_hour": null'
        end if
        write(u,'(A)') '}'
        close(u)
    end subroutine write_report

    subroutine usage()
        write(*,'(A)') 'Usage: run_srp_ads_history -opm INPUT.opm.json -o OUTPUT_PREFIX '// &
            '--eta-sigma VALUE [--days 15 | --hours 360] [--save-hours 1] '// &
            '[--da-order 4] [--max-depth 8] [--pos-tol-km 0.1] '// &
            '[--vel-tol-kms 1e-6] [-cfg config/config.txt]'
    end subroutine usage

    subroutine fatal(reason)
        character(len=*),intent(in) :: reason
        write(*,'(A)') 'ERROR: '//trim(reason)
        stop 1
    end subroutine fatal
end program run_srp_ads_history
