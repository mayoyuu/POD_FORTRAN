!> Paired ADS evaluation on a fixed Monte Carlo design for SRP uncertainty effects.
program run_srp_ads_mc_pair
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config
    use pod_data_format_module, only: load_initial_opm
    use pod_uq_ads_coordinates_module, only: ADS_COORD_WHITENED, ads_unit_to_physical
    use pod_uq_hfem_ads_module, only: hfem_ads_options_type, hfem_ads_history_type, &
        hfem_ads_split_event_type, hfem_ads_history_init, hfem_ads_history_advance, &
        hfem_ads_history_evaluate, hfem_ads_history_destroy
    use pod_srp_ads_history_module, only: SHAPE_COUNT, VALIDATION_COUNT, build_study_points, &
        assess_accuracy, advance_real_probes, validate_json_opm_fields
    use pod_srp_ads_snapshot_module, only: write_ads_patch_snapshot
    use pod_ads_split_module, only: sh_count
    use pod_srp_mc_io_module, only: read_srp_mc_samples, write_mc_point_header, &
        write_mc_point_rows
    use pod_srp_mc_statistics_module, only: srp_mc_group_stats_type, &
        srp_mc_pair_stats_type, compute_group_stats, compute_pair_stats
    implicit none

    integer, parameter :: PREFIX_COUNT=3
    integer, parameter :: PREFIX_LABELS(PREFIX_COUNT)=[10000,20000,50000]
    type(hfem_ads_options_type) :: treatment_options,control_options
    type(hfem_ads_history_type) :: treatment_history,control_history
    type(srp_mc_group_stats_type) :: treatment_stats,control_stats
    type(srp_mc_pair_stats_type) :: pair_stats
    type(hfem_ads_split_event_type) :: treatment_first_split,control_first_split
    type(hfem_ads_split_event_type) :: treatment_first_eta_split
    character(len=MAX_STRING_LEN) :: opm_file,prefix,config_file,sample_file,arg,message,path
    real(DP) :: epoch0,nominal0(6),cov(6,6),eta_sigma
    real(DP) :: pos_tolerance,vel_tolerance,duration_hours,save_hours
    real(DP) :: t_seconds,prev_seconds,final_seconds,current_hour,previous_hour,interval_count_real
    real(DP) :: treatment_center(6),control_center(6),mapped(7)
    real(DP) :: treatment_validation_pos,treatment_validation_vel
    real(DP) :: control_validation_pos,control_validation_vel
    real(DP) :: center_position_difference,center_velocity_difference
    real(DP) :: max_center_position_difference,max_center_velocity_difference
    real(DP),allocatable :: z(:,:),u(:,:),treatment_states(:,:),control_states(:,:)
    real(DP),allocatable :: treatment_errors(:,:),control_errors(:,:)
    real(DP),allocatable :: previous_treatment_errors(:,:),previous_control_errors(:,:)
    real(DP),allocatable :: probe_points(:,:),treatment_probe_ads(:,:),control_probe_ads(:,:)
    real(DP),allocatable :: treatment_shape(:,:),control_shape(:,:)
    real(DP),allocatable :: treatment_truth(:,:),control_truth(:,:)
    real(DP) :: treatment_shape_mean(6),control_shape_mean(6)
    real(DP) :: treatment_shape_pos_max,treatment_shape_vel_max
    real(DP) :: control_shape_pos_max,control_shape_vel_max
    real(DP),allocatable :: treatment_eta(:),control_eta(:)
    real(DP) :: treatment_zero(7,1),control_zero(6,1)
    logical,allocatable :: treatment_found(:),control_found(:),probe_found(:),shape_found(:)
    logical,allocatable :: keyframe_saved(:)
    integer,allocatable :: validation_ids(:),all_ids(:),visualization_ids(:)
    integer :: argc,i,j,k,n,actual_count,visualization_count
    integer :: status,step_index,steps,unit_hourly,unit_paired,unit_visual,unit_key
    integer :: unit_treatment_events,unit_control_events,treatment_worst,control_worst
    integer :: treatment_first_failure_step,control_first_failure_step
    integer :: treatment_previous_pass_step,control_previous_pass_step
    integer :: treatment_last_pass_step,control_last_pass_step
    integer :: treatment_first_unconfirmed_step,control_first_unconfirmed_step
    integer :: treatment_first_split_step,control_first_split_step,treatment_first_eta_step
    integer :: treatment_event_count,control_event_count
    integer :: local_event
    logical :: treatment_failed,control_failed,treatment_split_seen,control_split_seen
    logical :: treatment_eta_split_seen,keyframe,diagnostic_keyframe,pairing_complete
    logical :: has_opm,has_prefix,has_samples
    character(len=32) :: keyframe_kind

    opm_file=''
    prefix=''
    sample_file=''
    config_file='config/config.txt'
    eta_sigma=0.3_DP
    duration_hours=360.0_DP
    save_hours=1.0_DP
    pos_tolerance=0.1_DP
    vel_tolerance=1.0e-6_DP
    n=50000
    visualization_count=2048
    has_opm=.false.
    has_prefix=.false.
    has_samples=.false.

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
        case('--mc-samples')
            call get_command_argument(i+1,sample_file)
            has_samples=.true.
        case('--mc-count')
            call get_command_argument(i+1,arg)
            read(arg,*,iostat=status) n
            if(status/=0) call fatal('invalid --mc-count')
        case('--visualization-count')
            call get_command_argument(i+1,arg)
            read(arg,*,iostat=status) visualization_count
            if(status/=0) call fatal('invalid --visualization-count')
        case('--eta-sigma')
            call get_command_argument(i+1,arg)
            read(arg,*,iostat=status) eta_sigma
            if(status/=0) call fatal('invalid --eta-sigma')
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
            read(arg,*,iostat=status) treatment_options%da_order
            if(status/=0) call fatal('invalid --da-order')
            control_options%da_order=treatment_options%da_order
        case('--max-depth')
            call get_command_argument(i+1,arg)
            read(arg,*,iostat=status) treatment_options%max_split_depth
            if(status/=0) call fatal('invalid --max-depth')
            control_options%max_split_depth=treatment_options%max_split_depth
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
    if(.not.has_opm.or..not.has_prefix.or..not.has_samples) then
        call usage()
        call fatal('-opm, -o, and --mc-samples are required')
    end if
    if(n<2.or.visualization_count<1.or.duration_hours<0.0_DP.or.save_hours<=0.0_DP) &
        call fatal('MC count must be at least two; duration/cadence/count are invalid')
    if(.not.ieee_is_finite(eta_sigma).or.eta_sigma<=0.0_DP) &
        call fatal('treatment --eta-sigma must be finite and positive')
    if(pos_tolerance<=0.0_DP.or.vel_tolerance<=0.0_DP) call fatal('tolerances must be positive')

    call validate_json_opm_fields(trim(opm_file),status,message)
    if(status/=0) call fatal('invalid JSON OPM: '//trim(message))
    call read_srp_mc_samples(trim(sample_file),n,z,u,status,message,eta_sigma)
    if(status/=0) call fatal('MC sample input: '//trim(message))
    call pod_engine_init(trim(config_file))
    config%use_srp=.true.
    call load_initial_opm(trim(opm_file),epoch0,nominal0,cov)

    treatment_options%coordinate_mode=ADS_COORD_WHITENED
    treatment_options%srp_sigma=eta_sigma
    treatment_options%da_storage_variables=7
    treatment_options%error_tolerance=[pos_tolerance,pos_tolerance,pos_tolerance, &
        vel_tolerance,vel_tolerance,vel_tolerance]
    control_options=treatment_options
    control_options%srp_sigma=0.0_DP
    control_options%da_storage_variables=7

    call hfem_ads_history_init(treatment_history,nominal0,cov,epoch0, &
        treatment_options,status,message)
    if(status/=0) call fatal('treatment ADS init: '//trim(message))
    call hfem_ads_history_init(control_history,nominal0,cov,epoch0, &
        control_options,status,message)
    if(status/=0) call fatal('control ADS init: '//trim(message))
    if(treatment_history%coordinate_map%n_variables/=7.or. &
       control_history%coordinate_map%n_variables/=6.or. &
       treatment_history%da_variables/=7.or.control_history%da_variables/=7) &
        call fatal('paired ADS dimension contract was not established')

    allocate(treatment_states(6,n),control_states(6,n),treatment_errors(6,n), &
        control_errors(6,n),previous_treatment_errors(6,n), &
        previous_control_errors(6,n),treatment_found(n),control_found(n),all_ids(n))
    all_ids=[(i,i=1,n)]
    call choose_visualization_ids(n,min(visualization_count,n),visualization_ids)

    call build_study_points(probe_points,validation_ids)
    allocate(treatment_probe_ads(6,VALIDATION_COUNT),control_probe_ads(6,VALIDATION_COUNT), &
        treatment_shape(6,SHAPE_COUNT),control_shape(6,SHAPE_COUNT),shape_found(SHAPE_COUNT), &
        treatment_truth(6,VALIDATION_COUNT),control_truth(6,VALIDATION_COUNT), &
        treatment_eta(VALIDATION_COUNT),control_eta(VALIDATION_COUNT), &
        probe_found(VALIDATION_COUNT))
    do j=1,VALIDATION_COUNT
        call ads_unit_to_physical(treatment_history%coordinate_map, &
            probe_points(:,validation_ids(j)),mapped(1:6),mapped(7),status)
        if(status/=0) call fatal('treatment validation mapping failed')
        treatment_truth(:,j)=nominal0+mapped(1:6)
        treatment_eta(j)=mapped(7)
        call ads_unit_to_physical(control_history%coordinate_map, &
            probe_points(1:6,validation_ids(j)),mapped(1:6),mapped(7),status)
        if(status/=0) call fatal('control validation mapping failed')
        control_truth(:,j)=nominal0+mapped(1:6)
        control_eta(j)=0.0_DP
    end do

    call open_outputs()
    treatment_first_failure_step=-1
    control_first_failure_step=-1
    treatment_previous_pass_step=-1
    control_previous_pass_step=-1
    treatment_last_pass_step=-1
    control_last_pass_step=-1
    treatment_first_unconfirmed_step=-1
    control_first_unconfirmed_step=-1
    treatment_first_split_step=-1
    control_first_split_step=-1
    treatment_first_eta_step=-1
    treatment_event_count=0
    control_event_count=0
    treatment_split_seen=.false.
    control_split_seen=.false.
    treatment_eta_split_seen=.false.
    max_center_position_difference=0.0_DP
    max_center_velocity_difference=0.0_DP

    final_seconds=duration_hours*3600.0_DP
    interval_count_real=duration_hours/save_hours
    if(abs(interval_count_real-real(nint(interval_count_real),DP))<= &
       1.0e-12_DP*max(1.0_DP,abs(interval_count_real))) then
        steps=nint(interval_count_real)
    else
        steps=ceiling(interval_count_real)
    end if
    allocate(keyframe_saved(0:steps),source=.false.)
    previous_treatment_errors=0.0_DP
    previous_control_errors=0.0_DP
    prev_seconds=0.0_DP

    do step_index=0,steps
        t_seconds=min(final_seconds,real(step_index,DP)*save_hours*3600.0_DP)
        current_hour=t_seconds/3600.0_DP
        previous_hour=prev_seconds/3600.0_DP

        if(step_index>0) then
            call hfem_ads_history_advance(treatment_history,t_seconds,status,message)
            if(status/=0) call fatal('treatment ADS advance: '//trim(message))
            call capture_events('treatment',treatment_history,unit_treatment_events, &
                treatment_event_count,treatment_split_seen,treatment_first_split_step, &
                treatment_first_split,treatment_eta_split_seen,treatment_first_eta_step, &
                treatment_first_eta_split)
            call hfem_ads_history_advance(control_history,t_seconds,status,message)
            if(status/=0) call fatal('control ADS advance: '//trim(message))
            call capture_events('control',control_history,unit_control_events, &
                control_event_count,control_split_seen,control_first_split_step, &
                control_first_split)

            call advance_real_probes(treatment_truth,treatment_eta,epoch0,prev_seconds, &
                t_seconds,treatment_options,status,message)
            if(status/=0) call fatal('treatment real validation advance: '//trim(message))
            call advance_real_probes(control_truth,control_eta,epoch0,prev_seconds, &
                t_seconds,control_options,status,message)
            if(status/=0) call fatal('control real validation advance: '//trim(message))
        end if

        call hfem_ads_history_evaluate(treatment_history,u,treatment_states, &
            treatment_found,status)
        if(status/=0.or..not.all(treatment_found)) call fatal('treatment MC evaluation failed')
        call hfem_ads_history_evaluate(control_history,u(1:6,:),control_states, &
            control_found,status)
        if(status/=0.or..not.all(control_found)) call fatal('control MC evaluation failed')

        treatment_zero=0.0_DP
        call hfem_ads_history_evaluate(treatment_history,treatment_zero, &
            treatment_probe_ads(:,1:1),probe_found(1:1),status)
        if(status/=0.or..not.probe_found(1)) call fatal('treatment center evaluation failed')
        treatment_center=treatment_probe_ads(:,1)
        control_zero=0.0_DP
        call hfem_ads_history_evaluate(control_history,control_zero, &
            control_probe_ads(:,1:1),probe_found(1:1),status)
        if(status/=0.or..not.probe_found(1)) call fatal('control center evaluation failed')
        control_center=control_probe_ads(:,1)

        call hfem_ads_history_evaluate(treatment_history,probe_points,treatment_shape, &
            shape_found,status)
        if(status/=0.or..not.all(shape_found)) call fatal('treatment shape evaluation failed')
        treatment_probe_ads=treatment_shape(:,validation_ids)
        call assess_accuracy(treatment_probe_ads,treatment_truth,pos_tolerance, &
            vel_tolerance,treatment_validation_pos,treatment_validation_vel, &
            treatment_worst,treatment_failed)
        call hfem_ads_history_evaluate(control_history,probe_points(1:6,:),control_shape, &
            shape_found,status)
        if(status/=0.or..not.all(shape_found)) call fatal('control shape evaluation failed')
        control_probe_ads=control_shape(:,validation_ids)
        call assess_accuracy(control_probe_ads,control_truth,pos_tolerance, &
            vel_tolerance,control_validation_pos,control_validation_vel, &
            control_worst,control_failed)
        call compute_shape_summary(treatment_shape,treatment_center,treatment_shape_mean, &
            treatment_shape_pos_max,treatment_shape_vel_max)
        call compute_shape_summary(control_shape,control_center,control_shape_mean, &
            control_shape_pos_max,control_shape_vel_max)

        call update_failure(treatment_failed,step_index,treatment_last_pass_step, &
            treatment_first_failure_step,treatment_previous_pass_step)
        call update_failure(control_failed,step_index,control_last_pass_step, &
            control_first_failure_step,control_previous_pass_step)
        if(treatment_history%stats%depth_limited_patches>0.and. &
           treatment_first_unconfirmed_step<0) treatment_first_unconfirmed_step=step_index
        if(control_history%stats%depth_limited_patches>0.and. &
           control_first_unconfirmed_step<0) control_first_unconfirmed_step=step_index

        treatment_errors=treatment_states-spread(treatment_center,2,n)
        control_errors=control_states-spread(control_center,2,n)
        center_position_difference=sqrt(sum((treatment_center(1:3)-control_center(1:3))**2))
        center_velocity_difference=sqrt(sum((treatment_center(4:6)-control_center(4:6))**2))
        max_center_position_difference=max(max_center_position_difference,center_position_difference)
        max_center_velocity_difference=max(max_center_velocity_difference,center_velocity_difference)

        do k=1,PREFIX_COUNT
            actual_count=min(PREFIX_LABELS(k),n)
            call compute_group_stats(treatment_states(:,1:actual_count), &
                treatment_center,treatment_stats,status)
            if(status/=0) call fatal('treatment statistics failed')
            call write_group_row('treatment',PREFIX_LABELS(k),actual_count, &
                treatment_stats,treatment_history,treatment_validation_pos, &
                treatment_validation_vel,treatment_worst,treatment_failed)
            call compute_group_stats(control_states(:,1:actual_count), &
                control_center,control_stats,status)
            if(status/=0) call fatal('control statistics failed')
            call write_group_row('control',PREFIX_LABELS(k),actual_count, &
                control_stats,control_history,control_validation_pos, &
                control_validation_vel,control_worst,control_failed)
            call compute_pair_stats(treatment_states(:,1:actual_count), &
                control_states(:,1:actual_count),treatment_stats,control_stats, &
                pair_stats,status)
            if(status/=0) call fatal('paired statistics failed')
            call write_pair_row(PREFIX_LABELS(k),actual_count,pair_stats, &
                center_position_difference,center_velocity_difference)
        end do

        call write_mc_point_rows(unit_visual,current_hour,'hourly_visualization', &
            visualization_ids,treatment_errors(:,visualization_ids), &
            control_errors(:,visualization_ids),status,message)
        if(status/=0) call fatal('visualization point output: '//trim(message))

        diagnostic_keyframe=treatment_history%stats%split_occurred.or. &
            control_history%stats%split_occurred.or. &
            step_index==treatment_first_failure_step.or. &
            step_index==control_first_failure_step.or. &
            step_index==treatment_first_unconfirmed_step.or. &
            step_index==control_first_unconfirmed_step
        keyframe=is_keyframe(current_hour,step_index,steps).or.diagnostic_keyframe
        keyframe_kind=''
        if(step_index==0) then
            keyframe_kind='fixed_start'
        else if(diagnostic_keyframe) then
            keyframe_kind='diagnostic_event'
        else if(step_index==steps) then
            keyframe_kind='fixed_end'
        else
            keyframe_kind='fixed_hour'
        end if
        if(diagnostic_keyframe.and.step_index>0) then
            if(.not.keyframe_saved(step_index-1)) then
                call write_mc_point_rows(unit_key,previous_hour,'diagnostic_previous',all_ids, &
                    previous_treatment_errors,previous_control_errors,status,message)
                if(status/=0) call fatal('previous keyframe point output: '//trim(message))
                keyframe_saved(step_index-1)=.true.
            end if
        end if
        if(keyframe.and..not.keyframe_saved(step_index)) then
            call write_mc_point_rows(unit_key,current_hour,trim(keyframe_kind),all_ids, &
                treatment_errors,control_errors,status,message)
            if(status/=0) call fatal('keyframe point output: '//trim(message))
            keyframe_saved(step_index)=.true.
            call write_snapshots(.false.,nint(t_seconds))
            if(diagnostic_keyframe.and.step_index>0) then
                call write_snapshots(.true.,nint(prev_seconds))
            end if
        end if
        flush(unit_hourly)
        flush(unit_paired)
        flush(unit_visual)
        flush(unit_key)
        write(*,'(A,F8.2,A,I0,A,I0,A,ES10.3)') 'hour ',current_hour, &
            ' treatment patches ',treatment_history%domain%n_patches, &
            ' control patches ',control_history%domain%n_patches, &
            ' center delta km ',center_position_difference
        previous_treatment_errors=treatment_errors
        previous_control_errors=control_errors
        prev_seconds=t_seconds
    end do

    close(unit_hourly)
    close(unit_paired)
    close(unit_visual)
    close(unit_key)
    close(unit_treatment_events)
    close(unit_control_events)
    pairing_complete=max_center_position_difference<=1.0e-9_DP.and. &
        max_center_velocity_difference<=1.0e-12_DP
    call write_report()
    call hfem_ads_history_destroy(treatment_history)
    call hfem_ads_history_destroy(control_history)

contains

    subroutine choose_visualization_ids(total,count,ids)
        integer,intent(in) :: total,count
        integer,allocatable,intent(out) :: ids(:)
        integer :: q
        allocate(ids(count))
        if(count==1) then
            ids(1)=1
        else
            do q=1,count
                ids(q)=1+nint(real(q-1,DP)*real(total-1,DP)/real(count-1,DP))
            end do
        end if
    end subroutine choose_visualization_ids

    subroutine update_failure(failed,step,last_pass,first_failure,previous_pass)
        logical,intent(in) :: failed
        integer,intent(in) :: step
        integer,intent(inout) :: last_pass,first_failure,previous_pass
        if(failed.and.first_failure<0) then
            first_failure=step
            previous_pass=last_pass
        end if
        if(.not.failed) last_pass=step
    end subroutine update_failure

    logical function is_keyframe(hour,step,last_step)
        real(DP),intent(in) :: hour
        integer,intent(in) :: step,last_step
        real(DP),parameter :: fixed(5)=[24.0_DP,72.0_DP,168.0_DP,240.0_DP,360.0_DP]
        is_keyframe=step==0.or.step==last_step.or.any(abs(hour-fixed)<1.0e-9_DP)
    end function is_keyframe

    subroutine open_outputs()
        character(len=*),parameter :: EVENT_HEADER='interval_start_hour,interval_end_hour,'// &
            'event_id,step_event_index,bfs_iteration,parent_depth,component,direction,'// &
            'error1,error2,error3,error4,error5,error6,'// &
            'tolerance1,tolerance2,tolerance3,tolerance4,tolerance5,tolerance6,'// &
            'excess1,excess2,excess3,excess4,excess5,excess6,'// &
            'direction_score1,direction_score2,direction_score3,direction_score4,'// &
            'direction_score5,direction_score6,direction_score7,'// &
            'parent_center1,parent_center2,parent_center3,parent_center4,'// &
            'parent_center5,parent_center6,parent_center7,'// &
            'parent_width1,parent_width2,parent_width3,parent_width4,'// &
            'parent_width5,parent_width6,parent_width7'
        open(newunit=unit_hourly,file=trim(prefix)//'_mc_hourly_stats.csv', &
            status='replace',action='write',iostat=status)
        if(status/=0) call fatal('cannot open hourly statistics output')
        open(newunit=unit_paired,file=trim(prefix)//'_mc_paired_stats.csv', &
            status='replace',action='write',iostat=status)
        if(status/=0) call fatal('cannot open paired statistics output')
        open(newunit=unit_visual,file=trim(prefix)//'_mc_visualization_points.csv', &
            status='replace',action='write',iostat=status)
        if(status/=0) call fatal('cannot open visualization output')
        open(newunit=unit_key,file=trim(prefix)//'_mc_keyframe_points.csv', &
            status='replace',action='write',iostat=status)
        if(status/=0) call fatal('cannot open keyframe output')
        open(newunit=unit_treatment_events,file=trim(prefix)// &
            '_treatment_split_events.csv',status='replace',action='write',iostat=status)
        if(status/=0) call fatal('cannot open treatment split events')
        open(newunit=unit_control_events,file=trim(prefix)// &
            '_control_split_events.csv',status='replace',action='write',iostat=status)
        if(status/=0) call fatal('cannot open control split events')

        write(unit_hourly,'(A)') 'hour,epoch_et,group,sample_count,actual_sample_count,'// &
            'active_dimensions,patch_count,depth_limited_patches,'// &
            'split_dim1,split_dim2,split_dim3,split_dim4,split_dim5,split_dim6,split_dim7,'// &
            'bfs_iterations,max_queue_size,step_elapsed_seconds,'// &
            'validation_position_error_km,validation_velocity_error_kms,'// &
            'validation_worst_point,validation_failed,ads_accuracy_valid,'// &
            'precision_confirmed,maximum_patch_depth,'// &
            'center_x_km,center_y_km,center_z_km,center_vx_kms,center_vy_kms,center_vz_kms,'// &
            'mean_dx_km,mean_dy_km,mean_dz_km,mean_dvx_kms,mean_dvy_kms,mean_dvz_kms,'// &
            'std_dx_km,std_dy_km,std_dz_km,std_dvx_kms,std_dvy_kms,std_dvz_kms,'// &
            'min_dx_km,min_dy_km,min_dz_km,min_dvx_kms,min_dvy_kms,min_dvz_kms,'// &
            'max_dx_km,max_dy_km,max_dz_km,max_dvx_kms,max_dvy_kms,max_dvz_kms,'// &
            'skew_dx,skew_dy,skew_dz,skew_dvx,skew_dvy,skew_dvz,'// &
            'excess_kurt_dx,excess_kurt_dy,excess_kurt_dz,'// &
            'excess_kurt_dvx,excess_kurt_dvy,excess_kurt_dvz,'// &
            'position_axis1_rms_km,position_axis2_rms_km,position_axis3_rms_km,'// &
            'velocity_axis1_rms_kms,velocity_axis2_rms_kms,velocity_axis3_rms_kms,'// &
            'position_norm_rms_km,velocity_norm_rms_kms,'// &
            'position_p50_km,position_p95_km,position_p99_km,position_p999_km,'// &
            'velocity_p50_kms,velocity_p95_kms,velocity_p99_kms,velocity_p999_kms,'// &
            'position_max_km,velocity_max_kms,position_worst_id,velocity_worst_id,'// &
            'position_max_to_axis1,velocity_max_to_axis1,'// &
            'shape_mean_dx_km,shape_mean_dy_km,shape_mean_dz_km,'// &
            'shape_mean_dvx_kms,shape_mean_dvy_kms,shape_mean_dvz_kms,'// &
            'shape_position_max_km,shape_velocity_max_kms,'// &
            'shape_position_max_to_axis1,shape_velocity_max_to_axis1,'// &
            'max_abs_skewness,max_abs_skewness_component,'// &
            'max_abs_excess_kurtosis,max_abs_excess_kurtosis_component,'// &
            'cov_position_xx,cov_position_xy,cov_position_xz,cov_position_yy,'// &
            'cov_position_yz,cov_position_zz,cov_velocity_xx,cov_velocity_xy,'// &
            'cov_velocity_xz,cov_velocity_yy,cov_velocity_yz,cov_velocity_zz'
        write(unit_paired,'(A)') 'hour,epoch_et,sample_count,actual_sample_count,'// &
            'center_position_difference_km,center_velocity_difference_kms,'// &
            'mean_delta_x_km,mean_delta_y_km,mean_delta_z_km,'// &
            'mean_delta_vx_kms,mean_delta_vy_kms,mean_delta_vz_kms,'// &
            'std_delta_x_km,std_delta_y_km,std_delta_z_km,'// &
            'std_delta_vx_kms,std_delta_vy_kms,std_delta_vz_kms,'// &
            'min_delta_x_km,min_delta_y_km,min_delta_z_km,'// &
            'min_delta_vx_kms,min_delta_vy_kms,min_delta_vz_kms,'// &
            'max_delta_x_km,max_delta_y_km,max_delta_z_km,'// &
            'max_delta_vx_kms,max_delta_vy_kms,max_delta_vz_kms,'// &
            'delta_position_rms_km,delta_velocity_rms_kms,'// &
            'delta_position_p50_km,delta_position_p95_km,delta_position_p99_km,'// &
            'delta_position_p999_km,delta_velocity_p50_kms,delta_velocity_p95_kms,'// &
            'delta_velocity_p99_kms,delta_velocity_p999_kms,'// &
            'delta_position_max_km,delta_velocity_max_kms,'// &
            'delta_position_worst_id,delta_velocity_worst_id,'// &
            'position_max_increment_km,velocity_max_increment_kms,'// &
            'position_max_amplification,velocity_max_amplification,'// &
            'axis1_position_increment_km,axis1_velocity_increment_kms,'// &
            'axis1_position_amplification,axis1_velocity_amplification,'// &
            'delta_position_max_to_treatment_axis1,delta_position_max_to_control_axis1,'// &
            'delta_velocity_max_to_treatment_axis1,delta_velocity_max_to_control_axis1,'// &
            'position_p95_increment_km,position_p99_increment_km,'// &
            'velocity_p95_increment_kms,velocity_p99_increment_kms,'// &
            'position_p95_amplification,position_p99_amplification,'// &
            'velocity_p95_amplification,velocity_p99_amplification,'// &
            'cov_delta_11,cov_delta_12,cov_delta_13,cov_delta_14,cov_delta_15,cov_delta_16,'// &
            'cov_delta_21,cov_delta_22,cov_delta_23,cov_delta_24,cov_delta_25,cov_delta_26,'// &
            'cov_delta_31,cov_delta_32,cov_delta_33,cov_delta_34,cov_delta_35,cov_delta_36,'// &
            'cov_delta_41,cov_delta_42,cov_delta_43,cov_delta_44,cov_delta_45,cov_delta_46,'// &
            'cov_delta_51,cov_delta_52,cov_delta_53,cov_delta_54,cov_delta_55,cov_delta_56,'// &
            'cov_delta_61,cov_delta_62,cov_delta_63,cov_delta_64,cov_delta_65,cov_delta_66'
        call write_mc_point_header(unit_visual,status,message)
        if(status/=0) call fatal(trim(message))
        call write_mc_point_header(unit_key,status,message)
        if(status/=0) call fatal(trim(message))
        write(unit_treatment_events,'(A)') EVENT_HEADER
        write(unit_control_events,'(A)') EVENT_HEADER
    end subroutine open_outputs

    subroutine write_group_row(group,label,count_value,stats,history,validation_pos, &
                               validation_vel,worst,failed)
        character(len=*),intent(in) :: group
        integer,intent(in) :: label,count_value,worst
        type(srp_mc_group_stats_type),intent(in) :: stats
        type(hfem_ads_history_type),intent(in) :: history
        real(DP),intent(in) :: validation_pos,validation_vel
        logical,intent(in) :: failed
        integer :: split7(7),maximum_depth,accuracy_valid,precision_confirmed
        integer :: skew_component,kurtosis_component
        real(DP) :: shape_mean(6),shape_pos_max,shape_vel_max
        logical :: is_treatment
        is_treatment=trim(group)=='treatment'
        shape_mean=merge(treatment_shape_mean,control_shape_mean,is_treatment)
        shape_pos_max=merge(treatment_shape_pos_max,control_shape_pos_max,is_treatment)
        shape_vel_max=merge(treatment_shape_vel_max,control_shape_vel_max,is_treatment)
        if(is_treatment) then
            accuracy_valid=merge(1,0,treatment_first_failure_step<0.or. &
                step_index<treatment_first_failure_step)
            precision_confirmed=merge(1,0,treatment_first_unconfirmed_step<0.or. &
                step_index<treatment_first_unconfirmed_step)
        else
            accuracy_valid=merge(1,0,control_first_failure_step<0.or. &
                step_index<control_first_failure_step)
            precision_confirmed=merge(1,0,control_first_unconfirmed_step<0.or. &
                step_index<control_first_unconfirmed_step)
        end if
        maximum_depth=max_patch_depth(history)
        skew_component=maxloc(abs(stats%skewness),dim=1)
        kurtosis_component=maxloc(abs(stats%excess_kurtosis),dim=1)
        split7=0
        split7(1:size(history%stats%split_counts))=history%stats%split_counts
        write(unit_hourly,'(*(G0,:,","))') current_hour,epoch0+t_seconds,trim(group), &
            label,count_value,history%coordinate_map%n_variables,history%domain%n_patches, &
            history%stats%depth_limited_patches,split7,history%stats%bfs_iterations, &
            history%stats%max_queue_size,history%stats%elapsed_seconds,validation_pos, &
            validation_vel,worst,merge(1,0,failed),accuracy_valid,precision_confirmed, &
            maximum_depth,merge(treatment_center,control_center,is_treatment), &
            stats%mean_error,stats%standard_deviation,stats%component_min,stats%component_max, &
            stats%skewness,stats%excess_kurtosis,stats%position_axes_rms, &
            stats%velocity_axes_rms,stats%position_norm_rms,stats%velocity_norm_rms, &
            stats%position_norm_quantiles,stats%velocity_norm_quantiles, &
            stats%position_norm_max,stats%velocity_norm_max,stats%position_worst_id, &
            stats%velocity_worst_id,stats%position_max_to_longest_axis, &
            stats%velocity_max_to_longest_axis,shape_mean,shape_pos_max,shape_vel_max, &
            safe_ratio_local(shape_pos_max,stats%position_axes_rms(1)), &
            safe_ratio_local(shape_vel_max,stats%velocity_axes_rms(1)), &
            maxval(abs(stats%skewness)),skew_component, &
            maxval(abs(stats%excess_kurtosis)),kurtosis_component, &
            stats%covariance(1,1),stats%covariance(1,2),stats%covariance(1,3), &
            stats%covariance(2,2),stats%covariance(2,3),stats%covariance(3,3), &
            stats%covariance(4,4),stats%covariance(4,5),stats%covariance(4,6), &
            stats%covariance(5,5),stats%covariance(5,6),stats%covariance(6,6)
    end subroutine write_group_row

    subroutine write_pair_row(label,count_value,stats,center_pos,center_vel)
        integer,intent(in) :: label,count_value
        type(srp_mc_pair_stats_type),intent(in) :: stats
        real(DP),intent(in) :: center_pos,center_vel
        write(unit_paired,'(*(G0,:,","))') current_hour,epoch0+t_seconds,label,count_value, &
            center_pos,center_vel,stats%mean_delta,stats%standard_deviation, &
            stats%component_delta_min,stats%component_delta_max, &
            stats%position_delta_rms,stats%velocity_delta_rms, &
            stats%position_delta_quantiles,stats%velocity_delta_quantiles, &
            stats%position_delta_max,stats%velocity_delta_max, &
            stats%position_worst_id,stats%velocity_worst_id, &
            stats%position_max_increment,stats%velocity_max_increment, &
            stats%position_max_amplification,stats%velocity_max_amplification, &
            stats%longest_position_axis_increment,stats%longest_velocity_axis_increment, &
            stats%longest_position_axis_amplification, &
            stats%longest_velocity_axis_amplification, &
            stats%delta_position_max_to_treatment_axis, &
            stats%delta_position_max_to_control_axis, &
            stats%delta_velocity_max_to_treatment_axis, &
            stats%delta_velocity_max_to_control_axis, &
            treatment_stats%position_norm_quantiles(2)-control_stats%position_norm_quantiles(2), &
            treatment_stats%position_norm_quantiles(3)-control_stats%position_norm_quantiles(3), &
            treatment_stats%velocity_norm_quantiles(2)-control_stats%velocity_norm_quantiles(2), &
            treatment_stats%velocity_norm_quantiles(3)-control_stats%velocity_norm_quantiles(3), &
            safe_ratio_local(treatment_stats%position_norm_quantiles(2), &
                control_stats%position_norm_quantiles(2)), &
            safe_ratio_local(treatment_stats%position_norm_quantiles(3), &
                control_stats%position_norm_quantiles(3)), &
            safe_ratio_local(treatment_stats%velocity_norm_quantiles(2), &
                control_stats%velocity_norm_quantiles(2)), &
            safe_ratio_local(treatment_stats%velocity_norm_quantiles(3), &
                control_stats%velocity_norm_quantiles(3)), &
            ((stats%covariance(i,j),j=1,6),i=1,6)
    end subroutine write_pair_row

    subroutine compute_shape_summary(states,center,mean_error,pos_max,vel_max)
        real(DP),intent(in) :: states(:,:),center(6)
        real(DP),intent(out) :: mean_error(6),pos_max,vel_max
        real(DP),allocatable :: errors(:,:)
        allocate(errors(6,size(states,2)))
        errors=states-spread(center,dim=2,ncopies=size(states,2))
        mean_error=sum(errors,dim=2)/real(size(states,2),DP)
        pos_max=maxval(sqrt(sum(errors(1:3,:)**2,dim=1)))
        vel_max=maxval(sqrt(sum(errors(4:6,:)**2,dim=1)))
        deallocate(errors)
    end subroutine compute_shape_summary

    integer function max_patch_depth(history) result(depth)
        type(hfem_ads_history_type),intent(in) :: history
        integer :: q
        depth=0
        do q=1,history%domain%n_patches
            depth=max(depth,sh_count(history%domain%patches(q)%history,0))
        end do
    end function max_patch_depth

    real(DP) function safe_ratio_local(numerator,denominator) result(value)
        real(DP),intent(in) :: numerator,denominator
        if(abs(denominator)>100.0_DP*epsilon(1.0_DP)) then
            value=numerator/denominator
        else if(abs(numerator)<=100.0_DP*epsilon(1.0_DP)) then
            value=0.0_DP
        else
            value=sign(huge(1.0_DP),numerator)
        end if
    end function safe_ratio_local

    subroutine capture_events(group,history,output_unit,event_count,split_seen, &
                              first_step,first_event,eta_seen,eta_step,eta_event)
        character(len=*),intent(in) :: group
        type(hfem_ads_history_type),intent(in) :: history
        integer,intent(in) :: output_unit
        integer,intent(inout) :: event_count,first_step
        logical,intent(inout) :: split_seen
        type(hfem_ads_split_event_type),intent(inout) :: first_event
        logical,optional,intent(inout) :: eta_seen
        integer,optional,intent(inout) :: eta_step
        type(hfem_ads_split_event_type),optional,intent(inout) :: eta_event
        if(.not.allocated(history%stats%split_events)) return
        do local_event=1,size(history%stats%split_events)
            event_count=event_count+1
            call write_split_event(output_unit,history%stats%split_events(local_event), &
                local_event,event_count)
            if(.not.split_seen) then
                split_seen=.true.
                first_step=step_index
                first_event=history%stats%split_events(local_event)
            end if
            if(present(eta_seen).and.present(eta_step).and.present(eta_event)) then
                if(history%stats%split_events(local_event)%direction==7.and..not.eta_seen) then
                    eta_seen=.true.
                    eta_step=step_index
                    eta_event=history%stats%split_events(local_event)
                end if
            end if
        end do
    end subroutine capture_events

    subroutine write_split_event(output_unit,event,local_index,global_index)
        integer,intent(in) :: output_unit,local_index,global_index
        type(hfem_ads_split_event_type),intent(in) :: event
        write(output_unit,'(*(G0,:,","))') previous_hour,current_hour,global_index, &
            local_index,event%bfs_iteration,event%parent_depth,event%component, &
            event%direction,event%errors,event%tolerances,event%excess, &
            event%direction_scores,event%parent_center,event%parent_width
    end subroutine write_split_event

    subroutine write_snapshots(use_previous,time_seconds)
        logical,intent(in) :: use_previous
        integer,intent(in) :: time_seconds
        character(len=16) :: tag
        write(tag,'(I7.7)') time_seconds
        path=trim(prefix)//'_treatment_patch_s'//trim(tag)
        call write_ads_patch_snapshot(treatment_history,trim(path),use_previous,status,message)
        if(status/=0) call fatal('treatment patch snapshot: '//trim(message))
        path=trim(prefix)//'_control_patch_s'//trim(tag)
        call write_ads_patch_snapshot(control_history,trim(path),use_previous,status,message)
        if(status/=0) call fatal('control patch snapshot: '//trim(message))
    end subroutine write_snapshots

    subroutine write_report()
        integer :: report_unit,ios
        open(newunit=report_unit,file=trim(prefix)//'_mc_report.json', &
            status='replace',action='write',iostat=ios)
        if(ios/=0) call fatal('cannot open MC JSON report')
        write(report_unit,'(A)') '{'
        write(report_unit,'(A)') '  "output_schema_version": "srp_ads_mc_pair_v1",'
        write(report_unit,'(A)') '  "method": "paired fixed-design ADS Monte Carlo evaluation",'
        write(report_unit,'(A)') '  "nominal_srp_enabled": true,'
        write(report_unit,'(A)') '  "control_srp_uncertainty_enabled": false,'
        write(report_unit,'(A,I0,A)') '  "treatment_active_dimensions": ', &
            treatment_history%coordinate_map%n_variables,','
        write(report_unit,'(A,I0,A)') '  "control_active_dimensions": ', &
            control_history%coordinate_map%n_variables,','
        write(report_unit,'(A,I0,A)') '  "dace_storage_dimensions": ', &
            treatment_history%da_variables,','
        write(report_unit,'(A,I0,A)') '  "sample_count": ',n,','
        write(report_unit,'(A,I0,A)') '  "visualization_sample_count": ', &
            size(visualization_ids),','
        write(report_unit,'(A,ES25.16E3,A)') '  "duration_hours": ',duration_hours,','
        write(report_unit,'(A,ES25.16E3,A)') '  "save_hours": ',save_hours,','
        write(report_unit,'(A,ES25.16E3,A)') '  "treatment_eta_sigma": ',eta_sigma,','
        write(report_unit,'(A,ES25.16E3,A)') '  "control_eta_sigma": ',0.0_DP,','
        write(report_unit,'(A,I0,A)') '  "da_order": ',treatment_options%da_order,','
        write(report_unit,'(A,I0,A)') '  "max_split_depth": ', &
            treatment_options%max_split_depth,','
        write(report_unit,'(A,ES25.16E3,A)') '  "position_tolerance_km": ', &
            pos_tolerance,','
        write(report_unit,'(A,ES25.16E3,A)') '  "velocity_tolerance_kms": ', &
            vel_tolerance,','
        write(report_unit,'(A)') '  "sample_statistics_interpretation": '// &
            '"fixed truncated-normal design metrics, not probability confidence intervals",'
        write(report_unit,'(A)') '  "split_decision_rule": '// &
            '"split when any component truncation estimate exceeds its tolerance; '// &
            'component maximizes positive excess and direction maximizes top-order score",'
        write(report_unit,'(A)') '  "opm_file": "'//trim(opm_file)//'",'
        write(report_unit,'(A)') '  "sample_file": "'//trim(sample_file)//'",'
        write(report_unit,'(A,ES25.16E3,A)') '  "maximum_center_position_difference_km": ', &
            max_center_position_difference,','
        write(report_unit,'(A,ES25.16E3,A)') '  "maximum_center_velocity_difference_kms": ', &
            max_center_velocity_difference,','
        if(pairing_complete) then
            write(report_unit,'(A)') '  "pairing_complete": true,'
        else
            write(report_unit,'(A)') '  "pairing_complete": false,'
        end if
        call write_group_events(report_unit,'treatment_events',treatment_split_seen, &
            treatment_first_split_step,treatment_first_split, &
            treatment_first_unconfirmed_step,treatment_first_failure_step, &
            treatment_previous_pass_step,.true.)
        write(report_unit,'(A)') ','
        call write_group_events(report_unit,'control_events',control_split_seen, &
            control_first_split_step,control_first_split,control_first_unconfirmed_step, &
            control_first_failure_step,control_previous_pass_step,.false.)
        write(report_unit,'(A)') ''
        write(report_unit,'(A)') '}'
        close(report_unit)
    end subroutine write_report

    subroutine write_group_events(output_unit,key,split_seen,split_step,split_event, &
                                  unconfirmed_step,failure_step,previous_pass,include_eta)
        integer,intent(in) :: output_unit,split_step,unconfirmed_step,failure_step,previous_pass
        character(len=*),intent(in) :: key
        logical,intent(in) :: split_seen,include_eta
        type(hfem_ads_split_event_type),intent(in) :: split_event
        write(output_unit,'(A)') '  "'//trim(key)//'": {'
        call write_named_split(output_unit,'first_split',split_seen,split_step,split_event,.true.)
        if(include_eta) then
            call write_named_split(output_unit,'first_eta_split',treatment_eta_split_seen, &
                treatment_first_eta_step,treatment_first_eta_split,.true.)
        end if
        call write_step_event(output_unit,'first_unconfirmed',unconfirmed_step,.true.)
        call write_failure_event(output_unit,failure_step,previous_pass)
        write(output_unit,'(A)',advance='no') '  }'
    end subroutine write_group_events

    subroutine write_named_split(output_unit,key,seen,step,event,trailing_comma)
        integer,intent(in) :: output_unit,step
        character(len=*),intent(in) :: key
        logical,intent(in) :: seen,trailing_comma
        type(hfem_ads_split_event_type),intent(in) :: event
        character(len=1) :: suffix
        suffix=merge(',', ' ',trailing_comma)
        if(.not.seen) then
            write(output_unit,'(A)') '    "'//trim(key)//'": {'// &
                '"observed": false, "interval_start_hour": null, "hour": null, '// &
                '"full_criterion": null}'//suffix
            return
        end if
        write(output_unit,'(A)') '    "'//trim(key)//'": {'
        write(output_unit,'(A)') '      "observed": true,'
        write(output_unit,'(A,ES25.16E3,A)') '      "interval_start_hour": ', &
            step_hour(step-1),','
        write(output_unit,'(A,ES25.16E3,A)') '      "hour": ',step_hour(step),','
        write(output_unit,'(A)') '      "full_criterion": {'
        write(output_unit,'(A,I0,A)') '        "component": ',event%component,','
        write(output_unit,'(A,I0,A)') '        "direction": ',event%direction,','
        write(output_unit,'(A,I0,A)') '        "parent_depth": ',event%parent_depth,','
        write(output_unit,'(A,I0,A)') '        "bfs_iteration": ',event%bfs_iteration,','
        call write_json_array(output_unit,'truncation_estimates',event%errors,.true.,8)
        call write_json_array(output_unit,'component_tolerances',event%tolerances,.true.,8)
        call write_json_array(output_unit,'positive_excess',event%excess,.true.,8)
        call write_json_array(output_unit,'direction_scores',event%direction_scores,.true.,8)
        call write_json_array(output_unit,'parent_unit_center',event%parent_center,.true.,8)
        call write_json_array(output_unit,'parent_unit_width',event%parent_width,.false.,8)
        write(output_unit,'(A)') '      }'
        write(output_unit,'(A)') '    }'//suffix
    end subroutine write_named_split

    subroutine write_step_event(output_unit,key,step,trailing_comma)
        integer,intent(in) :: output_unit,step
        character(len=*),intent(in) :: key
        logical,intent(in) :: trailing_comma
        character(len=1) :: suffix
        suffix=merge(',', ' ',trailing_comma)
        if(step<0) then
            write(output_unit,'(A)') '    "'//trim(key)//'": {'// &
                '"observed": false, "hour": null}'//suffix
        else
            write(output_unit,'(A,ES25.16E3,A)') '    "'//trim(key)//'": {'// &
                '"observed": true, "hour": ',step_hour(step),'}'//suffix
        end if
    end subroutine write_step_event

    subroutine write_failure_event(output_unit,step,previous_pass)
        integer,intent(in) :: output_unit,step,previous_pass
        if(step<0) then
            write(output_unit,'(A)') '    "first_validation_failure": {'// &
                '"observed": false, "hour": null, "previous_passing_hour": null}'
        else if(previous_pass<0) then
            write(output_unit,'(A,ES25.16E3,A)') '    "first_validation_failure": {'// &
                '"observed": true, "hour": ',step_hour(step), &
                ', "previous_passing_hour": null}'
        else
            write(output_unit,'(A,ES25.16E3,A,ES25.16E3,A)') &
                '    "first_validation_failure": {"observed": true, "hour": ', &
                step_hour(step),', "previous_passing_hour": ',step_hour(previous_pass),'}'
        end if
    end subroutine write_failure_event

    subroutine write_json_array(output_unit,key,values,trailing_comma,indent)
        integer,intent(in) :: output_unit,indent
        character(len=*),intent(in) :: key
        real(DP),intent(in) :: values(:)
        logical,intent(in) :: trailing_comma
        integer :: q
        character(len=1) :: suffix
        suffix=merge(',', ' ',trailing_comma)
        write(output_unit,'(A)',advance='no') repeat(' ',indent)//'"'//trim(key)//'": ['
        do q=1,size(values)
            if(q>1) write(output_unit,'(A)',advance='no') ', '
            write(output_unit,'(ES25.16E3)',advance='no') values(q)
        end do
        write(output_unit,'(A)') ']'//suffix
    end subroutine write_json_array

    real(DP) function step_hour(step)
        integer,intent(in) :: step
        if(step<0) then
            step_hour=0.0_DP
        else
            step_hour=min(duration_hours,real(step,DP)*save_hours)
        end if
    end function step_hour

    subroutine usage()
        write(*,'(A)') 'Usage: run_srp_ads_mc_pair -opm INPUT.opm.json -o PREFIX '// &
            '--mc-samples SAMPLES.csv [--mc-count 50000] [--visualization-count 2048] '// &
            '[--eta-sigma 0.3] [--days 15|--hours 360] [--save-hours 1] '// &
            '[--da-order 4] [--max-depth 8] [--pos-tol-km 0.1] '// &
            '[--vel-tol-kms 1e-6] [-cfg config/config.txt]'
    end subroutine usage

    subroutine fatal(reason)
        character(len=*),intent(in) :: reason
        write(*,'(A)') 'ERROR: '//trim(reason)
        stop 1
    end subroutine fatal
end program run_srp_ads_mc_pair
