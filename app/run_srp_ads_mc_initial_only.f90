!> Initial-state-only ADS evaluation on a fixed Monte Carlo design for box-wing SRP.
program run_srp_ads_mc_initial_only
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
    use pod_srp_mc_io_module, only: read_srp_mc_samples, &
        write_initial_only_point_header, write_initial_only_point_rows
    use pod_srp_mc_statistics_module, only: srp_mc_group_stats_type, compute_group_stats
    use pod_spacecraft_geometry, only: srp_geometry_diagnostics_type, &
        reset_srp_geometry_diagnostics, get_srp_geometry_diagnostics
    implicit none

    integer, parameter :: PREFIX_COUNT=3
    integer, parameter :: PREFIX_LABELS(PREFIX_COUNT)=[10000,20000,50000]
    type(hfem_ads_options_type) :: options
    type(hfem_ads_history_type) :: history
    type(srp_mc_group_stats_type) :: stats
    type(hfem_ads_split_event_type) :: first_split
    type(srp_geometry_diagnostics_type) :: geometry_diagnostics
    character(len=MAX_STRING_LEN) :: opm_file,prefix,config_file,sample_file,arg,message,path
    real(DP) :: epoch0,nominal0(6),cov(6,6)
    real(DP) :: pos_tolerance,vel_tolerance,duration_hours,save_hours
    real(DP) :: t_seconds,prev_seconds,final_seconds,current_hour,previous_hour,interval_count_real
    real(DP) :: center(6),previous_center(6),mapped_state(6),mapped_eta
    real(DP) :: validation_pos,validation_vel,max_validation_pos,max_validation_vel,total_elapsed
    real(DP),allocatable :: z(:,:),u(:,:),states(:,:),previous_states(:,:)
    real(DP),allocatable :: probe_points(:,:),probe_ads(:,:),shape(:,:),truth(:,:),eta(:)
    real(DP) :: shape_mean(6),shape_pos_max,shape_vel_max,zero_point(6,1)
    logical,allocatable :: found(:),probe_found(:),shape_found(:),keyframe_saved(:)
    integer,allocatable :: validation_ids(:),all_ids(:),visualization_ids(:)
    integer :: argc,i,j,k,n,actual_count,visualization_count
    integer :: status,step_index,steps,unit_hourly,unit_visual,unit_key,unit_events
    integer :: worst,first_failure_step,previous_pass_step,last_pass_step
    integer :: first_unconfirmed_step,first_split_step,event_count,local_event
    logical :: failed,split_seen,keyframe,diagnostic_keyframe
    logical :: has_opm,has_prefix,has_samples
    character(len=32) :: keyframe_kind,attitude_mode

    opm_file=''
    prefix=''
    sample_file=''
    config_file=''
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
    if(.not.has_opm.or..not.has_prefix.or..not.has_samples.or.len_trim(config_file)==0) then
        call usage()
        call fatal('-opm, -o, --mc-samples, and -cfg are required')
    end if
    if(n<2.or.visualization_count<1.or.duration_hours<0.0_DP.or.save_hours<=0.0_DP) &
        call fatal('MC count must be at least two; duration/cadence/count are invalid')
    if(pos_tolerance<=0.0_DP.or.vel_tolerance<=0.0_DP) call fatal('tolerances must be positive')

    call validate_json_opm_fields(trim(opm_file),status,message)
    if(status/=0) call fatal('invalid JSON OPM: '//trim(message))
    call read_srp_mc_samples(trim(sample_file),n,z,u,status,message)
    if(status/=0) call fatal('MC sample input: '//trim(message))
    call pod_engine_init(trim(config_file))
    if(trim(config%srp_model)/='box_wing') call fatal('requires srp_model=box_wing')
    if(trim(config%srp_attitude_mode)/='earth'.and.trim(config%srp_attitude_mode)/='moon') &
        call fatal('requires earth or moon attitude')
    if(config%srp_scale_da_span/=0.0_DP) call fatal('SRP scale uncertainty must be zero')
    if(any(config%srp_attitude_bias_span_arcsec/=0.0_DP)) &
        call fatal('SRP attitude uncertainty must be zero')
    if(config%srp_array_angle_span_deg/=0.0_DP) &
        call fatal('SRP array-angle uncertainty must be zero')
    config%use_srp=.true.
    attitude_mode=trim(config%srp_attitude_mode)
    call reset_srp_geometry_diagnostics()
    call load_initial_opm(trim(opm_file),epoch0,nominal0,cov)

    options%coordinate_mode=ADS_COORD_WHITENED
    options%srp_sigma=0.0_DP
    options%da_storage_variables=6
    options%error_tolerance=[pos_tolerance,pos_tolerance,pos_tolerance, &
        vel_tolerance,vel_tolerance,vel_tolerance]
    call hfem_ads_history_init(history,nominal0,cov,epoch0,options,status,message)
    if(status/=0) call fatal('ADS init: '//trim(message))
    if(history%coordinate_map%n_variables/=6.or.history%da_variables/=6) &
        call fatal('six-dimensional ADS contract was not established')

    allocate(states(6,n),previous_states(6,n),found(n),all_ids(n))
    all_ids=[(i,i=1,n)]
    call choose_visualization_ids(n,min(visualization_count,n),visualization_ids)
    call build_study_points(probe_points,validation_ids)
    allocate(probe_ads(6,VALIDATION_COUNT),shape(6,SHAPE_COUNT),shape_found(SHAPE_COUNT), &
        truth(6,VALIDATION_COUNT),eta(VALIDATION_COUNT),probe_found(VALIDATION_COUNT))
    do j=1,VALIDATION_COUNT
        call ads_unit_to_physical(history%coordinate_map,probe_points(1:6,validation_ids(j)), &
            mapped_state,mapped_eta,status)
        if(status/=0) call fatal('validation mapping failed')
        truth(:,j)=nominal0+mapped_state
        eta(j)=0.0_DP
    end do

    call open_outputs()
    first_failure_step=-1
    previous_pass_step=-1
    last_pass_step=-1
    first_unconfirmed_step=-1
    first_split_step=-1
    event_count=0
    split_seen=.false.
    max_validation_pos=0.0_DP
    max_validation_vel=0.0_DP
    total_elapsed=0.0_DP

    final_seconds=duration_hours*3600.0_DP
    interval_count_real=duration_hours/save_hours
    if(abs(interval_count_real-real(nint(interval_count_real),DP))<= &
       1.0e-12_DP*max(1.0_DP,abs(interval_count_real))) then
        steps=nint(interval_count_real)
    else
        steps=ceiling(interval_count_real)
    end if
    allocate(keyframe_saved(0:steps),source=.false.)
    previous_states=0.0_DP
    previous_center=nominal0
    prev_seconds=0.0_DP

    do step_index=0,steps
        t_seconds=min(final_seconds,real(step_index,DP)*save_hours*3600.0_DP)
        current_hour=t_seconds/3600.0_DP
        previous_hour=prev_seconds/3600.0_DP
        if(step_index>0) then
            call hfem_ads_history_advance(history,t_seconds,status,message)
            if(status/=0) call fatal('ADS advance: '//trim(message))
            total_elapsed=total_elapsed+history%stats%elapsed_seconds
            call capture_events()
            call advance_real_probes(truth,eta,epoch0,prev_seconds,t_seconds,options,status,message)
            if(status/=0) call fatal('real validation advance: '//trim(message))
        end if

        call hfem_ads_history_evaluate(history,u(1:6,:),states,found,status)
        if(status/=0.or..not.all(found)) call fatal('MC evaluation failed')
        zero_point=0.0_DP
        call hfem_ads_history_evaluate(history,zero_point,probe_ads(:,1:1),probe_found(1:1),status)
        if(status/=0.or..not.probe_found(1)) call fatal('center evaluation failed')
        center=probe_ads(:,1)
        call hfem_ads_history_evaluate(history,probe_points(1:6,:),shape,shape_found,status)
        if(status/=0.or..not.all(shape_found)) call fatal('shape evaluation failed')
        probe_ads=shape(:,validation_ids)
        call assess_accuracy(probe_ads,truth,pos_tolerance,vel_tolerance,validation_pos, &
            validation_vel,worst,failed)
        call compute_shape_summary(shape,center,shape_mean,shape_pos_max,shape_vel_max)
        max_validation_pos=max(max_validation_pos,validation_pos)
        max_validation_vel=max(max_validation_vel,validation_vel)
        call update_failure(failed)
        if(history%stats%depth_limited_patches>0.and.first_unconfirmed_step<0) &
            first_unconfirmed_step=step_index

        do k=1,PREFIX_COUNT
            actual_count=min(PREFIX_LABELS(k),n)
            call compute_group_stats(states(:,1:actual_count),center,stats,status)
            if(status/=0) call fatal('statistics failed')
            call write_group_row(PREFIX_LABELS(k),actual_count)
        end do

        call write_initial_only_point_rows(unit_visual,current_hour,'hourly_visualization', &
            visualization_ids,states(:,visualization_ids),center,status,message)
        if(status/=0) call fatal('visualization point output: '//trim(message))

        diagnostic_keyframe=step_index==first_split_step.or. &
            step_index==first_unconfirmed_step.or.step_index==first_failure_step
        keyframe=is_keyframe(current_hour,step_index,steps).or.diagnostic_keyframe
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
                call write_initial_only_point_rows(unit_key,previous_hour,'diagnostic_previous', &
                    all_ids,previous_states,previous_center,status,message)
                if(status/=0) call fatal('previous keyframe point output: '//trim(message))
                keyframe_saved(step_index-1)=.true.
                call write_snapshot(.true.,nint(prev_seconds))
            end if
        end if
        if(keyframe.and..not.keyframe_saved(step_index)) then
            call write_initial_only_point_rows(unit_key,current_hour,trim(keyframe_kind), &
                all_ids,states,center,status,message)
            if(status/=0) call fatal('keyframe point output: '//trim(message))
            keyframe_saved(step_index)=.true.
            call write_snapshot(.false.,nint(t_seconds))
        end if
        flush(unit_hourly)
        flush(unit_visual)
        flush(unit_key)
        write(*,'(A,F8.2,A,I0,A,ES10.3,A,ES10.3)') 'hour ',current_hour, &
            ' patches ',history%domain%n_patches,' validation_r_km ',validation_pos, &
            ' validation_v_kms ',validation_vel
        previous_states=states
        previous_center=center
        prev_seconds=t_seconds
    end do

    close(unit_hourly)
    close(unit_visual)
    close(unit_key)
    close(unit_events)
    call write_report()
    call hfem_ads_history_destroy(history)

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

    subroutine update_failure(is_failed)
        logical,intent(in) :: is_failed
        if(is_failed.and.first_failure_step<0) then
            first_failure_step=step_index
            previous_pass_step=last_pass_step
        end if
        if(.not.is_failed) last_pass_step=step_index
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
        open(newunit=unit_visual,file=trim(prefix)//'_mc_visualization_points.csv', &
            status='replace',action='write',iostat=status)
        if(status/=0) call fatal('cannot open visualization output')
        open(newunit=unit_key,file=trim(prefix)//'_mc_keyframe_points.csv', &
            status='replace',action='write',iostat=status)
        if(status/=0) call fatal('cannot open keyframe output')
        open(newunit=unit_events,file=trim(prefix)//'_split_events.csv', &
            status='replace',action='write',iostat=status)
        if(status/=0) call fatal('cannot open split events output')

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
        call write_initial_only_point_header(unit_visual,status,message)
        if(status/=0) call fatal(trim(message))
        call write_initial_only_point_header(unit_key,status,message)
        if(status/=0) call fatal(trim(message))
        write(unit_events,'(A)') EVENT_HEADER
    end subroutine open_outputs

    subroutine write_group_row(label,count_value)
        integer,intent(in) :: label,count_value
        integer :: split7(7),maximum_depth,accuracy_valid,precision_confirmed
        integer :: skew_component,kurtosis_component
        accuracy_valid=merge(1,0,first_failure_step<0.or.step_index<first_failure_step)
        precision_confirmed=merge(1,0,first_unconfirmed_step<0.or. &
            step_index<first_unconfirmed_step)
        maximum_depth=max_patch_depth()
        skew_component=maxloc(abs(stats%skewness),dim=1)
        kurtosis_component=maxloc(abs(stats%excess_kurtosis),dim=1)
        split7=0
        split7(1:size(history%stats%split_counts))=history%stats%split_counts
        write(unit_hourly,'(*(G0,:,","))') current_hour,epoch0+t_seconds,trim(attitude_mode), &
            label,count_value,history%coordinate_map%n_variables,history%domain%n_patches, &
            history%stats%depth_limited_patches,split7,history%stats%bfs_iterations, &
            history%stats%max_queue_size,history%stats%elapsed_seconds,validation_pos, &
            validation_vel,worst,merge(1,0,failed),accuracy_valid,precision_confirmed, &
            maximum_depth,center,stats%mean_error,stats%standard_deviation, &
            stats%component_min,stats%component_max,stats%skewness,stats%excess_kurtosis, &
            stats%position_axes_rms,stats%velocity_axes_rms,stats%position_norm_rms, &
            stats%velocity_norm_rms,stats%position_norm_quantiles,stats%velocity_norm_quantiles, &
            stats%position_norm_max,stats%velocity_norm_max,stats%position_worst_id, &
            stats%velocity_worst_id,stats%position_max_to_longest_axis, &
            stats%velocity_max_to_longest_axis,shape_mean,shape_pos_max,shape_vel_max, &
            safe_ratio(shape_pos_max,stats%position_axes_rms(1)), &
            safe_ratio(shape_vel_max,stats%velocity_axes_rms(1)), &
            maxval(abs(stats%skewness)),skew_component, &
            maxval(abs(stats%excess_kurtosis)),kurtosis_component, &
            stats%covariance(1,1),stats%covariance(1,2),stats%covariance(1,3), &
            stats%covariance(2,2),stats%covariance(2,3),stats%covariance(3,3), &
            stats%covariance(4,4),stats%covariance(4,5),stats%covariance(4,6), &
            stats%covariance(5,5),stats%covariance(5,6),stats%covariance(6,6)
    end subroutine write_group_row

    subroutine compute_shape_summary(input_states,input_center,mean_error,pos_max,vel_max)
        real(DP),intent(in) :: input_states(:,:),input_center(6)
        real(DP),intent(out) :: mean_error(6),pos_max,vel_max
        real(DP),allocatable :: errors(:,:)
        allocate(errors(6,size(input_states,2)))
        errors=input_states-spread(input_center,dim=2,ncopies=size(input_states,2))
        mean_error=sum(errors,dim=2)/real(size(input_states,2),DP)
        pos_max=maxval(sqrt(sum(errors(1:3,:)**2,dim=1)))
        vel_max=maxval(sqrt(sum(errors(4:6,:)**2,dim=1)))
    end subroutine compute_shape_summary

    integer function max_patch_depth() result(depth)
        integer :: q
        depth=0
        do q=1,history%domain%n_patches
            depth=max(depth,sh_count(history%domain%patches(q)%history,0))
        end do
    end function max_patch_depth

    real(DP) function safe_ratio(numerator,denominator) result(value)
        real(DP),intent(in) :: numerator,denominator
        if(abs(denominator)>100.0_DP*epsilon(1.0_DP)) then
            value=numerator/denominator
        else if(abs(numerator)<=100.0_DP*epsilon(1.0_DP)) then
            value=0.0_DP
        else
            value=sign(huge(1.0_DP),numerator)
        end if
    end function safe_ratio

    subroutine capture_events()
        if(.not.allocated(history%stats%split_events)) return
        do local_event=1,size(history%stats%split_events)
            event_count=event_count+1
            call write_split_event(history%stats%split_events(local_event),local_event,event_count)
            if(.not.split_seen) then
                split_seen=.true.
                first_split_step=step_index
                first_split=history%stats%split_events(local_event)
            end if
        end do
    end subroutine capture_events

    subroutine write_split_event(event,local_index,global_index)
        type(hfem_ads_split_event_type),intent(in) :: event
        integer,intent(in) :: local_index,global_index
        write(unit_events,'(*(G0,:,","))') previous_hour,current_hour,global_index,local_index, &
            event%bfs_iteration,event%parent_depth,event%component,event%direction,event%errors, &
            event%tolerances,event%excess,event%direction_scores(1:7),event%parent_center(1:7),event%parent_width(1:7)
    end subroutine write_split_event

    subroutine write_snapshot(use_previous,time_seconds)
        logical,intent(in) :: use_previous
        integer,intent(in) :: time_seconds
        character(len=16) :: tag
        write(tag,'(I7.7)') time_seconds
        path=trim(prefix)//'_patch_s'//trim(tag)
        call write_ads_patch_snapshot(history,trim(path),use_previous,status,message)
        if(status/=0) call fatal('patch snapshot: '//trim(message))
    end subroutine write_snapshot

    subroutine write_report()
        integer :: report_unit,ios,q,written
        call get_srp_geometry_diagnostics(geometry_diagnostics)
        open(newunit=report_unit,file=trim(prefix)//'_mc_report.json', &
            status='replace',action='write',iostat=ios)
        if(ios/=0) call fatal('cannot open MC JSON report')
        write(report_unit,'(A)') '{'
        write(report_unit,'(A)') '  "output_schema_version": "srp_ads_mc_initial_only_v1",'
        write(report_unit,'(A)') '  "method": "fixed-design initial-error-only ADS Monte Carlo evaluation",'
        write(report_unit,'(A)') '  "nominal_srp_enabled": true,'
        write(report_unit,'(A)') '  "srp_parameter_uncertainty_enabled": false,'
        write(report_unit,'(A)') '  "srp_model": "'//trim(config%srp_model)//'",'
        write(report_unit,'(A)') '  "srp_attitude_mode": "'//trim(attitude_mode)//'",'
        write(report_unit,'(A,I0,A)') '  "active_dimensions": ', &
            history%coordinate_map%n_variables,','
        write(report_unit,'(A,I0,A)') '  "dace_storage_dimensions": ',history%da_variables,','
        write(report_unit,'(A,ES25.16E3,A)') '  "srp_eta": ',0.0_DP,','
        write(report_unit,'(A,I0,A)') '  "sample_count": ',n,','
        write(report_unit,'(A,I0,A)') '  "visualization_sample_count": ', &
            size(visualization_ids),','
        write(report_unit,'(A,I0,A)') '  "validation_point_count": ',VALIDATION_COUNT,','
        write(report_unit,'(A,I0,A)') '  "shape_point_count": ',SHAPE_COUNT,','
        write(report_unit,'(A,ES25.16E3,A)') '  "duration_hours": ',duration_hours,','
        write(report_unit,'(A,ES25.16E3,A)') '  "save_hours": ',save_hours,','
        write(report_unit,'(A,I0,A)') '  "da_order": ',options%da_order,','
        write(report_unit,'(A,I0,A)') '  "max_split_depth": ',options%max_split_depth,','
        write(report_unit,'(A,ES25.16E3,A)') '  "position_tolerance_km": ',pos_tolerance,','
        write(report_unit,'(A,ES25.16E3,A)') '  "velocity_tolerance_kms": ',vel_tolerance,','
        write(report_unit,'(A,ES25.16E3,A)') '  "maximum_validation_position_error_km": ', &
            max_validation_pos,','
        write(report_unit,'(A,ES25.16E3,A)') '  "maximum_validation_velocity_error_kms": ', &
            max_validation_vel,','
        write(report_unit,'(A,ES25.16E3,A)') '  "ads_propagation_elapsed_seconds": ', &
            total_elapsed,','
        write(report_unit,'(A)') '  "sample_statistics_interpretation": '// &
            '"fixed truncated-normal design metrics, not probability confidence intervals",'
        write(report_unit,'(A)') '  "split_decision_rule": '// &
            '"split when any component truncation estimate exceeds its tolerance; '// &
            'component maximizes positive excess and direction maximizes top-order score",'
        write(report_unit,'(A)') '  "opm_file": "'//trim(opm_file)//'",'
        write(report_unit,'(A)') '  "config_file": "'//trim(config_file)//'",'
        write(report_unit,'(A)') '  "sample_file": "'//trim(sample_file)//'",'
        write(report_unit,'(A)') '  "geometry_fallback_diagnostics": {'
        write(report_unit,'(A,I0,A)') '    "real_count": ', &
            geometry_diagnostics%real_fallback_count,','
        write(report_unit,'(A,I0,A)') '    "da_count": ', &
            geometry_diagnostics%da_fallback_count,','
        write(report_unit,'(A)') '    "first_real_message": "'// &
            trim(geometry_diagnostics%first_real_message)//'",'
        write(report_unit,'(A)') '    "first_da_message": "'// &
            trim(geometry_diagnostics%first_da_message)//'"'
        write(report_unit,'(A)') '  },'
        write(report_unit,'(A)',advance='no') '  "saved_keyframe_hours": ['
        written=0
        do q=0,steps
            if(.not.keyframe_saved(q)) cycle
            if(written>0) write(report_unit,'(A)',advance='no') ', '
            write(report_unit,'(ES25.16E3)',advance='no') step_hour(q)
            written=written+1
        end do
        write(report_unit,'(A)') '],'
        write(report_unit,'(A)') '  "events": {'
        call write_first_split(report_unit)
        write(report_unit,'(A)') ','
        call write_step_event(report_unit,'first_unconfirmed',first_unconfirmed_step,.true.)
        call write_failure_event(report_unit,first_failure_step,previous_pass_step)
        write(report_unit,'(A)') '  }'
        write(report_unit,'(A)') '}'
        close(report_unit)
    end subroutine write_report

    subroutine write_first_split(output_unit)
        integer,intent(in) :: output_unit
        if(.not.split_seen) then
            write(output_unit,'(A)',advance='no') '    "first_split": {'// &
                '"observed": false, "interval_start_hour": null, "hour": null, '// &
                '"full_criterion": null}'
            return
        end if
        write(output_unit,'(A)') '    "first_split": {'
        write(output_unit,'(A)') '      "observed": true,'
        write(output_unit,'(A,ES25.16E3,A)') '      "interval_start_hour": ', &
            step_hour(first_split_step-1),','
        write(output_unit,'(A,ES25.16E3,A)') '      "hour": ',step_hour(first_split_step),','
        write(output_unit,'(A)') '      "full_criterion": {'
        write(output_unit,'(A,I0,A)') '        "component": ',first_split%component,','
        write(output_unit,'(A,I0,A)') '        "direction": ',first_split%direction,','
        write(output_unit,'(A,I0,A)') '        "parent_depth": ',first_split%parent_depth,','
        write(output_unit,'(A,I0,A)') '        "bfs_iteration": ',first_split%bfs_iteration,','
        call write_json_array(output_unit,'truncation_estimates',first_split%errors,.true.,8)
        call write_json_array(output_unit,'component_tolerances',first_split%tolerances,.true.,8)
        call write_json_array(output_unit,'positive_excess',first_split%excess,.true.,8)
        call write_json_array(output_unit,'direction_scores',first_split%direction_scores(1:7),.true.,8)
        call write_json_array(output_unit,'parent_unit_center',first_split%parent_center(1:7),.true.,8)
        call write_json_array(output_unit,'parent_unit_width',first_split%parent_width(1:7),.false.,8)
        write(output_unit,'(A)') '      }'
        write(output_unit,'(A)',advance='no') '    }'
    end subroutine write_first_split

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
        write(*,'(A)') 'Usage: run_srp_ads_mc_initial_only -opm INPUT.opm.json -o PREFIX '// &
            '--mc-samples SAMPLES.csv -cfg BOX_WING_CONFIG [--mc-count 50000] '// &
            '[--visualization-count 2048] [--days 15|--hours 360] [--save-hours 1] '// &
            '[--da-order 4] [--max-depth 8] [--pos-tol-km 0.1] [--vel-tol-kms 1e-6]'
    end subroutine usage

    subroutine fatal(reason)
        character(len=*),intent(in) :: reason
        write(*,'(A)') 'ERROR: '//trim(reason)
        stop 1
    end subroutine fatal
end program run_srp_ads_mc_initial_only
