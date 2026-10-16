!> Separate C2 pipeline. Raw results are never labelled formally audited.
program run_srp_modes_ads
    use pod_global, only: DP,MAX_STRING_LEN
    use pod_config, only: config,validate_config
    use pod_engine_module, only: pod_engine_init
    use pod_data_format_module, only: load_initial_opm
    use pod_json_object_module, only: json_string
    use pod_srp_mode_module
    use pod_uq_ads_coordinates_module, only: ADS_COORD_WHITENED,ADS_PARAM_SRP_MODES, &
        ads_physical_parameters_type,ads_unit_to_physical
    use pod_uq_hfem_ads_module
    use pod_srp_ads_history_module, only: build_study_points,advance_real_mode_probes,advance_real_probes, &
        assess_accuracy,compute_shape_metrics,validate_json_opm_fields
    use pod_srp_ads_snapshot_module, only: write_ads_patch_snapshot
    use pod_srp_modes_probability_io_module
    use pod_ads_split_module, only: mf_find_patch
    use pod_force_model_module, only: clear_srp_scale_error,clear_srp_mode_coefficients_real
    use pod_da_force_model_module, only: clear_srp_scale_uncertainty,clear_srp_mode_uncertainty
    use pod_dace_classes, only: active_da_count
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    type(hfem_ads_history_type) :: history
    type(hfem_ads_options_type) :: options
    type(ads_physical_parameters_type) :: parameters
    type(probability_population_type) :: population
    type(probability_output_type) :: probability_output
    character(len=MAX_STRING_LEN) :: probability_samples=''
    integer :: probability_count=50000,probability_visual_count=2000,probability_probe_count=64,probability_probe_unit=0
    integer,allocatable :: probability_ids(:),probability_owners(:)
    real(DP),allocatable :: probability_states(:,:),probability_truth(:,:),probability_s(:),probability_u(:)
    logical,allocatable :: probability_found(:)
    logical :: probability_enabled=.false.,probability_failed=.false.,probability_step_failed
    real(DP) :: probability_pos=0,probability_vel=0,probability_max_pos=0,probability_max_vel=0
    character(len=MAX_STRING_LEN) :: opm_file='',cfg_file='',prefix='',arg,message,filename,revision='unrecorded'
    character(len=64) :: csv_hash,meta_hash='UNVERIFIED'
    character(len=:),allocatable :: metadata,header
    character(len=MAX_STRING_LEN),allocatable :: snapshots(:)
    real(DP) :: initial(6),cov(6,6),epoch,duration=360.0_DP,cadence=1.0_DP,pos_tol=0.1_DP,vel_tol=1e-6_DP
    real(DP) :: previous,t,nominal(6),deviation(6),extras(8,4),extra_states(6,4)
    real(DP) :: lower(6),upper(6),pos_max,vel_max,pos_axes(3),vel_axes(3)
    real(DP) :: max_pos=0,max_vel=0,step_pos,step_vel,zero_pos=0,zero_vel=0,zero_real_pos=0,zero_real_vel=0
    real(DP) :: c0(6,1),zero_real(6,1),offset_pos(2),offset_vel(2)
    real(DP),allocatable :: points(:,:),ads(:,:),truth(:,:),s(:),u(:)
    integer,allocatable :: validation_ids(:)
    logical,allocatable :: found(:)
    logical :: exists,failed,any_failed=.false.,unconfirmed=.false.,lock_open=.false.,extra_found(4)
    integer :: i,j,k,status,n_steps,n_validation,unit_initial,unit_history,unit_shape,unit_probes,unit_events,unit_lock
    integer :: before,after,unit,bytes,step,worst,total_events=0,snapshot_count=0,depth_count=0,patches=0
    character(len=16) :: tag
    i=1
    do while(i<=command_argument_count())
        call get_command_argument(i,arg,status=status)
        if(status/=0) call fatal('CLI argument too long')
        if(trim(arg)=='--help' .or. trim(arg)=='-h') then
            write(*,'(A)') '8D C2 modes ADS: -opm INPUT -cfg CONFIG -o PREFIX'
            write(*,'(A)') '--hours 360 --save-hours 1 --da-order 4 --max-depth 8 --pos-tol-km 0.1 --vel-tol-kms 1e-6'
            write(*,'(A)') '--rel-tol 1e-12 --abs-tol 1e-12'
            write(*,'(A)') '--probability-samples CSV --probability-count 50000 --probability-visualization-count 2000'
            write(*,'(A)') '--probability-validation-count 64 (independent Real probes; no change to original 337)'
            write(*,'(A)') 'Metadata SHA-256 uses sha256sum; optional --mode-meta-sha256 verifies the expected identity.'
            write(*,'(A)') 'Formal runs use tools/run_srp_modes_ads.py.'
            stop
        end if
        if(i==command_argument_count()) call fatal('option requires value: '//trim(arg))
        select case(trim(arg))
        case('-opm'); call get_command_argument(i+1,opm_file,status=status)
        case('-cfg'); call get_command_argument(i+1,cfg_file,status=status)
        case('-o'); call get_command_argument(i+1,prefix,status=status)
        case('--mode-meta-sha256'); call get_command_argument(i+1,meta_hash,status=status)
        case('--source-revision'); call get_command_argument(i+1,revision,status=status)
        case('--probability-samples'); call get_command_argument(i+1,probability_samples,status=status)
        case('--probability-count','--probability-visualization-count','--probability-validation-count')
            filename=arg
            call get_command_argument(i+1,arg,status=status)
            if(status/=0) call fatal('probability CLI value too long')
            select case(trim(filename))
            case('--probability-count'); read(arg,*,iostat=status) probability_count
            case('--probability-visualization-count'); read(arg,*,iostat=status) probability_visual_count
            case('--probability-validation-count'); read(arg,*,iostat=status) probability_probe_count
            end select
        case('--hours','--save-hours','--da-order','--max-depth','--pos-tol-km','--vel-tol-kms','--rel-tol','--abs-tol')
            filename=arg
            call get_command_argument(i+1,arg,status=status)
            if(status/=0) call fatal('CLI value too long')
            select case(trim(filename))
            case('--hours'); read(arg,*,iostat=status) duration
            case('--save-hours'); read(arg,*,iostat=status) cadence
            case('--da-order'); read(arg,*,iostat=status) options%da_order
            case('--max-depth'); read(arg,*,iostat=status) options%max_split_depth
            case('--pos-tol-km'); read(arg,*,iostat=status) pos_tol
            case('--vel-tol-kms'); read(arg,*,iostat=status) vel_tol
            case('--rel-tol'); read(arg,*,iostat=status) options%rel_tol
            case('--abs-tol'); read(arg,*,iostat=status) options%abs_tol
            end select
        case default
            call fatal('unknown option: '//trim(arg))
        end select
        if(status/=0) call fatal('invalid CLI value')
        i=i+2
    end do
    if(min(len_trim(opm_file),len_trim(cfg_file),len_trim(prefix))==0) call fatal('-opm, -cfg, -o are required')
    probability_enabled=len_trim(probability_samples)>0
    if(probability_enabled) then
        if(probability_count<2.or.min(probability_visual_count,probability_probe_count)<1) &
            call fatal('invalid probability sample/visualization/probe count')
        call load_probability_samples(trim(probability_samples),probability_count,population,status)
        if(status/=0) call fatal('invalid probability sample CSV/schema/count/mapping')
    end if
    if(.not.all(ieee_is_finite([duration,cadence,pos_tol,vel_tol,options%rel_tol,options%abs_tol]))) &
        call fatal('CLI numbers must be finite')
    if(min(duration,cadence,pos_tol,vel_tol,options%rel_tol,options%abs_tol)<=0 .or. &
        options%da_order<1 .or. options%max_split_depth<0) &
        call fatal('invalid positive limit')
    if(duration/cadence>10000 .or. len_trim(prefix)>MAX_STRING_LEN-45) call fatal('output size/path limit exceeded')
    n_steps=nint(duration/cadence)
    if(n_steps<1 .or. abs(n_steps*cadence-duration)>1e-10_DP) call fatal('duration must contain whole save intervals')
    call validate_json_opm_fields(trim(opm_file),status,message)
    if(status/=0) call fatal(trim(message))
    call pod_engine_init(trim(cfg_file))
    if(.not.validate_config()) call fatal('invalid configuration')
    if(trim(config%srp_model)/='time_varying_modes' .or. .not.config%use_srp) call fatal('C2 requires time_varying_modes SRP')
    call load_initial_opm(trim(opm_file),epoch,initial,cov)
    call load_srp_mode_table(trim(config%srp_mode_table_file),trim(config%srp_mode_metadata_file),epoch,status,message)
    if(status/=0) call fatal(trim(message))
    if(duration>360.0_DP) then
        block
            real(DP) :: eE(3),eM(3)
            call interpolate_srp_modes(epoch+duration*3600,eE,eM,status,message)
            if(status/=0) call fatal(trim(message))
        end block
    end if
    open(newunit=unit,file=trim(config%srp_mode_metadata_file),access='stream',form='unformatted',status='old')
    inquire(unit=unit,size=bytes); allocate(character(len=bytes) :: metadata)
    read(unit) metadata; close(unit)
    call json_string(metadata,'csv_sha256',header,exists)
    if(.not.exists) call fatal('metadata lacks CSV hash')
    csv_hash=header
    inquire(file=trim(prefix)//'_report.json',exist=exists)
    if(exists) call fatal('output already exists')
    open(newunit=unit_lock,file=trim(prefix)//'.lock',status='new',iostat=status)
    if(status/=0) call fatal('output already exists or prefix unavailable')
    lock_open=.true.
    call bind_metadata_hash()
    options%parameterization=ADS_PARAM_SRP_MODES; options%coordinate_mode=ADS_COORD_WHITENED
    options%error_tolerance=[pos_tol,pos_tol,pos_tol,vel_tol,vel_tol,vel_tol]
    before=active_da_count()
    call hfem_ads_history_init(history,initial,cov,epoch,options,status,message)
    if(status/=0) call fatal(trim(message))
    if(probability_enabled) then
        call probability_choose_ids(probability_count,probability_probe_count,probability_ids,status)
        if(status/=0) call fatal('invalid probability probe selection')
        probability_probe_count=size(probability_ids)
        allocate(probability_states(6,probability_count),probability_found(probability_count), &
            probability_owners(probability_count),probability_truth(6,probability_probe_count), &
            probability_s(probability_probe_count),probability_u(probability_probe_count))
        do j=1,probability_probe_count
            k=probability_ids(j)
            call ads_unit_to_physical(history%coordinate_map,population%q(:,k),deviation,parameters,status)
            if(status/=0) call fatal('probability initial map failed')
            probability_truth(:,j)=initial+deviation
            probability_s(j)=population%modes(1,k); probability_u(j)=population%modes(2,k)
        end do
        call probability_open(probability_output,trim(prefix),population,probability_visual_count,status)
        if(status/=0) call fatal('cannot open new probability outputs')
        open(newunit=probability_probe_unit,file=trim(prefix)//'_probability_probes.csv',status='new')
        write(probability_probe_unit,'(A)') 'hour,sample_id,s,u,ads_x,ads_y,ads_z,ads_vx,ads_vy,ads_vz,'// &
            'real_x,real_y,real_z,real_vx,real_vy,real_vz,ads_numerical_error_pos_km,ads_numerical_error_vel_kms'
    end if
    call build_study_points(points,validation_ids,8)
    n_validation=size(validation_ids)
    allocate(ads(6,size(points,2)),truth(6,n_validation),s(n_validation),u(n_validation),found(size(points,2)))
    allocate(snapshots(2*n_steps+2))
    open(newunit=unit_initial,file=trim(prefix)//'_initial_points.csv',status='new')
    write(unit_initial,'(A)') 'point_id,u1,u2,u3,u4,u5,u6,u7,u8,s,u,x0_km,y0_km,z0_km,vx0_kms,vy0_kms,vz0_kms,validation'
    do j=1,size(points,2)
        call ads_unit_to_physical(history%coordinate_map,points(:,j),deviation,parameters,status)
        if(status/=0) call fatal('study point mapping failed')
        write(unit_initial,'(I0,16(",",ES25.16E3),",",I0)') j,points(:,j),parameters%s,parameters%u, &
            initial+deviation,merge(1,0,j<=n_validation)
        if(j<=n_validation) then
            truth(:,j)=initial+deviation; s(j)=parameters%s; u(j)=parameters%u
        end if
    end do
    close(unit_initial)
    open(newunit=unit_history,file=trim(prefix)//'_history.csv',status='new')
    open(newunit=unit_shape,file=trim(prefix)//'_shape.csv',status='new')
    open(newunit=unit_probes,file=trim(prefix)//'_probes.csv',status='new')
    open(newunit=unit_events,file=trim(prefix)//'_split_events.csv',status='new')
    header='hour,epoch_et,nominal_x_km,nominal_y_km,nominal_z_km,nominal_vx_kms,nominal_vy_kms,nominal_vz_kms'
    do k=1,6
        write(tag,'(I0)') k; header=header//',lower'//trim(tag)
    end do
    do k=1,6
        write(tag,'(I0)') k; header=header//',upper'//trim(tag)
    end do
    header=header//',sampled_pos_max_km,sampled_vel_max_kms,patches,depth_limited'
    do k=1,8
        write(tag,'(I0)') k; header=header//',split_dim'//trim(tag)
    end do
    write(unit_history,'(A)') header//',ads_numerical_error_pos_km,ads_numerical_error_vel_kms,zero_mode_pos_km,zero_mode_vel_kms'
    write(unit_shape,'(A)') 'hour,point_id,dx_km,dy_km,dz_km,dvx_kms,dvy_kms,dvz_kms'
    write(unit_probes,'(A)') 'hour,point_id,s,u,ads_x,ads_y,ads_z,ads_vx,ads_vy,ads_vz,'// &
        'real_x,real_y,real_z,real_vx,real_vy,real_vz,ads_numerical_error_pos_km,ads_numerical_error_vel_kms'
    header='interval_start_hour,interval_end_hour,event_id,bfs_iteration,parent_depth,component,direction'
    do i=1,6
        write(tag,'(I0)') i; header=header//',error'//trim(tag)
    end do
    do i=1,6
        write(tag,'(I0)') i; header=header//',tolerance'//trim(tag)
    end do
    do i=1,6
        write(tag,'(I0)') i; header=header//',excess'//trim(tag)
    end do
    do i=1,8
        write(tag,'(I0)') i; header=header//',direction_score'//trim(tag)
    end do
    do i=1,8
        write(tag,'(I0)') i; header=header//',parent_center'//trim(tag)
    end do
    do i=1,8
        write(tag,'(I0)') i; header=header//',parent_width'//trim(tag)
    end do
    write(unit_events,'(A)') header
    c0(:,1)=initial; zero_real=c0; previous=0
    extras=0; extras(7:8,2)=-1; extras(7,3)=1; extras(7,4)=-1; extras(8,4)=1
    do step=0,n_steps
        t=real(step,DP)*cadence*3600
        if(step>0) then
            call hfem_ads_history_advance(history,t,status,message)
            if(status/=0) call fatal(trim(message))
            call advance_real_mode_probes(truth,s,u,epoch,previous,t,options,status,message)
            if(status/=0) call fatal(trim(message))
            call advance_real_mode_probes(zero_real,[0.0_DP],[0.0_DP],epoch,previous,t,options,status,message)
            if(status/=0) call fatal(trim(message))
            config%srp_model='cannonball'
            call advance_real_probes(c0,[0.0_DP],epoch,previous,t,options,status,message)
            config%srp_model='time_varying_modes'
            if(status/=0) call fatal(trim(message))
            if(probability_enabled) then
                call advance_real_mode_probes(probability_truth,probability_s,probability_u,epoch,previous,t,options,status,message)
                if(status/=0) call fatal(trim(message))
            end if
        end if
        nominal=c0(:,1)
        call hfem_ads_history_evaluate(history,points,ads,found,status)
        if(status/=0 .or. .not.all(found) .or. .not.all(ieee_is_finite(ads))) call fatal('invalid ADS study evaluation')
        call assess_accuracy(ads(:,validation_ids),truth,pos_tol,vel_tol,step_pos,step_vel,worst,failed)
        any_failed=any_failed.or.failed; max_pos=max(max_pos,step_pos); max_vel=max(max_vel,step_vel)
        call hfem_ads_history_evaluate(history,extras,extra_states,extra_found,status)
        if(status/=0 .or. .not.all(extra_found)) call fatal('invalid endpoint evaluation')
        zero_pos=max(zero_pos,norm2(extra_states(1:3,2)-nominal(1:3)))
        zero_vel=max(zero_vel,norm2(extra_states(4:6,2)-nominal(4:6)))
        zero_real_pos=max(zero_real_pos,norm2(zero_real(1:3,1)-nominal(1:3)))
        zero_real_vel=max(zero_real_vel,norm2(zero_real(4:6,1)-nominal(4:6)))
        do j=1,2
            offset_pos(j)=norm2(extra_states(1:3,j+2)-nominal(1:3))
            offset_vel(j)=norm2(extra_states(4:6,j+2)-nominal(4:6))
        end do
        call compute_shape_metrics(ads,nominal,lower,upper,pos_max,vel_max,pos_axes,vel_axes,status)
        if(status/=0) call fatal('invalid sampled envelope')
        patches=history%domain%n_patches
        depth_count=depth_count+history%stats%depth_limited_patches
        unconfirmed=unconfirmed.or.history%stats%depth_limited_patches>0
        if(probability_enabled) then
            call hfem_ads_history_evaluate(history,population%q,probability_states,probability_found,status)
            if(status/=0.or..not.all(probability_found).or..not.all(ieee_is_finite(probability_states))) &
                call fatal('invalid full probability population evaluation; no sample may be dropped')
            do j=1,probability_count
                call mf_find_patch(history%domain,population%q(:,j),probability_owners(j))
            end do
            if(any(probability_owners<1)) call fatal('probability point has no owning patch')
            call assess_accuracy(probability_states(:,probability_ids),probability_truth,pos_tol,vel_tol, &
                probability_pos,probability_vel,worst,probability_step_failed)
            probability_failed=probability_failed.or.probability_step_failed
            probability_max_pos=max(probability_max_pos,probability_pos)
            probability_max_vel=max(probability_max_vel,probability_vel)
            do j=1,probability_probe_count
                k=probability_ids(j)
                write(probability_probe_unit,'(ES25.16E3,",",I0,16(",",ES25.16E3))') t/3600,k, &
                    probability_s(j),probability_u(j),probability_states(:,k),probability_truth(:,j), &
                    norm2(probability_states(1:3,k)-probability_truth(1:3,j)), &
                    norm2(probability_states(4:6,k)-probability_truth(4:6,j))
            end do
            call probability_write_hour(probability_output,population,t/3600,epoch+t,probability_states,nominal, &
                probability_owners,.not.(any_failed.or.probability_failed.or.zero_pos>1e-6_DP.or.zero_vel>1e-9_DP), &
                .not.unconfirmed,history%stats%split_occurred,history%stats%depth_limited_patches>0, &
                any_failed.or.probability_failed.or.zero_pos>1e-6_DP.or.zero_vel>1e-9_DP,patches, &
                history%stats%depth_limited_patches,step_pos,step_vel,probability_pos,probability_vel,step==n_steps,status)
            if(status/=0) call fatal('probability full-population statistics/output failed')
        end if
        write(unit_history,'(22(ES25.16E3,","),10(I0,","),3(ES25.16E3,","),ES25.16E3)') &
            t/3600,epoch+t,nominal,lower,upper,pos_max,vel_max,patches,history%stats%depth_limited_patches, &
            history%stats%split_counts,step_pos,step_vel,norm2(extra_states(1:3,2)-nominal(1:3)), &
            norm2(extra_states(4:6,2)-nominal(4:6))
        do j=1,size(points,2)
            write(unit_shape,'(ES25.16E3,",",I0,6(",",ES25.16E3))') t/3600,j,ads(:,j)-nominal
        end do
        do j=1,n_validation
            write(unit_probes,'(ES25.16E3,",",I0,16(",",ES25.16E3))') t/3600,validation_ids(j),s(j),u(j), &
                ads(:,validation_ids(j)),truth(:,j),norm2(ads(1:3,validation_ids(j))-truth(1:3,j)), &
                norm2(ads(4:6,validation_ids(j))-truth(4:6,j))
        end do
        if(allocated(history%stats%split_events)) then
            do j=1,size(history%stats%split_events)
                total_events=total_events+1
                associate(event=>history%stats%split_events(j))
                    write(unit_events,'(2(ES25.16E3,","),5(I0,","),41(ES25.16E3,","),ES25.16E3)') &
                        previous/3600,t/3600,total_events,event%bfs_iteration,event%parent_depth,event%component,event%direction, &
                        event%errors,event%tolerances,event%excess,event%direction_scores,event%parent_center,event%parent_width
                end associate
            end do
        end if
        if(step==0 .or. step==n_steps .or. history%stats%split_occurred.or.probability_enabled) call snapshot(.false.)
        if(history%stats%split_occurred) call snapshot(.true.)
        write(*,'(A,F8.3,A,I0,A,ES12.4,A,ES12.4)') 'C2 hour=',t/3600,' patches=',patches, &
            ' ADS pos=',step_pos,' vel=',step_vel
        previous=t
    end do
    close(unit_history); close(unit_shape); close(unit_probes); close(unit_events)
    if(probability_enabled) then
        call probability_close(probability_output,status)
        if(status/=0) call fatal('probability output close failed')
        close(probability_probe_unit); probability_probe_unit=0
    end if
    call hfem_ads_history_destroy(history); call clear_states(); after=active_da_count()
    if(after/=before) call fatal('DA handle leak')
    open(newunit=unit,file=trim(prefix)//'_report.json',status='new')
    write(unit,'(A)') '{"schema_version":"srp_modes_ads_v1","formal_result":false,"n_variables":8,'
    write(unit,'(A)') '"parameterization":"SRP_MODES","variable_roles":["x","y","z","vx","vy","vz","s","u"],'
    call string_field('opm_file',trim(opm_file)); call string_field('config_file',trim(cfg_file))
    call string_field('mode_csv_file',trim(config%srp_mode_table_file))
    call string_field('mode_metadata_file',trim(config%srp_mode_metadata_file))
    call string_field('mode_csv_sha256',csv_hash); call string_field('mode_metadata_sha256',trim(meta_hash))
    call string_field('source_revision',trim(revision))
    call number_field('epoch0_tdb_s',epoch); call number_field('duration_hours',duration)
    call number_field('save_hours',cadence); call number_field('da_order',real(options%da_order,DP))
    call number_field('max_depth',real(options%max_split_depth,DP)); call number_field('domain_sigma',options%domain_sigma)
    call number_field('position_tolerance_km',pos_tol); call number_field('velocity_tolerance_kms',vel_tol)
    call number_field('rel_tol',options%rel_tol); call number_field('abs_tol',options%abs_tol)
    call number_field('dt_min_s',options%dt_min); call number_field('dt_max_s',options%dt_max)
    call number_field('shape_count',real(size(points,2),DP)); call number_field('validation_count',real(n_validation,DP))
    call number_field('split_event_count',real(total_events,DP)); call number_field('depth_limited_patches',real(depth_count,DP))
    call number_field('da_handles_baseline',real(before,DP)); call number_field('da_handles_final',real(after,DP))
    call bool_field('validation_failed',any_failed); call bool_field('precision_unconfirmed',unconfirmed)
    if(probability_enabled) then
        write(unit,'(A)') '"probability_outputs":{'
        call number_field('statistics_sample_count',real(probability_count,DP))
        call number_field('visualization_sample_count',real(size(probability_output%visual_ids),DP))
        call number_field('validation_count',real(probability_probe_count,DP))
        call bool_field('validation_failed',probability_failed)
        call number_field('max_position_error_km',probability_max_pos)
        write(unit,'(A,ES25.16E3,A)') '"max_velocity_error_kms":',probability_max_vel,'},'
    end if
    write(unit,'(A,ES25.16E3,A,ES25.16E3,A)') '"ads_numerical_error":{"max_position_km":',max_pos, &
        ',"max_velocity_kms":',max_vel,'},'
    write(unit,'(A,ES25.16E3,A,ES25.16E3,A)') '"zero_mode_error":{"max_position_km":',zero_pos, &
        ',"max_velocity_kms":',zero_vel,',"control":"cannonball_real_nominal"},'
    write(unit,'(A,ES25.16E3,A,ES25.16E3,A)') '"zero_mode_real_control_error":{"max_position_km":',zero_real_pos, &
        ',"max_velocity_kms":',zero_real_vel,'},'
    write(unit,'(A,ES25.16E3,A,ES25.16E3,A)') '"conditional_srp_offset":{"final_earth_position_km":',offset_pos(1), &
        ',"final_moon_position_km":',offset_pos(2),'},'
    write(unit,'(A)') '"surrogate_model_error":{"status":"not_evaluated_P4"},"snapshot_prefixes":['
    do j=1,snapshot_count
        if(j>1) write(unit,'(A)',advance='no') ','
        write(unit,'(A)') quoted(trim(snapshots(j)))
    end do
    write(unit,'(A)') ']}'
    close(unit); close(unit_lock,status='delete'); lock_open=.false.
    if(any_failed.or.probability_failed.or.zero_pos>1e-6_DP.or.zero_vel>1e-9_DP) &
        call fatal('C2 numerical gate failed; raw diagnostics preserved')
    write(*,'(A)') 'PASS: C2 raw run; formal publication requires Python audit'
contains
    subroutine bind_metadata_hash()
        ! Use the platform digest utility, not a second Fortran SHA implementation.
        ! The redirection target is exclusively created by this runner before use.
        character(len=:),allocatable :: hash_file,command
        character(len=64) :: actual_hash
        integer :: hash_unit,io,command_status,exit_status,k
        hash_file=trim(prefix)//'_mode_meta_hash.tmp'
        open(newunit=hash_unit,file=hash_file,status='new',iostat=io)
        if(io/=0) call fatal('cannot reserve metadata hash temporary file')
        close(hash_unit)
        command='sha256sum --zero -- '//shell_quoted(trim(config%srp_mode_metadata_file))//' > '//shell_quoted(hash_file)
        call execute_command_line(command,exitstat=exit_status,cmdstat=command_status)
        if(command_status/=0) call fatal('cannot execute sha256sum for metadata identity')
        if(exit_status/=0) call fatal('sha256sum failed for metadata identity')
        open(newunit=hash_unit,file=hash_file,status='old',access='stream',form='unformatted',iostat=io)
        if(io/=0) call fatal('cannot read metadata digest')
        read(hash_unit,iostat=io) actual_hash
        close(hash_unit,status='delete')
        if(io/=0) call fatal('incomplete metadata digest')
        do k=1,64
            if(index('0123456789abcdef',actual_hash(k:k))==0) call fatal('invalid metadata digest')
        end do
        if(trim(meta_hash)/='UNVERIFIED' .and. meta_hash/=actual_hash) call fatal('metadata hash differs from expected identity')
        meta_hash=actual_hash
    end subroutine
    function shell_quoted(text) result(value)
        character(len=*),intent(in) :: text
        character(len=:),allocatable :: value
        integer :: k
        value=achar(39)
        do k=1,len(text)
            if(text(k:k)==achar(39)) then
                value=value//achar(39)//achar(92)//achar(39)//achar(39)
            else
                value=value//text(k:k)
            end if
        end do
        value=value//achar(39)
    end function
    subroutine clear_states()
        call clear_srp_scale_error(); call clear_srp_mode_coefficients_real()
        call clear_srp_scale_uncertainty(); call clear_srp_mode_uncertainty(); call clear_srp_mode_table()
    end subroutine
    subroutine fatal(text)
        character(len=*),intent(in) :: text
        write(*,'(A)') 'C2 failed: '//text
        call hfem_ads_history_destroy(history); call clear_states()
        error stop 2
    end subroutine
    subroutine snapshot(use_previous)
        logical,intent(in) :: use_previous
        real(DP) :: time
        time=history%current_time
        if(use_previous) time=history%previous_time
        write(tag,'(I7.7)') nint(time)
        filename=trim(prefix)//'_patch_s'//trim(tag)
        inquire(file=trim(filename)//'_meta.csv',exist=exists)
        if(exists) return
        call write_ads_patch_snapshot(history,trim(filename),use_previous,status,message,schema_version=2, &
            mode_csv_sha256=csv_hash,mode_metadata_sha256=trim(meta_hash))
        if(status/=0) call fatal(trim(message))
        snapshot_count=snapshot_count+1; snapshots(snapshot_count)=filename
    end subroutine
    function quoted(text) result(value)
        character(len=*),intent(in) :: text
        character(len=:),allocatable :: value
        integer :: k
        value='"'
        do k=1,len_trim(text)
            if(text(k:k)=='"' .or. text(k:k)==achar(92)) value=value//achar(92)
            if(iachar(text(k:k))<32) then
                value=value//' '
            else
                value=value//text(k:k)
            end if
        end do
        value=value//'"'
    end function
    subroutine string_field(key,value)
        character(len=*),intent(in) :: key,value
        write(unit,'(A)') quoted(key)//':'//quoted(value)//','
    end subroutine
    subroutine number_field(key,value)
        character(len=*),intent(in) :: key
        real(DP),intent(in) :: value
        write(unit,'(A,ES25.16E3,A)') quoted(key)//':',value,','
    end subroutine
    subroutine bool_field(key,value)
        character(len=*),intent(in) :: key
        logical,intent(in) :: value
        if(value) then
            write(unit,'(A)') quoted(key)//':true,'
        else
            write(unit,'(A)') quoted(key)//':false,'
        end if
    end subroutine
end program
