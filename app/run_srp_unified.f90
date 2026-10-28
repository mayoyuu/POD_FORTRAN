! Shared, file-driven uncertainty propagation. All physical states are km, km/s.
program run_srp_unified
    use, intrinsic :: iso_fortran_env, only: int64, error_unit, output_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use pod_global, only: DP
    use pod_config, only: config
    use pod_engine_module, only: pod_engine_init
    use pod_spice, only: set_spice_kernel_path, str2et
    use pod_data_format_module, only: load_initial_opm
    use pod_srp_mode_module, only: load_srp_mode_table
    use pod_srp_ads_history_module, only: advance_real_probes, advance_real_mode_probes, &
        validate_json_opm_fields
    use pod_uq_hfem_ads_module, only: hfem_ads_options_type, hfem_ads_history_type, &
        hfem_ads_history_init, hfem_ads_history_advance, hfem_ads_history_evaluate, &
        hfem_ads_history_destroy
    use pod_uq_ads_coordinates_module, only: ads_physical_to_unit, ads_physical_parameters_type, &
        ADS_COORD_WHITENED, ADS_PARAM_NONE, ADS_PARAM_SRP_ETA, ADS_PARAM_SRP_MODES
    implicit none

    integer, parameter :: PATH_LEN=4096
    character(len=PATH_LEN) :: opm_file='', config_file='', points_file='', outdir='', kernel_dir=''
    character(len=16) :: method='', error_kind=''
    character(len=128) :: end_utc=''
    real(DP) :: end_et
    integer :: requested_saved_points=2000, requested_final_only=0, span_unit
    real(DP) :: hours=-1, save_hours=-1, eta_sigma=-1, state_scale=-1, domain_sigma=-1
    real(DP) :: pos_tol=-1, vel_tol=-1
    integer :: order=-1, max_order=-1, adaptive=-1, requested_probes=-1
    logical :: seen(21)=.false., is_map, is_modes, accuracy_valid=.true., precision_confirmed=.true.
    type(hfem_ads_options_type) :: options
    type(hfem_ads_history_type) :: history
    type(ads_physical_parameters_type) :: parameters
    real(DP) :: epoch0, nominal(6), initial_cov(6,6), reference(6,1), mean(6), covariance(6,6)
    real(DP), allocatable :: states(:,:), initial(:,:), eta(:), s(:), u(:), wm(:), wc(:)
    real(DP), allocatable :: unit_points(:,:), one_unit(:), probe_states(:,:), schedule(:)
    character(len=128), allocatable :: sample_ids(:)
    integer, allocatable :: saved_ids(:), probe_ids(:)
    logical, allocatable :: found(:), save_at(:)
    integer :: n, i, j, k, status, nm, np, step, replay, nsteps, depth_seen, patches, depth_limited
    integer :: moments_unit, points_unit, splits_unit, timing_unit, probe_unit, ios
    character(len=4096) :: message
    real(DP) :: start_wall, mark, current, previous, validation_pos, validation_vel
    real(DP) :: propagation_seconds=0, evaluation_seconds=0, statistics_seconds=0
    real(DP) :: validation_seconds=0, io_seconds=0, compute_seconds, total_seconds

    start_wall=wall_seconds()
    call parse_arguments()
    is_map=method=='ADS' .or. method=='DAMC'
    mark=wall_seconds()
    call read_points()
    io_seconds=io_seconds+wall_seconds()-mark
    if(len_trim(kernel_dir)>0) call set_spice_kernel_path(trim(kernel_dir))
    call pod_engine_init(trim(config_file))
    call validate_json_opm_fields(trim(opm_file),status,message)
    call require_ok('OPM preflight')
    call load_initial_opm(trim(opm_file),epoch0,nominal,initial_cov)
    if(.not.all(ieee_is_finite([epoch0,nominal]))) call fail('nonfinite OPM state/epoch')
    if(len_trim(end_utc)>0) then
        call str2et(trim(end_utc),end_et)
        if(.not.ieee_is_finite(end_et) .or. end_et<=epoch0) call fail('end UTC must be after OPM epoch')
        hours=(end_et-epoch0)/3600.0_DP
    else
        end_et=epoch0+3600.0_DP*hours
    end if
    call open_output('time_span.csv',span_unit)
    write(span_unit,'(a)') 'start_et,end_et,duration_hours'
    write(span_unit,'(3(es24.16e3,:,","))') epoch0,end_et,hours
    close(span_unit)
    initial_cov=initial_cov*state_scale**2
    is_modes=trim(config%srp_model)=='time_varying_modes'
    if(is_modes) then
        call load_srp_mode_table(trim(config%srp_mode_table_file), &
            trim(config%srp_mode_metadata_file),epoch0,status,message)
        call require_ok('mode table')
        if(error_kind=='eta') call fail('eta cannot be combined with time_varying_modes')
        if(any(abs(eta)>0.0_DP)) call fail('mode trajectories require eta=0')
    else if(error_kind=='attitude') then
        call fail('attitude requires time_varying_modes config')
    end if
    options%coordinate_mode=ADS_COORD_WHITENED
    options%parameterization=ADS_PARAM_NONE
    if(error_kind=='eta') then
        options%parameterization=ADS_PARAM_SRP_ETA
        options%srp_sigma=eta_sigma
    end if
    if(is_modes) options%parameterization=ADS_PARAM_SRP_MODES
    options%domain_sigma=domain_sigma
    options%da_order=order
    options%error_tolerance(1:3)=pos_tol
    options%error_tolerance(4:6)=vel_tol
    options%rel_tol=config%rkf78_rel_tol
    options%abs_tol=config%rkf78_abs_tol
    options%dt_min=config%rkf78_min_step
    options%dt_max=config%rkf78_max_step
    options%max_steps=config%max_propagation_steps
    if(method=='DAMC') options%max_split_depth=0
    ! Engine/bootstrap work is included in total_seconds, separate from propagation.
    mark=wall_seconds()
    reference(:,1)=nominal
    depth_seen=0
    if(is_map) then
        call init_history()
        allocate(unit_points(history%coordinate_map%n_variables,n), found(n))
        do i=1,n
            parameters%eta=eta(i); parameters%s=s(i); parameters%u=u(i)
            call ads_physical_to_unit(history%coordinate_map,initial(:,i)-nominal,parameters,one_unit,status)
            if(status/=0) call fail('sample cannot be mapped to ADS coordinates: '//trim(sample_ids(i)))
            if(any(abs(one_unit)>1.0_DP+1e-12_DP)) &
                call fail('sample outside ADS domain: '//trim(sample_ids(i)))
            unit_points(:,i)=max(-1.0_DP,min(1.0_DP,one_unit))
        end do
    end if
    nm=n
    if(requested_saved_points>0) nm=min(requested_saved_points,n)
    if(method=='UT') nm=n
    allocate(saved_ids(nm))
    call evenly_spaced(n,saved_ids)
    np=0
    if(is_map) np=min(requested_probes,n)
    allocate(probe_ids(np),probe_states(6,np))
    call evenly_spaced(n,probe_ids)
    if(np>0) probe_states=initial(:,probe_ids)
    call build_schedule()
    propagation_seconds=propagation_seconds+wall_seconds()-mark
    mark=wall_seconds()
    call open_output('hourly_moments.csv',moments_unit)
    call open_output('points.csv',points_unit)
    call open_output('splits.csv',splits_unit)
    call open_output('timing.csv',timing_unit)
    call open_output('validation_ids.csv',probe_unit)
    call write_headers()
    write(probe_unit,'(a)') 'probe_index,sample_id'
    do i=1,np
        write(probe_unit,'(i0,a,a)') i,',',trim(sample_ids(probe_ids(i)))
    end do
    close(probe_unit)
    io_seconds=io_seconds+wall_seconds()-mark
    validation_pos=0; validation_vel=0; current=0
    patches=0; depth_limited=0
    if(is_map) patches=history%stats%n_patches
    if(requested_final_only==0 .or. nsteps==0) then
        call save_output()
    else
        call write_progress()
    end if

    do step=1,nsteps
        previous=current
        current=schedule(step)
        mark=wall_seconds()
        if(is_map) then
            call hfem_ads_history_advance(history,current,status,message)
            call require_ok('history advance')
            depth_seen=max(depth_seen,history%stats%depth_limited_patches)
        else
            call advance_particles(states,eta,s,u,previous,current)
        end if
        call advance_reference(previous,current)
        propagation_seconds=propagation_seconds+wall_seconds()-mark
        if(is_map) then
            mark=wall_seconds()
            if(np>0) call advance_particles(probe_states,eta(probe_ids),s(probe_ids),u(probe_ids), &
                previous,current)
            validation_seconds=validation_seconds+wall_seconds()-mark
            call evaluate_history()
            call validate_history()
            if(method=='DAMC' .and. adaptive==1) then
                do while((validation_pos>pos_tol .or. validation_vel>vel_tol) .and. order<max_order)
                    ! A new algebra order cannot recover terms discarded by an old map.
                    ! Destroy every old DA handle BEFORE reinitializing DACE, then replay.
                    mark=wall_seconds()
                    call hfem_ads_history_destroy(history)
                    order=order+1; options%da_order=order; depth_seen=0
                    call init_history()
                    do replay=1,step
                        call hfem_ads_history_advance(history,schedule(replay),status,message)
                        call require_ok('DAMC initial-epoch replay')
                        depth_seen=max(depth_seen,history%stats%depth_limited_patches)
                    end do
                    propagation_seconds=propagation_seconds+wall_seconds()-mark
                    call evaluate_history()
                    call validate_history()
                end do
            end if
            if(np==0 .or. validation_pos>pos_tol .or. validation_vel>vel_tol) accuracy_valid=.false.
            patches=history%stats%n_patches
            depth_limited=history%stats%depth_limited_patches
            precision_confirmed=accuracy_valid .and. depth_seen==0
            if(method=='ADS') then
                mark=wall_seconds()
                call write_splits(previous/3600.0_DP,current/3600.0_DP)
                io_seconds=io_seconds+wall_seconds()-mark
            end if
        end if
        if(save_at(step) .and. (requested_final_only==0 .or. step==nsteps)) then
            call save_output()
        else if(requested_final_only==1) then
            ! Progress is independent of particle snapshots in terminal-only runs.
            call write_progress()
        end if
    end do
    if(is_map) call hfem_ads_history_destroy(history)
    mark=wall_seconds()
    close(moments_unit); close(points_unit); close(splits_unit)
    write(timing_unit,'(a)') 'propagation_seconds,evaluation_seconds,statistics_seconds,'// &
        'validation_seconds,io_seconds,compute_excluding_io_seconds,total_seconds'
    flush(timing_unit)
    io_seconds=io_seconds+wall_seconds()-mark
    compute_seconds=propagation_seconds+evaluation_seconds+statistics_seconds+validation_seconds
    total_seconds=wall_seconds()-start_wall
    write(timing_unit,'(*(es24.16e3,:,","))') propagation_seconds,evaluation_seconds,statistics_seconds, &
        validation_seconds,io_seconds,compute_seconds,total_seconds
    close(timing_unit)

contains

    real(DP) function wall_seconds() result(value)
        integer(int64) :: ticks,rate
        call system_clock(ticks,rate)
        value=real(ticks,DP)/real(rate,DP)
    end function

    subroutine fail(text)
        character(len=*), intent(in) :: text
        write(error_unit,'(a)') 'run_srp_unified: '//trim(text)
        stop 2
    end subroutine

    subroutine require_ok(context)
        character(len=*),intent(in) :: context
        if(status/=0) call fail(trim(context)//': '//trim(message))
    end subroutine

    subroutine help()
        write(*,'(a)') 'run_srp_unified --method ADS|UT|MC|DAMC --opm PATH --config PATH --points PATH --outdir PATH'
        write(*,'(a)') '  --hours REAL --save-hours REAL --error none|eta|attitude --eta-sigma REAL'
        write(*,'(a)') '  --state-scale REAL --domain-sigma REAL --order INT --max-order INT --adaptive 0|1'
        write(*,'(a)') '  --probes INT --pos-tol-km REAL --vel-tol-kms REAL [--kernel-dir PATH]'
        write(*,'(a)') '  [--saved-points INT (0=all; default 2000)] [--end-utc ISO-UTC] [--final-only 0|1]'
        write(*,'(a)') 'Base options except --kernel-dir are required; --end-utc overrides --hours.'
    end subroutine

    subroutine parse_arguments()
        integer :: arg, nargs, idx, rc
        character(len=PATH_LEN) :: key,value
        nargs=command_argument_count()
        if(nargs==1) then
            call get_command_argument(1,key)
            if(trim(key)=='--help') then
                call help()
                stop
            end if
        end if
        arg=1
        do while(arg<=nargs)
            call get_command_argument(arg,key)
            select case(trim(key))
            case('--method'); idx=1
            case('--opm'); idx=2
            case('--config'); idx=3
            case('--points'); idx=4
            case('--outdir'); idx=5
            case('--hours'); idx=6
            case('--save-hours'); idx=7
            case('--error'); idx=8
            case('--eta-sigma'); idx=9
            case('--state-scale'); idx=10
            case('--domain-sigma'); idx=11
            case('--order'); idx=12
            case('--max-order'); idx=13
            case('--adaptive'); idx=14
            case('--probes'); idx=15
            case('--pos-tol-km'); idx=16
            case('--vel-tol-kms'); idx=17
            case('--kernel-dir'); idx=18
            case('--saved-points'); idx=19
            case('--end-utc'); idx=20
            case('--final-only'); idx=21
            case default
                call fail('unknown option '//trim(key))
            end select
            if(seen(idx)) call fail('duplicate option '//trim(key))
            if(arg==nargs) call fail('missing value for '//trim(key))
            call get_command_argument(arg+1,value)
            if(len_trim(value)==0 .or. index(value,'--')==1) call fail('missing value for '//trim(key))
            seen(idx)=.true.; rc=0
            select case(idx)
            case(1); method=trim(value)
            case(2); opm_file=trim(value)
            case(3); config_file=trim(value)
            case(4); points_file=trim(value)
            case(5); outdir=trim(value)
            case(6); read(value,*,iostat=rc) hours
            case(7); read(value,*,iostat=rc) save_hours
            case(8); error_kind=trim(value)
            case(9); read(value,*,iostat=rc) eta_sigma
            case(10); read(value,*,iostat=rc) state_scale
            case(11); read(value,*,iostat=rc) domain_sigma
            case(12); read(value,*,iostat=rc) order
            case(13); read(value,*,iostat=rc) max_order
            case(14); read(value,*,iostat=rc) adaptive
            case(15); read(value,*,iostat=rc) requested_probes
            case(16); read(value,*,iostat=rc) pos_tol
            case(17); read(value,*,iostat=rc) vel_tol
            case(18); kernel_dir=trim(value)
            case(19); read(value,*,iostat=rc) requested_saved_points
            case(20); end_utc=trim(value)
            case(21); read(value,*,iostat=rc) requested_final_only
            end select
            if(rc/=0) call fail('invalid value for '//trim(key))
            arg=arg+2
        end do
        if(.not.all(seen(1:17))) call fail('missing required options; use --help')
        if(method/='ADS' .and. method/='UT' .and. method/='MC' .and. method/='DAMC') call fail('invalid method')
        if(error_kind/='none' .and. error_kind/='eta' .and. error_kind/='attitude') call fail('invalid error')
        if(.not.all(ieee_is_finite([hours,save_hours,eta_sigma,state_scale,domain_sigma,pos_tol,vel_tol]))) &
            call fail('numeric options must be finite')
        if(hours<0 .or. save_hours<=0 .or. eta_sigma<0 .or. state_scale<=0 .or. domain_sigma<=0 .or. &
            pos_tol<=0 .or. vel_tol<=0) call fail('invalid numeric option range')
        if(order<1 .or. max_order<order .or. requested_probes<0 .or. requested_saved_points<0 .or. &
            (adaptive/=0 .and. adaptive/=1)) call fail('invalid integer option range')
        if(requested_final_only/=0 .and. requested_final_only/=1) call fail('invalid final-only flag')
        if(error_kind=='eta' .and. eta_sigma<=0) call fail('eta error requires positive eta-sigma')
    end subroutine

    subroutine read_points()
        integer :: input_unit, rc, row, col, last, comma
        character(len=8192) :: line
        character(len=256) :: fields(12)
        real(DP) :: numeric(11)
        open(newunit=input_unit,file=trim(points_file),status='old',action='read',iostat=rc)
        if(rc/=0) call fail('cannot open points input')
        read(input_unit,'(a)',iostat=rc) line
        if(rc/=0 .or. trim(line)/='sample_id,x,y,z,vx,vy,vz,eta,s,u,wm,wc') call fail('incorrect points header')
        n=0
        do
            read(input_unit,'(a)',iostat=rc) line
            if(rc<0) exit
            if(rc/=0) call fail('cannot count points input')
            if(len_trim(line)==0) call fail('blank points record')
            n=n+1
        end do
        if(n<1 .or. (method/='UT' .and. n<2)) call fail('insufficient input samples')
        allocate(initial(6,n),states(6,n),eta(n),s(n),u(n),wm(n),wc(n),sample_ids(n))
        rewind(input_unit)
        read(input_unit,'(a)') line
        do row=1,n
            read(input_unit,'(a)',iostat=rc) line
            if(rc/=0) call fail('cannot read points record')
            last=1
            do col=1,11
                comma=index(line(last:),',')
                if(comma==0) call fail('points record must have 12 fields')
                if(comma==1 .or. comma-1>len(fields(col))) call fail('invalid points field')
                fields(col)=line(last:last+comma-2)
                last=last+comma
            end do
            if(index(line(last:),',')/=0 .or. len_trim(line(last:))==0) call fail('invalid final points field')
            fields(12)=trim(line(last:))
            if(len_trim(fields(1))>len(sample_ids(row))) call fail('sample id exceeds 128 characters')
            sample_ids(row)=trim(adjustl(fields(1)))
            do col=1,11
                read(fields(col+1),*,iostat=rc) numeric(col)
                if(rc/=0) call fail('invalid numeric sample field')
            end do
            if(.not.all(ieee_is_finite(numeric))) call fail('sample values must be finite')
            initial(:,row)=numeric(1:6)
            eta(row)=numeric(7); s(row)=numeric(8); u(row)=numeric(9)
            wm(row)=numeric(10); wc(row)=numeric(11)
        end do
        close(input_unit)
        states=initial
        if(abs(sum(wm)-1.0_DP)>1e-8_DP*max(1.0_DP,sum(abs(wm)))) call fail('mean weights must sum to one')
        if(method/='UT') then
            if(maxval(abs(wm-1.0_DP/real(n,DP)))>1e-12_DP .or. &
                maxval(abs(wc-1.0_DP/real(n-1,DP)))>1e-12_DP) call fail('MC/map weights must be 1/N and 1/(N-1)')
        end if
        if(error_kind=='none' .and. any(abs(eta)>0.0_DP)) call fail('none error requires eta=0')
    end subroutine

    subroutine evenly_spaced(count,ids)
        integer,intent(in) :: count
        integer,intent(out) :: ids(:)
        integer :: index_id, m
        m=size(ids)
        if(m==1) then
            ids(1)=1
        else
            do index_id=1,m
                ids(index_id)=1+nint(real(index_id-1,DP)*real(count-1,DP)/real(m-1,DP))
            end do
        end if
    end subroutine

    subroutine build_schedule()
        real(DP) :: next_hour,next_save,finish,next_time,tol
        integer :: capacity
        capacity=ceiling(hours)+ceiling(hours/save_hours)+2
        allocate(schedule(capacity),save_at(capacity))
        next_hour=3600.0_DP; next_save=3600.0_DP*save_hours; finish=3600.0_DP*hours
        tol=100*epsilon(1.0_DP)*max(1.0_DP,finish)
        nsteps=0
        do while(min(next_hour,next_save,finish)>0)
            next_time=min(next_hour,next_save,finish)
            if(abs(next_time-finish)<=tol) next_time=finish
            nsteps=nsteps+1; schedule(nsteps)=next_time
            save_at(nsteps)=abs(next_time-next_save)<=tol .or. abs(next_time-finish)<=tol
            if(abs(next_time-next_hour)<=tol) next_hour=next_hour+3600.0_DP
            if(abs(next_time-next_save)<=tol) next_save=next_save+3600.0_DP*save_hours
            if(abs(next_time-finish)<=tol) exit
        end do
    end subroutine

    subroutine init_history()
        call hfem_ads_history_init(history,nominal,initial_cov,epoch0,options,status,message)
        call require_ok('history initialization')
    end subroutine

    subroutine advance_particles(values,eta_values,s_values,u_values,t0,t1)
        real(DP),intent(inout) :: values(:,:)
        real(DP),intent(in) :: eta_values(:),s_values(:),u_values(:),t0,t1
        if(is_modes) then
            call advance_real_mode_probes(values,s_values,u_values,epoch0,t0,t1,options,status,message)
        else
            call advance_real_probes(values,eta_values,epoch0,t0,t1,options,status,message)
        end if
        call require_ok('Real propagation')
    end subroutine

    subroutine advance_reference(t0,t1)
        real(DP),intent(in) :: t0,t1
        logical :: saved_srp
        character(len=len(config%srp_model)) :: saved_model
        saved_srp=config%use_srp; saved_model=config%srp_model
        ! Serial call: no particle propagation can see this temporary model.
        config%use_srp=.true.; config%srp_model='cannonball'
        call advance_real_probes(reference,[0.0_DP],epoch0,t0,t1,options,status,message)
        config%use_srp=saved_srp; config%srp_model=saved_model
        call require_ok('cannonball reference propagation')
    end subroutine

    subroutine evaluate_history()
        real(DP) :: begin
        begin=wall_seconds()
        call hfem_ads_history_evaluate(history,unit_points,states,found,status)
        if(status/=0 .or. .not.all(found)) call fail('ADS evaluation missed input samples')
        if(.not.all(ieee_is_finite(states))) call fail('ADS evaluation produced nonfinite values')
        evaluation_seconds=evaluation_seconds+wall_seconds()-begin
    end subroutine

    subroutine validate_history()
        real(DP) :: begin, difference(6)
        integer :: probe
        begin=wall_seconds()
        validation_pos=0; validation_vel=0
        do probe=1,np
            difference=states(:,probe_ids(probe))-probe_states(:,probe)
            validation_pos=max(validation_pos,norm2(difference(1:3)))
            validation_vel=max(validation_vel,norm2(difference(4:6)))
        end do
        validation_seconds=validation_seconds+wall_seconds()-begin
    end subroutine

    subroutine open_output(name,output_unit)
        character(len=*),intent(in) :: name
        integer,intent(out) :: output_unit
        open(newunit=output_unit,file=trim(outdir)//'/'//name,status='replace',action='write',iostat=ios)
        if(ios/=0) call fail('cannot open output '//trim(name)//'; outdir must already exist')
    end subroutine

    subroutine write_headers()
        integer :: a,b
        write(moments_unit,'(a)',advance='no') 'hour'
        do a=1,6
            write(moments_unit,'(a,i0)',advance='no') ',mean',a
        end do
        do a=1,6
            do b=1,6
                write(moments_unit,'(a,i0,i0)',advance='no') ',cov',a,b
            end do
        end do
        do a=1,6
            write(moments_unit,'(a,i0)',advance='no') ',reference',a
        end do
        write(moments_unit,'(a)') ',patches,depth_limited,validation_pos_km,validation_vel_kms,'// &
            'accuracy_valid,precision_confirmed,order,propagation_seconds,evaluation_seconds,'// &
            'statistics_seconds,validation_seconds,io_seconds'
        write(points_unit,'(a)') 'hour,sample_id,x,y,z,vx,vy,vz,eta,s,u,wm,wc'
        write(splits_unit,'(a)',advance='no') 'interval_start_hour,interval_end_hour,event,'// &
            'direction,parent_depth,component'
        do a=1,6
            write(splits_unit,'(a,i0)',advance='no') ',error',a
        end do
        do a=1,6
            write(splits_unit,'(a,i0)',advance='no') ',tolerance',a
        end do
        do a=1,6
            write(splits_unit,'(a,i0)',advance='no') ',excess',a
        end do
        do a=1,8
            write(splits_unit,'(a,i0)',advance='no') ',score',a
        end do
        do a=1,8
            write(splits_unit,'(a,i0)',advance='no') ',center',a
        end do
        do a=1,8
            write(splits_unit,'(a,i0)',advance='no') ',width',a
        end do
        write(splits_unit,'(a)') ''
    end subroutine

    subroutine write_splits(t0,t1)
        real(DP),intent(in) :: t0,t1
        integer :: event_id
        if(.not.allocated(history%stats%split_events)) return
        do event_id=1,size(history%stats%split_events)
            associate(e=>history%stats%split_events(event_id))
                write(splits_unit,'(es24.16e3,",",es24.16e3,4(",",i0),42(",",es24.16e3))') &
                    t0,t1,e%bfs_iteration,e%direction,e%parent_depth,e%component, &
                    e%errors,e%tolerances,e%excess,e%direction_scores,e%parent_center,e%parent_width
            end associate
        end do
        flush(splits_unit)
    end subroutine

    subroutine write_progress()
        write(*,'(a,a,a,f10.4,a,i0,a,i0,a,i0)') 'unified ',trim(method),' hour=',current/3600.0_DP, &
            ' order=',order,' patches=',patches,' accuracy_valid=',merge(1,0,accuracy_valid)
        flush(output_unit)
    end subroutine

    subroutine save_output()
        integer :: member,a,b
        real(DP) :: begin,deviation(6),mean_correction(6),term(6),adjusted(6),next_mean(6),anchor(6)
        begin=wall_seconds()
        ! Compensated summation protects signed UT weights from unnecessary cancellation.
        anchor=states(:,1)
        mean=0; mean_correction=0
        do member=1,n
            term=wm(member)*(states(:,member)-anchor)
            adjusted=term-mean_correction; next_mean=mean+adjusted
            mean_correction=(next_mean-mean)-adjusted; mean=next_mean
        end do
        mean=anchor+mean
        covariance=0
        do member=1,n
            deviation=states(:,member)-mean
            do a=1,6
                do b=1,6
                    covariance(a,b)=covariance(a,b)+wc(member)*deviation(a)*deviation(b)
                end do
            end do
        end do
        statistics_seconds=statistics_seconds+wall_seconds()-begin
        begin=wall_seconds()
        write(moments_unit,'(es24.16e3,48(",",es24.16e3),2(",",i0),'// &
            '2(",",es24.16e3),3(",",i0),5(",",es24.16e3))') &
            current/3600.0_DP,mean,reshape(transpose(covariance),[36]),reference(:,1), &
            patches,depth_limited,validation_pos,validation_vel,merge(1,0,accuracy_valid), &
            merge(1,0,precision_confirmed),order,propagation_seconds,evaluation_seconds, &
            statistics_seconds,validation_seconds,io_seconds
        do member=1,nm
            a=saved_ids(member)
            write(points_unit,'(es24.16e3,a,a,11(",",es24.16e3))') current/3600.0_DP,',', &
                trim(sample_ids(a)),states(:,a),eta(a),s(a),u(a),wm(a),wc(a)
        end do
        flush(moments_unit); flush(points_unit)
        call write_progress()
        io_seconds=io_seconds+wall_seconds()-begin
    end subroutine
end program run_srp_unified
