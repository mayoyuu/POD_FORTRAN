!> Independent Real nominal endpoint trajectories; never creates a DA domain.
program run_srp_modes_truth
    use pod_global, only: DP,MAX_STRING_LEN
    use pod_config, only: config,validate_config
    use pod_engine_module, only: pod_engine_init
    use pod_data_format_module, only: load_initial_opm
    use pod_srp_ads_history_module, only: validate_json_opm_fields
    use pod_srp_mode_module, only: load_srp_mode_table,clear_srp_mode_table
    use pod_force_model_module, only: set_propagation_epoch,set_srp_scale_error,clear_srp_scale_error, &
        set_srp_mode_coefficients_real,clear_srp_mode_coefficients_real
    use pod_integrator_module, only: adaptive_step_integrate,METHOD_RKF78
    use pod_dace_classes, only: active_da_count
    use pod_spacecraft_geometry, only: reset_srp_geometry_diagnostics,get_srp_geometry_diagnostics, &
        srp_geometry_diagnostics_type
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    character(len=MAX_STRING_LEN) :: opm='',cfg='',prefix='',arg,value,message
    character(len=16),parameter :: scenarios(5)=[character(len=16) :: &
        'Deta','proxy_earth','proxy_moon','boxwing_earth','boxwing_moon']
    real(DP) :: hours=360,cadence=1,rel_tol=1e-13_DP,abs_tol=1e-13_DP,epoch,initial(6),cov(6,6)
    real(DP) :: state(6),nd(6),previous,t
    real(DP),allocatable :: times(:),integrated(:,:)
    integer :: i,status,profile,step,intervals,n_steps,unit,report_unit,lock_unit,before,after
    logical :: exists
    type(srp_geometry_diagnostics_type) :: diagnostics
    i=1
    do while(i<=command_argument_count())
        call get_command_argument(i,arg,status=status)
        if(status/=0) call fatal('CLI argument too long')
        if(trim(arg)=='--help') then
            write(*,'(A)') 'Real-only SRP truth: -opm INPUT -cfg MODES_CONFIG -o PREFIX'
            write(*,'(A)') '--hours 360 --save-hours 1 --rel-tol 1e-13 --abs-tol 1e-13'
            stop
        end if
        if(i==command_argument_count()) call fatal('option requires a value')
        call get_command_argument(i+1,value,status=status)
        if(status/=0) call fatal('CLI value too long')
        select case(trim(arg))
        case('-opm'); opm=value
        case('-cfg'); cfg=value
        case('-o'); prefix=value
        case('--hours'); read(value,*,iostat=status) hours
        case('--save-hours'); read(value,*,iostat=status) cadence
        case('--rel-tol'); read(value,*,iostat=status) rel_tol
        case('--abs-tol'); read(value,*,iostat=status) abs_tol
        case default; call fatal('unknown option: '//trim(arg))
        end select
        if(status/=0) call fatal('invalid CLI number')
        i=i+2
    end do
    if(min(len_trim(opm),len_trim(cfg),len_trim(prefix))==0) call fatal('-opm, -cfg and -o required')
    if(.not.all(ieee_is_finite([hours,cadence,rel_tol,abs_tol]))) call fatal('numbers must be finite')
    if(min(hours,cadence,rel_tol,abs_tol)<=0 .or. hours>360 .or. hours/cadence>10000) call fatal('invalid limits')
    if(len_trim(prefix)>MAX_STRING_LEN-24) call fatal('output path too long')
    intervals=nint(hours/cadence)
    if(intervals<1 .or. abs(intervals*cadence-hours)>1e-10_DP) call fatal('whole save intervals required')
    call validate_json_opm_fields(trim(opm),status,message)
    if(status/=0) call fatal(trim(message))
    before=active_da_count()
    call pod_engine_init(trim(cfg))
    if(.not.validate_config() .or. .not.config%use_srp) call fatal('invalid active-SRP configuration')
    call load_initial_opm(trim(opm),epoch,initial,cov)
    inquire(file=trim(prefix)//'_report.json',exist=exists)
    if(exists) call fatal('output exists')
    open(newunit=lock_unit,file=trim(prefix)//'.lock',status='new',iostat=status)
    if(status/=0) call fatal('output exists or prefix unavailable')
    open(newunit=unit,file=trim(prefix)//'_states.csv',status='new',iostat=status)
    if(status/=0) call fatal('output CSV exists or unavailable')
    write(unit,'(A)') 'scenario,hour,epoch_et,x_km,y_km,z_km,vx_kms,vy_kms,vz_kms'
    call reset_srp_geometry_diagnostics()
    do profile=1,5
        call pod_engine_init(trim(cfg))
        if(config%max_propagation_steps<2) call fatal('invalid max_propagation_steps')
        select case(profile)
        case(1)
            config%srp_model='cannonball'
            call set_srp_scale_error(0.3_DP)
        case(2,3)
            config%srp_model='time_varying_modes'
            call load_srp_mode_table(trim(config%srp_mode_table_file),trim(config%srp_mode_metadata_file),epoch,status,message)
            if(status/=0) call fatal(trim(message))
            if(profile==2) call set_srp_mode_coefficients_real(1.0_DP,0.5_DP)
            if(profile==3) call set_srp_mode_coefficients_real(0.0_DP,1.0_DP)
        case(4,5)
            config%srp_model='box_wing'
            config%srp_attitude_mode='earth'
            if(profile==5) config%srp_attitude_mode='moon'
            if(.not.validate_config()) call fatal('invalid physical Box-Wing geometry')
        end select
        call set_propagation_epoch(epoch)
        state=initial; previous=0
        do step=0,intervals
            t=real(step,DP)*cadence*3600
            if(step>0) then
                nd(1:3)=state(1:3)/config%LU; nd(4:6)=state(4:6)/config%VU
                call adaptive_step_integrate(nd,previous/config%TU,t/config%TU,METHOD_RKF78,times,integrated, &
                    n_steps,config%max_propagation_steps,rel_tol,abs_tol,1e-6_DP,3600.0_DP)
                if(.not.allocated(times) .or. .not.allocated(integrated) .or. n_steps<1) call fatal('did not reach end time')
                if(n_steps>size(times) .or. n_steps>size(integrated,1)) call fatal('did not reach end time')
                if(.not.all(ieee_is_finite(integrated(n_steps,:))) .or. &
                    .not.ieee_is_finite(times(n_steps))) call fatal('did not reach finite end state')
                if(abs(times(n_steps)-t/config%TU)>1000*epsilon(1.0_DP)*max(1.0_DP,abs(t/config%TU))) &
                    call fatal('did not reach end time')
                state(1:3)=integrated(n_steps,1:3)*config%LU
                state(4:6)=integrated(n_steps,4:6)*config%VU
            end if
            if(.not.all(ieee_is_finite(state))) call fatal('nonfinite physical state')
            write(unit,'(A,8(",",ES25.16E3))') trim(scenarios(profile)),t/3600,epoch+t,state
            previous=t
        end do
        call clear_srp_scale_error(); call clear_srp_mode_coefficients_real(); call clear_srp_mode_table()
        write(*,'(A,A,A,F8.3)') 'Real truth ',trim(scenarios(profile)),' completed hours=',hours
    end do
    close(unit)
    call get_srp_geometry_diagnostics(diagnostics)
    after=active_da_count()
    if(before/=after) call fatal('unexpected DA ownership')
    open(newunit=report_unit,file=trim(prefix)//'_report.json',status='new')
    write(report_unit,'(A)') '{"schema_version":"srp_modes_truth_v1","formal_result":false,'
    write(report_unit,'(A)') '"state_source":"independent_real_nominal",'
    write(report_unit,'(A,ES25.16E3,A)') '"duration_hours":',hours,','
    write(report_unit,'(A,ES25.16E3,A)') '"save_hours":',cadence,','
    write(report_unit,'(A,ES25.16E3,A)') '"epoch0_tdb_s":',epoch,','
    write(report_unit,'(A,ES25.16E3,A)') '"rel_tol":',rel_tol,','
    write(report_unit,'(A,ES25.16E3,A)') '"abs_tol":',abs_tol,','
    write(report_unit,'(A,I0,A)') '"row_count":',5*(intervals+1),','
    write(report_unit,'(A,I0,A)') '"geometry_fallback_count":',diagnostics%real_fallback_count,','
    write(report_unit,'(A,I0,A)') '"da_handles_baseline":',before,','
    write(report_unit,'(A,I0,A)') '"da_handles_final":',after,'}'
    close(report_unit); close(lock_unit,status='delete')
    write(*,'(A)') 'PASS: independent Real-only truth; formal publication requires Python audit'
contains
    subroutine fatal(text)
        character(len=*),intent(in) :: text
        call clear_srp_scale_error(); call clear_srp_mode_coefficients_real(); call clear_srp_mode_table()
        write(*,'(A)') 'Real truth failed: '//text
        error stop 2
    end subroutine
end program
