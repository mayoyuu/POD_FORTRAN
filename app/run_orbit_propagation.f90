! File-driven counterpart of test_orbit_output: hourly RKF78 ephemeris and final OPM.
program run_orbit_propagation
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use pod_global, only: DP
    use pod_config, only: config, validate_config
    use pod_engine_module, only: pod_engine_init
    use pod_spice, only: str2et, et2utc, set_spice_kernel_path
    use pod_json_object_module, only: json_member, json_number, json_string
    use pod_uq_transform_ut_module, only: generate_sigma_points, reconstruct_ut_moments
    use pod_uq_hfem_ads_module, only: hfem_ads_options_type
    use pod_srp_ads_history_module, only: advance_real_probes
    use pod_srp_mode_module, only: load_srp_mode_table
    implicit none
    character(len=4096) :: opm='', cfg='config/config.txt', prefix='', kernel_dir='', end_utc=''
    character(len=4096) :: arg, value, message
    character(len=:), allocatable :: source, raw, epoch_text, frame
    character(len=5), parameter :: keys(6)=[character(len=5) :: 'X','Y','Z','X_DOT','Y_DOT','Z_DOT']
    character(len=24), parameter :: metadata(*)=[character(len=24) :: &
        'ASL_CAT_ID','OBJECT_NAME','OBJECT_ID','DRAG_AREA','DRAG_CD','MASS', &
        'SOLAR_RAD_AREA','SOLAR_RAD_COEFF']
    character(len=32) :: utc
    real(DP) :: hours=-1, epoch0, end_et, duration, current, next_time
    real(DP) :: initial(6), covariance(6,6), trial(6,6), mean(6)
    real(DP), allocatable :: points(:,:), initial_points(:,:), weights(:), eta(:)
    real(DP), allocatable :: ut_mean(:), ut_cov(:,:), cross_cov(:,:)
    type(hfem_ads_options_type) :: options
    integer :: i,j,u,io,nbytes,ncov,step,idx
    logical :: has_cov, ok, exists, seen(7)=.false.
    interface
        subroutine dpotrf(uplo,n,a,lda,info)
            import DP
            character(len=1),intent(in) :: uplo
            integer,intent(in) :: n,lda
            real(DP),intent(inout) :: a(lda,*)
            integer,intent(out) :: info
        end subroutine
    end interface

    call arguments()
    ! Fail before engine initialization (load_config otherwise creates missing files).
    inquire(file=trim(cfg),exist=exists)
    if(.not.exists) call fatal('configuration file does not exist: '//trim(cfg))
    open(newunit=u,file=trim(opm),access='stream',form='unformatted',status='old',iostat=io)
    if(io/=0) call fatal('cannot open OPM: '//trim(opm))
    inquire(unit=u,size=nbytes)
    if(nbytes<2) call fatal('empty OPM')
    allocate(character(len=nbytes) :: source)
    read(u,iostat=io) source
    close(u)
    if(io/=0) call fatal('cannot read OPM')
    call json_string(source,'EPOCH',epoch_text,ok)
    if(.not.ok) call fatal('OPM requires a valid top-level EPOCH in JSON format')
    do i=1,6
        call json_number(source,trim(keys(i)),initial(i),ok)
        if(.not.ok) call fatal('OPM requires finite state field '//trim(keys(i)))
    end do
    call check_frame('REF_FRAME','J2000','GCRF')
    call check_frame('COV_REF_FRAME','J2000','GCRF')
    call check_frame('CENTER_NAME','EARTH','Earth')
    call check_frame('TIME_SYSTEM','UTC','UTC')
    ncov=0; covariance=0
    do i=1,6
        do j=1,i
            call json_member(source,cov_key(i,j),raw,ok)
            if(.not.ok) cycle
            ncov=ncov+1
            call json_number(source,cov_key(i,j),covariance(i,j),ok)
            if(.not.ok) call fatal('invalid covariance field '//cov_key(i,j))
            covariance(j,i)=covariance(i,j)
        end do
    end do
    if(ncov/=0 .and. ncov/=21) call fatal('incomplete covariance: all 21 lower-triangle fields are required')
    has_cov=ncov==21
    if(has_cov .and. any(covariance/=0)) then
        trial=covariance
        call dpotrf('L',6,trial,6,io)
        if(io/=0) call fatal('initial covariance must be positive definite (or identically zero)')
        call generate_sigma_points(initial,covariance,0.5_DP,points,weights)
    else
        allocate(points(6,1),weights(1))
        points(:,1)=initial; weights=1
    end if
    initial_points=points
    allocate(eta(size(weights))); eta=0
    if(len_trim(kernel_dir)>0) then
        if(kernel_dir(len_trim(kernel_dir):len_trim(kernel_dir))/='/') kernel_dir=trim(kernel_dir)//'/'
        call set_spice_kernel_path(trim(kernel_dir))
    end if
    call pod_engine_init(trim(cfg))
    if(.not.validate_config()) call fatal('invalid dynamics configuration')
    if(trim(config%central_body)/='EARTH' .or. trim(config%reference_frame)/='J2000') &
        call fatal('dynamics configuration must use EARTH / J2000')
    call str2et(epoch_text,epoch0)
    if(seen(4)) then
        call str2et(trim(end_utc),end_et)
        duration=end_et-epoch0
    else
        duration=3600*hours
        end_et=epoch0+duration
    end if
    if(.not.all(ieee_is_finite([epoch0,end_et,duration])) .or. duration<0) &
        call fatal('end time must be finite and at or after the OPM epoch')
    if(trim(config%srp_model)=='time_varying_modes' .and. config%use_srp) then
        call load_srp_mode_table(trim(config%srp_mode_table_file),trim(config%srp_mode_metadata_file),epoch0,io,message)
        if(io/=0) call fatal(trim(message))
    end if
    options%rel_tol=config%rkf78_rel_tol
    options%abs_tol=config%rkf78_abs_tol
    options%dt_min=config%rkf78_min_step
    options%dt_max=config%rkf78_max_step
    options%max_steps=config%max_propagation_steps
    if(options%max_steps<2) call fatal('max_propagation_steps must be at least 2')
    inquire(file=trim(prefix)//'_final.opm.json',exist=exists)
    if(exists) call fatal('final OPM already exists; choose another output prefix')
    open(newunit=u,file=trim(prefix)//'_ephemeris.csv',status='new',action='write',iostat=io)
    if(io/=0) call fatal('cannot create ephemeris; check output directory or existing files')
    write(u,'(a)',advance='no') 'utc,elapsed_seconds,x_km,y_km,z_km,vx_kms,vy_kms,vz_kms'
    if(has_cov) then
        do i=1,6
            do j=1,i
                write(u,'(a)',advance='no') ','//cov_key(i,j)
            end do
        end do
    end if
    write(u,*)
    mean=initial; current=0; step=0
    call write_row()
    do while(current<duration)
        next_time=min(current+3600.0_DP,duration)
        if(next_time<=current) call fatal('propagation duration exceeds time resolution')
        call advance_real_probes(points,eta,epoch0,current,next_time,options,io,message)
        if(io/=0) call fatal(trim(message))
        if(has_cov .and. size(weights)>1) then
            call reconstruct_ut_moments(initial_points,points,weights,ut_mean,ut_cov,cross_cov)
            mean=ut_mean; covariance=ut_cov
        else
            mean=points(:,1)
        end if
        if(.not.all(ieee_is_finite(mean)) .or. .not.all(ieee_is_finite(covariance))) &
            call fatal('propagation produced nonfinite state or covariance')
        current=next_time; step=step+1
        call write_row()
        if(mod(step,24)==0 .or. current==duration) then
            write(output_unit,'(a,f12.4)') 'run_orbit_propagation: hours=',current/3600
            flush(output_unit)
        end if
    end do
    close(u,iostat=io)
    if(io/=0) call fatal('cannot close ephemeris')
    call write_final_opm()
    write(*,'(a)') 'Ephemeris: '//trim(prefix)//'_ephemeris.csv'
    write(*,'(a)') 'Final OPM: '//trim(prefix)//'_final.opm.json'
contains
    subroutine fatal(text)
        character(len=*),intent(in) :: text
        write(error_unit,'(a)') 'run_orbit_propagation: '//text
        stop 2
    end subroutine

    subroutine arguments()
        integer :: k,rc
        k=1
        do while(k<=command_argument_count())
            call get_command_argument(k,arg,status=rc)
            if(rc/=0) call fatal('argument too long')
            if(trim(arg)=='--help' .or. trim(arg)=='-h') then
                write(*,'(a)') 'run_orbit_propagation --opm FILE (--hours H | --end-utc UTC) [options]'
                write(*,'(a)') '  --config FILE      default config/config.txt; controls all forces'
                write(*,'(a)') '  --output PREFIX    default input path without extension + _propagated'
                write(*,'(a)') '  --kernel-dir DIR   optional SPICE kernel root'
                write(*,'(a)') '  --seconds S        alternative duration in seconds'
                write(*,'(a)') 'JSON OPM: Earth-centered J2000/GCRF, UTC, km and km/s.'
                write(*,'(a)') 'Covariance present: UT (kappa=0.5); absent: deterministic RKF78.'
                write(*,'(a)') 'Outputs: PREFIX_ephemeris.csv and PREFIX_final.opm.json (no overwrite).'
                stop
            end if
            if(k==command_argument_count()) call fatal('missing option value')
            call get_command_argument(k+1,value,status=rc)
            if(rc/=0 .or. len_trim(value)==0) call fatal('invalid option value')
            select case(trim(arg))
            case('--opm','-opm'); idx=1; opm=value
            case('--config','-cfg'); idx=2; cfg=value
            case('--hours'); idx=3; read(value,*,iostat=rc) hours
            case('--seconds','-dt'); idx=3; read(value,*,iostat=rc) hours; hours=hours/3600
            case('--end-utc','-et'); idx=4; end_utc=value
            case('--output','-o'); idx=5; prefix=value
            case('--kernel-dir'); idx=6; kernel_dir=value
            case default; call fatal('unknown option '//trim(arg))
            end select
            if(rc/=0) call fatal('invalid numeric option')
            if(seen(idx)) call fatal('duplicate option '//trim(arg))
            seen(idx)=.true.; k=k+2
        end do
        if(.not.seen(1)) call fatal('--opm is required')
        if(seen(3) .eqv. seen(4)) call fatal('specify exactly one of --hours, --seconds or --end-utc')
        if(seen(3)) then
            if(.not.ieee_is_finite(hours) .or. hours<0) call fatal('duration must be finite and nonnegative')
        end if
        if(.not.seen(5)) then
            prefix=opm
            k=len_trim(prefix)
            if(k>5) then
                if(prefix(k-4:k)=='.json') prefix=prefix(:k-5)
            end if
            k=len_trim(prefix)
            if(k>4) then
                if(prefix(k-3:k)=='.opm') prefix=prefix(:k-4)
            end if
            prefix=trim(prefix)//'_propagated'
        end if
    end subroutine

    subroutine check_frame(key,allowed1,allowed2)
        character(len=*),intent(in) :: key,allowed1,allowed2
        call json_member(source,key,raw,ok)
        if(.not.ok) return
        call json_string(source,key,frame,ok)
        if(.not.ok) call fatal('invalid OPM '//key)
        if(frame/=allowed1 .and. frame/=allowed2) call fatal('unsupported OPM '//key//': '//frame)
    end subroutine

    function cov_key(row,col) result(key)
        integer,intent(in) :: row,col
        character(len=:),allocatable :: key
        key='C'//trim(keys(row))//'_'//trim(keys(col))
    end function

    subroutine write_row()
        call et2utc(epoch0+current,'ISOC',6,utc)
        write(u,'(a,7(",",es25.16e3))',advance='no') trim(utc),current,mean
        if(has_cov) then
            do i=1,6
                do j=1,i
                    write(u,'(",",es25.16e3)',advance='no') covariance(i,j)
                end do
            end do
        end if
        write(u,*,iostat=io)
        if(io/=0) call fatal('cannot write ephemeris')
        flush(u)
    end subroutine

    subroutine write_final_opm()
        open(newunit=u,file=trim(prefix)//'_final.opm.json',status='new',action='write',iostat=io)
        if(io/=0) call fatal('cannot create final OPM')
        write(u,'(a)') '{'
        ! Preserve identity and spacecraft metadata, never stale fit/GMM results.
        do i=1,size(metadata)
            call json_member(source,trim(metadata(i)),raw,ok)
            if(ok) write(u,'(a)') '  "'//trim(metadata(i))//'": '//raw//','
        end do
        write(u,'(a)') '  "CENTER_NAME": "EARTH",'
        write(u,'(a)') '  "REF_FRAME": "J2000",'
        write(u,'(a)') '  "TIME_SYSTEM": "UTC",'
        call et2utc(end_et,'ISOC',6,utc)
        write(u,'(a)') '  "EPOCH": "'//trim(utc)//'",'
        do i=1,6
            write(u,'(a,es25.16e3,a)') '  "'//trim(keys(i))//'": ',mean(i),','
        end do
        if(has_cov) then
            write(u,'(a)') '  "COV_REF_FRAME": "J2000",'
            do i=1,6
                do j=1,i
                    write(u,'(a,es25.16e3,a)') '  "'//cov_key(i,j)//'": ',covariance(i,j),','
                end do
            end do
            write(u,'(a)') '  "PROPAGATION_METHOD": "UT",'
            write(u,'(a)') '  "UT_KAPPA": 0.5,'
        else
            write(u,'(a)') '  "PROPAGATION_METHOD": "DETERMINISTIC",'
        end if
        write(u,'(a)') '  "STATUS": "PROPAGATED"'
        write(u,'(a)') '}'
        close(u,iostat=io)
        if(io/=0) call fatal('cannot close final OPM')
    end subroutine
end program run_orbit_propagation
