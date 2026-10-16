!> Audited full-sun mode table. Owns ordinary Real storage only, never DA handles.
module pod_srp_mode_module
    use pod_global, only: DP
    use pod_config, only: config
    use pod_json_object_module
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: iostat_eor
    implicit none
    private
    public :: load_srp_mode_table, clear_srp_mode_table, srp_mode_table_is_loaded, interpolate_srp_modes
    real(DP),allocatable,save :: epochs(:),earth(:,:),moon(:,:)
    real(DP),save :: physical(4)=0.0_DP
    real(DP),parameter :: time_tolerance=1e-6_DP, identity_tolerance=1e-12_DP
    character(len=*),parameter :: header= &
        'time_hours,epoch_tdb_s,ref_x_km,ref_y_km,ref_z_km,ref_vx_km_s,ref_vy_km_s,ref_vz_km_s,'// &
        'a_cb_x_km_s2,a_cb_y_km_s2,a_cb_z_km_s2,a_cb_norm_km_s2,'// &
        'a_bwe_x_km_s2,a_bwe_y_km_s2,a_bwe_z_km_s2,a_bwm_x_km_s2,a_bwm_y_km_s2,a_bwm_z_km_s2,'// &
        'eE_x,eE_y,eE_z,eE_norm,eM_x,eM_y,eM_z,eM_norm'
contains
    subroutine clear_srp_mode_table()
        if(allocated(epochs)) deallocate(epochs,earth,moon)
        physical=0.0_DP
    end subroutine

    logical function srp_mode_table_is_loaded()
        srp_mode_table_is_loaded=allocated(epochs)
    end function

    function current_physical() result(values)
        real(DP) :: values(4),pressure
        pressure=config%srp_pressure_1au_n_m2
        if(pressure<0.0_DP) pressure=1367.0_DP/299792458.0_DP
        values=[config%srp_cannonball_cr,config%srp_cannonball_area_m2,config%srp_mass_kg,pressure]
    end function

    logical function physical_matches(values)
        real(DP),intent(in) :: values(4)
        real(DP) :: actual(4)
        actual=current_physical()
        physical_matches=all(ieee_is_finite(actual)) .and. all(actual>0.0_DP)
        if(physical_matches) physical_matches=all(abs(actual-values)<=identity_tolerance*abs(values))
    end function

    subroutine load_srp_mode_table(csv_file,meta_file,expected_epoch0,status,message)
        character(len=*),intent(in) :: csv_file,meta_file
        real(DP),intent(in) :: expected_epoch0
        integer,intent(out) :: status
        character(len=*),intent(out) :: message
        character(len=:),allocatable :: text,units,params,str
        character(len=16384) :: line
        real(DP) :: epoch0,duration,cadence,count,values(4),row(26),cbnorm
        real(DP),allocatable :: times(:),e(:,:),m(:,:)
        integer :: unit,ios,n,i,j,bytes
        logical :: ok,flag
        call clear_srp_mode_table()
        status=1; message='invalid SRP mode metadata'
        if(.not.ieee_is_finite(expected_epoch0)) return
        open(newunit=unit,file=meta_file,status='old',action='read',access='stream',form='unformatted',iostat=ios)
        if(ios/=0) then
            message='cannot open SRP mode metadata'; return
        end if
        inquire(unit=unit,size=bytes,iostat=ios)
        if(ios/=0 .or. bytes<=0 .or. bytes>65536) then
            close(unit); return
        end if
        allocate(character(len=bytes) :: text)
        read(unit,iostat=ios) text
        close(unit)
        if(ios/=0) return
        if(.not.string_equals(text,'schema_version','boxwing_srp_modes_v1')) return
        call json_boolean(text,'audit_complete',flag,ok)
        if(.not.ok .or. .not.flag) return
        call json_boolean(text,'gate_passed',flag,ok)
        if(.not.ok .or. .not.flag) return
        call json_boolean(text,'opm_epoch_verified',flag,ok)
        if(.not.ok .or. .not.flag) return
        if(.not.string_equals(text,'frame','J2000')) return
        if(.not.string_equals(text,'center','EARTH')) return
        if(.not.string_equals(text,'timescale','TDB')) return
        if(.not.string_equals(text,'illumination','full_sun')) return
        if(.not.string_equals(text,'reference_model','cannonball')) return
        if(.not.string_equals(text,'evaluation_contract','common_reference_state')) return
        call json_member(text,'units',units,ok)
        if(.not.ok) return
        if(.not.string_equals(units,'position','km')) return
        if(.not.string_equals(units,'velocity','km/s')) return
        if(.not.string_equals(units,'acceleration','km/s2')) return
        if(.not.string_equals(units,'modes','dimensionless')) return
        call json_number(text,'epoch0_tdb_s',epoch0,ok)
        if(.not.ok) return
        if(abs(epoch0-expected_epoch0)>time_tolerance) then
            message='SRP mode OPM epoch mismatch'; return
        end if
        call json_number(text,'duration_hours',duration,ok)
        if(.not.ok) return
        call json_number(text,'sample_hours',cadence,ok)
        if(.not.ok) return
        call json_number(text,'row_count',count,ok)
        if(.not.ok) return
        if(duration<360.0_DP .or. cadence<=0.0_DP .or. count<2 .or. count>1000000) return
        n=int(count)
        if(count/=real(n,DP)) return
        if(abs(real(n-1,DP)*cadence-duration)*3600.0_DP>time_tolerance) return
        call json_member(text,'physical_parameters',params,ok)
        if(.not.ok) return
        call json_number(params,'srp_cannonball_cr',values(1),ok)
        if(.not.ok) return
        call json_number(params,'srp_cannonball_area_m2',values(2),ok)
        if(.not.ok) return
        call json_number(params,'srp_mass_kg',values(3),ok)
        if(.not.ok) return
        call json_number(params,'srp_pressure_1au_n_m2',values(4),ok)
        if(.not.ok) return
        if(any(values<=0.0_DP) .or. .not.physical_matches(values)) then
            message='SRP mode Cannonball physical parameters mismatch'; return
        end if
        ! SHA recomputation is the Python preflight's responsibility, not this numerical reader.
        do i=1,3
            select case(i)
            case(1); call json_string(text,'csv_sha256',str,ok)
            case(2); call json_string(text,'opm_sha256',str,ok)
            case(3); call json_string(text,'config_sha256',str,ok)
            end select
            if(.not.ok .or. len(str)/=64) return
            do j=1,len(str)
                if(index('0123456789abcdefABCDEF',str(j:j))==0) return
            end do
        end do
        message='invalid SRP mode CSV'
        open(newunit=unit,file=csv_file,status='old',action='read',iostat=ios)
        if(ios/=0) then
            message='cannot open SRP mode CSV'; return
        end if
        call read_csv_record(unit,line,ios)
        if(ios/=0) then
            close(unit); return
        end if
        if(trim(line)/=header) then
            close(unit); return
        end if
        allocate(times(n),e(3,n),m(3,n))
        do i=1,n
            call read_csv_record(unit,line,ios)
            if(ios/=0) exit
            if(count_commas(trim(line))/=25) then
                ios=1; exit
            end if
            ! List-directed null fields and slash termination are not valid numeric CSV.
            if(index(line,'/')/=0 .or. index(line,'*')/=0 .or. index(line,';')/=0) then
                ios=1; exit
            end if
            call parse_numeric_row(trim(line),row,ok)
            if(.not.ok .or. .not.all(ieee_is_finite(row))) then
                ios=1; exit
            end if
            if(abs(row(1)-real(i-1,DP)*cadence)*3600.0_DP>time_tolerance .or. &
               abs(row(2)-(epoch0+row(1)*3600.0_DP))>time_tolerance) then
                ios=1; exit
            end if
            if(i>1) then
                if(row(2)<=times(i-1)) then
                    ios=1; exit
                end if
            end if
            cbnorm=norm2(row(9:11))
            if(.not.ieee_is_finite(cbnorm) .or. cbnorm<=0.0_DP .or. &
               abs(row(12)-cbnorm)>identity_tolerance*cbnorm) then
                ios=1; exit
            end if
            if(norm2((row(13:15)-row(9:11))/cbnorm-row(19:21))>identity_tolerance .or. &
               norm2((row(16:18)-row(9:11))/cbnorm-row(23:25))>identity_tolerance .or. &
               abs(norm2(row(19:21))-row(22))>identity_tolerance .or. &
               abs(norm2(row(23:25))-row(26))>identity_tolerance .or. &
               row(22)>0.3_DP .or. row(26)>0.3_DP) then
                ios=1; exit
            end if
            times(i)=row(2); e(:,i)=row(19:21); m(:,i)=row(23:25)
        end do
        if(ios==0) then
            call read_csv_record(unit,line,ios)
            if(ios<0) then
                call move_alloc(times,epochs); call move_alloc(e,earth); call move_alloc(m,moon)
                physical=values; status=0; message=''
            end if
        end if
        close(unit)
    end subroutine

    logical function string_equals(text,key,expected)
        character(len=*),intent(in) :: text,key,expected
        character(len=:),allocatable :: actual
        logical :: ok
        call json_string(text,key,actual,ok)
        string_equals=ok .and. actual==expected
    end function

    ! Reject overlong records instead of silently discarding their unparsed suffix.
    subroutine read_csv_record(unit,line,status)
        integer,intent(in) :: unit
        character(len=*),intent(out) :: line
        integer,intent(out) :: status
        read(unit,'(A)',advance='no',iostat=status) line
        if(status==0) then
            status=1
        else if(status==iostat_eor) then
            status=0
        end if
    end subroutine

    subroutine parse_numeric_row(line,row,ok)
        character(len=*),intent(in) :: line
        real(DP),intent(out) :: row(26)
        logical,intent(out) :: ok
        integer :: first,last,i,comma
        first=1; row=0; ok=.false.
        do i=1,26
            comma=index(line(first:),',')
            last=len(line)
            if(comma>0) last=first+comma-2
            if(last<first) return
            call json_number('{"value":'//trim(adjustl(line(first:last)))//'}','value',row(i),ok)
            if(.not.ok) return
            first=last+2
        end do
    end subroutine

    integer function count_commas(text) result(n)
        character(len=*),intent(in) :: text
        integer :: i
        n=0
        do i=1,len(text)
            if(text(i:i)==',') n=n+1
        end do
    end function

    subroutine interpolate_srp_modes(epoch,eE,eM,status,message)
        real(DP),intent(in) :: epoch
        real(DP),intent(out) :: eE(3),eM(3)
        integer,intent(out) :: status
        character(len=*),intent(out) :: message
        integer :: lo,hi,mid,n
        real(DP) :: t,f
        eE=0.0_DP; eM=0.0_DP; status=1; message='SRP mode table not loaded'
        if(.not.allocated(epochs)) return
        message='SRP mode physical identity changed'
        if(.not.physical_matches(physical)) return
        message='SRP mode time outside table coverage'
        if(.not.ieee_is_finite(epoch)) return
        n=size(epochs)
        if(epoch<epochs(1)-time_tolerance .or. epoch>epochs(n)+time_tolerance) return
        t=max(epochs(1),min(epochs(n),epoch))
        lo=1; hi=n
        do while(hi-lo>1)
            mid=(lo+hi)/2
            if(t<epochs(mid)) then
                hi=mid
            else
                lo=mid
            end if
        end do
        f=(t-epochs(lo))/(epochs(hi)-epochs(lo))
        eE=(1.0_DP-f)*earth(:,lo)+f*earth(:,hi)
        eM=(1.0_DP-f)*moon(:,lo)+f*moon(:,hi)
        status=0; message=''
    end subroutine
end module pod_srp_mode_module
