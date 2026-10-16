program test_srp_mode_table
    use pod_global, only: DP
    use pod_config, only: config, set_default_config
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
    use pod_srp_mode_module, only: load_srp_mode_table, clear_srp_mode_table, &
        srp_mode_table_is_loaded, interpolate_srp_modes
    implicit none
    character(len=256) :: prefix, csv_file, meta_file, message
    real(DP) :: eE(3), eM(3)
    integer :: status, clock_value, u, variant,failures

    call set_default_config()
    config%srp_pressure_1au_n_m2=4.559821181419152e-6_DP
    call system_clock(clock_value)
    write(prefix,'(A,I0)') '/tmp/pod_srp_mode_table_',clock_value
    csv_file=trim(prefix)//'.csv'
    meta_file=trim(prefix)//'.json'
    call clear_srp_mode_table()
    call interpolate_srp_modes(1000.0_DP,eE,eM,status,message)
    call require(status/=0 .and. .not.srp_mode_table_is_loaded(),'unloaded interpolation must fail')
    call write_fixture(0)
    call load_srp_mode_table(trim(csv_file),trim(meta_file),1000.0_DP,status,message)
    call require(status==0 .and. srp_mode_table_is_loaded(),'load valid table: '//trim(message))
    call interpolate_srp_modes(1000.0_DP,eE,eM,status,message)
    call require(status==0 .and. abs(eE(1)-0.1_DP)<1e-15_DP,'first node')
    call require(abs(eM(2)-0.2_DP)<1e-15_DP,'Moon mode')
    call interpolate_srp_modes(2800.0_DP,eE,eM,status,message)
    call require(status==0 .and. abs(eE(1)-(0.1_DP+0.5_DP/3600.0_DP))<1e-15_DP,'midpoint')
    call interpolate_srp_modes(1000.0_DP-5e-7_DP,eE,eM,status,message)
    call require(status==0 .and. abs(eE(1)-0.1_DP)<1e-15_DP,'roundoff clamps first node')
    call interpolate_srp_modes(1000.0_DP-2e-6_DP,eE,eM,status,message)
    call require(status/=0,'left extrapolation rejected')
    call interpolate_srp_modes(1297000.0_DP+5e-7_DP,eE,eM,status,message)
    call require(status==0 .and. abs(eE(1)-0.2_DP)<1e-15_DP,'roundoff clamps last node')
    call interpolate_srp_modes(1297000.0_DP+2e-6_DP,eE,eM,status,message)
    call require(status/=0,'right extrapolation rejected')
    config%srp_mass_kg=1201.0_DP
    call interpolate_srp_modes(1000.0_DP,eE,eM,status,message)
    call require(status/=0,'runtime physical identity drift rejected')
    call load_srp_mode_table(trim(csv_file),trim(meta_file),1000.0_DP,status,message)
    call require(status/=0 .and. .not.srp_mode_table_is_loaded(),'failed reload clears old table')
    config%srp_mass_kg=1200.0_DP
    call load_srp_mode_table(trim(csv_file),trim(meta_file),1001.0_DP,status,message)
    call require(status/=0,'OPM epoch mismatch')
    failures=0
    do variant=1,25
        call write_fixture(variant)
        call load_srp_mode_table(trim(csv_file),trim(meta_file),1000.0_DP,status,message)
        if(status==0 .or. srp_mode_table_is_loaded()) then
            write(*,'(A,I0)') 'accepted invalid fixture=',variant
            failures=failures+1
        end if
    end do
    call require(failures==0,'all invalid fixtures rejected')
    call write_fixture(0)
    call load_srp_mode_table(trim(csv_file),trim(meta_file),1000.0_DP,status,message)
    call require(status==0,'valid reload after failures')
    call load_srp_mode_table('/nonexistent/pod_mode.csv',trim(meta_file),1000.0_DP,status,message)
    call require(status/=0 .and. .not.srp_mode_table_is_loaded(),'missing CSV invalidates loaded table')
    call load_srp_mode_table(trim(csv_file),'/nonexistent/pod_meta.json',1000.0_DP,status,message)
    call require(status/=0 .and. .not.srp_mode_table_is_loaded(),'missing metadata rejected')
    call clear_srp_mode_table()
    call clear_srp_mode_table()
    call require(.not.srp_mode_table_is_loaded(),'idempotent clear')
    open(newunit=u,file=trim(csv_file),status='old'); close(u,status='delete')
    open(newunit=u,file=trim(meta_file),status='old'); close(u,status='delete')
    write(*,'(A)') 'PASS: SRP mode table parsing, identity, interpolation and lifecycle'

contains
    subroutine require(ok,text)
        logical,intent(in) :: ok
        character(len=*),intent(in) :: text
        if(.not.ok) then
            write(*,'(A)') 'FAIL: '//trim(text)
            error stop 1
        end if
    end subroutine

    subroutine write_fixture(bad)
        integer,intent(in) :: bad
        integer :: unit, i
        real(DP) :: h, row(26), a(3), earth(3), moon(3)
        character(len=:),allocatable :: text
        character(len=*),parameter :: header= &
            'time_hours,epoch_tdb_s,ref_x_km,ref_y_km,ref_z_km,ref_vx_km_s,ref_vy_km_s,ref_vz_km_s,'// &
            'a_cb_x_km_s2,a_cb_y_km_s2,a_cb_z_km_s2,a_cb_norm_km_s2,'// &
            'a_bwe_x_km_s2,a_bwe_y_km_s2,a_bwe_z_km_s2,a_bwm_x_km_s2,a_bwm_y_km_s2,a_bwm_z_km_s2,'// &
            'eE_x,eE_y,eE_z,eE_norm,eM_x,eM_y,eM_z,eM_norm'
        open(newunit=unit,file=trim(csv_file),status='replace')
        if(bad==1) then
            write(unit,'(A)') 'wrong_header'
        else
            write(unit,'(A)') header
        end if
        a=[1e-10_DP,0.0_DP,0.0_DP]
        do i=0,360
            h=real(i,DP)
            if(bad==2 .and. i==2) h=1.0_DP
            earth=[0.1_DP+h/3600.0_DP,0.0_DP,0.0_DP]
            moon=[0.0_DP,0.2_DP,0.0_DP]
            row=[h,1000.0_DP+3600.0_DP*h,300000.0_DP,1000.0_DP,20000.0_DP, &
                 0.1_DP,0.5_DP,0.2_DP,a,norm2(a),a+norm2(a)*earth,a+norm2(a)*moon, &
                 earth,norm2(earth),moon,norm2(moon)]
            if(bad==3 .and. i==1) row(22)=row(22)+0.1_DP
            if(bad==4 .and. i==1) row(19)=0.3_DP
            if(bad==19 .and. i==1) row(2)=row(2)+1.0_DP
            if(bad==20 .and. i==1) then
                row(19:21)=[0.4_DP,0.0_DP,0.0_DP]; row(22)=0.4_DP
                row(13:15)=a+norm2(a)*row(19:21)
            end if
            if(bad==21 .and. i==1) row(12)=ieee_value(0.0_DP,ieee_quiet_nan)
            if(bad==22 .and. i==1) row(12)=0.0_DP
            if(bad==23 .and. i==1) row(26)=0.1_DP
            if(bad==24 .and. i==1) then
                row(9:18)=huge(1.0_DP)*0.75_DP; row(19:26)=0
            end if
            if(bad==9 .and. i==1) then
                write(unit,'(*(ES25.17E3,:,","))',advance='no') row(:25)
                write(unit,'(A)') ',0.2 999' ! Two numbers in the last field must not be silently accepted.
            else if(bad==25 .and. i==1) then
                write(unit,'(*(ES25.17E3,:,","))',advance='no') row
                write(unit,'(A)') repeat(' ',17000)//'garbage'
            else
                write(unit,'(*(ES25.17E3,:,","))') row
            end if
        end do
        close(unit)
        text='{"schema_version":"boxwing_srp_modes_v1","audit_complete":true,'// &
            '"gate_passed":true,"opm_epoch_verified":true,"frame":"J2000",'// &
            '"center":"EARTH","timescale":"TDB","illumination":"full_sun",'// &
            '"reference_model":"cannonball","evaluation_contract":"common_reference_state",'// &
            '"epoch0_tdb_s":1000,"duration_hours":360,"sample_hours":1,"row_count":361,'// &
            '"units":{"position":"km","velocity":"km/s","acceleration":"km/s2","modes":"dimensionless"},'// &
            '"physical_parameters":{"srp_cannonball_cr":1.2213,"srp_cannonball_area_m2":30,'// &
            '"srp_mass_kg":1200,"srp_pressure_1au_n_m2":4.559821181419152e-6},'// &
            '"csv_sha256":"'//repeat('a',64)//'","opm_sha256":"'//repeat('b',64)// &
            '","config_sha256":"'//repeat('c',64)//'"}'
        if(bad==5) text='{"schema_version":"boxwing_srp_modes_v1","audit_complete":false}'
        if(bad==6) text=text(:len(text)-1)//',"frame":"RTN"}'
        if(bad==7) text=text(:len(text)-1) ! Missing final brace.
        if(bad==8) text='{"unused":{"frame":"J2000"},"audit_complete":true}'
        if(bad==10) text=text//repeat(' ',17000)//'garbage'
        if(bad==11) text=replace_text(text,'"km/s2"','"m/s2"')
        if(bad==12) text=replace_text(text,'"sample_hours":1','"sample_hours":0.5')
        if(bad==13) text=replace_text(text,'"duration_hours":360','"duration_hours":359')
        if(bad==14) text=replace_text(text,'"row_count":361','"row_count":361.5')
        if(bad==15) text=replace_text(text,'"audit_complete":true','"audit_complete":"true"')
        if(bad==16) text=replace_text(text,'"epoch0_tdb_s":1000','"epoch0_tdb_s":NaN')
        if(bad==17) text=replace_text(text,'"schema_version":"boxwing_srp_modes_v1"', &
                                         '"schema_version":"wrong_schema"')
        if(bad==18) text=replace_text(text,repeat('a',64),repeat('g',64))
        open(newunit=unit,file=trim(meta_file),status='replace')
        write(unit,'(A)') text
        close(unit)
    end subroutine

    function replace_text(text,old,new) result(value)
        character(len=*),intent(in) :: text,old,new
        character(len=:),allocatable :: value
        integer :: p
        p=index(text,old)
        call require(p>0,'fixture mutation has target')
        value=text(:p-1)//new//text(p+len(old):)
    end function
end program test_srp_mode_table
