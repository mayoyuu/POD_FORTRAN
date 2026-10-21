program test_srp_modes_force_consistency
    use pod_global, only: DP
    use pod_config, only: config,load_config,set_default_config,validate_config
    use pod_spice, only: spice_init
    use pod_engine_module, only: pod_engine_init
    use pod_srp_mode_module
    use pod_force_model_module, only: compute_solar_radiation_pressure,set_srp_mode_coefficients_real, &
        clear_srp_mode_coefficients_real,set_srp_scale_error,clear_srp_scale_error, &
        set_real_ballistic=>set_srp_ballistic_parameters,clear_real_ballistic=>clear_srp_ballistic_parameters
    use pod_da_force_model_module, only: da_compute_solar_radiation_pressure,ForceModelTempPool, &
        set_srp_mode_uncertainty,clear_srp_mode_uncertainty,set_srp_scale_uncertainty,clear_srp_scale_uncertainty, &
        set_da_ballistic=>set_srp_ballistic_parameters,clear_da_ballistic=>clear_srp_ballistic_parameters
    use pod_dace_classes
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
    implicit none
    type(AlgebraicVector) :: position_da,acceleration_da
    type(ForceModelTempPool) :: pool
    type(DA) :: variable
    real(DP) :: row(26),position(3),acceleration(3),expected(3),eE(3),eM(3),coeffs(2,5),q(8)
    real(DP) :: epoch0,epoch,s,u,max_endpoint,max_constant,actual(3),peak,error
    integer :: unit,ios,i,j,k,status,before,warm,clock_value
    character(len=256) :: message,arg
    character(len=512) :: directory,csv_file,meta_file
    character(len=256) :: cfg_file
    logical :: exists
    call get_command_argument(1,arg)
    call spice_init()
    call load_config('config/config_box_wing_earth.txt')
    directory='SRP/boxwing_modes_DRO-2_init/'
    inquire(file=trim(directory)//'boxwing_modes.csv',exist=exists)
    if(.not.exists) directory='../../SRP/boxwing_modes_DRO-2_init/'
    csv_file=trim(directory)//'boxwing_modes.csv'; meta_file=trim(directory)//'boxwing_modes_meta.json'
    call system_clock(clock_value)
    write(cfg_file,'(A,I0,A)') '/tmp/pod_srp_modes_cfg_',clock_value,'.txt'
    open(newunit=unit,file=trim(cfg_file),status='new')
    write(unit,'(A)') 'srp_model = time_varying_modes'
    write(unit,'(A)') 'srp_mode_table_file = '//trim(csv_file)
    write(unit,'(A)') 'srp_mode_metadata_file = '//trim(meta_file)
    write(unit,'(A)') 'srp_pressure_1au_n_m2 = 4.559821181419152e-6'
    close(unit)
    call load_config(trim(cfg_file))
    open(newunit=unit,file=trim(cfg_file),status='old'); close(unit,status='delete')
    call require(trim(config%srp_mode_table_file)==trim(csv_file) .and. &
                 trim(config%srp_mode_metadata_file)==trim(meta_file),'mode config paths are parsed')
    call require(all(config%srp_box_dimensions_m<0),'minimal modes config has no Box-Wing geometry')
    call require(validate_config(),'mode config does not need live Box-Wing geometry')
    config%srp_mode_table_file='/nonexistent/pod_mode.csv'
    call require(.not.validate_config(),'mode config missing file')
    config%srp_mode_table_file=csv_file
    open(newunit=unit,file=trim(csv_file),status='old')
    read(unit,*) ! header
    read(unit,*) row; epoch0=row(2); position=row(3:5); epoch=epoch0
    call dace_initialize(3,8); before=active_da_count()
    call position_da%init(3); call acceleration_da%init(3); call pool%init(3)
    do i=1,3
        call variable%init_var(i)
        call da_mul(variable,100.0_DP,position_da%elements(i))
        call da_add(position_da%elements(i),position(i),position_da%elements(i))
    end do
    call variable%destroy()
    call clear_srp_mode_table(); call clear_srp_mode_coefficients_real(); call clear_srp_mode_uncertainty()
    select case(trim(arg))
    case('invalid-real'); call set_srp_mode_coefficients_real(-0.1_DP,0.5_DP)
    case('nan-real'); call set_srp_mode_coefficients_real(ieee_value(0.0_DP,ieee_quiet_nan),0.5_DP)
    case('invalid-index'); call set_srp_mode_uncertainty(7,0.5_DP,0.5_DP,7,0.5_DP,0.5_DP)
    case('out-of-range-index'); call set_srp_mode_uncertainty(7,0.5_DP,0.5_DP,9,0.5_DP,0.5_DP)
    case('invalid-span'); call set_srp_mode_uncertainty(7,0.1_DP,0.5_DP,8,0.5_DP,0.5_DP)
    case('nan-da'); call set_srp_mode_uncertainty(7,0.5_DP,ieee_value(0.0_DP,ieee_quiet_nan),8,0.5_DP,0.5_DP)
    case('unset-real'); call compute_solar_radiation_pressure(position,epoch,acceleration)
    case('unset-da'); call da_compute_solar_radiation_pressure(position_da,epoch,acceleration_da,pool)
    case('no-table')
        call set_srp_mode_coefficients_real(0.0_DP,0.0_DP)
        call compute_solar_radiation_pressure(position,epoch,acceleration)
    end select
    call load_srp_mode_table(trim(csv_file),trim(meta_file),epoch0,status,message)
    call require(status==0,'P1 audited table: '//trim(message))
    call set_srp_mode_coefficients_real(0.5_DP,0.5_DP)
    call set_srp_mode_uncertainty(7,0.5_DP,0.5_DP,8,0.5_DP,0.5_DP)
    select case(trim(arg))
    case('mixed-real')
        call set_srp_scale_error(0.1_DP)
        call compute_solar_radiation_pressure(position,epoch,acceleration)
    case('mixed-da')
        call set_srp_scale_uncertainty(7,0.1_DP,0.1_DP)
        call da_compute_solar_radiation_pressure(position_da,epoch,acceleration_da,pool)
    case('override-real')
        call set_real_ballistic(Cr=1.2213_DP)
        call compute_solar_radiation_pressure(position,epoch,acceleration)
    case('override-da')
        call set_da_ballistic(Cr=1.2213_DP)
        call da_compute_solar_radiation_pressure(position_da,epoch,acceleration_da,pool)
    case('extrapolation')
        call compute_solar_radiation_pressure(position,epoch0-1.0_DP,acceleration)
    case('init-clears-real','init-clears-da')
        call pod_engine_init('config/config_box_wing_earth.txt')
        call require(.not.srp_mode_table_is_loaded(),'engine clears mode table')
        config%srp_model='time_varying_modes'
        if(trim(arg)=='init-clears-real') call compute_solar_radiation_pressure(position,epoch,acceleration)
        if(trim(arg)=='init-clears-da') call da_compute_solar_radiation_pressure(position_da,epoch,acceleration_da,pool)
    end select
    call require(len_trim(arg)==0,'negative path unexpectedly returned')
    coeffs(:,1)=[0.0_DP,0.0_DP]; coeffs(:,2)=[1.0_DP,0.0_DP]; coeffs(:,3)=[0.0_DP,1.0_DP]
    coeffs(:,4)=[1.0_DP,1.0_DP]; coeffs(:,5)=[0.37_DP,0.81_DP]
    max_endpoint=0; max_constant=0; peak=0
    do k=1,361
        if(k>1) read(unit,*,iostat=ios) row
        position=row(3:5); epoch=row(2)
        do j=1,5
            s=coeffs(1,j); u=coeffs(2,j)
            call set_srp_mode_coefficients_real(s,u)
            call compute_solar_radiation_pressure(position,epoch,acceleration)
            expected=row(9:11)+row(12)*(s*row(19:21)+(1-s)*u*row(23:25))
            if(j==1) expected=row(9:11)
            if(j==2 .or. j==4) expected=row(13:15)
            if(j==3) expected=row(16:18)
            error=norm2(acceleration-expected)/norm2(expected); max_endpoint=max(max_endpoint,error)
            call require(error<1e-12_DP,'361 node endpoint/interior Real identity')
            if(mod(k-1,60)==0) then
                do i=1,3
                    position_da%elements(i)=position(i)
                end do
                call set_srp_mode_uncertainty(7,s,0.0_DP,8,u,0.0_DP)
                call da_compute_solar_radiation_pressure(position_da,epoch,acceleration_da,pool)
                actual=acceleration_da%cons()
                error=norm2(actual-acceleration)/norm2(acceleration); max_constant=max(max_constant,error)
                call require(error<1e-12_DP,'DA/Real constant endpoint/interior')
            end if
        end do
    end do
    close(unit)
    call interpolate_srp_modes(epoch,eE,eM,status,message)
    position=position+[120.0_DP,-80.0_DP,30.0_DP]
    do i=1,3
        call variable%init_var(i)
        call da_mul(variable,1.0_DP,position_da%elements(i))
        call da_add(position_da%elements(i),position(i),position_da%elements(i))
    end do
    call variable%destroy()
    call set_srp_mode_coefficients_real(0.25_DP,0.75_DP)
    call set_srp_mode_uncertainty(7,0.25_DP,0.25_DP,8,0.75_DP,0.25_DP)
    call compute_solar_radiation_pressure(position,epoch,acceleration)
    call da_compute_solar_radiation_pressure(position_da,epoch,acceleration_da,pool)
    actual=acceleration_da%cons()
    error=norm2(actual-acceleration)/norm2(acceleration); max_constant=max(max_constant,error)
    call require(error<1e-12_DP,'off-reference constant')
    q=0; q(7)=0.4_DP; q(8)=-0.6_DP
    actual=acceleration_da%eval(q)
    call set_srp_mode_coefficients_real(0.25_DP+0.25_DP*q(7),0.75_DP+0.25_DP*q(8))
    call compute_solar_radiation_pressure(position,epoch,acceleration)
    call require(norm2(actual-acceleration)/norm2(acceleration)<1e-12_DP,'local coefficients preserve global s/u')
    do j=1,3
        call require(abs(acceleration_da%elements(j)%get_deriv_value(7))>0,'DA s sensitivity')
        call require(abs(acceleration_da%elements(j)%get_deriv_value(8))>0,'DA u sensitivity')
    end do
    warm=active_da_count()
    do k=1,500
        call da_compute_solar_radiation_pressure(position_da,epoch,acceleration_da,pool)
        call require(active_da_count()==warm,'mode force does not leak DA handles')
    end do
    call pool%destroy(); call position_da%destroy(); call acceleration_da%destroy()
    call clear_srp_mode_coefficients_real(); call clear_srp_mode_uncertainty(); call clear_srp_mode_table()
    call require(active_da_count()==before,'all owned DA handles released')
    write(*,'(A,ES12.4)') 'P2_MAX_ENDPOINT_REL=',max_endpoint
    write(*,'(A,ES12.4)') 'P2_MAX_CONSTANT_REL=',max_constant
    write(*,'(A,I0,A,I0)') 'P2_DA_HANDLES: baseline=',before,', final=',active_da_count()
    write(*,'(A)') 'PASS: P1 nodes, 8D proxy forces, local sensitivities and handle ownership'
contains
    subroutine require(ok,label)
        logical,intent(in) :: ok
        character(len=*),intent(in) :: label
        if(.not.ok) then
            write(*,'(A)') 'FAIL: '//label
            error stop 1
        end if
    end subroutine
end program
