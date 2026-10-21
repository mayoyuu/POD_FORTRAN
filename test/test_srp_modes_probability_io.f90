program test_srp_modes_probability_io
    use pod_global, only: DP
    use pod_srp_modes_probability_io_module
    use pod_config, only: config
    use pod_engine_module, only: pod_engine_init
    use pod_srp_mode_module, only: load_srp_mode_table,clear_srp_mode_table
    use pod_uq_ads_coordinates_module, only: ADS_PARAM_SRP_MODES,ads_unit_to_physical,ads_physical_parameters_type
    use pod_uq_hfem_ads_module
    use pod_ads_split_module, only: patch_type,patch_split,patch_destroy,mf_destroy,mf_init,mf_push,mf_find_patch
    use pod_dace_classes, only: DA,da_add,da_mul,active_da_count
    implicit none
    type(probability_population_type) :: population,loaded
    type(probability_output_type) :: output
    real(DP) :: states(6,8),nominal(6),q(8),z(6)
    integer :: i,j,status,unit,clock_value
    integer :: owners(8)
    integer,allocatable :: ids(:)
    character(len=512) :: prefix
    character(len=32) :: scenario
    character(len=4096) :: line
    call get_command_argument(1,prefix)
    call get_command_argument(2,scenario)
    if(len_trim(prefix)==0) then
        call system_clock(clock_value)
        write(prefix,'(A,I0)') '/tmp/pod_probability_io_',clock_value
    end if
    call probability_choose_ids(8,3,ids,status)
    call require(status==0.and.all(ids==[1,5,8]),'fixed IDs')
    allocate(population%q(8,8),population%modes(2,8))
    population%q=0; population%modes=.5_DP
    open(newunit=unit,file=trim(prefix)//'_samples.csv',status='new')
    write(unit,'(A)') 'schema=srp_modes_probability_samples_v1'
    write(unit,'(A)') PROBABILITY_SAMPLE_HEADER
    do i=1,8
        q=real(i-4,DP)/4; z=3*q(1:6)
        write(unit,'(I0,18(",",ES25.16E3))') i,z,q,(1+q(7:8))/2,(1+q(7))/2, &
            (1-(1+q(7))/2)*(1+q(8))/2
    end do
    close(unit)
    call load_probability_samples(trim(prefix)//'_samples.csv',8,loaded,status)
    call require(status==0.and.size(loaded%q,2)==8,'strict 8D loader')
    call require(abs(loaded%q(8,8)-1.0_DP)<1e-14_DP,'global mode coordinate')
    call load_probability_samples(trim(prefix)//'_samples.csv',9,loaded,status)
    call require(status/=0,'missing sample rejected')
    nominal=[300000.0_DP,10000.0_DP,-20000.0_DP,.1_DP,.2_DP,.3_DP]
    do i=1,8
        states(:,i)=nominal+[real(i,DP),sin(real(i,DP)),cos(real(i,DP)), &
            .01_DP*real(i,DP),.02_DP*sin(real(i,DP)),0.0_DP]
    end do
    if(trim(scenario)=='nearzero') then
        nominal(6)=0
        states(6,:)=[(real(i,DP)*1e-16_DP,i=1,8)]
    end if
    owners=1
    if(trim(scenario)=='split') call split_population()
    call probability_open(output,trim(prefix),population,2,status)
    call require(status==0,'open output')
    do j=0,2
        call probability_write_hour(output,population,real(j,DP),1234.0_DP+3600*j,states,nominal,owners, &
            .not.(trim(scenario)=='events'.and.j==2),.not.(trim(scenario)=='events'.and.j==2), &
            trim(scenario)=='events'.and.j==1,trim(scenario)=='events'.and.j==2, &
            trim(scenario)=='events'.and.j==2,merge(4,1,trim(scenario)=='split'), &
            merge(1,0,trim(scenario)=='events'.and.j==2),0.0_DP,0.0_DP,0.0_DP,0.0_DP,j==2,status)
        call require(status==0,'write hourly moments')
    end do
    call probability_close(output,status)
    call require(status==0,'close output')
    open(newunit=unit,file=trim(prefix)//'_probability_hourly_stats.csv',status='old')
    read(unit,'(A)') line
    call require(index(line,'cov_14')>0.and.index(line,'cov_66')>0,'full 6x6 covariance header')
    call require(index(line,'skewness_valid6')>0,'degenerate moment flags')
    close(unit)
    write(*,'(A)') 'PASS: probability input, full moments and layered native output'
contains
    subroutine split_population()
        type(hfem_ads_history_type) :: history
        type(hfem_ads_options_type) :: options
        type(ads_physical_parameters_type) :: parameters
        type(patch_type) :: a,b,left,right
        type(DA) :: variable
        real(DP) :: cov(6,6),deviation(6),expected(6,8)
        integer :: baseline,i
        logical :: found(8)
        character(len=512) :: message
        call pod_engine_init('config/config_srp_modes_DRO-2.txt')
        call load_srp_mode_table(trim(config%srp_mode_table_file),trim(config%srp_mode_metadata_file), &
            852033729.1839125_DP,status,message)
        call require(status==0,'split mode table')
        options%parameterization=ADS_PARAM_SRP_MODES; cov=0
        baseline=active_da_count()
        call hfem_ads_history_init(history,nominal,cov,852033729.1839125_DP,options,status,message)
        call require(status==0,'split probability history')
        call variable%init_var(7)
        call da_add(history%domain%patches(1)%da_vec%elements(1),variable,history%domain%patches(1)%da_vec%elements(1))
        call variable%destroy(); call variable%init_var(8)
        call da_mul(variable,2.0_DP,variable)
        call da_add(history%domain%patches(1)%da_vec%elements(1),variable,history%domain%patches(1)%da_vec%elements(1))
        call variable%destroy()
        call patch_split(history%domain%patches(1),8,a,b)
        call mf_destroy(history%domain); call mf_init(history%domain)
        call patch_split(a,7,left,right); call mf_push(history%domain,left); call mf_push(history%domain,right)
        call patch_destroy(left); call patch_destroy(right)
        call patch_split(b,7,left,right); call mf_push(history%domain,left); call mf_push(history%domain,right)
        call patch_destroy(left); call patch_destroy(right); call patch_destroy(a); call patch_destroy(b)
        population%q(7,:)=[-.9_DP,.9_DP,-.7_DP,.7_DP,0.0_DP,0.0_DP,-.4_DP,.4_DP]
        population%q(8,:)=[-.8_DP,-.8_DP,.8_DP,.8_DP,0.0_DP,.5_DP,0.0_DP,0.0_DP]
        population%modes=(1+population%q(7:8,:))/2
        expected=spread(nominal,2,8)
        expected(1,:)=expected(1,:)+population%q(7,:)+2*population%q(8,:)
        call hfem_ads_history_evaluate(history,population%q,states,found,status)
        call require(status==0.and.all(found).and.maxval(abs(states-expected))<1e-9_DP,'split global-coordinate batch values')
        do i=1,8
            call mf_find_patch(history%domain,population%q(:,i),owners(i))
            call ads_unit_to_physical(history%coordinate_map,population%q(:,i),deviation,parameters,status)
            call require(status==0.and.abs(parameters%s-population%modes(1,i))<1e-14_DP.and. &
                abs(parameters%u-population%modes(2,i))<1e-14_DP,'split does not relabel global physical modes')
        end do
        call require(all([(any(owners==i),i=1,4)]),'all four split patches covered')
        call hfem_ads_history_destroy(history); call clear_srp_mode_table()
        call require(active_da_count()==baseline,'split probability evaluation has no DA leak')
    end subroutine
    subroutine require(ok,label)
        logical,intent(in) :: ok
        character(len=*),intent(in) :: label
        if(.not.ok) then
            write(*,'(A)') 'FAIL: '//label
            error stop 1
        end if
    end subroutine
end program
