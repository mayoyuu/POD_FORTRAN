program test_srp_modes_snapshot
    use pod_global, only: DP
    use pod_config, only: config,load_config
    use pod_srp_mode_module
    use pod_uq_hfem_ads_module
    use pod_uq_ads_coordinates_module, only: ADS_PARAM_SRP_MODES
    use pod_srp_ads_snapshot_module, only: write_ads_patch_snapshot
    use pod_dace_classes, only: DA,active_da_count,da_add,assignment(=)
    use pod_ads_split_module, only: patch_type,patch_split,patch_destroy,mf_destroy,mf_init,mf_push
    implicit none
    type(hfem_ads_history_type) :: history
    type(hfem_ads_options_type) :: options
    type(patch_type) :: left,right
    type(DA) :: variable
    real(DP) :: state(6),cov(6,6),epoch
    integer :: status,clock_value,before,after,unit,i
    character(len=256) :: prefix,message
    character(len=2048) :: line
    character(len=11),parameter :: suffixes(4)=[character(len=11) :: '_meta.csv','_coeff.csv','_splits.csv','_map.csv']
    call load_config('config/config_box_wing_earth.txt')
    config%srp_model='time_varying_modes'
    epoch=852033729.1839125_DP
    call load_srp_mode_table('../../SRP/boxwing_modes_DRO-2_init/boxwing_modes.csv', &
        '../../SRP/boxwing_modes_DRO-2_init/boxwing_modes_meta.json',epoch,status,message)
    call require(status==0,'P1 table')
    options%parameterization=ADS_PARAM_SRP_MODES
    state=[300000.0_DP,10000.0_DP,20000.0_DP,0.1_DP,0.2_DP,0.3_DP]; cov=0
    before=active_da_count()
    call hfem_ads_history_init(history,state,cov,epoch,options,status,message)
    call require(status==0,'8D history')
    call variable%init_var(8)
    call da_add(history%domain%patches(1)%da_vec%elements(1),variable,history%domain%patches(1)%da_vec%elements(1))
    call variable%destroy()
    call patch_split(history%domain%patches(1),8,left,right)
    call mf_destroy(history%domain); call mf_init(history%domain)
    call mf_push(history%domain,left); call mf_push(history%domain,right)
    call patch_destroy(left); call patch_destroy(right)
    after=active_da_count()
    call system_clock(clock_value)
    write(prefix,'(A,I0)') '/tmp/pod_srp_modes_snapshot_',clock_value
    call write_ads_patch_snapshot(history,trim(prefix),.false.,status,message,schema_version=2, &
        mode_csv_sha256=repeat('a',64),mode_metadata_sha256=repeat('b',64))
    call require(status==0,'v2 snapshot: '//trim(message))
    call require(active_da_count()==after,'snapshot has no DA ownership side effects')
    open(newunit=unit,file=trim(prefix)//'_meta.csv',status='old')
    read(unit,'(A)') line; close(unit)
    call require(index(line,'center8')>0 .and. index(line,'width8')>0,'8D metadata header')
    open(newunit=unit,file=trim(prefix)//'_coeff.csv',status='old')
    read(unit,'(A)') line; close(unit)
    call require(index(line,'power8')>0,'8D coefficient header')
    call hfem_ads_history_destroy(history)
    call require(active_da_count()==before,'snapshot and split cleanup')
    do i=1,4
        open(newunit=unit,file=trim(prefix)//trim(suffixes(i)),status='old'); close(unit,status='delete')
    end do
    call clear_srp_mode_table()
    write(*,'(A)') 'PASS: v2 eight-dimensional split snapshots and handle ownership'
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
