program test_hfem_ads_modes_8d
    use pod_global, only: DP
    use pod_config, only: config
    use pod_engine_module, only: pod_engine_init
    use pod_data_format_module, only: load_initial_opm
    use pod_srp_mode_module
    use pod_uq_ads_coordinates_module, only: ADS_PARAM_SRP_MODES
    use pod_uq_hfem_ads_module
    use pod_dace_classes, only: DA,active_da_count,da_exp_sub,da_mul,da_add,assignment(=)
    use pod_force_model_module, only: set_propagation_epoch,set_srp_mode_coefficients_real,clear_srp_mode_coefficients_real
    use pod_integrator_module, only: adaptive_step_integrate,METHOD_RKF78
    implicit none
    type(hfem_ads_history_type) :: history
    type(hfem_ads_options_type) :: options
    type(hfem_ads_stats_type) :: stats
    type(DA) :: variable,cubic
    real(DP) :: state(6),cov(6,6),epoch,q(8,4),input(8,4),truth(6),nondim(6),values(6,4),saved(6,4)
    real(DP),allocatable :: output(:,:),times(:),integrated(:,:)
    logical :: found(4)
    integer :: status,before,j,i,n_steps,direction,live_handles,saved_max_steps
    character(len=256) :: message
    call pod_engine_init('config/config_box_wing_earth.txt')
    call load_initial_opm('OPM/DRO-2/DRO-2_init.opm.json',epoch,state,cov)
    options%parameterization=ADS_PARAM_SRP_MODES
    options%da_order=4
    before=active_da_count()
    config%srp_model='cannonball'
    call hfem_ads_history_init(history,state,cov,epoch,options,status,message)
    call require(status/=0 .and. .not.history%initialized,'reject modes with wrong force model')
    config%srp_model='time_varying_modes'
    call hfem_ads_history_init(history,state,cov,epoch,options,status,message)
    call require(status/=0 .and. active_da_count()==before,'reject unloaded modes without handles')
    call load_srp_mode_table('../../SRP/boxwing_modes_DRO-2_init/boxwing_modes.csv', &
        '../../SRP/boxwing_modes_DRO-2_init/boxwing_modes_meta.json',epoch,status,message)
    call require(status==0,'load P1 table')
    call hfem_ads_history_init(history,state,cov,epoch,options,status,message)
    call require(status==0 .and. history%da_variables==8,'8D history init')
    call require(size(history%stats%first_split_direction_scores)==8,'8D scores')
    q=0; q(7,:)=[-1.0_DP,1.0_DP,-1.0_DP,0.4_DP]; q(8,:)=[-1.0_DP,0.3_DP,1.0_DP,-0.6_DP]
    call hfem_ads_history_evaluate(history,q,values,found,status)
    call require(status==0 .and. all(found),'8D root evaluation')
    call require(maxval(abs(values-spread(state,2,4)))<1e-9_DP,'initial state mapping')
    input(1:6,:)=values; input(7:8,:)=0.5_DP*(1+q(7:8,:))
    ! Force more than six saved DA nodes: DA states are (6,n_steps), not (n_steps,6).
    options%dt_max=1.0_DP
    history%options%dt_max=options%dt_max
    call hfem_ads_history_advance(history,30.0_DP,status,message)
    call require(status==0,'8D history short advance: '//trim(message))
    call hfem_ads_history_evaluate(history,q,values,found,status)
    do j=1,4
        call set_propagation_epoch(epoch)
        call set_srp_mode_coefficients_real(input(7,j),input(8,j))
        nondim(1:3)=state(1:3)/config%LU; nondim(4:6)=state(4:6)/config%VU
        call adaptive_step_integrate(nondim,0.0_DP,30.0_DP/config%TU,METHOD_RKF78, &
            times,integrated,n_steps,options%max_steps,options%rel_tol,options%abs_tol,options%dt_min,options%dt_max)
        truth(1:3)=integrated(n_steps,1:3)*config%LU; truth(4:6)=integrated(n_steps,4:6)*config%VU
        call require(norm2(values(1:3,j)-truth(1:3))<0.1_DP,'history matches Real position')
        call require(norm2(values(4:6,j)-truth(4:6))<1e-6_DP,'history matches Real velocity')
    end do
    call clear_srp_mode_coefficients_real()
    saved=values; live_handles=active_da_count()
    history%options%max_steps=2
    call hfem_ads_history_advance(history,60.0_DP,status,message)
    call require(status/=0 .and. history%current_time==30.0_DP,'incomplete history leaves old checkpoint intact')
    call require(active_da_count()==live_handles,'failed owned history integration cleanup')
    call hfem_ads_history_evaluate(history,q,values,found,status)
    call require(status==0 .and. maxval(abs(values-saved))==0,'failed history preserves old polynomial')
    call hfem_ads_history_destroy(history)
    call require(active_da_count()==before,'history cleanup')
    options%dt_max=0.1_DP
    call hfem_ads_propagate(state,cov,epoch,0.0_DP,1.0_DP,input,options,output,stats,status,message)
    call require(status==0 .and. size(output,1)==8 .and. size(output,2)==4,'8D one-shot propagation')
    call require(maxval(abs(output(7:8,:)-input(7:8,:)))==0,'preserve both physical parameters')
    call require(active_da_count()==before,'one-shot cleanup')
    saved_max_steps=options%max_steps; options%max_steps=2
    call hfem_ads_propagate(state,cov,epoch,0.0_DP,1.0_DP,input,options,output,stats,status,message)
    call require(status/=0 .and. active_da_count()==before,'failed owned one-shot integration cleanup')
    options%max_steps=saved_max_steps
    ! A synthetic inherited position polynomial forces each parameter-direction split
    ! through the real history BFS, with SRP off to isolate event/translation mechanics.
    config%use_srp=.false.
    options%dt_max=3600.0_DP
    cov=0.0_DP ! Isolate inherited parameter curvature from the six physical-state directions.
    options%max_split_depth=1; options%error_tolerance=1e-9_DP
    options%da_order=2
    do direction=7,8
        call hfem_ads_history_init(history,state,cov,epoch,options,status,message)
        call variable%init_var(direction)
        call da_exp_sub(variable,cubic)
        call da_mul(cubic,1.0_DP,cubic)
        call da_add(history%domain%patches(1)%da_vec%elements(1),cubic,history%domain%patches(1)%da_vec%elements(1))
        call variable%destroy(); call cubic%destroy()
        call hfem_ads_history_advance(history,1.0_DP,status,message)
        call require(status==0 .and. history%stats%split_occurred,'actual parameter split')
        call require(history%stats%first_split_direction==direction,'split searches all eight variables')
        call require(size(history%stats%split_events(1)%parent_center)==8,'event retains eight coordinates')
        call hfem_ads_history_destroy(history)
        call require(active_da_count()==before,'split cleanup')
    end do
    call clear_srp_mode_table()
    write(*,'(A)') 'PASS: HFEM 8D propagation, parameter retention, events and ownership'
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
