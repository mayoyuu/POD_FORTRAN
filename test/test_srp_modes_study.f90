program test_srp_modes_study
    use pod_global, only: DP
    use pod_config, only: config,load_config
    use pod_spice, only: spice_init
    use pod_data_format_module, only: load_initial_opm
    use pod_srp_mode_module
    use pod_srp_ads_history_module, only: build_study_points,advance_real_mode_probes
    use pod_uq_hfem_ads_module, only: hfem_ads_options_type
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
    implicit none
    real(DP),allocatable :: points(:,:),legacy(:,:),explicit7(:,:)
    integer,allocatable :: ids(:),oldids(:),explicitids(:)
    real(DP) :: states(6,4),initial(6),cov(6,6),epoch,s(4),u(4)
    type(hfem_ads_options_type) :: options
    character(len=256) :: message
    integer :: status,j,k
    call build_study_points(legacy,oldids)
    call build_study_points(explicit7,explicitids,7)
    call require(all(legacy==explicit7) .and. all(oldids==explicitids),'legacy grid unchanged')
    call build_study_points(points,ids,8)
    call require(size(points,1)==8 .and. size(points,2)==512 .and. size(ids)==337,'8D grid sizes')
    do j=0,255
        do k=1,8
            call require(points(k,j+1)==merge(1.0_DP,-1.0_DP,btest(j,k-1)),'all 256 corners')
        end do
    end do
    call require(all(points(:,273)==0),'8D center')
    call require(all(abs(points)<=1),'finite root box')
    call build_study_points(points,ids,6)
    call require(size(points,1)==6,'6D grid')
    call spice_init(); call load_config('config/config_box_wing_earth.txt')
    call load_initial_opm('OPM/DRO-2/DRO-2_init.opm.json',epoch,initial,cov)
    config%srp_model='time_varying_modes'
    call load_srp_mode_table('../../SRP/boxwing_modes_DRO-2_init/boxwing_modes.csv', &
        '../../SRP/boxwing_modes_DRO-2_init/boxwing_modes_meta.json',epoch,status,message)
    call require(status==0,'P1 table')
    states=spread(initial,2,4); s=[0.0_DP,1.0_DP,0.0_DP,0.5_DP]; u=[0.0_DP,0.2_DP,1.0_DP,0.5_DP]
    call advance_real_mode_probes(states,s,u,epoch,0.0_DP,30.0_DP,options,status,message)
    call require(status==0,'independent probes advance: '//trim(message))
    call require(norm2(states(:,2)-states(:,1))>1e-10_DP,'s affects Real trajectory')
    call require(norm2(states(:,3)-states(:,1))>1e-10_DP,'u affects Real trajectory')
    call advance_real_mode_probes(states,s(:3),u,epoch,30.0_DP,60.0_DP,options,status,message)
    call require(status/=0,'shape mismatch')
    u(4)=ieee_value(0.0_DP,ieee_quiet_nan)
    call advance_real_mode_probes(states,s,u,epoch,30.0_DP,60.0_DP,options,status,message)
    call require(status/=0,'NaN rejected without error stop')
    u(4)=1.1_DP
    call advance_real_mode_probes(states,s,u,epoch,30.0_DP,60.0_DP,options,status,message)
    call require(status/=0,'outside coefficients')
    u(4)=0.5_DP
    call advance_real_mode_probes(states,s,u,epoch,60.0_DP,30.0_DP,options,status,message)
    call require(status/=0,'backwards time')
    call advance_real_mode_probes(states,s,u,epoch,30.0_DP,1296001.0_DP,options,status,message)
    call require(status/=0,'coverage preflight')
    call clear_srp_mode_table()
    write(*,'(A)') 'PASS: 6/7/8D study points and independent Real mode probes'
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
