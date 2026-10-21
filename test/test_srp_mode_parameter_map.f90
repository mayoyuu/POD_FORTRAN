program test_srp_mode_parameter_map
    use pod_global, only: DP
    use pod_uq_ads_coordinates_module
    use pod_dace_classes, only: DA,AlgebraicVector,dace_initialize,active_da_count,da_mul,da_add,assignment(=)
    use pod_ads_split_module, only: patch_type,patch_init,patch_split,patch_destroy,sh_center,sh_width,sh_map_point
    implicit none
    type(ads_coordinate_map_type) :: map
    type(ads_physical_parameters_type) :: parameters,recovered
    type(AlgebraicVector) :: source
    type(DA) :: s,u
    type(patch_type) :: parent,left,right
    real(DP) :: mean(6),covariance(6,6),deviation(6),q(8),expected(8),s0,ds,u0,du
    real(DP),allocatable :: result(:),center(:),width(:),values(:)
    integer :: status,i,k,direction,side,before
    character(len=256) :: message
    mean=0.0_DP; covariance=0.0_DP
    parameters=ads_physical_parameters_type()
    call require(parameters%eta==0.0_DP .and. parameters%s==0.0_DP .and. parameters%u==0.0_DP,'parameter defaults')
    do i=1,6
        covariance(i,i)=real(i,DP)
    end do
    call ads_build_coordinate_map(mean,covariance,ADS_COORD_COMPONENT,3.0_DP,0.0_DP,map,status,message)
    call require(status==0 .and. map%n_variables==6 .and. map%parameterization==ADS_PARAM_NONE,'AUTO 6D')
    call ads_build_coordinate_map(mean,covariance,ADS_COORD_COMPONENT,3.0_DP,0.02_DP,map,status,message)
    call require(status==0 .and. map%n_variables==7 .and. map%parameterization==ADS_PARAM_SRP_ETA,'AUTO 7D')
    call ads_build_coordinate_map(mean,covariance,ADS_COORD_WHITENED,3.0_DP,0.0_DP,map,status,message,ADS_PARAM_SRP_MODES)
    call require(status==0 .and. map%n_variables==8,'explicit 8D')
    do k=0,256
        if(k<256) then
            do i=1,8
                q(i)=merge(1.0_DP,-1.0_DP,btest(k,i-1))
            end do
        else
            q=0.0_DP
        end if
        expected=q
        call ads_unit_to_physical(map,q,deviation,parameters,status)
        call require(status==0,'unit conversion')
        call require(parameters%s>=0 .and. parameters%s<=1 .and. parameters%u>=0 .and. parameters%u<=1,'domain')
        call require(parameters%s+(1.0_DP-parameters%s)*parameters%u<=1.0_DP,'triangle')
        call ads_physical_to_unit(map,deviation,parameters,result,status)
        call require(status==0 .and. maxval(abs(result-expected))<1e-13_DP,'corner/center round trip')
    end do
    do k=1,100
        call random_number(q); q=2*q-1
        call ads_unit_to_physical(map,q,deviation,parameters,status)
        call ads_physical_to_unit(map,deviation,parameters,result,status)
        call require(status==0 .and. maxval(abs(result-q))<1e-13_DP,'random round trip')
    end do
    q=0; q(7)=1+1e-8_DP
    call ads_unit_to_physical(map,q,deviation,parameters,status)
    call require(status/=0,'outside unit domain')
    q=0; q(8)=1+100*epsilon(1.0_DP)
    call ads_unit_to_physical(map,q,deviation,parameters,status)
    call require(status==0 .and. parameters%u==1,'roundoff clips physical boundary')
    parameters=ads_physical_parameters_type(); parameters%s=-1e-5_DP
    call ads_physical_to_unit(map,mean,parameters,result,status)
    call require(status/=0,'negative physical s')
    parameters%s=0.5_DP; parameters%eta=0.01_DP
    call ads_physical_to_unit(map,mean,parameters,result,status)
    call require(status/=0,'mixed physical eta/modes')
    call ads_unit_to_physical(map,q,deviation,s0,status)
    call require(status/=0,'legacy scalar wrapper rejects modes')
    call dace_initialize(2,8); before=active_da_count()
    call source%init(2); call s%init_var(7); call u%init_var(8)
    call da_mul(s,0.5_DP,source%elements(1)); call da_add(source%elements(1),0.5_DP,source%elements(1))
    call da_mul(u,0.5_DP,source%elements(2)); call da_add(source%elements(2),0.5_DP,source%elements(2))
    call patch_init(parent,source)
    do direction=7,8
        call patch_split(parent,direction,left,right)
        do side=1,2
            if(side==1) then
                center=sh_center(left%history); width=sh_width(left%history)
            else
                center=sh_center(right%history); width=sh_width(right%history)
            end if
            call ads_patch_mode_parameters(map,center,width,s0,ds,u0,du,status)
            call require(status==0,'patch mode coefficients')
            q=0; q(direction)=merge(-0.75_DP,0.75_DP,side==1)
            expected=q
            if(side==1) then
                call sh_map_point(left%history,q); values=left%da_vec%eval(q)
            else
                call sh_map_point(right%history,q); values=right%da_vec%eval(q)
            end if
            call require(abs(values(1)-(s0+ds*q(7)))<1e-14_DP,'split s global value')
            call require(abs(values(2)-(u0+du*q(8)))<1e-14_DP,'split u global value')
            call require(maxval(abs(values-0.5_DP*(1+expected(7:8))))<1e-14_DP,'split not renormalized')
        end do
        call patch_destroy(left); call patch_destroy(right)
    end do
    call patch_destroy(parent); call source%destroy(); call s%destroy(); call u%destroy()
    call require(active_da_count()==before,'all ADS handles released')
    call ads_build_coordinate_map(mean,covariance,ADS_COORD_COMPONENT,3.0_DP,0.02_DP,map,status,message,ADS_PARAM_SRP_MODES)
    call require(status/=0,'mixed map parameterization')
    write(*,'(A)') 'PASS: 8D modes coordinate maps and actual ADS s/u splits'
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
