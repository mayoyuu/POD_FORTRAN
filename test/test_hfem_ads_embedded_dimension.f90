program test_hfem_ads_embedded_dimension
    use pod_global, only: DP
    use pod_uq_ads_coordinates_module, only: ADS_COORD_WHITENED
    use pod_uq_hfem_ads_module, only: hfem_ads_options_type,hfem_ads_history_type, &
        hfem_ads_history_init,hfem_ads_history_evaluate,hfem_ads_history_destroy
    implicit none
    type(hfem_ads_options_type) :: options
    type(hfem_ads_history_type) :: history
    real(DP) :: nominal(6),covariance(6,6),points(6,2),values(6,2)
    logical :: found(2)
    character(len=256) :: message
    integer :: i,status

    nominal=[1.0_DP,2.0_DP,3.0_DP,4.0_DP,5.0_DP,6.0_DP]
    covariance=0.0_DP
    do i=1,6
        covariance(i,i)=real(i,DP)**2
    end do
    options%coordinate_mode=ADS_COORD_WHITENED
    options%domain_sigma=3.0_DP
    options%srp_sigma=0.0_DP
    options%da_storage_variables=7
    call hfem_ads_history_init(history,nominal,covariance,0.0_DP,options,status,message)
    call assert_true(status==0,'embedded init: '//trim(message))
    call assert_true(history%coordinate_map%n_variables==6,'six active variables')
    call assert_true(history%da_variables==7,'seven storage variables')
    points=0.0_DP
    points(:,2)=0.5_DP
    call hfem_ads_history_evaluate(history,points,values,found,status)
    call assert_true(status==0.and.all(found),'six-dimensional evaluation')
    call assert_true(maxval(abs(values(:,1)-nominal))<1.0e-13_DP,'center value')
    call assert_true(maxval(abs(values(:,2)-(nominal+1.5_DP*[(real(i,DP),i=1,6)])))< &
        1.0e-12_DP,'embedded affine value')
    call hfem_ads_history_destroy(history)

    write(*,'(a)') 'PASS: six-dimensional control embedded in seven-variable DACE context.'
contains
    subroutine assert_true(condition,label)
        logical,intent(in) :: condition
        character(len=*),intent(in) :: label
        if(.not.condition) then
            write(*,'(a)') 'FAIL: '//trim(label)
            error stop 1
        end if
    end subroutine assert_true
end program test_hfem_ads_embedded_dimension
