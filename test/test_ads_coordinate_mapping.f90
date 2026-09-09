program test_ads_coordinate_mapping
    use pod_global, only: DP
    use pod_uq_ads_coordinates_module, only: ADS_COORD_COMPONENT, &
        ADS_COORD_WHITENED, ads_coordinate_map_type, ads_build_coordinate_map, &
        ads_physical_to_unit, ads_unit_to_physical
    implicit none

    type(ads_coordinate_map_type) :: map
    real(DP) :: mean6(6), covariance(6,6), expected_basis(6,6)
    real(DP) :: lower(6,6), deviation(6), recovered(6), eta, recovered_eta
    real(DP), allocatable :: unit_point(:)
    character(len=256) :: message
    integer :: status, i

    mean6 = [1.0_DP,2.0_DP,3.0_DP,4.0_DP,5.0_DP,6.0_DP]
    covariance = 0.0_DP
    do i=1,6
        covariance(i,i)=real(i*i,DP)
    end do
    covariance(1,2)=0.5_DP
    covariance(2,1)=0.5_DP

    call ads_build_coordinate_map(mean6,covariance,ADS_COORD_COMPONENT, &
        3.0_DP,0.0_DP,map,status,message)
    call assert_true(status==0,'component map status')
    call assert_true(map%n_variables==6,'srp sigma zero selects 6D')
    call assert_true(map%mode==ADS_COORD_COMPONENT,'component mode')
    expected_basis=0.0_DP
    do i=1,6
        expected_basis(i,i)=3.0_DP*sqrt(covariance(i,i))
    end do
    call assert_close(maxval(abs(map%basis6-expected_basis)),0.0_DP, &
        1.0e-14_DP,'component basis')

    deviation=[0.75_DP,-1.5_DP,2.25_DP,-3.0_DP,3.75_DP,-4.5_DP]
    call ads_physical_to_unit(map,deviation,0.0_DP,unit_point,status)
    call assert_true(status==0,'component physical-to-unit status')
    call assert_true(size(unit_point)==6,'component unit dimension')
    call ads_unit_to_physical(map,unit_point,recovered,recovered_eta,status)
    call assert_true(status==0,'component unit-to-physical status')
    call assert_close(maxval(abs(recovered-deviation)),0.0_DP,1.0e-14_DP, &
        'component round trip')
    call assert_close(recovered_eta,0.0_DP,1.0e-14_DP,'disabled SRP eta')

    lower=0.0_DP
    lower(1,1)=2.0_DP
    lower(2,1)=0.4_DP; lower(2,2)=1.5_DP
    lower(3,1)=-0.2_DP; lower(3,2)=0.3_DP; lower(3,3)=1.2_DP
    lower(4,1)=0.1_DP; lower(4,4)=0.8_DP
    lower(5,2)=-0.1_DP; lower(5,4)=0.2_DP; lower(5,5)=0.6_DP
    lower(6,3)=0.05_DP; lower(6,5)=-0.08_DP; lower(6,6)=0.4_DP
    covariance=matmul(lower,transpose(lower))
    call ads_build_coordinate_map(mean6,covariance,ADS_COORD_WHITENED, &
        2.5_DP,0.02_DP,map,status,message)
    call assert_true(status==0,'whitened map status')
    call assert_true(map%n_variables==7,'positive srp sigma selects 7D')
    call assert_close(maxval(abs(matmul(map%basis6,transpose(map%basis6))- &
        6.25_DP*covariance)),0.0_DP,1.0e-12_DP,'whitened covariance basis')
    call assert_close(upper_triangle_max(map%basis6), &
        0.0_DP,1.0e-14_DP,'whitened upper triangle')

    unit_point=[0.2_DP,-0.3_DP,0.4_DP,-0.5_DP,0.6_DP,-0.7_DP,0.25_DP]
    call ads_unit_to_physical(map,unit_point,deviation,eta,status)
    call assert_true(status==0,'whitened unit-to-physical status')
    call assert_close(eta,2.5_DP*0.02_DP*0.25_DP,1.0e-14_DP, &
        'optional SRP coordinate scale')
    call ads_physical_to_unit(map,deviation,eta,unit_point,status)
    call assert_true(status==0,'whitened physical-to-unit status')
    call assert_close(maxval(abs(unit_point- &
        [0.2_DP,-0.3_DP,0.4_DP,-0.5_DP,0.6_DP,-0.7_DP,0.25_DP])), &
        0.0_DP,1.0e-13_DP,'whitened round trip')

    call ads_build_coordinate_map(mean6,covariance,ADS_COORD_COMPONENT, &
        0.0_DP,0.0_DP,map,status,message)
    call assert_true(status<0,'zero domain sigma rejected')
    call ads_build_coordinate_map(mean6,covariance,ADS_COORD_COMPONENT, &
        3.0_DP,-0.1_DP,map,status,message)
    call assert_true(status<0,'negative SRP sigma rejected')

    covariance=0.0_DP
    do i=1,5
        covariance(i,i)=1.0_DP
    end do
    covariance(6,6)=-1.0_DP
    call ads_build_coordinate_map(mean6,covariance,ADS_COORD_WHITENED, &
        3.0_DP,0.0_DP,map,status,message)
    call assert_true(status<0,'non-positive-definite covariance rejected')

    covariance=0.0_DP
    covariance(1,1)=1.0_DP
    call ads_build_coordinate_map(mean6,covariance,ADS_COORD_COMPONENT, &
        3.0_DP,0.0_DP,map,status,message)
    call assert_true(status==0,'component permits zero variance')
    deviation=0.0_DP
    deviation(2)=1.0_DP
    call ads_physical_to_unit(map,deviation,0.0_DP,unit_point,status)
    call assert_true(status<0,'nonzero deviation on zero scale rejected')
    deviation(2)=0.0_DP
    call ads_physical_to_unit(map,deviation,0.0_DP,unit_point,status)
    call assert_true(status==0,'zero deviation on zero scale accepted')

    write(*,'(a)') 'PASS: ADS coordinate mappings and optional SRP dimension.'

contains

    pure real(DP) function upper_triangle_max(matrix) result(value)
        real(DP), intent(in) :: matrix(6,6)
        integer :: row, column

        value = 0.0_DP
        do column = 2, 6
            do row = 1, column - 1
                value = max(value, abs(matrix(row,column)))
            end do
        end do
    end function upper_triangle_max

    subroutine assert_true(condition,label)
        logical,intent(in) :: condition
        character(len=*),intent(in) :: label
        if(.not.condition) then
            write(*,'(a)') 'FAIL: '//trim(label)
            error stop 1
        end if
    end subroutine assert_true

    subroutine assert_close(actual,expected,tolerance,label)
        real(DP),intent(in) :: actual,expected,tolerance
        character(len=*),intent(in) :: label
        if(abs(actual-expected)>tolerance) then
            write(*,'(a,3(1x,es24.16))') 'FAIL: '//trim(label), &
                actual,expected,tolerance
            error stop 1
        end if
    end subroutine assert_close

end program test_ads_coordinate_mapping
