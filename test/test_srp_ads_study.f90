program test_srp_ads_study
    use pod_global, only: DP
    use pod_srp_ads_history_module, only: build_study_points, assess_accuracy, &
        validate_json_opm_fields
    implicit none
    real(DP), allocatable :: points(:,:)
    integer, allocatable :: validation_ids(:)
    real(DP) :: ads(6,3), truth(6,3), max_pos, max_vel
    integer :: worst_id, status, i
    logical :: failed
    character(len=256) :: message

    call build_study_points(points,validation_ids)
    call require(size(points,1)==7 .and. size(points,2)==512, '512 seven-dimensional shape points')
    call require(size(validation_ids)==207, '207 independent truth probes')
    call require(all(points(:,1)==-1.0_DP), 'first corner is the negative corner')
    call require(all(points(:,128)==1.0_DP), 'last corner is the positive corner')
    call require(all(points(:,143)==0.0_DP), 'center point follows corners and axes')
    call require(all(validation_ids==[(status,status=1,143), &
        (status,status=144,207)]), 'validation IDs are stable shape IDs')

    ads=0.0_DP
    truth=0.0_DP
    call assess_accuracy(ads,truth,0.1_DP,1.0e-6_DP, &
        max_pos,max_vel,worst_id,failed)
    call require(.not.failed .and. max_pos==0.0_DP, 'zero mismatch passes')
    ads(1,3)=0.11_DP
    call assess_accuracy(ads,truth,0.1_DP,1.0e-6_DP, &
        max_pos,max_vel,worst_id,failed)
    call require(failed .and. worst_id==3, 'position threshold identifies failing probe')
    call require(abs(max_pos-0.11_DP)<1.0e-12_DP, 'maximum position mismatch is recorded')
    ads=0.0_DP
    ads(4,2)=1.1e-6_DP
    call assess_accuracy(ads,truth,0.1_DP,1.0e-6_DP, &
        max_pos,max_vel,worst_id,failed)
    call require(failed .and. worst_id==2, 'velocity threshold identifies failing probe')

    call validate_json_opm_fields('input/DROb_20251210_9.opm',status,message)
    call require(status==0,'existing JSON OPM has required state and covariance')
    write(*,'(A)') 'PASS: SRP ADS fixed grid, validation rule, and OPM checks'
contains
    subroutine require(ok,label)
        logical,intent(in) :: ok
        character(len=*),intent(in) :: label
        if(ok) return
        write(*,'(A)') 'FAIL: '//trim(label)
        stop 1
    end subroutine require
end program test_srp_ads_study
