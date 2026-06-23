program test_mahalanobis_prepost_diagnostics
    use pod_global, only: DP
    use pod_error_analysis_module, only: compute_orbit_error, compute_covariance_condition_number
    implicit none

    real(DP) :: ref_state(6), pred_state(6), upd_state(6)
    real(DP) :: pred_cov(6,6), upd_cov(6,6)
    real(DP) :: pos_err(3), vel_err(3), pos_rms, vel_rms
    real(DP) :: pred_mahal, upd_mahal
    real(DP) :: pred_cond, upd_cond
    real(DP) :: pred_lambda_min, pred_lambda_max
    real(DP) :: upd_lambda_min, upd_lambda_max
    integer :: pred_info, upd_info, i

    ref_state = 0.0_DP
    pred_state = 0.0_DP
    upd_state = 0.0_DP
    pred_state(1) = 2.0_DP
    upd_state(1) = 1.0_DP

    pred_cov = 0.0_DP
    upd_cov = 0.0_DP
    do i = 1, 6
        pred_cov(i,i) = 4.0_DP
        upd_cov(i,i) = 4.0_DP
    end do
    upd_cov(1,1) = 0.01_DP

    call compute_orbit_error(pred_state, pred_cov, ref_state, &
                             pos_err, vel_err, pos_rms, vel_rms, pred_mahal)
    if (abs(pos_rms - 2.0_DP) > 1.0e-12_DP) error stop "unexpected prediction position RMS"
    if (abs(pred_mahal - 1.0_DP) > 1.0e-12_DP) error stop "unexpected prediction Mahalanobis distance"

    call compute_orbit_error(upd_state, upd_cov, ref_state, &
                             pos_err, vel_err, pos_rms, vel_rms, upd_mahal)
    if (abs(pos_rms - 1.0_DP) > 1.0e-12_DP) error stop "unexpected update position RMS"
    if (abs(upd_mahal - 10.0_DP) > 1.0e-12_DP) error stop "unexpected update Mahalanobis distance"
    if (upd_mahal <= pred_mahal) error stop "update Mahalanobis should grow when covariance collapses"

    call compute_covariance_condition_number(pred_cov, pred_cond, pred_info, &
                                             lambda_min=pred_lambda_min, lambda_max=pred_lambda_max)
    call compute_covariance_condition_number(upd_cov, upd_cond, upd_info, &
                                             lambda_min=upd_lambda_min, lambda_max=upd_lambda_max)

    if (pred_info /= 0) error stop "prediction covariance should be valid"
    if (upd_info /= 0) error stop "update covariance should be valid"
    if (abs(pred_cond - 1.0_DP) > 1.0e-12_DP) error stop "unexpected prediction condition number"
    if (abs(upd_cond - 400.0_DP) > 1.0e-10_DP) error stop "unexpected update condition number"
    if (abs(pred_lambda_min - 4.0_DP) > 1.0e-12_DP) error stop "unexpected prediction lambda_min"
    if (abs(pred_lambda_max - 4.0_DP) > 1.0e-12_DP) error stop "unexpected prediction lambda_max"
    if (abs(upd_lambda_min - 0.01_DP) > 1.0e-12_DP) error stop "unexpected update lambda_min"
    if (abs(upd_lambda_max - 4.0_DP) > 1.0e-12_DP) error stop "unexpected update lambda_max"

    write(*,*) "test_mahalanobis_prepost_diagnostics passed"
end program test_mahalanobis_prepost_diagnostics
