program test_gmm_likelihood_diagnostics
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_spice, only: str2et
    use pod_uq_gmm_state_module, only: uq_gmm_state_type
    use pod_data_format_module, only: write_gmm_measurement_diag_lines
    implicit none

    type(uq_gmm_state_type) :: gmm
    real(DP) :: et, ref_state(6)
    character(len=*), parameter :: tmpfile = "test_gmm_like_diag_tmp"
    character(len=1024) :: line, data_line
    integer :: u, ios, comp_read, cov_info_read
    logical :: found_loglike, found_mahal, found_pzzcond
    real(DP) :: loglike_noprior(2), z_mahal_sq(2), pzz_cond(2)
    real(DP) :: innovation_ra_arcsec(2), innovation_dec_arcsec(2)
    real(DP) :: det_pzz(2), logweight_prior(2), logweight_posterior(2)
    character(len=32) :: utc_read, phase_read
    real(DP) :: et_read, weight_read, pos_err_read(3), vel_err_read(3)
    real(DP) :: pos_rms_read, vel_rms_read, mahal_read, cond_read
    real(DP) :: lambda_min_read, lambda_max_read
    real(DP) :: z_mahal_read, loglike_read, logprior_read, logpost_read
    real(DP) :: det_pzz_read, pzz_cond_read, innov_ra_read, innov_dec_read

    call pod_engine_init("test_config.txt")
    call str2et("2027-02-24T16:00:00.000", et)

    call init_test_gmm(gmm)
    ref_state = [384400.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP]

    loglike_noprior = [-10.0_DP, -5.0_DP]
    z_mahal_sq = [2.5_DP, 1.3_DP]
    pzz_cond = [10.0_DP, 8.0_DP]
    innovation_ra_arcsec = [0.1_DP, -0.05_DP]
    innovation_dec_arcsec = [0.02_DP, 0.03_DP]
    det_pzz = [1.0e-12_DP, 2.0e-12_DP]
    logweight_prior = [log(0.3_DP), log(0.7_DP)]
    logweight_posterior = [log(0.15_DP), log(0.85_DP)]

    call write_gmm_measurement_diag_lines(tmpfile, et, "TEST", gmm, ref_state, .true., &
        loglike_noprior, z_mahal_sq, pzz_cond, &
        innovation_ra_arcsec, innovation_dec_arcsec, &
        det_pzz, logweight_prior, logweight_posterior)

    open(newunit=u, file=trim(tmpfile)//".gmm_diag", status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'FAIL: diagnostic file was not created'
        error stop 1
    end if

    read(u, '(A)', iostat=ios) line
    read(u, '(A)', iostat=ios) data_line
    close(u)

    found_loglike = index(line, 'LogLike_NoPrior') > 0
    found_mahal = index(line, 'Z_Mahal_Sq') > 0
    found_pzzcond = index(line, 'Pzz_Cond') > 0

    if (.not. found_loglike) then
        write(*,*) 'FAIL: header missing LogLike_NoPrior'
        error stop 1
    end if
    if (.not. found_mahal) then
        write(*,*) 'FAIL: header missing Z_Mahal_Sq'
        error stop 1
    end if
    if (.not. found_pzzcond) then
        write(*,*) 'FAIL: header missing Pzz_Cond'
        error stop 1
    end if

    read(data_line, *, iostat=ios) utc_read, et_read, phase_read, comp_read, &
        weight_read, pos_err_read, vel_err_read, pos_rms_read, vel_rms_read, &
        mahal_read, cond_read, lambda_min_read, lambda_max_read, cov_info_read, &
        z_mahal_read, loglike_read, logprior_read, logpost_read, det_pzz_read, &
        pzz_cond_read, innov_ra_read, innov_dec_read
    if (ios /= 0) then
        write(*,*) 'FAIL: diagnostic data row column order is not parseable'
        error stop 1
    end if
    if (abs(z_mahal_read - z_mahal_sq(1)) > 1.0e-12_DP) then
        write(*,*) 'FAIL: Z_Mahal_Sq is not in the expected column'
        error stop 1
    end if
    if (abs(loglike_read - loglike_noprior(1)) > 1.0e-12_DP) then
        write(*,*) 'FAIL: LogLike_NoPrior is not in the expected column'
        error stop 1
    end if

    write(*,*) 'test_gmm_likelihood_diagnostics passed'

contains

    subroutine init_test_gmm(gmm)
        type(uq_gmm_state_type), intent(inout) :: gmm
        real(DP) :: mean1(6), mean2(6), cov(6,6)
        integer :: i

        call gmm%allocate_components(2, 6)
        mean1 = [384400.0_DP, 1000.0_DP, 5000.0_DP, 0.0_DP, 1.0_DP, 0.0_DP]
        mean2 = [384400.0_DP, -1000.0_DP, -5000.0_DP, 0.0_DP, -1.0_DP, 0.0_DP]
        cov = 0.0_DP
        do i = 1, 6
            cov(i,i) = 1.0_DP
        end do

        call gmm%components(1)%init(0.5_DP, mean1, cov)
        call gmm%components(2)%init(0.5_DP, mean2, cov)
    end subroutine init_test_gmm

end program test_gmm_likelihood_diagnostics
