program test_error_line_diagnostics
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_spice, only: str2et
    use pod_data_format_module, only: write_error_line
    implicit none

    character(len=*), parameter :: prefix = "/tmp/pod_error_line_diagnostics"
    character(len=512) :: header, data_line
    real(DP) :: et
    real(DP) :: pos_err(3), vel_err(3)
    integer :: unit, ios

    call pod_engine_init("config/config.txt")
    call str2et("2027-01-21T07:12:00.000", et)

    pos_err = [1.0_DP, 2.0_DP, 3.0_DP]
    vel_err = [0.1_DP, 0.2_DP, 0.3_DP]

    call write_error_line(prefix, et, pos_err, vel_err, &
                          3.7416573867739413_DP, 0.37416573867739417_DP, &
                          3.0_DP, .true., &
                          posterior_mahalanobis_d=4.5_DP, &
                          prior_cond_p=10.0_DP, &
                          prior_lambda_min_p=1.0e-12_DP, &
                          prior_lambda_max_p=0.10_DP, &
                          posterior_cond_p=20.0_DP, &
                          posterior_lambda_min_p=2.0e-12_DP, &
                          posterior_lambda_max_p=0.40_DP)

    open(newunit=unit, file=prefix // ".err", status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "failed to open diagnostic .err output"
    read(unit, '(A)') header
    read(unit, '(A)') data_line
    close(unit)

    if (index(header, "Prior_Mahal_D") == 0) error stop "missing Prior_Mahal_D column"
    if (index(header, "Posterior_Mahal_D") == 0) error stop "missing Posterior_Mahal_D column"
    if (index(header, "Prior_CondP") == 0) error stop "missing Prior_CondP column"
    if (index(header, "Prior_LambdaMinP") == 0) error stop "missing Prior_LambdaMinP column"
    if (index(header, "Prior_LambdaMaxP") == 0) error stop "missing Prior_LambdaMaxP column"
    if (index(header, "Posterior_CondP") == 0) error stop "missing Posterior_CondP column"
    if (index(header, "Posterior_LambdaMinP") == 0) error stop "missing Posterior_LambdaMinP column"
    if (index(header, "Posterior_LambdaMaxP") == 0) error stop "missing Posterior_LambdaMaxP column"

    if (index(data_line, "3.000000") == 0) error stop "missing prior Mahalanobis value"
    if (index(data_line, "4.500000") == 0) error stop "missing posterior Mahalanobis value"
    if (index(data_line, "10.000000") == 0) error stop "missing prior condition value"
    if (index(data_line, "20.000000") == 0) error stop "missing posterior condition value"
    if (index(data_line, "1.000000E-12") == 0) error stop "prior lambda min rounded to zero"
    if (index(data_line, "2.000000E-12") == 0) error stop "posterior lambda min rounded to zero"

    write(*,*) "test_error_line_diagnostics passed"
end program test_error_line_diagnostics
