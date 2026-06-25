program test_gmm_component_diagnostics
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_spice, only: str2et
    use pod_uq_gmm_state_module, only: uq_gmm_state_type
    use pod_data_format_module, only: write_gmm_diag_lines
    implicit none

    character(len=*), parameter :: prefix = "/tmp/pod_gmm_component_diagnostics"
    type(uq_gmm_state_type) :: gmm
    real(DP) :: et
    real(DP) :: ref_state(6), cov1(6,6), cov2(6,6)
    character(len=1024) :: header, line1, line2
    integer :: unit, ios, i

    call pod_engine_init("config/config.txt")
    call str2et("2027-02-24T16:00:59.904", et)

    ref_state = [0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP]
    cov1 = 0.0_DP
    cov2 = 0.0_DP
    do i = 1, 6
        cov1(i,i) = 1.0_DP
        cov2(i,i) = 2.0_DP
    end do

    call gmm%allocate_components(2, 6)
    call gmm%components(1)%init(0.25_DP, &
        [1.0_DP, 2.0_DP, 3.0_DP, 0.1_DP, 0.2_DP, 0.3_DP], cov1)
    call gmm%components(2)%init(0.75_DP, &
        [4.0_DP, 5.0_DP, 6.0_DP, 0.4_DP, 0.5_DP, 0.6_DP], cov2)

    call write_gmm_diag_lines(prefix, et, "POSTERIOR", gmm, ref_state, .true.)

    open(newunit=unit, file=prefix // ".gmm_diag", status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "failed to open .gmm_diag output"
    read(unit, '(A)') header
    read(unit, '(A)') line1
    read(unit, '(A)') line2
    close(unit)

    if (index(header, "Phase") == 0) error stop "missing Phase column"
    if (index(header, "Comp") == 0) error stop "missing Comp column"
    if (index(header, "Weight") == 0) error stop "missing Weight column"
    if (index(header, "Mean_Pos_RMS") == 0) error stop "missing Mean_Pos_RMS column"
    if (index(header, "CondP") == 0) error stop "missing CondP column"
    if (index(header, "LambdaMinP") == 0) error stop "missing LambdaMinP column"

    if (index(line1, "POSTERIOR") == 0) error stop "missing posterior phase"
    if (index(line1, " 1 ") == 0) error stop "missing component 1 row"
    if (index(line1, "2.5000000000E-001") == 0) error stop "missing component 1 weight"
    if (index(line2, " 2 ") == 0) error stop "missing component 2 row"
    if (index(line2, "7.5000000000E-001") == 0) error stop "missing component 2 weight"

    write(*,*) "test_gmm_component_diagnostics passed"
end program test_gmm_component_diagnostics