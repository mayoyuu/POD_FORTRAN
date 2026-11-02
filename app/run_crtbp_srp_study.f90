!> Single-case SRP sensitivity study over one CRTBP orbit+phase.
!!
!! The external batch (Python) reads crtbp_initial_conditions.json, enumerates
!! orbits/phases/models/domains, and calls this executable once per case in
!! parallel. This program runs the time-advancing ADS propagation for one
!! parameter box, computes endpoint mean/covariance from a fixed sample set,
!! and (when beta=0 lies in the box) the d(x)/d(beta) sensitivity at beta=0.
!!
!! Arguments (fixed order):
!!   1-6  x0(6) = x y z vx vy vz   (nondimensional)
!!   7    mu                       (CRTBP mass ratio)
!!   8    period_nd                (orbit period, nondimensional)
!!   9    alpha0_rad               (solar initial phase)
!!   10   model                    ("cannonball" | "flatplate")
!!   11   beta_lower               (m2/kg)
!!   12   beta_upper               (m2/kg)
!!   13   tilt_max_deg             (0 => u,v inactive)
!!   14   chi_lower
!!   15   chi_upper                (==chi_lower => chi inactive)
!!   16   da_order                 (e.g. 6)
!!   17   err_toll                 (truncation-error threshold, ND)
!!   18   n_split_max              (max split depth)
!!   19   output_json_path
program run_crtbp_srp_study
    use pod_global, only: DP
    use, intrinsic :: iso_fortran_env, only: int64
    use pod_crtbp_module, only: set_crtbp_mu
    use pod_srp_study_parameters_module
    use pod_crtbp_srp_study_module, only: CrtbpSrpStudyConfig, &
        STUDY_CANNONBALL, STUDY_FLAT_PLATE
    use pod_crtbp_srp_study_ads_module, only: study_ads_propagate, study_ads_log
    use pod_crtbp_srp_study_stats_module, only: study_endpoint_statistics, &
        study_beta_zero_sensitivity
    use pod_ads_split_module, only: manifold_type, mf_destroy
    implicit none

    real(DP), parameter :: TU_S = 375237.080249162_DP   ! Earth-Moon time unit (s)
    real(DP), parameter :: REL_TOL = 1.0e-10_DP, ABS_TOL = 1.0e-12_DP
    real(DP), parameter :: DT_MIN = 1.0e-13_DP
    integer,  parameter :: MAX_STEPS = 100000
    integer,  parameter :: N_SAMPLES = 10000, SEED = 20261102
    integer,  parameter :: N_SEGMENTS = 48
    real(DP), parameter :: BETA_EPS = 1.0e-4_DP

    real(DP) :: x0(6), mu, period_nd, alpha0_rad
    real(DP) :: beta_lower, beta_upper, tilt_max_deg, chi_lower, chi_upper
    real(DP) :: err_toll, t_end, dt_max
    integer  :: da_order, n_split_max
    character(len=64) :: model
    character(len=1024) :: out_path
    type(CrtbpSrpStudyConfig) :: cfg
    type(SrpStudyDomain) :: domain
    type(manifold_type) :: result_manifold
    type(study_ads_log) :: log
    integer :: index_map(4), nvar, status
    real(DP) :: err_toll_vec(6), mean_state(6), cov_state(6,6), beta_sens(6)
    real(DP) :: wall_t0, wall_t1
    integer(int64) :: count0, count1, count_rate
    real(DP) :: n_split_dir(4)
    integer :: i

    call parse_args(x0, mu, period_nd, alpha0_rad, model, beta_lower, beta_upper, &
        tilt_max_deg, chi_lower, chi_upper, da_order, err_toll, n_split_max, out_path)

    call set_crtbp_mu(mu)

    cfg%sun%alpha0_rad = alpha0_rad
    cfg%sun%tu_s = TU_S
    if (trim(model) == 'cannonball') then
        cfg%model = STUDY_CANNONBALL
    else if (trim(model) == 'flatplate') then
        cfg%model = STUDY_FLAT_PLATE
    else
        error stop 'run_crtbp_srp_study: model must be cannonball or flatplate'
    end if

    if (tilt_max_deg > 0.0_DP .and. chi_upper > chi_lower) then
        domain = srp_study_domain(beta_lower, beta_upper, tilt_max_deg=tilt_max_deg, &
            chi_lower=chi_lower, chi_upper=chi_upper)
    else if (tilt_max_deg > 0.0_DP) then
        domain = srp_study_domain(beta_lower, beta_upper, tilt_max_deg=tilt_max_deg)
    else if (chi_upper > chi_lower) then
        domain = srp_study_domain(beta_lower, beta_upper, chi_lower=chi_lower, chi_upper=chi_upper)
    else
        domain = srp_study_domain(beta_lower, beta_upper)
    end if

    t_end = 2.0_DP*period_nd
    if (t_end*TU_S < 15.0_DP*86400.0_DP) t_end = 3.0_DP*period_nd
    dt_max = t_end

    err_toll_vec = err_toll

    call system_clock(count0, count_rate)
    wall_t0 = real(count0, DP) / real(count_rate, DP)

    call study_ads_propagate(x0, cfg, domain, index_map, da_order, err_toll_vec, &
        n_split_max, t_end, N_SEGMENTS, REL_TOL, ABS_TOL, DT_MIN, dt_max, MAX_STEPS, &
        result_manifold, log, .false.)

    mean_state = 0.0_DP
    cov_state = 0.0_DP
    status = 0
    call study_endpoint_statistics(result_manifold, N_SAMPLES, SEED, mean_state, cov_state, status)

    beta_sens = 0.0_DP
    if (beta_lower == 0.0_DP) then
        call study_beta_zero_sensitivity(x0, cfg, 0.6_DP, BETA_EPS, t_end, &
            REL_TOL, ABS_TOL, DT_MIN, dt_max, MAX_STEPS, beta_sens)
    end if

    call system_clock(count1)
    wall_t1 = real(count1, DP) / real(count_rate, DP)

    n_split_dir = real(log%split_count_dir, DP)

    call write_json(out_path, x0, mu, period_nd, alpha0_rad, model, &
        beta_lower, beta_upper, tilt_max_deg, chi_lower, chi_upper, da_order, err_toll, &
        n_split_max, t_end, t_end*TU_S/86400.0_DP, result_manifold%n_patches, &
        log%n_splits_total, n_split_dir, log%n_max_active, log%first_split_time, &
        mean_state, cov_state, beta_sens, beta_lower == 0.0_DP, wall_t1 - wall_t0, status)

    call mf_destroy(result_manifold)

contains

    subroutine parse_args(x0, mu, period_nd, alpha0_rad, model, beta_lower, beta_upper, &
            tilt_max_deg, chi_lower, chi_upper, da_order, err_toll, n_split_max, out_path)
        real(DP), intent(out) :: x0(6), mu, period_nd, alpha0_rad
        real(DP), intent(out) :: beta_lower, beta_upper, tilt_max_deg, chi_lower, chi_upper
        real(DP), intent(out) :: err_toll
        integer, intent(out) :: da_order, n_split_max
        character(len=*), intent(out) :: model, out_path
        integer :: nargs, i
        character(len=1024) :: arg
        nargs = command_argument_count()
        if (nargs < 19) error stop 'run_crtbp_srp_study: expected 19 arguments'
        do i = 1, 6
            call get_command_argument(i, arg)
            read(arg, *) x0(i)
        end do
        call get_command_argument(7, arg);  read(arg, *) mu
        call get_command_argument(8, arg);  read(arg, *) period_nd
        call get_command_argument(9, arg);  read(arg, *) alpha0_rad
        call get_command_argument(10, model)
        call get_command_argument(11, arg); read(arg, *) beta_lower
        call get_command_argument(12, arg); read(arg, *) beta_upper
        call get_command_argument(13, arg); read(arg, *) tilt_max_deg
        call get_command_argument(14, arg); read(arg, *) chi_lower
        call get_command_argument(15, arg); read(arg, *) chi_upper
        call get_command_argument(16, arg); read(arg, *) da_order
        call get_command_argument(17, arg); read(arg, *) err_toll
        call get_command_argument(18, arg); read(arg, *) n_split_max
        call get_command_argument(19, out_path)
    end subroutine parse_args

    subroutine write_json(path, x0, mu, period_nd, alpha0_rad, model, &
            beta_lower, beta_upper, tilt_max_deg, chi_lower, chi_upper, da_order, err_toll, &
            n_split_max, t_end_nd, t_end_days, n_patches, n_splits, split_dir, n_max_active, &
            first_split_time, mean_state, cov_state, beta_sens, has_beta_sens, wall_s, status)
        character(len=*), intent(in) :: path, model
        real(DP), intent(in) :: x0(6), mu, period_nd, alpha0_rad
        real(DP), intent(in) :: beta_lower, beta_upper, tilt_max_deg, chi_lower, chi_upper
        real(DP), intent(in) :: err_toll, t_end_nd, t_end_days, first_split_time, wall_s
        integer, intent(in) :: da_order, n_split_max, n_patches, n_splits, n_max_active, status
        real(DP), intent(in) :: split_dir(4), mean_state(6), cov_state(6,6), beta_sens(6)
        logical, intent(in) :: has_beta_sens
        integer :: u, i, j

        open(newunit=u, file=trim(path), status='replace', action='write')
        write(u, '(A)') '{'
        write(u, '(A)') '  "model": "'//trim(model)//'",'
        write(u, '(A,ES20.12,A)') '  "mu": ', mu, ','
        write(u, '(A,ES20.12,A)') '  "period_nd": ', period_nd, ','
        write(u, '(A,ES20.12,A)') '  "alpha0_rad": ', alpha0_rad, ','
        write(u, '(A,ES20.12,A)') '  "beta_lower": ', beta_lower, ','
        write(u, '(A,ES20.12,A)') '  "beta_upper": ', beta_upper, ','
        write(u, '(A,ES20.12,A)') '  "tilt_max_deg": ', tilt_max_deg, ','
        write(u, '(A,ES20.12,A)') '  "chi_lower": ', chi_lower, ','
        write(u, '(A,ES20.12,A)') '  "chi_upper": ', chi_upper, ','
        write(u, '(A,I0,A)') '  "da_order": ', da_order, ','
        write(u, '(A,ES20.12,A)') '  "err_toll": ', err_toll, ','
        write(u, '(A,I0,A)') '  "n_split_max": ', n_split_max, ','
        write(u, '(A,ES20.12,A)') '  "t_end_nd": ', t_end_nd, ','
        write(u, '(A,ES20.12,A)') '  "t_end_days": ', t_end_days, ','
        write(u, '(A,I0,A)') '  "n_patches": ', n_patches, ','
        write(u, '(A,I0,A)') '  "n_splits_total": ', n_splits, ','
        write(u, '(A,I0,A)') '  "n_max_active": ', n_max_active, ','
        write(u, '(A,ES20.12,A)') '  "first_split_time": ', first_split_time, ','
        write(u, '(A,ES20.12,A)') '  "wall_time_s": ', wall_s, ','
        write(u, '(A,I0,A)') '  "status": ', status, ','
        write(u, '(A)', advance='no') '  "x0": ['
        call write_real_array(u, x0, 6)
        write(u, '(A)') '],'
        write(u, '(A)', advance='no') '  "split_count_dir": ['
        do i = 1, 4
            write(u, '(I0)', advance='no') nint(split_dir(i))
            if (i < 4) write(u, '(A)', advance='no') ', '
        end do
        write(u, '(A)') '],'
        write(u, '(A)', advance='no') '  "mean_state": ['
        call write_real_array(u, mean_state, 6)
        write(u, '(A)') '],'
        write(u, '(A)', advance='no') '  "cov_state": ['
        do i = 1, 6
            do j = 1, 6
                write(u, '(ES20.12)', advance='no') cov_state(i, j)
                if (.not. (i == 6 .and. j == 6)) write(u, '(A)', advance='no') ', '
            end do
        end do
        write(u, '(A)') '],'
        if (has_beta_sens) then
            write(u, '(A)', advance='no') '  "beta_sensitivity": ['
            call write_real_array(u, beta_sens, 6)
            write(u, '(A)') ']'
        else
            write(u, '(A)') '  "beta_sensitivity": null'
        end if
        write(u, '(A)') '}'
        close(u)
    end subroutine write_json

    subroutine write_real_array(u, a, n)
        integer, intent(in) :: u, n
        real(DP), intent(in) :: a(n)
        integer :: i
        do i = 1, n
            write(u, '(ES20.12)', advance='no') a(i)
            if (i < n) write(u, '(A)', advance='no') ', '
        end do
    end subroutine write_real_array

end program run_crtbp_srp_study
