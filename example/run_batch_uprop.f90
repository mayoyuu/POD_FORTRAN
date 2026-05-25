program run_batch_uprop
    use pod_global, only: DP
    use pod_uq_crtbp_mc_module, only: crtbp_mc_propagate
    use pod_uq_crtbp_da_module, only: crtbp_da_propagate
    use pod_dace_classes
    implicit none

    real(DP), parameter :: LU = 384400.0_DP
    real(DP), parameter :: TU = 375190.26189464843_DP
    real(DP), parameter :: MU = 0.012153614091892_DP
    real(DP), parameter :: T_END_DAYS = 4.0_DP
    real(DP), parameter :: DAY_SEC = 86400.0_DP

    character(len=32)  :: method
    character(len=512) :: input_csv, output_csv
    character(len=2048) :: line
    character(len=32)  :: center, family

    real(DP) :: nominal_state(6), cov(6,6), T_nd
    real(DP) :: final_mean(6), final_cov(6,6)
    real(DP), allocatable :: final_samples(:,:)
    real(DP) :: propagated_ref(6)

    real(DP) :: mean_dim(6), cov_dim(6,6), T_dim
    real(DP) :: sigma_pos(3), sigma_vel(3), sigma_pos_rms, sigma_vel_rms
    real(DP) :: t_end

    integer  :: n_samples, da_order, max_steps_int
    real(DP) :: rel_tol, abs_tol, dt_min, dt_max

    integer :: in_unit, out_unit, ios, n_rows, i

    if (command_argument_count() /= 2) then
        print *, 'Usage: run_batch_uprop <method> <input_csv>'
        print *, '  method   : MC or DA'
        print *, '  input_csv: path to input CSV file'
        stop
    end if

    call get_command_argument(1, method)
    call get_command_argument(2, input_csv)

    if (trim(method) /= 'MC' .and. trim(method) /= 'DA') then
        print *, 'Error: method must be "MC" or "DA", got: ', trim(method)
        stop
    end if

    t_end = DAY_SEC * T_END_DAYS / TU

    n_samples     = 10000
    da_order      = 4
    rel_tol       = 1.0e-12_DP
    abs_tol       = 1.0e-12_DP
    dt_min        = 1.0e-10_DP
    dt_max        = 0.1_DP
    max_steps_int = 100000

    ! ---- Open input ----
    open(newunit=in_unit, file=trim(input_csv), status='old', action='read')
    read(in_unit, '(A)') line  ! skip header

    ! ---- Open output ----
    call make_output_filename(input_csv, method, output_csv)
    open(newunit=out_unit, file=trim(output_csv), status='replace', action='write')
    write(out_unit, '(A)') &
        'Center,Family,T(day),x_mean(km),y_mean(km),z_mean(km),' // &
        'vx_mean(km/s),vy_mean(km/s),vz_mean(km/s),' // &
        'sigma_x(km),sigma_y(km),sigma_z(km),sigma_pos_RMS(km),' // &
        'sigma_vx(km/s),sigma_vy(km/s),sigma_vz(km/s),sigma_vel_RMS(km/s)'

    if (trim(method) == 'DA') call dace_initialize(da_order, 6)

    ! =========================================================
    !  Main loop over input rows
    ! =========================================================
    n_rows = 0
    do
        read(in_unit, '(A)', iostat=ios) line
        if (ios /= 0) exit

        call parse_csv_line(line, center, family, nominal_state, T_nd)
        n_rows = n_rows + 1

        ! Dimensional initial covariance (same as run_CRTBP_uprop)
        cov = 0.0_DP
        cov(1,1) = 3.0_DP
        cov(2,2) = 3.0_DP
        cov(3,3) = 3.0_DP
        cov(4,4) = 3.33e-5_DP*1e-6_DP
        cov(5,5) = 3.33e-5_DP*1e-6_DP
        cov(6,6) = 3.33e-5_DP*1e-6_DP

        call nondim_cov(cov, LU, TU)

        if (trim(method) == 'MC') then
            call crtbp_mc_propagate(nominal_state, cov, MU, t_end, &
                n_samples, rel_tol, abs_tol, dt_min, dt_max, max_steps_int, &
                final_samples, final_mean, final_cov, .false.)
        else
            call crtbp_da_propagate(nominal_state, cov, MU, t_end, &
                da_order, n_samples, rel_tol, abs_tol, dt_min, dt_max, max_steps_int, &
                final_samples, final_mean, final_cov, propagated_ref, .false.)
        end if

        if (allocated(final_samples)) deallocate(final_samples)

        call dim_mean_and_cov(final_mean, final_cov, LU, TU, mean_dim, cov_dim)

        do i = 1, 3
            sigma_pos(i) = sqrt(max(cov_dim(i,i), 0.0_DP))
            sigma_vel(i) = sqrt(max(cov_dim(i+3,i+3), 0.0_DP))
        end do
        sigma_pos_rms = sqrt(max(cov_dim(1,1) + cov_dim(2,2) + cov_dim(3,3), 0.0_DP))
        sigma_vel_rms = sqrt(max(cov_dim(4,4) + cov_dim(5,5) + cov_dim(6,6), 0.0_DP))

        T_dim = T_nd * TU / DAY_SEC

        write(out_unit, '(A,",",A,",",*(ES22.14,:,","))') &
            trim(center), trim(family), T_dim, &
            mean_dim(1:3), mean_dim(4:6), &
            sigma_pos(1:3), sigma_pos_rms, &
            sigma_vel(1:3), sigma_vel_rms

        write(*, '(A,I0,A,A,",",A)') '  Row ', n_rows, ': ', trim(center), trim(family)
    end do

    close(in_unit)
    close(out_unit)

    print *, '========================================'
    print '(A,I0,A)', 'Processed ', n_rows, ' rows'
    print '(A,A)',  'Results saved to: ', trim(output_csv)

contains

    subroutine make_output_filename(path, meth, out_path)
        character(len=*), intent(in)  :: path, meth
        character(len=*), intent(out) :: out_path
        character(len=512) :: basename
        integer :: dot_pos, sep_pos

        sep_pos = scan(trim(path), '/\', back=.true.)
        if (sep_pos > 0) then
            basename = path(sep_pos+1:)
        else
            basename = path
        end if

        dot_pos = scan(trim(basename), '.', back=.true.)
        if (dot_pos > 0) basename = basename(1:dot_pos-1)

        write(out_path, '(A,A,A,A,A)') 'output/', trim(basename), '_batch_uprop_', trim(meth), '.csv'
    end subroutine make_output_filename

    subroutine parse_csv_line(line, center, family, state, T_nd)
        character(len=*), intent(in)  :: line
        character(len=*), intent(out) :: center, family
        real(DP), intent(out) :: state(6), T_nd
        character(len=32) :: orient
        real(DP) :: C_val, nu_val, k_val

        read(line, *) center, family, orient, &
            state(1), state(2), state(3), state(4), state(5), state(6), &
            T_nd, C_val, nu_val, k_val
    end subroutine parse_csv_line

    subroutine nondim_cov(cov, LU, TU)
        real(DP), intent(inout) :: cov(6,6)
        real(DP), intent(in)    :: LU, TU
        real(DP) :: pp_factor, pv_factor, vv_factor
        integer :: i, j

        pp_factor = 1.0_DP / (LU * LU)
        pv_factor = TU / (LU * LU)
        vv_factor = (TU * TU) / (LU * LU)

        do i = 1, 3
            do j = 1, 3
                cov(i,j) = cov(i,j) * pp_factor
            end do
        end do
        do i = 1, 3
            do j = 4, 6
                cov(i,j) = cov(i,j) * pv_factor
                cov(j,i) = cov(j,i) * pv_factor
            end do
        end do
        do i = 4, 6
            do j = 4, 6
                cov(i,j) = cov(i,j) * vv_factor
            end do
        end do
    end subroutine nondim_cov

    subroutine dim_mean_and_cov(mean_nd, cov_nd, LU, TU, mean_dim, cov_dim)
        real(DP), intent(in)  :: mean_nd(6), cov_nd(6,6)
        real(DP), intent(in)  :: LU, TU
        real(DP), intent(out) :: mean_dim(6), cov_dim(6,6)
        real(DP) :: pp_factor, pv_factor, vv_factor
        integer :: i, j

        pp_factor = LU * LU
        pv_factor = LU * LU / TU
        vv_factor = LU * LU / (TU * TU)

        mean_dim(1:3) = mean_nd(1:3) * LU
        mean_dim(4:6) = mean_nd(4:6) * LU / TU

        do i = 1, 3
            do j = 1, 3
                cov_dim(i,j) = cov_nd(i,j) * pp_factor
            end do
        end do
        do i = 1, 3
            do j = 4, 6
                cov_dim(i,j) = cov_nd(i,j) * pv_factor
                cov_dim(j,i) = cov_nd(j,i) * pv_factor
            end do
        end do
        do i = 4, 6
            do j = 4, 6
                cov_dim(i,j) = cov_nd(i,j) * vv_factor
            end do
        end do
    end subroutine dim_mean_and_cov

end program run_batch_uprop
