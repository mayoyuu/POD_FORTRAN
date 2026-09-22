!> Fixed 7D study grid, independent SRP truth probes, and history diagnostics.
module pod_srp_ads_history_module
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use pod_global, only: DP
    use pod_config, only: config
    use pod_uq_hfem_ads_module, only: hfem_ads_options_type
    use pod_force_model_module, only: set_propagation_epoch, &
        set_srp_scale_error, clear_srp_scale_error
    use pod_integrator_module, only: adaptive_step_integrate, METHOD_RKF78
    use pod_uncertainty_diagnostics_module, only: compute_covariance_axes
    implicit none
    private
    integer, parameter, public :: SHAPE_COUNT = 512
    integer, parameter, public :: VALIDATION_COUNT = 207
    public :: build_study_points, assess_accuracy, advance_real_probes
    public :: compute_shape_metrics, validate_json_opm_fields

contains

    subroutine build_study_points(points, validation_ids)
        real(DP), allocatable, intent(out) :: points(:,:)
        integer, allocatable, intent(out) :: validation_ids(:)
        integer, parameter :: primes(7) = [2,3,5,7,11,13,17]
        integer :: j, k

        allocate(points(7,SHAPE_COUNT), validation_ids(VALIDATION_COUNT))
        points = 0.0_DP
        do j = 0, 127
            do k = 1, 7
                points(k,j+1) = merge(1.0_DP,-1.0_DP,btest(j,k-1))
            end do
        end do
        do k = 1, 7
            points(k,128+2*k-1) = -1.0_DP
            points(k,128+2*k) = 1.0_DP
        end do
        ! Point 143 is the center. Interior points start at 144.
        do j = 1, SHAPE_COUNT-143
            do k = 1, 7
                points(k,143+j) = 2.0_DP*halton(j,primes(k))-1.0_DP
            end do
        end do
        validation_ids = [(j,j=1,VALIDATION_COUNT)]
    end subroutine build_study_points

    pure real(DP) function halton(index_value, base) result(value)
        integer, intent(in) :: index_value, base
        integer :: n
        real(DP) :: denominator

        n = index_value
        value = 0.0_DP
        denominator = real(base,DP)
        do while (n > 0)
            value = value + real(mod(n,base),DP)/denominator
            n = n/base
            denominator = denominator*real(base,DP)
        end do
    end function halton

    subroutine assess_accuracy(ads, truth, pos_tolerance, vel_tolerance, &
            max_position, max_velocity, worst_id, failed)
        real(DP), intent(in) :: ads(:,:), truth(:,:)
        real(DP), intent(in) :: pos_tolerance, vel_tolerance
        real(DP), intent(out) :: max_position, max_velocity
        integer, intent(out) :: worst_id
        logical, intent(out) :: failed
        real(DP) :: position_error, velocity_error, ratio, worst_ratio
        integer :: j

        max_position = 0.0_DP
        max_velocity = 0.0_DP
        worst_ratio = 0.0_DP
        worst_id = 0
        failed = .false.
        do j = 1, size(ads,2)
            position_error = sqrt(sum((ads(1:3,j)-truth(1:3,j))**2))
            velocity_error = sqrt(sum((ads(4:6,j)-truth(4:6,j))**2))
            max_position = max(max_position,position_error)
            max_velocity = max(max_velocity,velocity_error)
            ratio = max(position_error/pos_tolerance,velocity_error/vel_tolerance)
            if (ratio > worst_ratio) then
                worst_ratio = ratio
                worst_id = j
            end if
        end do
        failed = max_position > pos_tolerance .or. max_velocity > vel_tolerance
    end subroutine assess_accuracy

    subroutine advance_real_probes(states, eta, epoch0, t_start, t_end, options, &
            status, message)
        real(DP), intent(inout) :: states(:,:)
        real(DP), intent(in) :: eta(:), epoch0, t_start, t_end
        type(hfem_ads_options_type), intent(in) :: options
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        real(DP) :: nondim(6)
        real(DP), allocatable :: times(:), integrated(:,:)
        integer :: j, n_steps

        status = 0
        message = ''
        if (size(states,1) /= 6 .or. size(states,2) /= size(eta)) then
            status = -1
            message = 'real SRP probe dimensions do not match'
            return
        end if
        if (t_end < t_start) then
            status = -2
            message = 'real SRP probe time goes backward'
            return
        end if
        call set_propagation_epoch(epoch0)
        do j = 1, size(eta)
            call set_srp_scale_error(eta(j))
            nondim(1:3) = states(1:3,j)/config%LU
            nondim(4:6) = states(4:6,j)/config%VU
            call adaptive_step_integrate(nondim,t_start/config%TU, &
                t_end/config%TU,METHOD_RKF78,times,integrated,n_steps, &
                options%max_steps,options%rel_tol,options%abs_tol, &
                options%dt_min,options%dt_max)
            if (n_steps < 1 .or. .not. allocated(integrated)) then
                call clear_srp_scale_error()
                status = -3
                message = 'real SRP probe integration produced no state'
                return
            end if
            if (.not. all(ieee_is_finite(integrated(n_steps,:)))) then
                call clear_srp_scale_error()
                status = -3
                message = 'real SRP probe integration produced a nonfinite state'
                return
            end if
            states(1:3,j) = integrated(n_steps,1:3)*config%LU
            states(4:6,j) = integrated(n_steps,4:6)*config%VU
        end do
        call clear_srp_scale_error()
    end subroutine advance_real_probes

    subroutine compute_shape_metrics(states, nominal, lower, upper, &
            max_position, max_velocity, position_axes, velocity_axes, status)
        real(DP), intent(in) :: states(:,:), nominal(6)
        real(DP), intent(out) :: lower(6), upper(6), max_position, max_velocity
        real(DP), intent(out) :: position_axes(3), velocity_axes(3)
        integer, intent(out) :: status
        real(DP) :: mean(6), cov3(3,3), eig(3), vec(3,3)
        real(DP) :: deviation(6)
        integer :: j, k, l, n, axis_status

        status = 0
        lower = 0.0_DP
        upper = 0.0_DP
        position_axes = 0.0_DP
        velocity_axes = 0.0_DP
        max_position = 0.0_DP
        max_velocity = 0.0_DP
        if (size(states,1) /= 6 .or. size(states,2) < 2) then
            status = -1
            return
        end if
        if (.not. all(ieee_is_finite(states))) then
            status = -2
            return
        end if
        n = size(states,2)
        do k = 1, 6
            lower(k) = minval(states(k,:)-nominal(k))
            upper(k) = maxval(states(k,:)-nominal(k))
        end do
        mean = sum(states,dim=2)/real(n,DP)
        do j = 1, n
            deviation = states(:,j)-nominal
            max_position = max(max_position,sqrt(sum(deviation(1:3)**2)))
            max_velocity = max(max_velocity,sqrt(sum(deviation(4:6)**2)))
        end do
        do l = 0, 1
            cov3 = 0.0_DP
            do j = 1, n
                do k = 1, 3
                    cov3(:,k) = cov3(:,k) + &
                        (states(1+3*l:3+3*l,j)-mean(1+3*l:3+3*l))* &
                        (states(k+3*l,j)-mean(k+3*l))
                end do
            end do
            cov3 = cov3/real(n-1,DP)
            if (l == 0) then
                call compute_covariance_axes(cov3,eig,vec,position_axes,axis_status)
            else
                call compute_covariance_axes(cov3,eig,vec,velocity_axes,axis_status)
            end if
            if (axis_status /= 0) then
                status = -3
                return
            end if
        end do
    end subroutine compute_shape_metrics

    subroutine validate_json_opm_fields(filename, status, message)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        character(len=8192) :: line
        character(len=5), parameter :: keys(6) = [character(len=5) :: &
            'X','Y','Z','X_DOT','Y_DOT','Z_DOT']
        logical :: state_seen(6), covariance_seen(6,6), epoch_seen
        integer :: unit, io, depth, i, j, ch

        status = 0
        message = ''
        state_seen = .false.
        covariance_seen = .false.
        epoch_seen = .false.
        open(newunit=unit,file=filename,status='old',action='read',iostat=io)
        if (io /= 0) then
            status = -1
            message = 'cannot open OPM file'
            return
        end if
        depth = 0
        do
            read(unit,'(A)',iostat=io) line
            if (io < 0) exit
            if (io /= 0) then
                status = -2
                message = 'cannot read OPM file'
                exit
            end if
            do ch = 1, len_trim(line)
                if (line(ch:ch) == '{') depth = depth+1
                if (line(ch:ch) == '}') depth = depth-1
            end do
            if (depth /= 1) cycle
            if (index(line,'"EPOCH"') > 0) epoch_seen = .true.
            if (index(line,'"REF_FRAME"') > 0) then
                if (index(line,'GCRF') == 0 .and. index(line,'J2000') == 0) then
                    status = -3
                    message = 'OPM REF_FRAME must be GCRF or J2000'
                    exit
                end if
            end if
            do i = 1, 6
                if (index(line,'"'//trim(keys(i))//'"') > 0) state_seen(i) = .true.
                do j = 1, i
                    if (index(line,'"C'//trim(keys(i))//'_'//trim(keys(j))//'"') > 0) &
                        covariance_seen(i,j) = .true.
                end do
            end do
        end do
        close(unit)
        if (status /= 0) return
        if (.not. epoch_seen .or. .not. all(state_seen)) then
            status = -4
            message = 'OPM lacks top-level epoch or Cartesian state'
            return
        end if
        do i = 1, 6
            do j = 1, i
                if (.not. covariance_seen(i,j)) then
                    status = -5
                    message = 'OPM lacks a covariance element'
                    return
                end if
            end do
        end do
    end subroutine validate_json_opm_fields

end module pod_srp_ads_history_module
