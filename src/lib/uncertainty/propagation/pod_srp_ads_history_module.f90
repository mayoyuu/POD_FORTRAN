!> Dimension-aware study grids, independent SRP probes, and history diagnostics.
module pod_srp_ads_history_module
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use pod_global, only: DP
    use pod_config, only: config
    use pod_uq_hfem_ads_module, only: hfem_ads_options_type
    use pod_force_model_module, only: set_propagation_epoch, &
        set_srp_scale_error, clear_srp_scale_error, &
        set_srp_mode_coefficients_real,clear_srp_mode_coefficients_real
    use pod_srp_mode_module, only: interpolate_srp_modes
    use pod_integrator_module, only: adaptive_step_integrate, METHOD_RKF78
    use pod_uncertainty_diagnostics_module, only: compute_covariance_axes
    implicit none
    private
    integer, parameter, public :: SHAPE_COUNT = 512
    integer, parameter, public :: VALIDATION_COUNT = 207
    public :: build_study_points, assess_accuracy, advance_real_probes
    public :: compute_shape_metrics, validate_json_opm_fields
    public :: advance_real_mode_probes

contains

    subroutine build_study_points(points, validation_ids,n_variables)
        real(DP), allocatable, intent(out) :: points(:,:)
        integer, allocatable, intent(out) :: validation_ids(:)
        integer,intent(in),optional :: n_variables
        integer, parameter :: primes(8) = [2,3,5,7,11,13,17,19]
        integer :: j,k,n,ncorners,nbase,nvalidation

        n=7
        if(present(n_variables)) n=n_variables
        if(n<6 .or. n>8) error stop 'study grid supports 6, 7 or 8 variables'
        ncorners=2**n; nbase=ncorners+2*n+1
        nvalidation=VALIDATION_COUNT
        if(n==8) nvalidation=nbase+64
        allocate(points(n,SHAPE_COUNT), validation_ids(nvalidation))
        points = 0.0_DP
        do j = 0, ncorners-1
            do k = 1, n
                points(k,j+1) = merge(1.0_DP,-1.0_DP,btest(j,k-1))
            end do
        end do
        do k = 1, n
            points(k,ncorners+2*k-1) = -1.0_DP
            points(k,ncorners+2*k) = 1.0_DP
        end do
        ! Legacy 7D retains center=143 and every existing Halton point.
        do j = 1, SHAPE_COUNT-nbase
            do k = 1, n
                points(k,nbase+j) = 2.0_DP*halton(j,primes(k))-1.0_DP
            end do
        end do
        validation_ids = [(j,j=1,nvalidation)]
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
            ! A finite last checkpoint is not sufficient after an aborted RK step.
            if (.not. allocated(times)) then
                call clear_srp_scale_error()
                status = -3
                message = 'real SRP probe integration produced no end time'
                return
            end if
            if (n_steps > size(times)) then
                call clear_srp_scale_error()
                status = -3
                message = 'real SRP probe integration produced an invalid end time'
                return
            end if
            if (.not. ieee_is_finite(times(n_steps))) then
                call clear_srp_scale_error()
                status = -3
                message = 'real SRP probe integration produced a nonfinite end time'
                return
            end if
            if (abs(times(n_steps)-t_end/config%TU) > &
                1000*epsilon(1.0_DP)*max(1.0_DP,abs(t_end/config%TU))) then
                call clear_srp_scale_error()
                status = -3
                message = 'real SRP probe integration did not reach requested end time'
                return
            end if
            states(1:3,j) = integrated(n_steps,1:3)*config%LU
            states(4:6,j) = integrated(n_steps,4:6)*config%VU
        end do
        call clear_srp_scale_error()
    end subroutine advance_real_probes

    subroutine advance_real_mode_probes(states,s_values,u_values,epoch0,t_start,t_end,options,status,message)
        real(DP),intent(inout) :: states(:,:)
        real(DP),intent(in) :: s_values(:),u_values(:),epoch0,t_start,t_end
        type(hfem_ads_options_type),intent(in) :: options
        integer,intent(out) :: status
        character(len=*),intent(out) :: message
        real(DP) :: nondim(6),eE(3),eM(3)
        real(DP),allocatable :: times(:),integrated(:,:),working(:,:)
        integer :: j,n_steps
        call clear_srp_scale_error(); call clear_srp_mode_coefficients_real()
        status=-1; message='real mode probe dimensions do not match'
        if(size(states,1)/=6 .or. size(states,2)/=size(s_values) .or. size(s_values)/=size(u_values)) return
        status=-2; message='real mode probes require finite state, coefficients and forward time'
        if(.not.all(ieee_is_finite(states)) .or. .not.all(ieee_is_finite(s_values)) .or. &
           .not.all(ieee_is_finite(u_values)) .or. .not.all(ieee_is_finite([epoch0,t_start,t_end]))) return
        if(t_end<t_start .or. any(s_values<0) .or. any(s_values>1) .or. any(u_values<0) .or. any(u_values>1)) return
        status=-3; message='invalid Real mode integration options or force model'
        if(trim(config%srp_model)/='time_varying_modes') return
        if(options%max_steps<2 .or. .not.all(ieee_is_finite([options%rel_tol,options%abs_tol, &
            options%dt_min,options%dt_max])) .or. min(options%rel_tol,options%abs_tol,options%dt_min)<=0 .or. &
            options%dt_max<options%dt_min) return
        call interpolate_srp_modes(epoch0+t_start,eE,eM,status,message)
        if(status/=0) return
        call interpolate_srp_modes(epoch0+t_end,eE,eM,status,message)
        if(status/=0 .or. t_start==t_end) return
        working=states
        call set_propagation_epoch(epoch0)
        do j=1,size(s_values)
            call set_srp_mode_coefficients_real(s_values(j),u_values(j))
            nondim(1:3)=states(1:3,j)/config%LU; nondim(4:6)=states(4:6,j)/config%VU
            call adaptive_step_integrate(nondim,t_start/config%TU,t_end/config%TU,METHOD_RKF78, &
                times,integrated,n_steps,options%max_steps,options%rel_tol,options%abs_tol,options%dt_min,options%dt_max)
            status=-4; message='Real mode probe integration did not reach a finite end state'
            if(n_steps<1 .or. .not.allocated(times) .or. .not.allocated(integrated)) exit
            if(n_steps>size(times) .or. n_steps>size(integrated,1)) exit
            if(.not.all(ieee_is_finite(integrated(n_steps,:)))) exit
            if(abs(times(n_steps)-t_end/config%TU)>1000*epsilon(1.0_DP)*max(1.0_DP,abs(t_end/config%TU))) exit
            working(1:3,j)=integrated(n_steps,1:3)*config%LU
            working(4:6,j)=integrated(n_steps,4:6)*config%VU
            status=0; message=''
        end do
        call clear_srp_mode_coefficients_real(); call clear_srp_scale_error()
        if(status==0) states=working
    end subroutine

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
