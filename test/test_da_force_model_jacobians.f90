program test_da_force_model_jacobians
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config
    use pod_spice, only: str2et
    use pod_force_model_module, only: compute_gravity_network, compute_atmospheric_drag, &
                                      compute_solar_radiation_pressure, compute_post_newtonian, &
                                      init_real_gravity => init_gravity_network
    use pod_da_force_model_module, only: da_compute_gravity_network, da_compute_atmospheric_drag, &
                                         da_compute_solar_radiation_pressure, da_compute_post_newtonian, &
                                         ForceModelTempPool, &
                                         init_da_gravity => init_gravity_network, &
                                         cleanup_da_gravity => cleanup_gravity_network
    use pod_dace_classes
    implicit none

    character(len=*), parameter :: CONFIG_FILE = 'config/dummy_test_config.txt'
    character(len=*), parameter :: TEST_EPOCH = '2024-03-09T12:00:00'
    real(DP), parameter :: H_POS = 1.0e-1_DP
    real(DP), parameter :: H_DRAG_POS = 1.0e-3_DP
    real(DP), parameter :: H_VEL = 1.0e-6_DP

    real(DP) :: et
    integer :: n_fail

    call pod_engine_init(CONFIG_FILE)
    call str2et(TEST_EPOCH, et)
    call dace_initialize(2, 6)

    n_fail = 0
    call test_gravity_position_jacobian(et, n_fail)
    call test_drag_jacobians(n_fail)
    call test_srp_position_jacobian(et, n_fail)
    call test_relativity_jacobians(et, n_fail)

    call cleanup_da_gravity()

    if (n_fail /= 0) then
        write(*,*) 'test_da_force_model_jacobians failed: ', n_fail
        stop 1
    end if
    write(*,*) 'test_da_force_model_jacobians passed'

contains

    subroutine test_gravity_position_jacobian(epoch, n_fail)
        real(DP), intent(in) :: epoch
        integer, intent(inout) :: n_fail
        real(DP) :: pos(3), vel(3), jac_da(3,3), jac_fd(3,3)
        type(AlgebraicVector) :: pos_da, vel_da, acc_da
        type(ForceModelTempPool) :: pool

        call configure_gravity([.false., .false., .true., .false., .false., .false., &
                                .false., .false., .false., .false., .false.])
        pos = [100000.0_DP, 50000.0_DP, 20000.0_DP]
        vel = [1.5_DP, 2.5_DP, 0.5_DP]
        call init_da_state(pos, vel, pos_da, vel_da, acc_da, pool)
        call da_compute_gravity_network(pos_da, epoch, acc_da, pool)
        call extract_jacobian(acc_da, jac_da, 1)
        call finite_difference_position_gravity(pos, epoch, H_POS, jac_fd)
        call assert_jacobian_close('gravity position', jac_da, jac_fd, 2.0e-7_DP, 1.0e-20_DP, n_fail)
        call destroy_da_state(pos_da, vel_da, acc_da, pool)
    end subroutine test_gravity_position_jacobian

    subroutine test_drag_jacobians(n_fail)
        integer, intent(inout) :: n_fail
        real(DP) :: pos(3), vel(3), jac_da_pos(3,3), jac_da_vel(3,3)
        real(DP) :: jac_fd_pos(3,3), jac_fd_vel(3,3)
        type(AlgebraicVector) :: pos_da, vel_da, acc_da
        type(ForceModelTempPool) :: pool

        pos = [6428.137_DP, 20.0_DP, -10.0_DP]
        vel = [0.1_DP, 7.5_DP, -0.2_DP]
        call init_da_state(pos, vel, pos_da, vel_da, acc_da, pool)
        call da_compute_atmospheric_drag(pos_da, vel_da, acc_da, pool)
        call extract_jacobian(acc_da, jac_da_pos, 1)
        call extract_jacobian(acc_da, jac_da_vel, 4)
        call finite_difference_position_drag(pos, vel, H_DRAG_POS, jac_fd_pos)
        call finite_difference_velocity_drag(pos, vel, H_VEL, jac_fd_vel)
        call assert_jacobian_close('drag position', jac_da_pos, jac_fd_pos, 2.0e-6_DP, 1.0e-20_DP, n_fail)
        call assert_jacobian_close('drag velocity', jac_da_vel, jac_fd_vel, 2.0e-6_DP, 1.0e-20_DP, n_fail)
        call destroy_da_state(pos_da, vel_da, acc_da, pool)
    end subroutine test_drag_jacobians

    subroutine test_srp_position_jacobian(epoch, n_fail)
        real(DP), intent(in) :: epoch
        integer, intent(inout) :: n_fail
        real(DP) :: pos(3), vel(3), jac_da(3,3), jac_fd(3,3)
        type(AlgebraicVector) :: pos_da, vel_da, acc_da
        type(ForceModelTempPool) :: pool

        pos = [100000.0_DP, 50000.0_DP, 20000.0_DP]
        vel = [1.5_DP, 2.5_DP, 0.5_DP]
        call init_da_state(pos, vel, pos_da, vel_da, acc_da, pool)
        call da_compute_solar_radiation_pressure(pos_da, epoch, acc_da, pool)
        call extract_jacobian(acc_da, jac_da, 1)
        call finite_difference_position_srp(pos, epoch, H_POS, jac_fd)
        call assert_jacobian_close('SRP position', jac_da, jac_fd, 2.0e-6_DP, 1.0e-28_DP, n_fail)
        call destroy_da_state(pos_da, vel_da, acc_da, pool)
    end subroutine test_srp_position_jacobian

    subroutine test_relativity_jacobians(epoch, n_fail)
        real(DP), intent(in) :: epoch
        integer, intent(inout) :: n_fail
        logical :: enabled_bodies(11)
        real(DP) :: pos(3), vel(3), jac_da_pos(3,3), jac_da_vel(3,3)
        real(DP) :: jac_fd_pos(3,3), jac_fd_vel(3,3)
        type(AlgebraicVector) :: pos_da, vel_da, acc_da
        type(ForceModelTempPool) :: pool

        enabled_bodies = .false.
        enabled_bodies(3) = .true.
        enabled_bodies(10) = .true.
        enabled_bodies(11) = .true.
        call configure_gravity(enabled_bodies)
        pos = [100000.0_DP, 50000.0_DP, 20000.0_DP]
        vel = [1.5_DP, 2.5_DP, 0.5_DP]
        call init_da_state(pos, vel, pos_da, vel_da, acc_da, pool)
        call da_compute_post_newtonian(pos_da, vel_da, epoch, acc_da, pool)
        call extract_jacobian(acc_da, jac_da_pos, 1)
        call extract_jacobian(acc_da, jac_da_vel, 4)
        call finite_difference_position_relativity(pos, vel, epoch, H_POS, jac_fd_pos)
        call finite_difference_velocity_relativity(pos, vel, epoch, H_VEL, jac_fd_vel)
        call assert_jacobian_close('relativity position', jac_da_pos, jac_fd_pos, 2.0e-5_DP, 1.0e-24_DP, n_fail)
        call assert_jacobian_close('relativity velocity', jac_da_vel, jac_fd_vel, 2.0e-5_DP, 1.0e-24_DP, n_fail)
        call destroy_da_state(pos_da, vel_da, acc_da, pool)
    end subroutine test_relativity_jacobians

    subroutine configure_gravity(enabled_bodies)
        logical, intent(in) :: enabled_bodies(11)
        call cleanup_da_gravity()
        config%use_planet = enabled_bodies
        config%use_earth_nspheric = .false.
        config%use_moon_nspheric = .false.
        config%use_third_body = any(enabled_bodies([1,2,4,5,6,7,8,9,10,11]))
        call init_real_gravity()
        call init_da_gravity()
    end subroutine configure_gravity

    subroutine init_da_state(pos, vel, pos_da, vel_da, acc_da, pool)
        real(DP), intent(in) :: pos(3), vel(3)
        type(AlgebraicVector), intent(inout) :: pos_da, vel_da, acc_da
        type(ForceModelTempPool), intent(inout) :: pool
        integer :: i

        call pos_da%init(3)
        call vel_da%init(3)
        call acc_da%init(3)
        call pool%init(3)
        do i = 1, 3
            pos_da%elements(i) = pos(i) + da_var(i)
            vel_da%elements(i) = vel(i) + da_var(i+3)
        end do
    end subroutine init_da_state

    subroutine destroy_da_state(pos_da, vel_da, acc_da, pool)
        type(AlgebraicVector), intent(inout) :: pos_da, vel_da, acc_da
        type(ForceModelTempPool), intent(inout) :: pool
        call pool%destroy()
        call pos_da%destroy()
        call vel_da%destroy()
        call acc_da%destroy()
    end subroutine destroy_da_state

    subroutine extract_jacobian(acc_da, jacobian, first_var)
        type(AlgebraicVector), intent(in) :: acc_da
        real(DP), intent(out) :: jacobian(3,3)
        integer, intent(in) :: first_var
        integer :: i, j
        do i = 1, 3
            do j = 1, 3
                jacobian(i,j) = acc_da%elements(i)%get_deriv_value(first_var + j - 1)
            end do
        end do
    end subroutine extract_jacobian

    subroutine assert_jacobian_close(label, actual, expected, rtol, atol, n_fail)
        character(len=*), intent(in) :: label
        real(DP), intent(in) :: actual(3,3), expected(3,3), rtol, atol
        integer, intent(inout) :: n_fail
        real(DP) :: max_error, scale, limit

        max_error = maxval(abs(actual - expected))
        scale = maxval(abs(expected))
        limit = atol + rtol * scale
        if (max_error > limit) then
            write(*,'(A,A)') 'FAIL: ', trim(label)
            write(*,'(A,ES14.6,A,ES14.6,A,ES14.6)') &
                '  max_error=', max_error, ' scale=', scale, ' limit=', limit
            n_fail = n_fail + 1
        else
            write(*,'(A,A,A,ES14.6)') 'PASS: ', trim(label), ' relative_error=', max_error / max(scale, atol)
        end if
    end subroutine assert_jacobian_close

    subroutine finite_difference_position_gravity(pos, epoch, h, jacobian)
        real(DP), intent(in) :: pos(3), epoch, h
        real(DP), intent(out) :: jacobian(3,3)
        real(DP) :: plus(3), minus(3), acc_plus(3), acc_minus(3)
        integer :: j
        do j = 1, 3
            plus = pos; minus = pos
            plus(j) = plus(j) + h; minus(j) = minus(j) - h
            call compute_gravity_network(plus, epoch, acc_plus)
            call compute_gravity_network(minus, epoch, acc_minus)
            jacobian(:,j) = (acc_plus - acc_minus) / (2.0_DP * h)
        end do
    end subroutine finite_difference_position_gravity

    subroutine finite_difference_position_drag(pos, vel, h, jacobian)
        real(DP), intent(in) :: pos(3), vel(3), h
        real(DP), intent(out) :: jacobian(3,3)
        real(DP) :: plus(3), minus(3), acc_plus(3), acc_minus(3)
        integer :: j
        do j = 1, 3
            plus = pos; minus = pos
            plus(j) = plus(j) + h; minus(j) = minus(j) - h
            call compute_atmospheric_drag(plus, vel, acc_plus)
            call compute_atmospheric_drag(minus, vel, acc_minus)
            jacobian(:,j) = (acc_plus - acc_minus) / (2.0_DP * h)
        end do
    end subroutine finite_difference_position_drag

    subroutine finite_difference_velocity_drag(pos, vel, h, jacobian)
        real(DP), intent(in) :: pos(3), vel(3), h
        real(DP), intent(out) :: jacobian(3,3)
        real(DP) :: plus(3), minus(3), acc_plus(3), acc_minus(3)
        integer :: j
        do j = 1, 3
            plus = vel; minus = vel
            plus(j) = plus(j) + h; minus(j) = minus(j) - h
            call compute_atmospheric_drag(pos, plus, acc_plus)
            call compute_atmospheric_drag(pos, minus, acc_minus)
            jacobian(:,j) = (acc_plus - acc_minus) / (2.0_DP * h)
        end do
    end subroutine finite_difference_velocity_drag

    subroutine finite_difference_position_srp(pos, epoch, h, jacobian)
        real(DP), intent(in) :: pos(3), epoch, h
        real(DP), intent(out) :: jacobian(3,3)
        real(DP) :: plus(3), minus(3), acc_plus(3), acc_minus(3)
        integer :: j
        do j = 1, 3
            plus = pos; minus = pos
            plus(j) = plus(j) + h; minus(j) = minus(j) - h
            call compute_solar_radiation_pressure(plus, epoch, acc_plus)
            call compute_solar_radiation_pressure(minus, epoch, acc_minus)
            jacobian(:,j) = (acc_plus - acc_minus) / (2.0_DP * h)
        end do
    end subroutine finite_difference_position_srp

    subroutine finite_difference_position_relativity(pos, vel, epoch, h, jacobian)
        real(DP), intent(in) :: pos(3), vel(3), epoch, h
        real(DP), intent(out) :: jacobian(3,3)
        real(DP) :: plus(3), minus(3), acc_plus(3), acc_minus(3)
        integer :: j
        do j = 1, 3
            plus = pos; minus = pos
            plus(j) = plus(j) + h; minus(j) = minus(j) - h
            call compute_post_newtonian(plus, vel, epoch, acc_plus)
            call compute_post_newtonian(minus, vel, epoch, acc_minus)
            jacobian(:,j) = (acc_plus - acc_minus) / (2.0_DP * h)
        end do
    end subroutine finite_difference_position_relativity

    subroutine finite_difference_velocity_relativity(pos, vel, epoch, h, jacobian)
        real(DP), intent(in) :: pos(3), vel(3), epoch, h
        real(DP), intent(out) :: jacobian(3,3)
        real(DP) :: plus(3), minus(3), acc_plus(3), acc_minus(3)
        integer :: j
        do j = 1, 3
            plus = vel; minus = vel
            plus(j) = plus(j) + h; minus(j) = minus(j) - h
            call compute_post_newtonian(pos, plus, epoch, acc_plus)
            call compute_post_newtonian(pos, minus, epoch, acc_minus)
            jacobian(:,j) = (acc_plus - acc_minus) / (2.0_DP * h)
        end do
    end subroutine finite_difference_velocity_relativity

end program test_da_force_model_jacobians
