program test_srp_param_da_propagation
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config
    use pod_spice, only: str2et
    use pod_uq_propagation, only: run_particle_propagation, METHOD_DA
    use pod_uq_state_module, only: uq_state_type
    use pod_dace_classes, only: dace_initialize
    implicit none

    character(len=*), parameter :: CONFIG_FILE = 'config/dummy_test_config.txt'
    character(len=*), parameter :: TEST_EPOCH = '2024-03-09T12:00:00'
    real(DP), parameter :: DT = 3600.0_DP
    integer, parameter :: DA_ORDER = 2
    integer, parameter :: N_PARTICLES = 3

    real(DP) :: epoch0
    real(DP) :: reference_orbit(6), reference_out(6)
    type(uq_state_type) :: state6, out6, state7, out7
    real(DP) :: eta_effect
    integer :: i

    call pod_engine_init(CONFIG_FILE)

    config%use_planet = .false.
    config%use_earth_nspheric = .false.
    config%use_moon_nspheric = .false.
    config%use_third_body = .false.
    config%use_relativity = .false.
    config%use_drag = .false.
    config%use_srp = .true.

    call str2et(TEST_EPOCH, epoch0)
    reference_orbit = [100000.0_DP, 50000.0_DP, 20000.0_DP, &
                       1.5_DP, 2.5_DP, 0.5_DP]

    call dace_initialize(DA_ORDER, 6)
    call state6%allocate_memory(6, N_PARTICLES)
    do i = 1, N_PARTICLES
        state6%samples(:, i) = reference_orbit
    end do

    call run_particle_propagation(state6, reference_orbit, epoch0, 0.0_DP, DT, &
                                  METHOD_DA, out6, da_order=DA_ORDER, &
                                  reference_orbit_out=reference_out)

    if (size(out6%samples, 1) /= 6) then
        write(*,*) 'Expected 6-D output for legacy 6-D propagation, got ', size(out6%samples, 1)
        stop 1
    end if

    call dace_initialize(DA_ORDER, 7)
    call state7%allocate_memory(7, N_PARTICLES)
    do i = 1, N_PARTICLES
        state7%samples(1:6, i) = reference_orbit
    end do
    state7%samples(7, 1) = -0.20_DP
    state7%samples(7, 2) =  0.00_DP
    state7%samples(7, 3) =  0.20_DP

    call run_particle_propagation(state7, reference_orbit, epoch0, 0.0_DP, DT, &
                                  METHOD_DA, out7, da_order=DA_ORDER, &
                                  reference_orbit_out=reference_out)

    if (size(out7%samples, 1) /= 7) then
        write(*,*) 'Expected 7-D output for SRP-parameter propagation, got ', size(out7%samples, 1)
        stop 1
    end if

    if (maxval(abs(out7%samples(7, :) - state7%samples(7, :))) > 1.0e-14_DP) then
        write(*,*) 'SRP parameter samples should be preserved as constant parameters.'
        stop 1
    end if

    eta_effect = maxval(abs(out7%samples(1:6, 3) - out7%samples(1:6, 1)))
    if (eta_effect <= 1.0e-12_DP) then
        write(*,*) 'Changing eta_srp did not affect the propagated orbit.'
        stop 1
    end if

    call state6%deallocate_memory()
    call out6%deallocate_memory()
    call state7%deallocate_memory()
    call out7%deallocate_memory()

    write(*,*) 'SRP DA parameter propagation test passed.'
end program test_srp_param_da_propagation
