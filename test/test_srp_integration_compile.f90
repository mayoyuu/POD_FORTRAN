program test_srp_integration_compile
    use pod_force_model_module, only: compute_solar_radiation_pressure
    use pod_da_force_model_module, only: build_srp_da_parameter_map
    use pod_uq_da_module, only: uq_da_propagator
    implicit none
    type(uq_da_propagator) :: propagator

    if (propagator%da_order <= 0) stop 1
    write(*,*) 'SRP real/DA/UQ integration compile test passed.'
end program test_srp_integration_compile
