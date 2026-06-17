program test_ads_domain_scale
    use pod_global, only: DP
    use pod_uq_crtbp_ads_module, only: crtbp_ads_propagate_deviates
    implicit none

    real(DP) :: nominal(6), deviates(6, 2), domain_scale(6)
    real(DP), allocatable :: final_samples(:,:)
    real(DP) :: final_mean(6), final_cov(6,6), err_toll(6)
    integer :: n_patches

    nominal = [1.1309265107780351_DP, 0.0_DP, 0.0_DP, &
               0.0_DP, -0.46540743845849059_DP, 0.0_DP]
    domain_scale = [2.0e-3_DP, 3.0e-3_DP, 4.0e-3_DP, &
                    5.0e-3_DP, 6.0e-3_DP, 7.0e-3_DP]

    deviates(:, 1) = 0.25_DP * domain_scale
    deviates(:, 2) = -0.50_DP * domain_scale
    err_toll = 1.0e30_DP

    call crtbp_ads_propagate_deviates(nominal, deviates, &
        0.012153614091892_DP, 0.0_DP, &
        2, 0, err_toll, &
        1.0e-14_DP, 1.0e-14_DP, 1.0e-10_DP, 0.01_DP, 10, &
        final_samples, final_mean, final_cov, n_patches, .false., &
        domain_scale=domain_scale)

    if (n_patches /= 1) then
        write(*,*) 'FAIL: expected one accepted ADS patch, got ', n_patches
        stop 1
    end if

    if (maxval(abs(final_samples(:, 1) - (nominal + deviates(:, 1)))) > 1.0e-12_DP) then
        write(*,*) 'FAIL: positive deviate was not evaluated in physical coordinates'
        write(*,*) 'got     ', final_samples(:, 1)
        write(*,*) 'expected', nominal + deviates(:, 1)
        stop 1
    end if

    if (maxval(abs(final_samples(:, 2) - (nominal + deviates(:, 2)))) > 1.0e-12_DP) then
        write(*,*) 'FAIL: negative deviate was not evaluated in physical coordinates'
        write(*,*) 'got     ', final_samples(:, 2)
        write(*,*) 'expected', nominal + deviates(:, 2)
        stop 1
    end if

    write(*,*) 'PASS: ADS domain_scale maps physical deviates through unit split coordinates'

    if (allocated(final_samples)) deallocate(final_samples)
end program test_ads_domain_scale
