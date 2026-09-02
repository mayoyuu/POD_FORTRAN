program test_srp_panel_real
    use pod_global, only: DP
    use pod_spacecraft_geometry, only: srp_optical_properties_type, compute_panel_force_real
    implicit none

    real(DP), parameter :: pressure = 4.0e-6_DP, area = 5.0_DP, tol = 1.0e-18_DP
    real(DP) :: e_s(3), normal(3), force(3)
    type(srp_optical_properties_type) :: optical

    e_s = [1.0_DP, 0.0_DP, 0.0_DP]
    normal = e_s

    optical%absorptivity = 1.0_DP
    optical%specular_reflectivity = 0.0_DP
    optical%diffuse_reflectivity = 0.0_DP
    call compute_panel_force_real(e_s, normal, area, pressure, optical, force)
    call assert_close(force, [-pressure*area, 0.0_DP, 0.0_DP], 'normal perfect absorption')

    optical%absorptivity = 0.0_DP
    optical%specular_reflectivity = 1.0_DP
    call compute_panel_force_real(e_s, normal, area, pressure, optical, force)
    call assert_close(force, [-2.0_DP*pressure*area, 0.0_DP, 0.0_DP], 'normal perfect specular')

    normal = [0.0_DP, 1.0_DP, 0.0_DP]
    call compute_panel_force_real(e_s, normal, area, pressure, optical, force)
    call assert_close(force, [0.0_DP, 0.0_DP, 0.0_DP], 'grazing incidence')

    normal = [-1.0_DP, 0.0_DP, 0.0_DP]
    call compute_panel_force_real(e_s, normal, area, pressure, optical, force)
    call assert_close(force, [0.0_DP, 0.0_DP, 0.0_DP], 'single-sided back face')

    write(*,*) 'Real SRP panel tests passed.'

contains
    subroutine assert_close(actual, expected, label)
        real(DP), intent(in) :: actual(3), expected(3)
        character(len=*), intent(in) :: label
        if (maxval(abs(actual - expected)) > tol) then
            write(*,*) 'FAILED: ', trim(label), actual, expected
            stop 1
        end if
    end subroutine assert_close
end program test_srp_panel_real
