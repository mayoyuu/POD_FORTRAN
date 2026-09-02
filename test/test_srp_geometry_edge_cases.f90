program test_srp_geometry_edge_cases
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use pod_global, only: DP
    use pod_spacecraft_geometry, only: compute_pointing_attitude_real, compute_array_normal_real, &
                                       spacecraft_geometry_type
    implicit none

    real(DP) :: position(3), velocity(3), sun_position(3), earth_position(3), moon_position(3)
    real(DP) :: c_i_b(3,3), normal(3), acceleration(3), invalid_c(3,3)
    type(spacecraft_geometry_type) :: geometry
    integer :: status
    character(len=256) :: message

    position = 0.0_DP
    velocity = [0.0_DP,1.0_DP,0.0_DP]
    sun_position = [1.0_DP,0.0_DP,0.0_DP]
    earth_position = [-1.0_DP,0.0_DP,0.0_DP]
    moon_position = [0.0_DP,1.0_DP,0.0_DP]
    call compute_pointing_attitude_real(position, velocity, sun_position, earth_position, moon_position, &
                                        'sun', 'sun', [1.0_DP,0.0_DP,0.0_DP], &
                                        [0.0_DP,1.0_DP,0.0_DP], c_i_b, status, message)
    call assert_true(status > 0, 'parallel roll reference uses fallback')
    call assert_true(all(ieee_is_finite(c_i_b)), 'fallback attitude is finite')

    call compute_pointing_attitude_real(position, velocity, sun_position, earth_position, moon_position, &
                                        'sun', 'orbit_normal', [1.0_DP,0.0_DP,0.0_DP], &
                                        [0.0_DP,0.0_DP,1.0_DP], c_i_b, status, message, &
                                        sun_velocity=[0.0_DP,2.0_DP,0.0_DP])
    call assert_true(status == 0, 'relative target velocity attitude status')
    call assert_true(dot_product(matmul(c_i_b,[0.0_DP,0.0_DP,1.0_DP]), &
                                 [0.0_DP,0.0_DP,1.0_DP]) > 1.0_DP-1.0e-13_DP, &
                     'orbit normal uses spacecraft minus target velocity')

    call compute_array_normal_real(c_i_b, matmul(c_i_b,[0.0_DP,1.0_DP,0.0_DP]), 'single_axis', &
                                   [0.0_DP,1.0_DP,0.0_DP], [1.0_DP,0.0_DP,0.0_DP], 0.0_DP, &
                                   normal, status, message)
    call assert_true(status > 0, 'Sun parallel to hinge uses reference normal')
    call assert_true(all(ieee_is_finite(normal)), 'fallback array normal is finite')

    call configure_geometry(geometry)
    invalid_c = 0.0_DP
    invalid_c(1,1) = 2.0_DP
    invalid_c(2,2) = 1.0_DP
    invalid_c(3,3) = 1.0_DP
    call geometry%compute_srp_with_attitude_real(position, sun_position, invalid_c, acceleration, &
                                                  status, message)
    call assert_true(status < 0, 'invalid external attitude rejected')

    geometry%box_optical%absorptivity = 0.8_DP
    call geometry%validate(status, message)
    call assert_true(status < 0, 'local geometry rejects invalid optical sum')
    call configure_geometry(geometry)
    geometry%pressure_1au_n_m2 = 0.0_DP
    call geometry%validate(status, message)
    call assert_true(status < 0, 'local geometry rejects zero pressure')
    write(*,*) 'SRP geometry edge-case tests passed.'

contains
    subroutine configure_geometry(g)
        type(spacecraft_geometry_type), intent(out) :: g
        g%mass_kg = 1000.0_DP
        g%box_dimensions_m = [2.0_DP,3.0_DP,4.0_DP]
        g%box_optical%absorptivity = 1.0_DP
        g%box_optical%specular_reflectivity = 0.0_DP
        g%box_optical%diffuse_reflectivity = 0.0_DP
        g%array_total_area_m2 = 0.0_DP
        g%pressure_1au_n_m2 = 4.56e-6_DP
    end subroutine configure_geometry
    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not. condition) then
            write(*,*) 'FAILED: ', trim(label)
            stop 1
        end if
    end subroutine assert_true
end program test_srp_geometry_edge_cases
