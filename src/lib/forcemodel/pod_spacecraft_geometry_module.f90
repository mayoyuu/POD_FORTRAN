!> Simplified box-wing spacecraft geometry and solar-radiation-pressure core.
!!
!! Coordinate convention: C_I_B maps a body-frame vector into the inertial frame,
!! v_I = C_I_B v_B.  The primary body axis points at exactly one selected target
!! (Sun, Earth, or Moon); the secondary body axis only resolves roll.
module pod_spacecraft_geometry
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_config, only: config_params
    use pod_dace_classes
    implicit none
    private

    real(DP), parameter :: DEFAULT_GEOMETRY_TOLERANCE = 1.0e-12_DP
    real(DP), parameter :: AU_KM = 149597870.7_DP
    real(DP), parameter :: DEFAULT_SOLAR_PRESSURE_1AU = 1367.0_DP / 299792458.0_DP

    type, public :: srp_optical_properties_type
        real(DP) :: absorptivity = 1.0_DP
        real(DP) :: specular_reflectivity = 0.0_DP
        real(DP) :: diffuse_reflectivity = 0.0_DP
    end type srp_optical_properties_type

    type, public :: srp_panel_type
        real(DP) :: area_m2 = 0.0_DP
        real(DP) :: normal_body(3) = 0.0_DP
        type(srp_optical_properties_type) :: optical
    end type srp_panel_type

    type, public :: srp_da_parameter_map_type
        integer :: global_scale_index = 0
        integer :: attitude_bias_index(3) = 0
        integer :: array_angle_index = 0
        real(DP) :: global_scale_nominal = 0.0_DP
        real(DP) :: global_scale_span = 0.0_DP
        real(DP) :: attitude_bias_span_rad(3) = 0.0_DP
        real(DP) :: array_angle_span_rad = 0.0_DP
    end type srp_da_parameter_map_type

    type, public :: da_attitude_type
        type(AlgebraicVector) :: x_body_in_inertial
        type(AlgebraicVector) :: y_body_in_inertial
        type(AlgebraicVector) :: z_body_in_inertial
    contains
        procedure :: init => da_attitude_init
        procedure :: destroy => da_attitude_destroy
    end type da_attitude_type

    type, public :: spacecraft_geometry_type
        real(DP) :: mass_kg = -1.0_DP
        real(DP) :: box_dimensions_m(3) = -1.0_DP
        type(srp_optical_properties_type) :: box_optical
        real(DP) :: array_total_area_m2 = 0.0_DP
        character(len=MAX_STRING_LEN) :: array_tracking_mode = 'single_axis'
        real(DP) :: array_hinge_axis_body(3) = [0.0_DP, 1.0_DP, 0.0_DP]
        real(DP) :: array_reference_normal_body(3) = [1.0_DP, 0.0_DP, 0.0_DP]
        type(srp_optical_properties_type) :: array_front_optical
        type(srp_optical_properties_type) :: array_back_optical
        character(len=MAX_STRING_LEN) :: attitude_mode = 'sun'
        character(len=MAX_STRING_LEN) :: roll_reference = 'orbit_normal'
        real(DP) :: primary_axis_body(3) = [0.0_DP, 0.0_DP, 1.0_DP]
        real(DP) :: secondary_axis_body(3) = [0.0_DP, 1.0_DP, 0.0_DP]
        real(DP) :: pressure_1au_n_m2 = DEFAULT_SOLAR_PRESSURE_1AU
        real(DP) :: geometry_tolerance = DEFAULT_GEOMETRY_TOLERANCE
    contains
        procedure :: initialize_from_config => initialize_geometry_from_config
        procedure :: validate => validate_spacecraft_geometry
        procedure :: compute_srp_from_ephemerides_real
        procedure :: compute_srp_with_attitude_real
        procedure :: compute_srp_from_ephemerides_da
        procedure :: compute_srp_with_attitude_da
        procedure :: compute_srp_with_real_attitude_da
    end type spacecraft_geometry_type

    public :: compute_pointing_attitude_real
    public :: compute_array_normal_real
    public :: compute_panel_force_real
    public :: compute_pointing_attitude_da
    public :: compute_array_normal_da
    public :: compute_panel_force_da
    public :: normalize_vector_real
    public :: cross_product_real

contains

    subroutine da_attitude_init(self)
        class(da_attitude_type), intent(inout) :: self
        call self%x_body_in_inertial%init(3)
        call self%y_body_in_inertial%init(3)
        call self%z_body_in_inertial%init(3)
    end subroutine da_attitude_init

    subroutine da_attitude_destroy(self)
        class(da_attitude_type), intent(inout) :: self
        call self%x_body_in_inertial%destroy()
        call self%y_body_in_inertial%destroy()
        call self%z_body_in_inertial%destroy()
    end subroutine da_attitude_destroy

    subroutine compute_pointing_attitude_da(position, velocity, sun_position, earth_position, &
                                             moon_position, attitude_mode, roll_reference, &
                                             primary_axis_body, secondary_axis_body, attitude, &
                                             status, message, tolerance, sun_velocity, earth_velocity, &
                                             moon_velocity)
        type(AlgebraicVector), intent(in) :: position, velocity
        real(DP), intent(in) :: sun_position(3), earth_position(3), moon_position(3)
        character(len=*), intent(in) :: attitude_mode, roll_reference
        real(DP), intent(in) :: primary_axis_body(3), secondary_axis_body(3)
        type(da_attitude_type), intent(inout) :: attitude
        integer, intent(out) :: status
        character(len=*), intent(out), optional :: message
        real(DP), intent(in), optional :: tolerance
        real(DP), intent(in), optional :: sun_velocity(3), earth_velocity(3), moon_velocity(3)

        type(AlgebraicVector) :: target_direction, roll_vector, relative_position, relative_velocity
        type(AlgebraicVector) :: i1, i2, i3, projected
        real(DP) :: target_position(3), target_velocity(3), fallback(3), tol
        real(DP) :: b1(3), b2(3), b3(3)
        logical :: ok

        tol = DEFAULT_GEOMETRY_TOLERANCE
        if (present(tolerance)) tol = tolerance
        status = 0
        call set_message(message, '')
        call attitude%init()
        target_velocity = 0.0_DP

        select case (trim(attitude_mode))
        case ('sun')
            target_position = sun_position
            if (present(sun_velocity)) target_velocity = sun_velocity
        case ('earth')
            target_position = earth_position
            if (present(earth_velocity)) target_velocity = earth_velocity
        case ('moon')
            target_position = moon_position
            if (present(moon_velocity)) target_velocity = moon_velocity
        case default
            status = -1
            call set_message(message, 'unknown attitude mode: '//trim(attitude_mode))
            goto 900
        end select

        call normalize_vector_real(primary_axis_body, b1, ok, tol)
        if (.not. ok) then
            status = -3
            call set_message(message, 'primary body axis is zero')
            goto 900
        end if
        b2 = secondary_axis_body - dot_product(secondary_axis_body, b1)*b1
        call normalize_vector_real_inplace(b2, ok, tol)
        if (.not. ok) then
            status = -4
            call set_message(message, 'primary and secondary body axes are parallel')
            goto 900
        end if
        b3 = cross_product_real(b1, b2)

        call da_real_minus_vector(target_position, position, target_direction)
        call da_normalize_vector(target_direction, i1, ok, tol)
        if (.not. ok) then
            status = -2
            call set_message(message, 'spacecraft coincides with attitude target')
            goto 900
        end if

        select case (trim(roll_reference))
        case ('orbit_normal')
            call vec_sub(position, target_position, relative_position)
            call vec_sub(velocity, target_velocity, relative_velocity)
            call da_cross_vector(relative_position, relative_velocity, roll_vector)
        case ('inertial_x')
            call da_set_real_vector(roll_vector, [1.0_DP, 0.0_DP, 0.0_DP])
        case ('inertial_z')
            call da_set_real_vector(roll_vector, [0.0_DP, 0.0_DP, 1.0_DP])
        case ('sun')
            call da_real_minus_vector(sun_position, position, roll_vector)
        case ('earth')
            call da_real_minus_vector(earth_position, position, roll_vector)
        case ('moon')
            call da_real_minus_vector(moon_position, position, roll_vector)
        case default
            status = -5
            call set_message(message, 'unknown roll reference: '//trim(roll_reference))
            goto 900
        end select

        call da_project_perpendicular(roll_vector, i1, projected)
        call da_normalize_vector(projected, i2, ok, tol)
        if (.not. ok) then
            call choose_fallback_axis(i1%cons(), fallback)
            call da_set_real_vector(roll_vector, fallback)
            call da_project_perpendicular(roll_vector, i1, projected)
            call da_normalize_vector(projected, i2, ok, tol)
            if (.not. ok) then
                status = -6
                call set_message(message, 'unable to construct DA roll reference')
                goto 900
            end if
            status = 1
            call set_message(message, 'DA roll reference degenerate; inertial fallback used')
        end if
        call da_cross_vector(i1, i2, i3)
        call da_linear_combo3_real(i1, b1(1), i2, b2(1), i3, b3(1), &
                                   attitude%x_body_in_inertial)
        call da_linear_combo3_real(i1, b1(2), i2, b2(2), i3, b3(2), &
                                   attitude%y_body_in_inertial)
        call da_linear_combo3_real(i1, b1(3), i2, b2(3), i3, b3(3), &
                                   attitude%z_body_in_inertial)

900     continue
        call target_direction%destroy()
        call roll_vector%destroy()
        call relative_position%destroy()
        call relative_velocity%destroy()
        call i1%destroy(); call i2%destroy(); call i3%destroy(); call projected%destroy()
    end subroutine compute_pointing_attitude_da

    subroutine compute_srp_from_ephemerides_da(self, position, velocity, sun_position, &
                                                earth_position, moon_position, da_map, acceleration, &
                                                status, message, sun_velocity, earth_velocity, moon_velocity)
        class(spacecraft_geometry_type), intent(in) :: self
        type(AlgebraicVector), intent(in) :: position, velocity
        real(DP), intent(in) :: sun_position(3), earth_position(3), moon_position(3)
        type(srp_da_parameter_map_type), intent(in) :: da_map
        type(AlgebraicVector), intent(inout) :: acceleration
        integer, intent(out) :: status
        character(len=*), intent(out), optional :: message
        real(DP), intent(in), optional :: sun_velocity(3), earth_velocity(3), moon_velocity(3)

        type(da_attitude_type) :: attitude
        character(len=256) :: attitude_message, force_message
        integer :: attitude_status
        real(DP) :: sun_vel(3), earth_vel(3), moon_vel(3)

        call acceleration%init(3)
        acceleration = 0.0_DP
        sun_vel = 0.0_DP; earth_vel = 0.0_DP; moon_vel = 0.0_DP
        if (present(sun_velocity)) sun_vel = sun_velocity
        if (present(earth_velocity)) earth_vel = earth_velocity
        if (present(moon_velocity)) moon_vel = moon_velocity
        call compute_pointing_attitude_da(position, velocity, sun_position, earth_position, &
                                          moon_position, self%attitude_mode, self%roll_reference, &
                                          self%primary_axis_body, self%secondary_axis_body, attitude, &
                                          attitude_status, attitude_message, self%geometry_tolerance, &
                                          sun_vel, earth_vel, moon_vel)
        if (attitude_status < 0) then
            status = attitude_status
            call set_message(message, attitude_message)
            call attitude%destroy()
            return
        end if
        call self%compute_srp_with_attitude_da(position, sun_position, attitude, da_map, &
                                               acceleration, status, force_message)
        if (status < 0) then
            call set_message(message, force_message)
        else if (attitude_status > 0) then
            status = attitude_status
            call set_message(message, attitude_message)
        else
            call set_message(message, force_message)
        end if
        call attitude%destroy()
    end subroutine compute_srp_from_ephemerides_da

    subroutine compute_srp_with_attitude_da(self, position, sun_position, attitude, da_map, &
                                             acceleration, status, message)
        class(spacecraft_geometry_type), intent(in) :: self
        type(AlgebraicVector), intent(in) :: position
        real(DP), intent(in) :: sun_position(3)
        type(da_attitude_type), intent(in) :: attitude
        type(srp_da_parameter_map_type), intent(in) :: da_map
        type(AlgebraicVector), intent(inout) :: acceleration
        integer, intent(out) :: status
        character(len=*), intent(out), optional :: message

        real(DP), parameter :: normals_body(3,6) = reshape([ &
            1.0_DP,0.0_DP,0.0_DP, -1.0_DP,0.0_DP,0.0_DP, &
            0.0_DP,1.0_DP,0.0_DP, 0.0_DP,-1.0_DP,0.0_DP, &
            0.0_DP,0.0_DP,1.0_DP, 0.0_DP,0.0_DP,-1.0_DP], [3,6])
        type(AlgebraicVector) :: sun_vector, sun_direction, normal_inertial
        type(AlgebraicVector) :: force, total_force, array_normal_body
        type(da_attitude_type) :: working_attitude
        type(DA) :: distance, distance2, pressure, scale, scale_var, factor
        real(DP) :: areas(6), c_constant(3,3), identity_error(3,3)
        integer :: i, array_status
        logical :: ok
        character(len=256) :: local_message

        call acceleration%init(3); acceleration = 0.0_DP
        status = 0
        call self%validate(status, local_message)
        if (status < 0) then
            call set_message(message, local_message)
            return
        end if
        c_constant(:,1) = attitude%x_body_in_inertial%cons()
        c_constant(:,2) = attitude%y_body_in_inertial%cons()
        c_constant(:,3) = attitude%z_body_in_inertial%cons()
        identity_error = matmul(transpose(c_constant),c_constant)
        identity_error(1,1)=identity_error(1,1)-1.0_DP
        identity_error(2,2)=identity_error(2,2)-1.0_DP
        identity_error(3,3)=identity_error(3,3)-1.0_DP
        if (maxval(abs(identity_error)) > 1.0e-10_DP .or. &
            abs(determinant3(c_constant)-1.0_DP) > 1.0e-10_DP) then
            status = -10
            call set_message(message, 'DA attitude constant is not a proper C_I_B matrix')
            return
        end if
        call working_attitude%init()
        call da_copy_vector(attitude%x_body_in_inertial, working_attitude%x_body_in_inertial)
        call da_copy_vector(attitude%y_body_in_inertial, working_attitude%y_body_in_inertial)
        call da_copy_vector(attitude%z_body_in_inertial, working_attitude%z_body_in_inertial)
        call apply_attitude_bias_da(working_attitude, da_map)
        call da_real_minus_vector(sun_position, position, sun_vector)
        call vector_norm2_sub(sun_vector, distance)
        if (distance%cons() <= self%geometry_tolerance) then
            status = -11
            call set_message(message, 'spacecraft coincides with Sun')
            goto 900
        end if
        call vec_div(sun_vector, distance, sun_direction)
        call da_mul(distance, distance, distance2)
        call real_div_da_sub(self%pressure_1au_n_m2*AU_KM**2, distance2, pressure)
        call total_force%init(3); total_force = 0.0_DP

        areas = [self%box_dimensions_m(2)*self%box_dimensions_m(3), &
                 self%box_dimensions_m(2)*self%box_dimensions_m(3), &
                 self%box_dimensions_m(1)*self%box_dimensions_m(3), &
                 self%box_dimensions_m(1)*self%box_dimensions_m(3), &
                 self%box_dimensions_m(1)*self%box_dimensions_m(2), &
                 self%box_dimensions_m(1)*self%box_dimensions_m(2)]
        do i = 1, 6
            call da_transform_body_real(working_attitude, normals_body(:,i), normal_inertial)
            call compute_panel_force_da(sun_direction, normal_inertial, areas(i), pressure, &
                                        self%box_optical, force)
            call da_accumulate_vector(total_force, force)
        end do

        if (self%array_total_area_m2 > 0.0_DP) then
            call compute_array_normal_da(working_attitude, sun_direction, self%array_tracking_mode, &
                                         self%array_hinge_axis_body, &
                                         self%array_reference_normal_body, da_map, array_normal_body, &
                                         array_status, local_message, self%geometry_tolerance)
            if (array_status < 0) then
                status = array_status - 20
                call set_message(message, local_message)
                goto 900
            end if
            call da_transform_body_da(working_attitude, array_normal_body, normal_inertial)
            call compute_panel_force_da(sun_direction, normal_inertial, self%array_total_area_m2, &
                                        pressure, self%array_front_optical, force)
            call da_accumulate_vector(total_force, force)
            call da_negate_vector(normal_inertial, array_normal_body)
            call compute_panel_force_da(sun_direction, array_normal_body, self%array_total_area_m2, &
                                        pressure, self%array_back_optical, force)
            call da_accumulate_vector(total_force, force)
            if (array_status > 0) then
                status = array_status
                call set_message(message, local_message)
            end if
        end if

        call scale%init(); scale = 1.0_DP + da_map%global_scale_nominal
        if (da_map%global_scale_index > 0) then
            call scale_var%init_var(da_map%global_scale_index)
            scale = 1.0_DP + da_map%global_scale_nominal + da_map%global_scale_span*scale_var
        end if
        factor = scale * (1.0e-3_DP/self%mass_kg)
        call vec_mul(factor, total_force, acceleration)
        if (status == 0) call set_message(message, '')

900     continue
        call sun_vector%destroy(); call sun_direction%destroy(); call normal_inertial%destroy()
        call force%destroy(); call total_force%destroy(); call array_normal_body%destroy()
        call distance%destroy(); call distance2%destroy(); call pressure%destroy()
        call scale%destroy(); call scale_var%destroy(); call factor%destroy()
        call working_attitude%destroy()
    end subroutine compute_srp_with_attitude_da

    subroutine compute_srp_with_real_attitude_da(self, position, sun_position, c_i_b, da_map, &
                                                  acceleration, status, message)
        class(spacecraft_geometry_type), intent(in) :: self
        type(AlgebraicVector), intent(in) :: position
        real(DP), intent(in) :: sun_position(3), c_i_b(3,3)
        type(srp_da_parameter_map_type), intent(in) :: da_map
        type(AlgebraicVector), intent(inout) :: acceleration
        integer, intent(out) :: status
        character(len=*), intent(out), optional :: message
        type(da_attitude_type) :: attitude
        character(len=256) :: local_message

        call attitude%init()
        call da_set_real_vector(attitude%x_body_in_inertial, c_i_b(:,1))
        call da_set_real_vector(attitude%y_body_in_inertial, c_i_b(:,2))
        call da_set_real_vector(attitude%z_body_in_inertial, c_i_b(:,3))
        call self%compute_srp_with_attitude_da(position, sun_position, attitude, da_map, &
                                               acceleration, status, local_message)
        call set_message(message, local_message)
        call attitude%destroy()
    end subroutine compute_srp_with_real_attitude_da

    subroutine compute_panel_force_da(sun_direction, panel_normal, area_m2, pressure, optical, force)
        type(AlgebraicVector), intent(in) :: sun_direction, panel_normal
        real(DP), intent(in) :: area_m2
        type(DA), intent(in) :: pressure
        type(srp_optical_properties_type), intent(in) :: optical
        type(AlgebraicVector), intent(inout) :: force

        type(DA) :: mu, reflection, factor
        type(AlgebraicVector) :: incident_term, normal_term, bracket

        call force%init(3); force = 0.0_DP
        call vector_dot_vector_sub(panel_normal, sun_direction, mu)
        if (mu%cons() <= 0.0_DP .or. area_m2 <= 0.0_DP) goto 900
        reflection = 2.0_DP*(optical%specular_reflectivity*mu + &
                             optical%diffuse_reflectivity/3.0_DP)
        call vec_mul(1.0_DP-optical%specular_reflectivity, sun_direction, incident_term)
        call vec_mul(reflection, panel_normal, normal_term)
        call vec_add(incident_term, normal_term, bracket)
        factor = -area_m2*mu*pressure
        call vec_mul(factor, bracket, force)
900     continue
        call mu%destroy(); call reflection%destroy(); call factor%destroy()
        call incident_term%destroy(); call normal_term%destroy(); call bracket%destroy()
    end subroutine compute_panel_force_da

    subroutine compute_array_normal_da(attitude, sun_direction_inertial, tracking_mode, &
                                       hinge_axis_body, reference_normal_body, da_map, normal_body, &
                                       status, message, tolerance)
        type(da_attitude_type), intent(in) :: attitude
        type(AlgebraicVector), intent(in) :: sun_direction_inertial
        character(len=*), intent(in) :: tracking_mode
        real(DP), intent(in) :: hinge_axis_body(3), reference_normal_body(3)
        type(srp_da_parameter_map_type), intent(in) :: da_map
        type(AlgebraicVector), intent(inout) :: normal_body
        integer, intent(out) :: status
        character(len=*), intent(out), optional :: message
        real(DP), intent(in), optional :: tolerance

        type(AlgebraicVector) :: sun_body, hinge_da, reference_da, base, scaled
        type(AlgebraicVector) :: cross_term, term1, term2, term3, sum12
        type(DA) :: projection, angle, angle_var, sine, cosine, hinge_component, coefficient
        real(DP) :: hinge(3), reference(3), tol
        logical :: ok

        tol = DEFAULT_GEOMETRY_TOLERANCE
        if (present(tolerance)) tol = tolerance
        status = 0
        call set_message(message, '')
        call normalize_vector_real(hinge_axis_body, hinge, ok, tol)
        if (.not. ok) then
            status = -1
            call set_message(message, 'solar-array hinge axis is zero')
            goto 900
        end if
        reference = reference_normal_body - dot_product(reference_normal_body, hinge)*hinge
        call normalize_vector_real_inplace(reference, ok, tol)
        if (.not. ok) then
            status = -2
            call set_message(message, 'solar-array reference normal is parallel to hinge')
            goto 900
        end if
        call da_set_real_vector(hinge_da, hinge)
        call da_set_real_vector(reference_da, reference)

        select case (trim(tracking_mode))
        case ('fixed')
            call vec_mul(1.0_DP, reference_da, base)
        case ('single_axis')
            call da_transform_inertial_to_body(attitude, sun_direction_inertial, sun_body)
            call vector_dot_vector_sub(sun_body, hinge_da, projection)
            call vec_mul(projection, hinge_da, scaled)
            call vec_sub(sun_body, scaled, base)
            call da_normalize_vector_inplace(base, ok, tol)
            if (.not. ok) then
                call vec_mul(1.0_DP, reference_da, base)
                status = 1
                call set_message(message, 'DA Sun direction parallel to hinge; reference normal used')
            end if
        case default
            status = -3
            call set_message(message, 'unknown solar-array tracking mode: '//trim(tracking_mode))
            goto 900
        end select

        call angle%init(); angle = 0.0_DP
        if (da_map%array_angle_index > 0) then
            call angle_var%init_var(da_map%array_angle_index)
            angle = da_map%array_angle_span_rad*angle_var
        end if
        call da_sin_sub(angle, sine)
        call da_cos_sub(angle, cosine)
        call da_cross_vector(hinge_da, base, cross_term)
        call vec_mul(cosine, base, term1)
        call vec_mul(sine, cross_term, term2)
        call vector_dot_vector_sub(hinge_da, base, hinge_component)
        coefficient = hinge_component*(1.0_DP-cosine)
        call vec_mul(coefficient, hinge_da, term3)
        call vec_add(term1, term2, sum12)
        call vec_add(sum12, term3, normal_body)
        call da_normalize_vector_inplace(normal_body, ok, tol)
        if (.not. ok) then
            status = -4
            call set_message(message, 'computed DA solar-array normal is zero')
        end if

900     continue
        call sun_body%destroy(); call hinge_da%destroy(); call reference_da%destroy()
        call base%destroy(); call scaled%destroy(); call cross_term%destroy()
        call term1%destroy(); call term2%destroy(); call term3%destroy(); call sum12%destroy()
        call projection%destroy(); call angle%destroy(); call angle_var%destroy()
        call sine%destroy(); call cosine%destroy(); call hinge_component%destroy(); call coefficient%destroy()
    end subroutine compute_array_normal_da

    subroutine apply_attitude_bias_da(attitude, da_map)
        type(da_attitude_type), intent(inout) :: attitude
        type(srp_da_parameter_map_type), intent(in) :: da_map
        type(da_attitude_type) :: nominal
        type(AlgebraicVector) :: delta, basis, k1, k2, t1, t2v, sum12, rotated_body
        type(AlgebraicVector) :: new_columns(3)
        type(DA) :: variable, theta2, theta4, theta6, theta8, sinc_theta, cosc_theta
        integer :: axis, component

        call delta%init(3); delta = 0.0_DP
        do component = 1, 3
            if (da_map%attitude_bias_index(component) <= 0) cycle
            call variable%init_var(da_map%attitude_bias_index(component))
            delta%elements(component) = da_map%attitude_bias_span_rad(component)*variable
            call variable%destroy()
        end do
        if (maxval(abs(delta%cons())) == 0.0_DP .and. all(da_map%attitude_bias_index <= 0)) then
            call delta%destroy()
            return
        end if
        call vector_dot_vector_sub(delta, delta, theta2)
        theta4 = theta2*theta2; theta6 = theta4*theta2; theta8 = theta4*theta4
        sinc_theta = 1.0_DP - theta2/6.0_DP + theta4/120.0_DP - &
                     theta6/5040.0_DP + theta8/362880.0_DP
        cosc_theta = 0.5_DP - theta2/24.0_DP + theta4/720.0_DP - &
                     theta6/40320.0_DP + theta8/3628800.0_DP
        call nominal%init()
        call da_copy_vector(attitude%x_body_in_inertial, nominal%x_body_in_inertial)
        call da_copy_vector(attitude%y_body_in_inertial, nominal%y_body_in_inertial)
        call da_copy_vector(attitude%z_body_in_inertial, nominal%z_body_in_inertial)
        do axis = 1, 3
            call da_set_real_vector(basis, unit_basis(axis))
            call da_cross_vector(delta, basis, k1)
            call da_cross_vector(delta, k1, k2)
            call vec_mul(sinc_theta, k1, t1)
            call vec_mul(cosc_theta, k2, t2v)
            call vec_add(basis, t1, sum12)
            call vec_add(sum12, t2v, rotated_body)
            call da_transform_body_da(nominal, rotated_body, new_columns(axis))
        end do
        call da_copy_vector(new_columns(1), attitude%x_body_in_inertial)
        call da_copy_vector(new_columns(2), attitude%y_body_in_inertial)
        call da_copy_vector(new_columns(3), attitude%z_body_in_inertial)
        call nominal%destroy(); call delta%destroy(); call basis%destroy(); call k1%destroy(); call k2%destroy()
        call t1%destroy(); call t2v%destroy(); call sum12%destroy(); call rotated_body%destroy()
        do axis=1,3; call new_columns(axis)%destroy(); end do
        call theta2%destroy(); call theta4%destroy(); call theta6%destroy(); call theta8%destroy()
        call sinc_theta%destroy(); call cosc_theta%destroy(); call variable%destroy()
    end subroutine apply_attitude_bias_da

    pure function unit_basis(axis) result(vector)
        integer, intent(in) :: axis
        real(DP) :: vector(3)
        vector = 0.0_DP
        vector(axis) = 1.0_DP
    end function unit_basis

    subroutine da_set_real_vector(vector, values)
        type(AlgebraicVector), intent(inout) :: vector
        real(DP), intent(in) :: values(3)
        integer :: i
        call vector%init(3)
        do i = 1, 3
            vector%elements(i) = values(i)
        end do
    end subroutine da_set_real_vector

    subroutine da_copy_vector(source, destination)
        type(AlgebraicVector), intent(in) :: source
        type(AlgebraicVector), intent(inout) :: destination
        integer :: i
        call destination%init(3)
        do i = 1, 3
            destination%elements(i) = source%elements(i)
        end do
    end subroutine da_copy_vector

    subroutine da_real_minus_vector(real_vector, vector, result_vector)
        real(DP), intent(in) :: real_vector(3)
        type(AlgebraicVector), intent(in) :: vector
        type(AlgebraicVector), intent(inout) :: result_vector
        integer :: i
        call result_vector%init(3)
        do i = 1, 3
            result_vector%elements(i) = real_vector(i) - vector%elements(i)
        end do
    end subroutine da_real_minus_vector

    subroutine da_normalize_vector(vector, unit_vector, success, tolerance)
        type(AlgebraicVector), intent(in) :: vector
        type(AlgebraicVector), intent(inout) :: unit_vector
        logical, intent(out) :: success
        real(DP), intent(in) :: tolerance
        type(DA) :: magnitude
        integer :: i

        call vector_norm2_sub(vector, magnitude)
        success = magnitude%cons() > tolerance
        call unit_vector%init(3)
        if (success) then
            do i = 1, 3
                unit_vector%elements(i) = vector%elements(i)/magnitude
            end do
        else
            unit_vector = 0.0_DP
        end if
        call magnitude%destroy()
    end subroutine da_normalize_vector

    subroutine da_normalize_vector_inplace(vector, success, tolerance)
        type(AlgebraicVector), intent(inout) :: vector
        logical, intent(out) :: success
        real(DP), intent(in) :: tolerance
        type(AlgebraicVector) :: normalized
        call da_normalize_vector(vector, normalized, success, tolerance)
        call da_copy_vector(normalized, vector)
        call normalized%destroy()
    end subroutine da_normalize_vector_inplace

    subroutine da_project_perpendicular(vector, unit_axis, projected)
        type(AlgebraicVector), intent(in) :: vector, unit_axis
        type(AlgebraicVector), intent(inout) :: projected
        type(AlgebraicVector) :: parallel
        type(DA) :: component
        call vector_dot_vector_sub(vector, unit_axis, component)
        call vec_mul(component, unit_axis, parallel)
        call vec_sub(vector, parallel, projected)
        call parallel%destroy(); call component%destroy()
    end subroutine da_project_perpendicular

    subroutine da_cross_vector(a, b, c)
        type(AlgebraicVector), intent(in) :: a, b
        type(AlgebraicVector), intent(inout) :: c
        call c%init(3)
        c%elements(1) = a%elements(2)*b%elements(3) - a%elements(3)*b%elements(2)
        c%elements(2) = a%elements(3)*b%elements(1) - a%elements(1)*b%elements(3)
        c%elements(3) = a%elements(1)*b%elements(2) - a%elements(2)*b%elements(1)
    end subroutine da_cross_vector

    subroutine da_linear_combo3_real(v1, a1, v2, a2, v3, a3, result_vector)
        type(AlgebraicVector), intent(in) :: v1, v2, v3
        real(DP), intent(in) :: a1, a2, a3
        type(AlgebraicVector), intent(inout) :: result_vector
        type(AlgebraicVector) :: t1, t2, t3, sum12
        call vec_mul(a1, v1, t1); call vec_mul(a2, v2, t2); call vec_mul(a3, v3, t3)
        call vec_add(t1, t2, sum12); call vec_add(sum12, t3, result_vector)
        call t1%destroy(); call t2%destroy(); call t3%destroy(); call sum12%destroy()
    end subroutine da_linear_combo3_real

    subroutine da_linear_combo2_da(v1, a1, v2, a2, result_vector)
        type(AlgebraicVector), intent(in) :: v1, v2
        type(DA), intent(in) :: a1, a2
        type(AlgebraicVector), intent(inout) :: result_vector
        type(AlgebraicVector) :: t1, t2
        call vec_mul(a1, v1, t1); call vec_mul(a2, v2, t2)
        call vec_add(t1, t2, result_vector)
        call t1%destroy(); call t2%destroy()
    end subroutine da_linear_combo2_da

    subroutine da_transform_body_real(attitude, body_vector, inertial_vector)
        type(da_attitude_type), intent(in) :: attitude
        real(DP), intent(in) :: body_vector(3)
        type(AlgebraicVector), intent(inout) :: inertial_vector
        call da_linear_combo3_real(attitude%x_body_in_inertial, body_vector(1), &
                                   attitude%y_body_in_inertial, body_vector(2), &
                                   attitude%z_body_in_inertial, body_vector(3), inertial_vector)
    end subroutine da_transform_body_real

    subroutine da_transform_body_da(attitude, body_vector, inertial_vector)
        type(da_attitude_type), intent(in) :: attitude
        type(AlgebraicVector), intent(in) :: body_vector
        type(AlgebraicVector), intent(inout) :: inertial_vector
        type(AlgebraicVector) :: t1, t2, t3, sum12
        call vec_mul(body_vector%elements(1), attitude%x_body_in_inertial, t1)
        call vec_mul(body_vector%elements(2), attitude%y_body_in_inertial, t2)
        call vec_mul(body_vector%elements(3), attitude%z_body_in_inertial, t3)
        call vec_add(t1, t2, sum12); call vec_add(sum12, t3, inertial_vector)
        call t1%destroy(); call t2%destroy(); call t3%destroy(); call sum12%destroy()
    end subroutine da_transform_body_da

    subroutine da_transform_inertial_to_body(attitude, inertial_vector, body_vector)
        type(da_attitude_type), intent(in) :: attitude
        type(AlgebraicVector), intent(in) :: inertial_vector
        type(AlgebraicVector), intent(inout) :: body_vector
        call body_vector%init(3)
        call vector_dot_vector_sub(attitude%x_body_in_inertial, inertial_vector, body_vector%elements(1))
        call vector_dot_vector_sub(attitude%y_body_in_inertial, inertial_vector, body_vector%elements(2))
        call vector_dot_vector_sub(attitude%z_body_in_inertial, inertial_vector, body_vector%elements(3))
    end subroutine da_transform_inertial_to_body

    subroutine da_negate_vector(vector, negative)
        type(AlgebraicVector), intent(in) :: vector
        type(AlgebraicVector), intent(inout) :: negative
        call vec_mul(-1.0_DP, vector, negative)
    end subroutine da_negate_vector

    subroutine da_accumulate_vector(total, addend)
        type(AlgebraicVector), intent(inout) :: total
        type(AlgebraicVector), intent(in) :: addend
        type(AlgebraicVector) :: temporary
        call vec_add(total, addend, temporary)
        call da_copy_vector(temporary, total)
        call temporary%destroy()
    end subroutine da_accumulate_vector

    subroutine initialize_geometry_from_config(self, cfg, status, message)
        class(spacecraft_geometry_type), intent(inout) :: self
        type(config_params), intent(in) :: cfg
        integer, intent(out) :: status
        character(len=*), intent(out), optional :: message

        self%mass_kg = cfg%srp_mass_kg
        self%box_dimensions_m = cfg%srp_box_dimensions_m
        call set_optical_from_array(self%box_optical, cfg%srp_box_optical)
        self%array_total_area_m2 = cfg%srp_array_total_area_m2
        self%array_tracking_mode = cfg%srp_array_tracking_mode
        self%array_hinge_axis_body = cfg%srp_array_hinge_axis_body
        self%array_reference_normal_body = cfg%srp_array_reference_normal_body
        call set_optical_from_array(self%array_front_optical, cfg%srp_array_front_optical)
        call set_optical_from_array(self%array_back_optical, cfg%srp_array_back_optical)
        self%attitude_mode = cfg%srp_attitude_mode
        self%roll_reference = cfg%srp_roll_reference
        self%primary_axis_body = cfg%srp_primary_axis_body
        self%secondary_axis_body = cfg%srp_secondary_axis_body
        if (cfg%srp_pressure_1au_n_m2 == -1.0_DP) then
            self%pressure_1au_n_m2 = DEFAULT_SOLAR_PRESSURE_1AU
        else
            self%pressure_1au_n_m2 = cfg%srp_pressure_1au_n_m2
        end if
        self%geometry_tolerance = cfg%srp_geometry_tolerance
        call self%validate(status, message)
    end subroutine initialize_geometry_from_config

    subroutine validate_spacecraft_geometry(self, status, message)
        class(spacecraft_geometry_type), intent(in) :: self
        integer, intent(out) :: status
        character(len=*), intent(out), optional :: message
        real(DP) :: primary_norm, secondary_norm, hinge_norm, reference_norm

        status = 0
        call set_message(message, '')
        if (self%mass_kg <= 0.0_DP) then
            status = -1
            call set_message(message, 'spacecraft mass must be positive')
        else if (any(self%box_dimensions_m <= 0.0_DP)) then
            status = -2
            call set_message(message, 'box dimensions must be positive')
        else if (self%array_total_area_m2 < 0.0_DP) then
            status = -3
            call set_message(message, 'solar-array area cannot be negative')
        else if (self%pressure_1au_n_m2 <= 0.0_DP) then
            status = -4
            call set_message(message, 'solar pressure must be positive')
        else if (self%geometry_tolerance <= 0.0_DP) then
            status = -5
            call set_message(message, 'geometry tolerance must be positive')
        else if (trim(self%attitude_mode) /= 'sun' .and. trim(self%attitude_mode) /= 'earth' .and. &
                 trim(self%attitude_mode) /= 'moon') then
            status = -6
            call set_message(message, 'attitude mode must be sun, earth, or moon')
        else if (trim(self%roll_reference) /= 'orbit_normal' .and. &
                 trim(self%roll_reference) /= 'inertial_x' .and. &
                 trim(self%roll_reference) /= 'inertial_z' .and. &
                 trim(self%roll_reference) /= 'sun' .and. trim(self%roll_reference) /= 'earth' .and. &
                 trim(self%roll_reference) /= 'moon') then
            status = -7
            call set_message(message, 'invalid roll reference')
        else if (.not. valid_optical_properties(self%box_optical)) then
            status = -8
            call set_message(message, 'invalid box optical properties')
        else if (trim(self%array_tracking_mode) /= 'fixed' .and. &
                 trim(self%array_tracking_mode) /= 'single_axis') then
            status = -9
            call set_message(message, 'array tracking mode must be fixed or single_axis')
        end if
        if (status < 0) return
        primary_norm = sqrt(sum(self%primary_axis_body**2))
        secondary_norm = sqrt(sum(self%secondary_axis_body**2))
        if (primary_norm <= self%geometry_tolerance .or. secondary_norm <= self%geometry_tolerance) then
            status = -10
            call set_message(message, 'primary and secondary body axes must be nonzero')
        else if (abs(dot_product(self%primary_axis_body, self%secondary_axis_body)) >= &
                 (1.0_DP-self%geometry_tolerance)*primary_norm*secondary_norm) then
            status = -11
            call set_message(message, 'primary and secondary body axes cannot be parallel')
        end if
        if (status < 0 .or. self%array_total_area_m2 == 0.0_DP) return
        hinge_norm = sqrt(sum(self%array_hinge_axis_body**2))
        reference_norm = sqrt(sum(self%array_reference_normal_body**2))
        if (hinge_norm <= self%geometry_tolerance .or. reference_norm <= self%geometry_tolerance) then
            status = -12
            call set_message(message, 'array hinge and reference normal must be nonzero')
        else if (abs(dot_product(self%array_hinge_axis_body, self%array_reference_normal_body)) >= &
                 (1.0_DP-self%geometry_tolerance)*hinge_norm*reference_norm) then
            status = -13
            call set_message(message, 'array hinge and reference normal cannot be parallel')
        else if (.not. valid_optical_properties(self%array_front_optical) .or. &
                 .not. valid_optical_properties(self%array_back_optical)) then
            status = -14
            call set_message(message, 'invalid solar-array optical properties')
        end if
    end subroutine validate_spacecraft_geometry

    subroutine compute_srp_from_ephemerides_real(self, position, velocity, sun_position, &
                                                  earth_position, moon_position, acceleration, &
                                                  status, message, global_scale_error, array_angle_bias_rad, &
                                                  sun_velocity, earth_velocity, moon_velocity)
        class(spacecraft_geometry_type), intent(in) :: self
        real(DP), intent(in) :: position(3), velocity(3)
        real(DP), intent(in) :: sun_position(3), earth_position(3), moon_position(3)
        real(DP), intent(out) :: acceleration(3)
        integer, intent(out) :: status
        character(len=*), intent(out), optional :: message
        real(DP), intent(in), optional :: global_scale_error, array_angle_bias_rad
        real(DP), intent(in), optional :: sun_velocity(3), earth_velocity(3), moon_velocity(3)

        real(DP) :: c_i_b(3,3)
        integer :: attitude_status
        character(len=256) :: attitude_message, force_message
        real(DP) :: sun_vel(3), earth_vel(3), moon_vel(3)

        acceleration = 0.0_DP
        sun_vel = 0.0_DP; earth_vel = 0.0_DP; moon_vel = 0.0_DP
        if (present(sun_velocity)) sun_vel = sun_velocity
        if (present(earth_velocity)) earth_vel = earth_velocity
        if (present(moon_velocity)) moon_vel = moon_velocity
        call compute_pointing_attitude_real(position, velocity, sun_position, earth_position, &
                                             moon_position, self%attitude_mode, self%roll_reference, &
                                             self%primary_axis_body, self%secondary_axis_body, c_i_b, &
                                             attitude_status, attitude_message, self%geometry_tolerance, &
                                             sun_vel, earth_vel, moon_vel)
        if (attitude_status < 0) then
            status = attitude_status
            call set_message(message, attitude_message)
            return
        end if

        call self%compute_srp_with_attitude_real(position, sun_position, c_i_b, acceleration, &
                                                 status, force_message, global_scale_error, &
                                                 array_angle_bias_rad)
        if (status < 0) then
            call set_message(message, force_message)
        else if (attitude_status > 0) then
            status = attitude_status
            call set_message(message, attitude_message)
        else
            call set_message(message, force_message)
        end if
    end subroutine compute_srp_from_ephemerides_real

    subroutine compute_srp_with_attitude_real(self, position, sun_position, c_i_b, acceleration, &
                                               status, message, global_scale_error, array_angle_bias_rad)
        class(spacecraft_geometry_type), intent(in) :: self
        real(DP), intent(in) :: position(3), sun_position(3), c_i_b(3,3)
        real(DP), intent(out) :: acceleration(3)
        integer, intent(out) :: status
        character(len=*), intent(out), optional :: message
        real(DP), intent(in), optional :: global_scale_error, array_angle_bias_rad

        real(DP), parameter :: normals_body(3,6) = reshape([ &
            1.0_DP, 0.0_DP, 0.0_DP, -1.0_DP, 0.0_DP, 0.0_DP, &
            0.0_DP, 1.0_DP, 0.0_DP,  0.0_DP,-1.0_DP, 0.0_DP, &
            0.0_DP, 0.0_DP, 1.0_DP,  0.0_DP, 0.0_DP,-1.0_DP], [3,6])
        real(DP) :: areas(6), sun_vector(3), sun_direction(3), solar_distance
        real(DP) :: pressure, force(3), total_force(3), normal_inertial(3), array_normal_body(3)
        real(DP) :: scale_factor, array_angle, identity_error(3,3)
        integer :: i, array_status
        logical :: ok
        character(len=256) :: local_message

        acceleration = 0.0_DP
        status = 0
        call set_message(message, '')
        call self%validate(status, local_message)
        if (status < 0) then
            call set_message(message, local_message)
            return
        end if

        identity_error = matmul(transpose(c_i_b), c_i_b)
        identity_error(1,1) = identity_error(1,1) - 1.0_DP
        identity_error(2,2) = identity_error(2,2) - 1.0_DP
        identity_error(3,3) = identity_error(3,3) - 1.0_DP
        if (maxval(abs(identity_error)) > 1.0e-10_DP .or. &
            abs(determinant3(c_i_b) - 1.0_DP) > 1.0e-10_DP) then
            status = -10
            call set_message(message, 'C_I_B is not a proper orthonormal attitude matrix')
            return
        end if

        sun_vector = sun_position - position
        solar_distance = sqrt(dot_product(sun_vector, sun_vector))
        call normalize_vector_real(sun_vector, sun_direction, ok, self%geometry_tolerance)
        if (.not. ok) then
            status = -11
            call set_message(message, 'spacecraft coincides with Sun')
            return
        end if
        pressure = self%pressure_1au_n_m2 * (AU_KM / solar_distance)**2
        total_force = 0.0_DP

        areas = [self%box_dimensions_m(2)*self%box_dimensions_m(3), &
                 self%box_dimensions_m(2)*self%box_dimensions_m(3), &
                 self%box_dimensions_m(1)*self%box_dimensions_m(3), &
                 self%box_dimensions_m(1)*self%box_dimensions_m(3), &
                 self%box_dimensions_m(1)*self%box_dimensions_m(2), &
                 self%box_dimensions_m(1)*self%box_dimensions_m(2)]
        do i = 1, 6
            normal_inertial = matmul(c_i_b, normals_body(:,i))
            call compute_panel_force_real(sun_direction, normal_inertial, areas(i), pressure, &
                                          self%box_optical, force)
            total_force = total_force + force
        end do

        if (self%array_total_area_m2 > 0.0_DP) then
            array_angle = 0.0_DP
            if (present(array_angle_bias_rad)) array_angle = array_angle_bias_rad
            call compute_array_normal_real(c_i_b, sun_direction, self%array_tracking_mode, &
                                           self%array_hinge_axis_body, &
                                           self%array_reference_normal_body, array_angle, &
                                           array_normal_body, array_status, local_message, &
                                           self%geometry_tolerance)
            if (array_status < 0) then
                status = array_status - 20
                call set_message(message, local_message)
                return
            end if
            normal_inertial = matmul(c_i_b, array_normal_body)
            call compute_panel_force_real(sun_direction, normal_inertial, self%array_total_area_m2, &
                                          pressure, self%array_front_optical, force)
            total_force = total_force + force
            call compute_panel_force_real(sun_direction, -normal_inertial, self%array_total_area_m2, &
                                          pressure, self%array_back_optical, force)
            total_force = total_force + force
            if (array_status > 0) then
                status = array_status
                call set_message(message, local_message)
            end if
        end if

        scale_factor = 1.0_DP
        if (present(global_scale_error)) scale_factor = 1.0_DP + global_scale_error
        acceleration = scale_factor * total_force / self%mass_kg * 1.0e-3_DP
    end subroutine compute_srp_with_attitude_real

    subroutine compute_pointing_attitude_real(position, velocity, sun_position, earth_position, &
                                               moon_position, attitude_mode, roll_reference, &
                                               primary_axis_body, secondary_axis_body, c_i_b, &
                                               status, message, tolerance, sun_velocity, earth_velocity, &
                                               moon_velocity)
        real(DP), intent(in) :: position(3), velocity(3)
        real(DP), intent(in) :: sun_position(3), earth_position(3), moon_position(3)
        character(len=*), intent(in) :: attitude_mode, roll_reference
        real(DP), intent(in) :: primary_axis_body(3), secondary_axis_body(3)
        real(DP), intent(out) :: c_i_b(3,3)
        integer, intent(out) :: status
        character(len=*), intent(out), optional :: message
        real(DP), intent(in), optional :: tolerance
        real(DP), intent(in), optional :: sun_velocity(3), earth_velocity(3), moon_velocity(3)

        real(DP) :: tol, target_position(3), target_velocity(3), target_direction(3), roll_vector(3)
        real(DP) :: b1(3), b2(3), b3(3), i1(3), i2(3), i3(3)
        logical :: ok

        tol = DEFAULT_GEOMETRY_TOLERANCE
        if (present(tolerance)) tol = tolerance
        c_i_b = 0.0_DP
        status = 0
        call set_message(message, '')
        target_velocity = 0.0_DP

        select case (trim(attitude_mode))
        case ('sun')
            target_position = sun_position
            if (present(sun_velocity)) target_velocity = sun_velocity
        case ('earth')
            target_position = earth_position
            if (present(earth_velocity)) target_velocity = earth_velocity
        case ('moon')
            target_position = moon_position
            if (present(moon_velocity)) target_velocity = moon_velocity
        case default
            status = -1
            call set_message(message, 'unknown attitude mode: '//trim(attitude_mode))
            return
        end select

        target_direction = target_position - position
        call normalize_vector_real(target_direction, i1, ok, tol)
        if (.not. ok) then
            status = -2
            call set_message(message, 'spacecraft coincides with attitude target')
            return
        end if

        call normalize_vector_real(primary_axis_body, b1, ok, tol)
        if (.not. ok) then
            status = -3
            call set_message(message, 'primary body axis is zero')
            return
        end if
        b2 = secondary_axis_body - dot_product(secondary_axis_body, b1) * b1
        call normalize_vector_real_inplace(b2, ok, tol)
        if (.not. ok) then
            status = -4
            call set_message(message, 'primary and secondary body axes are parallel')
            return
        end if
        b3 = cross_product_real(b1, b2)

        select case (trim(roll_reference))
        case ('orbit_normal')
            roll_vector = cross_product_real(position - target_position, velocity - target_velocity)
        case ('inertial_x')
            roll_vector = [1.0_DP, 0.0_DP, 0.0_DP]
        case ('inertial_z')
            roll_vector = [0.0_DP, 0.0_DP, 1.0_DP]
        case ('sun')
            roll_vector = sun_position - position
        case ('earth')
            roll_vector = earth_position - position
        case ('moon')
            roll_vector = moon_position - position
        case default
            status = -5
            call set_message(message, 'unknown roll reference: '//trim(roll_reference))
            return
        end select

        i2 = roll_vector - dot_product(roll_vector, i1) * i1
        call normalize_vector_real_inplace(i2, ok, tol)
        if (.not. ok) then
            call choose_fallback_axis(i1, roll_vector)
            i2 = roll_vector - dot_product(roll_vector, i1) * i1
            call normalize_vector_real_inplace(i2, ok, tol)
            if (.not. ok) then
                status = -6
                call set_message(message, 'unable to construct roll reference')
                return
            end if
            status = 1
            call set_message(message, 'roll reference degenerate; inertial fallback used')
        end if
        i3 = cross_product_real(i1, i2)

        c_i_b = outer_product(i1, b1) + outer_product(i2, b2) + outer_product(i3, b3)
    end subroutine compute_pointing_attitude_real

    subroutine compute_array_normal_real(c_i_b, sun_direction_inertial, tracking_mode, &
                                         hinge_axis_body, reference_normal_body, angle_bias_rad, &
                                         normal_body, status, message, tolerance)
        real(DP), intent(in) :: c_i_b(3,3), sun_direction_inertial(3)
        character(len=*), intent(in) :: tracking_mode
        real(DP), intent(in) :: hinge_axis_body(3), reference_normal_body(3)
        real(DP), intent(in) :: angle_bias_rad
        real(DP), intent(out) :: normal_body(3)
        integer, intent(out) :: status
        character(len=*), intent(out), optional :: message
        real(DP), intent(in), optional :: tolerance

        real(DP) :: tol, hinge(3), reference(3), sun_body(3), base_normal(3)
        logical :: ok

        tol = DEFAULT_GEOMETRY_TOLERANCE
        if (present(tolerance)) tol = tolerance
        normal_body = 0.0_DP
        status = 0
        call set_message(message, '')

        call normalize_vector_real(hinge_axis_body, hinge, ok, tol)
        if (.not. ok) then
            status = -1
            call set_message(message, 'solar-array hinge axis is zero')
            return
        end if
        reference = reference_normal_body - dot_product(reference_normal_body, hinge) * hinge
        call normalize_vector_real_inplace(reference, ok, tol)
        if (.not. ok) then
            status = -2
            call set_message(message, 'solar-array reference normal is parallel to hinge')
            return
        end if

        select case (trim(tracking_mode))
        case ('fixed')
            base_normal = reference
        case ('single_axis')
            sun_body = matmul(transpose(c_i_b), sun_direction_inertial)
            base_normal = sun_body - dot_product(sun_body, hinge) * hinge
            call normalize_vector_real_inplace(base_normal, ok, tol)
            if (.not. ok) then
                base_normal = reference
                status = 1
                call set_message(message, 'Sun direction parallel to hinge; reference normal used')
            end if
        case default
            status = -3
            call set_message(message, 'unknown solar-array tracking mode: '//trim(tracking_mode))
            return
        end select

        normal_body = rotate_about_axis(base_normal, hinge, angle_bias_rad)
        call normalize_vector_real_inplace(normal_body, ok, tol)
        if (.not. ok) then
            status = -4
            call set_message(message, 'computed solar-array normal is zero')
        end if
    end subroutine compute_array_normal_real

    subroutine compute_panel_force_real(sun_direction, panel_normal, area_m2, pressure_n_m2, &
                                        optical, force_newton)
        real(DP), intent(in) :: sun_direction(3), panel_normal(3)
        real(DP), intent(in) :: area_m2, pressure_n_m2
        type(srp_optical_properties_type), intent(in) :: optical
        real(DP), intent(out) :: force_newton(3)

        real(DP) :: mu

        force_newton = 0.0_DP
        mu = dot_product(panel_normal, sun_direction)
        if (mu <= 0.0_DP .or. area_m2 <= 0.0_DP .or. pressure_n_m2 <= 0.0_DP) return

        force_newton = -pressure_n_m2 * area_m2 * mu * &
            ((1.0_DP - optical%specular_reflectivity) * sun_direction + &
             2.0_DP * (optical%specular_reflectivity * mu + &
                       optical%diffuse_reflectivity / 3.0_DP) * panel_normal)
    end subroutine compute_panel_force_real

    subroutine normalize_vector_real(vector, unit_vector, success, tolerance)
        real(DP), intent(in) :: vector(3)
        real(DP), intent(out) :: unit_vector(3)
        logical, intent(out) :: success
        real(DP), intent(in), optional :: tolerance

        real(DP) :: norm_value, tol

        tol = DEFAULT_GEOMETRY_TOLERANCE
        if (present(tolerance)) tol = tolerance
        norm_value = sqrt(dot_product(vector, vector))
        success = norm_value > tol
        if (success) then
            unit_vector = vector / norm_value
        else
            unit_vector = 0.0_DP
        end if
    end subroutine normalize_vector_real

    subroutine normalize_vector_real_inplace(vector, success, tolerance)
        real(DP), intent(inout) :: vector(3)
        logical, intent(out) :: success
        real(DP), intent(in) :: tolerance
        real(DP) :: magnitude
        magnitude = sqrt(dot_product(vector, vector))
        success = magnitude > tolerance
        if (success) then
            vector = vector/magnitude
        else
            vector = 0.0_DP
        end if
    end subroutine normalize_vector_real_inplace

    pure function cross_product_real(a, b) result(c)
        real(DP), intent(in) :: a(3), b(3)
        real(DP) :: c(3)

        c = [a(2)*b(3) - a(3)*b(2), &
             a(3)*b(1) - a(1)*b(3), &
             a(1)*b(2) - a(2)*b(1)]
    end function cross_product_real

    pure function outer_product(a, b) result(matrix)
        real(DP), intent(in) :: a(3), b(3)
        real(DP) :: matrix(3,3)
        integer :: i, j

        do j = 1, 3
            do i = 1, 3
                matrix(i,j) = a(i) * b(j)
            end do
        end do
    end function outer_product

    pure function rotate_about_axis(vector, axis, angle) result(rotated)
        real(DP), intent(in) :: vector(3), axis(3), angle
        real(DP) :: rotated(3)

        rotated = vector * cos(angle) + cross_product_real(axis, vector) * sin(angle) + &
                  axis * dot_product(axis, vector) * (1.0_DP - cos(angle))
    end function rotate_about_axis

    subroutine choose_fallback_axis(primary, fallback)
        real(DP), intent(in) :: primary(3)
        real(DP), intent(out) :: fallback(3)
        integer :: index_min

        index_min = minloc(abs(primary), dim=1)
        fallback = 0.0_DP
        fallback(index_min) = 1.0_DP
    end subroutine choose_fallback_axis

    subroutine set_message(message, text)
        character(len=*), intent(out), optional :: message
        character(len=*), intent(in) :: text
        if (present(message)) message = text
    end subroutine set_message

    subroutine set_optical_from_array(optical, values)
        type(srp_optical_properties_type), intent(out) :: optical
        real(DP), intent(in) :: values(3)
        optical%absorptivity = values(1)
        optical%specular_reflectivity = values(2)
        optical%diffuse_reflectivity = values(3)
    end subroutine set_optical_from_array

    logical function valid_optical_properties(optical)
        type(srp_optical_properties_type), intent(in) :: optical
        real(DP) :: values(3)
        values = [optical%absorptivity, optical%specular_reflectivity, optical%diffuse_reflectivity]
        valid_optical_properties = all(values >= 0.0_DP) .and. all(values <= 1.0_DP) .and. &
                                   abs(sum(values)-1.0_DP) <= 1.0e-12_DP
    end function valid_optical_properties

    pure real(DP) function determinant3(matrix)
        real(DP), intent(in) :: matrix(3,3)
        determinant3 = matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2)) &
                     - matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1)) &
                     + matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
    end function determinant3

end module pod_spacecraft_geometry
