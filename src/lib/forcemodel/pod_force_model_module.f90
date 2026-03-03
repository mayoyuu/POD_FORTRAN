module pod_force_model_module
    use pod_global, only: DP
    use pod_config, only: config
    use pod_spice, only: get_sun_position, get_moon_position
    
    implicit none
    

    
    ! 力模型配置
    logical :: use_j2_perturbation = .true.
    logical :: use_nspheric_perturbation = .true.
    logical :: use_atmospheric_drag = .false.
    logical :: use_solar_radiation_pressure = .false.
    logical :: use_third_body_perturbation = .false.
    
    ! 地球物理常数
    real(DP), parameter :: EARTH_RADIUS = 6378.137_DP  ! km
    real(DP), parameter :: J2_COEFFICIENT = 1.08263e-3_DP
    real(DP), parameter :: EARTH_MU = 398600.4418_DP  ! km³/s²
    
    ! 太阳常数
    real(DP), parameter :: SOLAR_CONSTANT = 1367.0_DP  ! W/m²
    real(DP), parameter :: SPEED_OF_LIGHT = 299792458.0_DP  ! m/s
    
contains

    subroutine compute_acceleration(position, velocity, time, acceleration)
        real(DP), dimension(3), intent(in) :: position, velocity
        real(DP), intent(in) :: time
        real(DP), dimension(3), intent(out) :: acceleration
        
        real(DP), dimension(3) :: acc_central, acc_j2, acc_drag, acc_srp, acc_third_body
        real(DP), dimension(3) :: acc_total
        
        ! 中心引力加速度
        call compute_central_gravity(position, acc_central)
        
        ! J2摄动力加速度
        if (use_j2_perturbation) then
            call compute_j2_perturbation(position, acc_j2)
        else
            acc_j2 = 0.0_DP
        end if
        
        ! 大气阻力加速度
        if (use_atmospheric_drag) then
            call compute_atmospheric_drag(position, velocity, acc_drag)
        else
            acc_drag = 0.0_DP
        end if
        
        ! 太阳辐射压加速度
        if (use_solar_radiation_pressure) then
            call compute_solar_radiation_pressure(position, acc_srp)
        else
            acc_srp = 0.0_DP
        end if
        
        ! 第三体摄动力加速度
        if (use_third_body_perturbation) then
            call compute_third_body_perturbation(position, time, acc_third_body)
        else
            acc_third_body = 0.0_DP
        end if
        
        ! 总加速度
        acceleration = acc_central + acc_j2 + acc_drag + acc_srp + acc_third_body
    end subroutine compute_acceleration
    
    subroutine compute_central_gravity(position, acceleration)
        real(DP), dimension(3), intent(in) :: position
        real(DP), dimension(3), intent(out) :: acceleration
        
        real(DP) :: r_mag, r_mag3
        
        r_mag = sqrt(sum(position**2))
        r_mag3 = r_mag**3
        
        ! 中心引力: a = -μ * r / |r|³
        acceleration = -config%gravitational_constant * position / r_mag3
    end subroutine compute_central_gravity
    
    subroutine compute_j2_perturbation(position, acceleration)
        real(DP), dimension(3), intent(in) :: position
        real(DP), dimension(3), intent(out) :: acceleration
        
        real(DP) :: r_mag, r_mag2, r_mag5, z_r, z_r2
        real(DP) :: j2_factor, j2_x, j2_y, j2_z
        
        r_mag = sqrt(sum(position**2))
        r_mag2 = r_mag**2
        r_mag5 = r_mag**5
        
        z_r = position(3) / r_mag
        z_r2 = z_r**2
        
        ! J2摄动力系数
        j2_factor = 1.5_DP * J2_COEFFICIENT * config%gravitational_constant * EARTH_RADIUS**2 / r_mag5
        
        ! J2摄动力分量
        j2_x = j2_factor * position(1) * (5.0_DP * z_r2 - 1.0_DP)
        j2_y = j2_factor * position(2) * (5.0_DP * z_r2 - 1.0_DP)
        j2_z = j2_factor * position(3) * (5.0_DP * z_r2 - 3.0_DP)
        
        acceleration = [j2_x, j2_y, j2_z]
    end subroutine compute_j2_perturbation
    
    subroutine compute_atmospheric_drag(position, velocity, acceleration)
        real(DP), dimension(3), intent(in) :: position, velocity
        real(DP), dimension(3), intent(out) :: acceleration
        
        real(DP) :: altitude, density, velocity_mag, drag_coefficient, area_mass_ratio
        real(DP) :: drag_factor
        
        ! 计算高度
        altitude = sqrt(sum(position**2)) - EARTH_RADIUS
        
        ! 计算大气密度 (简化模型)
        call compute_atmospheric_density(altitude, density)
        
        ! 计算相对速度大小
        velocity_mag = sqrt(sum(velocity**2))
        
        ! 阻力系数和面积质量比 (需要根据具体卫星参数设置)
        drag_coefficient = 2.2_DP  ! 典型值
        area_mass_ratio = 0.01_DP  ! m²/kg
        
        ! 阻力因子
        drag_factor = -0.5_DP * density * drag_coefficient * area_mass_ratio * velocity_mag
        
        ! 阻力加速度
        acceleration = drag_factor * velocity
    end subroutine compute_atmospheric_drag
    
    subroutine compute_atmospheric_density(altitude, density)
        real(DP), intent(in) :: altitude
        real(DP), intent(out) :: density
        
        real(DP), parameter :: H0 = 8.5_DP  ! km
        real(DP), parameter :: RHO0 = 1.225e-3_DP  ! kg/km³ (海平面密度)
        
        ! 指数大气密度模型
        if (altitude < 100.0_DP) then
            density = RHO0 * exp(-altitude / H0)
        else
            density = 0.0_DP  ! 高空大气密度很小
        end if
    end subroutine compute_atmospheric_density
    
    subroutine compute_solar_radiation_pressure(position, acceleration)
        real(DP), dimension(3), intent(in) :: position
        real(DP), dimension(3), intent(out) :: acceleration
        
        real(DP) :: solar_distance, solar_factor, reflectivity, area_mass_ratio
        real(DP), dimension(3) :: solar_direction
        
        ! 太阳方向 (简化处理，假设太阳在X轴正方向)
        solar_direction = [1.0_DP, 0.0_DP, 0.0_DP]
        
        ! 太阳距离 (简化处理)
        solar_distance = 1.496e8_DP  ! km (1 AU)
        
        ! 太阳辐射压因子
        reflectivity = 1.0_DP  ! 完全反射
        area_mass_ratio = 0.01_DP  ! m²/kg
        
        solar_factor = SOLAR_CONSTANT / (SPEED_OF_LIGHT * 1000.0_DP) * &
                      reflectivity * area_mass_ratio / (solar_distance**2)
        
        ! 太阳辐射压加速度
        acceleration = solar_factor * solar_direction
    end subroutine compute_solar_radiation_pressure
    
    subroutine compute_third_body_perturbation(position, time, acceleration)
        real(DP), dimension(3), intent(in) :: position
        real(DP), intent(in) :: time
        real(DP), dimension(3), intent(out) :: acceleration
        
        real(DP), dimension(3) :: moon_position, sun_position, moon_velocity, sun_velocity
        real(DP), dimension(3) :: acc_moon, acc_sun
        
        ! 获取太阳位置
        call get_sun_position(time, 'EARTH', sun_position, sun_velocity)
        
        ! 获取月球位置
        call get_moon_position(time, 'EARTH', moon_position, moon_velocity)
        
        ! 月球摄动力
        call compute_third_body_gravity(position, moon_position, 4902.8_DP, acc_moon)
        
        ! 太阳摄动力
        call compute_third_body_gravity(position, sun_position, 1.327e11_DP, acc_sun)
        
        ! 总第三体摄动力
        acceleration = acc_moon + acc_sun
    end subroutine compute_third_body_perturbation
    
    subroutine compute_third_body_gravity(position, body_position, body_mu, acceleration)
        real(DP), dimension(3), intent(in) :: position, body_position
        real(DP), intent(in) :: body_mu
        real(DP), dimension(3), intent(out) :: acceleration
        
        real(DP), dimension(3) :: r_sat, r_body, r_rel
        real(DP) :: r_sat_mag, r_body_mag, r_rel_mag
        real(DP) :: r_sat_mag3, r_body_mag3, r_rel_mag3
        
        ! 卫星到地心的向量
        r_sat = position
        r_sat_mag = sqrt(sum(r_sat**2))
        r_sat_mag3 = r_sat_mag**3
        
        ! 第三体到地心的向量
        r_body = body_position
        r_body_mag = sqrt(sum(r_body**2))
        r_body_mag3 = r_body_mag**3
        
        ! 卫星到第三体的向量
        r_rel = r_body - r_sat
        r_rel_mag = sqrt(sum(r_rel**2))
        r_rel_mag3 = r_rel_mag**3
        
        ! 第三体摄动力: a = μ_body * (r_rel/|r_rel|³ - r_body/|r_body|³)
        acceleration = body_mu * (r_rel / r_rel_mag3 - r_body / r_body_mag3)
    end subroutine compute_third_body_gravity
    
    subroutine set_force_model_options(use_j2, use_drag, use_srp, use_third_body)
        logical, intent(in) :: use_j2, use_drag, use_srp, use_third_body
        
        use_j2_perturbation = use_j2
        use_atmospheric_drag = use_drag
        use_solar_radiation_pressure = use_srp
        use_third_body_perturbation = use_third_body
    end subroutine set_force_model_options
    
    subroutine print_force_model_status()
        write(*, *) '=== 力模型状态 ==='
        write(*, *) '中心引力: 启用'
        write(*, *) 'J2摄动力: ', merge('启用', '禁用', use_j2_perturbation)
        write(*, *) '大气阻力: ', merge('启用', '禁用', use_atmospheric_drag)
        write(*, *) '太阳辐射压: ', merge('启用', '禁用', use_solar_radiation_pressure)
        write(*, *) '第三体摄动力: ', merge('启用', '禁用', use_third_body_perturbation)
    end subroutine print_force_model_status

end module pod_force_model_module
