module cat_measurement_model_module
    use cat_global, only: DP, MAX_STRING_LEN
    use cat_frame_module, only: eci_to_geodetic, geodetic_to_eci
    
    implicit none
    
    ! 观测站类型
    type observation_station
        character(len=MAX_STRING_LEN) :: name
        real(DP) :: latitude, longitude, altitude
        real(DP), dimension(3) :: eci_position
        character(len=MAX_STRING_LEN) :: station_type  ! 'RADAR', 'OPTICAL', 'GPS'
    end type observation_station
    
contains

    subroutine compute_measurement(state, time, measurement_type, measurement)
        real(DP), dimension(6), intent(in) :: state
        real(DP), intent(in) :: time
        character(len=*), intent(in) :: measurement_type
        real(DP), dimension(3), intent(out) :: measurement
        
        real(DP), dimension(3) :: position, velocity
        type(observation_station) :: station
        
        ! 提取位置和速度
        position = state(1:3)
        velocity = state(4:6)
        
        ! 设置默认观测站 (简化处理)
        call set_default_station(station)
        
        select case (trim(measurement_type))
            case ('RADAR')
                call compute_radar_measurement(position, station, measurement)
            case ('OPTICAL')
                call compute_optical_measurement(position, station, measurement)
            case ('GPS')
                call compute_gps_measurement(position, measurement)
            case default
                write(*, *) '警告: 未知的测量类型: ', trim(measurement_type)
                measurement = 0.0_DP
        end select
    end subroutine compute_measurement
    
    subroutine compute_measurement_jacobian(state, time, measurement_type, jacobian)
        real(DP), dimension(6), intent(in) :: state
        real(DP), intent(in) :: time
        character(len=*), intent(in) :: measurement_type
        real(DP), dimension(3,6), intent(out) :: jacobian
        
        real(DP), dimension(3) :: position, velocity
        real(DP) :: h, perturbation
        real(DP), dimension(6) :: state_perturbed
        real(DP), dimension(3) :: measurement_original, measurement_perturbed
        integer :: i
        
        ! 提取位置和速度
        position = state(1:3)
        velocity = state(4:6)
        
        ! 计算原始测量值
        call compute_measurement(state, time, measurement_type, measurement_original)
        
        ! 数值微分计算雅可比矩阵
        h = 1.0e-6_DP  ! 扰动大小
        
        do i = 1, 6
            state_perturbed = state
            state_perturbed(i) = state_perturbed(i) + h
            
            call compute_measurement(state_perturbed, time, measurement_type, measurement_perturbed)
            
            jacobian(:, i) = (measurement_perturbed - measurement_original) / h
        end do
    end subroutine compute_measurement_jacobian
    
    subroutine compute_radar_measurement(position, station, measurement)
        real(DP), dimension(3), intent(in) :: position
        type(observation_station), intent(in) :: station
        real(DP), dimension(3), intent(out) :: measurement
        
        real(DP), dimension(3) :: range_vector, range_unit
        real(DP) :: range_mag, azimuth, elevation
        
        ! 计算距离向量
        range_vector = position - station%eci_position
        range_mag = sqrt(sum(range_vector**2))
        range_unit = range_vector / range_mag
        
        ! 计算距离
        measurement(1) = range_mag
        
        ! 计算方位角 (在站心坐标系中)
        azimuth = atan2(range_unit(2), range_unit(1)) * 180.0_DP / (4.0_DP * atan(1.0_DP))
        if (azimuth < 0.0_DP) azimuth = azimuth + 360.0_DP
        measurement(2) = azimuth
        
        ! 计算俯仰角
        elevation = asin(range_unit(3)) * 180.0_DP / (4.0_DP * atan(1.0_DP))
        measurement(3) = elevation
    end subroutine compute_radar_measurement
    
    subroutine compute_optical_measurement(position, station, measurement)
        real(DP), dimension(3), intent(in) :: position
        type(observation_station), intent(in) :: station
        real(DP), dimension(3), intent(out) :: measurement
        
        real(DP), dimension(3) :: range_vector, range_unit
        real(DP) :: right_ascension, declination
        
        ! 计算距离向量
        range_vector = position - station%eci_position
        range_unit = range_vector / sqrt(sum(range_vector**2))
        
        ! 计算赤经
        right_ascension = atan2(range_unit(2), range_unit(1)) * 180.0_DP / (4.0_DP * atan(1.0_DP))
        if (right_ascension < 0.0_DP) right_ascension = right_ascension + 360.0_DP
        measurement(1) = right_ascension
        
        ! 计算赤纬
        declination = asin(range_unit(3)) * 180.0_DP / (4.0_DP * atan(1.0_DP))
        measurement(2) = declination
        
        ! 第三个分量设为0 (光学观测通常只有两个角度)
        measurement(3) = 0.0_DP
    end subroutine compute_optical_measurement
    
    subroutine compute_gps_measurement(position, measurement)
        real(DP), dimension(3), intent(in) :: position
        real(DP), dimension(3), intent(out) :: measurement
        
        ! GPS观测直接给出ECI坐标系中的位置
        measurement = position
    end subroutine compute_gps_measurement
    
    subroutine set_default_station(station)
        type(observation_station), intent(out) :: station
        
        ! 设置默认观测站 (北京天文台)
        station%name = 'Beijing_Observatory'
        station%latitude = 39.9_DP  ! 度
        station%longitude = 116.4_DP  ! 度
        station%altitude = 0.05_DP  ! km
        station%station_type = 'RADAR'
        
        ! 计算ECI位置 (简化处理，假设在J2000历元)
        call geodetic_to_eci(station%latitude, station%longitude, station%altitude, &
                           station%eci_position)
    end subroutine set_default_station
    
    subroutine add_measurement_noise(measurement, noise_std)
        real(DP), dimension(3), intent(inout) :: measurement
        real(DP), dimension(3), intent(in) :: noise_std
        
        real(DP), dimension(3) :: noise
        integer :: i
        
        ! 添加高斯噪声
        do i = 1, 3
            call random_gaussian(0.0_DP, noise_std(i), noise(i))
        end do
        
        measurement = measurement + noise
    end subroutine add_measurement_noise
    
    subroutine random_gaussian(mean, std, value)
        real(DP), intent(in) :: mean, std
        real(DP), intent(out) :: value
        
        real(DP) :: u1, u2, z
        
        ! Box-Muller变换生成高斯随机数
        call random_number(u1)
        call random_number(u2)
        
        z = sqrt(-2.0_DP * log(u1)) * cos(2.0_DP * 4.0_DP * atan(1.0_DP) * u2)
        value = mean + std * z
    end subroutine random_gaussian
    
    subroutine compute_visibility(position, station, elevation_limit, is_visible)
        real(DP), dimension(3), intent(in) :: position
        type(observation_station), intent(in) :: station
        real(DP), intent(in) :: elevation_limit
        logical, intent(out) :: is_visible
        
        real(DP), dimension(3) :: range_vector, range_unit
        real(DP) :: elevation
        
        ! 计算距离向量
        range_vector = position - station%eci_position
        range_unit = range_vector / sqrt(sum(range_vector**2))
        
        ! 计算俯仰角
        elevation = asin(range_unit(3)) * 180.0_DP / (4.0_DP * atan(1.0_DP))
        
        ! 检查可见性
        is_visible = (elevation >= elevation_limit)
    end subroutine compute_visibility
    
    subroutine compute_range_rate(position, velocity, station, range_rate)
        real(DP), dimension(3), intent(in) :: position, velocity
        type(observation_station), intent(in) :: station
        real(DP), intent(out) :: range_rate
        
        real(DP), dimension(3) :: range_vector, range_unit
        
        ! 计算距离向量
        range_vector = position - station%eci_position
        range_unit = range_vector / sqrt(sum(range_vector**2))
        
        ! 计算距离变化率
        range_rate = sum(velocity * range_unit)
    end subroutine compute_range_rate
    
    subroutine compute_doppler_shift(range_rate, frequency, doppler_shift)
        real(DP), intent(in) :: range_rate, frequency
        real(DP), intent(out) :: doppler_shift
        
        real(DP), parameter :: SPEED_OF_LIGHT = 299792458.0_DP  ! m/s
        
        ! 多普勒频移: Δf = -f * v_r / c
        doppler_shift = -frequency * range_rate * 1000.0_DP / SPEED_OF_LIGHT
    end subroutine compute_doppler_shift

end module cat_measurement_model_module
