!--------------------------------------------------------------------------------------------------------------
!> # POD Measurement Model Module
!> 
!> Provides observation equations for transforming spacecraft state vectors 
!> into predicted station measurements (Optical RA/DEC and Radar Range/Az/El).
!--------------------------------------------------------------------------------------------------------------

module pod_measurement_model_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_spice, only: get_frame_transform
    use pod_basicmath_module, only: vector_magnitude, normalize_vector
    
    implicit none
    private
    
    ! 对外暴露的接口
    public :: observation_station
    public :: set_station_from_geodetic
    public :: compute_measurement
    
    !> 观测站类型重构
    type observation_station
        character(len=MAX_STRING_LEN) :: name
        real(DP) :: latitude, longitude, altitude
        real(DP), dimension(3) :: ecef_position  ! 保存地固系(ITRS/ECEF)坐标
        character(len=MAX_STRING_LEN) :: station_type  ! 'RADAR', 'OPTICAL', 'GPS'
    end type observation_station
    
    real(DP), parameter :: PI = 3.14159265358979323846_DP

contains

    !> 初始化测站：将经纬高转化为 ECEF
    subroutine set_station_from_geodetic(station, name, lat_deg, lon_deg, alt_km, type_str)
        type(observation_station), intent(out) :: station
        character(len=*), intent(in) :: name
        real(DP), intent(in) :: lat_deg, lon_deg, alt_km
        character(len=*), intent(in) :: type_str
        
        real(DP) :: lat_rad, lon_rad
        real(DP) :: ae, e2, N
        
        station%name = name
        station%latitude = lat_deg
        station%longitude = lon_deg
        station%altitude = alt_km
        station%station_type = type_str
        
        ! 转换为弧度
        lat_rad = lat_deg * PI / 180.0_DP
        lon_rad = lon_deg * PI / 180.0_DP
        
        ! WGS84 椭球参数
        ae = 6378.1370_DP
        e2 = 0.00669438_DP
        
        ! 计算曲率半径 N
        N = ae / sqrt(1.0_DP - e2 * sin(lat_rad)**2)
        
        ! 计算 ECEF 坐标 (X, Y, Z)
        station%ecef_position(1) = (N + alt_km) * cos(lat_rad) * cos(lon_rad)
        station%ecef_position(2) = (N + alt_km) * cos(lat_rad) * sin(lon_rad)
        station%ecef_position(3) = (N * (1.0_DP - e2) + alt_km) * sin(lat_rad)
        
    end subroutine set_station_from_geodetic

    !> 顶层测量派发器
    subroutine compute_measurement(state, et, station, measurement_type, measurement)
        real(DP), dimension(6), intent(in) :: state         ! J2000系状态 [x,y,z,vx,vy,vz]
        real(DP), intent(in)               :: et            ! SPICE 历元时间 (Ephemeris Time)
        type(observation_station), intent(in):: station     ! 测站信息
        character(len=*), intent(in)       :: measurement_type
        real(DP), dimension(:), allocatable, intent(out) :: measurement
        
        select case (trim(measurement_type))
            case ('OPTICAL')
                allocate(measurement(2)) ! RA, DEC
                call compute_optical_measurement(state(1:3), et, station, measurement)
                
            case ('RADAR')
                allocate(measurement(3)) ! Range, Azimuth, Elevation
                call compute_radar_measurement(state(1:3), et, station, measurement)
                
            case default
                write(*, *) '[ERROR] 未知的测量类型: ', trim(measurement_type)
        end select
    end subroutine compute_measurement

    !> =====================================================================
    !> 光学测量方程 (GCRS/J2000 -> 测站赤经赤纬)
    !> =====================================================================
    subroutine compute_optical_measurement(pos_j2000, et, station, measurement)
        real(DP), dimension(3), intent(in) :: pos_j2000     ! 目标在 J2000 下的位置
        real(DP), intent(in)               :: et            ! 星历时间
        type(observation_station), intent(in):: station     ! 测站
        real(DP), dimension(2), intent(out):: measurement   ! 输出 [RA, DEC] (弧度)
        
        real(DP), dimension(3,3) :: rot_itrf_to_j2000
        real(DP), dimension(3)   :: obs_j2000, rel_j2000, rel_unit
        real(DP) :: ra, dec
        
        ! 1. 获取动态坐标转换矩阵：ITRF93 (地固系) -> J2000 (惯性系)
        call get_frame_transform('ITRF93', 'J2000', et, rot_itrf_to_j2000)
        
        ! 2. 将测站坐标从地固系旋转到当前的 J2000 惯性系
        obs_j2000 = matmul(rot_itrf_to_j2000, station%ecef_position)
        
        ! 3. 计算相对视线向量 (Line of Sight) 并归一化
        rel_j2000 = pos_j2000 - obs_j2000
        rel_unit = normalize_vector(rel_j2000)
        
        ! 4. 计算赤经 (Right Ascension) [0, 2π]
        ra = atan2(rel_unit(2), rel_unit(1))
        if (ra < 0.0_DP) ra = ra + 2.0_DP * PI
        
        ! 5. 计算赤纬 (Declination) [-π/2, π/2]
        dec = asin(rel_unit(3))
        
        measurement(1) = ra
        measurement(2) = dec
    end subroutine compute_optical_measurement

    !> =====================================================================
    !> 雷达测量方程 (J2000 -> 测站极坐标: 距离, 方位角, 俯仰角)
    !> =====================================================================
    subroutine compute_radar_measurement(pos_j2000, et, station, measurement)
        real(DP), dimension(3), intent(in) :: pos_j2000
        real(DP), intent(in)               :: et
        type(observation_station), intent(in):: station
        real(DP), dimension(3), intent(out):: measurement
        
        real(DP), dimension(3,3) :: rot_j2000_to_itrf
        real(DP), dimension(3)   :: pos_itrf, rel_itrf, rel_enu
        real(DP) :: range_mag, azimuth, elevation
        real(DP) :: lat_rad, lon_rad
        real(DP) :: sin_lat, cos_lat, sin_lon, cos_lon
        
        ! 1. 获取 J2000 到 ITRF 的转换矩阵
        call get_frame_transform('J2000', 'ITRF93', et, rot_j2000_to_itrf)
        
        ! 2. 将目标位置转换到地固系 (ITRF)
        pos_itrf = matmul(rot_j2000_to_itrf, pos_j2000)
        
        ! 3. 计算地固系下的相对位置向量
        rel_itrf = pos_itrf - station%ecef_position
        
        ! 4. 计算斜距 (Range)
        range_mag = vector_magnitude(rel_itrf)
        
        ! 5. ITRF 到 ENU (East, North, Up) 坐标系转换
        lat_rad = station%latitude * PI / 180.0_DP
        lon_rad = station%longitude * PI / 180.0_DP
        
        sin_lat = sin(lat_rad)
        cos_lat = cos(lat_rad)
        sin_lon = sin(lon_rad)
        cos_lon = cos(lon_rad)
        
        ! 计算 ENU 坐标分量
        rel_enu(1) = -sin_lon * rel_itrf(1) + cos_lon * rel_itrf(2)
        rel_enu(2) = -sin_lat * cos_lon * rel_itrf(1) - sin_lat * sin_lon * rel_itrf(2) + cos_lat * rel_itrf(3)
        rel_enu(3) =  cos_lat * cos_lon * rel_itrf(1) + cos_lat * sin_lon * rel_itrf(2) + sin_lat * rel_itrf(3)
        
        ! 6. 计算方位角 (Azimuth) 和俯仰角 (Elevation)
        ! 方位角通常从正北(North)起算，顺时针为正，所以用 atan2(East, North)
        azimuth = atan2(rel_enu(1), rel_enu(2))
        if (azimuth < 0.0_DP) azimuth = azimuth + 2.0_DP * PI
        
        elevation = asin(rel_enu(3) / range_mag)
        
        ! 返回 [Range(km), Azimuth(rad), Elevation(rad)]
        measurement(1) = range_mag
        measurement(2) = azimuth
        measurement(3) = elevation
    end subroutine compute_radar_measurement

end module pod_measurement_model_module