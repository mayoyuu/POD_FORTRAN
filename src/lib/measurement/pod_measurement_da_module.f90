!--------------------------------------------------------------------------------------------------------------
!> # POD Measurement Model DA Module
!> 
!> Differential Algebra (DA) version of the measurement model.
!> Propagates state uncertainties (Taylor polynomials) directly into measurement space
!> (RA/DEC or Range/Az/El) using DACE engine overloading.
!--------------------------------------------------------------------------------------------------------------

module pod_measurement_da_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_spice, only: get_frame_transform
    use pod_basicmath_module, only: PI
    use pod_measurement_base_module, only: observation_station
    use pod_dace_classes
    
    implicit none
    private
    
    public :: compute_measurement_da


contains

    !> DA 顶层测量派发器
    subroutine compute_measurement_da(state_da, et, station, measurement_type, measurement_da)
        class(AlgebraicVector), intent(in)    :: state_da         ! J2000系 DA 状态向量 [x,y,z,vx,vy,vz]
        real(DP), intent(in)                  :: et               ! 星历时间
        type(observation_station), intent(in) :: station          ! 测站信息
        character(len=*), intent(in)          :: measurement_type
        type(AlgebraicVector), intent(inout)  :: measurement_da   ! 输出 DA 测量向量
        
        type(AlgebraicVector) :: pos_j2000
        
        ! 1. 从 6 维状态中提取前 3 维位置 DA 向量
        call pos_j2000%init(3)
        call pos_j2000%set(1, state_da%elements(1))
        call pos_j2000%set(2, state_da%elements(2))
        call pos_j2000%set(3, state_da%elements(3))
        
        select case (trim(measurement_type))
            case ('OPTICAL')
                call measurement_da%init(2) ! DA RA, DA DEC
                call compute_optical_measurement_da(pos_j2000, et, station, measurement_da)
                
            case ('RADAR')
                call measurement_da%init(3) ! DA Range, DA Azimuth, DA Elevation
                call compute_radar_measurement_da(pos_j2000, et, station, measurement_da)
                
            case default
                write(*, *) '[ERROR] 未知的测量类型: ', trim(measurement_type)
        end select
        
        ! 释放局部变量
        call pos_j2000%destroy()
    end subroutine compute_measurement_da


    !> =====================================================================
    !> DA 光学测量方程 (GCRS/J2000 -> DA 赤经赤纬)
    !> =====================================================================
    subroutine compute_optical_measurement_da(pos_j2000, et, station, measurement_da)
        class(AlgebraicVector), intent(in)    :: pos_j2000
        real(DP), intent(in)                  :: et
        type(observation_station), intent(in) :: station
        type(AlgebraicVector), intent(inout)  :: measurement_da
        
        real(DP), dimension(3,3) :: rot_itrf_to_j2000
        real(DP), dimension(3)   :: obs_j2000
        type(AlgebraicVector)    :: rel_j2000, rel_unit
        type(DA)                 :: ra, dec, range_mag
        
        ! 1. 获取转换矩阵并计算测站在 J2000 下的双精度常数位置
        call get_frame_transform('ITRF93', 'J2000', et, rot_itrf_to_j2000)
        obs_j2000 = matmul(rot_itrf_to_j2000, station%ecef_position)
        
        ! 2. 相对向量 (触发 DA 向量 - 实数数组 重载)
        rel_j2000 = pos_j2000 - obs_j2000
        
        ! 3. 归一化 (触发 Vector norm2 和 向量 / DA 重载)
        range_mag = rel_j2000%norm2()
        rel_unit = rel_j2000 / range_mag
        
        ! 4. 计算 DA 赤经赤纬
        ra = atan2(rel_unit%elements(2), rel_unit%elements(1))
        
        ! 使用常数项判断象限，并加上 2PI (触发 DA + 实数 重载)
        if (ra%cons() < 0.0_DP) then
            ra = ra + 2.0_DP * PI
        end if
        
        dec = asin(rel_unit%elements(3))
        
        ! 5. 存入结果
        call measurement_da%set(1, ra)
        call measurement_da%set(2, dec)
        
    end subroutine compute_optical_measurement_da


    !> =====================================================================
    !> DA 雷达测量方程 (J2000 -> DA 极坐标: 距离, 方位角, 俯仰角)
    !> =====================================================================
    subroutine compute_radar_measurement_da(pos_j2000, et, station, measurement_da)
        class(AlgebraicVector), intent(in)    :: pos_j2000
        real(DP), intent(in)                  :: et
        type(observation_station), intent(in) :: station
        type(AlgebraicVector), intent(inout)  :: measurement_da
        
        real(DP), dimension(3,3) :: rot_j2000_to_itrf, rot_itrf_to_enu
        type(AlgebraicVector)    :: pos_itrf, rel_itrf, rel_enu
        type(DA)                 :: range_mag, azimuth, elevation
        real(DP)                 :: lat_rad, lon_rad, sin_lat, cos_lat, sin_lon, cos_lon
        
        ! 1. 转换目标位置到地固系 (触发 实数矩阵 * DA向量 重载)
        call get_frame_transform('J2000', 'ITRF93', et, rot_j2000_to_itrf)
        pos_itrf = matmul(rot_j2000_to_itrf, pos_j2000)
        
        ! 2. 相对位置 (触发 DA向量 - 实数数组 重载)
        rel_itrf = pos_itrf - station%ecef_position
        
        ! 3. DA 斜距
        range_mag = rel_itrf%norm2()
        
        ! 4. 构建 ITRF 到 ENU 的实数旋转矩阵
        lat_rad = station%latitude * PI / 180.0_DP
        lon_rad = station%longitude * PI / 180.0_DP
        sin_lat = sin(lat_rad); cos_lat = cos(lat_rad)
        sin_lon = sin(lon_rad); cos_lon = cos(lon_rad)
        
        rot_itrf_to_enu(1,:) = [-sin_lon, cos_lon, 0.0_DP]
        rot_itrf_to_enu(2,:) = [-sin_lat*cos_lon, -sin_lat*sin_lon, cos_lat]
        rot_itrf_to_enu(3,:) = [cos_lat*cos_lon, cos_lat*sin_lon, sin_lat]
        
        ! 5. 旋转 DA 向量到 ENU 坐标系 (触发 实数矩阵 * DA向量 重载)
        rel_enu = matmul(rot_itrf_to_enu, rel_itrf)
        
        ! 6. 计算 DA 方位角和俯仰角
        azimuth = atan2(rel_enu%elements(1), rel_enu%elements(2))
        if (azimuth%cons() < 0.0_DP) then
            azimuth = azimuth + 2.0_DP * PI
        end if
        
        ! 注意这里使用了 rel_enu(3) / range_mag (DA / DA)
        elevation = asin(rel_enu%elements(3) / range_mag)
        
        ! 7. 存入结果
        call measurement_da%set(1, range_mag)
        call measurement_da%set(2, azimuth)
        call measurement_da%set(3, elevation)
        
    end subroutine compute_radar_measurement_da

end module pod_measurement_da_module