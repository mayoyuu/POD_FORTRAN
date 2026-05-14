!--------------------------------------------------------------------------------------------------------------
!> # POD Measurement Model DA Module
!> 
!> Differential Algebra (DA) version of the measurement model.
!> Propagates state uncertainties (Taylor polynomials) directly into measurement space
!> (RA/DEC or Range/Az/El) using DACE engine overloading.
!> 
!> 内存安全说明：
!> - 使用临时变量池 (MeasurementTempPool) 管理所有 DA/AlgebraicVector 临时变量
!> - 所有 AlgebraicVector 运算均使用显式 subroutine 版本 (vec_add, vec_sub, vec_div 等)
!> - DA 标量运算也使用显式 subroutine 版本 (da_add, da_sub, da_mul, da_div 等)
!> - 避免 sub(a,b,a) 等别名模式，始终使用临时变量中转
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

    ! =========================================================
    ! 临时变量池：集中管理 DA 运算所需的临时变量，
    ! 避免在循环中反复申请/释放 DA 句柄导致内存泄露。
    ! =========================================================
    type :: MeasurementTempPool
        ! --- AlgebraicVector 临时变量 (均为 3 维) ---
        type(AlgebraicVector) :: pos_j2000   ! 从状态向量提取的前 3 维位置
        type(AlgebraicVector) :: rel_j2000   ! 光学: J2000 相对向量
        type(AlgebraicVector) :: rel_unit    ! 光学: 归一化方向向量
        type(AlgebraicVector) :: pos_itrf    ! 雷达: 地固系位置
        type(AlgebraicVector) :: rel_itrf    ! 雷达: 地固系相对位置
        type(AlgebraicVector) :: rel_enu     ! 雷达: ENU 坐标系相对位置
        
        ! --- DA 标量临时变量 ---
        type(DA) :: range_mag   ! 斜距 (光学/雷达共用)
        type(DA) :: ra          ! 赤经
        type(DA) :: dec         ! 赤纬
        type(DA) :: azimuth     ! 方位角
        type(DA) :: elevation   ! 俯仰角
        type(DA) :: tmp_da1     ! 通用临时 DA (用于象限修正等)
    contains
        procedure :: init => pool_init
        procedure :: destroy => pool_destroy
    end type MeasurementTempPool

contains

    ! =========================================================
    ! 临时变量池初始化
    ! =========================================================
    subroutine pool_init(this)
        class(MeasurementTempPool), intent(inout) :: this
        
        ! AlgebraicVector 临时变量 (均为 3 维)
        call this%pos_j2000%init(3)
        call this%rel_j2000%init(3)
        call this%rel_unit%init(3)
        call this%pos_itrf%init(3)
        call this%rel_itrf%init(3)
        call this%rel_enu%init(3)
        
        ! DA 标量临时变量
        call this%range_mag%init()
        call this%ra%init()
        call this%dec%init()
        call this%azimuth%init()
        call this%elevation%init()
        call this%tmp_da1%init()
    end subroutine pool_init

    ! =========================================================
    ! 临时变量池销毁
    ! =========================================================
    subroutine pool_destroy(this)
        class(MeasurementTempPool), intent(inout) :: this
        
        call this%pos_j2000%destroy()
        call this%rel_j2000%destroy()
        call this%rel_unit%destroy()
        call this%pos_itrf%destroy()
        call this%rel_itrf%destroy()
        call this%rel_enu%destroy()
        
        call this%range_mag%destroy()
        call this%ra%destroy()
        call this%dec%destroy()
        call this%azimuth%destroy()
        call this%elevation%destroy()
        call this%tmp_da1%destroy()
    end subroutine pool_destroy

    !> DA 顶层测量派发器
    subroutine compute_measurement_da(state_da, et, station, measurement_type, measurement_da)
        class(AlgebraicVector), intent(in)    :: state_da         ! J2000系 DA 状态向量 [x,y,z,vx,vy,vz]
        real(DP), intent(in)                  :: et               ! 星历时间
        type(observation_station), intent(in) :: station          ! 测站信息
        character(len=*), intent(in)          :: measurement_type
        type(AlgebraicVector), intent(inout)  :: measurement_da   ! 输出 DA 测量向量
        
        type(MeasurementTempPool) :: pool
        
        ! 初始化临时变量池
        call pool%init()
        
        ! 1. 从 6 维状态中提取前 3 维位置 DA 向量
        call pool%pos_j2000%set(1, state_da%elements(1))
        call pool%pos_j2000%set(2, state_da%elements(2))
        call pool%pos_j2000%set(3, state_da%elements(3))
        
        select case (trim(measurement_type))
            case ('OPTICAL')
                call measurement_da%init(2) ! DA RA, DA DEC
                call compute_optical_measurement_da(pool, et, station, measurement_da)
                
            case ('RADAR')
                call measurement_da%init(3) ! DA Range, DA Azimuth, DA Elevation
                call compute_radar_measurement_da(pool, et, station, measurement_da)
                
            case default
                write(*, *) '[ERROR] 未知的测量类型: ', trim(measurement_type)
        end select
        
        ! 销毁临时变量池
        call pool%destroy()
    end subroutine compute_measurement_da


    !> =====================================================================
    !> DA 光学测量方程 (GCRS/J2000 -> DA 赤经赤纬)
    !> =====================================================================
    subroutine compute_optical_measurement_da(pool, et, station, measurement_da)
        type(MeasurementTempPool), intent(inout) :: pool
        real(DP), intent(in)                  :: et
        type(observation_station), intent(in) :: station
        type(AlgebraicVector), intent(inout)  :: measurement_da
        
        real(DP), dimension(3,3) :: rot_itrf_to_j2000
        real(DP), dimension(3)   :: obs_j2000
        
        ! 1. 获取转换矩阵并计算测站在 J2000 下的双精度常数位置
        call get_frame_transform('ITRF93', 'J2000', et, rot_itrf_to_j2000)
        obs_j2000 = matmul(rot_itrf_to_j2000, station%ecef_position)
        ! call vec_matmul(rot_itrf_to_j2000, station%ecef_position, obs_j2000)
        
        ! 2. 相对向量: rel_j2000 = pos_j2000 - obs_j2000
        call vec_sub(pool%pos_j2000, obs_j2000, pool%rel_j2000)
        
        ! 3. 归一化: range_mag = ||rel_j2000||, rel_unit = rel_j2000 / range_mag
        call vector_norm2_sub(pool%rel_j2000, pool%range_mag)
        call vec_div(pool%rel_j2000, pool%range_mag, pool%rel_unit)
        
        ! 4. 计算 DA 赤经赤纬
        ! ra = atan2(rel_unit(2), rel_unit(1))
        call da_atan2_sub(pool%rel_unit%elements(2), pool%rel_unit%elements(1), pool%ra)
        
        ! 使用常数项判断象限，并加上 2PI
        if (pool%ra%cons() < 0.0_DP) then
            ! ra = ra + 2.0_DP * PI  (使用临时变量避免别名)
            call da_add(pool%ra, 2.0_DP * PI, pool%tmp_da1)
            pool%ra = pool%tmp_da1
        end if
        
        ! dec = asin(rel_unit(3))
        call da_asin_sub(pool%rel_unit%elements(3), pool%dec)
        
        ! 5. 存入结果
        call measurement_da%set(1, pool%ra)
        call measurement_da%set(2, pool%dec)
        
    end subroutine compute_optical_measurement_da


    !> =====================================================================
    !> DA 雷达测量方程 (J2000 -> DA 极坐标: 距离, 方位角, 俯仰角)
    !> =====================================================================
    subroutine compute_radar_measurement_da(pool, et, station, measurement_da)
        type(MeasurementTempPool), intent(inout) :: pool
        real(DP), intent(in)                  :: et
        type(observation_station), intent(in) :: station
        type(AlgebraicVector), intent(inout)  :: measurement_da
        
        real(DP), dimension(3,3) :: rot_j2000_to_itrf, rot_itrf_to_enu
        real(DP)                 :: lat_rad, lon_rad, sin_lat, cos_lat, sin_lon, cos_lon
        
        ! 1. 转换目标位置到地固系: pos_itrf = matmul(rot_j2000_to_itrf, pos_j2000)
        call get_frame_transform('J2000', 'ITRF93', et, rot_j2000_to_itrf)
        call vec_matmul(rot_j2000_to_itrf, pool%pos_j2000, pool%pos_itrf)
        
        ! 2. 相对位置: rel_itrf = pos_itrf - station%ecef_position
        call vec_sub(pool%pos_itrf, station%ecef_position, pool%rel_itrf)
        
        ! 3. DA 斜距: range_mag = ||rel_itrf||
        call vector_norm2_sub(pool%rel_itrf, pool%range_mag)
        
        ! 4. 构建 ITRF 到 ENU 的实数旋转矩阵
        lat_rad = station%latitude * PI / 180.0_DP
        lon_rad = station%longitude * PI / 180.0_DP
        sin_lat = sin(lat_rad); cos_lat = cos(lat_rad)
        sin_lon = sin(lon_rad); cos_lon = cos(lon_rad)
        
        rot_itrf_to_enu(1,:) = [-sin_lon, cos_lon, 0.0_DP]
        rot_itrf_to_enu(2,:) = [-sin_lat*cos_lon, -sin_lat*sin_lon, cos_lat]
        rot_itrf_to_enu(3,:) = [cos_lat*cos_lon, cos_lat*sin_lon, sin_lat]
        
        ! 5. 旋转 DA 向量到 ENU 坐标系: rel_enu = matmul(rot_itrf_to_enu, rel_itrf)
        call vec_matmul(rot_itrf_to_enu, pool%rel_itrf, pool%rel_enu)
        
        ! 6. 计算 DA 方位角和俯仰角
        ! azimuth = atan2(rel_enu(1), rel_enu(2))
        call da_atan2_sub(pool%rel_enu%elements(1), pool%rel_enu%elements(2), pool%azimuth)
        if (pool%azimuth%cons() < 0.0_DP) then
            ! azimuth = azimuth + 2.0_DP * PI (使用临时变量避免别名)
            call da_add(pool%azimuth, 2.0_DP * PI, pool%tmp_da1)
            pool%azimuth = pool%tmp_da1
        end if
        
        ! elevation = asin(rel_enu(3) / range_mag)
        ! 先计算 rel_enu(3) / range_mag 存入 tmp_da1，再 asin
        call da_div(pool%rel_enu%elements(3), pool%range_mag, pool%tmp_da1)
        call da_asin_sub(pool%tmp_da1, pool%elevation)
        
        ! 7. 存入结果
        call measurement_da%set(1, pool%range_mag)
        call measurement_da%set(2, pool%azimuth)
        call measurement_da%set(3, pool%elevation)
        
    end subroutine compute_radar_measurement_da

end module pod_measurement_da_module
