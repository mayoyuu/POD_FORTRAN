!> @file pod_measurement_base_module.f90
!> @brief 测量系统的底层公共数据结构
module pod_measurement_base_module
    use pod_global, only: DP, MAX_STRING_LEN
    implicit none
    
    private
    public :: observation_station,PI
    
    !> 测站结构体 (唯一来源)
    type observation_station
        character(len=MAX_STRING_LEN) :: name
        real(DP) :: latitude, longitude, altitude
        real(DP), dimension(3) :: ecef_position  ! 保存地固系(ITRS/ECEF)坐标
        character(len=MAX_STRING_LEN) :: station_type  ! 'RADAR', 'OPTICAL', 'GPS'
    end type observation_station

    real(DP), parameter :: PI = 3.14159265358979323846_DP
end module pod_measurement_base_module