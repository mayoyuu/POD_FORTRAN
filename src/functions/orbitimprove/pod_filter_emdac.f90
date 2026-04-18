module pod_measurement_model_module
    use pod_global, only: DP
    use pod_basicmath_module, only: norm_vector
    implicit none
    private
    public :: compute_optical_measurement

contains

    !> 计算 GCRS 状态到测站视线的 赤经(RA) 和 赤纬(DEC)
    !> 对应 C++ 中的 GCRS_2_RADEC
    subroutine compute_optical_measurement(state_gcrs, obs_ecef, HG_t, z_out)
        real(DP), dimension(6), intent(in)  :: state_gcrs  ! [x, y, z, vx, vy, vz] in GCRS
        real(DP), dimension(3), intent(in)  :: obs_ecef    ! 测站在 ECEF 系下的坐标
        real(DP), dimension(3,3), intent(in):: HG_t        ! ECEF 到 GCRS 的转换矩阵 (考虑了岁差、章动、EOP)
        real(DP), dimension(2), intent(out) :: z_out       ! [RA, DEC] 弧度
        
        real(DP), dimension(3) :: pos_gcrs, obs_gcrs, rel_gcrs
        real(DP) :: royl, ra, dec
        real(DP) :: PI = 3.14159265358979323846_DP
        
        pos_gcrs = state_gcrs(1:3)
        
        ! 1. 将测站坐标转到 GCRS 系 (Rxp_out 操作)
        obs_gcrs = matmul(HG_t, obs_ecef)
        
        ! 2. 计算视线向量 (Line of Sight)
        rel_gcrs = pos_gcrs - obs_gcrs
        
        ! 3. 求解赤经 RA (0 到 2PI)
        ra = atan2(rel_gcrs(2), rel_gcrs(1))
        if (ra < 0.0_DP) ra = ra + 2.0_DP * PI
        
        ! 4. 求解赤纬 DEC
        royl = norm_vector(rel_gcrs)
        dec = asin(rel_gcrs(3) / royl)
        
        z_out(1) = ra
        z_out(2) = dec
    end subroutine compute_optical_measurement

end module pod_measurement_model_module