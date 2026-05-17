!> @file pod_error_analysis_module.f90
!> @brief 轨道误差分析：计算估计状态相对于参考轨道的各项误差指标
module pod_error_analysis_module
    use pod_global, only: DP
    use pod_basicmath_module, only: inverse_and_determinant

    implicit none
    private

    public :: compute_orbit_error

contains

    !> ======================================================================
    !> 计算估计状态相对于参考轨道（真值）的误差
    !>
    !> 计算时机：时间更新之后、测量更新之前
    !> ======================================================================
    subroutine compute_orbit_error(est_state, est_cov, ref_state, &
                                   pos_err, vel_err, pos_rms, vel_rms, mahalanobis_d)
        real(DP), intent(in)  :: est_state(6)   ! 滤波器估计 [X,Y,Z,Vx,Vy,Vz]
        real(DP), intent(in)  :: est_cov(6,6)   ! 滤波器协方差 P_{k|k-1}
        real(DP), intent(in)  :: ref_state(6)   ! 参考轨道（真值）
        real(DP), intent(out) :: pos_err(3)     ! dX, dY, dZ (km)
        real(DP), intent(out) :: vel_err(3)     ! dVx, dVy, dVz (km/s)
        real(DP), intent(out) :: pos_rms        ! sqrt(dX^2 + dY^2 + dZ^2) (km)
        real(DP), intent(out) :: vel_rms        ! sqrt(dVx^2 + dVy^2 + dVz^2) (km/s)
        real(DP), intent(out) :: mahalanobis_d  ! sqrt(dx^T * P^{-1} * dx)

        real(DP) :: dx(6), inv_cov(6,6), det_cov
        integer  :: info

        ! 1. 位置与速度误差
        pos_err = est_state(1:3) - ref_state(1:3)
        vel_err = est_state(4:6) - ref_state(4:6)

        ! 2. 位置/速度 RMS
        pos_rms = sqrt(sum(pos_err**2))
        vel_rms = sqrt(sum(vel_err**2))

        ! 3. 马氏距离 sqrt(dx^T * P^{-1} * dx)
        dx = est_state - ref_state
        call inverse_and_determinant(est_cov, inv_cov, det_cov, info)

        if (info == 0) then
            mahalanobis_d = sqrt(dot_product(dx, matmul(inv_cov, dx)))
        else
            mahalanobis_d = -1.0_DP  ! 标记：协方差奇异，无法计算马氏距离
        end if

    end subroutine compute_orbit_error

end module pod_error_analysis_module
