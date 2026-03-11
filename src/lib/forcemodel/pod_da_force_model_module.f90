module pod_da_force_model_module
    use pod_global, only: DP
    use pod_dace_classes
    implicit none

    ! 地球引力常数 (km^3/s^2)
    real(DP), parameter :: MU_EARTH = 398600.4415_DP

contains

    ! =========================================================
    ! 二体问题加速度计算 (DA版)
    ! =========================================================
    subroutine da_compute_acceleration(pos, vel, time, acc)
        type(AlgebraicVector), intent(in) :: pos, vel
        real(DP), intent(in) :: time
        type(AlgebraicVector), intent(inout) :: acc
        
        type(DA) :: r_sq, r_mag, r_cube, mu_over_r3
        integer :: i
        
        ! 1. 初始化局部 DA 标量
        call r_sq%init(); call r_mag%init()
        call r_cube%init(); call mu_over_r3%init()
        if (acc%size /= 3) call acc%init(3)
        
        ! 2. 计算距离的平方: r^2 = x^2 + y^2 + z^2
        ! 这里完全依赖于我们重载的 DA 乘法和加法
        r_sq = (pos%elements(1) * pos%elements(1)) + &
               (pos%elements(2) * pos%elements(2)) + &
               (pos%elements(3) * pos%elements(3))
        
        ! 3. 计算 r 和 r^3
        r_mag = sqrt(r_sq)
        r_cube = r_mag * r_mag * r_mag
        
        ! 4. 计算引力系数: -mu / r^3
        mu_over_r3 = -MU_EARTH / r_cube
        
        ! 5. 计算加速度向量: a = (-mu/r^3) * r
        do i = 1, 3
            acc%elements(i) = mu_over_r3 * pos%elements(i)
        end do
        
        ! 6. 安全释放内存
        call r_sq%destroy(); call r_mag%destroy()
        call r_cube%destroy(); call mu_over_r3%destroy()
        
    end subroutine da_compute_acceleration

end module pod_da_force_model_module