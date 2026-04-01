program calc_normalization_units
    implicit none
    
    integer, parameter :: DP = selected_real_kind(15, 307)
    
    real(DP) :: distance_unit, time_unit, mass_unit
    real(DP) :: velocity_unit, accel_norm
    real(DP) :: G, mu_earth, mu_moon, mu_sys
    
    ! 物理常数输入
    distance_unit = 384400.0_DP
    mu_earth      = 398600.43550702254_DP
    mu_moon       = 4902.800118_DP
    mu_sys        = mu_earth + mu_moon
    G             = 6.67430e-20_DP
    
    ! 派生单位计算
    time_unit     = sqrt((distance_unit**3) / mu_sys)
    mass_unit     = mu_sys / G
    velocity_unit = distance_unit / time_unit
    accel_norm    = distance_unit / (time_unit**2)
    
    ! 输出高精度结果
    write(*,*) '=== 地月空间归一化常数 ==='
    write(*, '(A, F25.15, A)') ' LU   = ', distance_unit, ' _DP'
    write(*, '(A, F25.15, A)') ' TU   = ', time_unit, ' _DP'
    write(*, '(A, ES25.15, A)') ' MU   = ', mass_unit, ' _DP'
    write(*, '(A, F25.15, A)') ' VU   = ', velocity_unit, ' _DP'
    write(*, '(A, ES25.15, A)') ' AccU = ', accel_norm, ' _DP'
    
end program calc_normalization_units