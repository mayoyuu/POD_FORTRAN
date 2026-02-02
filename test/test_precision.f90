program test_precision
    implicit none
    
    ! 定义精度类型
    integer, parameter :: DP = selected_real_kind(15, 307)
    integer, parameter :: SP = selected_real_kind(6, 37)
    
    ! 声明变量
    real(DP) :: dp_var
    real(SP) :: sp_var
    
    write(*, *) '=========================================='
    write(*, *) '    精度类型测试程序'
    write(*, *) '=========================================='
    write(*, *)
    
    ! 测试双精度
    write(*, *) '=== 双精度 (DP) 测试 ==='
    write(*, '(A, I0)') 'DP kind值: ', DP
    write(*, '(A, I0)') 'DP 存储大小 (字节): ', storage_size(dp_var) / 8
    write(*, '(A, I0)') 'DP 精度位数: ', precision(dp_var)
    write(*, '(A, I0)') 'DP 指数范围: ', range(dp_var)
    write(*, *)
    
    ! 测试单精度
    write(*, *) '=== 单精度 (SP) 测试 ==='
    write(*, '(A, I0)') 'SP kind值: ', SP
    write(*, '(A, I0)') 'SP 存储大小 (字节): ', storage_size(sp_var) / 8
    write(*, '(A, I0)') 'SP 精度位数: ', precision(sp_var)
    write(*, '(A, I0)') 'SP 指数范围: ', range(sp_var)
    write(*, *)
    
    ! 测试数值范围
    write(*, *) '=== 数值范围测试 ==='
    dp_var = 1.0_DP
    write(*, '(A, ES25.15)') 'DP 1.0: ', dp_var
    
    dp_var = 1.0e-300_DP
    write(*, '(A, ES25.15)') 'DP 1.0e-300: ', dp_var
    
    dp_var = 1.0e+300_DP
    write(*, '(A, ES25.15)') 'DP 1.0e+300: ', dp_var
    
    write(*, *)
    
    ! 测试精度
    write(*, *) '=== 精度测试 ==='
    dp_var = 3.141592653589793238462643383279_DP
    write(*, '(A, F30.20)') 'DP π (30位): ', dp_var
    
    sp_var = 3.141592653589793238462643383279_SP
    write(*, '(A, F30.20)') 'SP π (30位): ', sp_var
    
    write(*, *)
    
    ! 测试物理常数
    write(*, *) '=== 物理常数测试 ==='
    dp_var = 299792458.0_DP  ! 光速
    write(*, '(A, F20.6)') '光速 (m/s): ', dp_var
    
    dp_var = 6.67430e-11_DP  ! 万有引力常数
    write(*, '(A, ES20.10)') '万有引力常数: ', dp_var
    
    dp_var = 398600.5_DP     ! 地球引力参数
    write(*, '(A, F20.6)') '地球引力参数: ', dp_var
    
    write(*, *)
    write(*, *) '=========================================='
    write(*, *) '    测试完成'
    write(*, *) '=========================================='
    
end program test_precision
