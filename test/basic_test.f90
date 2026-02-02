program basic_test
    implicit none
    
    real(8) :: test_value
    real(8), dimension(3) :: test_vector
    
    write(*, *) '=== 基本Fortran功能测试 ==='
    
    ! 测试基本计算
    test_value = 3.14159_8
    test_vector = [1.0_8, 2.0_8, 3.0_8]
    
    write(*, *) '测试值: ', test_value
    write(*, *) '测试向量: ', test_vector
    write(*, *) '向量模长: ', sqrt(sum(test_vector**2))
    
    write(*, *) ''
    write(*, *) '=== 测试完成 ==='
    
end program basic_test
