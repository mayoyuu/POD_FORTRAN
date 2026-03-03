program time_test
    use pod_global, only: DP, pod_init, pod_cleanup
    use pod_spice, only: str2et, et2utc, spice_init, spice_cleanup
    
    implicit none
    
    real(DP) :: time_et
    character(len=100) :: utc_string, utc_output
    
    write(*, *) '=== SPICE时间转换测试 ==='
    
    ! 初始化系统
    call pod_init()
    call spice_init()
    
    ! 测试时间转换
    write(*, *) ''
    write(*, *) '1. 测试时间转换功能:'
    
    ! 测试2024年1月1日
    utc_string = '2024-01-01T12:00:00'
    write(*, *) '输入UTC时间: ', trim(utc_string)
    
    call str2et(utc_string, time_et)
    write(*, *) '转换到ET时间: ', time_et
    
    call et2utc(time_et, 'ISOC', 3, utc_output)
    write(*, *) '转换回UTC时间: ', trim(utc_output)
    
    ! 测试2024年6月1日
    utc_string = '2024-06-01T12:00:00'
    write(*, *) ''
    write(*, *) '输入UTC时间: ', trim(utc_string)
    
    call str2et(utc_string, time_et)
    write(*, *) '转换到ET时间: ', time_et
    
    call et2utc(time_et, 'ISOC', 3, utc_output)
    write(*, *) '转换回UTC时间: ', trim(utc_output)
    
    ! 测试2024年12月1日
    utc_string = '2024-12-01T12:00:00'
    write(*, *) ''
    write(*, *) '输入UTC时间: ', trim(utc_string)
    
    call str2et(utc_string, time_et)
    write(*, *) '转换到ET时间: ', time_et
    
    call et2utc(time_et, 'ISOC', 3, utc_output)
    write(*, *) '转换回UTC时间: ', trim(utc_output)
    
    ! 测试不同格式
    write(*, *) ''
    write(*, *) '2. 测试不同时间格式:'
    
    ! 测试ISO格式
    utc_string = '2024-01-01T00:00:00'
    call str2et(utc_string, time_et)
    write(*, *) 'ISO格式: ', trim(utc_string), ' -> ET: ', time_et
    
    ! 测试标准格式
    utc_string = '2024 JAN 01 12:00:00'
    call str2et(utc_string, time_et)
    write(*, *) '标准格式: ', trim(utc_string), ' -> ET: ', time_et
    
    write(*, *) ''
    write(*, *) '=== 时间转换测试完成 ==='
    
    ! 清理系统
    call spice_cleanup()
    call pod_cleanup()
    
end program time_test
