program test_time_system
    use pod_global, only: DP, pod_init, pod_cleanup
    use pod_time_module, only: time_transfer, utc_to_et, et_to_utc, utc_to_jd, jd_to_utc, &
                              utc_to_mjd, mjd_to_utc, validate_time_string, get_current_utc
    use pod_spice, only: spice_init, spice_cleanup, str2et, et2utc
    
    implicit none
    
    real(DP) :: et_time, jd_time, mjd_time, et_output
    character(len=100) :: utc_string, utc_output, current_time
    integer :: i, status
    logical :: is_valid
    
    write(*, *) '=== CAT 时间系统模块测试 ==='
    
    ! 初始化系统
    call pod_init()
    call spice_init()
    
    ! 测试1: UTC到ET转换
    write(*, *) ''
    write(*, *) '1. 测试UTC到ET转换:'
    
    utc_string = '2024-01-01T12:00:00'
    et_time = utc_to_et(utc_string)
    write(*, *) '  输入UTC: ', trim(utc_string)
    write(*, *) '  输出ET: ', et_time
    
    ! 验证转换
    call et2utc(et_time, 'ISOC', 3, utc_output)
    write(*, *) '  验证UTC: ', trim(utc_output)
    
    ! 测试2: ET到UTC转换
    write(*, *) ''
    write(*, *) '2. 测试ET到UTC转换:'
    
    utc_output = et_to_utc(et_time, 'ISOC', 3)
    write(*, *) '  ET输入: ', et_time
    write(*, *) '  UTC输出: ', trim(utc_output)
    
    ! 测试不同格式
    write(*, *) '  日历格式: ', trim(et_to_utc(et_time, 'C', 3))
    write(*, *) '  儒略日格式: ', trim(et_to_utc(et_time, 'J', 3))
    
    ! 测试3: UTC到儒略日转换
    write(*, *) ''
    write(*, *) '3. 测试UTC到儒略日转换:'
    
    jd_time = utc_to_jd(utc_string)
    write(*, *) '  UTC输入: ', trim(utc_string)
    write(*, *) '  JD输出: ', jd_time
    
    ! 验证转换
    utc_output = jd_to_utc(jd_time, 'ISOC')
    write(*, *) '  验证UTC: ', trim(utc_output)
    
    ! 测试4: 修正儒略日转换
    write(*, *) ''
    write(*, *) '4. 测试修正儒略日转换:'
    
    mjd_time = utc_to_mjd(utc_string)
    write(*, *) '  UTC输入: ', trim(utc_string)
    write(*, *) '  MJD输出: ', mjd_time
    
    ! 验证转换
    utc_output = mjd_to_utc(mjd_time, 'ISOC')
    write(*, *) '  验证UTC: ', trim(utc_output)
    
    ! 测试5: 时间系统转换
    write(*, *) ''
    write(*, *) '5. 测试时间系统转换:'
    
    ! TDB到TAI转换
    call time_transfer(et_time, et_output, 'TDB', 'TAI', status)
    if (status == 0) then
        write(*, *) '  TDB到TAI转换成功'
        write(*, *) '  TDB: ', et_time
        write(*, *) '  TAI: ', et_output
        write(*, *) '  差值: ', et_output - et_time, ' 秒'
    else
        write(*, *) '  TDB到TAI转换失败'
    end if
    
    ! TDB到TDT转换
    call time_transfer(et_time, et_output, 'TDB', 'TDT', status)
    if (status == 0) then
        write(*, *) '  TDB到TDT转换成功'
        write(*, *) '  TDB: ', et_time
        write(*, *) '  TDT: ', et_output
        write(*, *) '  差值: ', et_output - et_time, ' 秒'
    else
        write(*, *) '  TDB到TDT转换失败'
    end if
    
    ! 测试6: 时间字符串验证
    write(*, *) ''
    write(*, *) '6. 测试时间字符串验证:'
    
    ! 有效时间字符串
    is_valid = validate_time_string('2024-01-01T12:00:00')
    write(*, *) '  有效时间字符串: ', is_valid
    
    ! 空字符串
    is_valid = validate_time_string('')
    write(*, *) '  空字符串: ', is_valid
    
    write(*, *) '  注意: 无效时间字符串测试跳过（SPICE错误处理）'
    
    ! 测试7: 多个时间点测试
    write(*, *) ''
    write(*, *) '7. 测试多个时间点:'
    
    do i = 1, 5
        write(utc_string, '(A,I0,A)') '2024-01-', i, 'T12:00:00'
        et_time = utc_to_et(utc_string)
        jd_time = utc_to_jd(utc_string)
        mjd_time = utc_to_mjd(utc_string)
        write(*, *) '  ', trim(utc_string)
        write(*, *) '    ET: ', et_time
        write(*, *) '    JD: ', jd_time
        write(*, *) '    MJD: ', mjd_time
    end do
    
    ! 测试8: 边界条件测试
    write(*, *) ''
    write(*, *) '8. 测试边界条件:'
    
    ! 测试闰年
    utc_string = '2024-02-29T12:00:00'
    et_time = utc_to_et(utc_string)
    write(*, *) '  闰年2月29日: ', trim(utc_string)
    write(*, *) '  ET: ', et_time
    
    ! 测试年末
    utc_string = '2023-12-31T23:59:59'
    et_time = utc_to_et(utc_string)
    write(*, *) '  年末: ', trim(utc_string)
    write(*, *) '  ET: ', et_time
    
    ! 测试9: 当前时间
    write(*, *) ''
    write(*, *) '9. 测试当前时间:'
    
    current_time = get_current_utc('ISOC')
    write(*, *) '  当前UTC时间: ', trim(current_time)
    
    ! 测试10: 精度测试
    write(*, *) ''
    write(*, *) '10. 测试精度:'
    
    utc_string = '2024-01-01T12:00:00.123456'
    et_time = utc_to_et(utc_string)
    utc_output = et_to_utc(et_time, 'ISOC', 6)
    write(*, *) '  高精度输入: ', trim(utc_string)
    write(*, *) '  高精度输出: ', trim(utc_output)
    
    write(*, *) ''
    write(*, *) '=== 时间系统测试完成 ==='
    
    ! 清理系统
    call spice_cleanup()
    call pod_cleanup()
    
end program test_time_system