program full_spice_test
    use pod_global, only: DP, pod_init, pod_cleanup
    use pod_spice, only: spice_init, spice_cleanup, spkezr, str2et, et2utc
    use pod_basicmath, only: vector_magnitude
    
    implicit none
    
    real(DP) :: time_et
    real(DP), dimension(6) :: state
    real(DP), dimension(3) :: sun_pos, sun_vel, moon_pos, moon_vel, earth_pos, earth_vel
    character(len=100) :: utc_string, utc_output
    logical :: found
    
    write(*, *) '=== 完整SPICE功能测试 ==='
    
    ! 初始化系统
    call pod_init()
    call spice_init()
    
    ! 测试时间转换
    write(*, *) ''
    write(*, *) '1. 测试时间转换功能:'
    
    ! 测试2024年1月01日
    utc_string = '2024-01-01T12:00:00'
    write(*, *) '输入UTC时间: ', trim(utc_string)
    
    call str2et(utc_string, time_et)
    write(*, *) '转换到ET时间: ', time_et
    
    call et2utc(time_et, 'ISOC', 3, utc_output)
    write(*, *) '转换回UTC时间: ', trim(utc_output)
    
    ! 测试天体位置计算
    write(*, *) ''
    write(*, *) '2. 测试天体位置计算:'
    
    ! 太阳位置
    write(*, *) ''
    write(*, *) '太阳位置 (相对于地球):'
    call spkezr('SUN', time_et, 'J2000', 'NONE', 'EARTH', state, found)
    if (found) then
        sun_pos = state(1:3)
        sun_vel = state(4:6)
        write(*, *) '  位置 (km): ', sun_pos
        write(*, *) '  速度 (km/s): ', sun_vel
        write(*, *) '  距离 (km): ', vector_magnitude(sun_pos)
    else
        write(*, *) '  错误: 无法获取太阳位置'
    end if
    
    ! 月球位置
    write(*, *) ''
    write(*, *) '月球位置 (相对于地球):'
    call spkezr('MOON', time_et, 'J2000', 'NONE', 'EARTH', state, found)
    if (found) then
        moon_pos = state(1:3)
        moon_vel = state(4:6)
        write(*, *) '  位置 (km): ', moon_pos
        write(*, *) '  速度 (km/s): ', moon_vel
        write(*, *) '  距离 (km): ', vector_magnitude(moon_pos)
    else
        write(*, *) '  错误: 无法获取月球位置'
    end if
    
    ! 地球位置（相对于太阳）
    write(*, *) ''
    write(*, *) '地球位置 (相对于太阳):'
    call spkezr('EARTH', time_et, 'J2000', 'NONE', 'SUN', state, found)
    if (found) then
        earth_pos = state(1:3)
        earth_vel = state(4:6)
        write(*, *) '  位置 (km): ', earth_pos
        write(*, *) '  速度 (km/s): ', earth_vel
        write(*, *) '  距离 (km): ', vector_magnitude(earth_pos)
    else
        write(*, *) '  错误: 无法获取地球位置'
    end if
    
    ! 测试不同时间点
    write(*, *) ''
    write(*, *) '3. 测试不同时间点的天体位置:'
    
    ! 测试2024年6月1日
    utc_string = '2024-06-01T12:00:00'
    call str2et(utc_string, time_et)
    write(*, *) '时间: ', trim(utc_string)
    
    call spkezr('SUN', time_et, 'J2000', 'NONE', 'EARTH', state, found)
    if (found) then
        sun_pos = state(1:3)
        write(*, *) '太阳距离 (km): ', vector_magnitude(sun_pos)
    end if
    
    call spkezr('MOON', time_et, 'J2000', 'NONE', 'EARTH', state, found)
    if (found) then
        moon_pos = state(1:3)
        write(*, *) '月球距离 (km): ', vector_magnitude(moon_pos)
    end if
    
    ! 测试2024年12月1日
    utc_string = '2024-12-01T12:00:00'
    call str2et(utc_string, time_et)
    write(*, *) '时间: ', trim(utc_string)
    
    call spkezr('SUN', time_et, 'J2000', 'NONE', 'EARTH', state, found)
    if (found) then
        sun_pos = state(1:3)
        write(*, *) '太阳距离 (km): ', vector_magnitude(sun_pos)
    end if
    
    call spkezr('MOON', time_et, 'J2000', 'NONE', 'EARTH', state, found)
    if (found) then
        moon_pos = state(1:3)
        write(*, *) '月球距离 (km): ', vector_magnitude(moon_pos)
    end if
    
    write(*, *) ''
    write(*, *) '=== SPICE 完整测试完成 ==='
    
    ! 清理系统
    call spice_cleanup()
    call pod_cleanup()
    
end program full_spice_test