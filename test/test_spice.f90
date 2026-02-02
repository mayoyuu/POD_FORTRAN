program test_spice
    use cat_global, only: DP, cat_init, cat_cleanup
    use cat_spice, only: spice_init, spice_cleanup, spkezr, str2et, et2utc
    use cat_basicmath, only: vector_magnitude
    
    implicit none
    
    real(DP) :: time_et, time_utc
    real(DP), dimension(6) :: state
    real(DP), dimension(3) :: sun_pos, sun_vel, moon_pos, moon_vel, earth_pos, earth_vel
    character(len=100) :: utc_string
    logical :: found
    
    write(*, *) '=== SPICE 功能测试 ==='
    
    ! 初始化系统
    call cat_init()
    call spice_init()
    
    ! 测试时间转换
    write(*, *) ''
    write(*, *) '1. 测试时间转换功能:'
    utc_string = '2024-01-01T12:00:00'
    call str2et(utc_string, time_et)
    write(*, *) 'UTC时间: ', trim(utc_string)
    write(*, *) 'ET时间: ', time_et
    
    call et2utc(time_et, 'ISOC', 3, utc_string)
    write(*, *) '转换回UTC: ', trim(utc_string)
    
    ! 测试天体位置计算
    write(*, *) ''
    write(*, *) '2. 测试天体位置计算:'
    
    ! 太阳位置
    call spkezr('SUN', time_et, 'J2000', 'NONE', 'EARTH', state, found)
    if (found) then
        sun_pos = state(1:3)
        sun_vel = state(4:6)
        write(*, *) '太阳位置 (km): ', sun_pos
        write(*, *) '太阳速度 (km/s): ', sun_vel
        write(*, *) '太阳距离 (km): ', vector_magnitude(sun_pos)
    end if
    
    ! 月球位置
    call spkezr('MOON', time_et, 'J2000', 'NONE', 'EARTH', state, found)
    if (found) then
        moon_pos = state(1:3)
        moon_vel = state(4:6)
        write(*, *) '月球位置 (km): ', moon_pos
        write(*, *) '月球速度 (km/s): ', moon_vel
        write(*, *) '月球距离 (km): ', vector_magnitude(moon_pos)
    end if
    
    ! 地球位置（相对于太阳）
    call spkezr('EARTH', time_et, 'J2000', 'NONE', 'SUN', state, found)
    if (found) then
        earth_pos = state(1:3)
        earth_vel = state(4:6)
        write(*, *) '地球位置 (km): ', earth_pos
        write(*, *) '地球速度 (km/s): ', earth_vel
        write(*, *) '地球距离 (km): ', vector_magnitude(earth_pos)
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
    write(*, *) '=== SPICE 测试完成 ==='
    
    ! 清理系统
    call spice_cleanup()
    call cat_cleanup()
    
end program test_spice
