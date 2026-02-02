program direct_spice_test
    implicit none
    
    ! 直接声明SPICE函数接口
    interface
        subroutine furnsh(kernel)
            character(len=*), intent(in) :: kernel
        end subroutine furnsh
        
        subroutine unload(kernel)
            character(len=*), intent(in) :: kernel
        end subroutine unload
        
        subroutine str2et(time_string, et)
            character(len=*), intent(in) :: time_string
            real(8), intent(out) :: et
        end subroutine str2et
        
        subroutine et2utc(et, format, precision, utc_string)
            real(8), intent(in) :: et
            character(len=*), intent(in) :: format
            integer, intent(in) :: precision
            character(len=*), intent(out) :: utc_string
        end subroutine et2utc
        
        subroutine spkezr(target, et, ref, abcorr, observer, state, lt)
            character(len=*), intent(in) :: target
            real(8), intent(in) :: et
            character(len=*), intent(in) :: ref
            character(len=*), intent(in) :: abcorr
            character(len=*), intent(in) :: observer
            real(8), dimension(6), intent(out) :: state
            real(8), intent(out) :: lt
        end subroutine spkezr
    end interface
    
    ! 变量声明
    real(8) :: time_et, light_time
    real(8), dimension(6) :: state
    character(len=100) :: utc_string
    integer :: i
    
    write(*, *) '=== 直接SPICE库测试 ==='
    
    ! 1. 加载SPICE内核
    write(*, *) ''
    write(*, *) '1. 加载SPICE内核:'
    
    call furnsh('./kernels/lsk/naif0012.tls')
    write(*, *) '已加载闰秒内核: ./kernels/lsk/naif0012.tls'
    
    call furnsh('./kernels/pck/pck00010.tpc')
    write(*, *) '已加载行星内核: ./kernels/pck/pck00010.tpc'
    
    call furnsh('./kernels/spk/de421.bsp')
    write(*, *) '已加载星历内核: ./kernels/spk/de421.bsp'
    
    call furnsh('./kernels/pck/earth_000101_230801_230509.bpc')
    write(*, *) '已加载地球内核: ./kernels/pck/earth_000101_230801_230509.bpc'
    
    ! 2. 测试时间转换
    write(*, *) ''
    write(*, *) '2. 测试时间转换:'
    
    utc_string = '2024-01-01T12:00:00'
    call str2et(utc_string, time_et)
    write(*, *) 'UTC时间: ', trim(utc_string)
    write(*, *) 'ET时间: ', time_et
    
    ! 转换回UTC
    call et2utc(time_et, 'ISOC', 0, utc_string)
    write(*, *) '转换回UTC: ', trim(utc_string)
    
    ! 3. 测试天体位置计算
    write(*, *) ''
    write(*, *) '3. 测试天体位置计算:'
    
    ! 太阳位置
    call spkezr('SUN', time_et, 'J2000', 'NONE', 'EARTH', state, light_time)
    write(*, *) '太阳位置 (km): ', state(1:3)
    write(*, *) '太阳速度 (km/s): ', state(4:6)
    write(*, *) '太阳距离 (km): ', sqrt(sum(state(1:3)**2))
    
    ! 月球位置
    call spkezr('MOON', time_et, 'J2000', 'NONE', 'EARTH', state, light_time)
    write(*, *) '月球位置 (km): ', state(1:3)
    write(*, *) '月球速度 (km/s): ', state(4:6)
    write(*, *) '月球距离 (km): ', sqrt(sum(state(1:3)**2))
    
    ! 地球位置（相对于太阳）
    call spkezr('EARTH', time_et, 'J2000', 'NONE', 'SUN', state, light_time)
    write(*, *) '地球位置 (km): ', state(1:3)
    write(*, *) '地球速度 (km/s): ', state(4:6)
    write(*, *) '地球距离 (km): ', sqrt(sum(state(1:3)**2))
    
    ! 4. 测试不同时间点
    write(*, *) ''
    write(*, *) '4. 测试不同时间点:'
    
    ! 测试2024年6月1日
    utc_string = '2024-06-01T12:00:00'
    call str2et(utc_string, time_et)
    write(*, *) '时间: ', trim(utc_string)
    
    call spkezr('SUN', time_et, 'J2000', 'NONE', 'EARTH', state, light_time)
    write(*, *) '太阳距离 (km): ', sqrt(sum(state(1:3)**2))
    
    call spkezr('MOON', time_et, 'J2000', 'NONE', 'EARTH', state, light_time)
    write(*, *) '月球距离 (km): ', sqrt(sum(state(1:3)**2))
    
    ! 测试2024年12月1日
    utc_string = '2024-12-01T12:00:00'
    call str2et(utc_string, time_et)
    write(*, *) '时间: ', trim(utc_string)
    
    call spkezr('SUN', time_et, 'J2000', 'NONE', 'EARTH', state, light_time)
    write(*, *) '太阳距离 (km): ', sqrt(sum(state(1:3)**2))
    
    call spkezr('MOON', time_et, 'J2000', 'NONE', 'EARTH', state, light_time)
    write(*, *) '月球距离 (km): ', sqrt(sum(state(1:3)**2))
    
    write(*, *) ''
    write(*, *) '=== 直接SPICE测试完成 ==='
    
    ! 清理（可选，避免segfault）
    ! call unload('*')
    
end program direct_spice_test
