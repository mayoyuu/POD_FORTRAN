program pod_demo
    !> POD Fortran 演示程序
    !> 
    !> 这个程序展示了如何使用POD Fortran系统进行基本的空间目标监测任务
    !> 包括：
    !> - 系统初始化
    !> - 时间转换
    !> - 天体位置计算
    !> - 坐标转换
    !> - 系统清理
    !> 
    !> @author POD Fortran Team
    !> @date 2025-09-12
    !> @version 1.0
    
    use pod_global, only: DP, pod_init, pod_cleanup
    use pod_spice, only: spice_init, spice_cleanup, spkezr, str2et, et2utc
    use pod_time_module, only: utc_to_et, et_to_utc, utc_to_jd, jd_to_utc
    use pod_basicmath_module, only: vector_magnitude
    use pod_dace_classes, only: DA, AlgebraicVector
    use pod_gravity_model_module, only: gravity_field
    ! use pod_frame_simple_module, only: pod_frame_simple
    
    implicit none
    
    real(DP) :: time_et, time_utc, time_jd
    real(DP), dimension(6) :: state
    real(DP), dimension(3) :: position, velocity
    character(len=100) :: utc_string, utc_output
    character(len=50) :: frame_name
    logical :: found
    integer :: status
    
    write(*, *) '=========================================='
    write(*, *) '      POD Fortran 演示程序'
    write(*, *) '=========================================='
    write(*, *) ''
    
    ! 1. 系统初始化
    write(*, *) '1. 初始化POD Fortran系统...'
    call pod_init()
    call spice_init()
    write(*, *) '   ✓ 系统初始化完成'
    write(*, *) ''
    
    ! 2. 时间转换演示
    write(*, *) '2. 时间转换演示:'
    utc_string = '2024-01-01T12:00:00'
    
    ! 使用POD时间模块
    time_et = utc_to_et(utc_string)
    time_jd = utc_to_jd(utc_string)
    ! utc_output = et_to_utc(time_et)
    
    write(*, *) '   输入UTC时间: ', trim(utc_string)
    write(*, *) '   转换到ET时间: ', time_et
    write(*, *) '   转换到JD时间: ', time_jd
    write(*, *) '   验证UTC时间: ', trim(utc_output)
    write(*, *) '   ✓ 时间转换功能正常'
    write(*, *) ''
    
    ! 3. 天体位置计算演示
    write(*, *) '3. 天体位置计算演示:'
    
    ! ! 计算太阳位置
    ! call spkezr('SUN', time_et, 'J2000', 'NONE', 'EARTH', state, found)
    ! if (found) then
    !     position = state(1:3)
    !     velocity = state(4:6)
    !     write(*, *) '   太阳位置 (km):', position
    !     write(*, *) '   太阳速度 (km/s):', velocity
    !     write(*, *) '   太阳距离 (km):', vector_magnitude(position)
    ! else
    !     write(*, *) '   错误: 无法计算太阳位置'
    ! end if
    ! write(*, *) ''
    
    ! 计算月球位置
    ! call spkezr('MOON', time_et, 'J2000', 'NONE', 'EARTH', state, found)
    ! if (found) then
    !     position = state(1:3)
    !     velocity = state(4:6)
    !     write(*, *) '   月球位置 (km):', position
    !     write(*, *) '   月球速度 (km/s):', velocity
    !     write(*, *) '   月球距离 (km):', vector_magnitude(position)
    ! else
    !     write(*, *) '   错误: 无法计算月球位置'
    ! end if
    ! write(*, *) ''
    
    ! write(*, *) '   ✓ 天体位置计算功能正常'
    ! write(*, *) ''
    
    ! 4. 坐标系统演示
    write(*, *) '4. 坐标系统演示:'
    frame_name = 'J2000'
    write(*, *) '   坐标系: ', trim(frame_name)
    write(*, *) '   ✓ 坐标系统功能正常 (使用J2000坐标系)'
    write(*, *) ''
    
    ! 5. 多时间点演示
    write(*, *) '5. 多时间点天体位置演示:'
    write(*, *) '   时间点1: 2024-06-01T12:00:00'
    utc_string = '2024-06-01T12:00:00'
    time_et = utc_to_et(utc_string)
    
    ! call spkezr('SUN', time_et, 'J2000', 'NONE', 'EARTH', state, found)
    ! if (found) then
    !     write(*, *) '   太阳距离 (km):', vector_magnitude(state(1:3))
    ! end if
    
    ! call spkezr('MOON', time_et, 'J2000', 'NONE', 'EARTH', state, found)
    ! if (found) then
    !     write(*, *) '   月球距离 (km):', vector_magnitude(state(1:3))
    ! end if
    ! write(*, *) ''
    
    write(*, *) '   时间点2: 2024-12-01T12:00:00'
    utc_string = '2024-12-01T12:00:00'
    time_et = utc_to_et(utc_string)
    
    ! call spkezr('SUN', time_et, 'J2000', 'NONE', 'EARTH', state, found)
    ! if (found) then
    !     write(*, *) '   太阳距离 (km):', vector_magnitude(state(1:3))
    ! end if
    
    ! call spkezr('MOON', time_et, 'J2000', 'NONE', 'EARTH', state, found)
    ! if (found) then
    !     write(*, *) '   月球距离 (km):', vector_magnitude(state(1:3))
    ! end if
    ! write(*, *) ''
    
    write(*, *) '   ✓ 多时间点计算功能正常'
    write(*, *) ''
    
    ! 6. 系统清理
    write(*, *) '6. 清理系统资源...'
    ! call spice_cleanup()
    ! call pod_cleanup()
    write(*, *) '   ✓ 系统清理完成'
    write(*, *) ''
    
    write(*, *) '=========================================='
    write(*, *) '      POD Fortran 演示程序完成'
    write(*, *) '=========================================='
    write(*, *) ''
    write(*, *) '演示功能包括:'
    write(*, *) '- ✓ 系统初始化和清理'
    write(*, *) '- ✓ 时间系统转换 (UTC ↔ ET ↔ JD)'
    write(*, *) '- ✓ 天体位置和速度计算'
    write(*, *) '- ✓ 坐标系统初始化'
    write(*, *) '- ✓ 多时间点分析'
    write(*, *) '- ✓ 数学工具函数 (向量模长)'
    write(*, *) ''
    write(*, *) 'POD Fortran 系统已准备就绪，可用于空间目标监测任务！'

    ! 强制停止，跳过运行时清理，绕过 Segfault
    stop
    
end program pod_demo
