program test_frame_oop
    use pod_global, only: DP
    use pod_spice, only: spice_init, str2et  ! 主程序只需负责初始化和定个时间
    use pod_frame_simple_module, only: pod_frame_state
    use pod_frame_module, only: transform_frame_state  ! <--- 直接调用最高层的封装
    
    implicit none
    
    ! 声明面向对象的状态容器
    type(pod_frame_state) :: state_in, state_out
    real(DP) :: time_et
    
    write(*, *) '=== 终极架构验证: 坐标系高层封装测试 ==='
    
    ! 1. 初始化 SPICE 环境 (加载内核库)
    call spice_init()
    
    ! 2. 设定测试时间 (注意这里我们用了 2023 年以规避之前的报错)
    call str2et('2023-01-01T12:00:00', time_et)
    
    ! 3. 准备输入数据: GCRS (J2000) 下的卫星状态
    call state_in%set_frame_name('J2000')
    call state_in%set_epoch(time_et)
    call state_in%set_position([7000.0_DP, 0.0_DP, 0.0_DP])
    call state_in%set_velocity([0.0_DP, 7.5_DP, 0.0_DP])
    call state_in%set_valid(.true.)
    
    write(*, *) ''
    write(*, *) '--- 转换前 ---'
    write(*, *) '当前坐标系: ', trim(state_in%get_frame_name())
    write(*, '(A, 3F15.4)') ' 位置 (km)   : ', state_in%get_position()
    write(*, '(A, 3F15.4)') ' 速度 (km/s) : ', state_in%get_velocity()
    
    ! 4. 见证奇迹的时刻：一键调用高层转换接口！
    ! 没有任何复杂的矩阵乘法，没有任何底层的 SPICE 函数暴露在这里
    call transform_frame_state(state_in, 'ITRF93', state_out)
    
    ! 5. 打印结果
    write(*, *) ''
    if (state_out%is_valid()) then
        write(*, *) '--- 转换后 ---'
        write(*, *) '目标坐标系: ', trim(state_out%get_frame_name())
        write(*, '(A, 3F15.4)') ' 位置 (km)   : ', state_out%get_position()
        write(*, '(A, 3F15.4)') ' 速度 (km/s) : ', state_out%get_velocity()
    else
        write(*, *) '[ERROR] 转换失败！请检查内核文件或时间覆盖范围。'
    end if

end program test_frame_oop