program test_orbit_output
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_config, only: config
    use pod_engine_module, only: pod_engine_init
    use pod_spice, only: str2et, et2utc
    use pod_data_format_module, only: load_initial_opm
    use pod_orbit_propagation, only: orbit_state, propagation_result, &
                                   propagate_orbit, cleanup_propagation_result
    
    implicit none

    ! ---- 变量定义 ----
    character(len=*), parameter :: OPM_FILE = 'input/FAC_2604_1.opm'
    character(len=*), parameter :: OUTPUT_FILE = 'output/orbit_integration_results_4030_0808_without_srp.txt'
    character(len=*), parameter :: FILTER_TIME_STR = '2026-04-30T20:00:00.000'
    character(len=*), parameter :: FINAL_TIME_STR  = '2026-08-08T07:20:00.000'
    
    real(DP) :: et_start, et_final, et_filter, propagation_duration
    ! 【修复2：将循环内的变量声明移动到最上方】
    real(DP) :: current_et 
    
    real(DP) :: initial_state_vec(6), initial_cov(6,6)
    type(orbit_state) :: start_state
    type(propagation_result) :: result
    
    integer :: i, u_out, ios
    character(len=32) :: utc_str
    
    ! ==========================================
    ! 1. 系统与环境初始化
    ! ==========================================
    ! 初始化配置、SPICE 内核等 (如果你的配置文件名不同，请自行修改)
    call pod_engine_init('config/dummy_test_config.txt') 
    
    ! 解析时间阈值和最终时间
    call str2et(FILTER_TIME_STR, et_filter)
    call str2et(FINAL_TIME_STR, et_final)
    
    ! ==========================================
    ! 2. 读取 OPM 初始状态
    ! ==========================================
    write(*,*) ">>> 正在从 OPM 读取初始轨道数据..."
    call load_initial_opm(OPM_FILE, et_start, initial_state_vec, initial_cov)
    
    start_state%state = initial_state_vec
    start_state%epoch = et_start
    
    ! 计算需要传播的总时长 (秒)
    propagation_duration = et_final - et_start
    
    if (propagation_duration <= 0.0_DP) then
        write(*,*) "[ERROR] 积分终止时间早于或等于初始历元！"
        stop
    end if
    
    ! ==========================================
    ! 3. 执行轨道传播 (积分)
    ! ==========================================
    write(*,*) ">>> 开始积分传播至 ", FINAL_TIME_STR
    ! 使用推荐的 RKF78 积分器 (integrator_choice = 2)
    call propagate_orbit(start_state, propagation_duration, 2, result)
    write(*,*) ">>> 积分完成，总步数: ", result%n_steps
    
    ! ==========================================
    ! 4. 格式化输出到 TXT 文件
    ! ==========================================
    write(*,*) ">>> 正在保存符合时间阈值的结果到: ", OUTPUT_FILE
    open(newunit=u_out, file=OUTPUT_FILE, status='replace', action='write', iostat=ios)
    if (ios /= 0) stop "[ERROR] 无法创建输出文件"
    
    ! 【修复1：使用 // 字符串拼接符安全折行，避免长度超限报警】
    write(u_out, '(A)') '# UTC_Time                 X(km)                 ' // &
                        'Y(km)                 Z(km)               ' // &
                        'VX(km/s)              VY(km/s)              VZ(km/s)'
    
    do i = 1, result%n_steps
        ! 计算当前步的绝对历元 (ET) = 初始历元 + 相对积分时间
        current_et = et_start + result%times(i)
        
        ! 判定是否在 2026-04-30T20:00:00 之后
        if (current_et >= et_filter) then
            ! 将 ET 转换为 UTC 字符串 (SPICE 格式)
            call et2utc(current_et, 'ISOC', 3, utc_str)
            
            ! 按照指定格式写入：UTC时间(A24) + 6个状态分量(ES24.14)
            write(u_out, '(A24, 2X, 6(ES24.14, 2X))') &
                trim(utc_str), &
                result%states(i, 1), result%states(i, 2), result%states(i, 3), &
                result%states(i, 4), result%states(i, 5), result%states(i, 6)
        end if
    end do
    
    close(u_out)
    
    ! 清理内存
    call cleanup_propagation_result(result)
    
    write(*,*) ">>> 处理完毕！"

end program test_orbit_output