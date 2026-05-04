!> @file pod_ut_runner_module.f90
!> @brief 无迹变换（UT）轨道确定集成封装层
module pod_ut_runner_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_filter_ut_module, only: ut_filter          ! 你的 UT 滤波器模块
    use pod_obs_io_module, only: load_single_observation
    use pod_measurement_base_module, only: observation_station
    use pod_basicmath_module, only: PI
    use pod_data_format_module, only: load_initial_opm, write_json_opm

    implicit none
    private

    public :: run_ut_orbit_determination

contains

    !> ======================================================================
    !> 执行基于无迹卡尔曼滤波（UKF）的轨道确定流程
    !> ======================================================================
    subroutine run_ut_orbit_determination(obs_file, site_json_file, initial_json_file, &
                                          output_json_file, kappa)
        character(len=*), intent(in) :: obs_file           ! 观测文件 (.obs)
        character(len=*), intent(in) :: site_json_file     ! 测站配置文件 (.json)
        character(len=*), intent(in) :: initial_json_file  ! 初始先验状态文件 (.json)
        character(len=*), intent(in) :: output_json_file   ! 输出结果文件 (.json)
        real(DP),         intent(in) :: kappa              ! UT 参数，通常取值 3-n 或 0

        ! 核心对象
        type(ut_filter)       :: my_filter
        type(observation_station) :: current_station

        ! 状态与时间
        real(DP) :: initial_mean(6), final_mean(6)
        real(DP) :: initial_cov(6,6), final_cov(6,6), noise_Q(6,6)
        real(DP) :: y_meas(2), noise_R(2,2)
        real(DP) :: et_current, et_obs, dt
        integer  :: obs_count, i
        logical  :: is_eof

        ! ---- 1. 测量噪声协方差（光学赤经赤纬，0.1角秒精度）----
        noise_R = 0.0_DP
        noise_R(1,1) = (0.1_DP * PI / 180.0_DP / 3600.0_DP)**2
        noise_R(2,2) = noise_R(1,1)

        ! ---- 2. 从 JSON 加载初始状态（仅均值和协方差）----
        call load_initial_opm(initial_json_file, et_current, initial_mean, initial_cov)

        ! ---- 3. 初始化 UT 滤波器 ----
        call my_filter%filter_init(initial_mean, initial_cov, kappa)

        write(*,*) '  [UT Runner] 滤波器初始化完成'
        write(*,*) '    初始历元：', et_current
        write(*,*) '    初始均值：', initial_mean(1:3)

        ! ---- 4. 主循环：逐观测进行时间更新 + 测量更新 ----
        obs_count = 1
        do
            ! 读取单次观测
            call load_single_observation(obs_file, site_json_file, obs_count, &
                                         et_obs, y_meas(1), y_meas(2), current_station, is_eof)
            if (is_eof) exit

            ! 当前滤波器时刻
            call my_filter%get_current_epoch(et_current)
            dt = et_obs - et_current

            ! 构造与时间步长相关的过程噪声 Q（简单白噪声模型）
            noise_Q = 0.0_DP
            do i = 1, 3
                noise_Q(i,i)   = (0.5_DP * dt**2 * 1.0e-6_DP)**2   ! 位置噪声
                noise_Q(i+3,i+3) = (dt * 1.0e-6_DP)**2             ! 速度噪声
            end do

            write(*,'(A,I0,A,F10.2,A)') '  [UT Runner] 处理观测 #', obs_count, &
                  '  时间步长 dt =', dt, ' 秒'

            ! 时间更新（传播 + 过程噪声叠加）
            call my_filter%time_update(et_obs, noise_Q)

            ! 测量更新（利用实际观测修正）
            call my_filter%measurement_update(y_meas, noise_R, et_obs, current_station)

            obs_count = obs_count + 1
        end do

        write(*,*) '  [UT Runner] 滤波结束，有效观测数：', obs_count - 1

        ! ---- 5. 提取最终结果并写入 JSON 文件 ----
        call my_filter%get_current_epoch(et_current)
        call my_filter%get_current_state(final_mean)
        call my_filter%get_current_cov(final_cov)

        ! 使用不带 GMM 状态参数的 write_json_opm 重载版本（若库中提供）
        ! 若库中仅有带 gmm_state 的版本，可改为：
        !   call write_json_opm(output_json_file, final_mean, final_cov, et_current)
        call write_json_opm(output_json_file, final_mean, final_cov, et_current)

        write(*,*) '  [UT Runner] 结果已写入：', trim(output_json_file)

    end subroutine run_ut_orbit_determination

end module pod_ut_runner_module