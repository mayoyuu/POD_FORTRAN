! program test_da_twobody
!     use pod_global, only: DP
!     use pod_config, only: config
!     use pod_force_model_module, only: set_force_model_options
!     use pod_integrator_module, only: rk4_integrate, adaptive_step_integrate
!     use pod_dace_bridge, only: dace_initialize
!     use pod_dace_classes
!     use pod_da_integrator_module
!     implicit none
    
!     type(AlgebraicVector) :: state_da, new_state_da
!     type(DA) :: delta_var
!     real(DP) :: state_real(6)
!     real(DP) :: state_real_out(6)
!     real(DP) :: dt, current_time
!     integer :: i
    
!     ! --- RKF45 变步长相关变量 ---
!     real(DP) :: t_start, t_end, tolerance
!     integer :: max_steps, n_steps_da, n_steps_real
    
!     ! DA 版本的历史记录
!     real(DP), allocatable :: times_da(:)
!     type(AlgebraicVector), allocatable :: states_da(:)
    
!     ! 实数版本的历史记录
!     real(DP), allocatable :: times_real(:)
!     real(DP), allocatable :: states_real(:,:)
    
!     write(*,*) "================================================"
!     write(*,*) "       DA vs 实数 二体问题轨道传播基准测试"
!     write(*,*) "================================================"
    
!     ! [极其关键的公平性设置]
!     ! 1. 确保纯实数力模型只保留中心引力，关闭 J2、大气等摄动
!     call set_force_model_options(.false., .false., .false., .false.)
!     ! 2. 统一引力常数 (需确保 config 已分配，或你的 pod_config 有默认值)
!     config%gravitational_constant = 398600.4415_DP 
    
!     ! 初始化 DACE 引擎: 阶数 = 2, 变量数 = 6
!     call dace_initialize(2, 6)
    
!     ! 标称初始轨道
!     state_real = [7000.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 7.546_DP, 0.0_DP]
    
!     ! 组装 DA 初始状态
!     call state_da%init(6)
!     call new_state_da%init(6)
!     do i = 1, 6
!         call delta_var%init_var(i) 
!         state_da%elements(i) = state_real(i) + delta_var
!         call delta_var%destroy()
!     end do
    
!     ! ========================================================
!     ! 阶段一：定步长 RK4 对比测试 (单步 60 秒)
!     ! ========================================================
!     write(*,*) ">>> [阶段 1] 正在执行 RK4 积分对比 (步长: 60s)..."
!     dt = 60.0_DP
!     current_time = 0.0_DP
    
!     ! 跑纯实数 RK4
!     call rk4_integrate(state_real, dt, current_time, state_real_out)
    
!     ! 跑 DA RK4
!     call da_rk4_integrate(state_da, dt, current_time, new_state_da)
    
!     write(*,*) "------------------------------------------------"
!     write(*,*) "RK4 积分 60s 后的 X 坐标标称值对比："
!     write(*,*) "  纯实数 RK4 结果 : ", state_real_out(1)
!     write(*,*) "  DA RK4 常数项   : ", new_state_da%elements(1)%cons()
!     write(*,*) "  差异 (绝对误差) : ", abs(state_real_out(1) - new_state_da%elements(1)%cons())
    
!     ! ========================================================
!     ! 阶段二：自适应变步长 RKF45 对比测试
!     ! ========================================================
!     write(*,*) ""
!     write(*,*) ">>> [阶段 2] 正在执行 RKF45 变步长积分对比..."
!     t_start = 0.0_DP
!     t_end = 3600.0_DP       ! 测试跑 1 个小时
!     max_steps = 2000        ! 把 100000 改为 2000
!     tolerance = 1.0e-8_DP   ! 保持 1e-8 的合理容差
    
!     ! 跑纯实数 RKF45
!     call adaptive_step_integrate(state_real, t_start, t_end, max_steps, tolerance, &
!                                  times_real, states_real, n_steps_real)
                                 
!     ! 跑 DA RKF45
!     call da_adaptive_step_integrate(state_da, t_start, t_end, max_steps, tolerance, &
!                                     times_da, states_da, n_steps_da)
                                    
!     write(*,*) "------------------------------------------------"
!     write(*,*) "RKF45 积分 1 小时 统计数据对比："
!     write(*,*) "  纯实数 RKF45 步数 : ", n_steps_real
!     write(*,*) "  DA RKF45 步数     : ", n_steps_da
    
!     if (n_steps_real > 0 .and. n_steps_da > 0) then
!         write(*,*) "------------------------------------------------"
!         write(*,*) "RKF45 积分 1 小时 后的最终 X 坐标对比："
!         write(*,*) "  纯实数 RKF45 结果 : ", states_real(n_steps_real, 1)
!         write(*,*) "  DA RKF45 常数项   : ", states_da(n_steps_da)%elements(1)%cons()
!         write(*,*) "  差异 (绝对误差)   : ", abs(states_real(n_steps_real, 1) - states_da(n_steps_da)%elements(1)%cons())
!     end if
!     write(*,*) "================================================"

!     ! 内存清理
!     call state_da%destroy()
!     call new_state_da%destroy()
    
!     if (allocated(states_da)) then
!         do i = 1, size(states_da)
!             call states_da(i)%destroy()
!         end do
!         deallocate(states_da)
!     end if
!     if (allocated(times_da)) deallocate(times_da)
!     if (allocated(states_real)) deallocate(states_real)
!     if (allocated(times_real)) deallocate(times_real)

! end program test_da_twobody