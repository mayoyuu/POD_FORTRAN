module pod_uq_propagation
    use pod_global, only: DP, MAX_STRING_LEN, output_directory
    use pod_utils,  only: print_separator, print_vector, print_matrix
    
    ! 引入底层计算库
    use pod_uq_base_module,  only: uq_propagator_base
    use pod_uq_da_module,    only: uq_da_propagator
    use pod_uq_mc_module,    only: uq_mc_propagator
    use pod_uq_ut_module,    only: uq_ut_propagator
    use pod_uq_ads_module,   only: uq_ads_propagator
    use pod_uq_state_module, only: uq_state_type
    use pod_uq_hfem_ads_module, only: hfem_ads_options_type, hfem_ads_stats_type
    use pod_uq_ads_coordinates_module, only: ads_coordinate_map_type, &
        ads_build_coordinate_map, ads_physical_to_unit
    use pod_random_module, only: init_random_seed, generate_multivariate_normal, randn
    
    implicit none

    ! 传播算法常量定义 (对外暴露作为开关)
    integer, parameter, public :: METHOD_MC = 1
    integer, parameter, public :: METHOD_DA = 2
    integer, parameter, public :: METHOD_UT = 3
    integer, parameter, public :: METHOD_ADS = 4
    
contains

    ! ====================================================================
    ! 核心 API：执行不确定性传播
    ! 调用者只需传入标称状态、协方差矩阵以及控制开关，即可获得传播前后的完整分布状态
    ! ====================================================================
    subroutine run_uq_propagation(nominal_state, initial_cov, epoch0, t_start, t_end,&
                                  method_switch, n_particles, &
                                  save_results_to_file, initial_state_out, final_state_out, &
                                  da_order,integrator_switch, file_prefix, &
                                  ads_options, ads_stats, initial_perturbations)
        
        real(DP), intent(in) :: nominal_state(6)   ! 标称轨道 (作为均值)
        real(DP), intent(in) :: initial_cov(6,6)   ! 初始协方差矩阵
        real(DP), intent(in) :: epoch0             ! 【新增】物理历元基准 (TDB 秒)
        real(DP), intent(in) :: t_start, t_end     ! 积分起止相对时间
        integer,  intent(in) :: method_switch      ! 方法开关 (METHOD_MC 或 METHOD_DA)
        ! integer,  intent(in) :: integrator_switch  ! 积分器开关 (INTEG_RK4, INTEG_RKF45 等)
        integer,  intent(in) :: n_particles        ! 生成的粒子数量
        logical,  intent(in) :: save_results_to_file ! 是否将结果落盘
        integer,  intent(in), optional :: da_order  ! 可选的 DA 展开阶数 (仅对 DA 方法有效)
        integer,  intent(in), optional :: integrator_switch ! 可选的积分器开关参数
        character(len=*), intent(in), optional :: file_prefix ! 【新增】声明
        type(hfem_ads_options_type), intent(in), optional :: ads_options
        type(hfem_ads_stats_type), intent(out), optional :: ads_stats
        real(DP), intent(in), optional :: initial_perturbations(:,:)
        
        ! 输出参数 (将内存生命周期交还给上层，方便上层继续用于粒子滤波更新)
        type(uq_state_type), intent(out) :: initial_state_out
        type(uq_state_type), intent(out) :: final_state_out
        
        ! 局部变量
        class(uq_propagator_base), allocatable :: propagator
        character(len=MAX_STRING_LEN) :: actual_prefix
        real(DP) :: cpu_start, cpu_end
        integer :: ads_sampled_count, ads_rejected_count

        ads_sampled_count = 0
        ads_rejected_count = 0
        
        call print_separator('不确定性/误差传播 (UQ API)')
        write(*,*) '[INFO] 正在初始化分布参数, 粒子数: ', n_particles
        
        call init_random_seed(.true.) ! 如果你希望每次运行轨迹不同，可以改为 .false.
        if (method_switch == METHOD_ADS) then
            if (present(initial_perturbations)) then
                if (present(ads_options)) then
                    call set_ads_samples_from_perturbations(nominal_state, initial_cov, &
                        ads_options, initial_perturbations, initial_state_out)
                else
                    call set_ads_samples_from_perturbations(nominal_state, initial_cov, &
                        hfem_ads_options_type(), initial_perturbations, initial_state_out)
                end if
                ads_sampled_count = size(initial_perturbations,2)
            else if (present(ads_options)) then
                call generate_ads_samples(nominal_state, initial_cov, n_particles, &
                    ads_options, initial_state_out, ads_sampled_count, ads_rejected_count)
            else
                call generate_ads_samples(nominal_state, initial_cov, n_particles, &
                    hfem_ads_options_type(), initial_state_out, &
                    ads_sampled_count, ads_rejected_count)
            end if
        else
            call initial_state_out%allocate_memory(6, n_particles)
            call generate_multivariate_normal(nominal_state, initial_cov, &
                initial_state_out%samples)
            call initial_state_out%compute_moments()
        end if
        
        ! 2. 实例化对应的传播器
        select case(method_switch)
            case(METHOD_MC)
                allocate(uq_mc_propagator :: propagator)
            case(METHOD_DA)
                allocate(uq_da_propagator :: propagator)
            case(METHOD_UT)
                allocate(uq_ut_propagator :: propagator)
            case(METHOD_ADS)
                allocate(uq_ads_propagator :: propagator)
            case default
                write(*,*) '[ERROR] UQ API: 未知的传播方法开关!'
                return
        end select
        
        ! =======================================================
        ! 3. 配置传播器参数 (核心改动区域)
        ! =======================================================
        propagator%epoch0 = epoch0 ! 【新增】将历元基准注入给动力学传播引擎
        ! 如果用户显式传了积分器开关，才去覆盖默认值；否则使用基类默认的 METHOD_RKF78
        if (present(integrator_switch)) then
            call propagator%set_integrator(integrator_switch)
        end if
        call propagator%set_verbosity(.true.)

        if (present(da_order) .and. method_switch == METHOD_DA) then
            select type (propagator)
                type is (uq_da_propagator)
                    call propagator%set_da_order(da_order)
            end select
        end if
        if (method_switch == METHOD_ADS) then
            select type (propagator)
                type is (uq_ads_propagator)
                    if (present(ads_options)) call propagator%set_ads_options(ads_options)
                    if (present(da_order)) call propagator%set_ads_order(da_order)
            end select
        end if
        
        ! 4. 核心计算
        write(*,*) '--------------------------------------------------'
        write(*,*) ' 开始执行传播计算...'
        write(*,*) ' 使用算法: ', trim(propagator%get_method_name())
        
        call cpu_time(cpu_start)
        call propagator%propagate(t_start, t_end, initial_state_out, final_state_out)
        call cpu_time(cpu_end)
        if (method_switch == METHOD_ADS) then
            select type (propagator)
                type is (uq_ads_propagator)
                    propagator%last_stats%sampled_count = ads_sampled_count
                    propagator%last_stats%rejected_count = ads_rejected_count
                    if (present(ads_stats)) ads_stats = propagator%last_stats
                    if (propagator%last_status /= 0) then
                        write(*,*) '[ERROR] HFEM ADS: ', trim(propagator%last_message)
                        return
                    end if
            end select
        end if
        
        write(*,*) ' ✅ 传播计算完成！耗时: ', cpu_end - cpu_start, ' 秒'
        write(*,*) '--------------------------------------------------'
        
        ! 5. 计算并打印统计信息 (此时 final_state_out 内已经由传播器自动调用过了 compute_moments)
        call display_uq_results(initial_state_out, final_state_out)
        
        ! 6. 文件保存控制
        if (save_results_to_file) then
            ! 如果用户传了前缀就用用户的，没传就给个安全默认值
            if (present(file_prefix)) then
                actual_prefix = file_prefix
            else
                actual_prefix = './output/uq_default' 
            end if
            
            ! 传入前缀进行保存
            call save_uq_results(final_state_out, actual_prefix)
        end if
        
    end subroutine run_uq_propagation

    subroutine generate_ads_samples(nominal_state, covariance, n_particles, &
            options, state, sampled_count, rejected_count)
        real(DP), intent(in) :: nominal_state(6), covariance(6,6)
        integer, intent(in) :: n_particles
        type(hfem_ads_options_type), intent(in) :: options
        type(uq_state_type), intent(inout) :: state
        integer, intent(out) :: sampled_count, rejected_count
        type(ads_coordinate_map_type) :: coordinate_map
        real(DP), allocatable :: candidates(:,:), unit_point(:)
        real(DP) :: eta
        character(len=MAX_STRING_LEN) :: message
        integer :: status, filled, batch_size, j, n_variables

        call ads_build_coordinate_map(nominal_state, covariance, &
            options%coordinate_mode, options%domain_sigma, options%srp_sigma, &
            coordinate_map, status, message)
        if (status /= 0) error stop trim(message)
        if (n_particles < 2) error stop 'ADS requires at least two particles'
        n_variables = coordinate_map%n_variables
        call state%deallocate_memory()
        call state%allocate_memory(n_variables, n_particles)
        allocate(state%mean(n_variables), state%cov(n_variables,n_variables))
        state%mean = 0.0_DP
        state%mean(1:6) = nominal_state
        state%cov = 0.0_DP
        state%cov(1:6,1:6) = covariance
        if (n_variables == 7) state%cov(7,7) = options%srp_sigma**2

        filled = 0
        sampled_count = 0
        do while (filled < n_particles)
            batch_size = max(16, 2*(n_particles-filled))
            allocate(candidates(6,batch_size))
            call generate_multivariate_normal(nominal_state, covariance, candidates)
            do j = 1, batch_size
                sampled_count = sampled_count + 1
                eta = 0.0_DP
                if (n_variables == 7) eta = options%srp_sigma*randn()
                call ads_physical_to_unit(coordinate_map, &
                    candidates(:,j)-nominal_state, eta, unit_point, status)
                if (status /= 0 .or. maxval(abs(unit_point)) > 1.0_DP) cycle
                filled = filled + 1
                state%samples(1:6,filled) = candidates(:,j)
                if (n_variables == 7) state%samples(7,filled) = eta
                if (filled == n_particles) exit
            end do
            deallocate(candidates)
        end do
        rejected_count = sampled_count - n_particles
    end subroutine generate_ads_samples

    subroutine set_ads_samples_from_perturbations(nominal_state, covariance, &
            options, perturbations, state)
        real(DP), intent(in) :: nominal_state(6), covariance(6,6)
        type(hfem_ads_options_type), intent(in) :: options
        real(DP), intent(in) :: perturbations(:,:)
        type(uq_state_type), intent(inout) :: state
        integer :: n_variables

        n_variables = merge(7,6,options%srp_sigma > 0.0_DP)
        if (size(perturbations,1) /= n_variables) &
            error stop 'ADS perturbation file dimension conflicts with srp_sigma'
        if (size(perturbations,2) < 2) &
            error stop 'ADS perturbation file requires at least two samples'
        call state%deallocate_memory()
        call state%allocate_memory(n_variables,size(perturbations,2))
        state%samples = perturbations
        state%samples(1:6,:) = state%samples(1:6,:) + &
            spread(nominal_state,dim=2,ncopies=size(perturbations,2))
        allocate(state%mean(n_variables),state%cov(n_variables,n_variables))
        state%mean = 0.0_DP
        state%mean(1:6) = nominal_state
        state%cov = 0.0_DP
        state%cov(1:6,1:6) = covariance
        if (n_variables == 7) state%cov(7,7) = options%srp_sigma**2
    end subroutine set_ads_samples_from_perturbations
   ! ====================================================================
    ! 专门为粒子滤波 (Particle Filter) 优化的纯粒子传播 API
    ! 直接基于 uq_state_type 对象进行 Array -> Array 的高频调用映射
    ! ====================================================================
    subroutine run_particle_propagation(initial_state, reference_orbit, epoch0, t_start, t_end,&
                                        method_switch, final_state, integrator_switch, da_order,&
                                        reference_orbit_out)
        
        ! 输入/输出参数直接使用 OOP 对象
        type(uq_state_type), intent(inout) :: initial_state ! 传入时内部 samples 已分配并填充
        type(uq_state_type), intent(inout) :: final_state   ! 传播后的状态对象 (底层会自动分配内存)
        
        real(DP), intent(in)  :: reference_orbit(6)     ! DA 需要明确指出的参考轨道
        real(DP), intent(in)  :: epoch0                 ! 物理历元基准 (TDB 秒)
        real(DP), intent(in)  :: t_start, t_end         ! 积分起止相对时间
        integer,  intent(in)  :: method_switch          ! 方法开关 (METHOD_MC 或 METHOD_DA)
        
        integer,  intent(in), optional :: integrator_switch ! 积分器开关 (如 INTEG_RKF45 等)
        integer,  intent(in), optional :: da_order          ! 可选的 DA 阶数
        real(DP), intent(out), optional :: reference_orbit_out(6) ! 返回传播后的参考轨道常数项
        
        ! 局部变量
        class(uq_propagator_base), allocatable :: propagator
        
        ! 1. 实例化对应的传播器
        select case(method_switch)
            case(METHOD_MC)
                allocate(uq_mc_propagator :: propagator)
            case(METHOD_DA)
                allocate(uq_da_propagator :: propagator)
            case(METHOD_UT)
                allocate(uq_ut_propagator :: propagator)
            case default
                write(*,*) '[ERROR] UQ API: 未知的传播方法开关!'
                return
        end select

        ! 2. 配置 DA 阶数 (仅当使用 DA 方法且传入了 da_order 时有效)
        if (present(da_order) .and. method_switch == METHOD_DA) then
            select type (propagator)
                type is (uq_da_propagator)
                    call propagator%set_da_order(da_order)
            end select
        end if
        
        ! 3. 配置参考轨道与传播器参数
        ! 将 reference_orbit 作为均值注入，这是多数 DA 展开的默认基准点
        if (allocated(initial_state%mean)) deallocate(initial_state%mean)
        allocate(initial_state%mean(size(initial_state%samples, 1)))
        initial_state%mean = 0.0_DP
        initial_state%mean(1:6) = reference_orbit
        if (size(initial_state%samples, 1) > 6) then
            initial_state%mean(7:size(initial_state%samples, 1)) = &
                sum(initial_state%samples(7:size(initial_state%samples, 1), :), dim=2) / &
                real(size(initial_state%samples, 2), DP)
        end if
        
        propagator%epoch0 = epoch0
        if (present(integrator_switch)) then
            call propagator%set_integrator(integrator_switch)
        end if
        
        ! 在粒子滤波中高频调用，强制关闭打印输出以提升性能 (同时屏蔽不必要的矩计算)
        call propagator%set_verbosity(.false.)

        ! 对于 UT 方法: 确保协方差矩阵已从样本计算
        select type (prop => propagator)
            type is (uq_ut_propagator)
                if (.not. allocated(initial_state%cov)) then
                    call initial_state%compute_moments()
                end if
        end select

        ! 4. 核心传播计算 (这里会调用内部的 RKF 或其他积分器逻辑)
        ! 底层的 propagate 内部会自动调用 final_state%allocate_memory
        call propagator%propagate(t_start, t_end, initial_state, final_state)
        
        ! 5. 精准提取参考轨道的常数项
        if (present(reference_orbit_out)) then
            select type (prop => propagator)
                type is (uq_da_propagator)
                    ! DA 方法：直接拿刚才缓存好的精确常数项
                    reference_orbit_out = prop%propagated_ref_orbit
                type is (uq_mc_propagator)
                    ! MC 方法：退化为使用样本均值
                    call final_state%compute_moments()
                    reference_orbit_out = final_state%mean
                type is (uq_ut_propagator)
                    reference_orbit_out = final_state%mean
            end select
        end if
        
    end subroutine run_particle_propagation
    
    ! ====================================================================
    ! 以下为内部辅助例程 
    ! ====================================================================
    subroutine display_uq_results(initial_state, final_state)
        type(uq_state_type), intent(in) :: initial_state, final_state
        real(DP) :: std_dev(6)
        integer :: i
        
        call print_separator('传播结果统计摘要')
        call print_vector(initial_state%mean, '初始均值 (Mean_0):', '(6(ES14.5, 1X))')
        
        do i = 1, 6
            std_dev(i) = sqrt(initial_state%cov(i,i))
        end do
        call print_vector(std_dev, '初始标准差 (1-Sigma_0):', '(6(ES14.5, 1X))')
        
        write(*,*) '--------------------------------------------------'
        call print_vector(final_state%mean, '最终均值 (Mean_f):', '(6(ES14.5, 1X))')
        
        do i = 1, 6
            std_dev(i) = sqrt(final_state%cov(i,i))
        end do
        call print_vector(std_dev, '最终标准差 (1-Sigma_f):', '(6(ES14.5, 1X))')
    end subroutine display_uq_results

    subroutine save_uq_results(final_state, file_prefix)
        type(uq_state_type), intent(in) :: final_state
        character(len=*), intent(in) :: file_prefix ! 【新增】用户决定的文件前缀路径
        
        integer :: file_unit, i, dim, n_particles
        character(len=MAX_STRING_LEN) :: filepath_particles, filepath_stats
        
        dim = size(final_state%samples, 1)
        n_particles = size(final_state%samples, 2)
        
        ! 智能拼接路径 (使用用户传入的前缀，避开全局变量中的隐藏换行符)
        filepath_particles = trim(file_prefix) // '_particles.csv'
        filepath_stats = trim(file_prefix) // '_stats.csv'
        
        ! 1. 存粒子
        open(newunit=file_unit, file=trim(filepath_particles), status='replace', action='write')
        if (dim == 7) then
            write(file_unit, '(A)') 'x,y,z,vx,vy,vz,eta_srp'
        else
            write(file_unit, '(A)') 'x,y,z,vx,vy,vz'
        end if
        do i = 1, n_particles
            write(file_unit, '(*(ES22.14, :, ","))') final_state%samples(:, i)
        end do
        close(file_unit)
        
        ! 2. 存统计矩
        open(newunit=file_unit, file=trim(filepath_stats), status='replace', action='write')
        write(file_unit, '(A)') '# Mean'
        write(file_unit, '(*(ES22.14, :, ","))') final_state%mean(:)
        write(file_unit, '(A)') '# Covariance Matrix'
        do i = 1, dim
            write(file_unit, '(*(ES22.14, :, ","))') final_state%cov(i, :)
        end do
        close(file_unit)
        
        write(*,*) '✅ 结果已保存至:'
        write(*,*) '  ', trim(filepath_particles)
        write(*,*) '  ', trim(filepath_stats)
    end subroutine save_uq_results

end module pod_uq_propagation
