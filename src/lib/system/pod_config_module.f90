!> # POD Configuration Management Module
!> 
!> This module provides comprehensive configuration management for the POD Fortran
!> space object monitoring system. It handles loading, saving, and validation
!> of configuration parameters including numerical methods, output formats,
!> and system settings.
!> 
!> ## Features
!> 
!> - **Configuration Loading**: Load configuration from text files
!> - **Parameter Validation**: Comprehensive parameter validation and range checking
!> - **Default Configuration**: Automatic generation of default configuration templates
!> - **Configuration Display**: Formatted display of all configuration parameters
!> - **Parameter Access**: Get and set individual configuration parameters
!> 
!> ## Configuration PODegories
!> 
!> - **Orbit Propagation**: Step sizes, tolerances, integrator parameters
!> - **Orbit Improvement**: Convergence criteria, iteration limits, damping
!> - **Coordinate Systems**: Reference frames, time systems, corrections
!> - **Output Control**: Formats, precision, file loPODions
!> - **Numerical Methods**: Tolerances, matrix sizes, precision
!> - **Performance**: Threading, memory limits, GPU acceleration
!> - **Force Models**: Perturbation models, atmospheric drag, radiation pressure
!> - **Measurement Models**: Noise models, weighting, outlier detection
!> - **File I/O**: Directories, backup settings, file formats
!> - **Debugging**: Debug levels, intermediate results, logging
!> 
!> ## Dependencies
!> 
!> - `pod_global`: For data types and global state management
!> 
!> ## Input/Output Files
!> 
!> - **Input**: `config/pod_config.txt`
!> - **Output**: `config/pod_config.txt` (default template)
!> 
!> ## Author
!> 
!> Zhao Yuhui (PMO, zhaoyuhui@pmo.ac.cn)
!> 
!> ## Version
!> 
!> - **Created**: 2025-09-08
!> - **Ref**: [待添加参考文献]
!> 
!> @note This module provides 50+ configuration parameters covering all aspects
!>       of the POD Fortran system operation.
!> 
!> @warning Configuration validation is performed automatically when loading
!>          configuration files. Invalid parameters will be reported and
!>          automatically corrected where possible.

module pod_config
    use pod_global, only: DP, MAX_STRING_LEN, global_state
    
    implicit none
    
    !> Configuration parameters structure
    !!
    !! This derived type contains all configuration parameters for the POD Fortran
    !! system, organized into logical PODegories for easy management and validation.
    type config_params
        ! 轨道传播参数
        ! 归一化单位
        real(DP) :: LU = 384400.0_DP
        real(DP) :: TU = 375190.26189464843_DP
        real(DP) :: MU = 6.045626292270688E+24_DP
        real(DP) :: VU = 1.024546847401753_DP
        real(DP) :: AccU = 2.730739444643263E-06_DP
        !积分器参数
        real(DP) :: propagation_step = 60.0_DP  ! 传播步长（秒）
        real(DP) :: propagation_abs_tol = 1.0e-12_DP  ! 积分器容差
        real(DP) :: propagation_rel_tol = 1.0e-12_DP  ! 积分器相对容差
        integer :: max_propagation_steps = 10000  ! 最大传播步数
        real(DP) :: rk4_step_size = 60.0_DP  ! RK4积分器步长
        real(DP) :: rkf45_abs_tol = 1.0e-12_DP  ! RKF45积分器容差
        real(DP) :: rkf45_rel_tol = 1.0e-12_DP  ! RKF45积分器相对容差
        real(DP) :: rkf45_min_step = 1.0e-6_DP  ! RKF45最小步长
        real(DP) :: rkf45_max_step = 3600.0_DP  ! RKF45最大步长
        real(DP) :: rkf78_abs_tol = 1.0e-12_DP  ! RKF78积分器容差
        real(DP) :: rkf78_rel_tol = 1.0e-12_DP  ! RKF78积分器相对容差
        real(DP) :: rkf78_min_step = 1.0e-6_DP  ! RKF78最小步长
        real(DP) :: rkf78_max_step = 3600.0_DP  ! RKF78最大步长

        ! 轨道改进参数
        real(DP) :: convergence_tolerance = 1.0e-8_DP  ! 收敛容差
        integer :: max_iterations = 100  ! 最大迭代次数
        real(DP) :: damping_factor = 0.1_DP  ! 阻尼因子
        logical :: use_adaptive_damping = .true.  ! 使用自适应阻尼
        
        ! 坐标系统参数
        character(len=MAX_STRING_LEN) :: reference_frame = 'J2000'
        character(len=MAX_STRING_LEN) :: central_body = 'EARTH'
        character(len=MAX_STRING_LEN) :: time_system = 'UTC'
        logical :: use_precession = .true.  ! 使用岁差修正
        logical :: use_nutation = .true.  ! 使用章动修正
        
        ! 输出参数
        logical :: verbose_output = .true.
        logical :: save_results = .true.
        character(len=MAX_STRING_LEN) :: output_format = 'CSV'
        character(len=MAX_STRING_LEN) :: date_format = 'YYYY-MM-DD HH:MM:SS'
        character(len=MAX_STRING_LEN) :: number_format = 'F20.12'
        logical :: scientific_notation = .false.
        integer :: output_precision = 12  ! 输出精度位数
        
        ! 数值计算参数
        real(DP) :: machine_epsilon = epsilon(1.0_DP)
        real(DP) :: gravitational_constant = 398600.4418_DP  ! km³/s² (地球)
        real(DP) :: numerical_tolerance = 1.0e-15_DP  ! 数值计算容差
        integer :: max_matrix_size = 1000  ! 最大矩阵大小
        
        ! 性能参数
        integer :: max_threads = 1  ! 最大线程数
        logical :: parallel_computation = .false.  ! 并行计算
        real(DP) :: memory_limit = 1.0e9_DP  ! 内存限制 (bytes)
        logical :: use_gpu = .false.  ! 使用GPU加速
        
        ! 力模型参数
        ! =====================================
        ! N体与高阶引力场参数 (Cislunar 特化)
        ! =====================================
        logical :: use_earth_nspheric = .true.
        logical :: use_moon_nspheric  = .true.
        logical :: use_third_body     = .true.
        logical :: use_srp            = .true.
        logical :: use_drag           = .false.
        logical :: use_relativity     = .false.
        
        integer :: earth_degree = 10
        integer :: moon_degree  = 10
        
        ! 激活的天体列表 (1:水星, ..., 3:地球, ..., 10:月球, 11:太阳)
        logical, dimension(11) :: use_planet = .false.
        
        ! 测量模型参数
        real(DP) :: measurement_noise_std = 1.0e-3_DP  ! 测量噪声标准差
        logical :: use_measurement_weights = .true.  ! 使用测量权重
        real(DP) :: outlier_threshold = 3.0_DP  ! 异常值阈值
        
        ! 文件I/O参数
        character(len=MAX_STRING_LEN) :: data_directory = './data/'
        character(len=MAX_STRING_LEN) :: results_directory = './results/'
        logical :: auto_backup = .true.  ! 自动备份
        integer :: backup_interval = 3600  ! 备份间隔 (秒)
        
        ! 调试参数
        logical :: debug_mode = .false.  ! 调试模式
        integer :: debug_level = 1  ! 调试级别
        logical :: save_intermediate_results = .false.  ! 保存中间结果
        character(len=MAX_STRING_LEN) :: debug_directory = './debug/'
    end type config_params
    
    type(config_params) :: config
    
contains

    !> Load configuration from file
    !!
    !! This subroutine loads configuration parameters from a text file.
    !! The file should contain key-value pairs in the format "key = value".
    !! Comments (lines starting with #) and empty lines are ignored.
    !!
    !! @param[in] config_file Path to the configuration file
    !!
    !! @note If the configuration file does not exist, a default configuration
    !!       will be created and used.
    !!
    !! @warning Invalid parameter values will be reported but the loading
    !!          process will continue with default values.
    subroutine load_config(config_file)
        character(len=*), intent(in) :: config_file
        integer :: unit, ios, line_num
        character(len=MAX_STRING_LEN) :: line, key, value
        
        ! 设置默认配置
        call set_default_config()
        
        ! 检查配置文件是否存在
        if (.not. file_exists(config_file)) then
            write(*, *) '配置文件不存在，使用默认配置'
            call save_default_config(config_file)
            return
        end if
        
        ! 打开配置文件
        open(newunit=unit, file=config_file, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*, *) '无法打开配置文件: ', trim(config_file)
            return
        end if
        
        line_num = 0
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            
            line_num = line_num + 1
            
            ! 跳过空行和注释行
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            
            ! 解析键值对
            call parse_key_value(line, key, value)
            if (len_trim(key) > 0) then
                call set_config_value(key, value)
            end if
        end do
        
        close(unit)
        
        global_state%config_loaded = .true.
        write(*, *) '配置文件加载完成: ', trim(config_file)
    end subroutine load_config
    
    subroutine set_default_config()
        ! 轨道传播参数
        config%LU = 384400.0_DP
        config%TU = 375190.26189464843_DP
        config%MU = 6.045626292270688E+24_DP
        config%VU = 1.024546847401753_DP
        config%AccU = 2.730739444643263E-06_DP

        config%propagation_step = 60.0_DP
        config%propagation_abs_tol = 1.0e-12_DP
        config%propagation_rel_tol = 1.0e-12_DP
        config%max_propagation_steps = 10000
        config%rk4_step_size = 60.0_DP
        config%rkf45_abs_tol = 1.0e-12_DP
        config%rkf45_rel_tol = 1.0e-12_DP
        config%rkf45_min_step = 1.0e-6_DP
        config%rkf45_max_step = 3600.0_DP
        config%rkf78_abs_tol = 1.0e-12_DP
        config%rkf78_rel_tol = 1.0e-12_DP
        config%rkf78_min_step = 1.0e-6_DP
        config%rkf78_max_step = 3600.0_DP
        
        ! 轨道改进参数
        config%convergence_tolerance = 1.0e-8_DP
        config%max_iterations = 100
        config%damping_factor = 0.1_DP
        config%use_adaptive_damping = .true.
        
        ! 坐标系统参数
        config%reference_frame = 'J2000'
        config%central_body = 'EARTH'
        config%time_system = 'UTC'
        config%use_precession = .true.
        config%use_nutation = .true.
        
        ! 输出参数
        config%verbose_output = .true.
        config%save_results = .true.
        config%output_format = 'CSV'
        config%date_format = 'YYYY-MM-DD HH:MM:SS'
        config%number_format = 'F20.12'
        config%scientific_notation = .false.
        config%output_precision = 12
        
        ! 数值计算参数
        config%machine_epsilon = epsilon(1.0_DP)
        config%gravitational_constant = 398600.4418_DP
        config%numerical_tolerance = 1.0e-15_DP
        config%max_matrix_size = 1000
        
        ! 性能参数
        config%max_threads = 1
        config%parallel_computation = .false.
        config%memory_limit = 1.0e9_DP
        config%use_gpu = .false.
        
        ! 力模型参数
        config%use_earth_nspheric = .true.
        config%use_moon_nspheric = .true.
        config%use_third_body = .true.
        config%use_srp = .true.
        config%use_drag = .false.
        config%use_relativity = .false.
        config%earth_degree = 10
        config%moon_degree = 10
        config%use_planet = .false.
        config%use_planet(3)  = .true.  ! 默认给地球
        config%use_planet(10) = .true.  ! 默认给月球
        config%use_planet(11) = .true.  ! 默认给太阳
  
        ! 测量模型参数
        config%measurement_noise_std = 1.0e-3_DP
        config%use_measurement_weights = .true.
        config%outlier_threshold = 3.0_DP
        
        ! 文件I/O参数
        config%data_directory = './data/'
        config%results_directory = './results/'
        config%auto_backup = .true.
        config%backup_interval = 3600
        
        ! 调试参数
        config%debug_mode = .false.
        config%debug_level = 1
        config%save_intermediate_results = .false.
        config%debug_directory = './debug/'
        call resolve_config_dependencies()

        write(*, *) '默认配置已设置'
    end subroutine set_default_config

    ! 配置依赖决议与修正
    subroutine resolve_config_dependencies()
        ! 1. 月球高阶场依赖检查
        if (config%use_moon_nspheric) then
            ! 如果第三体总开关没开，或者月球(10)没有被激活
            if (.not. config%use_third_body .or. .not. config%use_planet(10)) then
                write(*, *) '⚠️ 警告: 启用了月球高阶场，但未激活第三体摄动或未包含月球(10)。'// &
                '已自动强制禁用月球高阶场！'
                config%use_moon_nspheric = .false.
            end if
        end if
        
        ! 2. 地球高阶场依赖检查 (以防万一有人把地球中心引力给关了)
        if (config%use_earth_nspheric) then
            if (.not. config%use_planet(3)) then
                write(*, *) '⚠️ 警告: 启用了地球高阶场，但未激活地球(3)。'// &
                '已自动强制禁用地球高阶场！'
                config%use_earth_nspheric = .false.
            end if
        end if
        
        ! 如果未来有别的依赖，比如用了大气阻力就必须指定迎风面积，也可以写在这里
    end subroutine resolve_config_dependencies
    
    subroutine save_default_config(config_file)
        character(len=*), intent(in) :: config_file
        integer :: unit, ios
        
        open(newunit=unit, file=config_file, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(*, *) '错误: 无法创建配置文件: ', trim(config_file)
            return
        end if
        
        write(unit, '(A)') '# POD Fortran 配置文件'
        write(unit, '(A)') '# 生成时间: ' // get_current_time()
        write(unit, '(A)') ''
        
        write(unit, '(A)') '# 轨道传播参数'
        write(unit, '(A)') 'LU = 384400.0'
        write(unit, '(A)') 'TU = 375190.26189464843'
        write(unit, '(A)') 'MU = 6.045626292270688E+24'
        write(unit, '(A)') 'VU = 1.024546847401753'
        write(unit, '(A)') 'AccU = 2.730739444643263E-06'
        write(unit, '(A)') 'propagation_step = 60.0'
        write(unit, '(A)') 'propagation_abs_tol = 1.0e-12'
        write(unit, '(A)') 'propagation_rel_tol = 1.0e-12'
        write(unit, '(A)') 'max_propagation_steps = 10000'
        write(unit, '(A)') 'rk4_step_size = 60.0'
        write(unit, '(A)') 'rkf45_abs_tol = 1.0e-12'
        write(unit, '(A)') 'rkf45_rel_tol = 1.0e-12'
        write(unit, '(A)') 'rkf45_min_step = 1.0e-6'
        write(unit, '(A)') 'rkf45_max_step = 3600.0'
        write(unit, '(A)') 'rkf78_abs_tol = 1.0e-12'
        write(unit, '(A)') 'rkf78_rel_tol = 1.0e-12'
        write(unit, '(A)') 'rkf78_min_step = 1.0e-6'
        write(unit, '(A)') 'rkf78_max_step = 3600.0'
        write(unit, '(A)') ''
        
        write(unit, '(A)') '# 轨道改进参数'
        write(unit, '(A)') 'convergence_tolerance = 1.0e-8'
        write(unit, '(A)') 'max_iterations = 100'
        write(unit, '(A)') 'damping_factor = 0.1'
        write(unit, '(A)') 'use_adaptive_damping = true'
        write(unit, '(A)') ''
        
        write(unit, '(A)') '# 坐标系统参数'
        write(unit, '(A)') 'reference_frame = J2000'
        write(unit, '(A)') 'central_body = EARTH'
        write(unit, '(A)') 'time_system = UTC'
        write(unit, '(A)') 'use_precession = true'
        write(unit, '(A)') 'use_nutation = true'
        write(unit, '(A)') ''
        
        write(unit, '(A)') '# 输出参数'
        write(unit, '(A)') 'verbose_output = true'
        write(unit, '(A)') 'save_results = true'
        write(unit, '(A)') 'output_format = CSV'
        write(unit, '(A)') 'date_format = YYYY-MM-DD HH:MM:SS'
        write(unit, '(A)') 'number_format = F20.12'
        write(unit, '(A)') 'scientific_notation = false'
        write(unit, '(A)') 'output_precision = 12'
        write(unit, '(A)') ''
        
        write(unit, '(A)') '# 数值计算参数'
        write(unit, '(A)') 'gravitational_constant = 398600.4418'
        write(unit, '(A)') 'numerical_tolerance = 1.0e-15'
        write(unit, '(A)') 'max_matrix_size = 1000'
        write(unit, '(A)') ''
        
        write(unit, '(A)') '# 性能参数'
        write(unit, '(A)') 'max_threads = 1'
        write(unit, '(A)') 'parallel_computation = false'
        write(unit, '(A)') 'memory_limit = 1.0e9'
        write(unit, '(A)') 'use_gpu = false'
        write(unit, '(A)') ''
        
        write(unit, '(A)') '# 力模型参数 (地月空间 Cislunar 特化)'
        write(unit, '(A)') 'use_earth_nspheric = true'
        write(unit, '(A)') 'use_moon_nspheric = true'
        write(unit, '(A)') 'use_third_body = true'
        write(unit, '(A)') 'use_srp = true'
        write(unit, '(A)') 'use_drag = false'
        write(unit, '(A)') 'use_relativity = false'
        write(unit, '(A)') 'earth_degree = 10'
        write(unit, '(A)') 'moon_degree = 10'
        write(unit, '(A)') 'active_planets = 3, 10, 11'
        write(unit, '(A)') ''
        
        write(unit, '(A)') '# 测量模型参数'
        write(unit, '(A)') 'measurement_noise_std = 1.0e-3'
        write(unit, '(A)') 'use_measurement_weights = true'
        write(unit, '(A)') 'outlier_threshold = 3.0'
        write(unit, '(A)') ''
        
        write(unit, '(A)') '# 文件I/O参数'
        write(unit, '(A)') 'data_directory = ./data/'
        write(unit, '(A)') 'results_directory = ./results/'
        write(unit, '(A)') 'auto_backup = true'
        write(unit, '(A)') 'backup_interval = 3600'
        write(unit, '(A)') ''
        
        write(unit, '(A)') '# 调试参数'
        write(unit, '(A)') 'debug_mode = false'
        write(unit, '(A)') 'debug_level = 1'
        write(unit, '(A)') 'save_intermediate_results = false'
        write(unit, '(A)') 'debug_directory = ./debug/'
        
        close(unit)
        
        write(*, *) '默认配置文件已创建: ', trim(config_file)
    end subroutine save_default_config
    
    function get_current_time() result(time_str)
        character(len=32) :: time_str
        character(len=8) :: date
        character(len=10) :: time
        character(len=5) :: zone
        
        call date_and_time(date, time, zone)
        time_str = date(1:4) // '-' // date(5:6) // '-' // date(7:8) // ' ' // &
                   time(1:2) // ':' // time(3:4) // ':' // time(5:6)
    end function get_current_time
    
    subroutine parse_key_value(line, key, value)
        character(len=*), intent(in) :: line
        character(len=*), intent(out) :: key, value
        integer :: pos
        
        pos = index(line, '=')
        if (pos > 0) then
            key = adjustl(line(1:pos-1))
            value = adjustl(line(pos+1:))
        else
            key = ''
            value = ''
        end if
    end subroutine parse_key_value
    
    subroutine set_config_value(key, value)
        character(len=*), intent(in) :: key, value
        integer :: ios
        
        select case (trim(key))
            ! 轨道传播参数
            case ('LU')
                read(value, *, iostat=ios) config%LU
            case ('TU')
                read(value, *, iostat=ios) config%TU
            case ('MU')
                read(value, *, iostat=ios) config%MU
            case ('VU')
                read(value, *, iostat=ios) config%VU
            case ('AccU')
                read(value, *, iostat=ios) config%AccU
            case ('propagation_step')
                read(value, *, iostat=ios) config%propagation_step
            case ('propagation_abs_tol')
                read(value, *, iostat=ios) config%propagation_abs_tol
            case ('propagation_rel_tol')
                read(value, *, iostat=ios) config%propagation_rel_tol
            case ('max_propagation_steps')
                read(value, *, iostat=ios) config%max_propagation_steps
            case ('rk4_step_size')
                read(value, *, iostat=ios) config%rk4_step_size
            case ('rkf45_abs_tol')
                read(value, *, iostat=ios) config%rkf45_abs_tol
            case ('rkf45_rel_tol')
                read(value, *, iostat=ios) config%rkf45_rel_tol
            case ('rkf45_min_step')
                read(value, *, iostat=ios) config%rkf45_min_step
            case ('rkf45_max_step')
                read(value, *, iostat=ios) config%rkf45_max_step
            case ('rkf78_abs_tol')
                read(value, *, iostat=ios) config%rkf78_abs_tol
            case ('rkf78_rel_tol')
                read(value, *, iostat=ios) config%rkf78_rel_tol
            case ('rkf78_min_step')
                read(value, *, iostat=ios) config%rkf78_min_step
            case ('rkf78_max_step')
                read(value, *, iostat=ios) config%rkf78_max_step
            
            ! 轨道改进参数
            case ('convergence_tolerance')
                read(value, *, iostat=ios) config%convergence_tolerance
            case ('max_iterations')
                read(value, *, iostat=ios) config%max_iterations
            case ('damping_factor')
                read(value, *, iostat=ios) config%damping_factor
            case ('use_adaptive_damping')
                config%use_adaptive_damping = (trim(value) == 'true')
            
            ! 坐标系统参数
            case ('reference_frame')
                config%reference_frame = trim(value)
            case ('central_body')
                config%central_body = trim(value)
            case ('time_system')
                config%time_system = trim(value)
            case ('use_precession')
                config%use_precession = (trim(value) == 'true')
            case ('use_nutation')
                config%use_nutation = (trim(value) == 'true')
            
            ! 输出参数
            case ('verbose_output')
                config%verbose_output = (trim(value) == 'true')
            case ('save_results')
                config%save_results = (trim(value) == 'true')
            case ('output_format')
                config%output_format = trim(value)
            case ('date_format')
                config%date_format = trim(value)
            case ('number_format')
                config%number_format = trim(value)
            case ('scientific_notation')
                config%scientific_notation = (trim(value) == 'true')
            case ('output_precision')
                read(value, *, iostat=ios) config%output_precision
            
            ! 数值计算参数
            case ('gravitational_constant')
                read(value, *, iostat=ios) config%gravitational_constant
            case ('numerical_tolerance')
                read(value, *, iostat=ios) config%numerical_tolerance
            case ('max_matrix_size')
                read(value, *, iostat=ios) config%max_matrix_size
            
            ! 性能参数
            case ('max_threads')
                read(value, *, iostat=ios) config%max_threads
            case ('parallel_computation')
                config%parallel_computation = (trim(value) == 'true')
            case ('memory_limit')
                read(value, *, iostat=ios) config%memory_limit
            case ('use_gpu')
                config%use_gpu = (trim(value) == 'true')
            
            ! 力模型参数
            ! =====================================
            ! 力模型参数 (地月空间 Cislunar 特化)
            ! =====================================
            case ('use_earth_nspheric')
                config%use_earth_nspheric = (trim(value) == 'true')
            case ('use_moon_nspheric')
                config%use_moon_nspheric = (trim(value) == 'true')
            case ('use_third_body')
                config%use_third_body = (trim(value) == 'true')
            case ('use_srp')
                config%use_srp = (trim(value) == 'true')
            case ('use_drag')
                config%use_drag = (trim(value) == 'true')
            case ('use_relativity')
                config%use_relativity = (trim(value) == 'true')
            case ('earth_degree')
                read(value, *, iostat=ios) config%earth_degree
            case ('moon_degree')
                read(value, *, iostat=ios) config%moon_degree
            case ('active_planets')
                ! 读取天体ID列表 (支持逗号或空格分隔，如 "3, 10, 11")
                call parse_active_planets(value)
            
            ! 测量模型参数
            case ('measurement_noise_std')
                read(value, *, iostat=ios) config%measurement_noise_std
            case ('use_measurement_weights')
                config%use_measurement_weights = (trim(value) == 'true')
            case ('outlier_threshold')
                read(value, *, iostat=ios) config%outlier_threshold
            
            ! 文件I/O参数
            case ('data_directory')
                config%data_directory = trim(value)
            case ('results_directory')
                config%results_directory = trim(value)
            case ('auto_backup')
                config%auto_backup = (trim(value) == 'true')
            case ('backup_interval')
                read(value, *, iostat=ios) config%backup_interval
            
            ! 调试参数
            case ('debug_mode')
                config%debug_mode = (trim(value) == 'true')
            case ('debug_level')
                read(value, *, iostat=ios) config%debug_level
            case ('save_intermediate_results')
                config%save_intermediate_results = (trim(value) == 'true')
            case ('debug_directory')
                config%debug_directory = trim(value)
            
            case default
                write(*, *) '警告: 未知配置参数: ', trim(key)
        end select
        
        if (ios /= 0) then
            write(*, *) '错误: 配置参数解析错误: ', trim(key), ' = ', trim(value)
        end if

        call resolve_config_dependencies()  ! 每次设置参数后检查依赖关系
    end subroutine set_config_value

    subroutine parse_active_planets(val_str)
        character(len=*), intent(in) :: val_str
        integer :: id, i, ios
        character(len=MAX_STRING_LEN) :: temp_str
        
        ! 先清空所有的激活状态
        config%use_planet = .false.
        temp_str = val_str
        
        ! 将所有逗号替换为空格，方便 list-directed read
        do i = 1, len_trim(temp_str)
            if (temp_str(i:i) == ',') temp_str(i:i) = ' '
        end do
        
        ! 循环读取 ID
        do
            if (len_trim(temp_str) == 0) exit
            read(temp_str, *, iostat=ios) id
            if (ios /= 0) exit
            
            if (id >= 1 .and. id <= 11) then
                config%use_planet(id) = .true.
            end if
            
            ! 截断已经读过的数字，准备读下一个
            i = scan(trim(temp_str), '0123456789') ! 找到第一个数字
            if (i > 0) then
                i = i + scan(trim(temp_str(i:)), ' ') - 1 ! 找到数字后的第一个空格
                if (i > 0) then
                    temp_str = adjustl(temp_str(i:))
                else
                    exit
                end if
            else
                exit
            end if
        end do
    end subroutine parse_active_planets
    
    logical function file_exists(filename)
        character(len=*), intent(in) :: filename
        inquire(file=filename, exist=file_exists)
    end function file_exists
    
    subroutine print_config()
        integer :: i

        write(*, *) '=== 当前配置 ==='
        write(*, *) ''
        write(*, *) '轨道传播参数:'
        write(*, *) '  长度单位 (LU): ', config%LU, ' km'
        write(*, *) '  时间单位 (TU): ', config%TU, ' 秒'
        write(*, *) '  质量单位 (MU): ', config%MU, ' kg'
        write(*, *) '  速度单位 (VU): ', config%VU, ' km/s'
        write(*, *) '  加速度单位 (AccU): ', config%AccU, ' km/s²'
        write(*, *) '  传播步长: ', config%propagation_step, ' 秒'
        write(*, *) '  传播绝对容差: ', config%propagation_abs_tol
        write(*, *) '  传播相对容差: ', config%propagation_rel_tol
        write(*, *) '  最大传播步数: ', config%max_propagation_steps
        write(*, *) '  RK4步长: ', config%rk4_step_size, ' 秒'
        write(*, *) '  RKF45绝对容差: ', config%rkf45_abs_tol
        write(*, *) '  RKF45相对容差: ', config%rkf45_rel_tol
        write(*, *) '  RKF45最小步长: ', config%rkf45_min_step, ' 秒'
        write(*, *) '  RKF45最大步长: ', config%rkf45_max_step, ' 秒'
        write(*, *) '  RKF78绝对容差: ', config%rkf78_abs_tol
        write(*, *) '  RKF78相对容差: ', config%rkf78_rel_tol
        write(*, *) '  RKF78最小步长: ', config%rkf78_min_step, ' 秒'
        write(*, *) '  RKF78最大步长: ', config%rkf78_max_step, ' 秒'
        write(*, *) ''
        write(*, *) '轨道改进参数:'
        write(*, *) '  收敛容差: ', config%convergence_tolerance
        write(*, *) '  最大迭代次数: ', config%max_iterations
        write(*, *) '  阻尼因子: ', config%damping_factor
        write(*, *) '  自适应阻尼: ', config%use_adaptive_damping
        write(*, *) ''
        write(*, *) '坐标系统参数:'
        write(*, *) '  参考坐标系: ', trim(config%reference_frame)
        write(*, *) '  中心天体: ', trim(config%central_body)
        write(*, *) '  时间系统: ', trim(config%time_system)
        write(*, *) '  岁差修正: ', config%use_precession
        write(*, *) '  章动修正: ', config%use_nutation
        write(*, *) ''
        write(*, *) '输出参数:'
        write(*, *) '  详细输出: ', config%verbose_output
        write(*, *) '  保存结果: ', config%save_results
        write(*, *) '  输出格式: ', trim(config%output_format)
        write(*, *) '  日期格式: ', trim(config%date_format)
        write(*, *) '  数字格式: ', trim(config%number_format)
        write(*, *) '  科学计数法: ', config%scientific_notation
        write(*, *) '  输出精度: ', config%output_precision, ' 位'
        write(*, *) ''
        write(*, *) '数值计算参数:'
        write(*, *) '  引力常数: ', config%gravitational_constant, ' km³/s²'
        write(*, *) '  数值容差: ', config%numerical_tolerance
        write(*, *) '  最大矩阵大小: ', config%max_matrix_size
        write(*, *) ''
        write(*, *) '性能参数:'
        write(*, *) '  最大线程数: ', config%max_threads
        write(*, *) '  并行计算: ', config%parallel_computation
        write(*, *) '  内存限制: ', config%memory_limit, ' bytes'
        write(*, *) '  GPU加速: ', config%use_gpu
        write(*, *) ''
        write(*, *) '力模型参数 :'
        ! 地球高阶场状态输出
        if (.not. config%use_planet(3)) then
            write(*, *) '  地球高阶引力场: [未激活地球本体, 强制禁用]'
        else if (config%use_earth_nspheric) then
            write(*, *) '  地球高阶引力场: 启用 (', config%earth_degree, '阶)'
        else
            write(*, *) '  地球高阶引力场: 禁用'
        end if
        ! 月球高阶场状态输出
        if (.not. config%use_third_body .or. .not. config%use_planet(10)) then
            write(*, *) '  月球高阶引力场: [未激活月球第三体, 强制禁用]'
        else if (config%use_moon_nspheric) then
            write(*, *) '  月球高阶引力场: 启用 (', config%moon_degree, '阶)'
        else
            write(*, *) '  月球高阶引力场: 禁用'
        end if
        write(*, *) '  第三体摄动: ', config%use_third_body
        write(*, *) '  太阳辐射压: ', config%use_srp
        write(*, *) '  大气阻力: ', config%use_drag
        write(*, *) '  相对论效应: ', config%use_relativity
        write(*, "(A)", advance='no') '  已激活的引力网络节点 (ID): '
        do i = 1, 11
            if (config%use_planet(i)) write(*, "(I3)", advance='no') i
        end do
        write(*, *) '' ! 换行
        write(*, *) ''
        write(*, *) '测量模型参数:'
        write(*, *) '  测量噪声标准差: ', config%measurement_noise_std
        write(*, *) '  使用测量权重: ', config%use_measurement_weights
        write(*, *) '  异常值阈值: ', config%outlier_threshold
        write(*, *) ''
        write(*, *) '文件I/O参数:'
        write(*, *) '  数据目录: ', trim(config%data_directory)
        write(*, *) '  结果目录: ', trim(config%results_directory)
        write(*, *) '  自动备份: ', config%auto_backup
        write(*, *) '  备份间隔: ', config%backup_interval, ' 秒'
        write(*, *) ''
        write(*, *) '调试参数:'
        write(*, *) '  调试模式: ', config%debug_mode
        write(*, *) '  调试级别: ', config%debug_level
        write(*, *) '  保存中间结果: ', config%save_intermediate_results
        write(*, *) '  调试目录: ', trim(config%debug_directory)
        write(*, *) ''
    end subroutine print_config
    
    ! 配置验证函数
    logical function validate_config()
        validate_config = .true.
        
        ! 验证轨道传播参数
        if (config%propagation_step <= 0.0_DP) then
            write(*, *) '错误: 传播步长必须大于0'
            validate_config = .false.
        end if
        
        if (config%propagation_abs_tol <= 0.0_DP) then
            write(*, *) '错误: 传播绝对容差必须大于0'
            validate_config = .false.
        end if

        if (config%propagation_rel_tol <= 0.0_DP) then
            write(*, *) '错误: 传播相对容差必须大于0'
            validate_config = .false.
        end if

        if (config%max_propagation_steps <= 0) then
            write(*, *) '错误: 最大传播步数必须大于0'
            validate_config = .false.
        end if
        
        ! 验证积分器参数
        if (config%rk4_step_size <= 0.0_DP) then
            write(*, *) '错误: RK4步长必须大于0'
            validate_config = .false.
        end if
        
        if (config%rkf45_abs_tol <= 0.0_DP) then
                        write(*, *) '错误: RKF45绝对容差必须大于0'
            validate_config = .false.
        end if

        if (config%rkf45_rel_tol <= 0.0_DP) then
                        write(*, *) '错误: RKF45相对容差必须大于0'
            validate_config = .false.
        end if
        
        if (config%rkf45_min_step >= config%rkf45_max_step) then
            write(*, *) '错误: RKF45最小步长必须小于最大步长'
            validate_config = .false.
        end if
        
        if (config%rkf78_abs_tol <= 0.0_DP) then
            write(*, *) '错误: RKF78绝对容差必须大于0'
            validate_config = .false.
        end if

        if (config%rkf78_rel_tol <= 0.0_DP) then
            write(*, *) '错误: RKF78相对容差必须大于0'
            validate_config = .false.
        end if
        
        if (config%rkf78_min_step >= config%rkf78_max_step) then
            write(*, *) '错误: RKF78最小步长必须小于最大步长'
            validate_config = .false.
        end if
        
        ! 验证轨道改进参数
        if (config%convergence_tolerance <= 0.0_DP) then
            write(*, *) '错误: 收敛容差必须大于0'
            validate_config = .false.
        end if
        
        if (config%max_iterations <= 0) then
            write(*, *) '错误: 最大迭代次数必须大于0'
            validate_config = .false.
        end if
        
        if (config%damping_factor < 0.0_DP .or. config%damping_factor > 1.0_DP) then
            write(*, *) '错误: 阻尼因子必须在0到1之间'
            validate_config = .false.
        end if
        
        ! 验证数值计算参数
        if (config%gravitational_constant <= 0.0_DP) then
            write(*, *) '错误: 引力常数必须大于0'
            validate_config = .false.
        end if
        
        if (config%numerical_tolerance <= 0.0_DP) then
            write(*, *) '错误: 数值容差必须大于0'
            validate_config = .false.
        end if
        
        if (config%max_matrix_size <= 0) then
            write(*, *) '错误: 最大矩阵大小必须大于0'
            validate_config = .false.
        end if
        
        ! 验证性能参数
        if (config%max_threads <= 0) then
            write(*, *) '错误: 最大线程数必须大于0'
            validate_config = .false.
        end if
        
        if (config%memory_limit <= 0.0_DP) then
            write(*, *) '错误: 内存限制必须大于0'
            validate_config = .false.
        end if
        
        ! ! 验证力模型参数
        if (config%use_earth_nspheric .and. (config%earth_degree < 2 .or. config%earth_degree > 300)) then
            write(*, *) '错误: 地球重力场阶数设置不合理 (通常在 2-300 之间)'
            validate_config = .false.
        end if
        
        if (config%use_moon_nspheric .and. (config%moon_degree < 2 .or. config%moon_degree > 300)) then
            write(*, *) '错误: 月球重力场阶数设置不合理 (通常在 2-300 之间)'
            validate_config = .false.
        end if
        
        ! 校验多体引力网
        if (.not. any(config%use_planet)) then
            write(*, *) '警告: 未激活任何天体的引力！系统处于无重力环境。'
            ! 这里给警告即可，因为有时候单纯测光压的时候可能会故意关掉引力
        end if
        
        ! 验证测量模型参数
        if (config%measurement_noise_std < 0.0_DP) then
            write(*, *) '错误: 测量噪声标准差不能为负'
            validate_config = .false.
        end if
        
        if (config%outlier_threshold <= 0.0_DP) then
            write(*, *) '错误: 异常值阈值必须大于0'
            validate_config = .false.
        end if
        
        ! 验证输出参数
        if (config%output_precision < 1 .or. config%output_precision > 20) then
            write(*, *) '错误: 输出精度必须在1到20之间'
            validate_config = .false.
        end if
        
        ! 验证调试参数
        if (config%debug_level < 0 .or. config%debug_level > 5) then
            write(*, *) '错误: 调试级别必须在0到5之间'
            validate_config = .false.
        end if
        
        if (config%backup_interval <= 0) then
            write(*, *) '错误: 备份间隔必须大于0'
            validate_config = .false.
        end if
        
        if (validate_config) then
            write(*, *) '配置验证通过'
        else
            write(*, *) '配置验证失败'
        end if
    end function validate_config
    
    ! 配置范围检查
    subroutine check_config_ranges()
        ! 检查并修正超出合理范围的参数
        
        ! 传播步长范围检查
        if (config%propagation_step < 1.0_DP) then
            config%propagation_step = 1.0_DP
            write(*, *) '警告: 传播步长已调整为最小值: 1.0秒'
        else if (config%propagation_step > 3600.0_DP) then
            config%propagation_step = 3600.0_DP
            write(*, *) '警告: 传播步长已调整为最大值: 3600.0秒'
        end if
        
        ! 容差范围检查
        if (config%propagation_abs_tol < 1.0e-15_DP) then
            config%propagation_abs_tol = 1.0e-15_DP
            write(*, *) '警告: 传播绝对容差已调整为最小值: 1.0e-15'
        else if (config%propagation_rel_tol > 1.0e-6_DP) then
            config%propagation_rel_tol = 1.0e-6_DP
            write(*, *) '警告: 传播相对容差已调整为最大值: 1.0e-6'
        end if
        
        ! 最大传播步数范围检查
        if (config%max_propagation_steps < 100) then
            config%max_propagation_steps = 100
            write(*, *) '警告: 最大传播步数已调整为最小值: 100'
        else if (config%max_propagation_steps > 1000000) then
            config%max_propagation_steps = 1000000
            write(*, *) '警告: 最大传播步数已调整为最大值: 1000000'
        end if
        
        ! 最大迭代次数范围检查
        if (config%max_iterations < 10) then
            config%max_iterations = 10
            write(*, *) '警告: 最大迭代次数已调整为最小值: 10'
        else if (config%max_iterations > 10000) then
            config%max_iterations = 10000
            write(*, *) '警告: 最大迭代次数已调整为最大值: 10000'
        end if
        
        ! 线程数范围检查
        if (config%max_threads < 1) then
            config%max_threads = 1
            write(*, *) '警告: 最大线程数已调整为最小值: 1'
        else if (config%max_threads > 64) then
            config%max_threads = 64
            write(*, *) '警告: 最大线程数已调整为最大值: 64'
        end if
        
        write(*, *) '配置范围检查完成'
    end subroutine check_config_ranges
    
    ! 获取配置值函数
    function get_config_value(key) result(value)
        character(len=*), intent(in) :: key
        character(len=MAX_STRING_LEN) :: value
        
        select case (trim(key))
            case ('propagation_step')
                write(value, '(F20.12)') config%propagation_step
            case ('propagation_abs_tol')
                write(value, '(E20.12)') config%propagation_abs_tol
            case ('propagation_rel_tol')
                write(value, '(E20.12)') config%propagation_rel_tol
            case ('max_propagation_steps')
                write(value, '(I10)') config%max_propagation_steps
            case ('convergence_tolerance')
                write(value, '(E20.12)') config%convergence_tolerance
            case ('max_iterations')
                write(value, '(I10)') config%max_iterations
            case ('reference_frame')
                value = trim(config%reference_frame)
            case ('central_body')
                value = trim(config%central_body)
            case ('verbose_output')
                value = merge('true ', 'false', config%verbose_output)
            case ('save_results')
                value = merge('true ', 'false', config%save_results)
            case ('output_format')
                value = trim(config%output_format)
            case ('gravitational_constant')
                write(value, '(F20.12)') config%gravitational_constant
            case default
                value = 'unknown'
                write(*, *) '警告: 未知配置参数: ', trim(key)
        end select
    end function get_config_value

end module pod_config
