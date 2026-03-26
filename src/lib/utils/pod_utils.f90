module pod_utils
    use pod_global, only: DP, MAX_STRING_LEN
    
    implicit none
    
contains

    subroutine print_banner()
        write(*, *) '=========================================='
        write(*, *) '    CAT Fortran - 地月空间轨道动力学工具'
        write(*, *) '    Cislunar-space Astrodynamics Tools'
        write(*, *) '=========================================='
        write(*, *)
    end subroutine print_banner
    
    subroutine print_help()
        write(*, *) 'CAT Fortran - 地月空间轨道动力学工具'
        write(*, *)
        write(*, *) '用法: pod_fortran [选项]'
        write(*, *)
        write(*, *) '选项:'
        write(*, *) '  -c, --config FILE    指定配置文件路径'
        write(*, *) '  -h, --help           显示此帮助信息'
        write(*, *)
        write(*, *) '功能模块:'
        write(*, *) '  OP - 轨道传播 (Orbit Propagation)'
        write(*, *) '  OM - 轨道改进 (Orbit Improvement)'
        write(*, *) '  GS - 观测几何仿真 (Geometry Simulation)'
        write(*, *) '  VC - 可见性计算 (Visibility Calculation)'
        write(*, *) '  OS - 观测仿真 (Observation Simulation)'
        write(*, *) '  IOD - 初轨确定 (Initial Orbit Determination)'
        write(*, *) '  ID - 目标识别 (Identification)'
        write(*, *) '  LK - 轨道关联 (Linkage)'
        write(*, *) '  MM - 机动建模 (Maneuver Modeling)'
        write(*, *) '  MDA - 任务设计与分析 (Mission Design & Analysis)'
        write(*, *)
    end subroutine print_help
    
    integer function get_user_choice(prompt, min_val, max_val)
        character(len=*), intent(in) :: prompt
        integer, intent(in) :: min_val, max_val
        integer :: choice, ios
        
        do
            write(*, '(A)', advance='no') prompt
            read(*, *, iostat=ios) choice
            
            if (ios == 0 .and. choice >= min_val .and. choice <= max_val) then
                get_user_choice = choice
                return
            else
                write(*, *) '无效输入，请输入 ', min_val, ' 到 ', max_val, ' 之间的数字'
            end if
        end do
    end function get_user_choice
    
    real(DP) function get_user_real(prompt, min_val, max_val)
        character(len=*), intent(in) :: prompt
        real(DP), intent(in) :: min_val, max_val
        real(DP) :: value
        integer :: ios
        
        do
            write(*, '(A)', advance='no') prompt
            read(*, *, iostat=ios) value
            
            if (ios == 0 .and. value >= min_val .and. value <= max_val) then
                get_user_real = value
                return
            else
                write(*, *) '无效输入，请输入 ', min_val, ' 到 ', max_val, ' 之间的数值'
            end if
        end do
    end function get_user_real
    
    subroutine get_user_string(prompt, result)
        character(len=*), intent(in) :: prompt
        character(len=*), intent(out) :: result
        
        write(*, '(A)', advance='no') prompt
        read(*, '(A)') result
        result = adjustl(result)
    end subroutine get_user_string
    
    subroutine print_separator(title)
        character(len=*), intent(in), optional :: title
        integer :: i, len_title
        
        if (present(title)) then
            len_title = len_trim(title)
            write(*, *) repeat('=', 20) // ' ' // trim(title) // ' ' // repeat('=', 20)
        else
            write(*, *) repeat('=', 60)
        end if
    end subroutine print_separator
    
    subroutine print_matrix(matrix, title, format_str)
        real(DP), dimension(:,:), intent(in) :: matrix
        character(len=*), intent(in), optional :: title
        character(len=*), intent(in), optional :: format_str
        integer :: i, j, rows, cols
        character(len=32) :: fmt
        
        rows = size(matrix, 1)
        cols = size(matrix, 2)
        
        if (present(title)) then
            write(*, *) trim(title)
        end if
        
        if (present(format_str)) then
            fmt = format_str
        else
            fmt = '(F12.6)'
        end if
        
        do i = 1, rows
            write(*, '(A,I3,A)', advance='no') 'Row ', i, ': '
            do j = 1, cols
                write(*, fmt, advance='no') matrix(i, j)
                if (j < cols) write(*, '(A)', advance='no') '  '
            end do
            write(*, *)
        end do
    end subroutine print_matrix
    
    subroutine print_vector(vector, title, format_str)
        real(DP), dimension(:), intent(in) :: vector
        character(len=*), intent(in), optional :: title
        character(len=*), intent(in), optional :: format_str
        integer :: i, n
        character(len=32) :: fmt
        
        n = size(vector)
        
        if (present(title)) then
            write(*, *) trim(title)
        end if
        
        if (present(format_str)) then
            fmt = format_str
        else
            fmt = '(F12.6)'
        end if
        
        write(*, '(A)', advance='no') 'Vector: '
        do i = 1, n
            write(*, fmt, advance='no') vector(i)
            if (i < n) write(*, '(A)', advance='no') '  '
        end do
        write(*, *)
    end subroutine print_vector
    
    subroutine save_matrix_to_file(matrix, filename, title)
        real(DP), dimension(:,:), intent(in) :: matrix
        character(len=*), intent(in) :: filename
        character(len=*), intent(in), optional :: title
        integer :: unit, i, j, rows, cols
        
        rows = size(matrix, 1)
        cols = size(matrix, 2)
        
        open(newunit=unit, file=filename, status='replace', action='write')
        
        if (present(title)) then
            write(unit, '(A)') '# ' // trim(title)
        end if
        
        write(unit, '(A,I6,A,I6)') '# Matrix size: ', rows, ' x ', cols
        
        do i = 1, rows
            do j = 1, cols
                write(unit, '(F20.12)', advance='no') matrix(i, j)
                if (j < cols) write(unit, '(A)', advance='no') ','
            end do
            write(unit, *)
        end do
        
        close(unit)
        write(*, *) '矩阵已保存到文件: ', trim(filename)
    end subroutine save_matrix_to_file
    
    subroutine save_vector_to_file(vector, filename, title)
        real(DP), dimension(:), intent(in) :: vector
        character(len=*), intent(in) :: filename
        character(len=*), intent(in), optional :: title
        integer :: unit, i, n
        
        n = size(vector)
        
        open(newunit=unit, file=filename, status='replace', action='write')
        
        if (present(title)) then
            write(unit, '(A)') '# ' // trim(title)
        end if
        
        write(unit, '(A,I6)') '# Vector size: ', n
        
        do i = 1, n
            write(unit, '(F20.12)') vector(i)
        end do
        
        close(unit)
        write(*, *) '向量已保存到文件: ', trim(filename)
    end subroutine save_vector_to_file
    
    logical function confirm_action(prompt)
        character(len=*), intent(in) :: prompt
        character(len=10) :: response
        
        write(*, '(A)', advance='no') prompt // ' (y/n): '
        read(*, '(A)') response
        
        confirm_action = (response(1:1) == 'y' .or. response(1:1) == 'Y')
    end function confirm_action

    subroutine save_uq_results(final_state, base_filename)
        use pod_global, only: DP, MAX_STRING_LEN, output_directory
        use pod_uq_state_module, only: uq_state_type
        implicit none
        
        ! 传入参数
        type(uq_state_type), intent(in) :: final_state
        character(len=*), intent(in)    :: base_filename
        
        ! 局部变量
        integer :: file_unit, i, dim, n_particles
        character(len=MAX_STRING_LEN) :: filepath_particles, filepath_stats
        
        ! 获取维度信息
        dim = size(final_state%samples, 1)
        n_particles = size(final_state%samples, 2)
        
        ! ==========================================================
        ! 任务 1：保存所有粒子 (Samples) 到 CSV
        ! ==========================================================
        ! 拼接完整文件路径 (例如: ./output/up_result_particles.csv)
        filepath_particles = trim(output_directory) // trim(base_filename) // '_particles.csv'
        
        ! 打开文件
        open(newunit=file_unit, file=trim(filepath_particles), status='replace', action='write')
        
        ! 写表头 (假设是 6 维)
        write(file_unit, '(A)') 'x,y,z,vx,vy,vz'
        
        ! 遍历所有粒子写入文件
        do i = 1, n_particles
            ! 使用科学计数法，保留 14 位小数，逗号分隔
            write(file_unit, '(ES22.14, 5(",", ES22.14))') final_state%samples(:, i)
        end do
        
        close(file_unit)
        
        ! ==========================================================
        ! 任务 2：保存统计矩 (Mean & Covariance) 到同一个 CSV
        ! ==========================================================
        filepath_stats = trim(output_directory) // trim(base_filename) // '_stats.csv'
        open(newunit=file_unit, file=trim(filepath_stats), status='replace', action='write')
        
        ! 写入均值
        write(file_unit, '(A)') '# Mean'
        write(file_unit, '(ES22.14, 5(",", ES22.14))') final_state%mean(:)
        
        ! 写入协方差矩阵
        write(file_unit, '(A)') '# Covariance Matrix'
        do i = 1, dim
            write(file_unit, '(ES22.14, 5(",", ES22.14))') final_state%cov(i, :)
        end do
        
        close(file_unit)
        
        write(*,*) '[INFO] 不确定性传播结果已成功保存至:'
        write(*,*) '       ', trim(filepath_particles)
        write(*,*) '       ', trim(filepath_stats)

    end subroutine save_uq_results
    
    subroutine pause_execution()
        character(len=1) :: dummy
        
        write(*, *) '按回车键继续...'
        read(*, '(A)') dummy
    end subroutine pause_execution
    
    real(DP) function norm_vector(vector)
        real(DP), dimension(:), intent(in) :: vector
        norm_vector = sqrt(sum(vector**2))
    end function norm_vector
    
    real(DP) function norm_matrix(matrix)
        real(DP), dimension(:,:), intent(in) :: matrix
        norm_matrix = sqrt(sum(matrix**2))
    end function norm_matrix
    

end module pod_utils
