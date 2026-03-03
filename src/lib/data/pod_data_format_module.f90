module pod_data_format_module
    use pod_global, only: DP, MAX_STRING_LEN
    implicit none
    
    private
    public :: read_tle_file, write_tle_file, read_observation_file, write_observation_file, &
              read_state_vector_file, write_state_vector_file, format_time_string, &
              parse_csv_line, format_csv_line
    
contains

    ! 读取TLE文件
    subroutine read_tle_file(filename, tle_data, n_lines, error_code)
        character(len=*), intent(in) :: filename
        character(len=*), dimension(:), intent(out) :: tle_data
        integer, intent(out) :: n_lines
        integer, intent(out) :: error_code
        
        integer :: unit, i, ios
        character(len=MAX_STRING_LEN) :: line
        
        error_code = 0
        n_lines = 0
        
        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            error_code = 1001  ! 文件打开错误
            return
        end if
        
        do i = 1, size(tle_data)
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            tle_data(i) = line
            n_lines = n_lines + 1
        end do
        
        close(unit)
        
    end subroutine read_tle_file
    
    ! 写入TLE文件
    subroutine write_tle_file(filename, tle_data, n_lines, error_code)
        character(len=*), intent(in) :: filename
        character(len=*), dimension(:), intent(in) :: tle_data
        integer, intent(in) :: n_lines
        integer, intent(out) :: error_code
        
        integer :: unit, i, ios
        
        error_code = 0
        
        open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            error_code = 1002  ! 文件创建错误
            return
        end if
        
        do i = 1, n_lines
            write(unit, '(A)', iostat=ios) trim(tle_data(i))
            if (ios /= 0) then
                error_code = 1003  ! 文件写入错误
                exit
            end if
        end do
        
        close(unit)
        
    end subroutine write_tle_file
    
    ! 读取观测数据文件
    subroutine read_observation_file(filename, observations, n_obs, error_code)
        character(len=*), intent(in) :: filename
        real(DP), dimension(:,:), intent(out) :: observations
        integer, intent(out) :: n_obs
        integer, intent(out) :: error_code
        
        integer :: unit, i, ios
        real(DP) :: time, range, az, el
        
        error_code = 0
        n_obs = 0
        
        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            error_code = 1001
            return
        end if
        
        do i = 1, size(observations, 1)
            read(unit, *, iostat=ios) time, range, az, el
            if (ios /= 0) exit
            observations(i, 1) = time
            observations(i, 2) = range
            observations(i, 3) = az
            observations(i, 4) = el
            n_obs = n_obs + 1
        end do
        
        close(unit)
        
    end subroutine read_observation_file
    
    ! 写入观测数据文件
    subroutine write_observation_file(filename, observations, n_obs, error_code)
        character(len=*), intent(in) :: filename
        real(DP), dimension(:,:), intent(in) :: observations
        integer, intent(in) :: n_obs
        integer, intent(out) :: error_code
        
        integer :: unit, i, ios
        
        error_code = 0
        
        open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            error_code = 1002
            return
        end if
        
        do i = 1, n_obs
            write(unit, '(F20.8,3F15.8)', iostat=ios) observations(i, 1), &
                  observations(i, 2), observations(i, 3), observations(i, 4)
            if (ios /= 0) then
                error_code = 1003
                exit
            end if
        end do
        
        close(unit)
        
    end subroutine write_observation_file
    
    ! 读取状态向量文件
    subroutine read_state_vector_file(filename, states, times, n_states, error_code)
        character(len=*), intent(in) :: filename
        real(DP), dimension(:,:), intent(out) :: states
        real(DP), dimension(:), intent(out) :: times
        integer, intent(out) :: n_states
        integer, intent(out) :: error_code
        
        integer :: unit, i, ios
        real(DP) :: time, x, y, z, vx, vy, vz
        
        error_code = 0
        n_states = 0
        
        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            error_code = 1001
            return
        end if
        
        do i = 1, size(states, 1)
            read(unit, *, iostat=ios) time, x, y, z, vx, vy, vz
            if (ios /= 0) exit
            times(i) = time
            states(i, 1) = x
            states(i, 2) = y
            states(i, 3) = z
            states(i, 4) = vx
            states(i, 5) = vy
            states(i, 6) = vz
            n_states = n_states + 1
        end do
        
        close(unit)
        
    end subroutine read_state_vector_file
    
    ! 写入状态向量文件
    subroutine write_state_vector_file(filename, states, times, n_states, error_code)
        character(len=*), intent(in) :: filename
        real(DP), dimension(:,:), intent(in) :: states
        real(DP), dimension(:), intent(in) :: times
        integer, intent(in) :: n_states
        integer, intent(out) :: error_code
        
        integer :: unit, i, ios
        
        error_code = 0
        
        open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            error_code = 1002
            return
        end if
        
        do i = 1, n_states
            write(unit, '(F20.8,6F15.8)', iostat=ios) times(i), &
                  states(i, 1), states(i, 2), states(i, 3), &
                  states(i, 4), states(i, 5), states(i, 6)
            if (ios /= 0) then
                error_code = 1003
                exit
            end if
        end do
        
        close(unit)
        
    end subroutine write_state_vector_file
    
    ! 格式化时间字符串
    subroutine format_time_string(time_et, format_type, time_string)
        real(DP), intent(in) :: time_et
        character(len=*), intent(in) :: format_type
        character(len=*), intent(out) :: time_string
        
        ! 这里应该调用SPICE的时间转换函数
        ! 暂时使用简单格式
        write(time_string, '(F20.8)') time_et
        
    end subroutine format_time_string
    
    ! 解析CSV行
    subroutine parse_csv_line(line, values, n_values, error_code)
        character(len=*), intent(in) :: line
        real(DP), dimension(:), intent(out) :: values
        integer, intent(out) :: n_values
        integer, intent(out) :: error_code
        
        integer :: i, pos, start_pos, ios
        character(len=MAX_STRING_LEN) :: token
        
        error_code = 0
        n_values = 0
        pos = 1
        
        do i = 1, size(values)
            start_pos = pos
            do while (pos <= len_trim(line) .and. line(pos:pos) /= ',')
                pos = pos + 1
            end do
            
            if (start_pos <= pos - 1) then
                token = line(start_pos:pos-1)
                read(token, *, iostat=ios) values(i)
                if (ios /= 0) then
                    error_code = 1004  ! 数据解析错误
                    return
                end if
                n_values = n_values + 1
            end if
            
            if (pos > len_trim(line)) exit
            pos = pos + 1  ! 跳过逗号
        end do
        
    end subroutine parse_csv_line
    
    ! 格式化CSV行
    subroutine format_csv_line(values, n_values, line)
        real(DP), dimension(:), intent(in) :: values
        integer, intent(in) :: n_values
        character(len=*), intent(out) :: line
        
        integer :: i
        character(len=MAX_STRING_LEN) :: temp_line
        
        line = ''
        do i = 1, n_values
            if (i > 1) then
                line = trim(line) // ','
            end if
            write(temp_line, '(F15.8)') values(i)
            line = trim(line) // trim(adjustl(temp_line))
        end do
        
    end subroutine format_csv_line

end module pod_data_format_module
