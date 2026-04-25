!> @file pod_data_format_module.f90
!> @brief 负责轨道参数消息 (OPM) 与 JSON 格式的序列化与反序列化
module pod_data_format_module
    use pod_global, only: DP, MAX_STRING_LEN
    ! 如果你封装了 SPICE 的时间转换，从这里引入：
    use pod_spice, only: str2et, et2utc

    use pod_uq_gmm_state_module, only: uq_gmm_state_type
     ! 其他必要的模块...
    
    implicit none
    private
    
    public :: load_initial_opm
    public :: write_json_opm

contains


   !> ======================================================================
    !> 读取标准的 OPM JSON 格式并初始化状态
    !> ======================================================================
    subroutine load_initial_opm(json_file, et, state, cov, gmm_state, has_gmm)
        character(len=*), intent(in) :: json_file
        real(DP), intent(out)        :: et
        real(DP), intent(out)        :: state(6)
        real(DP), intent(out)        :: cov(6,6)
        type(uq_gmm_state_type), intent(out), optional :: gmm_state
        logical, intent(out), optional :: has_gmm

        integer :: u_json, ios, i, j, k
        integer :: comp_idx, n_components
        real(DP) :: tmp_val
        character(len=MAX_STRING_LEN) :: line
        character(len=64) :: epoch_str
        logical :: found
        
        character(len=5), parameter :: s_keys(6) = ["X    ", "Y    ", "Z    ", "X_DOT", "Y_DOT", "Z_DOT"]
        
        state = 0.0_DP
        cov = 0.0_DP
        et = 0.0_DP
        epoch_str = ""
        comp_idx = 0
        if (present(has_gmm)) has_gmm = .false.
        
        open(newunit=u_json, file=json_file, status='old', iostat=ios)
        if (ios /= 0) stop "[ERROR] 无法打开初始 OPM 文件: " // trim(json_file)
        
        do
            read(u_json, '(A)', iostat=ios) line
            if (ios < 0) exit ! 文件结束
            
            ! 提取时间
            call extract_json_string(line, '"EPOCH"', epoch_str, found)
            
            ! 提取状态向量与主协方差矩阵 (自动处理对称性)
            do i = 1, 6
                call extract_json_value(line, '"'//trim(s_keys(i))//'"', state(i), found)
                do j = 1, i
                    call extract_json_value(line, '"C'//trim(s_keys(i))//'_'//trim(s_keys(j))//'"', cov(i,j), found)
                    if (found) cov(j,i) = cov(i,j) 
                end do
            end do
            
           ! 提取 GMM 状态
            if (present(gmm_state)) then
                ! 1. 解析总分量数并初始化
                call extract_json_value(line, '"GMM_N_COMPONENTS"', tmp_val, found)
                if (found) then
                    n_components = nint(tmp_val)
                    
                    ! [修正1] 调用模块中真实存在的方法分配主容器
                    call gmm_state%allocate_components(n_components, 6)
                    
                    ! [修正2] 极其重要！必须立即为每个分量的内部 allocatable 数组分配内存
                    do i = 1, n_components
                        allocate(gmm_state%components(i)%mean(6))
                        allocate(gmm_state%components(i)%cov(6,6))
                        gmm_state%components(i)%weight = 0.0_DP
                        gmm_state%components(i)%mean   = 0.0_DP
                        gmm_state%components(i)%cov    = 0.0_DP
                    end do
                    
                    if (present(has_gmm)) has_gmm = .true.
                end if
                
                ! 2. 解析当前块的 INDEX
                call extract_json_value(line, '"INDEX"', tmp_val, found)
                if (found) comp_idx = nint(tmp_val)
                
                ! 3. 填充对应 INDEX 的属性
                if (comp_idx > 0 .and. comp_idx <= gmm_state%n_components) then
                    call extract_json_value(line, '"WEIGHT"', gmm_state%components(comp_idx)%weight, found)
                    
                    ! 此时 mean 已经 allocate，可以直接安全地接收数据
                    call extract_json_array6(line, '"MEAN"', gmm_state%components(comp_idx)%mean, found)
                    
                    ! 提取 GMM 全协方差矩阵 (读取下三角，对称赋值)
                    do j = 1, 6
                        do k = 1, j
                            call extract_json_value(line, '"COV_' // char(48+j) // '_' // char(48+k) // '"', tmp_val, found)
                            if (found) then
                                ! 此时 cov 已经 allocate，直接按索引赋值是安全的
                                gmm_state%components(comp_idx)%cov(j,k) = tmp_val
                                gmm_state%components(comp_idx)%cov(k,j) = tmp_val ! 保证对称矩阵
                            end if
                        end do
                    end do
                end if
            end if
        end do
        
        close(u_json)
        
        if (trim(epoch_str) /= "") call str2et(trim(epoch_str), et)
        
    end subroutine load_initial_opm

    !> ======================================================================
    !> 将滤波后的最终结果输出为标准的 OPM JSON 格式
    !> ======================================================================
    subroutine write_json_opm(filename, final_state, final_cov, gmm_state, rms, obj_id, et_last)
        character(len=*), intent(in) :: filename, obj_id
        real(DP), intent(in)         :: final_state(6), final_cov(6,6), rms, et_last
        type(uq_gmm_state_type), intent(in), optional :: gmm_state
        
        integer :: u, i, j, k
        character(len=64) :: epoch_str
        
        ! 定义状态键名常量，用于自动生成状态名与协方差名
        character(len=5), parameter :: s_keys(6) = ["X    ", "Y    ", "Z    ", "X_DOT", "Y_DOT", "Z_DOT"]
        
        ! 将 TDB 秒转换回 UTC 字符串
        call et2utc(et_last, 'ISOC', 6, epoch_str) 
        
        open(newunit=u, file=filename, status='replace', action='write')
        
        write(u, '(A)') '{'
        write(u, '(A,A,A)') '    "ASL_CAT_ID": "', trim(obj_id), '",'
        
        ! 动力学参数标称配置
        write(u, '(A)') '    "DRAG_AREA": 1.0,'
        write(u, '(A)') '    "DRAG_CD": 2.2,'
        write(u, '(A)') '    "MASS": 1000.0,'
        write(u, '(A)') '    "SOLAR_RAD_AREA": 1.0,'
        write(u, '(A)') '    "SOLAR_RAD_COEFF": 1.3,'
        
        ! 历元与状态向量
        write(u, '(A,A,A)') '    "EPOCH": "', trim(epoch_str), '",'
        do i = 1, 6
            write(u, '(A,A,A,ES22.15,A)') '    "', trim(s_keys(i)), '": ', final_state(i), ','
        end do
        
        ! 协方差矩阵 (自动生成 OPM 标准的下三角键名: C + row_name + _ + col_name)
        do i = 1, 6
            do j = 1, i
                write(u, '(A,A,A,A,A,ES22.15,A)') '    "C', trim(s_keys(i)), '_', trim(s_keys(j)), '": ', final_cov(i,j), ','
            end do
        end do
        
       ! 写入 GMM 核心描述
        if (present(gmm_state)) then
            write(u, '(A,I0,A)') '    "GMM_N_COMPONENTS": ', gmm_state%n_components, ','
            write(u, '(A)') '    "GMM_COMPONENTS": ['
            
            do i = 1, gmm_state%n_components
                write(u, '(A)') '        {'
                write(u, '(A,I0,A)')          '            "INDEX": ', i, ','
                write(u, '(A,ES22.15,A)')     '            "WEIGHT": ', gmm_state%components(i)%weight, ','
                write(u, '(A, 5(ES22.15, ", "), ES22.15, A)') '            "MEAN": [', gmm_state%components(i)%mean, '],'
                
                ! 展开输出下三角协方差矩阵 (COV_1_1 到 COV_6_6)
                do j = 1, 6
                    do k = 1, j
                        if (j == 6 .and. k == 6) then
                            ! 最后一项不加逗号，遵循严格的 JSON 语法
                            write(u, '(A,A,A,A,A,ES22.15)') '            "COV_', char(48+j), '_', char(48+k), '": ',&
                             gmm_state%components(i)%cov(j,k)
                        else
                            write(u, '(A,A,A,A,A,ES22.15,A)') '            "COV_', char(48+j), '_', char(48+k), '": ', &
                            gmm_state%components(i)%cov(j,k), ','
                        end if
                    end do
                end do
                
                if (i < gmm_state%n_components) then
                    write(u, '(A)') '        },'
                else
                    write(u, '(A)') '        }'
                end if
            end do
            write(u, '(A)') '    ],'
        end if
        
        ! 统计信息
        write(u, '(A,ES22.15,A)') '    "RMS": ', rms, ','
        write(u, '(A,A,A)')       '    "LAST_OBS": "', trim(epoch_str), '",'
        write(u, '(A)')           '    "STATUS": "SUCCESS"'
        
        write(u, '(A)') '}'
        close(u)
        
    end subroutine write_json_opm

    !> ======================================================================
    !> 内部辅助工具：从 JSON 字符串行中提取实数值
    !> ======================================================================
    subroutine extract_json_value(line, key, val, found)
        character(len=*), intent(in)  :: line, key
        real(DP), intent(out)         :: val
        logical, intent(out)          :: found
        
        integer :: idx_key, idx_colon, idx_comma
        character(len=MAX_STRING_LEN) :: val_str
        
        found = .false.
        idx_key = index(line, trim(key))
        if (idx_key > 0) then
            idx_colon = index(line(idx_key:), ':') + idx_key - 1
            if (idx_colon >= idx_key) then
                idx_comma = index(line(idx_colon:), ',')
                if (idx_comma > 0) then
                    val_str = line(idx_colon+1 : idx_colon+idx_comma-2)
                else
                    val_str = line(idx_colon+1 :) ! 没有逗号，说明是最后一项
                end if
                read(val_str, *, err=10) val
                found = .true.
            end if
        end if
        return
        
10      found = .false. ! 读取浮点数失败的异常处理
    end subroutine extract_json_value

    !> ======================================================================
    !> 内部辅助工具：从 JSON 字符串行中提取字符串值
    !> ======================================================================
    subroutine extract_json_string(line, key, str_val, found)
        character(len=*), intent(in)  :: line, key
        character(len=*), intent(out) :: str_val
        logical, intent(out)          :: found
        
        integer :: idx_key, idx_colon, start_quote, end_quote
        
        found = .false.
        idx_key = index(line, trim(key))
        if (idx_key > 0) then
            idx_colon = index(line(idx_key:), ':') + idx_key - 1
            start_quote = index(line(idx_colon:), '"') + idx_colon - 1
            if (start_quote >= idx_colon) then
                end_quote = index(line(start_quote+1:), '"') + start_quote
                if (end_quote > start_quote) then
                    str_val = line(start_quote+1 : end_quote-1)
                    found = .true.
                end if
            end if
        end if
    end subroutine extract_json_string


    !> ======================================================================
    !> 内部辅助工具：从 JSON 数组字符串提取 6 个实数 (如 MEAN 或 COV_DIAG)
    !> ======================================================================
    subroutine extract_json_array6(line, key, vals, found)
        character(len=*), intent(in)  :: line, key
        real(DP), intent(out)         :: vals(6)
        logical, intent(out)          :: found
        
        integer :: idx_key, idx_bracket1, idx_bracket2
        character(len=MAX_STRING_LEN) :: arr_str
        
        found = .false.
        idx_key = index(line, trim(key))
        if (idx_key > 0) then
            idx_bracket1 = index(line(idx_key:), '[') + idx_key - 1
            idx_bracket2 = index(line(idx_bracket1:), ']') + idx_bracket1 - 1
            
            if (idx_bracket1 >= idx_key .and. idx_bracket2 > idx_bracket1) then
                arr_str = line(idx_bracket1+1 : idx_bracket2-1)
                ! Fortran 的表控格式输入(*) 会自动将逗号视为空白分隔符
                read(arr_str, *, err=20) vals
                found = .true.
            end if
        end if
        return
20      found = .false.
    end subroutine extract_json_array6


end module pod_data_format_module