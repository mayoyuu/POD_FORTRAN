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
    !> 读取初始 OPM JSON 文件，提取状态向量和协方差矩阵
    !> 假设 JSON 是格式化的 (Pretty-Printed)，即每行包含独立的键值对
    !> ======================================================================
    subroutine load_initial_opm(json_file, et, state, cov ,gmm_state, has_gmm)
        character(len=*), intent(in) :: json_file
        real(DP), intent(out)        :: et
        real(DP), intent(out)        :: state(6)
        real(DP), intent(out)        :: cov(6,6)
        type(uq_gmm_state_type), intent(out), optional :: gmm_state
        logical, intent(out), optional :: has_gmm

        integer :: u_json, ios
        character(len=MAX_STRING_LEN) :: line
        character(len=64) :: epoch_str
        logical :: found
        
        ! 初始全赋零/默认值
        state = 0.0_DP
        cov = 0.0_DP
        et = 0.0_DP
        epoch_str = ""
        
        open(newunit=u_json, file=json_file, status='old', iostat=ios)
        if (ios /= 0) stop "[ERROR] 无法打开初始 OPM 文件: " // trim(json_file)
        
        do
            read(u_json, '(A)', iostat=ios) line
            if (ios < 0) exit ! 文件结束
            
            ! 提取时间 (字符串)
            call extract_json_string(line, '"EPOCH"', epoch_str, found)
            
            ! 提取状态向量 (X, Y, Z, X_DOT, Y_DOT, Z_DOT)
            call extract_json_value(line, '"X"', state(1), found)
            call extract_json_value(line, '"Y"', state(2), found)
            call extract_json_value(line, '"Z"', state(3), found)
            call extract_json_value(line, '"X_DOT"', state(4), found)
            call extract_json_value(line, '"Y_DOT"', state(5), found)
            call extract_json_value(line, '"Z_DOT"', state(6), found)
            
            ! 提取协方差矩阵对角线 (示例：CX_X, CY_Y...)
            ! (如果初始文件包含完整的 21 个独立项，可按此模式继续补充提取)
            call extract_json_value(line, '"CX_X"', cov(1,1), found)
            call extract_json_value(line, '"CY_Y"', cov(2,2), found)
            call extract_json_value(line, '"CZ_Z"', cov(3,3), found)
            call extract_json_value(line, '"CX_DOT_X_DOT"', cov(4,4), found)
            call extract_json_value(line, '"CY_DOT_Y_DOT"', cov(5,5), found)
            call extract_json_value(line, '"CZ_DOT_Z_DOT"', cov(6,6), found)
        end do
        
        close(u_json)
        
        ! 调用 SPICE 将 ISO 8601 字符串转换为 TDB 秒
        call str2et(trim(epoch_str), et)
        
    end subroutine load_initial_opm

    !> ======================================================================
    !> 将滤波后的最终结果输出为标准的 OPM JSON 格式
    !> ======================================================================
    subroutine write_json_opm(filename, final_state, final_cov, rms, obj_id, et_last)
        character(len=*), intent(in) :: filename, obj_id
        real(DP), intent(in)         :: final_state(6), final_cov(6,6), rms, et_last
        
        integer :: u
        character(len=64) :: epoch_str
        
        ! 将 TDB 秒转换回 UTC 字符串
        call et2utc(et_last, 'ISOC', 6, epoch_str) 
        
        open(newunit=u, file=filename, status='replace', action='write')
        
        write(u, '(A)') '{'
        write(u, '(A,A,A)') '    "ASL_CAT_ID": "', trim(obj_id), '",'
        
        ! 动力学参数 (如果是估值可替换为动态变量，此处输出标称配置)
        write(u, '(A)') '    "DRAG_AREA": 1.0,'
        write(u, '(A)') '    "DRAG_CD": 2.2,'
        write(u, '(A)') '    "MASS": 1000.0,'
        write(u, '(A)') '    "SOLAR_RAD_AREA": 1.0,'
        write(u, '(A)') '    "SOLAR_RAD_COEFF": 1.3,'
        
        ! 历元与状态向量
        write(u, '(A,A,A)') '    "EPOCH": "', trim(epoch_str), '",'
        write(u, '(A,ES22.15,A)') '    "X": ', final_state(1), ','
        write(u, '(A,ES22.15,A)') '    "Y": ', final_state(2), ','
        write(u, '(A,ES22.15,A)') '    "Z": ', final_state(3), ','
        write(u, '(A,ES22.15,A)') '    "X_DOT": ', final_state(4), ','
        write(u, '(A,ES22.15,A)') '    "Y_DOT": ', final_state(5), ','
        write(u, '(A,ES22.15,A)') '    "Z_DOT": ', final_state(6), ','
        
        ! 协方差矩阵 (只输出下三角或全展开)
        write(u, '(A,ES22.15,A)') '    "CX_X": ', final_cov(1,1), ','
        write(u, '(A,ES22.15,A)') '    "CY_Y": ', final_cov(2,2), ','
        write(u, '(A,ES22.15,A)') '    "CZ_Z": ', final_cov(3,3), ','
        write(u, '(A,ES22.15,A)') '    "CX_DOT_X_DOT": ', final_cov(4,4), ','
        write(u, '(A,ES22.15,A)') '    "CY_DOT_Y_DOT": ', final_cov(5,5), ','
        write(u, '(A,ES22.15,A)') '    "CZ_DOT_Z_DOT": ', final_cov(6,6), ','
        ! 如果需要，可以继续写入交叉协方差如 "CX_DOT_X": final_cov(4,1) 等...
        
        ! 统计信息
        write(u, '(A,ES22.15,A)') '    "RMS": ', rms, ','
        write(u, '(A,A,A)') '    "LAST_OBS": "', trim(epoch_str), '",'
        write(u, '(A)') '    "STATUS": "SUCCESS"'
        
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

    !> 从 JSON 文件中探测并加载 GMM 结构
    subroutine load_gmm_state_from_json(filename, gmm_out, et, has_gmm)
        character(len=*), intent(in) :: filename
        type(uq_gmm_state_type), intent(out) :: gmm_out
        real(DP), intent(out) :: et
        logical, intent(out) :: has_gmm
        
        ! 逻辑：
        ! 1. 打开文件扫描 "GMM_N_COMPONENTS" 关键字
        ! 2. 如果没有，has_gmm = .false.
        ! 3. 如果有，按照上一轮讨论的数组解析逻辑，allocate 并填充 gmm_out
    end subroutine load_gmm_state_from_json

    !> 将当前的 GMM 状态完整存入 JSON (供后续任务断点续传)
    subroutine save_gmm_state_to_json(filename, gmm_state, et, obj_id)
        character(len=*), intent(in) :: filename, obj_id
        type(uq_gmm_state_type), intent(in) :: gmm_state
        real(DP), intent(in) :: et
        ! 逻辑：遍历 gmm_state%components 写入 JSON 数组
    end subroutine save_gmm_state_to_json

end module pod_data_format_module