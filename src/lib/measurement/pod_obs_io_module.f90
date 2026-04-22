!> @file pod_obs_io_module.f90
!> @brief 解析光学观测 .obs 文件和测站 site.json 文件
module pod_obs_io_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_measurement_base_module, only: observation_station
    use pod_basicmath_module, only: PI
    implicit none
    private
    
    public :: load_single_observation
    
contains

    !> 读取给定索引的一行观测，并查找对应的 JSON 测站信息
    subroutine load_single_observation(obs_file, json_file, line_num, et, ra_rad, dec_rad, station, is_eof)
        character(len=*), intent(in) :: obs_file, json_file
        integer, intent(in)          :: line_num
        real(DP), intent(out)        :: et, ra_rad, dec_rad
        type(observation_station), intent(out) :: station
        logical, intent(out)         :: is_eof
        
        integer :: u_obs, ios, i
        character(len=MAX_STRING_LEN) :: line, sys, site_id
        integer :: year, month, day, hour, min
        real(DP) :: sec, ra_deg, dec_deg, dummy1, dummy2
        
        is_eof = .false.
        
        ! 1. 扫描 OBS 文件到达指定行
        open(newunit=u_obs, file=obs_file, status='old', iostat=ios)
        if (ios /= 0) stop "[ERROR] 无法打开 OBS 文件"
        
        do i = 1, line_num
            read(u_obs, '(A)', iostat=ios) line
            if (ios < 0) then
                is_eof = .true.
                close(u_obs)
                return
            end if
        end do
        close(u_obs)
        
        ! 2. 解析 OBS 格式: UTC YYYY MM DD HH MM SS.SSS RA DEC SITE 0.0 0.0
        read(line, *, iostat=ios) sys, year, month, day, hour, min, sec, &
                                  ra_deg, dec_deg, site_id, dummy1, dummy2
        if (ios /= 0) stop "[ERROR] 解析 OBS 行失败"
        
        ! 转换角度为弧度
        ra_rad  = ra_deg * PI / 180.0_DP
        dec_rad = dec_deg * PI / 180.0_DP
        
        ! TODO: 调用你的时间转换库将 YYYY MM DD HH MM SS.SSS 转为历元 et (TDB秒)
        ! et = calendar_to_tdb(year, month, day, hour, min, sec)
        et = 0.0_DP ! 占位
        
        ! 3. 在 JSON 中查找对应的测站位置
        call parse_site_json(json_file, trim(site_id), station)
        
    end subroutine load_single_observation

    !> 轻量级手写 JSON 扫描器
    subroutine parse_site_json(json_file, target_id, station)
        character(len=*), intent(in) :: json_file, target_id
        type(observation_station), intent(out) :: station
        
        integer :: u_json, ios, idx_id, idx_lbh
        character(len=MAX_STRING_LEN) :: line, search_str
        character(len=100) :: coord_str
        real(DP) :: lon_deg, lat_deg, alt_m
        logical :: found
        
        station%name = target_id
        station%station_type = 'OPTICAL'
        found = .false.
        
        ! 目标搜索串： "id":"R92"
        search_str = '"id":"' // trim(target_id) // '"'
        
        open(newunit=u_json, file=json_file, status='old', iostat=ios)
        if (ios /= 0) stop "[ERROR] 无法打开 site.json 文件"
        
        do
            read(u_json, '(A)', iostat=ios) line
            if (ios < 0) exit
            
            idx_id = index(line, trim(search_str))
            if (idx_id > 0) then
                ! 找到了对应的 site, 提取 "lbh":[lon, lat, alt]
                idx_lbh = index(line, '"lbh":[')
                if (idx_lbh > 0) then
                    ! 截取中括号内的字符串
                    coord_str = line(idx_lbh + 7 : index(line, ']') - 1)
                    ! 从提取的字符串中读取以逗号分隔的三个浮点数
                    read(coord_str, *) lon_deg, lat_deg, alt_m
                    found = .true.
                    exit
                end if
            end if
        end do
        close(u_json)
        
        if (.not. found) then
            write(*,*) "[ERROR] 在 json 中找不到测站 ID: ", trim(target_id)
            stop
        end if
        
        ! 转为弧度并存入 station 结构体
        station%longitude = lon_deg * PI / 180.0_DP
        station%latitude  = lat_deg * PI / 180.0_DP
        station%altitude  = alt_m
        
        ! 将经纬高 (LBH) 转换为地固系坐标 (ECEF XYZ)
        call lbh_to_ecef(station%longitude, station%latitude, station%altitude, station%ecef_position)
        
    end subroutine parse_site_json

    subroutine lbh_to_ecef(lon_rad, lat_rad, alt_m, ecef_position)
        real(DP), intent(in) :: lon_rad, lat_rad, alt_m
        real(DP), dimension(3), intent(out) :: ecef_position
        
        real(DP) :: ae, e2, N
        
        ! WGS84 椭球参数
        ae = 6378137.0_DP
        e2 = 0.00669437999014_DP
        
        ! 计算曲率半径 N
        N = ae / sqrt(1.0_DP - e2 * sin(lat_rad)**2)
        
        ! 计算 ECEF 坐标 (X, Y, Z)
        ecef_position(1) = (N + alt_m) * cos(lat_rad) * cos(lon_rad)
        ecef_position(2) = (N + alt_m) * cos(lat_rad) * sin(lon_rad)
        ecef_position(3) = (N * (1.0_DP - e2) + alt_m) * sin(lat_rad)
    end subroutine lbh_to_ecef

end module pod_obs_io_module