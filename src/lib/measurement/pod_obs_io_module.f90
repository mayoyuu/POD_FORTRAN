!> @file pod_obs_io_module.f90
!> @brief 解析光学观测 .obs 文件和测站 site.json 文件
module pod_obs_io_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_measurement_base_module, only: observation_station
    use pod_basicmath_module, only: PI
    use pod_time_module, only: utc_to_et
    implicit none
    private
    
    public :: load_single_observation
    public :: obs_record, station_record, ref_orbit_record
    public :: preload_observations, preload_stations, find_station_by_id
    public :: preload_reference_orbits
    public :: find_reference_by_et

    !> 轻量级观测记录（仅含滤波所需字段）
    type :: obs_record
        real(DP) :: et          ! 观测历元（TDB 秒）
        real(DP) :: ra, dec     ! 光学测量 [rad]
        character(len=16) :: station_id  ! 测站标识符
    end type obs_record

    !> 测站字典条目：ID 与完整测站信息的捆绑
    type :: station_record
        character(len=16) :: id
        type(observation_station) :: station
    end type station_record

    !> 参考轨道记录（与 obs_record 对称）
    type :: ref_orbit_record
        real(DP) :: et              ! 历元 (TDB 秒)
        real(DP) :: state(6)        ! [X, Y, Z, Vx, Vy, Vz] (km, km/s)
    end type ref_orbit_record

contains

    !> 读取给定索引的一行观测，并查找对应的 JSON 测站信息
    subroutine load_single_observation(obs_file, json_file, line_num, et, ra_rad, dec_rad, station, is_eof)
        character(len=*), intent(in) :: obs_file, json_file
        integer, intent(in)          :: line_num
        real(DP), intent(out)        :: et, ra_rad, dec_rad
        type(observation_station), intent(out) :: station
        logical, intent(out)         :: is_eof
        
        integer :: u_obs, ios, i
        character(len=MAX_STRING_LEN) :: line, sys, time_sys, site_id
        character(len=40) :: utc_str
        character(len=20) :: sec_str
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
        read(line, *, iostat=ios) sys, time_sys, year, month, day, hour, min, sec, &
                                  ra_deg, dec_deg, site_id, dummy1, dummy2
        if (ios /= 0) stop "[ERROR] 解析 OBS 行失败"
        
        ! 转换角度为弧度
        ra_rad  = ra_deg * PI / 180.0_DP
        dec_rad = dec_deg * PI / 180.0_DP

        call normalize_utc_seconds(year, month, day, hour, min, sec)
        write(sec_str, '(F12.6)') sec
        write(utc_str, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",A)') &
            year, month, day, hour, min, trim(adjustl(sec_str))

        ! 调用 pod_time_module 转换历元 (UTC -> TDB)
        et = utc_to_et(trim(utc_str))
        
        ! 3. 在 JSON 中查找对应的测站位置
        call parse_site_json(json_file, trim(site_id), station)
        
    end subroutine load_single_observation

     !> ---------------------------------------------------------------
    !> 核心优化：一次性将 OBS 文件全部解析到内存
    !> ---------------------------------------------------------------
    subroutine preload_observations(obs_file, obs_list)
        character(len=*), intent(in) :: obs_file
        type(obs_record), allocatable, intent(out) :: obs_list(:)

        integer :: u, ios, n, i
        character(len=MAX_STRING_LEN) :: line
        character(len=8) :: sys, time_sys
        integer :: year, month, day, hour, min
        real(DP) :: sec, ra_deg, dec_deg
        character(len=16) :: site_id
        real(DP) :: dummy1, dummy2
        character(len=40) :: utc_str
        character(len=20) :: sec_str

        ! 1. 统计行数
        open(newunit=u, file=obs_file, status='old', iostat=ios)
        if (ios /= 0) stop "[ERROR] 无法打开 OBS 文件"
        n = 0
        do
            read(u, '(A)', iostat=ios) line
            if (ios < 0) exit  ! EOF
            if (ios > 0) stop "[ERROR] 读取 OBS 文件出错"
            n = n + 1
        end do
        rewind(u)

        allocate(obs_list(n))

        ! 2. 逐行解析，填充结构体数组
        do i = 1, n
            read(u, '(A)') line
            read(line, *, iostat=ios) sys, time_sys, year, month, day, hour, min, sec, &
                                       ra_deg, dec_deg, site_id, dummy1, dummy2
            if (ios /= 0) then
                write(*,*) "[ERROR] 解析 OBS 行失败，行号：", i
                stop
            end if

            ! 角度 -> 弧度
            obs_list(i)%ra  = ra_deg * PI / 180.0_DP
            obs_list(i)%dec = dec_deg * PI / 180.0_DP
            obs_list(i)%station_id = trim(site_id)

            ! 构造 UTC 字符串并转为 ET
            call normalize_utc_seconds(year, month, day, hour, min, sec)
            write(sec_str, '(F12.6)') sec
            write(utc_str, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",A)') &
                year, month, day, hour, min, trim(adjustl(sec_str))
            obs_list(i)%et = utc_to_et(trim(utc_str))
        end do
        close(u)

    end subroutine preload_observations

    !> ---------------------------------------------------------------
    !> 一次性解析 site.json，构建测站 ID -> 站信息的快速查找表
    !> ---------------------------------------------------------------
    subroutine preload_stations(json_file, station_list)
        character(len=*), intent(in) :: json_file
        type(station_record), allocatable, intent(out) :: station_list(:)

        integer :: u, ios, count, i
        character(len=MAX_STRING_LEN) :: line
        character(len=16) :: site_id
        integer :: idx_id, idx_lbh, idx_end
        character(len=100) :: coord_str
        real(DP) :: lon_deg, lat_deg, alt_m

        open(newunit=u, file=json_file, status='old', iostat=ios)
        if (ios /= 0) stop "[ERROR] 无法打开 site.json"

        count = 0
        ! 先计数有效测站行（同时包含 "id" 和 "lbh"）
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, '"id"') > 0 .and. index(line, '"lbh"') > 0) count = count + 1
        end do
        rewind(u)

        allocate(station_list(count))

        i = 0
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit

            idx_id = index(line, '"id"')
            idx_lbh = index(line, '"lbh"')
            if (idx_id > 0 .and. idx_lbh > 0) then
                i = i + 1

                ! 提取 ID
                idx_id = idx_id + 6
                site_id = line(idx_id : idx_id + index(line(idx_id:), '"') - 2)
                station_list(i)%id = trim(site_id)

                ! 提取 lbh 数组内容
                idx_lbh = idx_lbh + 7
                idx_end = index(line(idx_lbh:), ']') + idx_lbh - 2
                coord_str = line(idx_lbh:idx_end)
                read(coord_str, *) lon_deg, lat_deg, alt_m

                ! 构建测站对象
                station_list(i)%station%name       = trim(site_id)
                station_list(i)%station%station_type = 'OPTICAL'
                station_list(i)%station%longitude  = lon_deg * PI / 180.0_DP
                station_list(i)%station%latitude   = lat_deg * PI / 180.0_DP
                station_list(i)%station%altitude   = alt_m / 1000.0_DP   ! m -> km

                ! ECEF 坐标（km）
                call lbh_to_ecef(station_list(i)%station%longitude, &
                                 station_list(i)%station%latitude,  &
                                 station_list(i)%station%altitude,  &
                                 station_list(i)%station%ecef_position)
            end if
        end do
        close(u)

    end subroutine preload_stations

    !> ---------------------------------------------------------------
    !> 从预加载的测站列表中按 ID 查找（线性搜索，足够高效）
    !> ---------------------------------------------------------------
    function find_station_by_id(station_id, station_list) result(station)
        character(len=*), intent(in) :: station_id
        type(station_record), intent(in) :: station_list(:)
        type(observation_station) :: station
        integer :: j

        do j = 1, size(station_list)
            if (trim(station_list(j)%id) == trim(station_id)) then
                station = station_list(j)%station
                return
            end if
        end do

        write(*,*) "[ERROR] 在预加载列表中未找到测站 ID: ", trim(station_id)
        stop
    end function find_station_by_id


    !> 轻量级手写 JSON 扫描器
    subroutine find_reference_by_et(et_query, ref_list, ref_state, found, tolerance, matched_index)
        real(DP), intent(in) :: et_query
        type(ref_orbit_record), intent(in) :: ref_list(:)
        real(DP), intent(out) :: ref_state(6)
        logical, intent(out) :: found
        real(DP), intent(in), optional :: tolerance
        integer, intent(out), optional :: matched_index

        real(DP) :: tol, best_dt, dt
        integer :: i, best_i

        tol = 0.5_DP
        if (present(tolerance)) tol = tolerance

        found = .false.
        ref_state = 0.0_DP
        best_dt = huge(1.0_DP)
        best_i = 0

        do i = 1, size(ref_list)
            dt = abs(ref_list(i)%et - et_query)
            if (dt < best_dt) then
                best_dt = dt
                best_i = i
            end if
        end do

        if (best_i > 0 .and. best_dt <= tol) then
            found = .true.
            ref_state = ref_list(best_i)%state
        end if

        if (present(matched_index)) matched_index = best_i
    end subroutine find_reference_by_et
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
        station%altitude  = alt_m/1000.0_DP  ! 转换为 km
        
        ! 将经纬高 (LBH) 转换为地固系坐标 (ECEF XYZ)
        call lbh_to_ecef(station%longitude, station%latitude, station%altitude, station%ecef_position)

        
    end subroutine parse_site_json

    subroutine lbh_to_ecef(lon_rad, lat_rad, alt_m, ecef_position)
        real(DP), intent(in) :: lon_rad, lat_rad, alt_m
        real(DP), dimension(3), intent(out) :: ecef_position
        
        real(DP) :: ae, e2, N
        
        ! WGS84 椭球参数
        ae = 6378137.0_DP/1000.0_DP  ! 地球赤道半径，单位 km
        e2 = 0.00669437999014_DP
        
        ! 计算曲率半径 N
        N = ae / sqrt(1.0_DP - e2 * sin(lat_rad)**2)
        
        ! 计算 ECEF 坐标 (X, Y, Z)
        ecef_position(1) = (N + alt_m) * cos(lat_rad) * cos(lon_rad)
        ecef_position(2) = (N + alt_m) * cos(lat_rad) * sin(lon_rad)
        ecef_position(3) = (N * (1.0_DP - e2) + alt_m) * sin(lat_rad)
    end subroutine lbh_to_ecef

    !> ---------------------------------------------------------------
    !> 一次性将 ORBITS_REF 文件全部解析到内存
    !> 格式: Name Sys YYYY MM DD HH MM SS.SSSSSS X Y Z Vx Vy Vz
    !> ---------------------------------------------------------------
    subroutine preload_reference_orbits(ref_file, ref_list)
        character(len=*), intent(in) :: ref_file
        type(ref_orbit_record), allocatable, intent(out) :: ref_list(:)

        integer :: u, ios, n, i
        character(len=MAX_STRING_LEN) :: line
        character(len=24) :: obj_name, time_sys
        integer :: year, month, day, hour, min
        real(DP) :: sec, X, Y, Z, Vx, Vy, Vz
        character(len=40) :: utc_str
        character(len=20) :: sec_str

        ! 1. 统计行数
        open(newunit=u, file=ref_file, status='old', iostat=ios)
        if (ios /= 0) stop "[ERROR] 无法打开 REF 文件: " // trim(ref_file)
        n = 0
        do
            read(u, '(A)', iostat=ios) line
            if (ios < 0) exit
            if (ios > 0) stop "[ERROR] 读取 REF 文件出错"
            n = n + 1
        end do
        rewind(u)

        allocate(ref_list(n))

        ! 2. 逐行解析
        do i = 1, n
            read(u, '(A)') line
            read(line, *, iostat=ios) obj_name, time_sys, &
                                       year, month, day, hour, min, sec, &
                                       X, Y, Z, Vx, Vy, Vz
            if (ios /= 0) then
                write(*,*) "[ERROR] 解析 REF 行失败，行号：", i
                stop
            end if

            ref_list(i)%state(1) = X
            ref_list(i)%state(2) = Y
            ref_list(i)%state(3) = Z
            ref_list(i)%state(4) = Vx
            ref_list(i)%state(5) = Vy
            ref_list(i)%state(6) = Vz

            ! 替换为新代码：
            call normalize_utc_seconds(year, month, day, hour, min, sec)
            write(sec_str, '(F12.6)') sec
            write(utc_str, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",A)') &
                year, month, day, hour, min, trim(adjustl(sec_str))
            ref_list(i)%et = utc_to_et(trim(utc_str))

        end do
        close(u)

    end subroutine preload_reference_orbits

   !> 秒数接近 60 时进位归一化，避免 SPICE STR2ET 报错
    subroutine normalize_utc_seconds(year, month, day, hour, min, sec)
        integer, intent(inout)  :: year, month, day, hour, min
        real(DP), intent(inout) :: sec

        ! 检查是否大于或等于 59.9999995，防止写入字符串时四舍五入变为 60.0
        if (sec >= 59.9999995_DP) then
            sec = 0.0_DP
            min = min + 1
            if (min >= 60) then
                min = min - 60
                hour = hour + 1
                if (hour >= 24) then
                    hour = hour - 24
                    day = day + 1
                end if
            end if
        end if
    end subroutine normalize_utc_seconds

end module pod_obs_io_module