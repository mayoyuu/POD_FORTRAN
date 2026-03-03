!> # POD Frame Module
!>
!> 统一的坐标系/参考系（frame）转换与相关几何变换模块。
!>
!> 本模块参考时间系统模块的设计风格，提供对外一致的接口，并与
!> SPICE 工具包深度集成：所有帧标识统一使用 SPICE 帧名字符串
!>（不在程序内自定义/维护帧 ID），必要时通过 `namfrm` 动态查询
!> 帧 ID；具体帧间旋转由 `pxform`/`sxform` 获得。
!>
!> ## Features
!> - **Frame 转换**: 基于 SPICE `pxform`/`sxform` 的帧间旋转/状态转换
!> - **轨道要素 ↔ 笛卡尔**: `cartesian_to_keplerian` / `keplerian_to_cartesian`
!> - **几何工具**: 旋转矩阵构造、角度弧度转换等（按需扩展）
!> - **高精度**: 使用双精度实数计算
!>
!> ## Supported Frames (使用帧名，不维护本地ID)
!> - `J2000`                 : 地心惯性参考系（SPICE标准）
!> - `BCRS`                  : 需要在 Frames Kernel (FK) 中定义后使用
!> - `ITRF93`/`ITRF2000`     : 地固参考系（按所用 FK/内核选择命名）
!> - `TOPOCENTRIC_*`         : 台站本地地平系（建议在 FK 中按台站命名）
!> - `EARTH-MOON ROTATION`   : 需在 FK 中定义的自定义帧名（示例名）
!>
!> 以上帧名均作为字符串传入/传出；如确需帧 ID，使用 `namfrm` 由
!> 帧名动态获取，以避免本地常量失配风险。
!>
!> ## SPICE Integration & Kernel Requirements
!> - 函数/子程序: `pxform`, `sxform`, `namfrm`（通过 `pod_spice` 暴露）
!> - 必需内核: 
!>   - LSK（闰秒内核）: 时间一致性
!>   - PCK/SPK（视需求）: 天体/姿态数据
!>   - FK（Frames Kernel）: 自定义帧（如 BCRS、台站 TOPO、EARTH-MOON ROTATION）
!>
!> ## Version
!> - **Created**: 2025-09-12
!> - **Updated**: 2025-09-12 - 统一采用 SPICE 帧名；不新增本地帧ID
!>
!> @note 所有对外接口建议使用帧名字符串（如 `from_frame`, `to_frame`）。
!> @warning 使用前需确保相应 SPICE 内核（尤其 FK）已正确加载。
!> @todo 提供基于帧名的高层封装接口，如 `transform_frame(from,to,et,vec)`。
module pod_frame_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_config, only: config
    use pod_frame_simple_module, only: pod_frame_state  ! <--- 引入状态对象
    use pod_spice, only: sxform, pxform, namfrm
    
    implicit none
    
    ! 数学常数
    real(DP), parameter :: PI = 4.0_DP * atan(1.0_DP)
    real(DP), parameter :: DEG_TO_RAD = PI / 180.0_DP
    real(DP), parameter :: RAD_TO_DEG = 180.0_DP / PI
    
contains

    subroutine cartesian_to_keplerian(cartesian_state, keplerian_elements)
        real(DP), dimension(6), intent(in) :: cartesian_state
        real(DP), dimension(6), intent(out) :: keplerian_elements
        
        real(DP), dimension(3) :: r, v, h, e
        real(DP) :: r_mag, v_mag, h_mag, e_mag, energy, a, e_scalar
        real(DP) :: i, omega, w, nu, cos_nu, sin_nu, cos_w
        real(DP) :: n_x, n_y, n_mag, cos_i, sin_i
        
        ! 提取位置和速度向量
        r = cartesian_state(1:3)
        v = cartesian_state(4:6)
        
        ! 计算位置和速度大小
        r_mag = sqrt(sum(r**2))
        v_mag = sqrt(sum(v**2))
        
        ! 计算角动量向量
        h = cross_product(r, v)
        h_mag = sqrt(sum(h**2))
        
        ! 计算偏心率向量
        e = cross_product(v, h) / config%gravitational_constant - r / r_mag
        e_mag = sqrt(sum(e**2))
        
        ! 计算轨道能量
        energy = 0.5_DP * v_mag**2 - config%gravitational_constant / r_mag
        
        ! 计算半长轴
        if (abs(energy) > 1.0e-12_DP) then
            a = -config%gravitational_constant / (2.0_DP * energy)
        else
            a = huge(1.0_DP)  ! 抛物线轨道
        end if
        
        ! 计算偏心率
        e_scalar = e_mag
        
        ! 计算轨道倾角
        cos_i = h(3) / h_mag
        i = acos(cos_i) * RAD_TO_DEG
        
        ! 计算升交点赤经
        n_x = -h(2)
        n_y = h(1)
        n_mag = sqrt(n_x**2 + n_y**2)
        
        if (n_mag > 1.0e-12_DP) then
            omega = atan2(n_y, n_x) * RAD_TO_DEG
            if (omega < 0.0_DP) omega = omega + 360.0_DP
        else
            omega = 0.0_DP  ! 赤道轨道
        end if
        
        ! 计算近地点幅角
        if (n_mag > 1.0e-12_DP .and. e_scalar > 1.0e-12_DP) then
            cos_w = (n_x * e(1) + n_y * e(2)) / (n_mag * e_scalar)
            cos_w = max(-1.0_DP, min(1.0_DP, cos_w))  ! 限制在[-1, 1]范围内
            w = acos(cos_w) * RAD_TO_DEG
            
            if (e(3) < 0.0_DP) w = 360.0_DP - w
        else
            w = 0.0_DP
        end if
        
        ! 计算真近点角
        if (e_scalar > 1.0e-12_DP) then
            cos_nu = vector_dot_product(e, r) / (e_scalar * r_mag)
            cos_nu = max(-1.0_DP, min(1.0_DP, cos_nu))
            sin_nu = vector_dot_product(h, cross_product(e, r)) / (h_mag * e_scalar * r_mag)
            nu = atan2(sin_nu, cos_nu) * RAD_TO_DEG
            if (nu < 0.0_DP) nu = nu + 360.0_DP
        else
            nu = 0.0_DP
        end if
        
        ! 保存开普勒轨道根数
        keplerian_elements(1) = a
        keplerian_elements(2) = e_scalar
        keplerian_elements(3) = i
        keplerian_elements(4) = omega
        keplerian_elements(5) = w
        keplerian_elements(6) = nu
    end subroutine cartesian_to_keplerian
    
    subroutine keplerian_to_cartesian(keplerian_elements, cartesian_state)
        real(DP), dimension(6), intent(in) :: keplerian_elements
        real(DP), dimension(6), intent(out) :: cartesian_state
        
        real(DP) :: a, e, i, omega, w, nu
        real(DP) :: r_mag, v_mag, r_orbital, v_orbital
        real(DP), dimension(3) :: r_orb, v_orb, r_eci, v_eci
        real(DP), dimension(3,3) :: rotation_matrix
        
        ! 提取开普勒轨道根数
        a = keplerian_elements(1)
        e = keplerian_elements(2)
        i = keplerian_elements(3) * DEG_TO_RAD
        omega = keplerian_elements(4) * DEG_TO_RAD
        w = keplerian_elements(5) * DEG_TO_RAD
        nu = keplerian_elements(6) * DEG_TO_RAD
        
        ! 计算轨道平面内的位置和速度
        r_orbital = a * (1.0_DP - e**2) / (1.0_DP + e * cos(nu))
        v_orbital = sqrt(config%gravitational_constant / (a * (1.0_DP - e**2)))
        
        ! 轨道平面内的位置向量
        r_orb(1) = r_orbital * cos(nu)
        r_orb(2) = r_orbital * sin(nu)
        r_orb(3) = 0.0_DP
        
        ! 轨道平面内的速度向量
        v_orb(1) = -v_orbital * sin(nu)
        v_orb(2) = v_orbital * (e + cos(nu))
        v_orb(3) = 0.0_DP
        
        ! 构建旋转矩阵 (轨道平面到ECI坐标系)
        call build_rotation_matrix(omega, i, w, rotation_matrix)
        
        ! 转换到ECI坐标系
        r_eci = matmul(rotation_matrix, r_orb)
        v_eci = matmul(rotation_matrix, v_orb)
        
        ! 保存笛卡尔状态
        cartesian_state(1:3) = r_eci
        cartesian_state(4:6) = v_eci
    end subroutine keplerian_to_cartesian
    
    function cross_product(a, b) result(c)
        real(DP), dimension(3), intent(in) :: a, b
        real(DP), dimension(3) :: c
        
        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
    end function cross_product
    
    function vector_dot_product(a, b) result(c)
        real(DP), dimension(:), intent(in) :: a, b
        real(DP) :: c
        
        c = sum(a * b)
    end function vector_dot_product
    
    subroutine build_rotation_matrix(omega, i, w, rotation_matrix)
        real(DP), intent(in) :: omega, i, w
        real(DP), dimension(3,3), intent(out) :: rotation_matrix
        
        real(DP), dimension(3,3) :: R_omega, R_i, R_w
        
        ! 绕Z轴旋转 (升交点赤经)
        R_omega(1,1) = cos(omega)
        R_omega(1,2) = -sin(omega)
        R_omega(1,3) = 0.0_DP
        R_omega(2,1) = sin(omega)
        R_omega(2,2) = cos(omega)
        R_omega(2,3) = 0.0_DP
        R_omega(3,1) = 0.0_DP
        R_omega(3,2) = 0.0_DP
        R_omega(3,3) = 1.0_DP
        
        ! 绕X轴旋转 (轨道倾角)
        R_i(1,1) = 1.0_DP
        R_i(1,2) = 0.0_DP
        R_i(1,3) = 0.0_DP
        R_i(2,1) = 0.0_DP
        R_i(2,2) = cos(i)
        R_i(2,3) = -sin(i)
        R_i(3,1) = 0.0_DP
        R_i(3,2) = sin(i)
        R_i(3,3) = cos(i)
        
        ! 绕Z轴旋转 (近地点幅角)
        R_w(1,1) = cos(w)
        R_w(1,2) = -sin(w)
        R_w(1,3) = 0.0_DP
        R_w(2,1) = sin(w)
        R_w(2,2) = cos(w)
        R_w(2,3) = 0.0_DP
        R_w(3,1) = 0.0_DP
        R_w(3,2) = 0.0_DP
        R_w(3,3) = 1.0_DP
        
        ! 组合旋转矩阵: R = R_omega * R_i * R_w
        rotation_matrix = matmul(matmul(R_omega, R_i), R_w)
    end subroutine build_rotation_matrix
    
    subroutine eci_to_geodetic(eci_pos, lat, lon, alt)
        real(DP), dimension(3), intent(in) :: eci_pos
        real(DP), intent(out) :: lat, lon, alt
        
        real(DP), parameter :: EARTH_RADIUS = 6378.137_DP  ! km
        real(DP), parameter :: EARTH_FLATTENING = 1.0_DP / 298.257223563_DP
        real(DP), parameter :: EARTH_ECC2 = 2.0_DP * EARTH_FLATTENING - EARTH_FLATTENING**2
        
        real(DP) :: x, y, z, r_xy, r_mag
        real(DP) :: lat_old, lat_new, sin_lat, cos_lat, N
        integer :: i
        
        x = eci_pos(1)
        y = eci_pos(2)
        z = eci_pos(3)
        
        ! 计算经度
        lon = atan2(y, x) * RAD_TO_DEG
        
        ! 计算纬度 (迭代方法)
        r_xy = sqrt(x**2 + y**2)
        r_mag = sqrt(x**2 + y**2 + z**2)
        
        ! 初始猜测
        lat_old = asin(z / r_mag) * RAD_TO_DEG
        
        ! 迭代求解
        do i = 1, 10
            sin_lat = sin(lat_old * DEG_TO_RAD)
            cos_lat = cos(lat_old * DEG_TO_RAD)
            
            N = EARTH_RADIUS / sqrt(1.0_DP - EARTH_ECC2 * sin_lat**2)
            
            lat_new = atan2(z + N * EARTH_ECC2 * sin_lat, r_xy) * RAD_TO_DEG
            
            if (abs(lat_new - lat_old) < 1.0e-12_DP) exit
            lat_old = lat_new
        end do
        
        lat = lat_new
        
        ! 计算高度
        alt = r_xy / cos_lat - N
    end subroutine eci_to_geodetic
    
    subroutine geodetic_to_eci(lat, lon, alt, eci_pos)
        real(DP), intent(in) :: lat, lon, alt
        real(DP), dimension(3), intent(out) :: eci_pos
        
        real(DP), parameter :: EARTH_RADIUS = 6378.137_DP  ! km
        real(DP), parameter :: EARTH_FLATTENING = 1.0_DP / 298.257223563_DP
        real(DP), parameter :: EARTH_ECC2 = 2.0_DP * EARTH_FLATTENING - EARTH_FLATTENING**2
        
        real(DP) :: sin_lat, cos_lat, sin_lon, cos_lon, N
        
        sin_lat = sin(lat * DEG_TO_RAD)
        cos_lat = cos(lat * DEG_TO_RAD)
        sin_lon = sin(lon * DEG_TO_RAD)
        cos_lon = cos(lon * DEG_TO_RAD)
        
        N = EARTH_RADIUS / sqrt(1.0_DP - EARTH_ECC2 * sin_lat**2)
        
        eci_pos(1) = (N + alt) * cos_lat * cos_lon
        eci_pos(2) = (N + alt) * cos_lat * sin_lon
        eci_pos(3) = (N * (1.0_DP - EARTH_ECC2) + alt) * sin_lat
    end subroutine geodetic_to_eci

    !> 核心功能：基于 SPICE 的通用参考系状态转换
    !> @param state_in 输入的状态对象 (必须有效且包含历元)
    !> @param target_frame_name 目标参考系名称 (如 'ITRF93', 'J2000')
    !> @param state_out 输出的转换后状态对象
    subroutine transform_frame_state(state_in, target_frame_name, state_out)
        class(pod_frame_state), intent(in) :: state_in
        character(len=*), intent(in) :: target_frame_name
        class(pod_frame_state), intent(inout) :: state_out
        
        real(DP), dimension(6,6) :: xform
        real(DP), dimension(6) :: vec_in, vec_out
        character(len=MAX_STRING_LEN) :: source_frame
        
        ! 1. 安全检查：如果输入状态无效，直接返回无效的输出
        if (.not. state_in%is_valid()) then
            call state_out%clear_state()
            return
        end if
        
        source_frame = state_in%get_frame_name()
        
        ! 2. 性能优化：如果源参考系和目标参考系完全一样，直接复制对象并返回
        if (trim(source_frame) == trim(target_frame_name)) then
            call state_out%copy_state(state_in)
            return
        end if
        
        ! 3. 调用 SPICE 核心计算 6x6 状态转移矩阵 (包含自转对速度的科氏力等修正)
        call sxform(trim(source_frame), trim(target_frame_name), &
                    state_in%get_epoch(), xform)
        
        ! 4. 组装 6 维向量并执行矩阵乘法
        vec_in(1:3) = state_in%get_position()
        vec_in(4:6) = state_in%get_velocity()
        
        vec_out = matmul(xform, vec_in)
        
        ! 5. 将计算结果组装到新的目标对象中
        call state_out%set_position(vec_out(1:3))
        call state_out%set_velocity(vec_out(4:6))
        call state_out%set_epoch(state_in%get_epoch())
        call state_out%set_frame_name(trim(target_frame_name))
        call state_out%set_valid(.true.)
        
    end subroutine transform_frame_state

    !> 高级功能：将任意坐标系下的状态，安全地转换为地球经纬度高程
    !> 内部会自动先将其转换到地固系 (ITRF93) 以确保地球自转被正确处理
    subroutine state_to_geodetic(state_in, lat, lon, alt)
        class(pod_frame_state), intent(in) :: state_in
        real(DP), intent(out) :: lat, lon, alt
        
        type(pod_frame_state) :: state_itrs
        
        ! 1. 强制转换到高精度地固系 (ITRF93)
        call transform_frame_state(state_in, 'ITRF93', state_itrs)
        
        ! 2. 确保转换成功
        if (.not. state_itrs%is_valid()) then
            lat = 0.0_DP
            lon = 0.0_DP
            alt = 0.0_DP
            return
        end if
        
        ! 3. 调用 eci_to_geodetic (此时传入的已经是安全的 ITRS 坐标了)
        call eci_to_geodetic(state_itrs%get_position(), lat, lon, alt)
        
    end subroutine state_to_geodetic

end module pod_frame_module
