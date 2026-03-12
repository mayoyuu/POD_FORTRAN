module pod_force_model_module
    use pod_global, only: DP
    use pod_config, only: config
    use pod_spice, only: get_body_state, pxform, bodvrd
    use pod_gravity_model_module, only: gravity_field
    
    implicit none
    
! =========================================================
    ! N 体常量定义
    ! =========================================================
    integer, parameter :: MAX_BODIES = 11
    character(len=16), dimension(MAX_BODIES) :: body_names = &
        ['MERCURY', 'VENUS  ', 'EARTH  ', 'MARS   ', 'JUPITER', &
         'SATURN ', 'URANUS ', 'NEPTUNE', 'PLUTO  ', 'MOON   ', 'SUN    ']
    real(DP), dimension(MAX_BODIES) :: gm_planets

    ! =========================================================
    ! 力模型统一配置结构体 (Configuration Object)
    ! =========================================================
    type, public :: force_model_config_type
        ! 1. 摄动效应总开关 (Perturbation Switches)
        logical :: use_earth_nspheric = .true.
        logical :: use_moon_nspheric  = .true.
        logical :: use_sun_nspheric   = .true.
        logical :: use_third_body     = .true.
        logical :: use_srp            = .true.
        logical :: use_drag           = .false.
        logical :: use_relativity     = .false.
        
        ! 2. 高阶引力场精度配置 (Gravity Field Degrees)
        integer :: earth_degree = 10
        integer :: moon_degree  = 10
        integer :: sun_degree   = 0   ! 预留，太阳常作点质量处理
        
        ! 3. 多体引力网激活清单 (Active N-Body Network)
        logical, dimension(MAX_BODIES) :: use_planet = .false.
    end type force_model_config_type

    ! 实例化全局力模型配置对象
    type(force_model_config_type), public :: fm_config

    ! 地球与月球的独立高阶引力场对象
    type(gravity_field) :: earth_grav
    type(gravity_field) :: moon_grav
    type(gravity_field) :: sun_grav  ! 预留太阳高阶引力场对象
    logical :: is_gravity_network_loaded = .false.

    ! 地球物理常数
    real(DP), parameter :: EARTH_RADIUS = 6378.137_DP  ! km
    ! 太阳常数以及天文单位
    real(DP), parameter :: SOLAR_CONSTANT = 1367.0_DP       ! W/m^2
    real(DP), parameter :: SPEED_OF_LIGHT = 299792458.0_DP  ! m/s
    real(DP), parameter :: AU_KM = 149597870.7_DP           ! 1 AU (km)

contains

    subroutine init_gravity_network()
        integer :: i, dim
        real(8) :: gm_val(1)
        
        ! 1. 动态加载被激活天体的 GM 值
        do i = 1, MAX_BODIES
            if (fm_config%use_planet(i)) then
                call bodvrd(trim(body_names(i)), 'GM', 1, dim, gm_val)
                gm_planets(i) = gm_val(1)
            end if
        end do
        
        ! 2. 挂载地球高阶场
        if (fm_config%use_earth_nspheric) then
            earth_grav%cen_body = 3
            earth_grav%ncs = fm_config%earth_degree
            call earth_grav%read_gravity_field()
        end if
        
        ! 3. 挂载月球高阶场
        if (fm_config%use_moon_nspheric) then
            moon_grav%cen_body = 10
            moon_grav%ncs = fm_config%moon_degree
            call moon_grav%read_gravity_field()
        end if

        ! 4. 初始化太阳高阶引力场 (通常不考虑，但保留接口)
        if (fm_config%use_planet(11) .and. fm_config%sun_degree > 0) then
            sun_grav%cen_body = 11
            sun_grav%ncs = fm_config%sun_degree
            call sun_grav%read_gravity_field()
        end if
        
        is_gravity_network_loaded = .true.
    end subroutine init_gravity_network

    !> 计算总加速度的主函数
    subroutine compute_acceleration(position, velocity, time, acceleration)
        real(DP), dimension(3), intent(in) :: position, velocity
        real(DP), intent(in) :: time  ! 这里的 time 必须是绝对 TDB 秒数
        real(DP), dimension(3), intent(out) :: acceleration
        
        real(DP), dimension(3) :: acc_grav, acc_drag, acc_srp
        
        ! 安全检查
        if (.not. is_gravity_network_loaded) call init_gravity_network()
        
        ! 1. 计算多体统一引力 (包含中心点质量、第三体点质量、地月高阶非球形)
        call compute_gravity_network(position, time, acc_grav)
        
        ! 2. 算非引力项 (后续完善)
                ! 大气阻力加速度
        if (fm_config%use_drag) then
            call compute_atmospheric_drag(position, velocity, acc_drag)
        else
            acc_drag = 0.0_DP
        end if
        
        ! 太阳辐射压加速度
        if (fm_config%use_srp) then
            ! 【关键修改】把 time 传进去
            call compute_solar_radiation_pressure(position, time, acc_srp)
        else
            acc_srp = 0.0_DP
        end if
        
        ! 3. 总和
        acceleration = acc_grav + acc_drag + acc_srp
    end subroutine compute_acceleration

    !> 核心：多体统一引力网计算 (对标深空架构)
    subroutine compute_gravity_network(position, time, acc_grav)
        real(DP), dimension(3), intent(in) :: position
        real(DP), intent(in) :: time
        real(DP), dimension(3), intent(out) :: acc_grav
        
        integer :: i
        real(DP) :: r_mag, r_body_mag, r_rel_mag
        real(DP), dimension(3) :: body_pos, body_vel, r_rel
        real(DP), dimension(3) :: acc_z, acc_t
        real(DP), dimension(3,3) :: rot_to_body
        
        acc_grav = 0.0_DP
        
        do i = 1, MAX_BODIES
            if (.not. fm_config%use_planet(i)) cycle
            
            ! ==========================================
            ! 情形 A: 中心天体 (地球, i=3)
            ! ==========================================
            if (i == 3) then
                r_mag = norm2(position)
                
                ! A.1 中心点质量引力
                acc_grav = acc_grav - gm_planets(i) * position / r_mag**3
                
                ! A.2 地球高阶非球形引力
                if (fm_config%use_earth_nspheric) then
                    call pxform('J2000', 'IAU_EARTH', time, rot_to_body)
                    earth_grav%dr = matmul(rot_to_body, position)
                    call earth_grav%f_zonal(acc_z)
                    call earth_grav%f_tesseral(acc_t)
                    acc_grav = acc_grav + matmul(transpose(rot_to_body), acc_z + acc_t)
                end if
                
            ! ==========================================
            ! 情形 B: 第三体摄动 (如月球、太阳、木星)
            ! ==========================================
            else
                ! 获取天体相对于地球的位置
                call get_body_state(trim(body_names(i)), time, 'EARTH', body_pos, body_vel)
                
                ! 向量 r_rel = 卫星位置 - 天体位置 (从天体指向卫星)
                r_rel = position - body_pos
                r_body_mag = norm2(body_pos)
                r_rel_mag  = norm2(r_rel)
                
                ! B.1 第三体点质量引力 (直接项 + 间接项)
                acc_grav = acc_grav - gm_planets(i) * (r_rel / r_rel_mag**3 + body_pos / r_body_mag**3)
                
                ! B.2 [Cislunar 杀手锏] 月球高阶非球形引力
                if (i == 10 .and. fm_config%use_moon_nspheric) then
                    ! 对于月球，必须转到月固系 (IAU_MOON 或 MOON_PA)
                    ! 注意：送给重力模型的位置是卫星相对于月球的位置 (r_rel)
                    call pxform('J2000', 'IAU_MOON', time, rot_to_body)
                    moon_grav%dr = matmul(rot_to_body, r_rel)
                    call moon_grav%f_zonal(acc_z)
                    call moon_grav%f_tesseral(acc_t)
                    acc_grav = acc_grav + matmul(transpose(rot_to_body), acc_z + acc_t)
                end if
                ! 3. 太阳高阶引力 (如果启用且阶数 > 0)
                if (i == 11 .and. fm_config%use_sun_nspheric) then
                    call pxform('J2000', 'IAU_SUN', time, rot_to_body)
                    sun_grav%dr = matmul(rot_to_body, r_rel)
                    call sun_grav%f_zonal(acc_z)
                    call sun_grav%f_tesseral(acc_t)
                    acc_grav = acc_grav + matmul(transpose(rot_to_body), acc_z + acc_t)
                end if
                
            end if
        end do
    end subroutine compute_gravity_network

    ! ... (保留原来的 options/status 打印等辅助函数) ... 
    subroutine compute_atmospheric_drag(position, velocity, acceleration)
        real(DP), dimension(3), intent(in) :: position, velocity
        real(DP), dimension(3), intent(out) :: acceleration
        
        real(DP) :: altitude, density, velocity_mag, drag_coefficient, area_mass_ratio
        real(DP) :: drag_factor
        
        ! 计算高度
        altitude = sqrt(sum(position**2)) - EARTH_RADIUS
        
        ! 计算大气密度 (简化模型)
        call compute_atmospheric_density(altitude, density)
        
        ! 计算相对速度大小
        velocity_mag = sqrt(sum(velocity**2))
        
        ! 阻力系数和面积质量比 (需要根据具体卫星参数设置)
        drag_coefficient = 2.2_DP  ! 典型值
        area_mass_ratio = 0.01_DP  ! m²/kg
        
        ! 阻力因子
        drag_factor = -0.5_DP * density * drag_coefficient * area_mass_ratio * velocity_mag
        
        ! 阻力加速度
        acceleration = drag_factor * velocity
    end subroutine compute_atmospheric_drag
    
    subroutine compute_atmospheric_density(altitude, density)
        real(DP), intent(in) :: altitude
        real(DP), intent(out) :: density
        
        real(DP), parameter :: H0 = 8.5_DP  ! km
        real(DP), parameter :: RHO0 = 1.225e-3_DP  ! kg/km³ (海平面密度)
        
        ! 指数大气密度模型
        if (altitude < 100.0_DP) then
            density = RHO0 * exp(-altitude / H0)
        else
            density = 0.0_DP  ! 高空大气密度很小
        end if
    end subroutine compute_atmospheric_density

    !! 计算光压效应
    ! ======================================================================
    ! 计算太阳辐射压 (SRP) 加速度
    ! ======================================================================
    subroutine compute_solar_radiation_pressure(position, time, acceleration)
        real(DP), dimension(3), intent(in) :: position
        real(DP), intent(in) :: time  ! 必须传入时间以获取动态太阳位置
        real(DP), dimension(3), intent(out) :: acceleration
        
        real(DP) :: solar_distance, solar_factor, reflectivity, area_mass_ratio
        real(DP) :: nominal_srp_pressure
        real(DP), dimension(3) :: sun_position, sun_velocity
        real(DP), dimension(3) :: relative_pos, solar_direction
        
        ! 1 AU 的标准距离 (单位: km)
        real(DP), parameter :: AU_KM = 149597870.7_DP 
        
        ! 1. 动态获取当前时刻太阳相对于地球的位置
        call get_body_state('SUN', time, 'EARTH', sun_position, sun_velocity)
        
        ! 2. 计算从太阳指向卫星的相对位置向量
        ! 注意：位置相减 (Sat - Sun) 得到的是背离太阳的方向，正是光压的推力方向
        relative_pos = position - sun_position
        solar_distance = norm2(relative_pos)
        
        ! 提取单位方向向量
        solar_direction = relative_pos / solar_distance
        
        ! 3. 卫星属性参数设置 (可由外部传入，这里暂用默认值)
        reflectivity = 1.25_DP         ! 反射系数 (Cr)
        area_mass_ratio = 7.0e-3_DP    ! 面质比 (m²/kg)
        
        ! 4. 计算 1 AU 处的标准辐射压强 P_0 (单位: N/m² 或 kg/(m·s²))
        ! P_0 = 1367 (W/m²) / 299792458 (m/s) ≈ 4.56e-6 N/m²
        nominal_srp_pressure = SOLAR_CONSTANT / SPEED_OF_LIGHT
        
        ! 5. 核心物理公式计算 (得到加速度大小，单位: m/s²)
        ! 公式: a = P_0 * Cr * (A/M) * (AU / r)^2
        solar_factor = nominal_srp_pressure * reflectivity * area_mass_ratio * &
                      (AU_KM / solar_distance)**2
        
        ! 6. 附加上方向，并进行极其关键的单位转换 (m/s² -> km/s²)
        acceleration = (solar_factor * 1.0e-3_DP) * solar_direction

        ! ... 前面计算 solar_factor 和 solar_direction 的代码保持不变 ...
        
        ! 引入阴影模型
        ! real(DP) :: nu_earth, nu_moon, nu_total
        ! real(DP), dimension(3) :: moon_position, moon_velocity
        
        ! ! 7. 计算地球阴影因子 (传入相对于地球的坐标)
        ! call compute_illumination_factor(position, sun_position, EARTH_RADIUS, nu_earth)
        
        ! ! 8. 计算月球阴影因子 (必须传入相对于月球的坐标)
        ! call get_body_state('MOON', time, 'EARTH', moon_position, moon_velocity)
        ! call compute_illumination_factor(position - moon_position, &
        !                                  sun_position - moon_position, &
        !                                  1737.4_DP, nu_moon) ! 月球半径 1737.4 km
        
        ! ! 9. 综合光照因子 (取两者相乘是业界处理多重遮挡的标准近似做法)
        ! nu_total = nu_earth * nu_moon
        
        ! ! 10. 最终施加带有圆锥形阴影修正的加速度
        ! acceleration = nu_total * (solar_factor * 1.0e-3_DP) * solar_direction
        
        
    end subroutine compute_solar_radiation_pressure

    ! ======================================================================
    ! 计算极其严密的圆锥形光照因子 (视圆面相交法)
    ! 返回值 nu: 0.0 (本影/全食), 1.0 (全照), (0.0, 1.0) (半影/偏食)
    ! ======================================================================
    subroutine compute_illumination_factor(r_sat_wrt_body, r_sun_wrt_body, R_body, nu)
        real(DP), dimension(3), intent(in) :: r_sat_wrt_body  ! 卫星相对于遮挡天体的位置
        real(DP), dimension(3), intent(in) :: r_sun_wrt_body  ! 太阳相对于遮挡天体的位置
        real(DP), intent(in) :: R_body                        ! 遮挡天体的真实物理半径(km)
        real(DP), intent(out) :: nu                           ! 光照因子 (0 到 1)
        
        real(DP), parameter :: R_SUN = 695700.0_DP            ! 太阳真实半径 (km)
        real(DP), parameter :: PI = 3.14159265358979323846_DP
        
        real(DP), dimension(3) :: r_sat_to_sun, r_sat_to_body
        real(DP) :: dist_sat_to_sun, dist_sat_to_body
        real(DP) :: a, b, c, cos_c, x, y, area_overlap
        
        ! 1. 计算卫星到太阳、卫星到遮挡天体的向量和距离
        r_sat_to_sun = r_sun_wrt_body - r_sat_wrt_body
        r_sat_to_body = -r_sat_wrt_body  ! 遮挡天体在原点
        
        dist_sat_to_sun = norm2(r_sat_to_sun)
        dist_sat_to_body = norm2(r_sat_to_body)
        
        ! 极值保护：如果卫星已经砸在天体内部了，直接全黑
        if (dist_sat_to_body <= R_body) then
            nu = 0.0_DP
            return
        end if
        
        ! 2. 计算视半径 (Apparent Angular Radius)
        ! a = 太阳的视半径，b = 遮挡天体的视半径
        a = asin(R_SUN / dist_sat_to_sun)
        b = asin(R_body / dist_sat_to_body)
        
        ! 3. 计算视圆心角距 (Apparent Angular Separation) c
        cos_c = dot_product(r_sat_to_sun, r_sat_to_body) / (dist_sat_to_sun * dist_sat_to_body)
        cos_c = max(-1.0_DP, min(1.0_DP, cos_c)) ! 防止浮点截断导致 acos 越界
        c = acos(cos_c)
        
        ! 4. 核心几何判断：四个极其严密的物理边界
        if (c >= a + b) then
            ! 完全没有遮挡 (全照)
            nu = 1.0_DP
        else if (c <= b - a) then
            ! 天体视圆完全盖住太阳 (全食 / 本影)
            nu = 0.0_DP
        else if (c <= a - b) then
            ! 太阳视圆太大，天体完全在太阳内部 (环食 / 伪本影)
            ! 光照因子 = 1 - (天体视面积 / 太阳视面积)
            nu = 1.0_DP - (b / a)**2
        else
            ! 圆面部分相交 (偏食 / 半影)
            ! 利用余弦定理求解两个扇形的重叠面积
            x = (c**2 + a**2 - b**2) / (2.0_DP * c)
            y = sqrt(max(0.0_DP, a**2 - x**2))
            
            area_overlap = a**2 * acos(x / a) + b**2 * acos((c - x) / b) - c * y
            
            nu = 1.0_DP - area_overlap / (PI * a**2)
        end if
        
    end subroutine compute_illumination_factor
    
    
!> 接口 1：统一设置摄动总开关
    subroutine set_perturbation_switches(earth_grav, moon_grav, third_body, srp, drag, relativity)
        logical, intent(in) :: earth_grav, moon_grav, third_body, srp, drag, relativity
        
        fm_config%use_earth_nspheric = earth_grav
        fm_config%use_moon_nspheric  = moon_grav
        fm_config%use_third_body     = third_body
        fm_config%use_srp            = srp
        fm_config%use_drag           = drag
        fm_config%use_relativity     = relativity
        
        ! 任何开关的改动，都标记需要重新校验引力网
        is_gravity_network_loaded = .false. 
    end subroutine set_perturbation_switches

    !> 接口 2：统一设置高阶重力场截断阶数
    subroutine set_gravity_degrees(earth_deg, moon_deg)
        integer, intent(in) :: earth_deg, moon_deg
        
        fm_config%earth_degree = earth_deg
        fm_config%moon_degree  = moon_deg
        is_gravity_network_loaded = .false. 
    end subroutine set_gravity_degrees

    !> 接口 3：配置多体引力网的激活天体 (传入需要激活的天体编号数组)
    subroutine set_active_planets(active_body_ids)
        integer, dimension(:), intent(in) :: active_body_ids
        integer :: i, body_id
        
        fm_config%use_planet = .false. ! 先全量清空
        
        do i = 1, size(active_body_ids)
            body_id = active_body_ids(i)
            if (body_id >= 1 .and. body_id <= MAX_BODIES) then
                fm_config%use_planet(body_id) = .true.
            end if
        end do
        is_gravity_network_loaded = .false.
    end subroutine set_active_planets

    !> 打印当前力模型状态
    subroutine print_force_model_status()
        integer :: i
        character(len=64) :: status_str
        
        write(*, *) '=================================================='
        write(*, *) '             POD 力模型状态监控面板               '
        write(*, *) '=================================================='
        
        write(*, *) '[1. 引力与多体摄动]'
        write(*, *) '  中心引力 (点质量) : 默认启用'
        write(*, *) '  第三体摄动总开关  : ', merge('启用', '禁用', fm_config%use_third_body)
        
        ! 动态格式化地球高阶场状态
        if (fm_config%use_earth_nspheric) then
            write(status_str, "(A,I3,A)") '启用 (', fm_config%earth_degree, ' 阶)'
        else
            status_str = '禁用'
        end if
        write(*, *) '  地球高阶非球形    : ', trim(status_str)
        
        ! 动态格式化月球高阶场状态
        if (fm_config%use_moon_nspheric) then
            write(status_str, "(A,I3,A)") '启用 (', fm_config%moon_degree, ' 阶)'
        else
            status_str = '禁用'
        end if
        write(*, *) '  月球高阶非球形    : ', trim(status_str)
        
        write(*, *) ''
        write(*, *) '[2. 表面力与其他效应]'
        write(*, *) '  太阳辐射压 (SRP)  : ', merge('启用', '禁用', fm_config%use_srp)
        write(*, *) '  大气阻力 (Drag)   : ', merge('启用', '禁用', fm_config%use_drag)
        write(*, *) '  相对论效应 (PNE)  : ', merge('启用', '禁用', fm_config%use_relativity)
        
        write(*, *) ''
        write(*, *) '[3. 动态引力网激活清单]'
        write(*, "(A)", advance='no') '  包含天体: '
        do i = 1, MAX_BODIES
            if (fm_config%use_planet(i)) then
                write(*, "(A)", advance='no') trim(body_names(i)) // ' '
            end if
        end do
        write(*, *) '' ! 补一个换行
        
        write(*, *) '=================================================='
    end subroutine print_force_model_status

end module pod_force_model_module
