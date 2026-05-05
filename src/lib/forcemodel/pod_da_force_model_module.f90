
module pod_da_force_model_module
    use pod_global, only: DP
    use pod_config, only: config
    use pod_spice, only: get_body_state, pxform, bodvrd, bodvcd
    use pod_gravity_model_module, only: gravity_field
    use pod_dace_classes
    implicit none

    real(DP), public :: current_epoch0 = 0.0_DP
    
    public :: set_propagation_epoch

    ! =========================================================
    ! N 体常量定义
    ! =========================================================
    integer, parameter :: MAX_BODIES = 11
    ! 将 JUPITER 改为 JUPITER BARYCENTER，将其它外行星也做类似处理，以确保与 DE440 内核中的天体 ID 和命名完全对齐
    character(len=20), dimension(MAX_BODIES) :: body_names = &
        [character(len=20) :: 'MERCURY', 'VENUS', 'EARTH', 'MARS', &
         'JUPITER BARYCENTER', 'SATURN BARYCENTER', 'URANUS BARYCENTER', &
         'NEPTUNE BARYCENTER', 'PLUTO BARYCENTER', 'MOON', 'SUN']
    real(DP), dimension(MAX_BODIES) :: gm_planets

    ! 地球与月球的独立高阶引力场对象
    type(gravity_field) :: earth_grav
    type(gravity_field) :: moon_grav
    logical :: is_gravity_network_loaded = .false.

    ! 地球物理常数
    real(DP), parameter :: EARTH_RADIUS = 6378.137_DP  ! km
    ! 太阳常数以及天文单位
    real(DP), parameter :: SOLAR_CONSTANT = 1367.0_DP       ! W/m^2
    real(DP), parameter :: SPEED_OF_LIGHT = 299792458.0_DP  ! m/s
    real(DP), parameter :: AU_KM = 149597870.7_DP           ! 1 AU (km)
    ! 新增：用于相对论效应的光速 (km/s)
    real(DP), parameter :: C_LIGHT_KM = SPEED_OF_LIGHT * 1.0e-3_DP   ! 299792.458 km/s


contains
    ! 设置基准历元的接口 (顶层调用)
    subroutine set_propagation_epoch(epoch)
        real(DP), intent(in) :: epoch
        current_epoch0 = epoch
    end subroutine set_propagation_epoch
        !> 计算总加速度的主函数
    subroutine da_compute_acceleration(position, velocity, time, acceleration)
        type(AlgebraicVector), intent(in) :: position, velocity
        real(DP), intent(in) :: time  ! 这里的 time 必须是绝对 TDB 秒数
        type(AlgebraicVector), intent(out) :: acceleration
        
        type(AlgebraicVector) :: acc_grav, acc_drag, acc_srp, acc_rel

        ! 初始化AlgebraicVector
        if (acc_grav%size /= 3) call acc_grav%init(3)
        if (acc_drag%size /= 3) call acc_drag%init(3)
        if (acc_srp%size /= 3) call acc_srp%init(3)
        if (acc_rel%size /= 3) call acc_rel%init(3)
        if (acceleration%size /= 3) call acceleration%init(3)
        
        ! 安全检查
        if (.not. is_gravity_network_loaded) call init_gravity_network()
        
        ! 1. 计算多体统一引力 (包含中心点质量、第三体点质量、地月高阶非球形)
        call da_compute_gravity_network(position, time, acc_grav)
        
        ! 2. 算非引力项 (后续完善)
                ! 大气阻力加速度
        if (config%use_drag) then
            call da_compute_atmospheric_drag(position, velocity, acc_drag)
        else
            acc_drag = 0.0_DP
        end if

        ! 2.2 太阳辐射压 (SRP)
        if (config%use_srp) then
            call da_compute_solar_radiation_pressure(position, time, acc_srp)
        else
            acc_srp = 0.0_DP
        end if

        ! 相对论效应 (新增)
        if (config%use_relativity) then
            call da_compute_post_newtonian(position, velocity, time, acc_rel)
        else
            acc_rel = 0.0_DP
        end if
        
        ! 3. 总和
        acceleration = acc_grav + acc_drag + acc_srp + acc_rel
    end subroutine da_compute_acceleration

    !> 核心：多体统一引力网计算 (对标深空架构)
    subroutine da_compute_gravity_network(position, time, acc_grav)
        type(AlgebraicVector), intent(in) :: position
        real(DP), intent(in) :: time
        type(AlgebraicVector), intent(inout) :: acc_grav
        
        integer :: i
        real(DP) :: r_body_mag
        type(DA) :: r_mag, r_rel_mag
        real(DP), dimension(3) :: body_pos, body_vel
        type(AlgebraicVector) ::  r_rel
        type(AlgebraicVector) :: acc_z, acc_t
        real(DP), dimension(3,3) :: rot_to_body, rot_to_inertial

        acc_grav = 0.0_DP

        ! 必须分配内存，否则传给底层的引力场计算必崩
        call acc_z%init(3)
        call acc_t%init(3)
        
        do i = 1, MAX_BODIES
            if (.not. config%use_planet(i)) cycle
            
            ! ==========================================
            ! 情形 A: 中心天体 (地球, i=3)
            ! ==========================================
            if (i == 3) then
                r_mag = position%norm2()  ! 计算卫星位置的模长 (DA 类型)
                
                ! A.1 中心点质量引力
                acc_grav = acc_grav - gm_planets(i) * position / r_mag**3
                ! write(*,*) '1st step', time
                ! write(*,*) '>>> 地球 ', gm_planets(i), ' km^3/s^2; 卫星位置模长 ', r_mag, ' km'
                ! write(*,*) '>>> 此时的引力为 ', acc_grav%cons()
                ! A.2 地球高阶非球形引力
                if (config%use_earth_nspheric) then
                    call pxform('J2000', 'IAU_EARTH', time, rot_to_body)
                    earth_grav%dr_da = matmul(rot_to_body, position)

                    
                    call earth_grav%f_zonal_da(acc_z)
                    call earth_grav%f_tesseral_da(acc_t)

                    ! write(*,*) '>>> 地球非球形引力分量 (地固系), 带谐: ', acc_z%elements(1), acc_z%elements(2), &
                    ! acc_z%elements(3)
                    ! write(*,*) '>>> 地球非球形引力分量 (地固系), 田谐: ', acc_t%elements(1), acc_t%elements(2), &
                    ! acc_t%elements(3)
                    
                    rot_to_inertial = transpose(rot_to_body)
                    acc_grav = acc_grav + matmul(rot_to_inertial, acc_z + acc_t)
                end if
                
            ! ==========================================
            ! 情形 B: 第三体摄动 (如月球、太阳、木星)
            ! ==========================================
            else
                if (.not. config%use_third_body) cycle
                call get_body_state(body_names(i), time, 'EARTH', body_pos, body_vel)

                ! 向量 r_rel = 卫星位置 - 天体位置 (从天体指向卫星)
                r_rel = position - body_pos
                r_body_mag = norm2(body_pos)
                r_rel_mag  = r_rel%norm2()

                ! write(*,*) '>>> 天体位置为: ', r_body_mag, ' km; 卫星相对位置为: ', r_rel_mag, ' km'
                
                ! B.1 第三体点质量引力 (直接项 + 间接项)
                acc_grav = acc_grav - gm_planets(i) * (r_rel / r_rel_mag**3 + body_pos / r_body_mag**3)

                ! B.2 月球高阶非球形引力
                if (i == 10 .and. config%use_moon_nspheric) then
                    ! 对于月球，必须转到月固系 (IAU_MOON 或 MOON_PA)
                    ! 注意：送给重力模型的位置是卫星相对于月球的位置 (r_rel)
                    call pxform('J2000', 'IAU_MOON', time, rot_to_body)
                    moon_grav%dr_da = matmul(rot_to_body, r_rel)
                    call moon_grav%f_zonal_da(acc_z)
                    call moon_grav%f_tesseral_da(acc_t)
                    rot_to_inertial = transpose(rot_to_body)
                    acc_grav = acc_grav + matmul(rot_to_inertial, acc_z + acc_t)
                end if
                
            end if
        end do
    end subroutine da_compute_gravity_network

        !> 初始化多体引力网与高阶重力场
    subroutine init_gravity_network()
        integer :: i
        ! 抛弃 SPICE bodvcd 提取，硬编码 DE440 标准 GM 值 
        ! 单位: km^3/s^2。这保证了你的力学模型永远和 DE440 物理基准绝对对齐！
        ! 索引: 1=水星, 2=金星, 3=地球, 4=火星, 5=木星, 6=土星, 7=天王星, 8=海王星, 9=冥王星, 10=月球, 11=太阳
        real(DP), dimension(MAX_BODIES) :: de440_gms = [ &
            22032.080486418_DP,     & ! 1. Mercury
            324858.592_DP,          & ! 2. Venus
            398600.435507_DP,       & ! 3. Earth
            42828.375214_DP,        & ! 4. Mars
            126686531.900_DP,       & ! 5. Jupiter
            37931206.234_DP,        & ! 6. Saturn
            5793951.322_DP,         & ! 7. Uranus
            6836527.100580_DP,      & ! 8. Neptune
            975.5014_DP,            & ! 9. Pluto
            4902.80021852638_DP,    & ! 10. Moon 
            132712440041.279419_DP  ] ! 11. Sun
        
        write(*,*) '>>> 正在初始化多体统一引力网...'
        
        ! 直接从标准数组中映射激活天体的 GM 值
        do i = 1, MAX_BODIES
            if (config%use_planet(i)) then
                gm_planets(i) = de440_gms(i)
            end if
        end do
        
        ! 2. 挂载地球高阶场
        if (config%use_earth_nspheric) then
            earth_grav%cen_body = 3
            earth_grav%ncs = config%earth_degree
            call earth_grav%read_gravity_field()
        end if
        
        ! 3. 挂载月球高阶场
        if (config%use_moon_nspheric) then
            moon_grav%cen_body = 10
            moon_grav%ncs = config%moon_degree
            call moon_grav%read_gravity_field()
        end if

        is_gravity_network_loaded = .true.
    end subroutine init_gravity_network

    subroutine da_compute_atmospheric_drag(position, velocity, acceleration)
        type(AlgebraicVector), intent(in) :: position, velocity
        type(AlgebraicVector), intent(inout) :: acceleration
        
        real(DP) :: altitude, density, velocity_mag, drag_coefficient, area_mass_ratio
        real(DP) :: drag_factor

        
        ! 计算高度
        altitude = norm2(position%cons()) - EARTH_RADIUS
        
        ! 计算大气密度 (简化模型)
        call compute_atmospheric_density(altitude, density)
        
        ! 计算相对速度大小
        velocity_mag = norm2(velocity%cons())
        
        ! 阻力系数和面积质量比 (需要根据具体卫星参数设置)
        drag_coefficient = 2.2_DP  ! 典型值
        area_mass_ratio = 0.01_DP  ! m²/kg
        
        ! 阻力因子
        drag_factor = -0.5_DP * density * drag_coefficient * area_mass_ratio * velocity_mag
        
        ! 阻力加速度
        acceleration = drag_factor * velocity
    end subroutine da_compute_atmospheric_drag
    
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

    ! ======================================================================
    ! 太阳辐射压 (SRP) —— 标准炮弹球模型，支持可选参数
    ! ======================================================================
    subroutine da_compute_solar_radiation_pressure(position, time, acceleration, Cr, SMR, RP)
        type(AlgebraicVector), intent(in) :: position
        real(DP), intent(in) :: time
        type(AlgebraicVector), intent(inout) :: acceleration
        real(DP), intent(in), optional :: Cr, SMR, RP   ! 新增可选参数
        
        real(DP) :: solar_distance, reflectivity, area_mass_ratio, nominal_rp
        real(DP), dimension(3) :: sun_position, sun_velocity
        type(AlgebraicVector) :: relative_pos, solar_direction

        call relative_pos%init(3)
        call solar_direction%init(3)
        
        ! 默认值 (与 f_SRP 对齐)
        reflectivity = 1.25_DP
        area_mass_ratio = 7.5e-3_DP
        nominal_rp = SOLAR_CONSTANT / SPEED_OF_LIGHT   ! ≈ 4.56e-6 N/m²
        
        if (present(Cr)) reflectivity = Cr
        if (present(SMR)) area_mass_ratio = SMR
        if (present(RP))  nominal_rp = RP
        
        ! 动态获取太阳位置
        call get_body_state('SUN', time, 'EARTH', sun_position, sun_velocity)
        relative_pos = position - sun_position
        solar_distance = norm2(relative_pos%cons())
        solar_direction = relative_pos / solar_distance
        
        ! 核心炮弹球模型公式
        acceleration = reflectivity * area_mass_ratio * nominal_rp * &
                      (AU_KM / solar_distance)**2 * solar_direction
        ! 单位转换 m/s² -> km/s²
        acceleration = acceleration * 1.0e-3_DP
        
        ! 如需阴影模型，可在此调用 compute_illumination_factor (已实现，当前未激活)
    end subroutine da_compute_solar_radiation_pressure

    ! ======================================================================
    ! 后牛顿相对论效应 (一级)，只考虑日、地、月 (DA 版本)
    ! ======================================================================
    subroutine da_compute_post_newtonian(position, velocity, time, acc_rel)
        type(AlgebraicVector), intent(in) :: position, velocity
        real(DP), intent(in) :: time
        type(AlgebraicVector), intent(inout) :: acc_rel
        
        integer :: i
        integer, parameter :: rel_bodies(3) = [3, 10, 11]   ! 地球、月球、太阳
        real(DP), dimension(3) :: body_pos, body_vel
        type(AlgebraicVector) :: r_rel, v_rel
        type(DA) :: r, v2, r_dot_v

        acc_rel = 0.0_DP
        
        ! 需要临时变量
        call r_rel%init(3)
        call v_rel%init(3)
        
        do i = 1, 3
            call get_body_state(body_names(rel_bodies(i)), time, 'EARTH', body_pos, body_vel)
            r_rel = position - body_pos      ! DA - DP
            v_rel = velocity - body_vel
            r = r_rel%norm2()                ! DA 模长
            if (r%cons() < 1.0e-12_DP) cycle ! 避免奇异
            
            v2 = dot_product(v_rel, v_rel)   ! DA 标量
            r_dot_v = dot_product(r_rel, v_rel) / r   ! (DA⋅DA)/DA = DA
            
            ! 一阶后牛顿加速度 (与 f_PNE 完全一致)
            acc_rel = acc_rel + gm_planets(rel_bodies(i)) / (C_LIGHT_KM**2 * r**2) * ( &
                     (4.0_DP * gm_planets(rel_bodies(i)) / r - v2) * (r_rel / r) + &
                     4.0_DP * r_dot_v * v_rel )
        end do
    end subroutine da_compute_post_newtonian

    ! ======================================================================
    ! 光照因子 (圆锥形阴影) — 保留作为可选扩展，当前未在主 SRP 中调用
    ! ======================================================================
    subroutine compute_illumination_factor(r_sat_wrt_body, r_sun_wrt_body, R_body, nu)
        real(DP), dimension(3), intent(in) :: r_sat_wrt_body
        real(DP), dimension(3), intent(in) :: r_sun_wrt_body
        real(DP), intent(in) :: R_body
        real(DP), intent(out) :: nu
        
        real(DP), parameter :: R_SUN = 695700.0_DP
        real(DP), parameter :: PI = 3.14159265358979323846_DP
        
        real(DP), dimension(3) :: r_sat_to_sun, r_sat_to_body
        real(DP) :: dist_sat_to_sun, dist_sat_to_body
        real(DP) :: a, b, c, cos_c, x, y, area_overlap
        
        r_sat_to_sun = r_sun_wrt_body - r_sat_wrt_body
        r_sat_to_body = -r_sat_wrt_body
        dist_sat_to_sun = norm2(r_sat_to_sun)
        dist_sat_to_body = norm2(r_sat_to_body)
        
        if (dist_sat_to_body <= R_body) then
            nu = 0.0_DP
            return
        end if
        
        a = asin(R_SUN / dist_sat_to_sun)
        b = asin(R_body / dist_sat_to_body)
        cos_c = dot_product(r_sat_to_sun, r_sat_to_body) / (dist_sat_to_sun * dist_sat_to_body)
        cos_c = max(-1.0_DP, min(1.0_DP, cos_c))
        c = acos(cos_c)
        
        if (c >= a + b) then
            nu = 1.0_DP
        else if (c <= b - a) then
            nu = 0.0_DP
        else if (c <= a - b) then
            nu = 1.0_DP - (b / a)**2
        else
            x = (c**2 + a**2 - b**2) / (2.0_DP * c)
            y = sqrt(max(0.0_DP, a**2 - x**2))
            area_overlap = a**2 * acos(x / a) + b**2 * acos((c - x) / b) - c * y
            nu = 1.0_DP - area_overlap / (PI * a**2)
        end if
    end subroutine compute_illumination_factor
    

end module pod_da_force_model_module
   

