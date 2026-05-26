
module pod_da_force_model_module
    use pod_global, only: DP
    use pod_config, only: config
    use pod_spice, only: get_body_state, pxform, bodvrd, bodvcd
    use pod_gravity_model_module, only: gravity_field
    use pod_dace_classes
    implicit none

    real(DP), public :: current_epoch0 = 0.0_DP
    public :: set_propagation_epoch, cleanup_gravity_network

    ! =========================================================
    ! N 体常量定义
    ! =========================================================
    integer, parameter :: MAX_BODIES = 11
    ! 所有行星统一使用 BARYCENTER 名称，以确保与 DE440 内核中的天体 ID 和命名完全对齐
    character(len=20), dimension(MAX_BODIES) :: body_names = &
        [character(len=20) :: 'MERCURY BARYCENTER', 'VENUS BARYCENTER', 'EARTH', 'MARS BARYCENTER', &
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

    ! =========================================================
    ! 临时变量池：集中管理 DA 运算所需的临时变量，
    ! 避免在循环中反复申请/释放 DA 句柄导致内存泄露。
    ! =========================================================
    type :: ForceModelTempPool
        ! --- DA 标量临时变量 ---
        type(DA) :: r_mag, r_rel_mag  ! 模长计算临时变量
        type(DA) :: tmp_da1, tmp_da2, tmp_da3, tmp_da4, tmp_da5
        type(DA) :: tmp_da6, tmp_da7, tmp_da8
        type(DA) :: v2, r_dot_v       ! 相对论效应临时变量
        
        ! --- AlgebraicVector 临时变量 ---
        type(AlgebraicVector) :: tmp_vec1, tmp_vec2, tmp_vec3, tmp_vec4
        type(AlgebraicVector) :: tmp_vec5, tmp_vec6, tmp_vec7
        type(AlgebraicVector) :: r_rel, v_rel  ! 相对论效应的向量临时变量
        type(AlgebraicVector) :: relative_pos, solar_direction  ! SRP 临时变量
    contains
        procedure :: init => force_pool_init
        procedure :: destroy => force_pool_destroy
    end type ForceModelTempPool

contains

    ! =========================================================
    ! 临时变量池初始化
    ! =========================================================
    subroutine force_pool_init(this, n)
        class(ForceModelTempPool), intent(inout) :: this
        integer, intent(in) :: n
        call this%r_mag%init(); call this%r_rel_mag%init()
        call this%tmp_da1%init();  call this%tmp_da2%init();  call this%tmp_da3%init()
        call this%tmp_da4%init();  call this%tmp_da5%init();  call this%tmp_da6%init()
        call this%tmp_da7%init();  call this%tmp_da8%init()
        call this%v2%init(); call this%r_dot_v%init()
        call this%tmp_vec1%init(n); call this%tmp_vec2%init(n)
        call this%tmp_vec3%init(n); call this%tmp_vec4%init(n)
        call this%tmp_vec5%init(n); call this%tmp_vec6%init(n)
        call this%tmp_vec7%init(n)
        call this%r_rel%init(n); call this%v_rel%init(n)
        call this%relative_pos%init(n); call this%solar_direction%init(n)
    end subroutine force_pool_init

    ! =========================================================
    ! 临时变量池销毁
    ! =========================================================
    subroutine force_pool_destroy(this)
        class(ForceModelTempPool), intent(inout) :: this
        call this%r_mag%destroy(); call this%r_rel_mag%destroy()
        call this%tmp_da1%destroy();  call this%tmp_da2%destroy();  call this%tmp_da3%destroy()
        call this%tmp_da4%destroy();  call this%tmp_da5%destroy();  call this%tmp_da6%destroy()
        call this%tmp_da7%destroy();  call this%tmp_da8%destroy()
        call this%v2%destroy(); call this%r_dot_v%destroy()
        call this%tmp_vec1%destroy(); call this%tmp_vec2%destroy()
        call this%tmp_vec3%destroy(); call this%tmp_vec4%destroy()
        call this%tmp_vec5%destroy(); call this%tmp_vec6%destroy()
        call this%tmp_vec7%destroy()
        call this%r_rel%destroy(); call this%v_rel%destroy()
        call this%relative_pos%destroy(); call this%solar_direction%destroy()
    end subroutine force_pool_destroy

    ! 设置基准历元的接口 (顶层调用)
    subroutine set_propagation_epoch(epoch)
        real(DP), intent(in) :: epoch
        current_epoch0 = epoch
    end subroutine set_propagation_epoch
    
    !> 计算总加速度的主函数
    subroutine da_compute_acceleration(position, velocity, time, acceleration)
        type(AlgebraicVector), intent(in) :: position, velocity
        real(DP), intent(in) :: time  ! 这里的 time 必须是绝对 TDB 秒数
        type(AlgebraicVector), intent(inout) :: acceleration
        
        type(AlgebraicVector) :: acc_grav, acc_drag, acc_srp, acc_rel
        type(ForceModelTempPool) :: pool

        ! 初始化AlgebraicVector
        if (acc_grav%size /= 3) call acc_grav%init(3)
        if (acc_drag%size /= 3) call acc_drag%init(3)
        if (acc_srp%size /= 3) call acc_srp%init(3)
        if (acc_rel%size /= 3) call acc_rel%init(3)
        if (acceleration%size /= 3) call acceleration%init(3)
        
        ! 初始化临时变量池
        call pool%init(3)
        
        ! 安全检查
        if (.not. is_gravity_network_loaded) call init_gravity_network()
        
        ! 1. 计算多体统一引力 (包含中心点质量、第三体点质量、地月高阶非球形)
        call da_compute_gravity_network(position, time, acc_grav, pool)
        
        ! 2. 算非引力项 (后续完善)
        ! 大气阻力加速度
        if (config%use_drag) then
            call da_compute_atmospheric_drag(position, velocity, acc_drag, pool)
        else
            call vec_mul(0.0_DP, acc_drag, acc_drag)
        end if

        ! 2.2 太阳辐射压 (SRP)
        if (config%use_srp) then
            call da_compute_solar_radiation_pressure(position, time, acc_srp, pool)
        else
            call vec_mul(0.0_DP, acc_srp, acc_srp)
        end if

        ! 相对论效应 (新增)
        if (config%use_relativity) then
            call da_compute_post_newtonian(position, velocity, time, acc_rel, pool)
        else
            call vec_mul(0.0_DP, acc_rel, acc_rel)
        end if
        
        ! 3. 总和: acceleration = acc_grav + acc_drag + acc_srp + acc_rel
        ! 使用显式子程序版本避免临时句柄泄漏
        call vec_add(acc_grav, acc_drag, pool%tmp_vec1)
        call vec_add(pool%tmp_vec1, acc_srp, pool%tmp_vec2)
        call vec_add(pool%tmp_vec2, acc_rel, acceleration)
        
        ! 销毁临时变量池
        call pool%destroy()
        
        ! 销毁局部 DA 变量
        call acc_grav%destroy()
        call acc_drag%destroy()
        call acc_srp%destroy()
        call acc_rel%destroy()

    end subroutine da_compute_acceleration

    !> 核心：多体统一引力网计算 (对标深空架构)
    subroutine da_compute_gravity_network(position, time, acc_grav, pool)
        type(AlgebraicVector), intent(in) :: position
        real(DP), intent(in) :: time
        type(AlgebraicVector), intent(inout) :: acc_grav
        type(ForceModelTempPool), intent(inout) :: pool
        
        integer :: i
        real(DP) :: r_body_mag
        real(DP), dimension(3) :: body_pos, body_vel
        real(DP), dimension(3,3) :: rot_to_body, rot_to_inertial

        ! acc_grav = 0.0_DP
        call vec_mul(0.0_DP, acc_grav, acc_grav)
        
        do i = 1, MAX_BODIES
            if (.not. config%use_planet(i)) cycle
            
            ! ==========================================
            ! 情形 A: 中心天体 (地球, i=3)
            ! ==========================================
            if (i == 3) then
                ! r_mag = position%norm2()  (使用 subroutine 版本避免临时句柄泄漏)
                call vector_norm2_sub(position, pool%r_mag)
                
                ! A.1 中心点质量引力: acc_grav = acc_grav - gm_planets(i) * position / r_mag**3
                ! 先计算 r_mag**3
                call da_mul(pool%r_mag, pool%r_mag, pool%tmp_da1)
                call da_mul(pool%tmp_da1, pool%r_mag, pool%tmp_da2)
                ! tmp_da2 = r_mag**3
                ! 计算 position / r_mag**3
                call vec_div(position, pool%tmp_da2, pool%tmp_vec1)
                ! 乘以 gm_planets(i)
                call vec_mul(gm_planets(i), pool%tmp_vec1, pool%tmp_vec2)
                ! acc_grav = acc_grav - tmp_vec2
                call vec_sub(acc_grav, pool%tmp_vec2, pool%tmp_vec3)
                acc_grav = pool%tmp_vec3
                
               ! A.2 地球高阶非球形引力
                if (config%use_earth_nspheric) then
                    call pxform('J2000', 'IAU_EARTH', time, rot_to_body)
                    call vec_matmul(rot_to_body, position, earth_grav%dr_da)
                    
                    ! 1. 用临时向量接收地固系下的带谐和田谐加速度
                    call earth_grav%f_zonal_da(pool%tmp_vec6)
                    call earth_grav%f_tesseral_da(pool%tmp_vec7)
            
                    ! 2. 在地固系下相加：tmp_vec6 = acc_z + acc_t
                    call vec_add(pool%tmp_vec6, pool%tmp_vec7, pool%tmp_vec6)

                    ! 3. 转回 J2000 惯性系：tmp_vec7 = matmul(rot_to_inertial, tmp_vec6)
                    rot_to_inertial = transpose(rot_to_body)
                    call vec_matmul(rot_to_inertial, pool%tmp_vec6, pool%tmp_vec7)
                    
                    ! 4. 安全地累加到总加速度中
                    call vec_add(acc_grav, pool%tmp_vec7, pool%tmp_vec6)
                    acc_grav = pool%tmp_vec6

                    call earth_grav%dr_da%destroy()  ! 释放 3 个 DA 句柄
                end if
                
            ! ==========================================
            ! 情形 B: 第三体摄动 (如月球、太阳、木星)
            ! ==========================================
            else
                if (.not. config%use_third_body) cycle
                call get_body_state(body_names(i), time, 'EARTH', body_pos, body_vel)

                ! 向量 r_rel = 卫星位置 - 天体位置 (从天体指向卫星)
                ! r_rel = position - body_pos
                call vec_sub(position, body_pos, pool%r_rel)
                r_body_mag = norm2(body_pos)
                ! r_rel_mag = r_rel%norm2()
                call vector_norm2_sub(pool%r_rel, pool%r_rel_mag)

                ! B.1 第三体点质量引力 (直接项 + 间接项)
                ! acc_grav = acc_grav - gm_planets(i) * (r_rel / r_rel_mag**3 + body_pos / r_body_mag**3)
                ! 计算 r_rel_mag**3
                call da_mul(pool%r_rel_mag, pool%r_rel_mag, pool%tmp_da1)
                call da_mul(pool%tmp_da1, pool%r_rel_mag, pool%tmp_da2)
                ! tmp_da2 = r_rel_mag**3
                ! 计算 r_rel / r_rel_mag**3
                call vec_div(pool%r_rel, pool%tmp_da2, pool%tmp_vec1)
                ! 计算 body_pos / r_body_mag**3 (实数除法)
                ! 注意：body_pos / r_body_mag**3 是实数数组，使用 vec_add 将实数数组加到 DA 向量上
                call vec_add(pool%tmp_vec1, body_pos / r_body_mag**3, pool%tmp_vec3)
                ! 乘以 gm_planets(i)
                call vec_mul(gm_planets(i), pool%tmp_vec3, pool%tmp_vec4)
                ! acc_grav = acc_grav - tmp_vec4
                call vec_sub(acc_grav, pool%tmp_vec4, pool%tmp_vec5)
                acc_grav = pool%tmp_vec5

                ! B.2 月球高阶非球形引力
              ! B.2 月球高阶非球形引力
                if (i == 10 .and. config%use_moon_nspheric) then
                    ! 对于月球，必须转到月固系 (IAU_MOON 或 MOON_PA)
                    call pxform('J2000', 'IAU_MOON', time, rot_to_body)
                    call vec_matmul(rot_to_body, pool%r_rel, moon_grav%dr_da)

                    ! 1. 用临时向量接收月固系下的带谐和田谐加速度
                    call moon_grav%f_zonal_da(pool%tmp_vec6)
                    call moon_grav%f_tesseral_da(pool%tmp_vec7)

                    ! 2. 在月固系下相加：tmp_vec6 = acc_z + acc_t
                    call vec_add(pool%tmp_vec6, pool%tmp_vec7, pool%tmp_vec6)

                    ! 3. 转回 J2000 惯性系：tmp_vec7 = matmul(rot_to_inertial, tmp_vec6)
                    rot_to_inertial = transpose(rot_to_body)
                    call vec_matmul(rot_to_inertial, pool%tmp_vec6, pool%tmp_vec7)

                    ! 4. 安全地累加到总加速度中
                    call vec_add(acc_grav, pool%tmp_vec7, pool%tmp_vec6)
                    acc_grav = pool%tmp_vec6
                    call moon_grav%dr_da%destroy()  ! 释放 3 个 DA 句柄
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

    !> 清理多体引力网与高阶重力场释放内存
    subroutine cleanup_gravity_network()
        if (is_gravity_network_loaded) then
            if (earth_grav%dr_da%size > 0) call earth_grav%dr_da%destroy()
            if (moon_grav%dr_da%size > 0) call moon_grav%dr_da%destroy()
            is_gravity_network_loaded = .false.
            ! write(*,*) '>>> 多体统一引力网内存已安全释放。'
        end if
    end subroutine cleanup_gravity_network

    subroutine da_compute_atmospheric_drag(position, velocity, acceleration, pool)
        type(AlgebraicVector), intent(in) :: position, velocity
        type(AlgebraicVector), intent(inout) :: acceleration
        type(ForceModelTempPool), intent(inout) :: pool
        
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
        
        ! 阻力加速度: acceleration = drag_factor * velocity
        call vec_mul(drag_factor, velocity, acceleration)
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
    subroutine da_compute_solar_radiation_pressure(position, time, acceleration, pool, Cr, SMR, RP)
        type(AlgebraicVector), intent(in) :: position
        real(DP), intent(in) :: time
        type(AlgebraicVector), intent(inout) :: acceleration
        type(ForceModelTempPool), intent(inout) :: pool
        real(DP), intent(in), optional :: Cr, SMR, RP   ! 新增可选参数
        
        real(DP) :: solar_distance, reflectivity, area_mass_ratio, nominal_rp
        real(DP), dimension(3) :: sun_position, sun_velocity

        ! 默认值 (与 f_SRP 对齐)
        reflectivity = 1.25_DP
        area_mass_ratio = 7.5e-3_DP
        nominal_rp = SOLAR_CONSTANT / SPEED_OF_LIGHT   ! ≈ 4.56e-6 N/m²
        
        if (present(Cr)) reflectivity = Cr
        if (present(SMR)) area_mass_ratio = SMR
        if (present(RP))  nominal_rp = RP
        
        ! 动态获取太阳位置
        call get_body_state('SUN', time, 'EARTH', sun_position, sun_velocity)
        ! relative_pos = position - sun_position
        call vec_sub(position, sun_position, pool%relative_pos)
        solar_distance = norm2(pool%relative_pos%cons())
        ! solar_direction = relative_pos / solar_distance
        call vec_div(pool%relative_pos, solar_distance, pool%solar_direction)
        
        ! 核心炮弹球模型公式
        ! acceleration = reflectivity * area_mass_ratio * nominal_rp * (AU_KM / solar_distance)**2 * solar_direction
        call vec_mul(reflectivity * area_mass_ratio * nominal_rp * (AU_KM / solar_distance)**2, &
                     pool%solar_direction, acceleration)
        ! 单位转换 m/s² -> km/s²
        call vec_mul(1.0e-3_DP, acceleration, acceleration)
        
        ! 如需阴影模型，可在此调用 compute_illumination_factor (已实现，当前未激活)
    end subroutine da_compute_solar_radiation_pressure

    ! ======================================================================
    ! 后牛顿相对论效应 (一级)，只考虑日、地、月 (DA 版本)
    ! ======================================================================
    subroutine da_compute_post_newtonian(position, velocity, time, acc_rel, pool)
        type(AlgebraicVector), intent(in) :: position, velocity
        real(DP), intent(in) :: time
        type(AlgebraicVector), intent(inout) :: acc_rel
        type(ForceModelTempPool), intent(inout) :: pool
        
        integer :: i
        integer, parameter :: rel_bodies(3) = [3, 10, 11]   ! 地球、月球、太阳
        real(DP), dimension(3) :: body_pos, body_vel

        ! acc_rel = 0.0_DP
        call vec_mul(0.0_DP, acc_rel, acc_rel)
        
        do i = 1, 3
            call get_body_state(body_names(rel_bodies(i)), time, 'EARTH', body_pos, body_vel)
            ! r_rel = position - body_pos      ! DA - DP
            call vec_sub(position, body_pos, pool%r_rel)
            ! v_rel = velocity - body_vel
            call vec_sub(velocity, body_vel, pool%v_rel)
            ! r = r_rel%norm2()                ! DA 模长
            call vector_norm2_sub(pool%r_rel, pool%r_mag)
            if (pool%r_mag%cons() < 1.0e-12_DP) cycle ! 避免奇异
            
            ! v2 = v_rel*v_rel   ! DA 标量 (点积)
            call vector_dot_vector_sub(pool%v_rel, pool%v_rel, pool%v2)
            ! r_dot_v = r_rel*v_rel / r   ! (DA⋅DA)/DA = DA
            call vector_dot_vector_sub(pool%r_rel, pool%v_rel, pool%tmp_da1)
            call da_div(pool%tmp_da1, pool%r_mag, pool%r_dot_v)
            
            ! 一阶后牛顿加速度 (与 f_PNE 完全一致)
            ! acc_rel = acc_rel + gm_planets(rel_bodies(i)) / (C_LIGHT_KM**2 * r**2) * (
            !          (4.0_DP * gm_planets(rel_bodies(i)) / r - v2) * (r_rel / r) +
            !          4.0_DP * r_dot_v * v_rel )
            
            ! 计算公共因子: gm / (C_LIGHT_KM**2 * r**2)
            call da_mul(pool%r_mag, pool%r_mag, pool%tmp_da1)  ! r**2
            call da_mul(pool%tmp_da1, C_LIGHT_KM**2, pool%tmp_da2)  ! C_LIGHT_KM**2 * r**2
            ! da_div 不支持 real/DA，直接调用 real_div_da_sub
            call real_div_da_sub(gm_planets(rel_bodies(i)), pool%tmp_da2, pool%tmp_da3)  ! gm / (C_LIGHT_KM**2 * r**2)
            
            ! 计算第一项: (4.0_DP * gm / r - v2) * (r_rel / r)
            ! da_div 不支持 real/DA，直接调用 real_div_da_sub
            call real_div_da_sub(4.0_DP * gm_planets(rel_bodies(i)), pool%r_mag, pool%tmp_da4)
            call da_sub(pool%tmp_da4, pool%v2, pool%tmp_da5)  ! (4*gm/r - v2)
            call vec_div(pool%r_rel, pool%r_mag, pool%tmp_vec1)  ! r_rel / r
            call vec_mul(pool%tmp_da5, pool%tmp_vec1, pool%tmp_vec2)  ! (4*gm/r - v2) * (r_rel/r)
            
            ! 计算第二项: 4.0_DP * r_dot_v * v_rel
            call da_mul(pool%r_dot_v, 4.0_DP, pool%tmp_da6)
            call vec_mul(pool%tmp_da6, pool%v_rel, pool%tmp_vec3)
            
            ! 相加: tmp_vec2 + tmp_vec3
            call vec_add(pool%tmp_vec2, pool%tmp_vec3, pool%tmp_vec4)
            
            ! 乘以公共因子
            call vec_mul(pool%tmp_da3, pool%tmp_vec4, pool%tmp_vec5)
            
            ! acc_rel = acc_rel + tmp_vec5
            call vec_add(acc_rel, pool%tmp_vec5, pool%tmp_vec6)
            acc_rel = pool%tmp_vec6
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
   
