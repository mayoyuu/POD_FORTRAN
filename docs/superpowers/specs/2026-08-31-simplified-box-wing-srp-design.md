# POD_Fortran 简化 Box-Wing 太阳光压模型设计

日期：2026-08-31

状态：设计已确认，等待书面规格审阅

适用范围：POD_Fortran 的 real 与 DA 高保真动力学传播

## 1. 决策摘要

在现有炮弹球模型之外增加简化 box-wing SRP 模型。模型由六个箱体表面和一个代表两侧太阳翼总面积的等效双面平板组成，不计算地影、月影、自遮挡、多次反射、热辐射和光压转矩。

姿态采用混合接口：

1. 默认解析生成对日、对地、对月三种互斥姿态。
2. 提供给定实数姿态矩阵的低层入口，未来 CK 或外部姿态文件可以接入。
3. 解析姿态具有完整 real 和 DA 实现；外部实数姿态可叠加三个 DA 小角度误差。
4. 现有 kernels/ck 中的 Rosetta CK 不装载、不参与当前卫星计算。

炮弹球模型保留且仍为默认模型。只有显式设置 srp_model=box_wing 时才启用新模型。

## 2. 目标与非目标

目标：

- 计算对日、对地、对月三种定姿方式下的 box-wing 光压。
- 箱体面和太阳翼使用可配置光学参数。
- 同一物理模型具有 real 和 DA 两条路径。
- DA 保留状态—视线—姿态—入射角—光压的完整导数链。
- 保留全局 SRP 尺度、姿态小角度和太阳翼转角不确定性接口。
- 与现有 real/DA 力模型和 UQ 传播兼容。
- 形成物理极限、real/DA 一致性、有限差分导数和端到端传播测试闭环。

非目标：不计算进出影、半影、自遮挡、反照、红外、热反冲和 SRP 力矩；不积分姿态动力学；不实现太阳翼机械限位或柔性；不下载或修复任务内核；第一版时间不是 DA 变量。

## 3. 现有 CK 的适用性

- RATT_*.BC（对象 -226000）是 Rosetta 航天器姿态。
- ROS_SA_2015_V0013.BC（对象 -226015/-226025）是 Rosetta 太阳翼姿态，覆盖约 2015-01-01 至 2015-04-10。
- ROS_HGA_*.BC 是 Rosetta 天线姿态；CATT_*.BC 是 67P 彗星本体姿态。
- ckbrief 报告缺少匹配 SCLK；仓库也没有完整 Rosetta FK。
- 当前 pod_spice.f90 没有 CK/SCLK/FK 装载和 SCE2C、CKGPAV 封装。

所以第一版使用解析姿态。未来取得目标卫星匹配的 CK/SCLK/FK 后，通过外部姿态入口接入。CKGPAV 的参考系到本体系矩阵必须转成本文约定的本体系到惯性系矩阵，并以已知轴向算例验证。

## 4. 坐标、方向与单位

- 惯性系记为 I，第一版固定 J2000；卫星本体系记为 B。
- C_I_B 从本体系旋转到惯性系：v_I = C_I_B v_B。
- C_I_B 的三列是 +X_B、+Y_B、+Z_B 在惯性系中的坐标。
- srp_pointing_axis_body 指定对准目标的本体系轴，默认 +Z_B=[0,0,1]^T。
- 卫星到太阳单位方向：

~~~text
e_S = (r_S - r_sc) / norm(r_S - r_sc)
~~~

光子传播方向是 -e_S，所以吸收力远离太阳。

单位：位置 km，速度 km/s，时间 ET/TDB 秒，面积 m²，尺寸 m，质量 kg，光压 N/m²，公开加速度 km/s²。

## 5. 姿态模型

目标 T 按 sun、earth、moon 三种互斥模式选择：

~~~text
u_P = (r_T - r_sc) / norm(r_T - r_sc)
~~~

一根轴对准目标只约束两个自由度。默认用相对中心天体的轨道面法向作为滚转参考：

~~~text
q = r_rel × v_rel
z_I = u_P
y_I = normalize(q - dot(q,z_I) z_I)
x_I = y_I × z_I
~~~

对配置的本体系指向轴和滚转轴作相同正交化，得到 T_B；令 T_I=[x_I,y_I,z_I]：

~~~text
C_I_B = T_I transpose(T_B)
~~~

这允许指向轴不是 +Z_B，并保持右手系和 det(C_I_B)=+1。

退化处理：

- 目标距离小于 srp_geometry_tolerance：返回错误。
- 轨道法向或滚转投影退化：从 inertial_z、inertial_x、inertial_y 中选与主轴最不共线者。
- DA 路径以 DA 常数项决定分支，不宣称单个展开可以跨越分支边界。

外部姿态 API 接受实数 C_I_B。DA 路径可叠加小角度：

~~~text
C_I_B_DA = C_I_B_real exp(skew(delta_theta_DA))
~~~

实现采用 Rodrigues/稳定小角度函数。解析姿态则直接由 DA 状态构造，不冻结名义姿态。

## 6. 几何与太阳翼

箱体尺寸为 Lx、Ly、Lz，生成六个单面表面：

| 面 | 本体系外法向 | 面积 |
|---|---|---:|
| +X_B | [1,0,0]^T | Ly Lz |
| -X_B | [-1,0,0]^T | Ly Lz |
| +Y_B | [0,1,0]^T | Lx Lz |
| -Y_B | [0,-1,0]^T | Lx Lz |
| +Z_B | [0,0,1]^T | Lx Ly |
| -Z_B | [0,0,-1]^T | Lx Ly |

第一版六面共享一组箱体光学系数，但底层面板类型保留逐面参数能力。

两侧物理太阳翼合并成一个面积等于两翼总有效面积的双面板。因为不计算力矩和自遮挡，合并不损失净平动力。

太阳翼支持 fixed 和 single_axis。固定模式使用配置法向；单轴模式绕本体系固定铰链轴 a_h 理想跟踪太阳：

~~~text
e_S_B = transpose(C_I_B) e_S_I
p     = e_S_B - dot(e_S_B,a_h) a_h
n_SA_B = normalize(p)
~~~

法向符号选为正面朝向太阳。若投影退化，使用参考法向在铰链法平面内的投影。翼角偏差通过绕铰链轴的 Rodrigues 旋转加入。模型不含机械限位。

## 7. 面板光压公式

每个表面有吸收率 alpha、镜面反射率 rho_s、漫反射率 rho_d：

~~~text
alpha + rho_s + rho_d = 1
0 <= alpha, rho_s, rho_d <= 1
~~~

面板惯性系法向和入射余弦：

~~~text
n_i_I = C_I_B n_i_B
mu_i  = dot(n_i_I,e_S)
~~~

单面表面仅在 mu_i>0 时受光。双面太阳翼用法向相反的正反两个光学表面表示。距离 d 处的光压：

~~~text
P(d) = P_1AU (AU/d)^2
~~~

第 i 个受光面的力：

~~~text
F_i = -P(d) A_i mu_i [
        (1-rho_s_i) e_S
        + 2 (rho_s_i mu_i + rho_d_i/3) n_i_I
      ]
~~~

总加速度：

~~~text
a_SRP = (1+eta_SRP)/mass * sum(F_i)
~~~

结果由 m/s² 转成 km/s²，不乘阴影因子。物理极限：完全吸收法向入射为 -PA e_S；完全镜面为 -2PA e_S；掠入射和单面背光为零。

## 8. 软件结构与接口

新增 src/lib/forcemodel/pod_spacecraft_geometry_module.f90。主要类型：

~~~fortran
type :: srp_optical_properties_type
    real(DP) :: absorptivity
    real(DP) :: specular_reflectivity
    real(DP) :: diffuse_reflectivity
end type

type :: srp_panel_type
    real(DP) :: area_m2
    real(DP) :: normal_body(3)
    type(srp_optical_properties_type) :: optical
end type

type :: srp_da_parameter_map_type
    integer :: global_scale_index
    integer :: attitude_bias_index(3)
    integer :: array_angle_index
    real(DP) :: global_scale_span
    real(DP) :: attitude_bias_span_rad(3)
    real(DP) :: array_angle_span_rad
end type

type :: da_attitude_type
    type(AlgebraicVector) :: x_body_in_inertial
    type(AlgebraicVector) :: y_body_in_inertial
    type(AlgebraicVector) :: z_body_in_inertial
contains
    procedure :: init
    procedure :: destroy
end type

type :: spacecraft_geometry_type
contains
    procedure :: initialize_from_config
    procedure :: validate
    procedure :: compute_attitude_real
    procedure :: compute_attitude_da
    procedure :: compute_srp_real
    procedure :: compute_srp_da
    procedure :: compute_srp_with_attitude_real
    procedure :: compute_srp_with_attitude_da
end type
~~~

高层解析姿态调用：

~~~fortran
call geometry%compute_srp_real(position, velocity, et, acceleration, status, message)
call geometry%compute_srp_da(position_da, velocity_da, et, da_map, &
                             acceleration_da, status, message)
~~~

外部姿态调用：

~~~fortran
call geometry%compute_srp_with_attitude_real(position, et, c_i_b, &
                                              acceleration, status, message)
call geometry%compute_srp_with_attitude_da(position_da, et, attitude_da, da_map, &
                                            acceleration_da, status, message)
~~~

所有运行时状态通过参数传入。几何对象初始化后只读；DA 临时池为调用局部对象，不在全局保存位置、速度或 DA 句柄。pod_engine_init 读取配置后显式初始化默认几何，real 和 DA 力模型共享它；测试可构造局部对象。

## 9. 配置参数

~~~text
srp_model                       = cannonball | box_wing
srp_attitude_mode               = sun | earth | moon
srp_roll_reference             = orbit_normal | inertial_z | inertial_x | sun | earth | moon
srp_pointing_axis_body          = x y z
srp_roll_axis_body              = x y z
srp_mass_kg                     = positive
srp_box_dimensions_m            = Lx Ly Lz
srp_box_optical                 = alpha rho_specular rho_diffuse
srp_array_total_area_m2         = nonnegative
srp_array_tracking_mode         = fixed | single_axis
srp_array_hinge_axis_body       = x y z
srp_array_reference_normal_body = x y z
srp_array_front_optical         = alpha rho_specular rho_diffuse
srp_array_back_optical          = alpha rho_specular rho_diffuse
srp_pressure_1au_n_m2           = positive or negative-to-use-existing-default
srp_geometry_tolerance          = positive
srp_scale_da_span               = nonnegative
srp_attitude_bias_span_arcsec   = sx sy sz
srp_array_angle_span_deg        = nonnegative
~~~

默认保持 srp_model=cannonball；指向轴 [0,0,1]、滚转轴 [0,1,0]、滚转参考 orbit_normal、太阳翼 single_axis。质量、尺寸和光学系数不使用虚构卫星默认值，启用 box-wing 时必须显式给出。太阳翼面积可为零。DA 变量编号不进配置，只存物理尺度。

校验覆盖质量、尺寸、面积、轴向非零且不共线、铰链与参考法向不平行、光学系数范围及和为 1、枚举值和容差。炮弹球不要求 box-wing 参数；box-wing 缺参数不得回退。

## 10. DA 设计

解析姿态 DA 路径中，以下量必须保持为 DA：卫星位置/速度；相对目标向量与单位视线；轨道法向、滚转基和姿态三轴；太阳翼跟踪法向；面板惯性法向与入射余弦；日距光压和总加速度。

SPICE 日地月状态与时间保持 real：

~~~text
d_S_DA = r_S_real(t) - r_sc_DA
~~~

所以仍保留对卫星状态的完整导数。

推荐 DA 布局：

| DA 变量 | 含义 |
|---:|---|
| 1–3 | 初始位置 |
| 4–6 | 初始速度 |
| 7 | 全局 SRP 比例误差 eta_SRP |
| 8–10 | 本体系三个姿态小角度 |
| 11 | 太阳翼转角偏差 |

布局不在几何模块硬编码；传播器按实际维数构造 srp_da_parameter_map_type。6 维保持纯状态传播，7 维保持现有尺度参数行为，额外变量只在映射有效时启用。质量、面积、压力和光学系数第一版为 real，相关误差优先吸收到全局尺度中。

受光判断和退化后备由 DA 常数项决定。测试点远离 mu=0 和退化边界；边界测试只要求有限、无 NaN 和分支正确。

DACE 生命周期沿用 pod_gravity_model_module 的显式临时池：进入时初始化，面板循环中复用，所有退出路径销毁，避免不可控的重载表达式临时句柄。

## 11. 与现有模块集成

- pod_force_model_module::compute_solar_radiation_pressure 按 srp_model 分派。炮弹球保持原公式；box-wing 分支增加速度输入并调用 real 接口。
- pod_da_force_model_module::da_compute_solar_radiation_pressure 同样分派，并把现有尺度参数设置转换成 DA 参数映射。新模块使用自己的 SRP 临时池。
- pod_uq_prop_da_module 集中构造映射；额外参数作为常参数随样本原样保留，兼容现有第 7 维策略。
- 第一版不修改 pod_spice.f90 的 CK 接口；未来 CK 适配器输出符合约定的 C_I_B 后调用外部姿态入口。

## 12. 错误处理与诊断

初始化和计算接口返回 status 与可选 message：0 成功，正值表示启用后备几何，负值表示不可计算。生产传播遇到负值明确终止，不返回零加速度或静默切换模型。

外部姿态矩阵校验：

~~~text
maxabs(transpose(C) C - I) <= 1e-10
abs(det(C)-1)              <= 1e-10
~~~

计算核心不写文件。调试模式可选返回姿态矩阵、逐面 mu、逐面加速度和太阳翼法向，由调用者决定输出。

## 13. 测试闭环

闭环顺序是“需求 → 独立单元测试 → real/DA 对照 → 主模型集成 → 全量回归”。所有测试由 FPM 注册并可单独运行。

### 13.1 配置闭环

test_spacecraft_geometry_config.f90：

- 默认 srp_model=cannonball。
- 完整 box-wing 配置正确解析。
- 非法质量、面积、轴向、光学系数和枚举校验失败。
- 炮弹球不要求 box-wing 参数；box-wing 缺参数不得回退。

### 13.2 面板物理闭环

test_srp_panel_force_real.f90 在固定 1 AU 和单位面积下验证完全吸收、完全镜面、漫反射手算值、掠入射、单面背光、双面激活和 km/s² 单位。标量极限容差：

~~~text
abs(x-x_ref) <= 1e-15 + 1e-12 abs(x_ref)
~~~

### 13.3 姿态闭环

test_srp_attitude_modes.f90 注入合成日地月状态以隔离 SPICE：

- 三种模式的指定本体系轴分别对准各自目标。
- 不会同时约束其他轴对准另外两个天体。
- transpose(C)C=I、det(C)=+1。
- 非默认指向轴正确。
- 近共线时启用后备轴且无 NaN。

轴对准、正交和行列式误差均不超过 1e-12。

### 13.4 太阳翼闭环

test_srp_solar_array_tracking.f90 验证固定模式、单轴法向单位化及其与铰链正交、朝阳符号、姿态变化、翼角偏差正负一阶效应，以及太阳方向平行铰链时的后备法向。

### 13.5 Box-wing 合力闭环

test_box_wing_srp_real.f90 验证合力等于逐面显式和、单面退化、整体旋转协变、非对称几何下三姿态结果可区分，以及太阳翼面积为零时只剩箱体。

### 13.6 real/DA 常数项闭环

test_box_wing_srp_da_consistency.f90 覆盖三种姿态与 fixed/single_axis：

~~~text
norm(a_DA_constant-a_real)
  <= 1e-18 + 1e-11 max(norm(a_real),1e-15) km/s²
~~~

### 13.7 DA 导数闭环

test_box_wing_srp_da_jacobian.f90 分别检查位置、速度、全局尺度、三个姿态小角度和太阳翼转角。real 用中心有限差分，DA 提取一阶系数；测试点远离分支边界。建议步长：位置 1e2–1e3 km，速度 1e-5 km/s，角度和尺度 1e-6。方向增量验收：

~~~text
norm(delta_a_DA-delta_a_FD)
  <= 1e-16 + 1e-5 max(norm(delta_a_DA),norm(delta_a_FD)) km/s²
~~~

理论零交叉导数不超过 1e-18 km/s²/归一化变量。

### 13.8 集成与回归闭环

- 默认炮弹球结果与修改前基准一致。
- 选择 box-wing 后 real 和 DA 主力模型确实调用新模块。
- use_srp=false 时两条路径 SRP 为零。
- 保留并通过现有 test_srp_param_da_propagation：6 维兼容、7 维参数保留并影响轨道。
- 新增 10/11 维短弧传播：额外参数原样保留，姿态或翼角变化影响轨道。
- 最终运行完整 fpm test，新旧测试全部通过。

### 13.9 需求—测试追踪矩阵

| 需求 | 直接测试 | 集成验证 |
|---|---|---|
| 六面箱体公式 | panel、box-wing real | real 主力模型 |
| 双面太阳翼 | panel、solar-array | real/DA 主力模型 |
| 日/地/月三模式 | attitude modes | 三模式短弧传播 |
| 不计算地影 | 固定光照算例 | SRP 全程非零基准 |
| real/DA 同公式 | DA consistency | 6/7/10/11 维传播 |
| 状态—姿态导数 | DA Jacobian | DA 短弧敏感性 |
| 尺度/姿态/翼角不确定性 | DA Jacobian | 7/10/11 维传播 |
| 炮弹球兼容 | config、baseline | 完整 fpm test |
| 外部姿态接口 | identity/known rotation | 暂不依赖 CK |

只有每项至少一个直接测试通过、所有集成测试通过且完整测试套件无回归时，功能才算完成。

## 14. 性能与可维护性

- 几何和光学参数只初始化一次，传播循环内不分配面板数组。
- real 与 DA 使用相同面板顺序和物理分解，减少公式漂移。
- SPICE 每步只查询必需目标；对日模式复用太阳状态。
- 全局只保存初始化后只读的默认几何；运行时状态和 DA 临时量局部化。
- 不引入网格、多次反射等超出本任务的框架。

## 15. 预计修改范围

新增：

- src/lib/forcemodel/pod_spacecraft_geometry_module.f90
- 第 13 节列出的单元与集成测试。

修改：

- src/lib/system/pod_config_module.f90
- config/config.txt
- src/lib/api/pod_engine_module.f90
- src/lib/forcemodel/pod_force_model_module.f90
- src/lib/forcemodel/pod_da_force_model_module.f90
- src/lib/uncertainty/propagation/pod_uq_prop_da_module.f90
- 必要时更新 fpm.toml 显式测试条目。

## 16. 完成标准

1. 三种解析姿态均可在 real 和 DA 主力模型运行。
2. 箱体与太阳翼公式、单位和符号通过解析极限测试。
3. DA 常数项与 real 一致，关键一阶导数与中心有限差分一致。
4. 现有 6/7 维 UQ 行为保持，新 10/11 维闭环通过。
5. 默认炮弹球无回归，完整 fpm test 通过。
6. 无 DACE 句柄泄漏、NaN 或静默模型切换。
7. 当前 Rosetta CK 未被误装载，外部姿态接口具有明确方向约定和校验。
