# CRTBP 平板 SRP：Q1–Q4 的力模型与参数接口

本次完成 `src/lib/forcemodel/` 中的 real / DA 力模型和参数域，供后续单工况传播、ADS 和批量实验调用。实现依据三层实验设计及用户最后修订：取消 HZ；关注均值、协方差和 ADS 分裂方向；面积 180 m² 为笔误，参考面积采用 80 / 160 m²。Q1–Q4 是研究问题，三层参数维度依次为 1 / 3 / 4。

当前未接入 JSON 批处理、逐小时统计、首次/逐次分裂日志或最终分布保存。以下接口不把参数区间自动解释成某种概率分布，也不设置 ADS 阶数、阈值或分裂策略。不存在等待另一份 ADS 文档的前置条件；这些属于下一步基于现有 ADS 的实现与收敛验证。

## 文件与接口

| 文件 | 提供内容 |
|---|---|
| `src/lib/forcemodel/pod_srp_study_parameters_module.f90` | 物理参数索引、光学族、参考面积质量转换、区间/中点映射、DA 活跃变量对应关系 |
| `src/lib/forcemodel/pod_crtbp_srp_study_module.f90` | 绝对 β 炮弹球/平板加速度、CRTBP 组合 RHS、real / DA 版本和 DA 临时池 |
| `src/lib/forcemodel/pod_bicircular_sun_module.f90` | 复用确定性双圆太阳几何；没有太阳引力 |
| `test/test_crtbp_flatplate.f90` | 力、光学系数、参数域、导数和资源释放测试 |
| `test/test_crtbp_flatplate_propagation.f90` | 整段轨迹等效、坐标变换和参数 DA 传播回归 |

原 `pod_crtbp_srp_module` 的相对 `eta` 接口保留，见 `README_CRTBP_SRP.md`。新接口用绝对 β，不能将两个接口的参数数组互换。

## 坐标、单位与太阳

状态为无量纲地月质心会合系 `[x,y,z,vx,vy,vz]`；时间 `tau` 为相对同一物理历元的无量纲时间。组合 RHS 复用原 CRTBP 引力、科里奥利和离心项，仅向第 4–6 项加入 SRP。光压力已在会合系表达，无须再加旋转惯性项。

新 `CrtbpSrpStudyConfig` 默认使用初值 JSON 的尺度：

| 量 | 数值 |
|---|---:|
| `cfg%sun%du_m` | 384400000 m |
| `cfg%sun%tu_s` | 375237.080249162 s |
| `cfg%sun%au_m`、`radius_m` | 149597870700 m |
| `cfg%sun%period_s` | 31558149.8 s |
| `cfg%sun%alpha0_rad` | 0 rad |
| `cfg%pressure_1au` | 4.56e-6 N/m² |

在传播入口调用 `set_crtbp_mu(0.012153614091891635_DP)`；现有 CRTBP 默认约化值 `0.012153614091892` 也符合用户要求。`mu` 是既有模块的全局量，不由光压模块修改。旧 `BicircularSunConfig()` / `CrtbpSrpConfig()` 的 TU 默认仍为 375190.26189464843 s；新研究配置有意覆盖为目录 TU，跨接口对照时须统一太阳配置。

```text
alpha(tau) = alpha0 + (2*pi*TU/period_s - 1)*tau
r_sun      = (radius_m/DU) * [cos(alpha), sin(alpha), 0]
d          = r_sun - r
s          = d/|d|                       # 航天器指向太阳
q          = (AU/DU)^2 / |d|^2
P_nd       = pressure_1au * TU^2 / DU
```

太阳只提供方向、距离和随时间变化的相位，全日照，无太阳引力、阴影、力矩或热再辐射。炮弹球和平板统一采用完整日距因子 `q`；这是沿用此前完整日距实现的选择，与附件“约 1 AU 常压”建议有区别。`alpha0_rad` 是太阳初相，与 JSON 的轨道初相分开设置。ADS 分裂或续传播时保持连续的 `tau`。

## 平板参数和物理公式

物理数组始终是 `p(4)=[beta,chi,u,v]`，可以用 `SRP_BETA/SRP_CHI/SRP_U/SRP_V` 索引。

- `beta`：正入射等效 `Cr*A/m`，单位 m²/kg；取 0 代表关闭 SRP 的物理参考。
- `chi`：镜面/漫反射组成，范围 `[0,1]`；`rho_s=0.5*chi`，`rho_d=0.75*(1-chi)`，吸收份额 `0.25+0.25*chi`。始终满足 `C_eq=1+rho_s+2*rho_d/3=1.5`。
- `u,v`：无量纲太阳切向斜率。参数沿弧段不变，太阳切向基随位置和时间变化，因此这不是固定惯性法向的姿态模型。

默认 `chi=0.6`，得到 `rho_s=rho_d=0.3`。`srp_face_on_beta(area_m2,mass_kg)` 返回 `1.5*area_m2/mass_kg`：质量 1200 kg、面积 80 / 160 m² 对应 β=0.1 / 0.2。区间内变化的 β 是等效强度参数，不能解释成同一块固定面积平板仅靠改变材料就能从零光压变化到该上限。

切向基由会合系固定参考轴 `k=cfg%tangent_reference_axis` 定义，默认 z 轴：

```text
b1 = normalize(k cross s)
b2 = s cross b1
e1 = cos(delta)*b1 + sin(delta)*b2
e2 = -sin(delta)*b1 + cos(delta)*b2
c  = 1/sqrt(1+u*u+v*v)
n  = c*(s + u*e1 + v*e2)

a_CB = -P_nd*q*beta*s
a_FP = -P_nd*q*(beta/1.5)*c * [(1-rho_s)*s + 2*(rho_s*c+rho_d/3)*n]
```

`delta=cfg%tangent_rotation_rad` 默认为 0，供固定切向基旋转对照使用。改变 delta 时，同一物理法向要求同时反向旋转 `(u,v)`。不自动切换参考轴；当固定轴与太阳方向平行时明确报错，避免在 DA 域内引入不连续的坐标切换。应选择在整个研究弧段/参数域中均不接近极点的固定轴。

有限 `u,v` 保证 `c=n·s>0`，基线不涉及正背面切换。参数工厂只接受最大倾角 `[0,90)` 度。这里没有额外的 `Cr` 乘子。`u=v=0` 时，对所有允许的 χ，平板严格退化为炮弹球。DA 中位置、日距和切向基全部保留多项式依赖；不会因 β 的常数项为零而丢掉其非零导数。

## Q1–Q4 参数域

`SrpStudyDomain` 保存四个参数的 `lower/upper`，相等的上下界表示固定参数。`center()` 与 `halfwidth()` 返回中点和半宽。`srp_study_init_da` 按物理顺序压缩活跃自变量，返回 `index_map(physical_index)=DA_index`，固定参数对应 0。

| 研究 | 模型/参数 | 示例域 | `index_map` |
|---|---|---|---|
| Q1 强度 | Cannonball，β | `[0,0.1]` 或 `[0,0.2]`；χ=.6，u=v=0 | `[1,0,0,0]` |
| Q2 强度/姿态 | Flat plate，β,u,v | 局部 β `[.08,.12]` 或主/扩展域；χ=.6 | `[1,0,2,3]` |
| Q3 分裂方向 | 在 Q2/Q4 参数域上记录分裂 | 需同时保存参数轴定义、局部域和尺度 | 同相应参数域 |
| Q4 反射/姿态耦合 | Flat plate，β,χ,u,v | χ `[0,1]`；窄域对照 `[.4,.8]` | `[1,2,3,4]` |

例如 `tilt_max_deg=10` 产生 `u,v∈[-q_u,q_u]`，其中 `q_u=tan(10°)/sqrt(2)`。矩形角点倾角恰为 10°，并非整个均匀圆锥。5°/20°/30° 等域通过同一接口显式选择。也可直接设置域的上下界，构造固定 β 的纯姿态对照。

`chi∈[0,1]` 的 DA 中心是 0.5，物理基准 χ=0.6 需在多项式对应位置求值；不能混淆两者。主域 β `[0,.1]` 的中心是 .05；扩展域 `[0,.2]` 的中心是 .1；两者的物理参考仍是 β=0。

## 调用示例

real 组合 RHS：

```fortran
use pod_global, only: DP
use pod_crtbp_module, only: set_crtbp_mu
use pod_crtbp_srp_study_module
type(CrtbpSrpStudyConfig) :: cfg
real(DP) :: x(6), rhs(6), p(4), tau
! 由调用者填写 x、tau；默认已经选择 STUDY_FLAT_PLATE。
call set_crtbp_mu(0.012153614091891635_DP)
p = [0.1_DP, 0.6_DP, 0.02_DP, -0.03_DP]
call crtbp_srp_study_derivatives_real(x, rhs, tau, cfg, p)
! 炮弹球：cfg%model=STUDY_CANNONBALL；关闭光压：cfg%enabled=.false.
```

DA 初始化和 RHS（六维初始状态全部为常数）：

```fortran
use pod_dace_classes
use pod_srp_study_parameters_module
! 其余 use、cfg、x、tau 同上；声明放在可执行语句前。
type(SrpStudyDomain) :: domain
type(SrpStudyTempPool) :: pool
type(DA) :: xd(6), fd(6), pd(4)
integer :: i, index_map(4)
domain = srp_study_domain(0.0_DP, 0.1_DP, tilt_max_deg=10.0_DP)
! Q4 时再加 chi_lower=0.0_DP, chi_upper=1.0_DP。
! 6 仅是示例阶数，后续 ADS 实验仍需阶数/阈值收敛验证。
call dace_initialize(6, domain%n_variables())
call pool%init()
call srp_study_init_da(domain, pd, index_map)
do i=1,6
    xd(i)=x(i)
end do
call crtbp_srp_study_derivatives_da(xd, fd, tau, cfg, pd, pool)
! 传播结束后显式释放数组内的 DACE handles。
do i=1,6
    call xd(i)%destroy()
    call fd(i)%destroy()
end do
do i=1,4
    call pd(i)%destroy()
end do
call pool%destroy()
```

独立加速度接口为 `srp_study_acceleration_real(tau,r,p,cfg,a,a_si)` 和 `srp_study_acceleration_da(tau,r,p,cfg,a,pool,a_si)`，`a_si` 可省略，SI 单位为 m/s²，方向仍是会合系分量。泛型名 `srp_study_acceleration` 和 `crtbp_srp_study_derivatives` 按参数类型分派。调用前初始化池，输入、输出、池不可共用同一对象；DACE 重设变量数前释放旧 handles。

## 下一阶段的衔接约定

单工况入口接收轨道初值、传播时间、太阳初相、模型和参数域；JSON 的枚举与批处理留给外层调度。终点保持用户规则：先取 2T，仅当 `2*T*TU < 15*86400` 时改取 3T。不是补足任意整数周期，也不是固定 15 天。

ADS 子域必须对**状态多项式和参数多项式同时**作相同的仿射替换，时间连续。不能每次分裂后重新调用根域的初始化函数。分裂方向应记录活跃 DA 索引和对应物理名；三维中的第 2 轴是 u，四维中的第 2 轴是 χ。

逐小时的 `∂x/∂beta | beta=0` 要在包含 β=0 的局部映射中，对多项式求导再在相应坐标求值，并除以该局部 β 半宽。根域主/扩展实验的 β=0 都对应 `xi_beta=-1`；直接读一次项得到的是中心导数，不是 β=0 的导数。比较平板基准时明确其它参数，例如 χ=.6、u=v=0。

均值、方差/协方差必须基于明确的参数概率分布和分裂后权重计算；保存参数区间本身不等于定义先验。HZ 已取消。分裂次数最多的参数轴也不能直接等同于最重要的物理因素，仍需保留域尺度、判据与敏感度作对照。

## 验证

在 WSL/Linux 项目目录运行：

```bash
bash test/run_crtbp_flatplate_checks.sh
bash test/run_crtbp_srp_checks.sh
```

脚本默认 `gfortran-13`，可用 `FPM_FC` 覆盖；依赖环境由项目 `setup_env.sh` 加载。新测试日志在 `build/crtbp_flatplate_validation/`，旧模型回归日志在 `build/crtbp_srp_validation/`。

新增回归采用测试专用 RK4、4096 步、约 19.99347 天；包含整个弧段的 β=0、对日等效/χ 不变性和切向基旋转不变性。三/四维 DA 使用六阶、β `[.08,.12]`、最大倾角 3°，四维 χ `[.4,.8]`，检查中心、内部点与三个角点，并以半步长实数积分复核时间离散误差。这只验证当前力模型及参数传播接口，不能替代全域 10°/30°、全轨道族和 ADS 分裂后的收敛实验。

2026-11-02，GNU Fortran 13.3.0、debug 构建的实测结果：

| 检查 | 结果（状态差为无量纲） |
|---|---:|
| β=0 与关闭 SRP 的全弧段差 | 0 |
| 对日平板/炮弹球，χ=0/.6/1 的全弧段最大差 | 0 |
| 切向基旋转后的全弧段最大差 | 8.47e-22 |
| 3D DA 末端最大绝对误差 | 5.54e-12 |
| 4D DA 末端最大绝对误差 | 2.70e-11 |
| 3D / 4D 实数角点半步长最大差 | 2.76e-11 / 2.78e-11 |
| 力模型、光学极限、混合导数、DA 资源释放 | 通过 |
| 非法参数域/倾角、切向极点、DA 维度不匹配 | 均按预期拒绝 |
| 原 eta 力模型和 15 天 DA 传播回归 | 通过 |
