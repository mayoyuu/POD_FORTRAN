# CRTBP + bicircular Sun geometry + SRP

后续 Q1–Q4 的绝对 β 炮弹球/平板 real、DA 接口见 [README_CRTBP_FLATPLATE.md](README_CRTBP_FLATPLATE.md)。本文件继续描述原相对 eta 接口，其 TU 默认保持原值；新研究接口默认使用初值 JSON 的 TU=375237.080249162 s。

本模块位于项目现有的 `src/lib/forcemodel/`（不是另建 `force_model/`），包含：

- `pod_bicircular_sun_module.f90`：确定性的双圆太阳几何与时间换算。
- `pod_crtbp_srp_module.f90`：real / DA 光压加速度、组合 CRTBP RHS、DA 参数初始化及临时池。

太阳仅提供方向和距离；没有太阳直接/间接引力、阴影、姿态或面元模型。
默认全日照。地月引力和旋转惯性项直接调用原有 `pod_crtbp_module`。
这一轮不修改 `app/`、积分器或 ADS 分裂算法。

## 坐标、时间和单位

地月质心旋转系中，地球在 `(-mu,0,0)`，月球在 `(1-mu,0,0)`。
状态为 `[x,y,z,dx/dtau,dy/dtau,dz/dtau]`，全部无量纲。
`mu` 沿用 `pod_crtbp_module` 的当前值，使用原有 `set_crtbp_mu(mu)` 设置；
它是项目现有的全局量，不能在并发调用之间设置成不同值。

默认尺度严格沿用 `app/run_CRTBP_uprop.f90`：

| 字段 | 默认值 | 单位 |
|---|---:|---|
| `cfg%sun%du_m` | 384400000 | m |
| `cfg%sun%tu_s` | 375190.26189464843 | s |
| `cfg%sun%au_m` | 149597870700 | m |
| `cfg%sun%radius_m` | 149597870700 | m |
| `cfg%sun%period_s` | 31558149.8 | s |
| `cfg%sun%alpha0_rad` | 0 | rad |
| `cfg%pressure_1au` | 4.56e-6 | N/m² |
| `cfg%cr0` | 1.2 | 无量纲 |
| `cfg%area_mass` | 0.02 | m²/kg |

任务书示例 GM 与现有 CRTBP 常量略有差别，因此这里不从示例 GM 重新推导
mu/TU。DU、TU 为本接口的权威尺度；如果切换系统，调用者须同时提供一致的
mu、DU、TU（隐含 `GM_total=DU³/TU²`）。SRP 模块不改写全局 mu。

`tau=(t_seconds-t0_seconds)/TU`；可用 `sun_elapsed_tau(t_seconds,t0_seconds,cfg%sun)`
转换。`alpha0_rad` 是 tau=0 时刻的太阳相位。积分中传入连续、相对同一个 t0 的
tau，续积分或 ADS 分裂时不要把 tau 重置为零。

设 `n=1/TU`，`n_s=2*pi/period_s`，旋转系太阳角速度为 `omega=n_s*TU-1`。

```text
alpha(tau) = alpha0_rad + omega*tau
r_sun      = (radius_m/du_m) [cos(alpha), sin(alpha), 0]
```

太阳在默认地月会合系中反向转动；周期由 `2*pi/abs(omega)` 得出。
这等价于对惯性太阳方向施加 `Rz(-tau)`，再将惯性方向角 `alpha0+n_s*TU*tau`
代入。力矢量转到旋转系只需旋转分量；原 CRTBP RHS 已有离心和科里奥利项。

## 力模型

`eta=[eta_r,eta_t,eta_n]` 为弧段内固定的无量纲模型参数，名义值为零。

```text
d      = r - r_sun                     # 太阳指向航天器
er     = d / |d|
q      = (AU/DU)^2 / |d|^2
gamma0 = P0 * Cr0 * (A/m)0 * TU^2/DU
```

| `cfg%model` | 加速度 |
|---|---|
| `SRP_OFF` | 0 |
| `SRP_CANNONBALL_FULL`（默认） | `gamma0*q*(1+eta_r)*er` |
| `SRP_CANNONBALL_FARFIELD` | 同上，但令 `d=-r_sun`，忽略航天器视差 |
| `SRP_EMPIRICAL_THREE_COMPONENT` | `gamma0*q*((1+eta_r)*er+eta_t*et+eta_n*en)` |

M2 中 `et=normalize(k cross er)`，`k=(0,0,1)`，`en=er cross et`。
当 `|k cross er|²<1e-20` 时改用 x 轴作参考；`used_fallback` 返回是否触发。
M2 不重新归一化合成向量，横向/法向参数表示经验加速度分量，不是姿态角。
M1 忽略 `eta_t,eta_n`，只应激活 `eta_r`。

`enabled=.false.`、`SRP_OFF` 或任一名义光压因子为零时，光压严格为零。
这条快速路径不计算太阳几何。全零 eta 在启用模型时仍有名义 SRP。
SI 加速度为 `a_nd*DU/TU²`；默认 1 AU 远场幅值是 `1.0944e-7 m/s²`。
太阳位于 +x 时，质心附近航天器的光压指向 -x。
完整模型允许 z 不为零；远场默认几何的直接 z 加速度为零。

## 调用接口

```fortran
use pod_global, only: DP
use pod_crtbp_srp_module

type(CrtbpSrpConfig) :: cfg
real(DP) :: x(6), rhs(6), eta(3), a(3), a_si(3), tau
! 先设置 x 和 tau；x 是无量纲会合系状态。
eta = [0.05_DP, 0.0_DP, 0.0_DP]
call srp_acceleration_real(tau, x(1:3), eta, cfg, a, a_si)
call crtbp_srp_derivatives_real(x, rhs, tau, cfg, eta)
```

也可调用泛型 `srp_acceleration` 和 `crtbp_srp_derivatives`，由 real/DA 类型分派。
RHS 输出为 `dx/dtau`，仅在第 4–6 项加光压。
加速度接口还接受可选的 `prescribed_sun_nd(3)`，用于替换几何 provider；其语义
必须是同一时刻、同一会合系、相同 DU 下的“质心→太阳”位置。RHS 当前默认使用
双圆几何；更换 provider 时由调用者组合原 CRTBP RHS 和加速度接口。

DA 例子（所有声明应放在程序声明区）：

```fortran
use pod_dace_classes
type(DA) :: x_da(6), rhs_da(6), eta_da(3)
type(CrtbpSrpForceTempPool) :: pool
integer :: i

! 仅在创建本轮任何 DA 对象之前初始化；已有句柄必须先 destroy。
call dace_initialize(4, 1)
call pool%init()
do i=1,6
    x_da(i) = x(i)                      ! 六维初值均为常数
end do
call srp_init_da_parameters([.true.,.false.,.false.], &
    [0.30_DP,0.0_DP,0.0_DP], eta_da)     ! eta_r=0.30*xi_1
call crtbp_srp_derivatives_da(x_da, rhs_da, tau, cfg, eta_da, pool)

! 积分时始终复用 eta_da 和 pool；状态因参数影响逐步成为 DA 多项式。
! 结束后显式释放数组元素，符合项目 DACE 句柄生命周期约定。
do i=1,6
    call x_da(i)%destroy()
    call rhs_da(i)%destroy()
end do
do i=1,3
    call eta_da(i)%destroy()
end do
call pool%destroy()
```

DA 加速度签名为 `srp_acceleration_da(tau,r,eta,cfg,a,pool,a_si,used_fallback,prescribed_sun_nd)`。
`a_si` 也为 DA 数组，可选。输出数组与输入数组及 pool 必须是不同对象。
池应 init 一次、重复调用、最后 destroy；不得将 pool 的成员再作为输出实参传入。

三参数模型用 `dace_initialize(order,3)`、`active=[.true.,.true.,.true.]`。
辅助函数按 r、t、n 顺序紧凑分配编号，检查 DACE 变量数恰等于激活参数数。
例如只激活 n 时，`eta_n=s_n*xi_1`，其余为零。没有激活参数时直接使用 real。
辅助函数用于初始参数域；也允许调用者直接传入带偏置或重映射后的 eta DA 多项式。
不得在每次 RHS 调用中重新调用初始化或重置 eta。

若末端 `state_da(i)` 对 `xi_1` 的一次系数为 c，则 `dX_i/deta_r=c/scale_r`。
可用 `state_da(i)%get_coeff([1])/scale_r` 提取；多参数使用对应的单位指数向量。
尺度表示区间半宽，不自动代表概率分布或 1σ。

**后续 ADS 接口约定：** 分裂参数域时，对状态多项式和 eta 多项式应用同一个
变量仿射替换（例如 `xi -> center + halfwidth*xi_local`）。仅分裂六维初始状态
会丢失参数域，因为初值全部为常数。现有 CRTBP ADS/积分器硬连接纯 CRTBP RHS，
本次新接口尚未自动接入它们；下一阶段需要接入新 RHS、携带参数域并记录分裂时间。

## 验证

GNU Fortran 必须 >=13.3.0，遵循主 README 的 DA 临时对象析构要求。
不需要 SPICE 星历加载；构建仍沿用项目现有链接库。

一键验证：`bash test/run_crtbp_srp_checks.sh`。日志和可选时间序列写入
`build/crtbp_srp_validation/`；包含正常用例和四个预期失败的输入检查。
环境变量 `FPM_FC` 可指定符合版本要求的编译器。

```bash
cd POD_Fortran
source setup_env.sh
fpm test --compiler gfortran-13 --profile debug --flag '-ffree-line-length-none' --target test_crtbp_srp
fpm test --compiler gfortran-13 --profile debug --flag '-ffree-line-length-none' --target test_crtbp_srp_propagation
# 可选：输出时间序列，路径传给测试程序
fpm test --compiler gfortran-13 --profile debug --flag '-ffree-line-length-none' --target test_crtbp_srp_propagation -- /tmp/crtbp_srp_15day.csv
# 以下两条应非零退出，并明确报告单位或日距错误：
fpm test --compiler gfortran-13 --profile debug --flag '-ffree-line-length-none' --target test_crtbp_srp -- invalid_units
fpm test --compiler gfortran-13 --profile debug --flag '-ffree-line-length-none' --target test_crtbp_srp -- sun_collision
```

力测试覆盖符号、周期、SI/ND 换算、逆平方律、三维视差、M1/M2 退化、
开关退化、DA 零阶/一阶、传播状态的 DA 几何链式依赖以及句柄计数。
传播测试使用现有 CRTBP 测试中的初始状态，执行 15 天、2048 步的测试专用 RK4：
对照 SRP off、名义 SRP、四阶单参数 DA、eta=±0.05 的 real 传播；另用 ±0.001
核对敏感性。实数名义轨道以步长减半核对积分误差，SRP off 检查 Jacobi 漂移。
这不是生产 ADS 积分器，也不声称该初值已在本次工作中校正成周期轨道。

CSV 给出时间、太阳角/位置、名义 SRP、off/nominal/±0.05 轨道；
参数效应应相对名义 SRP 轨道计算，名义光压效应则相对 off 轨道计算。
测试配置在源码中固定且打印 DA 阶数/变量数/尺度；批处理元数据输出留给 app 阶段。

2026-11-01 在 GNU Fortran 13.3.0 debug 构建中运行上述验证脚本全部通过。
该固定测试用例的结果如下（混合位置/速度分量的最大值均在无量纲状态中计算）：

| 检查 | 结果 |
|---|---:|
| 15 天 DA 与同一步长 real 的末端最大差，eta=±0.05 | 4.2571e-14 |
| real 名义传播步长减半的末端最大差 | 1.1039e-10 |
| dX/deta 与有限差分的最大差 | 2.8562e-11 |
| SRP off 的 Jacobi 常数漂移 | 4.8272e-13 |

第一项检验 DA 参数映射的一致性，不代表轨道的绝对积分误差；积分误差另由
第二项检查。力模型用例、单参数/三参数 DA 用例、句柄释放检查，以及四个
预期失败输入用例均通过；独立只读代码审查未发现可操作的问题。

## 适用范围与错误

参数/单位无效和距离小于等于 `1e-12 DU` 时通过明确的 `error stop` 终止。
DA 对日距的检查在展开中心完成；它不能证明整个参数域不穿越奇异点或 M2
备用轴切换边界。DA 域必须留在同一光滑分支内，这一点需要后续 ADS 调用者保证。
该极端几何远离默认地月轨道范围。配置和太阳时间保持 real；距离、方向、q
在 DA 路径上完整保留多项式，不取常数项代替。

模型依据：任务书 `CRTBP_SRP_Force_Model_Codex_Spec.md` 的数学定义。
炮弹球光压方向与幅值亦可参见
[CannonballSRP 模型文档](https://ccar-orcca.github.io/scarabaeus-docs/api_ref/dynamics/cannonballsrp.html)。
