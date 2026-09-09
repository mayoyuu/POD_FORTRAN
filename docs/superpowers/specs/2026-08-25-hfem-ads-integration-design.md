# HFEM ADS 误差传播集成设计

## 1. 目标

将现有 Fortran ADS 从“CRTBP 原型 + HFEM 测试内局部实现”整理为正式的 HFEM 不确定性传播方法，并接入统一 UQ 调度与 `run_HFEM_uprop` 应用入口。

完成后应支持：

- 使用完整 OPM 初始协方差生成相关轨道粒子；
- 以 `component` 或 `whitened` 坐标构造 ADS 根域；
- 将 SRP 全局相对尺度误差作为独立第七维 ADS 参数；
- 从内部高斯采样或外部扰动文件获得初始误差云；
- 通过最终 ADS manifold 批量求值得到终端误差云；
- 保存粒子、最终统计矩和 ADS 分裂统计；
- 保持 CRTBP ADS 与 HFEM ADS 的动力学职责分离。

## 2. ADS 坐标约定

### 2.1 ADS 根域

ADS 的独立变量始终位于单位盒：

```text
xi(i) in [-1, 1]
```

`domain_sigma` 控制单位盒边界对应多少个标准差，默认值为 `3.0`。

### 2.2 Component 模式

默认模式为 `component`：

```text
B = domain_sigma * diag(sqrt(cov(i,i)))
x0 = mean + B * xi
```

物理偏差到 ADS 坐标的映射为：

```text
xi(i) = deviation(i) / B(i,i)
```

该模式与 C++ 参考实现的逐分量 `3-sigma` 标准化一致。在协方差为对角阵时，它等价于完整白化。复核 `ADS_cpp_Core_file/ads_prop_perturb.cpp` 可见，参考实现直接构造 `SV(i) + 3*sigma(i)*DA(i)`，并把物理扰动除以 `3*sigma(i)`；它没有读取完整协方差，也没有执行 Cholesky 白化。因此，“C++ 原实现已使用完整白化”不是当前代码事实。

选择它作为默认模式的原因是：

- 满足 ADS 对统一单位盒坐标的要求；
- 与当前 C++ 参考行为兼容；
- `split_count(1:6)` 可直接解释为 `x,y,z,vx,vy,vz` 的分裂次数；
- 完整协方差仍用于生成粒子，所以最终误差云保留相关性；
- 对零方差分量和半正定协方差更稳健。

Component 根域是协方差椭球的轴对齐包围盒，可能传播统计概率很低的盒角区域，因此 patch 数可能高于白化模式。

### 2.3 Whitened 模式

可选模式为 `whitened`。对六维轨道协方差作 Cholesky 分解：

```text
cov = L * transpose(L)
B = domain_sigma * L
x0 = mean + B * xi
```

物理偏差到 ADS 坐标的映射通过三角方程求解：

```text
B * xi = deviation
```

该模式使用完整协方差方向，通常产生体积更小、与概率分布更一致的根域。它要求六维轨道协方差正定；Cholesky 失败时程序明确报错并提示使用 `component`。

实现时调用 `src/lib/math/pod_basicmath_module.f90` 已公开的 LAPACK `dpotrf` 接口，在坐标模块中只封装尺寸、对称性、返回码和上三角清零，不重复实现 Cholesky 算法。

白化模式下 `xi(i)` 表示 `L` 第 `i` 列对应的联合方向，不再等同于单个物理状态分量。输出必须保存基矩阵 `B`，分裂标签使用 `xi1...xi6`，避免把白化方向分裂数误解释为物理分量分裂数。

### 2.4 分裂次数的科学含义

ADS 分裂方向由最终 DA 多项式最高阶项对超差输出分量的贡献决定。分裂数同时依赖：初始不确定度尺度、坐标参数化、传播时长、DA 阶数、截断误差阈值、最大分裂深度和动力学模型。

因此分裂数表示“在指定不确定度、参数化和误差容限下，各输入方向对非线性截断误差的贡献”，不是脱离实验条件的纯动力学不变量。

## 3. SRP 第七维

SRP 相对尺度误差记为 `eta_srp`，其一倍标准差由运行时参数 `srp_sigma` 输入：

```text
eta_srp = domain_sigma * srp_sigma * xi(7)
srp_multiplier = 1 + eta_srp
```

默认概率模型为 `eta_srp ~ N(0, srp_sigma^2)`。`srp_sigma = 0` 时禁用 SRP 不确定性并使用六维 ADS；`srp_sigma > 0` 时使用七维 ADS。

轨道初始误差与 SRP 参数默认相互独立。SRP 不是积分状态。每个 ADS patch 积分前必须利用 splitting history 重建该 patch 的 SRP 局部中心和跨度：

```text
eta_center = domain_sigma * srp_sigma * patch_center(7)
eta_span   = domain_sigma * srp_sigma * patch_width(7) / 2
srp_multiplier = 1 + eta_center + eta_span * xi_local(7)
```

第七维不参与六维轨道 Cholesky 白化，因此其分裂数在两种坐标模式下都可直接解释为 SRP 参数触发的分裂次数。

最新 box-wing SRP 代码已经约定 DA 参数编号：第 7 维为全局尺度，第 8--10 维为三轴姿态偏差，第 11 维为太阳翼角偏差。本轮 ADS 只实现用户已确认的“轨道 6 维 + 全局 SRP 尺度 1 维”，即最多初始化 7 个 DA 变量；不得因为配置文件中的 `srp_attitude_bias_span_arcsec` 或 `srp_array_angle_span_deg` 非零而隐式扩展到 10/11 维。姿态和太阳翼角不确定性留作后续独立扩展。

`srp_sigma` 是概率模型的一倍标准差，而 `config%srp_scale_da_span` 是现有普通 DA 传播器使用的 DA 系数跨度，二者语义不同。ADS 不从 `config%srp_scale_da_span` 推导概率标准差；它在每个 patch 积分前显式调用 `set_srp_scale_uncertainty(7, eta_center, eta_span)`。该接口同时作用于 cannonball 和 box-wing 模型中的全局尺度乘子。传播结束和所有错误退出路径都必须调用 `clear_srp_scale_uncertainty()`，避免模块级状态污染后续计算。

## 4. 目标架构

采用分层方案：

```text
pod_ads_split_module
    通用 Patch / Manifold / SplittingHistory 与批量求值
                |
pod_uq_ads_coordinates_module
    component / whitened 坐标映射
                |
pod_uq_hfem_ads_module
    HFEM ADS 建域、SRP 第七维、粒子求值与统计
                |
pod_uq_prop_ads_module
    uq_ads_propagator OOP 适配器
                |
pod_uq_propagation
    METHOD_ADS 统一调度
                |
run_HFEM_uprop
    CLI、输入选择和结果输出
```

CRTBP 路径继续由 `pod_uq_crtbp_ads_module` 独立负责。HFEM 模块不得依赖 CRTBP ADS 模块，CRTBP 模块也不得依赖 HFEM ADS 服务层。

## 5. 模块职责

### 5.1 通用 ADS core

修改 `src/lib/math/pod_ads_split_module.f90`，保留现有类型和分裂算法，新增 `mf_find_patch`、`mf_evaluate_point` 和 `mf_evaluate_points`。

批量求值接口负责多点定位、按 patch 分组、局部坐标映射、DA 求值和 `found` 状态返回。底层直接复用当前工作区已经新增的 `CompiledDA%eval_into` 和 `CompiledDA%eval_batch_into`；不在 ADS 任务中再次修改 C++ wrapper 或另建第二套 FFI。每个 patch 只编译一次，同一 patch 的输入整理为连续二维数组后批量求值。它不包含 HFEM、CRTBP、协方差或文件 I/O 知识。

### 5.2 ADS 坐标模块

新增 `src/lib/uncertainty/propagation/pod_uq_ads_coordinates_module.f90`，定义坐标模式常量、`ads_coordinate_map_type`、坐标映射构造、正向映射、逆向映射，以及协方差、零尺度和 Cholesky 错误检查。

### 5.3 HFEM ADS 服务层

新增 `src/lib/uncertainty/propagation/pod_uq_hfem_ads_module.f90`，定义：

- `hfem_ads_options_type`；
- `hfem_ads_stats_type`；
- `hfem_ads_build_domain`；
- `hfem_ads_propagate_samples`；
- `hfem_ads_propagate`。

该模块负责六维/七维初始 DA 映射、HFEM DA 积分、BFS 分裂、物理量纲截断误差、SRP patch-local 参数、最终 manifold 批量求值、统计计算和资源清理。

HFEM 服务层只选择并调用现有力模型，不复制 box-wing 几何、定姿或光压公式。`config%srp_model` 决定使用 cannonball 还是 box-wing；ADS 第七维对两者都只表示最终 SRP 加速度的全局相对尺度。

### 5.4 OOP ADS 传播器

重写 `src/lib/uncertainty/propagation/pod_uq_prop_ads_module.f90`：

- 模块名改为 `pod_uq_ads_module`，与同级 MC/DA/UT 模块一致；
- 保留公开类型 `uq_ads_propagator`；
- 删除 CRTBP 依赖、`mu`、`set_ads_mu` 和硬编码粒子数；
- 增加 ADS 选项设置和 `last_stats`；
- 调用 HFEM ADS 服务层；
- 将最终粒子、均值和协方差完整保存在 `output_state`；启用 SRP 时保留第七行 `eta_srp`，但动力学输出仍只有前六行轨道状态。

### 5.5 统一 UQ 调度

修改 `src/functions/orbitprop/pod_uq_propagation.f90`：

- 新增 `METHOD_ADS = 4`；
- 增加 `uq_ads_propagator` 分配分支；
- 接收可选 ADS 选项和统计输出；
- 支持内部样本或外部初始样本；
- 内部模式使用完整 OPM 协方差生成相关轨道粒子；
- `srp_sigma > 0` 时独立生成高斯 SRP 样本；
- 使用准确名义状态和输入协方差建域，不以有限样本估计量替代输入协方差。

现有 `run_particle_propagation` 和 `uq_da_propagator` 的 7D DA-MC 路径保持不变。ADS 不通过普通 DA 传播器间接实现，因为后者将 `config%srp_scale_da_span` 与样本偏差组合，其参数化不等同于 ADS 的 patch-local 第七维。

## 6. 粒子输入

### 6.1 内部生成

内部模式由 `n_particles` 指定最终粒子数。轨道样本来自 OPM 完整协方差，SRP 样本来自独立高斯分布。对单位域外样本执行拒绝重采样，直到获得请求数量，并记录总抽样次数和拒绝数。

Gaussian 分布无界而 ADS 根域有限，因此输出文档和统计中必须明确粒子云是由 `domain_sigma` 截断的 Gaussian 分布。

当前工作区新增的 `generate_antithetic_standard_normal_samples` 主要服务确定性诊断算例，要求偶数样本且使用总体矩匹配。本轮生产 ADS 默认仍沿用通用多元高斯采样和域外拒绝重采样，不把该诊断采样器设为隐式默认；否则会改变现有 MC/DA 的随机样本语义。

### 6.2 外部扰动文件

新增可选 `--perturb-file`：

- 六列表示六维物理轨道偏差；
- 七列的第七列表示物理 `eta_srp`；
- 六列文件配合 `srp_sigma > 0` 时，SRP 样本在内部独立生成；
- 域外样本不得写成零状态，而是排除出传播统计并计入 `outside_domain`；
- 若无有效样本，程序失败并返回清晰错误。

## 7. 文件 I/O

新增 `src/lib/data/pod_uq_sample_io_module.f90`，集中提供六/七列扰动 CSV 读取、动态维数粒子 CSV、moments JSON 和 ADS stats JSON 写出。

启用 SRP 时粒子 CSV 表头为 `x,y,z,vx,vy,vz,eta_srp`，否则为 `x,y,z,vx,vy,vz`。

ADS stats 至少保存：坐标模式、`domain_sigma`、`srp_sigma`、DA 参数、请求/输入/传播/域外/写出粒子数、patch 与 BFS 统计、分裂标签和次数、坐标基矩阵以及耗时。

最终无权重粒子的均值、协方差、偏度和峰度沿用 `uq_state_type%compute_moments` / `compute_higher_moments`，以保持 `run_HFEM_uprop` 现有结果口径，其中协方差分母为 `N-1`。新 `pod_uncertainty_diagnostics_module%compute_weighted_moments` 使用归一化权重的总体协方差并输出超额峰度，只用于加权诊断，不可直接替换生产 ADS moments 而不同时变更字段定义。`_ads_stats.json` 必须记录 `covariance_normalization = "sample_n_minus_1"`。

为避免影响近期 SRP 工作，本轮不重构 `app/run_srp_uq_propagation.f90`，也不改变 MC/DA/UT 既有文件名和字段。新动态维数 I/O 先服务 ADS，后续再单独统一其他入口。

## 8. CLI 集成

修改 `app/run_HFEM_uprop.f90`，正式支持 `-m MC|DA|UT|ADS`，并新增：

```text
--ads-coordinates component|whitened
--ads-domain-sigma <real>
--ads-max-depth <integer>
--ads-pos-tol <km>
--ads-vel-tol <km/s>
--srp-sigma <dimensionless 1-sigma>
--perturb-file <csv>
```

默认值为 component、3 sigma、最大深度 8、位置误差阈值 `1.0e-1 km`、速度误差阈值 `1.0e-6 km/s`、`srp_sigma=0`。

现有 `-da/--da-order` 同时表示 DA 和 ADS 多项式阶数，`-n/--n-particles` 控制内部粒子数或外部文件读取上限。

`--srp-sigma` 与现有 `run_srp_uq_propagation` 中同名参数保持“全局 SRP 相对尺度误差的 1-sigma”语义；它不等同于配置项 `srp_scale_da_span`。本轮 ADS 的 SRP 均值固定为零，名义 box-wing/cannonball 参数仍来自配置和现有力模型。

ADS 输出为：

```text
<prefix>_particles.csv
<prefix>_moments.json
<prefix>_ads_stats.json
```

## 9. 错误处理

- `domain_sigma` 必须大于零；
- `srp_sigma` 必须非负；
- whitened 模式要求六维协方差正定；
- component 模式允许零方差分量，但该分量的非零外部偏差视为域外；
- 外部文件每行必须包含六列或七列有限实数；
- 所有粒子均在域外时终止；
- 不允许以零向量代替无法定位 patch 的传播结果；
- 输出统计必须区分请求数、输入数、域内数和写出数。

## 10. 测试设计

### 10.1 坐标映射测试

新增 `test/test_ads_coordinate_mapping.f90`，检查 component/whitened 往返、对角等价性、非对角耦合、零方差和非正定错误路径。

### 10.2 ADS core 测试

扩展 `test/test_ads_core_dynamic_dimensions.f90`，检查六/七维批量求值、子 patch 定位、按 patch 分组、域外状态和 DA 资源释放。把当前工作区的 `test/test_compiled_da_eval_into.f90` 与 `test/test_compiled_da_batch_eval.f90` 作为底层接口前置回归，不在 ADS 测试中复制其覆盖内容。

### 10.3 HFEM 快速集成测试

新增 `test/test_hfem_ads_zero_duration_cloud.f90`，在零传播时间下验证两种坐标模式、SRP 第七列、粒子数、均值/协方差和域外计数。

### 10.4 HFEM 长传播测试

重构 `test/test_hfem_ads_propagation.f90`，删除内部重复 BFS 和 point evaluator，改用生产接口并继续验证 6D/7D、中心/内部/SRP 端点、SRP 响应、patch 覆盖和有限结果。至少各跑一个 cannonball 和 box-wing 的 7D 全局尺度算例，并断言活动 DA 变量数为 7、姿态/太阳翼角变量未被启用。

### 10.5 回归验证

至少运行：

```text
fpm build
fpm test test_compiled_da_eval_into
fpm test test_compiled_da_batch_eval
fpm test test_srp_da_parameter_map
fpm test test_ads_core_dynamic_dimensions
fpm test test_ads_coordinate_mapping
fpm test test_ads_domain_scale
fpm test test_hfem_ads_zero_duration_cloud
fpm test test_hfem_ads_propagation
fpm test test_uq_crtbp_comparison
```

最后使用小粒子数运行 `run_HFEM_uprop -m ADS` CLI 冒烟测试，再运行正式七天算例。

## 11. 文档

新增 `docs/hfem_ads_usage.md`，说明 ADS 根域与截断 Gaussian、两种坐标模式、分裂次数解释、SRP 一倍标准差参数、外部文件格式和输出字段。

## 12. 非目标和保护范围

本轮不做 CRTBP ADS 重写、2D 包络、manifold 长期序列化、全时序误差云保存、SRP 与轨道初始误差的相关建模，以及 box-wing 姿态偏差/太阳翼角偏差的 ADS 传播。

以下用户已有工作区改动不修改、不恢复：

- `external/dace_build/dace_wrapper.cpp`；
- `rename_to_pod.sh`；
- `setup_env.sh`；
- `src/lib/system/pod_dace_classes.f90`；
- `src/lib/statistics/pod_uncertainty_diagnostics_module.f90`；
- 当前未提交的 SRP、批量求值和 15 天诊断测试；
- 已删除的 `test/test_zero_init_cov_scripts.sh`。

## 13. 2026-09-09 最新代码冲突复核

复核范围包括设计提交 `2547c95` 之后的 `80af9d7`、`8570f7a`、`97d59e7`，以及当前未提交工作区。结论如下：

| 近期变化 | 与原设计的关系 | 修订后的处理 |
|---|---|---|
| box-wing SRP 与 7/8--10/11 维参数映射 | 有范围歧义，但无动力学接口冲突 | 本轮锁定 7D 全局尺度；复用现有力模型；8--11 维列为非目标 |
| `CompiledDA%eval_into` / `%eval_batch_into` | 与原计划的批量求值底层能力重叠 | ADS core 只做 manifold 定位和分组，直接调用现有批量接口，不修改 wrapper |
| `pod_uncertainty_diagnostics_module` | 与“统计计算”职责部分重叠，且协方差/峰度口径不同 | 生产粒子统计沿用 `uq_state_type`；诊断模块只用于加权分析 |
| 新的 15 天 SRP/初始状态研究测试 | 算例与 ADS 验证互补，不是生产入口 | 作为科学对照，不搬入 ADS 模块、不改其输出 |
| 测试程序从 `app/` 移到 `test/`、`fpm.toml` 清理 | 无接口冲突 | 新 ADS 测试继续放 `test/`；依赖 FPM 自动发现或按项目现行规则注册 |
| `run_HFEM_uprop` 与 `pod_uq_propagation` 尚未支持 ADS | 与原设计一致 | 仍按第 5、8 节补齐正式入口 |
| 旧 `pod_uq_prop_ads_module.f90` 仍依赖 CRTBP 且丢弃样本 | 与原设计一致，仍是必须修复的问题 | 重写为 HFEM OOP 适配器并保留最终粒子 |

因此当前设计不存在需要改换总体方向的冲突，但实现计划必须以上述修订边界为准，尤其不能覆盖尚未提交的 DACE 批量接口和 SRP 研究代码。
