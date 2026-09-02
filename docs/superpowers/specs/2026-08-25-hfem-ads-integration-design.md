# HFEM ADS 误差传播集成设计

## 1. 目标

将现有 Fortran ADS 从“CRTBP 原型 + HFEM 测试内局部实现”整理为正式的 HFEM 不确定性传播方法，并接入统一 UQ 调度与 `run_HFEM_uprop` 应用入口。

完成后应支持：

- 使用完整 OPM 初始协方差生成相关轨道粒子；
- 以 `component` 或 `whitened` 坐标构造 ADS 根域；
- 将 SRP 相对尺度误差作为独立第七维 ADS 参数；
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

该模式与 C++ 参考实现的逐分量 `3-sigma` 标准化一致。在协方差为对角阵时，它等价于完整白化。

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

实现时考虑依托已有的src/math文件夹下的已有的 Cholesky 算法，避免重复实现。

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

批量求值接口负责编译 patch、多点定位、局部坐标映射、DA 求值和 `found` 状态返回。它不包含 HFEM、CRTBP、协方差或文件 I/O 知识。

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

### 5.4 OOP ADS 传播器

重写 `src/lib/uncertainty/propagation/pod_uq_prop_ads_module.f90`：

- 模块名改为 `pod_uq_ads_module`，与同级 MC/DA/UT 模块一致；
- 保留公开类型 `uq_ads_propagator`；
- 删除 CRTBP 依赖、`mu`、`set_ads_mu` 和硬编码粒子数；
- 增加 ADS 选项设置和 `last_stats`；
- 调用 HFEM ADS 服务层；
- 将最终粒子、均值和协方差完整保存在 `output_state`。

### 5.5 统一 UQ 调度

修改 `src/functions/orbitprop/pod_uq_propagation.f90`：

- 新增 `METHOD_ADS = 4`；
- 增加 `uq_ads_propagator` 分配分支；
- 接收可选 ADS 选项和统计输出；
- 支持内部样本或外部初始样本；
- 内部模式使用完整 OPM 协方差生成相关轨道粒子；
- `srp_sigma > 0` 时独立生成高斯 SRP 样本；
- 使用准确名义状态和输入协方差建域，不以有限样本估计量替代输入协方差。

## 6. 粒子输入

### 6.1 内部生成

内部模式由 `n_particles` 指定最终粒子数。轨道样本来自 OPM 完整协方差，SRP 样本来自独立高斯分布。对单位域外样本执行拒绝重采样，直到获得请求数量，并记录总抽样次数和拒绝数。

Gaussian 分布无界而 ADS 根域有限，因此输出文档和统计中必须明确粒子云是由 `domain_sigma` 截断的 Gaussian 分布。

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

扩展 `test/test_ads_core_dynamic_dimensions.f90`，检查六/七维批量求值、子 patch 定位、域外状态和 DA 资源释放。

### 10.3 HFEM 快速集成测试

新增 `test/test_hfem_ads_zero_duration_cloud.f90`，在零传播时间下验证两种坐标模式、SRP 第七列、粒子数、均值/协方差和域外计数。

### 10.4 HFEM 长传播测试

重构 `test/test_hfem_ads_propagation.f90`，删除内部重复 BFS 和 point evaluator，改用生产接口并继续验证 6D/7D、中心/内部/SRP 端点、SRP 响应、patch 覆盖和有限结果。

### 10.5 回归验证

至少运行：

```text
fpm build
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

本轮不做 CRTBP ADS 重写、2D 包络、manifold 长期序列化、全时序误差云保存，以及 SRP 与轨道初始误差的相关建模。

以下用户已有工作区改动不修改、不恢复：

- `rename_to_pod.sh`；
- `setup_env.sh`；
- 已删除的 `test/test_zero_init_cov_scripts.sh`。
