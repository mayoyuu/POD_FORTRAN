# HFEM ADS Integration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 将 Fortran ADS 正式接入 HFEM 与 `run_HFEM_uprop`，支持可选第七维全局 SRP 尺度误差，并保存最终误差云、统计矩、粒子数和分裂统计。

**Architecture:** 通用 ADS core 只管理 patch/manifold 和批量求值；独立坐标模块完成 component/whitened 映射；HFEM 服务层负责建域、分裂、SRP patch-local 参数和粒子传播；OOP 适配器接入统一 UQ；应用层只解析参数和写文件。默认是 6D component，只有 `srp_sigma > 0` 时才创建 7D ADS。

**Tech Stack:** Fortran 2008、FPM、DACE/CompiledDA、C++ wrapper、OpenBLAS/LAPACK `dpotrf`、SPICE、Git。

## Global Constraints

- 在分支 `feature/hfem-ads-integration` 上工作，不修改或提交用户现有未提交文件。
- ADS 默认坐标为 `component`，`whitened` 仅为显式可选模式。
- ADS 根坐标始终是 `[-1,1]`；`domain_sigma` 默认 `3.0`。
- 第七维是可选项：`srp_sigma == 0` 使用 6D，`srp_sigma > 0` 使用 7D。
- 第七维只表示全局 SRP 相对尺度误差；box-wing 的姿态三轴和太阳翼角第 8--11 维不在本轮范围。
- `srp_sigma` 是概率 1-sigma，不读取 `config%srp_scale_da_span` 代替。
- 最终误差云必须保存在 `output_state%samples` 和 `<prefix>_particles.csv`。
- 最终协方差沿用 `uq_state_type%compute_moments` 的 `N-1` 口径。
- 复用工作区已有 `CompiledDA%eval_into` / `eval_batch_into`，不修改其 wrapper。
- 每个生产行为必须先有失败测试，再写最小实现。

---

### Task 1: ADS 坐标映射与可选第七维

**Files:**
- Create: `src/lib/uncertainty/propagation/pod_uq_ads_coordinates_module.f90`
- Create: `test/test_ads_coordinate_mapping.f90`

**Interfaces:**
- Produces: `ADS_COORD_COMPONENT`、`ADS_COORD_WHITENED`、`ads_coordinate_map_type`
- Produces: `ads_build_coordinate_map(mean6,cov6,mode,domain_sigma,srp_sigma,map,status,message)`
- Produces: `ads_physical_to_unit(map,deviation,eta_srp,unit_point,status)`
- Produces: `ads_unit_to_physical(map,unit_point,deviation,eta_srp,status)`

- [ ] **Step 1: 写坐标与维数失败测试**

测试必须断言：默认/显式 component 的 `B=domain_sigma*diag(sqrt(cov_ii))`；whitened 的 `B*B^T=domain_sigma^2*cov`；`srp_sigma=0` 时 `map%n_variables=6`，正值时为 7；第七维往返满足 `eta=domain_sigma*srp_sigma*xi7`；非正定 whitened、负 sigma、零 `domain_sigma` 返回负状态。

- [ ] **Step 2: 运行 RED**

Run: `fpm test test_ads_coordinate_mapping`

Expected: FAIL，原因是 `pod_uq_ads_coordinates_module` 尚不存在。

- [ ] **Step 3: 实现最小坐标模块**

核心类型和接口固定为：

```fortran
integer,parameter,public :: ADS_COORD_COMPONENT=1,ADS_COORD_WHITENED=2
type,public :: ads_coordinate_map_type
    integer :: mode=ADS_COORD_COMPONENT
    integer :: n_variables=6
    real(DP) :: domain_sigma=3.0_DP
    real(DP) :: srp_sigma=0.0_DP
    real(DP) :: mean6(6)=0.0_DP
    real(DP) :: basis6(6,6)=0.0_DP
end type
```

whitened 分解调用 `pod_basicmath_module::dpotrf('L',6,...)`，成功后显式清零上三角；逆映射用前代求解 `basis6*xi=deviation`。component 零方差分量只接受零偏差。

- [ ] **Step 4: 运行 GREEN**

Run: `fpm test test_ads_coordinate_mapping`

Expected: PASS。

- [ ] **Step 5: 提交**

```bash
git add src/lib/uncertainty/propagation/pod_uq_ads_coordinates_module.f90 test/test_ads_coordinate_mapping.f90
git commit -m "feat: add ADS coordinate mappings"
```

### Task 2: Manifold 定位与批量求值

**Files:**
- Modify: `src/lib/math/pod_ads_split_module.f90`
- Modify: `test/test_ads_core_dynamic_dimensions.f90`

**Interfaces:**
- Consumes: `CompiledDA%eval_into`、`CompiledDA%eval_batch_into`
- Produces: `mf_find_patch(manifold,point,patch_index)`
- Produces: `mf_evaluate_point(manifold,point,value,found,status)`
- Produces: `mf_evaluate_points(manifold,points,values,found,status)`

- [ ] **Step 1: 扩展失败测试**

构造两个 7D 子 patch，批量输入含左域、右域、边界和域外点；断言批量结果等于逐点 DA 求值、`found` 正确、域外输出不伪造为零传播状态、调用前后 `active_da_count` 不变。

- [ ] **Step 2: 运行 RED**

Run: `fpm test test_ads_core_dynamic_dimensions`

Expected: FAIL，原因是新 manifold 接口不存在。

- [ ] **Step 3: 实现定位和分组批量求值**

`mf_find_patch` 按现有 patch 顺序返回第一个包含点的索引；`mf_evaluate_points` 先定位所有列，再为每个 patch 分配连续局部输入矩阵，只编译一次并调用 `eval_batch_into`，最后销毁 compiled handle。未找到点的 `found=.false.`，其输出设 IEEE quiet NaN，调用方不得纳入统计。

- [ ] **Step 4: 运行 GREEN 与底层回归**

Run:

```bash
fpm test test_compiled_da_eval_into
fpm test test_compiled_da_batch_eval
fpm test test_ads_core_dynamic_dimensions
```

Expected: 全部 PASS。

- [ ] **Step 5: 提交**

```bash
git add src/lib/math/pod_ads_split_module.f90 test/test_ads_core_dynamic_dimensions.f90
git commit -m "feat: evaluate ADS manifolds in batches"
```

### Task 3: HFEM ADS 服务层

**Files:**
- Create: `src/lib/uncertainty/propagation/pod_uq_hfem_ads_module.f90`
- Create: `test/test_hfem_ads_zero_duration_cloud.f90`
- Modify: `test/test_hfem_ads_propagation.f90`

**Interfaces:**
- Consumes: `ads_coordinate_map_type`、`manifold_type`、HFEM DA 积分器、`set_srp_scale_uncertainty`
- Produces: `hfem_ads_options_type`
- Produces: `hfem_ads_stats_type`
- Produces: `hfem_ads_propagate(nominal_state,covariance,epoch0,t_start,t_end,input_samples,options,output_samples,stats,status,message)`

- [ ] **Step 1: 写零时长失败测试**

分别测试 6D component、6D whitened 和 7D component。输入粒子使用已知相关协方差内点；`t_end=t_start` 时前六维输出应等于初始轨道粒子，第七维启用时原样保留；输出列数等于域内粒子数；`stats%requested/input/inside/propagated/written` 一致。

- [ ] **Step 2: 运行 RED**

Run: `fpm test test_hfem_ads_zero_duration_cloud`

Expected: FAIL，原因是 HFEM ADS 服务模块不存在。

- [ ] **Step 3: 实现 options/stats 和建域**

```fortran
type,public :: hfem_ads_options_type
    integer :: coordinate_mode=ADS_COORD_COMPONENT
    real(DP) :: domain_sigma=3.0_DP
    real(DP) :: srp_sigma=0.0_DP
    integer :: da_order=4
    integer :: max_split_depth=8
    real(DP) :: error_tolerance(6)=[(1.0e-1_DP,i=1,3),(1.0e-6_DP,i=1,3)]
    real(DP) :: rel_tol=1.0e-12_DP,abs_tol=1.0e-12_DP
    real(DP) :: dt_min=1.0e-6_DP,dt_max=3600.0_DP
    integer :: max_steps=100000
end type
```

`stats` 至少含请求/输入/域内/域外/传播/写出数量、patch 数、BFS 次数、最大队列、深度受限 patch、`split_counts(:)` 和耗时。

- [ ] **Step 4: 实现可选 SRP patch-local 传播**

`options%srp_sigma==0` 时初始化 DACE 6 变量并始终 `clear_srp_scale_uncertainty`；大于零时初始化 7 变量，并对每个 patch 用 history center/width 计算 `eta_center`、`eta_span` 后调用 `set_srp_scale_uncertainty(7,...)`。无论成功失败都清除 SRP 模块状态。力模型由 `config%srp_model` 决定，不构造 8--11 维。

- [ ] **Step 5: 实现样本映射和统计计数**

把物理轨道偏差和可选 eta 映射到单位域；域外列标记并排除；调用 `mf_evaluate_points`；输出只包含有效列，且第七行复制输入 eta。零有效样本返回负状态。

- [ ] **Step 6: 运行零时长 GREEN**

Run: `fpm test test_hfem_ads_zero_duration_cloud`

Expected: PASS。

- [ ] **Step 7: 重构长传播测试并先观察失败**

删除测试内重复 BFS/evaluator，改调用生产接口；新增 cannonball/box-wing 7D 尺度响应断言，并确认 `dace_max_variables()==7`。

Run: `fpm test test_hfem_ads_propagation`

Expected: 在生产接口尚未补足长传播行为时 FAIL。

- [ ] **Step 8: 完成最小长传播实现并回归**

Run:

```bash
fpm test test_hfem_ads_zero_duration_cloud
fpm test test_hfem_ads_propagation
fpm test test_srp_da_parameter_map
```

Expected: 全部 PASS。

- [ ] **Step 9: 提交**

```bash
git add src/lib/uncertainty/propagation/pod_uq_hfem_ads_module.f90 test/test_hfem_ads_zero_duration_cloud.f90 test/test_hfem_ads_propagation.f90
git commit -m "feat: propagate HFEM uncertainty with ADS"
```

### Task 4: OOP 适配器与统一 UQ 调度

**Files:**
- Rewrite: `src/lib/uncertainty/propagation/pod_uq_prop_ads_module.f90`
- Modify: `src/functions/orbitprop/pod_uq_propagation.f90`
- Create: `test/test_uq_hfem_ads_dispatch.f90`

**Interfaces:**
- Produces: module `pod_uq_ads_module` and type `uq_ads_propagator`
- Produces: `METHOD_ADS=4`
- Extends: `run_uq_propagation(...,ads_options,ads_stats)`

- [ ] **Step 1: 写统一调度失败测试**

测试 `METHOD_ADS` 的零时长 6D 与可选 7D 调度，断言 `final_state%samples` 已分配、粒子数未丢失、均值/协方差已计算；6D 不含 eta，7D 含 eta。

- [ ] **Step 2: 运行 RED**

Run: `fpm test test_uq_hfem_ads_dispatch`

Expected: FAIL，原因是 `METHOD_ADS` 或新适配器不存在。

- [ ] **Step 3: 重写适配器**

删除 CRTBP、`mu`、`set_ads_mu` 和硬编码 `100000`；保存 `hfem_ads_options_type :: options` 与 `hfem_ads_stats_type :: last_stats`；`propagate` 将 `input_state%samples` 传给 HFEM 服务并调用 `output_state%compute_moments()`。

- [ ] **Step 4: 接入统一调度**

在 `run_uq_propagation` 增加 `METHOD_ADS=4` 和 ADS 分配分支。内部粒子生成仍从完整 OPM 协方差产生相关轨道粒子；只有 `ads_options%srp_sigma>0` 才把 state 分配为 7D 并独立生成 eta。为有限 ADS 域执行拒绝重采样，直到获得请求数，同时统计总抽样与拒绝数量。

- [ ] **Step 5: 运行 GREEN 与旧方法回归**

Run:

```bash
fpm test test_uq_hfem_ads_dispatch
fpm test test_uq
fpm test test_run_hfem_uprop_zero_cov_io
```

Expected: 全部 PASS。

- [ ] **Step 6: 提交**

```bash
git add src/lib/uncertainty/propagation/pod_uq_prop_ads_module.f90 src/functions/orbitprop/pod_uq_propagation.f90 test/test_uq_hfem_ads_dispatch.f90
git commit -m "feat: dispatch HFEM ADS through unified UQ"
```

### Task 5: 动态粒子 I/O 与 HFEM CLI

**Files:**
- Create: `src/lib/data/pod_uq_sample_io_module.f90`
- Modify: `app/run_HFEM_uprop.f90`
- Create: `test/test_uq_sample_io.f90`
- Create: `test/test_run_hfem_uprop_ads_io.f90`

**Interfaces:**
- Produces: `read_uq_perturbations_csv`
- Produces: `write_uq_particles_csv`
- Produces: `write_uq_moments_json`
- Produces: `write_ads_stats_json`

- [ ] **Step 1: 写动态 I/O 失败测试**

测试六列/七列扰动读取、混合列数拒绝、动态 CSV 表头、moments JSON 的 6D/7D 数组、stats JSON 粒子计数/坐标模式/基矩阵/分裂次数和 `sample_n_minus_1` 字段。

- [ ] **Step 2: 运行 RED**

Run: `fpm test test_uq_sample_io`

Expected: FAIL，原因是 I/O 模块不存在。

- [ ] **Step 3: 实现最小动态 I/O**

读取器只接受整文件统一 6 或 7 列有限实数；写入器按实际维数生成表头；JSON 使用动态循环，不用固定六元素格式。所有 open/read/write 错误返回 `status,message`。

- [ ] **Step 4: 运行 I/O GREEN**

Run: `fpm test test_uq_sample_io`

Expected: PASS。

- [ ] **Step 5: 写 CLI 失败测试**

以零时长/极短时长和小粒子数调用：

```text
fpm run run_HFEM_uprop -- -opm <fixture> -m ADS -dt 0 -n 8 -o /tmp/pod_ads_cli --srp-sigma 0
```

另跑 `--srp-sigma 0.02`；断言 6D/7D 表头、精确数据行数、moments 和 ads_stats 文件存在且计数一致。

- [ ] **Step 6: 运行 CLI RED**

Run: `fpm test test_run_hfem_uprop_ads_io`

Expected: FAIL，原因是 CLI 尚不接受 ADS。

- [ ] **Step 7: 接入 CLI**

增加 `-m ADS`、`--ads-coordinates`、`--ads-domain-sigma`、`--ads-max-depth`、`--ads-pos-tol`、`--ads-vel-tol`、`--srp-sigma`、`--perturb-file`。DACE 初始化变量数为 `merge(7,6,srp_sigma>0)`。ADS 输出固定为：

```text
<prefix>_particles.csv
<prefix>_moments.json
<prefix>_ads_stats.json
```

不改变 MC/DA/UT 和 `run_srp_uq_propagation` 的既有行为。

- [ ] **Step 8: 运行 CLI GREEN**

Run:

```bash
fpm test test_uq_sample_io
fpm test test_run_hfem_uprop_ads_io
fpm test test_run_hfem_uprop_zero_cov_io
fpm test test_run_srp_uq_propagation_io
```

Expected: 全部 PASS。

- [ ] **Step 9: 提交**

```bash
git add src/lib/data/pod_uq_sample_io_module.f90 app/run_HFEM_uprop.f90 test/test_uq_sample_io.f90 test/test_run_hfem_uprop_ads_io.f90
git commit -m "feat: expose ADS in HFEM CLI"
```

### Task 6: 使用文档与最终验证

**Files:**
- Create: `docs/hfem_ads_usage.md`
- Modify: `docs/superpowers/plans/2026-09-09-hfem-ads-integration.md`

**Interfaces:**
- Documents: 6D/7D 选择、component/whitened、截断 Gaussian、SRP sigma、输出文件和后台运行命令。

- [ ] **Step 1: 写使用文档**

至少给出 6D、7D component、6D whitened 和外部扰动文件四条命令；明确正确后台命令先创建目录，并把 stdout/stderr 与结果前缀分开：

```bash
mkdir -p output/ads_test
nohup fpm run run_HFEM_uprop -- <args> -o output/ads_test/hfem_propagation > output/ads_test/run.log 2>&1 &
```

- [ ] **Step 2: 执行完整定向回归**

Run:

```bash
fpm build
fpm test test_compiled_da_eval_into
fpm test test_compiled_da_batch_eval
fpm test test_srp_da_parameter_map
fpm test test_ads_core_dynamic_dimensions
fpm test test_ads_coordinate_mapping
fpm test test_ads_domain_scale
fpm test test_hfem_ads_zero_duration_cloud
fpm test test_hfem_ads_propagation
fpm test test_uq_hfem_ads_dispatch
fpm test test_uq_sample_io
fpm test test_run_hfem_uprop_ads_io
fpm test test_uq_crtbp_comparison
```

Expected: build 成功，全部测试 PASS。

- [ ] **Step 3: CLI 冒烟**

分别运行 `srp_sigma=0` 和正值的小粒子算例，核对三类输出文件、CSV 行数和 stats JSON 计数。

- [ ] **Step 4: 检查变更边界**

Run:

```bash
git diff --check
git status --short
git diff main...HEAD --stat
```

Expected: 无空白错误；所有提交只包含本计划文件；用户原有未提交修改仍存在且未进入 ADS 提交。

- [ ] **Step 5: 更新计划勾选并提交文档**

```bash
git add -f docs/hfem_ads_usage.md docs/superpowers/plans/2026-09-09-hfem-ads-integration.md
git commit -m "docs: explain HFEM ADS propagation"
```
