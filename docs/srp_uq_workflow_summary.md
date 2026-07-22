# SRP 参数误差传播实验工作流程总结

## 1. 研究目标

本轮工作的目标是评估太阳光压（Solar Radiation Pressure, SRP）模型参数误差对 HALO 轨道误差传播结果的影响。前期 10 day、`srp_sigma=0.05` 的结果显示 SRP 参数误差对最终状态分布的影响较小，因此进一步设计了多传播时间、多 SRP 不确定性幅值的批量实验，并将结果同时放在惯性状态空间（RV）和站心光学观测空间（RADEC）下比较。

实验只关注误差传播，不涉及轨道改进或滤波更新。对照组采用“有标称 SRP、但不把 SRP 参数作为不确定性状态”的 6 维传播结果；实验组采用“状态 + SRP 尺度误差”的 7 维 DA-MC 传播结果。

## 2. 光压模型和误差参数设计

SRP 采用炮弹球模型。基本加速度尺度由反射系数、面积质量比和太阳辐射压共同决定：

```text
a_srp ∝ Cr * (A / m) * P_sun * (AU / r_sun)^2
```

当前 SRP-UQ 入口为 `app/run_srp_uq_propagation.f90`，默认参数为：

```text
Cr       = 1.25
A / m    = 7.5e-3 m^2/kg
P_sun    = 1367 / c N/m^2
eta_mean = 0.0
```

SRP 误差参数记为 `eta_srp`。它不是一个绝对加速度误差，而是 SRP 加速度的乘性尺度误差：

```text
a_srp_with_error = a_srp_nominal * (1 + eta_srp)
eta_srp ~ N(eta_mean, eta_sigma^2)
```

因此：

```text
srp_sigma = 0.01  ->  约 1% SRP 尺度误差
srp_sigma = 0.05  ->  约 5% SRP 尺度误差
srp_sigma = 0.10  ->  约 10% SRP 尺度误差
srp_sigma = 0.20  ->  约 20% SRP 尺度误差
```

在 DA-MC 传播中，`eta_srp` 被作为第 7 维不确定性状态加入初始样本：

```text
[x, y, z, vx, vy, vz, eta_srp]
```

输出的 SRP-UQ 粒子文件表头为：

```text
x,y,z,vx,vy,vz,eta_srp
```

nominal 对照仍为 6 维状态传播，不输出 `eta_srp`。

## 3. 实验矩阵

为避免 `init2` 和 `init` 混用，本轮只使用每个 HALO 轨道文件夹中的 `*_init.opm.json`：

```text
OPM/L1Halo-1/L1Halo-1_init.opm.json
OPM/L1Halo-2/L1Halo-2_init.opm.json
OPM/L1Halo-3/L1Halo-3_init.opm.json
OPM/L2Halo-1/L2Halo-1_init.opm.json
OPM/L2Halo-2/L2Halo-2_init.opm.json
OPM/L2Halo-3/L2Halo-3_init.opm.json
```

传播时间设置为：

```text
5 day   = 432000 s
10 day  = 864000 s
15 day  = 1296000 s
20 day  = 1728000 s
```

SRP 尺度误差标准差设置为：

```text
eta_sigma = 0.01, 0.05, 0.10, 0.20
```

因此实验规模为：

```text
SRP-UQ:    6 orbits * 4 propagation times * 4 sigma values = 96 cases
Nominal:   6 orbits * 4 propagation times                  = 24 cases
RADEC:     96 + 24 = 120 converted particle files
```

默认粒子数和 DA 阶数为：

```text
PARTICLES = 100000
DA_ORDER  = 4
```

## 4. 传播和转换流程

### 4.1 SRP-UQ 传播

脚本：

```text
scripts/run_halo_srp_sweep_parallel.sh
```

该脚本调用：

```bash
fpm run run_srp_uq_propagation -- \
  -cfg config/config.txt \
  -opm <opm> \
  -dt <dt_seconds> \
  -o <output_base> \
  -n <particles> \
  -da <da_order> \
  -srp-mean 0.0 \
  -srp-sigma <sigma>
```

默认输出目录：

```text
output/2026_07_for_wy_srp_sweep
```

输出文件名前缀示例：

```text
L1Halo-1_init_srpMean_0p0_srpSigma_0p05_dt_432000_n100000_o4
```

每个 case 生成：

```text
*_particles.csv
*_moments.json
```

### 4.2 nominal SRP 对照传播

脚本：

```text
scripts/run_halo_nominal_srp_sweep_parallel.sh
```

该脚本复用已有 6 维传播入口：

```bash
fpm run run_HFEM_uprop -- \
  -cfg config/config.txt \
  -opm <opm> \
  -m DA \
  -dt <dt_seconds> \
  -o <output_base> \
  -n <particles> \
  -da <da_order>
```

默认输出目录：

```text
output/2026_07_for_wy_nominal_srp_sweep
```

输出文件名前缀示例：

```text
L1Halo-1_init_nominalSrp_dt_432000_n100000_o4
```

这组结果包含标称 SRP 力模型，但不传播 `eta_srp` 参数误差，是 SRP-UQ 结果的直接对照。

### 4.3 RADEC 坐标转换

脚本：

```text
scripts/run_halo_radec_batch.sh
```

该脚本复用已有程序：

```text
app/run_radec_convert.f90
```

默认输入根目录为：

```text
output/2026_07_for_wy_srp_sweep
output/2026_07_for_wy_nominal_srp_sweep
```

因此它会同时转换 SRP-UQ 和 nominal 两类 `*_particles.csv`。默认测站为 `R91`，默认输出角度单位为 degree。输出文件追加：

```text
_radec_R91_deg.csv
```

例如：

```text
L1Halo-1_init_srpMean_0p0_srpSigma_0p05_dt_432000_n100000_o4_radec_R91_deg.csv
L1Halo-1_init_nominalSrp_dt_432000_n100000_o4_radec_R91_deg.csv
```

转换时，脚本从对应 OPM 文件读取初始 `EPOCH`，再根据文件名中的 `dt_seconds` 计算最终历元：

```text
epoch_final = OPM_EPOCH + dt_seconds
```

随后将该历元传给 `run_radec_convert` 的 `-et` 参数。SRP-UQ 的粒子文件虽然有 7 列，但 `run_radec_convert` 只读取前 6 个状态量，因此可以直接处理。

## 5. 结果分析方法

分析脚本：

```text
analysis/analyze_srp_sweep_results.py
```

该脚本将同一轨道、同一传播时间下的 SRP-UQ 结果与 nominal 对照逐 case 配对比较。分析同时在两个空间中进行：

```text
RV:    x, y, z, vx, vy, vz
RADEC: ra, dec
```

核心量化指标为前 4 阶矩的最大相对误差：

```text
1阶矩: mean
2阶矩: standard deviation
3阶矩: skewness
4阶矩: kurtosis
```

对每个 case，计算：

```text
max_relative_error = max_i |moment_srp(i) - moment_nominal(i)| / max(|moment_nominal(i)|, eps)
```

并额外输出以下辅助指标：

```text
RV:
  - 位置协方差 trace 的相对变化
  - 速度协方差 trace 的相对变化
  - 位置均值偏移范数
  - 速度均值偏移范数

RADEC:
  - RADEC 协方差 trace 的相对变化
  - RADEC 协方差面积 sqrt(det(cov)) 的相对变化
  - RADEC 均值偏移范数
```

分析输出包括：

```text
metrics_by_case.csv
aggregate_by_dt_sigma.csv
plots/
```

其中：

```text
metrics_by_case.csv
```

保存每一个 `orbit / dt / sigma` 的详细指标；

```text
aggregate_by_dt_sigma.csv
```

将同一 `dt / sigma` 下 6 个 HALO 轨道的指标做平均和最大值汇总，用来观察 SRP 影响随传播时间和不确定性大小的变化趋势。

## 6. 图像展示设计

分析脚本会生成两类图。

第一类是量化指标热力图，用于展示趋势：

```text
heatmap_rv_moments_1to4_max_rel_err.png
heatmap_rv_pos_cov_trace_rel_delta.png
heatmap_radec_moments_1to4_max_rel_err.png
heatmap_radec_cov_area_rel_delta.png
```

这些图的横轴是传播时间，纵轴是 `eta_srp` 的标准差。颜色表示相应指标在 6 个轨道上的最大值或平均值。它们回答的问题是：

```text
SRP 参数误差的影响是否随时间增长？
SRP 参数误差的影响是否随 sigma 增大而增强？
这种影响在 RV 空间和 RADEC 空间是否一致？
```

第二类是分布对比图，用于展示具体 case 的粒子分布：

```text
distribution_rv_selected_case.png
distribution_radec_selected_case.png
distribution_grid_rv_position_<orbit>.png
distribution_grid_radec_<orbit>.png
```

这些图将 SRP-UQ 分布与 nominal 分布画在一起，直观比较分布宽度、偏移和形状变化。

## 7. 当前结果入口

本轮传播和 RADEC 转换的核心结果位于：

```text
output/2026_07_for_wy_srp_sweep
output/2026_07_for_wy_nominal_srp_sweep
output/2026_07_for_wy_radec_batch
```

已验证的队列规模为：

```text
SRP-UQ queue:    96 cases
Nominal queue:   24 cases
RADEC queue:    120 cases
```

分析脚本已在一个小范围 case 上验证通过：

```text
orbit      = L1Halo-1
dt         = 432000 s = 5 day
srp_sigma  = 0.05
```

该 case 的样例指标显示：

```text
RV 前 4 阶矩最大相对误差      ≈ 6.80e-3
RV 位置协方差 trace 相对变化  ≈ 6.50e-6
RADEC 前 4 阶矩最大相对误差   ≈ 1.43e-3
RADEC 协方差面积相对变化      ≈ -1.08e-5
```

这说明在 5 day、5% SRP 尺度误差下，SRP 参数误差对协方差尺度的影响仍然很小，但高阶分布形状已经能在前 4 阶矩指标中被捕捉到。完整结论需要基于 `aggregate_by_dt_sigma.csv` 和热力图观察所有 `dt / sigma` 组合；预期如果 SRP 影响确实存在，应在更长传播时间和更大 `srp_sigma` 下呈现更清楚的增长趋势。

完整分析命令为：

```bash
cd /home/songyu/fortran_test/POD_Fortran
python3 analysis/analyze_srp_sweep_results.py \
  --srp-root output/2026_07_for_wy_srp_sweep \
  --nominal-root output/2026_07_for_wy_nominal_srp_sweep \
  --out-dir output/2026_07_for_wy_analysis
```

如果只想先分析一个轨道或一个网格点，可以使用：

```bash
python3 analysis/analyze_srp_sweep_results.py \
  --case-orbit L1Halo-1 \
  --dt-list 432000 \
  --sigma-list 0.05 \
  --out-dir output/2026_07_for_wy_analysis_check_small \
  --no-plots
```

注意：画图需要 Python 环境中安装 `numpy` 和 `matplotlib`。如果只生成 CSV 指标，可以使用 `--no-plots`。

## 8. 结论性说明

目前的 SRP 参数误差建模方式是合理的乘性尺度误差设计，而不是把一个无量纲误差错误地当成绝对力误差。早期结果中 SRP 影响很小，更可能来自以下原因：

```text
1. 传播时间较短，SRP 小扰动尚未充分积累；
2. HALO 轨道动力学中主导误差仍由初始状态协方差贡献；
3. SRP 加速度本身相对主要引力项较小；
4. RADEC 投影下某些方向的状态差异不敏感；
5. sigma=0.05 对当前炮弹球参数而言仍属于温和模型误差。
```

因此，本轮多时间、多 `srp_sigma` 实验的核心作用，是把“SRP 影响小”从单一 case 观察扩展为系统性结论：如果随着 `dt` 和 `srp_sigma` 增大，RV 和 RADEC 指标均持续增加，则说明 SRP 参数误差的影响存在但在原始设置下较弱；如果即使在 20 day、`srp_sigma=0.20` 下仍不明显，则可进一步判断当前任务场景对 SRP 参数误差不敏感，或需要重新审查面积质量比、反射系数、阴影模型等 SRP 物理参数设置。
