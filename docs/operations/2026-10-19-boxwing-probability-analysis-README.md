# Box-Wing SRP 概率结果离线分析

本包只读取已有结果，不重新传播轨道，也不修改容差、原始数据或正式验收状态。
默认分析原有 15 条轨道的 0–300 h；NRHO 未纳入。

桌面目录：
`C:\Users\13113\Desktop\科研\误差传播\SRP影响\20261005SRP\boxwing_modes_probability_20261019`

## Windows 运行

在本目录打开 PowerShell；如交付包已有 `.venv`，运行：

```powershell
.\.venv\Scripts\python.exe .\analyze_boxwing_probability.py
```

也可以自行安装 Python 3.11 或更新版本，并创建环境：

```powershell
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r .\requirements-srp-probability-analysis.txt
.\.venv\Scripts\python.exe .\analyze_boxwing_probability.py --max-hour 300
```

无需 WSL、Fortran、DACE 或 SPICE。脚本和 `srp_probability_analysis_data.py` 必须放在一起。
默认输入根目录是脚本所在目录，默认输出为该目录的 `analysis`。
可指定其他位置（不改变原始输入）：

```powershell
.\.venv\Scripts\python.exe .\analyze_boxwing_probability.py --root 'D:\概率结果' --output 'D:\概率分析' --max-hour 300
```

`--max-hour` 为非负整数，要求从 0 到该小时的每小时数据完整。必需数据缺失、重复、非有限或格式错误时程序报错并返回非零状态；请保留错误中的轨道及文件路径。
只有明确位于所选窗口之后、且是文件最后一条的损坏记录才允许作为警告保存。

## 输入与输出

输入根目录保留三个原名目录：

- `boxwing_modes_DRO-2_init_c2_probability_15d_20261016`
- `boxwing_modes_six_orbit_probability_20261016/results`
- `boxwing_modes_remaining10_probability_20261017/results`

其中 `.轨道名.staging-*` 是保留的失败或中断结果。程序只识别已知 15 条轨道，发现重复候选即停止。
读取 `c2_probability_hourly_stats.csv`、`c2_history.csv` 和可选的 `c2_report.json`。
点云优先读取 `c2_probability_points.h5` 的 `visualization/error` 单小时切片；无 HDF5 时流式读取 `c2_probability_visualization.csv`。
不加载 50,000 点全量关键帧。

`analysis/campaign_summary.csv` 汇总原始路径、已保存最后小时、所选窗口、连续通过前缀终点、首个无效小时及直接原因、原始报告状态与警告。
`analysis/provenance.json` 保存输入统计/历史/报告的 SHA-256、原始源码版本（报告缺失时为空）、单位、解释约束及软件版本。
完整复制审计以包内复制清单为准。

每条轨道的子目录包含：

- `summary.csv` 和 `hourly_diagnostics.csv`：完整所选窗口、各项数值误差、深度/精度/有效性位、连续通过掩码及失败标记。
- `moments.npz`：`hours`、`mean_state`、`mean_error`、完整 `covariance`（默认 301×6×6）、`std`、布尔 `recorded_prefix`、分量次序及单位。
- `numerical_errors.png`：全窗口三类位置/速度数值误差及各自容差线。
- `mean_std.png`、`norm_quantiles.png`：只画从 0 开始连续通过的前缀。
- `last_prefix_cloud.png`：最后通过小时的 2,000 个固定样本投影和 50,000 样本协方差的 1/3-sigma 轮廓。

## 科学解释与限制

主要统计选用每小时 50,000 样本；10,000/20,000 前缀仅用于一致性检查，不当作额外小时。
位置单位 km，速度单位 km/s；协方差位置块为 km²、速度块为 km²/s²、交叉块为 km²/s。
物理误差向量的参考是同一运行的 cannonball Real nominal。

原始/概率探针容差为 0.1 km 和 1e-6 km/s；zero-mode 容差为 1e-6 km 和 1e-9 km/s。
任一直接超限、深度受限或累计有效性/精度位失败都会截断连续前缀；之后即使数值或标志恢复，也不重新接续。
无法从当前行直接解释的累计失败标为 `unknown_recorded_failure`。
缺少最终报告标为 `interrupted_unfinalized`；保留 `formal_result` 原值，不凭有效性位认定正式成功。

这是已有数值检查的记录前缀，不是新的正式验收，也不是对整个不确定域的精度证明。
50,000 样本是 ADS 映射求值，不是 50,000 次独立 Real 积分。
ADS 数值误差比较的是同一个 SRP 代理模型下 ADS 与 Real 的结果，不是相对于完整实时 Box-Wing 模型的误差。
点云子集不用于重算全体统计；协方差轮廓不是非高斯分布的概率保证，p999 也不是罕见事件保证。
无效段保留在 CSV/NPZ 和数值误差图中，但不作为可信分布绘制。
