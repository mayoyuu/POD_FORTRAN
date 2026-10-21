# Box-Wing 8D 概率输出：代码实现与两小时 smoke

日期：2026-10-16。范围：DRO-2_init、现有 C2 模式代理 ADS，新增显式启用的概率分支。**本报告不是新的 15 天概率结果报告；360 h 概率运行尚未启动。**

## 已实现

- 同一批固定 8D 样本逐时 ADS 求值。正式默认 N=50000，六个独立 ±3 截断白化正态，加独立均匀 s/u；参数整段固定，不把 50k ADS 样本称为 50k Real 真值轨道。
- 每小时全体 N 点及 10k/20k/50k 嵌套前缀计算均值、完整 6×6 协方差（ddof=1，含位置—速度交叉项）、各阶矩、主轴、范数分位数与有效性。统计计算不使用 V 点可视化子集。
- 每小时仅保存固定 V=min(2000,N) 个编号；全量关键帧为固定日期、起终点及首次分裂/深度受限/验证失败与前一小时。前一帧只在内存缓存；不先写全量逐时大 CSV。
- gzip4/shuffle 分组 HDF5，分别声明 N/V、角色、权重和单位。保留输入总体及逐小时完整 v2 ADS 快照。
- 原 512 个设计点、337 个 Real 探针和零模式门槛不变；独立增加 64 个概率 Real 探针。每小时快照还核对原 512 点与固定概率探针。
- 每个时刻/前缀的矩一致性与退化标志校验；全量关键帧另以 NumPy 独立重算。普通小时不声称能用 2k 点重算 50k 统计。
- 样本/schema/映射/哈希/点号/时间/补丁/PSD 审计，源码和输入前后冻结身份、磁盘预检与运行中 2 GB 余量检查、审计后原子发布。旧 C2 的普通路径仍仅依赖 stdlib。

代码在 `/home/songyu/fortran_test/POD_Fortran/.worktrees/boxwing-modes-p1`，保留既有 dirty P1–P4 改动，无提交、合并或推送。新增独立采样器、原生 IO/统计模块、HDF5 打包审计器、NumPy/h5py 最小依赖及测试；只在现有 C2 app/wrapper 接入可选分支，不修改积分器或 ADS 分裂数学。

## 实际 smoke 结果

输出：`/home/songyu/fortran_test/POD_Fortran/SRP/boxwing_modes_DRO-2_init_c2_probability_smoke_20261016`。

| 项目 | 已核验结果 |
| --- | --- |
| 时间 | 0、1、2 h，共 3 个唯一时刻 |
| 每小时计算总体 | 50000；共 150000 次概率样本状态求值 |
| 常规保存 | 3×2000 点；HDF5 state shape `(3,2000,6)` |
| 全量关键帧 | 0、2 h；state shape `(2,50000,6)`；1 h 无全量关键帧 |
| 统计行 | 3 小时×3 前缀=9 行；covariance shape `(3,3,6,6)` |
| 原检查 | 512 个设计点、337 个 Real 探针/小时 |
| 概率验证 | 64 个独立 Real 探针/小时；快照概率检查共 192 次 |
| 原 ADS/Real 最大误差 | 4.247564216879669e-10 km / 1.8742721836921976e-15 km/s |
| 概率 ADS/Real 最大误差 | 2.982255177909857e-10 km / 1.5827980514452614e-15 km/s |
| 零模式最大误差 | 8.231806349783991e-11 km / 4.710277376051325e-16 km/s |
| 精度未确认 | false；无分裂、无深度受限 |
| HDF5 实际大小 | 15780213 bytes，约 15.78 MB |
| 本结果目录文件总计 | 39822362 bytes，约 39.82 MB（含一份冻结样本和 executable） |

原生输出与 HDF5 的逐帧数据核验、关键帧全部前缀独立重算、全部文件哈希、当前源码身份与输入身份均通过。两份点 staging CSV 仅在已验证 HDF5 后删除；点数据保存在 HDF5，可恢复。未删除任何既有结果。

同配置历史结果 `SRP/boxwing_modes_DRO-2_init_c2_15d_20261014` 只读复用前 0/1/2 h：共同名义轨迹和 1536 条设计点记录的最大位置/速度差全部为 **0**，没有重跑基线。另一个 N=32 短例中，概率开/关的 history/shape 逐字节一致；同一 N 将 V 从 5 改为 11，协方差逐字节一致。

报告身份：

```text
c2_report.json SHA256 = ea3ab7b9ad845cb0a2545e8181fbe03cc246539b0ce9631010df47a00af5fc04
c2_probability_points.h5 SHA256 = 178f2767ffd6c48c39e42e8d41a049dcfc9d2793ee197781ece52c7b169184a9
probability_samples.csv SHA256 = dc39316538c76bcc5f290a90b71f514b461d132c4840543a5452d0eb991d456a
```

`formal_result=true` 表示这次 2 h 结果通过数值与身份审计。`sample_convergence_confirmed=true` / `probability_output_accepted=true` 仅表示这次所计算 **0–2 h** 上的 20k→50k 前缀稳定性通过，不表示 360 h 已验证，也不保证 P99.9 或更罕见尾部概率。小于 50k 的开发 smoke 明确标为不适用正式样本充分性。

## 验证与审查

按 TDD 实施，先观察采样器/新模块/打包器/接口缺失的失败，再实现；另以失败用例驱动修正原生 CSV 多余分隔符、近零方差高阶矩标记、普通小时矩一致性、发布后输入身份路径。

只读代码审查确认全量统计与少量保存的边界、前缀、事件缓存、独立 Real 探针及快照逻辑；提出普通小时统计审计不足，已补充每行代数/标记/范围校验，并加入非关键小时篡改回归。审查复核后无未处理 Critical/Important 项。

隔离环境：`.venv-probability`，NumPy 2.5.4，h5py 3.16.0；未安装进系统 Python。`python3 -S tools/run_srp_modes_ads.py --help` 通过，证明普通帮助路径不需要新包。

主要验证命令（工作目录为上述 worktree；先 `source setup_env.sh`）：

最终冻结源码版本上的结果：15 项新增概率 Python 测试、6 项原 wrapper 测试、4 项原 native runner Python 测试，共 25 项通过；下列 7 个原生 target 均通过；`git diff --check` 通过。

```bash
.venv-probability/bin/python -m unittest discover -s test -p 'test_srp_modes_probability_*.py' -v
.venv-probability/bin/python -m unittest discover -s test -p test_run_srp_modes_ads.py -v
.venv-probability/bin/python -m unittest discover -s test -p test_srp_modes_ads_runner.py -v
fpm test --target test_srp_modes_probability_io
fpm test --target test_srp_modes_snapshot
fpm test --target test_srp_modes_study
fpm test --target test_hfem_ads_modes_8d
fpm test --target test_srp_modes_force_consistency
fpm test --target test_srp_mode_parameter_map
fpm test --target test_srp_mode_table
git diff --check
```

新输出测试包含合成相关点云、完整交叉协方差、零/近零方差、事件前帧标签合并，以及实际 8D 四补丁（第 7/8 维分裂）的全局 s/u、补丁归属、批量求值和 DA 句柄清理。原统计模块的 array-temporary 调试警告仍存在；不代表计算失败，本期没有为消除此警告而改旧统计数学。

本功能的定向回归不等于原 P4 全仓回归转绿；历史缺夹具/超时等未在本期重新标为通过。

## 后续完整运行

复用已经生成并核验的 50k 样本：

```bash
cd /home/songyu/fortran_test/POD_Fortran/.worktrees/boxwing-modes-p1
source setup_env.sh
.venv-probability/bin/python tools/run_srp_modes_ads.py \
  --opm OPM/DRO-2/DRO-2_init.opm.json \
  --config config/config_srp_modes_DRO-2.txt \
  --output-dir /home/songyu/fortran_test/POD_Fortran/SRP/boxwing_modes_DRO-2_init_c2_probability_15d_20261016 \
  --hours 360 --save-hours 1 \
  --probability-samples /home/songyu/fortran_test/POD_Fortran/SRP/boxwing_modes_DRO-2_init_probability_samples_20261016/probability_samples.csv \
  --probability-manifest /home/songyu/fortran_test/POD_Fortran/SRP/boxwing_modes_DRO-2_init_probability_samples_20261016/probability_samples.json \
  --probability-count 50000 --probability-visualization-count 2000 --probability-validation-count 64
```

该命令尚未执行。完整验收还需要 361 小时、1083 行统计、361×2000 常规点、六个固定全量关键帧（另加实际事件帧）、全部快照与门禁、整弧 20k→50k 稳定性。P99.9 始终单列有限样本限制；不扩展到 DRO-3，不重跑 C0/C1。
