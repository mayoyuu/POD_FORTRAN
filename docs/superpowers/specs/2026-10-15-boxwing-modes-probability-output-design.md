# Box-Wing 8D 联合概率传播与逐小时分布输出设计

日期：2026-10-15。状态更新（2026-10-16）：本书面规格已获用户确认并授权实现。可选概率传播分支已接入现有 C2；已完成 50,000 点计算 / 2,000 点落盘的两小时 smoke 和输出审计。新的 360 h 概率算例尚未启动，不能用短例代替完整验收。见 `../reports/2026-10-16-boxwing-modes-probability-implementation-smoke.md`。

## 1. 目标与范围

在已经通过 15 天数值和收敛验证的 `DRO-2_init` C2 模式代理 ADS 上，增加可选的概率样本求值和分层输出：每小时用全部 50,000 点计算均值、完整 6×6 状态协方差等统计，常规点云只保存 2,000 个固定编号样本；关键时刻额外保存全量。8D 指六个初始状态变量和两个固定模式参数，实际输出轨道状态仍为六维。

参考现有设计和实现：

- `docs/superpowers/specs/2026-10-08-box-wing-initial-only-mc-design.md`：固定样本编号、ADS–MC、逐时统计、HDF5、验证与审计发布。
- `docs/superpowers/specs/2026-09-28-srp-cloud-statistics-design.md`：设计点不等于概率样本、完整快照决定能否离线重采样。
- `src/lib/statistics/pod_srp_mc_statistics_module.f90`：均值、完整协方差、标准差、偏度/峰度及范数分位数。
- `tools/pack_srp_mc_initial_points_hdf5.py`：流式严格校验、压缩、临时文件和原子发布。

当前旧输出仅保存 512 个确定性设计点，以及初/末完整多项式，无法从中恢复任意时刻的可靠 50,000 点概率分布。需要新增一次 C2 概率输出运行；不是把原点云重新标成 MC，也不是独立 Real 积分 50,000 条轨迹。

范围限于当前 `DRO-2_init`、原始完整 OPM 协方差、原模式表与动力学、360 h、每小时输出。**不自动扩展到 DRO-3，不重跑 C0/C1，也不重跑旧三组收敛实验。** 已发布基线、变体和 Real 真值保持只读。

## 2. 用户选择与概率测度

### 2.1 初始状态

设 OPM 均值为 `mu0`，正定协方差为 `P0=L L^T`，`L` 为与现有白化坐标一致的下三角 Cholesky 因子。

```text
z_i ~ N(0,1) conditioned on |z_i|<=3, i=1..6
z_1,...,z_6 mutually independent
q_i = z_i/3
x0 = mu0 + L z = mu0 + (3L) q_1:6
```

这是白化坐标逐分量截断，不是物理分量逐一裁剪，也不是三倍 Mahalanobis 半径内的球。不得把样本裁剪到边界；采用拒绝采样。

`P0` 是截断前母高斯的协方差。截断后的目标初始协方差为 `kappa*P0`，其中 `kappa=1-6*phi(3)/(2*Phi(3)-1)`，不能强行要求实际样本协方差精确等于 `P0`。

### 2.2 模式参数

按所选方案中的独立均匀先验确定：

```text
s ~ Uniform[0,1]
u ~ Uniform[0,1]
s independent of u and of z_1:6
q7=2s-1, q8=2u-1
xi_E=s, xi_M=(1-s)u
```

每个样本的 `s,u` 在整段传播中固定，不按小时重新抽取，不引入随机过程或过程噪声。此先验是明确的建模假设，不是由 P1 模式表或航天器物理自动推出的已知分布。

`s,u` 均匀**不等于** `(xi_E,xi_M)` 在混合系数三角形上均匀。不得改成三角形均匀采样或把两个 xi 都独立均匀。作为先验检查，`E[s]=E[u]=.5`、`Var(s)=Var(u)=1/12`，`E[xi_E]=.5`、`E[xi_M]=.25`。

概率结论仅对上述先验、±3σ 截断、该轨道和代理模型成立。原有“不赋予 s/u 概率测度”的结果继续解释为确定性包络，不追溯赋予新先验。

## 3. 实现路线与取舍

选择：**为现有 C2 runner 增加显式启用的概率输出分支**，复用同一条增量 ADS 历史，每小时批量求值同一组固定概率坐标。概率读写拆成小模块，原积分器、ADS 分裂核心、力公式不变。

采用“计算 50k、常规保存 2k、关键时刻保存 50k”的分层口径，沿用现有 ADS–MC 的统计/可视化/关键帧分离设计。相较原工程的 2048 默认值，本期按用户提出的 2000 设置，不能静默换成 2048。减少的是 I/O 和存储，不减少全部 50k 点的 ADS 求值、统计和覆盖检查。

另外两条路线不作为本期实现：

1. 只用旧快照离线求值：当前仅初末快照，不能满足逐小时输出。
2. 复制一个独立概率传播主程序：会重复 C2 的零模式、Real 探针、事件、快照和审计逻辑，容易产生两套数值行为。

默认关闭概率分支。未指定概率样本时，旧 C2 CLI/schema 和数据含义不变；旧 6D/7D runner 完全不改。概率开关不得反过来改变 ADS 分裂判据、步长或物理状态。

## 4. 固定样本

新增专用 8D 概率样本生成器，不改变旧 7D eta 样本格式。正式默认 `N=50000`，编号 `sample_id=1..N`，等权；前 10,000、20,000、50,000 为嵌套前缀。

常规保存点数 `V=min(2000,N)`，只生成一次固定编号列表，所有时刻不变。沿用已有按编号等距选点规则：V>1 时 `id_j=1+floor(j*(N-1)/(V-1)+.5)`，j=0..V-1；V=1 时取编号 1。选点只依赖编号，不依赖传播状态、误差大小或尾部位置，不能逐小时重新挑点。清单与其哈希写入 manifest。

统计总体的等权为 `1/N`，可视化子集若单独估计经验分布，权重应重新归一化为 `1/V`。该子集只是较低分辨率的分布样本，不能冒充完整 50k 人群，不能把它的均值/协方差覆盖正式全量统计。

使用 `numpy.random.PCG64DXSM`。主种子固定为 `20261015`，通过 `SeedSequence.spawn(2)` 派生初始状态和模式参数两个独立子流，记录 spawn key、NumPy 版本、算法和拒绝次数。改变 N 不得改变共同前缀。

样本 CSV 第一行为 `schema=srp_modes_probability_samples_v1`，第二行为字段头，字段固定为：

```text
sample_id,z1,z2,z3,z4,z5,z6,q1,q2,q3,q4,q5,q6,q7,q8,s,u,xi_E,xi_M
```

manifest 保存分布、相互独立假设、映射公式、截断规则、种子、样本数、前缀数和 CSV 哈希。各前缀使用归一化等权 `1/n`；不得将完整总体的 `1/50000` 权重直接用于较小前缀而不重归一化。

读取时严格检查 schema、顺序/唯一编号、全部有限值、z/q/s/u 的域和映射关系、xi 非负及其和不超过 1。禁止缺点补零、裁剪越界样本或悄悄重新抽样。

2026-10-15 实查主仓库 `20261005SRP/samples/` 为空，未找到旧 `mc_standard_normal_50000` 文件。因此本期生成新的有完整身份记录的样本，不假装复用了旧文件，也不声称与历史 MC 人群逐点配对；C0/C1 的历史结果仍保持原样。

## 5. 逐时求值、参考与验证

在 `0..360 h` 的全部 361 个时刻，对同一组全部 50,000 个 q 求值：

```text
x_i(t) = saved/current 8D ADS polynomial evaluated at q_i
e_i(t) = x_i(t) - x_CB_nominal(t)
```

共同参考是**当前这次 C2 内部**同配置 Cannonball Real 名义轨迹；不是 `q=0` 所对应的 `s=u=.5`，也不替换为均值或历史 C0 的近似匹配名义轨迹。

每个概率点必须找到所属 patch，采用现有全局到 patch-local 映射；分裂边界采用现有确定性归属约定，不能重复计权。内存中对全部 N 点检查归属；落盘的常规 V 点和关键帧 N 点各自保存 patch_id。覆盖完整性按全部 N 点检查，不只检查保存的少量点。

原有 512 个设计点、337 个 Real 探针、零模式控制、事件与数值门槛保持独立。新增从概率样本中等距选择固定 64 个编号（`round(1+i*(N-1)/63)`，i=0..63），独立 Real 累积传播并逐小时对比 ADS；小样本测试取 min(N,64) 并去重。

原探针和概率探针均执行位置 ≤0.1 km、速度 ≤1e-6 km/s；零模式仍为 1e-6 km / 1e-9 km/s。不降低现有门槛。概率探针保存独立的逐点文件和最大误差，不把验证点数 337 改写成 401 或把二者混为一组。

这些是经过独立 Real 抽样检查的 **ADS 概率映射样本**，不是 50,000 条独立高保真 Real 真值，也不是连续域精度证明。原有名义轨迹上的物理代理误差不能自动解释为整个概率样本域的代理误差上界。

## 6. 每小时统计

每个时刻分别对可用的 10k/20k/50k 前缀计算，正式运行应有 361×3=1083 行；较小 smoke 的前缀为标准前缀中小于 N 的项再加 N，必须如实记录实际数量。

这些统计直接使用内存中的全部对应前缀，不从 2,000 个落盘点反算。统计行明确区分 `statistics_sample_count` 与 `saved_visualization_count`。常规保存从 50k 降为 2k 不改变均值、完整协方差、偏度/峰度及分位数的计算精度；只降低今后仅依靠保存散点重新估计分布的分辨率。

保存：

- 当前名义参考六维状态、六维绝对样本均值和六维误差均值。
- 完整 6×6 协方差，所有 36 个元素，包含 3×3 的位置—速度交叉块。
- 六分量标准差、最小/最大值、偏度、超额峰度。
- 位置/速度各自的三条 RMS 主轴、范数 RMS。
- 位置/速度误差范数的 P50、P95、P99、P99.9、最大值及最坏样本号。
- 样本数、补丁与深度状态、原探针和概率探针误差、有效性和首次失败/未确认时刻。

复用 `compute_group_stats` 的统计定义：协方差分母 `n-1`（无偏样本协方差估计），分位数采用排序样本的线性插值；偏度、峰度保留现有有限样本修正。新文件中明确记录这些定义，不与旧 512 设计点按 N 计算的几何总体矩混用。

误差平移不改变协方差，因此绝对状态与 e 使用同一协方差。计算中心矩优先基于 e，避免大绝对坐标造成消减精度问题。零方差分量的标准化偏度/峰度标为不可计算，不用旧统计函数返回的 0 冒充已验证的高斯性。逐分量提供 `skewness_valid` / `kurtosis_valid`，分别要求有效方差及 n>2 / n>3；不可计算值在 CSV 留空、JSON 用 null 表示，不写 NaN/Infinity。

协方差 CSV 字段为 `cov_11...cov_66`，索引顺序严格为 `[dx,dy,dz,dvx,dvy,dvz]`；位置块单位 km²，交叉块 km²/s，速度块 km²/s²。位置/速度主轴分开计算，不能把未经尺度化的混合单位 6×6 矩阵拿来解释“最长物理轴”。

不得预设传播后仍为高斯，不输出未经依据拟合的“唯一高斯 PDF”。逐时 2k 点支持散点图、主要概率质量区的粗略边缘分布/经验 CDF/直方图/KDE，不保证复杂六维联合密度或稀有尾部的精度。2k 点在真实 P99/P99.9 以上平均仅约 20/2 点，因此正式分位数必须使用实时全量 50k 计算值。关键帧保留全量点；其他时刻需要更密集样本时可依靠逐小时完整 ADS 快照和完整初始样本定义重建，不必再次积分传播。

## 7. 样本量收敛与充分性

“5 万点”是默认计算规模，不是尾部概率精度的自动保证。逐小时保留各前缀结果，检查 20k→50k，10k 用于辅助趋势判断。

对每个物理位置/速度块分别计算：

1. 均值变化：全弧最大均值向量差，除以 50k 全弧最大 `sqrt(trace(C_block))`，要求 ≤5%。
2. 协方差变化：全弧最大块 Frobenius 差，除以 50k 对应块全弧最大 Frobenius 范数，要求 ≤5%。
3. 位置—速度交叉块变化：最大差除以 50k 全弧最大 `sqrt(||C_rr||F*||C_vv||F)`，要求 ≤5%。
4. P50/P95/P99 范数分位数变化：各指标全弧最大绝对差，除以 50k 该指标全弧最大值，要求 ≤5%。

分母为零且差也为零记为通过且说明退化；分母为零而差非零记为未确认，JSON 用 null 和原因表示，不写 NaN/Infinity。禁止对近零均值直接做相对百分比，或对混合单位全协方差直接做一个未尺度化范数。

P99.9 单独报告变化和有效样本数，不把 50 个量级的尾部样本自动当作可靠极罕见事件保证。上述前缀比较只是经验样本量稳定性检查，不是独立重复试验或置信区间证明。

保存独立状态：`numerical_validation_passed`、`sample_convergence_confirmed`、`probability_output_accepted`。后者只有完整性、数值门禁和上述主统计收敛全部通过时才能为 true。收敛不足时保留结果并明确未确认，不悄悄把 50k 称为“足够”，不自动启动更多长算。

本期正式验收固定 N=50000。小样本 smoke 不具备 20k/50k 前缀时标为 `not_applicable_smoke`，不能据此声称满足正式概率充分性。其他 N 必须记录实际前缀，不自动沿用本次 50k 验收结论。

## 8. 输出与 HDF5

概率分支新增文件，旧 C2 文件均保留：

| 文件 | 内容 |
| --- | --- |
| `probability_samples.csv` / `probability_samples.json` | 原始概率输入和完整先验/种子/哈希身份 |
| `c2_probability_hourly_stats.csv` | 每小时各前缀的全部统计与有效性 |
| `c2_probability_points.h5` | 常规逐时 V 点、关键帧 N 点、完整 N 点初始定义及全量统计，按组明确区分 |
| `c2_probability_probes.csv` | 固定概率验证样本的 ADS/独立 Real 状态和误差 |
| `c2_probability_convergence.json` | 各指标稳定性、失败原因和最终充分性状态 |
| `c2_report.json` 的可选 `probability_outputs` 块 | 与原 C2 数值/物理报告关联的概率身份与结果索引 |

HDF5 新 schema 为 `srp_modes_probability_hdf5_v1`，不能冒用原 visualization/keyframe schema。布局：

```text
/hour                         (T,)          float64
/epoch_et                     (T,)          float64
/reference_state              (T,6)         float64
/population/sample_id         (N,)          int64
/population/z_state           (N,6)         float64
/population/unit_coordinates  (N,8)         float64
/population/mode_parameters   (N,2)         float64, [s,u]
/population/physical_weights  (N,2)         float64, [xi_E,xi_M]
/population/weight            (N,)          float64, 1/N
/visualization/sample_id      (V,)          int64
/visualization/weight         (V,)          float64, 1/V
/visualization/state          (T,V,6)       float64
/visualization/error          (T,V,6)       float64
/visualization/patch_id       (T,V)         int32
/keyframes/hour_index         (F,)          int64, zero-based index into /hour
/keyframes/kind               (F,)          UTF-8, fixed/event labels
/keyframes/sample_id          (N,)          int64
/keyframes/weight             (N,)          float64, 1/N
/keyframes/state              (F,N,6)       float64
/keyframes/error              (F,N,6)       float64
/keyframes/patch_id           (F,N)         int32
/moments/sample_count         (K,)          int64, 10k/20k/50k
/moments/mean_state           (T,K,6)       float64
/moments/mean_error           (T,K,6)       float64
/moments/covariance           (T,K,6,6)     float64
/validity/ads_accuracy_valid  (T,)          bool
/validity/precision_confirmed (T,)          bool
```

属性记录轨道/历元、frame/center/time scale、状态和协方差单位表、先验及截断规则、参考中心定义、sample/OPM/config/mode/source/executable 哈希、schema、统计估计定义和软件依赖版本。另记录 `statistics_sample_count=50000`、`visualization_sample_count=2000` 和 `keyframe_sample_count=50000`；每组 point_role 必须明确，不能根据一个含糊的 sample_count 推断全时域保存了 50k 点。两类点群均验证 `state-error=reference_state`；同一关键时刻 visualization 中的点须与 keyframe 同编号状态一致。

沿用 gzip level 4 + shuffle、按小时分块，常规点向量 chunk 为 `(1,min(V,2048),6)`，关键帧为 `(1,min(N,2048),6)`；不把全时域全部点读入内存。

常规 361×2000 点的 state+error 两份 float64 场未压缩约 69.3 MB；六个固定全量关键帧另增约 28.8 MB，合计约 98.1 MB，相较逐小时全量约 1.73 GB 减少约 94.3%。这是点状态主体的十进制 MB 估计，不含 patch_id、完整初始样本、元数据、额外事件关键帧和快照；压缩率不预先保证。

固定全量关键帧为 0、24、72、168、240、360 h，并始终包含实际运行起/终点。第一次分裂、第一次深度受限、第一次验证失败的前一保存时刻和发生时刻额外保留全量；同一时刻只写一份，可累计多个事件标签。普通 hourly visualization 仍包含这些小时的固定 V 点。前一小时全量状态只在内存暂存，供事件触发保存，不事先把每小时 N 点写盘。

原生只将常规 V 点和关键帧 N 点写入两个专用 staging CSV，首行为 `schema=srp_modes_probability_points_v1`，第二行为 `hour,epoch_et,keyframe_kind,sample_id,patch_id,x_km,y_km,z_km,vx_kms,vy_kms,vz_kms,dx_km,dy_km,dz_km,dvx_kms,dvy_kms,dvz_kms`。普通流的 keyframe_kind 固定为 `hourly_visualization`，关键帧标签不用逗号分隔。Python 两遍流式校验/打包后按组写入同一 HDF5。不能先写每小时 50k 大 CSV，再事后删除 48k 来宣称减少了输出成本。

磁盘预检按预计常规点、固定关键帧、额外事件余量、临时 CSV 和逐小时快照综合计算，至少保留 2 GB 可用余量，运行中检查剩余空间。HDF5 完整审计并绑定哈希后，才可移除本次 staging 中可由 HDF5 完整恢复的点 CSV；可选保留 CSV 便于诊断。旧结果及非本次生成文件不得删除。

## 9. 逐小时完整快照

启用概率输出时，每小时保存完整 v2 8D patch/map/split/coeff 快照；包含所有 patch，而不是只保存抽样点。事件前后仍按现有规则记录，同一时刻去重。普通 C2 不启用概率输出时，快照规则不变。

每份快照在原 512 点上重建核验，同时用概率样本中的固定 64 点检查 snapshot/live ADS 逐点一致。其用途是为后续不同样本规模或先验的离线求值提供足够资料；本期不新建通用离线重采样产品，也不声称快照中的 ADS 多项式是独立 Real 真值。

## 10. 文件边界与兼容性

实现位于现有 `/home/songyu/fortran_test/POD_Fortran/.worktrees/boxwing-modes-p1`，保留未提交 P1–P4 改动，不创建新分支、不提交/合并/推送。

计划新增：

- `tools/generate_srp_modes_probability_samples.py`：独立先验采样和样本身份。
- `src/lib/uncertainty/io/pod_srp_modes_probability_io_module.f90`：严格样本读取、逐时统计/完整协方差及点记录；复用已有统计计算。
- `tools/pack_srp_modes_probability_hdf5.py`：新 HDF5 schema、流式打包与数值核验。
- 对应采样、原生 IO/统计、HDF5、runner 和负例测试。
- 专用最小 Python requirements：NumPy、h5py；不因本功能强制引入 SciPy/pandas。

修改范围仅为 `app/run_srp_modes_ads.f90` 与 `tools/run_srp_modes_ads.py` 的 opt-in 分支和相关测试。复用现有 history/evaluate/Real probes、统计模块和快照接口，不改其既有数学行为。

wrapper 增加 `--probability-samples`、`--probability-manifest`、`--probability-count`（默认 50000）、`--probability-visualization-count`（默认 2000）、`--probability-validation-count`（默认 64）及可选保留 staging CSV。样本/manifest 必须同时给出，启用时读取真实样本数而不是只相信 CLI。visualization_count 必须为正，实际 V=min(请求值,N) 并写入报告；不能据此改变 statistics_count。默认 C2 的原参数含义不变。

NumPy/h5py 只在概率路径按需导入，不能让原 stdlib-only C2 的帮助、preflight 或普通运行依赖新增包。当前 WSL 环境缺少这些依赖；实现期必须先在项目隔离环境验证依赖和打包短测试，不修改系统 Python，也不在依赖不可用时先启动 15 天长算。

## 11. 审计、失败与发布

继续采用独占输出锁、新 staging、冻结 executable、输入和源码运行前后哈希检查、审计后原子发布。新输出写到独立目录，例如 `SRP/boxwing_modes_DRO-2_init_c2_probability_<run-id>`；不覆盖原 `*_20261014` 或 `*_20261015`。

概率采样 manifest、生成器、原生模块、packer、依赖版本及所有最终文件进入 provenance。旧 P4 报告只对其冻结源码有效；新增功能造成源码身份变化是新的运行身份，不追改历史报告。

未找到 patch、非有限值、域/映射错误、遗漏/重复 sample_id、缺小时、哈希变动或打包错误：不得正式发布，保留诊断 staging。ADS/Real 精度失败或深度受限：标记首次失败起的区段不可用于可信概率统计，保留诊断数据，但本次正式 360 h 验收不通过。不得跳过失败样本后重新归一化以掩盖丢点。

HDF5 逐帧核对落盘点 CSV、实际点号/数量、量纲、参考平移；moments 组与全部逐时原生统计 CSV 一致。关键帧可从落盘完整 N 点独立重算所有前缀的均值、完整协方差及分位数并与实时统计比较；普通小时只有 V 点，不能声称据此重新计算并验证了 N 点统计。普通小时全量统计依赖经合成/独立 NumPy 测试核对的原生全量计算、覆盖计数及有效性检查，并明确记录审计层级；如需追查，可用完整逐时快照和输入样本重建全量。正定性检查使用 `D=diag(sqrt(diag(P0)))` 尺度化后的无量纲协方差，允许最小特征值不小于 `-1e-10*max(1,lambda_max)` 的舍入容差；不静默修正真实负特征值。

`formal_result=true` 继续表示输出完整且数值/身份审计通过；若样本量收敛未确认，必须同时 `probability_output_accepted=false`，不能凭 formal_result 把概率充分性判作已通过。

## 12. 验证与完成标准

按 TDD 先失败后实现；在现有 dirty worktree 内仅改本功能文件：

1. 先验/样本：固定种子可复现、前缀稳定、截断而非裁剪、全部映射成立；检查截断正态和 Uniform 的已知一/二阶矩。样本文件缺失时明确失败，不伪造历史复用。
2. 统计：用已知六维相关点云验证完整协方差，特别是非零位置—速度交叉项；与独立 NumPy 计算交叉核验 `ddof=1`、均值、各阶矩及分位数。零方差和单位混合有专门负例。
3. HDF5：检查全部时刻/点号/坐标/权重、两种点群 shape、有限值、状态—偏差—参考关系；确认常规固定 V 点、关键帧全量 N 点和全量统计没有混淆。用重复、缺帧、错单位、错误协方差、改变统计总体、篡改哈希等反例确认拒绝路径。
4. 原生短例：概率开/关不改变原 history/shape 的数值；在实际 8D 分裂 patch 上保持样本全局 s/u 和正确归属；原 512/337 检查保留，新增概率点独立 Real 检查通过；句柄清理一致。
5. 2 h smoke：每小时全部 50k 求值和统计，0/1/2 h 各有 2k 常规点，0/2 h 各有 50k 全量关键帧，1 h 不应被当作全量关键帧；三个前缀统计、完整 6×6 协方差及可解析 HDF5/快照齐全。关闭点输出或改变 V，不得改变同一 N 点的统计结果。通过原有门槛和新增审计，才进入长算。
6. 新的一次完整 DRO-2_init C2：361 个唯一小时、每小时用完整 50k 计算统计（共 18,050,000 次样本状态求值），但只保存 361×2000=722,000 条常规点状态；六个固定关键帧另有 300,000 条全量状态，事件帧按实际数量另计。逐小时完整快照和有效性记录齐全，统计中的 N 始终为 50k，不因 V=2k 缩减。所有数值门禁、重建和输入/源码身份检查通过。
7. 同配置新旧 C2 的共同名义轨迹与 512 点云做兼容性比较，位置最大差 ≤1e-6 km、速度 ≤1e-9 km/s；若不满足先诊断，不把新增概率支路当作改变原传播的理由。
8. 主概率统计的 20k→50k 稳定性通过才称默认样本规模满足本规格；P99.9 与其他尾部结论仍按有限样本限制说明。
9. 运行相关原生/Python 回归并如实区分本功能失败、历史缺夹具、依赖缺失和超时；不把原 P4 尚未全绿的全仓回归重标为通过。

本期完成的含义是：在明确的新先验下，获得一条经过验证、可重建和可审计的 15 天逐小时概率样本与矩输出；不表示任意先验、任意轨道或极罕见概率都已被验证。
