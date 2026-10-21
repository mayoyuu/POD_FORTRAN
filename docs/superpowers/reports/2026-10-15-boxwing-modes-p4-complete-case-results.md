# P4：DRO-2_init 完整 15 天案例核验

核验日期：2026-10-15。结论：**一个完整 360 h 的 C2 轨道传播案例已经跑通，三组全弧收敛检查通过，独立 Real 物理比较完成。全仓回归尚未全绿，因此不宣称设计 §17 的全部完成标准已满足。**

遵循用户要求，复用原 C0/C1 历史结果，不重新传播 6D/7D 对照域；未重跑已完成 C2 基线。源码保留在原隔离工作区，未提交、合并或推送；既有结果未覆盖。

## 1. 完整性与正式结果

所有 C2 运行均为 8 个活动变量、DA/ADS，domain_sigma=3、最大深度 8；512 个相同研究点、337 个相同独立 Real 探针；每小时保存，覆盖 0…360 h，共 361 个时刻。起始/终止历元为 `852033729.1839125` / `853329729.1839125` TDB s，J2000、EARTH，full-sun。s/u 是有界确定性参数，不赋予概率测度。

正式目录位于主仓库 `/home/songyu/fortran_test/POD_Fortran/SRP/`：

| 运行 | 正式目录名 | DA 阶 | 模式表间隔 | RKF78 rel/abs |
| --- | --- | --- | --- | --- |
| 原基线，直接复用 | `boxwing_modes_DRO-2_init_c2_15d_20261014` | 4 | 1 h | 1e-12 |
| 模式表加密 | `boxwing_modes_DRO-2_init_c2_half_hour_20261015` | 4 | 0.5 h | 1e-12 |
| 积分容差收紧 | `boxwing_modes_DRO-2_init_c2_tight_rkf_20261015` | 4 | 1 h | 1e-13 |
| DA 升阶 | `boxwing_modes_DRO-2_init_c2_da5_20261015` | 5 | 1 h | 1e-12 |

四份 `c2_report.json` 均为 `formal_result=true`。每份实际 361×512 点云、361×337 探针记录已重新读取，检查唯一性、连续时刻、有限数值、实际状态差及报告一致性。全部 artifact/input 哈希重核通过。四组均为 1 patch、没有自然分裂、未达到深度上限，`validation_failed=false`、`precision_unconfirmed=false`，DA 句柄 0→0。分裂行为由已有专门测试覆盖，不从本案例无分裂推导分裂路径已被自然激活。

汇总分析结果：[case_report.json](/home/songyu/fortran_test/POD_Fortran/SRP/boxwing_modes_DRO-2_init_p4_analysis_20261015/case_report.json)；逐小时可检查指标：[hourly_metrics.csv](/home/songyu/fortran_test/POD_Fortran/SRP/boxwing_modes_DRO-2_init_p4_analysis_20261015/hourly_metrics.csv)。

## 2. 数值门禁

门禁未降低：ADS 对同代理 Real 最大位置/速度误差 ≤0.1 km / 1e-6 km/s；零模式 ADS 对 Cannonball Real 名义控制 ≤1e-6 km / 1e-9 km/s。

| 运行 | ADS 位置 km | ADS 速度 km/s | 零模式位置 km | 零模式速度 km/s |
| --- | ---: | ---: | ---: | ---: |
| 原基线 | 1.863267429652309e-5 | 2.3448859556127377e-10 | 1.2523157749779443e-7 | 1.2160771899899936e-12 |
| 0.5 h 表 | 1.8635681765278537e-5 | 2.34493571576413e-10 | 1.1782729162669417e-7 | 1.2074073343769302e-12 |
| RKF 1e-13 | 1.863267429652309e-5 | 2.3448859556127377e-10 | 1.2523157749779443e-7 | 1.2160771899899936e-12 |
| DA 5 | 2.4769989226455524e-7 | 3.4955820411178488e-12 | 5.1314289318810894e-8 | 3.2925004047334924e-13 |

基线最大 ADS 位置误差约 1.86 cm。所有运行通过原数值门禁。基线零模式独立 Real 对 Cannonball Real 的控制误差为 0。

四组初/末快照均在内存中重新求值于保存的全部 512 个坐标，没有改写原快照、重建 CSV 或原审计报告。最大重建位置差 `5.820766091346741e-11 km`，最大速度差 `5.551115123125783e-17 km/s`，通过 1e-8 km / 1e-11 km/s 的重建门槛。

基线末时刻 ADS 对新独立 proxy Real 的位置/速度误差：Earth `2.740205621635553e-8 km / 1.5735428708053844e-14 km/s`，Moon `1.236350080809748e-7 km / 1.3862745025930546e-12 km/s`。

## 3. 收敛性

使用同一 512 点网格，将每组保存偏差还原为绝对状态，再减去共同的历史 C0 名义轨迹。主指标为每小时位置/速度采样包络最大范数；相对变化定义为 `max_hour(abs(variant_metric-baseline_metric)) / max_hour(abs(baseline_metric))`。同时检查 Earth/Moon 末时刻条件 SRP 位移相对变化。两类均须 ≤5%；不对接近 0 的数值残差计算误导性百分比。

| 变体 | 位置包络相对变化 | 速度包络相对变化 | Earth/Moon 位移相对变化 | 结论 |
| --- | ---: | ---: | --- | --- |
| 0.5 h 表 | 7.801362054403109e-7 | 6.897182842361699e-7 | 3.554400426222593e-6 / 2.429201794102125e-5 | 通过 |
| RKF 1e-13 | 0 | 0 | 0 / 0 | 通过 |
| DA 5 | 3.885553323013091e-9 | 1.695518609037026e-8 | 0 / 1.564883225144107e-9 | 通过 |

表中是无量纲比例，不是百分数。最大片刻/研究点状态差分别为：0.5 h 表 `9.363543721171273e-4 km / 2.8407963909431817e-9 km/s`；RKF 加严 0；DA 5 `1.8587671414504303e-5 km / 2.334499077079892e-10 km/s`。加严 RKF 在本案例产生完全相同的保存状态，不能据此声称所有轨道都与积分容差无关。

## 4. 物理比较：三种误差/效应分开报告

独立 Real-only 工具采用原 Real 积分器与物理力模型，不构建 DA 域。五条名义轨迹：Deta η=+.3、proxy Earth (s=1,u=.5)、proxy Moon (s=0,u=1)、完整 Box-Wing Earth/Moon。每条均实际到达 360 h，共 1805 条记录；rel/abs=1e-13，几何 fallback 次数 0，DA 句柄 0→0。

正式结果：[truth_report.json](/home/songyu/fortran_test/POD_Fortran/SRP/boxwing_modes_DRO-2_init_truth_20261015/truth_report.json)。

| 物理指标 | Earth | Moon |
| --- | ---: | ---: |
| proxy Real − 完整 Box-Wing Real 全弧最大位置差 km | 0.00029510227003063196 | 0.0006986788344201328 |
| 同上最大速度差 km/s | 1.3310323333614272e-9 | 4.1978687232857605e-9 |
| 完整 Box-Wing − C0 名义轨迹末位置位移 km | 17.063230563825265 | 37.95227139032189 |
| proxy − C0 名义轨迹末位置位移 km | 17.062937241268198 | 37.95292844681379 |

代理模型位置差约 0.30 m / 0.70 m，最大值均出现在 360 h。它是沿这两条名义轨迹量化的模型近似误差，不是 ADS 数值误差，也不是全不确定域代理误差上界。17/38 km 则是物理模型相对炮弹球的条件轨道位移，不是算法失败。Deta η=+.3 的末位置位移为 `45.12815041026905 km`。

## 5. C0/C1 复用范围

C0：复用 `SRP/SRP_UCERTAINTY_PROP_WITHINIT_ETA0/DRO-2`。OPM 哈希、原配置哈希、起始历元、初始均值与 6×6 状态基、非 SRP 物理配置、Cr/面积/质量相符；361 个小时的数值/精度状态通过。

不得隐藏的小差异：C0 默认光压为 `1367/c=4.559821181358738e-6 N/m²`，C2 显式值为 `4.559821181419152e-6`，相对差 `1.3249e-11`。实际全弧名义轨迹最大差为 `9.808390277792056e-9 km / 2.753873468939917e-14 km/s`。因此 C0 是已记录差异的共同参考，不标作重新运行的完全同配置对照。

C1：现有 `SRP/SRP_UCERTAINTY_PROP_WITHINIT/DRO-2` 的物理域是 η±.9，必须称 C1-H。其旧报告缺少 OPM/config/executable 身份哈希。仅按 `q7_new/3` 对末时刻多项式做 η±.3 子域重建，不声称获得全新、完全溯源匹配的 C1-30 或其逐小时多项式。该限制已写入机器可读结果。子域零状态 η+.3 终点对独立 Deta Real 的差为 `5.725767217017969e-9 km / 5.1791785103661725e-14 km/s`，支持该端点的一致性，但不补全缺失的历史身份信息。

共同初始状态方向下，末时刻位置/速度包络：C0 `1159.2410572572726 km / .0035641132893550493 km/s`；历史 C1-H 的 ±.3 子域 `1201.2768345297486 / .003686039882178524`；完整 C1-H `1285.9629225981728 / .003936928554479495`。基线 C2 原保存偏差包络 `1195.9532322727653 / .0036720744098945675`。这些是采样盒域包络，不是连续域证明或概率置信区间。

详见隔离工作区的 [复用审计](/home/songyu/fortran_test/POD_Fortran/.worktrees/boxwing-modes-p1/docs/superpowers/reports/2026-10-14-p4-reuse-audit.md)。

## 6. 回归与尚未满足的整体验收项

全套原生 Fortran 测试在私有副本 `/home/songyu/fortran_test/POD_Fortran/.worktrees/p4-native-20261015.NwRuiu` 执行，OPM/config/SRP 输出均与历史目录隔离，每项限时 300 s。101 个目标实际全部尝试：**91 通过，7 超时，3 缺输入文件失败**。不得将退出码 124 或缺文件称为数值通过。

7 个超时目标：`test_da_orbit_propagation`、`test_hfem_ads_propagation`、`test_halo_srp_initial_state_uncertainty_15day`、`test_uq`、`test_hfem_uq_timing_table`、`test_halo_srp_scale_uncertainty_15day`、`test_da_integrator`。其覆盖的长时传播/性能实验未在该时限内完成；本次不降低实验规模或断言门槛以制造通过，也不能从超时断言没有回归。

3 个缺文件失败：

- `test_ut_gap_propagation` 缺 `output/L1Halo-2_emdac_n5_before_gap.txt`。
- `test_stt_crtbp` 缺 `../rand_list300km_1.0ms.txt`。
- `test_uq_crtbp_comparison` 缺 `../rand_list200km_0.7ms.txt`。

8D 模式表/参数映射/代理力、嵌入维度、实际 8D HFEM、事件/研究点/v2 快照、旧 6D/7D 映射和相关 SRP/重建原生目标通过。全部结果逐项保存在 [原生回归索引](/home/songyu/fortran_test/POD_Fortran/.worktrees/p4-native-20261015.NwRuiu/regression_logs/results.jsonl) 与对应日志。

Python 全发现初次运行暴露 6 个模块依赖不可加载（numpy/pandas/h5py 缺失），以及测试环境误把 `POD_SRP_ADS_EXE` 指向另一种 CLI。已核对并改正该环境变量为 `run_srp_ads_history`，未修改旧接口。当前可加载的全部模块重新执行 **67 tests，全部通过**，包含新增分析 6 项、Real truth 4 项和原生/正式 C2 runner 回归。6 个依赖不满足的模块未声称运行通过，没有安装新依赖或改写无关 MC 代码。

日志：[Python 可运行模块复验](/home/songyu/fortran_test/POD_Fortran/.worktrees/p4-native-20261015.NwRuiu/python-available-modules.log)、[全发现诊断](/home/songyu/fortran_test/POD_Fortran/.worktrees/p4-native-20261015.NwRuiu/python-test.log)、[完整 FPM 日志](/home/songyu/fortran_test/POD_Fortran/.worktrees/p4-native-20261015.NwRuiu/fpm-test.log)。

本轮新增容差接口、Real truth 与分析逻辑经过 TDD 红/绿测试及主 agent 逐项复核；隔离回归 launcher 已实际执行全部 101 个目标。独立子 agent 审查因额度/容量限制未能完成，不作独立复审已通过的声明。上述回归与审查差项意味着设计 §17.9 及完整工程验收尚未全部关闭；**这不改变已实际完成并核验的完整 C2 轨道案例及其三组收敛结果。**

## 7. 可复现命令与身份

以下命令从 `/home/songyu/fortran_test/POD_Fortran/.worktrees/boxwing-modes-p1` 执行。既有正式输出是只读输入；分析输出必须使用不存在的新目录，不能覆盖本文结果。

```bash
source setup_env.sh
python3 tools/analyze_srp_modes_case.py \
  --baseline ../../SRP/boxwing_modes_DRO-2_init_c2_15d_20261014 \
  --half-hour ../../SRP/boxwing_modes_DRO-2_init_c2_half_hour_20261015 \
  --tight-rkf ../../SRP/boxwing_modes_DRO-2_init_c2_tight_rkf_20261015 \
  --da5 ../../SRP/boxwing_modes_DRO-2_init_c2_da5_20261015 \
  --truth ../../SRP/boxwing_modes_DRO-2_init_truth_20261015 \
  --c0 ../../SRP/SRP_UCERTAINTY_PROP_WITHINIT_ETA0/DRO-2 \
  --c1 ../../SRP/SRP_UCERTAINTY_PROP_WITHINIT/DRO-2 \
  --output ../../SRP/NEW_EXCLUSIVE_ANALYSIS_DIRECTORY
```

三组变体均调用 `tools/run_srp_modes_ads.py --hours 360 --save-hours 1 --max-depth 8 --pos-tol-km .1 --vel-tol-kms 1e-6`，只按 §1 改表路径、DA 阶或 rel/abs，使用同一冻结 runner。新 Real truth 调用 `tools/run_srp_modes_truth.py`，同 OPM/config，360 h，每小时保存，rel/abs=1e-13。原始 `run.log`、正式报告与冻结 runner 均保留。

```text
OPM: 982f389d31af0c6a986bdf0ea87ecc7e4d4b09c468bfc402ab3bfef2781d49d8
Hourly mode CSV: f82cd81ee5aab0e833c2ec4d9015d7df80ff63d6a78d9c8b808c2207e9f69f5e
Hourly metadata: ee42d4d41084788fdb4e4ec5723e890f55c961507c08ca0804c63f35d9da2d12
Original C2 runner: 89c03449d2d5d72da311cf2a989a721488cf5181299ba50b7ec6770c89a0ef62
Variant runner: 7d8cfd752136cff5d352edb48be7f0ee60d81f6b3bba064ac60129e385ce35a2
Real truth runner: 6bff2afa76142b19ecd50cb3e5f77f7ec1a4b9f1aeb12e2cefbc270805cef9a2
Analysis source: 9f468bbe7d1b955341a3ea4f5da5df0aec6da91918123e3712f50255fe83bc67
Case report: 8509cc64a31203a44cead17e83a2d84af7174bca83a2ac030ed170edd314e37c
Hourly metrics: 3011bceed4c93d51ec4350d3277ad35e19dd36cebd915da6b388986494b607fc
```

基线 HEAD `8d999bfc317eee44e22f36e56e0a2e2e87f493e1` 不包含未提交实现；正式报告逐文件 `source_files_sha256` 才是本轮源码身份，不把基线 revision 冒充实现提交。旧结果只读复用，未改核心积分器、ADS 分裂核心或既有 runner 的 CLI/schema。
