# Box-Wing 时变 SRP 误差模式 DA/ADS 传播设计

原始日期：2026-10-08

修订日期：2026-10-16

状态：P1/P2/P3 完成；P4 完整 C2 360 h 案例、独立 Real 比较及三组收敛通过；全仓回归与整体工程验收仍有差项。概率输出扩展已实现并通过 50k/2k 两小时 smoke，完整概率运行须单独验收。用户已授权先本地提交、再启动 360 h 概率传播；实际状态以新运行日志和报告为准。

适用范围：POD_Fortran 的 Real 与 DA 高保真动力学传播；算例 `DRO-2_init`，15 天。

## 1. 决策摘要

在现有恒定标量 SRP 误差 `eta` 模型之外，新增时变双模式 SRP 代理模型。从同一条 Cannonball 名义参考轨道上提取 Box-wing 对地、对月相对 Cannonball 的三维时变误差模式 `e_E(t)`、`e_M(t)`，再用两个固定但不确定的系数 `(s,u)` 构成 8 维 DA/ADS：

~~~text
a_SRP(t,x; xi_E,xi_M)
  = a_CB(t,x) + ||a_CB(t,x)|| [xi_E e_E(t) + xi_M e_M(t)]

xi_E = s
xi_M = (1-s)u
s,u in [0,1]
~~~

方形参数域经上述映射覆盖三角形系数域：

~~~text
xi_E >= 0
xi_M >= 0
xi_E + xi_M <= 1
~~~

工程路线确定为：

1. 同状态模式提取是独立 go/no-go 数据门禁；未通过 30% 矢量上界时，不进入 8D 实现。
2. ADS 坐标映射增加显式参数化类型，不再用变量数猜测物理含义。
3. 新增独立 C2 runner；现有 6D/7D runner 的 CLI 与输出语义保持不变。
4. 只泛化 C2 必需的共享基础设施：HFEM 事件、研究点、快照和离线重建。
5. Real 与 DA 使用同一模式表和同一代理公式，并分别提供显式 `(s,u)` 接口。

本设计不对 `(s,u)` 赋予概率测度，只做有界域包络、确定性端点和代理误差验证。

上述限制继续适用于本设计既有包络实验和已发布结果。用户于 2026-10-15 选择新增 8D 联合概率输出：初态白化 ±3σ 截断高斯、s/u 独立均匀先验；它作为显式启用的新实验，按 [概率输出扩展规格](/home/songyu/fortran_test/POD_Fortran/docs/superpowers/specs/2026-10-15-boxwing-modes-probability-output-design.md) 单独复核与实现，不追溯改变旧结果的解释。

## 2. 目标与非目标

### 2.1 目标

- 在同一条确定性 Cannonball 参考轨道上重算 Cannonball、Box-wing 对地、Box-wing 对月瞬时 SRP 加速度。
- 生成带完整物理身份和来源信息的可复现模式表。
- 建立 6 个初始状态变量加 2 个模式变量的 8D DA/ADS 域。
- 保证 ADS 分裂后 `(s,u)` 仍表示正确的全局物理系数。
- 完成 C2 代理模型的 Real/DA 一致性验证和 15 天 ADS 传播。
- 分开报告 DA/ADS 数值误差、时变代理误差和 SRP 参数导致的条件轨道偏移。
- 保持现有 6D/7D C0/C1 行为、CLI、文件格式和测试不回归。

### 2.2 非目标

- 不引入三角形均匀先验，不计算模式参数概率矩。
- 不重写 RKF 积分器、`pod_ads_split_module` 或 DACE 包装层核心算法。
- 不修改 Cannonball 和 Box-wing 原始物理公式。
- 不引入阴影或半影；MVP 沿用当前 full-sun 口径。
- 不把结论推广到其他轨道、历元或航天器构型。
- 不把现有 7D `eta` runner 原地改造成同时承担 8D 实验的多用途程序。

## 3. 实施前代码约束与设计响应

以下是设计修订时的实施前基线；实施后的当前状态见 §16.2：

- `pod_uq_ads_coordinates_module` 只能通过 `srp_sigma` 在 6D/7D 之间选择，第 7 维固定解释为 `eta`。
- ADS 全局根域和补丁局部域均使用 `[-1,1]` 坐标。
- HFEM 对 `eta` 使用补丁中心和宽度重新构造局部 DA 参数；模式参数必须采用同样机制。
- HFEM 事件、研究点、快照、离线重建和两个现有 runner 含固定 7 维结构。
- 当前 Real 验证探针只接受 `eta`。
- `pod_engine_init` 尚不加载 SRP 模式表。

因此，只把数组长度从 7 改为 8 不是完整实现。本设计明确增加参数化枚举、8D 坐标映射、patch-local 参数注入、模式表生命周期和独立 C2 runner。

## 4. 阶段与门禁

| 阶段 | 内容 | 进入下一阶段的门禁 |
|---|---|---|
| P1 | 同状态模式提取与数据审计 | `max ||e_E|| <= 0.3` 且 `max ||e_M|| <= 0.3` |
| P2 | 模式表、Real/DA 力模型与 8D 坐标基础设施 | 瞬时端点、坐标往返、分裂后参数映射通过 |
| P3 | 8D runner、快照和独立 Real 探针 | 1 小时 smoke 满足 ADS 数值阈值 |
| P4 | 15 天实验、代理误差与收敛性报告 | 完整输出、无非法求值、回归通过 |

P1 是正式 go/no-go 门禁。任一模式超过 0.3 时，保留提取结果并报告实测上界，但停止使用“30% 矢量上界”和后续固定 30% 对照口径；由用户决定扩大参数口径或修改代理模型后再继续。

## 5. P1：同状态模式提取

### 5.1 独立提取程序

新增：

~~~text
app/run_boxwing_mode_extract.f90
~~~

长时数据生成不放入默认 `fpm test`。程序执行：

1. `pod_engine_init(config/config_box_wing_earth.txt)`，读取完整 Box-wing 几何与共同动力学配置。当前 `config/config.txt` 的几何字段为 `-1` 占位值，不能直接用于提取；程序随后强制 Cannonball 参考传播，输入配置的模型选择不改变参考轨道。
2. 读取 `OPM/DRO-2/DRO-2_init.opm.json` 的 `epoch0,state0`。
3. 强制 Cannonball、`eta=0`，用 RKF78 传播 15 天。
4. 使用 `REL_TOL=ABS_TOL=1e-12`、`DT_MIN=1e-6 s`、`DT_MAX=3600 s`。默认 `--hours 360 --sample-hours 1`，保存 0..360 h 共 361 个整点状态；允许 `--sample-hours 0.5` 生成 721 个节点。时长必须包含整数个采样区间。
5. 在每个 `t_k,x_ref(t_k)` 上，用同一位置、速度和绝对 TDB 历元重算 `a_CB,a_BWE,a_BWM`。
6. 计算：

~~~text
e_E(t_k) = [a_BWE(t_k) - a_CB(t_k)] / ||a_CB(t_k)||
e_M(t_k) = [a_BWM(t_k) - a_CB(t_k)] / ||a_CB(t_k)||
~~~

所有状态、加速度和模式分量必须有限，`||a_CB||>0`，整点时间误差不超过 `1e-6 s`。

### 5.2 输出

~~~text
SRP/boxwing_modes_DRO-2_init/
  boxwing_modes.csv
  boxwing_modes_meta.json
~~~

CSV 固定列为：

~~~text
time_hours,epoch_tdb_s,
ref_x_km,ref_y_km,ref_z_km,ref_vx_km_s,ref_vy_km_s,ref_vz_km_s,
a_cb_x_km_s2,a_cb_y_km_s2,a_cb_z_km_s2,a_cb_norm_km_s2,
a_bwe_x_km_s2,a_bwe_y_km_s2,a_bwe_z_km_s2,
a_bwm_x_km_s2,a_bwm_y_km_s2,a_bwm_z_km_s2,
eE_x,eE_y,eE_z,eE_norm,
eM_x,eM_y,eM_z,eM_norm
~~~

metadata 至少记录：

- schema 版本；
- OPM 路径、SHA-256、`epoch0`；
- 配置路径、SHA-256、源码 revision；
- J2000、EARTH、TDB 和单位；
- full-sun 标识、实际采样间隔（默认 1 h）；
- Cannonball/Box-wing 物理参数；
- 行数、起止历元、两个最大模式范数及发生时刻；
- 数据 CSV 的 SHA-256。

提取程序先写 staging CSV 和基础 metadata；`tools/audit_boxwing_mode_table.py --finalize` 重算数据、补齐 SHA-256，并在全部审计通过后原子发布上述两个正式文件。正式 C2 运行只接受已完成审计的 metadata。

正式目录额外保存冻结的提取器 executable，metadata 使用相对文件名 `extractor` 及其 SHA-256。审计通过 `--utc-to-tdb` 接口调用项目 SPICE，独立比较输入 OPM EPOCH 与表的 TDB 起点，容差 `1e-6 s`；冻结副本不受后续编译影响。

### 5.3 P1 验收

- 默认恰有 361 行，`time_hours=0..360` 严格递增且无重复；0.5 h 加密采样恰有 721 行。行数、时间网格和 metadata 必须一致。
- `epoch_tdb_s=epoch0+3600*time_hours`，误差不超过 `1e-6 s`。
- 差分、分量范数和 metadata 最大值可从 CSV 重算。
- J2000 与同一参考状态构造的 RTN 旋转保持向量范数，相对误差不超过 `1e-12`。
- `max ||e_E||<=0.3` 且 `max ||e_M||<=0.3`。

分段线性插值的向量范数是区间参数的凸函数，因此每段最大值不超过两个端点的较大值；节点门禁同时保证插值段内不突破 0.3。

### 5.4 P1 实测结果（2026-10-13）

使用仓库 DRO-2 初值、完整 Earth 配置中的公共物理参数和 Cannonball 参考轨道，完成 15 天同状态提取：

| 采样间隔 | 节点数 | 最大 Earth 模式范数 | 时刻 / h | 最大 Moon 模式范数 | 时刻 / h |
|---|---:|---:|---:|---:|---:|
| 1 h | 361 | 0.2728419545768155 | 325 | 0.2733193844007026 | 153 |
| 0.5 h | 721 | 0.2728419545768162 | 325 | 0.2733193844007047 | 153 |

两组均通过 30% 门禁，metadata 已 finalized。1 h 独立重复运行的 CSV 字节完全一致。全部数据通过初值身份、时间、有限性、模式差分、范数、RTN 旋转和哈希审计，几何 fallback 计数均为 0。加密采样确认本算例的模式最大值稳定，但不替代 P4 的代理传播、DA 阶数和积分容差收敛验证。

正式数据位于 `SRP/boxwing_modes_DRO-2_init/`，0.5 h 对照位于 `SRP/boxwing_modes_DRO-2_init_half_hour/`。P1 代码与测试位于分支 `codex/boxwing-modes-p1` 的 `.worktrees/boxwing-modes-p1/` 工作区，尚未合并到主分支。

## 6. 模式表与初始化生命周期

新增 `src/lib/forcemodel/pod_srp_mode_module.f90`，公开接口：

~~~fortran
load_srp_mode_table(csv_file, meta_file, expected_epoch0, status, message)
clear_srp_mode_table()
srp_mode_table_is_loaded()
interpolate_srp_modes(epoch_tdb_s, e_E, e_M, status, message)
~~~

模块只持有普通 `real(DP)` 节点和向量，不持有 DA 句柄。加载时验证：

- CSV header、metadata schema、行数、单调时间、有限值和存储范数；
- metadata 已标记为审计完成；
- metadata 历元与 `expected_epoch0` 在 `1e-6 s` 内一致；
- frame、center、timescale、单位和 full-sun 标识；
- metadata 的 Cannonball 公共物理参数（Cr、面积、质量、有效 1 AU 光压）与当前配置一致，相对容差 `1e-12`；Box-wing 几何是 P1 冻结来源身份，由 Python 审计其源配置，不要求模式运行配置重复提供有效几何；
- 覆盖区间至少为 `[epoch0,epoch0+15 days]`。

Fortran 运行时不自行实现 SHA-256；`tools/audit_boxwing_mode_table.py` 重算并验证 OPM、配置和 CSV 哈希，再原子写入带 `audit_complete=true` 的最终 metadata。Fortran loader 负责检查 finalized 标识、语义、历元和物理参数，两者共同构成数据身份门禁。

P3 原生 C2 runner 为可重建 v2 快照调用 Linux GNU `sha256sum` 计算 metadata 实际摘要；可选 `--mode-meta-sha256` 仅用于核对预期身份。此摘要不等于完整数据审计，原生输出仍标记 `formal_result=false`；正式发布必须经过 Python wrapper。

正式 runner 每次运行前仍须执行 Python 哈希审计，不能仅凭 loader 接受 `audit_complete=true` 推断数据未被改动。公共参数在加载后被修改时，插值也必须失败。`srp_pressure_1au_n_m2=-1` 按 `1367/c` 比较；P1 实测表使用显式值 `4.559821181419152e-6`，运行配置必须保持该值。

`pod_engine_init` 每次读取配置时先清理旧模式表，但不加载新表，因为此时尚不知道 OPM 历元。C2 runner 的固定顺序为：

1. `pod_engine_init(config_file)`；
2. `load_initial_opm` 得到 `epoch0`；
3. `load_srp_mode_table(...,expected_epoch0=epoch0)` 验证 finalized metadata；
4. 初始化 ADS。

任何 `time_varying_modes` 力计算在模式表未加载时都必须失败。其他模型保持模式表为空，避免同进程切换配置时复用旧数据。

插值输入统一为绝对 TDB 秒。首末节点允许 `1e-6 s` 舍入容差并夹到端点；真正越界时报错，不外推。

## 7. 配置契约

新增：

~~~text
srp_model = cannonball | box_wing | time_varying_modes
srp_mode_table_file = SRP/boxwing_modes_DRO-2_init/boxwing_modes.csv
srp_mode_metadata_file = SRP/boxwing_modes_DRO-2_init/boxwing_modes_meta.json
~~~

只支持分段线性插值，不增加可配置插值算法字段。

- `time_varying_modes` 必须提供两个存在的文件。
- 它使用 Cannonball 公共物理参数，但不要求 Box-wing 几何字段有效。
- `box_wing` 继续执行现有几何验证。
- 模式系数不存入配置；Real 探针和 ADS 映射通过显式接口传入。

## 8. Real/DA 力模型接口

### 8.1 Cannonball helper

Real 与 DA 力模型内部各提取一个不读取 `srp_model` 的 Cannonball full-sun helper。Cannonball 分支和模式分支均调用它，禁止模式分支重新调用顶层 SRP 分派函数。

### 8.2 Real 系数

新增：

~~~fortran
set_srp_mode_coefficients_real(s, u)
clear_srp_mode_coefficients_real()
~~~

`s,u` 必须有限且位于 `[0,1]`。`clear_srp_mode_coefficients_real` 将状态标记为“未设置”；`time_varying_modes` Real 分支在未显式设置时必须失败。即使零模式也必须显式调用 `set_srp_mode_coefficients_real(0,0)`，结束和错误出口再清理。

### 8.3 DA patch-local 系数

新增：

~~~fortran
set_srp_mode_uncertainty(s_index, s_nominal, s_span, &
                         u_index, u_nominal, u_span)
clear_srp_mode_uncertainty()
~~~

ADS 全局坐标与物理系数的关系为：

~~~text
s = 0.5(1+q7)
u = 0.5(1+q8)
q7,q8 in [-1,1]
~~~

若补丁历史给出全局坐标中心 `c` 和宽度 `w`，则局部 DA 参数必须设置为：

~~~text
s_nominal = 0.5(1+c7)
s_span    = 0.25 w7
u_nominal = 0.5(1+c8)
u_span    = 0.25 w8

s_da = s_nominal + s_span da_var(7)
u_da = u_nominal + u_span da_var(8)
~~~

不得直接把补丁局部 `da_var(7/8)` 当作物理 `s/u`。

### 8.4 代理力

~~~fortran
a_cb = cannonball_fullsun(t, x)
call interpolate_srp_modes(t, e_E, e_M, status, message)
xi_E = s
xi_M = (1-s)*u
a_srp = a_cb + norm(a_cb)*(xi_E*e_E + xi_M*e_M)
~~~

Real 与 DA 使用同一模式表、绝对历元和公式。时间函数为普通浮点量，状态与 `(s,u)` 为 DA 量。

## 9. ADS 参数化与坐标映射

新增公开枚举：

~~~fortran
ADS_PARAM_AUTO      ! 仅供旧调用兼容
ADS_PARAM_NONE
ADS_PARAM_SRP_ETA
ADS_PARAM_SRP_MODES
~~~

`hfem_ads_options_type` 与 `ads_coordinate_map_type` 增加 `parameterization`：

| 参数化 | 活跃维数 | 参数维 |
|---|---:|---|
| `NONE` | 6 | 无 |
| `SRP_ETA` | 7 | `eta` |
| `SRP_MODES` | 8 | `s,u` |

`hfem_ads_options_type` 默认使用 `ADS_PARAM_AUTO`：`srp_sigma>0` 解析为 `SRP_ETA`，否则解析为 `NONE`；解析后的 `ads_coordinate_map_type` 只保存 `NONE/ETA/MODES`，不保存 `AUTO`。C2 runner 必须显式选择 `SRP_MODES`。

新增统一物理参数对象：

~~~fortran
type ads_physical_parameters_type
    real(DP) :: eta = 0
    real(DP) :: s = 0
    real(DP) :: u = 0
end type
~~~

坐标转换根据参数化处理：

~~~text
s = 0.5(1+q7), q7 = 2s-1
u = 0.5(1+q8), q8 = 2u-1
~~~

保留现有 6D/7D wrapper，避免同时破坏旧调用方。物理参数超出 `[0,1]` 超过 `1000*epsilon` 时失败；只在该容差内夹到端点。

## 10. HFEM 与配套链路的 8D 扩展

### 10.1 HFEM

以下结构全部支持 8 维：

- `direction_scores`；
- `parent_center,parent_width`；
- `first_split_direction_scores`；
- `first_split_parent_center,first_split_parent_width`；
- 两条传播路径的局部评分数组；
- `record_split_event` 有效范围；
- `da_storage_variables` 上限。

每个补丁传播前按参数化三分派：

- `NONE`：清理 `eta` 和模式 DA 状态；
- `SRP_ETA`：保持现有第 7 维 patch-local 行为；
- `SRP_MODES`：按 §8.3 设置第 7/8 维。

输出样本必须保留全部物理参数列；8D 输出不能只复制第 7 维。

### 10.2 研究点与 Real 探针

`build_study_points` 增加 `n_variables` 参数并支持 6/7/8 维；Halton 素数表扩为 8 个。旧调用显式传 7，保持原点集不变。

新增：

~~~fortran
advance_real_mode_probes(states, s_values, u_values, epoch0, &
                         t_start, t_end, options, status, message)
~~~

它逐点设置 `(s,u)`、传播同一代理模型，并在所有出口清理参数状态。现有 `advance_real_probes(...,eta,...)` 保持不变。

### 10.3 快照与离线重建

快照按 `history%da_variables` 动态生成：

- `center1..centerN,width1..widthN`；
- `power1..powerN`；
- 参数化类型和变量角色；
- 模式 CSV/metadata SHA-256。

新 schema 为 v2。`tools/reconstruct_srp_ads_snapshot.py` 同时支持旧 6D/7D v1 和新 6D/7D/8D v2；既有快照不得失效。增加 8D 多项式、分裂历史和坐标往返测试。

## 11. C2 独立 runner

新增：

~~~text
app/run_srp_modes_ads.f90
~~~

必需参数：

~~~text
-opm INPUT.opm.json
-o OUTPUT_PREFIX
-cfg CONFIG.txt
~~~

配置必须为 `srp_model=time_varying_modes`。runner 固定建立 `ADS_PARAM_SRP_MODES` 的 8D 域，并输出：

- 每小时名义状态、域包络、补丁数和验证误差；
- 8 个方向的分裂统计和完整事件；
- 8D 研究点与独立 Real 探针结果；
- 起点、终点和事件时刻的 v2 快照；
- 原始 JSON 报告记录所有输入路径、metadata 身份和源码 revision；
- 正式运行 wrapper 重算并记录 OPM/config/mode CSV/metadata/executable 的 SHA-256。

默认值沿用现有实验：15 天、1 h 保存、DA 4 阶、最大深度 8、位置阈值 0.1 km、速度阈值 `1e-6 km/s`、RKF78 `1e-12/1e-12`。

正式实验通过 `tools/run_srp_modes_ads.py` 启动。wrapper 在运行前复核 finalized mode metadata 和全部输入哈希，将 Fortran 输出写入 staging，补齐报告哈希并通过结构审计后原子发布；失败结果留在 staging。Fortran 可执行程序仍可用于 focused smoke，但其原始输出不得直接标记为正式结果。

现有 `run_srp_ads_history` 和 `run_srp_ads_mc_pair` 的 CLI、默认行为和 schema 不变。

## 12. 实验组与标量口径

| 组 | DA 维 | 物理域 | 用途 |
|---|---:|---|---|
| C0 | 6 | Cannonball，无额外参数 | 控制组 |
| C1-30 | 7 | `eta in [-0.3,0.3]` | 与 30% 口径匹配的标量域 |
| C1-H | 7 | 现有 `eta_sigma=0.3,domain_sigma=3`，即 `[-0.9,0.9]` | 历史背景，不解释为等宽域 |
| C2 | 8 | `s,u in [0,1]` | 时变双模式处理组 |
| D0 | Real | CB，`eta=0` | 确定性端点 |
| Deta | Real | CB，`eta=+0.3` | 确定性 30% 标量端点 |
| D1 | Real | 完整 Box-wing Earth | 代理真值端点 |
| D2 | Real | 完整 Box-wing Moon | 代理真值端点 |

C1-30 在现有映射下使用 `domain_sigma=3,srp_sigma=0.1`。C1-H 可复用已有结果，但所有报告必须标注其物理域为 ±0.9。

C0、C1-30、C2 使用相同 OPM、协方差、历元、非 SRP 动力学、积分器、DA 阶数、ADS 容差和输出时刻，并记录所有输入哈希。

## 13. 验证与数值验收

### 13.1 单元与接口

- **T1 表解析**：拒绝错误 header、时间、有限性、哈希和历元。
- **T2 插值**：节点精确返回，中点为相邻节点平均，端点容差和越界失败正确。
- **T3 参数域**：坐标往返覆盖中心、角点、边界和随机点，三角形约束成立。
- **T4 补丁映射**：根、左、右补丁中的全局 `(s,u)` 物理值一致。
- **T5 三端点**：361 个参考状态上，`(0,0)`、`(1,任意u)`、`(0,1)` 分别复现 CB、Earth、Moon；加速度相对误差 `<=1e-12`。
- **T6 Real/DA 常数项**：代表性历元和参数点相对误差 `<=1e-12`。
- **T7 生命周期**：8D 构造、分裂、快照、销毁后无 DACE 句柄增长。

### 13.2 传播

1. C2 `(0,0)` 与 C0 每小时状态差：
   - 位置范数 `<=1e-6 km`；
   - 速度范数 `<=1e-9 km/s`。
2. 指定 8D 点的 ADS 求值与同一代理模型 Real 积分：
   - 位置 `<=0.1 km`；
   - 速度 `<=1e-6 km/s`。
3. C2 Earth/Moon 顶点代理传播分别与 D1/D2 比较，报告每小时和全程最大误差。代理误差是模型适用性指标，不作为软件正确性门禁。
4. C2 参数变化相对 C0 的状态差单独报告，不与 ADS 数值误差混名。
5. 至少执行 0.5 h 模式采样、更严积分容差和更高 DA 阶三组收敛性检查。关键指标相对基准变化 `<=5%` 记为“已确认收敛”，否则报告“未确认”。

输出名称统一使用 `ads_numerical_error`、`surrogate_model_error`、`conditional_srp_offset`。

## 14. 错误处理

- 表未加载、身份不匹配、时间越界或插值失败时立即报错，不回退 Cannonball。
- setter 拒绝非有限值、非法索引和越界系数。
- runner 在初始化、结束和所有错误出口清理 `eta` 与模式参数状态。
- 传播失败时不发布正式结果；临时输出保留诊断。
- 达到最大分裂深度可保存结果，但标记 `precision_unconfirmed=true`。
- 代理误差不得标成 ADS 数值误差，反之亦然。

## 15. 文件范围

### 15.1 新增

- `src/lib/forcemodel/pod_srp_mode_module.f90`
- `app/run_boxwing_mode_extract.f90`
- `app/run_srp_modes_ads.f90`
- `tools/audit_boxwing_mode_table.py`
- `tools/run_srp_modes_ads.py`
- `test/test_audit_boxwing_mode_table.py`
- `test/test_boxwing_mode_extract.py`
- `test/test_run_srp_modes_ads.py`
- `test/test_srp_mode_table.f90`
- `test/test_srp_mode_parameter_map.f90`
- `test/test_srp_modes_force_consistency.f90`
- `test/test_hfem_ads_modes_8d.f90`
- `test/test_srp_modes_ads_runner.py`

### 15.2 修改

- `src/lib/api/pod_engine_module.f90`：模式表生命周期。
- `src/lib/system/pod_config_module.f90`：新模型、路径和验证。
- Real/DA force modules：Cannonball helper、setter 和模式分派。
- `pod_uq_ads_coordinates_module.f90`：参数化枚举、8D 映射、兼容 wrapper。
- `pod_uq_hfem_ads_module.f90`：8D 事件、参数注入和输出。
- `pod_srp_ads_history_module.f90`：可指定维数的研究点和模式 Real 探针。
- `pod_srp_ads_snapshot_module.f90`：动态维度 v2 快照。
- `tools/reconstruct_srp_ads_snapshot.py`：兼容 v1 的 8D v2 重建。
- 对应现有 Fortran/Python 测试。

`fpm.toml` 已启用自动发现；仅在目标名冲突时增加显式条目。

### 15.3 保持不动

- RKF 积分器；
- `pod_ads_split_module` 核心；
- `pod_dace_classes`；
- 现有两个 runner 的 CLI/schema；
- 已有 C0/C1 结果文件。

## 16. 实施顺序

1. P1 extractor、审计和 30% go/no-go。
2. 模式表解析、metadata 校验和初始化生命周期。
3. 参数化枚举、8D 坐标往返和 patch-local 映射。
4. Real/DA 代理力与瞬时端点。
5. HFEM 事件、研究点、快照和离线重建。
6. C2 runner 与 1 小时 smoke。
7. C0/C1-30/C2 15 天传播和 D0/Deta/D1/D2 确定性传播。
8. 收敛性、完整回归和最终分类报告。

任何阶段失败时停在当前阶段修复，不绕过门禁启动正式运行。

### 16.1 P2 实施状态（2026-10-13）

模式表、typed metadata、配置/engine 生命周期、Real/DA 代理 SRP、8D 坐标及 patch-local 系数已实现于现有 `codex/boxwing-modes-p1` 隔离工作区。361 个 P1 节点上五种参数组合的 Real 最大相对误差为 0；代表节点及离参考状态的 DA/Real 常数项最大相对误差为 `6.8643e-16`，小于 `1e-12`。重复 500 次代理力计算不增加 DA 句柄，清理后回到 0。实际 ADS 第 7、8 维左右分裂保持全局 s/u。

17 个 Fortran focused targets、Python P1 审计/提取器/旧快照测试及 16 条独立进程错误路径通过；两份 P1 正式表重审通过。尚未运行完整 `fpm test`，未声称 1 h/15 天 C2 精度通过。

在 P2 结束时，`hfem_ads_options_type` 的 MODES 入口、HFEM 8D 事件/评分/参数注入、研究点、v2 快照与 C2 runner 仍未实施；该历史状态已由下述 P3 覆盖。详见 [P2 结果报告](/home/songyu/fortran_test/POD_Fortran/docs/superpowers/reports/2026-10-13-boxwing-modes-p2-results.md)。

### 16.2 P3 实施状态（2026-10-14）

完整 8D HFEM/ADS、512 个研究点、337 个独立 Real 探针、8D 事件与 v2 快照/离线重建、独立 C2 runner 和审计发布 wrapper 已实现。旧 6D/7D runner CLI/schema 与 v1 读取保持兼容；RKF 和 ADS 分裂核心不变。

最终 1 h smoke 位于 `SRP/boxwing_modes_DRO-2_init_p3_smoke_20261014_verified`，DA 阶 4、深度 8、RKF78 容差 `1e-12`。ADS 对同代理 Real 最大位置/速度误差为 `3.519638739074529e-10 km` / `1.6541950504657498e-15 km/s`；零模式 ADS 对同动力学 Cannonball Real 名义控制为 `6.66852530405717e-11 km` / `2.482534153247273e-16 km/s`，均通过 §13 门槛。句柄数 0→0，未达到深度上限；初始/终点 v2 重构对保存点云误差为 0。

21 个 Fortran focused targets、40 个 Python tests 通过，两份 P1 表重审通过。独立审查发现的 DA 完成检查转置、默认原生 metadata 哈希和 6/7D 嵌入 8D 重构列数问题已由红/绿回归修复，复审无剩余实质问题。源码仍未提交或合并，尚未运行完整 `fpm test`。

P3 的零模式比较使用 Cannonball Real 名义控制，不代表已运行完整 C0/C1-30 对照域。本次 smoke 未自然分裂，分裂路径由单元测试覆盖。P4 的 15 天传播、真实 Box-Wing 代理误差和三组收敛检查尚未执行，整体 §17 完成标准尚未达到。详见 [P3 验收报告](/home/songyu/fortran_test/POD_Fortran/docs/superpowers/reports/2026-10-14-boxwing-modes-p3-results.md)。

### 16.3 P4 C2 启动历史（2026-10-14，已由 16.4 更新）

按用户指示，优先复用已有 C0/C1 结果，不重新运行对照组，直接启动 C2 的 360 h 基准。P1 数据重新审计通过，源码与 P3 验收身份一致；使用 P3 冻结 executable，DA 阶 4、深度 8、每小时保存。独立后台 wrapper/native 已确认跨启动 shell 存活并完成实际 checkpoint；本次最终目录为 `SRP/boxwing_modes_DRO-2_init_c2_15d_20261014`，仅在运行与审计成功后发布。

已有 6D C0 的 OPM 哈希、配置来源及非 SRP 设置匹配，可优先复用。当前找到的标量历史组为 eta±0.9（C1-H），不得直接标成 C1-30；未因此启动替代对照任务。C2 内部 Real 探针与零模式名义控制仍按验收要求执行。15 天最终数值门禁、真值代理误差、收敛检查与整个 P4 尚未通过。进程、日志、输入身份和初步复用核对见 [P4 C2 运行记录](/home/songyu/fortran_test/POD_Fortran/docs/superpowers/reports/2026-10-14-boxwing-modes-p4-c2-run.md)。

### 16.4 P4 完整轨道案例核验（2026-10-15）

原 C2 基线及模式表 0.5 h、RKF rel/abs=1e-13、DA 阶 5 三组变体均已实际完成 360 h、361 个时刻并正式发布。四组相同 512 个研究点、337 个独立 Real 探针；ADS 与零模式门禁均通过，句柄 0→0，没有深度受限或精度未确认。基线最大 ADS 误差 `1.863267429652309e-5 km / 2.3448859556127377e-10 km/s`。逐小时共同 C0 参考下的采样包络及 Earth/Moon 末时刻物理位移，三组变化均低于 5% 收敛门槛。四组初/末 v2 快照已只读重新求值验证。

新增五条完整独立 Real 名义轨迹覆盖 Deta、proxy Earth/Moon、完整 Box-Wing Earth/Moon；总计 1805 条状态，rel/abs=1e-13，几何 fallback 为 0。Earth/Moon proxy 对完整模型的全弧最大位置差约 0.30 m / 0.70 m；完整物理模型对炮弹球 C0 参考的末位置位移约 17.063 km / 37.952 km。代理误差、ADS 数值误差与条件物理位移保持分类。

按用户要求不重跑原 C0/C1 域或 C2 基线。C0 初态/基/非 SRP 配置可复用，但默认与显式光压有 `1.3249e-11` 相对差，已记录而非称完全一致；实际名义轨迹差约 `9.81e-9 km`。历史 C1 是 eta±.9 的 C1-H，仅对其末时刻多项式做 eta±.3 子域重建；旧报告缺失来源哈希的限制没有被掩盖，也没有声称获得完整匹配的逐小时 C1-30。

完整原生回归 101 项已尝试：91 通过、7 达到单项 300 s 时限、3 缺历史输入文件失败；可加载的全部 Python 模块 67 tests 通过，另有 6 个模块因 numpy/pandas/h5py 缺失未能加载。因此**完整 15 天案例与收敛已验证，但 §17.9 全仓回归通过及整体完成标准尚未全部满足**。独立子 agent 复审受额度/容量限制，本轮只声明主 agent 复核。源码仍未提交或合并。详见 [P4 完整案例核验报告](/home/songyu/fortran_test/POD_Fortran/docs/superpowers/reports/2026-10-15-boxwing-modes-p4-complete-case-results.md)。

## 17. 完成标准

1. P1 模式表可复现、身份完整并通过 30% 门禁。
2. 8D 根域和分裂补丁保持 `s,u in [0,1]` 及正确全局映射。
3. Real/DA 端点和常数项达到 §13 阈值。
4. 1 小时及 15 天 C2 无非法 DA 求值、NaN、越界插值或句柄泄漏。
5. v2 8D 快照可重建，既有 v1 仍可读取。
6. 零模式、ADS 数值误差和收敛状态按 §13 报告。
7. 代理误差独立量化，不与 ADS 数值误差混淆。
8. C1-30 与 C1-H 物理域明确区分。
9. 现有 6D/7D focused tests、完整 `fpm test` 和相关 Python 测试通过。
10. 正式输出记录 OPM、config、模式数据、可执行文件和源码 revision 身份。
