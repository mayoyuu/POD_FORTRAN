# CRTBP 平板 θ/φ DA 传播

2026-11-05：按用户确认的太阳切向基接入 θ/φ。此次只实现接口和编译，不进行 u/v、实数轨道或其他数值对照验证，也不启动批量传播。

## 最小改动

保留四个参数槽位、CRTBP RHS、RKF78 和 ADS 传播接口。新增姿态坐标选择：

- 默认 SRP_ATTITUDE_UV：p=[beta,chi,u,v]，兼容旧调用。
- SRP_ATTITUDE_THETA_PHI：p=[beta,chi,theta,phi]，内部角度全部为弧度。

初始六维状态仍为 DA 常数。上下界相等的参数固定，其余参数在区间中点展开：

| 活跃变量 | 设置 | 物理槽位到 DA 自变量的 index_map |
|---|---|---|
| θ、φ | β、χ 上下界相等 | [0,0,1,2] |
| β、θ、φ | χ 上下界相等 | [1,0,2,3] |
| β、χ、θ、φ | 四个参数都有非零区间宽度 | [1,2,3,4] |

对应源码：

- src/lib/forcemodel/pod_srp_study_parameters_module.f90：角度域工厂、模式与参数名称。
- src/lib/forcemodel/pod_crtbp_srp_study_module.f90：real/DA 法向构造分支。
- src/lib/uncertainty/propagation/pod_crtbp_srp_study_ads_module.f90：检查模式一致、分裂轴名称。
- app/run_crtbp_srp_study.f90：可选命令行参数及 JSON 元数据。

积分器、ADS 分裂策略和统计计算本体没有修改。

## 几何与单位

s 从航天器指向太阳，默认参考轴 k=(0,0,1) 位于会合系中：

~~~text
e1 = normalize(k cross s)
e2 = s cross e1
n  = cos(theta)*s + sin(theta)*(cos(phi)*e1 + sin(phi)*e2)
cos_incidence = cos(theta)
~~~

保留现有 tangent_rotation_rad，默认 0；θ/φ 相对于最终 e1/e2 定义。太阳视线和切向基仍包含传播状态的 DA 依赖。θ/φ 的三角函数直接用 DA 运算，不读其常数项替代，不经 u/v 转换。

θ 范围为 [0,pi/2)，避免进入正背面切换。φ 使用连续角区间，不执行取模，例如输入 [350,370] 度表示跨越 360° 的连续范围。θ=0 时 φ 的一阶作用消失，但保留 θ–φ 的混合项。

角度工厂输入度、内部区间存弧度。直接设置 domain%lower/upper 或调用力模型时必须使用弧度。均值/协方差沿用现有参数矩形域的均匀采样；角度模式意味着 θ、φ 各自在所设区间均匀，并非圆锥内等立体角分布。

## Fortran 接入

使用现有模块声明 cfg、domain、params 等对象后，只需更改域和模式：

~~~fortran
use pod_srp_study_parameters_module
use pod_crtbp_srp_study_module

! 示例：beta=.1、chi=.6 固定，仅 theta/phi 为 DA 变量。
! 角度范围只是接口示例，并非已选定的全量实验配置。
domain = srp_study_angle_domain( &
    beta_lower=0.1_DP, beta_upper=0.1_DP, &
    theta_lower_deg=4.0_DP, theta_upper_deg=6.0_DP, &
    phi_lower_deg=25.0_DP, phi_upper_deg=35.0_DP, &
    chi_lower=0.6_DP, chi_upper=0.6_DP)
cfg%model = STUDY_FLAT_PLATE
cfg%attitude_coordinates = domain%attitude_coordinates
~~~

直接 DA/RKF78 路径继续使用：

~~~fortran
! 初始化 DACE 前必须释放上一套 DACE 对象。
call dace_initialize(6, domain%n_variables())
call srp_study_init_da(domain, params, index_map)
! 后续既有 RHS/RKF78 调用不变。
~~~

若调用 study_ads_propagate，该函数已经负责 DACE 和根域初始化，应直接传入 cfg、domain，不提前创建另一套带活动句柄的 DA 对象。ADS 原来对状态和四个参数共同进行仿射替换的流程继续使用。

令 beta_lower/beta_upper 不同即可加入 β；再令 chi_lower/chi_upper 不同即可加入 χ。domain%parameter_names() 返回本次所用的物理参数名；旧常量 SRP_PARAMETER_NAMES 仍专用于 u/v 兼容接口。

## 单工况命令行

原来 19 个位置参数保持原意。启用角度模式时，在输出文件路径后追加五项：

~~~text
theta_phi theta_lower_deg theta_upper_deg phi_lower_deg phi_upper_deg
~~~

例如仅 θ/φ 不确定的调用（本次未执行）：

~~~bash
source setup_env.sh
fpm run --compiler gfortran-13 --profile debug --flag -ffree-line-length-none \
  --target run_crtbp_srp_study -- \
  1.1309265107780351 0 0 0 -0.46540743845849059 0 \
  0.012153614091891635 2.3017923284002024 0 \
  flatplate 0.1 0.1 0 0.6 0.6 6 1e-8 8 \
  build/theta_phi_case.json theta_phi 4 6 25 35
~~~

这里第 13 项 tilt_max_deg=0 是旧接口位置占位，在角度模式下忽略；角度以最后四个界限为准。太阳初相为第 9 项，示例取 0 rad。传播时间仍为 2T；若不足 15 天则为 3T。

输出增加：

- attitude_coordinates：uv 或 theta_phi。
- parameter_names、parameter_units：例如 [beta,chi,theta,phi] 与 [m2/kg,1,rad,rad]。
- parameter_lower/upper：实际根域物理界限，角度为弧度。
- parameter_da_index：固定参数为 0，活跃参数为对应 DA 索引。
- parameter_prior：现有参数矩形域均匀采样的明确标记。
- 每个 split event 的 axis_name：与已有 ax 对应。ax、split_count_dir 仍按四个物理槽位排序。

事件中的 c/w/call/wall 保留原有归一化坐标含义，不直接作为角度值。物理值由根域中点和半宽换算。后处理应读取 parameter_names，不能在角度模式中继续把第 3/4 个物理槽位标为 u/v。

原入口在 beta_lower=0 时输出的 beta_sensitivity 仍是正对日（θ=0）、χ=.6 的既有单边差分参考，不是该角度区间的平均敏感度。本次未改这条诊断路径。

## Python 批量运行 Q2 / Q4（2026-11-06）

tools/run_crtbp_srp_study_batch.py 现已接入角度模式。Q2/Q4 保留各自的 beta、chi 区间，使用以下参数选择姿态：

| 选项 | 含义 |
|---|---|
| --attitude theta_phi | 平板配置使用 theta/phi；默认仍为 uv |
| --theta-deg LOWER UPPER | 可选；默认 theta=[0,配置原 tilt_deg] |
| --phi-deg LOWER UPPER | 角度模式必填，连续方位角区间，单位度 |

例如在 POD_Fortran 目录运行以下命令，仅预览 Q2/Q4 各一个工况：

~~~bash
python3 tools/run_crtbp_srp_study_batch.py \
  --configs Q2_L2_t10 Q4_L3_chi01 \
  --attitude theta_phi --theta-deg 0 10 --phi-deg -45 45 \
  --solar-phases 0 --limit 2 --dry-run \
  --outdir build/q2_q4_angles_preview
~~~

dry-run 只打印命令，不创建结果目录或改写 manifest，也不编译或调用 Fortran。未指定 --executable 时打印的 run_crtbp_srp_study 是占位名称。

正式批次可用同样选项，去掉 --dry-run 和 --limit，并指定 --jobs 及新的结果目录：

~~~bash
source setup_env.sh
FPM_FC=gfortran-13 FPM_FFLAGS=-ffree-line-length-none \
python3 tools/run_crtbp_srp_study_batch.py \
  --configs Q2_L2_t10 Q4_L3_chi01 \
  --attitude theta_phi --theta-deg 0 10 --phi-deg -45 45 \
  --solar-phases 0 --jobs 4 \
  --outdir results/q2_q4_theta_phi
~~~

以上角度只是调用示例。也可以用 --orbits / --phases 筛选指定轨道与轨道初相；--solar-phases 始终为太阳初相（弧度），与姿态 phi 无关。

| 配置 | beta (m2/kg) | chi | DA 活跃变量 |
|---|---|---|---|
| Q2_L2_t10 | [0.08,0.12] | 固定 0.6 | beta,theta,phi |
| Q4_L3_chi01 | [0.08,0.12] | [0,1] | beta,chi,theta,phi |

Q2_L2_ext_main 与 Q4_L3_ext 对应 beta=[0,0.1]；Q2_L2_ext_ext 对应 beta=[0,0.2]。原 --beta 选项专用于自定义 Q1 炮弹球，不能用它覆盖 Q2/Q4。更改 beta/chi 配置时可直接修改 CONFIGS 中相应条目。

不提供 --theta-deg 时，Q2_L2_t5/t10/t20 分别使用 [0,5]/[0,10]/[0,20] 度。提供时则统一覆盖所有选中平板配置。phi 允许例如 [350,370] 度，不取模；theta 必须满足 0<=lower<=upper<90。

角度工况 ID 会包含 theta_phi 与实际角度区间，例如 Q2_L2_t10__theta_phi_t0_10_p-45_45。配置、角度区间和结果坐标元数据保存到 manifest/campaign 中，防止与旧 u/v 结果混用。若 --executable 指向不支持角度模式的旧程序，结果会标为 configuration_mismatch，需重新编译。

本次完成 Python 语法、参数校验与 dry-run 命令生成检查；检查中禁止子进程，没有启动 fpm、Fortran、传播或任何数值对照。
