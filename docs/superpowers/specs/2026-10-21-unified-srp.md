# 统一 SRP 联合不确定性传播

用户于 2026-10-21 确认统一入口方案，并补充：除 UT 外，每种需要 MC 的方法每小时保存固定 2000 点；ADS 保留分裂信息。

## 代码复用与交付

初版在本文件所在独立目录暂存；用户随后要求合入原 POD_Fortran 的 main。以原仓库为最终运行位置，不替换其 fpm.toml、不删除 Git 库或 worktree。复用 boxwing-modes-p1 的动力学、积分器、6/7/8D ADS、模式表加载器。一个 app/run.sh 调用 tools/run_srp_unified.py，后者组织样本、构建配置、批处理和 HZ；一个 app/run_srp_unified.f90 负责全部传播方法。旧入口的几何研究/历史兼容功能不删除，不继续复制其逐案例运行循环。

案例：C0 无 SRP；C1 标称 cannonball；C2 对月 box-wing；C3 对地 box-wing；T1 cannonball 加 eta 误差；T2 现有 time_varying_modes 姿态模式代理。T2 与实时 box-wing 模型区别必须写入结果元数据。

## 方法

UT 使用同一输入样本的一二阶矩构造联合 6/7/8D scaled UT（alpha=0.5,beta=2,kappa=0，可覆盖），保留状态与 SRP 参数相关性。sigma 点参数在整个轨迹上固定；逐小时接续原点轨迹，不重新高斯化。越界 sigma 点报错，不裁剪。UT 记录 sigma 点和权重，HZ 不适用于带权确定性 sigma 点。

ADS 复用现有增量历史和分裂模块。MC 对每个真实样本独立实数积分。DAMC 使用相同域映射，但禁用分裂，以 k 阶 DA 映射大量样本。自适应 DAMC 在固定独立 Real 探针的误差超过阈值时增加阶数；因低阶丢失的系数不能恢复，升阶必须从初始历元重建至当前小时。达到最高阶仍超限时保留结果并标无效。

## 样本与输出

一次生成统一 z1..z7,s,u，所有方法/案例共享前六维初始样本，默认 N=50000，逐维正态截断 3 sigma。支持已有标准化 sample CSV 和物理联合样本。初始协方差可以用 state-scale 缩放；OPM 完整 6x6 协方差保留。

每小时保存全体样本均值、完整 6x6 协方差、共同 cannonball nominal 参考、ADS/DA 的独立验证误差和有效性。非 UT 保存固定最多 2000 个样本，数量不足时明确记录 actual_saved_count。UT 保存其原生 sigma 点及 Wm/Wc，禁止用 sigma 点当作随机样本做 HZ。

HZ 默认使用固定保存的最多 2000 点，记录统计量、p 值、alpha、样本数、秩和可用性；0 h 同样检验以揭示输入截断的影响。均值/协方差始终来自完整 N，不从 2000 点重算。HZ 拒绝正态不自动解释成动力学新产生非高斯性。

分项 wall-clock：传播、映射求值、统计、独立验证、HZ、输出、端到端。compute_excluding_io 是实际分项之和；不从不同进程/不同阶段不兼容的计时简单相减。编译和准备时间单列。

ADS splits.csv 保存每个事件的区间、方向、父深度、六分量截断估计/容差/超限量、方向评分、中心和宽度；hourly_moments.csv 保存补丁数与精度有效位。

## 迁移与验收

路径从脚本位置解析，星历和外部库通过 --kernel-dir、DACE_ROOT、OPENBLAS_ROOT、SPICE_LIB_DIR 指定。默认单线程外层作业，可用 --jobs 控制进程并发；进程内状态不跨作业共享。构建一次，案例批量复用同一可执行文件。保存 resolved_config、输入哈希、源码来源、运行命令和逐案例状态。

首轮验收 L1Halo-1，UT/ADS × C0 C1 C2 C3 T1 T2，0..2 h，每小时保存，N=50000、2000 点输出。用户要求测试简短有效并明确取消额外压力测试：不再运行一阶/极严格阈值压力测试，不扩展 24 h 测试。主分支仅补必要短时入口检查；360 h 脚本与 GNU parallel/内置进程池命令由用户自行启动。不得将短时冒烟通过称为 15 天数值精度认证。

## 来源

- Julier & Uhlmann, 2004, Unscented Filtering and Nonlinear Estimation, DOI 10.1109/JPROC.2003.823141。
- Henze & Zirkler, 1990, A class of invariant consistent tests for multivariate normality, DOI 10.1080/03610929008830400。
