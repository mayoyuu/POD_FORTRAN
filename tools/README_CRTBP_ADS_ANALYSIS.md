# CRTBP + SRP：首次分裂时间与初始参数域盒子

所有命令在 POD_Fortran/ 目录执行，只读取已保存 JSON，不编译、不传播。
图中的盒子是**初始参数域的 ADS 划分**，不是航天器位置分布。
依赖：Python 3.9+；PNG 需要 matplotlib，可执行 python3 -m pip install matplotlib。
仅导出 CSV/JSON 可加 --no-png，不需要第三方包。

## 1. H1S 首次分裂汇总

以已运行的 Q2 theta=[0,10] deg、phi=[-45,45] deg、太阳初相0为例：

~~~bash
RESULTS=results/srp_sensitivity/results/q2_theta_phi/results
CONFIG=Q2_L2_t10__theta_phi_t0_10_p-45_45
ANALYSIS=results/srp_sensitivity/analysis/H1S

python3 tools/analyze_srp_study_progress.py \
  --results-root "$RESULTS" \
  --configs "$CONFIG" --orbit-prefix H1S --solar-phase 0 \
  --outdir "$ANALYSIS"
~~~

批量传播的 --outdir 下还有一层 results/，里面才是配置目录。
本脚本也接受批次目录并自动进入其 results/。新工况替换 RESULTS、CONFIG 即可。
用 ls "$RESULTS" 查看实际配置 ID；省略 --configs 汇总目录内所有配置。

输出：

- split_stats_H1S.csv：逐工况首次分裂无量纲时间、天数、分裂次数、分裂方向、发散数、终止天数、源文件路径。
- heatmap_first_split_time_<CONFIG>.png：振幅×轨道初相的首次分裂天数。
- heatmap_n_splits_total_<CONFIG>.png 和各方向热图。
- direction_summary_H1S.csv：每配置/太阳初相/轨道族的方向汇总。

原始 first_split_time=-1 表示该传播弧段内未记录分裂；CSV天数字段留空，不能当0天。
保留 app_status 和 n_diverged，发生发散的工况也保留其已记录的首次分裂。
theta_phi 模式输出 split_theta/split_phi；旧 split_u/split_v 列留空。
角度配置目录中残留的旧 u/v 输出会打印 configuration mismatch 警告并排除。
--solar-phase 是太阳初相，单位弧度，与姿态 phi 无关。省略时按太阳相位分别生成热图，不混为一个单元格。
PNG 灰色包括无分裂和缺失单元格，请以逐工况 CSV 判断。
新分析会覆盖同名汇总文件；不同批次需分别保存时，将 ANALYSIS 设置为不同子目录。

## 2. 第14天指定工况的 theta–phi 盒子

以下选择 H1S_A005000/p00__a0.json，即轨道初相索引0、太阳初相0；可换为实际工况。
Q2_L2_t10 的 beta 根域为 [0.08,0.12] m2/kg。

~~~bash
CASE="$RESULTS/$CONFIG/H1S_A005000/p00__a0.json"

python3 tools/plot_srp_ads_partition.py \
  --json "$CASE" --time-days 14 \
  --axes theta phi --slice beta=0.1 \
  --outdir "$ANALYSIS"
~~~

自动创建 $ANALYSIS/$CONFIG/H1S_A005000/。
文件名包含源文件名、天数、坐标轴和切片值：

- partition_*_boxes.csv：该时刻全部完整参数盒子的边界、归一化中心/完整宽度、深度和 in_slice。
- partition_*_summary.json：源路径、查询时间、已应用事件数、切片参数、分区性质。
- partition_*.png：所选二维切片盒子边界。

CSV 的 theta/phi 物理边界用弧度（列名带 _rad）；CLI --slice theta=... / phi=... 以及图的角度轴用度。
未指定的隐藏参数固定为根域中点。落在隐藏维度的分裂边界时取右侧盒子；根域最右端点包含在内，避免重复切片。
CSV 所有盒子均保留，in_slice=1 标记图中显示的盒子。

查看 beta–theta 切片（phi固定0度）：

~~~bash
python3 tools/plot_srp_ads_partition.py \
  --json "$CASE" --time-days 14 \
  --axes beta theta --slice phi=0 \
  --outdir "$ANALYSIS"
~~~

Q4 的 theta–phi 切片可显式固定两个隐藏参数：
--axes theta phi --slice beta=0.1 chi=0.6。
两个绘图轴必须具有非零不确定区间。

查看全部高维盒子的二维投影，可用 --projection 替代 --slice ...。
投影会重叠，不能把重叠矩形当成互不相交的二维划分。

## 3. 可恢复范围

- 由完整 split_events 的 call/wall 从根域回放二叉均分；wall 是完整宽度，不能误当半宽。
- 使用不晚于查询时刻的分裂记录；检查点之间划分保持不变。这是记录的 ADS 检查点时间，不是连续时间的精确阈值穿越时刻。
- 事件时间只写出7位有效数字，边界比较容差写入 summary.json。
- 未发生盒子删除时，得到准确的活跃参数域划分。
- 查询时刻已发生发散删除时，现有 JSON 没有删除盒子的边界。默认仍绘制历史分裂分区，但打印警告，PNG/CSV/JSON 均标记 historical，不能判断每个盒子是否仍存活。加 --require-active 可直接报错。
- 首次发散时间缺失而 n_diverged>0 时，也保守标记 historical。
- 缺少 call/wall、事件数与总数不符、父盒无法匹配等情况报错，不猜测盒子。
- 查询14天超出该工况 t_end_days 时直接报错，不会外推。
- 不能恢复沿轨道传播后的状态分布；那需要相应时刻的 DA 映射或状态样本。

## 验证

python3 -m unittest discover -s test -p test_crtbp_ads_partition.py

仅使用人工构造的日志，检查时间切片、根域体积守恒、缺失日志、发散标记、Q2/Q4隐藏维度、角度转换与多太阳相位分组，不执行 Fortran 或数值传播对照。
