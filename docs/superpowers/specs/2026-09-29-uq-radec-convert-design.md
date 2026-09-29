# run_uq_radec 传播结果到 RA/DEC 平面后处理程序设计

日期：2026-09-29

## 1. 目标与范围

新增后处理程序 `run_uq_radec`，把 `run_HFEM_uprop` 的传播结果统一转换到站心 RA/DEC 平面：

- 粒子 CSV（`x,y,z,vx,vy,vz` 或 7 列 `…,eta_srp`）→ `ra,dec` CSV；
- 矩 JSON（`mean`[6] + `covariance`[6×6]）→ RA/DEC 的 `mean`[2] + `covariance`[2×2]，用无迹变换（UT）完成。

参考 `app/run_radec_convert.f90`（已能把粒子 CSV 转 RA/DEC，但不处理矩 JSON，且读 7 列文件会失败）。本程序是它的超集，覆盖 MC/DA/UT/ADS 四种方法的输出。

本程序不改动 `run_HFEM_uprop` 的传播逻辑，只新增一个独立可执行程序与一个 JSON 读取子程序。

## 2. 命令行接口

```
fpm run run_uq_radec -- -i <particles.csv> -j <moments.json> -s <station_id> -et <epoch> \
    [-o <out_prefix>] [--site <site.json>] [--deg] [-cfg <config.txt>]
```

- `-i` 粒子 CSV（可选）；`-j` 矩 JSON（可选）。**至少给一个**。
- `-s` 测站 ID（必填）；`-et` 观测历元（必填，字符串，经 `str2et` 解析）。
- `--site` 测站 JSON（默认 `config/site-used.json`）；`-cfg` 配置文件（默认 `config/config.txt`）。
- `--deg` 角度输出为度；默认输出为弧度。
- 只给 `-i` → 仅 CSV 路径；只给 `-j` → 仅 JSON 路径；都给 → 两条路径都执行。

输出命名（未给 `-o` 时从各输入名派生；给 `-o <prefix>` 时统一用 prefix）：

- CSV 输出：`<csv_stem>_radec.csv` 或 `<prefix>_radec.csv`；
- JSON 输出：`<json_stem>_radec_moments.json` 或 `<prefix>_radec_moments.json`。

## 3. 数据流

### 3.1 CSV 路径（粒子 → RA/DEC）

1. 读表头后逐行解析。每行按空格/逗号分隔读取数值；列数必须为 6 或 7（取前 3 列作为 J2000 位置，忽略速度与 `eta_srp`）。列数非 6/7 的行跳过并计数告警。
2. 用 `compute_measurement(state, et, station, 'OPTICAL', measurement)` 得到 `ra,dec`（弧度）。
3. 按 `--deg` 决定是否乘 `180/π`，逐行写入 `ra,dec`。

### 3.2 JSON 路径（矩 → RA/DEC，无迹变换）

1. `read_uq_moments_json` 读出 `mean`[6] 与 `covariance`[6×6]。
2. `generate_sigma_points(mean, cov, kappa=0.5, sigmas_x, weights)`：6 维生成 `2·6+1 = 13` 个 sigma 点。
3. 每个 sigma 点（前 3 维位置）过 `compute_measurement` 得到对应的 `ra,dec`，组成 `sigmas_y(2, 13)`。
4. `reconstruct_ut_moments(sigmas_x, sigmas_y, weights, mean_y, P_yy, P_xy, is_angle=[.true.,.false.])`：
   - `is_angle(1)=.true.` 对 RA 用圆均值（`atan2(Σw·sin, Σw·cos)`）并修正 0/2π 跨界；`is_angle(2)=.false.` 对 DEC 用线性均值。
   - `P_xy`（交叉协方差）本程序不输出，但调用接口保留。
5. 输出 `mean_radec`[2] 与 `cov_radec`[2×2]。`--deg` 时均值乘 `180/π`，协方差乘 `(180/π)²`。

## 4. 模块改动

在 `src/lib/data/pod_uq_sample_io_module.f90` 新增一个公开子程序，与已有的 `write_uq_moments_json` 对称：

```fortran
subroutine read_uq_moments_json(filename, mean, cov, status, message)
    character(len=*), intent(in) :: filename
    real(DP), allocatable, intent(out) :: mean(:)   ! 6
    real(DP), allocatable, intent(out) :: cov(:,:)  ! 6x6
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
```

行式解析（本库三种 `_moments.json` 形状的公共字段 `mean`/`covariance` 格式一致）：

- 含 `"mean"` 的行：取出 `[` 与 `]` 之间的 6 个数（表控读入，逗号视为分隔符）。
- 含 `"covariance"` 键的行：其后连续 6 行各含一组 `[ ... ]`，逐行取 6 个数填入协方差行。

`status=0` 成功；非 0 时 `message` 说明失败原因（打不开、缺 `mean`/`covariance`、数值不足、非有限值等）。

## 5. 输出格式

CSV 表头与 `run_radec_convert` 一致：

- 弧度：`ra_rad,dec_rad`；度：`ra_deg,dec_deg`。每行 `ES22.14,ES22.14`。

JSON 输出：

```json
{
  "method": "radec_ut",
  "mean": [ra, dec],
  "covariance": [[c_ra_ra, c_ra_dec], [c_dec_ra, c_dec_dec]]
}
```

单位：弧度（`--deg` 时为度）。数值用 `ES24.16E3`，与 `write_uq_moments_json` 一致。

## 6. 错误处理

- `-i` 与 `-j` 均未提供 → 报错退出。
- `-s` / `-et` 缺失 → 报错退出。
- 输入文件打不开 / 为空 / 无有效数据行 → 报错退出。
- CSV 列数非 6/7 → 跳过该行，累计告警；全跳过则报错。
- `read_uq_moments_json` 失败 → 报错退出。
- Cholesky 分解失败（协方差非正定，如近零协方差）→ 退化处理：仅用均值点求 `ra,dec`，协方差记为零，并打印告警。

## 7. 验证

新增 `test/test_uq_radec_io.f90`，并注册到 `fpm.toml`：

1. **JSON 读取往返**：构造已知 `mean`/`cov` 写临时 `_moments.json`，`read_uq_moments_json` 读回，逐元素比对一致。
2. **CSV 6/7 列读取**：6 列与 7 列粒子文件均能解析，输出行数与输入行数一致。
3. **UT 一致性**：对 `output/TD1_2604_2_TO_260612_times_100_*` 样例，JSON 路径得到的 `mean_radec` 与粒子 CSV 逐点 `compute_measurement` 后的经验均值近似吻合（允许 UT 与经验矩的正常差异）；`cov_radec` 半正定。
4. **角度单位**：`--deg` 与默认弧度的换算关系正确（含协方差 `(180/π)²` 缩放）。
5. **RA 跨界**：构造 RA 跨 0/2π 的合成 sigma 点，验证圆均值不出现 180° 错误。

手动验收：对 `output/TD1_2604_2_TO_260612_times_100_particles.csv` + `_moments.json` 用 `-s R91 -et <该历元>` 跑通，产出 CSV 与 JSON 且可解析。

## 8. 已知限制

- `-et` 需用户显式给出传播终态历元；`_moments.json` 当前不记录历元，程序无法自行推断。
- UT 仅传播到二阶矩（`mean`+`cov`），不产出 RA/DEC 的偏度/峰度。
- `--deg` 输出时 RA 单位为度且约定在 `[0, 360)`，DEC 在 `[-90, 90]`；圆均值修正只在 RA 维度启用。
- 测站模型沿用 `run_radec_convert` 的 `preload_stations` / `find_station_by_id` / `compute_measurement`，不新增测量类型。
