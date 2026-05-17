# Orbit Error Analysis — 设计规格

**日期**: 2026-05-16
**目标**: 为 EMDAC 和 UT runner 增加参考轨道误差分析功能

## 1. 功能概述

根据输入的 ORBITS_REF 参考轨道文件（真值），在每次测量更新之前，计算滤波器估计状态相对于真实轨道的误差，输出到 `.err` 文件。

### 计算指标

| 指标 | 说明 |
|------|------|
| dX, dY, dZ | 各方向位置误差 (km) |
| dVx, dVy, dVz | 各方向速度误差 (km/s) |
| Pos_RMS | sqrt(dX² + dY² + dZ²) (km) |
| Vel_RMS | sqrt(dVx² + dVy² + dVz²) (km/s) |
| Mahalanobis_D | sqrt(dxᵀ · P⁻¹ · dx)，使用 P_{k|k-1} |

### 计算时机

时间更新之后、测量更新之前（与残差写入位置一致）。

---

## 2. 架构概览

```
ORBITS_REF 文件 ──preload──→ ref_orbit_record 数组
OBS 文件 ────────preload──→ obs_record 数组（已有）

Runner 主循环:
  for i = 1, N:
    1. 从 obs_list(i) 取观测
    2. 从 ref_list(i) 取参考轨道
    3. time_update
    4. compute_orbit_error()  ← 新增
    5. write_error_line()     ← 新增
    6. measurement_update
    7. write_residual_line()  ← 已有
```

### I/O 与计算职责分离

| 职责 | 位置 | 方式 |
|------|------|------|
| 参考轨道数据结构 + 预加载 | `pod_obs_io_module` | 新增 `ref_orbit_record`, `preload_reference_orbits` |
| 误差计算 | `pod_error_analysis_module` (新建) | `compute_orbit_error` |
| .err 文件写入 | `pod_data_format_module` | 新增 `write_error_line` |
| Runner 集成 | `pod_emdac_runner_module`, `pod_ut_runner_module` | 预加载 + 循环内调用 |

---

## 3. 模块详细设计

### 3.1 pod_obs_io_module — 新增类型与读取

```fortran
type :: ref_orbit_record
    real(DP) :: et
    real(DP) :: state(6)  ! X, Y, Z, Vx, Vy, Vz (km, km/s)
end type ref_orbit_record

subroutine preload_reference_orbits(ref_file, ref_list)
    ! 读入 .ref 文件全部行，解析为 ref_orbit_record 数组
    ! 格式: Name Sys YYYY MM DD HH MM SS.SSSS X Y Z Vx Vy Vz
    ! 时间通过 utc_to_et 转为 TDB 秒
```

- 时间匹配容差：不对参考轨道做时间插值，按索引一一对应
- 调用方保证 OBS 与 REF 行数一致；不做行数强制校验

### 3.2 pod_error_analysis_module (新建)

路径：`POD_Fortran/src/lib/math/pod_error_analysis_module.f90`

```fortran
subroutine compute_orbit_error(est_state, est_cov, ref_state, &
                                pos_err, vel_err, pos_rms, vel_rms, mahalanobis_d)
    real(DP), intent(in)  :: est_state(6), est_cov(6,6), ref_state(6)
    real(DP), intent(out) :: pos_err(3), vel_err(3)
    real(DP), intent(out) :: pos_rms, vel_rms, mahalanobis_d
```

计算步骤：
1. `pos_err = est_state(1:3) - ref_state(1:3)`
2. `vel_err = est_state(4:6) - ref_state(4:6)`
3. `pos_rms = sqrt(sum(pos_err**2))`
4. `vel_rms = sqrt(sum(vel_err**2))`
5. `dx = est_state - ref_state`；对 `est_cov` 求逆 (调用 `inverse_and_determinant`)；`mahalanobis_d = sqrt(dot_product(dx, matmul(inv_cov, dx)))`

依赖：`pod_basicmath_module`（求逆）

### 3.3 pod_data_format_module — 新增写入

```fortran
subroutine write_error_line(filename, et, pos_err, vel_err, &
                             pos_rms, vel_rms, mahalanobis_d, is_first)
```

- 首次 (`is_first=.true.`)：`status='replace'`，写表头
- 后续：`status='old', position='append'`
- 输出文件：`trim(filename) // ".err"`
- 格式：

```
# UTC_Time                 ET_Seconds     dX      dY      dZ     dVx     dVy     dVz   Pos_RMS  Vel_RMS  Mahal_D
2027-01-21T07:12:00.000    123456789.000  -1.234   0.567  ...    ...     ...     ...    1.234    0.567    2.345
```

### 3.4 Runner 接口变更

EMDAC：

```fortran
subroutine run_emdac_orbit_determination(          &
    obs_file,           &  ! 观测数据 (.obs)
    site_json_file,     &  ! 测站配置 (.json)
    ref_orbit_file,     &  ! 参考轨道 (.ref)        —— 新增
    initial_json_file,  &  ! 初始先验状态
    output_opm_file,    &  ! 输出定轨结果           —— 原名 output_file_name
    output_residual_file,& ! 输出残差 (.residual)   —— 新增
    output_error_file,  &  ! 输出误差 (.err)        —— 新增
    gmm_in_switch, n_components, max_da_order, &
    opt_particles, opt_em_max_iter, opt_em_tol)
```

UT：

```fortran
subroutine run_ut_orbit_determination(          &
    obs_file,           &  ! 观测数据 (.obs)
    site_json_file,     &  ! 测站配置 (.json)
    ref_orbit_file,     &  ! 参考轨道 (.ref)        —— 新增
    initial_json_file,  &  ! 初始先验状态
    output_opm_file,    &  ! 输出定轨结果
    output_residual_file,& ! 输出残差 (.residual)   —— 新增
    output_error_file)    ! 输出误差 (.err)         —— 新增
```

**关键变更：预加载替代逐行扫描**

循环前：
```fortran
call preload_observations(obs_file, obs_list)
call preload_stations(site_json_file, station_list)
call preload_reference_orbits(ref_orbit_file, ref_list)  ! 新增
```

循环内不再调用 `load_single_observation`（逐行打开扫描），直接从 `obs_list(i)` 取数据，用 `find_station_by_id` 查测站。

utc_to_et 调用移至预加载内部完成，循环内不再重复转换。

---

## 4. 数据流

```
初始: OBS_File ──preload──→ obs_list[]
       Site_Json ──preload──→ station_list[]
       REF_File ──preload──→ ref_list[]      ← 新增

循环 (i=1..N):
  obs_list(i) ──→ et_obs, y_meas(ra,dec), station_id
  find_station_by_id(station_list, station_id) ──→ current_station
  ref_list(i) ──→ ref_state(6)                 ← 新增

  time_update(et_obs, Q) → P_{k|k-1}, state_mean
  compute_orbit_error(state_mean, P_{k|k-1}, ref_state) → err metrics  ← 新增
  write_error_line(...)                                                ← 新增

  measurement_update(y_meas, R, ...) → corrected state
  write_residual_line(...)  (已有)
```

---

## 5. 文件变更汇总

| 文件 | 操作 | 说明 |
|------|------|------|
| `POD_Fortran/src/lib/math/pod_error_analysis_module.f90` | **新建** | 误差计算模块 |
| `POD_Fortran/src/lib/measurement/pod_obs_io_module.f90` | 修改 | 新增 `ref_orbit_record`, `preload_reference_orbits`，公开 `preload_observations`, `obs_record`, `preload_stations`, `station_record`, `find_station_by_id` |
| `POD_Fortran/src/lib/data/pod_data_format_module.f90` | 修改 | 新增 `write_error_line`，公开 |
| `POD_Fortran/src/functions/orbitimprove/pod_emdac_runner_module.f90` | 修改 | 新接口参数 + 预加载 + 误差计算调用 |
| `POD_Fortran/src/functions/orbitimprove/pod_ut_runner_module.f90` | 修改 | 新接口参数 + 预加载 + 误差计算调用 |

---

## 6. 错误处理

- `preload_reference_orbits`：文件打开失败则 `stop`，与 `preload_observations` 风格一致
- `compute_orbit_error`：协方差矩阵求逆异常由 `inverse_and_determinant` 内部处理
- `write_error_line`：文件写入失败静默退出（`iostat /= 0` 时 `return`），与 `write_residual_line` 一致
- OBS/REF 行数不匹配：不做强制校验，由调用方确保一致
