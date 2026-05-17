# Orbit Error Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add reference orbit error analysis to EMDAC and UT runners, computing position/velocity errors, RMS, and Mahalanobis distance at each measurement epoch and writing results to `.err` files.

**Architecture:** I/O (reader + writer) stays in existing modules (`pod_obs_io_module`, `pod_data_format_module`). Error computation goes into a new `pod_error_analysis_module` under `src/lib/math/`. Both runners switch from per-observation file scanning to bulk preloading for OBS, stations, and reference orbits, then compute errors inline between time update and measurement update.

**Tech Stack:** Fortran 2003+, SPICE (time conversion), LAPACK (matrix inverse via `pod_basicmath_module`), DACE (DA propagation in EMDAC)

---

## File Structure

| File | Responsibility |
|------|---------------|
| `src/lib/measurement/pod_obs_io_module.f90` | OBS/REF record types, preload subroutines, station lookup |
| `src/lib/math/pod_error_analysis_module.f90` (new) | `compute_orbit_error` — position/velocity error, RMS, Mahalanobis distance |
| `src/lib/data/pod_data_format_module.f90` | `write_error_line` — `.err` file output with `is_first` pattern |
| `src/functions/orbitimprove/pod_emdac_runner_module.f90` | EMDAC runner — new interface, preloading, error computation |
| `src/functions/orbitimprove/pod_ut_runner_module.f90` | UT runner — new interface, preloading, error computation |

---

### Task 1: Add `ref_orbit_record` type and `preload_reference_orbits` to `pod_obs_io_module`

**Files:**
- Modify: `POD_Fortran/src/lib/measurement/pod_obs_io_module.f90`

- [ ] **Step 1: Add `ref_orbit_record` type and public declaration**

Add after line 12 (`public :: obs_record, station_record`), change the public list to include the new type and subroutine:

```fortran
    public :: obs_record, station_record, ref_orbit_record
```

Add after line 14 (`public :: preload_observations, preload_stations, find_station_by_id`):

```fortran
    public :: preload_observations, preload_stations, find_station_by_id
    public :: preload_reference_orbits
```

Add the new type after the `station_record` type block (after line 26):

```fortran
    !> 参考轨道记录（与 obs_record 对称）
    type :: ref_orbit_record
        real(DP) :: et              ! 历元 (TDB 秒)
        real(DP) :: state(6)        ! [X, Y, Z, Vx, Vy, Vz] (km, km/s)
    end type ref_orbit_record
```

- [ ] **Step 2: Add `preload_reference_orbits` subroutine at end of module (before `end module`)**

```fortran
    !> ---------------------------------------------------------------
    !> 一次性将 ORBITS_REF 文件全部解析到内存
    !> 格式: Idx Name Sys YYYY MM DD HH MM SS.SSSSSS X Y Z Vx Vy Vz
    !> ---------------------------------------------------------------
    subroutine preload_reference_orbits(ref_file, ref_list)
        character(len=*), intent(in) :: ref_file
        type(ref_orbit_record), allocatable, intent(out) :: ref_list(:)

        integer :: u, ios, n, i, idx
        character(len=MAX_STRING_LEN) :: line
        character(len=24) :: obj_name, time_sys
        integer :: year, month, day, hour, min
        real(DP) :: sec, X, Y, Z, Vx, Vy, Vz
        character(len=40) :: utc_str

        ! 1. 统计行数
        open(newunit=u, file=ref_file, status='old', iostat=ios)
        if (ios /= 0) stop "[ERROR] 无法打开 REF 文件: " // trim(ref_file)
        n = 0
        do
            read(u, '(A)', iostat=ios) line
            if (ios < 0) exit
            if (ios > 0) stop "[ERROR] 读取 REF 文件出错"
            n = n + 1
        end do
        rewind(u)

        allocate(ref_list(n))

        ! 2. 逐行解析
        do i = 1, n
            read(u, '(A)') line
            read(line, *, iostat=ios) idx, obj_name, time_sys, &
                                       year, month, day, hour, min, sec, &
                                       X, Y, Z, Vx, Vy, Vz
            if (ios /= 0) then
                write(*,*) "[ERROR] 解析 REF 行失败，行号：", i
                stop
            end if

            ref_list(i)%state(1) = X
            ref_list(i)%state(2) = Y
            ref_list(i)%state(3) = Z
            ref_list(i)%state(4) = Vx
            ref_list(i)%state(5) = Vy
            ref_list(i)%state(6) = Vz

            write(utc_str, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",F12.6)') &
                  year, month, day, hour, min, sec
            ref_list(i)%et = utc_to_et(trim(utc_str))
        end do
        close(u)

    end subroutine preload_reference_orbits
```

- [ ] **Step 3: Commit**

```bash
git add POD_Fortran/src/lib/measurement/pod_obs_io_module.f90
git commit -m "feat: add ref_orbit_record type and preload_reference_orbits"
```

---

### Task 2: Create `pod_error_analysis_module` with `compute_orbit_error`

**Files:**
- Create: `POD_Fortran/src/lib/math/pod_error_analysis_module.f90`

- [ ] **Step 1: Create the module file**

```fortran
!> @file pod_error_analysis_module.f90
!> @brief 轨道误差分析：计算估计状态相对于参考轨道的各项误差指标
module pod_error_analysis_module
    use pod_global, only: DP
    use pod_basicmath_module, only: inverse_and_determinant

    implicit none
    private

    public :: compute_orbit_error

contains

    !> ======================================================================
    !> 计算估计状态相对于参考轨道（真值）的误差
    !>
    !> 计算时机：时间更新之后、测量更新之前
    !> ======================================================================
    subroutine compute_orbit_error(est_state, est_cov, ref_state, &
                                   pos_err, vel_err, pos_rms, vel_rms, mahalanobis_d)
        real(DP), intent(in)  :: est_state(6)   ! 滤波器估计 [X,Y,Z,Vx,Vy,Vz]
        real(DP), intent(in)  :: est_cov(6,6)   ! 滤波器协方差 P_{k|k-1}
        real(DP), intent(in)  :: ref_state(6)   ! 参考轨道（真值）
        real(DP), intent(out) :: pos_err(3)     ! dX, dY, dZ (km)
        real(DP), intent(out) :: vel_err(3)     ! dVx, dVy, dVz (km/s)
        real(DP), intent(out) :: pos_rms        ! sqrt(dX^2 + dY^2 + dZ^2) (km)
        real(DP), intent(out) :: vel_rms        ! sqrt(dVx^2 + dVy^2 + dVz^2) (km/s)
        real(DP), intent(out) :: mahalanobis_d  ! sqrt(dx^T * P^{-1} * dx)

        real(DP) :: dx(6), inv_cov(6,6), det_cov
        integer  :: info

        ! 1. 位置与速度误差
        pos_err = est_state(1:3) - ref_state(1:3)
        vel_err = est_state(4:6) - ref_state(4:6)

        ! 2. 位置/速度 RMS
        pos_rms = sqrt(sum(pos_err**2))
        vel_rms = sqrt(sum(vel_err**2))

        ! 3. 马氏距离 sqrt(dx^T * P^{-1} * dx)
        dx = est_state - ref_state
        call inverse_and_determinant(est_cov, inv_cov, det_cov, info)

        if (info == 0) then
            mahalanobis_d = sqrt(dot_product(dx, matmul(inv_cov, dx)))
        else
            mahalanobis_d = -1.0_DP  ! 标记：协方差奇异，无法计算马氏距离
        end if

    end subroutine compute_orbit_error

end module pod_error_analysis_module
```

- [ ] **Step 2: Commit**

```bash
git add POD_Fortran/src/lib/math/pod_error_analysis_module.f90
git commit -m "feat: add pod_error_analysis_module with compute_orbit_error"
```

---

### Task 3: Add `write_error_line` to `pod_data_format_module`

**Files:**
- Modify: `POD_Fortran/src/lib/data/pod_data_format_module.f90`

- [ ] **Step 1: Add `write_error_line` to the public list**

Change line 15 from:
```fortran
    public :: load_initial_opm
    public :: write_json_opm, write_residual_line
```
To:
```fortran
    public :: load_initial_opm
    public :: write_json_opm, write_residual_line, write_error_line
```

- [ ] **Step 2: Add `write_error_line` subroutine after `write_residual_line` (after line 276)**

```fortran
    !> ======================================================================
    !> 写入单行轨道误差数据 (对称于 write_residual_line)
    !> ======================================================================
    subroutine write_error_line(filename, et, pos_err, vel_err, &
                                 pos_rms, vel_rms, mahalanobis_d, is_first)
        character(len=*), intent(in) :: filename
        real(DP), intent(in)         :: et
        real(DP), intent(in)         :: pos_err(3), vel_err(3)
        real(DP), intent(in)         :: pos_rms, vel_rms, mahalanobis_d
        logical, intent(in)          :: is_first

        integer :: u, ios
        character(len=64) :: utc_str
        character(len=256) :: full_filename

        full_filename = trim(filename) // ".err"

        call et2utc(et, 'ISOC', 3, utc_str)

        if (is_first) then
            open(newunit=u, file=full_filename, status='replace', action='write', iostat=ios)
            write(u, '(A)') "# UTC_Time                 ET_Seconds     dX(km)       dY(km)       dZ(km)     &
             dVx(km/s)     dVy(km/s)     dVz(km/s)    Pos_RMS(km)  Vel_RMS(km/s)  Mahal_D"
        else
            open(newunit=u, file=full_filename, status='old', position='append', action='write', iostat=ios)
        end if

        if (ios /= 0) return

        write(u, '(A24, F16.4, 3F14.6, 3F14.8, 2F14.6, F14.6)') &
            utc_str, et, &
            pos_err(1), pos_err(2), pos_err(3), &
            vel_err(1), vel_err(2), vel_err(3), &
            pos_rms, vel_rms, mahalanobis_d

        close(u)
    end subroutine write_error_line
```

- [ ] **Step 3: Commit**

```bash
git add POD_Fortran/src/lib/data/pod_data_format_module.f90
git commit -m "feat: add write_error_line for .err file output"
```

---

### Task 4: Update `pod_emdac_runner_module` — new interface, preloading, error computation

**Files:**
- Modify: `POD_Fortran/src/functions/orbitimprove/pod_emdac_runner_module.f90`

- [ ] **Step 1: Update `use` statements at top of module**

Change lines 6-11:
```fortran
    use pod_filter_emdac_module, only: emdac_filter
    use pod_obs_io_module, only: load_single_observation
    use pod_measurement_base_module, only: observation_station
    use pod_basicmath_module, only: PI
    ! 假设引入了 JSON 读写模块
    use pod_data_format_module, only: load_initial_opm, write_json_opm, write_residual_line
```

To:
```fortran
    use pod_filter_emdac_module, only: emdac_filter
    use pod_obs_io_module, only: obs_record, preload_observations, &
                                  station_record, preload_stations, find_station_by_id, &
                                  ref_orbit_record, preload_reference_orbits
    use pod_measurement_base_module, only: observation_station
    use pod_basicmath_module, only: PI
    use pod_data_format_module, only: load_initial_opm, write_json_opm, &
                                       write_residual_line, write_error_line
    use pod_error_analysis_module, only: compute_orbit_error
```

- [ ] **Step 2: Update subroutine signature**

Replace lines 24-27:
```fortran
    subroutine run_emdac_orbit_determination(obs_file, site_json_file, gmm_in_switch, &
                                             initial_json_file, output_file_name, n_components,&
                                             max_da_order,opt_particles, &
                                             opt_em_max_iter, opt_em_tol)
```

With:
```fortran
    subroutine run_emdac_orbit_determination(obs_file, site_json_file, ref_orbit_file, &
                                             initial_json_file, output_opm_file, &
                                             output_residual_file, output_error_file, &
                                             gmm_in_switch, n_components, max_da_order, &
                                             opt_particles, opt_em_max_iter, opt_em_tol)
```

- [ ] **Step 3: Update intent declarations**

Replace lines 29-34:
```fortran
        character(len=*), intent(in) :: obs_file           ! 观测文件路径 (.obs)
        character(len=*), intent(in) :: site_json_file     ! 测站配置文件路径 (.json)
        character(len=*), intent(in) :: initial_json_file  ! 初始先验状态文件路径 (.opm/.json)
        character(len=*), intent(in) :: output_file_name   ! 输出定轨结果文件路径 (.opm/.json)
        logical, intent(in) :: gmm_in_switch               ! GMM 初始化开关
        integer,  intent(in) :: n_components       ! GMM 分量数量
        integer,  intent(in) :: max_da_order       ! DA 阶数
```

With:
```fortran
        character(len=*), intent(in) :: obs_file             ! 观测数据文件 (.obs)
        character(len=*), intent(in) :: site_json_file       ! 测站配置文件 (.json)
        character(len=*), intent(in) :: ref_orbit_file       ! 参考轨道文件 (.ref)
        character(len=*), intent(in) :: initial_json_file    ! 初始先验状态文件 (.json)
        character(len=*), intent(in) :: output_opm_file      ! 输出定轨结果文件
        character(len=*), intent(in) :: output_residual_file ! 输出残差文件
        character(len=*), intent(in) :: output_error_file    ! 输出误差文件 (.err)
        logical, intent(in) :: gmm_in_switch
        integer,  intent(in) :: n_components
        integer,  intent(in) :: max_da_order
```

- [ ] **Step 4: Add preload declarations and remove `load_single_observation` variables**

In the local variable section, replace:
```fortran
        type(observation_station) :: current_station
```
With:
```fortran
        type(observation_station) :: current_station
        type(obs_record), allocatable :: obs_list(:)
        type(station_record), allocatable :: station_list(:)
        type(ref_orbit_record), allocatable :: ref_list(:)
```

- [ ] **Step 5: Add error metrics local variables**

Add after existing locals:
```fortran
        real(DP) :: pos_err(3), vel_err(3)
        real(DP) :: pos_rms, vel_rms, mahalanobis_d
```

- [ ] **Step 6: Replace `is_eof` with loop counter and add preload calls**

Remove `is_eof` from the logical declaration. Then replace the loop section (after filter init, around line 97):
```fortran
        ! 4. 核心数据同化流 (Time Update + Measurement Update)
        obs_count = 1
        do
            call load_single_observation(obs_file, site_json_file, obs_count, &
                                         et_obs, y_meas(1), y_meas(2), current_station, is_eof)
            if (is_eof) exit
```

With:
```fortran
        ! 4. 预加载全部观测、测站与参考轨道
        call preload_observations(obs_file, obs_list)
        call preload_stations(site_json_file, station_list)
        call preload_reference_orbits(ref_orbit_file, ref_list)

        write(*,*) '  [Runner] 预加载完成：观测 ', size(obs_list), &
                   ' 条，测站 ', size(station_list), ' 个，参考轨道 ', size(ref_list), ' 条'

        ! 5. 核心数据同化流 (Time Update + Measurement Update)
        do obs_count = 1, size(obs_list)
            et_obs = obs_list(obs_count)%et
            y_meas(1) = obs_list(obs_count)%ra
            y_meas(2) = obs_list(obs_count)%dec
            current_station = find_station_by_id(obs_list(obs_count)%station_id, station_list)
```

- [ ] **Step 7: Close the loop properly**

Change `end do` to match the `do obs_count = 1, size(obs_list)` loop. The remaining body (time update, measurement update, residual writing) stays the same except for the additions in Step 8.

- [ ] **Step 8: Insert error computation between time update and measurement update**

After `call my_filter%get_current_state(final_mean)` (post-time-update) and before `call my_filter%measurement_update`, add:

```fortran
            call my_filter%get_current_cov(final_cov)

            call compute_orbit_error(final_mean, final_cov, ref_list(obs_count)%state, &
                                      pos_err, vel_err, pos_rms, vel_rms, mahalanobis_d)

            call write_error_line(output_error_file, et_obs, pos_err, vel_err, &
                                   pos_rms, vel_rms, mahalanobis_d, (obs_count == 1))
```

- [ ] **Step 9: Update `write_residual_line` call to use new filename parameters**

Change:
```fortran
            call write_residual_line(output_file_name, et_obs, y_meas, step_comp, &
                                    step_res, trim(current_station%name), (obs_count == 1))
```
To:
```fortran
            call write_residual_line(output_residual_file, et_obs, y_meas, step_comp, &
                                    step_res, trim(current_station%name), (obs_count == 1))
```

- [ ] **Step 10: Update `write_json_opm` call**

Change:
```fortran
        call write_json_opm(output_file_name, final_mean, final_cov, my_filter%gmm_state, 0.0_DP, "DRO", et_current)
```
To:
```fortran
        call write_json_opm(output_opm_file, final_mean, final_cov, my_filter%gmm_state, 0.0_DP, "DRO", et_current)
```

- [ ] **Step 11: Remove obsolete code**

Remove `is_eof` from the logical declaration line and any remaining references to it.

- [ ] **Step 12: Commit**

```bash
git add POD_Fortran/src/functions/orbitimprove/pod_emdac_runner_module.f90
git commit -m "feat: update EMDAC runner with preloading, error analysis, new interface"
```

---

### Task 5: Update `pod_ut_runner_module` — new interface, preloading, error computation

**Files:**
- Modify: `POD_Fortran/src/functions/orbitimprove/pod_ut_runner_module.f90`

- [ ] **Step 1: Update `use` statements**

Change lines 5-9:
```fortran
    use pod_filter_ut_module, only: ut_filter          ! 你的 UT 滤波器模块
    use pod_obs_io_module, only: load_single_observation
    use pod_measurement_base_module, only: observation_station
    use pod_basicmath_module, only: PI
    use pod_data_format_module, only: load_initial_opm, write_json_opm, write_residual_line
```

To:
```fortran
    use pod_filter_ut_module, only: ut_filter
    use pod_obs_io_module, only: obs_record, preload_observations, &
                                  station_record, preload_stations, find_station_by_id, &
                                  ref_orbit_record, preload_reference_orbits
    use pod_measurement_base_module, only: observation_station
    use pod_basicmath_module, only: PI
    use pod_data_format_module, only: load_initial_opm, write_json_opm, &
                                       write_residual_line, write_error_line
    use pod_error_analysis_module, only: compute_orbit_error
```

- [ ] **Step 2: Update subroutine signature**

Replace:
```fortran
    subroutine run_ut_orbit_determination(obs_file, site_json_file, initial_json_file, &
                                          output_file_name)
        character(len=*), intent(in) :: obs_file           ! 观测文件 (.obs)
        character(len=*), intent(in) :: site_json_file     ! 测站配置文件 (.json)
        character(len=*), intent(in) :: initial_json_file  ! 初始先验状态文件 (.json)
        character(len=*), intent(in) :: output_file_name   ! 输出结果文件 (.json)
```

With:
```fortran
    subroutine run_ut_orbit_determination(obs_file, site_json_file, ref_orbit_file, &
                                          initial_json_file, output_opm_file, &
                                          output_residual_file, output_error_file)
        character(len=*), intent(in) :: obs_file
        character(len=*), intent(in) :: site_json_file
        character(len=*), intent(in) :: ref_orbit_file
        character(len=*), intent(in) :: initial_json_file
        character(len=*), intent(in) :: output_opm_file
        character(len=*), intent(in) :: output_residual_file
        character(len=*), intent(in) :: output_error_file
```

- [ ] **Step 3: Add preload data structures and error variables to locals**

Add after `type(observation_station) :: current_station`:
```fortran
        type(obs_record), allocatable :: obs_list(:)
        type(station_record), allocatable :: station_list(:)
        type(ref_orbit_record), allocatable :: ref_list(:)
```

Add after `real(DP) :: step_res(6)`:
```fortran
        real(DP) :: pos_err(3), vel_err(3)
        real(DP) :: pos_rms, vel_rms, mahalanobis_d
```

Remove `logical :: is_eof` from declarations.

- [ ] **Step 4: Replace the main loop — preload + iterate**

Replace lines 61-109 (from `obs_count = 1` to `end do` of the main read loop):

```fortran
        ! ---- 4. 预加载全部观测、测站与参考轨道 ----
        call preload_observations(obs_file, obs_list)
        call preload_stations(site_json_file, station_list)
        call preload_reference_orbits(ref_orbit_file, ref_list)

        write(*,*) '  [UT Runner] 预加载完成：观测 ', size(obs_list), &
                   ' 条，测站 ', size(station_list), ' 个，参考轨道 ', size(ref_list), ' 条'

        ! ---- 5. 主循环：逐观测进行时间更新 + 测量更新 ----
        do obs_count = 1, size(obs_list)
            et_obs = obs_list(obs_count)%et
            y_meas(1) = obs_list(obs_count)%ra
            y_meas(2) = obs_list(obs_count)%dec
            current_station = find_station_by_id(obs_list(obs_count)%station_id, station_list)

            ! 当前滤波器时刻
            call my_filter%get_current_epoch(et_current)
            dt = et_obs - et_current

            ! 构造与时间步长相关的过程噪声 Q
            noise_Q = 0.0_DP
            do i = 1, 3
                noise_Q(i,i)       = (dt**4 / 4.0_DP) * sigma_a**2
                noise_Q(i+3,i+3)   = dt**2 * sigma_a**2
            end do

            write(*,*) '  [UT Runner] 处理观测 #', obs_count, &
                    '  观测时刻 et_obs = ', et_obs, ' 秒', &
                    '  时间步长 dt =', dt, ' 秒'

            call my_filter%time_update(et_obs, noise_Q)
            call my_filter%get_current_epoch(et_current)
            call my_filter%get_current_state(final_mean)
            call my_filter%get_current_cov(final_cov)

            call compute_orbit_error(final_mean, final_cov, ref_list(obs_count)%state, &
                                      pos_err, vel_err, pos_rms, vel_rms, mahalanobis_d)

            call write_error_line(output_error_file, et_obs, pos_err, vel_err, &
                                   pos_rms, vel_rms, mahalanobis_d, (obs_count == 1))

            call my_filter%measurement_update(y_meas, noise_R, et_obs, current_station)
            call my_filter%get_current_epoch(et_current)
            call my_filter%get_current_state(final_mean)

            call my_filter%get_last_residual(step_res, step_comp)

            call write_residual_line(output_residual_file, et_obs, y_meas, step_comp, &
                                    step_res, trim(current_station%name), (obs_count == 1))

        end do
```

- [ ] **Step 5: Update `write_json_opm` call**

Change:
```fortran
        call write_json_opm(output_file_name, final_mean, final_cov, &
                    rms = 0.0_DP, obj_id = "DRO", et_last = et_current)

        write(*,*) '  [UT Runner] 结果已写入：', trim(output_file_name)
```
To:
```fortran
        call write_json_opm(output_opm_file, final_mean, final_cov, &
                    rms = 0.0_DP, obj_id = "DRO", et_last = et_current)

        write(*,*) '  [UT Runner] 结果已写入：', trim(output_opm_file)
```

- [ ] **Step 6: Commit**

```bash
git add POD_Fortran/src/functions/orbitimprove/pod_ut_runner_module.f90
git commit -m "feat: update UT runner with preloading, error analysis, new interface"
```

---

### Task 6: Update `pod_emdac_test.f90` to match new runner interface

**Files:**
- Modify: `POD_Fortran/test/pod_emdac_test.f90`

- [ ] **Step 1: Add new file path variables**

After `character(len=MAX_STRING_LEN) :: obs_file, initial_json_file, output_file_name, site_json_file`, add:

```fortran
    character(len=MAX_STRING_LEN) :: ref_orbit_file, output_residual_file, output_error_file
```

- [ ] **Step 2: Add new command-line argument parsing**

After the existing `select case` block (after the `-gmm` case), add new cases inside the loop:

```fortran
            case ('-ref', '--ref_orbit')
                call get_command_argument(i+1, arg_str); ref_orbit_file = trim(arg_str)
                i = i + 1
            case ('-res', '--residual')
                call get_command_argument(i+1, arg_str); output_residual_file = trim(arg_str)
                i = i + 1
            case ('-err', '--error')
                call get_command_argument(i+1, arg_str); output_error_file = trim(arg_str)
                i = i + 1
```

- [ ] **Step 3: Set default values for new parameters**

After `output_file_name = 'output/DROb_202601_1_emdac_result_with_process_noise'`, add:

```fortran
    ref_orbit_file      = 'ORBITS_REF/DRO/DRO_single_R91_1h.ref'
    output_residual_file = 'output/DROb_202601_1_emdac_residual'
    output_error_file    = 'output/DROb_202601_1_emdac_error'
```

- [ ] **Step 4: Update the runner call to match new interface**

Replace the existing `call run_emdac_orbit_determination(...)` block with:

```fortran
    call run_emdac_orbit_determination( &
        obs_file             = obs_file, &
        site_json_file       = site_json_file, &
        ref_orbit_file       = ref_orbit_file, &
        initial_json_file    = initial_json_file, &
        output_opm_file      = output_file_name, &
        output_residual_file = output_residual_file, &
        output_error_file    = output_error_file, &
        gmm_in_switch        = gmm_in_switch, &
        opt_particles        = opt_particles, &
        max_da_order         = opt_da_order, &
        opt_em_max_iter      = opt_em_max_iter, &
        opt_em_tol           = opt_em_tol, &
        n_components         = n_components)
```

- [ ] **Step 5: Commit**

```bash
git add POD_Fortran/test/pod_emdac_test.f90
git commit -m "test: update EMDAC test for new runner interface with error analysis"
```

---

### Task 7: Update `pod_ut_test.f90` to match new runner interface

**Files:**
- Modify: `POD_Fortran/test/pod_ut_test.f90`

- [ ] **Step 1: Add new file path variables**

After `character(len=MAX_STRING_LEN) :: obs_file, initial_json_file, output_file_name, site_json_file`, add:

```fortran
    character(len=MAX_STRING_LEN) :: ref_orbit_file, output_residual_file, output_error_file
```

- [ ] **Step 2: Set default values for new parameters**

After `output_file_name  = 'output/DROb_202601_ut_result'`, add:

```fortran
    ref_orbit_file      = 'ORBITS_REF/DRO/DRO_single_R91_1h.ref'
    output_residual_file = 'output/DROb_202601_ut_residual'
    output_error_file    = 'output/DROb_202601_ut_error'
```

- [ ] **Step 3: Update the runner call to match new interface**

Replace the existing `call run_ut_orbit_determination(...)` block with:

```fortran
    call run_ut_orbit_determination( &
        obs_file             = obs_file, &
        site_json_file       = site_json_file, &
        ref_orbit_file       = ref_orbit_file, &
        initial_json_file    = initial_json_file, &
        output_opm_file      = output_file_name, &
        output_residual_file = output_residual_file, &
        output_error_file    = output_error_file)
```

- [ ] **Step 4: Commit**

```bash
git add POD_Fortran/test/pod_ut_test.f90
git commit -m "test: update UT test for new runner interface with error analysis"
```

---

### Task 8: Build and test with fpm

- [ ] **Step 1: Clean previous build artifacts**

```bash
cd POD_Fortran
fpm clean
```

Expected: Clean build directory.

- [ ] **Step 2: Source environment and build**

```bash
source ./setup_env.sh
fpm build
```

Expected: All modules compile successfully. Look for errors in `pod_error_analysis_module`, `pod_obs_io_module`, `pod_data_format_module`, `pod_emdac_runner_module`, `pod_ut_runner_module`.

- [ ] **Step 3: Run EMDAC test**

```bash
fpm test pod_emdac_test
```

Expected: Runner initializes, preloads OBS/REF/station data, runs filter loop, writes `.opm.json`, `.residual`, and `.err` files to `output/`.

- [ ] **Step 4: Verify output files**

```bash
head -5 output/DROb_202601_1_emdac_error.err
```

Expected: Header line + 4 data rows showing dX, dY, dZ, dVx, dVy, dVz, Pos_RMS, Vel_RMS, Mahal_D columns.

- [ ] **Step 5: Run UT test**

```bash
fpm test pod_ut_test
```

Expected: Same as EMDAC — `.opm.json`, `.residual`, `.err` files generated.

- [ ] **Step 6: Verify UT output**

```bash
head -5 output/DROb_202601_ut_error.err
```

Expected: Same format as EMDAC error file.

- [ ] **Step 7: Commit**

```bash
git add .
git commit -m "test: verify build and test pass with error analysis feature"
```
