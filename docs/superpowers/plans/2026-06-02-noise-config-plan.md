# Noise Configuration for EMDAC & UT Runners — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make measurement noise (RA/Dec) and process noise (sigma_a, enable switch) configurable via config.txt for both emdac_runner and ut_runner.

**Architecture:** Add 4 fields to the existing `config_params` type in `pod_config_module.f90`. Runners import `pod_config` and read values from the global `config` variable directly — no new function parameters. Existing defaults match current hardcoded values for zero behavioral change.

**Tech Stack:** Fortran 2003+, fpm build system

---

## File Map

| File | Role |
|------|------|
| `src/lib/system/pod_config_module.f90` | Config type definition, parsing, serialization, printing |
| `src/functions/orbitimprove/pod_emdac_runner_module.f90` | EMDAC-N runner — reads noise from config |
| `src/functions/orbitimprove/pod_ut_runner_module.f90` | UT runner — reads noise from config |
| `config/config.txt` | EMDAC config file — new keys added |
| `config/ut_config.txt` | UT config file — new keys added |

---

### Task 1: Add noise fields to config infrastructure

**File:** `src/lib/system/pod_config_module.f90`

- [ ] **Step 1: Add 4 fields to `config_params` type**

After line 140 (`real(DP) :: outlier_threshold = 3.0_DP`), insert:

```fortran
        ! 滤波器噪声参数
        real(DP) :: measurement_noise_ra_arcsec = 0.1_DP     ! RA 测量噪声标准差 (角秒)
        real(DP) :: measurement_noise_dec_arcsec = 0.1_DP    ! Dec 测量噪声标准差 (角秒)
        real(DP) :: process_noise_sigma_a = 1.0e-11_DP       ! 过程噪声加速度标准差 (km/s²)
        logical  :: use_process_noise = .true.                ! 是否启用过程噪声
```

- [ ] **Step 2: Add defaults in `set_default_config`**

After `config%outlier_threshold = 3.0_DP` (line 291 in `set_default_config`), insert:

```fortran
        config%measurement_noise_ra_arcsec  = 0.1_DP
        config%measurement_noise_dec_arcsec = 0.1_DP
        config%process_noise_sigma_a        = 1.0e-11_DP
        config%use_process_noise            = .true.
```

- [ ] **Step 3: Add case branches in `set_config_value`**

After `case ('outlier_threshold')` block (line 597), insert:

```fortran
            case ('measurement_noise_ra_arcsec')
                read(value, *, iostat=ios) config%measurement_noise_ra_arcsec
            case ('measurement_noise_dec_arcsec')
                read(value, *, iostat=ios) config%measurement_noise_dec_arcsec
            case ('process_noise_sigma_a')
                read(value, *, iostat=ios) config%process_noise_sigma_a
            case ('use_process_noise')
                config%use_process_noise = (trim(value) == 'true')
```

- [ ] **Step 4: Add write lines in `save_default_config`**

After `write(unit, '(A)') 'outlier_threshold = 3.0'` (line 421), insert:

```fortran
        write(unit, '(A)') 'measurement_noise_ra_arcsec = 0.1'
        write(unit, '(A)') 'measurement_noise_dec_arcsec = 0.1'
        write(unit, '(A)') 'process_noise_sigma_a = 1.0e-11'
        write(unit, '(A)') 'use_process_noise = true'
        write(unit, '(A)') ''
```

- [ ] **Step 5: Add print statements in `print_config`**

Replace the existing measurement model print block (lines 759-762):
```fortran
        write(*, *) '测量模型参数:'
        write(*, *) '  测量噪声标准差: ', config%measurement_noise_std
        write(*, *) '  使用测量权重: ', config%use_measurement_weights
        write(*, *) '  异常值阈值: ', config%outlier_threshold
```

With:
```fortran
        write(*, *) '测量模型参数:'
        write(*, *) '  测量噪声标准差(旧): ', config%measurement_noise_std
        write(*, *) '  使用测量权重: ', config%use_measurement_weights
        write(*, *) '  异常值阈值: ', config%outlier_threshold
        write(*, *) ''
        write(*, *) '滤波器噪声参数:'
        write(*, *) '  RA 测量噪声 (角秒): ', config%measurement_noise_ra_arcsec
        write(*, *) '  Dec 测量噪声 (角秒): ', config%measurement_noise_dec_arcsec
        write(*, *) '  过程噪声 sigma_a (km/s²): ', config%process_noise_sigma_a
        write(*, *) '  启用过程噪声: ', config%use_process_noise
```

- [ ] **Step 6: Build and verify compilation**

Run: `cd POD_Fortran && fpm build`
Expected: compilation succeeds with no errors.

- [ ] **Step 7: Commit**

```bash
git add src/lib/system/pod_config_module.f90
git commit -m "feat: add noise config fields to pod_config module

Add measurement_noise_ra_arcsec, measurement_noise_dec_arcsec,
process_noise_sigma_a, and use_process_noise to config_params type.
Include parsing, serialization, defaults, and print support."
```

---

### Task 2: Update EMDAC runner to read noise from config

**File:** `src/functions/orbitimprove/pod_emdac_runner_module.f90`

- [ ] **Step 1: Add `use pod_config` import**

Replace line 4:
```fortran
    use pod_global, only: DP, MAX_STRING_LEN
```

With:
```fortran
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_config, only: config
```

- [ ] **Step 2: Remove hardcoded `sigma_a` parameter**

Delete line 62:
```fortran
        real(DP), parameter :: sigma_a = 1.0e-11_DP   ! km/s²
```

- [ ] **Step 3: Replace hardcoded `noise_R` with config values**

Replace lines 84-86:
```fortran
        noise_R = 0.0_DP
        noise_R(1,1) = (0.1_DP * PI / 180.0_DP / 3600.0_DP)**2 
        noise_R(2,2) = noise_R(1,1)
```

With:
```fortran
        noise_R = 0.0_DP
        noise_R(1,1) = (config%measurement_noise_ra_arcsec * PI / 180.0_DP / 3600.0_DP)**2
        noise_R(2,2) = (config%measurement_noise_dec_arcsec * PI / 180.0_DP / 3600.0_DP)**2
```

- [ ] **Step 4: Gate Q computation with `use_process_noise`**

Replace lines 158-164:
```fortran
            noise_Q = 0.0_DP
            do i = 1, 3
                noise_Q(i,i)       = (dt**4 / 4.0_DP) * sigma_a**2
                noise_Q(i+3,i+3)   = dt**2 * sigma_a**2
                ! noise_Q(i,i+3)     = (dt**3 / 2.0_DP) * sigma_a**2
                ! noise_Q(i+3,i)     = noise_Q(i,i+3)       ! 对称
            end do
```

With:
```fortran
            noise_Q = 0.0_DP
            if (config%use_process_noise) then
                do i = 1, 3
                    noise_Q(i,i)     = (dt**4 / 4.0_DP) * config%process_noise_sigma_a**2
                    noise_Q(i+3,i+3) = dt**2 * config%process_noise_sigma_a**2
                end do
            end if
```

- [ ] **Step 5: Build and verify compilation**

Run: `cd POD_Fortran && fpm build`
Expected: compilation succeeds with no errors.

- [ ] **Step 6: Commit**

```bash
git add src/functions/orbitimprove/pod_emdac_runner_module.f90
git commit -m "feat: emdac_runner reads measurement and process noise from config

Replace hardcoded sigma_a and noise_R with config% values.
Add use_process_noise gate for Q matrix computation."
```

---

### Task 3: Update UT runner to read noise from config

**File:** `src/functions/orbitimprove/pod_ut_runner_module.f90`

- [ ] **Step 1: Add `use pod_config` import**

Replace line 4:
```fortran
    use pod_global, only: DP, MAX_STRING_LEN
```

With:
```fortran
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_config, only: config
```

- [ ] **Step 2: Remove hardcoded `sigma_a` parameter**

Delete line 49:
```fortran
        real(DP), parameter :: sigma_a = 1.0e-11_DP   ! km/s²
```

- [ ] **Step 3: Replace hardcoded `noise_R` with config values**

Replace lines 58-60:
```fortran
        noise_R = 0.0_DP
        noise_R(1,1) = (0.1_DP * PI / 180.0_DP / 3600.0_DP)**2
        noise_R(2,2) = noise_R(1,1)
```

With:
```fortran
        noise_R = 0.0_DP
        noise_R(1,1) = (config%measurement_noise_ra_arcsec * PI / 180.0_DP / 3600.0_DP)**2
        noise_R(2,2) = (config%measurement_noise_dec_arcsec * PI / 180.0_DP / 3600.0_DP)**2
```

- [ ] **Step 4: Gate Q computation with `use_process_noise`**

Replace lines 99-103:
```fortran
            noise_Q = 0.0_DP
            do i = 1, 3
                noise_Q(i,i)       = (dt**4 / 4.0_DP) * sigma_a**2
                noise_Q(i+3,i+3)   = dt**2 * sigma_a**2
            end do
```

With:
```fortran
            noise_Q = 0.0_DP
            if (config%use_process_noise) then
                do i = 1, 3
                    noise_Q(i,i)     = (dt**4 / 4.0_DP) * config%process_noise_sigma_a**2
                    noise_Q(i+3,i+3) = dt**2 * config%process_noise_sigma_a**2
                end do
            end if
```

- [ ] **Step 5: Build and verify compilation**

Run: `cd POD_Fortran && fpm build`
Expected: compilation succeeds with no errors.

- [ ] **Step 6: Commit**

```bash
git add src/functions/orbitimprove/pod_ut_runner_module.f90
git commit -m "feat: ut_runner reads measurement and process noise from config

Replace hardcoded sigma_a and noise_R with config% values.
Add use_process_noise gate for Q matrix computation."
```

---

### Task 4: Add new noise keys to config files

**Files:** `config/config.txt`, `config/ut_config.txt`

- [ ] **Step 1: Update `config/config.txt`**

Replace the `# 测量模型参数` section:
```
# 测量模型参数
measurement_noise_std = 1.0e-3
use_measurement_weights = true
outlier_threshold = 3.0
```

With:
```
# 测量模型参数
measurement_noise_std = 1.0e-3
use_measurement_weights = true
outlier_threshold = 3.0

# 滤波器噪声参数
measurement_noise_ra_arcsec = 0.1
measurement_noise_dec_arcsec = 0.1
process_noise_sigma_a = 1.0e-11
use_process_noise = true
```

- [ ] **Step 2: Update `config/ut_config.txt`**

Same replacement as Step 1.

- [ ] **Step 3: Commit**

```bash
git add config/config.txt config/ut_config.txt
git commit -m "feat: add noise config keys to config files

Add measurement_noise_ra_arcsec, measurement_noise_dec_arcsec,
process_noise_sigma_a, and use_process_noise to both config files."
```

---

### Task 5: Verification — build and smoke test

- [ ] **Step 1: Full clean build**

```bash
cd POD_Fortran && fpm build
```
Expected: clean build, no errors.

- [ ] **Step 2: Run existing tests**

```bash
cd POD_Fortran && fpm test
```
Expected: all existing tests pass (default noise values match old hardcoded values).

- [ ] **Step 3: Quick smoke test with EMDAC runner**

```bash
cd POD_Fortran && fpm run run_pod_emdac -- \
  -obs example/example_emdac_orbit_determination.f90 \
  -init example/example_emdac_orbit_determination.f90 \
  -cfg config/config.txt 2>&1 | head -5
```
The runner should load the config and print the noise parameters. (It will fail on invalid obs/init files, but the config loading should succeed.)
Expected: config is loaded, noise parameters appear in printed config.

- [ ] **Step 4: Quick smoke test with UT runner**

```bash
cd POD_Fortran && fpm run run_pod_ut -- \
  -obs example/example_ut_orbit_determination.f90 \
  -init example/example_ut_orbit_determination.f90 \
  -cfg config/ut_config.txt 2>&1 | head -5
```
Expected: config loads, noise parameters appear in printed config.
```
