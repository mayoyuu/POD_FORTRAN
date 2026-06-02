# Noise Configuration Spec: EMDAC & UT Runners

**Date:** 2026-06-02
**Status:** approved

## Motivation

Both `pod_emdac_runner_module.f90` and `pod_ut_runner_module.f90` hardcode measurement noise (R matrix) and process noise (Q matrix scale `sigma_a`). Changing these requires source edits and recompilation. The goal is to make them configurable through the existing `config.txt` key-value system.

## New Config Fields

Four fields added to `config_params` in `pod_config_module.f90`:

| Field                          | Type      | Default     | Description                              |
|--------------------------------|-----------|-------------|------------------------------------------|
| `measurement_noise_ra_arcsec`  | real(DP)  | 0.1         | RA measurement noise std (arcseconds)    |
| `measurement_noise_dec_arcsec` | real(DP)  | 0.1         | Dec measurement noise std (arcseconds)   |
| `process_noise_sigma_a`        | real(DP)  | 1.0e-11     | Process noise acceleration std (km/s^2)  |
| `use_process_noise`            | logical   | .true.      | Enable/disable process noise Q matrix    |

Measurement noise has no separate enable switch — set both values to 0.0 to effectively disable it.

## Calculation

Measurement noise R (2x2 diagonal):
```
noise_R(1,1) = (meas_noise_ra_arcsec  * PI / 180.0 / 3600.0)^2
noise_R(2,2) = (meas_noise_dec_arcsec * PI / 180.0 / 3600.0)^2
```

Process noise Q (6x6 diagonal):
```
if use_process_noise:
    noise_Q(i,i)     = (dt^4 / 4) * sigma_a^2     for i=1,2,3 (position)
    noise_Q(i+3,i+3) = dt^2 * sigma_a^2           for i=1,2,3 (velocity)
else:
    noise_Q = 0
```

## Files Changed

| File | Change |
|------|--------|
| `src/lib/system/pod_config_module.f90` | Add 4 fields to `config_params` type, defaults, parsing, serialization, printing |
| `src/functions/orbitimprove/pod_emdac_runner_module.f90` | `use pod_config`; read noise values from `config%` instead of hardcoded |
| `src/functions/orbitimprove/pod_ut_runner_module.f90` | Same as emdac_runner |
| `config/config.txt` | Add new keys under `# 测量模型参数` |
| `config/ut_config.txt` | Same as config.txt |

## Data Flow

```
config.txt  -->  pod_config%load_config  -->  config (global, type config_params)
                                                    |
                    emdac_runner / ut_runner  <------'  (use pod_config, only: config)
```

No new function parameters. Runners access `config%xxx` directly.

## Backward Compatibility

- Existing `measurement_noise_std` field in config is retained but unused by runners (deprecated).
- Default values match the current hardcoded values (0.1 arcsec, 1.0e-11 km/s^2, process noise enabled).
- Config files without the new keys will use defaults — existing configs continue to work unchanged.
