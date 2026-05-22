# CRTBP Uncertainty Propagation Design

**Date:** 2026-05-22
**Project:** POD_Fortran

## Overview

Add a Circular Restricted Three-Body Problem (CRTBP) force model and uncertainty propagation capability to POD_Fortran. The CRTBP is a distinct dynamical system — rotating frame, two primaries, prescaled non-dimensional — that does not interact with the existing Earth-centered force model. All new code lives in independent files; zero changes to existing source files.

## Requirements

1. CRTBP equations of motion (non-dimensional, from C++ reference implementation)
2. RKF45 adaptive integrator dedicated to CRTBP (completely independent from `pod_integrator_module`)
3. Monte Carlo (MC) uncertainty propagation
4. Differential Algebra (DA) uncertainty propagation
5. Independent executable `run_CRTBP_uprop` with `--in <json>` CLI
6. All inputs/outputs are CRTBP non-dimensional units; dimensional inputs are rejected
7. DA AlgebraicVector operations use explicit `vec_add/sub/mul/div` subroutines (no operator overloading for vectors) to prevent memory leaks

## New Files

| File | Purpose |
|------|---------|
| `src/lib/forcemodel/pod_crtbp_module.f90` | CRTBP EOM — real and DA derivative subroutines |
| `src/lib/integrator/pod_crtbp_integrator_module.f90` | Standalone RKF45 adaptive integrator (real + DA) |
| `src/lib/uncertainty/propagation/pod_uq_crtbp_mc_module.f90` | MC propagator: sample → integrate → statistics |
| `src/lib/uncertainty/propagation/pod_uq_crtbp_da_module.f90` | DA propagator: single DA integration → polynomial eval |
| `app/run_CRTBP_uprop.f90` | Executable: CLI, JSON parsing, method routing, output |
| `config/crtbp_uprop_example.json` | Example configuration file |

## Data Flow

```
config.json  ──►  run_CRTBP_uprop  ──►  MC or DA propagator
                                              │
                    ┌─────────────────────────┤
                    ▼                         ▼
              MC Propagator              DA Propagator
              (N particles)              (1 DA integration)
                    │                         │
                    ▼                         ▼
         crtbp_integrator_module     crtbp_integrator_module
         (adaptive RKF45, real)      (adaptive RKF45, DA)
                    │                         │
                    ▼                         ▼
           crtbp_module              crtbp_module
           (real EOM)                (DA EOM)
                    │                         │
                    ▼                         ▼
              final samples            compiled polynomial
                                         → eval all samples
                                              │
                    ┌─────────────────────────┘
                    ▼
             CSV output (particles.csv + stats.csv)
```

## Architecture Decisions

- **No reuse of existing integrator module**: The CRTBP EOM are non-dimensional by nature and have no physical unit conversion chain. A standalone integrator with inlined Butcher coefficients is cleaner than parameterizing the existing one.
- **No reuse of existing UQ base class**: CRTBP has no `config%LU/VU/TU`, no SPICE, no `set_propagation_epoch`. The existing `uq_propagator_base` is tightly coupled to these concepts. The CRTBP propagators are self-contained modules, not subclasses.
- **Reuse of existing utilities**: `pod_random_module` for `generate_multivariate_normal`, `pod_dace_classes` for DA types, `pod_global` for `DP`.

## Module Details

### pod_crtbp_module

```
module pod_crtbp_module
  use pod_global, only: DP
  use pod_dace_classes  (for DA version)

  real(DP), public :: crtbp_mu = 0.01_DP
  public :: set_crtbp_mu
  public :: crtbp_derivatives_real
  public :: crtbp_derivatives_da

contains
  subroutine set_crtbp_mu(mu)
  subroutine crtbp_derivatives_real(x, dxdt, t)   ! x(6), dxdt(6), t in
  subroutine crtbp_derivatives_da(x_da, dxdt_da, t, pool)  ! AlgebraicVector in/out
end module
```

The real version is a direct translation of the C++ CRTBP `operator()`. The DA version uses explicit subroutine calls for all AlgebraicVector operations, with a temp pool passed in to avoid allocations inside the hot loop.

### pod_crtbp_integrator_module

```
module pod_crtbp_integrator_module
  use pod_global, only: DP

  public :: rkf45_crtbp_step
  public :: adaptive_integrate_crtbp

contains
  subroutine rkf45_crtbp_step(state, dt, t, state_5th, state_4th, error_est)
    ! Inlined RKF45 Butcher tableau, calls crtbp_derivatives_real
  end subroutine

  subroutine adaptive_integrate_crtbp(state, t_start, t_end, &
      rel_tol, abs_tol, dt_min, dt_max, max_steps, &
      times, states, n_steps, nfe)
    ! Adaptive step control with WRMS error, same algorithm as pod_integrator_module
  end subroutine
end module
```

### pod_uq_crtbp_mc_module

```
module pod_uq_crtbp_mc_module
  use pod_global, only: DP
  public :: crtbp_mc_propagate

contains
  subroutine crtbp_mc_propagate(nominal_state, cov, mu, t_end, &
      n_particles, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
      final_samples, final_mean, final_cov, verbose)
    ! 1. set_crtbp_mu(mu)
    ! 2. Generate particles via generate_multivariate_normal
    ! 3. OpenMP parallel: adaptive_integrate_crtbp per particle
    ! 4. Compute moments (mean, covariance)
  end subroutine
end module
```

### pod_uq_crtbp_da_module

```
module pod_uq_crtbp_da_module
  use pod_global, only: DP
  use pod_dace_classes
  public :: crtbp_da_propagate

contains
  subroutine crtbp_da_propagate(nominal_state, cov, mu, t_end, &
      da_order, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
      final_samples, final_mean, final_cov, propagated_ref, verbose)
    ! 1. set_crtbp_mu(mu)
    ! 2. Build DA initial state: nominal + da_var(i)
    ! 3. Single DA adaptive integration
    ! 4. Compile polynomial, evaluate all deviations
    ! 5. Compute moments
  end subroutine
end module
```

### run_CRTBP_uprop

```
program run_CRTBP_uprop
  ! CLI: --in <json_path>
  ! 1. Parse command-line for --in
  ! 2. Read and parse JSON config file
  ! 3. Validate: reject inputs with non-CRTBP scales (|x| > 100 or |v| > 10 → warning)
  ! 4. Route to crtbp_mc_propagate or crtbp_da_propagate
  ! 5. Write output CSV files
end program
```

## JSON Config Schema

```json
{
    "mu": 0.01,
    "initial_state": [x, y, z, vx, vy, vz],
    "initial_covariance": [[...6x6...]],
    "propagation_time": 10.0,
    "method": "MC",
    "num_samples": 10000,
    "da_order": 2,
    "integrator": {
        "rel_tol": 1.0e-12,
        "abs_tol": 1.0e-12,
        "min_step": 1.0e-10,
        "max_step": 0.1,
        "max_steps": 100000
    },
    "output": {
        "prefix": "./output/crtbp_uprop"
    }
}
```

Fields:
- `mu`: CRTBP mass ratio (required, must be in (0, 1))
- `initial_state`: 6-element non-dimensional state [x,y,z,vx,vy,vz] (required)
- `initial_covariance`: 6×6 non-dimensional covariance matrix (required)
- `propagation_time`: non-dimensional propagation duration (required, > 0)
- `method`: "MC" or "DA" (required)
- `num_samples`: particle count for MC (required)
- `da_order`: DA expansion order, only used when method="DA" (optional, default 2)
- `integrator.rel_tol/abs_tol`: RKF45 tolerances (optional, defaults provided)
- `integrator.min_step/max_step`: step size bounds (optional, defaults provided)
- `integrator.max_steps`: max adaptive steps (optional, default 100000)
- `output.prefix`: path prefix for output CSV files (required)

## Output

Two CSV files per run:
- `<prefix>_particles.csv`: header row + N rows of final state particles
- `<prefix>_stats.csv`: mean vector (1 row) + covariance matrix (6 rows)

## Non-dimensional Validation

If `|initial_state(1:3)| > 100.0` or `|initial_state(4:6)| > 10.0`, the program errors with a message suggesting the user check their units. This is a heuristic guard since CRTBP states are typically O(1).

## Files NOT Modified

- `fpm.toml` — user will add new executable entry
- All existing `src/`, `app/`, `test/` files — zero changes
- `pod_integrator_module`, `pod_force_model_module`, `pod_uq_*` — untouched
