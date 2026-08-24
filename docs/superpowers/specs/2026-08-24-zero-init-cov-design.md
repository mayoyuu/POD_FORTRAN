# ZERO_INIT_COV Design

## Goal

Allow both HALO SRP sweep scripts to set the six-dimensional orbital initial
covariance to exactly zero while leaving the nominal orbital state and the
`eta_srp` uncertainty sweep unchanged.

## User Interface

Both sweep scripts accept the environment variable `ZERO_INIT_COV`:

```bash
ZERO_INIT_COV=1 bash scripts/run_halo_srp_sweep_parallel.sh
ZERO_INIT_COV=1 bash scripts/run_halo_nominal_srp_sweep_parallel.sh
```

The values `1`, `true`, and `TRUE` enable the option. An unset variable or any
other value leaves the existing OPM covariance behavior unchanged.

Both Fortran executables accept a consistent command-line flag:

```text
--zero-init-cov
```

The scripts translate the enabled environment variable into that CLI flag.

## Semantics and Data Flow

1. Each executable loads the epoch, nominal six-dimensional state, and
   covariance from the selected OPM as before.
2. When `--zero-init-cov` is present, only the loaded 6-by-6 orbital covariance
   is replaced by zero.
3. The OPM state mean is not modified.
4. In `run_srp_uq_propagation`, `srp_eta_mean` and `srp_eta_sigma` remain
   independent of this option. The first six sample components are initially
   deterministic, while the seventh component continues to follow the requested
   `eta_srp` distribution.
5. In `run_HFEM_uprop`, all initial samples equal the OPM nominal state when the
   option is enabled.

## Zero-Covariance Sampling

`generate_multivariate_normal` will explicitly recognize an exactly all-zero
covariance matrix. In that case it fills every sample column with the supplied
mean and returns without calling Cholesky factorization. Nonzero covariance
matrices retain the current Cholesky-based behavior, including the existing
diagnostics for invalid non-positive-definite inputs.

Handling the zero case in the shared sampler prevents duplicated special-case
sampling logic and also makes the existing `run_HFEM_uprop -pr 0 -vr 0` input
well-defined.

## Output Isolation and Traceability

When enabled, each script adds `_zeroInitCov` to its output base and run-log
name so that zero-covariance runs cannot overwrite ordinary runs in the same
output root. Existing output names remain unchanged when the option is disabled.

Each batch log records `zero_init_cov`, and each `summary.tsv` row contains a
`zero_init_cov` column. Queue entries use the isolated output paths.

## Error Handling

The scripts follow the existing boolean convention: `1`, `true`, and `TRUE`
enable the mode; all other values disable it. The Fortran flag has no argument,
so it cannot consume or misinterpret the next CLI option.

Zero covariance is treated as a valid deterministic distribution. Other
singular or indefinite covariance matrices are outside this feature's scope and
continue through the current error-diagnostic path.

## Tests

Automated tests will verify:

- the shared multivariate-normal sampler returns the exact mean in every sample
  for an all-zero covariance;
- its existing nonzero-covariance behavior remains available;
- both executables accept `--zero-init-cov` and complete a short propagation;
- the SRP executable retains the requested nonzero `eta_srp` sigma while the
  orbital initial covariance is zero;
- both scripts translate `ZERO_INIT_COV=1` into `--zero-init-cov`, isolate output
  names, and record the mode in summary/log metadata;
- both scripts remain syntactically valid Bash programs.

The implementation will follow a red-green-refactor sequence: add and observe
failing tests first, then make the minimum production changes required to pass.

## Scope

This change does not modify OPM files, change nominal states, change the SRP
force model, alter the default covariance, or add support for arbitrary
positive-semidefinite singular covariance matrices.
