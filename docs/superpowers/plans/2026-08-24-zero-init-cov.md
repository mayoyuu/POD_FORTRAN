# ZERO_INIT_COV Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an exact deterministic six-dimensional initial-orbit mode to both HALO SRP sweep workflows without changing the independent `eta_srp` uncertainty.

**Architecture:** Both Fortran CLIs accept `--zero-init-cov` and clear only the loaded 6-by-6 OPM covariance. The shared multivariate-normal sampler treats an exactly zero covariance as a deterministic distribution, while both Bash drivers translate `ZERO_INIT_COV=1` into the CLI flag and isolate output names.

**Tech Stack:** Fortran 2008, FPM, LAPACK-backed sampling, Bash, Git.

## Global Constraints

- `ZERO_INIT_COV` clears only the six-dimensional orbital covariance.
- `srp_eta_mean` and `SRP_SIGMA_LIST` behavior must remain unchanged.
- Existing behavior and output names remain unchanged when the option is disabled.
- Exact zero covariance must bypass Cholesky; other singular covariance behavior is unchanged.
- Preserve the user's unrelated modifications to `rename_to_pod.sh` and `setup_env.sh`.

---

### Task 1: Make zero covariance a valid deterministic distribution

**Files:**
- Create: `test/test_zero_covariance_sampling.f90`
- Modify: `src/lib/statistics/pod_random_module.f90:46-122`

**Interfaces:**
- Consumes: `generate_multivariate_normal(mean, cov, samples)`
- Produces: The same procedure, extended so `cov == 0` fills every sample with `mean`.

- [ ] **Step 1: Write the failing sampler test**

Create `test/test_zero_covariance_sampling.f90`:

```fortran
program test_zero_covariance_sampling
    use pod_global, only: DP
    use pod_random_module, only: init_random_seed, generate_multivariate_normal
    implicit none

    integer, parameter :: dim = 6, n_samples = 4
    real(DP) :: mean(dim), cov(dim, dim), samples(dim, n_samples)
    real(DP) :: expected(dim, n_samples)

    mean = [1.0_DP, -2.0_DP, 3.0_DP, -4.0_DP, 5.0_DP, -6.0_DP]
    expected = spread(mean, dim=2, ncopies=n_samples)
    cov = 0.0_DP
    samples = huge(1.0_DP)

    call init_random_seed(.true.)
    call generate_multivariate_normal(mean, cov, samples)

    if (any(samples /= expected)) then
        write(*,*) 'FAIL: zero covariance did not return deterministic mean samples'
        error stop 1
    end if

    cov = 0.0_DP
    cov(1,1) = 1.0_DP
    cov(2,2) = 1.0_DP
    cov(3,3) = 1.0_DP
    cov(4,4) = 1.0_DP
    cov(5,5) = 1.0_DP
    cov(6,6) = 1.0_DP
    samples = expected

    call generate_multivariate_normal(mean, cov, samples)

    if (.not. any(samples /= expected)) then
        write(*,*) 'FAIL: positive-definite covariance did not sample deviations'
        error stop 1
    end if

    write(*,*) 'test_zero_covariance_sampling passed'
end program test_zero_covariance_sampling
```

- [ ] **Step 2: Run the test and verify RED**

Run:

```bash
source setup_env.sh
fpm test test_zero_covariance_sampling
```

Expected: FAIL because Cholesky rejects the zero matrix and the sentinel sample values remain.

- [ ] **Step 3: Implement the zero-covariance fast path**

In `generate_multivariate_normal`, immediately after computing `dim` and `n_particles`, add:

```fortran
        if (all(cov == 0.0_DP)) then
            do j = 1, n_particles
                samples(:, j) = mean
            end do
            return
        end if
```

Do not change the existing Cholesky path.

- [ ] **Step 4: Run the sampler test and verify GREEN**

Run:

```bash
source setup_env.sh
fpm test test_zero_covariance_sampling
```

Expected: PASS and no Cholesky failure message.

- [ ] **Step 5: Commit the sampler behavior**

```bash
git add test/test_zero_covariance_sampling.f90 src/lib/statistics/pod_random_module.f90
git commit -m "feat: support deterministic zero covariance sampling"
```

### Task 2: Add a consistent CLI flag to both propagation programs

**Files:**
- Modify: `test/test_run_srp_uq_propagation_io.f90`
- Create: `test/test_run_hfem_uprop_zero_cov_io.f90`
- Modify: `app/run_srp_uq_propagation.f90:13-225`
- Modify: `app/run_HFEM_uprop.f90:27-235`

**Interfaces:**
- Consumes: CLI token `--zero-init-cov`
- Produces: `cov6 = 0.0_DP` or `cov = 0.0_DP` after OPM loading, plus an explicit log line.

- [ ] **Step 1: Extend the SRP CLI integration test**

In the valid invocation in `test/test_run_srp_uq_propagation_io.f90`, add:

```fortran
                  '--zero-init-cov ' // &
                  '-srp-sigma 0.05 ' // &
```

Replace its existing `-srp-sigma 0.0` argument, then assert:

```fortran
        call assert_file_contains(SUCCESS_LOG, 'Initial orbit covariance: ZERO', &
                                  'zero initial covariance log', n_fail)
        call assert_file_contains(PREFIX // '_moments.json', &
                                  '"eta_sigma_input"', 'eta sigma retained', n_fail)
```

- [ ] **Step 2: Add the HFEM CLI integration test**

Create `test/test_run_hfem_uprop_zero_cov_io.f90` with a short real invocation:

```fortran
program test_run_hfem_uprop_zero_cov_io
    implicit none
    character(len=*), parameter :: prefix = '/tmp/test_hfem_zero_cov'
    character(len=*), parameter :: log_file = '/tmp/test_hfem_zero_cov.log'
    character(len=2048) :: command
    integer :: cmd_status, exit_status, u, ios
    character(len=512) :: line
    logical :: found

    call execute_command_line('rm -f ' // prefix // '_particles.csv ' // &
                              prefix // '_moments.json ' // log_file)

    command = 'fpm run run_HFEM_uprop -- ' // &
              '-cfg config/config.txt ' // &
              '-opm OPM/L1Halo-1/L1Halo-1_init.opm.json ' // &
              '-m DA -dt 60 -o ' // prefix // ' -n 4 -da 2 ' // &
              '--zero-init-cov > ' // log_file // ' 2>&1'
    call execute_command_line(command, cmdstat=cmd_status, exitstat=exit_status)

    if (cmd_status /= 0 .or. exit_status /= 0) error stop 'HFEM zero-cov command failed'

    inquire(file=prefix // '_moments.json', exist=found)
    if (.not. found) error stop 'HFEM zero-cov moments output missing'

    found = .false.
    open(newunit=u, file=log_file, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'HFEM zero-cov log missing'
    do
        read(u, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (index(line, 'Initial orbit covariance: ZERO') > 0) found = .true.
    end do
    close(u)
    if (.not. found) error stop 'HFEM zero-cov log marker missing'

    call execute_command_line('rm -f ' // prefix // '_particles.csv ' // &
                              prefix // '_moments.json ' // log_file)
    write(*,*) 'test_run_hfem_uprop_zero_cov_io passed'
end program test_run_hfem_uprop_zero_cov_io
```

- [ ] **Step 3: Run both CLI tests and verify RED**

Run:

```bash
source setup_env.sh
fpm test test_run_srp_uq_propagation_io
fpm test test_run_hfem_uprop_zero_cov_io
```

Expected: FAIL because the flag is ignored and the required ZERO log marker is absent.

- [ ] **Step 4: Implement the SRP CLI flag**

In `app/run_srp_uq_propagation.f90`:

```fortran
    logical :: has_opm, has_dt, has_et, has_output, zero_init_cov
```

Initialize it:

```fortran
    zero_init_cov = .false.
```

Parse it:

```fortran
        case ('--zero-init-cov')
            zero_init_cov = .true.
```

Immediately after `load_initial_opm`:

```fortran
    if (zero_init_cov) then
        cov6 = 0.0_DP
        write(*,*) 'Initial orbit covariance: ZERO'
    end if
```

Add `--zero-init-cov` to `print_usage`. Do not modify `srp_eta_mean` or `srp_eta_sigma`.

- [ ] **Step 5: Implement the HFEM CLI flag**

In `app/run_HFEM_uprop.f90`, add and initialize `zero_init_cov`, parse the same flag, and after the existing `-pr/-vr` override block add:

```fortran
    if (zero_init_cov) then
        cov = 0.0_DP
        write(*,*) 'Initial orbit covariance: ZERO'
    end if
```

This makes `--zero-init-cov` take precedence if it is combined with `-pr/-vr`. Add the flag to the file's usage comments.

- [ ] **Step 6: Run both CLI tests and verify GREEN**

Run:

```bash
source setup_env.sh
fpm test test_run_srp_uq_propagation_io
fpm test test_run_hfem_uprop_zero_cov_io
```

Expected: both PASS; the SRP moments JSON still reports `eta_sigma_input = 0.05`.

- [ ] **Step 7: Commit both CLI changes**

```bash
git add app/run_srp_uq_propagation.f90 app/run_HFEM_uprop.f90 \
        test/test_run_srp_uq_propagation_io.f90 \
        test/test_run_hfem_uprop_zero_cov_io.f90
git commit -m "feat: add zero initial covariance CLI option"
```

### Task 3: Expose ZERO_INIT_COV in both HALO sweep scripts

**Files:**
- Create: `test/test_zero_init_cov_scripts.sh`
- Modify: `scripts/run_halo_srp_sweep_parallel.sh`
- Modify: `scripts/run_halo_nominal_srp_sweep_parallel.sh`

**Interfaces:**
- Consumes: `ZERO_INIT_COV=1|true|TRUE`
- Produces: `--zero-init-cov` CLI argument, `_zeroInitCov` output suffix, and metadata columns.

- [ ] **Step 1: Write the failing Bash regression test**

Create `test/test_zero_init_cov_scripts.sh`:

```bash
#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

tmp_root="$(mktemp -d)"
trap 'rm -rf "$tmp_root"' EXIT

scripts=(
    scripts/run_halo_srp_sweep_parallel.sh
    scripts/run_halo_nominal_srp_sweep_parallel.sh
)

for script in "${scripts[@]}"; do
    bash -n "$script"
    grep -q -- '--zero-init-cov' "$script"
done

ZERO_INIT_COV=1 DRY_RUN=1 DT_LIST=60 \
OUT_ROOT="$tmp_root/srp" bash scripts/run_halo_srp_sweep_parallel.sh

grep -q '^zero_init_cov=1$' "$tmp_root/srp/batch.log"
grep -q '_zeroInitCov' "$tmp_root/srp/queue.tsv"
grep -q 'zero_init_cov' "$tmp_root/srp/summary.tsv"

ZERO_INIT_COV=TRUE DRY_RUN=1 DT_LIST=60 \
OUT_ROOT="$tmp_root/nominal" bash scripts/run_halo_nominal_srp_sweep_parallel.sh

grep -q '^zero_init_cov=1$' "$tmp_root/nominal/batch.log"
grep -q '_zeroInitCov' "$tmp_root/nominal/queue.tsv"
grep -q 'zero_init_cov' "$tmp_root/nominal/summary.tsv"

echo 'test_zero_init_cov_scripts passed'
```

- [ ] **Step 2: Run the Bash test and verify RED**

Run:

```bash
bash test/test_zero_init_cov_scripts.sh
```

Expected: FAIL because neither script contains `--zero-init-cov`.

- [ ] **Step 3: Implement shared behavior in each script**

Near the other environment defaults, add:

```bash
zero_init_cov_raw="${ZERO_INIT_COV:-0}"
zero_init_cov=0
zero_init_cov_tag=""
zero_init_cov_args=()
if [[ "$zero_init_cov_raw" == "1" || "$zero_init_cov_raw" == "true" || "$zero_init_cov_raw" == "TRUE" ]]; then
    zero_init_cov=1
    zero_init_cov_tag="_zeroInitCov"
    zero_init_cov_args=(--zero-init-cov)
fi
```

Append `$zero_init_cov_tag` to each output base and run-log stem. Add:

```bash
        "${zero_init_cov_args[@]}" \
```

to each `fpm run ... --` argument list.

Add `zero_init_cov` to the summary header and each summary row, and add:

```bash
    printf 'zero_init_cov=%s\n' "$zero_init_cov"
```

to each batch configuration log. Queue paths will inherit the new suffix.

- [ ] **Step 4: Run the Bash test and verify GREEN**

Run:

```bash
bash test/test_zero_init_cov_scripts.sh
```

Expected: PASS for both `1` and `TRUE`.

- [ ] **Step 5: Verify disabled-mode compatibility**

Run:

```bash
tmp_root="$(mktemp -d)"
DRY_RUN=1 DT_LIST=60 OUT_ROOT="$tmp_root/default" \
  bash scripts/run_halo_nominal_srp_sweep_parallel.sh
if grep -q '_zeroInitCov' "$tmp_root/default/queue.tsv"; then exit 1; fi
rm -rf "$tmp_root"
```

Expected: exit 0 and no zero-covariance suffix.

- [ ] **Step 6: Commit the script integration**

```bash
git add scripts/run_halo_srp_sweep_parallel.sh \
        scripts/run_halo_nominal_srp_sweep_parallel.sh \
        test/test_zero_init_cov_scripts.sh
git commit -m "feat: expose zero covariance in halo sweeps"
```

### Task 4: Full verification and documentation check

**Files:**
- Verify: all files from Tasks 1-3
- Verify: `docs/superpowers/specs/2026-08-24-zero-init-cov-design.md`

**Interfaces:**
- Consumes: completed feature
- Produces: fresh evidence that build, targeted tests, and scripts pass.

- [ ] **Step 1: Build both executables**

```bash
source setup_env.sh
fpm build run_HFEM_uprop
fpm build run_srp_uq_propagation
```

Expected: both commands exit 0.

- [ ] **Step 2: Run all feature tests**

```bash
source setup_env.sh
fpm test test_zero_covariance_sampling
fpm test test_run_srp_uq_propagation_io
fpm test test_run_hfem_uprop_zero_cov_io
bash test/test_zero_init_cov_scripts.sh
```

Expected: all four tests PASS.

- [ ] **Step 3: Run the full FPM test suite**

```bash
source setup_env.sh
fpm test
```

Expected: exit 0 with zero failed tests. If unrelated pre-existing failures occur, record their exact names and outputs separately.

- [ ] **Step 4: Check diffs and repository cleanliness**

```bash
git diff --check
git status --short
```

Expected: no whitespace errors; only the user's pre-existing `rename_to_pod.sh` and `setup_env.sh` modifications may remain uncommitted.

- [ ] **Step 5: Confirm requirements line by line**

Verify from fresh command output:

- `ZERO_INIT_COV=1` reaches both CLIs.
- orbital covariance becomes exactly zero;
- SRP sigma remains independently configurable;
- no Cholesky call occurs for the zero matrix;
- output names and metadata distinguish zero-covariance runs;
- default-mode naming and covariance behavior remain unchanged.
