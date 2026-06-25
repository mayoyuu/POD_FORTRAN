# EMDAC GMM Stabilization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stabilize EMDAc-N for the L1Halo-1 optical-only case by first matching the Servadio 2021 algorithm where appropriate, then testing controlled mixture-management changes against the current `.err` and `.gmm_diag` evidence.

**Architecture:** Keep the core EMDAc flow intact: DA/particle time update, K-means/EM GMM refit, Gaussian-sum measurement update, then global mean/covariance recomposition. Add diagnostics and small policy switches around measurement likelihood, weight management, and covariance regularization so each hypothesis can be tested independently.

**Tech Stack:** Fortran/FPM, DACE, OpenBLAS, existing POD_Fortran test programs, `.err`, `.residual`, `.gmm_diag`, and Python-only post-processing scripts for offline analysis.

## Global Constraints

- Do not use a git worktree unless the user is told in advance what will happen in the worktree.
- Do not delete or revert user changes.
- Use TDD for production behavior changes: write a focused failing test, run it, implement, then rerun.
- Run `source setup_env.sh && fpm build` after implementation tasks.
- For orbit-level validation, use `run_pod_emdac` and inspect `.err`, `.residual`, and `.gmm_diag`.
- Baseline configuration for current comparisons: `measurement_noise_ra_arcsec = 1`, `measurement_noise_dec_arcsec = 1`, `process_noise_sigma_a = 1.0e-11`, `use_process_noise = true`.
- Primary screening case: `n=3`, `p=10000`, `o=4`, L1Halo-1 single-station R91 1h data.
- Do not trust final RMS alone; report final/median/max position error, posterior Mahalanobis distance, residual RMS, max component weight, effective component count, min `lambda_min(P)`, max `cond(P)`, and first collapse epoch.

---

## Current Evidence

The new 1 arcsec runs show:

- `n1_p{5,10,50}k` stays stable with final position error around 40 km and posterior Mahalanobis distance around 10.
- `n3/n5` look acceptable through much of the run but collapse in the tail.
- The first serious posterior Mahalanobis excursions occur around 2027-02-13 to 2027-02-15.
- All `n3/n5` cases cross posterior Mahalanobis distance greater than 1000 around 2027-02-24.
- `.gmm_diag` shows that the measurement update often changes balanced prior weights into a single surviving posterior component, for example `0.31,0.37,0.33 -> 0,0,1`.
- Median `lambda_min(P_post) / lambda_min(P_prior)` is approximately 1, so the dominant immediate failure is not per-component covariance shrinkage during the measurement update. The dominant failure is mixture weight collapse.

## Paper Alignment Notes

Servadio 2021 EMDAc-N uses:

- K-means to initialize EM, optionally K-means++ as a more robust alternative.
- EM outputs particle responsibilities `w_i,j` and sums `W_i`.
- Measurement update uses the whole ensemble for each component through `w_i,j / W_i`, not hard cluster membership.
- Measurement likelihood is computed from predicted measurement mean and covariance:
  `P(y | i, Y_k) = N(y; z_i^-, P_ZZ,i)`.
- Weight update is Eq. (34):
  `mu_i^+ = mu_i^- P(y | i, Y_k) / sum_h mu_h^- P(y | h, Y_k)`.
- The posterior global state is recomposed with Eq. (35)-(36).
- The paper does not describe a weight floor, pruning, likelihood tempering, or merge policy. The examples keep a fixed small number of components and use scenarios where the measurement likelihood does not repeatedly annihilate all but one mode.

Current code status:

- `src/functions/orbitimprove/pod_filter_emdac.f90` measurement weight update currently matches Eq. (34) in log-softmax form.
- `src/lib/statistics/pod_gmm_math.f90` EM responsibilities include `log(weights(i))` in `log_P`. This is a standard GMM EM responsibility formula and is consistent with the Bayes form discussed before Algorithm 3, but the extracted Algorithm 3 text is ambiguous because it presents a simplified log-density ratio. Keep this as a verification item rather than changing it blindly.
- Current `sample_particles_from_gmm` can assign zero particles to tiny-weight components. That makes weight collapse irreversible in the next time update.

## Test Scope Answer

Use `p=10000` as the default screening case. It is enough for most implementation iterations because current evidence shows `p=10000` reproduces the failure while running faster than `p=50000`.

Do not rely only on 10k for final conclusions. After a candidate fix passes `n3_p10k`, run:

- `n1_p10k` to ensure the stable single-Gaussian-equivalent case is not degraded.
- `n3_p10k` as the primary failing case.
- `n5_p10k` to verify the fix generalizes to more components.
- One stress check at `n5_p50k` before trusting the conclusion.

## Floor vs Pruning Decision

For this codebase, test floor first, pruning second.

Reason:

- The current failure is premature component death. Pruning intentionally removes low-weight components, so it can make the observed failure irreversible even faster.
- A small posterior weight floor preserves hypotheses long enough for later observations and EM refitting to recover.
- Pruning is useful after there is an adaptive birth/splitting/merge mechanism or if the number of components grows. Here `N` is fixed and small, so pruning without replacement is not the right first defense.

Recommended order:

1. Add diagnostics for measurement likelihood spread.
2. Add optional likelihood tempering.
3. Add optional posterior weight floor.
4. Only then evaluate pruning/merge logic, and only if the number of components becomes adaptive.

## Reference Software

- DACE is the DA toolchain used by the paper; the public DACE repository describes it as the Differential Algebra Computational Toolbox and provides build/use instructions: https://github.com/dacelib/dace
- Stone Soup is a mature Python tracking/state-estimation framework. Its `GaussianMixtureReducer` documents pruning, merging, and truncating for Gaussian mixtures. Use it as a reference for mixture-management ideas, not as a direct EMDAc implementation: https://stonesoup.readthedocs.io/en/latest/stonesoup.mixturereducer.html
- The paper references adaptive Gaussian-sum filtering and Gaussian mixture PHD filtering literature. Those are mature algorithm families for pruning/merging, but they solve a different mixture-management problem than this fixed-N EMDAc implementation.

## Task 1: Lock Down Paper-Formula Regression Tests

**Files:**
- Create: `test/test_gmm_weight_update_formula.f90`
- Modify: `fpm.toml`
- Inspect: `src/functions/orbitimprove/pod_filter_emdac.f90`

**Interfaces:**
- Produces a tiny deterministic test for Eq. (34)-style normalization using known prior weights and likelihoods.
- The test should not require SPICE, DACE propagation, or observation files.

- [ ] Write a pure helper or test-local calculation for:
  `post_i = prior_i * likelihood_i / sum(prior * likelihood)`.
- [ ] Use priors `[0.3, 0.4, 0.3]` and likelihoods `[1.0, 10.0, 2.0]`.
- [ ] Expected posterior:
  `[0.06122448979591837, 0.8163265306122449, 0.12244897959183673]`.
- [ ] Run:
  `source setup_env.sh && fpm test test_gmm_weight_update_formula`
- [ ] Expected before helper extraction: compile failure or missing public helper.
- [ ] Implement the smallest public/testable helper in a statistics module if needed.
- [ ] Rerun:
  `source setup_env.sh && fpm test test_gmm_weight_update_formula`
- [ ] Expected: pass.
- [ ] Rerun:
  `source setup_env.sh && fpm build`
- [ ] Expected: project compiled successfully.

## Task 2: Add Measurement Likelihood Diagnostics

**Files:**
- Modify: `src/functions/orbitimprove/pod_filter_emdac.f90`
- Modify: `src/lib/data/pod_data_format_module.f90`
- Create: `test/test_gmm_likelihood_diagnostics.f90`
- Modify: `fpm.toml`

**Interfaces:**
- Extend `.gmm_diag` rows with measurement-update internals, or add `.gmm_like_diag` if line width becomes unwieldy.
- Required fields: `Z_Mahal_Sq`, `LogLike_NoPrior`, `LogWeightPrior`, `LogWeightPosterior`, `DetPzz`, `Pzz_Cond`, `Innovation_RA_arcsec`, `Innovation_Dec_arcsec`.

- [ ] Write a failing test that calls the diagnostic writer with two components and checks that the header contains `LogLike_NoPrior`, `Z_Mahal_Sq`, and `Pzz_Cond`.
- [ ] Run:
  `source setup_env.sh && fpm test test_gmm_likelihood_diagnostics`
- [ ] Expected: fail because the fields are not written.
- [ ] Store `mahalanobis_sq(i)`, `log_likelihood(i) - log(weight_i)`, `log(weight_i)`, posterior log weight, `det_Pzz(i)`, and a 2x2 condition number in `filter_measurement_update`.
- [ ] Write the diagnostics after the measurement likelihood computation, before local arrays are deallocated.
- [ ] Rerun:
  `source setup_env.sh && fpm test test_gmm_likelihood_diagnostics`
- [ ] Expected: pass.
- [ ] Smoke run `n3_p10k` for only the first 20 observations and confirm the file has finite values.

## Task 3: Verify Angle Innovation Handling

**Files:**
- Modify: `src/lib/measurement` or `src/functions/orbitimprove/pod_filter_emdac.f90`
- Create: `test/test_angle_innovation_wrap.f90`
- Modify: `fpm.toml`

**Interfaces:**
- Produce a helper such as `wrap_angle_rad(delta)` or `compute_optical_innovation(y_meas, z_pred, innovation)`.
- RA innovation must wrap to `[-pi, pi]` before likelihood and Kalman update.

- [ ] Write a failing test:
  measured RA `359.999 deg`, predicted RA `0.001 deg`, expected RA residual about `-0.002 deg`, not `359.998 deg`.
- [ ] Run:
  `source setup_env.sh && fpm test test_angle_innovation_wrap`
- [ ] Expected: fail if current code uses direct subtraction.
- [ ] Use the helper in:
  `innovation(:, i) = y_meas - means_z(:, i)`
  and
  component mean update `y_meas - means_z(:, i)`.
- [ ] Keep Dec innovation unwrapped.
- [ ] Rerun:
  `source setup_env.sh && fpm test test_angle_innovation_wrap`
- [ ] Expected: pass.
- [ ] Rerun:
  `source setup_env.sh && fpm build`
- [ ] Expected: pass.

## Task 4: Test Likelihood Tempering

**Files:**
- Modify: `src/lib/config/pod_config_module.f90` or the existing config module.
- Modify: `config/config.txt`
- Modify: `src/functions/orbitimprove/pod_filter_emdac.f90`
- Create: `test/test_likelihood_tempering.f90`
- Modify: `fpm.toml`

**Interfaces:**
- Config parameter: `gmm_likelihood_temperature`.
- Default value: `1.0`.
- Weight update uses:
  `tempered_log_likelihood = log(prior_weight) + (log_like_no_prior / temperature)`.

- [ ] Write a failing unit test with log likelihoods `[0, -100, -200]`, equal priors, and `T=10`; assert the second component is not rounded to zero.
- [ ] Run:
  `source setup_env.sh && fpm test test_likelihood_tempering`
- [ ] Expected: fail until helper/config exists.
- [ ] Add config parsing with default `1.0`.
- [ ] Apply temperature only to measurement likelihood, not to prior weight.
- [ ] Rerun:
  `source setup_env.sh && fpm test test_likelihood_tempering`
- [ ] Expected: pass.
- [ ] Orbit tests:
  `n3_p10k` with `T=1,2,5,10`.
- [ ] Accept candidate only if `post_max_weight` no longer becomes 1.0 at 2027-02-11/2027-02-15 and final position/MD improve.

## Task 5: Test Posterior Weight Floor

**Files:**
- Modify: config module and `config/config.txt`
- Modify: `src/functions/orbitimprove/pod_filter_emdac.f90`
- Create: `test/test_gmm_weight_floor.f90`
- Modify: `fpm.toml`

**Interfaces:**
- Config parameter: `gmm_weight_floor`.
- Default value: `0.0`.
- After softmax weight update, apply:
  `weight_i = max(weight_i, floor)`, then renormalize.

- [ ] Write a failing test with posterior weights `[1.0, 0.0, 0.0]` and floor `1e-4`.
- [ ] Expected normalized weights:
  approximately `[0.99980004, 0.00009998, 0.00009998]`.
- [ ] Run:
  `source setup_env.sh && fpm test test_gmm_weight_floor`
- [ ] Expected: fail until helper/config exists.
- [ ] Implement floor as a separate helper so it can be tested without full EMDAC.
- [ ] Rerun:
  `source setup_env.sh && fpm test test_gmm_weight_floor`
- [ ] Expected: pass.
- [ ] Orbit tests:
  `n3_p10k` with floor `0, 1e-8, 1e-6, 1e-4`.
- [ ] Accept only if effective component count stays above 1.2 through 2027-02-24 without worsening n1_p10k.

## Task 6: Protect Sampling From Zero-Particle Components

**Files:**
- Modify: `src/functions/orbitimprove/pod_filter_emdac.f90`
- Create: `test/test_gmm_min_particles_per_component.f90`
- Modify: `fpm.toml`

**Interfaces:**
- Config parameter: `gmm_min_particles_per_component`.
- Default value: `0`.
- Allocation must ensure any component with positive retained weight receives at least this many particles, then distribute the remainder by normalized weights.

- [ ] Write a failing test for weights `[0.9998, 1e-4, 1e-4]`, total particles `10000`, min per component `10`.
- [ ] Expected counts all at least 10 and sum exactly 10000.
- [ ] Run:
  `source setup_env.sh && fpm test test_gmm_min_particles_per_component`
- [ ] Expected: fail until allocation helper exists.
- [ ] Implement a helper for counts; call it from `sample_particles_from_gmm`.
- [ ] Rerun:
  `source setup_env.sh && fpm test test_gmm_min_particles_per_component`
- [ ] Expected: pass.
- [ ] Orbit test:
  `n3_p10k` with best floor/temperature from earlier and min particles `10`.
- [ ] Inspect `.gmm_diag` for components returning from tiny weights after later observations.

## Task 7: Covariance Regularization and Inflation

**Files:**
- Modify: covariance helper module or `pod_filter_emdac.f90`
- Create: `test/test_covariance_regularization.f90`
- Modify: `fpm.toml`

**Interfaces:**
- Config parameters:
  `gmm_cov_lambda_min_floor`
  `gmm_cov_inflation`
- Default values:
  `0.0`
  `1.0`

- [ ] Write a failing test with a covariance whose smallest eigenvalue is below `1e-12`.
- [ ] Expected regularized covariance has `lambda_min >= floor`.
- [ ] Run:
  `source setup_env.sh && fpm test test_covariance_regularization`
- [ ] Expected: fail until regularizer exists.
- [ ] Apply regularization after EM covariance construction and after measurement covariance update.
- [ ] Rerun:
  `source setup_env.sh && fpm test test_covariance_regularization`
- [ ] Expected: pass.
- [ ] Orbit test:
  `n3_p10k` with and without regularization, holding floor/temperature fixed.
- [ ] Accept only if Mahalanobis distance improves without increasing residual RMS substantially.

## Task 8: Evaluate Pruning/Merging Only After Stabilization

**Files:**
- Create: `test/test_gmm_prune_merge_policy.f90`
- Modify: GMM utility module if adaptive mixture management is pursued.

**Interfaces:**
- Keep fixed `N` by default.
- If pruning is enabled, do not reduce active components without a replacement/splitting rule.
- Merge only components close in squared Mahalanobis distance and recompute merged covariance with moment matching.

- [ ] Do not implement pruning before Tasks 1-7 are evaluated.
- [ ] If needed, write a test for moment-matched merge of two components.
- [ ] Compare to Stone Soup's documented mixture reducer concepts: pruning removes low-weight components, merging combines nearby components, truncating caps component count.
- [ ] Only enable pruning for adaptive-N experiments, not for the current fixed-N baseline.

## Task 9: Standard Orbit Validation Matrix

**Files:**
- Create or update: `scripts/analyze_emdac_diag.py`
- Create or update: `docs/emdac_run_matrix.md`

**Interfaces:**
- Script reads `.err`, `.residual`, `.gmm_diag`.
- Script prints one row per run with:
  `final_pos`, `median_pos`, `max_pos`, `final_post_md`, `median_post_md`, `max_post_md`, `residual_rms`, `max_weight_median`, `first_weight_collapse_epoch`, `first_post_md_gt_1000`, `min_lmin`, `max_cond`.

- [ ] Write parser against current fixed-width `.err` and whitespace `.gmm_diag`.
- [ ] Include support for Fortran overflow-like numbers such as `1.7976931349+308`.
- [ ] Run on current 9 cases and save a baseline table.
- [ ] For every candidate fix, run:
  `n1_p10k`, `n3_p10k`, `n5_p10k`.
- [ ] For final acceptance, run:
  `n1_p5k`, `n1_p10k`, `n1_p50k`,
  `n3_p5k`, `n3_p10k`, `n3_p50k`,
  `n5_p5k`, `n5_p10k`, `n5_p50k`.

## Acceptance Criteria

A candidate change is useful only if all of the following hold:

- `n1_p10k` remains stable and final position error is not materially worse than the current approximately 40 km baseline.
- `n3_p10k` no longer collapses to a single posterior weight before 2027-02-24.
- `n3_p10k` final position error improves by at least 10x from the current approximately 153,525 km.
- `n3_p10k` posterior Mahalanobis distance no longer exceeds 1000 around 2027-02-24.
- `n5_p10k` shows the same qualitative improvement.
- `.gmm_diag` effective component count remains above 1.2 through the critical 2027-02-13 to 2027-02-24 interval.

## Recommended Execution Order

1. Task 1: formula regression tests.
2. Task 2: likelihood diagnostics.
3. Task 3: angle innovation wrap.
4. Task 4: likelihood tempering.
5. Task 5: posterior weight floor.
6. Task 6: minimum particles per retained component.
7. Task 7: covariance regularization/inflation.
8. Task 9: validation matrix script.
9. Task 8: pruning/merging only if adaptive-N becomes necessary.

## Notes For Interpreting Future Runs

- If posterior weights collapse but likelihood diagnostics show one component genuinely has much smaller innovation, tempering/floor is a robustness choice, not a correction of the paper formula.
- If all components have similar innovation but weights still collapse, inspect determinant and `Pzz` conditioning; the likelihood may be dominated by an over-small `Pzz`.
- If `Pzz` determinant is tiny and `Z_Mahal_Sq` is huge, inspect measurement units and angle wrapping first.
- If weight floor preserves bad branches that later pollute the global mean, reduce the floor or pair it with likelihood tempering instead of raising the floor further.
- Pruning is not a fix for current fixed-N collapse. It is a mixture-size control tool.
