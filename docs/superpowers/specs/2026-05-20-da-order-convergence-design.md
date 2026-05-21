# DA Order Convergence Test Design

## Goal

Determine the minimum DA (Differential Algebra) expansion order needed for 15-day
orbit uncertainty propagation from the L1 Halo OPM, such that every particle's
DA-predicted state matches the MC-propagated state within 1% relative error in
all six components.

## Parameters

| Parameter | Value |
|-----------|-------|
| OPM | `OPM/L1Halo/L1Halo_init.opm.json` |
| Particles | 1000 |
| Duration | 15 days (1,296,000 s) |
| DA orders tested | 3, 4, 5, 6 |
| Integrator | RKF45 |
| Force model | Earth 10×10, Moon 10×10, Sun third-body, SRP |

## Error metric

For each direction component \(c \in \{x, y, z, v_x, v_y, v_z\}\) and each particle \(i\):

\[
\text{rel\_err}(c, i) = \frac{|\text{DA}_c(i) - \text{MC}_c(i)|}{|\text{MC}_c(i)|}
\]

An order **passes** if \(\max_{c,i}\ \text{rel\_err}(c, i) < 0.01\).

## Architecture

New file: `test/test_da_convergence.f90`

### Phase 0 — System init
- Load config, init SPICE kernels, init gravity network
- Init DACE engine once with `dace_initialize(max_order=6, max_vars=6)` — the max
  order we'll test; individual DA runs use `dace_push_to(order)` internally

### Phase 1 — Read OPM
- Parse `OPM/L1Halo/L1Halo_init.opm.json`
- Extract: epoch (convert UTC → TDB), nominal state (6), covariance (6×6)

### Phase 2 — Generate particles
- Call `run_uq_propagation` with DA method just to generate initial particles
  into a `uq_state_type` object
- This object is reused for both MC and DA phases (same random seeds)

### Phase 3 — MC golden reference
- Call `run_particle_propagation` with `method_switch=METHOD_MC`, `integrator_switch=RKF45`
- Propagate all 1000 particles through full 15-day arc
- Save `final_mc%samples(6, 1000)` as golden reference

### Phase 4 — DA scan
- For each order in [3, 4, 5, 6]:
  - Call `run_particle_propagation` with `method_switch=METHOD_DA`, `da_order=order`, `integrator_switch=RKF45`
  - The DA propagator's internal `dace_push_to(order)` handles order switching
  - Compute per-component max relative error vs MC golden
  - Print per-order results immediately

### Phase 5 — Summary
- Table: order | max_rel_err per component | pass/fail
- Conclusion: minimum order that satisfies the 1% threshold

## Output example

```
============================================================
  DA Order Convergence Test — L1Halo / 15-day / RKF45
============================================================
Particles:    1000
DA orders:    3, 4, 5, 6

Order  |   X(%)  |   Y(%)  |   Z(%)  |  VX(%)  |  VY(%)  |  VZ(%)  | Status
  3    |  x.xx   |  x.xx   |  x.xx   |  x.xx   |  x.xx   |  x.xx   | FAIL/PASS
  4    |  x.xx   |  x.xx   |  x.xx   |  x.xx   |  x.xx   |  x.xx   | FAIL/PASS
  5    |  x.xx   |  x.xx   |  x.xx   |  x.xx   |  x.xx   |  x.xx   | FAIL/PASS
  6    |  x.xx   |  x.xx   |  x.xx   |  x.xx   |  x.xx   |  x.xx   | FAIL/PASS

CONCLUSION: Minimum DA order = X
```

## Key implementation notes

- Use `run_particle_propagation` (not `run_uq_propagation`) for phases 3-4 to
  share pre-generated initial particles and bypass the built-in random sampling
- DACE initialized once at max order 6; `dace_push_to(order)` inside
  `da_propagate` handles per-run order switching
- MC saves results to memory (not file) for comparison — use a local copy of
  `final_mc%samples`
- RKF45 is passed via `integrator_switch` parameter in both DA and MC calls
- Verbosity should be enabled for MC, disabled for DA during the scan loop
