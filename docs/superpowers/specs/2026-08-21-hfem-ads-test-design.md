# HFEM ADS propagation test design

## Goal

Add a manually runnable Fortran test that propagates the existing `L1Halo-1`
case for seven days with the HFEM high-fidelity force model and ADS. The test
runs two domains:

1. Six orbital uncertainty variables.
2. Six orbital uncertainty variables plus the SRP scale parameter `eta_srp`.

The implementation is written and registered with fpm, but is not compiled or
run by Codex, as requested by the user.

## Default case

- Configuration: `config/config.txt`.
- Initial state: `OPM/L1Halo-1/L1Halo-1_init.opm.json`.
- Propagation duration: 604800 seconds (seven days).
- Integrator: HFEM DA RKF78.
- DA order: 4.
- Orbital 1-sigma uncertainty: 100 km in each position component and 0.3 m/s
  in each velocity component.
- ADS orbital domain: plus or minus 3 sigma.
- SRP parameter: `eta_srp`, with nominal value zero and bounded ADS half-width
  0.1. The physical SRP multiplier therefore spans `[0.9, 1.1]`.
- Maximum splitting-history depth: 8 per accepted Patch.

## Normalized-domain mapping

ADS always splits normalized independent variables
`xi(i) in [-1, 1]`. Physical scaling is embedded in the DA maps before each
propagation.

For the orbital variables:

```
r_da(i) = (r_mean(i) + r_half_width(i) * xi(i)) / LU
v_da(i) = (v_mean(i) + v_half_width(i) * xi(i+3)) / VU
```

For the seventh SRP variable:

```
eta_srp = eta_mean + eta_half_width * xi(7)
a_srp   = a_srp_nominal * (1 + eta_srp)
```

Particle or deterministic point evaluation performs the inverse normalization
first, selects the containing patch in the global unit domain, maps that point
to the patch-local unit domain, and evaluates the accepted DA flow map.

## Critical parameter-splitting rule

The six orbital DA variables are embedded in the propagated state, so
`patch_split` automatically composes their physical scaling with each child
domain transformation. `eta_srp` is different: it is created by the HFEM force
model and is not a seventh integrated state component.

Consequently, every HFEM patch integration must reconstruct the local SRP
parameter map from the patch splitting history. If the seventh normalized
coordinate of a patch has center `c7` and width `w7`, the force model must use:

```
eta_patch_center = eta_mean + eta_half_width * c7
eta_patch_span   = eta_half_width * w7 / 2
eta_srp          = eta_patch_center + eta_patch_span * xi_local(7)
```

For the root patch, `c7 = 0` and `w7 = 2`, which recovers the full physical
interval. For a left/right child, the center and span become the corresponding
physical half-interval. Merely recording a split in direction 7 while continuing
to use `eta_mean + eta_half_width * xi(7)` for both children is incorrect.

## ADS core generalization

The existing Fortran ADS core conflates two dimensions:

- Flow-map component count: six propagated Cartesian state components.
- Independent-variable count: currently hard-coded to six, but seven for the
  SRP case.

Restore the same separation used by the C++ reference. The active DACE engine
is the authority for the independent-variable count, exposed through
`dace_max_variables()`. Center, width, containment, point mapping, replay, and
split-direction selection query that value. Patch construction, translation,
and truncation-error estimation use `patch%da_vec%size` for the output-component
count. No dimension is stored redundantly in the splitting history. Existing
CRTBP callers remain compatible because they initialize DACE with six
independent variables and propagate six output components.

## SRP force-model interface

Extend `set_srp_scale_uncertainty` with an optional DA span argument. Its
default is one, preserving the existing DA-MC convention in which the DA
evaluation coordinate is the physical `eta_srp` deviation.

For HFEM ADS, the driver supplies the patch-local center and span described
above before every patch integration. `clear_srp_scale_uncertainty` resets the
index, nominal value, and span after each scenario and on normal completion.

## HFEM ADS test driver

Add `test/test_hfem_ads_propagation.f90` with a local HFEM-specific BFS driver.
The driver reuses the production ADS types and split operations but calls
`da_adaptive_step_integrate` directly for HFEM propagation.

For every queued patch it:

1. Reconstructs the patch-local SRP map when running the 7D case.
2. Integrates the six-component nondimensional HFEM DA state with RKF78.
3. Converts the final DA state to km and km/s before truncation-error testing.
4. Accepts the flow map if all component tolerances pass or the split-depth
   limit is reached.
5. Otherwise chooses the largest offending output component, selects the
   independent variable with the largest top-order contribution, and splits
   the input patch.

The same executable runs the 6D case first and the 7D SRP case second.

## Test checks and output

The program reports for each scenario:

- DA order and normalized/physical domain half-widths.
- Number of accepted patches, BFS iterations, and maximum queue size.
- Split counts for every independent variable, including `eta_srp`.
- Wall time for domain construction and deterministic point evaluation.
- Final states evaluated at the domain center and selected boundary/interior
  points.

The program fails if:

- No patch is accepted.
- A deterministic normalized point cannot be assigned to a patch.
- Any evaluated state component is non-finite.
- Zero-duration mapping does not recover the requested physical initial state.
- In the 7D case, `eta_srp = -0.1` and `eta_srp = +0.1` produce no resolvable
  difference after propagation.

Exact patch counts and timings are diagnostic rather than assertions because
they depend strongly on force-model and integrator tolerances.

## Files

Modify:

- `src/lib/math/pod_ads_split_module.f90`
- `src/lib/system/pod_dace_classes.f90`
- `src/lib/forcemodel/pod_da_force_model_module.f90`
- `fpm.toml`

Add:

- `test/test_hfem_ads_propagation.f90`

No existing user-owned setup-script changes are touched.

## Manual verification

Codex does not execute these commands. The user can build and run the new test
after reviewing the changes:

```
fpm build
fpm test test_hfem_ads_propagation
```

The existing ADS domain-scale and CRTBP tests should also be run to check that
the six-variable default remains backward compatible.
