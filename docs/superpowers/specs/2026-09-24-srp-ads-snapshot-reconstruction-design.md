# SRP ADS Snapshot Reconstruction Design

## Objective

Add a standalone offline program that reads one saved ADS patch snapshot and
densely samples the union of its DA polynomial patches. It must reconstruct
both 7D SRP-uncertain snapshots and padded 6D eta-zero control snapshots
without rerunning orbit propagation.

## Approaches considered

1. **Patch-local deterministic resampling — selected.** Read each patch's
   local polynomial and sample its local unit box. This preserves every saved
   patch, scales with patch count, and needs no DACE runtime.
2. **One dense global design evaluated against all patches.** This is useful
   for strict treatment/control pairing, but inefficient for small patches and
   can undersample heavily split regions.
3. **Boundary optimization or explicit surface meshing.** This could estimate
   projection boundaries more sharply, but it is substantially more complex
   and cannot represent the complete 6D image as one surface.

The selected utility focuses on inspecting one reconstructed ADS domain.
Strict high-density treatment/control differences remain a separate future
operation that must use one common global point set.

## Inputs and interface

The command is:

```sh
python3 tools/reconstruct_srp_ads_snapshot.py \
  --prefix PATH/ORBIT_patch_s0003600 \
  --points-per-patch 4096 \
  --output PATH/reconstructed.csv
```

The prefix identifies:

- `PREFIX_meta.csv`
- `PREFIX_coeff.csv`
- `PREFIX_splits.csv`
- `PREFIX_map.csv`

The program infers the parent run prefix by removing
`_patch_s<seven digits>`. It then reads `RUN_history.csv`,
`RUN_initial_points.csv`, and `RUN_shape.csv` when available. An explicit
`--history` path overrides inferred history. `--absolute-only` allows
reconstruction when no history exists. `--no-verify` skips comparison with
the saved 512-point cloud.

Defaults are 4096 points per patch and a boundary fraction of 0.5. The program
uses only Python's standard library.

## Snapshot model

For active dimension `d`, each patch metadata row defines global center
`c` and full width `w`. Local coordinates `v in [-1,1]^d` map to global
normalized coordinates as:

```text
u_i = c_i + 0.5 * w_i * v_i
```

A global point maps back to one patch by containment
`abs(u_i-c_i) <= 0.5*w_i`, then:

```text
v_i = 2 * (u_i-c_i) / w_i
```

Each sparse coefficient row contributes:

```text
coefficient * product(v_i ** power_i)
```

to one of six Cartesian state components. Components 1–3 are km and 4–6 are
km/s. A 6D control snapshot has `center7=width7=power7=0`; reconstructed
`local7`, `u7`, and physical eta are zero.

## Sampling design

Every patch receives exactly `points_per_patch` deterministic points:

1. all `2^d` local corners;
2. all `2d` local axis endpoints;
3. the local center;
4. additional face Halton points assigned cyclically across the `2d` faces;
5. additional interior Halton points.

The requested boundary fraction controls the fraction of additional points
placed on faces. Halton bases are the first seven primes. The program rejects
a requested count smaller than `2^d + 2d + 1`.

These points reconstruct the saved ADS polynomial image at controllable
density. They do not turn the sampled image into a rigorous continuous-domain
or physical-truth guarantee.

## Outputs

The reconstructed CSV contains:

- `patch_id`, per-patch `sample_id`, and `sample_kind`;
- `local1...local7`;
- global `u1...u7`;
- physical relative SRP `eta`;
- absolute `x_km,y_km,z_km,vx_kms,vy_kms,vz_kms`;
- when history is available,
  `dx_km,dy_km,dz_km,dvx_kms,dvy_kms,dvz_kms`.

The JSON report contains schema version, source paths, hashes, active
dimension, DA order, snapshot time, patch/sample counts, sampling settings,
state/deviation extrema, structural validation, and saved-cloud verification
errors.

## Structural validation and failure handling

Before sampling, the program rejects:

- missing or duplicate patch IDs;
- inconsistent snapshot times;
- missing state components;
- invalid powers, components, widths, or non-finite values;
- nonzero seventh powers in a 6D snapshot;
- split histories that do not reproduce metadata center and width;
- coefficients referring to unknown patches.

Errors identify the file, patch, and offending field.

## Saved-cloud verification

When the parent initial-points, history, and shape files exist, the program:

1. finds the history row at `time_seconds/3600`;
2. evaluates each of the 512 saved global unit points with the first containing
   patch, matching the Fortran manifold rule;
3. compares absolute reconstructed state against
   `nominal + saved shape deviation`;
4. records maximum position, velocity, and component differences;
5. fails if position exceeds `1e-8 km` or velocity exceeds
   `1e-11 km/s`.

This verifies coefficient parsing, patch selection, coordinate mapping, and
polynomial evaluation against the original runner outputs.

## Testing

Synthetic tests cover:

- one unsplit 7D polynomial with known values;
- two split patches and signed split-history reconstruction;
- padded 6D control coefficients and zero eta;
- deterministic sample counts and boundary/interior classifications;
- saved-cloud verification;
- rejection of malformed metadata, powers, and split history.

A real runner snapshot is reconstructed in a short integration check and must
match its saved 512-point cloud within the stated tolerances.
