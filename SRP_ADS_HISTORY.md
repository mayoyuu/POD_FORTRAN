# SRP–ADS uncertainty-domain history

## Run one orbit

From the project root, after building with the existing FPM environment:

```sh
fpm run --target run_srp_ads_history -- -opm input/DROb_20251210_9.opm -o output/DROb_srp_ads --eta-sigma 0.02
```

The runner accepts the repository's JSON .opm or .opm.json files. The OPM must contain a complete 6×6 Cartesian covariance. Defaults: 15 days, one checkpoint per hour (hours 0–360), DA order 4, maximum split depth 8, position tolerance 0.1 km, velocity tolerance 1e-6 km/s. Optional arguments are --days, --hours, --save-hours, --da-order, --max-depth, --pos-tol-km, --vel-tol-kms, and -cfg.

The initial domain is the whitened 6D ±3σ covariance box plus a seventh ±3σ SRP relative multiplier η. Its mean is zero; --eta-sigma supplies its standard deviation. The OPM state follows the project's GCRF/J2000 convention. Positions and velocities are in km and km/s.

## Outputs

For output prefix PREFIX:

- PREFIX_history.csv: checkpoint epoch, nominal real-SRP state, **sampled** component minima/maxima relative to nominal, maximum sampled position/velocity norms, three position and velocity principal **1σ semiaxis** lengths from the 512-point cloud, patch count, depth-limited patch count, split counts by variable, queue diagnostics, and validation discrepancies.
- PREFIX_shape.csv: 512 fixed unit-box points per checkpoint, with stable IDs and six Cartesian errors relative to nominal.
- PREFIX_initial_points.csv: each point's seven unit-box coordinates, physical initial state, η, and validation membership. Point 143 is the domain center.
- PREFIX_report.json: first saved hour failing the independent real-SRP comparison, previous passing hour, worst point/error, first hour where splitting reached maximum depth, thresholds, and validation scope.
- PREFIX_patch_sXXXXXXX_{meta,coeff,splits,map}.csv: sparse DA polynomials and split history at the start, end, and first failed or first unconfirmed checkpoint with its predecessor. Suffix is elapsed seconds. Each coefficient uses seven **local patch coordinates**; splits.csv and meta.csv locate the patch in the original unit box. map.csv stores the whitening basis and units.

The 207 validation probes are 128 seven-dimensional box corners, 14 axial points, the center, and 64 fixed interior Halton points. They are propagated by an independent real-valued SRP integrator using the same initial states and η values. A checkpoint fails when any probe differs from ADS by **more than 0.1 km or 1e-6 km/s** by default. A depth-limited patch is reported separately as an unconfirmed ADS truncation tolerance. Processing continues after first failure. If no failure occurs, the report states only that none was observed at the saved times and sampled probes. Component ranges and principal axes are sampled estimates, not rigorous interval bounds.

To make one SVG showing selected saved hours and the size history:

```sh
python3 tools/plot_srp_ads_history.py --prefix output/DROb_srp_ads --hours 0 24 360 --output output/DROb_srp_ads_domain.svg
```

The script uses only Python's standard library and renders selected checkpoints on demand.
