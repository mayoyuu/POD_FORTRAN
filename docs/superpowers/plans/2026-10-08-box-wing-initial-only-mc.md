# Box-Wing Initial-Only ADS–MC Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Earth-pointing and Moon-pointing box-wing 6D ADS propagation that reuses the current 50,000 fixed initial-error samples and saves 2,048 points hourly plus 50,000 points at fixed and diagnostic keyframes.

**Architecture:** Add two immutable config files and a standalone single-group Fortran runner so the current paired cannonball runner and its running binary remain unchanged. Extend the existing MC I/O module with a single-group point schema, then add a Python campaign wrapper that stages the new executable, packs point CSV files, audits each orbit, and atomically publishes Earth/Moon outputs.

**Tech Stack:** Fortran 2008, FPM, DACE/ADS, Python 3, NumPy, h5py, unittest.

## Global Constraints

- Perform implementation in an isolated Git worktree and use its own FPM build directory; do not rebuild the main checkout while the existing campaign is running.
- Do not modify `config/config.txt`.
- Use OPM full 6×6 covariance and exactly the existing sample file's `sample_id,u1..u6`.
- Keep SRP scale error, attitude bias, and array-angle error equal to zero.
- Defaults: 50,000 samples, 2,048 hourly saved points, 360 h, 1 h cadence, DA order 4, maximum ADS depth 8, 0.1 km and 1e-6 km/s validation thresholds.
- Evaluate all 50,000 points every hour; save all points only at 0, 24, 72, 168, 240, 360 h and first split/depth-limit/validation-failure event pairs.
- Preserve model-induced center displacement by saving nominal centers and absolute sample states.
- A failed orbit must not terminate the remaining batch.

---

### Task 1: Add and verify the two box-wing configurations

**Files:**
- Create: `config/config_box_wing_earth.txt`
- Create: `config/config_box_wing_moon.txt`
- Create: `test/test_box_wing_campaign_configs.f90`

**Interfaces:**
- Consumes: `pod_config::set_default_config,load_config,validate_config`.
- Produces: two complete config paths accepted by `pod_engine_init`.

- [ ] **Step 1: Write the failing configuration test**

Create a Fortran test that loads each file after `set_default_config()` and checks:

```fortran
call load_and_check('config/config_box_wing_earth.txt','earth')
call load_and_check('config/config_box_wing_moon.txt','moon')
...
call assert_true(trim(config%srp_model)=='box_wing','box-wing model')
call assert_true(trim(config%srp_attitude_mode)==expected_mode,'attitude')
call assert_close(config%srp_box_dimensions_m,[2.0_DP,3.0_DP,4.0_DP],'box')
call assert_close(config%srp_box_optical,[0.30_DP,0.40_DP,0.30_DP],'box optical')
call assert_scalar(config%srp_array_total_area_m2,24.0_DP,'array area')
call assert_close(config%srp_array_front_optical,[0.85_DP,0.08_DP,0.07_DP],'front')
call assert_close(config%srp_array_back_optical,[0.60_DP,0.20_DP,0.20_DP],'back')
call assert_scalar(config%srp_scale_da_span,0.0_DP,'scale span')
call assert_close(config%srp_attitude_bias_span_arcsec,[0.0_DP,0.0_DP,0.0_DP],'attitude span')
call assert_scalar(config%srp_array_angle_span_deg,0.0_DP,'array span')
call assert_true(validate_config(),'valid campaign config')
```

- [ ] **Step 2: Run the test and verify it fails because the files do not exist**

Run: `fpm test --target test_box_wing_campaign_configs`
Expected: nonzero exit while loading `config_box_wing_earth.txt`.

- [ ] **Step 3: Create complete configs from `config/config.txt`**

Copy every existing non-SRP setting verbatim. Replace the SRP block with the documented 1200 kg, 2×3×4 m, 24 m² and optical values. Use `srp_pressure_1au_n_m2 = 4.559820616707957e-6`, zero all uncertainty spans, and vary only `srp_attitude_mode`.

- [ ] **Step 4: Run the configuration tests**

Run: `fpm test --target test_box_wing_campaign_configs`
Expected: PASS for both modes.

- [ ] **Step 5: Commit**

```bash
git add config/config_box_wing_earth.txt config/config_box_wing_moon.txt test/test_box_wing_campaign_configs.f90
git commit -m "feat: add box-wing campaign configs"
```

### Task 2: Add single-group point CSV I/O

**Files:**
- Modify: `src/lib/uncertainty/io/pod_srp_mc_io_module.f90`
- Modify: `test/test_srp_mc_io.f90`

**Interfaces:**
- Produces:
  - `write_initial_only_point_header(unit_value,status,message)`
  - `write_initial_only_point_rows(unit_value,hour,keyframe_kind,sample_ids,states,center,status,message)`
  - schema constant `SRP_MC_INITIAL_POINTS_SCHEMA='srp_mc_initial_points_v1'`.

- [ ] **Step 1: Extend the I/O test before implementation**

Add two finite states, a six-component center and IDs 11/12. Assert the emitted header contains:

```text
hour,keyframe_kind,sample_id,x_km,y_km,z_km,vx_kms,vy_kms,vz_kms,dx_km,dy_km,dz_km,dvx_kms,dvy_kms,dvz_kms
```

Assert the first row contains the absolute state and `state-center`. Add failure cases for mismatched state shape, mismatched ID count, nonfinite state and nonfinite center.

- [ ] **Step 2: Run the focused test and verify missing symbols fail compilation**

Run: `fpm test --target test_srp_mc_io`
Expected: compile failure naming the new procedures.

- [ ] **Step 3: Implement the schema and procedures**

Validate `size(states,1)==6`, `size(states,2)==size(sample_ids)`, finite hour/states/center and a comma-free nonempty keyframe kind. Write each row with:

```fortran
error_state=states(:,row)-center
write(unit_value,'(ES25.16E3,",",A,",",I0,12(",",ES25.16E3))') &
    hour,trim(keyframe_kind),sample_ids(row),states(:,row),error_state
```

- [ ] **Step 4: Run the focused test**

Run: `fpm test --target test_srp_mc_io`
Expected: `PASS: SRP MC sample and point I/O.`

- [ ] **Step 5: Commit**

```bash
git add src/lib/uncertainty/io/pod_srp_mc_io_module.f90 test/test_srp_mc_io.f90
git commit -m "feat: write single-group SRP MC points"
```

### Task 3: Implement the 6D box-wing ADS–MC runner

**Files:**
- Create: `app/run_srp_ads_mc_initial_only.f90`
- Create: `test/test_srp_ads_mc_initial_runner.py`

**Interfaces:**
- CLI: `run_srp_ads_mc_initial_only -opm FILE -o PREFIX --mc-samples FILE -cfg FILE [options]`.
- Outputs: `_mc_hourly_stats.csv`, `_mc_visualization_points.csv`, `_mc_keyframe_points.csv`, `_split_events.csv`, `_mc_report.json`, and snapshot files.

- [ ] **Step 1: Add CLI contract tests**

Use a temporary directory and a mock executable invocation where possible. Check `--help`, missing required arguments, rejection of cannonball config, rejection of a positive SRP uncertainty option, and JSON report fields:

```python
self.assertEqual(report["output_schema_version"], "srp_ads_mc_initial_only_v1")
self.assertEqual(report["active_dimensions"], 6)
self.assertEqual(report["srp_eta"], 0.0)
self.assertIn(report["srp_attitude_mode"], ("earth", "moon"))
```

- [ ] **Step 2: Run the new test and verify the executable is absent**

Run: `python3 -m unittest test/test_srp_ads_mc_initial_runner.py -v`
Expected: FAIL because the application is not built.

- [ ] **Step 3: Implement argument parsing and model validation**

Mirror the paired runner options except `--eta-sigma`. After `pod_engine_init(config_file)`, require:

```fortran
if(trim(config%srp_model)/='box_wing') call fatal('requires srp_model=box_wing')
if(trim(config%srp_attitude_mode)/='earth'.and. &
   trim(config%srp_attitude_mode)/='moon') call fatal('requires earth or moon attitude')
if(config%srp_scale_da_span/=0.0_DP) call fatal('SRP scale uncertainty must be zero')
```

- [ ] **Step 4: Initialize one 6D ADS history and matching Real probes**

Set `coordinate_mode=ADS_COORD_WHITENED`, `srp_sigma=0`, `da_storage_variables=6`, and the existing error tolerances. Require both `coordinate_map%n_variables==6` and `history%da_variables==6`. Evaluate MC inputs with `u(1:6,:)`; map all validation probes with `probe_points(1:6,:)`; set every Real-probe eta to zero.

- [ ] **Step 5: Implement the hourly loop and statistics**

Advance one history and one Real validation set. Evaluate 50,000 points, the zero point and 512 shape points. Reuse `compute_group_stats`, `compute_shape_summary`, `assess_accuracy`, the existing split-event CSV fields and failure tracking. Write one row for each 10k/20k/50k prefix with `group=<attitude mode>`.

- [ ] **Step 6: Implement point and snapshot saving**

Use `choose_visualization_ids` unchanged so the same 2,048 IDs are selected. Write absolute states and model-relative errors through Task 2. Fixed keyframes are 0/24/72/168/240/360 h; diagnostic keyframes include first split, first depth-limited hour and first validation failure, plus their prior hour. Save snapshots through `write_ads_patch_snapshot`.

- [ ] **Step 7: Implement the report**

Include model, attitude, active/storage dimensions, OPM/config/sample paths, counts, tolerances, first split criterion, first unconfirmed hour, first failure and previous pass, maximum validation errors, elapsed time and keyframe hours. Explicitly state that point statistics are fixed truncated-normal design metrics.

- [ ] **Step 8: Build and run contract tests**

Run: `fpm build --target run_srp_ads_mc_initial_only`
Run: `python3 -m unittest test/test_srp_ads_mc_initial_runner.py -v`
Expected: all tests PASS.

- [ ] **Step 9: Commit**

```bash
git add app/run_srp_ads_mc_initial_only.f90 test/test_srp_ads_mc_initial_runner.py
git commit -m "feat: propagate initial-only box-wing ADS clouds"
```

### Task 4: Pack and audit single-group point outputs

**Files:**
- Create: `tools/pack_srp_mc_initial_points_hdf5.py`
- Create: `test/test_pack_srp_mc_initial_points_hdf5.py`

**Interfaces:**
- `pack_initial_points(csv_path,h5_path,expected_points,kind,report_path) -> dict`.
- HDF5 datasets: `hour`, `keyframe_kind`, `sample_id`, `state`, `error`.

- [ ] **Step 1: Write synthetic CSV packing tests**

Build two hours with stable IDs. Assert HDF5 shapes, units, exact absolute state/error values, unique IDs, hour count and schema attributes. Add rejection tests for missing IDs, duplicate `(hour,sample_id)`, nonfinite values and wrong row count.

- [ ] **Step 2: Run the test and verify import failure**

Run: `python3 -m unittest test/test_pack_srp_mc_initial_points_hdf5.py -v`
Expected: FAIL because the module does not exist.

- [ ] **Step 3: Implement streaming packing**

Read CSV in chunks, validate the schema line and header, create compressed chunked float64 datasets, and atomically replace the final HDF5 only after validation. Store source/report SHA-256, model, attitude, units and schema version as attributes.

- [ ] **Step 4: Run tests**

Run: `python3 -m unittest test/test_pack_srp_mc_initial_points_hdf5.py -v`
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```bash
git add tools/pack_srp_mc_initial_points_hdf5.py test/test_pack_srp_mc_initial_points_hdf5.py
git commit -m "feat: pack initial-only SRP point clouds"
```

### Task 5: Add the Earth/Moon batch campaign wrapper

**Files:**
- Create: `tools/run_box_wing_initial_mc_campaign.py`
- Create: `test/test_run_box_wing_initial_mc_campaign.py`
- Create: `20261005SRP/RUN_BOX_WING_INITIAL_CAMPAIGN.md`

**Interfaces:**
- CLI supports `--attitudes earth moon`, `--jobs`, `--dry-run`, `--audit-only`, `--force`, `--executable`, and all propagation defaults.
- Results root: `20261005SRP/results_box_wing_initial_only/<attitude>/<orbit>`.

- [ ] **Step 1: Write batch tests with a fake propagator**

Verify command construction omits `--eta-sigma`, chooses the correct config, retains the sample path, runs both attitudes, keeps failures isolated, rejects incomplete outputs, and atomically replaces only validated results.

- [ ] **Step 2: Run tests and verify module absence**

Run: `python3 -m unittest test/test_run_box_wing_initial_mc_campaign.py -v`
Expected: FAIL because the script does not exist.

- [ ] **Step 3: Implement preflight and execution**

Reuse the existing 15-orbit whitelist and OPM covariance audit by importing side-effect-free helpers from `run_srp_mc_campaign.py`. Stage the executable by SHA-256. For each attitude/orbit, create a UUID staging directory, run the new executable with the common sample file, inject OPM/config/sample/executable hashes, pack both point files, audit 361 hours and expected point counts, then publish atomically.

- [ ] **Step 4: Implement resumability and batch reports**

Skip complete directories unless `--force`. Write `audit/input_audit.csv`, `audit/campaign_results.csv`, per-orbit manifests and `population/batch_summary.json`. Return nonzero only after all submitted jobs finish if any orbit failed.

- [ ] **Step 5: Write the operation guide**

Document isolated build/staged executable use, dry run, one-orbit smoke invocation, full Earth/Moon launch, background launch, process inspection and stop commands. State that the current cannonball campaign must not be rebuilt or replaced.

- [ ] **Step 6: Run batch tests**

Run: `python3 -m unittest test/test_run_box_wing_initial_mc_campaign.py -v`
Expected: all tests PASS.

- [ ] **Step 7: Commit**

```bash
git add tools/run_box_wing_initial_mc_campaign.py test/test_run_box_wing_initial_mc_campaign.py 20261005SRP/RUN_BOX_WING_INITIAL_CAMPAIGN.md
git commit -m "feat: batch box-wing initial-only propagation"
```

### Task 6: Verify physics, output contracts and a real short run

**Files:**
- Modify only if a concrete failure identifies a root cause in files from Tasks 1–5.

**Interfaces:**
- Produces a staged executable and short-run artifacts outside tracked source directories.

- [ ] **Step 1: Run focused Fortran tests**

```bash
fpm test --target test_box_wing_campaign_configs
fpm test --target test_srp_mc_io
fpm test --target test_halo_srp_real_da_10day
```

Expected: all PASS.

- [ ] **Step 2: Run all new Python tests**

```bash
python3 -m unittest   test/test_srp_ads_mc_initial_runner.py   test/test_pack_srp_mc_initial_points_hdf5.py   test/test_run_box_wing_initial_mc_campaign.py -v
```

Expected: all PASS.

- [ ] **Step 3: Run one-hour Earth and Moon smoke propagations**

Use one eligible OPM, 50,000 MC points, 2,048 visualization points and `--hours 1`. Expected for each mode: two unique hours, finite 10k/20k/50k statistics, 4,096 visualization rows, 100,000 fixed start/end keyframe rows, 6 active dimensions, zero eta and parseable events/report/snapshots.

- [ ] **Step 4: Confirm strict sample pairing**

Compare the box-wing point IDs and initial absolute states against the current cannonball control group. Expected: IDs and initial physical states match within output precision; only propagated dynamics differ after hour 0.

- [ ] **Step 5: Run regression tests**

Run: `fpm test` and the existing SRP Python test suite. Expected: no new failures.

- [ ] **Step 6: Commit any evidence-only documentation update**

If commands and output schemas match the plan, update the operation guide with measured smoke runtime and file sizes, then commit only that document.

### Task 7: Full 15-day gate before batch launch

**Files:**
- Output only under `20261005SRP/results_box_wing_initial_only/.qualification/`.

**Interfaces:**
- Produces qualification reports for one Earth and one Moon run.

- [ ] **Step 1: Run one eligible orbit for 360 h in Earth mode**

Expected: 361 unique hours and complete hourly/keyframe/event/snapshot artifacts.

- [ ] **Step 2: Audit Earth output**

Require stable 2,048 hourly IDs, 50,000 unique IDs at each keyframe, finite values, matching hashes, explicit ADS validity flags and a model/attitude report of `box_wing/earth`.

- [ ] **Step 3: Repeat for Moon mode**

Apply the identical run and audit with `box_wing/moon`.

- [ ] **Step 4: Record qualification**

Write a JSON summary containing runtimes, sizes, first split/depth/failure events, maximum validation errors and readiness. Mark `batch_ready=true` only if both modes pass structural audits; preserve scientific failure times rather than treating them as missing output.

- [ ] **Step 5: Stop for user review before launching the all-orbit batch**

Present qualification paths and measured resource use. Do not launch the full Earth/Moon batch until the user explicitly requests it.
